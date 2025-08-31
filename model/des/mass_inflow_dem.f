MODULE MASS_INFLOW_DEM_MOD

   USE bc
   USE constant, only: PI
   USE des_psd, only: psd_radius
   USE derived_types, only: dg_pic
   USE derived_types, only: pic
   USE des_allocate, only: particle_grow, particle_grow_gluedspheredem
   USE des_bc
   USE des_physical_prop_mod, only: des_physical_prop
   USE des_rxns, only: des_x_s, des_min_pmass
   USE des_thermo, only: des_t_s
   USE desgrid
   USE discretelement
   USE functions, only: funijk, is_ghost, is_entering, is_exiting
   USE functions, only: is_entering_ghost, is_exiting_ghost, is_nonexistent
   USE functions, only: set_entering, set_entering_ghost, set_ghost, set_normal
   USE functions, only: set_entering_gsp, set_normal_gsp
   use gsp_math_mod, only: SET_MARKED_GSP, SET_MARKED_GSP_STATE, findIndex_1i
   USE geometry
   USE GSP_MATH_MOD, only: gsp_quat_to_exyz
   USE indices
   USE mpi_utility, only: BCAST, global_all_min
   USE param1, only: half, undefined, zero
   USE physprop, only: d_p0, ro_s0, c_ps0
   USE physprop, only: nmax
   USE randomno
   USE resize, only: integer_grow
   USE run, only: ANY_SOLIDS_SPECIES_EQ, ENERGY_EQ, TSTOP
   use sq_properties_mod

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_MASS_INLET                                          !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!  Revised: Hang zhou, Renjie Ke                      Date: 14-Mar-24  !
!                                                                      !
!  Purpose:  This routine fills in the necessary information for new   !
!  particles entering the system.                                      !
!  Revised:  Provide mass inflow support for gsp                       !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE MASS_INFLOW_DEM(lDTSOLID)

      use desmpi_wrapper, only: DES_MPI_STOP
      implicit none

! Dummy Arguments
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: lDTSOLID

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: IP, LS, M, NP, IJK, LC, PP, gid
      INTEGER :: BCV, BCV_I
! number of spheres needed to seed: for GluedSphereDEM, it is the number
! of spheres in each non-spherical particle; for all other models, it is
! one.
      INTEGER :: numSpherSeed
      LOGICAL :: CHECK_FOR_ERRORS, OWNS
! I/J/K index of fluid cell containing the new particle.
      INTEGER :: IJKP(3)
      DOUBLE PRECISION :: DIST, POS(3), qnorm, dist2
! Random numbers -shared by all ranks-
      DOUBLE PRECISION :: RAND(3)
! Max number of glued sphere after pLoop
      INTEGER :: dummy_index, old_NGluedParticles, old_iMAX_GLOBAL_ID, tmp_gid_size
      DOUBLE PRECISION :: TMP_R
      DOUBLE PRECISION, PARAMETER  :: THIRD = (1.0d0/3.0d0)

!......................................................................!

      CHECK_FOR_ERRORS = .FALSE.
      ! all ranks save a copy of current NGluedParticles before going to seed the inflow gsps
      IF(gsp_explicit) THEN
         old_NGluedParticles = NGluedparticles
         old_iMAX_GLOBAL_ID  = iMAX_GLOBAL_ID
      ENDIF

      DO BCV_I = 1, DEM_BCMI
         BCV = DEM_BCMI_MAP(BCV_I)

         DO LC=DEM_BCMI_IJKSTART(BCV_I), DEM_BCMI_IJKEND(BCV_I)
            IJK = DEM_BCMI_IJK(LC)

            DO LS= 1,DG_PIC(IJK)%ISIZE
               NP = DG_PIC(IJK)%P(LS)
               IF(IS_EXITING(NP) .or. IS_EXITING_GHOST(NP)) CYCLE
	            IF (GSP_EXPLICIT) THEN
                  ! note that for mpi, gid_list is global and should be located through iglobal_id
	               gid = gid_list(NP)

                  ! check if gsp is fully inside the domain
               	SELECT CASE (BC_PLANE(BCV))
                  CASE('N'); DIST = gluePos(gid,2) - YN(BC_J_s(BCV)) ! BC_J_s(BCV) is the boundary value
                  CASE('S'); DIST = YN(BC_J_s(BCV)-1) - gluePos(gid,2)
                  CASE('E'); DIST = gluePos(gid,1) - XE(BC_I_w(BCV))
                  CASE('W'); DIST = XE(BC_I_w(BCV)-1) - gluePos(gid,1)
                  CASE('T'); DIST = gluePos(gid,3) - ZT(BC_K_b(BCV))
                  CASE('B'); DIST = ZT(BC_K_b(BCV)-1) - gluePos(gid,3)
                  END SELECT
! The particle is still inside the domain
                  IF(DIST > glueBounding(gid)/2.0D0) THEN
		               IF(IS_ENTERING(NP)) CALL SET_NORMAL(NP)
		               IF(IS_ENTERING_GHOST(NP)) CALL SET_GHOST(NP)
                     CALL SET_NORMAL_GSP(gid)
                     ! save gsp whose gsp_state is changed into marked array
                     ! later, marked array will be used to built the global particle_state_gsp in PE_IO
                     IF(numPEs > 1) THEN
                        ! avoid the same gid being pushed into set marked gsp multiple times
                        dummy_index = findIndex_1i(MARKED_GSP,gid)
                        if(dummy_index /= -1) then
                           marked_state(dummy_index) = particle_state_gsp(gid)
                        else
                           CALL SET_MARKED_GSP(gid)
                           CALL SET_MARKED_GSP_STATE(particle_state_gsp(gid))
                        endif
                     ENDIF
		            ENDIF
	            ELSE
	               SELECT CASE (BC_PLANE(BCV))
                  CASE('N'); DIST = DES_POS_NEW(NP,2) - YN(BC_J_s(BCV))
                  CASE('S'); DIST = YN(BC_J_s(BCV)-1) - DES_POS_NEW(NP,2)
                  CASE('E'); DIST = DES_POS_NEW(NP,1) - XE(BC_I_w(BCV))
                  CASE('W'); DIST = XE(BC_I_w(BCV)-1) - DES_POS_NEW(NP,1)
                  CASE('T'); DIST = DES_POS_NEW(NP,3) - ZT(BC_K_b(BCV))
                  CASE('B'); DIST = ZT(BC_K_b(BCV)-1) - DES_POS_NEW(NP,3)
                  END SELECT
! The particle is still inside the domain
                  TMP_R = DES_RADIUS(NP)
                  !@renjieke, in gsp implicit, use bounding diameter to determine if particle is still in domain
                  !IF(GSP_IMPLICIT) TMP_R = (6*PVOL(NP)/PI)**(THIRD)
                  IF(DIST > TMP_R) THEN
                     IF(IS_ENTERING(NP)) CALL SET_NORMAL(NP)
                     IF(IS_ENTERING_GHOST(NP)) CALL SET_GHOST(NP)
                  ENDIF
	            ENDIF ! end of IF (GluedSphereDEM)
            ENDDO ! end of LS= 1,DG_PIC(IJK)%ISIZE
         ENDDO ! end of LC=DEM_BCMI_IJKSTART(BCV_I), DEM_BCMI_IJKEND(BCV_I)

! All ranks generate random numbers but PE_IO BCASTS its values so that
! the calculated position (with random wiggles) is the same on all ranks.
         CALL RANDOM_NUMBER(RAND)
         CALL BCAST(RAND)

! Check if any particles need seeded this time step.
         IF(DEM_MI_TIME(BCV_I) > S_TIME) CYCLE
         IF((S_TIME<BC_MI_START_TIME(BCV)).OR.(S_TIME>BC_MI_END_TIME(BCV))) CYCLE

         LS = 1
! Loop over the particles being injected
         PLoop: DO IP = 1, PI_COUNT(BCV_I)
            ! before doing the operation, make sure all ranks know current max
            ! also a specific rank might already add 1 to NGluedParticles and imax_global_id
            ! so need to unify their values
            IF(GSP_EXPLICIT) THEN
               CALL global_all_max(NGluedParticles)
            ENDIF
            CALL global_all_max(imax_global_id)
! Determine the location and phase of the incoming particle.
            CALL SEED_NEW_PARTICLE(BCV, BCV_I, RAND, M, POS, IJKP, OWNS)
! Only the rank receiving the new particle needs to continue
            IF(.NOT.OWNS) CYCLE PLoop
	         IF(GSP_EXPLICIT) THEN
	    	      numSpherSeed = sphereNumber(M)
		         NGluedParticles = NGluedParticles + 1
		         CALL PARTICLE_GROW_GluedSphereDEM(NGluedParticles)
               CALL SET_ENTERING_GSP(NGluedParticles)
               IF(numPEs > 1) THEN
                  dummy_index = findIndex_1i(MARKED_GSP,NGluedParticles)
                  if(dummy_index /= -1) then
                        marked_state(dummy_index) = particle_state_gsp(NGluedParticles)
                  else
                        CALL SET_MARKED_GSP(NGluedParticles)
                        CALL SET_MARKED_GSP_STATE(particle_state_gsp(NGluedParticles))
                  endif
               ENDIF

		         IF(BC_GSP_RANDOM_Q(BCV_I,M)) THEN
         	      qnorm = 0
                  DO WHILE (qnorm .gt. 1 .or. qnorm .lt. 1.0d-8)
! If point is outside unit 4-sphere, reject it and try again. This assures uniform distribution.
                     CALL RANDOM_NUMBER(glueQuat(NGluedParticles,:)) ! random number from [0,1]
                     glueQuat(NGluedParticles,:) = 2.0D0*glueQuat(NGluedParticles,:) - 1.0D0 ! shift range to [-1,1]
                     qnorm = glueQuat(NGluedParticles,1)**2 + glueQuat(NGluedParticles,2)**2 + &
                     glueQuat(NGluedParticles,3)**2 + glueQuat(NGluedParticles,4)**2
                  ENDDO
                  glueQuat(NGluedParticles,:) = glueQuat(NGluedParticles,:) /sqrt(qnorm)
               ELSE
                  glueQuat(NGluedParticles,1)=GSP_q1(M)
                  glueQuat(NGluedParticles,2)=GSP_q2(M)
                  glueQuat(NGluedParticles,3)=GSP_q3(M)
                  glueQuat(NGluedParticles,4)=GSP_q4(M)
               ENDIF
		         CALL gsp_quat_to_exyz(glueQuat(NGluedParticles,:), glueEX(NGluedParticles,:), glueEY(NGluedParticles,:), glueEZ(NGluedParticles,:))
		         gluePos(NGluedParticles,:) = POS(:)
		         glueInertia(NGluedParticles, :) = zero
		         if(BC_PSD_TYPE(BCV_I,M) .eq. 'MONO') then
                  glueBounding(NGluedParticles) =  D_P0(M)
               else
                  glueBounding(NGluedParticles) = psd_radius(2, BCV_I,M)*2.0d0    ! 2 indicates BC (1 to indicate IC)
               end if
	         ELSE  !normal dem
	    	      numSpherSeed = 1
	         ENDIF !IF(GluedSphereDEM)

	         DO PP = 1, numSpherSeed
! Increment the global max particle ID (all ranks).
            	iMAX_GLOBAL_ID = iMAX_GLOBAL_ID + 1
! Increment the number of particle on the processor by one. If the max
! number of particles is exceeded, set the error flag and cycle.
            	PIP = PIP + 1
            	CALL PARTICLE_GROW(PIP)
            	MAX_PIP = max(PIP,MAX_PIP)

! Find the first free space in the particle existence array.
            	NP_LP: DO NP = LS, MAX_PIP
               	IF(IS_NONEXISTENT(NP)) THEN
                  	LS = NP
                  	EXIT NP_LP
               	ENDIF
            	ENDDO NP_LP

! Set the particle's global ID.
            	iGLOBAL_ID(LS) = iMAX_GLOBAL_ID
! Set the properties of the new particle.
		      IF(GSP_EXPLICIT) THEN
            	   CALL SET_NEW_PARTICLE_PROPS_GSP(BCV, M, LS, POS, IJKP, PP)
! Update glueVolume
                     glueVolume(NGluedParticles) = glueVolume(NGluedParticles) + PVOL(LS)
! Update glueMass
                     glueMass(NGluedParticles) = glueMass(NGluedParticles) + PMASS(LS)
! Update glueInertia
                     IF (sphereNumber(M) .eq. 1) THEN
	                  glueInertia(NGluedParticles, :) = 2.0d0/5.0d0*PMASS(LS)*DES_RADIUS(LS)**2
                     ELSE
                        dist2 = sc2gpc_vec(LS,2)*sc2gpc_vec(LS,2) + sc2gpc_vec(LS,3)*sc2gpc_vec(LS,3)
                        glueInertia(NGluedParticles,1) = glueInertia(NGluedParticles,1) + PMASS(LS) * dist2+ (2.0d0/5.0d0)*PMASS(LS)*DES_RADIUS(LS)**2

                        dist2 = sc2gpc_vec(LS,1)*sc2gpc_vec(LS,1) + sc2gpc_vec(LS,3)*sc2gpc_vec(LS,3)
	                  glueInertia(NGluedParticles,2) = glueInertia(NGluedParticles,2) + PMASS(LS) * dist2 + (2.0d0/5.0d0)*PMASS(LS)*DES_RADIUS(LS)**2

                        dist2 = sc2gpc_vec(LS,1)*sc2gpc_vec(LS,1) + sc2gpc_vec(LS,2)*sc2gpc_vec(LS,2)
                        glueInertia(NGluedParticles,3) = glueInertia(NGluedParticles,3) + PMASS(LS) * dist2 + (2.0d0/5.0d0)*PMASS(LS)*DES_RADIUS(LS)**2
                     ENDIF
! Update glueVel
                     glueVel(NGluedParticles,1) = BC_U_s(BCV,M)
                     glueVel(NGluedParticles,2) = BC_V_s(BCV,M)
                     glueVel(NGluedParticles,3) = BC_W_s(BCV,M)
	            ELSE
                     CALL SET_NEW_PARTICLE_PROPS(BCV, M, LS, POS, IJKP)
		      ENDIF ! end of IF(GluedSphereDEM)
! Update the minimum particle mass for each solid phases
                  IF(ANY_SOLIDS_SPECIES_EQ) THEN
                     IF(PMASS(LS) .LT. DES_MIN_PMASS(M)) DES_MIN_PMASS(M)=PMASS(LS)
                  ENDIF
	         ENDDO ! end of PP = 1, numSpherSeed
            IF(GSP_EXPLICIT) glueDiameter(NGluedParticles) = 2.0D0*((3.0D0/4.0D0)*glueVolume(NGluedParticles)/PI)**(1.0D0/3.0D0)
         ENDDO PLoop

! after Ploop, all core should have same NGluedParticles and imax_global_id
! but the glued related array might have different size

! Update the time for seeding the next particle.
         DEM_MI_TIME(BCV_I) =  DEM_MI_TIME(BCV_I) + PI_FACTOR(BCV_I)*lDTSOLID
! Set the flag for error checking.
         CHECK_FOR_ERRORS = .TRUE.
      ENDDO ! end of DO BCV_I = 1, DEM_BCMI

      IF(ANY_SOLIDS_SPECIES_EQ) CALL global_all_min(des_min_pmass)

! @renjieke 5-12-2025, global_all_sum birngs trouble when there are two mass inflow BC.
      IF(GSP_EXPLICIT .and. numPEs > 1) THEN
         call global_all_max(NGluedParticles)
         call global_all_max(imax_global_id)
         ! unify the size of glued related variables in each ranks
         CALL PARTICLE_GROW_GluedSphereDEM(NGluedParticles)

         CALL GLOBAL_ALL_SUM(gluePos(old_NGluedParticles+1:NGluedParticles,:))
         CALL GLOBAL_ALL_SUM(glueVel(old_NGluedParticles+1:NGluedParticles,:))
         CALL GLOBAL_ALL_SUM(glueInertia(old_NGluedParticles+1:NGluedParticles,:))
         CALL GLOBAL_ALL_SUM(glueQuat(old_NGluedParticles+1:NGluedParticles,:))
         CALL GLOBAL_ALL_SUM(glueVolume(old_NGluedParticles+1:NGluedParticles))
         CALL GLOBAL_ALL_SUM(glueDiameter(old_NGluedParticles+1:NGluedParticles))
         CALL GLOBAL_ALL_SUM(glueBounding(old_NGluedParticles+1:NGluedParticles))
      ENDIF

      IF(CHECK_FOR_ERRORS) THEN
      ENDIF

      RETURN
   END SUBROUTINE MASS_INFLOW_DEM

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SEED_NEW_PARTICLE                                       !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Purpose:  This routine uses the classification information to place !
!  a new particle in the proper location.                              !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE SEED_NEW_PARTICLE(lBCV, lBCV_I, lRAND, lM, lPOS, &
         lIJKP, lOWNS)
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! The associated bc no.
      INTEGER, INTENT(IN) :: lBCV, lBCV_I
! Random numbers
      DOUBLE PRECISION, INTENT(IN) :: lRAND(3)
! Phase of incoming particle.
      INTEGER, INTENT(OUT) :: lM
! Position of incoming particle.
      DOUBLE PRECISION, INTENT(OUT) :: lPOS(3)
! I/J/K index of fluid cell containing the new particle.
      INTEGER, INTENT(OUT) :: lIJKP(3)
! Logical indicating that the current rank is the owner
      LOGICAL, INTENT(OUT) :: lOWNS

! Local variables
!---------------------------------------------------------------------//
! the associated bc no.
!     INTEGER :: BCV
! Random position
      DOUBLE PRECISION RAND1, RAND2
! Array index of vacant position
      INTEGER :: VACANCY
! Number of array position.
      INTEGER :: OCCUPANTS
      INTEGER :: RAND_I
      INTEGER :: lI, lJ, lK
      DOUBLE PRECISION :: WINDOW
!......................................................................!

!      IF(PARTICLE_PLCMNT(lBCV_I) == 'ORDR')THEN
      VACANCY = DEM_MI(lBCV_I)%VACANCY
      OCCUPANTS = DEM_MI(lBCV_I)%OCCUPANTS
      DEM_MI(lBCV_I)%VACANCY = MOD(VACANCY,OCCUPANTS) + 1
!      ELSE
!         call bcast(lpar_rad)
!         call bcast(lpar_pos)
!         call bcast(m)
!      ENDIF
! Assign the phase of the incoming particle.
      IF(DEM_MI(lBCV_I)%POLYDISPERSE) THEN
         RAND_I = ceiling(dble(NUMFRAC_LIMIT)*lRAND(1))
         lM = DEM_BC_POLY_LAYOUT(lBCV_I, RAND_I)
      ELSE
         lM = DEM_BC_POLY_LAYOUT(lBCV_I,1)
      ENDIF

      WINDOW = DEM_MI(lBCV_I)%WINDOW
      RAND1 = HALF*D_p0(lM) + (WINDOW - D_p0(lM))*lRAND(1)
      RAND2 = HALF*D_p0(lM) + (WINDOW - D_p0(lM))*lRAND(2)

! Set the physical location and I/J/K location of the particle.
      SELECT CASE(BC_PLANE(lBCV))
      CASE('N','S')

         lPOS(1) = DEM_MI(lBCV_I)%P(VACANCY) + RAND1
         lPOS(3) = DEM_MI(lBCV_I)%Q(VACANCY) + RAND2
         lPOS(2) = DEM_MI(lBCV_I)%OFFSET

         lIJKP(1) = DEM_MI(lBCV_I)%W(VACANCY)
         lIJKP(3) = DEM_MI(lBCV_I)%H(VACANCY)
         lIJKP(2) = DEM_MI(lBCV_I)%L

      CASE('E','W')

         lPOS(2) = DEM_MI(lBCV_I)%P(VACANCY) + RAND1
         lPOS(3) = DEM_MI(lBCV_I)%Q(VACANCY) + RAND2
         lPOS(1) = DEM_MI(lBCV_I)%OFFSET

         lIJKP(2) = DEM_MI(lBCV_I)%W(VACANCY)
         lIJKP(3) = DEM_MI(lBCV_I)%H(VACANCY)
         lIJKP(1) = DEM_MI(lBCV_I)%L

      CASE('T','B')

         lPOS(1) = DEM_MI(lBCV_I)%P(VACANCY) + RAND1
         lPOS(2) = DEM_MI(lBCV_I)%Q(VACANCY) + RAND2
         lPOS(3) = DEM_MI(lBCV_I)%OFFSET

         lIJKP(1) = DEM_MI(lBCV_I)%W(VACANCY)
         lIJKP(2) = DEM_MI(lBCV_I)%H(VACANCY)
         lIJKP(3) = DEM_MI(lBCV_I)%L

      END SELECT

! Easier and cleaner to clear out the third component at the end.
      IF(NO_K) lPOS(3) = ZERO

         lI = IofPOS(lPOS(1))
         lJ = JofPOS(lPOS(2))
         lK = KofPOS(lPOS(3))

         lOWNS = ((DG_ISTART <= lI) .AND. (lI <= DG_IEND) .AND. &
            (DG_JSTART <= lJ) .AND. (lJ <= DG_JEND))

         IF(DO_K) lOWNS = lOWNS .AND. (DG_KSTART<=lK) .AND. (lK<=DG_KEND)

      RETURN
   END SUBROUTINE SEED_NEW_PARTICLE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_NEW_PARTICLE_PROPS                                  !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Purpose:  Set the properties of the new particle.                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE SET_NEW_PARTICLE_PROPS(lBCV, lM, lNP, lPOS, lIJKP)

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! The associated bc no.
      INTEGER, INTENT(IN) :: lBCV
! Phase of incoming particle.
      INTEGER, INTENT(IN) :: lM
! Index of new particle
      INTEGER, INTENT(IN) :: lNP
! Position of incoming particle.
      DOUBLE PRECISION, INTENT(IN) :: lPOS(3)
! I/J/K index of fluid cell containing the new particle.
      INTEGER, INTENT(IN) :: lIJKP(3)

! Local variables
!---------------------------------------------------------------------//
! I/J/K index of DES grid cell
      INTEGER :: lI, lJ, lK, n_rows, PP
      DOUBLE PRECISION :: axi(3), m,n, IXX, IYY, IZZ, SVOLUME, qnorm
      DOUBLE PRECISION :: actual_b, sphere_radius, xloc, yloc, zloc
      DOUBLE PRECISION :: sphere_mass_tmp, dist2, Itemp
!......................................................................!

! The particle exists and is entering, not exiting nor a ghost particle
      IF (IS_GHOST(lNP)) THEN
         CALL SET_ENTERING_GHOST(lNP)
      ELSE
         CALL SET_ENTERING(lNP)
      ENDIF

! Set the initial position values based on mass inlet class
      DES_POS_NEW(lNP,:) = lPOS(:)
      PPOS(lNP,:) = lPOS(:)
      DES_VEL_NEW(lNP,1) = BC_U_s(lBCV,lM)
      DES_VEL_NEW(lNP,2) = BC_V_s(lBCV,lM)
      DES_VEL_NEW(lNP,3) = BC_W_s(lBCV,lM)

! Set the initial velocity values
      IF (DO_OLD) THEN
         DES_POS_OLD(lNP,:) = lPOS(:)
         DES_VEL_OLD(lNP,:) = DES_VEL_NEW(lNP,:)
         OMEGA_OLD(lNP,:) = 0
      ENDIF

! Set the initial angular velocity values
      OMEGA_NEW(lNP,:) = 0

! Set the initial angular position values
      IF(PARTICLE_ORIENTATION) ORIENTATION(lNP,1:3) = INIT_ORIENTATION

! Set the particle radius value
      if(BC_PSD_TYPE(lBCV,lM) .eq. 'MONO') then
        DES_RADIUS(lNP) = D_P0(lM)*HALF
	! d_p0 was already scaled to parcel size
        IF(CGP_SCALING_METHOD(lM)==2) then
           DES_RADIUS(lNP) = DES_RADIUS(lNP)/CGP_STAT_WT(lM)**(1.0d0/3.0d0)
        endif
      else
        DES_RADIUS(lNP) = psd_radius(2, lBCV, lM)    ! 2 indicates BC (1 to indicate IC)
      end if


      if(CGDEM) then
         IF(CGP_SCALING_METHOD(lM)==1) THEN
            DES_CGP_STW(lNP) = CGP_STAT_WT(lM)
            DES_CGP_RPR(lNP) = DES_RADIUS(lNP)/DES_CGP_STW(lNP)**(1.0d0/3.0d0)
         ELSEIF(CGP_SCALING_METHOD(lM)==2) THEN
            DES_CGP_STW(lNP) = (HALF*CGP_D_P0(lM)/DES_RADIUS(lNP))**3
            DES_CGP_RPR(lNP) = DES_RADIUS(lNP)
            DES_RADIUS(lNP) = HALF*CGP_D_P0(lM)
         ENDIF
      endif

! Set the particle density value
      RO_Sol(lNP) = RO_S0(lM)

! Store the I/J/K indices of the particle.
      PIJK(lNP,1:3) = lIJKP(:)
      PIJK(lNP,4) = FUNIJK(lIJKP(1), lIJKP(2), lIJKP(3))

! Set the particle mass phase
      PIJK(lNP,5) = lM

! Calculate the DES grid cell indices.
      lI = min(DG_IEND2,max(DG_ISTART2,IOFPOS(DES_POS_NEW(lNP,1))))
      lJ = min(DG_JEND2,max(DG_JSTART2,JOFPOS(DES_POS_NEW(lNP,2))))
      IF(NO_K) THEN
         lK = 1
      ELSE
         lK = min(DG_KEND2,max(DG_KSTART2,KOFPOS(DES_POS_NEW(lNP,3))))
      ENDIF
! Store the triple
      DG_PIJK(lNP) = DG_FUNIJK(lI,lJ,lK)

! Calculate the new particle's Volume, Mass, OMOI3
! SuperDEM
         IF (SuperDEM)   then
            super_r(lNP,1)=SQP_a(lM) * 2.0D0 * DES_RADIUS(lNP) / D_P0(lM) ! scale a,b,c for polydisperse systems
            super_r(lNP,2)=SQP_b(lM) * 2.0D0 * DES_RADIUS(lNP) / D_P0(lM)
            super_r(lNP,3)=SQP_c(lM) * 2.0D0 * DES_RADIUS(lNP) / D_P0(lM)
            super_mn(lNP,1)=SQP_m(lM)
            super_mn(lNP,2)=SQP_n(lM)

            IF(BC_SQP_RANDOM_Q(lBCV,lM)) THEN
               qnorm = 0
               DO WHILE (qnorm .gt. 1 .or. qnorm .lt. 1.0d-8)
                   ! If point is outside unit 4-sphere, reject it and try again. This assures uniform distribution.
                   CALL RANDOM_NUMBER(super_q(lNP,:))            ! random number from [0,1]
                   super_q(lNP,:) = 2.0D0*super_q(lNP,:) - 1.0D0 ! shift range to [-1,1]
                   qnorm = super_q(lNP,1)**2 + super_q(lNP,2)**2 + super_q(lNP,3)**2 + super_q(lNP,4)**2
               END DO
               super_q(lNP,:) = super_q(lNP,:) / SQRT(qnorm)
            ELSE
               super_q(lNP,1)=SQP_q1(lM)
               super_q(lNP,2)=SQP_q2(lM)
               super_q(lNP,3)=SQP_q3(lM)
               super_q(lNP,4)=SQP_q4(lM)
            ENDIF

            axi(1)=super_r(lNP,1)
            axi(2)=super_r(lNP,2)
            axi(3)=super_r(lNP,3)
            m=super_mn(lNP,1)
            n=super_mn(lNP,2)
            CALL SQ_VOLUME (axi,m,n,svolume)
            PVOL(lNP)=svolume
            PMASS(lNP) = PVOL(lNP)*RO_SOL(lNP)
            call SQ_INERTIA(axi,m,n,IXX,IYY,IZZ)
            OMOI3(lNP,1)=1.0/(IXX*RO_SOL(lNP))
            OMOI3(lNP,2)=1.0/(IYY*RO_SOL(lNP))
            OMOI3(lNP,3)=1.0/(IZZ*RO_SOL(lNP))
! GSP implicit modeling
         ELSEIF (GSP_IMPLICIT) THEN
            IF(BC_GSP_RANDOM_Q(lBCV,lM)) THEN
               qnorm = 0
               DO WHILE (qnorm .gt. 1 .or. qnorm .lt. 1.0d-8)
                  CALL RANDOM_NUMBER(glueQuat(lNP,:))
                  glueQuat(lNP,:) = 2.0D0*glueQuat(lNP,:) - 1.0D0
                  qnorm = glueQuat(lNP,1)**2 + glueQuat(lNP,2)**2 + glueQuat(lNP,3)**2 + glueQuat(lNP,4)**2
               ENDDO
               glueQuat(lNP,:) = glueQuat(lNP,:)/sqrt(qnorm)
            ELSE
               glueQuat(lNP,1)=GSP_q1(lM)
               glueQuat(lNP,2)=GSP_q2(lM)
               glueQuat(lNP,3)=GSP_q3(lM)
               glueQuat(lNP,4)=GSP_q4(lM)
            ENDIF
            IF (BC_PSD_TYPE(lBCV,lM) .eq. 'MONO') THEN
               PVOL(lNP) = gsp_vol(lM)
            ELSE
               PVOL(lNP) = (4.0D0/3.0D0)*PI*DES_RADIUS(lNP)**3
               PVOL(lNP) = PVOL(lNP)*gsp_vol(lM)/((pi*d_p0(lM)**3.0D0)/6.0d0)
            ENDIF
            PMASS(lNP) = PVOL(lNP)*RO_SOL(lNP)
            glueSurface = 0.0D0
            CALL gsp_quat_to_exyz(glueQuat(lNP, :), glueEX(lNP,:), glueEY(lNP,:), glueEZ(lNP,:))

            IF(sphereNumber(lM) .eq. 1) THEN
               OMOI3(lNP, :) = 2.5D0/(PMASS(lNP)*DES_RADIUS(lNP)**2)
            ELSE
               OMOI3(lNP, :) = 0.0d0
               n_rows = size(cspart_data,1)
               actual_b = 2.0 * des_radius(lNP)
               DO PP = 1, n_rows
                  if(cspart_data(PP,4,lM) == 0.0) cycle
                  sphere_radius = cspart_data(PP,4,lM) / 2.0 * actual_b
                  xloc = cspart_data(PP,1,lM) * actual_b
                  yloc = cspart_data(PP,2,lM) * actual_b
                  zloc = cspart_data(PP,3,lM) * actual_b
                  glueSurface(lNP) = glueSurface(lNP) + cspart_data(PP,5,lM) * actual_b * actual_b

                  sphere_mass_tmp = RO_SOL(lNP) * PI / 6.0D0 * sphere_radius ** 3.0D0
                  Itemp = (2.0/5.0)*sphere_mass_tmp*sphere_radius**2.0

                  dist2 = yloc ** 2.0 + zloc ** 2.0
                  OMOI3(lNP,1) = OMOI3(lNP,1) + sphere_mass_tmp * dist2 + Itemp

                  dist2 = xloc ** 2.0 + zloc ** 2.0
                  OMOI3(lNP,2) = OMOI3(lNP,2) + sphere_mass_tmp * dist2 + Itemp

                  dist2 = xloc ** 2.0 + yloc ** 2.0
                  OMOI3(lNP,3) = OMOI3(lNP,3) + sphere_mass_tmp * dist2 + Itemp
               ENDDO
               OMOI3(lNP, 1) = 1.0d0/OMOI3(lNP, 1)
               OMOI3(lNP, 2) = 1.0d0/OMOI3(lNP, 2)
               OMOI3(lNP, 3) = 1.0d0/OMOI3(lNP, 3)

               glueAngMom(lNP, 1) = OMEGA_NEW(lNP,1)/OMOI3(lNP,1)
               glueAngMom(lNP, 2) = OMEGA_NEW(lNP,2)/OMOI3(lNP,2)
               glueAngMom(lNP, 3) = OMEGA_NEW(lNP,3)/OMOI3(lNP,3)
            ENDIF
         ELSE
            PVOL(lNP) = (4.0D0/3.0D0)*PI*DES_RADIUS(lNP)**3
            PMASS(lNP) = PVOL(lNP)*RO_SOL(lNP)
! for spherical particle, IXX=IYY=IZZ
            OMOI(lNP) = 2.5D0/(PMASS(lNP)*DES_RADIUS(lNP)**2) !ONE OVER MOI
         ENDIF


! Clear the drag force
      IF(DES_EXPLICITLY_COUPLED) then
         F_GP(lNP) = ZERO
         DRAG_FC(lNP,:) = ZERO
      ENDIF

! If solving the energy equations, set the temperature
      IF(ANY_SOLIDS_SPECIES_EQ .OR. ENERGY_EQ ) DES_T_s(lNP) = BC_T_s(lBCV,lM)

! Set species mass fractions
      IF((ENERGY_EQ .AND. C_PS0(lM)==UNDEFINED) .OR. ANY_SOLIDS_SPECIES_EQ)&
         DES_X_s(lNP,1:NMAX(lM)) = BC_X_s(lBCV,lM,1:NMAX(lM))

! Residence time
      RESIDENCE_TIME(lNP) = ZERO

! Calculate time dependent physical properties
      CALL DES_PHYSICAL_PROP


      RETURN
   END SUBROUTINE SET_NEW_PARTICLE_PROPS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_NEW_PARTICLE_PROPS_GSP                              !
!  Author: Hang Zhou                                   Date: 16-Jan-24 !
!                                                                      !
!  Purpose:  Set the properties of the new particle for GSP.           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE SET_NEW_PARTICLE_PROPS_GSP(lBCV, lM, lNP, lPOS, lIJKP, PP)

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! The associated bc no.
      INTEGER, INTENT(IN) :: lBCV
! Phase of incoming particle.
      INTEGER, INTENT(IN) :: lM
! Index of new particle
      INTEGER, INTENT(IN) :: lNP
! Position of incoming particle.
      DOUBLE PRECISION, INTENT(IN) :: lPOS(3)
! I/J/K index of fluid cell containing the new particle.
      INTEGER, INTENT(IN) :: lIJKP(3)
! Index of spheres in the list of spheres in a non-spherical particle
      INTEGER, INTENT(IN) :: PP

! Local variables
!---------------------------------------------------------------------//
! I/J/K index of DES grid cell
      INTEGER :: lI, lJ, lK, PIP_START
      DOUBLE PRECISION :: axi(3), m,n, IXX, IYY, IZZ, SVOLUME, qnorm
!......................................................................!
! Dist2, for inertia calculation
      DOUBLE PRECISION :: dist2

! The particle exists and is entering, not exiting nor a ghost particle
      IF (IS_GHOST(lNP)) THEN
         CALL SET_ENTERING_GHOST(lNP)
      ELSE
         CALL SET_ENTERING(lNP)
      ENDIF
      gid_list(lNP) = NGluedParticles

! read heat conduction related parameters
      gp_sa(lNP) = tmp_gp_sa(start_index_sphere(lM)+PP-1) * glueBounding(NGluedParticles) * glueBounding(NGluedParticles)
! already reset to relative id in read_particle_input
      gp_neighsid(lNP,:) = tmp_gp_neighsid(start_index_sphere(lM)+PP-1,:)
      PIP_START = lNP + 1
      where (gp_neighsid(lNP,:) /= -1) gp_neighsid(lNP,:) = gp_neighsid(lNP,:) + PIP_START
      gp_neighsa(lNP,:) = tmp_gp_neighsa(start_index_sphere(lM)+PP-1,:) * glueBounding(NGluedParticles) * glueBounding(NGluedParticles)

! Set the initial angular position values
      IF(PARTICLE_ORIENTATION) ORIENTATION(lNP,1:3) = INIT_ORIENTATION

! Set the particle radius value
      DES_RADIUS(lNP) = glueBounding(NGluedParticles)*HALF*sphereRelDia(start_index_sphere(lM)+PP-1)
! relative position of spheres to glued particle
      sc2gpc_vec(lNP,:) = glueBounding(NGluedParticles)*sphereRelPos(start_index_sphere(lM)+PP-1,:)
      ! Set the initial position values based on mass inlet class
      DES_POS_NEW(lNP,:) = lPOS(:)&
      			  +sc2gpc_vec(lNP,1)*glueEX(NGluedParticles,:) &
                          +sc2gpc_vec(lNP,2)*glueEY(NGluedParticles,:) &
                          +sc2gpc_vec(lNP,3)*glueEZ(NGluedParticles,:)
      PPOS(lNP,:) = lPOS(:)&
      		    +sc2gpc_vec(lNP,1)*glueEX(NGluedParticles,:) &
                    +sc2gpc_vec(lNP,2)*glueEY(NGluedParticles,:) &
                    +sc2gpc_vec(lNP,3)*glueEZ(NGluedParticles,:)
      ! Set the initial velocity values
      DES_VEL_NEW(lNP,1) = BC_U_s(lBCV,lM)
      DES_VEL_NEW(lNP,2) = BC_V_s(lBCV,lM)
      DES_VEL_NEW(lNP,3) = BC_W_s(lBCV,lM)

      ! Set the initial angular velocity values
      OMEGA_NEW(lNP,:) = 0.0d0


      IF (DO_OLD) THEN
         DES_POS_OLD(lNP,:) = DES_POS_NEW(lNP,:)
         DES_VEL_OLD(lNP,:) = DES_VEL_NEW(lNP,:)
         OMEGA_OLD(lNP,:) = 0.0d0
      ENDIF


! Set the particle density value
      RO_Sol(lNP) = RO_S0(lM)

! Store the I/J/K indices of the particle.
      PIJK(lNP,1:3) = lIJKP(:)
      PIJK(lNP,4) = FUNIJK(lIJKP(1), lIJKP(2), lIJKP(3))

! Set the particle mass phase
      PIJK(lNP,5) = lM

! Calculate the DES grid cell indices.
      lI = min(DG_IEND2,max(DG_ISTART2,IOFPOS(DES_POS_NEW(lNP,1))))
      lJ = min(DG_JEND2,max(DG_JSTART2,JOFPOS(DES_POS_NEW(lNP,2))))
      IF(NO_K) THEN
         lK = 1
      ELSE
         lK = min(DG_KEND2,max(DG_KSTART2,KOFPOS(DES_POS_NEW(lNP,3))))
      ENDIF
! Store the triple
      DG_PIJK(lNP) = DG_FUNIJK(lI,lJ,lK)

! Calculate the new particle's Volume, Mass, OMOI3
      PVOL(lNP) = (4.0D0/3.0D0)*PI*DES_RADIUS(lNP)**3
      PMASS(lNP) = PVOL(lNP)*RO_SOL(lNP)
! for spherical particle, IXX=IYY=IZZ
      OMOI(lNP) = 2.5D0/(PMASS(lNP)*DES_RADIUS(lNP)**2) !ONE OVER MOI

! Clear the drag force
      IF(DES_EXPLICITLY_COUPLED) then
         F_GP(lNP) = ZERO
         DRAG_FC(lNP,:) = ZERO
      ENDIF

! If solving the energy equations, set the temperature
      IF(ANY_SOLIDS_SPECIES_EQ .OR. ENERGY_EQ ) DES_T_s(lNP) = BC_T_s(lBCV,lM)

! Set species mass fractions
      IF((ENERGY_EQ .AND. C_PS0(lM)==UNDEFINED) .OR. ANY_SOLIDS_SPECIES_EQ)&
         DES_X_s(lNP,1:NMAX(lM)) = BC_X_s(lBCV,lM,1:NMAX(lM))

! Residence time
      RESIDENCE_TIME(lNP) = ZERO

! Calculate time dependent physical properties
      CALL DES_PHYSICAL_PROP


      RETURN
   END SUBROUTINE SET_NEW_PARTICLE_PROPS_GSP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine:  DES_NEW_PARTICLE_TEST                                  !
!                                                                      !
!  Purpose:  This routine checks if a new particle placed using the    !
!  random inlet was placed in contact with an existing particle.  If   !
!  so a flag is set indicating contact, and the new particle is        !
!  repositioned within the inlet domain.                               !
!                                                                      !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Purpose: This routine has to be modified for parallel version       !
!           the parameter now accepts the lpar_rad and lpar_pos and    !
!           tests if it touches any particles                          !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE DES_NEW_PARTICLE_TEST(BCV_I,ppar_rad,ppar_pos,TOUCHING)

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! index of boundary condition
      INTEGER, INTENT(IN) :: BCV_I
      DOUBLE PRECISION, INTENT(IN) :: ppar_pos(DIMN)
      DOUBLE PRECISION, INTENT(IN) :: ppar_rad
      LOGICAL, INTENT(INOUT) :: TOUCHING

! Local variables
!---------------------------------------------------------------------//
! particle number id of a potential overlapping/contacting particle
      INTEGER NP2
! total number of particles in current ijk cell and loop counter
      INTEGER NPG, LL
! i, j, k indices along boundary used for loop counters
      INTEGER I, J, K, IJK
! for parallel processing
      integer listart,liend,ljstart,ljend,lkstart,lkend

      DOUBLE PRECISION  DISTVEC(DIMN), DIST, R_LM
!......................................................................!

      TOUCHING = .FALSE.

! For parallel processing the arrays has to be limited
!      select case (des_mi_class(bcv_i))
!      case ('XW','XE', 'YZw','YZe')
!         listart = gs_array(bcv_i,1)
!         liend = gs_array(bcv_i,2)
!         ljstart = max(gs_array(bcv_i,3),jstart)
!         ljend = min(gs_array(bcv_i,4),jend)
!         lkstart = max(gs_array(bcv_i,5),jstart)
!         lkend = min(gs_array(bcv_i,6),jend)
!      case ('YN','YS', 'XZn','XZs')
!         listart = max(gs_array(bcv_i,1),istart)
!         liend = min(gs_array(bcv_i,2),iend)
!         ljstart = gs_array(bcv_i,3)
!         ljend = gs_array(bcv_i,4)
!         lkstart = max(gs_array(bcv_i,5),jstart)
!         lkend = min(gs_array(bcv_i,6),jend)
!      case ('ZT','ZB', 'XYt','XYb')
!         listart = max(gs_array(bcv_i,1),istart)
!         liend = min(gs_array(bcv_i,2),iend)
!         ljstart = max(gs_array(bcv_i,3),jstart)
!         ljend = min(gs_array(bcv_i,4),jend)
!         lkstart = gs_array(bcv_i,5)
!         lkend = gs_array(bcv_i,6)
!      end select

      listart = 1
      ljstart = 1
      lkstart = 1
      liend = 1
      ljend = 1
      lkend = 1
      DO k = lkstart,lkend
      DO j = ljstart,ljend
      DO i = listart,liend
!      DO K = GS_ARRAY(BCV_I,5), GS_ARRAY(BCV_I,6)
!         DO J = GS_ARRAY(BCV_I,3), GS_ARRAY(BCV_I,4)
!           DO I =  GS_ARRAY(BCV_I,1), GS_ARRAY(BCV_I,2)
             IJK = FUNIJK(I,J,K)
             IF(ASSOCIATED(PIC(IJK)%P)) THEN
               NPG =  SIZE(PIC(IJK)%P)
               DO LL = 1, NPG
                  NP2 = PIC(IJK)%P(LL)
                  DISTVEC(:) = ppar_pos(:) - DES_POS_NEW(NP2,:)
                  DIST = DOT_PRODUCT(DISTVEC,DISTVEC)
                  R_LM = ppar_rad + DES_RADIUS(NP2)
                  IF(DIST .LE. R_LM*R_LM) TOUCHING = .TRUE.
               ENDDO
             ENDIF
           ENDDO
         ENDDO
       ENDDO

      RETURN

   END SUBROUTINE DES_NEW_PARTICLE_TEST

END MODULE MASS_INFLOW_DEM_MOD
