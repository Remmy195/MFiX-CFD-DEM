
MODULE CFNEWVALUES_MOD
      USE SQ_OBB_MOD
      USE SQ_CONTACT_NEWTON_DPmethod_MOD
      USE SQ_ROTATION_MOD
      USE SQ_PROPERTIES_MOD
      USE Sq_math_mod
      USE SQ_CONTACT_WALL
      USE GSP_MATH_MOD

CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: UPDATE_GSP_PARTICLE_STATE                                !
! Author: Renjie Ke                                   Date Mar-15-2024 !
! Purpose: synchronize the gsp particle state before doing the         !
! half time step update                                                !
! TBD:                                                                 !
!     Even only the rank has gsp state change send the date to pe_io   !
!     communication in each dem time step is still consider costly     !
!     consider build the gsp_particle_state as local array             !
!     in the future                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

   SUBROUTINE UPDATE_GSP_PARTICLE_STATE
   ! use des_mpi
      use cdist
      use mpi_comm_des
      use mpi_utility
      use discretelement
   ! Number of TFM solid phases
      use physprop, only:  MMAX
   ! Error manager
      use error_manager

      IMPLICIT NONE

      INTEGER :: LOCAL_CNT, GLOBAL_CNT
      INTEGER :: L
      INTEGER, ALLOCATABLE :: ltemp_array(:)  ! local
      INTEGER, ALLOCATABLE :: gtemp_gid_array(:)  ! global
      INTEGER, ALLOCATABLE :: gtemp_state_array(:)  ! global
   ! Variables related to gather
      integer :: lgathercnts(0:numpes-1), lproc
   ! Logical
      logical :: isChanged

      LOCAL_CNT = 0

      if(NumPEs<=1) return

   ! if a rank does not have particle_state_gsp change, then its marked_gsp(:) = 0
      IF(ALL(marked_gsp == 0)) THEN
         isChanged = .False.
      ELSE
         isChanged = .True.
      ENDIF

   ! LBOUND is always 1, so find the first available space
      DO L = 1, size(MARKED_GSP)
         IF(MARKED_GSP(L) .eq. 0) EXIT
      ENDDO

   ! marked_state and marked_gsp have the same local_cnt
      LOCAL_CNT = L - 1
      call global_sum(LOCAL_CNT, GLOBAL_CNT)
      igath_sendcnt = LOCAL_CNT
      if(.NOT. isChanged) igath_sendcnt = 0

      lgathercnts = 0
      lgathercnts(myPE) = LOCAL_CNT
      if(.NOT. isChanged) lgathercnts(myPE) = 0
      call global_sum(lgathercnts,igathercnts)
   ! Calculate the rank displacements.
      idispls(0) = 0
      DO lPROC = 1,NUMPEs-1
         idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)
      ENDDO

      ALLOCATE (iprocbuf(LOCAL_CNT))
      ALLOCATE (ltemp_array(LOCAL_CNT))
      if(mype == pe_io) then
         allocate(gtemp_gid_array(GLOBAL_CNT))
         allocate(gtemp_state_array(GLOBAL_CNT))
         allocate(irootbuf(GLOBAL_CNT))
      else
         allocate(gtemp_gid_array(10))
         allocate(gtemp_state_array(10))
         allocate(irootbuf(10))
      endif

   ! Pack marked_gsp in a temporary local array
   ! even if nothing changed in current rank, still assign a valid iprocbuf, but it does not matter
      ltemp_array(1:LOCAL_CNT) = marked_gsp(1:LOCAL_CNT)
      iprocbuf(1:LOCAL_CNT)=ltemp_array(1:LOCAL_CNT)

      CALL desmpi_gatherv(ptype=1)
      IF(MYPE == PE_IO) gtemp_gid_array(:) = irootbuf(:)

      ltemp_array(1:LOCAL_CNT) = marked_state(1:LOCAL_CNT)
      iprocbuf(1:LOCAL_CNT)=ltemp_array(1:LOCAL_CNT)

      CALL desmpi_gatherv(ptype=1)
      IF(MYPE == PE_IO) gtemp_state_array(:) = irootbuf(:)

   ! unpack gtemp_gid_array(:) and gtemp_state_array(:) to update the particle_state_gsp
      if(mype == pe_io) then
         do L = 1,size(gtemp_gid_array)
               ! even if only the non-zero marked_gsp elements are sent, still put a safe check
               if(gtemp_gid_array(L) .lt. 1) cycle
               particle_state_gsp(gtemp_gid_array(L)) = gtemp_state_array(L)
         enddo
      endif

   ! after pe_io finish its update for particle_state_gsp, send back to other ranks
      CALL bcast(particle_state_gsp, pe_io)

   ! then reset the local marked_gsp and marked_gsp_state in each rank
      marked_gsp(:) = 0
      marked_state(:) = -1

      deallocate (iProcBuf, iRootBuf, ltemp_array, gtemp_gid_array, gtemp_state_array)
      ! to this end, global-wise particle_state_gsp has been successfully updated!
   END SUBROUTINE UPDATE_GSP_PARTICLE_STATE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: GSP_CFNEWVALUES_vhalfxfull                             !
!  Reviewer: Renjie Ke                                Date: Mar-18-24  !
!  Purpose: Verlet algorithm, update velocity by half time step        !
!           then update the position by a full step                    !
!           update on glued particle itself                            !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

   SUBROUTINE GSP_CFNEWVALUES_vhalfxfull

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE constant
      USE des_bc
      USE discretelement
      USE fldvar
      USE functions
      USE mfix_pic
      USE mpi_utility
      USE parallel
      USE run
      USE param
      USE param1
      USE physprop
      use mpi_funs_des, only: DES_PAR_EXCHANGE
      USE desgrid, only: desgrid_pic

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER ::L
      DOUBLE PRECISION :: quat(4), exone(3), eyone(3), ezone(3), torloc(3), tbody(3), mbody(3),rot(3,3),local_omg(3)
      DOUBLE PRECISION :: conjqm(4), fquat(4), inertia(3), angmom(3), omega(3)
      double precision, allocatable, dimension(:) :: tmpgluemass
      INTEGER :: gid, TEMP_PIP

!-----------------------------------------------
! step 0, set glued particle mass to 0
      glueMass(:) = 0.d0
      temp_pip = pip
! step 1, recompute glueMass due to the mass outflow condition - should always reset before!
      DO L = 1, MAX_PIP
         IF(IS_NONEXISTENT(L)) CYCLE
         IF(IS_ANY_GHOST(L)) CYCLE
         gid = gid_list(L)
         IF(IS_NORMAL(L) .AND. (.NOT. (PARTICLE_STATE_GSP(gid)==NORMAL_PARTICLE))) THEN
            PARTICLE_STATE(L)=PARTICLE_STATE_GSP(gid)
            ! @renjieke 2025-06-11 a component particle state is changed based on the GSP particle state.
            ! update pip here to make sure the restart file gets the correct number of particles.
            IF(PARTICLE_STATE(L) == 0) TEMP_PIP = TEMP_PIP - 1
         ENDIF
         IF(gid .lt. 1) CYCLE
         glueMass(gid) = glueMass(gid) + PMASS(L)
      ENDDO

! Mass is changed so need to communicate with other PE
      IF(numPEs > 1) THEN
         CALL GLOBAL_ALL_SUM(glueMass)
         PIP = TEMP_PIP
      ENDIF

      IF(mype == pe_io) THEN
! Note: "force & torque are already summarized in GSP_CFNEWVALUES_vfull"
! step 2, calculate velocity $ position of glued particle
         DO gid = 1, NGluedParticles
            IF(PARTICLE_STATE_GSP(gid)== NONEXISTENT) CYCLE
            ! acceleration
            IF(glueMass(gid) .ne. 0) THEN
               ! half dt velocity and full dt gluedparticle position approximation
               if(PARTICLE_STATE_GSP(gid) == NORMAL_PARTICLE) then
                  glueAcc(gid,:) = glueForce(gid,:)/glueMass(gid) + GRAV(:)
                  glueVel(gid,:) = glueVel(gid,:) + glueAcc(gid,:)*dtsolid/2
               endif
               gluePos(gid,:) = gluePos(gid,:) + glueVel(gid,:)*dtsolid
            ENDIF
         ENDDO

! step 3, update rotation
         DO gid = 1, NGluedParticles
            IF(PARTICLE_STATE_GSP(gid)== NONEXISTENT) CYCLE
            angmom = glueAngMom(gid,:)
            quat   = glueQuat(gid,:)
            inertia= glueInertia(gid,:)

            CALL gsp_quat_to_mat(quat,rot)
            CALL gsp_transpose_matvec(rot,angmom,mbody)
            CALL gsp_quatvec(quat,mbody,conjqm)
            conjqm = 2 * conjqm
            torloc(:) = glueTorque(gid,:)

            CALL gsp_transpose_matvec(rot,torloc,tbody)
            CALL gsp_quatvec(quat,tbody,fquat)
            conjqm = conjqm + dtsolid*fquat

            CALL gsp_no_squish_rotate(3,conjqm,quat,inertia,dtsolid/2)
            CALL gsp_no_squish_rotate(2,conjqm,quat,inertia,dtsolid/2)
            CALL gsp_no_squish_rotate(1,conjqm,quat,inertia,dtsolid)
            CALL gsp_no_squish_rotate(2,conjqm,quat,inertia,dtsolid/2)
            CALL gsp_no_squish_rotate(3,conjqm,quat,inertia,dtsolid/2)

            CALL gsp_qnormalize(quat)
            CALL gsp_quat_to_mat(quat,rot)
            CALL gsp_invquatvec(quat,conjqm,mbody)
            mbody = mbody * 0.5
            local_omg = mbody / inertia

            exone(:) = glueEX(gid,:)
            eyone(:) = glueEY(gid,:)
            ezone(:) = glueEZ(gid,:)
            CALL gsp_quat_to_exyz(quat,exone,eyone,ezone)
            glueEX(gid,:) = exone(:)
            glueEY(gid,:) = eyone(:)
            glueEZ(gid,:) = ezone(:)

            CALL gsp_matvec(rot,mbody,angmom)
            CALL gsp_matvec(rot,local_omg,omega)

            glueAngMom(gid,:) = angmom
            glueQuat(gid,:) = quat
            glueOmg(gid,:) = omega
         ENDDO
      ENDIF !(mype == pe_io)

! after updating glue-related array in pe_io, broadcast back to each processor
      IF(numPEs > 1) THEN
         CALL BCAST(gluePos,PE_IO)
         CALL BCAST(glueVel,PE_IO)

         CALL BCAST(glueEX, PE_IO)
         CALL BCAST(glueEY, PE_IO)
         CALL BCAST(glueEZ, PE_IO)

         CALL BCAST(glueAngMom, PE_IO)
         CALL BCAST(glueOmg, PE_IO)
      ENDIF

! step 4, update individual spheres inside the glued particle
! each pe only update normal component sphere, ignore ghost and nonexistent component spheres
      DO L = 1, MAX_PIP
         IF(IS_NONEXISTENT(L)) CYCLE
         IF(IS_ANY_GHOST(L)) CYCLE
         gid = gid_list(L)
         IF(gid .lt. 1) CYCLE
         DES_VEL_NEW(L,:) = glueVel(gid,:)
         DES_POS_NEW(L,:) = gluePos(gid,:) + sc2gpc_vec(L,1)*glueEX(gid,:) &
                     + sc2gpc_vec(L,2)*glueEY(gid,:) + sc2gpc_vec(L,3)*glueEZ(gid,:)
         gp_squat(L,1:4) = glueQuat(gid,:) !TBD only used for output
      ENDDO

      ! IF(numPEs > 1) THEN
      ! ! rebin to desgrid is needed before calling des_par_exchange
      ! CALL DESGRID_PIC(.TRUE.)
      ! CALL DES_PAR_EXCHANGE
      ! ENDIF

      RETURN

   END SUBROUTINE GSP_CFNEWVALUES_vhalfxfull

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: GSP_CFNEWVALUES_vfull                                  !
!  Reviewer: Renjie Ke                                Date: Mar-18-24  !
!                                                                      !
!  Purpose: Verlet algorithm,                                          !
!           update another half step velocity                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

   SUBROUTINE GSP_CFNEWVALUES_vfull

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE constant
      USE des_bc
      USE discretelement
      USE fldvar
      USE mfix_pic
      USE mpi_utility
      USE parallel
      USE run
      USE param
      USE param1
      USE physprop
      USE functions, only: is_nonexistent, is_ghost, is_entering_ghost, is_any_ghost

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER ::L, tempPhaseNumber
      DOUBLE PRECISION :: quat(4),exone(3),eyone(3),ezone(3),torloc(3),tbody(3),mbody(3),rot(3,3),local_omg(3)
      DOUBLE PRECISION :: conjqm(4),fquat(4),inertia(3),angmom(3),omega(3),dist2
      INTEGER :: gid
      double precision, allocatable, dimension(:) :: tmpgluemass
      double precision, allocatable, dimension(:,:) :: tmpglueInertia, tmpglueForce, tmpglueTorque

      integer, save :: totalcnt = 0
!-----------------------------------------------

      !step0, set to 0
      glueMass(:) = 0.d0
      glueInertia(:,:) =0.0D0

      ! step1, sum gas-solid force from individual particles to glued particles
      ! glueForce already includes p-w and p-p collisions
      ! mass can change due to chemical reaction, so recalculate glueMass
      DO L = 1, MAX_PIP
         IF(IS_NONEXISTENT(L)) CYCLE
         IF(IS_ANY_GHOST(L)) CYCLE
         gid = gid_list(L)
         IF(gid .lt. 1) CYCLE
         glueMass(gid) = glueMass(gid) + PMASS(L)
      ENDDO

      if(numPEs > 1) CALL GLOBAL_ALL_SUM(glueMass)

      DO L = 1, MAX_PIP
         IF(IS_NONEXISTENT(L)) CYCLE
         IF(IS_ANY_GHOST(L)) CYCLE
         gid = gid_list(L)
         IF(gid .lt. 1) CYCLE
         tempPhaseNumber = pijk(L, 5)
         ! dist2: distance square
         ! update the inertia
         IF (allocated(sphereNumber) .and. sphereNumber(tempPhaseNumber) .eq. 1) THEN
            ! special case when glued particles only have 1 component sphere
            glueInertia(gid, :) = 2.0d0/5.0d0*PMASS(L)*DES_RADIUS(L)**2
         ELSE
            dist2 = sc2gpc_vec(L,2)*sc2gpc_vec(L,2) + sc2gpc_vec(L,3)*sc2gpc_vec(L,3)
            glueInertia(gid,1) = glueInertia(gid,1) + PMASS(L) * dist2 + (2.0d0/5.0d0)*PMASS(L)*DES_RADIUS(L)**2

            dist2 = sc2gpc_vec(L,1)*sc2gpc_vec(L,1) + sc2gpc_vec(L,3)*sc2gpc_vec(L,3)
            glueInertia(gid,2) = glueInertia(gid,2) + PMASS(L) * dist2 + (2.0d0/5.0d0)*PMASS(L)*DES_RADIUS(L)**2

            dist2 = sc2gpc_vec(L,2)*sc2gpc_vec(L,2) + sc2gpc_vec(L,1)*sc2gpc_vec(L,1)
            glueInertia(gid,3) = glueInertia(gid,3) + PMASS(L) * dist2 + (2.0d0/5.0d0)*PMASS(L)*DES_RADIUS(L)**2
         ENDIF
         ! skip non-exsitent and ghost component spheres, for drag force, mpi all_reduce should be fine
         glueForce(gid,:) = glueForce(gid,:) + drag_fc(L,:)
      ENDDO

      if(numPEs > 1) THEN
         CALL GLOBAL_ALL_SUM(glueInertia)
         CALL GLOBAL_ALL_SUM(glueForce)
         CALL GLOBAL_ALL_SUM(glueTorque)
      endif

      IF(mype == pe_io) THEN
!step2, update velocity & pos of glued sphere
         DO gid = 1, NGluedParticles
         IF(PARTICLE_STATE_GSP(gid)== NONEXISTENT) CYCLE
            IF(glueMass(gid) .ne. 0) THEN
               if(PARTICLE_STATE_GSP(gid) == NORMAL_PARTICLE) then
                  glueAcc(gid,:) = glueForce(gid,:)/glueMass(gid) + GRAV(:)
                  IF (isNaN(glueAcc(gid,1))) &
                     WRITE(*,"(A,I4,7ES15.8)") 'gidgvel=',gid,glueAcc(gid,:),glueForce(gid,:),glueMass(gid)
                  glueVel(gid,:) = glueVel(gid,:) + glueAcc(gid,:)*dtsolid/2
               endif
            ENDIF
         ENDDO

         !step3, update rotation
         DO gid = 1, NGluedParticles
            IF(PARTICLE_STATE_GSP(gid)== NONEXISTENT) CYCLE
            angmom = glueAngMom(gid,:)
            quat = glueQuat(gid,:)
            inertia = glueInertia(gid,:)

            CALL gsp_quat_to_mat(quat,rot)
            torloc(:) = glueTorque(gid,:)
            CALL gsp_transpose_matvec(rot,torloc,tbody)
            CALL gsp_transpose_matvec(rot,angmom,mbody)
            CALL gsp_quatvec(quat,mbody,conjqm)
            conjqm = 2*conjqm

            CALL gsp_quatvec(quat,tbody,fquat)
            conjqm = conjqm + dtsolid*fquat

            CALL gsp_invquatvec(quat,conjqm,mbody)
            mbody = mbody*0.5
            local_omg = mbody/inertia
            CALL gsp_matvec(rot,mbody,angmom)
            CALL gsp_matvec(rot,local_omg,omega)

            glueAngMom(gid,:) = angmom
            glueOmg(gid,:) = omega
         ENDDO
      ENDIF !(mype == pe_io)

      if(numPEs > 1) THEN
         CALL BCAST(gluePos,PE_IO)
         CALL BCAST(glueVel,PE_IO)

         CALL BCAST(glueAngMom, PE_IO)
         CALL BCAST(glueOmg, PE_IO)
      ENDIF

      DO L = 1, MAX_PIP
         IF(IS_NONEXISTENT(L)) CYCLE
         IF(IS_ANY_GHOST(L)) CYCLE
         gid = gid_list(L)
         IF(gid .lt. 1) CYCLE
         DES_VEL_NEW(L,:) = glueVel(gid,:) ! + CROSS(glueOmg,dist2centerVEC) --> this will give total velocity instead of translation velocity
         OMEGA_NEW(L,:) = glueOmg(gid,:) ! used to write output of angular velocity of glued particle
      ENDDO
      ! if(mype == pe_io) then
      !    ! force report
      !    OPEN(UNIT=2856, FILE='output1/gp_force.txt', STATUS='UNKNOWN', ACTION='READWRITE', POSITION='APPEND')
      !       DO gid = 1, NGluedParticles
      !          write(2856,*) "glueForce:",gid,glueForce(gid,1),glueForce(gid,2),glueForce(gid,3)
      !       ENDDO
      !       write(2856,*) "---------------------------------------------->"
      !    CLOSE(UNIT=2856)
      RETURN

      END SUBROUTINE  GSP_CFNEWVALUES_vfull

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: GSP_CFNEWVALUES_vhalfxfull_implicit                    !
!  Author     : Hang Zhou                              Date: Apr-23-24 !
!  Revised    : Renjie Ke                              Date: Oct-31-24 !
!  Purpose: Verlet algorithm, update velocity by half time step        !
!           then update the position by a full step                    !
!           update on glued particle itself                            !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

   SUBROUTINE GSP_CFNEWVALUES_vhalfxfull_implicit

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE constant
      USE des_bc
      USE discretelement
      USE fldvar
      USE functions
      USE mfix_pic
      USE mpi_utility
      USE parallel
      USE run
      USE param
      USE param1
      USE physprop
      use mpi_funs_des, only: DES_PAR_EXCHANGE
      USE desgrid, only: desgrid_pic

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER ::L, pp, mm
      DOUBLE PRECISION :: quat(4), exone(3), eyone(3), ezone(3), torloc(3), tbody(3), mbody(3),rot(3,3),local_omg(3)
      DOUBLE PRECISION :: conjqm(4), fquat(4), inertia(3), angmom(3), omega(3), tow_part(3)
      DO L = 1, MAX_PIP
         IF(IS_NONEXISTENT(L) .OR. IS_ANY_GHOST(L)) CYCLE
         ! Only update velocity if it is a normal particle
         IF(IS_NORMAL(L)) &
            DES_VEL_NEW(L,:) = DES_VEL_NEW(L,:) + (FC(L,:)/PMASS(L)+GRAV(:))*DTSOLID/2.0d0
         DES_POS_NEW(L,:) = DES_POS_NEW(L,:) + DES_VEL_NEW(L,:) * DTSOLID
         FC(L, :) = ZERO

         angmom = glueAngMom(L,:)
         quat   = glueQuat(L,:)
         inertia= 1.0d0/OMOI3(L,:)
         tow_part = TOW(L,:)

         CALL gsp_quat_to_mat(quat,rot)
         CALL gsp_transpose_matvec(rot,angmom,mbody)
         CALL gsp_quatvec(quat,mbody,conjqm)
         conjqm = 2 * conjqm

         CALL gsp_transpose_matvec(rot,tow_part,tbody)
         CALL gsp_quatvec(quat,tbody,fquat)
         conjqm = conjqm + dtsolid*fquat

         CALL gsp_no_squish_rotate(3,conjqm,quat,inertia,dtsolid/2)
         CALL gsp_no_squish_rotate(2,conjqm,quat,inertia,dtsolid/2)
         CALL gsp_no_squish_rotate(1,conjqm,quat,inertia,dtsolid)
         CALL gsp_no_squish_rotate(2,conjqm,quat,inertia,dtsolid/2)
         CALL gsp_no_squish_rotate(3,conjqm,quat,inertia,dtsolid/2)

         CALL gsp_qnormalize(quat)
         CALL gsp_quat_to_mat(quat,rot)
         CALL gsp_invquatvec(quat,conjqm,mbody)
         mbody = mbody * 0.5
         local_omg = mbody / inertia

         CALL gsp_quat_to_exyz(quat,exone,eyone,ezone)
         glueEX(L,:) = exone(:)
         glueEY(L,:) = eyone(:)
         glueEZ(L,:) = ezone(:)

         CALL gsp_matvec(rot,mbody,angmom)
         CALL gsp_matvec(rot,local_omg,omega)

         glueAngMom(L,:) = angmom
         glueQuat(L,:) = quat
         omega_new(L,1:3) = omega
         IF(PARTICLE_ORIENTATION) orientation(L,1:3) = eyone(1:3)

         TOW(L,:) = ZERO
      ENDDO

      IF(numPEs > 1) THEN
         ! rebin to desgrid is needed before calling des_par_exchange
         CALL DESGRID_PIC(.TRUE.)
         CALL DES_PAR_EXCHANGE
      ENDIF

      RETURN

   END SUBROUTINE GSP_CFNEWVALUES_vhalfxfull_implicit

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: GSP_CFNEWVALUES_vfull_implicit                         !
!  Author     : Hang Zhou                              Date: Apr-23-24 !
!  Revised    : Renjie Ke                              Date: Oct-31-24 !
!                                                                      !
!  Purpose: Verlet algorithm,                                          !
!           update another half step velocity                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

   SUBROUTINE GSP_CFNEWVALUES_vfull_implicit

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE constant
      USE des_bc
      USE discretelement
      USE fldvar
      USE mfix_pic
      USE mpi_utility
      USE parallel
      USE run
      USE param
      USE param1
      USE physprop
      USE functions, only: is_nonexistent, is_ghost, is_entering_ghost, is_exiting_ghost, is_normal, is_any_ghost

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER ::L, tempPhaseNumber
      DOUBLE PRECISION :: quat(4),exone(3),eyone(3),ezone(3),torloc(3),tbody(3),mbody(3),rot(3,3),local_omg(3)
      DOUBLE PRECISION :: conjqm(4),fquat(4),inertia(3),angmom(3),omega(3),dist2, tow_part(3)
      double precision, allocatable, dimension(:) :: tmpgluemass
      double precision, allocatable, dimension(:,:) :: tmpglueInertia, tmpglueForce, tmpglueTorque
      integer, save :: totalcnt = 0
!-----------------------------------------------

      DO L = 1, MAX_PIP
         IF(IS_NONEXISTENT(L) .or. is_any_ghost(L)) CYCLE
         !update velocity from t+1/2dt to d+dt
         IF(IS_NORMAL(L)) &
            DES_VEL_NEW(L,:) = DES_VEL_NEW(L,:) + (FC(L,:)/PMASS(L) + GRAV(:))*DTSOLID/2.0d0
         angmom = glueAngMom(L,:)
         quat = glueQuat(L,:)
         ! @renjieke 10-31-2024 if chemical flow, inertia requires the update
         inertia= 1.0d0/OMOI3(L,:)
         tow_part = TOW(L,:)

         CALL gsp_quat_to_mat(quat,rot)
         CALL gsp_transpose_matvec(rot,tow_part,tbody)
         CALL gsp_transpose_matvec(rot,angmom,mbody)
         CALL gsp_quatvec(quat,mbody,conjqm)
         conjqm = 2*conjqm

         CALL gsp_quatvec(quat,tbody,fquat)
         conjqm = conjqm + dtsolid*fquat

         CALL gsp_invquatvec(quat,conjqm,mbody)
         mbody = mbody*0.5
         local_omg = mbody/inertia
         CALL gsp_matvec(rot,mbody,angmom)
         CALL gsp_matvec(rot,local_omg,omega)

         glueAngMom(L,:) = angmom
         OMEGA_NEW(L,1:3)= omega
      ENDDO
      RETURN

      END SUBROUTINE  GSP_CFNEWVALUES_vfull_implicit

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES                                            C
!
!  Purpose: DES - Calculate the new values of particle velocity,
!           position, angular velocity etc
!
!                                                                      C
!  Comments: Implements Eqns 1, 2, 3, 4 & 5  from the following paper:
!    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical
!    simulation of plug glow of cohesionless particles in a
!    horizontal pipe", Powder technology, 71, 239-250, 1992
!
!  pradeep : changes for parallel processing
!          1. periodic boundaries might lie in different proc. so adjust
!             particle position for periodic removed
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CFNEWVALUES

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE constant
      USE des_bc
      USE discretelement
      USE fldvar
      USE mfix_pic
      USE mpi_utility
      USE parallel
      USE run
      USE param
      USE param1
      USE physprop
      use geometry, only: DO_K, NO_K

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: L
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.
      DOUBLE PRECISION :: OMEGA_MAG,OMEGA_UNIT(3),ROT_ANGLE

!-----------------------------------------------

! Apply prescribed motion to particles if needed
      CALL DEM_SET_RIGID_MOTION

! Adams-Bashforth defaults to Euler for the first time step.
      IF(FIRST_PASS .AND. INTG_ADAMS_BASHFORTH) THEN
         DO L =1, MAX_PIP
            IF(IS_NONEXISTENT(L)) CYCLE                       ! Only real particles
            IF(IS_ENTERING(L).or.IS_ENTERING_GHOST(L)) CYCLE  ! Only non-entering
            IF(IS_GHOST(L)) CYCLE                             ! Skip ghost particles
            IF(IS_RIGID_MOTION_GHOST(L)) CYCLE
            DES_ACC_OLD(L,:) = FC(L,:)/PMASS(L) + GRAV(:)
            ROT_ACC_OLD(L,:) = TOW(L,:)
         ENDDO
      ENDIF


!!$omp parallel default(none)                    &
!!$omp shared(MAX_PIP,INTG_EULER,INTG_ADAMS_BASHFORTH,fc,tow,do_nsearch,   &
!!$omp       omega_new,omega_old,pmass,grav,des_vel_new,des_pos_new,       &
!!$omp       des_vel_old,des_pos_old,dtsolid,omoi,des_acc_old,rot_acc_old, &
!!$omp       ppos,neighbor_search_rad_ratio,des_radius,DO_OLD, iGlobal_ID, &
!!$omp       residence_time, &
!!$omp       particle_orientation,orientation,particle_state) &
!!$omp private(l,rot_angle,omega_mag,omega_unit,aabb)

! If a particle is classified as new, then forces are ignored.
! Classification from new to existing is performed in routine
! des_check_new_particle.f

! Advance particle position, velocity
! first-order method
      IF (INTG_EULER) THEN
!!$omp sections
!!$omp section


         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,1) =   &
            DES_VEL_NEW(:,1) + DTSOLID*(FC(:,1)/PMASS(:) + GRAV(1))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,1) =       &
            DES_POS_NEW(:,1) + DES_VEL_NEW(:,1)*DTSOLID
         FC(:,1) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,2) =   &
            DES_VEL_NEW(:,2) + DTSOLID*(FC(:,2)/PMASS(:) + GRAV(2))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,2) =       &
            DES_POS_NEW(:,2) + DES_VEL_NEW(:,2)*DTSOLID
         FC(:,2) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,3) =   &
            DES_VEL_NEW(:,3) + DTSOLID*(FC(:,3)/PMASS(:) + GRAV(3))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,3) =       &
            DES_POS_NEW(:,3) + DES_VEL_NEW(:,3)*DTSOLID
         FC(:,3) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) OMEGA_NEW(:,1) =     &
            OMEGA_NEW(:,1) + TOW(:,1)*OMOI(:)*DTSOLID
         TOW(:,1) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE)OMEGA_NEW(:,2) =      &
            OMEGA_NEW(:,2) + TOW(:,2)*OMOI(:)*DTSOLID
         TOW(:,2) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE)OMEGA_NEW(:,3) =      &
            OMEGA_NEW(:,3) + TOW(:,3)*OMOI(:)*DTSOLID
         TOW(:,3) = ZERO
!!$omp end sections

! Second-order Adams-Bashforth/Trapezoidal scheme
      ELSEIF (INTG_ADAMS_BASHFORTH) THEN

!!$omp sections
!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            FC(:MAX_PIP,1) = FC(:MAX_PIP,1)/PMASS(:MAX_PIP) + GRAV(1)
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            DES_VEL_NEW(:MAX_PIP,1) = DES_VEL_OLD(:MAX_PIP,1) + 0.5d0* &
               (3.d0*FC(:MAX_PIP,1)-DES_ACC_OLD(:MAX_PIP,1) )*DTSOLID
            DES_ACC_OLD(:MAX_PIP,1) = FC(:MAX_PIP,1)

         WHERE(PARTICLE_STATE < NORMAL_GHOST)                          &
            DES_POS_NEW(:MAX_PIP,1) = DES_POS_OLD(:MAX_PIP,1) + 0.5d0* &
               (DES_VEL_OLD(:MAX_PIP,1)+DES_VEL_NEW(:MAX_PIP,1))*DTSOLID
         FC(:MAX_PIP,1) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            FC(:MAX_PIP,2) = FC(:MAX_PIP,2)/PMASS(:MAX_PIP) + GRAV(2)
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            DES_VEL_NEW(:MAX_PIP,2) = DES_VEL_OLD(:MAX_PIP,2) + 0.5d0* &
               (3.d0*FC(:MAX_PIP,2)-DES_ACC_OLD(:MAX_PIP,2) )*DTSOLID
            DES_ACC_OLD(:MAX_PIP,2) = FC(:MAX_PIP,2)

         WHERE(PARTICLE_STATE < NORMAL_GHOST)                          &
            DES_POS_NEW(:MAX_PIP,2) = DES_POS_OLD(:MAX_PIP,2) + 0.5d0* &
               (DES_VEL_OLD(:MAX_PIP,2)+DES_VEL_NEW(:MAX_PIP,2))*DTSOLID
         FC(:MAX_PIP,2) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            FC(:MAX_PIP,3) = FC(:MAX_PIP,3)/PMASS(:MAX_PIP) + GRAV(3)
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            DES_VEL_NEW(:MAX_PIP,3) = DES_VEL_OLD(:MAX_PIP,3) + 0.5d0* &
                 (3.d0*FC(:MAX_PIP,3)-DES_ACC_OLD(:MAX_PIP,3) )*DTSOLID
            DES_ACC_OLD(:MAX_PIP,3) = FC(:MAX_PIP,3)

         WHERE(PARTICLE_STATE < NORMAL_GHOST)                          &
            DES_POS_NEW(:MAX_PIP,3) = DES_POS_OLD(:MAX_PIP,3) + 0.5d0* &
               (DES_VEL_OLD(:MAX_PIP,3)+DES_VEL_NEW(:MAX_PIP,3))*DTSOLID
         FC(:MAX_PIP,3) = ZERO

!!$omp section
        WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE) &
           OMEGA_NEW(:MAX_PIP,1) = OMEGA_OLD(:MAX_PIP,1) + 0.5d0*     &
               (3.d0*TOW(:MAX_PIP,1)*OMOI(:MAX_PIP) -                  &
               ROT_ACC_OLD(:MAX_PIP,1))*DTSOLID
        WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE) &
            ROT_ACC_OLD(:MAX_PIP,1) = TOW(:MAX_PIP,1)*OMOI(:MAX_PIP)
         TOW(:MAX_PIP,1) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            OMEGA_NEW(:MAX_PIP,2) = OMEGA_OLD(:MAX_PIP,2) + 0.5d0*     &
                 (3.d0*TOW(:MAX_PIP,2)*OMOI(:MAX_PIP)-                 &
                 ROT_ACC_OLD(:MAX_PIP,2) )*DTSOLID
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            ROT_ACC_OLD(:MAX_PIP,2) = TOW(:MAX_PIP,2)*OMOI(:MAX_PIP)
         TOW(:MAX_PIP,2) = ZERO

!!$omp section
        WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            OMEGA_NEW(:MAX_PIP,3) = OMEGA_OLD(:MAX_PIP,3) + 0.5d0*     &
               (3.d0*TOW(:MAX_PIP,3)*OMOI(:MAX_PIP)-                   &
               ROT_ACC_OLD(:MAX_PIP,3) )*DTSOLID
        WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            ROT_ACC_OLD(:MAX_PIP,3) = TOW(:MAX_PIP,3)*OMOI(:MAX_PIP)
         TOW(:MAX_PIP,3) = ZERO

!!$omp end sections
      ENDIF

! Update particle orientation - Always first order
! When omega is non-zero, compute the rotation angle, and apply the
! Rodrigues' rotation formula

      IF(PARTICLE_ORIENTATION) THEN
         DO L = 1, MAX_PIP
            OMEGA_MAG = OMEGA_NEW(L,1)**2 + OMEGA_NEW(L,2)**2 + OMEGA_NEW(L,3)**2

            IF(OMEGA_MAG>ZERO) THEN
               OMEGA_MAG=DSQRT(OMEGA_MAG)
               OMEGA_UNIT(:) = OMEGA_NEW(L,:)/OMEGA_MAG
               ROT_ANGLE = OMEGA_MAG * DTSOLID

               ORIENTATION(L,:) = ORIENTATION(L,:)*DCOS(ROT_ANGLE) + &
                  CROSS(OMEGA_UNIT,ORIENTATION(L,:))*DSIN(ROT_ANGLE) + &
                  OMEGA_UNIT(:)*DOT_PRODUCT(OMEGA_UNIT,ORIENTATION(L,:))*&
                  (ONE-DCOS(ROT_ANGLE))
            ENDIF
         ENDDO
      ENDIF

! Residence time
!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) RESIDENCE_TIME(:) = RESIDENCE_TIME(:) + DTSOLID
!!$omp end sections

! Check if the particle has moved a distance greater than or equal to
! its radius since the last time a neighbor search was called. if so,
! make sure that neighbor is called in des_time_march
      IF(.NOT.DO_NSEARCH) THEN
!!$omp do reduction (.or.:do_nsearch)
         DO L = 1, MAX_PIP
            DO_NSEARCH = DO_NSEARCH .or. &
               (DES_POS_NEW(L,1) - PPOS(L,1))**2+              &
               (DES_POS_NEW(L,2) - PPOS(L,2))**2+              &
               (DES_POS_NEW(L,3) - PPOS(L,3))**2  >=           &
               (NEIGHBOR_SEARCH_RAD_RATIO*DES_RADIUS(L))**2
         ENDDO
      ENDIF

!!$omp end parallel

      FIRST_PASS = .FALSE.

      RETURN

      contains

      include 'functions.inc'

   END SUBROUTINE CFNEWVALUES

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SuperDEM_CFNEWVALUES                                   C
!                                                                      C
!  Author: Xi Gao                                 Date: 25-Feb-2019    C
!                                                                      C
!  Purpose: SuperDEM - Calculate the new values of particle velocity,  C
!           position, angular velocity etc                             C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE SuperDEM_CFNEWVALUES

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE constant
      USE des_bc
      USE discretelement
      USE fldvar
      USE mfix_pic
      USE mpi_utility
      USE parallel
      USE run
      USE param
      USE param1
      USE physprop
      use geometry, only: DO_K, NO_K

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: L, ixg
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.
      DOUBLE PRECISION :: OMEGA_MAG,OMEGA_UNIT(3),ROT_ANGLE !! nb this is not the keyword omega_unit
      DOUBLE PRECISION :: Qc(4),Qc_new(4),DW(3),DW_local(3)
      DOUBLE PRECISION :: OMEGA_NEW_global(3),OMEGA_NEW_LOCAL(3)
      DOUBLE PRECISION :: tow_global(3),tow_local(3)
      LOGICAL, PARAMETER :: ldebug = .false.

!-----------------------------------------------

! Adams-Bashforth defaults to Euler for the first time step.
      IF(FIRST_PASS .AND. INTG_ADAMS_BASHFORTH) THEN
         DO L =1, MAX_PIP
            IF(IS_NONEXISTENT(L)) CYCLE                       ! Only real particles
            IF(IS_ENTERING(L).or.IS_ENTERING_GHOST(L)) CYCLE  ! Only non-entering
            IF(IS_GHOST(L)) CYCLE                             ! Skip ghost particles
            IF(IS_RIGID_MOTION_GHOST(L)) CYCLE                ! Skip rigid motion ghost particles
            DES_ACC_OLD(L,:) = FC(L,:)/PMASS(L) + GRAV(:)     ! should this 'FC' and 'TOW'
            ROT_ACC_OLD(L,:) = TOW(L,:)                       ! have the dimension (3, PARTICLES)?
         ENDDO
      ENDIF

! If a particle is classified as new, then forces are ignored.
! Classification from new to existing is performed in routine
! des_check_new_particle.f

! Advance particle position, velocity

         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,1) =   &
            DES_VEL_NEW(:,1) + DTSOLID*(FC(:,1)/PMASS(:) + GRAV(1))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,1) =       &
            DES_POS_NEW(:,1) + DES_VEL_NEW(:,1)*DTSOLID
         FC(:,1) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,2) =   &
            DES_VEL_NEW(:,2) + DTSOLID*(FC(:,2)/PMASS(:) + GRAV(2))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,2) =       &
            DES_POS_NEW(:,2) + DES_VEL_NEW(:,2)*DTSOLID
         FC(:,2) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,3) =   &
            DES_VEL_NEW(:,3) + DTSOLID*(FC(:,3)/PMASS(:) + GRAV(3))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,3) =       &
            DES_POS_NEW(:,3) + DES_VEL_NEW(:,3)*DTSOLID
         FC(:,3) = ZERO


! Update of angular velocity
      IF (SuperDEM) THEN
        do L = 1, MAX_PIP
! Normal_particle
           if (PARTICLE_STATE(L) == normal_particle) then
           Qc(1)= super_q(L,1)
           Qc(2)= super_q(L,2)
           Qc(3)= super_q(L,3)
           Qc(4)= super_q(L,4)
           OMEGA_NEW_global(1)=OMEGA_NEW(L,1)
           OMEGA_NEW_global(2)=OMEGA_NEW(L,2)
           OMEGA_NEW_global(3)=OMEGA_NEW(L,3)
! Rotate angular velocity from global to local
           ixg=1
           CALL QROTATE(Qc,OMEGA_NEW_global,OMEGA_NEW_local,Ixg)
           TOW_global(1)=TOW(L,1)
           TOW_global(2)=TOW(L,2)
           TOW_global(3)=TOW(L,3)
           ixg=1
           CALL QROTATE(Qc,TOW_global,TOW_local,Ixg)

            OMEGA_NEW_LOCAL(1) = OMEGA_NEW_LOCAL(1) + &
                        ((1.0/OMOI3(L,2)-1.0/OMOI3(L,3))*OMEGA_NEW_LOCAL(2)*&
                        OMEGA_NEW_LOCAL(3)+TOW_LOCAL(1))*DTSOLID*OMOI3(L,1)

            OMEGA_NEW_LOCAL(2) = OMEGA_NEW_LOCAL(2) + &
                        ((1.0/OMOI3(L,3)-1.0/OMOI3(L,1))*OMEGA_NEW_LOCAL(3)*&
                        OMEGA_NEW_LOCAL(1)+TOW_LOCAL(2))*DTSOLID*OMOI3(L,2)

            OMEGA_NEW_LOCAL(3) = OMEGA_NEW_LOCAL(3) + &
                        ((1.0/OMOI3(L,1)-1.0/OMOI3(L,2))*OMEGA_NEW_LOCAL(1)*&
                        OMEGA_NEW_LOCAL(2)+TOW_LOCAL(3))*DTSOLID*OMOI3(L,3)
! Rotate angular velocity from local to global
            IXG=2
            CALL QROTATE(Qc,OMEGA_NEW_global,OMEGA_NEW_local,Ixg)
            OMEGA_NEW(L,1)= OMEGA_NEW_global(1)
            OMEGA_NEW(L,2)= OMEGA_NEW_global(2)
            OMEGA_NEW(L,3)= OMEGA_NEW_global(3)
            tow(L,1)=0
            tow(L,2)=0
            tow(L,3)=0
            endif
        enddo
      ENDIF !SuperDEM

!Update particle quaternion
     if (SuperDEM) then
        do L = 1, MAX_PIP
           if (PARTICLE_STATE(L) == normal_particle) then
              Qc(1)= super_q(L,1)
              Qc(2)= super_q(L,2)
              Qc(3)= super_q(L,3)
              Qc(4)= super_q(L,4)
! Incremental of angular velocity in the global frame
              DW(1)= OMEGA_NEW(L,1)*DTSOLID
              DW(2)= OMEGA_NEW(L,2)*DTSOLID
              DW(3)= OMEGA_NEW(L,3)*DTSOLID
! Apply the incremental rotation to the orientation quaternion
              ixg=1
              CALL QROTATE(Qc,DW,DW_local,Ixg)
              CALL QINCROTATE2(Qc,Qc_new,Dw_local)
              super_q(L,1)= Qc_new(1)
              super_q(L,2)= Qc_new(2)
              super_q(L,3)= Qc_new(3)
              super_q(L,4)= Qc_new(4)

              IF(PARTICLE_ORIENTATION) THEN
                 ! Convert quaternion to orientation vector
                 ! Input: quaternion Qc_new
                 !        initial orientation vector INIT_ORIENTATION
                 ! Output: orientation vector (global frame) ORIENTATION
                 ! Note the last argument is 2 to apply the rotation on the
                 ! initial orientation
                 CALL QROTATE(Qc_new,ORIENTATION(L,:),INIT_ORIENTATION,2)
              ENDIF
           endif
        ENDDO
     endif !SuperDEM

! Residence time
!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) RESIDENCE_TIME(:) = RESIDENCE_TIME(:) + DTSOLID
!!$omp end sections

! Check if the particle has moved a distance greater than or equal to
! its radius since the last time a neighbor search was called. if so,
! make sure that neighbor is called in des_time_march
      IF(.NOT.DO_NSEARCH) THEN
!!$omp do reduction (.or.:do_nsearch)
         DO L = 1, MAX_PIP
            DO_NSEARCH = DO_NSEARCH .or. &
               (DES_POS_NEW(L,1) - PPOS(L,1))**2+              &
               (DES_POS_NEW(L,2) - PPOS(L,2))**2+              &
               (DES_POS_NEW(L,3) - PPOS(L,3))**2  >=           &
               (NEIGHBOR_SEARCH_RAD_RATIO*DES_RADIUS(L))**2
         ENDDO
      ENDIF

      FIRST_PASS = .FALSE.

      RETURN

      contains

      include 'functions.inc'

   END SUBROUTINE SuperDEM_CFNEWVALUES

   SUBROUTINE DEM_SET_RIGID_MOTION
   USE discretelement
   USE functions
   USE keyframe_mod, only:omega_conv_fact
   use get_stl_data_mod

   IMPLICIT NONE
   INTEGER :: L,M
! Line joining center of rotation to particle center
   double precision, dimension(dimn) :: OP, OP_NEW
! Angular velocity (rad/sec)
   double precision, dimension(3) :: omega

! Rotation matrices
   DOUBLE PRECISION :: Theta_deg(3)
   DOUBLE PRECISION, DIMENSION(3,3) :: Rx, Ry, Rz
   DOUBLE PRECISION, DIMENSION(DIM_M,3,3) :: R

   DO M=1,DIM_M
! OMEGA_CONV_FACT converts the angle to radians.
! The routines building the rotation matrices needs degrees.
      Theta_deg(:) = DES_RIGID_MOTION_OMEGA(M,:) * DTSOLID * OMEGA_CONV_FACT * 180.0/PI
      CALL BUILD_X_ROTATION_MATRIX(Theta_deg(1), Rx)
      CALL BUILD_Y_ROTATION_MATRIX(Theta_deg(2), Ry)
      CALL BUILD_Z_ROTATION_MATRIX(Theta_deg(3), Rz)
      R(M,:,:) = MATMUL(Rz,MATMUL(Ry,Rx)) ! incremental rotation matrix
   ENDDO

   DO L =1, MAX_PIP
      IF(.NOT.IS_RIGID_MOTION(L)) CYCLE                       ! Only real particles
      M = PIJK(L,5)

! It is better to apply the rotation and compute the velocity from the
! change in position rather than compute the velocity from omega because this
! can lead to large errors if velocity is large (large omega or large time
! step), where particles do not stay along a circle.
      OP = DES_POS_NEW(L,:) - DES_RIGID_MOTION_ROT_CENTER(M, :)
      OP_NEW = MATMUL(R(M,:,:),OP(:))
      DES_VEL_NEW(L,:) = des_rigid_motion_vel(M,:) + (OP_NEW(:)- OP(:))/DTSOLID

      omega = des_rigid_motion_omega(m, :) * omega_conv_fact
      OMEGA_NEW(L,:) = omega(:)

   ENDDO
   RETURN
   END SUBROUTINE DEM_SET_RIGID_MOTION

END MODULE CFNEWVALUES_MOD
