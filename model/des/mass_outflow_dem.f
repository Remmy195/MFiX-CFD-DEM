MODULE MASS_OUTFLOW_DEM_MOD
CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MASS_OUTFLOW_DEM                                        !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose:  This routine fills in the necessary information for new   !
!  particles entering the system.                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MASS_OUTFLOW_DEM(FORCE_NSEARCH)

      use compar, only: pe_io, mype, numPEs
      use derived_types, only: dg_pic
      use discretelement
      use des_bc
      use bc
      use functions
      use param1, only: zero
      use gsp_math_mod, only: SET_MARKED_GSP, SET_MARKED_GSP_STATE, findIndex_1i

      use mpi_utility, only: GLOBAL_ALL_OR

      implicit none

      LOGICAL, INTENT(INOUT) :: FORCE_NSEARCH

      INTEGER :: IJK
      INTEGER :: LC, LP, NP, M, gid, dummy_index
      INTEGER :: BCV, BCV_I, IDX

      DOUBLE PRECISION :: SGN
      DOUBLE PRECISION :: DIST

      LOGICAL :: FREEZE_VEL
      DOUBLE PRECISION :: FREEZE(3)

! distanc used to check if particles are exiting:
! for gsp_explicit, it is the bounding diameter of the gluedSphere
! for all other models, it is the particle diameter
      DOUBLE PRECISION :: compareDist

      DO BCV_I = 1, DEM_BCMO
         BCV = DEM_BCMO_MAP(BCV_I)
         FREEZE_VEL = (BC_TYPE_ENUM(BCV) /= MASS_OUTFLOW)

         SELECT CASE (BC_PLANE(BCV))
         CASE('N'); FREEZE = (/0.0d0, 1.0d0, 0.0d0/); IDX=2; SGN=-1.0d0
         CASE('S'); FREEZE = (/0.0d0, 1.0d0, 0.0d0/); IDX=2; SGN= 1.0d0
         CASE('E'); FREEZE = (/1.0d0, 0.0d0, 0.0d0/); IDX=1; SGN=-1.0d0
         CASE('W'); FREEZE = (/1.0d0, 0.0d0, 0.0d0/); IDX=1; SGN= 1.0d0
         CASE('T'); FREEZE = (/0.0d0, 0.0d0, 1.0d0/); IDX=3; SGN=-1.0d0
         CASE('B'); FREEZE = (/0.0d0, 0.0d0, 1.0d0/); IDX=3; SGN= 1.0d0
         END SELECT

         DO LC=DEM_BCMO_IJKSTART(BCV_I), DEM_BCMO_IJKEND(BCV_I)
            IJK = DEM_BCMO_IJK(LC)
            DO LP= 1,DG_PIC(IJK)%ISIZE

               NP = DG_PIC(IJK)%P(LP) ! NP for sure is local
               IF(IS_NONEXISTENT(NP)) CYCLE
               IF(IS_ANY_GHOST(NP)) CYCLE
               IF(IS_ENTERING(NP)) CYCLE

               IF (gsp_explicit) then
                  ! find the gsp this NP belongs to
                  gid = gid_list(NP)
                  ! use the gluePos of this gsp to determine the state
                  SELECT CASE (BC_PLANE(BCV))
                  CASE('N'); DIST = gluePos(gid,2) - YN(BC_J_s(BCV))
                  CASE('S'); DIST = YN(BC_J_s(BCV)-1) - gluePos(gid,2)
                  CASE('E'); DIST = gluePos(gid,1) - XE(BC_I_w(BCV))
                  CASE('W'); DIST = XE(BC_I_w(BCV)-1) - gluePos(gid,1)
                  CASE('T'); DIST = gluePos(gid,3) - ZT(BC_K_b(BCV))
                  CASE('B'); DIST = ZT(BC_K_b(BCV)-1) - gluePos(gid,3)
                  END SELECT
                  compareDist = glueBounding(gid)/2.0D0
               ELSE
                  SELECT CASE (BC_PLANE(BCV))
                  CASE('S'); DIST = YN(BC_J_s(BCV)-1) - DES_POS_NEW(NP,2)
                  CASE('N'); DIST = DES_POS_NEW(NP,2) - YN(BC_J_s(BCV))
                  CASE('W'); DIST = XE(BC_I_w(BCV)-1) - DES_POS_NEW(NP,1)
                  CASE('E'); DIST = DES_POS_NEW(NP,1) - XE(BC_I_w(BCV))
                  CASE('B'); DIST = ZT(BC_K_b(BCV)-1) - DES_POS_NEW(NP,3)
                  CASE('T'); DIST = DES_POS_NEW(NP,3) - ZT(BC_K_b(BCV))
                  END SELECT
                  compareDist = DES_RADIUS(NP)
               ENDIF

! The particle is still inside the domain
               IF(DIST > compareDist) THEN
                  CALL SET_NORMAL(NP)
                  ! this gsp is still inside the domain
                  IF (gsp_explicit) THEN
                     CALL SET_NORMAL_GSP(gid)
                     IF(numPEs > 1) THEN
                        ! I dont want the same gid being pushed into set marked gsp multiple times
                        dummy_index = findIndex_1i(MARKED_GSP,gid)
                        if(dummy_index /= -1) then
                           marked_state(dummy_index) = particle_state_gsp(gid)
                        else
                           CALL SET_MARKED_GSP(gid)
                           CALL SET_MARKED_GSP_STATE(particle_state_gsp(gid))
                        endif
                     ENDIF
                  ENDIF

! Check if the particle is crossing over the outlet plane.
               ELSEIF(DIST > ZERO) THEN
! The velocity is 'frozen' normal to the outflow plane. This approach
! is strict because complex BCs (via STLs) can let particles pop through
! the wall along the outlet.
                  IF(FREEZE_VEL) THEN
! Only 'freeze' a particle's velocy if it has it moving out of the
! domain. Otherwise, particles flagged as exiting but moving away from
! the BC appear to moon-walk through the domain until it crashes.
                     IF(DES_VEL_NEW(NP,IDX)*SGN > 0.0d0) THEN
                        DES_VEL_NEW(NP,:) = DES_VEL_NEW(NP,:)*FREEZE(:)
! Set the flags for an exiting particle.
                        IF (IS_GHOST(NP)) THEN
                           CALL SET_EXITING_GHOST(NP)
                        ELSE
                           CALL SET_EXITING(NP)
                        ENDIF
                        IF (gsp_explicit) THEN
                           glueVel(gid_list(NP),:) = glueVel(gid_list(NP),:)*FREEZE(:)
                           CALL SET_EXITING_GSP(gid_list(NP))
                           IF(numPEs > 1) THEN
                              dummy_index = findIndex_1i(MARKED_GSP,gid)
                              if(dummy_index /= -1) then
                                 marked_state(dummy_index) = particle_state_gsp(gid)
                              else
                                 CALL SET_MARKED_GSP(gid)
                                 CALL SET_MARKED_GSP_STATE(particle_state_gsp(gid))
                              endif
                           ENDIF
                        ENDIF
                     ENDIF

! The user specified velocity is applied to the exiting particle. This
! only applies to mass outflows where the speed at which particles
! exit needs to be controlled.
                  ELSE
                     M = PIJK(NP,5)
                     DES_VEL_NEW(NP,1) = BC_U_s(BCV,M)
                     DES_VEL_NEW(NP,2) = BC_V_s(BCV,M)
                     DES_VEL_NEW(NP,3) = BC_W_s(BCV,M)
! Set the flags for an exiting particle.
                     IF (IS_GHOST(NP)) THEN
                        CALL SET_EXITING_GHOST(NP)
                     ELSE
                        CALL SET_EXITING(NP)
                     ENDIF
                     IF (gsp_explicit) THEN
                        glueVel(gid_list(NP), :) = DES_VEL_NEW(NP, :)
                        CALL SET_EXITING_GSP(gid_list(NP))
                        IF(numPEs > 1) THEN
                           dummy_index = findIndex_1i(MARKED_GSP,gid)
                           if(dummy_index /= -1) then
                              marked_state(dummy_index) = particle_state_gsp(gid)
                           else
                              CALL SET_MARKED_GSP(gid)
                              CALL SET_MARKED_GSP_STATE(particle_state_gsp(gid))
                           endif
                        ENDIF
                     ENDIF
                  ENDIF

! Ladies and gentlemen, the particle has left the building.
               ELSE
                  ! this gsp is crossing the boundary and call to delete the particle associate with it
                  IF (gsp_explicit) THEN
                     CALL DELETE_PARTICLE_GSP(NP)
                  ELSE
                     CALL DELETE_PARTICLE(NP)
                  ENDIF
                  FORCE_NSEARCH = .TRUE.
               ENDIF
            ENDDO
         ENDDO
      ENDDO
! Sync the search flag across all processes.
      CALL GLOBAL_ALL_OR(FORCE_NSEARCH)

      RETURN
      END SUBROUTINE MASS_OUTFLOW_DEM

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DELETE_PARTICLE                                         !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose:  This routine is used to check if a new particle has fully !
!  entered the domain.  If so, the flag classifying the particle as new!
!  is removed, allowing the particle to respond to contact forces from !
!  walls and other particles.                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DELETE_PARTICLE(NP)

      USE compar
      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE indices
      USE param1
      USE physprop
      USE functions

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NP

      INTEGER :: num_delete, gid, M, PP, id_delete
      INTEGER, ALLOCATABLE :: delete_list(:)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      IF(gsp_explicit) THEN
        M = PIJK(NP,5)
        gid = gid_list(NP)
        num_delete = 1
        allocate(delete_list(sphereNumber(M))); delete_list = 0
        if(gid == -1) return ! once a gsp got deleted, its gid reset to -1
        DO PP= max(1, NP-sphereNumber(M)+1),min(NP+sphereNumber(M)-1, size(gid_list))
          IF(gid_list(PP) .eq. gid) THEN
             delete_list(num_delete) = PP
             num_delete = num_delete + 1
          ENDIF
        ENDDO
        CALL SET_NONEXISTENT_GSP(gid)
      ELSE
         allocate(delete_list(1))
         delete_list(1) = NP
      ENDIF

      DO PP=1, size(delete_list)
         id_delete = delete_list(PP)
         IF(id_delete == 0) CYCLE
         IF(IS_NONEXISTENT(id_delete)) CYCLE
         iGLOBAL_ID(id_delete) = -1
         IF(gsp_explicit) gid_list(id_delete) = -1

         DES_POS_NEW(id_delete,:) = ZERO
         DES_VEL_NEW(id_delete,:) = ZERO
         OMEGA_NEW(id_delete,:) = ZERO

         IF(PARTICLE_ORIENTATION) ORIENTATION(id_delete,1:3) = INIT_ORIENTATION

         IF (DO_OLD) THEN
            DES_POS_OLD(id_delete,:) = ZERO
            DES_VEL_OLD(id_delete,:) = ZERO
            OMEGA_OLD(id_delete,:) = ZERO
         ENDIF

         DES_RADIUS(id_delete) = ZERO
         PMASS(id_delete) = HUGE(0.0)
         PVOL(id_delete) = HUGE(0.0)
         RO_Sol(id_delete) = ZERO
         OMOI(id_delete) = ZERO
!SuperDEM
         IF(SuperDEM)OMOI3(id_delete,:) = ZERO

         FC(id_delete,:) = ZERO
         TOW(id_delete,:) = ZERO

         PPOS(id_delete,:) = ZERO

         WALL_COLLISION_FACET_ID(:,id_delete) = -1
!gsp_explicit
!TBD, all glue variables related to PIP should reset to 0
!     should also reset variables related to NGluedParticles by using gid_list
         IF(gsp_explicit) THEN
            gp_sa(id_delete) = ZERO
            gp_neighsid(id_delete,:) = ZERO
            gp_neighsa(id_delete,:) = ZERO
            gp_kps(id_delete) = ZERO
            !gp_sdrag_mass(id_delete) = ZERO
            gp_squat(id_delete,:) = ZERO
            sc2gpc_vec(id_delete,:) = ZERO
         ENDIF

         PIP = PIP - 1

!JFD: If the particle was a ghost particle, we also need to decrement ighost_cnt,
!     otherwise this inconsistency propagates and failure occurred while writing
!     the DES restart file (The number of particle PIP-ighost_cnt becomes
!     negative). This occurred with non-uniform fluid grid, where the des grid and
!     the fluid grid are not aligned.
! @RENJIEKE 5-20-2025, change from is ghost to is any ghost for fixing gsp implicit restart issue.
         IF(IS_ANY_GHOST(id_delete)) ighost_cnt = ighost_cnt - 1

         CALL SET_NONEXISTENT(id_delete)
      ENDDO
      deallocate(delete_list)
      RETURN

   END SUBROUTINE DELETE_PARTICLE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DELETE_PARTICLE_GSP                                     !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!  revised: Renjie Ke                                 Date: 15-Mar-24  !
!                                                                      !
!  Purpose:  This routine is called in mass_outflow_dem for gsp model  !
!            only. It is almost the same as original DELETE_PARTICLE,  !
!            however, iglobal_id is not yet filled during auto-seeding !
!            for gsp. So we used different delete particle subroutine  !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DELETE_PARTICLE_GSP(NP)

      USE compar
      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE indices
      USE param1
      USE physprop
      USE functions
      use gsp_math_mod, only: SET_MARKED_GSP, SET_MARKED_GSP_STATE, findIndex_1i

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NP

      INTEGER :: num_delete, gid, M, PP, id_delete, dummy_index
      INTEGER, ALLOCATABLE :: delete_list(:)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      IF(gsp_explicit) THEN
        M = PIJK(NP,5)
        gid = gid_list(NP)
        num_delete = 1
        allocate(delete_list(sphereNumber(M))); delete_list = 0
        if(gid .lt. 1) return ! once a gsp got deleted, its gid reset to -1
        DO PP = 1, size(iglobal_id)
          if(iglobal_id(PP) .lt. 1) cycle
          IF(gid_list(PP) .eq. gid) THEN
            if(particle_state(PP) == 1 .or. particle_state(PP) == 2 .or. particle_state(PP) == 3) then
             delete_list(num_delete) = PP
             num_delete = num_delete + 1
            endif
          ENDIF
        ENDDO

        CALL SET_NONEXISTENT_GSP(gid)
        IF(numPEs > 1) THEN
         dummy_index = findIndex_1i(MARKED_GSP,gid)
         if(dummy_index /= -1) then
            marked_state(dummy_index) = particle_state_gsp(gid)
         else
            CALL SET_MARKED_GSP(gid)
            CALL SET_MARKED_GSP_STATE(particle_state_gsp(gid))
         endif
        ENDIF
      ELSE
         allocate(delete_list(1))
         delete_list(1) = NP
      ENDIF

      DO PP=1, size(delete_list)
         id_delete = delete_list(PP)
         ! a gsp can be half normal half ghost in current rank, only normal will be considered
         ! ghost leads to a zero in id_delete
         IF(id_delete == 0) CYCLE
	      IF(IS_NONEXISTENT(id_delete)) CYCLE
         iGLOBAL_ID(id_delete) = -1
         IF(gsp_explicit) gid_list(id_delete) = -1

         DES_POS_NEW(id_delete,:) = ZERO
         DES_VEL_NEW(id_delete,:) = ZERO
         OMEGA_NEW(id_delete,:) = ZERO

         IF(PARTICLE_ORIENTATION) ORIENTATION(id_delete,1:3) = INIT_ORIENTATION

         IF (DO_OLD) THEN
            DES_POS_OLD(id_delete,:) = ZERO
            DES_VEL_OLD(id_delete,:) = ZERO
            OMEGA_OLD(id_delete,:) = ZERO
         ENDIF

         DES_RADIUS(id_delete) = ZERO
         PMASS(id_delete) = HUGE(0.0)
         PVOL(id_delete) = HUGE(0.0)
         RO_Sol(id_delete) = ZERO
         OMOI(id_delete) = ZERO

         FC(id_delete,:) = ZERO
         TOW(id_delete,:) = ZERO

         PPOS(id_delete,:) = ZERO

         WALL_COLLISION_FACET_ID(:,id_delete) = -1
!gsp_explicit
!TBD, all glue variables related to PIP should reset to 0
!     should also reset variables related to NGluedParticles by using gid_list
         IF(gsp_explicit) THEN
            gp_sa(id_delete) = ZERO
            gp_neighsid(id_delete,:) = ZERO
            gp_neighsa(id_delete,:) = ZERO
            gp_kps(id_delete) = ZERO
            !gp_sdrag_mass(id_delete) = ZERO
            gp_squat(id_delete,:) = ZERO
            sc2gpc_vec(id_delete,:) = ZERO
         ENDIF

         PIP = PIP - 1

!JFD: If the particle was a ghost particle, we also need to decrement ighost_cnt,
!     otherwise this inconsistency propagates and failure occurred while writing
!     the DES restart file (The number of particle PIP-ighost_cnt becomes
!     negative). This occurred with non-uniform fluid grid, where the des grid and
!     the fluid grid are not aligned.
         IF(IS_GHOST(id_delete)) ighost_cnt = ighost_cnt - 1
         CALL SET_NONEXISTENT(id_delete)
      ENDDO

      deallocate(delete_list)
      RETURN

   END SUBROUTINE DELETE_PARTICLE_GSP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: UPDATE_EXIT_DES_VEL                                     !
!                                                                      !
!  Purpose:  This routine updates the velocity of particles that are   !
!  already identified as EXITING particles associated with             !
!  Non-MASS_OUTFLOW boundaries.                                        !
!  The variable FC utilized in the update has the contact force of ONLY!
!  NORMAL particles on each EXITING particle.                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE UPDATE_EXIT_DES_VEL()

      use derived_types, only: dg_pic
      use discretelement, only: DES_VEL_NEW, DES_VEL_OLD, DES_ACC_OLD
      use discretelement, only: DTSOLID, GRAV, FC, PMASS
      use discretelement, only: INTG_EULER, INTG_ADAMS_BASHFORTH
      use des_bc, only: DEM_BCMO, DEM_BCMO_MAP
      use des_bc, only: DEM_BCMO_IJK, DEM_BCMO_IJKSTART, DEM_BCMO_IJKEND
      use bc, only: BC_TYPE_ENUM, MASS_OUTFLOW, BC_PLANE
      use functions, only: IS_EXITING

      implicit none

      INTEGER :: IJK
      INTEGER :: LC, LP, NP
      INTEGER :: BCV, BCV_I
      DOUBLE PRECISION :: FREEZE(3)
      LOGICAL, SAVE :: FIRST_PASS = .TRUE. ! Used with Adams-Bashforth

      DO BCV_I = 1, DEM_BCMO

         BCV = DEM_BCMO_MAP(BCV_I)

         ! Only Non-MASS_OUTFLOW boundaries are required
         IF(BC_TYPE_ENUM(BCV) == MASS_OUTFLOW) cycle

         SELECT CASE (BC_PLANE(BCV))
         CASE('N'); FREEZE = (/0.0d0, 1.0d0, 0.0d0/)
         CASE('S'); FREEZE = (/0.0d0, 1.0d0, 0.0d0/)
         CASE('E'); FREEZE = (/1.0d0, 0.0d0, 0.0d0/)
         CASE('W'); FREEZE = (/1.0d0, 0.0d0, 0.0d0/)
         CASE('T'); FREEZE = (/0.0d0, 0.0d0, 1.0d0/)
         CASE('B'); FREEZE = (/0.0d0, 0.0d0, 1.0d0/)
         END SELECT

         DO LC=DEM_BCMO_IJKSTART(BCV_I), DEM_BCMO_IJKEND(BCV_I)
            IJK = DEM_BCMO_IJK(LC)
            DO LP= 1,DG_PIC(IJK)%ISIZE

               NP = DG_PIC(IJK)%P(LP)
               ! Only EXITING particles are required
               IF(.NOT. IS_EXITING(NP)) CYCLE
               ! Update EXITING particle velocity based on scheme
               ! Written following CFNEWVALUES
               ! Adams-Bashforth defaults to Euler
               ! for the first time step.
               IF(FIRST_PASS .AND. INTG_ADAMS_BASHFORTH) THEN
                  DES_ACC_OLD(NP,:) = FC(NP,:)/PMASS(NP) + GRAV(:)
               ENDIF

               IF(INTG_EULER) THEN
                  DES_VEL_NEW(NP,:) = DES_VEL_NEW(NP,:) &
                          + DTSOLID*(FC(NP,:)/PMASS(NP) + GRAV(:))
               ELSEIF(INTG_ADAMS_BASHFORTH) THEN
                  FC(NP,:) = FC(NP,:)/PMASS(NP) + GRAV(:)
                  DES_VEL_NEW(NP,:) = DES_VEL_OLD(NP,:) + 0.5d0* &
                     (3.0d0*FC(NP,:)-DES_ACC_OLD(NP,:))*DTSOLID
                  DES_ACC_OLD(NP,:) = FC(NP,:)
               ENDIF

               ! Zero out the non-normal component for the particle
               DES_VEL_NEW(NP,:) = DES_VEL_NEW(NP,:)*FREEZE(:)

            ENDDO
         ENDDO
      ENDDO

      FIRST_PASS = .FALSE.

      RETURN
      END SUBROUTINE UPDATE_EXIT_DES_VEL

END MODULE MASS_OUTFLOW_DEM_MOD
