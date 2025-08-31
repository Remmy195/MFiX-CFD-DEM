#include "error.inc"
MODULE DES_REACTION_MODEL_MOD
CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_REACTION_MODEL                                     !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE DES_REACTION_MODEL

      USE constant, only: pi
      USE compar, only: PE_IO, myPE, numPEs
      USE des_rxns, only: des_r_s, des_x_s
      USE discretelement, only: des_radius, pvol, ro_sol, pijk, pmass, omoi,omoi3, intg_euler, dtsolid
      USE discretelement, only: max_pip, des_explicitly_coupled, particle_state, normal_particle
      USE error_manager
      USE discretelement, only:  glueDiameter, glueinertia, gluevolume, sc2gpc_vec, gp_neighsa, gp_neighsid, iglobal_id
      USE discretelement, only: gp_sa, gid_list, des_pos_new
      USE functions, only: is_normal
      USE param, only: dimension_n_s
      USE param1, only: zero
      USE run, only: ANY_SOLIDS_SPECIES_EQ, SPECIES_EQ
      USE run, only: DT
      USE run, only: SOLVE_ROs
      USE discretelement, only: CGDEM, DES_CGP_RPR, DES_CGP_STW
      USE discretelement, only: super_r, super_mn,SuperDEM
      USE discretelement, only: GluedSphereDEM, gsp_explicit, gsp_implicit, glueMass
      USE sq_properties_mod
      USE mpi_utility, only: allgather_1i, global_sum, global_all_sum, global_all_max, bcast
      IMPLICIT NONE

! Loop counter
      INTEGER :: NN, LL
! total rate of consumption/production of species (g/sec)
      DOUBLE PRECISION  :: SUM_DES_Rs(1:MAX_PIP)

      DOUBLE PRECISION :: PIx4o3, MOI, dist2
      DOUBLE PRECISION :: o3 = 1.0d0/3.0d0

      DOUBLE PRECISION :: lDT, lOoDT

      DOUBLE PRECISION :: axi(3),m,n,IXX, IYY, IZZ, SVOLUME, scale_factor
      ! GSP Variables
      DOUBLE PRECISION :: rad_p_old
      INTEGER :: nbIDX, lcurpos, new_nbidx, locpp, J, gid
      DOUBLE PRECISION, allocatable, dimension(:) :: new_total_mass
      LOGICAL :: is_constant_density = .False.
      LOGICAL :: TMP_NSEARCH
      DOUBLE PRECISION :: cur_nsa, next_nsa
! Logical for Adams-Bashfort integration.
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
! summation of mass fractions of species in each partilce
      DOUBLE PRECISION :: sum_X_s(MAX_PIP)

!---------------------------------------------------------------------//
      IF(.NOT.ANY_SOLIDS_SPECIES_EQ) RETURN
      PIx4o3 = Pi*4.0d0/3.0d0

      lDT = merge(DT, DTSOLID, DES_EXPLICITLY_COUPLED)
      lOoDT = -1.0d0/lDT

! Bound the amount of mass loss. The definition of lOodT as a negative
! value allows the max function evaluation here
      WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE)            &
         SUM_DES_Rs(:MAX_PIP) = sum(DES_R_s(:MAX_PIP,:),DIM=2)

      DO LL=1,MAX_PIP
         IF (PARTICLE_STATE(LL) .ne. NORMAL_PARTICLE) CYCLE
         DO NN=1, DIMENSION_N_S
            IF(DES_R_s(LL,NN) .lt. DES_X_s(LL,NN)*lOoDT) THEN
               DES_R_s(LL,NN) = DES_X_s(LL,NN)*lOoDT
               WRITE(ERR_MSG,1000) LL, NN
               CALL LOG_WARNING()
            ENDIF
         ENDDO
      ENDDO

 1000 FORMAT('Warning 1000: Overconsumption of solid species occurs in the particle.',   &
            /'Its source term will be set to its amount in the particle.',&
            /'To resolve this issue, please consider:',                  &
            /' a) Using stiff solver in the Chemistry pane ',      &
            /' b) Reducing the reaction rates ',&
            /' c) Reducing the time step ',&
            /'Global particle ID:',I8, &
            /'Species index:',I8)


! Always use First-order method: Euler
! when integrating species mass fractions
      ! IF(INTG_EULER) THEN

! sum_des_rs is the rate (unit amount of sustance per time) divided
! by pmass
         WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE)            &
            PMASS(:MAX_PIP) = PMASS(:MAX_PIP) + lDT*                   &
               SUM_DES_Rs(:MAX_PIP)*PMASS(:MAX_PIP)

         FORALL(NN=1:DIMENSION_N_S)
            WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE .AND.    &
               SPECIES_EQ(PIJK(:MAX_PIP,5)))         &
               DES_X_s(:MAX_PIP,NN) = max(DES_X_s(:MAX_PIP,NN) + lDT*  &
                  (DES_R_s(:MAX_PIP,NN) - DES_X_s(:MAX_PIP,NN)*        &
                  SUM_DES_Rs(:MAX_PIP)), ZERO)
         END FORALL

         WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE .AND.    &
               SPECIES_EQ(PIJK(:MAX_PIP,5)))             &
            sum_X_s(:MAX_PIP) = sum(DES_X_s(:MAX_PIP,:), DIM=2)

! renormalize the species mass fractions
         FORALL(NN=1:DIMENSION_N_S)
            WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE .AND.    &
               SPECIES_EQ(PIJK(:MAX_PIP,5)))  &
               DES_X_s(:MAX_PIP,NN) = DES_X_s(:MAX_PIP,NN)  / sum_X_s(:MAX_PIP)
         END FORALL

      ! ELSE
         ! IF(FIRST_PASS) THEN
         ! ENDIF
      ! ENDIF

      IF(SuperDEM) THEN
         DO NN=1,MAX_PIP
            IF(IS_NORMAL(NN)) THEN

               axi(1)=super_r(NN,1)
               axi(2)=super_r(NN,2)
               axi(3)=super_r(NN,3)
               m=super_mn(NN,1)
               n=super_mn(NN,2)

               IF(SOLVE_ROs(PIJK(NN,5))) THEN     ! variable density, constant size
                  RO_Sol(NN)= PMASS(NN)/PVOL(NN)
                  call SQ_INERTIA(axi,m,n,IXX,IYY,IZZ)
                  OMOI3(NN,1)=1.0/(IXX*RO_SOL(NN))
                  OMOI3(NN,2)=1.0/(IYY*RO_SOL(NN))
                  OMOI3(NN,3)=1.0/(IZZ*RO_SOL(NN))
               ELSE                               ! constant density, variable size
! Assume the particle shape remains the same (same aspect ratio, same m,n values)
! Then we can scale a,b,c by multiplying them by (new mass/ old mass)^(1/3)
! Here the old mass is still PVOL(NN)*RO_SOL(NN) because the chemical reaction
! has only changed PMASS at this point.
                  scale_factor = (PMASS(NN)/(PVOL(NN)*RO_SOL(NN)))**(o3)
                  axi(:)       = scale_factor * axi(:)

                  CALL SQ_VOLUME (axi,m,n,svolume)
                  PVOL(NN)=svolume
                  call SQ_INERTIA(axi,m,n,IXX,IYY,IZZ)
                  OMOI3(NN,1)=1.0/(IXX*RO_SOL(NN))
                  OMOI3(NN,2)=1.0/(IYY*RO_SOL(NN))
                  OMOI3(NN,3)=1.0/(IZZ*RO_SOL(NN))

                  super_r(NN,1)  = axi(1)
                  super_r(NN,2)  = axi(2)
                  super_r(NN,3)  = axi(3)
                  DES_RADIUS(NN) = scale_factor * DES_RADIUS(NN)
               ENDIF
            ENDIF
         ENDDO
      ELSEIF (GSP_EXPLICIT) THEN !Glued sphere particles
         glueVolume(:) = ZERO
         glueInertia(:,:) = ZERO
         allocate(new_total_mass(size(glueMass)))
         new_total_mass(:) = ZERO
         DO NN=1,MAX_PIP
            IF(.NOT. IS_NORMAL(NN)) CYCLE
            gid = gid_list(NN)
            IF(SOLVE_ROs(PIJK(NN,5))) THEN     ! variable density, constant size
               RO_Sol(NN)= PMASS(NN)/PVOL(NN)
               MOI = (2.0d0/5.0d0)*PMASS(NN)*DES_RADIUS(NN)**2
               dist2 = sc2gpc_vec(NN,2)*sc2gpc_vec(NN,2) + sc2gpc_vec(NN,3)*sc2gpc_vec(NN,3)
               glueInertia(gid,1) = glueInertia(gid,1) + PMASS(NN) * dist2 + MOI

               dist2 = sc2gpc_vec(NN,1)*sc2gpc_vec(NN,1) + sc2gpc_vec(NN,3)*sc2gpc_vec(NN,3)
               glueInertia(gid,2) = glueInertia(gid,2) + PMASS(NN) * dist2 + MOI

               dist2 = sc2gpc_vec(NN,2)*sc2gpc_vec(NN,2) + sc2gpc_vec(NN,1)*sc2gpc_vec(NN,1)
               glueInertia(gid,3) = glueInertia(gid,3) + PMASS(NN) * dist2 + MOI

               glueVolume(gid) = glueVolume(gid) + PVOL(NN)
            ELSE                               ! constant density, variable size
               is_constant_density = .True.
               new_total_mass(gid) = new_total_mass(gid) + pmass(NN)
            ENDIF
            OMOI(NN) = 2.5D0/(PMASS(NN)*DES_RADIUS(NN)**2) !UPDATE ONE OVER MOI
         ENDDO ! end of particle loop

         IF(numPEs > 1 .and. is_constant_density) THEN
            CALL GLOBAL_ALL_SUM(new_total_mass)
         ENDIF

         IF(is_constant_density) THEN
            DO NN = 1, MAX_PIP
               IF(.NOT. IS_NORMAL(NN)) CYCLE
               gid = gid_list(NN)
               IF(.NOT. SOLVE_ROS(PIJK(NN,5))) THEN
                  ! if not expose surface area, calculate size change based on total mass change
                  scale_factor = (new_total_mass(gid) / glueMass(gid)) ** (1.0/3.0)

                  des_radius(NN) = des_radius(NN) * scale_factor
                  pvol(NN) = (4.0/3.0)*pi*des_radius(NN) ** 3.0
                  gp_neighsa(NN,:) = gp_neighsa(NN,:) * scale_factor ** 2.0
                  gp_sa(NN) = gp_sa(NN) * scale_factor ** 2.0

                  sc2gpc_vec(NN,:) = sc2gpc_vec(NN,:) * scale_factor
                  moi = (2.0d0/5.0d0)*pmass(NN)*des_radius(NN)**2
                  dist2 = sc2gpc_vec(NN,2)*sc2gpc_vec(NN,2) + sc2gpc_vec(NN,3)*sc2gpc_vec(NN,3)
                  glueinertia(gid,1) = glueinertia(gid,1) + pmass(NN) * dist2 + moi

                  dist2 = sc2gpc_vec(NN,1)*sc2gpc_vec(NN,1) + sc2gpc_vec(NN,3)*sc2gpc_vec(NN,3)
                  glueinertia(gid,2) = glueinertia(gid,2) + pmass(NN) * dist2 + moi

                  dist2 = sc2gpc_vec(NN,2)*sc2gpc_vec(NN,2) + sc2gpc_vec(NN,1)*sc2gpc_vec(NN,1)
                  glueinertia(gid,3) = glueinertia(gid,3) + pmass(NN) * dist2 + moi

                  glueVolume(gid) = glueVolume(gid) + pvol(NN)
               ENDIF ! end of constant density reaction
            ENDDO
         ENDIF

         IF(numPEs > 1) THEN
            CALL GLOBAL_ALL_SUM(glueVolume)
            CALL GLOBAL_ALL_SUM(glueInertia)
         ENDIF

         glueDiameter = (6.0*glueVolume/PI) ** (1.0/3.0)
      ELSE
         DO NN=1,MAX_PIP
            IF(IS_NORMAL(NN)) THEN
               IF(SOLVE_ROs(PIJK(NN,5))) THEN     ! variable density, constant size
                  RO_Sol(NN)= PMASS(NN)/PVOL(NN)
               ELSE                               ! constant density, variable size
                  IF(GSP_IMPLICIT) THEN
                     WRITE(ERR_MSG, 1099)
                     CALL LOG_ERROR()
                  ENDIF
 1099 FORMAT('FATAL - variable size chemical reaction is not supported in gsp implicit model.')
                  DES_RADIUS(NN) = (PMASS(NN)/(Pix4o3*RO_SOL(NN)))**o3
                  PVOL(NN) = PMASS(NN)/RO_SOL(NN)
!llu, 11/7/18, update real particle size
                  if(CGDEM) DES_CGP_RPR(NN) = DES_RADIUS(NN) / DES_CGP_STW(NN)**o3
               ENDIF
               OMOI(NN) = 2.5D0/(PMASS(NN)*DES_RADIUS(NN)**2) !UPDATE ONE OVER MOI
            ENDIF
         ENDDO
      ENDIF

! Clear the necessary variables.
      DES_R_s = ZERO

! Flag that the first pass is over
      FIRST_PASS = .FALSE.

      RETURN

   END SUBROUTINE DES_REACTION_MODEL

END MODULE DES_REACTION_MODEL_MOD
