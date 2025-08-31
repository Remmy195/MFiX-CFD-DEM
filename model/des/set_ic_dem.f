#include "error.inc"

MODULE SET_IC_DEM_MOD

   use error_manager

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE Name: DES_SET_IC                                         !
!                                                                      !
!  Purpose: Assign initial conditions to particles based on their      !
!  location within the domain.                                         !
!                                                                      !
!  Author: J.Musser                                   Date: 15-Feb-11  !
!  Revised: Renjie Ke                                 Data: 09-Aug-24  !
!                                                                      !
!  Comments:                                                           !
!  Revised Comments:                                                   !
!     Providing the flexibility for user to using IC values to         !
!     values given in particle_input files                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE SET_IC_DEM

      use run, only: ENERGY_EQ, SPECIES_EQ

      use ic

      use des_thermo

      use derived_types, only: PIC
      use discretelement, only : gener_part_config, P_INPUT_DAT_VERSION
      use discretelement, only: MAX_PIP
      use discretelement, only: PINC
      use discretelement, only: PIJK
      use discretelement, only: particle_state, normal_particle, particle_state_gsp
      use discretelement, only: gid_list, gluePos,GluedSphereDEM, gp_component_IC, gsp_explicit
      use discretelement, only: des_vel_new

      USE des_rxns, only: DES_X_s

      use mfix_pic, only: MPPIC

      use param1, only: undefined, zero

      use physprop, only: C_PS0
      use physprop, only: NMAX

      USE compar
      use indices
      use geometry

      use functions
      use toleranc

      IMPLICIT NONE

! Dummy indices
      INTEGER :: ICV, lICV
      INTEGER :: I, J, K, IJK
      INTEGER :: M, NN
      INTEGER :: NP
      INTEGER :: NINDX
      INTEGER :: gid
      DOUBLE PRECISION :: gPos(3)
      LOGICAL :: pIver_1, pIver_4, pVer

! Particle temperature and species only need to be set here
! from IC region settings if using particle_input.dat version 1.0
! If using particle_input.dat version 2.0 and above, temperature
! and species are set from data in particle_input.dat.
! If using IC region settings (gener_part_config=.true.), temperature
! and species are set in ADD_PARTICLE subroutine (generate_particles_mod.f).
! Temperature and species are also set here when using MPPIC

! It must not auto seeding!
      pIver_1 = trim(P_INPUT_DAT_VERSION)=='1.0'
      pIver_4 = trim(P_INPUT_DAT_VERSION)=='4.0'
      pVer = (.NOT. GENER_PART_CONFIG) .AND. (pIver_4 .OR. pIver_1)

      IF (pVer .OR. MPPIC) THEN
         DO ICV = 1, DIMENSION_IC
            IF(.NOT.IC_DEFINED(ICV)) CYCLE
            IF(gsp_explicit .and. .not.(allocated(gp_component_IC))) CYCLE
! If all IC choose not to overwrite particle_input data, then CYCLE
! If at least one of IC choose to overwrite the data, executing the code below
            IF(.NOT. ANY(IC_OVERRIDE_PART_INPUT(ICV,:)) .AND. pIver_4) CYCLE

            DO K = IC_K_B(ICV), IC_K_T(ICV)
            DO J = IC_J_S(ICV), IC_J_N(ICV)
            DO I = IC_I_W(ICV), IC_I_E(ICV)

! Set the initial conditions for particles in cells that are
! not dead and that this rank owns.
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
               IF (DEAD_CELL_AT(I,J,K)) CYCLE

               IJK = FUNIJK(I,J,K)

! Loop through particles in cell IJK.
               DO NINDX = 1,PINC(IJK)
                  NP = PIC(IJK)%P(NINDX)
! Shift the phase index to the absolute phase index.
                  M = PIJK(NP,5)
! If this phase choose not to overwrite, then cycle to the next particle
                  IF(.NOT.IC_OVERRIDE_PART_INPUT(ICV,M) .AND. pIver_4) CYCLE

! Note GSP is a little bit tricky, components might be in a different IC
! Therefore, for GSP model, set the IC based on the gluePos
                  lICV = ICV
                  IF(gsp_explicit) THEN
                     ! In MPI, only overwrite values for normal particles, the update of ghost is through par_exchange
                     IF(PARTICLE_STATE(NP) .NE. NORMAL_PARTICLE) CYCLE
                     IF(PARTICLE_STATE_GSP(gid_list(NP)) .NE. NORMAL_PARTICLE) CYCLE
                     ! for GSP, this lICV might not equal to ICV
                     lICV = gp_component_IC(NP)
                     ! so need to check if this ICV determine to overwrite or not
                     IF(.NOT.IC_OVERRIDE_PART_INPUT(lICV,M)) CYCLE
                  ENDIF

! Set the initial particle velocities to overwrite what were read from particle_input data files
                  IF(pIver_4) THEN
                     des_vel_new(NP,1) = IC_U_s(lICV,M)
                     des_vel_new(NP,2) = IC_V_s(lICV,M)
                     des_vel_new(NP,3) = IC_W_s(lICV,M)
                  ENDIF

! Set the initial particle temperature.
                  IF(ENERGY_EQ) THEN
                     DES_T_s(NP) = IC_T_s(lICV,M)
                  ENDIF

! Set the initial species composition.
                  IF((ENERGY_EQ .AND. C_Ps0(M) == UNDEFINED) .OR. SPECIES_EQ(M)) THEN
                     DES_X_s(NP,:) = ZERO
                     DO NN = 1, NMAX(M)
                        DES_X_s(NP,NN) = IC_X_s(lICV,M,NN)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         ! to this end, gp_component_IC is not used any more
         IF(ALLOCATED(gp_component_IC)) DEALLOCATE(gp_component_IC)
      ENDIF

! Verify that all particles have a specified temperature and species
! mass fractions that sum to one. These checks are needed as the
! basic checks for IC_T_s and IC_X_s can be skipped if the IC region
! is specified with EPg = 1.
      DO NP = 1, MAX_PIP
! skipping non-existent particles
         IF(IS_NONEXISTENT(NP)) CYCLE
! skipping ghost particles
         IF(IS_ANY_GHOST(NP)) CYCLE

         M = PIJK(NP,5)

! Check that the temperature is specified.
         IF(ENERGY_EQ) THEN
            IF(DES_T_s(NP) == ZERO) THEN
               WRITE(ERR_MSG, 2000) trim(iVal(NP)), trim(iVal(lICV)), trim(iVal(M))
               CALL LOG_ERROR()
            ENDIF
         ENDIF

 2000 FORMAT('Error 2000: Particle ',A,' does not have a specified ',  &
         'initial',/'temperature. Verify that the IC region ',         &
         'containing this particle',/'has a solids temperature ',      &
         'defined: IC_T_s(',A,',',A,').')

! Check that the species mass fractions are specified.
         IF((ENERGY_EQ .AND. C_Ps0(M) == UNDEFINED) .OR.               &
            SPECIES_EQ(M)) THEN
            IF(.NOT.COMPARE(sum(DES_X_s(NP,1:NMAX(M))),ONE)) THEN
               WRITE(ERR_MSG, 2001) trim(iVal(NP)), trim(iVal(lICV)), trim(iVal(M))
               CALL LOG_ERROR()
            ENDIF
         ENDIF

 2001 FORMAT('Error 2001: The initial species mass fractions for ',     &
         'particle ',A,/'do not sum to 1. Verify that the IC ',    &
         'region containing this particle',/'has the solids species ', &
         'mass fractions defined: IC_X_s(',A,',',A,',:).')


      ENDDO

! Initialize specific heat
      IF(ENERGY_EQ) THEN
         FORALL(NP=1:MAX_PIP, PARTICLE_STATE(NP)==NORMAL_PARTICLE) &
            DES_C_PS(NP) = CALC_CP_DES(NP)
      ENDIF
      RETURN

   END SUBROUTINE SET_IC_DEM

END MODULE SET_IC_DEM_MOD
