#include "error.inc"
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_ARRHENIUS_RRATES_DES(NP, pM, IJK, DES_RATES)       C
!  Purpose: Calculate reaction rates for various reactions in cell ijk C
!           using information from the data file. Unit: kmol/s.        C
!                                                                      C
!  Author: Hang Zhou                                   Date: July 2024 C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_ARRHENIUS_RRATES_DES(NP, pM, IJK, DES_RATES)

! Global domain parameters.
!`````````````````````````````````````````````````````````````````````//

      USE constant, ONLY: GAS_CONST_cal
      USE des_rxns, only: DES_X_s, NO_OF_DES_RXNS, DES_Reaction
      use des_thermo, only: des_t_s
      use discretelement, only: ro_sol, pvol
      USE fldvar, ONLY: EP_g, P_g, Ro_g, T_g, X_g
      ! Number of TFM reactions.
      USE parse
      USE param1, ONLY: undefined_c
      USE physprop, ONLY: Mw_g, Mw_s
      USE toleranc, ONLY: chem_min_species_solid, chem_min_species_fluid, tmax, tmin


      IMPLICIT NONE

! Local variables:
!`````````````````````````````````````````````````````````````````````//

      INTEGER, INTENT(IN) :: NP  ! Global index of particle
      INTEGER, INTENT(IN) :: pM  ! Solid phase index of particle NP
      INTEGER, INTENT(IN) :: IJK ! Fluid cell index containing NP

! Calculated reaction rates. (reacted moles per sec)
      DOUBLE PRECISION, INTENT(OUT) :: DES_RATES(NO_OF_DES_RXNS)

! index of reactions, species, phase
      INTEGER :: rr, lN, NN, M, PP, LL
      DOUBLE PRECISION :: gas_constant, exp_factor, P_atm
! summation of species stoichiometric coefficients
      DOUBLE PRECISION :: sum_sto
! forward or reverse rate constants
      DOUBLE PRECISION :: k_rate
! concentration of gas species
      DOUBLE PRECISION, DIMENSION(DIMENSION_N_G) :: Cgi
! concentration of solid species
      DOUBLE PRECISION, DIMENSION(DIMENSION_M, DIMENSION_N_S) :: Csi

      ! Local representation of the number of solids phases.
      INTEGER :: lMMAX

! Flag indicating mass fractions of all reactants are greater than chem_min_species_fluid or chem_min_species_solid
      LOGICAL :: reactant_sufficient
! Flag indicating if there is gas phase in the reaction
      LOGICAL :: has_gas
! Flag indicating if there is solid phase in the reaction
      LOGICAL :: has_solid

! temperature used to calculate the reaction rates
      DOUBLE PRECISION :: Tgs

      gas_constant = 4.183925d0 * GAS_CONST_cal !(J/mol.K)
      P_atm = P_g(ijk)/101325.0d0  ! atm

      DES_RATES(:) = 0.0d0

      Do LL = 1, DIMENSION_N_G
         Cgi(LL) = RO_g(IJK)*X_g(IJK,LL)/Mw_g(LL)  ! kmol/m3
      ENDDO
      DO M=1, DIMENSION_M
         DO LL = 1, DIMENSION_N_S
            Csi(M, LL) = RO_SOL(NP)*DES_X_s(NP,LL)/Mw_s(pM,LL)  ! kmol/m3
         ENDDO
      ENDDO

! Calculate reaction rates based on arrhenius coefficients.
! Loop over reactions.
      RXN_LP: DO rr = 1, NO_OF_DES_RXNS
         ! check if there are enough reactants for the reactions
         reactant_sufficient = .TRUE.
         has_gas = .FALSE.
         has_solid = .FALSE.
         DO lN = 1, DES_Reaction(rr)%nSpecies
            ! check for reactants
            IF(DES_Reaction(rr)%Species(lN)%Coeff .LT. 0) THEN
               ! Global phase index.
               M = DES_Reaction(rr)%Species(lN)%pMap
               ! Global species index.
               NN = DES_Reaction(rr)%Species(lN)%sMap
               IF(M .EQ. 0) THEN
                   has_gas = .TRUE.
                   IF(X_g(IJK, NN) .LT. chem_min_species_fluid) THEN
                     reactant_sufficient = .FALSE.
                     EXIT
                   ENDIF
               ELSE
                   has_solid = .TRUE.
                   IF(DES_X_s(NP, NN) .LT. chem_min_species_solid) THEN
                     reactant_sufficient = .FALSE.
                     EXIT
                   ENDIF
               ENDIF
            ENDIF
         ENDDO

         IF(reactant_sufficient) THEN
             ! here, we assume only one solid phase will appear in a reaction
             IF(has_solid .and. has_gas) THEN
                 Tgs = (max(min(TMAX, T_g(IJK)), TMIN) + max(min(TMAX, DES_T_s(NP)), TMIN)) * HALF
             ELSEIF(has_gas) THEN
                 Tgs = max(min(TMAX, T_g(IJK)), TMIN)
             ELSEIF(has_solid) THEN
                 Tgs = max(min(TMAX, DES_T_s(NP)), TMIN)
             ENDIF
             exp_factor = 1/Tgs/(gas_constant*1000.d0) !1/RT, kmol/J
            ! Skip empty reactions
            IF(DES_Reaction(rr)%nSpecies == 0) CYCLE RXN_LP
            ! getting the forward rate constant
            IF(des_arrhenius_coeff(rr,2)==0.0d0) THEN
                k_rate = des_arrhenius_coeff(rr,1)* exp(-des_arrhenius_coeff(rr,3)*exp_factor)
            ELSEIF(des_arrhenius_coeff(rr,2)==1.0d0) THEN
                k_rate = des_arrhenius_coeff(rr,1) * T_g(IJK) &
                       * exp(-des_arrhenius_coeff(rr,3)*exp_factor)
            ELSEIF(des_arrhenius_coeff(rr,2)==-1.0d0) THEN
                k_rate = des_arrhenius_coeff(rr,1) / T_g(IJK) &
                       * exp(-des_arrhenius_coeff(rr,3)*exp_factor)
             ELSE
               k_rate = des_arrhenius_coeff(rr,1) * T_g(IJK)**des_arrhenius_coeff(rr,2) &
                    * exp(-des_arrhenius_coeff(rr,3)*exp_factor)
            ENDIF

       ! getting the rate constant for reverse reactions
        ! reverse_calc is UNDEFINED_C means the reaction is forward reaction
        ! reverse_calc is fromForwardRateConstant means that
        !  1. the reaction is the reverse reaction
        !  2. the rate constant needs to be calculated based on the forward rate constant
        ! reverse_calc is fromArrheniusCoeff or defined means
        !  1. the reaction is the reverse reaction
        !  2. arrhenius coefficients are given for the calculation of the rate constant,
        !     which has been calculated in the above calculation of k_rate
            IF(INDEX(des_reverse_calc(rr), "FROMFORWARDRATECONSTANT") .NE. 0) THEN
               WRITE(ERR_MSG,"(/1X,70('*')/' From: From: rrates_arrhenius_des:',/   &
               ' Message: FROMFORWARDRATECONSTANT cannot be used for the reverse DES reactions. ',/  &
                A,'. The Arrhenius coefficients must be given by arrhenius_coeff!')") trim(DES_Reaction(rr)%ChemEq)
               CALL log_error()
            ENDIF


            DES_RATES(rr) = PVOL(NP)* k_rate
            DO lN = 1, DES_Reaction(rr)%nSpecies
               ! check for reactants
               IF(DES_Reaction(rr)%Species(lN)%Coeff .LT. 0) THEN
                  ! Global phase index.
                  M = DES_Reaction(rr)%Species(lN)%pMap
                  ! Global species index.
                  NN = DES_Reaction(rr)%Species(lN)%sMap
                  IF(M .EQ. 0) THEN
                     DES_RATES(rr) = DES_RATES(rr) * Cgi(NN)**DES_Reaction(rr)%Species(lN)%RxnOrder
                  ELSE
                     DES_RATES(rr) = DES_RATES(rr) * Csi(M, NN)**DES_Reaction(rr)%Species(lN)%RxnOrder
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO RXN_LP ! Loop over reactions.

      RETURN

      END SUBROUTINE CALC_ARRHENIUS_RRATES_DES
