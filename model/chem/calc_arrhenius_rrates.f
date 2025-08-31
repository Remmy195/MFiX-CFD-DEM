#include "error.inc"
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_ARRHENIUS_RRATES(IJK, RATES)                       C
!  Purpose: Calculate reaction rates for various reactions in cell ijk C
!           using information from the data file. Unit: kmol/m3/s.     C
!                                                                      C
!  Author: Hang Zhou                                   Date: July 2024 C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_ARRHENIUS_RRATES(IJK, RATES)

! Global domain parameters.
!`````````````````````````````````````````````````````````````````````//

      USE constant, ONLY: GAS_CONST_cal
      USE fldvar, ONLY: EP_g, P_g, Ro_g, T_g, X_g, RO_s, T_s, X_s
      ! Number of TFM reactions.
      USE parse
      USE param1, ONLY: one, undefined_c
      USE physprop, ONLY: Mw_g, Mw_s
      USE read_thermochemical, ONLY:calc_H0oRT, calc_S0oR
      USE run, ONLY: tfm_solids
      USE rxns, ONLY: NO_OF_RXNS, Reaction, SPECIES_ALIAS_g
      USE toleranc, ONLY: chem_min_species_solid, chem_min_species_fluid, tmax, tmin


      IMPLICIT NONE

! Local variables:
!`````````````````````````````````````````````````````````````````````//

      INTEGER, INTENT(IN) :: IJK

      DOUBLE PRECISION, DIMENSION(NO_OF_RXNS), INTENT(OUT) :: RATES
! index of reactions, species, phase
      INTEGER :: rr, lN, NN, M, PP, LL
      DOUBLE PRECISION :: gas_constant, exp_factor, P_atm
! summation of species stoichiometric coefficients
      DOUBLE PRECISION :: sum_sto
! forward or reverse rate constants
      DOUBLE PRECISION :: k_rate
! equilibrium constant in pressure and concentration units
      DOUBLE PRECISION :: k_equ_p, k_equ_c
! concentration of gas species
      DOUBLE PRECISION, DIMENSION(DIMENSION_N_G) :: Cgi
! concentration of solid species
      DOUBLE PRECISION, DIMENSION(DIMENSION_M, DIMENSION_N_S) :: Csi
! entropy of all gas species over gas constant
      DOUBLE PRECISION, DIMENSION(DIMENSION_N_G) :: SioR_g
! enthalpy of all gas species over (gas constant * Tg)
      DOUBLE PRECISION, DIMENSION(DIMENSION_N_G) :: HioRT_g
! change of entropy over gas consntat and ehthalpy over (gas_constant*Tg) for each reaction
      DOUBLE PRECISION :: dS0, dH0
! concentration of third body
      DOUBLE PRECISION :: M_Body
! specific species used as third body
      CHARACTER(len=32):: M_species
! parameters for pressure-dependent reactions
      DOUBLE PRECISION :: k_extrapre, Pr, f_press
! for TROE model
      DOUBLE PRECISION :: f_cent, f_press_log, c_troe, n_troe, d_troe, alpha_troe, T3_troe, T1_troe, T2_troe
! for SRI model
      DOUBLE PRECISION :: a_sri, b_sri, c_sri, d_sri, e_sri, X_sri
! for PLOG model
      DOUBLE PRECISION, ALLOCATABLE :: k_plog(:), press_plog(:)
      INTEGER, ALLOCATABLE :: index_press_plog(:)
      INTEGER :: num_press_plog !number of different pressures for PLOG model
      INTEGER :: end_index_plog

      ! Local representation of the number of solids phases.
      INTEGER :: lMMAX

! Flag indicating mass fractions of all reactants are greater than chem_min_species_fluid or chem_min_species_solid
      LOGICAL :: reactant_sufficient
! Flag indicating if there is gas phase in the reaction
      LOGICAL :: has_gas
! Flag indicating if there is solid phase in the reaction
      LOGICAL :: has_solid

      CHARACTER(LEN=32) species_tmp, species_tmp2
! temperature used to calculate the reaction rates
      DOUBLE PRECISION :: Tgs

      gas_constant = 4.183925d0 * GAS_CONST_cal !(J/mol.K)
      P_atm = P_g(ijk)/101325.0d0  ! atm

      RATES(:) = 0.0d0

      Do LL = 1, DIMENSION_N_G
         Cgi(LL) = RO_g(IJK)*X_g(IJK,LL)/Mw_g(LL)  ! kmol/m3
         HioRT_g(LL) = calc_H0oRT(T_g(ijk), 0, LL)
         SioR_g(LL) = calc_S0oR(T_g(ijk), 0, LL)
      ENDDO
      IF(TFM_SOLIDS) THEN
         DO M=1, DIMENSION_M
            DO LL = 1, DIMENSION_N_S
               Csi(M, LL) = RO_s(IJK, M)*X_s(IJK,M,LL)/Mw_s(M,LL)  ! kmol/m3
	    ENDDO
         ENDDO
      ENDIF

! Calculate reaction rates based on arrhenius coefficients.
! Loop over reactions.
      RXN_LP: DO rr = 1, NO_OF_RXNS
         ! check if there are enough reactants for the reactions
         reactant_sufficient = .TRUE.
         has_gas = .FALSE.
         has_solid = .FALSE.
         DO lN = 1, Reaction(rr)%nSpecies
            ! check for reactants
            IF(Reaction(rr)%Species(lN)%Coeff .LT. 0) THEN
               ! Global phase index.
               M = Reaction(rr)%Species(lN)%pMap
               ! Global species index.
               NN = Reaction(rr)%Species(lN)%sMap
               IF(M .EQ. 0) THEN
                  has_gas = .TRUE.
	          IF(X_g(IJK, NN) .LT. chem_min_species_fluid) THEN
                     reactant_sufficient = .FALSE.
                     EXIT
		  ENDIF
               ELSE
                  has_solid = .TRUE.
                  Tgs = max(min(TMAX, T_s(IJK,M)), TMIN)
	          IF(X_s(IJK, M, NN) .LT. chem_min_species_solid) THEN
                     reactant_sufficient = .FALSE.
                     EXIT
		  ENDIF
               ENDIF
            ENDIF
         ENDDO

         IF(reactant_sufficient) THEN
             IF(has_solid .and. has_gas) THEN
                 Tgs = (max(min(TMAX, T_g(IJK)), TMIN) + Tgs) * HALF
             ELSEIF(has_gas) THEN
                 Tgs = max(min(TMAX, T_g(IJK)), TMIN)
             ENDIF
             exp_factor = 1/Tgs/(gas_constant*1000.d0) !1/RT, kmol/J (1000 is used to convert mol to kmol)
            ! Skip empty reactions
            IF(Reaction(rr)%nSpecies == 0) CYCLE RXN_LP
            ! getting the forward rate constant
            IF(arrhenius_coeff(rr,2)==0.0d0) THEN
	        k_rate = arrhenius_coeff(rr,1)* exp(-arrhenius_coeff(rr,3)*exp_factor)
            ELSEIF(arrhenius_coeff(rr,2)==1.0d0) THEN
	        k_rate = arrhenius_coeff(rr,1) * T_g(IJK) &
	                   * exp(-arrhenius_coeff(rr,3)*exp_factor)
            ELSEIF(arrhenius_coeff(rr,2)==-1.0d0) THEN
	        k_rate = arrhenius_coeff(rr,1) / T_g(IJK) &
	                   * exp(-arrhenius_coeff(rr,3)*exp_factor)
	    ELSE
               k_rate = arrhenius_coeff(rr,1) * T_g(IJK)**arrhenius_coeff(rr,2) &
	                * exp(-arrhenius_coeff(rr,3)*exp_factor)
            ENDIF
	        ! check if there are third body in the reaction
            IF(has_solid .AND. (third_body_model(rr) .NE. UNDEFINED_C)) THEN
               WRITE(ERR_MSG,"(/1X,70('*')/' From: rrates_arrhenius:',/   &
               ' Message: Third body information is given for reaction with solid phase, ',/  &
                 A,', which is not supported!')") trim(Reaction(rr)%ChemEq)
               CALL log_error()
            ENDIF
            IF(third_body_model(rr) .EQ. UNDEFINED_C) THEN
                M_body = 1
            ELSEIF(third_body_model(rr) .EQ. 'M_ALL') THEN
                M_body = 0
                Do LL = 1, DIMENSION_N_G
                    M_body = M_body + Cgi(LL)*third_body_coeff(rr,LL)
                ENDDO
            ELSE
                M_species = trim(third_body_model(rr)(3:))
                DO LL = 1, DIMENSION_N_G
                    species_tmp = trim(SPECIES_ALIAS_g(LL))
                    species_tmp2 =trim(M_species)
                    ! Remove case sensitivity.
                    CALL MAKE_UPPER_CASE (species_tmp,32)
                    CALL MAKE_UPPER_CASE (species_tmp2,32)
                    IF(trim(species_tmp) .EQ. trim(species_tmp2)) THEN
                        M_body = Cgi(LL)*third_body_coeff(rr,LL)
                        EXIT
                    ENDIF
                 ENDDO
	     ENDIF
            ! check for pressure-related parameters
            IF(has_solid .AND. (press_rxn_model(rr) .NE. UNDEFINED_C)) THEN
               WRITE(ERR_MSG,"(/1X,70('*')/' From: rrates_arrhenius:',/   &
               ' Message: Pressure-related parameters are given for reaction with solid phase, ',/  &
                 A,', which is not supported!')") trim(Reaction(rr)%ChemEq)
               CALL log_error()
            ENDIF
            IF(press_rxn_model(rr) .NE. UNDEFINED_C) THEN
               IF(INDEX(press_rxn_model(rr), "PLOG") .NE. 0) THEN
                  ALLOCATE(press_plog(SIZE(press_rxn_coeff(rr,:))))
                  ALLOCATE(index_press_plog(SIZE(press_rxn_coeff(rr,:))))
                  press_plog = UNDEFINED
                  index_press_plog = UNDEFINED_I

                  num_press_plog = 1
                  index_press_plog(num_press_plog) = 1
                  press_plog(num_press_plog) = press_rxn_coeff(rr,1)
                  DO LL = 2, SIZE(press_rxn_coeff(rr,:))
                     IF((press_rxn_coeff(rr,LL) .NE. UNDEFINED) .AND. &
                        (press_rxn_coeff(rr,LL) .NE. press_plog(num_press_plog))) THEN
                        num_press_plog = num_press_plog + 1
                        index_press_plog(num_press_plog) = LL
                        press_plog(num_press_plog) = press_rxn_coeff(rr,LL)
                     ENDIF
                  ENDDO
                  ALLOCATE(k_plog(num_press_plog))
                  k_plog = 0.d0
                  DO LL = 1, num_press_plog
                     IF(LL .NE. num_press_plog) THEN
                        end_index_plog = index_press_plog(LL+1)-1
                     ELSE
                        end_index_plog = SIZE(press_rxn_coeff(rr,:))
                        DO WHILE (press_rxn_coeff(rr,end_index_plog) .EQ. UNDEFINED)
                           end_index_plog = end_index_plog -1
                        ENDDO
                     ENDIF
                     DO PP = index_press_plog(LL), end_index_plog
                        IF(arrhenius_coeff_press(rr,PP,2)==0.0d0) THEN
                           k_plog(LL) = k_plog(LL) + arrhenius_coeff_press(rr,PP,1) &
                                      * exp(-arrhenius_coeff_press(rr,PP,3)*exp_factor)
                        ELSEIF(arrhenius_coeff_press(rr,PP,2)==1.0d0) THEN
                           k_plog(LL) = k_plog(LL) + arrhenius_coeff_press(rr,PP,1) * T_g(IJK) &
                                      * exp(-arrhenius_coeff_press(rr,PP,3)*exp_factor)
                        ELSEIF(arrhenius_coeff_press(rr,PP,2)==-1.0d0) THEN
                           k_plog(LL) = k_plog(LL) + arrhenius_coeff_press(rr,PP,1) / T_g(IJK) &
                                      * exp(-arrhenius_coeff_press(rr,PP,3)*exp_factor)
                        ELSE
                           k_plog(LL) = k_plog(LL) + arrhenius_coeff_press(rr,PP,1) * T_g(IJK)**arrhenius_coeff_press(rr,PP,2) &
	                              * exp(-arrhenius_coeff_press(rr,PP,3)*exp_factor)
                        ENDIF
                     ENDDO
                  ENDDO
                  IF(P_atm .LE. press_plog(1)) THEN
                      k_rate = k_plog(1)
                  ELSEIF(P_atm .GT. press_plog(num_press_plog)) THEN
                      k_rate = k_plog(num_press_plog)
                  ELSE
                      DO LL = 2, num_press_plog
                         IF(P_atm .LE. press_plog(LL)) THEN
                            k_rate = exp(LOG(k_plog(LL-1))+(LOG(P_atm)-LOG(press_plog(LL-1)))/(LOG(press_plog(LL))-LOG(press_plog(LL-1))) &
                                                            *(LOG(k_plog(LL))-LOG(k_plog(LL-1))))
                            EXIT
                         ENDIF
                      ENDDO
                  ENDIF
                  DEALLOCATE(press_plog)
                  DEALLOCATE(index_press_plog)
                  DEALLOCATE(k_plog)
               ELSE
                  ! get the reaction rate constant at low (for fall-off reactions)
                  ! or high pressure (for chemically activated bimolecular reactions)
                  IF(arrhenius_coeff_press(rr,1,2)==0.0d0) THEN
                      k_extrapre = arrhenius_coeff_press(rr,1,1) &
                                * exp(-arrhenius_coeff_press(rr,1,3)*exp_factor)
                  ELSEIF(arrhenius_coeff_press(rr,1,2)==1.0d0) THEN
                      k_extrapre = arrhenius_coeff_press(rr,1,1) * T_g(IJK) &
	                          * exp(-arrhenius_coeff_press(rr,1,3)*exp_factor)
                  ELSEIF(arrhenius_coeff_press(rr,1,2)==-1.0d0) THEN
                      k_extrapre = arrhenius_coeff_press(rr,1,1) / T_g(IJK) &
	                          * exp(-arrhenius_coeff_press(rr,1,3)*exp_factor)
                  ELSE
                      k_extrapre = arrhenius_coeff_press(rr,1,1) * T_g(IJK)**arrhenius_coeff_press(rr,1,2) &
	                          * exp(-arrhenius_coeff_press(rr,1,3)*exp_factor)
                  ENDIF

                  IF(INDEX(press_rxn_model(rr), 'FALLOFF') .NE. 0) THEN
                     ! get reduced pressure Pr
                     Pr = k_extrapre/k_rate * M_body
                  ELSEIF(INDEX(press_rxn_model(rr), 'BIMO') .NE. 0) THEN
                     ! get reduced pressure Pr
                     Pr = k_rate/k_extrapre * M_body
                  ENDIF
                  ! Lindemann model
                  IF(INDEX(press_rxn_model(rr), 'LINDEMANN') .NE. 0) THEN
                     f_press = 1
                  ! TROE model
                  ELSEIF(INDEX(press_rxn_model(rr), 'TROE') .NE. 0) THEN
                     alpha_troe = press_rxn_coeff(rr,1)
                     T3_troe = press_rxn_coeff(rr,2)
                     T1_troe = press_rxn_coeff(rr,3)
                     T2_troe = press_rxn_coeff(rr,4)
                     f_cent = (1-alpha_troe)*exp(-T_g(IJK)/T3_troe)+alpha_troe*exp(-T_g(IJK)/T1_troe)
                     IF(T2_troe .NE. UNDEFINED) f_cent = f_cent + exp(-T2_troe/T_g(IJK))
                     c_troe = -0.4d0-0.67d0*LOG10(f_cent)
                     n_troe = 0.75d0-1.27d0*LOG10(f_cent)
                     d_troe = 0.14d0
                     f_press_log = LOG10(f_cent)/(1+((LOG10(Pr)+c_troe)/(n_troe-d_troe*(LOG10(Pr)+c_troe)))**2.0d0)
                     f_press = 10**(f_press_log)
                  ! SRI model
                  ELSEIF(INDEX(press_rxn_model(rr), 'SRI') .NE. 0) THEN
                     a_sri = press_rxn_coeff(rr,1)
                     b_sri = press_rxn_coeff(rr,2)
                     c_sri = press_rxn_coeff(rr,3)
                     d_sri = press_rxn_coeff(rr,4)
                     e_sri = press_rxn_coeff(rr,5)
                     X_sri = 1/(1+LOG10(Pr)*LOG10(Pr))
                     f_press = a_sri*exp(-b_sri/T_g(IJK))+exp(-T_g(IJK)/c_sri)**X_sri
                     IF(d_sri .NE. UNDEFINED .AND. e_sri .NE. UNDEFINED) THEN
                        f_press = f_press*d_sri*T_g(IJK)**e_sri
                     ENDIF
                  ENDIF

                  IF(INDEX(press_rxn_model(rr), 'FALLOFF') .NE. 0) THEN
                     k_rate = k_rate*Pr/(1+Pr)*f_press
                  ELSEIF(INDEX(press_rxn_model(rr), 'BIMO') .NE. 0) THEN
                     k_rate = k_rate/(1+Pr)*f_press
                  ENDIF
               ENDIF
            ENDIF
            ! check for Landau-Teller reactions
            IF(LT_coeff(rr,1) .NE. UNDEFINED) THEN
               IF(ANY( LT_coeff(rr,:)== UNDEFINED)) THEN
                  WRITE(ERR_MSG,"(/1X,70('*')/' From: rrates_arrhenius:',/   &
                          ' Message: 2 parameters must be givn for Landau-Teller type of reactions.',/  &
                           A)") trim(Reaction(rr)%ChemEq)
                  CALL log_error()
               ENDIF
               k_rate = arrhenius_coeff(rr,1) * T_g(IJK)**arrhenius_coeff(rr,2) * &
                        exp(-arrhenius_coeff(rr,3)*exp_factor + &
			LT_coeff(rr,1)/T_g(IJK)**(1.0d0/3.0d0) + LT_coeff(rr,2)/T_g(IJK)**(2.0d0/3.0d0))
            ENDIF
            ! check for other rate constant fitting options
            IF(rate_fit_model(rr) .NE. UNDEFINED_C) THEN
               IF(rate_fit_model(rr) .EQ. "JAN") THEN
                  IF( ANY( rate_fit_coeff(rr,:)== UNDEFINED)) THEN
                     WRITE(ERR_MSG,"(/1X,70('*')/' From: rrates_arrhenius:',/   &
                          ' Message: 9 parameters must be givn for JAN type of reactions.',/  &
                           A)") trim(Reaction(rr)%ChemEq)
                     CALL log_error()
                  ENDIF
                  k_rate = arrhenius_coeff(rr,1) * T_g(IJK)**arrhenius_coeff(rr,2) &
		           * exp(arrhenius_coeff(rr,3)/T_g(IJK)+ &
                                rate_fit_coeff(rr,1) + &
                                rate_fit_coeff(rr,2)*LOG(T_g(IJK)) + &
                                rate_fit_coeff(rr,3)*(LOG(T_g(IJK)))**2.0d0 + &
                                rate_fit_coeff(rr,4)*(LOG(T_g(IJK)))**3.0d0 + &
                                rate_fit_coeff(rr,5)*(LOG(T_g(IJK)))**4.0d0 + &
                                rate_fit_coeff(rr,6)*(LOG(T_g(IJK)))**5.0d0 + &
                                rate_fit_coeff(rr,7)*(LOG(T_g(IJK)))**6.0d0 + &
                                rate_fit_coeff(rr,8)*(LOG(T_g(IJK)))**7.0d0 + &
                                rate_fit_coeff(rr,9)*(LOG(T_g(IJK)))**8.0d0)

               ELSEIF(rate_fit_model(rr) .EQ. "FIT1") THEN
                  IF( ANY( rate_fit_coeff(rr,1:4)== UNDEFINED)) THEN
                     WRITE(ERR_MSG,"(/1X,70('*')/' From: rrates_arrhenius:',/   &
                          ' Message: 4 parameters must be givn for FIT1 type of reactions.',/  &
                           A)") trim(Reaction(rr)%ChemEq)
                     CALL log_error()
                  ENDIF
                  k_rate = arrhenius_coeff(rr,1) * T_g(IJK)**arrhenius_coeff(rr,2) &
                           *exp(rate_fit_coeff(rr,1)/T_g(IJK) + &
                                rate_fit_coeff(rr,2)/(T_g(IJK)**2.0d0) + &
                                rate_fit_coeff(rr,3)/(T_g(IJK)**3.0d0) + &
                                rate_fit_coeff(rr,4)/(T_g(IJK)**4.0d0))
               ENDIF
               IF(arrhenius_coeff(rr,2)==0.0d0) THEN
	               k_rate = k_rate
               ELSEIF(arrhenius_coeff(rr,2)==1.0d0) THEN
	               k_rate = k_rate * T_g(IJK)
               ELSEIF(arrhenius_coeff(rr,2)==-1.0d0) THEN
	               k_rate = k_rate / T_g(IJK)
	            ELSE
                  k_rate = k_rate * T_g(IJK)**arrhenius_coeff(rr,2)
               ENDIF
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
            IF(INDEX(reverse_calc(rr), "FROMFORWARDRATECONSTANT") .NE. 0) THEN
               IF(has_solid) THEN
                  WRITE(ERR_MSG,"(/1X,70('*')/' From: rrates_arrhenius:',/   &
                  ' Message: FROMFORWARDRATECONSTANT cannot be used for the reverse reaction with solid phase. ',/  &
                    A,'. The Arrhenius coefficients must be given by arrhenius_coeff!')") trim(Reaction(rr)%ChemEq)
                  CALL log_error()
               ENDIF
	       dS0 = 0.0d0
	       dH0 = 0.0d0
	       sum_sto  = 0.0d0
               ! sum_sto, dS0 and dH0 should be calculated based on the forward reaction
	       DO lN = 1, Reaction(rr)%nSpecies
	          ! Global species index.
                  NN = Reaction(rr)%Species(lN)%sMap
	          sum_sto = sum_sto - Reaction(rr)%Species(lN)%Coeff
	          dS0 = dS0 - Reaction(rr)%Species(lN)%Coeff*SioR_g(NN)
                  dH0 = dH0 - Reaction(rr)%Species(lN)%Coeff*HioRT_g(NN)
	       ENDDO
               k_equ_p = exp(dS0-dH0)
	       IF(sum_sto .EQ. 0.0d0) THEN
	          k_equ_c = k_equ_p
	       ELSEIF(sum_sto .EQ. 1.0d0) THEN
	          k_equ_c = k_equ_p*(101.325/gas_constant/T_g(ijk))
	       ELSEIF(sum_sto .EQ. -1.0d0) THEN
	          k_equ_c = k_equ_p/(101.325/gas_constant/T_g(ijk))
	       ELSE
	          k_equ_c = k_equ_p*(101.325/gas_constant/T_g(ijk))**sum_sto
	       ENDIF
	       k_rate = k_rate/k_equ_c
	    ENDIF

            IF(press_rxn_model(rr) .NE. UNDEFINED_C .AND. INDEX(press_rxn_model(rr), "PLOG") .EQ. 0) THEN
               RATES(rr) = EP_g(IJK)* k_rate
            ELSE
               RATES(rr) = M_body*EP_g(IJK)* k_rate
            ENDIF
            DO lN = 1, Reaction(rr)%nSpecies
               ! check for reactants
               IF(Reaction(rr)%Species(lN)%Coeff .LT. 0) THEN
                  ! Global phase index.
                  M = Reaction(rr)%Species(lN)%pMap
                  ! Global species index.
                  NN = Reaction(rr)%Species(lN)%sMap
                  IF(M .EQ. 0) THEN
                     RATES(rr) = RATES(rr) * Cgi(NN)**Reaction(rr)%Species(lN)%RxnOrder
                  ELSE
                     RATES(rr) = RATES(rr) * Csi(M, NN)**Reaction(rr)%Species(lN)%RxnOrder
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO RXN_LP ! Loop over reactions.

      RETURN

      END SUBROUTINE CALC_ARRHENIUS_RRATES
