#include "error.inc"

MODULE parse

   USE compar
   USE error_manager
   USE funits
   USE make_upper_case_mod, only: make_upper_case
   USE param
   USE param1

      IMPLICIT NONE

! Strings indicating arithmetic operation and reaction blocks.
      CHARACTER(LEN=2), PARAMETER :: START_STR = '@('  ! start
      CHARACTER(LEN=1), PARAMETER :: END_STR = ')'     ! end

! Strings indicating reaction blocks.
      CHARACTER(LEN=4), PARAMETER :: RXN_BLK     = 'RXNS'      ! start block
      CHARACTER(LEN=8), PARAMETER :: DES_RXN_BLK = 'DES_RXNS'  ! start block
      CHARACTER(LEN=3), PARAMETER :: END_BLK     = 'END'       ! end block

      LOGICAL READING_RXN
      LOGICAL READING_RATE

      LOGICAL DES_RXN
      LOGICAL TFM_RXN

! Logical indicating that the start of a reaction construct has
! been identified.
      LOGICAL IN_CONSTRUCT

! Logical indicating that the chemical equation spans additional lines.
      LOGICAL MORE_ChemEq
! Logical indicating that the rxn orders spans additional lines.
      LOGICAL MORE_RxnOrder
! Logical indicating that the parameters for pressure-related reactions
! spans additional lines.
      LOGICAL MORE_PressParam
! Logical indicating that the parameters for optional rate-constant fit
! expression spans additional lines.
      LOGICAL MORE_fitParam
! Logical indicating that the coefficients of third body spans additional lines.
      LOGICAL MORE_thirdBody

! Reaction names
      CHARACTER(len=32),  DIMENSION(:),   ALLOCATABLE :: RXN_NAME
! chemical Equations
      CHARACTER(len=512), DIMENSION(:),   ALLOCATABLE :: RXN_CHEM_EQ
! User defined heat of reaction
      DOUBLE PRECISION,   DIMENSION(:),   ALLOCATABLE :: usrDH
! User defined heat of reaction partitions.
      DOUBLE PRECISION,   DIMENSION(:,:), ALLOCATABLE :: usrfDH

! flag to use Arrhenius model to calculate reaction rates for fluid phase
      LOGICAL :: ARRHENIUS_RRATES_FLUID
! flag to use Arrhenius model to calculate reaction rates for particle phase
      LOGICAL :: ARRHENIUS_RRATES_DES
! Arrhenius parameters for each reaction.
      DOUBLE PRECISION,   DIMENSION(:,:), ALLOCATABLE :: arrhenius_coeff

! Tag indicating if calculate the reaction rates for reverse reactions
! based on the forward rate constant and equilbrium constant or the arrhenius coefficients.
! For the reverse reaction rates, there are two methods to calculate the rate constant.
! (1) fromArrheniusCoeff: Similar to the forward reaction rates:
!     Giving the Arrhenius parameters and calculating the rate constant accordingly.
! (2) fromForwardRateConstant: Calculating the forward rate constant Kf and the equilibrium constant Kc.
!     Then the reverse rate constant is equal to Kf/Kc.
      CHARACTER(len=32),  DIMENSION(:),   ALLOCATABLE :: reverse_calc

! species reaction orders when different from its stoichiometric coefficient
      CHARACTER(len=512), DIMENSION(:),   ALLOCATABLE :: rxn_order_array

! Tag indicating if the third body is the mixture or a specific species.
! It can be "M_all" or "M_species name".
! A specific species can be used as third body in the fall-off region
! when the reaction rates depends on pressure
      CHARACTER(len=32),  DIMENSION(:),   ALLOCATABLE :: third_body_model
! Third body coefficients. It can be different from 1
      DOUBLE PRECISION,   DIMENSION(:,:), ALLOCATABLE :: third_body_coeff

! Model used to calculate the reaction rates when reaction is pressure-depedent
! it can be "Lindemann_falloff", "Troe_falloff" or "SRI_falloff"
! or "Lindemann_bimo", "Troe_bimo" or "SRI_bimo",
! or "PLOG".
! "falloff" or "bimo" represent if the reaction is unimolecular/recombination fall-off reactions
! or chemically activated bimolecular reactions. This will determine if the arrhenius_coeff and
! arrhenius_coeff_pres are for high pressure or low pressure.
! For unimolecular/recombination fall-off reactions, arrhenius_coeff is for high-pressure limit
! and arrhenius_coeff_pres is for low-pressure limit.
! For chemically activated bimolecular reactions, arrhenius_coeff is for low-pressure limit
! and arrhenius_coeff_pres is for high-pressure limit.
! For "PLOG" model, several pressures (saved in press_rxn_coeff) and the corresponding
! arrhenius coefficients (saved in arrhenius_coeff_press) are required.
      CHARACTER(len=32),  DIMENSION(:),   ALLOCATABLE :: press_rxn_model
! Coefficients for the pressure-depedent reactions
      DOUBLE PRECISION,   DIMENSION(:,:), ALLOCATABLE :: press_rxn_coeff
! Arrhenius parameters for high (chemically activated bimolecilar reactions)
! or low (fall-off region) pressure, or plog models.
      DOUBLE PRECISION,   DIMENSION(:,:,:), ALLOCATABLE :: arrhenius_coeff_press

! Coefficients for Landau_Teller reactions, for which Arrhenius form cannot
! be used for the rate expression
      DOUBLE PRECISION,   DIMENSION(:,:), ALLOCATABLE :: LT_coeff

! Model to specify the optional rate-constant fit expression
! it can be "Jan" or "Fit1"
      CHARACTER(len=32),  DIMENSION(:),   ALLOCATABLE :: rate_fit_model
! Coefficients for rate fit expressions
      DOUBLE PRECISION,   DIMENSION(:,:), ALLOCATABLE :: rate_fit_coeff

! Logical indicating that the start of a reaction construct has
! been identified.
      LOGICAL IN_DES_CONSTRUCT

! Reaction names
      CHARACTER(len=32),  DIMENSION(:),   ALLOCATABLE :: DES_RXN_NAME
! chemical Equations
      CHARACTER(len=512), DIMENSION(:),   ALLOCATABLE :: DES_RXN_CHEM_EQ
! User defined heat of reaction
      DOUBLE PRECISION,   DIMENSION(:),   ALLOCATABLE :: DES_usrDH
! User defined heat of reaction partitions.
      DOUBLE PRECISION,   DIMENSION(:,:), ALLOCATABLE :: DES_usrfDH
! Arrhenius parameters for each reaction.
      DOUBLE PRECISION,   DIMENSION(:,:), ALLOCATABLE :: DES_arrhenius_coeff
! Tag indicate if calculate the reaction rates
! based on the forward rate constant and equilbrium constant or the arrhenius coefficients
! For the reverse reaction rates, there are two methods to calculate the rate constant.
! (1) fromArrheniusCoeff: Similar to the forward reaction rates:
!     Giving the Arrhenius parameters and calculating the rate constant accordingly.
! (2) fromForwardRateConstant: Calculating the forward rate constant Kf and the equilibrium constant Kc.
!     Then the reverse rate constant is equal to Kf/Kc.
      CHARACTER(len=32),  DIMENSION(:),   ALLOCATABLE :: DES_reverse_calc
! species reaction orders when different from its stoichiometric coefficient
      CHARACTER(len=512), DIMENSION(:),   ALLOCATABLE :: DES_rxn_order_array

      CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: setReaction                                          !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE setReaction(RxN, lNg, lSAg, lM, lNs, lSAs, lDH, lfDH, lReverse, lrxnOrder)

      use rxn_com
      use toleranc

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
! Data structure for storing reaction data.
      TYPE(REACTION_BLOCK), POINTER, INTENT(INOUT) :: RxN
! Number of gas species
      INTEGER, INTENT(IN) :: lNg
! Gas phase species aliases
      CHARACTER(len=32), DIMENSION(DIM_N_g), INTENT(IN) :: lSAg
! Number of solids phases
      INTEGER, INTENT(IN) :: lM
! Number of species in each solids phase.
      INTEGER, DIMENSION(DIM_M), INTENT(IN) :: lNs
! Solids phase species aliases.
      CHARACTER(len=32), DIMENSION(DIM_M, DIM_N_s), INTENT(IN) :: lSAs
! User defined heat of reaction.
      DOUBLE PRECISION, INTENT(IN) :: lDH
! User defined heat of reaction partition.
      DOUBLE PRECISION, DIMENSION(0:DIM_M), INTENT(IN) :: lfDH
! Tag indicating how to calculate reaction rate for reverse reaction
      CHARACTER(len=*), INTENT(IN) :: lReverse
! User defined reaction orders of each species
      CHARACTER(len=*), INTENT(IN) :: lrxnOrder

! Local Variables:
!---------------------------------------------------------------------//
! Alias, phase, species, stoich coeff :: Reactants    Products
      CHARACTER(LEN=32), DIMENSION(50)     :: rAlias    ,  pAlias
      DOUBLE PRECISION, DIMENSION(50) :: rCoeff    ,  pCoeff

! Number of products and reactants
      INTEGER rNo, pNo
! Positions in ChemEq delineating the reactants and products.
      INTEGER rEnd, pStart
! Loop counters
      INTEGER L, LL, M, lN

! Sum of user specified heat of reaction partitions. If fracDH is set
! by the user, they must sum to one over all phases.
      DOUBLE PRECISION sumFDH
! Local storage for chemical equations, left-adjusted and trimmed
      CHARACTER(LEN=512) lChemEq
! Local storage for reaction name, left-adjusted and trimmed
      CHARACTER(LEN=32)  lName
! Logical indicating the reaction is skipped.
      LOGICAL Skip

      LOGICAL pMap(0:lM)

      INTEGER nSpecies, nPhases

      LOGICAL blankAlias(0:(DIM_N_g + lM*DIM_N_s))

! Initialize local reaction name and chemical equation variables.
      lName = trim(adjustl(RxN%Name))
      lChemEq = trim(adjustl(RxN%ChemEq))

      RxN%Classification = "Undefined"
      RxN%Calc_DH = .TRUE.
      RxN%nSpecies = 0
      RxN%nPhases = 0
! Verify that the reactants are separated by --> or = signs. If the
! chemical equation is NONE, the reaction is skipped.
      CALL checkSplit(lName, lChemEq, rEnd, pStart, Skip)
      IF(Skip) THEN
         RxN%nSpecies = 0
         RETURN
      ENDIF
! Set the flag to calculate heat of reaction.
      RxN%Calc_DH = .TRUE.
      IF(lDH /= UNDEFINED) RxN%Calc_DH = .FALSE.

! Pull off the reactants from the chemical equations.
      CALL splitEntries(lName, lChemEq, 1, rEnd, rNo, rAlias, rCoeff)
! Pull off the products from the chemical equations.
      CALL splitEntries(lName, lChemEq, pStart, len_trim(lChemEq),     &
         pNo, pAlias, pCoeff)

      nSpecies = rNo + pNo
      RxN%nSpecies = nSpecies
      Allocate( RxN%Species( nSpecies ))

      CALL checkBlankAliases(lNg, lSAg, lM, lNs, lSAs, blankAlias)

! Check that species in the chemical equation match a species alias
! in one of the phases.

      CALL mapAliases(lName, lChemEq, lNg, lSAg, lM, lNs, lSAs, rNo,   &
         rAlias, rCoeff, -ONE, 0, blankAlias, lrxnOrder, RxN)

! Check that species in the chemical equation match a species alias
! in one of the phases.
      CALL mapAliases(lName, lChemEq, lNg, lSAg, lM, lNs, lSAs, pNo,   &
         pAlias, pCoeff, ONE, rNo, blankAlias, lrxnOrder, RxN)

! All the relevant data has been collected at this point. Build the
! reaction block data structure.
      L = max(1,lM)
      LL = (L * (L-1)/2)
      Allocate( RxN%rPhase( LL+L ))


! Initialize local map and global values
      pMap(:) = .FALSE.
      nPhases = 0
      DO lN = 1, nSpecies
         M = RxN%Species(lN)%pMap

         RxN%Species(lN)%mXfr = M
         RxN%Species(lN)%xXfr = ZERO

         IF(.NOT.pMap(M)) THEN
            pMap(M) = .TRUE.
            nPhases = nPhases + 1
         ENDIF
      ENDDO
      RxN%nPhases = nPhases

! Initialize sum of heat of reaction partitions.
      sumFDH = ZERO
! The user specified the heat of reaction.
      IF(.NOT.RxN%Calc_DH) THEN
! Allocate and initialize the heat of reaction storage array.
         Allocate( RxN%HoR( 0:lM ))
         RxN%HoR(:) = ZERO
         DO M=0,lM
! The phase is referenced by the reaction and heat of reaction is
! allocated (in part or fully) this this phase.
            IF(pMap(M) .AND. lFDH(M) .NE. UNDEFINED) THEN
! Store the heat of reaction.
               RxN%HoR(M) = lFDH(M) * lDH
               sumFDH = sumFDH + lFDH(M)
! The phase is not referenced by the reaction, but the heat of reaction
! is allocated (in part or fully) to this phase. Flag error and exit.
            ELSEIF(.NOT.pMap(M) .AND. lFDH(M) .NE. UNDEFINED) THEN
                  write(err_msg, "(/1X,70('*')/' From: setReaction:',/       &
                     ' Message: Heat of reaction is proportioned to a phase not',  &
                     ' referenced',/' by the chemical equation for reaction ',A,'.', &
                     /' If this is a catalytic reaction, reference one of the', &
                     ' species of the',/' catalyst phase within the chemical',  &
                     ' equation with a stoichiometric',/' coefficient of zero.')") trim(lName)
               call log_error()
            ENDIF
         ENDDO
! Logical check: No partition was assigned to an undefined phase.
         DO M=lM+1,DIM_M
            IF(.NOT.RxN%Calc_DH .AND. lFDH(M) .NE. UNDEFINED) THEN
               write(err_msg, "(/1X,70('*')/' From: From: setReaction:',/       &
                  ' Message: Heat of reaction is proportioned to a phase not',  &
                  ' referenced',/' by the chemical equation for reaction ',A,'.', &
                  /' If this is a catalytic reaction, reference one of the', &
                  ' species of the',/' catalyst phase within the chemical', &
                  ' equation with a stoichiometric',/' coefficient of zero.')") trim(lName)
               call log_error()
            ENDIF
         ENDDO
      ENDIF

! Verify that the heat of reaction partitions sum to one.
      IF(.NOT.RxN%Calc_DH .AND. .NOT. COMPARE(sumFDH, ONE)) THEN
            write(ERR_MSG,"(/1X,70('*')/' From: From: setReaction:',/   &
               ' Message: The heat of reaction partitions (fracDH) to all',  &
               ' phases do',/' not sum to 1 for reaction ',A,'.')") trim(lName)
            call log_error()
      ENDIF

      RETURN

      END SUBROUTINE setReaction


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: checkSplit                                           !
!                                                                      !
!  Purpose: Determine the location of reactatns and products within    !
!  the chemical equation. If the entry is NONE, flag that the reaction !
!  is to be skipped for further processing.                            !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE checkSplit(lName, lChemEq, lrEnd, lpStart, lSkip)

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
! Chemical reaction name.
      CHARACTER(len=*), INTENT(IN) :: lName
! Chemical equation from input file.
      CHARACTER(len=*), INTENT(IN) :: lChemEq
! Position specifying the end of the reactants in lChemEq
      INTEGER, INTENT(OUT) :: lrEnd
! Position specifying the start of the products in lChemEq
      INTEGER, INTENT(OUT) :: lpStart
! If the chemical equation is NONE, the the skip flag is set to avoid
! further processing.
      LOGICAL, INTENT(OUT) :: lSkip

! Local Variables:
!---------------------------------------------------------------------//
! Position of the head (>) and tail (-) of an arror (-->)
      INTEGER hArr, tArr
! Position of the head and tail of equal signs (=, ==, ===, ...)
      INTEGER hEqs, tEqs
! Position of the head/tail of a reverse arrow (<--)
      INTEGER hRArr, tRArr
! A flag generated to point out the location of the entry error.
      CHARACTER(LEN=512) FLAG
      FLAG = ''

! If the chemical equation is set to 'none', then the reaction is
! skipped for the simulation.
      lSkip = .FALSE.
      IF(INDEX(lChemEq,'NONE') > 0) THEN
         lSkip = .TRUE.
         lrEnd = UNDEFINED_I
         lpStart = UNDEFINED_I
         RETURN
      ENDIF

! Search for > (part of -->) and search for < (part of <--)
      tArr = INDEX(lChemEq,'-', BACK=.FALSE.)
      hArr = INDEX(lChemEq,">", BACK=.TRUE.)
! Search for the first and last instances of equal signs.
      tEqs = INDEX(lChemEq,"=", BACK=.FALSE.)
      hEqs = INDEX(lChemEq,"=", BACK=.TRUE.)
! Search for < (as part of <-- or <-->). Illegal chemical equation.
      hRArr = INDEX(lChemEq,"<", BACK=.FALSE.)
      tRArr = INDEX(lChemEq,"-", BACK=.TRUE.)

! An illegal arrow was found. Flag error and exit.
      IF(hRArr > 0) THEN
! Construct the error flag.
            IF(hArr > 0) THEN
               FLAG = setFlag(20, hRArr, hArr)
            ELSEIF(tRArr > 0) THEN
               FLAG = setFlag(20, hRArr, tRArr)
            ELSE
               FLAG = setFlag(20, hRArr)
            ENDIF
            write(ERR_MSG,"(/1X,70('*')/' From: From: setReaction --> checkSplit',/   &
               ' Message: Error in determining the reactants and products',  &
               ' in the',/' chemical equation for reaction ',A,'.', &
               /1X,A,' operators were found.', &
               /' Chemical equation: ',A,/1X, A)") trim(lName), 'Illegal', trim(lChemEq), trim(Flag)
            call log_error()
      ENDIF
! If there are more than one operator, flag error and exit.
      IF(hArr /= 0 .AND. hEqs /= 0) THEN
            FLAG = setFlag(20, hArr, hEqs)
            write(ERR_MSG, "(/1X,70('*')/' From: setReaction --> checkSplit',/   &
               ' Message: Error in determining the reactants and products',  &
               ' in the',/' chemical equation for reaction ',A,'.', &
               /1X,A,' operators were found.', &
               /' Chemical equation: ',A,/1X, A)") trim(lName), 'Too many', trim(lChemEq), trim(Flag)
            call log_error()
! If there is no operator (--> or =), flag error and exit.
      ELSEIF(hArr == 0 .AND. hEqs == 0) THEN
            write(err_msg, "(/1X,70('*')/'From: setReaction --> checkSplit',/   &
               ' Message: Error in determining the reactants and products',  &
               ' in the',/' chemical equation for reaction ',A,'.', &
               /1X,A,' operators were found.', &
               /' Chemical equation: ',A)") trim(lName), 'No', trim(lChemEq)
            call log_error()
! The head of an arrow was found.
      ELSEIF(hArr /= 0) THEN
! Verify that a tail was found.
         IF(tArr == 0) THEN
! Construct the error flag.
            FLAG = setFlag(20, hArr)
               write(err_msg,"(/1X,70('*')/' From: From: setReaction --> checkSplit',/   &
                  ' Message: Error in determining the reactants and products',  &
                  ' in the',/' chemical equation for reaction ',A,'.', &
                  ' Incorrect operator format! ',A, &
                  /' Chemical equation: ',A,/1X, A)") trim(lName), 'Missing the tail; -->', trim(lChemEq), trim(Flag)
               call log_error()
         ELSEIF(tArr > hArr) THEN
               FLAG = setFlag(20, hArr, INDEX(lChemEq,'-',BACK=.TRUE.))
               write(err_msg,"(/1X,70('*')/' From: From: setReaction --> checkSplit',/   &
                  ' Message: Error in determining the reactants and products',  &
                  ' in the',/' chemical equation for reaction ',A,'.', &
                  ' Incorrect operator format  ',A, &
                  /' Chemical equation: ',A,/1X, A/)") trim(lName), 'Arror head precedes the tail; -->', trim(lChemEq), trim(Flag)
               call log_error()
         ELSE
! An arror was used to separate reactants and products. Send back the
! ending index of reactants and the starting index for products.
            lrEnd = tArr - 1
            lpStart = hArr + 1
         ENDIF
! Equals sign(s) were used to specify the reaction. Send back the ending
! index of reactants and the starting index for products.
      ELSEIF(hEqs /= 0) THEN
         lrEnd = tEqs - 1
         lpStart = hEqs + 1
! Fatal Error. One of the above checks should have caught any problems
! and sent out an error message.
      ELSE

         write(err_msg,"(/1X,70('*')/' From: From: setReaction --> checkSplit',/  &
            ' Message: Error in determining the reactants and products',  &
            ' in the',/' chemical equation for reaction ',A,'.', &
            ' FATAL ERROR: All logical checks failed.', &
            /' Chemical equation: ',A,/1X, A/)") trim(lName), trim(lChemEq), trim(Flag)

            call log_error()
      ENDIF

      RETURN

      END SUBROUTINE checkSplit


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: splitEntries()                                       !
!                                                                      !
!  Purpose: Takes a string of either reactants or products and splits  !
!  the string into individual species/stoichiometric entries.          !
!                                                                      !
!  A call to splitAliasAndCoeff is made to further split the entries   !
!  into species aliases and matching stoichiometric coefficients.      !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE splitEntries(lName, lChemEq, lStart, lEnd, lNo,       &
         lAlias, lCoeff)

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
! Chemical reaction name.
      CHARACTER(len=*), INTENT(IN) :: lName
! Chemical equation.
      CHARACTER(len=*), INTENT(IN) :: lChemEq
! Starting position for substring analysis.
      INTEGER, INTENT(IN) :: lStart
! Ending position for substring analysis.
      INTEGER, INTENT(IN) :: lEnd
! The number of individual species found in lSpecies.
      INTEGER, INTENT(OUT) :: lNo
! Species Aliases from the chemical equation.
      CHARACTER(LEN=32), DIMENSION(50), INTENT(OUT) :: lAlias
! Stoichiometric coefficient pulled from the chemical equation.
      DOUBLE PRECISION, DIMENSION(50), INTENT(OUT) :: lCoeff

! Local Variables:
!---------------------------------------------------------------------//
! Flag indicating that there are more entries to process.
      LOGICAL MORE
! Starting position for left-to-right search.
      INTEGER lPOS
! Position of plus sign character found in search.
      INTEGER rPOS

! Initialize storage variables.
      lNo = 0
      lAlias(:) = ''
      lCoeff(:) = UNDEFINED

! Initialize local variables.
      lPOS = lStart
      MORE = .TRUE.
! Loop through the string, splitting entries separated by a plus sign.
      DO WHILE(MORE)
! Increment the species counter.
         lNo = lNo + 1
! Locate the plus sign. (Left to right)
         rPOS = (lPOS-1) + INDEX(lChemEq(lPOS:lEnd),"+", BACK=.FALSE.)
! A plus sign was found.
         IF(rPOS .GT. lPOS) THEN
! Extract the entry and split it into the species alias and
! stoichiometric coefficient.
            CALL splitAliasAndCoeff(lName, lChemEq, lPOS, rPOS-1,      &
               lAlias(lNo), lCoeff(lNo))
! Indicate that there are more entries to process.
            MORE = .TRUE.
! No plus sign was found. This is the last entry.
         ELSE
! Extract the entry and split it into the species alias and
! stoichiometric coefficient.
            CALL splitAliasAndCoeff(lName, lChemEq, lPOS, lEnd,        &
               lAlias(lNo), lCoeff(lNo))
! Indicate that there are no more entries to process.
            MORE = .FALSE.
         ENDIF
! Move past the found plus sign for next search.
         lPOS = rPOS + 1
      ENDDO

      RETURN
      END SUBROUTINE splitEntries

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: splitAliasAndCoeff()                                 !
!                                                                      !
!  Purpose: Take a string containing a species alias and stoichio-     !
!  metric coefficient and splits them into their respective parts.     !
!                                                                      !
!  If no numerical coefficient is found, it is set to one.             !
!                                                                      !
!  If present, asterisks (*) are assumed to separate numerical         !
!  coefficients and the species alias. If more than one asterisk is    !
!  found, and error is reported and MFIX_EXIT is called.               !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE splitAliasAndCoeff(lName, lChemEq, lStart, lEnd,      &
         lAlias, lCoeff)

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
! Chemical reaction name.
      CHARACTER(len=*), INTENT(IN) :: lName
! Chemical equation.
      CHARACTER(len=*), INTENT(IN) :: lChemEq
! Starting position for substring analysis.
      INTEGER, INTENT(IN) :: lStart
! Ending position for substring analysis.
      INTEGER, INTENT(IN) :: lEnd
! Species Aliases from the chemical equation.
      CHARACTER(LEN=32), INTENT(OUT) :: lAlias
! Stoichiometric coefficient pulled from the chemical equation.
      DOUBLE PRECISION, INTENT(OUT) :: lCoeff

! Local Variables:
!---------------------------------------------------------------------//
! Flag
      LOGICAL MATCH
      INTEGER nPOS

      INTEGER L, N, IOS, aPOS, a2POS

      CHARACTER(LEN=12), PARAMETER :: Numbers = '.0123456789'

! A flag generated to point out the location of the entry error.
      CHARACTER(LEN=512) FLAG
      FLAG = ''

! Locate the first asterisk (if any). Search left-to-right.
      aPOS = INDEX(lChemEq(lStart:lEnd),"*", BACK=.FALSE.)
! An asterisk was found.
      IF(aPOS .GT. ZERO) THEN
! Make sure that there isn't another asterisk further down the string.
         a2POS = INDEX(lChemEq(lStart:lEnd),"*", BACK=.TRUE.)
         IF(aPOS /= a2POS) THEN
               FLAG = setFlag(20, lStart+aPOS, lStart+a2POS)
               write(err_msg,"(/1X,70('*')/' From: From: setReaction -->',               &
                  ' splitAliasAndCoeff',/' Message: Error determining the',     &
                  ' stoichiometric coefficient in the',/' chemical equation',   &
                  ' for reaction ',A,'.', &
                  /1X,A,' operators were found.', &
                  /' Chemical equation: ',A,/1X, A/)") trim(lName), 'Too many', trim(lChemEq), trim(Flag)

               call log_error()
         ELSE
! Store left-of-asterisk as the coefficient. If an error occurs in
! converting the string to double precision, flag the problem and
! call MFIX_EXIT.
            READ(lChemEq(lStart:(lStart+aPOS-2)),*,IOSTAT=IOS) lCoeff
            IF(IOS .NE. 0) THEN
! Construct the error flag.
               FLAG = setFlag(20, lStart + int(aPOS/2))
               write(err_msg,"(/1X,70('*')/' From: From: setReaction -->',      &
                  ' splitAliasAndCoeff',/' Message: Error determining the',     &
                  ' stoichiometric coefficient in the',/' chemical equation',   &
                  ' for reaction ',A,'.', /' Chemical equation: ',A,/1X, A)") trim(lName), trim(lChemEq), trim(Flag)
               call log_error()
            ENDIF
! Store right-of-asterisk as the species alias.
            WRITE(lAlias,"(A)") &
               trim(adjustl(lChemEq((lStart+aPOS):lEnd)))
         ENDIF
! If no asterisk was found, search for numbers and spaces.
      ELSE
! Initialize the position of last consecutive number.
         nPOS = 0
! In a left-to-right search, check if the characters in the entry are
! numbers or punctuation.
         DO L=lStart,lEnd
            MATCH = .FALSE.
            DO N=1,12
               IF(lChemEq(L:L) /= Numbers(N:N)) CYCLE
! Note the position of the number.
               nPOS = L
! Flag that a match was made.
               MATCH = .TRUE.
            ENDDO
! If no match, assume the end of the coefficient was found.
            IF(.NOT.MATCH) EXIT
         ENDDO
! If no numbers or punctuation was found, assumed the stoichiometric
! coefficient is one.
         IF(trim(lChemEq(lStart:nPOS)) =='') THEN
            lCoeff = 1.0d0
         ELSE
! If leading numbers were found, store as the stoich-coeff.
            READ(lChemEq(lStart:nPOS),*,IOSTAT=IOS) lCoeff
! Report any problems in converting the string to double precision.
            IF(IOS .NE. 0) THEN
! Construct the error flag.
               FLAG = setFlag(20, &
                  lStart+int(len_trim(lChemEq(lStart:nPOS))/2))
               write(err_msg,"(/1X,70('*')/' From: From: setReaction -->',      &
                  ' splitAliasAndCoeff',/' Message: Error determining the',     &
                  ' stoichiometric coefficient in the',/' chemical equation',   &
                  ' for reaction ',A,'.', /' Chemical equation: ',A,/1X, A)") trim(lName), trim(lChemEq), trim(Flag)
               call log_error()
            ENDIF
         ENDIF
! Store right-of-coefficient as the species alias.
         READ(lChemEq(nPOS+1:lEnd),*,IOSTAT=IOS) lAlias
      ENDIF
! Quick check to make sure that the species alias is not empty.
      IF(LEN_TRIM(lAlias) == 0) THEN
! Construct the error flag.
         FLAG = setFlag(20, lStart + int(lEnd/2))
         write(err_msg,"(/1X,70('*')/' From: From: setReaction -->',      &
            ' splitAliasAndCoeff',/' Message: Error determining the',     &
            ' species in the chemical equation for',/' reaction ',A,'.'/, &
            /' Chemical equation: ',A,/1X, A)") trim(lName), trim(lChemEq), trim(Flag)

         call log_error()
      ENDIF

      RETURN

      END SUBROUTINE splitAliasAndCoeff

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: checkBlankAliases()                                  !
!                                                                      !
!  Purpose: Take a string containing a species alias and stoichio-     !
!  metric coefficient and splits them into their respective parts.     !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE checkBlankAliases(lNg, lSAg, lM, lNs, lSAs, lBA)

      IMPLICIT NONE

! Number of gas species
      INTEGER, INTENT(IN) :: lNg
! Gas phase species aliases
      CHARACTER(len=32), DIMENSION(DIM_N_g), INTENT(IN) :: lSAg
! Number of solids phases
      INTEGER, INTENT(IN) :: lM
! Number of species in each solids phase.
      INTEGER, DIMENSION(DIM_M), INTENT(IN) :: lNs
! Solids phase species aliases.
      CHARACTER(len=32), DIMENSION(DIM_M, DIM_N_s), INTENT(IN) :: lSAs

      LOGICAL, INTENT(OUT) :: lBA(0:(DIM_N_g + lM*DIM_N_s))

      INTEGER M, N

! Loop counter for continuum and discrete species
      INTEGER Nsp

! Initialize counters
      Nsp = 0

      lBA(0) = .FALSE.
      DO N = 1, lNg
         Nsp = Nsp + 1
         lBA(Nsp) = .FALSE.
         IF(len_trim(lSAg(N)) == 0) THEN
            lBA(Nsp) = .TRUE.
            lBA(0) = .TRUE.
         ENDIF
      ENDDO

! Compare aliaes between solids phases
      DO M = 1, lM
         DO N = 1, lNs(M)
            Nsp = Nsp + 1
            lBA(Nsp) = .FALSE.
            IF(len_trim(lSAs(M,N)) == 0) THEN
               lBA(Nsp) = .TRUE.
               lBA(0) = .TRUE.
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE checkBlankAliases

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: mapAliases()                                         !
!                                                                      !
!  Purpose: Take a string containing a species alias and stoichio-     !
!  metric coefficient and splits them into their respective parts.     !
!                                                                      !
!  If no numerical coefficient is found, it is set to one.             !
!                                                                      !
!  If present, asterisks (*) are assumed to separate numerical         !
!  coefficients and the species alias. If more than one asterisk is    !
!  found, and error is reported and MFIX_EXIT is called.               !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE mapAliases(lName, lChemEq, lNg, lSAg, lM, lNs, lSAs,  &
         lNo, lAlias, lCoeff, lSgn, lStart, lBA, lOrder, lRxN)

      use rxn_com

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
! Chemical reaction name.
      CHARACTER(len=*), INTENT(IN) :: lName
! Chemical equation.
      CHARACTER(len=*), INTENT(IN) :: lChemEq
! Number of gas species
      INTEGER, INTENT(IN) :: lNg
! Gas phase species aliases
      CHARACTER(len=32), DIMENSION(DIM_N_g), INTENT(IN) :: lSAg
! Number of solids phases
      INTEGER, INTENT(IN) :: lM
! Number of species in each solids phase.
      INTEGER, DIMENSION(DIM_M), INTENT(IN) :: lNs
! Solids phase species aliases.
      CHARACTER(len=32), DIMENSION(DIM_M, DIM_N_s), INTENT(IN) :: lSAs
! Number of products (or reactants)
      INTEGER, INTENT(IN) :: lNo
! Species Alaises pulled from the chemical equation.
      CHARACTER(LEN=32), DIMENSION(50), INTENT(IN) :: lAlias

      DOUBLE PRECISION, DIMENSION(50), INTENT(IN) :: lCoeff

      DOUBLE PRECISION, INTENT(IN) :: lSgn

      INTEGER, INTENT(IN) :: lStart

      LOGICAL, INTENT(IN) :: lBA(0:(DIM_N_g + lM*DIM_N_s))

      CHARACTER(len=*), INTENT(IN) :: lOrder

! Data structure for storing reaction data.
      TYPE(REACTION_BLOCK), POINTER, INTENT(INOUT) :: lRxN


! Local Variables:
!---------------------------------------------------------------------//
! Loop counters.
      INTEGER L, M, NN, CC
! A flag generated to point out the location of the entry error.
      CHARACTER(LEN=512) FLAG
! Location in string to locate error.
      INTEGER lPOS, rPOS

! array of reaction orders of each reactants in the reaction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::rxnOrder_coeff
! array of reactants given for the reaction orders by users
      CHARACTER(len=32), DIMENSION(:), ALLOCATABLE :: rxnOrder_species
! number of species whose reaction order is given by user
      INTEGER num_coeff
      INTEGER index_c(lNo), index_s
      CHARACTER(len=512) :: char_tmp

! Initialize the error flag.
      FLAG = ''

! If lSgn is negative, the species are reactants
      IF((lSgn .LT. 0) .AND. (lOrder .NE. UNDEFINED_C)) THEN
! split the user defined reaction orders
! get the number of reaction orders, which can be less than the number of reactants,
! and the index of colons for each reactants
	 num_coeff = 0
	 DO CC=1, LEN_TRIM(lOrder)
	    IF(lOrder(CC:CC).EQ. ':') THEN
	       num_coeff = num_coeff+1
	       IF(num_coeff .LE. lNo) THEN
	          index_c(num_coeff) = CC
	       ELSE
	          WRITE(ERR_MSG, 1002) num_coeff, lNo, trim(adjustl(lChemEq))
                call log_error()
	       ENDIF
            ENDIF
	 ENDDO

	 ALLOCATE(rxnOrder_species(num_coeff))
	 ALLOCATE(rxnOrder_coeff(num_coeff))

	 DO CC = 1, num_coeff
	    IF(CC .EQ. 1) rxnOrder_species(CC) = trim(adjustl(lOrder(1:index_c(cc)-1)))
	    IF(CC .LT. num_coeff) THEN
	       char_tmp = adjustl(lOrder(index_c(cc)+1:index_c(cc+1)-1))
	       ! get the index of the first space in char_tmp
	       index_s = INDEX(trim(char_tmp), ' ')
	       IF(index_s .EQ. 0) THEN
	       	  WRITE(ERR_MSG, 1003) trim(adjustl(lChemEq))
                  call log_error()
	       ENDIF
	       READ(char_tmp(1:index_s-1), *) rxnOrder_coeff(CC)
	       rxnOrder_species(CC+1) = trim(adjustl(char_tmp(index_s+1:)))
	    ELSE
	       char_tmp = lOrder(index_c(cc)+1:)
	       READ(char_tmp, *) rxnOrder_coeff(CC)
	    ENDIF
	 ENDDO
      ENDIF

! Loop over the number of (reactants/products)
      ALOOP : DO L=1, lNo
! initialize the reaction orders for all species to be its stoichimetric.
! reaction orders of products will not be used
         lRxn%Species(lStart + L)%RxnOrder = lCoeff(L)
! Compare entry with gas phase species.
         DO NN = 1, lNg
            IF( checkMatch(lSAg(NN), lAlias(L))) THEN
               lRxN%Species(lStart + L)%pMap = 0
               lRxN%Species(lStart + L)%sMap = NN
               lRxN%Species(lStart + L)%Coeff = lSgn * lCoeff(L)
               CYCLE ALOOP
            ENDIF
         ENDDO
! Compare entry with solids phase species.
         DO M = 1, lM
            DO NN = 1, lNs(M)
               IF(checkMatch(lSAs(M,NN),lAlias(L))) THEN
                  lRxN%Species(lStart + L)%pMap = M
                  lRxN%Species(lStart + L)%sMap = NN
                  lRxN%Species(lStart + L)%Coeff = lSgn * lCoeff(L)
                  CYCLE ALOOP
               ENDIF
            ENDDO
         ENDDO
! No matching species was located. Flag an error and exit.

         lPOS = INDEX(lChemEq,trim(lAlias(L)), BACK=.FALSE.)
         rPOS = INDEX(lChemEq,trim(lAlias(L)), BACK=.TRUE.)
         FLAG = setFlag(20, 1 + int((lPOS + rPOS)/2))

         write(ERR_MSG,1000) trim(lAlias(L)), trim(lName), trim(lChemEq), trim(Flag)
         IF(lBA(0)) CALL writeBA()
         call log_error()

      ENDDO ALOOP

! reset the reaction order of reactants
      IF(lSgn .LT. 0.0) THEN
         IF(lOrder .NE. UNDEFINED_C) THEN
            DO CC = 1, num_coeff
	       DO L=1, lNo
                  IF( checkMatch(rxnOrder_species(CC), lAlias(L))) THEN
                     lRxn%Species(lStart + L)%RxnOrder = rxnOrder_coeff(CC)
                     EXIT
                  ENDIF
                  IF(L .EQ. lNo) THEN
                     write(ERR_MSG,1004) trim(rxnOrder_species(CC)), trim(lChemEq)
                     call log_error()
                  ENDIF

               ENDDO
	     ENDDO
            ENDIF
      ENDIF

      RETURN

 1000 FORMAT(/1X,/' From: setReaction --> mapAliases',/                &
         ' Message: Unable to match species ',A,' in the chemical',    &
         ' equation for ',/' reaction ',A,'. Chemical equation: ',A,/1X, A)


 1002 FORMAT(/1X, i4, i4, 70('*')/' From: setReaction --> mapAliases',/&
         ' Message: Number of reaction orders given by rxn_oder, ', A, &
	 ' , is greater than the number of reactants, ', A,/         &
	 ' in chemical equation ', A, '.')

 1003 FORMAT(/1X,70('*')/' From:setReaction --> mapAliases',/         &
         ' Error 1003: Space must be added between coefficient and ',  &
         ' next species name in rxn_order for chemical equation ',/    &
         A, '.')

 1004 FORMAT(/1X,70('*')/' From: setReaction --> mapAliases',/   &
         ' Error 1004: Reaction order of, ', A, ' is given, but not a reactant ',/ &
	 ' in chemical equation, ', A, '.')

 1010 FORMAT(/' Chemical equation: ',A,/1X, A/)

      contains

!......................................................................!
!  Function name: checkMatch                                           !
!                                                                      !
!  Purpose: Takes two species aliases as arguments, converts them to   !
!  uppercase and checks if they match.                                 !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!......................................................................!
      LOGICAL FUNCTION checkMatch(lSA, ceSA)

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
      CHARACTER(LEN=32), INTENT(IN) :: lSA, ceSA

! Local Variables:
!---------------------------------------------------------------------//
      CHARACTER(LEN=32) tlSA

! Copy species alias.
      tlSA = lSA
! Remove case sensitivity.
      CALL MAKE_UPPER_CASE (tlSA,32)
! Compare the two strings.
      checkMatch = .FALSE.
      IF(trim(tlSA) == trim(ceSA)) checkMatch = .TRUE.
      RETURN
      END FUNCTION checkMatch

!......................................................................!
!  Function name: updateMap                                            !
!                                                                      !
!  Purpose: Flags that the passed phase is part of the chemical        !
!  reaction. If the phase was not already noted, the number of phases  !
!  in the reaction is increased and the flag set true. Additionally,   !
!  The number of species (either product or reactant) for the phase    !
!  is incremented.                                                     !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!......................................................................!
      SUBROUTINE updateMap(lnP, lpMap, llNoP)

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
      INTEGER, INTENT(INOUT) :: lnP
! Map of phases for the reaction.
      LOGICAL, INTENT(INOUT) :: lpMap
      INTEGER, INTENT(INOUT) :: llNoP

! Local Variables:
!---------------------------------------------------------------------//
! None

! Increment the number of reactants/products this phase has involved in
! the current reaction.
      llNoP = llNoP + 1
! If the phase was already identified, return.
      IF(lpMap) RETURN
! If this is the first time the phase is identififed, set the flag to
! true and increment the total number of phases in the reaction.
      lnP = lnP + 1
      lpMap = .TRUE.

      RETURN
      END SUBROUTINE updateMap


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: writeBA()                                            !
!                                                                      !
!  Purpose: Print out which species were not given aliases.            !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE writeBA

      IMPLICIT NONE

      INTEGER M, N

! Loop counter for continuum and discrete species
      INTEGER Nsp

      write(ERR_MSG,1000)

! Initialize counters
      Nsp = 0

      DO N = 1, lNg
         Nsp = Nsp + 1
         IF(lBA(Nsp)) THEN
            write(ERR_MSG,1001) N
         ENDIF
      ENDDO

! Compare aliaes between solids phases
      DO M = 1, lM
         DO N = 1, lNs(M)
            Nsp = Nsp + 1

            IF(lBA(Nsp)) THEN
               write(ERR_MSG,1002)M, N
            ENDIF
         ENDDO
      ENDDO

      CALL LOG_WARNING()

      RETURN

 1000 FORMAT(' Species aliases were not provided for the following:')
 1001 FORMAT(3X, ' Gas phase species ',I2)
 1002 FORMAT(3X, ' Solid phase ',I2,' specie ',I2)

      END SUBROUTINE writeBA

      END SUBROUTINE mapAliases


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: setFlag                                              !
!                                                                      !
!  Purpose: Creates a flag pointing to a particular string location.   !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      CHARACTER(len=512) FUNCTION setFlag(fill, flg1, flg2) RESULT(OUT)

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
! Number of leading spaces to fill with dashes. This is used to jump
! past any lead-in text.
      INTEGER, INTENT(IN) :: fill
! Location in string to place the pointer.
      INTEGER, INTENT(IN) :: flg1
! Optional - location of second pointer in string
      INTEGER, INTENT(IN), OPTIONAL :: flg2

! Local Variables:
!---------------------------------------------------------------------//
      INTEGER L, FILL1, FILL2

! Create a string with "FILL" dash characters.
      OUT = ''
      DO L = 1, FILL-1
         WRITE(OUT,"(A,A)") trim(OUT), '-'
      ENDDO

! If a second pointer is present, determined the larger of the two.
      IF(PRESENT(flg2)) THEN
         IF(flg1 < flg2) THEN
            FILL1 = flg1 - 1
            FILL2 = (flg2-flg1) - 1
         ELSE
            FILL1 = flg2 - 1
            FILL2 = (flg1-flg2) - 1
         ENDIF
      ELSE
         FILL1 = flg1 - 1
         FILL2 = 0
      ENDIF

! Fill with dashes up to the the first pointer. ----^
      DO L = 1, FILL1
         WRITE(OUT,"(A,A)") trim(OUT), '-'
      ENDDO
      WRITE(OUT,"(A,A)") trim(OUT), '^'
! Fill with dashes up to the second pointer. ----^---^
      IF(FILL2 > 0) THEN
         DO L = 1, FILL2
            WRITE(OUT,"(A,A)") trim(OUT), '-'
         ENDDO
         WRITE(OUT,"(A,A)") trim(OUT), '^'
      ENDIF

      END FUNCTION setFlag

END MODULE parse
