#include "error.inc"

MODULE PARSE_RXN_MOD

   use error_manager

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARSE_RXN(LINE, LMAX)                                  C
!  Purpose: Parse input line                                           C
!                                                                      C
!  Author: P. Nicoletti                               Date: 30-JUN-97  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: This routine was complete rewritten as part of the effort  C
!  to simplify reaction inputs in MFiX.                                C
!  Author: J. Musser                                  Date: 01-Oct-12  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: Add the parameters used to calculate the reaction rates.   C
!           For example, the arrhenius parameters can be given and the C
!  reaction rates will be calculated accordingly without giving UDFs.  C
!  More parameters, like user defined reaction orders, third body      C
!  reactions and pressure-depedent reactions can be given.             C
!  Please refer to CHEMKIN manual for more details.                    C
!  NOTE: Parameters about third-body reactions, pressure-dependent     C
!       reactions, Landau_Teller reactions, and optional rate-constant C
!       fit expression are only used for gas-phase reactons. If one of C
!       them is given, the rest must be given as well.                 C
!  Author: Hang Zhou                                  Date: 23-May-24  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PARSE_RXN(LINE, lNoOfRxns, lName, lChemEq, lDH, lFDH, &
                           lArrhen, lReverse, lrxnOrder, lMmodel, lMcoeff, &
                           lpreModel, lpreCoeff, lArrhenPre, lLTCoeff, &
                           lFitModel, lFitCoeff)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE funits
      USE param
      USE param1
      USE parse

      IMPLICIT NONE

! Input line from input file
      CHARACTER(len=*), INTENT(IN) :: LINE
! Array of reaction names.
      INTEGER, INTENT(INOUT) :: lNoOfRxns
! Array of reaction names.
      CHARACTER(len=*), INTENT(INOUT), DIMENSION(DIMENSION_RXN) :: lName
! Array of Chemical reaction equations.
      CHARACTER(len=*), INTENT(INOUT), DIMENSION(DIMENSION_RXN) :: lChemEq
! Array of User defined heat of reactions.
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(DIMENSION_RXN) :: lDH
! Array of User defined heat of reaction phase partitions.
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(DIMENSION_RXN, 0:DIM_M) :: lFDH
! Array of User defined Arrhenius parameters for each reaction.
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(DIMENSION_RXN, 3) :: lArrhen
! Array about how to calculate the rate constant for reverse reactions.
      CHARACTER(len=*), INTENT(INOUT), DIMENSION(DIMENSION_RXN) :: lReverse
! Array of user defined species reaction order.
      CHARACTER(len=*), INTENT(INOUT), DIMENSION(DIMENSION_RXN) :: lrxnOrder
! Array about species considered as third body.
      CHARACTER(len=*), OPTIONAL, INTENT(INOUT), DIMENSION(DIMENSION_RXN) :: lMmodel
! Array of third body coefficients. It can be different from 1
      DOUBLE PRECISION, OPTIONAL, INTENT(INOUT), DIMENSION(:, :) :: lMcoeff
! Array of models used for presure-dependent reactions
      CHARACTER(len=*), OPTIONAL, INTENT(INOUT), DIMENSION(:) :: lpreModel
! Array of coefficients for the pressure-depedent reactions
      DOUBLE PRECISION, ALLOCATABLE, OPTIONAL, INTENT(INOUT), DIMENSION(:,:) :: lpreCoeff
! Array of user defined Arrhenius parameters for pressure-dependent reactions
      DOUBLE PRECISION, ALLOCATABLE, OPTIONAL, INTENT(INOUT), DIMENSION(:,:,:) :: lArrhenPre
! Array of coefficients for Landau Teller reactions
      DOUBLE PRECISION, OPTIONAL, INTENT(INOUT), DIMENSION(:,:) :: lLTCoeff
! Model to specify the optional rate-constant fit expression
! it can be "Jan" or "Fit1"
      CHARACTER(len=*), OPTIONAL, INTENT(INOUT), DIMENSION(:) :: lFitModel
! Coefficients for rate fit expressions
      DOUBLE PRECISION, OPTIONAL, INTENT(INOUT), DIMENSION(:,:) :: lFitCoeff
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      CHARACTER, PARAMETER :: CT_BEG = '{'
      CHARACTER, PARAMETER :: CT_END = '}'

! Positions of braces {...}
      INTEGER bIDX, eIDX
! Reaction Name
      CHARACTER(LEN=1024) :: INPUT
! Index of reaction.
      INTEGER IDX,ll
! string of the parameters for pressure-related reaction rates
      CHARACTER(LEN=1024) :: lpre_param_str
! string of the coefficients for third-body species
      CHARACTER(LEN=1024) :: lMcoeff_str
! string of the parameters for optional rate-constant fit expression
      CHARACTER(LEN=512) :: lfit_param_str

! Skip this routine for blank lines
      if(len_trim(input) == 0) return

! Copy line to input for processing.
      INPUT = TRIM(ADJUSTL(LINE))

! Look for the start and end of a reaction construct by checking for
! left and right braces.
      bIDX = INDEX(INPUT,CT_BEG)
      eIDX = INDEX(INPUT,CT_END)

! If not already inside/reading from a reaction construct, check to see
! if this is the start of a construct.
      IF(.NOT.IN_CONSTRUCT) THEN

! An identifier for the end of a construct was found.
         IF(eIDX .GT. 0) THEN
            IF(eIDX .GT. bIDX) THEN
! The reaction construct is specified in a single line.
! rxn_name { A + B --> AB }
               IDX = getReactionIndex(lNoOfRxns, 'NEW')
! Pull off the reaction construct name.
               CALL getName(INPUT,(bIDX-1), lNAME(IDX))
! Process the rest of the line.
               IF(isFracDH(INPUT(bIDX+1:eIDX-1)))THEN
                  WRITE(ERR_MSG, 1002) 'FracDH', trim(adjustl(INPUT))
                  call log_error()
               ELSEIF(isDH(INPUT(bIDX+1:eIDX-1)))THEN
                  WRITE (ERR_MSG, 1002) 'DH', trim(adjustl(INPUT))
                  call log_error()
               ELSEIF(isArrhenius(INPUT(bIDX+1:eIDX-1)))THEN
                  WRITE (ERR_MSG, 1002) 'arrhenius_coeff', trim(adjustl(INPUT))
                  call log_error()
               ELSEIF(isReverseCalc(INPUT(bIDX+1:eIDX-1)))THEN
                  WRITE (ERR_MSG, 1002) 'reverse_calc', trim(adjustl(INPUT))
                  call log_error()
               ELSEIF(isRxnOrder(INPUT(bIDX+1:eIDX-1)))THEN
                  WRITE (ERR_MSG, 1002) 'rxn_order', trim(adjustl(INPUT))
                  call log_error()
               ELSEIF(isThirdBody(INPUT(bIDX+1:eIDX-1)))THEN
                  WRITE (ERR_MSG, 1002) 'third_body_param', trim(adjustl(INPUT))
                  call log_error()
               ELSEIF(isPressRxn(INPUT(bIDX+1:eIDX-1)))THEN
                  WRITE (ERR_MSG, 1002) 'press_rxn_param', trim(adjustl(INPUT))
                  call log_error()
               ELSEIF(isLTCoeff(INPUT(bIDX+1:eIDX-1)))THEN
                  WRITE (ERR_MSG, 1002) 'LT_coeff', trim(adjustl(INPUT))
                  call log_error()
               ELSEIF(isFITCoeff(INPUT(bIDX+1:eIDX-1)))THEN
                  WRITE (ERR_MSG, 1002) 'FIT_coeff', trim(adjustl(INPUT))
                  call log_error()
               ELSEIF(.NOT.isChemEq(INPUT(bIDX+1:eIDX-1)))THEN
                  WRITE(ERR_MSG, 1003) trim(adjustl(INPUT))
                  call log_error()
               ENDIF
               CALL getChemEq(INPUT(bIDX+1:eIDX-1), lChemEq(IDX))
            ELSE
! The format given in the input file is incorrect. Brace mismatch.
               WRITE(ERR_MSG, 1001) trim(adjustl(LINE))
               call log_error()
            ENDIF
         ELSE
! This is the start of a reaction construct.
            lpre_param_str = ""
            lfit_param_str = ""
            lMcoeff_str = ""
            IF(bIDX .GT. 0) THEN
! Get the reaction index.
               IDX = getReactionIndex(lNoOfRxns, 'NEW')
! Extract the reaction name.
               CALL getName(INPUT,(bIDX-1), lNAME(IDX))
! Process any data.
               IF(LEN_TRIM(ADJUSTL(INPUT(bIDX+1:eIDX-1))) .GT. 0) THEN
                IF (present(lMmodel)) THEN
                     CALL readConstruct(INPUT(bIDX+1:eIDX-1), lChemEq(IDX), lDH(IDX), lFDH(IDX,:), &
                                        lArrhen(IDX,:), lReverse(IDX), lrxnOrder(IDX),&
                                        lMmodel(IDX), lMcoeff(IDX,:), lMcoeff_str, lpreModel(IDX), lpreCoeff, lArrhenPre, &
                                        lpre_param_str, lLTCoeff(IDX,:), lFitModel(IDX), lFitCoeff(IDX,:), &
                                        lfit_param_str, IDX)
                ELSE
                 CALL readConstruct(INPUT(bIDX+1:eIDX-1), lChemEq(IDX), lDH(IDX), lFDH(IDX,:), &
                                        lArrhen(IDX,:), lReverse(IDX), lrxnOrder(IDX))
              ENDIF
             ENDIF
               IN_CONSTRUCT = .TRUE.
            ELSE
! Format Error.
               WRITE(ERR_MSG, 1004) trim(adjustl(INPUT))
               call log_error()
            ENDIF
         ENDIF
      ELSE

         IF(bIDX .GT. 0) THEN
! Format Error.
            WRITE(ERR_MSG, 1005) trim(adjustl(INPUT))
            call log_error()
! This is the last line of the reaction construct which may or may not
! contain additional data.
         ELSEIF(eIDX .GT. 0) THEN
           IDX = getReactionIndex(lNoOfRxns)
           IF (present(lMmodel)) THEN
              CALL readConstruct(INPUT(bIDX+1:eIDX-1), lChemEq(IDX), lDH(IDX), lFDH(IDX,:), &
                                 lArrhen(IDX,:), lReverse(IDX), lrxnOrder(IDX),&
                                 lMmodel(IDX), lMcoeff(IDX,:), lMcoeff_str, lpreModel(IDX), lpreCoeff, lArrhenPre, &
                                 lpre_param_str, lLTCoeff(IDX,:), lFitModel(IDX), lFitCoeff(IDX,:), &
                                 lfit_param_str, IDX)
         ELSE
            CALL readConstruct(INPUT(bIDX+1:eIDX-1), lChemEq(IDX), lDH(IDX), lFDH(IDX,:), &
                                 lArrhen(IDX,:), lReverse(IDX), lrxnOrder(IDX))
           ENDIF
           IN_CONSTRUCT = .FALSE.

! Reading from somewhere inside of a reaction construct.
         ELSE
            IDX = getReactionIndex(lNoOfRxns)

          IF (present(lMmodel)) THEN
               CALL readConstruct(INPUT(bIDX+1:), lChemEq(IDX), lDH(IDX), lFDH(IDX,:), &
                                  lArrhen(IDX,:), lReverse(IDX), lrxnOrder(IDX), &
                                  lMmodel(IDX), lMcoeff(IDX,:), lMcoeff_str, lpreModel(IDX), lpreCoeff, lArrhenPre, &
                                  lpre_param_str, lLTCoeff(IDX,:), lFitModel(IDX), lFitCoeff(IDX,:), &
                                  lfit_param_str, IDX)
          ELSE
             CALL readConstruct(INPUT(bIDX+1:), lChemEq(IDX), lDH(IDX), lFDH(IDX,:), &
                                  lArrhen(IDX,:), lReverse(IDX), lrxnOrder(IDX))
            ENDIF
         ENDIF
      ENDIF
      RETURN

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN',/                         &
         ' Error 1001: Mismatch of braces "{...}" in reaction ',       &
         ' construct.',/                                               &
     /' INPUT: ',A)

 1002 FORMAT(//1X,70('*')/' From: PARSE_RXN',/                         &
         ' Error 1002: Input format error in reaction construct.',     &
         ' Opening and',/' closing braces were found on the same line',&
         ' along with the',/' keyword ',A,'.',/' Single line',         &
         ' constructs can only contain a chemical equation.',//        &
         ' INPUT: ',A,//                                               &
         ' Example 1: RXN_NAME { chem_eq = "A + B --> AB" }',//        &
         ' Example 2: RXN_NAME {',/14X,'chem_eq = "A + B --> AB"',/14X,&
         'DH = 2.5d4',/14X,'fracDH(0) = 1.0',/12X,'}')

 1003 FORMAT(//1X,70('*')/' From: PARSE_RXN',/                         &
         ' Error 1003: Input format error in reaction construct.',     &
         ' Opening and',/' closing braces were found on the same line',&
         ' and chem_eq was NOT found.',/' Single line constructs can', &
         ' only contain a chemical equation.',//' INPUT: ',A,//        &
         ' Example 1: RXN_NAME { chem_eq = "A + B --> AB" }',//        &
         ' Example 2: RXN_NAME {',/14X,'chem_eq = "A + B --> AB"',/14X,&
         'DH = 2.5d4',/14X,'fracDH(0) = 1.0',/12X,'}')

 1004 FORMAT(//1X,70('*')/' From: PARSE_RXN',/                         &
         ' Error 1004: Data within the reaction block was identified', &
         ' outside of a',/' reaction construct. ',                     &
     /' INPUT: ',A)

 1005 FORMAT(//1X,70('*')/' From: PARSE_RXN',/                         &
         ' Error 1005: The start of a new reaction construct was',     &
         ' found before the',/' closing of the previous construct.',/  &
         /' INPUT: ',A)

      CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: getReactionIndex()                                   !
!                                                                      !
!  Purpose: Extract the reaction name from a construct.                !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      INTEGER FUNCTION getReactionIndex(lNOR, STAT)

      use rxns

      IMPLICIT NONE
! Number of reactions.
      INTEGER, INTENT(INOUT) :: lNOR
! Status
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: STAT

      IF(.NOT.PRESENT(STAT)) THEN
         getReactionIndex = lNOR

      ELSE
         IF(STAT == 'NEW') THEN
! Increment the number of reactions processed from the data file and
! return the new value as the index of the reaction being processed.
             lNOR = lNOR + 1
             if (lnor > DIMENSION_RXN ) then
                 write(ERR_MSG,*) "Error, too many reactions! Maximum allowed: ", DIMENSION_RXN
                 call log_error()
             endif
            getReactionIndex = lNOR
         ELSE
            WRITE(ERR_MSG,*) ' Unknown status'
            call log_error()
         ENDIF
      ENDIF

      RETURN
      END FUNCTION getReactionIndex



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: readConstruct(IN, ChemEq, uDH, uFDH, uArrhen,        !
!                               uReverse, urxnOrder, uMmodel, uMcoeff, !
!                               uMcoeff_str, upreModel, upreCoeff,     !
!                               uArrhenPre, upre_param_str, uLTCoeff,  !
!                               uFitModel, uFitCoeff, ufit_param_str,  !
!                               id_rxn)                                !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Revision Number: 1                                                  !
!  Purpose: Add the parameters used to calculate the reaction rates.   !
!           For example, the arrhenius parameters can be given and the !
!           reaction rates will be calculated accordingly without      !
!           giving UDFs. More parameters, like user defined reaction   !
!           orders, third body reactions and pressure-depedent         !
!           reactions can be given.                                    !
!           Please refer to CHEMKIN manual for more details.           !
!  NOTE: Parameters about third-body reactions, pressure-dependent     !
!        reactions, Landau Teller reactions, and optional rate-constant!
!        fit expression are only used for gas-phase reactons.          !
!        These are optional arguments. You can give ALL or NONE of the !
!        optional arguments. ou cannot give only part of the optional  !
!        arguments.                                                    !
!  Author: Hang Zhou                                  Date: 23-May-24  !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE readConstruct(IN, ChemEq, uDH, uFDH, uArrhen, uReverse, &
                               urxnOrder, uMmodel, uMcoeff, uMcoeff_str, upreModel, &
                               upreCoeff, uArrhenPre, upre_param_str, uLTCoeff, &
                               uFitModel, uFitCoeff, ufit_param_str, id_rxn)

      IMPLICIT NONE

! Input string being parsed.
      CHARACTER(len=*), INTENT(IN) :: IN
! Chemical equation.
      CHARACTER(len=*), INTENT(OUT) :: ChemEq
! User defined heat of reaction.
      DOUBLE PRECISION, INTENT(INOUT) :: uDH
! User defined splitting of heat of reaction
      DOUBLE PRECISION, INTENT(INOUT) :: uFDH(0:DIM_M)
! Array of User defined Arrhenius parameters
      DOUBLE PRECISION, INTENT(INOUT) :: uArrhen(:)
! Array about how to calculate the rate constant for reverse reactions.
      CHARACTER(len=*), INTENT(INOUT) :: uReverse
! Array of user defined species reaction order.
      CHARACTER(len=*), INTENT(INOUT) :: urxnOrder
! Array about species considered as third body.
      CHARACTER(len=*), OPTIONAL, INTENT(INOUT) :: uMmodel
! Array of third body coefficients. It can be different from 1
      DOUBLE PRECISION, OPTIONAL, INTENT(INOUT) :: uMcoeff(:)
! Array of coefficients used for third-body species
      CHARACTER(len=*), OPTIONAL, INTENT(INOUT) :: uMcoeff_str
! Array of models used for presure-dependent reactions
      CHARACTER(len=*), OPTIONAL, INTENT(INOUT) :: upreModel
! Array of coefficients for the pressure-depedent reactions
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: upreCoeff
! Array of user defined Arrhenius parameters for pressure-dependent reactions
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: uArrhenPre
! Array of parameters used for presure-dependent reactions
      CHARACTER(len=*), OPTIONAL, INTENT(INOUT) :: upre_param_str
! Input string being parsed.
      INTEGER, OPTIONAL, INTENT(IN) :: id_rxn
! Array of coefficients for Landau_Teller reactions
      DOUBLE PRECISION, OPTIONAL, INTENT(INOUT) :: uLTCoeff(:)
! Model to specify the optional rate-constant fit expression
! it can be "Jan" or "Fit1"
      CHARACTER(len=*), OPTIONAL, INTENT(INOUT) :: uFitModel
! Coefficients for rate fit expressions
      DOUBLE PRECISION, OPTIONAL, INTENT(INOUT) :: uFitCoeff(:)
! Array of parameters used for presure-dependent reactions
      CHARACTER(len=*), OPTIONAL, INTENT(INOUT) :: ufit_param_str
! The input line contains no additional data.
      IF(LEN_TRIM(ADJUSTL(IN)) == 0) RETURN
! The input contains chemical equation data.
      IF(MORE_ChemEq .OR. isChemEq(IN)) THEN
         CALL getChemEq(IN, ChemEq)
! The input contains heat of reaction parsing data.
      ELSEIF(isFracDH(IN)) THEN
         CALL getFracDH(IN, uFDH(:))
! The input contains heat of reaction data.
      ELSEIF(isDH(IN)) THEN
         CALL getDH(IN, uDH)
! The input contains Arrhenius parameters.
      ELSEIF(isArrhenius(IN)) THEN
         CALL getArrhenius(IN, uArrhen(:))
! The input contains infos about how to calculate the rate constant for reverse reactions.
      ELSEIF(isReverseCalc(IN)) THEN
         CALL getReverseCalc(IN, uReverse)
! The input contains user defined reaction orders.
      ELSEIF(MORE_RxnOrder .OR. isRxnOrder(IN)) THEN
         CALL getRxnOrder(IN, urxnOrder)
! The input contains third body information.
      ELSEIF(MORE_thirdBody .OR. isThirdBody(IN)) THEN
         CALL getThirdBody(IN, uMmodel, uMcoeff, uMcoeff_str)
! The reaction is pressure-dependent.
      ELSEIF(MORE_PressParam .OR. isPressRxn(IN)) THEN
         CALL getPressRxn(IN, upreModel, upreCoeff, uArrhenPre, id_rxn, upre_param_str)
! The input contains coefficients for Landau_Teller reactions
      ELSEIF(isLTCoeff(IN)) THEN
         CALL getLTCoeff(IN, uLTCoeff(:))
! The input contains coefficients for Landau_Teller reactions
      ELSEIF(MORE_fitParam .OR. isFITCoeff(IN)) THEN
         CALL getFITCoeff(IN, uFitModel, uFITCoeff(:), ufit_param_str)
! The entry doesn't match any of the keywords.
      ELSE
! Unidentified keyword.
         WRITE(ERR_MSG, 1001) trim(adjustl(IN))
         call log_error()
      ENDIF

      RETURN

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> readConstruct',/       &
         ' Error 1001: Unidentified keyword in reaction construct.'/,  &
         /' INPUT: ',A)

      END SUBROUTINE readConstruct



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: isChemEq(INPUT)                                      !
!                                                                      !
!  Purpose: Checks if the line contains the chemical Eq.               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION isChemEq(INPUT)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT

! Check to see if the line contains 'END'
      IF (INDEX(LINE(1:),"CHEM_EQ") == 0) THEN
! 'CHEM_EQ' was not found. This line does not contains a chemical eq.
         isChemEq = .FALSE.
      ELSE
! 'CHEM_EQ' was found. This line contains all or part of a chemical eq.
         isChemEq = .TRUE.
      ENDIF

      END FUNCTION isChemEq


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: isDH(INPUT)                                          !
!                                                                      !
!  Purpose: Checks if the line contains user defined heat of reaction. !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION isDH(INPUT)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! Check to see if the line contains 'END'
      IF (INDEX(LINE(1:),"DH") == 0) THEN
! 'DH' was not found. This line does not contains a heat of reaction.
         isDH = .FALSE.
      ELSE
! 'DH' was found. This line contains the heat of reaction
         isDH = .TRUE.
      ENDIF
      END FUNCTION isDH

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: isFracDH(INPUT)                                          !
!                                                                      !
!  Purpose: Checks if the line contains user defined heat of reaction. !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION isFracDH(INPUT)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT

! Check to see if the line contains 'END'
      IF (INDEX(LINE(1:),"FRACDH") == 0) THEN
! 'FRACDH' was not found.
         isFracDH = .FALSE.
      ELSE
! 'FRACDH' was found.
         isFracDH = .TRUE.
      ENDIF

      END FUNCTION isFracDH

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: isArrhenius(INPUT)                                   !
!                                                                      !
!  Purpose: Checks if the line contains user defined arrhenius         !
!           parameters.                                                !
!                                                                      !
!  Author: Hang Zhou                                Date: May 23, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION isArrhenius(INPUT)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! Check to see if the line contains 'END'
      IF (INDEX(LINE(1:),"ARRHENIUS_COEFF") == 0) THEN
! 'arrhenius_coeff' was not found.
         isArrhenius = .FALSE.
      ELSE
! 'arrhenius_coeff' was found.
         isArrhenius = .TRUE.
      ENDIF
      END FUNCTION isArrhenius

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: isReverseCalc(INPUT)                                 !
!                                                                      !
!  Purpose: Checks if the line contains information about the way to   !
!           to calculate reverse reaction rate                         !
!                                                                      !
!  Author: Hang Zhou                                Date: May 23, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION isReverseCalc(INPUT)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT

! Check to see if the line contains 'END'
      IF (INDEX(LINE(1:),"REVERSE_CALC") == 0) THEN
! 'reverse_calc' was not found.
         isReverseCalc = .FALSE.
      ELSE
! 'reverse_calc' was found.
         isReverseCalc = .TRUE.
      ENDIF

      END FUNCTION isReverseCalc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: isRxnOrder(INPUT)                                    !
!                                                                      !
!  Purpose: Checks if the line contains user defined reaction order.   !
!                                                                      !
!  Author: Hang Zhou                                Date: May 23, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION isRxnOrder(INPUT)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT

! Check to see if the line contains 'END'
      IF (INDEX(LINE(1:),"RXN_ORDER") == 0) THEN
! 'rxn_order' was not found.
         isRxnOrder = .FALSE.
      ELSE
! 'rxn_order' was found.
         isRxnOrder = .TRUE.
      ENDIF

      END FUNCTION isRxnOrder

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: isThirdBody(INPUT)                                   !
!                                                                      !
!  Purpose: Checks if the line contains user defined                   !
!           third body coefficients.                                   !
!                                                                      !
!  Author: Hang Zhou                                Date: May 23, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION isThirdBody(INPUT)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT

! Check to see if the line contains 'END'
      IF (INDEX(LINE(1:),"THIRD_BODY_PARAM") == 0) THEN
! 'third_body_param' was not found.
         isThirdBody = .FALSE.
      ELSE
! 'third_body_param' was found.
         isThirdBody = .TRUE.
      ENDIF

      END FUNCTION isThirdBody

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: isPressRxn(INPUT)                                     !
!                                                                      !
!  Purpose: Checks if the line contains user defined                   !
!           coefficients under various pressures.                      !
!                                                                      !
!  Author: Hang Zhou                                Date: May 23, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION isPressRxn(INPUT)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! Check to see if the line contains 'END'
      IF (INDEX(LINE(1:),"PRESS_RXN_PARAM") == 0) THEN
! 'press_rxn_param' was not found.
         isPressRxn = .FALSE.
      ELSE
! 'press_rxn_param' was found.
         isPressRxn = .TRUE.
      ENDIF

      END FUNCTION isPressRxn


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: isLTCoeff(INPUT)                                     !
!                                                                      !
!  Purpose: Checks if the line contains user defined coefficients      !
!           for Landau Teller reactions.                               !
!                                                                      !
!  Author: Hang Zhou                                Date: May 23, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION isLTCoeff(INPUT)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT

! Check to see if the line contains 'END'
      IF (INDEX(LINE(1:),"LT_COEFF") == 0) THEN
! 'LT_coeff' was not found.
         isLTCoeff = .FALSE.
      ELSE
! 'LT_coeff' was found.
         isLTCoeff = .TRUE.
      ENDIF

      END FUNCTION isLTCoeff

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: isFITCoeff(INPUT)                                    !
!                                                                      !
!  Purpose: Checks if the line contains user defined coefficients      !
!           to specify the optional rate-constant fit expression.      !
!                                                                      !
!  Author: Hang Zhou                                Date: May 23, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION isFITCoeff(INPUT)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT

! Check to see if the line contains 'END'
      IF (INDEX(LINE(1:),"FIT_COEFF") == 0) THEN
! 'FIT_coeff' was not found.
         isFITCoeff = .FALSE.
      ELSE
! 'FIT_coeff' was found.
         isFITCoeff = .TRUE.
      ENDIF

      END FUNCTION isFITCoeff

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: get_ChemEq(INPUT, lNAME, IER)                      !
!                                                                      !
!  Purpose: Extract the reaction name from a construct.                !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getName(INPUT, rPOS, lNAME)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! End of search location for reaction name.
      INTEGER, INTENT(IN) :: rPOS
! Name of reaction pulled from input.
      CHARACTER(LEN=32) , INTENT(OUT) :: lNAME

      INTEGER NAME_LEN

! Initialize the return value.
      lNAME = ''
! Verify that the name is not too long. This should be caught by
! preprocessing of the data file. However, if the user changed the
! reaction name after compiling (an error check for later) this check
! prevents and overflow error.
      NAME_LEN = len_trim(adjustl(INPUT(1:rPOS)))
      IF(NAME_LEN .GT. 32) THEN
         WRITE(ERR_MSG, 1001) trim(adjustl(INPUT))
         call log_error()
! Verify that the name was not deleted after compiling.
! prevents and overflow error.
      ELSEIF(NAME_LEN .EQ. 0) THEN
         WRITE(ERR_MSG, 1002) trim(adjustl(INPUT))
         call log_error()
      ELSE
         lNAME = trim(adjustl(INPUT(1:rPOS)))
      ENDIF

! There shouldn't be any crazy characters at this point because the
! code should fail to compile if the reaction names are not defined
! or contain invalid characters.

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> get_ChemEq',/          &
         ' Error 1001: Reaction name too long! Reaction names are',    &
         ' limited to 32',/' characters.',/ &
        /' Reaction Name: ',A)

 1002 FORMAT(//1X,70('*')/' From: PARSE_RXN --> get_ChemEq',/          &
         ' Error 1002: Unable to determine reaction name.',/          &
         /' INPUT: ',A)

      RETURN
      END SUBROUTINE getName


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: getDH(INPUT, lDH)                                  !
!                                                                      !
!  Purpose: Extract the reaction name from a construct.                !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getDH(INPUT, lDH)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! Name of reaction pulled from input.
      DOUBLE PRECISION, INTENT(OUT) :: lDH

      INTEGER lQ
      INTEGER lLMAX
! read/write output status
      INTEGER IOS

      lLMAX = LEN_TRIM(INPUT)

      IF(INDEX(INPUT,"DH") .EQ. 0) THEN
         WRITE (ERR_MSG, 1100) trim(adjustl(INPUT))
         call log_error()
      ENDIF

      lQ = INDEX(INPUT(:lLMAX),'=')

      IF(lQ .EQ. 0) THEN
         WRITE (ERR_MSG, 1001) trim(adjustl(INPUT))
         call log_error()
      ENDIF

! Convert the entrying into an double precision value.
      READ(INPUT(lQ+1:),*,IOSTAT=IOS) lDH
      IF(IOS .NE. 0) THEN
         WRITE(ERR_MSG, 1002) trim(adjustl(INPUT))
         call log_error()
      ENDIF


 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getDH',/               &
         ' Error 1001: Input format error for DH.',/                   &
        /' INPUT: ',A)

 1002 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getDH',/               &
         ' Error 1002: Unable to determine DH value from input.',/     &
         ' Cannot convert specified value to double precision value.',/&
         /' INPUT: ',A)

 1100 FORMAT(//1X,70('*')/' From: PARSE_RXN',/                         &
         ' Error 1105: DH was initially located within the input line',&
         /' however its location cannot be determined.',/&
         /' INPUT: ',A)

      RETURN
      END SUBROUTINE getDH


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: getfracDH(INPUT, lChemEq)                          !
!                                                                      !
!  Purpose: Extract the reaction name from a construct.                !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getFracDH(INPUT, lFracDH)

      USE param
      USE param1

      IMPLICIT NONE

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! Name of reaction pulled from input.
      DOUBLE PRECISION, INTENT(OUT) :: lFracDH(0:DIM_M)


      INTEGER POS, lP, rP, lQ
      INTEGER lLMAX
! read/write output status
      INTEGER IOS
! Phase Index
      INTEGER pIDX

      lLMAX = LEN_TRIM(INPUT)
      POS = INDEX(INPUT,"FRACDH")

      IF(POS == 0) THEN
         WRITE (ERR_MSG, 1100) trim(adjustl(INPUT))
         call log_error()
      ENDIF

      lP = INDEX(INPUT(:lLMAX),'(')
      rP = INDEX(INPUT(:lLMAX),')')
      lQ = INDEX(INPUT(:lLMAX),'=')

      IF(lP .EQ. rP .AND. lP .EQ. ZERO) THEN
         WRITE(ERR_MSG, 1001) trim(adjustl(INPUT))
         call log_error()
      ELSEIF(lP .GE. rP) THEN
         WRITE(ERR_MSG, 1002) trim(adjustl(INPUT))
         call log_error()
      ELSEIF(rP .GE. lQ) THEN
         WRITE(ERR_MSG, 1002) trim(adjustl(INPUT))
         call log_error()
      ENDIF
! Convert the entrying into an integer value.
      READ(INPUT(lP+1:rP-1),*,IOSTAT=IOS) pIDX
      IF(IOS .NE. 0) THEN
         WRITE(ERR_MSG, 1003) trim(adjustl(INPUT))
         call log_error()
      ELSEIF(pIDX .LT. 0 .OR. pIDX .GT. DIM_M)THEN
         WRITE(ERR_MSG, 1004) trim(adjustl(INPUT))
         call log_error()
      ENDIF

! Convert the entrying into an double precision value.
      READ(INPUT(lQ+1:),*,IOSTAT=IOS) lFracDH(pIDX)
      IF(IOS .NE. 0) THEN
         WRITE(ERR_MSG, 1005)trim(adjustl(INPUT))
         call log_error()
      ENDIF

      RETURN

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFracDH',/           &
         ' Error 1001: Unable to determine phase association for',     &
         ' fracDH. When',/' specifying heat of reaction (DH), the',    &
         ' fraction of DH assigned to',/' each phase must be',         &
         ' given explicitly.',//' Example: fracDH(0) = 0.25  ! 25% of',&
         ' DH is assigned to gas phase',/'          fracDH(1) = 0.75 ',&
         ' ! 75% of DH is assigned to solids phase 1',//' Note:',      &
         ' fracDH(0) + fracDH(1) + ... + frachDH(M) == 1.0',/          &
         /' INPUT: ',A)

 1002 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFracDH',/           &
         ' Error 1002: Input format error for fracDH.',/               &
         /' INPUT: ',A)

 1003 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFracDH',/           &
         ' Error 1003: Unable to determine phase index for fracDH.',// &
         /' INPUT: ',A)

 1004 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFracDH',/           &
         ' Error 1004: Phase index for fracDH exceeds DIM_M!',//       &
         /' INPUT: ',A)

 1005 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFracDH',/           &
         ' Error 1005: Unable to determine fracDH value from input.',/ &
         ' Cannot convert specified value to double precision value.',/&
         /' INPUT: ',A)

 1100 FORMAT(//1X,70('*')/' From: PARSE_RXN',/                         &
         ' Error 1105: fracDH was initially located within the',       &
         ' input line,',/' however its location cannot be determined.',&
         /' INPUT: ',A)



      END SUBROUTINE getFracDH

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: getChemEq(INPUT, lChemEq)                          !
!                                                                      !
!  Purpose: Extract the reaction name from a construct.                !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getChemEq(IN, lChemEq)

      IMPLICIT NONE

! Input line.
      CHARACTER(len=*), INTENT(IN) :: IN
! Name of reaction pulled from input.
      CHARACTER(len=*), INTENT(OUT) :: lChemEq

! read/write output status
      INTEGER IOS

      INTEGER POS, lPOS, rPOS, ldP, lsP, aPOS
      INTEGER lLMAX

! Chemical equations start with the keyword: CHEM_EQ. If this is not a
! continuation of the previous line, search for the keyword and
! flag an error if not found.
      IF(.NOT.MORE_ChemEq) THEN
         lLMAX = LEN_TRIM(IN)
         POS = INDEX(IN,"CHEM_EQ")

         IF(POS == 0) THEN
            WRITE(ERR_MSG, 1105) 'Chem_Eq'
            call log_error()
         ENDIF
! Initialize
         lChemEq = ''
! Update POS to skip over the keyword: CHEM_EQ
         POS = POS+7
      ELSE
         POS = 1
      ENDIF

! Search for quote marks bounding the chemical equation.
      ldP = POS + INDEX(IN(POS:),'"')  ! double quote "
      lsP = POS + INDEX(IN(POS:),"'")  ! single quote '

      IF(ldP .GT. POS .AND. lsP .EQ. POS) THEN
! The chemical equation is bounded by double quotes
         lPOS = ldP
! Search for the second quote mark.
         rPOS = lPOS + INDEX(IN(lPOS+1:),'"')
      ELSEIF(ldP .EQ. POS .AND. lsP .GT. POS) THEN
! The chemical equation is bounded by single quotes
         lPOS = lsP
! Search for the second quote mark.
         rPOS = lPOS + INDEX(IN(lPOS+1:),"'")
      ELSE
! Different errors are thrown depending if this is a continuation
! (MORE_ChemEq) or the start of a chemical equation.
         IF(.NOT.MORE_ChemEq) THEN
            WRITE(ERR_MSG, 1001) trim(adjustl(IN))
         ELSE
            IF(isFracDH(IN)) THEN
              WRITE(ERR_MSG, 1002) 'Keyword fracDH was found inside',        &
                 ' the chemical equation!', trim(adjustl(IN))
            ELSEIF(isDH(IN)) THEN
              WRITE(ERR_MSG, 1002) 'Keyword DH was found inside',            &
                 ' the chemical equation!', trim(adjustl(IN))
            ELSEIF(isArrhenius(IN)) THEN
              WRITE(ERR_MSG, 1002) 'Keyword arrhenius_coeff was found inside',&
                 ' the chemical equation!', trim(adjustl(IN))
            ELSEIF(isReverseCalc(IN)) THEN
              WRITE(ERR_MSG, 1002) 'Keyword reverse_calc was found inside',   &
                 ' the chemical equation!', trim(adjustl(IN))
            ELSEIF(isRxnOrder(IN)) THEN
              WRITE(ERR_MSG, 1002) 'Keyword rxn_order was found inside',      &
                 ' the chemical equation!', trim(adjustl(IN))
            ELSEIF(isThirdBody(IN)) THEN
              WRITE(ERR_MSG, 1002) 'Keyword third_body_param was found inside',            &
                 ' the chemical equation!', trim(adjustl(IN))
            ELSEIF(isPressRxn(IN)) THEN
              WRITE(ERR_MSG, 1002) 'Keyword press_rxn_param was found inside',            &
                 ' the chemical equation!', trim(adjustl(IN))
            ELSEIF(isLTCoeff(IN)) THEN
              WRITE(ERR_MSG, 1002) 'Keyword LT_coeff was found inside',            &
                 ' the chemical equation!', trim(adjustl(IN))
            ELSEIF(isFITCoeff(IN)) THEN
              WRITE(ERR_MSG, 1002) 'Keyword FIT_coeff was found inside',            &
                 ' the chemical equation!', trim(adjustl(IN))
! The entry doesn't match any of the keywords.
            ELSE
              WRITE(ERR_MSG, 1002) 'Unbalanced or missing parentheses', '',  &
                 trim(adjustl(IN))
            ENDIF
         ENDIF
         call log_error()
      ENDIF

! Mismatch/Unbalanced parentheses
      IF(lPOS .EQ. rPOS) THEN
! Different errors are thrown depending if this is a continuation
! (MORE_ChemEq) or the start of a chemical equation.
         IF(.NOT.MORE_ChemEq) THEN
            WRITE(ERR_MSG, 1001) trim(adjustl(IN))
         ELSE
           WRITE(ERR_MSG, 1002) 'Unbalanced or missing parentheses', '',  &
              trim(adjustl(IN))
         ENDIF
         call log_error()
      ENDIF

! Search for an ampersand.
      aPOS = lPOS + INDEX(IN(lPOS+1:),'&')
! An ampersand was found.
      IF(aPOS .GT. lPOS) THEN
         MORE_ChemEq = .TRUE.
! The ampersand should be further to the right than the last quote mark.
         IF(aPOS .LE. rPOS) THEN
            WRITE(ERR_MSG, 1003) trim(adjustl(IN))
            call log_error()
         ENDIF
      ELSE
         MORE_ChemEq = .FALSE.
      ENDIF

! Store the chemical equation.
      WRITE(lChemEq,"(A,1X,A)",IOSTAT=IOS) trim(lChemEq), &
         trim(adjustl(IN(lPOS:rPOS-1)))
      IF(IOS .NE. 0) THEN
         WRITE(ERR_MSG, 1004) trim(lChemEq), trim(adjustl(IN))
         call log_error()
      ENDIF

      RETURN

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getChemEq',/           &
         ' Error 1001: Unbalanced or missing parentheses for chem_eq.',&
         //' INPUT: ',A,//' Example 1: RXN_NAME { chem_eq = ',         &
         '"A + B --> AB" }',//' Example 2: RXN_NAME {',/14X,           &
         'chem_eq = "A + B --> AB"',/14X, 'DH = 2.5d4',/14X,           &
         'fracDH(0) = 1.0',/12X,'}')

 1002 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getChemEq',/           &
         ' Error 1002: Chemical equation continuation input error.',   &
         //'  > ',2A//' INPUT: ',A)

 1003 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getChemEq',/           &
         ' Error 1003: Input format error for chem_eq. An amperand',   &
         ' (&)',/' was located within the parentheses.',//' INPUT: ',A)

 1004 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getChemEq',/           &
         ' Error 1004: Unable to process chemical equation input.',/   &
         ' A possible error is variable overflow as the total length', &
         ' is limited',/' to 1024 characters.',//' lChemEq: ',A,/      &
         ' INPUT: ',A)

 1105 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getChemEq',/           &
         ' Error 1105: chem_eq was initially located within the',      &
         ' input line,',/' however its location cannot be determined.',&
         /' INPUT: ',A)

      END SUBROUTINE getChemEq

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: getArrhenius(INPUT, lArrhen)                       !
!                                                                      !
!  Purpose: Extract the arrhenius parameters from a construct.         !
!                                                                      !
!  Author: Hang Zhou                                Date: May 23, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getArrhenius(INPUT, lArrhen)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! Arrhenius parameters pulled from input.
      DOUBLE PRECISION, INTENT(OUT) :: lArrhen(:)

      INTEGER lQ
      INTEGER lLMAX
! read/write output status
      INTEGER IOS

      lLMAX = LEN_TRIM(INPUT)

      IF(INDEX(INPUT,"ARRHENIUS_COEFF") .EQ. 0) THEN
         WRITE (ERR_MSG, 1100) trim(adjustl(INPUT))
         call log_error()
      ENDIF

      lQ = INDEX(INPUT(:lLMAX),'=')

      IF(lQ .EQ. 0) THEN
         WRITE (ERR_MSG, 1001) trim(adjustl(INPUT))
         call log_error()
      ENDIF

! Convert the entrying into an double precision value.
      READ(INPUT(lQ+1:),*,IOSTAT=IOS) lArrhen(1:3)
      IF(IOS .NE. 0) THEN
         WRITE(ERR_MSG, 1002) trim(adjustl(INPUT))
         call log_error()
      ENDIF

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getArrhenius',/      &
         ' Error 1001: Input format error for arrhenius_coeff.',/    &
     /' INPUT: ',A)

 1002 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getArrhenius',/      &
         ' Error 1002: Unable to determine arrhenius_coeff value ',  &
         ' from input.',/                                            &
      ' Cannot convert specified value to double precision value.',/ &
         /' INPUT: ',A)

 1100 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getArrhenius',/      &
         ' Error 1100: arrhenius_coeff was initially located ',/     &
         ' within the input line, however its location cannot ',/    &
         ' be determined.', /                                        &
         ' INPUT: ',A)

      RETURN
      END SUBROUTINE getArrhenius

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: getReverseCalc(INPUT, lReverse)                    !
!                                                                      !
!  Purpose: Extract the infos about calculating reverse reaction rate  !
!           from a construct.                                          !
!                                                                      !
!  Author: Hang Zhou                                Date: May 23, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getReverseCalc(INPUT, lReverse)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! Infos about calculating reverse reaction rate  pulled from input.
      CHARACTER(len=*), INTENT(OUT) :: lReverse

      INTEGER lQ
      INTEGER lLMAX

      lLMAX = LEN_TRIM(INPUT)

      IF(INDEX(INPUT,"REVERSE_CALC") .EQ. 0) THEN
         WRITE (ERR_MSG, 1100) trim(adjustl(INPUT))
         call log_error()
      ENDIF

      lQ = INDEX(INPUT(:lLMAX),'=')

      IF(lQ .EQ. 0) THEN
         WRITE (ERR_MSG, 1001) trim(adjustl(INPUT))
         call log_error()
      ENDIF

      lReverse = trim(INPUT(lQ+1:))

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getReverseCalc',/      &
         ' Error 1001: Input format error for reverse_calc.',/         &
     /' INPUT: ',A)

 1100 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getReverseCalc',/    &
         ' Error 1100: reverse_calc was initially located within ',/ &
         'the input line, however its location cannot be determined.',/ &
         /' INPUT: ',A)
      RETURN
      END SUBROUTINE getReverseCalc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: getRxnOrder(INPUT, lrxnOrder)                      !
!                                                                      !
!  Purpose: Extract the reaction orders from a construct.              !
!                                                                      !
!  Author: Hang Zhou                                Date: May 23, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getRxnOrder(INPUT, lrxnOrder)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! Reaction orders pulled from input.
      CHARACTER(len=*), INTENT(OUT) :: lrxnOrder

      INTEGER lQ
      INTEGER lLMAX
      ! read/write output status
      INTEGER IOS
      INTEGER POS, aPOS

      IF(.NOT. MORE_RxnOrder) THEN
         lLMAX = LEN_TRIM(INPUT)

         IF(INDEX(INPUT,"RXN_ORDER") .EQ. 0) THEN
            WRITE (ERR_MSG, 1100) trim(adjustl(INPUT))
            call log_error()
         ENDIF

         lQ = INDEX(INPUT(:lLMAX),'=')
         IF(lQ .EQ. 0) THEN
            WRITE (ERR_MSG, 1001) trim(adjustl(INPUT))
            call log_error()
         ENDIF

         POS = lQ + 1
      ELSE
         POS = 1
      ENDIF

      ! Search for an ampersand.
      aPOS = INDEX(INPUT,"&")
      ! An ampersand was found.
      IF(aPOS .GT. POS) THEN
         MORE_RxnOrder = .TRUE.
      ! The ampersand should be at the end of the line.
         IF(LEN(trim(adjustl(INPUT(aPOS+1:)))) .NE. 0) THEN
            WRITE(ERR_MSG, 1002) trim(adjustl(INPUT))
            call log_error()
         ENDIF
         WRITE(lrxnOrder,"(A,1X,A)",IOSTAT=IOS) trim(lrxnOrder), &
         trim(adjustl(INPUT(POS:aPOS-1)))
      ELSE
         MORE_RxnOrder = .FALSE.
         WRITE(lrxnOrder,"(A,1X,A)",IOSTAT=IOS) trim(lrxnOrder), &
         trim(adjustl(INPUT(POS:)))
      ENDIF

      IF(IOS .NE. 0) THEN
         WRITE(ERR_MSG, 1003) trim(lrxnOrder), trim(adjustl(INPUT))
         call log_error()
      ENDIF

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getRxnOrder',/         &
         ' Error 1001: Input format error for rxn_order.',/            &
        /' INPUT: ',A)

 1002 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getRxnOrder',/           &
         ' Error 1002: Input format error for rxn_order. An amperand',   &
         ' (&)',/' was located within the parameters.',//' INPUT: ',A)

 1003 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getRxnOrder',/         &
         ' Error 1003: Unable to process reaction order input.',/      &
         ' A possible error is variable overflow as the total length', &
         ' is limited',/' to 1024 characters.',//' lRxnOrder: ',A,/    &
         ' INPUT: ',A)

 1100 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getRxnOrder',/       &
         ' Error 1100: rxn_order was initially located within ',/    &
         'the input line, however its location cannot be determined.',/&
         /' INPUT: ',A)

      RETURN
      END SUBROUTINE getRxnOrder

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: getThirdBody(INPUT, lMmodel, lMcoeff, lMcoeff_str) !
!                                                                      !
!  Purpose: Extract the infos about third body from a construct.       !
!                                                                      !
!  Author: Hang Zhou                                Date: May 23, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getThirdBody(INPUT, lMmodel, lMcoeff, lMcoeff_str)

      use rxns, only: SPECIES_ALIAS_g

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! Infos about third body pulled from input.
      CHARACTER(len=*), INTENT(OUT) :: lMmodel
! Third body coefficients. It can be different from 1
      DOUBLE PRECISION, INTENT(OUT) :: lMcoeff(:)
      CHARACTER(len=*), INTENT(INOUT) :: lMcoeff_str
      CHARACTER(len=512) :: char_tmp
! array of species names and coefficient given as third-body coefficient
      CHARACTER(len=32),  DIMENSION(:), ALLOCATABLE :: species_M
      DOUBLE PRECISION,   DIMENSION(:), ALLOCATABLE :: coeff_M
      CHARACTER(LEN=32) species_tmp, species_tmp2

      INTEGER lQ
      INTEGER lLMAX
! read/write output status
      INTEGER IOS
      INTEGER LL, CC
! index of space
      INTEGER index_s
! number of coefficients given in the third body block
      INTEGER num_coeff
! index of colons in the characters for coefficient
      INTEGER index_c(DIM_N_G)
      INTEGER POS, rPOS, ePOS, aPOS

      ! Search for an ampersand.
      aPOS = INDEX(INPUT,"&")
      IF(.NOT. MORE_thirdBody) THEN
         lLMAX = LEN_TRIM(INPUT)
         IF(INDEX(INPUT,"THIRD_BODY_PARAM") .EQ. 0) THEN
            WRITE (ERR_MSG, 1100) trim(adjustl(INPUT))
            call log_error()
         ENDIF

         lQ = INDEX(INPUT(:lLMAX),'=')
         IF(lQ .EQ. 0) THEN
            WRITE (ERR_MSG, 1001) trim(adjustl(INPUT))
            call log_error()
         ENDIF
         POS = lQ + 1
         rPOS = INDEX(INPUT,'"')

         lMcoeff = 1.0d0
         IF(rPOS .EQ. 0 .AND. aPOS .EQ. 0) THEN
            READ(INPUT(POS:),*,IOSTAT=IOS) lMmodel
         ELSEIF(rPOS .EQ. 0 .AND. aPOS .NE. 0) THEN
            READ(INPUT(POS:aPOS-1),*,IOSTAT=IOS) lMmodel
         ELSE
            ! Search for the second quote mark.
            ePOS = rPOS + INDEX(INPUT(rPOS+1:),'"')
            READ(INPUT(POS:ePOS),*,IOSTAT=IOS) lMmodel, lMcoeff_str
         ENDIF
      ELSE
         POS = 1
         rPOS = INDEX(INPUT,'"')
         ! Search for the second quote mark.
         ePOS = rPOS + INDEX(INPUT(rPOS+1:),'"')
         WRITE(lMcoeff_str,"(A,1X,A)",IOSTAT=IOS) trim(lMcoeff_str), &
            trim(adjustl(INPUT(rPOS+POS:ePOS-1)))
         IF(IOS .NE. 0) THEN
            WRITE(ERR_MSG, 1009) lMmodel, trim(adjustl(lMcoeff_str))
            call log_error()
         ENDIF
      ENDIF

    IF(trim(lMmodel) .NE. 'M_ALL') THEN
       species_tmp = trim(lMmodel(3:))
         ! Remove case sensitivity.
         CALL MAKE_UPPER_CASE (species_tmp,32)
       IF(.NOT. ANY(SPECIES_ALIAS_g == trim(species_tmp))) THEN
          WRITE(ERR_MSG, 1004) trim(adjustl(INPUT))
          CALL LOG_ERROR()
       ENDIF
    ENDIF

      ! An ampersand was found.
      IF(aPOS .GT. POS) THEN
         MORE_thirdBody = .TRUE.
         ! The ampersand should be further to the end of the line.
         IF(LEN(trim(adjustl(INPUT(aPOS+1:)))) .NE. 0) THEN
            WRITE(ERR_MSG, 1010) trim(adjustl(INPUT))
            call log_error()
         ENDIF
      ELSE
         MORE_thirdBody = .FALSE.
         IF(LEN_TRIM(lMcoeff_str) .NE. 0) THEN
          ! split the coefficients of third body species
          ! get the number of coefficients given as third-body species,
          ! and the index of colons for each third-body species
            num_coeff = 0
            DO LL=1, LEN_TRIM(lMcoeff_str)
               IF(lMcoeff_str(LL:LL).EQ. ':') THEN
                  num_coeff = num_coeff+1
                  index_c(num_coeff) = LL
               ENDIF
            ENDDO

            ALLOCATE(species_M(num_coeff))
            ALLOCATE(coeff_M(num_coeff))

	    DO CC = 1, num_coeff
	       IF(CC .EQ. 1) species_M(CC) = trim(adjustl(lMcoeff_str(1:index_c(cc)-1)))
	       IF(CC .LT. num_coeff) THEN
	          char_tmp = adjustl(lMcoeff_str(index_c(cc)+1:index_c(cc+1)-1))
                  ! get the index of the first space in char_tmp
                  ! there must be space between the coefficient of the previous species and
                  ! the name of the next species

	          index_s = INDEX(trim(char_tmp), ' ')
	          IF(index_s .EQ. 0) THEN
	       	     WRITE(ERR_MSG, 1005) trim(adjustl(INPUT))
                     call log_error()
                  ENDIF
                  READ(char_tmp(1:index_s-1), *) coeff_M(CC)
                  species_M(CC+1) = trim(adjustl(char_tmp(index_s+1:)))
               ELSE
                  char_tmp = lMcoeff_str(index_c(cc)+1:)
                  READ(char_tmp, *) coeff_M(CC)
               ENDIF
               DO LL = 1, DIM_N_G
                  species_tmp = trim(species_M(CC))
                  ! Remove case sensitivity.
                  CALL MAKE_UPPER_CASE (species_tmp,32)
                  IF(SPECIES_ALIAS_g(LL) .EQ. trim(species_tmp)) THEN
                     IF(lMcoeff(LL) .NE. 1.0d0) THEN
                        WRITE(ERR_MSG, 1007) SPECIES_ALIAS_g(LL), trim(adjustl(INPUT))
                        call log_error()
                     ENDIF
                    lMcoeff(LL) = coeff_M(CC)
                    EXIT
                  ENDIF
                  IF(LL .EQ. DIM_N_G) THEN
                     WRITE(ERR_MSG, 1006) trim(adjustl(INPUT))
                     call log_error()
                  ENDIF
               ENDDO
            ENDDO

            IF((trim(lMmodel) .NE. 'M_ALL')) THEN
               IF(num_coeff .EQ. one) THEN
                  species_tmp = trim(species_M(1))
                  species_tmp2 =trim(lMmodel(3:))
                  ! Remove case sensitivity.
                  CALL MAKE_UPPER_CASE (species_tmp,32)
                  CALL MAKE_UPPER_CASE (species_tmp2,32)
                  IF(trim(species_tmp) .NE. trim(species_tmp2)) THEN
                     WRITE(ERR_MSG, 1003) trim(adjustl(INPUT))
                     CALL LOG_WARNING()
                  ENDIF
               ELSE
                  WRITE(ERR_MSG, 1003) trim(adjustl(INPUT))
                  CALL LOG_WARNING()
               ENDIF
            ENDIF
            DEALLOCATE(species_M)
            DEALLOCATE(coeff_M)

         ENDIF
      ENDIF

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getThirdBody',/        &
         ' Error 1001: Input format error for third_body_param.',/     &
        /' INPUT: ',A)

 1003 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getThirdBody',/        &
         ' Warning 1003: A specific species is used as third body.',/  &
         ' Only coefficient for this species given in the third body',/&
         ' coefficients list will be used! If it is not given in the',/&
         ' list, one will be used!',/                                  &
         /' INPUT: ',A)

 1004 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getThirdBody',/        &
         ' Error 1004: The specific species given in third body model ',&
         ' is not in species list. The third body format needs to be ',/&
     ' "M_all" for all species used as third body or ',/           &
     ' "M_speciesnName" when one specific species is used as third'/&
     ' body, and the species must be given in the species list',/  &
         /' INPUT: ',A)

 1005 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getThirdBody',/       &
         ' Error 1005: Space must be added between coefficient of ', &
         ' species and next species name.',/                         &
         /' INPUT: ',A)

 1006 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getThirdBody',/        &
         ' Error 1006: The species given in the third body coefficients',&
         ' list is not in species list.',/                             &
         /' INPUT: ',A)

 1007 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getThirdBody',/        &
         ' Error 1007: The coefficient for species ', A, 'has already'/&
         ' been modified. This may happen from duplicate coefficients',&
         ' given for this species. ',/                           &
         /' INPUT: ',A)

 1009 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getThirdBody',/       &
         ' Error 1009: Unable to process third_body_param input.',/   &
         ' A possible error is variable overflow as the total length',&
         ' is limited',/' to 1024 characters.',//' model: ',A,//       &
         ' parameters: ',A)


 1010 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getThirdBody',/           &
         ' Error 1010: Input format error for third_body_param. An amperand',   &
         ' (&)',/' was located within the parentheses.',//' INPUT: ',A)

 1100 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getThirdBody',/      &
         ' Error 1100: third_body_param was initially located ',/    &
         ' within the input line however its location cannot ',/     &
         ' be determined.',/                                         &
         ' INPUT: ',A)

      RETURN
      END SUBROUTINE getThirdBody


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: getPressRxn(INPUT, lpreModel, lpreCoefflList,      !
!                               lArrhenPreList, lpre_param_str)        !
!                                                                      !
!  Purpose: Extract the infos about parameters for pressure-dependent  !
!           reacitons from a construct.                                !
!                                                                      !
!  Author: Hang Zhou                                Date: May 23, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getPressRxn(INPUT, lpreModel, lpreCoeffl, lArrhenPre, lid_rxn, lpre_param_str)

      USE resize, only: real_grow, real_grow2, real_grow2_reverse, real_grow3_reverse2

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! Infos about model used for pressure-dependent reactions pulled from input.
      CHARACTER(len=*), INTENT(OUT) :: lpreModel
! index of the reaction
      INTEGER, INTENT(IN) :: lid_rxn
! Coefficients in the model
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: lpreCoeffl
! Arrhenius coefficients for low or high pressure limits
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: lArrhenPre
      CHARACTER(len=*), INTENT(INOUT) :: lpre_param_str

      INTEGER lQ
      INTEGER lLMAX
! read/write output status
      INTEGER IOS
      INTEGER POS, rPOS, aPOS, ePOS
      INTEGER index_arrhen, index_coeff, num_coeff, ii
      INTEGER end_press, end_arrhen, num_comma, old_size, index_tmp
      INTEGER, DIMENSION(:), ALLOCATABLE :: index_comma

      ! Search for an ampersand.
      aPOS = INDEX(INPUT,"&")

      IF(.NOT. MORE_PressParam) THEN
         lLMAX = LEN_TRIM(INPUT)
         IF(INDEX(INPUT,"PRESS_RXN_PARAM") .EQ. 0) THEN
            WRITE (ERR_MSG, 1100) trim(adjustl(INPUT))
            call log_error()
         ENDIF

         lQ = INDEX(INPUT(:lLMAX),'=')
         IF(lQ .EQ. 0) THEN
            WRITE (ERR_MSG, 1001) trim(adjustl(INPUT))
            call log_error()
         ENDIF
         POS = lQ + 1
         rPOS = INDEX(INPUT,'"')

         IF(rPOS .EQ. 0 .AND. aPOS .EQ. 0) THEN
            WRITE(ERR_MSG, 1003) trim(adjustl(INPUT))
            call log_error()
         ELSEIF(rPOS .EQ. 0 .AND. aPOS .NE. 0) THEN
            READ(INPUT(POS:aPOS-1),*,IOSTAT=IOS) lpreModel
         ELSE
            ! Search for the second quote mark.
            ePOS = rPOS + INDEX(INPUT(rPOS+1:),'"')
            READ(INPUT(POS:ePOS),*,IOSTAT=IOS) lpreModel, lpre_param_str
         ENDIF
      ELSE
         POS = 1
         rPOS = INDEX(INPUT,'"')
         ! Search for the second quote mark.
         ePOS = rPOS + INDEX(INPUT(rPOS+1:),'"')
         WRITE(lpre_param_str,"(A,1X,A)",IOSTAT=IOS) trim(lpre_param_str), &
            trim(adjustl(INPUT(rPOS+POS:ePOS-1)))
         IF(IOS .NE. 0) THEN
            WRITE(ERR_MSG, 1009) lPreModel, trim(adjustl(lpre_param_str))
            call log_error()
         ENDIF
      ENDIF

      ! An ampersand was found.
      IF(aPOS .GT. POS) THEN
         MORE_PressParam = .TRUE.
      ! The ampersand should be further to the end of the line.
         IF(LEN(trim(adjustl(INPUT(aPOS+1:)))) .NE. 0) THEN
            WRITE(ERR_MSG, 1010) trim(adjustl(INPUT))
            call log_error()
         ENDIF
      ELSE
         MORE_PressParam = .FALSE.
         IF(INDEX(lpre_param_str,'ARRHENIUS_PRESS') .EQ. 0) THEN
            WRITE (ERR_MSG, 1003) lPreModel, trim(adjustl(lpre_param_str))
            call log_error()
         ENDIF
         IF((INDEX(lpreModel, 'LINDEMANN') .EQ. 0) .AND. (INDEX(lpre_param_str, 'PRESS_COEFF') .EQ. 0)) THEN
            WRITE(ERR_MSG, 1004) lpreModel, trim(adjustl(lpre_param_str))
            call log_error()
         ELSEIF((INDEX(lpreModel, 'LINDEMANN') .NE. 0) .AND. (INDEX(lpre_param_str, 'PRESS_COEFF') .NE. 0)) THEN
            WRITE(ERR_MSG, 1005) lpreModel, trim(adjustl(lpre_param_str))
          call log_warning()
         ENDIF

         index_arrhen = INDEX(lpre_param_str, 'ARRHENIUS_PRESS')
         index_coeff = INDEX(lpre_param_str, 'PRESS_COEFF')

         IF(index_coeff .EQ. 0) THEN ! for LINDEMANN model
            READ(lpre_param_str(index_arrhen+16:), *) lArrhenPre(lid_rxn,1,1:3)
         ELSE
            num_coeff = 0
          IF(index_arrhen .LT. index_coeff) THEN
             ! Get the number of coefficient given in press_coeff
             DO ii=index_coeff+12, LEN_TRIM(lpre_param_str)
                IF((lpre_param_str(ii-1:ii-1) .EQ. ' ') .AND. lpre_param_str(ii:ii) .NE. ' ') &
                   num_coeff = num_coeff+1
             ENDDO
             end_press = LEN(lpre_param_str)
             end_arrhen = index_coeff-1
          ELSE
             DO ii=index_coeff+12, index_arrhen-1
                IF((lpre_param_str(ii-1:ii-1) .EQ. ' ') .AND. lpre_param_str(ii:ii) .NE. ' ') &
                   num_coeff = num_coeff+1
             ENDDO
             end_press = index_arrhen-1
             end_arrhen = LEN(lpre_param_str)
          ENDIF

          IF((INDEX(lpreModel, 'TROE') .NE. 0) .AND. (num_coeff .GT. 4)) THEN
             WRITE(ERR_MSG, 1006) 'TROE', trim(adjustl(lpre_param_str))
               call log_error()
          ENDIF

          IF((INDEX(lpreModel, 'SRI') .NE. 0) .AND. (num_coeff .GT. 5)) THEN
             WRITE(ERR_MSG, 1006) 'SRI', trim(adjustl(lpre_param_str))
              call log_error()
          ENDIF

          IF(INDEX(lpreModel, 'PLOG') .EQ. 0) THEN ! for LINDEMANN, Troe or SRI model
             READ(lpre_param_str(index_coeff+12:end_press), *) lpreCoeffl(lid_rxn, 1:num_coeff)
             READ(lpre_param_str(index_arrhen+16:end_arrhen), *) lArrhenPre(lid_rxn, 1,1:3)
             index_tmp = INDEX(lpre_param_str,',')
             IF((index_tmp .NE. 0) .and. (trim(lpre_param_str(index_tmp+1:end_arrhen)).NE. '')) THEN
                WRITE(ERR_MSG, 1008) lpreModel, trim(adjustl(lpre_param_str))
                  call log_error()
             ENDIF
          ELSE ! for PLOG model
             old_size = size(lpreCoeffl, 2)
             IF(num_coeff .GT. old_size) THEN
                CALL real_grow2(lpreCoeffl,num_coeff)
                lpreCoeffl(:, old_size+1:) = UNDEFINED
             ENDIF
             old_size = size(lArrhenPre, 2)
             IF(num_coeff .GT. old_size) THEN
                CALL real_grow3_reverse2(lArrhenPre,num_coeff)
                lArrhenPre(:, old_size+1:,:) = UNDEFINED
             ENDIF
             READ(lpre_param_str(index_coeff+12:end_press), *) lpreCoeffl(lid_rxn,1:num_coeff)
             DO ii=2, num_coeff
                 IF(lpreCoeffl(lid_rxn,ii) .LT. lpreCoeffl(lid_rxn,ii-1)) THEN
                   WRITE(ERR_MSG, 1011) trim(adjustl(lpre_param_str))
                   call log_error()
                 ENDIF
             ENDDO
             ALLOCATE(index_comma(num_coeff-1))
             ! the Arrhenius coefficients for each pressure is separated by comma.
             ! get the number of commas
             num_comma = 0
             DO ii=index_arrhen+16, end_arrhen
                IF(lpre_param_str(ii-1:ii-1) .EQ. ',') THEN
                   num_comma = num_comma+1
                 IF(num_comma .GT. num_coeff-1) THEN
                      WRITE(ERR_MSG, 1007) trim(adjustl(lpre_param_str))
                      call log_error()
                 ENDIF
                 index_comma(num_comma) = ii
                ENDIF
             ENDDO

             IF(num_comma .LT. num_coeff-1) THEN
                WRITE(ERR_MSG, 1007) trim(adjustl(lpre_param_str))
                  call log_error()
             ENDIF
             IF(num_comma .EQ. 1) THEN
                READ(lpre_param_str(index_arrhen+16:index_comma(1)), *) lArrhenPre(lid_rxn,1, 1:3)
                READ(lpre_param_str(index_comma(1):end_arrhen), *) lArrhenPre(lid_rxn,2, 1:3)
             ELSE
                DO ii = 1, num_comma
                   IF(ii .EQ. 1) THEN
                      READ(lpre_param_str(index_arrhen+16:index_comma(ii)), *) lArrhenPre(lid_rxn,ii, 1:3)
                   ELSE IF(ii .LT. num_comma) THEN
                      READ(lpre_param_str(index_comma(ii-1):index_comma(ii)), *) lArrhenPre(lid_rxn,ii, 1:3)
                   ELSE
                      READ(lpre_param_str(index_comma(ii-1):index_comma(ii)), *) lArrhenPre(lid_rxn,ii, 1:3)
                      READ(lpre_param_str(index_comma(ii):end_arrhen), *) lArrhenPre(lid_rxn,ii+1, 1:3)
                   ENDIF
                ENDDO
             ENDIF
          ENDIF
         ENDIF
      ENDIF

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getPressRxn',/           &
         ' Error 1001: Input format error for press_rxn_param.',/        &
      /' INPUT: ',A)

 1003 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getPressRxn',/               &
         ' Error 1003: Unable to determine arrhenius_press in press_rxn_param value from input.',/     &
         ' An extra set of Arrhenius coefficients must be given for pressure-dependent reactions.',/&
         /' INPUT: ',A)

 1004 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getPressRxn',/               &
         ' Error 1004: press_coeff must be given in press_rxn_param of the input for model Troe, SRI or PLOG.',/ &
         /' INPUT: ',A)

 1005 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getPressRxn',/               &
         ' Warning 1005: Lindemann model is used for the calculation of the reaction rates of pressure-dependent reactions.',/ &
         ' The parameters given by press_coeff in press_rxn_param will not be used.',/&
         /' INPUT: ',A)

 1006 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getPressRxn',/                       &
         ' Error 1006: Too many coefficients are given in press_coeff.', / &
         ' 3 or 4 coefficients for Troe model and 3 or 5 coefficients for SRI model are required.' / &
         /' INPUT: ',A)

 1007 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getPressRxn',/                             &
         ' Error 1007: Number of pressures given in press_coeff is different from the ', / &
         ' number of the groups of arrhenius coefficients give in arrhenius_press for ', / &
         ' PLOG model. In PLOG model, each pressure is given by press_coeff. And the',/    &
         ' corresponding arrhenius coefficients are given in arrhenius_press separated ',/&
         ' by comma. INPUT: ',A)

 1008 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getPressRxn',/                             &
         ' Error 1008: More then one group of Arrhenius coefficients are given in ', /     &
         ' arrhenius_press for model ', A, '.' / &
         ' INPUT: ',A)

 1009 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getPressRxn',/           &
         ' Error 1009: Unable to process press_rxn_param input.',/   &
         ' A possible error is variable overflow as the total length', &
         ' is limited',/' to 1024 characters.',//' model: ',A,//      &
         ' parameters: ',A)

 1010 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getPressRxn',/           &
         ' Error 1010: Input format error for press_rxn_param. An amperand',   &
         ' (&)',/' was located within the parentheses.',//' INPUT: ',A)


 1011 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getPressRxn',/                             &
         ' Error 1011: Coefficients given in press_coeff for PLOG model is not in ascending', / &
         ' order of pressure and need to be fixed.', / &
         ' INPUT: ',A)

 1100 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getPressRxn',/            &
         ' Error 1100: press_rxn_param was initially located within the ',&
         ' input line, however its location cannot be determined.',/        &
         /' INPUT: ',A)

      RETURN

      END SUBROUTINE getPressRxn

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: resizePressCoeff(uArrhenPre, upreCoeffl, new_size) !
!                                                                      !
!  Purpose: resize arrhenius_coeff_press and press_rxn_coeff           !
!           based on parameters given in PLOG model.                   !
!                                                                      !
!  Author: Hang Zhou                               Date: June 06, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!     SUBROUTINE resizePressCoeff(uArrhenPre, upreCoeffl, new_size)

         ! Coefficients in the model.
!         DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: upreCoeffl
!        ! Arrhenius coefficients for low or high pressure limits
!         DOUBLE PRECISION, DIMENSION(:,:,:),INTENT(INOUT) :: uArrhenPre
!! new size: number of pressures
    ! INTEGER, INTENT(IN) :: new_size

!     INTEGER size_press, size_arrhen, num_rxn
!     DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: press_coeff_tmp
!     DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: arrhen_press_tmp

!     num_rxn = size(upreCoeffl, 1)
!     size_press = size(upreCoeffl, 2)
!     size_arrhen = size(uArrhenPre, 2)

!     IF(new_size .GT. size_press) THEN
!        allocate(press_coeff_tmp(num_rxn, new_size))
!        press_coeff_tmp(:, 1:size_press) = upreCoeffl(:, 1:size_press)
!        call move_alloc(press_coeff_tmp,upreCoeffl)
!        upreCoeffl(:, size_press:) = UNDEFINED
!     ENDIF

!     IF(new_size .GT. size_arrhen) THEN
!        allocate(arrhen_press_tmp(num_rxn, new_size, 3))
!        arrhen_press_tmp(:, 1:size_arrhen, :) = uArrhenPre(:, 1:size_arrhen, :)
!        call move_alloc(arrhen_press_tmp,uArrhenPre)
!        upreCoeffl(:, size_arrhen:,:) = UNDEFINED
!     ENDIF

 !     END SUBROUTINE resizePressCoeff


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: getLTCoeff(INPUT, lLTCoeff)                        !
!                                                                      !
!  Purpose: Extract the coefficients for Landau_Teller reactions       !
!           from a construct.                                          !
!                                                                      !
!  Author: Hang Zhou                                Date: May 23, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getLTCoeff(INPUT, lLTCoeff)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! Third body coefficients. It can be different from 1
      DOUBLE PRECISION, INTENT(OUT) :: lLTCoeff(:)

      INTEGER lQ
      INTEGER lLMAX
! read/write output status
      INTEGER IOS

      lLMAX = LEN_TRIM(INPUT)

      IF(INDEX(INPUT,'LT_COEFF') .EQ. 0) THEN
         WRITE (ERR_MSG, 1100) trim(adjustl(INPUT))
         call log_error()
      ENDIF

      lQ = INDEX(INPUT(:lLMAX),'=')

      IF(lQ .EQ. 0) THEN
         WRITE (ERR_MSG, 1001) trim(adjustl(INPUT))
         call log_error()
      ENDIF

! Convert the entrying into an double precision value.
      READ(INPUT(lQ+1:),*,IOSTAT=IOS) lLTCoeff(1:2)
      IF(IOS .NE. 0) THEN
         WRITE(ERR_MSG, 1002) trim(adjustl(INPUT))
         call log_error()
      ENDIF

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getLTCoeff',/          &
         ' Error 1001: Input format error for LT_coeff.',/             &
         /' INPUT: ',A)

 1002 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getLTCoeff',/          &
         ' Error 1002: Unable to determine LT_coeff value from input.',/&
         ' Cannot convert specified value to double precision value.',/&
         /' INPUT: ',A)

 1100 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getLTCoeff',/                         &
         ' Error 1100: LT_coeff was initially located within the input line',&
         /' however its location cannot be determined.',/              &
         /' INPUT: ',A)

      RETURN
      END SUBROUTINE getLTCoeff


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: getFITCoeff(INPUT, lFITModel, lFITCoeff,           !
!                               lfit_param_str)                        !
!                                                                      !
!  Purpose: Extract the infos for the optional rate-constant fit       !
!           expression from a construct.                               !
!                                                                      !
!  Author: Hang Zhou                                Date: May 23, 2024 !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getFITCoeff(INPUT, lFITModel, lFITCoeff, lfit_param_str)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! Infos about third body pulled from input.
      CHARACTER(len=*), INTENT(OUT) :: lFITModel
! Third body coefficients. It can be different from 1
      DOUBLE PRECISION, INTENT(OUT) :: lFITCoeff(:)
      CHARACTER(len=*), INTENT(INOUT) :: lfit_param_str

      INTEGER lQ
      INTEGER lLMAX
! read/write output status
      INTEGER IOS
      INTEGER POS, rPOS, ePOS, aPOS
      INTEGER index_coeff, num_coeff, ii

      ! Search for an ampersand.
      aPOS = INDEX(INPUT,"&")

      IF(.NOT. MORE_fitParam) THEN
         lLMAX = LEN_TRIM(INPUT)

         IF(INDEX(INPUT,'FIT_COEFF') .EQ. 0) THEN
            WRITE (ERR_MSG, 1100) trim(adjustl(INPUT))
            call log_error()
         ENDIF

         lQ = INDEX(INPUT(:lLMAX),'=')
         IF(lQ .EQ. 0) THEN
            WRITE (ERR_MSG, 1001) trim(adjustl(INPUT))
            call log_error()
         ENDIF
         POS = lQ + 1
         rPOS = INDEX(INPUT,'"')

         IF(rPOS .EQ. 0 .AND. aPOS .EQ. 0) THEN
            WRITE(ERR_MSG, 1002) trim(adjustl(INPUT))
            call log_error()
         ELSEIF(rPOS .EQ. 0 .AND. aPOS .NE. 0) THEN
            READ(INPUT(POS:aPOS-1),*,IOSTAT=IOS) lFITModel
         ELSE
            ! Search for the second quote mark.
            ePOS = rPOS + INDEX(INPUT(rPOS+1:),'"')
            READ(INPUT(POS:ePOS),*,IOSTAT=IOS) lFITModel, lfit_param_str
         ENDIF
      ELSE
         POS = 1
         rPOS = INDEX(INPUT,'"')
         ! Search for the second quote mark.
         ePOS = rPOS + INDEX(INPUT(rPOS+1:),'"')
         WRITE(lfit_param_str,"(A,1X,A)",IOSTAT=IOS) trim(lfit_param_str), &
            trim(adjustl(INPUT(rPOS+POS:ePOS-1)))
         IF(IOS .NE. 0) THEN
            WRITE(ERR_MSG, 1009) lFITModel, trim(adjustl(lfit_param_str))
            call log_error()
         ENDIF
      ENDIF

      ! An ampersand was found.
      IF(aPOS .GT. POS) THEN
         MORE_fitParam = .TRUE.
         ! The ampersand should be further to the right than the last quote mark.
         IF(LEN(trim(adjustl(INPUT(aPOS+1:)))) .NE. 0) THEN
            WRITE(ERR_MSG, 1003) trim(adjustl(INPUT))
            call log_error()
         ENDIF
      ELSE
         MORE_fitParam = .FALSE.
         index_coeff = INDEX(lfit_param_str, ':')
         num_coeff = 0
         DO ii=index_coeff+1, LEN_TRIM(lfit_param_str)
            IF((lfit_param_str(ii-1:ii-1) .EQ. ' ') .AND. lfit_param_str(ii:ii) .NE. ' ') &
             num_coeff = num_coeff+1
         ENDDO
         READ(lfit_param_str(index_coeff+1:), *) lFITCoeff(1:num_coeff)
      ENDIF

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFITCoeff',/         &
         ' Error 1001: Input format error for FIT_coeff.',/            &
        /' INPUT: ',A)

 1002 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFITCoeff',/          &
         ' Error 1002: Unable to determine the coefficients in ',/      &
         ' FIT_coeff value from input.',/                               &
         /' INPUT: ',A)

 1003 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFITCoeff',/           &
         ' Error 1003: Input format error for FIT_Coeff. An amperand',   &
         ' (&)',/' was located within the parentheses.',//' INPUT: ',A,// &
         ' to the data file.')

 1009 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFITCoeff',/           &
         ' Error 1009: Unable to process fit_coeff input.',/   &
         ' A possible error is variable overflow as the total length', &
         ' is limited',/' to 1024 characters.',//' model: ',A,//      &
         ' parameters: ',A)

 1100 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFITCoeff',/                         &
         ' Error 1105: FIT_coeff was initially located within the input line',&
         /' however its location cannot be determined.'/,              &
         /' INPUT: ',A)

      RETURN
      END SUBROUTINE getFITCoeff

   END SUBROUTINE PARSE_RXN

END MODULE PARSE_RXN_MOD
