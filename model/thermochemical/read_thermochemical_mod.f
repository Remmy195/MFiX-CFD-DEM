!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: read_thermochemical_mod.f                              C
!  Purpose: Read thermochemical data                                   C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: None                                C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

#include "error.inc"

MODULE read_thermochemical

   use error_manager
   USE make_upper_case_mod, only: make_upper_case
   use param1, only: zero, undefined
   use resize, only: real_grow3, real_grow4

  IMPLICIT NONE

! Atomic weight data structure to compute molecular weight
  TYPE ATOMIC
     INTEGER :: ID
     CHARACTER(LEN=2) :: SYMBOL
     CHARACTER(LEN=16) :: NAME
     CHARACTER(LEN=16) :: WEIGHT_STRING
     DOUBLE PRECISION :: WEIGHT
     LOGICAL :: RADIOACTIVE
  END TYPE ATOMIC

  TYPE(ATOMIC), DIMENSION(118) :: ELEMENT
  DATA ELEMENT /&
#include "atomic.inc"
      /

  INTEGER :: ATOMIC_MAX = 118

CONTAINS

!      Program Test; CALL Read_Therm_tester; END Program Test

      SUBROUTINE READ_Therm_tester(THERM)
      Implicit none
      CHARACTER(LEN=*), INTENT(IN) :: THERM
      DOUBLE PRECISION Ahigh(7), Alow(7)
      DOUBLE PRECISION Thigh, Tlow, Tcom, MW
      DOUBLE PRECISION Cp1, h1, h2, Hf298oR
      CHARACTER(LEN=132) :: PATH
      CHARACTER(LEN=18) :: SPECIES
      integer funit, IER
      CHARACTER(len=255) FILENAME
      LOGICAL LocalCopy

      SPECIES = 'CH4'
      PATH = '.'
      funit = 5

      INQUIRE(FILE=TRIM(THERM),EXIST=LocalCopy)
      IF(LocalCopy)Then
        OPEN(UNIT=funit,FILE=TRIM(THERM))
      ELSE
        FILENAME = './BURCAT.THR'
        OPEN(UNIT=funit,FILE=TRIM(FILENAME), ERR=500)
      ENDIF

 !      Call Read_Therm(PATH, 'N2', Thigh, Tlow, Tcom, Ahigh, Alow, Hf298oR)
      Call Read_Therm('BURCAT',PATH, SPECIES, Thigh, Tlow, Tcom, MW, Ahigh, &
         Alow, Hf298oR, IER)
      IF(IER /= 0) GOTO 200

      print *, SPECIES
      print *, Thigh, Tlow, Tcom, MW, Hf298oR*1.987207

!      print *, Hf298oR
!      T = 300
!      DO i = 1, 12
!        Cp1 = calc_CpoR(T, Thigh, Tlow, Tcom, Ahigh, Alow)*1.987207
!        T = T + 100
!        print *, T, Cp1
!      ENDDO

!      Cp1 = calc_CpoR(8D2, Thigh, Tlow, Tcom, Ahigh, Alow)*1.987207
!      h1 = calc_H0oR(4D2, Thigh, Tlow, Tcom, Ahigh, Alow)*1.987207
!      h2 = calc_H0oR(12D2, Thigh, Tlow, Tcom, Ahigh, Alow)*1.987207
      print *, Cp1, h1, h2
      CLOSE(UNIT=funit)
      ERROR STOP
200   PRINT *, 'READ_Therm_tester: Species ', &
         TRIM(SPECIES), ' not found in Database!'
      ERROR STOP
500   PRINT *, 'READ_Therm_tester: Cannot Open file ', TRIM(THERM), '!'
      PRINT *, 'Check path or copy mfix/model/thermochemical/', &
         TRIM(THERM), ' into run directory'
      ERROR STOP
      END Subroutine READ_Therm_tester


!     GET_LINE Subroutine
!
!      Sets line_string to the substring of data_str beginning at
!      start_index, and ending at the next newline (or end of string)

      SUBROUTINE get_line(line_string, data_str, start_index)

        CHARACTER(len=*), intent(out) :: line_string
        CHARACTER(len=*), intent(in) :: data_str
        INTEGER, intent(inout) :: start_index

        INTEGER :: END_INDEX

        END_INDEX = INDEX(DATA_STR(START_INDEX:), NEW_LINE(''))
        IF (END_INDEX == 0) THEN
           LINE_STRING = DATA_STR(START_INDEX:)
        ELSE
           LINE_STRING = DATA_STR(START_INDEX:START_INDEX+END_INDEX-2)
        END IF
        START_INDEX = START_INDEX + END_INDEX
        IF (END_INDEX == 0) START_INDEX = -1

      END SUBROUTINE get_line


      !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv  C
      !                                                                        C
      !     Module name: READ_Therm()                                          C
      !     Purpose: Read Thermo coefficients from Burcat and Ruscic           C
      !     Author: M. Syamlal                                 Date: 30-SEP-05 C
      !                                                                        C
      !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  C
      !
      SUBROUTINE READ_Therm(data_format,data_str, Sp, Thigh, Tlow, Tcom, MW, Ahigh, &
         Alow, Hf298oR, IER, lM, lN, Tmiddle, Amiddle)

      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: data_format
      CHARACTER(*), INTENT(IN) :: data_str, Sp

!     holds one line in the input file
      CHARACTER(len=100) :: LINE_STRING
      CHARACTER(len=24) :: THERMO_SOURCE
      CHARACTER(len=18) :: SPECIES, ss
      CHARACTER(len=15) :: Hf298oR_string, blank15
      INTEGER :: LINE_START
      INTEGER IER
      DOUBLE PRECISION Ahigh(7), Alow(7), Hf298oR
      DOUBLE PRECISION, INTENT(OUT) :: Thigh, Tlow, Tcom, MW
      INTEGER, OPTIONAL, INTENT(IN) :: lM, lN ! phase id and index of species
      DOUBLE PRECISION, ALLOCATABLE, OPTIONAL, INTENT(INOUT), DIMENSION(:,:,:) :: Tmiddle
      DOUBLE PRECISION, ALLOCATABLE, OPTIONAL, INTENT(INOUT), DIMENSION(:,:,:,:) :: Amiddle
      LOGICAL :: moreElement
      INTEGER :: id_space, num_element, j, num_temp, old_size
      CHARACTER(len=2) :: element_name
      CHARACTER(len=100) :: sub_line

      IER = 0
      SPECIES = SP

      LINE_STRING = '                '
      BLANK15 = '               '

      line_start = 1

! Supported Thermo blocks:
! 1) One Burcat block:  The block starts with "THERMO DATA"
! 2) One Chemkin block: The block starts with "THERMO DATA (CHEMKIN)"
! 3) One Burcat block followed by a Chemkin block
!    The Burcat block starts with "THERMO DATA"
!    The Chemkin block starts with "THERMO DATA (CHEMKIN)"

      IF(trim(data_format)=='BURCAT') THEN

         DO WHILE(LINE_STRING(1:11) /= 'THERMO DATA')
            CALL get_line(line_string, data_str, line_start)
            if (line_start < 0) then
               ier = 1
               return
            end if
         END DO

! If we want to read burcat format and it is actually chemkin format, don't read it.
! Flag species as not found (ier=1) and exit
         IF(trim(data_format)=='BURCAT'.AND.LINE_STRING(13:21)=='(CHEMKIN)') then
            ier = 1
            RETURN
         endif
! Otherwise, it is burcat data format
         THERMO_SOURCE = 'BURCAT'

      ELSEIF(trim(data_format)=='CHEMKIN') THEN

         DO WHILE(LINE_STRING(1:21) /= 'THERMO DATA (CHEMKIN)')
            CALL get_line(line_string, data_str, line_start)
            if (line_start < 0) then
               ier = 1
               return
            end if
         END DO

         THERMO_SOURCE = 'CHEMKIN'
      ENDIF

! We have identified the start of a THERMO block (either BURCAT or CHEMKIN), now look for the species.
      ss = '                 '
      call trimTab(SPECIES)
      DO WHILE(TRIM(ss) /= TRIM(SPECIES))
        CALL get_line(line_string, data_str, line_start)
        IF(LINE_STRING(1:21) == 'THERMO DATA (CHEMKIN)') then
           IER = 1
           RETURN !Stop looking for species if we found another chemkin block
        endif

        if (line_start < 0) then
           ier = 1
           return
        end if
        ss = LINE_STRING(1:18)
        call trimTab(ss)

      END DO


      SELECT CASE(TRIM(THERMO_SOURCE))

         CASE('BURCAT')
! Line 1
            call get_burcat_line_1(LINE_STRING, Tlow, Thigh, MW, Tcom)

! Tcom is almost always 1000K, however there are a few species where
! this value is too high and causes a problem (e.g., liquid water).
! Therefore, set Tcom = Thigh when Thigh < Tcom.
            Tcom = min(Tcom, Thigh)
! Line 2
            CALL get_line(line_string, data_str, line_start)
            READ(line_string, FMT='(5E15.0)',ERR=300,END=300)Ahigh(1:5)

! Line 3
            CALL get_line(line_string, data_str, line_start)
            READ(line_string, FMT='(5E15.0)',ERR=300,END=300)Ahigh(6:7), Alow(1:3)

! Line 4
            CALL get_line(line_string, data_str, line_start)
            READ(line_string, FMT='(4E15.0,A15)',ERR=300,END=300)Alow(4:7), Hf298oR_string
            IF(Hf298oR_string==blank15) THEN
! Heat of formation was left blank, compute it from polynomial coefficients
               Hf298oR = calc_Hf298oR(Thigh, Tlow, Tcom, Ahigh, Alow)
            ELSE ! convert string to float
               READ(Hf298oR_string, FMT='(E15.0)',ERR=300,END=300) Hf298oR
            ENDIF


         CASE('CHEMKIN')
! Line 1
            ! Flag indicating if element information is given in the lines after the first line
            moreElement = .False.
            call get_chemkin_line_1(LINE_STRING,Tlow, Thigh, Tcom, MW, moreElement)
            ! check the element information given in multiple lines after the first line
            DO WHILE (moreElement)
               CALL get_line(line_string, data_str, line_start)
               sub_line = trim(adjustl(line_string))
               DO WHILE((sub_line .NE. "") .AND. (trim(sub_line) .NE. "&"))
                  id_space = index(sub_line, " ")
                  element_name = sub_line(1:id_space-1)
                  sub_line = adjustl(sub_line(id_space:))
                  id_space = index(sub_line, " ")
                  read(sub_line(1:id_space), *) num_element
                  sub_line = adjustl(sub_line(id_space:))

                  IF(num_element>0) THEN
                     CALL MAKE_UPPER_CASE(element_name,2)
                     ! `E` represents electron, which does not affect the MW of species. MW(E) = 0.0
                     IF(element_name .NE. 'E') THEN
                        DO j=1,ATOMIC_MAX
                           IF(element_name==ELEMENT(j)%SYMBOL) THEN
                              IF(ELEMENT(j)%RADIOACTIVE) THEN
                                 WRITE (ERR_MSG, 1100) TRIM(ELEMENT(j)%SYMBOL), TRIM(ELEMENT(j)%NAME), TRIM(line_string)
                                 CALL LOG_ERROR()
                              ENDIF
                              MW = MW + num_element*ELEMENT(j)%WEIGHT
                              EXIT
                           ENDIF
                           IF(j == ATOMIC_MAX) THEN
                              WRITE (ERR_MSG, 1200) TRIM(element_name), TRIM(line_string)
                              CALL LOG_ERROR()
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
               IF(INDEX(line_string, '&') .EQ. 0) moreElement = .False.
            END DO

1100     FORMAT(//1X,/1x,'From: read_thermochemical -> read_therm',/1x,'Error 1100: ', &
            'ATTEMPTING TO USE A RADIOACTIVE ELEMENT, WHOSE ATOMIC WEIGHT IS NOT DEFINED.',/ &
            ' ELEMENT SYMBOL: ', A, '. ELEMENT NAME  : ', A,/ &
            ' IN LINE: ', A)

1200     FORMAT(//1X,/1x,'From: read_thermochemical -> read_therm',/1x,'Error 1200: ', &
            'ELEMENT ', A, ', GIVEN IN LINE (', A, '), IS NOT SUPPORTED IN MFIX.')

! Tcom is almost always 1000K, however there are a few species where
! this value is too high and causes a problem (e.g., liquid water).
! Therefore, set Tcom = Thigh when Thigh < Tcom.
            Tcom = min(Tcom, Thigh)

! Line 2
            CALL get_line(line_string, data_str, line_start)
! Check if more than two temperature range is given.
! If so, the second line will start with `TEMP` followed by the temperatures
            IF(INDEX(line_string, 'TEMP') .NE. 0) THEN
                num_temp = 0
                sub_line = trim(adjustl(line_string))
                sub_line = trim(adjustl(sub_line(5:)))
                DO WHILE(sub_line .NE. "")
                    num_temp = num_temp + 1
                    id_space = index(sub_line, " ")
                    sub_line = adjustl(sub_line(id_space:))
                ENDDO
                sub_line = trim(adjustl(line_string))
                sub_line = sub_line(5:)
                IF(num_temp .EQ. 3) THEN  ! Only two temperature ranges are given
                    READ(sub_line, *) Tlow, Tcom, Thigh
                    IF(Tcom .LE. Tlow) THEN
                        WRITE (ERR_MSG, 1300) TRIM(line_string)
                        CALL LOG_ERROR()
                    ENDIF
                    IF(Thigh .LE. Tcom) THEN
                        WRITE (ERR_MSG, 1300) TRIM(line_string)
                        CALL LOG_ERROR()
                    ENDIF
                    ! Lines for the coefficient
                    CALL get_line(line_string, data_str, line_start)
                    READ(line_string, FMT='(5E15.0)',ERR=300,END=300)Ahigh(1:5)
                    CALL get_line(line_string, data_str, line_start)
                    READ(line_string, FMT='(2E15.0)',ERR=300,END=300)Ahigh(6:7)
                    CALL get_line(line_string, data_str, line_start)
                    READ(line_string, FMT='(5E15.0)',ERR=300,END=300)Alow(1:5)
                    CALL get_line(line_string, data_str, line_start)
                    READ(line_string, FMT='(2E15.0)',ERR=300,END=300)Alow(6:7)

                    ! Compute it from polynomial coefficients
                    Hf298oR = calc_Hf298oR(Thigh, Tlow, Tcom, Ahigh, Alow)
                ELSE
                    old_size = SIZE(Tmiddle, 3)
                    IF(num_temp .GT. (old_size+2)) THEN
                        CALL REAL_GROW3(Tmiddle, num_temp-2)
                        CALL REAL_GROW4(Amiddle, num_temp-3)
                        Tmiddle(:,:,old_size+1:) = UNDEFINED
                        Amiddle(:,:,:,old_size+1:) = UNDEFINED
                    ENDIF
                    READ(sub_line, *) Tlow, Tmiddle(lM,lN,1:(num_temp-2)), Thigh
                    IF(Tmiddle(lM,lN,1) .LE. Tlow) THEN
                        WRITE (ERR_MSG, 1300) TRIM(line_string)
                        CALL LOG_ERROR()
                    ENDIF
                    DO j=2,num_temp-3
                        IF(Tmiddle(lM,lN,j) .LT. Tmiddle(lM,lN,j-1)) THEN
                            WRITE (ERR_MSG, 1300) TRIM(line_string)
                            CALL LOG_ERROR()
                        ENDIF
                    ENDDO
                    IF(Thigh .LE. Tmiddle(lM,lN,num_temp-2)) THEN
                        WRITE (ERR_MSG, 1300) TRIM(line_string)
                        CALL LOG_ERROR()
                    ENDIF
                    ! Lines for the coefficient
                    CALL get_line(line_string, data_str, line_start)
                    READ(line_string, FMT='(5E15.0)',ERR=300,END=300)Ahigh(1:5)
                    CALL get_line(line_string, data_str, line_start)
                    READ(line_string, FMT='(2E15.0)',ERR=300,END=300)Ahigh(6:7)
                    DO j=2,num_temp-2
                        CALL get_line(line_string, data_str, line_start)
                        READ(line_string, FMT='(5E15.0)',ERR=300,END=300)Amiddle(1:5,lM,lN,j-1)
                        CALL get_line(line_string, data_str, line_start)
                        READ(line_string, FMT='(2E15.0)',ERR=300,END=300)Amiddle(6:7,lM,lN,j-1)
                    ENDDO
                    CALL get_line(line_string, data_str, line_start)
                    READ(line_string, FMT='(5E15.0)',ERR=300,END=300)Alow(1:5)
                    CALL get_line(line_string, data_str, line_start)
                    READ(line_string, FMT='(2E15.0)',ERR=300,END=300)Alow(6:7)

                    ! Compute it from polynomial coefficients
                    Hf298oR = calc_Hf298oR(Thigh, Tlow, Tcom, Ahigh, Alow, Tmiddle(lM,lN,:), Amiddle(:,lM,lN,:))
                ENDIF
            ELSE
                READ(line_string, FMT='(5E15.0)',ERR=300,END=300)Ahigh(1:5)
                ! Line 3
                CALL get_line(line_string, data_str, line_start)
                READ(line_string, FMT='(5E15.0)',ERR=300,END=300)Ahigh(6:7), Alow(1:3)
                ! Line 4
                CALL get_line(line_string, data_str, line_start)
                READ(line_string, FMT='(4E15.0,A15)',ERR=300,END=300)Alow(4:7), Hf298oR_string
                IF(Hf298oR_string==blank15) THEN
                ! Heat of formation was left blank, compute it from polynomial coefficients
                    Hf298oR = calc_Hf298oR(Thigh, Tlow, Tcom, Ahigh, Alow)
                ELSE ! convert string to float
                    READ(Hf298oR_string, FMT='(E15.0)',ERR=300,END=300) Hf298oR
                ENDIF
            ENDIF

         CASE DEFAULT

            ! TODO: add error message

      END SELECT

1300     FORMAT(//1X,/1x,'From: read_thermochemical -> read_therm',/1x,'Error 1300: ', &
            'Temperatures given for species with multiple temperature ranges in line ('/, &
            A, ') is not in ascending order.')

      RETURN

300   PRINT *, 'READ_Therm: Error reading coefficients for species ', &
         TRIM(LINE_STRING(1:18))
      ERROR STOP

      END SUBROUTINE READ_Therm



!**********************************************************************!
! Function: calc_CpoR                                                  !
! Purpose: Evaluate the polynomial form of the specific heat.          !
!                                                                      !
!**********************************************************************!
      PURE DOUBLE PRECISION FUNCTION  calc_CpoR(T, M, N)

! Polynomial coefficients
      use physprop, only: Ahigh  ! for T in [Tcom, Thigh]
      use physprop, only: Alow   ! for T in [Tlow, Tcom)
      use physprop, only: Thigh  ! Upper bound of use
      use physprop, only: Tlow   ! Lower bound of use
      use physprop, only: Tcom   ! Switch from low to high coeffs
      use physprop, only: Amiddle ! for temperatures between Tlow and Thigh
      use physprop, only: Tmiddle ! more than two temperature ranges

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Evaluation temperaure (K)
      DOUBLE PRECISION, intent(in) :: T
! Phase index.
      INTEGER, intent(in) :: M
! Species index.
      INTEGER, intent(in) :: N

! Local Variables:
!-----------------------------------------------------------------------
! Bounded temperature.
      DOUBLE PRECISION :: xT
      INTEGER :: i

! Bound the temperature to the valid region.
      xT = max(min(T,Thigh(M,N)),Tlow(M,N))

! Evaluate the polynomial form.
     if(Tmiddle(M,N,1) .NE. UNDEFINED) then
          if(xT<=Tmiddle(M,N,1)) then
              calc_CpoR = calc_CpoR0(xT, Alow(1:5,M,N))
          else
              DO i = 2, size(Tmiddle(M,N,:))
                  if(xT <= Tmiddle(M,N,i)) then
                     calc_CpoR = calc_CpoR0(xT, Amiddle(1:5,M,N,i-1))
                     exit
                  endif
                  calc_CpoR = calc_CpoR0(xT, Ahigh(1:5,M,N))
              ENDDO
          endif
      else
          IF(xT <= Tcom(M,N))THEN
              calc_CpoR = calc_CpoR0(xT, Alow(1:5,M,N))
          ELSE
              calc_CpoR = calc_CpoR0(xT, Ahigh(1:5,M,N))
          ENDIF
      endif

      RETURN
      END Function calc_CpoR


!**********************************************************************!
! Function: calc_CpoR0                                                 !
! Purpose: Evaluate the polynomial form of the specific heat.          !
!                                                                      !
!**********************************************************************!
      PURE DOUBLE PRECISION FUNCTION  calc_CpoR0(T, A)

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Evaluation temperaure (K)
      DOUBLE PRECISION, intent(in) :: T
! Polynomial coefficients.
      DOUBLE PRECISION, intent(in) :: A(1:5)

! Evaluate the polynomial.
      calc_CpoR0 = (((A(5)*T +A(4))*T + A(3))*T + A(2))*T + A(1)

      RETURN
      END Function calc_CpoR0


!**********************************************************************!
! Function: calc_ICpoR                                                 !
! Purpose: Integrate the polynomial form of the specific heat.         !
!                                                                      !
!**********************************************************************!
      DOUBLE PRECISION FUNCTION calc_ICpoR(T, M, N, IER)

      use physprop, only: Ahigh
      use physprop, only: Thigh
      use physprop, only: ICpoR_h
      use physprop, only: Alow
      use physprop, only: Tlow
      use physprop, only: ICpoR_l
      use physprop, only: Tcom
      use physprop, only: Amiddle ! for temperatures between Tlow and Thigh
      use physprop, only: Tmiddle ! more than two temperature ranges
      use physprop, only: ICpoR_m

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Evaluation temperaure (K)
      DOUBLE PRECISION, intent(in) :: T
! Phase index.
      INTEGER, intent(in) :: M
! Species index.
      INTEGER, intent(in) :: N
! Error Flag.
      INTEGER, intent(inout) :: IER

! Local Variables.
!-----------------------------------------------------------------------
      DOUBLE PRECISION :: xT
      INTEGER :: i
!-----------------------------------------------------------------------

! Initialize the bounded temperature and error flag.
      xT = T
      IER = 0

! Verify that the temperature is in a valid range.
      if(T > Thigh(M,N)) THEN
        xT = Thigh(M,N)
      elseif(T < Tlow(M,N)) THEN
        xT = Tlow(M,N)
      endif

! Integrate the polynomial from 0.0 to T.
     if(Tmiddle(M,N,1) .NE. UNDEFINED) then
          if(xT<=Tmiddle(M,N,1)) then
              calc_ICpoR = calc_ICpoR0(xT, Alow(1:5,M,N),  ICpoR_l(M,N))
          else
              DO i = 2, size(Tmiddle,3)
                  if(xT <= Tmiddle(M,N,i)) then
                     calc_ICpoR = calc_ICpoR0(xT, Amiddle(1:5,M,N,i-1),  ICpoR_m(M,N,i-1))
                     exit
                  endif
                  calc_ICpoR = calc_ICpoR0(xT, Ahigh(1:5,M,N),  ICpoR_h(M,N))
              ENDDO
          endif
      else
          IF(xT < Tcom(M,N))THEN
              calc_ICpoR = calc_ICpoR0(xT, Alow(1:5,M,N),  ICpoR_l(M,N))
          ELSE
              calc_ICpoR = calc_ICpoR0(xT, Ahigh(1:5,M,N),  ICpoR_h(M,N))
          ENDIF
      endif

      RETURN
      END FUNCTION calc_ICpoR


!**********************************************************************!
! Function: calc_ICpoR                                                 !
! Purpose: Integrate the polynomial form of the specific heat.         !
!                                                                      !
!**********************************************************************!
      PURE DOUBLE PRECISION FUNCTION calc_ICpoR0(T, A, REF_ICpoR)

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Evaluation temperaure (K)
      DOUBLE PRECISION, intent(in) :: T
! Polynomial coefficients.
      DOUBLE PRECISION, intent(in) :: A(1:5)
! Reference Integral
      DOUBLE PRECISION, intent(in) :: REF_ICpoR

! Local Variables.
!-----------------------------------------------------------------------
! Integral of specific heat polynomial (from 0 to T) over T
      DOUBLE PRECISION ICpoRT

!-----------------------------------------------------------------------

      ICpoRT = (((A(5)*T/5.0d0 + A(4)/4.0d0)*T + A(3)/3.0d0)*T +       &
         A(2)/2.0d0)*T + A(1)

      calc_ICpoR0 = T*ICpoRT - REF_ICpoR

      RETURN
      END FUNCTION calc_ICpoR0



!**********************************************************************!
! Function: calc_H0oR                                                  !
! Purpose: Calculate the heat of formation from the first six poly-    !
!          nomial coefficients.                                        !
!                                                                      !
! >>> This function is currently unused.                               !
!                                                                      !
!**********************************************************************!
      PURE DOUBLE PRECISION FUNCTION calc_H0oR(T, Th, Tl, Tc, Ah, Al)

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Polynomial coefficients. (High/Low)
      DOUBLE PRECISION, intent(in) :: Ah(7), Al(7)
! Temperature ranges of polynomials.
      DOUBLE PRECISION, intent(in) :: Th   ! Max temp (for Ahigh)
      DOUBLE PRECISION, intent(in) :: Tl   ! Min temp (for Alow)
      DOUBLE PRECISION, intent(in) :: Tc   ! switch from low to high
! Evaluation temperaure (K)
      DOUBLE PRECISION, intent(in) :: T

! Local Variables.
!-----------------------------------------------------------------------

      !ICp = calc_ICpoR(T, Th, Tl, Tc, Ah, Al)
      !If (T < Tc) then
      !  calc_H0oR = ICp + Al(6)
      !else
      !  calc_H0oR = ICp + Ah(6)
      !endif

      calc_H0oR = 0.0

      return
      END FUNCTION calc_H0oR

!**********************************************************************!
! Function: calc_H0oRT                                                 !
! Author: Hang Zhou                                   Date: April 2024 !
! Purpose: Evaluate the polynomial form of the enthalpy.               !
!                                                                      !
!**********************************************************************!
      DOUBLE PRECISION FUNCTION  calc_H0oRT(T, M, N)

! Polynomial coefficients
      use physprop, only: Ahigh  ! for T in [Tcom, Thigh]
      use physprop, only: Alow   ! for T in [Tlow, Tcom)
      use physprop, only: Thigh  ! Upper bound of use
      use physprop, only: Tlow   ! Lower bound of use
      use physprop, only: Tcom   ! Switch from low to high coeffs
      use physprop, only: Amiddle ! for temperatures between Tlow and Thigh
      use physprop, only: Tmiddle ! more than two temperature ranges

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Evaluation temperaure (K)
      DOUBLE PRECISION, intent(in) :: T
! Phase index.
      INTEGER, intent(in) :: M
! Species index.
      INTEGER, intent(in) :: N

! Local Variables:
!-----------------------------------------------------------------------
! Bounded temperature.
      DOUBLE PRECISION :: xT

      INTEGER :: i

! Bound the temperature to the valid region.
      xT = max(min(T,Thigh(M,N)),Tlow(M,N))
! Evaluate the polynomial form.
      if(Tmiddle(M,N,1) .NE. UNDEFINED) then
          if(xT<=Tmiddle(M,N,1)) then
              calc_H0oRT = (calc_ICpoR0(xT, Alow(1:5,M,N), ZERO) +  Alow(6, M, N))/xT
          else
              DO i = 2, size(Tmiddle,3)
                  if(xT <= Tmiddle(M,N,i)) then
                     calc_H0oRT = (calc_ICpoR0(xT, Amiddle(1:5,M,N,i-1), ZERO) +  Amiddle(6,M,N,i-1))/xT
                     exit
                  endif
                  calc_H0oRT = (calc_ICpoR0(xT, Ahigh(1:5,M,N), ZERO) +  Ahigh(6, M, N))/xT
              ENDDO
          endif
      else
          IF(xT <= Tcom(M,N))THEN
              calc_H0oRT = (calc_ICpoR0(xT, Alow(1:5,M,N), ZERO) +  Alow(6, M, N))/xT
          ELSE
              calc_H0oRT = (calc_ICpoR0(xT, Ahigh(1:5,M,N), ZERO) +  Ahigh(6, M, N))/xT
          ENDIF
      endif

      RETURN
      END Function calc_H0oRT

!**********************************************************************!
! Function: calc_Hf298oR                                               !
! Purpose: Calculate the heat of formation from the first six poly-    !
!          nomial coefficients.                                        !
!                                                                      !
!**********************************************************************!
      DOUBLE PRECISION FUNCTION calc_Hf298oR(Th, Tl, Tc, Ah, Al, Tm, Am)

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Temperature ranges of polynomials.
      DOUBLE PRECISION, intent(in) :: Th   ! Max temp (for Ahigh)
      DOUBLE PRECISION, intent(in) :: Tl   ! Min temp (for Alow)
      DOUBLE PRECISION, intent(in) :: Tc   ! switch from low to high
      DOUBLE PRECISION, OPTIONAL, intent(in), dimension(:) :: Tm
! Polynomial coefficients. (High/Low)
      DOUBLE PRECISION, intent(in) :: Ah(7), Al(7)
      DOUBLE PRECISION, OPTIONAL, intent(in), dimension(:,:) :: Am
! Evaluation temperature (298K)
      DOUBLE PRECISION :: T
      INTEGER :: i

! Local Variables.
!-----------------------------------------------------------------------

      T = 298.15D0
      if (present(Tm)) then
          if(T<=Tm(1)) then
              calc_Hf298oR = calc_ICpoR0(T, Al(1:5), ZERO) +  Al(6)
          else
              DO i = 2, size(Tm)
                  if(T <= Tm(i)) then
                     calc_Hf298oR = calc_ICpoR0(T, Am(1:5,i-1), ZERO) +  Am(6,i-1)
                     exit
                  endif
                  calc_Hf298oR = calc_ICpoR0(T, Ah(1:5), ZERO) +  Ah(6)
              ENDDO
          endif
      else
          If (T < Tc) then
             calc_Hf298oR = calc_ICpoR0(T, Al(1:5), ZERO) +  Al(6)
          else
             calc_Hf298oR = calc_ICpoR0(T, Ah(1:5), ZERO) +  Ah(6)
          endif
      endif

      return
      END FUNCTION calc_Hf298oR

!**********************************************************************!
! Function: calc_S0oR                                                  !
! Author: Hang Zhou                                   Date: April 2024 !
! Purpose: Evaluate the polynomial form of the entropy.                !
!                                                                      !
!**********************************************************************!
      PURE DOUBLE PRECISION FUNCTION  calc_S0oR(T, M, N)

! Polynomial coefficients
      use physprop, only: Ahigh  ! for T in [Tcom, Thigh]
      use physprop, only: Alow   ! for T in [Tlow, Tcom)
      use physprop, only: Thigh  ! Upper bound of use
      use physprop, only: Tlow   ! Lower bound of use
      use physprop, only: Tcom   ! Switch from low to high coeffs
      use physprop, only: Amiddle ! for temperatures between Tlow and Thigh
      use physprop, only: Tmiddle ! more than two temperature ranges

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Evaluation temperaure (K)
      DOUBLE PRECISION, intent(in) :: T
! Phase index.
      INTEGER, intent(in) :: M
! Species index.
      INTEGER, intent(in) :: N

! Local Variables:
!-----------------------------------------------------------------------
! Bounded temperature.
      DOUBLE PRECISION :: xT

      INTEGER :: i

! Bound the temperature to the valid region.
      xT = max(min(T,Thigh(M,N)),Tlow(M,N))

! Evaluate the polynomial form.
     if(Tmiddle(M,N,1) .NE. UNDEFINED) then
          if(xT<=Tmiddle(M,N,1)) then
              calc_S0oR = Alow(1, M, N)*LOG(xT) + Alow(2, M, N)*xT + Alow(3, M, N)/2.0d0*xT*xT &
	           +Alow(4, M, N)/3.0d0*xT*xT*xT + Alow(5, M, N)/4.0d0*xT*xT*xT*xT + Alow(7, M, N)
          else
              DO i = 2, size(Tmiddle(M,N,:))
                  if(xT <= Tmiddle(M,N,i)) then
                     calc_S0oR = Amiddle(1, M, N,i-1)*LOG(xT) + Amiddle(2, M, N,i-1)*xT + Amiddle(3, M, N,i-1)/2.0d0*xT*xT &
	                  +Amiddle(4, M, N,i-1)/3.0d0*xT*xT*xT + Amiddle(5, M, N,i-1)/4.0d0*xT*xT*xT*xT + Amiddle(7, M, N,i-1)
                     exit
                  endif
                  calc_S0oR = Ahigh(1, M, N)*LOG(xT) + Ahigh(2, M, N)*xT + Ahigh(3, M, N)/2.0d0*xT*xT &
	           +Ahigh(4, M, N)/3.0d0*xT*xT*xT + Ahigh(5, M, N)/4.0d0*xT*xT*xT*xT + Ahigh(7, M, N)
              ENDDO
          endif
      else
          IF(xT <= Tcom(M,N))THEN
              calc_S0oR = Alow(1, M, N)*LOG(xT) + Alow(2, M, N)*xT + Alow(3, M, N)/2.0d0*xT*xT &
	           +Alow(4, M, N)/3.0d0*xT*xT*xT + Alow(5, M, N)/4.0d0*xT*xT*xT*xT + Alow(7, M, N)
          ELSE
              calc_S0oR = Ahigh(1, M, N)*LOG(xT) + Ahigh(2, M, N)*xT + Ahigh(3, M, N)/2.0d0*xT*xT &
	           +Ahigh(4, M, N)/3.0d0*xT*xT*xT + Ahigh(5, M, N)/4.0d0*xT*xT*xT*xT + Ahigh(7, M, N)
          ENDIF
      endif

      RETURN
      END Function calc_S0oR



!**********************************************************************!
! SUBROUTINE: replaceTab                                               !
! Purpose: Replace all instances of a tab with a single space.         !
!**********************************************************************!
      SUBROUTINE replaceTab(C)

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Incoming string that will have tabs removed.
      CHARACTER(len=*) :: C

! Local Variables:
!-----------------------------------------------------------------------
! Loop counter
      INTEGER :: I

      DO I = 1, len(C)
        IF(C(I:I) == CHAR(9)) C(I:I)=' '
      ENDDO

      RETURN
      END SUBROUTINE replaceTab


!**********************************************************************!
! SUBROUTINE: trimTab                                                  !
! Purpose: Search a string for the first instance of a tab. The        !
!          location of the tab and all remaining string entries are    !
!          replaced with blank spaces.                                 !
!**********************************************************************!
      SUBROUTINE trimTab(C)

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Incoming string that will have tabs removed.
      CHARACTER(len=*) :: C

! Local Variables:
!-----------------------------------------------------------------------
! Loop counter
      INTEGER :: I
! Logical indicating that a tab was located.
      LOGICAL :: tabFound

! Initialize flag
      tabFound = .FALSE.

! Look at each entry of the string. Once a tab is located, the rest of
! the string is replaced by blank spaces.
      DO I = 1, len(C)
        IF(C(I:I) == CHAR(9) ) tabFound = .TRUE.
        if(tabFound) C(I:I)=' '
      ENDDO

      RETURN
      END SUBROUTINE trimTab

    END MODULE read_thermochemical
