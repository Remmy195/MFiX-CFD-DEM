#include "error.inc"

!      program test; call tester; end program test

      subroutine get_burcat_line_1_tester
      implicit none

      integer          :: i
      double precision :: value1 , value2 , value3, value4
       character(len=81)     :: line

      open (unit=10,file='BURCAT.THR',status='old')
      do i = 1,9765
        read (10,'(a80)') line
      end do

      write (*,'(1x,a80)') line

      call get_burcat_line_1(line,value1,value2,value3, value4)

      write (*,*) value1
      write (*,*) value2
      write (*,*) value3
      write (*,*) value4

      stop
      end subroutine get_burcat_line_1_tester


      subroutine get_burcat_line_1(line,value1,value2,value3, value4)

      implicit none

      character(*)    :: line
      character(len=80)     :: tokens(80)
      logical          :: bSpace
      integer          :: nTokens , start , i, tcomIndex, equalIndex
      double precision :: value1 , value2 , value3, value4

      bSpace = .false.
      nTokens = 0
      start = 1
      tcomIndex = 0
      value4 = 1000.d0
      equalIndex = 0

      do i = 1,len(line)
        if (line(i:i).eq.char(9) .or. line(i:i).eq.char(32)) then
           if (.not.bSpace) then
              bSpace = .true.
               nTokens = nTokens + 1
              tokens(nTokens) = line(start:i-1)
           if(start.eq.19 .and. INDEX(tokens(nTokens),"=").ne.0) tcomIndex=nTokens
           end if
        else
           if (bSpace) start = i
           bSpace = .false.
        end if
      end do

      if (start .ne. len(line)) nTokens = nTokens - 1
      if (tcomIndex .ne. 0) then
         equalIndex = INDEX(tokens(tcomIndex),"=")
         if (tokens(tcomIndex)(equalIndex+1:equalIndex+1) .ne. char(32) ) then
            read (tokens(tcomIndex)(equalIndex+1:), '(f16.8)') value4
         else
            read (tokens(tcomIndex+1), '(f16.8)') value4
         end if
      end if

       if (nTokens .ge. 4) then
         read (tokens(nTokens  ),'(f16.8)') value3
         read (tokens(nTokens-1),'(f16.8)',err=111) value2
         read (tokens(nTokens-2),'(f16.8)') value1
         return
 111     continue
         read (tokens(nTokens-3),'(f16.8)') value1
         read (tokens(nTokens-2),'(f16.8)') value2
       end if

       return
       end


      subroutine get_chemkin_line_1(line,Tlow, Thigh, Tcomm, MW, moreElement)

! Chemkin format taken from Appendix B of
! THE CHEMKIN THERMODYNAMIC DATA BASE
! R. J. Kee, F. M. Rupley, J. A. Miller
! SANDIA REPORT
! SANDIA-8215B. UC-4
! Unliimited Release
! Printed March 1990
! Supersedes SANDIA-8215
! Sandia National Laboratories
! Albuquerque, New Mexico 87185 and Livermore. California 94550

! APPENDIX B. DATA BASE LISTING IN THE CHEMKIN FORMAT
! Chemkin expects the thermodynamic data in a specific format, which is essentially the
! same as that used by the NASA Complex Chemical Equilibrium Program. 13 In addition
! to the fourteen fit coefficients, the data base also contains the species' name, its elemental
! composition, its electronic charge, and an indication of its phase (gas, liquid or solid). The
! data for each species requires four formatted 80 column card images. Table Bl provides
! the data format specification.
!
! TABLE B1.
! FORMAT FOR THERMODYNAMIC DATA
!
!
! Card											Card
! Number	Contents						Format		Column
!
! 1		Species name (must start in Column 1) 			18A1		1 to 18
! 		Date ( not used in the code) 				6A1		19 to 24
! 		Atomic symbols and formula 				4(2A1,I3)	25 to 44
! 		Phase of species (S, L, or G for solid, 		A1		45
! 		liquid, or gas, respectively)
! 		Low temperature 					E1O.O		46 to 55
! 		High temperature 					E1O.O		56 to 65
! 		Common temperature (if needed) 				E8.0		66 to 73
! 		(blank for default)
! 		Atomic symbols and formula (if needed) 			2A1,I3		74 to 78
! 		(blank for default)
! 		The integer "1" 					I1		80
!
! 2		Coefficients a 1 - a 5 in Eqs. (1-3), 			5(E15.0)	1 to 75
! 		for upper temperature interval
! 		The integer "2" 					I1		80
!
! 3		Coefficients a6, a7 for upper temperature interval, 	5(E15.0)	1 to 75
! 		and a1, a2, and a3 for lower
! 		The integer "3" 					I1		80
!
! 4		Coefficients a4, as, a6, a7				4(E15.0)	1 to 60
! 		for lower temperature interval
! 		The integer "4" 					I1		80

      use error_manager
      USE make_upper_case_mod, only: make_upper_case
      use param1, only: zero, undefined
      use read_thermochemical, only: ELEMENT, ATOMIC_MAX

      implicit none

      character(*)    :: line
      integer          :: i, j
      character(len=18):: species_name
      character(len=6) :: date
      character(len=2), dimension(9) :: AS
      INTEGER, dimension(9) :: AF
      double precision :: Tlow, Thigh, Tcomm
      double precision :: MW
      character(len=5) :: ASF
      integer :: numElement  ! number of element
      logical :: moreElement


! We only need to read:
! - The atomic symbols and formula to compute the molecular weight
! - The low, high and common temperatures

      AS = UNDEFINED_C
      AF = 0
      numElement = 4
! Read atomic symbols and formula (4 symbols)
      read (LINE(25:44),'(4(A2,I3))') (AS(i),AF(i), i=1,4)
      IF (len(trim(LINE(74:78))) .NE. 0) THEN
         read (LINE(74:78),'(A2,I3)') AS(5),AF(5)
         numElement = numElement + 1
      ENDIF
      IF (len(trim(LINE(81:100))) .NE. 0) THEN
         IF(INDEX(LINE(81:100), '&') .EQ. 0) THEN
            read (LINE(81:100),'(4(A2,I3))') (AS(i),AF(i), i=numElement+1,numElement+4)
            numElement = numElement + 4
         ELSE
            moreElement = .True.
         ENDIF
      ENDIF

! Compute the molecular weight, by looking up the atomic weights of each element
! Reject radioactive elements
      MW = 0.0
      DO i = 1,9
         IF(AF(i)>0) THEN
            CALL MAKE_UPPER_CASE(AS(i),2)
! `E` represents electron, which does not affect the MW of species. MW(E) = 0.0
            IF(AS(i) .NE. 'E') THEN
               DO j=1,ATOMIC_MAX
                  IF(AS(i)==ELEMENT(j)%SYMBOL) THEN
                     IF(ELEMENT(j)%RADIOACTIVE) THEN
                        WRITE (ERR_MSG, 1100) TRIM(ELEMENT(j)%SYMBOL), TRIM(ELEMENT(j)%NAME), TRIM(LINE)
                        CALL LOG_ERROR()
                     ENDIF
                     MW = MW + AF(i)*ELEMENT(j)%WEIGHT
                     EXIT
                  ENDIF
                  IF(j == ATOMIC_MAX) THEN
                     WRITE (ERR_MSG, 1200) TRIM(AS(i)), TRIM(LINE)
                     CALL LOG_ERROR()
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDDO

! Low temperature
      read (LINE(46:55),'(E10.0)') Tlow
! High temperature
      read (LINE(56:65),'(E10.0)') Thigh
! Common temperature
      read (LINE(66:73),'(E8.0)') Tcomm
      if(Tcomm==ZERO) Tcomm=1000.0D0  ! If Tcom field was left blank, use
                                      ! default value of 1000 K

1100     FORMAT(//1X,/1x,'From: get_values -> get_chemkin_line_1',/1x,'Error 1100: ', &
            'ATTEMPTING TO USE A RADIOACTIVE ELEMENT, WHOSE ATOMIC WEIGHT IS NOT DEFINED.',/ &
            ' ELEMENT SYMBOL: ', A, '. ELEMENT NAME  : ', A,/ &
            ' IN LINE: ', A)

1200     FORMAT(//1X,/1x,'From: get_values -> get_chemkin_line_1',/1x,'Error 1200: ', &
            'ELEMENT ', A, ', GIVEN IN LINE (', A, '), IS NOT SUPPORTED IN MFIX.')

      return
      end
