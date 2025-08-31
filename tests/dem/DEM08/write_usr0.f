!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: WRITE_USR0                                             !
!  Author:                                            Date: dd-mmm-yy  !
!                                                                      !
!  Purpose: Write initial part of user-defined output                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_USR0

      use compar, only: myPE, PE_IO
      use run, only: DESCRIPTION

      IMPLICIT NONE

      INTEGER :: M, N

      CHARACTER(len=6) :: FNAME='KE.log'

! logical used for testing is the data file already exists
      LOGICAL :: EXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: fUNIT = 2030

      IF(myPE /= PE_IO) RETURN

      INQUIRE(FILE=FNAME,EXIST=EXISTS)
      IF (.NOT.EXISTS) THEN
         OPEN(UNIT=fUNIT,FILE=FNAME,STATUS='NEW')
         !WRITE(fUNIT,"(2/,25x,A)") trim(DESCRIPTION)
         WRITE(fUNIT,"(2/,7x,'TIME',5x,'No. above midplane'&
                ,5x,'Avg. KE')")
      ELSE
         OPEN(UNIT=fUNIT,FILE=FNAME,POSITION="APPEND",STATUS='OLD')
         WRITE(fUNIT,"(2/,70('-'),/,35x,'RESTART',/,70('-'),2/)")
      ENDIF

      CLOSE(fUNIT)

      RETURN
      END SUBROUTINE WRITE_USR0

