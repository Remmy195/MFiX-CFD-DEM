!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR1 (L)                                         C
!  Purpose: Write user-defined output                                  C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_USR1(L)

      use discretelement
      use run
      use usr
      use compar

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: L
      INTEGER, PARAMETER :: uPos = 2030

! Open the files.
      OPEN(UNIT=uPOS,FILE='POST_POS_VEL.dat', &
        POSITION="APPEND",STATUS='OLD')

! Write the results to a file.
      WRITE(uPos,"(F15.8,5x,F15.8,5x,F15.8)") Time, &
        DES_POS_new(1,1), DES_VEL_new(1,1)
      CLOSE(uPos)

      RETURN
      END SUBROUTINE WRITE_USR1
