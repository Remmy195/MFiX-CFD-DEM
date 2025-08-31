!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR0                                                   C
!  Purpose: This routine is called before the time loop starts and is  C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting constants and checking errors in   C
!           data.  This routine is not called from an IJK loop, hence  C
!           all indices are undefined.                                 C
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
!
      SUBROUTINE USR0

      use discretelement, only: DTSOLID
      use constant, only: USR_C

      use param1, only: UNDEFINED
      use error_manager

      IMPLICIT NONE

      !CALL INIT_ERR_MSG('USR0')

! Fix the DEM step size
      DTSOLID = 5.0d-06

! Verify that the velocity ramp value is defined.
      IF(USR_C(1) == UNDEFINED) THEN
         WRITE(ERR_MSG,1000)
         !CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1000 FORMAT("FATAL ERR: Constant USR_C(1) is undefined.")

      !CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE USR0
