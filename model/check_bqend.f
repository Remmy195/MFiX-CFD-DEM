#include "error.inc"

MODULE CHECK_BATCH_QUEUE_END_MOD

   use compar, only: PE_IO
   use error_manager, only: err_msg, loglevel_info, log_message, loglevel_status
   use mpi_utility, only: BCAST
   use run, only: BATCH_WALLCLOCK, CHK_BATCHQ_END, TERM_BUFFER, GET_TUNIT
   use time_mod, only: get_elapsed_time

CONTAINS

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CHECK_BATCH_QUEUE_END                                   !
!  Author: A.Gel                                                       !
!                                                                      !
!----------------------------------------------------------------------!

   SUBROUTINE CHECK_BATCH_QUEUE_END(pEXIT_FLAG)

      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: pEXIT_FLAG

      pEXIT_FLAG = pEXIT_FLAG .OR. MFIX_STOP_EXISTS() .OR. MONITOR_STEADY_STATE()

      IF (CHK_BATCHQ_END) THEN
         pEXIT_FLAG = pEXIT_FLAG .OR. BATCH_WALLTIME_REACHED()
      END IF

      call bcast (pEXIT_FLAG,PE_IO)

   END SUBROUTINE CHECK_BATCH_QUEUE_END

!----------------------------------------------------------------------!
!                                                                      !
!  Function: MFIX_STOP_EXISTS                                          !
!  Author: Mark Meredith                                               !
!                                                                      !
!  Purpose:  Check if MFIX.STOP exists                                 !
!                                                                      !
!----------------------------------------------------------------------!


   LOGICAL FUNCTION MFIX_STOP_EXISTS()

! Logical flags for halt cases.
      LOGICAL :: USER_HALT

      INQUIRE(file="MFIX.STOP", exist=USER_HALT)
      MFIX_STOP_EXISTS = USER_HALT

! Report that the user signaled halt by creating MFIX.STOP
      IF(MFIX_STOP_EXISTS) THEN
         WRITE(ERR_MSG, 1200)
         CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

1200  FORMAT(2/,19('='),' MFIX STOP SIGNAL DETECTED ',19('='),/'MFIX.',&
         'STOP file detected in run directory. Terminating MFIX.',/,   &
         /19('='),' MFIX STOP SIGNAL DETECTED ',19('='))

   END FUNCTION MFIX_STOP_EXISTS

!----------------------------------------------------------------------!
!                                                                      !
!  Function: BATCH_WALLTIME_REACHED                                    !
!  Author: Mark Meredith                                               !
!                                                                      !
!  Purpose:  Check if batch walltime limit is reached                  !
!                                                                      !
!----------------------------------------------------------------------!


   LOGICAL FUNCTION BATCH_WALLTIME_REACHED()

      IMPLICIT NONE

! Elapsed wall time, and fancy formatted buffer/batch queue times.
      DOUBLE PRECISION :: ELAPSED_TIME, FANCY_BUFF, FANCY_BATCH
! Time units for formatted output.
      CHARACTER(LEN=4) :: WT_UNIT, BF_UNIT, BC_UNIT

      ELAPSED_TIME = get_elapsed_time()

! Set flags for wall time exceeded
      BATCH_WALLTIME_REACHED = ((ELAPSED_TIME+TERM_BUFFER) >= BATCH_WALLCLOCK)

! Report that the max user wall time was reached and exit.
      IF(BATCH_WALLTIME_REACHED) THEN
         CALL GET_TUNIT(ELAPSED_TIME,WT_UNIT)
         FANCY_BUFF = TERM_BUFFER
         CALL GET_TUNIT(FANCY_BUFF, BF_UNIT)
         FANCY_BATCH = BATCH_WALLCLOCK
         CALL GET_TUNIT(FANCY_BATCH, BC_UNIT)
         WRITE(ERR_MSG, 1100) FANCY_BATCH, BC_UNIT, ELAPSED_TIME, WT_UNIT, FANCY_BUFF, BF_UNIT
         CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

1100  FORMAT(2/,15('='),' Requested time limit reached ',('='),/ &
         'Batch wall time:',3X,F9.2,1X,A,/ &
         'Elapsed wall time: ',F9.2, 1X,A,/ &
         'Term buffer:',7X,F9.2,1X,A)

   END FUNCTION BATCH_WALLTIME_REACHED

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: MONITOR_STEADY_STATE                                    !
!  Author: Jeff Dietiker                              Date: 6/10/2022  !
!                                                                      !
!  Purpose:  Check if steady state was reached in monitor data         !
!                                                                      !
!----------------------------------------------------------------------!
   LOGICAL FUNCTION MONITOR_STEADY_STATE()

   use monitor, only: monitor_ss_exit

   MONITOR_STEADY_STATE = monitor_ss_exit

! Report that the user signaled halt by creating MFIX.STOP
      IF(MONITOR_STEADY_STATE) THEN
         WRITE(ERR_MSG, 1200)
         ! CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)
         CALL LOG_STATUS()
      ENDIF

1200  FORMAT(2/,19('='),' STEADY STATE DETECTED ',19('='), &
              /'Steady state was reached based on monitor data.',&
              /'Terminating MFIX.',/,   &
         /19('='),' STEADY STATE DETECTED ',19('='))

   END FUNCTION MONITOR_STEADY_STATE

END MODULE CHECK_BATCH_QUEUE_END_MOD
