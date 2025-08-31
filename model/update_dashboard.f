module update_dashboard_mod
   use get_smass_mod, only: get_smass
contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: UPDATE_DASHBOARD                                       C
!  Purpose: Updates and writes dashboard file                          C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 30-JAN-09  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE UPDATE_DASHBOARD(NIT)

!-----------------------------------------------
!   Modules
!-----------------------------------------------
      USE compar
      USE dashboard
      USE leqsol
      USE machine
      USE parallel
      USE run, ONLY: get_tunit, description, dt, dt_dir, run_name, time, tstop
      USE sendrecv
      USE time_mod, only: get_elapsed_time, get_remaining_time, get_io_time
      USE vtk
      Use residual

      IMPLICIT NONE

      DOUBLE PRECISION :: Smass
      DOUBLE PRECISION :: elapsed_time, remaining_time, io_time
      INTEGER :: NIT
      CHARACTER(LEN=4)  :: elapsed_tunit, remaining_tunit, io_tunit
!     temporary array to hold time data
      INTEGER DAT(8)
      CHARACTER(LEN=10) DATE, TIME_VALUE, ZONE
      LOGICAL :: Sm_flag

      IF(myPE /= PE_IO) RETURN

! Intel Linux compiler supports this function through its portability library
      CALL DATE_AND_TIME(DATE, TIME_VALUE, ZONE, DAT)
      ID_YEAR   = DAT(1)
      ID_MONTH  = DAT(2)
      ID_DAY    = DAT(3)
      ID_HOUR   = DAT(5)
      ID_MINUTE = DAT(6)
      ID_SECOND = DAT(7)

      elapsed_time = get_elapsed_time()
      call get_tunit(elapsed_time, elapsed_tunit)

      remaining_time = get_remaining_time()
      call get_tunit(remaining_time, remaining_tunit)

      io_time = get_io_time()
      call get_tunit(io_time, io_tunit)

      Sm_flag = (TRIM(RUN_STATUS)=='In Progress...'.OR.TRIM(RUN_STATUS)=='Complete.')

      IF(Sm_flag) THEN
         CALL GET_SMASS (SMASS)
         SMMIN = DMIN1(SMASS,SMMIN)
         SMMAX = DMAX1(SMASS,SMMAX)
      ENDIF

      DTMIN = DMIN1(DT,DTMIN)
      DTMAX = DMAX1(DT,DTMAX)

!      IF(NIT==0) NIT=NIT_MIN

      NIT_MIN = MIN0(NIT,NIT_MIN)
      NIT_MAX = MAX0(NIT,NIT_MAX)

      OPEN(UNIT     =  111           , &
           FILE     = 'DASHBOARD.TXT', &
           FORM     = 'FORMATTED'    , &
           ACCESS   = 'SEQUENTIAL'   , &
           STATUS   = 'REPLACE'      , &
           ACTION   = 'WRITE')

      WRITE(111,30) '________________________________________________________________________________'
      WRITE (111,10)'| RUN_NAME             = ',trim(RUN_NAME),                                    '|'
      WRITE (111,10)'| Description          = ',trim(DESCRIPTION),                                 '|'
      WRITE (111,10)'| Run status           = ',RUN_STATUS,                                        '|'
      WRITE (111,15)'| Elapsed time         = ',ELAPSED_TIME, ELAPSED_TUNIT,                       '|'
      WRITE (111,15)'| I/O time             = ',IO_TIME, IO_TUNIT,                                 '|'
    IF(RUN_STATUS/='Complete.') THEN
      WRITE (111,15)'| Remaining time       = ',REMAINING_TIME, REMAINING_TUNIT,                   '|'
    END IF
    IF(WRITE_VTK_FILES) THEN
      WRITE(111,10) '| Latest vtu file      = ',trim(VTU_FILENAME),                                '|'
    END IF
    IF(IS_SERIAL) THEN
      WRITE (111,30)'| Serial run                                                                   |'
    ELSE IF (IS_SMP) THEN
      WRITE (111,25)'| SMP run, num_threads = ',num_threads,                                       '|'
    ELSE
      WRITE (111,25)'| DMP run, num_PEs = ', numPEs,                                               '|'
    END IF
      WRITE(111,30) '|______________________________________________________________________________|'
      WRITE(111,30) '|         |         |         |         |         |                            |'
      WRITE(111,30) '|  Name   |  Value  |   Min   |   Max   | % of max|0%       Progress       100%|'
      WRITE(111,30) '|_________|_________|_________|_________|_________|____________________________|'
      WRITE(111,40,ADVANCE='NO')' Time    ',Time,Init_Time,Tstop
    CALL WRITE_SIMPLE_PROGRESS_BAR(Time,TStop-Init_Time)
    IF(DT_DIR>0) THEN
      WRITE(111,40,ADVANCE='NO')' DT (+)  ',DT,DTMIN,DTMAX
!     WRITE(111,40,ADVANCE='NO')' DT      ',DT,DTMIN,DTMAX
    ELSE
      WRITE(111,40,ADVANCE='NO')' DT (-)  ',DT,DTMIN,DTMAX
    END IF
    CALL WRITE_SIMPLE_PROGRESS_BAR(DT,DTMAX)
    IF(Sm_flag) THEN
      WRITE(111,40,ADVANCE='NO')' Sm      ',SMASS,SMMIN,SMMAX
      CALL WRITE_SIMPLE_PROGRESS_BAR(SMASS,SMMAX)
    ELSE
      WRITE(111,30) '| Sm      |         |         |         |         |                            |'
    END IF
      WRITE(111,50,ADVANCE='NO')' NIT     ',NIT,NIT_MIN,NIT_MAX
    CALL WRITE_SIMPLE_PROGRESS_BAR(dble(NIT),dble(NIT_MAX))
    IF (RESID_INDEX(8,1) == UNDEFINED_I) THEN
         WRITE (111,55) ' Max res ',RESID_STRING(8)
    END IF
      WRITE(111,30) '|_________|_________|_________|_________|_________|____________________________|'
      WRITE(111,60) '  Last updated at: ',ID_HOUR,ID_MINUTE,ID_SECOND,' on: ',ID_MONTH,ID_DAY,ID_YEAR
10    FORMAT(A,A,T80,A)
15    FORMAT(A, F6.1, 1X, A,T80,A)
25    FORMAT(A,I6,T80,A)
30    FORMAT(A)
40    FORMAT('|',A,'|',3(E9.2,'|'))
50    FORMAT('|',A,'|',3(I9,'|'))
55    FORMAT('|',A,'|',A7,'  |         |         |         |                            |')
60    FORMAT(A,I2.2,':',I2.2,':',I2.2,A,I2.2,'/',I2.2,'/',I4)
      CLOSE(111)

      RETURN
      END SUBROUTINE UPDATE_DASHBOARD

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_SIMPLE_PROGRESS_BAR                              C
!  Purpose: Displays a progress bar on the screen                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 30-JAN-09  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_SIMPLE_PROGRESS_BAR(x,x_MAX)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param1, only: zero


      IMPLICIT NONE

      INTEGER            :: BAR_WIDTH
      CHARACTER (LEN=1)  :: BAR_CHAR
      DOUBLE PRECISION   :: BAR_RESOLUTION

      INTEGER :: PROGRESS
      INTEGER :: P
      CHARACTER (LEN=9) :: TEXT
      CHARACTER (LEN=70) :: PROGRESSBAR
      DOUBLE PRECISION :: PERCENT,PTEST
      DOUBLE PRECISION :: x,x_max


      IF(X_MAX==ZERO) THEN
         WRITE(111,5)
5        FORMAT(9X,'|',28X,'|')
         RETURN
      ENDIF

      BAR_WIDTH = 28
      BAR_CHAR = '='
      BAR_RESOLUTION = 1.0

      PERCENT  = x/x_MAX * 100.0
      PROGRESS = INT(PERCENT * BAR_WIDTH)

      WRITE(TEXT,10) PERCENT
10    FORMAT(' ',F5.1,' % ')

      DO P = 1, BAR_WIDTH
         PTEST = FLOAT(P)/FLOAT(BAR_WIDTH) * 100.0
         IF(PERCENT<PTEST-BAR_RESOLUTION) THEN
            PROGRESSBAR(P:P)= ' '
         ELSE
            PROGRESSBAR(P:P)= BAR_CHAR
         ENDIF
      ENDDO

      WRITE(111,15)TEXT,'|',TRIM(PROGRESSBAR),'|'

15    FORMAT(A,A,A28,A)

      RETURN
      END SUBROUTINE WRITE_SIMPLE_PROGRESS_BAR

end module update_dashboard_mod
