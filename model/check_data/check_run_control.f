#include "error.inc"

MODULE CHECK_RUN_CONTROL_MOD
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_RUN_CONTROL                                       !
!  Purpose: Check the run control namelist section                     !
!                                                                      !
!  Author: P.Nicoletti                                Date: 27-NOV-91  !
!          J.Musser                                   Date: 31-JAN-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_RUN_CONTROL


! Global Variables:
!---------------------------------------------------------------------//
! New or restart
      USE run, only: RUN_TYPE
! Brief description of simulation.
      USE run, only: DESCRIPTION
! Simulation units: SI, CGS
      USE run, only: UNITS
! Simulation start/stop times.
      USE run, only: TIME, TSTOP, optflag1
! Time step size, one over time step size.
      USE run, only: DT, ODT, STEADY_STATE
      USE run, only: ishii, jackson
! Maximum, minimum time step and nice dt list
      USE run, only: DT_MAX, DT_MIN, NICE_DT, NICE_DT_LIST, DT_INDEX

! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: UNDEFINED, UNDEFINED_C
      USE param1, only: ONE, ZERO

! Global Module procedures:
!---------------------------------------------------------------------//
      USE error_manager

! Skip data check when doing preprocessing only
      USE run, only:ppo
! Check that re-indexing is called only when cut cells are used
      USE cutcell, only: cartesian_grid, re_indexing

      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//

!......................................................................!

      IF(PPO) THEN
         TSTOP = ZERO
         RETURN
      ENDIF

! Clear out the run description if not specified.
      IF (DESCRIPTION == UNDEFINED_C) DESCRIPTION = ' '

! Verify UNITS input.
      IF(UNITS == UNDEFINED_C) THEN
         WRITE(ERR_MSG,1000) 'UNITS'
         CALL LOG_ERROR()
      ELSEIF((UNITS /= 'CGS') .AND. (UNITS /= 'SI')) THEN
         WRITE(ERR_MSG,1001) 'UNITS', UNITS
         CALL LOG_ERROR()
      ENDIF

! Verify that DT is valid.
      IF (DT < ZERO) THEN
         WRITE(ERR_MSG,1002) 'DT', DT
         CALL LOG_ERROR()

! Steady-state simulation.
      ELSEIF(DT == UNDEFINED .OR. DT == ZERO) THEN
         STEADY_STATE = .TRUE.
         DT = ZERO
         ODT = ZERO
         TIME = ZERO

! Transient simulation.
      ELSE
         STEADY_STATE = .FALSE.
! Calculate reciprocal of initial timestep.
         ODT = ONE/DT
! Verify the remaining time settings.
         IF (TIME == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) 'TIME'
            CALL LOG_ERROR()

         ELSEIF (TSTOP == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) 'TSTOP'
            CALL LOG_ERROR()

         ELSEIF (TIME < ZERO) THEN
            WRITE(ERR_MSG,1002)'TIME', TIME
            CALL LOG_ERROR()

         ELSEIF (TSTOP < ZERO) THEN
            WRITE(ERR_MSG,1002) 'TSTOP', TSTOP
            CALL LOG_ERROR()
         ENDIF
      ENDIF

! Generate the initial nice dt list, calculate the dt_index
      IF(NICE_DT) THEN
         IF (DT_MAX < DT_MIN) THEN
            WRITE(ERR_MSG, 1003) 'checking MFIX run!' 
            CALL LOG_ERROR()
            1003 FORMAT('Error 1003: DT_MAX should be larger than DT_MIN,',A)
         ELSE  
            CALL CREATE_NICE_LIST(DT_MIN, DT_MAX, NICE_DT_LIST)
            CALL CREATE_DT_INDEX(DT, NICE_DT_LIST, DT_INDEX)
         ENDIF
      ENDIF

! Verify the run type.
      IF(.NOT.(RUN_TYPE=='NEW' .OR. RUN_TYPE=='RESTART_1'              &
         .OR. RUN_TYPE=='RESTART_2')) THEN
         WRITE(ERR_MSG,1001) 'RUN_TYPE', RUN_TYPE
         CALL LOG_ERROR()
      ENDIF

      CALL CHECK_TURBULENCE_MODEL


! Ishii and jackson form of governing equations cannot both be invoked
      IF (ISHII .AND. JACKSON) THEN
         WRITE(ERR_MSG,2002)
         CALL LOG_ERROR()
 2002 FORMAT('Error 2002: Cannot set both ISHII = .T. and JACKSON = ',&
             '.T.')
      ENDIF

      IF (RE_INDEXING .AND. (.NOT.CARTESIAN_GRID)) THEN
         RE_INDEXING = .FALSE.
         WRITE(ERR_MSG,3000)
         CALL LOG_WARNING()
 3000 FORMAT('Warning: Re-indexing can only be used when CARTESIAN_GRID (cut cells) is used. Re-indexing will be turned off.')
      ENDIF
      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A)

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A)

 1002 FORMAT('Error 1002: Illegal or unknown input: ',A,' = ',G14.4)

      END SUBROUTINE CHECK_RUN_CONTROL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_RUN_CONTROL                                       !
!  Purpose: Check the run control namelist section                     !
!                                                                      !
!  Author: P.Nicoletti                                Date: 27-NOV-91  !
!          J.Musser                                   Date: 31-JAN-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_TURBULENCE_MODEL

! Global Variables:
!---------------------------------------------------------------------//
! Flag: for turbulence model.
      use derived_types, only: TURBULENCE_MODEL_ENUM
      use derived_types, only: MIXING_LENGTH_ENUM
      use derived_types, only: K_EPSILON_ENUM
      use derived_types, only: NO_TURBULENCE_ENUM

      use turb, only: TURBULENCE_MODEL
      use turb, only: K_EPSILON
! Viscosity bound.
      use visc_g, only: MU_GMAX

! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: UNDEFINED, UNDEFINED_C
      USE param1, only: ONE, ZERO

! Global Module procedures:
!---------------------------------------------------------------------//
      USE error_manager

      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//

!......................................................................!

      K_EPSILON = .FALSE.

      SELECT CASE(TRIM(TURBULENCE_MODEL))
      CASE('NONE')
         TURBULENCE_MODEL_ENUM = NO_TURBULENCE_ENUM
      CASE('MIXING_LENGTH')
         TURBULENCE_MODEL_ENUM = MIXING_LENGTH_ENUM
!  Check whether MU_gmax is specified for turbulence
         IF( MU_GMAX==UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) 'MU_GMAX'
            CALL LOG_ERROR()
         ENDIF
      CASE('K_EPSILON')
         TURBULENCE_MODEL_ENUM = K_EPSILON_ENUM
         K_Epsilon = .TRUE.
!  Check whether MU_gmax is specified for turbulence
         IF( MU_GMAX==UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) 'MU_GMAX'
            CALL LOG_ERROR()
         ENDIF
      CASE DEFAULT
            WRITE(ERR_MSG, 1150)
            CALL LOG_ERROR()
      END SELECT

 1150 FORMAT('Error 1150: Unknown TURBULENCE_MODEL')

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A)

      END SUBROUTINE CHECK_TURBULENCE_MODEL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CREATE_NICE_LIST(),CREATE_DT_INDEX()                    !
!  Purpose: create a nice DT list between DT_MIN and DT_MAX            !                                             
!           create DT_INDEX                                            !
!  Author: R.Ke                                       Date: 02-NOV-23  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE CREATE_NICE_LIST(MIN, MAX, ARRAY)
! Dummy arguments
            DOUBLE PRECISION, INTENT(IN) :: MIN, MAX
            DOUBLE PRECISION, ALLOCATABLE, INTENT(OUT) :: ARRAY(:)
! Local Variables
            DOUBLE PRECISION :: magnitude
            INTEGER :: i, j, n, cnt, order_of_mag_min, order_of_mag_max

            order_of_mag_min = int(floor(dlog10(MIN)))
            order_of_mag_max = int(floor(dlog10(MAX)))

! Calculate the size of array needed based on the order of MAX and MIN
            n = 0
            DO i = order_of_mag_min, order_of_mag_max
                  magnitude = 10.0d0 ** i
                  DO j = 1, 9
                  IF (j * magnitude <= MAX + 1.0d-12 .and. j * magnitude >= MIN) THEN
                        n = n + 1
                  ENDIF
                  END DO
            END DO
            ALLOCATE(ARRAY(n))

! Create the list of nice DT
            cnt = 1
            DO i = order_of_mag_min, order_of_mag_max
                  magnitude = 10.0d0 ** i
                  DO j = 1, 9
                        IF (j * magnitude <= MAX + 1.0d-12 .and. j * magnitude >= MIN) THEN
                              ARRAY(cnt) = j * magnitude
                        cnt = cnt + 1
                        ENDIF
                  END DO
            END DO
      END SUBROUTINE CREATE_NICE_LIST

      SUBROUTINE CREATE_DT_INDEX(TIMESTEP,LIST,INDEX)
! Dummy arguments
            DOUBLE PRECISION, INTENT(IN) :: TIMESTEP
            INTEGER, INTENT(OUT) :: INDEX
            DOUBLE PRECISION, ALLOCATABLE, INTENT(IN) :: LIST(:)
! Find the closest value of DT from the list, this is the target value
            INDEX = MINLOC(ABS(LIST-TIMESTEP), DIM=1)
      END SUBROUTINE CREATE_DT_INDEX
      
END MODULE CHECK_RUN_CONTROL_MOD
