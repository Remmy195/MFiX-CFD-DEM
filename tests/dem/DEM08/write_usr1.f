!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: WRITE_USR1 (L)                                         !
!  Author: J.Musser                                   Date: 14-JLY-15  !
!                                                                      !
!  Purpose: Write user-defined output                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_USR1(L)

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: L

      SELECT CASE(L)
      CASE(1); CALL COLLECT_KE_DATA()
      END SELECT

      RETURN
      END SUBROUTINE WRITE_USR1



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: COLLECT_KE_DATA                                    !
!  Author: A. Vaidheeswaran                           Date: 07-JUN-21  !
!                                                                      !
!  Purpose: Routine to collect particle KE data for pile formation     !
!           to test rolling friction model in DEM.                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COLLECT_KE_DATA

      use run, only: TIME

      use output, only: USR_Y_n
      use discretelement, only: MAX_PIP
      use discretelement, only: PMASS
      use discretelement, only: DES_POS_NEW, DES_VEL_NEW
      use mpi_utility, only: GLOBAL_SUM
      use compar, only: myPE, PE_IO
      use functions, only: is_normal

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
      CHARACTER(len=6) :: FNAME = 'KE.log'

! Local Variables
!---------------------------------------------------------------------//
      INTEGER :: LC

! Number of particles above midplane
      INTEGER :: lNUM, gNUM
! KE of particles above midplane
      DOUBLE PRECISION :: lKE, gKE
! file unit for mass data
      INTEGER, PARAMETER :: fUNIT = 2030
! formatting
      CHARACTER(LEN=30) :: FORMAT_MOD
!......................................................................!

      lNUM = 0
      lKE  = 0.0

! Count number of particles below hopper outlet
      DO LC=1, MAX_PIP
         IF(IS_NORMAL(LC) .AND. DES_POS_NEW(LC,2) >  0.15) THEN
            lNUM = lNUM + 1
            lKE = lKE + 0.5*pMASS(LC)*(DES_VEL_NEW(LC,1)**2.0 + &
                  DES_VEL_NEW(LC,2)**2.0 + DES_VEL_NEW(LC,3)**2.0)
         END IF
      ENDDO

      CALL GLOBAL_SUM(lNUM, gNUM)
      CALL GLOBAL_SUM(lKE, gKE)

      IF(myPE /= PE_IO) RETURN

! Write out the total.
      OPEN(UNIT=fUNIT, FILE=FNAME, STATUS='OLD', POSITION='APPEND')
      FORMAT_MOD = "(5x,F8.5,5x,I5,5x,F16.14)" 
      WRITE(fUNIT,FORMAT_MOD) TIME, gNUM, gKE/gNUM
      CLOSE(fUNIT)

      RETURN
      END SUBROUTINE COLLECT_KE_DATA
