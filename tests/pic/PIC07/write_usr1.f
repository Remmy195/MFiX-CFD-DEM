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
      CASE(1); CALL GET_PRESSURE_DATA('Pg_HGH.csv',L)
      END SELECT

      RETURN
      END SUBROUTINE WRITE_USR1



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: SSCP_PRESSURE_DATA                                     !
!  Author: J.Musser                                   Date: 14-JLY-15  !
!                                                                      !
!  Purpose: Routine to collect particle data for SSCP setup.           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GET_PRESSURE_DATA(FNAME,LC)

      use output, only: USR_J_s, USR_J_n
      use run, only: TIME
      use compar, only: myPE, PE_IO
      use bc, only: BC_V_g

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
      CHARACTER(len=*), INTENT(IN) :: FNAME
      INTEGER, INTENT(IN) :: LC

! Local Variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: PgS, PgN
! file unit for heat transfer data
      INTEGER, PARAMETER :: fUNIT = 2030
!......................................................................!

      PgS = avgPg(LC, USR_J_s(LC))
      PgN = avgPg(LC, USR_J_n(LC))

      IF(myPE /= PE_IO) RETURN

! Write out the total.
      OPEN(UNIT=fUNIT, FILE=FNAME, STATUS='OLD', POSITION='APPEND')
      WRITE(fUNIT,"(4(3x,F12.4))") TIME, BC_V_g(10), PgS-PgN
      CLOSE(fUNIT)

      RETURN

      contains

!```````````````````````````````````````````````````````````````````````!
!                                                                       !
!                                                                       !
!```````````````````````````````````````````````````````````````````````!
      DOUBLE PRECISION FUNCTION avgPg(pLC, pJ)

      use fldvar, only: P_G
      use geometry, only: VOL
      use functions, only: IS_ON_myPE_OWNs
      use functions, only: FUNIJK
      use functions, only: FLUID_AT
      use mpi_utility, only: GLOBAL_SUM

      use output, only: USR_I_w, USR_I_e
      use output, only: USR_K_b, USR_K_t

      INTEGER, INTENT(IN) :: pLC, pJ

      INTEGER :: lI, lK, lIJK
      DOUBLE PRECISION :: lDATA(2), rDATA(2)

      lDATA=0.0d0
      DO lK=USR_K_b(pLC),USR_K_t(pLC)
      DO lI=USR_I_w(pLC),USR_I_e(pLC)
         IF(IS_ON_myPE_OWNS(lI,pJ,lK)) THEN
            lIJK=FUNIJK(lI,pJ,lK)
            IF(FLUID_AT(lIJK)) THEN
               lDATA(1) = lData(1) + VOL(lIJK)*P_G(lIJK)
               lDATA(2) = lData(2) + VOL(lIJK)
            ENDIF
         ENDIF
      ENDDO
      ENDDO

      rDATA=0.0d0
      CALL GLOBAL_SUM(lDATA, rDATA)

      IF(rDATA(2) > 0.0d0) THEN
         avgPg = rDATA(1)/rDATA(2)
      ELSE
         avgPg = HUGE(0.0)
      ENDIF

      END FUNCTION AvgPg

      END SUBROUTINE GET_PRESSURE_DATA
