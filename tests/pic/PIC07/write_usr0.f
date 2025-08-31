!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: WRITE_USR0                                             !
!  Author:                                            Date: dd-mmm-yy  !
!                                                                      !
!  Purpose: Write initial part of user-defined output                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_USR0

      use fs_util, only: CREATE_DIR
      use run, only: DEM_SOLIDS, PIC_SOLIDS, TFM_SOLIDS

      IMPLICIT NONE

      CALL WRITE_PRESSURE_HEADER('Pg_HGH.csv',1)

      RETURN
      END SUBROUTINE WRITE_USR0


!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE WRITE_PRESSURE_HEADER(FNAME,LC)

      use run, only: DESCRIPTION

      use output, only: USR_X_w, USR_X_e, USR_I_w, USR_I_e
      use output, only: USR_Z_b, USR_Z_t, USR_K_b, USR_K_t
      use output, only: USR_Y_s, USR_Y_n, USR_J_s, USR_J_n

      use param1, only: ZERO, UNDEFINED

      use geometry, only: DX, DY, DZ, DO_K
      use geometry, only: XMIN, KMIN1
      use geometry, only: IMAX, JMAX, KMAX, XMIN
      use compar, only: myPE, PE_IO
      use calc_cell_mod

      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LC

      INTEGER :: lIJK(4,9)
      DOUBLE PRECISION :: lWT(4,9), lLOC(4,9)

      CHARACTER(len=*) :: FNAME

! logical used for testing is the data file already exists
      LOGICAL :: EXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: fUNIT = 2030

      DOUBLE PRECISION :: LOC(4)
      INTEGER :: LL
      CHARACTER :: LLC

      !CALL INIT_ERR_MSG("WRITE_PRESSURE_HEADER")

! Set the I/J/K values used to extract data from the simulation.
      CALL CALC_CELL (XMIN, USR_X_w(LC), DX, IMAX, USR_I_w(LC))
      CALL CALC_CELL (XMIN, USR_X_e(LC), DX, IMAX, USR_I_e(LC))
      USR_I_w(LC) = USR_I_w(LC) + 1

      CALL CALC_CELL (ZERO, USR_Y_s(LC), DY, JMAX, USR_J_s(LC))
      CALL CALC_CELL (ZERO, USR_Y_n(LC), DY, JMAX, USR_J_n(LC))
      USR_J_s(LC) = USR_J_s(LC) + 1

      IF(do_K) THEN
         CALL CALC_CELL (ZERO, USR_Z_b(LC), DZ, KMAX, USR_K_b(LC))
         CALL CALC_CELL (ZERO, USR_Z_t(LC), DZ, KMAX, USR_K_t(LC))
         USR_K_b(LC) = USR_K_b(LC) + 1
      ELSE
         USR_K_b(LC) = KMIN1
         USR_K_t(LC) = KMIN1
      ENDIF

      IF(myPE /= PE_IO) RETURN

      INQUIRE(FILE=FNAME,EXIST=EXISTS)
      IF (.NOT.EXISTS) THEN

         OPEN(UNIT=fUNIT,FILE=FNAME,STATUS='NEW')
         WRITE(fUNIT,"(2/,25x,A,/)") trim(DESCRIPTION)

         WRITE(fUNIT,"('USR_X:',2(2x,f12.4),2(2x,I4))")     &
            USR_X_w(LC), USR_X_e(LC), USR_I_w(LC), USR_I_e(LC)
         WRITE(fUNIT,"('USR_Y:',2(2x,f12.4),2(2x,I4))")     &
            USR_Y_s(LC), USR_Y_n(LC), USR_J_s(LC), USR_J_n(LC)
         WRITE(fUNIT,"('USR_Z:',2(2x,f12.4),2(2x,I4))")     &
            USR_Z_b(LC), USR_Z_t(LC), USR_K_b(LC), USR_K_t(LC)

         WRITE(fUNIT,"(2/,5x,'TIME',2(11x,A4))")'BC Vg','dPg'

      ELSE
         OPEN(UNIT=fUNIT,FILE=FNAME,POSITION="APPEND",STATUS='OLD')
         WRITE(fUNIT,"(2/,70('-'),/,35x,'RESTART',/,70('-'),2/)")
      ENDIF

      CLOSE(fUNIT)

      !CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE WRITE_PRESSURE_HEADER
