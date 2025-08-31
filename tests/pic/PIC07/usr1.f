!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR1                                                   C
!  Author:                                            Date: dd-mmm-yy  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR1

      use run,only: TIME

      use bc, only: BC_J_s
      use bc, only: BC_I_w, BC_I_e
      use bc, only: BC_K_b, BC_K_t

      use bc, only: BC_V_G
      use bc, only: BC_TYPE_ENUM, MASS_INFLOW
      use fldvar, only: V_G

      use constant, only: USR_C

      use functions, only: is_on_mype_plus2layers, JM_OF
      use functions, only: bound_funijk
      use compar, only: dead_cell_at, myPE, PE_IO

      use error_manager

      IMPLICIT NONE

      INTEGER :: BCV
      INTEGER :: lI, lJ, lK

      IF(myPE /= PE_IO) RETURN

!  Set the boundary conditions.
      BCV = 10
      IF(BC_TYPE_ENUM(BCV) /= MASS_INFLOW) THEN
         !CALL INIT_ERR_MSG('USR1')
         WRITE(ERR_MSG,"('BC ',I3,' is not a mass inflow.')") BCV
         !CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Set the BC gas velocity
      BC_V_G(BCV) = MAX(TIME*USR_C(1),0.02)

      lJ=BC_J_S(BCV)
      DO lK=BC_K_B(BCV), BC_K_T(BCV)
      DO lI=BC_I_W(BCV), BC_I_E(BCV)
         IF(.NOT.IS_ON_myPE_plus2layers(lI,lJ,lK)) CYCLE
         IF(DEAD_CELL_AT(lI,lJ,lK)) CYCLE  ! skip dead cells
         V_G(JM_OF(BOUND_FUNIJK(lI,lJ,lK))) = BC_V_G(BCV)
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE USR1
