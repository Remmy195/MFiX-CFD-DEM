!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: WRITE_USR1 (L)                                         !
!  Purpose: Write user-defined output                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_USR1(L)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: L

      SELECT CASE(L)
      CASE(1); CALL WRITE_LINF_NORM
      CASE(2); CALL WRITE_INT_VOLS
      END SELECT

      RETURN
      END SUBROUTINE WRITE_USR1


!......................................................................!
!  Subroutine: WRITE_L1_NORMS                                          !
!                                                                      !
!  Purpose: Compare the current position of each particle to its       !
!  starting position. Calculate the difference and report the L-Inf    !
!  norm.                                                               !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!......................................................................!
      SUBROUTINE WRITE_LINF_NORM

      use constant
      use discretelement
      use output
      use run

      use compar, only: myPE, PE_IO
      use mpi_utility, only: GLOBAL_ALL_SUM
      use mpi_utility, only: GLOBAL_ALL_MAX

      IMPLICIT NONE

! Loop index
      INTEGER :: NP

      CHARACTER*13, PARAMETER :: FNAME = 'POST_LINF.dat'
      INTEGER, PARAMETER :: FUNIT = 2030

      DOUBLE PRECISION :: X_ERR, Y_ERR, Z_ERR, lPHI
      DOUBLE PRECISION :: L_INF

      L_INF = 0.0d0
      DO NP = 1, MAX_PIP
         !IF(PEA(NP,1) .AND. .NOT.PEA(NP,4)) THEN

            X_ERR = DES_USR_VAR(NP,1) - DES_POS_NEW(NP,1)
            Y_ERR = DES_USR_VAR(NP,2) - DES_POS_NEW(NP,2)
            Z_ERR = DES_USR_VAR(NP,3) - DES_POS_NEW(NP,3)

! Shift to account for particles on the edge of the periodic BC.
            IF(abs(X_ERR) > 0.90) X_ERR = abs(X_ERR) - 1.0
            IF(abs(Y_ERR) > 0.90) Y_ERR = abs(Y_ERR) - 1.0
            IF(abs(Z_ERR) > 0.90) Z_ERR = abs(Z_ERR) - 1.0

            L_INF = max(abs(X_ERR),abs(Y_ERR),abs(Z_ERR))

         !ENDIF
      ENDDO

      write(myPE+3000,*) L_INF

      CALL GLOBAL_ALL_MAX(L_INF)

      IF(myPE /= PE_IO) RETURN

      OPEN(UNIT=FUNIT,FILE=FNAME,POSITION="APPEND",STATUS='OLD')

      WRITE(FUNIT,"(3X,f7.4,5X,f3.0,2X,3x,g11.4)") &
         TIME, TIME/USR_DT(1), L_INF

      CLOSE(FUNIT)

      RETURN
      END SUBROUTINE WRITE_LINF_NORM

!......................................................................!
!  Subroutine: WRITE_L1_NORMS                                          !
!                                                                      !
!  Purpose: Compare the current position of each particle to its       !
!  starting position. Calculate the difference and report the max L1   !
!  norm and the sum of all L1 norms.                                   !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!  Modified: A.Vaidheeswaran (for PIC)                Date:  Aug-19    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!......................................................................!
      SUBROUTINE WRITE_INT_VOLS

      use discretelement
      use output
      use run
      use mfix_pic

      use physprop, only: D_P0

      use constant, only: PI
      use fldvar, only: EP_G
      use geometry, only: VOL, XLENGTH, YLENGTH, ZLENGTH

      use indices, only: I_OF, J_OF, K_OF

      use compar, only: myPE, PE_IO
      use compar, only: IJKStart3, IJKEnd3

      use mpi_utility, only: GLOBAL_ALL_SUM
      use functions, only: FLUID_AT, IS_ON_myPE_WOBND

      IMPLICIT NONE


      CHARACTER*12, PARAMETER :: FNAME = 'POST_VOL.dat'
      INTEGER, PARAMETER :: FUNIT = 2030


      INTEGER :: IJK, PCNT
      DOUBLE PRECISION :: VOL_SYS, VOL_FLD, VOL_SLD, ERR, STAT_WT


      VOL_FLD = 0.0d0
      DO IJK = IJKStart3, IJKEnd3
         IF(IS_ON_myPE_wobnd(I_OF(IJK), J_OF(IJK), K_OF(IJK))) THEN
            IF (FLUID_AT(IJK)) THEN
               VOL_FLD = VOL_FLD + EP_G(IJK)*VOL(IJK)
            ENDIF
         ENDIF
      ENDDO

      CALL GLOBAL_ALL_SUM(VOL_SYS)
      CALL GLOBAL_ALL_SUM(VOL_FLD)

      PCNT = PIP - IGHOST_CNT
      CALL GLOBAL_ALL_SUM(PCNT)

      IF(myPE /= PE_IO) RETURN

      VOL_SYS = XLENGTH*YLENGTH*ZLENGTH

      !STAT_WT = DES_STAT_WT(MAX_PIP)
      STAT_WT = 1.0
      VOL_SLD = PCNT*STAT_WT*(PI/6.0d0)*(D_P0(1)**3)

      ERR = abs(VOL_SYS - (VOL_FLD + VOL_SLD))/VOL_SYS

      OPEN(UNIT=FUNIT,FILE=FNAME,POSITION="APPEND",STATUS='OLD')
      WRITE(FUNIT,"(3X,F7.4,3(3x,g15.8))") TIME, VOL_SLD, VOL_FLD, ERR
      CLOSE(FUNIT)

      RETURN
      END SUBROUTINE WRITE_INT_VOLS


