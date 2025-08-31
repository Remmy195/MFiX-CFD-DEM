!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR0                                                   !
!  Purpose: This routine is called before the time loop starts and is  !
!           user-definable.  The user may insert code in this routine  !
!           or call appropriate user defined subroutines.  This        !
!           can be used for setting constants and checking errors in   !
!           data.  This routine is not called from an IJK loop, hence  !
!           all indices are undefined.                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!
      SUBROUTINE USR0

      use usr

      use mfix_pic
      use discretelement, only: MAX_PIP
      use discretelement, only: DES_POS_NEW
      use discretelement, only: DES_USR_VAR

      IMPLICIT NONE
      INTEGER :: NP
      DOUBLE PRECISION :: STAT_WT

      STAT_WT = 1.0

      DO NP=1, MAX_PIP
         !IF(PEA(NP,1) .AND. .NOT. PEA(NP,4)) THEN
            DES_USR_VAR(NP,:) = DES_POS_NEW(NP,:)
            DES_STAT_WT(NP) = STAT_WT    
         !ELSE
         !   DES_USR_VAR(NP,:) = 0.0d0
         !ENDIF

      ENDDO

      RETURN
      END SUBROUTINE USR0
