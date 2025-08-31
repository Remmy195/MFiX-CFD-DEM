! -*- f90 -*-
MODULE residual_pub

      Use param, only: dim_n, dim_m
      use residual


!                      Residual sum within a group of equations
      LOGICAL          :: GROUP_RESID
   CONTAINS

      FUNCTION GET_RESID_STRING_COUNT()
         IMPLICIT NONE
         INTEGER :: GET_RESID_STRING_COUNT

         GET_RESID_STRING_COUNT = MAX_RESID_INDEX

      END FUNCTION GET_RESID_STRING_COUNT

      RECURSIVE SUBROUTINE GET_RESID_STRING(INDEX, RES)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(OUT) :: RES
         INTEGER, INTENT(IN) :: INDEX
         RES = RESID_STRING(INDEX)
      END SUBROUTINE GET_RESID_STRING

      RECURSIVE FUNCTION GET_RESID_GRP_STRING_COUNT()
         IMPLICIT NONE
         INTEGER :: GET_RESID_GRP_STRING_COUNT

         GET_RESID_GRP_STRING_COUNT = 6     ! ?

      END FUNCTION GET_RESID_GRP_STRING_COUNT

      RECURSIVE SUBROUTINE GET_RESID_GRP_STRING(INDEX, RES)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(OUT) :: RES
         INTEGER, INTENT(IN) :: INDEX

         RES = RESID_GRP_STRING(INDEX)

      END SUBROUTINE GET_RESID_GRP_STRING

      RECURSIVE FUNCTION GET_RESID(INDEX)
         IMPLICIT NONE
         DOUBLE PRECISION :: GET_RESID
         INTEGER, INTENT(IN) :: INDEX
         INTEGER :: RI, RI2

         if (.not. allocated(resid)) then
            get_resid = 0.0
            return
         endif

         IF (INDEX < 1 .OR. INDEX > SIZE(RESID_INDEX,1)) THEN
            ! PRINT *,__FILE__," INVALID VALUE FOR INDEX ",INDEX
            GET_RESID = 0.0
            RETURN
         ENDIF
         RI = RESID_INDEX(INDEX,1)
         RI2 = RESID_INDEX(INDEX,2)
         IF (RI < 1 .OR. RI > SIZE(RESID,1)) THEN
            ! PRINT *,__FILE__," INVALID VALUE FOR RESID_INDEX 1 ",RI
            GET_RESID = 0.0
            RETURN
         ENDIF
         IF (RI < 0 .OR. RI2 > SIZE(RESID,1)) THEN
            ! PRINT *,__FILE__," INVALID VALUE FOR RESID_INDEX 2 ",RI2
            GET_RESID = 0.0
            RETURN
         ENDIF
         GET_RESID = RESID(RI,RI2)

      END FUNCTION GET_RESID

      RECURSIVE FUNCTION GET_RESID_GRP(INDEX)
         IMPLICIT NONE
         DOUBLE PRECISION :: GET_RESID_GRP
         INTEGER, INTENT(IN) :: INDEX

          GET_RESID_GRP = RESID_GRP(INDEX)

      END FUNCTION GET_RESID_GRP

   END MODULE residual_pub
