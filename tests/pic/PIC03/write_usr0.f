!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR0                                             C
!  Purpose: Write initial part of user-defined output                  C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_USR0

      use compar, only: myPE, PE_IO

      IMPLICIT NONE

      IF(myPE /= PE_IO) RETURN

      CALL WRITE_DAT_HEADER('POST_LINF.dat', 1)
      CALL WRITE_DAT_HEADER('POST_VOL.dat', 2)

      RETURN

      CONTAINS

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE WRITE_DAT_HEADER(FNAME, FID)

      use particle_filter, only: DES_INTERP_ON
      use particle_filter, only: DES_INTERP_SCHEME
      use particle_filter, only: DES_INTERP_MEAN_FIELDS
      use particle_filter, only: DES_INTERP_WIDTH

      use run, only: DESCRIPTION

      use param1, only: UNDEFINED

      IMPLICIT NONE

      CHARACTER(len=*) :: FNAME
      INTEGER, INTENT(IN) :: FID ! File ID

! logical used for testing is the data file already exists
      LOGICAL :: EXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: fUNIT = 2030

      INQUIRE(FILE=FNAME,EXIST=EXISTS)
      IF (.NOT.EXISTS) THEN
         OPEN(UNIT=fUNIT,FILE=FNAME,STATUS='NEW')
         WRITE(fUNIT, 1000) trim(DESCRIPTION)
      ELSE
         OPEN(UNIT=fUNIT,FILE=FNAME,POSITION="APPEND",STATUS='OLD')
      ENDIF

      WRITE(fUNIT, 1110, ADVANCE='NO') DES_INTERP_ON
      IF(DES_INTERP_ON) THEN
         WRITE(fUNIT, 1111, ADVANCE='YES') DES_INTERP_SCHEME
      ELSE
         WRITE(fUNIT, *) '   '
      ENDIF

      WRITE(fUNIT, 1120, ADVANCE='NO') DES_INTERP_MEAN_FIELDS
      IF(DES_INTERP_WIDTH /= UNDEFINED) THEN
         WRITE(fUNIT, 1121, ADVANCE='YES') DES_INTERP_WIDTH
      ELSE
         WRITE(fUNIT, *) '   '
      ENDIF


      IF(FID == 1) WRITE(fUNIT, 1251)
      IF(FID == 2) WRITE(fUNIT, 1252)

 1000 FORMAT(2/,25x,A)

 1110 FORMAT(2/,7x,'DES_INTERP_ON =',11x,L1)
 1111 FORMAT(5x,'DES_INTERP_SCHEME = ',A)

 1120 FORMAT(7x,'DES_INTERP_MEAN_FIELDS =',2x,L1)
 1121 FORMAT(5x,'DES_INTERP_WIDTH =  ',F7.2)

 1251 FORMAT(/5x,'TIME',5X,'CYCLE',5x,'L-Inf NORM')
 1252 FORMAT(/5x,'TIME',8X,'Solid Vol',7x,'Fluid Vol',11x,'Abs Error')


      CLOSE(fUNIT)
      RETURN
      END SUBROUTINE WRITE_DAT_HEADER

      END SUBROUTINE WRITE_USR0
