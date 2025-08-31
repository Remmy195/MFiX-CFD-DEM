#include "error.inc"

MODULE SET_BC_DEM_MO_MOD

   use bc, only: bc_plane
   use bc, only: bc_x_w, bc_x_e, bc_y_s, bc_y_n, bc_z_b, bc_z_t
   use des_bc, only: dem_bcmo, dem_bcmo_map, dem_bcmo_ijk
   use des_bc, only: dem_bcmo_ijkstart, dem_bcmo_ijkend
   use desgrid, only: dg_funijk
   use desgrid, only: dg_is_on_mype_plus1layers
   use desgrid, only: iofbcpos, jofbcpos, kofbcpos
   use discretelement, only: des_neighbor_search, des_le_bc, des_le_shear_dir, des_continuum_coupled, des_periodic_walls
   use discretelement, only: des_periodic_walls_x, des_periodic_walls_y, des_periodic_walls_z
   use discretelement, only: dimn
   use error_manager, only: log_message, loglevel_error
   use funits, only: dmp_log, unit_log
   use geometry, only: do_k


CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_DEM_MO                                           !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 23-Nov-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_DEM_MO

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: BCV, BCV_I      ! BC loop counter

      INTEGER :: LC

      LOGICAL, parameter :: setDBG = .FALSE.
      LOGICAL :: dFlag

      INTEGER :: MAX_CELLS, BND1, BND2

      INTEGER, ALLOCATABLE :: LOC_DEM_BCMO_IJK(:)

      INTEGER :: I,J,K,IJK
      INTEGER :: I_w, I_e, J_s, J_n, K_b, K_t

! Initialize the data structures:
      allocate( DEM_BCMO_IJKSTART(DEM_BCMO) )
      allocate( DEM_BCMO_IJKEND(DEM_BCMO) )

      DEM_BCMO_IJKSTART = -1
      DEM_BCMO_IJKEND   = -1

      dFlag = (DMP_LOG .AND. setDBG)
      if(dFlag) write(*,"(2/,2x,'DEM outlet count: ',I4)") DEM_BCMO

! Loop over the outflow BCs to get an approximate count of the number
! of fluid cells that are adjacent to the outlet.
      MAX_CELLS = 0
      DO BCV_I = 1, DEM_BCMO
         BCV = DEM_BCMO_MAP(BCV_I)

! Set the search area to the dimensions of the inlet.
         if(dFlag) WRITE(*,"(/2x,'Adding cells for BCV: ',I3)") BCV
         SELECT CASE (BC_PLANE(BCV))
         CASE('N','S')
            BND1 = IofbcPOS(BC_X_e(BCV)) - IofbcPOS(BC_X_w(BCV))
            BND2 = KofbcPOS(BC_Z_t(BCV)) - KofbcPOS(BC_Z_b(BCV))

         CASE('E','W')
            BND1 = JofbcPOS(BC_Y_n(BCV)) - JofbcPOS(BC_Y_s(BCV))
            BND2 = KofbcPOS(BC_Z_t(BCV)) - KofbcPOS(BC_Z_b(BCV))

         CASE('T','B')
            BND1 = IofbcPOS(BC_X_e(BCV)) - IofbcPOS(BC_X_w(BCV))
            BND2 = JofbcPOS(BC_Y_n(BCV)) - JofbcPOS(BC_Y_s(BCV))
         END SELECT

         MAX_CELLS = MAX_CELLS +                                      &
            2*(BND1+1)*(BND2+1) + 2*(BND1+2) + 2*(BND2+2)

         if(dFlag) WRITE(*,"(4x,'Plane:   ',A)") BC_PLANE(BCV)
         if(dFlag) WRITE(*,"(4x,'Cells: ', I8)") (BND1+1)*(BND2+1)
      ENDDO

      if(dFlag) write(*,"(2x,'Max Cells: ',I8)") MAX_CELLS

! Allocate an array to hold the IJK values. This array should be
! more than enough to store all the IJKs.
      if(allocated(LOC_DEM_BCMO_IJK)) deallocate(LOC_DEM_BCMO_IJK)
      allocate( LOC_DEM_BCMO_IJK(MAX_CELLS) )

! Loop over the IJKs for each BC and store only the IJKs that you
! own as well as the start/end array positions for each BC.
      LC = 1
      DO BCV_I = 1, DEM_BCMO

         DEM_BCMO_IJKSTART(BCV_I) = LC
         BCV = DEM_BCMO_MAP(BCV_I)


! JFD: The logic of searching des grid cells to apply the MO BC must be reviewed.
!      The original code gives the wrong set of cells and therefore
!      particles are not deleted when the leave the MO BC plane.
!      New code uses NINT in {I,J,K}ofbcpos instead of FLOOR in
!      {I,J,K}ofpos. Code was simplified to include on layer of desgrid
!      cells below the BC plane and one layer of desgrid cells above the
!      BC plane.
!      For example the argument to FLOOR could be  36.000000000000021 and it
!      would become 36, or it could be 23.999999999999996 and become 23,
!      which leads to inconsistencies. Using NINT would give 36 and 24.        


! JFD: These are des grid cells, not fluid cells
         if(dFlag) write(*,"(/2x,'Searching for des grid cells:',I3)") BCV

         I_w = IofbcPOS(BC_X_w(BCV)); I_e = IofbcPOS(BC_X_e(BCV)) + 1
         J_s = JofbcPOS(BC_Y_s(BCV)); J_n = JofbcPOS(BC_Y_n(BCV)) + 1
         IF(DO_K) THEN
            K_b = KofbcPOS(BC_Z_b(BCV)); K_t = KofbcPOS(BC_Z_t(BCV)) + 1
         ELSE
            K_b = 1; K_t = 1
         ENDIF


         if(dFlag) then
            write(*,"(4x,'Search bounds: ')")
            write(*,"(6x,'I_w/I_e:',2(2x,I6))") I_w, I_e
            write(*,"(6x,'J_s/J_n:',2(2x,I6))") J_s, J_n
            write(*,"(6x,'K_b/K_t:',2(2x,I6))") K_b, K_t
         endif

! Store the IJKs.
         DO K = K_b, K_t
         DO J = J_s, J_n
         DO I = I_w, I_e
! Skip cells that this rank does not own or are considered dead.
            IF(.NOT.dg_is_ON_myPE_plus1layers(I,J,K))CYCLE

            IJK = DG_FUNIJK(I,J,K)
            LOC_DEM_BCMO_IJK(LC) = IJK
            LC = LC+1
         ENDDO
         ENDDO
         ENDDO


         DEM_BCMO_IJKEND(BCV_I) = LC-1

         if(dFLAG) write(*,1111) BCV, BCV_I,                           &
            DEM_BCMO_IJKSTART(BCV_I),DEM_BCMO_IJKEND(BCV_I)

      ENDDO

 1111 FORMAT(/2x,'DEM Mass Outflow:',/4x,'BC:',I4,3x,'MAP:',I4,&
         /4x,'IJKSTART:',I6,/4x,'IJKEND:  ',I6)

! Allocate the global store array. This changes across MPI ranks.
      if(allocated(DEM_BCMO_IJK)) deallocate(DEM_BCMO_IJK)
      IF(LC > 1) THEN
         allocate( DEM_BCMO_IJK(LC-1) )
         DEM_BCMO_IJK(1:LC-1) = LOC_DEM_BCMO_IJK(1:LC-1)
      ELSE
         allocate( DEM_BCMO_IJK(1) )
         DEM_BCMO_IJK(1) = LOC_DEM_BCMO_IJK(1)
      ENDIF

      deallocate(LOC_DEM_BCMO_IJK)

      RETURN
      END SUBROUTINE SET_BC_DEM_MO


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_DES_LE_BC                                        !
!                                                                      !
!  Purpose: Check/set parameters for DES Lees Edeards BC.              !
!                                                                      !
!  Author: J.Musser                                   Date: 11-DEC-13  !
!                                                                      !
!  Comments: *** DES Lees Edwards BC functionality has been lost. ***  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_LE_BC

      IMPLICIT NONE

! Lees Edwards BC functionality has been lost in current DEM code
      IF(DES_LE_BC) THEN
         IF (DES_CONTINUUM_COUPLED) THEN
            WRITE(UNIT_LOG, 1064)
             call log_error()
         ENDIF
         IF (DES_NEIGHBOR_SEARCH .NE. 4) THEN
            WRITE(UNIT_LOG, 1060)
            call log_error()
         ENDIF
! not all possible shear directions are fully coded
         IF (DIMN .EQ. 2) THEN
            IF(TRIM(DES_LE_SHEAR_DIR) .NE. 'DUDY' .AND. &
               TRIM(DES_LE_SHEAR_DIR) .NE. 'DVDX') THEN
               WRITE(UNIT_LOG, 1061)
               call log_error()
            ENDIF
         ELSEIF(DIMN.EQ.3) THEN
            IF(TRIM(DES_LE_SHEAR_DIR) .NE. 'DUDY') THEN ! .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DUDZ' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DVDX' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DVDZ' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DWDX' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DWDY') THEN
               WRITE(UNIT_LOG, 1062)
               call log_error()
            ENDIF
         ENDIF
         IF (DES_PERIODIC_WALLS) THEN
            DES_PERIODIC_WALLS = .FALSE.
            DES_PERIODIC_WALLS_X = .FALSE.
            DES_PERIODIC_WALLS_Y = .FALSE.
            DES_PERIODIC_WALLS_Z = .FALSE.
            WRITE(UNIT_LOG, 1063)
            WRITE(*,1063)
         ENDIF
      ENDIF

      RETURN

 1060 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Only the grid based search option is allowed when using',&
         'using',/10X,'Lees & Edwards BC.',/1X,70('*')/)

 1061 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Invalid option for keyword DES_LE_SHEAR_DIR. When',/10X,&
         'DIMN=2 shear options are DUDY or DVDX',/1X,70('*')/)

1062  FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Invalid option for keyword DES_LE_SHEAR_DIR. When',/10X,&
         'DIMN=3 shear options are DUDY, DUDZ, DVDX, DVDZ, DWDX or',&
         'DWDY.',/1X,70('*')/)

 1063 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: DES_PERIODIC_WALLS set to false when DES_LE_BC.',&
         /10X,'DES_LE_BC implies periodic walls, however, the ',&
         'periodicity is handled',/10X, 'independently of ',&
         'DES_PERIODIC_WALLS option and so it is shut off.',&
         /1X,70('*')/)

 1064 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES_CONTINUUM_COUPLED cannot be true when using ',&
         'DES_LE_BC.',/1X,70('*')/)

   END SUBROUTINE CHECK_DES_LE_BC

END MODULE SET_BC_DEM_MO_MOD
