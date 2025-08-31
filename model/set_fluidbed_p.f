#include "error.inc"

MODULE SET_FLUIDBED_P_MOD
CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_FLUIDBED_P                                          C
!  Purpose: Set the pressure field inside the bed assuming a fluidized C
!           bed with gravity acting the -ve y-direction                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Modifications for including cylindrical geometry           C
!  Author: M. Syamlal                                 Date:  6-MAR-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!  Revision Number: 2                                                  C
!  Purpose: Set pressure drop for cyclic boundary condition w/         C
!           pressure drop                                              C
!  Author: M. Syamlal                                 Date: 29-APR-94  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: BC_DEFINED, BC_TYPE, IC_P_g, BC_P_g           C
!                        EP_g, MW_MIX_G, RO_g0, T_g,                   C
!                        SMAX, ROP_s,                                  C
!                        DX, DY, DZ, BFY_G, DELP_X, DELP_Y, DELP_Z,    C
!                        DO_I, DO_J, DO_K, IMIN1, KMIN1, JMIN1, IMAX1, C
!                        IMAX2, JMAX1, JMAX2, KMAX1, KMAX2             C
!  Variables modified: P_g                                             C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

   SUBROUTINE SET_FLUIDBED_P

      USE bc
      USE bodyforce
      USE compar
      USE constant
      USE discretelement
      USE eos, ONLY: EOSG
      USE fldvar
      USE functions
      USE funits
      USE geometry
      USE ic
      USE indices
      USE machine, only: start_log, end_log
      USE mpi_utility
      USE param
      USE param1
      USE physprop
      USE scales
      USE sendrecv

      IMPLICIT NONE

! indices
      INTEGER :: I, J, K, IJK, M
! Local loop counter
      INTEGER :: L
! Gas pressure at the axial location j
      DOUBLE PRECISION :: PJ
! Bed weight per unit area
      DOUBLE PRECISION :: BED_WEIGHT
! Total area of a x-z plane
      DOUBLE PRECISION :: AREA
! x-z plane area of one cell
      DOUBLE PRECISION :: dAREA
! Average pressure drop per unit length
      DOUBLE PRECISION :: DPoDX, DPoDY, DPoDZ
!-----------------------------------------------

! If any initial pressures are unspecified skip next section
! calculations.
      DO L = 1, DIMENSION_IC
         IF (IC_DEFINED(L)) THEN
            IF (IC_P_G(L) == UNDEFINED) GOTO 60
            PJ = IC_P_G(L)
         ENDIF
      ENDDO

! Here the pressure in each cell is determined from a specified pressure
! drop across the domain length. This section requires that the pressure
! is already defined in all initial condition regions (otherwise this
! section would be skipped)
! ---------------------------------------------------------------->>>
      IF (DO_I .AND. DELP_X/=UNDEFINED) THEN
         DPODX = DELP_X/XLENGTH
         PJ = PJ - DPODX*HALF*(DX(IMAX1)+DX(IMAX2))
         DO I = IMAX1, IMIN1, -1
            PJ = PJ + DPODX*HALF*(DX(I)+DX(I+1))
            DO K = KMIN1, KMAX1
               DO J = JMIN1, JMAX1
! Bound Checking
                  IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
                  IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = FUNIJK(I,J,K)
                  IF (FLUID_AT(IJK)) P_G(IJK) = SCALE_PRESSURE(PJ)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF (DO_J .AND. DELP_Y/=UNDEFINED) THEN
         DPODY = DELP_Y/YLENGTH
         PJ = PJ - DPODY*HALF*(DY(JMAX1)+DY(JMAX2))
         DO J = JMAX1, JMIN1, -1
            PJ = PJ + DPODY*HALF*(DY(J)+DY(J+1))
            DO K = KMIN1, KMAX1
               DO I = IMIN1, IMAX1
! Bound Checking
                  IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
                  IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = FUNIJK(I,J,K)
                  IF (FLUID_AT(IJK)) P_G(IJK) = SCALE_PRESSURE(PJ)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF (DO_K .AND. DELP_Z/=UNDEFINED) THEN
         DPODZ = DELP_Z/ZLENGTH
         PJ = PJ - DPODZ*HALF*(DZ(KMAX1)+DZ(KMAX2))
         DO K = KMAX1, KMIN1, -1
            PJ = PJ + DPODZ*HALF*(DZ(K)+DZ(K+1))
            DO J = JMIN1, JMAX1
               DO I = IMIN1, IMAX1
! Bound Checking
                  IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
                  IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = FUNIJK(I,J,K)
                  IF (FLUID_AT(IJK)) P_G(IJK) = SCALE_PRESSURE(PJ)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
! ----------------------------------------------------------------<<<
      GOTO 100   ! pressure in all initial condition region cells was defined

   60 CONTINUE   ! pressure in an initial condition region cell was undefined

! ---------------------------------------------------------------->>>
! Search for an outflow boundary condition where pressure is specified
      PJ = UNDEFINED
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L) .AND. BC_TYPE_ENUM(L)==P_OUTFLOW) PJ = BC_P_G(L)
      ENDDO

      IF (PJ == UNDEFINED) THEN
! either a PO was not specified and/or a PO was specified but not the
! pressure at the outlet
         IF (RO_G0 /= UNDEFINED) THEN
! If incompressible flow set P_g to zero
            DO IJK = IJKSTART3, IJKEND3
               IF (FLUID_AT(IJK)) P_G(IJK) = ZERO
            ENDDO
            GOTO 100

         ELSE   ! compressible case

! Error condition -- no pressure outflow boundary condition is specified
! if a case is compressible and pressure in any of the initial
! conditions regions is unspecified, then a PO is effectively required
! (i.e., is specifies a bc_p_g).
            WRITE (err_msg, 1000)
            call log_error()
         ENDIF
      ENDIF


! Set an approximate pressure field assuming that the pressure drop
! balances the weight of the bed, if the initial pressure-field is not
! specified
      DO J = JMAX2, JMIN1, -1

! Find the average weight per unit area over an x-z slice
         BED_WEIGHT = 0.0
         AREA = 0.0
         DO K = KMIN1, KMAX1
            DO I = IMIN1, IMAX1
               IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I,J,K)
               IF (FLUID_AT(IJK)) THEN
                  IF (COORDINATES == 'CARTESIAN') THEN
                     DAREA = DX(I)*DZ(K)
                  ELSE IF (CYLINDRICAL) THEN
                     DAREA = DX(I)*X(I)*DZ(K)
                  ENDIF
                  AREA = AREA + DAREA
                  IF (RO_G0 == UNDEFINED) THEN
                     BED_WEIGHT = BED_WEIGHT - DY(J)*BFY_G(IJK)*EP_G(IJK)*EOSG(&
                        MW_MIX_G(IJK),PJ,T_G(IJK))*DAREA
                  ELSE
                     BED_WEIGHT = BED_WEIGHT - DY(J)*BFY_G(IJK)*EP_G(IJK)*RO_G0&
                        *DAREA
                  ENDIF
! This code is turned off for DEM runs until the value of rop_s can be
! ensured valid values for a DEM run at this point in the code.
                  IF (.NOT.DISCRETE_ELEMENT) THEN
                     DO M = 1, SMAX
                        BED_WEIGHT = BED_WEIGHT - DY(J)*BFY_S(IJK,M)*ROP_S(IJK,M)*&
                           DAREA
                     ENDDO
                  ENDIF  ! end if (.not.discrete_element)
               ENDIF  ! end if (fluid_at(ijk))
            ENDDO    ! end do loop (i=imin1,imax1)
         ENDDO    ! end do loop (k=kmin1,kmax1)

! Global Sum
         call global_all_sum(bed_weight)
         call global_all_sum(area)
         IF (AREA /= 0.0) BED_WEIGHT = BED_WEIGHT/AREA

         PJ = PJ + BED_WEIGHT
         DO K = KMIN1, KMAX1
            DO I = IMIN1, IMAX1
               IF(.NOT.IS_ON_MYPE_OWNS(I,J,K)) CYCLE
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I,J,K)
               IF(FLUID_AT(IJK).AND.P_G(IJK)==UNDEFINED)P_G(IJK)=SCALE_PRESSURE(PJ)
            ENDDO    ! end do (i=imin1,imax1)
         ENDDO   ! end do (k = kmin1,kmax1)
      ENDDO   ! end do (j=jmax2,jimn1, -1)
! end setting an undefined pressure in an initial condition region
! ----------------------------------------------------------------<<<

  100 CONTINUE

      call send_recv(P_G,2)

      CALL SET_BC0_OUTFLOW_P

      RETURN

 1000 FORMAT('From: SET_FLUIDBED_P'/'Message: Outflow ',&
         'pressure boundary condition (P_OUTFLOW) not found.',/&
         'All the initial pressures (IC_P_g) or at least one P_OUTFLOW',/&
         'condition need to be specified',/1X,/)

   END SUBROUTINE SET_FLUIDBED_P

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: set_bc0_outflow_p                                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

   SUBROUTINE SET_BC0_OUTFLOW_P

   use bc, only: bc_plane
   use bc, only: bc_k_b, bc_k_t
   use bc, only: bc_j_s, bc_j_n
   use bc, only: bc_i_w, bc_i_e
   use bc, only: bc_p_g, bc_type_enum, bc_defined
   use bc, only: p_outflow, mass_outflow, outflow
   use bc, only: bc_massflow_g

   use compar, only: dead_cell_at, mype, numpes

   use functions, only: is_on_mype_plus2layers, bound_funijk
   use fldvar, only: p_g
   use indices, only: im1, ip1, jm1, jp1, km1, kp1

   use param, only: DIMENSION_BC, dim_m
   use param1, only: undefined, zero
   use physprop, only: smax, mmax

      IMPLICIT NONE

! Local variables
!--------------------------------------------------------------------//
! Local index for boundary condition
      INTEGER ::  BCV
! indices
      INTEGER :: I, J, K, IJK, FIJK
! flags whether mo has non-zero massflow/volflow
      LOGICAL :: BC_FLAGG

!--------------------------------------------------------------------//
      DO BCV = 1, DIMENSION_BC
        IF (BC_DEFINED(BCV)) THEN
          IF(BC_TYPE_ENUM(BCV)==MASS_OUTFLOW .OR. &
             BC_TYPE_ENUM(BCV)==OUTFLOW) THEN

             BC_FLAGG = .FALSE.
             IF(BC_TYPE_ENUM(BCV) == MASS_OUTFLOW) THEN
               IF(BC_MASSFLOW_G(BCV) /= UNDEFINED .AND. &
                  BC_MASSFLOW_G(BCV) /= ZERO) THEN
                  BC_FLAGG = .TRUE.
               ENDIF
             ENDIF

             DO K = BC_K_B(BCV), BC_K_T(BCV)
             DO J = BC_J_S(BCV), BC_J_N(BCV)
             DO I = BC_I_W(BCV), BC_I_E(BCV)
                IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                IJK = BOUND_FUNIJK(I,J,K)

                SELECT CASE (TRIM(BC_PLANE(BCV)))
                CASE ('W')   ! fluid cell at west
                  FIJK = BOUND_FUNIJK(IM1(I),J,K)
                CASE ('E')
                  FIJK = BOUND_FUNIJK(IP1(I),J,K)
                CASE ('S')   ! fluid cell at south
                  FIJK = BOUND_FUNIJK(I,JM1(J),K)
                CASE ('N')
                  FIJK = BOUND_FUNIJK(I,JP1(J),K)
                CASE ('B')   ! fluid cell at bottom
                  FIJK = BOUND_FUNIJK(I,J,KM1(K))
                CASE ('T')
                  FIJK = BOUND_FUNIJK(I,J,KP1(K))
                END SELECT

! do not set pressure if MO with mass flow specified and user specified
! pressure was given
                IF (.NOT.(BC_FLAGG .AND. BC_P_G(BCV) /= UNDEFINED)) &
                   P_G(IJK) = P_G(FIJK)

             ENDDO
             ENDDO
             ENDDO
          ENDIF
        ENDIF
      ENDDO


   END SUBROUTINE SET_BC0_OUTFLOW_P


END MODULE SET_FLUIDBED_P_MOD
