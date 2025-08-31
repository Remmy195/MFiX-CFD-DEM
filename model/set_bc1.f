MODULE SET_BC1_MOD


CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_BC1                                                 C
!  Purpose: Set transient flow boundary conditions                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC1

! Modules
!---------------------------------------------------------------------//

      use bc
      use param, only: dimension_bc
      use set_outflow_mod, only: set_outflow

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! index for boundary condition
      INTEGER :: L


! Set the boundary conditions
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

            SELECT CASE(BC_TYPE_ENUM(L))
            CASE (P_OUTFLOW)
               CALL SET_OUTFLOW(L)
            CASE (MASS_OUTFLOW)
               CALL SET_OUTFLOW(L)
            CASE (MASS_INFLOW)

            CASE (P_INFLOW)
               CALL SET_OUTFLOW(L)
            CASE (OUTFLOW)
               CALL SET_OUTFLOW(L)
            END SELECT
         ENDIF   ! end if (bc_defined(l))
      ENDDO    ! end do loop (l=1,dimension_bc)

      RETURN
      END SUBROUTINE SET_BC1

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_BC1_DT                                              C
!  Purpose: Unlike set_bc1 this routine is called only once every time C
!  step.  It is currently used for the following:                      C
!                                                                      C
!  1) Setting ramp flow boundary conditions for mass_inflow boundaries C
!     no other bc types are supported for this.                        C
!                                                                      C
!  2) Calculate a running sum of the current mass /volume in outflow   C
!     boundaries based on stable time step values (previous version    C
!     summed values multiple times every iteration). Note that this    C
!     sum is really off by a time step as it should technically be     C
!     done at the end of a successful iteration.                       C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC1_DT

! Modules
!---------------------------------------------------------------------//

  use bc
  use param, only: dimension_bc
  use param1, only: undefined
  use calc_outflow_mod, only: calc_outflow


      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! index for boundary condition
      INTEGER :: L


! Set the boundary conditions
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

            SELECT CASE(BC_TYPE_ENUM(L))
            CASE (P_OUTFLOW)
               IF (BC_DT_0(L) /= UNDEFINED) THEN
                  CALL CALC_OUTFLOW(L)
                  CALL SET_BC1_REPORT_OUTFLOW(L)
               ENDIF

            CASE (MASS_OUTFLOW)
               CALL CALC_OUTFLOW(L)
               CALL SET_BC1_ADJUST_OUTFLOW(L)
            CASE (MASS_INFLOW)
               CALL SET_BC1_JET(L)
               CALL SET_BC1_RAMP(L)

            CASE (P_INFLOW)

            CASE (OUTFLOW)
               IF (BC_DT_0(L) /= UNDEFINED) THEN
                  CALL CALC_OUTFLOW(L)
                  CALL SET_BC1_REPORT_OUTFLOW(L)
               ENDIF
            END SELECT
         ENDIF   ! end if (bc_defined(l))
      ENDDO    ! end do loop (l=1,dimension_bc)

      RETURN
      END SUBROUTINE SET_BC1_DT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc1_jet                                             C
!  Purpose: update transient jet conditions                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC1_JET(BCV)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_plane
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_i_w, bc_i_e
      use bc, only: bc_jet
      use bc, only: bc_jet_g
      use bc, only: bc_time
      use bc, only: bc_jet_gh, bc_dt_h
      use bc, only: bc_jet_gl, bc_dt_l

      use fldvar, only: u_g, v_g, w_g

      use run, only: time, dt

      use param1, only: undefined

      use functions, only: im_of, jm_of, km_of
      use functions, only: is_on_mype_plus2layers
      use functions, only: funijk
      use compar, only: dead_cell_at

      IMPLICIT NONE
! Dummy arguments
!---------------------------------------------------------------------//
! index for boundary condition
      INTEGER, INTENT(IN) :: BCV

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: I, J, K, IJK
! IJK index for setting velocity bc
      INTEGER :: IJK2
!---------------------------------------------------------------------//

      IF(.NOT.BC_JET(BCV)) RETURN

      IF (TIME + 0.1d0*DT>=BC_TIME(BCV) .AND. &
          BC_JET_G(BCV)/=UNDEFINED) THEN

         IF (BC_JET_G(BCV) == BC_JET_GH(BCV)) THEN
            BC_JET_G(BCV) = BC_JET_GL(BCV)
            BC_TIME(BCV) = TIME + BC_DT_L(BCV)
         ELSEIF (BC_JET_G(BCV) == BC_JET_GL(BCV)) THEN
            BC_JET_G(BCV) = BC_JET_GH(BCV)
            BC_TIME(BCV) = TIME + BC_DT_H(BCV)
         ELSE
            BC_JET_G(BCV) = BC_JET_GH(BCV)
            BC_TIME(BCV) = TIME + BC_DT_H(BCV)
         ENDIF

         DO K = BC_K_B(BCV), BC_K_T(BCV)
         DO J = BC_J_S(BCV), BC_J_N(BCV)
         DO I = BC_I_W(BCV), BC_I_E(BCV)
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I,J,K)
! Why is the velocity of the boundary cell not always set (in case of
! w, s, or b plane then velocity of the adjacent fluid cell is set)?
! It should not really matter for MI...
            SELECT CASE (TRIM(BC_PLANE(BCV)))
            CASE ('W')
               IJK2 = IM_OF(IJK)
               U_G(IJK2) = BC_JET_G(BCV)
            CASE ('E')
               U_G(IJK) = BC_JET_G(BCV)
            CASE ('S')
               IJK2 = JM_OF(IJK)
               V_G(IJK2) = BC_JET_G(BCV)
            CASE ('N')
               V_G(IJK) = BC_JET_G(BCV)
            CASE ('B')
               IJK2 = KM_OF(IJK)
               W_G(IJK2) = BC_JET_G(BCV)
            CASE ('T')
               W_G(IJK) = BC_JET_G(BCV)
            END SELECT
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SET_BC1_JET


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc1_report_outflow                                  C
!  Purpose: print out outflow conditions                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC1_REPORT_OUTFLOW(BCV)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_dt_0, bc_time
      use bc, only: bc_out_n
      use bc, only: bc_mout_g, bc_mout_s
      use bc, only: bc_vout_g, bc_vout_s
      use funits, only: dmp_log, unit_log
      use param1, only: undefined, zero
      use physprop, only: smax
      use machine, only: start_log, end_log
      use run, only: time, dt, tstop
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! index for boundary condition
      INTEGER, INTENT(IN) :: BCV

! Local variables
!---------------------------------------------------------------------//
! Solids phase index
      INTEGER :: M

      IF (BC_DT_0(BCV) == UNDEFINED) RETURN

! Calculate and accumulate the actual mass and volume outflow
      IF (TIME + 0.1d0*DT>=BC_TIME(BCV) .OR. &
          TIME+0.1d0*DT>=TSTOP) THEN
         BC_TIME(BCV) = TIME + BC_DT_0(BCV)

! Average and print out the flow rates
         BC_MOUT_G(BCV) = ABS(BC_MOUT_G(BCV))/BC_OUT_N(BCV)
         BC_VOUT_G(BCV) = ABS(BC_VOUT_G(BCV))/BC_OUT_N(BCV)
         CALL START_LOG
         IF(DMP_LOG)WRITE (UNIT_LOG, 1000) BCV, TIME
         IF(DMP_LOG)WRITE (UNIT_LOG, 1100) BC_MOUT_G(BCV), &
            BC_VOUT_G(BCV)
         BC_MOUT_G(BCV) = ZERO
         BC_VOUT_G(BCV) = ZERO
         DO M = 1, SMAX
            BC_MOUT_S(BCV,M) = ABS(BC_MOUT_S(BCV,M))/BC_OUT_N(BCV)
            BC_VOUT_S(BCV,M) = ABS(BC_VOUT_S(BCV,M))/BC_OUT_N(BCV)
            IF(DMP_LOG)WRITE (UNIT_LOG, 1200) M, &
               BC_MOUT_S(BCV,M), BC_VOUT_S(BCV,M)
            BC_MOUT_S(BCV,M) = ZERO
            BC_VOUT_S(BCV,M) = ZERO
         ENDDO
         CALL END_LOG
         BC_OUT_N(BCV) = 0
      ENDIF

 1000 FORMAT(/,1X,'Average outflow rates at BC No. ',I2,'  At t = ',G12.5)
 1100 FORMAT(3X,'Gas : Mass flow = ',G12.5,'     Volumetric flow = ',G12.5)
 1200 FORMAT(3X,'Solids-',I1,' : Mass flow = ',G12.5,&
         '     Volumetric flow = ',G12.5)

      RETURN
      END SUBROUTINE SET_BC1_REPORT_OUTFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc1_adjust_outflow                                  C
!  Purpose: Adjust velocities to get specified mass or volumetric      C
!  flow rate based on average outflow rate                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC1_ADJUST_OUTFLOW(BCV)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_dt_0, bc_time
      use bc, only: bc_i_w, bc_i_e
      use bc, only: bc_j_s, bc_j_n
      use bc, only: bc_k_b, bc_k_t
      use bc, only: bc_massflow_g, bc_massflow_s
      use bc, only: bc_mout_g, bc_mout_s
      use bc, only: bc_out_n
      use bc, only: bc_plane
      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use bc, only: bc_u_s, bc_v_s, bc_w_s
      use bc, only: bc_volflow_g, bc_volflow_s
      use bc, only: bc_vout_g, bc_vout_s
      use compar, only: dead_cell_at
      use fldvar, only: u_g, v_g, w_g
      use fldvar, only: u_s, v_s, w_s, rop_s
      use functions, only: funijk
      use functions, only: im_of, jm_of, km_of
      use functions, only: is_on_mype_plus2layers
      use funits, only: dmp_log, unit_log
      use machine, only: start_log, end_log
      use param1, only: undefined, zero, one, small_number
      use physprop, only: smax, mmax
      use run, only: kt_type_enum, ghd_2007
      use run, only: time, dt, tstop
      use utilities, only: cap_bc_vel
      use toleranc, only: max_inlet_vel

      IMPLICIT NONE
! Dummy arguments
!---------------------------------------------------------------------//
! index for boundary condition
      INTEGER, INTENT(IN) :: BCV

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: I, J, K, IJK
! IJK index for setting velocity bc
      INTEGER :: IJK2
! Solids phase index
      INTEGER :: M
! local solids velocity for mixture (for ghd)
      DOUBLE PRECISION :: lvel_s

      DOUBLE PRECISION :: SGN
!---------------------------------------------------------------------//

! Calculate and accumulate the actual mass and volume outflow
      IF (TIME + 0.1d0*DT>=BC_TIME(BCV) .OR. &
          TIME+0.1d0*DT>=TSTOP) THEN
         BC_TIME(BCV) = TIME + BC_DT_0(BCV)

! Average and print out the flow rates
         BC_MOUT_G(BCV) = ABS(BC_MOUT_G(BCV))/BC_OUT_N(BCV)
         BC_VOUT_G(BCV) = ABS(BC_VOUT_G(BCV))/BC_OUT_N(BCV)
         CALL START_LOG
         IF(DMP_LOG)WRITE (UNIT_LOG, 1000) BCV, TIME
         IF(DMP_LOG)WRITE (UNIT_LOG, 1100) BC_MOUT_G(BCV), &
            BC_VOUT_G(BCV)
         DO M = 1, SMAX
            BC_MOUT_S(BCV,M) = ABS(BC_MOUT_S(BCV,M))/BC_OUT_N(BCV)
            BC_VOUT_S(BCV,M) = ABS(BC_VOUT_S(BCV,M))/BC_OUT_N(BCV)
            IF(DMP_LOG)WRITE (UNIT_LOG, 1200) M, &
               BC_MOUT_S(BCV,M), BC_VOUT_S(BCV,M)
         ENDDO
         CALL END_LOG
         BC_OUT_N(BCV) = 0

         SGN = ONE
! fluid at east, north, top (boundary w, s, b)
         SELECT CASE (BC_PLANE(BCV))
         CASE ('E'); SGN = -SGN
         CASE ('N'); SGN = -SGN
         CASE ('T'); SGN = -SGN
         END SELECT

! Now that we know the mass and volume outflow update the bc velocities
! (gas phase)
         IF (BC_MASSFLOW_G(BCV) /= UNDEFINED) THEN
            IF (BC_MOUT_G(BCV) > SMALL_NUMBER) THEN
               SELECT CASE (TRIM(BC_PLANE(BCV)))
               CASE ('W', 'E')
                  BC_U_G(BCV) = BC_U_G(BCV)*BC_MASSFLOW_G(BCV)/&
                     BC_MOUT_G(BCV)
                  BC_U_G(BCV) = CAP_BC_VEL(BC_U_G(BCV), SGN)
               CASE ('S', 'N')
                  BC_V_G(BCV) = BC_V_G(BCV)*BC_MASSFLOW_G(BCV)/&
                     BC_MOUT_G(BCV)
                  BC_V_G(BCV) = CAP_BC_VEL(BC_V_G(BCV), SGN)
               CASE ('B', 'T')
                  BC_W_G(BCV) = BC_W_G(BCV)*BC_MASSFLOW_G(BCV)/&
                     BC_MOUT_G(BCV)
                  BC_W_G(BCV) = CAP_BC_VEL(BC_W_G(BCV), SGN)
               END SELECT
            ENDIF
         ELSEIF (BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN
            IF (BC_VOUT_G(BCV) > SMALL_NUMBER) THEN
               SELECT CASE (TRIM(BC_PLANE(BCV)))
               CASE ('W', 'E')
                  BC_U_G(BCV) = BC_U_G(BCV)*BC_VOLFLOW_G(BCV)/&
                     BC_VOUT_G(BCV)
                  BC_U_G(BCV) = CAP_BC_VEL(BC_U_G(BCV), SGN)
               CASE ('S', 'N')
                  BC_V_G(BCV) = BC_V_G(BCV)*BC_VOLFLOW_G(BCV)/&
                     BC_VOUT_G(BCV)
                  BC_V_G(BCV) = CAP_BC_VEL(BC_V_G(BCV), SGN)
               CASE ('B', 'T')
                  BC_W_G(BCV) = BC_W_G(BCV)*BC_VOLFLOW_G(BCV)/&
                     BC_VOUT_G(BCV)
                  BC_W_G(BCV) = CAP_BC_VEL(BC_W_G(BCV), SGN)
               END SELECT
            ENDIF
         ENDIF
! zero out counter for new cycle
         BC_MOUT_G(BCV) = zero
         BC_VOUT_G(BCV) = zero

! Now that we know the mass and volume outflow update the bc velocities
! (solids phase)
         DO M = 1, SMAX
            IF (BC_MASSFLOW_S(BCV,M) /= UNDEFINED) THEN
               IF (BC_MOUT_S(BCV,M) > SMALL_NUMBER) THEN
                  SELECT CASE (TRIM(BC_PLANE(BCV)))
                  CASE ('W', 'E')
                     BC_U_S(BCV,M) = BC_U_S(BCV,M)*&
                        BC_MASSFLOW_S(BCV,M)/BC_MOUT_S(BCV,M)
                     BC_U_S(BCV,M) = CAP_BC_VEL(BC_U_S(BCV,M), SGN)
                  CASE ('S', 'N')
                     BC_V_S(BCV,M) = BC_V_S(BCV,M)*&
                        BC_MASSFLOW_S(BCV,M)/BC_MOUT_S(BCV,M)
                     BC_V_S(BCV,M) = CAP_BC_VEL(BC_V_S(BCV,M), SGN)
                  CASE ('B', 'T')
                     BC_W_S(BCV,M) = BC_W_S(BCV,M)*&
                        BC_MASSFLOW_S(BCV,M)/BC_MOUT_S(BCV,M)
                     BC_W_S(BCV,M) = CAP_BC_VEL(BC_W_S(BCV,M), SGN)
                  END SELECT
               ENDIF
            ELSEIF (BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN
               IF (BC_VOUT_S(BCV,M) > SMALL_NUMBER) THEN
                  SELECT CASE (TRIM(BC_PLANE(BCV)))
                  CASE ('W', 'E')
                     BC_U_S(BCV,M) = BC_U_S(BCV,M)*&
                        BC_VOLFLOW_S(BCV,M)/BC_VOUT_S(BCV,M)
                     BC_U_S(BCV,M) = CAP_BC_VEL(BC_U_S(BCV,M), SGN)
                  CASE ('S', 'N')
                     BC_V_S(BCV,M) = BC_V_S(BCV,M)*&
                        BC_VOLFLOW_S(BCV,M)/BC_VOUT_S(BCV,M)
                     BC_V_S(BCV,M) = CAP_BC_VEL(BC_V_S(BCV,M), SGN)
                  CASE ('B', 'T')
                     BC_W_S(BCV,M) = BC_W_S(BCV,M)*&
                        BC_VOLFLOW_S(BCV,M)/BC_VOUT_S(BCV,M)
                     BC_W_S(BCV,M) = CAP_BC_VEL(BC_W_S(BCV,M), SGN)
                 END SELECT
               ENDIF
            ENDIF
! zero out counter for new cycle
            BC_MOUT_S(BCV,M) = zero
            BC_VOUT_S(BCV,M) = zero
         ENDDO

! Apply updated boundary velocities - Define the field variables at the
! boundaries according to user specifications with modifications from
! the above calculations.
! If the boundary plane is W, S, or B (i.e., the fluid cell is on the
! west, south or bottom of the boundary cell) then define the velocity
! of the adjacent fluid cell according to the boundary velocity rather
! than the velocity of the boundary cell.
! Why not set the velocity in the boundary cell itself?  Based on the
! momentum bc routine it should not really matter for MO.
         DO K = BC_K_B(BCV), BC_K_T(BCV)
         DO J = BC_J_S(BCV), BC_J_N(BCV)
         DO I = BC_I_W(BCV), BC_I_E(BCV)
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I,J,K)
            SELECT CASE (TRIM(BC_PLANE(BCV)))
            CASE ('W')
               IJK2 = IM_OF(IJK)
               U_G(IJK2) = BC_U_G(BCV)
            CASE ('E')
               U_G(IJK) = BC_U_G(BCV)
            CASE ('S')
               IJK2 = JM_OF(IJK)
               V_G(IJK2) = BC_V_G(BCV)
            CASE ('N')
               V_G(IJK) = BC_V_G(BCV)
            CASE ('B')
               IJK2 = KM_OF(IJK)
               W_G(IJK2) = BC_W_G(BCV)
            CASE ('T')
               W_G(IJK) = BC_W_G(BCV)
            END SELECT
            DO M = 1, SMAX
               SELECT CASE (TRIM(BC_PLANE(BCV)))
               CASE ('W')
                  IJK2 = IM_OF(IJK)
                  U_S(IJK2,M) = BC_U_S(BCV,M)
               CASE ('E')
                  U_S(IJK,M) = BC_U_S(BCV,M)
               CASE ('S')
                  IJK2 = JM_OF(IJK)
                  V_S(IJK2,M) = BC_V_S(BCV,M)
               CASE ('N')
                  V_S(IJK,M) = BC_V_S(BCV,M)
               CASE ('B')
                  IJK2 = KM_OF(IJK)
                  W_S(IJK2,M) = BC_W_S(BCV,M)
               CASE ('T')
                  W_S(IJK,M) = BC_W_S(BCV,M)
               END SELECT
            ENDDO

! compute mixutre velocity BC for GHD theory
            IF(KT_TYPE_ENUM == GHD_2007) THEN
               lvel_s = zero

! bulk density is set by set_outflow and is set in bc according to
! neighboring fluid appropriately. so we don't need to reference
! bulk density with ijk vs ijk2 index (ijk value will = ijk2 value).
               DO M = 1, SMAX
                  SELECT CASE (TRIM(BC_PLANE(BCV)))
                  CASE ('W', 'E')
                     lvel_s = lvel_s + &
                        BC_U_S(BCV,M)*ROP_S(IJK,M)
                  CASE ('S', 'N')
                     lvel_s = lvel_s + &
                        BC_V_S(BCV,M)*ROP_S(IJK,M)
                  CASE ('B', 'T')
                     lvel_s = lvel_s + &
                        BC_W_S(BCV,M)*ROP_S(IJK,M)
                  END SELECT
               ENDDO

               IF (ROP_S(IJK,MMAX) > 0) THEN
                  lvel_s = lvel_s /ROP_S(IJK,MMAX)
               ELSE
                  lvel_s = ZERO
               ENDIF

               SELECT CASE (TRIM(BC_PLANE(BCV)))
               CASE ('W')
                  IJK2 = IM_OF(IJK)
                  U_S(IJK2,MMAX) =  lvel_s
               CASE ('E')
                  U_S(IJK,MMAX) = lvel_s
               CASE ('S')
                  IJK2 = JM_OF(IJK)
                  V_S(IJK2,MMAX) = lvel_s
               CASE ('N')
                  V_S(IJK,MMAX) = lvel_s
               CASE ('B')
                  IJK2 = KM_OF(IJK)
                  W_S(IJK2,MMAX) = lvel_s
               CASE ('T')
                  W_S(IJK,MMAX) = lvel_s
               END SELECT

            ENDIF   ! end if (kt_type_enum==ghd_2007)
         ENDDO
         ENDDO
         ENDDO
      ENDIF   ! if time to update outflow condition

 1000 FORMAT(/,1X,'Average outflow rates at BC No. ',I2,'  At t = ',G12.5)
 1100 FORMAT(3X,'Gas : Mass flow = ',G12.5,'     Volumetric flow = ',G12.5)
 1200 FORMAT(3X,'Solids-',I1,' : Mass flow = ',G12.5,&
         '     Volumetric flow = ',G12.5)

      RETURN
      END SUBROUTINE SET_BC1_ADJUST_OUTFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc1_ramp                                            C
!  Purpose: ramp BC parameters                                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_BC1_RAMP(BCV)

      use set_bc0_mod
      use flow_to_vel_mod, only: flow_to_vel
      use check_bc0_flow_mod, only: check_bc_vel_inflow
      use set_bc0_flow_mod

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_jet_g
      use bc, only: bc_time
      use bc, only: bc_jet_gh, bc_dt_h
      use bc, only: bc_jet_gl, bc_dt_l


      use run, only: time, dt

      use param1, only: undefined


! Total number of (actual) continuum solids.
      use physprop, only: SMAX
! Total number of discrete solids.
      use discretelement, only: DES_MMAX
! Maximum number of disperse phases
      use param, only: DIM_M

      use param, only: dimension_bc
      use param1, only: undefined_i
      use param1, only: zero, one, undefined

      use bc
! User specified BC solids bulk density
      use bc, only: BC_ROP_s
! Solids volume fraction at BC
      use bc, only: BC_EP_s
      use bc, only: BC_EP_g
      use bc, only: cg_mi_converted_to_ps

      use cutcell, only: cartesian_grid

!
      IMPLICIT NONE
! Dummy arguments

!---------------------------------------------------------------------//
! index for boundary condition
      INTEGER, INTENT(IN) :: BCV

! Local variables
!---------------------------------------------------------------------//
! Total number of solids phases (continuum + discrete)
      INTEGER :: MMAX_TOT
      INTEGER :: M
! Flag to skip checks on indexed solid phase.
      LOGICAL :: SKIP(1:DIM_M)
! Logical on whether to perform a check on the direction (+/-; in/out)
! of any specified bc velocity for the given bc plane
      LOGICAL :: CHECK_VEL
!......................................................................!

      IF(.NOT.BC_RAMP(BCV)) RETURN

! Total number of solids.
      MMAX_TOT = SMAX! + DES_MMAX

      CHECK_VEL = .TRUE.
      IF (CARTESIAN_GRID) CHECK_VEL = .FALSE.


! ramp BC

! Gas pressure
      if(BC_P_g_ramp_t1(BCV)/=UNDEFINED)  &
         BC_P_g(BCV) = set_bc_value (time, BC_P_g_ramp_t0(BCV), BC_P_g_ramp_t1(BCV),&
                       BC_P_g_ramp_v0(BCV), BC_P_g_ramp_v1(BCV))
!
! Gas phase temperature
      if(BC_T_g_ramp_t1(BCV)/=UNDEFINED)  &
         BC_T_g(BCV) = set_bc_value (time, BC_T_g_ramp_t0(BCV), BC_T_g_ramp_t1(BCV),&
                       BC_T_g_ramp_v0(BCV), BC_T_g_ramp_v1(BCV))
!
! Solids phase temperature
      DO M = 1, SMAX
         if(BC_T_s_ramp_t1(BCV,M)/=UNDEFINED)  &
            BC_T_s(BCV,M) = set_bc_value (time, BC_T_s_ramp_t0(BCV,M), BC_T_s_ramp_t1(BCV,M),&
                          BC_T_s_ramp_v0(BCV,M), BC_T_s_ramp_v1(BCV,M))
      ENDDO

! x-component of gas velocity
      if(BC_U_g_ramp_t1(BCV)/=UNDEFINED)  &
         BC_U_g(BCV) = set_bc_value (time, BC_U_g_ramp_t0(BCV), BC_U_g_ramp_t1(BCV),&
                       BC_U_g_ramp_v0(BCV), BC_U_g_ramp_v1(BCV))
!
! x-component of solids phase velocity
      DO M = 1, SMAX
         if(BC_U_s_ramp_t1(BCV,M)/=UNDEFINED)  &
            BC_U_s(BCV,M) = set_bc_value (time, BC_U_s_ramp_t0(BCV,M), BC_U_s_ramp_t1(BCV,M),&
                          BC_U_s_ramp_v0(BCV,M), BC_U_s_ramp_v1(BCV,M))
      ENDDO
!
! y-component of gas velocity
      if(BC_V_g_ramp_t1(BCV)/=UNDEFINED)  &
         BC_V_g(BCV) = set_bc_value (time, BC_V_g_ramp_t0(BCV), BC_V_g_ramp_t1(BCV),&
                       BC_V_g_ramp_v0(BCV), BC_V_g_ramp_v1(BCV))
!
! y-component of solids phase velocity
      DO M = 1, SMAX
         if(BC_V_s_ramp_t1(BCV,M)/=UNDEFINED)  &
            BC_V_s(BCV,M) = set_bc_value (time, BC_V_s_ramp_t0(BCV,M), BC_V_s_ramp_t1(BCV,M),&
                          BC_V_s_ramp_v0(BCV,M), BC_V_s_ramp_v1(BCV,M))
      ENDDO
!
! z-component of gas velocity
      if(BC_W_g_ramp_t1(BCV)/=UNDEFINED)  &
         BC_W_g(BCV) = set_bc_value (time, BC_W_g_ramp_t0(BCV), BC_W_g_ramp_t1(BCV),&
                       BC_W_g_ramp_v0(BCV), BC_W_g_ramp_v1(BCV))
!
! z-component of solids phase velocity
      DO M = 1, SMAX
         if(BC_W_s_ramp_t1(BCV,M)/=UNDEFINED)  &
            BC_W_s(BCV,M) = set_bc_value (time, BC_W_s_ramp_t0(BCV,M), BC_W_s_ramp_t1(BCV,M),&
                          BC_W_s_ramp_v0(BCV,M), BC_W_s_ramp_v1(BCV,M))
      ENDDO
!
! Determine which solids phases are present.
      SKIP=(BC_ROP_S(BCV,:)==UNDEFINED.OR.BC_ROP_S(BCV,:)==ZERO) &
         .AND.(BC_EP_S(BCV,:)==UNDEFINED.OR.BC_EP_S(BCV,:)==ZERO)

         IF(MMAX_TOT == 1 .AND. BC_EP_g(BCV)/=ONE) SKIP(1) = .FALSE.

! set bc scalar values to corresponding field values; this this would be
! a problem for bc_types other than mass inflow (i.e., any outflow)
      CALL SET_BC0_INFLOW(BCV)

! set bc velocity values to corresponding field values
      CALL SET_BC0_VEL_INFLOW(BCV)

      RETURN
      END SUBROUTINE SET_BC1_RAMP


      DOUBLE PRECISION FUNCTION set_bc_value (time, bc_t0, bc_t1, bc_v0, bc_v1)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: time, bc_t0, bc_t1, bc_v0, bc_v1

      if(time<bc_t0) then
         set_bc_value = bc_v0
      elseif(time>=bc_t0.and.time<=bc_t1) then
         set_bc_value = bc_v0 + (bc_v1 - bc_v0)/(bc_t1 - bc_t0) * (time - bc_t0)
      else
         set_bc_value = bc_v1
      endif


      RETURN
      END FUNCTION set_bc_value


END MODULE SET_BC1_MOD
