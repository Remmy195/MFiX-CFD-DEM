#include "error.inc"

MODULE CHECK_BOUNDARY_CONDITIONS_MOD

   use bc, only: bc_type_enum, mass_inflow, p_inflow, outflow, mass_outflow, p_outflow, free_slip_wall, no_slip_wall, par_slip_wall, bc_defined, dummy
   use bc, only: is_cg, cg_nsw, cg_fsw, cg_psw
   use check_bc_dem_mod, only:  check_bc_dem
   use check_bc_pic_mod, only:  check_bc_pic
   use check_bc_walls_mod, only:  check_bc_walls
   use check_bc_geometry_mod, only: check_bc_geometry_flow, check_bc_geometry_wall, check_bc_geometry
   use check_bc_inflow_mod, only: check_bc_inflow, check_bc_pressure_inflow, check_bc_mass_inflow
   use check_bc_outflow_mod, only: check_bc_pressure_flow, check_bc_mass_outflowa, check_bc_mass_outflowb
   use cutcell, only: CARTESIAN_GRID
   use error_manager

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BOUNDARY_CONDITIONS                               !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Check boundary condition specifications                    !
!     - convert physical locations to i, j, k's (GET_FLOW_BC)          !
!     - compute area of boundary surfaces (GET_BC_AREA)                !
!     - convert mass and volumetric flows to velocities (FLOW_TO_VEL)  !
!     - check specification of physical quantities                     !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE CHECK_BOUNDARY_CONDITIONS

! Global Variables:
!---------------------------------------------------------------------//
! Total number of (actual) continuum solids.
      use physprop, only: SMAX, MMAX
! Total number of discrete solids.
      use discretelement, only: DES_MMAX
! Special case to handle mixed mi/po for dem
      use discretelement, only: DISCRETE_ELEMENT, CGDEM
      use bc, only: bc_po_apply_to_des, bc_mi_apply_to_des
! Flag: BC dimensions or Type is specified
      use bc, only: BC_DEFINED
! User specified BC solids bulk density
      use bc, only: BC_ROP_s
! Solids volume fraction at BC
      use bc, only: BC_EP_s
      use bc, only: BC_EP_g
! Allow Pressure outflow to inject DEM particles
      use bc, only: bc_mi_apply_to_des
! Run-time flag for DEM solids
      use run, only: DEM_SOLIDS
! Run-time flag for PIC solids
      use run, only: PIC_SOLIDS

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants
      use param1, only: ZERO, ONE, UNDEFINED
! Maximum number of BCs
      use param, only: DIMENSION_BC
! Maximum number of disperse phases
      use param, only: DIM_M

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

! Skip data check when doing preprocessing only
      USE run, only:ppo
      USE run, only: GENERATE_MESH

      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Loop counter for BCs
      INTEGER :: BCV
! Total number of solids phases (continuum + discrete)
      INTEGER :: MMAX_TOT
! Flag to skip checks on indexed solid phase.
      LOGICAL :: SKIP(1:DIM_M)
!......................................................................!

! Determine which BCs are DEFINED
      CALL CHECK_BC_GEOMETRY

! Only check geometry when generating mesh
      IF(GENERATE_MESH) THEN
         DO BCV = 1, DIMENSION_BC

            IF (BC_DEFINED(BCV)) THEN

               SELECT CASE (BC_TYPE_ENUM(BCV))

               CASE (MASS_INFLOW)
                  CALL CHECK_BC_GEOMETRY_FLOW(BCV)

               CASE (P_INFLOW)
                  CALL CHECK_BC_GEOMETRY_FLOW(BCV)

               CASE (OUTFLOW)
                  CALL CHECK_BC_GEOMETRY_FLOW(BCV)

               CASE (MASS_OUTFLOW)
                  CALL CHECK_BC_GEOMETRY_FLOW(BCV)

               CASE (P_OUTFLOW)
                  CALL CHECK_BC_GEOMETRY_FLOW(BCV)

               CASE (FREE_SLIP_WALL)
                  CALL CHECK_BC_GEOMETRY_WALL(BCV)

               CASE (NO_SLIP_WALL)
                  CALL CHECK_BC_GEOMETRY_WALL(BCV)

               CASE (PAR_SLIP_WALL)
                  CALL CHECK_BC_GEOMETRY_WALL(BCV)

               END SELECT
            ENDIF
         ENDDO
         ! Skip remaining BC checks if only doing preprocessing
         IF(PPO) RETURN
      ENDIF

! Total number of solids. (this won't work for GHD/hybrid)
      MMAX_TOT = SMAX + DES_MMAX

! Loop over each defined BC and check the user data.
      DO BCV = 1, DIMENSION_BC

         IF (BC_DEFINED(BCV)) THEN

! Determine which solids phases are present.
            SKIP=(BC_ROP_S(BCV,:)==UNDEFINED.OR.BC_ROP_S(BCV,:)==ZERO) &
               .AND.(BC_EP_S(BCV,:)==UNDEFINED.OR.BC_EP_S(BCV,:)==ZERO)

            IF(MMAX_TOT == 1 .AND. BC_EP_g(BCV)/=ONE) SKIP(1) = .FALSE.

            IF(IS_CG(BC_TYPE_ENUM(BCV)) .AND. .NOT. cartesian_grid) THEN
               WRITE(ERR_MSG,1200)  trim(iVal(BCV))
               CALL LOG_ERROR()
            ENDIF

            SELECT CASE (BC_TYPE_ENUM(BCV))

            CASE (MASS_INFLOW)
               CALL CHECK_BC_GEOMETRY_FLOW(BCV)
               CALL CHECK_BC_RAMP_JET(BCV)

               IF(DISCRETE_ELEMENT.OR.CGDEM) THEN
                  IF(bc_mi_apply_to_des(BCV).AND.bc_po_apply_to_des(BCV)) THEN
                     bc_po_apply_to_des(BCV) = .FALSE.
                     ! JFD: turning off this warning message. This is not needed.
                     ! WRITE(ERR_MSG,1201)  trim(iVal(BCV))
                     ! CALL LOG_WARNING()
                  ENDIF
                  IF(BC_EP_G(BCV)==UNDEFINED.OR.BC_EP_G(BCV)==ONE) bc_mi_apply_to_des(BCV) = .False.
               ENDIF

               CALL CHECK_BC_MASS_INFLOW(MMAX_TOT, SKIP, BCV)
               CALL CHECK_BC_INFLOW(MMAX_TOT,SKIP,BCV)

            CASE (P_INFLOW)
               CALL CHECK_BC_GEOMETRY_FLOW(BCV)
               CALL CHECK_BC_PRESSURE_FLOW(MMAX_TOT, BCV)
               CALL CHECK_BC_PRESSURE_INFLOW(MMAX_TOT, SKIP, BCV)
               CALL CHECK_BC_INFLOW(MMAX_TOT, SKIP, BCV)

            CASE (OUTFLOW)
               CALL CHECK_BC_GEOMETRY_FLOW(BCV)

            CASE (MASS_OUTFLOW)
               CALL CHECK_BC_GEOMETRY_FLOW(BCV)
               CALL CHECK_BC_MASS_OUTFLOWA(MMAX_TOT, BCV)
               CALL CHECK_BC_MASS_OUTFLOWB(MMAX_TOT, SKIP, BCV)

            CASE (P_OUTFLOW)
               IF(DISCRETE_ELEMENT.AND.bc_po_apply_to_des(BCV).AND.bc_mi_apply_to_des(BCV)) THEN
                  bc_mi_apply_to_des(BCV) = .FALSE.
                  ! JFD: turning off this warning message. This is not needed.
                  ! WRITE(ERR_MSG,1202)  trim(iVal(BCV))
                  ! CALL LOG_WARNING()
               ENDIF
               CALL CHECK_BC_GEOMETRY_FLOW(BCV)
               CALL CHECK_BC_PRESSURE_FLOW(MMAX_TOT, BCV)

               IF(DISCRETE_ELEMENT.OR.CGDEM) THEN
                  IF(BC_EP_G(BCV)==UNDEFINED.OR.BC_EP_G(BCV)==ONE) bc_mi_apply_to_des(BCV) = .False.
                  IF(BC_MI_APPLY_TO_DES(BCV)) THEN
                     CALL CHECK_BC_MASS_INFLOW(MMAX_TOT, SKIP, BCV)
                     CALL CHECK_BC_INFLOW(MMAX_TOT,SKIP,BCV)
                  ENDIF
               ENDIF

            CASE (FREE_SLIP_WALL, CG_FSW)
               if(.NOT.IS_CG(BC_TYPE_ENUM(BCV))) CALL CHECK_BC_GEOMETRY_WALL(BCV)
               CALL CHECK_BC_WALLS(MMAX_TOT, SKIP, BCV)

            CASE (NO_SLIP_WALL, CG_NSW)
               if(.NOT.IS_CG(BC_TYPE_ENUM(BCV))) CALL CHECK_BC_GEOMETRY_WALL(BCV)
               CALL CHECK_BC_WALLS(MMAX_TOT, SKIP, BCV)

            CASE (PAR_SLIP_WALL, CG_PSW)
               if(.NOT.IS_CG(BC_TYPE_ENUM(BCV))) CALL CHECK_BC_GEOMETRY_WALL(BCV)
               CALL CHECK_BC_WALLS(MMAX_TOT, SKIP, BCV)

            END SELECT

! Check whether BC values are specified for undefined BC locations
         ELSEIF(BC_TYPE_ENUM(BCV) /= DUMMY .AND.                          &
            .NOT.IS_CG(BC_TYPE_ENUM(BCV))) THEN

            CALL CHECK_BC_RANGE(BCV)

         ENDIF
      ENDDO

      RETURN

 1200 FORMAT('Error 1200: Cartesian grid boundary condition is specified in BC# ',A, &
         /'without turning on Cartesian grid method.')

 1201 FORMAT('Warning 1201: Mass inflow BC# ',A, &
         /'Cannot apply both mi and po to des. Turning off bc_po_apply_to_des.')

 1202 FORMAT('Warning 1202: Pressure or Mass outflow BC# ',A, &
         /'Cannot apply both mi and po to des. Turning off bc_mi_apply_to_des.')
      END SUBROUTINE CHECK_BOUNDARY_CONDITIONS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BC_RANGE                                          !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Verify that data was not given for undefined BC regions.   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_RANGE(BCV)

! Global Variables:
!---------------------------------------------------------------------//
! Gas phase BC variables
      use bc, only: BC_EP_g, BC_T_g, BC_X_g, BC_P_g
      use bc, only: BC_U_g, BC_V_g, BC_W_g
! Solids phase BC variables.
      USE bc, only: BC_EP_s, BC_ROP_s, BC_T_s, BC_X_s
      use bc, only: BC_U_s, BC_V_s, BC_W_s
! Scalar equation BC variables.
      USE bc, only: BC_SCALAR


! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constant for unspecified values.
      use param1, only: UNDEFINED
! Maximum number of disperse phases.
      use param, only: DIM_M
! Maximum number of species gas/solids
      use param, only: DIMENSION_N_G, DIMENSION_N_S
! Maximum number of scalar equations.
      use param, only: DIM_SCALAR


! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------/
! Boundary condition index.
      INTEGER, INTENT(in) :: BCV

! Local Variables:
!---------------------------------------------------------------------//
! Generic loop variables.
      INTEGER :: M, N
!......................................................................!

! Check gas phase variables.
      IF(BC_U_G(BCV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_U_g',BCV))
         CALL LOG_ERROR()
      ENDIF
      IF(BC_V_G(BCV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_V_g',BCV))
         CALL LOG_ERROR()
      ENDIF
      IF (BC_W_G(BCV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_W_g',BCV))
         CALL LOG_ERROR()
      ENDIF
      IF (BC_EP_G(BCV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_EP_g',BCV))
         CALL LOG_ERROR()
      ENDIF
      IF (BC_P_G(BCV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_P_g',BCV))
         CALL LOG_ERROR()
      ENDIF
      IF (BC_T_G(BCV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_T_g',BCV))
         CALL LOG_ERROR()
      ENDIF

      DO N = 1, DIMENSION_N_G
         IF(BC_X_G(BCV,N) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_X_g',BCV,N))
            CALL LOG_ERROR()
         ENDIF
      ENDDO

! Check solids phase variables.
      DO M = 1, DIM_M
         IF(BC_ROP_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_ROP_s',BCV,M))
            CALL LOG_ERROR()
         ENDIF
         IF(BC_EP_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_EP_s',BCV,M))
            CALL LOG_ERROR()
         ENDIF
         IF(BC_U_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_U_s',BCV,M))
            CALL LOG_ERROR()
         ENDIF
         IF(BC_V_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_V_s',BCV,M))
            CALL LOG_ERROR()
         ENDIF

         IF(BC_W_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_W_s',BCV,M))
            CALL LOG_ERROR()
         ENDIF
         IF(BC_T_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_T_s',BCV,M))
            CALL LOG_ERROR()
         ENDIF

         DO N = 1, DIMENSION_N_S
            IF(BC_X_S(BCV,M,N) /= UNDEFINED) THEN
               WRITE(ERR_MSG,1100) trim(iVar('BC_X_s',BCV,M,N))
               CALL LOG_ERROR()
            ENDIF
         ENDDO

      ENDDO

! Check scalar equation variables.
      DO N = 1, DIM_SCALAR
         IF(BC_Scalar(BCV,N) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_Scalar',BCV))
            CALL LOG_ERROR()
         ENDIF
      ENDDO

      RETURN

 1100 FORMAT('Error 1100:',A,' specified for an undefined BC location')

   END SUBROUTINE CHECK_BC_RANGE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BC_RAMP_JET                                       !
!  Author: Jeff Dietiker                              Date: 30-May-25  !
!                                                                      !
!  Purpose: Verify that ramp bc and jet bc are not defined             !
!  simultaneously                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_RAMP_JET(BCV)

! Global Variables:
!---------------------------------------------------------------------//
      USE bc

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constant for unspecified values.
      use param1, only: UNDEFINED
! Maximum number of disperse phases.
      use param, only: DIM_M

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------/
! Boundary condition index.
      INTEGER, INTENT(in) :: BCV

! Local Variables:
!---------------------------------------------------------------------//
! Generic loop variables.
      INTEGER :: M
!......................................................................!

! RAMP BC
      BC_RAMP(BCV) = .FALSE.

      CALL CHECK_BC_RAMP('P_G',BCV, 0, &
                          bc_P_g_ramp_t0(BCV),bc_P_g_ramp_t1(BCV), &
                          bc_P_g_ramp_v0(BCV),bc_P_g_ramp_v1(BCV))

      CALL CHECK_BC_RAMP('T_G',BCV, 0, &
                          bc_T_g_ramp_t0(BCV),bc_T_g_ramp_t1(BCV), &
                          bc_T_g_ramp_v0(BCV),bc_T_g_ramp_v1(BCV))

      CALL CHECK_BC_RAMP('U_G',BCV, 0, &
                          bc_U_g_ramp_t0(BCV),bc_U_g_ramp_t1(BCV), &
                          bc_U_g_ramp_v0(BCV),bc_U_g_ramp_v1(BCV))

      CALL CHECK_BC_RAMP('V_G',BCV, 0, &
                          bc_V_g_ramp_t0(BCV),bc_V_g_ramp_t1(BCV), &
                          bc_V_g_ramp_v0(BCV),bc_V_g_ramp_v1(BCV))

      CALL CHECK_BC_RAMP('W_G',BCV, 0, &
                          bc_W_g_ramp_t0(BCV),bc_W_g_ramp_t1(BCV), &
                          bc_W_g_ramp_v0(BCV),bc_W_g_ramp_v1(BCV))

      CALL CHECK_BC_RAMP('VOLFLOW_G',BCV, 0, &
                          bc_VOLFLOW_g_ramp_t0(BCV),bc_VOLFLOW_g_ramp_t1(BCV), &
                          bc_VOLFLOW_g_ramp_v0(BCV),bc_VOLFLOW_g_ramp_v1(BCV))

      CALL CHECK_BC_RAMP('MASSFLOW_G',BCV, 0, &
                          bc_MASSFLOW_g_ramp_t0(BCV),bc_MASSFLOW_g_ramp_t1(BCV), &
                          bc_MASSFLOW_g_ramp_v0(BCV),bc_MASSFLOW_g_ramp_v1(BCV))

      DO M = 1, DIM_M
         CALL CHECK_BC_RAMP('T_S',BCV, M, &
                             bc_T_s_ramp_t0(BCV, M),bc_T_s_ramp_t1(BCV, M), &
                             bc_T_s_ramp_v0(BCV, M),bc_T_s_ramp_v1(BCV, M))

         CALL CHECK_BC_RAMP('THETA_M',BCV, M, &
                             bc_Theta_m_ramp_t0(BCV, M),bc_Theta_m_ramp_t1(BCV, M), &
                             bc_Theta_m_ramp_v0(BCV, M),bc_Theta_m_ramp_v1(BCV, M))

         CALL CHECK_BC_RAMP('U_S',BCV, M, &
                             bc_U_s_ramp_t0(BCV, M),bc_U_s_ramp_t1(BCV, M), &
                             bc_U_s_ramp_v0(BCV, M),bc_U_s_ramp_v1(BCV, M))

         CALL CHECK_BC_RAMP('V_S',BCV, M, &
                             bc_V_s_ramp_t0(BCV, M),bc_V_s_ramp_t1(BCV, M), &
                             bc_V_s_ramp_v0(BCV, M),bc_V_s_ramp_v1(BCV, M))

         CALL CHECK_BC_RAMP('W_S',BCV, M, &
                             bc_W_s_ramp_t0(BCV, M),bc_W_s_ramp_t1(BCV, M), &
                             bc_W_s_ramp_v0(BCV, M),bc_W_s_ramp_v1(BCV, M))

         CALL CHECK_BC_RAMP('VOLFLOW_S',BCV, M, &
                             bc_VOLFLOW_s_ramp_t0(BCV, M),bc_VOLFLOW_s_ramp_t1(BCV, M), &
                             bc_VOLFLOW_s_ramp_v0(BCV, M),bc_VOLFLOW_s_ramp_v1(BCV, M))

         CALL CHECK_BC_RAMP('MASSFLOW_S',BCV, M, &
                             bc_MASSFLOW_s_ramp_t0(BCV, M),bc_MASSFLOW_s_ramp_t1(BCV, M), &
                             bc_MASSFLOW_s_ramp_v0(BCV, M),bc_MASSFLOW_s_ramp_v1(BCV, M))
      ENDDO


! JET BC
      BC_JET(BCV) = .FALSE.

      CALL CHECK_BC_JET(BCV)

! Check that only BC_RAMP OR BC_JET is defined. Having both BC_RAMP and
! BC_JET defined simultaneously is not supported.

      IF(BC_RAMP(BCV).AND.BC_JET(BCV)) THEN
         WRITE(ERR_MSG,1100) BCV
         CALL LOG_ERROR()
      ENDIF

      RETURN

 1100 FORMAT('Error 1200: Both ramp and jet BCs are defined for BC#:',I4, &
            /'This is not supported.')

   END SUBROUTINE CHECK_BC_RAMP_JET

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BC_RAMP                                           !
!  Author: Jeff Dietiker                              Date: 30-May-25  !
!                                                                      !
!  Purpose: Verify that ramp bc parameters are all defined.            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE CHECK_BC_RAMP(key, BCV, M, t0, t1, v0, v1)

! Global Variables:
!---------------------------------------------------------------------//
      USE bc

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constant for unspecified values.
      use param1, only: UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------/
      CHARACTER(LEN=*), INTENT(in) :: key
      INTEGER, INTENT(in) :: BCV, M
      DOUBLE PRECISION, INTENT(in):: t0,t1,v0,v1

! Local Variables:
!---------------------------------------------------------------------//
      CHARACTER(LEN=32) :: keychar, kt0, kt1, kv0, kv1
      LOGICAL :: BCRAMP
!......................................................................!

      BCRAMP = .FALSE.

      IF(t1/=UNDEFINED) BCRAMP = .TRUE.
      IF(v0/=UNDEFINED) BCRAMP = .TRUE.
      IF(v1/=UNDEFINED) BCRAMP = .TRUE.

      IF(BCRAMP) THEN
         IF(t1==UNDEFINED.OR.v0==UNDEFINED.OR.v1==UNDEFINED) THEN
            keychar = 'BC_'//KEY//'_RAMP'
            kt0 = trim(keychar)//'_T0'
            kt1 = trim(keychar)//'_T1'
            kv0 = trim(keychar)//'_V0'
            kv1 = trim(keychar)//'_V1'
            WRITE(ERR_MSG,1100) trim(keychar), &
                                BCV, M, &
                                trim(kt0),t0, &
                                trim(kt1),t1, &
                                trim(kv0),v0, &
                                trim(kv1),v1
            CALL LOG_ERROR()
         ENDIF

         BC_RAMP(BCV) = .TRUE.

      ENDIF

 1100 FORMAT('Error 1200: All bc ramp parameters must be defined for : ',A, &
            /'BC#   =',I4, &
            /'Phase =',I4, &
            4(/A,' = ',G0))

      RETURN
   END SUBROUTINE CHECK_BC_RAMP
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BC_JET                                            !
!  Author: Jeff Dietiker                              Date: 30-May-25  !
!                                                                      !
!  Purpose: Verify that jet bc parameters are all defined.             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE CHECK_BC_JET(BCV)

! Global Variables:
!---------------------------------------------------------------------//
      USE bc

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constant for unspecified values.
      use param1, only: UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------/
      INTEGER, INTENT(in) :: BCV

! Local Variables:
!---------------------------------------------------------------------//
      LOGICAL :: BCJET
!......................................................................!

      BCJET = .FALSE.

      IF(BC_JET_G0(BCV)/=UNDEFINED) BCJET = .TRUE.
      IF(BC_JET_GH(BCV)/=UNDEFINED) BCJET = .TRUE.
      IF(BC_JET_GL(BCV)/=UNDEFINED) BCJET = .TRUE.
      IF(BC_DT_0(BCV)/=UNDEFINED) BCJET = .TRUE.
      IF(BC_DT_H(BCV)/=UNDEFINED) BCJET = .TRUE.
      IF(BC_DT_L(BCV)/=UNDEFINED) BCJET = .TRUE.

      IF(BCJET) THEN
         IF(BC_JET_G0(BCV)==UNDEFINED.OR. &
            BC_JET_GH(BCV)==UNDEFINED.OR. &
            BC_JET_GL(BCV)==UNDEFINED.OR. &
            BC_DT_0(BCV)==UNDEFINED.OR.   &
            BC_DT_H(BCV)==UNDEFINED.OR.   &
            BC_DT_L(BCV)==UNDEFINED) THEN

            WRITE(ERR_MSG,1100) BCV, &
                                trim(iVar('BC_JET_GO',BCV)), BC_JET_G0(BCV),  &
                                trim(iVar('BC_JET_GH',BCV)), BC_JET_GH(BCV),  &
                                trim(iVar('BC_JET_GL',BCV)), BC_JET_GL(BCV),  &
                                trim(iVar('BC_DT_O',BCV)), BC_DT_0(BCV),  &
                                trim(iVar('BC_DT_H',BCV)), BC_DT_H(BCV),  &
                                trim(iVar('BC_DT_L',BCV)), BC_DT_L(BCV)
            CALL LOG_ERROR()
         ENDIF

         BC_JET(BCV) = .TRUE.

      ENDIF

 1100 FORMAT('Error 1210: All bc jet parameters must be defined for BC#: ',I4, &
            6(/A,' = ',G0))

      RETURN
   END SUBROUTINE CHECK_BC_JET

END MODULE CHECK_BOUNDARY_CONDITIONS_MOD
