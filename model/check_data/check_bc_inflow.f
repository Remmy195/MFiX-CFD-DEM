#include "error.inc"

MODULE CHECK_BC_INFLOW_MOD

   use error_manager

CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_PRESSURE_INFLOW                                 !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided detailed error message on bc                       !
!                                                                      !
! Comments:                                                            !
!     Unlike the MI boundary, for the PI boundary the velocities at    !
!     the inflow face are calculated by solving the momentum eqns      !
!     and are not fixed. In this way, the PI is similar to the PO      !
!     except that the flow is into the domain and hence all other      !
!     scalars (e.g., mass fractions, void fraction, temperature,       !
!     etc.,) at the inflow cells need to be specified. To satisfy      !
!     the error routines at the start of the simulation, both the      !
!     tangential and normal components at the inflow also need to      !
!     be specified. The velocities values essentially serve as IC.     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_PRESSURE_INFLOW(M_TOT, SKIP, BCV)

! Modules
!---------------------------------------------------------------------//
      use bc, only: bc_p_g, bc_rop_s
      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use bc, only: bc_u_s, bc_v_s, bc_w_s
      use geometry, only: no_i, no_j, no_k
      use param, only: dim_m
      use param1, only: undefined, zero
      use physprop, only: ro_g0

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: M_TOT
      LOGICAL, INTENT(in) :: SKIP(DIM_M)
      INTEGER, INTENT(in) :: BCV

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: M

!---------------------------------------------------------------------//

! Check that velocities are also specified. These are essentially used
! as initial conditions for the boundary region. If they are not
! specified then a default value is set here otherwise check_data_20
! will complain and cause MFIX to exit.
      IF(BC_U_G(BCV) == UNDEFINED) THEN
         BC_U_G(BCV) = ZERO
         IF(.NOT.NO_I) THEN
            WRITE(ERR_MSG, 1300) trim(iVar('BC_U_g',BCV))
            CALL LOG_WARNING()
         ENDIF
      ENDIF

      IF(BC_V_G(BCV) == UNDEFINED) THEN
         BC_V_G(BCV) = ZERO
         IF(.NOT.NO_J) THEN
            WRITE(ERR_MSG, 1300) trim(iVar('BC_V_g',BCV))
            CALL LOG_WARNING()
         ENDIF
      ENDIF

      IF(BC_W_G(BCV) == UNDEFINED) THEN
         BC_W_G(BCV) = ZERO
         IF(.NOT.NO_K) THEN
            WRITE(ERR_MSG, 1300) trim(iVar('BC_W_g',BCV))
            CALL LOG_WARNING()
         ENDIF
      ENDIF

      DO M = 1, M_TOT
         IF (SKIP(M)) THEN
            BC_U_S(BCV,M) = ZERO
            BC_V_S(BCV,M) = ZERO
            BC_W_S(BCV,M) = ZERO
         ELSE
            IF(BC_U_S(BCV,M) == UNDEFINED) THEN
               BC_U_S(BCV,M) = ZERO
               IF(BC_ROP_S(BCV,M) /= ZERO .AND. .NOT.NO_I) THEN
                  WRITE(ERR_MSG, 1300) trim(iVar('BC_U_s',BCV,M))
                  CALL LOG_WARNING()
               ENDIF
            ENDIF

            IF(BC_V_S(BCV,M) == UNDEFINED) THEN
               BC_V_S(BCV,M) = ZERO
               IF(BC_ROP_S(BCV,M) /= ZERO .AND. .NOT.NO_J) THEN
                  WRITE(ERR_MSG, 1300) trim(iVar('BC_V_s',BCV,M))
                  CALL LOG_WARNING()
               ENDIF
            ENDIF

            IF(BC_W_S(BCV,M) == UNDEFINED) THEN
               BC_W_S(BCV,M) = ZERO
               IF(BC_ROP_S(BCV,M) /= ZERO .AND. .NOT.NO_K) THEN
                  WRITE(ERR_MSG, 1300) trim(iVar('BC_W_s',BCV,M))
                  CALL LOG_WARNING()
               ENDIF
            ENDIF
         ENDIF
      ENDDO

 1300 FORMAT('Warning 1300: ',A,' was undefined. This variable was ', &
         'set ',/ 'to zero to be used as the initial value in the BC ',&
         'region.')

      RETURN

   END SUBROUTINE CHECK_BC_PRESSURE_INFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_INFLOW                                          !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message on common inflow BC       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_INFLOW(M_TOT, SKIP, BCV)

! Modules
!---------------------------------------------------------------------
      use bc, only: bc_type_enum, mass_inflow
      use bc, only: bc_x_g, bc_x_s
      use bc, only: bc_t_g, bc_t_s, bc_theta_m
      use bc, only: bc_k_turb_g, bc_e_turb_g
      use bc, only: bc_scalar
      use bc, only: bc_u_s, bc_v_s, bc_w_s
      use bc, only: bc_t_g_ramp_v0, bc_t_s_ramp_v0
      use param, only: dim_m
      use param1, only: undefined, one, zero
      use physprop, only: inert_species
      use physprop, only: mu_g0
      use physprop, only: mw_avg, nmax
      use physprop, only: ro_g0
      use run, only: energy_eq, granular_energy
      use run, only: solids_model, solve_ros, species_eq
      use scalars, only: nscalar
      use toleranc, only: compare
      use turb, only: k_epsilon

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: BCV
      INTEGER, INTENT(in) :: M_TOT
      LOGICAL, INTENT(in) :: SKIP(DIM_M)

! Local variables
!---------------------------------------------------------------------
! loop/variable indices
      INTEGER :: M, N
      DOUBLE PRECISION SUM_X
! Index of inert species
      INTEGER :: INERT
!---------------------------------------------------------------------

! Check temperature dependency.
      IF(BC_T_G_RAMP_V0(BCV)/=UNDEFINED.AND.BC_T_G(BCV)==UNDEFINED) &
         BC_T_G(BCV) = BC_T_G_RAMP_V0(BCV)

      IF((ENERGY_EQ .OR. RO_G0==UNDEFINED .OR.MU_G0==UNDEFINED) .AND. &
         BC_T_G(BCV)==UNDEFINED) THEN
         WRITE(ERR_MSG, 1000) trim(iVar('BC_T_g',BCV))
         CALL LOG_ERROR()
      ENDIF

! Sum together defined gas phase species mass fractions.
      SUM_X = ZERO
      DO N = 1, NMAX(0)
         IF(BC_X_G(BCV,N) /= UNDEFINED) THEN
            SUM_X = SUM_X + BC_X_G(BCV,N)
         ELSE
            BC_X_G(BCV,N) = ZERO
         ENDIF
      ENDDO

! Enforce that the species mass fractions must sum to one.
      IF(.NOT.COMPARE(ONE,SUM_X)) THEN

         IF(SPECIES_EQ(0)) THEN
            WRITE(ERR_MSG, 1110) BCV
            CALL LOG_ERROR()

 1110 FORMAT('Error 1110: BC_X_g(',I3,',:) do not sum to 1 and the ',&
        'gas phase',/'species equations are solved.')

         ELSEIF(RO_G0 == UNDEFINED .AND. MW_AVG == UNDEFINED) THEN
            WRITE(ERR_MSG, 1111) BCV
            CALL LOG_ERROR()

 1111 FORMAT('Error 1111: BC_X_g(',I3,',:) do not sum to 1 and the ',&
         'gas phase',/'is compressible and MW_AVG is UNDEFINED.')

         ELSEIF(.NOT.COMPARE(SUM_X,ZERO)) THEN
            WRITE(ERR_MSG, 1112) BCV
            CALL LOG_ERROR()

 1112 FORMAT('Error 1112: BC_X_g(',I3,',:) do not sum to 1 or 0 ',&
         'and they',/'are not needed.')

         ELSE
            BC_X_G(BCV,:) = ZERO
            BC_X_G(BCV,1) = ONE
         ENDIF
      ENDIF


! Verify that species mass fractions are defined for mass flow BCs when
! using variable solids density. Needed to calculation RO_s
      DO M = 1, M_TOT

! If this phase is not present, clear out x_s for the BC and
! cycle the solids loop. No need to continue checks.
         IF(SKIP(M)) THEN
            IF(SPECIES_EQ(M))THEN
               BC_X_S(BCV,M,:) = ZERO
               BC_X_S(BCV,M,1) = ONE
            ENDIF
            CYCLE
         ENDIF
! Sum together defined species mass fractions.
         SUM_X = ZERO
         DO N = 1, NMAX(M)
            IF(BC_X_S(BCV,M,N) /= UNDEFINED) THEN
               SUM_X = SUM_X + BC_X_S(BCV,M,N)
            ELSE
               BC_X_S(BCV,M,N) = ZERO
            ENDIF
         ENDDO

! Enforce that the species mass fractions must sum to one.
         IF(.NOT.COMPARE(ONE,SUM_X)) THEN
            IF(SPECIES_EQ(M)) THEN
               WRITE(ERR_MSG, 1210) BCV, M
               CALL LOG_ERROR()

 1210 FORMAT('Error 1210: BC_X_s(',I3,',',I2,':) do not sum to 1 ',  &
         'and the solids phase',/'species equations are solved. ')

            ELSEIF(SOLVE_ROS(M)) THEN
               WRITE(ERR_MSG, 1211) BCV, M
               CALL LOG_ERROR()

 1211 FORMAT('Error 1211: BC_X_s(',I3,',',I2,':) do not sum to 1 ',  &
         'and the solids phase',/'density is calculated.')

            ELSEIF(.NOT.COMPARE(SUM_X,ZERO)) THEN
               WRITE(ERR_MSG, 1212) BCV, M
               CALL LOG_ERROR()

 1212 FORMAT('Error 1212: BC_X_s(',I3,',',I2,':) do not sum to 1 ',  &
         'or ZERO and',/'they are not needed.')

            ELSE
               BC_X_S(BCV,M,:) = ZERO
               BC_X_S(BCV,M,1) = ONE
            ENDIF
         ENDIF

! Set the solids density for the BC region.
         IF(SOLVE_ROs(M)) THEN
! Verify that the species mass fraction for the inert material is not
! zero in the IC region when the solids is present.
            INERT = INERT_SPECIES(M)
            IF(BC_X_S(BCV,M,INERT) == ZERO) THEN
               WRITE(ERR_MSG,1213) M, BCV
               CALL LOG_ERROR()
 1213 FORMAT('Error 1213: No inert species for phase ',I2,' in BC ',   &
         'region',I3,'.',/'Unable to calculate solids phase density.')
            ENDIF
         ENDIF
      ENDDO


      DO M = 1, M_TOT
! Check solids phase temperature dependency.
         IF(BC_T_S_RAMP_V0(BCV,M)/=UNDEFINED.AND.BC_T_S(BCV,M)==UNDEFINED) &
            BC_T_S(BCV,M) = BC_T_S_RAMP_V0(BCV,M)

         IF(ENERGY_EQ .AND. BC_T_S(BCV,M)==UNDEFINED) THEN
            IF(SKIP(M)) THEN
               BC_T_S(BCV,M) = BC_T_G(BCV)
            ELSE
               WRITE(ERR_MSG, 1000) trim(iVar('BC_T_s',BCV,M))
               CALL LOG_ERROR()
            ENDIF
         ENDIF

! Check granular energy dependency
         IF(GRANULAR_ENERGY) THEN
            IF(BC_THETA_M(BCV,M) == UNDEFINED) THEN
               IF(SKIP(M) .OR. SOLIDS_MODEL(M) /= 'TFM') THEN
                  BC_THETA_M(BCV,M) = ZERO
               ELSE
                  WRITE(ERR_MSG,1000) trim(iVar('BC_Theta_m',BCV,M))
                  CALL LOG_ERROR()
               ENDIF
            ENDIF
         ENDIF
      ENDDO

! Check K-Epsilon BCs.
      IF(K_Epsilon) THEN
         IF(BC_K_Turb_G(BCV) == UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) trim(iVar('BC_K_Turb_G',BCV))
            CALL LOG_ERROR()
         ENDIF
         IF (BC_K_Turb_G(BCV) < ZERO) THEN
            WRITE(ERR_MSG, 1002) iVar('BC_K_Turb_G',BCV),ival(BC_K_Turb_G(BCV))
            CALL LOG_ERROR()
         ENDIF

         IF(BC_E_Turb_G(BCV) == UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) trim(iVar('BC_E_Turb_G',BCV))
            CALL LOG_ERROR()
         ENDIF
         IF (BC_E_Turb_G(BCV) <= ZERO) THEN
            WRITE(ERR_MSG, 1003) iVar('BC_E_Turb_G',BCV),ival(BC_E_Turb_G(BCV))
            CALL LOG_ERROR()
         ENDIF
      ENDIF

! Check scalar equation BCs.
      DO N = 1, NScalar
         IF(BC_Scalar(BCV,N) == UNDEFINED) THEN
            WRITE(ERR_MSG, 1001) trim(iVar('BC_Scalar',BCV,N))
            CALL LOG_ERROR()
         ENDIF
      ENDDO

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A)

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A)

 1002 FORMAT('Error 1002: Turbulent kinetic energy must be greater ',&
        'than or equal to 0.0 ', /'Illegal value: ',A,' = ',A)

 1003 FORMAT('Error 1003: Turbulent dissipation rate must be ',&
        'greater than 0.0 ',     &
         /'Illegal value: ',A,' = ',A)
      END SUBROUTINE CHECK_BC_INFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_MASS_INFLOW                                     !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message on BC                     !
!                                                                      !
! Comments:                                                            !
!     The velocities at the inflow/outflow face are fixed; the         !
!     momentum equations are not solved in the MI/MO cells. Depending  !
!     on user specification (mass flow or volume flow or velocity)     !
!     different checks need to be made. Consistently on the sum of     !
!     the volume fractions at the boundary are also made for MI.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_MASS_INFLOW(M_TOT, SKIP, BCV)

! Modules
!---------------------------------------------------------------------
      use bc, only: bc_ep_g, bc_p_g
      use bc, only: bc_rop_s, bc_ep_s, bc_x_s
      use bc, only: bc_massflow_g, bc_massflow_s
      use bc, only: bc_volflow_g, bc_volflow_s
      use bc, only: bc_p_g_ramp_v0
      use eos, only: eoss
      use param, only: dim_m, dimension_bc
      use param1, only: undefined, one, zero
      use physprop, only: ro_g0, ro_s0
      use physprop, only: inert_species, x_s0
      use run, only: solve_ros
      use toleranc, only: compare
      use usr_prop, only: usr_ros, usr_rog

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: M_TOT
      LOGICAL, INTENT(in) :: SKIP(DIM_M)
      INTEGER, INTENT(in) :: BCV

! Local variables
!---------------------------------------------------------------------
! loop/variable indices
      INTEGER :: M
      DOUBLE PRECISION :: SUM_EP
! Solids phase density in BC region.
      DOUBLE PRECISION :: BC_ROs(DIM_M)
! Index of inert species
      INTEGER :: INERT
! check routine looks that user has not specified both rop_s and ep_s
! (see routine below for details)
      LOGICAL, SAVE :: FIRST_PASS(DIMENSION_BC) = .TRUE.
!---------------------------------------------------------------------
! The check on velocities, among other checks, is made slightly later
! in set_bc0_flow -> check_bc_vel_inflow. Massflow or volflow needs
! to be converted to velocity using known density/volume fraction.
! So such checks are delayed until the calculations can be done.

! Check gas phase volume fraction.
! BC_EP_G is needed not only to convert massflow and volflow to a
! velocity but also in defining the flux of material through the
! boundary plane (rop_g*v_g)
      IF(BC_EP_G(BCV) == UNDEFINED) THEN
         WRITE(ERR_MSG, 1096) trim(iVar('BC_EP_G',BCV))
         CALL LOG_ERROR()
 1096 FORMAT('Error 1096: ',A,' must be defined when requesting',/, &
         'a MASS_INFLOW boundary.')
      ENDIF

! Verify compressible boundary condition variables.
      IF(RO_G0 == UNDEFINED) THEN
         IF(BC_P_G_RAMP_V0(BCV)/=UNDEFINED.AND.BC_P_G(BCV)==UNDEFINED) &
            BC_P_G(BCV) = BC_P_G_RAMP_V0(BCV)
         IF(BC_P_G(BCV) == UNDEFINED) THEN
            IF(BC_MASSFLOW_G(BCV) /= UNDEFINED .AND.                   &
               BC_MASSFLOW_G(BCV) /= ZERO) THEN
               WRITE(ERR_MSG, 1100) trim(iVar('BC_P_g',BCV))
               CALL LOG_ERROR()
            ENDIF
 1100 FORMAT('Error 1100: ',A,' must be specified for compressible ',  &
         'flows',/'when specifying BC_MASSFLOW_g to make the ',        &
         'conversion to velocity.')


         ELSEIF(BC_P_G(BCV) <= ZERO) THEN
            WRITE(ERR_MSG, 1101) BCV, trim(iVal(BC_P_G(BCV)))
            CALL LOG_ERROR()
         ENDIF
 1101 FORMAT('Error 1101: Pressure must be greater than zero for ',    &
         'compressible flow',/' >>>  BC_P_g(',I3,') = ',A)
      ENDIF

! Calculate the solids volume fraction from the gas phase if there is
! only one solids phase. Do this check regardless of whether it needs
! to be set or not.
      IF(M_TOT == 1 .AND. BC_EP_S(BCV,1) == UNDEFINED) THEN
         IF (BC_EP_G(BCV) /= UNDEFINED) BC_EP_S(BCV,1) = ONE - BC_EP_g(BCV)
      ENDIF

! Define both as zero since ep_g=1 (overwrite any other user setting)
      IF(COMPARE(BC_EP_g(BCV),ONE)) THEN
         DO M = 1, M_TOT
            BC_EP_S(BCV,M)  = ZERO
            BC_ROP_S(BCV,M) = ZERO
         ENDDO
      ENDIF

! At this point, either one or both ep_s/rop_s must be defined as zero;
! So make sure both are set to zero.
      DO M = 1, M_TOT
         IF(SKIP(M)) THEN
            BC_EP_S(BCV,M)  = ZERO
            BC_ROP_S(BCV,M) = ZERO
            CYCLE
         ENDIF
      ENDDO

! Bulk density or solids volume fraction should be explicitly defined
! but not both.
      DO M = 1, M_TOT
         IF(BC_ROP_S(BCV,M) == UNDEFINED .AND. &
            BC_EP_S(BCV,M) == UNDEFINED) THEN
            WRITE(ERR_MSG, 1205) BCV, M, 'BC_ROP_s or BC_EP_s'
            CALL LOG_ERROR()
 1205 FORMAT('Error 1205: Incomplete specification for BC', I3,'.',/&
         'BC_EP_G is defined but insufficient corresponding solids ',&
         'phase ',I2,' data:',/,A,' not specified.')
         ENDIF

         IF (FIRST_PASS(BCV)) THEN
! JEC: currently a necessary work-around due to a very generic second
! call to check_boundary_conditions following get_stl_data rather than
! a targeted call invoking only which is needed. Ideally, the second
! call would be more targeted or a top-level skip would be enforced for
! bc_types
         IF((BC_ROP_S(BCV,M) /= UNDEFINED .AND. &
            BC_ROP_S(BCV,M) /= ZERO) .AND. &
            (BC_EP_S(BCV,M) /= UNDEFINED .AND. &
            BC_EP_S(BCV,m) /= ZERO)) THEN
            WRITE(ERR_MSG, 1206) BCV, M
            CALL LOG_ERROR()
 1206 FORMAT('Error 1206: Over specification for BC', I3,'.',/&
         'Define either BC_ROP_S or BC_EP_s for phase ', I2,' but ',&
         ' not both ',/,'to ensure consistency.')
         ENDIF
         FIRST_PASS(BCV) = .FALSE.
         ENDIF

         IF (USR_ROS(M)) THEN
            IF(BC_ROP_S(BCV,M) /= UNDEFINED .AND. &
               BC_ROP_S(BCV,M) /= ZERO) THEN
               WRITE(ERR_MSG, 1207) trim(iVar('BC_ROP_S',BCV,M))
               CALL LOG_ERROR()
 1207 FORMAT('Error 1207: Specification of ',A,' not allowed ',&
         'in combination ',/,'with usr_ros. Specify BC_EP_S ',&
         'instead.')
            ENDIF
            IF(BC_EP_S(BCV,M) == UNDEFINED) THEN
               WRITE(ERR_MSG, 1208) trim(iVar('BC_EP_S',BCV,M))
               CALL LOG_ERROR()
 1208 FORMAT('Error 1208: ', A,' must be defined when requesting',/, &
          'a MASS_INFLOW boundary with usr_ros.')
            ENDIF
         ENDIF
      ENDDO

! For usr defined density bc_ep_s must be specified. This permits the
! check on volume fractions here. Otherwise it would have to wait until
! after set_ro_s
      DO M = 1, M_TOT

! Set the solids density for the BC region.
         IF(SOLVE_ROs(M)) THEN
! presence of non-zero inert species is checked by bc_inflow
            INERT = INERT_SPECIES(M)
            BC_ROs(M) = EOSS(RO_s0(M), X_s0(M,INERT),&
               BC_X_S(BCV,M,INERT))
         ELSEIF (USR_ROs(M)) THEN
! Set an BC density for this material. The actual value is not overly
! important. This facilitates the call to set_ro_s wherein the phase
! density is called in all cells. At that point, the value of rop_s
! will be updated based on ro_s and enforcing any ic_ep_s and bc_ep_s.
            BC_ROs(M) = ONE

         ELSE
            BC_ROs(M) = RO_s0(M)
         ENDIF
      ENDDO

! Now check the sum of the total volume fraction; ensure any user
! settings of volume/void fraction are self-consistent
      SUM_EP = BC_EP_G(BCV)
      DO M=1, M_TOT
         IF(BC_ROP_S(BCV,M) == UNDEFINED) THEN
            BC_ROP_S(BCV,M) = BC_EP_S(BCV,M) *BC_ROs(M)
         ELSEIF(BC_EP_S(BCV,M) == UNDEFINED) THEN
            BC_EP_S(BCV,M) = BC_ROP_S(BCV,M) / BC_ROs(M)
         ENDIF
! Add this phase to the total volume fraction.
         SUM_EP = SUM_EP + BC_EP_S(BCV,M)
      ENDDO

! Verify that the volume fractions sum to one.
      IF(.NOT.COMPARE(SUM_EP,ONE)) THEN
         WRITE(ERR_MSG,1215) BCV, trim(iVal(SUM_EP))
         CALL LOG_ERROR()
      ENDIF

 1215 FORMAT('Error 1215: Illegal boundary condition region: ',I3,'. ',&
         'Sum of volume',/'fractions does NOT equal ONE. (SUM = ',A,   &
         ')')

      RETURN
      END SUBROUTINE CHECK_BC_MASS_INFLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_USR_RO                                          !
!                                                                      !
! Purpose: Provided a detailed error message on BC                     !
!                                                                      !
! Comments:                                                            !
!   The following restrictions are currently difficult to lift. To do  !
!   so requires knowing ro in the boundary which is based on user      !
!   defined function.  Therefore all field (i.e. ijk) quantities must  !
!   be specified in the boundary cells in order to evaluate ro.  Also, !
!   the subroutine flow_to_vel would need to be modified to rely on    !
!   such ijk field values rather than singular bc values.              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_USR_RO(M_TOT, BCV)

! Modules
!---------------------------------------------------------------------
      use bc, only: bc_rop_s
      use bc, only: bc_massflow_g, bc_massflow_s
      use param, only: dim_m
      use param1, only: undefined, zero
      use usr_prop, only: usr_ros, usr_rog

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: BCV
      INTEGER, INTENT(in) :: M_TOT

! Local variables
!---------------------------------------------------------------------
! loop/variable indices
      INTEGER :: M

!---------------------------------------------------------------------


      IF(USR_ROG) THEN
         IF(BC_MASSFLOW_G(BCV) /= UNDEFINED .AND.                   &
            BC_MASSFLOW_G(BCV) /= ZERO) THEN
            WRITE(ERR_MSG, 1098) trim(iVar('BC_MASSFLOW_G',BCV))
            CALL LOG_ERROR()
         ENDIF
 1098 FORMAT('Error 1098: ',A,' cannot be used when invoking',/, &
         'a user specified density function via USR_ROg.',&
         'The conversion',/,'to volumetric flow rate is not readily',&
         'possible.')
      ENDIF

      DO M = 1, M_TOT
         IF(USR_ROS(M)) THEN
            IF(BC_MASSFLOW_S(BCV,M) /= UNDEFINED .AND.  &
               BC_MASSFLOW_S(BCV,M) /= ZERO) THEN
               WRITE(ERR_MSG, 1199) trim(iVar('BC_MASSFLOW_S',BCV,M)), &
                  trim(iVAR('USR_ROS',M))
               CALL LOG_ERROR()
            ENDIF
 1199 FORMAT('Error 1199: ',A,' cannot be used when invoking a',/, &
         'user specified density function via', A,&
         '. The conversion',/,'to volumetric flow rate is not readily',&
         'possible.')

            IF(BC_ROP_s(BCV,M) /= UNDEFINED .AND. &
               BC_ROP_S(BCV,M) /= ZERO) THEN
! allow for zero specification since some checks may result in setting
! bc_rop_s to zero (only if bc_ep_s is also zero)
               WRITE(ERR_MSG,1200) trim(iVar('BC_ROP_s',BCV,M)), &
                  trim(iVar('USR_ROs',M))
               CALL LOG_ERROR()
            ENDIF
 1200 FORMAT('Error 1200: ',A,' must not be specified when also ',&
         'requesting ',/,A,' for the same phase.')
         ENDIF
     ENDDO
     RETURN
     END SUBROUTINE CHECK_BC_USR_RO



END MODULE CHECK_BC_INFLOW_MOD
