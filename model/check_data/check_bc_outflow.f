#include "error.inc"

MODULE CHECK_BC_OUTFLOW_MOD

   use error_manager

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_PRESSURE_FLOW                                   !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message on bc                     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_PRESSURE_FLOW(M_TOT, BCV)

! Modules
! --------------------------------------------------------------------
      USE param1, only: UNDEFINED
      USE param1, only: ZERO
      use physprop, only: RO_g0
      use bc, only: BC_P_g

      IMPLICIT NONE

! Dummy arguments
! --------------------------------------------------------------------
! loop/variable indices
      INTEGER, INTENT(in) :: BCV
      INTEGER, INTENT(in) :: M_TOT
! --------------------------------------------------------------------

      IF (BC_P_G(BCV) == UNDEFINED.AND.RO_G0>ZERO) THEN
         WRITE(ERR_MSG,1000) trim(iVar('BC_P_g',BCV))
         CALL LOG_ERROR()

      ELSEIF (BC_P_G(BCV)<=ZERO .AND. RO_G0==UNDEFINED) THEN
         WRITE(ERR_MSG, 1100) BCV, trim(iVal(BC_P_G(BCV)))
         CALL LOG_ERROR()
      ENDIF

 1100 FORMAT('Error 1100: Pressure must be greater than zero for ',    &
         'compressible flow',/3x,'BC_P_g(',I3,') = ',A)

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A)

      END SUBROUTINE CHECK_BC_PRESSURE_FLOW


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_MASS_OUTFLOWA                                   !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Comments:                                                            !
!     The velocities at the outflow face are fixed and the momentum    !
!     equations are not solved in the outflow cells. Since the flow    !
!     is out of the domain none of the other scalars should need to    !
!     be specified (e.g., mass fractions, void fraction, etc.,).       !
!     Such values will become defined according to their adjacent      !
!     fluid cell. Checks are made on massflow/volflow/velocity in      !
!     a separate routine.                                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_MASS_OUTFLOWA(M_TOT, BCV)

! Modules
! --------------------------------------------------------------------

      use bc, only: bc_plane
      use bc, only: bc_dt_0
      use bc, only: bc_massflow_g, bc_volflow_g
      use bc, only: bc_massflow_s, bc_volflow_s
      use bc, only: bc_ep_g, bc_rop_s, bc_ep_s

      use bc, only: bc_u_g, bc_v_g, bc_w_g
      use param, only: dimension_bc
      use param1, only: undefined, zero, one
      use toleranc, only: compare
      IMPLICIT NONE

! Dummy arguments
! --------------------------------------------------------------------
! loop/variable indices
      INTEGER, intent(in) :: BCV
      INTEGER, intent(in) :: M_TOT
! Local variables
! --------------------------------------------------------------------
      DOUBLE PRECISION :: sum_ep
      INTEGER :: M
! check routine looks that user has not specified both rop_s and ep_s
! (see routine below for details)
      LOGICAL, SAVE :: FIRST_PASS(DIMENSION_BC) = .TRUE.
! --------------------------------------------------------------------
      IF (FIRST_PASS(BCV)) THEN
! JEC: currently a necessary work-around due to a very generic second
! call to check_boundary_conditions following get_stl_data rather than
! a targeted call invoking only which is needed. Ideally, the second
! call would be more targeted or a top-level skip would be enforced for
! bc_types
         FIRST_PASS(BCV) = .FALSE.
      ELSE
         RETURN
      ENDIF

      IF(BC_DT_0(BCV) == UNDEFINED) THEN
         WRITE(ERR_MSG, 1000) trim(iVar('BC_DT_0',BCV))
         CALL LOG_ERROR()
      ENDIF

      DO M = 1, M_TOT
         IF(BC_MASSFLOW_S(BCV,M) /= UNDEFINED .OR. &
            BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1102) trim(iVar('BC_MASSFLOW_S',BCV,M)), &
               trim(iVar('BC_VOLFLOW_S',BCV,M))
            CALL LOG_WARNING()
         ENDIF
 1102 FORMAT('Warning 1102: ', A,' and/or ', A,' have been',/ 'defined',&
         ' at a mass outflow boundary. A specified solids flow ',&
         'rate may',/'not be physically achievable depending on the ',&
         'system and simulation ',/'setup.')
      ENDDO

! If massflow or volflow is specified then require the user set ep_g
! at the boundary. For MO, BC_EP_G is used in the sequence of converting
! gas massflow or volflow to a velocity. Similarly, BC_EP_S is needed
! for converting solids massflow or volflow to a velocity. For an outflow
! condition it should generally be found from the neighboring fluid cell
! and not defined. However, it becomes too challenging for the solids
! phase to accommodate all the possibilities while enforcing consistency.
! The easiest route is simply request the phase volume fraction be set
! for both gas and solids.

      IF(BC_MASSFLOW_G(BCV) /= UNDEFINED .OR. &
         BC_VOLFLOW_G(BCV) /= UNDEFINED ) THEN

         IF(BC_MASSFLOW_G(BCV) == ZERO .OR. &
            BC_VOLFLOW_G(BCV) == ZERO) THEN
            IF(BC_EP_G(BCV) == UNDEFINED) BC_EP_G(BCV) = ZERO
! consistency between massflow/volflow in this case will be checked in
! set_flow_vel
         ELSE
            IF(BC_EP_G(BCV) == UNDEFINED) THEN
               WRITE(ERR_MSG, 1103) trim(iVar('BC_EP_G',BCV))
               CALL LOG_ERROR()
 1103 FORMAT('ERROR 1103: ',A,' was undefined for a MASS_OUTFLOW',/,&
        'boundary with nonzero MASSFLOW_G or VOLFLOW_G. This ',&
        'value ',/'is needed only to determine an IC setting for ',&
        'the conversion to velocity.')
            ENDIF
         ENDIF
      ENDIF ! end if bc_ep_g /= undefined

! if BC_EP_G is defined then set bc_ep_s to ensure consistency at the
! boundary (bc_ep_g must be defined if non zero mass or vol flow)
      IF(BC_EP_G(BCV) /= UNDEFINED) THEN
! If ep_g=1, define bc_eps/bc_rop_s as zero
! (overwrite any other user setting)
         IF(COMPARE(BC_EP_g(BCV),ONE)) THEN
            DO M = 1, M_TOT
               BC_EP_S(BCV,M)  = ZERO
               BC_ROP_S(BCV,M) = ZERO
            ENDDO
         ENDIF
      ENDIF

      DO M = 1,M_TOT
         IF (BC_ROP_S(BCV,M) /= UNDEFINED .AND. &
            BC_ROP_S(BCV,M) /= ZERO) THEN
! Do not allow bc_rop_s in outflow as it is too involved to ensure
! consistency depending on the boundary settings
            WRITE(ERR_MSG, 1106) trim(iVar('BC_ROP_S',BCV,M))
            CALL LOG_ERROR()
         ENDIF
 1106 FORMAT('Error 1106: ',A,' should not be defined for a ',/&
        'MASS_OUTFLOW boundary with nonzero MASSFLOW_S or ',&
        'VOLFLOW_S.',/,'Set BC_EP_S instead.')
      ENDDO

! Calculate the solids volume fraction from the gas phase if there is
! only one solids phase.
      IF(M_TOT == 1 .AND. BC_EP_G(BCV) /= UNDEFINED) THEN
         BC_EP_S(BCV,1) = ONE - BC_EP_g(BCV)
      ENDIF

      DO M = 1,M_TOT
         IF(BC_MASSFLOW_S(BCV,M) /= UNDEFINED .OR. &
            BC_VOLFLOW_S(BCV,M) /= UNDEFINED ) THEN

            IF(BC_MASSFLOW_S(BCV,M) == ZERO .OR. &
               BC_VOLFLOW_S(BCV,M) == ZERO) THEN
               IF(BC_EP_S(BCV,M) == UNDEFINED) BC_EP_S(BCV,M) = ZERO
            ELSE
               IF(BC_EP_S(BCV,M) == UNDEFINED) THEN
                  WRITE(ERR_MSG, 1104) trim(iVar('BC_EP_S',BCV,M))
                  CALL LOG_ERROR()
 1104 FORMAT('ERROR 1104: ',A,' was undefined for a MASS_OUTFLOW',/,&
        'boundary with nonzero MASSFLOWS or VOLFLOW_S. This ',&
        'value ',/'is needed only to determine an IC setting for ',&
        'the conversion to velocity.')
               ENDIF
            ENDIF !non-zero mass/vol flow.
! assign a value to bc_rop_s to satisfy checks; it will get corrected
            BC_ROP_S(BCV,M) = BC_EP_S(BCV,M)

            IF(M_TOT == 1 .AND. BC_EP_G(BCV) == UNDEFINED) THEN
               BC_EP_G(BCV) = ONE - BC_EP_s(BCV,1)
            ENDIF
         ENDIF ! defined mass or vol flow.
      ENDDO ! end do

! It is possible to have bc_ep_g set and not bc_ep_s or vice versa;
! If one is set then the other should be set...
! This does not check if one bc_ep_s is set that all are set; however
! this should be caught in the check on sum of volume fractions
      DO M = 1,M_TOT
         IF ((BC_EP_G(BCV) /= UNDEFINED .AND. &
            BC_EP_S(BCV,M) == UNDEFINED) .OR. &
            (BC_EP_S(BCV,M) /= UNDEFINED .AND. &
            BC_EP_G(BCV) == UNDEFINED)) THEN
            WRITE(ERR_MSG, 1105) BCV
            CALL LOG_ERROR()
         ENDIF
 1105 FORMAT('Error 1105: Incomplete specification for BC', I3,'.',/&
         'BC_EP_G or BC_EP_S is defined but not both. This may lead',/&
         'to inconsistency at the boundary. Please correct.')
      ENDDO

! At this stage if a phase has nonzero massflow or volflow then the
! respective volume fraction should be set

! Now check the sum of the total volume fraction; ensure any user
! settings of volume/void fraction are self-consistent
      IF (BC_EP_G(BCV) /= UNDEFINED) THEN
         SUM_EP = BC_EP_G(BCV)
         DO M = 1, M_TOT

! Add this phase to the total volume fraction.
            SUM_EP = SUM_EP + BC_EP_S(BCV,M)
         ENDDO

! Verify that the volume fractions sum to one.
         IF(.NOT.COMPARE(SUM_EP,ONE)) THEN
            WRITE(ERR_MSG,1107) BCV, trim(iVal(SUM_EP))
            CALL LOG_ERROR()
         ENDIF

 1107 FORMAT('Error 1107: Illegal boundary condition region: ',I3,'. ',&
         'Sum of volume',/'fractions does NOT equal ONE. (SUM = ',A,   &
         ')')
      ENDIF

! Density is needed to convert massflow to velocity. Calculation of
! density in a mass outflow may rely on various dependent variables that
! the user may or may not specify. If they are specified then they
! generally serve only as an IC for density. If they are not specified
! then these values should generally come from the neighboring fluid
! cell.
      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A)

   END SUBROUTINE CHECK_BC_MASS_OUTFLOWA


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_MASS_OUTFLOWB                                   !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message on BC                     !
!                                                                      !
! Comments:                                                            !
!     The velocities at the outflow face are fixed and the momentum    !
!     equations are not solved in the MO cells. Depending on user      !
!     specification (mass flow or volume flow or velocity)             !
!     different checks need to be made. Consistently on the sum of      !
!     the volume fractions at the boundary are made as needed.         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_MASS_OUTFLOWB(M_TOT, SKIP, BCV)

! Modules
!---------------------------------------------------------------------
      use bc, only: bc_ep_g, bc_p_g, bc_t_g, bc_x_g
      use bc, only: bc_rop_s, bc_ep_s, bc_x_s
      use bc, only: bc_massflow_g, bc_massflow_s
      use bc, only: bc_volflow_g, bc_volflow_s
      use param, only: dim_m
      use param1, only: undefined, one, zero
      use physprop, only: ro_g0, ro_s0
      use physprop, only: mw_avg, nmax
      use physprop, only: inert_species
      use run, only: solve_ros, species_eq
      use toleranc, only: compare

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: M_TOT
      LOGICAL, INTENT(in) :: SKIP(DIM_M)
      INTEGER, INTENT(IN) :: BCV

! Local variables
!---------------------------------------------------------------------
! loop/variable indices
      INTEGER :: M, N
      DOUBLE PRECISION SUM_X
! Index of inert species
      INTEGER :: INERT

!---------------------------------------------------------------------

! The check on velocities, among other checks, is made slightly later
! in set_bc0_flow -> check_bc_vel_outflow. Massflow or volflow needs
! to be converted to velocity using known density/volume fraction.
! So such checks are delayed until the calculations can be done.

! For mass outflow with mass or vol flow it is not necessary to set any
! associated boundary condition quantities that are typically needed to
! define a gas density/void fraction as these can be lifted from the
! neighboring fluid cell. However, if these quantities are defined then
! they can supersede the fluid cell value. Note in this case they are
! essentially used as an IC for setting the velocity and are otherwise
! unnecessary. If they are set this routine makes sure they are valid.
      IF(RO_G0 == UNDEFINED) THEN
         IF(BC_MASSFLOW_G(BCV) /= UNDEFINED .AND. &
            BC_MASSFLOW_G(BCV) /= ZERO) THEN
            IF(BC_P_G(BCV) /= UNDEFINED) THEN
               WRITE(ERR_MSG, 1102) trim(iVar('BC_P_g',BCV))
               CALL LOG_WARNING()

               IF(BC_P_G(BCV) <= ZERO) THEN
                  WRITE(ERR_MSG, 1101) trim(iVar('BC_P_G',BCV)),&
                     trim(iVal(BC_P_G(BCV)))
                  CALL LOG_ERROR()
               ENDIF
            ENDIF

            IF(BC_T_G(BCV) /= UNDEFINED) THEN
               WRITE(ERR_MSG, 1102) trim(iVar('BC_T_g',BCV))
               CALL LOG_WARNING()

               IF(BC_T_G(BCV) <= ZERO) THEN
                  WRITE(ERR_MSG, 1101) trim(iVar('BC_T_G',BCV)),&
                     trim(iVal(BC_T_G(BCV)))
                  CALL LOG_ERROR()
               ENDIF
            ENDIF
 1101 FORMAT('Error 1101: If defined, then it must be greater ',&
        'than zero ',/,A,': ',A)

            IF(ANY(BC_X_G(BCV,1:NMAX(0)) /= UNDEFINED)) THEN
               WRITE(ERR_MSG, 1102) trim(iVar('BC_X_g',BCV,N))
               CALL LOG_WARNING()

 1102 FORMAT('Warning 1102: ',A,' was defined for a MASS_OUTFLOW',/,&
        'boundary with nonzero MASSFLOW_G in a compressible system.',/&
        'This value is only used to determine an initial density ',/&
        'in the flow boundary that is then needed in the conversion ',/&
        'of flow to an initial velocity. Velocity is later adjusted ',/&
        'based on the target mass flowrate.')

! if any bc_x_g are defined then ensure consistency
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
            CALL          LOG_ERROR()
 1112 FORMAT('Error 1112: BC_X_g(',I3,',:) do not sum to 1 or 0 ',&
         'and they',/'are not needed.')

                  ELSE
                     BC_X_G(BCV,:) = ZERO
                     BC_X_G(BCV,1) = ONE
                  ENDIF
               ENDIF
            ENDIF !if any bc_x_g are defined

         ENDIF ! if gas massflow rate is set
      ENDIF ! if compressible (variable density)


      DO M = 1,M_TOT
         IF(RO_S0(M) == UNDEFINED) THEN
            IF(BC_MASSFLOW_S(BCV,M) /= UNDEFINED .AND. &
               BC_MASSFLOW_S(BCV,M) /= ZERO) THEN

               IF(ANY(BC_X_S(BCV,M,1:NMAX(M)) /= UNDEFINED)) THEN
                  WRITE(ERR_MSG, 1122) trim(iVar('BC_X_s',BCV,M))
                  CALL LOG_WARNING()

 1122 FORMAT('Warning 1122: ',A,' was defined for a MASS_OUTFLOW',/,&
        'boundary with nonzero MASSFLOW_S in a varying density system.',/&
        'This value may be used to determine an initial density ',/&
        'in the flow boundary that is then needed in the conversion ',/&
        'of flow to an initial velocity. Velocity is later adjusted ',/&
        'based on the target mass flow rate.')

! if any bc_x_s are defined now check consistency
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
               ENDIF ! if any bc_x_s are defined

! check settings for solve_ros
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
               ENDIF ! endif solve_ros

            ENDIF ! endif bc_massflow_s is defined and not zero
         ENDIF ! endif ro_s0 is undefined (variable density)
      ENDDO !end loop over solids



      RETURN
      END SUBROUTINE CHECK_BC_MASS_OUTFLOWB

END MODULE CHECK_BC_OUTFLOW_MOD
