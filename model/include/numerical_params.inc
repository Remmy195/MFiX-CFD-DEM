!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: NUMERICAL_PARAM_LOCKED                                    !
!  Author: J.Musser                                   Date: 03-APR-15  !
!    Charles G. Waldman                                     23-APR-25  !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
       NAMELIST / NUMERICAL_PARAM_LOCKED /&
       REPEATABLE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: NUMERICAL_PARAM_UNLOCKED                                  !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / NUMERICAL_PARAM_UNLOCKED /&

      MAX_NIT, DETECT_STALL, NORM_G, NORM_S, MAX_INLET_VEL_FAC, &
      TOL_RESID, TOL_RESID_T, TOL_RESID_X, TOL_RESID_TH, &
      TOL_RESID_SCALAR, TOL_DIVERGE, TOL_SUM_RESID_ABS, &
      LEQ_IT, LEQ_METHOD, LEQ_SWEEP, LEQ_TOL, LEQ_PC, DO_TRANSPOSE, &
      UR_FAC, UR_F_GS, UR_KTH_SML, &
      DISCRETIZE, DEF_COR, CHI_SCHEME, CN_ON, &
      DIL_EP_S, ZERO_EP_s, ZERO_X_gs, &
      EPP_DEN, PPG_DEN, TMIN, TMAX, &
      FPE_TRAP_INVALID, FPE_TRAP_ZERO, FPE_TRAP_OVERFLOW, &
      REPEATABLE
