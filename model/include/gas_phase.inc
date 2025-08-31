!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: GAS_PHASE_LOCKED                                          !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / GAS_PHASE_LOCKED /&

      NMAX_g, SPECIES_g, SPECIES_ALIAS_g

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: GAS_PHASE_UNLOCKED                                        !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / GAS_PHASE_UNLOCKED /&

      RO_g0, MU_G_MODEL, MU_g0, &
      SL_muref, SL_Tref, SL_S, &
      HB_tau0, HB_k0, HB_n, HB_gama_c, &
      MW_AVG, MW_g, KG_MODEL, K_g0, C_pg0, DIF_g0, &
      multi_component_diffusion, dif_coeff_kt, dif_thermal, &
      LJsig, LJeps, Dabg, &
      config_transport, mu_transport, alpha_transport, zrot_transport, &
      poly_mu_g, poly_kg
