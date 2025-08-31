!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: CHEM_EQUATIONS_LOCKED                                     !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
       NAMELIST / CHEM_EQUATIONS_LOCKED /&
       STIFF_CHEMISTRY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: CHEM_EQUATIONS_UNLOCKED                                   !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / CHEM_EQUATIONS_UNLOCKED /&

      ARRHENIUS_RRATES_FLUID, ARRHENIUS_RRATES_DES, &
      STIFF_CHEMISTRY, STIFF_CHEM_MAX_STEPS, &
      STIFF_CHEM_ABS_TOL, STIFF_CHEM_REL_TOL, &
      STIFF_CHEM_MIN_RATE, STIFF_CHEM_MIN_RATE_DPM, &
      CHEM_MIN_SPECIES_FLUID, CHEM_MIN_SPECIES_SOLID, &
      DES_MIN_PMASS_FRAC
