MODULE CALC_MU_G_MOD

   use calc_vol_fr_mod, only: set_ep_factors
   use calc_trd_g_mod, only: calc_deriv_vel_gas
   use compar, only: ijkstart3, ijkend3
   use constant, only: to_si
   use derived_types, only: K_EPSILON_ENUM, MIXING_LENGTH_ENUM, TURBULENCE_MODEL_ENUM
   use derived_types, only: MU_G_MODEL
   use derived_types, only: MU_G_MODEL_ENUM, MU_G_CONSTANT, MU_G_USR
   use derived_types, only: MU_G_SUTHERLAND, MU_G_HERSCHEL_BULKLEY
   use derived_types, only: MU_G_LENNARD_JONES, MU_G_CANTERA_POLY
   use drag, only: f_gs
   use fldvar, only: T_g, epg_ifac, k_turb_g, e_turb_g
   use fldvar, only: ro_g, ep_s, ro_s
   use fldvar, only: x_g
   use functions, only: fluid_at
   use mms, only: use_mms
   use param1, only: zero, one, small_number, undefined
   use physprop, only: mu_g, mu_g0, nmax, mw_g, mw_mix_g
   use physprop, only: SL_muref, SL_Tref, SL_S
   use physprop, only: HB_tau0, HB_k0, HB_n, HB_gama_c
   use physprop, only: LJeps, LJsig, poly_mu_g, mu_species_g
   use run, only: kt_type_enum, ahmadi_1995
   use turb, only: l_scale, tau_1, turb_c_mu
   use usr_prop, only: usr_mug, calc_usr_prop, gas_viscosity
   use visc_g, only: lambda_gt, mu_gmax
   use visc_g, only: mu_gt, epmu_gt, lambda_gt, eplambda_gt
   use visc_s, only: ep_star_array



CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_MU_g                                               C
!  Purpose: Calculate the effective viscosity for a turbulent flow,    C
!           which is the sum of molecular and eddy viscosities         C
!                                                                      C
!                                                                      C
!  Comments: This routine is called even if mu_g0 is defined           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_MU_G()

      IMPLICIT NONE

      SELECT CASE( MU_G_MODEL_ENUM)
      CASE(MU_G_CONSTANT)
      CASE(MU_G_SUTHERLAND)
         CALL CALC_MU_G_SUTHERLAND
      CASE(MU_G_USR)
         CALL CALC_USR_PROP(Gas_Viscosity,lm=0)
      CASE(MU_G_HERSCHEL_BULKLEY)
         CALL CALC_MU_G_HERSCHEL_BULKLEY
      CASE(MU_G_LENNARD_JONES)
         CALL CALC_MU_G_LENNARD_JONES
      CASE(MU_G_CANTERA_POLY)
         CALL CALC_MU_G_CANTERA_POLY
      END SELECT


! adjust viscosity for tubulence
      SELECT CASE( TURBULENCE_MODEL_ENUM)
      CASE(K_EPSILON_ENUM)
         CALL CALC_K_EPSILON_MU
      CASE(MIXING_LENGTH_ENUM)
         CALL CALC_LSCALE_MU
      END SELECT

      CALL SET_EPMUG_VALUES

      RETURN
   END SUBROUTINE CALC_MU_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_EPMUG_VALUES                                        C
!  Purpose: This routine sets the internal variables epmu_g and        C
!  eplambda_g that are used in the stress calculations. If the         C
!  keyword Ishii is invoked then these quantities represent the        C
!  viscosity and second viscosity multiplied by the volume fraction    C
!  otherwise they are simply viscosity/second viscosity (i.e. are      C
!  multiplied by one).                                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE SET_EPMUG_VALUES

! Local variables
!---------------------------------------------------------------------//
! cell index
      INTEGER :: IJK
!---------------------------------------------------------------------//

!      EPMU_GT(:) = EPG_IFAC(:)*MU_gt(:)
!      EPLAMBDA_GT(:) = EPG_IFAC(:)*LAMBDA_gt(:)

! Assign and update ep_g*mu_gt and ep_g*lambda_gt. This is needed even
! for a constant viscosity case

! JFD: Need additional call to SET_EP_FACTORS here
      CALL SET_EP_FACTORS

      DO IJK = ijkstart3, ijkend3
! MMS: Force constant gas viscosity at all cells including ghost cells.
         IF (FLUID_AT(IJK) .OR. USE_MMS) THEN

! if ishii then multiply by void fraction otherwise multiply by 1
            EPMU_GT(IJK)= EPG_IFAC(IJK)*MU_GT(IJK)
            EPLAMBDA_GT(IJK) = EPG_IFAC(IJK)*LAMBDA_GT(IJK)
         ELSE
            EPMU_GT(IJK) = ZERO
            EPLAMBDA_GT(IJK) = ZERO
         ENDIF   ! end if (fluid_at(ijk) .or. use_mms)
      ENDDO   ! end do (ijk=ijkstart3,ijkend3)

      RETURN
   END SUBROUTINE SET_EPMUG_VALUES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Compute a default value of gas viscosity where gas is      C
!  assumed to be air                                                   C
!  Author: W. Sams/M. Syamlal                         Date: 18-JUL-94  C
!                                                                      C
!  Literature/Document References:                                     C
!     Sutherland, W., The Viscosity of Gases and Molecular Force,      C
!        Phil. Mag. 5:507-531, 1893. Eq 4.                             C
!     Dry Air (Ideal Gas State), The NBS-NACA Tables of Thermal        C
!        Properties of Gases, U.S. Department of Commerce, National    C
!        Bureau of Standards, Table 2.39, Morey F. C., 1950.           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_DEFAULT_MUG

      IMPLICIT NONE

! Local parameters
!---------------------------------------------------------------------//
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0

! Local variables
!---------------------------------------------------------------------//
! Cell indices
      INTEGER :: IJK
!---------------------------------------------------------------------//

!!$omp parallel do private(ijk) schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN

! Gas viscosity   (in Poise or Pa.s)
! Calculating gas viscosity using Sutherland's formula with
! Sutherland's constant (C) as 110.0 and reference viscosity of
! air at 273K: mu = 1.71*10-4 poise at T = 273K
            MU_G(IJK) = to_SI*1.7D-4 * &
               (T_G(IJK)/273.0D0)**1.5D0 * (383.D0/(T_G(IJK)+110.D0))

! assign values to quantities used throughout code
            MU_GT(IJK) = MU_G(IJK)
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
         ELSE
            MU_G(IJK) = ZERO
            MU_GT(IJK) = ZERO
            LAMBDA_GT(IJK) = ZERO
         ENDIF
      ENDDO

      RETURN
   END SUBROUTINE CALC_DEFAULT_MUG

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Compute gas viscosity based on Sutherland's law            C
!  Author: Jeff Dietiker                              Date: 11-JAN-24  C
!                                                                      C
!  Literature/Document References:                                     C
!     Sutherland, W., The Viscosity of Gases and Molecular Force,      C
!        Phil. Mag. 5:507-531, 1893. Eq 4.                             C
!     Dry Air (Ideal Gas State), The NBS-NACA Tables of Thermal        C
!        Properties of Gases, U.S. Department of Commerce, National    C
!        Bureau of Standards, Table 2.39, Morey F. C., 1950.           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_MU_G_SUTHERLAND

      IMPLICIT NONE

! Local parameters
!---------------------------------------------------------------------//
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0

! Local variables
!---------------------------------------------------------------------//
! Cell indices
      INTEGER :: IJK
!---------------------------------------------------------------------//

!!$omp parallel do private(ijk) schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN

! Gas viscosity   (in Poise or Pa.s)
! Calculating gas viscosity using Sutherland's law:
            MU_G(IJK) = SL_MUREF * &
               (T_G(IJK)/SL_TREF)**1.5D0 * ((SL_TREF + SL_S)/(T_G(IJK) + SL_S))

! assign values to quantities used throughout code
            MU_GT(IJK) = MU_G(IJK)
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
         ELSE
            MU_G(IJK) = ZERO
            MU_GT(IJK) = ZERO
            LAMBDA_GT(IJK) = ZERO
         ENDIF
      ENDDO

      RETURN
   END SUBROUTINE CALC_MU_G_SUTHERLAND

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: compute gas viscosity based on non-Newtonian               C
!  Herschel-Bulkley model                                              C
!  Author: S. Benyahia                                Date: JAN-11-24  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_MU_G_HERSCHEL_BULKLEY

      IMPLICIT NONE

! Local parameters
!---------------------------------------------------------------------//
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0

! Local variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: D_g(3,3)
! Gas velocity gradient
      DOUBLE PRECISION :: DelV_g(3,3)
! Second invariant of the deviator of D_g
      DOUBLE PRECISION :: I2_devD_g
! cell index
      INTEGER :: IJK
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
! Calculate the rate of strain tensor D_g
            CALL CALC_DERIV_VEL_GAS(ijk, DelV_G, D_G)

! Calculate the second invariant of the deviator of D_g
            I2_DEVD_G = ((D_G(1,1)-D_G(2,2))**2+(D_G(2,2)-D_G(3,3))**2+&
                         (D_G(3,3)-D_G(1,1))**2)/6.D0 + &
                        D_G(1,2)**2 + D_G(2,3)**2 + D_G(3,1)**2
            I2_DEVD_G = 2.0d0*dsqrt(I2_DEVD_G) + small_number

!!	    MU_G(IJK) = K0/I2_DEVD_G**(one-n)  ! singularity at very small I2...

! model from "Progress in Computational Fluid Dynamics, Vol. 9, Nos. 6/7, 2009"
! works for tauy > zero

! Model below hacked from FLUENT 13 Users guide downloaded from internet.
! The model below works for all cases

            if(I2_DEVD_G > HB_gama_c) then  ! tauy/muy = 1 sec^-1
               MU_G(IJK) = HB_Tau0 / I2_DEVD_G + HB_K0*(HB_gama_c/I2_DEVD_G)**(one-HB_n)
            else
               MU_G(IJK) =  HB_Tau0*(2d0-(I2_DEVD_G/HB_gama_c))/HB_gama_c + &
                            HB_K0*((2d0-HB_n)+(HB_n-1d0)*I2_DEVD_G/HB_gama_c)
            endif


! assign values to quantities used throughout code
            MU_GT(IJK) = MU_G(IJK)
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
         ELSE
            MU_G(IJK) = ZERO
            MU_GT(IJK) = ZERO
            LAMBDA_GT(IJK) = ZERO
         ENDIF

      ENDDO

      RETURN
   END SUBROUTINE CALC_MU_G_HERSCHEL_BULKLEY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Compute gas viscosity based on Lennard-Jones model.        C
!  Author: Hang Zhou                                  Date: 20-Apr-25  C
!                                                                      C
!  Literature/Document References:                                     C
!   Bird, R., Stewart, W. and Lightfoot, E. (2002) Transport Phenomena.C
!   2nd Edition, John Wiley and Sons, New York.                        C
!    - Eq.(1.4-14) on Page 26 for the calculation of viscosity         C
!    - Eq.(E.2-1) on Page 866 for the collision integral.              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_MU_G_LENNARD_JONES

      IMPLICIT NONE

! Local parameters
!---------------------------------------------------------------------//
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0

! Local variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: omega_mu, Tstar
      DOUBLE PRECISION :: phi_ij(NMAX(0), NMAX(0)), sum_yphi(NMAX(0))
! Cell indices
      INTEGER :: IJK
! Species indices
      INTEGER :: M, N
!---------------------------------------------------------------------//

!!$omp parallel do private(ijk) schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN

            DO N = 1, nmax(0)
               Tstar = T_g(ijk) / LJeps(N)
               omega_mu = 1.16145D0/Tstar**0.14874D0 + 0.52487D0/exp(0.77320D0*Tstar) &
                        + 2.16178D0/exp(2.43787D0*Tstar)
               ! mu = 5/16*sqrt(pi*m*Kb*T)/(pi*sigma**2*omega)
               ! Kb is the Boltzmann's constant = 1.38066e-23 J/K
               ! m is the molecular mass in kg. m = M/6.02214e23 with M for molecular weight in kg/mol.
               ! The following equation has used all the constant above.
               ! NOTE: LJsig (sigma) is in unit of Angstrom = 1e-10 m. The conversion has been added as well.
               mu_species_g(N) = 2.6696D-6*sqrt(Mw_g(N)*T_g(ijk))/(LJsig(N)*LJsig(N)*omega_mu)
            ENDDO

            CALL CALC_MU_G_MIXTURE(mu_species_g, ijk)

! assign values to quantities used throughout code
            MU_GT(IJK) = MU_G(IJK)
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
         ELSE
            MU_G(IJK) = ZERO
            MU_GT(IJK) = ZERO
            LAMBDA_GT(IJK) = ZERO
         ENDIF
      ENDDO

      RETURN
   END SUBROUTINE CALC_MU_G_LENNARD_JONES

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Compute gas viscosity based on fitted polynomial from      C
!           Cantera based on species transport data.                   C
!  Author: Hang Zhou                                  Date: 20-Apr-25  C
!                                                                      C
!  Literature/Document References:                                     C
!   https://cantera.org/dev/cxx/d8/d58/classCantera_1_1GasTransport.html#a267d20cdea7486fe84e722ee82d84b20  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_MU_G_CANTERA_POLY

      IMPLICIT NONE

! Local parameters
!---------------------------------------------------------------------//
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0

! Local variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: poly_value
! Cell indices
      INTEGER :: IJK
! Species indices
      INTEGER :: N
! Polynominal indices
      INTEGER :: M
!---------------------------------------------------------------------//

!!$omp parallel do private(ijk) schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN

            DO N = 1, nmax(0)
               poly_value = 0.0D0
               DO M = 1,5
                  poly_value = poly_value + poly_mu_g(N,M)*(LOG(T_g(ijk)))**(M-1)
               ENDDO
               mu_species_g(N) = sqrt(T_g(ijk))*(poly_value)**2.0D0
            ENDDO

            CALL CALC_MU_G_MIXTURE(mu_species_g, ijk)

! assign values to quantities used throughout code
            MU_GT(IJK) = MU_G(IJK)
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
         ELSE
            MU_G(IJK) = ZERO
            MU_GT(IJK) = ZERO
            LAMBDA_GT(IJK) = ZERO
         ENDIF
      ENDDO

      RETURN
   END SUBROUTINE CALC_MU_G_CANTERA_POLY

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Compute gas viscosity from species mixture.                C
!  Author: Hang Zhou                                  Date: 20-Apr-25  C
!                                                                      C
! Method of Wilke is used to calculate the mixture viscosity           C
! Reference: Poling, Bruce E., John M. Prausnitz, and John P. O’Connell.C
!            2001. Properties of Gases and Liquids. 5th ed. New York:  C
!            McGraw-Hill Education.                                    C
!            Eqs. 9-5.13-9-5.15.                                       C
!            This is consistent with what is used is BSL.              C
!            Bird, R., Stewart, W. and Lightfoot, E. (2002) Transport  C
!            Phenomena. 2nd Edition, John Wiley and Sons, New York.    C
!            Eqs. 1.4-15 and 1.4-16.                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_MU_G_MIXTURE(mu_species, ijk)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN), DIMENSION(*) :: mu_species
      INTEGER, INTENT(IN) :: ijk

! Local parameters
!---------------------------------------------------------------------//
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0

! Local variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: phi_ij(NMAX(0),NMAX(0)), sum_yphi(NMAX(0))
! Species indices
      INTEGER :: M, N
!---------------------------------------------------------------------//

      sum_yphi = 0.0D0
      DO M = 1, nmax(0)
         DO N = M, nmax(0)
            IF(M==N) THEN
               phi_ij(M,N) = 1.0D0
            ELSE
               phi_ij(M,N) = (1+sqrt(mu_species(M)/mu_species(N))*sqrt(sqrt(mw_g(N)/mw_g(M))))**2 &
                              /sqrt(8*(1+mw_g(M)/mw_g(N)))
            ENDIF
            sum_yphi(M) = sum_yphi(M) + X_g(ijk,N)/MW_G(N)*MW_MIX_G(ijk)*phi_ij(M,N)
         ENDDO
         IF(M .GT. ONE) THEN
            DO N = 1, M-1
               phi_ij(M,N) = phi_ij(N,M)*mu_species(M)/mu_species(N)*mw_g(N)/mw_g(M)
               sum_yphi(M) = sum_yphi(M) + X_g(ijk,N)/MW_G(N)*MW_MIX_G(ijk)*phi_ij(M,N)
            ENDDO
         ENDIF
      ENDDO

      MU_G(ijk) = 0.0D0
      DO M = 1, nmax(0)
         MU_G(ijk) = MU_G(ijk) + X_g(ijk,M)/MW_G(M)*MW_MIX_G(ijk)*mu_species(M)/sum_yphi(M)
      ENDDO

      RETURN
   END SUBROUTINE CALC_MU_G_MIXTURE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: compute turbulent eddy viscosity                           C
!  Author: S. Benyahia                                Date: May-13-04  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!   Cao, J. and Ahmadi, G., 1995, Gas-particle two-phase turbulent     C
!      flow in a vertical duct. Int. J. Multiphase Flow, vol. 21,      C
!      No. 6, pp. 1203-1228.                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_K_EPSILON_MU

      IMPLICIT NONE

! Local parameters
!---------------------------------------------------------------------//
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0

! Local variables
!---------------------------------------------------------------------//
! Solids phase index
      INTEGER :: M
! Constant in turbulent viscosity formulation
      DOUBLE PRECISION :: C_MU
! particle relaxation time
      DOUBLE PRECISION :: Tau_12_st
! cell index
      INTEGER :: IJK
!---------------------------------------------------------------------//

! initialize
      C_MU = turb_c_mu

      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN

! Correction in Ahmadi paper (Cao and Ahmadi)
            IF(KT_TYPE_ENUM == AHMADI_1995 .AND.&
               F_GS(IJK,1) > SMALL_NUMBER) THEN
! solids phase index used throughout routine...
               M = 1 ! for solids phase
               Tau_12_st = Ep_s(IJK,M)*RO_S(IJK,M)/F_GS(IJK,1)
               C_MU = turb_C_MU/(ONE+ Tau_12_st/Tau_1(IJK) * &
                  (EP_s(IJK,M)/(ONE-EP_star_array(IJK)))**3)
            ENDIF

! I'm not very confident about this correction in Peirano paper,
! but it's made available here, uncomment to use it.
! sof@fluent.com --> 02/01/05
!            IF(KT_TYPE_ENUM==SIMONIN_1996 .AND.&
!               F_GS(IJK,1) > SMALL_NUMBER) THEN
! solids phase index used throughout routine...
!               M = 1 ! for solids phase
!               Tau_12_st = Ep_s(IJK,M)*RO_S(IJK,M)/F_GS(IJK,1)
!               X_21 = Ep_s(IJK,M)*RO_S(IJK,M)/(EP_g(IJK)*RO_g(IJK))
! new definition of C_mu (equation A.12, Peirano et al. (2002),
! Powder tech. 122,69-82)
!               IF( K_12(IJK)/(2.0D0*K_Turb_G(IJK)) < ONE) &
!                  C_MU = C_MU/(ONE+ 0.314D0*X_21*Tau_12_st / Tau_1(IJK) * &
!                         (ONE - K_12(IJK)/(2.0D0*K_Turb_G(IJK))) )
!            ENDIF

! Definition of the turbulent viscosity
            MU_GT(IJK) = MU_G(IJK) + RO_G(IJK)*C_MU*&
               K_Turb_G(IJK)**2 / (E_Turb_G(IJK) + SMALL_NUMBER)

            MU_GT(IJK) = MIN(MU_GMAX, MU_GT(IJK))
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
         ELSE
            MU_GT(IJK) = ZERO
            LAMBDA_GT(IJK) = ZERO
         ENDIF
      ENDDO

      RETURN
   END SUBROUTINE CALC_K_EPSILON_MU


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: compute mixing length model                                C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_LSCALE_MU

      IMPLICIT NONE

! Local parameters
!---------------------------------------------------------------------//
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0

! Local variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: D_g(3,3)
! Gas velocity gradient
      DOUBLE PRECISION :: DelV_g(3,3)
! Second invariant of the deviator of D_g
      DOUBLE PRECISION :: I2_devD_g
! cell index
      INTEGER :: IJK
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
! Calculate the rate of strain tensor D_g
            CALL CALC_DERIV_VEL_GAS(ijk, DelV_G, D_G)

! Calculate the second invariant of the deviator of D_g
            I2_DEVD_G = ((D_G(1,1)-D_G(2,2))**2+(D_G(2,2)-D_G(3,3))**2+&
                         (D_G(3,3)-D_G(1,1))**2)/6.D0 + &
                        D_G(1,2)**2 + D_G(2,3)**2 + D_G(3,1)**2

            MU_GT(IJK) =  MU_G(IJK)+2.0*L_SCALE(IJK)*L_SCALE(IJK)*&
                                     RO_G(IJK)*SQRT(I2_DEVD_G)

            MU_GT(IJK) = MIN(MU_GMAX, MU_GT(IJK))
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
         ELSE
            MU_GT(IJK) = ZERO
            LAMBDA_GT(IJK) = ZERO
         ENDIF
      ENDDO

      RETURN
   END SUBROUTINE CALC_LSCALE_MU

END MODULE CALC_MU_G_MOD
