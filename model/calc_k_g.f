MODULE CALC_K_G_MOD

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE constant, only: GAS_CONST_cal
      USE derived_types, only: KG_MODEL, KG_MODEL_ENUM
      USE derived_types, only: KG_CONSTANT, KG_AIR, KG_USR
      USE derived_types, only: KG_LENNARD_JONES, KG_CANTERA_POLY
      USE fldvar, only: T_g, x_g
      USE functions, only: fluid_at
      USE param1, only: undefined, one, zero
      USE physprop, only: k_g0, k_g, nmax, mw_g, mw_mix_g
      USE physprop, only: LJeps, LJsig, poly_kg, mu_species_g, k_species_g
      USE read_thermochemical, only: calc_cpoR
      USE sendrecv, only: send_recv
! invoke user defined quantity
      USE usr_prop, only: usr_kg, calc_usr_prop
      USE usr_prop, only: gas_conductivity

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_K_g                                                C
!  Purpose: Calculate the effective conductivity of fluid phase        C
!                                                                      C
!                                                                      C
!  Comments:                                                           C
!  This routine will not be called if k_g0 is defined                  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_K_G()

      IMPLICIT NONE
!---------------------------------------------------------------------//

      SELECT CASE( KG_MODEL_ENUM)
      CASE(KG_CONSTANT)
      CASE(KG_AIR)
         CALL CALC_DEFAULT_Kg
      CASE(KG_USR)
         CALL CALC_USR_PROP(Gas_Conductivity,lm=0)
      CASE(KG_LENNARD_JONES)
         CALL CALC_K_G_LENNARD_JONES
      CASE(KG_CANTERA_POLY)
         CALL CALC_K_G_CANTERA_POLY
      END SELECT

      CALL send_recv(K_G, 2)

      RETURN
   END SUBROUTINE CALC_K_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Compute the default value for gas conductivity where the   C
!  gas phase is assumed to be air                                      C
!  Author:M. Syamlal                                  Date: 24-APR-96  C
!                                                                      C
!  Literature/Document References:                                     C
!  Bird, Stewart, and Lightfoot (1960) --                              C
!    Temperature dependence from formula 8.3-12 on p. 255 and          C
!    conductivity value at 300 K from p. 263                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_DEFAULT_Kg

! Modules
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
      USE fldvar, only: T_g
      USE functions, only: fluid_at
      USE param1, only: zero
      USE physprop, only: K_g
      USE run, only: units
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
!---------------------------------------------------------------------//

!!$omp parallel do private(ijk) &
!!$omp& schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
! Gas conductivity (air) in cal/(s.cm.K)
            K_G(IJK) = 6.02D-5*SQRT(T_G(IJK)/300.D0)
         ELSE
            K_G(IJK) = ZERO
         ENDIF

! 1 cal = 4.183925D0 J
         IF (UNITS == 'SI') K_G(IJK) = 418.3925D0*K_G(IJK)      !J/s.m.K

      ENDDO

      RETURN
   END SUBROUTINE CALC_DEFAULT_Kg

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Compute gas thermal conductivity based on                  C
!           Lennard-Jones model.                                       C
!  Author: Hang Zhou                                  Date: 20-Apr-25  C
!                                                                      C
!  Literature/Document References:                                     C
!   Bird, R., Stewart, W. and Lightfoot, E. (2002) Transport Phenomena.C
!   2nd Edition, John Wiley and Sons, New York.                        C
!- Eq.(9.3-13) on Page 275 for the calculation of thermal conductivity C
!- Eq.(E.2-1) on Page 866 for the collision integral.                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_K_G_LENNARD_JONES

      IMPLICIT NONE

! Local parameters
!---------------------------------------------------------------------//

! Local variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: omega_k, Tstar, mu
      DOUBLE PRECISION :: phi_ij(NMAX(0), NMAX(0)), sum_yphi(NMAX(0))
! Cell indices
      INTEGER :: IJK
! Species indices
      INTEGER :: M, N
! Gas constant
      DOUBLE PRECISION :: GAS_CONSTANT
!---------------------------------------------------------------------//

      GAS_CONSTANT = 4.183925d3 * GAS_CONST_cal !(J/kmol.K)

!!$omp parallel do private(ijk) schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN

            DO N = 1, nmax(0)
               Tstar = T_g(ijk) / LJeps(N)
               omega_k = 1.16145D0/Tstar**0.14874D0 + 0.52487D0/exp(0.77320D0*Tstar) &
                        + 2.16178D0/exp(2.43787D0*Tstar)
               ! k = 25/32*sqrt(pi*m*Kb*T)/(pi*sigma**2*omega)*Cv = 15/4*R/M*mu = 5/2*Cv*mu for monatomic gas
               !     Kb is the Boltzmann's constant = 1.38066e-23 J/K
               !     m is the molecular mass in kg. m = M/6.02214e23 with M for molecular weight in kg/mol.
               !     Cv = 3/2*R/M for monatomic gas. R is the gas constant = 8.31451 J/mol/K.
               !     It can give k_species(N) = 8.3236D-2*sqrt(T_g(ijk)/Mw_g(N))/(LJsig(N)*LJsig(N)*omega_k)
               ! For polyatomic gas, k = (Cp+5/4*R/M)*mu.
               !     This works for monatomic gas as well because Cp = 5/2* R/M for monatomic gases.
               ! Therefore, we are using the format for polyatomic gas here for k.
               ! For each species, we only have LJ model or cantera-based polynomial models for viscosity (mu).
               ! So, to be consistent, we replace mu in the equation with the LJ model for viscosity
               !     when the species viscosity is not calculated.
               IF(mu_species_g(N) .EQ. undefined) THEN
                  mu = 2.6696D-6*sqrt(Mw_g(N)*T_g(ijk))/(LJsig(N)*LJsig(N)*omega_k)
                  k_species_g(N) = (calc_cpoR(T_g(IJK),0,N)+1.25d0)*GAS_CONSTANT/Mw_g(N) * mu
               ELSE
                  k_species_g(N) = (calc_cpoR(T_g(IJK),0,N)+1.25d0)*GAS_CONSTANT/Mw_g(N) * mu_species_g(N)
               ENDIF
            ENDDO

            CALL CALC_K_G_MIXTURE_LJ(k_species_g, ijk)

         ELSE
            K_G(IJK) = ZERO
         ENDIF

      ENDDO

      RETURN
   END SUBROUTINE CALC_K_G_LENNARD_JONES

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Compute gas thermal conductivity based on fitted polynomialC
!           from Cantera based on species transport data.              C
!  Author: Hang Zhou                                  Date: 20-Apr-25  C
!                                                                      C
!  Literature/Document References:                                     C
!   https://cantera.org/dev/cxx/d8/d58/classCantera_1_1GasTransport.html#a267d20cdea7486fe84e722ee82d84b20  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_K_G_CANTERA_POLY

      IMPLICIT NONE

! Local parameters
!---------------------------------------------------------------------//

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
                  poly_value = poly_value + poly_kg(N,M)*(LOG(T_g(ijk)))**(M-1)
               ENDDO
               k_species_g(N) = sqrt(T_g(ijk))*poly_value
            ENDDO

            CALL CALC_K_G_MIXTURE_CANTERA(k_species_g, ijk)
         ELSE
            K_G(IJK) = ZERO
         ENDIF

      ENDDO

      RETURN
   END SUBROUTINE CALC_K_G_CANTERA_POLY

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Compute gas thermal conductivity from species mixture.     C
!  Author: Hang Zhou                                  Date: 20-Apr-25  C
!                                                                      C
! Method of Wilke is used to calculate the mixture viscosity           C
! Reference: Poling, Bruce E., John M. Prausnitz, and John P. O’Connell.C
!            2001. Properties of Gases and Liquids. 5th ed. New York:  C
!            McGraw-Hill Education.                                    C
!            Eqs. 10-6.1, 10-6.2.                                      C
!                                                                      C
!            This is consistent with what is used is BSL.              C
!            Bird, R., Stewart, W. and Lightfoot, E. (2002) Transport  C
!            Phenomena. 2nd Edition, John Wiley and Sons, New York.    C
!            Eq. 9.3-18                                                C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_K_G_MIXTURE_LJ(k_species, ijk)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN), DIMENSION(*) :: k_species
      INTEGER, INTENT(IN) :: ijk

! Local parameters
!---------------------------------------------------------------------//

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
               phi_ij(M,N) = (1+sqrt(k_species(M)/k_species(N))*sqrt(sqrt(mw_g(N)/mw_g(M))))**2 &
                              /sqrt(8*(1+mw_g(M)/mw_g(N)))
            ENDIF
            sum_yphi(M) = sum_yphi(M) + X_g(ijk,N)/MW_G(N)*MW_MIX_G(ijk)*phi_ij(M,N)
         ENDDO
         IF(M .GT. ONE) THEN
            DO N = 1, M-1
               phi_ij(M,N) = phi_ij(N,M)*k_species(M)/k_species(N)*mw_g(N)/mw_g(M)
               sum_yphi(M) = sum_yphi(M) + X_g(ijk,N)/MW_G(N)*MW_MIX_G(ijk)*phi_ij(M,N)
            ENDDO
         ENDIF
      ENDDO

      K_G(ijk) = 0.0D0
      DO M = 1, nmax(0)
         K_G(ijk) = K_G(ijk) + X_g(ijk,M)/MW_G(M)*MW_MIX_G(ijk)*k_species(M)/sum_yphi(M)
      ENDDO

      RETURN
   END SUBROUTINE CALC_K_G_MIXTURE_LJ

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Compute gas thermal conductivity from species mixture.     C
!  Author: Hang Zhou                                  Date: 20-Apr-25  C
!                                                                      C
! Symmetric mean mixing rule based on mole fraction weighting and      C
! harmonic averaging is used to calculate mixture thermal condictiity. C
! This is consistent with what Cantera used for MiXTransport.         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_K_G_MIXTURE_CANTERA(k_species, ijk)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN), DIMENSION(*) :: k_species
      INTEGER, INTENT(IN) :: ijk

! Local parameters
!---------------------------------------------------------------------//

! Local variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: sum1, sum2
! Species indices
      INTEGER :: N
!---------------------------------------------------------------------//

      sum1 = 0.0d0
      sum2 = 0.0d0

      DO N = 1, nmax(0)
         sum1 = sum1 + X_g(ijk,N)/MW_G(N)*MW_MIX_G(ijk) * k_species(N)
         sum2 = sum2 + X_g(ijk,N)/MW_G(N)*MW_MIX_G(ijk) / k_species(N)
      ENDDO

      K_G(ijk) = 0.5D0*(sum1+1.0D0/sum2)

      RETURN
   END SUBROUTINE CALC_K_G_MIXTURE_CANTERA

END MODULE CALC_K_G_MOD
