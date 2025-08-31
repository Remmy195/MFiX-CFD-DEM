MODULE CALC_DIF_g_MOD
CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DIF_g                                              C
!  Purpose: Calculate the effective diffusivity of fluid phase in      C
!  dimensions length^2/time.                                           C
!                                                                      C
!                                                                      C
!  Comments:                                                           C
!  This routine will not be called if dif_g is defined                 C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_DIF_G()

! Modules
!---------------------------------------------------------------------//
      USE param1, only: undefined
      USE physprop, only: dif_g0, dif_g, divJi, Multi_Component_Diffusion
      USE sendrecv, only: send_recv
! invoke user defined quantity
      USE usr_prop, only: usr_difg, calc_usr_prop
      USE usr_prop, only: gas_diffusivity
      IMPLICIT NONE
!---------------------------------------------------------------------//

      IF (USR_Difg) THEN
         CALL CALC_USR_PROP(Gas_Diffusivity,lm=0)
      ELSEIF (Multi_Component_Diffusion) THEN
         CALL MultiComponentDiffusion
      ELSEIF (Dif_g0 == UNDEFINED) THEN
! unnecessary check but included for clarity
         CALL CALC_DEFAULT_DIF_GAS
      ENDIF

      CALL send_recv(DIF_G, 2)
      IF(Multi_Component_Diffusion) CALL send_recv(divJi, 2)

      RETURN
   END SUBROUTINE CALC_DIF_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DEFAULT_DIF_GAS                                    C
!  Purpose: Compute the default value for diffusivity of each gas      C
!  species; each species is assigned the same value with the base      C
!  value corresponding to CO2 in N2 at T=298K and P~=1atm.             C
!                                                                      C
!  Author:M. Syamlal                                  Date: 13-FEB-98  C
!                                                                      C
!  Revision: Include dilute mixture approximation for calculation of   C
!  multicomponent diffusion coefficients                               C
!  Author:N. Reuge                                    Date: 11-APR-07  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Fuller relation - to account for influence of gas temperature       C
!     and pressure                                                     C
!  Curtiss-Hirschfelder, Wilke & Blanc - dilute mixture approximation  C
!     for multicomponent diffusion. Valid if the mass fraction of the  C
!     carrier species > 0.9                                            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_DEFAULT_DIF_GAS

! Modules
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
      use fldvar, only: T_g, P_g, X_g
      use functions, only: fluid_at
      use param1, only: zero, one
      use physprop, only: NMAX, Dif_g
      use scales, only: unscale_pressure
      use toleranc, only: zero_x_gs
      IMPLICIT NONE

! Local Variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
! Species index
      INTEGER :: N, L
! Binary diffusion coefficient
      DOUBLE PRECISION :: Dab(NMAX(0),NMAX(0))
! Reference temperature and pressure for each species diffusion
! coefficient
      DOUBLE PRECISION :: Tg_ref(NMAX(0))
      DOUBLE PRECISION :: Pg_ref(NMAX(0))
! Intermediate calculation to determine weighted average diffusion
! coefficient
      DOUBLE PRECISION :: lDab, Sum_XgLoDab, Sum_XgL, lXgN
!---------------------------------------------------------------------//

      CALL SET_BINARY_DAB_GAS(Dab, Tg_ref, Pg_ref)

!!$omp  parallel do private(ijk) &
!!$omp& schedule(dynamic,chunk_size)
      DO N = 1, NMAX(0)
         DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN

               IF (NMAX(0) == 1) THEN
! Only 1 species present in gas phase
                  lDab = Dab(N,N)

               ELSE
! Dilute mixture approximation for multi-component diffusion

                  SUM_XgLoDab = ZERO
                  SUM_XgL = ZERO
                  DO L = 1, NMAX(0)
                     IF (L /= N) THEN
                        SUM_XgLoDab = SUM_XgLoDab+&
                                     X_g(IJK,L)/Dab(N,L)
! Sum_XgL = 1-XgN
                        SUM_XgL = SUM_XgL + X_g(IJK,L)
                     ENDIF
                  ENDDO
! it is possible for lXgN to evaluate <0
                  lXgN = ONE-SUM_XgL

                  IF (lXgN > ZERO_X_GS .AND. &
                     SUM_XgL > ZERO_X_GS) THEN
! i.e. when cell is not only species N
! If this criteria is too strict (i.e. when XgN->1), then this section
! may be evaluated when the other XgN do not carry significant value
! (i.e. are effectively zero). As a result, lDab may be evaluated and
! become unrealistically large. This may happen when happen when
! 1-XgN != sum_XgL. Generally, sum_XgL should equal 1-XgN but they may
! differ due to the numerical evaluation of each Xg. If XgN->1, then
! use both sum_XgL and 1-XgN to determine whether calculations should
! proceed. Given the noted limitation of this approximation then an
! additional criteria should probably be added to only evaluate when
! sum_xgL <0.1.
                     IF (SUM_XgLoDab > ZERO) THEN
! for numerical reasons use Sum_XgL rather than ONE-X_g(IJK,N)
                        lDab = SUM_XgL/SUM_XgLoDab
                     ELSE
! this should not occur...
                        lDab = Dab(N,N)
                     ENDIF
                  ELSE
! Address case when the mass fraction of the Nth species is nearly 1.
                     lDab = Dab(N,N)
                  ENDIF
               ENDIF

! Influence of gas temperature and gas pressure from Fuller relation
               DIF_G(IJK,N) = lDab* (T_g(IJK)/Tg_ref(N))**1.75* &
                           Pg_ref(N)/UNSCALE_PRESSURE(P_g(IJK))
            ELSE
               DIF_G(IJK,N) = ZERO
            ENDIF
         ENDDO
      ENDDO

      RETURN
   END SUBROUTINE CALC_DEFAULT_DIF_GAS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_BINARY_DAB_GAS                                      C
!  Purpose: Set binary diffusion coefficients for all cross species    C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Bird, Stewart and lightfoot (1960)                                  C
!  Reid, Prausnitz and Poling (1987)                                   C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE SET_BINARY_DAB_GAS(lDab, lTg_ref, lPg_ref)

! Modules
!---------------------------------------------------------------------//
      use physprop, only: NMAX
      use run, only: units
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Binary diffusion coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDab(NMAX(0),NMAX(0))
      DOUBLE PRECISION, INTENT(OUT) :: lTg_ref(NMAX(0))
      DOUBLE PRECISION, INTENT(OUT) :: lPg_ref(NMAX(0))

! Local Variables
!---------------------------------------------------------------------//
! Species index
      INTEGER :: N, L
!---------------------------------------------------------------------//


! Default gas diffusion coefficient
! Bird, Stewart, and Lightfoot (1960) -- CO2--N2 at 298.2 K
      DO N = 1,NMAX(0)
         DO L = N,NMAX(0)
!         DO L = N+1, NMAX(0)
!            IF (N /= L) THEN
! assign the diagonal case for when nmax = 1 or when only species N present
! in given cell.
               lDab(N,L) = 0.165D0       !cm^2/s
               lDab(L,N) = lDab(N,L)
!            ENDIF
         ENDDO
         lTg_ref(N) = 298.2d0
         lPg_ref(N) = 1.01D6             !dyne
      ENDDO

! Gas diffusion coefficients for 3 species system: SiH4, H2 & N2
! Calculated using relation derived from Chapman and Enskog's theory
! of gases - Reid, Prausnitz and Poling (1987)
! Binary diffusion coefficient SiH4/H2 at 873 K
!      lDab(1,2) = 3.78      ! cm^2/s
!      lDab(2,1) = lDab(1,2)
! Binary diffusion coefficient SiH4/N2 at 873 K
!      lDab(1,3) = 1.02      ! cm^2/s
!      lDab(3,1) = lDab(1,3)
! Binary diffusion coefficient H2/N2 at 873 K
!      lDab(2,3) = 4.52      ! cm^2/s
!      lDab(3,2) = lDab(2,3)
!      DO N = 1,NMAX(0)
!         lTg_ref(N) = 873.0
!         lPg_ref(N) = 1.01e6         ! dyne
!      ENDDO

      IF(UNITS == 'SI') THEN
         DO N = 1,NMAX(0)
            DO L = N,NMAX(0)
!               IF (N /= L) THEN
! assign the diagonal case nmax = 1 or when only species N present
! in given cell.
                  lDab(N,L) = lDab(N,L)*0.0001D0   !m^2/s
                  lDab(L,N) = lDab(N,L)
!               ENDIF
            ENDDO
            lPg_ref(N) = lPg_ref(N)/10.D0          !Pa
         ENDDO
      ENDIF

      RETURN
   END SUBROUTINE SET_BINARY_DAB_GAS

END MODULE CALC_DIF_g_MOD

!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: MultiComponentDiffusion                                 C
!  Purpose: Compute the Fickian coefficients based on computed Stefan- C
!           Maxwell (ordinary) diffusion coef. based on kinetic theory C
!  Author: S. Benyahia                                 Date: 26-OCT-24 C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  BSL, Bird Stewart and Lightfoot  Transport phenomena, 2002          C
!  2nd ed.  Wiley.                                                     C
!                                                                      C
!  Taylor, R & Krishna, R - Multicomponent mass transfer, 1993, Wiley  C
!     for multicomponent diffusion expressed in mole fraction.         C
!                                                                      C
!  Fluent 13 manual for mass flux expression in mass fraction. Those   C
!   expressions were rederived in the context of this work and can     C
!   also be found in BSL textbook.                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE MultiComponentDiffusion

! Modules
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
      use constant, only: to_SI
      use cooling_rate_mod, only: ludcmp, lubksb
      use cutcell, only: cut_cell_at, cut_treatment_at
      use fldvar, only: T_g, P_g, X_g, RO_g, ROP_g
      use fun_avg, only: avg_x, avg_y, avg_z, avg_x_s, avg_y_s, avg_z_s
      use fun_avg, only: avg_x_h, avg_y_h, avg_z_h
      use functions, only: fluid_at, im_of, jm_of, km_of, ip_of, jp_of, kp_of
      use functions, only: east_of, west_of, north_of, south_of, top_of, bottom_of, im_of, jm_of, km_of
      use geometry
      use indices, only: i_of, j_of, k_of, im1, jm1, km1
      use param1, only: zero, half, one, small_number
      use param, only: dimension_3
      use physprop, only: NMAX, Dif_g, MW_g, mw_mix_g
      use physprop, only: divJi, dif_coeff_kt, dif_thermal, LJsig, LJeps, Dabg
      use run, only: energy_eq
      use scales, only: unscale_pressure
      use toleranc, only: zero_x_gs
      IMPLICIT NONE

! Local Variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, IJKE, IJKN, IJKT, IJKW, IJKS, IJKB
      INTEGER :: IJKP, IJPK, IPJK
      INTEGER :: IM, JM, KM, IMJK, IJMK, IJKM
! Species index
      INTEGER :: i, j, k, kk, ss, s, indx((NMAX(0)-1))
! Binary diffusion coefficient Fickian (Dab) and Stefan-Maxwell (lDab)
      DOUBLE PRECISION :: Dab(NMAX(0),NMAX(0)), lDab(NMAX(0),NMAX(0))
      DOUBLE PRECISION :: Aij((NMAX(0)-1),(NMAX(0)-1)), Bij((NMAX(0)-1),(NMAX(0)-1))
      DOUBLE PRECISION :: Amat0((NMAX(0)-1),(NMAX(0)-1)), Bmat0((NMAX(0)-1))
      DOUBLE PRECISION :: fDif_g(DIMENSION_3, (NMAX(0)-1), (NMAX(0)-1))
      DOUBLE PRECISION :: DTi(DIMENSION_3, (NMAX(0)-1))  ! thermal diffusion coefficient
      DOUBLE PRECISION :: C_AE, C_AW, C_AN, C_AS, C_AT, C_AB
! miscel
      DOUBLE PRECISION :: Mt, SumAii, SumFluxX, SumFluxY, SumFluxZ, sumMx1, sumMx2

      DOUBLE PRECISION :: D_Fe, D_Fw, D_Fn, D_Fs, D_Ft, D_Fb, xg(NMAX(0))
!
! variables to compute diffusion coefficients from KT
      DOUBLE PRECISION :: Pg, sig_ab, eps_ab, Tstar, Ohmega_ab

      double precision :: d  ! most of these parameters will come out of LU
      INTEGER :: NP          !max no. of linear equations to solve for species.
!
!---------------------------------------------------------------------//
!
! Initialize & set self diffusion (diagonal) to zero

      lDab = zero

! Set binary diffusion coefficient read from input file
      if(.not.dif_coeff_kt) then !
	 do I = 1, nmax(0)
	    dO J = (I+1), nmax(0)
	       lDab(i,j) = Dabg(i,j)
	       lDab(j,i) = lDab(i,j)  ! must be symmetric
	    enddO
	 enddo
      endif
!
      S = NMAX(0)-1   ! solving for (n-1) species
      DTi(:,:) = zero
      DIF_G(:,:) = zero
!
!!$omp  parallel do private(ijk)
!!$omp& schedule(dynamic,chunk_size)
      DO IJK = IJKSTART3, IJKEND3
         IF (FLUID_AT(IJK)) THEN
	    Mt = zero
	    do I = 1, nmax(0)
	      Mt = Mt + X_G(ijk,i)/MW_G(i)
	    enddo
	    Mt = one/Mt

	    SumMx1 = zero
	    SumMx2 = zero
	    do I = 1, nmax(0)
	      xg(I) = X_G(ijk,i)*Mt/MW_G(i)  ! calculate mole fraction
	      SumMx1 = SumMx1 + xg(I)*MW_G(i)**0.511d0  ! sums needed for thermal diffusion
	      SumMx2 = SumMx2 + xg(I)*MW_G(i)**0.489d0
	    enddo
!
! compute thermal diffusion coefficients
!
            if(dif_thermal) THEN
	       do i = 1, s
	          DTi(ijk,i) = -2.59d-7/T_g(ijk)**0.341d0 *   &   ! (0.659-1=-0.341) includes T in dT/T
	                      (xg(i)*MW_G(i)**0.511d0/SumMx1 - X_G(ijk,i)) *(SumMx1/SumMx2)
	       enddo
            endif
!
! compute diffusion coefficient from KT of low density gases
	    if(dif_coeff_kt) then !
	    do I = 1, nmax(0)
	       dO J = (I+1), nmax(0)
! When using thermodynamic pressure, as in below, we must use a compressible flow (i.e. ro_g .ne. const)
		 Pg = UNSCALE_PRESSURE(P_g(ijk))/(1.013d6*to_SI)  ! from Pa/dyne to atm
		 sig_ab = half*(LJsig(i) + LJsig(j))
		 eps_ab = dsqrt(LJeps(i) * LJeps(j))
		 Tstar = T_g(ijk) / eps_ab
! collision integral for diffusion Eq. E 2.2 Page 866 in BSL
		 Ohmega_ab = 1.06036d0/Tstar**0.1561d0 + 0.193d0/dEXP(0.47635d0*Tstar) + &
		             1.03587d0/dEXP(1.52996d0*Tstar) + 1.76474d0/dEXP(3.89411d0*Tstar)
! compute ordinary diffusion coefficient from KT, eq. 17.3-12 Page 526 in BSL
! as mentioned in the BSL textbook, the Dab obtained below is within 6-10% of measured values.
		 lDab(i,j) = 1.8583d-3*dsqrt(T_g(ijk)**3 * (one/Mw_g(i)+one/Mw_g(j))) / &
		             (Pg * sig_ab**2 * Ohmega_ab)
		 lDab(i,j) = 1d-4*lDab(i,j)  ! from cm^2/sec to m^2/sec
		 lDab(j,i) = lDab(i,j) ! symmetric matrix for Stefan-Maxwell coefficients
!
	       enddO
	    enddo
            endif ! change later this keyword to Dif_Coeff_KT
! end of Dab calculations from kinetic theory
!
! Now setting up the matrices Aij and Bij to calculate the Fickian
! diffusion coefficients (non-symmetric matrix fDif_g)
!
	    DO I = 1, s
	       SumAii = zero
	       DO J = 1, NMAX(0)
	         IF(J/=I) SumAii = SumAii + (xg(J)*Mt/(lDab(I,J)*Mw_g(I)))
	       ENDDO
	       DO J = 1, s
	         IF(J==I) Aij(I,J) = (xg(I)*Mt/(lDab(I,NMAX(0))*Mw_g(NMAX(0))) &
		                      + SumAii)
	         IF(J/=I) Aij(I,J) = xg(I)*(Mt/(lDab(I,NMAX(0))*Mw_g(NMAX(0)))-&
		                                 Mt/(lDab(I,J)*Mw_g(J)))
	         IF(J==I) Bij(I,J) = (xg(I)*Mt/Mw_g(NMAX(0)) + (one-xg(I))&
		                                *Mt/Mw_g(I))
	         IF(J/=I) Bij(I,J) = xg(I)*(Mt/Mw_g(NMAX(0)) - Mt/Mw_g(J))
	       ENDDO
            ENDDO
!
! Solving the system of linear equations. The details of the LU decomposition
! routines can be found in cooling_rate.f GHD theory in MFIX.
!
      do kk=1,s
         do i=1,s
            do j=1,s
                Amat0(i,j) = Aij(i,j)
            enddo
         enddo
         do i=1,s
            bmat0(i) = Bij(i,kk)
         enddo
         NP = NMAX(0)  ! max number of linear equations solved below
         CALL LUDCMP(Amat0, s, NP, indx, d, 'multiCompDiff') ! solve system of (nmax(0)-1) linear
         CALL LUBKSB(Amat0, s, NP, indx, bmat0)              ! equations using LU decomposition

         do i=1,s
            Dab(i,kk) = bmat0(i)   !! just (n-1) as we do not solve for n species equation !!
	    fDif_g(IJK,i,kk) = Dab(i,kk)*ROP_g(IJK)  ! include gas density as done in solve_species_eq
	    if(dabs(fDif_g(IJK,i,kk)) < small_number) fDif_g(IJK,i,kk) = zero  ! avoid small numbers
         enddo
      enddo
!
!  End of call to LU decomposition routines...
!
           do i=1,s
	     DIF_G(IJK,i) = Dab(i,i)  ! used in conv_dif routine; important for convergence
	     fDif_g(IJK,i,i) = zero   ! only off-diagonal values used in flux calculation below
           enddo
         ELSE  ! if not fluid at ijk
	   DIF_G(IJK,:) = zero
	   fDif_g(IJK,:,:) = zero
	   DTi(ijk,:) = zero
         ENDIF
      ENDDO   ! end do ijk
!
! end of Fickian diffusion coef. calculation
!
      divJi(:,:) = zero      ! species off-diagonal fluxes
!
! compute diffusion fluxes at cell faces, code hacked from conv_dif_phi
!
!!$omp  parallel do private(ijk)
      DO IJK = IJKSTART3, IJKEND3
	IF (FLUID_AT(IJK)) THEN

         IMJK = IM_OF(IJK)
         IJMK = JM_OF(IJK)
         IJKM = KM_OF(IJK)


         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IM = IM1(I)
         JM = JM1(J)
         KM = KM1(K)

         IJKE = EAST_OF(IJK)
         IJKN = NORTH_OF(IJK)
         IJKS = SOUTH_OF(IJK)
         IJKW = WEST_OF(IJK)

         C_AE = ODX_E(I)*AYZ(IJK)
         C_AW = ODX_E(IM)*AYZ(IMJK)
         C_AN = ODY_N(J)*AXZ(IJK)
         C_AS = ODY_N(JM)*AXZ(IJMK)
         C_AT = OX(I)*ODZ_T(K)*AXY(IJK)
         C_AB = OX(I)*ODZ_T(KM)*AXY(IJKM)

         IF(CUT_TREATMENT_AT(IJK).AND.CUT_CELL_AT(IJK)) THEN
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKP = KP_OF(IJK)

            IF (.NOT.FLUID_AT(IPJK)) C_AE = ODX_E(I)*DY(J)*DZ(K)
            IF (.NOT.FLUID_AT(IMJK)) C_AW = ODX_E(IM)*DY(J)*DZ(K)
            IF (.NOT.FLUID_AT(IJPK)) C_AN = ODY_N(J)*DX(I)*DZ(K)
            IF (.NOT.FLUID_AT(IJMK)) C_AS = ODY_N(JM)*DX(I)*DZ(K)
            IF (.NOT.FLUID_AT(IJKP)) C_AT = OX(I)*ODZ_T(K)*DX(I)*DY(J)
            IF (.NOT.FLUID_AT(IJKM)) C_AB = OX(I)*ODZ_T(KM)*DX(I)*DY(J)
         ENDIF

	 Do ss = 1, s
	    SumFluxX = zero
	    SumFluxY = zero
	    SumFluxZ = zero
	    Do kk = 1, s
	      if(kk /= ss) then  ! only for off-diagonal fluxes
! note off-diagonal diffusion coefficients can be negative, thus using AVG_X_S

                D_Fe = AVG_X_S(fDif_g(IJK,ss,kk),fDif_g(IJKE,ss,kk),I)*C_AE*( X_g(IJKE, kk) - X_g(IJK, kk) )

                D_FW = AVG_X_S(fDif_g(IJKW,ss,kk),fDif_g(IJK,ss,kk),IM)*C_AW*( X_g(IJK, kk) - X_g(IJKW, kk) )

                D_FN = AVG_Y_S(fDif_g(IJK,ss,kk),fDif_g(IJKN,ss,kk),J)*C_AN*( X_g(IJKN, kk) - X_g(IJK, kk) )

                D_FS = AVG_Y_S(fDif_g(IJKS,ss,kk),fDif_g(IJK,ss,kk),JM)*C_AS*( X_g(IJK, kk) - X_g(IJKS, kk) )

                IF (DO_K) THEN

		  IJKT = TOP_OF(IJK)
                  IJKB = BOTTOM_OF(IJK)

                  D_FT = AVG_Z_S(fDif_g(IJK,ss,kk),fDif_g(IJKT,ss,kk),K)*C_AT*( X_g(IJKT, kk) - X_g(IJK, kk) )
                  D_FB = AVG_Z_S(fDif_g(IJKB,ss,kk),fDif_g(IJK,ss,kk),KM)*C_AB*( X_g(IJK, kk) - X_g(IJKB, kk) )

		ENDIF

	        SumFluxX = SumFluxX + (D_Fe - D_FW)
	        SumFluxY = SumFluxY + (D_FN - D_FS)
	        IF (DO_K) SumFluxZ = SumFluxZ + (D_FT - D_FB)

	      endif ! only off diagonal fluxes computed here
	    endDo  ! for kk loop

	    if(dif_thermal) then  ! include thermal diffusion contribution to Ji
                D_Fe = AVG_X_S(DTi(IJK,ss),DTi(IJKE,ss),I)*C_AE*( T_g(IJKE) - T_g(IJK) ) ! /T included in DTi
                D_FW = AVG_X_S(DTi(IJKW,ss),DTi(IJK,ss),IM)*C_AW*( T_g(IJK) - T_g(IJKW) )
                D_FN = AVG_Y_S(DTi(IJK,ss),DTi(IJKN,ss),J)*C_AN*( T_g(IJKN) - T_g(IJK) )
                D_FS = AVG_Y_S(DTi(IJKS,ss),DTi(IJK,ss),JM)*C_AS*( T_g(IJK) - T_g(IJKS) )
		IF (DO_K) THEN
                  IJKT = TOP_OF(IJK)
                  IJKB = BOTTOM_OF(IJK)
                  D_FT = AVG_Z_S(DTi(IJK,ss),DTi(IJKT,ss),K)*C_AT*( T_g(IJKT) - T_g(IJK) )
                  D_FB = AVG_Z_S(DTi(IJKB,ss),DTi(IJK,ss),KM)*C_AB*( T_g(IJK) - T_g(IJKB) )
                ENDIF

	        SumFluxX = SumFluxX + (D_Fe - D_FW)
	        SumFluxY = SumFluxY + (D_FN - D_FS)
	        IF (DO_K) SumFluxZ = SumFluxZ + (D_FT - D_FB)
	    endif

	    divJi(ijk,ss) = divJi(ijk,ss) + SumFluxX
	    divJi(ijk,ss) = divJi(ijk,ss) + SumFluxY
	    IF (DO_K) divJi(ijk,ss) = divJi(ijk,ss) + SumFluxZ

	    if(X_g(IJK, ss) > zero_x_gs) then  ! to amplify central coefficient
	      if(divJi(ijk,ss) < zero) divJi(ijk,ss) = divJi(ijk,ss) / X_g(IJK, ss)
	     ! if(divJi(ijk,ss) > zero) do nothing, term will go to right hand side in Bm coef.
	    else
	      divJi(ijk,ss) = divJi(ijk,ss)  ! or zero
	    endif
	 endDo ! for ss loop
	ELSE  ! not fluid at
	 divJi(ijk,:) = zero
	ENDIF ! fluid at
      ENDDO ! ijk

      RETURN
   END SUBROUTINE MultiComponentDiffusion
