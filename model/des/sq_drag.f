
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Module name: Sq_drag                                             !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Provide a hook for non-spherical drag law implementation.         !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     Module Sq_drag
     contains
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: DRAG_SQP_DiFelice_Ganser                                  !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Gas-solids drag for non-spherical particles in SuperDEM           !
!                                                                             !
!  Reference:                                                                 !
! 1.Xi Gao, Jia Yu, Liqiang Lu, Cheng Li and William A Rogers,                !
!   Development and validation of SuperDEM-CFD coupled model for simulating   !
!	 non-spherical particles hydrodynamics in fluidized beds,                 !
!	 Chemical Engineering Journal, 2021,420: 127654." and	                  !
! 2. G.H. Ganser, A rational approach to drag prediction of spherical         !
!	 and nonspherical particles, Powder Technol. 77 (1993) 143-152.           !
!  This model requires specification of a sphericity and a reference length   !
! (typically the bed diameter)"/>                                             !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DRAG_SQP_DiFelice_Ganser(NP,lDgA,EPg,Mug,ROg,VREL, DPM,&
                lUg, lVg, lWg)
!-----------------------------------------------
! Modules
!-----------------------------------------------
      use error_manager
      USE param1, only:LARGE_NUMBER, ONE, ZERO, HALF
      use discretelement, only: ro_sol
      use discretelement, only: super_q,super_r,super_mn
      USE run
      USE parallel
      USE geometry
      USE indices
      USE constant
      USE fldvar
      USE rxns
      USE physprop
      USE functions
      USE discretelement, only: des_usr_var
      use SQ_PROPERTIES_MOD
      USE run, only: REF_LENGTH_DG

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM
! fluid velocity components:
! o TFM: Averaged from faces to cell center
! o DES: Interpolated to the particle's position
      DOUBLE PRECISION, INTENT(IN) :: lUg, lVg, lWg
! particle number id.
      INTEGER , INTENT(IN) :: NP
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      DOUBLE PRECISION :: RE
! Single sphere drag coefficient
      DOUBLE PRECISION :: C_d,ROPg
!-----------------------------------------------
! Superquadric particle properties
      DOUBLE PRECISION :: ugc(3),ugcc(3),Qc(4), axi(3),m,n
      DOUBLE PRECISION :: cross_sphericity_GANSER,cross_sphericity_HOLZER
      DOUBLE PRECISION :: project_area
      DOUBLE PRECISION :: f_u,f_v,f_w
! sphericity, radius of equivalent volume, diameter of equivalent volume
      DOUBLE PRECISION :: fai,eq_r_v,DPM_V,ex
! Parameters in Ganser non-spherical drag coefficient
      DOUBLE PRECISION :: dvdD,k1,k2
      INTEGER:: shapes
!----------------------------------------------

      ROPg=ROg*EPg

      IF(Mug > 0.0) THEN
         RE = DPM*VREL*ROg*EPg/Mug
      ELSE
         RE = LARGE_NUMBER
      ENDIF

      IF (RE == ZERO) THEN
         lDgA = ZERO
         RETURN
      ENDIF

      Qc(1)= super_q(NP,1)
      Qc(2)= super_q(NP,2)
      Qc(3)= super_q(NP,3)
      Qc(4)= super_q(NP,4)
      axi(1)=super_r(NP,1)
      axi(2)=super_r(NP,2)
      axi(3)=super_r(NP,3)
      m=super_mn(NP,1)
      n=super_mn(NP,2)

! Shape of a particle
     call Sq_shapes(axi,m,n,shapes)
! Superquadric particle sphericity-fai
      call Sq_sphericity(axi,m,n,fai,eq_r_v,shapes)
      dpm_v=2.0*eq_r_v
      dvdD = dpm_v/REF_LENGTH_DG
!superquadric cross sphericity
      ugcc(1)=lug
      ugcc(2)=lvg
      ugcc(3)=lwg
      call Sq_cross_sphericity(ugcc, Qc, axi,m,n, &
          project_area,cross_sphericity_ganser,cross_sphericity_holzer,shapes)
! Difelice-Ganser hybrid drag model
      RE=RE+0.000001d0
      k1=(1.0d0/3.0d0*cross_sphericity_ganser+2.0d0/3.0d0*fai**(-0.5d0))**(-1.0d0)-2.25d0*dvdD
      k2= 10.0d0**(1.8148d0*(-log10(fai))**0.5743d0)
      C_d = (24.D0/k1/(RE)) * (1.0 + 0.1118D0*(RE*k1*k2)**0.6567D0)+&
                     0.4305d0*k2/(1+3305d0/(RE*k1*k2))
      ex=3.70-0.65*exp(-0.5*(1.5-log10(RE+0.00000001D0))**2.0)
      lDgA = 0.75D0*C_d*VREL*ROPg*EPg**(1.0-ex) / DPM_v

      RETURN
      END SUBROUTINE DRAG_SQP_DiFelice_Ganser


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: DRAG_SQP_DiFelice_Holzer_Sommerfeld                       !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Gas-solids drag for non-spherical particles in SuperDEM           !
!                                                                             !
!  Reference:                                                                 !
! 1.Xi Gao, Jia Yu, Liqiang Lu, Cheng Li and William A Rogers,                !
!     Development and validation of SuperDEM-CFD coupled model for simulating !
!	  non-spherical particles hydrodynamics in fluidized beds,                !
!	  Chemical Engineering Journal, 2021,420: 127654.   	                  !
! 2. A. Hölzer, M. Sommerfeld, New simple correlation formula for the drag    !
!	  coefficient of nonspherical particles,                                  !
!	  Powder Technology 184 (2008) 361-365.                                   !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_SQP_DiFelice_Holzer_Sommerfeld(NP,lDgA,EPg,Mug,ROg,VREL, DPM,&
                lUg, lVg, lWg)
!-----------------------------------------------
! Modules
!-----------------------------------------------
      use error_manager
      USE param1, only:LARGE_NUMBER, ONE, ZERO, HALF
      use discretelement, only: ro_sol
      use discretelement, only: super_q,super_r,super_mn
      USE run
      USE parallel
      USE geometry
      USE indices
      USE constant
      USE fldvar
      USE rxns
      USE physprop
      USE functions
      USE discretelement, only: des_usr_var
      use SQ_PROPERTIES_MOD
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM
! fluid velocity components:
! o TFM: Averaged from faces to cell center
! o DES: Interpolated to the particle's position
      DOUBLE PRECISION, INTENT(IN) :: lUg, lVg, lWg
! particle number id.
      INTEGER , INTENT(IN) :: NP
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Reynolds number
      DOUBLE PRECISION :: RE
! Single sphere drag coefficient
      DOUBLE PRECISION :: C_d,ROPg
!-----------------------------------------------
! Superquadric particle properties
      DOUBLE PRECISION :: ugc(3),ugcc(3),Qc(4), axi(3),m,n
      DOUBLE PRECISION :: cross_sphericity_GANSOR,cross_sphericity_HOLZER
      DOUBLE PRECISION :: project_area
      DOUBLE PRECISION :: f_u,f_v,f_w
! sphericity, radius of equivalent volume, diameter of equivalent volume
      DOUBLE PRECISION :: fai,eq_r_v,DPM_V,ex
      INTEGER:: shapes

      ROPg=ROg*EPg
      IF(Mug > 0.0) THEN
         RE = DPM*VREL*ROg*EPg/Mug
      ELSE
         RE = LARGE_NUMBER
      ENDIF

      IF (RE == ZERO) THEN
         lDgA = ZERO
         RETURN
      ENDIF

      Qc(1)= super_q(NP,1)
      Qc(2)= super_q(NP,2)
      Qc(3)= super_q(NP,3)
      Qc(4)= super_q(NP,4)
      axi(1)=super_r(NP,1)
      axi(2)=super_r(NP,2)
      axi(3)=super_r(NP,3)
      m=super_mn(NP,1)
      n=super_mn(NP,2)

! Shape of a particle
     call Sq_shapes(axi,m,n,shapes)
! Superquadric particle sphericity-fai
      call Sq_sphericity(axi,m,n,fai,eq_r_v,shapes)
      dpm_v=2.0*eq_r_v
!superquadric cross sphericity
      ugcc(1)=lug
      ugcc(2)=lvg
      ugcc(3)=lwg
      call Sq_cross_sphericity(ugcc, Qc, axi,m,n, &
          project_area,cross_sphericity_gansor,cross_sphericity_holzer,shapes)
! Holzer and sommerfeld drag, note not 0.4210, its 0.42*10^
      C_d =  8.0/RE/cross_sphericity_holzer**0.5D0 +&
            16.0/RE/fai**0.5D0 +&
             3.0/RE**0.5D0/fai**0.75d0 + &
            0.42*10**(0.4d0*(-log10(fai))**0.2d0)/cross_sphericity_holzer
      ex=3.70-0.65*exp(-0.5*(1.5-log10(RE))**2.0)
      lDgA= 0.5d0*C_d*VREL*ROPg*EPg**(1.0-ex)*project_area/(4.0/3.0*PI*eq_r_v**3.0)
      RETURN

      END SUBROUTINE DRAG_SQP_DiFelice_Holzer_Sommerfeld


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: DRAG_SQP_GIDASPOW_ganser                                  !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Gas-solids drag for non-spherical particles in SuperDEM           !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DRAG_SQP_GIDASPOW_ganser(NP,lDgA,EPg,Mug,ROg,VREL, DPM,&
                lUg, lVg, lWg)
!-----------------------------------------------
! Modules
!-----------------------------------------------
      use error_manager
      USE param1, only:LARGE_NUMBER, ONE, ZERO, HALF
      use discretelement, only: ro_sol
      use discretelement, only: super_q,super_r,super_mn
      USE run
      USE parallel
      USE geometry
      USE indices
      USE constant
      USE fldvar
      USE rxns
      USE physprop
      USE functions
      USE discretelement, only: des_usr_var
      use SQ_PROPERTIES_MOD
      USE run, only: REF_LENGTH_DG
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM
! fluid velocity components:
! o TFM: Averaged from faces to cell center
! o DES: Interpolated to the particle's position
      DOUBLE PRECISION, INTENT(IN) :: lUg, lVg, lWg
! particle number id.
      INTEGER , INTENT(IN) :: NP
! Local variables
!-----------------------------------------------
! Reynolds number
      DOUBLE PRECISION :: RE
! Single sphere drag coefficient
      DOUBLE PRECISION :: C_d,ROPg
!-----------------------------------------------
! Superquadric particle properties
      DOUBLE PRECISION :: ugc(3),ugcc(3),Qc(4), axi(3),m,n
      DOUBLE PRECISION :: cross_sphericity_GANSER,cross_sphericity_HOLZER
      DOUBLE PRECISION :: project_area
      DOUBLE PRECISION :: f_u,f_v,f_w
! sphericity, radius of equivalent volume, diameter of equivalent volume
      DOUBLE PRECISION :: fai,eq_r_v,DPM_V
      DOUBLE PRECISION :: bed_diameter,dvdD,k1,k2
      INTEGER:: shapes
!----------------------------------------------

      IF(Mug > 0.0) THEN
         RE = DPM*VREL*ROg*EPg/Mug
      ELSE
         RE = LARGE_NUMBER
      ENDIF

      IF (RE == ZERO) THEN
         lDgA = ZERO
         RETURN
      ENDIF

      Qc(1)= super_q(NP,1)
      Qc(2)= super_q(NP,2)
      Qc(3)= super_q(NP,3)
      Qc(4)= super_q(NP,4)
      axi(1)=super_r(NP,1)
      axi(2)=super_r(NP,2)
      axi(3)=super_r(NP,3)
      m=super_mn(NP,1)
      n=super_mn(NP,2)

! Shape of a particle
     call Sq_shapes(axi,m,n,shapes)
! Superquadric sphericity
      call Sq_sphericity(axi,m,n,fai,eq_r_v,shapes)
      dpm_v=2.0*eq_r_v
      dvdD = dpm_v/REF_LENGTH_DG
!Superquadric cross sphericity
      ugcc(1)=lug
      ugcc(2)=lvg
      ugcc(3)=lwg
      call Sq_cross_sphericity(ugcc, Qc, axi,m,n, &
            project_area,cross_sphericity_ganser,cross_sphericity_holzer,shapes)
! Ganser drag for non-spherical
      k1=(1.0d0/3.0d0*cross_sphericity_ganser+2.0d0/3.0d0*fai**(-0.5d0))**(-1.0d0)-2.25d0*dvdD
      k2= 10.0d0**(1.8148d0*(-log10(fai))**0.5743d0)
! Dense phase
         IF(EPg <= 0.80D0) THEN
            lDgA = 150D0*(ONE-EPg)*Mug / (EPg*DPM_v**2*fai**2) + &
                   1.75D0*ROg*VREL/DPM_v/fai
         ELSE
! Dilute phase - EP_g >= 0.8
            IF(RE <= 1000D0)THEN
! this could be replaced with the function C_DS_SN
               C_d = (24.D0/k1/(RE+0.00000001D0)) * (1.0 + 0.1118D0*(RE*k1*k2)**0.6567D0)+&
                     0.4305d0*k2/(1+3305d0/(RE*k1*k2))
            ELSE
               C_d = 0.44D0
            ENDIF
        lDgA= 0.5d0*C_d*VREL*ROg*EPg**(-2.65d0)*project_area/(4.0/3.0*PI*eq_r_v**3.0)
        ENDIF

      RETURN

      END SUBROUTINE DRAG_SQP_GIDASPOW_Ganser

      END Module sq_drag
