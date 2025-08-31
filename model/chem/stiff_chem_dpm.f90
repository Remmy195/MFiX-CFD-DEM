#include "error.inc"
module stiff_chem_dpm

  public :: stiff_chem_rrates_dpm
  public :: mapMFIXtoODE_dpm
  public :: mapODEtoMFIX_dpm
  public :: reportODEVar_dpm

contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Module name: stiff_chem_rrates_dpm                                !
  !  Author: J. Musser                             Date: May 05, 2023  !
  !                                                                    !
  !  Purpose: Calculate reaction rates for various reactions for DES.  !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine stiff_chem_rrates_dpm(a_NEq, a_Y, a_norm, a_YDOT)

    use fldvar,         only: T_g, X_g, rop_g

    use des_thermo,     only: T_p => des_t_s
    use des_rxns,       only: X_p => des_x_s
    use discretelement, only: m_p => pmass

    use physprop, only: mmax, nmax

    use discretelement, only: pinc, pijk
    use derived_types,  only: pic
    use des_rxns,       only: des_min_pmass

    use mfix_pic, only: mppic, des_stat_wt

    use des_rxns, only: no_of_des_rxns

    use des_rxns, only: des_reaction

    use rxns, only: reaction
    use rxns, only: no_of_rxns

    use run, only: units

    use param,  only: dimension_n_g, dimension_n_s
    use param1, only: zero, small_number

    use parse, only: ARRHENIUS_RRATES_FLUID, ARRHENIUS_RRATES_DES

    use toleranc, only: chem_min_species_fluid, chem_min_species_solid, des_min_pmass_frac

    use calc_h_mod, only: calc_h
    use toleranc,   only: compare
    use vtk

    use cutcell

    implicit none

    ! Passed Variables: Dummy argument format required by ODEPACK.
    !---------------------------------------------------------------------//
    ! (1) Number of ODEs to be solve
    ! (2) Fluid cell index
    ! (3) Number of dpm solids
    ! (3) Start of solve/don't solve flags for solids
    integer         , intent(in   ) :: a_NEq(:)
    ! Array of dependent variable initial values.
    double precision, intent(in   ) :: a_Y(a_NEq(1))
    ! Array of particle mass before calling stiff solver,
    ! used to normalize particle mass and species mass in particles
    double precision, intent(in   ) :: a_norm(:)
    ! Rate of change of dependent variables.
    double precision, intent(  out) :: a_YDOT(a_NEq(1))

    ! Local variables
    !---------------------------------------------------------------------//
    integer :: ijk  ! fluid cell index
    integer :: h    ! Reaction loop counter
    integer :: p    ! Local loop index for particle
    integer :: pid  ! index of particle in global des arrays

    ! > Cumulative over all reactions.
    !...................................................................
    ! Rate of formation (+) or consumption (-) of fluid species
    double precision :: cRg(dimension_n_g)
    ! Rate of formation (+) or consumption (-) of fluid species from homogeneous reactions
    double precision :: cRg_homo(dimension_n_g)
    ! Rate of formation (+) or consumption (-) of dpm species.
    double precision :: cRp(dimension_n_s)
    ! Heat of reaction for a reaction.
    double precision :: cHoRg, cHoRp, CHoRg_homo

    ! User-defined reaction rates returned from DES_USR_RATES
    double precision :: des_rates(no_of_des_rxns)
    ! User-defined reaction rates returned from USR_RATES
    double precision :: rates(no_of_rxns)
    ! Weight of dpm particles (1 for DEM)
    double precision :: stat_wt
    ! Scale factor for energy units
    double precision :: lto_SI

! Keep a copy of last node value at the end of the p-loop
    integer :: last_node
! index of partilces which have reactions
    integer, dimension (:), allocatable :: p_rxn
    integer :: ss

    lto_SI = merge(4.183925d3, 1.0d0, UNITS == 'SI')

    ijk = a_NEq(2)
    a_YDOT = 0.0d0

    ! Map the current ODE independent variables to MFIX variables.
    call mapODEtoMFIX_dpm(size(a_NEq), a_NEq, size(a_Y), a_norm, a_Y)

! last_node is initially set to the number of nodes for the gas phase
! (density + temperature + number of gas species)
! Each particle p will start incrementing the node value from last_node
! Each particle p will add a node for temperature.
! If reactions occur for particle p, (1 + nmax_s(m)) nodes will be added
! (mass + number of solids species)for this particle.
    last_node = 2 + nmax(0)

! calculate reaction rates of homogeneous reactions
    ! Calculate the changes of gas phase from homogeneous reactions
    cRg_homo = zero; CHoRg_homo = zero
    RATES  = ZERO

    ! Calculate user defined reaction rates.
    IF(ARRHENIUS_RRATES_FLUID) THEN
       CALL CALC_ARRHENIUS_RRATES(IJK, RATES)
    ELSE
       CALL USR_RATES(IJK, RATES)
    ENDIF

    ! Loop over reactions.
    rxn_lp_homo: do h = 1, no_of_rxns
      ! Skip empty reactions
      if (Reaction(H)%nSpecies == 0) cycle rxn_lp_homo
      rxn_bk_homo: block
      ! Rate of formation (+) or consumption (-) of a species
      double precision :: rRate
      ! Rate of formation (+) or consumption (-) of fluid species from jomogeneour reactions
      double precision :: rRg_homo(dimension_n_g)
      ! Heat of reaction for a reaction.
      double precision :: rHoRg_homo
      integer :: lN   ! Local reaction species index/loop counter
      integer :: m    ! Global Phase index loop counter
      integer :: n    ! Global species inde

      ! Initialize local loop arrays
      rRg_homo = zero; rHoRg_homo = zero

      ! Check if there are enough reactants for this reaction
      do lN = 1,Reaction(H)%nSpecies
        ! Global species index.
        N = Reaction(H)%Species(lN)%sMap

        if (Reaction(H)%Species(lN)%MWxStoich .lt. zero) then
          if(X_g(ijk,n) .le. chem_min_species_fluid) then
            rates(h) = zero
            cycle rxn_lp_homo
          endif
       endif
      enddo

      ! Calculate the rate of formation/consumption for each species.
      !-------------------------------------------------------------//
      do lN = 1, Reaction(H)%nSpecies

        ! Global phase index. M should be 0 here for dpm.
        M = Reaction(H)%Species(lN)%pMap
        ! Global species index.
        N = Reaction(H)%Species(lN)%sMap
        ! Convert rate from moles to kg (g for csg units)
        rRate = rates(h) * Reaction(H)%Species(lN)%MWxStoich
        ! Consumption/Formation of gas phase species.
        rRg_homo(n) = rRg_homo(n) + rRate
      enddo ! Loop of species
      ! Copy the single reaction rates of formation/consumption to
      ! the total (accumulative) rates of formation/consumption
      ! arrays. This ensures that any reactions without sufficient
      ! reactants are not included.
      cRg_homo(:) = cRg_homo(:) + rRg_homo(:)

      ! Calculate and store the heat of reaction.
      !-------------------------------------------------------------//
      ! Automated heat of reaction calculations
      if (Reaction(H)%calc_dH) then
         ! Loop over reaction species.
         do lN = 1, Reaction(H)%nSpecies

           ! Global phase index. m should be 0 here.
           m = Reaction(h)%Species(lN)%pMap
           ! Global species index.
           n = Reaction(h)%Species(lN)%sMap
           ! Rate of formation/consumption for species N
           rRate = rates(h) * Reaction(h)%Species(lN)%MWxStoich

           ! Gas phase enthalpy change
           rHoRg_homo = rHoRg_homo + rRate * calc_H(T_g(ijk),0,n)
         enddo

         ! Convert the heat of reaction to the appropriate units
         ! (if SI), and store in the global array.
         cHORg_homo = cHORg_homo + rHoRg_homo*lto_SI

      else ! User-defined heat of reaction.
        cHoRg_homo = cHoRg_homo + Reaction(H)%HoR(0) * RATES(H)
      endif
      end block rxn_bk_homo
    enddo rxn_lp_homo ! loop over homogeneous reactions.

    ! Save rates so they can be written in vtu files and/or monitors
    call save_fluid_rrates_cell_data(ijk, rates)

    fill_ydot_homo: block
    !---------------------------------------------------------------//
    ! Map the reaction rates into the ODE solver

      use fldvar, only: rop_g
      use physprop,   only: C_pg
      integer :: node, n ! YDOT and species loop counters

      node = 1

      ! Bulk density:  rop_g
      a_YDOT(node) = a_YDOT(node) + sum(cRg_homo);       node = node + 1

      ! Temperature: T_g
      cHoRg_homo = cHoRg_homo/(rop_g(ijk)*c_pg(ijk))
      a_YDOT(node) = a_YDOT(node) - cHoRg_homo;              node = node + 1

      ! Bulk species density: rop_g * X_g
      do n=1,nmax(0)
        a_YDOT(node) = a_YDOT(node) + cRg_homo(n);       node = node + 1
      enddo

    end block fill_ydot_homo

! Return if there are no particles in the cell
    if(a_NEQ(3)==0) return
! Get the list of particles having reactions
    ss = COUNT(a_NEQ == 1)
    if(ss .ne. 0) then
      allocate(p_rxn(ss))
      p_rxn = pack([(ss-3, ss=1, size(a_NEQ))], a_NEQ==1)
    else
      return
    endif

    do ss=1, size(p_rxn)
      if(p_rxn(ss) .le. 0) cycle
      p = p_rxn(ss)

      ! clear cumulative arrays for each particle
      cRg = zero; cHORg = zero
      cRp = zero; cHORp = zero

      pid = pic(ijk)%p(p)

      ! Parcel mass scaling multiplier
      if (mppic) then
         stat_wt = des_stat_wt(pid)
      else
         stat_wt = 1.0d0
      endif

      ! Calculate user defined reaction rates.
      des_rates(:) = zero
      if(arrhenius_rrates_des) then
         call calc_arrhenius_rrates_des(pid, pijk(pid,5), ijk, des_rates)
      else
         call usr_rates_des(pid, pijk(pid,5), ijk, des_rates)
      endif

      ! Loop over reactions.
      rxn_lp: do h = 1, no_of_des_rxns

        ! Skip empty reactions
        if (des_Reaction(H)%nSpecies == 0) cycle rxn_lp

        ! scale des_rates by mass to facilitate calculations involving
        ! very small/slow rates that would otherwise be zeroed; this is
        ! exacerbated by working in SI units where mass is in unit kg and
        ! unit substances is kmol with very small/light particles
        if (compare(des_rates(h)/m_p(pid),zero)) cycle rxn_lp
        rxn_bk: block
        !-------------------------------------------------------------//
        ! Rate of formation (+) or consumption (-) of a species
        double precision :: rRate
        ! Rate of formation (+) or consumption (-) of fluid species
        double precision :: rRg(dimension_n_g)
        ! Rate of formation (+) or consumption (-) of dpm species.
        double precision :: rRp(dimension_n_s)
        ! Heat of reaction for a reaction.
        double precision :: rHoRg, rHoRp
        ! Rate of interphase enthalpy transfer due to mass transfer.
        double precision :: rRxH

        integer :: lN   ! Local reaction species index/loop counter
        integer :: m    ! Global Phase index loop counter
        integer :: n    ! Global species index
        integer :: mXfr ! Global phase index for mass transfer

        ! Initialize local loop arrays
        rRg = zero; rHoRg = zero
        rRp = zero; rHoRp = zero
        rRxH = zero

        ! Check if there are enough reactants for this reaction
        do lN = 1, des_reaction(H)%nSpecies
          ! Global phase index.
          m = des_reaction(H)%Species(lN)%pMap
          ! Global species index.
          n = des_reaction(H)%Species(lN)%sMap

          if(des_reaction(H)%Species(lN)%MWxStoich .lt. zero) then
            if(m == 0) then
              if (X_g(ijk,n) .le. chem_min_species_fluid) then
                des_rates(h) = zero
                cycle rxn_lp
              endif
            else
              if (X_p(pid,n) .le. chem_min_species_solid) then
                des_rates(h) = zero
                cycle rxn_lp
              endif
              ! check if the particle mass is big enough to consider reactions.
              if(m_p(pid) .le. des_min_pmass(m)*des_min_pmass_frac) then
                des_rates(h) = zero
                cycle rxn_lp
              endif
            endif
          endif
        enddo

        ! Calculate the rate of formation/consumption for each species.
        !-------------------------------------------------------------//
        do lN = 1, des_reaction(H)%nSpecies

          ! Global phase index.
          m = des_reaction(H)%Species(lN)%pMap
          ! Global species index.
          n = des_reaction(H)%Species(lN)%sMap

          ! Index for interphase mass transfer. For a gas/solid reaction,
          ! the index is stored with the gas phase.
          mXfr = des_reaction(H)%Species(lN)%mXfr

          rRate = des_rates(h) * des_reaction(H)%Species(lN)%MWxStoich

          ! Gas Phase:
          if (m == 0) then
            ! Consumption/Formation of gas phase species.
            rRg(n) = rRg(n) + rRate
            ! Enthalpy transfer associated with mass transfer. (gas/solid)
            if (m /= mXfr) rRxH = rRxH + rRate*calc_H(T_g(ijk),0,n)
          else ! discrete particle model
            ! Formation/consumption of solids phase species.
            rRp(n) = rRp(n) + rRate
          endif

        enddo ! Loop of species


        ! Copy the single reaction rates of formation/consumption to
        ! the total (accumulative) rates of formation/consumption
        ! arrays. This ensures that any reactions without sufficient
        ! reactants are not included.
        cRg(:) = cRg(:) + rRg(:)*stat_wt
        cRp(:) = cRp(:) + rRp(:)

        ! Calculate and store the heat of reaction.
        !-------------------------------------------------------------//
        rHoRg = zero
        rHoRp = zero

        ! Automated heat of reaction calculations
        if (des_reaction(H)%calc_DH) then
          ! Loop over reaction species.
          do lN = 1, des_reaction(H)%nSpecies
            ! Global phase index.
            M = des_reaction(H)%Species(lN)%pMap
            ! Global species index.
            N = des_reaction(H)%Species(lN)%sMap

            ! Rate of formation/consumption for species N
            rRate = des_rates(h) * DES_Reaction(H)%Species(lN)%MWxStoich

            if(m == 0) then
              ! Gas phase enthalpy change from energy equation derivation.
              rHORg = rHORg + calc_H(T_g(ijk),0,n) * rRate
            else
              ! Particle enthalpy change from energy equation derivation.
              rHORp = rHORp + calc_H(T_p(pid),m,n) * rRate
            endif
          enddo

          ! Apply enthalpy transfer associated with mass transfer to get the
          ! complete heat of reaction for Reaction H.
          rHORg = rHORg - rRxH
          rHORp = rHORp + rRxH

          ! Convert the heat of reaction to the appropriate units (if SI),
          ! and store in the global array.
          cHORg = cHORg + lto_SI*rHORg * stat_wt
          cHORp = cHORp + lto_SI*rHORp

        else

          ! User-defined heat of reaction.
          rHoRg = des_rates(h) * des_reaction(H)%HoR( 0)
          rHoRp = des_rates(h) * des_reaction(H)%HoR(pijk(pid,5))

          cHoRg = cHoRg + rHoRg * stat_wt
          cHORp = cHORp + rHoRp

        endif

        end block rxn_bk

      enddo rxn_lp ! Loop over reactions.

      call save_des_rrates_data(pid, ijk, des_rates)


      fill_ydot: block
      !---------------------------------------------------------------//
      ! Map the reaction rates into the ODE solver

        use fldvar, only: rop_g

        use physprop,   only: C_pg
        use des_thermo, only: C_ps => des_C_ps

        use geometry, only: vol

        integer :: node, n ! YDOT and species loop counters
        double precision :: mpCps, ivol

        ivol = 1.0d0 / vol(ijk)

        node = 1

        ! Bulk density:  rop_g
        a_YDOT(node) = a_YDOT(node) + sum(cRg)*ivol;       node = node + 1

        ! Temperature: T_g
        cHORg = (cHORg/(rop_g(ijk)*c_pg(ijk)))*ivol
        a_YDOT(node) = a_YDOT(node)  - cHoRg;              node = node + 1

        ! Bulk species density: rop_g * X_g
        do n=1,nmax(0)
          a_YDOT(node) = a_YDOT(node) + cRg(n)*ivol;       node = node + 1
        enddo


! Reset the node so it starts at the first slot for particle p
! Temperature node is always added
        node = last_node + 1

        ! Solids temperature
        mpCps = m_p(pid)*C_ps(pid)

        if (mpCps > small_number) then
           a_YDOT(node) = -cHoRp/mpCps
        else
           a_YDOT(node) = zero
        endif

        node = node + 1
          ! particle mass: mp
          a_YDOT(node) = sum(cRp)/a_norm(p);                            node = node + 1

          ! species mass mp * X_pn
          do n=1, nmax(pijk(pid,5))
            a_YDOT(node) = cRp(n)/a_norm(p);                          node = node + 1
          enddo

! Save a copy of the last node after filling YDOT for the current particle
! The next particle will start from last_node + 1 after its contribution
! to the gas phase have been added
        last_node = node - 1


      end block fill_ydot

    enddo ! loop over particle p


    !---------------------------------------------------------------//

  end subroutine stiff_chem_rrates_dpm


  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Module name: mapMFIXtoODE_dpm                                     !
  !  Author: J. Musser                             Date: May 05, 2023  !
  !                                                                    !
  !  Purpose: This routine maps MFIX variables into the ODE array.     !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine mapMFIXtoODE_dpm (a_NEq_d, a_NEq, a_ODE_Vars_d, a_ODE_norm, a_ODE_Vars)

    ! Global Variables:
    !-----------------------------------------------------------------//
    use fldvar, only: ROP_g, T_g, X_g

    use physprop, only: nmax

    use discretelement, only: pinc, pijk
    use derived_types,  only: pic

    use des_rxns,       only: X_p => des_x_s
    use des_thermo,     only: T_p => des_T_s
    use discretelement, only: m_p => pmass

    implicit none

    ! Passed Variables:
    !-----------------------------------------------------------------//
    ! (1) Number of ODEs to be solve
    ! (2) Fluid cell index
    integer         , intent(in   ) :: a_NEq_d
    integer         , intent(in   ) :: a_NEq(a_NEq_d)

    ! Array of dependent variable initial values.
    integer         , intent(in   ) :: a_ODE_Vars_d
    ! Values used to normalize the particle mass and species mass in particles
    double precision, intent(in   ) :: a_ODE_norm(:)
    double precision, intent(  out) :: a_ODE_Vars(a_ODE_Vars_d)

    ! Local Variables:
    !-----------------------------------------------------------------//
    ! Local indices
    integer :: node, ijk, n, p, pid
    !...................................................................

    ! Initialize.
    ijk = a_NEq(2)
    a_ODE_Vars = 0.0d0

    node = 1

    ! Gas phase bulk density.
    a_ODE_Vars(node) = rop_g(ijk);                       node = node + 1

    ! Gas phase temperature.
    a_ODE_Vars(node) = T_g(ijk);                         node = node + 1

    ! Gas phase species mass.
    do n=1,nmax(0)
      a_ODE_Vars(node) = rop_g(ijk)*X_g(ijk,n);          node = node + 1
    enddo

    do p=1,a_NEq(3)

      pid = pic(ijk)%p(p)

      if (a_NEq(3 + p) == 1) then

        ! Solids temperature.
        a_ODE_Vars(node) = T_p(pid);                                   node = node + 1

        ! Total solids mass
        a_ODE_Vars(node) = m_p(pid)/a_ODE_norm(p);                     node = node + 1

        ! Species mass
        do n=1,nmax(pijk(pid,5))
          a_ODE_Vars(node) = m_p(pid)*X_p(pid,n)/a_ODE_norm(p);        node = node + 1
        enddo

      endif
    enddo

  end subroutine mapMFIXtoODE_dpm


  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Module name: mapODEtoMFIX_dpm                                     !
  !  Author: J. Musser                             Date: May 05, 2023  !
  !                                                                    !
  !  Purpose: This is a driver routine for mapping variables stored in !
  !  the ODE array back to MFIX field variables.                       !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine mapODEtoMFIX_dpm(a_NEq_d, a_NEq, a_ODE_Vars_d, a_ODE_norm, a_ODE_Vars)

    ! Global Variables:
    !-----------------------------------------------------------------//
    use fldvar, only: ROP_g, RO_g, T_g, X_g, P_g, EP_g

    use des_rxns,       only: X_p  => des_x_s
    use des_thermo,     only: T_p  => des_T_s
    use discretelement, only: m_p  => pmass

    use discretelement, only: rho_p => ro_sol
    use discretelement, only: rad_p  => des_radius
    use discretelement, only: vol_p  => pvol
    use discretelement, only: omoi, omoi3
    use discretelement, only: CGDEM, DES_CGP_RPR, DES_CGP_STW
    use discretelement, only: super_r, super_mn,SuperDEM
    use discretelement, only: gsp_explicit, gsp_implicit, glueInertia
    use discretelement, only: glueVolume, glueDiameter, glueMass
    use discretelement, only: gid_list, particle_state, iglobal_id
    use discretelement, only: gp_sa, gp_neighsa, gp_neighsid, sc2gpc_vec
    use discretelement, only: max_pip
    use error_manager
    use sq_properties_mod
    use physprop, only: mmax, nmax

    ! Gas constant (cal/g.K)
    use constant, only: GAS_CONST

    use constant, only: pi

    ! MPI
    use compar, only: PE_IO, myPE, numPEs
    use mpi_utility, only: bcast, global_all_sum

    ! Mixture and species molecular weights
    use physprop, only: MW_MIX_g, MW_g

    use run, only: solve_ros

    use discretelement, only: pinc, pijk
    use derived_types,  only: pic

    ! Variable solids density calculation.
    use eos, only: eoss

    ! Global Parameters:
    !-----------------------------------------------------------------//
    use param1,   only : zero, small_number
    use param1,   only :  one, large_number

    use scales, only: scale_pressure

    implicit none

    ! Passed Variables:
    !-----------------------------------------------------------------//
    ! (1) Number of ODEs to be solve
    ! (2) Fluid cell index
    integer         , intent(in   ) :: a_NEq_d
    integer         , intent(in   ) :: a_NEq(a_NEq_d)

    ! Array of dependent variable initial values.
    integer         , intent(in   ) :: a_ODE_Vars_d
    ! Values used to normalize the particle mass and species mass in particles
    double precision, intent(in   ) :: a_ODE_norm(:)
    double precision, intent(in   ) :: a_ODE_Vars(a_ODE_Vars_d)

    ! Local Variables:
    !-----------------------------------------------------------------//
    double precision, parameter :: pix4o3 = pi*(4.0d0/3.0d0)
    double precision, parameter :: third = (1.0d0/3.0d0)

    double precision :: im_p, mp_normalize
    double precision :: axi(3),mm,nn,IXX, IYY, IZZ, SVOLUME, scale_factor
    integer :: node, ijk, n, p, pid
    ! GSP Variables
    double precision :: rad_p_old
    !...................................................................

    ! Initialize.
    ijk = a_NEq(2)
    node = 1

    ! Gas phase density.
    rop_g(ijk) = a_ODE_Vars(node);                       node = node + 1

    ! Gas phase temperature.
    T_g(ijk) = a_ODE_Vars(node);                         node = node + 1

    ! Gas phase species mass fractions.
    do n=1,nmax(0)
      X_g(ijk,n) = a_ODE_Vars(node)/rop_g(ijk);          node = node + 1
      X_g(ijk,n) = max(0.0, X_g(ijk,n))
    enddo

    do p=1, a_NEq(3)

      pid = pic(ijk)%p(p)

      if (a_NEq(p + 3) == 1) then
        ! Solids temperature.
        T_p(pid) = a_ODE_Vars(node);                       node = node + 1
        ! Total solids mass
        mp_normalize = max(a_ODE_Vars(node), zero);          node = node + 1
        m_p(pid) = mp_normalize * a_ODE_norm(p)
        if (mp_normalize > small_number) then
          im_p = one / mp_normalize
        else
          im_p = zero
        endif

        ! Species mass
        do n=1,nmax(pijk(pid,5))
          X_p(pid,n) = im_p*a_ODE_Vars(node);            node = node + 1
          X_p(pid,n) = max(0.0, X_p(pid,n))
        enddo

        if(superDEM) then
          axi(1)=super_r(pid,1)
          axi(2)=super_r(pid,2)
          axi(3)=super_r(pid,3)
          mm=super_mn(pid,1)
          nn=super_mn(pid,2)

          if(solve_ros(pijk(pid,5))) then     ! variable density, constant size
            rho_p(pid)= m_p(pid)/vol_p(pid)
            call SQ_INERTIA(axi,mm,nn,IXX,IYY,IZZ)
            omoi3(pid,1)=1.0/(IXX*rho_p(pid))
            omoi3(pid,2)=1.0/(IYY*rho_p(pid))
            omoi3(pid,3)=1.0/(IZZ*rho_p(pid))
          else                               ! constant density, variable size
! Assume the particle shape remains the same (same aspect ratio, same m,n values)
! Then we can scale a,b,c by multiplying them by (new mass/ old mass)^(1/3)
! Here the old mass is still PVOL(NN)*RO_SOL(NN) because the chemical reaction
! has only changed PMASS at this point.
            scale_factor = (m_p(pid)/(vol_p(pid)*rho_p(pid)))**third
            axi(:)       = scale_factor * axi(:)

            CALL SQ_VOLUME (axi,mm,nn,svolume)
            vol_p(pid)=svolume
            call SQ_INERTIA(axi,mm,nn,IXX,IYY,IZZ)
            omoi3(pid,1)=1.0/(IXX*rho_p(pid))
            omoi3(pid,2)=1.0/(IYY*rho_p(pid))
            omoi3(pid,3)=1.0/(IZZ*rho_p(pid))

            super_r(pid,1)  = axi(1)
            super_r(pid,2)  = axi(2)
            super_r(pid,3)  = axi(3)
            rad_p(pid) = scale_factor * rad_p(pid)
          endif
        elseif(gsp_explicit) then
          if(particle_state(pid) /= 1) cycle
          if(solve_ros(pijk(pid,5))) then ! variable density, constant size
            rho_p(pid)= m_p(pid)/vol_p(pid)
          else                           ! constant density, variable size
            ! No operation for now, scale everything based total mass change later
            continue
          endif ! end of constant density
          ! for gsp model, omoi is not used anywhere
          omoi(pid) = 2.5d0/(m_p(pid)*rad_p(pid)**2) ! update one over moi
        else
          if(solve_ros(pijk(pid,5))) then
            rho_p(pid)= m_p(pid)/vol_p(pid)
          else
            if(gsp_implicit) then
                write(err_msg, 1099)
                call log_error()
            endif
 1099 FORMAT('FATAL - variable size chemical reaction is not supported in gsp implicit model.')
            rad_p(pid) = ( m_p(pid) / (Pix4o3*rho_p(pid)) )**third
            vol_p(pid) = m_p(pid) / rho_p(pid)
            ! update real particle size
            if(CGDEM) DES_CGP_RPR(pid) = rad_p(pid) / DES_CGP_STW(pid)**third
          endif
          omoi(pid) = 2.5d0/(m_p(pid)*rad_p(pid)**2)
        endif ! end check for spherical particles
      endif
    enddo

    ! Unlike the tfm implementation, we do not update volume fraction
    ! because that would require a deposition step.
    ! ep_g != (one - sum(ep_s))

    if (ep_g(ijk) > small_number) then

      ! Gas phase bulk density is updated within the stiff solver (a_Var(1)).
      ! Now that the gas phase volume fraction is updated, the gas phase
      ! density can be backed out. RO_g * EP_g = ROP_g
      ro_g(ijk) = rop_g(ijk) / ep_g(ijk)

    else

      ! This case shouldn't happen, however 'LARGE_NUMBER' is used to aid
      ! in tracking errors should this somehow become and issue.
      ro_g(ijk) = large_number

    endif

    ! Calculate the mixture molecular weight.
    mw_mix_g(ijk) = sum(X_g(ijk,1:nmax(0))/MW_g(1:nmax(0)))
    mw_mix_g(ijk) = one/mw_mix_g(ijk)

    ! Calculate the gas phase pressure.
    P_g(ijk) = scale_pressure((ro_g(ijk)*gas_const*T_g(ijk))/mw_mix_g(ijk))

    return
  end subroutine mapODEtoMFIX_dpm


  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Module name: reportODEVar_dpm                                      !
  !  Author: J. Musser                             Date: May 05, 2023  !
  !                                                                    !
  !  Purpose: Checks the values in the ODE array. Currently not used.  !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine reportODEVar_dpm (a_NEq_d, a_NEq, a_ODE_Vars_d, a_norm, a_ODE_Vars)

    use run, only: species_eq

    use discretelement, only: pijk
    use derived_types,  only: pic

    use rxns, only: species_alias_g, species_alias_s
    use physprop, only: nmax
    use param1,   only : small_number

    implicit none

    ! Passed Variables:
    !-----------------------------------------------------------------//
    ! (1) number of ODEs to be solve (required by ODEPACK)
    ! (2) fluid cell index
    ! (3) number of solids
    ! (4) start of solids flags
    integer         , intent(in   ) :: a_NEq_d
    integer         , intent(in   ) :: a_NEq(a_NEq_d)

    ! Array of dependent variable initial values.
    integer         , intent(in   ) :: a_ODE_Vars_d
    double precision, intent(in   ) :: a_ODE_Vars(a_ODE_Vars_d)
    double precision, intent(in   ) :: a_norm(:)
    double precision :: ep_rho, inv_mp, mp, Xs_mp, Xs
    double precision :: sum_Xg, sum_Xs, sum_Xs_mp

    integer :: node, m, n, p, pid

    node = 1

    ! Gas phase bulk density.
    ep_rho = a_ODE_Vars(node)
    node = node + 1
    write(*,"(I3,2x,A24,F12.4)") node, "ROP_g:", ep_rho

    ! Skip checks on temperature. The ODE solver has no means to
    ! enforce these bounds. Therefore, like the energy equations
    ! linear solvers, we 'clip' temperature when mapping back to
    ! global field arrays.
    write(*,"(I3,2x,A24,F12.4)") node, "T_g:", a_ODE_Vars(node)
    node = node + 1

    ! Gas phase species mass fractions.
    sum_Xg = 0.
    do n=1,nmax(0)
      sum_Xg = sum_Xg + a_ODE_Vars(node)/ep_rho
      write(*,"(I3,2x,A24,3x,E18.8)") node, "X_g("//trim(species_alias_g(n))//"):", a_ODE_Vars(node)/ep_rho
      node = node + 1
    enddo


    do p=1, a_NEq(3)

      pid = pic(a_NEq(2))%p(p)
      m = pijk(pid,5)

      ! Particle temperature
      write(*,"(I3,2x,A24,F12.4)") node, "T_p:", a_ODE_Vars(node)
      node = node + 1

      if (a_NEq(3 + p) == 1) then

        ! Particle mass
        mp = a_ODE_Vars(Node)*a_norm(p)
        node = node + 1
        write(*,"(I3,2x,A24,3x,E18.8)") node, "m_p:", mp

        if (mp > small_number) then
           inv_mp = 1.0d0 / mp
        else
           inv_mp = 0.0d0
        endif

        ! Solids phase species mass fractions.
        sum_Xs = 0.
        sum_Xs_mp = 0.
        do n=1,nmax(m)

          Xs_mp = a_ODE_Vars(node)*a_norm(p)
          node = node + 1

          Xs = Xs_mp*inv_mp
          write(*,"(I3,2x,A24,3x,E18.8)") node, "X_s("//trim(species_alias_s(m,n))//"):", Xs

          sum_Xs_mp = sum_Xs_mp + Xs_mp
          sum_Xs    = sum_Xs + Xs

        enddo

      endif
    enddo

    write(*,"(2x,A16,F12.4)") "   sum Xg:", sum_Xg
    write(*,"(2x,A16,F12.4)") "   sum Xs:", sum_Xs
    write(*,"(2x,A16,F12.4)") "sum Xs_mp:", sum_Xs_mp

  end subroutine reportODEVar_dpm


end module stiff_chem_dpm
