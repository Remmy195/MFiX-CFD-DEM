module stiff_chem_tfm

  public :: stiff_chem_rrates_tfm
  public :: mapMFIXtoODE_tfm
  public :: mapODEtoMFIX_tfm
  public :: reportODEVar_tfm

contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Subroutine: stiff_chem_rrates_tfm                                 !
  !  Author: J. Musser                                Date: 12-Feb-13  !
  !                                                                    !
  !  Purpose: Calculate reaction rates for chemical reactions in ijk.  !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine stiff_chem_rrates_tfm (a_NEq, a_Y, a_norm, a_YDOT)

    ! Global Parameters:
    !-----------------------------------------------------------------//
    ! Maximum number of solids phases
    use param, only: dimension_m, dimension_n_g, dimension_n_s
    use param1, only: zero

    ! Global Variables:
    !-----------------------------------------------------------------//
    ! Number of solids phases
    use physprop, only: mmax, nmax

    use fldvar, only: T_g, X_g
    use fldvar, only: T_s, X_s

    use rxns, only: reaction
    use rxns, only: no_of_rxns

    use run, only: units

    use toleranc, only: chem_min_species_fluid, chem_min_species_solid

    use calc_h_mod, only: calc_h

    ! External Function for comparing two numbers.
    use toleranc, only: compare
    use toleranc, only: zero_ep_s

    ! External routine for saving reaction rates in vtu files
    use vtk, only: save_fluid_rrates_cell_data

    use parse, only: ARRHENIUS_RRATES_FLUID

    implicit none

    ! Passed Variables: Dummy argument format required by ODEPACK.
    !-----------------------------------------------------------------//
    ! (1) Number of ODEs to be solve
    ! (2) Fluid cell index
    ! (3) Number of solids (tfm or dpm)
    ! (4) Start of solve/don't solve flags for solids
    integer         , intent(in   ) :: a_NEq(:)
    ! Array of dependent variable initial values.
    double precision, intent(in   ) :: a_Y(:)
    ! Array of particle mass before calling stiff solver,
    ! used to normalize particle mass and species mass in particles
    double precision, intent(in   ) :: a_norm(:)
    ! Rate of change of dependent variables.
    double precision, intent(  out) :: a_YDOT(:)

    ! Local Variables:
    !-----------------------------------------------------------------//
    integer :: IJK  ! Fluid Cell Index
    integer :: H    ! Reaction loop counter
    integer :: L, M ! Global Phase index loop counters
    integer :: N    ! Global species index
    integer :: lN   ! Local reaction species index/loop counter
    integer :: mXfr ! Global phase index for mass transfer
    integer :: node ! Counter for YDOT

    ! User-defined reaction rates returned from USR_RATES
    double precision :: rates(no_of_rxns)

    ! > Single reaction.
    !...................................................................
    ! Rate of formation (+) or consumption (-) of a species (GENERIC).
    double precision :: rRate
    ! Rate of formation (+) or consumption (-) of each gas phase species.
    double precision :: rRg(dimension_n_g)
    ! Rate of formation (+) or consumption (-) of each solids phase species.
    double precision :: rRs(dimension_m, dimension_n_s)
    ! Heat of reaction for a reaction.
    double precision :: rHoR(0:dimension_m)
    ! Rate of interphase enthalpy transfer due to mass transfer.
    double precision :: rRxH(0:dimension_m, 0:dimension_m)

    ! > Cumulative over all reactions.
    !...................................................................
    ! Rate of formation (+) or consumption (-) of each gas phase species.
    double precision :: cRg(dimension_n_g)
    ! Rate of formation (+) or consumption (-) of each solids phase species.
    double precision :: cRs(dimension_m, dimension_n_s)
    ! Total heat of reaction
    double precision :: cHoR(0:dimension_m)

    ! UDF for Reaction Rates:
    external usr_rates

    double precision :: stat_wt, lto_SI

    lto_SI = merge(4.183925d3, 1.0d0, UNITS == 'SI')

    ! Initialize variables:
    ijk = a_NEq(2)

    a_YDOT   = ZERO

    cRg = ZERO;
    cRs = ZERO;

    cHoR = ZERO

    ! Map the current ODE independent variables to MFIX variables.
    call mapODEtoMFIX_tfm(size(a_NEq), a_NEq, size(a_Y), a_norm, a_Y)

    RATES  = ZERO

    ! Calculate user defined reaction rates.
    IF(ARRHENIUS_RRATES_FLUID) THEN
       CALL CALC_ARRHENIUS_RRATES(IJK, RATES)
    ELSE
       CALL USR_RATES(IJK, RATES)
    ENDIF

    ! Loop over reactions.
    rxn_lp: do h = 1, no_of_rxns

      ! Skip empty reactions
      if (Reaction(H)%nSpecies == 0) cycle rxn_lp

      ! Initialize local loop arrays
      rRg  = zero
      rRs  = zero
      rHOR = zero
      rRxH = zero

      ! Check if there are enough reactants for this reaction
      do lN = 1,Reaction(H)%nSpecies
        ! Global phase index.
        M = Reaction(H)%Species(lN)%pMap
        ! Global species index.
        N = Reaction(H)%Species(lN)%sMap

        if (Reaction(H)%Species(lN)%MWxStoich .lt. zero) then
          if (m == 0) then ! gas phase:
            if(X_g(ijk,n) .le. chem_min_species_fluid) then
              rates(h) = zero
              cycle rxn_lp
            endif
          else ! tfm solids m
            if (X_s(ijk,m,n) .le. chem_min_species_solid) then
              rates(h) = zero
              cycle rxn_lp
            endif
          endif
       endif
      enddo

      ! Calculate the rate of formation/consumption for each species.
      !-------------------------------------------------------------//
      do lN = 1, Reaction(H)%nSpecies

        ! Global phase index.
        M = Reaction(H)%Species(lN)%pMap
        ! Global species index.
        N = Reaction(H)%Species(lN)%sMap

        ! Index for interphase mass transfer. For a gas/solid reaction,
        ! the index is stored with the gas phase. For solid/solid mass
        ! transfer the index is stored with the source phase.
        mXfr = Reaction(H)%Species(lN)%mXfr

        ! Convert rate from moles to kg (g for csg units)
        rRate = rates(h) * Reaction(H)%Species(lN)%MWxStoich

        if (m == 0) then ! gas phase:

          ! Consumption/Formation of gas phase species.
          rRg(n) = rRg(n) + rRate
          ! Enthalpy transfer associated with mass transfer. (gas/solid)
          if (m /= mXfr) rRxH(m,mXfr) =  rRxH(m,mXfr) + &
                    rRate * calc_H(T_g(ijk),0,n)

        else ! tfm solids m

          ! Consumption of solids phase species.
          if (rRate < zero) then
            rRs(m,n) = rRs(m,n) + rRate
            ! Enthalpy transfer associated with mass transfer. (solid/solid) This
            ! is only calculated from the source (reactant) material.
            if (m /= mXfr) then
              if (m < mXfr) then
                rRxH(M,mXfr) =  rRxH(M,mXfr) +               &
                  rRate * calc_H(T_s(ijk,m),m,n) *         &
                  Reaction(H)%Species(lN)%xXfr
              else
                rRxH(mXfr,m) =  rRxH(mXfr,m) -                &
                  rRate * calc_H(T_s(ijk,m),m,n) *          &
                  Reaction(H)%Species(lN)%xXfr
              endif
            endif
          else
            ! Formation of solids phase species.
            rRs(M,N) = rRs(M,N) + rRate
          endif
        endif
      enddo ! Loop of species

      ! Copy the single reaction rates of formation/consumption to
      ! the total (accumulative) rates of formation/consumption
      ! arrays. This ensures that any reactions without sufficient
      ! reactants are not included.
      cRg = cRg + rRg
      cRs = cRs + rRs

      ! Calculate and store the heat of reaction.
      !-------------------------------------------------------------//
      ! Automated heat of reaction calculations
       if (Reaction(H)%calc_dH) then

         ! Loop over reaction species.
         do lN = 1, Reaction(H)%nSpecies

           ! Global phase index.
           m = Reaction(h)%Species(lN)%pMap
           ! Global species index.
           n = Reaction(h)%Species(lN)%sMap
           ! Rate of formation/consumption for species N
           rRate = rates(h) * Reaction(h)%Species(lN)%MWxStoich

           if (m == 0) then
             ! Gas phase enthalpy change
             rHoR(0) = rHoR(0) + rRate * calc_H(T_g(ijk),0,n)
           else
             ! Solid phase enthalpy change
             rHoR(M) = rHoR(M) + rRate * calc_H(T_s(ijk,m),m,n)
           endif
         enddo

         ! Complete the skew-symmetric for enthalpy transfer
         ! with mass transfer
         do m=1, mmax
           do l=0, m-1
             rRxH(m,l) = - rRxH(l,m)
           enddo
         enddo

         ! Apply enthalpy transfer associated with mass transfer
         ! to get the complete heat of reaction for Reaction H.
         do l=0, mmax
           do m = 0, mmax
             if(l /= m) rHOR(l) = rHOR(l) - rRxH(l,m)
           enddo
         enddo

         ! Convert the heat of reaction to the appropriate units
         ! (if SI), and store in the global array.
         cHOR(:) = cHOR(:) + rHOR(:)*lto_SI

      else ! User-defined heat of reaction.
        cHoR(:) = cHoR(:) + Reaction(H)%HoR(:) * RATES(H)
      endif

    enddo rxn_lp ! loop over reactions.

    ! Save rates so they can be written in vtu files and/or monitors
    call save_fluid_rrates_cell_data(ijk, rates)

    fill_ydot: block
    !---------------------------------------------------------------//
    ! Map the reaction rates into the ODE solver

    use fldvar, only: rop_g
    use fldvar, only: rop_s, ro_s

    use physprop, only: C_pg, C_ps

    node = 1

    ! Bulk density:  rop_g
    a_YDOT(node) = sum(cRg);                           node = node + 1

    ! Temperature: T_g
    a_YDOT(node) = -cHoR(0)/(rop_g(ijk)*C_pg(ijk));    node = node + 1

    ! Bulk species density: rop_g * X_g
    do n=1,nmax(0)
      a_YDOT(node) = cRg(n);                           node = node + 1
    enddo

    do m = 1, a_NEq(3)

      if (a_NEq(3 + m) == 1) THEN

        ! Solids temperature: T_s
        if (rop_s(ijk,m) > zero_ep_s*ro_s(ijk,m)) then
           a_YDOT(node) = -cHoR(m)/(rop_s(ijk,m)*C_ps(ijk,m))
        else
           a_YDOT(node) = 0.0d0
        endif;                                           node = node + 1

        ! Solids bulk density: rop_s
        a_YDOT(node) = sum(cRs(m,:))/a_norm(m);                  node = node + 1

        ! Bulk species mass rop_s * X_s
        do n=1, nmax(m)
          a_YDOT(node) = cRs(m,n)/a_norm(m);                     node = node + 1
        enddo

      endif
    enddo
  end block fill_ydot

    return
  end subroutine stiff_chem_rrates_tfm


  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Module name: mapMFIXtoODE                                         !
  !  Author: J.Musser                                 Date: 07-Feb-13  !
  !                                                                    !
  !  Purpose: This routine maps MFIX variables into the ODE array.     !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine mapMFIXtoODE_tfm(a_NEq_d, a_NEq, a_ODE_Vars_d, a_ODE_norm, a_ODE_Vars)

    ! Global Variables:
    !-----------------------------------------------------------------//
    use fldvar, only: ROP_g, T_g, X_g
    use fldvar, only: ROP_s, T_s, X_s

    use physprop, only: nmax

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
    integer :: node, ijk, m, n
    !...................................................................

    ! Initialize.
    ijk = a_NEq(2)
    a_ODE_Vars = 0.0d0

    node = 1

    ! Gas phase bulk density.
    a_ODE_Vars(node) = rop_g(ijk);               node = node + 1

    ! Gas phase temperature.
    a_ODE_Vars(node) = T_g(ijk);                 node = node + 1

    ! Gas phase species mass.
    do n=1,nmax(0)
      a_ODE_Vars(node) = rop_g(ijk)*X_g(ijk,n);  node = node + 1
    enddo
    do m = 1, a_NEq(3)
        if (a_NEq(3 + m) == 1) then
            ! Solids temperature.
            a_ODE_Vars(Node) = T_s(ijk,m);           node = node + 1

            ! Solids bulk density.
            a_ODE_Vars(node) = rop_s(ijk,m)/a_ODE_norm(m);       node = node + 1

            ! Solids species mass.
            do n=1,nmax(m)
                a_ODE_Vars(node) = &
                    rop_s(ijk,m)*X_s(ijk,m,n)/a_ODE_norm(m);     node = Node + 1
            enddo

        endif
    enddo

  end subroutine mapMFIXtoODE_tfm

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Module name: mapODEtoMFIX                                         !
  !  Author: J.Musser                                 Date: 07-Feb-13  !
  !                                                                    !
  !  Purpose: This is a driver routine for mapping variables stored in !
  !  the ODE array back to MFIX field variables.                       !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine mapODEtoMFIX_tfm(a_NEq_d, a_NEq, a_ODE_Vars_d, a_ODE_norm, a_ODE_Vars)

    ! Global Variables:
    !-----------------------------------------------------------------//
    use fldvar, only: rop_g, ro_g, t_g, X_g, ep_g, P_g
    use fldvar, only: rop_s, ro_s, t_s, X_s, ep_s

    use physprop, only: mmax, nmax

    ! Gas constant (cal/g.K)
    use constant, only: GAS_CONST

    ! Mixture and species molecular weights
    use physprop, only: MW_MIX_g, MW_g
    ! Baseline/Unreacted solids density.
    use physprop, only: RO_s0, X_s0, inert_species

    use run, only: solve_ros

    ! Variable solids density calculation.
    use eos, only: eoss

    ! Global Parameters:
    !-----------------------------------------------------------------//
    use param1,   only : zero, small_number
    use param1,   only :  one, large_number

    use toleranc, only: zero_ep_s

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
    integer :: node, ijk, m, n, p, pid
    double precision :: sum_eps, inv_rops, rops_normalize
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
      X_g(ijk,n) = max(0.0,X_g(ijk,n))
    enddo

    sum_eps = zero

    do m = 1, a_NEq(3)

      if (a_NEq(3 + m) == 1) then

        ! Solids temperature.
        T_s(ijk,m) = a_ODE_Vars(node);                     node = node + 1

        ! Solids bulk density.
        rops_normalize = a_ODE_Vars(node)
        rop_s(ijk,m) = a_ODE_Vars(node)*a_ODE_norm(m);                 node = node + 1

        if (rops_normalize > small_number) then
          inv_rops = one/rops_normalize
        else
          inv_rops = zero
        endif

        ! Solids phase species mass fractions.
        do n=1,nmax(m)
          X_s(ijk,m,n) = inv_rops*a_ODE_Vars(node);      Node = Node + 1
          X_s(ijk,m,n) = max(0.0,X_s(ijk,m,n))
        enddo

        ! Update the solids density from species composition. This update must
        ! come after the species data is mapped back into the MFIX variables.
        if(solve_ros(m)) ro_s(ijk,m) = eoss(ro_s0(m),              &
          x_s0(m,inert_species(m)), x_s(ijk,m,inert_species(m)))
      endif

      sum_eps = sum_eps + ep_s(ijk,m)
    enddo

    ! Calculate the gas volume fraction from solids volume fractions.
    ep_g(ijk) = min(max(one - sum_eps, zero), one)

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
  end subroutine mapODEtoMFIX_tfm

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Module name: reportODEVar_tfm                                     !
  !  Author: J.Musser                                 Date: 22-May-23  !
  !                                                                    !
  !  Purpose: Checks the values in the ODE array. Currently not used.  !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine reportODEVar_tfm (a_NEq_d, a_NEq, a_ODE_Vars_d, a_norm, a_ODE_Vars)

    use run, only: tfm_solids, dem_solids, pic_solids
    use run, only: species_eq

    use rxns, only: species_alias_g, species_alias_s
    use physprop, only: nmax

    use fldvar, only: ro_s
    use toleranc, only: zero_ep_s

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

    double precision :: ep_rho, inv_rho_eps, rho_eps, Xs_rho_eps, Xs
    double precision :: sum_Xg, sum_Xs, sum_Xs_rho_eps

    integer :: node, m, n, ijk

    node = 1
    ijk = a_NEq(2)

    ! Gas phase density.
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

    do m=1, a_NEq(3)

      ! Solids temperature
      write(*,"(I3,2x,A24,F12.4)") node, "T_s:", a_ODE_Vars(node)
      node = node + 1

      if (a_NEq(3 + m) == 1) then

        ! Solids bulk density fraction.
        rho_eps = a_ODE_Vars(Node)*a_norm(m)
        node = node + 1
        write(*,"(I3,2x,A24,3x,E18.8)") node, "ROP_s:", rho_eps

        if (rho_eps > zero_ep_s*ro_s(ijk,m)) then
           inv_rho_eps = 1.0d0 / rho_eps
        else
           inv_rho_eps = 0.0d0
        endif

        ! Solids phase species mass fractions.
        sum_Xs = 0.
        sum_Xs_rho_eps = 0.
        do n=1,nmax(m)

          Xs_rho_eps = a_ODE_Vars(node)*a_norm(m)
          node = node + 1

          Xs = Xs_rho_eps*inv_rho_eps
          write(*,"(I3,2x,A24,3x,E18.8)") node, "X_s("//trim(species_alias_s(m,n))//"):", Xs

          sum_Xs_rho_eps = sum_Xs_rho_eps + Xs_rho_eps
          sum_Xs    = sum_Xs + Xs

        enddo

      endif
    enddo

    write(*,"(2x,A16,F12.4)") "   sum Xg:", sum_Xg
    write(*,"(2x,A16,F12.4)") "   sum Xs:", sum_Xs
    write(*,"(2x,A16,F12.4)") "sum Xs_mp:", sum_Xs_rho_eps

  end subroutine reportODEVar_tfm

end module stiff_chem_tfm
