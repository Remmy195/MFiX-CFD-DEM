module stiff_chem

! Runtime Flags:
!---------------------------------------------------------------------//
! Flag to invoke stiff chemistry solver.
  logical :: stiff_chemistry

! ODEPACK Controlling parameters:
!---------------------------------------------------------------------//
! The maximum number of steps ODEPACK may use to integrate.
  integer :: stiff_chem_max_steps

! ! Absolute tolerance
!   double precision :: stiff_chem_abs_tol
!
! ! Relative tolerance
!   double precision :: stiff_chem_rel_tol
!
! ! Minimum reaction rate. If all reaction rates are below this threshold, a
! ! cell or particle is said to be non-reactive and no ODE is set for this
! ! cell/particle.
!   double precision :: stiff_chem_min_rate
!


contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  ! Module name: stiff_chem_solver                                     !
  ! Purpose: Called in time_march.f to do rxns calcs                   !
  !                                                                    !
  ! Author: J.Musser                                  Date: 21-May-23  !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine stiff_chem_solver (ODE_dt, iErr)

    use output, only: full_log

    use compar, only: istart, jstart, kstart
    use compar, only: iend,   jend,   kend

    use functions, only: funijk
    use functions, only: fluid_at

    use compar,    only: dead_cell_at

    use StiffChemClass, only: TSolver
    use StiffChemClass, only: clearStats
    use StiffChemClass, only: reportStats

    implicit none

    ! Passed Variables:
    !------------------------------------------------------------------!
    double precision, intent(in   ) :: ODE_dt  ! Time integral length.
    integer         , intent(  out) :: iErr    ! Error Flag


    ! Local Variables:
    !------------------------------------------------------------------!

    ! Fluid Cell Index
    integer :: i,j,k,ijk

    Type(TSolver) :: StiffSolver

    double precision :: timeIn, timeOut
    if (full_log) call clearStats()

    do k=kstart, kend
      do i=istart, iend
        do j=jstart, jend

          if (dead_cell_at(i,j,k)) cycle

          ijk = funijk(i,j,k)

          if( fluid_at(ijk)) then

            timeIn  = 0.d0
            timeOut = ODE_dt

            StiffSolver = TSolver(ijk, stiff_chem_max_steps)
            if ( StiffSolver%isActive() ) then

              call StiffSolver%PrepareForSolve()

              call StiffSolver%Solve(timeIn, timeOut)

              call StiffSolver%Finalize()

            endif
          endif

        enddo ! j-loop
      enddo ! i-loop
    enddo ! k-loop

    call finalize_stiff_solver()

    if (full_log) call reportStats()

    iErr = 0

    return
  end subroutine stiff_chem_solver


  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Module name: FINALIZE_STIFF_SOLVER                                !
  !                                                                    !
  !  Purpose:                                                          !
  !                                                                    !
  !  Author: J.Musser                                 Date: 07-Feb-13  !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine finalize_stiff_solver

    use fldvar, only: ro_g, rop_g, T_g, X_g, ep_g, p_g
    use fldvar, only: ro_s, rop_s, T_s, X_s

    use physprop, only: MW_mix_g

    use geometry, only: ijkmax2
    use physprop, only: mmax, nmax

    use utilities, only: bound_x
    use sendrecv, only: send_recv

    implicit none

    ! Local loop indices.
    INTEGER :: M    ! Solids phase index
    INTEGER :: NN   ! Species index

    call send_recv(ep_g,2)
    call send_recv(ro_g,2)
    call send_recv(rop_g,2)

    call send_recv(T_g,2)
    call send_recv(P_g,2)
    call send_recv(MW_mix_g,2)

    do nn=1,nmax(0)
       call send_recv(X_g(:,nn),2)
       call bound_x (X_g(1,nn), ijkmax2)
    enddo

    DO M = 1, MMAX
      ! Solids volume fraction. (Constant Solids Density)
      call send_recv(ro_s(:,m),2)
      call send_recv(rop_s(:,m),2)
      ! Solids temperature.
      call send_recv(T_s(:,m),2)
      ! Solids phase species mass fractions.
      do nn=1,nmax(m)
        call send_recv(X_s(:,m,nn),2)
        call bound_X (X_s(1,m,nn), ijkmax2)
      enddo
    enddo

    return
  end subroutine finalize_stiff_solver

end module stiff_chem
