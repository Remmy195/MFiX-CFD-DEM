#include "error.inc"

module StiffChemClass

  implicit none
  ! Declare the values used to normalize the particle mass and
  ! species mass in solid. Here, the particle mass obtained before
  ! calling stiff solver is used..
  double precision, allocatable, public :: m_norm(:)

  ! Counters to track status of solves
  integer, private :: m_failed
  integer, private :: m_successful
  integer, private :: m_incomplete

  ! Counters to track number of hetero/homogeneous reactions
  integer, private :: m_hetrgns
  integer, private :: m_homogns

  double precision, private :: m_minTolSF
  double precision, private :: m_maxTolSF

  double precision, private :: m_startTime

  integer, private :: m_minSteps
  integer, private :: m_maxSteps

  integer, private :: m_binnedStepCounts(5)
  integer, private :: m_binnedSolveTimes(7)

  ! Frequency to report the number of steps distribution.
  integer, parameter, private :: m_reportFreq = 10
  integer, private :: m_report = 0

  ! Solver type
  type, public :: TSolver

    ! The maximum number of internal steps ODEPack may use to
    ! integrate over the time interval. The default is 500.
    integer, private :: m_max_steps

    ! Flag: cell has one or more reactions
    logical, private :: m_active = .False.

    ! (1) :: Number of ODEs
    ! (2) :: Fluid cell index (IJK) passed into ODEPack
    ! (:) :: Flag for solving solids
    integer, allocatable, private :: m_NEq(:)

    ! ODE variables
    double precision, allocatable :: m_ODE_Vars(:)
    double precision, allocatable :: m_ODE_Vars_backup(:)

    ! ODEPack settings .......................................

    ! Jacobian type indicator specifies how the Jacobian matrix
    ! df/dy will be treated, if and when DLSODA requires this matrix.
    ! 2 :: internally generated full Jacobian
    !      (using NEQ extra calls to F per df/dy value).
    integer, private :: m_jacobian_type = 2

    ! Type of Error control.
    ! 1 :: EWT(i) = RTOL*ABS(Y(i)) + ATOL
    ! 2 :: EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)
    integer, private :: m_ODE_iTol = 2

    ! Relative error tolerance parameter.
    double precision, private :: m_ODE_rTol

    ! Absolute error tolerance parameter.
    double precision, allocatable, private :: m_ODE_aTol(:)

    ! Declared length of rWork.
    integer, private :: m_ODE_LRW
    double precision, allocatable, private :: m_ODE_rWork(:)

    ! Declared length of iWork.
    integer, private :: m_ODE_LIW
    integer, allocatable, private :: m_ODE_iWork(:)

    ! Stats ..................................................

    ! local timer
    double precision, private :: m_startTime

    contains

      procedure, public :: isActive
      procedure, public :: PrepareForSolve
      procedure, public :: Solve
      procedure, public :: Finalize

  end type TSolver

  ! This constructor allows us to create the data object
  ! on the fly, much like c++.
  interface TSolver
    module procedure TSolverConstructor
  end interface TSolver

  public :: clearStats
  public :: reportStats

contains

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Subroutine: TSolverConstructor                                    !
  !  Author: J. Musser                             Date: June 6, 2023  !
  !                                                                    !
  !  Purpose: Creates and returns a instance of the stiff solver.      !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  function TSolverConstructor(a_ijk, a_max_steps) &
  result(NewSolver)

    use debug, only: assert
    use derived_types,  only: pic
    use des_rxns, only: no_of_des_rxns
    use discretelement, only: m_p => pmass
    use discretelement, only: pinc, pijk
    use fldvar, only: ep_s, rop_s
    use parse, only: arrhenius_rrates_fluid, arrhenius_rrates_des
    use physprop, only: mmax, nmax
    use run, only: tfm_solids, dem_solids, pic_solids, species_eq
    use rxns, only: no_of_rxns
    use time_mod, only: get_elapsed_time
    use toleranc, only: stiff_chem_abs_tol, stiff_chem_rel_tol, stiff_chem_min_rate
    use toleranc, only: stiff_chem_min_rate_dpm, zero_ep_s

    integer, intent(in   ) :: a_ijk
    integer, intent(in   ) :: a_max_steps

    double precision :: rates(no_of_rxns)
    double precision :: rates_dpm(no_of_des_rxns)

    type(TSolver) :: NewSolver

    integer :: NEq_d, ODE_d
    integer :: m, pid, p
    integer :: hetrgns

    NewSolver%m_startTime = get_elapsed_time()

    ! Store the number of max steps
    NewSolver%m_max_steps = a_max_steps

    ! Compute the dimension of the NEq array
    NEq_d = 3 ! Number of ODEs, fluid cell IJK, solids count

    if (tfm_solids) then
      ! Flags for continuum solids
      NEq_d = NEq_d + mmax

    else if (dem_solids .or. pic_solids) then
      ! Flags for each dpm entity in cell ijk
      NEq_d = NEq_d + pinc(a_ijk)

    endif

    hetrgns = 0

    allocate(NewSolver%m_NEq(NEq_d), source=0)

    ! Compute the number of ODEs to solve
    ! Fluid density, temperature, species
    ODE_d = 2 + nmax(0)

    if (allocated(m_norm)) deallocate(m_norm)
    if (tfm_solids) then
      NewSolver%m_NEq(3) = mmax
      ! temperature for all solids phases
      ODE_d = ODE_d + mmax
      allocate(m_norm(mmax), source=0D0)

      do m=1, mmax
        NewSolver%m_NEq(3+m) = 0
        m_norm(m) = rop_s(a_ijk, m)
        if(species_eq(m) .and. (ep_s(a_ijk,m) > zero_ep_s)) then
          ! flag heterogeneous reaction
          hetrgns = 1
          ! flag to solve this solid density and species
          NewSolver%m_NEq(3+m) = 1
          ! density and species mass
          ODE_d = ODE_d + 1 + nmax(m)
        endif
      enddo

    else if (dem_solids .or. pic_solids) then
      NewSolver%m_NEq(3) = pinc(a_ijk)
      if(pinc(a_ijk)>0) then
         allocate(m_norm(pinc(a_ijk)), source=0D0)
      else
         allocate(m_norm(1), source=0D0)
      endif

      do p=1, pinc(a_ijk)
        pid = pic(a_ijk)%p(p)
        m = pijk(pid,5)
        call assert(a_ijk .eq. pijk(pid,4),                       &
          & "** StiffChem::StiffChem() Initialization error: " // &
          & "Cell and particle ijk values do not match!",         &
          & a_ijk, pijk(pid,4))

        NewSolver%m_NEq(3+p) = 0
        m_norm(p) = m_p(pid)

        if (species_eq(m)) then
          ! Dry run- compute rates for particle.
          rates_dpm = 0.0d0
          if(arrhenius_rrates_des) then
             call calc_arrhenius_rrates_des(pid, m, a_ijk, rates_dpm)
          else
             call usr_rates_des(pid, m, a_ijk, rates_dpm)
          endif

          if (maxval(rates_dpm) > stiff_chem_min_rate_dpm) then
            NewSolver%m_active = .True.
            ! flag heterogeneous reaction
            hetrgns = 1
            ! flag to solve this solid
            NewSolver%m_NEq(3+p) = 1
            ! temperature, total mass and species of this particle
            ODE_d = ODE_d + 2 + nmax(m)
          endif
        endif
      enddo
    endif

    ! Increment homogeneous / heterogeneous reaction counters
    m_homogns = m_homogns + (1-hetrgns)
    m_hetrgns = m_hetrgns + hetrgns

    call assert(ODE_d .gt. 0,                                 &
      & "** StiffChem::StiffChem() Initialization error: " // &
      & "Number of ODEs to be solved is <= 0!",               &
      & ODE_d, 0)

    ! If not active, test if there are continuum reactions
    if ( .not.NewSolver%m_active ) then
      rates = 0.0d0
      if(arrhenius_rrates_fluid) then
         call calc_arrhenius_rrates(a_ijk, rates)
      else
         call usr_rates(a_ijk, rates)
      endif

      if (maxval(rates) > stiff_chem_min_rate) NewSolver%m_active = .True.

    endif

    NewSolver%m_NEq(1) = ODE_d ! number of ODEs to solve
    NewSolver%m_NEq(2) = a_ijk ! fluid cell index

    if ( NewSolver%m_active ) then
      ! Array for ODE variables
      allocate(NewSolver%m_ODE_Vars(ODE_d), source=0D0)
      ! Declared length of rWork.
      NewSolver%m_ODE_LRW = 22 + ODE_d * max(16, ODE_d + 9)
      allocate(NewSolver%m_ODE_rWork(NewSolver%m_ODE_LRW ), source=0D0)
      ! Declared length of iWork.
      NewSolver%m_ODE_LIW = 20 + ODE_d
      allocate(NewSolver%m_ODE_iWork(NewSolver%m_ODE_LIW ), source=0)
      ! Relative error tolerance parameter.
      NewSolver%m_ODE_rTol = stiff_chem_rel_tol
      ! Absolute error tolerance parameter.
      allocate(NewSolver%m_ODE_aTol(ODE_d), source=0D0)
      NewSolver%m_ODE_aTOL(:) = stiff_chem_abs_tol ! All Equations
    endif
  end function TSolverConstructor

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Subroutine: isActive                                              !
  !  Author: J. Musser                             Date: June 6, 2023  !
  !                                                                    !
  !  Purpose: Access to member variable m_active. This lets the IJK    !
  !  loop skip cells that have no reactions.                           !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  logical function isActive (a_this)
    class(TSolver), intent(in   ) :: a_this
    isActive = a_this%m_active
  end function isActive

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Subroutine: PrepareForSolve                                       !
  !  Author: J. Musser                             Date: June 6, 2023  !
  !                                                                    !
  !  Purpose: Map MFIX field variables to the ODE array, initialize    !
  !  the ODEPack work arrays, and set various ODEPack settings.        !
  !                                                                    !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine PrepareForSolve (a_this)
    class(TSolver), intent(inout) :: a_this
    ! Map MFIX variables to ODE variables.
    call mapMFIXtoODE(size(a_this%m_NEq), a_this%m_NEq, &
    &  a_this%m_NEq(1), m_norm, a_this%m_ODE_Vars)

    ! Store a copy of the original field variables so we can
    ! reset them if the stiff solver fails.
    a_this%m_ODE_Vars_backup = a_this%m_ODE_Vars

    ! Initialize the real and integer work arrays
    a_this%m_ODE_rWork(:) = 0.
    a_this%m_ODE_iWork(:) = 0.

    ! Optional inputs for ODEPack

    ! (5)  flag to generate extra printing at method switches.
    !      IXPR = 0 means no extra printing (the default).
    !      IXPR = 1 means print data on each switch.
    !      T, H, and NST will be printed on the same logical
    !      unit as used for error messages.
    !a_this%m_ODE_iWork(5) = 0

    ! (6)  maximum number of (internally defined) steps
    !      allowed during one call to the solver.
    !      The default value is 500.
    a_this%m_ODE_iWork(6) = a_this%m_max_steps

    ! (7)  maximum number of messages printed (per problem)
    !      warning that T + H = T on a step (H = step size).
    !      This must be positive to result in a non-default
    !      value.  The default value is 10.
    !a_this%m_ODE_iWork(7) = 0

    ! (8)  the maximum order to be allowed for the nonstiff
    !      (Adams) method.  the default value is 12.
    !      if MXORDN exceeds the default value, it will
    !      be reduced to the default value.
    !      MXORDN is held constant during the problem.
    !a_this%m_ODE_iWork(8) = 0

    ! (9)  the maximum order to be allowed for the stiff
    !      (BDF) method.  The default value is 5.
    !      If MXORDS exceeds the default value, it will
    !      be reduced to the default value.
    !      MXORDS is held constant during the problem.
    !a_this%m_ODE_iWork(9) = 0

  end subroutine PrepareForSolve

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Subroutine: Solve                                                 !
  !  Author: J. Musser                             Date: June 6, 2023  !
  !                                                                    !
  !  Purpose: Drives the call to ODEPack and checks the solve status.  !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine Solve(a_this, a_timeIn, a_timeOut)

    use error_manager, only: err_msg, log_message, loglevel_error

    class(TSolver)  , intent(inout) :: a_this
    double precision, intent(inout) :: a_timeIn, a_timeOut

    ! Explicit interface for ODEPack
    !-----------------------------------------------------------------//
    interface
      subroutine dlsoda (F, NEQ, Y, timeIn, timeOut, iTol, rTol, aTol, &
        iTask,iState, iOpt, rWork, LRW, iWork, LIW, Jac, JT)
        external F
        integer :: iTol, iTask, iState, iOpt, LRW, LIW, JT
        integer, dimension(2) :: NEQ
        integer, dimension(LIW) :: IWORK
        double precision :: timeIn, timeOut, JAC, rTol
        double precision, dimension(LRW) :: RWORK
        double precision, dimension(NEQ(1)) :: Y, ATOL
      end subroutine dlsoda
    end interface

    ! Jacobian Routine (not used)
    double precision :: jacobian_routine

    ! Index specifying the ODEPack task.
    integer, parameter :: iTask = 1

    ! Flag indicating optional inputs are used.
    integer, parameter :: iOpT = 1

    ! Specifies the state of ODEPack
    integer :: iState

    integer, parameter :: max_iter = 3
    integer :: iter

    double precision :: tolSF

    tolSF = 1.0d0

    iter = 0
    iState = 1

    do while (iState > 0 .and. iter < max_iter)
      iter = iter + 1
      ! Integrate flow field variables to incorporate reactions.
      call dlsoda(stiff_chem_rrates, a_this%m_NEq, a_this%m_ODE_Vars, &
        a_timeIn, a_timeOut, &
        a_this%m_ODE_iTol, a_this%m_ODE_rTol, a_this%m_ODE_aTol, &
        iTask, iState, iOpt,  &
        a_this%m_ODE_rWork, a_this%m_ODE_LRW, &
        a_this%m_ODE_iWork, a_this%m_ODE_LIW, &
        jacobian_routine, a_this%m_jacobian_type)

      if ( ODEVarHasNaNs(a_this%m_NEq(1), a_this%m_ODE_Vars) ) then
          m_failed = m_failed + 1
          a_this%m_ODE_Vars = a_this%m_ODE_Vars_backup
      else
        ! On output, iSTATE has the following values and meanings.
        select case (iState)
        case(1,2)
          ! 2  means the integration was performed successfully.
          ! 1  means nothing was done; TOUT = T and ISTATE = 1 on input.
          iState = 0 ! Successful integration
          m_successful = m_successful + 1

        case(-1)
          !-1  an excessive amount of work (more than MXSTEP steps)
          !    was done on this call, before completing the requested
          !    task, but the integration was otherwise successful as
          !    as far as time T (saved in timeIn).
          iState = 2 ! Continue the calculation.

        case(-2)
          !-2  too much accuracy was requested for the precision
          !    of the machine being used. This was detected before
          !    completing the requested task, but the integration
          !    was successful as far as T.

          ! ODEPack provides a tolerance scale factor when a too much
          ! accuracy was requested. The solver is deemed likely to
          ! succeed if rTol and aTol are scaled by this factor for the
          ! next call.
          tolSF = a_this%m_ODE_rWork(14)
          a_this%m_ODE_rTol = tolSF*a_this%m_ODE_rTol
          a_this%m_ODE_aTol = tolSF*a_this%m_ODE_aTol
          iState = 3 ! Try again with scaled tolerances

        case(-3)
          !-3  illegal input was detected.
          ! The interface to ODEPack is static and assumed correct so
          ! we treat this error as fatal.
          write(err_msg, 1100)
          call log_error()

        case(-4,-5,-6)
          !-4  there were repeated error test failures on one attempted
          !    step, before completing the requested task. The problem
          !    may have a singularity, or the input may be inappropriate.
          !
          !-5  there were repeated convergence test failures on
          !    one attempted step, before completing the requested
          !    task. This may be caused by an inaccurate Jacobian
          !    matrix, if one is being used.

          !-6  EWT(i) became zero for some i during the integration.
          !    Pure relative error control (ATOL(i)=0.0) was requested
          !    on a variable which has now vanished.

          ! In all cases, the integration was successful as far as
          ! time T, but integration is incomplete.

          m_incomplete = m_incomplete + 1

          call reportODEVar(size(a_this%m_NEq), a_this%m_NEq, &
          &  a_this%m_NEq(1), m_norm, a_this%m_ODE_Vars)

        case(-7)
          !-7  the length of RWORK and/or IWORK was too small

          ! The interface to ODEPack is static and assumed correct so
          ! we treat this error as fatal.
          write(err_msg, 1110)
          call log_error()
        end select
      endif
    enddo ! while iState > 0

    m_minTolSF = min(tolSF, m_minTolSF)
    m_maxTolSF = max(tolSF, m_maxTolSF)

    if (iter >= max_iter .and. iState .ne. 0) then
      m_incomplete = m_incomplete + 1
    endif

 1100 format('Error 1100: ODEPack reported a fatal error:',/,&
             'ISTATE = -3: Illegal input detected.')

 1110 format('Error 1110: ODEPack reported a fatal error:',/,&
             'ISTATE = -7: Length of work arrays too small!')

  end subroutine Solve

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Subroutine: Finalize                                              !
  !  Author: J. Musser                             Date: June 6, 2023  !
  !                                                                    !
  !  Purpose: Maps the ODE array back to MFIX field arrays and bins    !
  !  the number of steps needed to solve this cell as well as the time !
  !  needed for the solve.                                             !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine Finalize(a_this)
    use time_mod, only: get_elapsed_time
    class(TSolver), intent(inout) :: a_this
    ! IWORK(11) the number of steps taken for the problem so far.
    integer :: NST
    double precision :: solveTime

    call mapODEtoMFIX(size(a_this%m_NEq), a_this%m_NEq, &
    &  a_this%m_NEq(1), m_norm, a_this%m_ODE_Vars)

    NST = a_this%m_ODE_iWork(11)

    m_minSteps = min(m_minSteps, NST)
    m_maxSteps = max(m_maxSteps, NST)

    if (NST < 10) then
       m_binnedStepCounts(1) = m_binnedStepCounts(1) + 1
    else if (NST < 100) then
       m_binnedStepCounts(2) = m_binnedStepCounts(2) + 1
    else if (NST < 1000) then
       m_binnedStepCounts(3) = m_binnedStepCounts(3) + 1
    else if (NST < 10000) then
       m_binnedStepCounts(4) = m_binnedStepCounts(4) + 1
    else
       m_binnedStepCounts(5) = m_binnedStepCounts(5) + 1
    endif

    solveTime = get_elapsed_time() - a_this%m_startTime

    if (solveTime < 1.0d-3) then
       m_binnedSolveTimes(1) = m_binnedSolveTimes(1) + 1
    else if (solveTime < 1.0d-2) then
       m_binnedSolveTimes(2) = m_binnedSolveTimes(2) + 1
    else if (solveTime < 1.0d-1) then
       m_binnedSolveTimes(3) = m_binnedSolveTimes(3) + 1
    else if (solveTime < 1.0d0 ) then
       m_binnedSolveTimes(4) = m_binnedSolveTimes(4) + 1
    else if (solveTime < 1.0d+1) then
       m_binnedSolveTimes(5) = m_binnedSolveTimes(5) + 1
    else if (solveTime < 1.0d+2) then
       m_binnedSolveTimes(6) = m_binnedSolveTimes(6) + 1
    else
       m_binnedSolveTimes(7) = m_binnedSolveTimes(7) + 1
    endif

  end subroutine Finalize

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Subroutine: ODEVarHasNaNs                                     !
  !  Author: J. Musser                             Date: June 6, 2023  !
  !                                                                    !
  !  Purpose: This function checks the ODE array for NaNs.             !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  function ODEVarHasNaNs (a_ODE_d, a_ODE_Vars) &
  result(a_NaN)

    use utilities, only: mfix_isNaN

    integer,          intent(in   ) :: a_ODE_d
    double precision, intent(in   ) :: a_ODE_Vars(a_ODE_d)

    logical :: a_NaN

    integer :: node

    ! Check the results for NaNs.
    do node=1, a_ODE_d
      if (mfix_isNaN(a_ODE_Vars(node)) ) then
        a_NaN = .True.
        return
      endif
    enddo

    a_NaN = .False.
    return

  end function ODEVarHasNaNs

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Subroutine: stiff_chem_rrates                                     !
  !  Author: J. Musser                             Date: June 6, 2023  !
  !                                                                    !
  !  Purpose: This is a driver routine to computing reaction rates     !
  !  needed for ODEPack. Different routines are used for tfm and dpm.  !
  !  Fluid-only runs use the TFM interface.                            !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine stiff_chem_rrates (a_NEq, a_Time, a_Y, a_YDOT)

    use run, only: dem_solids, pic_solids

    use stiff_chem_dpm, only: stiff_chem_rrates_dpm
    use stiff_chem_tfm, only: stiff_chem_rrates_tfm

    ! Passed Variables: Dummy argument format required by ODEPack.
    !-----------------------------------------------------------------//
    ! (1) Number of ODEs to be solve
    ! (2) Fluid cell index
    ! (3) Number of solids (tfm or dpm)
    ! (4) Start of solve/don't solve flags for solids
    integer         , intent(in   ) :: a_NEq(*)
    ! Independent variable (not used)
    double precision, intent(in   ) :: a_Time
    ! Array of dependent variable initial values.
    double precision, intent(in   ) ::a_Y(*)
    ! Rate of change of dependent variables.
    double precision, intent(  out) ::a_YDOT(*)

    double precision :: d_Time
    integer :: NEq_d, ODE_d

    NEq_d = 3 + a_NEq(3)
    ODE_d = a_NEq(1)

    if (dem_solids .or. pic_solids) then
      call stiff_chem_rrates_dpm(a_NEq(1:NEq_d), &
             a_Y(1:ODE_d), m_norm(:), a_YDOT(1:ODE_d))
    else
      call stiff_chem_rrates_tfm(a_NEq(1:NEq_d), &
             a_Y(1:ODE_d), m_norm(:), a_YDOT(1:ODE_d))
    endif

    ! This is inserted to suppress 'unused variable warnings'
    ! from the compiler. We know we don't use it but it is
    ! part of the function signature for ODEPack so we have
    ! to include it.
    if (.False.) d_Time = a_Time

  end subroutine stiff_chem_rrates

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Subroutine: mapMFIXtoODE                                          !
  !  Author: J. Musser                             Date: June 6, 2023  !
  !                                                                    !
  !  Purpose: This is a driver routine to copy MFIX field variables to !
  !  the ODE array. Different routines are used for tfm and dpm.       !
  !  Fluid-only runs use the tfm interface.                            !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine mapMFIXtoODE (a_NEq_d, a_NEq, a_ODE_d, a_ODE_norm, a_ODE_Vars)

    use stiff_chem_dpm, only: mapMFIXtoODE_dpm
    use stiff_chem_tfm, only: mapMFIXtoODE_tfm

    use run, only: dem_solids, pic_solids

    implicit none

    ! Passed Variables:
    !-----------------------------------------------------------------//
    ! (1) Number of ODEs to be solve
    ! (2) Fluid cell index
    integer         , intent(in   ) :: a_NEq_d
    integer         , intent(in   ) :: a_NEq(a_NEq_d)

    ! Array of dependent variable initial values.
    integer         , intent(in   ) :: a_ODE_d
    ! Values used to normalize the particle mass and species mass in particles
    double precision, intent(in   ) :: a_ODE_norm(:)
    double precision, intent(  out) :: a_ODE_Vars(a_ODE_d)

    if (dem_solids .or. pic_solids) then
      call mapMFIXtoODE_dpm(a_NEq_d, a_NEq, a_ODE_d, a_ODE_norm, a_ODE_Vars)
    else
      call mapMFIXtoODE_tfm(a_NEq_d, a_NEq, a_ODE_d, a_ODE_norm, a_ODE_Vars)
    endif
  end subroutine mapMFIXtoODE

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Subroutine: mapODEtoMFIX                                          !
  !  Author: J. Musser                             Date: June 6, 2023  !
  !                                                                    !
  !  Purpose: This is a driver routine for mapping variables stored in !
  !  the ODE array back to MFIX field variables. Different routines    !
  !  are used for tfm and dpm runs. Fluid-only runs use the tfm        !
  !  interface.                                                        !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine mapODEtoMFIX (a_NEq_d, a_NEq, a_ODE_d, a_ODE_norm, a_ODE_Vars)

    use stiff_chem_dpm, only: mapODEtoMFIX_dpm
    use stiff_chem_tfm, only: mapODEtoMFIX_tfm

    use run, only: dem_solids, pic_solids

    implicit none

    ! Passed Variables:
    !-----------------------------------------------------------------//
    ! (1) Number of ODEs to be solve
    ! (2) Fluid cell index
    integer         , intent(in   ) :: a_NEq_d
    integer         , intent(in   ) :: a_NEq(a_NEq_d)

    ! Array of dependent variable initial values.
    integer         , intent(in   ) :: a_ODE_d
    ! Values used to normalize the particle mass and species mass in particles
    double precision, intent(in   ) :: a_ODE_norm(:)
    double precision, intent(  out) :: a_ODE_Vars(a_ODE_d)

    if (dem_solids .or. pic_solids) then
      call mapODEtoMFIX_dpm(a_NEq_d, a_NEq, a_ODE_d, a_ODE_norm, a_ODE_Vars)
    else
      call mapODEtoMFIX_tfm(a_NEq_d, a_NEq, a_ODE_d, a_ODE_norm, a_ODE_Vars)
    endif

  end subroutine mapODEtoMFIX

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Subroutine: reportODEVar                                          !
  !  Author: J. Musser                             Date: June 6, 2023  !
  !                                                                    !
  !  Purpose: This is a driver routine for checking variables stored   !
  !  in the ODE array. This is only used for debugging and is not      !
  !  not actively called by the stiff solver.                          !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine reportODEVar (a_NEq_d, a_NEq, a_ODE_d, a_ODE_norm, a_ODE_Vars)

    use stiff_chem_dpm, only: reportODEVar_dpm
    use stiff_chem_tfm, only: reportODEVar_tfm

    use run, only: dem_solids, pic_solids

    implicit none

    ! Passed variables:
    !-----------------------------------------------------------------//
    ! (1) Number of ODEs to be solve
    ! (2) Fluid cell index
    integer         , intent(in   ) :: a_NEq_d
    integer         , intent(in   ) :: a_NEq(a_NEq_d)

    ! Array of dependent variable initial values.
    integer         , intent(in   ) :: a_ODE_d
    ! Values used to normalize the particle mass and species mass in particles
    double precision, intent(in   ) :: a_ODE_norm(:)
    double precision, intent(in   ) :: a_ODE_Vars(a_ODE_d)

    if (dem_solids .or. pic_solids) then
      call reportODEVar_dpm(a_NEq_d, a_NEq, a_ODE_d, a_ODE_norm, a_ODE_Vars)
    else
      call reportODEVar_tfm(a_NEq_d, a_NEq, a_ODE_d, a_ODE_norm, a_ODE_Vars)
    endif
  end subroutine reportODEVar

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Subroutine: ClearStats                                            !
  !  Author: J. Musser                             Date: June 6, 2023  !
  !                                                                    !
  !  Purpose: Clear a few variables that we'll use to store solver     !
  !  information to report at the end of the routine.                  !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine clearStats ( )

    use error_manager, only: err_msg, log_message, loglevel_status
    use time_mod, only: get_elapsed_time
    implicit none

    m_startTime = get_elapsed_time()

    write(err_msg,"(/3x,'Integrating stiff chemistry...')")
    call log_status()

    m_failed = 0
    m_successful = 0
    m_incomplete = 0
    m_hetrgns = 0
    m_homogns = 0

    m_minTolSF =  1.0d8
    m_maxTolSF = -1.0d8

    m_minSteps =  100000000
    m_maxSteps = -100000000

    if (m_report == 0) then
      m_binnedStepCounts = 0
      m_binnedSolveTimes = 0
      m_report = m_reportFreq
    endif

  end subroutine clearStats

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  !                                                                    !
  !  Subroutine: ReportStats                                           !
  !  Author: J. Musser                             Date: June 6, 2023  !
  !                                                                    !
  !  Purpose: Print out solver information to give the user an idea    !
  !  of how the solver is performing.                                  !
  !                                                                    !
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine reportStats ( )

    use compar, only: PE_IO, myPE
    use error_manager, only: err_msg, log_message, loglevel_status
    use mpi_utility, only: global_sum, global_min, global_max
    use time_mod, only: get_elapsed_time

    implicit none

    ! Message buffer.
    character(len=64) :: lMsg0, lMsg1, lMsg2

    double precision :: solveTime

    ! array to collect global data
    integer :: int_global(3 + 2 + 5 + 7)
    integer :: int_local( 3 + 2 + 5 + 7)

    double precision :: dbl_gmax(3), dbl_gmin(3)
    double precision :: dbl_lmax(3), dbl_lmin(3)

    m_report = m_report - 1

    int_local(1) = m_successful
    int_local(2) = m_incomplete
    int_local(3) = m_failed

    int_local(4) = m_homogns
    int_local(5) = m_hetrgns

    int_local( 6:10) = m_binnedStepCounts(1:5)
    int_local(11:17) = m_binnedSolveTimes(1:7)

    call global_sum(int_local, int_global)

    solveTime = get_elapsed_time() - m_startTime

    dbl_lmin(1) = dble(m_minSteps)
    dbl_lmin(2) = m_minTolSF
    dbl_lmin(3) = solveTime

    call global_min(dbl_lmin, dbl_gmin)

    dbl_lmax(1) = dble(m_maxSteps)
    dbl_lmax(2) = m_maxTolSF
    dbl_lmax(3) = solveTime
    call global_max(dbl_lmax, dbl_gmax)

    ! Only the IO processor needs to continue.
    ! -!- NO MPI CALLS AFTER THIS POINT -!-
    if(myPE .ne. PE_IO) return

    if ( dbl_gmin(1) ==  100000000 .and. &
      &  dbl_gmax(1) == -100000000) then

      write(err_msg,"(A)") "No reaction information to report."
      call log_status()

      return
    endif

    ! Report Min/Max steps:
    lMsg0=''; write(lMsg0,*) int(dbl_gmin(1))
    lMsg1=''; write(lMsg1,*) int(dbl_gmax(1))
    write(err_msg,1000)  trim(adjustl(lMsg0)), trim(adjustl(lMsg1))
    call log_status()

    ! Report Homogeneous/Heterogeneous reactions:
    lMsg0=''; write(lMsg0,*) int_global(4)
    lMsg1=''; write(lMsg1,*) int_global(5)
    write(err_msg,1001) trim(adjustl(lMsg0)), trim(adjustl(lMsg1))
    call log_status()

    lMsg0=''; write(lMsg0,*) int_global(1)
    lMsg1=''; write(lMsg1,*) int_global(2)
    lMsg2=''; write(lMsg2,*) int_global(3)
    write(ERR_MSG,1002) trim(adjustl(lMsg0)), &
      &                 trim(adjustl(lMsg1)), &
      &                 trim(adjustl(lMsg2))
    call log_status()

    solveTime = dbl_gmin(3)
    if (solveTime > 3.6d3) then
       lMsg0=''; write(lMsg0,"(f10.6,' hrs')") solveTime/3.6d3
    elseif(solveTime > 6.0d1) then
       lMsg0=''; write(lMsg0,"(f10.6,' min')") solveTime/6.0d1
    else
       lMsg0=''; write(lMsg0,"(f10.6,' sec')") solveTime
    endif

    solveTime = dbl_gmax(3)
    if (solveTime > 3.6d3) then
       lMsg1=''; write(lMsg1,"(f10.6,' hrs')") solveTime/3.6d3
    elseif(solveTime > 6.0d1) then
       lMsg1=''; write(lMsg1,"(f10.6,' min')") solveTime/6.0d1
    else
       lMsg1=''; write(lMsg1,"(f10.6,' sec')") solveTime
    endif

    write(err_msg,1003) trim(adjustl(lMsg0)), trim(adjustl(lMsg1))
    call log_status()

    ! Report Homogeneous/Heterogeneous reactions:
    lMsg0=''; write(lMsg0,"(f6.4)") dbl_gmin(2)
    lMsg1=''; write(lMsg1,"(f6.4)") dbl_gmax(2)
    write(err_msg,1004) trim(adjustl(lMsg0)), trim(adjustl(lMsg1))
    call log_status()

    if (m_report <= 0) then

      write(err_msg,"(5x,'Average integration distribution:')")
      call log_status()
      write(err_msg,1100) "NST < 10^1", int_global( 6)/m_reportFreq
      call log_status()
      write(err_msg,1100) "NST < 10^2", int_global( 7)/m_reportFreq
      call log_status()
      write(err_msg,1100) "NST < 10^3", int_global( 8)/m_reportFreq
      call log_status()
      write(err_msg,1100) "NST < 10^4", int_global( 9)/m_reportFreq
      call log_status()
      write(err_msg,1100) "NST > 10^4", int_global(10)/m_reportFreq
      call log_status()

      write(err_msg,"(5x,'Average time per solve:')")
      call log_status()
      write(err_msg,1100)"time < 10^-3 sec", int_global(11)/m_reportFreq
      call log_status()
      write(err_msg,1100)"time < 10^-2 sec", int_global(12)/m_reportFreq
      call log_status()
      write(err_msg,1100)"time < 10^-1 sec", int_global(13)/m_reportFreq
      call log_status()
      write(err_msg,1100)"time < 10^+0 sec", int_global(14)/m_reportFreq
      call log_status()
      write(err_msg,1100)"time < 10^+1 sec", int_global(15)/m_reportFreq
      call log_status()
      write(err_msg,1100)"time < 10^+2 sec", int_global(16)/m_reportFreq
      call log_status()
      write(err_msg,1100)"time > 10^+2 sec", int_global(17)/m_reportFreq
      call log_status()

    endif

 1000 Format(5x,'Minimum/maximum number of steps over all cells: ',A,'/',A)
 1001 Format(5x,'Number of cells with Homogeneous/Heterogeneous reactions: ',A,'/',A)
 1002 Format(5x,'Number of successful/incomplete/failed integrations: ',&
             A,'/',A,'/',A)
 1003 Format(5x,'Minimum/maximum time Used: ',A,' / ',A)
 1004 Format(5x,'Minimum/maximum tolarance scale factor: ',A,'/',A)
 1100 Format(7x,A,': ',I6)
  end subroutine reportStats

end module StiffChemClass
