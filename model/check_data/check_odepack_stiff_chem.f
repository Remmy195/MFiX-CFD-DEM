#include "error.inc"

module check_odepack_stiff_chem_mod

  use error_manager

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_DATA_CHEM                                        !
!  Author: J.Musser                                   Date: 02-Aug-13  !
!                                                                      !
!  Purpose: Check the chemical rxns namelist variables.                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine check_odepack_stiff_chem

    ! Global Variables:
    !-----------------------------------------------------------------//
    ! Skip data check when doing preprocessing only
    use run, only: ppo
    ! Constant gas phase density
    use physprop, only : RO_G0
    ! Runtime logical for solving energy and species equations.
    use run, only: energy_eq, species_eq
    ! Net rate of gas/solids phase production/consumption
    use rxns, only: SUM_R_g, SUM_R_s
    ! Run time logical for using stiff chemistry solver
    use stiff_chem, only: stiff_chemistry
    use stiff_chem, only: stiff_chem_max_steps


    ! Double precision zero.
    use param1, only: zero, undefined, undefined_i

    implicit none

    if(ppo) return

    ! Verify that there is sufficient run complexity to use the stiff solver
    if(stiff_chemistry) then

    ! Energy equations must be solved.
     if(.not.energy_eq) then
        write(err_msg,1004)'ENERGY_EQ = .FALSE.'
        call log_error()
     endif

     ! Cannot have constant fluid density
     if(ro_g0 /= undefined) then
       write(err_msg,1004)'RO_G0 /= UNDEFINED'
       call log_error()
     endif

     ! Must solve fluid species equations
     if(.not.species_eq(0)) then
        write(err_msg,1004)'SPECIES_EQ(0) = .FALSE.'
        call log_error()
     endif

     if(stiff_chem_max_steps == undefined_i) then
        stiff_chem_max_steps = 500000
     endif

! Clear the interphase mass transfer terms as the stiff solver
! does no use them.
     if(allocated(sum_r_g)) SUM_R_g = zero
     if(allocated(sum_r_s)) SUM_R_s = zero

   endif

 1004 FORMAT('Error 1004: ',                                           &
         'Invalid parameters for stiff chemistry solver!',//           &
         ' The following criteria must be satisfied:',/                &
         '   > Solving the energy equations.',/                        &
         '   > Compressible gas phase.',/                              &
         '   > Solving gas phase species equations.',//                &
         ' >>> Invalid Parameter: ',A,//                               &
         'Check the user documentation and correct the project settings.')

  return
  end subroutine check_odepack_stiff_chem

end module check_odepack_stiff_chem_mod
