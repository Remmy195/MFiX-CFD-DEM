module pause

    use exit, only : exit_flag, check_exit_flag
    use compar, only : mype
    use mpi_utility
    use reinit, only: reinitialize
    use error_manager, only: reinit_error
    use reset_new_mod, only: reset_new
    use time_mod, only: start_pause, end_pause
    use, intrinsic:: iso_c_binding, only: c_int

    implicit none

    interface
       subroutine usleep(us) bind (C)
           import c_int
           integer(c_int), value :: us
       end subroutine usleep
    end interface


    logical :: pause_flag = .false.
    logical :: reinit_flag = .false.
    logical :: autostart_flag = .false.  ! Automatically restart after successful reinit
    integer :: ier = 0

! TODO dynamically allocate this, instead of fixed-length
    character(100000):: reinit_data

contains
    subroutine wait_while_paused
#ifdef PYMFIX
        double precision :: wall_paused_prev

        if (check_exit_flag()) return

        if (paused_flag_bcast()) then
            call start_pause
            if (mype .eq. 0) then
                print *, "Paused"
!           else
!               print *, "Paused", mype
            end if

            do while (paused_flag_bcast() .and. .not. check_exit_flag())
                if (check_reinit_flag()) then
                    call do_reinit
                    reinit_flag = .false.
                else
                    call usleep(100000) ! 10 hz
!                   call sleep(1) ! 1 hz
                end if
            end do

            if (check_reinit_flag()) then
                call do_reinit
                reinit_flag = .false.
            end if

            if (mype .eq. 0) then
                print *, "resuming"
!           else
!               print *, "resuming", mype
            end if
            call end_pause

        end if
#endif
    end subroutine wait_while_paused

#ifdef PYMFIX
    logical function paused_flag_bcast()
        logical :: local
        local = pause_flag ! Assuming atomic assignment of boolean
        call bcast(local)
        paused_flag_bcast = local
    end function paused_flag_bcast


    logical function check_reinit_flag()
        logical :: local
        local = reinit_flag ! Assuming atomic assignment of boolean
        call bcast(local)
        check_reinit_flag = local

        if(check_reinit_flag) call bcast(reinit_data)  ! Should use mutex for here and set_reinit_data

    end function check_reinit_flag


    subroutine set_reinit_data(data, length)
        integer, intent(in):: length
        character(len=length), intent(in):: data
        reinit_data = data
    end subroutine set_reinit_data


    subroutine do_reinit

        if (mype .eq. 0) then
            print *, "Reinitializing"
        end if

        call bcast(reinit_data)
        call reset_new(reinit_data) ! we might be paused in the middle of a time step

        call reinitialize(reinit_data, ier)

        if (check_exit_flag() .or. reinit_error()) then
            if (mype .eq. 0) then
                print *, "Reinitialization failed"
            end if
            exit_flag = .false.  !! allow user to retry.
        else
            if (autostart_flag) then
                autostart_flag = .false.
                pause_flag = .false.
            end if
        end if

    end subroutine do_reinit


#endif
end module pause
