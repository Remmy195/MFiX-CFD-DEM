module time_mod
!! All times are wall time.
!! Do not use CPU time for anything. See https://mfix.netl.doe.gov/gitlab/develop/mfix/-/issues/48

   use run, only: time, tstop, tstart

   implicit none

   double precision:: paused_time =  0.0d0
   double precision:: io_time = 0.0d0
   double precision:: pause_start_time = 0.0d0
   double precision:: io_start_time = 0.0d0
   integer:: io_depth = 0      ! start_io/end_io calls may be nested
   integer:: pause_depth = 0   ! start_pause/end_pause calls may be nested

   interface
      double precision function call_etime() bind(C,name="etime")
          use, intrinsic :: iso_c_binding, only : c_double
          real(kind=c_double):: etime
      end function call_etime
   end interface

contains

    double precision function get_elapsed_time()
        get_elapsed_time = call_etime()
    end function get_elapsed_time

    double precision function get_paused_time()
        if (pause_depth > 0) then
            get_paused_time = paused_time + get_elapsed_time() - pause_start_time
        else
            get_paused_time = paused_time
        endif
    end function get_paused_time

    double precision function get_io_time()
        if (io_depth > 0) then
            get_io_time = io_time + get_elapsed_time() - io_start_time
        else
            get_io_time = io_time
        endif
    end function get_io_time

    double precision function get_remaining_time()
        get_remaining_time = 0.0d0
        if (time > tstart .and. time <= tstop) then
            get_remaining_time = (get_elapsed_time() - get_paused_time()) &
                                   *  (tstop-time)/(time-tstart)
        end if
    end function get_remaining_time

    subroutine start_io
        if (io_depth == 0) then
            io_start_time = get_elapsed_time()
        end if
        io_depth = io_depth + 1
    end subroutine start_io

    subroutine end_io
        if (io_depth <= 0) then
            write(*,*) "Error, end_io called before start_io"
            error stop
        end if
        io_depth = io_depth - 1
        if (io_depth == 0) then
            io_time = io_time + get_elapsed_time() - io_start_time
        end if
    end subroutine end_io

    subroutine start_pause
        if (pause_depth == 0) then
            pause_start_time = get_elapsed_time()
        end if
        pause_depth = pause_depth + 1
    end subroutine start_pause

    subroutine end_pause
        if (pause_depth <= 0) then
            write(*,*) "Error, end_io called before start_io"
            error stop
        end if
        pause_depth = pause_depth - 1
        if (pause_depth == 0) then
            paused_time = paused_time + get_elapsed_time() - pause_start_time
        end if
    end subroutine end_pause

end module time_mod
