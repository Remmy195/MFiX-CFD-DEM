#include "error.inc"

! FIXME main.f says this writes to the .OUT file but it does not.

module write_out3_mod
contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: write_out3                                             !
!  Author: M. Syamlal                                 Date: 10-Jan-92  !
!                                                                      !
!  Purpose: to write time used by the simulation                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine write_out3(elapsed_time, io_time, paused_time)

      use error_manager
      use run, only: time, tstart, get_tunit

      implicit none

      double precision, intent(in) :: elapsed_time
      double precision, intent(in) :: io_time
      double precision, intent(in) :: paused_time

      double precision:: time1
      double precision:: tstart1
      double precision:: elapsed_time1
      double precision:: io_time1
      double precision:: paused_time1

      character(len=4) :: unit

! copy values since get_tuinit modifies them and arguments are intent(in)
      tstart1 = tstart
      time1 = time
      elapsed_time1 = elapsed_time
      io_time1 = io_time
      paused_time1 = paused_time

      call get_tunit(tstart1, unit)
      write(err_msg, 1000) "Simulation start time", trim(ival(tstart1)), unit
      call log_status()

      call get_tunit(time1, unit)
      write(err_msg, 1000) "Simulation time reached", trim(ival(time)), unit
      call log_status()

      call get_tunit(elapsed_time1, unit)
      write(err_msg, 1000) 'Elapsed real time', trim(ival(elapsed_time1)), unit
      call log_status()

      call get_tunit(io_time1, unit)
      write(err_msg, 1000) 'Time spent in I/O', trim(ival(io_time1)), unit
      call log_status()

      if (paused_time1 > 0) then
          call get_tunit(paused_time1, unit)
          write(err_msg, 1000) 'Time spent paused', trim(ival(paused_time1)), unit
          call log_status()
      end if

 1000 format(a,' = ',a, a)

      return

   end subroutine write_out3

end module write_out3_mod
