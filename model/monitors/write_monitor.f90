module write_monitor_mod

   use calc_cell_mod, only: calc_cell

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Author: J.Musser                                     Date: 11-2017  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  subroutine write_monitor(id, tavg)

    use monitor, only: monitor_type, monitor_var_count
    use monitor, only: monitor_x_w, monitor_x_e
    use monitor, only: monitor_y_s, monitor_y_n
    use monitor, only: monitor_z_b, monitor_z_t

    use monitor, only: monitor_i_w, monitor_i_e
    use monitor, only: monitor_j_s, monitor_j_n
    use monitor, only: monitor_k_b, monitor_k_t

    use monitor, only: monitor_time, monitor_vars
    use monitor, only: monitor_dt, monitor_dvars_odt
    use monitor, only: monitor_tavg, monitor_tavg_vars
    use monitor, only: monitor_tavg_dvars_odt

    use geometry, only: x_min, imax, dx
    use geometry, only: y_min, jmax, dy
    use geometry, only: z_min, kmax, dz, do_k

    use compar, only: mype,pe_io

    use monitor_functions, only: monitor_open_file
    use monitor_functions, only: monitor_write_values
    use monitor_functions, only: monitor_write_header
    use monitor_functions, only: monitor_calc_area
    use monitor_functions, only: monitor_calc_volume

    use monitor_average_mod,       only: monitor_average
    use monitor_minmax_mod,        only: monitor_minmax
    use monitor_area_avg_mod,      only: monitor_area_avg
    use monitor_volume_avg_mod,    only: monitor_volume_avg

    use monitor_average_des_mod,   only: monitor_average_des
    use monitor_minmax_des_mod,    only: monitor_minmax_des
    use monitor_flowrate_des_mod,  only: monitor_flowrate_des

    use monitor, only: monitor_tavg_dt, monitor_tavg_count, monitor_tavg_vars
    use monitor, only: monitor_tavg_start, monitor_tavg_reset_dt
    use monitor, only: monitor_tavg_ss_tol, monitor_tavg_ss_duration
    use monitor, only: monitor_tavg_ss_span, monitor_ss
    use monitor, only: monitor_header_written

    USE mpi_utility, only: bcast
    use param1, only: undefined, zero
    use run, only: time

    implicit none


    integer, intent(in   ) :: id
    logical, intent(in   ) :: tavg  ! set tavg to .true. to write time-averaged data
                                    ! set tavg to .false. to write instantaneous data

    integer          :: i_w, i_e, j_s, j_n, k_b, k_t
    double precision :: x_w, x_e, y_s, y_n, z_b, z_t

    integer                       :: funit, funit_tavg
    double precision, allocatable :: data(:)


    integer :: var, ss_count
    double precision :: n  ! Warning!! Here n is a float so it can handle large number of samples
    double precision :: previous  ! Used to compute rate of change of time-averaged data

! When time-averaging, quit if sampling interval is not defined or the
! averaging start time has not been reached

    if(tavg.and.monitor_tavg_dt(id) == undefined) return
    if(tavg.and.time<monitor_tavg_start(id)) return

! Monitor cell indices
    i_w = monitor_i_w(id)
    i_e = monitor_i_e(id)
    j_s = monitor_j_s(id)
    j_n = monitor_j_n(id)
    k_b = monitor_k_b(id)
    k_t = monitor_k_t(id)

    if(monitor_var_count(id) <= 0) then
       return ! Error
    else
       allocate(data(monitor_var_count(id)))
    endif

    if(.not.monitor_header_written(id)) then
       call monitor_calc_area(id)
       call monitor_calc_volume(id)
    endif

    if(mype == pe_io) then
       if(.not.tavg) then
          call monitor_open_file(id, funit,.false.)
          call monitor_write_header(id, funit)
          if (monitor_tavg(id)) then
             call monitor_open_file(id, funit_tavg,.true.)
             call monitor_write_header(id, funit_tavg)
          endif
       endif
       monitor_header_written(id) = .true.
    endif

    call bcast(monitor_header_written(id))

    select case(monitor_type(id))

    case(0,1,4,5)
       !  0: Point value
       !  1: Sum over region
       !  4: Arithmetic average over region
       !  5: Standard deviation over region
       call monitor_average(id, data, i_w, i_e, j_s, j_n, k_b, k_t)

    case(2,3)
       !  2: Minimum over region
       !  3: Maximum over region
       call monitor_minmax(id, data, i_w, i_e, j_s, j_n, k_b, k_t)

    case(6,7,8,9,10)
       !  6: Area-weighted average over surface
       !  7: Flow rate across surface
       !  8: Mass flow rate across surface
       !  9: Mass-weighted average over surface
       ! 10: Volumetric flow rate over surface
       call monitor_area_avg(id, data, i_w, i_e, j_s, j_n, k_b, k_t)

    case(11,12,13,14)
       ! 11: Volume integral
       ! 12: Volume-weighted average
       ! 13: Mass-weighted volume integral
       ! 14: Mass-weighted volume average
       call monitor_volume_avg(id, data, i_w, i_e, j_s, j_n, k_b, k_t)

    case(101,104,105,106,107)
       ! 101: Sum over region (DEM/PIC)
       ! 104: Arithmetic average over region (DEM/PIC)
       ! 105: Standard deviation over region (DEM/PIC)
       ! 106: Mass-weighted average over region (DEM/PIC)
       ! 107: Volume-weighted average over region (DEM/PIC)
       call monitor_average_des(id, data, i_w, i_e, j_s, j_n, k_b, k_t)

    case(102,103)
       ! 102: Minimum over region (DEM/PIC)
       ! 103: Maximum over region (DEM/PIC)
       call monitor_minmax_des(id, data, i_w, i_e, j_s, j_n, k_b, k_t)

    case(108,109,110)
       ! 108: Flow rate (DEM/PIC)
       ! 109: Mass-weighted flow rate (DEM/PIC)
       ! 110: Volume-weighted flow rate (DEM/PIC)
       call monitor_flowrate_des(id, data)

    end select


    if(mype == pe_io) then
! Time-average
! if monitor_tavg_reset_dt is defined, the time average is reset every
! monitor_tavg_reset_dt (similar to a running average)
!
       if(tavg) then

          n = monitor_tavg_count(id)
          monitor_tavg_count(id) = n + 1

          if(n==0) then
             monitor_tavg_vars(id, 1:monitor_var_count(id)) = data(:)
          else
             ss_count = 0
             do var = 1, monitor_var_count(id)
                previous = monitor_tavg_vars(id, var)
! New values
                monitor_tavg_vars(id, var) = ( monitor_tavg_vars(id, var) * n + data(var) ) / (n + 1)
! Rate of change
                monitor_tavg_dvars_odt(id, var) = (monitor_tavg_vars(id, var) - previous) / monitor_tavg_dt(id)

! Check if steady state was reached
                if(dabs(monitor_tavg_dvars_odt(id, var))<monitor_tavg_ss_tol(id)) then
                   monitor_tavg_ss_duration(id,var) = monitor_tavg_ss_duration(id,var) + monitor_tavg_dt(id)
                else
                   monitor_tavg_ss_duration(id,var) = zero
                endif
                if(monitor_tavg_ss_duration(id,var)>monitor_tavg_ss_span(id)) ss_count = ss_count + 1
             enddo
             if(ss_count==monitor_var_count(id)) monitor_ss(id) = .true.

          endif

          if(monitor_tavg_reset_dt(id) /= undefined) then
             if(monitor_tavg_count(id)>monitor_tavg_reset_dt(id)/monitor_tavg_dt(id) - 1.0)   monitor_tavg_count(id) = 0.0
          endif

       else
          monitor_time(id) = time
! Rate of change
          monitor_dvars_odt(id,1:monitor_var_count(id)) = (data(1:monitor_var_count(id)) - monitor_vars(id,1:monitor_var_count(id)) )/monitor_dt(id)
! New values
          monitor_vars(id,1:monitor_var_count(id)) = data(1:monitor_var_count(id))
          call monitor_write_values(funit, data)
          close(funit)
          if (monitor_tavg(id)) then
             if(time>=monitor_tavg_start(id)) call monitor_write_values(funit_tavg, monitor_tavg_vars(id,1:monitor_var_count(id)))
             close(funit_tavg)
          endif
       endif

    endif


    if(allocated(data)) deallocate(data)

    return

  end subroutine write_monitor



end module write_monitor_mod
