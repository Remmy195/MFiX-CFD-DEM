MODULE monitor

  Use param1
! Maximum number of solids phases.
  use param, only: DIM_M
! Maximum number of gas phase species
  use param, only: DIM_N_g
! Maximum number of solids phase species
  use param, only: DIM_N_s
! Maximum number of scalar equations
  use param, only: DIM_Scalar
! Maximum number of DEM solids phase species
  use param, only: DIMENSION_N_S
! Maximum number of reaction rates saved in monitors
  use rxns, only: nrrmax

  integer, parameter :: dimension_monitor = 100

  character(len=256) :: monitor_name(dimension_monitor)
  integer :: monitor_type(dimension_monitor)
  double precision :: monitor_dt(dimension_monitor)

  logical :: monitor_defined (dimension_monitor)
  double precision :: monitor_time(dimension_monitor)
  double precision :: monitor_next_time(dimension_monitor)
  integer :: monitor_var_count(dimension_monitor)

! Spatial location of monitors
  double precision :: monitor_x_w(dimension_monitor)
  double precision :: monitor_x_e(dimension_monitor)
  double precision :: monitor_y_s(dimension_monitor)
  double precision :: monitor_y_n(dimension_monitor)
  double precision :: monitor_z_b(dimension_monitor)
  double precision :: monitor_z_t(dimension_monitor)

! Indices location of monitors
  integer :: monitor_i_w(dimension_monitor)
  integer :: monitor_i_e(dimension_monitor)
  integer :: monitor_j_s(dimension_monitor)
  integer :: monitor_j_n(dimension_monitor)
  integer :: monitor_k_b(dimension_monitor)
  integer :: monitor_k_t(dimension_monitor)

! Geometric type (point, plane, volume)
  character(len=2) :: monitor_geotype(dimension_monitor)

! Area and volume
  double precision :: monitor_area(dimension_monitor)
  double precision :: monitor_volume(dimension_monitor)

! Flag to know if header was written or not
  logical :: monitor_header_written(dimension_monitor) = .false.

! Fluid phase variables
  logical :: monitor_ep_g(dimension_monitor)
  logical :: monitor_ro_g(dimension_monitor)
  logical :: monitor_p_g(dimension_monitor)
  logical :: monitor_u_g(dimension_monitor)
  logical :: monitor_v_g(dimension_monitor)
  logical :: monitor_w_g(dimension_monitor)
  logical :: monitor_t_g(dimension_monitor)
  logical :: monitor_mw_mix_g(dimension_monitor)
  logical :: monitor_x_g(dimension_monitor, dim_n_g)
  logical :: monitor_y_g(dimension_monitor, dim_n_g)

  logical :: monitor_k_turb_g(dimension_monitor)
  logical :: monitor_e_turb_g(dimension_monitor)

! Eulerian solids variables
  logical :: monitor_ep_s (dimension_monitor, dim_m)
  logical :: monitor_rop_s(dimension_monitor, dim_m)
  logical :: monitor_ro_s(dimension_monitor, dim_m)
  logical :: monitor_p_star(dimension_monitor)
  logical :: monitor_p_s(dimension_monitor, dim_m)
  logical :: monitor_u_s(dimension_monitor, dim_m)
  logical :: monitor_v_s(dimension_monitor, dim_m)
  logical :: monitor_w_s(dimension_monitor, dim_m)
  logical :: monitor_t_s(dimension_monitor, dim_m)
  logical :: monitor_x_s(dimension_monitor, dim_m, dim_n_s)
  logical :: monitor_theta_m(dimension_monitor, dim_m)

  ! Other Field variables
  logical :: monitor_scalar(dimension_monitor, dim_scalar)
  logical :: monitor_rrate(dimension_monitor,nrrmax)
  logical :: monitor_fluid_rrate(dimension_monitor,nrrmax)
  logical :: monitor_des_rrate(dimension_monitor,nrrmax)

  ! DES specific variables
  logical :: monitor_part_phase(dimension_monitor, dim_m)

  logical :: monitor_radius(dimension_monitor)
  logical :: monitor_pmass(dimension_monitor)
  logical :: monitor_pcount(dimension_monitor)
  logical :: monitor_pvol(dimension_monitor)
  logical :: monitor_ro_p(dimension_monitor)

  logical :: monitor_pos_x(dimension_monitor)
  logical :: monitor_pos_y(dimension_monitor)
  logical :: monitor_pos_z(dimension_monitor)

  logical :: monitor_vel_x(dimension_monitor)
  logical :: monitor_vel_y(dimension_monitor)
  logical :: monitor_vel_z(dimension_monitor)

  logical :: monitor_rot_x(dimension_monitor)
  logical :: monitor_rot_y(dimension_monitor)
  logical :: monitor_rot_z(dimension_monitor)

  logical :: monitor_t_p(dimension_monitor)
  logical :: monitor_x_p(dimension_monitor, dim_n_s)

  logical :: monitor_des_usr_var(dimension_monitor, 100)
  logical :: monitor_part_rrate(dimension_monitor, 2000)

! Residence time
  logical :: monitor_part_residence_time(dimension_monitor)

! Monitor type name
  character(len=64) :: monitor_typename(dimension_monitor)
! Monitor variable name list
  character(len=128) :: monitor_varname_list(dimension_monitor,256)

! Monitor variable values
  double precision :: monitor_vars(dimension_monitor,256)
  double precision :: monitor_dvars_odt(dimension_monitor,256)

! Monitor time-averaged data
  double precision :: monitor_tavg_dt(dimension_monitor)
  double precision :: monitor_tavg_start(dimension_monitor)
  double precision :: monitor_tavg_reset_dt(dimension_monitor)
  double precision :: monitor_tavg_ss_tol(dimension_monitor)
  double precision :: monitor_tavg_ss_span(dimension_monitor)

  logical :: monitor_tavg(dimension_monitor) = .false.
  double precision :: monitor_next_tavg_time(dimension_monitor)
  double precision :: monitor_tavg_next_time(dimension_monitor)
  double precision :: monitor_tavg_count(dimension_monitor)
  double precision :: monitor_tavg_vars(dimension_monitor,256)
  double precision :: monitor_tavg_dvars_odt(dimension_monitor,256)
  double precision :: monitor_tavg_ss_duration(dimension_monitor,256)=0.0D0
  logical :: monitor_ss(dimension_monitor) = .false.
  logical :: monitor_ss_exit = .false.

end module monitor
