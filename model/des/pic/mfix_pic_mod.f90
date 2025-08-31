!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: MFIX_PIC                                                    !
!  Purpose: MP-PIC related data                                        !
!  Author: R. Garg                                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
MODULE MFIX_PIC

  use param, only: dim_m
  use param, only: dimension_bc
!......................................................................!


! User inputs: .........................................................


  double precision :: psfac_fric_pic
  double precision :: fric_exp_pic
  double precision :: fric_non_sing_fac

  double precision :: mppic_coeff_en1
  double precision :: mppic_coeff_en_wall
  double precision :: mppic_coeff_et_wall
  double precision :: mppic_velfac_coeff

  double precision :: pic_cfl
  double precision :: pic_cfl_parcel_fraction
  character(3) :: pic_cfl_control


! Run time variables ....................................................


  ! Flag for turning on MP-PIC method
  logical :: mppic

  ! grid dimension minimums
  double precision :: dxyz_min

  ! Count and map between PIC MI IDs and user BC index.
  integer :: pic_bcmi, pic_bcmi_map(dimension_bc)

  ! Mass inflow data structure.
  integer, parameter:: dim_pic_mi_bin=1001 ! maximum bin size (used for custom distribution)
  integer, parameter:: pic_numbin=100      ! bin size for Normal and Log-Normal distributions
                                           ! Note: It would be better to
                                           ! dynamically allocate arrays but
                                           ! some compilers have issues with
                                           ! allocatable arrays in derived type
                                           ! structures.
                                           ! Must have pic_numbin <  dim_pic_mi_bin

  type pic_mi_
     double precision :: pflow(dim_m)
     double precision :: remdr(dim_m)

     double precision :: offset(3)
     double precision :: length(3)

     logical, allocatable :: blocked(:,:,:)

  !incorporate polydispersity  (set at dim_pic_mi_bin to accommodate custom
  !                             distribution limits; normal and log-normal
  !                             distributions require pic_numbin+1)
     double precision :: polypflow(dim_pic_mi_bin,dim_m)
     double precision :: median_dia(dim_pic_mi_bin,dim_m)
     double precision :: polyremdr(dim_pic_mi_bin,dim_m)
     integer :: NUMBIN(dim_m)

  end type pic_mi_

  ! Mass inflow data.
  type(pic_mi_), target, allocatable :: pic_mi(:)

  ! Mass inflow maximum skipped seeding
  integer :: pic_mi_max_skip = 1000

  ! Parcel Statistical weight
  double precision, allocatable :: des_stat_wt(:)

  ! solid pressure gradient interpolated to parcel location
  double precision, allocatable :: ps_grad(:,:)

  ! Volume fraction interpolated to parcel location (not used)
  double precision, allocatable :: epg_p(:)

  ! Flag for turning on collision damping
  logical :: PIC_COLLISION_DAMPING

  ! Restitution coefficient used with collision damping
  double precision :: PIC_CD_E
end module mfix_pic
