#include "error.inc"

MODULE READ_PART_INPUT_MOD

   USE error_manager
   USE read_database_mod, only: read_database

CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: READ_PART_INPUT                                           !
!                                                                      !
! Revised: Renjie Ke                                  Date: July, 2024 !
! Purpose: Read the particle input and broadcasts the particle data to !
! respective processors.                                               !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_PART_INPUT
!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      use cdist
      use compar
      use desmpi
      use functions
      use funits
      use geometry, only: NO_K
      use mpi_init_des, only: des_scatter_particle
      use mpi_utility

      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: k
! index of particle
      INTEGER :: lcurpar
! local unit
      INTEGER, PARAMETER :: lunit=10
! local filename
      character(255) :: lfilename, lfilename_dat, lfilename_csv
! IO Status:
      INTEGER :: IOS
! Flag to indicate if file exists.
      LOGICAL :: lEXISTS_csv, lEXISTS_dat
! Read dimension: 2D vs 3D data
      integer :: RDMN
! local filename
      character(11) Version_string
!-----------------------------------------------

      IOS = 0
      RDMN = merge(2,3,NO_K)

! Setup the file name based on distributed or serial IO.
      IF(bDIST_IO) THEN
         lFILENAME = ''
         WRITE(lFILENAME,'("particle_input_",I4.4,".dat")') myPE
      ELSE
         if(mype==pe_io) WRITE(*,*) "Loading particle_input.csv ..."
         lFILENAME_csv = "particle_input.csv"
      ENDIF

      IF(myPE == PE_IO) THEN
         INQUIRE(FILE=lFILENAME_csv, EXIST=lEXISTS_csv)
! First check if csv file exists.
         IF(lEXISTS_csv) THEN
            Version_string = 'Version 4.0'
         ELSE
! Second check if dat file exists.
            lfilename_dat= "particle_input.dat"
            INQUIRE(FILE=lfilename_dat, EXIST=lEXISTS_dat)

            IF(.NOT. lEXISTS_dat) THEN
! If BOTH csv and dat file are non-existent, then report error
               WRITE(ERR_MSG, 1100)
               CALL LOG_ERROR()
               IOS = 1
            ELSE
               OPEN(UNIT=lUNIT, FILE=lfilename_dat, FORM="FORMATTED")
               read(lUNIT,'(A11)') Version_string
               backspace(lUNIT)
            ENDIF
         ENDIF
      ENDIF

      call bcast(VERSION_STRING, PE_IO)

! Collect the error message and quit.
      CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) call LOG_ERROR()

 1100 FORMAT('FATAL - DEM particle input file not found.')

! Read the file
!----------------------------------------------------------------->>>
      SELECT CASE (Version_string)

         CASE('Version 4.0')
! Version 4.0 is a unified version, we expect using this to deal with ALL solid models
! TODO deprecate all older versions

            P_INPUT_DAT_VERSION = '4.0'
            WRITE(ERR_MSG, 1200) P_INPUT_DAT_VERSION
            CALL LOG_STATUS()
            CALL READ_PART_INPUT_V4P0

         CASE('Version 3.0')

            P_INPUT_DAT_VERSION = '3.0'
            WRITE(ERR_MSG, 1200) P_INPUT_DAT_VERSION
            CALL LOG_STATUS()
            CALL READ_PART_INPUT_V3P0

         CASE('Version 2.0')

            P_INPUT_DAT_VERSION = '2.0'
            WRITE(ERR_MSG, 1200) P_INPUT_DAT_VERSION
            CALL LOG_STATUS()
            IF(GluedsphereDEM) THEN
               CALL READ_PART_INPUT_GluedSphereDEM
            ELSE
               CALL READ_PART_INPUT_V2P0
            ENDIF

         CASE DEFAULT

            P_INPUT_DAT_VERSION = '1.0'
            IF(PARTICLES>0) THEN
               WRITE(ERR_MSG, 1200) P_INPUT_DAT_VERSION
               !CALL LOG_STATUS()
               CALL CHECK_PARTICLE_INPUT_FILTER_WITH_V1
               IF (SuperDEM) THEN
                   CALL READ_PART_INPUT_SuperDEM
               ELSE
                   CALL READ_PART_INPUT_V1P0
               ENDIF
            ELSE
               WRITE(ERR_MSG, 1210) P_INPUT_DAT_VERSION
               CALL LOG_INFO()
            ENDIF

      END SELECT
 1200 FORMAT('Info: Reading DEM particle input file version: ', A)
 1210 FORMAT('Info: PARTICLES is set to 0: Skipping reading DEM particle input file version: ', A)
      RETURN
      END SUBROUTINE READ_PART_INPUT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: READ_PART_INPUT_V2P0                                      !
!                                                                      !
! Purpose: Read the particle input and broadcasts the particle data to !
! respective processors.                                               !
! Version 2.0                                                          !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_PART_INPUT_V2P0

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      use cdist
      use compar
      use des_allocate, only: particle_grow
      use desmpi
      use desmpi_wrapper, only: DES_MPI_STOP
      use des_thermo, only: des_t_s
      use discretelement, only: des_usr_var_size
      use functions
      use funits
      use geometry, only: NO_K
      use mass_outflow_dem_mod, only: delete_particle
      use mpi_comm_des, only: desmpi_scatterv
      use mpi_init_des, only: des_scatter_particle
      use mpi_utility
      use parallel_mpi
      use param, only: dimension_n_s
      use param1, only: UNDEFINED_I
      use run, only: any_solids_species_eq, energy_eq
      use des_rxns, only: des_x_s

      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: i,k,ilines,m,n
! index of particle
      INTEGER :: lcurpar, LC1
! local unit
      INTEGER, PARAMETER :: lunit=10
! local filename
      character(255) lfilename
! IO Status:
      INTEGER :: IOS
! Flag to indicate if file exists.
      LOGICAL :: lEXISTS
! Read dimension: 2D vs 3D data
      integer :: RDIMN
! Buffer
      character(32) Buffer
! Option to read or use default value (T/F)
      character(1) option
! Integer buffer
      integer :: int_buffer
! Double precision buffer, can be a list of up to 100 numbers
      integer, dimension(100) :: dp_buffer
! Length of list of variables to read (number of columns)
      integer :: n_vars, n_cols
      integer :: XYZ_start
      double precision, dimension(3) :: Dummy_xyz
      integer :: Phase_ID_start
      double precision :: Uniform_Phase_ID
      integer :: Diameter_start
      double precision :: Uniform_Diameter
      integer :: Density_start
      double precision :: Uniform_Density
      integer :: Velocity_start
      double precision, dimension(3) :: Uniform_Velocity
      integer :: Temperature_start
      double precision :: Uniform_Temperature
      integer :: Species_start
      double precision, dimension(100) :: Uniform_Species
      integer :: User_scalar_start
      double precision, dimension(100) :: Uniform_User_scalar

      double precision, allocatable, dimension(:,:) :: part_data
      integer :: lproc,lbuf,lpacketsize
      integer :: lproc_parcnt(0:numpes-1)
      integer, allocatable, dimension(:) :: lpar_proc
      logical :: x_test, y_test, z_test, xyz_test
! Particle filter
      LOGICAL, allocatable, dimension(:) :: keep_particle
      LOGICAL :: keep_me
      character(32) :: label

      character(1024) :: first_line
      Integer :: CC1,CC2

      character(7) :: filter_type
      integer :: file_particles ! number of particles in the file
      Integer :: fpip ! Number of particles after filtering the data
!-----------------------------------------------

      n_vars = 0

      IOS = 0
      RDIMN = merge(2,3,NO_K)

      IF (myPE .eq. PE_IO) THEN
! Read Header
! Skip Instructions
         ilines = 18
         do i = 1, ilines
            read(lunit,*) buffer
         enddo

! Dimension
         read(lunit,*) buffer,int_buffer
         if(RDIMN/=int_buffer) then
            WRITE(ERR_MSG,100) int_buffer, RDIMN
            CALL LOG_ERROR()
         endif

! Number of particles
         read(lunit,*) buffer, file_particles

         if (particles .eq. UNDEFINED_I) particles = file_particles

! Skip header
         do i = 1, 3
            read(lunit,*) buffer
         enddo

! Coordinates option
         call parse_vector_option(lunit, 'Coordinates', .TRUE., rdimn, Dummy_xyz(1:rdimn), XYZ_start, n_cols, n_vars)

! Phase ID option
         call parse_scalar_option(lunit, 'Phase_ID', .TRUE., Uniform_Phase_ID, Phase_ID_start,  n_cols, n_vars)

! Diameter option
         call parse_scalar_option(lunit, 'Diameter', .TRUE., Uniform_Diameter, Diameter_start,  n_cols, n_vars)

! Density option
         call parse_scalar_option(lunit, 'Density', .TRUE., Uniform_Density, Density_start,  n_cols, n_vars)

! Velocity option
          call parse_vector_option(lunit, 'Velocity', .TRUE., rdimn, Uniform_Velocity(1:rdimn), Velocity_start, n_cols, n_vars)

! Temperature option
         call parse_scalar_option(lunit, 'Temperature', ENERGY_EQ, Uniform_Temperature, Temperature_start, n_cols, n_vars)

! Species option
          call parse_vector_option(lunit, 'Species', ANY_SOLIDS_SPECIES_EQ, DIMENSION_N_S, Uniform_Species(1:DIMENSION_N_S), Species_start, n_cols, n_vars)

! User Scalar option
          call parse_vector_option(lunit, 'User scalar', (DES_USR_VAR_SIZE>0), DES_USR_VAR_SIZE, Uniform_User_scalar(1:DES_USR_VAR_SIZE), User_scalar_start, n_cols, n_vars)


 100 FORMAT('Error reading particle input file.',&
           /'Dimension (2D or 3D) does not match simulation setup file.',&
           /'Dimension in particle_input.dat    = ',I3,&
           /'Dimension in simulation setup file = ',I3)


! Skip the next three lines. This should get us to line 35 where the data
! starts
! Skip header
         do i = 1, 3
            read(lunit,*) buffer
         enddo


      ENDIF !  IF (myPE .eq. PE_IO)

! Broadcast column indices and default values to all PEs
      call bcast(n_vars,pe_io)

      call bcast(XYZ_start,pe_io)
      call bcast(Phase_ID_start,pe_io)
      call bcast(Uniform_Phase_ID,pe_io)
      call bcast(Diameter_start,pe_io)
      call bcast(Uniform_Diameter,pe_io)
      call bcast(Density_start,pe_io)
      call bcast(Uniform_Density,pe_io)
      call bcast(Velocity_start,pe_io)
      call bcast(Uniform_Velocity,pe_io)
      call bcast(Temperature_start,pe_io)
      call bcast(Uniform_Temperature,pe_io)
      call bcast(Species_start,pe_io)
      call bcast(Uniform_Species,pe_io)
      call bcast(User_scalar_start,pe_io)
      call bcast(Uniform_User_scalar,pe_io)

! Read the file. Distributed I/O is not supported
! Read into temporary variable and scatter
      IF (myPE .eq. PE_IO) THEN

! Allocate and initialize temporary variables.
         ALLOCATE (part_data(particles,n_vars)); part_data=0.0D0
         ALLOCATE (lpar_proc(particles))

! Detecting number of columns from the first line in the data
         WRITE(ERR_MSG,'(A,/,A,I2)')'Detecting number of columns from the first line in the data.', &
                                    'Expected number of columns = ',n_vars
         CALL LOG_STATUS()

         read (lunit,'(A)',IOSTAT=IOS) first_line

         call check_number_of_columns(first_line, n_vars,&
              XYZ_start, Phase_ID_start,Diameter_start,Density_start,Velocity_start, &
              Temperature_start,Species_start,User_scalar_start)

! Echo filtering options
         call echo_all_filter_options( &
         PART_IN_X_MIN,PART_IN_X_MAX,PART_IN_X_EXCLUDE, &
         PART_IN_Y_MIN,PART_IN_Y_MAX,PART_IN_Y_EXCLUDE, &
         PART_IN_Z_MIN,PART_IN_Z_MAX, PART_IN_Z_EXCLUDE, &
         PART_IN_PHASE, &
         PART_IN_DIAMETER_MIN, PART_IN_DIAMETER_MAX,PART_IN_DIAMETER_EXCLUDE, &
         PART_IN_DENSITY_MIN, PART_IN_DENSITY_MAX,PART_IN_DENSITY_EXCLUDE, &
         PART_IN_U_MIN, PART_IN_U_MAX,PART_IN_U_EXCLUDE, &
         PART_IN_V_MIN, PART_IN_V_MAX,PART_IN_V_EXCLUDE, &
         PART_IN_W_MIN, PART_IN_W_MAX,PART_IN_W_EXCLUDE, &
         PART_IN_TEMP_MIN, PART_IN_TEMP_MAX,PART_IN_TEMP_EXCLUDE, &
         PART_IN_X_S_MIN,PART_IN_X_S_MAX,PART_IN_X_S_EXCLUDE, &
         PART_IN_USR_VAR_MIN,PART_IN_USR_VAR_MAX,PART_IN_USR_VAR_EXCLUDE)

! If the correct number of columns was detected, go back one line and read the
! data
         backspace(lunit)
! Loop through the input file.
         DO lcurpar = 1, particles
            read (lunit,*,IOSTAT=IOS) (part_data(lcurpar,k),k=1,n_vars)

! Report read errors.
            IF(IOS > 0) THEN
               WRITE(ERR_MSG,1200)
               CALL LOG_ERROR()
               EXIT
 1200 FORMAT('Error reading particle input ',&
         'file.',/'A common error is 2D input for 3D cases.')

! Report End of file errors.
            ELSEIF(IOS < 0) THEN
               WRITE(ERR_MSG,1201) &
                     trim(iVal(lcurpar)), trim(iVal(Particles))
               CALL LOG_ERROR()
               EXIT
 1201 FORMAT('Error reading particle input ',&
         'file.',/'End of file found for particle ',A,' and ',A,1X,    &
         'entries are expected.')

            ENDIF

         ENDDO ! Loop over particles

      ENDIF ! (myPE == PE_IO)

      CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) CALL LOG_ERROR()

! Scatter particles

! set the packet size for transfer
      lpacketsize = n_vars


! build the send buffer in PE_IO proc
! first pass to get the count of particles
      if (particles .eq. UNDEFINED_I) particles = file_particles
      lproc_parcnt(:) = 0
      if(myPE.eq.pe_io) then
         lpar_proc(:) =-1
         do lcurpar = 1,particles
            do lproc= 0,numpes-1
! test if particles belongs to processor's domain
! Do individual tests in each direction (x,y, and z)
               x_test = (part_data(lcurpar,1).ge.xe(istart1_all(lproc)-1).and. &
                         part_data(lcurpar,1).lt.xe(iend1_all(lproc)))
               y_test = (part_data(lcurpar,2).ge.yn(jstart1_all(lproc)-1).and. &
                         part_data(lcurpar,2).lt.yn(jend1_all(lproc)))
               xyz_test = x_test.and.y_test
               if(do_k) then
                  z_test = (part_data(lcurpar,3).ge.zt(kstart1_all(lproc)-1).and. &
                            part_data(lcurpar,3).lt.zt(kend1_all(lproc)))
                  xyz_test = xyz_test.and.z_test
               endif

               if ( xyz_test ) then
                  lpar_proc(lcurpar) = lproc
                  lproc_parcnt(lproc) = lproc_parcnt(lproc) + 1
                  exit
               end if
            end do ! (lproc= 0,numpes-1)
            if (lpar_proc(lcurpar).eq.-1) then
               WRITE(*,501) lcurpar
               call des_mpi_stop
            endif
         enddo ! (lcurpar = 1,particles)
      endif ! if (mype.eq.pe_io)

      call bcast(lproc_parcnt(0:numpes-1),pe_io)

! second pass: set and allocate scatter related variables
      pip = lproc_parcnt(mype)
      call PARTICLE_GROW(pip)
      max_pip = max(pip,max_pip)
      iscr_recvcnt = pip*lpacketsize
      allocate (dprocbuf(iscr_recvcnt))
      if (mype.eq.pe_io) then
         allocate (drootbuf(particles*lpacketsize))
      else
         allocate (drootbuf(10))
      endif

! in the IO processor build the drootbuffer and idispls required
! for mpi communication
      if(mype.eq.pe_io) then
         idispls(0) = 0
         iscattercnts(0) = lproc_parcnt(0)*lpacketsize
         do lproc = 1,numpes-1
            idispls(lproc) = idispls(lproc-1) + iscattercnts(lproc-1)
            iscattercnts(lproc) = lproc_parcnt(lproc)*lpacketsize
         end do
         lproc_parcnt(:) = 0
         do lcurpar = 1,particles
            lproc = lpar_proc(lcurpar)
            lbuf = idispls(lproc)+lproc_parcnt(lproc)*lpacketsize+1
            drootbuf(lbuf:lbuf+lpacketsize-1) = part_data(lcurpar,:)
            lproc_parcnt(lproc) = lproc_parcnt(lproc) + 1
         enddo
      endif

      call desmpi_scatterv(ptype=2)

! unpack the particles in each processor and set the pip
      lc1 = 0
      do lcurpar = 1,pip
         lbuf = (lcurpar-1)*lpacketsize

         lc1 = lc1 + 1 ! initialize potential new particle ID
         keep_me = .true.  ! used to filter particles based on variable range

! Particle position (x,y,z, coordinates)
         des_pos_new(lc1,1:rdimn) = dprocbuf(lbuf + XYZ_start:lbuf + XYZ_start + rdimn - 1 )

         call filter_single_particle_based_on_min_max_var(des_pos_new(lc1,1),'x-coordinate',part_in_x_min,part_in_x_max,part_in_x_exclude,keep_me)
         call filter_single_particle_based_on_min_max_var(des_pos_new(lc1,2),'y-coordinate',part_in_y_min,part_in_y_max,part_in_y_exclude,keep_me)
         if(rdimn==3)call filter_single_particle_based_on_min_max_var(des_pos_new(lc1,3),'z-coordinate',part_in_z_min,part_in_z_max,part_in_z_exclude,keep_me)

! Particle Phase ID
         if(Phase_ID_start>0) then
            PIJK(lc1,5) = INT(dprocbuf(lbuf + Phase_ID_start))
         else
            PIJK(lc1,5) = INT(Uniform_Phase_ID)
         endif

         IF(.NOT.part_in_phase(PIJK(lc1,5))) keep_me = .FALSE.

! Particle Radius
         if(Diameter_start>0) then
            des_radius(lc1) = 0.5D0*dprocbuf(lbuf + Diameter_start)
         else
            des_radius(lc1) = 0.5D0*Uniform_Diameter
         endif

         call filter_single_particle_based_on_min_max_var(2.0D0*des_radius(lc1),'Diameter',part_in_diameter_min,part_in_diameter_max,part_in_diameter_exclude,keep_me)

! Statistical weight, radius and physical radius (CGDEM only)
! The diameter in particle_input.dat is the coarse-grained (not physical) diameter
         if(CGDEM) then
            m = PIJK(lc1,5)
            des_cgp_stw(lc1) = cgp_stat_wt(m)
            des_cgp_rpr(lc1) = des_radius(lc1)/(des_cgp_stw(lc1)**(1.0d0/3.0d0))
         endif

! Particle Density
         if(Density_start>0) then
            ro_sol(lc1) = dprocbuf(lbuf + Density_Start)
         else
            ro_sol(lc1) = Uniform_Density
         endif

         call filter_single_particle_based_on_min_max_var(ro_sol(lc1),'Density',part_in_density_min,part_in_density_max,part_in_density_exclude,keep_me)

! Particle Velocity
         if(Velocity_start>0) then
            des_vel_new(lc1,1:rdimn) = dprocbuf(lbuf + Velocity_start:lbuf + Velocity_start + rdimn - 1)
         else
            des_vel_new(lc1,1:rdimn) = Uniform_Velocity(1:rdimn)
         endif

         call filter_single_particle_based_on_min_max_var(des_vel_new(lc1,1),'x-velocity',part_in_u_min,part_in_u_max,part_in_u_exclude,keep_me)
         call filter_single_particle_based_on_min_max_var(des_vel_new(lc1,2),'y-velocity',part_in_v_min,part_in_v_max,part_in_v_exclude,keep_me)
         if(rdimn==3)call filter_single_particle_based_on_min_max_var(des_vel_new(lc1,3),'z-velocity',part_in_w_min,part_in_w_max,part_in_w_exclude,keep_me)

! Particle Temperature
         if(ENERGY_EQ) then
            if(Temperature_start>0) then
               des_t_s(lc1) = dprocbuf(lbuf + Temperature_Start)
            else
               des_t_s(lc1) = Uniform_Temperature
            endif

            call filter_single_particle_based_on_min_max_var(des_t_s(lc1),'Temperature',part_in_temp_min,part_in_temp_max,part_in_temp_exclude,keep_me)

         endif

! Particle Species: Always need DIMENSION_N_S values for all particles
         if(ANY_SOLIDS_SPECIES_EQ) then
            if(Species_start>0) then
               des_x_s(lc1,1:DIMENSION_N_S) = dprocbuf(lbuf + Species_start:lbuf + Species_start + DIMENSION_N_S - 1)
            else
               des_x_s(lc1,1:DIMENSION_N_S) = Uniform_Species(1:DIMENSION_N_S)
            endif

            do n = 1,DIMENSION_N_S
               write(label,'("Species",I2)') n
               call filter_single_particle_based_on_min_max_var(des_x_s(lc1,n),label,part_in_x_s_min(n),part_in_x_s_max(n),part_in_x_s_exclude(n),keep_me)
            enddo

         endif

! Particle Scalar: CAUTION THIS ARRAY IS TRANSPOSED, FIRST INDEX IS USR_VAR ID, SECOND INDEX IS PARTICLE ID
         if(DES_USR_VAR_SIZE>0) then
            if(User_scalar_start>0) then
               des_usr_var(1:DES_USR_VAR_SIZE,lc1) = dprocbuf(lbuf + User_scalar_start:lbuf + User_scalar_start + DES_USR_VAR_SIZE - 1)
            else
               des_usr_var(1:DES_USR_VAR_SIZE,lc1) = Uniform_user_scalar(1:DES_USR_VAR_SIZE)
            endif

            do n = 1,DES_USR_VAR_SIZE
                write(label,'("User scalar # ",I2)') n
                call filter_single_particle_based_on_min_max_var(des_usr_var(n,lc1),label,part_in_usr_var_min(n),part_in_usr_var_max(n),part_in_usr_var_exclude(n),keep_me)
            enddo

         endif

! Keep or reject particle based on filtering
         if(keep_me) then
! Set particle status to "normal"
            m = PIJK(lc1,5)
            if(des_rigid_motion(m)) then
               call set_rigid_motion(lc1)
            else
               call set_normal(lc1)
            endif
         else
! Reject current particle and go back in the list of particles
            lc1 = lc1 - 1
         endif

      enddo! end of do lcurpar = 1,pip
      deallocate (dprocbuf,drootbuf)

! Set particle count after filtering
      pip = lc1
      max_pip = pip

      CALL  GLOBAL_ALL_SUM(PIP,FPIP)

      WRITE(ERR_MSG, 490)PARTICLES,FPIP
      CALL LOG_INFO()

! Abort if no particles are left after filtering the data
      IF(FPIP==0) THEN
         WRITE(ERR_MSG, 495)
         CALL LOG_ERROR()
      ENDIF


 490  FORMAT('Number of particles read from particle_input.dat  = ', I9 ,&
            /'Number of particles left after filtering the data = ', I9)

 495  FORMAT('No particles left in the system after filtering the data.', &
            /'Please verify the filter settings.')

 500  FORMAT(/2X,'From: DES_SCATTER_PARTICLE: (0)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')
 501  FORMAT(/2X,'From: DES_SCATTER_PARTICLE: (1)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')

      IF(myPE == PE_IO) deallocate (part_data,lpar_proc)


      IF(myPE == PE_IO) CLOSE(lUNIT)
      RETURN
   END SUBROUTINE READ_PART_INPUT_V2P0

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: READ_PART_INPUT_V3P0                                      !
!                                                                      !
! Purpose: Read the particle input and broadcasts the particle data to !
! respective processors.                                               !
! Version 3.0                                                          !
! Added support for SuperDEM                                           !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_PART_INPUT_V3P0

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      use cdist
      use compar
      use des_allocate, only: particle_grow
      use desmpi
      use desmpi_wrapper, only: DES_MPI_STOP
      use des_thermo, only: des_t_s
      use discretelement, only: des_usr_var_size
      use functions
      use funits
      use geometry, only: NO_K
      use mpi_comm_des, only: desmpi_scatterv
      use mpi_init_des, only: des_scatter_particle
      use mpi_utility
      use parallel_mpi
      use param, only: dimension_n_s
      use param1, only: zero, UNDEFINED_I
      use run, only: any_solids_species_eq, energy_eq
      use des_rxns, only: des_x_s
      use sq_rotation_mod

      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: i,k,ilines,m,n
! index of particle
      INTEGER :: lcurpar, LC1
! local unit
      INTEGER, PARAMETER :: lunit=10
! local filename
      character(255) lfilename
! IO Status:
      INTEGER :: IOS
! Flag to indicate if file exists.
      LOGICAL :: lEXISTS
! Read dimension: 2D vs 3D data
      integer :: RDIMN
! Buffer
      character(32) Buffer
! Option to read or use default value (T/F)
      character(1) option
! Integer buffer
      integer :: int_buffer
! Double precision buffer, can be a list of up to 100 numbers
      integer, dimension(100) :: dp_buffer
! Length of list of variables to read (number of columns)
      integer :: n_vars, n_cols
      integer :: XYZ_start
      double precision, dimension(3) :: Dummy_xyz
      integer :: Phase_ID_start
      double precision :: Uniform_Phase_ID
      integer :: Diameter_start
      double precision :: Uniform_Diameter
      integer :: Density_start
      double precision :: Uniform_Density
      integer :: Velocity_start
      double precision, dimension(3) :: Uniform_Velocity
      integer :: Temperature_start
      double precision :: Uniform_Temperature
      integer :: Species_start
      double precision, dimension(100) :: Uniform_Species
      integer :: User_scalar_start
      double precision, dimension(100) :: Uniform_User_scalar
      integer :: SQP_semi_axis_start
      double precision, dimension(3) :: Uniform_SQP_semi_axis
      integer :: SQP_roundness_start
      double precision, dimension(2) :: Uniform_SQP_roundness
      integer :: SQP_quaternion_start
      double precision, dimension(4) :: Uniform_SQP_quaternion

      double precision, allocatable, dimension(:,:) :: part_data
      integer :: lproc,lbuf,lpacketsize
      integer :: lproc_parcnt(0:numpes-1)
      integer, allocatable, dimension(:) :: lpar_proc
      logical :: x_test, y_test, z_test, xyz_test
      double precision :: q_norm

! Particle filter
      LOGICAL, allocatable, dimension(:) :: keep_particle
      LOGICAL :: keep_me
      character(32) :: label

      character(1024) :: first_line
      Integer :: CC1,CC2

     character(7) :: filter_type
     integer :: file_particles
     Integer :: fpip ! Number of particles after filtering the data

!-----------------------------------------------

      n_vars = 0

      IOS = 0
      RDIMN = merge(2,3,NO_K)

      IF (myPE .eq. PE_IO) THEN
! Read Header
! Skip Instructions
         ilines = 19
         do i = 1, ilines
            read(lunit,*) buffer
         enddo

! Dimension
         read(lunit,*) buffer,int_buffer
         if(RDIMN/=int_buffer) then
            WRITE(ERR_MSG,100) int_buffer, RDIMN
            CALL LOG_ERROR()
         endif

! Number of particles
         read(lunit,*) buffer, file_particles
         if (particles .eq. UNDEFINED_I) particles = file_particles
! Skip header
         do i = 1, 3
            read(lunit,*) buffer
         enddo

! Coordinates option
         call parse_vector_option(lunit, 'Coordinates', .TRUE., rdimn, Dummy_xyz(1:rdimn), XYZ_start, n_cols, n_vars)

! Phase ID option
         call parse_scalar_option(lunit, 'Phase_ID', .TRUE., Uniform_Phase_ID, Phase_ID_start, n_cols, n_vars)

! Diameter option
         call parse_scalar_option(lunit, 'Diameter', .TRUE., Uniform_Diameter, Diameter_start, n_cols, n_vars)

! Density option
         call parse_scalar_option(lunit, 'Density', .TRUE., Uniform_Density, Density_start, n_cols, n_vars)

! Velocity option
          call parse_vector_option(lunit, 'Velocity', .TRUE., rdimn, Uniform_Velocity(1:rdimn), Velocity_start, n_cols, n_vars)

! Temperature option
         call parse_scalar_option(lunit, 'Temperature', ENERGY_EQ, Uniform_Temperature, Temperature_start, n_cols, n_vars)

! Species option
          call parse_vector_option(lunit, 'Species', ANY_SOLIDS_SPECIES_EQ, DIMENSION_N_S, Uniform_Species(1:DIMENSION_N_S), Species_start, n_cols, n_vars)

! User Scalar option
          call parse_vector_option(lunit, 'User scalar', (DES_USR_VAR_SIZE>0), DES_USR_VAR_SIZE, Uniform_User_scalar(1:DES_USR_VAR_SIZE), User_scalar_start, n_cols, n_vars)

! SuperDEM semiaxis
          call parse_vector_option(lunit, 'SuperDEM semiaxis', SuperDEM, 3, Uniform_SQP_semi_axis(1:3), SQP_semi_axis_start, n_cols, n_vars)

! SuperDEM roundness
          call parse_vector_option(lunit, 'SuperDEM roundness', SuperDEM, 2, Uniform_SQP_roundness(1:2), SQP_roundness_start, n_cols, n_vars)

! SuperDEM quaternion
          call parse_vector_option(lunit, 'SuperDEM quaternion', SuperDEM, 4, Uniform_SQP_quaternion(1:4), SQP_quaternion_start, n_cols, n_vars)

 100 FORMAT('Error reading particle input file.',&
           /'Dimension (2D or 3D) does not match simulation setup file.',&
           /'Dimension in particle_input.dat    = ',I3&
           /'Dimension in simulation setup file = ',I3)

! Skip the next three lines. This should get us to line 39 where the data
! starts
! Skip header
         do i = 1, 3
            read(lunit,*) buffer
         enddo
      ENDIF !  IF (myPE .eq. PE_IO)

       ! IF(myPE == PE_IO) THEN
       !    print*,'n_vars            = ',n_vars
       !    print*,'Phase_ID_start    = ',Phase_ID_start
       !    print*,'Diameter_start    = ',Diameter_start
       !    print*,'Density_start     = ',Density_start
       !    print*,'Velocity_start    = ',Velocity_start
       !    print*,'Temperature_start = ',Temperature_start
       !    print*,'Species_start     = ',Species_start
       !    print*,'User_scalar_start = ',User_scalar_start
       !    print*,'SQP_semi_axis     = ',SQP_semi_axis_start
       !    print*,'SQP_roundness     = ',SQP_roundness_start
       !    print*,'SQP_quaternion    = ',SQP_quaternion_start
       ! ENDIF

! Broadcast column indices and default values to all PEs
      call bcast(n_vars,pe_io)

      call bcast(XYZ_start,pe_io)
      call bcast(Phase_ID_start,pe_io)
      call bcast(Uniform_Phase_ID,pe_io)
      call bcast(Diameter_start,pe_io)
      call bcast(Uniform_Diameter,pe_io)
      call bcast(Density_start,pe_io)
      call bcast(Uniform_Density,pe_io)
      call bcast(Velocity_start,pe_io)
      call bcast(Uniform_Velocity,pe_io)
      call bcast(Temperature_start,pe_io)
      call bcast(Uniform_Temperature,pe_io)
      call bcast(Species_start,pe_io)
      call bcast(Uniform_Species,pe_io)
      call bcast(User_scalar_start,pe_io)
      call bcast(Uniform_User_scalar,pe_io)
      call bcast(SQP_semi_axis_start,pe_io)
      call bcast(Uniform_SQP_semi_axis,pe_io)
      call bcast(SQP_roundness_start,pe_io)
      call bcast(Uniform_SQP_roundness,pe_io)
      call bcast(SQP_quaternion_start,pe_io)
      call bcast(Uniform_SQP_quaternion,pe_io)

! Read the file. Distributed I/O is not supported
! Read into temporary variable and scatter
      IF (myPE .eq. PE_IO) THEN
! Allocate and initialize temporary variables.
         ALLOCATE (part_data(particles,n_vars)); part_data=0.0D0
         ALLOCATE (lpar_proc(particles))

! Detecting number of columns from the first line in the data
         WRITE(ERR_MSG,'(A,/,A,I2)')'Detecting number of columns from the first line in the data.', &
                                    'Expected number of columns = ',n_vars
         CALL LOG_STATUS()

         read (lunit,'(A)',IOSTAT=IOS) first_line

         call check_number_of_columns(first_line, n_vars,&
              XYZ_start, Phase_ID_start,Diameter_start,Density_start,Velocity_start, &
              Temperature_start,Species_start,User_scalar_start)

! Echo filtering options
         call echo_all_filter_options( &
         PART_IN_X_MIN,PART_IN_X_MAX,PART_IN_X_EXCLUDE, &
         PART_IN_Y_MIN,PART_IN_Y_MAX,PART_IN_Y_EXCLUDE, &
         PART_IN_Z_MIN,PART_IN_Z_MAX, PART_IN_Z_EXCLUDE, &
         PART_IN_PHASE, &
         PART_IN_DIAMETER_MIN, PART_IN_DIAMETER_MAX,PART_IN_DIAMETER_EXCLUDE, &
         PART_IN_DENSITY_MIN, PART_IN_DENSITY_MAX,PART_IN_DENSITY_EXCLUDE, &
         PART_IN_U_MIN, PART_IN_U_MAX,PART_IN_U_EXCLUDE, &
         PART_IN_V_MIN, PART_IN_V_MAX,PART_IN_V_EXCLUDE, &
         PART_IN_W_MIN, PART_IN_W_MAX,PART_IN_W_EXCLUDE, &
         PART_IN_TEMP_MIN, PART_IN_TEMP_MAX,PART_IN_TEMP_EXCLUDE, &
         PART_IN_X_S_MIN,PART_IN_X_S_MAX,PART_IN_X_S_EXCLUDE, &
         PART_IN_USR_VAR_MIN,PART_IN_USR_VAR_MAX,PART_IN_USR_VAR_EXCLUDE)

! If the correct number of columns was detected, go back one line and read the
! data
         backspace(lunit)
! Loop through the input file.
         DO lcurpar = 1, particles
            read (lunit,*,IOSTAT=IOS) (part_data(lcurpar,k),k=1,n_vars)

! Report read errors.
            IF(IOS > 0) THEN
               WRITE(ERR_MSG,1200)
               CALL LOG_ERROR()
               EXIT
 1200 FORMAT('Error reading particle input ',&
         'file.',/'A common error is 2D input for 3D cases.')

! Report End of file errors.
            ELSEIF(IOS < 0) THEN
               WRITE(ERR_MSG,1201) &
                     trim(iVal(lcurpar)), trim(iVal(Particles))
               CALL LOG_ERROR()
               EXIT
 1201 FORMAT('Error reading particle input ',&
         'file.',/'End of file found for particle ',A,' and ',A,1X,    &
         'entries are expected.')

            ENDIF

         ENDDO ! Loop over particles

      ENDIF ! (myPE == PE_IO)

      CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) CALL LOG_ERROR()

! Scatter particles

! set the packet size for transfer
      lpacketsize = n_vars

! build the send buffer in PE_IO proc
! first pass to get the count of particles

      lproc_parcnt(:) = 0
      if(myPE.eq.pe_io) then
         lpar_proc(:) =-1
         do lcurpar = 1,particles
            do lproc= 0,numpes-1
! test if particles belongs to processor's domain
! Do individual tests in each direction (x,y, and z)
               x_test = (part_data(lcurpar,1).ge.xe(istart1_all(lproc)-1).and. &
                         part_data(lcurpar,1).lt.xe(iend1_all(lproc)))
               y_test = (part_data(lcurpar,2).ge.yn(jstart1_all(lproc)-1).and. &
                         part_data(lcurpar,2).lt.yn(jend1_all(lproc)))
               xyz_test = x_test.and.y_test
               if(do_k) then
                  z_test = (part_data(lcurpar,3).ge.zt(kstart1_all(lproc)-1).and. &
                            part_data(lcurpar,3).lt.zt(kend1_all(lproc)))
                  xyz_test = xyz_test.and.z_test
               endif

               if ( xyz_test ) then
                  lpar_proc(lcurpar) = lproc
                  lproc_parcnt(lproc) = lproc_parcnt(lproc) + 1
                  exit
               end if
            end do ! (lproc= 0,numpes-1)
            if (lpar_proc(lcurpar).eq.-1) then
               WRITE(*,501) lcurpar
               call des_mpi_stop
            endif
         enddo ! (lcurpar = 1,particles)
      endif ! if (my_pe.eq.pe_io)

      call bcast(lproc_parcnt(0:numpes-1),pe_io)

! second pass: set and allocate scatter related variables
      pip = lproc_parcnt(mype)
      call PARTICLE_GROW(pip)
      max_pip = max(pip,max_pip)
      iscr_recvcnt = pip*lpacketsize
      allocate (dprocbuf(iscr_recvcnt))
      if (mype.eq.pe_io) then
         allocate (drootbuf(particles*lpacketsize))
      else
         allocate (drootbuf(10))
      endif

! in the IO processor build the drootbuffer and idispls required
! for mpi communication
      if(mype.eq.pe_io) then
         idispls(0) = 0
         iscattercnts(0) = lproc_parcnt(0)*lpacketsize
         do lproc = 1,numpes-1
            idispls(lproc) = idispls(lproc-1) + iscattercnts(lproc-1)
            iscattercnts(lproc) = lproc_parcnt(lproc)*lpacketsize
         end do
         lproc_parcnt(:) = 0
         do lcurpar = 1,particles
            lproc = lpar_proc(lcurpar)
            lbuf = idispls(lproc)+lproc_parcnt(lproc)*lpacketsize+1
            drootbuf(lbuf:lbuf+lpacketsize-1) = part_data(lcurpar,:)
            lproc_parcnt(lproc) = lproc_parcnt(lproc) + 1
         enddo
      endif

      call desmpi_scatterv(ptype=2)

! unpack the particles in each processor and set the pip
      lc1 = 0
      do lcurpar = 1,pip
         lbuf = (lcurpar-1)*lpacketsize

         lc1 = lc1 + 1 ! initialize potential new particle ID
         keep_me = .true.  ! used to filter particles based on variable range

! Particle position (x,y,z, coordinates)
         des_pos_new(lc1,1:rdimn) = dprocbuf(lbuf + XYZ_start:lbuf + XYZ_start + rdimn - 1 )

         call filter_single_particle_based_on_min_max_var(des_pos_new(lc1,1),'x-coordinate',part_in_x_min,part_in_x_max,part_in_x_exclude,keep_me)
         call filter_single_particle_based_on_min_max_var(des_pos_new(lc1,2),'y-coordinate',part_in_y_min,part_in_y_max,part_in_y_exclude,keep_me)
         if(rdimn==3)call filter_single_particle_based_on_min_max_var(des_pos_new(lc1,3),'z-coordinate',part_in_z_min,part_in_z_max,part_in_z_exclude,keep_me)

! Particle Phase ID
         if(Phase_ID_start>0) then
            PIJK(lc1,5) = INT(dprocbuf(lbuf + Phase_ID_start))
         else
            PIJK(lc1,5) = INT(Uniform_Phase_ID)
         endif

         IF(.NOT.part_in_phase(PIJK(lc1,5))) keep_me = .FALSE.

! Particle Radius
         if(Diameter_start>0) then
            des_radius(lc1) = 0.5D0*dprocbuf(lbuf + Diameter_start)
         else
            des_radius(lc1) = 0.5D0*Uniform_Diameter
         endif

         call filter_single_particle_based_on_min_max_var(2.0D0*des_radius(lc1),'Diameter',part_in_diameter_min,part_in_diameter_max,part_in_diameter_exclude,keep_me)

! Statistical weight, radius and physical radius (CGDEM only)
! The diameter in particle_input.dat is the coarse-grained (not physical) diameter
         if(CGDEM) then
            m = PIJK(lc1,5)
            des_cgp_stw(lc1) = cgp_stat_wt(m)
            des_cgp_rpr(lc1) = des_radius(lc1)/(des_cgp_stw(lc1)**(1.0d0/3.0d0))
         endif

! Particle Density
         if(Density_start>0) then
            ro_sol(lc1) = dprocbuf(lbuf + Density_Start)
         else
            ro_sol(lc1) = Uniform_Density
         endif

         call filter_single_particle_based_on_min_max_var(ro_sol(lc1),'Density',part_in_density_min,part_in_density_max,part_in_density_exclude,keep_me)

! Particle Velocity
         if(Velocity_start>0) then
            des_vel_new(lc1,1:rdimn) = dprocbuf(lbuf + Velocity_start:lbuf + Velocity_start + rdimn - 1)
         else
            des_vel_new(lc1,1:rdimn) = Uniform_Velocity(1:rdimn)
         endif

         call filter_single_particle_based_on_min_max_var(des_vel_new(lc1,1),'x-velocity',part_in_u_min,part_in_u_max,part_in_u_exclude,keep_me)
         call filter_single_particle_based_on_min_max_var(des_vel_new(lc1,2),'y-velocity',part_in_v_min,part_in_v_max,part_in_v_exclude,keep_me)
         if(rdimn==3)call filter_single_particle_based_on_min_max_var(des_vel_new(lc1,3),'z-velocity',part_in_w_min,part_in_w_max,part_in_w_exclude,keep_me)

! Particle Temperature
         if(ENERGY_EQ) then
            if(Temperature_start>0) then
               des_t_s(lc1) = dprocbuf(lbuf + Temperature_Start)
            else
               des_t_s(lc1) = Uniform_Temperature
            endif

            call filter_single_particle_based_on_min_max_var(des_t_s(lc1),'Temperature',part_in_temp_min,part_in_temp_max,part_in_temp_exclude,keep_me)

         endif

! Particle Species: Always need DIMENSION_N_S values for all particles
         if(ANY_SOLIDS_SPECIES_EQ) then
            if(Species_start>0) then
               des_x_s(lc1,1:DIMENSION_N_S) = dprocbuf(lbuf + Species_start:lbuf + Species_start + DIMENSION_N_S - 1)
            else
               des_x_s(lc1,1:DIMENSION_N_S) = Uniform_Species(1:DIMENSION_N_S)
            endif

            do n = 1,DIMENSION_N_S
               write(label,'("Species",I2)') n
               call filter_single_particle_based_on_min_max_var(des_x_s(lc1,n),label,part_in_x_s_min(n),part_in_x_s_max(n),part_in_x_s_exclude(n),keep_me)
            enddo

         endif

! Particle Scalar: CAUTION THIS ARRAY IS TRANSPOSED, FIRST INDEX IS USR_VAR ID, SECOND INDEX IS PARTICLE ID
         if(DES_USR_VAR_SIZE>0) then
            if(User_scalar_start>0) then
               des_usr_var(1:DES_USR_VAR_SIZE,lc1) = dprocbuf(lbuf + User_scalar_start:lbuf + User_scalar_start + DES_USR_VAR_SIZE - 1)
            else
               des_usr_var(1:DES_USR_VAR_SIZE,lc1) = Uniform_user_scalar(1:DES_USR_VAR_SIZE)
            endif

            do n = 1,DES_USR_VAR_SIZE
                write(label,'("User scalar # ",I2)') n
                call filter_single_particle_based_on_min_max_var(des_usr_var(n,lc1),label,part_in_usr_var_min(n),part_in_usr_var_max(n),part_in_usr_var_exclude(n),keep_me)
            enddo

         endif

! Superquadic parameters
         if(SuperDEM) then
            if(SQP_semi_axis_start>0) then
               super_r(lc1,1:3) = dprocbuf(lbuf + SQP_semi_axis_start:lbuf + SQP_semi_axis_start + 3 - 1)
            else
               super_r(lc1,1:3) = Uniform_SQP_semi_axis(1:3)
            endif

            if(SQP_roundness_start>0) then
               super_mn(lc1,1:2) = dprocbuf(lbuf + SQP_roundness_start:lbuf + SQP_roundness_start + 2 - 1)
            else
               super_mn(lc1,1:2) = Uniform_SQP_roundness(1:2)
            endif

            if(SQP_quaternion_start>0) then
               super_q(lc1,1:4) = dprocbuf(lbuf + SQP_quaternion_start:lbuf + SQP_quaternion_start + 4 - 1)
            else
               super_q(lc1,1:4) = Uniform_SQP_quaternion(1:4)
            endif

            ! Normalize quaternion
            q_norm = dsqrt(super_q(lc1,1)**2 + super_q(lc1,2)**2 &
                          +super_q(lc1,3)**2 + super_q(lc1,4)**2)

            if(q_norm == zero) then
               WRITE(ERR_MSG,110)
               CALL LOG_ERROR()
            endif

 110 FORMAT('Error reading particle input file.',&
           /'Quaternion has invalid value of zero.')

            super_q(lc1,1:4) = super_q(lc1,1:4) / q_norm

            IF(PARTICLE_ORIENTATION) THEN
            ! Convert quaternion to orientation vector
               CALL QROTATE(super_q(lc1,1:4),ORIENTATION(lc1,:),INIT_ORIENTATION,2)
            ENDIF
         endif

! Keep or reject particle based on filtering
         if(keep_me) then
! Set particle status to "normal"
            call set_normal(lc1)
         else
! Reject current particle and go back in the list of particles
            lc1 = lc1 - 1
         endif

      enddo

      deallocate (dprocbuf,drootbuf)

! Set particle count after filtering
      pip = lc1
      max_pip = pip

      CALL  GLOBAL_ALL_SUM(PIP,FPIP)

      WRITE(ERR_MSG, 490)PARTICLES,FPIP
      CALL LOG_INFO()

! Abort if no particles are left after filtering the data
      IF(FPIP==0) THEN
         WRITE(ERR_MSG, 495)
         CALL LOG_ERROR()
      ENDIF


 490  FORMAT('Number of particles read from particle_input.dat  = ', I9 ,&
            /'Number of particles left after filtering the data = ', I9)

 495  FORMAT('No particles left in the system after filtering the data.', &
            /'Please verify the filter settings.')

 500  FORMAT(/2X,'From: DES_SCATTER_PARTICLE: (0)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')
 501  FORMAT(/2X,'From: DES_SCATTER_PARTICLE: (1)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')

      IF(myPE == PE_IO) deallocate (part_data,lpar_proc)
      IF(myPE == PE_IO) CLOSE(lUNIT)
      RETURN

   END SUBROUTINE READ_PART_INPUT_V3P0

!----------------------------------------------------------------------!
! Subroutine: READ_PART_INPUT_GluedSphereDEM                            !
!                                                                      !
! Author: Renjie Ke                                    Date: Feb, 2024 !
! Purpose:                                                             !
! rewrite based on read_part_input_v2p0 and for gluedsphereDEM only     !
! read the particle_input.dat and seed the glued-sphere into predefined!
! locations, automatically build the connectivity and global id.       !
! after reading, bcast all glued related arrays to each pe             !
!----------------------------------------------------------------------!
   SUBROUTINE READ_PART_INPUT_GluedSphereDEM

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use cdist
      use compar
      use des_allocate, only: particle_grow
      use desmpi
      use desmpi_wrapper, only: DES_MPI_STOP
      use des_thermo, only: des_t_s
      use discretelement
      use functions
      use funits
      use geometry, only: NO_K
      USE GSP_MATH_MOD, only: gsp_quat_to_exyz
      use mpi_comm_des, only: desmpi_scatterv
      use mpi_init_des, only: des_scatter_particle
      use mpi_utility
      use parallel_mpi
      use param, only: dimension_n_s
      use param1, only: zero, UNDEFINED_I
      use run, only: any_solids_species_eq, energy_eq
      use des_rxns, only: des_x_s
      use sq_rotation_mod
      use constant, only: pi
      use gsp_math_mod, only: SET_MARKED_GSP, SET_MARKED_GSP_STATE

      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: i,k,ilines,m,n
! index of particle
      INTEGER :: lcurpar, LC1, L
! local unit
      INTEGER, PARAMETER :: lunit=10
! local filename
      character(255) lfilename
! IO Status:
      INTEGER :: IOS
! Flag to indicate if file exists.
      LOGICAL :: lEXISTS
! Read dimension: 2D vs 3D data
      integer :: RDIMN
! Buffer
      character(32) Buffer
! Option to read or use default value (T/F)
      character(1) option
! Integer buffer
      integer :: int_buffer
! Double precision buffer, can be a list of up to 100 numbers
      integer, dimension(100) :: dp_buffer
! Length of list of variables to read (number of columns)
      integer :: n_vars, n_cols
      integer :: XYZ_start
      double precision, dimension(3) :: Dummy_xyz
      integer :: Phase_ID_start
      double precision :: Uniform_Phase_ID
      integer :: Diameter_start
      double precision :: Uniform_Diameter
      integer :: Density_start
      double precision :: Uniform_Density
      integer :: Velocity_start
      double precision, dimension(3) :: Uniform_Velocity
      integer :: quat_start
      double precision, dimension(4) :: uniform_gsp_quat
      integer :: Temperature_start
      double precision :: Uniform_Temperature
      integer :: Species_start
      double precision, dimension(100) :: Uniform_Species
      integer :: User_scalar_start
      double precision, dimension(100) :: Uniform_User_scalar
      integer :: SQP_semi_axis_start
      double precision, dimension(3) :: Uniform_SQP_semi_axis
      integer :: SQP_roundness_start
      double precision, dimension(2) :: Uniform_SQP_roundness
      integer :: SQP_quaternion_start
      double precision, dimension(4) :: Uniform_SQP_quaternion

      double precision, allocatable, dimension(:,:) :: part_data
      double precision, allocatable, dimension(:,:) :: cpart_data ! component sphere data
      integer :: lproc,lbuf,lpacketsize
      integer :: lproc_parcnt(0:numpes-1)
      integer, allocatable, dimension(:) :: lpar_proc
      logical :: x_test, y_test, z_test, xyz_test
      double precision :: q_norm

! Particle filter
      LOGICAL, allocatable, dimension(:) :: keep_particle
      LOGICAL :: keep_me
      character(32) :: label

      character(1024) :: first_line
      Integer :: CC1,CC2

      character(7) :: filter_type

      Integer :: fpip ! Number of particles after filtering the data
! GluedSphere use
     integer :: cur_phaseid,cur_nsphere,lstart,lend,pip_start,cur_global, cpart_nvars
     integer :: cphase_id_start,cdensity_start,cvel_start,ctemp_start,cspecies_start,cusr_var_start
     double precision :: Bounding_diameter, cur_rho, cur_mass, cur_pmass, cur_vol, cur_pvol, cur_radius, dist2
     double precision :: exone(3), eyone(3), ezone(3), temp_s2cvec(3), q(4)
     integer :: dummy_gid
     integer :: file_particles
     integer :: tmp_size
!-----------------------------------------------

      n_vars = 0
      IOS = 0
      RDIMN = merge(2,3,NO_K)

! glue variables are reset in des_allocate_mod, no need to do it here
      ! glueMass(:) = 0.d0
      ! gluePos(:,:) = 0.0D0
      ! glueVel(:,:) = 0.0D0
      ! glueOmg(:,:) = 0.0D0
      ! glueAngMom(:,:) = 0.0D0
      ! glueInertia(:,:) = 0.0D0
      ! glueForce(:,:) = 0.0D0
      ! glueTorque(:,:)= 0.0D0
      ! glueVolume(:) = 0.0D0

      cpart_nvars = 22
      ! cpart_data at least have 22 columns
      !  X  Y  Z  Relative Diameter  surfaceArea  W  E  S  N  Top  Bot  Wa  Ea  Sa  Na  Ta  Ba
      !  sc2gpc_vec(L,1),sc2gpc_vec(L,2),sc2gpc_vec(L,3), cur_global(this is ordered by gsp_config.dat), gid
      cphase_id_start=0; cdensity_start=0; cvel_start=0; ctemp_start=0; cspecies_start=0; cusr_var_start=0

      IF (myPE .eq. PE_IO) THEN
! Read Header
! Skip Instructions
         ilines = 18
         do i = 1, ilines
            read(lunit,*) buffer
         enddo

! Dimension
         read(lunit,*) buffer,int_buffer
         if(RDIMN/=int_buffer) then
            WRITE(ERR_MSG,100) int_buffer, RDIMN
            CALL LOG_ERROR()
         endif

! Number of particles
         read(lunit,*) buffer, file_particles
         if (particles .eq. UNDEFINED_I) particles = file_particles

! Skip header
         do i = 1, 3
            read(lunit,*) buffer
         enddo

! Coordinates option
         call parse_vector_option(lunit, 'Coordinates', .TRUE., rdimn, Dummy_xyz(1:rdimn), XYZ_start, n_cols, n_vars)

! Phase ID option
         call parse_scalar_option(lunit, 'Phase_ID', .TRUE., Uniform_Phase_ID, Phase_ID_start, n_cols, n_vars)

! Diameter option
         call parse_scalar_option(lunit, 'Diameter', .TRUE., Uniform_Diameter, Diameter_start, n_cols, n_vars)

! Density option
         call parse_scalar_option(lunit, 'Density', .TRUE., Uniform_Density, Density_start, n_cols, n_vars)

! Velocity option
          call parse_vector_option(lunit, 'Velocity', .TRUE., rdimn, Uniform_Velocity(1:rdimn), Velocity_start, n_cols, n_vars)

! Quat option
          call parse_vector_option(lunit, 'Quaternion', .TRUE., rdimn+1, uniform_gsp_quat(1:rdimn+1), quat_start, n_cols, n_vars)

! Temperature option
         call parse_scalar_option(lunit, 'Temperature', ENERGY_EQ, Uniform_Temperature, Temperature_start, n_cols, n_vars)

! Species option
          call parse_vector_option(lunit, 'Species', ANY_SOLIDS_SPECIES_EQ, DIMENSION_N_S, Uniform_Species(1:DIMENSION_N_S), Species_start, n_cols, n_vars)

! User Scalar option
          call parse_vector_option(lunit, 'User scalar', (DES_USR_VAR_SIZE>0), DES_USR_VAR_SIZE, Uniform_User_scalar(1:DES_USR_VAR_SIZE), User_scalar_start, n_cols, n_vars)

 100 FORMAT('Error reading particle input file.',&
           /'Dimension (2D or 3D) does not match simulation setup file.',&
           /'Dimension in particle_input.dat    = ',I3&
           /'Dimension in simulation setup file = ',I3)

! Skip the next three lines. This should get us to line 39 where the data
! starts
! Skip header
         do i = 1, 3
            read(lunit,*) buffer
         enddo

         ! print*,'n_vars            = ',n_vars
         ! print*,'Phase_ID_start    = ',Phase_ID_start
         ! print*,'Diameter_start    = ',Diameter_start
         ! print*,'Density_start     = ',Density_start
         ! print*,'Velocity_start    = ',Velocity_start
         ! print*,'quat_start        = ',quat_start
         ! print*,'Temperature_start = ',Temperature_start
         ! print*,'Species_start     = ',Species_start
         ! print*,'User_scalar_start = ',User_scalar_start

         ! only those in particle_input.dat but also key for component spheres and component spheres do not have it
         ! then we include it in cpart_data

         if(Phase_ID_start > 0) then
            cpart_nvars = cpart_nvars + 1; cphase_id_start = cpart_nvars
         endif
         if(Density_start > 0) then
            cpart_nvars = cpart_nvars + 1; cdensity_start = cpart_nvars
         endif
         if(Velocity_start > 0) then
            cpart_nvars = cpart_nvars + rdimn; cvel_start = cpart_nvars-rdimn+1
         endif
         if(Temperature_start > 0) then
            cpart_nvars = cpart_nvars + 1; ctemp_start = cpart_nvars
         endif
         if(Species_start > 0) then
            cpart_nvars = cpart_nvars + DIMENSION_N_S; cspecies_start = cpart_nvars-DIMENSION_N_S+1
         endif
         if(User_scalar_start > 0) then
            cpart_nvars = cpart_nvars + DES_USR_VAR_SIZE; cusr_var_start = cpart_nvars-DES_USR_VAR_SIZE+1
         endif

         ! print*,'cphase_id_start   = ',cphase_id_start
         ! print*,'cdensity_start    = ',cdensity_start
         ! print*,'cvel_start        = ',cvel_start
         ! print*,'ctemp_start       = ',ctemp_start
         ! print*,'cusr_var_start    = ',cusr_var_start
         ! print*,'cusr_var_start    = ',cusr_var_start

      ENDIF !IF (myPE .eq. PE_IO)

! Broadcast column indices and default values to all PEs
      call bcast(n_vars,pe_io)
      call bcast(cpart_nvars,pe_io)

      call bcast(XYZ_start,pe_io)
      call bcast(Phase_ID_start,pe_io)
      call bcast(Uniform_Phase_ID,pe_io)
      call bcast(Diameter_start,pe_io)
      call bcast(Uniform_Diameter,pe_io)
      call bcast(Density_start,pe_io)
      call bcast(Uniform_Density,pe_io)
      call bcast(Velocity_start,pe_io)
      call bcast(Uniform_Velocity,pe_io)
      call bcast(quat_start,pe_io)
      call bcast(uniform_gsp_quat,pe_io)
      call bcast(Temperature_start,pe_io)
      call bcast(Uniform_Temperature,pe_io)
      call bcast(Species_start,pe_io)
      call bcast(Uniform_Species,pe_io)
      call bcast(User_scalar_start,pe_io)
      call bcast(Uniform_User_scalar,pe_io)

      call bcast(cphase_id_start,pe_io)
      call bcast(cdensity_start,pe_io)
      call bcast(cvel_start,pe_io)
      call bcast(ctemp_start,pe_io)
      call bcast(cspecies_start,pe_io)
      call bcast(cusr_var_start,pe_io)

! particles is provided in mfx file when using manually seeding
! so it is safe to allocate just use 'particles' on each processor
      ! if(allocated(gid_list)) then
      !    deallocate(gid_list)
      !    allocate(gid_list(particles))
      !    gid_list(:) = 0
      ! endif

! Read the file. Distributed I/O is not supported
! Read into temporary variable and scatter
      IF (myPE .eq. PE_IO) THEN
! Allocate and initialize temporary variables.
         ! particles is read from dat file, replaced by NGluedparticles now
         ! n_vars is set before
         ! note: they are actually all for glued-sphere instead of component sphere
         ALLOCATE (part_data(NGluedParticles, n_vars)); part_data=0.0D0
         ALLOCATE (cpart_data(particles, cpart_nvars)); cpart_data=0.0D0
         ALLOCATE (lpar_proc(particles))
         cur_global = 0
         pip_start = 1

! Detecting number of columns from the first line in the data
         WRITE(ERR_MSG,'(A,/,A,I2)')'Detecting number of columns from the first line in the data.', &
                                    'Expected number of columns = ',n_vars
         CALL LOG_STATUS()

         read (lunit,'(A)',IOSTAT=IOS) first_line

         call check_number_of_columns(first_line, n_vars,&
              XYZ_start, Phase_ID_start,Diameter_start,Density_start,Velocity_start, &
              Temperature_start,Species_start,User_scalar_start)

! Echo filtering options
         call echo_all_filter_options( &
         PART_IN_X_MIN,PART_IN_X_MAX,PART_IN_X_EXCLUDE, &
         PART_IN_Y_MIN,PART_IN_Y_MAX,PART_IN_Y_EXCLUDE, &
         PART_IN_Z_MIN,PART_IN_Z_MAX, PART_IN_Z_EXCLUDE, &
         PART_IN_PHASE, &
         PART_IN_DIAMETER_MIN, PART_IN_DIAMETER_MAX,PART_IN_DIAMETER_EXCLUDE, &
         PART_IN_DENSITY_MIN, PART_IN_DENSITY_MAX,PART_IN_DENSITY_EXCLUDE, &
         PART_IN_U_MIN, PART_IN_U_MAX,PART_IN_U_EXCLUDE, &
         PART_IN_V_MIN, PART_IN_V_MAX,PART_IN_V_EXCLUDE, &
         PART_IN_W_MIN, PART_IN_W_MAX,PART_IN_W_EXCLUDE, &
         PART_IN_TEMP_MIN, PART_IN_TEMP_MAX,PART_IN_TEMP_EXCLUDE, &
         PART_IN_X_S_MIN,PART_IN_X_S_MAX,PART_IN_X_S_EXCLUDE, &
         PART_IN_USR_VAR_MIN,PART_IN_USR_VAR_MAX,PART_IN_USR_VAR_EXCLUDE)

! If the correct number of columns was detected, go back one line and read the
! data
         backspace(lunit)
! Loop through the input file.
         DO lcurpar = 1, NGluedParticles
            read (lunit,*,IOSTAT=IOS) (part_data(lcurpar,k),k=1,n_vars)

! Report read errors.
            IF(IOS > 0) THEN
               WRITE(ERR_MSG,1200)
               CALL LOG_ERROR()
               EXIT
 1200 FORMAT('Error reading particle input ',&
         'file.',/'A common error is 2D input for 3D cases.')

! Report End of file errors.
            ELSEIF(IOS < 0) THEN
               WRITE(ERR_MSG,1201) &
                     trim(iVal(lcurpar)), trim(iVal(NGluedParticles))
               CALL LOG_ERROR()
               EXIT
 1201 FORMAT('Error reading particle input ',&
         'file.',/'End of file found for particle ',A,' and ',A,1X,    &
         'entries are expected.')

            ENDIF
            if(Diameter_start > 0) then
               Bounding_diameter = part_data(lcurpar,Diameter_start)
            else
               ! use the default value!
               Bounding_diameter = Uniform_Diameter
            endif
            gluePos(lcurpar,:) = part_data(lcurpar,1:3)
            if(Velocity_start>0) glueVel(lcurpar,:) = part_data(lcurpar,Velocity_start:Velocity_start+2)
            ! glueQuat could also be read directly
            exone=(/1.0D0, 0.0D0, 0.0D0/)
            eyone=(/0.0D0, 1.0D0, 0.0D0/)
            ezone=(/0.0D0, 0.0D0, 1.0D0/)
            if(quat_start>0) then
               glueQuat(lcurpar,:) = part_data(lcurpar,quat_start:quat_start+rdimn)
            else
               ! use the default value!
               glueQuat(lcurpar,:) = uniform_gsp_quat(:)
            endif
            q(:) = glueQuat(lcurpar,:)
            CALL gsp_quat_to_exyz(q(:),exone,eyone,ezone)
! random q renjieke 3-1-2024 add later

            ! calculate glueMass from sphereRelDia
            if(Phase_ID_start>0) then
               cur_phaseid = part_data(lcurpar,Phase_ID_start)
            else
               cur_phaseid = Uniform_Phase_ID
            endif

            cur_nsphere = sphereNumber(cur_phaseid)

            if(Density_start>0) then
               cur_rho = part_data(lcurpar,Density_start)
            else
               cur_rho = Uniform_Density
            endif

            ! lstart and lend refers to the first and last component sphere belong to this glued-sphere
            ! they are the location in the gsp_config
            lstart = start_index_sphere(cur_phaseid)
            lend = start_index_sphere(cur_phaseid)+cur_nsphere-1
            cur_mass = 0.0D00
            cur_vol = 0.0D0

            DO L = lstart, lend
               cur_global = cur_global + 1
               cur_radius = 0.5*sphereRelDia(L)*Bounding_diameter
               cur_pvol = (4.0D0/3.0D0)*PI*(cur_radius)**3
               cur_pmass = cur_rho * cur_pvol
               cur_vol = cur_vol + cur_pvol
               cur_mass = cur_mass + cur_pmass

               dist2 = sphereRelPos(L,2)*sphereRelPos(L,2) + sphereRelPos(L,3)*sphereRelPos(L,3)
               dist2 = dist2 * Bounding_diameter * Bounding_diameter
               glueInertia(lcurpar,1) = glueInertia(lcurpar,1) + cur_pmass * dist2 + (2.0d0/5.0d0)*cur_pmass*cur_radius**2

               dist2 = sphereRelPos(L,1)*sphereRelPos(L,1) + sphereRelPos(L,3)*sphereRelPos(L,3)
               dist2 = dist2 * Bounding_diameter * Bounding_diameter
               glueInertia(lcurpar,2) = glueInertia(lcurpar,2) + cur_pmass * dist2 + (2.0d0/5.0d0)*cur_pmass*cur_radius**2

               dist2 = sphereRelPos(L,1)*sphereRelPos(L,1) + sphereRelPos(L,2)*sphereRelPos(L,2)
               dist2 = dist2 * Bounding_diameter * Bounding_diameter
               glueInertia(lcurpar,3) = glueInertia(lcurpar,3) + cur_pmass * dist2 + (2.0d0/5.0d0)*cur_pmass*cur_radius**2

               temp_s2cvec(:) = Bounding_diameter * sphereRelPos(L,:)
               ! adding rotation to determine the exact location of component spheres
               cpart_data(cur_global,1:3) = gluePos(lcurpar,:) + temp_s2cvec(1) * exone &
                                                            + temp_s2cvec(2) * eyone &
                                                            + temp_s2cvec(3) * ezone
               cpart_data(cur_global,4) = Bounding_diameter * sphereRelDia(L)
               cpart_data(cur_global,5) = tmp_gp_sa(L) * Bounding_diameter * Bounding_diameter
               cpart_data(cur_global,6:11) = tmp_gp_neighsid(L,:)
               where (cpart_data(cur_global,6:11) /= -1) cpart_data(cur_global,6:11) = cpart_data(cur_global,6:11) + pip_start
               cpart_data(cur_global,12:17) = tmp_gp_neighsa(L,:) * Bounding_diameter
               cpart_data(cur_global,18:20) = Bounding_diameter * sphereRelPos(L,:)
               cpart_data(cur_global,21) = cur_global
               ! gid_list is set here, in cpart_data, it is still global but will localize in each pe after mpi_scatter
               cpart_data(cur_global,22) = lcurpar
               if(cphase_id_start>0) cpart_data(cur_global,cphase_id_start) = cur_phaseid
               if(cdensity_start>0) cpart_data(cur_global,cdensity_start) = cur_rho
               if(cvel_start>0) cpart_data(cur_global,cvel_start:cvel_start+rdimn-1) = glueVel(lcurpar,:)
               if(ctemp_start>0) cpart_data(cur_global,ctemp_start) = part_data(lcurpar,Temperature_start)
               if(cspecies_start>0) cpart_data(cur_global,cspecies_start:cspecies_start+DIMENSION_N_S-1) &
                                        = part_data(lcurpar,Species_start:Species_start+DIMENSION_N_S-1)
               if(cusr_var_start>0) cpart_data(cur_global,cusr_var_start:cusr_var_start+DES_USR_VAR_SIZE-1) &
                                           = part_data(lcurpar,User_scalar_start:User_scalar_start+DES_USR_VAR_SIZE-1)
            ENDDO ! L = lstart, lend
            pip_start = pip_start + cur_nsphere
            glueMass(lcurpar) = cur_mass
	         glueVolume(lcurpar) = cur_vol
            glueDiameter(lcurpar) = 2.0D0 * (cur_vol*(3.0D0/4.0D0)/PI)**(1.0D0/3.0D0)
            glueBounding(lcurpar) = Bounding_diameter
            ! OOPS, preferring all data in particle_input.dat is valid, so 1 can be used to initialize the particle_state_gsp
            particle_state_gsp(lcurpar) = 1
         ENDDO ! Loop over NGluedParticles
      ENDIF ! (myPE == PE_IO)

      CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) CALL LOG_ERROR()

! broadcast all those back to all PE
      call bcast(gluePos,pe_io)
      call bcast(glueVel,pe_io)
      call bcast(glueMass,pe_io)
      call bcast(glueVolume, pe_io)
      call bcast(glueInertia,pe_io)
      call bcast(glueQuat,pe_io)
      call bcast(glueDiameter,pe_io)
      call bcast(glueBounding,pe_io)
      call bcast(particle_state_gsp,pe_io)

      glueForce(:,:) = 0.0D0
      glueTorque(:,:) = 0.0D0
      glueOmg(:,:) = 0.0D0
      glueAngMom(:,:) =0.0D0
! Scatter particles
! set the packet size for transfer
      lpacketsize = cpart_nvars ! at least 22
! build the send buffer in PE_IO proc
! first pass to get the count of particles
      lproc_parcnt(:) = 0
      if(myPE.eq.pe_io) then
         lpar_proc(:) =-1
         do lcurpar = 1,particles ! components loop, using cpart_data instead
            do lproc= 0,numpes-1
! test if particles belongs to processor's domain
! Do individual tests in each direction (x,y, and z)
               x_test = (cpart_data(lcurpar,1).ge.xe(istart1_all(lproc)-1).and. &
                         cpart_data(lcurpar,1).lt.xe(iend1_all(lproc)))
               y_test = (cpart_data(lcurpar,2).ge.yn(jstart1_all(lproc)-1).and. &
                         cpart_data(lcurpar,2).lt.yn(jend1_all(lproc)))
               xyz_test = x_test.and.y_test
               if(do_k) then
                  z_test = (cpart_data(lcurpar,3).ge.zt(kstart1_all(lproc)-1).and. &
                            cpart_data(lcurpar,3).lt.zt(kend1_all(lproc)))
                  xyz_test = xyz_test.and.z_test
               endif

               if ( xyz_test ) then
                  lpar_proc(lcurpar) = lproc
                  lproc_parcnt(lproc) = lproc_parcnt(lproc) + 1
                  exit
               end if
            end do ! (lproc= 0,numpes-1)
            if (lpar_proc(lcurpar).eq.-1) then
               WRITE(*,501) lcurpar
               call des_mpi_stop
            endif
         enddo ! (lcurpar = 1,particles)
      endif ! if (my_pe.eq.pe_io)

      call bcast(lproc_parcnt(0:numpes-1),pe_io)

! second pass: set and allocate scatter related variables
      pip = lproc_parcnt(mype)
      ! MPI initialize pip by particles/numpes
      ! might not large enough and call particle grow
      call PARTICLE_GROW(pip)
      max_pip = max(pip,max_pip)
      iscr_recvcnt = pip*lpacketsize
      allocate (dprocbuf(iscr_recvcnt))
      if (mype.eq.pe_io) then
         allocate (drootbuf(particles*lpacketsize))
      else
         allocate (drootbuf(10))
      endif

! in the IO processor build the drootbuffer and idispls required
! for mpi communication
      if(mype.eq.pe_io) then
         idispls(0) = 0
         iscattercnts(0) = lproc_parcnt(0)*lpacketsize
         do lproc = 1,numpes-1
            idispls(lproc) = idispls(lproc-1) + iscattercnts(lproc-1)
            iscattercnts(lproc) = lproc_parcnt(lproc)*lpacketsize
         end do
         lproc_parcnt(:) = 0
         do lcurpar = 1,particles
            lproc = lpar_proc(lcurpar)
            lbuf = idispls(lproc)+lproc_parcnt(lproc)*lpacketsize+1
            drootbuf(lbuf:lbuf+lpacketsize-1) = cpart_data(lcurpar,:)
            lproc_parcnt(lproc) = lproc_parcnt(lproc) + 1
         enddo
      endif

      call desmpi_scatterv(ptype=2) ! scatter to each pe now!
      gp_sa(:) = 0.0D0
      sc2gpc_vec(:,:) = 0.0D0
      cglobal_id(:) = 0
      if(ENERGY_EQ) then
         gp_neighsid(:,:) = 0
         gp_neighsa(:,:) = 0.0D0
      endif

! unpack the particles in each processor and set the pip
!  1 X 2 Y 3 Z 4 Relative Diameter 5 surfaceArea 6 W 7 E 8 S 9 N 10 Top 11 Bot 12 Wa 13 Ea 14 Sa 15 Na 16 Ta 17 Ba
!  18 sc2gpc_vec(L,1) 19 sc2gpc_vec(L,2) 20 sc2gpc_vec(L,3) 21 cglobal_id(L) 22 gid_list(L)
      lc1 = 0
      do lcurpar = 1,pip
         lbuf = (lcurpar-1)*lpacketsize

         lc1 = lc1 + 1 ! initialize potential new particle ID
         keep_me = .true.  ! used to filter particles based on variable range

! Particle position (x,y,z, coordinates)
         ! XYZ_start is a constant and it will not change
         des_pos_new(lc1,1:rdimn) = dprocbuf(lbuf + 1:lbuf + 1 + rdimn - 1)

         call filter_single_particle_based_on_min_max_var(des_pos_new(lc1,1),'x-coordinate',part_in_x_min,part_in_x_max,part_in_x_exclude,keep_me)
         call filter_single_particle_based_on_min_max_var(des_pos_new(lc1,2),'y-coordinate',part_in_y_min,part_in_y_max,part_in_y_exclude,keep_me)
         if(rdimn==3)call filter_single_particle_based_on_min_max_var(des_pos_new(lc1,3),'z-coordinate',part_in_z_min,part_in_z_max,part_in_z_exclude,keep_me)
         ! unpack sc2gpc_vec, gp_sa, cglobal_id, gid_list
         gp_sa(lc1) = dprocbuf(lbuf+5)
         sc2gpc_vec(lc1,:) = dprocbuf(lbuf+18:lbuf+20)
         cglobal_id(lc1) = dprocbuf(lbuf+21)
         gid_list(lc1) = dprocbuf(lbuf+22) ! index is local pip, value is global gid

! Particle Phase ID
         if(Phase_ID_start>0) then
            PIJK(lc1,5) = INT(dprocbuf(lbuf + cphase_id_start))
         else
            PIJK(lc1,5) = INT(Uniform_Phase_ID)
         endif

         IF(.NOT.part_in_phase(PIJK(lc1,5))) keep_me = .FALSE.

! Particle Radius
         ! not matter bounding using default value or not
         des_radius(lc1) = 0.5D0*dprocbuf(lbuf + 4)

         call filter_single_particle_based_on_min_max_var(2.0D0*des_radius(lc1),'Diameter',part_in_diameter_min,part_in_diameter_max,part_in_diameter_exclude,keep_me)

! Statistical weight, radius and physical radius (CGDEM only)
! The diameter in particle_input.dat is the coarse-grained (not physical) diameter
         if(CGDEM) then
            m = PIJK(lc1,5)
            des_cgp_stw(lc1) = cgp_stat_wt(m)
            des_cgp_rpr(lc1) = des_radius(lc1)/(des_cgp_stw(lc1)**(1.0d0/3.0d0))
         endif

! Particle Density
         if(Density_start>0) then
            ro_sol(lc1) = dprocbuf(lbuf + cdensity_start)
         else
            ro_sol(lc1) = Uniform_Density
         endif

         call filter_single_particle_based_on_min_max_var(ro_sol(lc1),'Density',part_in_density_min,part_in_density_max,part_in_density_exclude,keep_me)

! Particle Velocity
         if(Velocity_start>0) then
            des_vel_new(lc1,1:rdimn) = dprocbuf(lbuf + cvel_start:lbuf + cvel_start + rdimn - 1)
         else
            des_vel_new(lc1,1:rdimn) = Uniform_Velocity(1:rdimn)
         endif

         call filter_single_particle_based_on_min_max_var(des_vel_new(lc1,1),'x-velocity',part_in_u_min,part_in_u_max,part_in_u_exclude,keep_me)
         call filter_single_particle_based_on_min_max_var(des_vel_new(lc1,2),'y-velocity',part_in_v_min,part_in_v_max,part_in_v_exclude,keep_me)
         if(rdimn==3)call filter_single_particle_based_on_min_max_var(des_vel_new(lc1,3),'z-velocity',part_in_w_min,part_in_w_max,part_in_w_exclude,keep_me)

! Particle Temperature
         if(ENERGY_EQ) then
            if(Temperature_start>0) then
               des_t_s(lc1) = dprocbuf(lbuf + ctemp_start)
            else
               des_t_s(lc1) = Uniform_Temperature
            endif
            ! when energy_eq, unpack gp_neighsid and gp_neighsa as well
            gp_neighsid(lc1,:) = dprocbuf(lbuf+6:lbuf+11)
            gp_neighsa(lc1,:) = dprocbuf(lbuf+12:lbuf+17)
            call filter_single_particle_based_on_min_max_var(des_t_s(lc1),'Temperature',part_in_temp_min,part_in_temp_max,part_in_temp_exclude,keep_me)
         endif

! Particle Species: Always need DIMENSION_N_S values for all particles
         if(ANY_SOLIDS_SPECIES_EQ) then
            if(Species_start>0) then
               des_x_s(lc1,1:DIMENSION_N_S) = dprocbuf(lbuf + cspecies_start:lbuf + cspecies_start + DIMENSION_N_S - 1)
            else
               des_x_s(lc1,1:DIMENSION_N_S) = Uniform_Species(1:DIMENSION_N_S)
            endif

            do n = 1,DIMENSION_N_S
               write(label,'("Species",I2)') n
               call filter_single_particle_based_on_min_max_var(des_x_s(lc1,n),label,part_in_x_s_min(n),part_in_x_s_max(n),part_in_x_s_exclude(n),keep_me)
            enddo

         endif

! Particle Scalar: CAUTION THIS ARRAY IS TRANSPOSED, FIRST INDEX IS USR_VAR ID, SECOND INDEX IS PARTICLE ID
         if(DES_USR_VAR_SIZE>0) then
            if(User_scalar_start>0) then
               des_usr_var(1:DES_USR_VAR_SIZE,lc1) = dprocbuf(lbuf + cusr_var_start:lbuf + cusr_var_start + DES_USR_VAR_SIZE - 1)
            else
               des_usr_var(1:DES_USR_VAR_SIZE,lc1) = Uniform_user_scalar(1:DES_USR_VAR_SIZE)
            endif

            do n = 1,DES_USR_VAR_SIZE
                write(label,'("User scalar # ",I2)') n
                call filter_single_particle_based_on_min_max_var(des_usr_var(n,lc1),label,part_in_usr_var_min(n),part_in_usr_var_max(n),part_in_usr_var_exclude(n),keep_me)
            enddo

         endif

! Keep or reject particle based on filtering
         if(keep_me) then
! Set particle status to "normal"
            call set_normal(lc1)
            !call set_normal_gsp(gid_list(lc1))
         else
            ! since gid_list is now local array, below is unneeded
            ! dummy_gid = gid_list(lc1)
            ! CALL SET_NONEXISTENT_GSP(dummy_gid)
            ! IF(numPEs > 1) THEN
            !    tmp_size = NGluedParticles/numpes
            !    tmp_size = MAX(tmp_size,4)
            !    allocate(marked_gsp(tmp_size));marked_gsp(:) = 0
            !    allocate(marked_state(tmp_size));marked_state(:) = -1
            !    CALL SET_MARKED_GSP(dummy_gid)
            !    CALL SET_MARKED_GSP_STATE(particle_state_gsp(dummy_gid))
            ! ENDIF
            ! CALL SET_NONEXISTENT(lc1)
            lc1 = lc1 - 1
         endif

      enddo

      deallocate (dprocbuf,drootbuf)

! Set particle count after filtering
      pip = lc1
      max_pip = pip

      CALL  GLOBAL_ALL_SUM(PIP,FPIP)

      WRITE(ERR_MSG, 490)PARTICLES,FPIP
      CALL LOG_INFO()

! Abort if no particles are left after filtering the data
      IF(FPIP==0) THEN
         WRITE(ERR_MSG, 495)
         CALL LOG_ERROR()
      ENDIF

 490  FORMAT('Number of particles read from particle_input.dat  = ', I9 ,&
            /'Number of particles left after filtering the data = ', I9)

 495  FORMAT('No particles left in the system after filtering the data.', &
            /'Please verify the filter settings.')

 500  FORMAT(/2X,'From: DES_SCATTER_PARTICLE: (0)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')
 501  FORMAT(/2X,'From: DES_SCATTER_PARTICLE: (1)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')

      IF(myPE == PE_IO) deallocate (part_data, cpart_data, lpar_proc)
      IF(myPE == PE_IO) CLOSE(lUNIT)
      RETURN
   END SUBROUTINE READ_PART_INPUT_GluedSphereDEM

   SUBROUTINE HANDLE_HEADER_LINE(first_line, n_vars, header_array)
! Dump header names in the csv file into header_array
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(IN) :: first_line
      INTEGER, INTENT(OUT) :: n_vars
      CHARACTER(LEN=25), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: header_array

      CHARACTER(9999) :: temp_str
      CHARACTER(len=1) :: delimiter
      INTEGER :: str_length, count
      INTEGER :: i,j, start,end

      temp_str = TRIM(first_line)
      str_length = LEN_TRIM(temp_str)
      delimiter = ','
      count = 1
      DO i = 1, str_length
         IF (temp_str(i:i) == delimiter) THEN
            count = count + 1
         ENDIF
      ENDDO

      n_vars = count
      ALLOCATE(header_array(n_vars))

      start = 1
      j = 1
      DO i = 1, str_length
         IF (temp_str(i:i) == delimiter .OR. i == str_length) THEN
            end = i
            IF (i == str_length) end = end + 1
            header_array(j) = to_lowercase(TRIM(ADJUSTL(temp_str(start:end-1))))
            start = i + 1
            j = j + 1
         ENDIF
      ENDDO
   END SUBROUTINE HANDLE_HEADER_LINE

   SUBROUTINE COUNT_PARTICLES_IN_FILE(is_cnt_gsp, lfilename, expected_lpar_num, lpar_cnt, lgsp_cnt, expected_lcomponent_num, PNMAX_LOC)
! Count the number of particles based on the number of rows in the csv file
      USE discretelement, only: sphereNumber, gluedSphereDEM
      USE geometry, only: NO_K
      IMPLICIT NONE
      CHARACTER(255),INTENT(IN) :: lfilename
      LOGICAL,INTENT(IN) :: is_cnt_gsp
      INTEGER,INTENT(IN) :: expected_lpar_num
      INTEGER,INTENT(OUT) :: lpar_cnt, lgsp_cnt
      INTEGER, INTENT(OUT), optional :: expected_lcomponent_num
      INTEGER, INTENT(INOUT), optional :: PNMAX_LOC

      INTEGER, PARAMETER :: lunit_cnt = 12
      INTEGER :: nvars
      INTEGER :: IOS
      INTEGER :: line_cnt, component_cnt, expect_component_cnt
      DOUBLE PRECISION :: col1
      CHARACTER(len=25), dimension(:), allocatable :: header_array
      CHARACTER(1024) :: first_line
      INTEGER :: position,i
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: data_line
      INTEGER :: ub_sphereNumber

      line_cnt = 0
      component_cnt = 0
      expect_component_cnt = 0
      nvars = 0
      position = -1
      PNMAX_LOC = 0
      IF(gluedSphereDEM) ub_sphereNumber = ubound(sphereNumber,1)
      OPEN(lunit_cnt, file=lfilename, status='old', action='READ')
! skip the first line
      READ(lunit_cnt,'(A)',iostat=ios) first_line
      CALL HANDLE_HEADER_LINE(first_line, nvars, header_array)

      do i = 1, size(header_array)
         if(trim(header_array(i)) == 'phase_id') then
            position = i
            exit
         endif
      enddo

      if(position .ne. -1) allocate(data_line(position))

      DO
         IF(is_cnt_gsp) THEN
            IF(position == -1) THEN
               ! just need count, values don't actually matter here.
               READ(lunit_cnt,*,iostat=ios) col1
               PNMAX_LOC = 1
            ELSE
               READ(lunit_cnt,*,iostat=ios) data_line(1:position)
               PNMAX_LOC = max(PNMAX_LOC, int(data_line(position)))
            ENDIF
         ELSE
            IF(position == -1) THEN
               READ(lunit_cnt,*,iostat=ios) col1
               PNMAX_LOC = 1
            ELSE
               READ(lunit_cnt,*,iostat=ios) data_line(1:position)
               PNMAX_LOC = max(PNMAX_LOC, int(data_line(position)))
            ENDIF
         ENDIF

         IF(ios /= 0) exit ! end of file
         IF(is_cnt_gsp) THEN

            IF(PNMAX_LOC .gt. ub_sphereNumber) THEN
               WRITE(ERR_MSG,12000) trim(ival(PNMAX_LOC))
               CALL LOG_ERROR()
   12000 FORMAT('Phase ', A,' is undefined for GSP model, please check your solid materials.')
            ENDIF

            IF(position == -1) THEN
               component_cnt = component_cnt + sphereNumber(int(1))
            ELSE
               component_cnt = component_cnt + sphereNumber(int(data_line(position)))
            ENDIF
         ENDIF
         line_cnt = line_cnt + 1
         IF((line_cnt == expected_lpar_num) .AND. is_cnt_gsp) THEN
            expect_component_cnt = component_cnt
         ENDIF
      ENDDO

      IF(is_cnt_gsp) THEN
         lgsp_cnt = line_cnt
         lpar_cnt = component_cnt
         ! if not count gsps, expect_component_cnt = 0
         IF(expect_component_cnt == 0) expect_component_cnt = component_cnt
         ! if count gsps, assign expect_component_cnt to output directly
         expected_lcomponent_num = expect_component_cnt
      ELSE
         lgsp_cnt = 0
         lpar_cnt = line_cnt
      ENDIF
      CLOSE(lunit_cnt)
      IF(allocated(data_line)) deallocate(data_line)

   END SUBROUTINE COUNT_PARTICLES_IN_FILE

   FUNCTION to_lowercase(str) RESULT(lower_str)
      CHARACTER(len=*), INTENT(IN) :: str
      CHARACTER(len=LEN(str)) :: lower_str
      INTEGER :: i

      lower_str = str
      DO i = 1, LEN(str)
         IF (ICHAR(str(i:i)) >= ICHAR('A') .AND. ICHAR(str(i:i)) <= ICHAR('Z')) THEN
            lower_str(i:i) = CHAR(ICHAR(str(i:i)) + 32)
         END IF
      END DO
   END FUNCTION to_lowercase

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: READ_PART_INPUT_V4P0                                      !
!                                                                       !
! Revised: Renjie Ke                                   Date: July, 2024 !
!                                                                       !
! Purpose: Read the particle input and broadcasts the particle data to  !
! respective processors.                                                !
! Version 4.0                                                           !
! support for all DEM submodels, 2d or 3d                               !
! the input file is particle_input.csv                                  !
! the first row in csv file is header line                              !
! for GSP, particle_input.csv and gsp_config.dat are needed             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE READ_PART_INPUT_V4P0

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      USE cdist
      USE CFNEWVALUES_MOD, only: UPDATE_GSP_PARTICLE_STATE
      USE compar
      USE constant, only: pi
      USE des_allocate, only: particle_grow, particle_grow_gluedspheredem
      USE desmpi
      USE desmpi_wrapper, only: DES_MPI_STOP
      USE des_thermo, only: des_t_s
      USE discretelement, only: des_usr_var_size
      USE functions
      USE funits
      USE geometry, only: NO_K, XLENGTH, YLENGTH, ZLENGTH
      USE mass_outflow_dem_mod, only: delete_particle
      USE mpi_comm_des, only: desmpi_scatterv
      USE mpi_init_des, only: des_scatter_particle
      USE mpi_utility
      USE parallel_mpi
      USE param, only: dimension_n_s
      USE param1, only: zero, undefined, UNDEFINED_I
      USE run, only: any_solids_species_eq, energy_eq
      USE des_rxns, only: des_x_s
      USE sq_rotation_mod
      USE gsp_math_mod, only: gsp_quat_to_exyz
      USE gsp_math_mod, only: set_marked_gsp, set_marked_gsp_state, findIndex_1i
      USE ic, only: DIMENSION_IC, IC_DEFINED
      USE ic, only: ic_x_w, ic_x_e, ic_y_s, ic_y_n, ic_z_b, ic_z_t
      USE ic, only: ic_override_part_input
      USE physprop, only: d_p0, ro_s0, x_s0, RO_Xs0

      IMPLICIT NONE
   !-----------------------------------------------
   ! Local variables
   !-----------------------------------------------
   ! Indices
      INTEGER :: i, j, k, lcurpar, lstart, lend, m, n, lc1, L
   ! IC indices
      INTEGER :: ICV
   ! Local unit
      INTEGER, PARAMETER :: lunit=10
      INTEGER, PARAMETER :: lunit_gsp_distribution=11
   ! Local filename
      CHARACTER(255) lfilename
   ! IO Status:
      INTEGER :: IOS
   ! Flag to indicate if file exists.
      LOGICAL :: lEXISTS
   ! Read dimension: 2D vs 3D data
      INTEGER :: RDIMN
      DOUBLE PRECISION :: MAX_D
   ! Number of variables (columns in data file)
      INTEGER :: n_vars, cpart_nvars, expect_n_vars, distribution_n_vars
   ! Number of particles (for GluedSphereDEM it will be NGluedParticles)
      INTEGER :: n_par
   ! Header line of csv file
      CHARACTER(len=25), dimension(:), allocatable :: header_array
      CHARACTER(len=25), dimension(:), allocatable :: intra_header_array
      CHARACTER(len=50), dimension(:), allocatable :: not_in_header_array
   ! Must have variables locations in the array
      INTEGER :: XYZ_start, Phase_ID_start, Diameter_start, Density_start, Velocity_start
      LOGICAL :: FOUND_XYZ, FOUND_PID, FOUND_DIA, FOUND_ROL, FOUND_VEL
   ! Extra variables locations in the array
      INTEGER :: temperature_start, species_start, species_end, usr_var_start, usr_var_end
      INTEGER :: ctemperature_start, cspecies_start, cusr_var_start
   ! Extra variables flag
      LOGICAL :: FOUND_TEMPERATURE, FOUND_SPECIES, FOUND_USR_VAR
   ! Distribution to ranks in DMP
      DOUBLE PRECISION, allocatable, dimension(:,:) :: part_data
      INTEGER, allocatable, dimension(:) :: lpar_proc
      INTEGER :: lproc, lbuf, lpacketsize
      INTEGER :: lproc_parcnt(0:numpes-1)
      LOGICAL :: x_test, y_test, z_test, xyz_test
      DOUBLE PRECISION :: xloc, yloc, zloc
   ! SuperDEM specific variables
      INTEGER :: superdem_semiaxis_start, superdem_semiaxis_end
      INTEGER :: superdem_roundness_start, superdem_roundness_end
      INTEGER :: superdem_quat_start, superdem_quat_end
      LOGICAL :: found_quat, found_semiaxis, found_roundness
      DOUBLE PRECISION :: q_norm
      INTEGER :: sqp_pid
   ! GluedSphereDEM specific variables
      DOUBLE PRECISION, allocatable, dimension(:,:) :: cpart_data
      DOUBLE PRECISION, allocatable, dimension(:,:) :: distribution_data
      INTEGER :: cur_phaseid, cur_nsphere, cur_global, cstart, pip_start, gid
      INTEGER :: quat_start, quat_end, gp_neighsa_start, gp_neighsid_start
      DOUBLE PRECISION :: Bounding_diameter, cur_rho, cur_mass, cur_pmass, cur_vol, cur_radius, cur_pvol, dist2
      DOUBLE PRECISION :: exone(3), eyone(3), ezone(3), temp_s2cvec(3), gPos(3)
      INTEGER :: dummy_index, temp_size
      INTEGER :: intra_density_start,intra_temperature_start, intra_species_start, intra_User_scalar_start
   ! CGDEM specific variables
      INTEGER :: cgp_stat_wt_start
      LOGICAL :: found_cgp_stat_wt

   ! Particle filter
      LOGICAL, allocatable, dimension(:) :: keep_particle
      LOGICAL :: keep_me
      CHARACTER(32) :: label
      CHARACTER(9999) :: first_line,intra_first_line
      CHARACTER(7) :: filter_type
      INTEGER :: fpip ! Number of particles after filtering the data
      INTEGER :: N_DELETED

   ! GluedSphereDEM intra distribution override
      LOGICAL :: usr_defined_distributions
      CHARACTER(255) :: lfilename_distribution
      LOGICAL :: lexists_distribution
   ! Count the number of particles/component spheres from the csv file
      INTEGER :: par_cnt
   ! Count the number of GSPs from the csv file and the number of GSPs after filtering
      INTEGER :: gsp_cnt
      INTEGER :: expected_component_cnt
      INTEGER :: no_found_columns
      INTEGER :: no_found_cnt
      INTEGER :: PNMAX

! Distributed I/O is not supported by current subroutine
      par_cnt = 0
      usr_defined_distributions = .True.
      IF(GSP_IMPLICIT) usr_defined_distributions = .False.
      IOS = 0
      RDIMN = MERGE(2,3,NO_K)
      MAX_D = 0.0
      PNMAX = 0

      n_vars = 0
      cpart_nvars = 15
      expect_n_vars = 2 ! X, Y in 2D
      IF(do_k) expect_n_vars = 3 ! X, Y, Z coordinates in 3D

! Count the number of particles provided in particle_input.csv
      lfilename = 'particle_input.csv'
      IF(mype == pe_io) THEN
         par_cnt = 0
         gsp_cnt = 0
         expected_component_cnt = 0
         IF(GSP_EXPLICIT) THEN
             CALL COUNT_PARTICLES_IN_FILE(.TRUE., lfilename, particles, par_cnt, gsp_cnt, expected_component_cnt, PNMAX_LOC=PNMAX)
             if (NGluedParticles .eq. UNDEFINED_I) NGluedParticles = gsp_cnt
             call PARTICLE_GROW_GluedSphereDEM(NGluedParticles)
         ELSE
            ! gsp_implicit also goes into here
             CALL COUNT_PARTICLES_IN_FILE(.FALSE., lfilename, particles, par_cnt, gsp_cnt, PNMAX_LOC = PNMAX)
             if (particles .eq. UNDEFINED_I) particles = par_cnt
         ENDIF

         ! TODO consolidate these messages into a single message
         IF(GSP_EXPLICIT) THEN
            IF(particles .gt. gsp_cnt) THEN
               WRITE(ERR_MSG,1200) trim(ival(particles)),'glued-sphere particles',trim(ival(gsp_cnt))
               CALL LOG_ERROR()
            ELSEIF(particles .lt. gsp_cnt) THEN
               NGluedParticles = particles
               WRITE(ERR_MSG,1201) 'glued-sphere particles',trim(ival(particles)),trim(ival(gsp_cnt))
               CALL LOG_MESSAGE(__FILE__,__LINE__,LOGLEVEL_INFO,footer=.false.)
               WRITE(ERR_MSG,1202) 'the first '//trim(ival(particles)),'glued-sphere particles'
               CALL LOG_MESSAGE(__FILE__,__LINE__,LOGLEVEL_INFO,header=.false.,footer=.false.)
               particles = expected_component_cnt
               WRITE(ERR_MSG,1202) trim(ival(particles)),'component spheres'
               CALL LOG_MESSAGE(__FILE__,__LINE__,LOGLEVEL_INFO,header=.false.)
            ELSE
               particles = expected_component_cnt
            ENDIF
        ELSE
            if (particles .eq. UNDEFINED_I) particles = par_cnt
            IF(particles .gt. par_cnt) THEN
               WRITE(ERR_MSG,1200) trim(ival(particles)),'particles',trim(ival(par_cnt))
               CALL LOG_ERROR()
            ELSEIF(particles .lt. par_cnt) THEN
               WRITE(ERR_MSG,1201) 'particles',trim(ival(particles)),trim(ival(par_cnt))
               CALL LOG_MESSAGE(__FILE__,__LINE__,LOGLEVEL_INFO,header=.false.,footer=.false.)
               IF(GSP_IMPLICIT) THEN
                  WRITE(ERR_MSG,1202) 'the first '//trim(ival(particles)),'glued-sphere particles'
                  CALL LOG_MESSAGE(__FILE__,__LINE__,LOGLEVEL_INFO,header=.false.)
               ELSE
                  WRITE(ERR_MSG,1202) 'the first '//trim(ival(particles)),'particles'
                  CALL LOG_MESSAGE(__FILE__,__LINE__,LOGLEVEL_INFO,header=.false.)
               ENDIF
            ENDIF
         ENDIF
 1200 FORMAT('Expecting ',A,' ',A,' in particle_input.csv, but there are only ',A,'.')
 1201 FORMAT('Number of ',A,' requested is ',A,'.' &
            /'Particle_input.csv contains ',A,'.')
 1202 FORMAT('Only read ',A,' ',A,' from the file.')

      ENDIF

      CALL BCAST(NGLUEDPARTICLES)
      CALL BCAST(PARTICLES)

! Read Header line and save a header array
      IF (myPE .eq. PE_IO) THEN
         OPEN(lunit, file='particle_input.csv', status='old', action='READ')
         READ(lunit, '(A)', IOSTAT=IOS) first_line
         IF (IOS /= 0) CALL LOG_ERROR()
         first_line = trim(first_line)
! Header_array (case insensitive)
         CALL HANDLE_HEADER_LINE(first_line, n_vars, header_array)
      ENDIF !  IF (myPE .eq. PE_IO)

! Broadcast column indices and default values to all PEs
      CALL bcast(n_vars,pe_io)

! Read the file. Distributed I/O is not supported
      IF (myPE .eq. PE_IO) THEN
! Only x,y,z and phase_id are must have variables
         XYZ_start = 0
         Phase_ID_start = 0
         Diameter_start = 0
         Density_start = 0
         Velocity_start = 0

         FOUND_XYZ = .false.
         FOUND_PID = .false.
         FOUND_DIA = .false.
         FOUND_ROL = .false.
         FOUND_VEL = .false.

! Extra variables -- temperature, species and user scalar
         temperature_start = 0
         species_start = 0
         usr_var_start = 0
         species_end = 0
         usr_var_end = 0

         IF(GSP_EXPLICIT) THEN
            ctemperature_start = 0
            cspecies_start = 0
            cusr_var_start = 0
         ENDIF

         FOUND_TEMPERATURE = .false.
         FOUND_SPECIES = .false.
         FOUND_USR_VAR = .false.

! Allocate and initialize temporary variables.
         IF(GSP_EXPLICIT) THEN
            ALLOCATE (part_data(NGluedParticles, n_vars)); part_data=0.0D0
            ALLOCATE (lpar_proc(particles))
            cur_global = 0
            pip_start = 1
            n_par = NGluedParticles
            ! GluedsphereDEM model specific variables
            quat_start = 0
            quat_end = 0
            found_quat = .false.
         ELSE
            ALLOCATE (part_data(particles,n_vars)); part_data=0.0D0
            ALLOCATE (lpar_proc(particles))
            n_par = particles
            IF(GSP_IMPLICIT) THEN
               found_quat = .false.
               quat_start = 0
               quat_end = 0
            ENDIF
            IF(SuperDEM) THEN
               ! SuperDEM model specific variables
               superdem_semiaxis_start = 0
               superdem_roundness_start = 0
               superdem_quat_start = 0
               superdem_semiaxis_end = 0
               superdem_roundness_end = 0
               superdem_quat_end = 0
               found_quat = .false.
               found_semiaxis = .false.
               found_roundness = .false.
            ENDIF
            IF(CGDEM) THEN
               ! CGDEM model specific variable --> cgp statistical weight
               cgp_stat_wt_start = 0
               found_cgp_stat_wt = .false.
            ENDIF
         ENDIF

! Perform a read and check on the header array to determine if a specific header name is existed
         DO i = 1, size(header_array)
            IF(TRIM(header_array(i)) == 'x') THEN
               FOUND_XYZ = .true.
               XYZ_start = i
            ELSE IF(index(TRIM(header_array(i)), 'phase_id') /= 0) THEN
               FOUND_PID = .true.
               Phase_ID_start = i
               expect_n_vars = expect_n_vars + 1
            ELSE IF(index(TRIM(header_array(i)), 'diameter') /= 0) THEN
               FOUND_DIA = .true.
               Diameter_start = i
               expect_n_vars = expect_n_vars + 1
            ELSE IF(index(TRIM(header_array(i)), 'density') /= 0) THEN
               FOUND_ROL = .true.
               Density_start = i
               expect_n_vars = expect_n_vars + 1
            ELSE IF(TRIM(header_array(i)) == 'u') THEN
               FOUND_VEL = .true.
               Velocity_start = i
               expect_n_vars = expect_n_vars + 2 ! u, v
               IF(do_k) expect_n_vars = expect_n_vars + 1 ! w
            ELSE IF(index(TRIM(header_array(i)), 't_s') /= 0 .or. index(TRIM(header_array(i)), 'temperature') /= 0) THEN
               FOUND_TEMPERATURE = .true.
               temperature_start = i
               cpart_nvars = cpart_nvars + 1 ! t_s
               cpart_nvars = cpart_nvars + 12 ! neighsid(6), neighsa(6)
            ELSE IF(index(TRIM(header_array(i)), 'species') /= 0) THEN
               FOUND_SPECIES = .true.
               species_end = max(species_end, i)
               cpart_nvars = cpart_nvars + 1
            ELSE IF(index(TRIM(header_array(i)), 'usr_var') /= 0) THEN
               FOUND_USR_VAR = .true.
               usr_var_end = max(usr_var_end, i)
               cpart_nvars = cpart_nvars + 1
! For CGDEM, look for cgp_stat_wt
            ELSE IF(index(TRIM(header_array(i)), 'cgp_stat_wt') /= 0 .AND. CGDEM) THEN
               found_cgp_stat_wt = .true.
               cgp_stat_wt_start = i
               expect_n_vars = expect_n_vars + 1
! For GluedSphereDEM, look for quaternion
            ELSE IF(index(TRIM(header_array(i)), 'gsp_q1') /= 0 .AND. GluedsphereDEM) THEN
               found_quat = .true.
               quat_start = i
               expect_n_vars = expect_n_vars + 4
! For SuperDEM model, look for quaternion, semiaxis, roundness
            ELSE IF(index(TRIM(header_array(i)), 'sqp_a') /= 0 .AND. SuperDEM) THEN
               found_semiaxis = .true.
               superdem_semiaxis_start = i
               expect_n_vars = expect_n_vars + 3
            ELSE IF(index(TRIM(header_array(i)), 'sqp_m') /= 0 .AND. SuperDEM) THEN
               found_roundness = .true.
               superdem_roundness_start = i
               expect_n_vars = expect_n_vars + 2
            ELSE IF(index(TRIM(header_array(i)), 'sqp_q1') /= 0 .AND. SuperDEM) THEN
               found_quat = .true.
               superdem_quat_start = i
               expect_n_vars = expect_n_vars + 4
            ENDIF
         ENDDO ! end of DO i = 1, size(header_array)

         IF(ENERGY_EQ .and. FOUND_TEMPERATURE) expect_n_vars = expect_n_vars + 1
         IF(any_solids_species_eq .and. FOUND_SPECIES) expect_n_vars = expect_n_vars + DIMENSION_N_S
         IF(DES_USR_VAR_SIZE>0 .and. FOUND_USR_VAR) expect_n_vars = expect_n_vars + des_usr_var_size

         IF(FOUND_SPECIES) species_start = species_end - DIMENSION_N_S + 1
         IF(FOUND_USR_VAR) usr_var_start = usr_var_end - DES_USR_VAR_SIZE + 1

! Check if expected number of columns is equal to n_vars
         IF(expect_n_vars .NE. n_vars) THEN
            WRITE(ERR_MSG,1204) expect_n_vars, n_vars
            CALL LOG_ERROR()
         ENDIF

 1204 FORMAT('The expected number of columns does not match the actual number of columns in particle_input.csv,',&
            /'expected number of columns = ', I4,',' &
            /'number of columns in file = ', I4,'.')

! Calculate the number of variables read from keyword or default values
         ! dimension_n_s at least is one, so no_found_columns will be at least one
         no_found_columns = 9 + 1 + dimension_n_s + des_usr_var_size - expect_n_vars
         IF(SuperDEM) no_found_columns = 9 + 1 + dimension_n_s + des_usr_var_size + 3 + 2 + 4 - expect_n_vars
         IF(GluedSphereDEM) no_found_columns = 9 + 1 + dimension_n_s + des_usr_var_size + 4 - expect_n_vars
         IF(CGDEM) no_found_columns = 9 + 1 + dimension_n_s + des_usr_var_size + 1 - expect_n_vars

         allocate(not_in_header_array(no_found_columns))
         not_in_header_array(:) = ""
         no_found_cnt = 0

! Data check if some columns are not presented
         IF(.NOT. FOUND_XYZ) THEN
            WRITE(ERR_MSG,1205) 'particle coordinates'
            CALL LOG_ERROR()
         ENDIF

         IF(.NOT. FOUND_PID) THEN
            IF(DES_MMAX .gt. 1) THEN
               WRITE(ERR_MSG,1206) 'phase'
               CALL LOG_ERROR()
            ELSE
               no_found_cnt = no_found_cnt + 1
               not_in_header_array(no_found_cnt) = 'phase_id (1)'
            ENDIF
         ENDIF
 1205 FORMAT('Missing ', A, ' info from particle_input.csv')
 1206 FORMAT('Missing ', A, ' info from particle_input.csv when number of discrete solids phases is larger than 1.')

         ! diameter and density can be read from GUI "solids" panel.
         ! If not finding "diameter" column, d_p0 should be defined.
         IF (.NOT. FOUND_DIA) THEN
            IF(ALL(d_p0 == undefined)) THEN
               WRITE(ERR_MSG,1208) 'diameter','d_p0 is'
               CALL LOG_ERROR()
            ELSE
               no_found_cnt = no_found_cnt + 1
               not_in_header_array(no_found_cnt) = 'diameter (d_p0)'
            ENDIF

            IF(FOUND_PID .and. d_p0(PNMAX) == undefined) THEN
               WRITE(ERR_MSG,1207) trim(ival(PNMAX))
               CALL LOG_ERROR()
            ENDIF
 1207 FORMAT('Diameter for phase ', A,' is undefined in solver, please check solid materials.')
         ENDIF

         ! If not found "density" column, either RO_S0 or RO_Xs0 should be defined.
         IF (.NOT. FOUND_ROL) THEN
            IF(ALL(RO_S0 == undefined) .AND. ALL(RO_Xs0 == undefined)) THEN
               WRITE(ERR_MSG,1208) 'density','ro_s0/ro_xs0 is'
               CALL LOG_ERROR()
            ELSE
               no_found_cnt = no_found_cnt + 1
               not_in_header_array(no_found_cnt) = 'density (RO_S0/RO_Xs0)'
            ENDIF
         ENDIF

         ! velocity can be read from IC, but if not found "U, V, W" columns, they can just be [0, 0, 0] or [0, 0] in 2d
         ! they also can be overwritten by ic_u_s/ic_v_s/ic_w_s if ic_overwrite flag is turned on
         ! so they can't be error, because at least they can have zero velocities.
         IF (.NOT. FOUND_VEL) THEN
            no_found_cnt = no_found_cnt + 1
            not_in_header_array(no_found_cnt) = 'velocity (/0.0, 0.0, 0.0/)'
         ENDIF

         ! so temperature can at least have default 273.15
         IF(.NOT. FOUND_TEMPERATURE .AND. ENERGY_EQ) THEN
            cpart_nvars = cpart_nvars + 1
            cpart_nvars = cpart_nvars + 12
            no_found_cnt = no_found_cnt + 1
            not_in_header_array(no_found_cnt) = 't_s/temperature (273.15)'
         ENDIF

         ! If not found "species1,species2..." columns, X_S0 should be defined.
         IF(.NOT. FOUND_SPECIES .AND. ANY_SOLIDS_SPECIES_EQ) THEN
            IF(ALL(X_S0 == undefined)) THEN
            ! case 1: if constant density is used, then x_s0 will be undefined. So 'ALL(X_S0 == undefined)' cannot be an error here.
            ! case 2: if variable density is used, then undefined x_s0 will be an error.
            ! but if it is case 2, simulation will break before enter par input reading
            ! so case 1 will be the only case to here !
               WRITE(ERR_MSG,12081)
               CALL LOG_INFO()
               X_s0(:,:) = 0.0
               X_s0(:,1) = 1.0
            ENDIF
            cpart_nvars = cpart_nvars + DIMENSION_N_S
            no_found_cnt = no_found_cnt + 1
            not_in_header_array(no_found_cnt) = 'species (X_s0)'
         ENDIF

         ! There are no keywords or initial default values for des_usr_var
         ! But currently if the user does not provide usr_var column, set them all to default 1.0 @renjieke 2024-10-16
         IF(.NOT. FOUND_USR_VAR .AND. DES_USR_VAR_SIZE>0) THEN
            cpart_nvars = cpart_nvars + des_usr_var_size
            no_found_cnt = no_found_cnt + 1
            not_in_header_array(no_found_cnt) = 'user_scalar (1.0)'
         ENDIF

! CGP model, if not found statistical weight, using keyword cgp_stat_wt, GUI already forced it to >= 1
         IF(.NOT. found_cgp_stat_wt .and. CGDEM) THEN
            no_found_cnt = no_found_cnt + 1
            not_in_header_array(no_found_cnt) = 'cgp statistical weight (cgp_stat_wt)'
         ENDIF
! In GSP model, quaternion is required.
         IF(.NOT. found_quat .and. GluedSphereDEM) THEN
            no_found_cnt = no_found_cnt + 1
            not_in_header_array(no_found_cnt) = 'gsp quaternion (gsp_q1,q2,q3,q4)'
         ENDIF

! In SQP model, quat and superquadric parameters are required.
         IF(.NOT. found_semiaxis .and. SuperDEM) THEN
            no_found_cnt = no_found_cnt + 1
            not_in_header_array(no_found_cnt) = 'sqp semiaxis (sqp_a,b,c)'
         ENDIF

         IF(.NOT. found_roundness .and. SuperDEM) THEN
            no_found_cnt = no_found_cnt + 1
            not_in_header_array(no_found_cnt) = 'sqp roundness (sqp_m,n)'
         ENDIF

         IF(.NOT. found_quat .and. SuperDEM) THEN
            no_found_cnt = no_found_cnt + 1
            not_in_header_array(no_found_cnt) = 'sqp quaternion (sqp_q1,q2,q3,q4)'
         ENDIF

 1208 FORMAT('Missing ', A, ' info from particle_input.csv and the keyword ',A,' undefined.')
 12081 FORMAT('Warning: Constant density particle is used in the reacting flow,',&
             /'species mass fractions are set to 1.0, 0.0, 0.0, ...,',&
             /'please turn on ic_override_part_input or provide them in csv.')
! Output to let users know which variables are read from file and which variables are read from keyword/default values
! it's awkward to do this with error_manager so we print the headers and footers ourselves
#ifdef PYMFIX
         write(*,*) ">>>>>================================================================="
         write(*,'(A,I4)') "MESSAGE from des/read_particle_input.f:", __LINE__
#else
         print*,"=========================== Particle input ==========================="

         write(*,*)
#endif
         write(*,*) "Variables read from the csv file:"
         write(*,'(A)',advance='no') '    | '
         DO i = 1, size(header_array)
            if(mod(i,6) == 0) then
               write(*,*)
               if(trim(header_array(i)) == '') exit
               write(*,'(A)',advance='no') '    | '
            endif
            write(*,'(A)',advance='no') trim(header_array(i))
            if(trim(header_array(i)) == '') exit
            write(*,'(A)',advance='no') ' | '
         ENDDO
         write(*,*)
         write(*,*)
         write(*,*) "variables read from keywords/set to default:"
         if(all(not_in_header_array=="")) then
            write(*,'(A)') '    None'
         else
            write(*,'(A)',advance='no') '    | '
            DO i = 1, size(not_in_header_array)
               if(mod(i,4) == 0) then
                  write(*,*)
                  if(trim(not_in_header_array(i)) == '') exit
                  write(*,'(A)',advance='no') '    | '
               endif
               write(*,'(A)',advance='no') trim(not_in_header_array(i))
               if(trim(not_in_header_array(i)) == '') exit
               write(*,'(A)',advance='no') ' | '
            ENDDO
            write(*,*)
         endif ! end of all(not_in_header_array=="")
         write(*,*)
         write(*,*) "IC variables such as temperature, velocities and species will be"
         write(*,*) "overridden if keyword IC_OVERRIDE_PART_INPUT is enabled."
#ifdef PYMFIX
         write(*,*) "<<<<<================================================================="
#else
         write(*,*)
         print*,"======================================================================"
#endif

! Model specific variables start and end position in particle_input.csv
         IF(GluedSphereDEM .AND. found_quat) quat_end = quat_start + 3
         IF(GSP_EXPLICIT) THEN
            ALLOCATE (cpart_data(particles, cpart_nvars))
            cpart_data = ZERO
         ENDIF
         IF(SuperDEM) THEN
            IF(found_semiaxis) superdem_semiaxis_end = superdem_semiaxis_start + 2
            IF(found_roundness) superdem_roundness_end = superdem_roundness_start + 1
            IF(found_quat) superdem_quat_end = superdem_quat_start + 3
         ENDIF

! Echo filtering options
         CALL echo_all_filter_options( &
         PART_IN_X_MIN,PART_IN_X_MAX,PART_IN_X_EXCLUDE, &
         PART_IN_Y_MIN,PART_IN_Y_MAX,PART_IN_Y_EXCLUDE, &
         PART_IN_Z_MIN,PART_IN_Z_MAX, PART_IN_Z_EXCLUDE, &
         PART_IN_PHASE, &
         PART_IN_DIAMETER_MIN, PART_IN_DIAMETER_MAX,PART_IN_DIAMETER_EXCLUDE, &
         PART_IN_DENSITY_MIN, PART_IN_DENSITY_MAX,PART_IN_DENSITY_EXCLUDE, &
         PART_IN_U_MIN, PART_IN_U_MAX,PART_IN_U_EXCLUDE, &
         PART_IN_V_MIN, PART_IN_V_MAX,PART_IN_V_EXCLUDE, &
         PART_IN_W_MIN, PART_IN_W_MAX,PART_IN_W_EXCLUDE, &
         PART_IN_TEMP_MIN, PART_IN_TEMP_MAX,PART_IN_TEMP_EXCLUDE, &
         PART_IN_X_S_MIN,PART_IN_X_S_MAX,PART_IN_X_S_EXCLUDE, &
         PART_IN_USR_VAR_MIN,PART_IN_USR_VAR_MAX,PART_IN_USR_VAR_EXCLUDE)

! Override the intra particle distributions
! If ic_overwrite flag is turned on, vel/temperature/composition will be overwritten by ic values later
         IF(usr_defined_distributions .and. gluedsphereDEM) THEN
! Check if intra_distribution.csv exists
            lfilename_distribution = 'intra_distribution.csv'
            INQUIRE(FILE=lfilename_distribution, EXIST=lexists_distribution)

            IF(lexists_distribution) THEN
               gsp_cnt = 0
               par_cnt = 0

               WRITE(ERR_MSG,1211)
               CALL LOG_INFO()

               CALL COUNT_PARTICLES_IN_FILE(.FALSE., lfilename_distribution, particles, par_cnt, gsp_cnt)

! Particle number counted from intra_distribution.csv can be smaller or larger than 'particles' in system
               IF(par_cnt .gt. particles) par_cnt = particles
               ! global id, temperature, number of species, des_usr_var size
               intra_density_start = 0
               intra_temperature_start = 0
               intra_species_start = 0
               intra_User_scalar_start = 0

               OPEN(lunit_gsp_distribution, file=lfilename_distribution, status='old', action='READ')
               ! Process the first line
               READ (lunit_gsp_distribution,'(A)',IOSTAT=IOS) intra_first_line
               intra_first_line = TRIM(intra_first_line)

               CALL HANDLE_HEADER_LINE(intra_first_line, distribution_n_vars, intra_header_array)

               ALLOCATE (distribution_data(par_cnt, distribution_n_vars))
               distribution_data = ZERO

               DO i = 1, size(intra_header_array)
                  IF(index(TRIM(intra_header_array(i)), 'density') /= 0) THEN
                     intra_density_start = i
                  ELSE IF(index(TRIM(intra_header_array(i)), 't_s') /= 0) THEN
                     intra_temperature_start = i
                  ELSE IF(index(TRIM(intra_header_array(i)), 'species1') /= 0) THEN
                     intra_species_start = i
                  ELSE IF(index(TRIM(intra_header_array(i)), 'usr_var1') /= 0) THEN
                     intra_User_scalar_start = i
                  ENDIF
               ENDDO

               DO lcurpar = 1, par_cnt
                  READ (lunit_gsp_distribution,*,IOSTAT=IOS) (distribution_data(lcurpar,k),k=1,distribution_n_vars)
                  IF(IOS > 0) THEN
                     WRITE(ERR_MSG,1209)
                     CALL LOG_ERROR()
                  ELSEIF(IOS < 0) THEN
                     WRITE(ERR_MSG,1210) TRIM(iVal(lcurpar)),TRIM(iVal(par_cnt))
                     CALL LOG_ERROR()
                  ENDIF
               ENDDO
               CLOSE(lunit_gsp_distribution)
            ELSE
               WRITE(ERR_MSG,1212)
               CALL LOG_INFO()
               usr_defined_distributions = .False.
            ENDIF
         ENDIF

 1209 FORMAT('Error reading the intra distribution ',&
         'file.',/'A common error is given too many variables in file.')
 1210 FORMAT('Error reading the intra distribution ',&
         'file.',/'End of file found for particle ',A,' and ',A,1X,    &
         'entries are expected.')

 1211 FORMAT('Using user-defined intra-particle distributions instead of uniform distributions.')
 1212 FORMAT('Could not find intra_distribution.csv. Defaulting to uniform intra-particle distribution instead.')
      ENDIF ! end of PE_IO

      IF(myPE == PE_IO) THEN
   ! Loop through the input file.
         DO lcurpar = 1, n_par
            READ (lunit,*,IOSTAT=IOS) (part_data(lcurpar,k),k=1,n_vars)
   ! Report read errors.
            IF(IOS > 0) THEN
               WRITE(ERR_MSG,1213)
               CALL LOG_ERROR()
               EXIT
   ! Report End of file errors.
            ELSEIF(IOS < 0) THEN
               WRITE(ERR_MSG,1214) &
                     TRIM(iVal(lcurpar)), TRIM(iVal(n_par))
               CALL LOG_ERROR()
               EXIT
            ENDIF

 1213 FORMAT('Error reading particle input ',&
        'file.',/'A common error is 2D input for 3D cases.')
 1214 FORMAT('Error reading particle input ',&
        'file.',/'End of file found for particle ',A,' and ',A,1X,    &
        'entries are expected.')

! Append mother glued-sphere particle info into its component spheres
            IF(GSP_EXPLICIT) THEN
               gluePos(lcurpar,:) = part_data(lcurpar,XYZ_start:XYZ_start+2)

               cur_phaseid = 1
               IF(FOUND_PID) cur_phaseid = part_data(lcurpar,Phase_ID_start)

               cur_nsphere = sphereNumber(cur_phaseid)

               IF(FOUND_DIA) THEN ! this should be here or on the distribution?
                  Bounding_diameter = part_data(lcurpar, Diameter_start)
               ELSE
                  Bounding_diameter = d_p0(cur_phaseid)
               ENDIF

               IF(FOUND_ROL) THEN
                  cur_rho = part_data(lcurpar,Density_Start)
               ELSE
                  cur_rho = ro_s0(cur_phaseid)
               ENDIF

               IF(FOUND_VEL) THEN
                  glueVel(lcurpar,:) = part_data(lcurpar,Velocity_start:Velocity_start+2)
               ELSE
                  glueVel(lcurpar,:) = (/0.0, 0.0, 0.0/)
               ENDIF

               lstart = start_index_sphere(cur_phaseid)
               lend = start_index_sphere(cur_phaseid)+cur_nsphere-1
               cur_mass = 0.0D00
               cur_vol = 0.0D0

               exone=(/1.0D0, 0.0D0, 0.0D0/)
               eyone=(/0.0D0, 1.0D0, 0.0D0/)
               ezone=(/0.0D0, 0.0D0, 1.0D0/)

               IF(found_quat) THEN
                  glueQuat(lcurpar,:) = part_data(lcurpar, quat_start:quat_end)
                  CALL gsp_quat_to_exyz(glueQuat(lcurpar,:),exone,eyone,ezone)
               ELSE
                  glueQuat(lcurpar,:) = (/gsp_q1(cur_phaseid), gsp_q2(cur_phaseid), gsp_q3(cur_phaseid), gsp_q4(cur_phaseid)/)
               ENDIF

               DO L = lstart, lend
                  cur_global = cur_global + 1
                  cur_radius = 0.5*sphereRelDia(L)*Bounding_diameter
                  cur_pvol = (4.0D0/3.0D0)*PI*(cur_radius)**3
                  cur_pmass = cur_rho * cur_pvol
                  cur_vol = cur_vol + cur_pvol
                  cur_mass = cur_mass + cur_pmass

                  DO i = 1, 3
                     j = mod(i, 3) + 1
                     k = mod(i + 1, 3) + 1
                     dist2 = sphereRelPos(L,j) * sphereRelPos(L,j) + sphereRelPos(L,k) * sphereRelPos(L,k)
                     dist2 = dist2 * Bounding_diameter ** (2.0)
                     glueInertia(lcurpar,i) = glueInertia(lcurpar,i) + cur_pmass * dist2 &
                           + (2.0d0 / 5.0d0) * cur_pmass * cur_radius ** 2
                  ENDDO

                  temp_s2cvec(:) = Bounding_diameter * sphereRelPos(L,:)

                  cpart_data(cur_global,1:3) = gluePos(lcurpar,:) + temp_s2cvec(1) * exone &
                                                + temp_s2cvec(2) * eyone &
                                                + temp_s2cvec(3) * ezone

                  cpart_data(cur_global,4) = Bounding_diameter * sphereRelDia(L)
                  cpart_data(cur_global,5) = tmp_gp_sa(L) * Bounding_diameter * Bounding_diameter
                  cpart_data(cur_global,6:8) = temp_s2cvec(:)
                  cpart_data(cur_global,9) = cur_global
                  cpart_data(cur_global,10) = lcurpar
                  cpart_data(cur_global,11) = cur_phaseid
                  cpart_data(cur_global,12) = cur_rho
                  cpart_data(cur_global,13:15) = glueVel(lcurpar,:)

                  cstart = 15

                  IF(energy_eq) THEN
                     cstart = cstart + 1 ! cstart = 16
                     IF(FOUND_TEMPERATURE) THEN
                        cpart_data(cur_global,cstart) = part_data(lcurpar,temperature_start)
                     ELSE
                        cpart_data(cur_global,cstart) = 273.15
                     ENDIF
                     ctemperature_start = cstart
                     cstart = cstart + 1 ! cstart = 17
                     cpart_data(cur_global,cstart:cstart+5) = tmp_gp_neighsid(L,:)
                     where (cpart_data(cur_global,cstart:cstart+5) /= -1) cpart_data(cur_global,cstart:cstart+5) = cpart_data(cur_global,cstart:cstart+5) + pip_start
                     gp_neighsid_start = cstart
                     cstart = cstart + 5 + 1 ! cstart = 23
                     cpart_data(cur_global,cstart:cstart+5) = tmp_gp_neighsa(L,:) * Bounding_diameter * Bounding_diameter
                     gp_neighsa_start = cstart
                     cstart = cstart + 5 ! cstart = 28
                  ENDIF

                  IF(any_solids_species_eq) THEN
                     cstart = cstart + 1
                     IF(FOUND_SPECIES) THEN
                        cpart_data(cur_global,cstart:cstart+DIMENSION_N_S-1) &
                                          = part_data(lcurpar,species_start:species_end)
                     ELSE
                        cpart_data(cur_global,cstart:cstart+DIMENSION_N_S-1) = X_S0(cur_phaseid,1:DIMENSION_N_S)
                     ENDIF
                     cspecies_start = cstart
                     cstart = cstart+DIMENSION_N_S-1
                  ENDIF

                  IF(DES_USR_VAR_SIZE>0) THEN
                     cstart = cstart + 1
                     IF(FOUND_USR_VAR) THEN
                        cpart_data(cur_global,cstart:cstart+DES_USR_VAR_SIZE-1) &
                                             = part_data(lcurpar,usr_var_start:usr_var_end)
                     ELSE
                        cpart_data(cur_global,cstart:cstart+DES_USR_VAR_SIZE-1) &
                                             = 1.0
                     ENDIF
                     cusr_var_start = cstart
                     cstart = cstart+DES_USR_VAR_SIZE-1
                  ENDIF
               ENDDO ! L = lstart, lend
               pip_start = pip_start + cur_nsphere
               glueMass(lcurpar) = cur_mass
               glueVolume(lcurpar) = cur_vol
               glueDiameter(lcurpar) = 2.0D0 * (cur_vol*(3.0D0/4.0D0)/PI)**(1.0D0/3.0D0)
               glueBounding(lcurpar) = Bounding_diameter
               particle_state_gsp(lcurpar) = 1
            ENDIF ! end of IF(GSP_EXPLICIT)
         ENDDO ! lcurpar = 1, n_par

         IF(GSP_EXPLICIT) THEN
            ! Unlike other solid models, GSP distributes cpart_data to different ranks
            ! start position for variables might change, manually set to correct locations
            Diameter_start = 4
            Phase_ID_start = 11
            Density_start = 12
            Velocity_start = 13
            IF(FOUND_TEMPERATURE .AND. energy_eq) temperature_start = ctemperature_start
            IF(FOUND_SPECIES .AND. any_solids_species_eq) species_start = cspecies_start
            IF(FOUND_USR_VAR .AND. DES_USR_VAR_SIZE>0) usr_var_start = cusr_var_start
         ENDIF

         IF(usr_defined_distributions .and. gluedsphereDEM) THEN
            ! replace values in cpart_data by user-defined distribution data
            ! 1:par_cnt only since the row number of intra_distribution.csv can be smaller than 'particles'
            IF(intra_density_start > 0)cpart_data(1:par_cnt,density_start) = distribution_data(1:par_cnt,intra_density_start) ! density distribution might also be non-uniform
            IF(intra_temperature_start > 0) cpart_data(1:par_cnt,ctemperature_start) = distribution_data(1:par_cnt,intra_temperature_start)
            IF(intra_species_start > 0) cpart_data(1:par_cnt,cspecies_start:cspecies_start+DIMENSION_N_S-1) = &
                                                      distribution_data(1:par_cnt,intra_species_start:intra_species_start+DIMENSION_N_S-1)
            IF(intra_User_scalar_start > 0) cpart_data(1:par_cnt,cusr_var_start:cusr_var_start+DES_USR_VAR_SIZE-1) = &
                              distribution_data(1:par_cnt,intra_User_scalar_start:intra_User_scalar_start+DES_USR_VAR_SIZE-1)
         ENDIF
      ENDIF ! (myPE == PE_IO)

      CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) CALL LOG_ERROR()

! Broadcast column indices and default values to all PEs
      CALL BCAST(XYZ_START,PE_IO)
      CALL BCAST(PHASE_ID_START,PE_IO)
      CALL BCAST(DIAMETER_START,PE_IO)
      CALL BCAST(DENSITY_START,PE_IO)
      CALL BCAST(VELOCITY_START,PE_IO)
      CALL BCAST(TEMPERATURE_START,PE_IO)
      CALL BCAST(SPECIES_START,PE_IO)
      CALL BCAST(USR_VAR_START,PE_IO)
      CALL BCAST(FOUND_XYZ,PE_IO)
      CALL BCAST(FOUND_PID,PE_IO)
      CALL BCAST(FOUND_DIA,PE_IO)
      CALL BCAST(FOUND_ROL,PE_IO)
      CALL BCAST(FOUND_VEL,PE_IO)
      CALL BCAST(FOUND_TEMPERATURE,PE_IO)
      CALL BCAST(FOUND_USR_VAR,PE_IO)
      CALL BCAST(FOUND_SPECIES,PE_IO)

! Broadcast all those back to all PE and they are exclusive for GluedSphereDEM model
      IF(GluedsphereDEM) THEN
         IF(GSP_EXPLICIT) THEN
            CALL bcast(gluePos,pe_io)
            CALL bcast(glueVel,pe_io)
            CALL bcast(glueMass,pe_io)
            CALL bcast(glueVolume, pe_io)
            CALL bcast(glueInertia,pe_io)
            CALL bcast(glueQuat,pe_io)
            CALL bcast(glueDiameter,pe_io)
            CALL bcast(glueBounding,pe_io)
            CALL bcast(particle_state_gsp,pe_io)
            CALL bcast(ctemperature_start,pe_io)
            CALL bcast(cspecies_start,pe_io)
            CALL bcast(cusr_var_start,pe_io)
            CALL bcast(cpart_nvars,pe_io)
            CALL bcast(gp_neighsid_start,pe_io)
            CALL bcast(gp_neighsa_start,pe_io)

            glueForce(:,:) = 0.0D0
            glueTorque(:,:) = 0.0D0
            glueOmg(:,:) = 0.0D0
         ELSEIF(GSP_IMPLICIT) THEN
            CALL bcast(quat_start, PE_IO)
         ENDIF
         glueAngMom(:,:) =0.0D0
      ENDIF

      IF(SuperDEM) THEN
         CALL bcast(found_semiaxis, pe_io)
         CALL bcast(found_roundness, pe_io)
         CALL bcast(found_quat, pe_io)
         CALL bcast(superdem_semiaxis_start,pe_io)
         CALL bcast(superdem_roundness_start,pe_io)
         IF(found_quat) CALL bcast(superdem_quat_start,pe_io)
      ENDIF

      IF(CGDEM) THEN
         CALL bcast(found_cgp_stat_wt, pe_io)
         CALL bcast(cgp_stat_wt_start, pe_io)
      ENDIF

! Scatter particles, set the packet size for transfer
      lpacketsize = n_vars
      IF(GSP_EXPLICIT) lpacketsize = cpart_nvars

   ! Build the send buffer in PE_IO proc
   ! First pass to get the count of particles
      lproc_parcnt(:) = 0
      IF(myPE .eq. pe_io) THEN
         lpar_proc(:) = -1
   ! Loop on the particle level, for gsp its component spheres
         DO lcurpar = 1,particles
            DO lproc= 0,numpes-1
   ! Test if particles belongs to processor's domain
   ! Do individual tests in each direction (x,y, and z)
               IF(GSP_EXPLICIT) THEN
                  xloc = cpart_data(lcurpar,1)
                  yloc = cpart_data(lcurpar,2)
                  IF(do_k) zloc = cpart_data(lcurpar,3)
               ELSE
                  xloc = part_data(lcurpar,1)
                  yloc = part_data(lcurpar,2)
                  IF(do_k) zloc = part_data(lcurpar,3)
               ENDIF

               x_test = (xloc .ge. xe(istart1_all(lproc)-1) .AND. &
                        xloc .lt. xe(iend1_all(lproc)))
               y_test = (yloc .ge. yn(jstart1_all(lproc)-1) .AND. &
                        yloc .lt. yn(jend1_all(lproc)))
               xyz_test = x_test.AND.y_test
               IF(do_k) THEN
                  z_test = (zloc .ge. zt(kstart1_all(lproc)-1) .AND. &
                           zloc .lt. zt(kend1_all(lproc)))
                  xyz_test = xyz_test .AND. z_test
               ENDIF

               IF ( xyz_test ) THEN
                  lpar_proc(lcurpar) = lproc
                  lproc_parcnt(lproc) = lproc_parcnt(lproc) + 1
                  exit
               ENDIF
            ENDDO ! (lproc= 0,numpes-1)
            IF (lpar_proc(lcurpar).eq.-1) THEN
               !WRITE(*,501) lcurpar
               !CALL des_mpi_stop
!JFD: Trick: assign out of bounds particles to process 0. These will be
!removed when applying filters (default filter)
               lpar_proc(lcurpar) = 0
               lproc_parcnt(0) = lproc_parcnt(0) + 1
            ENDIF
         ENDDO ! (lcurpar = 1,particles)
      ENDIF ! if (mype .eq. pe_io)

      CALL bcast(lproc_parcnt(0:numpes-1),pe_io)

! Second pass: set and allocate scatter related variables
      pip = lproc_parcnt(mype)
      CALL PARTICLE_GROW(pip)
      max_pip = max(pip,max_pip)
      iscr_recvcnt = pip*lpacketsize
      allocate (dprocbuf(iscr_recvcnt))
      IF (mype.eq.pe_io) THEN
         allocate (drootbuf(particles*lpacketsize))
      ELSE
         allocate (drootbuf(10))
      ENDIF

! In the IO processor build the drootbuffer and idispls required
! For mpi communication
      IF(mype.eq.pe_io) THEN
         idispls(0) = 0
         iscattercnts(0) = lproc_parcnt(0)*lpacketsize
         DO lproc = 1,numpes-1
            idispls(lproc) = idispls(lproc-1) + iscattercnts(lproc-1)
            iscattercnts(lproc) = lproc_parcnt(lproc)*lpacketsize
         ENDDO
         lproc_parcnt(:) = 0
         DO lcurpar = 1,particles
            lproc = lpar_proc(lcurpar)
            lbuf = idispls(lproc)+lproc_parcnt(lproc)*lpacketsize+1
            IF(GSP_EXPLICIT) THEN
               drootbuf(lbuf:lbuf+lpacketsize-1) = cpart_data(lcurpar,:)
            ELSE
               drootbuf(lbuf:lbuf+lpacketsize-1) = part_data(lcurpar,:)
            ENDIF
            lproc_parcnt(lproc) = lproc_parcnt(lproc) + 1
         ENDDO
      ENDIF

      CALL desmpi_scatterv(ptype=2)

! Exclusive for GSP, set the container to zero in case some random values errors
      IF(GSP_EXPLICIT) THEN
         gp_sa(:) = 0.0D0
         sc2gpc_vec(:,:) = 0.0D0
         cglobal_id(:) = 0
         IF(ENERGY_EQ) THEN
            gp_neighsid(:,:) = 0
            gp_neighsa(:,:) = 0.0D0
         ENDIF
      ENDIF

      IF(GSP_EXPLICIT) THEN
         IF(.NOT. ALLOCATED(marked_gsp)) THEN
            temp_size = MAX(NGluedParticles/numpes, 4)
            ALLOCATE(marked_gsp(temp_size)); marked_gsp(:) = 0
            ALLOCATE(marked_state(temp_size)); marked_state(:) = -1
         ENDIF
      ENDIF

! Unpack the particles in each processor and set the pip
      lc1 = 0
      N_DELETED = 0
      DO lcurpar = 1,pip
         lbuf = (lcurpar-1)*lpacketsize

         lc1 = lc1 + 1 ! initialize potential new particle ID
         keep_me = .true.  ! used to filter particles based on variable range

! Particle position (x,y,z, coordinates)
         des_pos_new(lc1,1:rdimn) = dprocbuf(lbuf + 1:lbuf + 1 + rdimn - 1)

         CALL filter_single_particle_based_on_min_max_var(des_pos_new(lc1,1),'x-coordinate',part_in_x_min,part_in_x_max,part_in_x_exclude,keep_me)
         CALL filter_single_particle_based_on_min_max_var(des_pos_new(lc1,2),'y-coordinate',part_in_y_min,part_in_y_max,part_in_y_exclude,keep_me)
         IF(rdimn==3) CALL filter_single_particle_based_on_min_max_var(des_pos_new(lc1,3),'z-coordinate',part_in_z_min,part_in_z_max,part_in_z_exclude,keep_me)
         IF(GSP_EXPLICIT) THEN
            gp_sa(lc1) = dprocbuf(lbuf+5)
            sc2gpc_vec(lc1,:) = dprocbuf(lbuf+6:lbuf+8)
            cglobal_id(lc1) = dprocbuf(lbuf+9)
            gid_list(lc1) = dprocbuf(lbuf+10)
            gPos(:) = gluePos(gid_list(lc1),:)
         ENDIF

! Particle Phase ID
         IF(FOUND_PID) THEN
            PIJK(lc1,5) = INT(dprocbuf(lbuf + Phase_ID_start))
         ELSE
            PIJK(lc1,5) = 1
         ENDIF
         ! also determine the IC after knowing phase id
         ! @renjieke no need to formulate gp_component_ic for gsp_implicit
         IF(GSP_EXPLICIT .AND. ANY(IC_OVERRIDE_PART_INPUT)) THEN
            DO ICV = 1, DIMENSION_IC
               IF (.NOT. IC_DEFINED(ICV)) CYCLE
               IF (.NOT. IC_OVERRIDE_PART_INPUT(ICV,PIJK(lc1,5))) CYCLE
               IF ((ic_x_w(ICV) <= gPos(1) .AND. gPos(1) <= ic_x_e(ICV)) .AND. &
                  (ic_y_s(ICV) <= gPos(2) .AND. gPos(2) <= ic_y_n(ICV)) .AND. &
                  (ic_z_b(ICV) <= gPos(3) .AND. gPos(3) <= ic_z_t(ICV))) THEN
                  gp_component_IC(lc1) = ICV
               ENDIF
            ENDDO
         ENDIF
         IF(.NOT.part_in_phase(PIJK(lc1,5))) keep_me = .FALSE.

! Particle Radius
         ! @renjieke for gsp_implicit, des_radius refers to bounding diameter
         IF(FOUND_DIA .or. (GSP_EXPLICIT)) THEN
            des_radius(lc1) = 0.5D0*dprocbuf(lbuf + Diameter_start)
         ELSE
            des_radius(lc1) = 0.5D0*d_p0(PIJK(LC1,5))
         ENDIF
         ! for CGDEM, this value is used to replace D_P0, not CG_D_P0
         CALL filter_single_particle_based_on_min_max_var(2.0D0*des_radius(lc1),'Diameter',part_in_diameter_min,part_in_diameter_max,part_in_diameter_exclude,keep_me)

! Remove particles that are outside the MFiX box
         CALL filter_single_particle_based_on_min_max_var(des_pos_new(lc1,1),'x-coordinate',x_min,x_max,.false.,keep_me)
         CALL filter_single_particle_based_on_min_max_var(des_pos_new(lc1,2),'y-coordinate',y_min,y_max,.false.,keep_me)
         IF(rdimn==3) CALL filter_single_particle_based_on_min_max_var(des_pos_new(lc1,3),'z-coordinate',z_min,z_max,.false.,keep_me)

! Statistical weight, radius and physical radius (CGDEM only)
! The diameter in particle_input.csv is the coarse-grained (not physical) diameter
         IF(CGDEM) THEN
            ! des_radius is the coarse grain radius, it's larger than physical radius
            IF(found_cgp_stat_wt) THEN
               des_cgp_stw(lc1) = dprocbuf(lbuf + cgp_stat_wt_start)
            ELSE
               des_cgp_stw(lc1) = CGP_STAT_WT(PIJK(LC1,5))
            ENDIF
            des_cgp_rpr(lc1) = des_radius(lc1)/(des_cgp_stw(lc1)**(1.0d0/3.0d0))
         ENDIF

! Particle Density
         IF(FOUND_ROL) THEN
            ro_sol(lc1) = dprocbuf(lbuf + Density_Start)
         ELSE
            ro_sol(lc1) = ro_s0(PIJK(LC1, 5))
         ENDIF
         CALL filter_single_particle_based_on_min_max_var(ro_sol(lc1),'Density',part_in_density_min,part_in_density_max,part_in_density_exclude,keep_me)

   ! Particle Velocity
         IF(FOUND_VEL) THEN
            des_vel_new(lc1,1:rdimn) = dprocbuf(lbuf + Velocity_start:lbuf + Velocity_start + rdimn - 1)
         ELSE
            des_vel_new(lc1,1:rdimn) = (/0.0, 0.0, 0.0/)
         ENDIF

         CALL filter_single_particle_based_on_min_max_var(des_vel_new(lc1,1),'x-velocity',part_in_u_min,part_in_u_max,part_in_u_exclude,keep_me)
         CALL filter_single_particle_based_on_min_max_var(des_vel_new(lc1,2),'y-velocity',part_in_v_min,part_in_v_max,part_in_v_exclude,keep_me)
         IF(rdimn==3) CALL filter_single_particle_based_on_min_max_var(des_vel_new(lc1,3),'z-velocity',part_in_w_min,part_in_w_max,part_in_w_exclude,keep_me)

   ! Particle Temperature
         IF(ENERGY_EQ) THEN
            IF(FOUND_TEMPERATURE) THEN
               des_t_s(lc1) = dprocbuf(lbuf + temperature_Start)
            ELSE
               des_t_s(lc1) = 273.15
            ENDIF
            IF(GSP_EXPLICIT) THEN
               gp_neighsid(lc1,:) = dprocbuf(lbuf+gp_neighsid_start:lbuf+gp_neighsid_start+5)
               gp_neighsa(lc1,:) = dprocbuf(lbuf+gp_neighsa_start:lbuf+gp_neighsa_start+5)
            ENDIF
            CALL filter_single_particle_based_on_min_max_var(des_t_s(lc1),'Temperature',part_in_temp_min,part_in_temp_max,part_in_temp_exclude,keep_me)
         ENDIF

   ! Particle Species: Always need DIMENSION_N_S values for all particles
         IF(ANY_SOLIDS_SPECIES_EQ) THEN
            IF(FOUND_SPECIES) THEN
               des_x_s(lc1,1:DIMENSION_N_S) = dprocbuf(lbuf + species_start:lbuf + species_start + DIMENSION_N_S - 1)
            ELSE
               des_x_s(lc1,1:DIMENSION_N_S) = X_S0(PIJK(LC1,5),1:DIMENSION_N_S)
            ENDIF
            DO n = 1,DIMENSION_N_S
               WRITE(label,'("Species",I2)') n
               CALL filter_single_particle_based_on_min_max_var(des_x_s(lc1,n),label,part_in_x_s_min(n),part_in_x_s_max(n),part_in_x_s_exclude(n),keep_me)
            ENDDO
         ENDIF

   ! Particle Scalar: CAUTION THIS ARRAY IS TRANSPOSED, FIRST INDEX IS USR_VAR ID, SECOND INDEX IS PARTICLE ID
         IF(DES_USR_VAR_SIZE>0) THEN
            IF(FOUND_USR_VAR) THEN
               des_usr_var(1:DES_USR_VAR_SIZE,lc1) = dprocbuf(lbuf + usr_var_start:lbuf + usr_var_start + DES_USR_VAR_SIZE - 1)
            ELSE
               des_usr_var(1:DES_USR_VAR_SIZE,lc1) = 1.0
            ENDIF
            DO n = 1,DES_USR_VAR_SIZE
               WRITE(label,'("User scalar # ",I2)') n
               CALL filter_single_particle_based_on_min_max_var(des_usr_var(n,lc1),label,part_in_usr_var_min(n),part_in_usr_var_max(n),part_in_usr_var_exclude(n),keep_me)
            ENDDO

         ENDIF

! Glued-sphere Dem parameters
         !@renjieke, for gsp_implicit need to save and distribute glueQuat and glueEX,EY,EZ
         IF(GSP_IMPLICIT) THEN
            IF(found_quat) THEN
               glueQuat(lc1,1:4) = dprocbuf(lbuf + quat_start:lbuf + quat_start + 4 - 1)
            ELSE
               glueQuat(lc1,1:4) = (/gsp_q1(PIJK(lc1,5)), gsp_q2(PIJK(lc1,5)), gsp_q3(PIJK(lc1,5)), gsp_q4(PIJK(lc1,5))/)
            ENDIF
            CALL gsp_quat_to_exyz(glueQuat(lc1,:),exone,eyone,ezone)
            glueEX(lc1,:) = exone
            glueEY(lc1,:) = eyone
            glueEZ(lc1,:) = ezone
         ENDIF

! Superquadric parameters
         IF(SuperDEM) THEN
            sqp_pid = PIJK(LC1,5)
            ! those sqp parameters won't have the "undefined" problems, they have default values.
            IF(found_semiaxis) THEN
               super_r(lc1,1:3) = dprocbuf(lbuf + superdem_semiaxis_start:lbuf + superdem_semiaxis_start + 3 - 1)
            ELSE
               super_r(lc1,1:3) = (/sqp_a(sqp_pid), sqp_b(sqp_pid), sqp_c(sqp_pid)/)
            ENDIF

            IF(found_roundness) THEN
               super_mn(lc1,1:2) = dprocbuf(lbuf + superdem_roundness_start:lbuf + superdem_roundness_start + 2 - 1)
            ELSE
               super_mn(lc1,1:2) = (/sqp_m(sqp_pid), sqp_n(sqp_pid)/)
            ENDIF

            IF(found_quat) THEN
               super_q(lc1,1:4) = dprocbuf(lbuf + superdem_quat_start:lbuf + superdem_quat_start + 4 - 1)
            ELSE
               super_q(lc1,1:4) = (/sqp_q1(sqp_pid), sqp_q2(sqp_pid), sqp_q3(sqp_pid), sqp_q4(sqp_pid)/)
            ENDIF

            ! Normalize quaternion
            q_norm = dsqrt(super_q(lc1,1)**2 + super_q(lc1,2)**2 &
                           +super_q(lc1,3)**2 + super_q(lc1,4)**2)

            IF(q_norm == zero) THEN
               WRITE(ERR_MSG,1215)
               CALL LOG_ERROR()
            ENDIF

 1215 FORMAT('Error reading particle input file.',&
           /'Quaternion has invalid value of zero.')

            super_q(lc1,1:4) = super_q(lc1,1:4) / q_norm

            IF(PARTICLE_ORIENTATION) THEN
            ! Convert quaternion to orientation vector
               CALL QROTATE(super_q(lc1,1:4),ORIENTATION(lc1,:),INIT_ORIENTATION,2)
            ENDIF
         ENDIF

   ! Keep or reject particle based on filtering
         IF(keep_me) THEN
            MAX_D = max(MAX_D,2*des_radius(lc1))
   ! Set particle status to "normal"
            m = PIJK(lc1,5)
            IF(des_rigid_motion(m)) THEN
               CALL set_rigid_motion(lc1)
            ELSE
               CALL set_normal(lc1)
            ENDIF
         ELSE ! end of keep_me
   ! Reject current particle and go back in the list of particles
            IF(GSP_EXPLICIT) THEN
               ! for GluedSphereDEM, do not go back and set lc1 and corresponding gid as non-existence
               CALL set_nonexistent(lc1)
               gid = gid_list(lc1)
               CALL set_nonexistent_gsp(gid)
               N_DELETED = N_DELETED + 1

               IF(numpes .gt. 1) THEN
                  dummy_index = findIndex_1i(MARKED_GSP,gid)
                  IF(dummy_index /= -1) THEN
                     marked_state(dummy_index) = particle_state_gsp(gid)
                  ELSE
                     CALL SET_MARKED_GSP(gid)
                     CALL SET_MARKED_GSP_STATE(particle_state_gsp(gid))
                  ENDIF
               ENDIF

            ELSE ! end of GluedSphereDEM
               lc1 = lc1 - 1
            ENDIF ! end of other DEM
         ENDIF ! end of not keep me

      ENDDO! end of do lcurpar = 1,pip
      deallocate(dprocbuf,drootbuf)

      call global_all_max(max_d)
      IF(mype == pe_io) THEN
         IF((XLENGTH/dble(DESGRIDSEARCH_IMAX)) < MAX_D) THEN
            WRITE(ERR_MSG, 1216) 'X', MAX_D, &
               XLENGTH/dble(DESGRIDSEARCH_IMAX)
            CALL LOG_INFO()
         ELSEIF((YLENGTH/dble(DESGRIDSEARCH_JMAX)) < MAX_D) THEN
            WRITE(ERR_MSG, 1216) 'Y', MAX_D, &
               YLENGTH/dble(DESGRIDSEARCH_JMAX)
            CALL LOG_INFO()
         ELSEIF((RDIMN == 3) .AND. ((ZLENGTH/dble(DESGRIDSEARCH_KMAX)) < MAX_D)) THEN
            WRITE(ERR_MSG, 1216) 'Z', MAX_D, &
               ZLENGTH/dble(DESGRIDSEARCH_KMAX)
            CALL LOG_INFO()
         ENDIF
      ENDIF
 1216 FORMAT('Warning: The des search grid is too fine in the ',A1, &
         '-direction. The',/'maximum particle diameter from CSV is larger ',&
         'than the cell width:',/2x,'MAX DIAM:   ',g12.5,/2x,'CELL ',      &
         'WIDTH: ',g12.5,/'Decrease the values for DESGRIDSEARCH in ', &
         'the project settings.')

! Set particle count after filtering
      pip = lc1
      max_pip = pip

      IF(GSP_EXPLICIT) THEN
         gsp_cnt = 0
         ! only exchange gsp_particle_state when it is DMP
         IF(NUMPES > 1) THEN
            CALL GLOBAL_ALL_SUM(N_DELETED)
            IF(N_DELETED>0) CALL UPDATE_GSP_PARTICLE_STATE
         ENDIF
         ! in DMP, through update_gsp_particle_state, pe_io should already have correct info about all the gsp states
         IF(myPE == PE_IO) THEN
            DO L = 1, NGluedParticles
               IF(PARTICLE_STATE_GSP(L) == 1) gsp_cnt = gsp_cnt + 1
            ENDDO
         ENDIF
         CALL BCAST(gsp_cnt)
      ENDIF

      CALL GLOBAL_ALL_SUM(PIP,FPIP)

   ! Print out the number of particles after filtering the data
      IF(GSP_EXPLICIT) THEN
         WRITE(ERR_MSG, 1217) NGluedParticles, gsp_cnt
         CALL LOG_INFO()
   ! Abort if no particles are left after filtering the data
         IF(gsp_cnt==0 .and. particles .ne. 0) THEN
            WRITE(ERR_MSG, 1218)
            CALL LOG_ERROR()
         ENDIF
      ELSE ! end of GSP_EXPLICIT
         WRITE(ERR_MSG, 1219)PARTICLES,FPIP
         CALL LOG_INFO()
   ! Abort if no particles are left after filtering the data
         IF(FPIP==0 .and. particles .ne. 0) THEN
            WRITE(ERR_MSG, 1220)
            CALL LOG_ERROR()
         ENDIF
      ENDIF ! end of other DEM

 1217  FORMAT('Number of glued-sphere particles from input  = ', I9 ,&
            /'Number of glued-sphere particles left after filtering the data = ', I9)

 1218  FORMAT('No glued-sphere particles left in the system after filtering the data.', &
            /'Please verify the filter settings.')

 1219  FORMAT('Number of particles read from particle_input.csv  = ', I9 ,&
            /'Number of particles left after filtering the data = ', I9)

 1220  FORMAT('No particles left in the system after filtering the data.', &
            /'Please verify the filter settings.')

 500  FORMAT(/2X,'From: DES_SCATTER_PARTICLE: (0)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')
 501  FORMAT(/2X,'From: DES_SCATTER_PARTICLE: (1)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')

      IF(myPE == PE_IO) deallocate(part_data,lpar_proc)
      IF(myPE == PE_IO .AND. GSP_EXPLICIT) deallocate(cpart_data)
      IF(myPE == PE_IO .AND. allocated(distribution_data)) deallocate(distribution_data)
      IF(myPE == PE_IO) CLOSE(lUNIT)

      IF(myPE == PE_IO) THEN
         IF(allocated(header_array)) deallocate(header_array)
         IF(allocated(intra_header_array)) deallocate(intra_header_array)
         IF(allocated(not_in_header_array)) deallocate(not_in_header_array)
      ENDIF

      RETURN
   END SUBROUTINE READ_PART_INPUT_V4P0

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
   Subroutine parse_scalar_option(lunit, Scalar_name, Scalar_condition, Default_scalar, Starting_index, n_cols, n_vars)

        USE make_upper_case_mod, only: make_upper_case

        IMPLICIT NONE
! Local unit
        integer, intent(in) :: lunit
! Scalar name
        character(LEN=*), intent(in) :: scalar_name
! Condition to read the option
        logical, intent(in) :: Scalar_condition
! Option to read or use default value (T/F)
        character(1) :: scalar_option
! Default value
        double precision, intent(out) :: Default_scalar
! Column index
        integer, intent(out) :: Starting_index
! Total number of columns
        integer, intent(inout) :: n_cols
! Total number of variables
        integer, intent(inout) :: n_vars
! Buffer
        character(32) Buffer
! IO Status:
        INTEGER :: IOS

        IOS            = 0
        Default_scalar = 0.0D0
        Starting_index = 0

        read(lunit,*) buffer, Scalar_option

        if(Scalar_condition) then

            CALL MAKE_UPPER_CASE (scalar_option,1)

            if(scalar_option=='F') then
                backspace(lUNIT)
                read(lunit,*,IOSTAT=IOS) buffer, buffer, Default_scalar
! Report read errors.
                IF(IOS > 0) THEN
                    WRITE(ERR_MSG,110) TRIM(Scalar_name)
                    CALL LOG_ERROR()
                ENDIF
                Starting_index   = 0
            elseif(scalar_option == 'T') then
                starting_index  = n_vars + 1
                n_vars          = starting_index
            else
                WRITE(ERR_MSG,120) trim(Scalar_name), Scalar_option
                CALL LOG_ERROR()
            endif

        endif

        return

 110 FORMAT('Error reading particle input file.',&
           /'Invalid or missing default value for : ',A)

 120 FORMAT('Error reading particle input file.',&
           /'Invalid option for : ',A,' = ',A,&
           /'Valid options are "T" or "F" ')
   end Subroutine parse_scalar_option

   Subroutine parse_vector_option(lunit, Vector_name, Vector_condition, Vector_size, Default_vector, Starting_index, n_cols, n_vars)
      USE make_upper_case_mod, only: make_upper_case

      IMPLICIT NONE

! Local unit
      integer, intent(in) :: lunit
! Vector name
      character(LEN=*), intent(in) :: Vector_name
! Condition to read the option
      logical, intent(in) :: Vector_condition
! Vector size
      integer, intent(in) :: Vector_size
! Option to read or use default value (T/F)
      character(1) :: vector_option
! Default value
      double precision, intent(out), dimension(Vector_size) :: Default_vector
! Column index
      integer, intent(out) :: Starting_index
! Total number of columns
      integer, intent(inout) :: n_cols
! Total number of variables
      integer, intent(inout) :: n_vars
! Buffer
      character(32) Buffer
! IO Status:
      INTEGER :: IOS

      IOS               = 0
      Default_vector(:) = 0.0D0
      Starting_index    = 0

      read(lunit,*) buffer, Vector_option

      if(Vector_name=='Coordinates'.AND.Vector_option/='T') then
         WRITE(ERR_MSG,110) trim(Vector_name), Vector_option
         CALL LOG_ERROR()
      endif


      if(Vector_condition) then

         CALL MAKE_UPPER_CASE (Vector_option,1)

         if(Vector_option=='F') then
            backspace(lUNIT)
            read(lunit,*,IOSTAT=IOS) buffer, buffer, Default_vector(1:Vector_size)
! Report read errors.
               IF(IOS > 0) THEN
                  WRITE(ERR_MSG,120) TRIM(Vector_name),Vector_size
                  CALL LOG_ERROR()
               ENDIF
            Starting_index   = 0
         elseif(Vector_option == 'T') then
            starting_index  = n_vars + 1
            n_vars          = starting_index + Vector_size - 1
         else
            WRITE(ERR_MSG,130) trim(Vector_name), Vector_option
            CALL LOG_ERROR()
         endif

      endif

      return

 110 FORMAT('Error reading particle input file.',&
           /'Invalid option for : ',A,' = ',A,&
           /'Only valid option is "T" ')

 120 FORMAT('Error reading particle input file.',&
           /'Invalid or missing default value(s) for : ',A,&
           /'Expected number of components = ',I2)

 130 FORMAT('Error reading particle input file.',&
           /'Invalid option for : ',A,' = ',A,&
           /'Valid options are "T" or "F" ')
      end Subroutine parse_vector_option


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: READ_PART_INPUT_V1P0                                      !
!                                                                      !
! Purpose: Read the particle input and broadcasts the particle data to !
! respective processors.                                               !
! Version 1.0 (legacy)                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_PART_INPUT_V1P0

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      use cdist
      use compar
      use desmpi
      use functions
      use funits
      use geometry, only: NO_K
      use mpi_init_des, only: des_scatter_particle
      use mpi_utility

      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: k
! index of particle
      INTEGER :: lcurpar
! local unit
      INTEGER, PARAMETER :: lunit=10
! local filename
      character(255) lfilename
! IO Status:
      INTEGER :: IOS
! Flag to indicate if file exists.
      LOGICAL :: lEXISTS
! Read dimension: 2D vs 3D data
      integer :: RDMN
!-----------------------------------------------

      IOS = 0
      RDMN = merge(2,3,NO_K)


! Read the file
!----------------------------------------------------------------->>>
! In distributed IO the first line of the file will be number of
! particles in that processor
      IF (bdist_io) then
         read(lunit,*) pip
         DO lcurpar = 1,pip
            call set_normal(lcurpar)
            read (lunit,*) (des_pos_new(lcurpar,k),k=1,RDMN),&
               des_radius(lcurpar), ro_sol(lcurpar),&
               (des_vel_new(lcurpar,k),k=1,RDMN)
         ENDDO

! Serial IO (not bDIST_IO)
      ELSE
!----------------------------------------------------------------->>>

! Read into temporary variable and scatter
         IF (myPE .eq. PE_IO) THEN

! Allocate and initialize temporary variables.
            ALLOCATE (dpar_pos(particles,3)); dpar_pos=0.0
            ALLOCATE (dpar_vel(particles,3)); dpar_vel=0.0
            ALLOCATE (dpar_rad(particles));   dpar_rad=0.0
            ALLOCATE (dpar_den(particles));   dpar_den = 0.0
! Loop through the input file.
            DO lcurpar = 1, particles
               read (lunit,*,IOSTAT=IOS)                               &
               (dpar_pos(lcurpar,k),k=1,RDMN),dpar_rad(lcurpar),       &
               dpar_den(lcurpar),(dpar_vel(lcurpar,k),k=1,RDMN)

! Report read errors.
               IF(IOS > 0) THEN
                  WRITE(ERR_MSG,1200)
                  CALL LOG_ERROR()
                  EXIT
 1200 FORMAT('Error reading particle input ',&
         'file.',/'A common error is 2D input for 3D cases.')

! Report End of file errors.
               ELSEIF(IOS < 0) THEN
                  WRITE(ERR_MSG,1201) &
                     trim(iVal(lcurpar)), trim(iVal(Particles))
                  CALL LOG_ERROR()
                  EXIT
 1201 FORMAT('Error reading particle input ',&
         'file.',/'End of file found for particle ',A,' and ',A,1X,    &
         'entries are expected.')

               ENDIF

            ENDDO

         ENDIF

         CALL GLOBAL_ALL_SUM(IOS)
         IF(IOS /= 0) CALL LOG_ERROR()

         CALL DES_SCATTER_PARTICLE

         IF(myPE == PE_IO) &
            deallocate (dpar_pos,dpar_vel,dpar_rad,dpar_den)

      ENDIF   ! end if/else bdist_io


      IF(bDIST_IO .OR. myPE == PE_IO) CLOSE(lUNIT)
      RETURN
 END SUBROUTINE READ_PART_INPUT_V1P0

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: READ_PART_INPUT_SuperDEM                                  !
!                                                                      !
! Purpose: Read the particle input and broadcasts the particle data to !
! respective processors.                                               !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_PART_INPUT_SuperDEM

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      use cdist
      use compar
      use desmpi
      use functions
      use funits
      use geometry, only: NO_K
      use mpi_init_des, only: des_scatter_particle
      use mpi_utility

      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      integer :: k
! index of particle
      INTEGER :: lcurpar
! local unit
      INTEGER, PARAMETER :: lunit=10
! local filename
      character(255) lfilename
! IO Status:
      INTEGER :: IOS
! Flag to indicate if file exists.
      LOGICAL :: lEXISTS
! Read dimension: 2D vs 3D data
      integer :: RDMN
!-----------------------------------------------

      IOS = 0
      RDMN = merge(2,3,NO_K)

! Setup the file name based on distributed or serial IO.
      IF(bDIST_IO) THEN
         lFILENAME = ''
         WRITE(lFILENAME,'("particle_input_",I4.4,".dat")') myPE
      ELSE
         lFILENAME= "particle_input.dat"
      ENDIF

! Check the the file exists and open it.
      IF(bDIST_IO .OR. myPE == PE_IO) THEN
         INQUIRE(FILE=lFILENAME, EXIST=lEXISTS)
         IF(.NOT.LEXISTS) THEN
            WRITE(ERR_MSG, 1100)
            CALL LOG_ERROR()
            IOS = 1
         ELSE
            OPEN(UNIT=lUNIT, FILE=lFILENAME, FORM="FORMATTED")
         ENDIF
      ENDIF

! Collect the error message and quit.
      CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) call LOG_ERROR()

 1100 FORMAT('FATAL - DEM particle input file not found.')

! Read the file
!----------------------------------------------------------------->>>
! In distributed IO the first line of the file will be number of
! particles in that processor
      IF (bdist_io) then
         read(lunit,*) pip
         DO lcurpar = 1,pip
            call set_normal(lcurpar)
            read (lunit,*) (des_pos_new(lcurpar,k),k=1,RDMN),&
               des_radius(lcurpar), ro_sol(lcurpar),&
               (des_vel_new(lcurpar,k),k=1,RDMN)
         ENDDO

! Serial IO (not bDIST_IO)
      ELSE
!----------------------------------------------------------------->>>

! Read into temporary variable and scatter
         IF (myPE .eq. PE_IO) THEN

! Allocate and initialize temporary variables.
            ALLOCATE (dpar_pos(particles,3)); dpar_pos=0.0
            ALLOCATE (dpar_vel(particles,3)); dpar_vel=0.0
            ALLOCATE (dpar_rad(particles));   dpar_rad=0.0
            ALLOCATE (dpar_den(particles));   dpar_den = 0.0
!SuperDEM
            IF(SuperDEM) THEN
               ALLOCATE (dpar_super_r(particles,3));   dpar_super_r = 0.0
               ALLOCATE (dpar_super_mn(particles,2));   dpar_super_mn = 0.0
               ALLOCATE (dpar_super_q(particles,4));   dpar_super_q = 0.0
               ALLOCATE (dpar_var(des_usr_var_size,particles));   dpar_var = 0.0
            ENDIF
! Loop through the input file.
            DO lcurpar = 1, particles
               read (lunit,*,IOSTAT=IOS)                               &
               (dpar_pos(lcurpar,k),k=1,RDMN),dpar_rad(lcurpar),       &
               dpar_den(lcurpar),(dpar_vel(lcurpar,k),k=1,RDMN)

! Report read errors.
               IF(IOS > 0) THEN
                  WRITE(ERR_MSG,1200)
                  CALL LOG_ERROR()
                  EXIT
 1200 FORMAT('Error reading particle input ',&
         'file.',/'A common error is 2D input for 3D cases.')

! Report End of file errors.
               ELSEIF(IOS < 0) THEN
                  WRITE(ERR_MSG,1201) &
                     trim(iVal(lcurpar)), trim(iVal(Particles))
                  CALL LOG_ERROR()
                  EXIT
 1201 FORMAT('Error reading particle input ',&
         'file.',/'End of file found for particle ',A,' and ',A,1X,    &
         'entries are expected.')

               ENDIF

            ENDDO



            open(unit=1091,file='sq_input.dat', action='read')
            do lcurpar = 1, particles
               read(1091,*) dpar_super_r(lcurpar,1),dpar_super_r(lcurpar,2),&
                          dpar_super_r(lcurpar,3),dpar_super_mn(lcurpar,1),&
                          dpar_super_mn(lcurpar,2), dpar_super_q(lcurpar,1),&
                          dpar_super_q(lcurpar,2),dpar_super_q(lcurpar,3),&
                          dpar_super_q(lcurpar,4)
              enddo
              close(1091)

       IF(myPE == PE_IO)  WRITE(*,'(/," Finished reading superquadric particle data.")')

         ENDIF  !mype_io

         CALL GLOBAL_ALL_SUM(IOS)
         IF(IOS /= 0) CALL LOG_ERROR()

         CALL DES_SCATTER_PARTICLE

         IF(myPE == PE_IO) &
            deallocate (dpar_pos,dpar_vel,dpar_rad,dpar_den,dpar_super_r,dpar_super_mn,dpar_super_q,dpar_var)



      ENDIF   ! end if/else bdist_io

      IF(bDIST_IO .OR. myPE == PE_IO) CLOSE(lUNIT)
      RETURN
   END SUBROUTINE READ_PART_INPUT_SuperDEM

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: READ_SPHERE_INPUT_GluedSphereDEM                         !
!                                                                      !
! Author: Hang Zhou                                      Date: Dec 2023!
! Revised: Renjie Ke                                     Date: Mar 2024!
! Purpose: Read the sphere particle input                              !
!               and broadcasts the particle data to                    !
!               respective processors for glued-sphere models          !
!               only read in pe_io, scatter later in READ_PART_INPUT    !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_SPHERE_INPUT_GluedSphereDEM

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use discretelement, only: sphereNumber, sphereRelDia, sphereRelPos, start_index_sphere, gsp_vol
      use discretelement, only: tmp_gp_sa, tmp_gp_neighsid, tmp_gp_neighsa, gp_sa, gp_neighsid, gp_neighsa
      use discretelement, only: des_mmax, pbb
      use geometry, only: NO_K
      USE param1, only: undefined, zero
      use physprop, only: d_p0
      use mpi_utility, only: bcast
      use constant, only: pi

      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------

! local filename
      character(255) lfilename
      character(255) lfilename1
! Flag to indicate if file exists.
      LOGICAL :: lEXISTS
      LOGICAL :: lEXISTS1

! indices
      integer :: i,j,k,ilines
! index of spheres
      integer :: lcurpar
! index of solid phase
      integer :: MM
! IO Status:
      integer :: IOS, IOS1
! local unit
      integer, parameter :: lunit=10
      integer, parameter :: lunit1=11
! Read dimension: 2D vs 3D data
      integer :: RDIMN
! Buffer
      character(32) Buffer
! Integer buffer
      integer :: int_buffer
! Line
      character(200) Line
      integer :: ppos,linepos
! Length of list of variables to read (number of columns)
      integer :: n_vars
! Index of solid phase
      integer :: phase_index
! Number of solid phases using glued-sphere model. All others are assumed sphere.
      integer, dimension(:), allocatable :: idGluedPhase ! EX idGluedPhase = /1,3,5/
! Number of spheres in each non-spherical particle of each phase
      integer, dimension(:), allocatable :: sphereNumber_tmp
! Phase number
      integer :: PN
! Sphere index while reading file
      integer :: sphere_index
! GSP extent (including sphere size)
      double precision :: gsp_min(3), gsp_max(3)
!-----------------------------------------------
      IOS = 0
      RDIMN = merge(2,3,NO_K)
! Check the the file exists and open it.
      lfilename = "gsp_config.dat"
      lfilename1 = "gsp_config.csv"
      IF(myPE == PE_IO) THEN
         INQUIRE(FILE=LFILENAME, EXIST=lEXISTS)
         IF(.NOT.LEXISTS) THEN
            WRITE(ERR_MSG, 1100)
            CALL LOG_ERROR()
            IOS = 1
         ELSE
            INQUIRE(FILE=LFILENAME1, EXIST=lEXISTS1)
            IF(.NOT. lEXISTS1) CALL GSP_CONFIG_DAT_TO_CSV
            OPEN(UNIT=lUNIT, FILE=lFILENAME, FORM="FORMATTED")
         ENDIF
      ENDIF

! Collect the error message and quit.
      CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) call LOG_ERROR()

 1100 FORMAT('FATAL - gsp_config.dat not found.')

      IF (myPE .eq. PE_IO) THEN
! Skip Instructions
         ilines = 14
         do i = 1, ilines
            read(lunit,*) buffer
         enddo

! Dimension
         read(lunit,*) buffer,int_buffer
         if(RDIMN/=int_buffer) then
            WRITE(ERR_MSG,100) int_buffer, RDIMN
            CALL LOG_ERROR()
         endif
         n_vars = int_buffer+1
! ID of Solid Phases using glued-sphere model
         read(lunit,'(a)') line
         linepos = index(line,':')
         line = line(linepos+1:)
         ! Number of solid phases given in the input file
         i = count(transfer(trim(adjustl(line)), 'a', len(trim(adjustl(line))))==" ")
         ALLOCATE(idGluedPhase(i+1))
         ALLOCATE(sphereNumber_tmp(i+1))
         read(line, *) idGluedPhase

! Number of Spheres in One Non-Spherical Particle in EachPhase
         read(lunit,*) buffer,sphereNumber_tmp

! Skip header
         do i = 1, 3
            read(lunit,*) buffer
         enddo

 100 FORMAT('Error reading sphere particle input file.',&
           /'Dimension (2D or 3D) does not match simulation setup file.',&
           /'Dimension in sphere_particle_input.dat    = ',I3,&
           /'Dimension in simulation setup file = ',I3)
      ENDIF !  IF (myPE .eq. PE_IO)

      ! even use read particle_input.dat, still need to specify each phase in GUI to get the proper DES_MMAX
      Allocate(sphereNumber(DES_MMAX)); sphereNumber = 1

! Read the file. Distributed I/O is not supported
      IF (myPE .eq. PE_IO) THEN
! Set the number of sphere in each solid phase (The initial value is 1, assuming the particle is sphere)
         do MM = 1, size(idGluedPhase) ! EX if idGluedPhase = (1,3,5), so then MM iterates from 1 to 3
            sphereNumber(idGluedPhase(MM)) = sphereNumber_tmp(MM) !
         enddo
      ENDIF

! Allocate the relateive position and diameter of spheres to its host non-spherical particle in each phase
      CALL BCAST(sphereNumber,pe_io)
      ALLOCATE(sphereRelPos(sum(sphereNumber),RDIMN)); sphereRelPos=0.0D0
      ALLOCATE(sphereRelDia(sum(sphereNumber))); sphereRelDia=1.0D0
! Also allocate the surface area and neigh component info
         ! they are also relative values
         ! sa * bounding diameter^2
         ! neighsid will relative to the first component sphere in each solid phase
         ! neighsa * bounding diameter (it include eff_dist)
      ALLOCATE(tmp_gp_sa(sum(sphereNumber))); tmp_gp_sa=0.0D0
      ALLOCATE(tmp_gp_neighsid(sum(sphereNumber),6)); tmp_gp_neighsid=0
      ALLOCATE(tmp_gp_neighsa(sum(sphereNumber),6)); tmp_gp_neighsa=0.0D0

      ALLOCATE(gsp_vol(DES_MMAX)); gsp_vol = zero

      IF (myPE .eq. PE_IO) THEN
! Get the start index of each solid phase in sphereRelPos and sphereRelDia
! DES_MMAX: number of discrete 'solids phases'
! Initialize gsp_vol to that of a single sphere of diameter d_p0,
! so that the GSP model can handle GSPs made up of one sphere
! (GSPs that are not included in gsp_config.dat are considered
! single sphere GSPs).

         DO MM = 1, DES_MMAX
            gsp_vol(MM) = (pi*d_p0(MM)**3.0D0)/6.0d0
            IF (MM .eq. 1) then
               start_index_sphere(MM) = 1
            ELSE
               start_index_sphere(MM) = start_index_sphere(MM-1) + sphereNumber(MM-1)
            ENDIF
         ENDDO

! Loop through the input file.
         DO MM = 1, size(idGluedPhase)
            phase_index = idGluedPhase(MM)
            gsp_vol(MM) = zero
            gsp_min(:) =   UNDEFINED
            gsp_max(:) = - UNDEFINED
            DO lcurpar = 1, sphereNumber(phase_index)
               sphere_index = start_index_sphere(phase_index) + lcurpar - 1
               read (lunit,*,IOSTAT=IOS) (sphereRelPos(sphere_index, ppos), ppos=1,RDIMN), &
                                          sphereRelDia(sphere_index), &
                                          tmp_gp_sa(sphere_index), &
                                          (tmp_gp_neighsid(sphere_index,i),i=1,6), &
                                          (tmp_gp_neighsa(sphere_index,j),j=1,6)

               gsp_min(:) = min(gsp_min(:),sphereRelPos(sphere_index,:) - 0.5*sphereRelDia(sphere_index))
               gsp_max(:) = max(gsp_max(:),sphereRelPos(sphere_index,:) + 0.5*sphereRelDia(sphere_index))
               gsp_vol(MM) = gsp_vol(MM)+pi/6.0D0*(d_p0(MM)*sphereRelDia(sphere_index))**3.0D0

! Report read errors.
               IF(IOS > 0) THEN
                  WRITE(ERR_MSG,1200)
                  CALL LOG_ERROR()
                  EXIT
 1200 FORMAT('Error reading gsp_config.dat.',&
         /'A common error is 2D input for 3D cases.')
! Report End of file errors.
               ELSEIF(IOS < 0) THEN
                  WRITE(ERR_MSG,1202) &
                        trim(iVal(lcurpar)), trim(iVal(sum(sphereNumber_tmp)))
                  CALL LOG_ERROR()
                  EXIT
 1202 FORMAT('Error reading gsp_config.dat.',&
         /'End of file found for particle ',A,' and ',A,1X,    &
         'entries are expected.')
               ENDIF

            ENDDO ! END of lcurpar
            ! particle bounding box
            pbb(phase_index,:) = gsp_max(:) - gsp_min(:)
         ENDDO ! END of MM

         ! set tmp_gp_neighsid relative to the first component spheres for auto seeding
         where (tmp_gp_neighsid /= -1) tmp_gp_neighsid = tmp_gp_neighsid - 1
      ENDIF ! (myPE == PE_IO)

      CALL BCAST(start_index_sphere,pe_io)
      CALL BCAST(sphereRelPos,pe_io)
      CALL BCAST(sphereRelDia,pe_io)
      CALL BCAST(tmp_gp_sa, pe_io)
      CALL BCAST(tmp_gp_neighsid, pe_io)
      CALL BCAST(tmp_gp_neighsa, pe_io)
      CALL BCAST(gsp_vol, pe_io)
      CALL BCAST(pbb, pe_io)

      CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) CALL LOG_ERROR()

      IF(myPE == PE_IO) THEN
         CLOSE(lUNIT)
      ENDIF
      RETURN
   END SUBROUTINE READ_SPHERE_INPUT_GluedSphereDEM

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: READ_SPHERE_INPUT_GSP_NO_INTRA                           !
!                                                                      !
! Author: Renjie Ke                                      Date: Nov 2024!
! Purpose: Read the sphere particle input                              !
!               and broadcasts the particle data to                    !
!               respective processors for glued-sphere models          !
!               only read in pe_io, scatter later in READ_PART_INPUT    !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_SPHERE_INPUT_GSP_NO_INTRA

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use discretelement, only: sphereNumber, sphereRelDia, sphereRelPos, start_index_sphere, gsp_vol
      use discretelement, only: tmp_gp_sa, tmp_gp_neighsid, tmp_gp_neighsa, gp_sa, gp_neighsid, gp_neighsa
      use discretelement, only: des_mmax, pbb, cspart_data
      use geometry, only: NO_K
      USE param1, only: undefined, zero
      use physprop, only: d_p0
      use mpi_utility, only: bcast
      use constant, only: pi

      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------

! local filename
      character(255) lfilename
      character(255) lfilename1
! Flag to indicate if file exists.
      LOGICAL :: lEXISTS
      LOGICAL :: lEXISTS1

! indices
      integer :: i,j,k,ilines
! index of spheres
      integer :: lcurpar
! index of solid phase
      integer :: MM
! IO Status:
      integer :: IOS, IOS1
! local unit
      integer, parameter :: lunit=10
      integer, parameter :: lunit1=11
! Read dimension: 2D vs 3D data
      integer :: RDIMN
! Buffer
      character(32) Buffer
! Integer buffer
      integer :: int_buffer
! Line
      character(200) Line
      integer :: ppos,linepos
! Length of list of variables to read (number of columns)
      integer :: n_vars
! Index of solid phase
      integer :: phase_index
! Number of solid phases using glued-sphere model. All others are assumed sphere.
      integer, dimension(:), allocatable :: idGluedPhase ! EX idGluedPhase = /1,3,5/
! Number of spheres in each non-spherical particle of each phase
      integer, dimension(:), allocatable :: sphereNumber_tmp
! Phase number
      integer :: PN
! GSP extent (including sphere size)
      double precision :: gsp_min(3), gsp_max(3)
! The max phaseid in idGluedPhase
      integer :: maxval_idGluedPhase
!-----------------------------------------------
      IOS = 0
      RDIMN = merge(2,3,NO_K)
! Check the the file exists and open it.
      lfilename = "gsp_config.dat"
      lfilename1 = "gsp_config.csv"
      IF(myPE == PE_IO) THEN
         INQUIRE(FILE=LFILENAME, EXIST=lEXISTS)
         IF(.NOT.LEXISTS) THEN
            WRITE(ERR_MSG, 1100)
            CALL LOG_ERROR()
            IOS = 1
         ELSE
            INQUIRE(FILE=LFILENAME1, EXIST=lEXISTS1)
            IF(.NOT. lEXISTS1) CALL GSP_CONFIG_DAT_TO_CSV
            OPEN(UNIT=lUNIT, FILE=lFILENAME, FORM="FORMATTED")
         ENDIF
      ENDIF

! Collect the error message and quit.
      CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) call LOG_ERROR()

 1100 FORMAT('FATAL - gsp_config.dat not found.')

      IF (myPE .eq. PE_IO) THEN
! Skip Instructions
         ilines = 14
         do i = 1, ilines
            read(lunit,*) buffer
         enddo

! Dimension
         read(lunit,*) buffer,int_buffer
         if(RDIMN/=int_buffer) then
            WRITE(ERR_MSG,100) int_buffer, RDIMN
            CALL LOG_ERROR()
         endif
         n_vars = int_buffer+1
! ID of Solid Phases using glued-sphere model
         read(lunit,'(a)') line
         linepos = index(line,':')
         line = line(linepos+1:)
         ! Number of solid phases given in the input file
         i = count(transfer(trim(adjustl(line)), 'a', len(trim(adjustl(line))))==" ")
         ALLOCATE(idGluedPhase(i+1))
         ALLOCATE(sphereNumber_tmp(i+1))
         read(line, *) idGluedPhase

! Number of Spheres in One Non-Spherical Particle in EachPhase
         read(lunit,*) buffer,sphereNumber_tmp

! Skip header
         do i = 1, 3
            read(lunit,*) buffer
         enddo

 100 FORMAT('Error reading sphere particle input file.',&
           /'Dimension (2D or 3D) does not match simulation setup file.',&
           /'Dimension in sphere_particle_input.dat    = ',I3,&
           /'Dimension in simulation setup file = ',I3)
      ENDIF !  IF (myPE .eq. PE_IO)

      ! even use read particle_input.dat, still need to specify each phase in GUI to get the proper DES_MMAX
      Allocate(sphereNumber(DES_MMAX)); sphereNumber = 1

! Read the file. Distributed I/O is not supported
      IF (myPE .eq. PE_IO) THEN
! Set the number of sphere in each solid phase (The initial value is 1, assuming the particle is sphere)
         do MM = 1, size(idGluedPhase) ! EX if idGluedPhase = (1,3,5), so then MM iterates from 1 to 3
            sphereNumber(idGluedPhase(MM)) = sphereNumber_tmp(MM) !
         enddo
         maxval_idGluedPhase = maxval(idGluedPhase)
      ENDIF

! Allocate the relateive position and diameter of spheres to its host non-spherical particle in each phase
      CALL BCAST(sphereNumber,pe_io)
      ! ALLOCATE(sphereRelPos(sum(sphereNumber),RDIMN)); sphereRelPos=0.0D0
      ! ALLOCATE(sphereRelDia(sum(sphereNumber))); sphereRelDia=1.0D0

      ! CALL BCAST(n_var,pe_io)
      CALL BCAST(maxval_idGluedPhase, pe_io)
      ! first 5 column should be enough for no intra heat transfer case
      Allocate(cspart_data(maxval(sphereNumber),5,maxval_idGluedPhase))
      cspart_data(:,:,:) = 0.0

      ALLOCATE(gsp_vol(DES_MMAX)); gsp_vol = zero

      IF (myPE .eq. PE_IO) THEN
! start_index_sphere(M) tells start position of a phase
         DO MM = 1, DES_MMAX
            gsp_vol(MM) = (pi*d_p0(MM)**3.0D0)/6.0d0
         ENDDO

! Loop through the input file.
         DO MM = 1, size(idGluedPhase)
            phase_index = idGluedPhase(MM)
            gsp_vol(MM) = zero
            gsp_min(:) =   UNDEFINED
            gsp_max(:) = - UNDEFINED
            DO lcurpar = 1, sphereNumber(phase_index)
               read (lunit,*,IOSTAT=IOS) cspart_data(lcurpar,:,phase_index)
               gsp_min(:) = min(gsp_min(:),cspart_data(lcurpar,1:3,phase_index) - 0.5*cspart_data(lcurpar,4,phase_index))
               gsp_max(:) = max(gsp_max(:),cspart_data(lcurpar,1:3,phase_index) + 0.5*cspart_data(lcurpar,4,phase_index))
               gsp_vol(MM) = gsp_vol(MM)+pi/6.0D0*(d_p0(MM)*cspart_data(lcurpar,4,phase_index))**3.0D0

! Report read errors.
               IF(IOS > 0) THEN
                  WRITE(ERR_MSG,1200)
                  CALL LOG_ERROR()
                  EXIT
 1200 FORMAT('Error reading gsp_config.dat.',&
         /'A common error is 2D input for 3D cases.')
! Report End of file errors.
               ELSEIF(IOS < 0) THEN
                  WRITE(ERR_MSG,1202) &
                        trim(iVal(lcurpar)), trim(iVal(sum(sphereNumber_tmp)))
                  CALL LOG_ERROR()
                  EXIT
 1202 FORMAT('Error reading gsp_config.dat.',&
         /'End of file found for particle ',A,' and ',A,1X,    &
         'entries are expected.')
               ENDIF

            ENDDO ! END of lcurpar
            ! particle bounding box
            pbb(phase_index,:) = gsp_max(:) - gsp_min(:)
         ENDDO ! END of MM
      ENDIF ! (myPE == PE_IO)

      CALL BCAST(cspart_data, pe_io)
      CALL BCAST(gsp_vol, pe_io)
      CALL BCAST(pbb, pe_io)

      CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) CALL LOG_ERROR()

      IF(myPE == PE_IO) THEN
         CLOSE(lUNIT)
      ENDIF
      RETURN
   END SUBROUTINE READ_SPHERE_INPUT_GSP_NO_INTRA

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: READ_SPHERE_INPUT_GSP                                    !
!                                                                      !
! Author: Renjie Ke                                      Date: Jan 2025!
! Purpose: One routine work for gsp model implicit and explicit        !
!          Also read csv format of gsp_config only                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE READ_SPHERE_INPUT_GSP
!-----------------------------------------------
! Modules
!-----------------------------------------------
      use discretelement, only: sphereNumber, sphereRelDia, sphereRelPos, start_index_sphere, gsp_vol
      use discretelement, only: tmp_gp_sa, tmp_gp_neighsid, tmp_gp_neighsa, gp_sa, gp_neighsid, gp_neighsa
      use discretelement, only: des_mmax, pbb, cspart_data, GSP_EXPLICIT, GSP_IMPLICIT
      use geometry, only: NO_K
      USE param1, only: undefined, zero
      use physprop, only: d_p0
      use mpi_utility, only: bcast
      use constant, only: pi

      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local filename
      character(255) lfilename
      character(255) lfilename1
! Flag to indicate if file exists.
      LOGICAL :: lEXISTS
      LOGICAL :: lEXISTS1
! indices
      integer :: i,j,lcurpar
! index of solid phase
      integer :: MM
! IO Status:
      integer :: IOS
! local unit
      integer, parameter :: lunit1=11
! Read dimension: 2D vs 3D data
      integer :: RDIMN
! Line
      integer :: ppos
      CHARACTER(1024) :: vline
      CHARACTER(1024) :: hline
! Length of list of variables to read (number of columns)
      integer :: n_vars
! Sphere index while reading file
      integer :: sphere_index
! GSP extent (including sphere size)
      double precision :: gsp_min(3), gsp_max(3)
! Count the csv file
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: data_line
      integer :: old_pid, new_pid, phase_cnt, maxval_idGluedPhase

!-----------------------------------------------

      IOS = 0
      RDIMN = merge(2,3,NO_K)
! Check the the file exists and open it.
      lfilename = "gsp_config.dat"
      lfilename1 = "gsp_config.csv"
      IF(myPE == PE_IO) THEN
         INQUIRE(FILE=LFILENAME, EXIST=lEXISTS)
         INQUIRE(FILE=LFILENAME1, EXIST=lEXISTS1)

         IF(.NOT. LEXISTS .AND. .NOT. LEXISTS1) THEN
            WRITE(ERR_MSG, 1100)
            CALL LOG_ERROR()
            IOS = 1
         ELSEIF(LEXISTS .AND. .NOT. LEXISTS1) THEN
            ! The size of array such as
            ! sphereNumber, sphereRelDia, sphereRelPos, start_index_sphere, gsp_vol
            ! should also be set in this migration code
            CALL GSP_CONFIG_DAT_TO_CSV
            INQUIRE(FILE=LFILENAME1, EXIST=lEXISTS1)
            IF(.NOT. lEXISTS1) THEN
               WRITE(ERR_MSG, 1101)
               CALL LOG_ERROR()
            ENDIF
         ENDIF

         OPEN(UNIT=lUNIT1, FILE=lFILENAME1, status='old', action='READ')
         ! gsp_config.csv has the fixed number of columns
         n_vars = 18
         READ(LUNIT1, '(A)', IOSTAT=IOS) vline
         READ(LUNIT1, '(A)', IOSTAT=IOS) hline
         allocate(data_line(n_vars))
      ENDIF ! end of mype == pe_io

! Collect the error message and quit.
      CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) call LOG_ERROR()

 1100 FORMAT('FATAL - neither gsp_config.dat nor gsp_config.csv found.')
 1101 FORMAT('FATAL - migration from dat to csv failed.')

      allocate(sphereNumber(DES_MMAX)); sphereNumber = 0

      IF(myPE == PE_IO) THEN
         ! read the whole gsp_config.csv and count the sphereNumber
         DO
            READ(LUNIT1,*,IOSTAT = IOS) data_line(1:n_vars)
            IF(ios /= 0) exit ! end of file
            new_pid = data_line(n_vars)
            sphereNumber(new_pid) = sphereNumber(new_pid) + 1
         ENDDO

         CLOSE(lUNIT1)

         phase_cnt = 0
         maxval_idGluedPhase = -1
         DO MM = 1, size(sphereNumber)
            IF(sphereNumber(MM) == 0) CYCLE
            phase_cnt = phase_cnt + 1
            maxval_idGluedPhase = max(phase_cnt, maxval_idGluedPhase)
         ENDDO
      ENDIF ! end of mype == pe_io

      call bcast(sphereNumber, pe_io)

      IF(GSP_EXPLICIT) THEN
         ALLOCATE(sphereRelPos(sum(sphereNumber),RDIMN)); sphereRelPos=0.0D0
         ALLOCATE(sphereRelDia(sum(sphereNumber))); sphereRelDia=1.0D0
         ALLOCATE(tmp_gp_sa(sum(sphereNumber))); tmp_gp_sa=0.0D0
         ALLOCATE(tmp_gp_neighsid(sum(sphereNumber),6)); tmp_gp_neighsid=0
         ALLOCATE(tmp_gp_neighsa(sum(sphereNumber),6)); tmp_gp_neighsa=0.0D0
      ELSE
         call bcast(maxval_idGluedPhase, pe_io)
         Allocate(cspart_data(maxval(sphereNumber),5,maxval_idGluedPhase))
         cspart_data(:,:,:) = 0.0
      ENDIF
      ALLOCATE(gsp_vol(DES_MMAX)); gsp_vol = zero

      ! At this point, we are not thinking mix of gsp model with other model
      IF(myPE == PE_IO) THEN
         OPEN(UNIT=lUNIT1, FILE=lFILENAME1, status='old', action='READ')
         READ(LUNIT1, '(A)', IOSTAT=IOS) vline
         READ(LUNIT1, '(A)', IOSTAT=IOS) hline
         IF(GSP_EXPLICIT) THEN
            DO MM = 1, DES_MMAX
               IF (MM .eq. 1) then
                  start_index_sphere(MM) = 1
               ELSE
                  start_index_sphere(MM) = start_index_sphere(MM-1) + sphereNumber(MM-1)
               ENDIF
            ENDDO
         ENDIF
         DO MM = 1, phase_cnt
            gsp_vol(MM) = zero
            gsp_min(:) =   UNDEFINED
            gsp_max(:) = - UNDEFINED

            DO lcurpar = 1, sphereNumber(MM)
               IF(GSP_EXPLICIT) THEN
                  sphere_index = start_index_sphere(MM) + lcurpar - 1
                  read (lunit1,*,IOSTAT=IOS) (sphereRelPos(sphere_index, ppos), ppos=1,RDIMN), &
                                             sphereRelDia(sphere_index), &
                                             tmp_gp_sa(sphere_index), &
                                             (tmp_gp_neighsid(sphere_index,i),i=1,6), &
                                             (tmp_gp_neighsa(sphere_index,j),j=1,6)
                  gsp_min(:) = min(gsp_min(:),sphereRelPos(sphere_index,:) - 0.5*sphereRelDia(sphere_index))
                  gsp_max(:) = max(gsp_max(:),sphereRelPos(sphere_index,:) + 0.5*sphereRelDia(sphere_index))
                  gsp_vol(MM) = gsp_vol(MM)+pi/6.0D0*(d_p0(MM)*sphereRelDia(sphere_index))**3.0D0
               ELSE
                  read (lunit1,*,IOSTAT=IOS) cspart_data(lcurpar,:,MM)
                  gsp_min(:) = min(gsp_min(:),cspart_data(lcurpar,1:3,MM) - 0.5*cspart_data(lcurpar,4,MM))
                  gsp_max(:) = max(gsp_max(:),cspart_data(lcurpar,1:3,MM) + 0.5*cspart_data(lcurpar,4,MM))
                  gsp_vol(MM) = gsp_vol(MM)+pi/6.0D0*(d_p0(MM)*cspart_data(lcurpar,4,MM))**3.0D0
               ENDIF

               ! Report read errors.
               IF(IOS > 0) THEN
                  WRITE(ERR_MSG,1200)
                  CALL LOG_ERROR()
                  EXIT
 1200 FORMAT('Error reading gsp_config.dat.',&
         /'A common error is the wrong column number of gsp_config.csv.')
               ENDIF
            ENDDO ! END of lcurpar
            ! particle bounding box
            pbb(MM,:) = gsp_max(:) - gsp_min(:)
         ENDDO
         ! set tmp_gp_neighsid relative to the first component spheres for auto seeding
         IF(GSP_EXPLICIT) where (tmp_gp_neighsid /= -1) tmp_gp_neighsid = tmp_gp_neighsid - 1

      ENDIF
      IF(GSP_EXPLICIT) THEN
         CALL BCAST(start_index_sphere,pe_io)
         CALL BCAST(sphereRelPos,pe_io)
         CALL BCAST(sphereRelDia,pe_io)
         CALL BCAST(tmp_gp_sa, pe_io)
         CALL BCAST(tmp_gp_neighsid, pe_io)
         CALL BCAST(tmp_gp_neighsa, pe_io)
         CALL BCAST(gsp_vol, pe_io)
         CALL BCAST(pbb, pe_io)
      ELSE
         CALL BCAST(cspart_data, pe_io)
         CALL BCAST(gsp_vol, pe_io)
         CALL BCAST(pbb, pe_io)
      ENDIF

      CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) CALL LOG_ERROR()

      IF(myPE == PE_IO) THEN
         CLOSE(lUNIT1)
      ENDIF

      RETURN
   END SUBROUTINE READ_SPHERE_INPUT_GSP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: GSP_CONFIG_DAT_TO_CSV                                    !
!                                                                      !
! Author: Renjie Ke                                      Date: Jan 2025!
! Purpose: Migration function for users who still use gsp_config.dat   !
!          This subroutine shall be removed after some time            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE GSP_CONFIG_DAT_TO_CSV
      use geometry, only: NO_K
      USE physprop, only: d_p0
      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------

! local filename
      character(255) lfilename_dat
      character(255) lfilename_csv
! Flag to indicate if file exists.
      LOGICAL :: lEXISTS_DAT
      LOGICAL :: lEXISTS_CSV
! Read dimension: 2D vs 3D data
      integer :: RDIMN
! Buffer
      character(32) Buffer
! local unit
      integer, parameter :: lunit=10
      integer, parameter :: lunit1=11
! local format string
      character(255) FORMAT_STRING
! Line
      character(200) Line
      integer :: ppos,linepos
! Integer buffer
      integer :: int_buffer, DIM_NUM
! Number of solid phases using glued-sphere model. All others are assumed sphere.
      integer, dimension(:), allocatable :: idGluedPhase ! EX idGluedPhase = /1,3,5/
! Number of spheres in each non-spherical particle of each phase
      integer, dimension(:), allocatable :: sphereNumber_tmp
! Number of components per phase
      integer :: NM, TOTALNC, lcurpar, NC, PID, i, ilines, cur_pos, nbidx, neighbor_id
      double precision :: dist(3)
! Temp array to store the data
      double precision, dimension(:,:), allocatable :: darray1, darray2
      integer, dimension(:,:), allocatable :: iarray

! IO Status:
      integer :: IOS

      IOS = 0
      RDIMN = merge(2,3,NO_K)
! Check the the file exists and open it.
      LFILENAME_DAT = "gsp_config.dat"
      LFILENAME_CSV = "gsp_config.csv"

      IF(myPE == PE_IO) OPEN(UNIT=lUNIT, FILE=LFILENAME_DAT, FORM="FORMATTED")


 1101 FORMAT('FATAL - both gsp_config.dat and gsp_config.csv are not found.')
 1102 FORMAT('FATAL - gsp_config.dat not found, stop transferring to csv.')

! Collect the error message and quit.
      CALL GLOBAL_ALL_SUM(IOS)
      IF(IOS /= 0) call LOG_ERROR()

      IF (myPE .eq. PE_IO) THEN
! Skip Instructions
         ilines = 14
         do i = 1, ilines
            read(lunit,*) buffer
         enddo
! Dimension
         read(lunit,*) buffer,int_buffer
         DIM_NUM = int_buffer
! ID of Solid Phases using glued-sphere model
         read(lunit,'(a)') line
         linepos = index(line,':')
         line = line(linepos+1:)
         ! Number of solid phases given in the input file
         i = count(transfer(trim(adjustl(line)), 'a', len(trim(adjustl(line))))==" ")
         ALLOCATE(idGluedPhase(i+1))
         ALLOCATE(sphereNumber_tmp(i+1))
         read(line, *) idGluedPhase
! Number of Spheres in One Non-Spherical Particle in EachPhase
         read(lunit,*) buffer,sphereNumber_tmp
! Skip header
         do i = 1, 3
            read(lunit,*) buffer
         enddo
      ENDIF ! end of mype == pe_io

      IF (myPE .eq. PE_IO) THEN
         TOTALNC = 0
         DO NM = 1, SIZE(sphereNumber_tmp)
            TOTALNC = TOTALNC + sphereNumber_tmp(NM)
         ENDDO

         ALLOCATE(darray1(TOTALNC,5))
         ALLOCATE(iarray(TOTALNC,6))
         ALLOCATE(darray2(TOTALNC,6))
         DO lcurpar = 1, TOTALNC
            read (lunit,*,IOSTAT=IOS) (darray1(lcurpar, ppos), ppos=1,5), &
                        (iarray(lcurpar,ppos),ppos=1,6), &
                        (darray2(lcurpar,ppos),ppos=1,6)
         ENDDO
         CLOSE(lUNIT)

         OPEN(UNIT=lunit1, FILE=lfilename_csv, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND', IOSTAT=ios)
         WRITE(lunit1, '(A)') '#version = 1'
         ! write the header line in gsp_config.csv
         WRITE(lunit1, '(A)') 'x,y,z,rel_d,rel_sa,w,e,s,n,t,b,w_na,e_na,s_na,n_na,t_na,b_na,pid'

         FORMAT_STRING = '(G0,",",G0,",",G0,",",G0,",",G0,",",I0,",",I0,",",I0,",",I0,",",I0,",",' // &
                        'I0,",",G0,",",G0,",",G0,",",G0,",",G0,",",G0,",",I0)'

         ! In gsp_config.dat, the Wa, Ea, Sa, Na, Ta, Ba columns are neighbor surface area dividied by distance
         ! Scale back to area only
         cur_pos = 1
         DO NM = 1,SIZE(sphereNumber_tmp)
            NC = sphereNumber_tmp(NM)
            PID = idGluedPhase(NM)
            DO lcurpar = cur_pos, cur_pos + NC - 1
               DO nbidx = 1,6
                  IF(iarray(lcurpar,nbidx) == -1) CYCLE
                  neighbor_id = iarray(lcurpar,nbidx)
                  ! x,y,z are scaled by bounding diameter in gsp_config.dat, scale them back
                  dist(:) = abs((darray1(lcurpar,1:3) - darray1(cur_pos -1 + neighbor_id,1:3)) * d_p0(pid))
                  ! neighbor surface area is scaled by bounding diameter as well, scale them back and after calculation, divided by bounding diameter ** 2.0
                  IF(nbidx == 1 .or. nbidx == 2) THEN
                     ! darray2(lcurpar,nbidx) = darray2(lcurpar,nbidx) * d_p0(pid) * dist(1) / d_p0(pid) / d_p0(pid)
                     darray2(lcurpar,nbidx) = darray2(lcurpar,nbidx) * dist(1) / d_p0(pid)
                  ELSEIF(nbidx == 3 .or. nbidx == 4) THEN
                     darray2(lcurpar,nbidx) = darray2(lcurpar,nbidx) * dist(2) / d_p0(pid)
                  ELSE
                     darray2(lcurpar,nbidx) = darray2(lcurpar,nbidx) * dist(3) / d_p0(pid)
                  ENDIF
               ENDDO
               ! write the data line in gsp_config.csv
               WRITE(lunit1, FORMAT_STRING) (darray1(lcurpar, ppos), ppos=1,5), &
                                             (iarray(lcurpar,ppos),ppos=1,6), &
                                             (darray2(lcurpar,ppos),ppos=1,6), PID
            ENDDO
            cur_pos = cur_pos + NC
         ENDDO

         CLOSE(LUNIT1)

         WRITE(ERR_MSG, 99)
         CALL LOG_INFO()
 99 FORMAT('Info: Migration from gsp_config.dat to CSV is done')
      ENDIF ! end of mype == pe_io

      RETURN
   END SUBROUTINE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: WRITE_PART_OUTPUT                                        !
!                                                                      !
! Purpose: Write a particle output file at the end of a simulation.    !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE WRITE_PART_OUTPUT
!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      use cdist
      use compar
      use des_allocate, only: particle_grow
      use desmpi
      use desmpi_wrapper, only: DES_MPI_STOP
      use des_thermo, only: des_t_s
      use discretelement, only: des_usr_var_size
      use functions
      use funits
      use geometry, only: NO_K
      use mpi_comm_des, only: desmpi_scatterv
      use mpi_comm_des, only: desmpi_gatherv, des_gather
      use mpi_init_des, only: des_scatter_particle
      use mpi_utility
      use parallel_mpi
      use param, only: dimension_n_s
      use run, only: any_solids_species_eq, energy_eq
      use des_rxns, only: des_x_s

      implicit none

       CALL WRITE_PART_OUTPUT_V4P0
       return
   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: WRITE_PART_OUTPUT_V4P0                                   !
!                                                                      !
! Author: Renjie Ke                                   Date: 11-04-2024 !
! Revise: #1                                          Date: 11-08-2024 !
! Purpose:                                                             !
!   write a particle_output.csv.                                       !
!   also, when using gsp model and energy equation is solved,          !
!   write a intra_distribution.csv.                                    !
! Revise #1 content:                                                   !
!   add support for gsp model without intra heat transfer              !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE WRITE_PART_OUTPUT_V4P0

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      USE cdist
      USE compar
      USE des_allocate, only: particle_grow
      USE desmpi
      USE desmpi_wrapper, only: DES_MPI_STOP
      USE des_thermo, only: des_t_s
      USE discretelement, only: des_usr_var_size
      USE functions
      USE funits
      USE geometry, only: NO_K
      USE mpi_comm_des, only: desmpi_scatterv
      USE mpi_comm_des, only: desmpi_gatherv, des_gather
      USE mpi_init_des, only: des_scatter_particle
      USE mpi_utility
      USE parallel_mpi
      USE param, only: dimension_n_s
      USE run, only: any_solids_species_eq, energy_eq
      USE des_rxns, only: des_x_s

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: n
! Local unit
      INTEGER, PARAMETER :: lunit=10
      INTEGER, PARAMETER :: lunit_intra=11
! Local filename
      CHARACTER(255) lfilename
! IO Status:
      INTEGER :: IOS
! Read dimension: 2D vs 3D data
      INTEGER :: RDMN
! File version
      character(11) Version_string
! Global Particle count before filter is applied
      INTEGER :: NP
! Local and global Particle count after filter is applied
      INTEGER :: local_cnt, global_cnt
! Length of list of variables to read (number of columns)
      INTEGER :: n_vars
      INTEGER :: Phase_ID_start
      INTEGER :: Diameter_start
      INTEGER :: Density_start
      INTEGER :: Velocity_start
      INTEGER :: Temperature_start
      INTEGER :: Species_start
      INTEGER :: User_scalar_start
      INTEGER :: SQP_semi_axis_start
      INTEGER :: SQP_roundness_start
      INTEGER :: SQP_quaternion_start
      INTEGER :: GSP_quaternion_start
      INTEGER :: CGP_stat_wt_start
! Rank communication
      INTEGER :: lproc
! Header line
      CHARACTER(9999) :: column_header
! Variables related to gather
      INTEGER :: lgathercnts(0:numpes-1)
      INTEGER :: PC, LC1, LC2
      DOUBLE PRECISION, ALLOCATABLE :: ltemp_array(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: gtemp_array(:,:)
! Particle filter
      LOGICAL, allocatable, dimension(:) :: keep_particle
! Intra distribution
      INTEGER :: intra_n_vars, intra_temperature_start, intra_species_start, intra_User_scalar_start
      LOGICAL :: is_intra_distribution
      CHARACTER(9999) :: intra_column_header
      CHARACTER(255) :: lFILENAME_INTRA
      DOUBLE PRECISION, ALLOCATABLE :: intra_gtemp_array(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: intra_ltemp_array(:,:)
      CHARACTER(32) :: label
! Extra info to pass
      DOUBLE PRECISION, allocatable, dimension(:) :: gsp_temperature
      DOUBLE PRECISION, allocatable, dimension(:,:) :: gsp_species
      DOUBLE PRECISION, allocatable, dimension(:,:) :: gsp_user_scalar
      DOUBLE PRECISION, allocatable, dimension(:) :: gsp_phase
      INTEGER :: NS, gid, pid

      IF(.NOT.WRITE_PART_OUT)  RETURN

      WRITE(ERR_MSG, 10)
      CALL LOG_INFO()
      IF(GSP_EXPLICIT) THEN
         WRITE(ERR_MSG, 11)
         CALL LOG_INFO()
      ENDIF
 10   FORMAT('Writing DEM particle_output.csv file, version: 4.0.', &
      /'Please rename the file as "particle_input.csv" to use for', &
      /'another simulation initial condition.')

 11   FORMAT('GSP model is used: ', &
      /'an extra file "intra_distribution_output.csv" is generated,', &
      /'when using energy_equation/species/user_scalar.', &
      /'Rename this file to "intra_distribution.csv"', &
      /'and use the same "gsp_config.csv" and "GSP-(phase_id)" folders for another simulation.')

      n_vars = 0

!-----------------------------------------------

      IOS = 0
      RDMN = merge(2,3,NO_K)
      Version_string = 'Version 4.0'

! Initial the global count.
      GLOBAL_CNT = 10
! Calculate the number of 'real' particles on the local process.
      ! LOCAL_CNT = PIP - iGHOST_CNT
      LOCAL_CNT = 0

! Filter particles
      call echo_all_filter_options( &
      PART_OUT_X_MIN,PART_OUT_X_MAX,PART_OUT_X_EXCLUDE, &
      PART_OUT_Y_MIN,PART_OUT_Y_MAX,PART_OUT_Y_EXCLUDE, &
      PART_OUT_Z_MIN,PART_OUT_Z_MAX, PART_OUT_Z_EXCLUDE, &
      PART_OUT_PHASE, &
      PART_OUT_DIAMETER_MIN, PART_OUT_DIAMETER_MAX,PART_OUT_DIAMETER_EXCLUDE, &
      PART_OUT_DENSITY_MIN, PART_OUT_DENSITY_MAX,PART_OUT_DENSITY_EXCLUDE, &
      PART_OUT_U_MIN, PART_OUT_U_MAX,PART_OUT_U_EXCLUDE, &
      PART_OUT_V_MIN, PART_OUT_V_MAX,PART_OUT_V_EXCLUDE, &
      PART_OUT_W_MIN, PART_OUT_W_MAX,PART_OUT_W_EXCLUDE, &
      PART_OUT_TEMP_MIN, PART_OUT_TEMP_MAX,PART_OUT_TEMP_EXCLUDE, &
      PART_OUT_X_S_MIN,PART_OUT_X_S_MAX,PART_OUT_X_S_EXCLUDE, &
      PART_OUT_USR_VAR_MIN,PART_OUT_USR_VAR_MAX,PART_OUT_USR_VAR_EXCLUDE)

      if(allocated(keep_particle)) deallocate (keep_particle)
      allocate (keep_particle(MAX_PIP))
      keep_particle(:) = .TRUE.

      DO LC1 = 1, MAX_PIP
         IF(.NOT.IS_NORMAL(LC1)) keep_particle(LC1) = .FALSE.
      ENDDO

      call filter_particle_output(keep_particle)

      DO LC1 = 1, MAX_PIP
         IF(keep_particle(LC1)) LOCAL_CNT = LOCAL_CNT + 1
      ENDDO

! Calculate the total number of particles system-wide.
      NP = PIP - IGHOST_CNT
      CALL GLOBAL_ALL_SUM(NP)

      GLOBAL_CNT=10
      call global_all_sum(LOCAL_CNT, GLOBAL_CNT)

      WRITE(ERR_MSG, 20) NP,GLOBAL_CNT
      CALL LOG_INFO()

 20   FORMAT('Number of particles in the system                 = ', I9 ,&
            /'Number of particles left after filtering the data = ', I9)
! Set the send count from the local process.
      igath_sendcnt = LOCAL_CNT

! Collect the number of particles on each rank.all ranks.
      lgathercnts = 0
      lgathercnts(myPE) = LOCAL_CNT
      call global_sum(lgathercnts,igathercnts)

! Calculate the rank displacements.
      idispls(0) = 0
      DO lPROC = 1,NUMPEs-1
         idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)
      ENDDO

! Column header
      column_header = "X,Y,phase_id,diameter,density,U,V"
      n_vars = 7
      intra_n_vars = 0
      Phase_ID_start = 3
      Diameter_start = 4
      Density_start = 5
      Velocity_start = 6

      IF(RDMN==3) THEN
         column_header = "X,Y,Z,phase_id,diameter,density,U,V,W"
         n_vars = 9
         Phase_ID_start = 4
         Diameter_start = 5
         Density_start = 6
         Velocity_start = 7
      ENDIF

      is_intra_distribution = .False.

      IF(GSP_EXPLICIT) THEN
         is_intra_distribution = energy_eq .or. any_solids_species_eq .or. (DES_USR_VAR_SIZE>0)
         intra_column_header = "ID,Density"
         intra_n_vars = 2
      ENDIF

! Temperature option
      IF(ENERGY_EQ) THEN
         column_header = trim(column_header)//",t_s"
         Temperature_start   = n_vars + 1
         n_vars              = n_vars + 1
         intra_temperature_start = -1
         IF(is_intra_distribution) THEN
            intra_column_header = trim(intra_column_header)//",t_s"
            intra_temperature_start = intra_n_vars + 1
            intra_n_vars = intra_temperature_start
         ENDIF
      ELSE
         Temperature_start   = -1
         intra_temperature_start = -1
      ENDIF

! Species option
      IF(ANY_SOLIDS_SPECIES_EQ) THEN
         DO n = 1,DIMENSION_N_S
            write(label,'(",species", I0)') n
            column_header = trim(column_header)//trim(label)
            IF(is_intra_distribution) THEN
               intra_column_header = trim(intra_column_header)//trim(label)
            ENDIF
         ENDDO
         Species_start   = n_vars + 1
         n_vars          = n_vars + DIMENSION_N_S
         IF(is_intra_distribution) THEN
            intra_species_start = intra_n_vars + 1
            intra_n_vars = intra_n_vars + DIMENSION_N_S
         ENDIF
      ELSE
         Species_start   = -1
         intra_species_start = -1
      ENDIF

! User Scalar option
      IF(DES_USR_VAR_SIZE>0) THEN
         DO n = 1,DES_USR_VAR_SIZE
            write(label,'(",usr_var", I0)') n
            column_header = trim(column_header)//trim(label)
            IF(is_intra_distribution) THEN
               intra_column_header = trim(intra_column_header)//trim(label)
            ENDIF
         ENDDO
         User_scalar_start   = n_vars + 1
         n_vars              = n_vars + DES_USR_VAR_SIZE
         IF(is_intra_distribution) THEN
            intra_User_scalar_start = intra_n_vars + 1
            intra_n_vars = intra_n_vars + DES_USR_VAR_SIZE
         ENDIF
      ELSE
         User_scalar_start   = -1
         intra_User_scalar_start = -1
      ENDIF

! CGDEM specific
      IF(CGDEM) THEN
         column_header = trim(column_header)//",cgp_stat_wt"
         CGP_stat_wt_start = n_vars + 1
         n_vars            = n_vars + 1
      ELSE
         CGP_stat_wt_start = -1
      ENDIF

! GluedSphereDEM specific
      IF(GluedSphereDEM) THEN
         column_header = trim(column_header)//",gsp_q1,gsp_q2,gsp_q3,gsp_q4"
         GSP_quaternion_start = n_vars + 1
         n_vars               = n_vars + 4
      ELSE
         GSP_quaternion_start = -1
      ENDIF

      IF(SuperDEM) THEN
! SQP semi-axis option
         column_header = trim(column_header)//",sqp_a,sqp_b,sqp_c"
         SQP_semi_axis_start = n_vars + 1
         n_vars              = n_vars + 3

! SQP roundness option
         column_header = trim(column_header)//",sqp_m,sqp_n"
         SQP_roundness_start = n_vars + 1
         n_vars              = n_vars + 2

! SQP quaternion option
         column_header = trim(column_header)//",sqp_q1,sqp_q2,sqp_q3,sqp_q4"
         SQP_quaternion_start = n_vars + 1
         n_vars               = n_vars + 4
      ELSE
         SQP_semi_axis_start = -1
         SQP_roundness_start = -1
         SQP_quaternion_start = -1
      ENDIF

! Setup the file name
      lFILENAME= "particle_output.csv"
      ! exclusive for gsp model
      IF(is_intra_distribution) lFILENAME_INTRA = "intra_distribution_output.csv"

! Write header of into the csv file
      IF(myPE == PE_IO) THEN
         OPEN(UNIT=lUNIT, FILE=lFILENAME, FORM="FORMATTED")
         write(lUNIT,'(A)') TRIM(column_header)

         IF(is_intra_distribution) THEN
            OPEN(UNIT=lunit_intra, FILE=lFILENAME_INTRA, FORM="FORMATTED")
            write(lUNIT_intra,'(A)') TRIM(intra_column_header)
         ENDIF
      ENDIF

      ALLOCATE (dprocBuf(LOCAL_CNT))

      IF(is_intra_distribution) THEN
         ALLOCATE (intra_ltemp_array(intra_n_vars,LOCAL_CNT))
      ELSE
         ALLOCATE (ltemp_array(n_vars,LOCAL_CNT))
      ENDIF

      IF(myPE == PE_IO) THEN
         ALLOCATE (dRootBuf(GLOBAL_CNT))
         IF(is_intra_distribution) THEN
            ALLOCATE (intra_gtemp_array(intra_n_vars,GLOBAL_CNT))
         ELSE
            ALLOCATE (gtemp_array(n_vars,GLOBAL_CNT))
         ENDIF
      ELSE
         ALLOCATE (dRootBuf(10))
         IF(is_intra_distribution) THEN
            ALLOCATE (intra_gtemp_array(intra_n_vars,10))
         ELSE
            ALLOCATE (gtemp_array(n_vars,10))
         ENDIF
      ENDIF

! @Renjieke, MAX_PIP works fine for all other model except GSP
! For each vector component, pack component list in a local array

      IF(GSP_EXPLICIT) THEN
         IF(Temperature_Start>0) allocate(gsp_temperature(NGluedParticles))
         IF(Species_Start>0) allocate(gsp_species(NGluedParticles, DIMENSION_N_S))
         IF(User_scalar_start>0) allocate(gsp_user_scalar(NGluedParticles, DES_USR_VAR_SIZE))
         allocate(gsp_phase(NGluedParticles))

         IF(Temperature_Start>0) gsp_temperature = 0.0
         IF(Species_Start>0) gsp_species = 0.0
         IF(User_scalar_start>0) gsp_user_scalar = 0.0
         gsp_phase = 0.0
         DO LC1 = 1, MAX_PIP
            IF(PARTICLE_STATE(LC1) .ne. NORMAL_PARTICLE) CYCLE
            gid = gid_list(LC1)
            pid = PIJK(LC1,5)
            ns = sphereNumber(pid)
            ! @renjieke, we cannot simply make gsp_phase(gid) = pid, because for a gsp
            ! it's component spheres can be half normal and half ghost in a rank
            ! so it is possible for a specific gid, gsp_phase(gid) in two different ranks can have the same non-zero values
            ! so it will cause trouble if GLOBAL_ALL_SUM(gsp_phase) is called
            ! thus I use this DBLE(pid) / ns
            gsp_phase(gid) = gsp_phase(gid) + DBLE(pid) / DBLE(ns)
            IF(Temperature_Start>0) gsp_temperature(gid) = gsp_temperature(gid) + des_t_s(LC1) / DBLE(ns)
            IF(Species_Start>0) THEN
               DO LC2 = 1, DIMENSION_N_S
                  gsp_species(gid,LC2) = gsp_species(gid,LC2) + des_x_s(LC1,LC2) / DBLE(ns)
               ENDDO
            ENDIF
            IF(User_scalar_start>0) THEN
               DO LC2 = 1, DES_USR_VAR_SIZE
                  gsp_user_scalar(gid,LC2) = gsp_user_scalar(gid,LC2) + des_usr_var(LC2,LC1) / DBLE(ns)
               ENDDO
            ENDIF
         ENDDO
         CALL GLOBAL_ALL_SUM(gsp_phase)
         IF(Temperature_Start>0) call GLOBAL_ALL_SUM(gsp_temperature)
         IF(Species_Start>0) call GLOBAL_ALL_SUM(gsp_species)
         IF(User_scalar_start>0) call GLOBAL_ALL_SUM(gsp_user_scalar)
      ELSE
         PC = 0
         DO LC1 = 1, MAX_PIP
            IF(.NOT.KEEP_PARTICLE(LC1)) CYCLE
            PC =PC + 1
! Coordinates (x,y) in 2D, (x,y,z) in 3D
            DO LC2=1,RDMN
               ltemp_array(LC2,PC) = DES_POS_NEW(LC1,LC2)
            ENDDO
! Phase ID
            IF(Phase_ID_start>0) THEN
               ltemp_array(Phase_ID_start,PC) = dfloat(PIJK(LC1,5))
            ENDIF

! Diameter
            IF(Diameter_start>0) THEN
               ltemp_array(Diameter_start,PC) = 2.0D0*DES_RADIUS(LC1)
            ENDIF

! Density
            IF(Density_start>0) THEN
               ltemp_array(Density_start,PC) = ro_sol(LC1)
            ENDIF

! Velocity
            IF(Velocity_start>0) THEN
               DO LC2=1,RDMN
                  ltemp_array(Velocity_start+LC2-1,PC) = DES_VEL_NEW(LC1,LC2)
               ENDDO
            ENDIF
! Temperature
            IF(Temperature_start>0) THEN
               ltemp_array(Temperature_start,PC) = des_t_s(LC1)
            ENDIF

! Species
            IF(Species_start>0) THEN
               DO LC2=1,DIMENSION_N_S
                  ltemp_array(Species_start+LC2-1,PC) = des_x_s(LC1,LC2)
               ENDDO
            ENDIF

! Particle Scalar: CAUTION THIS ARRAY IS TRANSPOSED, FIRST INDEX IS USR_VAR ID, SECOND INDEX IS PARTICLE ID
            IF(User_scalar_start>0) THEN
               DO LC2=1,DES_USR_VAR_SIZE
                  ltemp_array(User_scalar_start+LC2-1,PC) = des_usr_var(LC2,LC1)
               ENDDO
            ENDIF

! SQP_semi-axes:
            IF(SQP_semi_axis_start>0) THEN
               DO LC2=1,3
                  ltemp_array(SQP_semi_axis_start+LC2-1,PC) = super_r(LC1,LC2)
               ENDDO
            ENDIF

! SQP_roundness:
            IF(SQP_roundness_start>0) THEN
               DO LC2=1,2
                  ltemp_array(SQP_roundness_start+LC2-1,PC) = super_mn(LC1,LC2)
               ENDDO
            ENDIF

! SQP_quaternion:
            IF(SQP_quaternion_start>0) THEN
               DO LC2=1,4
                  ltemp_array(SQP_quaternion_start+LC2-1,PC) = super_q(LC1,LC2)
               ENDDO
            ENDIF

! CGP_stat_wt:
            IF(CGP_stat_wt_start>0) THEN
               ltemp_array(CGP_stat_wt_start,PC) = des_cgp_stw(LC1)
            ENDIF

! GSP_quaternion:
            IF(GSP_quaternion_start>0) THEN
               DO LC2=1,4
                  ltemp_array(GSP_quaternion_start+LC2-1,PC) = glueQuat(LC1,LC2)
               ENDDO
            ENDIF

         ENDDO ! particle loop

! For each component, gather the local list to global temporary array
         DO LC1 = 1,n_vars
               dprocbuf(1:LOCAL_CNT)=ltemp_array(LC1,1:LOCAL_CNT)
               CALL desmpi_gatherv(ptype=2)
               gtemp_array(LC1,:) = drootbuf(:)
         ENDDO
      ENDIF ! end of non GSP and GSP without intra heat transfer

      IF(is_intra_distribution) THEN
         PC = 0
         DO LC1 = 1, MAX_PIP
            IF(.NOT.KEEP_PARTICLE(LC1)) CYCLE
            ! IF(particle_state .ne. NORMAL_PARTICLE) cycle
            PC =PC + 1

            intra_ltemp_array(1,PC) = iglobal_id(LC1)
! intra density
            intra_ltemp_array(2,PC) = ro_sol(LC1)
! intra temperature
            IF(intra_temperature_start>0) THEN
               intra_ltemp_array(intra_temperature_start,PC) = des_t_s(LC1)
            ENDIF
! intra Species
            IF(intra_species_start>0) THEN
               DO LC2=1,DIMENSION_N_S
                  intra_ltemp_array(intra_species_start+LC2-1,PC) = des_x_s(LC1,LC2)
               ENDDO
            ENDIF
! Intra Particle Scalar: CAUTION THIS ARRAY IS TRANSPOSED, FIRST INDEX IS USR_VAR ID, SECOND INDEX IS PARTICLE ID
            IF(intra_User_scalar_start>0) THEN
               DO LC2=1,DES_USR_VAR_SIZE
                  intra_ltemp_array(intra_User_scalar_start+LC2-1,PC) = des_usr_var(LC2,LC1)
               ENDDO
            ENDIF
         ENDDO

         DO LC1 = 1,intra_n_vars
            dprocbuf(1:LOCAL_CNT)=intra_ltemp_array(LC1,1:LOCAL_CNT)
            CALL desmpi_gatherv(ptype=2)
            intra_gtemp_array(LC1,:) = drootbuf(:)
         ENDDO
      ENDIF

! Write the data
      IF(mype == pe_io) THEN
! If below arrays are not allocated before, meaning they will not be used
! allocate to a small scale in case of memory error
         if(.not. allocated(gtemp_array)) allocate(gtemp_array(1,1))
         if(.not. allocated(intra_gtemp_array)) allocate(intra_gtemp_array(1,1))
         if(.not. allocated(gsp_temperature)) allocate(gsp_temperature(1))
         if(.not. allocated(gsp_species)) allocate(gsp_species(1,1))
         if(.not. allocated(gsp_user_scalar)) allocate(gsp_user_scalar(1,1))

         call write_part_data2csv(gtemp_array, intra_gtemp_array, &
            Temperature_start, Species_start, User_scalar_start, &
            SQP_semi_axis_start, SQP_roundness_start, SQP_quaternion_start, &
            intra_Temperature_start, intra_species_start, intra_User_scalar_start, &
            GSP_quaternion_start, &
            CGP_stat_wt_start, &
            Phase_ID_start,Diameter_start,Density_start,Velocity_start, &
            is_intra_distribution, &
            RDMN, GLOBAL_CNT, &
            lunit,lunit_intra, &
            gsp_temperature, gsp_phase, gsp_species, gsp_user_scalar, &
            KEEP_PARTICLE)
      ENDIF ! end of mype == pe_io

      if(allocated(dProcBuf)) deallocate(dProcBuf)
      if(allocated(dRootBuf)) deallocate(dRootBuf)
      if(allocated(ltemp_array)) deallocate(ltemp_array)
      if(allocated(gtemp_array)) deallocate(gtemp_array)
      if(allocated(intra_ltemp_array)) deallocate(intra_ltemp_array)
      if(allocated(intra_gtemp_array)) deallocate(intra_gtemp_array)

      if(allocated(gsp_temperature))deallocate(gsp_temperature)
      if(allocated(gsp_phase)) deallocate(gsp_phase)
      if(allocated(gsp_species)) deallocate(gsp_species)
      if(allocated(gsp_user_scalar)) deallocate(gsp_user_scalar)

      IF(myPE == PE_IO) CLOSE(lUNIT)
      IF(myPE == PE_IO) CLOSE(lunit_intra)

      RETURN

   END SUBROUTINE WRITE_PART_OUTPUT_V4P0

!  contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
! Subroutine: write_part_data2csv                                      !
!   dump particle data to string                                       !
!   and then dump string to csv files                                  !
!   why using string and character --> because comma                   !
!   there might be a better way to do that                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine write_part_data2csv(gtemp_array, intra_gtemp_array, &
               Temperature_start, Species_start, User_scalar_start, &
               SQP_semi_axis_start, SQP_roundness_start, SQP_quaternion_start, &
               intra_Temperature_start, intra_species_start, intra_User_scalar_start, &
               GSP_quaternion_start, &
               CGP_stat_wt_start, &
               Phase_ID_start,Diameter_start,Density_start,Velocity_start, &
               is_intra_distribution, &
               RDMN, GLOBAL_CNT, &
               lunit,lunit_intra, &
               gsp_temperature, gsp_phase, gsp_species, gsp_user_scalar, &
               KEEP_PARTICLE)
!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      USE param

      IMPLICIT NONE

      DOUBLE PRECISION, intent(in) :: gtemp_array(:,:)
      DOUBLE PRECISION, intent(in) :: intra_gtemp_array(:,:)
      INTEGER, intent(in) :: Temperature_start, Species_start, User_scalar_start
      INTEGER, intent(in) :: SQP_semi_axis_start, SQP_roundness_start, SQP_quaternion_start
      INTEGER, intent(in) :: intra_Temperature_start, intra_species_start, intra_User_scalar_start
      INTEGER, intent(in) :: Phase_ID_start, Diameter_start, Density_start, Velocity_start
      INTEGER, intent(in) :: GSP_quaternion_start
      INTEGER, intent(in) :: CGP_stat_wt_start
      INTEGER, intent(in) :: GLOBAL_CNT, RDMN
      INTEGER, intent(in) :: lunit, lunit_intra
      LOGICAL, intent(in) :: is_intra_distribution
      LOGICAL, intent(in) :: keep_particle(:)

      DOUBLE PRECISION, intent(in) :: gsp_temperature(:), gsp_phase(:)
      DOUBLE PRECISION, intent(in) :: gsp_species(:,:), gsp_user_scalar(:,:)

! local variables
      CHARACTER(255) :: temp_t_s_string, temp_pos_string, temp_dia_rho_vel_string, temp_phase_string
      CHARACTER(225) :: temp_gsp_quat_string
      CHARACTER(:), ALLOCATABLE :: temp_species_string, temp_user_scalar_string
      CHARACTER(:), ALLOCATABLE :: temp_row_string_gsp
      CHARACTER(:), ALLOCATABLE :: temp_row_string
      CHARACTER(255) :: temp_id_string
      CHARACTER(255) :: temp_sqp_abc_string
      CHARACTER(255) :: temp_sqp_mn_string
      CHARACTER(255) :: temp_sqp_quat_string
      CHARACTER(255) :: temp_rho_string
      CHARACTER(255) :: temp_cgp_stat_wt_string
      CHARACTER(32) :: label

      INTEGER :: LC1, LC2, r1, c1, global_id, i, j
      INTEGER :: row_string_size
      DOUBLE PRECISION :: dummy_zero(3)
      DOUBLE PRECISION, allocatable, dimension(:,:) :: intra_gtemp_array_sorted, temp_sorted
      LOGICAL :: rk_debug = .False.
! --------------------------------------------------------------------------------------------------->
      IF(rk_debug) THEN
         ! DEBUG CHECK
         print*,"debug check in part_data2string"
         print*,"gtemp_array size:",size(gtemp_array,1),size(gtemp_array,2)
         print*,"intra_gtemp_array size:",size(intra_gtemp_array,1),size(intra_gtemp_array,2)
         print*,"Temperature_start, Species_start, User_scalar_start:",Temperature_start, Species_start, User_scalar_start
         print*,"SQP_semi_axis_start, SQP_roundness_start, SQP_quaternion_start:",SQP_semi_axis_start, SQP_roundness_start, SQP_quaternion_start
         print*,"intra_Temperature_start, intra_species_start, intra_User_scalar_start:",intra_Temperature_start, intra_species_start, intra_User_scalar_start
         print*,"GSP_quaternion_start:",GSP_quaternion_start
         print*,"CGP_stat_wt_start:",CGP_stat_wt_start
         print*,"Phase_ID_start,Diameter_start,Density_start,Velocity_start:",Phase_ID_start,Diameter_start,Density_start,Velocity_start
         print*,"is SuperDEM?  ",SuperDEM
         print*,"is GluedSphereDEM?  ",GluedSphereDEM
         print*,"is is_intra_distribution?  ",is_intra_distribution
         print*,"GLOBAL_CNT,NGluedParticles:",GLOBAL_CNT,NGluedParticles
         print*,"lunit,lunit_intra:",lunit,lunit_intra
         print*,"DIMENSION_N_S, DES_USR_VAR_SIZE:",DIMENSION_N_S, DES_USR_VAR_SIZE
         print*,"keep_particle size: ",size(keep_particle,1)
         PRINT*,'============================================================================'
      ENDIF

      dummy_zero(:) = 0.0

      IF(gsp_explicit) THEN
         DO LC1 = 1, NGluedParticles
            if(particle_state_gsp(LC1) .lt. 1) CYCLE
            WRITE(label,*) glueBounding(LC1)
            temp_dia_rho_vel_string = ","//trim(adjustl(label))
            WRITE(label,*) (glueMass(LC1)/glueVolume(LC1))
            temp_dia_rho_vel_string = trim(temp_dia_rho_vel_string)//","//trim(adjustl(label))
            temp_pos_string = ""
            DO LC2 = 1,3
               WRITE(label, *) gluePos(LC1,LC2)
               temp_pos_string = trim(temp_pos_string)//","//trim(adjustl(label))
               IF(PART_OUT_ZERO_VEL) THEN
                  WRITE(label, '(F18.1)') dummy_zero(LC2)
               ELSE
                  WRITE(label, *) glueVel(LC1,LC2)
               ENDIF
               temp_dia_rho_vel_string = trim(temp_dia_rho_vel_string)//","//trim(adjustl(label))
            ENDDO
            temp_pos_string = temp_pos_string(2:)
            temp_phase_string = ","//trim(iVal(NINT(gsp_phase(LC1))))
! Take average of temperature, species, user_scalar on component spheres
! However, they will be overridden by intra_distribution.csv anyway
            IF(Temperature_start>0) THEN
               WRITE(temp_t_s_string,*) gsp_temperature(LC1)
               temp_t_s_string = ","//trim(adjustl(temp_t_s_string))
            ENDIF

            IF(Species_start>0) THEN
               IF(.NOT. ALLOCATED(temp_species_string)) ALLOCATE(CHARACTER(len=DIMENSION_N_S*19) :: temp_species_string)
               temp_species_string = ""
               DO LC2=1,DIMENSION_N_S
                  WRITE(label, *) gsp_species(LC1,LC2)
                  temp_species_string = trim(temp_species_string)//","//trim(adjustl(label))
               ENDDO
            ENDIF

            IF(User_scalar_start>0) THEN
               IF(.NOT. ALLOCATED(temp_user_scalar_string)) ALLOCATE(CHARACTER(len=DES_USR_VAR_SIZE*19) :: temp_user_scalar_string)
               temp_user_scalar_string = ""
               DO LC2=1,DES_USR_VAR_SIZE
                  WRITE(label, *) gsp_user_scalar(LC1,LC2)
                  temp_user_scalar_string = trim(temp_user_scalar_string)//","//trim(adjustl(label))
               ENDDO
            ENDIF

            temp_gsp_quat_string = ""
            IF(GSP_quaternion_start>0) THEN
               DO LC2=1,4
                  WRITE(label, *) glueQuat(LC1,LC2)
                  temp_gsp_quat_string = trim(temp_gsp_quat_string)//","//trim(adjustl(label))
               ENDDO
            ENDIF

            ! Calculate the size of the row string
            row_string_size = LEN_TRIM(temp_pos_string) +  LEN_TRIM(temp_phase_string) + LEN_TRIM(temp_dia_rho_vel_string)
            IF(Temperature_start>0) row_string_size = row_string_size + LEN_TRIM(temp_t_s_string)
            IF(Species_start>0) row_string_size = row_string_size + LEN_TRIM(temp_species_string)
            IF(User_scalar_start>0) row_string_size = row_string_size + LEN_TRIM(temp_user_scalar_string)
            IF(GSP_quaternion_start>0) row_string_size = row_string_size + LEN_TRIM(temp_gsp_quat_string)
            ! Allocate the array with a slightly larger size to ensure it can hold all the data.
            row_string_size = row_string_size + 10

            ! Dynamically allocate the size of the temp_row_string_gsp
            IF(.NOT. ALLOCATED(temp_row_string_gsp)) ALLOCATE(CHARACTER(len=row_string_size) :: temp_row_string_gsp)
            temp_row_string_gsp = ""
            temp_row_string_gsp = trim(temp_pos_string)//trim(temp_phase_string)//trim(temp_dia_rho_vel_string)
            IF(Temperature_start>0) temp_row_string_gsp = trim(temp_row_string_gsp) // trim(temp_t_s_string)
            IF(Species_start>0) temp_row_string_gsp = trim(temp_row_string_gsp) // trim(temp_species_string)
            IF(User_scalar_start>0) temp_row_string_gsp = trim(temp_row_string_gsp) // trim(temp_user_scalar_string)
            IF(GSP_quaternion_start>0) temp_row_string_gsp = trim(temp_row_string_gsp) // trim(temp_gsp_quat_string)
            WRITE(lUNIT,'(A)') trim(temp_row_string_gsp)
         ENDDO

         IF(is_intra_distribution) THEN
! Reorder the intra_gtemp_array to intra_gtemp_array_sorted
            r1 = size(intra_gtemp_array,1)
            c1 = size(intra_gtemp_array,2)
            allocate(intra_gtemp_array_sorted(r1,c1))
            allocate(temp_sorted(r1,particles)) ! XXX
            intra_gtemp_array_sorted(:,:) = 0.0D0
            temp_sorted(:,:) = 0.0D0
            DO LC1 = 1, GLOBAL_CNT
               global_id = intra_gtemp_array(1,LC1)
               temp_sorted(1, global_id) = global_id
               temp_sorted(2:, global_id) = intra_gtemp_array(2:, LC1)
            ENDDO

            j = 1
            DO i = 1, particles
               if(temp_sorted(1, i) .ne. 0.0) then
                  intra_gtemp_array_sorted(:,j) = temp_sorted(:, i)
                  j = j + 1
               endif
            ENDDO

            DO LC1 = 1, GLOBAL_CNT
               temp_id_string = ""
               temp_id_string = trim(iVal(INT(intra_gtemp_array_sorted(1,LC1))))

               temp_rho_string = ""
               WRITE(temp_rho_string,*) intra_gtemp_array_sorted(2,LC1)
               temp_rho_string = ","//trim(adjustl(temp_rho_string))

               IF(intra_temperature_start>0) THEN
                  temp_t_s_string = ""
                  WRITE(temp_t_s_string,*) intra_gtemp_array_sorted(intra_temperature_start,LC1)
                  temp_t_s_string = ","//trim(adjustl(temp_t_s_string))
               ENDIF

               IF(intra_species_start>0) THEN
                  IF(.NOT. ALLOCATED(temp_species_string)) ALLOCATE(CHARACTER(len=DIMENSION_N_S*19) :: temp_species_string)
                  temp_species_string = ""
                  DO LC2=1,DIMENSION_N_S
                     WRITE(label, *) intra_gtemp_array_sorted(intra_species_start+LC2-1,LC1)
                     temp_species_string = trim(temp_species_string)//","//trim(adjustl(label))
                  ENDDO
               ENDIF

               IF(intra_User_scalar_start > 0) THEN
                  IF(.NOT. ALLOCATED(temp_user_scalar_string)) ALLOCATE(CHARACTER(len=DES_USR_VAR_SIZE*19) :: temp_user_scalar_string)
                  temp_user_scalar_string = ""
                  DO LC2=1,DES_USR_VAR_SIZE
                     WRITE(label, *) intra_gtemp_array_sorted(intra_User_scalar_start+LC2-1,LC1)
                     temp_user_scalar_string = trim(temp_user_scalar_string)//","//trim(adjustl(label))
                  ENDDO
               ENDIF

               row_string_size = LEN_TRIM(temp_id_string) +  LEN_TRIM(temp_rho_string)
               IF(intra_Temperature_start>0) row_string_size = row_string_size + LEN_TRIM(temp_t_s_string)
               IF(intra_species_start>0) row_string_size = row_string_size + LEN_TRIM(temp_species_string)
               IF(intra_User_scalar_start>0) row_string_size = row_string_size + LEN_TRIM(temp_user_scalar_string)
               row_string_size = row_string_size + 10

               IF(.NOT. ALLOCATED(temp_row_string)) ALLOCATE(CHARACTER(len=row_string_size) :: temp_row_string)
               temp_row_string = ""
               temp_row_string = trim(temp_id_string)//trim(temp_rho_string)
               IF(intra_Temperature_start>0) temp_row_string = trim(temp_row_string) // trim(temp_t_s_string)
               IF(intra_species_start>0) temp_row_string = trim(temp_row_string) // trim(temp_species_string)
               IF(intra_User_scalar_start>0) temp_row_string = trim(temp_row_string) // trim(temp_user_scalar_string)
               WRITE(lUNIT_intra,'(A)') trim(temp_row_string)
            ENDDO
         ENDIF
      ELSE ! end of gsp for intra_distribution.csv
         DO LC1 = 1, GLOBAL_CNT
            WRITE(label,*) gtemp_array(Diameter_start,LC1)
            temp_dia_rho_vel_string = ","//trim(adjustl(label))
            WRITE(label,*) gtemp_array(Density_start,LC1)
            temp_dia_rho_vel_string = trim(temp_dia_rho_vel_string)//","//trim(adjustl(label))
            temp_pos_string = ""
            DO LC2 = 1,RDMN
               WRITE(label, *) gtemp_array(LC2,LC1)
               temp_pos_string = trim(temp_pos_string)//","//trim(adjustl(label))
               IF(PART_OUT_ZERO_VEL) THEN
                  WRITE(label, '(F18.1)') dummy_zero(LC2)
               ELSE
                  WRITE(label, *) gtemp_array(Velocity_start+LC2-1,LC1)
               ENDIF
               temp_dia_rho_vel_string = trim(temp_dia_rho_vel_string)//","//trim(adjustl(label))
            ENDDO
            temp_pos_string = temp_pos_string(2:)
            temp_phase_string = ","//trim(iVal(NINT(gtemp_array(Phase_ID_start,LC1))))

            IF(Temperature_start>0) THEN
               WRITE(temp_t_s_string,*) gtemp_array(Temperature_start,LC1)
               temp_t_s_string = ","//trim(adjustl(temp_t_s_string))
            ENDIF

            IF(Species_start>0) THEN
               IF(.NOT. ALLOCATED(temp_species_string)) ALLOCATE(CHARACTER(len=DIMENSION_N_S*19) :: temp_species_string)
               temp_species_string = ""
               DO LC2=1,DIMENSION_N_S
                  WRITE(label, *) gtemp_array(Species_start+LC2-1,LC1)
                  temp_species_string = trim(temp_species_string)//","//trim(adjustl(label))
               ENDDO
            ENDIF

            IF(User_scalar_start>0) THEN
               IF(.NOT. ALLOCATED(temp_user_scalar_string)) ALLOCATE(CHARACTER(len=DES_USR_VAR_SIZE*19) :: temp_user_scalar_string)
               temp_user_scalar_string = ""
               DO LC2=1,DES_USR_VAR_SIZE
                  WRITE(label, *) gtemp_array(User_scalar_start+LC2-1,LC1)
                  temp_user_scalar_string = trim(temp_user_scalar_string)//","//trim(adjustl(label))
               ENDDO
            ENDIF

            IF(SQP_semi_axis_start>0) THEN
               temp_sqp_abc_string = ""
               DO LC2=1,3
                  WRITE(label, *) gtemp_array(SQP_semi_axis_start+LC2-1,LC1)
                  temp_sqp_abc_string = trim(temp_sqp_abc_string)//","//trim(adjustl(label))
               ENDDO
            ENDIF

            IF(SQP_roundness_start>0) THEN
               temp_sqp_mn_string = ""
               DO LC2=1,2
                  WRITE(label, *) gtemp_array(SQP_roundness_start+LC2-1,LC1)
                  temp_sqp_mn_string = trim(temp_sqp_mn_string)//","//trim(adjustl(label))
               ENDDO
            ENDIF

            IF(SQP_quaternion_start>0) THEN
               temp_sqp_quat_string = ""
               DO LC2=1,4
                  WRITE(label, *) gtemp_array(SQP_quaternion_start+LC2-1,LC1)
                  temp_sqp_quat_string = trim(temp_sqp_quat_string)//","//trim(adjustl(label))
               ENDDO
            ENDIF

            IF(CGP_stat_wt_start>0) THEN
               temp_cgp_stat_wt_string = ","//trim(iVal(gtemp_array(CGP_stat_wt_start,LC1)))
            ENDIF

            IF(GSP_quaternion_start>0) THEN
            ! @renjieke, glueQuat is now attached to pip in GSP project2
            ! so they need to be put in gtemp_array as well
               temp_gsp_quat_string = ""
               DO LC2=1,4
                  WRITE(label, *) gtemp_array(GSP_quaternion_start+LC2-1,LC1)
                  temp_gsp_quat_string = trim(temp_gsp_quat_string)//","//trim(adjustl(label))
               ENDDO
            ENDIF

            row_string_size = LEN_TRIM(temp_pos_string) + LEN_TRIM(temp_phase_string) + LEN_TRIM(temp_dia_rho_vel_string)
            IF(Temperature_start>0) row_string_size = row_string_size + LEN_TRIM(temp_t_s_string)
            IF(Species_start>0) row_string_size = row_string_size + LEN_TRIM(temp_species_string)
            IF(User_scalar_start>0) row_string_size = row_string_size + LEN_TRIM(temp_user_scalar_string)
            IF(SQP_semi_axis_start>0) row_string_size = row_string_size + LEN_TRIM(temp_sqp_abc_string)
            IF(SQP_roundness_start>0) row_string_size = row_string_size + LEN_TRIM(temp_sqp_mn_string)
            IF(SQP_quaternion_start>0) row_string_size = row_string_size + LEN_TRIM(temp_sqp_quat_string)
            IF(CGP_stat_wt_start>0) row_string_size = row_string_size + LEN_TRIM(temp_cgp_stat_wt_string)
            IF(GSP_quaternion_start>0) row_string_size = row_string_size + LEN_TRIM(temp_gsp_quat_string)
            ! Allocate the array with a slightly larger size to ensure it can hold all the data.
            row_string_size = row_string_size + 10

            ! Dynamically allocate the size of the temp_row_string
            IF(.NOT. ALLOCATED(temp_row_string)) ALLOCATE(CHARACTER(len=row_string_size) :: temp_row_string)
            temp_row_string = ""

            temp_row_string = trim(temp_pos_string)//trim(temp_phase_string)//trim(temp_dia_rho_vel_string)
            IF(Temperature_start>0) temp_row_string = trim(temp_row_string) // trim(temp_t_s_string)
            IF(Species_start>0) temp_row_string = trim(temp_row_string) // trim(temp_species_string)
            IF(User_scalar_start>0) temp_row_string = trim(temp_row_string) // trim(temp_user_scalar_string)
            IF(SQP_semi_axis_start>0) temp_row_string = trim(temp_row_string) // trim(temp_sqp_abc_string)
            IF(SQP_roundness_start>0) temp_row_string = trim(temp_row_string) // trim(temp_sqp_mn_string)
            IF(SQP_quaternion_start>0) temp_row_string = trim(temp_row_string) // trim(temp_sqp_quat_string)
            IF(CGP_stat_wt_start>0) temp_row_string = trim(temp_row_string) // trim(temp_cgp_stat_wt_string)
            IF(GSP_quaternion_start>0) temp_row_string = trim(temp_row_string) // trim(temp_gsp_quat_string)
            WRITE(lUNIT,'(A)') trim(temp_row_string)
         ENDDO
         IF(allocated(intra_gtemp_array_sorted)) deallocate(intra_gtemp_array_sorted)
         IF(allocated(temp_sorted)) deallocate(temp_sorted)
      ENDIF ! end of non gsp
      IF(allocated(temp_row_string_gsp)) deallocate(temp_row_string_gsp)
      IF(allocated(temp_row_string)) deallocate(temp_row_string)
      IF(allocated(temp_species_string)) deallocate(temp_species_string)
      IF(allocated(temp_user_scalar_string)) deallocate(temp_user_scalar_string)
   end subroutine write_part_data2csv

      subroutine check_if_uniform(data,uniform, uniform_value)

           USE param1, only: undefined

           IMPLICIT NONE

           DOUBLE PRECISION, INTENT(in) :: DATA(:)
           LOGICAL, INTENT(OUT) :: uniform
           DOUBLE PRECISION, INTENT(OUT) :: uniform_value
           DOUBLE PRECISION :: lmin,lmax,gmin,gmax
           DOUBLE PRECISION, PARAMETER :: tol = 1.0D-6
           INTEGER :: NP

           lmin =  UNDEFINED
           lmax = -UNDEFINED
           DO NP=1,MAX_PIP
               IF(.NOT.IS_NORMAL(NP)) CYCLE
               lmin = MIN(lmin,DATA(NP))
               lmax = MAX(lmax,DATA(NP))
           ENDDO

           call global_all_min(lmin, gmin)
           call global_all_max(lmax, gmax)

           if(dabs(gmax-gmin)<tol) then
               uniform = .True.
               uniform_value = gmin
           else
               uniform = .false.
               uniform_value = undefined
           endif

           return

      end subroutine check_if_uniform


      subroutine filter_particle_output(keep_particle)

           USE param1, only: undefined

           IMPLICIT NONE

           LOGICAL, INTENT(inout) :: KEEP_PARTICLE(:)
           INTEGER :: n,NP,M
           character(32) :: label

           DO NP=1,MAX_PIP
               IF(.NOT.IS_NORMAL(NP)) THEN
                   KEEP_PARTICLE(NP)=.FALSE.
                   CYCLE
               ENDIF
               ! phase id
               M = PIJK(NP,5)
               IF(.NOT.part_out_phase(M)) KEEP_PARTICLE(NP)=.FALSE.
           ENDDO

! x, y and z-coordinate
           call filter_based_on_min_max_var(des_pos_new(:,1),'x-coordinate',part_out_x_min,part_out_x_max,part_out_x_exclude,keep_particle)
           call filter_based_on_min_max_var(des_pos_new(:,2),'y-coordinate',part_out_y_min,part_out_y_max,part_out_y_exclude,keep_particle)
           call filter_based_on_min_max_var(des_pos_new(:,3),'z-coordinate',part_out_z_min,part_out_z_max,part_out_z_exclude,keep_particle)

           call filter_based_on_min_max_var(2.0D0*des_radius,'Diameter',part_out_diameter_min,part_out_diameter_max,part_out_diameter_exclude,keep_particle)
           call filter_based_on_min_max_var(ro_sol,'Density',part_out_density_min,part_out_density_max,part_out_density_exclude,keep_particle)

           call filter_based_on_min_max_var(des_vel_new(:,1),'x-velocity',part_out_u_min,part_out_u_max,part_out_u_exclude,keep_particle)
           call filter_based_on_min_max_var(des_vel_new(:,2),'y-velocity',part_out_v_min,part_out_v_max,part_out_v_exclude,keep_particle)
           call filter_based_on_min_max_var(des_vel_new(:,3),'z-velocity',part_out_w_min,part_out_w_max,part_out_w_exclude,keep_particle)

           IF(ENERGY_EQ) THEN
               call filter_based_on_min_max_var(des_t_s,'Temperature',part_out_temp_min,part_out_temp_max,part_out_temp_exclude,keep_particle)
           ENDIF

           IF(ANY_SOLIDS_SPECIES_EQ) THEN
               do n = 1,DIMENSION_N_S
                   write(label,'("Species",I2)') n
                   call filter_based_on_min_max_var(des_x_s(:,n),label,part_out_x_s_min(n),part_out_x_s_max(n),part_out_x_s_exclude(n),keep_particle)
          enddo
      ENDIF

! Particle Scalar: CAUTION THIS ARRAY IS TRANSPOSED, FIRST INDEX IS USR_VAR ID, SECOND INDEX IS PARTICLE ID
      IF(DES_USR_VAR_SIZE>0) THEN
          do n = 1,DES_USR_VAR_SIZE
              write(label,'("User scalar # ",I2)') n
              call filter_based_on_min_max_var(des_usr_var(n,:),label,part_out_usr_var_min(n),part_out_usr_var_max(n),part_out_usr_var_exclude(n),keep_particle)
          enddo
      ENDIF
      return
 end subroutine filter_particle_output

END SUBROUTINE WRITE_PART_OUTPUT



subroutine filter_based_on_min_max_var(data,label,min_value,max_value,exclude_within_range,keep_particle)

      USE discretelement
      USE param1, only: undefined

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in) :: DATA(:)
      CHARACTER(LEN=*), INTENT(in) :: label
      DOUBLE PRECISION, INTENT(in) :: min_value,max_value
      LOGICAL, INTENT(in) :: exclude_within_range
      LOGICAL, INTENT(inout) :: KEEP_PARTICLE(:)
      INTEGER :: NP

! Skip filtering altogether if min and max values are undefined.
! This leaves KEEP_PARTICLE unchanged

      IF(min_value==-UNDEFINED.and.max_value==UNDEFINED) RETURN

!      WRITE(ERR_MSG, 100) trim(label)
!      CALL LOG_INFO()

! 100 FORMAT('Info: Filtering DEM particle_[in/out]put.dat file based on: ', A)

      IF(exclude_within_range) then
         DO NP=1,MAX_PIP
            IF(.NOT.KEEP_PARTICLE(NP)) CYCLE
            IF(min_value<=DATA(NP).AND.DATA(NP)<=max_value) KEEP_PARTICLE(NP)=.FALSE.
         ENDDO
      ELSE
         DO NP=1,MAX_PIP
            IF(.NOT.KEEP_PARTICLE(NP)) CYCLE
            IF(DATA(NP)<=min_value.OR.DATA(NP)>=max_value) KEEP_PARTICLE(NP)=.FALSE.
         ENDDO
      ENDIF

      return

 end subroutine filter_based_on_min_max_var

 subroutine filter_single_particle_based_on_min_max_var(data,label,min_value,max_value,exclude_within_range,keep_particle)

      USE discretelement
      USE param1, only: undefined

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in) :: DATA
      CHARACTER(LEN=*), INTENT(in) :: label
      DOUBLE PRECISION, INTENT(in) :: min_value,max_value
      LOGICAL, INTENT(in) :: exclude_within_range
      LOGICAL, INTENT(inout) :: KEEP_PARTICLE
      INTEGER :: NP

! Single particle version, used when reading particle_input.dat file
! Skip filtering altogether if min and max values are undefined.
! This leaves KEEP_PARTICLE unchanged

      IF(min_value==-UNDEFINED.and.max_value==UNDEFINED) RETURN
      IF(.NOT.KEEP_PARTICLE) RETURN

!      WRITE(ERR_MSG, 100) trim(label)
!      CALL LOG_INFO()

! 100 FORMAT('Info: Filtering DEM particle_[in/out]put.dat file based on: ', A)
! if one of min_value or max_value is defined, filter the data.
! If keep_within_range is .True., we keep particles if   min_value<=data<=max_value
! If keep_within_range is .False., we keep particles if  data<=min_value  and  data>=max_value
! Note how this is coded differently since we use the flag KEEP_PARTICLE.

      IF(exclude_within_range) then
          IF(min_value<=data .AND. data<=max_value) KEEP_PARTICLE=.FALSE.
      ELSE
          IF(data<=min_value .OR. data>=max_value) KEEP_PARTICLE=.FALSE.
      ENDIF

      return

 end subroutine filter_single_particle_based_on_min_max_var

 subroutine echo_filter_option(label,min_value,max_value,exclude_within_range,n_filters)

      USE param1, only: undefined

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(in) :: label
      DOUBLE PRECISION, INTENT(in) :: min_value,max_value
      LOGICAL, INTENT(in) :: exclude_within_range
      INTEGER, INTENT(inout) :: n_filters
      CHARACTER(7) :: filter_type
      ! Echo the filtering option for a given variable.
      ! If no filter is applied, nothing is displayed

      IF(min_value==-UNDEFINED.and.max_value==UNDEFINED) RETURN

      n_filters = n_filters + 1

      filter_type = merge('EXCLUDE','KEEP   ',exclude_within_range)
      WRITE(ERR_MSG, 101) trim(label), filter_type,min_value,max_value
      CALL LOG_MESSAGE(__FILE__,__LINE__,LOGLEVEL_INFO,header=.false.,footer=.false.)

 101 FORMAT(A12,': ',A7, ' data between ',G0, ' and ',G0)


      return

 end subroutine echo_filter_option

 subroutine check_particle_input_filter_with_V1

      USE discretelement
      USE param1, only: undefined
      use param, only: dimension_n_s
      use run, only: energy_eq

      IMPLICIT NONE

      INTEGER :: M,N
      LOGICAL :: PART_IN_FILTER

! If PART_IN_* are not the default values and we use version 1.0 of
! particle_input.dat, flag this as a fatal error.

      PART_IN_FILTER = .FALSE.
      IF(PART_IN_X_MIN/=-UNDEFINED.OR.PART_IN_X_MAX/=UNDEFINED.OR.PART_IN_X_EXCLUDE) PART_IN_FILTER = .TRUE.
      IF(PART_IN_Y_MIN/=-UNDEFINED.OR.PART_IN_Y_MAX/=UNDEFINED.OR.PART_IN_Y_EXCLUDE) PART_IN_FILTER = .TRUE.
      IF(PART_IN_Z_MIN/=-UNDEFINED.OR.PART_IN_Z_MAX/=UNDEFINED.OR.PART_IN_Z_EXCLUDE) PART_IN_FILTER = .TRUE.
      DO M = 1,DES_MMAX
         IF(.NOT.PART_IN_PHASE(M)) PART_IN_FILTER = .TRUE.
      END DO
      IF(PART_IN_DIAMETER_MIN/=-UNDEFINED.OR.PART_IN_DIAMETER_MAX/=UNDEFINED.OR.PART_IN_DIAMETER_EXCLUDE) PART_IN_FILTER = .TRUE.
      IF(PART_IN_DENSITY_MIN/=-UNDEFINED.OR.PART_IN_DENSITY_MAX/=UNDEFINED.OR.PART_IN_DENSITY_EXCLUDE) PART_IN_FILTER = .TRUE.
      IF(PART_IN_U_MIN/=-UNDEFINED.OR.PART_IN_U_MAX/=UNDEFINED.OR.PART_IN_U_EXCLUDE) PART_IN_FILTER = .TRUE.
      IF(PART_IN_V_MIN/=-UNDEFINED.OR.PART_IN_V_MAX/=UNDEFINED.OR.PART_IN_V_EXCLUDE) PART_IN_FILTER = .TRUE.
      IF(PART_IN_W_MIN/=-UNDEFINED.OR.PART_IN_W_MAX/=UNDEFINED.OR.PART_IN_W_EXCLUDE) PART_IN_FILTER = .TRUE.
      IF(ENERGY_EQ) THEN
         IF(PART_IN_TEMP_MIN/=-UNDEFINED.OR.PART_IN_TEMP_MAX/=UNDEFINED.OR.PART_IN_TEMP_EXCLUDE) PART_IN_FILTER = .TRUE.
      ENDIF
      DO N = 1,DIMENSION_N_S
         IF(PART_IN_X_S_MIN(N)/=-UNDEFINED.OR.PART_IN_X_S_MAX(N)/=UNDEFINED.OR.PART_IN_X_S_EXCLUDE(N)) PART_IN_FILTER = .TRUE.
      ENDDO
      DO N = 1,DES_USR_VAR_SIZE
         IF(PART_IN_USR_VAR_MIN(N)/=-UNDEFINED.OR.PART_IN_USR_VAR_MAX(N)/=UNDEFINED.OR.PART_IN_USR_VAR_EXCLUDE(N)) PART_IN_FILTER = .TRUE.
      ENDDO

      IF(PART_IN_FILTER.AND.P_INPUT_DAT_VERSION == '1.0') THEN
         WRITE(ERR_MSG, 100)
         CALL LOG_ERROR()
      ENDIF

100 FORMAT('Filtering particle_input.dat is not supported for version 1.0.', &
          /'Please use particle_input.dat version 2.0 or above', &
          /'or remove all filtering parameters to keep using version 1.0.')

      RETURN

 end subroutine check_particle_input_filter_with_V1


 subroutine check_number_of_columns(first_line, n_vars,&
     XYZ_start, Phase_ID_start,Diameter_start,Density_start,Velocity_start, &
     Temperature_start,Species_start,User_scalar_start)

      use discretelement, only: des_usr_var_size
      use geometry, only: NO_K
      use param, only: dimension_n_s
      use run, only: any_solids_species_eq, energy_eq

      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Read dimension: 2D vs 3D data
      integer :: RDIMN
!-----------------------------------------------
      ! character(1024), intent(in) :: first_line
      character(1024) :: first_line
      integer, intent(in) :: n_vars
      integer, intent(in) :: XYZ_start
      integer, intent(in) :: Phase_ID_start
      integer, intent(in) :: Diameter_start
      integer, intent(in) :: Density_start
      integer, intent(in) :: Velocity_start
      integer, intent(in) :: Temperature_start
      integer, intent(in) :: Species_start
      integer, intent(in) :: User_scalar_start
      Integer :: i,CC1,CC2
      integer :: n_cols

      RDIMN = merge(2,3,NO_K)

! Number of columns is equal to the number of delimiters + 1,
! so n_cols is initialized at one. Delimiters are assumed to be either:
! tabs    (ASCII code = 9)
! spaces  (ASCII code = 32)
! commas  (ASCII code = 44).

! First, replace commas and tabs by spaces
      DO I = 1,len(trim(first_line))
          CC1 = IACHAR(first_line(I:I))
          IF (CC1== 9) first_line(I:I) = ' '
          IF (CC1==44) first_line(I:I) = ' '
      ENDDO

! Remove leading spaces
      first_line = adjustl(first_line)

! A new column is detected if we find a space followed by a non-space
! character. Repeating spaces are ignored.
      n_cols = 1
      DO I = 1,len(trim(first_line))-1
         CC1 = IACHAR(first_line(I:I))
         CC2 = IACHAR(first_line(I+1:I+1))
         IF(CC1==32.AND.CC2==32) CYCLE  ! ignore 2 consecutive spaces
         IF(CC1==32.AND.CC2/=32) n_cols = n_cols + 1 ! A space followed by a non-space
      ENDDO

! Abort if we don't detect the expected number of columns
      IF(n_cols/=n_vars) THEN
         WRITE(ERR_MSG,'(A,I2/,A/,A)')      'Detected number of columns = ',n_cols, &
                                            'Incorrect number of columns in particle_input.dat file.',&
                                            'Please verify the data matches the following column layout:'

         WRITE(ERR_MSG(4),'(A)')            'Variable:             Column(s)'
! Coordinates option
         IF(rdimn==2) THEN
            WRITE(ERR_MSG(5),'(A)')         'Coordinates (x,y):    1 - 2'
         ELSE
            WRITE(ERR_MSG(5),'(A)')         'Coordinates (x,y,z):  1 - 3'
         ENDIF

! Phase ID option
         IF(Phase_ID_start>0) THEN
            WRITE(ERR_MSG(6),'(A,I2)')      'Phase ID:            ',Phase_ID_start
         ELSE
            WRITE(ERR_MSG(6),'(A)')         'Phase ID:             Not needed (uniform value)'
         ENDIF

! Diameter option
         IF(Diameter_start>0) THEN
            WRITE(ERR_MSG(7),'(A,I2)')      'Diameter:            ',Diameter_start
         ELSE
            WRITE(ERR_MSG(7),'(A)')         'Diameter:             Not needed (uniform value)'
         ENDIF

! Density option
         IF(Density_start>0) THEN
            WRITE(ERR_MSG(8),'(A,I2)')      'Density:             ',Density_start
         ELSE
            WRITE(ERR_MSG(8),'(A)')         'Density:              Not needed (uniform value)'
         ENDIF

! Velocity option
         IF(rdimn==2.AND.Velocity_start>0) THEN
            WRITE(ERR_MSG(9),'(A,I2,A,I2)') 'Velocity (u,v):      ',Velocity_start,' -',Velocity_start + 1
         ELSEIF(rdimn==3.AND.Velocity_start>0) THEN
            WRITE(ERR_MSG(9),'(A,I2,A,I2)') 'Velocity (u,v,w):    ',Velocity_start,' -',Velocity_start + 2
         ELSEIF(rdimn==2.AND.Velocity_start==0) THEN
            WRITE(ERR_MSG(9),'(A)')         'Velocity (u,v):       Not needed (uniform value)'
         ELSEIF(rdimn==3.AND.Velocity_start==0) THEN
            WRITE(ERR_MSG(9),'(A)')         'Velocity (u,v,w):     Not needed (uniform value)'
         ENDIF

! Temperature option
         IF(ENERGY_EQ.AND.Temperature_start>0) THEN
            WRITE(ERR_MSG(10),'(A,I2)')     'Temperature:         ',Temperature_start
         ELSE IF(ENERGY_EQ.AND.Temperature_start==0) THEN
            WRITE(ERR_MSG(10),'(A)')        'Temperature:          Not needed (uniform value)'
         ELSE
            WRITE(ERR_MSG(10),'(A)')        'Temperature:          Not needed (Energy equation is not solved)'
         ENDIF

! Species option
         IF(ANY_SOLIDS_SPECIES_EQ.AND.Species_start>0) THEN
            WRITE(ERR_MSG(11),'(A,I2,A,I2)') 'Species:             ',Species_start,' -',Species_start + DIMENSION_N_S - 1
         ELSEIF(ANY_SOLIDS_SPECIES_EQ.AND.Species_start==0) THEN
            WRITE(ERR_MSG(11),'(A)')        'Species:              Not needed (uniform value)'
         ELSE
            WRITE(ERR_MSG(11),'(A)')        'Species:              Not needed (Species equation(s) are not solved)'
         ENDIF

! User Scalar option
         IF(DES_USR_VAR_SIZE>0.AND.User_scalar_start>0) THEN
            WRITE(ERR_MSG(12),'(A,I2,A,I2)')'User scalars:        ',User_scalar_start,' -',User_scalar_start + DES_USR_VAR_SIZE - 1
         ELSEIF(DES_USR_VAR_SIZE>0.AND.User_scalar_start==0) THEN
            WRITE(ERR_MSG(12),'(A)')        'User scalars:         Not needed (uniform value)'
         ELSE
            WRITE(ERR_MSG(12),'(A)')        'User scalars:         Not needed (No user scalars defined)'
         ENDIF

         CALL LOG_ERROR()
      ELSE
         WRITE(ERR_MSG,'(A,I2)')'Detected number of columns = ',n_cols
         CALL LOG_STATUS()
      ENDIF

      RETURN

      end subroutine check_number_of_columns



      subroutine echo_all_filter_options( &
      PART_X_MIN,PART_X_MAX,PART_X_EXCLUDE, &
      PART_Y_MIN,PART_Y_MAX,PART_Y_EXCLUDE, &
      PART_Z_MIN,PART_Z_MAX, PART_Z_EXCLUDE, &
      PART_PHASE, &
      PART_DIAMETER_MIN, PART_DIAMETER_MAX,PART_DIAMETER_EXCLUDE, &
      PART_DENSITY_MIN, PART_DENSITY_MAX,PART_DENSITY_EXCLUDE, &
      PART_U_MIN, PART_U_MAX,PART_U_EXCLUDE, &
      PART_V_MIN, PART_V_MAX,PART_V_EXCLUDE, &
      PART_W_MIN, PART_W_MAX,PART_W_EXCLUDE, &
      PART_TEMP_MIN, PART_TEMP_MAX,PART_TEMP_EXCLUDE, &
      PART_X_S_MIN,PART_X_S_MAX,PART_X_S_EXCLUDE, &
      PART_USR_VAR_MIN,PART_USR_VAR_MAX,PART_USR_VAR_EXCLUDE)

      use discretelement
      use geometry, only: NO_K
      use param, only: dimension_n_s
      use run, only: any_solids_species_eq, energy_eq

      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      DOUBLE PRECISION, INTENT(IN) :: PART_X_MIN,PART_X_MAX
      LOGICAL, INTENT(IN) :: PART_X_EXCLUDE
      DOUBLE PRECISION, INTENT(IN) :: PART_Y_MIN,PART_Y_MAX
      LOGICAL, INTENT(IN) :: PART_Y_EXCLUDE
      DOUBLE PRECISION, INTENT(IN) :: PART_Z_MIN,PART_Z_MAX
      LOGICAL, INTENT(IN) :: PART_Z_EXCLUDE
      LOGICAL, INTENT(IN) :: PART_PHASE(DIM_M)
      DOUBLE PRECISION, INTENT(IN) :: PART_DIAMETER_MIN, PART_DIAMETER_MAX
      LOGICAL, INTENT(IN) :: PART_DIAMETER_EXCLUDE
      DOUBLE PRECISION, INTENT(IN) :: PART_DENSITY_MIN, PART_DENSITY_MAX
      LOGICAL, INTENT(IN) :: PART_DENSITY_EXCLUDE
      DOUBLE PRECISION, INTENT(IN) :: PART_U_MIN, PART_U_MAX
      LOGICAL, INTENT(IN) :: PART_U_EXCLUDE
      DOUBLE PRECISION, INTENT(IN) :: PART_V_MIN, PART_V_MAX
      LOGICAL, INTENT(IN) :: PART_V_EXCLUDE
      DOUBLE PRECISION, INTENT(IN) :: PART_W_MIN, PART_W_MAX
      LOGICAL, INTENT(IN) :: PART_W_EXCLUDE
      DOUBLE PRECISION, INTENT(IN) :: PART_TEMP_MIN, PART_TEMP_MAX
      LOGICAL, INTENT(IN) :: PART_TEMP_EXCLUDE
! Max size for PART species
      INTEGER, PARAMETER :: PART_X_size_max = 100
      DOUBLE PRECISION, INTENT(IN) :: PART_X_S_MIN(PART_X_size_max), PART_X_S_MAX(PART_X_size_max)
      LOGICAL, INTENT(IN) :: PART_X_S_EXCLUDE(PART_X_size_max)
      INTEGER, PARAMETER :: PART_USR_size_max = 100
      DOUBLE PRECISION, INTENT(IN) :: PART_USR_VAR_MIN(PART_USR_size_max), PART_USR_VAR_MAX(PART_USR_size_max)
      LOGICAL, INTENT(IN) :: PART_USR_VAR_EXCLUDE(PART_USR_size_max)

!-----------------------------------------------
! Read dimension: 2D vs 3D data
      integer :: RDIMN
      Integer :: i,m,n
      integer :: n_filters
      character(32) :: label
!-----------------------------------------------

      RDIMN = merge(2,3,NO_K)

         WRITE(ERR_MSG,'(A/,A)')'The following filters are applied to particle data: ','Default filter: Out of bound particles are automatically removed'
         CALL LOG_MESSAGE(__FILE__,__LINE__,LOGLEVEL_INFO, footer=.false.)

         n_filters = 0
! Particle position (x,y,z, coordinates)
         call echo_filter_option('x-coordinate',part_x_min,part_x_max,part_x_exclude,n_filters)
         call echo_filter_option('y-coordinate',part_y_min,part_y_max,part_y_exclude,n_filters)
         if(rdimn==3)call echo_filter_option('z-coordinate',part_z_min,part_z_max,part_z_exclude,n_filters)

! Particle Phase ID
         do M = 1,DES_MMAX
            if(.NOT.PART_PHASE(M)) then
               write(label,'("Phase ",I2)') M
               WRITE(ERR_MSG,'(A12,A,A)')trim(label),': EXCLUDE particles belonging to ', trim(label)
               CALL LOG_MESSAGE(__FILE__,__LINE__,LOGLEVEL_INFO, header=.false.,footer=.false.)
            endif
         end do

! Particle Diameter
         call echo_filter_option('Diameter',part_diameter_min,part_diameter_max,part_diameter_exclude,n_filters)

! Particle Density
         call echo_filter_option('Density',part_density_min,part_density_max,part_density_exclude,n_filters)

! Particle Velocity
         call echo_filter_option('x-velocity',part_u_min,part_u_max,part_u_exclude,n_filters)
         call echo_filter_option('y-velocity',part_v_min,part_v_max,part_v_exclude,n_filters)
         if(rdimn==3)call echo_filter_option('z-velocity',part_w_min,part_w_max,part_w_exclude,n_filters)

! Particle Temperature
         if(ENERGY_EQ) then
            call echo_filter_option('Temperature',part_temp_min,part_temp_max,part_temp_exclude,n_filters)
         endif

! Particle Species: Always need DIMENSION_N_S values for all particles
         if(ANY_SOLIDS_SPECIES_EQ) then
            do n = 1,DIMENSION_N_S
               write(label,'("Species",I2)') n
               call echo_filter_option(label,part_x_s_min(n),part_x_s_max(n),part_x_s_exclude(n),n_filters)
            enddo
         endif

! Particle Scalar: CAUTION THIS ARRAY IS TRANSPOSED, FIRST INDEX IS USR_VAR ID, SECOND INDEX IS PARTICLE ID
         if(DES_USR_VAR_SIZE>0) then
            do n = 1,DES_USR_VAR_SIZE
                write(label,'("User scalar # ",I2)') n
                call echo_filter_option(label,part_usr_var_min(n),part_usr_var_max(n),part_usr_var_exclude(n),n_filters)
            enddo
         endif

         if(n_filters==0) then
            WRITE(ERR_MSG,'(A)')'No other filter is applied.'
         else
            write(ERR_MSG, *)"" ! print final footer
         endif
         CALL LOG_MESSAGE(__FILE__,__LINE__,LOGLEVEL_INFO, header=.false.)


      RETURN

      end subroutine echo_all_filter_options

END MODULE READ_PART_INPUT_MOD
