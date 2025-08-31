#include "error.inc"

    MODULE keyframe_mod
!
!     Declare the user-defined namelist variables (usrnlst.inc) in this module.
!     Also Include user-defined variables in this module.  To access the
!     variables from a subroutine add the statement "Use usr".  If allocatable
!     arrays are defined in this module allocate them in usr0.  To turn on the
!     user defined subroutines (usr0, usr1, and usr2) set keyword CALL_USR to true.

      use error_manager
      use compar, only: mype, pe_io
      use mpi_utility, only: bcast
      use discretelement, only: dimn, cross
      use fs_util, only: file_exists
      use compar, only: mype, pe_io
      use exit, only: mfix_exit

      IMPLICIT NONE

! Maximum number of keyframes files
      integer, parameter :: max_kf = 100
! Maximum number of keyframes dependent variables (per keyframe file)
      integer, parameter :: max_kfvar = 20
! Maximum number of indices for a keyframe keyword
      integer, parameter :: max_kf_kw_indices = 3
! Flag to read keyframe file
      logical, dimension(max_kf) ::  keyframe
! Keyframe keyword
      character(len=64), dimension(max_kf, max_kfvar) :: kf_keyword
! Keyframe independent variable identifier (0=time, 1=x-coord., 2=y-coord., 3=z-coord)
      integer, dimension(max_kf) :: kf_iv
! Keyframe interpolation method
      character(len=10),dimension(max_kf) :: kf_interp
! Flag for dynamic wall boundaries
      logical, dimension(max_kf) ::  dynamic_bc
! Total number of keyframe tables
      integer :: kf_count
! Independent variable description (map integer to string)
     character(len=64), dimension(-1:3) :: kf_iv_description
! Rotational velocity unit
      CHARACTER(LEN=7) :: OMEGA_UNIT
! Conversion factor to express OMEGA in rad/sec
      DOUBLE PRECISION :: OMEGA_CONV_FACT

     data kf_iv_description / 'None', 'Time', 'x-coordinate', 'y-coordinate', 'z-coordinate' /

! Struct for storing keyframe data
      type kf_struct
        integer :: ID
        character(len=512) :: filename
        integer :: nrows, nvars
        character(len=10) :: interp_method
        character(len=64) :: x_name
        character(len=64), allocatable :: y_name(:)
        integer :: kw_id
        character(len=64),allocatable :: kw_name(:)
        integer, allocatable :: kw_index(:,:)
        integer :: iv
        integer :: x_index
        double precision, allocatable :: x(:)
        double precision, allocatable :: y(:,:)
        double precision, allocatable :: y_interp(:)
      end type kf_struct


      type(kf_struct), allocatable,  dimension(:) :: kf
      double precision, allocatable, dimension(:) :: interpolated_var
      integer, dimension(max_kf) :: list
      integer, dimension(max_kf) :: rlist

contains



!----------------------------------------------------------------------!
! Below is the code for keyframe data management.                      !
! It is anticipated that users won't need to modify it.                !
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: read_keyframe_data                                      !
!  Author: Sathish Sanjeevi                         Date:23-Oct-2019   !
!                                                                      !
!  Purpose: Keyframes are data tables provided by user to perform user-!
!           defined operations such as ramping or other time-dependent !
!           operations. This functions reads keyframe files provided   !
!           by the user.                                               !
!                                                                      !
!----------------------------------------------------------------------!

      subroutine read_keyframe_data


      IMPLICIT NONE
! Dummy arguments:
!---------------------------------------------------------------------
      integer i, j, k, eml, kwi
      integer id, fileunit, allocstatus, io, obi
      character(len = 8) :: xi
      character(len = 64) :: kwname
      double precision :: x_init

      i = 0
      list = 0
      rlist = 0
! Dry run to count actual number of keyframes
      do id = 1, max_kf
        if(keyframe(id) .eqv. .true.) then
          i = i+1
! Keyframe IDs can be non-contiguous. Efficient storage needs contiguous
! structs. Store indices based on KF IDs and vice versa.
          list(id) = i
          rlist(i) = id
        end if
      end do
      kf_count = i

      allocate(kf(kf_count), stat = allocstatus)


! Loop to read each keyframe file
      do i = 1, kf_count
         id = rlist(i)
         write(xi, '(I0.4)') id
         kf(i)%filename = '.KF/keyframe_'//trim(xi)//'.csv'

         if(myPE == pe_io) then
! Report error if file does not exist
            if (.not. file_exists(kf(i)%filename)) then
               write(err_msg, 9999)  trim(kf(i)%filename)
9999           format(2/,1x,70('*')/' from: keyframe_mod.f',/' error 9999: ',    &
                     'input file does not exist: ',a/,'aborting.',/1x,   &
                     70('*'),2/)
               call log_error()
            endif
            call get_number_of_rows_and_columns(kf(i)%filename, ',', kf(i)%nrows, kf(i)%nvars)
         endif

         call bcast(kf(i)%nrows, pe_io)
         call bcast(kf(i)%nvars, pe_io)

         allocate(kf(i)%y_name(kf(i)%nvars))
         allocate(kf(i)%x(kf(i)%nrows))
         allocate(kf(i)%y(kf(i)%nrows, kf(i)%nvars))
         allocate(kf(i)%y_interp(kf(i)%nvars))
         allocate(kf(i)%kw_name(kf(i)%nvars))
         allocate(kf(i)%kw_index(kf(i)%nvars, max_kf_kw_indices))


         if(myPE == pe_io) then
! Read keyframe file
            open(newunit = fileunit, file = kf(i)%filename)
            kf(i)%interp_method = kf_interp(id)
            read(fileunit, *, iostat = io) kf(i)%x_name,(kf(i)%y_name(k), k = 1, kf(i)%nvars)
            kf(i)%x_index = 1
! Actual reading of keyframe table
            do j = 1, kf(i)%nrows
               read(fileunit, *, iostat = io) kf(i)%x(j), &
                 (kf(i)%y(j,k), k = 1, kf(i)%nvars)
               ! check if file is in proper condition
               if (io>0) then
                  write(err_msg, 1111) j+3, kf(i)%filename
1111              format('Check if data is a valid floating point at line : ', i5/,&
                        'in file : ', a)
                  call log_error()
               else if (io<0) then
                  write(err_msg, 2222) kf(i)%filename
2222              format('End of file reached earlier in file : ', a)
                  call log_error()
               end if
            end do
            close(fileunit)

! Important: check if keyframe independent variable is sorted in ascending order
            x_init = 0.D0
            do j = 1, kf(i)%nrows
               if(kf(i)%x(j) .ge. x_init) then
                  x_init = kf(i)%x(j)
               else
                  write(err_msg, 9555) trim(kf(i)%filename), j+2
9555              format('Time cannot be decreasing. Check keyframe file : ',a/,&
                       'first column (time) at line : ', i5)
                  call log_error()
               end if
            end do

! Initialize starting keyframe variable to the first one
            kf(i)%y_interp(:) = kf(i)%y(1,:)
            kf(i)%x_index = 1

! Keyword name (defined in the .mfx file)
            kf(i)%kw_name(1:kf(i)%nvars) = kf_keyword(id,1:kf(i)%nvars)
! Keyword indices
            !!kf(i)%kw_index(1:kf(i)%nvars,1:max_kf_kw_indices) = kf_kw_index(id,1:kf(i)%nvars,1:max_kf_kw_indices) ! cgw

            kf(i)%iv = kf_iv(id)

            WRITE(ERR_MSG, 100) trim(kf(i)%filename),              &
                                ival(kf(i)%nrows),                 &
                                ival(kf(i)%nvars + 1),             &
                                kf(i)%interp_method,               &
                                trim(kf(i)%x_name),                &
                                trim(kf_iv_description(kf(i)%iv)), &
                                ival(minval(kf(i)%x(1:kf(i)%nrows))), &
                                ival(maxval(kf(i)%x(1:kf(i)%nrows)))
            EML=11

            do k = 1, kf(i)%nvars
               WRITE(ERR_MSG(EML), '("Dependent variable # ",A)') ival(k)   ; EML = EML + 1
               WRITE(ERR_MSG(EML), '(" - Name    = ",A)') trim(kf(i)%y_name(k))   ; EML = EML + 1
               kwname = kf(i)%kw_name(k)
               obi = index(kf(i)%kw_name(k),'(')
               if(obi>0) kwname = kf(i)%kw_name(k)(1:obi - 1)
               WRITE(ERR_MSG(EML), '(" - Keyword = ",A)') trim(kwname)   ; EML = EML + 1
! Get keyword indices, if any, and store them
               IF(trim(kf(i)%kw_name(k))/='None') THEN
                  call get_keyword_indices(trim(kf(i)%kw_name(k)), kf(i)%kw_index(k,:))
! Remove indices from the keyword name
                  kf(i)%kw_name(k) = kwname
                  DO kwi = 1,max_kf_kw_indices
                     IF(kf(i)%kw_index(k,kwi)>0) THEN
                        WRITE(ERR_MSG(EML), '(" - Index ",A," = ", A)') trim(ival(kwi)), ival(kf(i)%kw_index(k,kwi))   ; EML = EML + 1
                     ENDIF
                  ENDDO
               ENDIF
               WRITE(ERR_MSG(EML), '(" - Min     = ",A)')ival(minval(kf(i)%y(1:kf(i)%nrows,k))) ; EML = EML + 1
               WRITE(ERR_MSG(EML), '(" - Max     = ",A)')ival(maxval(kf(i)%y(1:kf(i)%nrows,k))) ; EML = EML + 1

            enddo


 100 FORMAT(/'Found keyframe file  : ',A,      &
            /'Number of rows       = ',A,      &
            /'Number of columns    = ',A,      &
            /'Interpolation method = ',A,      &
            /'Independent variable:',          &
            /' - Name = ',A,                   &
            /' - Type = ',A,                   &
            /' - Min  = ',A,                   &
            /' - Max  = ',A,                   &
            /)

            CALL LOG_INFO()
         end if ! myPE == pe_io
      end do !  keyframe file

      ! Broadcast keyframe data so each rank can do the interpolation
      do i = 1, kf_count
         call bcast(KF(i)%ID, pe_io)
         call bcast(KF(i)%filename, pe_io)
         call bcast(KF(i)%NROWS, pe_io)
         call bcast(KF(i)%NVARS, pe_io)
         call bcast(KF(i)%interp_method, pe_io)
         call bcast(KF(i)%x_name, pe_io)
         call bcast(KF(i)%y_name, pe_io)
         call bcast(KF(i)%kw_id, pe_io)
         call bcast(KF(i)%kw_name, pe_io)
         call bcast(KF(i)%kw_index, pe_io)
         call bcast(KF(i)%iv, pe_io)
         call bcast(KF(i)%x_index, pe_io)
         call bcast(KF(i)%x, pe_io)
         call bcast(KF(i)%y, pe_io)
      enddo

      call check_omega_unit

      return

      end subroutine read_keyframe_data



!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: interpolate_keyframe_data                               !
!  Author: Sathish Sanjeevi                         Date:23-Oct-2019   !
!                                                                      !
!  Purpose: Interpolate keyframe data based on current simulation      !
!           timestep                                                   !
!                                                                      !
!----------------------------------------------------------------------!


      subroutine interpolate_keyframe_data(x)
      IMPLICIT NONE
      double precision, intent(in) :: x
      integer :: i
      integer :: x_index, nrows
      character(len=10) :: interp
      double precision :: begin_x, end_x
      double precision :: kf_x, kf_next_x
      double precision :: denom
      double precision, allocatable :: kf_y(:), kf_next_y(:)
      logical :: flag


! Interpolate keyframe value (interp_y) based on current value of x
      do i = 1, kf_count
         allocate(kf_y(kf(i)%nvars))
         allocate(kf_next_y(kf(i)%nvars))
         nrows        = kf(i)%nrows
         interp       = kf(i)%interp_method
         begin_x      = kf(i)%x(1)
         end_x        = kf(i)%x(nrows)
         x_index      = kf(i)%x_index
         kf_x         = kf(i)%x(x_index)
         kf_next_x    = kf(i)%x(x_index + 1)
         kf_y(:)      = kf(i)%y(x_index,:)
         kf_next_y(:) = kf(i)%y(x_index + 1,:)

! Interpolate only if x lies within bounds of keyframe x-data
         if(x.gt.begin_x.and.x.lt.end_x) then
            flag = x.gt.kf(i)%x(x_index).and.&
                   x.le.kf(i)%x(x_index + 1)

! Loop repeatedly until x falls within the
! corresponding keyframe x intervals. This can happen if keyframe
! x increments are smaller than the simulation x-increments
            do while(.not.flag)
               x_index            = x_index + 1
               kf(i)%x_index = x_index
               kf_x               = kf(i)%x(x_index)
               kf_next_x          = kf(i)%x(x_index + 1)
               kf_y(:)            = kf(i)%y(x_index,:)
               kf_next_y(:)       = kf(i)%y(x_index + 1,:)

               flag = x.gt.kf(i)%x(x_index).and.&
                      x.le.kf(i)%x(x_index + 1)
            end do

! Interpolation can be a `linear` one or `step` (floor) function
            if(x.gt.kf_x.and.x.le.kf_next_x) then
               if(interp.eq.'LINEAR') then
                  denom = kf_next_x - kf_x
! Check if denominator is zero
                  if(denom.gt.0.D0) then
                     kf(i)%y_interp(:) = kf_y(:) + (x - kf_x)/denom  &
                                         *(kf_next_y(:) - kf_y(:))
                  else
                     kf(i)%y_interp(:) = kf_y(:)
                  end if
               elseif(interp.eq.'STEP') then
                  kf(i)%y_interp(:) = kf_y(:)
               else
                  write(*,5555)
 5555             FORMAT('Enter a valid interpolation scheme: LINEAR or STEP')
               endif
            endif
         endif

! Use the first or last keyframe value when x is outside keyframe bounds
         if(x.le.begin_x) kf(i)%y_interp(:) = kf(i)%y(1,:)
         if(x.ge.end_x)  kf(i)%y_interp(:) = kf(i)%y(nrows,:)


         deallocate(kf_y)
         deallocate(kf_next_y)
      enddo


      return
      end subroutine interpolate_keyframe_data





      subroutine assign_kf_value_to_keyword(i)
      use bc
      use is
      use discretelement, only: des_rigid_motion_vel, des_rigid_motion_omega, des_rigid_motion_rot_center

      IMPLICIT NONE
      integer, intent(in) :: i ! keyframe ID
      integer :: j ! keyframe dependent var ID
      integer :: bc_id ! Boundary Condition ID
      integer :: is_id ! Internal Surface ID
      integer :: xyz_id ! Coordinate ID (1=x, 2=y, 3=z)
      integer :: m ! Phase ID

      DO J = 1, KF(I)%NVARS ! loop over keyframe dependent variables

! Keyword names must be upper case in the CASE("***") statement

         SELECT CASE(TRIM(KF(I)%KW_NAME(J)))
         CASE('NONE')

         CASE('BC_U_G')
            bc_id = kf(i)%kw_index(j,1)
            bc_u_g(bc_id) = kf(i)%y_interp(j)

         CASE('BC_V_G')
            bc_id = kf(i)%kw_index(j,1)
            bc_v_g(bc_id) = kf(i)%y_interp(j)

         CASE('BC_W_G')
            bc_id = kf(i)%kw_index(j,1)
            bc_w_g(bc_id) = kf(i)%y_interp(j)

         CASE('BC_U_S')
            bc_id = kf(i)%kw_index(j,1)
            m = kf(i)%kw_index(j,2)
            bc_u_s(bc_id,m) = kf(i)%y_interp(j)

         CASE('BC_V_S')
            bc_id = kf(i)%kw_index(j,1)
            m = kf(i)%kw_index(j,2)
            bc_v_s(bc_id,m) = kf(i)%y_interp(j)

         CASE('BC_W_S')
            bc_id = kf(i)%kw_index(j,1)
            m = kf(i)%kw_index(j,2)
            bc_w_s(bc_id,m) = kf(i)%y_interp(j)

         CASE('BC_WALL_VEL')
            bc_id = kf(i)%kw_index(j,1)
            xyz_id = kf(i)%kw_index(j,2)
            bc_wall_vel(bc_id, xyz_id) = kf(i)%y_interp(j)

         CASE('BC_WALL_ROT_CENTER')
            bc_id = kf(i)%kw_index(j,1)
            xyz_id = kf(i)%kw_index(j,2)
            bc_wall_rot_center(bc_id, xyz_id) = kf(i)%y_interp(j)

         CASE('BC_WALL_OMEGA')
            bc_id = kf(i)%kw_index(j,1)
            xyz_id = kf(i)%kw_index(j,2)
            bc_wall_omega(bc_id, xyz_id) = kf(i)%y_interp(j)

         CASE('IS_WALL_VEL')
            is_id = kf(i)%kw_index(j,1)
            xyz_id = kf(i)%kw_index(j,2)
            is_wall_vel(is_id, xyz_id) = kf(i)%y_interp(j)

         CASE('IS_WALL_ROT_CENTER')
            is_id = kf(i)%kw_index(j,1)
            xyz_id = kf(i)%kw_index(j,2)
            is_wall_rot_center(bc_id, xyz_id) = kf(i)%y_interp(j)

         CASE('IS_WALL_OMEGA')
            is_id = kf(i)%kw_index(j,1)
            xyz_id = kf(i)%kw_index(j,2)
            is_wall_omega(is_id, xyz_id) = kf(i)%y_interp(j)

         CASE('DES_RIGID_MOTION_VEL')
            m = kf(i)%kw_index(j,1)
            xyz_id = kf(i)%kw_index(j,2)
            des_rigid_motion_vel(m, xyz_id) = kf(i)%y_interp(j)

         CASE('DES_RIGID_MOTION_OMEGA')
            m = kf(i)%kw_index(j,1)
            xyz_id = kf(i)%kw_index(j,2)
            des_rigid_motion_omega(m, xyz_id) = kf(i)%y_interp(j)

         CASE('DES_RIGID_MOTION_ROT_CENTER')
            m = kf(i)%kw_index(j,1)
            xyz_id = kf(i)%kw_index(j,2)
            des_rigid_motion_rot_center(m, xyz_id) = kf(i)%y_interp(j)


         CASE DEFAULT
            WRITE(ERR_MSG, 1000) TRIM(KF(I)%KW_NAME(J))
            CALL LOG_ERROR()
         END SELECT

      ENDDO
      RETURN

1000  FORMAT(2/,1X,70('*')/' FROM: KEYFRAME_MOD.F',/' ERROR 1000: ',&
      'UNKNOWN OR UNSUPPORTED KEYWORD: ', A/,&
      1x,'ABORTING.',/1X,70('*'),2/)
      end subroutine assign_kf_value_to_keyword

      subroutine update_keyframe_data
      use bc
      use run, only: time

      IMPLICIT NONE
      integer :: i ! keyframe ID

      call interpolate_keyframe_data(time)
      do i = 1, kf_count
         call assign_kf_value_to_keyword(i)
      end do

      return
      end subroutine update_keyframe_data


      subroutine get_keyword_indices(kw_name, kw_index)
      IMPLICIT NONE

      character(len=*), intent(in) :: kw_name
      integer, intent(out) :: kw_index(max_kf_kw_indices)
      integer :: obi ! open bracket index
      integer :: cbi ! close bracket index
      integer :: ci ! coma index
      integer :: ipos ! index position
      character(len=64) :: substring

      kw_index(:) = -1

      obi = index(kw_name,'(')
      cbi = index(kw_name,')')

      if(obi==0) then ! There is no opened bracket in the keyword name
                      ! Assume it means this is a scalar (no indices)
         return
      endif

! There is at least one index
! Get characters between the brackets
      substring = kw_name(obi + 1:cbi -1)

      ipos = 1
      ci   = index(substring,',')

! Get each index. Look at characters until a comma is found,
! and convert this substring to an integer.
! Repeat until there are no commas in substring.
      do while (ci>0)
         read(substring(1:ci - 1),*) kw_index(ipos)
         ipos = ipos + 1
! The index is the text between commas
         substring = substring(ci +1:)
         ci   = index(substring,',')
      enddo

! Last index
      read(substring,*) kw_index(ipos)

      return
      end subroutine get_keyword_indices


      subroutine get_number_of_rows_and_columns(filename, delimiter, nr, nc)
      IMPLICIT NONE
      CHARACTER(LEN=*), intent(in) :: filename
      CHARACTER(LEN=1), intent(in) :: delimiter
      integer, intent(out) :: nr
      integer, intent(out) :: nc
      integer :: fileunit, fs, io, i, j, ncp
      integer :: rl ! record length
      character :: C, CR, LF
      character(len=100) :: ioerrmsg
      logical :: blank_line

      ! if(mype/=pe_io) return

      LF = char(10) ! Line Feed
      CR = char(13) ! Carriage Return

! Report error if file does not exist
      if (.not. file_exists(filename)) then
         write(err_msg, 9999)  trim(filename)
9999     format(2/,1x,70('*')/' from: keyframe_mod.f',/' error 9999: ',    &
                     'input file does not exist: ',a/,'aborting.',/1x,   &
                     70('*'),2/)
         call log_error()
      endif
! open file
      open (newunit = fileunit, file=filename, form='unformatted', access='direct', status='old', recl = 1, iostat=io)
! fs is the file size, i.e. the number of characters in the file (record  length is one)
! fs is used to determine if we have reached the end of the file
      inquire(unit=fileunit, recl=rl, size=fs )
      !print*,'record, file size= ',rl,fs

      if(fs < 1) then
         write(err_msg, 900)  trim(filename), ival(fs)
900     format(2/,1x,70('*')/' from: keyframe_mod.f',/' error 1010: ',    &
                'Keyframe data file: ',a/,1x,&
                'File size: ',a/,1x,&
                'File size must be larger than 1',/1x,  &
                 70('*'),2/)
         call log_error()
      endif



! i keeps track of the file line number
! j keeps track of the number of characters read
      i  = 0
      j  = 0
      ioerrmsg = ''
! Read file one character at a time until EOF (j>fs)
      do
         i  = i + 1
         blank_line = .true.
         nc = 1
         C  = ''
! On each line, count the number of delimiters
         do while ((C .ne. CR) .and. (C .ne. LF))
            j = j + 1
            read(fileunit, rec=j, iostat=io, iomsg=ioerrmsg) C
            if(iachar(C)>32) blank_line = .false.
            !print*,C, io,j,fs, trim(ioerrmsg), iachar(C)
            if(j>fs) exit
            if (C .eq. delimiter) nc = nc + 1
         enddo
         if (C .eq. CR) then
             !write(*,*) "CR detected"
             j = j + 1
             read(fileunit, rec=j, iostat=io, iomsg=ioerrmsg) C
             if(iachar(C)>32) blank_line = .false.
             !write(*,*) "Read", ichar(C)
             if (C .ne. LF) then
                 write(err_msg, 1001) trim(filename), ival(i)
1001             format("Keyframe data file: ", a/ "Line: ", a/, "Invalid line ending")
                 call log_error()
             endif
         endif
         if(j>fs) exit
         if(i>1) then ! verify each line has the same number of columns
            if((nc .ne. ncp).and.(.not.blank_line)) then
                write(err_msg, 1000)  trim(filename), ival(i), ival(nc), ival(ncp)
1000            format(2/,1x,70('*')/' from: keyframe_mod.f',/' error 1000: ',    &
                     'Keyframe data file: ',a/,1x,&
                     'Line: ',a/,1x, &
                     'Number of detected columns: ',a/,1x,&
                     'Number of expected columns: ',a,/1x,   &
                     70('*'),2/)
                call log_error()
            endif
         endif
         if(.not.blank_line) ncp = nc
         !print*, 'Line', i, ' contains ', nc, ' columns ', 'blank_line=',blank_line
! Decrement number of rows when a blank line is encountered
! Note that only blank lines at the end of the file is supported
         if(blank_line) i = i - 1
      enddo

      close(fileunit)

      nr = i - 2

      if(iachar(C)>32) nr = nr + 1 ! correct number of rows if last line
                                   ! doesn't end with a NL
      nc = ncp - 1

! Data check

      if(nr < 1) then
         write(err_msg, 1010)  trim(filename), ival(nr)
1010     format(2/,1x,70('*')/' from: keyframe_mod.f',/' error 1010: ',    &
                'Keyframe data file: ',a/,1x,&
                'Number of detected rows: ',a/,1x,&
                'Must have at least one row (below header)',/1x,  &
                 70('*'),2/)
         call log_error()
      endif

      if((ncp < 2).and.(.not.blank_line)) then
         write(err_msg, 1020)  trim(filename), ival(ncp)
1020     format(2/,1x,70('*')/' from: keyframe_mod.f',/' error 1020: ',    &
                'Keyframe data file: ',a/,1x,&
                'Number of detected columns: ',a/,1x,&
                'Must have at least two columns.',/1x,  &
                 70('*'),2/)
         call log_error()
      endif


      return
      end subroutine get_number_of_rows_and_columns

      SUBROUTINE CHECK_OMEGA_UNIT
      USE CONSTANT, ONLY: PI
      USE PARAM1, ONLY: ZERO, ONE

      IMPLICIT NONE

      SELECT CASE(TRIM(OMEGA_UNIT))

         CASE('RAD/S')
            OMEGA_CONV_FACT = ONE

         CASE('DEG/S')
            OMEGA_CONV_FACT = ONE/180.0*PI

         CASE('REV/S')
            OMEGA_CONV_FACT = 2.0*PI

         CASE('RAD/MIN')
            OMEGA_CONV_FACT = ONE/60.0

         CASE('DEG/MIN')
            OMEGA_CONV_FACT = PI/180.0/60.0

         CASE('REV/MIN')
            OMEGA_CONV_FACT = 2.0*PI/60.0

         CASE DEFAULT
            WRITE(ERR_MSG, 1000) trim(OMEGA_UNIT)
            CALL LOG_ERROR()

      END SELECT

      RETURN

 1000 FORMAT('Error 1000: Unknown value for OMEGA_UNIT: ',A)

      END SUBROUTINE CHECK_OMEGA_UNIT

      END MODULE keyframe_mod
