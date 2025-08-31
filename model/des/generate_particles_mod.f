#include "error.inc"

MODULE GENERATE_PARTICLES

   USE error_manager

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PARTICLE_COUNT

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                C
!                                                                      C
!  Purpose: Generate particle configuration based on maximum particle  C
!           radius and filling from top to bottom within specified     C
!           bounds                                                     C
!                                                                      C
!                                                                      C
!  Authors: Rahul Garg                                Date: 19-Mar-14  C
!  Revised: Renjie Ke                                 Date: 11-Mar-24  C
!          add mpi support for GSP model when auto seeding             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE GENERATE_PARTICLE_CONFIG

      use mfix_pic, only: MPPIC
      use discretelement, only: PIP, MAX_PIP, PARTICLES, NGluedParticles, SuperDEM
! Flag indicating that the IC region is defined.
      use ic, only: IC_DEFINED
! Parameter for detecting unspecified values, zero, and one
      use param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE, Half
! Maximum number of initial conditions
      use param, only: DIMENSION_IC, dimension_n_s
! IC Region gas volume fraction.
      use ic, only: IC_EP_G
      use discretelement
      use mpi_utility
! grow the local arrays
      use des_allocate, only: PARTICLE_GROW, PARTICLE_GROW_GluedSphereDEM
      use resize
      !use desmpi_wrapper, only: DES_MPI_STOP
      use desgrid, only: dg_xstart, dg_ystart, dg_zstart
      use desgrid, only: dg_xend, dg_yend, dg_zend

      IMPLICIT NONE

      INTEGER :: ICV
      INTEGER :: totalNGSP
      integer :: lproc_parcnt(0:numpes-1) ! local component spheres count
      integer :: lproc_gspcnt(0:numpes-1) ! local gsp count
      integer :: lproc_gidsize(1:numpes) ! local gid size
      integer :: redundant
      integer :: bcast_pe
      integer :: lproc,lcurpar
      integer :: start_index,end_index,start_gid
      integer, allocatable, dimension(:) :: global_cid
      double precision, allocatable, dimension(:) :: global_mass
      double precision, allocatable, dimension(:) :: global_volume
      double precision, allocatable, dimension(:) :: global_dia
      double precision, allocatable, dimension(:) :: global_bound
      double precision, allocatable, dimension(:,:) :: global_pos
      double precision, allocatable, dimension(:,:) :: global_vel
      double precision, allocatable, dimension(:,:) :: global_inertia
      double precision, allocatable, dimension(:,:) :: global_quat
      integer, allocatable, dimension(:) :: global_gsp_state

! all glue variables are initialized to 0 in des_allocate_mod, not need to do it here.

      DO ICV = 1, DIMENSION_IC

         IF(.NOT.IC_DEFINED(ICV)) CYCLE
         IF(IC_EP_G(ICV) == ONE) CYCLE

         IF(MPPIC) THEN
            CALL GENERATE_PARTICLE_CONFIG_MPPIC(ICV)
         ELSE
            CALL GENERATE_PARTICLE_CONFIG_DEM(ICV)
         ENDIF
      ENDDO

      CALL GLOBAL_SUM(PIP,PARTICLES)

      WRITE(ERR_MSG, 1004) PARTICLES
 1004 FORMAT(/,'Total number of particles in the system: ',I15)

      CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)

! after auto seeding, it is highly possible each processor has different size of glued variables
      IF(GSP_EXPLICIT) THEN
         CALL GLOBAL_ALL_SUM(NGluedParticles,totalNGSP)
         WRITE(ERR_MSG, 1005) totalNGSP
         1005 FORMAT(/,'Total number of glued sphere in the system: ',I15)
         CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)

         if(numPEs > 1) then
            lproc_parcnt(:) = 0
            lproc_gspcnt(:) = 0
            lproc_parcnt(mype) = pip
            lproc_gspcnt(mype) = NGluedParticles
            ! might have use in the future
            call global_all_sum(lproc_parcnt)
            call global_all_sum(lproc_gspcnt)

            call bcast(particles, pe_io)
            call PARTICLE_GROW_GluedSphereDEM(totalNGSP)

            ! reset cglobal_id to zero in order to sort components global id
            if(allocated(cglobal_id)) cglobal_id(:) = ZERO

            start_index = 0
            start_gid = 0
            do lproc = 0,numPEs-1
               if(mype == lproc) exit
               start_index = start_index + lproc_parcnt(lproc)
               start_gid = start_gid + lproc_gspcnt(lproc)
            enddo

            end_index = start_index + lproc_parcnt(mype)
            start_index = start_index + 1

            ! during auto seeding, each pe have a copy of gid_list starting from 1, so reorder
            gid_list(1:PIP) = gid_list(1:PIP) + start_gid

            allocate(global_cid(particles))
            do lcurpar = 1, particles
               global_cid(lcurpar) = lcurpar
            enddo

            ! component global id should be continuous within a gsp for internal component spheres heat transfer
            cglobal_id(1:pip) = global_cid(start_index:end_index)

            do lcurpar = 1, pip
              where (gp_neighsid(lcurpar,:) /= -1) gp_neighsid(lcurpar,:) = gp_neighsid(lcurpar,:)+start_index-1
            enddo

            ! similar to manually seeding, build following global gsp variables
            ! call bcast(gluePos,pe_io)
            allocate(global_pos(totalNGSP,3)); global_pos(:,:) = 0.0D0
            global_pos(start_gid+1:start_gid+lproc_gspcnt(lproc),:) = gluePos(1:NGluedParticles,:)
            call global_all_sum(global_pos)
            gluePos(1:totalNGSP,:) = global_pos(:,:)

            ! call bcast(glueVel,pe_io)
            allocate(global_vel(totalNGSP,3)); global_vel(:,:) = 0.0D0
            global_vel(start_gid+1:start_gid+lproc_gspcnt(lproc),:) = glueVel(1:NGluedParticles,:)
            call global_all_sum(global_vel)
            glueVel(1:totalNGSP,:) = global_vel(:,:)

            ! call bcast(glueMass,pe_io)
            allocate(global_mass(totalNGSP)); global_mass(:) = 0.0D0
            global_mass(start_gid+1:start_gid+lproc_gspcnt(lproc)) = glueMass(1:NGluedParticles)
            call global_all_sum(global_mass)
            glueMass(1:totalNGSP) = global_mass(:)

	    ! call bcast(glueVolume,pe_io)
            allocate(global_volume(totalNGSP)); global_volume(:) = 0.0D0
            global_volume(start_gid+1:start_gid+lproc_gspcnt(lproc)) = glueVolume(1:NGluedParticles)
            call global_all_sum(global_volume)
            glueVolume(1:totalNGSP) = global_volume(:)

            ! call bcast(glueInertia,pe_io)
            allocate(global_inertia(totalNGSP,3)); global_inertia(:,:) = 0.0D0
            global_inertia(start_gid+1:start_gid+lproc_gspcnt(lproc),:) = glueInertia(1:NGluedParticles,:)
            call global_all_sum(global_inertia)
            glueInertia(1:totalNGSP,:) = global_inertia(:,:)

            ! call bcast(glueQuat,pe_io)
            allocate(global_quat(totalNGSP,4)); global_quat(:,:) = 0.0D0
            global_quat(start_gid+1:start_gid+lproc_gspcnt(lproc),:) = glueQuat(1:NGluedParticles,:)
            call global_all_sum(global_quat)
            glueQuat(1:totalNGSP,:) = global_quat(:,:)

            ! call bcast(glueDiameter,pe_io)
            allocate(global_dia(totalNGSP)); global_dia(:) = 0.0D0
            global_dia(start_gid+1:start_gid+lproc_gspcnt(lproc)) = glueDiameter(1:NGluedParticles)
            call global_all_sum(global_dia)
            glueDiameter(1:totalNGSP) = global_dia(:)

            ! call bcast(glueBounding,pe_io)
            allocate(global_bound(totalNGSP)); global_bound(:) = 0.0D0
            global_bound(start_gid+1:start_gid+lproc_gspcnt(lproc)) = glueBounding(1:NGluedParticles)
            call global_all_sum(global_bound)
            glueBounding(1:totalNGSP) = global_bound(:)

            ! make particle_state_gsp also global-wise
            ! since its index is the global gid, it still has to be a global array
            allocate(global_gsp_state(totalNGSP)); global_gsp_state(:) = 0
            global_gsp_state(start_gid+1:start_gid+lproc_gspcnt(lproc)) = particle_state_gsp(1:NGluedParticles)
            call global_all_sum(global_gsp_state)
            particle_state_gsp(1:totalNGSP) = global_gsp_state(:)

            ! those should be reset to zero as well
            glueForce(:,:) = 0.0D0
            glueTorque(:,:) = 0.0D0
            glueOmg(:,:) = 0.0D0
            glueAngMom(:,:) =0.0D0

            NGluedParticles = totalNGSP

            deallocate(global_cid)
            deallocate(global_mass)
            deallocate(global_dia)
            deallocate(global_bound)
            deallocate(global_pos)
            deallocate(global_vel)
            deallocate(global_inertia)
            deallocate(global_quat)
            deallocate(global_gsp_state)

            ! for auto seeding in mpi, call the following subroutine to correct it's mother PE
            CALL CORRECT_COMPONENTS_MOTHER_PE
         endif ! end of numPEs > 1
      ENDIF ! end of GluedSphereDEM

      RETURN
      END SUBROUTINE GENERATE_PARTICLE_CONFIG

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: CORRECT_COMPONENTS_MOTHER_PE                            C
!                                                                      C
!  Purpose: allow auto seeding to put GSP across the mpi boundary      C
!           assign component out of current PE bound a new PE          C
!           set this component as non-exsitent in current PE           C
!                                                                      C
!  Author: Renjie Ke                                 Date: 04-June-24  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CORRECT_COMPONENTS_MOTHER_PE

      use mfix_pic, only: MPPIC
      use discretelement, only: PIP, MAX_PIP, PARTICLES, NGluedParticles, SuperDEM
! Flag indicating that the IC region is defined.
      use ic, only: IC_DEFINED
! Parameter for detecting unspecified values, zero, and one
      use param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE, Half
! Maximum number of initial conditions
      use param, only: DIMENSION_IC, dimension_n_s
! IC Region gas volume fraction.
      use ic, only: IC_EP_G
      use discretelement
      use mpi_utility
! grow the local arrays
      use des_allocate, only: PARTICLE_GROW, PARTICLE_GROW_GluedSphereDEM
      use resize

      use desgrid, only: dg_xstart, dg_ystart, dg_zstart
      use desgrid, only: dg_xend, dg_yend, dg_zend

      use run, only: any_solids_species_eq, energy_eq
      use des_rxns, only: des_x_s
      use des_thermo, only: des_t_s
      use desmpi_wrapper, only: des_mpi_barrier, des_mpi_stop
      use functions, only: is_nonexistent, set_nonexistent, SET_NORMAL

      IMPLICIT NONE

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: correction_table
      INTEGER :: component_variable_size
      logical :: x_test, y_test, z_test, xyz_test
      INTEGER :: cstart, new_pip, locspot, gcid, lcurpar,new_max_pip
      INTEGER :: lrow, nrows
      DOUBLE PRECISION :: temp_pos(3)

      ! if_change_pe_flag, component_original_rank_id(1), xyz(3), phase_id(1), radius(1), density(1), vel(3), sc2gpc_vec(3), cglobal_id(1), gid_list(1)
      ! if_change_pe_flag == 1 means this component is outside of current PE boundary, it should change its rank id
      ! total 16 columns initially
      component_variable_size = 16

      IF(ENERGY_EQ) THEN
      ! Temperature
      component_variable_size = component_variable_size + 1
      ! gp_neighsid(6),gp_neighsa(6),gp_sa(1)
      component_variable_size = component_variable_size + 13
      ENDIF

      IF(ANY_SOLIDS_SPECIES_EQ) THEN
      ! species
      component_variable_size = component_variable_size + dimension_n_s
      ENDIF

      IF(DES_USR_VAR_SIZE>0) THEN
      ! des_usr_var
      component_variable_size = component_variable_size + des_usr_var_size
      ENDIF

      ! initialize the correction_table and set to ZERO
      ALLOCATE(correction_table(particles, component_variable_size))
      correction_table(:,:) = ZERO

      ! save the original PIP and MAX_PIP in each rank
      new_pip = pip
      new_max_pip = max_pip

      ! first loop, iterate components in current PE to see if any component is outside of current PE boundary
      DO lcurpar = 1, PIP
         IF(particle_state(lcurpar) .lt. 1) CYCLE

         x_test = (des_pos_new(lcurpar,1) .ge. xe(istart1_all(mype)-1) .and. &
               des_pos_new(lcurpar,1) .lt. xe(iend1_all(mype)))
         y_test = (des_pos_new(lcurpar,2) .ge. yn(jstart1_all(mype)-1) .and. &
               des_pos_new(lcurpar,2) .lt. yn(jend1_all(mype)))
         xyz_test = x_test .and. y_test

         IF(do_k) THEN
               z_test = (des_pos_new(lcurpar,3).ge.zt(kstart1_all(mype)-1).and. &
                           des_pos_new(lcurpar,3).lt.zt(kend1_all(mype)))
               xyz_test = xyz_test.and.z_test
         ENDIF

         IF(xyz_test) CYCLE

         CALL set_nonexistent(lcurpar)
         new_pip = new_pip - 1

         ! filling the correction_table with components info
         gcid = cglobal_id(lcurpar)
         correction_table(gcid,1) = 1
         correction_table(gcid,2) = mype
         correction_table(gcid,3:5) = des_pos_new(lcurpar,:)
         correction_table(gcid,6) = PIJK(lcurpar,5)
         correction_table(gcid,7) = des_radius(lcurpar)
         correction_table(gcid,8) = RO_SOL(lcurpar)
         correction_table(gcid,9:11) = des_vel_new(lcurpar,:)
         correction_table(gcid,12:14) = sc2gpc_vec(lcurpar,:)
         ! global id
         correction_table(gcid,15) = cglobal_id(lcurpar)
         ! global glued sphere id
         correction_table(gcid,16) = gid_list(lcurpar)

         ! after saving, reset variables related to local non-exsitence
         des_pos_new(lcurpar,:) = ZERO
         PIJK(lcurpar,5) = ZERO
         des_radius(lcurpar) = ZERO
         RO_SOL(lcurpar) = ZERO
         des_vel_new(lcurpar,:) = ZERO
         sc2gpc_vec(lcurpar,:) = ZERO
         sc2gpc_vec(lcurpar,:) = ZERO
         cglobal_id(lcurpar) = ZERO
         gid_list(lcurpar) = ZERO

         cstart = 17

         IF(energy_eq) THEN
               correction_table(gcid,cstart) = des_t_s(lcurpar)
               correction_table(gcid,cstart+1:cstart+6) = gp_neighsid(lcurpar,:)
               correction_table(gcid,cstart+7:cstart+12) = gp_neighsa(lcurpar,:)
               correction_table(gcid,cstart+13) = gp_sa(lcurpar)
               cstart = cstart + 14

               des_t_s(lcurpar) = ZERO
               gp_neighsid(lcurpar,:) = ZERO
               gp_neighsa(lcurpar,:) = ZERO
               gp_sa(lcurpar) = ZERO
         ENDIF

         IF(any_solids_species_eq) THEN
               correction_table(gcid,cstart:cstart+DIMENSION_N_S-1) = des_x_s(lcurpar,1:DIMENSION_N_S)
               cstart = cstart + DIMENSION_N_S

               des_x_s(lcurpar,1:DIMENSION_N_S) = ZERO
         ENDIF

         IF(des_usr_var_size>0) THEN
               correction_table(gcid,cstart:cstart+des_usr_var_size-1) = des_usr_var(1:des_usr_var_size,lcurpar)
               cstart = cstart + des_usr_var_size

               des_usr_var(1:des_usr_var_size,lcurpar) = ZERO
         ENDIF
      ENDDO ! end of lcurpar = 1, PIP - first loop

      CALL GLOBAL_ALL_SUM(correction_table)

      nrows = SIZE(correction_table,1)
      locspot = 1

      ! second loop, every rank iterate rows of correction table to find a component within its boundary box
      DO lrow = 1, nrows ! 1, particles

         ! if_change_pe_flag is false or mype == component_original_rank_id CYCLE
         IF(correction_table(lrow,1) .lt. 1.0 .or. correction_table(lrow,2) == mype) CYCLE
         temp_pos(:) = correction_table(lrow,3:5)

         ! check if this component is within current rank boundary
         x_test = (temp_pos(1) .ge. xe(istart1_all(mype)-1) .and. &
               temp_pos(1) .lt. xe(iend1_all(mype)))
         y_test = (temp_pos(2) .ge. yn(jstart1_all(mype)-1) .and. &
               temp_pos(2) .lt. yn(jend1_all(mype)))
         xyz_test = x_test .and. y_test

         IF(do_k) THEN
               z_test = (temp_pos(3).ge.zt(kstart1_all(mype)-1).and. &
                           temp_pos(3).lt.zt(kend1_all(mype)))
               xyz_test = xyz_test.and.z_test
         ENDIF

         IF(.NOT. xyz_test) CYCLE

         ! if true, add pip count
         new_pip = new_pip + 1

         ! grow local pip arrays in case not enough space
         CALL PARTICLE_GROW(new_pip)

         ! find the first location with non-exsitence
         DO WHILE(.not.is_nonexistent(locspot))
            locspot = locspot + 1
         ENDDO

         ! set this locspot as a normal particle
         CALL SET_NORMAL(locspot)

         ! copy the value from correction_table and fill the local arrays
         des_pos_new(locspot,:) = correction_table(lrow,3:5) !xyz
         PIJK(locspot,5) = int(correction_table(lrow,6)) ! phase id
         des_radius(locspot) = correction_table(lrow,7) ! component radius
         RO_SOL(locspot) = correction_table(lrow,8) ! density
         des_vel_new(locspot,:) = correction_table(lrow,9:11) ! velocity
         sc2gpc_vec(locspot,:) = correction_table(lrow,12:14) ! componet2gpc vector
         cglobal_id(locspot) = correction_table(lrow,15) ! component global id
         gid_list(locspot) = int(correction_table(lrow,16)) ! global gid

         cstart = 17
         IF(energy_eq) THEN
               des_t_s(locspot) = correction_table(lrow,cstart)
               gp_neighsid(locspot,:) = int(correction_table(lrow,cstart+1:cstart+6))
               gp_neighsa(locspot,:) = correction_table(lrow,cstart+7:cstart+12)
               gp_sa(locspot) = correction_table(lrow,cstart+13)
               cstart = cstart + 14
         ENDIF

         IF(any_solids_species_eq) THEN
               des_x_s(locspot,1:dimension_n_s) = correction_table(lrow,cstart:cstart+DIMENSION_N_S-1) ! species
               cstart = cstart + DIMENSION_N_S
         ENDIF

         IF(des_usr_var_size>0) THEN
               des_usr_var(1:des_usr_var_size,locspot) = correction_table(lrow,cstart:cstart+des_usr_var_size-1) ! des_usr_var
               cstart = cstart + des_usr_var_size
         ENDIF

      ENDDO ! end of lcurpar = 1, nrows

      ! deallocate the correction_table
      IF(allocated(correction_table)) deallocate(correction_table)
      ! reset the local PIP and MAX_PIP
      PIP = NEW_PIP
      MAX_PIP = MAX(PIP, MAX_PIP)

      END SUBROUTINE CORRECT_COMPONENTS_MOTHER_PE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG_DEM                            !
!  Authors: Rahul Garg, Jeff Dietiker                Date: 21-Mar-2014 !
!                                                                      !
!  Purpose: Generate particle configuration for DEM solids for each IC !
!           region. Now using the particle linked lists for initial    !
!           build                                                      !
!           This routine will ultimately supersede the older routine   !
!           that has not been deleted yet                              !
!                                                                      !
!----------------------------------------------------------------------!
!   Revision Date    : 09/2016                                         !
!   Revision Purpose : Make data structure changes for polydispersity  !
!                                                                      !
!   Revision By      : ASU MFIX-DEM Phi Team (Shaohua Chen, Yang Jiao, !
!                      Aytekin Gel, Manogna Adepu, Heather Emady)      !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GENERATE_PARTICLE_CONFIG_DEM(ICV)

! Global Variables:
!---------------------------------------------------------------------//
! particle radius and density
      use discretelement, only: DES_RADIUS, RO_Sol, PVOL
! SQP particle half-axis
      use discretelement, only: SQP_a, SQP_b, SQP_c
! particle position new and old
      use discretelement, only: DES_POS_NEW, DES_POS_OLD
! particle velocity new and old
      use discretelement, only: DES_VEL_NEW, DES_VEL_OLD
! Number of particles in the system (current)
      use discretelement, only: PIP
! Number of DEM solids phases.
      use discretelement, only: DES_MMAX
! Flag to use _OLD variables
      use discretelement, only: DO_OLD
! Angular velocity
      use discretelement, only: OMEGA_OLD, OMEGA_NEW, PIJK
! solid phase diameters and densities.
      use physprop, only: D_p0, RO_s0, MMAX
! IC Region solids volume fraction.
      use ic, only: IC_EP_s, IC_RO_s
! particle size distribution:
      use ic, only: IC_PSD_TYPE
      use ic, only: IC_PSD_MEAN_DP, IC_PSD_STDEV
      use ic, only: IC_PSD_MIN_DP
      use ic, only: IC_PSD_MAX_DP

! Constant: 3.14159...
      use constant, only: PI
! min and max physical coordinates of IC regions in each direction
      use ic, only: IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t
! initially specified velocity field and granular temperature
      use ic, only: IC_U_s, IC_V_s, IC_W_s, IC_Theta_M
      use ic, only: IC_DES_LATTICE, IC_DES_SPACING,IC_DES_SM,IC_DES_NP
      use ic, only: IC_DES_ABORT_ON_FAILED_SEEDING
      use ic, only: IC_DES_SPACE_FACTOR_X, IC_DES_SPACE_FACTOR_Y, IC_DES_SPACE_FACTOR_Z
      use ic, only: IC_DES_RAND, IC_DES_RAND_FACTOR_X, IC_DES_RAND_FACTOR_Y
      use ic, only: IC_DES_RAND_FACTOR_Z
      use ic, only: IC_DES_CHECK_STL_OVERLAP
      use ic, only: IC_WITH_DES_SOLID_PHASE
      use ic, only: DES_IC_POS_TMP
      use ic, only: IC_SQP_RANDOM_Q, IC_GSP_RANDOM_Q
! Parameter for detecting unspecified values, zero, and one
      use param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE, Half
! Parameter for small numbers
      use param1, only: SMALL_NUMBER

! to access random number generator subroutines
      use randomno
      use mpi_utility

      use desgrid, only: dg_xstart, dg_ystart, dg_zstart
      use desgrid, only: dg_xend, dg_yend, dg_zend

! direction wise spans of the domain and grid spacing in each direction
      use geometry, only: xlength, ylength, zlength

      use cutcell, only: CARTESIAN_GRID
      use stl_functions_des, only: CHECK_IF_PARTICLE_OVERLAPS_STL
      use run, only: solids_model
      use des_allocate, only: PARTICLE_GROW, PARTICLE_GROW_GluedSphereDEM

      use desgrid, only: IofPOS, JofPOS, KofPOS
      use desgrid, only: dg_is_ON_myPE_OWNs
      use toleranc, only: compare

      use discretelement, only: max_pip, max_radius, xe, yn, zt, pbb
      use discretelement, only: cgp_scaling_method, cgp_d_p0
      use discretelement, only: sphereNumber, sphereRelDia, sphereRelPos, NGluedParticles, gid_list
      use discretelement, only: gp_sa, gp_neighsid, gp_neighsa, gsp_vol
      use discretelement, only: gsp_q1, gsp_q2, gsp_q3, gsp_q4
      use discretelement, only: glueMass, glueVolume, glueDiameter, glueBounding, gluePos, glueVel, glueOmg, glueAngMom, glueInertia
      use discretelement, only: glueForce, glueTorque, glueNumber, gluedSphereDEM
      use discretelement, only: PARTICLE_STATE
      use discretelement, only: MIN_OVERLAP, MAX_OVERLAP, MEAN_OVERLAP, COORDINATION
      use discretelement, only: gsp_explicit, gsp_implicit
      use functions
      use param, only: dim_m
      use param, only: dimension_i, dimension_j, dimension_k
      use param1, only: undefined
      use mass_outflow_dem_mod, only: delete_particle
      use des_psd

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: ICV

! Local variables
!---------------------------------------------------------------------//
! Starting positions in the axial directions
      DOUBLE PRECISION :: xINIT, yINIT, zINIT
! Fractor used to scale particle diameter
      DOUBLE PRECISION :: lFAC,lFAC_X,lFAC_Y,lFAC_Z
! Particle position and velocity
      DOUBLE PRECISION :: POS(3), VEL(3)
! Number of particles in the lattice
      INTEGER :: SEED_X, SEED_Y, SEED_Z
! Loop indices phase/fluid cell
      INTEGER :: M, MM, I, J, K, IJK
! Loop indices for seeding
      INTEGER :: II, JJ, KK, LL
! Start and end bound for IC region.
      DOUBLE PRECISION :: IC_START(3), IC_END(3)
! Volume and lengths of the IC Region
      DOUBLE PRECISION :: DOM_VOL, DOML(3)
! Flag to skip the particle
      LOGICAL :: SKIP
! Diameter,radius adjusted for space padding
      DOUBLE PRECISION :: ADJ_DIA, ADJ_RADIUS
      DOUBLE PRECISION :: ADJ_DIA_X, ADJ_RADIUS_X
      DOUBLE PRECISION :: ADJ_DIA_Y, ADJ_RADIUS_Y
      DOUBLE PRECISION :: ADJ_DIA_Z, ADJ_RADIUS_Z
! Number of particles calculated from volume fracton
      INTEGER :: rPARTS(DIM_M), tPARTS
! Spacing between particles.
      DOUBLE PRECISION :: lDEL, lDX, lDY, lDZ
! Number of seeded particles
      INTEGER :: pCOUNT(DIM_M), tCOUNT, pSEEDED(DIM_M)


      DOUBLE PRECISION :: SOLIDS_DATA(0:DIM_M), lRO(DIM_M), GLOBAL_SOLIDS_DATA !(0:DIM_M)

      LOGICAL :: VEL_FLUCT, XYZ_FLUCT
      DOUBLE PRECISION, ALLOCATABLE :: randVEL(:,:), randXYZ(:,:)

      !DOUBLE PRECISION :: BC_MAX_RADIUS

      DOUBLE PRECISION :: rfac, rfac_base, rfac_x, rfac_y, rfac_z
      DOUBLE PRECISION :: IC_VOL, P_VOL

      INTEGER :: LATTICE
      INTEGER, PARAMETER :: LATTICE_CUBIC = 1 , LATTICE_HEXA = 2

      LOGICAL :: NARROW_Z

      INTEGER :: ALL_PART_COUNT_INIT(0:numPEs-1)
      INTEGER :: ALL_PART_COUNT_ADJ(0:numPEs-1)

      INTEGER :: iproc, IERR
      INTEGER :: NP_to_remove, NR
      INTEGER :: OLD_PIP,PIP_INCREMENT, NEW_PIP
      DOUBLE PRECISION :: ECHO_EPS, ECHO_SM
      INTEGER :: ECHO_NP, ECHO_NP_SPHERE
      CHARACTER(LEN=11) :: ECHO_INPUT_VAR, ECHO_INPUT_VALUE

! A variable for temporary storing particle diameter
      DOUBLE PRECISION, DIMENSION(1) :: TMP_DIAMETER

      DOUBLE PRECISION, DIMENSION(1) :: tmp_rand_number
      DOUBLE PRECISION :: tmp_prob_sum
! Variable to store volumetric mean diameter
      DOUBLE PRECISION :: vol_mean_diameter
      INTEGER :: phase_number
      INTEGER, ALLOCATABLE :: tmp_array(:)
!      DOUBLE PRECISION, ALLOCATABLE :: prob_for_each_phase(:)
! Max diameter for a given phase M for lattice spacing
      double precision :: phase_max_dia
      double precision :: total_volume
      integer :: allocstatus
! PIP form previous IC
      integer :: PIP_PREVIOUS
! Particle bounding box in the current IC
      double precision :: icpbb(3)
! Is the GSP baseline shape rotated
      logical :: gsp_rotated



      total_volume = 0.D0
!      totaltime = 0.D0
!......................................................................!

      WRITE(ERR_MSG,"(2/,'Generating particle configuration for initial condition', I3)") ICV

      CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.TRUE., FOOTER=.FALSE.)

! Keep a copy of PIP from the previous IC seeding
      PIP_PREVIOUS = PIP

      SOLIDS_DATA = ZERO
      CALL GET_IC_VOLUME(ICV, SOLIDS_DATA(0))
      IC_VOL = SOLIDS_DATA(0)
      CALL GLOBAL_ALL_SUM(IC_VOL)

! setting particle seed spacing grid to be slightly greater than
! the maximum particle diameter. seed at ~particle radii
      lFAC = 1.05D0

      ! Verify that the IC region volume is not zero.
      IF(IC_VOL <= 0.0d0) THEN
         WRITE(ERR_MSG,1000) ICV, IC_VOL
         CALL LOG_ERROR()
      ENDIF

1000  FORMAT('Error 1000: Invalid IC region volume: IC=',I3,' VOL=', ES15.4)

! Setup local arrays with IC region bounds.
      IC_START(1)=IC_X_W(ICV);   IC_END(1)=IC_X_E(ICV)
      IC_START(2)=IC_Y_S(ICV);   IC_END(2)=IC_Y_N(ICV)
      IC_START(3)=IC_Z_B(ICV);   IC_END(3)=IC_Z_T(ICV)

      DOML = IC_END-IC_START
      IF(NO_K) DOML(3)=DZ(1)

      ! Create a local array of solids densities so we don't have to
      ! reevaluate the function for every particle. Doing this for each
      ! phase is a bit overkill as DEM says that each IC region can only
      ! have on solids.
      DO M=MMAX+1,MMAX+DES_MMAX
         IF((SOLIDS_MODEL(M) == 'DEM'.OR.SOLIDS_MODEL(M) == 'CGP'.OR.SOLIDS_MODEL(M) == 'SQP'.OR.SOLIDS_MODEL(M) == 'GSP') &
            .and. IC_EP_S(ICV,M) > ZERO) THEN
            lRO(M) = IC_RO_S(ICV, M)
         ELSE
            lRO(M) = -UNDEFINED
         END IF
      END DO

      ! Volume of the IC region
      DOM_VOL = DOML(1)*DOML(2)*DOML(3)

      phase_number = 0
      rPARTS=0
      !BC_MAX_RADIUS = ZERO

      DO M=MMAX+1,MMAX+DES_MMAX
! HangZ: add glued sphere check
         IF(SOLIDS_MODEL(M) == 'DEM' .OR. SOLIDS_MODEL(M) == 'CGP' .OR. SOLIDS_MODEL(M) == 'SQP' .OR. SOLIDS_MODEL(M) == 'GSP') THEN
            phase_number = phase_number + 1

! Mean volume of particles given a distribution is not same as volume based on
! distribution's mean diameter. mean_volume() computes actual mean volume for a given distribution
            if(IC_PSD_TYPE(ICV,M) .eq. 'MONO') then
              p_vol = (pi*d_p0(m)**3.0D0)/6.0d0
              if(SuperDEM) p_vol = sqp_vol(m)
	      if(GluedSphereDEM) p_vol = gsp_vol(m)
            elseif(CGDEM.and.cgp_scaling_method(m)==2) then
              p_vol = (pi*cgp_d_p0(m)**3.0D0)/6.0d0
            else
              p_vol = mean_volume(1, ICV, M)  ! 1 selects IC; 2 selects BC
! Particle size distribution is given based on bounding diameter of glued particle.
! To get the actual volume of gsp, mean_volume needs to be converted based on the ratio
! between actual volume of gsp (gsp_vol) and volume from bounding diameter.
	      if(GluedSphereDEM) p_vol = p_vol*gsp_vol(m)/((pi*d_p0(m)**3.0D0)/6.0d0)
            end if
            call global_all_sum(p_vol)
            p_vol = p_vol/dble(numPEs)

! Use the computed particle volume to compute number of particles
! if solids volume fraction is given
            IF(IC_EP_S(ICV,M)>ZERO.AND.IC_EP_S(ICV,M)/=UNDEFINED) THEN
               rPARTS(M) = floor((IC_EP_S(ICV,M)*IC_VOL)/P_VOL)
! if solids mass is given
            ELSEIF(IC_DES_SM(ICV,M)>ZERO.AND.IC_DES_SM(ICV,M)/=UNDEFINED) THEN
               rPARTS(M) = floor(IC_DES_SM(ICV,M)/(RO_S0(M)*P_VOL))
! if number of particles are specified explicitly
            ELSEIF(IC_DES_NP(ICV,M)>0.AND.IC_DES_NP(ICV,M)/=UNDEFINED_I) THEN
               rPARTS(M) = IC_DES_NP(ICV,M)
            ELSE
               rPARTS(M) = 0
            ENDIF

! Negative values of rParts indicate overflow. Treat this as a fatal error.
            IF(rPARTS(M)<0) THEN
               WRITE(ERR_MSG,100) ICV, M, DOM_VOL,D_P0(M),&
                    (6.0d0*IC_EP_S(ICV,M)*DOM_VOL)/(PI*(D_P0(M)**3))
               CALL LOG_ERROR()
            ENDIF

! @JFD: max_radius is not used here in this function. Should be removed ?!?
! already, max_radius is searched and updated in
            ! IF(IC_EP_S(ICV,M) > ZERO .AND. D_P0(M) > (2.D0*MAX_RADIUS)) then
            !    MAX_RADIUS = HALF * D_P0(M)
            !    print *, "higher max_radius : ", max_radius
            !    print *, "ICV and phase M : ", ICV, M
            !    print *, "diameter d_p0 : ", d_p0(m)
            !    read (*,*)
            ! ENDIF
         ENDIF
      ENDDO


100 FORMAT('Error 100: Overflow in IC region, IC=',I3,' , M=',I3, &
      /'IC volume=',ES15.4,' , Particle diameter=', ES15.4, ' , Number of particles = ',ES15.4, &
      /'This error is usually triggered by attempting to initialize with a very large number of particles.' &
      /'Please verify your settings (Particle size and region extents).', &
      /'This is a fatal error. MFiX will exit now.')


      allocate(tmp_array(phase_number))
!      allocate(prob_for_each_phase(phase_number))

! Total number of particles in this IC region.
      tPARTS = sum(rPARTS)
      IF(tPARTS == 0) RETURN


      xINIT = IC_START(1)
      yINIT = IC_START(2)
      zINIT = IC_START(3)

      pCOUNT = 0
      pSEEDED = 0
      tCOUNT = 0

!     get particle number fraction for each phase
      phase_number = 0
      DO M=MMAX+1,MMAX+DES_MMAX
! HangZ: add glued sphere check
         IF(SOLIDS_MODEL(M) == 'DEM'.OR.SOLIDS_MODEL(M) == 'CGP' &
            .OR. SOLIDS_MODEL(M) == 'SQP' .OR. SOLIDS_MODEL(M) == 'GSP') THEN
               phase_number = phase_number + 1
               tmp_array(phase_number) = M
!               prob_for_each_phase(phase_number) = dble(rPARTS(M))/dble(tPARTS)
         ENDIF
      ENDDO

! Fill particles for each solids phase
      DO MM=1, PHASE_NUMBER
         M = tmp_array(MM)
! Only seed particles if number of particles is non-zero
         if(rparts(M)==0) cycle

! Seed particles from cubic of hexagonal arrangement
            WRITE(ERR_MSG,"(2/,'Phase',I3,': Number of particles to seed = ',I12)") M,rpARTS(M)
            CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)

            GLOBAL_SOLIDS_DATA = ZERO

            if(IC_PSD_TYPE(icv, m) == 'LOG_NORMAL' .or.&
               IC_PSD_TYPE(icv, m) == 'NORMAL'.or. &
               IC_PSD_TYPE(icv, m) == 'UNIFORM') then
               phase_max_dia = IC_PSD_MAX_DP(icv, m)
            else if(IC_PSD_TYPE(icv, m) == 'CUSTOM') then
! Either upper bound from the "CUSTOM" distribution itself or user-specified max diameter
               phase_max_dia = min(IC_PSD_MAX_DP(icv, m), find_max_dia(1, icv, m))
               max_radius = dmax1(max_radius,HALF*phase_max_dia)
            else
               phase_max_dia = d_p0(m)
            endif

            SELECT CASE(SOLIDS_MODEL(M))

            CASE('DEM', 'CGP')
               icpbb(:) = ONE

            CASE('SQP')
               ! icpbb(M,1) = 2.0 * SQP_a(M) / phase_max_dia
               ! icpbb(M,2) = 2.0 * SQP_b(M) / phase_max_dia
               ! icpbb(M,3) = 2.0 * SQP_c(M) / phase_max_dia
               ! The above leads to fpe. Don't use the spacing correction for now
               icpbb(:) = ONE

            CASE('GSP')
! Only use the particle bounding box (pbb) when the gsp is not rotated.
! In the future, we should recompute the pbb when the gsp is rotated.
! This can lead to initial overlap with hexa lattice. Only use pbb with cubic
! lattice
               gsp_rotated = .NOT.(gsp_q1(m)==one.and.gsp_q2(m)==zero.and.gsp_q3(m)==zero.and.gsp_q4(m)==zero)
               if(ic_gsp_random_q(1,1).or.gsp_rotated.or.TRIM(IC_DES_LATTICE(ICV,M))=='HEXA') then
                  icpbb(:) = one
               else
                  icpbb(:) = pbb(M,:)
               endif

            CASE DEFAULT

            END SELECT
! Use the CGP diameter if needed
            if(CGDEM.and.cgp_scaling_method(m)==2) phase_max_dia = cgp_d_p0(m)
! setting particle seed spacing grid to be slightly greater than
! the particle diameter. Spacing can vary in each direction
            lFAC         = ONE + IC_DES_SPACING(ICV,M)
            lFAC_X       = ONE + IC_DES_SPACING(ICV,M)*IC_DES_SPACE_FACTOR_X(ICV,M)
            lFAC_Y       = ONE + IC_DES_SPACING(ICV,M)*IC_DES_SPACE_FACTOR_Y(ICV,M)
            lFAC_Z       = ONE + IC_DES_SPACING(ICV,M)*IC_DES_SPACE_FACTOR_Z(ICV,M)
            ADJ_RADIUS   = 0.5d0*phase_max_dia*lFAC
            ADJ_DIA      = 2.0d0*ADJ_RADIUS
            ADJ_RADIUS_X = 0.5d0*phase_max_dia*lFAC_X*icpbb(1)
            ADJ_RADIUS_Y = 0.5d0*phase_max_dia*lFAC_Y*icpbb(2)
            ADJ_RADIUS_Z = 0.5d0*phase_max_dia*lFAC_Z*icpbb(3)
            ADJ_DIA_X    = 2.0d0*ADJ_RADIUS_X
            ADJ_DIA_Y    = 2.0d0*ADJ_RADIUS_Y
            ADJ_DIA_Z    = 2.0d0*ADJ_RADIUS_Z

! Compute seeding info based on lattice type
! Convert lattice keyword into integer for faster processing at the particle level (see triple nested loop below)
! Narrow z-domains with hexa lattice degenerate to triangular lattice

            NARROW_Z = .FALSE.

            SELECT CASE(TRIM(IC_DES_LATTICE(ICV,M)))

               CASE('CUBIC')
                  LATTICE = LATTICE_CUBIC
                  SEED_X = max(1,floor((IC_END(1)-IC_START(1)-ADJ_DIA_X)/ADJ_DIA_X))
! yinit instead of IC_START(2) to consider cases with multiple solids phases
! Note yInit is updated at the end after initializing each solids phase to next solids phase
                  !SEED_Y = max(1,floor((IC_END(2)-IC_START(2)-ADJ_DIA_Y)/ADJ_DIA_Y))
                  SEED_Y = max(1,floor((IC_END(2)-yInit-ADJ_DIA_Y)/ADJ_DIA_Y))
                  SEED_Z = max(1,floor((IC_END(3)-IC_START(3)-ADJ_DIA_Z)/ADJ_DIA_Z))
! JFD: Do not convert SEED_{X,Y,Z} back to spacing, use the adjusted
! diameter, so the spacing defined with advanced dem seeding options are
! honored.
!                  lDX = DOML(1)/dble(SEED_X)
                  !lDY = DOML(2)/dble(SEED_Y)
!                  lDY = (IC_END(2)-yInit)/dble(SEED_Y)
                  lDX = ADJ_DIA_X
                  lDY = ADJ_DIA_Y
                  IF(DO_K) THEN
!                     lDZ = DOML(3)/dble(SEED_Z)
                     lDZ = ADJ_DIA_Z
                  ELSE
                     lDZ = 0.0d0
                  ENDIF
                  rfac_base = ONE

               CASE('HEXA')
                  LATTICE = LATTICE_HEXA
                  IF(DO_K) THEN
                     SEED_X = (IC_END(1)-IC_START(1)-ADJ_RADIUS_X)/ADJ_RADIUS_X
                     SEED_X = floor(0.5d0*(SEED_X) - 1.5D0)

! yinit instead of IC_START(2) to consider cases with multiple solids phases
! Note yInit is updated at the end after initializing each solids phase to next solids phase
                     !SEED_Y = (IC_END(2)-IC_START(2)-ADJ_RADIUS_Y)/ADJ_RADIUS_Y
                     SEED_Y = (IC_END(2)-yInit-ADJ_RADIUS_Y)/ADJ_RADIUS_Y
                     SEED_Y = floor(1.5d0*(SEED_Y - 1.0D0)/dsqrt(6.0d0))

                     SEED_Z = (IC_END(3)-IC_START(3)-ADJ_RADIUS_Z)/ADJ_RADIUS_Z
                     SEED_Z = floor(1.0d0*(SEED_Z - 1.0D0)/dsqrt(3.0d0) - 1.0d0/3.0d0)
                     IF(SEED_Z<1) THEN
                        SEED_Z   = 1
                        NARROW_Z = .TRUE.
                        ZINIT    = 0.5D0*(IC_START(3)+IC_END(3))
                     ENDIF
                  ELSE
                     SEED_X = (IC_END(1)-IC_START(1)-ADJ_RADIUS_X)/ADJ_RADIUS_X
                     SEED_X = floor(0.5d0*(SEED_X) - 0.5D0)

! yinit instead of IC_START(2) to consider cases with multiple solids phases
! Note yInit is updated at the end after initializing each solids phase to next solids phase
                     !SEED_Y = (IC_END(2)-IC_START(2)-ADJ_RADIUS_X)/ADJ_RADIUS_X
                     SEED_Y = (IC_END(2)-yInit-ADJ_RADIUS_X)/ADJ_RADIUS_X
                     SEED_Y = floor(1.5d0*(SEED_Y - 1.0D0)/dsqrt(3.0d0))

                     SEED_Z = 1
                  ENDIF
                  rfac_base = HALF

            END SELECT

            SEED_X = MAX(SEED_X,1)
            SEED_Y = MAX(SEED_Y,1)
            SEED_Z = MAX(SEED_Z,1)

! Set velocity and position fluctuations if needed
            VEL_FLUCT = SET_VEL_FLUCT(ICV,M)
            XYZ_FLUCT = SET_XYZ_FLUCT(ICV,M)

! Go through the lattice vertically (JJ loop).
! Find position of next particle, and add it.
! The whole lattice will seed more than needed (rPARTS(M) particles)
! After each xz plane is seeded, count total number of particles
! and stop seeding if particle count is greater than rPARTS(M)
! Then each rank will remove particles until the exact count
! of rPARTS(M) is obtained.
            JJ_LP: DO JJ=1, SEED_Y
            KK_LP: DO KK=1, SEED_Z
            II_LP: DO II=1, SEED_X

               SELECT CASE(LATTICE)

                  CASE(LATTICE_CUBIC)
                     POS(1) = xINIT + (II-1)*lDX + HALF*lDX
                     POS(2) = YINIT + (JJ-1)*lDY + HALF*lDY
                     POS(3) = ZINIT + (KK-1)*lDZ + HALF*lDZ

                  CASE(LATTICE_HEXA)
                     IF(DO_K.AND..NOT.NARROW_Z) THEN
                        POS(1) = xINIT + ( 2*II + mod(JJ+KK,2) ) * ADJ_RADIUS_X
                        POS(2) = YINIT + ( 2*JJ/3.0  ) * dsqrt(6.0d0) * ADJ_RADIUS_Y
                        POS(3) = ZINIT + ( KK + mod(JJ,2)/3.0 ) * dsqrt(3.0d0) * ADJ_RADIUS_Z
                     ELSE
                        POS(1) = xINIT + ( 2*II + mod(JJ,2) ) * ADJ_RADIUS_X
                        POS(2) = YINIT + ( 2*JJ/2.0  ) * dsqrt(3.0d0) * ADJ_RADIUS_Y
                        POS(3) = ZINIT
                     ENDIF
               END SELECT

! HangZ: instead of adding 1, we should add number of spheres in that non-spherical particles
! @renjieke, in the implicit model, gsp works like normal dem for seeding
               IF (SOLIDS_MODEL(M) == 'GSP' .AND. GSP_EXPLICIT) THEN
                  pCOUNT(M) = pCOUNT(M) + 1
                  tCOUNT = tCOUNT + sphereNumber(M)
               ELSE
                  pCOUNT(M) = pCOUNT(M) + 1
                  tCOUNT = tCOUNT + 1
               ENDIF

! Add random position fluctuation if needed
! Not available for GSP
               IF(XYZ_FLUCT) THEN
                  rfac   = rfac_base * IC_DES_RAND(ICV,M) * IC_DES_SPACING(ICV,M) * phase_max_dia
                  rfac_x = rfac * IC_DES_SPACE_FACTOR_X(ICV,M)
                  rfac_y = rfac * IC_DES_SPACE_FACTOR_Y(ICV,M)
                  rfac_z = rfac * IC_DES_SPACE_FACTOR_Z(ICV,M)

                  POS(1) = POS(1) + (randXYZ(pCOUNT(M),1)-half) * rfac * IC_DES_RAND_FACTOR_X(ICV,M)
                  POS(2) = POS(2) + (randXYZ(pCOUNT(M),2)-half) * rfac * IC_DES_RAND_FACTOR_Y(ICV,M)
                  POS(3) = POS(3) + (randXYZ(pCOUNT(M),3)-half) * rfac * IC_DES_RAND_FACTOR_Z(ICV,M)
               ENDIF

! Nudge particles that fall right on top of a des grid line
               IF(compare(POS(1),dg_xstart) .OR. compare(POS(1),dg_xend))    &
                  POS(1) = POS(1) + SMALL_NUMBER
               IF(compare(POS(2),dg_ystart) .OR. compare(POS(2),dg_yend))    &
                  POS(2) = POS(2) + SMALL_NUMBER
               IF(DO_K) THEN
                  IF(compare(POS(3),dg_zstart) .OR. compare(POS(3),dg_zend)) &
                     POS(3) = POS(3) + SMALL_NUMBER
               ENDIF
! Nudge gluePos that fall right on top of a des grid line, which normally means it lies on boundary
               ! d_p0(M)
               ! IF(POS(1)<dg_xstart+0.5*d_p0(M) .or. POS(1)>dg_xend-0.5*d_p0(M)) CYCLE II_LP
               ! IF(POS(1)<dg_ystart+0.5*d_p0(M) .or. POS(1)>dg_yend-0.5*d_p0(M)) CYCLE JJ_LP
               ! IF(DO_K) THEN
               !    IF(POS(1)<dg_zstart+0.5*d_p0(M) .or. POS(1)>dg_zend-0.5*d_p0(M)) CYCLE KK_LP
               ! ENDIF

! HangZ: getting positions for each sphere ==> add a new subroutine called ADD_PARTICLE_GSP
! Add the particle at coordinates POS. This excludes particles that are not in a fluid cell
! or overlap with walls. All variables are set inside ADD_PARTICLE
               IF (SOLIDS_MODEL(M) == 'GSP') THEN
                  IF (GSP_EXPLICIT) THEN
                     CALL ADD_PARTICLE_GSP(POS,PSEEDED(M),IC_DES_CHECK_STL_OVERLAP(ICV,M),IC_GSP_RANDOM_Q(ICV,M),ADJ_RADIUS_X,ADJ_RADIUS_Y,ADJ_RADIUS_Z)
                  ELSE
                     CALL ADD_PARTICLE(POS,PSEEDED(M),IC_DES_CHECK_STL_OVERLAP(ICV,M),IC_GSP_RANDOM_Q(ICV,M))
                  ENDIF
               ELSE
                  CALL ADD_PARTICLE(POS,PSEEDED(M),IC_DES_CHECK_STL_OVERLAP(ICV,M),IC_SQP_RANDOM_Q(ICV,M))
               ENDIF
               ! Quit seeding along the xz plane if rank has already seeded more
               ! than the total number of particles of phase M
               IF(pSEEDED(M)>rPARTS(M)) EXIT KK_LP
            ENDDO II_LP
            ENDDO KK_LP
! Exit if all particles were seeded.
! The test is done at the JJ_LP loop level to minimize communication
! among processors when doing the global sum on SOLIDS_DATA
            GLOBAL_SOLIDS_DATA = SOLIDS_DATA(M)
            CALL GLOBAL_ALL_SUM(GLOBAL_SOLIDS_DATA)
            IF(GLOBAL_SOLIDS_DATA >= rPARTS(M)) THEN
! Get ready for the next phase. Will start just above the current xz plane
               YINIT = POS(2) + ADJ_RADIUS
               EXIT JJ_LP
            ENDIF

            ENDDO JJ_LP
! Since global sum was done at the JJ_LP loop level, the particle count will
! most likely exceed the desired particle count.
! It is adjusted below by removing particles belonging to each processors
! 1) Gather PIP from each PE onto head node
! 2) On head node, loop through PIP, and decrease count by one if PIP>0
! 3) Continue until sum of PIP matched the desire count
! 4) Scatter new PIP values to processors
            ALL_PART_COUNT_INIT(MyPE) = INT(SOLIDS_DATA(M))
            if(numPEs>1) then
              CALL allgather_1i (ALL_PART_COUNT_INIT(MyPE),ALL_PART_COUNT_ADJ,IERR)
            else
              all_part_count_adj = all_part_count_init
            end if
            IF (myPE == PE_IO) THEN
               NP_to_remove = SUM(ALL_PART_COUNT_ADJ) - rPARTS(M)
               NR_LP: DO NR = 1, NP_to_remove
                  DO iproc = 0,NumPEs-1
                     IF(ALL_PART_COUNT_ADJ(iproc)>0) ALL_PART_COUNT_ADJ(iproc) = ALL_PART_COUNT_ADJ(iproc) - 1
                     IF(SUM(ALL_PART_COUNT_ADJ)<=rPARTS(M) ) EXIT NR_LP
                  ENDDO
               ENDDO NR_LP
            ENDIF

            CALL BCAST(ALL_PART_COUNT_ADJ)

! Now each rank knows how many particles it needs to remove from its local count
! Reset PIP and set type to NONEXISTENT for particles that are removed
            SOLIDS_DATA(M) = SOLIDS_DATA(M) - FLOAT(ALL_PART_COUNT_INIT(MyPE)-ALL_PART_COUNT_ADJ(MyPE))
            OLD_PIP = PIP
            IF (SOLIDS_MODEL(M) == 'GSP') THEN
               IF(GSP_EXPLICIT) THEN
                  NGluedParticles = NGluedParticles - FLOAT(ALL_PART_COUNT_INIT(MyPE)-ALL_PART_COUNT_ADJ(MyPE))
                  NEW_PIP = PIP + (ALL_PART_COUNT_ADJ(MyPE)- ALL_PART_COUNT_INIT(MyPE))*sphereNumber(M)
               ELSE
                  ! @renjieke there is no need to use extra NSphereGSP to check the number of components
                  ! NSphereGSP = NSphereGSP + (ALL_PART_COUNT_ADJ(MyPE)- ALL_PART_COUNT_INIT(MyPE))*sphereNumber(M)
                  NEW_PIP = PIP - ALL_PART_COUNT_INIT(MyPE) + ALL_PART_COUNT_ADJ(MyPE)
               ENDIF
            ELSE
               NEW_PIP = PIP - ALL_PART_COUNT_INIT(MyPE) + ALL_PART_COUNT_ADJ(MyPE)
            ENDIF
            IF(NEW_PIP<OLD_PIP) THEN
               DO LL=NEW_PIP+1,OLD_PIP
                  CALL DELETE_PARTICLE(LL)
               ENDDO
            ENDIF
      ENDDO   ! Loop over phase M

! Collect the data
      CALL GLOBAL_ALL_SUM(SOLIDS_DATA)

! Verify that the IC region volume is not zero.
      IF(SOLIDS_DATA(0) <= 0.0d0) THEN
         WRITE(ERR_MSG,1001) ICV, SOLIDS_DATA(0)
         CALL LOG_ERROR()
      ENDIF

1001 FORMAT('Error 1001: Invalid IC region volume: IC=',I3,' VOL=', ES15.4)

      IF (SOLIDS_MODEL(M) == 'GSP') THEN
         WRITE(ERR_MSG,2020) ICV
      ELSE
         WRITE(ERR_MSG,2000) ICV
      ENDIF

      CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)

      DO M=MMAX+1, MMAX+DES_MMAX
         IF(SOLIDS_DATA(M) < SMALL_NUMBER) CYCLE
         total_volume = 0.D0

         DO II = PIP_PREVIOUS + 1, MAX_PIP
           IF(.NOT.IS_NORMAL(II)) CYCLE
           IF(PIJK(II,5) == M) THEN
             TOTAL_VOLUME = TOTAL_VOLUME + PVOL(II)
           END IF
         END DO

         CALL GLOBAL_ALL_SUM(total_volume)

         ECHO_EPS = total_volume/SOLIDS_DATA(0)

         ECHO_SM = total_volume*RO_S0(M)
         ECHO_NP = INT(SOLIDS_DATA(M))

         IF(IC_EP_S(ICV,M)>ZERO.AND.IC_EP_S(ICV,M)/=UNDEFINED) THEN
            ECHO_INPUT_VAR = '   EP_S    '
            WRITE(ECHO_INPUT_VALUE,"(1x,ES9.2,1x)") IC_EP_S(ICV,M)

         ELSEIF(IC_DES_SM(ICV,M)>ZERO.AND.IC_DES_SM(ICV,M)/=UNDEFINED) THEN
            ECHO_INPUT_VAR = ' SOL. MASS '
            WRITE(ECHO_INPUT_VALUE,"(1x,ES9.2,1x)") IC_DES_SM(ICV,M)

         ELSEIF(IC_DES_NP(ICV,M)>0.AND.IC_DES_NP(ICV,M)/=UNDEFINED_I) THEN
            ECHO_INPUT_VAR = '  PART NUM '
            WRITE(ECHO_INPUT_VALUE,"(1x,I9,1x)") IC_DES_NP(ICV,M)

         ELSE
            ECHO_INPUT_VAR   = '  UNKNOWN  '
            ECHO_INPUT_VALUE = '  UNKNOWN  '

         ENDIF

         IF (SOLIDS_MODEL(M) == 'GSP') THEN
            ECHO_NP_SPHERE = INT(SOLIDS_DATA(M)*sphereNumber(M))
            WRITE(ERR_MSG,2030) M, ECHO_INPUT_VAR, ECHO_INPUT_VALUE, &
                                ECHO_NP, ECHO_NP_SPHERE, ECHO_EPS, ECHO_SM
         ELSE
            WRITE(ERR_MSG,2010) M, ECHO_INPUT_VAR, ECHO_INPUT_VALUE, &
                                ECHO_NP, ECHO_EPS, ECHO_SM
         ENDIF

         CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)
      ENDDO
      write(err_msg,*) ''
      CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.TRUE.)

! Verif the seeding was successful
      DO M=MMAX+1, MMAX+DES_MMAX
         IF(SOLIDS_DATA(M) < SMALL_NUMBER) CYCLE
            IF(SOLIDS_DATA(M) == rPARTS(M)) THEN
               WRITE(ERR_MSG,"(2/,'Initial condition',I3,', Phase',I3,': Successful seeding')") ICV, M
               CALL LOG_INFO()
            ELSE
               WRITE(ERR_MSG,"(2/,'Initial condition',I3,', Phase',I3,': Unsuccessful seeding'  &
               /,'Please adjust seeding volume fraction, or advanced seeding options (lattice, spacing).'  &
               /,'To bypass this error (not recommended), set IC_DES_ABORT_ON_FAILED_SEEDING(ICV,M)=.FALSE.' &
               )") ICV, M
!#  FIXMEif IC_DES_ABORT is set, do not print the "bypass" text
                 IF(IC_DES_ABORT_ON_FAILED_SEEDING(ICV,M)) THEN
                    CALL LOG_ERROR()
                 ELSE
                    CALL LOG_INFO()
                 ENDIF
            ENDIF

      ENDDO
!      IF(allocated(kf_data)) deallocate(kf_data)
      IF(allocated(tmp_array)) deallocate(tmp_array)
!      IF(allocated(prob_for_each_phase)) deallocate(prob_for_each_phase)

      IF(allocated(randVEL)) deallocate(randVEL)

! Initialize arrays (in case they were grown)
      MIN_OVERLAP(:) = ZERO
      MAX_OVERLAP(:) = ZERO
      MEAN_OVERLAP(:) = ZERO
      COORDINATION(:) = ZERO



      RETURN

! 2000 FORMAT(/2x,'|',43('-'),'|',/2x,'| IC region: ',I3,28x,'|',/2x,   &
!         '|',43('-'),'|',/2x,'| Phase | Number of |    EPs    |    EP',&
!         's    |',/2x,'|   ID  | Particles | Specified |   Actual  |', &
!         /2x,'|-------|',3(11('-'),'|'))
!
! 2010 FORMAT(2x,'|  ',I3,'  |',1x,I9,1x,'|',2(1x,ES9.2,1x,'|'),/2x,    &
!         '|-------|',3(11('-'),'|'))
 2000 FORMAT(/2x,'|',67('-'),'|',              &
             /2x,'| IC region: ',I3,52x,'|',   &
             /2x,'|',67('-'),'|',              &
             /2x,'| Phase |   Input   |   Input   | Number of |    EPs    |   Solid   |', &
             /2x,'|   ID  | specified |   value   | particles |           |   mass    |', &
             /2x,'|-------|',5(11('-'),'|'))
 2010 FORMAT(2x,'|  ',I3,'  |',A11,'|',A11,'|',1x,I9,1x,'|',2(1x,ES9.2,1x,'|'),/2x,    &
         '|-------|',5(11('-'),'|'))

 2020 FORMAT(/2x,'|',85('-'),'|',              &
             /2x,'| IC region: ',I3,70x,'|',   &
             /2x,'|',85('-'),'|',              &
             /2x,'| Phase |   Input   |   Input   |    Number of    | Number of |    EPs    |   Solid   |', &
             /2x,'|   ID  | specified |   value   | glued particles |  spheres  |           |   mass    |', &
             /2x,'|-------|',2(11('-'),'|'), 17('-'),'|', 3(11('-'),'|'))

 2030 FORMAT(2x,'|  ',I3,'  |',A11,'|',A11,'|',1x,I15,1x,'|',1x,I9,1x,'|',2(1x,ES9.2,1x,'|'),/2x,    &
         '|-------|',2(11('-'),'|'), 17('-'),'|', 3(11('-'),'|'))


      CONTAINS

!......................................................................!
! Function: SET_VEL_FLUCT                                              !
! Purpose: Set the flag for random velocity fluctuations. If needed    !
! the random velocities are calculated.                                !
!......................................................................!
      LOGICAL FUNCTION SET_VEL_FLUCT(lICV, lM)
      INTEGER, INTENT(IN) :: lICV, lM
      DOUBLE PRECISION :: VEL_SIG
      SET_VEL_FLUCT=(IC_Theta_M(lICV,lM) /= 0.0d0)
      IF(SET_VEL_FLUCT) THEN
         if(allocated(randVEL)) deallocate(randVEL)
         allocate(randVEL(100+int(rPARTS(lM)),3))
         VEL_SIG = sqrt(IC_Theta_M(lICV,lM))
         CALL NOR_RNO(randVEL(:,1), IC_U_s(lICV,lM),VEL_SIG)
         CALL NOR_RNO(randVEL(:,2), IC_V_s(lICV,lM),VEL_SIG)
         IF(DO_K) CALL NOR_RNO(randVEL(:,3),IC_W_s(lICV,lM),VEL_SIG)
      ELSE
         if(allocated(randVEL)) deallocate(randVEL)
         allocate(randVEL(1,3))
      ENDIF

      END FUNCTION SET_VEL_FLUCT


!......................................................................!
! Function: SET_XYZ_FLUCT                                              !
! Purpose: Set the flag for random position fluctuations. If needed    !
! the random positions  are calculated.                                !
!......................................................................!
      LOGICAL FUNCTION SET_XYZ_FLUCT(lICV, lM)
      use randomno
      INTEGER, INTENT(IN) :: lICV, lM
      DOUBLE PRECISION :: XYZ_SIG
      SET_XYZ_FLUCT=(IC_DES_RAND(lICV,lM) /= 0.0d0)
      IF(SET_XYZ_FLUCT) THEN
         if(allocated(randXYZ)) deallocate(randXYZ)
         allocate(randXYZ(10*int(rPARTS(lM)),3))
         XYZ_SIG = IC_DES_RAND(ICV,M) * 0.5D0 * D_P0(lM)
         CALL UNI_RNO(randXYZ(:,1)) !, 0.0D0,XYZ_SIG)
         CALL UNI_RNO(randXYZ(:,2)) !, 0.0D0,XYZ_SIG)
         IF(DO_K) CALL UNI_RNO(randXYZ(:,3)) !, 0.0D0,XYZ_SIG)
      ENDIF

      END FUNCTION SET_XYZ_FLUCT

!......................................................................!
! Function: SET_VEL_FLUCT_mod                                          !
! Purpose: Faster/optimized version of SET_VEL_FLUCT                   !
!......................................................................!
      LOGICAL FUNCTION SET_VEL_FLUCT_mod(lICV, lM)
      INTEGER, INTENT(IN) :: lICV, lM
      DOUBLE PRECISION :: VEL_SIG
      SET_VEL_FLUCT_mod=(IC_Theta_M(lICV,lM) /= 0.0d0)

      END FUNCTION SET_VEL_FLUCT_mod


!......................................................................!
! subroutine: ADD_PARTICLE                                             !
! Purpose: Add particle at position POS with additional checks         !
!......................................................................!

      SUBROUTINE ADD_PARTICLE(POS,PSEEDED,CHECK_STL_OVERLAP,RANDOM_QUATERNION)

      use cutcell, only: SMALL_CELL_AT
      USE des_rxns, only: DES_X_s
      use des_thermo, only: DES_T_s
      use discretelement
      use ic
      use param1, only: undefined, zero
      use physprop, only: C_PS0
      use physprop, only: NMAX
      use run, only: ENERGY_EQ, SPECIES_EQ
      use sq_rotation_mod
      use sq_properties_mod
      use gsp_math_mod, only: gsp_quat_to_exyz

      IMPLICIT NONE
! Particle position and velocity
      DOUBLE PRECISION :: POS(3), VEL(3)
! Number of seeded particle, per solids phase
      INTEGER :: PSEEDED
! Flag to check if particle overlaps with STL
      LOGICAL :: CHECK_STL_OVERLAP
! Flag to randomize initial sqp or gsp orientation
      LOGICAL :: RANDOM_QUATERNION
! Species loop index
      INTEGER :: NN

      DOUBLE PRECISION :: SVOLUME
      DOUBLE PRECISION :: qnorm

! Keep only particles that belong to this process.
      IF(.NOT.dg_is_ON_myPE_OWNs(IofPOS(POS(1)), &
         JofPOS(POS(2)),KofPOS(POS(3)))) RETURN

! Bin the parcel to the fluid grid.
      CALL PIC_SEARCH(I, POS(1), XE, DIMENSION_I, IMIN2, IMAX2)
      CALL PIC_SEARCH(J, POS(2), YN, DIMENSION_J, JMIN2, JMAX2)
      K=1
      IF(DO_K) CALL PIC_SEARCH(K, POS(3), ZT, DIMENSION_K, KMIN2, KMAX2)

! Skip cells that return invalid IJKs.
      IF(DEAD_CELL_AT(I,J,K)) RETURN

! Skip cells that are not part of the local fluid domain.
! Small cut cells are allowed here, if particles intersect the STL geometry,
! they can be removed below
      IJK = FUNIJK(I,J,K)
      IF(.NOT.FLUID_AT(IJK).AND..NOT.SMALL_CELL_AT(IJK)) RETURN

      IF(CARTESIAN_GRID.AND.CHECK_STL_OVERLAP) THEN
         CALL CHECK_IF_PARTICLE_OVERLAPS_STL(POS, I, J, K, SKIP)
         IF(SKIP) RETURN
      ENDIF

      PIP = PIP + 1
      PSEEDED = PSEEDED + 1
      CALL PARTICLE_GROW(PIP)
      MAX_PIP = max(PIP,MAX_PIP)

      CALL SET_NORMAL(PIP)

      IF(VEL_FLUCT.and.pip<=size(randvel,1)) THEN  ! avoid getting out of bounds
         VEL(1) = randVEL(PIP,1)
         VEL(2) = randVEL(PIP,2)
         VEL(3) = randVEL(PIP,3)
      ELSE
         VEL(1) = IC_U_s(ICV,M)
         VEL(2) = IC_V_s(ICV,M)
         VEL(3) = IC_W_s(ICV,M)
      ENDIF
      IF(NO_K) VEL(3) = 0.0d0


!----------------------
      if(IC_PSD_TYPE(ICV,M) .eq. 'MONO') then
        DES_RADIUS(PIP) = D_P0(M)*HALF
        ! llu, d_p0 was already scaled to parcel size
        IF(CGP_SCALING_METHOD(M)==2) then
           DES_RADIUS(PIP) = DES_RADIUS(PIP)/CGP_STAT_WT(M)**(1.0d0/3.0d0)
        endif
      else
        DES_RADIUS(PIP) = psd_radius(1, ICV, M)
      end if

!----------------------


      DES_POS_NEW(PIP,:) = POS(:)
      DES_VEL_NEW(PIP,:) = VEL(:)
      OMEGA_NEW(PIP,:) = 0.0d0

      RO_SOL(PIP) =  RO_S0(M)

      PIJK(PIP,1) = I
      PIJK(PIP,2) = J
      PIJK(PIP,3) = K
      PIJK(PIP,4) = IJK
      PIJK(PIP,5) = M

      IF(DO_OLD) THEN
         DES_VEL_OLD(PIP,:) = DES_VEL_NEW(PIP,:)
         DES_POS_OLD(PIP,:) = DES_POS_NEW(PIP,:)
         OMEGA_OLD(PIP,:)   = ZERO
      ENDIF

      IF(CGDEM) THEN
         IF(CGP_SCALING_METHOD(M)==1) THEN
            DES_CGP_STW(PIP) = CGP_STAT_WT(M)
            DES_CGP_RPR(PIP) = DES_RADIUS(PIP)/DES_CGP_STW(PIP)**(1.0d0/3.0d0)
         ELSEIF(CGP_SCALING_METHOD(M)==2) THEN
            DES_CGP_STW(PIP) = (HALF*CGP_D_P0(M)/DES_RADIUS(PIP))**3
            DES_CGP_RPR(PIP) = DES_RADIUS(PIP)
            DES_RADIUS(PIP) = HALF*CGP_D_P0(M)
         ENDIF
      ENDIF

! Update volume
      PVOL(PIP) = (4.0D0/3.0D0)*PI*DES_RADIUS(PIP)**3

      IF(SuperDEM) THEN
         super_mn(PIP,1)=SQP_m(M)
         super_mn(PIP,2)=SQP_n(M)
         super_r(PIP,1)=SQP_a(M) * 2.0D0 * DES_RADIUS(PIP) / D_P0(M) ! scale a,b,c for polydisperse systems
         super_r(PIP,2)=SQP_b(M) * 2.0D0 * DES_RADIUS(PIP) / D_P0(M)
         super_r(PIP,3)=SQP_c(M) * 2.0D0 * DES_RADIUS(PIP) / D_P0(M)

         IF(RANDOM_QUATERNION) THEN
            qnorm = 0
            DO WHILE (qnorm .gt. 1 .or. qnorm .lt. 1.0d-8)
                ! If point is outside unit 4-sphere, reject it and try again. This assures uniform distribution.
                CALL RANDOM_NUMBER(super_q(PIP,:))            ! random number from [0,1]
                super_q(PIP,:) = 2.0D0*super_q(PIP,:) - 1.0D0 ! shift range to [-1,1]
                qnorm = super_q(PIP,1)**2 + super_q(PIP,2)**2 + super_q(PIP,3)**2 + super_q(PIP,4)**2
            END DO
            super_q(PIP,:) = super_q(PIP,:) / sqrt(qnorm)
         ELSE
            super_q(PIP,1)=SQP_q1(M)
            super_q(PIP,2)=SQP_q2(M)
            super_q(PIP,3)=SQP_q3(M)
            super_q(PIP,4)=SQP_q4(M)
         ENDIF

         CALL SQ_VOLUME (super_r(PIP,1:3),super_mn(PIP,1),super_mn(PIP,2),svolume)
         PVOL(PIP)=svolume

         IF(PARTICLE_ORIENTATION) THEN
            ! Convert quaternion to orientation vector
            CALL QROTATE(super_q(PIP,:),ORIENTATION(PIP,:),INIT_ORIENTATION,2)
         ENDIF
      ELSE IF(GluedSphereDEM) THEN
         IF (IC_PSD_TYPE(ICV,M) .eq. 'MONO') THEN
            PVOL(PIP) = gsp_vol(m)
         ELSE
         ! Particle size distribution is given based on bounding diameter of glued particle.
         ! To get the actual volume of gsp, mean_volume needs to be converted based on the ratio
         ! between actual volume of gsp (gsp_vol) and volume from bounding diameter.
            PVOL(PIP) = PVOL(PIP)*gsp_vol(m)/((pi*d_p0(m)**3.0D0)/6.0d0)
         ENDIF
         PMASS(PIP) = PVOL(PIP)*RO_SOL(PIP)
         IF(RANDOM_QUATERNION) THEN
            qnorm = 0
            DO WHILE (qnorm .gt. 1 .or. qnorm .lt. 1.0d-8)
            ! If point is outside unit 4-sphere, reject it and try again. This assures uniform distribution.
               CALL RANDOM_NUMBER(glueQuat(PIP,:)) ! random number from [0,1]
               glueQuat(PIP,:) = 2.0D0*glueQuat(PIP,:) - 1.0D0 ! shift range to [-1,1]
               qnorm = glueQuat(PIP,1)**2 + glueQuat(PIP,2)**2 + &
                       glueQuat(PIP,3)**2 + glueQuat(PIP,4)**2
            ENDDO
            glueQuat(PIP,:) = glueQuat(PIP,:) / sqrt(qnorm)
         ELSE
            glueQuat(PIP,1)=GSP_q1(M)
            glueQuat(PIP,2)=GSP_q2(M)
            glueQuat(PIP,3)=GSP_q3(M)
            glueQuat(PIP,4)=GSP_q4(M)
         ENDIF
         IF(PARTICLE_ORIENTATION) THEN
            ! Convert quaternion to orientation vector
            CALL QROTATE(glueQuat(PIP,:),ORIENTATION(PIP,:),INIT_ORIENTATION,2)
         ENDIF

         glueAngMom(PIP, :) = 0.0d0
         glueSurface(PIP) = 0.0d0
         call gsp_quat_to_exyz(glueQuat(PIP, :), glueEX(PIP,:), glueEY(PIP,:), glueEZ(PIP,:))

      ENDIF

! Set the initial particle temperature.
      IF(ENERGY_EQ) THEN
         DES_T_s(PIP) = IC_T_s(ICV,M)
      ENDIF

! Set the initial species composition.
      IF((ENERGY_EQ .AND. C_Ps0(M) == UNDEFINED) .OR.         &
         SPECIES_EQ(M)) THEN
         DES_X_s(PIP,:) = ZERO
         DO NN = 1, NMAX(M)
            DES_X_s(PIP,NN) = IC_X_s(ICV,M,NN)
         ENDDO
      ENDIF

      SOLIDS_DATA(M) = SOLIDS_DATA(M) + 1.0

! Residence time
      RESIDENCE_TIME(PIP) = ZERO

      RETURN
      END SUBROUTINE ADD_PARTICLE

!......................................................................!
! subroutine: ADD_PARTICLE_GSP                                         !
! Author: Hang Zhou                                     Date: Dec 2023 !
! Revised: Renjie Ke                                    Date: Mar 2024 !
! Purpose: Add sphere particle at position POS with additional checks  !
!          Also include surface area, neighbor info array construction !
!          Revised: all variables now are relative values based on     !
!                   bounding diameter                                  !
!......................................................................!
      SUBROUTINE ADD_PARTICLE_GSP(POS,PSEEDED,CHECK_STL_OVERLAP,RANDOM_QUATERNION,ADJ_RADIUS_X,ADJ_RADIUS_Y,ADJ_RADIUS_Z)

      use cutcell, only: SMALL_CELL_AT
      USE des_rxns, only: DES_X_s
      use des_thermo, only: DES_T_s
      use discretelement
      use ic
      use param1, only: undefined, zero
      use physprop, only: C_PS0
      use physprop, only: NMAX
      use run, only: ENERGY_EQ, SPECIES_EQ
      use GSP_MATH_MOD, only: gsp_quat_to_exyz

      IMPLICIT NONE
! Particle position and velocity
      DOUBLE PRECISION :: POS(3), VEL(3)
! Distance vector, current sphere points to neighbor sphere
      DOUBLE PRECISION :: distvec(3)
! Number of seeded particle, per solids phase
      INTEGER :: PSEEDED
! Flag to check if particle overlaps with STL
      LOGICAL :: CHECK_STL_OVERLAP
! Flag to randomize initial glued sphere orientation
      LOGICAL :: RANDOM_QUATERNION
! Species loop index
      INTEGER :: NN, nbID
! index of spheres
      INTEGER :: lcurpar, PIP_START, Jid
      DOUBLE PRECISION :: SVOLUME, eff_dist
      DOUBLE PRECISION :: qnorm
! Radius of non-spherical particle if particle size distribution is used
      DOUBLE PRECISION :: psd_radius_tmp
! Dist2, for inertia calculation
      DOUBLE PRECISION :: dist2
      DOUBLE PRECISION, INTENT(IN) :: ADJ_RADIUS_X
      DOUBLE PRECISION, INTENT(IN) :: ADJ_RADIUS_Y
      DOUBLE PRECISION, INTENT(IN) :: ADJ_RADIUS_Z
! ---------------------------------------------------------------------->

! HangZ: set the process based on the positions of the whole gluded particles
! Keep only particles that belong to this process.
      IF(.NOT.dg_is_ON_myPE_OWNs(IofPOS(POS(1)), &
         JofPOS(POS(2)),KofPOS(POS(3)))) RETURN
! HangZ: set the fluid grid based on the positions of the whole gluded particles
! Bin the parcel to the fluid grid.
      CALL PIC_SEARCH(I, POS(1), XE, DIMENSION_I, IMIN2, IMAX2)
      CALL PIC_SEARCH(J, POS(2), YN, DIMENSION_J, JMIN2, JMAX2)
      K=1
      IF(DO_K) CALL PIC_SEARCH(K, POS(3), ZT, DIMENSION_K, KMIN2, KMAX2)

! Skip cells that return invalid IJKs.
      IF(DEAD_CELL_AT(I,J,K)) RETURN

! Skip cells that are not part of the local fluid domain.
! Small cut cells are allowed here, if particles intersect the STL geometry,
! they can be removed below
      IJK = FUNIJK(I,J,K)
      IF(.NOT.FLUID_AT(IJK).AND..NOT.SMALL_CELL_AT(IJK)) RETURN

      IF(CARTESIAN_GRID.AND.CHECK_STL_OVERLAP) THEN
         CALL CHECK_IF_PARTICLE_OVERLAPS_STL(POS, I, J, K, SKIP)
         IF(SKIP) RETURN
      ENDIF

      PSEEDED = PSEEDED + 1
      NGluedParticles = NGluedParticles+1
      CALL PARTICLE_GROW_GluedSphereDEM(NGluedParticles)
      GlueInertia(NGluedParticles,:) = ZERO
      GlueMass(NGluedParticles) = ZERO
      GlueVolume(NGluedParticles) = ZERO

      IF(RANDOM_QUATERNION) THEN
         qnorm = 0
         DO WHILE (qnorm .gt. 1 .or. qnorm .lt. 1.0d-8)
             ! If point is outside unit 4-sphere, reject it and try again. This assures uniform distribution.
             CALL RANDOM_NUMBER(glueQuat(NGluedParticles,:))! random number from [0,1]
             glueQuat(NGluedParticles,:) = 2.0D0*glueQuat(NGluedParticles,:) - 1.0D0 ! shift range to [-1,1]
             qnorm = glueQuat(NGluedParticles,1)**2 + glueQuat(NGluedParticles,2)**2 + &
                     glueQuat(NGluedParticles,3)**2 + glueQuat(NGluedParticles,4)**2
         ENDDO
         glueQuat(NGluedParticles,:) = glueQuat(NGluedParticles,:) / sqrt(qnorm)
      ELSE
         glueQuat(NGluedParticles,1)=GSP_q1(M)
         glueQuat(NGluedParticles,2)=GSP_q2(M)
         glueQuat(NGluedParticles,3)=GSP_q3(M)
         glueQuat(NGluedParticles,4)=GSP_q4(M)
      ENDIF

      CALL gsp_quat_to_exyz(glueQuat(NGluedParticles,:), glueEX(NGluedParticles,:), glueEY(NGluedParticles,:), glueEZ(NGluedParticles,:))
      IF(VEL_FLUCT) THEN  ! avoid getting out of bounds
         IF(PSEEDED<=size(randvel,1)) THEN
            VEL(1) = randVEL(PSEEDED,1)
            VEL(2) = randVEL(PSEEDED,2)
            VEL(3) = randVEL(PSEEDED,3)
         ENDIF
      ELSE
         VEL(1) = IC_U_s(ICV,M)
         VEL(2) = IC_V_s(ICV,M)
         VEL(3) = IC_W_s(ICV,M)
      ENDIF
      IF(NO_K) VEL(3) = 0.0d0
!HangZ: getting information of each sphere, including position, diameter, velocity, density, species, and temperature
!HangZ: for now, the density, species and temperatures of sphere in one non-spherical particle are same
      ! monodisperse versus polydisverse
      if(IC_PSD_TYPE(ICV,M) .eq. 'MONO') then
         glueBounding(NGluedParticles) = D_P0(M)
      else
         psd_radius_tmp = psd_radius(1, ICV, M)
         glueBounding(NGluedParticles) = psd_radius_tmp*2.0d0
      endif
      ! this subroutine will be called under 1:M loop
      ! PIP_START is the first global id of component sphere in a glued sphere
      PIP_START = PIP + 1
      DO lcurpar =1, sphereNumber(M)!HangZ: for each sphere
         ! grow the array related to component sphere size
         ! REAL_GROW2_reverse change the first dimension
         ! REAL_GROW2 change the second dimension
         PIP = PIP + 1
         CALL PARTICLE_GROW(PIP)
         MAX_PIP = max(PIP,MAX_PIP)
         CALL SET_NORMAL(PIP)

!HangZ: getting the diameter and position of each sphere
         DES_RADIUS(PIP) = glueBounding(NGluedParticles)*HALF*sphereRelDia(start_index_sphere(M)+lcurpar-1)
         sc2gpc_vec(PIP,:) = glueBounding(NGluedParticles)*sphereRelPos(start_index_sphere(M)+lcurpar-1,:)
         DES_POS_NEW(PIP,:) = POS(:)+sc2gpc_vec(PIP,1)*glueEX(NGluedParticles,:) &
                                    +sc2gpc_vec(PIP,2)*glueEY(NGluedParticles,:) &
                                    +sc2gpc_vec(PIP,3)*glueEZ(NGluedParticles,:)
! read heat conduction related parameters
         gp_sa(PIP) = tmp_gp_sa(start_index_sphere(M)+lcurpar-1) * glueBounding(NGluedParticles) * glueBounding(NGluedParticles)
! already reset to relative id in read_particle_input
         gp_neighsid(PIP,:) = tmp_gp_neighsid(start_index_sphere(M)+lcurpar-1,:)
         where (gp_neighsid(PIP,:) /= -1) gp_neighsid(PIP,:) = gp_neighsid(PIP,:) + PIP_START
         gp_neighsa(PIP,:) = tmp_gp_neighsa(start_index_sphere(M)+lcurpar-1,:) * glueBounding(NGluedParticles) * glueBounding(NGluedParticles)

         DES_VEL_NEW(PIP,:) = VEL(:)
         OMEGA_NEW(PIP,:) = 0.0d0

         RO_SOL(PIP) =  RO_S0(M)

         PIJK(PIP,1) = I
         PIJK(PIP,2) = J
         PIJK(PIP,3) = K
         PIJK(PIP,4) = IJK
         PIJK(PIP,5) = M
         gid_list(PIP) = NGluedParticles

         IF(DO_OLD) THEN
            DES_VEL_OLD(PIP,:) = DES_VEL_NEW(PIP,:)
            DES_POS_OLD(PIP,:) = DES_POS_NEW(PIP,:)
            OMEGA_OLD(PIP,:)   = ZERO
         ENDIF

! Update volume
         PVOL(PIP) = (4.0D0/3.0D0)*PI*DES_RADIUS(PIP)**3
! Update glued diameter
         glueVolume(NGluedParticles) = glueVolume(NGluedParticles) + PVOL(PIP)
! Update Mass
         PMASS(PIP) = PVOL(PIP) * RO_SOL(PIP)
         glueMass(NGluedParticles) = glueMass(NGluedParticles) + PMASS(PIP)
! Update glueInertia
! for a solid sphere , if sphereNumber(phase) == 1 then gluedInertia == 2/5*m*r^2
          if (sphereNumber(M) .eq. 1) then
             glueInertia(NGluedParticles, :) = 2.0d0/5.0d0*PMASS(PIP)*DES_RADIUS(PIP)**2
          else
             dist2 = sc2gpc_vec(PIP,2)*sc2gpc_vec(PIP,2) + sc2gpc_vec(PIP,3)*sc2gpc_vec(PIP,3)
             glueInertia(NGluedParticles,1) = glueInertia(NGluedParticles,1) + PMASS(PIP) * dist2 + (2.0d0/5.0d0)*PMASS(PIP)*DES_RADIUS(PIP)**2

             dist2 = sc2gpc_vec(PIP,1)*sc2gpc_vec(PIP,1) + sc2gpc_vec(PIP,3)*sc2gpc_vec(PIP,3)
             glueInertia(NGluedParticles,2) = glueInertia(NGluedParticles,2) + PMASS(PIP) * dist2 + (2.0d0/5.0d0)*PMASS(PIP)*DES_RADIUS(PIP)**2

             dist2 = sc2gpc_vec(PIP,1)*sc2gpc_vec(PIP,1) + sc2gpc_vec(PIP,2)*sc2gpc_vec(PIP,2)
             glueInertia(NGluedParticles,3) = glueInertia(NGluedParticles,3) + PMASS(PIP) * dist2 + (2.0d0/5.0d0)*PMASS(PIP)*DES_RADIUS(PIP)**2
          endif
! Set the initial particle temperature.
         IF(ENERGY_EQ) THEN
            DES_T_s(PIP) = IC_T_s(ICV,M)
         ENDIF

! Set the initial species composition.
         IF((ENERGY_EQ .AND. C_Ps0(M) == UNDEFINED) .OR. &
            SPECIES_EQ(M)) THEN
            DES_X_s(PIP,:) = ZERO

! HangZ: confirm where IC_X_s is defined
            DO NN = 1, NMAX(M)
               DES_X_s(PIP,NN) = IC_X_s(ICV,M,NN)
            ENDDO
         ENDIF

! Residence time
         RESIDENCE_TIME(PIP) = ZERO
      ENDDO ! HangZ: end loop for all spheres in one non-spherical particle

      gluePos(NGluedParticles,:) = POS(:)
      glueVel(NGluedParticles,:) = VEL(:)
      glueDiameter(NGluedParticles) = 2.0D0*((3.0D0/4.0D0)*glueVolume(NGluedParticles)/PI)**(1.0D0/3.0D0)
      CALL SET_NORMAL_GSP(NGluedParticles)
      SOLIDS_DATA(M) = SOLIDS_DATA(M) + 1.0

      RETURN
      END SUBROUTINE ADD_PARTICLE_GSP

      END SUBROUTINE GENERATE_PARTICLE_CONFIG_DEM


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_PARTICLE_XYZ_IN_CSV                              C
!  Purpose: Write particles coordinates in a csv file                  C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 16-Feb-18  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_PARTICLE_XYZ_IN_CSV

!      use des_mpi
      use cdist
      use mpi_comm_des
      use mpi_utility
! Number of DEM solids phases.
      use discretelement, only: DES_MMAX
! Number of TFM solid phases
      use physprop, only:  MMAX
! Find available file unit
      use funits, only: newunit
! Error manager
      use error_manager

      IMPLICIT NONE


      INTEGER :: LOCAL_CNT, GLOBAL_CNT


      DOUBLE PRECISION, ALLOCATABLE :: ltemp_array(:,:)  ! local
      DOUBLE PRECISION, ALLOCATABLE :: gtemp_array(:,:)  ! global

      LOGICAL, ALLOCATABLE :: WRITE_TO_CSV(:)

      INTEGER :: M,MM
      INTEGER :: LB, UB
      INTEGER :: PC, LC1, LC2
! Variables related to gather
      integer :: lgathercnts(0:numpes-1), lproc
! local unit
      INTEGER :: lunit
! local filename
       character(255) :: lfilename
! IO Status:
       INTEGER :: IOS
! Flag to indicate if file exists.
       LOGICAL :: lEXISTS


! Loop through all particles and kee a list of particles belonging to a VTK region

! Since the data is appended (i.e., written after all tags), the
! offset, in number of bytes must be specified.  The offset includes
! the size of the data for each field, plus the size of the integer
! that stores the number of bytes.  this is why the offset of a field
! equals the offset of the previous field plus sizeof(int) plus the
! number of bytes of the field.

! Next, the actual data is written for the geometry (PASS=WRITE_DATA)
! The DATA is converted to single precision to save memory.


!      CALL INIT_ERR_MSG("WRITE_PARTICLE_XYZ_IN_CSV")

      IF(ALLOCATED(WRITE_TO_CSV)) DEALLOCATE(WRITE_TO_CSV)
      ALLOCATE(WRITE_TO_CSV(MAX_PIP))

      DO MM=MMAX+1,MMAX+DES_MMAX

         WRITE_TO_CSV(:) = .FALSE.

         IF (.NOT.BDIST_IO) THEN
! The number of particles

            LOCAL_CNT = 0
            GLOBAL_CNT = 10
            PC = 1
            DO LC1 = 1, MAX_PIP
               IF(PC > PIP) EXIT
               IF(IS_NONEXISTENT(LC1)) CYCLE
               PC = PC+1
               IF(IS_ANY_GHOST(LC1)) CYCLE

               M = PIJK(LC1,5)
               IF(M==MM) THEN
                  WRITE_TO_CSV(LC1) = .TRUE.
                  LOCAL_CNT = LOCAL_CNT + 1
               ENDIF
            ENDDO ! particle loop

! Calculate the total number of particles system-wide.
            call global_sum(LOCAL_CNT, GLOBAL_CNT)

! 1788 ! No need to set the send/reccv when using distributed IO
! 1789       IF (BDIST_IO) RETURN
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

            LB = LBOUND(DES_POS_NEW,2)
            UB = UBOUND(DES_POS_NEW,2)

            ALLOCATE (dProcBuf(LOCAL_CNT) )
            ALLOCATE (ltemp_array((UB-LB)+1,LOCAL_CNT))
            IF(MYPE==PE_IO) THEN
               ALLOCATE (gtemp_array((UB-LB)+1,GLOBAL_CNT))
               ALLOCATE (dRootBuf(GLOBAL_CNT))
            ELSE
               ALLOCATE (gtemp_array((UB-LB)+1,10))
               ALLOCATE (dRootBuf(10))
            ENDIF

! Pack particle coordinates in a temporary local array
            PC = 0
            DO LC1 = 1, MAX_PIP
               IF(WRITE_TO_CSV(LC1)) THEN
                  PC =PC + 1
                  DO LC2=LB, UB
                     ltemp_array(LC2,PC) = DES_POS_NEW(LC1,LC2)
                  ENDDO
               ENDIF
               IF(PC==LOCAL_CNT) EXIT
            ENDDO

! For each coordinate (x,y, and z), gather the local list to global temporary array
            DO LC1 = LB, UB
               dprocbuf(1:LOCAL_CNT)=ltemp_array(LC1,1:LOCAL_CNT)
               CALL desmpi_gatherv(ptype=2)
               IF(MYPE == PE_IO) gtemp_array(LC1,:) = drootbuf(:)
            ENDDO

! Write the list of coordinates
            IF(myPE == PE_IO) THEN
               lFILENAME = ''
               WRITE(lFILENAME,'("PART_XYZ_AT_TSTOP_M=",I2.2,".CSV")') MM
               INQUIRE(FILE=lFILENAME, EXIST=lEXISTS)
               IF(LEXISTS) THEN
                  WRITE(ERR_MSG, 1000) TRIM(lFILENAME)
                  CALL LOG_ERROR()
               ENDIF
               lUNIT = newunit()
               OPEN(UNIT=lUNIT, FILE=lFILENAME, FORM="FORMATTED")

1000 FORMAT('Warning 1100: ',A,' already exists and will be overwritten.')

               WRITE(lUNIT,*)  "x, y, z"
               DO LC1=1, GLOBAL_CNT
                  WRITE(lUNIT,1100)  gtemp_array(LB:UB,LC1)
               ENDDO
               CLOSE(lUNIT)
            ENDIF

1100 FORMAT(E14.8,' , ',E14.8,' , ',E14.8)

            deallocate (dProcBuf, dRootBuf, ltemp_array,gtemp_array)




         ELSEIF(BDIST_IO.AND.LOCAL_CNT>0) THEN

            IF(LOCAL_CNT==0) RETURN

         ENDIF

      ENDDO ! MM Solids phase loop

!      CALL FINL_ERR_MSG
      RETURN

      END SUBROUTINE WRITE_PARTICLE_XYZ_IN_CSV





!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG (OLD, monodisperse version)    !
!  Authors: Rahul Garg                               Date: 21-Mar-2014 !
!                                                                      !
!  Purpose: Generate particle configuration for DEM solids for each IC !
!           region. Now using the particle linked lists for initial    !
!           build                                                      !
!           This routine will ultimately supersede the older routine   !
!           that has not been deleted yet                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GENERATE_PARTICLE_CONFIG_DEM_0(ICV)

! Global Variables:
!---------------------------------------------------------------------//
! particle radius and density
      use discretelement, only: DES_RADIUS, RO_Sol
! particle position new and old
      use discretelement, only: DES_POS_NEW, DES_POS_OLD
! particle velocity new and old
      use discretelement, only: DES_VEL_NEW, DES_VEL_OLD
! Number of particles in the system (current)
      use discretelement, only: PIP
! Number of DEM solids phases.
      use discretelement, only: DES_MMAX
! Flag to use _OLD variables
      use discretelement, only: DO_OLD
! Angular velocity
      use discretelement, only: OMEGA_OLD, OMEGA_NEW, PIJK
! solid phase diameters and densities.
      use physprop, only: D_p0, RO_s0, MMAX
! IC Region solids volume fraction.
      use ic, only: IC_EP_s, IC_RO_s

! Constant: 3.14159...
      use constant, only: PI
! min and max physical coordinates of IC regions in each direction
      use ic, only: IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t
! initially specified velocity field and granular temperature
      use ic, only: IC_U_s, IC_V_s, IC_W_s, IC_Theta_M
! Parameter for detecting unspecified values, zero, and one
      use param1, only: ZERO, Half
! Parameter for small numbers
      use param1, only: SMALL_NUMBER

! to access random number generator subroutines
      use randomno
      use mpi_utility

      use desgrid, only: dg_xstart, dg_ystart, dg_zstart
      use desgrid, only: dg_xend, dg_yend, dg_zend

! direction wise spans of the domain and grid spacing in each direction
      use geometry, only: xlength, ylength, zlength

      use cutcell, only: CARTESIAN_GRID
      use stl_functions_des, only: CHECK_IF_PARTICLE_OVERLAPS_STL
      use run, only: solids_model
      use des_allocate, only: PARTICLE_GROW

      use desgrid, only: IofPOS, JofPOS, KofPOS
      use desgrid, only: dg_is_ON_myPE_OWNs
      use toleranc, only: compare

      use discretelement, only: max_pip, xe, yn, zt
      use functions
      use param, only: dim_m
      use param, only: dimension_i, dimension_j, dimension_k
      use param1, only: undefined
      use ic, only: IC_PIC_CONST_STATWT
      use discretelement, only: cgdem, des_cgp_stw, des_cgp_rpr
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: ICV

! Local variables
!---------------------------------------------------------------------//
! Starting positions in the axial directions
      DOUBLE PRECISION :: xINIT, yINIT, zINIT
! Fractor used to scale particle diameter
      DOUBLE PRECISION :: lFAC
! Particle position and velocity
      DOUBLE PRECISION :: POS(3), VEL(3)
! Number of particles in the lattice
      INTEGER :: SEED_X, SEED_Y, SEED_Z
! Loop indices phase/fluid cell
      INTEGER :: M, MM, I, J, K, IJK
! Loop indices for seeding
      INTEGER :: II, JJ, KK
! Start and end bound for IC region.
      DOUBLE PRECISION :: IC_START(3), IC_END(3)
! Volume and lengths of the IC Region
      DOUBLE PRECISION :: DOM_VOL, DOML(3)
! Flag to skip the particle
      LOGICAL :: SKIP
! Diameter adjusted for space padding
      DOUBLE PRECISION :: ADJ_DIA
! Number of particles calculated from volume fracton
      INTEGER :: rPARTS(DIM_M), tPARTS
! Spacing between particles.
      DOUBLE PRECISION :: lDEL, lDX, lDY, lDZ
! Number of seeded particles
      INTEGER :: pCOUNT(DIM_M), tCOUNT

      DOUBLE PRECISION :: SOLIDS_DATA(0:DIM_M), lRO(DIM_M)

      LOGICAL :: VEL_FLUCT
      DOUBLE PRECISION, ALLOCATABLE :: randVEL(:,:)

      DOUBLE PRECISION :: BC_MAX_RADIUS

!......................................................................!

      WRITE(ERR_MSG,"(2/,'Generating initial particle configuration:')")
      CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)

      SOLIDS_DATA = ZERO
      CALL GET_IC_VOLUME(ICV, SOLIDS_DATA(0))

! setting particle seed spacing grid to be slightly greater than
! the maximum particle diameter. seed at ~particle radii
      lFAC = 1.05D0

! Setup local arrays with IC region bounds.
      IC_START(1)=IC_X_W(ICV);   IC_END(1)=IC_X_E(ICV)
      IC_START(2)=IC_Y_S(ICV);   IC_END(2)=IC_Y_N(ICV)
      IC_START(3)=IC_Z_B(ICV);   IC_END(3)=IC_Z_T(ICV)

      DOML = IC_END-IC_START
      IF(NO_K) DOML(3)=DZ(1)

      ! Create a local array of solids densities so we don't have to
      ! reevaluate the function for every particle. Doing this for each
      ! phase is a bit overkill as DEM says that each IC region can only
      ! have on solids.
      DO M=MMAX+1,MMAX+DES_MMAX
         IF(SOLIDS_MODEL(M) == 'DEM' .and. IC_EP_S(ICV,M) > ZERO) THEN
            lRO(M) = IC_RO_S(ICV, M)
         ELSE
            lRO(M) = -UNDEFINED
         END IF
      END DO

! Volume of the IC region
      DOM_VOL = DOML(1)*DOML(2)*DOML(3)

      rPARTS=0
      BC_MAX_RADIUS = ZERO
      DO M=MMAX+1,MMAX+DES_MMAX
         IF(SOLIDS_MODEL(M) == 'DEM') THEN
            ! Number of particles for phase M
            rPARTS(M) = &
                 floor((6.0d0*IC_EP_S(ICV,M)*DOM_VOL)/(PI*(D_P0(M)**3)))
! Negative values of rParts indicate overflow. Treat this as a fatal error.
      IF(rPARTS(M)<0) THEN
         WRITE(ERR_MSG,100) ICV, M, DOM_VOL,D_P0(M),(6.0d0*IC_EP_S(ICV,M)*DOM_VOL)/(PI*(D_P0(M)**3))
         CALL LOG_ERROR()
      ENDIF

100 FORMAT('Error 100: Overflow in IC region, IC=',I3,' , M=',I3, &
      /'IC Volume=',ES15.4,' , Particle Diameter=', ES15.4, ' , Number of particles = ',ES15.4, &
      /'This error is usually triggered by attempting to initialize with a very large number of particles.' &
      /'Please verify your settings (Particle size and region extents).', &
      /'This is a fatal error. MFiX will exit now.')
            IF(IC_EP_S(ICV,M) > ZERO .AND. D_P0(M) > BC_MAX_RADIUS)&
                 BC_MAX_RADIUS = HALF * D_P0(M)
         ENDIF
      ENDDO

! Total number of particles in this IC region.
      tPARTS = sum(rPARTS)
      IF(tPARTS == 0) RETURN

      ADJ_DIA = 2.0d0*BC_MAX_RADIUS*lFAC

! Generic filling
      SEED_X = max(1,floor((IC_END(1)-IC_START(1)-ADJ_DIA)/ADJ_DIA))
      SEED_Y = max(1,floor((IC_END(2)-IC_START(2)-ADJ_DIA)/ADJ_DIA))
      SEED_Z = max(1,floor((IC_END(3)-IC_START(3)-ADJ_DIA)/ADJ_DIA))

      lDX = DOML(1)/dble(SEED_X)
      lDY = DOML(2)/dble(SEED_Y)
      IF(DO_K) THEN
         lDZ = DOML(3)/dble(SEED_Z)
      ELSE
         lDZ = 0.0d0
      ENDIF

      xINIT = IC_START(1)+HALF*lDX
      yINIT = IC_START(2)+HALF*lDY
      zINIT = IC_START(3)+HALF*lDZ

      M=1
      pCOUNT = 0
      tCOUNT = 0

      VEL_FLUCT = SET_VEL_FLUCT(ICV,M)

      JJ_LP: DO JJ=1, SEED_Y
         POS(2) = YINIT + (JJ-1)*lDY
         IF(compare(POS(2),dg_ystart) .OR. compare(POS(2),dg_yend))    &
            POS(2) = POS(2) + SMALL_NUMBER

      KK_LP: DO KK=1, SEED_Z
         POS(3) = ZINIT + (KK-1)*lDZ
         IF(DO_K) THEN
            IF(compare(POS(3),dg_zstart) .OR. compare(POS(3),dg_zend)) &
               POS(3) = POS(3) + SMALL_NUMBER
         ENDIF

      II_LP: DO II=1, SEED_X
         POS(1) = xINIT + (II-1)*lDX
         IF(compare(POS(1),dg_xstart) .OR. compare(POS(1),dg_xend))    &
            POS(1) = POS(1) + SMALL_NUMBER

! Exit if all particles were seeded.
         IF(tCOUNT > int(tPARTS)) THEN
            EXIT JJ_LP
! Find the next phase that needs to be seeded
         ELSEIF((pCOUNT(M) > int(rPARTS(M))).OR.(rParts(M)==ZERO)) THEN
            MM_LP: DO MM=M+1,MMAX+DES_MMAX
               IF(rPARTS(MM) > 0.0) THEN
                  M=MM
                  EXIT MM_LP
               ENDIF
            ENDDO MM_LP
            IF(M > MMAX+DES_MMAX) EXIT JJ_LP
            VEL_FLUCT = SET_VEL_FLUCT(ICV,M)
         ENDIF

         pCOUNT(M) = pCOUNT(M) + 1
         tCOUNT = tCOUNT + 1

! Keep only particles that belong to this process.
         IF(.NOT.dg_is_ON_myPE_OWNs(IofPOS(POS(1)), &
            JofPOS(POS(2)),KofPOS(POS(3)))) CYCLE

! Bin the parcel to the fluid grid.
         K=1
         IF(DO_K) CALL PIC_SEARCH(K, POS(3), ZT, DIMENSION_K, KMIN2, KMAX2)
         CALL PIC_SEARCH(J, POS(2), YN, DIMENSION_J, JMIN2, JMAX2)
         CALL PIC_SEARCH(I, POS(1), XE, DIMENSION_I, IMIN2, IMAX2)

! Skip cells that return invalid IJKs.
         IF(DEAD_CELL_AT(I,J,K)) CYCLE

! Skip cells that are not part of the local fluid domain.
         IJK = FUNIJK(I,J,K)
         IF(.NOT.FLUID_AT(IJK)) CYCLE

         IF(CARTESIAN_GRID) THEN
            CALL CHECK_IF_PARTICLE_OVERLAPS_STL(POS, I, J, K, SKIP)
            IF(SKIP) CYCLE
         ENDIF

         PIP = PIP + 1
         CALL PARTICLE_GROW(PIP)
         MAX_PIP = max(PIP,MAX_PIP)

         CALL SET_NORMAL(PIP)

         IF(VEL_FLUCT) THEN
            VEL(1) = randVEL(pCOUNT(M),1)
            VEL(2) = randVEL(pCOUNT(M),2)
            VEL(3) = randVEL(pCOUNT(M),3)
         ELSE
            VEL(1) = IC_U_s(ICV,M)
            VEL(2) = IC_V_s(ICV,M)
            VEL(3) = IC_W_s(ICV,M)
         ENDIF
         IF(NO_K) VEL(3) = 0.0d0


         DES_POS_NEW(PIP,:) = POS(:)
         DES_VEL_NEW(PIP,:) = VEL(:)
         OMEGA_NEW(PIP,:) = 0.0d0

         DES_RADIUS(PIP) = D_P0(M)*HALF
         RO_SOL(PIP) =  lRO(M)

         PIJK(PIP,1) = I
         PIJK(PIP,2) = J
         PIJK(PIP,3) = K
         PIJK(PIP,4) = IJK
         PIJK(PIP,5) = M

         IF(DO_OLD) THEN
            DES_VEL_OLD(PIP,:) = DES_VEL_NEW(PIP,:)
            DES_POS_OLD(PIP,:) = DES_POS_NEW(PIP,:)
            OMEGA_OLD(PIP,:) = ZERO
         ENDIF

         IF(CGDEM) THEN
            DES_CGP_STW(PIP) = IC_PIC_CONST_STATWT(ICV,M)
            DES_CGP_RPR(PIP) = DES_RADIUS(PIP)/DES_CGP_STW(PIP)**(1.0d0/3.0d0)
            print*,'PIP,real radiu=',PIP,DES_CGP_RPR(PIP),DES_RADIUS(PIP),DES_CGP_STW(PIP)
      if(DES_CGP_STW(PIP) < 0.99) write(*,*) 'Error! Check IC_PIC_CONST_STATWT in initial domain: ', ICV, 'Phase: ', M
         ENDIF


         SOLIDS_DATA(M) = SOLIDS_DATA(M) + 1.0

      ENDDO II_LP
      ENDDO KK_LP
      ENDDO JJ_LP

! Collect the data
      CALL GLOBAL_ALL_SUM(SOLIDS_DATA)

! Verify that the IC region volume is not zero.
      IF(SOLIDS_DATA(0) <= 0.0d0) THEN
         WRITE(ERR_MSG,1000) ICV, SOLIDS_DATA(0)
         CALL LOG_ERROR()
      ENDIF

1000 FORMAT('Error 1000: Invalid IC region volume: IC=',I3,' VOL=', ES15.4)

      WRITE(ERR_MSG,2000) ICV
      CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)

      DO M=MMAX+1, MMAX+DES_MMAX
         IF(SOLIDS_DATA(M) < SMALL_NUMBER) CYCLE
         WRITE(ERR_MSG,2010) M, int(SOLIDS_DATA(M)), IC_EP_S(ICV,M),   &
            (dble(SOLIDS_DATA(M))*(Pi/6.0d0)*D_P0(M)**3)/SOLIDS_DATA(0)
         CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)
      ENDDO

      IF(allocated(randVEL)) deallocate(randVEL)

      RETURN

 2000 FORMAT(/2x,'|',43('-'),'|',/2x,'| IC region: ',I3,28x,'|',/2x,   &
         '|',43('-'),'|',/2x,'| Phase | Number of |    EPs    |    EP',&
         's    |',/2x,'|   ID  | Particles | Specified |   Actual  |', &
         /2x,'|-------|',3(11('-'),'|'))

 2010 FORMAT(2x,'|  ',I3,'  |',1x,I9,1x,'|',2(1x,ES9.2,1x,'|'),/2x,    &
         '|-------|',3(11('-'),'|'))


      CONTAINS

!......................................................................!
! Function: SET_VEL_FLUCT                                              !
! Purpose: Set the flag for random velocity fluctuations. If needed    !
! the random velocities are calculated.                                !
!......................................................................!
      LOGICAL FUNCTION SET_VEL_FLUCT(lICV, lM)
      INTEGER, INTENT(IN) :: lICV, lM
      DOUBLE PRECISION :: VEL_SIG
      SET_VEL_FLUCT=(IC_Theta_M(lICV,lM) /= 0.0d0)
      IF(SET_VEL_FLUCT) THEN
         if(allocated(randVEL)) deallocate(randVEL)
         allocate(randVEL(100+int(rPARTS(lM)),3))
         VEL_SIG = sqrt(IC_Theta_M(lICV,lM))
         CALL NOR_RNO(randVEL(:,1), IC_U_s(lICV,lM),VEL_SIG)
         CALL NOR_RNO(randVEL(:,2), IC_V_s(lICV,lM),VEL_SIG)
         IF(DO_K) CALL NOR_RNO(randVEL(:,3),IC_W_s(lICV,lM),VEL_SIG)
      ENDIF

      END FUNCTION SET_VEL_FLUCT

      END SUBROUTINE GENERATE_PARTICLE_CONFIG_DEM_0




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GENERATE_PARTICLE_CONFIG_MMPPIC                         !
!  Author: Rahul Garg                                 Date: 3-May-2011 !
!                                                                      !
!  Purpose: Generates particle position distribution for MPPIC.        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GENERATE_PARTICLE_CONFIG_MPPIC(ICV)

! Global variables
!---------------------------------------------------------------------//
! Number of DES solids phases.
      use discretelement, only: DES_MMAX
! IC Region bulk density (RO_s * EP_s)
      use ic, only: IC_ROP_s
      use ic, only: IC_EP_s
      use ic, only: IC_PIC_CONST_STATWT

      use param1, only: ZERO
      use param, only: DIM_M
      use physprop, only: mmax

!(new) Polydispersity information
      use ic, only: IC_PSD_TYPE
      use ic, only: IC_PSD_MEAN_DP, IC_PSD_STDEV
      use ic, only: IC_PSD_MIN_DP
      use ic, only: IC_PSD_MAX_DP

! The accumulated number of particles in each IJK.
      use mpi_utility, only: GLOBAL_ALL_SUM

      use run, only: solids_model
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: ICV

! Local variables
!---------------------------------------------------------------------//
! Generic loop counters
      INTEGER :: M
! Solids data in IC Region by phase:
      double precision :: solids_data(2*dim_m), ic_vol

      integer :: parcels
      double precision :: stat_wt, particles, ic_eps, act_eps

!......................................................................!

      SOLIDS_DATA = 0.0d0

      WRITE(ERR_MSG,"(2/,'Generating initial parcel configuration:')")
      CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)

      CALL GET_IC_VOLUME(ICV, IC_VOL)

      ! Collect the data
      CALL GLOBAL_ALL_SUM(IC_VOL)

      ! Verify that the IC region volume is not zero.
      IF(IC_VOL <= 0.0d0) THEN
         WRITE(ERR_MSG,1000) ICV, IC_VOL
         CALL LOG_ERROR()
      ENDIF

1000  FORMAT('Error 1000: Invalid IC region volume: IC=',I3,' VOL=', ES15.4)

      WRITE(ERR_MSG,2000) ICV
      CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_WARNING, HEADER=.FALSE., FOOTER=.FALSE.)

! Set up the individual solids phases.
      DO M=MMAX+1, DES_MMAX+MMAX
         IF(SOLIDS_MODEL(M) == 'PIC') THEN
            IF(IC_ROP_S(ICV,M) == ZERO) CYCLE
!     (NEW)
            IF(IC_PSD_TYPE(ICV,M).eq.'MONO') THEN
!              Seed parcels with constant dp. Using old routines. No change.
               CALL GPC_MPPIC_CONST_STATWT(ICV, M, SOLIDS_DATA((2*M-1):(2*M)))
            ELSEIF(IC_PSD_TYPE(ICV,M).eq.'NORMAL'.OR. &
                   IC_PSD_TYPE(ICV,M).eq.'LOG_NORMAL'.OR. &
                   IC_PSD_TYPE(ICV,M).eq.'CUSTOM') THEN
!              Seed parcels with poly dp.  New routine
               CALL POLYPIC_GPC_MPPIC_CONST_STATWT(ICV, M, &
                                                   SOLIDS_DATA((2*M-1):(2*M)))
            END IF

         ENDIF
      ENDDO

      ! Collect the data
      CALL GLOBAL_ALL_SUM(solids_data)

! Report solids information for the IC region.
      DO M=MMAX+1, DES_MMAX+MMAX
         IF(SOLIDS_MODEL(M) == 'PIC') THEN

            stat_wt   = ic_pic_const_statwt(icv,m)
            ic_eps    = ic_ep_s(icv,m)
            parcels   = int(SOLIDS_DATA(2*M-1))
            particles = parcels*stat_wt
            act_eps   = solids_data(2*m)/ic_vol

            WRITE(ERR_MSG,2010) M, parcels, particles, &
                 stat_wt, ic_eps, act_eps
            CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF
      ENDDO

      RETURN

 2000 FORMAT(/2x,'|',67('-'),'|',/2x,'| IC region: ',I3,52x,'|',/2x,   &
         '|',67('-'),'|',/2x,'| Phase | Num. Comp | Num. Real ',       &
         '|Statistical|    EPs    |    EPs    |',/2x,'|   ID  |  ',    &
         'Parcels  | Particles |   Weight  | Specified |   Actual  |', &
         /2x,'|-------|',5(11('-'),'|'))

2010  FORMAT(2x,'|  ',I3,'  |',1x,I9,1x,'|',es10.1,1x,'|',3(1x,ES9.2,1x,'|'),/2x,&
         '|-------|',5(11('-'),'|'))

      END SUBROUTINE GENERATE_PARTICLE_CONFIG_MPPIC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GENERATE_PARTICLE_CONFIG_MMPPIC                         !
!  Author: Rahul Garg                                 Date: 3-May-2011 !
!                                                                      !
!  Purpose: generates particle position distribution for MPPIC         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GPC_MPPIC_CONST_STATWT(ICV, M, sDATA)

! Global variables
!---------------------------------------------------------------------//
! Constant: 3.14159...
      use constant, only: PI
! Cut_cell identifier array
      use discretelement, only: XE, YN, ZT

      use des_allocate, only: PARTICLE_GROW
      use discretelement, only: PIJK, PIP, MAX_PIP
      use discretelement, only:  DES_VEL_NEW, DES_POS_NEW, DES_RADIUS, RO_SOL
      use discretelement, only:  RESIDENCE_TIME

! IC Region bulk density (RO_s * EP_s)
      use ic, only: IC_EP_s, IC_RO_s
      use ic, only: IC_I_w, IC_I_e, IC_J_s, IC_J_n, IC_K_b, IC_K_t

! initially specified velocity field and granular temperature
      use ic, only: IC_U_s, IC_V_s, IC_W_s, IC_Theta_M
      use ic, only: IC_PIC_CONST_STATWT

      use mfix_pic, only: des_stat_wt
      use mpi_utility

      use param1, only: ZERO, HALF

! solid phase diameters and densities.
      use physprop, only: D_p0

      use stl_functions_des, only: picSTLoverlap

      use desgrid, only: dg_funijk
      use desgrid, only: dg_istart2, dg_jstart2, dg_kstart2
      use desgrid, only: dg_iend2,   dg_jend2,   dg_kend2
      use desgrid, only: iofpos,     jofpos,     kofpos

      use randomno
      use functions
      use des_allocate, only: PARTICLE_GROW

      IMPLICIT NONE

! Dummy Arguments
!----------------------------------------------------------------------//
! Index of IC region and solids phase
      INTEGER, INTENT(IN) :: ICV, M
! Data about solids in the IC region.
      DOUBLE PRECISION, INTENT(OUT) :: sDATA(2)

! Local variables
!----------------------------------------------------------------------//

! Number of real and comp. particles in a cell.
      DOUBLE PRECISION ::  rPARTS
      INTEGER :: maxPARTS
      DOUBLE PRECISION :: DOML(3), IC_START(3)
! Parcel position with random
      DOUBLE PRECISION :: POS(3)
! Solids density for phase in IC region
      DOUBLE PRECISION :: lRO
! Average velocity and standard deivation
      DOUBLE PRECISION :: IC_VEL(3), VEL_SIG
! Arrasy for assigning random position and velocities
      DOUBLE PRECISION, ALLOCATABLE :: randVEL(:,:)
      DOUBLE PRECISION :: RAND(3)
! Statistical weights
      DOUBLE PRECISION :: STAT_WT
! Volume of a parcel and total solids volume
      DOUBLE PRECISION :: sVOL, sVOL_TOT, EFF_RAD, remainder
! Counter for seeded parcels.
      INTEGER :: SEEDED
! Generic loop indices and loop counters
      INTEGER :: I, J, K, IJK, LC, LC_MAX
      INTEGER :: dg_I, dg_J, dg_K, dg_IJK, skipped
!......................................................................!

      maxPARTS=25
      allocate(randVEL(maxPARTS,3))

      IC_VEL(1) = IC_U_s(ICV,M)
      IC_VEL(2) = IC_V_s(ICV,M)
      IC_VEL(3) = merge(IC_W_s(ICV,M),0.0d0,DO_K)

      VEL_SIG = sqrt(IC_Theta_M(ICV,M))

! Volume occupied by one particle
      sVOL = (Pi/6.0d0)*(D_P0(M)**3.d0)

      SEEDED = 0
      sVOL_TOT = 0.0d0

      remainder = 0.0d0

      STAT_WT = IC_PIC_CONST_STATWT(ICV,M)

      lRO   = IC_RO_S(ICV, M)

      DO K = IC_K_B(ICV), IC_K_T(ICV)
      DO J = IC_J_S(ICV), IC_J_N(ICV)
      DO I = IC_I_W(ICV), IC_I_E(ICV)

         IF(.not.IS_ON_myPE_wobnd(I,J,K)) cycle
         IF(DEAD_CELL_AT(I,J,K)) cycle

         IJK = FUNIJK(I,J,K)
         IF(.not.FLUID_AT(IJK)) cycle

         rPARTS = IC_EP_s(ICV,M)*VOL(IJK)/sVOL + remainder

! Seed parcels with a constant statistical weight
         LC_MAX = floor(rPARTS/STAT_WT)

         remainder = max(rPARTS - dble(LC_MAX)*STAT_WT, 0.0d0)

         IF(LC_MAX == 0) cycle

! Increase particle buffer
         IF(LC_MAX > maxPARTS) THEN
            maxPARTS = 2*LC_MAX
            if(allocated(randVEL)) deallocate(randVEL)
            allocate(randVEL(maxPARTS,3))
         ENDIF

         DO LC=1, merge(2,3,NO_K)
            IF(VEL_SIG > ZERO) THEN
               CALL NOR_RNO(randVEL(1:LC_MAX,LC), IC_VEL(LC), VEL_SIG)
            ELSE
               randVEL(1:LC_MAX,LC) = IC_VEL(LC)
            ENDIF
         ENDDO
         IF(NO_K) randVEL(1:LC_MAX,3) = 0.0d0

         IC_START(1) = XE(I-1)
         IC_START(2) = YN(J-1)
         IC_START(3) = ZERO;  IF(DO_K) IC_START(3) = ZT(K-1)

         DOML(1) = DX(I)
         DOML(2) = DY(J)
         DOML(3) = ZERO;  IF(DO_K) DOML(3) = DZ(K)

         EFF_RAD = D_P0(M)*HALF*(STAT_WT)**(1.0/3.0)

         PIC_LP: DO LC=1,LC_MAX

            CALL RANDOM_NUMBER(RAND)
            POS(:) = IC_START + DOML*RAND

            dg_i = min(dg_iend2,max(dg_istart2, iofpos(pos(1))))
            dg_j = min(dg_jend2,max(dg_jstart2, jofpos(pos(2))))
            dg_k = min(dg_kend2,max(dg_kstart2, kofpos(pos(3))))
            dg_ijk = dg_funijk(dg_i,dg_j,dg_k)

            skipped = 0
            do while(picSTLoverlap(dg_ijk, pos, eff_rad))

               skipped = skipped + 1
               call random_number(rand)
               pos(:) = ic_start + doml*rand

               dg_i = min(dg_iend2,max(dg_istart2, iofpos(pos(1))))
               dg_j = min(dg_jend2,max(dg_jstart2, jofpos(pos(2))))
               dg_k = min(dg_kend2,max(dg_kstart2, kofpos(pos(3))))
               dg_ijk = dg_funijk(dg_i,dg_j,dg_k)

               if(skipped >= 1000 ) cycle pic_lp
            enddo


            PIP = PIP + 1
            CALL PARTICLE_GROW(PIP)
            MAX_PIP = max(PIP,MAX_PIP)

            DES_POS_NEW(PIP,:) = POS(:)
            DES_VEL_NEW(PIP,:) = 1.0e-4*rand!randVEL(LC,:)

            DES_RADIUS(PIP) = D_P0(M)*HALF
            RO_SOL(PIP) =  lRO

            PIJK(PIP,1) = I
            PIJK(PIP,2) = J
            PIJK(PIP,3) = K
            PIJK(PIP,4) = IJK
            PIJK(PIP,5) = M

            DES_STAT_WT(PIP) = STAT_WT
            sVOL_TOT = sVOL_TOT + sVOL*STAT_WT

            RESIDENCE_TIME(PIP) = ZERO

            CALL SET_NORMAL(PIP)

            SEEDED = SEEDED + 1

         ENDDO PIC_LP

      ENDDO
      ENDDO
      ENDDO

      IF(allocated(randVEL)) deallocate(randVEL)

      sDATA(1) = dble(SEEDED)
      sDATA(2) = sVOL_TOT

      RETURN
      END SUBROUTINE GPC_MPPIC_CONST_STATWT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GET_IC_VOLUME                                           !
!  Author: J.Musser                                 Date: 26-Aug-2015  !
!                                                                      !
!  Purpose: Calculate the actual volume of the IC region.              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GET_IC_VOLUME(ICV, IC_VOL)

! IC region index bounds
      use ic, only: IC_I_w, IC_I_e
      use ic, only: IC_J_s, IC_J_n
      use ic, only: IC_K_b, IC_K_t

! Volume of computational cells.
      use geometry, only: VOL

      use functions
      use compar, only: dead_cell_at

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Index of IC region
      INTEGER, INTENT(IN) :: ICV
! Total calculated volume of IC region
      DOUBLE PRECISION, INTENT(OUT) :: IC_VOL

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: I, J, K, IJK
!......................................................................!


      IC_VOL = 0.0d0
      DO K = IC_K_B(ICV), IC_K_T(ICV)
      DO J = IC_J_S(ICV), IC_J_N(ICV)
      DO I = IC_I_W(ICV), IC_I_E(ICV)

         IF(.NOT.IS_ON_MYPE_WOBND(I,J,K)) CYCLE
         IF(DEAD_CELL_AT(I,J,K)) CYCLE

         IJK = FUNIJK(I,J,K)
         IF(FLUID_AT(IJK)) IC_VOL = IC_VOL + VOL(IJK)

      ENDDO
      ENDDO
      ENDDO

      RETURN

   END SUBROUTINE GET_IC_VOLUME

! NEW
! ---------------------------------------------------------------------//
! Subroutine: GENERATE_PARTICLE_CONFIG_MPPIC  PolyVersion
! Author:  M. Clarke
! October 2021
! Purpose:  Incorporate polydispersity during particle generation
! Note: Reusing as much of original non-poly version as possible
!---------------------------------------------------------------------//
      SUBROUTINE POLYPIC_GPC_MPPIC_CONST_STATWT(ICV,M,sDATA)

! Global Variables
!--------------------------------------------------------------------//
! Constants
      use constant, only: PI
      use param1, only: ZERO, HALF, ONE
! Managing cut cell identifiers
      use discretelement, only: XE, YN, ZT
! Use discrete parameters
      use des_allocate, only: PARTICLE_GROW
      use discretelement
! IC Region information
      use ic, only: IC_EP_s, IC_RO_s
      use ic, only: IC_I_w, IC_I_e, IC_J_s, IC_J_n, IC_K_b, IC_K_t
      use ic, only: IC_U_s, IC_V_s, IC_W_s, IC_Theta_M
      use ic, only: IC_PIC_CONST_STATWT
! PIC information
      use mfix_pic, only: des_stat_wt, dim_pic_mi_bin, pic_numbin
! MPI
      use mpi_utility
! Physical properties
      use physprop, only: D_p0
! STL files
      use stl_functions_des, only: picSTLoverlap
! functions for des grid
      use desgrid, only: dg_funijk
      use desgrid, only: dg_istart2, dg_jstart2, dg_kstart2
      use desgrid, only: dg_iend2, dg_jend2, dg_kend2
      use desgrid, only: iofpos, jofpos, kofpos
! polydisperse function information
      use ic, only: IC_PSD_TYPE
      use ic, only: IC_PSD_MEAN_DP
      use ic, only: IC_PSD_STDEV
      use ic, only: IC_PSD_MIN_DP
      use ic, only: IC_PSD_MAX_DP

      use des_psd   !for custom distributions
      use randomno
      use functions

      IMPLICIT NONE

!     index of IC region and solids phase
      INTEGER, INTENT(IN) :: ICV, M
!     solids data in the IC region
      DOUBLE PRECISION, INTENT(OUT) :: sDATA(2)

! Local Variables
!-------------------------------------------------------------------//
! Number of particles in a cell
      DOUBLE PRECISION :: rPARTS
      INTEGER ::  maxPARTS
      DOUBLE PRECISION :: DOML(3), IC_START(3)
! Parcel position by direction
      DOUBLE PRECISION :: POS(3)
! Solids density or phase in IC region
      DOUBLE PRECISION :: localRO
! Average velocity and standard deviation
      DOUBLE PRECISION :: IC_VEL(3), VEL_SIG
! Placeholder for random position and velocity
      DOUBLE PRECISION, ALLOCATABLE :: randVEL(:,:)
      DOUBLE PRECISION :: RAND(3)
! Local stat wt
      DOUBLE PRECISION :: STAT_WT
! Volume of a parcel and total solids volume plus related variables
      DOUBLE PRECISION :: sVOL, SVOL_TOT, EFF_RAD, remainder
! Counter for seeded parcels
      INTEGER :: SEEDED
! Generic loop indices and counters
      INTEGER :: I, J, K, IJK, LC, LC_MAX, N
      INTEGER :: dg_I, dg_J, dg_K, dg_IJK, skipped
! Left and right hand endpoints of bins
      DOUBLE PRECISION :: LEFT_END, RIGHT_END
! Incremental distance between the Z-bins
      DOUBLE PRECISION :: Z_INCREMENT
! Array that holds CDF of polydisperse function
      DOUBLE PRECISION :: ZCDF(dim_pic_mi_bin)
! Particle count in each bin (based on average diameter of the bin)
      DOUBLE PRECISION :: PARTICLES_IN_BIN_byvol
      DOUBLE PRECISION :: PARTICLES_IN_BIN_bymass
      DOUBLE PRECISION :: PARTICLES_IN_BIN
! Median diameter of a perticle in bin
      DOUBLE PRECISION :: MEDIAN_DIA(dim_pic_mi_bin)
      DOUBLE PRECISION :: PSD_MEAN_DP, PSD_STDEV
      DOUBLE PRECISION :: PSD_MIN_DP, PSD_MAX_DP
! Conversion of PSD for Log-normal distribution
      DOUBLE PRECISION :: PSD_MEAN_DP_LN, PSD_STDEV_LN
! Volume variables used to make particle counts
      DOUBLE PRECISION :: VOL_OF_PHASE
      DOUBLE PRECISION :: VOL_PARTICLE
! Number variables used to make particle counts
      DOUBLE PRECISION :: PART_ESTIMATE
      DOUBLE PRECISION :: NUM_ERROR
      DOUBLE PRECISION :: VOL_MEAN
! Volume of the IC
      DOUBLE PRECISION :: IC_VOL
! Number of bins
      INTEGER :: NUMBIN
! Keyframe management for custom distributions
      INTEGER :: KF, NKF
!      DOUBLE PRECISION :: TEST_DIA
      INTEGER :: localM

      LOGICAL :: ldebug=.FALSE.

!---------------------------------------------------------------------//

      maxPARTS = 25
      allocate(randVEL(maxPARTS,3))

!     set initial velocity for all parcels in phase M
      IC_VEL(1) = IC_U_s(ICV,M)
      IC_VEL(2) = IC_V_s(ICV,M)
      IC_VEL(3) = merge(IC_W_s(ICV,M),0.0d0,DO_K)

      VEL_SIG = sqrt(IC_Theta_M(ICV,M))

      STAT_WT = IC_PIC_CONST_STATWT(ICV,M)
      localRO = IC_RO_S(ICV,M)
      SEEDED = 0
      sVOL_TOT = 0.0d0

!     incorporate known information from polydisperse functions
      PSD_MEAN_DP = IC_PSD_MEAN_DP(ICV,M)
      PSD_STDEV   = IC_PSD_STDEV(ICV,M)
      PSD_MIN_DP  = IC_PSD_MIN_DP(ICV,M) !Lower Clip-off
      PSD_MAX_DP  = IC_PSD_MAX_DP(ICV,M) !Upper Clip-off

      NUMBIN = PIC_NUMBIN ! Histogram setting for NORMAL and LOG-NORMAL
                          ! PIC_NUMBIN is set to 100 in mfix_pic_mod.f90
                          ! We could let the user set this.
!     There are NUMBIN bins, and (NUMBIN+1) control points.
!     Bin I is bounded by control points I and (I+1).
!     Control point 1 correspond to PSD_MIN_DP (minimum diameter).
!     Control point (NUMBIN+1) correspond to PSD_MAX_DP (maximum diameter).
!     MEDIAN_DIAM(I) the diameter at the center of bin I.

!     Manage normal and log_normal distribution types
!     added custom
      IF (IC_PSD_TYPE(ICV,M).EQ.'NORMAL'.OR.  &
          IC_PSD_TYPE(ICV,M).EQ.'LOG_NORMAL'.OR. &
          IC_PSD_TYPE(ICV,M).EQ.'CUSTOM') THEN

!        Create increment for the bins
         Z_INCREMENT = (PSD_MAX_DP-PSD_MIN_DP)/REAL(NUMBIN)

         if(IC_PSD_TYPE(ICV,M).EQ.'NORMAL') then
            ! Create the CDF at each control point and diameter at center of each bin.
            ! Note that ZCDF(1) is not zero and ZCDF(NUMBIN+1) is not one due to clipping.
            ZCDF(1)  = HALF * (1.0d0 + ERF((PSD_MIN_DP - PSD_MEAN_DP)/ &
                             (SQRT(2.0d0)*PSD_STDEV)))
            DO I = 1, NUMBIN, 1
               MEDIAN_DIA(I) = PSD_MIN_DP + (FLOAT(I)-HALF)*Z_INCREMENT
               RIGHT_END     = PSD_MIN_DP +  FLOAT(I)      *Z_INCREMENT
               ZCDF(I+1)     = HALF * (1.0d0 + ERF((RIGHT_END-PSD_MEAN_DP)/ &
                                      (SQRT(2.0d0)*PSD_STDEV)))
            END DO

         endif !normal


         if(IC_PSD_TYPE(ICV,M).EQ.'LOG_NORMAL')then
         !convert mean and stdev to log normal form
            PSD_MEAN_DP_LN = LOG(PSD_MEAN_DP**2/ &
                                 sqrt(PSD_MEAN_DP**2 + PSD_STDEV**2))
            PSD_STDEV_LN = sqrt(LOG(ONE+PSD_STDEV**2/PSD_MEAN_DP**2))
         ! Create the CDF at each control point and diameter at center of each bin.
         ! Note that ZCDF(1) is not zero and ZCDF(NUMBIN+1) is not one due to
         ! clipping
            ZCDF(1)     = HALF + HALF*ERF((LOG(PSD_MIN_DP) - PSD_MEAN_DP_LN)/ &
                                         (SQRT(2.0d0)*PSD_STDEV_LN))
            DO I = 1, NUMBIN, 1
               MEDIAN_DIA(I) = PSD_MIN_DP + (FLOAT(I)-HALF)*Z_INCREMENT
               RIGHT_END     = PSD_MIN_DP +  FLOAT(I)      *Z_INCREMENT
               ZCDF(I+1)     = HALF + HALF*ERF((LOG(RIGHT_END) - PSD_MEAN_DP_LN)/ &
                                               (SQRT(2.0d0)*PSD_STDEV_LN))
            END DO
         endif  !lognormal

         IF (IC_PSD_TYPE(ICV,M).EQ.'CUSTOM') THEN
           !pick up the keyframe count for the distribution
           KF = kf_count(1) !1 means IC
           !isolate which keyframe carries the IC/M combo
           if (KF > 0) then
           !pick which KF is associated with the current IC
              DO NKF = 1, kf_count(1), 1 !1 means IC
                 if(rlist_CV(NKF,1).EQ.ICV)then
                   KF = NKF
                 end if
              END DO
           !check which phase is associated with that KF
              localM = rlist_M(KF,1)
           endif !kf
           !check if custom distribution for the ICV and phase match
           if (localM.eq.M) then
           ! The custom PSD is written at each control point (one control point per line).
           ! The number of bins will be equal to the number of control points minus one.
           ! Note that here ZCDF(1) = 0 and ZCDF(NUMBIN+1) = 1.0
              NUMBIN = dist_data(KF,1)%nrows - 1
              ZCDF(1) = dist_data(KF,1)%y(1) ! This should always be zero
              DO I = 1, NUMBIN, 1
                 !data is stored under dist_data
                 !pick up the median diameters and associate CDF
                 MEDIAN_DIA(I) = HALF * (dist_data(KF,1)%x(I) + dist_data(KF,1)%x(I+1))
                 ZCDF(I+1) = dist_data(KF,1)%y(I+1)
              END DO
           else
              !there is no associated custom distribution and error out
              write(*,*) "Something went wrong with CUSTOM ICV=",ICV,&
                     "and PHASE=",localM
              ! FIXME make this a LOG_ERROR message
           endif !phase check
           if (ldebug) then
              write(*,*) "NUMBIN = ", NUMBIN
              write(*,*) "Local M = ", localM, "M = ", M
              write(*,*) "KF = ", KF
              do I=1,NUMBIN,1
                 write(*,*) "ZCDF(",I,")=", ZCDF(I), &
                            "    MEDIAN_DIA(",I,")=", MEDIAN_DIA(I)
              end do
           end if
         END IF !custom


!     Calculate volume of solids phase in a particular region
         CALL GET_IC_VOLUME(ICV,IC_VOL)
         VOL_OF_PHASE = IC_EP_S(ICV,M)*IC_VOL
!     Estimate number of particles using mean particle volume of system
!     Estimate must be scaled by CDF for non-symmetric systems
         VOL_MEAN = mean_volume(1,ICV,M)
         PART_ESTIMATE = (VOL_OF_PHASE/VOL_MEAN)/(ZCDF(NUMBIN + 1)-ZCDF(1))
         if (ldebug) write(*,*) "PART_ESTIMATE = ", PART_ESTIMATE

         NUM_ERROR = 0.0d0
         DO N = 1, NUMBIN, 1
!           volume of single particle in bin
            VOL_PARTICLE = (PI/6.0d0)*(MEDIAN_DIA(N)**3.0d0) !vol approach
!           number of particles in the bin including fractional part
!           (number based)
            PARTICLES_IN_BIN = (ZCDF(N+1)-ZCDF(N))*PART_ESTIMATE + NUM_ERROR

            if(ldebug) then
               write(*,*) "MEDIAN_DIA(",N,")=", MEDIAN_DIA(N), &
                     "   CDF(",N,")=", ZCDF(N), "   CDF(",N+1,")=", ZCDF(N+1),&
                       "   NUMBER OF PARTICLES = ", PARTICLES_IN_BIN
            end if
!           check if any particles in the bin
            IF(PARTICLES_IN_BIN<1.0d0) cycle

!           now we allocate these particles (as parcels) over the region before we
!           calculate the next set

!           recycling original seeding logic but employing poly properties
            remainder = 0.0d0

!           Loop over the IC region and populate with PARTICLES_IN_BIN(N)
            DO K = IC_K_B(ICV), IC_K_T(ICV)
            DO J = IC_J_S(ICV), IC_J_N(ICV)
            DO I = IC_I_W(ICV), IC_I_E(ICV)

!           Cycle out of non-productive cells
               IF(.not.IS_ON_myPE_wobnd(I,J,K)) cycle
               IF(DEAD_CELL_AT(I,J,K)) cycle

               IJK = FUNIJK(I,J,K)
               IF(.not.FLUID_AT(IJK)) cycle

!           recast as number based distribution
                 rPARTS = (VOL(IJK)/IC_VOL)* &
                          PARTICLES_IN_BIN + remainder

!           Increase particle buffer if needed

!           Number of parcels in a particular cell.
               LC_MAX = floor(rPARTS/STAT_WT)
!           Keep track of leftovers in remainder and add back in on next sweep
               remainder = max(rPARTS-dble(LC_MAX)*STAT_WT,0.0d0)

               IF(LC_MAX == 0) cycle

               IF(LC_MAX > maxPARTS) THEN
                   maxPARTS = 2*LC_MAX
                   if(allocated(randVEL)) deallocate(randVEL)
                   allocate(randVEL(maxPARTS,3))
               END IF

               DO LC=1, merge(2,3,NO_K)
                  IF(VEL_SIG > ZERO) THEN
                     CALL NOR_RNO(randVEL(1:LC_MAX,LC),IC_VEL(LC), VEL_SIG)
                  ELSE
                     randvel(1:LC_MAX,3) = IC_VEL(LC)
                  END IF
               END DO

               IF(NO_K) randVEL(1:LC_MAX,3) = 0.0d0

               IC_START(1) = XE(I-1)
               IC_START(2) = YN(J-1)
               IC_START(3) = ZERO;  IF(DO_K) IC_START(3) = ZT(K-1)

               DOML(1) = DX(I)
               DOML(2) = DY(J)
               DOML(3) = ZERO;  IF(DO_K) DOML(3) = DZ(K)

               EFF_RAD = MEDIAN_DIA(N)*HALF*(STAT_WT)**(1.0/3.0)

               PIC_LP:  DO LC=1, LC_MAX

                  CALL RANDOM_NUMBER(RAND)
                  POS(:) = IC_START + DOML*RAND

                  dg_i = min(dg_iend2,max(dg_istart2, iofpos(pos(1))))
                  dg_j = min(dg_jend2,max(dg_jstart2, jofpos(pos(2))))
                  dg_k = min(dg_kend2,max(dg_kstart2, kofpos(pos(3))))
                  dg_ijk = dg_funijk(dg_i,dg_j,dg_k)

                  skipped = 0
                  do while(picSTLoverlap(dg_ijk, pos, eff_rad))

                     skipped = skipped + 1
                     CALL RANDOM_NUMBER(RAND)
                     POS(:) = IC_START + DOML*RAND

                     dg_i = min(dg_iend2,max(dg_istart2, iofpos(pos(1))))
                     dg_j = min(dg_jend2,max(dg_jstart2, jofpos(pos(2))))
                     dg_k = min(dg_kend2,max(dg_kstart2, kofpos(pos(3))))
                     dg_ijk = dg_funijk(dg_i,dg_j,dg_k)

                     if(skipped >= 1000) cycle PIC_LP
                  end do

                  PIP = PIP + 1
                  CALL PARTICLE_GROW(PIP)
                  MAX_PIP = max(PIP,MAX_PIP)

                  DES_POS_NEW(PIP,:) = POS(:)
                  DES_VEL_NEW(PIP,:) = 1.0e-4*rand

                  DES_RADIUS(PIP) = MEDIAN_DIA(N)*HALF
                  RO_SOL(PIP) = localRO

                  PIJK(PIP,1) = I
                  PIJK(PIP,2) = J
                  PIJK(PIP,3) = K
                  PIJK(PIP,4) = IJK
                  PIJK(PIP,5) = M

                  DES_STAT_WT(PIP) = STAT_WT
                  sVOL_TOT = sVOL_TOT + VOL_PARTICLE*STAT_WT

                  RESIDENCE_TIME(PIP) = ZERO

                  CALL SET_NORMAL(PIP)

                  SEEDED = SEEDED + 1

               END DO PIC_LP

            END DO !I
            END DO !J
            END DO !K

            ! apply strict inequality-  error lost in last bin
            if (N < NUMBIN)  NUM_ERROR = remainder * (MEDIAN_DIA(N)**3/MEDIAN_DIA(N+1)**3)

         END DO !N

      END IF !NORMAL

      IF(allocated(randvel)) deallocate(randVEL)
      sDATA(1) = dble(SEEDED)
      sDATA(2) = sVOL_TOT

      RETURN




      END SUBROUTINE POLYPIC_GPC_MPPIC_CONST_STATWT

END MODULE GENERATE_PARTICLES
