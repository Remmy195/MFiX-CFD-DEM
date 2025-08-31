#include "error.inc"

MODULE MAKE_ARRAYS_DES_MOD

   use calc_interp_weights_mod, only: calc_interp_weights
   use cfassign_mod, only: cfassign
   use error_manager
   use init_settling_dem_mod, only: init_settling_dem
   use neighbour_mod, only: neighbour
   use particles_in_cell_mod, only: particles_in_cell, init_particles_in_cell
   use read_part_input_mod, only: read_part_input
   use read_res0_des_mod, only: read_res0_des
   use set_filter_des_mod, only: set_filter_des
   use set_ic_dem_mod, only: set_ic_dem
   use set_phase_index_mod, only: set_phase_index
   use write_des_data_mod, only: write_des_data
   use sq_properties_mod

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: MAKE_ARRAYS_DES                                        !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: DES - allocating DES arrays                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MAKE_ARRAYS_DES

      use calc_collision_wall
      use comp_mean_fields_mod, only: comp_mean_fields
      use compar
      use constant, only: pi
      use cutcell
      use des_functions_mod, only: des_sort_particle_arrays, des_getvalid_fluid_cells
      use des_rxns
      use des_thermo
      use desgrid
      use discretelement
      use functions
      use funits
      use generate_particles, only: generate_particle_config
      use geometry
      USE mfix_pic, only: MPPIC
      use mpi_funs_des, only: des_par_exchange
      use mpi_funs_pic, only: pic_par_exchange
      use mpi_utility
      use param, only: dimension_3, dimension_3_alloc, DIMENSION_IC
      use param1
      use run
      use stl
      use stl_functions_des
      use stl_preproc_des, only: add_facet
      use des_rxns, only: SAVE_PART_RRATES, Part_RRates_out
      use read_part_input_mod, only: READ_SPHERE_INPUT_GSP
      use check_solids_dem_mod, only: check_solids_dem
      use ic,only: IC_OVERRIDE_PART_INPUT
      !use GSP_MATH_MOD, only: gsp_result_type

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: I, J, K, L, IJK
      INTEGER :: I1, I2, J1, J2, K1, K2, II, JJ, KK, IJK2
      INTEGER :: lcurpar, lpip_all(0:numpes-1), lglobal_id, gid
      DOUBLE PRECISION :: axi(3), m,n, IXX, IYY, IZZ, SVOLUME,NS
      DOUBLE PRECISION :: sphere_radius,xloc,yloc,zloc,Itemp,actual_b,dist2,sphere_mass_tmp
      INTEGER :: rktest(4)
      integer :: tmp_size, n_rows,PP,MM, gsp_implicit_output_size

! Check interpolation input.
      CALL SET_FILTER_DES

! cfassign and des_init_bc called before reading the particle info
      CALL CFASSIGN

      VOL_SURR(:) = ZERO

      ! initialize VOL_SURR array
      DO K = KSTART2, KEND1
         DO J = JSTART2, JEND1
            DO I = ISTART2, IEND1
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = funijk(I,J,K)
               I1 = I
               I2 = I+1
               J1 = J
               J2 = J+1
               K1 = K
               K2 = merge(K, K+1, NO_K)

! looping over stencil points (node values)
               DO KK = K1, K2
                  DO JJ = J1, J2
                     DO II = I1, I2
                        IF (DEAD_CELL_AT(II,JJ,KK)) CYCLE  ! skip dead cells
                        IJK2 = funijk_map_c(II, JJ, KK)
                        IF(FLUID_AT(IJK2)) VOL_SURR(IJK) = &
                        VOL_SURR(IJK)+VOL(IJK2)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

! Set the initial particle data.
      IF(RUN_TYPE == 'NEW') THEN
         IF (GluedSphereDEM) THEN
            ! This is now reading the gsp_config.csv instead of dat
            ! A migration code is attached to convert dat to csv if csv is not found
            CALL READ_SPHERE_INPUT_GSP

            IF (GENER_PART_CONFIG) THEN
               if(mype==pe_io) write(*,*) "Enable GSP model and perform auto-seeding"
               CALL GENERATE_PARTICLE_CONFIG
            ELSE
               if(mype==pe_io) write(*,*) "Enable GSP model and read the particle input file"
               CALL READ_PART_INPUT
            ENDIF
         ELSE
            IF (GENER_PART_CONFIG) THEN
               CALL GENERATE_PARTICLE_CONFIG
            ELSE
               CALL READ_PART_INPUT
            ENDIF
         ENDIF

! Set the global ID for the particles and set the ghost cnt
         ighost_cnt = 0
         lpip_all = 0
         lpip_all(mype) = pip
         CALL global_all_sum(lpip_all)
         lglobal_id = sum(lpip_all(0:mype-1))
         imax_global_id = 0
! DEM set the global id based on the order of PE in MPI
         DO lcurpar  = 1,pip
            lglobal_id = lglobal_id + 1
            iglobal_id(lcurpar) = lglobal_id
            imax_global_id = iglobal_id(pip)
         ENDDO
! Sort the components global id
         IF(gsp_explicit .AND. numpes > 1) THEN
            DO lcurpar = 1, MAX_PIP
               IF(PARTICLE_STATE(lcurpar) .NE. NORMAL_PARTICLE) CYCLE
               iglobal_id(lcurpar) = cglobal_id(lcurpar)
               imax_global_id = MAX(iglobal_id(lcurpar),imax_global_id)
            ENDDO
         ENDIF

         CALL GLOBAL_ALL_MAX(imax_global_id)

! deallocate cglobal_id and do not use anymore
         IF(ALLOCATED(cglobal_id)) DEALLOCATE(cglobal_id)
         ! only when mpirun > 1 & gsp model & marked_gsp is not yet allocated
         ! marked_gsp and marked_state should only allocated at the same time
         IF(gsp_explicit .AND. numpes > 1 .AND. (.NOT. ALLOCATED(marked_gsp))) THEN
            tmp_size = NGluedParticles/numpes
            tmp_size = MAX(tmp_size,4)
            ALLOCATE(marked_gsp(tmp_size)); marked_gsp(:) = 0
            ALLOCATE(marked_state(tmp_size)); marked_state(:) = -1
         ENDIF
! Initialize old values
         omega_new(:,:) = zero

! Particle orientation
         IF(PARTICLE_ORIENTATION) THEN
            IF(.NOT.SuperDEM) THEN
               ORIENTATION(:,1) = INIT_ORIENTATION(1)
               ORIENTATION(:,2) = INIT_ORIENTATION(2)
               ORIENTATION(:,3) = INIT_ORIENTATION(3)
            ENDIF
         ENDIF

         IF (DO_OLD) THEN
            omega_old(:,:)   = zero
            des_pos_old(:,:) = des_pos_new(:,:)
            des_vel_old(:,:) = des_vel_new(:,:)
         ENDIF

         IF(SAVE_PART_RRATES) Part_RRates_out(:,:) = ZERO

! Read the restart file.
      ELSEIF(RUN_TYPE == 'RESTART_1' .OR. RUN_TYPE == 'RESTART_2') THEN
         IF (GluedSphereDEM) CALL READ_SPHERE_INPUT_GSP
         CALL READ_RES0_DES
         imax_global_id = maxval(iglobal_id(1:pip))
         call global_all_max(imax_global_id)

         if(GluedSphereDEM .and. numpes > 1) then
            allocate(marked_gsp(NGluedParticles/numpes))
            allocate(marked_state(NGluedParticles/numpes))
            marked_gsp(:) = 0
            marked_state(:) = -1
         endif

! Initialize the old values.
         IF (DO_OLD) THEN
            omega_old(:,:)   = omega_new(:,:)
            des_pos_old(:,:) = des_pos_new(:,:)
            des_vel_old(:,:) = des_vel_new(:,:)
         ENDIF
      ELSE
         WRITE(ERR_MSG, 1100)
         CALL LOG_ERROR()
 1100 FORMAT('Error 1100: Unsupported RUN_TYPE for DES.')
      ENDIF

      IF(RUN_TYPE == 'RESTART_2') VTP_FINDEX=0

! setting additional particle properties now that the particles
! have been identified
      NSphereGSP = 0 ! in each rank calculate the number of normal component spheres
                     ! only for the gsp implicit model
! Setting the minimum particle mass of each phase
      IF(ANY_SOLIDS_SPECIES_EQ) DES_MIN_PMASS = UNDEFINED
      DO L = 1, MAX_PIP
! Skip 'empty' locations when populating the particle property arrays.
         IF(IS_NONEXISTENT(L)) CYCLE
         IF(IS_ANY_GHOST(L)) CYCLE

! particle volume, mass and interia for superquadric particle
         IF (SuperDEM)   then
            axi(1)=super_r(L,1)
            axi(2)=super_r(L,2)
            axi(3)=super_r(L,3)
            m=super_mn(L,1)
            n=super_mn(L,2)
            CALL SQ_VOLUME (axi,m,n,svolume)
            PVOL(L)=svolume
            PMASS(L) = PVOL(L)*RO_SOL(L)
            call SQ_INERTIA(axi,m,n,IXX,IYY,IZZ)
            OMOI3(L,1)=1.0/(IXX*RO_SOL(L))
            OMOI3(L,2)=1.0/(IYY*RO_SOL(L))
            OMOI3(L,3)=1.0/(IZZ*RO_SOL(L))
         ELSEIF(GSP_IMPLICIT) THEN
            MM = PIJK(L,5)
            IF(sphereNumber(MM) .eq. 1) THEN
               ! if only has 1 components, bounding diameter == actual diameter
               PVOL(L) = (4.0D0/3.0D0)*PI*DES_RADIUS(L)**3.0d0
               PMASS(L) = PVOL(L)*RO_SOL(L)
               OMOI3(L, 1) = 2.5D0/(PMASS(L)*DES_RADIUS(L)**2)
               OMOI3(L, 2) = 2.5D0/(PMASS(L)*DES_RADIUS(L)**2)
               OMOI3(L, 3) = 2.5D0/(PMASS(L)*DES_RADIUS(L)**2)
               glueSurface(L) = 4*PI*DES_RADIUS(L)**2.0
               NSphereGSP = NSphereGSP + 1
            ELSE
               ! @ renjieke need to recompute pvol and pmass?
               PVOL(L) = 0.0D0
               OMOI3(L,:) = 0.0D0
               glueSurface(L) = 0.0D0
               n_rows = size(cspart_data,1)
               ! This des_radius already consider the difference between poly-disperse and mono-disperse cases
               actual_b = 2.0 * des_radius(L)
               NSphereGSP = NSphereGSP + sphereNumber(MM)
               DO PP = 1, n_rows
                  if(cspart_data(PP,4,MM) == 0.0) cycle
                  sphere_radius = cspart_data(PP,4,MM) * actual_b / 2.0
                  xloc = cspart_data(PP,1,MM) * actual_b
                  yloc = cspart_data(PP,2,MM) * actual_b
                  zloc = cspart_data(PP,3,MM) * actual_b
                  glueSurface(L) = glueSurface(L) + cspart_data(PP,5,MM) * actual_b * actual_b
                  PVOL(L) = PVOL(L) + 4.0/3.0*PI*sphere_radius**(3.0)

                  sphere_mass_tmp = RO_SOL(L) * 4.0D0 * PI / 3.0D0 * sphere_radius ** 3.0D0
                  Itemp = (2.0/5.0)*sphere_mass_tmp*sphere_radius**2.0

                  dist2 = yloc ** 2.0 + zloc ** 2.0
                  OMOI3(L,1) = OMOI3(L,1) + sphere_mass_tmp * dist2 + Itemp

                  dist2 = xloc ** 2.0 + zloc ** 2.0
                  OMOI3(L,2) = OMOI3(L,2) + sphere_mass_tmp * dist2 + Itemp

                  dist2 = xloc ** 2.0 + yloc ** 2.0
                  OMOI3(L,3) = OMOI3(L,3) + sphere_mass_tmp * dist2 + Itemp
               ENDDO
               OMOI3(L, 1) = 1.0d0/OMOI3(L, 1)
               OMOI3(L, 2) = 1.0d0/OMOI3(L, 2)
               OMOI3(L, 3) = 1.0d0/OMOI3(L, 3)
               PMASS(L) = PVOL(L)*RO_SOL(L)
               glueDiameter(L) = (6.0D0 * PVOL(L)/PI) ** (1.0D0/3.0D0)
               glueBounding(L) = des_radius(L) * 2.0D0
               IF(PARTICLE_ORIENTATION) &
                 CALL QROTATE(glueQuat(L,:),ORIENTATION(L,:),INIT_ORIENTATION,2)
            ENDIF
         ELSE
           PVOL(L) = (4.0D0/3.0D0)*PI*DES_RADIUS(L)**3
           PMASS(L) = PVOL(L)*RO_SOL(L)
! for spherical particle, IXX=IYY=IZZ
           OMOI(L) = 2.5D0/(PMASS(L)*DES_RADIUS(L)**2) !ONE OVER MOI
         ENDIF

         MM = PIJK(L,5)
         IF(ANY_SOLIDS_SPECIES_EQ) THEN
            IF(PMASS(L) .LT. DES_MIN_PMASS(MM)) DES_MIN_PMASS(MM)=PMASS(L)
         ENDIF
      ENDDO
      IF(ANY_SOLIDS_SPECIES_EQ) CALL global_all_min(des_min_pmass)

      IF(GSP_IMPLICIT) THEN
      ! @renjieke only for saving purpose, no calculation, no info exchange
      ! lets say if max_pip increased, NSphereGSP need to be changed as well, so the size of below needs to grow
         gsp_implicit_output_size = max(NSphereGSP, 4)
         allocate( sphere_radius_out(gsp_implicit_output_size) )
         allocate( sphere_rho_out(gsp_implicit_output_size) )
         allocate( sphere_mass_out(gsp_implicit_output_size) )
         allocate( sphere_vol_out(gsp_implicit_output_size) )
         allocate( sphere_rank_out(gsp_implicit_output_size) )
         allocate( sphere_pid_out(gsp_implicit_output_size) )
         allocate( sphere_id_out(gsp_implicit_output_size) )
         allocate( sphere_res_time_out(gsp_implicit_output_size) )
         allocate( sphere_coordination_out(gsp_implicit_output_size) )
         allocate( sphere_MIN_OVERLAP_out(gsp_implicit_output_size) )
         allocate( sphere_MAX_OVERLAP_out(gsp_implicit_output_size) )
         allocate( sphere_MEAN_OVERLAP_out(gsp_implicit_output_size) )

         allocate( sphere_rot_out(gsp_implicit_output_size, 3) )
         allocate( sphere_vel_out(gsp_implicit_output_size, 3) )
         allocate( sphere_orient_out(gsp_implicit_output_size, 3) )

         allocate( sphere_col_force_out(gsp_implicit_output_size, 3) )
         allocate( sphere_drag_out(gsp_implicit_output_size) )

         ! Dump gsp total mass and bounding diameter to each component sphere in implicit model
         allocate( sphere_total_mass_out(gsp_implicit_output_size) )
         allocate( sphere_bounding_out(gsp_implicit_output_size) )
         allocate( sphere_eqv_dia_out(gsp_implicit_output_size) )

         IF(ENERGY_EQ) THEN
            allocate( sphere_temp_out(gsp_implicit_output_size) )
            allocate( sphere_cps_out(gsp_implicit_output_size) )
         ENDIF
         if(DES_USR_VAR_SIZE > 0) &
            allocate( sphere_usr_var_out(DES_USR_VAR_SIZE, gsp_implicit_output_size) )
         IF(ANY_SOLIDS_SPECIES_EQ) &
            allocate( sphere_xs_out(gsp_implicit_output_size, DIMENSION_N_S) )
         if (SAVE_PART_RRATES) &
            allocate( sphere_prates_out(gsp_implicit_output_size, NO_OF_DES_RXNS) )
      ENDIF

      CALL SET_PHASE_INDEX
      CALL INIT_PARTICLES_IN_CELL

! do_nsearch should be set before calling particle in cell
      DO_NSEARCH =.TRUE.
! Bin the particles to the DES grid.
      CALL DESGRID_PIC(PLOCATE=.TRUE.)
      IF(MPPIC) THEN
         CALL PIC_PAR_EXCHANGE
      ELSE
         CALL DES_PAR_EXCHANGE
      ENDIF

      CALL PARTICLES_IN_CELL
!Sitaraman======================================
if (optflag1.eq.1) then
      CALL DES_SORT_PARTICLE_ARRAYS
      CALL DES_GETVALID_FLUID_CELLS
!Initialize vol_surr_inv
      vol_surr_inv(:)=ZERO
      do i=1,DIMENSION_3_ALLOC
         if(vol_surr(i) .gt. ZERO) then
            vol_surr_inv(i)=ONE/vol_surr(i)
         endif
      enddo

!Initialize fluid_at_mask1
      fluid_at_mask1=ZERO
      do ijk=1,DIMENSION_3_ALLOC
         if(fluid_at(ijk)) fluid_at_mask1(ijk)=ONE
      enddo

!Initialize fluid_at_mask2
      fluid_at_mask2=ZERO
      DO K = KSTART2, KEND2
         DO J = JSTART2, JEND2
            DO I = ISTART2, IEND2

               ijk=funijk_map_c(i,j,k)
               if(fluid_at(ijk) .and. IS_ON_myPE_wobnd(I,J,K)) then
                  fluid_at_mask2(i-ISTART2+1,j-JSTART2+1,k-KSTART2+1)=ONE
               endif

            enddo
         enddo
      enddo

!Initialize fluid_at_mask3
      fluid_at_mask3=ZERO
      do ijk=1,dimension_3_alloc
         if(vol_surr(ijk) .gt. ZERO) fluid_at_mask3(ijk)=ONE
      enddo

!Initialize do_k mask (vaidhynathan)
      IF(DO_K) THEN
        do_k_mask = ONE
      ELSE
        do_k_mask = ZERO
      ENDIF
endif
!==================================================

      IF(DEM_SOLIDS) THEN
         CALL NEIGHBOUR
         CALL INIT_SETTLING_DEM
      ENDIF

      IF(RUN_TYPE == 'NEW') CALL SET_IC_DEM
! After SET_IC_DEM, des_vel_new might changed, so perform the info exchange again
! to make sure ghost info is properly updated
! this only happens if you have at least one IC choose to overwrite its initial values
      IF(ANY(IC_OVERRIDE_PART_INPUT)) THEN
         IF(numPEs > 1) THEN
            CALL DESGRID_PIC(.TRUE.)
            CALL DES_PAR_EXCHANGE
         ENDIF
         ! @renjieke 8-9-2024, gsp_explicit requires updating of glueVel in DMP
         IF(gsp_explicit) THEN
            glueVel(:,:) = ZERO
            DO L = 1,MAX_PIP
               IF(PARTICLE_STATE(L) .NE. NORMAL_PARTICLE) CYCLE
               glueVel(gid_list(L),:) = des_vel_new(L,:)
            ENDDO
            IF(numPEs > 1) CALL GLOBAL_ALL_SUM(glueVel)
         ENDIF
      ENDIF

! Calculate interpolation weights
      CALL CALC_INTERP_WEIGHTS
! Calculate mean fields using either interpolation or cell averaging.
      CALL COMP_MEAN_FIELDS

      IF(RUN_TYPE /= 'RESTART_1' .AND. PRINT_DES_DATA) THEN
         S_TIME = TIME
         CALL WRITE_DES_DATA
      ENDIF

! Set time when DTSOLID will be updated next
      DTSOLID_UPDATE_TIME = TIME + DTSOLID_UPDATE_DT

! Check GSP solids here to update DTSOLID for GSP after seeding
      IF(GluedSphereDEM) CALL CHECK_SOLIDS_DEM

      RETURN

   END SUBROUTINE MAKE_ARRAYS_DES

END MODULE MAKE_ARRAYS_DES_MOD
