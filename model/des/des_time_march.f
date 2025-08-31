#include "error.inc"

MODULE DES_TIME_MARCH

      use bc

      use calc_collision_wall, only: calc_dem_thermo_with_wall_stl
      use calc_drag_des_mod, only: calc_drag_des
      use calc_force_dem_mod, only: calc_force_dem, calc_coordination
      use SQ_CALC_FORCE_MOD, only: CALC_FORCE_SuperDEM
      use calc_interp_weights_mod, only: calc_interp_weights
      use calc_pg_grad_mod, only: calc_pg_grad
      use calc_rrates_des_mod, only: zero_rrate_des
      use calc_thermo_des_mod, only: calc_thermo_des
      use cfnewvalues_mod, only: cfnewvalues, SuperDEM_cfnewvalues
      use cfupdateold_mod, only: cfupdateold
      use comp_mean_fields_mod, only: comp_mean_fields
      use compar, only: ADJUST_PARTITION
      use conv_gs_des1_mod, only: zero_energy_source, conv_gs_des1
      use des_bc, only: DEM_BCMI, DEM_BCMO
      use des_functions_mod, only: des_sort_particles_spatially, des_sort_particle_arrays, des_getvalid_fluid_cells
      use des_radiation_mod, only: calc_avgts
      use des_reaction_model_mod, only: des_reaction_model
      use des_thermo, only: CALC_RADT_DES
      use des_thermo_newvalues_mod, only: des_thermo_newvalues, GSP_HEAT_COND_EXPLICIT
      use desgrid, only: desgrid_pic
      use CHECK_SOLIDS_DEM_MOD, only: update_dtsolid
      use discretelement
      use drag_gs_des1_mod, only: drag_gs_des1
      use error_manager
      USE get_stl_data_mod, only: move_is_stl
      use keyframe_mod
      use mass_inflow_dem_mod, only: mass_inflow_dem
      use mass_outflow_dem_mod, only: mass_outflow_dem, update_exit_des_vel
      use mpi_funs_des, only: DESMPI_SEND_RECV_FIELD_VARS
      use mpi_funs_des, only: DES_PAR_EXCHANGE
      use mpi_utility
      use neighbour_mod, only: neighbour
      use output, only: DLB,DLB_TIME
      use output_manager, only: write_outputs
      use output_manager, only: UPDATE_TIME_STEP_INFO_INT
      use output_manager, only: PRINT_TIME_STEP_INFO
      use particles_in_cell_mod, only: particles_in_cell
      use run , only: optflag1
      use run, only: ANY_SOLIDS_SPECIES_EQ
      use run, only: CALL_USR
      use run, only: ENERGY_EQ
      use run, only: NSTEP
      use run, only: TIME, TSTOP, DT, DT_PREV
      use run, only: FREEZE_EXIT_DES_VEL
      use run, only: USE_DT_PREV
      use sendrecv
      use time_mod, only: get_elapsed_time
      use vtk, only: vtk_time, vtk_dt, dimension_vtk, write_vtk_files
      use vtp, only: write_vtp_file, vtk_data
      use CFNEWVALUES_MOD, only: GSP_CFNEWVALUES_vhalfxfull, GSP_CFNEWVALUES_vfull
      use CFNEWVALUES_MOD, only: GSP_CFNEWVALUES_vhalfxfull_implicit, GSP_CFNEWVALUES_vfull_implicit
      use CFNEWVALUES_MOD, only: UPDATE_GSP_PARTICLE_STATE
      use stiff_chem, only : STIFF_CHEMISTRY
!---------------------------------------------------------------------//
! Total number of particles
      INTEGER, SAVE :: NP=0

! loop counter index for any initial particle settling in coupled cases
      INTEGER :: FACTOR
! Temporary variables when des_continuum_coupled is T to track
! changes in solid time step
      DOUBLE PRECISION :: DTSOLID_TMP
! Numbers to calculate wall time spent in DEM calculations.
      DOUBLE PRECISION :: TIME_TMP

! Temporary fluid time step.
      DOUBLE PRECISION :: DT_TMP

      LOGICAL :: EXIT_LOOP

!......................................................................!

   CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Subroutine: DES_TIME_INIT                                        !
!     Author: Jay Boyalakuntla                        Date: 21-Jun-04  !
!                                                                      !
!     Purpose: Main DEM driver routine                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_TIME_INIT
         use run, only: solve_ros
         use constant, only: pi
         USE des_thermo, only: des_t_s
         IMPLICIT NONE

         INTEGER :: lcurpar, gid
         DOUBLE PRECISION :: moi, dist2, scale_factor
         DOUBLE PRECISION, allocatable, dimension(:) :: new_total_mass
         INTEGER :: lcurpos, nbIDX, J, locpp, new_nbidx
         DOUBLE PRECISION :: cur_nsa, next_nsa
         LOGICAL :: TMP_NSEARCH
         LOGICAL :: IS_CONSTANT_DENSITY = .FALSE.
         EXIT_LOOP = .FALSE.

! In case of restarts assign S_TIME from MFIX TIME
      S_TIME = TIME

! JFD: Update DTSOLID
      CALL UPDATE_DTSOLID

      DTSOLID_TMP = DTSOLID
      TIME_TMP = GET_ELAPSED_TIME()

!hangZ: DT was modified after fluid iteration. Use the stored DT for DEM update._
      IF(USE_DT_PREV) THEN
         DT_TMP = DT_PREV
      ELSE
         DT_TMP = DT
      ENDIF
! Initialize time stepping variables for coupled gas/solids simulations.
      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DT_TMP.GE.DTSOLID) THEN
            FACTOR = CEILING(real(DT_TMP/DTSOLID))
         ELSE
            FACTOR = 1
            DTSOLID = DT_TMP
         ENDIF

! Initialize time stepping variable for pure granular simulations.
      ELSE
         FACTOR = CEILING(real((TSTOP-TIME)/DTSOLID))
         DT = DTSOLID
         CALL WRITE_OUTPUTS(.FALSE., .FALSE.)
      ENDIF   ! end if/else (des_continuum_coupled)

      NP = PIP - IGHOST_CNT
      CALL GLOBAL_ALL_SUM(NP)
#ifndef QUIET
      IF(DES_CONTINUUM_COUPLED) THEN
          WRITE(ERR_MSG, 1000) trim(iVal(factor)), trim(iVAL(NP))
          CALL LOG_STATUS()
      ELSE
          WRITE(ERR_MSG, 1100) TIME, DTSOLID, trim(iVal(factor))
          CALL LOG_STATUS()
      ENDIF
1000  FORMAT(/'DEM NITs: ',A,3x,'Total PIP: ', A)
1100  FORMAT(/'Simulation time: ',g12.5,3x,'DT: ',g12.5,3x,'DEM NITs: ',A)
#endif
      IF(DES_CONTINUUM_COUPLED) CALL UPDATE_TIME_STEP_INFO_INT(15,factor)
      CALL UPDATE_TIME_STEP_INFO_INT(16,NP)

      IF(CALL_USR) CALL USR0_DES

!@renjieke, for gsp model, only do explicitly coupled for drag force
      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DES_EXPLICITLY_COUPLED) THEN
            CALL DRAG_GS_DES1
            IF(ENERGY_EQ) CALL CONV_GS_DES1
            IF(ANY_SOLIDS_SPECIES_EQ .AND. (.NOT. stiff_chemistry)) CALL DES_REACTION_MODEL
            ! if stiff solver have done its mapping
            IF(stiff_chemistry .and. gsp_explicit) THEN
               glueVolume(:) = 0.0
               glueInertia(:,:) = 0.0

               allocate(new_total_mass(size(glueMass)))
               new_total_mass(:) = 0.0

               DO lcurpar = 1, max_pip
                  gid = gid_list(lcurpar)
                  if(particle_state(lcurpar) /= 1) cycle
                  if(.not. solve_ros(pijk(lcurpar,5))) then
                     is_constant_density = .true.
                     new_total_mass(gid) = new_total_mass(gid) + pmass(lcurpar)
                  else
                     moi = (2.0d0/5.0d0)*pmass(lcurpar)*des_radius(lcurpar)**2

                     dist2 = sc2gpc_vec(lcurpar,2)*sc2gpc_vec(lcurpar,2) + sc2gpc_vec(lcurpar,3)*sc2gpc_vec(lcurpar,3)
                     glueinertia(gid,1) = glueinertia(gid,1) + pmass(lcurpar) * dist2 + moi

                     dist2 = sc2gpc_vec(lcurpar,1)*sc2gpc_vec(lcurpar,1) + sc2gpc_vec(lcurpar,3)*sc2gpc_vec(lcurpar,3)
                     glueinertia(gid,2) = glueinertia(gid,2) + pmass(lcurpar) * dist2 + moi

                     dist2 = sc2gpc_vec(lcurpar,2)*sc2gpc_vec(lcurpar,2) + sc2gpc_vec(lcurpar,1)*sc2gpc_vec(lcurpar,1)
                     glueinertia(gid,3) = glueinertia(gid,3) + pmass(lcurpar) * dist2 + moi
                     glueVolume(gid) = glueVolume(gid) + pvol(lcurpar)
                  endif
               ENDDO

               IF(numpes > 1 .and. is_constant_density) THEN
                  CALL GLOBAL_ALL_SUM(new_total_mass)
               ENDIF

               IF(is_constant_density) THEN
                  DO lcurpar = 1, max_pip
                     IF(PARTICLE_STATE(lcurpar) /= 1) cycle
                     gid = gid_list(lcurpar)
                     if(.not. solve_ros(pijk(lcurpar,5))) then
                        scale_factor = (new_total_mass(gid)/glueMass(gid)) ** (1.0/3.0)

                        des_radius(lcurpar) = des_radius(lcurpar) * scale_factor
                        pvol(lcurpar) = (4.0/3.0)*pi*des_radius(lcurpar) ** 3.0
                        gp_neighsa(lcurpar,:) = gp_neighsa(lcurpar,:) * scale_factor ** 2.0
                        gp_sa(lcurpar) = gp_sa(lcurpar) * scale_factor ** 2.0

                        sc2gpc_vec(lcurpar,:) = sc2gpc_vec(lcurpar,:) * scale_factor
                        moi = (2.0d0/5.0d0)*pmass(lcurpar)*des_radius(lcurpar)**2

                        dist2 = sc2gpc_vec(lcurpar,2)*sc2gpc_vec(lcurpar,2) + sc2gpc_vec(lcurpar,3)*sc2gpc_vec(lcurpar,3)
                        glueinertia(gid,1) = glueinertia(gid,1) + pmass(lcurpar) * dist2 + moi

                        dist2 = sc2gpc_vec(lcurpar,1)*sc2gpc_vec(lcurpar,1) + sc2gpc_vec(lcurpar,3)*sc2gpc_vec(lcurpar,3)
                        glueinertia(gid,2) = glueinertia(gid,2) + pmass(lcurpar) * dist2 + moi

                        dist2 = sc2gpc_vec(lcurpar,2)*sc2gpc_vec(lcurpar,2) + sc2gpc_vec(lcurpar,1)*sc2gpc_vec(lcurpar,1)
                        glueinertia(gid,3) = glueinertia(gid,3) + pmass(lcurpar) * dist2 + moi

                        glueVolume(gid) = glueVolume(gid) + pvol(lcurpar)

                        ! @renjieke 3-15-2025, keep the commented code for now
                        ! this code block is used to transfer neighbor surface area to particle surface area
                        ! in the release 25.1, this method is not adpoted due to lack of the experimental verification.
                        ! if(gp_sa(lcurpar) /= 0.0) then
                        !    lcurpos = 0
                        !    do nbIDX = 1, 6
                        !       J = gp_neighsid(lcurpar, nbIDX) ! current PIP's neighbor pip id
                        !       if(J < 0) cycle
                        !       do locpp = 1, size(iglobal_id)
                        !          if(iglobal_id(locpp) == J) then
                        !          lcurpos = locpp
                        !          exit
                        !          endif
                        !       enddo

                        !       if(lcurpos .le. 0) cycle

                        !       if(nbidx == 1 .or. nbidx == 3 .or. nbidx == 5) then
                        !       new_nbidx = nbidx + 1
                        !       else
                        !       new_nbidx = nbidx - 1
                        !       endif

                        !       ! here shift the difference in gp_neighsa between two neigh spheres to gp_sa
                        !       cur_nsa = gp_neighsa(lcurpar,nbidx)
                        !       next_nsa = gp_neighsa(lcurpos,new_nbidx)
                        !       if (cur_nsa == next_nsa) cycle ! this guarantee the same pair will not be updated twice
                        !       if (cur_nsa .lt. next_nsa) then
                        !          gp_sa(lcurpos) = gp_sa(lcurpos) + abs(cur_nsa - next_nsa)
                        !          gp_neighsa(lcurpos,new_nbidx) = cur_nsa
                        !       else
                        !          gp_sa(lcurpar) = gp_sa(lcurpar) + abs(cur_nsa - next_nsa)
                        !          gp_neighsa(lcurpar,nbidx) = next_nsa
                        !       endif
                        !    enddo ! end of nbIDX
                        ! endif ! end of gp_sa /= 0.0
                     endif ! end of not solve_ros
                  ENDDO
               ENDIF ! end of is_constant_density

               IF(numPEs > 1) THEN
                  CALL GLOBAL_ALL_SUM(glueVolume)
                  CALL GLOBAL_ALL_SUM(glueInertia)
               ENDIF

               glueDiameter = (6.0*glueVolume/PI) ** (1.0/3.0)

               deallocate(new_total_mass)
            ENDIF
         ELSE
            IF(ANY_SOLIDS_SPECIES_EQ) CALL ZERO_RRATE_DES
            IF(ENERGY_EQ) CALL ZERO_ENERGY_SOURCE
         ENDIF
         CALL CALC_PG_GRAD
      ENDIF

      IF(any(CALC_RADT_DES)) CALL CALC_avgTs

      !Hari Sitaraman	(particle sorting)===================
      if (optflag1.eq.1) then
        print *,"call spatial sort"
        CALL DES_SORT_PARTICLES_SPATIALLY()
      endif
      !======================================================

   END SUBROUTINE DES_TIME_INIT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Subroutine: DES_TIME_STEP                                        !
!     Author: Jay Boyalakuntla                        Date: 21-Jun-04  !
!                                                                      !
!     Purpose: Main DEM driver routine                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE DES_TIME_STEP(NN)

! Modules
      use compar, only: pe_io, mype, numPEs
!---------------------------------------------------------------------//
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NN
      logical :: mod_assertion
      INTEGER :: LC
! Wall time at the start of IO operations.
      DOUBLE PRECISION :: WALL_START_IO

      mod_assertion = (NN == 1 .OR. MOD(NN,NEIGHBOR_SEARCH_N) == 0)

         IF(DES_CONTINUUM_COUPLED) THEN
! If the current time in the discrete loop exceeds the current time in
! the continuum simulation, exit the discrete loop
            IF(S_TIME.GE.(TIME+DT_TMP)) THEN
               EXIT_LOOP = .TRUE.
               RETURN
            ENDIF
! If next time step in the discrete loop will exceed the current time
! in the continuum simulation, modify the discrete time step so final
! time will match
            IF((S_TIME+DTSOLID).GT.(TIME+DT_TMP)) &
               DTSOLID = TIME + DT_TMP - S_TIME
         ENDIF

 ! NREL_CU_OPT sorting,swappping routines
         if (optflag1.eq.1) then
            CALL DES_SORT_PARTICLE_ARRAYS
            CALL DES_GETVALID_FLUID_CELLS
         endif

         IF(GSP_EXPLICIT .and. numPEs > 1) THEN
            IF(DEM_BCMI > 0 .OR. DEM_BCMO > 0 .OR. N_GSP_DELETED>0) THEN
               CALL UPDATE_GSP_PARTICLE_STATE
               N_GSP_DELETED = 0
            ENDIF
         ENDIF

! Verlet velocity update for Glued particle DEM
         IF(GSP_EXPLICIT) THEN
            CALL GSP_CFNEWVALUES_vhalfxfull
         ELSEIF(GSP_IMPLICIT) THEN
            CALL GSP_CFNEWVALUES_vhalfxfull_implicit
         ENDIF

! Calculate inter particle forces acting (collisional, cohesion)
         IF(SuperDEM) THEN
            CALL CALC_FORCE_SuperDEM
         ELSE
            CALL CALC_FORCE_DEM
         ENDIF

! Update coordination number
         CALL CALC_COORDINATION

! Calculate or distribute fluid-particle drag force.
! gsp_implicit also includes drag force here
         IF (.NOT. GSP_EXPLICIT) CALL CALC_DRAG_DES

! Calculate heat conduction to/from wall
         ! accurate heat conduction is not expected in the gsp implicit
         IF(ENERGY_EQ) CALL CALC_DEM_THERMO_WITH_WALL_STL

! Update the old values of particle position and velocity with the new
! values computed
         IF (DO_OLD) CALL CFUPDATEOLD
! Calculate thermochemical sources (energy and  rates of formation).
         CALL CALC_THERMO_DES
! Update keyframe data
         CALL UPDATE_KEYFRAME_DATA
! Move internal surfaces
         CALL MOVE_IS_STL
! Call user functions.
         IF(CALL_USR) CALL USR1_DES
! Call component spheres conduction, set conduction kp in GSP_HEAT_COND_EXPLICIT
         IF(GSP_EXPLICIT .AND. ENERGY_EQ) CALL GSP_HEAT_COND_EXPLICIT
! Update velocity of ONLY EXITING particles associated with
! Non-MASS_OUTFLOW boundaries, when a flag criterion is met
         IF((DEM_BCMO > 0) .AND.         &
            (.NOT. FREEZE_EXIT_DES_VEL)) CALL UPDATE_EXIT_DES_VEL

! Update position and velocities
         IF(SuperDEM) THEN
            CALL SuperDEM_CFNEWVALUES
         ELSEIF (GSP_EXPLICIT) THEN
            CALL GSP_CFNEWVALUES_vfull
! reset FC for GSP explicit model !
            FC(:,1) = 0.d0
            FC(:,2) = 0.d0
            FC(:,3) = 0.d0
            TOW(:,1) = 0.d0
            TOW(:,2) = 0.d0
            TOW(:,3) = 0.d0
         ELSEIF(GSP_IMPLICIT) THEN
            CALL GSP_CFNEWVALUES_vfull_implicit
         ELSE
            CALL CFNEWVALUES
         ENDIF

! Update particle temperatures
         CALL DES_THERMO_NEWVALUES

         DO_NSEARCH = (NN == 1 .OR. MOD(NN,NEIGHBOR_SEARCH_N) == 0)
         if (DO_NSEARCH .neqv. mod_assertion) THEN
            WRITE(ERR_MSG, *) "Failed assertion; problem in neighbor search algorithm"
            CALL LOG_ERROR()
         ENDIF
! Add/Remove particles to the system via flow BCs.
         IF(DEM_BCMI > 0) CALL MASS_INFLOW_DEM(DTSOLID_TMP)
         IF(DEM_BCMO > 0) CALL MASS_OUTFLOW_DEM(DO_NSEARCH)

! Call exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
         IF (DO_NSEARCH .OR. (numPEs>1) .OR. DES_PERIODIC_WALLS) THEN
            CALL DESGRID_PIC(.TRUE.)
            CALL DES_PAR_EXCHANGE
         ENDIF

         IF(DO_NSEARCH) CALL NEIGHBOUR

! Explicitly coupled simulations do not need to rebin particles to
! the fluid grid every time step. However, this implies that the
! fluid cell information and interpolation weights become stale.
         IF(DES_CONTINUUM_COUPLED .AND. &
            .NOT.DES_EXPLICITLY_COUPLED) THEN
! Bin particles to fluid grid.
            CALL PARTICLES_IN_CELL
! Calculate interpolation weights
            CALL CALC_INTERP_WEIGHTS
! Calculate mean fields (EPg).
            CALL COMP_MEAN_FIELDS
         ENDIF

! Update time to reflect changes
         S_TIME = S_TIME + DTSOLID

      IF(WRITE_VTK_FILES.AND.DISCRETE_ELEMENT) THEN
         DO LC = 1, DIMENSION_VTK
            IF(VTK_DATA(LC)/='C'.AND.VTK_DATA(LC)/='G') THEN
               IF(S_TIME+0.1d0*DTSOLID>=VTK_TIME(LC)) THEN
                  VTK_TIME(LC) = (INT((S_TIME + 0.1d0*DTSOLID)/VTK_DT(LC))+1)*VTK_DT(LC)
                  CALL WRITE_VTP_FILE(LC,0,S_TIME)
               ENDIF
            ENDIF
         ENDDO
      ENDIF

! The following section targets data writes for DEM only cases:
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN
! Keep track of TIME and number of steps for DEM simulations
            TIME = S_TIME
            NSTEP = NSTEP + 1
! Call the output manager to write RES and SPx data.
            DLB = .TRUE.
            CALL WRITE_OUTPUTS(.FALSE., .FALSE.)
         ENDIF  ! end if (.not.des_continuum_coupled)

         IF(CALL_USR) CALL USR2_DES

         IF(ADJUST_PARTITION) THEN
            EXIT_LOOP = .TRUE.
            RETURN
         ENDIF

   END SUBROUTINE DES_TIME_STEP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Subroutine: DES_TIME_END                                         !
!     Author: Jay Boyalakuntla                        Date: 21-Jun-04  !
!                                                                      !
!     Purpose: Main DEM driver routine                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE DES_TIME_END

      IMPLICIT NONE

      IF(CALL_USR) CALL USR3_DES

! Reset the discrete time step to original value.
      DTSOLID = DTSOLID_TMP

      IF(DES_CONTINUUM_COUPLED) CALL DESMPI_SEND_RECV_FIELD_VARS
#ifndef QUIET
      TIME_TMP = GET_ELAPSED_TIME() - TIME_TMP
      IF(TIME_TMP > 1.0d-10) THEN
          WRITE(ERR_MSG, 9000) trim(iVal(dble(FACTOR)/TIME_TMP))
      ELSE
          WRITE(ERR_MSG, 9000) '+Inf'
      ENDIF
      CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)
 9000 FORMAT('    NITs/SEC = ',A)
#endif
      RETURN

   END SUBROUTINE DES_TIME_END

END MODULE DES_TIME_MARCH
