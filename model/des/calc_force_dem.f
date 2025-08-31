MODULE CALC_FORCE_DEM_MOD

   USE CFRELVEL_MOD, ONLY: CFRELVEL

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_FORCE_DEM                                          !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Calculate contact force and torque on particle from        !
!           particle-particle and particle-wall collisions. Treats     !
!           wall interaction also as a two-particle interaction but    !
!           accounting for the wall properties                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE CALC_FORCE_DEM

! Modules
!---------------------------------------------------------------------//
      USE calc_collision_wall
      use compar, only: pe_io, mype, numPEs
      USE constant, ONLY: Pi
      use desmpi_wrapper, only: DES_MPI_STOP,des_mpi_barrier
      USE des_thermo
      USE des_thermo_cond
      USE discretelement
      USE resize
      USE run
      use param1, only: one, small_number, zero, undefined
      use mpi_utility, only: global_all_sum

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! percent of particle radius when excess overlap will be flagged
      DOUBLE PRECISION, PARAMETER :: flag_overlap = 0.20d0 ! what is excess overlap?
! particle no. indices
      INTEGER :: I, LL, cc, CC_START, CC_END
! the overlap occurring between particle-particle or particle-wall
! collision in the normal direction
      DOUBLE PRECISION :: OVERLAP_N, OVERLAP_T(3)
! square root of the overlap
      DOUBLE PRECISION :: SQRT_OVERLAP
! distance vector between two particle centers or between a particle
! center and wall when the two surfaces are just at contact (i.e. no
! overlap)
      DOUBLE PRECISION :: R_LM,DIST_CI,DIST_CL
! the normal and tangential components of the translational relative
! velocity
      DOUBLE PRECISION :: V_REL_TRANS_NORM, rad
! distance vector between two particle centers or between a particle
! center and wall at current and previous time steps
      DOUBLE PRECISION :: DIST(3), NORMAL(3), DIST_MAG, POS(3)
! tangent to the plane of contact at current time step
      DOUBLE PRECISION :: VREL_T(3)
! normal and tangential forces
      DOUBLE PRECISION :: FN(3), FT(3)
! temporary storage of force
      DOUBLE PRECISION :: FC_TMP(3)
! temporary storage of force for torque
      DOUBLE PRECISION :: TOW_FORCE(3)
! temporary storage of torque
      DOUBLE PRECISION :: TOW_TMP(3,2)
! temporary storage of conduction/radiation
      DOUBLE PRECISION :: QQ_TMP

! store solids phase index of particle (i.e. pijk(np,5))
      INTEGER :: PHASEI, PHASELL
! local values used spring constants and damping coefficients
      DOUBLE PRECISION :: ETAN_DES, ETAT_DES
      DOUBLE PRECISION :: KN_DES, KT_DES
! local values used for calculating cohesive forces
      DOUBLE PRECISION :: FORCE_COH, EQ_RADIUS, DistApart

      LOGICAL, PARAMETER :: report_excess_overlap = .FALSE.

      DOUBLE PRECISION :: FNMD, FTMD, MAG_OVERLAP_T, TANGENT(3)
      integer :: start_index, end_index, do_index

!{SBAN
! Effective particle radius and mass to recalculate collision parameters
      DOUBLE PRECISION :: rad_eff, m_eff
!}

! Glued particle DEM model
! the glued particle id which sphere LL and sphere I belong to, LL-I is a cntact sphere pair
      INTEGER :: gidI, gidLL, gid, io_err, ppdim1, ppdim2
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ppforce !(NGluedParticles, 3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: pptorque !(NGluedParticles, 3)
! contact center, contact center to glued particle LL belonging to, contact center to glued particle I belonging to
! translational relative velocity
      DOUBLE PRECISION :: Contact_PT(3),dist2centerVECLL(3),dist2centerVECI(3),VRELTRANS(3)

! Force chain data
      INTEGER :: FCHAINC_OLD           ! old value of FCHAINC
      INTEGER :: old_size,new_size ! New force chain array size when resizing (growing) arrays
      DOUBLE PRECISION :: FC_MAG_MAX, MAG_ORIENT
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Coh_force

      INTEGER :: FC_MAX_INDEX         ! FC index where max occurs (within LL-I pair)

! Rolling friction
      DOUBLE PRECISION :: RFTOW(3), omega_mag
!------------------calculations start here!--------------------------------------------------->

! Initialize cohesive forces
      IF(USE_COHESION) PostCohesive(:) = ZERO

! calculate particle wall contact
      CALL CALC_DEM_FORCE_WITH_WALL_STL

! at the beginning each time step, reset those arrays to zero
      IF (GSP_EXPLICIT) THEN
         allocate(ppforce(NGluedParticles,3)); ppforce(:,:) = ZERO
         allocate(pptorque(NGluedParticles,3)); pptorque(:,:) = ZERO
         des_col_force(:,:) = ZERO
      ENDIF

      if (optflag1.eq.1) then
        start_index = 1
        end_index = MAX_PIP_EXIST
      else
        start_index = 1
        end_index = MAX_PIP
      endif


      IF(USE_COHESION) THEN
         allocate(Coh_force(MAX_PIP,3))
         Coh_force(:,:) = ZERO
      ENDIF

! Check particle LL neighbor contacts
!---------------------------------------------------------------------//

!$omp parallel default(none) private(pos,rad,cc,cc_start,cc_end,ll,i,  &
!$omp    overlap_n,vrel_t,v_rel_trans_norm,sqrt_overlap,dist,r_lm,     &
!$omp    kn_des,kt_des,phasell,phasei,etan_des,rad_eff,m_eff,          &
!$omp    gidLL,gidi,contact_pt,dist2centervecll,dist2centerveci,       &
!$omp    VRELTRANS,                                                    &
!$omp    etat_des,fn,ft,overlap_t,tangent,mag_overlap_t,               &
!$omp    eq_radius,distapart,force_coh,dist_mag,NORMAL,ftmd,fnmd,      &
!$omp    dist_cl, dist_ci, fc_tmp, tow_tmp, tow_force, qq_tmp,         &
!$omp    rftow,omega_mag )                                             &
!$omp shared(max_pip,neighbors,neighbor_index,des_pos_new,des_radius,  &
!$omp    omega_new, mew_r,                                             &
!$omp    des_coll_model_enum,kn,kt,pft_neighbor,pijk,hert_kn,hert_kt,  &
!$omp    des_etan,des_etat,mew,use_cohesion, calc_cond_des, dtsolid,   &
!$omp    van_der_waals,vdw_outer_cutoff,vdw_inner_cutoff,              &
!$omp    hamaker_constant,asperities,surface_energy, optflag1,         &
!$omp    start_index, end_index,rk_global_cnt_debug,                   &
!$omp    des_usr_var,GluedSphereDEM,gluepos,des_vel_new,glueomg,       &
!$omp    gsp_implicit, gsp_explicit,                                   &
!$omp    ppforce,pptorque,gid_list,                                    &
!$omp    tow, fc, energy_eq, grav_mag, postcohesive, pmass, q_source,  &
!$omp    write_force_chain, des_col_force, iglobal_id,particle_state,  &
!$omp    old_size,new_size,io_err,                                     &
!$omp    fchain_fcmax,fc_mag_max,fc_max_index, mag_orient,coh_force,   &
!$omp    FCHAINC,FCHAINC_OLD,                                          &
!$omp    FCHAIN_MIDPOINT,FCHAIN_CONTACT_POINT,                         &
!$omp    FCHAIN_ORIENT, FCHAIN_FN,FCHAIN_FT,FCHAIN_FC,                 &
!$omp    FCHAIN_LENGTH,FCHAIN_FN_MAG,FCHAIN_LOCAL_ID1,FCHAIN_LOCAL_ID2,&
!$omp    FCHAIN_GLOBAL_ID1,FCHAIN_GLOBAL_ID2, FCHAIN_OVERLAP) &
!$omp shared(max_pip_exist)

! Force chain initialization
      IF(WRITE_FORCE_CHAIN) THEN
         FCHAINC_OLD = min(FCHAINC, size(FCHAIN_MIDPOINT,1))
         FCHAINC = 0
         FCHAIN_MIDPOINT(:,:) = UNDEFINED
         FCHAIN_CONTACT_POINT(:,:) = UNDEFINED
         FCHAIN_ORIENT(:,:) = UNDEFINED
         FCHAIN_FCMAX(:) = .FALSE.
      ENDIF

!$omp do reduction (+: FCHAINC)
      DO do_index = start_index, end_index
         if (optflag1.eq.1) then
            ll = do_index
         else
            ll = do_index
            IF(IS_NONEXISTENT(LL)) CYCLE
         endif

         IF(WRITE_FORCE_CHAIN) THEN
            FC_MAG_MAX = ZERO
            FC_MAX_INDEX = -1
         ENDIF


         POS = DES_POS_NEW(LL,:)
         rad = DES_RADIUS(LL)

! read the gid of the current sphere LL!
         IF (GSP_EXPLICIT) gidLL = gid_list(LL)
!-------------------------------------------------------------------------------------
! NEIGHBORS: store the neighbor component spheres in order
! NEIGHBOR_INDEX: the range index of neighbor component sphere
!-------------------------------------------------------------------------------------
         CC_START = 1
         IF (LL.gt.1) CC_START = NEIGHBOR_INDEX(LL-1)
         CC_END   = NEIGHBOR_INDEX(LL)

         DO CC = CC_START, CC_END-1
! read the gid of the current sphere I!
            I  = NEIGHBORS(CC)
            IF(GSP_EXPLICIT) gidI = gid_list(I)
            IF(IS_NONEXISTENT(I)) CYCLE

! distance between LL and I center just about to contact, no overlap
            R_LM = rad + DES_RADIUS(I)

! distance between current and previous time steps
            DIST(:) = DES_POS_NEW(I,:) - POS(:)
            DIST_MAG = dot_product(DIST,DIST)

! reset temporary force of collision
            FC_TMP(:) = ZERO

! Compute particle-particle VDW cohesive short-range forces
            IF(USE_COHESION .AND. VAN_DER_WAALS) THEN
               IF(DIST_MAG < (R_LM+VDW_OUTER_CUTOFF)**2) THEN
                  EQ_RADIUS = 2d0 * DES_RADIUS(LL)*DES_RADIUS(I) /     &
                    (DES_RADIUS(LL)+DES_RADIUS(I))
                  IF(DIST_MAG > (VDW_INNER_CUTOFF+R_LM)**2) THEN
                     DistApart = (SQRT(DIST_MAG)-R_LM)
                     FORCE_COH = HAMAKER_CONSTANT * EQ_RADIUS /           &
                        (12d0*DistApart**2) * (Asperities/(Asperities+    &
                        EQ_RADIUS) + ONE/(ONE+Asperities/DistApart)**2 )
                  ELSE
                     FORCE_COH = 2d0 * PI * SURFACE_ENERGY * EQ_RADIUS *  &
                       (Asperities/(Asperities+EQ_RADIUS) + ONE/          &
                       (ONE+Asperities/VDW_INNER_CUTOFF)**2 )
                  ENDIF
                  FC_TMP(:) = DIST(:)*FORCE_COH/SQRT(DIST_MAG)
                  Coh_Force(LL,:) = Coh_Force(LL,:) + FC_TMP(:)
                  Coh_Force(I,:)  = Coh_Force(I ,:) - FC_TMP(:)
                  TOW_TMP(:,:) = ZERO

               ENDIF
            ENDIF

            IF (ENERGY_EQ) THEN
               IF(gsp_implicit) THEN
                  ! No action needed, project2 method does not support for now
                  CONTINUE
               ELSE
! Calculate conduction and radiation for thermodynamic neighbors
                  IF(CALC_COND_DES(PIJK(LL,5))) THEN
                     QQ_TMP = DES_CONDUCTION(LL, I, sqrt(DIST_MAG), PIJK(LL,5), PIJK(LL,4))
!$omp atomic
                     Q_Source(LL) = Q_Source(LL) + QQ_TMP
!$omp atomic
                     Q_Source(I) = Q_Source(I) - QQ_TMP
                  ENDIF
               ENDIF
            ENDIF

            IF(DIST_MAG > (R_LM - SMALL_NUMBER)**2) THEN
               PFT_NEIGHBOR(:,CC) = 0.0
               CYCLE
            ENDIF

            IF(DIST_MAG == 0) THEN
               WRITE(*,8550) LL, I
               ERROR STOP "division by zero"
 8550 FORMAT('distance between particles is zero:',2(2x,I10))
            ENDIF

            DIST_MAG = SQRT(DIST_MAG)
            NORMAL(:)= DIST(:)/DIST_MAG
! Calculate the normal overlap
            OVERLAP_N = R_LM-DIST_MAG

            IF(REPORT_EXCESS_OVERLAP) CALL PRINT_EXCESS_OVERLAP

            IF (GSP_EXPLICIT) THEN
! translational relative velocity for glued particles
               Contact_PT = POS(:) +  NORMAL(:)*(rad-OVERLAP_N/2)
               dist2centerVECLL = Contact_PT(:) - gluePos(gidLL,:)
               dist2centerVECI  = Contact_PT(:) - gluePos(gidI,:)
               VRELTRANS(:) = (DES_VEL_NEW(LL,:) - DES_VEL_NEW(I,:))
               VRELTRANS(:) =  VRELTRANS(:) + CROSS(glueOmg(gidLL,:),dist2centerVECLL)
               VRELTRANS(:) =  VRELTRANS(:) - CROSS(glueOmg(gidI,:),dist2centerVECI)
               V_REL_TRANS_NORM = DOT_PRODUCT(VRELTRANS,NORMAL)
               VREL_T(:) =  VRELTRANS(:) - V_REL_TRANS_NORM*NORMAL(:)
            ELSEIF (GSP_IMPLICIT) THEN
               CALL gsp_gsp_collison(LL, I, CC, V_REL_TRANS_NORM, VREL_T, NORMAL, OVERLAP_N, &
                                    dist2centerVECI, dist2centerVECLL)
               IF(OVERLAP_N == 0.0) CYCLE
            ELSE
! Calculate the components of translational relative velocity for a
! contacting particle pair and the tangent to the plane of contact
               CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, VREL_T, NORMAL(:), DIST_MAG)
            ENDIF
            phaseLL = PIJK(LL,5)
            phaseI = PIJK(I,5)
!{SBAN
            rad_eff = (rad * DES_RADIUS(I)) / (rad + DES_RADIUS(I))
            m_eff = (PMASS(LL) * PMASS(I)) / (PMASS(LL) + PMASS(I))
! Hertz spring-dashpot contact model
            IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
               sqrt_overlap = SQRT(OVERLAP_N)
               KN_DES = hert_kn(phaseLL,phaseI)*sqrt_overlap * sqrt(rad_eff)
               KT_DES = hert_kt(phaseLL,phaseI)*sqrt_overlap * sqrt(rad_eff)
               sqrt_overlap = SQRT(sqrt_overlap)
               ETAN_DES = DES_ETAN(phaseLL,phaseI)*sqrt_overlap * sqrt(hert_kn(phaseLL,phaseI)*sqrt(rad_eff)*m_eff)
               ETAT_DES = DES_ETAT(phaseLL,phaseI)*sqrt_overlap * sqrt(hert_kt(phaseLL,phaseI)*sqrt(rad_eff)*m_eff)
! Linear spring-dashpot contact model
            ELSE
               KN_DES = KN
               KT_DES = KT
               ETAN_DES = DES_ETAN(phaseLL,phaseI) * sqrt(m_eff)
               ETAT_DES = DES_ETAT(phaseLL,phaseI) * sqrt(m_eff)
!}
            ENDIF

! Calculate the normal contact force
            FN(:) =  -(KN_DES * OVERLAP_N * NORMAL(:) + &
               ETAN_DES * V_REL_TRANS_NORM * NORMAL(:))

! Calculate the tangential overlap
            OVERLAP_T(:) = DTSOLID*VREL_T(:) + PFT_NEIGHBOR(:,CC)
            MAG_OVERLAP_T = sqrt(dot_product(OVERLAP_T,OVERLAP_T))

! Calculate the tangential contact force.
            IF(MAG_OVERLAP_T > 0.0) THEN
! Calculate the tangential contact force.
               FT = -KT_DES*OVERLAP_T - ETAT_DES*VREL_T
               FTMD = sqrt(dot_product(FT,FT))
! Max force before the on set of frictional slip.
               FNMD = MEW*sqrt(dot_product(FN,FN))
! Frictional slip
               IF(FTMD > FNMD) THEN
! Direction of tangential force.
                  TANGENT = OVERLAP_T/MAG_OVERLAP_T
                  FT = -FNMD * TANGENT
                  OVERLAP_T = (FNMD/KT_DES) * TANGENT
               ENDIF
            ELSE
               FT = 0.0
            ENDIF

! Save tangential displacement history
            PFT_NEIGHBOR(:,CC) = OVERLAP_T(:)
            IF (GSP_EXPLICIT .OR. GSP_IMPLICIT) THEN
               TOW_TMP(:,1) = CROSS(dist2centerVECLL,FN) + CROSS(dist2centerVECLL,FT)
               TOW_TMP(:,2) = CROSS(dist2centerVECI ,-FN) + CROSS(dist2centerVECI ,-FT)
            ELSE
! calculate the distance from the particles' centers to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
               DIST_CL = DIST_MAG/2.d0 + (DES_RADIUS(LL)**2 - &
                  DES_RADIUS(I)**2)/(2.d0*DIST_MAG)

               DIST_CI = DIST_MAG - DIST_CL

               TOW_force(:) = CROSS(NORMAL(:), FT(:))
               TOW_TMP(:,1) = DIST_CL*TOW_force(:)
               TOW_TMP(:,2) = DIST_CI*TOW_force(:)
            ENDIF
! Calculate the total force FC of a collision pair
! total contact force ( FC_TMP may already include cohesive force)
            FC_TMP(:) = FC_TMP(:) + FN(:) + FT(:)

! force and torque summation, glued sphere can use the same format of summation

!$omp atomic
            FC(LL,1) = FC(LL,1) + FC_TMP(1)
!$omp atomic
            FC(LL,2) = FC(LL,2) + FC_TMP(2)
!$omp atomic
            FC(LL,3) = FC(LL,3) + FC_TMP(3)
!$omp atomic
            FC(I,1) = FC(I,1) - FC_TMP(1)
!$omp atomic
            FC(I,2) = FC(I,2) - FC_TMP(2)
!$omp atomic
            FC(I,3) = FC(I,3) - FC_TMP(3)

            IF(GSP_EXPLICIT .AND. (particle_state(LL) .EQ. 1)) THEN
!$omp atomic
               ppforce(gidLL,1) = ppforce(gidLL,1) + FC_TMP(1)
!$omp atomic
               ppforce(gidLL,2) = ppforce(gidLL,2) + FC_TMP(2)
!$omp atomic
               ppforce(gidLL,3) = ppforce(gidLL,3) + FC_TMP(3)
!$omp atomic
               ppforce(gidI,1) = ppforce(gidI,1) - FC_TMP(1)
!$omp atomic
               ppforce(gidI,2) = ppforce(gidI,2) - FC_TMP(2)
!$omp atomic
               ppforce(gidI,3) = ppforce(gidI,3) - FC_TMP(3)
            ENDIF

! for each particle the signs of norm and ft both flip, so add the same torque

!$omp atomic
            TOW(LL,1) = TOW(LL,1) + TOW_TMP(1,1)
!$omp atomic
            TOW(LL,2) = TOW(LL,2) + TOW_TMP(2,1)
!$omp atomic
            TOW(LL,3) = TOW(LL,3) + TOW_TMP(3,1)
!$omp atomic
            TOW(I,1)  = TOW(I,1)  + TOW_TMP(1,2)
!$omp atomic
            TOW(I,2)  = TOW(I,2)  + TOW_TMP(2,2)
!$omp atomic
            TOW(I,3)  = TOW(I,3)  + TOW_TMP(3,2)

            IF(GSP_EXPLICIT .AND. (particle_state(LL) .EQ. 1)) THEN
!$omp atomic
               pptorque(gidLL,1) = pptorque(gidLL,1) + TOW_TMP(1,1)
!$omp atomic
               pptorque(gidLL,2) = pptorque(gidLL,2) + TOW_TMP(2,1)
!$omp atomic
               pptorque(gidLL,3) = pptorque(gidLL,3) + TOW_TMP(3,1)
!$omp atomic
               pptorque(gidI,1) = pptorque(gidI,1) + TOW_TMP(1,2)
!$omp atomic
               pptorque(gidI,2) = pptorque(gidI,2) + TOW_TMP(2,2)
!$omp atomic
               pptorque(gidI,3) = pptorque(gidI,3) + TOW_TMP(3,2)

            ENDIF

            DES_COL_FORCE(I,:) = FC(I,:)

!######################################################################
! Force chain data - Implemented by Eric Breard (U. of Oregon)
!                    Reviewed by Jeff Dietiker
!######################################################################
! Note: Force Chain Data does not work in SMP
!######################################################################

            IF(WRITE_FORCE_CHAIN) THEN
! Resize arrays if needed
               IF(FCHAINC>=FCHAINC_OLD.AND.MOD(FCHAINC,100)==0) THEN
                  old_size = size(FCHAIN_MIDPOINT,1)
                  IF(FCHAINC>=INT(0.9*old_size)) THEN
                     new_size = INT(1.5*old_size)
                     call real_grow2_reverse(FCHAIN_MIDPOINT,new_size)
                     FCHAIN_MIDPOINT(old_size+1:new_size,:) = UNDEFINED
                     call real_grow2_reverse(FCHAIN_CONTACT_POINT,new_size)
                     FCHAIN_CONTACT_POINT(old_size+1:new_size,:) = UNDEFINED
                     call real_grow2_reverse(FCHAIN_ORIENT,new_size)
                     call real_grow2_reverse(FCHAIN_FN,new_size)
                     call real_grow2_reverse(FCHAIN_FT,new_size)
                     call real_grow2_reverse(FCHAIN_FC,new_size)
                     call real_grow(FCHAIN_LENGTH,new_size)
                     call real_grow(FCHAIN_FN_MAG,new_size)
                     call logical_grow(FCHAIN_FCMAX,new_size)
                     call integer_grow(FCHAIN_LOCAL_ID1,new_size)
                     call integer_grow(FCHAIN_LOCAL_ID2,new_size)
                     call integer_grow(FCHAIN_GLOBAL_ID1,new_size)
                     call integer_grow(FCHAIN_GLOBAL_ID2,new_size)
                     call real_grow(FCHAIN_OVERLAP,new_size)
                  ENDIF
               ENDIF
! Save data
               FCHAINC                     = FCHAINC + 1  ! Force chain counter
               FCHAIN_MIDPOINT(FCHAINC,:)  = 0.5D0 * (DES_POS_NEW(LL,:) + DES_POS_NEW(I,:)) ! Mid point coordinates
               FCHAIN_CONTACT_POINT(FCHAINC,:)  = DES_POS_NEW(LL,:) + NORMAL(:) * (DES_RADIUS(LL) + 0.5*OVERLAP_N) ! Contact point
               MAG_ORIENT                  = dsqrt(dot_product(FN,FN))
               if(MAG_ORIENT/=ZERO) FCHAIN_ORIENT(FCHAINC,:)    = FN(:)/MAG_ORIENT ! Orientation (same as Normal force)
               FCHAIN_FN(FCHAINC,:)        = FN(:)    ! Normal force
               FCHAIN_FT(FCHAINC,:)        = FT(:)    ! Tangential force
               FCHAIN_FC(FCHAINC,:)        = FC_TMP(:)! Total force (including cohesion)
               FCHAIN_LENGTH(FCHAINC)      = DIST_MAG !distance between LL and I
               FCHAIN_FN_MAG(FCHAINC)      = DSQRT(FN(1)**2 + FN(2)**2 + FN(3)**2) ! Normal force magnitude
               FCHAIN_LOCAL_ID1(FCHAINC)   = LL
               FCHAIN_LOCAL_ID2(FCHAINC)   = I
               FCHAIN_GLOBAL_ID1(FCHAINC)  = iglobal_id(LL)
               FCHAIN_GLOBAL_ID2(FCHAINC)  = iglobal_id(I)
               FCHAIN_OVERLAP(FCHAINC)     = OVERLAP_N
               IF(FCHAIN_FN_MAG(FCHAINC)>FC_MAG_MAX) THEN
                  FC_MAX_INDEX = FCHAINC  ! Keep track of max magnitude for pair LL-I
                  FC_MAG_MAX = FCHAIN_FN_MAG(FCHAINC)
               ENDIF

            ENDIF
!Rolling friction
! JFD: Here the rolling friction coefficient is non-dimensional. This is the equivalent to
! the Zhou paper (Physica A 269 (1999) 536-553) 's definition (Model A) divided by the
! particle diameter. Therefore here mew_r is multiplied by the diameter.

! TBD: implementation of rolling friction for glued particles @renjieke 1-1-2024
            IF (.not. GluedSphereDEM) THEN
            if(MEW_R>ZERO) then
               ! Particle LL
               OMEGA_MAG = sqrt(dot_product(OMEGA_NEW(LL,:),OMEGA_NEW(LL,:)))
               if(OMEGA_MAG > ZERO) then
                  RFTOW = -MEW_R*(2.0*DES_RADIUS(LL))*KN_DES * OVERLAP_N*OMEGA_NEW(LL,:)/OMEGA_MAG
                  TOW(LL,:) = TOW(LL,:) + RFTOW
               endif
               !particle i
               OMEGA_MAG = sqrt(dot_product(OMEGA_NEW(I,:),OMEGA_NEW(I,:)))
               if(OMEGA_MAG > ZERO) then
                  RFTOW = -MEW_R*(2.0*DES_RADIUS(I))*KN_DES * OVERLAP_N*OMEGA_NEW(I,:)/OMEGA_MAG
                  TOW(I,:) = TOW(I,:) + RFTOW
               endif
            endif
            ENDIF
         ENDDO ! Neighbor CC i.e. particle I loop
         DES_COL_FORCE(LL,:) = FC(LL,:)


         IF(WRITE_FORCE_CHAIN.AND.FC_MAX_INDEX>0) FCHAIN_FCMAX(FC_MAX_INDEX) = .TRUE. ! Record where max magnitude occurs (for a given particle LL
      ENDDO ! Particle LL loop
!$omp end do

!$omp end parallel

! just for post-processing mag. of cohesive forces on each particle
! Divide by weight or mass if gravity is turned off
      IF(USE_COHESION .AND. VAN_DER_WAALS) THEN
         IF(GRAV_MAG > ZERO) THEN
            DO LL = 1, MAX_PIP
               PostCohesive(LL) = dsqrt(Coh_Force(LL,1)**2+Coh_Force(LL,2)**2+Coh_Force(LL,3)**2)/(PMASS(LL)*GRAV_MAG)
            ENDDO
         ELSE
            DO LL = 1, MAX_PIP
               PostCohesive(LL) = dsqrt(Coh_Force(LL,1)**2+Coh_Force(LL,2)**2+Coh_Force(LL,3)**2)/PMASS(LL)
            ENDDO
         ENDIF
         deallocate(Coh_force)
      ENDIF

      IF(GSP_EXPLICIT) THEN
         DO gid = 1, NGluedParticles
            ! ppforce(gid) should be zero if this gsp does not contact, so remove if(glueColNum(gid)==0) cycle should be safe
            if(all(ppforce(gid,:) == 0.0)) CYCLE
            glueForce(gid,:)  = glueForce(gid,:)  + ppforce(gid,:)
            if(all(pptorque(gid,:) == 0.0)) CYCLE
            glueTorque(gid,:) = glueTorque(gid,:) + pptorque(gid,:)
         ENDDO
         IF(allocated(ppforce)) deallocate(ppforce)
         IF(allocated(pptorque)) deallocate(pptorque)
      ENDIF

      RETURN

      contains

        include 'functions.inc'




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: print_excess_overlap                                    !
!                                                                      !
!  Purpose: Print overlap warning messages.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PRINT_EXCESS_OVERLAP

      use error_manager

      IF(OVERLAP_N > flag_overlap*DES_RADIUS(LL) .OR.                  &
         OVERLAP_N > flag_overlap*DES_RADIUS(I)) THEN

         WRITE(ERR_MSG,1000) trim(iVAL(LL)), trim(iVAL(I)), S_TIME,    &
            DES_RADIUS(LL), DES_RADIUS(I), OVERLAP_N

         CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 1000 FORMAT('WARNING: Excessive overlap detected between ',          &
         'particles ',A,' and ',/A,' at time ',g11.4,'.',/             &
         'RADII:  ',g11.4,' and ',g11.4,4x,'OVERLAP: ',g11.4)

      END SUBROUTINE PRINT_EXCESS_OVERLAP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: gsp_gsp_collison                                        !
!                                                                      !
!  Purpose: Calculate the collision between glued particles            !
!           when internal property distributions are not considered.   !
!  Author: Hang Zhou                                  Date: April 2024 !
!  Revised: Renjie Ke                                 Date: 11/01/2024 !
!  Revised: Renjie Ke                                 Date: 03/15/2025 !
!           Force the contact normal in the implicit method the same   !
!           as in the explicit method; fix the MPI issue.              !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE gsp_gsp_collison(index_ll, index_i, index_CC, VRN, VSLIP, NORM, overlapn, dist2I, dist2LL)


      use GSP_MATH_MOD, only: gsp_quat_to_exyz
      use error_manager, only: err_msg

      IMPLICIT NONE

      INTEGER, INTENT(IN)::index_i, index_ll, index_CC !particle id
! slip velocity at point of contact
      DOUBLE PRECISION, INTENT(OUT) :: VSLIP(:)
! normal component of relative contact velocity (scalar)
      DOUBLE PRECISION, INTENT(OUT) :: VRN
! unit normal vector along the line of contact pointing from
! particle L to particle II
      DOUBLE PRECISION, INTENT(OUT) :: NORM(3)
! overlap between particles
      DOUBLE PRECISIOn, INTENT(OUT) :: OVERLAPN
! temporary storage of torque
      DOUBLE PRECISIOn, INTENT(OUT) :: dist2I(3), dist2LL(3)

      INTEGER :: ii, ll, mi, mll, ss
      INTEGER :: num_contact, n_rows
      INTEGER, DIMENSION (:), ALLOCATABLE :: si, sll, si_tmp, sll_tmp ! list of spheres
      DOUBLE PRECISION :: rel_omgrad(3)
      double precision :: startt

      double precision :: sphere_radius_i, sphere_radius_ll
      double precision :: sphere_pos_i(3), sphere_pos_ll(3)
! bounding diameter of particle I and LL, they are not always d_p0(phase)
      double precision :: real_b_ii, real_b_ll, xloc, yloc, zloc, r_sum
      double precision :: exone_i(3),eyone_i(3),ezone_i(3)
      double precision :: exone_ll(3),eyone_ll(3),ezone_ll(3)
      double precision :: distll2i(3),distll2i_mag,dist_tmp(3)
      double precision :: lli_contact_point(3)
      double precision :: rvel(3)

      mi = PIJK(index_i, 5) ! phase id of particle i
      mll = PIJK(index_ll, 5) ! phase id of partice j

! n_rows should be the same for all phase
      n_rows = size(cspart_data,1)


      NORM(:) = (/0.D0,0.D0,0.D0/)
      lli_contact_point(:) = (/0.D0,0.D0,0.D0/)
      OVERLAPN = 0.d0
      num_contact = 0

      real_b_ii = 2 * des_radius(index_i)
      real_b_ll = 2 * des_radius(index_ll)

      exone_i=(/1.0D0, 0.0D0, 0.0D0/)
      eyone_i=(/0.0D0, 1.0D0, 0.0D0/)
      ezone_i=(/0.0D0, 0.0D0, 1.0D0/)
      CALL gsp_quat_to_exyz(glueQuat(index_i,:),exone_i,eyone_i,ezone_i)

      exone_ll=(/1.0D0, 0.0D0, 0.0D0/)
      eyone_ll=(/0.0D0, 1.0D0, 0.0D0/)
      ezone_ll=(/0.0D0, 0.0D0, 1.0D0/)
      CALL gsp_quat_to_exyz(glueQuat(index_ll,:),exone_ll,eyone_ll,ezone_ll)
! loop over spheres in gsp index_ll
   ! loop over spheres in gsp index_i
      do ll=1, n_rows
         if(cspart_data(ll,4,mll) == 0.0) exit
         if(cspart_data(ll,5,mll) == 0.0) cycle

         sphere_radius_ll = cspart_data(ll,4,mll) * real_b_ll / 2.0
         xloc = cspart_data(ll,1,mll) * real_b_ll
         yloc = cspart_data(ll,2,mll) * real_b_ll
         zloc = cspart_data(ll,3,mll) * real_b_ll

         sphere_pos_ll(:) = des_pos_new(INDEX_LL,:) + xloc * exone_ll &
                  + yloc * eyone_ll + zloc * ezone_ll

         do ii = 1, n_rows
            if(cspart_data(ii,4,mi) == 0.0) exit
            if(cspart_data(ii,5,mi) == 0.0) cycle

            sphere_radius_i = cspart_data(ii,4,mi) * real_b_ii / 2.0
            xloc = cspart_data(ii,1,mi) * real_b_ii
            yloc = cspart_data(ii,2,mi) * real_b_ii
            zloc = cspart_data(ii,3,mi) * real_b_ii

            sphere_pos_i(:) = des_pos_new(INDEX_I,:) + xloc * exone_i &
                     + yloc * eyone_i + zloc * ezone_i

            r_sum = sphere_radius_i + sphere_radius_ll
            ! distll2i vector points to i
            distll2i(:) =  sphere_pos_i(:) - sphere_pos_ll(:)
            distll2i_mag = dot_product(distll2i,distll2i)
            ! if two components do not contact
            IF(distll2i_mag > (r_sum - SMALL_NUMBER)**2) CYCLE

            IF (ENERGY_EQ) THEN
            ! Calculate conduction and radiation for thermodynamic neighbors
               IF(CALC_COND_DES(mll)) THEN
                  QQ_TMP = DES_CONDUCTION(index_ll, index_i, sqrt(distll2i_mag), PIJK(index_ll,5), PIJK(index_ll,4), &
                     sphere_radius_i, sphere_radius_ll)
!$omp atomic
                  Q_Source(index_ll) = Q_Source(index_ll) + QQ_TMP
!$omp atomic
                  Q_Source(index_i) = Q_Source(index_i) - QQ_TMP
               ENDIF
            ENDIF

            distll2i_mag = sqrt(distll2i_mag)
            NORM(:) = NORM(:) + distll2i(:)/distll2i_mag
            OVERLAPN = OVERLAPN + r_sum - distll2i_mag
            ! the location of contact point should also be examed
            lli_contact_point(:)  = lli_contact_point(:) + sphere_pos_ll &
               + (sphere_radius_ll - (r_sum - distll2i_mag)/2)*distll2i(:)/distll2i_mag
            num_contact = num_contact + 1
         enddo ! end of ii
      enddo ! end of ll

      if(num_contact==0) then
         PFT_NEIGHBOR(:,index_CC) = 0.0
         return

      endif
      OVERLAPN = OVERLAPN/num_contact
      ! norm direction needs more carefully consideration!
      norm(:) = norm(:)/num_contact
      lli_contact_point(:) = lli_contact_point(:)/num_contact
      r_sum = DES_RADIUS(index_i) + DES_RADIUS(index_ll)

      dist2LL = lli_contact_point- des_pos_new(index_ll,:)
      dist2I = lli_contact_point- des_pos_new(index_i,:)
      rel_omgrad = CROSS(OMEGA_NEW(index_ll,:),dist2ll) - CROSS(OMEGA_NEW(index_i,:),dist2i)
      rvel(:) = des_vel_new(index_ll,:) - des_vel_new(index_i,:) + rel_omgrad

      VRN = DOT_PRODUCT(rvel,norm)
      VSLIP(:) = rvel(:) - VRN*norm

      END SUBROUTINE gsp_gsp_collison

    END SUBROUTINE CALC_FORCE_DEM

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_COORDINATION                                       C
!  Purpose: Compute coordination number and overlap statistics         C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 17-SEP-24  C
!  Reviewer:                                                           C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_COORDINATION

! Modules
!---------------------------------------------------------------------//
      USE discretelement
      USE param
      USE param1
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: I, LID1, LID2
      double precision :: OVERLAP
!---------------------------------------------------------------------//

! Coordination number can be extracted from force chain data
! A force chain links two particles (Local IDs LID1 and LID2), so the
! coordination in incremented for each force chain. We can keep track of
! the min, max and mean overlap for each particle.

      IF(.NOT.WRITE_FORCE_CHAIN) RETURN

      COORDINATION(:) = ZERO
      MIN_OVERLAP(:) = UNDEFINED
      MAX_OVERLAP(:) = -UNDEFINED
      MEAN_OVERLAP(:) = ZERO


      DO I = 1, FCHAINC

         LID1 = FCHAIN_LOCAL_ID1(I)
         LID2 = FCHAIN_LOCAL_ID2(I)

         OVERLAP = FCHAIN_OVERLAP(I)

         COORDINATION(LID1) = COORDINATION(LID1) + 1
         COORDINATION(LID2) = COORDINATION(LID2) + 1

         MIN_OVERLAP(LID1) = DMIN1(MIN_OVERLAP(LID1),OVERLAP)
         MAX_OVERLAP(LID1) = DMAX1(MAX_OVERLAP(LID1),OVERLAP)
         MEAN_OVERLAP(LID1) = MEAN_OVERLAP(LID1) + OVERLAP

         MIN_OVERLAP(LID2) = DMIN1(MIN_OVERLAP(LID2),OVERLAP)
         MAX_OVERLAP(LID2) = DMAX1(MAX_OVERLAP(LID2),OVERLAP)
         MEAN_OVERLAP(LID2) = MEAN_OVERLAP(LID2) + OVERLAP

      ENDDO

      WHERE(COORDINATION>ZERO)
         MEAN_OVERLAP(:) = MEAN_OVERLAP(:) / COORDINATION(:)
      ELSEWHERE
         MIN_OVERLAP(:) = ZERO
         MAX_OVERLAP(:) = ZERO
      ENDWHERE

! SuperDEM force chain is made up of two segments, so we need to divide
! coordination by two (min, max, mean are the same)
      IF(SUPERDEM) COORDINATION(:) = COORDINATION(:) / 2.0D0

      RETURN

      END SUBROUTINE CALC_COORDINATION

 END MODULE CALC_FORCE_DEM_MOD
