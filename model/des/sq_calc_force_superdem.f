MODULE SQ_CALC_FORCE_MOD

   USE CFRELVEL_MOD, ONLY: CFRELVEL, Sq_CFRELVEL
   USE SQ_OBB_MOD
   USE SQ_CONTACT_Newton_DPmethod_MOD
   USE SQ_CONTACT_EVOLUTION_MOD
   USE SQ_EQUIVALENT_RADIUS_MOD
   USE SQ_ROTATION_MOD
   USE SQ_PROPERTIES_MOD
CONTAINS

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: calculation of force for SuperDEM                                 !
!                                                                             !
! Reference:                                                                  !
!   Xi Gao, Jia Yu, Ricardo JF Portal, Jean-Francois Dietiker,                !
!   Mehrdad Shahnam and William A Rogers, Development and validation          !
!   of SuperDEM for non-spherical particulate systems using a                 !
!   superquadric particle method, Particuology,                               !
!   2021: doi.org/10.1016/j.partic.2020.1011.1007.                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!


   SUBROUTINE CALC_FORCE_SuperDEM

! Modules
!---------------------------------------------------------------------//
      USE calc_collision_wall
      USE constant, ONLY: Pi
      USE des_thermo
      USE des_thermo_cond
      USE discretelement
      USE resize
      USE run
      use param1, only: one, small_number, zero, undefined
      use usr

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! percent of particle radius when excess overlap will be flagged
      DOUBLE PRECISION, PARAMETER :: flag_overlap = 0.20d0
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
      DOUBLE PRECISION :: R_LM,DIST_CI,DIST_CL,DIST_CII(3),DIST_CLL(3)
! the normal and tangential components of the translational relative
! velocity
      DOUBLE PRECISION :: V_REL_TRANS_NORM, rad
! distance vector between two particle centers or between a particle
! center and wall at current and previous time steps
      DOUBLE PRECISION :: DIST(3), NORMAL(3), DIST_MAG, POS(3),dist_li
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

! Implicit indexes of superquadric surfaces i and j (roundness)
      DOUBLE PRECISION :: m1,n1,m2,n2
! Euler parameters of superquadric surfaces, also called quaternions
      DOUBLE PRECISION :: euler0A,euler1A, euler2A, euler3A,tra(3,3)
      DOUBLE PRECISION :: euler0B,euler1B, euler2B, euler3B,trb(3,3)
      DOUBLE PRECISION :: Qc(4),PFT_NEIGHBOR_Global(3)
! P is the location of superquadric particle center, axi is the sem-axes
      DOUBLE PRECISION :: p(6), axi(6),X(6), Rn1n(3), Rn2n(3)
      DOUBLE PRECISION :: contacttest, objreport,OBB_TEST
      DOUBLE PRECISION :: ContactTest_A, contactTest_B
      DOUBLE PRECISION :: X0_A(3),X0_B(3),XCON_A(3), XCON_B(3)
      DOUBLE PRECISION :: LAMBDA0_A,LAMBDA0_B,LAMBDA_A, LAMBDA_B
      DOUBLE PRECISION :: Rn1n_A(3), Rn2n_A(3),Rn1n_B(3), Rn2n_B(3)
      DOUBLE PRECISION :: pa(3),pb(3),xa(3),xb(3),axia(3),axib(3)
      DOUBLE PRECISION :: eq_ra,eq_rb,Ixg,dp
      DOUBLE PRECISION :: HERT_KN_SUPER(DIM_M,DIM_M),HERT_KT_SUPER(DIM_M,DIM_M)
      DOUBLE PRECISION :: DES_ETAN_SUPER(DIM_M,DIM_M),DES_ETAT_SUPER(DIM_M,DIM_M)
      DOUBLE PRECISION :: x_middle(3)
      DOUBLE PRECISION :: R_EFF_sq, R_EFF_bs, M_EFF
      LOGICAL, PARAMETER :: ldebug = .false.
      LOGICAL:: CONVERGE_A, CONVERGE_B
      INTEGER :: shape_LL, shape_I

      integer :: start_index, end_index, do_index

! Force chain data
      INTEGER :: FCHAINC_OLD           ! old value of FCHAINC
      INTEGER :: old_size,new_size ! New force chain array size when resizing (growing) arrays
      DOUBLE PRECISION :: FC_MAG_MAX, MAG_ORIENT
      INTEGER :: FC_MAX_INDEX         ! FC index where max occurs (within LL-I pair)

! Rolling friction
      DOUBLE PRECISION :: RFTOW(3), omega_mag
!......................................................................!

! This code often causes numeric overflows, which seem not to affect the
! simulation.  TODO:  Disable FPE traps in this subroutine (?)

! Initialize cohesive forces
      IF(USE_COHESION) PostCohesive(:) = ZERO

      CALL CALC_DEM_FORCE_WITH_WALL_STL

      if (optflag1.eq.1) then
        start_index = 1
        end_index = MAX_PIP_EXIST
      else
        start_index = 1
        end_index = MAX_PIP
      endif
! Check particle LL neighbor contacts
!---------------------------------------------------------------------//

!$omp parallel default(none) private(pos,rad,cc,cc_start,cc_end,ll,i,  &
!$omp    overlap_n,vrel_t,v_rel_trans_norm,sqrt_overlap,dist,r_lm,     &
!$omp    kn_des,kt_des,phasell,phasei,etan_des,                        &
!$omp    etat_des,fn,ft,overlap_t,tangent,mag_overlap_t,               &
!$omp    eq_radius,distapart,force_coh,dist_mag,NORMAL,ftmd,fnmd,      &
!$omp    dist_cl, dist_ci, fc_tmp, tow_tmp, tow_force, qq_tmp,         &
!$omp    mew_r,rftow,omega_new,omega_mag,                              &
!$omp    euler0a,euler1a,euler2a,euler3a,                              &
!$omp    euler0b,euler1b,euler2b,euler3b,                              &
!$omp    p,axi,axia,m1,m2,n1,n2,shape_LL,shape_I,OBB_test,X0_A,X0_B,   &
!$omp    lambda_a,lambda0_a,lambda_b,lambda0_b,converge_a,converge_b,  &
!$omp    x,xcon_a,xcon_b,contacttest_a,contacttest_b,rn1n_a,rn1n_b,    &
!$omp    rn2n_a,rn2n_b,objreport,dist_li,pa,xa,Ixg,pb,xb,axib,dp,      &
!$omp    eq_ra,eq_rb,TRa,TRb,x_middle,R_EFF_sq,R_EFF_bs,M_EFF,         &
!$omp    des_etan_super,des_etat_super,hert_kn_super,hert_kt_super,    &
!$omp    qc,pft_neighbor_global,dist_cll,dist_cii)                     &
!$omp shared(max_pip,neighbors,neighbor_index,des_pos_new,des_radius,  &
!$omp    des_coll_model_enum,kn,kt,pft_neighbor,pijk,hert_kn,hert_kt,  &
!$omp    des_etan,des_etat,mew,use_cohesion, calc_cond_des, dtsolid,dt,&
!$omp    van_der_waals,vdw_outer_cutoff,vdw_inner_cutoff,              &
!$omp    hamaker_constant,asperities,surface_energy, optflag1,         &
!$omp    start_index, end_index,time,                                  &
!$omp    tow, fc, energy_eq, grav_mag, postcohesive, pmass, q_source,  &
!$omp    write_force_chain, des_col_force, iglobal_id,                 &
!$omp    old_size,new_size,                                            &
!$omp    fchain_fcmax,fc_mag_max, fc_max_index, mag_orient,            &
!$omp    FCHAINC,FCHAINC_OLD,                                          &
!$omp    FCHAIN_MIDPOINT,FCHAIN_CONTACT_POINT,                         &
!$omp    FCHAIN_ORIENT, FCHAIN_FN,FCHAIN_FT,FCHAIN_FC,                 &
!$omp    FCHAIN_LENGTH,FCHAIN_FN_MAG,FCHAIN_LOCAL_ID1,FCHAIN_LOCAL_ID2,&
!$omp    FCHAIN_GLOBAL_ID1,FCHAIN_GLOBAL_ID2, FCHAIN_OVERLAP,          &
!$omp    max_pip_exist,contact_point_a,contact_point_b,                &
!$omp    contact_lambda_a,contact_lambda_b,                            &
!$omp    super_q,super_r,super_mn,SuperDEM,SQP_contact_evolution,      &
!$omp    SQP_init_contact_evolution)

!!$omp do

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


          pos = DES_POS_NEW(LL,:)
          rad = DES_RADIUS(LL)

! feed the quaternions for LL particle
          euler0A = super_q(ll,1)
          euler1A = super_q(ll,2)
          euler2A = super_q(ll,3)
          euler3A = super_q(ll,4)

          CC_START = 1
          IF (LL.gt.1) CC_START = NEIGHBOR_INDEX(LL-1)
          CC_END = NEIGHBOR_INDEX(LL)

          DO CC = CC_START, CC_END-1
              I  = NEIGHBORS(CC)
              IF(IS_NONEXISTENT(I)) CYCLE

              R_LM = rad + DES_RADIUS(I)
              DIST(:) = DES_POS_NEW(I,:) - POS(:)
              DIST_MAG = dot_product(DIST,DIST)

              FC_TMP(:) = ZERO
! feed the quaternions for i particle
              euler0B= super_q(i,1)
              euler1B= super_q(i,2)
              euler2B= super_q(i,3)
              euler3B= super_q(i,4)

! superquadric for cohesive force and energy are not implemented at present
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
                      TOW_TMP(:,:) = ZERO

! just for post-processing mag. of cohesive forces on each particle
                      PostCohesive(LL) = PostCohesive(LL) + FORCE_COH / PMASS(LL)
                  ENDIF
              ENDIF

              IF (ENERGY_EQ) THEN
            ! Calculate conduction and radiation for thermodynamic neighbors
                  IF(CALC_COND_DES(PIJK(LL,5))) THEN
                      QQ_TMP = DES_CONDUCTION(LL, I, sqrt(DIST_MAG), PIJK(LL,5), PIJK(LL,4))

!$omp atomic
                      Q_Source(LL) = Q_Source(LL) + QQ_TMP

!$omp atomic
                      Q_Source(I) = Q_Source(I) - QQ_TMP
                  ENDIF
              ENDIF

              IF(DIST_MAG == 0) THEN
                  WRITE(*,8550) LL, I
                  ERROR STOP "division by zero"
8550              FORMAT('distance between particles is zero:',2(2x,I10))
              ENDIF

! Superquadric broad phase contact detection, orented bounding sphere and OBB
! if the distance between two superquadric particle is even larger than that of two
! bounding sphere, then cycle, else enter OBB contact detection
              IF(DIST_MAG > (R_LM*5.0d0)**2) THEN
                  PFT_NEIGHBOR(1:3,CC) =0.0
                  contact_point_A(:,CC)=0.0d0
                  contact_point_B(:,CC)=0.0d0
                  contact_lambda_A(CC)=0.0d0
                  contact_lambda_B(CC)=0.0d0
                  CYCLE
              Endif


              IF(DIST_MAG > (R_LM*1.0d0)**2) THEN
                  PFT_NEIGHBOR(1:3,CC) =0.0
                  CYCLE
              ELSE
                  p(1)=DES_POS_NEW(LL,1)
                  p(2)=DES_POS_NEW(LL,2)
                  p(3)=DES_POS_NEW(LL,3)
                  p(4)=DES_POS_NEW(I,1)
                  p(5)=DES_POS_NEW(I,2)
                  p(6)=DES_POS_NEW(I,3)
                  axi(1)=super_r(LL,1)
                  axi(2)=super_r(LL,2)
                  axi(3)=super_r(LL,3)
                  axi(4)=super_r(i,1)
                  axi(5)=super_r(i,2)
                  axi(6)=super_r(i,3)
                  m1=super_mn(LL,1)
                  n1=super_mn(LL,2)
                  m2=super_mn(I,1)
                  n2=super_mn(I,2)
! Check the shape of two particles
                  call Sq_shapes(axi(1:3),m1,n1,shape_LL)
                  call Sq_shapes(axi(4:6),m2,n2,shape_I)
! Only non-sphere vs. non-sphere, sphere vs. non-sphere
                  IF(SuperDEM .and. (shape_LL .gt. 0 .or. shape_I .gt. 0) ) then

! Superquadric broad phase contact detection---OBB
                      CALL SQ_OBB(p,axi,euler0A,euler1A,euler2A,euler3A,&
                          euler0B,euler1B,euler2B,euler3B,OBB_test)
                      IF (OBB_TEST .GT. 0.0) THEN
! Distance is positive, contact not detected!
!	                  PFT_NEIGHBOR(1:3,CC) =0.0
                          contact_point_A(:,CC)=0.0
                          contact_point_B(:,CC)=0.0
                          contact_lambda_A(CC)=0.0
                          contact_lambda_B(CC)=0.0
                          CYCLE
                      ELSE

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! SUPERQUADRIC contact detection
! deepest point method
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                          if (time<DT.and.SQP_init_contact_evolution) then
                              converge_a = .False.
                              goto 500
                          endif
! Previous contact point is used to speed up
                          X0_A(1:3)=contact_point_a(:,CC)
                          LAMBDA0_A=contact_lambda_a(CC)
                          if (ldebug)  then
                              WRITE(*,'(/," contact_point_a_old")')
                              WRITE(*,'(4(f18.8))') X0_A(1),x0_a(2),x0_a(3),LAMBDA0_A
                          endif
! First test whether a point on superquadric surface A is inside B
! Method: Newton-Raphson
                          CALL SQ_CONTACT_NEWTON_DP_A(p,axi,m1,n1,m2,n2, &
                              euler0A,euler1A,euler2A,euler3A, &
                              euler0B,euler1B,euler2B,euler3B,&
                              Rn1n_A,Rn2n_A,X0_A,LAMBDA0_A,ContactTest_A,XCON_A,LAMBDA_A,converge_a)

                          if(SQP_contact_evolution) converge_a=.false.


500                       dp = 0
                          if (converge_a) dp = dot_product(Rn1n_A,Rn2n_A)
                          if (.not. converge_a .or. dp > 1e-12) then
                              call sq_evolution_newton_method_A(p,axi,m1,n1,m2,n2,&
                                  euler0A,euler1A,euler2A,euler3A,&
                                  euler0B,euler1B,euler2B,euler3B,&
                                  Rn1n_A,Rn2n_A,XCON_A,lambda_a,contacttest_A)
                          endif
                          X(1)=XCON_A(1)
                          X(2)=XCON_A(2)
                          X(3)=XCON_A(3)
                          contact_point_A(1:3,CC)=X(1:3)
                          contact_lambda_a(CC)=LAMBDA_A
                          if (ldebug)  then
                              WRITE(*,'(/," contact_point_a_NEW")')
                              WRITE(*,'(4(f18.8))') contact_point_A(1:3,CC),LAMBDA_A
                          endif
                          IF (ContactTest_A > 0) then
! contact not detected
                              CYCLE
                          ELSE
                              if (time<DT.and.SQP_init_contact_evolution) then
                                  converge_b = .False.
                                  goto 600
                              endif
                              X0_B(:)=contact_point_b(:,CC)
                              LAMBDA0_B=contact_lambda_B(CC)

! Further find a point on B
                              CALL SQ_CONTACT_NEWTON_DP_B(p,axi,m1,n1,m2,n2, &
                                  euler0A,euler1A,euler2A,euler3A, &
                                  euler0B,euler1B,euler2B,euler3B,&
                                  Rn1n_B, Rn2n_B, X0_B,LAMBDA0_B,ContactTest_B,XCON_B,LAMBDA_B,converge_b)

                              if(SQP_contact_evolution) converge_b=.false.
!                              rn1n_a and rn1n_b seem to be negatives of each other
600                           dp = 0
                              if (converge_b) dp = dot_product(Rn1n_B,Rn2n_B)
                              if(.not. converge_b .or. dp>1e-12) then
                                  call sq_evolution_newton_method_B(p,axi,m1,n1,m2,n2,&
                                      euler0A,euler1A,euler2A,euler3A,&
                                      euler0B,euler1B,euler2B,euler3B,&
                                      Rn1n_B,Rn2n_B,XCON_B,lambda_b,contacttest_B)
                              endif

                              X(4)=XCON_B(1)
                              X(5)=XCON_B(2)
                              X(6)=XCON_B(3)
                              contact_point_B(1:3,CC)=X(4:6)
                              contact_lambda_B(CC)=LAMBDA_B
                              if (ldebug)  then
                                  WRITE(*,'(/," contact_point_B_NEW")')
                                  WRITE(*,'(4(f18.8))') contact_point_B(1:3,CC),LAMBDA_B
                              endif
! penetration between two superquadric particles
                              OBJreport = -DSQRT((XCON_B(1)-XCON_A(1))**2+&
                                  (XCON_B(2)-XCON_A(2))**2+&
                                  (XCON_B(3)-XCON_A(3))**2)

! OBJreport is negative,contact detected between two superquadrics
                              OVERLAP_N = dabs(OBJreport)

! Prevent large penetration
                              if (OVERLAP_N> 0.8*min(axi(1),axi(2),axi(3))) then
                                  OVERLAP_N = 0.8*min(axi(1),axi(2),axi(3))
                                  !WRITE(*,*) 'overlap between two superquadric particles is too large!'
                              endif
! Distance between two contact points
                              DIST(1)=X(1)-X(4)
                              DIST(2)=X(2)-X(5)
                              DIST(3)=X(3)-X(6)

! Set a limiter for the extremely small overlap
                              if (OVERLAP_N < 1.0e-10*Min(axi(1),axi(2), axi(3))) then
!                             PFT_NEIGHBOR(1:3,CC) = 0.0
                                  CYCLE
                              else
                                  NORMAL(:)=  Rn1n_A(:)
! Calculate relative velocity
                                  CALL Sq_CFRELVEL(LL,I,DES_POS_NEW(LL,:),DES_POS_NEW(I,:),X, &
                                      V_REL_TRANS_NORM, VREL_T,NORMAL,dist_li,shape_LL, shape_I)
                              endif
                          ENDIF  !(ContactTest_A > 0)
                      ENDIF  ! (OBB_TEST .GT. 0.0)
! For spherical particle
                  ELSE
                      DIST_MAG = SQRT(DIST_MAG)
                      NORMAL(:)= DIST(:)/DIST_MAG
! Calculate the normal overlap
                      OVERLAP_N = R_LM-DIST_MAG

                      IF(REPORT_EXCESS_OVERLAP) CALL PRINT_EXCESS_OVERLAP

! Calculate the components of translational relative velocity for a
! contacting particle pair and the tangent to the plane of contact
                      CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, VREL_T,NORMAL(:), DIST_MAG)
                  ENDIF ! SuperDEM

              ENDIF  ! (DIST_MAG > (R_LM - SMALL_NUMBER)**2)

              phaseLL = PIJK(LL,5)
              phaseI = PIJK(I,5)

! Hertz spring-dashpot contact model
              IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
                  sqrt_overlap = SQRT(OVERLAP_N)

                  IF (SuperDEM .and. (shape_LL .gt. 0 .or. shape_I .gt. 0)) THEN
! Equivalent radius of particle A at contact point
                      pa(1)=p(1)
                      pa(2)=p(2)
                      pa(3)=p(3)
                      xa(1)=xcon_a(1)
                      xa(2)=xcon_a(2)
                      xa(3)=xcon_a(3)
                      axia(1)=axi(1)
                      axia(2)=axi(2)
                      axia(3)=axi(3)
                      Ixg=2
                      call EQUIVALENT_RADIUS(pa,xa,axia,m1,n1,&
                          euler0A,euler1A,euler2A,euler3A, &
                          eq_ra)
                      if (ldebug) then
                          WRITE(*,'(/," equivalent radius of particle A : ", f8.6, &
                              " m ",/)') eq_ra
                      endif
! Equivalent radius of particle B at contact point
                      pb(1)=p(4)
                      pb(2)=p(5)
                      pb(3)=p(6)
                      xb(1)=xcon_b(1)
                      xb(2)=xcon_b(2)
                      xb(3)=xcon_b(3)
                      axib(1)=axi(4)
                      axib(2)=axi(5)
                      axib(3)=axi(6)
                      Ixg=2
                      call EQUIVALENT_RADIUS(pb,xb,axib,m2,n2,&
                          euler0B,euler1B,euler2B,euler3B, &
                          eq_rb)
                      if (ldebug) then
                          WRITE(*,'(/," equivalent radius of particle B : ", f8.6, &
                              " m ",/)') eq_rb
                      endif
! find the middle point between two contact points

                      call trans_matrix(euler0a, euler1a, euler2a, euler3a, TRa)
                      call trans_matrix(euler0b, euler1b, euler2b, euler3b, TRb)

                      call Find_middle_point(xa,xb,P,axi,m1,n1,m2,n2,&
                          Tra,TrB,X_MIDDLE)
                      if (ldebug) then
                          WRITE(*,*) 'contact point a and b:'
                          WRITE(*,'(6(f18.8))')   xa,xb
                          WRITE(*,*) 'middle point:'
                          WRITE(*,'(3(f18.8))')   x_middle
                      endif
! Update the value of R_EFF for superquadric particle
                      R_EFF_sq = eq_ra*eq_rb/(eq_ra + eq_rb)
                      R_EFF_bs=(DES_RADIUS(LL)*DES_RADIUS(I))/&
                          (DES_RADIUS(LL) + DES_RADIUS(I))
                      M_EFF=(PMASS(LL)*PMASS(I))/(PMASS(LL)+PMASS(I))
! Update some parameters for superquadric particle
!{SBAN
                      hert_kn_SUPER(phaseLL,phaseI)=hert_kn(phaseLL,phaseI)*SQRT(R_eff_sq)
                      hert_kt_SUPER(phaseLL,phaseI)=hert_kt(phaseLL,phaseI)*SQRT(R_eff_sq)

                      HERT_KN_SUPER(phaseI,phaseLL) = HERT_KN_SUPER(phaseLL,phaseI)
                      HERT_KT_SUPER(phaseI,phaseLL) = HERT_KT_SUPER(phaseLL,phaseI)

                      DES_ETAN_SUPER(phaseLL,phaseI)=DES_ETAN(phaseLL,phaseI)*&
                          SQRT(HERT_KN_SUPER(phaseLL,phaseI))*SQRT(M_EFF)
                      DES_ETAT_SUPER(phaseLL,phaseI)=DES_ETAT(phaseLL,phaseI)*&
                          SQRT(HERT_KT_SUPER(phaseLL,phaseI))*SQRT(M_EFF)
!}
                      DES_ETAN_SUPER(phaseI,phaseLL) = DES_ETAN_SUPER(phaseLL,phaseI)
                      DES_ETAT_SUPER(phaseI,phaseLL) = DES_ETAT_SUPER(phaseLL,phaseI)

                      KN_DES = hert_kn_super(phaseLL,phaseI)*sqrt_overlap
                      KT_DES = hert_kt_super(phaseLL,phaseI)*sqrt_overlap
                      sqrt_overlap = SQRT(sqrt_overlap)
                      ETAN_DES = DES_ETAN_super(phaseLL,phaseI)*sqrt_overlap
                      ETAT_DES = DES_ETAT_super(phaseLL,phaseI)*sqrt_overlap
                  ELSE ! for sphere
!{SBAN
                      R_EFF_bs=(DES_RADIUS(LL)*DES_RADIUS(I))/&
                          (DES_RADIUS(LL) + DES_RADIUS(I))
                      M_EFF=(PMASS(LL)*PMASS(I))/(PMASS(LL)+PMASS(I))
                      KN_DES = hert_kn(phaseLL,phaseI)*sqrt_overlap * SQRT(R_EFF_bs)
                      KT_DES = hert_kt(phaseLL,phaseI)*sqrt_overlap * SQRT(R_EFF_bs)
                      sqrt_overlap = SQRT(sqrt_overlap)
                      ETAN_DES = DES_ETAN(phaseLL,phaseI)*sqrt_overlap * sqrt(hert_kn(phaseLL,phaseI)*sqrt(R_EFF_bs)*M_EFF)
                      ETAT_DES = DES_ETAT(phaseLL,phaseI)*sqrt_overlap * sqrt(hert_kt(phaseLL,phaseI)*sqrt(R_EFF_bs)*M_EFF)
!}
                  ENDIF !SuperDEM
! Linear spring-dashpot contact model
              ELSE
                  KN_DES = KN
                  KT_DES = KT
                  ETAN_DES = DES_ETAN(phaseLL,phaseI)
                  ETAT_DES = DES_ETAT(phaseLL,phaseI)
              ENDIF

! Calculate the normal contact force
              FN(:) =  -(KN_DES * OVERLAP_N * NORMAL(:) + &
                  ETAN_DES * V_REL_TRANS_NORM * NORMAL(:))
              if (ldebug) then
                  WRITE(*,*) 'forces:'
                  WRITE(*,'(3(f18.5))')    fn(1), fn(2), fn(3)
              endif
! Calculate the tangential overlap
              if (ldebug) then
                  WRITE(*,*) 'pft_neighbor_old:'
                  WRITE(*,'(1(f18.8))')    PFT_NEIGHBOR(1:3,CC)
              endif


! Rotate the PFT from local to global
              qc(1)=euler0A
              qc(2)=euler1A
              qc(3)=euler2A
              qc(4)=euler3A
              call QROTATE(Qc,PFT_NEIGHBOR_Global(1:3),PFT_NEIGHBOR(1:3,CC),2)

              OVERLAP_T(1:3) = DTSOLID*VREL_T(1:3) + PFT_NEIGHBOR_Global(1:3)
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
              PFT_NEIGHBOR_global(1:3) = OVERLAP_T(1:3)
! rotate the pft from global to local
              call QROTATE(Qc,PFT_NEIGHBOR_Global(1:3),PFT_NEIGHBOR(1:3,CC),1)

              if (ldebug) then
                  WRITE(*,*) 'pft_neighbor_new:'
                  WRITE(*,'(1(f18.8))')  PFT_NEIGHBOR(1:3,CC)
              endif

              IF(SuperDEM .and. (shape_LL .gt. 0 .or. shape_I .gt. 0)) THEN
                  DIST_CLL(1)=X_middle(1)-DES_POS_NEW(LL,1)
                  DIST_CLL(2)=X_middle(2)-DES_POS_NEW(LL,2)
                  DIST_CLL(3)=X_middle(3)-DES_POS_NEW(LL,3)
                  DIST_CII(1)=X_middle(1)-DES_POS_NEW(I,1)
                  DIST_CII(2)=X_middle(2)-DES_POS_NEW(I,2)
                  DIST_CII(3)=X_middle(3)-DES_POS_NEW(I,3)
                  TOW_TMP(:,1)=CROSS(DIST_CLL(:),(FN(:)+FT(:)))
                  TOW_TMP(:,2)=CROSS(DIST_CII(:),-(FN(:)+FT(:)))
              ELSE
! calculate the distance from the particles' centers to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
!                 WRITE(*,*) 'contact for sphere particles!!!'
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
              FC(LL,:) = FC(LL,:) + FC_TMP(:)

!$omp atomic
              FC(I,1) = FC(I,1) - FC_TMP(1)
!$omp atomic
              FC(I,2) = FC(I,2) - FC_TMP(2)
!$omp atomic
              FC(I,3) = FC(I,3) - FC_TMP(3)

! for each particle the signs of norm and ft both flip, so add the same torque
              TOW(LL,:) = TOW(LL,:) + TOW_TMP(:,1)

!$omp atomic
              TOW(I,1)  = TOW(I,1)  + TOW_TMP(1,2)
!$omp atomic
              TOW(I,2)  = TOW(I,2)  + TOW_TMP(2,2)
!$omp atomic
              TOW(I,3)  = TOW(I,3)  + TOW_TMP(3,2)

              DES_COL_FORCE(I,:) = FC(I,:)

!######################################################################
! Force chain data - We cannot use the same approach as for spheres.
!                    With spheres, the line from one particle center to
!                    it's neighbor's center also passes through the
!                    contact point. Therefore, one segment is sufficient.
!                    With non-spherical particles, we need to use two
!                    segments: one from the particle center to the
!                    contact point, and one from the neighbor's center
!                    to the contact point.
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
               IF(SuperDEM .and. (shape_LL .gt. 0 .or. shape_I .gt. 0)) THEN
! Need to save two data points (one for each particle in contact)
! The midpoint is half-way between the particle center and the contact point
                  FCHAINC                     = FCHAINC + 1  ! Force chain counter
                  FCHAIN_MIDPOINT(FCHAINC,:)  = 0.5D0 * (DES_POS_NEW(LL,:) + X_middle(:)) ! Mid point coordinates
                  FCHAIN_CONTACT_POINT(FCHAINC,:)  = X_middle(:) !  Contact point coordinates
                  MAG_ORIENT                  = dsqrt(dot_product(DIST_CLL,DIST_CLL))
                  if(MAG_ORIENT/=ZERO) FCHAIN_ORIENT(FCHAINC,:)    = DIST_CLL(:)/MAG_ORIENT ! Orientation
                  FCHAIN_FN(FCHAINC,:)        = FN(:)    ! Normal force
                  FCHAIN_FT(FCHAINC,:)        = FT(:)    ! Tangential force
                  FCHAIN_FC(FCHAINC,:)        = FC_TMP(:)! Total force (including cohesion)
                  FCHAIN_LENGTH(FCHAINC)      = DSQRT(DIST_CLL(1)**2 + DIST_CLL(2)**2 + DIST_CLL(3)**2) !distance between center and contact point
                  FCHAIN_FN_MAG(FCHAINC)      = DSQRT(FN(1)**2 + FN(2)**2 + FN(3)**2) ! Normal force magnitude
                  FCHAIN_LOCAL_ID1(FCHAINC)   = LL
                  FCHAIN_LOCAL_ID2(FCHAINC)   = I
                  FCHAIN_GLOBAL_ID1(FCHAINC)  = iglobal_id(LL)
                  FCHAIN_GLOBAL_ID2(FCHAINC)  = iglobal_id(I)
                  FCHAIN_OVERLAP(FCHAINC)     = OVERLAP_N

                  FCHAINC                     = FCHAINC + 1  ! Force chain counter
                  FCHAIN_MIDPOINT(FCHAINC,:)  = 0.5D0 * (DES_POS_NEW(I,:) + X_middle(:)) ! Mid point coordinates
                  FCHAIN_CONTACT_POINT(FCHAINC,:)  = X_middle(:) !  Contact point coordinates
                  MAG_ORIENT                  = dsqrt(dot_product(DIST_CII,DIST_CII))
                  if(MAG_ORIENT/=ZERO) FCHAIN_ORIENT(FCHAINC,:)    = DIST_CII(:)/MAG_ORIENT ! Orientation
                  FCHAIN_FN(FCHAINC,:)        = FN(:)    ! Normal force
                  FCHAIN_FT(FCHAINC,:)        = FT(:)    ! Tangential force
                  FCHAIN_FC(FCHAINC,:)        = FC_TMP(:)! Total force (including cohesion)
                  FCHAIN_LENGTH(FCHAINC)      = DSQRT(DIST_CII(1)**2 + DIST_CII(2)**2 + DIST_CII(3)**2) !distance between center and contact point
                  FCHAIN_FN_MAG(FCHAINC)      = DSQRT(FN(1)**2 + FN(2)**2 + FN(3)**2) ! Normal force magnitude
                  FCHAIN_LOCAL_ID1(FCHAINC)   = LL
                  FCHAIN_LOCAL_ID2(FCHAINC)   = I
                  FCHAIN_GLOBAL_ID1(FCHAINC)  = iglobal_id(LL)
                  FCHAIN_GLOBAL_ID2(FCHAINC)  = iglobal_id(I)
                  FCHAIN_OVERLAP(FCHAINC)     = OVERLAP_N

               ELSE ! Contact between two spheres

                  FCHAINC                     = FCHAINC + 1  ! Force chain counter
                  FCHAIN_MIDPOINT(FCHAINC,:)  = 0.5D0 * (DES_POS_NEW(LL,:) + DES_POS_NEW(I,:)) ! Mid point coordinates
                  FCHAIN_CONTACT_POINT(FCHAINC,:)  = DES_POS_NEW(LL,:) + NORMAL(:) * (DES_RADIUS(LL) + 0.5*OVERLAP_N) ! Contact point
                  FCHAIN_ORIENT(FCHAINC,:)    = FN(:)    ! Orientation (same as Normal force)
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

               ENDIF


               IF(FCHAIN_FN_MAG(FCHAINC)>FC_MAG_MAX) THEN
                  FC_MAX_INDEX = FCHAINC  ! Keep track of max magnitude for pair LL-I
                  FC_MAG_MAX = FCHAIN_FN_MAG(FCHAINC)
               ENDIF

            ENDIF


!Rolling friction
! JFD: Assume we can use the same formulation as for spheres.
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

          ENDDO ! Neighbor CC i.e. particle I loop
          ! Update DES_COL_FORCE
          DES_COL_FORCE(LL,:) = FC(LL,:)

          IF(WRITE_FORCE_CHAIN.AND.FC_MAX_INDEX>0) FCHAIN_FCMAX(FC_MAX_INDEX) = .TRUE. ! Record where max magnitude occurs (for a given particle LL)

      ENDDO ! Particle LL loop
!$omp end do

!$omp end parallel

! just for post-processing mag. of cohesive forces on each particle
      IF(USE_COHESION .AND. VAN_DER_WAALS .AND. GRAV_MAG > ZERO) THEN
          PostCohesive(:) = PostCohesive(:)/GRAV_MAG
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
      SUBROUTINE PRINT_EXCESS_OVERLAP_SUPERDEM

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

END SUBROUTINE PRINT_EXCESS_OVERLAP_SUPERDEM

END SUBROUTINE CALC_FORCE_SUPERDEM

END MODULE SQ_CALC_FORCE_MOD
