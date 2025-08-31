#include "error.inc"

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Module name: SQ_CONTACT_EVOLUTION_MOD                            !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Contact detection between superquadric surfaces using             !
!            deeptest point method with Newton-Raphson                        !
!                                                                             !
!  Reference:                                                                 !
! Reference:                                                                  !
!   Xi Gao, Jia Yu, Ricardo JF Portal, Jean-François Dietiker,                !
!   Mehrdad Shahnam and William A Rogers,Development and validation           !
!   of SuperDEM for non-spherical particulate systems using a                 !
!   superquadric particle method, Particuology,                               !
!   2021: doi.org/10.1016/j.partic.2020.1011.1007.                            !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
MODULE SQ_CONTACT_EVOLUTION_MOD
      USE SQ_CONTACT_NEWTON_DPmethod_MOD
      USE SQ_PROPERTIES_MOD
CONTAINS

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Module name: sq_evolution_method                                 !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: volume equivalent sphere contact as the initial state, evolute to  !
!  the final state in N steps                                                 !
!                                                                             !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE sq_evolution_newton_method_A(p,axi,m1,n1,m2,n2,&
                                  euler0A,euler1A,euler2A,euler3A,&
                                  euler0B,euler1B,euler2B,euler3B,&
                                  rn1n,rn2n,X_A,lambda_a,contact_test_A)

     IMPLICIT none
     INTEGER :: steps_base,steps,k,I
     DOUBLE PRECISION :: PN
     DOUBLE PRECISION :: P(6),AXI(6),X_A0(3),X_A(3),XLB_A(3),XUB_A(3)
     DOUBLE PRECISION :: X_A00(3),LAMBDA0_A0,LAMBDA0_A
     DOUBLE PRECISION :: m1,n1,m2,n2
! equivalent volume radius
     DOUBLE PRECISION :: eq_vol_rA,eq_vol_rb,svolume_a,svolume_b
     DOUBLE PRECISION :: del_axi(6),del_e1,del_e2,del_e3,del_e4
     DOUBLE PRECISION :: axi_temp(6),e1_temp,e2_temp,e3_temp,e4_temp
     DOUBLE PRECISION :: euler0A,euler1A,euler2A,euler3A
     DOUBLE PRECISION :: euler0B,euler1B,euler2B,euler3B
     DOUBLE PRECISION :: tra(3,3),trb(3,3),Rn1n(3),Rn2n(3)

     DOUBLE PRECISION :: norm_dist,dist(3)
     DOUBLE PRECISION :: contact_test_A,lambda_a
     DOUBLE PRECISION :: xcon_a(3)
     DOUBLE PRECISION :: e1,e2,e3,e4
     logical :: converge_a
     converge_a = .false.
! number of steps to evolute from sphere to superquadric particle
     STEPS_base=10
! calculate the volume of the superquadric particle
     call SQ_VOLUME(axi(1:3),m1,n1,svolume_A)
     call SQ_VOLUME(axi(4:6),m2,n2,svolume_B)
! the max curvature is the 10 times the volume equivalent sphere (?)
     eq_vol_rA=1.0*(3.0*svolume_A/4.0/3.1415926)**(1.0/3.0)
     eq_vol_rB=1.0*(3.0*svolume_B/4.0/3.1415926)**(1.0/3.0)
! calculate increasement of each step
     del_axi(1:3)=(axi(1:3)-eq_vol_ra)/steps_base
     del_axi(4:6)=(axi(4:6)-eq_vol_rb)/steps_base
     e1 = 2.0D0 / n1
     e2 = 2.0D0 / m1
     e3 = 2.0D0 / n2
     e4 = 2.0D0 / m2
     del_e1=(e1-1)/steps_base
     del_e2=(e2-1)/steps_base
     del_e3=(e3-1)/steps_base
     del_e4=(e4-1)/steps_base

!initial guess based on sphere
     dist(1)=p(4)-p(1)
     dist(2)=p(5)-p(2)
     dist(3)=p(6)-p(3)
     norm_dist=dist(1)**2+dist(2)**2+dist(3)**2
     norm_dist=sqrt(norm_dist)
     dist(:)=dist(:)/norm_dist

     X_A00(1) = p(1) + dist(1)*eq_vol_rA
     X_A00(2) = p(2) + dist(2)*eq_vol_rA
     X_A00(3) = p(3) + dist(3)*eq_vol_rA

     Do i=1,5
     X_A0(:)=X_A00(:)
     LAMBDA0_A=0.0
        steps=steps_base+(i-1)*10
     DO k=1,steps
        axi_temp(1:3)=eq_vol_ra+ k*del_axi(1:3)*steps_base/steps
        axi_temp(4:6)=eq_vol_rb+ k*del_axi(4:6)*steps_base/steps
        e1_temp=1.0 + k*del_e1*steps_base/steps
        e2_temp=1.0 + k*del_e2*steps_base/steps
        e3_temp=1.0 + k*del_e3*steps_base/steps
        e4_temp=1.0 + k*del_e4*steps_base/steps

! method: using newton-Raphson
        converge_a= .false.
              CALL SQ_CONTACT_NEWTON_DP_A(p,axi_temp,         &
                              2.0D0/e2_temp,2.0D0/e1_temp,    &
                              2.0D0/e4_temp,2.0D0/e3_temp,    &
                              euler0A,euler1A,euler2A,euler3A,&
                              euler0B,euler1B,euler2B,euler3B,&
                              Rn1n,Rn2n,X_A0,LAMBDA0_A,        &
                              Contact_Test_A,X_A,LAMBDA_A,     &
                              converge_a)

          X_A0=X_A
          LAMBDA0_A=LAMBDA_A

     ENDDO

     IF (CONTACT_TEST_A .NE. 0) THEN
        GOTO 100
     ENDIF

     ENDDO

     IF (CONTACT_TEST_A .EQ. 0) THEN

     ENDIF

     100 CONTINUE

     END SUBROUTINE sq_evolution_newton_method_a

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Module name: sq_evolution_method                                           !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: volume equivalent sphere contact as the initial state, evolute to  !
!  the final state in N steps                                                 !
!                                                                             !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE sq_evolution_newton_method_B(p,axi,m1,n1,m2,n2,&
                                  euler0A,euler1A,euler2A,euler3A,&
                                  euler0B,euler1B,euler2B,euler3B,&
                                  Rn1n,Rn2n,X_B,lambda_b,contact_test_B)

     IMPLICIT none
     INTEGER :: steps_base,steps,k,m,n,i
     DOUBLE PRECISION :: PN
     DOUBLE PRECISION :: P(6),AXI(6),X_B0(3),X_B(3),XLB_B(3),XUB_B(3)
     DOUBLE PRECISION :: LAMBDA0_b0,LAMBDA0_b
     DOUBLE PRECISION :: m1,n1,m2,n2
! equivalent volume radius
     DOUBLE PRECISION :: eq_vol_rA,eq_vol_rb,svolume_a,svolume_b
     DOUBLE PRECISION :: del_axi(6),del_e1,del_e2,del_e3,del_e4
     DOUBLE PRECISION :: axi_temp(6),e1_temp,e2_temp,e3_temp,e4_temp
     DOUBLE PRECISION :: euler0A,euler1A,euler2A,euler3A
     DOUBLE PRECISION :: euler0B,euler1B,euler2B,euler3B
     DOUBLE PRECISION :: tra(3,3),trb(3,3),Rn1n(3),Rn2n(3)
     DOUBLE PRECISION :: e1,e2,e3,e4
     DOUBLE PRECISION :: norm_dist,dist(3)
     DOUBLE PRECISION :: contact_test_B,lambda_b
     DOUBLE PRECISION :: xcon_b(3)
     logical :: converge_b
     converge_b = .false.

! number of steps to evolute from sphere to superquadric particle
     STEPS_base=10
! calculate the volume of the superquadric particle
     call SQ_VOLUME(axi(1:3),m1,n1,svolume_A)
     call SQ_VOLUME(axi(4:6),m2,n2,svolume_B)
! the max curvature is the 10 times the volume equivalent sphere
     eq_vol_rA=1.0*(3.0*svolume_A/4.0/3.1415926)**(1.0/3.0)
     eq_vol_rB=1.0*(3.0*svolume_B/4.0/3.1415926)**(1.0/3.0)
! calculate increasement of each step
     del_axi(1:3)=(axi(1:3)-eq_vol_ra)/steps_base
     del_axi(4:6)=(axi(4:6)-eq_vol_rb)/steps_base


     e1 = 2.0d0/n1
     e2 = 2.0d0/m1
     e3 = 2.0d0/n2
     e4 = 2.0d0/m2
     del_e1=(e1-1)/steps_base
     del_e2=(e2-1)/steps_base
     del_e3=(e3-1)/steps_base
     del_e4=(e4-1)/steps_base

! initial guess based on sphere
     dist(1)=p(1)-p(4)
     dist(2)=p(2)-p(5)
     dist(3)=p(3)-p(6)
     norm_dist=dist(1)**2+dist(2)**2+dist(3)**2
     norm_dist=sqrt(norm_dist)
     dist(:)=dist(:)/norm_dist

     X_B0(1) = p(4) + dist(1)*eq_vol_rB
     X_B0(2) = p(5) + dist(2)*eq_vol_rB
     X_B0(3) = p(6) + dist(3)*eq_vol_rB

     Do i=1,5
     X_B0(:)=X_B0(:)
     lambda0_b=0.0

        steps=steps_base+(i-1)*10

     DO k=1,steps

        axi_temp(1:3)=eq_vol_rA+ k*del_axi(1:3)*steps_base/steps
        axi_temp(4:6)=eq_vol_rB+ k*del_axi(4:6)*steps_base/steps
        e1_temp=1.0d0 + k*del_e1*steps_base/steps
        e2_temp=1.0d0 + k*del_e2*steps_base/steps
        e3_temp=1.0d0 + k*del_e3*steps_base/steps
        e4_temp=1.0d0 + k*del_e4*steps_base/steps

! method: using newton-Raphson
     converge_b = .false.

                CALL SQ_CONTACT_NEWTON_DP_B(p,axi_temp,                &
                                      2.0D0/e2_temp,2.0D0/e1_temp,     &
                                      2.0D0/e4_temp,2.0D0/e3_temp,     &
                                      euler0A,euler1A,euler2A,euler3A, &
                                      euler0B,euler1B,euler2B,euler3B, &
            Rn1n, Rn2n, X_B0,lambda0_b,Contact_Test_B,X_B,lambda_b,converge_b)

          X_B0=X_B
          LAMBDA0_B=LAMBDA_B

     ENDDO

     IF (CONTACT_TEST_B .NE. 0) THEN
        GOTO 100
     ENDIF

     ENDDO


     IF (CONTACT_TEST_B .EQ. 0) THEN

     ENDIF

     100 CONTINUE

     END SUBROUTINE sq_evolution_newton_method_B

END MODULE SQ_CONTACT_EVOLUTION_MOD
