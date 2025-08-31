#include "error.inc"

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Module name: SQ_CONTACT_NEWTON_DPmethod_MOD                                !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Contact detection between superquadric surfaces using             !
!            deeptest point method with Newton-Raphson                        !
!                                                                             !
!                                                                             !
! Reference:                                                                  !
!   Xi Gao, Jia Yu, Ricardo JF Portal, Jean-Francois Dietiker,                !
!   Mehrdad Shahnam and William A Rogers, Development and validation          !
!   of SuperDEM for non-spherical particulate systems using a                 !
!   superquadric particle method, Particuology,                               !
!   2021: doi.org/10.1016/j.partic.2020.1011.1007.                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
MODULE SQ_CONTACT_NEWTON_DPmethod_MOD
     USE Quadric, only: cross_product
     USE SQ_PROPERTIES_MOD
     USE SQ_EQUIVALENT_RADIUS_MOD
     USE Sq_math_mod
     USE Gmres_general

CONTAINS

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: SQ_CONTACT_NEWTON_DP_A                                    !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Give the properties (location, semi-axes, roundness,              !
!            euler parameters) of two potential contacting superquadrics      !
!            return the contact states (penetration or distance)              !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE SQ_CONTACT_NEWTON_DP_A(p, axi, m1, n1, m2, n2, &
         euler0A, euler1A, euler2A, euler3A, &
         euler0B, euler1B, euler2B, euler3B,&
         Rn1n_A, Rn2n_A, X0_A, LAMBDA0_A, ContactTest_A, XCON_A, LAMBDAA, converge_a)
          IMPLICIT none

          INTEGER :: s, ntrial, NP, Ixg
          INTEGER :: I, K, INDX(4)
          INTEGER :: N, itr_max
          INTEGER :: m
          DOUBLE PRECISION :: tol_abs, tol_rel, threshold

          DOUBLE PRECISION :: tolx, tolf, tolj
          DOUBLE PRECISION :: D, ERRF, ERRX
          DOUBLE PRECISION :: FJACOBIA(4,4), FUNCA(4), LAMBDAA, FA, FB
          DOUBLE PRECISION :: B(4), RHS_REC(4)
          DOUBLE PRECISION :: sscale
          ! Shape exponents of superquadric surfaces i and j
          DOUBLE PRECISION :: m1, n1, m2, n2
          ! Euler parameters of superquadric surfaces, also called quaternions
          DOUBLE PRECISION :: euler0A, euler1A, euler2A, euler3A
          DOUBLE PRECISION :: euler0B, euler1B, euler2B, euler3B
          DOUBLE PRECISION :: gap, gap2, XBOUND1, XBOUND2
          ! P is the location of superquadric center, axi is the sem-axes
          DOUBLE PRECISION :: p(6), axi(6)
          DOUBLE PRECISION :: X_A(4), X0_A(3), XLB_A(4), XUB_A(4), XCON_A(3)
          DOUBLE PRECISION :: LAMBDA0_A
          ! Transformation matrix
          DOUBLE PRECISION,DIMENSION (3,3) :: TRa(3,3), TRb(3,3)

          DOUBLE PRECISION :: grad_globalA(3), grad_localA(3)
          DOUBLE PRECISION :: grad_globalB(3), grad_localB(3)
          DOUBLE PRECISION :: grad_globalA01(3), grad_localA01(3)
          DOUBLE PRECISION :: grad_globalA02(3), grad_localA02(3)
          DOUBLE PRECISION :: grad_globalA03(3), grad_localA03(3)
          DOUBLE PRECISION :: grad_globalB01(3), grad_localB01(3)
          DOUBLE PRECISION :: grad_globalB02(3), grad_localB02(3)
          DOUBLE PRECISION :: grad_globalB03(3), grad_localB03(3)

          DOUBLE PRECISION :: hess_globalA(3,3), hess_localA(3,3)
          DOUBLE PRECISION :: hess_globalB(3,3), hess_localB(3,3)
          DOUBLE PRECISION :: ContactTest_A, Rn1n_A(3), Rn2n_A(3)
          DOUBLE PRECISION :: detA, FJACOBIA_INV(4,4), X_A_DELTA(4), X_AA(4)
          DOUBLE PRECISION :: c(4,1)

          INTEGER,PARAMETER :: MAX_ITER = 100
          DOUBLE PRECISION,PARAMETER :: g = (sqrt(5.0) - 1.0) / 2.0 ! golden ratio
          DOUBLE PRECISION,PARAMETER :: tol2 = 1e-7**2  !squared tolerance for golden mean
          DOUBLE PRECISION :: AAA(4), BBB(4), X1(4), X2(4), F1, F2, DELX(4)
          DOUBLE PRECISION :: funca1(4), funca2(4), funca01(4), funca02(4), FUNCA03(4)
          DOUBLE PRECISION :: X_A01(4), X_A02(4), X_A03(4)
          DOUBLE PRECISION :: meritA, resA, tolmerit, tolres, tol, smag2
          DOUBLE PRECISION :: meritA01, resA01, meritA02, resA02, meritA03, resA03
          DOUBLE PRECISION :: meritA1, resA1, meritA2, resA2, dist(3), norm_dist
          Logical :: converge_a
          converge_a = .true.

          X_A(:) = 0.0
          ! the aim of scale is for future chemical reaction
          sscale = 1
          axi(:) = axi(:)*sscale
          p(:) = p(:)*sscale

          ! Values for Newton-Raphson
          NP = 4
          S = 4

          ! Transformation Matrix A - should we keep these quaternionic?
          ! 1st column = TRA
          TRa(1,1) = euler0A**2 + euler1A**2 - euler2A**2 - euler3A**2
          TRa(2,1) = 2*(euler1A*euler2A + euler0A*euler3A)
          TRa(3,1) = 2*(euler1A*euler3A - euler0A*euler2A)
          ! 2nd column = TRA
          TRa(1,2) = 2*(euler1A*euler2A - euler0A*euler3A)
          TRa(2,2) = euler0A**2 - euler1A**2 + euler2A**2 - euler3A**2
          TRa(3,2) = 2*(euler2A*euler3A + euler0A*euler1A)
          ! 3rd column = TRA
          TRa(1,3) = 2*(euler1A*euler3A + euler0A*euler2A)
          TRa(2,3) = 2*(euler2A*euler3A - euler0A*euler1A)
          TRa(3,3) = euler0A**2 - euler1A**2 - euler2A**2 + euler3A**2

          ! Transformation Matrix B
          ! 1st column = TRb
          TRb(1,1) = euler0B**2 + euler1B**2 - euler2B**2 - euler3B**2
          TRb(2,1) = 2*(euler1B*euler2B + euler0B*euler3B)
          TRb(3,1) = 2*(euler1B*euler3B - euler0B*euler2B)
          ! 2nd column = TRb
          TRb(1,2) = 2*(euler1B*euler2B - euler0B*euler3B)
          TRb(2,2) = euler0B**2 - euler1B**2 + euler2B**2 - euler3B**2
          TRb(3,2) = 2*(euler2B*euler3B + euler0B*euler1B)
          ! 3rd column = TRb
          TRb(1,3) = 2*(euler1B*euler3B + euler0B*euler2B)
          TRb(2,3) = 2*(euler2B*euler3B - euler0B*euler1B)
          TRb(3,3) = euler0B**2 - euler1B**2 - euler2B**2 + euler3B**2

          ! initial values
          X_A_DELTA(:) = 0.0
          X_AA(:) = 0.0
          lambdaa = 0.0
          meritA = 0.0d0
          resA = 0.0d0
          !Ixg = 2, in the global
          Ixg = 2
          ! number of tests for Newton-Raphson
          NTRIAL = 20
          ! tolerance for convergence of contact
          tolmerit = 1E-8
          tolres = 1E-8

          ! step 1
          !  test whether the previous contact point is Ok, if ok then cycle, else
          !  it will serve as the initial condition
          IF (dabs(X0_A(1)) .GT. 1.0e-10) THEN
              X_A01(1:3) = X0_A(1:3)
              X_A01(4) = LAMBDA0_A
          else
              X_A01(:) = 0.0
              meritA01 = 1.0e30
              resA01 = 1.0e30
              GOTO 500
          ENDIF

          !test whether the initial or previous point is OK

          CALL FUNC_DP_A(p, X_A01(1:3), axi, m1, n1, m2, n2, &
              euler0A, euler1A, euler2A, euler3A, &
              euler0B, euler1B, euler2B, euler3B, &
              TRA, trb, X_A01(4), FUNCA01, meritA01, resA01)


          meritA = meritA01
          resA = resA01
          X_A(:) = X_A01(:)
          funca(:) = funca01(:)

          IF(meritA01 <= tolmerit .and. resA01 <= tolres) GOTO 100

          ! step 2
          !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          !   when the previous contact point is not available, guess the contact point
          !   base on sphere
500       CONTINUE
          dist(1) = p(4) - p(1)
          dist(2) = p(5) - p(2)
          dist(3) = p(6) - p(3)
          norm_dist = sqrt(dist(1)**2 + dist(2)**2 + dist(3)**2)
          dist(:) = dist(:) / norm_dist

          X_A02(1) = p(1) + dist(1) * axi(1)
          X_A02(2) = p(2) + dist(2) * axi(2)
          X_A02(3) = p(3) + dist(3) * axi(3)

          CALL SQ_GRADIENT(p(1:3), X_A02(1:3), axi(1:3), m1, n1, &
              euler0A, euler1A, euler2A, euler3A, &
              grad_globalA02, grad_localA02, 2)
          CALL SQ_GRADIENT(p(4:6), X_A02(1:3), axi(4:6), m2, n2, &
              euler0B, euler1B, euler2B, euler3B, &
              grad_globalB02, grad_localB02, 2)

          lambdaa=dot_product(grad_globalB02, grad_globalB02) / dot_product(grad_globalA02, grad_globalA02)

          X_A02(4) = sqrt(lambdaa)

          CALL FUNC_DP_A(p, X_A02(1:3), axi, m1, n1, m2, n2, &
              euler0A, euler1A, euler2A, euler3A, &
              euler0B, euler1B, euler2B, euler3B, &
              TRA, trb, X_A02(4), FUNCA02, meritA02, resA02)

          IF(meritA02 <= tolmerit .and. resA02<=tolres) then
              meritA=meritA02
              resA=resA02
              X_A(:)=X_A02(:)
              funca(:)=funca02(:)
              GOTO 100
          endif

          ! the previous contact point may not better than the guess based on sphere
          IF (resA02 < resA01 .and. meritA02 < meritA01) then
              meritA=meritA02
              resA=resA02
              X_A(:)=X_A02(:)
              funca(:)=funca02(:)
          endif

          ! step 3
          !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          Do k=1, NTRIAL

              CALL JACOBI_DP_A(p, X_A(1:3), axi, m1, n1, m2, n2, &
                  euler0A, euler1A, euler2A, euler3A, &
                  euler0B, euler1B, euler2B, euler3B, &
                  TRA, TRB, X_A(4), FJACOBIA)

              CALL FUNC_DP_A(p, X_A(1:3), axi, m1, n1, m2, n2, &
                  euler0A, euler1A, euler2A, euler3A, &
                  euler0B, euler1B, euler2B, euler3B, &
                  TRA, trb, X_A(4), FUNCA, meritA, resA)
              ERRF = 0d0

              DO I = 1, s
                  B(I) = -FUNCA(I) ! RHS of linear equations.
              ENDDO

              call matdet4(FJACOBIA, detA)

              ! cgw 2022-06-11 why is 1 the cutoff on determinant? (see also detB below)
              if (dabs(detA) > 1.0d0)  then
                  call matinv4(FJACOBIA, FJACOBIA_inv)
                  call matvec4(FJACOBIA_inv, B, X_AA)
              else
                  ! solve linear equation using gmres_general
                  m = 1000
                  threshold = 1E-10
                  X_A_DELTA(:) = 0.0
                  call gmres(FJACOBIA, B, X_A_DELTA, 4,m, threshold, X_AA,c)
              endif

              ! first check if without Golden mean is OK
              !initial bounds

              aaa(:) = X_A(:)
              bbb(:) = X_A(:) + X_AA(:)
              X_A03(:) = X_A(:) + X_AA(:)

              CALL FUNC_DP_A(p, X_A03(1:3), axi, m1, n1, m2, n2, &
                  euler0A, euler1A, euler2A, euler3A, &
                  euler0B, euler1B, euler2B, euler3B, &
                  TRA, TRB, X_A03(4), FUNCA03, meritA03, resA03)

              IF(meritA03 <= tolmerit .and. resA03<=tolres) then
                  meritA = meritA03
                  resA = resA03
                  X_A(:) = X_A03(:)
                  GOTO 100
              ENDIF

              X_A(:) = X_A03(:)

              ! Golden mean
              smag2 = dot_product(x_aa, x_aa)
              !bracket the minimum function value:
              DO i=1, max_iter
                  x1 = aaa + G*(bbb-aaa)
                  x2 = bbb - G*(bbb-aaa)

                  CALL FUNC_DP_A(p, x1(1:3), axi, m1, n1, m2, n2, euler0A, euler1A, euler2A, euler3A, &
                      euler0B, euler1B, euler2B, euler3B, TRA, TRB, X1(4), FUNCA1, meritA1, resA1)

                  f1=resA1

                  CALL FUNC_DP_A(p, x2(1:3), axi, m1, n1, m2, n2, euler0A, euler1A, euler2A, euler3A, &
                      euler0B, euler1B, euler2B, euler3B, TRA, TRB, X2(4), FUNCA2, meritA2, resA2)

                  f2=resA2

                  IF (F1<F2) THEN
                      aaa = x2; x2 = x1; F2 = F1
                      x1 = aaa + G*(bbb-aaa)
                      CALL FUNC_DP_A(p, x1(1:3), axi, m1, n1, m2, n2, euler0A, euler1A, euler2A, euler3A, &
                          euler0B, euler1B, euler2B, euler3B, TRA, TRB, X1(4), FUNCA1, meritA1, resA1)
                      f1 = resA1
                  ELSE
                      bbb = x1; x1 = x2; f1 = f2
                      x2 = bbb - G*(bbb-aaa)
                      CALL FUNC_DP_A(p, x2(1:3), axi, m1, n1, m2, n2, euler0A, euler1A, euler2A, euler3A, &
                          euler0B, euler1B, euler2B, euler3B, TRA, TRB, X2(4), FUNCA2, meritA2, resA2)
                      f2 = resA2
                  ENDIF
                  !check convergence
                  delx = x1 - x2
                  if(dot_product(delx, delx) <= tol2*smag2) then
                      X_A = 0.5D0 * (X1 + X2)
                      goto 200
                  endif
              ENDDO  ! golden_mean

200           continue
              delx = 0.50 * (x1 + x2) - X_A

              meritA = min(meritA1, meritA2)
              resA = min(resA1, resA2)
              IF(resA <= TOLRES .and. meritA < tolmerit) then
                  GOTO 100
              endif

              DO I = 1, s
                  B(I) = -FUNCA1(I) ! RHS of linear equations.
              ENDDO

              ERRX = 0d0
              DO I = 1, s
                  ERRX = ERRX + DABS(X_AA(I))
              ENDDO

          ENDDO   !ntrial

          ! check convergence
          IF(resA > TOLRES .or. meritA > tolmerit) then
              converge_a = .false.
              GOTO 600
          endif

100       CONTINUE

          ! save the point on A
          XCON_A(1) = X_A(1)
          XCON_A(2) = X_A(2)
          XCON_A(3) = X_A(3)
          LAMBDAA = X_A(4)
          CALL SQ_INOUT_F(P(4:6), XCON_A, axi(4:6), m2, n2, Trb, Fb)

          CALL SQ_NORMAL(P(1:3), XCON_A, axi(1:3), m1, n1, Tra, Rn1n_A)
          CALL SQ_NORMAL(P(4:6), XCON_A, axi(4:6), m2, n2, Trb, Rn2n_A)

          ! IF FB > 0, Point on superquadric particle A is not in particle B
          IF (fb .GT. 0.0) THEN
              ContactTEST_A = 1
              ! Then contact detected, find the point on A, that makes OBJ_A min.
          ELSE
              ContactTEST_A = -1
          END IF

600       continue

          RETURN

     END SUBROUTINE SQ_CONTACT_NEWTON_DP_A

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: SQ_CONTACT_NEWTON_DP_B                                    !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Give the properties (location, semi-axes, roundness,              !
!            euler parameters) of two potential contacting superquadrics      !
!            return the contact states (penetration or distance)              !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE SQ_CONTACT_NEWTON_DP_B(p, axi, m1, n1, m2, n2, &
         euler0A, euler1A, euler2A, euler3A, &
         euler0B, euler1B, euler2B, euler3B, &
         Rn1n_B, Rn2n_B, X0_B, LAMBDA0_B, ContactTest_B, XCON_B, LAMBDAB, converge_b)
          IMPLICIT none

          INTEGER :: s, ntrial, NP, Ixg
          INTEGER :: I, K, INDX(4)
          INTEGER :: N, itr_max
          INTEGER :: m
          DOUBLE PRECISION :: tol_abs, tol_rel, threshold

          DOUBLE PRECISION :: tolx, tolf, tolj
          DOUBLE PRECISION :: D, ERRF, ERRX
          DOUBLE PRECISION :: FJACOBIB(4,4), FUNCB(4),LAMBDAB,FA,FB
          DOUBLE PRECISION :: B(4),rhs_rec(4)
          DOUBLE PRECISION :: sscale
          ! Shape exponents indexes of superquadric surfaces i and j
          DOUBLE PRECISION :: m1, n1, m2, n2
          ! Euler parameters of superquadric surfaces, also called quaternions
          DOUBLE PRECISION :: euler0A, euler1A, euler2A, euler3A
          DOUBLE PRECISION :: euler0B, euler1B, euler2B, euler3B
          DOUBLE PRECISION :: gap, gap2, XBOUND1, XBOUND2
          ! P is the location of superquadric center, axi is the sem-axes
          DOUBLE PRECISION :: p(6), axi(6)
          DOUBLE PRECISION :: X_B(4), X0_B(3), XLB_B(4), XUB_B(4), XCON_B(3)
          DOUBLE PRECISION :: LAMBDA0_B
          ! Transformation matrix
          DOUBLE PRECISION,DIMENSION (3,3) :: TRa(3,3),TRb(3,3)

          DOUBLE PRECISION :: grad_globalA(3), grad_localA(3)
          DOUBLE PRECISION :: grad_globalB(3), grad_localB(3)
          DOUBLE PRECISION :: grad_globalA01(3), grad_localA01(3)
          DOUBLE PRECISION :: grad_globalA02(3), grad_localA02(3)
          DOUBLE PRECISION :: grad_globalA03(3), grad_localA03(3)
          DOUBLE PRECISION :: grad_globalB01(3), grad_localB01(3)
          DOUBLE PRECISION :: grad_globalB02(3), grad_localB02(3)
          DOUBLE PRECISION :: grad_globalB03(3), grad_localB03(3)

          DOUBLE PRECISION :: hess_globalA(3,3), hess_localA(3,3)
          DOUBLE PRECISION :: hess_globalB(3,3), hess_localB(3,3)
          DOUBLE PRECISION :: ContactTest_B, Rn1n_B(3), Rn2n_B(3)
          DOUBLE PRECISION :: detB, FJACOBIB_INV(4,4), X_B_DELTA(4), X_BB(4)
          DOUBLE PRECISION :: c(4,1)

          INTEGER,PARAMETER :: MAX_ITER=100
          DOUBLE PRECISION,PARAMETER :: G = (sqrt(5.0) - 1.0) / 2.0 ! golden ratio
          DOUBLE PRECISION,PARAMETER :: TOL2 = 1E-7 ** 2 !squared tolerance for golden mean
          DOUBLE PRECISION :: AAA(4), BBB(4), X1(4), X2(4), F1, F2, DELX(4)
          DOUBLE PRECISION :: funcb1(4), funcb2(4), FUNCB01(4), FUNCB02(4), funcb03(4)
          DOUBLE PRECISION :: X_B01(4), X_B02(4), X_B03(4)
          DOUBLE PRECISION :: meritB, resB, tolmerit, tolres, smag2
          DOUBLE PRECISION :: meritB01, resB01, meritB02, resB02, meritB03, resB03
          DOUBLE PRECISION :: meritB1, resB1, meritB2, resB2, dist(3), norm_dist
          Logical :: converge_b
          converge_b = .true.

          ! the aim of scale is for future chemical reaction
          sscale = 1
          axi(:) = axi(:)*sscale
          P(:) = P(:)*sscale
          ! Values for Newton_Raphson
          NP = 4
          S = 4

          ! Transformation Matrix A
          ! 1st column = TRA
          TRa(1,1) = euler0A**2 + euler1A**2 - euler2A**2 - euler3A**2
          TRa(2,1) = 2*(euler1A*euler2A + euler0A*euler3A)
          TRa(3,1) = 2*(euler1A*euler3A - euler0A*euler2A)
          ! 2nd column = TRA
          TRa(1,2) = 2*(euler1A*euler2A - euler0A*euler3A)
          TRa(2,2) = euler0A**2 - euler1A**2 + euler2A**2 - euler3A**2
          TRa(3,2) = 2*(euler2A*euler3A + euler0A*euler1A)
          ! 3rd column = TRA
          TRa(1,3) = 2*(euler1A*euler3A + euler0A*euler2A)
          TRa(2,3) = 2*(euler2A*euler3A - euler0A*euler1A)
          TRa(3,3) = euler0A**2 - euler1A**2 - euler2A**2 + euler3A**2

          ! Transformation Matrix B
          ! 1st column = TRb
          TRb(1,1) = euler0B**2 + euler1B**2 - euler2B**2 - euler3B**2
          TRb(2,1) = 2*(euler1B*euler2B + euler0B*euler3B)
          TRb(3,1) = 2*(euler1B*euler3B - euler0B*euler2B)
          ! 2nd column = TRb
          TRb(1,2) = 2*(euler1B*euler2B - euler0B*euler3B)
          TRb(2,2) = euler0B**2 - euler1B**2 + euler2B**2 - euler3B**2
          TRb(3,2) = 2*(euler2B*euler3B + euler0B*euler1B)
          ! 3rd column = TRb
          TRb(1,3) = 2*(euler1B*euler3B + euler0B*euler2B)
          TRb(2,3) = 2*(euler2B*euler3B - euler0B*euler1B)
          TRb(3,3) = euler0B**2 - euler1B**2 - euler2B**2 + euler3B**2

          ! initial
          X_B_DELTA(:) = 0.0
          X_BB(:) = 0.0
          lambdab = 0.0
          meritB = 0.0d0
          resB = 0.0d0
          !Ixg = 2, in the global
          Ixg = 2
          NTRIAL = 20
          tolmerit = 1E-8
          tolres = 1E-8

          ! step 1
          !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          !  test whether the previous contact point is Ok, if ok then cycle, else
          !  it will serve as the initial condition
          IF (dabs(X0_B(1)) .GT. 1.0e-10) THEN
              X_B01(1:3) = X0_B(1:3)
              X_B01(4) = LAMBDA0_B
          else
              X_B01(:) = 0.0
              meritB01 = 1.0e30
              resB01 = 1.0e30
              GOTO 500
          ENDIF

          !test whether the initial or previous point is OK

          CALL FUNC_DP_B(p, X_B01(1:3), axi, m1, n1, m2, n2, &
              euler0A, euler1A, euler2A, euler3A, &
              euler0B, euler1B, euler2B, euler3B, &
              TRA, trb, X_B01(4), FUNCB01, meritB01, resB01)

          meritB = meritB01
          resB = resB01
          X_B(:) = X_B01(:)
          FUNCB(:) = FUNCB01(:)

          IF(meritB01 <= tolmerit .AND. resB01 <= tolres) GOTO 100

          ! step 2
          !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          !   when the previous contact point is not available, guess the contact point
          !   base on sphere
500       CONTINUE
          dist(1) = p(1)-p(4)
          dist(2) = p(2)-p(5)
          dist(3) = p(3)-p(6)
          norm_dist = dist(1)**2+dist(2)**2+dist(3)**2
          norm_dist = sqrt(norm_dist)
          dist(:) = dist(:)/norm_dist

          X_B02(1) = p(4) + dist(1)*axi(4)
          X_B02(2) = p(5) + dist(2)*axi(5)
          X_B02(3) = p(6) + dist(3)*axi(6)

          CALL SQ_GRADIENT(p(1:3), X_B02(1:3), axi(1:3), m1, n1, &
              euler0A, euler1A, euler2A, euler3A, &
              grad_globalA02, grad_localA02, 2)
          CALL SQ_GRADIENT(p(4:6), X_B02(1:3), axi(4:6), m2, n2, &
              euler0B, euler1B, euler2B, euler3B, &
              grad_globalB02, grad_localB02, 2)

          lambdaB=dot_product(grad_globalA02, grad_globalA02) / dot_product(grad_globalB02, grad_globalB02)

          X_B02(4) = sqrt(lambdaB)

          CALL FUNC_DP_B(p, X_B02(1:3), axi, m1, n1, m2, n2, &
              euler0A, euler1A, euler2A, euler3A, &
              euler0B, euler1B, euler2B, euler3B, &
              TRA, trb, X_B02(4), FUNCB02, meritB02, resB02)

          IF(meritB02 <= tolmerit .AND. resB02<=tolres) then
              meritb =meritb02
              resb = resb02
              X_b(:) =X_b02(:)
              FUNCB(:) = FUNCB02(:)
              GOTO 100
          endif


          IF (resB02 < resB01 .and. meritB02 < meritB01) then
              meritb = meritb02
              resb = resb02
              X_b(:) = X_b02(:)
              FUNCB(:) = FUNCB02(:)
          endif


          ! step 3
          !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          Do k=1, NTRIAL

              CALL JACOBI_DP_B(p, x_b(1:3), axi, m1, n1, m2, n2, &
                  euler0A, euler1A, euler2A, euler3A, &
                  euler0B, euler1B, euler2B, euler3B, &
                  TRA, TRB, x_b(4), FJACOBIB)
              CALL FUNC_DP_B(p, x_b(1:3), axi, m1, n1, m2, n2, euler0A, euler1A, euler2A, euler3A, &
                  euler0B, euler1B, euler2B, euler3B, TRA, TRB, x_b(4), FUNCB, meritB, resB)

              ERRF = 0d0

              DO I = 1, s
                  B(I) = -FUNCB(I) ! RHS of linear equations.
              ENDDO

              call matdet4(FJACOBIB, detB)

              if (dabs(detB) > 1.0d0)  then
                  call matinv4(FJACOBIB, FJACOBIB_inv)
                  call matvec4(FJACOBIB_inv, B, X_BB)
              else
                  ! solve linear equation using gmres_general
                  m = 1000
                  threshold = 1E-10
                  X_B_DELTA(:) = 0.0
                  call gmres(FJACOBIB, B, X_B_DELTA, 4, m, threshold, X_BB, c)
              endif

              ! first check if without Golden mean is OK
              !initial bounds
              aaa(:) = X_B(:)
              bbb(:) = X_B(:) + X_BB(:)
              X_B03(:) = X_B(:) + X_BB(:)

              CALL FUNC_DP_B(p, X_B03(1:3), axi, m1, n1, m2, n2, &
                  euler0A, euler1A, euler2A, euler3A, &
                  euler0B, euler1B, euler2B, euler3B, &
                  TRA, TRB, X_B03(4), FUNCB03, meritB03, resB03)

              IF(meritB03 <= tolmerit .and. resB03<=tolres) then
                  meritB = meritB03
                  resB = resB03
                  X_B(:) = X_B03(:)
                  GOTO 100
              ENDIF

              X_B(:) = X_B03(:)

              ! Golden mean
              smag2 = dot_product(x_bb, x_bb)
              !bracket the minimum function value:
              DO i=1, max_iter
                  x1 = aaa + G*(bbb-aaa)
                  x2 = bbb - G*(bbb-aaa)

                  CALL FUNC_DP_B(p, x1(1:3), axi, m1, n1, m2, n2, euler0A, euler1A, euler2A, euler3A, &
                      euler0B, euler1B, euler2B, euler3B, TRA, TRB, X1(4), FUNCB1, meritB1, resB1)

                  f1=resB1

                  CALL FUNC_DP_B(p, x2(1:3), axi, m1, n1, m2, n2, euler0A, euler1A, euler2A, euler3A, &
                      euler0B, euler1B, euler2B, euler3B, TRA, TRB, X2(4), FUNCB2, meritB2, resB2)

                  f2=resB2

                  IF (F1<F2) THEN
                      aaa = x2; x2 = x1; F2 = F1
                      x1 = aaa + G*(bbb-aaa)

                      CALL FUNC_DP_B(p, x1(1:3), axi, m1, n1, m2, n2, euler0A, euler1A, euler2A, euler3A, &
                          euler0B, euler1B, euler2B, euler3B, TRA, TRB, X1(4), FUNCB1, meritB1, resB1)
                      f1 = resB1
                  ELSE
                      bbb = x1; x1 = x2; f1 = f2
                      x2 = bbb - G*(bbb-aaa)
                      CALL FUNC_DP_B(p, x2(1:3), axi, m1, n1, m2, n2, euler0A, euler1A, euler2A, euler3A, &
                          euler0B, euler1B, euler2B, euler3B, TRA, TRB, X2(4), FUNCB2, meritB2, resB2)
                      f2 = resB2
                  ENDIF
                  !check convergence
                  delx = x1-x2
                  if(dot_product(delx, delx) <= tol2*smag2) then
                      X_B = 0.5D0 * (X1 + X2)
                      goto 200
                  endif
              ENDDO  ! golden_mean

200           continue
              delx = 0.50 * (x1 + x2) - X_B

              meritB = min(meritB1, meritB2)
              resB = min(resB1, resB2)
              IF(resB <= TOLRES .and. meritB < tolmerit) then
                  GOTO 100
              endif

              DO I = 1, s
                  B(I) = -FUNCB1(I) ! RHS of linear equations.
              ENDDO

              ERRX = 0d0
              DO I = 1, s
                  ERRX = ERRX + DABS(X_BB(I))
              ENDDO

          ENDDO   !ntrial

          ! check convergence
          IF(resB > TOLRES .or. meritB > tolmerit) then
              converge_b = .false.
              GOTO 600
          endif

100       CONTINUE

          ! save the point on B
          XCON_B(1) = X_B(1)
          XCON_B(2) = X_B(2)
          XCON_B(3) = X_B(3)
          LAMBDAB = X_B(4)
          CALL SQ_INOUT_F(P(1:3), XCON_B, axi(1:3), m1, n1, Tra, Fa)

          CALL SQ_NORMAL(P(1:3), XCON_B, axi(1:3), m1, n1, Tra, Rn1n_B)
          CALL SQ_NORMAL(P(4:6), XCON_B, axi(4:6), m2, n2, Trb, Rn2n_B)

          ! IF FA > 0, Point on superquadric particle B is not in particle A
          IF (fa .GT. 0.0) THEN
              ContactTEST_B = 1
              ! Then contact detected, find the point on B, that makes OBJ_B min.
          ELSE
              ContactTEST_B = -1
          END IF

600       continue

          RETURN

     END SUBROUTINE SQ_CONTACT_NEWTON_DP_B



!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: FUNC_JACOBI_DP_A                                          !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Compute values of functions and Jacobians for Superquadrics A in  !
!  the deepest point method                                                   !
!   min. FA(X)                                                                !
!   subject to FB(X)=0                                                        !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE JACOBI_DP_B(p, x, axi, m1, n1, m2, n2, &
         euler0A, euler1A, euler2A, euler3A, &
         euler0B, euler1B, euler2B, euler3B, &
         TRA, TRB, LAMBDAB, FJACOBIB)
          IMPLICIT none
          INTEGER :: IXG
          DOUBLE PRECISION :: p(6), x(3), axi(6), m1, n1, m2, n2
          DOUBLE PRECISION :: euler0A, euler1A, euler2A, euler3A
          DOUBLE PRECISION :: euler0B, euler1B, euler2B, euler3B, TRA(3, 3), TRB(3, 3)
          DOUBLE PRECISION :: grad_globalA(3), grad_globalB(3)
          DOUBLE PRECISION :: grad_localA(3), grad_localB(3)
          DOUBLE PRECISION :: hess_globalA(3, 3), hess_globalB(3, 3)
          DOUBLE PRECISION :: hess_localA(3, 3), hess_localB(3, 3)
          DOUBLE PRECISION :: LAMBDAB, FA, FB, FJACOBIB(4, 4)

          ixg=2
          CALL SQ_GRADIENT(p(1:3), X(1:3), axi(1:3), m1, n1, &
              euler0A, euler1A, euler2A, euler3A, &
              grad_globalA, grad_localA, Ixg)
          CALL SQ_GRADIENT(p(4:6), X(1:3), axi(4:6), m2, n2, &
              euler0B, euler1B, euler2B, euler3B, &
              grad_globalB, grad_localB, Ixg)

          CALL SQ_Hessian(p(1:3), X(1:3), axi(1:3), m1, n1, &
              euler0A, euler1A, euler2A, euler3A, &
              hess_globalA, hess_localA, Ixg)
          CALL SQ_Hessian(p(4:6), X(1:3), axi(4:6), m2, n2, &
              euler0B, euler1B, euler2B, euler3B, &
              hess_globalB, hess_localB, Ixg)

          ! calculate the jacobian matrix
          FJACOBIB(:, :)=0.0
          FJACOBIB(1, 1)=hess_globalA(1, 1)+lambdaB*hess_globalB(1, 1)
          FJACOBIB(1, 2)=hess_globalA(1, 2)+lambdaB*hess_globalB(1, 2)
          FJACOBIB(1, 3)=hess_globalA(1, 3)+lambdaB*hess_globalB(1, 3)
          FJACOBIB(1, 4)=grad_globalB(1)

          FJACOBIB(2, 1)=hess_globalA(2, 1)+lambdaB*hess_globalB(2, 1)
          FJACOBIB(2, 2)=hess_globalA(2, 2)+lambdaB*hess_globalB(2, 2)
          FJACOBIB(2, 3)=hess_globalA(2, 3)+lambdaB*hess_globalB(2, 3)
          FJACOBIB(2, 4)=grad_globalB(2)

          FJACOBIB(3, 1)=hess_globalA(3, 1)+lambdaB*hess_globalB(3, 1)
          FJACOBIB(3, 2)=hess_globalA(3, 2)+lambdaB*hess_globalB(3, 2)
          FJACOBIB(3, 3)=hess_globalA(3, 3)+lambdaB*hess_globalB(3, 3)
          FJACOBIB(3, 4)=grad_globalB(3)

          FJACOBIB(4, 1)=grad_globalB(1)
          FJACOBIB(4, 2)=grad_globalB(2)
          FJACOBIB(4, 3)=grad_globalB(3)
          FJACOBIB(4, 4)=0.0

     END SUBROUTINE JACOBI_DP_B

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: JACOBI_DP_A                                               !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Compute values of functions and Jacobians for Superquadrics A in  !
!  the deepest point method                                                   !
!   min. FB(X)                                                                !
!   subject to FA(X)=0                                                        !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE JACOBI_DP_A(p, x, axi, m1, n1, m2, n2, &
         euler0A, euler1A, euler2A, euler3A, &
         euler0B, euler1B, euler2B, euler3B, &
         TRA, TRB, LAMBDAA, FJACOBIA)

          IMPLICIT none
          INTEGER :: IXG
          DOUBLE PRECISION :: p(6), x(3), axi(6), m1, n1, m2, n2
          DOUBLE PRECISION :: euler0A, euler1A, euler2A, euler3A
          DOUBLE PRECISION :: euler0B, euler1B, euler2B, euler3B
          DOUBLE PRECISION :: TRA(3, 3), TRB(3, 3)
          DOUBLE PRECISION :: grad_globalA(3), grad_globalB(3)
          DOUBLE PRECISION :: hess_globalA(3, 3), hess_globalB(3, 3)
          DOUBLE PRECISION :: hess_localA(3, 3), hess_localB(3, 3)
          DOUBLE PRECISION :: grad_localA(3), grad_localB(3)
          DOUBLE PRECISION :: LAMBDAA, FJACOBIA(4, 4)
          logical, parameter :: ldebug = .false.
          ixg=2

          CALL SQ_GRADIENT(p(1:3), X(1:3), axi(1:3), m1, n1, &
              euler0A, euler1A, euler2A, euler3A, &
              grad_globalA, grad_localA, Ixg)
          CALL SQ_GRADIENT(p(4:6), X(1:3), axi(4:6), m2, n2, &
              euler0B, euler1B, euler2B, euler3B, &
              grad_globalB, grad_localB, Ixg)

          CALL SQ_Hessian(p(1:3), X(1:3), axi(1:3), m1, n1, &
              euler0A, euler1A, euler2A, euler3A, &
              hess_globalA, hess_localA, Ixg)
          CALL SQ_Hessian(p(4:6), X(1:3), axi(4:6), m2, n2, &
              euler0B, euler1B, euler2B, euler3B, &
              hess_globalB, hess_localB, Ixg)
          if (ldebug) then
              WRITE(*, '(/, " hessian of particle A")')
              WRITE(*, '(3(f12.2))') hess_globalA(1, 1), hess_globalA(1, 2), hess_globalA(1, 3)
              WRITE(*, '(3(f12.2))') hess_globalA(2, 1), hess_globalA(2, 2), hess_globalA(2, 3)
              WRITE(*, '(3(f12.2))') hess_globalA(3, 1), hess_globalA(3, 2), hess_globalA(3, 3)
              WRITE(*, '(/, " hessian of particle b")')
              WRITE(*, '(3(f12.2))') hess_globalb(1, 1), hess_globalb(1, 2), hess_globalb(1, 3)
              WRITE(*, '(3(f12.2))') hess_globalb(2, 1), hess_globalb(2, 2), hess_globalb(2, 3)
              WRITE(*, '(3(f12.2))') hess_globalb(3, 1), hess_globalb(3, 2), hess_globalb(3, 3)
              WRITE(*, '(/, " gradient of particle A")')
              WRITE(*, '(3(f12.2))') grad_globalA(1), grad_globalA(2), grad_globalA(3)
              WRITE(*, '(/, " gradient of particle b")')
              WRITE(*, '(3(f12.2))') grad_globalB(1), grad_globalB(2), grad_globalB(3)
          endif
          ! calculate the jacobian matrix
          FJACOBIA(:,:) = 0.0
          FJACOBIA(1,1) = hess_globalB(1,1) + lambdaA*hess_globalA(1,1)
          FJACOBIA(1,2) = hess_globalB(1,2) + lambdaA*hess_globalA(1,2)
          FJACOBIA(1,3) = hess_globalB(1,3) + lambdaA*hess_globalA(1,3)
          FJACOBIA(1,4) = grad_globalA(1)

          FJACOBIA(2,1) = hess_globalB(2,1) + lambdaA*hess_globalA(2,1)
          FJACOBIA(2,2) = hess_globalB(2,2) + lambdaA*hess_globalA(2,2)
          FJACOBIA(2,3) = hess_globalB(2,3) + lambdaA*hess_globalA(2,3)
          FJACOBIA(2,4) = grad_globalA(2)

          FJACOBIA(3,1) = hess_globalB(3,1) + lambdaA*hess_globalA(3,1)
          FJACOBIA(3,2) = hess_globalB(3,2) + lambdaA*hess_globalA(3,2)
          FJACOBIA(3,3) = hess_globalB(3,3) + lambdaA*hess_globalA(3,3)
          FJACOBIA(3,4) = grad_globalA(3)

          FJACOBIA(4,1) = grad_globalA(1)
          FJACOBIA(4,2) = grad_globalA(2)
          FJACOBIA(4,3) = grad_globalA(3)
          FJACOBIA(4,4) = 0.0

     END SUBROUTINE JACOBI_DP_A

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: FUNC_DP_B                                                 !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Comopute values of functions and Jacobians for Superquadrics A in !
!  the deepest point method                                                   !
!   min. FA(X)                                                                !
!   subject to FB(X)=0                                                        !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE FUNC_DP_B(p, x, axi, m1, n1, m2, n2, euler0A, euler1A, euler2A, euler3A, &
         euler0B, euler1B, euler2B, euler3B, TRA, TRB, LAMBDAB, FUNCB, &
         meritB, resB)

          IMPLICIT none
          INTEGER :: IXG
          DOUBLE PRECISION :: p(6), x(3), axi(6), m1, n1, m2, n2
          DOUBLE PRECISION :: euler0A, euler1A, euler2A, euler3A
          DOUBLE PRECISION :: euler0B, euler1B, euler2B, euler3B
          DOUBLE PRECISION :: TRA(3, 3), TRB(3, 3)
          DOUBLE PRECISION :: grad_globalA(3), grad_globalB(3)
          DOUBLE PRECISION :: grad_localA(3), grad_localB(3)
          DOUBLE PRECISION :: LAMBDAB, FA, FB, FUNCB(4)
          DOUBLE PRECISION :: meritb, resb, v(3), VV, V1, V2, MM1, MM2
          logical, parameter :: ldebug = .false.

          ixg=2
          CALL SQ_GRADIENT(p(1:3), X(1:3), axi(1:3), m1, n1, &
              euler0A, euler1A, euler2A, euler3A, &
              grad_globalA, grad_localA, Ixg)
          CALL SQ_GRADIENT(p(4:6), X(1:3), axi(4:6), m2, n2, &
              euler0B, euler1B, euler2B, euler3B, &
              grad_globalB, grad_localB, Ixg)
          CALL SQ_INOUT_F(P(1:3), X(1:3), axi(1:3), m1, n1, Tra, Fa)
          CALL SQ_INOUT_F(P(4:6), X(1:3), axi(4:6), m2, n2, TrB, FB)

          ! RHS functions
          FUNCB(1)=  grad_globalA(1) + lambdaB*grad_globalB(1)
          FUNCB(2)=  grad_globalA(2) + lambdaB*grad_globalB(2)
          FUNCB(3)=  grad_globalA(3) + lambdaB*grad_globalB(3)
          FUNCB(4)=  fb

          CALL crossvec3(grad_globalA, grad_globalB, V)

          VV=dot_product(V, V)
          v1=dot_product(grad_globalA, grad_globalA)
          v2=dot_product(grad_globalB, grad_globalB)

          MM1=DABS(VV)/MAX(DABS(V1), 1e-16)/MAX(DABS(V2), 1E-16)
          MM2=DABS(Fb)/(MAX(dabs(FA), dabs(FB))+1E-16)
          meritb=max(MM1, MM2**2)
          resb=dot_product(FUncb, FUNCb)
          if (ldebug) then
              WRITE(*, '(/, " value of FB and fa")')
              WRITE(*, '(1(f16.8))') fb, fa
              WRITE(*, '(/, " value of M1, M2")')
              WRITE(*, '(1(f16.8))') MM1, MM2
          endif

     END SUBROUTINE FUNC_DP_B

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: FUNC_DP_A                                                 !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Compute values of functions and Jacobians for Superquadrics A in  !
!  the deepest point method                                                   !
!   min. FB(X)                                                                !
!   subject to FA(X)=0                                                        !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE FUNC_DP_A(p, x, axi, m1, n1, m2, n2, euler0A, euler1A, euler2A, euler3A, &
         euler0B, euler1B, euler2B, euler3B, TRA, TRB, LAMBDAA, FUNCA, &
         meritA, resA)

          IMPLICIT none
          INTEGER :: IXG
          DOUBLE PRECISION :: p(6), x(3), axi(6), m1, n1, m2, n2
          DOUBLE PRECISION :: euler0A, euler1A, euler2A, euler3A
          DOUBLE PRECISION :: euler0B, euler1B, euler2B, euler3B
          DOUBLE PRECISION :: grad_globalA(3), grad_globalB(3)
          DOUBLE PRECISION :: TRA(3, 3), TRB(3, 3)
          DOUBLE PRECISION :: grad_localA(3), grad_localB(3)
          DOUBLE PRECISION :: LAMBDAA, FA, FB, FUNCA(4)
          DOUBLE PRECISION :: meritA, resA, v(3), VV, V1, V2, MM1, MM2
          ixg=2
          CALL SQ_GRADIENT(p(1:3), X(1:3), axi(1:3), m1, n1, &
              euler0A, euler1A, euler2A, euler3A, &
              grad_globalA, grad_localA, Ixg)
          CALL SQ_GRADIENT(p(4:6), X(1:3), axi(4:6), m2, n2, &
              euler0B, euler1B, euler2B, euler3B, &
              grad_globalB, grad_localB, Ixg)

          CALL SQ_INOUT_F(P(1:3), X(1:3), axi(1:3), m1, n1, Tra, Fa)
          CALL SQ_INOUT_F(P(4:6), X(1:3), axi(4:6), m2, n2, Trb, Fb)
          ! RHS functions
          FUNCA(1) =  grad_globalB(1) + lambdaA*grad_globalA(1)
          FUNCA(2) =  grad_globalB(2) + lambdaA*grad_globalA(2)
          FUNCA(3) =  grad_globalB(3) + lambdaA*grad_globalA(3)
          FUNCA(4) =  fa
          !     fb=fb

          CALL crossvec3(grad_globalA, grad_globalB, V)
          VV = dot_product(V, V)
          v1 = dot_product(grad_globalA, grad_globalA)
          v2 = dot_product(grad_globalB, grad_globalB)
          MM1 = DABS(VV)/MAX(DABS(V1), 1e-16)/MAX(DABS(V2), 1E-16)
          MM2 = DABS(Fa)/(MAX(dabs(FA), dabs(FB))+1E-16)
          meritA = max(MM1, MM2**2)
          resA = dot_product(FUNCA, FUNCA)

     END SUBROUTINE FUNC_DP_A


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: Find_middle_point                                         !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: find the middle point between two contact points, that make       !
!            min(FA+FB)                                                       !
!                                                                             !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE Find_middle_point(XCON_A, XCON_B, P, axi, m1, n1, m2, n2, Tra, TrB, &
         X_MIDDLE)
          IMPLICIT NONE
          INTEGER :: max_iter, i
          DOUBLE PRECISION, intent(in) :: XCON_A(3), XCON_B(3), P(6), axi(6)
          DOUBLE PRECISION, intent(in) :: m1, n1, m2, n2, Tra(3, 3), Trb(3, 3)
          DOUBLE PRECISION :: a(3), b(3), s(3), smag2, x1(3), x2(3)
          DOUBLE PRECISION :: F1A, F2A, F1B, F2B, F1, F2, DELX(3), resid2, X_middle(3)
          DOUBLE PRECISION,PARAMETER :: G = (sqrt(5.0) - 1.0) / 2.0 ! golden ratio
          DOUBLE PRECISION,PARAMETER :: TOL2 = 1.0e-5 ** 2
          max_iter = 100

          a(:) = XCON_A(:)
          b(:) = XCON_B(:)
          s(:) = b(:) - a(:)
          smag2 = dot_product(s, s)

          DO i=1, max_iter
              x1 = a + G*(b-a)
              x2 = b - G*(b-a)
              CALL SQ_INOUT_F(P(1:3), X1, axi(1:3), m1, n1, Tra, F1A)
              CALL SQ_INOUT_F(P(4:6), X1, axi(4:6), m2, n2, Trb, F1B)
              F1=F1A + F1B
              CALL SQ_INOUT_F(P(1:3), X2, axi(1:3), m1, n1, Tra, F2A)
              CALL SQ_INOUT_F(P(4:6), X2, axi(4:6), m2, n2, Trb, F2B)
              F2=F2A + F2B

              IF (F1 < F2) THEN
                  a=x2; x2=x1; F2=F1
                  x1=a + G*(b-a)
                  CALL SQ_INOUT_F(P(1:3), X1, axi(1:3), m1, n1, Tra, F1A)
                  CALL SQ_INOUT_F(P(4:6), X1, axi(4:6), m2, n2, Trb, F1B)
                  F1 = F1A + F1B
              ELSE
                  b=x1; x1=x2; f1=f2
                  x2 = b - G*(b-a)
                  CALL SQ_INOUT_F(P(1:3), X2, axi(1:3), m1, n1, Tra, F2A)
                  CALL SQ_INOUT_F(P(4:6), X2, axi(4:6), m2, n2, Trb, F2B)
                  F2 = F2A + F2B
              ENDIF
              ! check convergence
              delx = x1 - x2
              resid2 = dot_product(delx, delx)
              if (resid2 <= TOL2*smag2) then
                  go to 100
              endif

          ENDDO
100       continue
          X_middle= (x1 + x2) / 2.0
          IF (resid2 > Tol2*smag2) then
              X_middle = (XCON_A + XCON_B) / 2.0
          ENDIF

          RETURN
     END SUBROUTINE Find_middle_point

END MODULE SQ_CONTACT_NEWTON_DPmethod_MOD
