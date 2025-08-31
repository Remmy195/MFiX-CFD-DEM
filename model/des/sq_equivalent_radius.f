#include "error.inc"

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Module name: SQ_EQUIVALENT_RADIUS_MOD                                      !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Calculate equivalent radius at the contact point                  !
!                                                                             !
! Reference:                                                                  !
!   Xi Gao, Jia Yu, Ricardo JF Portal, Jean-Francois Dietiker,                !
!   Mehrdad Shahnam and William A Rogers,Development and validation           !
!   of SuperDEM for non-spherical particulate systems using a                 !
!   superquadric particle method, Particuology,                               !
!   2021: doi.org/10.1016/j.partic.2020.1011.1007.                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
MODULE SQ_EQUIVALENT_RADIUS_MOD
      USE SQ_ROTATION_MOD
      USE SQ_PROPERTIES_MOD
CONTAINS

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: Sq_gradient                                               !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Calculate superquadric gradient vector                            !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SQ_GRADIENT(p,x,axi,m,n,&
                                       euler0,euler1,euler2,euler3, &
                                       gradient_global,gradient_local,Ixg)
      IMPLICIT NONE

      INTEGER*4 :: Ixg
      DOUBLE PRECISION :: m,n,Rn1p,euler0,euler1,euler2,euler3
      DOUBLE PRECISION :: axi(3),x(3),X_global(3),x_local(3),p(3)
      DOUBLE PRECISION :: Qc(4),Qrot(3,3),gradient_local(3),gradient_global(3)
!
      X_global(1)= X(1)-p(1)
      X_global(2)= X(2)-p(2)
      X_global(3)= X(3)-p(3)
      Qc(1)=euler0
      Qc(2)=euler1
      Qc(3)=euler2
      Qc(4)=euler3

! global position to local position
      CALL QROTATE(Qc,X_global,X_local,1)

      Rn1p = DABS(X_LOCAL(1)/axi(1))**m + DABS(X_LOCAL(2)/AXI(2))**m

! superquadric gradient,f'_x

      gradient_local(1)=n/axi(1)*DABS(X_LOCAL(1)/axi(1))**(m-1.0d0)*&
                        Rn1p**(n/m-1.0d0)*DSIGN(1.0d0,X_LOCAL(1))

! superquadric gradient,f'_y

      gradient_local(2)=n/axi(2)*DABS(X_LOCAL(2)/axi(2))**(m-1.0d0)* &
                        Rn1p**(n/m-1.0d0)*DSIGN(1.0d0,X_LOCAL(2))

! superquadric gradient,f'_z

      gradient_local(3)=n/axi(3)*DABS(X_LOCAL(3)/axi(3))**(n-1.0d0)* &
                        DSIGN(1.0d0,X_LOCAL(3))

! in the local
      IF (Ixg.EQ.1) then
         gradient_local(1)=gradient_local(1)
         gradient_local(2)=gradient_local(2)
         gradient_local(3)=gradient_local(3)
         ! cgw - the above is suspect.  but it's not used.
         ERROR STOP 1
      ELSEIF (Ixg.EQ.2) then
      CALL QROTATIONMATRIX(Qc,Qrot)
! from local to global
      CALL QMATRIXROTATEVECTOR(Qrot,gradient_global,gradient_local,Ixg)
      ENDIF
      !!write (*,*) gradient_local

      END SUBROUTINE SQ_GRADIENT

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: Sq_Hessian                                                !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Calculate sq gradient vector                                      !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SQ_Hessian(p,x,axi,m,n,&
                                        euler0,euler1,euler2,euler3, &
                                      hessian_global,hessian_local,Ixg)
      IMPLICIT NONE

      INTEGER*4 :: Ixg
      DOUBLE PRECISION :: m,n,Rn1p,euler0,euler1,euler2,euler3
      DOUBLE PRECISION :: axi(3),x(3),X_global(3),x_local(3),p(3)
      DOUBLE PRECISION :: Qc(4),Qrot(3,3)
      DOUBLE PRECISION :: xoa,yob,zoc
      DOUBLE PRECISION, PARAMETER :: tol=1.0D-10  ! Tolerance to bypass known singularity
      DOUBLE PRECISION :: hessian_local(3,3),hessian_global(3,3)
!

      X_global(1)= X(1)-p(1)
      X_global(2)= X(2)-p(2)
      X_global(3)= X(3)-p(3)
      Qc(1)=euler0
      Qc(2)=euler1
      Qc(3)=euler2
      Qc(4)=euler3
      hessian_local(:,:)=0.0d0

! global position to local position
      CALL QROTATE(Qc,X_global,X_local,1)

      xoa= DABS(X_LOCAL(1)/axi(1))
      yob= DABS(X_LOCAL(2)/axi(2))
      zoc= DABS(X_LOCAL(3)/axi(3))

! The Hessian matrix is singular at six points along the surface.
! To avoid the singularity, the matrix is evaluated slightly away
! from these points. The tolerance is hard-coded for now and may
! need to be adjusted (current cutoff is 1E-10).
      xoa = dmax1(xoa,tol)
      yob = dmax1(yob,tol)
      zoc = dmax1(zoc,tol)

      Rn1p = xoa**m + yob**m

      hessian_local(1,1)=1.0d0/axi(1)**2*n*(m-1)*xoa**(m-2)*Rn1p**(n/m-1)+&
                     1.0d0/axi(1)**2*n*(n-m)*xoa**(2*m-2)*Rn1p**(n/m-2)
! superquadric hessian,f''(xy)
      hessian_local(1,2)=1.0d0/axi(1)/axi(2)*n*(n-m)*xoa**(m-1)*&
                        yob**(m-1)*Rn1p**(n/m-2)*&
                   DSIGN(1.0d0,X_LOCAL(1))*DSIGN(1.0d0,X_LOCAL(2))
! superquadric hessian,f''(xz)
     hessian_local(1,3)= 0.0
! superquadric hessian,f''(yx)
     hessian_local(2,1)= hessian_local(1,2)
! superquadric hessian,f''(yy)
     hessian_local(2,2)= 1.0d0/axi(2)**2*n*(m-1)*yob**(m-2)*Rn1p**(n/m-1)&
                     +1.0d0/axi(2)**2*n*(n-m)*yob**(2*m-2)*Rn1p**(n/m-2)
! superquadric hessian,f''(yz)
     hessian_local(2,3)= 0.0
! superquadric hessian,f''(zx)
     hessian_local(3,1)= 0.0
! superquadric hessian,f''(zy)
     hessian_local(3,2)= 0.0
! superquadric hessian,f''(zy)
     hessian_local(3,3)= 1.0d0/axi(3)**2*n*(n-1)*zoc**(n-2)
! in the local
      IF (Ixg.EQ.1) then
        hessian_local(:,:)= hessian_local(:,:)
      ELSEIF (Ixg.EQ.2) then
      CALL QROTATIONMATRIX(Qc,Qrot)
! from local to global
      CALL QMATRIXROTATETENSOR(Qrot,hessian_global,hessian_local,Ixg)

      ENDIF
      END SUBROUTINE SQ_Hessian

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: EQUIVALENT_RADIUS                                         !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Only the rotation matrix for quaternions                          !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE EQUIVALENT_RADIUS(p,x,axi,m,n,&
                                        euler0,euler1,euler2,euler3, &
                                        eq_r)
      IMPLICIT NONE
      INTEGER :: Ixg
      DOUBLE PRECISION :: axi(3),x(3),p(3),m,n
      DOUBLE PRECISION :: euler0,euler1,euler2,euler3
      DOUBLE PRECISION :: gradient_global(3),gradient_local(3)
      DOUBLE PRECISION :: hessian_global(3,3),hessian_local(3,3)
      DOUBLE PRECISION :: gg_norm,part1,part2,kmean,eq_r,eq_r_max,svolume,eq_r_v

      Ixg=2
      eq_r =0
      gradient_global(:)=0.0
      gradient_local(:)=0.0

      call SQ_GRADIENT(p,x,axi,m,n,&
                                        euler0,euler1,euler2,euler3, &
                                      gradient_global,gradient_local,Ixg)
      hessian_global(:,:)=0.0
      hessian_local(:,:)=0.0

      call SQ_Hessian(p,x,axi,m,n,&
                                        euler0,euler1,euler2,euler3, &
                                      hessian_global,hessian_local,Ixg)

     gg_norm= DSQRT(gradient_global(1)**2 &
                  + gradient_global(2)**2 &
                  + gradient_global(3)**2)
     part1=(gradient_global(1)*hessian_global(1,1)+&
            gradient_global(2)*hessian_global(2,1)+&
            gradient_global(3)*hessian_global(3,1))*gradient_global(1) +&
           (gradient_global(1)*hessian_global(1,2)+&
            gradient_global(2)*hessian_global(2,2)+&
            gradient_global(3)*hessian_global(3,2))*gradient_global(2) +&
           (gradient_global(1)*hessian_global(1,3)+&
            gradient_global(2)*hessian_global(2,3)+&
            gradient_global(3)*hessian_global(3,3))*gradient_global(3)

     part2=hessian_global(1,1)+hessian_global(2,2)+hessian_global(3,3)
     kmean=-0.5/(gg_norm**3)*(part1-part2*gg_norm**2)
! calculate the volume of the superquadric particle
     call SQ_VOLUME(axi,m,n,svolume)
! the max curvature is the 10 times the volume equivalent sphere
     eq_r_v=1.0*(3.0*svolume/4.0/3.1415926)**(1.0/3.0)
     eq_r_max=10.0*eq_r_v

     if (kmean .eq. 0.0d0) then
       eq_r=eq_r_max
     else
       eq_r=1/dabs(kmean)
     endif
! set the minimum boundary for the particle to avoid large overlap
       if (eq_r .lt. 0.3*eq_r_v) eq_r=0.3*eq_r_v
       eq_r=min(eq_r,eq_r_max)

     END SUBROUTINE EQUIVALENT_RADIUS

END MODULE SQ_EQUIVALENT_RADIUS_MOD
