#include "error.inc"

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Module name: SQ_PROPERTIES_MOD                                             !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose:  calculate the properties of superquadric, e.g. moments of inertia!
!            volume                                                           !
!                                                                             !
! Reference:                                                                  !
!   Xi Gao, Jia Yu, Ricardo JF Portal, Jean-Francois Dietiker,                !
!   Mehrdad Shahnam and William A Rogers,Development and validation           !
!   of SuperDEM for non-spherical particulate systems using a                 !
!   superquadric particle method, Particuology,                               !
!   2021: doi.org/10.1016/j.partic.2020.1011.1007.                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
MODULE SQ_PROPERTIES_MOD
      use constant, only: PI
      use SQ_ROTATION_MOD

CONTAINS

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: TEST_SQ_PROPERTIES                                        !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: test the superquadric_properties                                  !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE TEST_SQ_PROPERTIES
     IMPLICIT NONE

     DOUBLE PRECISION :: axi(3),m,n,IXX,IYY,IZZ,SVOLUME

     axi(1)=1.0
     axi(2)=2.0
     axi(3)=3.0

     m=1.0e30  ! 2/e2=0.0
     n=1.0e30  ! 2/e1=0.0
     call sq_inertia (axi,m,n,IXX,IYY,IZZ)
     CALL SQ_VOLUME(axi,m,n,svolume)
     WRITE(*,*)'    IXX       IYY     IZZ   VOLUME of plate'
     WRITE(*,*)'  ---------------------------------'
     WRITE(*,'(4(f9.3))') IXX, IYY, IZZ, SVOLUME

     m=1.0e30   !2/e2=0.0
     n=2.0
     call sq_inertia (axi,m,n,IXX,IYY,IZZ)
     CALL SQ_VOLUME(axi,m,n,svolume)
     WRITE(*,*)'    IXX       IYY      IZZ  VOLUME of elliptical cylinder'
     WRITE(*,*)'  ---------------------------------'
     WRITE(*,'(4(f9.3))') IXX, IYY, IZZ, SVOLUME

     m=2.0
     n=2.0
     call sq_inertia (axi,m,n,IXX,IYY,IZZ)
     CALL SQ_VOLUME(axi,m,n,svolume)
     WRITE(*,*)'    IXX       IYY      IZZ    VOLUME of ellipsoid'
     WRITE(*,*)'  ---------------------------------'
     WRITE(*,'(4(f9.3))') IXX, IYY, IZZ, SVOLUME

     RETURN
     END SUBROUTINE TEST_SQ_PROPERTIES

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: SQ_min_bounding_sphere                                    !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: calculate the minimum bounding sphere diameter for                !
!  the superquadric particle                                                  !
!                                                                             !
!  Reference: Podlozhnyuk, A., Pirker, S. & Kloss, C. , "Efficient            !
!  implementation of superquadric particles in Discrete Element Method within !
!  an open-source framework. Comp. Part. Mech. 4, 101-118 (2017).             !
!  https://doi.org/10.1007/s40571-016-0131-6                                  !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE SQ_min_bounding_sphere(sqp_a, sqp_b, sqp_c, sqp_m, sqp_n, sqp_d0)
     IMPLICIT NONE
     DOUBLE PRECISION, INTENT(IN) :: sqp_a, sqp_b, sqp_c, sqp_m, sqp_n
     DOUBLE PRECISION, INTENT(OUT) :: sqp_d0
     DOUBLE PRECISION :: a, b, c, m, n, x, y, z, temp
     DOUBLE PRECISION :: alfa, beta, gama, zeta

! When n = 2 we need to reorder a and b. Since sqp_a and sqp_b have INTENT(IN)
! we need separate variables (c is used for consistent notation)
! m and n can also be overwritten to avoid a singularity
     a = sqp_a
     b = sqp_b
     c = sqp_c
     m = sqp_m
     n = sqp_n

! test case1, ellipsoidal, sqp_d0=2*a=0.002
!    a=0.001; b=0.0005; c=0.0002; n=2.0; m=2.0
! test case2, cylindrical, sqp_d0=0.0020396078 (this is approximately a cylinder)
!    a=0.001; b=0.0005; c=0.0002; n=200000; m=2.0
! test case3, cubic, sqp_d0=2*(0.0005^2+0.001^2+0.0002^2)^0.5=0.0022715634 (this is approximately a cube)
!     a=0.001; b=0.0005; c=0.0002; n=200000; m=200000
     IF(a<b) THEN
         temp = a
         a    = b
         b    = temp
     ENDIF

     IF(m .eq.2.0 .and. n .eq. 2.0) THEN ! spheroid
         sqp_d0 = 2.0*max(a,b,c)
         RETURN
     ENDIF

     IF (m .eq. 2.0) THEN  ! special case
         sqp_d0 = 2.0*sqrt(a*a + c*c)
         RETURN
     ENDIF

     ! Avoid singularity near n = 2
     IF (dabs(n - 2.0) .lt. 0.1) THEN
         IF (n.lt.2) THEN
             n = 1.9
         ELSE
             n = 2.1
         ENDIF
     ENDIF
     IF (m.eq. 2.0) THEN
         alfa = 0.0
     ELSE
         ! Avoid singularity near m = 2
         IF (dabs(m - 2.0) .lt. 0.1) THEN
             IF (m.lt.2) THEN
                 m = 1.9
             ELSE
                 m = 2.1
             ENDIF
          ENDIF
         alfa = (b/a)**(2.0/(m-2))
     ENDIF

     gama = (1.0+alfa**m)**(n/m-1.0)
     beta = (gama*c**2/a**2)**(1.0/(n-2))
     zeta = 1.0/((1.0+alfa**m)**(n/m)+beta**n)**(1.0/n)

     x      = a*zeta
     y      = b*alfa*zeta
     z      = c*beta*zeta
     sqp_d0 = 2.0*sqrt(x**2 + y**2 + z**2)

     ! WRITE(*,*)'    a,  b, c, m, n, sqp_d0'
     ! WRITE(*,*)'  ---------------------------------'
     ! WRITE(*,'(6(f25.10))') a,b,c,m,n,sqp_d0

     RETURN

     END SUBROUTINE SQ_min_bounding_sphere




!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: Sq_shapes                                                 !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: determine the shapes of superquadric particles, sphere, ellipsoid,!
!  cylinder, cuboid or a general superquadric particle                        !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE Sq_shapes(axi,m,n,shapes)
     IMPLICIT NONE
! shapes, sphere=0, cylinder=1, spheroid=2, cuboid=3, other=4
     Integer :: shapes
     DOUBLE PRECISION :: axi(3),m,n, e1, e2
     e1 = 2.0D0/n
     e2 = 2.0D0/m

     if( abs(axi(1)-axi(2)) .le. 1.0e-5 .and. abs(axi(1)-axi(3)) .le. 1.0e-5 &
         .and.  abs(e1-1.0) .le. 1.0e-5 .and. abs(e2-1.0) .le. 1.0e-5) then
         shapes=0  !sphere
     elseif(abs(e1-1.0) .le. 1.0e-5 .and. abs(e2-1.0) .le. 1.0e-5) then
         shapes=2  !ellipsoid
     elseif(abs(e1-1.0) .gt. 0.5 .and. abs(e2-1.0) .le. 1.0e-5) then
         shapes=1  ! cylinder-like
     elseif(abs(e1-1.0) .gt. 0.5 .and. abs(e2-1.0) .gt. 0.5) then
         shapes=3 ! cuboid-like
     else
         shapes=4 ! a general superquadric particle
     endif
     END SUBROUTINE Sq_shapes


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                               !
!  Subroutine name: Sq_cross_sphericity                                         !
!  Author: Xi Gao                                  Date: 25-Feb-2019            !
!                                                                               !
!  Purpose: calculate the cross_sphericity,projected area of some common shapes !
!                                                                               !
!  Reference:                                                                   !
!                                                                               !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE Sq_cross_sphericity(ugc, Qc, axi,m,n, &
             projected_area,cross_sphericity_gansor,cross_sphericity_Holzer, shapes)

     IMPLICIT NONE
! shapes, sphere=0, cylinder=1, spheroid=2, cuboid=3, other=4
     Integer :: shapes
     DOUBLE PRECISION :: L, D, LL, MM, NN

     DOUBLE PRECISION :: C_L, C_W, C_T,Ax,Ay,Az
     DOUBLE PRECISION :: Qc(4),axi(3),m,n,p
     DOUBLE PRECISION :: svolume,eq_r_v,eq_r_s
     DOUBLE PRECISION :: projected_area, cross_sphericity_gansor, cross_sphericity_holzer
     DOUBLE PRECISION :: ugcc(3),ugc(3),ugcc_norm,ugcc_local(3),ugcc_local_norm
     DOUBLE PRECISION :: ugcc_xy(3),ugcc_xy_norm
     DOUBLE PRECISION :: z_axis(3),x_axis(3),y_axis(3)
     DOUBLE PRECISION :: cos_ugcc_z,sin_ugcc_z
     DOUBLE PRECISION :: cos_ugcc_xy_x,sin_ugcc_xy_x

! unit vectors
     x_axis=(/1,0,0/)
     y_axis=(/0,1,0/)
     z_axis=(/0,0,1/)
! calculate the volume of the superquadric particle
     call SQ_VOLUME(axi,m,n,svolume)
! the max curvature is the 10 times the volume equivalent sphere
     eq_r_v=1.0*(3.0*svolume/4.0/PI)**(1.0/3.0)
! calculate normalized node center ugc
     ugcc=ugc
     ugcc_norm=(ugcc(1)**2.0+ugcc(2)**2.0+ugcc(3)**2.0)**0.5d0
     if (ugcc_norm .le. 1.0e-5) then
       cross_sphericity_gansor=1.0
       cross_sphericity_holzer=1.0
       projected_area=PI*eq_r_v**2.0d0
       goto 999
     else
        ugcc(:)=ugcc(:)/ugcc_norm
     endif

! rotate the ugc from global to local
     call QROTATE(Qc,ugcc,ugcc_local,1)

! ugc projected to the xy plane, set z component to zero
     ugcc_xy(1)=ugcc_local(1)
     ugcc_xy(2)=ugcc_local(2)
     ugcc_xy(3)=0.0

     ugcc_xy_norm=(ugcc_xy(1)**2.0+ugcc_xy(2)**2.0+ugcc_xy(3)**2.0)**0.5d0
     if (ugcc_xy_norm .le. 1.0e-5) then
      cross_sphericity_gansor=1.0
       cross_sphericity_holzer=1.0
       projected_area=PI*eq_r_v**2.0d0
       goto 999
     else
     ugcc_xy(:)=ugcc_xy(:)/ugcc_xy_norm
     endif

! angle with z axis (0,0,1), angle(1)=fai
      cos_ugcc_z=dot_product(ugcc_local,z_axis)
      sin_ugcc_z=(1.0-cos_ugcc_z**2.0)**0.5d0

! angle between x axis (1,0,0) and ugcc_xy, angle(2)=seta
      cos_ugcc_xy_x=dot_product(ugcc_xy,x_axis)
      sin_ugcc_xy_x=(1.0-cos_ugcc_xy_x**2.0)**0.5d0

! projected areas for cylinder
! L,D
     if(shapes .eq.1) then
       L=2.0*axi(3)
       D=2.0*axi(1)
       projected_area=L*D*abs(sin_ugcc_z)+PI/4.0*D**2.0*abs(cos_ugcc_z)
! projected area of a sphere=pi*r^2
      eq_r_s= (projected_area/PI)**0.5D0
      cross_sphericity_gansor=eq_r_s/eq_r_v
      cross_sphericity_holzer=(PI*eq_r_v**2.0/projected_area)

      else if (shapes .eq. 2) then
! projected areas for spheroid
      ll=sin_ugcc_z*cos_ugcc_xy_x
      mm=sin_ugcc_z*sin_ugcc_xy_x
      nn=cos_ugcc_z
      projected_area=PI*((ll*axi(2)*axi(3))**2.0+&
                                 (mm*axi(3)*axi(1))**2.0+&
                                  (nn*axi(1)*axi(2))**2.0)**0.5d0
! for ganser drag model
! projected area of a sphere=pi*r^2
      eq_r_s= (projected_area/PI)**0.5D0
! cross sphericity_gansor is the ratio between (diameter of sphere with equivalent projected area)
!   / (diameter if sphere with equivalent volume)
      cross_sphericity_Gansor=eq_r_s/eq_r_v

! for holzer drag model
! cross sphericity_holzer is the ratio between (cross sectional area of the volume equivalent sphere)
!   / (cross-sectional area of the non-spherical particle)
       cross_sphericity_holzer=(PI*eq_r_v**2.0/projected_area)

      else if (shapes .eq. 3 ) then
! projected area for cuboid,Ax=LT, Ay=LW, Az=WT
      C_L=2.0*axi(3)  ! z in the local
      C_T=2.0*axi(2)  ! y in the local
      C_W=2.0*axi(1)  ! x in the local
      Ax=C_L*C_T
      Ay=C_L*C_W
      Az=C_W*C_T
      ugcc_local_norm=(ugcc_local(1)**2.0+ugcc_local(2)**2.0+ugcc_local(3)**2.0 )**0.50
      ugcc_local=ugcc_local/ugcc_local_norm

      projected_area=(Ax+Ay+Az)/3.0
      projected_area=max(projected_area,Ax*abs(ugcc_local(1)) + &
                                    Ay*abs(ugcc_local(2)) + &
                                    Az*abs(ugcc_local(3)))

      eq_r_s= (projected_area/PI)**0.5D0
      cross_sphericity_Gansor=eq_r_s/eq_r_v
      cross_sphericity_holzer=(PI*eq_r_v**2.0/projected_area)

      else if (shapes .eq. 0 .or. shapes .eq. 4) then
       cross_sphericity_gansor=1.0
       cross_sphericity_holzer=1.0
       projected_area=PI*eq_r_v**2.0d0
      endif
  999 continue
      END SUBROUTINE Sq_cross_sphericity

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: Sq_angle_with_gravity                                     !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: calculate the superquadric particle angle with gravity            !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE Sq_angle_with_gravity(Qc, axi,m,n, cos_angle)
     USE constant, only : GRAVITY_X, GRAVITY_Y, GRAVITY_Z
     IMPLICIT NONE

     DOUBLE PRECISION :: Qc(4),axi(3),m,n,p, cos_angle,z_axis(3),cos_z
     ! Note that the shape factors m and n are unused
     DOUBLE PRECISION :: gravity(3),gravity_local(3),gravity_norm

     gravity=(/GRAVITY_X, GRAVITY_Y, GRAVITY_Z/)
     gravity_norm=(gravity(1)**2.0+gravity(2)**2.0+gravity(3)**2.0)**0.5d0
     if (gravity_norm == 0.0) then
         cos_angle = 0
     else
         gravity=gravity/gravity_norm
! particle local coordinate
         z_axis=(/0,0,1/)
! rotate the ugc from global to local
         call QROTATE(Qc,gravity,gravity_local,1)

! angle with gravity
         cos_angle=dot_product(gravity_local,z_axis)
      endif
      END SUBROUTINE Sq_angle_with_gravity

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: Sq_sphericity                                             !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: calculate the sphericity of some common shapes                    !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE Sq_sphericity(axi,m,n,sphericity,eq_r_v,shapes)
     IMPLICIT NONE
     Integer :: shapes
     DOUBLE PRECISION :: axi(3),m,n,sarea,svolume,eq_r_v,eq_r_s,sphericity

     call SQ_VOLUME(axi,m,n,svolume)
     if (shapes .eq. 1) then
      sarea=2.0*PI*axi(1)**2.0+2.0*PI*axi(1)*(2.0*axi(3))
      eq_r_v=1.0*(3.0*svolume/4.0/PI)**(1.0/3.0)
      eq_r_s=1.0*(sarea/4.0/PI)**(1.0/2.0)
      sphericity=(eq_r_v/eq_r_s)**2.0

     else if (shapes .eq. 2) then
      call SQ_area_ellipsoid(axi,sarea)
! volume of a sphere, 4/3*pi*r^3
      eq_r_v=1.0*(3.0*svolume/4.0/PI)**(1.0/3.0)
! surface area of a sphere, 4*pi*r^2
      eq_r_s=1.0*(sarea/4.0/PI)**(1.0/2.0)
      sphericity=(eq_r_v/eq_r_s)**2.0

     else if (shapes .eq. 3) then
! cuboid
      sarea=2.0*(2.0*axi(1)*2.0*axi(2)+ 2.0*axi(2)*2.0*axi(3)+ 2.0*axi(3)*2.0*axi(1))
      eq_r_v=1.0*(3.0*svolume/4.0/PI)**(1.0/3.0)
      eq_r_s=1.0*(sarea/4.0/PI)**(1.0/2.0)
      sphericity=(eq_r_v/eq_r_s)**2.0
      else if (shapes .eq. 0 .or. shapes .eq. 4) then
       eq_r_v=1.0*(3.0*svolume/4.0/PI)**(1.0/3.0)
       sphericity =1.0
      endif
      END SUBROUTINE Sq_sphericity

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: SQ_INOUT_F_local                                          !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: calculate the inside-outside function of superquadric             !
!               F(X)<0, X is inside the surface                               !
!               F(X)>0, X is outside the surface                              !
!               F(X)=0, X is inside the surface                               !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SQ_INOUT_F_local(X,axi,m,n,F)
      IMPLICIT NONE

      DOUBLE PRECISION :: x(3),axi(3),m,n,Fa
      DOUBLE PRECISION :: AA,BB,CC,Rn1p,F,EPS
      EPS=1D-12
      AA=X(1)/axi(1)
      BB=X(2)/axi(2)
      CC=X(3)/axi(3)
      IF ((DABS(AA) .LT. EPS) .and. (DABS(BB) .LT. EPS)) THEN
        f = dabs(cc)**n - 1.0
      ELSE
        Rn1p = dabs(AA)**m + dabs(BB)**m
        f= Rn1p**(n/m) + dabs(CC)**n - 1.0
      ENDIF
      RETURN
      END SUBROUTINE SQ_INOUT_F_local
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: SQ_INOUT_F                                                !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: calculate the inside-outside function of superquadric             !
!               F(X)<0, X is inside the surface                               !
!               F(X)>0, X is outside the surface                              !
!               F(X)=0, X is inside the surface                               !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SQ_INOUT_F(P,X,axi,m,n,Tra,F)
      IMPLICIT NONE

      DOUBLE PRECISION :: p(3),x(3),axi(3),m,n,Tra(3,3),Fa
      DOUBLE PRECISION :: AA,BB,CC,Rn1p,F,EPS
      intent(in):: p, x, axi, m, n, tra
      intent(out):: f
      EPS=1D-12
      AA=DABS(((X(1)-p(1))*TRa(1,1)+(X(2)-p(2))*TRa(2,1)+(X(3)-p(3))*TRa(3,1))/axi(1))
      BB=DABS(((X(1)-p(1))*TRa(1,2)+(X(2)-p(2))*TRa(2,2)+(X(3)-p(3))*TRa(3,2))/axi(2))
      CC=DABS(((X(1)-p(1))*TRa(1,3)+(X(2)-p(2))*TRa(2,3)+(X(3)-p(3))*TRa(3,3))/axi(3))

      IF ((AA .LT. EPS) .and. (BB .LT. EPS)) THEN
        f = cc**m - 1.0
      ELSE
        Rn1p = dabs(AA)**m + dabs(BB)**m
        f= Rn1p**(n/m) + dabs(CC)**n - 1.0
      ENDIF

      RETURN
      END SUBROUTINE SQ_INOUT_F

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: SQ_NORMAL                                                 !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: calculate the normals of superquadric                             !
                                                                              !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SQ_NORMAL(P,X,axi,m,n,Tra,Rn1n)
           IMPLICIT NONE
           DOUBLE PRECISION :: p(3),x(3),axi(3),m,n,Tra(3,3),Fa
           DOUBLE PRECISION :: AA,BB,CC,EPS
           DOUBLE PRECISION :: Rn1(3), Rn1n(3)
           DOUBLE PRECISION :: Rn1p, Rn1Aaux, Rn1Baux, Rn1Caux,Rn1norm
           EPS=1D-12
           AA=((X(1)-p(1))*TRa(1,1)+(X(2)-p(2))*TRa(2,1)+(X(3)-p(3))*TRa(3,1))/axi(1)
           BB=((X(1)-p(1))*TRa(1,2)+(X(2)-p(2))*TRa(2,2)+(X(3)-p(3))*TRa(3,2))/axi(2)
           CC=((X(1)-p(1))*TRa(1,3)+(X(2)-p(2))*TRa(2,3)+(X(3)-p(3))*TRa(3,3))/axi(3)
! normals calculation
           IF ((DABS(AA) .LT. EPS) .and. (DABS(BB) .LT. EPS)) THEN
               Rn1p = 0.0
               Rn1Caux= DSIGN(1.0d0,CC)*(dabs(CC)**(n-1.0))/axi(3)
               Rn1(1) = Rn1Caux*TRa(1,3)
               Rn1(2) = Rn1Caux*TRa(2,3)
               Rn1(3) = Rn1Caux*TRa(3,3)
           ELSE
               Rn1p = dabs(AA)**m + dabs(BB)**m
               Rn1Aaux = DSIGN(1.0d0,AA)*(dabs(AA)**(m-1.0))/axi(1)
               Rn1Baux = DSIGN(1.0d0,BB)*(dabs(BB)**(m-1.0))/axi(2)
               Rn1Caux = DSIGN(1.0d0,CC)*(dabs(CC)**(n-1.0))/axi(3)
               Rn1(1) = n*((Rn1p)**(n/m-1.0) * (Rn1Aaux*TRa(1,1) + &
                   Rn1Baux*TRa(1,2)) + Rn1Caux*TRa(1,3))
               Rn1(2) = n*((Rn1p)**(n/m-1.0) * (Rn1Aaux*TRa(2,1) + &
                   Rn1Baux*TRa(2,2)) + Rn1Caux*TRa(2,3))
               Rn1(3) = n*((Rn1p)**(n/m-1.0) * (Rn1Aaux*TRa(3,1) + &
                   Rn1Baux*TRa(3,2)) + Rn1Caux*TRa(3,3))
           ENDIF

           if (dabs(rn1(1)) > 1e10 .or. dabs(rn1(2)) > 1e10 .or. dabs(rn1(3)) > 1e10) then
               write(*,*) "Whoa ", rn1
     endif

     Rn1norm= DSQRT(Rn1(1)*RN1(1) + Rn1(2)*Rn1(2) + Rn1(3)*Rn1(3))

     if (rn1norm < 1e-10) then
         write(*,*) "What?", rn1norm
     endif

     IF (Rn1norm .GT. EPS) THEN
         Rn1n(1) = Rn1(1)/Rn1norm
         Rn1n(2) = Rn1(2)/Rn1norm
         Rn1n(3) = Rn1(3)/Rn1norm
     ENDIF

     RETURN
END SUBROUTINE SQ_NORMAL

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: SQ_NORMAL2                                                !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: calculate the normals of superquadric                             !
                                                                              !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SQ_NORMAL2(P,X,axi,m,n,&
                          euler0a,euler1a,euler2a,euler3a,Rn1n)
      IMPLICIT NONE

      DOUBLE PRECISION :: p(3),x(3),axi(3),m,n,Tra(3,3),Fa
      DOUBLE PRECISION :: euler0a,euler1a,euler2a,euler3a
      DOUBLE PRECISION :: AA,BB,CC,EPS
      DOUBLE PRECISION :: Rn1(3), Rn1n(3)
      DOUBLE PRECISION :: Rn1p, Rn1Aaux, Rn1Baux, Rn1Caux,Rn1norm

     EPS=1D-12

! Transformation Matrix A

     call trans_matrix(euler0a, euler1a, euler2a, euler3a, TRa)
     AA=((X(1)-p(1))*TRa(1,1)+(X(2)-p(2))*TRa(2,1)+(X(3)-p(3))*TRa(3,1))/axi(1)
     BB=((X(1)-p(1))*TRa(1,2)+(X(2)-p(2))*TRa(2,2)+(X(3)-p(3))*TRa(3,2))/axi(2)
     CC=((X(1)-p(1))*TRa(1,3)+(X(2)-p(2))*TRa(2,3)+(X(3)-p(3))*TRa(3,3))/axi(3)
! normals calculation
     IF ((DABS(AA) .LT. EPS) .and. (DABS(BB) .LT. EPS)) THEN
        Rn1p = 0.0
        Rn1Caux= DSIGN(1.0d0,CC)*(dabs(CC)**(n-1.0))/axi(3)
        Rn1(1) = Rn1Caux*TRa(1,3)
        Rn1(2) = Rn1Caux*TRa(2,3)
        Rn1(3) = Rn1Caux*TRa(3,3)
     ELSE
        Rn1p = dabs(AA)**m + dabs(BB)**m
        Rn1Aaux = DSIGN(1.0d0,AA)*(dabs(AA)**(n-1.0))/axi(1)
        Rn1Baux = DSIGN(1.0d0,BB)*(dabs(BB)**(n-1.0))/axi(2)
        Rn1Caux = DSIGN(1.0d0,CC)*(dabs(CC)**(m-1.0))/axi(3)
        Rn1(1) = n*((Rn1p)**(n/m-1.0) * (Rn1Aaux*TRa(1,1) + &
                     Rn1Baux*TRa(1,2)) + Rn1Caux*TRa(1,3))
        Rn1(2) = n*((Rn1p)**(n/m-1.0) * (Rn1Aaux*TRa(2,1) + &
                     Rn1Baux*TRa(2,2)) + Rn1Caux*TRa(2,3))
        Rn1(3) = n*((Rn1p)**(n/m-1.0) * (Rn1Aaux*TRa(3,1) + &
                     Rn1Baux*TRa(3,2)) + Rn1Caux*TRa(3,3))
     ENDIF

     Rn1norm= DSQRT(Rn1(1)**2 + Rn1(2)**2 + Rn1(3)**2)
     IF (Rn1norm .GT. EPS) THEN
	    Rn1n(1) = Rn1(1)/Rn1norm
	    Rn1n(2) = Rn1(2)/Rn1norm
	    Rn1n(3) = Rn1(3)/Rn1norm
     ENDIF
      RETURN
      END SUBROUTINE SQ_NORMAL2

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: Trans_matrix                                              !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: calculate the transformation matrix for a quaternion              !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine trans_matrix(qa0, qa1, qa2, qa3, m)
      implicit none
      double precision, intent(in):: qa0, qa1, qa2, qa3
      double precision :: q0, q1, q2, q3
      double precision, intent(out) :: m(3,3)
      double precision:: qnorm

      qnorm = dsqrt(qa0*qa0 + qa1*qa1 + qa2*qa2 + qa3*qa3)
      q0 = qa0/qnorm
      q1 = qa1/qnorm
      q2 = qa2/qnorm
      q3 = qa3/qnorm

      m(1,1) = 2*(q0*q0 + q1*q1) - 1
      m(2,2) = 2*(q0*q0 + q2*q2) - 1
      m(3,3) = 2*(q0*q0 + q3*q3) - 1
      m(2,1) = 2*(q1*q2 + q0*q3)
      m(1,2) = 2*(q1*q2 - q0*q3)
      m(1,3) = 2*(q1*q3 + q0*q2)
      m(3,1) = 2*(q1*q3 - q0*q2)
      m(3,2) = 2*(q2*q3 + q0*q1)
      m(2,3) = 2*(q2*q3 - q0*q1)

      end subroutine trans_matrix

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: SQ_INERTIA                                                !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: calculate moments of inertia of superquadric                      !
!           divided by solids density                                         !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SQ_INERTIA(axi, m, n, IXX, IYY, IZZ)
      IMPLICIT NONE

      DOUBLE PRECISION :: axi(3), m, n, IXX, IYY, IZZ, bt1, bt2, bt3, bt4
      DOUBLE PRECISION :: e1, e2


! make sure beta function returns a value for some special cases
      e1 = 2.0D0 / n
      e2 = 2.0D0 / m
!     if (e1 .eq. 0.0)  then
!         e1 = 1.0e-10
!     endif
!
!     if (e2 .eq. 0.0)  then
!         e2 = 1.0e-10
!     endif

! beta function
      CALL SQ_BETA(1.5*e2, 0.5*e2, bt1)
      CALL SQ_BETA(0.5*e1, 2.0*e1+1.0, bt2)
      CALL SQ_BETA(0.5*e2, 0.5*e2+1.0, bt3)
      CALL SQ_BETA(1.5*e1, e1+1.0, bt4)

! moments of inertia (IXX/Ro_sol) of x principal axis
      IXX=0.5 * axi(1) * axi(2) * axi(3) * e1 * e2 * (axi(2)**2 * bt1 * bt2 + &
          4.0 * axi(3)**2 * bt3 * bt4)
! moments of inertia (IYY/Ro_sol) of y principal axis
      IYY=0.5 * axi(1) * axi(2) * axi(3) * e1 * e2 * (axi(1)**2 * bt1 * bt2 + &
          4.0 * axi(3)**2 * bt3 * bt4)
! moments of inertia (IZZ/Ro_sol) of z principal axis
      IZZ=0.5 * axi(1) * axi(2) * axi(3) * e1 * e2 * (axi(1)**2 + axi(2)**2) * bt1 * bt2

      RETURN
      END SUBROUTINE SQ_INERTIA

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: SQ_VOLUME                                                 !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: calculate volume of superquadric                                  !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SQ_VOLUME(axi,m,n,svolume)
      IMPLICIT NONE

      DOUBLE PRECISION :: axi(3),m,n,svolume,bt5, bt6
      DOUBLE PRECISION :: e1, e2
      e1 = 2.0D0/n
      e2 = 2.0D0/m
      CALL SQ_BETA(0.5*e1+1,e1,bt5)
      CALL SQ_BETA(0.5*e2,0.5*e2,bt6)
! superquadric volume
      svolume=2.0*axi(1)*axi(2)*axi(3)*e1*e2*bt5*bt6
      RETURN
      END SUBROUTINE SQ_VOLUME

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: SQ_area_ellipsoid                                         !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: calculate area of superquadric ellipsoid  (e1=e2=1.0)             !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SQ_area_ellipsoid(axi,sarea)
      IMPLICIT NONE

      DOUBLE PRECISION :: axi(3),sarea,svolume,p, a, b, c
! superquadric area
      a=axi(1)
      b=axi(2)
      c=axi(3)
      p=1.6075
      sarea=4.0*PI*(((a*b)**p + (a*c)**p + (b*c)**p)/3.0)**(1.0/p)
      RETURN
      END SUBROUTINE SQ_area_ellipsoid

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: test_Sq_BETA                                              !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: This subroutine computes the beta function                        !
!           B(p,q) for p > 0 and q > 0, using subroutine Sq_BETA              !
!                                                                             !
!          Input :  p  --- Parameter  ( p > 0 )                               !
!                   q  --- Parameter  ( q > 0 )                               !
!          Output:  BT --- β(p,q)                                             !
!          Examples:                                                          !
!               p       q           β(p,q)                                    !
!              ---------------------------------                              !
!               1.5     2.0     .2666666667D+00                               !
!               2.5     2.0     .1142857143D+00                               !
!               1.5     3.0     .1523809524D+00                               !
! --------------------------------------------------------------              !
! REFERENCE: "Fortran Routines for Computation of Special                     !
!             Functions jin.ece.uiuc.edu/routines/routines.html"              !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE test_Sq_BETA
     IMPLICIT NONE

     DOUBLE PRECISION :: P, Q, BT

     WRITE(*,*)
     WRITE(*,*)'    p       q           B(p,q)'
     WRITE(*,*)'  ---------------------------------'
     P=1.5D0
     Q=2.0D0
     CALL SQ_BETA(P,Q,BT)
     WRITE(*,'(3(f9.3))') P,Q,BT
     P=2.5D0
     Q=2.0D0
     CALL SQ_BETA(P,Q,BT)
     WRITE(*,'(3(f9.3))') P,Q,BT
     P=1.5D0
     Q=3.0D0
     CALL SQ_BETA(P,Q,BT)
     WRITE(*,'(3(f9.3))') P,Q,BT

     END SUBROUTINE test_Sq_BETA

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: Sq_BETA                                                   !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Compute the beta function β(p,q)                                  !
!           Input :  p  --- Parameter  ( p > 0 )                              !
!                    q  --- Parameter  ( q > 0 )                              !
!           Output:  BT --- B(p,q)                                            !
!           Routine called: SQ_GAMMA for computing Γ(x)                       !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE SQ_BETA(P,Q,BT)
     IMPLICIT NONE

     DOUBLE PRECISION :: P, Q, BT, GP, GQ, GPQ
#ifdef __GFORTRAN__
     DOUBLE PRECISION :: lgp, lgq, lgpq
     lgp = log_gamma(p)
     lgq = log_gamma(q)
     lgpq = log_gamma(p+q)
     BT = exp(lgp + lgq - lgpq)
#else
     CALL SQ_GAMMA(P,GP)
     CALL SQ_GAMMA(Q,GQ)
     CALL SQ_GAMMA(P+Q,GPQ)
     BT=GP*GQ/GPQ
#endif
     RETURN
     END SUBROUTINE SQ_BETA

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: SQ_GAMMA                                                  !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: This subroutine computes the gamma function                       !
!       Purpose: Compute gamma function Γ(x)                                  !
!       Input :  x  --- Argument of Γ(x)                                      !
!                       ( x is not equal to 0,-1,-2,...)                      !
!       Output:  GA --- Γ(x)                                                  !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE SQ_GAMMA(X,GA)
          IMPLICIT NONE

          INTEGER :: M1, K, M
          DOUBLE PRECISION :: G(26), X, Z, R, GA, GR
          DOUBLE PRECISION, PARAMETER ::PI=3.141592653589793D0

#ifdef __GFORTRAN__
          GA = gamma(X)
          return
#else
          !!WRITE(*,*) "SQ_GAMMA", X

          IF (X.EQ.INT(X)) THEN
              IF (X.GT.0.0D0) THEN
                  GA=1.0D0
                  M1=X-1
                  DO  K=2,M1
                      GA=GA*K
                  ENDDO
              ELSE
                  GA=1.0D+300
              ENDIF
          ELSE
!              GA=GAMMA(X)

              IF (DABS(X).GT.1.0D0) THEN
                  Z=DABS(X)
                  M=INT(Z)
                  R=1.0D0
                  DO  K=1,M
                      R=R*(Z-K)
                  ENDDO
                  Z=Z-M
              ELSE
                  Z=X
              ENDIF
              DATA G/1.0D0, 0.5772156649015329D0,  &
                  -0.6558780715202538D0, -0.420026350340952D-1, &
                  0.1665386113822915D0, -.421977345555443D-1, &
                  -.96219715278770D-2, .72189432466630D-2, &
                  -.11651675918591D-2, -.2152416741149D-3, &
                  .1280502823882D-3, -.201348547807D-4, &
                  -.12504934821D-5, .11330272320D-5, &
                  -.2056338417D-6, .61160950D-8, &
                  .50020075D-8, -.11812746D-8, &
                  .1043427D-9, .77823D-11, &
                  -.36968D-11, .51D-12, &
                  -.206D-13, -.54D-14, .14D-14, .1D-15/
              GR = G(26)
              DO K = 25,1,-1
                  GR = GR*Z+G(K)
              ENDDO
              GA=1.0D0/(GR*Z)
              IF (DABS(X).GT.1.0D0) THEN
                  GA=GA*R
                  IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
              ENDIF
         ENDIF
         RETURN
#endif
     END SUBROUTINE SQ_GAMMA

END MODULE SQ_PROPERTIES_MOD
