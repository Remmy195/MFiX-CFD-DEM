!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: sq_contact_wall                                        !
!  Author: Xi Gao                                 Date: 25-Feb-2019    !
!                                                                      !
!  Purpose: This module containd routines for superquadric particle    !
! contact detection with wall                                          !
!                                                                      !
! Reference:                                                           !
!   Xi Gao, Jia Yu, Ricardo JF Portal, Jean-François Dietiker,         !
!   Mehrdad Shahnam and William A Rogers,Development and validation    !
!   of SuperDEM for non-spherical particulate systems using a          !
!   superquadric particle method, Particuology,                        !
!   2021: doi.org/10.1016/j.partic.2020.1011.1007.                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE sq_contact_wall

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: sq_contact_with_stl                                     !
!  Authors: Xi Gao                               Date: 25-Feb-2019     !
!                                                                      !
!  Purpose: superquadric particle contact with wall                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE SQ_CONTACT_WITH_STL(p,axi,m,n,&
                        euler0,euler1, euler2, euler3,wall_normal_global,&
                            vertex_global,contact_point_global,&
                     contact_point_wall_global,contact_test_wall)
     USE SQ_PROPERTIES_MOD
     USE SQ_ROTATION_MOD
     IMPLICIT NONE
     INTEGER :: CONTACT_TEST_WALL
! P is the location of superquadric center, axi is the sem-axes
     DOUBLE PRECISION :: p(3), axi(3),a,b,c
! Implicit indexes of superquadric surfaces i and j (roundness)
     DOUBLE PRECISION :: m, n
! Euler parameters of superquadric surfaces, also called quaternions
     DOUBLE PRECISION :: euler0,euler1, euler2, euler3,QC(4)
! contact point on superquadric surface with maximum pennetration
     DOUBLE PRECISION :: contact_point_global(3),contact_point_local(3)
! point on wall/stl with the lowest value of shape function of superquadric particle
     DOUBLE PRECISION :: contact_point_wall_global(3)
     DOUBLE PRECISION :: contact_point_wall_local(3)
! wall normal in the global and local coordinates
     DOUBLE PRECISION :: wall_normal_global(3),wall_normal_local(3)
     DOUBLE PRECISION :: nx,ny,nz
     DOUBLE PRECISION :: vertex_global(3),vertex_local(3), nn(3)
     DOUBLE PRECISION :: dist,f,dist2
     DOUBLE PRECISION :: alpha,beta,gamma

     a=axi(1)
     b=axi(2)
     c=axi(3)

     vertex_global(:)=vertex_global(:)-p(:)
     Qc(1)=euler0
     Qc(2)=euler1
     Qc(3)=euler2
     Qc(4)=euler3

     call QROTATE(Qc,wall_normal_global,wall_normal_local,1)
     call QROTATE(Qc,vertex_global,vertex_local,1)

! nx,ny,nz point to the non fluid domain side
     nx=-wall_normal_local(1)
     ny=-wall_normal_local(2)
     nz=-wall_normal_local(3)
     nn(1)=nx
     nn(2)=ny
     nn(3)=nz

     dist=dot_product(-1.0*wall_normal_local,vertex_local)

     if (dabs(nx) <= 1e-10 .and. dabs(ny) <= 1e-10) then
        contact_point_local(1)=0.0
        contact_point_local(2)=0.0
        contact_point_local(3)=1.0
        contact_point_wall_local(1)=0.0
        contact_point_wall_local(2)=0.0
        contact_point_wall_local(3)=dist/nz
     else
        if (dabs(nx)>dabs(ny)) then
           alpha=dabs(b*ny/(a*nx))**(1.0/(m-1.0))
           gamma=(1.0+alpha**m)**(n/m-1.0)
           beta=dabs(gamma*nz*c/(nx*a))**(1.0/(n-1.0))
           contact_point_local(1)=1.0/((1.0+alpha**m)**(n/m)+&
                                  beta**n)**(1/n)

           contact_point_local(2)=alpha*contact_point_local(1)
           contact_point_local(3)=beta*contact_point_local(1)
           contact_point_wall_local(1)=dabs(dist)/(dabs(nx)*a+&
                               alpha*dabs(ny)*b+beta*dabs(nz)*c)
           contact_point_wall_local(2)=alpha*contact_point_wall_local(1)
           contact_point_wall_local(3)=beta*contact_point_wall_local(1)

        else
           alpha=dabs(a*nx/(b*ny))**(1.0/(m-1.0))
           gamma=(1.0+alpha**m)**(n/m-1.0)
           beta=dabs(gamma*nz*c/(ny*b))**(1.0/(n-1.0))
           contact_point_local(2)=1.0/((1.0+alpha**m)**(n/m)+&
                                  beta**n)**(1.0/n)


           contact_point_local(1)=alpha*contact_point_local(2)
           contact_point_local(3)=beta*contact_point_local(2)
           contact_point_wall_local(2)=dabs(dist)/(dabs(ny)*b+&
                               alpha*dabs(nx)*a+beta*dabs(nz)*c)
           contact_point_wall_local(1)=alpha*contact_point_wall_local(2)
           contact_point_wall_local(3)=beta*contact_point_wall_local(2)

         endif
         contact_point_wall_local(1)=a*contact_point_wall_local(1)*dsign(1.0d0,nx)
         contact_point_wall_local(2)=b*contact_point_wall_local(2)*dsign(1.0d0,ny)
         contact_point_wall_local(3)=c*contact_point_wall_local(3)*dsign(1.0d0,nz)

      endif

         contact_point_local(1)=a*contact_point_local(1)*dsign(1.0d0,nx)
         contact_point_local(2)=b*contact_point_local(2)*dsign(1.0d0,ny)
         contact_point_local(3)=c*contact_point_local(3)*dsign(1.0d0,nz)
         dist=dot_product(vertex_local(:)-contact_point_local(:), nn(:))
         contact_point_wall_local(:)=contact_point_local(:)+dist*nn(:)
         dist2= dot_product(vertex_local-contact_point_local(:),nn(:))

         CALL SQ_INOUT_F_local(contact_point_wall_local,axi,m,n,f)
         if (f<= 0.0d0) then
            contact_test_wall = -1 ! contact with wall
         else
            contact_test_wall =  1 ! not contact with wall
         endif
     call QROTATE(Qc,contact_point_global,contact_point_local,2)
     call QROTATE(Qc,contact_point_wall_global,contact_point_wall_local,2)

     contact_point_wall_global(:)=contact_point_wall_global(:)+p(:)
     contact_point_global(:)=contact_point_global(:)+p(:)

      END SUBROUTINE SQ_CONTACT_WITH_STL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: point_wall_projection                                   !
!  Authors: Xi Gao                               Date: 25-Feb-2019     !
!                                                                      !
!  Purpose: project a point onto the wall                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE point_wall_projection(wall_normal,point_on_wall,input_point,&
                 output_point,dFace)
     IMPLICIT NONE
     DOUBLE PRECISION :: wall_normal(3),point_on_wall(3)
     DOUBLE PRECISION :: input_point(3),output_point(3)
     DOUBLE PRECISION :: wall_point_to_input_point(3)
     DOUBLE PRECISION :: alpha,normsq,norm2sq,dFace
     wall_point_to_input_point(:) = input_point(:)-point_on_wall
     normsq=dot_product(wall_normal,wall_normal)
     norm2sq=dot_product(wall_normal,wall_point_to_input_point)
     alpha=-norm2sq/normsq
     output_point(1)=input_point(1)+alpha*wall_normal(1)
     output_point(2)=input_point(2)+alpha*wall_normal(2)
     output_point(3)=input_point(3)+alpha*wall_normal(3)
     dFace = dabs(norm2sq)/sqrt(normsq)

     Return
     END SUBROUTINE point_wall_projection

     END MODULE sq_contact_wall
