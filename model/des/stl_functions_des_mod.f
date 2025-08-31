!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: stl_functions_des                                      !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: This module contains routines for geometric interaction    !
!  required for STL files.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE stl_functions_des

      IMPLICIT NONE

! Use this module only to define functions and subroutines.
      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ClosestPtPointTriangle                                  !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine ClosestPtPointTriangle(NF, pointp, points, &
                               closest_point, point_or_edge_int, edge)
      USE param1, only: zero, one
      USE discretelement, only: dimn
      USE geometry, only: NO_K
      USE stl, only: ignore_stl_edge_intersection
      USE param1, only: undefined


      IMPLICIT NONE

      integer, intent(in) :: NF
! points are the three nodes of the triangle
! point p is the sphere center
      double precision, intent(in), dimension(3,3) :: points
      double precision, intent(in), dimension(dimn) :: pointp
      double precision, intent(out), dimension(dimn) ::  closest_point
! Flag to distinguish point or edge intersection
      logical, intent(out) :: point_or_edge_int
      integer, intent(out) :: edge
! Local variables
      double precision, dimension(3) :: pointa, pointb, pointc
      double precision, dimension(dimn) :: ab, ac, ap, bp,cp
      double precision :: d1, d2, d3, d4, vc, v, d5, d6, vb, w, va, denom

! edge keeps track of the edge id when intersection occurs along an edge
! or a point. When occurring along a point, one of the two edges passing
! through the point is arbitrarily used. When intersection occurs truly
! inside the facet, point_or_edge_int is set to .false., and edge is set
! to 0. Additionally, the x-coordinate of closest_point is set to
! UNDEFINED when edge or point intersection is detected and we want to
! ignore such edge or point intersection. This will lead to a large
! distance between the particle and the intersection point and collision
! will not take place.

      point_or_edge_int = .false.
      edge = 0

      pointa = points(1,:)
      pointb = points(2,:)
      pointc = points(3,:)

      ab = pointb - pointa
      ac = pointc - pointa
      ap = pointp - pointa
      d1 = DOT_PRODUCT(ab, ap)
      d2 = DOT_PRODUCT(ac, ap)

! JFD: All tests <= and >= have been changed to < or > so that when the
! closest point is exactly along an edge, it is marked as belonging to
! the face. Only points that are outside the facet will be marked as
! belonging to and edge or a corner point. This helps in avoiding
! unwanted edge collisions when several triangles are coplanar.
! The last argument point_or_edge_int allows to distinguish between
! points that are inside and truly outside the facet.

! Check if P in vertex region outside A
      IF(d1 < Zero .AND. d2 < zero) then
         closest_point = pointa
         point_or_edge_int = .true.
         if(ignore_stl_edge_intersection(NF,1)) closest_point(1) = UNDEFINED
         edge = 1
         return
      end if

! Check if P in vertex region outside B
      bp = pointp - pointb
      d3 = DOT_PRODUCT(ab, bp);
      d4 = DOT_PRODUCT(ac, bp);
      if (d3 > zero .and. d4 < d3) then
         closest_point = pointb
         point_or_edge_int = .true.
         if(ignore_stl_edge_intersection(NF,1)) closest_point(1) = UNDEFINED
         edge = 1
         return
      endif

! Check if P in edge region of AB, if so return projection of P onto AB
      vc = d1*d4 - d3*d2;
      if (vc < zero .and. d1 > zero .and. d3 < zero) then
         v = d1 / (d1 - d3);
         closest_point =  pointa + v * ab;
         point_or_edge_int = .true.
         if(ignore_stl_edge_intersection(NF,1)) closest_point(1) = UNDEFINED
         edge = 1
         return
      end if

! Check if P in vertex region outside C
      cp = pointp - pointc
      d5 = DOT_PRODUCT(ab, cp)
      d6 = DOT_PRODUCT(ac, cp)
      if (d6 > zero .and. d5 < d6) then
         closest_point  = pointc
         point_or_edge_int = .true.
         if(ignore_stl_edge_intersection(NF,2)) closest_point(1) = UNDEFINED
         edge = 2
         return
      endif

! Check if P in edge region of AC, if so return projection of P onto AC
      vb = d5*d2 - d1*d6

      if (vb < zero .and. d2 > zero .and. d6 < zero) then
         w = d2 / (d2 - d6)
         closest_point = pointa + w * ac
         point_or_edge_int = .true.
         if(ignore_stl_edge_intersection(NF,2)) closest_point(1) = UNDEFINED
         edge = 2
         return
      end if

! Check if P in edge region of BC, if so return projection of P onto BC
      va = d3*d6 - d5*d4
      if (va < zero .and.(d4 - d3) > zero .and. (d5 - d6) > zero) then
         w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
         closest_point = pointb + w * (pointc - pointb)
         point_or_edge_int = .true.
         if(ignore_stl_edge_intersection(NF,3)) closest_point(1) = UNDEFINED
         edge = 3
         return
      end if


! P inside face region. Compute Q through its barycentric coordinates (u,v,w)
      denom = one / (va + vb + vc)
      v = vb * denom
      w = vc * denom
      closest_point = pointa + ab * v + ac * w;
! In 2D, force the z-coordinate to be the same as that of point P
      IF(NO_K) closest_point(3) = pointp(3)
      point_or_edge_int = .false.
      edge = 0
      ignore_stl_edge_intersection(NF,0) = .false.
      return
      end Subroutine ClosestPtPointTriangle

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: RayTriangle                                             !
!                                                                      !
!  Purpose: Return true if ray (origin + direction) intersect triangle.!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
pure function RayTriangle(origin, direction, triangle) result(intersect)

  implicit none

! Ray = origin + direction, triangle 3x3 vertices
  double precision, intent(in) :: origin(3), direction(3)
  double precision, intent(in) :: triangle(3,3)

! Local variables
  double precision :: vertex0(3), vertex1(3), vertex2(3)
  double precision :: edge1(3), edge2(3), h(3), s(3), q(3)

  double precision :: a, f, u, v

  double precision, parameter :: eps = epsilon(0.0)

  logical :: intersect

  intersect = .false.

  vertex0 = triangle(1,:)
  vertex1 = triangle(2,:)
  vertex2 = triangle(3,:)

  edge1 = vertex1 - vertex0
  edge2 = vertex2 - vertex0

  h(1) = direction(2) * edge2(3) - direction(3) * edge2(2)
  h(2) = direction(3) * edge2(1) - direction(1) * edge2(3)
  h(3) = direction(1) * edge2(2) - direction(2) * edge2(1)

  a = dot_product(edge1,h)

  if(a > -eps .and. a < eps) return

  f = 1.0d0/a

  s = origin - vertex0
  u = f*dot_product(s,h)

  if(u < 0.0d0 .or. u > 1.0d0) return

  q(1) = s(2) * edge1(3) - s(3) * edge1(2)
  q(2) = s(3) * edge1(1) - s(1) * edge1(3)
  q(3) = s(1) * edge1(2) - s(2) * edge1(1)

  v = f*dot_product(direction,q)

  if(v < 0.0d0 .or. u+v > 1.0) return

  intersect = .true.

  return

  ! if(t < eps) return
  ! t = f*dot_product(edge2,q)
  ! intersection point :: origin + t*direction

end function RayTriangle


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: intersectLnPlane                                        !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine intersectLnPlane(p_ref, dir_line, ref_plane,       &
         norm_plane, line_param)

      USE discretelement, only: dimn
      USE param1, only: zero

      IMPLICIT NONE

! Reference point and direction of the line
      DOUBLE PRECISION, INTENT(IN) :: P_REF(3),  DIR_LINE(3)
! reference point and normal of the plane
      DOUBLE PRECISION, INTENT(IN) :: REF_PLANE(3), NORM_PLANE(3)

! line is parameterized as p = p_ref + t * dir_line, t is line_param
      double precision, intent(out) :: line_param

      !local vars
      double precision :: denom

      denom = DOT_PRODUCT(dir_line, norm_plane)

      if(denom*denom.gt.zero) then
         line_param = DOT_PRODUCT(ref_plane-p_ref, norm_plane)
         line_param = line_param/denom
      endif

      return
      end subroutine intersectLnPlane

!......................................................................!
! Subroutine TRI_BOX_OVERLAP                                           !
! Author: J.Musser                                   Date: 10-22-2015  !
!                                                                      !
! Purpose: Determine if a box (DES grid cell) intersects the triangle  !
!    (SLT). Note that the DES grid size is slightly increased to       !
!    capture STLs near the boarder of the cell. Otherwise, some        !
!    collisions could ve over looked.                                  !
!                                                                      !
! Author: Tomas Akenine-Moller                   Accessed: 10-22-2015  !
! REF: http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/    !
!         code/tribox2.txt                                             !
!......................................................................!
      SUBROUTINE TRI_BOX_OVERLAP(pCENTER, pHALFSIZE, pVERTS, pOVERLAP)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: pCENTER(3), pHALFSIZE(3)
      DOUBLE PRECISION, INTENT(IN) :: pVERTS(3,3)
      LOGICAL, INTENT(OUT) :: pOVERLAP

      DOUBLE PRECISION :: v0(3), v1(3), v2(3)
      DOUBLE PRECISION :: fex, fey, fez
      DOUBLE PRECISION :: normal(3), e0(3), e1(3), e2(3)

      pOVERLAP = .FALSE.

      v0 = pVERTS(1,:) - pCENTER
      v1 = pVERTS(2,:) - pCENTER
      v2 = pVERTS(3,:) - pCENTER

      e0 = v1-v0
      e1 = v2-v1
      e2 = v0-v2

      fex = abs(e0(1))
      fey = abs(e0(2))
      fez = abs(e0(3))

      if(ATEST_X01(e0(3),e0(2),fez,fey)) return
      if(ATEST_Y02(e0(3),e0(1),fez,fex)) return
      if(ATEST_Z12(e0(2),e0(1),fey,fex)) return

      fex = abs(e1(1))
      fey = abs(e1(2))
      fez = abs(e1(3))

      if(ATEST_X01(e1(3),e1(2),fez,fey)) return
      if(ATEST_Y02(e1(3),e1(1),fez,fex)) return
      if(ATEST_Z0 (e1(2),e1(1),fey,fex)) return

      fex = abs(e2(1))
      fey = abs(e2(2))
      fez = abs(e2(3))

      if(ATEST_X2 (e2(3),e2(2),fez,fey)) return
      if(ATEST_Y1 (e2(3),e2(1),fez,fex)) return
      if(ATEST_Z12(e2(2),e2(1),fey,fex)) return

      if(findMin(v0(1),v1(1),v2(1)) > phalfsize(1) .OR. &
         findMax(v0(1),v1(1),v2(1)) <-phalfsize(1)) return

      if(findMin(v0(2),v1(2),v2(2)) > phalfsize(2) .OR. &
         findMax(v0(2),v1(2),v2(2)) <-phalfsize(2)) return

      if(findMin(v0(3),v1(3),v2(3)) > phalfsize(3) .OR. &
         findMax(v0(3),v1(3),v2(3)) <-phalfsize(3)) return


      normal(1) = e0(2)*e1(3)-e0(3)*e1(2)
      normal(2) = e0(3)*e1(1)-e0(1)*e1(3)
      normal(3) = e0(1)*e1(2)-e0(2)*e1(1)

      if(.NOT.planeBoxOverlap(normal,v0,phalfsize)) return

      pOVERLAP = .TRUE.

      RETURN

      CONTAINS

!``````````````````````````````````````````````````````````````````````!
! Function: planeBoxOverlap                                            !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION planeBoxOverlap(norm, vert, maxbox)

      double precision :: norm(3), vert(3), maxbox(3)

      integer :: lc
      double precision :: vmin(3), vmax(3), v

      do lc=1,3
         v=vert(lc)
         if(norm(lc) > 0.0d0) then
            vmin(lc) = -maxbox(lc) - v
            vmax(lc) =  maxbox(lc) - v
         else
            vmin(lc) = maxbox(lc) - v
            vmax(lc) =-maxbox(lc) - v
         endif
      enddo

      if(dot_product(norm,vmin) > 0.0d0) then
         planeBoxOverlap = .false.
         return
      elseif(dot_product(norm,vmax) >= 0.0d0) then
         planeBoxOverlap = .true.
         return
      endif

      planeBoxOverlap = .false.
      return

      RETURN
      END FUNCTION planeBoxOverlap

!``````````````````````````````````````````````````````````````````````!
! Function: findMin                                                    !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      DOUBLE PRECISION FUNCTION findMin(x0,x1,x2)

      double precision :: x0,x1,x2

      findMin = x0

      if(x1<findMin) findMin=x1
      if(x2<findMin) findMin=x2

      RETURN
      END FUNCTION findMin

!``````````````````````````````````````````````````````````````````````!
! Function: findMax                                                    !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      DOUBLE PRECISION FUNCTION findMax(x0,x1,x2)

      double precision :: x0,x1,x2

      findMax = x0

      if(x1>findMax) findMax=x1
      if(x2>findMax) findMax=x2

      RETURN
      END FUNCTION findMax

!``````````````````````````````````````````````````````````````````````!
! Function: ATEST_X01                                                  !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_X01(a,b,fa,fb)

      double precision :: a, b, fa, fb
      double precision :: lMin, lMax, p0, p2, rad

      p0 = a*v0(2) - b*v0(3)
      p2 = a*v2(2) - b*v2(3)

      if(p0<p2) then; lMIN=p0; lMAX=p2
      else; lMIN=p2; lMAX=p0; endif

      rad=fa*phalfsize(2) + fb*phalfsize(3)
      ATEST_X01=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_X01

!``````````````````````````````````````````````````````````````````````!
! Function: ATEST_X2                                                   !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_X2(a,b,fa,fb)

      double precision :: a, b, fa, fb
      double precision :: lMin, lMax, p0, p1, rad

      p0 = a*v0(2) - b*v0(3)
      p1 = a*v1(2) - b*v1(3)

      if(p0<p1) then; lMIN=p0; lMAX=p1
      else; lMIN=p1; lMAX=p0; endif

      rad=fa*phalfsize(2) + fb*phalfsize(3)
      ATEST_X2=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_X2

!``````````````````````````````````````````````````````````````````````!
! Function: ATEST_Y02                                                  !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_Y02(a,b,fa,fb)

      double precision :: a, b, fa, fb
      double precision :: lMin, lMax, p0, p2, rad

      p0 = -a*v0(1) + b*v0(3)
      p2 = -a*v2(1) + b*v2(3)

      if(p0<p2) then; lMIN=p0; lMAX=p2
      else; lMIN=p2; lMAX=p0; endif

      rad=fa*phalfsize(1) + fb*phalfsize(3)
      ATEST_Y02=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_Y02

!``````````````````````````````````````````````````````````````````````!
! Function: ATEST_Y1                                                   !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_Y1(a,b,fa,fb)

      double precision :: a, b, fa, fb
      double precision :: lMin, lMax, p0, p1, rad

      p0 = -a*v0(1) + b*v0(3)
      p1 = -a*v1(1) + b*v1(3)

      if(p0<p1) then; lMIN=p0; lMAX=p1
      else; lMIN=p1; lMAX=p0; endif

      rad=fa*phalfsize(1) + fb*phalfsize(3)
      ATEST_Y1=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_Y1

!``````````````````````````````````````````````````````````````````````!
! Function: ATEST_Z12                                                  !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_Z12(a,b,fa,fb)

      double precision :: a, b, fa, fb
      double precision :: lMin, lMax, p1, p2, rad

      p1 = a*v1(1) - b*v1(2)
      p2 = a*v2(1) - b*v2(2)

      if(p2<p1) then; lMIN=p2; lMAX=p1
      else; lMIN=p1; lMAX=p2; endif

      rad=fa*phalfsize(1) + fb*phalfsize(2)
      ATEST_Z12=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_Z12

!``````````````````````````````````````````````````````````````````````!
! Function: ATEST_Z0                                                   !
! Support function for TRI_BOX_OVERLAP                                 !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_Z0(a,b,fa,fb)

      double precision :: a, b, fa, fb
      double precision :: lMin, lMax, p0, p1, rad

      p0 = a*v0(1) - b*v0(2)
      p1 = a*v1(1) - b*v1(2)

      if(p0<p1) then; lMIN=p0; lMAX=p1
      else; lMIN=p1; lMAX=p0; endif

      rad=fa*phalfsize(1) + fb*phalfsize(2)
      ATEST_Z0=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_Z0

      END SUBROUTINE TRI_BOX_OVERLAP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_IF_PARTICLE_OVERLAPS_STL                          !
!  Authors: Rahul Garg                               Date: 21-Mar-2014 !
!                                                                      !
!  Purpose: This subroutine is special written to check if a particle  !
!          overlaps any of the STL faces. The routine exits on         !
!          detecting an overlap. It is called after initial            !
!          generation of lattice configuration to remove out of domain !
!          particles                                                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_IF_PARTICLE_OVERLAPS_STL(POS, fI, fJ, fK, REMOVE)

      USE compar
      USE constant
      USE cutcell
      USE desgrid
      USE discretelement, ONLY: dimn, xe, yn, zt, max_radius
      USE functions
      USE geometry
      USE indices
      USE param1
      USE run
      USE stl

      Implicit none

      DOUBLE PRECISION, INTENT(IN) :: POS(DIMN)
      INTEGER, INTENT(IN) :: fI, fJ, fK
      LOGICAL, INTENT(OUT) :: REMOVE

! Integers mapping the fluid cell corners to DES Grid indices.
      INTEGER :: I1, I2, J1, J2, K1, K2

      INTEGER I, J, K, IJK, NF, LC

      DOUBLE PRECISION :: LINE_T
      DOUBLE PRECISION :: RADSQ, DIST(3)

      double precision :: vv(3)

      REMOVE = .TRUE.

      I1 = IofPOS(XE(fI-1))
      I2 = IofPOS(XE( fI ))

      J1 = JofPOS(YN(fJ-1))
      J2 = JofPOS(YN( fJ ))

      K1 = KofPOS(ZT(fK-1))
      K2 = KofPOS(ZT( fK ))

      RADSQ = (1.05d0*MAX_RADIUS)**2

      DO K = K1, K2
      DO J = J1, J2
      DO I = I1, I2

         IF(.NOT.DG_is_ON_myPE_plus1layers(I,J,K)) CYCLE

         IJK = DG_FUNIJK(I,J,K)

! The point is on the non-fluid side of the plane if t>0
         DO LC = 1, FACETS_AT_DG(IJK)%COUNT
            NF = FACETS_AT_DG(IJK)%ID(LC)

            vv = VERTEX(1,:,NF)
            CALL INTERSECTLNPLANE(POS, NORM_FACE(:,NF), &
               vv, NORM_FACE(:,NF), LINE_T)

! Orthogonal projection puts the point outside of the domain or less than
! one particle radius to the facet.
            DIST = LINE_T*NORM_FACE(:,NF)
            IF(LINE_T > ZERO .OR. dot_product(DIST,DIST)<=RADSQ) RETURN

         ENDDO
      ENDDO
      ENDDO
      ENDDO

      REMOVE = .FALSE.

      RETURN
      END SUBROUTINE CHECK_IF_PARTICLE_OVERLAPS_STL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_IF_PARCEL_OVERLAPS_STL                            !
!                                                                      !
!  Purpose: This subroutine is special written to check if a particle  !
!          overlaps any of the STL faces. The routine exits on         !
!          detecting an overlap. It is called after initial            !
!          generation of lattice configuration to remove out of domain !
!          particles                                                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
function picSTLoverlap(ijk, pos, radius) result(overlap_exists)

  use cutcell, only: cartesian_grid
  use stl, only: facets_at_dg, vertex, norm_face

  Implicit none

  integer, intent(in) :: ijk
  double precision, intent(in) :: pos(3), radius

  logical :: overlap_exists, point_or_edge_int
  integer :: edge

  integer :: facet, lc

! line is parameterized as p = p_ref + t * dir_line, t is line_param
  double precision :: dist(3), closest_pt(3), NW_dist



  if(cartesian_grid) then ! this method works best with regular STL files
                          ! (not aligned with the background grid).

  if(facets_at_dg(ijk)%count > 0) then

     overlap_exists = .true.

     do lc = 1, facets_at_dg(ijk)%count
        facet = facets_at_dg(ijk)%id(lc)

        call ClosestPtPointTriangle(facet, pos, vertex(:,:,facet), closest_pt, point_or_edge_int, edge)

        ! dist = closest_pt - pos
        dist = pos - closest_pt
        Nw_dist = dot_product(dist,norm_face(:,facet))

        if(Nw_dist < radius) return

     enddo
  endif
  overlap_exists = .false.

  else  ! This method works better with rectangular geometry where the STL
        ! triangles are aligned with the grid.

! Assume there is no overlap. We loop over all facets (stl triangles), and we
! will return if overlap is detected with a triangle.
  overlap_exists = .false.

  if(facets_at_dg(ijk)%count > 0) then
fcl: do lc = 1, facets_at_dg(ijk)%count
        facet = facets_at_dg(ijk)%id(lc)

        call ClosestPtPointTriangle(facet, pos, vertex(:,:,facet), closest_pt, point_or_edge_int, edge)

! JFD: It is possible that a triangle is in contact with a cell only
! through an edge. In that case, the projected point (closest_pt) would be
! outside of the cell and along the shared edge and we should not use the triangle.
! This could falsely make the code skip the parcel if the dot product below is negative
! (it happens if we use a rectangular block of NSW cells).
! For the PIC silane pyrolysis tutorial, it appears sufficient to cycle on point
! or edge intersections.
! If the intersection is outside the cell and not along a triangle edge, then we
! would need to add a test. This will require adding the two cell corner
! coordinates as arguments to the function, and cycling if the closest point
! is outside these bounds.

        if (point_or_edge_int) cycle fcl

! pos is the location of the parcel,
! closest_pt is the projected point onto the triangle.
! dist goes from the triangle to the parcel.
! norm_face points towards the fluid region.
! If the dot product is positive, the parcel is in the fluid region,
! if the dot product is negative, the parcel is outside the fluid region and
! should be skipped.
! Nw_dist is the signed distance from the wall. Is Nw_dist is less than the
! parcel's radius, then there is overlap with the stl triangle and the parcel is
! skipped.
! Testing (Nw_dist < radius) is sufficient to rule out parcels that are either
! outside the fluid region or overlapping with a triangle.
        dist = pos - closest_pt
        Nw_dist = dot_product(dist,norm_face(:,facet))

        if(Nw_dist < radius) then
           overlap_exists = .true.
           return ! No need to check other facets
        endif

     enddo fcl

  endif

  endif  ! Cartesian grid
  RETURN
end function picSTLoverlap

      END MODULE STL_FUNCTIONS_DES
