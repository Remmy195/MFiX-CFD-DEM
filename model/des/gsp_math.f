#include "error.inc"
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: GSP_MATH_MOD                                           !
!  Author: Renjie Ke                                 Date: 18-Mar-2024 !
!                                                                      !
!  Purpose: Mathematical calculation for gsp model                     !
!           Rotational update subroutines for gsp model                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
MODULE GSP_MATH_MOD
    use error_manager
    type :: gsp_result_type
      double precision :: coords(3)
      double precision :: sphere_radius
    end type gsp_result_type

    CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: gsp_no_squish_rotate(3,conjqm,q,inertia,dtq)           !
!  Reviewer: Renjie Ke                               Date: 18-Mar-2024 !
!                                                                      !
!  Purpose: rigid body rotation,                                       !
!           distance between spheres inside a glued particle           !
!           remain constant                                            !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

    SUBROUTINE gsp_no_squish_rotate(k, p, q, inertia, dt)

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Variables
!-----------------------------------------------
      INTEGER, INTENT(IN)::k
      DOUBLE PRECISION, INTENT (IN) :: inertia(3),dt
      DOUBLE PRECISION, INTENT (INOUT) :: p(4),q(4)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      DOUBLE PRECISION :: phi, c_phi, s_phi, kp(4), kq(4)

!-----------------------------------------------
! apply permutation operator on p and q, get kp and kq
! see reference: "TF Miller III (2002), symplectic quaternion scheme for biophysics molecular dynamics"
! kq and kp refers to columns of S(q) in the paper
! instead of formulated a matrix, LLU did a column by column vector instead
! here is the S(q) =
!                   | q0 -q1 -q2 -q3 |
!                   | q1  q0 -q3  q2 |
!                   | q2  q3  q0 -q1 |
!                   | q3 -q2  q1  q0 |
! q^(dot) = 0.5 * S(q) * omega^(4)
!
! In Eqn. 2.18 it explains the a sum over permutation matrices
! That is the operation LLU doing here "IMPORTANT!"
      IF (k == 1) THEN
         kq(1) = -q(2);  kp(1) = -p(2)
         kq(2) =  q(1);  kp(2) =  p(1)
         kq(3) =  q(4);  kp(3) =  p(4)
         kq(4) = -q(3);  kp(4) = -p(3)
      ELSEIF (k == 2) THEN
         kq(1) = -q(3);  kp(1) = -p(3)
         kq(2) = -q(4);  kp(2) = -p(4)
         kq(3) =  q(1);  kp(3) =  p(1)
         kq(4) =  q(2);  kp(4) =  p(2)
      ELSEIF (k == 3) THEN
         kq(1) = -q(4);  kp(1) = -p(4)
         kq(2) =  q(3);  kp(2) =  p(3)
         kq(3) = -q(2);  kp(3) = -p(2)
         kq(4) =  q(1);  kp(4) =  p(1)
      ENDIF

! obtain phi, cosines and sines
      phi = p(1) * kq(1) + p(2) * kq(2) + p(3) * kq(3) + p(4) * kq(4)
      phi = phi / (4.0 * inertia(k))
      c_phi = cos(dt * phi)
      s_phi = sin(dt * phi)

! update p and q, output as conjqm and q
      p = c_phi * p + s_phi * kp
      q = c_phi * q + s_phi * kq

      RETURN

    END SUBROUTINE gsp_no_squish_rotate

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: gsp_invquatvec(q, conjqm, mbody)                       !
!  Reviewer: Renjie Ke                               Date: 18-Mar-2024 !
!  Purpose: inverse of quaternions vector q                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

   SUBROUTINE gsp_invquatvec(a, b, c)

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Variables
!-----------------------------------------------
      DOUBLE PRECISION, INTENT (IN) :: a(4), b(4)
      DOUBLE PRECISION, INTENT (OUT) :: c(3)

!-----------------------------------------------
      c(1) = -a(2)*b(1) + a(1)*b(2) + a(4)*b(3) - a(3)*b(4)
      c(2) = -a(3)*b(1) - a(4)*b(2) + a(1)*b(3) + a(2)*b(4)
      c(3) = -a(4)*b(1) + a(3)*b(2) - a(2)*b(3) + a(1)*b(4)

      RETURN

   END SUBROUTINE gsp_invquatvec

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: gsp_quat_to_exyz(q,ex,ey,ez)                           !
!  Reviewer: Renjie Ke                               Date: 18-Mar-2024 !
!                                                                      !
!  Purpose: compute space-frame ex,ey,ez from current quaternion q     !
!           ex,ey,ez = space-frame coords of 1st,2nd,3rd principal axis!
!           For example,                                               !
!           operation is ex = q' d q = Q d,                            !
!           where d is (1,0,0) = 1st axis in body frame                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE gsp_quat_to_exyz(q, ex, ey, ez)

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Variables
!-----------------------------------------------
      DOUBLE PRECISION, INTENT (IN) :: q(:)
      DOUBLE PRECISION, INTENT (INOUT) :: ex(:),ey(:),ez(:)

!-----------------------------------------------
      ex(1) = q(1) * q(1) + q(2) * q(2) - q(3) * q(3) - q(4) * q(4)
      ex(2) = 2.d0 * (q(2) * q(3) + q(1) * q(4))
      ex(3) = 2.d0 * (q(2) * q(4) - q(1) * q(3))

      ey(1) = 2.d0 * (q(2) * q(3) - q(1) * q(4))
      ey(2) = q(1) * q(1) - q(2) * q(2) + q(3) * q(3) - q(4) * q(4)
      ey(3) = 2.d0 * (q(3) * q(4) + q(1) * q(2))

      ez(1) = 2.d0 * (q(2) * q(4) + q(1) * q(3))
      ez(2) = 2.d0 * (q(3) * q(4) - q(1) * q(2))
      ez(3) = q(1) * q(1) - q(2) * q(2) - q(3) * q(3) + q(4) * q(4)

! normalize the ex, ey, ez to unit vector in space frame
      ex = ex / sqrt(dot_product(ex, ex))
      ey = ey / sqrt(dot_product(ey, ey))
      ez = ez / sqrt(dot_product(ez, ez))

      RETURN

   END SUBROUTINE gsp_quat_to_exyz

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: gsp_q_from_exyz(ex,ey,ez,q)                            !
!  Reviewer: Renjie Ke                               Date: 18-Mar-2024 !
!                                                                      !
!  Purpose: create unit quaternion from space-frame ex,ey,ez           !
!           ex,ey,ez are columns of a rotation matrix                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE gsp_q_from_exyz(ex, ey, ez, q)

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Variables
!-----------------------------------------------
      DOUBLE PRECISION, INTENT(IN)::ex(3), ey(3), ez(3)
      DOUBLE PRECISION, INTENT(OUT):: q(4)

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      DOUBLE PRECISION :: q0sq, q1sq, q2sq, q3sq, norm

!-----------------------------------------------
! squares of quaternion components
! q0^2 + q1^2 + q2^2 + q3^2 = 1 --> q0^2 = 1 - (q1^2 + q2^2 + q3^2)
! For rotation matrix R
! ex(1) is R_11, ey(2) = R_22, ez(3) = R_33
! R_11 = 1 - 2(q2^2 + q3^2), R_22 = 1 - 2(q1^2 + q3^2), R_33 = 1 - 2(q3^2 + q3^2)
! So q0^2 is obtained as following:
      q0sq = 0.25d0 * (ex(1) + ey(2) + ez(3) + 1.0)
      q1sq = q0sq - 0.5 * (ey(2) + ez(3))
      q2sq = q0sq - 0.5 * (ex(1) + ez(3))
      q3sq = q0sq - 0.5 * (ex(1) + ey(2))

!some component must be greater than 1/4 since they sum to
!compute other components from it
      IF (q0sq >= 0.25) THEN
         q(1) = sqrt(q0sq)
         q(2) = (ey(3) - ez(2)) / (4.0 * q(1))
         q(3) = (ez(1) - ex(3)) / (4.0 * q(1))
         q(4) = (ex(2) - ey(1)) / (4.0 * q(1))
      ELSEIF (q1sq >= 0.25) THEN
         q(2) = sqrt(q1sq)
         q(1) = (ey(3) - ez(2)) / (4.0 * q(2))
         q(3) = (ey(1) + ex(2)) / (4.0 * q(2))
         q(4) = (ex(3) + ez(1)) / (4.0 * q(2))
      ELSEIF (q2sq >= 0.25) THEN
         q(3) = sqrt(q2sq)
         q(1) = (ez(1) - ex(3)) / (4.0 * q(3))
         q(2) = (ey(1) + ex(2)) / (4.0 * q(3))
         q(4) = (ez(2) + ey(3)) / (4.0 * q(3))
      ELSEIF (q3sq >= 0.25) THEN
         q(4) = sqrt(q3sq)
         q(1) = (ex(2) - ey(1)) / (4.0 * q(4))
         q(2) = (ez(1) + ex(3)) / (4.0 * q(4))
         q(3) = (ez(2) + ey(3)) / (4.0 * q(4))
      ENDIF

      norm = 1.d0 / sqrt(q(1) * q(1) + q(2) * q(2) + q(3) * q(3) + q(4) * q(4))
      q(1:4) = q(1:4) * norm

      RETURN

   END SUBROUTINE gsp_q_from_exyz

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: gsp_omega_to_angmom(quat, omega, inertia, angmom)      !
!  Reviewer: Renjie Ke                               Date: 18-Mar-2024 !
!                                                                      !
!  Purpose: from rotational velocity omege,                            !
!           calculate the angular momentum                             !
!           not used here                                              !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE gsp_omega_to_angmom(quat, omega, inertia, angmom)

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Variables
!-----------------------------------------------
      DOUBLE PRECISION, INTENT (IN) :: quat(4), omega(3), inertia(3)
      DOUBLE PRECISION, INTENT (OUT) ::angmom(3)

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      DOUBLE PRECISION :: wbody(3), angmombody(3), rot(3,3)

!-----------------------------------------------
      CALL gsp_quat_to_mat(quat,rot)
      CALL gsp_transpose_matvec(rot,omega,wbody)

      angmombody(1:3) = wbody(1:3) * inertia(1:3)
      CALL gsp_matvec(rot,angmombody,angmom)

      RETURN

   END SUBROUTINE gsp_omega_to_angmom

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: gsp_quat_to_mat(quat,mat)                              !
!  Reviewer: Renjie Ke                               Date: 18-Mar-2024 !
!                                                                      !
!  Purpose: from quaternion to calculate the rotational matrix,        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE gsp_quat_to_mat(quat, mat)

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Variables
!-----------------------------------------------
      DOUBLE PRECISION, INTENT (IN) :: quat(4)
      DOUBLE PRECISION, INTENT (OUT) :: mat(3,3)

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      DOUBLE PRECISION :: w2, i2, j2, k2
      DOUBLE PRECISION :: twoij, twoik, twojk, twoiw, twojw, twokw

!-----------------------------------------------
! For example,
! R(1,1) = 1 - 2*(q2^2 + q3^3) --> R(1,1) = q0^2 + q1^2 - q2^2 - q3^2 for summation of qi^2 = 1
! individual term formulation
      w2 = quat(1)*quat(1);
      i2 = quat(2)*quat(2);
      j2 = quat(3)*quat(3);
      k2 = quat(4)*quat(4);
      twoij = 2.0*quat(2)*quat(3);
      twoik = 2.0*quat(2)*quat(4);
      twojk = 2.0*quat(3)*quat(4);
      twoiw = 2.0*quat(2)*quat(1);
      twojw = 2.0*quat(3)*quat(1);
      twokw = 2.0*quat(4)*quat(1);

! formulate the rotational matrix by assigning values
      mat(1,1) = w2+i2-j2-k2;
      mat(1,2) = twoij-twokw;
      mat(1,3) = twojw+twoik;

      mat(2,1) = twoij+twokw;
      mat(2,2) = w2-i2+j2-k2;
      mat(2,3) = twojk-twoiw;

      mat(3,1) = twoik-twojw;
      mat(3,2) = twojk+twoiw;
      mat(3,3) = w2-i2-j2+k2;

      RETURN

   END SUBROUTINE gsp_quat_to_mat

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: gsp_quatvec(a,b,c)                                     !
!               gsp_qnormalize(q)                                      !
!               gsp_matvec(m,v,answer)                                 !
!               gsp_transpose_matvec(m,v,answer)                       !
!  Reviewer: Renjie Ke                               Date: 18-Mar-2024 !
!                                                                      !
!  Purpose: vectors and matrices operation                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

   SUBROUTINE gsp_quatvec(a,b,c)
! purpose: calculate Q 'dot_product' Omega
! | q0, -q1, -q2, -q3 |         | 0  |
! | q1,  q0, -q3,  q2 |   dot   | Wx |
! | q2,  q3,  q0, -q1 | product | Wy |
! | q3, -q2,  q1,  q0 |         | Wz |
! b(3) refers to Wx, Wy, Wz
! a(4) refers to q0, q1, q2, q3
      IMPLICIT NONE
         ! Dummy variables
      DOUBLE PRECISION, INTENT (IN) :: a(4), b(3)
      DOUBLE PRECISION, INTENT (OUT) :: c(4)

         !------------------------------------------
         c(1) = -a(2)*b(1) - a(3)*b(2) - a(4)*b(3)
         c(2) = a(1)*b(1) + a(3)*b(3) - a(4)*b(2)
         c(3) = a(1)*b(2) + a(4)*b(1) - a(2)*b(3)
         c(4) = a(1)*b(3) + a(2)*b(2) - a(3)*b(1)

      RETURN

   END SUBROUTINE gsp_quatvec

   SUBROUTINE gsp_qnormalize(q)

      IMPLICIT NONE
      ! Dummy variables
      DOUBLE PRECISION, INTENT (INOUT) :: q(4)
      ! Local variables
      DOUBLE PRECISION:: invnorm

      !------------------------------------------
      invnorm = 1.0/sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3)+q(4)*q(4))
      q(1:4) = q(1:4)*invnorm

      RETURN

   END SUBROUTINE gsp_qnormalize


   SUBROUTINE gsp_matvec(m,v,answer)
! purpose: a matrix dot product a vector
      IMPLICIT NONE
      ! Dummy variables
      DOUBLE PRECISION, INTENT (IN) :: m(3,3), v(3)
      DOUBLE PRECISION, INTENT (OUT) :: answer(3)

      !------------------------------------------
      answer(1) = m(1,1)*v(1) + m(1,2)*v(2) + m(1,3)*v(3)
      answer(2) = m(2,1)*v(1) + m(2,2)*v(2) + m(2,3)*v(3)
      answer(3) = m(3,1)*v(1) + m(3,2)*v(2) + m(3,3)*v(3)

      RETURN

   END SUBROUTINE gsp_matvec

   SUBROUTINE gsp_transpose_matvec(m,v,answer)
! purpose: a transpose of a matrix dot product a vector
      IMPLICIT NONE
      ! Dummy variables
      DOUBLE PRECISION, INTENT (IN) :: m(3,3), v(3)
      DOUBLE PRECISION, INTENT (OUT) :: answer(3)

      !------------------------------------------
      answer(1) = m(1,1)*v(1) + m(2,1)*v(2) + m(3,1)*v(3)
      answer(2) = m(1,2)*v(1) + m(2,2)*v(2) + m(3,2)*v(3)
      answer(3) = m(1,3)*v(1) + m(2,3)*v(2) + m(3,3)*v(3)

      RETURN

   END SUBROUTINE gsp_transpose_matvec

   SUBROUTINE SET_MARKED_GSP(PP)
        USE resize, ONLY: INTEGER_GROW
        USE discretelement, ONLY: MARKED_GSP, NGluedParticles
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: PP
        INTEGER :: L, SIZE_OLD
        DO L = 1, SIZE(MARKED_GSP)
            IF(MARKED_GSP(L) == 0) EXIT
        ENDDO
        ! non zero found
        IF(L .gt. SIZE(MARKED_GSP)) THEN
            SIZE_OLD = SIZE(MARKED_GSP)
            IF(SIZE_OLD < NGluedParticles) CALL INTEGER_GROW(MARKED_GSP,2*SIZE_OLD)
            MARKED_GSP(SIZE_OLD+1:2*SIZE_OLD) = 0
            MARKED_GSP(SIZE_OLD+1) = PP
        ELSE
        ! zero found
            MARKED_GSP(L) = PP
        ENDIF
   END SUBROUTINE SET_MARKED_GSP

   SUBROUTINE SET_MARKED_GSP_STATE(PP)
        USE resize, ONLY: INTEGER_GROW
        USE discretelement, ONLY: marked_state, NGluedParticles
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: PP
        INTEGER :: L, SIZE_OLD
        DO L = 1, SIZE(marked_state)
            IF(marked_state(L) == -1) EXIT
        ENDDO
        ! non negative one found
        IF(L .gt. SIZE(marked_state)) THEN
            SIZE_OLD = SIZE(marked_state)
            IF(SIZE_OLD < NGluedParticles) CALL INTEGER_GROW(marked_state,2*SIZE_OLD)
            marked_state(SIZE_OLD+1:2*SIZE_OLD) = -1
            marked_state(SIZE_OLD+1) = PP
        ELSE
        ! negative one found
            marked_state(L) = PP
        ENDIF
   END SUBROUTINE SET_MARKED_GSP_STATE

   function findIndex_1i(array, value) result(index)
   ! 1i refers to 1d integer array
      integer, intent(in) :: array(:)
      integer, intent(in) :: value
      integer :: index
      integer :: i

      index = -1
      do i = 1, size(array)
         if (array(i) == value) then
            index = i
            exit
         endif
      end do
   end function findIndex_1i

! @renjieke, a special data type to pack coordinates and radius together.
   function calc_coords_radius_component(PP, MM, L) result(res)
      use discretelement, only: cspart_data, des_radius

      integer, intent(in) :: MM
      integer, intent(in) :: L
      integer, intent(in) :: PP

      type(gsp_result_type) :: res
      double precision :: actual_bounding

      actual_bounding = 2.0 * des_radius(L)
      res%sphere_radius = cspart_data(PP,4,MM) / 2.0 * actual_bounding
      res%coords(1) = cspart_data(PP,1,MM) * actual_bounding
      res%coords(2) = cspart_data(PP,2,MM) * actual_bounding
      res%coords(3) = cspart_data(PP,3,MM) * actual_bounding

   end function calc_coords_radius_component

END MODULE GSP_MATH_MOD