#include "error.inc"

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Module name: SQ_ROTATION_MOD                                     !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Rotate vector between global and local coordinates                 !
!            and renew quaternion                                             !
!                                                                             !
! Reference:                                                                  !
!   Xi Gao, Jia Yu, Ricardo JF Portal, Jean-François Dietiker,                !
!   Mehrdad Shahnam and William A Rogers,Development and validation           !
!   of SuperDEM for non-spherical particulate systems using a                 !
!   superquadric particle method, Particuology,                               !
!   2021: doi.org/10.1016/j.partic.2020.1011.1007.                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
MODULE SQ_ROTATION_MOD



CONTAINS

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: QROTATE                                                   !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Rotate vector between global and local coordinates.               !
!      This subroutine uses a quaternion Qc to rotate a vector "X" between    !
!       a global and local frame.                                             !
!                                                                             !
!      When ixg=1, then the quaternion Qc is used to rotate a vector from     !
!       the global frame to the local frame: X_global to X_local              !
!                                                                             !
!       When ixg=2, the inverse of Qc is used to rotate a vector from         !
!       the local frame to the global frame: X_local to X_global              !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE QROTATE(Qc,X_global,X_local,Ixg)
      IMPLICIT NONE

      INTEGER*4 :: Ixg
      DOUBLE PRECISION :: Qc, qrot, qirot, tt2, tt3, tt4, tt5, tt6,&
                          tt7, tt8, tt9, tt10, X_global, X_local
      DIMENSION :: Qc(4), qrot(3,3), qirot(3,3), X_global(3), X_local(3)
! Products for finding the rotation matrix of Qc
      tt2 = Qc(1)*Qc(2)
      tt3 = Qc(1)*Qc(3)
      tt4 = Qc(1)*Qc(4)
      tt5 = -Qc(2)*Qc(2)
      tt6 = Qc(2)*Qc(3)
      tt7 = Qc(2)*Qc(4)
      tt8 = -Qc(3)*Qc(3)
      tt9 = Qc(3)*Qc(4)
      tt10 = -Qc(4)*Qc(4)
! Rotation matrix of quaternion Qc, which rotates vectors from the global
! frame to the local frame
      qrot(1,1) = tt8 + tt10 + 0.5D0
      qrot(1,2) = tt6 + tt4
      qrot(1,3) = tt7 - tt3
      qrot(2,1) = tt6 - tt4
      qrot(2,2) = tt5 + tt10 + 0.5D0
      qrot(2,3) = tt9 + tt2
      qrot(3,1) = tt7 + tt3
      qrot(3,2) = tt9 - tt2
      qrot(3,3) = tt5 + tt8 + 0.5D0
!
      IF ( Ixg.EQ.1 ) THEN
! Rotate the vector "X_global" from the global frame to the
!         local frame
        X_local(1) = 2.D0*(qrot(1,1)*X_global(1)+qrot(1,2)*X_global(2) &
                    & +qrot(1,3)*X_global(3))
         X_local(2) = 2.D0*(qrot(2,1)*X_global(1)+qrot(2,2)*X_global(2) &
                    & +qrot(2,3)*X_global(3))
         X_local(3) = 2.D0*(qrot(3,1)*X_global(1)+qrot(3,2)*X_global(2) &
                    & +qrot(3,3)*X_global(3))
      ELSEIF ( Ixg.EQ.2 ) THEN
! Rotate the local vector "X_local" from the local frame to the
!         global frame
!
! Inverse rotation matrix QiRot: the rotation matrix of the
!         conjugate (inverse) of quaternion Qci, and it rotates
!         vectors from the local frame to the (global) frame
         qirot(1,1) = qrot(1,1)
         qirot(1,2) = qrot(2,1)
         qirot(1,3) = qrot(3,1)
         qirot(2,1) = qrot(1,2)
         qirot(2,2) = qrot(2,2)
         qirot(2,3) = qrot(3,2)
         qirot(3,1) = qrot(1,3)
         qirot(3,2) = qrot(2,3)
         qirot(3,3) = qrot(3,3)
! Rotate the vector "X_local" from the local frame to the
!         global frame
         X_global(1) = 2.D0*(qirot(1,1)*X_local(1)+qirot(1,2)*X_local(2)&
                     & +qirot(1,3)*X_local(3))
         X_global(2) = 2.D0*(qirot(2,1)*X_local(1)+qirot(2,2)*X_local(2)&
                     & +qirot(2,3)*X_local(3))
         X_global(3) = 2.D0*(qirot(3,1)*X_local(1)+qirot(3,2)*X_local(2)&
                     & +qirot(3,3)*X_local(3))
      ENDIF
      END SUBROUTINE QROTATE

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: QROTATE                                                   !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Rotate vector between global and local coordinates.               !
!      This subroutine uses a quaternion Qc to rotate a vector "X" between    !
!       a global and local frame.                                             !
!                                                                             !
!      When ixg=1, then the quaternion Qc is used to rotate a vector from     !
!       the global frame to the local frame: X_global to X_local              !
!                                                                             !
!       When ixg=2, the inverse of Qc is used to rotate a vector from         !
!       the local frame to the global frame: X_local to X_global              !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      !! FIXME use quaternion multiplication here instead of converting to matrix!
      SUBROUTINE QMultiPointsROTATE(Qc,X_global,X_local,Ixg,pointsnumber)
      IMPLICIT NONE

      INTEGER*4 :: I,Ixg,pointsnumber
      DOUBLE PRECISION :: Qc, qrot, qirot, tt2, tt3, tt4, tt5, tt6,&
                          tt7, tt8, tt9, tt10, X_global, X_local
      DIMENSION :: Qc(4), qrot(3,3), qirot(3,3)
      DIMENSION :: X_global(3*pointsnumber), X_local(3*pointsnumber)
! Products for finding the rotation matrix of Qc
      tt2 = Qc(1)*Qc(2)
      tt3 = Qc(1)*Qc(3)
      tt4 = Qc(1)*Qc(4)
      tt5 = -Qc(2)*Qc(2)
      tt6 = Qc(2)*Qc(3)
      tt7 = Qc(2)*Qc(4)
      tt8 = -Qc(3)*Qc(3)
      tt9 = Qc(3)*Qc(4)
      tt10 = -Qc(4)*Qc(4)

! Rotation matrix of quaternion Qc, which rotates vectors from the global
! frame to the local frame
      qrot(1,1) = tt8 + tt10 + 0.5D0
      qrot(1,2) = tt6 + tt4
      qrot(1,3) = tt7 - tt3
      qrot(2,1) = tt6 - tt4
      qrot(2,2) = tt5 + tt10 + 0.5D0
      qrot(2,3) = tt9 + tt2
      qrot(3,1) = tt7 + tt3
      qrot(3,2) = tt9 - tt2
      qrot(3,3) = tt5 + tt8 + 0.5D0
!
      IF ( Ixg.EQ.1 ) THEN
! Rotate the vector "X_global" from the global frame to the
!         local frame
        DO i=1,pointsnumber
        X_local((i-1)*3+1) = 2.D0*(qrot(1,1)*X_global((i-1)*3+1)+&
                                   qrot(1,2)*X_global((i-1)*3+2)+&
                                   qrot(1,3)*X_global((i-1)*3+3))
        X_local((i-1)*3+2) = 2.D0*(qrot(2,1)*X_global((i-1)*3+1)+&
                                   qrot(2,2)*X_global((i-1)*3+2)+&
                                   qrot(2,3)*X_global((i-1)*3+3))
        X_local((i-1)*3+3) = 2.D0*(qrot(3,1)*X_global((i-1)*3+1)+&
                                   qrot(3,2)*X_global((i-1)*3+2)+&
                                   qrot(3,3)*X_global((i-1)*3+3))
        ENDDO
      ELSEIF ( Ixg.EQ.2 ) THEN
! Rotate the local vector "X_local" from the local frame to the
!         global frame
!
! Inverse rotation matrix QiRot: the rotation matrix of the
!         conjugate (inverse) of quaternion Qci, and it rotates
!         vectors from the local frame to the (global) frame
         qirot(1,1) = qrot(1,1)
         qirot(1,2) = qrot(2,1)
         qirot(1,3) = qrot(3,1)
         qirot(2,1) = qrot(1,2)
         qirot(2,2) = qrot(2,2)
         qirot(2,3) = qrot(3,2)
         qirot(3,1) = qrot(1,3)
         qirot(3,2) = qrot(2,3)
         qirot(3,3) = qrot(3,3)
! Rotate the vector "X_local" from the local frame to the
!         global frame
         DO i=1,pointsnumber
         X_global((i-1)*3+1) = 2.D0*(qirot(1,1)*X_local((i-1)*3+1)+&
                                     qirot(1,2)*X_local((i-1)*3+2)+&
                                     qirot(1,3)*X_local((i-1)*3+3))
         X_global((i-1)*3+2) = 2.D0*(qirot(2,1)*X_local((i-1)*3+1)+&
                                     qirot(2,2)*X_local((i-1)*3+2)+&
                                     qirot(2,3)*X_local((i-1)*3+3))
         X_global((i-1)*3+3) = 2.D0*(qirot(3,1)*X_local((i-1)*3+1)+&
                                     qirot(3,2)*X_local((i-1)*3+2)+&
                                     qirot(3,3)*X_local((i-1)*3+3))
      ENDDO
      ENDIF
      END SUBROUTINE QMultiPointsROTATE

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: QROTATIONMATRIX                                           !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Only the rotation matrix for quaternions                          !
!    Rotation matrix of quaternion Qc, which rotates vectors from the global  !
!    frame to the local frame (half)                                          !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE QROTATIONMATRIX(Qc,Qrot)
      IMPLICIT NONE

      DOUBLE PRECISION :: Qc , Qrot , tt2 , tt3 , tt4 , tt5 , tt6 , tt7 ,  &
                     & tt8 , tt9 , tt10
      DIMENSION :: Qc(4) , Qrot(3,3)

!  QRot is actually one-half of the
!       true rotation matrix
!
!  Products for finding the rotation matrix of Qc
      tt2 = Qc(1)*Qc(2)
      tt3 = Qc(1)*Qc(3)
      tt4 = Qc(1)*Qc(4)
      tt5 = -Qc(2)*Qc(2)
      tt6 = Qc(2)*Qc(3)
      tt7 = Qc(2)*Qc(4)
      tt8 = -Qc(3)*Qc(3)
      tt9 = Qc(3)*Qc(4)
      tt10 = -Qc(4)*Qc(4)
!
! Rotation matrix of quaternion Qc, which rotates
!       vectors from the (global) frame of the local frame
      Qrot(1,1) = tt8 + tt10 + 0.5D0
      Qrot(1,2) = tt6 + tt4
      Qrot(1,3) = tt7 - tt3
      Qrot(2,1) = tt6 - tt4
      Qrot(2,2) = tt5 + tt10 + 0.5D0
      Qrot(2,3) = tt9 + tt2
      Qrot(3,1) = tt7 + tt3
      Qrot(3,2) = tt9 - tt2
      Qrot(3,3) = tt5 + tt8 + 0.5D0
!
      END SUBROUTINE QROTATIONMATRIX

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: QROTATIONMATRIX_INVERSE                                   !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Only the rotation matrix for quaternions                          !
!      Rotation matrix of quaternion Qc, which rotates vectors from the local !
!      frame to the global frame (half)                                       !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE QROTATIONMATRIX_INVERSE(Qc,Qirot)
      IMPLICIT NONE
      DOUBLE PRECISION :: Qc , Qrot , tt2 , tt3 , tt4 , tt5 , tt6 , tt7 ,  &
                     & tt8 , tt9 , tt10,Qirot
      DIMENSION :: Qc(4) , Qrot(3,3),Qirot(3,3)

!  Products for finding the rotation matrix of Qc
      tt2 = Qc(1)*Qc(2)
      tt3 = Qc(1)*Qc(3)
      tt4 = Qc(1)*Qc(4)
      tt5 = -Qc(2)*Qc(2)
      tt6 = Qc(2)*Qc(3)
      tt7 = Qc(2)*Qc(4)
      tt8 = -Qc(3)*Qc(3)
      tt9 = Qc(3)*Qc(4)
      tt10 = -Qc(4)*Qc(4)
!
! Rotation matrix of quaternion Qc, which rotates
!       vectors from the (global) frame of the local frame
      Qrot(1,1) = tt8 + tt10 + 0.5D0
      Qrot(1,2) = tt6 + tt4
      Qrot(1,3) = tt7 - tt3
      Qrot(2,1) = tt6 - tt4
      Qrot(2,2) = tt5 + tt10 + 0.5D0
      Qrot(2,3) = tt9 + tt2
      Qrot(3,1) = tt7 + tt3
      Qrot(3,2) = tt9 - tt2
      Qrot(3,3) = tt5 + tt8 + 0.5D0

      qirot(1,1) = qrot(1,1)
      qirot(1,2) = qrot(2,1)
      qirot(1,3) = qrot(3,1)
      qirot(2,1) = qrot(1,2)
      qirot(2,2) = qrot(2,2)
      qirot(2,3) = qrot(3,2)
      qirot(3,1) = qrot(1,3)
      qirot(3,2) = qrot(2,3)
      qirot(3,3) = qrot(3,3)
      END SUBROUTINE QROTATIONMATRIX_INVERSE

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: QMATRIXROTATEVECTOR                                       !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Only the rotation vector                                          !
!  This subroutine uses the rotation matrix QRot to rotate a                  !
!       vector "X" between the global frame (1) and local frame (2).          !
!                                                                             !
! NOTE: the matrix QRot is computed from a Quaternion, so QRot is             !
!       actual one-half of the true rotation matrix.  Hence, we multiply      !
!       by a factor of 2.d0 below.                                            !
!                                                                             !
! When ixg=1, then the matrix QRot is used to rotate a vector from            !
!       the global frame (1) to the local frame (2): X_frame_1 to X_frame_2   !
!                                                                             !
! When ixg=2, the inverse of QRot is used to rotate a vector from             !
!       the local frame (2) to the global frame (1): X_frame_2 to X_frame_1   !                                                                         !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE QMATRIXROTATEVECTOR(Qrot,X_frame_1,X_frame_2,Ixg)
      IMPLICIT NONE

      INTEGER*4 :: Ixg
      DOUBLE PRECISION :: Qrot, X_frame_2, X_frame_1
      DIMENSION :: Qrot(3,3), X_frame_2(3), X_frame_1(3)
!
      IF ( Ixg.EQ.1 ) THEN
! Rotate the vector "X_frame_1" from the global frame (1) to the local frame (2)
         X_frame_2(1) = 2.D0*(Qrot(1,1)*X_frame_1(1)+Qrot(1,2)          &
                      & *X_frame_1(2)+Qrot(1,3)*X_frame_1(3))
         X_frame_2(2) = 2.D0*(Qrot(2,1)*X_frame_1(1)+Qrot(2,2)          &
                      & *X_frame_1(2)+Qrot(2,3)*X_frame_1(3))
         X_frame_2(3) = 2.D0*(Qrot(3,1)*X_frame_1(1)+Qrot(3,2)          &
                      & *X_frame_1(2)+Qrot(3,3)*X_frame_1(3))
      ELSEIF ( Ixg.EQ.2 ) THEN
! Rotate the local vector "X_frame_2" from the local frame (2) to the
!         global frame (1)
!
! We use the inverse of the rotation matrix QRot, which is
!         its transpose
!
! Rotate the vector "X_frame_2" from the local frame (2) to the
!         global frame (1)
         X_frame_1(1) = 2.D0*(Qrot(1,1)*X_frame_2(1)+Qrot(2,1)          &
                      & *X_frame_2(2)+Qrot(3,1)*X_frame_2(3))
         X_frame_1(2) = 2.D0*(Qrot(1,2)*X_frame_2(1)+Qrot(2,2)          &
                      & *X_frame_2(2)+Qrot(3,2)*X_frame_2(3))
         X_frame_1(3) = 2.D0*(Qrot(1,3)*X_frame_2(1)+Qrot(2,3)          &
                      & *X_frame_2(2)+Qrot(3,3)*X_frame_2(3))
      ENDIF
      END SUBROUTINE QMATRIXROTATEVECTOR

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !                                                                           !
!  Subroutine name: QMATRIXROTATETENSOR                                       !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Only the rotation tensor                                          !
!  This subroutine uses the rotation matrix QRot to rotate a                  !
!       tensor "X" between the global frame (1) and local frame (2).          !
!                                                                             !
! NOTE: the matrix QRot is computed from a Quaternion, so QRot is             !
!       actual one-half of the true rotation matrix.  Hence, we multiply      !
!       by a factor of 2.d0 below.                                            !
!                                                                             !
! When ixg=1, then the matrix QRot is used to rotate a tensor from            !
!       the global frame (1) to the local frame (2): X_frame_1 to X_frame_2   !
!                                                                             !
! When ixg=2, the inverse of QRot is used to rotate a tensor from             !
!       the local frame (2) to the global frame (1): X_frame_2 to X_frame_1   !                                                                         !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE QMATRIXROTATETENSOR(Qrot,X_frame_1,X_frame_2,Ixg)
      IMPLICIT NONE
      INTEGER*4 :: Ixg
!Qrot is the rotation matrix from global to local
      DOUBLE PRECISION :: Qrot(3,3),Qirot(3,3)
      DOUBLE PRECISION :: X_frame_2(3,3), X_frame_1(3,3),M_1(3,3),M_2(3,3)
! Inverse rotation matrix QiRot: the rotation matrix of the
!         conjugate (inverse) of quaternion Qci, and it rotates
!         vectors from the local frame to the (global) frame
         qirot(1,1) = qrot(1,1)
         qirot(1,2) = qrot(2,1)
         qirot(1,3) = qrot(3,1)
         qirot(2,1) = qrot(1,2)
         qirot(2,2) = qrot(2,2)
         qirot(2,3) = qrot(3,2)
         qirot(3,1) = qrot(1,3)
         qirot(3,2) = qrot(2,3)
         qirot(3,3) = qrot(3,3)
      IF ( Ixg.EQ.1 ) THEN
! Rotate the tensor "X_frame_1" from the global frame (1) to the local frame (2)

         M_2(1,1) = 2.D0*(Qrot(1,1)*X_frame_1(1,1)   &
                         +Qrot(1,2)*X_frame_1(2,1)   &
                         +Qrot(1,3)*X_frame_1(3,1))
         M_2(1,2) = 2.D0*(Qrot(1,1)*X_frame_1(1,2)   &
                         +Qrot(1,2)*X_frame_1(2,2)   &
                         +Qrot(1,3)*X_frame_1(3,2))
         M_2(1,3) = 2.D0*(Qrot(1,1)*X_frame_1(1,3)   &
                         +Qrot(1,2)*X_frame_1(2,3)   &
                         +Qrot(1,3)*X_frame_1(3,3))
         M_2(2,1) = 2.D0*(Qrot(2,1)*X_frame_1(1,1)   &
                         +Qrot(2,2)*X_frame_1(2,1)   &
                         +Qrot(2,3)*X_frame_1(3,1))
         M_2(2,2) = 2.D0*(Qrot(2,1)*X_frame_1(1,2)   &
                         +Qrot(2,2)*X_frame_1(2,2)   &
                         +Qrot(2,3)*X_frame_1(3,2))
         M_2(2,3) = 2.D0*(Qrot(2,1)*X_frame_1(1,3)   &
                         +Qrot(2,2)*X_frame_1(2,3)   &
                         +Qrot(2,3)*X_frame_1(3,3))
         M_2(3,1) = 2.D0*(Qrot(3,1)*X_frame_1(1,1)   &
                         +Qrot(3,2)*X_frame_1(2,1)   &
                         +Qrot(3,3)*X_frame_1(3,1))
         M_2(3,2) = 2.D0*(Qrot(3,1)*X_frame_1(1,2)   &
                         +Qrot(3,2)*X_frame_1(2,2)   &
                         +Qrot(3,3)*X_frame_1(3,2))
         M_2(3,3) = 2.D0*(Qrot(3,1)*X_frame_1(1,3)   &
                         +Qrot(3,2)*X_frame_1(2,3)   &
                         +Qrot(3,3)*X_frame_1(3,3))
         X_frame_2(1,1) = 2.D0*(M_2(1,1)*Qirot(1,1)   &
                               +M_2(1,2)*Qirot(2,1)   &
                               +M_2(1,3)*Qirot(3,1))
         X_frame_2(1,2) = 2.D0*(M_2(1,1)*Qirot(1,2)   &
                               +M_2(1,2)*Qirot(2,2)   &
                               +M_2(1,3)*Qirot(3,2))
         X_frame_2(1,3) = 2.D0*(M_2(1,1)*Qirot(1,3)   &
                               +M_2(1,2)*Qirot(2,3)   &
                               +M_2(1,3)*Qirot(3,3))
         X_frame_2(2,1) = 2.D0*(M_2(2,1)*Qirot(1,1)   &
                               +M_2(2,2)*Qirot(2,1)   &
                               +M_2(2,3)*Qirot(3,1))
         X_frame_2(2,2) = 2.D0*(M_2(2,1)*Qirot(1,2)   &
                               +M_2(2,2)*Qirot(2,2)   &
                               +M_2(2,3)*Qirot(3,2))
         X_frame_2(2,3) = 2.D0*(M_2(2,1)*Qirot(1,3)   &
                               +M_2(2,2)*Qirot(2,3)   &
                               +M_2(2,3)*Qirot(3,3))
         X_frame_2(3,1) = 2.D0*(M_2(3,1)*Qirot(1,1)   &
                               +M_2(3,2)*Qirot(2,1)   &
                               +M_2(3,3)*Qirot(3,1))
         X_frame_2(3,2) = 2.D0*(M_2(3,1)*Qirot(1,2)   &
                               +M_2(3,2)*Qirot(2,2)   &
                               +M_2(3,3)*Qirot(3,2))
         X_frame_2(3,3) = 2.D0*(M_2(3,1)*Qirot(1,3)   &
                               +M_2(3,2)*Qirot(2,3)   &
                               +M_2(3,3)*Qirot(3,3))
      ELSEIF ( Ixg.EQ.2 ) THEN
! Rotate the local tensor "X_frame_2" from the local frame (2) to the
!         global frame (1)
!
! We use the inverse of the rotation matrix QRot, which is
!         its transpose
!
! Rotate the vector "X_frame_2" from the local frame (2) to the
!         global frame (1)
         M_1(1,1) = 2.D0*(Qirot(1,1)*X_frame_2(1,1)   &
                         +Qirot(1,2)*X_frame_2(2,1)   &
                         +Qirot(1,3)*X_frame_2(3,1))
         M_1(1,2) = 2.D0*(Qirot(1,1)*X_frame_2(1,2)   &
                         +Qirot(1,2)*X_frame_2(2,2)   &
                         +Qirot(1,3)*X_frame_2(3,2))
         M_1(1,3) = 2.D0*(Qirot(1,1)*X_frame_2(1,3)   &
                         +Qirot(1,2)*X_frame_2(2,3)   &
                         +Qirot(1,3)*X_frame_2(3,3))
         M_1(2,1) = 2.D0*(Qirot(2,1)*X_frame_2(1,1)   &
                         +Qirot(2,2)*X_frame_2(2,1)   &
                         +Qirot(2,3)*X_frame_2(3,1))
         M_1(2,2) = 2.D0*(Qirot(2,1)*X_frame_2(1,2)   &
                         +Qirot(2,2)*X_frame_2(2,2)   &
                         +Qirot(2,3)*X_frame_2(3,2))
         M_1(2,3) = 2.D0*(Qirot(2,1)*X_frame_2(1,3)   &
                         +Qirot(2,2)*X_frame_2(2,3)   &
                         +Qirot(2,3)*X_frame_2(3,3))
         M_1(3,1) = 2.D0*(Qirot(3,1)*X_frame_2(1,1)   &
                         +Qirot(3,2)*X_frame_2(2,1)   &
                         +Qirot(3,3)*X_frame_2(3,1))
         M_1(3,2) = 2.D0*(Qirot(3,1)*X_frame_2(1,2)   &
                         +Qirot(3,2)*X_frame_2(2,2)   &
                         +Qirot(3,3)*X_frame_2(3,2))
         M_1(3,3) = 2.D0*(Qirot(3,1)*X_frame_2(1,3)   &
                         +Qirot(3,2)*X_frame_2(2,3)   &
                         +Qirot(3,3)*X_frame_2(3,3))

         X_frame_1(1,1) = 2.D0*(M_1(1,1)*Qrot(1,1)   &
                               +M_1(1,2)*Qrot(2,1)   &
                               +M_1(1,3)*Qrot(3,1))
         X_frame_1(1,2) = 2.D0*(M_1(1,1)*Qrot(1,2)   &
                               +M_1(1,2)*Qrot(2,2)   &
                               +M_1(1,3)*Qrot(3,2))
         X_frame_1(1,3) = 2.D0*(M_1(1,1)*Qrot(1,3)   &
                               +M_1(1,2)*Qrot(2,3)   &
                               +M_1(1,3)*Qrot(3,3))
         X_frame_1(2,1) = 2.D0*(M_1(2,1)*Qrot(1,1)   &
                               +M_1(2,2)*Qrot(2,1)   &
                               +M_1(2,3)*Qrot(3,1))
         X_frame_1(2,2) = 2.D0*(M_1(2,1)*Qrot(1,2)   &
                               +M_1(2,2)*Qrot(2,2)   &
                               +M_1(2,3)*Qrot(3,2))
         X_frame_1(2,3) = 2.D0*(M_1(2,1)*Qrot(1,3)   &
                               +M_1(2,2)*Qrot(2,3)   &
                               +M_1(2,3)*Qrot(3,3))
         X_frame_1(3,1) = 2.D0*(M_1(3,1)*Qrot(1,1)   &
                               +M_1(3,2)*Qrot(2,1)   &
                               +M_1(3,3)*Qrot(3,1))
         X_frame_1(3,2) = 2.D0*(M_1(3,1)*Qrot(1,2)   &
                               +M_1(3,2)*Qrot(2,2)   &
                               +M_1(3,3)*Qrot(3,2))
         X_frame_1(3,3) = 2.D0*(M_1(3,1)*Qrot(1,3)   &
                               +M_1(3,2)*Qrot(2,3)   &
                               +M_1(3,3)*Qrot(3,3))
      ENDIF
      END SUBROUTINE QMATRIXROTATETENSOR

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !                                                                            !
!  Subroutine name: QINROTATE                                                 !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Compute the new quaternion of the contact frame:                  !
!          Qc_new -> Qc + dQc = Qc + (1/2) dw_global*Qc,                      !
!   or     Qc_new -> Qc + dQc = Qc + (1/2) Qc*dw_local                        !
!           where "*" is the quaternion composition operator.                 !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE QINCROTATE(Qc,Qc_new,Dw)
      IMPLICIT NONE

      INTEGER*4 :: Inorm
! Dw (=w*dtsolid)is the angular velocity incremental in the global frame
      DOUBLE PRECISION :: Dw, fqc, Qc, Qc_new, wq21, wq31, wq41,  &
                     & wq22, wq32, wq42, wq23, wq33, wq43, wq24, &
                     & wq34, wq44
      DIMENSION :: Qc(4), Qc_new(4), Dw(3)
      wq21 = Dw(1)*Qc(1)
      wq31 = Dw(2)*Qc(1)
      wq41 = Dw(3)*Qc(1)
      wq22 = Dw(1)*Qc(2)
      wq32 = Dw(2)*Qc(2)
      wq42 = Dw(3)*Qc(2)
      wq23 = Dw(1)*Qc(3)
      wq33 = Dw(2)*Qc(3)
      wq43 = Dw(3)*Qc(3)
      wq24 = Dw(1)*Qc(4)
      wq34 = Dw(2)*Qc(4)
      wq44 = Dw(3)*Qc(4)
!
      Qc_new(1) = Qc(1) + 0.5D0*(-wq22-wq33-wq44)
      Qc_new(2) = Qc(2) + 0.5D0*(wq21-wq43+wq34)
      Qc_new(3) = Qc(3) + 0.5D0*(wq31+wq42-wq24)
      Qc_new(4) = Qc(4) + 0.5D0*(wq41-wq32+wq23)
! Inorm = 2, default
     Inorm = 2
! Normalize the quaternion Qc_new, but only when "inorm" is not 0
      IF ( Inorm.NE.0 ) THEN
         IF ( Inorm.EQ.1 ) THEN
! The Katz correction factor
         fqc = 1.5D0 - 0.5D0*(Qc_new(1)**2+Qc_new(2)**2+Qc_new(3)**2+    &
                Qc_new(4)**2)
         ELSEIF ( Inorm.EQ.2 ) THEN
! A precise adjustment is, of course,
          fqc = 1.D0/SQRT(Qc_new(1)**2+Qc_new(2)**2+Qc_new(3)**2+Qc_new(4)**2)
         ENDIF
! Now apply the correction factor
         Qc_new(1) = fqc*Qc_new(1)
         Qc_new(2) = fqc*Qc_new(2)
         Qc_new(3) = fqc*Qc_new(3)
         Qc_new(4) = fqc*Qc_new(4)
      ENDIF

      END SUBROUTINE QINCROTATE
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !                                                                            !
!  Subroutine name: QINROTATE                                                 !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Compute the new quaternion of the contact frame:                  !
!          Qc_new -> Qc + dQc = Qc + (1/2) dw_global*Qc,                      !
!   or     Qc_new -> Qc + dQc = Qc + (1/2) Qc*dw_local  (this method is used) !
!           where "*" is the quaternion composition operator.                 !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE QINCROTATE2(Qc,Qc_new,Dw)
      IMPLICIT NONE
      INTEGER*4 :: Inorm
! Dw (=w*dtsolid)is the angular velocity incremental in the global frame
      DOUBLE PRECISION :: Dw, fqc, Qc, Qc_new, wq21, wq31, wq41,  &
                     & wq22, wq32, wq42, wq23, wq33, wq43, wq24, &
                     & wq34, wq44
      DIMENSION :: Qc(4), Qc_new(4), Dw(3)

      wq21 = Dw(1)*Qc(1)
      wq31 = Dw(2)*Qc(1)
      wq41 = Dw(3)*Qc(1)
      wq22 = Dw(1)*Qc(2)
      wq32 = Dw(2)*Qc(2)
      wq42 = Dw(3)*Qc(2)
      wq23 = Dw(1)*Qc(3)
      wq33 = Dw(2)*Qc(3)
      wq43 = Dw(3)*Qc(3)
      wq24 = Dw(1)*Qc(4)
      wq34 = Dw(2)*Qc(4)
      wq44 = Dw(3)*Qc(4)
!
      Qc_new(1) = Qc(1) + 0.5D0*(-wq22-wq33-wq44)
      Qc_new(2) = Qc(2) + 0.5D0*(wq21+wq43-wq34)
      Qc_new(3) = Qc(3) + 0.5D0*(wq31+wq24-wq42)
      Qc_new(4) = Qc(4) + 0.5D0*(wq41+wq32-wq23)
! Inorm = 2, default
     Inorm = 2
! Normalize the quaternion Qc_new, but only when "inorm" is not 0
      IF ( Inorm.NE.0 ) THEN
         IF ( Inorm.EQ.1 ) THEN
! The Katz correction factor
         fqc = 1.5D0 - 0.5D0*(Qc_new(1)**2+Qc_new(2)**2+Qc_new(3)**2+    &
                Qc_new(4)**2)
         ELSEIF ( Inorm.EQ.2 ) THEN
! A precise adjustment is, of course,
          fqc = 1.D0/SQRT(Qc_new(1)**2+Qc_new(2)**2+Qc_new(3)**2+Qc_new(4)**2)
         ENDIF
! Now apply the correction factor
         Qc_new(1) = fqc*Qc_new(1)
         Qc_new(2) = fqc*Qc_new(2)
         Qc_new(3) = fqc*Qc_new(3)
         Qc_new(4) = fqc*Qc_new(4)
      ENDIF

      END SUBROUTINE QINCROTATE2
END MODULE SQ_ROTATION_MOD
