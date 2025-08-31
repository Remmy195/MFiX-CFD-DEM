#include "error.inc"
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: sq_math_mod                                            !
!  Author: Xi Gao                                 Date: 25-Feb-2019    !
!                                                                      !
!  Purpose: Mathematical calculation                                   !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     MODULE Sq_math_mod
     use error_manager
     CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: matinv2(A,B)                                       !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:  Performs a direct calculation of the inverse              !
!	 of a 2x2 matrix                                                  !
!                                                                      !
!  References:                                                         !
!  http://fortranwiki.org/fortran/show/Matrix+inversion                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     Subroutine matinv2(A,B)
     IMPLICIT none
     DOUBLE PRECISION, intent(in) :: A(2,2)   ! Matrix
     DOUBLE PRECISION             :: B(2,2)   ! Inverse matrix
     DOUBLE PRECISION             :: detinv
! Calculate the inverse determinant of the matrix
     detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
! Calculate the inverse of the matrix
     B(1,1) = +detinv * A(2,2)
     B(2,1) = -detinv * A(2,1)
     B(1,2) = -detinv * A(1,2)
     B(2,2) = +detinv * A(1,1)
     END Subroutine matinv2

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: matdet2(A,det)                                     !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:  calculate the determinant of the matrix                   !
!                                                                      !
!  References:                                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     Subroutine matdet2(A,det)
     IMPLICIT none
     DOUBLE PRECISION, intent(in) :: A(2,2)   !! Matrix
     DOUBLE PRECISION             :: det

! Calculate the determinant of the matrix
     det = A(1,1)*A(2,2) - A(1,2)*A(2,1)

     END Subroutine matdet2

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: matinv3(A,B)                                       !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:  Performs a direct calculation of the inverse              !
!	 of a 3x3 matrix                                                  !
!                                                                      !
!  References:                                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     Subroutine matinv3(A,B)
     IMPLICIT none
     DOUBLE PRECISION, intent(in) :: A(3,3)   !! Matrix
     DOUBLE PRECISION             :: B(3,3)   !! Inverse matrix
     DOUBLE PRECISION             :: detinv
! Calculate the inverse determinant of the matrix
     detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

! Calculate the inverse of the matrix
     B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
     B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
     B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
     B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
     B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
     B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
     B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
     B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
     B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
     END Subroutine matinv3

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: matdet3(A,det)                                     !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:  calculate the determinant of the matrix                   !
!                                                                      !
!  References:                                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     Subroutine matdet3(A,det)
     IMPLICIT none
     DOUBLE PRECISION, intent(in) :: A(3,3)   !! Matrix
     DOUBLE PRECISION             :: det

! Calculate the determinant of the matrix
     det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)
     END Subroutine matdet3

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: matinv4(A,B)                                       !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:  Performs a direct calculation of the inverse              !
!	 of a 4x4 matrix                                               !
!                                                                      !
!  References:                                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     Subroutine matinv4(A,B)
     IMPLICIT none
     DOUBLE PRECISION, intent(in) :: A(4,4)   !! Matrix
     DOUBLE PRECISION             :: B(4,4)   !! Inverse matrix
     DOUBLE PRECISION             :: detinv

! Calculate the inverse determinant of the matrix
     detinv = &
      1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-&
                        A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
       - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-&
                        A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
       + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-&
                        A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
       - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-&
                        A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

! Calculate the inverse of the matrix
     B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+&
                      A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+&
                      A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
     B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+&
                      A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+&
                      A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
     B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+&
                      A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+&
                      A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
     B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+&
                      A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+&
                      A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
     B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+&
                      A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+&
                      A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
     B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+&
                      A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+&
                      A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
     B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+&
                      A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+&
                      A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
     B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+&
                      A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+&
                      A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
     B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+&
                      A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+&
                      A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
     B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+&
                      A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+&
                      A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
     B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+&
                      A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+&
                      A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
     B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+&
                      A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+&
                      A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
     B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+&
                      A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+&
                      A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
     B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+&
                      A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+&
                      A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
     B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+&
                      A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+&
                      A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
     B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+&
                      A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+&
                      A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
     END Subroutine matinv4

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: matdet4(A,det)                                     !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:  calculate the determinant of the matrix                   !
!                                                                      !
!  References:                                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     Subroutine matdet4(A,det)
     IMPLICIT none
     DOUBLE PRECISION, intent(in) :: A(4,4)   !! Matrix
     DOUBLE PRECISION             :: det

! Calculate the inverse determinant of the matrix
     det = &
      A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-&
                        A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
       - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-&
                        A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
       + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-&
                        A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
       - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-&
                        A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))

     END Subroutine matdet4

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: matvecN(A,vec,N)                                   !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:  calculate the matrix-vector multiplication                !
!                                                                      !
!  References:                                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     Subroutine matvec4(A,B,C)
     IMPLICIT none
     INTEGER                      :: i,j
     DOUBLE PRECISION, intent(in) :: A(4,4)   !! Matrix
     DOUBLE PRECISION             :: B(4),C(4)

     DO i=1,4
        C(i)=0.0d0
        DO j=1,4
!              if (A(i,j) .ne. 0.0d0 .and. B(j) .ne. 0.0d0)
                 C(i) = C(i) + A(i,j)*B(j)
!              endif
       ENDDO
     ENDDO

     END Subroutine matvec4

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: dotvecN(A,B,C,N)                                   !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:  calculate the vector-vector dot                           !
!                                                                      !
!  References:                                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     Subroutine dotvecN(A,B,C,N)
     IMPLICIT none
     INTEGER                      :: n,i,j
     DOUBLE PRECISION             :: A(n),B(n),C
     C=0.0d0
     DO i=1,n
        C= C + A(i)*B(i)
     ENDDO

     END Subroutine dotvecN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: normvecNsq(A,N,NORMsq,NORM)                        !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:                                                            !
!                                                                      !
!  References:                                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     Subroutine normvecN(A,N,NORMsq,NORM)
     IMPLICIT none
     INTEGER                      :: n,i,j
     DOUBLE PRECISION             :: A(n),NORMsq,NORM

     NORMsq=0.0d0
     DO i=1,n
        NORMsq = NORMsq + A(i)*A(i)
     ENDDO

     NORM=DSQRT(NORMsq)

     END Subroutine normvecN

     Subroutine crossvec3(A,B,C)
     IMPLICIT none
     INTEGER                      :: n,i,j
     DOUBLE PRECISION             :: A(3),B(3),C(3)

     C(1)= A(2)*B(3)-A(3)*B(2)
     C(2)= A(3)*B(1)-A(1)*B(3)
     C(3)= A(1)*B(2)-A(2)*B(1)

     END Subroutine crossvec3

     END MODULE Sq_math_mod
