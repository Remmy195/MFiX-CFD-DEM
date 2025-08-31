#include "error.inc"
     Module gmres_general
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: gmres_general                                          !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:  solve linear equations using general minimized            !
!	 residue method                                                !
!                                                                      !
!  References:                                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: gmres(A,vec,N)                                     !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:  solve linear equations using general minimized            !
!	  residue method                                               !
!                                                                      !
!  References:                                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine gmres( A, bb, xx0, n,m, threshold, xx,c)
      Implicit none
      INTENT(out) c
      Integer :: n,m,i,J,k

      Double Precision :: A(n,n),b(n,1),x0(n,1),c(n,1),x(n,1),r(n,1),threshold
      Double Precision :: normbsq,normb,normrsq,normr,normqq
      Double Precision :: sn(m+1,1),cs(m+1,1),e1(n,1),e(m+1,1)
      Double Precision :: error,r_norm
      Double Precision :: Q(n,m+1),beta(m+1,1),H(m+1,m),hh(m+1,1)
      Double Precision :: inv_h(m,m+1)
      Double Precision :: cs_k,sn_k,hhh
      Double Precision :: qq(n,1),qq_inv(1,n)
      Double Precision :: y(m+1,1)
      Double Precision :: dd
      Double Precision :: bb(n),xx0(n),xx(n)
!
      do i=1,n
         b(i,1)=bb(i)
         x0(i,1)=xx0(i)
      enddo

      call matvecN(A(1:n,1:n),x0(1:n,1),C(1:n,1),n)

      r(1:n,1)=b(1:n,1)-c(1:n,1)

      normb=0.0d0
      normr=0.0d0
      Do i=1,n
          normb=normb+b(i,1)*b(i,1)
          normR=normR+r(i,1)*r(i,1)
      ENDDO
      normb=DSQRT(normb)
      normr=DSQRT(normr)
      error= normr/normb

!initialize the 1D vectors
      xx(1:n) = 0.0d0
      sn(1:m,1) = 0.0d0
      cs(1:m,1) = 0.0d0
      e1(1:n,1) = 0.0d0
      e1(1,1) = 1.0d0

      e(1,1)=error
      r_norm=normr

      Q(1:n,1) = r(1:n,1)/r_norm
      beta(1:n,1) = r_norm*e1(1:n,1);

      DO k = 1,m
        call matvecN(A(1:n,1:n),Q(1:n,k),qq(1:n,1),n)
        qq_inv(1,1:n)=0.0d0
        qq_inv(1,1:n)=qq(1:n,1)

        DO i = 1,k
              dd=0.0d0
              DO j=1,n
                 dd= dd + qq_inv(1,j)*Q(j,i)
              ENDDO

           hh(i,1)=dd
           qq(1:n,1)=qq(1:n,1)-hh(i,1)*Q(1:n,i)
        enddo
          normqq=0.0d0
            Do i=1,n
               normqq=normqq+qq(i,1)*qq(i,1)
            ENDDO
               normqq=DSQRT(normqq)

        hh(k+1,1) = normqq
        if (normqq.ne.0.0) then
            qq(1:n,1) = qq(1:n,1) / normqq
        endif ! what to do if normqq is zero?

        Q(1:n,k+1) =qq(1:n,1)
        H(1:k,k)=hh(1:k,1)
        H(k+1,k)=hh(k+1,1)
!eliminate the last element in H ith row and update the rotation matrix
       call apply_givens_rotation(k,m,H(1:k+1,k), cs, sn,cs_k,sn_k);
       cs(k,1)=cs_k
       sn(k,1)=sn_k
!update the residual vector
       beta(k+1,1) = -sn(k,1)*beta(k,1)
       beta(k,1)   = cs(k,1)*beta(k,1)
       error  = dabs(beta(k+1,1)) / normb

!save the error
       e(k+1,1)=error

       if ( error .lt. threshold) then
           goto 100
       endif
       enddo
       k = m

!calculate the result
100    continue
       call inverse_uptri_matrix(h(1:k,1:k),k,inv_h(1:k,1:k))
       call matvecN(inv_h(1:k,1:k),beta(1:k,1),y(1:k,1),k)
       x(:,1)=x(:,1)+matmul(Q(1:n,1:k),y(1:k,1))

      do i=1,n
         xx(i)=x(i,1)
      enddo
      c(1:n,1)=0.0d0
      call matvecN(A(1:n,1:n),x(1:n,1),c(1:n,1),n)
      return
      end subroutine gmres

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: apply_givens_rotation                              !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:  Applying Givens Rotation to H col                         !
!                                                                      !
!  References:                                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     Subroutine apply_givens_rotation(k,m,h, cs, sn,cs_k, sn_k)
     Implicit none
     Integer :: i,k,m
     Double Precision :: cs(m,1),sn(m,1),h(k+1,1),temp
     Double Precision :: cs_k,sn_k
!apply for ith column
     Do i = 1,k-1
         temp   =  cs(i,1)*h(i,1) + sn(i,1)*h(i+1,1)
         h(i+1,1) = -sn(i,1)*h(i,1) + cs(i,1)*h(i+1,1)
         h(i,1)   = temp
     enddo
!update the next sin cos values for rotation
     call givens_rotation(h(k,1), h(k+1,1),cs_k,sn_k)
!eliminate H(i+1,i)
     h(k,1) = cs_k*h(k,1) + sn_k*h(k+1,1)
     h(k+1,1) = 0.0d0
     end subroutine apply_givens_rotation

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: givens_rotation(V1,V2,CS,SN)                       !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:  Calculate the Given rotation matrix                       !
!                                                                      !
!  References:                                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     Subroutine givens_rotation(v1, v2,cs,sn)
     Implicit none
     Double Precision :: cs,sn,v1,v2,t

     if (v1 .eq. 0d0) then
          cs = 0.0d0
          sn = 1.0d0
     else
          t=sqrt(v1**2.0d0+v2**2.0d0)
          cs = dabs(v1) / t
          sn = cs * v2 / v1
     endif
     end subroutine givens_rotation

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: inverse_uptri_matrix                               !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:                                                            !
!                                                                      !
!  References:                                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     subroutine inverse_uptri_matrix(A,n,inverse_A)
     implicit none
     integer :: i,j,k,n
     double precision :: summ
     double precision :: A(n,n),inverse_A(n,n)
     do i=1,n
        if (A(i,i)>= 0.0d0 .and. dabs(A(i,i))<1e-30) then
            inverse_A(i,i)=1e30
        elseif (A(i,i)<= 0.0d0 .and. dabs(A(i,i))<1e-30) then
            inverse_A(i,i)=-1e30
        else
        inverse_A(i,i)=1/A(i,i)
        endif
     enddo
     do i=1,n
        do j=1,n
           if(i<j) then
           summ=0.0d0
           do k=i,j-1
           summ=summ+inverse_A(i,k)*A(k,j)
           enddo
           if (A(j,j)<=0.0d0 .and. dabs(A(j,j))<1e-30 ) then

           inverse_A(j,j)=-summ/(-1e-30)
           elseif (A(j,j)>=0.0d0 .and. dabs(A(j,j))<1e-30 ) then

           inverse_A(j,j)=-summ/1e-30
           else
           inverse_A(i,j)=-summ/A(j,j)
           endif
           endif
        enddo
     enddo
     return
     end subroutine inverse_uptri_matrix

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: matvecN(A,vec,N)                                   !
!  Author: Xi Gao                                  Date: 25-Feb-2019   !
!  Purpose:  calculate the matrix-vector multiplication                !
!                                                                      !
!  References:                                                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     Subroutine matvecN(A,B,C,N)
     IMPLICIT none
     INTEGER                      :: n,i,j
     DOUBLE PRECISION, intent(in) :: A(n,n)
     DOUBLE PRECISION             :: B(n,1),C(n,1)

     DO i=1,n
        C(i,1)=0.0d0
        DO j=1,n
!              if (A(i,j) .ne. 0.0d0 .and. B(j) .ne. 0.0d0)
                 C(i,1) = C(i,1) + A(i,j)*B(j,1)
!              endif
       ENDDO
     ENDDO

     END Subroutine matvecN

     END MODULE GMRES_general
