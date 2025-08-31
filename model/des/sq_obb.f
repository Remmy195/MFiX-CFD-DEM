#include "error.inc"

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Module name: SQ_OBB_MOD                                          !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Contact detection between two Orentied Bounding Box (OBB-OBB)     !
!          broad phase contact before enter narrow phase S-S contact detection!
!                                                                             !
! Reference:                                                                  !
!   Xi Gao, Jia Yu, Ricardo JF Portal, Jean-François Dietiker,                !
!   Mehrdad Shahnam and William A Rogers,Development and validation           !
!   of SuperDEM for non-spherical particulate systems using a                 !
!   superquadric particle method, Particuology,                               !
!   2021: doi.org/10.1016/j.partic.2020.1011.1007.                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
MODULE SQ_OBB_MOD

CONTAINS

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Module name: SQ_OBB_CONTACT                                      !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Feed the coordinates of OBBs, Quaternions                         !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE SQ_OBB(p,axi,euler0A,euler1A,euler2A,euler3A, &
                                 euler0B,euler1B,euler2B,euler3B,&
                                 OBB_test)
     IMPLICIT none
     INTEGER,PARAMETER ::N=6
! Euler parameters of OBB, also called quaternions
     DOUBLE PRECISION :: euler0A,euler1A, euler2A, euler3A
     DOUBLE PRECISION :: euler0B,euler1B, euler2B, euler3B
! P is the location of OBB center, axi is the sem-axes
     DOUBLE PRECISION,DIMENSION (N) :: p, axi
! vector between locationS of two OBBs
     DOUBLE PRECISION :: pab(3),xab,yab,zab,HA(15),HB(15),HC(15)
! Transformation matrix
     DOUBLE PRECISION,DIMENSION (3,3) :: TRa(3,3),TRb(3,3),TRc(3,3)
     DOUBLE PRECISION :: OBB_TEST
     INTEGER :: I

! Transformation Matrix A
! 1st column = TRA
     TRa(1,1)= euler0A**2 + euler1A**2 - euler2A**2 - euler3A**2
     TRa(2,1)= 2*(euler1A*euler2A + euler0A*euler3A)
     TRa(3,1)= 2*(euler1A*euler3A - euler0A*euler2A)
! 2nd column = TRA
     TRa(1,2)= 2*(euler1A*euler2A - euler0A*euler3A)
     TRa(2,2)= euler0A**2 - euler1A**2 + euler2A**2 - euler3A**2
     TRa(3,2)= 2*(euler2A*euler3A + euler0A*euler1A)
! 3rd column = TRA
     TRa(1,3)= 2*(euler1A*euler3A + euler0A*euler2A)
     TRa(2,3)= 2*(euler2A*euler3A - euler0A*euler1A)
     TRa(3,3)= euler0A**2 - euler1A**2 - euler2A**2 + euler3A**2

! Transformation Matrix B
! 1st column = TRb
     TRb(1,1)= euler0B**2 + euler1B**2 - euler2B**2 - euler3B**2
     TRb(2,1)= 2*(euler1B*euler2B + euler0B*euler3B)
     TRb(3,1)= 2*(euler1B*euler3B - euler0B*euler2B)
! 2nd column = TRb
     TRb(1,2)= 2*(euler1B*euler2B - euler0B*euler3B)
     TRb(2,2)= euler0B**2 - euler1B**2 + euler2B**2 - euler3B**2
     TRb(3,2)= 2*(euler2B*euler3B + euler0B*euler1B)
! 3rd column = TRb
     TRb(1,3)= 2*(euler1B*euler3B + euler0B*euler2B)
     TRb(2,3)= 2*(euler2B*euler3B - euler0B*euler1B)
     TRb(3,3)= euler0B**2 - euler1B**2 - euler2B**2 + euler3B**2

! Transformation Matrix C=A^T*B
! 1st column = TRc
     TRc(1,1)=TRa(1,1)*TRb(1,1)+TRa(2,1)*TRb(2,1)+TRa(3,1)*TRb(3,1)
     TRc(1,2)=TRa(1,1)*TRb(1,2)+TRa(2,1)*TRb(2,2)+TRa(3,1)*TRb(3,2)
     TRc(1,3)=TRa(1,1)*TRb(1,3)+TRa(2,1)*TRb(2,3)+TRa(3,1)*TRb(3,3)
! 2nd column = TRc
     TRc(2,1)=TRa(1,2)*TRb(1,1)+TRa(2,2)*TRb(2,1)+TRa(3,2)*TRb(3,1)
     TRc(2,2)=TRa(1,2)*TRb(1,2)+TRa(2,2)*TRb(2,2)+TRa(3,2)*TRb(3,2)
     TRc(2,3)=TRa(1,2)*TRb(1,3)+TRa(2,2)*TRb(2,3)+TRa(3,2)*TRb(3,3)
! 3rd column = TRc
     TRc(3,1)=TRa(1,3)*TRb(1,1)+TRa(2,3)*TRb(2,1)+TRa(3,3)*TRb(3,1)
     TRc(3,2)=TRa(1,3)*TRb(1,2)+TRa(2,3)*TRb(2,2)+TRa(3,3)*TRb(3,2)
     TRc(3,3)=TRa(1,3)*TRb(1,3)+TRa(2,3)*TRb(2,3)+TRa(3,3)*TRb(3,3)

! vector between two OBBS in the global coordinat
     pab(1)=p(4)-p(1)
     pab(2)=p(5)-p(2)
     pab(3)=p(6)-p(3)
! pab in the local coordinates of OBB A,=A^T*pab
     xab= TRa(1,1)*pab(1)+TRa(2,1)*pab(2)+TRa(3,1)*pab(3)
     yab= TRa(1,2)*pab(1)+TRa(2,2)*pab(2)+TRa(3,2)*pab(3)
     zab= TRa(1,3)*pab(1)+TRa(2,3)*pab(2)+TRa(3,3)*pab(3)

! 15 tests

! case 1:a1
     HA(1)=axi(1)
     HB(1)=axi(4)*dabs(TRc(1,1))+axi(5)*dabs(TRc(1,2))+axi(6)*dabs(TRc(1,3))
     HC(1)=dabs(xab)
! case 2:a2
     HA(2)=axi(2)
     HB(2)=axi(4)*dabs(TRc(2,1))+axi(5)*dabs(TRc(2,2))+axi(6)*dabs(TRc(2,3))
     HC(2)=dabs(yab)
! case 3:a3
     HA(3)=axi(3)
     HB(3)=axi(4)*dabs(TRc(3,1))+axi(5)*dabs(TRc(3,2))+axi(6)*dabs(TRc(3,3))
     HC(3)=dabs(zab)

! case 4:b1
     HA(4)=axi(1)*dabs(TRc(1,1))+axi(2)*dabs(TRc(2,1))+axi(3)*dabs(TRc(3,1))
     HB(4)=axi(4)
     HC(4)=dabs(TRb(1,1)*pab(1)+TRb(2,1)*pab(2)+TRb(3,1)*pab(3))
! case 5:b2
     HA(5)=axi(1)*dabs(TRc(1,2))+axi(2)*dabs(TRc(2,2))+axi(3)*dabs(TRc(3,2))
     HB(5)=axi(5)
     HC(5)=dabs(TRb(1,2)*pab(1)+TRb(2,2)*pab(2)+TRb(3,2)*pab(3))
! case 6:b3
     HA(6)=axi(1)*dabs(TRc(1,3))+axi(2)*dabs(TRc(2,3))+axi(3)*dabs(TRc(3,3))
     HB(6)=axi(6)
     HC(6)=dabs(TRb(1,3)*pab(1)+TRb(2,3)*pab(2)+TRb(3,3)*pab(3))
! case 7:a1xb1
     HA(7)=axi(2)*dabs(TRc(3,1))+axi(3)*dabs(TRc(2,1))
     HB(7)=axi(5)*dabs(TRc(1,3))+axi(6)*dabs(TRc(1,2))
     HC(7)=dabs(zab*TRc(2,1)-yab*TRc(3,1))
! case 8:a1xb2
     HA(8)=axi(2)*dabs(TRc(3,2))+axi(3)*dabs(TRc(2,2))
     HB(8)=axi(4)*dabs(TRc(1,3))+axi(6)*dabs(TRc(1,1))
     HC(8)=dabs(zab*TRc(2,2)-yab*TRc(3,2))
! case 9:a1xb3
     HA(9)=axi(2)*dabs(TRc(3,3))+axi(3)*dabs(TRc(2,3))
     HB(9)=axi(4)*dabs(TRc(1,2))+axi(5)*dabs(TRc(1,1))
     HC(9)=dabs(zab*TRc(2,3)-yab*TRc(3,3))

! case 10:a2xb1
     HA(10)=axi(1)*dabs(TRc(3,1))+axi(3)*dabs(TRc(1,1))
     HB(10)=axi(5)*dabs(TRc(2,3))+axi(6)*dabs(TRc(2,2))
     HC(10)=dabs(xab*TRc(3,1)-zab*TRc(1,1))
! case 11:a2xb2
     HA(11)=axi(1)*dabs(TRc(3,2))+axi(3)*dabs(TRc(1,2))
     HB(11)=axi(4)*dabs(TRc(2,3))+axi(6)*dabs(TRc(2,1))
     HC(11)=dabs(xab*TRc(3,2)-zab*TRc(1,2))
! case 12:a2xb3
     HA(12)=axi(1)*dabs(TRc(3,3))+axi(3)*dabs(TRc(1,3))
     HB(12)=axi(4)*dabs(TRc(2,2))+axi(5)*dabs(TRc(2,1))
     HC(12)=dabs(xab*TRc(3,3)-zab*TRc(1,3))

! case 13:a3xb1
     HA(13)=axi(1)*dabs(TRc(2,1))+axi(2)*dabs(TRc(1,1))
     HB(13)=axi(5)*dabs(TRc(3,3))+axi(6)*dabs(TRc(3,2))
     HC(13)=dabs(yab*TRc(1,1)-xab*TRc(2,1))
! case 14:a3xb2
     HA(14)=axi(1)*dabs(TRc(2,2))+axi(2)*dabs(TRc(1,2))
     HB(14)=axi(4)*dabs(TRc(3,3))+axi(6)*dabs(TRc(3,1))
     HC(14)=dabs(yab*TRc(1,2)-xab*TRc(2,2))
! case 15:a3xb3
     HA(15)=axi(1)*dabs(TRc(2,3))+axi(2)*dabs(TRc(1,3))
     HB(15)=axi(4)*dabs(TRc(3,2))+axi(5)*dabs(TRc(3,1))
     HC(15)=dabs(yab*TRc(1,3)-xab*TRc(2,3))

! test OBB contacts
     DO  I=1,15
! if there is a separating axie exists, then not contact
         IF(HC(I)-HA(I)-HB(I)>0.0) then

            OBB_test= 1.0d0

             goto 100
          ELSE
!  a separating axis can not be found, then contacted
            OBB_test= -1.0d0
!     WRITE(*,'(" Contact detected with OBBs: ", f8.4,/)') OBB_TEST
          ENDIF

     ENDDO
     100 continue
!     WRITE(*,'(" Contact detected with OBBs: ", f8.4,/)') OBB_TEST
     RETURN

     END subroutine sq_OBB

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Module name: SQ_OBB_CONTACT_TEST                                 !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Feed the coordinates of OBBs, Quaternions                         !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE SQ_OBB_TEST

     IMPLICIT none
     INTEGER,PARAMETER ::N=6

! Euler parameters of superquadric surfaces, also called quaternions
     DOUBLE PRECISION :: euler0A,euler1A, euler2A, euler3A
     DOUBLE PRECISION :: euler0B,euler1B, euler2B, euler3B
     DOUBLE PRECISION :: m1,n1,m2,n2
! P is the location of superquadric center, axi is the sem-axes
     DOUBLE PRECISION,DIMENSION (N) :: p, axi

! Transformation matrix
     DOUBLE PRECISION,DIMENSION (3,3) :: TRa(3,3),TRb(3,3)
     INTEGER :: I, NDAT
     DOUBLE PRECISION :: OBB_test
! input 6 OBB pairs properties for test
     OPEN(10,file= 'input_data0001.txt')
     OPEN(20,file= 'input_data0002.txt')
     OPEN(30,file= 'input_data0003.txt')
     OPEN(40,file= 'input_data0004.txt')
     OPEN(50,file= 'input_data0005.txt')
     OPEN(60,file= 'input_data0006.txt')
     WRITE(*,*) 'Theoretical contact condition:'
     WRITE(*,*) 'case1: distance    = 0.0002m'
     WRITE(*,*) 'case2: penetration =-0.00002m'
     WRITE(*,*) 'case3: distance    = 0.00012m'
     WRITE(*,*) 'case4: penetration =-0.0000381m'
     WRITE(*,*) 'case5: distance    = 0.000159m'
     WRITE(*,*) 'case6: distance    = 0.0001821m'

     do i=1, 6
        NDAT=10*i
        READ (NDAT,*) p(1), p(2), p(3)
        READ (NDAT,*) p(4), p(5), p(6)
        READ (NDAT,*)
        READ (NDAT,*) axi(1), axi(2), axi(3)
        READ (NDAT,*) axi(4), axi(5), axi(6)
        READ (NDAT,*)
        READ (NDAT,*) m1, n1
        READ (NDAT,*) m2, n2
        READ (NDAT,*)
        READ (NDAT,*) euler1A, euler2A, euler3A
        READ (NDAT,*) euler1B, euler2B, euler3B
        READ (NDAT,*)
        CLOSE(NDAT)

        euler0A=DSQRT(1.0-euler1A**2-euler2A**2-euler3A**2)
        euler0B=DSQRT(1.0-euler1B**2-euler2B**2-euler3B**2)


        CALL SQ_OBB(p,axi,euler0A,euler1A,euler2A,euler3A, &
                               euler0B,euler1B,euler2B,euler3B,&
                               OBB_test)

        IF (OBB_TEST .LT. 0.0) THEN
           WRITE(*,*) 'OBBs contact detected'
           WRITE(*,'(1(f8.3))')   OBB_test
        ELSE
           WRITE(*,*) 'OBBs contact NOT detected'
           WRITE(*,'(1(f8.3))')   OBB_test
        ENDIF
     ENDDO

     RETURN
     END subroutine sq_OBB_TEST

END MODULE SQ_OBB_MOD
