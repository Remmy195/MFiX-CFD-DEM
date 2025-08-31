
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Module name: SQ_satellites                                                 !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Calculate satellites points                                       !
!                                                                             !
!  Reference:                                                                 !
!     Xi Gao, Jia Yu, Liqiang Lu, Cheng Li and William A Rogers,              !
!     Development and validation of SuperDEM-CFD coupled model for simulating !
!     non-spherical particles hydrodynamics in fluidized beds,                !
!     Chemical Engineering Journal, 2021,420: 127654.                         !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE SQ_SATELLITES
      use SQ_PROPERTIES_MOD
      contains

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: Sq_satellites_local                                       !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: satellites in the local coordinate                                !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE Sq_satellites_local(nx,ny,nz,nt,axi_x,axi_y,axi_z,&
                                                 m,n,s_pos,s_weights)
      IMPLICIT NONE
! number of satalites in the x, y, z direaction of the OBB
      INTEGER :: iii,jjj,kkk,ttt
! number of satellites in the superquadric
      INTEGER :: s_ttt

      INTEGER :: Nx,Ny,Nz,Nt
      DOUBLE PRECISION :: axi(3),axi_x,axi_y,axi_z,m,n,d_x,d_y,d_z
      DOUBLE PRECISION :: pos_x, pos_y,pos_z,s_pos(Nt,3)
      DOUBLE PRECISION :: F, s_weights(Nt)

      s_ttt=0
      pos_x=0.0
      pos_y=0.0
      pos_z=0.0
      s_pos=0.0


      axi(1)=axi_x
      axi(2)=axi_y
      axi(3)=axi_z

      d_x=2.0*axi_x/nx
      d_y=2.0*axi_y/ny
      d_z=2.0*axi_z/nz
      ttt=0


      DO iii=0,nx-1
          pos_x = (iii+0.5)*d_x-axi_x;
          DO jjj=0,ny-1
              pos_y = (jjj+0.5)*d_y-axi_y;
              DO kkk=0,nz-1
                  pos_z = (kkk+0.5)*d_z-axi_z;
                  ttt=ttt+1
                  s_pos(ttt,1)=pos_x
                  s_pos(ttt,2)=pos_y
                  s_pos(ttt,3)=pos_z

              ENDDO
          ENDDO
      ENDDO

! test whether the satellites points are inside the superquadric
      DO ttt=1, Nt
          call SQ_INOUT_F_local(s_pos(ttt,:),axi,m,n,F)
          IF (f>= 0.0 ) THEN
              s_weights(ttt)=0.0
          else
! give a initial value for s_weight when the exact value is not known
              s_weights(ttt)=2.0
              s_ttt=s_ttt + 1
          ENDIF
      ENDDO
! calculate the real weight

      DO TTT=1,Nt
          IF (s_weights(ttt)>1.0) then
              s_weights(ttt)=1.0/s_ttt
       ENDIF
   ENDDO
END SUBROUTINE Sq_satellites_local

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: Sq_satellites_global                                      !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: satellites in the global coordinate                               !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     SUBROUTINE Sq_satellites_global(Nt, Qc,sat_pos, sat_pos_global)
     USE USR
     use SQ_ROTATION_MOD
     IMPLICIT NONE
     Integer :: nt, t,L
     DOUBLE PRECISION :: Qc(4),sat_pos, sat_pos_global,s_local,s_global
     Dimension:: sat_pos(nt,3),sat_pos_global(nt,3)
     Dimension:: s_local(3*nt),s_global(3*nt)

! satellites in the local coordinate
       DO t=1,nt
           s_local((t-1)*3+1)=sat_pos(t,1)
           s_local((t-1)*3+2)=sat_pos(t,2)
           s_local((t-1)*3+3)=sat_pos(t,3)
       ENDDO
! rotate the sat_pos_local to global
       call QMultiPointsROTATE(Qc,s_global,s_local,2,nt)

       DO t=1,nt
           sat_pos_global(t,1)=s_global((t-1)*3+1)
           sat_pos_global(t,2)=s_global((t-1)*3+2)
           sat_pos_global(t,3)=s_global((t-1)*3+3)
       ENDDO

       END SUBROUTINE Sq_satellites_global

      END MODULE SQ_SATELLITES
