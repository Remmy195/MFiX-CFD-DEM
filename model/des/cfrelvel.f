MODULE CFRELVEL_MOD
CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFRELVEL
!  Purpose: Calculate the normal and tangential components of the
!           relative velocity between contacting particles
!
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04
!  Reviewer: Rahul Garg                               Date: 01-Aug-07
!
!  Comments: Relative (translational) velocity required for eqn 6
!  from the following paper:
!    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical
!    simulation of plug glow of cohesionless particles in a
!    horizontal pipe", Powder technology, 71, 239-250, 1992
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   SUBROUTINE CFRELVEL(L, II, VRN, VSLIP, NORM, DIST_LI)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement, only: DES_VEL_NEW, DES_RADIUS, OMEGA_NEW, CROSS
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! indices of particle-particle contact pair
      INTEGER, INTENT(IN) :: L, II
! distance between particle centers
      DOUBLE PRECISION, INTENT(IN) :: DIST_LI
! unit normal vector along the line of contact pointing from
! particle L to particle II
      DOUBLE PRECISION, INTENT(IN) :: NORM(3)
! slip velocity at point of contact
      DOUBLE PRECISION, INTENT(OUT) :: VSLIP(3)
! normal component of relative contact velocity (scalar)
      DOUBLE PRECISION, INTENT(OUT) :: VRN
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! translational relative velocity
      DOUBLE PRECISION :: VRELTRANS(3)
! rotational velocity at point of contact
      DOUBLE PRECISION :: V_ROT(3), OMEGA_SUM(3)
! distance from the contact point to the particle centers
      DOUBLE PRECISION :: DIST_CL, DIST_CI
!-----------------------------------------------

! translational relative velocity
         VRELTRANS(:) = (DES_VEL_NEW(L,:) - DES_VEL_NEW(II,:))

! calculate the distance from the particle center to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
         DIST_CL = (DIST_LI**2 + DES_RADIUS(L)**2 - DES_RADIUS(II)**2)/&
            (2.d0*DIST_LI)
         DIST_CI = DIST_LI - DIST_CL

         OMEGA_SUM(:) = OMEGA_NEW(L,:)*DIST_CL + &
              OMEGA_NEW(II,:)*DIST_CI

! calculate the rotational relative velocity
      V_ROT = CROSS(OMEGA_SUM, NORM)

! total relative velocity
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

! normal component of relative velocity (scalar)
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

! slip velocity of the contact point
! Equation (8) in Tsuji et al. 1992
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)

      RETURN

   END SUBROUTINE CFRELVEL



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: Sq_CFRELVEL
!  Purpose: Calculate the normal and tangential components of the
!           relative velocity between contacting superquadric particles
!
!  Author: Xi Gao                                     Date: 25-Feb-2019
!
!  Purpose: Calculate the normal and tangential components of the
!           relative velocity between two superquadric particles
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   SUBROUTINE Sq_CFRELVEL(L, II,POS_L, POS_II,X, VRN, VSLIP, NORM,dist_li,shape_L,shape_II)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement, only: DES_VEL_NEW, DES_RADIUS, OMEGA_NEW, CROSS
	  USE discretelement, only: SuperDEM
      use usr
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! indices of particle-particle contact pair
      INTEGER, INTENT(IN) :: L, II
! distance between particle centers
      DOUBLE PRECISION, INTENT(IN) :: DIST_LI
! unit normal vector along the line of contact pointing from
! particle L to particle II
      DOUBLE PRECISION, INTENT(IN) :: NORM(3)
! slip velocity at point of contact
      DOUBLE PRECISION, INTENT(OUT) :: VSLIP(3)
! normal component of relative contact velocity (scalar)
      DOUBLE PRECISION, INTENT(OUT) :: VRN
! particle psotion
      DOUBLE PRECISION :: POS_L(3),POS_II(3),X(6)
! momentum arms
      DOUBLE PRECISION :: MOM_L(3),MOM_II(3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! translational relative velocity
      DOUBLE PRECISION :: VRELTRANS(3)
! rotational velocity at point of contact
      DOUBLE PRECISION :: V_ROT(3), OMEGA_SUM(3)
! distance from the contact point to the particle centers
      DOUBLE PRECISION :: DIST_CL, DIST_CI

      INTEGER :: shape_L, shape_II
      LOGICAL,PARAMETER :: ldebug = .false.
!-----------------------------------------------

! translational relative velocity
      VRELTRANS(:) = (DES_VEL_NEW(L,:) - DES_VEL_NEW(II,:))
         if(ldebug) then
             WRITE(*,'(/," V_reltrans(:)")')
             WRITE(*,'(3(f8.5))') vreltrans(:)
         endif
      If (SuperDEM .and. (shape_L .gt. 0 .or. shape_II .gt. 0)) then
! Momentum of L particle  (center to contact point)
        MOM_L(1)=X(1)-POS_L(1)
        MOM_L(2)=X(2)-POS_L(2)
        MOM_L(3)=X(3)-POS_L(3)

! Momentum arm of II particle  (center to contact point)
        MOM_II(1)=X(4)-POS_II(1)
        MOM_II(2)=X(5)-POS_II(2)
        MOM_II(3)=X(6)-POS_II(3)
! calculate the rotational relative velocity
        V_ROT(:)= CROSS(OMEGA_NEW(L,:), MOM_L(:))- &
                  CROSS(OMEGA_NEW(II,:), MOM_II(:))
       if(ldebug) then
         WRITE(*,'(/," MOM_L(:),MOM_II(:),V_ROT(:)")')
         WRITE(*,'(6(f8.5))') MOM_L(:),MOM_II(:)
         WRITE(*,'(3(f8.5))') V_ROT(:)
       endif
      else
! calculate the distance from the particle center to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
         DIST_CL = (DIST_LI**2 + DES_RADIUS(L)**2 - DES_RADIUS(II)**2)/&
            (2.d0*DIST_LI)
         DIST_CI = DIST_LI - DIST_CL

         OMEGA_SUM(:) = OMEGA_NEW(L,:)*DIST_CL + &
              OMEGA_NEW(II,:)*DIST_CI

! calculate the rotational relative velocity
         V_ROT = CROSS(OMEGA_SUM, NORM)

     ENDIF

! total relative velocity
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

! normal component of relative velocity (scalar)
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

! slip velocity of the contact point
! Equation (8) in Tsuji et al. 1992
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)
     if(ldebug) then
       WRITE(*,'(/," VRELTRANS(:),VRN,VSLIP(:)")')
       WRITE(*,'(3(f8.5))') VRELTRANS(:)
       WRITE(*,'(4(f8.5))') VRN,VSLIP(:)
     endif
      RETURN

   END SUBROUTINE Sq_CFRELVEL

END MODULE CFRELVEL_MOD
