!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: solids_pressure                                        C
!  Purpose: To compute solids pressure and its inverse                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-FEB-93  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Comments:                                                           C
!        see set_constants.f and constant_mod.f                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE solids_pressure

  USE constant
  USE param1


CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Neg_H = P_star                                                       !
!                                                                      !
! Returns the solids plastic-flow pressure at the given void fraction  !
! and maximum packing void fraction. Is the MFIX default.              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  DOUBLE PRECISION FUNCTION Neg_H(epg,epg_star)
    IMPLICIT NONE
! void fraction
    DOUBLE PRECISION :: epg
! void fraction at maximum packing
    DOUBLE PRECISION :: epg_star
!----------------------------------------------------------------------
! P_star = a_ps*(EP_star - EP_g)**b_ps
    Neg_H = a_ps * (MAX(ZERO, (epg_star - epg )))**b_ps

  END FUNCTION Neg_H

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! inverse of Neg_H = ep_g                                              !
!                                                                      !
! Returns the void fraction corresponding to the given solids          !
! plastic-flow pressure and maximum packing void fraction              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  DOUBLE PRECISION FUNCTION INV_H(ps_star,epg_star)
    IMPLICIT NONE
! p_star
    DOUBLE PRECISION :: ps_star
! void fraction at maximum packing
    DOUBLE PRECISION :: epg_star
!----------------------------------------------------------------------
! epg = ep_star - (p_star/a_ps)**(1/b_ps)
    INV_H = epg_star - (ps_star/a_ps)**(ONE/dble(b_ps))

  END FUNCTION INV_H

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! dP_s/dEP_s                                                           !
!                                                                      !
! Differentiates the solids plastic-flow pressure model with respect   !
! to solids volume fraction for MFIX default model                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
  DOUBLE PRECISION FUNCTION dPodEP_s(eps,epg_star)
    IMPLICIT NONE
! solids volume fraction
    DOUBLE PRECISION :: eps
! void fraction at maximum packing
    DOUBLE PRECISION :: epg_star
!----------------------------------------------------------------------
! Differentiate P_s w.r.t. EP_s.
! d(p_star)/d(eps)=-d(p_star)/d(epg)
    dPodEP_s = a_ps * dble(b_ps)*&
       (MAX(ZERO, (epg_star - (ONE - eps) )))**(b_ps-1)
  END FUNCTION dPodEP_s


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Neg_H = P_star                                                       !
!                                                                      !
! Returns the solids pressure in plastic-flow stress formulation at    !
! the given void fraction and void fraction at maximum packing. Is     !
! Jackson's (1998) model.                                              !
!                                                                      !
! References:                                                          !
! R. Jackson, in: L.-S. Fan, T.M. Knowlton (Eds.), Fluidization IX,    !
!    Engineering Foundation Publication, New York, 1998, p. 1-13.      !
!    Eqn 5 & 9.                                                        !
! A. Srivastava and S. Sundaresan, Role of wall friction in            !
!    fluidization and standpipe flow, Powder Technology, 124, 2002     !
!    45-54.  Eqn 4 & table 2                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!    IMPLICIT NONE
! coefficient in Jackson's model (barye=dyne/cm^2=g/(cm.s^2))
!    DOUBLE PRECISION, PARAMETER :: A_ps_j = 0.5D0*100D0*2D0*981D0*0.4D0
! voidage at random close-packing
!    DOUBLE PRECISION, PARAMETER :: epg_cp = 0.35D0
!----------------------------------------------------------------------
! P_star = a_ps*(EP_star - EP_g)/(EP_g - EP_CP)
!    Neg_H(epg,epg_star) = to_SI*a_ps_jackson * &
!                         (MAX(ZERO, (epg_star-epg)/(epg-epg_cp)))

! epg = (ep_star*a_ps + epg_cp*ps_star)/(ps_star+a_ps)
!    INV_H(ps_star,epg_star) = (epg_star* to_SI*a_ps_jackson + ps_star*epg_cp)/&
!         (ps_star + to_SI*a_ps_jackson)

! dP_s/dEP_s (EP_s)
! Differentiate P_s w.r.t. EP_s.
!    dPodEP_s(XXX,YYY) = to_SI*a_ps_jackson * (YYY - epg_cp)/(XXX - epg_cp)**2

END MODULE solids_pressure
