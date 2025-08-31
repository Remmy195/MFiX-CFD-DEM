MODULE DES_THERMO_NEWVALUES_MOD

   use des_reaction_model_mod, only: des_reaction_model

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_THERMO_NEWVALUES                                   !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE DES_THERMO_NEWVALUES

      USE stiff_chem, only: stiff_chemistry
      USE des_thermo, only: q_source, q_source0, des_t_s, des_c_ps
      USE des_thermo_cond, only: DES_QW_cond
      USE discretelement, only: intg_euler, intg_adams_bashforth, max_pip, des_explicitly_coupled
      USE discretelement, only: normal_particle, particle_state, pmass, dtsolid
      USE param1, only: zero
      USE run, only: ANY_SOLIDS_SPECIES_EQ, ENERGY_EQ

      IMPLICIT NONE

! Passed variables
!-----------------------------------------------
! NONE

! Local variables
!---------------------------------------------------------------------//
! Logical for Adams-Bashfort integration.
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
!---------------------------------------------------------------------//

      IF(ENERGY_EQ) THEN

! Second-order Adams-Bashforth scheme defaults to Euler on first pass.
         IF(FIRST_PASS .AND. INTG_ADAMS_BASHFORTH) THEN
            WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE) &
                 Q_Source0(:MAX_PIP) = Q_Source(:MAX_PIP)/       &
                 (PMASS(:MAX_PIP)*DES_C_ps(:MAX_PIP))
         ENDIF
         FIRST_PASS = .FALSE.

! First-order method
         IF (INTG_EULER) THEN
            WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE)   &
                 DES_T_s(:MAX_PIP) = DES_T_s(:MAX_PIP) +   &
                 DTSOLID*(Q_Source(:MAX_PIP)/(PMASS(:MAX_PIP)*     &
                 DES_C_ps(:MAX_PIP)))

! Second-order Adams-Bashforth scheme
         ELSE
            WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE)
               DES_T_s(:MAX_PIP) = DES_T_s(:MAX_PIP) + DTSOLID *         &
                    (1.5d0*Q_Source(:MAX_PIP) -0.5d0*Q_Source0(:MAX_PIP))/ &
                    (PMASS(:MAX_PIP)*DES_C_ps(:MAX_PIP))
               Q_Source0(:MAX_PIP) = Q_Source(:MAX_PIP)
            ENDWHERE
         ENDIF

         Q_Source(:) = ZERO
         IF(ALLOCATED(DES_QW_Cond)) DES_QW_Cond(:,:) = ZERO

      ENDIF

! Update particle from reactive chemistry process.
      IF(ANY_SOLIDS_SPECIES_EQ .AND. .NOT.DES_EXPLICITLY_COUPLED .AND. .NOT.STIFF_CHEMISTRY)&
         CALL DES_REACTION_MODEL

      RETURN

   END SUBROUTINE DES_THERMO_NEWVALUES
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine name: GSP_HEAT_COND_EXPLICIT                             !
!                                                                      !
!  Purpose:                                                            !
!     calculate internal heat transfer for the GSP                     !
!                                                                      !
!  Author: R.Ke                                       Date: 13-Feb-24  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE GSP_HEAT_COND_EXPLICIT

! Module
!---------------------------------------------------------------------//
      Use des_rxns
      Use des_thermo
      Use discretelement
      Use run
      Use usr
      Use physprop, only: k_s0
      use compar, only: pe_io, mype, numPEs

      IMPLICIT NONE
! Local variables
!---------------------------------------------------------------------//
      double precision :: qcond, kp, kp_wood, kp_char, T_L
      integer :: L, J, nb, nbIDX, lcurpos, locpp
      double precision :: Kp_1, Kp_2, kps_eff

      DOUBLE PRECISION,save :: outTime=0
      DOUBLE PRECISION,parameter :: outDT = 1.0
      CHARACTER(LEN=64) :: FNAME
      LOGICAL :: F_EXISTS
      INTEGER, PARAMETER :: lUNIT = 2030
      DOUBLE PRECISION :: dist(3), eff_dist

! calculate kp for heat conduction in gsp explicit model
      do L=1, PIP
         IF (PARTICLE_STATE(L) .ne. NORMAL_PARTICLE) CYCLE
	      T_L = des_t_s(L)

         !llu's example, calculate kp from species
         !kp_wood=0.13+0.0003*(T_L-273.15)
         !kp_char=0.08-0.0001*(T_L-273.15)
         !Kp = kp_wood*DES_X_s(L,Biomass)+kp_char*( DES_X_s(L,char1)+ DES_X_s(L,char2))
         !gp_kps(L) = Kp

! phase-dependent constant conductivity coefficient phase dependent
         gp_kps(L) = k_s0(PIJK(L,5))
      enddo

! heat conduction at the level of component spheres
      do L=1, PIP
         IF (PARTICLE_STATE(L) .ne. NORMAL_PARTICLE) CYCLE
         qcond = 0
	      T_L = des_t_s(L)

         do nbIDX=1,6
            J = gp_neighsid(L,nbIDX)
            if(J>0) then

! this is based on assumption:
! if a PE have some components sphere of a gsp marking as normal, then it should also have all other components sphere marked as ghost
! so lcurpos can always be found in iglobal_id
! so the only outlier is a real long gsp that no only across boundary, but also across the two ghost layer
! in this outlier case, a PE know some normal, some ghost, the rest just nonexist to this pe
! so the FINDLOC will lead to 0
               !lcurpos = FINDLOC(iglobal_id,J,dim=1)
               DO locpp = 1, size(iglobal_id)
                  IF(iglobal_id(locpp) == J) THEN
                     lcurpos = locpp
                     EXIT
                  ENDIF
               ENDDO

	            if(lcurpos > 0) then
                  ! effective kps for conguate contact heat transfer
                  kps_eff = 0.0
                  Kp_1 = gp_kps(L)
                  Kp_2 = gp_kps(lcurpos)
                  kps_eff = 2.0d0 * (Kp_1*Kp_2)/(Kp_1 + Kp_2)
                  eff_dist = 0.0
                  dist(:) = 0.0
                  dist(:) = abs(sc2gpc_vec(L,:) - sc2gpc_vec(lcurpos,:))

                  ! heat transfer distance is perpendicular to the contact area
                  if(nbIDX == 1 .or. nbIDX == 2) then
                     eff_dist = dist(1)
                  elseif(nbIDX == 3 .or. nbIDX == 4) then
                     eff_dist = dist(2)
                  elseif(nbIDX == 5 .or. nbIDX == 6) then
                     eff_dist = dist(3)
                  endif

                  qcond = qcond + kps_eff*(des_t_s(lcurpos)-T_L)*gp_neighsa(L,nbIDX)/eff_dist
               endif
	         endif
         enddo
	   Q_Source(L) = Q_Source(L) + qcond
      enddo

      RETURN

   END SUBROUTINE GSP_HEAT_COND_EXPLICIT

END MODULE DES_THERMO_NEWVALUES_MOD
