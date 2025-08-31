#include "error.inc"

MODULE CHECK_GAS_PHASE_MOD

   USE error_manager
   USE read_database_mod, only: read_database

CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_GAS_PHASE                                        !
!  Purpose: Check the gas phase input section                          !
!                                                                      !
!  Author: P.Nicoletti                                Date: 02-DEC-91  !
!          J.Musser                                   Date: 01-FEB-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_GAS_PHASE(MFIX_DAT)


! Global Variables:
!---------------------------------------------------------------------//
! Flag: Solve species equations.
      use run, only: SPECIES_EQ, ENERGY_EQ
! Flag: Use legacy reaction rates implementation
      use rxns, only: USE_RRATES
! User specified: Constant gas viscosity
      use physprop, only: MU_G0
! User specified: Constant gas thermal conductivity
      use physprop, only: K_G0
! User specified: Constant gas mixture diffusion coefficient
      use physprop, only: DIF_G0
! User specified: Constant gas specific heat
      use physprop, only: C_PG0
! User specified: Constant gas density
      use physprop, only: RO_G0
! User specified: Constant gas mixture molecular weight
      use physprop, only: MW_AVG
! Flag: Use Multicomponent diffusion model
      use physprop, only: Multi_Component_Diffusion
! Flag: to compute diffusion coefficients from kinetic theory
      use physprop, only: dif_coeff_kt
! Flag: to include the effect of thermal diffusion
      use physprop, only: dif_thermal
! User specified: Parameter in Lennard-Jones intermolecular potential model
      use physprop, only: LJsig
! User specified: Parameter in Lennard-Jones intermolecular potential model
      use physprop, only: LJeps
! User specified: Multicomponent diffusion coefficients
      use physprop, only: Dabg
! User specified: total number of gas species
      use physprop, only: NMAX

      use mms, only: use_mms
! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants
      use param1, only: UNDEFINED, ZERO

! Skip data check when doing preprocessing only
      USE run, only:ppo

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: MFIX_DAT

! Local Variables:
!---------------------------------------------------------------------//
! NONE
      INTEGER :: I, J

!......................................................................!

      IF(PPO) RETURN

! CHECK C_pg0
      IF (C_PG0 < ZERO) THEN
         WRITE(ERR_MSG,1001) 'C_PG0', iVal(C_PG0)
         CALL LOG_ERROR()
      ENDIF

! CHECK DIF_g0
      IF (DIF_G0 < ZERO) THEN
         WRITE(ERR_MSG,1001) 'DIF_g0', iVal(DIF_g0)
         CALL LOG_ERROR()
      ENDIF

! Check the input specifications for gas species.
      IF(USE_RRATES)THEN
         CALL CHECK_GAS_SPECIES_LEGACY
      ELSE
         CALL CHECK_GAS_SPECIES(MFIX_DAT)
      ENDIF

! Currently MMS uses constant properties. These are in place simply
! to give the developer a heads-up that the code/setup may not fully
! encompass the use of non-constant properties
      IF (USE_MMS) THEN
         IF (MU_G0 == UNDEFINED) THEN
            WRITE(ERR_MSG, 1200) 'MU_G0'
            CALL LOG_ERROR()
         ENDIF
         IF (K_G0 == UNDEFINED .AND. ENERGY_EQ) THEN
            WRITE(ERR_MSG, 1200) 'K_G0'
            CALL LOG_ERROR()
         ENDIF
         IF (DIF_G0 == UNDEFINED .AND. SPECIES_EQ(0)) THEN
            WRITE(ERR_MSG, 1200) 'DIF_G0'
            CALL LOG_ERROR()
         ENDIF

 1200 FORMAT('Error 1200: ',A,' must be defined when USE_MMS is T.')

      ENDIF

! CHECK MW_AVG
      IF (SPECIES_EQ(0)) THEN
! MW_AVG is defined and the gas phase species equations are solved, then
! the user specified average molecular weight is ignored. The gas phase
! mixture molecular weight (MW_MIX_g) is used instead.
         IF (MW_AVG /= UNDEFINED) THEN
            WRITE (ERR_MSG, 1100) 'solving species equations'
            CALL LOG_WARNING()
            MW_AVG = UNDEFINED
         ENDIF
! Checks for multicomponent diffusion model
         IF (Multi_Component_Diffusion) THEN
	    IF (dif_coeff_kt) THEN
	       IF (RO_G0 /= UNDEFINED) THEN
                  WRITE(ERR_MSG, *) 'Constant gas density detected. Ideal gas law must be used &
                                    to compute diffusion coefficients from kinetic theory.'
                  CALL LOG_ERROR()
	       ENDIF
	       DO I = 1, NMAX(0)
		  IF (LJsig(I) == UNDEFINED) THEN
                     WRITE(ERR_MSG, 1002) 'LJsig', I
                     CALL LOG_ERROR()
		  ENDIF
		  IF (LJeps(I) == UNDEFINED) THEN
                     WRITE(ERR_MSG, 1002) 'LJeps', I
                     CALL LOG_ERROR()
		  ENDIF
	       ENDDO
	    ELSE
	       DO I = 1, NMAX(0)
	          DO J = (I+1), NMAX(0)
		     IF (Dabg(I,J) == UNDEFINED) THEN
                        WRITE(ERR_MSG, 1003) 'Dabg', I, J
                        CALL LOG_ERROR()
		     ENDIF
	          ENDDO
	       ENDDO
	    ENDIF
	    IF (dif_thermal .and. (.not.energy_eq)) THEN
               WRITE (ERR_MSG, *) 'Thermal diffusion is ineffective without solving for energy'
               CALL LOG_WARNING()
	    ENDIF
         ENDIF
      ELSE
! When the species equations are not solved and the gas phase is
! compressible, verify that the user provided average molecular weight
! has a physical value. (This does not include the case where MW_AVG
! is UNDEFINED.)
         IF (RO_G0 == UNDEFINED) THEN
            IF (MW_AVG <= ZERO) THEN
               WRITE(ERR_MSG, 1001) 'MW_AVG', iVal(MW_AVG)
               CALL LOG_ERROR()
            ENDIF
         ELSE
! Gas density for incompressible flows must be positive.
            IF (RO_G0 < ZERO) THEN
               WRITE(ERR_MSG, 1001) 'RO_G0', iVal(RO_G0)
               CALL LOG_ERROR()
            ENDIF
! Incompressible simulations do not need MW_AVG. Notify the user that
! the provided data is ignored.
            IF (MW_AVG /= UNDEFINED)THEN
               WRITE(ERR_MSG, 1100) 'RO_g0 is specified'
               CALL LOG_WARNING()
            ENDIF

         ENDIF
      ENDIF

      ! CHECK Gas viscosity
      CALL CHECK_GAS_VISCOSITY
      ! IF (MU_G0 < ZERO) THEN
      !    WRITE(ERR_MSG,1001) 'MU_G0', iVal(MU_G0)
      !    CALL LOG_ERROR()
      ! ENDIF

! CHECK Gas thermoal conductivity
      CALL CHECK_GAS_THERMAL_CONDUCTIVITY

      RETURN

 1001 FORMAT('Error 1001: Invalid or unknown input: ',A,' = ',A)
 1002 FORMAT('Error 1002: Input ',A,'(', I3, ') is undefined.')
 1003 FORMAT('Error 1003: Input ',A,'(', I3, ',', I3, ') is undefined.')
 1100 FORMAT('Message 2000: MW_AVG is not needed when ',A,'.')

      END SUBROUTINE CHECK_GAS_PHASE



!----------------------------------------------------------------------!
! Subroutine: CHECK_GAS_SPECIES                                        !
! Purpose: Gas phase species checks.                                   !
!                                                                      !
! Author: J. Musser                                  Date: 07-FEB-14   !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_GAS_SPECIES(MFIX_DAT)


! Global Variables:
!---------------------------------------------------------------------//
! Flag: Solve energy equations
      use run, only: ENERGY_EQ
! Flag: Solve species equations
      use run, only: SPECIES_EQ
! Flag: Database for phase X was read for species Y
      use rxns, only: rDatabase
! Flag: Code is reinitializing
      use run, only: REINITIALIZING
! Gas phase species database names.
      use rxns, only: SPECIES_g
! Gas phase molecular weights.
      use physprop, only: MW_g
! Number of gas phase species.
      use physprop, only: NMAX, NMAX_g
! User specified: Constant gas phase specific heat
      use physprop, only: C_PG0
! User specified: Constant gas density
      use physprop, only: RO_G0
! User specified: Constant gas phase mixture molecular weight
      use physprop, only: MW_AVG
! Check for granular flow
      use discretelement, only: DISCRETE_ELEMENT

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of gas phase species.
      USE param, only: DIM_N_g
! Constants.
      USE param1, only: UNDEFINED, UNDEFINED_I, UNDEFINED_C
      USE param1, only: ZERO

      implicit none

      CHARACTER(LEN=*), INTENT(IN) :: MFIX_DAT

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter.
      INTEGER :: N

! Flag that the energy equations are solved and constant gas phase
! specific heat is undefined.
! If true, a call to the thermochemical database is made.
      LOGICAL EEQ_CPG

! Flag that the average molecular weight (MW_AVG) and constant gas
! phase density are undefined.
! If true, a call to the thermochemical database is made.
      LOGICAL MWg_ROg

! Flag that the gas phase species equations are solved and the
! molecular weight for a species is not given in the data file.
! If true, a call to the thermochemical database is made.
      LOGICAL SEQ_MWg


!......................................................................!

! Reconcile the new species input method with the legacy input method.
      IF(SPECIES_EQ(0)) THEN

         IF(NMAX_g == UNDEFINED_I) THEN
            WRITE(ERR_MSG,1000) 'NMAX_g'
            CALL LOG_ERROR()
         ELSEIF(NMAX_g > DIM_N_G) THEN
            WRITE(ERR_MSG,1001) 'NMAX_g', iVal(NMAX_g)
            CALL LOG_ERROR()
         ELSE
            NMAX(0) = NMAX_g
         ENDIF
! Set the number of species to one if the species equations are not solved and
! the number of species is not specified.
      ELSE
         NMAX(0) = merge(1, NMAX_g, NMAX_g == UNDEFINED_I)
      ENDIF

! Set a default value of C_PGO, where running a granular Lagrangian
! simulation with energy equation if no gas species are set.
! The data check will fail if C_PG0 is not set and it is not convenient
! to set in the GUI because the Fluid pane is disabled.

      IF(DISCRETE_ELEMENT.AND.RO_G0==ZERO) THEN ! Granular flow
         IF(ENERGY_EQ .AND. C_PG0 == UNDEFINED.AND.(.NOT.SPECIES_EQ(0))) THEN
            C_PG0 = 1005.0D0
            WRITE(ERR_MSG,1500)
            CALL LOG_WARNING()
         ENDIF
      ENDIF

 1500 FORMAT(&
      'Message: 1500 This is a granular flow. The energy equations is being solved (ENERGY_EQ) and'     /&
      'the constant gas specific heat is undefined (C_PG0). The thermochemical' /&
      'A default value of C_PG0 = 1005.0 J/(kg.K) has been set.')


! Flag that the energy equations are solved and specified gas phase
! specific heat is undefined.
      EEQ_CPG = (ENERGY_EQ .AND. C_PG0 == UNDEFINED)
      IF(EEQ_CPG .AND. .NOT.REINITIALIZING) THEN
         WRITE(ERR_MSG,2000)
         CALL LOG_INFO()
      ENDIF

 2000 FORMAT(&
      'Message: 2000 The energy equations are being solved (ENERGY_EQ) and'     /&
      'the constant gas specific heat is undefined (C_PG0). The thermochemical' /&
      'database will be used to gather specific heat data on the individual gas'/&
      'phase species.')

      MWg_ROg = .FALSE.
      SEQ_MWg = .FALSE.
      IF(MW_AVG == UNDEFINED) THEN
         DO N=1,NMAX(0)
            IF(MW_g(N) == UNDEFINED) THEN
               IF(RO_G0 == UNDEFINED) MWg_ROg = .TRUE.
               IF(SPECIES_EQ(0)) SEQ_MWg = .TRUE.
            ENDIF
         ENDDO
      ENDIF

      IF(MWg_ROg .AND. REINITIALIZING) THEN
         WRITE(ERR_MSG, 2001)
         CALL LOG_INFO()
      ENDIF

 2001 FORMAT(&
      'Message 2001: MW_AVG and RO_G0 are undefined and one or more species'/&
      'molecular weights are undefined. The thermochemical database will be'/&
      'used to gather missing molecular weight data.')

      IF(SEQ_MWg .AND. REINITIALIZING) THEN
         WRITE(ERR_MSG, 2002)
         CALL LOG_INFO()
      ENDIF

 2002 FORMAT(&
      'Message 2002: One or more species molecular weights are undefined and' /&
      'the gas phase species equations are being solved (SPECIES_EQ(0)). The' /&
      'thermochemical database will be used to gather missing molecular weight'/&
      'data.')

! Initialize flag indicating the database was read for a species.
      rDatabase(0,:) = .FALSE.

      IF(EEQ_CPG .OR. SEQ_MWg .OR. MWg_ROg) THEN

         IF(.NOT.REINITIALIZING) THEN
            WRITE(ERR_MSG, 3000)
            CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, FOOTER=.FALSE.)
         ENDIF

 3000    FORMAT('Message 3000: Searching thermochemical databases for gas phase'/&
         'species data.',/'  ')

         DO N = 1, NMAX(0)
            IF(EEQ_CPG .OR. MW_g(N) == UNDEFINED) THEN
! Notify the user of the reason the thermochemical database is used.
! Flag that the species name is not provided.
               IF(SPECIES_g(N) == UNDEFINED_C) THEN
                  WRITE(ERR_MSG,1000) trim(iVar('SPECIES_g',N))
                  CALL LOG_ERROR()
               ENDIF
! Update the log files.
               IF(.NOT.REINITIALIZING) THEN
                  WRITE(ERR_MSG, 3001) N, trim(SPECIES_g(N))
                  CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)
               ENDIF
               3001 FORMAT(/2x,'>',I3,': Species: ',A)
! Read the database.
               CALL READ_DATABASE(MFIX_DAT, 0, N, SPECIES_g(N), MW_g(N))
! Flag variable to stating that the database was read.
               rDatabase(0,N) = .TRUE.
            ENDIF

         ENDDO ! Loop over species

         IF(.NOT.REINITIALIZING) THEN
            CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE.)
         ENDIF

      ENDIF

! Verify that no additional species information was given.
      DO N = NMAX(0) + 1, DIM_N_G
         IF(MW_G(N) /= UNDEFINED) THEN
            WRITE(ERR_MSG, 1002) trim(iVar('MW_g',N))
            CALL LOG_ERROR()
         ENDIF
      ENDDO

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A)

 1001 FORMAT('Error 1001: Invalid or unphysical input: ',A,' = ',A)

 1002 FORMAT('Error 1002: Invalid input: ',A,' specified out of range.')

      END SUBROUTINE CHECK_GAS_SPECIES


!----------------------------------------------------------------------!
! Subroutine: CHECK_GAS_SPECIES_LEGACY                                 !
! Purpose: These are legacy checks for using rrates.f to specify       !
! chemical reactions.                                                  !
!                                                                      !
! Author: J. Musser                                  Date: 03-FEB-14   !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_GAS_SPECIES_LEGACY


! Global Variables:
!---------------------------------------------------------------------//
! Flag: Solve species equations
      use run, only: SPECIES_EQ
! Gas phase molecular weights.
      use physprop, only: MW_g
! Number of gas phase species.
      use physprop, only: NMAX, NMAX_g
! Flag: Database was read. (legacy)
      use physprop, only: DATABASE_READ


! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of gas phase species.
      USE param, only: DIM_N_g
! Constants.
      USE param1, only: UNDEFINED_I, UNDEFINED, ZERO

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter.
      INTEGER :: N

!......................................................................!

! Reconcile the new species input method with the legacy input method.
      IF(SPECIES_EQ(0)) THEN
! Legacy checks for species equations.
         IF(NMAX_g /= UNDEFINED_I) THEN
            WRITE(ERR_MSG,2000) 'NMAX_g', 'undefined'
            CALL LOG_ERROR()
         ELSEIF(NMAX(0) == UNDEFINED_I) THEN
            WRITE(ERR_MSG,2000) trim(iVar('NMAX',0)), 'specified'
            CALL LOG_ERROR()
         ELSEIF(NMAX(0) > DIM_N_G) THEN
            WRITE(ERR_MSG,1001) trim(iVar('NMAX',0)), iVal(NMAX(0))
            CALL LOG_ERROR()
         ENDIF
! Set the number of species to one if the species equations are not
! solved and the number of species is not specified.
      ELSE
         IF(NMAX(0) == UNDEFINED_I) NMAX(0) = 1
      ENDIF

! Check MW_g if solids species are present
      DO N = 1, NMAX(0)
         IF(MW_G(N) == UNDEFINED) THEN
            WRITE(ERR_MSG,2000)trim(iVar('MW_g',N)), 'specified'
            CALL LOG_ERROR()
         ELSEIF(MW_G(N) <= ZERO) THEN
            WRITE(ERR_MSG,1001)trim(iVar('MW_g',N)), iVal(MW_G(N))
            CALL LOG_ERROR()
         ENDIF
      ENDDO ! Loop over species
      DO N = NMAX(0) + 1, DIM_N_G
         IF(MW_G(N) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1001)trim(iVar('MW_g',N)), iVal(MW_G(N))
            CALL LOG_ERROR()
         ENDIF
      ENDDO

! Set the legacy database flag. (Also in check_solids_common_all)
      DATABASE_READ = .FALSE.

      RETURN

 1001 FORMAT('Error 1001: Invalid or unphysical input: ',A,' = ',A)

 2000 FORMAT('Error 2000: Invalid input. ',A,' must be ',A,/'when ',    &
         'USE_RRATES is .TRUE.')

   END SUBROUTINE CHECK_GAS_SPECIES_LEGACY

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_GAS_VISOCITY                                      !
!  Purpose: Check the gas visocity model options                       !
!                                                                      !
!  Author: Jeff Dietiker                              Date: 11-JAN-24  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_GAS_VISCOSITY

! Global Variables:
!---------------------------------------------------------------------//
  use derived_types, only: MU_G_MODEL
  use derived_types, only: MU_G_MODEL_ENUM
  use derived_types, only: MU_G_CONSTANT
  use derived_types, only: MU_G_USR
  use derived_types, only: MU_G_SUTHERLAND
  use derived_types, only: MU_G_HERSCHEL_BULKLEY
  use derived_types, only: MU_G_LENNARD_JONES
  use derived_types, only: MU_G_CANTERA_POLY

  use physprop, only: MU_G0, NMAX
  use physprop, only: SL_muref, SL_Tref, SL_S
  use physprop, only: HB_tau0, HB_k0, HB_n, HB_gama_c
  use physprop, only: LJeps, LJsig, poly_mu_g

  use run, only: units, species_eq

  use usr_prop, only: usr_mug


! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: UNDEFINED, UNDEFINED_C
      USE param1, only: ONE, ZERO

! Global Module procedures:
!---------------------------------------------------------------------//
      USE error_manager

      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
      INTEGER :: I, J

!......................................................................!

! For backward compatibility, set the viscosity model to:
! - 'CONSTANT' if MU_G0 is NOT UNDEFINED
! - 'USR' if USR_MUG is TRUE
!
! Default model is Sutherland model


      if(MU_G0/=UNDEFINED) MU_G_MODEL='CONSTANT'
      if(USR_MUG) MU_G_MODEL='USR'

      SELECT CASE(TRIM(MU_G_MODEL))
      CASE('CONSTANT')
         MU_G_MODEL_ENUM = MU_G_CONSTANT
         IF (MU_G0 < ZERO.OR.MU_G0==UNDEFINED) THEN
            WRITE(ERR_MSG,1001) 'MU_G0', iVal(MU_G0)
            CALL LOG_ERROR()
         ENDIF

      CASE('SUTHERLAND')
         MU_G_MODEL_ENUM = MU_G_SUTHERLAND

         IF(SL_MUREF==UNDEFINED) SL_MUREF = merge(1.715D-5, 1.715D-4, UNITS =='SI')

         IF (SL_MUREF < ZERO) THEN
            WRITE(ERR_MSG,1001) 'SL_MUREF', iVal(SL_MUREF)
            CALL LOG_ERROR()
         ENDIF
         IF (SL_TREF < ZERO) THEN
            WRITE(ERR_MSG,1001) 'SL_TREF', iVal(SL_TREF)
            CALL LOG_ERROR()
         ENDIF
         IF (SL_S < ZERO) THEN
            WRITE(ERR_MSG,1001) 'SL_S', iVal(SL_S)
            CALL LOG_ERROR()
         ENDIF

      CASE('HERSCHEL-BULKLEY', 'HERSCHEL_BULKLEY')
         MU_G_MODEL_ENUM = MU_G_HERSCHEL_BULKLEY

         IF (HB_tau0 < ZERO) THEN
            WRITE(ERR_MSG,1001) 'HB_tau0', iVal(HB_TAU0)
            CALL LOG_ERROR()
         ENDIF

         IF (HB_k0 < ZERO) THEN
            WRITE(ERR_MSG,1001) 'HB_k0', iVal(HB_k0)
            CALL LOG_ERROR()
         ENDIF

         IF (HB_n < ZERO) THEN
            WRITE(ERR_MSG,1001) 'HB_n', iVal(HB_n)
            CALL LOG_ERROR()
         ENDIF

         IF (HB_gama_c < ZERO) THEN
            WRITE(ERR_MSG,1001) 'HB_gama_c', iVal(HB_gama_c)
            CALL LOG_ERROR()
         ENDIF
      CASE('LENNARD_JONES')
         MU_G_MODEL_ENUM = MU_G_LENNARD_JONES
         IF(.NOT. SPECIES_EQ(0)) THEN
            WRITE(ERR_MSG, 1004) 'LENNARD_JONES'
            CALL LOG_ERROR()
         ENDIF
         DO I = 1, NMAX(0)
            IF (LJsig(I) == UNDEFINED) THEN
               WRITE(ERR_MSG, 1002) 'LJsig', I
               CALL LOG_ERROR()
            ENDIF
            IF (LJeps(I) == UNDEFINED) THEN
               WRITE(ERR_MSG, 1002) 'LJeps', I
               CALL LOG_ERROR()
            ENDIF
         ENDDO
      CASE('CANTERA_POLY')
         MU_G_MODEL_ENUM = MU_G_CANTERA_POLY
         IF(.NOT. SPECIES_EQ(0)) THEN
            WRITE(ERR_MSG, 1004) 'CANTERA_POLY'
            CALL LOG_ERROR()
         ENDIF
         ! We only check the polynomial coefficients here.
         ! We do check the parameters from transport, like LJsig, LJeps, mu_transport, et al. here
         ! because they have been checked in GUI when getting the fitted polynominals from Cantera.
         DO I = 1, NMAX(0)
            DO J = 1, 5
               IF (poly_mu_g(I,J) == UNDEFINED) THEN
                  WRITE(ERR_MSG, 1003) 'poly_mu_g', I, J, ival(poly_mu_g(I,J))
                  CALL LOG_ERROR()
               ENDIF
            ENDDO
         ENDDO

      CASE('USR', 'MU_G_USR')
         MU_G_MODEL_ENUM = MU_G_USR

      CASE DEFAULT
            WRITE(ERR_MSG, 1100)
            CALL LOG_ERROR()
      END SELECT

 1001 FORMAT('Error 1001: Invalid or unknown input: ',A,' = ',A)
 1002 FORMAT('Error 1002: Input ',A,'(', I3, ') is undefined.')
 1003 FORMAT('Error 1003: Input ',A,'(', I3, ',', I1, ') is undefined.')
 1004 FORMAT('Error 1004: Species equation for gas phase needs to be turned on ', / &
             'to use', A,' model for viscosity.')
 1100 FORMAT('Error 1100: Unknown MU_G_MODEL')

      RETURN
      END SUBROUTINE CHECK_GAS_VISCOSITY

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_GAS_THERMAL_CONDUCTIVITY                          !
!  Purpose: Check the gas thermal conductivity model options           !
!                                                                      !
!  Author: Hang Zhou                                  Date: 25-Apr-25  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_GAS_THERMAL_CONDUCTIVITY

! Global Variables:
!---------------------------------------------------------------------//
  use derived_types, only: KG_MODEL
  use derived_types, only: KG_MODEL_ENUM
  use derived_types, only: KG_CONSTANT
  use derived_types, only: KG_USR
  use derived_types, only: KG_AIR
  use derived_types, only: KG_LENNARD_JONES
  use derived_types, only: KG_CANTERA_POLY

  use physprop, only: K_G0, NMAX
  use physprop, only: LJeps, LJsig, poly_kg

  use run, only: species_eq

  use usr_prop, only: usr_kg


! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: UNDEFINED
      USE param1, only: ZERO

! Global Module procedures:
!---------------------------------------------------------------------//
      USE error_manager

      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
      INTEGER :: I, J

!......................................................................!

! For backward compatibility, set the thermal conductivity model to:
! - 'CONSTANT' if K_G0 is NOT UNDEFINED
! - 'USR' if USR_KG is TRUE
!
! Default model is temperature-dependent model for air

      if(K_G0/=UNDEFINED) KG_MODEL='CONSTANT'
      if(USR_KG) KG_MODEL='USR'

      SELECT CASE(TRIM(KG_MODEL))
      CASE('CONSTANT')
         KG_MODEL_ENUM = KG_CONSTANT
         IF (K_G0 < ZERO.OR.K_G0==UNDEFINED) THEN
            WRITE(ERR_MSG,1001) 'K_G0', iVal(K_G0)
            CALL LOG_ERROR()
         ENDIF

      CASE('AIR')
         KG_MODEL_ENUM = KG_AIR

      CASE('LENNARD_JONES')
         KG_MODEL_ENUM = KG_LENNARD_JONES
         IF(.NOT. SPECIES_EQ(0)) THEN
            WRITE(ERR_MSG, 1004) 'LENNARD_JONES'
            CALL LOG_ERROR()
         ENDIF
         DO I = 1, NMAX(0)
            IF (LJsig(I) == UNDEFINED) THEN
               WRITE(ERR_MSG, 1002) 'LJsig', I
               CALL LOG_ERROR()
            ENDIF
            IF (LJeps(I) == UNDEFINED) THEN
               WRITE(ERR_MSG, 1002) 'LJeps', I
               CALL LOG_ERROR()
            ENDIF
         ENDDO
      CASE('CANTERA_POLY')
         KG_MODEL_ENUM = KG_CANTERA_POLY
         IF(.NOT. SPECIES_EQ(0)) THEN
            WRITE(ERR_MSG, 1004) 'CANTERA_POLY'
            CALL LOG_ERROR()
         ENDIF
         ! We only check the polynomial coefficients here.
         ! We do check the parameters from transport, like LJsig, LJeps, mu_transport, et al. here
         ! because they have been checked in GUI when getting the fitted polynominals from Cantera.
         DO I = 1, NMAX(0)
            DO J = 1, 5
               IF (poly_kg(I,J) == UNDEFINED) THEN
                  WRITE(ERR_MSG, 1003) 'poly_kg', I, J, ival(poly_kg(I,J))
                  CALL LOG_ERROR()
               ENDIF
            ENDDO
         ENDDO

      CASE('USR', 'KG_USR')
         KG_MODEL_ENUM = KG_USR


      CASE DEFAULT
            WRITE(ERR_MSG, 1100)
            CALL LOG_ERROR()
      END SELECT

 1001 FORMAT('Error 1001: Invalid or unknown input: ',A,' = ',A)
 1002 FORMAT('Error 1002: Input ',A,'(', I3, ') is undefined.')
 1003 FORMAT('Error 1003: Input ',A,'(', I3, ',', I1, ') is undefined.')
 1004 FORMAT('Error 1004: Species equation for gas phase needs to be turned on ', / &
             'to use', A,' model for thermal conductivity.')
 1100 FORMAT('Error 1100: Unknown KG_MODEL')

      RETURN

      END SUBROUTINE CHECK_GAS_THERMAL_CONDUCTIVITY
END MODULE CHECK_GAS_PHASE_MOD
