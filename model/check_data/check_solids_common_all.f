#include "error.inc"

MODULE CHECK_SOLIDS_COMMON_ALL_MOD

  use read_database_mod, only: read_database
  use error_manager

  CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_COMMON_ALL                                 !
!  Purpose: Check the solid phase input that is common to all solids   !
!  phase models.                                                       !
!                                                                      !
!    ****** DO NOT PUT MODEL SPECIFIC CHECKS IN THIS ROUTINE ******    !
!                                                                      !
!  Use the companion routines for checks specific to a particular      !
!  solids phase model:                                                 !
!                                                                      !
!    > CHECK_SOLIDS_CONTINUUM :: TFM solids phase model                !
!    > CHECK_SOLIDS_DEM       :: DEM solids phase model                !
!    > CHECK_SOLIDS_MPPIC     :: MPPIC solids phase model              !
!                                                                      !
!  Author: J.Musser                                  Date: 03-FEB-14   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_ALL(MFIX_DAT)

! Global Variables:
!---------------------------------------------------------------------//
      use run, only: energy_eq, species_eq, ks_model
! Flag: Use legacy rrates implementation.
      use rxns, only: USE_RRATES
! Number of continuum solids phases.
      use physprop, only: SMAX
! Number of discrete (DEM/MPPIC) solids.
      use discretelement, only: DES_MMAX
! User specified: Constant solids specific heat.
      use physprop, only: C_PS0
! User specified: Constant solids thermal conductivity.
      use physprop, only: K_S0
! User specified: Initial solids diameter.
      use physprop, only: D_P0
      use physprop, only: mu_s0, dif_s0
      use mms, only: use_mms

! Global Parameters:
!---------------------------------------------------------------------//
      use discretelement, only: SuperDEM, SQP_a, SQP_b, SQP_c, SQP_m, SQP_n
      use run, only: Solids_model
! Maximum number of solids phases.
      use param, only: DIM_M
! Parameter constants
      use param1, only: UNDEFINED, ZERO, UNDEFINED_C
      use sq_properties_mod


      implicit none

      CHARACTER(LEN=*), INTENT(IN) :: MFIX_DAT

! Local Variables:
!---------------------------------------------------------------------//
! Loop counters.
      INTEGER :: M
! Total number of all solids
      INTEGER :: MMAX_L

!......................................................................!

! Set the number of solids phases to be checked.
      MMAX_L = SMAX + DES_MMAX

! Check D_p0
      DO M = 1, MMAX_L
!        print*,'checking dp_0',M,D_p0(M), SQP_A(M),SQP_B(M),SQP_C(M)
         IF(D_P0(M) == UNDEFINED) THEN
            IF(SOLIDS_MODEL(M)=='SQP') THEN
! For SQP, d_p0 is the minimum bounding sphere diameter and can be computed based on superquadic
! parameters a, b, c, m, and n.
               CALL SQ_min_bounding_sphere(sqp_a(M), sqp_b(M), sqp_c(M), sqp_m(M), sqp_n(M), d_p0(M))
               WRITE(ERR_MSG, 1100) M, trim(iVar('D_p0',M)), d_p0(M)
               CALL LOG_INFO()
            ELSE
               WRITE(ERR_MSG, 1000) trim(iVar('D_p0',M))
               CALL LOG_ERROR()
            ENDIF
         ELSEIF(D_P0(M) <= ZERO)THEN
            WRITE(ERR_MSG, 1001) trim(iVar('D_p0',M)), iVal(D_P0(M))
            CALL LOG_ERROR()
         ENDIF
      ENDDO

      DO M = MMAX_L+1, DIM_M
         IF(D_P0(M) /= UNDEFINED)THEN
            WRITE(ERR_MSG,1002) trim(iVar('D_p0',M))
            CALL LOG_ERROR()
         ENDIF
      ENDDO


! Check K_s0
      DO M=1, MMAX_L
         IF (K_S0(M) < ZERO) THEN
! undefined will pass this check...
            WRITE(ERR_MSG, 1001) trim(iVar('K_s0',M)), iVal(K_s0(M))
            CALL LOG_ERROR()
         ENDIF
      ENDDO

      DO M = MMAX_L+1, DIM_M
         IF(K_s0(M) /= UNDEFINED)THEN
            WRITE(ERR_MSG,1002) trim(iVar('K_s0',M))
            CALL LOG_ERROR()
         ENDIF
      ENDDO

! if energy_eq, then ks_model needs to be defined but
! only for tfm and dem (not pic). further checks are
! made in corresponding routines.
      IF (ENERGY_EQ) THEN
         DO M = MMAX_L+1, DIM_M
            IF(Ks_MODEL(M) /= UNDEFINED_C)THEN
               WRITE(ERR_MSG,1002) trim(iVar('Ks_model',M))
               CALL LOG_ERROR()
            ENDIF
         ENDDO
      ENDIF


! Check C_ps0
      DO M=1, MMAX_L
         IF (C_PS0(M) < ZERO) THEN
            WRITE(ERR_MSG, 1001) trim(iVar('C_ps0',M)), iVal(C_ps0(M))
            CALL LOG_ERROR()
         ENDIF
      ENDDO

      DO M = MMAX_L+1, DIM_M
         IF(C_ps0(M) /= UNDEFINED)THEN
            WRITE(ERR_MSG,1002) trim(iVar('C_ps0',M))
            CALL LOG_ERROR()
         ENDIF
      ENDDO

! Check the input specifications for solids species.
      IF(USE_RRATES)THEN
         CALL CHECK_SOLIDS_SPECIES_LEGACY(MMAX_L)
      ELSE
         CALL CHECK_SOLIDS_SPECIES(MFIX_DAT, MMAX_L)
      ENDIF

! Currently MMS uses constant properties. These are in place simply
! to give the developer a heads-up that the code/setup may not fully
! encompass the use of non-constant properties
      IF (USE_MMS) THEN
         DO M = 1, MMAX_L
            IF (MU_S0(M) == UNDEFINED) THEN
               WRITE(ERR_MSG, 1200) trim(ivar('MU_s0',M))
               CALL LOG_ERROR()
            ENDIF
            IF (K_S0(M) == UNDEFINED .AND. ENERGY_EQ) THEN
               WRITE(ERR_MSG, 1200) trim(ivar('K_s0',M))
               CALL LOG_ERROR()
            ENDIF
            IF (DIF_S0(M) == UNDEFINED .AND. SPECIES_EQ(M)) THEN
               WRITE(ERR_MSG, 1200) trim(ivar('DIF_S0',M))
               CALL LOG_ERROR()
            ENDIF
         ENDDO
 1200 FORMAT('Error 1200: ',A,' must be defined when USE_MMS is set.')

      ENDIF


! Check solids drag model selection.
      CALL CHECK_SOLIDS_DRAG


! Check the solids density input parameters.
      CALL CHECK_SOLIDS_DENSITY(MMAX_L)

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A)

 1001 FORMAT('Error 1001: Invalid or unphysical input: ',A,' = ',A)

 1002 FORMAT('Error 1002: Invalid input: ',A,' specified out of range.')

 1100 FORMAT('Info: Solids phase',I3, ': SQP minimum bounding sphere diameter set to: ',A, ' = ',E12.4)
      END SUBROUTINE CHECK_SOLIDS_COMMON_ALL



!----------------------------------------------------------------------!
! Subroutine: CHECK_SOLIDS_DRAG                                        !
! Purpose: Check solids species input.                                 !
!                                                                      !
! Author: J. Musser                                  Date: 07-FEB-14   !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_SOLIDS_DRAG

! Global Variables:
!---------------------------------------------------------------------//
! User specified drag type, as string and enum
      use run, only: DRAG_TYPE
      use run, only: DRAG_TYPE_ENUM
! Possible DRAG_TYPE_ENUM values:
      use run, only: SYAM_OBRIEN
      use run, only: GIDASPOW
      use run, only: GIDASPOW_PCF
      use run, only: GIDASPOW_BLEND
      use run, only: GIDASPOW_BLEND_PCF
      use run, only: WEN_YU
      use run, only: WEN_YU_PCF
      use run, only: KOCH_HILL
      use run, only: KOCH_HILL_PCF
      use run, only: BVK, TPKKV
      use run, only: HYS
      use run, only: GAO
      use run, only: SARKAR
      use run, only: RADL
      use run, only: TGS,DIFELICE
      use run, only: DIFELICE
      use run, only: DIFELICE_GANSER

	  USE DERIVED_TYPES

!	  use run, only: SQP_DIFELICE_GANSER
!	  use run, only: SQP_DIFELICE_HOLZER_SOMMERFELD


      use run, only: USER_DRAG
      use constant, only : drag_c1, drag_d1
      use param1, only: zero


! Global Parameters:
!---------------------------------------------------------------------//

      implicit none


! Local Variables:
!---------------------------------------------------------------------//
!     NONE

!......................................................................!


      SELECT CASE(trim(adjustl(DRAG_TYPE)))

      CASE ('SYAM_OBRIEN'); DRAG_TYPE_ENUM                    = SYAM_OBRIEN
      CASE ('GIDASPOW'); DRAG_TYPE_ENUM                       = GIDASPOW
      CASE ('GIDASPOW_PCF'); DRAG_TYPE_ENUM                   = GIDASPOW_PCF
      CASE ('GIDASPOW_BLEND'); DRAG_TYPE_ENUM                 = GIDASPOW_BLEND
      CASE ('GIDASPOW_BLEND_PCF'); DRAG_TYPE_ENUM             = GIDASPOW_BLEND_PCF
      CASE ('WEN_YU'); DRAG_TYPE_ENUM                         = WEN_YU
      CASE ('WEN_YU_PCF'); DRAG_TYPE_ENUM                     = WEN_YU_PCF
      CASE ('KOCH_HILL'); DRAG_TYPE_ENUM                      = KOCH_HILL
      CASE ('KOCH_HILL_PCF'); DRAG_TYPE_ENUM                  = KOCH_HILL_PCF
      CASE ('BVK'); DRAG_TYPE_ENUM                            = BVK
      CASE ('TPKKV'); DRAG_TYPE_ENUM                          = TPKKV
      CASE ('HYS'); DRAG_TYPE_ENUM                            = HYS
      CASE ('GAO'); DRAG_TYPE_ENUM                            = GAO
      CASE ('SARKAR'); DRAG_TYPE_ENUM                         = SARKAR
      CASE ('RADL'); DRAG_TYPE_ENUM                           = RADL
      CASE ('TGS'); DRAG_TYPE_ENUM                            = TGS
      CASE ('DIFELICE'); DRAG_TYPE_ENUM                       = DIFELICE
      CASE ('DIFELICE_GANSER'); DRAG_TYPE_ENUM                = DIFELICE_GANSER
      CASE ('SQP_DIFELICE_GANSER'); DRAG_TYPE_ENUM            = SQP_DIFELICE_GANSER
      CASE ('SQP_DIFELICE_HOLZER_SOMMERFELD'); DRAG_TYPE_ENUM = SQP_DIFELICE_HOLZER_SOMMERFELD
      CASE ('QKQY'); DRAG_TYPE_ENUM                           = QKQY
      CASE ('USER_DRAG','USR_DRAG'); DRAG_TYPE_ENUM           = USER_DRAG

      CASE DEFAULT
         WRITE(ERR_MSG,1001)'DRAG_TYPE', trim(adjustl(DRAG_TYPE))
         CALL LOG_ERROR()
      END SELECT

! Check that drag_c1 and drag_d1 (Syam-O'Brien drag law) are not zero
      IF(DRAG_TYPE_ENUM == SYAM_OBRIEN) THEN
         IF(DRAG_C1<=ZERO) THEN
            WRITE(ERR_MSG,2001)'DRAG_C1 must be strictly &
            positive. Current value is ', DRAG_C1
            CALL LOG_ERROR()
         ENDIF
         IF(DRAG_D1<=ZERO) THEN
            WRITE(ERR_MSG,2001)' DRAG_D1 must be strictly &
            positive. Current value is ', DRAG_D1
            CALL LOG_ERROR()
         ENDIF
      ENDIF

      RETURN



 1001 FORMAT('Error 1001: Invalid or unknown input: ',A,' = ',A)
 2001 FORMAT('Error 2001: Invalid input: ',A,G12.5)

      END SUBROUTINE CHECK_SOLIDS_DRAG

!----------------------------------------------------------------------!
! Subroutine: CHECK_SOLIDS_SPECIES                                     !
! Purpose: Check solids species input.                                 !
!                                                                      !
! Author: J. Musser                                  Date: 07-FEB-14   !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_SOLIDS_SPECIES(MFIX_DAT, MMAX_LL)

! Global Variables:
!---------------------------------------------------------------------//
! Flag: Solve energy equations
      use run, only: ENERGY_EQ
! Flag: Solve species equations
      use run, only: SPECIES_EQ
! FLag: Reinitializing the code
      use run, only: REINITIALIZING
! Flag: Database for phase X was read for species Y
      use rxns, only: rDatabase
! Solids phase species database names.
      use rxns, only: SPECIES_s
! Solids phase molecular weights.
      use physprop, only: MW_s
! Number of solids phase species.
      use physprop, only: NMAX, NMAX_s
! User specified: Constant solids specific heat
      use physprop, only: C_PS0

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phase species.
      USE param, only: DIM_N_s
! Parameter constants.
      use param1, only: UNDEFINED, UNDEFINED_I, UNDEFINED_C

      implicit none

! Subroutine Arguments:
!---------------------------------------------------------------------//

      CHARACTER(LEN=*), INTENT(IN) :: MFIX_DAT

! Total number of solids phases
      INTEGER, intent(in) :: MMAX_LL

! Local Variables:
!---------------------------------------------------------------------//

! Flag that the energy equations are solved and specified solids phase
! specific heat is undefined.
! If true, a call to the thermochemical database is made.
      LOGICAL EEQ_CPS

! Flag that the solids phase species equations are solved and the
! molecular weight for a species are not given in the data file.
! If true, a call to the thermochemical database is made.
      LOGICAL SEQ_MWs

! Loop counters.
      INTEGER :: M, N

!......................................................................!

! Reconcile the new species input method with the legacy input method.
      DO M=1, MMAX_LL
         IF(SPECIES_EQ(M)) THEN
            IF(NMAX_s(M) == UNDEFINED_I) THEN
               WRITE(ERR_MSG,1000) iVar('NMAX_s',M)
               CALL LOG_ERROR()
            ELSEIF(NMAX_s(M) > DIM_N_S) THEN
               WRITE(ERR_MSG,1001) trim(iVar('NMAX_s',M)),            &
                  trim(iVal(NMAX_s(M)))
               CALL LOG_ERROR()
            ELSEIF(NMAX_s(M)==0) THEN
               WRITE(ERR_MSG, 2010) M, '(Species equation is solved for this phase).'
               CALL LOG_ERROR()
            ELSE
               NMAX(M) = NMAX_s(M)
            ENDIF

! Set the number of species to one if the species equations are not solved and
! the number of species is not specified.
         ELSE
            NMAX(M) = merge(1, NMAX_s(M), NMAX_s(M) == UNDEFINED_I)
         ENDIF

! Flag that the energy equations are solved and specified solids phase
! specific heat is undefined.
         EEQ_CPS = (ENERGY_EQ .AND. C_PS0(M) == UNDEFINED)
         IF(EEQ_CPS .AND. .NOT.REINITIALIZING)THEN
            WRITE(ERR_MSG,2000)
            CALL LOG_INFO()
            IF(NMAX(M)==0) THEN
               WRITE(ERR_MSG, 2010) M, '(needed to compute specific heat).'
               CALL LOG_ERROR()
            ENDIF
         ENDIF

 2000 FORMAT('Message: 2000 The energy equations are being solved ',   &
         '(ENERGY_EQ) and',/'the constant solids specific heat is ',   &
         'undefined (C_PS0). The',/'thermochemical database ',   &
         'will be used to gather specific heat data on',/'the ',       &
         'individual gas phase species.')

         SEQ_MWs = .FALSE.
         DO N=1,NMAX(M)
            IF(MW_s(M,N) == UNDEFINED) THEN
               IF(SPECIES_EQ(M)) SEQ_MWs = .TRUE.
            ENDIF
         ENDDO

         IF(SEQ_MWs .AND. .NOT.REINITIALIZING) THEN
            WRITE(ERR_MSG, 2001) M
            CALL LOG_INFO()
            IF(NMAX(M)==0) THEN
               WRITE(ERR_MSG, 2010) M, '(needed to compute molecular weight).'
               CALL LOG_ERROR()
            ENDIF
         ENDIF

 2001 FORMAT('Message 2001: One or more species molecular weights are',&
         ' undefined and',/'the solids phase species equations are ',  &
         'solved (SPECIES_EQ(',I2,')). The',/'thermochemical database ', &
         'will be searched for molecular weight data.')


 2010 FORMAT('Error 2010: There are no species defined for solids phase:',I3, &
            /'At least one species must be defined ',A)

! Initialize flag indicating the database was read for a species.
         rDatabase(M,:) = .FALSE.

         IF(EEQ_CPS .OR. SEQ_MWs)THEN

            IF(.NOT.REINITIALIZING) THEN
               WRITE(ERR_MSG, 3000) M
               CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, FOOTER=.FALSE.)
            ENDIF

 3000 FORMAT('Message 3000: Searching thermochemical databases for ',  &
         'solids phase',I3,/'species data.',/'  ')

            DO N = 1, NMAX(M)

! Notify the user of the reason the thermochemical database is used.
               IF(EEQ_CPS .OR. MW_s(M,N) == UNDEFINED) THEN

! Flag that the species name is not provided.
                  IF(SPECIES_s(M,N) == UNDEFINED_C) THEN
                     WRITE(ERR_MSG,1000) trim(iVar('SPECIES_s',M,N))
                     CALL LOG_ERROR()
                  ENDIF
! Update the log files.
                  IF(.NOT.REINITIALIZING) THEN
                     WRITE(ERR_MSG, 3001) N, trim(SPECIES_s(M,N))
                     CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)
                  ENDIF
                  3001 FORMAT(/2x,'>',I3,': Species: ',A)

                  CALL READ_DATABASE(MFIX_DAT, M, N, SPECIES_s(M,N), MW_S(M,N))
! Flag variable to stating that the database was read.
                  rDatabase(M,N) = .TRUE.
               ENDIF
            ENDDO ! Loop over species
            IF(.NOT.REINITIALIZING) CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE.)
         ENDIF

! Verify that no additional species information was given.
         DO N = NMAX(M) + 1, DIM_N_S
            IF(MW_S(M,N) /= UNDEFINED) THEN
               WRITE(ERR_MSG, 1002) trim(iVar('MW_s',M,N))
               CALL LOG_ERROR()
            ENDIF
         ENDDO
      ENDDO ! Loop over solids phases

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A)

 1001 FORMAT('Error 1001: Invalid or unphysical input: ',A,' = ',A)

 1002 FORMAT('Error 1002: Invalid input: ',A,' specified out of range.')

      END SUBROUTINE CHECK_SOLIDS_SPECIES


!----------------------------------------------------------------------!
! Subroutine: CHECK_SOLIDS_SPECIES_LEGACY                              !
! Purpose: These are legacy checks for using rrates.f to specify       !
! chemical reactions.                                                  !
!                                                                      !
! Author: J. Musser                                  Date: 03-FEB-14   !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_SOLIDS_SPECIES_LEGACY(MMAX_LL)


! Global Variables:
!---------------------------------------------------------------------//
! Flag: Solve species equations
      use run, only: SPECIES_EQ
! Solids phase molecular weights.
      use physprop, only: MW_s
! Number of solids phase species.
      use physprop, only: NMAX, NMAX_s
! Flag: Database was read. (legacy)
      use physprop, only: DATABASE_READ

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of gas phase species.
      USE param, only: DIM_N_s
! Constants.
      USE param1, only: UNDEFINED_I, UNDEFINED

      implicit none


! Arguments:
!---------------------------------------------------------------------//
! Total number of solids phases
      INTEGER, intent(in) :: MMAX_LL

! Local Variables:
!---------------------------------------------------------------------//
! Loop counters.
      INTEGER :: M, N


!......................................................................!


! Reconcile the new species input method with the legacy input method.
      DO M=1, MMAX_LL
         IF(SPECIES_EQ(M)) THEN

! Legacy checks for species equations.
            IF(NMAX_s(M) /= UNDEFINED_I) THEN
               WRITE(ERR_MSG,2000) trim(iVar('NMAX_s',M)), 'undefined'
               CALL LOG_ERROR()
            ELSEIF(NMAX(M) == UNDEFINED_I) THEN
               WRITE(ERR_MSG,2000) trim(iVar('NMAX',M)), 'specified'
               CALL LOG_ERROR()
            ELSEIF(NMAX(M) > DIM_N_S) THEN
               WRITE(ERR_MSG,1002) trim(iVar('NMAX',M))
               CALL LOG_ERROR()
            ENDIF

! Set the number of species to one if the species equations are not solved and
! the number of species is not specified.
         ELSE
            IF(NMAX(M) == UNDEFINED) NMAX(M) = 1
         ENDIF
      ENDDO

! Check MW_s if solids species are present
      DO M = 1, MMAX_LL
! Initialize flag indicating the database was read for a species.
         DO N = 1, NMAX(M)
            IF(MW_S(M,N) == UNDEFINED) THEN
               WRITE(ERR_MSG, 2000) trim(iVar('MW_s',M,N)), 'specified'
               CALL LOG_ERROR()
            ENDIF
         ENDDO ! Loop over species
         DO N = NMAX(M) + 1, DIM_N_S
            IF(MW_S(M,N) /= UNDEFINED) THEN
               WRITE(ERR_MSG,1002) trim(iVar('MW_s',M,N))
            ENDIF
         ENDDO
      ENDDO ! Loop over solids phases

! Set the legacy database flag. (Also in check_gas_phase.f)
      DATABASE_READ = .FALSE.

      RETURN

 1002 FORMAT('Error 1002: Invalid input: ',A,' specified out of range.')

 2000 FORMAT('Error 2000: Invalid input. ',A,' must be ',A,/'when ',   &
         'USE_RRATES is set.')

      END SUBROUTINE CHECK_SOLIDS_SPECIES_LEGACY



!----------------------------------------------------------------------!
!  Subroutine: CHECK_SOLIDS_DENSITY                                    !
!  Purpose: check the solid phase density input                        !
!                                                                      !
!  Author: J.Musser                                  Date: 03-FEB-14   !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_SOLIDS_DENSITY(MMAX_LL)


! Global Variables:
!---------------------------------------------------------------------//
! Flag: Solve variable solids density.
      use run, only: SOLVE_ROs
! User specified: constant solids density
      use physprop, only: RO_s0
! Baseline species densities
      use physprop, only: RO_Xs0
! Baseline species mass fractions.
      use physprop, only: X_s0
! Index of inert solids species
      use physprop, only: INERT_SPECIES
! Inert species mass fraction in dilute region
      use physprop, only: DIL_INERT_X_VSD
! Number of solids phase species.
      use physprop, only: NMAX

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants
      use param1, only: ZERO, ONE
      use param1, only: UNDEFINED, UNDEFINED_I
! Maximum number of solids phases.
      use param, only: DIM_M
! Maximum number of solids species.
      use param, only: DIM_N_s

! function to compare numbers
      use toleranc, only: compare
! Solids phase density caMulation
      use eos, ONLY: EOSS0

! Flag to indicate user function for density
      use usr_prop, only: usr_ros

      implicit none


! Arguments:
!---------------------------------------------------------------------//
! Total number of solids phases
      INTEGER, intent(in) :: MMAX_LL

! Local Variables:
!---------------------------------------------------------------------//
      INTEGER :: M, N
! local variable indicating that some variable density model is being used
! either a user defined function or mfix built in variable density model.
      LOGICAL :: SolveAnyRos
! local variable to sum species mass fraction
      DOUBLE PRECISION :: sum_xs
!......................................................................!


! Initialize the flag for variable solids density.
      SOLVE_ROs = .FALSE.
      SolveAnyRos = .FALSE.

! Check each solids phase.
      DO M = 1, MMAX_LL

! Set the flag for variable solids density if any of the input parameters
! are specified.
         IF(INERT_SPECIES(M) /= UNDEFINED_I) SOLVE_ROs(M) = .TRUE.
         DO N=1, DIM_N_s
            IF(RO_Xs0(M,N) /= UNDEFINED) SOLVE_ROs(M) = .TRUE.
            IF(X_S0(M,N) /= UNDEFINED) SOLVE_ROs(M) = .TRUE.
         ENDDO

         SolveAnyROs = (SOLVE_ROS(M) .OR. USR_ROS(M))
! Verify that one -and only one- solids density model is in use.
         IF(RO_S0(M) == UNDEFINED .AND. .NOT.SolveAnyROs) THEN
            WRITE(ERR_MSG, 1100) M
            CALL LOG_ERROR()

 1100 FORMAT('Error 1101: No solids density information for phase ',  &
         I2,'.')

! Check if the constant solids phase density is physical.
         ELSEIF(RO_S0(M) /= UNDEFINED .AND. SolveAnyROs) THEN
            WRITE(ERR_MSG, 1101) M
            CALL LOG_ERROR()

 1101 FORMAT('Error 1101: Conflicting solids density input specified ',&
         'for solids',/'phase ',I2,'. Constant solids density ',       &
         'specified (RO_s0) along with one',/'or more of the variable',&
         ' solids density parameters:',/'RO_Xs0, X_s0, INERT_SPECIES,',&
         ' or USR_ROS.')

         ENDIF

! Check physical restrictions on variable solids density input.
         IF(SOLVE_ROs(M)) THEN
! Check INERT_SPECIES
            IF(INERT_SPECIES(M) == UNDEFINED_I) THEN
               WRITE(ERR_MSG,1000) trim(iVar('INERT_SPECIES',M))
               CALL LOG_ERROR()
            ELSEIF(INERT_SPECIES(M) < 1 .OR.                          &
               INERT_SPECIES(M) > NMAX(M)) THEN
               WRITE(ERR_MSG, 1001) trim(iVar('INERT_SPECIES',M)),    &
                  trim(iVal(INERT_SPECIES(M)))
               CALL LOG_ERROR()
            ENDIF
            IF(DIL_INERT_X_VSD(M)<=ZERO) THEN
               WRITE(ERR_MSG,1103) M, DIL_INERT_X_VSD(M)
               CALL LOG_ERROR()
            ELSEIF(DIL_INERT_X_VSD(M)>ONE) THEN
               WRITE(ERR_MSG,1104) M, DIL_INERT_X_VSD(M)
               CALL LOG_ERROR()
            ENDIF

            sum_xs = zero
            DO N=1, NMAX(M)
! Check RO_Xs0
               IF(RO_Xs0(M,N) == UNDEFINED) THEN
                  WRITE(ERR_MSG,1000) trim(iVar('RO_Xs0',M,N))
                  CALL LOG_ERROR()
               ELSEIF(RO_Xs0(M,N) < ZERO) THEN
                  WRITE(ERR_MSG,1001) trim(iVar('RO_Xs0',M,N)),        &
                     trim(iVal(RO_xs0(M,N)))
                  CALL LOG_ERROR()
               ENDIF
! Check X_s0
               IF(X_s0(M,N) == UNDEFINED) THEN
                  WRITE(ERR_MSG,1000) trim(iVar('X_s0',M,N))
                  CALL LOG_ERROR()
               ELSEIF(X_s0(M,N) < ZERO .OR. X_s0(M,N) > ONE) THEN
                  WRITE(ERR_MSG,1001) trim(iVar('X_s0',M,N)),        &
                     trim(iVal(X_s0(M,N)))
                  CALL LOG_ERROR()
               ENDIF
! Sum of species mass fraction
               sum_xs = sum_xs + x_s0(M,N)
            ENDDO

! Enforce that the species mass fractions must sum to one.
            IF(.NOT.COMPARE(ONE,SUM_xs)) THEN
               WRITE(ERR_MSG, 1110) M
               CALL LOG_ERROR()
 1110 FORMAT('Error 1110: Xs_0(',I3,',:) do not sum to 1 and ',&
         'variable',/'solids density is requested.')
            ENDIF

! Check for input overflow.
            DO N=NMAX(M)+1, DIM_N_s
               IF(RO_Xs0(M,N) /= UNDEFINED) THEN
                  WRITE(ERR_MSG,1002) trim(iVar('RO_Xs0',M,N))
                  CALL LOG_ERROR()
               ENDIF
               IF(X_s0(M,N) /= UNDEFINED) THEN
                  WRITE(ERR_MSG,1002) trim(iVar('X_s0',M,N))
                  CALL LOG_ERROR()
               ENDIF
            ENDDO

! Check X_s0(Inert)
            IF(X_s0(M,INERT_SPECIES(M)) == ZERO) THEN
               WRITE(ERR_MSG,1102) M, INERT_SPECIES(M)
               CALL LOG_ERROR()
            ENDIF

 1102 FORMAT('Error 1102: Invalid baseline inert species mass',/       &
         'fraction. The inert species mass fraction must be greater ',  &
         'than zero.',/' Phase ',I2,' Inert Species: ',I3,' X_s0 = 0.0')

 1103 FORMAT('Error 1103: Invalid dilute region inert species mass',/   &
         'fraction. The inert species mass fraction must be greater ',  &
         'than zero.',/' Phase ',I2,/ &
         ' Please check the value of DIL_INERT_X_VSD:',G14.4)

 1104 FORMAT('Error 1104: Invalid dilute region inert species mass',/   &
         'fraction. The inert species mass fraction must be less ',  &
         'than or equal to one.',/' Phase ',I2,/ &
         ' Please check the value of DIL_INERT_X_VSD:',G14.4)

! All of the information for variable solids density has been verified
! as of this point. Calculate and store the baseline density.
            RO_s0(M) = EOSS0(M)

! Check physical restrictions on constant solids density input.
         ELSEIF(RO_S0(M) <= ZERO) THEN
            WRITE(ERR_MSG,1101) iVar('RO_s0',M), iVal(RO_s0(M))
            CALL LOG_ERROR()
         ENDIF

      ENDDO

! Check for input overflow.
      DO M = MMAX_LL+1, DIM_M
         IF(RO_S0(M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1002) trim(iVar('RO_s0',M))
            CALL LOG_ERROR()
         ENDIF
         IF(INERT_SPECIES(M) /= UNDEFINED_I) THEN
            WRITE(ERR_MSG,1002) trim(iVar('INERT_SPECIES',M))
            CALL LOG_ERROR()
         ENDIF
         DO N=1, DIM_N_s
            IF(RO_Xs0(M,N) /= UNDEFINED) THEN
               WRITE(ERR_MSG,1002) trim(iVar('RO_Xs0',M,N))
               CALL LOG_ERROR()
            ENDIF
            IF(X_s0(M,N) /= UNDEFINED) THEN
               WRITE(ERR_MSG,1002) trim(iVar('X_s0',M,N))
               CALL LOG_ERROR()
            ENDIF
         ENDDO
      ENDDO

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A)

 1001 FORMAT('Error 1001: Invalid or unphysical input: ',A,' = ',A)

 1002 FORMAT('Error 1002: Invalid input: ',A,' specified out of range.')

      END SUBROUTINE CHECK_SOLIDS_DENSITY

END MODULE CHECK_SOLIDS_COMMON_ALL_MOD
