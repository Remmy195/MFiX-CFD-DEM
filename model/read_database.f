#include "error.inc"

MODULE READ_DATABASE_MOD

  USE discretelement, only: discrete_element, des_continuum_hybrid
  USE error_manager
  USE funits, only: unit_dat, unit_log, dmp_log, file_size
  Use param, only: dim_m, dim_n
  USE param1, only: zero, undefined, small_number, undefined_c
  USE physprop, only: t_ref, icpor_h, icpor_l, tcom, ahigh, alow, thigh, tlow, hfrefor, tmiddle, amiddle,icpor_m
  USE physprop, only: database_read, nmax, mmax, mw_g, mw_s
  USE read_thermochemical, only: read_therm, calc_ICpoR, calc_CpoR
  use resize, only: real_grow3
  USE run, only: REINITIALIZING, mfix_input_filename
  USE rxns, only: species_name

  CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: read_database(Ier)                                     C
!  Purpose: read thermochemical database                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-OCT-05  C
!                                                                      C
!  Modification 1: J. Musser                          Date: 02-May-11  C
!  Purpose: Provided support for DEM access to database.               C
!                                                                      C
!  Modification 2: J. Musser                          Date: 02-Oct-12  C
!  Purpose: Calls to READ_DATABASE were moved to check_gas_phase and   C
!  check_solids_common_all during input data integrity checks for the  C
!  gas and solids phases.                                              C
!  Rather than looping through all species for each phase, the model   C
!  (TFM/DEM), phase index, species index, species name, and molecular  C
!  weight are passed as dummy arguments so that only information for   C
!  referenced species (lName) is obtained.                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE READ_DATABASE(MFIX_DAT, lM, lN, lName, lMW)

      IMPLICIT NONE

! project filename
      CHARACTER(LEN=*), INTENT(IN) :: MFIX_DAT

! Phase and species indices
      INTEGER, INTENT(IN) :: lM, lN
! Species name from data file.
      CHARACTER(len=*), INTENT(IN) :: lName
! Species molecular weight from the data file (if any)
      DOUBLE PRECISION, INTENT(INOUT) :: lMW

! Molecular weight read from database.
      DOUBLE PRECISION dbMW

! Error message returned from Read_Therm and sent to calling routine
      INTEGER IER

! Identifies if the species was found in any of the databases.
      LOGICAL SPECIES_FOUND

! Tcom +/- SMALL_NUMBER: This is done so that the specific heat
! polynomial can be evaluated at Tcom with the high and low
! coefficients.
      DOUBLE PRECISION :: xTc
! Tmiddle +/- SMALL_NUMBER: This is done so that the specific heat
! polynomial can be evaluated at Tmiddle with coefficients for the
! temperature ranges lower and higher than Tmiddle
      DOUBLE PRECISION :: xTm

! Various integrations of the specific heat polynomials:
      DOUBLE PRECISION :: ICpoR_TrL  ! 0.0 --> Tref using Alow
      DOUBLE PRECISION :: ICpoR_TcL  ! 0.0 --> Tcom using Alow
      DOUBLE PRECISION :: ICpoR_TcH  ! 0.0 --> Tcom using Ahigh
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ICpoR_TmL  ! 0.0 --> Tm using Amiddle-1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ICpoR_TmH  ! 0.0 --> Tm using Amiddle

      LOGICAL :: testCp = .FALSE.

      ! Full path to Burcat and Ruscic database
      CHARACTER(len=142) :: THERM = 'BURCAT.THR'
      INTEGER :: THERMO_SIZE
      INTEGER :: old_size, new_size, i

#ifdef BURCAT_THR
      THERM = BURCAT_THR
#endif


! Check for thermochemical data in input file.
      CALL READ_A_DATABASE(mfix_input_filename, MFIX_DAT)

      IF(.NOT.SPECIES_FOUND) THEN
! Check for thermochemical data from the BURCAT.THR database in the local
! run directory.
         THERMO_SIZE = FILE_SIZE(THERM)
         CALL READ_THERMO_DATABASE(THERMO_SIZE)
      ENDIF

! Error message control.
!-----------------------------------------------------------------------
! Write remaining error message if needed.
      IF(.NOT.SPECIES_FOUND) THEN
#ifdef PYMFIX
         WRITE(ERR_MSG,'(A,A)') 'Unable to find species: ', trim(lName)
#else
         WRITE(ERR_MSG,1010) trim(lName), trim(THERM)
#endif
         CALL LOG_ERROR()
      ENDIF

      RETURN

! Error Flags
!-----------------------------------------------------------------------
 1010 FORMAT('Message 1010: Species "',A,'" was not matched to any ',  &
         'entry in the',/'thermochemical databases.',2/,'SUGGESTION: ',&
         'Search the database for the exact species name. The ',/      &
         'species names are case sensitive and should match the names',&
         ' in',/'BURCAT.THR exactly excluding trailing blanks and ',   &
         'tabs. Also verify',/'that the data section in the project ', &
         'file (if any) is below a line',/'that starts with THERMO ',  &
         'DATA.',2/'Database location:', /A)

      CONTAINS

! Avoid using Fortran 2003 deferred-length string by passing in size to a subroutine
      SUBROUTINE READ_THERMO_DATABASE(THERM_SIZE)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: THERM_SIZE
        CHARACTER(THERM_SIZE) :: THERM_DB

        OPEN(UNIT=UNIT_DAT, FILE=TRIM(THERM), ACTION='read', FORM="unformatted", ACCESS="stream")
        READ(UNIT_DAT) THERM_DB
        CLOSE(UNIT_DAT)

        CALL READ_A_DATABASE(THERM, THERM_DB)

      END SUBROUTINE READ_THERMO_DATABASE


      SUBROUTINE READ_A_DATABASE(THERMO_DB_LABEL, THERMO_DB)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN) :: THERMO_DB_LABEL
         CHARACTER(LEN=*), INTENT(IN) :: THERMO_DB

         IF(.NOT. ALLOCATED(Tmiddle)) ALLOCATE(Tmiddle(0:DIM_M, DIM_N, 1))
         IF(.NOT. ALLOCATED(Amiddle)) ALLOCATE(Amiddle(7, 0:DIM_M, DIM_N, 1))
         IF(.NOT. ALLOCATED(ICpoR_m)) ALLOCATE(ICpoR_m(0:DIM_M, DIM_N, 1))
! Initialize the coefficient for different temperature ranges as zero
         Tmiddle(lM,lN,:) = UNDEFINED
         Amiddle(:,lM,lN,:) = UNDEFINED

         IER = 0
! First, try to find the species in a Burcat block
         CALL READ_THERM('BURCAT', &
            THERMO_DB, lName, Thigh(lM,lN), Tlow(lM,lN),      &
            Tcom(lM,lN), dbMW, Ahigh(:,lM,lN), Alow(:,lM,lN), &
            HfrefoR(lM,lN), IER)

! If the species was not found in a Burcat block, try to find it in a
! Chemkin block.
         IF(IER/=0) THEN
            CALL READ_THERM('CHEMKIN', &
               THERMO_DB, lName, Thigh(lM,lN), Tlow(lM,lN),      &
               Tcom(lM,lN), dbMW, Ahigh(:,lM,lN), Alow(:,lM,lN), &
               HfrefoR(lM,lN), IER, lM, lN, Tmiddle, Amiddle)
         ENDIF

         IF(IER == 0) THEN
! If the user did not supply a value for the gas phase molecular weight,
! use the value from the database.
            IF(lMW == UNDEFINED) THEN
               lMW = dbMW
               ! WRITE(ERR_MSG,2000) trim(lName),trim(lName),lMW
               ! CALL LOG_INFO()
            ENDIF
! There are a number of species with Tlow as 300, for which the
! following calculation will produce an error because T_ref = 298.  So
! slightly extend validity of the correlation.
            IF(ABS(Tlow(lM,lN)-T_ref)<=2.0D0 .AND. &
               Tlow(lM,lN) > T_ref) Tlow(lM,lN) = T_ref

! Some species (such as H2O_L) have Tcom<1000K and have zero Ahigh coefficients.
! Their Thigh was set to Tcom in read_thermochemical_mod.f so we only check when Thigh>Tcom
            IF(Thigh(lM,lN)>Tcom(lM,lN)) THEN
               IF(MAXVAL(DABS(Ahigh(:,lM,lN)))==ZERO) THEN
                  WRITE(ERR_MSG,1010) 'high', trim(lName)
                  CALL LOG_ERROR()
               ENDIF
            ENDIF

            IF(MAXVAL(DABS(Alow(:,lM,lN)))==ZERO) THEN
               WRITE(ERR_MSG,1010) 'low', trim(lName)
               CALL LOG_ERROR()
            ENDIF

! Initialize the reference integrals.
            ICpoR_l(lM,lN) = ZERO
            ICpoR_h(lM,lN) = ZERO
            ICpoR_m(lM,lN,:) = ZERO

! Calculate the integral of specific heat from zero to Tref using the
! Alow coefficients. Store the integrals in global variables.
            ICpoR_TrL = calc_ICpoR(T_ref, lM, lN, IER)
! Only two temperature ranges are given
            IF(Tmiddle(lM, lN, 1) .EQ. UNDEFINED) THEN
! Calculate the integral of specific heat from zero to Tcom using the
! Alow coefficients.
                xTc = Tcom(lM,lN)-Tcom(lM,lN)*1.0d-6
                ICpoR_TcL = calc_ICpoR(xTc, lM, lN, IER)
! Calculate the integral of specific heat from zero to Tcom using the
! Ahigh coefficients.
                IF(abs(Tcom(lM,lN) - Thigh(lM,lN)) > epsilon(0.0)) THEN
                    xTc = Tcom(lM,lN)+Tcom(lM,lN)*1.0d-6
                    ICpoR_TcH = calc_ICpoR(xTc, lM, lN, IER)
                ELSE
                    ICpoR_TcH = 0.0d0
                ENDIF
! Store the integrals in global variables.
                ICpoR_h(lM,lN) = ICpoR_TcH - (ICpoR_TcL - ICpoR_TrL)
! More than two temperature ranges are given
            ELSE
                old_size = SIZE(ICpoR_m, 3)
                new_size = SIZE(Tmiddle, 3)
                IF(new_size .GT. (old_size+1)) THEN
                    CALL REAL_GROW3(ICpoR_m, new_size-1)
                    ICpoR_m(:,:,old_size+1:) = ZERO
                ENDIF
                ALLOCATE(ICpoR_TmL(new_size))
                ALLOCATE(ICpoR_TmH(new_size))
                DO i=1, SIZE(Tmiddle, 3)
                   IF(Tmiddle(lM, lN, i) .NE. UNDEFINED) THEN
                       xTm = Tmiddle(lM,lN,i)-Tmiddle(lM,lN,i)*1.0d-6
                       ICpoR_TmL(i) = calc_ICpoR(xTm, lM, lN, IER)
                       xTm = Tmiddle(lM,lN,i)+Tmiddle(lM,lN,i)*1.0d-6
                       ICpoR_TmH(i) = calc_ICpoR(xTm, lM, lN, IER)

                   ENDIF
                ENDDO
!! Need to assign ICpoR_l, ICpoR_m, ICpoR_h after calculating all ICpoR_TrL, ICpoR_TmL and ICpoR_TmH
                DO i=1, SIZE(Tmiddle, 3)
                   IF(i .EQ. 1) THEN
                       ICpoR_m(lM,lN,i) = ICpoR_TmH(i) - (ICpoR_TmL(i) - ICpoR_TrL)
                   ELSEIF(i .LT. SIZE(Tmiddle, 3)) THEN
                       ICpoR_m(lM,lN,i) = ICpoR_TmH(i) - ICpoR_TmL(i) + ICpoR_m(lM,lN,i-1)
                   ELSE
                       ICpoR_h(lM,lN) = ICpoR_TmH(i) - ICpoR_TmL(i) + ICpoR_m(lM,lN,i-1)
                   ENDIF
                ENDDO
                DEALLOCATE(ICpoR_TmL)
                DEALLOCATE(ICpoR_TmH)
            ENDIF
            ICpoR_l(lM,lN) = ICpoR_TrL
         ENDIF

         SPECIES_FOUND = (IER == 0)
         IF(.NOT.SPECIES_FOUND) THEN
            IF(.NOT.REINITIALIZING)THEN
               WRITE(ERR_MSG,1001) trim(adjustl(THERMO_DB_LABEL)), 'Not found.'
               CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)
            ENDIF

         ELSE
            IF(.NOT.REINITIALIZING)THEN
               WRITE(ERR_MSG,1001) "Checking for ", trim(lName), " in ", trim(adjustl(THERMO_DB_LABEL)), 'Found.'
               CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)
            ENDIF

            if(testCP) CALL writeCp(lM, lN, lName, lMW)

! Check Cp(T) is defined for all temperatures
            CALL Check_Cp_vs_T(lM, lN, lName, lMW)

            CALL Echo_database_entry(lM, lN, lName, lMW)
         ENDIF

! Messages
!-----------------------------------------------------------------------
 1001 FORMAT(8X,A,1X,A,A,A,' : ',A)
 1010 FORMAT('Error 1010: All Cp coefficients for ',A, ' temperature range', &
       /,'are zero for species: ',A)
 2000 FORMAT('Info: Molecular weight keyword was not set for species: ',A, &
            /'Its value was read from the database instead.', &
            /A,' molecular weight =',G12.5)

      END SUBROUTINE READ_A_DATABASE

      END SUBROUTINE READ_DATABASE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: READ_DATABASE0(IER)                                    C
!  Purpose: Provides legacy support for rrates files.                  C
!                                                                      C
!  Author: J. Musser                                  Date: 02-Oct-12  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE READ_DATABASE0(MFIX_DAT)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: MFIX_DAT

! Loop indices for mass phase and species
      INTEGER M, NN
! Loop counter for continuum and discrete species
      INTEGER Nsp, DES_Nsp

! Return to the calling routine if the database has already been called.
      IF(database_read)RETURN

! Set the flag identifying that the database has been read.
      database_read = .TRUE.

! Initialize counters
      Nsp = 0
      DES_Nsp = 0

! Read species data for the gas phase.
!-----------------------------------------------------------------------
      DO NN = 1, NMAX(0)
         Nsp = Nsp + 1
! If a species name was not specified, flag error and exit.
          IF(SPECIES_NAME(Nsp) == UNDEFINED_C) THEN
            WRITE(*,1010) NN         ! screen
            IF(DMP_LOG) WRITE(UNIT_LOG,1010) NN  ! log file
             call log_error()
          ENDIF
! Read the database.
         CALL READ_DATABASE(MFIX_DAT, 0, NN, SPECIES_NAME(Nsp), MW_g(NN))
       ENDDO

! Read species data for the continuum solids phases.
!-----------------------------------------------------------------------
! Skip reading the database for the continuum solids phase if the
! simulation is only employing discrete solids.
      IF(.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID)THEN
          DO M = 1, MMAX
            DO NN = 1, NMAX(M)
                Nsp = Nsp + 1
! If a species name was not specified, flag error and exit.
                IF(SPECIES_NAME(Nsp) == UNDEFINED_C)THEN
                  WRITE(*,1011)'continuum', M, NN ! screen
                  IF(DMP_LOG) WRITE(UNIT_LOG,1011)'continuum', M, NN
                   call log_error()
                ENDIF
               CALL READ_DATABASE(MFIX_DAT, M, NN, SPECIES_NAME(Nsp), MW_s(M,NN))
             ENDDO   ! N=1, NMAX(M)
          ENDDO   ! M=1, MMAX
      ENDIF

      RETURN

! Error Messages
!-----------------------------------------------------------------------
 1010 FORMAT(/1X,70('*')/, ' From: READ_DATABASE0',/, ' Message: ',    &
         'No SPECIES_NAME provided for gas phase species ',I3,'.',/' ',&
         /1X,70('*')/)
 1011 FORMAT(/1X,70('*')/, ' From: READ_DATABASE0',/, ' Message: ',    &
         'No SPECIES_NAME provided for ',A,' solids phase ',I2,', ',/  &
         ' species ',I3,'.',/1X,70('*')/)

      END SUBROUTINE READ_DATABASE0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: read_database(Ier)                                     C
!  Purpose: read thermochemical database                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-OCT-05  C
!                                                                      C
!  Modification 1: J. Musser                          Date: 02-May-11  C
!  Purpose: Provided support for DEM access to database.               C
!                                                                      C
!  Modification 2: J. Musser                          Date: 02-Oct-12  C
!  Purpose: Calls to READ_DATABASE were moved to CHECK_DATA_04/05      C
!  duing input data integrity checks for the gas and solids phases.    C
!  Rather than looping through all species for each phase, the model   C
!  (TFM/DEM), phase index, species index, species name, and molecular  C
!  weight are passed as dummy arguments so that only information for   C
!  referenced species (lName) is obtained.                             C                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE writeCp(lM, lN, lName, lMW)

! Universal gas constant in cal/mol.K
      use constant, only: RGAS => GAS_CONST_cal

      IMPLICIT NONE

! Phase and species indices
      INTEGER, INTENT(IN) :: lM, lN
! Species name from data file.
      CHARACTER(len=*), INTENT(IN) :: lName
! Species molecular weight from the data file (if any)
      DOUBLE PRECISION, INTENT(in) :: lMW

      INTEGER :: IER1, IER2, lc

      DOUBLE PRECISION :: T
      DOUBLE PRECISION :: lCP, lICP

      CHARACTER(len=40) :: fmt
      INTEGER num_middleT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: T_compare

      write(*,"(2/,3x,'Specific Heat report for ',A)")trim(lName)

      IF(Tmiddle(lM, lN, 1) .EQ. UNDEFINED) THEN
         write(*,"(/12x,'Low',12x,'High')")
         write(*,"(6x,'T',3x,g12.5,2x,g12.5)") Tlow(lM,lN), Thigh(lM,lN)
         DO lc=1,5
            write(*,"(4x,'A(',I1,')',2(2x,g12.5))") lc, &
                 Alow(lc,lM,lN), Ahigh(lc,lM,lN)
         ENDDO
      ELSE
         DO lc=1,SIZE(Tmiddle(lM,lN,:))
            IF(Tmiddle(lM,lN,lc) .NE. UNDEFINED) num_middleT = lc
         ENDDO
         write(fmt,'(A, I3, A)') "(/12x,'Low',15x,'Middle',",15*num_middleT,"x,'High')"
         write(*,fmt)
         write(fmt,'(A, I3, A)') "(6x,'T',3x,",lc+2,'(g12.5,6x))'
         write(*,fmt) Tlow(lM,lN), Tmiddle(lM,lN,:), Thigh(lM,lN)
         DO lc=1,5
            write(fmt,'(A, I3, A)') "(5x,'A(',I1,')',",num_middleT+1,"(8x,g12.5))"
            write(*,fmt) lc, Alow(lc,lM,lN), Amiddle(lc,lM,lN,:), Ahigh(lc,lM,lN)
         ENDDO
      ENDIF

      write(*,"('')")
      IF(Tmiddle(lM, lN, 1) .EQ. UNDEFINED) write(*,"(5x,'Tcom: ',g12.5)")Tcom(lM,lN)
      write(*,"('')")

      write(*,"(5x,'Temperature',8x,'Cp',11x,'ICp')")
      IF(Tmiddle(lM, lN, 1) .EQ. UNDEFINED) THEN
         ALLOCATE(T_compare(1))
         T_compare = Tcom(lM,lN)
      ELSE
         ALLOCATE(T_compare(num_middleT))
         T_compare = Tmiddle(lM,lN,:)
      ENDIF

      DO lc=1,SIZE(T_compare)
         T = T_compare(lc) - 100.0
         DO WHILE(T <= T_compare(lc) - SMALL_NUMBER)

            IER1 = 0
            IER2 = 0

            write(*,"(7x,g12.5)",ADVANCE="NO") T
            lCP  = calc_CpoR(T, lM, lN) * RGAS / lMW
            lICP = calc_ICpoR(T, lM, lN, IER2) * RGAS / lMW
            write(*,"(2(3x,g12.5))",ADVANCE="NO")lCP, lICP

            IF(IER1 /= 0) write(*,"(3x,'Cp Error!')",ADVANCE="NO")
            IF(IER2 /= 0) write(*,"(3x,'ICp Error!')",ADVANCE="NO")
            write(*,"('')")

            T = T + 5.0
         ENDDO

         T = T_compare(lc) + SMALL_NUMBER
         DO WHILE(T <= T_compare(lc) + 100.0)

            IER1 = 0
            IER2 = 0

            write(*,"(7x,g12.5)",ADVANCE="NO") T
            lCP  = calc_CpoR(T, lM, lN) * RGAS / lMW
            lICP = calc_ICpoR(T, lM, lN, IER2) * RGAS / lMW
            write(*,"(2(3x,g12.5))",ADVANCE="NO")lCP, lICP

            IF(IER1 /= 0) write(*,"(3x,'Cp Error!')",ADVANCE="NO")
            IF(IER2 /= 0) write(*,"(3x,'ICp Error!')",ADVANCE="NO")
            write(*,"('')")

            T = T + 5.0
         ENDDO

         write(*,"('')")
      ENDDO

      write(*,"('')")

      DEALLOCATE(T_compare)
      END SUBROUTINE writeCp

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Check_Cp_vs_T                                          C
!  Purpose: Verify Cp is defined for all temperatures                  C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 22-APR-22  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE Check_Cp_vs_T(lM, lN, lName, lMW)

! Universal gas constant in cal/mol.K
      use constant, only: GAS_CONST_cal
      use param1, only: ZERO
      use toleranc, only: TMIN, TMAX
      use run, only: UNITS
      use read_thermochemical, only: calc_CpoR0

      IMPLICIT NONE

! Phase and species indices
      INTEGER, INTENT(IN) :: lM, lN
! Species name from data file.
      CHARACTER(len=*), INTENT(IN) :: lName
! Species molecular weight from the data file (if any)
      DOUBLE PRECISION, INTENT(in) :: lMW ! g/mol in cgs, kg/kmol in mks (SI)

      INTEGER :: IER1, IER2, lc, j, num_middleT

      DOUBLE PRECISION :: T
      INTEGER :: TLISTSIZE,I,NT
      DOUBLE PRECISION, ALLOCATABLE :: TLIST(:)  ! List of temperatures where Cp(T) is computed
      DOUBLE PRECISION, PARAMETER :: DT = 10.0D0 ! Temperature increment in list
      DOUBLE PRECISION :: lCP, lICP
      DOUBLE PRECISION :: GAS_CONSTANT
      CHARACTER(LEN=32):: CP_UNIT
      DOUBLE PRECISION :: TRANGE

      CHARACTER(len=100) :: fmt

! Setup units for Cp in case we need to display an error
      IF(UNITS == 'SI') THEN
         GAS_CONSTANT = 4.183925d3 * GAS_CONST_cal !(J/kmol.K)
         CP_UNIT      = 'J/(kg.K)' ! Universal gas constant divided by molecular weight (kg/kmol)
      ELSE
         GAS_CONSTANT = GAS_CONST_cal ! (cal/mol.K)
         CP_UNIT      = 'cal/(g.K)' ! Universal gas constant divided by molecular weight (g/mol)
      ENDIF

! Build a list of temperatures at which we compute Cp
! The list spans roughly from TMIN to TMAX with "nice" values with increment DT
! (default is DT = 10 Kelvins)
! Then we add known temperatures (TMIN, TMAX, Tlow, Tcom, Thigh)

      TLISTSIZE = INT((TMAX - TMIN) / DT) + 10 + SIZE(Tmiddle(lM,lN,:)) ! Set rough size of list + few extra slots

      allocate( TLIST  (TLISTSIZE) )

      TLIST(:) = ZERO

! Nicely incremented values.
      I = 0
      T = ZERO
      DO WHILE(T <= TMAX)
         IF(T >= TMIN) THEN
            I = I + 1
            TLIST(I) = T
         ENDIF
         T = T + DT
      ENDDO

! Add known "important" temperatures.There could be duplicates with
! the above if these temperatures are multiples of DT
      TLIST(I + 1) = Tlow(lM,lN)
      TLIST(I + 2) = TMIN
      IF(Tmiddle(lM,lN,1) .EQ. UNDEFINED) THEN
         TLIST(I + 3) = Tcom(lM,lN)
         TLIST(I + 4) = TMAX
         TLIST(I + 5) = Thigh(lM,lN)
! Keep last index in the list so we only loop through defined temperatures below.
         NT = I + 5
      ELSE
         num_middleT = 0
         DO j=1,SIZE(Tmiddle(lM,lN,:))
            IF(Tmiddle(lM,lN,j) .NE. UNDEFINED) THEN
               TLIST(I+3+j) = Tmiddle(lM,lN,j)
               num_middleT = num_middleT + 1
            ELSE
               exit
            ENDIF
         ENDDO
         TLIST(I+num_middleT+3) = TMAX
         TLIST(I+num_middleT+4) = Thigh(lM,lN)
! Keep last index in the list so we only loop through defined temperatures below.
         NT = I+num_middleT+4
      ENDIF

! Now go through the list of temperatures, compute Cp and verify it is not zero
! nor negative
      DO I = 1,NT

         T = TLIST(I)
         lCP  = calc_CpoR(T, lM, lN) * GAS_CONSTANT / lMW   !J/(kg.K) or cal/(g.K)

         IF(lCP<=ZERO) THEN
            IF(Tmiddle(lM,lN,1) .EQ. UNDEFINED) THEN
               WRITE(ERR_MSG,1020) trim(lName), &
                                   lM,lN, &
                                   T,lCP,trim(CP_UNIT), &
                                   TMIN, &
                                   TMAX, &
                                   Tlow(lM,lN), Tcom(lM,lN), &
                                   Tcom(lM,lN), Thigh(lM,lN), &
                                   Alow(1,lM,lN), Ahigh(1,lM,lN) , &
                                   Alow(2,lM,lN), Ahigh(2,lM,lN) , &
                                   Alow(3,lM,lN), Ahigh(3,lM,lN) , &
                                   Alow(4,lM,lN), Ahigh(4,lM,lN) , &
                                   Alow(5,lM,lN), Ahigh(5,lM,lN)
            ELSE
               WRITE(*,1030) trim(lName), &
                                   lM,lN, &
                                   T,lCP,trim(CP_UNIT), &
                                   TMIN, &
                                   TMAX, &
                                   Tlow(lM,lN)
               WRITE(fmt,'(A, I3, A)') "('Middle temperautres = ',", num_middleT,"(g12.5), 'K')"
               WRITE(*,fmt) Tmiddle(lM,lN,:)
               WRITE(fmt,'(A)') "('High temperature = ',g12.5,'K', /'        Low range           Middle range           High range')"
               WRITE(*,fmt) Thigh(lM,lN)
               DO j=1,4
                  WRITE(fmt,'(A, I3, A, I3, A)') "(' A(', I1, ')  ',g12.5, 6x,", num_middleT,"(g12.5, 6x), g12.5, 'K')"
                  WRITE(*,fmt) j, Alow(j,lM,lN), Amiddle(j,lM,lN,:), Ahigh(j,lM,lN)
               ENDDO
               WRITE(fmt,'(A, I3, A)') "(' A(5)  ', g12.5, 6x,", num_middleT,"(g12.5, 6x), g12.5, 'K')"
               WRITE(ERR_MSG,fmt) Alow(5,lM,lN), Amiddle(5,lM,lN,:), Ahigh(5,lM,lN)
            ENDIF

            CALL LOG_ERROR()  !Cp<=zero is a fatal error
         ENDIF

      ENDDO

! If there was no error, we are good to go
      WRITE(ERR_MSG,1010)
      CALL LOG_STATUS()

      deallocate( TLIST )

      RETURN


1010 FORMAT(/8x,'Checking Specific Heat is defined along entire temperature range: Yes')

1020 FORMAT(/8x,'Checking Specific Heat is defined along entire temperature range: No', &
           /'Error: Zero of negative Cp for Species: ',A      , &
           /'       Phase index :',I2,', species index:',I2   , &
           /'       Temperature = ',g12.5,', Cp = ',g12.5,A   , &
           /'Please verify the Cp coefficients and temperature range.', &
           /'Minimum allowed temperature (TMIN) = ',g12.5,'K', &
           /'Maximum allowed temperature (TMAX) = ',g12.5,'K', &
           /'Low  temperature range (Tlow to Tcom)  = ',g12.5,'K to',g12.5,'K', &
           /'High temperature range (Tcom to Thigh) = ',g12.5,'K to',g12.5,'K', &
           /'        Low range     High range',&
           /' A(1)  ', g12.5, 2x,g12.5, &
           /' A(2)  ', g12.5, 2x,g12.5, &
           /' A(3)  ', g12.5, 2x,g12.5, &
           /' A(4)  ', g12.5, 2x,g12.5, &
           /' A(5)  ', g12.5, 2x,g12.5 &
     )

1030 FORMAT(/8x,'Checking Specific Heat is defined along entire temperature range: No', &
           /'Error: Zero of negative Cp for Species: ',A      , &
           /'       Phase index :',I2,', species index:',I2   , &
           /'       Temperature = ',g12.5,', Cp = ',g12.5,A   , &
           /'Please verify the Cp coefficients and temperature range.', &
           /'Minimum allowed temperature (TMIN) = ',g12.5,'K', &
           /'Maximum allowed temperature (TMAX) = ',g12.5,'K', &
           /'Low  temperature  = ',g12.5,'K'&
     )


      END SUBROUTINE Check_Cp_vs_T

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ECHO_DATABASE_ENTRY                                    C
!  Purpose: Print a summary of the database entry                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 22-MAR-23  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE Echo_database_entry(lM, lN, lName, lMW)

! Universal gas constant in cal/mol.K
      use constant, only: GAS_CONST_cal
      use param1, only: ZERO
      use toleranc, only: TMIN, TMAX
      use run, only: UNITS

      IMPLICIT NONE

! Phase and species indices
      INTEGER, INTENT(IN) :: lM, lN
! Species name from data file.
      CHARACTER(len=*), INTENT(IN) :: lName
! Species molecular weight from the data file (if any)
      DOUBLE PRECISION, INTENT(in) :: lMW ! g/mol in cgs, kg/kmol in mks (SI)

      INTEGER :: IER1, IER2, lc

      DOUBLE PRECISION :: GAS_CONSTANT
      CHARACTER(LEN=32):: CP_UNIT
      DOUBLE PRECISION :: TRANGE

      CHARACTER(len=150) :: fmt
      INTEGER :: num_middleT, j

! Setup units for Cp in case we need to display an error
      IF(UNITS == 'SI') THEN
         GAS_CONSTANT = 4.183925d3 * GAS_CONST_cal !(J/kmol.K)
         CP_UNIT      = 'J/(kg.K)' ! Universal gas constant divided by molecular weight (kg/kmol)
      ELSE
         GAS_CONSTANT = GAS_CONST_cal ! (cal/mol.K)
         CP_UNIT      = 'cal/(g.K)' ! Universal gas constant divided by molecular weight (g/mol)
      ENDIF

      IF(Tmiddle(lM,lN,1) .EQ. UNDEFINED) THEN
         WRITE(ERR_MSG,1000) trim(lName), &
                             lM,lN, &
                             TMIN, &
                             TMAX, &
                             Tlow(lM,lN), Tcom(lM,lN), &
                             Tcom(lM,lN), Thigh(lM,lN), &
                             Alow(1,lM,lN), Ahigh(1,lM,lN) , &
                             Alow(2,lM,lN), Ahigh(2,lM,lN) , &
                             Alow(3,lM,lN), Ahigh(3,lM,lN) , &
                             Alow(4,lM,lN), Ahigh(4,lM,lN) , &
                             Alow(5,lM,lN), Ahigh(5,lM,lN) , &
                             Alow(6,lM,lN), Ahigh(6,lM,lN) , &
                             Alow(7,lM,lN), Ahigh(7,lM,lN) , &
                             lMW, HfrefoR(lM,lN)
      ELSE
         num_middleT = 0
         DO j=1,SIZE(Tmiddle(lM,lN,:))
            IF(Tmiddle(lM,lN,j) .NE. UNDEFINED) THEN
               num_middleT = num_middleT + 1
            ELSE
               exit
            ENDIF
         ENDDO
         WRITE(*,1010) trim(lName), &
                             lM,lN, &
                             TMIN, &
                             TMAX, &
                             Tlow(lM,lN)

         WRITE(fmt,'(A, I3, A)') "('Middle temperautres = ',", num_middleT,"(g12.5), 'K')"
         WRITE(*,fmt) Tmiddle(lM,lN,:)
         WRITE(fmt,'(A)') "('High temperature = ',g12.5,'K', /'Polynomial coefficients:', /'        Low range           Middle range           High range')"
         WRITE(*,fmt) Thigh(lM,lN)
         DO j=1,7
            WRITE(fmt,'(A, I3, A, I3, A)') "(' A(', I1, ')  ',g12.5, 6x,", num_middleT,"(g12.5, 6x), g12.5, 'K')"
            WRITE(*,fmt) j, Alow(j,lM,lN), Amiddle(j,lM,lN,:), Ahigh(j,lM,lN)
         ENDDO
         WRITE(fmt,'(A, I3, A)') "('Molecular weight = ', g12.5)"
         WRITE(*,fmt) lMW
         WRITE(fmt,'(A, I3, A)') "('Heat of formation/R = ', g12.5)"
         WRITE(ERR_MSG,fmt) HfrefoR(lM,lN)
      ENDIF


      CALL log_message(__FILE__, __LINE__, LOGLEVEL_INFO, header=.false., footer=.false. )


      RETURN



1000 FORMAT(/8x,'Summary of thermochemical data read from the database:', &
           /'Species: ',A      , &
           /'Phase index :',I2,', species index:',I2   , &
           /'Minimum allowed temperature (TMIN) = ',g12.5,'K', &
           /'Maximum allowed temperature (TMAX) = ',g12.5,'K', &
           /'Low  temperature range (Tlow to Tcom)  = ',g12.5,'K to',g12.5,'K', &
           /'High temperature range (Tcom to Thigh) = ',g12.5,'K to',g12.5,'K', &
           /'Polynomial coefficients:',&
           /'        Low range     High range',&
           /' A(1)  ', g12.5, 2x,g12.5, &
           /' A(2)  ', g12.5, 2x,g12.5, &
           /' A(3)  ', g12.5, 2x,g12.5, &
           /' A(4)  ', g12.5, 2x,g12.5, &
           /' A(5)  ', g12.5, 2x,g12.5, &
           /' A(6)  ', g12.5, 2x,g12.5, &
           /' A(7)  ', g12.5, 2x,g12.5, &
           /' Molecular weight = ',g12.5, &
           /' Heat of formation/R = ',g12.5 &
     )

1010 FORMAT(/8x,'Summary of thermochemical data read from the database:', &
           /'Species: ',A      , &
           /'Phase index :',I2,', species index:',I2   , &
           /'Minimum allowed temperature (TMIN) = ',g12.5,'K', &
           /'Maximum allowed temperature (TMAX) = ',g12.5,'K', &
           /'Low  temperature  = ',g12.5,'K' &
     )

      END SUBROUTINE Echo_database_entry

END MODULE READ_DATABASE_MOD
