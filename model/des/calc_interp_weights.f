MODULE CALC_INTERP_WEIGHTS_MOD

   use compar, only: FUNIJK_MAP_C
   use discretelement, only: DES_POS_NEW, MAX_PIP, PIJK, des_usr_var
   use discretelement, only: XE, YN, ZT
   use discretelement, only: SUPER_R, SUPER_MN, SUPER_Q, SuperDEM
   use geometry, only: DX, DY, DZ
   use geometry, only: oDX_E, oDY_N, oDZ_T, DO_K
   use param1, only: ZERO, HALF, ONE
   use particle_filter, only: DES_INTERP_DPVM
   use particle_filter, only: DES_INTERP_GAUSS
   use particle_filter, only: DES_INTERP_LHAT
   use particle_filter, only: DES_INTERP_SATELLITE

   use particle_filter, only: DES_INTERP_SCHEME_ENUM
   use particle_filter, only: FILTER_CELL
   use particle_filter, only: FILTER_WEIGHT
   use particle_filter, only: FILTER_WIDTH_INTERP
   use particle_filter, only: FILTER_WIDTH_INTERPx3
   use particle_filter, only: OoFILTER_VOL
   use SQ_SATELLITES
   use SQ_PROPERTIES_MOD
CONTAINS


   INCLUDE 'functions.inc'


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_INTERP_WEIGHTS                                     !
!                                                                      !
!  Purpose: Calculate weights used for interpolation.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE CALC_INTERP_WEIGHTS

      SELECT CASE(DES_INTERP_SCHEME_ENUM)
      CASE(DES_INTERP_DPVM);  CALL CALC_INTERP_WEIGHTS1
      CASE(DES_INTERP_GAUSS); CALL CALC_INTERP_WEIGHTS1
      CASE(DES_INTERP_LHAT);  CALL CALC_INTERP_WEIGHTS2
      CASE(DES_INTERP_SATELLITE); CALL CALC_INTERP_WEIGHTS3
      END SELECT

      RETURN
   END SUBROUTINE CALC_INTERP_WEIGHTS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_INTERP_WEIGHTS1                                    !
!                                                                      !
!  Purpose: Calculate weights used for interpolation.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE CALC_INTERP_WEIGHTS1

      IMPLICIT NONE

      INTEGER :: L, IDX
      INTEGER :: IJK, IJKt

      INTEGER :: I, J, K
      INTEGER :: IC, JC, KC
      INTEGER :: Km, Kp

      DOUBLE PRECISION :: WEIGHT
      DOUBLE PRECISION :: WEIGHT_I(-1:1)
      DOUBLE PRECISION :: WEIGHT_J(-1:1)
      DOUBLE PRECISION :: WEIGHT_K(-1:1)

!$omp parallel default(none) private(L, WEIGHT_I, WEIGHT_J, WEIGHT_K,  &
!$omp    KM, KP, IDX, IJK, I, J, K, WEIGHT, IJKT)                      &
!$omp shared(MAX_PIP, PIJK, DES_POS_NEW, XE, YN, ZT, DO_K,        &
!$omp    FILTER_CELL, FILTER_WEIGHT)
!$omp do
      DO L = 1, MAX_PIP

         IF(.NOT.IS_NORMAL(L)) CYCLE

         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

! Tentative weights for I indices to the West and East.
         WEIGHT_I(-1) = CALC_FILTER_WEIGHTS(DES_POS_NEW(L,1), XE(I-1))
         WEIGHT_I( 1) = CALC_FILTER_WEIGHTS(XE(I), DES_POS_NEW(L,1))
         WEIGHT_I( 0) = ONE - WEIGHT_I(-1) - WEIGHT_I(1)

! Tentative weights for J indices to the South and North.
         WEIGHT_J(-1) = CALC_FILTER_WEIGHTS(DES_POS_NEW(L,2), YN(J-1))
         WEIGHT_J( 1) = CALC_FILTER_WEIGHTS(YN(J), DES_POS_NEW(L,2))
         WEIGHT_J( 0) = ONE - WEIGHT_J(-1) - WEIGHT_J(1)

! Tentative weights for K indices to the Top and Bottom.
         IF(DO_K) THEN
            Km=-1;  Kp=1
            WEIGHT_K(-1) = CALC_FILTER_WEIGHTS(DES_POS_NEW(L,3),ZT(K-1))
            WEIGHT_K( 1) = CALC_FILTER_WEIGHTS(ZT(K), DES_POS_NEW(L,3))
            WEIGHT_K( 0) = ONE - WEIGHT_K(-1) - WEIGHT_K(1)
         ELSE
            Km= 0; Kp=0
            WEIGHT_K( 0) = ONE
         ENDIF

! Set the default fluid cell index and loop counter.
         IJK = PIJK(L,4)
         IDX=0

! Calculate weights for ghost particles. Only store weights that the
! current process owns.
         DO KC=Km,Kp
         DO IC=-1,+1
         DO JC=-1,+1
            IDX=IDX+1
            WEIGHT = WEIGHT_I(IC)*WEIGHT_J(JC)*WEIGHT_K(KC)
            IF(IS_ON_MYPE_PLUS2LAYERS(I+IC,J+JC,K+KC)) THEN
               IJKt = FUNIJK(I+IC,J+JC,K+KC)
               IF(FLUID_AT(IJKt) .OR. CYCLIC_AT(IJKt)) THEN
                  FILTER_CELL(IDX,L) = IJKt
                  FILTER_WEIGHT(IDX,L) = WEIGHT
               ELSE
                  FILTER_CELL(IDX,L) = IJK
                  FILTER_WEIGHT(IDX,L) = WEIGHT
               ENDIF
            ELSE
               FILTER_CELL(IDX,L) = IJK
               FILTER_WEIGHT(IDX,L) = ZERO
            ENDIF
         ENDDO
         ENDDO
         ENDDO

      ENDDO
!$omp end do
!$omp end parallel

      CONTAINS

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      DOUBLE PRECISION FUNCTION CALC_FILTER_WEIGHTS(POS1, POS2)

      DOUBLE PRECISION, INTENT(IN) :: POS1, POS2
      DOUBLE PRECISION :: OVERLAP, HEIGHT

      OVERLAP = POS1 - POS2
      IF(OVERLAP < FILTER_WIDTH_INTERP) THEN
         HEIGHT = FILTER_WIDTH_INTERP - OVERLAP
         CALC_FILTER_WEIGHTS = HEIGHT**2 * &
            (FILTER_WIDTH_INTERPx3 - HEIGHT)*OoFILTER_VOL
      ELSE
         CALC_FILTER_WEIGHTS = ZERO
      ENDIF

      END FUNCTION CALC_FILTER_WEIGHTS

   END SUBROUTINE CALC_INTERP_WEIGHTS1



!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!                                                                             !
!  Subroutine name: CALC_INTERP_WEIGHTS3                                      !
!  Author: Xi Gao                                  Date: 25-Feb-2019          !
!                                                                             !
!  Purpose: Calculate the weights for superquadric particle                    !
!                                                                             !
!  Reference:                                                                 !
!                                                                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

   SUBROUTINE CALC_INTERP_WEIGHTS3


      IMPLICIT NONE

      INTEGER :: L, IDX
      INTEGER :: IJK, IJKt

      INTEGER :: I, J, K
      INTEGER :: IC, JC, KC
      INTEGER :: Km, Kp

      DOUBLE PRECISION :: WEIGHT
      DOUBLE PRECISION :: WEIGHT_I(-1:1)
      DOUBLE PRECISION :: WEIGHT_J(-1:1)
      DOUBLE PRECISION :: WEIGHT_K(-1:1)
      INTEGER :: t, n_left, n_middle, n_right
      INTEGER, parameter:: N_x=3,N_y=3,N_z=3,N_t=27
      DOUBLE PRECISION ::  s_pos(N_t,3),s_weights(N_t)
      DOUBLE PRECISION ::  axi(3), m,n, qc(4),s_pos_global(n_t,3)
      INTEGER :: shape_L
      LOGICAL, PARAMETER :: ldebug = .false.

!$omp parallel default(none) private(L, WEIGHT_I, WEIGHT_J, WEIGHT_K,  &
!$omp    KM, KP, IDX, IJK, I, J, K, WEIGHT, IJKT, AXI,M,N,Qc,          &
!$omp    n_left,n_middle,n_right)                                      &
!$omp shared(MAX_PIP, PIJK, DES_POS_NEW, XE, YN, ZT, DO_K,             &
!$omp    FILTER_CELL, FILTER_WEIGHT,SUPER_R,SUPER_MN,SUPER_Q,S_POS,    &
!$omp    S_WEIGHTS,S_POS_GLOBAL)
!$omp do
      DO L = 1, MAX_PIP

         IF(.NOT.IS_NORMAL(L)) CYCLE

         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)
         IJK = PIJK(L,4)

! Superquadric particle properties
         axi(1)=super_r(L,1)
         axi(2)=super_r(L,2)
         axi(3)=super_r(L,3)
         m=super_mn(L,1)
         n=super_mn(L,2)
         Qc(1)= super_q(L,1)
         Qc(2)= super_q(L,2)
         Qc(3)= super_q(L,3)
         Qc(4)= super_q(L,4)
! Calculate the satellites points in the local coordinates
      call Sq_satellites_local(n_x,n_y,n_z,n_t,axi(1),&
                                          axi(2),axi(3),&
                                          m, n,s_pos,s_weights)
! Calculate the satellites in the global coorrdinate
      call Sq_satellites_global(N_t, Qc,s_pos, s_pos_global)
      DO t=1, n_t
          s_pos_global(t,:)=s_pos_global(t,:)+DES_POS_NEW(L,:)
      ENDDO

! Tentative weights for I indices to the West and East.
      N_left=0
      N_middle=0
      N_right=0
      DO t=1,n_t
          if (s_weights(t)>0.0) then
              if (s_pos_global(t,1)<=XE(I-1)) then
                  n_left=n_left+1
              elseif (s_pos_global(t,1)>XE(I-1) .and. s_pos_global(t,1)<XE(I)) then
                  n_middle=n_middle + 1
              else
                  n_right = n_right + 1
              endif
          endif
      ENDDO
! 1.0 is used to convert integer to double precision
      WEIGHT_I(-1)=1.0d0*n_left/(n_left+n_middle+n_right)
      WEIGHT_I( 0)=1.0d0*n_middle/(n_left+n_middle+n_right)
      WEIGHT_I( 1)=ONE-WEIGHT_I(-1)-WEIGHT_I(0)
      if (ldebug) then
          WRITE(*,*) 'weight_i'
          WRITE(*,'(3(f8.3))')  n_left*1.0, n_middle*1.0,n_right*1.0
          WRITE(*,'(3(f8.3))')  weight_i(-1), weight_i(0),weight_i(1)
      endif
! Tentative weights for J indices to the South and North.
      N_left=0
      N_middle=0
      N_right=0
      DO t=1,n_t
          if (s_weights(t)>0.0) then
              if (s_pos_global(t,2)<=YN(J-1)) then
                  n_left=n_left+1
              else if (s_pos_global(t,2)>YN(J-1) .and. s_pos_global(t,2)<YN(J)) then
                  n_middle=n_middle + 1
              else
                  n_right = n_right + 1
              endif
          endif
      ENDDO
      WEIGHT_J(-1)=n_left*1.0d0/(n_left+n_middle+n_right)
      WEIGHT_J( 0)=n_middle*1.0d0/(n_left+n_middle+n_right)
      WEIGHT_J( 1)=ONE-WEIGHT_J(-1)-WEIGHT_J(0)
      if (ldebug) then
          WRITE(*,*) 'weight_j'
          WRITE(*,'(3(f8.3))')  n_left*1.0, n_middle*1.0,n_right*1.0
          WRITE(*,'(3(f8.3))')  weight_j(-1), weight_j(0),weight_j(1)
      endif
! Tentative weights for K indices to the Top and Bottom.
      IF(DO_K) THEN
          Km=-1;  Kp=1

          N_left=0
          N_middle=0
          N_right=0
          DO t=1,n_t
              if (s_weights(t)>0.0) then
                  if (s_pos_global(t,3)<=ZT(K-1)) then
                      n_left=n_left+1
                  else if (s_pos_global(t,3)>ZT(K-1) .and. s_pos_global(t,3)<ZT(K)) then
                      n_middle=n_middle + 1
                  else
                      n_right = n_right + 1
                  endif
              endif
          ENDDO
          WEIGHT_K(-1)=1.0d0*n_left/(n_left+n_middle+n_right)
          WEIGHT_K( 0)=1.0d0*n_middle/(n_left+n_middle+n_right)
          WEIGHT_K( 1)=ONE-WEIGHT_K(-1)-WEIGHT_K(0)
          if (ldebug) then
              WRITE(*,*) 'weight_k'
              WRITE(*,'(3(f8.3))')  n_left*1.0, n_middle*1.0,n_right*1.0
              WRITE(*,'(3(f8.3))')  weight_k(-1), weight_k(0),weight_k(1)
          endif
      ELSE
          Km= 0; Kp=0
          WEIGHT_K( 0) = ONE
      ENDIF

! Set the default fluid cell index and loop counter.
      IJK = PIJK(L,4)
      IDX=0

! Calculate weights for ghost particles. Only store weights that the
! current process owns.
      DO KC=Km,Kp
          DO IC=-1,+1
              DO JC=-1,+1
                  IDX=IDX+1
                  WEIGHT = WEIGHT_I(IC)*WEIGHT_J(JC)*WEIGHT_K(KC)
                  if (ldebug) then
                      WRITE(*,*) 'weight_ijk'
                      WRITE(*,'(1(f8.3))')  weight
                  endif
                  IF(IS_ON_MYPE_PLUS2LAYERS(I+IC,J+JC,K+KC)) THEN
                      IJKt = FUNIJK(I+IC,J+JC,K+KC)
                      IF(FLUID_AT(IJKt) .OR. CYCLIC_AT(IJKt)) THEN
                          FILTER_CELL(IDX,L) = IJKt
                          FILTER_WEIGHT(IDX,L) = WEIGHT
                      ELSE
                          FILTER_CELL(IDX,L) = IJK
                          FILTER_WEIGHT(IDX,L) = WEIGHT
                      ENDIF
                  ELSE
                      FILTER_CELL(IDX,L) = IJK
                      FILTER_WEIGHT(IDX,L) = ZERO
                  ENDIF
              ENDDO
          ENDDO
      ENDDO
  ENDDO
!$omp end do
!$omp end parallel

END SUBROUTINE CALC_INTERP_WEIGHTS3

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_INTERP_WEIGHTS2                                    !
!                                                                      !
!  Purpose: Calculate weights used for interpolation.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
SUBROUTINE CALC_INTERP_WEIGHTS2

      IMPLICIT NONE

      INTEGER :: IDX
      INTEGER :: IJKt

! Interpolation weight
      DOUBLE PRECISION :: WEIGHT
! temporary variables used to calculate pressure at scalar cell edge
      INTEGER :: IJK, NP

      INTEGER :: I, J, K
      INTEGER :: lI, lJ, lK
      INTEGER :: IC, JC, KC
      INTEGER :: Kub
      DOUBLE PRECISION :: lOoDX, lOoDY, lOoDZ

      DOUBLE PRECISION :: XC, YC, ZC
!      DOUBLE PRECISION :: Sx(0:1), Sy(0:1), Sz(0:1)
!     new - adding 2 to keep seamless merge command for (27, 8) = (3D, 2D)
      DOUBLE PRECISION :: Sx(0:2), Sy(0:2), Sz(0:2)
      DOUBLE PRECISION, PARAMETER :: SGN(0:1) = (/-1.0d0, 1.0d0/)
!......................................................................!


      DO NP = 1, MAX_PIP

          IF(.NOT.IS_NORMAL(NP)) CYCLE

          I = PIJK(NP,1)
          J = PIJK(NP,2)
          K = PIJK(NP,3)

          IJK = PIJK(NP,4)

! Tentative weights for I indices to the West and East.
          XC = XE(I-1) + HALF*DX(I)
          IF(DES_POS_NEW(NP,1) < XC) THEN
              IF(FLUID_AT(WEST_OF(IJK))) THEN
                  lI = I-1
                  lOoDX = oDX_E(lI)
                  Sx(0) = (XC - DES_POS_NEW(NP,1))*lOoDX
                  Sx(1) = ONE - Sx(0)
                  Sx(2) = 0.0d0 !new
              ELSE
                  lI = I
                  lOoDX = oDX_E(lI)
                  !Sx(0:1)=(/1.0d0,0.0d0/)
                  Sx(0:2)=(/0.d0,1.0d0,0.0d0/)     !updated
            ENDIF
         ELSE
             IF(FLUID_AT(EAST_OF(IJK))) THEN
                 lI = I
                 lOoDX = oDX_E(lI)
                 Sx(0) = 0.0d0                    !updated
                 Sx(2) = (DES_POS_NEW(NP,1) - XC)*lOoDX
                 Sx(1) = ONE - Sx(2)              !updated
             ELSE
                 lI = I-1
                 lOoDX = oDX_E(lI)
                 !Sx(0:1)=(/0.0d0,1.0d0/)
                 Sx(0:2)=(/0.0d0,1.0d0,0.0d0/)
             ENDIF
         ENDIF

! Tentative weights for J indices to the South and North.
         YC = YN(J-1) + HALF*DY(J)
         IF(DES_POS_NEW(NP,2) < YC) THEN
             IF(FLUID_AT(SOUTH_OF(IJK))) THEN
                 lJ = J-1
                 lOoDY = oDY_N(lJ)
                 Sy(0) = (YC - DES_POS_NEW(NP,2))*lOoDY
                 Sy(1) = ONE - Sy(0)
                 Sy(2) = 0.0d0  !new
             ELSE
                 lJ = J
                 lOoDY = oDY_N(lJ)
                 !Sy(0:1)=(/1.0d0,0.0d0/)
                 Sy(0:2)=(/0.0d0,1.0d0,0.0d0/)     !updated
             ENDIF
         ELSE
             IF(FLUID_AT(NORTH_OF(IJK))) THEN
                 lJ = J
                 lOoDY = oDY_N(lJ)
                 Sy(0) = 0.0d0                     !updated
                 Sy(2) = (DES_POS_NEW(NP,2) - YC)*lOoDY
                 Sy(1) = ONE - Sy(2)               !updated
             ELSE
                 lJ = J-1
                 lOoDY = oDY_N(lJ)
                 !Sy(0:1)=(/0.0d0,1.0d0/)
               Sy(0:2)=(/0.0d0,1.0d0,0.0d0/)
           ENDIF
       ENDIF

! Tentative weights for K indices to the Top and Bottom.
       IF(DO_K) THEN
           ZC = ZT(K-1) + HALF*DZ(K)
           IF(DES_POS_NEW(NP,3) < ZC) THEN
               IF(FLUID_AT(BOTTOM_OF(IJK))) THEN
                   lK = K-1
                   lOoDZ = oDZ_T(lK)
                   Sz(0) = (ZC - DES_POS_NEW(NP,3))*lOoDZ
                   Sz(1) = ONE - Sz(0)
                   Sz(2) = 0.0d0  !new
               ELSE
                   lK = K
                   lOoDZ = oDZ_T(lK)
                   !Sz(0:1)=(/1.0d0,0.0d0/)
                   Sz(0:2)=(/0.0d0,1.0d0,0.0d0/)  !updated
               ENDIF
           ELSE
               IF(FLUID_AT(TOP_OF(IJK))) THEN
                   lK = K
                   lOoDZ = oDZ_T(lK)
                   Sz(0) = 0.0d0                  !updated
                   Sz(2) = (DES_POS_NEW(NP,3) - ZC)*lOoDZ
                   Sz(1) = ONE - Sz(2)            !updated
               ELSE
                   lK = K-1
                   lOoDZ = oDZ_T(lK)
                   !Sz(0:1)=(/0.0d0,1.0d0/)
                  Sz(0:2)=(/0.0d0,1.0d0,0.0d0/)
               ENDIF
            ENDIF
            !Kub = 1
            Kub = 2
         ELSE
            lK = K
            lOoDZ = 0.0d0
            Sz(:) = ONE
            Kub = 0
         ENDIF

! Set the default fluid cell index and loop counter.
         IDX=0

! Calculate weights for ghost particles. Only store weights that the
! current process owns.

!         write(*,*) "kub=", Kub

         DO KC=0,Kub
             DO JC=0,2   !changed from 1 to accommodate merge
                 DO IC=0,2   !changed from 1 to accommodate merge
                     IDX=IDX+1
                     WEIGHT = Sx(IC)*Sy(JC)*Sz(KC)
                     IF(IS_ON_MYPE_PLUS2LAYERS(I-1+IC,J-1+JC,K-1+KC)) THEN
                         IJKt = FUNIJK_MAP_C(I-1+IC,J-1+JC,K-1+KC)
                         IF(FLUID_AT(IJKt)) THEN
                             FILTER_CELL(IDX,NP) = IJKt
                             FILTER_WEIGHT(IDX,NP) = WEIGHT
                         ELSE
                             FILTER_CELL(IDX,NP) = IJK
                             FILTER_WEIGHT(IDX,NP) = ZERO
                         ENDIF
                     ELSE
                         FILTER_CELL(IDX,NP) = IJK
                         FILTER_WEIGHT(IDX,NP) = ZERO
                     ENDIF
!            if (IDX.LE.10) then
!              if (IS_ON_MYPE_PLUS2LAYERS(lI+IC,lJ+JC,lK+KC)) then
!                 write(*,*) "TRUE"
!              else
!                 write(*,*) "FALSE"
!              endif
!              if (FLUID_AT(IJKt)) write(*,*) "There is fluid"
!              WRITE(*,*) "(IC,JC,KC)=",IC,",",JC,",",KC
!              WRITE(*,*) "IJKt = ", IJKt
!              WRITE(*,*) "IJK = ", IJK
!            end if
                 ENDDO
             ENDDO
         ENDDO
     ENDDO

END SUBROUTINE CALC_INTERP_WEIGHTS2

END MODULE CALC_INTERP_WEIGHTS_MOD
