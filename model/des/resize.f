module resize

  public :: byte_grow
  public :: integer_grow
  public :: integer_grow2_reverse
  public :: integer_grow2
  public :: logical_grow
  public :: logical_grow2
  public :: real_grow
  public :: real_grow2
  public :: real_grow2_reverse
  public :: real_grow3
  public :: logical_grow2_reverse

  contains

      SUBROUTINE BYTE_GROW(byte_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: byte_array
        INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE :: byte_tmp
        INTEGER lSIZE
        INTEGER lBOUND1

        lSIZE = size(byte_array,1)
        lBOUND1 = lbound(byte_array,1)
        allocate(byte_tmp(lBOUND1:new_size+lBOUND1-1))
        byte_tmp(lBOUND1:lSIZE+lBOUND1-1) = byte_array(lBOUND1:lSIZE+lBOUND1-1)
        call move_alloc(byte_tmp,byte_array)

      END SUBROUTINE BYTE_GROW

      SUBROUTINE INTEGER_GROW(integer_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:), ALLOCATABLE :: integer_tmp
        INTEGER lSIZE
        INTEGER lBOUND1

        lSIZE = size(integer_array,1)
        lBOUND1 = lbound(integer_array,1)
        allocate(integer_tmp(lBOUND1:new_size+lBOUND1-1))
        integer_tmp(lBOUND1:lSIZE+lBOUND1-1) = integer_array(lBOUND1:lSIZE+lBOUND1-1)
        call move_alloc(integer_tmp,integer_array)

      END SUBROUTINE INTEGER_GROW

      SUBROUTINE INTEGER_GROW2_reverse(integer_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: integer_tmp
        INTEGER lSIZE, lSIZE2
        INTEGER lBOUND1, lBOUND2

        lSIZE = size(integer_array,1)
        lSIZE2 = size(integer_array,2)
        lBOUND1 = lbound(integer_array,1)
        lBOUND2 = lbound(integer_array,2)
        allocate(integer_tmp(lBOUND1:new_size+lBOUND1-1,lBOUND2:lSIZE2+lBOUND2-1))
        integer_tmp(lBOUND1:lSIZE+lBOUND1-1,:) = integer_array(lBOUND1:lSIZE+lBOUND1-1,:)
        call move_alloc(integer_tmp,integer_array)

      END SUBROUTINE INTEGER_GROW2_reverse

      SUBROUTINE INTEGER_GROW2(integer_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: integer_array
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: integer_tmp
        INTEGER lSIZE, lSIZE2
        INTEGER lBOUND1, lBOUND2

        lSIZE = size(integer_array,1)
        lSIZE2 = size(integer_array,2)
        lBOUND1 = lbound(integer_array,1)
        lBOUND2 = lbound(integer_array,2)
        allocate(integer_tmp(lBOUND1:lSIZE+lBOUND1-1,lBOUND2:new_size+lBOUND2-1))
        integer_tmp(:,lBOUND2:lSIZE2+lBOUND2-1) = integer_array(:,lBOUND2:lSIZE2+lBOUND2-1)
        call move_alloc(integer_tmp,integer_array)

      END SUBROUTINE INTEGER_GROW2

      SUBROUTINE LOGICAL_GROW(logical_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        LOGICAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: logical_array
        LOGICAL, DIMENSION(:), ALLOCATABLE :: logical_tmp
        INTEGER lSIZE
        INTEGER lBOUND1

        lSIZE = size(logical_array,1)
        lBOUND1 = lbound(logical_array,1)
        allocate(logical_tmp(lBOUND1:new_size+lBOUND1-1))
        logical_tmp(lBOUND1:lSIZE+lBOUND1-1) = logical_array(lBOUND1:lSIZE+lBOUND1-1)
        call move_alloc(logical_tmp,logical_array)

      END SUBROUTINE LOGICAL_GROW

      SUBROUTINE LOGICAL_GROW2(logical_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: logical_array
        LOGICAL, DIMENSION(:,:), ALLOCATABLE :: logical_tmp
        INTEGER lSIZE, lSIZE2
        INTEGER lBOUND1, lBOUND2

        lSIZE = size(logical_array,1)
        lSIZE2 = size(logical_array,2)
        lBOUND1 = lbound(logical_array, 1)
        lBOUND2 = lbound(logical_array, 2)
        allocate(logical_tmp(lBOUND1:lSIZE+lBOUND1-1,lBOUND2:new_size+lBOUND2-1))
        logical_tmp(:,lBOUND2:lSIZE2+lBOUND2-1) = logical_array(:,lBOUND2:lSIZE2+lBOUND2-1)
        call move_alloc(logical_tmp,logical_array)

      END SUBROUTINE LOGICAL_GROW2

      SUBROUTINE REAL_GROW(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE
        INTEGER lBOUND1

        lSIZE = size(real_array,1)
        lBOUND1 = lbound(real_array, 1)
        allocate(real_tmp(lBOUND1:new_size+lBOUND1-1))
        real_tmp(lBOUND1:lSIZE+lBOUND1-1) = real_array(lBOUND1:lSIZE+lBOUND1-1)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW

      SUBROUTINE REAL_GROW2(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2
        INTEGER lBOUND1, lBOUND2

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        lBOUND1 = lbound(real_array, 1)
        lBOUND2 = lbound(real_array, 2)
        allocate(real_tmp(lBOUND1:lSIZE+lBOUND1-1,lBOUND2:new_size+lBOUND2-1))
        real_tmp(:,lBOUND2:lSIZE2+lBOUND2-1) = real_array(:,lBOUND2:lSIZE2+lBOUND2-1)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW2

      SUBROUTINE REAL_GROW2_reverse(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2
        INTEGER lBOUND1, lBOUND2

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        lBOUND1 = lbound(real_array, 1)
        lBOUND2 = lbound(real_array, 2)
        allocate(real_tmp(lBOUND1:new_size+lBOUND1-1,lBOUND2:lSIZE2+lBOUND2-1))
        real_tmp(lBOUND1:lSIZE+lBOUND1-1,:) = real_array(lBOUND1:lSIZE+lBOUND1-1,:)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW2_REVERSE

      SUBROUTINE REAL_GROW3(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2, lSIZE3
        INTEGER lBOUND1, lBOUND2, lBOUND3

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        lSIZE3 = size(real_array,3)
        lBOUND1 = lbound(real_array, 1)
        lBOUND2 = lbound(real_array, 2)
        lBOUND3 = lbound(real_array, 3)
        allocate(real_tmp(lBOUND1:lSIZE+lBOUND1-1,lBOUND2:lSIZE2+lBOUND2-1,lBOUND3:new_size+lBOUND3-1))
        real_tmp(:,:,lBOUND3:lSIZE3+lBOUND3-1) = real_array(:,:,lBOUND3:lSIZE3+lBOUND3-1)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW3

      SUBROUTINE REAL_GROW3_REVERSE2(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2, lSIZE3
        INTEGER lBOUND1, lBOUND2, lBOUND3

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        lSIZE3 = size(real_array,3)
        lBOUND1 = lbound(real_array, 1)
        lBOUND2 = lbound(real_array, 2)
        lBOUND3 = lbound(real_array, 3)
        allocate(real_tmp(lBOUND1:lSIZE+lBOUND1-1,lBOUND2:new_size+lBOUND2-1,lBOUND3:lSIZE3+lBOUND3-1))
        real_tmp(:,lBOUND2:lSIZE2+lBOUND2-1,:) = real_array(:,lBOUND2:lSIZE2+lBOUND2-1,:)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW3_REVERSE2

      SUBROUTINE REAL_GROW4(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2, lSIZE3, lSIZE4
        INTEGER lBOUND1, lBOUND2, lBOUND3, lBOUND4

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        lSIZE3 = size(real_array,3)
        lSIZE4 = size(real_array,4)
        lBOUND1 = lbound(real_array, 1)
        lBOUND2 = lbound(real_array, 2)
        lBOUND3 = lbound(real_array, 3)
        lBOUND4 = lbound(real_array, 4)
        allocate(real_tmp(lBOUND1:lSIZE+lBOUND1-1,lBOUND2:lSIZE2+lBOUND2-1,lBOUND3:lSIZE3+lBOUND3-1,lBOUND4:new_size+lBOUND4-1))
        real_tmp(:,:,:,lBOUND4:lSIZE4+lBOUND4-1) = real_array(:,:,:,lBOUND4:lSIZE4+lBOUND4-1)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE REAL_GROW4

      SUBROUTINE LOGICAL_GROW2_REVERSE(real_array,new_size)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: new_size
        LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: real_array
        LOGICAL, DIMENSION(:,:), ALLOCATABLE :: real_tmp
        INTEGER lSIZE, lSIZE2
        INTEGER lBOUND1, lBOUND2

        lSIZE = size(real_array,1)
        lSIZE2 = size(real_array,2)
        lBOUND1 = lbound(real_array, 1)
        lBOUND2 = lbound(real_array, 2)
        allocate(real_tmp(lBOUND1:new_size+lBOUND1-1,lBOUND2:lSIZE2+lBOUND2-1))
        real_tmp(lBOUND1:lSIZE+lBOUND1-1,:) = real_array(lBOUND1:lSIZE+lBOUND1-1,:)
        call move_alloc(real_tmp,real_array)

      END SUBROUTINE LOGICAL_GROW2_REVERSE

end module resize
