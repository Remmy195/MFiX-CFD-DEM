! Definitions for parallel programming

      MODULE parallel
      IMPLICIT NONE

! Chunks of arrays used by each processor
      INTEGER :: Chunk_size
      PARAMETER (chunk_size=8)

      LOGICAL :: USE_DOLOOP
      LOGICAL :: IS_SERIAL
      LOGICAL :: IS_DMP
      LOGICAL :: IS_SMP

      END MODULE parallel
