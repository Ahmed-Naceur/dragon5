*DECK USSSEK
      SUBROUTINE USSSEK(NBNRS,NQT,LMOD,SIGR,CONRL,WEIGH,SIGL,PIJK,DIL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the dilution matrix preserving the non-correlated collision
* probability matrix in each subgroup. Use a fixed point iteration.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NBNRS   number of correlated fuel regions.
* NQT     number of subgroups in admixed resonant isotope.
* LMOD    moderator flag (=.TRUE. if all regions are containing the
*         resonant isotopes; =.FALSE. if a moderator region exists).
* SIGR    macroscopic total xs of the other isotopes.
* CONRL   number density of the admixed resonant isotope.
* WEIGH   multiband weights for the admixed resonant isotope.
* SIGL    microscopic total xs of the admixed resonant isotope.
* PIJK    non-correlated collision probability matrix.
*
*Parameters: input/output
* DIL     estimate and converged value of the dilution matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      LOGICAL LMOD
      INTEGER NBNRS,NQT
      REAL    SIGR(NBNRS),CONRL(NBNRS),WEIGH(NQT),SIGL(NQT),
     1        PIJK(0:NBNRS,0:NBNRS),DIL(0:NBNRS,0:NBNRS)
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: WORK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK(0:NBNRS,0:NBNRS,3))
*
      DEN=0.0
      DO 20 I=0,NBNRS
      DO 10 J=0,NBNRS
      WORK(I,J,3)=PIJK(I,J)
      DEN=MAX(DEN,ABS(PIJK(I,J)))
   10 CONTINUE
   20 CONTINUE
      IF(LMOD) THEN
         CALL ALINV(NBNRS,WORK(1,1,3),NBNRS+1,IER)
      ELSE
         CALL ALINV(NBNRS+1,WORK(0,0,3),NBNRS+1,IER)
      ENDIF
      IF(IER.NE.0) CALL XABORT('USSSEK: SINGULAR MATRIX(1).')
      ITER=0
   30 ITER=ITER+1
      IF(ITER.GT.50) CALL XABORT('USSSEK: MAXIMUM NB. OF ITERATIONS.')
      DO 45 I=0,NBNRS
      DO 40 J=0,NBNRS
      WORK(I,J,1)=0.0
   40 CONTINUE
   45 CONTINUE
      DO 72 L=1,NQT
      DO 55 I=0,NBNRS
      DO 50 J=0,NBNRS
      WORK(I,J,2)=DIL(I,J)
   50 CONTINUE
   55 CONTINUE
      DO 60 I=1,NBNRS
      WORK(I,I,2)=WORK(I,I,2)+SIGR(I)+CONRL(I)*SIGL(L)
   60 CONTINUE
      IF(LMOD) THEN
         CALL ALINV(NBNRS,WORK(1,1,2),NBNRS+1,IER)
      ELSE
         CALL ALINV(NBNRS+1,WORK(0,0,2),NBNRS+1,IER)
      ENDIF
      IF(IER.NE.0) CALL XABORT('USSSEK: SINGULAR MATRIX(2).')
      DO 71 I=0,NBNRS
      DO 70 J=0,NBNRS
      WORK(I,J,1)=WORK(I,J,1)+WEIGH(L)*WORK(I,J,2)
   70 CONTINUE
   71 CONTINUE
   72 CONTINUE
      ERR=0.0
      DO 85 I=0,NBNRS
      DO 80 J=0,NBNRS
      ERR=MAX(ERR,ABS(PIJK(I,J)-WORK(I,J,1)))
   80 CONTINUE
   85 CONTINUE
      IF(ERR.LT.1.0E-4*DEN) GO TO 110
      IF(LMOD) THEN
         CALL ALINV(NBNRS,WORK(1,1,1),NBNRS+1,IER)
      ELSE
         CALL ALINV(NBNRS+1,WORK(0,0,1),NBNRS+1,IER)
      ENDIF
      IF(IER.NE.0) CALL XABORT('USSSEK: SINGULAR MATRIX(3).')
      DO 105 I=0,NBNRS
      DO 100 J=0,NBNRS
      DIL(I,J)=DIL(I,J)+WORK(I,J,3)-WORK(I,J,1)
  100 CONTINUE
  105 CONTINUE
      GO TO 30
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  110 DEALLOCATE(WORK)
      RETURN
      END
