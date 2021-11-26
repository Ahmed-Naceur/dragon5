*DECK ALVDLF
      SUBROUTINE ALVDLF (ASS,MU1,ISEG,LON,NBL,LBL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* in-place L-D-L(T) factorization of a symmetric positive definite
* matrix in compressed diagonal storage mode. Supervectorial version
*
*Copyright:
* Copyright (C) 1993 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* ASS     coefficient matrix in compressed diagonal storage mode.
*         DIMENSION ASS(ISEG,MU1(L4))
* MU1     position of each diagonal element in vector ASS.
*         DIMENSION MU1(L4) where L4=SUM(LBL(I))
* ISEG    number of elements in a vector register.
* LON     number of groups of linear systems.
* NBL     number of linear systems in each group. DIMENSION NBL(LON)
* LBL     number of unknowns in each group. DIMENSION LBL(LON)
*
*Parameters: output
* ASS     LDL(T) factors of the coefficient matrix in compressed
*         diagonal storage mode.
*
*-----------------------------------------------------------------------
*
      INTEGER ISEG,LON,MU1(*),NBL(LON),LBL(LON)
      REAL ASS(ISEG,*)
      REAL, DIMENSION(:), ALLOCATABLE :: T
*
      ALLOCATE(T(ISEG))
      LBL0=0
      DO ILON=1,LON
        L4=LBL(ILON)
        NBANC=NBL(ILON)
        DO IB=1,NBANC
           ASS(IB,MU1(LBL0+1))=1.0/ASS(IB,MU1(LBL0+1))
        ENDDO
        IF (L4.NE.1) THEN
          DO K=LBL0+2,LBL0+L4
            K1=MU1(K)-K
            KM=MU1(K-1)+1-K1
            IF(KM+1-K .LE. 0) THEN
              IF(KM+1-K .LT. 0) THEN
                DO  I=KM+1,K-1
                  DO IB=1,NBANC
                    T(IB)=ASS(IB,K1+I)
                    ASS(IB,K1+I)=0.0
                  ENDDO
                  I1=MU1(I)-I
                  IM=MU1(I-1)+1-I1
                  IMIN=MAX0(IM,KM)
                  DO J=IMIN,I
                    DO IB=1,NBANC
                      T(IB)=T(IB)-ASS(IB,K1+J)*ASS(IB,I1+J)
                    ENDDO
                  ENDDO
                  DO IB=1,NBANC
                    ASS(IB,K1+I)=T(IB)
                  ENDDO
                ENDDO
              ENDIF
              DO IB=1,NBANC
                T(IB)=0.0
              ENDDO
              DO I=KM,K-1
                DO IB=1,NBANC
                  GAR=ASS(IB,K1+I)
                  ASS(IB,K1+I)=GAR*ASS(IB,MU1(I))
                  T(IB)=T(IB)+GAR*ASS(IB,K1+I)
                ENDDO
              ENDDO
              DO IB=1,NBANC
                ASS(IB,MU1(K))=ASS(IB,MU1(K))-T(IB)
              ENDDO
            ENDIF
            DO IB=1,NBANC
              ASS(IB,MU1(K))=1.0/ASS(IB,MU1(K))
            ENDDO
          ENDDO
        ENDIF
        LBL0=LBL0+L4
      ENDDO
      DEALLOCATE(T)
      RETURN
      END
