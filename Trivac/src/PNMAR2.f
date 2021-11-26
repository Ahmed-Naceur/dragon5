*DECK PNMAR2
      FUNCTION PNMAR2(NGPT,L1,L2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Return the dual Marshak boundary coefficients in plane geometry.
* These coefficients are specific to the left boundary.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NGPT    number of Gauss-Legendre base points for the integration of
*         the direction cosine. Set to 65 for exact integration.
* L1      first Legendre order (even number in mixed dual cases).
* L2      second Legendre order (odd number in mixed dual cases).

*Parameters: output
* PNMAR2  Marshak coefficient.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGPT,L1,L2
      REAL PNMAR2
*----
*  LOCAL VARIABLES
*----
      PARAMETER(MAXGPT=64)
      REAL ZGKSI(MAXGPT),WGKSI(MAXGPT)
      DOUBLE PRECISION SUM,PNL1,PNL2,P1,P2
*
      IF(MOD(L1,2).EQ.0) THEN
         CALL XABORT('PNMAR2: ODD FIRST INDEX EXPECTED.')
      ENDIF
      PNL1=0.0D0
      PNL2=0.0D0
      IF(NGPT.LE.64) THEN
*        USE A GAUSS-LEGENDRE QUADRATURE.
         CALL ALGPT(NGPT,-1.0,1.0,ZGKSI,WGKSI)
         SUM=0.0
         DO 30 I=NGPT/2+1,NGPT
            P1=1.0D0
            P2=ZGKSI(I)
            IF(L1.EQ.0) THEN
              PNL1=1.0D0
            ELSE IF(L1.EQ.1) THEN
              PNL1=P2
            ELSE
              DO 10 LL=2,L1
              PNL1=(ZGKSI(I)*REAL(2*LL-1)*P2-REAL(LL-1)*P1)/REAL(LL)
              P1=P2
              P2=PNL1
   10         CONTINUE
            ENDIF
            P1=1.0D0
            P2=ZGKSI(I)
            IF(L2.EQ.0) THEN
              PNL2=1.0D0
            ELSE IF(L2.EQ.1) THEN
              PNL2=P2
            ELSE
              DO 20 LL=2,L2
              PNL2=(ZGKSI(I)*REAL(2*LL-1)*P2-REAL(LL-1)*P1)/REAL(LL)
              P1=P2
              P2=PNL2
   20         CONTINUE
            ENDIF
            SUM=SUM+WGKSI(I)*ZGKSI(I)*(PNL1*PNL2)
   30    CONTINUE
         PNMAR2=REAL(SUM*REAL(2*L1+1))
      ELSE
*        USE EXACT INTEGRATION.
         NGPTE=16
         CALL ALGPT(NGPTE,0.0,1.0,ZGKSI,WGKSI)
         SUM=0.0D0
         DO 60 I=1,NGPTE
            P1=1.0D0
            P2=ZGKSI(I)
            IF(L1.EQ.0) THEN
              PNL1=1.0D0
            ELSE IF(L1.EQ.1) THEN
              PNL1=P2
            ELSE
              DO 40 LL=2,L1
              PNL1=(ZGKSI(I)*REAL(2*LL-1)*P2-REAL(LL-1)*P1)/REAL(LL)
              P1=P2
              P2=PNL1
   40         CONTINUE
            ENDIF
            P1=1.0D0
            P2=ZGKSI(I)
            IF(L2.EQ.0) THEN
              PNL2=1.0D0
            ELSE IF(L2.EQ.1) THEN
              PNL2=P2
            ELSE
              DO 50 LL=2,L2
              PNL2=(ZGKSI(I)*REAL(2*LL-1)*P2-REAL(LL-1)*P1)/REAL(LL)
              P1=P2
              P2=PNL2
   50         CONTINUE
            ENDIF
            SUM=SUM+WGKSI(I)*ZGKSI(I)*(PNL1*PNL2)
   60    CONTINUE
         PNMAR2=REAL(SUM*REAL(2*L1+1))
      ENDIF
      RETURN
      END
