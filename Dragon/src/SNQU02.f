*DECK SNQU02
      SUBROUTINE SNQU02(NLF,JOP,U,W,TPQ,UPQ,VPQ,WPQ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Set the level-symmetric (type 2) quadratures.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NLF     order of the SN approximation (even number).
*
*Parameters: output
* JOP     number of base points per axial level in one octant.
* U       base points in $\\xi$ of the axial quadrature. Used with
*         zero-weight points.
* W       weights for the axial quadrature in $\\xi$.
* TPQ     base points in $\\xi$ of the 2D SN quadrature.
* UPQ     base points in $\\mu$ of the 2D SN quadrature.
* VPQ     base points in $\\eta$ of the 2D SN quadrature.
* WPQ     weights of the 2D SN quadrature.
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NLF,JOP(NLF/2)
      REAL U(NLF/2),W(NLF/2),TPQ(NLF*(NLF/2+1)/4),UPQ(NLF*(NLF/2+1)/4),
     1 VPQ(NLF*(NLF/2+1)/4),WPQ(NLF*(NLF/2+1)/4)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(PI=3.141592654,MAXNLF=24,MAXEQ=64,MAXNBA=78,MAXW=16)
      INTEGER INWEI(MAXNBA)
      DOUBLE PRECISION WSUM2,WEI(MAXW),ZMAT(MAXEQ,MAXW+1),UD(MAXW)
*----
*  SET THE UNIQUE QUADRATURE VALUES.
*----
      IF(NLF.GT.MAXNLF) CALL XABORT('SNQU02: MAXNLF OVERFLOW.')
      M2=NLF/2
      NPQ=M2*(M2+1)/2
      ZMU1=1.0D0/(3.0D0*DBLE(NLF-1))
      NW=1+(NLF*(NLF+8)-1)/48
      IF(NW.GT.MAXW) CALL XABORT('SNQU02: MAXW OVERFLOW.')
      IF(NLF.EQ.2) THEN
         ZMU1=0.33333333
      ELSE IF(NLF.EQ.4) THEN
         ZMU1=0.12251480
      ELSE IF(NLF.EQ.6) THEN
         ZMU1=0.07109447
      ELSE IF(NLF.EQ.8) THEN
         ZMU1=0.04761903
      ELSE IF(NLF.EQ.10) THEN
         ZMU1=0.03584310
      ELSE IF(NLF.EQ.12) THEN
         ZMU1=0.02796615
      ELSE IF(NLF.EQ.14) THEN
         ZMU1=0.02310250
      ELSE IF(NLF.EQ.16) THEN
         ZMU1=0.01931398
      ELSE IF(NLF.EQ.18) THEN
         ZMU1=0.01692067
      ELSE IF(NLF.EQ.20) THEN
         ZMU1=0.01455253
      ELSE
         CALL XABORT('SNQU02: ORDER NOT AVAILABLE.')
      ENDIF
      U(1)=REAL(SQRT(ZMU1))
      DO I=2,M2
         ZMU2=ZMU1+2.0D0*DBLE(I-1)*(1.0D0-3.0D0*ZMU1)/DBLE(NLF-2)
         U(I)=REAL(SQRT(ZMU2))
      ENDDO
*----
*  COMPUTE THE POSITION OF WEIGHTS.
*----
      IPR=0
      INMAX=0
      DO IP=1,M2
         JOP(IP)=M2-IP+1
         DO IQ=1,JOP(IP)
            IPR=IPR+1
            IF(IPR.GT.MAXNBA) CALL XABORT('SNQU02: MAXNBA OVERFLOW.')
            TPQ(IPR)=U(IP)
            UPQ(IPR)=U(M2+2-IP-IQ)
            VPQ(IPR)=U(IQ)
            IS=MIN(IP,IQ,M2+2-IP-IQ)
            NW0=0
            DO II=1,IS-1
               NW0=NW0+(M2-3*(II-1)+1)/2
            ENDDO
            KK=IP-IS+1
            LL=IQ-IS+1
            IF(KK.EQ.1)THEN
               INWEI(IPR)=NW0+MIN(LL,M2-3*(IS-1)+1-LL)
            ELSEIF(LL.EQ.1)THEN
               INWEI(IPR)=NW0+MIN(KK,M2-3*(IS-1)+1-KK)
            ELSE
               INWEI(IPR)=NW0+MIN(KK,LL)
            ENDIF
            INMAX=MAX(INMAX,INWEI(IPR))
         ENDDO
      ENDDO
      IF(INMAX.NE.NW) CALL XABORT('SNQU02: INVALID VALUE OF NW.')
      IF(IPR.NE.NPQ) CALL XABORT('SNQU02: BAD VALUE ON NPQ.')
*----
*  SET THE RECTANGULAR SYSTEM AND SOLVE IT USING THE QR METHOD.
*----
      NEQ=0
      DO IPL=0,NLF,2
        DO IPK=IPL,NLF-IPL,2
          IF(MOD(IPL+IPK,2).EQ.1) CYCLE
          NEQ=NEQ+1
          IF(NEQ.GT.MAXEQ) CALL XABORT('SNQU02: MAXEQ OVERFLOW.')
          DO IW=1,NW
             ZMAT(NEQ,IW)=0.0D0
          ENDDO
          DO IPQ=1,NPQ
             ZMU=TPQ(IPQ)
             ZETA=UPQ(IPQ)
             IW=INWEI(IPQ)
             ZMAT(NEQ,IW)=ZMAT(NEQ,IW)+(ZMU**IPK)*(ZETA**IPL)
          ENDDO
          REF=1.0D0/DBLE(IPK+IPL+1)
          DO I=1,IPL-1,2
             REF=REF*DBLE(I)/DBLE(IPK+I)
          ENDDO
          ZMAT(NEQ,NW+1)=REF
        ENDDO
      ENDDO
      CALL ALST2F(MAXEQ,NEQ,NW,ZMAT,UD)
      CALL ALST2S(MAXEQ,NEQ,NW,ZMAT,UD,ZMAT(1,NW+1),WEI)
*----
*  SET THE LEVEL-SYMMETRIC QUADRATURES.
*----
      IPQ=0
      WSUM=0.0
      DO IP=1,M2
         WSUM2=0.0D0
         DO IQ=1,JOP(IP)
            IPQ=IPQ+1
            WPQ(IPQ)=REAL(WEI(INWEI(IPQ))*PI/2.0)
            WSUM2=WSUM2+WEI(INWEI(IPQ))
         ENDDO
         W(IP)=REAL(WSUM2)
         WSUM=WSUM+REAL(WSUM2*PI/2.0)
      ENDDO
      RETURN
      END
