*DECK NXTLSN
      SUBROUTINE NXTLSN(NDIM  ,ORDRE ,NQUAD ,NBANGL,DQUAD ,DANGLT,
     >                  DDENWT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To define level-symmetric (type 2) quadrature angles.
* 
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): 
*  A. Hebert and R. Le Tellier
*
*Parameters: input
* NDIM    number of dimensions for geometry.
* ORDRE   quadrature order.
* NQUAD   number of quadrant (in 3-D) and quarter (in 2-D).
* NBANGL  number of angles.
* DQUAD   relative density of each quadrant.
*
*Parameters: output
* DANGLT  director cosines of angles.
* DDENWT  angular density for each angle.
*
*----------
*
      IMPLICIT   NONE
*----
*  Subroutine arguments
*----
      INTEGER          NDIM,ORDRE,NQUAD,NBANGL
      DOUBLE PRECISION DQUAD(NQUAD),DANGLT(NDIM,NQUAD,NBANGL),
     >                 DDENWT(NQUAD,NBANGL)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTLSN')
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
      INTEGER          MAXORD,MAXNBA,MAXEQ,MAXW
      PARAMETER       (MAXORD=20,MAXNBA=55,MAXEQ=64,MAXW=12)
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI
*----
*  Local variables
*----
      INTEGER          I,IP,IQ,IS,IW,IANG,JOP(MAXORD/2),M2,NPQ,NW,IPR,
     >                 IPL,IPK,INMAX,IPQ,NW0,II,KK,LL,NEQ
      REAL             U(MAXORD/2),TPQ(MAXNBA),UPQ(MAXNBA),
     >                 VPQ(MAXNBA),WPQ(MAXNBA)
      DOUBLE PRECISION FAC,ZMU,ZETA,ZMU1,ZMU2,DDA,X,Y,Z,REF
      INTEGER          INWEI(MAXNBA)
      DOUBLE PRECISION WEI(MAXW),ZMAT(MAXEQ,MAXW+1),UD(MAXW)
*----
*  Set the unique quadrature values
*----
      IF(ORDRE.GT.MAXORD) CALL XABORT(NAMSBR//': MAXORD OVERFLOW.')
      M2=ORDRE/2
      NPQ=M2*(M2+1)/2
      ZMU1=1.0D0/(3.0D0*DBLE(ORDRE-1))
      NW=1+(ORDRE*(ORDRE+8)-1)/48
      IF(NW.GT.MAXW) CALL XABORT('NXTLSN: MAXW OVERFLOW.')
      IF(ORDRE.EQ.2) THEN
         ZMU1=0.33333333
      ELSE IF(ORDRE.EQ.4) THEN
         ZMU1=0.12251480
      ELSE IF(ORDRE.EQ.6) THEN
         ZMU1=0.07109447
      ELSE IF(ORDRE.EQ.8) THEN
         ZMU1=0.04761903
      ELSE IF(ORDRE.EQ.10) THEN
         ZMU1=0.03584310
      ELSE IF(ORDRE.EQ.12) THEN
         ZMU1=0.02796615
      ELSE IF(ORDRE.EQ.14) THEN
         ZMU1=0.02310250
      ELSE IF(ORDRE.EQ.16) THEN
         ZMU1=0.01931398
      ELSE IF(ORDRE.EQ.18) THEN
         ZMU1=0.01692067
      ELSE IF(ORDRE.EQ.20) THEN
         ZMU1=0.01455253
      ELSE
         CALL XABORT(NAMSBR//': ORDER NOT AVAILABLE.')
      ENDIF
      U(1)=REAL(SQRT(ZMU1))
      DO I=2,M2
         ZMU2=ZMU1+2.0D0*DBLE(I-1)*(1.0D0-3.0D0*ZMU1)/DBLE(ORDRE-2)
         U(I)=REAL(SQRT(ZMU2))
      ENDDO
*----
*  Compute the position of weights
*----
      IPR=0
      INMAX=0
      DO IP=1,M2
         JOP(IP)=M2-IP+1
         DO IQ=1,JOP(IP)
            IPR=IPR+1
            IF(IPR.GT.MAXNBA) CALL XABORT('NXTLSN: MAXNBA OVERFLOW.')
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
      IF(INMAX.NE.NW) CALL XABORT(NAMSBR//': INVALID VALUE OD NW.')
      IF(IPR.NE.NPQ) CALL XABORT(NAMSBR//': BAD VALUE ON NPQ.')
*----
*  Set the rectangular system and solve it using the QR method
*----
      NEQ=0
      DO IPL=0,ORDRE,2
        DO IPK=IPL,ORDRE-IPL,2
          IF(MOD(IPL+IPK,2).EQ.1) CYCLE
          NEQ=NEQ+1
          IF(NEQ.GT.MAXEQ) CALL XABORT(NAMSBR//': MAXEQ OVERFLOW.')
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
*  Set the level-symmetric quadratures
*----
      PI=XDRCST('Pi',' ')
      IPQ=0
      DO IP=1,M2
         DO IQ=1,JOP(IP)
            IPQ=IPQ+1
            WPQ(IPQ)=REAL(WEI(INWEI(IPQ))*PI)/2.0
         ENDDO
      ENDDO
*----
*  Fill-in DANGLT and DDENWT array
*----
      FAC=4.0*PI
      DO IANG=1,NBANGL
        X  = DBLE(TPQ(IANG))
        Y  = DBLE(UPQ(IANG))
        Z  = DBLE(VPQ(IANG))
        DDA=DBLE(FAC/WPQ(IANG))
        DANGLT(1,1,IANG)=X
        DANGLT(2,1,IANG)=Y
        DANGLT(3,1,IANG)=Z
        DDENWT(1,IANG)=DQUAD(1)*DDA             
        DANGLT(1,2,IANG)=-X
        DANGLT(2,2,IANG)=Y
        DANGLT(3,2,IANG)=Z
        DDENWT(2,IANG)=DQUAD(2)*DDA 
        DANGLT(1,3,IANG)=X
        DANGLT(2,3,IANG)=-Y
        DANGLT(3,3,IANG)=Z
        DDENWT(3,IANG)=DQUAD(3)*DDA 
        DANGLT(1,4,IANG)=-X
        DANGLT(2,4,IANG)=-Y
        DANGLT(3,4,IANG)=Z
        DDENWT(4,IANG)=DQUAD(4)*DDA 
      ENDDO            
*
      RETURN
      END
