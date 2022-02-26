*DECK SNQU03
      SUBROUTINE SNQU03(NLF,JOP,U,W,TPQ,UPQ,VPQ,WPQ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Set the level-symmetric quadratures of code SNOW.
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
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NLF,JOP(NLF/2)
      REAL U(NLF/2),W(NLF/2),TPQ(NLF*(NLF/2+1)/4),UPQ(NLF*(NLF/2+1)/4),
     1 VPQ(NLF*(NLF/2+1)/4),WPQ(NLF*(NLF/2+1)/4)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(PI=3.141592654,MAXNLF=20,MAXNBA=55,MAXW=12)
      INTEGER INWEI(MAXNBA)
      DOUBLE PRECISION ZMU1,ZMU2,WSUM2,WEI(MAXW)
*----
*  SET THE UNIQUE QUADRATURE VALUES.
*----
      IF(NLF.GT.MAXNLF) CALL XABORT('SNQU03: MAXNLF OVERFLOW.')
      M2=NLF/2
      NPQ=M2*(M2+1)/2
      ZMU1=1.0D0/(3.0D0*DBLE(NLF-1))
      NW=1+(NLF*(NLF+8)-1)/48
      IF(NW.GT.MAXW) CALL XABORT('SNQU03: MAXW OVERFLOW.')
      IF(NLF.EQ.2) THEN
         WEI(1)=1.0D0
      ELSE IF(NLF.EQ.4) THEN
         WEI(1)=0.3333333
      ELSE IF(NLF.EQ.6) THEN
         WEI(1)=0.1608612
         WEI(2)=0.1724721
      ELSE IF(NLF.EQ.8) THEN
         WEI(1)=0.1066008
         WEI(2)=0.1011727
         WEI(3)=0.0731616
      ELSE
         CALL XABORT('SNQU03: ORDER NOT AVAILABLE.')
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
            IF(IPR.GT.MAXNBA) CALL XABORT('SNQU03: MAXNBA OVERFLOW.')
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
      IF(INMAX.NE.NW) CALL XABORT('SNQU03: INVALID VALUE OD NW.')
      IF(IPR.NE.NPQ) CALL XABORT('SNQU03: BAD VALUE ON NPQ.')
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
