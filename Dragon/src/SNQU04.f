*DECK SNQU04
      SUBROUTINE SNQU04(NLF,JOP,U,W,TPQ,UPQ,VPQ,WPQ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Set the Gauss-Legendre, Gauss-Chebyshev quadratures.
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
      PARAMETER(PI=3.141592654,MAXNLF=64)
      REAL U2(MAXNLF),W2(MAXNLF)
*----
*  SET THE QUARATURE VALUES.
*----
      M2=NLF/2
      NPQ=M2*(M2+1)/2
      IF(NLF.EQ.2) THEN
      	U(1)=1/SQRT(3.0)
      	W(1)=1.0
      ELSE
        IF(NLF.GT.MAXNLF) CALL XABORT('SNQU04: TOO MANY GAUSS POINTS.')
        CALL ALGPT(NLF,-1.0,1.0,U2,W2)
        DO 45 M=1,M2
          U(M)=U2(M2+M)
          W(M)=W2(M2+M)
   45   CONTINUE
      ENDIF
      IPQ=0
      WSUM=0.0
      DO IP=1,M2
         JOP(IP)=M2-IP+1
         DO IQ=1,JOP(IP)
            IPQ=IPQ+1
            OMEGA=0.5*PI*(1.0-REAL(NLF-2*IP-2*IQ+3)/REAL(NLF-2*IP+2))
            TPQ(IPQ)=U(IP)
            UPQ(IPQ)=SQRT(1.0-U(IP)*U(IP))*COS(OMEGA)
            VPQ(IPQ)=SQRT(1.0-U(IP)*U(IP))*SIN(OMEGA)
            WPQ(IPQ)=PI*W(IP)/REAL(NLF-2*IP+2)
            WSUM=WSUM+WPQ(IPQ)
         ENDDO
      ENDDO
      IF(IPQ.NE.NPQ) CALL XABORT('SNQU04: BAD VALUE ON NPQ.')
      RETURN
      END
