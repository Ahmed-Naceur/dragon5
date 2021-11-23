*DECK MCGSCES
      SUBROUTINE MCGSCES(N,K,M,NOM,NZON,H,XST,B)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculate coefficients of a track for the characteristics integration.
* Step-Characteristics scheme with exact exponential calls with
* 'source term isolation' option turned off.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* N       number of elements in the current track.
* K       total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* M       number of material mixtures.
* NOM     vector containing the region number of the different segments
*         of this track.
* NZON    index-number of the mixture type assigned to each volume.
* H       vector containing the lenght of the different segments of this
*         track.
* XST     macroscopic total cross section.
*
*Parameters: output
* B       undefined.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N,K,M,NOM(N),NZON(K)
      REAL XST(0:M)
      DOUBLE PRECISION H(N),B(2,N)
*---
* LOCAL VARIABLES
*---
      INTEGER I,NOMI,NZI
      DOUBLE PRECISION TAUDMIN,TAUD,HID,HID2,TAUD3,TAUD4,TAUD5
*     tabulated exponential common block
      REAL             E0, E1, PAS1, DX1, XLIM1
      INTEGER          MEX1, LAU
      PARAMETER      ( MEX1=7936, TAUDMIN=2.D-2 )
      COMMON /EXP1/ E0(0:MEX1),E1(0:MEX1),PAS1,DX1,XLIM1
*
      DO I=2,N-1
         NOMI=NOM(I)
         NZI=NZON(NOMI)
         HID=H(I)
         TAUD=HID*XST(NZI)
         IF(TAUD.LE.TAUDMIN) THEN  
*        Linear interpolation in table of (1-exp(-x))/x
            LAU=INT(TAUD*PAS1)
            B(1,I)=HID*(E0(LAU)+E1(LAU)*TAUD)
*           and expansion in Taylor serie in O(TAUD^3)
            TAUD3=TAUD/3.D0
            TAUD4=0.125D0*TAUD
            TAUD5=0.2D0*TAUD
            HID2=HID*HID
            B(2,I)=HID2*(0.5D0-TAUD3*(0.5D0-TAUD4*(1.D0-TAUD5)))
         ELSE
*        Exact exponential
            B(1,I)=(1.D0-DEXP(-TAUD))/XST(NZI)
            B(2,I)=(HID-B(1,I))/XST(NZI)
         ENDIF
      ENDDO
*
      RETURN
      END
