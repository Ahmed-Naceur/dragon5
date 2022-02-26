*DECK MCGSCEL
      SUBROUTINE MCGSCEL(N,K,M,NOM,NZON,H,XST,B)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculate coefficients of a track for the characteristics 
* integration. Linear-Discontinuous-Characteristics scheme with 
* exact exponential calls. Source term isolation option turned off.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
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
* B       LDC coefficients.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N,K,M,NOM(N),NZON(K)
      REAL XST(0:M)
      DOUBLE PRECISION H(N),B(0:5,N)
*---
* LOCAL VARIABLES
*---
      INTEGER I,NOMI,NZI
      DOUBLE PRECISION TAUDMIN,SQ3,TAUD,HID,TEMP,DSIG,C1,C2,H2,H3,TEMP1
*     tabulated exponential common block
      PARAMETER      ( TAUDMIN=2.0D-2)
*
      SQ3=SQRT(3.0D0)
      DO I=2,N-1
         NOMI=NOM(I)
         NZI=NZON(NOMI)
         HID=H(I)
         TAUD=HID*XST(NZI)
         IF(TAUD.LE.TAUDMIN) THEN 
*        Use Taylor series expansions
            H2=HID*HID
            H3=H2*HID
            B(0,I)=TAUD*(0.5D0*TAUD-1.0D0)+1.0D0
            B(1,I)=HID*(TAUD*(TAUD/6.0D0-0.5D0)+1.0D0)
            B(2,I)=H2*(TAUD*(TAUD-4.0D0)+12.0D0)/24.0D0
            B(3,I)=-SQ3*H3*(TAUD*(TAUD-2.0D0)+4.0D0)/24.0D0
            B(4,I)=-SQ3*TAUD*(TAUD*(TAUD-2.0D0)+4.0D0)/24.0D0
            B(5,I)=H3*(TAUD*TAUD-TAUD+4.0D0)/40.0D0
         ELSE
*        Use exact exponential
            TEMP1=EXP(-TAUD)
            DSIG=DBLE(XST(NZI))
            B(0,I)=TEMP1
            TEMP=(1.D0-TEMP1)/TAUD
            B(1,I)=TEMP*HID
            B(2,I)=HID*(1.D0-TEMP)/DSIG
            B(3,I)=-SQ3*HID*(2.0D0-(TAUD+2.0D0)*TEMP)/DSIG**2
            B(4,I)=-SQ3*(2.0D0-(TAUD+2.0D0)*TEMP)/TAUD
            C1=TAUD*(TAUD-6.0D0)-12.0D0
            C2=TAUD*(3.0D0*TAUD+12.0D0)+12.0D0 
            B(5,I)=(C1+C2*TEMP)/(TAUD*DSIG**3)
         ENDIF
      ENDDO
*
      RETURN
      END
