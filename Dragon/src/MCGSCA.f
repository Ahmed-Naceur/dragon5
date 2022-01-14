*DECK MCGSCA
      SUBROUTINE MCGSCA(N,K,M,NOM,NZON,H,XST,B)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculate coefficients of a track for the characteristics integration.
* Step-Characteristics scheme with tabulated exponential calls with
* 'source term isolation' option turned on.
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
*         and step characteristics (SC).
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
      DOUBLE PRECISION H(N),B(N)
*---
* LOCAL VARIABLES
*---
      INTEGER I,NOMI,NZI
      REAL TAU
      DOUBLE PRECISION HID,TAUD
*     tabulated exponential common block
      REAL             E0, E1, PAS1, DX1, XLIM1
      INTEGER          MEX1, LAU
      PARAMETER      ( MEX1=7936 )
      COMMON /EXP1/ E0(0:MEX1),E1(0:MEX1),PAS1,DX1,XLIM1
*
      DO I=2,N-1
         NOMI=NOM(I)
         NZI=NZON(NOMI)
         HID=H(I)
         TAUD=HID*XST(NZI)
         TAU=REAL(TAUD)
         IF(TAU.GE.XLIM1) THEN
*        Out of the table range
            B(I)=1.D0/XST(NZI)
         ELSE
*        Linear interpolation in table of (1-exp(-x))/x
            LAU=INT(TAUD*PAS1)
            B(I)=HID*(E0(LAU)+E1(LAU)*TAUD)
         ENDIF
      ENDDO
*
      RETURN
      END
