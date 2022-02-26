*DECK MCGDS2E
      SUBROUTINE MCGDS2E(N,M,NFI,NOM,NZON,H,XST,XSW,DINV,B,A)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of ACA coefficients for this track
* (exact exponentials version).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s) : R. Le Tellier
*
*Parameters: input
* N       number of elements for this track.
* M       number of material mixtures.
* NFI     total number of volumes and surfaces for which specific values
*         of the neutron flux and reactions rates are required.
* NOM     integer tracking elements.
* NZON    zone number.
* H       tracking widths.
* XST     total cross sections array.
* XSW     scattering cross sections array.
*
*Parameters: output
* DINV    undefined.
* B       undefined.
* A       undefined.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*---
* SUBROUTINE ARGUMENTS
*---
      INTEGER N,M,NFI,NOM(N),NZON(NFI)
      REAL XST(0:M),XSW(0:M)
      DOUBLE PRECISION H(N),DINV(N),B(N),A(N)
*---
* LOCAL VARIABLES
*---
      DOUBLE PRECISION TAUDMIN
      PARAMETER(TAUDMIN=1.D-3)
      INTEGER I,NOMI,NZI
      DOUBLE PRECISION TAUD,ALPHA,TEMP
*
      DO I=1,N
         NOMI=NOM(I)
         NZI=NZON(NOMI)
         IF (NZI.GE.0) THEN
            TAUD=H(I)*DBLE(XST(NZI))
            IF (TAUD.GT.TAUDMIN) THEN
               TEMP=1.D0-DEXP(-TAUD)
               ALPHA=2.D0/TEMP-2.D0/TAUD-1.D0
               DINV(I)=TAUD/(2.D0+TAUD*ALPHA)
               B(I)=0.5D0*TAUD*(DINV(I)-ALPHA)
               A(I)=1.D0-B(I)
               B(I)=B(I)/XST(NZI)
               A(I)=A(I)+B(I)*XSW(NZI)
            ELSE
               DINV(I)=0.5D0*TAUD
               B(I)=0.D0
               A(I)=1.D0
            ENDIF
         ENDIF
      ENDDO
*
      RETURN
      END
