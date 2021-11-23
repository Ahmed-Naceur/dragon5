*DECK MCGDDFS
      SUBROUTINE MCGDDFS(N,K,M,NOM,NZON,H,XST,B)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculate coefficients of a track for the characteristics integration.
* Diamond-Differencing scheme without fix-up with
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
* B       DD1 coefficients.
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
      DOUBLE PRECISION TAUD,HID,TEMP
*
      DO I=2,N-1
         NOMI=NOM(I)
         NZI=NZON(NOMI)
         HID=H(I)
         TAUD=HID*XST(NZI)
         TEMP=HID/(2.D0+TAUD)
         B(1,I)=TEMP+TEMP
         B(2,I)=HID*TEMP
      ENDDO
*
      RETURN
      END
