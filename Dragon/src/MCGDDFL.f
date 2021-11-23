*DECK MCGDDFL
      SUBROUTINE MCGDDFL(N,K,M,NOM,NZON,H,XST,B)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculate coefficients of a track for the characteristics 
* integration.
* DD1 Diamond-Differencing scheme without fix-up with
* 'source term isolation' option turned off.
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
      DOUBLE PRECISION H(N),B(0:5,N)
*---
* LOCAL VARIABLES
*---
      INTEGER I,NOMI,NZI
      DOUBLE PRECISION SQ3,TAUD,HID,H2,H3,DEN
*
      SQ3=SQRT(3.0D0)
      DO I=2,N-1
         NOMI=NOM(I)
         NZI=NZON(NOMI)
         HID=H(I)
         TAUD=HID*XST(NZI)
         H2=HID*HID
         H3=H2*HID
         DEN=(TAUD+6.0D0)*TAUD+12.0D0
         B(0,I)=((TAUD-6.0D0)*TAUD+12.0D0)/DEN
         B(1,I)=12.0D0*HID/DEN
         B(2,I)=H2*(TAUD+6.0D0)/DEN
         B(3,I)=-2.0D0*SQ3*H3/DEN
         B(4,I)=-2.0D0*SQ3*TAUD/DEN
         B(5,I)=H3/DEN
      ENDDO
*
      RETURN
      END
