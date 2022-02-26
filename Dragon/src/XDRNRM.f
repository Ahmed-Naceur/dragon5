*DECK XDRNRM
      SUBROUTINE XDRNRM(NREGIO,NBMIX,MATCOD,VOLUME,XSSIGT,XSSIGW,
     > PIJSYM,PIS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Normalisation of the scattering-reduced cp matrix to force neutron
* conservation (no leakage).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NREGIO  number of regions considered.
* NBMIX   number of mixtures.
* MATCOD  mixture code in each region.
* VOLUME  volume of each region.
* XSSIGT  total macroscopic cross sections.
* XSSIGW  P0 within-group scattering cross sections.
*
*Parameters: input/output
* PIJSYM  group condensed reduce/symmetric scattering-reduced pij
*         matrix.
* PIS     escape probabilities.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER    NREGIO,NBMIX,MATCOD(NREGIO)
      REAL       VOLUME(NREGIO),XSSIGT(NBMIX),XSSIGW(NBMIX),
     >           PIJSYM(NREGIO*(NREGIO+1)/2),PIS(NREGIO)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION SUM,DD
*----
*  INTRINSIC FUNCTION FOR POSITION IN CONDENSED PIJ MATRIX
*----
      INDPOS(I,J)=MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
*
      DD=0.0D0
      DO 20 I=1,NREGIO
      SUM=0.0D0
      DO 10 J=1,NREGIO
      MATNUM=MATCOD(J)
      IF(MATNUM.GT.0) THEN
         SUM=SUM+(XSSIGT(MATNUM)-XSSIGW(MATNUM))*PIJSYM(INDPOS(I,J))/
     1   VOLUME(I)
      ENDIF
   10 CONTINUE
      PIS(I)=REAL(1.0D0-SUM)
      MATNUM=MATCOD(I)
      IF(MATNUM.GT.0) THEN
         DD=DD+(XSSIGT(MATNUM)-XSSIGW(MATNUM))*VOLUME(I)*PIS(I)
      ENDIF
   20 CONTINUE
      DO 40 I=1,NREGIO
      DO 30 J=1,I
      INDPIJ=INDPOS(I,J)
      PIJSYM(INDPIJ)=PIJSYM(INDPIJ)+PIS(I)*PIS(J)*VOLUME(I)*VOLUME(J)/
     > REAL(DD)
   30 CONTINUE
   40 CONTINUE
      RETURN
      END
