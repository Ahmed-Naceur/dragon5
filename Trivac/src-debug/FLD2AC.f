*DECK FLD2AC
      SUBROUTINE FLD2AC(NG,NUN,IG0,FLUX,ZMU)
*
*-----------------------------------------------------------------------
*
*Purpose:
* one-factor variationnal acceleration of the flux.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* NG      number of energy groups.
* NUN     number of unknowns per energy group.
* IG0     first group to accelerate.
*
*Parameters: input/output
* FLUX    neutron flux:
*         FLUX(:,:,1) <=old;
*         FLUX(:,:,2) <=present;
*         FLUX(:,:,3) <=new.
*
*Parameters: output
* ZMU     acceleration factor.
*
*-----------------------------------------------------------------------
*
      IMPLICIT  NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   NG, NUN, IG0
      REAL      FLUX(NUN,NG,3), ZMU
*----
*  LOCAL VARIABLES
*----
      INTEGER           IG, IR
      DOUBLE PRECISION  DMU, R1, R2
      DOUBLE PRECISION  DONE, DZERO, NOM, DENOM
      PARAMETER ( DONE=1.0D0, DZERO=0.0D0 )
*----
*  ZMU CALCULATION
*----
      NOM   = DZERO
      DENOM = DZERO
      DO  3 IG= IG0,NG
          DO  2 IR=1,NUN
             R1 = FLUX(IR,IG,2) - FLUX(IR,IG,1)
             R2 = FLUX(IR,IG,3) - FLUX(IR,IG,2)
             NOM = NOM + R1*(R2-R1)
             DENOM = DENOM + (R2-R1)*(R2-R1)
    2     CONTINUE
    3 CONTINUE
*
      DMU = - NOM / DENOM
      ZMU  = REAL(DMU)
      IF( DMU.GT.DZERO )THEN
       DO  13 IG= IG0,NG
          DO  12 IR=1,NUN
*
*           ACCELERATED VALUES FOR PHI(2) ET PHI(3)
            FLUX(IR,IG,3) = FLUX(IR,IG,2) + REAL(DMU) *
     >        (FLUX(IR,IG,3) - FLUX(IR,IG,2))
            FLUX(IR,IG,2) = FLUX(IR,IG,1) + REAL(DMU) *
     >        (FLUX(IR,IG,2) - FLUX(IR,IG,1))
   12     CONTINUE
   13  CONTINUE
      ELSE
       ZMU= 1.0
      ENDIF
      RETURN
      END
