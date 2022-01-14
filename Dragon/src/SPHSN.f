*DECK SPHSN
      SUBROUTINE SPHSN(IPTRK2,NPSYS,NREG,NUN,NMERGE,NALBP,NGCOND,NBMIX,
     1 MAT,VOL,KEY,MERG,SUNMER,FLXMER,SPH,IEX,FUNKNO,SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* SN source calculation over the macro-geometry.
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
* IPTRK2  pointer to the tracking of the macro-geometry
*         (L_TRACK signature).
* NPSYS   group masks.
* NREG    number of macro-regions (in the macro calculation).
* NUN     number of unknowns in the macro-calculation.
* NMERGE  number of merged regions.
* NALBP   number of physical albedos.
* NGCOND  number of condensed groups.
* NBMIX   number of macro-mixtures.
* MAT     mixture index per macro-region.
* VOL     volume of macro-regions.
* KEY     position of the flux components associated with each volume.
* MERG    index of merged macro-regions per macro-mixture.
* SUNMER  incoming source (scattering+fission) cross sections.
* FLXMER  flux estimate per mixture.
* SPH     SPH factors.
* IEX     iteration index.
* FUNKNO  neutron flux.
*
*Parameters: output
* SUNKNO  neutron sources.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK2
      INTEGER NPSYS(NGCOND),NREG,NUN,NMERGE,NALBP,NGCOND,NBMIX,
     1 MAT(NREG),KEY(NREG),MERG(NBMIX),IEX
      REAL VOL(NREG),SUNMER(NMERGE,NGCOND,NGCOND),FLXMER(NMERGE,NGCOND),
     1 SPH(NMERGE+NALBP,NGCOND),FUNKNO(NUN,NGCOND),SUNKNO(NUN,NGCOND)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (NSTATE=40)
      INTEGER     IPAR(NSTATE)
*----
*  RECOVER SN SPECIFIC PARAMETERS.
*----
      CALL LCMGET(IPTRK2,'STATE-VECTOR',IPAR)
      ISCAT=IPAR(16)
      NSCT=ISCAT*(ISCAT+1)/2
*----
*  COMPUTE THE NEUTRON SOURCES IN GROUP IGR
*----
      DO 20 IGR=1,NGCOND
      CALL XDRSET(SUNKNO(1,IGR),NUN,0.0)
      IF(NPSYS(IGR).EQ.0) GO TO 20
      COUR=0.0
      DO 10 I=1,NREG
      IF(MAT(I).GT.NBMIX) CALL XABORT('SPHSN: NBMIX OVERFLOW.')
   10 CONTINUE
      CALL SPHSNO(IEX,NGCOND,IGR,NMERGE,NALBP,FLXMER,SPH,
     1 SUNMER(1,1,IGR),NREG,NUN,NBMIX,NSCT,MAT,VOL,KEY,MERG,COUR,
     2 FUNKNO,SUNKNO(1,IGR))
   20 CONTINUE
      RETURN
      END
