*DECK SPHSNO
      SUBROUTINE SPHSNO (IEX,NGCOND,IGR,NMERGE,NALBP,FLXMER,SPH,SUNMER,
     1 NREG,NUN,NBMIX,NSCT,MAT,VOL,KEY,MERG,COUR,FUNKNO,SOURCE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* SPH source term calculation for a PIJ, MOC or SN discretization.
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
* IEX     iteration number.
* NGCOND  number of groups condensed.
* IGR     energy group index.
* NMERGE  number of merged regions.
* NALBP   number of physical albedos.
* FLXMER  flux estimate per mixture.
* SPH     SPH factors.
* SUNMER  incoming source (scattering+fission) cross sections.
* NREG    number of regions in macro-calculation.
* NUN     number of unknowns per group in macro-calculation.
* NBMIX   number of macro-mixtures.
* NSCT    maximum number of spherical harmonics moments of the flux.
* MAT     mixture index per region.
* VOL     volume of regions.
* KEY     position of the flux components associated with each volume.
* MERG    index of merged regions.
* COUR    four times the incoming current per unit surface.
* FUNKNO  previously calculated fluxes.
*
*Parameters: output
* SOURCE  fission and diffusion sources.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IEX,NGCOND,IGR,NMERGE,NALBP,NREG,NUN,NBMIX,MAT(NREG),
     1 KEY(NREG),MERG(NBMIX)
      REAL FLXMER(NMERGE,NGCOND),SPH(NMERGE+NALBP,NGCOND),
     1 SUNMER(NMERGE,NGCOND),VOL(NREG),COUR,FUNKNO(NUN,NGCOND),
     2 SOURCE(NUN)
*----
*  INCOMING CURRENT SOURCE.
*----
      IF(COUR.NE.0.0) THEN
         CALL XABORT('SPHSNO: INCOMING CURRENT NOT IMPLEMENTED.')
      ENDIF
*----
*  SCATTERING AND FISSION SOURCE.
*----
      CALL XDRSET(SOURCE,NUN,0.0)
      IGR2=NGCOND
      IF(IEX.EQ.1) IGR2=IGR-1
      DO 115 JGR=1,IGR2
      DO 110 K=1,NREG
      L=MAT(K)
      IF(L.EQ.0) GO TO 110
      IF(VOL(K).EQ.0.0) GO TO 110
      GARS=SUNMER(MERG(L),JGR)*SPH(MERG(L),JGR)
      DO 90 IL=1,NSCT
      IND1=(K-1)*NSCT+IL
      SOURCE(IND1)=SOURCE(IND1)+FUNKNO(IND1,JGR)*GARS
   90 CONTINUE
  110 CONTINUE
  115 CONTINUE
      IF(IEX.GT.1) RETURN
*----
*  INITIAL ITERATION.
*----
      DO 140 K=1,NREG
      L=MAT(K)
      IF(L.EQ.0) GO TO 140
      IF(VOL(K).EQ.0.0) GO TO 140
      IND1=KEY(K)
      PV=0.0
      DO 120 JGR=IGR,NGCOND
      PV=PV+SUNMER(MERG(L),JGR)*FLXMER(MERG(L),JGR)
  120 CONTINUE
      SOURCE(IND1)=SOURCE(IND1)+PV
  140 CONTINUE
      RETURN
      END
