*DECK TRIGSO
      SUBROUTINE TRIGSO (IEX,NGCOND,NMERGE,NALBP,FLXMER,SPH,SUNMER,
     1 CYLIND,IELEM,NREG,NUN,NBMIX,MAT,XX,DD,KN,QFR,VOL,MERG,COUR,
     2 FUNKNO,SOURCE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Source term calculation for a mixed-dual formulation of the finite
* element technique in a 3-D Cartesian geometry. Used in DRAGON.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
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
* NMERGE  number of merged regions.
* NALBP   number of physical albedos.
* FLXMER  flux estimate per mixture.
* SPH     SPH factors.
* SUNMER  incoming source (scattering+fission) cross sections.
* CYLIND  cylinderization flag (=.TRUE. for cylindrical geometry)
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* NREG    number of elements in BIVAC.
* NUN     number of unknowns per group in BIVAC.
* NBMIX   number of macro-mixtures.
* MAT     mixture index per region.
* XX      X-directed mesh spacings.
* DD      values used with a cylindrical geometry.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* VOL     volume of regions.
* MERG    index of merged regions per macro-mixture.
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
      INTEGER IEX,NGCOND,NMERGE,NALBP,IELEM,NREG,NUN,NBMIX,MAT(NREG),
     1 KN(NREG*(1+6*IELEM**2)),MERG(NBMIX)
      REAL FLXMER(NMERGE,NGCOND),SPH(NMERGE+NALBP,NGCOND),
     1 SUNMER(NMERGE,NGCOND),XX(NREG),DD(NREG),QFR(4*NREG),VOL(NREG),
     2 COUR,FUNKNO(NUN,NGCOND),SOURCE(NUN)
      LOGICAL CYLIND
*
      CALL XDRSET(SOURCE,NUN,0.0)
*----
*  INCOMING CURRENT SOURCE.
*----
      IF(COUR.NE.0.0) THEN
         NUM1=0
         NUM2=0
         DO 20 K=1,NREG
         IF(MAT(K).EQ.0) GO TO 20
         IF(VOL(K).EQ.0.0) GO TO 10
         IF(CYLIND) THEN
            DIN=1.0-0.5*XX(K)/DD(K)
            DOT=1.0+0.5*XX(K)/DD(K)
         ELSE
            DIN=1.0
            DOT=1.0
         ENDIF
*
         IF(KN(NUM1+2).NE.0) THEN
            IND1=ABS(KN(NUM1+2))
            SG=REAL(SIGN(1,KN(NUM1+2)))
            SOURCE(IND1)=SOURCE(IND1)-SG*0.5*COUR*DIN*QFR(NUM2+1)
         ENDIF
         IF(KN(NUM1+3).NE.0) THEN
            IND1=ABS(KN(NUM1+2+IELEM**2))
            SG=REAL(SIGN(1,KN(NUM1+2+IELEM**2)))
            SOURCE(IND1)=SOURCE(IND1)+SG*0.5*COUR*DOT*QFR(NUM2+2)
         ENDIF
         IF(KN(NUM1+4).NE.0) THEN
            IND1=ABS(KN(NUM1+2+2*IELEM**2))
            SG=REAL(SIGN(1,KN(NUM1+2+2*IELEM**2)))
            SOURCE(IND1)=SOURCE(IND1)-SG*0.5*COUR*QFR(NUM2+3)
         ENDIF
         IF(KN(NUM1+5).NE.0) THEN
            IND1=ABS(KN(NUM1+2+3*IELEM**2))
            SG=REAL(SIGN(1,KN(NUM1+2+3*IELEM**2)))
            SOURCE(IND1)=SOURCE(IND1)+SG*0.5*COUR*QFR(NUM2+4)
         ENDIF
         IF(KN(NUM1+6).NE.0) THEN
            IND1=ABS(KN(NUM1+2+4*IELEM**2))
            SG=REAL(SIGN(1,KN(NUM1+2+4*IELEM**2)))
            SOURCE(IND1)=SOURCE(IND1)-SG*0.5*COUR*QFR(NUM2+5)
         ENDIF
         IF(KN(NUM1+7).NE.0) THEN
            IND1=ABS(KN(NUM1+2+5*IELEM**2))
            SG=REAL(SIGN(1,KN(NUM1+2+5*IELEM**2)))
            SOURCE(IND1)=SOURCE(IND1)+SG*0.5*COUR*QFR(NUM2+6)
         ENDIF
   10    NUM1=NUM1+1+6*IELEM**2
         NUM2=NUM2+6
   20    CONTINUE
      ENDIF
*----
*  DIFFUSION AND FISSION SOURCE.
*----
      IF(IEX.EQ.1) GO TO 80
      DO 70 JGR=1,NGCOND
      NUM1=0
      DO 50 K=1,NREG
      L=MAT(K)
      IF(L.EQ.0) GO TO 50
      IF(VOL(K).EQ.0.0) GO TO 40
      GARS=SUNMER(MERG(L),JGR)*SPH(MERG(L),JGR)
      DO 30 I0=1,IELEM**3
      IND1=KN(NUM1+1)+I0-1
      SOURCE(IND1)=SOURCE(IND1)+FUNKNO(IND1,JGR)*VOL(K)*GARS
   30 CONTINUE
   40 NUM1=NUM1+1+6*IELEM**2
   50 CONTINUE
   70 CONTINUE
      RETURN
*----
*  INITIAL ITERATION.
*----
   80 NUM1=0
      DO 140 K=1,NREG
      L=MAT(K)
      IF(L.EQ.0) GO TO 140
      IF(VOL(K).EQ.0.0) GO TO 130
      IND1=KN(NUM1+1)
      PV=0.0
      DO 120 JGR=1,NGCOND
      PV=PV+SUNMER(MERG(L),JGR)*FLXMER(MERG(L),JGR)
  120 CONTINUE
      SOURCE(IND1)=SOURCE(IND1)+VOL(K)*PV
  130 NUM1=NUM1+1+6*IELEM**2
  140 CONTINUE
      RETURN
      END
