*DECK BIVGSR
      SUBROUTINE BIVGSR (IEX,NGCOND,NMERGE,NALBP,FLXMER,SPH,SUNMER,
     1 IELEM,NBLOS,NUN,NBMIX,SIDE,MAT,IPERT,KN,QFR,MERG,COUR,FUNKNO,
     2 SOURCE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Source term calculation for a Thomas-Raviart-Schneider formulation of
* the finite element technique in a 2-D hexagonal geometry.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
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
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* NBLOS   number of lozenges per direction, taking into account
*         mesh-splitting.
* NUN     number of unknowns per group in BIVAC.
* NBMIX   number of macro-mixtures.
* SIDE    side of the hexagons.
* MAT     mixture index per region.
* IPERT   mixture permutation index.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* MERG    index of merged regions prt macro-mixture.
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
      INTEGER IEX,NGCOND,NMERGE,NALBP,IELEM,NBLOS,NUN,NBMIX,
     1 MAT(3,NBLOS),IPERT(NBLOS),KN(NBLOS,4+6*IELEM*(IELEM+1)),
     2 MERG(NBMIX)
      REAL FLXMER(NMERGE,NGCOND),SPH(NMERGE+NALBP,NGCOND),
     1 SUNMER(NMERGE,NGCOND),SIDE,QFR(NBLOS,6),COUR,
     2 FUNKNO(NUN,NGCOND),SOURCE(NUN)
*
      CALL XDRSET(SOURCE,NUN,0.0)
      NELEM=IELEM*(IELEM+1)
      TTTT=0.5*SQRT(3.0)*SIDE*SIDE
*----
*  INCOMING CURRENT SOURCE.
*----
      IF(COUR.NE.0.0) THEN
         NUM=0
         DO 20 KEL=1,NBLOS
         IF(IPERT(KEL).EQ.0) GO TO 20
         NUM=NUM+1
         DO 15 K4=0,1
         DO 10 K2=1,IELEM+1,IELEM
         KNW1=KN(NUM,4+K4*NELEM+K2)
         KNX1=KN(NUM,4+(K4+2)*NELEM+K2)
         KNY1=KN(NUM,4+(K4+4)*NELEM+K2)
         IF((KNW1.NE.0).AND.(K2.EQ.1).AND.(K4.EQ.0)) THEN
            IND1=ABS(KNW1)
            SG=REAL(SIGN(1,KNW1))
            SOURCE(IND1)=SOURCE(IND1)-SG*0.5*COUR*QFR(NUM,1)
         ENDIF
         IF((KNW1.NE.0).AND.(K2.EQ.IELEM+1).AND.(K4.EQ.1)) THEN
            IND1=ABS(KNW1)
            SG=REAL(SIGN(1,KNW1))
            SOURCE(IND1)=SOURCE(IND1)+SG*0.5*COUR*QFR(NUM,2)
         ENDIF
         IF((KNX1.NE.0).AND.(K2.EQ.1).AND.(K4.EQ.0)) THEN
            IND1=ABS(KNX1)
            SG=REAL(SIGN(1,KNX1))
            SOURCE(IND1)=SOURCE(IND1)-SG*0.5*COUR*QFR(NUM,3)
         ENDIF
         IF((KNX1.NE.0).AND.(K2.EQ.IELEM+1).AND.(K4.EQ.1)) THEN
            IND1=ABS(KNW1)
            SG=REAL(SIGN(1,KNW1))
            SOURCE(IND1)=SOURCE(IND1)+SG*0.5*COUR*QFR(NUM,4)
         ENDIF
         IF((KNY1.NE.0).AND.(K2.EQ.1).AND.(K4.EQ.0)) THEN
            IND1=ABS(KNY1)
            SG=REAL(SIGN(1,KNY1))
            SOURCE(IND1)=SOURCE(IND1)-SG*0.5*COUR*QFR(NUM,5)
         ENDIF
         IF((KNY1.NE.0).AND.(K2.EQ.IELEM+1).AND.(K4.EQ.1)) THEN
            IND1=ABS(KNY1)
            SG=REAL(SIGN(1,KNY1))
            SOURCE(IND1)=SOURCE(IND1)+SG*0.5*COUR*QFR(NUM,6)
         ENDIF
   10    CONTINUE
   15    CONTINUE
   20    CONTINUE      
      ENDIF
*----
*  DIFFUSION AND FISSION SOURCE.
*----
      IF(IEX.EQ.1) GO TO 60
      DO 55 JGR=1,NGCOND
      NUM=0
      DO 50 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 50
      NUM=NUM+1
      L=MAT(1,IPERT(KEL))
      IF(L.EQ.0) GO TO 50
      GARS=SUNMER(MERG(L),JGR)*SPH(MERG(L),JGR)
      DO 45 K2=0,IELEM-1
      DO 40 K1=0,IELEM-1
      JND1=KN(NUM,1)+K2*IELEM+K1
      JND2=KN(NUM,2)+K2*IELEM+K1
      JND3=KN(NUM,3)+K2*IELEM+K1
      SOURCE(JND1)=SOURCE(JND1)+FUNKNO(JND1,JGR)*TTTT*GARS
      SOURCE(JND2)=SOURCE(JND2)+FUNKNO(JND2,JGR)*TTTT*GARS
      SOURCE(JND3)=SOURCE(JND3)+FUNKNO(JND3,JGR)*TTTT*GARS
   40 CONTINUE
   45 CONTINUE
   50 CONTINUE
   55 CONTINUE
      RETURN
*----
*  INITIAL ITERATION.
*----
   60 NUM=0
      DO 80 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 80
      NUM=NUM+1
      L=MAT(1,IPERT(KEL))
      IF(L.EQ.0) GO TO 80
      JND1=KN(NUM,1)
      JND2=KN(NUM,2)
      JND3=KN(NUM,3)
      PV=0.0
      DO 70 JGR=1,NGCOND
      PV=PV+SUNMER(MERG(L),JGR)*FLXMER(MERG(L),JGR)
   70 CONTINUE
      SOURCE(JND1)=SOURCE(JND1)+TTTT*PV
      SOURCE(JND2)=SOURCE(JND2)+TTTT*PV
      SOURCE(JND3)=SOURCE(JND3)+TTTT*PV
   80 CONTINUE
      RETURN
      END
