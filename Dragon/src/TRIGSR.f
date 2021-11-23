*DECK TRIGSR
      SUBROUTINE TRIGSR (IEX,NGCOND,NMERGE,NALBP,FLXMER,SPH,SUNMER,
     1 IELEM,NBLOS,NUN,NBMIX,SIDE,ZZ,FRZ,MAT,IPERT,KN,QFR,MERG,COUR,
     2 FUNKNO,SOURCE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Source term calculation for a Thomas-Raviart-Schneider formulation of
* the finite element technique in a 3-D hexagonal geometry. Used in
* DRAGON.
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
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* NBLOS   number of lozenges per direction, taking into account
*         mesh-splitting.
* NUN     number of unknowns per group in BIVAC.
* NBMIX   number of macro-mixtures.
* SIDE    side of the hexagons.
* ZZ      Z-directed mesh spacings.
* FRZ     volume fractions for the axial SYME boundary condition.
* MAT     mixture index per lozenge.
* IPERT   mixture permutation index.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
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
      INTEGER IEX,NGCOND,NMERGE,NALBP,IELEM,NBLOS,NUN,NBMIX,
     1 MAT(3,NBLOS),IPERT(NBLOS),MERG(NBMIX),
     2 KN(NBLOS,3+6*(IELEM+2)*IELEM**2)
      REAL FLXMER(NMERGE,NGCOND),SPH(NMERGE+NALBP,NGCOND),
     1 SUNMER(NMERGE,NGCOND),SIDE,ZZ(3,NBLOS),FRZ(NBLOS),
     2 QFR(NBLOS,8),COUR,FUNKNO(NUN,NGCOND),SOURCE(NUN)
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
         DO 16 K4=0,1
         DO 10 K2=1,IELEM+1,IELEM
         KNW1=KN(NUM,3+(K4+0)*NELEM*IELEM+K2)
         KNX1=KN(NUM,3+(K4+2)*NELEM*IELEM+K2)
         KNY1=KN(NUM,3+(K4+4)*NELEM*IELEM+K2)
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
         DO 15 K2=1,6,5
         KNZ1=KN(NUM,3+(K4+6)*NELEM*IELEM+K2)
         IF((KNZ1.NE.0).AND.(K2.EQ.1).AND.(K4.EQ.0)) THEN
            IND1=ABS(KNZ1)
            SG=REAL(SIGN(1,KNZ1))
            SOURCE(IND1)=SOURCE(IND1)-SG*0.5*COUR*QFR(NUM,7)
         ENDIF
         IF((KNZ1.NE.0).AND.(K2.EQ.6).AND.(K4.EQ.1)) THEN
            IND1=ABS(KNZ1)
            SG=REAL(SIGN(1,KNZ1))
            SOURCE(IND1)=SOURCE(IND1)+SG*0.5*COUR*QFR(NUM,8)
         ENDIF
   15    CONTINUE
   16    CONTINUE
   20    CONTINUE      
      ENDIF
*----
*  DIFFUSION AND FISSION SOURCE.
*----
      IF(IEX.EQ.1) GO TO 60
      DO 55 JGR=1,NGCOND
      NUM=0
      DO 40 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 40
      NUM=NUM+1
      L=MAT(1,IPERT(KEL))
      IF(L.EQ.0) GO TO 40
      VOL0=TTTT*ZZ(1,IPERT(KEL))*FRZ(KEL)
      GARS=SUNMER(MERG(L),JGR)*SPH(MERG(L),JGR)
      DO 32 K3=0,IELEM-1
      DO 31 K2=0,IELEM-1
      DO 30 K1=0,IELEM-1
      JND1=(NUM-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
      JND2=(KN(NUM,1)-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
      JND3=(KN(NUM,2)-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
      SOURCE(JND1)=SOURCE(JND1)+FUNKNO(JND1,JGR)*VOL0*GARS
      SOURCE(JND2)=SOURCE(JND2)+FUNKNO(JND2,JGR)*VOL0*GARS
      SOURCE(JND3)=SOURCE(JND3)+FUNKNO(JND3,JGR)*VOL0*GARS
   30 CONTINUE
   31 CONTINUE
   32 CONTINUE
   40 CONTINUE
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
      VOL0=TTTT*ZZ(1,IPERT(KEL))*FRZ(KEL)
      JND1=(NUM-1)*IELEM**3+1
      JND2=(KN(NUM,1)-1)*IELEM**3+1
      JND3=(KN(NUM,2)-1)*IELEM**3+1
      PV=0.0
      DO 70 JGR=1,NGCOND
      PV=PV+SUNMER(MERG(L),JGR)*FLXMER(MERG(L),JGR)
   70 CONTINUE
      SOURCE(JND1)=SOURCE(JND1)+VOL0*PV
      SOURCE(JND2)=SOURCE(JND2)+VOL0*PV
      SOURCE(JND3)=SOURCE(JND3)+VOL0*PV
   80 CONTINUE
      RETURN
      END
