*DECK BIVFSH
      SUBROUTINE BIVFSH (IEX,MAXKN,MAXQF,NGCOND,NMERGE,NALBP,FLXMER,
     1 SPH,SUNMER,NREG,NUN,ISPLH,NELEM,NBMIX,SIDE,MAT,KN,QFR,VOL,MERG,
     2 COUR,FUNKNO,T,RH,RT,SOURCE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Source term calculation for finite element or mesh corner finite
* differences in hexagonal geometry. 
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
* IEX     iteration number.
* MAXKN   dimension of array KN.
* MAXQF   dimension of array QFR.
* NGCOND  number of groups condensed.
* NMERGE  number of mixtures.
* NALBP   number of physical albedos.
* FLXMER  flux estimate per mixture.
* SPH     SPH factors.
* SUNMER  incoming source (scattering+fission) cross sections.
* NREG    number of hexagons in BIVAC.
* NUN     number of unknown per group in BIVAC.
* ISPLH   type of hexagonal mesh-splitting:
*         =1: hexagonal elements; >1: triangular elements.
* NELEM   number of finite elements (hexagons or triangles) excluding
*         the virtual elements.
* NBMIX   number of macro-mixtures.
* SIDE    side of the hexagons.
* MAT     mixture index per hexagon.
* KN      element-ordered unknown list.
* QFR     element-ordered information.
* VOL     volume of hexagons.
* MERG    index of merged hexagons per macro-mixture.
* COUR    four times the incoming current per unit surface.
* FUNKNO  previously calculated fluxes.
* T       unit matrix.
* RH      unit matrix.
* RT      unit matrix.
*
*Parameters: output
* SOURCE  fission and diffusion sources.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IEX,MAXKN,MAXQF,NGCOND,NMERGE,NALBP,NREG,NUN,ISPLH,NELEM,
     1 NBMIX,MAT(NREG),KN(MAXKN),MERG(NBMIX)
      REAL FLXMER(NMERGE,NGCOND),SPH(NMERGE+NALBP,NGCOND),
     1 SUNMER(NMERGE,NGCOND),SIDE,QFR(MAXQF),VOL(NREG),COUR,
     2 FUNKNO(NUN,NGCOND),T(2),RH(6,6),RT(3,3),SOURCE(NUN)
*----
*  LOCAL VARIABLES
*----
      INTEGER ISR(6,2),ISRH(6,2),ISRT(3,2)
      REAL TH(6),RH2(6,6)
      DATA ISRH/2,1,4,5,6,3,1,4,5,6,3,2/
      DATA ISRT/1,2,3,2,3,1/
*----
*  RECOVER THE HEXAGONAL MASS (RH2) MATRICES.
*----
      IF(ISPLH.EQ.1) THEN
*        HEXAGONAL BASIS.
         LH=6
         DO 15 I=1,6
         DO 10 J=1,2
         ISR(I,J)=ISRH(I,J)
   10    CONTINUE
   15    CONTINUE
         DO 25 I=1,6
         TH(I)=0.0
         DO 20 J=1,6
         RH2(I,J)=RH(I,J)
         TH(I)=TH(I)+RH(I,J)
   20    CONTINUE
   25    CONTINUE
         CONST=1.5*SQRT(3.0)
         AA=SIDE
      ELSE
*        TRIANGULAR BASIS.
         LH=3
         DO 35 I=1,3
         DO 30 J=1,2
         ISR(I,J)=ISRT(I,J)
   30    CONTINUE
   35    CONTINUE
         DO 45 I=1,3
         TH(I)=0.0
         DO 40 J=1,3
         RH2(I,J)=RT(I,J)
         TH(I)=TH(I)+RT(I,J)
   40    CONTINUE
   45    CONTINUE
         CONST=0.25*SQRT(3.0)
         AA=SIDE/REAL(ISPLH-1)
      ENDIF
*
      CALL XDRSET(SOURCE,NUN,0.0)
*----
*  INCOMING CURRENT SOURCE.
*----
      IF(COUR.NE.0.0) THEN
         NUM1=0
         DO 95 K=1,NELEM
         KHEX=KN(NUM1+LH+1)
         IF(VOL(KHEX).EQ.0.0) GO TO 90
         DO 80 IC=1,LH
         QFR1=QFR(NUM1+IC)
         IF(QFR1.EQ.0.0) GO TO 80
         DO 70 I1=1,2
         IND1=KN(NUM1+ISR(IC,I1))
         IF(IND1.NE.0) SOURCE(IND1)=SOURCE(IND1)+COUR*QFR1*T(I1)
   70    CONTINUE
   80    CONTINUE
   90    NUM1=NUM1+LH+1
   95    CONTINUE
      ENDIF
*----
*  DIFFUSION AND FISSION SOURCE.
*----
      IF(IEX.EQ.1) GO TO 160
      DO 156 JGR=1,NGCOND
      NUM1=0
      DO 155 K=1,NELEM
      KHEX=KN(NUM1+LH+1)
      IF(VOL(KHEX).EQ.0.0) GO TO 150
      L=MAT(KHEX)
      VOL0=QFR(NUM1+LH+1)
      GARS=SUNMER(MERG(L),JGR)*SPH(MERG(L),JGR)
      DO 140 I=1,LH
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 140
      DO 130 J=1,LH
      IND2=KN(NUM1+J)
      IF(IND2.EQ.0) GO TO 130
      SOURCE(IND1)=SOURCE(IND1)+RH2(I,J)*FUNKNO(IND2,JGR)*VOL0*GARS
  130 CONTINUE
  140 CONTINUE
  150 NUM1=NUM1+LH+1
  155 CONTINUE
  156 CONTINUE
      RETURN
*----
*  INITIAL ITERATION.
*----
  160 NUM1=0
      DO 195 K=1,NELEM
      KHEX=KN(NUM1+LH+1)
      IF(VOL(KHEX).EQ.0.0) GO TO 190
      L=MAT(KHEX)
      VOL0=QFR(NUM1+LH+1)
      PV=0.0
      DO 170 JGR=1,NGCOND
      PV=PV+SUNMER(MERG(L),JGR)*FLXMER(MERG(L),JGR)
  170 CONTINUE
      DO 180 I=1,LH
      IND1=KN(NUM1+I)
      IF(IND1.NE.0) SOURCE(IND1)=SOURCE(IND1)+TH(I)*VOL0*PV
  180 CONTINUE
  190 NUM1=NUM1+LH+1
  195 CONTINUE
      RETURN
      END
