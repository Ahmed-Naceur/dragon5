*DECK BIVFSO
      SUBROUTINE BIVFSO (IEX,MAXKN,NGCOND,NMERGE,NALBP,FLXMER,SPH,
     1 SUNMER,CYLIND,NREG,NUN,NBMIX,XX,DD,MAT,KN,QFR,VOL,MERG,COUR,
     2 FUNKNO,LC,T,TS,R,RS,SOURCE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Source term calculation for finite element or mesh corner finite
* differences in Cartesian geometry. 
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
* NGCOND  number of groups condensed.
* NMERGE  number of mixtures.
* NALBP   number of physical albedos.
* FLXMER  flux estimate per mixture.
* SPH     SPH factors.
* SUNMER  incoming source (scattering+fission) cross sections.
* CYLIND  cylinderization flag (=.TRUE. for cylindrical geometry).
* NREG    number of regions in BIVAC.
* NUN     number of unknown per group in BIVAC.
* NBMIX   number of macro-mixtures.
* XX      X-directed mesh spacings.
* DD      value used with a cylindrical geometry.
* MAT     mixture index per region.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* VOL     volume of regions.
* MERG    index of merged regions per macro-mixture.
* COUR    four times the incoming current per unit surface.
* FUNKNO  previously calculated fluxes.
* LC      order of the unit matrices.
* T       unit matrix.
* TS      unit matrix.
* R       unit matrix.
* RS      unit matrix.
*
*Parameters: output
* SOURCE  fission and diffusion sources.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IEX,MAXKN,NGCOND,NMERGE,NALBP,NREG,NUN,NBMIX,MAT(NREG),
     1 KN(MAXKN),MERG(NBMIX),LC
      REAL FLXMER(NMERGE,NGCOND),SPH(NMERGE+NALBP,NGCOND),
     1 SUNMER(NMERGE,NGCOND),XX(NREG),DD(NREG),QFR(4*NREG),VOL(NREG),
     2 COUR,FUNKNO(NUN,NGCOND),T(LC),TS(LC),R(LC,LC),RS(LC,LC),
     3 SOURCE(NUN)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      PARAMETER (DPI=6.283185307)
      INTEGER IJ1(25),IJ2(25),ISR(4,5)
*----
*  COMPUTE VECTORS IJ1, IJ2 AND MATRIX ISR.
*----
      CRZ=0.0
      LL=LC*LC
      DO 10 I=1,LL
      IJ1(I)=1+MOD(I-1,LC)
      IJ2(I)=1+(I-IJ1(I))/LC
   10 CONTINUE
      DO 20 I=1,LC
      ISR(1,I)=(I-1)*LC+1
      ISR(2,I)=I*LC
      ISR(3,I)=I
      ISR(4,I)=LL-LC+I
   20 CONTINUE
*
      CALL XDRSET(SOURCE,NUN,0.0)
*----
*  INCOMING CURRENT SOURCE.
*----
      IF(COUR.NE.0.0) THEN
         NUM1=0
         NUM2=0
         DO 80 K=1,NREG
         IF(MAT(K).EQ.0) GO TO 80
         IF(VOL(K).EQ.0.0) GO TO 70
         DO 60 IC=1,4
         QFR1=QFR(NUM2+IC)
         IF(QFR1.EQ.0.0) GO TO 60
         DO 50 I1=1,LC
         IND1=KN(NUM1+ISR(IC,I1))
         IF(IND1.EQ.0) GO TO 50
         IF(CYLIND) THEN
            IF(IC.EQ.1) THEN
               CRZ=-0.5*T(I1)
            ELSE IF(IC.EQ.2) THEN
               CRZ=0.5*T(I1)
            ELSE IF(IC.EQ.3) THEN
               CRZ=TS(I1)
            ELSE IF(IC.EQ.4) THEN
               CRZ=TS(I1)
            ENDIF
            RR=DPI*(XX(K)*CRZ+DD(K)*T(I1))
         ELSE
            RR=T(I1)
         ENDIF
         SOURCE(IND1)=SOURCE(IND1)+COUR*QFR1*RR
   50    CONTINUE
   60    CONTINUE
   70    NUM1=NUM1+LL
         NUM2=NUM2+4
   80    CONTINUE
      ENDIF
*----
*  DIFFUSION AND FISSION SOURCE.
*----
      IF(IEX.EQ.1) GO TO 165
      DO 161 JGR=1,NGCOND
      NUM1=0
      DO 160 K=1,NREG
      L=MAT(K)
      IF(L.EQ.0) GO TO 160
      IF(VOL(K).EQ.0.0) GO TO 150
      GARS=SUNMER(MERG(L),JGR)*SPH(MERG(L),JGR)
      DO 140 I=1,LL
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 140
      I1=IJ1(I)
      I2=IJ2(I)
      DO 130 J=1,LL
      IND2=KN(NUM1+J)
      IF(IND2.EQ.0) GO TO 130
      IF(CYLIND) THEN
         J1=IJ1(J)
         RR=DPI*(R(I1,J1)*DD(K)+RS(I1,J1)*XX(K))*R(I2,IJ2(J))
      ELSE
         RR=R(I1,IJ1(J))*R(I2,IJ2(J))
      ENDIF
      SOURCE(IND1)=SOURCE(IND1)+RR*FUNKNO(IND2,JGR)*VOL(K)*GARS
  130 CONTINUE
  140 CONTINUE
  150 NUM1=NUM1+LL
  160 CONTINUE
  161 CONTINUE
      RETURN
*----
*  INITIAL ITERATION.
*----
  165 NUM1=0
      DO 200 K=1,NREG
      L=MAT(K)
      IF(L.EQ.0) GO TO 200
      IF(VOL(K).EQ.0.0) GO TO 190
      PV=0.0
      DO 170 JGR=1,NGCOND
      PV=PV+SUNMER(MERG(L),JGR)*FLXMER(MERG(L),JGR)
  170 CONTINUE
      DO 180 I=1,LL
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 180
      SS=T(IJ1(I))*T(IJ2(I))*VOL(K)
      IF(CYLIND) THEN
         SS=DPI*(T(IJ1(I))*DD(K)+TS(IJ1(I))*XX(K))*T(IJ2(I))*VOL(K)
      ENDIF
      SOURCE(IND1)=SOURCE(IND1)+SS*PV
  180 CONTINUE
  190 NUM1=NUM1+LL
  200 CONTINUE
      RETURN
      END
