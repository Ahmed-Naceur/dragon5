*DECK BIVS01
      SUBROUTINE BIVS01(MAXKN,CYLIND,NREG,NUN,LL4,IIMAX,XX,DD,MAT,KN,
     1 BFR,VOL,IDL,MU,SOURCE,LC,T,TS,SYS,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* One-speed flux calculation in mesh corner finite difference or finite
* element approximation (Cartesian geometry). 
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
* MAXKN   dimension of array KN.
* CYLIND  cylinderization flag (=.TRUE. for cylindrical geometry).
* NREG    number of elements in BIVAC.
* NUN     dimension of vector FUNKNO.
* LL4     order of matrix SYS.
* IIMAX   allocated dimension of array SYS.
* XX      X-directed mesh spacings.
* DD      value used with a cylindrical geometry.
* MAT     mixture index per region.
* KN      element-ordered unknown list.
* BFR     element-ordered surface fractions.
* VOL     volume of regions.
* IDL     position of integrated fluxes into unknown vector.
* MU      indices used with compressed diagonal storage mode matrix SYS.
* SOURCE  fission and diffusion sources.
* LC      number of polynomials in a complete 1-D basis.
* T       Cartesian linear product vector.
* TS      cylindrical linear product vector.
* SYS     factorized system matrix.
*
*Parameters: output
* FUNKNO  unknown array. The first LL4 values contains the finite
*         element unknowns; the next NREG values contains element
*         averaged fluxes. The surface-averaged flux is located in
*         position FUNKNO(NUN).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXKN,NREG,NUN,LL4,IIMAX,MAT(NREG),KN(MAXKN),IDL(NREG),
     1 MU(LL4),LC
      REAL XX(NREG),DD(NREG),BFR(4*NREG),SOURCE(LL4),VOL(NREG),T(LC),
     1 TS(LC),SYS(IIMAX),FUNKNO(NUN)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      INTEGER IJ1(25),IJ2(25),ISR(4,5)
*----
*  COMPUTE VECTORS IJ1, IJ2 AND MATRIX ISR.
*----
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
*----
*  RESOLUTION.
*----
      DO 30 I=1,LL4
      FUNKNO(I)=SOURCE(I)
   30 CONTINUE
      CALL ALLDLS(LL4,MU,SYS,FUNKNO)
*----
*  CALCULATION OF ELEMENT-AVERAGED AND SURFACE-AVERAGED FLUXES.
*----
      FUNKNO(NUN)=0.0
      NUM1=0
      NUM2=0
      DO 170 K=1,NREG
      IF(MAT(K).EQ.0) GO TO 170
      FUNKNO(IDL(K))=0.0
      IF(VOL(K).EQ.0.0) GO TO 160
*
      DO 130 I=1,LL
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 130
      IF(CYLIND) THEN
         SS=(T(IJ1(I))+TS(IJ1(I))*XX(K)/DD(K))*T(IJ2(I))
      ELSE
         SS=T(IJ1(I))*T(IJ2(I))
      ENDIF
      FUNKNO(IDL(K))=FUNKNO(IDL(K))+SS*FUNKNO(IND1)
  130 CONTINUE
*
      DO 150 IC=1,4
      BFR1=BFR(NUM2+IC)
      IF(BFR1.EQ.0.0) GO TO 150
      DO 140 I1=1,LC
      IND1=KN(NUM1+ISR(IC,I1))
      IF(IND1.GT.0) FUNKNO(NUN)=FUNKNO(NUN)+T(I1)*FUNKNO(IND1)*BFR1
  140 CONTINUE
  150 CONTINUE
*
  160 NUM1=NUM1+LL
      NUM2=NUM2+4
  170 CONTINUE
      RETURN
      END
