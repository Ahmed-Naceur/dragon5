*DECK BIVS03
      SUBROUTINE BIVS03(MAXKN,MAXQF,NREG,NUN,LL4,ISPLH,NELEM,IIMAX,SIDE,
     1 KN,QFR,BFR,VOL,IDL,MU,SOURCE,RH,RT,SYS,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* One-speed flux calculation in mesh corner finite difference or finite
* element approximation (hexagonal geometry).
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
* MAXQF   dimension of arrays QFR and BFR.
* NREG    number of hexagons in BIVAC.
* NUN     dimension of vector FUNKNO.
* LL4     order of the matrix SYS.
* ISPLH   hexagonal geometry flag:
*         =1: hexagonal elements; >1: triangular elements.
* NELEM   number of finite elements (hexagons or triangles) excluding
*         the virtual elements.
* IIMAX   allocated dimension of array SYS.
* SIDE    side of the hexagons.
* KN      element-ordered unknown list (dimensionned to KN((LH+1)*NELEM)
*         where LH=6 (hexagons) or 3 (triangles)).
* QFR     element-ordered information.
* BFR     element-ordered surface fractions.
* VOL     volume of the hexagons.
* IDL     position of the average flux component associated with each
*         hexagon.
* MU      indices used with the compressed diagonal storage mode matrix
*         SYS.
* SOURCE  fission and diffusion sources.
* RH      unit matrix.
* RT      unit matrix.
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
      INTEGER MAXKN,MAXQF,NREG,NUN,LL4,ISPLH,NELEM,IIMAX,KN(MAXKN),
     1 IDL(NREG),MU(LL4)
      REAL SIDE,QFR(MAXQF),BFR(MAXQF),VOL(NREG),SOURCE(LL4),RH(6,6),
     1 RT(3,3),SYS(IIMAX),FUNKNO(NUN)
*----
*  LOCAL VARIABLES
*----
      INTEGER ISR(6,2),ISRH(6,2),ISRT(3,2)
      REAL TH(6)
      DATA ISRH/2,1,4,5,6,3,1,4,5,6,3,2/
      DATA ISRT/1,2,3,2,3,1/
*----
*  RECOVER THE HEXAGONAL SVALAR PRODUCT (TH) VECTOR.
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
         TH(I)=TH(I)+RH(I,J)
   20    CONTINUE
   25    CONTINUE
         CONST=1.5*SQRT(3.0)
         CONSB=2.0*SQRT(3.0)/3.0
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
         TH(I)=TH(I)+RT(I,J)
   40    CONTINUE
   45    CONTINUE
         CONST=0.25*SQRT(3.0)
         CONSB=2.0*SQRT(3.0)
         AA=SIDE/REAL(ISPLH-1)
      ENDIF
*----
*  RESOLUTION.
*----
      DO 120 I=1,LL4
      FUNKNO(I)=SOURCE(I)
  120 CONTINUE
      CALL ALLDLS(LL4,MU,SYS,FUNKNO)
*----
*  VOLUME-AVERAGED FLUXES.
*----
      DO 130 KHEX=1,NREG
      IF(IDL(KHEX).NE.0) FUNKNO(IDL(KHEX))=0.0
  130 CONTINUE
      FUNKNO(NUN)=0.0
      NUM1=0
      DO 180 K=1,NELEM
      KHEX=KN(NUM1+LH+1)
      IF(VOL(KHEX).EQ.0.0) GO TO 170
      DO 140 I=1,LH
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 140
      SS=TH(I)*QFR(NUM1+LH+1)/(CONST*VOL(KHEX))
      FUNKNO(IDL(KHEX))=FUNKNO(IDL(KHEX))+SS*FUNKNO(IND1)
  140 CONTINUE
*----
*  SURFACE-AVERAGED FLUX.
*----
      DO 160 IC=1,LH
      BFR1=BFR(NUM1+IC)*CONSB
      IF(BFR1.EQ.0.0) GO TO 160
      DO 150 I1=1,2
      IND1=KN(NUM1+ISR(IC,I1))
      IF(IND1.GT.0) FUNKNO(NUN)=FUNKNO(NUN)+TH(I1)*FUNKNO(IND1)*BFR1
  150 CONTINUE
  160 CONTINUE
*
  170 NUM1=NUM1+LH+1
  180 CONTINUE
      RETURN
      END
