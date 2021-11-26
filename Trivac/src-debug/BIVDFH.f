*DECK BIVDFH
      SUBROUTINE BIVDFH (MAXEV,MAXKN,IMPX,ISPLH,LX,SIDE,NELEM,NUN,IHEX,
     1 NCODE,ICODE,ZCODE,MAT,VOL,IDL,KN,QFR,IQFR,BFR,MUW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a mesh centered finite difference
* discretization of a 2-D hexagonal geometry.
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
* MAXEV   maximum number of unknowns.
* MAXKN   dimension of arrays KN, QFR and BFR.
* IMPX    print parameter.
* ISPLH   hexagonal mesh-splitting flag:
*         =1 for complete hexagons; >1 for triangular mesh-splitting
*         into 6*(ISPLH-1)**2 triangles.
* LX      number of hexagons.
* SIDE    side of an hexagon.
* NCODE   type of boundary condition applied on each side
*         (i=1: X-  i=2: X+  i=3: Y-  i=4: Y+):
*         NCODE(I)=1: VOID;   NCODE(I)=2: REFL;   NCODE(I)=5: SYME;
*         NCODE(I)=7: ZERO.
* ICODE   physical albedo index on each side of the domain.
* ZCODE   albedo corresponding to boundary condition 'VOID' on each
*         side (ZCODE(I)=0.0 by default).
* MAT     mixture index assigned to each hexagon.
* IHEX    type of hexagonal boundary condition.
*
*Parameters: output
* NELEM   order of the system matrices (number of elements).
* NUN     number of unknowns per energy group.
* VOL     volume of each hexagon.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* IQFR    element-ordered physical albedo indices.
* BFR     element-ordered surface fractions.
* MUW     compressed storage mode indices.
* IDL     position of the average flux component associated with each
*         hexagon.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXEV,MAXKN,IMPX,ISPLH,LX,NELEM,NUN,IHEX,NCODE(4),
     1 ICODE(4),MAT(LX),IDL(LX),KN(MAXKN),IQFR(MAXKN),MUW(NELEM)
      REAL SIDE,ZCODE(4),VOL(LX),QFR(MAXKN),BFR(MAXKN)
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*
      IF(IMPX.GT.0) WRITE(6,500)
      CALL BIVSBH (MAXEV,MAXKN,IMPX,ISPLH,LX,SIDE,NELEM,IHEX,NCODE,MAT,
     1 VOL,KN,QFR)
*----
*  PRODUCE STANDARD MESH CENTERED FINITE DIFFERENCE NUMBERING.
*----
      IF(ISPLH.EQ.1) THEN
         NSURF=6
      ELSE
         NSURF=3
      ENDIF
      SURFTOT=0.0
      NUM1=0
      DO 200 KX=1,NELEM
      DO 190 IC=1,NSURF
      N1=ABS(KN(NUM1+IC))
      IF(N1.GT.NELEM) THEN
         IF(NCODE(1).EQ.1) THEN
            N1=-1
         ELSE IF(NCODE(1).EQ.2) THEN
            N1=-2
         ELSE IF(NCODE(1).EQ.7) THEN
            N1=-3
         ENDIF
      ELSE IF(N1.EQ.KX) THEN
         N1=-2
      ENDIF
      KN(NUM1+IC)=N1
*----
*  PROCESS BOUNDARY CONDITIONS.
*----
      IF(NSURF.EQ.6) THEN
         BFR(NUM1+IC)=QFR(NUM1+IC)*QFR(NUM1+7)/(1.5*SQRT(3.0)*SIDE)
         IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0))THEN
            QFR(NUM1+IC)=ALB(ZCODE(1))*QFR(NUM1+IC)
         ELSE IF(NCODE(1).NE.1) THEN
            QFR(NUM1+IC)=0.0
         ENDIF
      ELSE
         AA=SIDE/REAL(ISPLH-1)
         BFR(NUM1+IC)=QFR(NUM1+IC)*QFR(NUM1+4)/(0.25*SQRT(3.0)*AA)
         IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0))THEN
            QFR(NUM1+IC)=ALB(ZCODE(1))*QFR(NUM1+IC)
         ELSE IF(NCODE(1).NE.1) THEN
            QFR(NUM1+IC)=0.0
         ENDIF
      ENDIF
      IQFR(NUM1+IC)=ICODE(1)
      SURFTOT=SURFTOT+BFR(NUM1+IC)
  190 CONTINUE
      NUM1=NUM1+NSURF+1
  200 CONTINUE
*----
*  COMPUTE THE SURFACE FRACTIONS.
*----
      IF(SURFTOT.GT.0.0) THEN
         DO 210 I=1,NUM1
         BFR(I)=BFR(I)/SURFTOT
  210    CONTINUE
      ENDIF
*
      IF((IMPX.GT.1).AND.(NSURF.EQ.6)) THEN
         WRITE(6,510)
         NUM1=0
         DO 220 I=1,NELEM
         WRITE(6,520) I,KN(NUM1+7),(KN(NUM1+J),J=1,6),(QFR(NUM1+J),
     1   J=1,7)
         NUM1=NUM1+7
  220    CONTINUE
         NUM1=0
         WRITE (6,580)
         DO 225 I=1,NELEM
         IF(MAT(I).LE.0) GO TO 225
         WRITE (6,590) I,(BFR(NUM1+J),J=1,6)
         NUM1=NUM1+7
  225    CONTINUE
      ELSE IF((IMPX.GT.1).AND.(NSURF.EQ.3)) THEN
         WRITE(6,530)
         NUM1=0
         DO 230 I=1,NELEM
         WRITE(6,540) I,KN(NUM1+4),(KN(NUM1+J),J=1,3),(QFR(NUM1+J),
     1   J=1,4),(BFR(NUM1+J),J=1,3)
         NUM1=NUM1+4
  230    CONTINUE
      ENDIF
      IF(IMPX.GT.0) WRITE(6,570) NELEM
*----
*  COMPUTE THE SYSTEM MATRIX BANDWIDTH.
*----
      DO 240 I=1,NELEM
      MUW(I)=1
  240 CONTINUE
      NUM1=0
      DO 260 INW1=1,NELEM
      DO 250 I=1,NSURF
      IF(KN(NUM1+I).GT.0) THEN
         INW2=KN(NUM1+I)
         IF(INW2.LT.INW1) THEN
            MUW(INW1)=MAX(MUW(INW1),INW1-INW2+1)
         ENDIF
      ENDIF
  250 CONTINUE
      NUM1=NUM1+NSURF+1
  260 CONTINUE
      IIMAX=0
      DO 270 I=1,NELEM
      IIMAX=IIMAX+MUW(I)
      MUW(I)=IIMAX
  270 CONTINUE
      IF(IMPX.GT.6) WRITE(6,550) 'MUW :',(MUW(I),I=1,NELEM)
      IF(IMPX.GT.2) WRITE(6,560) IIMAX
*----
*  APPEND THE AVERAGED FLUXES AT THE END OF UNKNOWN VECTOR.
*----
      NUN=0
      IF(ISPLH.GT.1) NUN=NELEM
      DO 280 I=1,LX
      IF(MAT(I).EQ.0) THEN
         IDL(I)=0
      ELSE
         NUN=NUN+1
         IDL(I)=NUN
      ENDIF
  280 CONTINUE
      RETURN
*
  500 FORMAT(//52H BIVDFH: NUMBERING FOR A MESH CENTERED FINITE DIFFER,
     1 42HENCE DISCRETIZATION IN HEXAGONAL GEOMETRY.)
  510 FORMAT(/31H BIVDFH: NUMBERING OF UNKNOWNS./1X,30(1H-)/9X,
     1 7HHEXAGON,3X,9HNEIGHBOUR,28X,23HVOID BOUNDARY CONDITION,15X,
     2 6HVOLUME)
  520 FORMAT (1X,2I6,2X,6I6,2X,6F6.2,5X,1P,E13.6)
  530 FORMAT(/31H BIVDFH: NUMBERING OF UNKNOWNS./1X,30(1H-)/9X,
     1 7HHEXAGON,3X,8HUNKNOWNS,11X,23HVOID BOUNDARY CONDITION,12X,
     2 6HVOLUME,13X,16HSURFACE FRACTION)
  540 FORMAT (1X,2I6,2X,3I6,2X,1P,3E11.2,5X,E13.6,5X,3E10.2)
  550 FORMAT(/1X,A5/(1X,20I6))
  560 FORMAT(/52H NUMBER OF TERMS IN THE COMPRESSED SYSTEM MATRICES =,
     > I6)
  570 FORMAT(/39H BIVDFH: NUMBER OF UNKNOWNS PER GROUP =,I6/)
  580 FORMAT (//17H SURFACE FRACTION//8H HEXAGON,5X,3HBFR)
  590 FORMAT (3X,I4,4X,1P,6E10.2)
      END
