*DECK TRISFH
      SUBROUTINE TRISFH (IMPX,MAXKN,MAXIP,NBLOS,ISPLH,IELEM,LXH,LZ,MAT,
     1 SIDE,ZZZ,NCODE,ICODE,ZCODE,LL4,LL4F,LL4W,LL4X,LL4Y,LL4Z,VOL,
     2 IDL,IPERT,ZZ,FRZ,KN,QFR,IQFR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a Thomas-Raviart-Schneider finite element
* discretization of a 3-D hexagonal geometry.
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
* IMPX    print parameter.
* MAXKN   number of components in KN.
* MAXIP   maximum number of currents
* NBLOS   number of lozenges per direction in 3D with mesh-splitting.
* ISPLH   mesh-splitting in 3*ISPLH**2 lozenges per hexagon.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic).
* LXH     number of hexagons in a plane.
* LZ      number of axial planes.
* MAT     mixture index assigned to each lozenge.
* SIDE    side of a lozenge.
* ZZZ     Z-coordinates of the axial planes.
* NCODE   type of boundary condition applied on each side (I=1: hbc):
*         NCODE(I)=1: VOID;          =2: REFL;       =6: ALBE;
*                 =5: SYME;          =7: ZERO.
* ICODE   physical albedo index on each side of the domain.
* ZCODE   albedo corresponding to boundary condition 'VOID' on each
*         side (ZCODE(I)=0.0 by default).
*
*Parameters: output
* LL4     order of the system matrices.
* LL4F    number of flux unknowns.
* LL4W    number of W-directed currents
* LL4X    number of X-directed currents
* LL4Y    number of Y-directed currents
* LL4Z    number of Z-directed currents
* ZZ      Z-sides of each hexagon.
* FRZ     volume fractions for the axial SYME boundary condition.
* VOL     volume of each lozenge.
* IDL     position of the average flux component associated with each
*         lozenge.
* IPERT   mixture permutation index.
* KN      ADI permutation indices for the volumes and currents.
* QFR     element-ordered boundary conditions.
* IQFR    element-ordered physical albedo indices.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,MAXKN,MAXIP,NBLOS,ISPLH,IELEM,LXH,LZ,
     1 MAT(3,ISPLH**2,LXH*LZ),NCODE(6),ICODE(6),LL4,LL4F,LL4W,LL4X,
     2 LL4Y,LL4Z,IDL(3,NBLOS),IPERT(NBLOS),KN(NBLOS,MAXKN/NBLOS),
     3 IQFR(NBLOS,8)
      REAL SIDE,ZZZ(LZ+1),ZCODE(6),VOL(3,NBLOS),ZZ(3,NBLOS),
     1 FRZ(NBLOS),QFR(NBLOS,8)
*----
*  LOCAL VARIABLES
*----
      LOGICAL COND,LL1,LL2
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IJP
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IZGLOB
      INTEGER, DIMENSION(:), ALLOCATABLE :: IP,I1,I3,I4,I5
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJP(LXH,ISPLH,ISPLH),IP(MAXIP),IZGLOB(NBLOS,3))
*----
*  THOMAS-RAVIART-SCHNEIDER SPECIFIC NUMEROTATION
*----
      NBC=INT((SQRT(REAL((4*LXH-1)/3))+1.)/2.)
      IF(LXH.NE.1+3*NBC*(NBC-1)) CALL XABORT('TRISFH: INVALID VALUE OF '
     1 //'LXH(1).')
      IF(ISPLH.EQ.1) THEN
         DO 10 I=1,LXH
         IJP(I,1,1)=I
   10    CONTINUE
      ELSE
         I=0
         DO 23 I0=1,2*NBC-1
         JMAX=NBC+I0-1
         IF(I0.GE.NBC) JMAX=3*NBC-I0-1
         IKEEP=I
         DO 22 J0=1,JMAX
         I=I+1
         DO 21 IM=1,ISPLH
         DO 20 JM=1,ISPLH
         IJP(I,IM,JM)=ISPLH*(IKEEP*ISPLH+(IM-1)*JMAX+J0-1)+JM
   20    CONTINUE
   21    CONTINUE
   22    CONTINUE
   23    CONTINUE
         IF(I.NE.LXH) CALL XABORT('TRISFH: INVALID VALUE OF LXH(2)')
      ENDIF
      ALLOCATE(I1(3*LXH),I3(2*LXH),I4(NBLOS),I5(NBLOS))
      DO 25 I=1,LXH
      I3(I)=I
   25 CONTINUE
      DO 30 I=1,LXH*LZ
      I4(I)=0
      IF(MAT(1,1,I).GT.0) I4(I)=I
   30 CONTINUE
      CALL XDISET(IZGLOB,3*NBLOS,0)
      J1=2+3*(NBC-1)*(NBC-2)
      IF(NBC.EQ.1) J1=1
      J3=J1+2*NBC-2
      J5=J3+2*NBC-2
      CALL BIVPER(J1,1,LXH,LXH,I1(1),I3)
      CALL BIVPER(J3,3,LXH,LXH,I1(LXH+1),I3)
      CALL BIVPER(J5,5,LXH,LXH,I1(2*LXH+1),I3)
      I=0
      DO 43 IZ=1,LZ
      DO 42 IX=1,LXH
      I=I+1
      IOFW=I1(IX)
      IOFX=I1(LXH+IX)
      IOFY=I1(2*LXH+IX)
      DO 41 IM=1,ISPLH
      DO 40 JM=1,ISPLH
      IZGLOB((IZ-1)*LXH*ISPLH**2+IJP(IOFW,IM,JM),1)=I4(I)
      IZGLOB((IZ-1)*LXH*ISPLH**2+IJP(IOFX,IM,JM),2)=I4(I)
      IZGLOB((IZ-1)*LXH*ISPLH**2+IJP(IOFY,IM,JM),3)=I4(I)
   40 CONTINUE
   41 CONTINUE
   42 CONTINUE
   43 CONTINUE
      DO 50 I=1,LXH
      II1=I1(I)
      II2=I1(LXH+I)
      II3=I1(2*LXH+I)
      I3(II1)=II2
      I3(LXH+II1)=II3
   50 CONTINUE
*----
*  COMPUTE THE FLUX PERMUTATION PART OF MATRIX KN (W <--> X)
*----
      CALL XDISET(KN,NBLOS*(3+6*IELEM*IELEM*(IELEM+2)),0)
      LT4=0
      DO 70 II2=1,NBLOS
      I=IZGLOB(II2,1)
      I4(II2)=0
      IF(I.NE.0) THEN
         LT4=LT4+1
         I4(II2)=LT4
      ENDIF
   70 CONTINUE
      LT4=0
      DO 80 II2=1,NBLOS
      I=IZGLOB(II2,2)
      I5(II2)=0
      IF(I.NE.0) THEN
         LT4=LT4+1
         I5(II2)=LT4
      ENDIF
   80 CONTINUE
      IF(ISPLH.EQ.1) THEN
         I=0
         DO 95 IZ=1,LZ
         DO 90 IX=1,LXH
         I=I+1
         IF(IZGLOB(I,1).EQ.0) GO TO 90
         IOF=(IZ-1)*LXH+I3(IX)
         KN(I4(I),1)=I5(IOF)+LT4
   90    CONTINUE
   95    CONTINUE
      ELSE
         I=0
         DO 105 I0=1,2*NBC-1
         JMAX=NBC+I0-1
         IF(I0.GE.NBC) JMAX=3*NBC-I0-1
         IKEEP=I
         DO 100 J0=1,JMAX
         I=I+1
         I1(I)=JMAX
         I1(LXH+I)=IKEEP
         I1(2*LXH+I)=J0
  100    CONTINUE
  105    CONTINUE
         DO 125 IZ=1,LZ
         DO 120 I=1,LXH
         JMAX=I1(I)
         IKEEP=I1(LXH+I)
         J00=I1(2*LXH+I)
         KMAX=I1(I3(I))
         JKEEP=I1(LXH+I3(I))
         K0=I1(2*LXH+I3(I))
         DO 115 IM=1,ISPLH
         DO 110 JM=1,ISPLH
         II1=ISPLH*(IKEEP*ISPLH+(IM-1)*JMAX+J00-1)+JM
         IOF1=(IZ-1)*LXH*ISPLH**2+II1
         IF(IZGLOB(IOF1,1).EQ.0) GO TO 120
         II2=ISPLH*(JKEEP*ISPLH+(ISPLH-JM)*KMAX+K0-1)+IM
         IOF2=(IZ-1)*LXH*ISPLH**2+II2
         KN(I4(IOF1),1)=I5(IOF2)+LT4
  110    CONTINUE
  115    CONTINUE
  120    CONTINUE
  125    CONTINUE
      ENDIF
*----
*  COMPUTE THE FLUX PERMUTATION PART OF MATRIX KN (X <--> Y)
*----
      LT4=0
      DO 130 II2=1,NBLOS
      I=IZGLOB(II2,3)
      I5(II2)=0
      IF(I.NE.0) THEN
         LT4=LT4+1
         I5(II2)=LT4
      ENDIF
  130 CONTINUE
      IF(ISPLH.EQ.1) THEN
         I=0
         DO 145 IZ=1,LZ
         DO 140 IX=1,LXH
         I=I+1
         IF(IZGLOB(I,1).EQ.0) GO TO 140
         IOF=(IZ-1)*LXH+I3(LXH+IX)
         KN(I4(I),2)=I5(IOF)+2*LT4
  140    CONTINUE
  145    CONTINUE
      ELSE
         I=0
         DO 155 I0=1,2*NBC-1
         JMAX=NBC+I0-1
         IF(I0.GE.NBC) JMAX=3*NBC-I0-1
         IKEEP=I
         DO 150 J0=1,JMAX
         I=I+1
         I1(I)=JMAX
         I1(LXH+I)=IKEEP
         I1(2*LXH+I)=J0
  150    CONTINUE
  155    CONTINUE
         DO 175 IZ=1,LZ
         DO 170 I=1,LXH
         JMAX=I1(I)
         IKEEP=I1(LXH+I)
         J00=I1(2*LXH+I)
         KMAX=I1(I3(LXH+I))
         JKEEP=I1(LXH+I3(LXH+I))
         K0=I1(2*LXH+I3(LXH+I))
         DO 165 IM=1,ISPLH
         DO 160 JM=1,ISPLH
         II1=ISPLH*(IKEEP*ISPLH+(IM-1)*JMAX+J00-1)+JM
         IOF1=(IZ-1)*LXH*ISPLH**2+II1
         IF(IZGLOB(IOF1,1).EQ.0) GO TO 170
         II2=ISPLH*(JKEEP*ISPLH+(ISPLH-IM)*KMAX+K0-1)+(ISPLH-JM+1)
         IOF2=(IZ-1)*LXH*ISPLH**2+II2
         KN(I4(IOF1),2)=I5(IOF2)+2*LT4
  160    CONTINUE
  165    CONTINUE
  170    CONTINUE
  175    CONTINUE
      ENDIF
*----
*  COMPUTE THE FLUX PERMUTATION PART OF MATRIX KN (Y <--> W)
*----
      IF(ISPLH.EQ.1) THEN
         DO 180 I=1,LXH*LZ
         IF(IZGLOB(I,1).EQ.0) GO TO 180
         KN(I4(I),3)=I4(I)
  180    CONTINUE
      ELSE
         I=0
         DO 195 I0=1,2*NBC-1
         JMAX=NBC+I0-1
         IF(I0.GE.NBC) JMAX=3*NBC-I0-1
         IKEEP=I
         DO 190 J0=1,JMAX
         I=I+1
         I1(I)=JMAX
         I1(LXH+I)=IKEEP
         I1(2*LXH+I)=J0
  190    CONTINUE
  195    CONTINUE
         DO 215 IZ=1,LZ
         DO 210 I=1,LXH
         JMAX=I1(I)
         IKEEP=I1(LXH+I)
         J00=I1(2*LXH+I)
         DO 205 IM=1,ISPLH
         DO 200 JM=1,ISPLH
         II1=ISPLH*(IKEEP*ISPLH+(IM-1)*JMAX+J00-1)+JM
         IOF1=(IZ-1)*LXH*ISPLH**2+II1
         IF(IZGLOB(IOF1,1).EQ.0) GO TO 210
         II2=ISPLH*(IKEEP*ISPLH+(JM-1)*JMAX+J00-1)+(ISPLH-IM+1)
         IOF2=(IZ-1)*LXH*ISPLH**2+II2
         KN(I4(IOF1),3)=I4(IOF2)
  200    CONTINUE
  205    CONTINUE
  210    CONTINUE
  215    CONTINUE
      ENDIF
      DEALLOCATE(I5,I4,I3,I1)
*----
*  SET THE CURRENT NUMBERING PART OF MATRIX KN AND MATRIX QFR (W-AXIS)
*----
      LL4W0=(2*LXH*ISPLH*IELEM+2*NBC-1)*ISPLH*LZ*IELEM**2
      LL4Z0=3*LXH*(LZ+1)*(ISPLH**2)*IELEM**2
      LL4F=3*LT4*IELEM**3
      CALL XDRSET(QFR,8*NBLOS,0.0)
      CALL XDISET(IQFR,8*NBLOS,0)
      ALBEDO=0.5*(1.0-ZCODE(1))/(1.0+ZCODE(1))
      NELEH=(IELEM+1)*IELEM**2
      NELEZ=6*IELEM**2
      NB1=2*NBC*ISPLH*IELEM+1
      NB2=2*(2*NBC-1)*ISPLH*IELEM+1
      KEL=0
      NDDIR=0
      NUM=0
      DO 345 IZ=1,LZ
      FRACT=1.0
      IF((NCODE(5).EQ.5).AND.(IZ.EQ.1)) FRACT=0.5
      IF((NCODE(6).EQ.5).AND.(IZ.EQ.LZ)) FRACT=0.5
      DZZ=ZZZ(IZ+1)-ZZZ(IZ)
      DO 290 JSTAGE=1,NBC
      DO 282 JEL=1,ISPLH
      DO 281 IRANG=1,NBC+JSTAGE-1
      DO 280 IEL=1,ISPLH
      KEL=KEL+1
      IF(IZGLOB(KEL,1).EQ.0) GO TO 280
      NUM=NUM+1
      IF((IRANG.EQ.1).AND.(IEL.EQ.1)) THEN
         LL1=.TRUE.
      ELSE
         LL1=(IZGLOB(KEL-1,1).EQ.0)
      ENDIF
      IF((IRANG.EQ.NBC+JSTAGE-1).AND.(IEL.EQ.ISPLH)) THEN
         LL2=.TRUE.
      ELSE
         LL2=(IZGLOB(KEL+1,1).EQ.0)
      ENDIF
      LCOUR=0
      DO 255 J=1,IELEM**2
      DO 250 I=1,IELEM+1
      LCOUR=LCOUR+1
      ITEMP = NDDIR
     >      + (JEL-1)*(2*(NBC+JSTAGE-1)*IELEM*ISPLH+1)*IELEM**2
     >      + (IRANG-1)*(2*IELEM*ISPLH)
     >      + (IEL-1)*IELEM
     >      + (J-1)*(2*(NBC+JSTAGE-1)*IELEM*ISPLH+1) + I
      IF(LCOUR.GT.NELEH) CALL XABORT('TRISFH: BUG1')
      IF(KEL.GT.NBLOS) CALL XABORT('TRISFH: BUG2')
      KN(NUM,3+LCOUR)=ITEMP
      KN(NUM,3+NELEH+LCOUR)=ITEMP+IELEM*ISPLH
  250 CONTINUE
  255 CONTINUE
      IF(LL1) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 260 I=1,IELEM**2
            KN(NUM,3+(I-1)*(IELEM+1)+1)=0
  260       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(NUM,1)=SIDE*DZZ*FRACT/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(NUM,1)=SIDE*DZZ*FRACT
            IQFR(NUM,1)=ICODE(1)
         ENDIF
      ENDIF
      IF(LL2) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 270 I=1,IELEM**2
            KN(NUM,3+NELEH+I*(IELEM+1))=0
  270       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(NUM,2)=SIDE*DZZ*FRACT/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(NUM,2)=SIDE*DZZ*FRACT
            IQFR(NUM,1)=ICODE(1)
         ENDIF
      ENDIF
  280 CONTINUE
  281 CONTINUE
  282 CONTINUE
      NDDIR=NDDIR+(NB1+2*(JSTAGE-1)*ISPLH*IELEM)*ISPLH*IELEM**2
  290 CONTINUE
*
      DO 340 JSTAGE=NBC+1,2*NBC-1
      DO 332 JEL=1,ISPLH
      DO 331 IRANG=1,(2*NBC-2)-(JSTAGE-NBC-1)
      DO 330 IEL=1,ISPLH
      KEL=KEL+1
      IF(IZGLOB(KEL,1).EQ.0) GO TO 330
      NUM=NUM+1
      IF((IRANG.EQ.1).AND.(IEL.EQ.1)) THEN
         LL1=.TRUE.
      ELSE
         LL1=(IZGLOB(KEL-1,1).EQ.0)
      ENDIF
      IF((IRANG.EQ.(2*NBC-2)-(JSTAGE-NBC-1)).AND.(IEL.EQ.ISPLH)) THEN
         LL2=.TRUE.
      ELSE
         LL2=(IZGLOB(KEL+1,1).EQ.0)
      ENDIF
      LCOUR=0
      DO 305 J=1,IELEM**2
      DO 300 I=1,IELEM+1
      LCOUR=LCOUR+1
      ITEMP = NDDIR
     >      + (JEL-1)*(2*(2*NBC-1+NBC-JSTAGE)*IELEM*ISPLH+1)*IELEM**2
     >      + (IRANG-1)*(2*IELEM*ISPLH)
     >      + (IEL-1)*IELEM
     >      + (J-1)*(2*(2*NBC-1+NBC-JSTAGE)*IELEM*ISPLH+1) + I
      IF(LCOUR.GT.NELEH) CALL XABORT('TRISFH: BUG3')
      IF(KEL.GT.NBLOS) CALL XABORT('TRISFH: BUG4')
      KN(NUM,3+LCOUR)=ITEMP
      KN(NUM,3+NELEH+LCOUR)=ITEMP+IELEM*ISPLH
  300 CONTINUE
  305 CONTINUE
      IF(LL1) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 310 I=1,IELEM**2
            KN(NUM,3+(I-1)*(IELEM+1)+1)=0
  310       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(NUM,1)=SIDE*DZZ*FRACT/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(NUM,1)=SIDE*DZZ*FRACT
            IQFR(NUM,1)=ICODE(1)
         ENDIF
      ENDIF
      IF(LL2) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 320 I=1,IELEM**2
            KN(NUM,3+NELEH+I*(IELEM+1))=0
  320       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(NUM,2)=SIDE*DZZ*FRACT/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(NUM,2)=SIDE*DZZ*FRACT
            IQFR(NUM,2)=ICODE(1)
         ENDIF
      ENDIF
  330 CONTINUE
  331 CONTINUE
  332 CONTINUE
      NDDIR=NDDIR+(NB2-2*(JSTAGE-NBC)*ISPLH*IELEM)*ISPLH*IELEM**2
  340 CONTINUE
  345 CONTINUE
*----
*  SET THE CURRENT NUMBERING PART OF MATRIX KN AND MATRIX QFR (X-AXIS)
*----
      CALL XDISET(IP,NBLOS,0)
      DO 350 NUM=1,LT4
      IP(KN(NUM,1)-LT4)=NUM
  350 CONTINUE
      KEL=0
      NUM=0
      DO 455 IZ=1,LZ
      FRACT=1.0
      IF((NCODE(5).EQ.5).AND.(IZ.EQ.1)) FRACT=0.5
      IF((NCODE(6).EQ.5).AND.(IZ.EQ.LZ)) FRACT=0.5
      DZZ=ZZZ(IZ+1)-ZZZ(IZ)
      DO 400 JSTAGE=1,NBC
      DO 392 JEL=1,ISPLH
      DO 391 IRANG=1,NBC+JSTAGE-1
      DO 390 IEL=1,ISPLH
      KEL=KEL+1
      IF(IZGLOB(KEL,2).EQ.0) GO TO 390
      NUM=NUM+1
      IF((IRANG.EQ.1).AND.(IEL.EQ.1)) THEN
         LL1=.TRUE.
      ELSE
         LL1=(IZGLOB(KEL-1,2).EQ.0)
      ENDIF
      IF((IRANG.EQ.NBC+JSTAGE-1).AND.(IEL.EQ.ISPLH)) THEN
         LL2=.TRUE.
      ELSE
         LL2=(IZGLOB(KEL+1,2).EQ.0)
      ENDIF
      LCOUR=0
      DO 365 J=1,IELEM**2
      DO 360 I=1,IELEM+1
      LCOUR=LCOUR+1
      ITEMP = NDDIR
     >      + (JEL-1)*(2*(NBC+JSTAGE-1)*IELEM*ISPLH+1)*IELEM**2
     >      + (IRANG-1)*(2*IELEM*ISPLH)
     >      + (IEL-1)*IELEM
     >      + (J-1)*(2*(NBC+JSTAGE-1)*IELEM*ISPLH+1) + I
      IF(LCOUR.GT.NELEH) CALL XABORT('TRISFH: BUG5')
      IF(KEL.GT.NBLOS) CALL XABORT('TRISFH: BUG6')
      KN(IP(NUM),3+2*NELEH+LCOUR)=ITEMP
      KN(IP(NUM),3+3*NELEH+LCOUR)=ITEMP+IELEM*ISPLH
  360 CONTINUE
  365 CONTINUE
      IF(LL1) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 370 I=1,IELEM**2
            KN(IP(NUM),3+2*NELEH+(I-1)*(IELEM+1)+1)=0
  370       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),3)=SIDE*DZZ*FRACT/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),3)=SIDE*DZZ*FRACT
            IQFR(IP(NUM),3)=ICODE(1)
         ENDIF
      ENDIF
      IF(LL2) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 380 I=1,IELEM**2
            KN(IP(NUM),3+3*NELEH+I*(IELEM+1))=0
  380       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),4)=SIDE*DZZ*FRACT/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),4)=SIDE*DZZ*FRACT
            IQFR(IP(NUM),4)=ICODE(1)
         ENDIF
      ENDIF
  390 CONTINUE
  391 CONTINUE
  392 CONTINUE
      NDDIR=NDDIR+(NB1+2*(JSTAGE-1)*ISPLH*IELEM)*ISPLH*IELEM**2
  400 CONTINUE
*
      DO 450 JSTAGE=NBC+1,2*NBC-1
      DO 442 JEL=1,ISPLH
      DO 441 IRANG=1,(2*NBC-2)-(JSTAGE-NBC-1)
      DO 440 IEL=1,ISPLH
      KEL=KEL+1
      IF(IZGLOB(KEL,2).EQ.0) GO TO 440
      NUM=NUM+1
      IF((IRANG.EQ.1).AND.(IEL.EQ.1)) THEN
         LL1=.TRUE.
      ELSE
         LL1=(IZGLOB(KEL-1,2).EQ.0)
      ENDIF
      IF((IRANG.EQ.(2*NBC-2)-(JSTAGE-NBC-1)).AND.(IEL.EQ.ISPLH)) THEN
         LL2=.TRUE.
      ELSE
         LL2=(IZGLOB(KEL+1,2).EQ.0)
      ENDIF
      LCOUR=0
      DO 415 J=1,IELEM**2
      DO 410 I=1,IELEM+1
      LCOUR=LCOUR+1
      ITEMP = NDDIR
     >      + (JEL-1)*(2*(2*NBC-1+NBC-JSTAGE)*IELEM*ISPLH+1)*IELEM**2
     >      + (IRANG-1)*(2*IELEM*ISPLH)
     >      + (IEL-1)*IELEM
     >      + (J-1)*(2*(2*NBC-1+NBC-JSTAGE)*IELEM*ISPLH+1) + I
      IF(LCOUR.GT.NELEH) CALL XABORT('TRISFH: BUG7')
      IF(KEL.GT.NBLOS) CALL XABORT('TRISFH: BUG8')
      KN(IP(NUM),3+2*NELEH+LCOUR)=ITEMP
      KN(IP(NUM),3+3*NELEH+LCOUR)=ITEMP+IELEM*ISPLH
  410 CONTINUE
  415 CONTINUE
      IF(LL1) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 420 I=1,IELEM**2
            KN(IP(NUM),3+2*NELEH+(I-1)*(IELEM+1)+1)=0
  420       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),3)=SIDE*DZZ*FRACT/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),3)=SIDE*DZZ*FRACT
            IQFR(IP(NUM),3)=ICODE(1)
         ENDIF
      ENDIF
      IF(LL2) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 430 I=1,IELEM**2
            KN(IP(NUM),3+3*NELEH+I*(IELEM+1))=0
  430       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),4)=SIDE*DZZ*FRACT/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),4)=SIDE*DZZ*FRACT
            IQFR(IP(NUM),4)=ICODE(1)
         ENDIF
      ENDIF
  440 CONTINUE
  441 CONTINUE
  442 CONTINUE
      NDDIR=NDDIR+(NB2-2*(JSTAGE-NBC)*ISPLH*IELEM)*ISPLH*IELEM**2
  450 CONTINUE
  455 CONTINUE
*----
*  SET THE CURRENT NUMBERING PART OF MATRIX KN AND MATRIX QFR (Y-AXIS)
*----
      CALL XDISET(IP,NBLOS,0)
      DO 460 NUM=1,LT4
      IP(KN(NUM,2)-2*LT4)=NUM
  460 CONTINUE
      KEL=0
      NUM=0
      DO 565 IZ=1,LZ
      FRACT=1.0
      IF((NCODE(5).EQ.5).AND.(IZ.EQ.1)) FRACT=0.5
      IF((NCODE(6).EQ.5).AND.(IZ.EQ.LZ)) FRACT=0.5
      DZZ=ZZZ(IZ+1)-ZZZ(IZ)
      DO 510 JSTAGE=1,NBC
      DO 502 JEL=1,ISPLH
      DO 501 IRANG=1,NBC+JSTAGE-1
      DO 500 IEL=1,ISPLH
      KEL=KEL+1
      IF(IZGLOB(KEL,3).EQ.0) GO TO 500
      NUM=NUM+1
      IF((IRANG.EQ.1).AND.(IEL.EQ.1)) THEN
         LL1=.TRUE.
      ELSE
         LL1=(IZGLOB(KEL-1,3).EQ.0)
      ENDIF
      IF((IRANG.EQ.NBC+JSTAGE-1).AND.(IEL.EQ.ISPLH)) THEN
         LL2=.TRUE.
      ELSE
         LL2=(IZGLOB(KEL+1,3).EQ.0)
      ENDIF
      LCOUR=0
      DO 475 J=1,IELEM**2
      DO 470 I=1,IELEM+1
      LCOUR=LCOUR+1
      ITEMP = NDDIR
     >      + (JEL-1)*(2*(NBC+JSTAGE-1)*IELEM*ISPLH+1)*IELEM**2
     >      + (IRANG-1)*(2*IELEM*ISPLH)
     >      + (IEL-1)*IELEM
     >      + (J-1)*(2*(NBC+JSTAGE-1)*IELEM*ISPLH+1) + I
      IF(LCOUR.GT.NELEH) CALL XABORT('TRISFH: BUG9')
      IF(KEL.GT.NBLOS) CALL XABORT('TRISFH: BUG10')
      KN(IP(NUM),3+4*NELEH+LCOUR)=ITEMP
      KN(IP(NUM),3+5*NELEH+LCOUR)=ITEMP+IELEM*ISPLH
  470 CONTINUE
  475 CONTINUE
      IF(LL1) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 480 I=1,IELEM**2
            KN(IP(NUM),3+4*NELEH+(I-1)*(IELEM+1)+1)=0
  480       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),5)=SIDE*DZZ*FRACT/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),5)=SIDE*DZZ*FRACT
            IQFR(IP(NUM),5)=ICODE(1)
         ENDIF
      ENDIF
      IF(LL2) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 490 I=1,IELEM**2
            KN(IP(NUM),3+5*NELEH+I*(IELEM+1))=0
  490       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),6)=SIDE*DZZ*FRACT/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),6)=SIDE*DZZ*FRACT
            IQFR(IP(NUM),6)=ICODE(1)
         ENDIF
      ENDIF
  500 CONTINUE
  501 CONTINUE
  502 CONTINUE
      NDDIR=NDDIR+(NB1+2*(JSTAGE-1)*ISPLH*IELEM)*ISPLH*IELEM**2
  510 CONTINUE
*
      DO 560 JSTAGE=NBC+1,2*NBC-1
      DO 552 JEL=1,ISPLH
      DO 551 IRANG=1,(2*NBC-2)-(JSTAGE-NBC-1)
      DO 550 IEL=1,ISPLH
      KEL=KEL+1
      IF(IZGLOB(KEL,3).EQ.0) GO TO 550
      NUM=NUM+1
      IF((IRANG.EQ.1).AND.(IEL.EQ.1)) THEN
         LL1=.TRUE.
      ELSE
         LL1=(IZGLOB(KEL-1,3).EQ.0)
      ENDIF
      IF((IRANG.EQ.(2*NBC-2)-(JSTAGE-NBC-1)).AND.(IEL.EQ.ISPLH)) THEN
         LL2=.TRUE.
      ELSE
         LL2=(IZGLOB(KEL+1,3).EQ.0)
      ENDIF
      LCOUR=0
      DO 525 J=1,IELEM**2
      DO 520 I=1,IELEM+1
      LCOUR=LCOUR+1
      ITEMP = NDDIR
     >      + (JEL-1)*(2*(2*NBC-1+NBC-JSTAGE)*IELEM*ISPLH+1)*IELEM**2
     >      + (IRANG-1)*(2*IELEM*ISPLH)
     >      + (IEL-1)*IELEM
     >      + (J-1)*(2*(2*NBC-1+NBC-JSTAGE)*IELEM*ISPLH+1) + I
      IF(LCOUR.GT.NELEH) CALL XABORT('TRISFH: BUG11')
      IF(KEL.GT.NBLOS) CALL XABORT('TRISFH: BUG12')
      KN(IP(NUM),3+4*NELEH+LCOUR)=ITEMP
      KN(IP(NUM),3+5*NELEH+LCOUR)=ITEMP+IELEM*ISPLH
  520 CONTINUE
  525 CONTINUE
      IF(LL1) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 530 I=1,IELEM**2
            KN(IP(NUM),3+4*NELEH+(I-1)*(IELEM+1)+1)=0
  530       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),5)=SIDE*DZZ*FRACT/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),5)=SIDE*DZZ*FRACT
            IQFR(IP(NUM),5)=ICODE(1)
         ENDIF
      ENDIF
      IF(LL2) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 540 I=1,IELEM**2
            KN(IP(NUM),3+5*NELEH+I*(IELEM+1))=0
  540       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),6)=SIDE*DZZ*FRACT/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),6)=SIDE*DZZ*FRACT
            IQFR(IP(NUM),6)=ICODE(1)
         ENDIF
      ENDIF
  550 CONTINUE
  551 CONTINUE
  552 CONTINUE
      NDDIR=NDDIR+(NB2-2*(JSTAGE-NBC)*ISPLH*IELEM)*ISPLH*IELEM**2
  560 CONTINUE
  565 CONTINUE
*----
*  SET THE CURRENT NUMBERING PART OF MATRIX KN AND MATRIX QFR (Z-AXIS)
*----
      KEL=0
      NUM=0
      DO 635 IZ=1,LZ
      DO 630 IX=1,LXH*ISPLH**2
      KEL=KEL+1
      IF(IZGLOB(KEL,1).EQ.0) GO TO 630
      NUM=NUM+1
      IF(IZ.EQ.1) THEN
         LL1=.TRUE.
      ELSE
         LL1=(IZGLOB(KEL-LXH*ISPLH**2,1).EQ.0)
      ENDIF
      IF(IZ.EQ.LZ) THEN
         LL2=.TRUE.
      ELSE
         LL2=(IZGLOB(KEL+LXH*ISPLH**2,1).EQ.0)
      ENDIF
      DO 572 K=0,2 ! THREE LOZENGES PER HEXAGON
      DO 571 I=0,1 ! FACE ZINF/ZSUP
      DO 570 J=1,IELEM**2
      LCOUR=(2*K+I)*IELEM**2+J
      IF(LCOUR.GT.NELEZ) CALL XABORT('TRISFH: BUG11')
      IF(KEL.GT.NBLOS) CALL XABORT('TRISFH: BUG12')
      ITEMP = NDDIR
     >      + 3*(IX-1)*(LZ+1)*IELEM**2
     >      + K*(LZ+1)*IELEM**2
     >      + (J-1)*(LZ+1) + IZ + I
      KN(NUM,3+6*NELEH+LCOUR)=ITEMP
  570 CONTINUE
  571 CONTINUE
  572 CONTINUE
*
*     REFL OR ALBE BOUNDARY CONDITION
      IF(LL1) THEN
         COND=(NCODE(5).EQ.2).OR.((NCODE(5).EQ.1).AND.(ZCODE(5).EQ.1.0))
         IF(COND) THEN
            DO 585 K=0,2
            DO 580 J=1,IELEM**2
            LCINF=2*K*IELEM**2+J
            KN(NUM,3+6*NELEH+LCINF)=0
  580       CONTINUE
  585       CONTINUE
         ELSE IF((NCODE(5).EQ.1).AND.(ICODE(5).EQ.0)) THEN
            ALBEDO=0.5*(1.0-ZCODE(5))/(1.0+ZCODE(5))
            QFR(NUM,7)=0.8660254038*SIDE*SIDE/ALBEDO
         ELSE IF(NCODE(5).EQ.1) THEN
            QFR(NUM,7)=0.8660254038*SIDE*SIDE
            IQFR(NUM,7)=ICODE(5)
         ENDIF
      ENDIF
      IF(LL2) THEN
         COND=(NCODE(6).EQ.2).OR.((NCODE(6).EQ.1).AND.(ZCODE(6).EQ.1.0))
         IF(COND) THEN
            DO 595 K=0,2
            DO 590 J=1,IELEM**2
            LCSUP=(2*K+1)*IELEM**2+J
            KN(NUM,3+6*NELEH+LCSUP)=0
  590       CONTINUE
  595       CONTINUE
         ELSE IF((NCODE(6).EQ.1).AND.(ICODE(6).EQ.0)) THEN
            ALBEDO=0.5*(1.0-ZCODE(6))/(1.0+ZCODE(6))
            QFR(NUM,8)=0.8660254038*SIDE*SIDE/ALBEDO
         ELSE IF(NCODE(6).EQ.1) THEN
            QFR(NUM,8)=0.8660254038*SIDE*SIDE
            IQFR(NUM,8)=ICODE(6)
         ENDIF
      ENDIF
*     TRAN BOUNDARY CONDITION
      IF((IZ.EQ.LZ).AND.(NCODE(6).EQ.4)) THEN
         DO 605 K=0,2
         DO 600 J=1,IELEM**2
         LCSUP=(2*K+1)*IELEM**2+J
         KN(NUM,3+6*NELEH+LCSUP)=KN(NUM,3+6*NELEH+LCSUP)-LZ
  600    CONTINUE
  605    CONTINUE
      ENDIF
*     SYME BOUNDARY CONDITION
      IF((NCODE(5).EQ.5).AND.(IZ.EQ.1)) THEN
         QFR(NUM,7)=QFR(NUM,8)
         IQFR(NUM,7)=IQFR(NUM,8)
         DO 615 K=0,2
         DO 610 J=1,IELEM**2
         LCINF=2*K*IELEM**2+J
         LCSUP=(2*K+1)*IELEM**2+J
         KN(NUM,3+6*NELEH+LCINF)=-KN(NUM,3+6*NELEH+LCSUP)
  610    CONTINUE
  615    CONTINUE
      ELSE IF((NCODE(6).EQ.5).AND.(IZ.EQ.LZ)) THEN
         QFR(NUM,8)=QFR(NUM,7)
         IQFR(NUM,8)=IQFR(NUM,7)
         DO 625 K=0,2
         DO 620 J=1,IELEM**2
         LCINF=2*K*IELEM**2+J
         LCSUP=(2*K+1)*IELEM**2+J
         KN(NUM,3+6*NELEH+LCSUP)=-KN(NUM,3+6*NELEH+LCINF)
  620    CONTINUE
  625    CONTINUE
      ENDIF
  630 CONTINUE
  635 CONTINUE
*----
*  REMOVING THE UNUSED UNKNOWNS INDICES FROM KN
*----
      CALL XDISET(IP,3*LL4W0+LL4Z0,0)
      DO 645 KEL=1,LT4
      DO 640 ICOUR=1,6*NELEH+NELEZ
      IND=ABS(KN(KEL,3+ICOUR))
      IF(IND.GT.MAXIP) CALL XABORT('TRISFH: MAXIP OVERFLOW.')
      IF(IND.NE.0) IP(IND)=1
  640 CONTINUE
  645 CONTINUE
      LL4W=0
      DO 650 IND=1,LL4W0
      IF(IP(IND).EQ.1) THEN
         LL4W=LL4W+1
         IP(IND)=LL4W
      ENDIF
  650 CONTINUE
      LL4X=0
      DO 660 IND=1,LL4W0
      IF(IP(LL4W0+IND).EQ.1) THEN
         LL4X=LL4X+1
         IP(LL4W0+IND)=LL4W+LL4X
      ENDIF
  660 CONTINUE
      LL4Y=0
      DO 670 IND=1,LL4W0
      IF(IP(2*LL4W0+IND).EQ.1) THEN
         LL4Y=LL4Y+1
         IP(2*LL4W0+IND)=LL4W+LL4X+LL4Y
      ENDIF
  670 CONTINUE
      LL4Z=0
      DO 680 IND=1,LL4Z0
      IF(IP(3*LL4W0+IND).EQ.1) THEN
         LL4Z=LL4Z+1
         IP(3*LL4W0+IND)=LL4W+LL4X+LL4Y+LL4Z
      ENDIF
  680 CONTINUE
      DO 695 KEL=1,LT4
      DO 690 ICOUR=1,6*NELEH+NELEZ
      IF(KN(KEL,3+ICOUR).NE.0) THEN
         IND=KN(KEL,3+ICOUR)
         KN(KEL,3+ICOUR)=SIGN(IP(ABS(IND)),IND)
      ENDIF
  690 CONTINUE
  695 CONTINUE
      LL4=LL4F+LL4W+LL4X+LL4Y+LL4Z
*----
*  PRINT A FEW GEOMETRY CHARACTERISTICS
*----
      IF(IMPX.GT.0) THEN
         write(6,*) ' '
         write(6,*) 'ISPLH =',ISPLH
         write(6,*) 'IELEM =',IELEM
         write(6,*) 'NELEH =',NELEH
         write(6,*) 'NELEZ =',NELEZ
         write(6,*) 'NBLOS =',NBLOS
         write(6,*) 'LL4F  =',LL4F
         write(6,*) 'LL4W  =',LL4W
         write(6,*) 'LL4X  =',LL4X
         write(6,*) 'LL4Y  =',LL4Y
         write(6,*) 'LL4Z  =',LL4Z
         write(6,*) 'NBC   =',NBC
      ENDIF
*----
*  SET IPERT
*----
      KEL=0
      DO 714 IZ=1,LZ
      DO 703 JSTAGE=1,NBC
      DO 702 JEL=1,ISPLH
      DO 701 IRANG=1,NBC+JSTAGE-1
      DO 700 IEL=1,ISPLH
      KEL=KEL+1
      IHEX=IZGLOB(KEL,1)
      IF(IHEX.EQ.0) THEN
        IPERT(KEL)=0
      ELSE
        IPERT(KEL)=(IHEX-1)*ISPLH**2+(IEL-1)*ISPLH+JEL
      ENDIF
      IF(IPERT(KEL).GT.NBLOS) call XABORT('TRISFH: NBLOS OVERFLOW(1)')
  700 CONTINUE
  701 CONTINUE
  702 CONTINUE
  703 CONTINUE
      DO 713 JSTAGE=NBC+1,2*NBC-1
      DO 712 JEL=1,ISPLH
      DO 711 IRANG=1,(2*NBC-2)-(JSTAGE-NBC-1)
      DO 710 IEL=1,ISPLH
      KEL=KEL+1
      IHEX=IZGLOB(KEL,1)
      IF(IHEX.EQ.0) THEN
        IPERT(KEL)=0
      ELSE
        IPERT(KEL)=(IHEX-1)*ISPLH**2+(IEL-1)*ISPLH+JEL
      ENDIF
      IF(IPERT(KEL).GT.NBLOS) call XABORT('TRISFH: NBLOS OVERFLOW(2)')
  710 CONTINUE
  711 CONTINUE
  712 CONTINUE
  713 CONTINUE
  714 CONTINUE
      IF(KEL.NE.NBLOS) CALL XABORT('TRISFH: IPERT FAILURE.')
*----
*  SET IDL, VOL, FRZ AND ZZ
*----
      NUM=0
      IDL(:3,:NBLOS)=0
      VOL(:3,:NBLOS)=0.0
      FRZ(:NBLOS)=0.0
      ZZ(:3,:NBLOS)=0.0
      DO 725 IZ=1,LZ
      FRACT=1.0
      IF((NCODE(5).EQ.5).AND.(IZ.EQ.1)) FRACT=0.5
      IF((NCODE(6).EQ.5).AND.(IZ.EQ.LZ)) FRACT=0.5
      DZ=ZZZ(IZ+1)-ZZZ(IZ)
      DO 720 J=1,LXH*ISPLH**2
      KEL=(IZ-1)*LXH*ISPLH**2+J
      KEL2=IPERT(KEL)
      IF(KEL2.EQ.0) GO TO 720
      NUM=NUM+1
      IDL(1,KEL2)=(NUM-1)*IELEM**3+1
      IDL(2,KEL2)=(KN(NUM,1)-1)*IELEM**3+1
      IDL(3,KEL2)=(KN(NUM,2)-1)*IELEM**3+1
      VOL(:3,KEL2)=2.59807587*SIDE*SIDE*DZ*FRACT/REAL(3)
      FRZ(KEL)=FRACT
      ZZ(:3,KEL2)=DZ
  720 CONTINUE
  725 CONTINUE
      IF(IMPX.GT.2) THEN
         WRITE(6,790) 'MAT',(((MAT(I,J,K),I=1,3),J=1,ISPLH**2),
     1   K=1,LXH*LZ)
         WRITE(6,790) 'IDL',((IDL(I,J),I=1,3),J=1,NBLOS)
         WRITE(6,800) 'ZZ ',((ZZ(I,J),I=1,3),J=1,NBLOS)
         WRITE(6,800) 'VOL',((VOL(I,J),I=1,3),J=1,NBLOS)
      ENDIF
*
      IF(IMPX.GT.0) WRITE(6,810) LL4
      IF(IMPX.GT.2) THEN
         WRITE (6,830)
         DO 730 K=1,NBLOS
         WRITE (6,840) K,(IZGLOB(K,I),I=1,3)
  730    CONTINUE
         WRITE (6,850)
         DO 740 K=1,LT4
         WRITE (6,860) K,(KN(K,I),I=1,3+2*NELEH)
         WRITE (6,870) 'X',(KN(K,I),I=3+2*NELEH+1,3+4*NELEH)
         WRITE (6,870) 'Y',(KN(K,I),I=3+4*NELEH+1,3+6*NELEH)
         IF(LL4Z.GT.0) THEN
            WRITE (6,870) 'Z',(KN(K,I),I=3+6*NELEH+1,3+6*NELEH+NELEZ)
         ENDIF
  740    CONTINUE
         WRITE (6,880)
         DO 750 K=1,LT4
         WRITE (6,890) K,(QFR(K,I),I=1,8)
  750    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IZGLOB,IP,IJP)
      RETURN
*
  790 FORMAT(1X,A3/14(2X,I6))
  800 FORMAT(1X,A3/7(2X,E12.5))
  810 FORMAT(31H NUMBER OF UNKNOWNS PER GROUP =,I8)
  830 FORMAT(/22H NUMBERING OF HEXAGONS/1X,21(1H-)//8H ELEMENT,4X,
     1 24H W ----- X ----- Y -----)
  840 FORMAT(1X,I6,5X,3I8)
  850 FORMAT(/22H NUMBERING OF UNKNOWNS/1X,21(1H-)//8H ELEMENT,5X,
     1 20H---> X ---> Y ---> W,4X,8HCURRENTS,89(1H.))
  860 FORMAT(1X,I6,5X,3I7,4X,1HW,12I8:/(38X,12I8))
  870 FORMAT(37X,A1,12I8:/(38X,12I8))
  880 FORMAT(/8H ELEMENT,3X,23HVOID BOUNDARY CONDITION/15X,7(1H-),
     1 3H W ,7(1H-),3X,7(1H-),3H X ,7(1H-),3X,7(1H-),3H Y ,7(1H-),
     2 3X,7(1H-),3H Z ,7(1H-))
  890 FORMAT(1X,I6,5X,1P,10E10.1/(12X,1P,10E10.1))
      END
