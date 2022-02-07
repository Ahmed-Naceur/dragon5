*DECK BIVSFH
      SUBROUTINE BIVSFH (MAXEV,NBLOS,IMPX,ISPLH,IELEM,LXH,MAT,SIDE,
     1 NCODE,ICODE,ZCODE,LL4,VOL,IDL,IPERT,KN,QFR,IQFR,BFR,MU)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a Thomas-Raviart-Schneider finite element
* discretization of a 2-D hexagonal geometry.
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
* MAXEV   allocated storage for vector MU.
* NBLOS   number of lozenges per direction, taking into account
*         mesh-splitting.
* IMPX    print parameter.
* ISPLH   mesh-splitting in 3*ISPLH**2 lozenges per hexagon.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic).
* LXH     number of hexagons.
* MAT     mixture index assigned to each lozenge.
* SIDE    side of a lozenge.
* NCODE   type of boundary condition applied on each side (I=1: hbc):
*         NCODE(I)=1: VOID;          =2: REFL;       =6: ALBE;
*                 =5: SYME;          =7: ZERO.
* ICODE   physical albedo index on each side of the domain.
* ZCODE   albedo corresponding to boundary condition 'VOID' on each
*         side (ZCODE(I)=0.0 by default).
*
*Parameters: output
* LL4     order of the system matrices.
* VOL     volume of each lozenge.
* IDL     position of the average flux component associated with each
*         lozenge.
* IPERT   mixture permutation index.
* KN      ADI permutation indices for the volumes and currents.
* QFR     element-ordered boundary conditions.
* IQFR    element-ordered physical albedo indices.
* BFR     element-ordered surface fractions.
* MU      compressed storage mode indices.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXEV,NBLOS,IMPX,ISPLH,IELEM,LXH,MAT(3,ISPLH**2,LXH),
     1 NCODE(4),ICODE(4),LL4,IDL(3,NBLOS),IPERT(NBLOS),
     2 KN(NBLOS,4+6*IELEM*(IELEM+1)),IQFR(NBLOS,6),MU(MAXEV)
      REAL SIDE,ZCODE(4),VOL(3,NBLOS),QFR(NBLOS,6),BFR(NBLOS,6)
*----
*  LOCAL VARIABLES
*----
      LOGICAL COND,LL1,LL2
      INTEGER, DIMENSION(:),ALLOCATABLE :: IP,I1,I3,I4,I5
      INTEGER, DIMENSION(:,:),ALLOCATABLE :: IZGLOB
      INTEGER, DIMENSION(:,:,:),ALLOCATABLE :: IJP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJP(LXH,ISPLH,ISPLH),IP(MAXEV),IZGLOB(NBLOS,3))
*----
*  THOMAS-RAVIART-SCHNEIDER SPECIFIC NUMEROTATION
*----
      NBC=INT((SQRT(REAL((4*LXH-1)/3))+1.)/2.)
      IF(LXH.NE.1+3*NBC*(NBC-1)) CALL XABORT('BIVSFH: INVALID VALUE OF'
     1 //' LXH(1).')
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
         IF(I.NE.LXH) CALL XABORT('BIVSFH: INVALID VALUE OF LXH(2)')
      ENDIF
      ALLOCATE(I1(3*LXH),I3(2*LXH),I4(NBLOS),I5(NBLOS))
      DO 30 I=1,LXH
      I3(I)=I
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
      DO 42 I=1,LXH
      IOFW=I1(I)
      IOFX=I1(LXH+I)
      IOFY=I1(2*LXH+I)
      DO 41 IM=1,ISPLH
      DO 40 JM=1,ISPLH
      IZGLOB(IJP(IOFW,IM,JM),1)=I4(I)
      IZGLOB(IJP(IOFX,IM,JM),2)=I4(I)
      IZGLOB(IJP(IOFY,IM,JM),3)=I4(I)
   40 CONTINUE
   41 CONTINUE
   42 CONTINUE
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
      CALL XDISET(KN,NBLOS*(4+6*IELEM*(IELEM+1)),0)
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
         DO 90 I=1,LXH
         IF(IZGLOB(I,1).EQ.0) GO TO 90
         KN(I4(I),2)=I5(I3(I))+LT4
   90    CONTINUE
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
         II2=ISPLH*(JKEEP*ISPLH+(ISPLH-JM)*KMAX+K0-1)+IM
         IF(IZGLOB(II1,1).EQ.0) GO TO 120
         KN(I4(II1),2)=I5(II2)+LT4
  110    CONTINUE
  115    CONTINUE
  120    CONTINUE
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
         DO 140 I=1,LXH
         IF(IZGLOB(I,1).EQ.0) GO TO 140
         KN(I4(I),3)=I5(I3(LXH+I))+2*LT4
  140    CONTINUE
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
         II2=ISPLH*(JKEEP*ISPLH+(ISPLH-IM)*KMAX+K0-1)+(ISPLH-JM+1)
         IF(IZGLOB(II1,1).EQ.0) GO TO 170
         KN(I4(II1),3)=I5(II2)+2*LT4
  160    CONTINUE
  165    CONTINUE
  170    CONTINUE
      ENDIF
*----
*  COMPUTE THE FLUX PERMUTATION PART OF MATRIX KN (Y <--> W)
*----
      IF(ISPLH.EQ.1) THEN
         DO 180 I=1,LXH
         IF(IZGLOB(I,1).EQ.0) GO TO 180
         KN(I4(I),4)=I4(I)
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
         DO 210 I=1,LXH
         JMAX=I1(I)
         IKEEP=I1(LXH+I)
         J00=I1(2*LXH+I)
         DO 205 IM=1,ISPLH
         DO 200 JM=1,ISPLH
         II1=ISPLH*(IKEEP*ISPLH+(IM-1)*JMAX+J00-1)+JM
         II2=ISPLH*(IKEEP*ISPLH+(JM-1)*JMAX+J00-1)+(ISPLH-IM+1)
         IF(IZGLOB(II1,1).EQ.0) GO TO 210
         KN(I4(II1),4)=I4(II2)
  200    CONTINUE
  205    CONTINUE
  210    CONTINUE
      ENDIF
      DEALLOCATE(I5,I4,I3,I1)
*----
*  SET THE CURRENT NUMBERING PART OF MATRIX KN AND MATRIX QFR (W-AXIS)
*----
      LL4W0=(2*NBLOS*IELEM+(2*NBC-1)*ISPLH)*IELEM
      LL4F=3*LT4*IELEM*IELEM
      CALL XDRSET(QFR,6*NBLOS,0.0)
      CALL XDISET(IQFR,6*NBLOS,0)
      CALL XDRSET(BFR,6*NBLOS,0.0)
      ALBEDO=0.5*(1.0-ZCODE(1))/(1.0+ZCODE(1))
      NELEM=IELEM*(IELEM+1)
      NB1=(2*NBC*IELEM*ISPLH+1)*IELEM*ISPLH
      KEL=0
      NDDIR=LL4F
      NUM=0
      DO 290 JSTAGE=1,NBC
      DO 282 JEL=1,ISPLH
      DO 281 IRANG=1,NBC+JSTAGE-1
      DO 280 IEL=1,ISPLH
      KEL=KEL+1
      IF(IZGLOB(KEL,1).EQ.0) GO TO 280
      NUM=NUM+1
      KN(NUM,1)=NUM
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
      DO 255 J=1,IELEM
      DO 250 I=1,IELEM+1
      LCOUR=LCOUR+1
      ITEMP = NDDIR
     >      + (JEL-1)*(2*(NBC+JSTAGE-1)*IELEM*ISPLH+1)*IELEM
     >      + (IRANG-1)*(2*IELEM*ISPLH)
     >      + (IEL-1)*IELEM
     >      + (J-1)*(2*(NBC+JSTAGE-1)*IELEM*ISPLH+1) + I
      IF(LCOUR.GT.NELEM) CALL XABORT('BIVSFH: bug1')
      IF(KEL.GT.NBLOS) CALL XABORT('BIVSFH: bug2')
      KN(NUM,4+LCOUR)=ITEMP
      KN(NUM,4+NELEM+LCOUR)=ITEMP+IELEM*ISPLH
  250 CONTINUE
  255 CONTINUE
      IF(LL1) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 260 I=1,IELEM
            KN(NUM,4+(I-1)*(IELEM+1)+1)=0
  260       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(NUM,1)=SIDE/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(NUM,1)=SIDE
            IQFR(NUM,1)=ICODE(1)
         ENDIF
         IF((NCODE(1).EQ.1).OR.(NCODE(1).EQ.7)) BFR(NUM,1)=SIDE
      ENDIF
      IF(LL2) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 270 I=1,IELEM
            KN(NUM,4+NELEM+I*(IELEM+1))=0
  270       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(NUM,2)=SIDE/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(NUM,2)=SIDE
            IQFR(NUM,2)=ICODE(1)
         ENDIF
         IF((NCODE(1).EQ.1).OR.(NCODE(1).EQ.7)) BFR(NUM,2)=SIDE
      ENDIF
  280 CONTINUE
  281 CONTINUE
  282 CONTINUE
      NDDIR=NDDIR+NB1+(2*(JSTAGE-1)*IELEM*ISPLH)*IELEM*ISPLH
  290 CONTINUE
*
      DO 340 JSTAGE=NBC+1,2*NBC-1
      DO 332 JEL=1,ISPLH
      DO 331 IRANG=1,(2*NBC-2)-(JSTAGE-NBC-1)
      DO 330 IEL=1,ISPLH
      KEL=KEL+1
      IF(IZGLOB(KEL,1).EQ.0) GO TO 330
      NUM=NUM+1
      KN(NUM,1)=NUM
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
      DO 305 J=1,IELEM
      DO 300 I=1,IELEM+1
      LCOUR=LCOUR+1
      ITEMP = NDDIR
     >      + (JEL-1)*(2*(2*NBC-1+NBC-JSTAGE)*IELEM*ISPLH+1)*IELEM
     >      + (IRANG-1)*(2*IELEM*ISPLH)
     >      + (IEL-1)*IELEM
     >      + (J-1)*(2*(2*NBC-1+NBC-JSTAGE)*IELEM*ISPLH+1) + I
      IF(LCOUR.GT.NELEM) CALL XABORT('BIVSFH: bug3')
      IF(KEL.GT.NBLOS) CALL XABORT('BIVSFH: bug4')
      KN(NUM,4+LCOUR)=ITEMP
      KN(NUM,4+NELEM+LCOUR)=ITEMP+IELEM*ISPLH
  300 CONTINUE
  305 CONTINUE
      IF(LL1) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 310 I=1,IELEM
            KN(NUM,4+(I-1)*(IELEM+1)+1)=0
  310       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(NUM,1)=SIDE/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(NUM,1)=SIDE
            IQFR(NUM,1)=ICODE(1)
         ENDIF
         IF((NCODE(1).EQ.1).OR.(NCODE(1).EQ.7)) BFR(NUM,1)=SIDE
      ENDIF
      IF(LL2) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 320 I=1,IELEM
            KN(NUM,4+NELEM+I*(IELEM+1))=0
  320       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(NUM,2)=SIDE/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(NUM,2)=SIDE
            IQFR(NUM,2)=ICODE(1)
         ENDIF
         IF((NCODE(1).EQ.1).OR.(NCODE(1).EQ.7)) BFR(NUM,2)=SIDE
      ENDIF
  330 CONTINUE
  331 CONTINUE
  332 CONTINUE
      NDDIR=NDDIR+(2*(2*NBC-1)*IELEM*ISPLH+1)*IELEM*ISPLH
     >           -(2*(JSTAGE-NBC)*IELEM*ISPLH)*IELEM*ISPLH
  340 CONTINUE
*----
*  SET THE CURRENT NUMBERING PART OF MATRIX KN AND MATRIX QFR (X-AXIS)
*----
      CALL XDISET(IP,NBLOS,0)
      DO 350 NUM=1,LT4
      IP(KN(NUM,2)-LT4)=NUM
  350 CONTINUE
      KEL=0
      NUM=0
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
      DO 365 J=1,IELEM
      DO 360 I=1,IELEM+1
      LCOUR=LCOUR+1
      ITEMP = NDDIR
     >      + (JEL-1)*(2*(NBC+JSTAGE-1)*IELEM*ISPLH+1)*IELEM
     >      + (IRANG-1)*(2*IELEM*ISPLH)
     >      + (IEL-1)*IELEM
     >      + (J-1)*(2*(NBC+JSTAGE-1)*IELEM*ISPLH+1) + I
      IF(LCOUR.GT.NELEM) CALL XABORT('BIVSFH: bug5')
      IF(KEL.GT.NBLOS) CALL XABORT('BIVSFH: bug6')
      KN(IP(NUM),4+2*NELEM+LCOUR)=ITEMP
      KN(IP(NUM),4+3*NELEM+LCOUR)=ITEMP+IELEM*ISPLH
  360 CONTINUE
  365 CONTINUE
      IF(LL1) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 370 I=1,IELEM
            KN(IP(NUM),4+2*NELEM+(I-1)*(IELEM+1)+1)=0
  370       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),3)=SIDE/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),3)=SIDE
            IQFR(NUM,3)=ICODE(1)
         ENDIF
         IF((NCODE(1).EQ.1).OR.(NCODE(1).EQ.7)) BFR(NUM,3)=SIDE
      ENDIF
      IF(LL2) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 380 I=1,IELEM
            KN(IP(NUM),4+3*NELEM+I*(IELEM+1))=0
  380       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),4)=SIDE/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),4)=SIDE
            IQFR(NUM,4)=ICODE(1)
         ENDIF
         IF((NCODE(1).EQ.1).OR.(NCODE(1).EQ.7)) BFR(NUM,4)=SIDE
      ENDIF
  390 CONTINUE
  391 CONTINUE
  392 CONTINUE
      NDDIR=NDDIR+NB1+(2*(JSTAGE-1)*IELEM*ISPLH)*IELEM*ISPLH
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
      DO 415 J=1,IELEM
      DO 410 I=1,IELEM+1
      LCOUR=LCOUR+1
      ITEMP = NDDIR
     >      + (JEL-1)*(2*(2*NBC-1+NBC-JSTAGE)*IELEM*ISPLH+1)*IELEM
     >      + (IRANG-1)*(2*IELEM*ISPLH)
     >      + (IEL-1)*IELEM
     >      + (J-1)*(2*(2*NBC-1+NBC-JSTAGE)*IELEM*ISPLH+1) + I
      IF(LCOUR.GT.NELEM) CALL XABORT('BIVSFH: bug7')
      IF(KEL.GT.NBLOS) CALL XABORT('BIVSFH: bug8')
      KN(IP(NUM),4+2*NELEM+LCOUR)=ITEMP
      KN(IP(NUM),4+3*NELEM+LCOUR)=ITEMP+IELEM*ISPLH
  410 CONTINUE
  415 CONTINUE
      IF(LL1) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 420 I=1,IELEM
            KN(IP(NUM),4+2*NELEM+(I-1)*(IELEM+1)+1)=0
  420       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),3)=SIDE/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),3)=SIDE
            IQFR(NUM,3)=ICODE(1)
         ENDIF
         IF((NCODE(1).EQ.1).OR.(NCODE(1).EQ.7)) BFR(NUM,3)=SIDE
      ENDIF
      IF(LL2) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 430 I=1,IELEM
            KN(IP(NUM),4+3*NELEM+I*(IELEM+1))=0
  430       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),4)=SIDE/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),4)=SIDE
            IQFR(NUM,4)=ICODE(1)
         ENDIF
         IF((NCODE(1).EQ.1).OR.(NCODE(1).EQ.7)) BFR(NUM,4)=SIDE
      ENDIF
  440 CONTINUE
  441 CONTINUE
  442 CONTINUE
      NDDIR=NDDIR+(2*(2*NBC-1)*IELEM*ISPLH+1)*IELEM*ISPLH
     >           -(2*(JSTAGE-NBC)*IELEM*ISPLH)*IELEM*ISPLH
  450 CONTINUE
*----
*  SET THE CURRENT NUMBERING PART OF MATRIX KN AND MATRIX QFR (Y-AXIS)
*----
      CALL XDISET(IP,NBLOS,0)
      DO 460 NUM=1,LT4
      IP(KN(NUM,3)-2*LT4)=NUM
  460 CONTINUE
      KEL=0
      NUM=0
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
      DO 475 J=1,IELEM
      DO 470 I=1,IELEM+1
      LCOUR=LCOUR+1
      ITEMP = NDDIR
     >      + (JEL-1)*(2*(NBC+JSTAGE-1)*IELEM*ISPLH+1)*IELEM
     >      + (IRANG-1)*(2*IELEM*ISPLH)
     >      + (IEL-1)*IELEM
     >      + (J-1)*(2*(NBC+JSTAGE-1)*IELEM*ISPLH+1) + I
      IF(LCOUR.GT.NELEM) CALL XABORT('BIVSFH: bug9')
      IF(KEL.GT.NBLOS) CALL XABORT('BIVSFH: bug10')
      KN(IP(NUM),4+4*NELEM+LCOUR)=ITEMP
      KN(IP(NUM),4+5*NELEM+LCOUR)=ITEMP+IELEM*ISPLH
  470 CONTINUE
  475 CONTINUE
      IF(LL1) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 480 I=1,IELEM
            KN(IP(NUM),4+4*NELEM+(I-1)*(IELEM+1)+1)=0
  480       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),5)=SIDE/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),5)=SIDE
            IQFR(NUM,5)=ICODE(1)
         ENDIF
         IF((NCODE(1).EQ.1).OR.(NCODE(1).EQ.7)) BFR(NUM,5)=SIDE
      ENDIF
      IF(LL2) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 490 I=1,IELEM
            KN(IP(NUM),4+5*NELEM+I*(IELEM+1))=0
  490       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),6)=SIDE/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),6)=SIDE
            IQFR(NUM,6)=ICODE(1)
         ENDIF
         IF((NCODE(1).EQ.1).OR.(NCODE(1).EQ.7)) BFR(NUM,6)=SIDE
      ENDIF
  500 CONTINUE
  501 CONTINUE
  502 CONTINUE
      NDDIR=NDDIR+NB1+(2*(JSTAGE-1)*IELEM*ISPLH)*IELEM*ISPLH
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
      DO 525 J=1,IELEM
      DO 520 I=1,IELEM+1
      LCOUR=LCOUR+1
      ITEMP = NDDIR
     >      + (JEL-1)*(2*(2*NBC-1+NBC-JSTAGE)*IELEM*ISPLH+1)*IELEM
     >      + (IRANG-1)*(2*IELEM*ISPLH)
     >      + (IEL-1)*IELEM
     >      + (J-1)*(2*(2*NBC-1+NBC-JSTAGE)*IELEM*ISPLH+1) + I
      IF(LCOUR.GT.NELEM) CALL XABORT('BIVSFH: bug11')
      IF(KEL.GT.NBLOS) CALL XABORT('BIVSFH: bug12')
      KN(IP(NUM),4+4*NELEM+LCOUR)=ITEMP
      KN(IP(NUM),4+5*NELEM+LCOUR)=ITEMP+IELEM*ISPLH
  520 CONTINUE
  525 CONTINUE
      IF(LL1) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 530 I=1,IELEM
            KN(IP(NUM),4+4*NELEM+(I-1)*(IELEM+1)+1)=0
  530       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),5)=SIDE/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),5)=SIDE
            IQFR(NUM,5)=ICODE(1)
         ENDIF
         IF((NCODE(1).EQ.1).OR.(NCODE(1).EQ.7)) BFR(NUM,5)=SIDE
      ENDIF
      IF(LL2) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 540 I=1,IELEM
            KN(IP(NUM),4+5*NELEM+I*(IELEM+1))=0
  540       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(IP(NUM),6)=SIDE/ALBEDO
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(IP(NUM),6)=SIDE
            IQFR(NUM,6)=ICODE(1)
         ENDIF
         IF((NCODE(1).EQ.1).OR.(NCODE(1).EQ.7)) BFR(NUM,6)=SIDE
      ENDIF
  550 CONTINUE
  551 CONTINUE
  552 CONTINUE
      NDDIR=NDDIR+(2*(2*NBC-1)*IELEM*ISPLH+1)*IELEM*ISPLH
     >           -(2*(JSTAGE-NBC)*IELEM*ISPLH)*IELEM*ISPLH
  560 CONTINUE
*----
*  COMPUTE THE SURFACE FRACTIONS
*----
      SURFTOT=0.0
      DO 566 I=1,NBLOS
      DO 565 J=1,6
      SURFTOT=SURFTOT+BFR(I,J)
  565 CONTINUE
  566 CONTINUE
      IF(SURFTOT.GT.0.0) THEN
         DO 575 I=1,NBLOS
         DO 570 J=1,6
         BFR(I,J)=BFR(I,J)/SURFTOT
  570    CONTINUE
  575    CONTINUE
      ENDIF
*----
*  REORDER THE UNKNOWNS AND REMOVE THE UNUSED UNKNOWNS INDICES FROM KN
*----
      CALL XDISET(IP,LL4F+3*LL4W0,0)
      LL4=0
      DO 591 KEL=1,LT4
      DO 582 IFLUX=1,4
      NUM=KN(KEL,IFLUX)
      DO 581 K2=1,IELEM
      DO 580 K1=1,IELEM
      JND1=(NUM-1)*IELEM**2+(K2-1)*IELEM+K1
      IF(JND1.GT.MAXEV) CALL XABORT('BIVSFH: MAXEV OVERFLOW(1).')
      IF(IP(JND1).EQ.0) THEN
         LL4=LL4+1
         IP(JND1)=LL4
      ENDIF
  580 CONTINUE
  581 CONTINUE
  582 CONTINUE
      DO 590 ICOUR=1,6*NELEM
      IND=ABS(KN(KEL,4+ICOUR))
      IF(IND.GT.MAXEV) CALL XABORT('BIVSFH: MAXEV OVERFLOW(2).')
      IF(IND.NE.0) THEN
         IF(IP(IND).EQ.0) THEN
            LL4=LL4+1
            IP(IND)=LL4
         ENDIF
      ENDIF
  590 CONTINUE
  591 CONTINUE
      DO 605 KEL=1,LT4
      DO 595 IFLUX=1,4
      NUM=KN(KEL,IFLUX)
      KN(KEL,IFLUX)=IP((NUM-1)*IELEM**2+1)
  595 CONTINUE
      DO 600 ICOUR=1,6*NELEM
      IF(KN(KEL,4+ICOUR).NE.0) THEN
         IND=KN(KEL,4+ICOUR)
         KN(KEL,4+ICOUR)=SIGN(IP(ABS(IND)),IND)
      ENDIF
  600 CONTINUE
  605 CONTINUE
*----
*  PRINT A FEW GEOMETRY CHARACTERISTICS
*----
      IF(IMPX.GT.0) THEN
         write(6,*) ' '
         write(6,*) 'ISPLH =',ISPLH
         write(6,*) 'IELEM =',IELEM
         write(6,*) 'NELEM =',NELEM
         write(6,*) 'NBLOS =',NBLOS
         write(6,*) 'LL4F  =',LL4F
         write(6,*) 'LL4   =',LL4
         write(6,*) 'NBC   =',NBC
      ENDIF
*----
*  SET IPERT
*----
      KEL=0
      DO 613 JSTAGE=1,NBC
      DO 612 JEL=1,ISPLH
      DO 611 IRANG=1,NBC+JSTAGE-1
      DO 610 IEL=1,ISPLH
      KEL=KEL+1
      IHEX=IZGLOB(KEL,1)
      IF(IHEX.EQ.0) THEN
        IPERT(KEL)=0
      ELSE
        IPERT(KEL)=(IHEX-1)*ISPLH**2+(IEL-1)*ISPLH+JEL
      ENDIF
  610 CONTINUE
  611 CONTINUE
  612 CONTINUE
  613 CONTINUE
      DO 623 JSTAGE=NBC+1,2*NBC-1
      DO 622 JEL=1,ISPLH
      DO 621 IRANG=1,(2*NBC-2)-(JSTAGE-NBC-1)
      DO 620 IEL=1,ISPLH
      KEL=KEL+1
      IHEX=IZGLOB(KEL,1)
      IF(IHEX.EQ.0) THEN
        IPERT(KEL)=0
      ELSE
        IPERT(KEL)=(IHEX-1)*ISPLH**2+(IEL-1)*ISPLH+JEL
      ENDIF
  620 CONTINUE
  621 CONTINUE
  622 CONTINUE
  623 CONTINUE
      IF(KEL.NE.NBLOS) CALL XABORT('BIVSFH: IPERT FAILURE.')
*----
*  SET IDL AND VOL
*----
      NUM=0
      IDL(:3,:NBLOS)=0
      VOL(:3,:NBLOS)=0.0
      DO 630 KEL=1,NBLOS
      KEL2=IPERT(KEL)
      IF(KEL2.EQ.0) GO TO 630
      NUM=NUM+1
      IDL(:3,KEL2)=KN(NUM,:3)
      VOL(:3,KEL2)=2.59807587*SIDE*SIDE/REAL(3)
  630 CONTINUE
      IF(IMPX.GT.2) THEN
         WRITE(6,800) 'MAT',(((MAT(I,J,K),I=1,3),J=1,ISPLH**2),K=1,LXH)
         WRITE(6,800) 'IDL',((IDL(I,J),I=1,3),J=1,NBLOS)
         WRITE(6,810) 'VOL',((VOL(I,J),I=1,3),J=1,NBLOS)
      ENDIF
*----
*  COMPUTE THE SYSTEM MATRIX BANDWIDTH.
*----
      CALL XDISET(MU,LL4,1)
      NUM=0
      DO 690 KEL=1,NBLOS
      IF(IZGLOB(KEL,1).EQ.0) GO TO 690
      NUM=NUM+1
      DO 663 K4=0,1
      DO 662 K3=0,IELEM-1
      DO 661 K2=1,IELEM+1
      INW1=ABS(KN(NUM,4+K4*NELEM+K3*(IELEM+1)+K2))
      INX1=ABS(KN(NUM,4+(K4+2)*NELEM+K3*(IELEM+1)+K2))
      INY1=ABS(KN(NUM,4+(K4+4)*NELEM+K3*(IELEM+1)+K2))
      DO 650 K1=1,IELEM+1
      INW2=ABS(KN(NUM,4+K4*NELEM+K3*(IELEM+1)+K1))
      INX2=ABS(KN(NUM,4+(K4+2)*NELEM+K3*(IELEM+1)+K1))
      INY2=ABS(KN(NUM,4+(K4+4)*NELEM+K3*(IELEM+1)+K1))
      IF((INW2.NE.0).AND.(INW1.NE.0)) THEN
         MU(INW1)=MAX(MU(INW1),INW1-INW2+1)
         MU(INW2)=MAX(MU(INW2),INW2-INW1+1)
      ENDIF
      IF((INX2.NE.0).AND.(INX1.NE.0)) THEN
         MU(INX1)=MAX(MU(INX1),INX1-INX2+1)
         MU(INX2)=MAX(MU(INX2),INX2-INX1+1)
      ENDIF
      IF((INY2.NE.0).AND.(INY1.NE.0)) THEN
         MU(INY1)=MAX(MU(INY1),INY1-INY2+1)
         MU(INY2)=MAX(MU(INY2),INY2-INY1+1)
      ENDIF
  650 CONTINUE
      DO 660 K1=0,IELEM-1
      IF(K4.EQ.0) THEN
         JND1=KN(NUM,1)+K3*IELEM+K1
         JND2=KN(NUM,2)+K3*IELEM+K1
         JND3=KN(NUM,3)+K3*IELEM+K1
      ELSE
         JND1=KN(NUM,2)+K1*IELEM+K3
         JND2=KN(NUM,3)+K1*IELEM+K3
         JND3=KN(NUM,4)+K1*IELEM+K3
      ENDIF
      IF(INW1.NE.0) THEN
         MU(JND1)=MAX(MU(JND1),JND1-INW1+1)
         MU(INW1)=MAX(MU(INW1),INW1-JND1+1)
      ENDIF
      IF(INX1.NE.0) THEN
         MU(JND2)=MAX(MU(JND2),JND2-INX1+1)
         MU(INX1)=MAX(MU(INX1),INX1-JND2+1)
      ENDIF
      IF(INY1.NE.0) THEN
         MU(JND3)=MAX(MU(JND3),JND3-INY1+1)
         MU(INY1)=MAX(MU(INY1),INY1-JND3+1)
      ENDIF
  660 CONTINUE
  661 CONTINUE
  662 CONTINUE
  663 CONTINUE
      ITRS=0
      DO I=1,LT4
         IF(KN(I,1).EQ.KN(NUM,4)) THEN
            ITRS=I
            GO TO 670
         ENDIF
      ENDDO
      CALL XABORT('BIVSFH: ITRS FAILURE.')
  670 DO 685 I=1,NELEM
      INW1=ABS(KN(ITRS,4+I))
      INX1=ABS(KN(NUM,4+2*NELEM+I))
      INY1=ABS(KN(NUM,4+4*NELEM+I))
      DO 680 J=1,NELEM
      INW2=ABS(KN(NUM,4+NELEM+J))
      INX2=ABS(KN(NUM,4+3*NELEM+J))
      INY2=ABS(KN(NUM,4+5*NELEM+J))
      IF((INY2.NE.0).AND.(INW1.NE.0)) THEN
         MU(INW1)=MAX(MU(INW1),INW1-INY2+1)
         MU(INY2)=MAX(MU(INY2),INY2-INW1+1)
      ENDIF
      IF((INW2.NE.0).AND.(INX1.NE.0)) THEN
         MU(INX1)=MAX(MU(INX1),INX1-INW2+1)
         MU(INW2)=MAX(MU(INW2),INW2-INX1+1)
      ENDIF
      IF((INX2.NE.0).AND.(INY1.NE.0)) THEN
         MU(INY1)=MAX(MU(INY1),INY1-INX2+1)
         MU(INX2)=MAX(MU(INX2),INX2-INY1+1)
      ENDIF
  680 CONTINUE
  685 CONTINUE
  690 CONTINUE
      MUMAX=0
      IIMAX=0
      DO 700 I=1,LL4
      MUMAX=MAX(MUMAX,MU(I))
      IIMAX=IIMAX+MU(I)
      MU(I)=IIMAX
  700 CONTINUE
*
      IF(IMPX.GT.0) WRITE(6,820) LL4
      IF(IMPX.GT.2) THEN
         WRITE (6,830) MUMAX,IIMAX
         WRITE (6,840)
         DO 710 K=1,LXH*ISPLH**2
         WRITE (6,850) K,(IZGLOB(K,I),I=1,3)
  710    CONTINUE
         WRITE (6,860)
         DO 720 K=1,LT4
         WRITE (6,870) K,(KN(K,I),I=1,4+2*NELEM)
         WRITE (6,880) 'X',(KN(K,I),I=4+2*NELEM+1,4+4*NELEM)
         WRITE (6,880) 'Y',(KN(K,I),I=4+4*NELEM+1,4+6*NELEM)
  720    CONTINUE
         WRITE (6,890)
         DO 730 K=1,LXH*ISPLH**2
         WRITE (6,900) K,(QFR(K,I),I=1,6)
  730    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IZGLOB,IJP,IP)
      RETURN
*
  800 FORMAT(1X,A3/14(2X,I6))
  810 FORMAT(1X,A3/7(2X,E12.5))
  820 FORMAT(31H NUMBER OF UNKNOWNS PER GROUP =,I6)
  830 FORMAT(/41H BIVSFH: MAXIMUM BANDWIDTH FOR MATRICES =,I6/9X,
     1 51HNUMBER OF TERMS IN THE COMPRESSED SYSTEM MATRICES =,I10)
  840 FORMAT(/22H NUMBERING OF HEXAGONS/1X,21(1H-)//8H ELEMENT,4X,
     1 24H W ----- X ----- Y -----)
  850 FORMAT(1X,I6,5X,3I8)
  860 FORMAT(/22H NUMBERING OF UNKNOWNS/1X,21(1H-)//8H ELEMENT,5X,
     1 27H---> W ---> X ---> Y ---> W,4X,8HCURRENTS,89(1H.))
  870 FORMAT(1X,I6,5X,4I7,4X,1HW,12I8:/(45X,12I8))
  880 FORMAT(44X,A1,12I8:/(45X,12I8))
  890 FORMAT(/8H ELEMENT,3X,23HVOID BOUNDARY CONDITION/15X,7(1H-),
     1 3H W ,7(1H-),3X,7(1H-),3H X ,7(1H-),3X,7(1H-),3H Y ,7(1H-))
  900 FORMAT(1X,I6,5X,1P,10E10.1/(12X,1P,10E10.1))
      END
