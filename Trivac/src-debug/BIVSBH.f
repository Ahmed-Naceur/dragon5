*DECK BIVSBH
      SUBROUTINE BIVSBH (MAXEV,MAXKN,IMPX,ISPLH,LX,SIDE,LL4,IHEX,NCODE,
     1 MAT,VOL,KN,QFR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering of an hexagonal 2-D geometry with or without triangular
* mesh-splitting.
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
* MAXEV   dimension for array IGAR:
*         if ISPLH=1: number of hexagons;
*         if ISPLH>1: (6*(ISPLH-1)**2)*LX where LX is the number of
*         hexagons.
* MAXKN   dimension for arrays KN and QFR.
* IMPX    print parameter.
* ISPLH   type of hexagonal mesh-splitting:
*         =1: no mesh splitting (complete hexagons);
*         =K: 6*(K-1)*(K-1) triangles per hexagon.
* LX      number of hexagons.
* SIDE    side of an hexagon.
* NCODE   type of boundary condition applied on each side
*         (i=1: X-  i=2: X+  i=3: Y-  i=4: Y+):
*         NCODE(I)=1: VOID;   NCODE(I)=2: REFL;   NCODE(I)=5: SYME;
*         NCODE(I)=7: ZERO.
* MAT     mixture index assigned to each hexagon.
* IHEX    type of hexagonal boundary condition.
*
*Parameters: output
* LL4     number of elements after mesh-splitting.
* VOL     volume of each hexagon.
* KN      element-ordered unknown list.
* QFR     element-ordered external surfaces: =1.0 on external surfaces;
*         =0.0 on internal surfaces.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXEV,MAXKN,IMPX,ISPLH,LX,LL4,IHEX,NCODE(6),MAT(LX),
     1 KN(MAXKN)
      REAL SIDE,VOL(LX),QFR(7*LX)
*----
*  LOCAL VARIABLES
*----
      INTEGER KK(6)
      CHARACTER HSMG*131
      LOGICAL LOGSUR
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGAR,KN2
      REAL, DIMENSION(:), ALLOCATABLE :: QFR2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IGAR(MAXEV),KN2(MAXKN),QFR2(MAXKN))
*
      IF(LX.GT.MAXEV) THEN
         WRITE(HSMG,'(30HBIVSBH: 1 INSUFFICIENT MAXEV (,I7,7H). SHOU,
     1   18HLD BE INCREASED TO,I7,1H.)') MAXEV,LX
         CALL XABORT(HSMG)
      ENDIF
      LL4=0
      DO 10 KX=1,LX
      IGAR(KX)=0
      IF(MAT(KX).LE.0) GO TO 10
      LL4=LL4+1
      IGAR(KX)=LL4
   10 CONTINUE
      NSURF=6
      NUM1=0
      DO 30 KX=1,LX
      VOL(KX)=0.0
      IF(MAT(KX).LE.0) GO TO 30
      IF(NUM1+7.GT.MAXKN) THEN
         WRITE(HSMG,'(30HBIVSBH: 1 INSUFFICIENT MAXKN (,I7,2H).)') MAXKN
         CALL XABORT(HSMG)
      ENDIF
      LOGSUR=(NCODE(1).EQ.1).OR.(NCODE(1).EQ.7)
      DO 20 IX=1,6
      N1=NEIGHB(KX,IX,IHEX,LX,POIDS)
      IF(N1.EQ.0) CALL XABORT('BIVSBH: NEIGHB FAILURE.')
      QFR(NUM1+IX)=0.0
      IF(ABS(N1).GT.LX) THEN
         IF(LOGSUR) QFR(NUM1+IX)=1.0
         KN(NUM1+IX)=SIGN(LX+1,N1)
      ELSE IF(MAT(ABS(N1)).LE.0) THEN
         IF(LOGSUR) QFR(NUM1+IX)=1.0
         KN(NUM1+IX)=SIGN(LX+1,N1)
         IF((IHEX.EQ.5).OR.(IHEX.EQ.6)) KN(NUM1+IX)=LX+1
      ELSE
         KN(NUM1+IX)=SIGN(IGAR(ABS(N1)),N1)
      ENDIF
   20 CONTINUE
      KN(NUM1+7)=KX
      VOL(KX)=2.59807587*SIDE*SIDE*POIDS
      QFR(NUM1+7)=VOL(KX)
      NUM1=NUM1+7
   30 CONTINUE
      MAXMAX=LX
      IF(IMPX.GT.4) THEN
         WRITE(6,510) 1
         NUM1=0
         DO 40 I=1,LL4
         WRITE(6,520) I,KN(NUM1+7),(KN(NUM1+J),J=1,6),(QFR(NUM1+J),
     1   J=1,7)
         NUM1=NUM1+7
   40    CONTINUE
      ENDIF
      IF(ISPLH.GE.2) THEN
*        HEXAGON TO TRIANGLE.
         NSURF=3
         IF(LL4*24.GT.MAXKN) THEN
            WRITE(HSMG,'(30HBIVSBH: 2 INSUFFICIENT MAXKN (,I7,7H). SHOU,
     1      18HLD BE INCREASED TO,I7,1H.)') MAXKN,LL4*24
            CALL XABORT(HSMG)
         ENDIF
         NUM1=0
         DO 60 KX=1,LL4
         IOF2=(KX-1)*24
         DO 50 IT=1,6
         KK(IT)=KN(NUM1+IT)
         KN2(IOF2+(IT-1)*4+4)=KN(NUM1+7)
         QFR2(IOF2+(IT-1)*4+1)=QFR(NUM1+IT)
         QFR2(IOF2+(IT-1)*4+2)=0.0
         QFR2(IOF2+(IT-1)*4+3)=0.0
         IF(IT.NE.6) KN2(IOF2+(IT-1)*4+2)=(KX-1)*6+IT+1
         IF(IT.EQ.6) KN2(IOF2+(IT-1)*4+2)=(KX-1)*6+1
         IF(IT.NE.1) KN2(IOF2+(IT-1)*4+3)=(KX-1)*6+IT-1
         IF(IT.EQ.1) KN2(IOF2+(IT-1)*4+3)=(KX-1)*6+6
         QFR2(IOF2+(IT-1)*4+4)=QFR(NUM1+7)/6.0
         IF((KK(IT).LT.0).AND.((IHEX.EQ.5).OR.(IHEX.EQ.6)).AND.
     1   (KX.GT.1)) THEN
            IC=0
            DO 45 I=1,6
            IF(-KN((-KK(IT)-1)*7+I).EQ.KX) IC=I
   45       CONTINUE
            IF(IC.EQ.0) CALL XABORT('BIVSBH: ALGORITHM FAILURE 1.')
            KN2(IOF2+(IT-1)*4+1)=-((-KK(IT)-1)*6+IC)
         ELSE IF(KK(IT).LT.0) THEN
            KN2(IOF2+(IT-1)*4+1)=0
         ELSE IF(KK(IT).EQ.KX) THEN
            KN2(IOF2+(IT-1)*4+1)=((KX-1)*6+IT)
         ELSE IF(ABS(KK(IT)).GT.MAXMAX) THEN
            KN2(IOF2+(IT-1)*4+1)=SIGN(LL4*6+1,KK(IT))
         ELSE
            KN2(IOF2+(IT-1)*4+1)=(KK(IT)-1)*6+IT+3-(IT/4)*6
         ENDIF
   50    CONTINUE
*        CHECK SYMMETRIES.
         IF((KX.EQ.1).AND.((IHEX.EQ.1).OR.(IHEX.EQ.10))) THEN
            KN2(2)=-1
            KN2(3)=1
            QFR2(4)=QFR(7)
         ELSE IF((KX.EQ.1).AND.((IHEX.EQ.2).OR.(IHEX.EQ.11))) THEN
            KN2((1-1)*4+2)=-KN2((1-1)*4+3)
            KN2((6-1)*4+3)=-KN2((6-1)*4+2)
            QFR2((1-1)*4+4)=QFR(7)/2.0
            QFR2((6-1)*4+4)=QFR(7)/2.0
         ELSE IF((KX.EQ.1).AND.(IHEX.EQ.3)) THEN
            KN2(2)=1
            KN2(3)=1
            QFR2(4)=QFR(7)
         ELSE IF((KX.EQ.1).AND.(IHEX.EQ.4)) THEN
            KN2((1-1)*4+3)=1
            KN2((2-1)*4+2)=-KN2((2-1)*4+3)
            QFR2((1-1)*4+4)=2.0*QFR(7)/3.0
            QFR2((2-1)*4+4)=QFR(7)/3.0
         ELSE IF((KX.EQ.1).AND.(IHEX.EQ.5)) THEN
            KN2((1-1)*4+3)=-KN2((1-1)*4+2)
            KN2((2-1)*4+2)=-KN2((2-1)*4+3)
            QFR2((1-1)*4+4)=QFR(7)/2.0
            QFR2((2-1)*4+4)=QFR(7)/2.0
         ELSE IF((KX.EQ.1).AND.(IHEX.EQ.6)) THEN
            KN2((2-1)*4+2)=-6
            KN2((6-1)*4+3)=-2
            QFR2((1-1)*4+4)=QFR(7)/3.0
            QFR2((2-1)*4+4)=QFR(7)/3.0
            QFR2((6-1)*4+4)=QFR(7)/3.0
         ELSE IF((KK(1).EQ.-KK(5)).AND.(KK(2).EQ.-KK(4)).AND.
     1   (KK(3).EQ.-KK(5)).AND.(KK(6).EQ.-KK(4))) THEN
            KN2(IOF2+(4-1)*4+3)=-KN2(IOF2+(4-1)*4+2)
            KN2(IOF2+(5-1)*4+2)=-KN2(IOF2+(5-1)*4+3)
            QFR2(IOF2+(4-1)*4+4)=QFR(NUM1+7)/2.0
            QFR2(IOF2+(5-1)*4+4)=QFR(NUM1+7)/2.0
         ELSE IF((KK(1).EQ.-KK(3)).AND.(KK(4).EQ.-KK(2)).AND.
     1   (KK(5).EQ.-KK(3)).AND.(KK(6).EQ.-KK(2))) THEN
            KN2(IOF2+(2-1)*4+3)=-KN2(IOF2+(2-1)*4+2)
            KN2(IOF2+(3-1)*4+2)=-KN2(IOF2+(3-1)*4+3)
            QFR2(IOF2+(2-1)*4+4)=QFR(NUM1+7)/2.0
            QFR2(IOF2+(3-1)*4+4)=QFR(NUM1+7)/2.0
         ELSE IF((KK(1).EQ.-KK(6)).AND.(KK(2).EQ.-KK(5)).AND.
     1   (KK(3).EQ.-KK(4))) THEN
            KN2(IOF2+(1-1)*4+3)=((KX-1)*6+1)
            KN2(IOF2+(3-1)*4+2)=((KX-1)*6+3)
            QFR2(IOF2+(1-1)*4+4)=QFR(NUM1+7)/3.0
            QFR2(IOF2+(2-1)*4+4)=QFR(NUM1+7)/3.0
            QFR2(IOF2+(3-1)*4+4)=QFR(NUM1+7)/3.0
         ELSE IF((KK(5).EQ.-KK(4)).AND.(KK(6).EQ.-KK(3)).AND.
     1   (KK(1).EQ.-KK(2))) THEN
            KN2(IOF2+(5-1)*4+3)=((KX-1)*6+5)
            KN2(IOF2+(1-1)*4+2)=((KX-1)*6+1)
            QFR2(IOF2+(5-1)*4+4)=QFR(NUM1+7)/3.0
            QFR2(IOF2+(6-1)*4+4)=QFR(NUM1+7)/3.0
            QFR2(IOF2+(1-1)*4+4)=QFR(NUM1+7)/3.0
         ELSE IF((KK(1).EQ.-KK(3)).AND.(KK(4).EQ.-KK(3)).AND.
     1   (KK(5).EQ.-KK(2)).AND.(KK(6).EQ.-KK(3))) THEN
            KN2(IOF2+(2-1)*4+3)=-KN2(IOF2+(2-1)*4+2)
            KN2(IOF2+(3-1)*4+2)=((KX-1)*6+3)
            QFR2(IOF2+(2-1)*4+4)=QFR(NUM1+7)/3.0
            QFR2(IOF2+(3-1)*4+4)=2.0*QFR(NUM1+7)/3.0
         ELSE IF((KK(2).EQ.-KK(6)).AND.(KK(3).EQ.-KK(5)).AND.
     1   (KK(2).LT.0)) THEN
            KN2(IOF2+(1-1)*4+2)=-KN2(IOF2+(1-1)*4+3)
            KN2(IOF2+(4-1)*4+3)=-KN2(IOF2+(4-1)*4+2)
            QFR2(IOF2+(5-1)*4+4)=2.0*QFR2(IOF2+(5-1)*4+4)
            QFR2(IOF2+(6-1)*4+4)=2.0*QFR2(IOF2+(6-1)*4+4)
         ELSE IF((KK(1).EQ.-KK(3)).AND.(KK(6).EQ.-KK(4)).AND.
     1   (KK(1).LT.0)) THEN
            KN2(IOF2+(2-1)*4+3)=-KN2(IOF2+(2-1)*4+2)
            KN2(IOF2+(5-1)*4+2)=-KN2(IOF2+(5-1)*4+3)
            QFR2(IOF2+(3-1)*4+4)=2.0*QFR2(IOF2+(3-1)*4+4)
            QFR2(IOF2+(4-1)*4+4)=2.0*QFR2(IOF2+(4-1)*4+4)
         ELSE IF((KK(4).EQ.-KK(2)).AND.(KK(5).EQ.-KK(1)).AND.
     1   (KK(4).LT.0)) THEN
            KN2(IOF2+(3-1)*4+2)=-KN2(IOF2+(3-1)*4+3)
            KN2(IOF2+(6-1)*4+3)=-KN2(IOF2+(6-1)*4+2)
            QFR2(IOF2+(1-1)*4+4)=2.0*QFR2(IOF2+(1-1)*4+4)
            QFR2(IOF2+(2-1)*4+4)=2.0*QFR2(IOF2+(2-1)*4+4)
         ELSE IF((KK(3).EQ.-KK(1)).AND.(KK(4).EQ.-KK(6)).AND.
     1   (KK(3).LT.0)) THEN
            KN2(IOF2+(2-1)*4+2)=-KN2(IOF2+(2-1)*4+3)
            KN2(IOF2+(5-1)*4+3)=-KN2(IOF2+(5-1)*4+2)
            QFR2(IOF2+(1-1)*4+4)=2.0*QFR2(IOF2+(1-1)*4+4)
            QFR2(IOF2+(6-1)*4+4)=2.0*QFR2(IOF2+(6-1)*4+4)
         ENDIF
         NUM1=NUM1+7
   60    CONTINUE
         MAXMAX=LL4*6
         IF(LL4*6.GT.MAXEV) THEN
            WRITE(HSMG,'(30HBIVSBH: 2 INSUFFICIENT MAXEV (,I7,7H). SHOU,
     1      18HLD BE INCREASED TO,I7,1H.)') MAXEV,LL4*6
            CALL XABORT(HSMG)
         ENDIF
         LL5=0
         NUM1=0
         NUM2=0
         DO 85 I=1,LL4*6
         IGAR(I)=0
         IF(KN2(NUM2+1).EQ.0) GO TO 80
         LL5=LL5+1
         IGAR(I)=LL5
         DO 70 J=1,4
         KN(NUM1+J)=KN2(NUM2+J)
         QFR(NUM1+J)=QFR2(NUM2+J)
   70    CONTINUE
         NUM1=NUM1+4
   80    NUM2=NUM2+4
   85    CONTINUE
         NUM1=0
         DO 100 I=1,LL5
         DO 90 K=1,3
         IF(ABS(KN(NUM1+K)).LE.LL4*6) THEN
            IF(IGAR(ABS(KN(NUM1+K))).EQ.0) CALL XABORT('BIVSBH: ALGORIT'
     1      //'HM FAILURE 2.')
            KN(NUM1+K)=SIGN(IGAR(ABS(KN(NUM1+K))),KN(NUM1+K))
         ENDIF
   90    CONTINUE
         NUM1=NUM1+4
  100    CONTINUE
         LL4=LL5
         IF(IMPX.GT.4) THEN
            WRITE(6,530) 2
            NUM1=0
            DO 110 I=1,LL4
            WRITE(6,540) I,KN(NUM1+4),(KN(NUM1+J),J=1,3),(QFR(NUM1+J),
     1      J=1,4)
            NUM1=NUM1+4
  110       CONTINUE
         ENDIF
*
*        TRIANGLE TO TRIANGLE.
         KSPLH=0
         IF(ISPLH.EQ.2) THEN
*           MESH-SPLITTING INTO 6 TRIANGLES.
            KSPLH=2
         ELSE IF(ISPLH.EQ.3) THEN
*           MESH-SPLITTING INTO 24 TRIANGLES.
            KSPLH=3
         ELSE IF(ISPLH.EQ.5) THEN
*           MESH-SPLITTING INTO 96 TRIANGLES.
            KSPLH=4
         ELSE IF(ISPLH.EQ.9) THEN
*           MESH-SPLITTING INTO 384 TRIANGLES.
            KSPLH=5
         ELSE IF(ISPLH.EQ.17) THEN
*           MESH-SPLITTING INTO 1536 TRIANGLES.
            KSPLH=6
         ELSE
            WRITE(HSMG,'(36HBIVSBH: UNABLE TO SPLIT WITH ISPLH =,I5,
     1      38H ISPLH = 1, 2, 3, 5, 9 AND 17 ALLOWED.)') ISPLH
            CALL XABORT(HSMG)
         ENDIF
         DO 230 JSPLH=3,KSPLH
         IF(LL4*16.GT.MAXKN) THEN
            WRITE(HSMG,'(30HBIVSBH: 3 INSUFFICIENT MAXKN (,I7,7H). SHOU,
     1      18HLD BE INCREASED TO,I7,1H.)') MAXKN,LL4*16
            CALL XABORT(HSMG)
         ENDIF
         NUM1=0
         DO 170 KX=1,LL4
         IOF2=(KX-1)*16
         DO 120 IT=1,3
         KK(IT)=KN(NUM1+IT)
  120    CONTINUE
         DO 130 IT=1,4
         KN2(IOF2+(IT-1)*4+4)=KN(NUM1+4)
         QFR2(IOF2+(IT-1)*4+1)=0.0
         QFR2(IOF2+(IT-1)*4+2)=0.0
         QFR2(IOF2+(IT-1)*4+3)=0.0
         QFR2(IOF2+(IT-1)*4+4)=QFR(NUM1+4)/4.0
  130    CONTINUE
         QFR2(IOF2+(1-1)*4+3)=QFR(NUM1+1)
         QFR2(IOF2+(3-1)*4+2)=QFR(NUM1+1)
         QFR2(IOF2+(1-1)*4+1)=QFR(NUM1+2)
         QFR2(IOF2+(4-1)*4+2)=QFR(NUM1+2)
         QFR2(IOF2+(3-1)*4+1)=QFR(NUM1+3)
         QFR2(IOF2+(4-1)*4+3)=QFR(NUM1+3)
         KN2(IOF2+(1-1)*4+1)=(KK(2)-1)*4+3
         KN2(IOF2+(1-1)*4+2)=(KX-1)*4+2
         KN2(IOF2+(1-1)*4+3)=(KK(1)-1)*4+3
         KN2(IOF2+(2-1)*4+1)=(KX-1)*4+4
         KN2(IOF2+(2-1)*4+2)=(KX-1)*4+3
         KN2(IOF2+(2-1)*4+3)=(KX-1)*4+1
         KN2(IOF2+(3-1)*4+1)=(KK(3)-1)*4+1
         KN2(IOF2+(3-1)*4+2)=(KK(1)-1)*4+1
         KN2(IOF2+(3-1)*4+3)=(KX-1)*4+2
         KN2(IOF2+(4-1)*4+1)=(KX-1)*4+2
         KN2(IOF2+(4-1)*4+2)=(KK(2)-1)*4+4
         KN2(IOF2+(4-1)*4+3)=(KK(3)-1)*4+4
         IF(ABS(KK(1)).GT.MAXMAX) THEN
            KN2(IOF2+(1-1)*4+3)=SIGN(LL4*4+1,KK(1))
            KN2(IOF2+(3-1)*4+2)=SIGN(LL4*4+1,KK(1))
         ENDIF
         IF(ABS(KK(2)).GT.MAXMAX) THEN
            KN2(IOF2+(1-1)*4+1)=SIGN(LL4*4+1,KK(2))
            KN2(IOF2+(4-1)*4+2)=SIGN(LL4*4+1,KK(2))
         ENDIF
         IF(ABS(KK(3)).GT.MAXMAX) THEN
            KN2(IOF2+(3-1)*4+1)=SIGN(LL4*4+1,KK(3))
            KN2(IOF2+(4-1)*4+3)=SIGN(LL4*4+1,KK(3))
         ENDIF
         IF((KK(1).LT.0).AND.((IHEX.EQ.5).OR.(IHEX.EQ.6))) THEN
            IC=0
            DO 140 I=1,3
            IF(-KN((-KK(1)-1)*4+I).EQ.KX) IC=I
  140       CONTINUE
            IF(IC.EQ.0) CALL XABORT('BIVSBH: ALGORITHM FAILURE 3.')
            KN2(IOF2+(1-1)*4+3)=-((-KK(1)-1)*4+3)
            KN2(IOF2+(3-1)*4+2)=-((-KK(1)-1)*4+1)
         ELSE IF((KK(2).LT.0).AND.((IHEX.EQ.5).OR.(IHEX.EQ.6))) THEN
            IC=0
            DO 150 I=1,3
            IF(-KN((-KK(2)-1)*4+I).EQ.KX) IC=I
  150       CONTINUE
            IF(IC.EQ.0) CALL XABORT('BIVSBH: ALGORITHM FAILURE 4.')
            KN2(IOF2+(1-1)*4+1)=-((-KK(2)-1)*4+3)
            KN2(IOF2+(4-1)*4+2)=-((-KK(2)-1)*4+4)
         ELSE IF((KK(3).LT.0).AND.((IHEX.EQ.5).OR.(IHEX.EQ.6))) THEN
            IC=0
            DO 160 I=1,3
            IF(-KN((-KK(3)-1)*4+I).EQ.KX) IC=I
  160       CONTINUE
            IF(IC.EQ.0) CALL XABORT('BIVSBH: ALGORITHM FAILURE 5.')
            KN2(IOF2+(3-1)*4+1)=-((-KK(3)-1)*4+1)
            KN2(IOF2+(4-1)*4+3)=-((-KK(3)-1)*4+4)
         ELSE IF((KK(1).EQ.-KK(2)).AND.(KK(1).LT.0)) THEN
            KN2(IOF2+(1-1)*4+3)=-KN2(IOF2+(1-1)*4+1)
            KN2(IOF2+(2-1)*4+2)=-KN2(IOF2+(2-1)*4+1)
            KN2(IOF2+(3-1)*4+1)=0
            QFR2(IOF2+(4-1)*4+4)=2.0*QFR2(IOF2+(4-1)*4+4)
         ELSE IF((KK(1).EQ.-KK(3)).AND.(KK(1).LT.0)) THEN
            KN2(IOF2+(2-1)*4+3)=-KN2(IOF2+(2-1)*4+1)
            KN2(IOF2+(3-1)*4+2)=-KN2(IOF2+(3-1)*4+1)
            KN2(IOF2+(1-1)*4+1)=0
            QFR2(IOF2+(4-1)*4+4)=2.0*QFR2(IOF2+(4-1)*4+4)
         ELSE IF((KK(2).EQ.-KK(3)).AND.(KK(2).LT.0)) THEN
            KN2(IOF2+(2-1)*4+3)=-KN2(IOF2+(2-1)*4+2)
            KN2(IOF2+(4-1)*4+2)=-KN2(IOF2+(4-1)*4+3)
            KN2(IOF2+(1-1)*4+1)=0
            QFR2(IOF2+(3-1)*4+4)=2.0*QFR2(IOF2+(3-1)*4+4)
         ELSE IF((KK(2).EQ.-KK(1)).AND.(KK(2).LT.0)) THEN
            KN2(IOF2+(1-1)*4+1)=-KN2(IOF2+(1-1)*4+3)
            KN2(IOF2+(2-1)*4+1)=-KN2(IOF2+(2-1)*4+2)
            KN2(IOF2+(4-1)*4+1)=0
            QFR2(IOF2+(3-1)*4+4)=2.0*QFR2(IOF2+(3-1)*4+4)
         ELSE IF((KK(3).EQ.-KK(1)).AND.(KK(3).LT.0)) THEN
            KN2(IOF2+(2-1)*4+1)=-KN2(IOF2+(2-1)*4+3)
            KN2(IOF2+(3-1)*4+1)=-KN2(IOF2+(3-1)*4+2)
            KN2(IOF2+(4-1)*4+1)=0
            QFR2(IOF2+(1-1)*4+4)=2.0*QFR2(IOF2+(1-1)*4+4)
         ELSE IF((KK(3).EQ.-KK(2)).AND.(KK(3).LT.0)) THEN
            KN2(IOF2+(2-1)*4+2)=-KN2(IOF2+(2-1)*4+3)
            KN2(IOF2+(4-1)*4+3)=-KN2(IOF2+(4-1)*4+2)
            KN2(IOF2+(3-1)*4+1)=0
            QFR2(IOF2+(1-1)*4+4)=2.0*QFR2(IOF2+(1-1)*4+4)
         ENDIF
         IF(KK(1).EQ.KX) THEN
           IF(KN2(IOF2+(1-1)*4+3).NE.0) KN2(IOF2+(1-1)*4+3)=((KX-1)*4+1)
           IF(KN2(IOF2+(3-1)*4+2).NE.0) KN2(IOF2+(3-1)*4+2)=((KX-1)*4+3)
         ENDIF
         IF(KK(2).EQ.KX) THEN
           IF(KN2(IOF2+(1-1)*4+1).NE.0) KN2(IOF2+(1-1)*4+1)=((KX-1)*4+1)
           IF(KN2(IOF2+(4-1)*4+2).NE.0) KN2(IOF2+(4-1)*4+2)=((KX-1)*4+4)
         ENDIF
         IF(KK(3).EQ.KX) THEN
           IF(KN2(IOF2+(3-1)*4+1).NE.0) KN2(IOF2+(3-1)*4+1)=((KX-1)*4+3)
           IF(KN2(IOF2+(4-1)*4+3).NE.0) KN2(IOF2+(4-1)*4+3)=((KX-1)*4+4)
         ENDIF
         NUM1=NUM1+4
  170    CONTINUE
         MAXMAX=LL4*4
         IF(LL4*4.GT.MAXEV) THEN
            WRITE(HSMG,'(30HBIVSBH: 3 INSUFFICIENT MAXEV (,I7,7H). SHOU,
     1      18HLD BE INCREASED TO,I7,1H.)') MAXEV,LL4*4
            CALL XABORT(HSMG)
         ENDIF
         LL5=0
         NUM1=0
         NUM2=0
         DO 195 I=1,LL4*4
         IGAR(I)=0
         IF(KN2(NUM2+1).EQ.0) GO TO 190
         LL5=LL5+1
         IGAR(I)=LL5
         DO 180 J=1,4
         KN(NUM1+J)=KN2(NUM2+J)
         QFR(NUM1+J)=QFR2(NUM2+J)
  180    CONTINUE
         NUM1=NUM1+4
  190    NUM2=NUM2+4
  195    CONTINUE
         NUM1=0
         DO 210 I=1,LL5
         DO 200 K=1,3
         IF(ABS(KN(NUM1+K)).LE.LL4*4) THEN
            IF(IGAR(ABS(KN(NUM1+K))).EQ.0) CALL XABORT('BIVSBH: ALGORIT'
     1      //'HM FAILURE 6.')
            KN(NUM1+K)=SIGN(IGAR(ABS(KN(NUM1+K))),KN(NUM1+K))
         ENDIF
  200    CONTINUE
         NUM1=NUM1+4
  210    CONTINUE
         LL4=LL5
         IF(IMPX.GT.4) THEN
            WRITE(6,530) JSPLH
            NUM1=0
            DO 220 I=1,LL4
            WRITE(6,540) I,KN(NUM1+4),(KN(NUM1+J),J=1,3),(QFR(NUM1+J),
     1      J=1,4)
            NUM1=NUM1+4
  220       CONTINUE
         ENDIF
  230    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IGAR,KN2,QFR2)
      RETURN
*
  510 FORMAT(/36H BIVSBH: NUMBERING OF UNKNOWNS. STEP,I3,1H./1X,40(1H-)/
     1 9X,7HHEXAGON,3X,9HNEIGHBOUR,27X,17HEXTERNAL SURFACES,22X,
     2 6HVOLUME)
  520 FORMAT (1X,2I6,2X,6I6,2X,6F6.2,5X,1P,E13.6)
  530 FORMAT(/36H BIVSBH: NUMBERING OF UNKNOWNS. STEP,I3,1H./1X,40(1H-)/
     1 9X,7HHEXAGON,3X,9HNEIGHBOUR,9X,17HEXTERNAL SURFACES,11X,
     2 6HVOLUME)
  540 FORMAT (1X,2I6,2X,3I6,2X,3F6.2,12X,1P,E13.6)
      END
