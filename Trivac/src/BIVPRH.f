*DECK BIVPRH
      SUBROUTINE BIVPRH (MAXEV,MAXKN,IMPX,ISPLH,LX,IHEX,NCODE,ICODE,
     1 ZCODE,MAT,SIDE,LL4,NELEM,VOL,KN,QFR,IQFR,BFR,MUW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a mesh corner finite difference or linear
* Lagrangian finite element discretization of a 2-D hexagonal geometry.
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
* MAXKN   dimension for arrays KN, QFR and BFR.
* IMPX    print parameter.
* ISPLH   type of hexagonal mesh-splitting:
*         =1: hexagonal elements; >1: 6*(ISPLH-1)**2 triangular elements
*         per hexagon.
* LX      number of hexagons.
* IHEX    type of hexagonal boundary condition.
* NCODE   type of boundary condition applied on each side
*         (i=1: X-  i=2: X+  i=3: Y-  i=4: Y+):
*         NCODE(I)=1: VOID;   NCODE(I)=2: REFL;   NCODE(I)=5: SYME;
*         NCODE(I)=7: ZERO.
* ICODE   physical albedo index on each side of the domain.
* ZCODE   ZCODE(I) is the albedo corresponding to boundary condition
*         'VOID' on each side (ZCODE(I)=0.0 by default).
* MAT     mixture index assigned to each hexagon.
* SIDE    side of the hexagon.
*
*Parameters: output
* LL4     order of system matrices.
* NELEM   number of finite elements (hexagons or triangles) excluding
*         the virtual elements.
* VOL     volume of each hexagon.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* IQFR    element-ordered physical albedo indices.
* BFR     element-ordered surface fractions.
* MUW     compressed storage mode indices.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXEV,MAXKN,IMPX,ISPLH,LX,IHEX,NCODE(4),ICODE(4),MAT(LX),
     1 LL4,NELEM,KN(MAXKN),IQFR(MAXKN),MUW(LL4)
      REAL ZCODE(4),SIDE,VOL(LX),QFR(MAXKN),BFR(MAXKN)
*----
*  LOCAL VARIABLES
*----
      INTEGER ISR(6,2),JCR(6),KK(6),ISRH(6,2),ISRT(3,2),ISRT2(3,2)
      CHARACTER HSMG*131
      LOGICAL LOG
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGAR,KN1
      ALB(X)=0.5*(1.0-X)/(1.0+X)
      DATA ISRH/2,1,4,5,6,3,1,4,5,6,3,2/
      DATA ISRT/1,2,3,2,3,1/
      DATA ISRT2/3,1,2,2,3,1/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IGAR(MAXEV),KN1(MAXKN))
*
      IF(IMPX.GT.0) WRITE(6,500)
      IF(ISPLH.EQ.1) THEN
         NSURF=6
         MAXNH=LX
      ELSE
         NSURF=3
         MAXNH=(6*(ISPLH-1)**2)*LX
      ENDIF
      CALL BIVSBH (MAXNH,MAXKN,IMPX,ISPLH,LX,SIDE,NELEM,IHEX,NCODE,
     1 MAT,VOL,KN1,QFR)
      IF(NELEM*NSURF.GT.MAXKN) THEN
         WRITE(HSMG,'(28HBIVPRH: INSUFFICIENT MAXKN (,I7,10H). SHOULD ,
     1   15HBE INCREASED TO,I7,1H.)') MAXKN,NELEM*NSURF
         CALL XABORT(HSMG)
      ENDIF
*----
*  PRODUCE STANDARD MESH CORNER FINITE DIFFERENCE NUMBERING.
*----
      DO 10 I=1,NELEM*(NSURF+1)
      KN(I)=-98
   10 CONTINUE
      DO 20 IC=1,NSURF
      IF(ISPLH.EQ.1) THEN
         JCR(IC)=IC+3-(IC/4)*6
         ISR(IC,1)=ISRH(IC,1)
         ISR(IC,2)=ISRH(IC,2)
      ELSE
         IF(IC.EQ.1) JCR(1)=1
         IF(IC.EQ.2) JCR(2)=3
         IF(IC.EQ.3) JCR(3)=2
         ISR(IC,1)=ISRT(IC,1)
         ISR(IC,2)=ISRT(IC,2)
      ENDIF
   20 CONTINUE
*----
*  SET ZERO BOUNDARY CONDITIONS
*----
      NUM1=0
      DO 30 KX=1,NELEM
      DO 26 IC=1,NSURF
      KY=ABS(KN1(NUM1+IC))
      DO 25 I1=1,2
      IF((KY.GT.NELEM).AND.(NCODE(1).EQ.7)) KN(NUM1+ISR(IC,I1))=-99
   25 CONTINUE
   26 CONTINUE
      NUM1=NUM1+NSURF+1
   30 CONTINUE
*
      SURFTOT=0.0
      LL4=0
      NUM1=0
      DO 50 KX=1,NELEM
      DO 40 IC=1,NSURF
      IF(NSURF.EQ.6) THEN
         BFR(NUM1+IC)=QFR(NUM1+IC)*QFR(NUM1+7)/(1.5*SQRT(3.0)*SIDE)
         IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(NUM1+IC)=ALB(ZCODE(1))*QFR(NUM1+IC)*QFR(NUM1+7)/
     1      (1.5*SQRT(3.0)*SIDE)
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(NUM1+IC)=QFR(NUM1+IC)*QFR(NUM1+7)/(1.5*SQRT(3.0)*SIDE)
         ELSE
            QFR(NUM1+IC)=0.0
         ENDIF
      ELSE
         AA=SIDE/REAL(ISPLH-1)
         BFR(NUM1+IC)=QFR(NUM1+IC)*QFR(NUM1+4)/(0.25*SQRT(3.0)*AA)
         IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(NUM1+IC)=ALB(ZCODE(1))*QFR(NUM1+IC)*QFR(NUM1+4)/
     1      (0.25*SQRT(3.0)*AA)
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(NUM1+IC)=QFR(NUM1+IC)*QFR(NUM1+4)/(0.25*SQRT(3.0)*AA)
         ELSE
            QFR(NUM1+IC)=0.0
         ENDIF
      ENDIF
      IQFR(NUM1+IC)=ICODE(1)
      SURFTOT=SURFTOT+BFR(NUM1+IC)
      KY0=KN1(NUM1+IC)
      IF((KY0.LT.0).AND.(IHEX.NE.5).AND.(IHEX.NE.6)) GO TO 40
      KY=ABS(KY0)
      DO 35 I1=1,2
      IND=ISR(IC,I1)
      IF((KY.GT.NELEM).OR.(KY0.LT.0)) THEN
         IF(KN(NUM1+IND).EQ.-98) THEN
            LL4=LL4+1
            KN(NUM1+IND)=LL4
         ENDIF
      ELSE
         JND=ISR(JCR(IC),I1+1-(I1/2)*2)
         IOF2=(KY-1)*(NSURF+1)+JND
         IF(IOF2.GT.MAXKN) CALL XABORT('BIVPRH: ALGORITHM FAILURE 2.')
         LOG=.FALSE.
         IF(KN(IOF2).EQ.-99) THEN
            KN(NUM1+IND)=-99
         ELSE IF(KN(NUM1+IND).EQ.-98) THEN
            LL4=LL4+1
            KN(NUM1+IND)=LL4
            IF(KY.NE.KX) KN(IOF2)=LL4
            IF((KY.NE.KX).AND.(ISPLH.GT.1)) LOG=.TRUE.
         ELSE IF(KN(NUM1+IND).EQ.-99) THEN
            GO TO 35
         ELSE IF((KN(IOF2).EQ.-98).AND.(KY.NE.KX)) THEN
            KN(IOF2)=KN(NUM1+IND)
            IF((KY.NE.KX).AND.(ISPLH.GT.1)) LOG=.TRUE.
         ELSE IF((KN(NUM1+IND).NE.KN(IOF2)).AND.(KY.NE.KX)) THEN
            CALL XABORT('BIVPRH: ALGORITHM FAILURE 3.')
         ELSE IF((KY.NE.KX).AND.(ISPLH.GT.1)) THEN
            LOG=.TRUE.
         ENDIF
         IF(LOG) THEN
            KND=0
            IF(JND.EQ.1) KND=2
            IF(JND.EQ.2) KND=1
            IF(JND.EQ.3) KND=3
            KZ=KN1((KY-1)*4+ISRT2(JCR(IC),I1+1-(I1/2)*2))
            IF((KZ.GT.0).AND.(KZ.LE.NELEM).AND.(KZ.NE.KY)) THEN
               IF(KN((KZ-1)*4+KND).EQ.-99) THEN
                  KN(NUM1+IND)=-99
               ELSE
                  KN((KZ-1)*4+KND)=KN(NUM1+IND)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
   35 CONTINUE
   40 CONTINUE
      KN(NUM1+NSURF+1)=KN1(KX*(NSURF+1))
      NUM1=NUM1+NSURF+1
   50 CONTINUE
*----
*  COMPUTE THE SURFACE FRACTIONS.
*----
      IF(SURFTOT.GT.0.0) THEN
         DO 55 I=1,NUM1
         BFR(I)=BFR(I)/SURFTOT
   55    CONTINUE
      ENDIF
*
      NUM1=0
      DO 150 KX=1,NELEM
      DO 60 IC=1,NSURF
      KK(IC)=KN1(NUM1+IC)
   60 CONTINUE
      IF(ISPLH.EQ.1) THEN
         IF((KX.EQ.1).AND.((IHEX.EQ.1).OR.(IHEX.EQ.10))) THEN
            DO 70 I=1,6
            KN(I)=KN(2)
   70       CONTINUE
         ELSE IF((KX.EQ.1).AND.((IHEX.EQ.2).OR.(IHEX.EQ.11))) THEN
            DO 80 I=1,6
            KN(I)=KN(2)
   80       CONTINUE
         ELSE IF((KX.EQ.1).AND.(IHEX.EQ.3)) THEN
            KN(3)=KN(1)
            KN(4)=KN(2)
            KN(5)=KN(1)
            KN(6)=KN(2)
         ELSE IF((KX.EQ.1).AND.(IHEX.EQ.4)) THEN
            KN(3)=KN(1)
            KN(4)=KN(1)
            KN(5)=KN(2)
            KN(6)=KN(1)
         ELSE IF((KX.EQ.1).AND.(IHEX.EQ.5)) THEN
            KN(3)=KN(1)
            KN(4)=KN(2)
            KN(5)=KN(1)
            KN(6)=KN(2)
         ELSE IF((KX.EQ.1).AND.(IHEX.EQ.6)) THEN
            KN(4)=KN(3)
            KN(5)=KN(2)
            KN(6)=KN(1)
         ELSE IF((KK(1).EQ.-KK(5)).AND.(KK(2).EQ.-KK(4)).AND.
     1   (KK(3).EQ.-KK(5)).AND.(KK(6).EQ.-KK(4))) THEN
            DO 90 I=1,6
            KN(NUM1+I)=KN(NUM1+6)
   90       CONTINUE
         ELSE IF((KK(1).EQ.-KK(3)).AND.(KK(4).EQ.-KK(2)).AND.
     1   (KK(5).EQ.-KK(3)).AND.(KK(6).EQ.-KK(2))) THEN
            DO 100 I=1,6
            KN(NUM1+I)=KN(NUM1+4)
  100       CONTINUE
         ELSE IF((KK(1).EQ.-KK(6)).AND.(KK(2).EQ.-KK(5)).AND.
     1   (KK(3).EQ.-KK(4))) THEN
            KN(NUM1+3)=KN(NUM1+1)
            KN(NUM1+6)=KN(NUM1+4)
         ELSE IF((KK(5).EQ.-KK(4)).AND.(KK(6).EQ.-KK(3)).AND.
     1   (KK(1).EQ.-KK(2))) THEN
            KN(NUM1+4)=KN(NUM1+2)
            KN(NUM1+5)=KN(NUM1+3)
         ELSE IF((KK(1).EQ.-KK(3)).AND.(KK(4).EQ.-KK(3)).AND.
     1   (KK(5).EQ.-KK(2)).AND.(KK(6).EQ.-KK(3))) THEN
            KN(NUM1+1)=KN(NUM1+4)
            KN(NUM1+2)=KN(NUM1+5)
            KN(NUM1+3)=KN(NUM1+4)
            KN(NUM1+6)=KN(NUM1+4)
         ELSE IF((KK(2).EQ.-KK(6)).AND.(KK(3).EQ.-KK(5)).AND.
     1   (KK(2).LT.0)) THEN
            KN(NUM1+1)=KN(NUM1+2)
            KN(NUM1+4)=KN(NUM1+3)
            KN(NUM1+5)=KN(NUM1+6)
         ELSE IF((KK(1).EQ.-KK(3)).AND.(KK(6).EQ.-KK(4)).AND.
     1   (KK(1).LT.0)) THEN
            KN(NUM1+1)=KN(NUM1+4)
            KN(NUM1+2)=KN(NUM1+5)
            KN(NUM1+3)=KN(NUM1+6)
         ELSE IF((KK(4).EQ.-KK(2)).AND.(KK(5).EQ.-KK(1)).AND.
     1   (KK(4).LT.0)) THEN
            KN(NUM1+3)=KN(NUM1+2)
            KN(NUM1+5)=KN(NUM1+4)
            KN(NUM1+6)=KN(NUM1+1)
         ELSE IF((KK(3).EQ.-KK(1)).AND.(KK(4).EQ.-KK(6)).AND.
     1   (KK(3).LT.0)) THEN
            KN(NUM1+4)=KN(NUM1+1)
            KN(NUM1+5)=KN(NUM1+2)
            KN(NUM1+6)=KN(NUM1+3)
         ENDIF
         DO 120 IC=1,NSURF
         IF((KK(IC).LT.0).AND.((IHEX.EQ.5).OR.(IHEX.EQ.6)).AND.
     1   (KX.GT.1)) THEN
            IS=0
            DO 110 I=1,6
            IF(-KN1((-KK(IC)-1)*7+I).EQ.KX) IS=I
  110       CONTINUE
            IF(IS.EQ.0) CALL XABORT('BIVPRH: ALGORITHM FAILURE 4.')
            KN((-KK(IC)-1)*7+ISRH(IS,2))=KN(NUM1+ISRH(IC,1))
            KN((-KK(IC)-1)*7+ISRH(IS,1))=KN(NUM1+ISRH(IC,2))
         ENDIF
  120    CONTINUE
      ELSE
         IF((IHEX.NE.5).AND.(IHEX.NE.6)) THEN
            IF((KK(1).EQ.-KK(2)).AND.(KK(1).LT.0)) THEN
               KN(NUM1+1)=KN(NUM1+3)
            ELSE IF((KK(1).EQ.-KK(3)).AND.(KK(1).LT.0)) THEN
               KN(NUM1+2)=KN(NUM1+3)
            ELSE IF((KK(2).EQ.-KK(3)).AND.(KK(2).LT.0)) THEN
               KN(NUM1+2)=KN(NUM1+1)
            ELSE IF((KK(2).EQ.-KK(1)).AND.(KK(2).LT.0)) THEN
               KN(NUM1+3)=KN(NUM1+1)
            ELSE IF((KK(3).EQ.-KK(1)).AND.(KK(3).LT.0)) THEN
               KN(NUM1+3)=KN(NUM1+2)
            ELSE IF((KK(3).EQ.-KK(2)).AND.(KK(3).LT.0)) THEN
               KN(NUM1+1)=KN(NUM1+2)
            ENDIF
         ENDIF
         DO 140 IC=1,NSURF
         IF((KK(IC).LT.0).AND.((IHEX.EQ.5).OR.(IHEX.EQ.6))) THEN
            IS=0
            DO 130 I=1,3
            IF(-KN1((-KK(IC)-1)*4+I).EQ.KX) IS=I
  130       CONTINUE
            IF(IS.EQ.0) CALL XABORT('BIVPRH: ALGORITHM FAILURE 5.')
            KY=-KK(IC)
            DO 135 I1=1,2
            J1=I1+1-(I1/2)*2
            IND=ISRT(IC,I1)
            JND=ISRT(IS,J1)
            KN((-KK(IC)-1)*4+JND)=KN(NUM1+IND)
            KND=0
            IF(JND.EQ.1) KND=2
            IF(JND.EQ.2) KND=1
            IF(JND.EQ.3) KND=3
            KZ=KN1((-KK(IC)-1)*4+ISRT2(IS,J1))
            IF((KZ.GT.0).AND.(KZ.LE.NELEM).AND.(KZ.NE.KY)) THEN
               KNZ=KN((KZ-1)*4+KND)
               KN((KZ-1)*4+KND)=KN(NUM1+IND)
               IF(IHEX.EQ.6) THEN
                  DO 132 L=1,NELEM
                  DO 131 LC=1,NSURF
                  IF(KN((L-1)*4+LC).EQ.KNZ) KN((L-1)*4+LC)=KN(NUM1+IND)
  131             CONTINUE
  132             CONTINUE
               ENDIF
            ENDIF
  135       CONTINUE
         ENDIF
  140    CONTINUE
      ENDIF
      NUM1=NUM1+NSURF+1
  150 CONTINUE
      LL5=0
      DO 170 I=1,MAXEV
      IGAR(I)=0
  170 CONTINUE
      NUM1=0
      DO 190 I=1,NELEM
      DO 180 IC=1,NSURF
      IND=KN(NUM1+IC)
      IF(IND.GT.MAXEV) THEN
         WRITE(HSMG,'(28HBIVPRH: INSUFFICIENT MAXEV (,I7,10H). SHOULD ,
     1   15HBE INCREASED TO,I7,1H.)') MAXEV,IND
         CALL XABORT(HSMG)
      ELSE IF(IND.EQ.-98) THEN
         CALL XABORT('BIVPRH: ALGORITHM FAILURE 6.')
      ELSE IF(IND.EQ.-99) THEN
         KN(NUM1+IC)=0
      ELSE IF(IGAR(IND).EQ.0) THEN
         LL5=LL5+1
         IGAR(IND)=LL5
      ENDIF
  180 CONTINUE
      NUM1=NUM1+NSURF+1
  190 CONTINUE
      NUM1=0
      DO 210 I=1,NELEM
      DO 200 IC=1,NSURF
      IF(KN(NUM1+IC).NE.0) THEN
         IF(IGAR(KN(NUM1+IC)).EQ.0) CALL XABORT('BIVPRH: ALGORITHM FAI'
     1   //'LURE 7.')
         KN(NUM1+IC)=IGAR(KN(NUM1+IC))
      ENDIF
  200 CONTINUE
      NUM1=NUM1+NSURF+1
  210 CONTINUE
      LL4=LL5
      IF(IMPX.GT.0) WRITE(6,570) LL4
      IF(LL4.GT.MAXEV) THEN
         WRITE(HSMG,'(28HBIVPRH: INSUFFICIENT MAXEV (,I7,10H). SHOULD ,
     1   15HBE INCREASED TO,I7,1H.)') MAXEV,LL4
         CALL XABORT(HSMG)
      ENDIF
      IF(LL4.GT.MAXEV) CALL XABORT('BIVPRH: INSUFFICIENT MAXEV.')
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
*----
*  COMPUTE THE SYSTEM MATRIX BANDWIDTH.
*----
      DO 240 I=1,LL4
      MUW(I)=1
  240 CONTINUE
      NUM1=0
      DO 270 K=1,NELEM
      DO 260 I=1,NSURF
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 260
      DO 250 J=1,NSURF
      IND2=KN(NUM1+J)
      IF(IND2.EQ.0) GO TO 250
      MUW(IND1)=MAX0(MUW(IND1),IND1-IND2+1)
  250 CONTINUE
  260 CONTINUE
      NUM1=NUM1+NSURF+1
  270 CONTINUE
      IIMAX=0
      DO 280 I=1,LL4
      IIMAX=IIMAX+MUW(I)
      MUW(I)=IIMAX
  280 CONTINUE
      IF(IMPX.GT.6) WRITE(6,550) 'MUW :',(MUW(I),I=1,LL4)
      IF(IMPX.GT.2) WRITE(6,560) IIMAX
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(KN1,IGAR)
      RETURN
*
  500 FORMAT(//52H BIVPRH: NUMBERING FOR A MESH CORNER FINITE DIFFEREN,
     1 60HCE OR LINEAR LAGRANGIAN FINITE ELEMENT DISCRETIZATION IN HEX,
     2 16HAGONAL GEOMETRY.)
  510 FORMAT(/31H BIVPRH: NUMBERING OF UNKNOWNS./1X,30(1H-)/9X,
     1 7HHEXAGON,3X,8HUNKNOWNS,29X,23HVOID BOUNDARY CONDITION,45X,
     2 6HVOLUME)
  520 FORMAT (1X,2I6,2X,6I6,2X,1P,6E11.2,5X,E13.6)
  530 FORMAT(/31H BIVPRH: NUMBERING OF UNKNOWNS./1X,30(1H-)/9X,
     1 7HHEXAGON,3X,8HUNKNOWNS,11X,23HVOID BOUNDARY CONDITION,12X,
     2 6HVOLUME,13X,16HSURFACE FRACTION)
  540 FORMAT (1X,2I6,2X,3I6,2X,1P,3E11.2,5X,E13.6,5X,3E10.2)
  550 FORMAT(/1X,A5/(1X,20I6))
  560 FORMAT(/52H NUMBER OF TERMS IN THE COMPRESSED SYSTEM MATRICES =,
     > I6)
  570 FORMAT(/39H BIVPRH: NUMBER OF UNKNOWNS PER GROUP =,I6/)
  580 FORMAT (//17H SURFACE FRACTION//8H HEXAGON,5X,3HBFR)
  590 FORMAT (3X,I4,4X,1P,6E10.2)
      END
