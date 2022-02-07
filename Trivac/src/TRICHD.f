*DECK TRICHD
      SUBROUTINE TRICHD(IMPX,LX,LY,LZ,CYLIND,IELEM,L4,LL4F,LL4X,
     1 LL4Y,LL4Z,MAT,VOL,XX,YY,ZZ,DD,KN,V,MUX,MUY,MUZ,IPBBX,IPBBY,IPBBZ,
     2 BBX,BBY,BBZ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Thomas-Raviart (dual) finite element unknown numbering for ADI
* solution in a 3D domain. Compute the storage info for ADI matrices
* in compressed diagonal storage mode. Compute the ADI permutation
* vectors. Compute the group-independent XB, YB and ZB matrices.
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
* IMPX    print parameter.
* LX      number of elements along the X axis.
* LY      number of elements along the Y axis.
* LZ      number of elements along the Z axis.
* CYLIND  cylindrical geometry flag (set with CYLIND=.true.).
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic).
* L4      total number of unknown (variational coefficients) per
*         energy group (order of system matrices).
* LL4F    exact number of flux unknowns.
* LL4X    exact number of X-directed current unknowns.
* LL4Y    exact number of Y-directed current unknowns.
* LL4Z    exact number of Z-directed current unknowns.
* MAT     mixture index assigned to each element.
* VOL     volume of each element.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* ZZ      Z-directed mesh spacings.
* DD      used with cylindrical geometry.
* KN      element-ordered unknown list.
* V       finite element unit matrix.
*
*Parameters: output
* MUX     X-directed compressed diagonal mode indices.
* MUY     Y-directed compressed diagonal mode indices.
* MUZ     Z-directed compressed diagonal mode indices.
* IPBBX   X-directed perdue storage indices.
* IPBBY   Y-directed perdue storage indices.
* IPBBZ   Z-directed perdue storage indices.
* BBX     X-directed flux-current matrices.
* BBY     Y-directed flux-current matrices.
* BBZ     Z-directed flux-current matrices.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      LOGICAL CYLIND
      INTEGER IMPX,LX,LY,LZ,IELEM,L4,LL4F,LL4X,LL4Y,LL4Z,
     1 MAT(LX*LY*LZ),KN(LX*LY*LZ*(1+6*IELEM**2)),MUX(L4),MUY(L4),
     2 MUZ(L4),IPBBX(2*IELEM,LL4X),IPBBY(2*IELEM,LL4Y),
     3 IPBBZ(2*IELEM,LL4Z)
      REAL VOL(LX*LY*LZ),XX(LX*LY*LZ),YY(LX*LY*LZ),ZZ(LX*LY*LZ),
     1 DD(LX*LY*LZ),V(IELEM+1,IELEM),BBX(2*IELEM,LL4X),
     2 BBY(2*IELEM,LL4Y),BBZ(2*IELEM,LL4Z)
*
      IF(IELEM.GT.4) CALL XABORT('TRICHD: 1 .LE. IELEM .LE. 3.')
      IF(L4.NE.LL4F+LL4X+LL4Y+LL4Z) CALL XABORT('TRICHD: INVALID L4.')
*----
*  COMPUTE THE X-ORIENTED SYSTEM BANDWIDTH VECTOR
*----
      CALL XDISET(MUX,L4,1)
      CALL XDISET(IPBBX,2*IELEM*LL4X,0)
      NUM1=0
      DO 20 KEL=1,LX*LY*LZ
      IF(MAT(KEL).EQ.0) GO TO 20
      DO 12 K3=0,IELEM-1
      DO 11 K2=0,IELEM-1
      KN1=KN(NUM1+2+K3*IELEM+K2)
      KN2=KN(NUM1+2+IELEM**2+K3*IELEM+K2)
      INX1=ABS(KN1)-LL4F
      INX2=ABS(KN2)-LL4F
      IF((KN1.NE.0).AND.(KN2.NE.0)) THEN
         MUX(INX2)=MAX(MUX(INX2),INX2-INX1+1)
         MUX(INX1)=MAX(MUX(INX1),INX1-INX2+1)
      ENDIF
      DO 10 K1=0,IELEM-1
      JND1=KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
      IF(KN1.NE.0) CALL TRINDX(JND1,IPBBX(1,INX1),2*IELEM)
      IF(KN2.NE.0) CALL TRINDX(JND1,IPBBX(1,INX2),2*IELEM)
   10 CONTINUE
   11 CONTINUE
   12 CONTINUE
      NUM1=NUM1+1+6*IELEM**2
   20 CONTINUE
*----
*  COMPUTE THE Y-ORIENTED SYSTEM BANDWIDTH VECTOR
*----
      CALL XDISET(MUY,L4,1)
      CALL XDISET(IPBBY,2*IELEM*LL4Y,0)
      NUM1=0
      DO 50 KEL=1,LX*LY*LZ
      IF(MAT(KEL).EQ.0) GO TO 50
      DO 42 K3=0,IELEM-1
      DO 41 K1=0,IELEM-1
      KN1=KN(NUM1+2+2*IELEM**2+K3*IELEM+K1)
      KN2=KN(NUM1+2+3*IELEM**2+K3*IELEM+K1)
      INY1=ABS(KN1)-LL4F-LL4X
      INY2=ABS(KN2)-LL4F-LL4X
      IF((KN1.NE.0).AND.(KN2.NE.0)) THEN
         MUY(INY2)=MAX(MUY(INY2),INY2-INY1+1)
         MUY(INY1)=MAX(MUY(INY1),INY1-INY2+1)
      ENDIF
      DO 40 K2=0,IELEM-1
      JND1=KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
      IF(KN1.NE.0) CALL TRINDX(JND1,IPBBY(1,INY1),2*IELEM)
      IF(KN2.NE.0) CALL TRINDX(JND1,IPBBY(1,INY2),2*IELEM)
   40 CONTINUE
   41 CONTINUE
   42 CONTINUE
      NUM1=NUM1+1+6*IELEM**2
   50 CONTINUE
*----
*  COMPUTE THE Z-ORIENTED SYSTEM BANDWIDTH VECTOR
*----
      CALL XDISET(MUZ,L4,1)
      CALL XDISET(IPBBZ,2*IELEM*LL4Z,0)
      NUM1=0
      DO 70 KEL=1,LX*LY*LZ
      IF(MAT(KEL).EQ.0) GO TO 70
      DO 62 K2=0,IELEM-1
      DO 61 K1=0,IELEM-1
      KN1=KN(NUM1+2+4*IELEM**2+K2*IELEM+K1)
      KN2=KN(NUM1+2+5*IELEM**2+K2*IELEM+K1)
      INZ1=ABS(KN1)-LL4F-LL4X-LL4Y
      INZ2=ABS(KN2)-LL4F-LL4X-LL4Y
      IF((KN1.NE.0).AND.(KN2.NE.0)) THEN
         MUZ(INZ2)=MAX(MUZ(INZ2),INZ2-INZ1+1)
         MUZ(INZ1)=MAX(MUZ(INZ1),INZ1-INZ2+1)
      ENDIF
      DO 60 K3=0,IELEM-1
      JND1=KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
      IF(KN1.NE.0) CALL TRINDX(JND1,IPBBZ(1,INZ1),2*IELEM)
      IF(KN2.NE.0) CALL TRINDX(JND1,IPBBZ(1,INZ2),2*IELEM)
   60 CONTINUE
   61 CONTINUE
   62 CONTINUE
      NUM1=NUM1+1+6*IELEM**2
   70 CONTINUE
*
      MUXMAX=0
      IIMAXX=0
      DO 80 I=1,LL4X
      MUXMAX=MAX(MUXMAX,MUX(I))
      IIMAXX=IIMAXX+MUX(I)
      MUX(I)=IIMAXX
   80 CONTINUE
*
      MUYMAX=0
      IIMAXY=0
      DO 90 I=1,LL4Y
      MUYMAX=MAX(MUYMAX,MUY(I))
      IIMAXY=IIMAXY+MUY(I)
      MUY(I)=IIMAXY
   90 CONTINUE
*
      MUZMAX=0
      IIMAXZ=0
      DO 100 I=1,LL4Z
      MUZMAX=MAX(MUZMAX,MUZ(I))
      IIMAXZ=IIMAXZ+MUZ(I)
      MUZ(I)=IIMAXZ
  100 CONTINUE
      IF(IMPX.GT.0) THEN
         WRITE (6,600) MUXMAX,MUYMAX,MUZMAX
         WRITE (6,610) IIMAXX,IIMAXY,IIMAXZ
      ENDIF
*----
*  COMPUTE THE FLUX-CURRENT COUPLING MATRICES XB, YB AND ZB.
*----
      CALL XDRSET(BBX,2*IELEM*LL4X,0.0)
      CALL XDRSET(BBY,2*IELEM*LL4Y,0.0)
      CALL XDRSET(BBZ,2*IELEM*LL4Z,0.0)
      NUM1=0
      DO 270 IE=1,LX*LY*LZ
      L=MAT(IE)
      IF(L.EQ.0) GO TO 270
      VOL0=VOL(IE)
      IF(VOL0.EQ.0.0) GO TO 260
      DX=XX(IE)
      DY=YY(IE)
      DZ=ZZ(IE)
      IF(CYLIND) THEN
         DIN=1.0-0.5*DX/DD(IE)
         DOT=1.0+0.5*DX/DD(IE)
      ELSE
         DIN=1.0
         DOT=1.0
      ENDIF
*
      DO 152 K3=0,IELEM-1
      DO 151 K2=0,IELEM-1
      INX1=ABS(KN(NUM1+2+K3*IELEM+K2))-LL4F
      INX2=ABS(KN(NUM1+2+IELEM**2+K3*IELEM+K2))-LL4F
      DO 150 K1=0,IELEM-1
      JND1=KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
      IF(KN(NUM1+2+K3*IELEM+K2).NE.0) THEN
         KK=0
         DO 110 I=1,2*IELEM
         IF(IPBBX(I,INX1).EQ.JND1) THEN
            KK=I
            GO TO 120
         ENDIF
  110    CONTINUE
         CALL XABORT('TRICHD: BUG1.')
  120    SG=REAL(SIGN(1,KN(NUM1+2+K3*IELEM+K2)))
         BBX(KK,INX1)=BBX(KK,INX1)+SG*(VOL0/DX)*DIN*V(1,K1+1)
      ENDIF
      IF(KN(NUM1+2+IELEM**2+K3*IELEM+K2).NE.0) THEN
         KK=0
         DO 130 I=1,2*IELEM
         IF(IPBBX(I,INX2).EQ.JND1) THEN
            KK=I
            GO TO 140
         ENDIF
  130    CONTINUE
         CALL XABORT('TRICHD: BUG2.')
  140    SG=REAL(SIGN(1,KN(NUM1+2+IELEM**2+K3*IELEM+K2)))
         BBX(KK,INX2)=BBX(KK,INX2)+SG*(VOL0/DX)*DOT*V(IELEM+1,K1+1)
      ENDIF
  150 CONTINUE
  151 CONTINUE
  152 CONTINUE
*
      DO 202 K3=0,IELEM-1
      DO 201 K1=0,IELEM-1
      INY1=ABS(KN(NUM1+2+2*IELEM**2+K3*IELEM+K1))-LL4F-LL4X
      INY2=ABS(KN(NUM1+2+3*IELEM**2+K3*IELEM+K1))-LL4F-LL4X
      DO 200 K2=0,IELEM-1
      JND1=KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
      IF(KN(NUM1+2+2*IELEM**2+K3*IELEM+K1).NE.0) THEN
         KK=0
         DO 160 I=1,2*IELEM
         IF(IPBBY(I,INY1).EQ.JND1) THEN
            KK=I
            GO TO 170
         ENDIF
  160    CONTINUE
         CALL XABORT('TRICHD: BUG3.')
  170    SG=REAL(SIGN(1,KN(NUM1+2+2*IELEM**2+K3*IELEM+K1)))
         BBY(KK,INY1)=BBY(KK,INY1)+SG*(VOL0/DY)*V(1,K2+1)
      ENDIF
      IF(KN(NUM1+2+3*IELEM**2+K3*IELEM+K1).NE.0) THEN
         KK=0
         DO 180 I=1,2*IELEM
         IF(IPBBY(I,INY2).EQ.JND1) THEN
            KK=I
            GO TO 190
         ENDIF
  180    CONTINUE
         CALL XABORT('TRICHD: BUG4.')
  190    SG=REAL(SIGN(1,KN(NUM1+2+3*IELEM**2+K3*IELEM+K1)))
         BBY(KK,INY2)=BBY(KK,INY2)+SG*(VOL0/DY)*V(IELEM+1,K2+1)
      ENDIF
  200 CONTINUE
  201 CONTINUE
  202 CONTINUE
*
      DO 252 K2=0,IELEM-1
      DO 251 K1=0,IELEM-1
      INZ1=ABS(KN(NUM1+2+4*IELEM**2+K2*IELEM+K1))-LL4F-LL4X-LL4Y
      INZ2=ABS(KN(NUM1+2+5*IELEM**2+K2*IELEM+K1))-LL4F-LL4X-LL4Y
      DO 250 K3=0,IELEM-1
      JND1=KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
      IF(KN(NUM1+2+4*IELEM**2+K2*IELEM+K1).NE.0) THEN
         KK=0
         DO 210 I=1,2*IELEM
         IF(IPBBZ(I,INZ1).EQ.JND1) THEN
            KK=I
            GO TO 220
         ENDIF
  210    CONTINUE
         CALL XABORT('TRICHD: BUG5.')
  220    SG=REAL(SIGN(1,KN(NUM1+2+4*IELEM**2+K2*IELEM+K1)))
         BBZ(KK,INZ1)=BBZ(KK,INZ1)+SG*(VOL0/DZ)*V(1,K3+1)
      ENDIF
      IF(KN(NUM1+2+5*IELEM**2+K2*IELEM+K1).NE.0) THEN
         KK=0
         DO 230 I=1,2*IELEM
         IF(IPBBZ(I,INZ2).EQ.JND1) THEN
            KK=I
            GO TO 240
         ENDIF
  230    CONTINUE
         CALL XABORT('TRICHD: BUG6.')
  240    SG=REAL(SIGN(1,KN(NUM1+2+5*IELEM**2+K2*IELEM+K1)))
         BBZ(KK,INZ2)=BBZ(KK,INZ2)+SG*(VOL0/DZ)*V(IELEM+1,K3+1)
      ENDIF
  250 CONTINUE
  251 CONTINUE
  252 CONTINUE
  260 NUM1=NUM1+1+6*IELEM**2
  270 CONTINUE
      RETURN
*
  600 FORMAT(/52H TRICHD: MAXIMUM BANDWIDTH FOR X-ORIENTED MATRICES =,
     1 I4/27X,25HFOR Y-ORIENTED MATRICES =,I4/27X,16HFOR Z-ORIENTED M,
     2 9HATRICES =,I4)
  610 FORMAT(/40H TRICHD: LENGTH OF X-ORIENTED MATRICES =,I10/16X,
     1 24HOF Y-ORIENTED MATRICES =,I10/16X,24HOF Z-ORIENTED MATRICES =,
     2 I10)
      END
