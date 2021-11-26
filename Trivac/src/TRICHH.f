*DECK TRICHH
      SUBROUTINE TRICHH(IMPX,MAXKN,NBLOS,LXH,LZ,IELEM,ISPLH,L4,LL4F,
     1 LL4W,LL4X,LL4Y,LL4Z,SIDE,ZZ,FRZ,IPERT,KN,V,H,MUW,MUX,MUY,MUZ,
     2 IPBBW,IPBBX,IPBBY,IPBBZ,BBW,BBX,BBY,BBZ,CTRAN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Thomas-Raviart-Schneider (dual) finite element unknown numbering for
* ADI solution in a 3D hexagonal domain. Compute the storage info for
* ADI matrices in compressed diagonal storage mode. Compute the ADI
* permutation vectors. Compute the group-independent WB, XB, YB and ZB
* matrices.
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
* NBLOS   number of lozenges per direction in 3D with mesh-splitting.
* LXH     number of hexagons in a plane.
* LZ      number of elements along the Z axis.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic).
* ISPLH   mesh-splitting in 3*ISPLH**2 lozenges per hexagon.
* L4      total number of unknown (variational coefficients) per
*         energy group (order of system matrices).
* LL4F    exact number of flux unknowns.
* LL4W    exact number of W-directed current unknowns.
* LL4X    exact number of X-directed current unknowns.
* LL4Y    exact number of Y-directed current unknowns.
* LL4Z    exact number of Z-directed current unknowns.
* SIDE    side of an hexagon.
* ZZ      Z-directed mesh spacings.
* FRZ     volume fractions for the axial SYME boundary condition.
* IPERT   mixture permutation index.
* KN      ADI permutation indices for the volumes and currents.
* V       nodal coupling matrix matrix.
* H       Piolat (hexagonal) coupling matrix.
*
*Parameters: output
* MUW     W-directed compressed diagonal mode indices.
* MUX     X-directed compressed diagonal mode indices.
* MUY     Y-directed compressed diagonal mode indices.
* MUZ     Z-directed compressed diagonal mode indices.
* IPBBW   W-directed perdue storage indices.
* IPBBX   X-directed perdue storage indices.
* IPBBY   Y-directed perdue storage indices.
* IPBBZ   Z-directed perdue storage indices.
* BBW     W-directed flux-current matrices.
* BBX     X-directed flux-current matrices.
* BBY     Y-directed flux-current matrices.
* BBZ     Z-directed flux-current matrices.
* CTRAN   tranverse coupling Piolat unit matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,MAXKN,NBLOS,LXH,LZ,IELEM,ISPLH,L4,IPERT(NBLOS),
     1 KN(NBLOS,MAXKN/NBLOS),LL4F,LL4W,LL4X,LL4Y,LL4Z,MUW(L4),
     2 MUX(L4),MUY(L4),MUZ(L4),IPBBW(2*IELEM,LL4W),IPBBX(2*IELEM,LL4W),
     3 IPBBY(2*IELEM,LL4W),IPBBZ(2*IELEM,LL4Z)
      REAL SIDE,ZZ(3,NBLOS),FRZ(NBLOS),V(IELEM+1,IELEM),
     1 H(IELEM+1,IELEM),BBW(2*IELEM,LL4W),BBX(2*IELEM,LL4W),
     2 BBY(2*IELEM,LL4W),BBZ(2*IELEM,LL4Z)
      DOUBLE PRECISION CTRAN((IELEM+1)*IELEM,(IELEM+1)*IELEM)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION TTTT,DENOM,VOL0
*
      NELEH=(IELEM+1)*IELEM**2
      NELEZ=6*IELEM**2
      NBC=INT((SQRT(REAL((4*LXH-1)/3))+1.)/2.)
      IF(LL4F.GT.3*NBLOS*IELEM**3) CALL XABORT('TRICHH: BUG1.')
      IF(LL4W.GT.(2*NBLOS*IELEM+(2*NBC-1)*ISPLH*LZ)*IELEM**2)
     1 CALL XABORT('TRICHH: BUG2.')
*----
*  COMPUTE THE TRANVERSE COUPLING PIOLAT UNIT MATRIX
*----
      CALL XDDSET(CTRAN,((IELEM+1)*IELEM)**2,0.0D0)
      CNORM=SIDE*SIDE/SQRT(3.)
      I=0
      DO 22 JS=1,IELEM
      DO 21 JT=1,IELEM+1
      J=0
      I=I+1
      SSS=1.0
      DO 20 IT=1,IELEM
      DO 10 IS=1,IELEM+1
      J=J+1
      CTRAN(I,J)=SSS*CNORM*H(IS,JS)*H(JT,IT)
   10 CONTINUE
      SSS=-SSS
   20 CONTINUE
   21 CONTINUE
   22 CONTINUE
      IF(IMPX.GT.1) THEN
         WRITE(6,*) 'TRICHH: MATRIX CTRAN'
         DO 30 I=1,(IELEM+1)*IELEM
         WRITE(6,'(10(1X,1P,E12.4))') (CTRAN(I,J),J=1,(IELEM+1)*IELEM)
   30    CONTINUE
         WRITE(6,*) '  '
      ENDIF
*----
*  COMPUTE THE W-, X- ,Y- AND Z-ORIENTED SYSTEM BANDWIDTH VECTORS
*----
      CALL XDISET(MUW,L4,1)
      CALL XDISET(MUX,L4,1)
      CALL XDISET(MUY,L4,1)
      CALL XDISET(MUZ,L4,1)
      CALL XDISET(IPBBW,2*IELEM*LL4W,0)
      CALL XDISET(IPBBX,2*IELEM*LL4X,0)
      CALL XDISET(IPBBY,2*IELEM*LL4Y,0)
      CALL XDISET(IPBBZ,2*IELEM*LL4Z,0)
      NUM=0
      DO 80 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 80
      NUM=NUM+1
      DO 64 K5=0,1 ! TWO LOZENGES PER HEXAGON
      DO 63 K4=0,IELEM-1
      DO 62 K3=0,IELEM-1
      DO 61 K2=1,IELEM+1
      KNW1=KN(NUM,3+K5*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
      KNX1=KN(NUM,3+(K5+2)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
      KNY1=KN(NUM,3+(K5+4)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
      INW1=ABS(KNW1)
      INX1=ABS(KNX1)-LL4W
      INY1=ABS(KNY1)-LL4W-LL4X
      DO 40 K1=1,IELEM+1
      KNW2=KN(NUM,3+K5*NELEH+(K4*IELEM+K3)*(IELEM+1)+K1)
      KNX2=KN(NUM,3+(K5+2)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K1)
      KNY2=KN(NUM,3+(K5+4)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K1)
      INW2=ABS(KNW2)
      INX2=ABS(KNX2)-LL4W
      INY2=ABS(KNY2)-LL4W-LL4X
      IF((KNW2.NE.0).AND.(KNW1.NE.0)) THEN
         MUW(INW1)=MAX(MUW(INW1),INW1-INW2+1)
         MUW(INW2)=MAX(MUW(INW2),INW2-INW1+1)
      ENDIF
      IF((KNX2.NE.0).AND.(KNX1.NE.0)) THEN
         MUX(INX1)=MAX(MUX(INX1),INX1-INX2+1)
         MUX(INX2)=MAX(MUX(INX2),INX2-INX1+1)
      ENDIF
      IF((KNY2.NE.0).AND.(KNY1.NE.0)) THEN
         MUY(INY1)=MAX(MUY(INY1),INY1-INY2+1)
         MUY(INY2)=MAX(MUY(INY2),INY2-INY1+1)
      ENDIF
   40 CONTINUE
      DO 60 K1=0,IELEM-1
      IF(V(K2,K1+1).EQ.0.0) GO TO 60
      IF(K5.EQ.0) THEN
         JND1=(NUM-1)*IELEM**3+K4*IELEM**2+K3*IELEM+K1+1
         JND2=(KN(NUM,1)-1)*IELEM**3+K4*IELEM**2+K3*IELEM+K1+1
         JND3=(KN(NUM,2)-1)*IELEM**3+K4*IELEM**2+K3*IELEM+K1+1
      ELSE
         JND1=(KN(NUM,1)-1)*IELEM**3+K4*IELEM**2+K1*IELEM+K3+1
         JND2=(KN(NUM,2)-1)*IELEM**3+K4*IELEM**2+K1*IELEM+K3+1
         JND3=(KN(NUM,3)-1)*IELEM**3+K4*IELEM**2+K1*IELEM+K3+1
      ENDIF
      IF(KNW1.NE.0) CALL TRINDX(JND1,IPBBW(1,INW1),2*IELEM)
      IF(KNX1.NE.0) CALL TRINDX(JND2,IPBBX(1,INX1),2*IELEM)
      IF(KNY1.NE.0) CALL TRINDX(JND3,IPBBY(1,INY1),2*IELEM)
   60 CONTINUE
   61 CONTINUE
   62 CONTINUE
   63 CONTINUE
   64 CONTINUE
      DO 73 K5=0,2 ! THREE LOZENGES PER HEXAGON
      DO 72 K2=0,IELEM-1
      DO 71 K1=0,IELEM-1
      KNZ1=KN(NUM,3+6*NELEH+2*K5*IELEM**2+K2*IELEM+K1+1)
      KNZ2=KN(NUM,3+6*NELEH+(2*K5+1)*IELEM**2+K2*IELEM+K1+1)
      INZ1=ABS(KNZ1)-LL4W-LL4X-LL4Y
      INZ2=ABS(KNZ2)-LL4W-LL4X-LL4Y
      IF((KNZ1.NE.0).AND.(KNZ2.NE.0)) THEN
         MUZ(INZ2)=MAX(MUZ(INZ2),INZ2-INZ1+1)
         MUZ(INZ1)=MAX(MUZ(INZ1),INZ1-INZ2+1)
      ENDIF
      DO 70 K3=0,IELEM-1
      IF(K5.EQ.0) THEN
         JND1=(NUM-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
      ELSE
         JND1=(KN(NUM,K5)-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
      ENDIF
      IF(KNZ1.NE.0) CALL TRINDX(JND1,IPBBZ(1,INZ1),2*IELEM)
      IF(KNZ2.NE.0) CALL TRINDX(JND1,IPBBZ(1,INZ2),2*IELEM)
   70 CONTINUE
   71 CONTINUE
   72 CONTINUE
   73 CONTINUE
   80 CONTINUE
*
      MUWMAX=0
      IIMAXW=0
      DO 90 I=1,LL4W
      MUWMAX=MAX(MUWMAX,MUW(I))
      IIMAXW=IIMAXW+MUW(I)
      MUW(I)=IIMAXW
   90 CONTINUE
      MUXMAX=0
      IIMAXX=0
      DO 100 I=1,LL4X
      MUXMAX=MAX(MUXMAX,MUX(I))
      IIMAXX=IIMAXX+MUX(I)
      MUX(I)=IIMAXX
  100 CONTINUE
      MUYMAX=0
      IIMAXY=0
      DO 110 I=1,LL4Y
      MUYMAX=MAX(MUYMAX,MUY(I))
      IIMAXY=IIMAXY+MUY(I)
      MUY(I)=IIMAXY
  110 CONTINUE
      MUZMAX=0
      IIMAXZ=0
      DO 120 I=1,LL4Z
      MUZMAX=MAX(MUZMAX,MUZ(I))
      IIMAXZ=IIMAXZ+MUZ(I)
      MUZ(I)=IIMAXZ
  120 CONTINUE
      IF(IMPX.GT.0) THEN
         WRITE (6,600) MUWMAX,MUXMAX,MUYMAX,MUZMAX
         WRITE (6,610) IIMAXW,IIMAXX,IIMAXY,IIMAXZ
      ENDIF
*----
*  COMPUTE THE FLUX-CURRENT COUPLING MATRICES WB, XB, YB AND ZB.
*----
      CALL XDRSET(BBW,2*IELEM*LL4W,0.0)
      CALL XDRSET(BBX,2*IELEM*LL4X,0.0)
      CALL XDRSET(BBY,2*IELEM*LL4Y,0.0)
      CALL XDRSET(BBZ,2*IELEM*LL4Z,0.0)
      TTTT=0.5D0*SQRT(3.D00)*SIDE*SIDE
      DENOM=0.5D0*SQRT(3.D00)*SIDE
      NUM=0
      DO 260 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 260
      NUM=NUM+1
      DZ=ZZ(1,IPERT(KEL))
      VOL0=TTTT*DZ*FRZ(KEL)
      DO 194 K5=0,1
      DO 193 K4=0,IELEM-1
      DO 192 K3=0,IELEM-1
      DO 191 K2=1,IELEM+1
      KNW1=KN(NUM,3+K5*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
      KNX1=KN(NUM,3+(K5+2)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
      KNY1=KN(NUM,3+(K5+4)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
      INW1=ABS(KNW1)
      INX1=ABS(KNX1)-LL4W
      INY1=ABS(KNY1)-LL4W-LL4X
      DO 190 K1=0,IELEM-1
      IF(V(K2,K1+1).EQ.0.0) GO TO 190
      IF(K5.EQ.0) THEN
         SSS=(-1.0)**K1
         JND1=(NUM-1)*IELEM**3+K4*IELEM**2+K3*IELEM+K1+1
         JND2=(KN(NUM,1)-1)*IELEM**3+K4*IELEM**2+K3*IELEM+K1+1
         JND3=(KN(NUM,2)-1)*IELEM**3+K4*IELEM**2+K3*IELEM+K1+1
      ELSE
         SSS=1.0
         JND1=(KN(NUM,1)-1)*IELEM**3+K4*IELEM**2+K1*IELEM+K3+1
         JND2=(KN(NUM,2)-1)*IELEM**3+K4*IELEM**2+K1*IELEM+K3+1
         JND3=(KN(NUM,3)-1)*IELEM**3+K4*IELEM**2+K1*IELEM+K3+1
      ENDIF
      IF(KNW1.NE.0.0) THEN
         KK=0
         DO 130 I=1,2*IELEM
         IF(IPBBW(I,INW1).EQ.JND1) THEN
            KK=I
            GO TO 140
         ENDIF
  130    CONTINUE
         CALL XABORT('TRICHH: BUG3.')
  140    SG=REAL(SIGN(1,KNW1))
         BBW(KK,INW1)=BBW(KK,INW1)+SG*SSS*REAL(VOL0/DENOM)*V(K2,K1+1)
      ENDIF
      IF(KNX1.NE.0.0) THEN
         KK=0
         DO 150 I=1,2*IELEM
         IF(IPBBX(I,INX1).EQ.JND2) THEN
            KK=I
            GO TO 160
         ENDIF
  150    CONTINUE
         CALL XABORT('TRICHH: BUG4.')
  160    SG=REAL(SIGN(1,KNX1))
         BBX(KK,INX1)=BBX(KK,INX1)+SG*SSS*REAL(VOL0/DENOM)*V(K2,K1+1)
      ENDIF
      IF(KNY1.NE.0.0) THEN
         KK=0
         DO 170 I=1,2*IELEM
         IF(IPBBY(I,INY1).EQ.JND3) THEN
            KK=I
            GO TO 180
         ENDIF
  170    CONTINUE
         CALL XABORT('TRICHH: BUG5.')
  180    SG=REAL(SIGN(1,KNY1))
         BBY(KK,INY1)=BBY(KK,INY1)+SG*SSS*REAL(VOL0/DENOM)*V(K2,K1+1)
      ENDIF
  190 CONTINUE
  191 CONTINUE
  192 CONTINUE
  193 CONTINUE
  194 CONTINUE
      DO 253 K5=0,2 ! THREE LOZENGES PER HEXAGON
      DO 252 K2=0,IELEM-1
      DO 251 K1=0,IELEM-1
      KNZ1=KN(NUM,3+6*NELEH+2*K5*IELEM**2+K2*IELEM+K1+1)
      KNZ2=KN(NUM,3+6*NELEH+(2*K5+1)*IELEM**2+K2*IELEM+K1+1)
      INZ1=ABS(KNZ1)-LL4W-LL4X-LL4Y
      INZ2=ABS(KNZ2)-LL4W-LL4X-LL4Y
      DO 250 K3=0,IELEM-1
      IF(K5.EQ.0) THEN
         JND1=(NUM-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
      ELSE
         JND1=(KN(NUM,K5)-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
      ENDIF
      IF(KNZ1.NE.0) THEN
         KK=0
         DO 210 I=1,2*IELEM
         IF(IPBBZ(I,INZ1).EQ.JND1) THEN
            KK=I
            GO TO 220
         ENDIF
  210    CONTINUE
         CALL XABORT('TRICHH: BUG6.')
  220    SG=REAL(SIGN(1,KNZ1))
         BBZ(KK,INZ1)=BBZ(KK,INZ1)+SG*REAL(VOL0/DZ)*V(1,K3+1)
      ENDIF
      IF(KNZ2.NE.0) THEN
         KK=0
         DO 230 I=1,2*IELEM
         IF(IPBBZ(I,INZ2).EQ.JND1) THEN
            KK=I
            GO TO 240
         ENDIF
  230    CONTINUE
         CALL XABORT('TRICHH: BUG7.')
  240    SG=REAL(SIGN(1,KNZ2))
         BBZ(KK,INZ2)=BBZ(KK,INZ2)+SG*REAL(VOL0/DZ)*V(IELEM+1,K3+1)
      ENDIF
  250 CONTINUE
  251 CONTINUE
  252 CONTINUE
  253 CONTINUE
  260 CONTINUE
      RETURN
*
  600 FORMAT(/52H TRICHH: MAXIMUM BANDWIDTH FOR W-ORIENTED MATRICES =,
     1 I4/27X,25HFOR X-ORIENTED MATRICES =,I4/27X,16HFOR Y-ORIENTED M,
     2 9HATRICES =,I4/27X,25HFOR Z-ORIENTED MATRICES =,I4)
  610 FORMAT(/40H TRICHH: LENGTH OF W-ORIENTED MATRICES =,I10/16X,
     1 24HOF X-ORIENTED MATRICES =,I10/16X,24HOF Y-ORIENTED MATRICES =,
     2 I10/16X,24HOF Z-ORIENTED MATRICES =,I10)
      END
