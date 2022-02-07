*DECK TRICHP
      SUBROUTINE TRICHP(IEL,LX,LY,LZ,L4,MAT,KN,MUX,MUY,MUZ,IPY,IPZ,
     1 IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Primal finite element unknown numbering for ADI solution in a 3D
* domain. Compute the storage info for ADI matrices in compressed
* diagona storage mode. Compute the ADI permutation vectors.
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
* IEL     degree of the Lagrangian finite elements. =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* L4      total number of unknown (variational coefficients) per
*         energy group (order of system matrices).
* MAT     mixture index assigned to each element.
* KN      element-ordered unknown list.
*
*Parameters: output
* MUX     X-directed compressed diagonal storage mode indices.
* MUY     Y-directed compressed diagonal storage mode indices.
* MUZ     Z-directed compressed diagonal storage mode indices.
* IPY     Y-directed permutation vectors.
* IPZ     Z-directed permutation vectors.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IEL,LX,LY,LZ,L4,MAT(LX*LY*LZ),KN(LX*LY*LZ*(IEL+1)**3),
     1 MUX(L4),MUY(L4),MUZ(L4),IPY(L4),IPZ(L4),IMPX
*----
*  LOCAL VARIABLES
*----
      INTEGER IJ1(125),IJ2(125),IJ3(125)
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IWRK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IWRK(LX*IEL+1,LY*IEL+1,LZ*IEL+1))
*
      LC=IEL+1
      LL=LC*LC*LC
      DO 5 L=1,LL
      L1=1+MOD(L-1,LC)
      L2=1+(L-L1)/LC
      L3=1+MOD(L2-1,LC)
      IJ1(L)=L1
      IJ2(L)=L3
      IJ3(L)=1+(L2-L3)/LC
    5 CONTINUE
*----
*  JUXTAPOSITION OF A CHECKERBOARD OVER A PLANE IN THE REACTOR
*----
      L2M=0
      LZTOT=LZ*(LC-1)+1
      LYTOT=LY*(LC-1)+1
      LXTOT=LX*(LC-1)+1
      DO 12 K=1,LZTOT
      DO 11 J=1,LYTOT
      DO 10 I=1,LXTOT
      IWRK(I,J,K)=0
   10 CONTINUE
   11 CONTINUE
   12 CONTINUE
      NUM1=0
      KEL=0
      DO 32 K0=1,LZ
      LK0=(K0-1)*(LC-1)
      DO 31 K1=1,LY
      LK1=(K1-1)*(LC-1)
      DO 30 K2=1,LX
      KEL=KEL+1
      IF(MAT(KEL).EQ.0) GO TO 30
      L2M=L2M+1
      LK2=(K2-1)*(LC-1)
      L=0
      DO 22 IK0=LK0+1,LK0+LC
      DO 21 IK1=LK1+1,LK1+LC
      DO 20 IK2=LK2+1,LK2+LC
      L=L+1
      IND1=KN(NUM1+L)
      IF(IND1.EQ.0) GO TO 20
      IF(IWRK(IK2,IK1,IK0).EQ.0) THEN
         IWRK(IK2,IK1,IK0)=IND1
      ELSE IF(IWRK(IK2,IK1,IK0).NE.IND1) THEN
         CALL XABORT('TRICHP: FAILURE OF THE RENUMBERING ALGORITHM(1).')
      ENDIF
   20 CONTINUE
   21 CONTINUE
   22 CONTINUE
      NUM1=NUM1+LL
   30 CONTINUE
   31 CONTINUE
   32 CONTINUE
*----
*  CALCULATION OF PERMUTATION VECTORS IPY AND IPZ
*----
      DO 40 I=1,L4
      IPY(I)=0
      IPZ(I)=0
   40 CONTINUE
      INEW=0
      DO 52 K0=1,LZTOT
      DO 51 K2=1,LXTOT
      IF(IWRK(K2,1,K0).EQ.IWRK(K2,LC,K0)) THEN
         K1MIN=1+LC/2
      ELSE
         K1MIN=1
      ENDIF
      DO 50 K1=K1MIN,LYTOT
      I=IWRK(K2,K1,K0)
      IF(I.EQ.0) GO TO 50
      IF(IPY(I).EQ.0) THEN
         INEW=INEW+1
         IPY(I)=INEW
      ENDIF
   50 CONTINUE
   51 CONTINUE
   52 CONTINUE
      IF(INEW.NE.L4) THEN
         CALL XABORT('TRICHP: FAILURE OF THE RENUMBERING ALGORITHM(2).')
      ENDIF
      INEW=0
      DO 72 K1=1,LYTOT
      DO 71 K2=1,LXTOT
      IF(IWRK(K2,K1,1).EQ.IWRK(K2,K1,LC)) THEN
         K0MIN=1+LC/2
      ELSE
         K0MIN=1
      ENDIF
      DO 70 K0=K0MIN,LZTOT
      I=IWRK(K2,K1,K0)
      IF(I.EQ.0) GO TO 70
      IF(IPZ(I).EQ.0) THEN
         INEW=INEW+1
         IPZ(I)=INEW
      ENDIF
   70 CONTINUE
   71 CONTINUE
   72 CONTINUE
      IF(INEW.NE.L4) THEN
         CALL XABORT('TRICHP: FAILURE OF THE RENUMBERING ALGORITHM(3).')
      ENDIF
*----
*  CALCULATION OF VECTORS MUX, MUY AND MUZ
*----
      DO 100 I=1,L4
      MUX(I)=1
      MUY(I)=1
      MUZ(I)=1
  100 CONTINUE
      NUM1=0
      DO 130 K=1,L2M
      DO 120 I=1,LL
      INX1=KN(NUM1+I)
      IF(INX1.EQ.0) GO TO 120
      INY1=IPY(INX1)
      INZ1=IPZ(INX1)
      DO 110 J=1,LL
      INX2=KN(NUM1+J)
      IF(INX2.EQ.0) GO TO 110
      INY2=IPY(INX2)
      INZ2=IPZ(INX2)
      IF((IJ2(I).EQ.IJ2(J)).AND.(IJ3(I).EQ.IJ3(J)).AND.(INX2.LT.INX1))
     1 THEN
         MUX(INX1)=MAX0(MUX(INX1),INX1-INX2+1)
      ELSE IF((IJ1(I).EQ.IJ1(J)).AND.(IJ3(I).EQ.IJ3(J)).AND.
     1 (INY2.LT.INY1)) THEN
         MUY(INY1)=MAX0(MUY(INY1),INY1-INY2+1)
      ELSE IF((IJ1(I).EQ.IJ1(J)).AND.(IJ2(I).EQ.IJ2(J)).AND.
     1 (INZ2.LT.INZ1)) THEN
         MUZ(INZ1)=MAX0(MUZ(INZ1),INZ1-INZ2+1)
      ENDIF
  110 CONTINUE
  120 CONTINUE
      NUM1=NUM1+LL
  130 CONTINUE
*
      MUXMAX=0
      MUYMAX=0
      MUZMAX=0
      IIMAXX=0
      IIMAXY=0
      IIMAXZ=0
      DO 140 I=1,L4
      MUXMAX=MAX(MUXMAX,MUX(I))
      MUYMAX=MAX(MUYMAX,MUY(I))
      MUZMAX=MAX(MUZMAX,MUZ(I))
      IIMAXX=IIMAXX+MUX(I)
      MUX(I)=IIMAXX
      IIMAXY=IIMAXY+MUY(I)
      MUY(I)=IIMAXY
      IIMAXZ=IIMAXZ+MUZ(I)
      MUZ(I)=IIMAXZ
  140 CONTINUE
      IF(IMPX.GT.0) WRITE (6,500) MUXMAX,MUYMAX,MUZMAX
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IWRK)
      RETURN
*
  500 FORMAT(/52H TRICHP: MAXIMUM BANDWIDTH FOR X-ORIENTED MATRICES =,
     1 I4/27X,25HFOR Y-ORIENTED MATRICES =,I4/27X,16HFOR Z-ORIENTED M,
     2 9HATRICES =,I4)
      END
