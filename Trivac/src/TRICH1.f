*DECK TRICH1
      SUBROUTINE TRICH1(IELEM,IDIM,LX,LY,LZ,L4,MAT,KN,MUX,MUY,MUZ,IPY,
     1 IPZ,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute the compressed diagonal storage indices (MUX, MUY and MUZ)
* and the permutation vectors (IPY and IPZ) for an ADI splitting of
* the nodal collocation leakage matrices.
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
* IELEM   degree of the polynomial expansion: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* IDIM    number of dimensions.
* LX      number of mesh along the X axis.
* LY      number of mesh along the Y axis.
* LZ      number of mesh along the Z axis.
* L4      order of system matrices
* MAT     mixture index assigned to each element.
* KN      element-ordered unknown list.
* IMPX    print parameter (equal to zero for no print).
*
*Parameters: output
* MUX     X-oriented compressed storage mode indices.
* MUY     Y-oriented compressed storage mode indices.
* MUZ     Z-oriented compressed storage mode indices.
* IPY     Y-oriented permutation matrices.
* IPZ     Z-oriented permutation matrices.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IELEM,IDIM,LX,LY,LZ,L4,MAT(LX*LY*LZ),KN(7*LX*LY*LZ),
     1 MUX(L4),MUY(L4),MUZ(L4),IPY(L4),IPZ(L4),IMPX
*----
*  LOCAL VARIABLES
*----
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWRK
*
      IORD(J,K,L,LL,IEL,IW)=(IEL*L+K)*LL*IEL+(1+IEL*(IW-1))+J
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IWRK(LX*LY*LZ))
*
      CALL XDISET(IWRK,LX*LY*LZ,0)
      LL4=0
      KEL=0
      DO 22 K0=1,LZ
      DO 21 K1=1,LY
      DO 20 K2=1,LX
      KEL=KEL+1
      IF(MAT(KEL).EQ.0) GO TO 20
      LL4=LL4+1
      IWRK((K0-1)*LX*LY+(K1-1)*LX+K2)=LL4
   20 CONTINUE
   21 CONTINUE
   22 CONTINUE
*----
*  COMPUTE THE PERMUTATION VECTORS IPY AND IPZ
*----
      IF(IDIM.GE.2) THEN
         INX1=0
         DO 52 K0=1,LZ
         DO 51 K2=1,LX
         DO 50 K1=1,LY
         INX2=IWRK((K0-1)*LX*LY+(K1-1)*LX+K2)
         IF(INX2.LE.0) GO TO 50
         INX1=INX1+1
         IF(IDIM.EQ.2) THEN
            DO 31 K=0,IELEM-1
            DO 30 J=0,IELEM-1
            I=IORD(J,K,0,LL4,IELEM,INX1)
            IPY(IORD(K,J,0,LL4,IELEM,INX2))=I
   30       CONTINUE
   31       CONTINUE
         ELSE IF(IDIM.EQ.3) THEN
            DO 42 L=0,IELEM-1
            DO 41 K=0,IELEM-1
            DO 40 J=0,IELEM-1
            I=IORD(J,K,L,LL4,IELEM,INX1)
            IPY(IORD(K,J,L,LL4,IELEM,INX2))=I
   40       CONTINUE
   41       CONTINUE
   42       CONTINUE
         ENDIF
   50    CONTINUE
   51    CONTINUE
   52    CONTINUE
         IF(INX1.NE.LL4) CALL XABORT('TRICH1: FAILURE OF THE RENUMBERI'
     1   //'NG ALGORITHM(1)')
         IF(IDIM.EQ.3) THEN
            INX1=0
            DO 72 K1=1,LY
            DO 71 K2=1,LX
            DO 70 K0=1,LZ
            INX2=IWRK((K0-1)*LX*LY+(K1-1)*LX+K2)
            IF(INX2.LE.0) GO TO 70
            INX1=INX1+1
            DO 62 L=0,IELEM-1
            DO 61 K=0,IELEM-1
            DO 60 J=0,IELEM-1
            I=IORD(J,K,L,LL4,IELEM,INX1)
            IPZ(IORD(K,L,J,LL4,IELEM,INX2))=I
   60       CONTINUE
   61       CONTINUE
   62       CONTINUE
   70       CONTINUE
   71       CONTINUE
   72       CONTINUE
            IF(INX1.NE.LL4) CALL XABORT('TRICH1: FAILURE OF THE RENUMB'
     1      //'ERING ALGORITHM(2)')
         ENDIF
      ENDIF
*
      L2M=0
      DO 80 KEL=1,LX*LY*LZ
      IF(MAT(KEL).EQ.0) GO TO 80
      L2M=L2M+1
      IWRK(KEL)=L2M
   80 CONTINUE
      DO 90 I=1,L4
      MUY(I)=0
      MUZ(I)=0
   90 CONTINUE
      LL5=L4/IELEM**(IDIM-1)
*----
*  COMPUTE VECTOR MUX
*----
      NUM1=0
      DO 130 KEL=1,LL4
      KK1=KN(NUM1+1)
      KK2=KN(NUM1+2)
      DO 100 J=0,IELEM-1
      INX1=IORD(J,0,0,LL4,IELEM,KEL)
      MUX(INX1)=J+1
*     X- SIDE:
      IF(KK1.GT.0) THEN
         INX2=IORD(0,0,0,LL4,IELEM,IWRK(KK1))
         MUX(INX1)=MAX0(MUX(INX1),INX1-INX2+1)
      ENDIF
*     X+ SIDE:
      IF(KK2.GT.0) THEN
         INX2=IORD(0,0,0,LL4,IELEM,IWRK(KK2))
         MUX(INX1)=MAX0(MUX(INX1),INX1-INX2+1)
      ENDIF
  100 CONTINUE
      NUM1=NUM1+6
  130 CONTINUE
*----
*  COMPUTE VECTOR MUY
*----
      IF(IDIM.GE.2) THEN
         NUM1=0
         DO 160 KEL=1,LL4
         KK3=KN(NUM1+3)
         KK4=KN(NUM1+4)
         DO 140 K=0,IELEM-1
         INY1=IPY(IORD(0,K,0,LL4,IELEM,KEL))
         MUY(INY1)=K+1
*        Y- SIDE:
         IF(KK3.GT.0) THEN
            INY2=IPY(IORD(0,0,0,LL4,IELEM,IWRK(KK3)))
            MUY(INY1)=MAX0(MUY(INY1),INY1-INY2+1)
         ENDIF
*        Y+ SIDE:
         IF(KK4.GT.0) THEN
            INY2=IPY(IORD(0,0,0,LL4,IELEM,IWRK(KK4)))
            MUY(INY1)=MAX0(MUY(INY1),INY1-INY2+1)
         ENDIF
  140    CONTINUE
         NUM1=NUM1+6
  160    CONTINUE
*----
*  COMPUTE VECTOR MUZ
*----
         IF(IDIM.EQ.3) THEN
            NUM1=0
            DO 180 KEL=1,LL4
            KK5=KN(NUM1+5)
            KK6=KN(NUM1+6)
            DO 170 L=0,IELEM-1
            INZ1=IPZ(IORD(0,0,L,LL4,IELEM,KEL))
            MUZ(INZ1)=L+1
*           Z- SIDE:
            IF(KK5.GT.0) THEN
               INZ2=IPZ(IORD(0,0,0,LL4,IELEM,IWRK(KK5)))
               MUZ(INZ1)=MAX0(MUZ(INZ1),INZ1-INZ2+1)
            ENDIF
*           Z+ SIDE:
            IF(KK6.GT.0) THEN
               INZ2=IPZ(IORD(0,0,0,LL4,IELEM,IWRK(KK6)))
               MUZ(INZ1)=MAX0(MUZ(INZ1),INZ1-INZ2+1)
            ENDIF
  170       CONTINUE
            NUM1=NUM1+6
  180       CONTINUE
            DO 195 J=1,IELEM-1
            DO 190 I=1,LL5
            MUX(I+J*LL5)=MUX(I)
            MUY(I+J*LL5)=MUY(I)
            MUZ(I+J*LL5)=MUZ(I)
  190       CONTINUE
  195       CONTINUE
            LL5=IELEM*LL5
         ENDIF
         DO 205 J=1,IELEM-1
         DO 200 I=1,LL5
         MUX(I+J*LL5)=MUX(I)
         MUY(I+J*LL5)=MUY(I)
         MUZ(I+J*LL5)=MUZ(I)
  200    CONTINUE
  205    CONTINUE
      ENDIF
*
      MUXMAX=0
      MUYMAX=0
      MUZMAX=0
      IIMAXX=0
      IIMAXY=0
      IIMAXZ=0
      DO 210 I=1,L4
      MUXMAX=MAX(MUXMAX,MUX(I))
      MUYMAX=MAX(MUYMAX,MUY(I))
      MUZMAX=MAX(MUZMAX,MUZ(I))
      IIMAXX=IIMAXX+MUX(I)
      MUX(I)=IIMAXX
      IIMAXY=IIMAXY+MUY(I)
      MUY(I)=IIMAXY
      IIMAXZ=IIMAXZ+MUZ(I)
      MUZ(I)=IIMAXZ
  210 CONTINUE
      IF(IMPX.GT.0) WRITE (6,230) MUXMAX,MUYMAX,MUZMAX
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IWRK)
      RETURN
*
  230 FORMAT(/41H TRICH1: MAXIMUM BANDWIDTH ALONG X AXIS =,I5/
     1 27X,14HALONG Y AXIS =,I5/27X,14HALONG Z AXIS =,I5)
      END
