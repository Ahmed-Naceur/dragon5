*DECK TRICH4
      SUBROUTINE TRICH4(ISPLH,IPTRK,IDIM,LX,LZ,L4,MAT,KN,MUW,MUX,MUY,
     1 MUZ,IPW,IPX,IPY,IPZ,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the compressed diagonal storage indices (MUW, MUX, MUY and
* MUZ) an the permutation vectors (IPW, IPX, IPY and IPZ) for an ADI
* splitting of a mesh centered finite difference discretization in
* hexagonal geometry.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Benaboud
*
*Parameters: input
* ISPLH   type of mesh-splitting: =1 for complete hexagons; =2 for
*         triangular mesh-splitting.
* IPTRK   L_TRACK pointer to the Trivac tracking information.
* IDIM    number of dimensions (2 or 3).
* LX      number of hexagons in a plane.
* LZ      number of axial planes.
* L4      order of system matrices.
* MAT     mixture index assigned to each element.
* KN      element-ordered unknown list. Dimensionned to 8*L4
*         for hexagons and to (18*(ISPLH-1)**2+3)*LX*LZ for
*         triangular mesh-splitting.
* IMPX    print parameter (equal to zero for no print).
*
*Parameters: output
* MUW     W-oriented compressed storage mode indices.
* MUX     X-oriented compressed storage mode indices.
* MUY     Y-oriented compressed storage mode indices.
* MUZ     Z-oriented compressed storage mode indices.
* IPW     W-oriented permutation matrices.
* IPX     X-oriented permutation matrices.
* IPY     Y-oriented permutation matrices.
* IPZ     Z-oriented permutation matrices.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER ISPLH,IDIM,LX,LZ,L4,MAT(LX*LZ),KN(*),MUW(L4),MUX(L4),
     1 MUY(L4),MUZ(L4),IPW(L4),IPX(L4),IPY(L4),IPZ(L4),IMPX
*----
*  LOCAL VARIABLES
*----
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWRK,I1,I2,I3,I4,IDX,IDY
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IWRK(LX*LZ))
*
      IF(ISPLH.EQ.1) THEN
         ALLOCATE(I1(LX),I2(LX),I3(LX),I4(LX))
         LT4=L4/LZ
         MEL = 0
         DO 250 KEL=1,LX
         I4(KEL) = 0
         IF(MAT(KEL).GT.0) THEN
            MEL = MEL + 1
            I4(KEL) = MEL
         ENDIF
  250    CONTINUE
*----
*  COMPUTE THE PERMUTATION VECTORS
*----
         DO 260 I=1,L4
         IPW(I)=0
         IPX(I)=0
         IPY(I)=0
         IPZ(I)=0
  260    CONTINUE
         NC = INT((SQRT(REAL((4*LX-1)/3))+1.)/2.)
         J1 = 2 + 3*(NC-1)*(NC-2)
         IF(NC.EQ.1) J1=1
         J2 = J1 + NC - 1
         J3 = J2 + NC - 1
         CALL BIVPER(J1,1,LX,LT4,I1,I4)
         CALL BIVPER(J2,2,LX,LT4,I2,I4)
         CALL BIVPER(J3,3,LX,LT4,I3,I4)
         KEL = 0
         DO 275 K0 = 1,LZ
         DO 270 K1 = 1,LT4
         KEL = KEL + 1
         IV = (K0-1)*LT4
         IPW(KEL) = I1(K1)+IV
         IPX(KEL) = I2(K1)+IV
         IPY(KEL) = I3(K1)+IV
  270    CONTINUE
  275    CONTINUE
         IF(IDIM.EQ.3) THEN
            JEL = 0
            DO 285 K1=1,LT4
            DO 280 K0=1,LZ
            JEL = JEL + 1
            IPZ((K0-1)*LT4+I1(K1)) = JEL
  280       CONTINUE
  285       CONTINUE
         ENDIF
         DEALLOCATE(I4,I3,I2,I1)
*
         DO 300 I=1,L4
         MUW(I)=0
         MUX(I)=0
         MUY(I)=0
         MUZ(I)=0
  300    CONTINUE
*----
*  COMPUTE THE COMPRESSED DIAGONAL STORAGE INDICES
*----
         NUM1=0
         DO 320 KEL=1,L4
         KK1=KN(NUM1+6)
         KK2=KN(NUM1+3)
         INW1=IPW(KEL)
         MUW(INW1)=1
         IF(KK1.GT.0) THEN
            INW2=IPW(KK1)
            MUW(INW1)=MAX0(MUW(INW1),INW1-INW2+1)
         ENDIF
         IF(KK2.GT.0) THEN
            INW2=IPW(KK2)
            MUW(INW1)=MAX0(MUW(INW1),INW1-INW2+1)
         ENDIF
         NUM1=NUM1+8
  320    CONTINUE
*
         NUM1=0
         DO 330 KEL=1,L4
         KK3=KN(NUM1+1)
         KK4=KN(NUM1+4)
         INX1=IPX(KEL)
         MUX(INX1)=1
         IF(KK3.GT.0) THEN
            INX2=IPX(KK3)
            MUX(INX1)=MAX0(MUX(INX1),INX1-INX2+1)
         ENDIF
         IF(KK4.GT.0) THEN
            INX2=IPX(KK4)
            MUX(INX1)=MAX0(MUX(INX1),INX1-INX2+1)
         ENDIF
         NUM1=NUM1+8
  330    CONTINUE
*
         NUM1=0
         DO 340 KEL=1,L4
         KK5=KN(NUM1+2)
         KK6=KN(NUM1+5)
         INY1=IPY(KEL)
         MUY(INY1)=1
         IF(KK5.GT.0) THEN
            INY2=IPY(KK5)
            MUY(INY1)=MAX0(MUY(INY1),INY1-INY2+1)
         ENDIF
         IF(KK6.GT.0) THEN
            INY2=IPY(KK6)
            MUY(INY1)=MAX0(MUY(INY1),INY1-INY2+1)
         ENDIF
         NUM1=NUM1+8
  340    CONTINUE
         IF(IDIM.EQ.3) THEN
            NUM1=0
            DO 350 KEL=1,L4
            KK7=KN(NUM1+7)
            KK8=KN(NUM1+8)
            INZ1=IPZ(KEL)
            MUZ(INZ1)=1
            IF(KK7.GT.0) THEN
               INZ2=IPZ(KK7)
               MUZ(INZ1)=MAX0(MUZ(INZ1),INZ1-INZ2+1)
            ENDIF
            IF(KK8.GT.0) THEN
               INZ2=IPZ(KK8)
               MUZ(INZ1)=MAX0(MUZ(INZ1),INZ1-INZ2+1)
            ENDIF
            NUM1=NUM1+8
  350       CONTINUE
         ENDIF
*
      ELSE IF(ISPLH.GE.2) THEN
*
         NTPH = 6*(ISPLH-1)**2
         NTPL = 1+2*(ISPLH-1)
         NVT1 = NTPL + 2 * (ISPLH-2) + NTPH / 2
         NVT2 = NTPH - NTPL - (ISPLH-4) * (NTPL+2)
         NVT3 = NTPH - (ISPLH-4) * NTPL
         IVAL = 3*NTPH + 8
         IF(ISPLH.EQ.3) NVT2 = NTPH
         IF(ISPLH.LE.3) ISAU =   2*(ISPLH-2)
         IF(ISPLH.GE.4) ISAU =   6*(ISPLH-3)
         ICR = ISAU*(1+2*(ISPLH-2))
*----
*  COMPUTE THE PERMUTATION VECTORS.
*----
         DO 400 I=1,L4
         IPW(I)=I
         IPX(I)=0
         IPY(I)=0
         IPZ(I)=0
  400    CONTINUE
         LI4 = L4/LZ
         ALLOCATE(IDX(LI4),IDY(LI4))
         CALL LCMGET(IPTRK,'ILX',IDX)
         CALL LCMGET(IPTRK,'ILY',IDY)
         DO 415 KZ=1,LZ
         DO 410 KI=1,LI4
         IPX(KI+(KZ-1)*LI4) = IDX(KI) + (KZ-1)*LI4
         IPY(KI+(KZ-1)*LI4) = IDY(KI) + (KZ-1)*LI4
  410    CONTINUE
  415    CONTINUE
         DEALLOCATE(IDY,IDX)
         IF(IDIM.EQ.3) THEN
            DO 425 K1=1,LI4
            DO 420 K0=1,LZ
            IPZ((K0-1)*LI4+K1) = K0 + (K1-1)*LZ
  420       CONTINUE
  425       CONTINUE
         ENDIF
*
         DO 500 I=1,L4
         MUW(I)=0
         MUX(I)=0
         MUY(I)=0
         MUZ(I)=0
  500    CONTINUE
*----
*  COMPUTE THE COMPRESSED DIAGONAL STORAGE INDICES
*----
         NUM1=0
         DO 520 K0=1,LX*LZ
         IF(MAT(K0).LE.0) GO TO 520
         DO 510 I = 1,NTPH
         CALL TRINEI(3,1,2,ISPLH,ICR,I,KK1,KK2,KK3,KEL,IQF,NUM1,NTPH,
     >   NTPL,NVT1,NVT2,NVT3,IVAL,KN)
         INW1=IPW(KEL)
         MUW(INW1)=1
         IF(KK1.GT.0) THEN
            INW2=IPW(KK1)
            MUW(INW1)=MAX0(MUW(INW1),INW1-INW2+1)
         ENDIF
         IF(KK2.GT.0) THEN
            INW2=IPW(KK2)
            MUW(INW1)=MAX0(MUW(INW1),INW1-INW2+1)
         ENDIF
  510    CONTINUE
         NUM1=NUM1+IVAL
  520    CONTINUE
*
         NUM1=0
         DO 540 K0=1,LX*LZ
         IF(MAT(K0).LE.0) GO TO 540
         DO 530 I = 1,NTPH
         CALL TRINEI(3,2,2,ISPLH,ICR,I,KK1,KK2,KK3,KEL,IQF,NUM1,NTPH,
     >   NTPL,NVT1,NVT2,NVT3,IVAL,KN)
         INX1=IPX(KEL)
         MUX(INX1)=1
         IF(KK1.GT.0) THEN
            INX2=IPX(KK1)
            MUX(INX1)=MAX0(MUX(INX1),INX1-INX2+1)
         ENDIF
         IF(KK2.GT.0) THEN
            INX2=IPX(KK2)
            MUX(INX1)=MAX0(MUX(INX1),INX1-INX2+1)
         ENDIF
  530    CONTINUE
         NUM1=NUM1+IVAL
  540    CONTINUE
*
         NUM1=0
         DO 560 K0=1,LX*LZ
         IF(MAT(K0).LE.0) GO TO 560
         DO 550 I = 1,NTPH
         CALL TRINEI(3,3,2,ISPLH,ICR,I,KK1,KK2,KK3,KEL,IQF,NUM1,NTPH,
     >   NTPL,NVT1,NVT2,NVT3,IVAL,KN)
         INY1=IPY(KEL)
         MUY(INY1)=1
         IF(KK1.GT.0) THEN
            INY2=IPY(KK1)
            MUY(INY1)=MAX0(MUY(INY1),INY1-INY2+1)
         ENDIF
         IF(KK2.GT.0) THEN
            INY2=IPY(KK2)
            MUY(INY1)=MAX0(MUY(INY1),INY1-INY2+1)
         ENDIF
  550    CONTINUE
         NUM1=NUM1+IVAL
  560    CONTINUE
         IF(IDIM.EQ.3) THEN
*
            NUM1=0
            DO 580 K0=1,LX*LZ
            IF(MAT(K0).LE.0) GO TO 580
            DO 570 I = 1,NTPH
            KK1 = KN(NUM1+NTPH+I)
            KK2 = KN(NUM1+2*NTPH+I)
            KEL = KN(NUM1+I)
            INZ1=IPZ(KEL)
            MUZ(INZ1)=1
            IF(KK1.GT.0) THEN
               INZ2=IPZ(KK1)
               MUZ(INZ1)=MAX0(MUZ(INZ1),INZ1-INZ2+1)
            ENDIF
            IF(KK2.GT.0) THEN
               INZ2=IPZ(KK2)
               MUZ(INZ1)=MAX0(MUZ(INZ1),INZ1-INZ2+1)
            ENDIF
  570       CONTINUE
            NUM1=NUM1+IVAL
  580       CONTINUE
         ENDIF
      ENDIF
      IF(IMPX.GE.4) THEN
         WRITE(6,710) 'IPW :',(IPW(I),I=1,L4)
         WRITE(6,710) 'MUW :',(MUW(I),I=1,L4)
         WRITE(6,710) 'IPX :',(IPX(I),I=1,L4)
         WRITE(6,710) 'MUX :',(MUX(I),I=1,L4)
         WRITE(6,710) 'IPY :',(IPY(I),I=1,L4)
         WRITE(6,710) 'MUY :',(MUY(I),I=1,L4)
         IF(IDIM.EQ.3) THEN
            WRITE(6,710) 'IPZ :',(IPZ(I),I=1,L4)
            WRITE(6,710) 'MUZ :',(MUZ(I),I=1,L4)
         ENDIF
      ENDIF
*
      MUWMAX=0
      MUXMAX=0
      MUYMAX=0
      MUZMAX=0
      IIMAWW=0
      IIMAWX=0
      IIMAWY=0
      IIMAWZ=0
      DO 590 I=1,L4
      MUWMAX=MAX(MUWMAX,MUW(I))
      MUXMAX=MAX(MUXMAX,MUX(I))
      MUYMAX=MAX(MUYMAX,MUY(I))
      MUZMAX=MAX(MUZMAX,MUZ(I))
      IIMAWW=IIMAWW+MUW(I)
      MUW(I)=IIMAWW
      IIMAWX=IIMAWX+MUX(I)
      MUX(I)=IIMAWX
      IIMAWY=IIMAWY+MUY(I)
      MUY(I)=IIMAWY
      IIMAWZ=IIMAWZ+MUZ(I)
      MUZ(I)=IIMAWZ
 590  CONTINUE
      IF(IMPX.GE.0) WRITE (6,720) MUWMAX,MUXMAX,MUYMAX,MUZMAX
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IWRK)
      RETURN
*
  710 FORMAT(/1X,A5/(1X,20I6))
  720 FORMAT(/41H TRICH4: MAXIMUM BANDWIDTH ALONG W AXIS =,I5/
     1 27X,14HALONG X AXIS =,I5/27X,14HALONG Y AXIS =,I5/27X,
     2 14HALONG Z AXIS =,I5)
      END
