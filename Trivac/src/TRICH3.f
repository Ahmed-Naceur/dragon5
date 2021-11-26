*DECK TRICH3
      SUBROUTINE TRICH3(ISPLH,IPTRK,LX,LZ,L4,MAT,KN,MUW,MUX,MUY,MUZ,
     1 IPW,IPX,IPY,IPZ,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the compressed diagonal storage indices (MUW, MUX, MUY and
* MUZ) an the permutation vectors (IPW, IPX, IPY and IPZ) for an ADI
* splitting of a mesh corner finite difference discretization in
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
* LX      number of hexagons in a plane.
* LZ      number of axial planes.
* L4      order of system matrices.
* MAT     mixture index assigned to each element.
* KN      element-ordered unknown list. Dimensionned to LL*LX*LZ
*         where LL=12 (hexagons) or 14 (triangles).
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
      INTEGER ISPLH,LX,LZ,L4,MAT(LX*LZ),KN(*),MUW(L4),MUX(L4),MUY(L4),
     1 MUZ(L4),IPW(L4),IPX(L4),IPY(L4),IPZ(L4),IMPX
*----
*  LOCAL VARIABLES
*----
      REAL HW(14,14),HX(14,14),HY(14,14),HZ(14,14),HL(2,2),RFAC(28,7),
     1 RF6(24,6),RF7(28,7)
      INTEGER NCODE(6),IJ1(14),IJ2(14),IJ27(14),IJ16(12),IJ26(12),
     1 IJ17(14)
      INTEGER, DIMENSION(:), ALLOCATABLE :: IDX,IDY
      COMMON /ELEMB/LC,T(5),TS(5),R(5,5),RS(5,5),Q(5,5),QS(5,5),V(5,4),
     1 E(5,5),RH(7,7),QH(7,7),RT(3,3),QT(3,3)
      DATA HL / 1.0,2*0.0,1.0/
      DATA IJ16,IJ26 /1,2,3,4,5,6,1,2,3,4,5,6,6*1,6*2/
      DATA IJ17,IJ27 /1,2,3,4,5,6,7,1,2,3,4,5,6,7,7*1,7*2/
      DATA RF6/
     >1.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,1.0,0.5,
     >1.0,0.0,1.0,1.0,0.0,0.5,1.0,0.0,0.0,0.0,0.0,0.0,
     >0.0,1.0,1.0,1.0,0.5,0.0,1.0,1.0,0.0,0.0,0.5,1.0,
     >0.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,
     >0.0,1.0,1.0,0.5,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,
     >1.0,0.0,1.0,0.5,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,
     >0.0,1.0,0.5,1.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,
     >1.0,0.0,0.5,1.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,
     >0.0,0.5,1.0,1.0,1.0,0.0,1.0,0.5,0.0,0.0,1.0,1.0,
     >0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,
     >0.0,0.0,0.0,0.0,0.0,1.0,0.5,1.0,0.0,0.0,1.0,1.0,
     >0.5,0.0,1.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0/
      DATA RF7/
     >1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.5,0.0,1.0,0.5,
     >1.0,0.0,1.0,0.5,1.0,0.0,0.5,1.0,0.0,0.0,0.0,0.0,0.0,0.0,
     >0.0,1.0,1.0,0.5,1.0,0.5,0.0,1.0,1.0,0.0,0.5,0.0,0.5,1.0,
     >0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,
     >0.0,1.0,1.0,0.5,0.5,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,
     >1.0,0.0,1.0,0.5,0.5,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,
     >0.0,0.5,0.5,1.0,0.5,0.5,0.0,0.5,0.5,0.0,1.0,0.0,0.5,0.5,
     >0.5,0.0,0.5,1.0,0.5,0.0,0.5,0.0,0.0,0.0,1.0,0.0,0.0,0.0,
     >0.0,1.0,0.5,0.5,1.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,
     >1.0,0.0,0.5,0.5,1.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,
     >0.0,0.5,1.0,0.5,1.0,1.0,0.0,1.0,0.5,0.0,0.5,0.0,1.0,1.0,
     >0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,
     >0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.5,1.0,0.0,0.5,0.0,1.0,1.0,
     >0.5,0.0,1.0,0.5,1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0/
*
      IF(ISPLH.EQ.1) THEN
         LC=6
         DO 10 I=1,2*LC
         IJ1(I)=IJ16(I)
         IJ2(I)=IJ26(I)
   10    CONTINUE
         DO 25 I=1,4*LC
         DO 20 J=1,LC
         RFAC(I,J)=RF6(I,J)
   20    CONTINUE
   25    CONTINUE
      ELSE
         LC=7
         DO 30 I=1,2*LC
         IJ1(I)=IJ17(I)
         IJ2(I)=IJ27(I)
   30    CONTINUE
         DO 45 I=1,4*LC
         DO 40 J=1,LC
         RFAC(I,J)=RF7(I,J)
   40    CONTINUE
   45    CONTINUE
      ENDIF
      LL=2*LC
      DO 55 I=1,LL
      I1=IJ1(I)
      I2=IJ2(I)
      DO 50 J=1,LL
      J1=IJ1(J)
      J2=IJ2(J)
      HW(I,J)   = RFAC(I1     ,J1) * HL(I2,J2)
      HX(I,J)   = RFAC(I1+LC  ,J1) * HL(I2,J2)
      HY(I,J)   = RFAC(I1+2*LC,J1) * HL(I2,J2)
      HZ(I,J)   = RFAC(I1+3*LC,J1)
   50 CONTINUE
   55 CONTINUE
*
      DO 65 I=1,LL
         I1 = IJ1(I)
         I2 = IJ2(I)
         DO 60 J=1,LL
            J1 = IJ1(J)
            J2 = IJ2(J)
            HW(I,J) = RFAC(I1     ,J1) * HL(I2,J2)
            HX(I,J) = RFAC(I1+LC  ,J1) * HL(I2,J2)
            HY(I,J) = RFAC(I1+2*LC,J1) * HL(I2,J2)
            HZ(I,J) = RFAC(I1+3*LC,J1)
  60     CONTINUE
  65  CONTINUE
*----
*  COMPUTE THE PERMUTATION VECTORS
*----
      DO 70 I=1,L4
      IPW(I)=I
      IPX(I)=0
      IPY(I)=0
      IPZ(I)=0
   70 CONTINUE
      LT4 = L4
      LPZ = LZ
      IF(LZ.GT.1) THEN
         LPZ = LZ+1
         CALL LCMGET (IPTRK,'NCODE',NCODE)
         IF((NCODE(5).EQ.7).OR.(NCODE(6).EQ.7))   LPZ = LZ
         IF((NCODE(5).EQ.7).AND.(NCODE(6).EQ.7)) LPZ = LZ-1
         LT4 = L4/LPZ
      ENDIF
      ALLOCATE(IDX(LT4),IDY(LT4))
      CALL LCMGET (IPTRK,'ILX',IDX)
      CALL LCMGET (IPTRK,'ILY',IDY)
      DO 85 KZ = 1, LPZ
      DO 80 KX = 1, LT4
         IPX(KX+(KZ-1)*LT4) = IDX(KX) + (KZ-1)*LT4
         IPY(KX+(KZ-1)*LT4) = IDY(KX) + (KZ-1)*LT4
   80 CONTINUE
   85 CONTINUE
      DEALLOCATE(IDY,IDX)
      KEL = 0
      DO 95 KX = 1, LT4
      DO 90 KZ = 1, LPZ
         KEL = KEL + 1
         IPZ(KX+(KZ-1)*LT4) = KEL
   90 CONTINUE
   95 CONTINUE
*----
*  COMPUTE THE COMPRESSED DIAGONAL STORAGE INDICES
*----
      DO 100 I=1,L4
      MUW(I)=1
      MUX(I)=1
      MUY(I)=1
      MUZ(I)=1
  100 CONTINUE
      NUM1=0
      DO 130 K=1,LX*LZ
         IF(MAT(K).LE.0) GO TO 130
         DO 120 I=1,LL
            INW1=KN(NUM1+I)
            IF(INW1.EQ.0) GO TO 120
            INX1=IPX(INW1)
            INY1=IPY(INW1)
            INZ1=IPZ(INW1)
            DO 110 J=1,LL
               INW2=KN(NUM1+J)
               IF(INW2.EQ.0) GO TO 110
               INX2=IPX(INW2)
               INY2=IPY(INW2)
               INZ2=IPZ(INW2)
               IF((HW(I,J).NE.0.0).AND.(INW2.LT.INW1))
     >            MUW(INW1)=MAX0(MUW(INW1),INW1-INW2+1)
               IF((HX(I,J).NE.0.0).AND.(INX2.LT.INX1))
     >            MUX(INX1)=MAX0(MUX(INX1),INX1-INX2+1)
               IF((HY(I,J).NE.0.0).AND.(INY2.LT.INY1))
     >            MUY(INY1)=MAX0(MUY(INY1),INY1-INY2+1)
               IF((HZ(I,J).NE.0.0).AND.(INZ2.LT.INZ1))
     >            MUZ(INZ1)=MAX0(MUZ(INZ1),INZ1-INZ2+1)
  110       CONTINUE
  120    CONTINUE
         NUM1=NUM1+LL
  130 CONTINUE
      IF(IMPX.GE.5) THEN
          WRITE(6,510) 'IPW :',(IPW(I),I=1,L4)
          WRITE(6,510) 'MUW :',(MUW(I),I=1,L4)
          WRITE(6,510) 'IPX :',(IPX(I),I=1,L4)
          WRITE(6,510) 'MUX :',(MUX(I),I=1,L4)
          WRITE(6,510) 'IPY :',(IPY(I),I=1,L4)
          WRITE(6,510) 'MUY :',(MUY(I),I=1,L4)
          IF(LZ.GT.1) THEN
             WRITE(6,510) 'IPZ :',(IPZ(I),I=1,L4)
             WRITE(6,510) 'MUZ :',(MUZ(I),I=1,L4)
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
      DO 140 I=1,L4
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
  140 CONTINUE
      IF(IMPX.GT.0) WRITE (6,500) MUWMAX,MUXMAX,MUYMAX,MUZMAX
      RETURN
*
  500 FORMAT(/41H TRICH3: MAXIMUM BANDWIDTH ALONG W AXIS =,I5/
     1 27X,14HALONG X AXIS =,I5/27X,14HALONG Y AXIS =,I5/27X,
     2 14HALONG Z AXIS =,I5)
  510 FORMAT(/1X,A5/(1X,20I6))
      END
