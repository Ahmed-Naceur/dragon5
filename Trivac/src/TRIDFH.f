*DECK TRIDFH
      SUBROUTINE TRIDFH (ISPLH,IPTRK,IDIM,LX,LZ,LL4,NUN,SIDE,ZZZ,ZZ,
     1 KN,QFR,IQFR,VOL,MAT,IDL,NCODE,ICODE,ZCODE,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a mesh centered finite difference
* discretization of a 3-D hexagonal geometry.
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
* IMPX    print parameter.
* ISPLH   type of mesh-splitting: =1 for complete hexagons; =2 for
*         triangular mesh-splitting.
* IPTRK   L_TRACK pointer to the tracking information.
* IDIM    number of dimensions (2 or 3).
* LX      number of hexagons.
* LZ      number of axial planes.
* SIDE    side of an hexagon.
* ZZZ     Z-coordinates of the axial planes.
* NCODE   type of boundary condition applied on each side (I=1: hbc):
*         NCODE(I)=1: VOID;          =2: REFL;       =6: ALBE;
*                 =5: SYME;          =7: ZERO.
* ICODE   physical albedo index on each side of the domain.
* ZCODE   albedo corresponding to boundary condition 'VOID' on each
*         side (ZCODE(I)=0.0 by default).
* MAT     mixture index assigned to each element.
*
*Parameters: output
* LL4     order of the system matrices.
* NUN     number of unknowns per energy group.
* VOL     volume of each element.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* IQFR    element-ordered physical albedo indices.
* IDL     position of the average flux component associated with each
*         volume.
* ZZ      Z-sides of each hexagon.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER ISPLH,IDIM,LX,LZ,LL4,NUN,KN(8*LX*LZ),IQFR(8*LX*LZ),
     1 MAT(LX*LZ),IDL(LX*LZ),NCODE(6),ICODE(6),IMPX
      REAL SIDE,ZZZ(LZ+1),ZZ(LX*LZ),QFR(8*LX*LZ),VOL(LX*LZ),ZCODE(6)
*----
*  LOCAL VARIABLES
*----
      LOGICAL LL1,LL2
      INTEGER, DIMENSION(:), ALLOCATABLE :: I1,KN1,KN2,KN3,KN4
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*----
*  MAIN LOOP OVER THE FINITE ELEMENTS
*----
      ALLOCATE(I1(LX*LZ),KN4(8*LX*LZ))
      MEL=0
      DO 10 KXZ=1,LX*LZ
         I1(KXZ) = 0
         IF(MAT(KXZ).LE.0) GO TO 10
         MEL=MEL+1
         I1(KXZ) = MEL
   10 CONTINUE
      IDEB = 0
      IFIN = 0
      IVAL=0
      NUM1=0
      MEL = MEL/LZ
      KEL =0
      DO 45 KZ=1,LZ
         LL1 = .FALSE.
         LL2 = .FALSE.
      DO 40 KX=1,LX
         KEL = KEL + 1
         ZZ(KEL) = 0.0
         IF(MAT(KEL).LE.0) GO TO 40
         ZZ(KEL) = ZZZ(KZ+1)-ZZZ(KZ)
         DO 20 IC=1,8
            QFR(NUM1+IC) = 0.0
            IQFR(NUM1+IC) = 0
   20    CONTINUE
         DO 30 IX=1,6
            KN4(NUM1+IX) = 0
            N1 = NEIGHB (KX,IX,9,LX,POIDS)
            IF(N1.GT.0) N1 = N1+(KZ-1)*LX
            IF(N1.GT.(LX+(KZ-1)*LX)) THEN
               IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
                  N1 = -1
                  QFR(NUM1+IX)=ALB(ZCODE(1))
               ELSE IF(NCODE(1).EQ.1) THEN
                  N1 = -1
                  QFR(NUM1+IX)=1.0
                  IQFR(NUM1+IX)=ICODE(1)
               ELSE IF(NCODE(1).EQ.2) THEN
                  N1 = -2
               ELSE IF(NCODE(1).EQ.7) THEN
                  N1 = -3
               ENDIF
            ELSE IF(MAT(N1).LE.0) THEN
               IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
                  N1 = -1
                  QFR(NUM1+IX)=ALB(ZCODE(1))
               ELSE IF(NCODE(1).EQ.1) THEN
                  N1 = -1
                  QFR(NUM1+IX)=1.0
                  IQFR(NUM1+IX)=ICODE(1)
               ELSE IF(NCODE(1).EQ.2) THEN
                  N1 = -2
               ELSE IF(NCODE(1).EQ.7) THEN
                  N1 = -3
               ENDIF
            ENDIF
            IF(N1.GT.0)   N1 = I1(N1)
            KN4(NUM1+IX) = N1
30       CONTINUE
         KK7 = I1(KEL) - MEL
         KK8 = I1(KEL) + MEL
*
*     VOID, REFL OR ZERO BOUNDARY CONDITIONS.
         IF(KZ.EQ.1) THEN
            LL1 = .TRUE.
         ENDIF
         IF(LL1) THEN
            IF((NCODE(5).EQ.1).AND.(ICODE(5).EQ.0)) THEN
               KK7=-1
               QFR(NUM1+7)=ALB(ZCODE(5))
            ELSE IF(NCODE(5).EQ.1) THEN
               KK7=-1
               QFR(NUM1+7)=1.0
               IQFR(NUM1+7)=ICODE(5)
            ELSE IF(NCODE(5).EQ.2) THEN
               KK7=-2
            ELSE IF(NCODE(5).EQ.7) THEN
               KK7=-3
            ENDIF
         ENDIF
*
         IF(KZ.EQ.LZ) THEN
            LL2 = .TRUE.
         ENDIF
         IF(LL2) THEN
            IF((NCODE(6).EQ.1).AND.(ICODE(6).EQ.0)) THEN
               KK8=-1
               QFR(NUM1+8)=ALB(ZCODE(6))
            ELSE IF(NCODE(6).EQ.1) THEN
               KK8=-1
               QFR(NUM1+8)=1.0
               IQFR(NUM1+8)=ICODE(6)
            ELSE IF(NCODE(6).EQ.2) THEN
               KK8=-2
            ELSE IF(NCODE(6).EQ.7) THEN
               KK8=-3
            ENDIF
         ENDIF
*
*     TRAN BOUNDARY CONDITION.
         IF((KZ.EQ.1).AND.(NCODE(5).EQ.4)) THEN
            KK7=-2
         ELSE IF((KZ.EQ.LZ).AND.(NCODE(6).EQ.4)) THEN
            KK8=-2
         ENDIF
*
*     SYME BOUNDARY CONDITION.
         IF((KZ.EQ.1).AND.(NCODE(5).EQ.5)) THEN
            KK7=-2
            ZZ(KEL)=0.5*ZZ(KEL)
         ELSE IF((KZ.EQ.LZ).AND.(NCODE(6).EQ.5)) THEN
            KK8=-2
            ZZ(KEL)=0.5*ZZ(KEL)
         ENDIF
*
         IF(KZ.EQ.1)  IDEB = KK7
         IF(KZ.EQ.LZ) IFIN = KK8
         KN4(NUM1+7)  = KK7
         KN4(NUM1+8)  = KK8
         NUM1=NUM1+8
40    CONTINUE
45    CONTINUE
* END OF THE MAIN LOOP OVER FINITE ELEMENTS.
*
*----
*  VOLUME CALCULATION
*----
      DO 55 KZ=1,LZ
      DO 50 KX=1,LX
         KEL = KX+(KZ-1)*LX
         IF(MAT(KEL).EQ.0) THEN
            VOL0 = 0.0
         ELSE
            VOL0 = 2.59807621*SIDE*SIDE*ZZ(KEL)
         ENDIF
         VOL(KEL) = VOL0
   50 CONTINUE
   55 CONTINUE
      IF(IMPX.NE.0) THEN
         WRITE(6,222) 'ZZ ',(ZZ(I),I=1,LX*LZ)
         WRITE(6,222) 'VOL',(VOL(I),I=1,LX*LZ)
      ENDIF
      LL4=0
      DO 60 KXZ=1,LX*LZ
         IF(MAT(KXZ).GT.0) LL4=LL4+1
   60 CONTINUE
*
      IF(ISPLH.EQ.1) THEN
         IVAL = 8
         DO 70 I=1,8*LL4
            KN(I) = 0
            KN(I) = KN4(I)
 70      CONTINUE
         DEALLOCATE(KN4)
      ELSE IF(ISPLH.GE.2) THEN
         IVAL =18*(ISPLH-1)**2+8
         CALL TRINTR(ISPLH,IPTRK,LX,LL4,9,MAT)
         ALLOCATE(KN1(LL4),KN2(LL4),KN3(LL4))
         CALL LCMGET (IPTRK,'IKN',KN1)
         NUM1 = 0
         NUM3 = 0
         DO 80 I=1,LZ*LX*IVAL
            KN(I) = 0
 80      CONTINUE
         DO 90 I=1,LL4
            KN2(I) = IDEB
            KN3(I) = IFIN
 90      CONTINUE
         DO 115 K=1,LZ
            NUM2 = 0
         DO 110 I=1,LX
            IF(MAT(I+(K-1)*LX).LE.0) GO TO 110
            DO 100 J=1,6*(ISPLH-1)**2
               IVAL1 = KN1(NUM2+J) + (K-1)*LL4
               IVAL2 = KN2(NUM2+J)
               IVAL3 = KN3(NUM2+J)
               IF(IDIM.EQ.3.AND.K.GT.1)   IVAL2 = IVAL1 - LL4
               IF(IDIM.EQ.3.AND.K.LT.LZ) IVAL3 = IVAL1 + LL4
               KN(NUM1+J                ) = IVAL1
               KN(NUM1+J+ 6*(ISPLH-1)**2) = IVAL2
               KN(NUM1+J+12*(ISPLH-1)**2) = IVAL3
 100        CONTINUE
            DO 105 KX=1,6
               KN(NUM1+KX+18*(ISPLH-1)**2) = KN4(NUM3+KX)
 105        CONTINUE
            KN(NUM1+IVAL-1) = KN4(NUM3+7)
            KN(NUM1+IVAL  ) = KN4(NUM3+8)
            NUM1 = NUM1 + IVAL
            NUM2 = NUM2 + 6*(ISPLH-1)**2
            NUM3 = NUM3 + 8
 110     CONTINUE
 115     CONTINUE
         LL4 = LZ * LL4
         DEALLOCATE(KN3,KN2,KN1,KN4)
      ENDIF
      IF(IMPX.GE.1) THEN
         NUM1=0
         NUM2=0
         IF(ISPLH.EQ.1) THEN
            WRITE (6,570)
            DO 130 KZ=1,LZ
               WRITE(6,'(/13H PLANE NUMBER,I6)') KZ
               WRITE (6,520)
               DO 120 KX=1,LX
                  KEL = KX+(KZ-1)*LX
                  IF(MAT(KEL).LE.0) GO TO 120
                  WRITE (6,530) I1(KEL),(KN(NUM1+I),I=1,IVAL),
     >                              (QFR(NUM2+I),I=1,8),VOL(KEL)
                  NUM1 = NUM1 + IVAL
                  NUM2 = NUM2 + 8
120            CONTINUE
130         CONTINUE
         ELSE
            WRITE (6,570)
            DO 160 KZ=1,LZ
               WRITE(6,'(/13H PLANE NUMBER,I6)') KZ
               WRITE (6,575)
               DO 140 KX=1,LX
                  KEL = KX+(KZ-1)*LX
                  IF(MAT(KEL).LE.0) GO TO 140
                  WRITE (6,580) I1(KEL),(KN(NUM1+I),I=1,IVAL-8)
                  NUM1 = NUM1 + IVAL
140            CONTINUE
               WRITE (6,585)
               DO 150 KX=1,LX
                  KEL = KX+(KZ-1)*LX
                  IF(MAT(KEL).LE.0) GO TO 150
                  WRITE (6,590) I1(KEL),(QFR(NUM2+I),I=1,8),
     *                          VOL(KEL)
                  NUM2 = NUM2 + 8
150            CONTINUE
160         CONTINUE
         ENDIF
         WRITE (6,560) LL4
      ENDIF
      DEALLOCATE(I1)
*----
*  APPEND THE AVERAGED FLUXES AT THE END OF UNKNOWN VECTOR
*----
      NUN=0
      IF(ISPLH.GT.1) NUN=LL4
      DO 190 I=1,LX*LZ
      IF(MAT(I).EQ.0) THEN
         IDL(I)=0
      ELSE
         NUN=NUN+1
         IDL(I)=NUN
      ENDIF
190   CONTINUE
      RETURN
*
222   FORMAT(1X,A3,/,7(2X,E12.5))
520   FORMAT (/8H ELEMENT,6X,10HNEIGHBOURS,37X,20HVOID BOUNDARY CONDIT,
     1 3HION,28X,6HVOLUME)
530   FORMAT (1X,I6,2X,8I6,2X,8F6.2,5X,E13.6)
560   FORMAT (/40H NUMBER OF NON VIRTUAL FINITE ELEMENTS =,I6/)
570   FORMAT (/22H NUMBERING OF UNKNOWNS/1X,21(1H-))
575   FORMAT (/8H ELEMENT,44X,10HNEIGHBOURS)
580   FORMAT (1X,I6,2X,20I6/(9X,20I6))
585   FORMAT (/8H ELEMENT,3X,23HVOID BOUNDARY CONDITION,28X,6HVOLUME)
590   FORMAT (1X,I6,2X,8F6.2,5X,E13.6)
      END
