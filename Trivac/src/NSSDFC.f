*DECK NSSDFC
      SUBROUTINE NSSDFC(IMPX,LX,LY,LZ,NCODE,ICODE,ZCODE,MAT,XXX,YYY,
     1 ZZZ,LL4F,LL4X,LL4Y,LL4Z,VOL,XX,YY,ZZ,IDL,KN,QFR,IQFR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a coarse mesh finite difference (NEM
* type) in a 3-D geometry.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
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
* NCODE   type of boundary condition applied on each side:
*         I=1: X-; I=2: X+; I=3: Y-; I=4: Y+; I=5: Z-; I=6: Z+;
*         NCODE(I)=1: VOID;  NCODE(I)=2: REFL;  NCODE(I)=4: TRAN;
*         NCODE(I)=7: ZERO.
* ICODE   physical albedo index on each side of the domain.
* ZCODE   albedo corresponding to boundary condition 'VOID' on each
*         side (ZCODE(i)=0.0 by default).
* MAT     mixture index assigned to each element.
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
* ZZZ     Cartesian coordinates along the Z axis.
*
*Parameters: output
* LL4F    total number of averaged flux unknown per energy group.
* LL4X    total number of X-direccted interface net currents.
* LL4Y    total number of Y-direccted interface net currents.
* LL4Z    total number of Z-direccted interface net currents.
* VOL     volume of each element.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* ZZ      Z-directed mesh spacings.
* IDL     position of averaged fluxes in unknown vector.
* KN      element-ordered interface net current unknown list.
* QFR     element-ordered boundary conditions.
* IQFR    element-ordered physical albedo indices.
*
*-----------------------------------------------------------------------
*
      INTEGER IMPX,LX,LY,LZ,NCODE(6),ICODE(6),MAT(LX,LY,LZ),LL4F,
     1 LL4X,LL4Y,LL4Z,IDL(LX,LY,LZ),KN(6,LX,LY,LZ),IQFR(6,LX,LY,LZ)
      REAL ZCODE(6),XXX(LX+1),YYY(LY+1),ZZZ(LZ+1),VOL(LX,LY,LZ),
     1 XX(LX,LY,LZ),YY(LX,LY,LZ),ZZ(LX,LY,LZ),QFR(6,LX,LY,LZ)
*----
*  LOCAL VARIABLES
*----
      LOGICAL LL1,LALB
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPX,IPY,IPZ
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*----
*  IDENTIFICATION OF THE NON VIRTUAL NODES
*----
      IF(IMPX.GT.0) WRITE(6,700) LX,LY,LZ
      ALLOCATE(IPX((LX+1)*LY*LZ),IPY((LY+1)*LX*LZ),IPZ((LZ+1)*LX*LY))
      IPX(:)=0
      IPY(:)=0
      IPZ(:)=0
      LL4F=0
      DO K0=1,LZ
        DO K1=1,LY
          DO K2=1,LX
            IDL(K2,K1,K0)=0
            KN(:6,K2,K1,K0)=0
            IF(MAT(K2,K1,K0).EQ.0) CYCLE
            LL4F=LL4F+1
            IDL(K2,K1,K0)=LL4F
            KN(1,K2,K1,K0)=K2+(LX+1)*(K1-1)+(LX+1)*LY*(K0-1)
            KN(2,K2,K1,K0)=(K2+1)+(LX+1)*(K1-1)+(LX+1)*LY*(K0-1)
            KN(3,K2,K1,K0)=K1+(LY+1)*(K2-1)+(LY+1)*LX*(K0-1)
            KN(4,K2,K1,K0)=(K1+1)+(LY+1)*(K2-1)+(LY+1)*LX*(K0-1)
            KN(5,K2,K1,K0)=K2+LX*(K1-1)+LX*LY*(K0-1)
            KN(6,K2,K1,K0)=K2+LX*(K1-1)+LX*LY*K0
            IPX(KN(1:2,K2,K1,K0))=1
            IPY(KN(3:4,K2,K1,K0))=1
            IPZ(KN(5:6,K2,K1,K0))=1
          ENDDO
        ENDDO
      ENDDO
      LL4X=0
      DO I=1,(LX+1)*LY*LZ
        IF(IPX(I).EQ.1) THEN
          LL4X=LL4X+1
          IPX(I)=LL4X
        ENDIF
      ENDDO
      LL4Y=0
      DO I=1,(LY+1)*LX*LZ
        IF(IPY(I).EQ.1) THEN
          LL4Y=LL4Y+1
          IPY(I)=LL4Y
        ENDIF
      ENDDO
      LL4Z=0
      DO I=1,(LZ+1)*LX*LY
        IF(IPZ(I).EQ.1) THEN
          LL4Z=LL4Z+1
          IPZ(I)=LL4Z
        ENDIF
      ENDDO
      DO K0=1,LZ
        DO K1=1,LY
          DO K2=1,LX
            IF(MAT(K2,K1,K0).EQ.0) CYCLE
            KN(1:2,K2,K1,K0)=IPX(KN(1:2,K2,K1,K0))
            KN(3:4,K2,K1,K0)=IPY(KN(3:4,K2,K1,K0))
            KN(5:6,K2,K1,K0)=IPZ(KN(5:6,K2,K1,K0))
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(IPZ,IPY,IPX)
*----
*  IDENTIFICATION OF THE GEOMETRY. MAIN LOOP OVER THE NODES
*----
      DO K0=1,LZ
        DO K1=1,LY
          DO K2=1,LX
            KEL=KEL+1
            XX(K2,K1,K0)=0.0
            YY(K2,K1,K0)=0.0
            ZZ(K2,K1,K0)=0.0
            VOL(K2,K1,K0)=0.0
            IF(MAT(K2,K1,K0).LE.0) CYCLE
            XX(K2,K1,K0)=XXX(K2+1)-XXX(K2)
            YY(K2,K1,K0)=YYY(K1+1)-YYY(K1)
            ZZ(K2,K1,K0)=ZZZ(K0+1)-ZZZ(K0)
            QFR(:6,K2,K1,K0)=0.0
            IQFR(:6,K2,K1,K0)=0
*----
*  VOID, REFL OR ZERO BOUNDARY CONTITION
*----
            IF(K2.EQ.1) THEN
              LL1=.TRUE.
            ELSE
              LL1=(MAT(K2-1,K1,K0).EQ.0)
            ENDIF
            IF(LL1) THEN
              LALB=(NCODE(1).EQ.1).OR.(NCODE(1).EQ.6)
              IF(LALB.AND.(ICODE(1).EQ.0)) THEN
                QFR(1,K2,K1,K0)=ALB(ZCODE(1))
                IQFR(1,K2,K1,K0)=-1
              ELSE IF(LALB) THEN
                QFR(1,K2,K1,K0)=1.0
                IQFR(1,K2,K1,K0)=ICODE(1)
              ELSE IF(NCODE(1).EQ.2) THEN
                IQFR(1,K2,K1,K0)=-2
              ELSE IF(NCODE(1).EQ.7) THEN
                IQFR(1,K2,K1,K0)=-3
              ENDIF
            ENDIF
*
            IF(K2.EQ.LX) THEN
              LL1=.TRUE.
            ELSE
              LL1=(MAT(K2+1,K1,K0).EQ.0)
            ENDIF
            IF(LL1) THEN
              LALB=(NCODE(2).EQ.1).OR.(NCODE(2).EQ.6)
              IF(LALB.AND.(ICODE(2).EQ.0)) THEN
                QFR(2,K2,K1,K0)=ALB(ZCODE(2))
                IQFR(2,K2,K1,K0)=-1
              ELSE IF(LALB) THEN
                QFR(2,K2,K1,K0)=1.0
                IQFR(2,K2,K1,K0)=ICODE(2)
              ELSE IF(NCODE(2).EQ.2) THEN
                IQFR(2,K2,K1,K0)=-2
              ELSE IF(NCODE(2).EQ.7) THEN
                IQFR(2,K2,K1,K0)=-3
              ENDIF
            ENDIF
*
            IF(K1.EQ.1) THEN
              LL1=.TRUE.
            ELSE
              LL1=(MAT(K2,K1-1,K0).EQ.0)
            ENDIF
            IF(LL1) THEN
              LALB=(NCODE(3).EQ.1).OR.(NCODE(3).EQ.6)
              IF(LALB.AND.(ICODE(3).EQ.0)) THEN
                QFR(3,K2,K1,K0)=ALB(ZCODE(3))
                IQFR(3,K2,K1,K0)=-1
              ELSE IF(LALB) THEN
                QFR(3,K2,K1,K0)=1.0
                IQFR(3,K2,K1,K0)=ICODE(3)
              ELSE IF(NCODE(3).EQ.2) THEN
                IQFR(3,K2,K1,K0)=-2
              ELSE IF(NCODE(3).EQ.7) THEN
                IQFR(3,K2,K1,K0)=-3
              ENDIF
            ENDIF
*
            IF(K1.EQ.LY) THEN
              LL1=.TRUE.
            ELSE
              LL1=(MAT(K2,K1+1,K0).EQ.0)
            ENDIF
            IF(LL1) THEN
              LALB=(NCODE(4).EQ.1).OR.(NCODE(4).EQ.6)
              IF(LALB.AND.(ICODE(4).EQ.0)) THEN
                QFR(4,K2,K1,K0)=ALB(ZCODE(4))
                IQFR(4,K2,K1,K0)=-1
              ELSE IF(LALB) THEN
                QFR(4,K2,K1,K0)=1.0
                IQFR(4,K2,K1,K0)=ICODE(4)
              ELSE IF(NCODE(4).EQ.2) THEN
                IQFR(4,K2,K1,K0)=-2
              ELSE IF(NCODE(4).EQ.7) THEN
                IQFR(4,K2,K1,K0)=-3
              ENDIF
            ENDIF
*
            IF(K0.EQ.1) THEN
              LL1=.TRUE.
            ELSE
              LL1=(MAT(K2,K1,K0-1).EQ.0)
            ENDIF
            IF(LL1) THEN
              LALB=(NCODE(5).EQ.1).OR.(NCODE(5).EQ.6)
              IF(LALB.AND.(ICODE(5).EQ.0)) THEN
                QFR(5,K2,K1,K0)=ALB(ZCODE(5))
                IQFR(5,K2,K1,K0)=-1
              ELSE IF(LALB) THEN
                QFR(5,K2,K1,K0)=1.0
                IQFR(5,K2,K1,K0)=ICODE(5)
              ELSE IF(NCODE(5).EQ.2) THEN
                IQFR(5,K2,K1,K0)=-2
              ELSE IF(NCODE(5).EQ.7) THEN
                IQFR(5,K2,K1,K0)=-3
              ENDIF
            ENDIF
*
            IF(K0.EQ.LZ) THEN
              LL1=.TRUE.
            ELSE
              LL1=(MAT(K2,K1,K0+1).EQ.0)
            ENDIF
            IF(LL1) THEN
              LALB=(NCODE(6).EQ.1).OR.(NCODE(6).EQ.6)
              IF(LALB.AND.(ICODE(6).EQ.0)) THEN
                QFR(6,K2,K1,K0)=ALB(ZCODE(6))
                IQFR(6,K2,K1,K0)=-1
              ELSE IF(LALB) THEN
                QFR(6,K2,K1,K0)=1.0
                IQFR(6,K2,K1,K0)=ICODE(6)
              ELSE IF(NCODE(6).EQ.2) THEN
                IQFR(6,K2,K1,K0)=-2
              ELSE IF(NCODE(6).EQ.7) THEN
                IQFR(6,K2,K1,K0)=-3
              ENDIF
            ENDIF
*----
*  TRAN BOUNDARY CONDITION
*----
            IF((K2.EQ.1).AND.(NCODE(1).EQ.4)) THEN
              KN(1,K2,K1,K0)=KN(2,LX,K1,K0)
            ENDIF
            IF((K2.EQ.LX).AND.(NCODE(2).EQ.4)) THEN
              KN(2,K2,K1,K0)=KN(1,1,K1,K0)
            ENDIF
            IF((K1.EQ.1).AND.(NCODE(3).EQ.4)) THEN
              KN(3,K2,K1,K0)=KN(2,K2,LY,K0)
            ENDIF
            IF((K1.EQ.LY).AND.(NCODE(4).EQ.4)) THEN
              KN(4,K2,K1,K0)=KN(1,K2,1,K0)
            ENDIF
            IF((K0.EQ.1).AND.(NCODE(5).EQ.4)) THEN
              KN(5,K2,K1,K0)=KN(6,K2,K1,LZ)
            ENDIF
            IF((K0.EQ.LZ).AND.(NCODE(6).EQ.4)) THEN
              KN(6,K2,K1,K0)=KN(5,K2,K1,1)
            ENDIF
*
            VOL(K2,K1,K0)=XX(K2,K1,K0)*YY(K2,K1,K0)*ZZ(K2,K1,K0)
          ENDDO
        ENDDO
      ENDDO
* END OF THE MAIN LOOP OVER NODES.
*
      IF(IMPX.GE.2) THEN
         WRITE(6,720) VOL(:LX,:LY,:LZ)
         WRITE(6,750)
         DO K0=1,LZ
           DO K1=1,LY
             DO K2=1,LX
               IF(MAT(K2,K1,K0).LE.0) CYCLE
               KEL=(K0-1)*LX*LY+(K1-1)*LX+K2
               WRITE (6,760) KEL,(KN(I,K2,K1,K0),I=1,6),
     1         (QFR(I,K2,K1,K0),I=1,6),(IQFR(I,K2,K1,K0),I=1,6)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      RETURN
*
  700 FORMAT(/46H NSSDFC: COARSE MESH FINITE DIFFERENCE METHOD.//3H NU,
     1 28HMBER OF NODES ALONG X AXIS =,I3/17X,14HALONG Y AXIS =,I3/
     2 17X,14HALONG Z AXIS =,I3)
  720 FORMAT(/17H VOLUMES PER NODE/(1X,1P,10E13.4))
  750 FORMAT(/22H NUMBERING OF UNKNOWNS/1X,21(1H-)//4X,4HNODE,5X,3HINT,
     1 26HERFACE NET CURRENT INDICES,28X,23HVOID BOUNDARY CONDITION)
  760 FORMAT(1X,I6,7X,6I8,6X,6F9.2/68X,6I9)
      END
