*DECK TRIHEX
      SUBROUTINE TRIHEX (IOPT,LX,LZ,LL4,MAT,KN,NCODE,IPTRK)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a mesh corner finite difference or
* Lagrangian finite element discretization of hexagonal geometry.
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
* IOPT    type of hexagonal Lagrangian finite element:
*         = 1 for hexagonal element with 6 points;
*         = 2 for hexagonal element with 7 points and triangular element
*         (IOPT.le.2) for Bivac;
*         = 3 for hexagonal element with 6 points;
*         = 4 for hexagonal element with 7 points and triangular element
*         (IOPT.gt.2) for Trivac.
* LX      number of elements in the XY plane.
* LZ      number of elements along Z axis.
* MAT     mixture index assigned to each element.
*
*Parameters: output
* LL4     total number of unknown (variational coefficients) per
*         energy group (order of system matrices).
* KN      element-ordered unknown list. Dimensionned to LC*LX*LZ
*         where LC= 7 for triangle/Bivac, 6 for hexagon/Bivac,
*         14 for triangle/Trivac and 12 for hexagon/Trivac.
* IPTRK   L_TRACK pointer to the tracking information.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IOPT,LX,LZ,LL4,MAT(LX*LZ),KN(*),NCODE(6)
      TYPE(C_PTR) IPTRK
*----
*  LOCAL VARIABLES
*----
      LOGICAL LC1,LMD,LMG,LED,LEG,LIV,CADI,LNC1,LNC2
      INTEGER IRO(2,7),ITAB(14),ICO(6,6)
      INTEGER, DIMENSION(:), ALLOCATABLE :: IISS,IPER,NWXY,IENV,IDX,IDY
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICG
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: NIK
      DATA ITAB /0,2*1,0,2*-1,0,1,0,-1,0,1,0,-1/
      DATA IRO  /2,3,3,7,7,6,4,4,1,2,5,1,6,5/
      DATA ICO  /2,1,5,6,7,3,1,5,6,7,3,2,3,2,1,5,6,7,2,1,5,6,7,3,7,3,2,
     >           1,5,6,3,2,1,5,6,7/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IISS(LX),IPER(LX),NWXY(LX),IENV(LX))
      ALLOCATE(ICG(14,LX),NIK(3,7,LX))
*
      IZI  = 0
      IZE  = 0
      IZS  = 0
      KEL  = 0
      IPAR = IOPT
      NC = INT((SQRT(REAL((4*LX-1)/3))+1.)/2.)
      IF(IPAR.LT.1.OR.IPAR.GT.4) CALL XABORT('TRIHEX : INVALID DATA.')
      M1 = 2 + 3*(NC-1)*(NC-2)
      IF(NC.EQ.1) M1=1
      DO 10 KX = 1,LX
         IISS(KX) = 0
         IENV(KX) = KX
         IF(MAT(KX).LE.0) GO TO 10
         KEL = KEL + 1
         IISS(KEL) = KX
  10  CONTINUE
      DO 20 IY = 1,14
         ICG(IY,1) = ITAB(IY)
  20  CONTINUE
      DO 45 J = 1,LX
      DO 40 KF = 1, 6
         N = NEIGHB(J,KF,9,LX,POIDS)
         IF(N.GT.LX) GOTO 40
         DO 30 I = 1,14
            IF(I.LE.7) THEN
               IF((KF.EQ.1).OR.(KF.EQ.5)) ICG(I,N)=ICG(I,J)+1
               IF((KF.EQ.2).OR.(KF.EQ.4)) ICG(I,N)=ICG(I,J)-1
               IF(KF.EQ.3)                ICG(I,N)=ICG(I,J)-2
               IF(KF.EQ.6)                ICG(I,N)=ICG(I,J)+2
            ELSE
               IF((KF.EQ.1).OR.(KF.EQ.3)) ICG(I,N)=ICG(I,J)+1
               IF((KF.EQ.4).OR.(KF.EQ.6)) ICG(I,N)=ICG(I,J)-1
               IF(KF.EQ.5)                ICG(I,N)=ICG(I,J)-2
               IF(KF.EQ.2)                ICG(I,N)=ICG(I,J)+2
            ENDIF
  30     CONTINUE
  40  CONTINUE
  45  CONTINUE
      ICHOIX = 1
      CADI = .FALSE.
      IF(IPAR.GT.2) THEN
         CADI = .TRUE.
         ICHOIX = 3
         IPAR = IPAR - 2
      ENDIF
       NELEMW = -1
       NELEMX = -1
       NELEMY = -1
       DO 190 IWXY = 1,ICHOIX
         DO 51 IY = 1,7
         DO 50 IX = 1,LX
            NIK(IWXY,IY,IX) = 0
  50     CONTINUE
  51     CONTINUE
         IF(NCODE(1).EQ.7) THEN
            DO 60 J  = 1,LX
               IF(MAT(J).LE.0) GO TO 60
               DO 55 KF = 1,6
                  N = NEIGHB(J,KF,9,LX,POIDS)
                  IF((N.GT.LX).OR.(MAT(N).LE.0)) THEN
                     NIK(IWXY,ICO(KF,2*IWXY-1),J) = -1
                     NIK(IWXY,ICO(KF,2*IWXY),J)   = -1
                  ENDIF
  55           CONTINUE
  60        CONTINUE
         ENDIF
         IF(IWXY.EQ.2.OR.IWXY.EQ.3) M1 = M1 + NC - 1
         CALL BIVPER(M1,IWXY,LX,LX,IPER,IENV)
         DO 65 I = 1,LX
            NWXY(IPER(I)) = I
  65     CONTINUE
         NELEM = 1
         DO 70 I = 1,NC
            IV = NWXY(I)
            IF(MAT(IV).LE.0) GO TO 70
            IF(NIK(IWXY,1,IV).EQ.-1) GO TO 70
            NIK(IWXY,1,IV) = NELEM
            NELEM = NELEM + 1
  70     CONTINUE
         DO 140 I = 1,LX
            LIV = (I.GT.1)
            IV1 = NWXY(I)
            IF(MAT(IV1).LE.0) THEN
               IFACE1 = IWXY + 5
               IFACE2 = IWXY + 2
               IF(IFACE1.GT.6) IFACE1 = IFACE1 - 6
               IF(IFACE2.GT.6) IFACE2 = IFACE2 - 6
               IND = NEIGHB(IV1,IFACE1,9,LX,P)
               ING = NEIGHB(IV1,IFACE2,9,LX,P)
               LED = (IND.GT.LX)
               LEG = (ING.GT.LX)
               LMD=.FALSE.
               LMG=.FALSE.
               IF(.NOT.LED) LMD = (MAT(IND).LE.0)
               IF(.NOT.LEG) LMG = (MAT(ING).LE.0)
               IF(LED.OR.LEG.OR.LMD.OR.LMG) THEN
                  IVK1 = -1
                  IVK2 = -1
                  IF((LMD.OR.LED).AND.(LMG.OR.LEG)) THEN
                     IVK1 = 1
                     IVK2 = 0
                  ELSE IF(LMD.OR.LED) THEN
                     IVK1 = 1
                     IVK2 = 3
                  ELSE IF(LMG.OR.LEG) THEN
                     IVK1 = 0
                     IVK2 = 0
                  ENDIF
                  IFACE3 = IWXY+IVK1+3
                  IFACE4 = IWXY+IVK2+2
                  IFACE5 = IWXY-IVK1+1
                  IF(IFACE3.GT.6) IFACE3 = IFACE3 - 6
                  IF(IFACE4.GT.6) IFACE4 = IFACE4 - 6
                  IF(IFACE5.GT.6) IFACE5 = IFACE5 - 6
                  IV3 = NEIGHB(IV1,IFACE3,9,LX,P)
                  IF((IV3.GT.LX).OR.(MAT(IV3).LE.0)) GOTO 90
  80              IF(NIK(IWXY,1,IV3).EQ.0)   THEN
                     NIK(IWXY,1,IV3) = NELEM
                     NELEM = NELEM + 1
                  ENDIF
                  IV3 = NEIGHB(IV3,IFACE4,9,LX,P)
                  IF((IV3.LE.LX).AND.(MAT(IV3).GT.0)) THEN
                     IV5 = NEIGHB(IV3,IFACE4-1,9,LX,P)
                     IF((IV5.GT.LX).OR.(MAT(IV5).GT.0)) GO TO 140
                     GO TO 80
                  ENDIF
  90              IV4 = NEIGHB(IV1,IFACE5,9,LX,P)
                  IF((IV4.GT.LX).OR.(MAT(IV4).LE.0)) GOTO 140
 100              IF(NIK(IWXY,7,IV4).EQ.0)   THEN
                     NIK(IWXY,7,IV4) = NELEM
                     NELEM = NELEM + 1
                  ENDIF
                  IV4 = NEIGHB(IV4,IFACE4,9,LX,P)
                  IF((IV4.LE.LX).AND.(MAT(IV4).GT.0)) GOTO 100
               ENDIF
            ELSE
               IF(IWXY.EQ.1) THEN
                  DO 110 K = 0,4
                     IF((IPAR.EQ.1).AND.(K.EQ.2)) GO TO 110
                     IF(LIV) THEN
                        IV2 = NWXY(I-1)
                        IF((ICG(K+2,IV1).EQ.ICG(5,IV2))
     >                     .AND.(MAT(IV2).GT.0)) THEN
                           NIK(1,2,IV1) = NIK(1,5,IV2)
                           NIK(1,3,IV1) = NIK(1,6,IV2)
                           GO TO 110
                        ENDIF
                     ENDIF
                     IF(NIK(1,K+2,IV1).EQ.-1) GO TO 110
                     NIK(1,K+2,IV1) = NELEM
                     NELEM = NELEM + 1
 110              CONTINUE
               ELSE IF(IWXY.EQ.2) THEN
                  DO 120 K = 0,4
                     IF((IPAR.EQ.1).AND.(K.EQ.2)) GO TO 120
                     IF(LIV) THEN
                        IV2 = NWXY(I-1)
                        IF((ICG(K+1,IV1).EQ.ICG(6,IV2))
     >                     .AND.(ICG(K+8,IV1).EQ.ICG(13,IV2))
     >                     .AND.(MAT(IV2).GT.0)) THEN
                           NIK(2,2,IV1) = NIK(2,5,IV2)
                           GO TO 120
                        ELSE IF((ICG(K+1,IV1).EQ.ICG(7,IV2))
     >                         .AND.(ICG(K+8,IV1).EQ.ICG(14,IV2))
     >                         .AND.(MAT(IV2).GT.0)) THEN
                           NIK(2,3,IV1) = NIK(2,6,IV2)
                           GO TO 120
                        ENDIF
                     ENDIF
                     IF(NIK(2,K+2,IV1).EQ.-1) GO TO 120
                     NIK(2,K+2,IV1) = NELEM
                     NELEM = NELEM + 1
 120              CONTINUE
               ELSE IF(IWXY.EQ.3) THEN
                  LK = 1
                  DO 130 K = 5,1,-1
                     LK = LK + 1
                     IF((IPAR.EQ.1).AND.(LK.EQ.4)) GO TO 130
                     IF(LIV) THEN
                        IV2 = NWXY(I-1)
                        IF((ICG(K,IV1).EQ.ICG(7,IV2))
     >                     .AND.(ICG(K+7,IV1).EQ.ICG(14,IV2))
     >                     .AND.(MAT(IV2).GT.0)) THEN
                           NIK(3,2,IV1) = NIK(3,5,IV2)
                           GO TO 130
                        ELSE IF((ICG(K-3,IV1).EQ.ICG(3,IV2))
     >                          .AND.(ICG(K+4,IV1).EQ.ICG(10,IV2))
     >                          .AND.(MAT(IV2).GT.0)) THEN
                           NIK(3,3,IV1) = NIK(3,6,IV2)
                           GO TO 130
                        ENDIF
                     ENDIF
                     IF(NIK(3,LK,IV1).EQ.-1) GO TO 130
                     NIK(3,LK,IV1) = NELEM
                     NELEM = NELEM + 1
 130              CONTINUE
               ENDIF
            ENDIF
 140     CONTINUE
         DO 160 I = 1,LX
            IV = NWXY(I)
            IF(MAT(IV).LE.0) GO TO 160
            DO 150 K = 1,2
               IFACE = K+IWXY-1
               IF(IFACE.GT.6) IFACE = IFACE - 6
               INE = NEIGHB(IV,IFACE,9,LX,P)
               IF((INE.GT.LX).OR.(MAT(INE).LE.0)) GO TO 150
               IF(K.EQ.1) NIK(IWXY,1,IV) = NIK(IWXY,6,INE)
               IF(K.EQ.2) NIK(IWXY,1,IV) = NIK(IWXY,3,INE)
 150        CONTINUE
 160     CONTINUE
         DO 180 I = 1,LX
            IV = NWXY(I)
            IF(MAT(IV).LE.0) GO TO 180
            LC1 = .FALSE.
            LIV = .TRUE.
            DO 170 K = 4,5
               IFACE = K+IWXY-1
               IF(IFACE.GT.6) IFACE = IFACE - 6
               INE = NEIGHB(IV,IFACE,9,LX,P)
               IF((INE.GT.LX).OR.(MAT(INE).LE.0)) THEN
                  IF(K.EQ.4) LC1 = .TRUE.
                  IF(NIK(IWXY,7,IV).EQ.-1) GO TO 170
                  IF(LC1.AND.K.EQ.5.AND.(NIK(IWXY,7,IV).EQ.0)) THEN
                     NIK(IWXY,7,IV) = NELEM
                     NELEM = NELEM + 1
                  ENDIF
               ELSE
                  IF(K.EQ.4) THEN
                     LIV = .FALSE.
                     NIK(IWXY,7,IV) = NIK(IWXY,2,INE)
                  ENDIF
                  IF(K.EQ.5.AND.LIV) NIK(IWXY,7,IV)=NIK(IWXY,5,INE)
               ENDIF
 170        CONTINUE
 180     CONTINUE
         IF(IWXY.EQ.1) NELEMW = NELEM - 1
         IF(IWXY.EQ.2) NELEMX = NELEM - 1
         IF(IWXY.EQ.3) NELEMY = NELEM - 1
 190  CONTINUE
      MEL = 0
      NPT = NELEM -1
      IF(ICHOIX.EQ.3) THEN
         IF((NELEMW.NE.NELEMX).OR.(NELEMW.NE.NELEMY)) THEN
            CALL XABORT('TRIHEX: ECHEC DE LA NUMEROTATION.')
         ENDIF
         ALLOCATE(IDX(NPT),IDY(NPT))
         DO 200 I = 1,LX
            IF(MAT(I).LE.0) GO TO 200
            DO 195 J = 1,7
               IF((IPAR.EQ.1).AND.(J.EQ.4)) GO TO 195
               IF(NIK(1,J,I).EQ.-1) GO TO 195
               IDX(NIK(1,J,I)) = NIK(2,IRO(1,J),I)
               IDY(NIK(1,J,I)) = NIK(3,IRO(2,J),I)
 195        CONTINUE
 200     CONTINUE
         CALL LCMPUT(IPTRK,'ILX',NPT,1,IDX)
         CALL LCMPUT(IPTRK,'ILY',NPT,1,IDY)
         DEALLOCATE(IDY,IDX)
      ENDIF
      NINC =  6
      IF(CADI) NINC = 12
      IF(IPAR.EQ.2) NINC =   7
      IF((IPAR.EQ.2).AND.CADI) NINC = 14
      LNC1 = .FALSE.
      LNC2 = .FALSE.
      KNZ = 0
      IF(LZ.GT.1)   KNZ = NPT
      IF(CADI.AND.(LZ.GT.1)) CALL LCMPUT(IPTRK,'NCODE',6,1,NCODE)
      IF((LZ.GT.1).AND.((NCODE(5).EQ.7).OR.(NCODE(6).EQ.7))) THEN
         LNC1 = (NCODE(5).EQ.7)
         LNC2 = (NCODE(6).EQ.7)
         IF(LNC1) IZI = 1
         IF(LNC2) IZS = 1
         IF(LNC1.AND.LNC2) THEN
            IZE = 1
            IZS = 1
            IZI = 0
         ENDIF
      ENDIF
      DO 230 MZ = 1, LZ
         KEL = 0
         DO 220 MX = 1, LX
            IF(MAT(MX).LE.0) GO TO 220
            KEL = KEL + 1
            ML = IISS(KEL)
            IJX = 1
            DO 210 JX = 1,7
               KN(MEL+IJX) = 0
               IF(CADI) KN(MEL+IJX+NINC/2) = 0
               IF((IPAR.EQ.1).AND.(JX.EQ.4)) GO TO 210
               IF(NIK(1,JX,ML).EQ.-1) GO TO 205
               IF(MZ.EQ.1.AND.LNC1) THEN
                  KN(MEL+IJX+NINC/2) = NIK(1,JX,ML)
                  GO TO 205
               ELSE IF(MZ.EQ.LZ.AND.LNC2) THEN
                  KN(MEL+IJX) = NIK(1,JX,ML)+(MZ-IZS-IZE)*KNZ
                  GO TO 205
               ENDIF
               KN(MEL+IJX) = NIK(1,JX,ML)+(MZ-1-IZE-IZI)*KNZ
               IF(CADI) KN(MEL+IJX+NINC/2)=NIK(1,JX,ML)+(MZ-IZE-IZI)*KNZ
 205           IJX = IJX + 1
 210        CONTINUE
            MEL = MEL + NINC
 220     CONTINUE
 230  CONTINUE
      LL4 = NPT
      IF(LZ.GT.1) THEN
         LL4 = NPT*(LZ+1)
         IF(LNC1.OR.LNC2)   LL4 = NPT*LZ
         IF(LNC1.AND.LNC2) LL4 = NPT*(LZ-1)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(NIK,ICG,IENV,NWXY,IPER,IISS)
      RETURN
      END
