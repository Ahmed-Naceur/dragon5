*DECK TRINTR
      SUBROUTINE TRINTR (ISPLH,IPTRK,LX,LI4,IHEX,MAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a mesh centred finite difference for
* hexagonal geometry (each hexagon represented by 6 triangles).
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
* ISPLH  used to compute the number of triangles per hexagon
*        (6*(ISPLH-1)**2).
* IPTRK  L_TRACK pointer to the tracking information.
* LX     number of elements.
* IHEX   type of hexagonal boundary condition.
* MAT    mixture index assigned to each element.
*
*Parameters: output
* LI4    total number of unknown (variational coefficients) per
*        energy group per plan.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER ISPLH,LX,LI4,IHEX,MAT(LX)
*----
*  LOCAL VARIABLES
*----
      LOGICAL LPAIR
      INTEGER IRO(180,2),NBL(20,2)
      INTEGER, DIMENSION(:), ALLOCATABLE :: IW,IY,IPO,IXN,IDX,IDY
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IC1,IC2
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: NIK
      DATA NBL  /3,3,5,7,7,5,7,9,11,11,9,7,9,11,13,15,15,13,11,9,
     >           3,3,7,5,5,7,11,9,7,7,9,11,15,13,11,9,9,11,13,15/
      DATA IRO  /4,1,2,5,6,3, 13,6,7,1,2,20,14,15,8,9,3,4,21,22,16,
     >  17,10,11,5,23,24,18,19,12, 28,17,18,8,9,1,2,39,29,30,19,20,
     >  10,11,3,4,48,40,41,31,32,21,22,12,13,5,6,49,50,42,43,33,34,
     >  23,24,14,15,7,51,52,44,45,35,36,25,26,16,53,54,46,47,37,38,
     >  27, 49,34,35,21,22,10,11,1,2,64,50,51,36,37,23,24,12,13,3,
     >  4,77,65,66,52,53,38,39,25,26,14,15,5,6,88,78,79,67,68,54,55,
     >  40,41,27,28,16,17,7,8,89,90,80,81,69,70,56,57,42,43,29,30,18,
     >  19,9,91,92,82,83,71,72,58,59,44,45,31,32,20,93,94,84,85,73,
     >  74,60,61,46,47,33,95,96,86,87,75,76,62,63,48,
     >  5,4,1,6,3,2, 21,20,14,13,6,23,22,16,15,8,7,1,24,18,17,10,9,
     >  3,2,19,12,11,5,4, 49,48,40,39,29,28,17,51,50,42,41,31,30,19,
     >  18,8,53,52,44,43,33,32,21,20,10,9,1,54,46,45,35,34,23,22,12,
     >  11,3,2,47,37,36,25,24,14,13,5,4,38,27,26,16,15,7,6,
     >  89,88,78,77,65,64,50,49,34,91,90,80,79,67,66,52,51,36,35,21,
     >  93,92,82,81,69,68,54,53,38,37,23,22,10,95,94,84,83,71,70,56,
     >  55,40,39,25,24,12,11,1,96,86,85,73,72,58,57,42,41,27,26,14,
     >  13,3,2,87,75,74,60,59,44,43,29,28,16,15,5,4,76,62,61,46,45,
     >  31,30,18,17,7,6,63,48,47,33,32,20,19,9,8/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IC1(3,2*LX),IC2(3,2*LX*(ISPLH-1)))
      ALLOCATE(NIK(3,6*(ISPLH-1)**2,LX))
*
      NBE = 0
      NC = INT((SQRT(REAL((4*LX-1)/3))+1.)/2.)
      L1 = 3*NC - 2
      COURS2 = REAL(L1)/2.
      LPAIR = (AINT(COURS2).EQ.COURS2)
      IF(ISPLH.LE.3) ISAU =   2*(ISPLH-2)
      IF(ISPLH.GE.4) ISAU =   6*(ISPLH-3)
      ALLOCATE(IW(3*L1),IY(L1))
      IW(1) = 2+3*(NC-1)*(NC-2)
      DO 10 I = 1,L1
         IF(I.LT.L1)     IW(I+1)    = 2+3*NC*(NC-1)-I
         IF(I.LE.NC)     IW(I+L1)   = 3+(3*NC-5)*(NC-1)-I
         IF(I.GT.NC)     IW(I+L1)   = 2+3*NC*(NC-1)-I+NC
         IF(I.LE.2*NC-1) IW(I+2*L1) = 3+(3*NC-4)*(NC-1)-I
         IF(I.GT.2*NC-1) IW(I+2*L1) = 2+3*NC*(NC-1)-I+2*NC-1
  10  CONTINUE
      IF(LPAIR) THEN
         DO 20 I = 1,L1/2
            IF(I.LE.NC) IY(I) = 1+2*(I-1)
            IF(I.GT.NC) IY(I) = IY(NC)
  20     CONTINUE
         KEL = 1
         DO 30 I = L1,L1/2,-1
            IF(I.GE.(L1-NC-1)) IY(I) = IY(KEL)
            IF(I.LT.(L1-NC-1)) IY(I) = IY(NC)
            KEL = KEL + 1
  30     CONTINUE
      ELSE
         DO 40 I = 1,(L1+1)/2
            IF(I.LE.NC) IY(I) = 1+2*(I-1)
            IF(I.GT.NC) IY(I) = IY(NC)
  40     CONTINUE
         KEL = 1
         DO 50 I = L1,(L1-1)/2,-1
            IF(I.GE.(L1-NC-1)) IY(I) = IY(KEL)
            IF(I.LT.(L1-NC-1)) IY(I) = IY(NC)
            KEL = KEL + 1
  50     CONTINUE
      ENDIF
      ICAS = 3
      DO 90 K = 1,ICAS
         KEL = 1
         DO 80 I = 1,L1
            IPAR = IW(I+(K-1)*L1)
            NPAR = IPAR
            IC1(K,KEL) = NPAR
            KEL = KEL + 1
            IF(I.GT.(2*NC-1)) GO TO 70
  60        NPAR = ABS(NEIGHB(NPAR,K+1,IHEX,LX,P))
            IF(NPAR.GT.LX) THEN
               IF(I.LT.NC.OR.I.GT.(2*NC-1)) GO TO 80
               IF(I.GE.NC.AND.I.LE.(2*NC-1)) NPAR = IPAR
            ENDIF
            IC1(K,KEL) = NPAR
            KEL = KEL + 1
  70        NPAR = ABS(NEIGHB(NPAR,K+2,IHEX,LX,P))
            IF(NPAR.GT.LX) GO TO 80
            IC1(K,KEL) = NPAR
            KEL = KEL + 1
            GO TO 60
  80     CONTINUE
  90  CONTINUE
      DO 140 K=1,ICAS
         IF(ISPLH.EQ.2) THEN
            DO 100 JX = 1,2*LX
               IC2(K,JX) = IC1(K,JX)
 100        CONTINUE
         ELSE
            JEL = 1
            IEL = 1
            KEL = 1
            MEL = 0
 110        IF(IEL.LE.2*LX) THEN
               IF(IC1(K,IEL).EQ.MEL) NBE = IY(KEL-1)
               IF(IC1(K,IEL).EQ.IW(KEL+(K-1)*L1)) THEN
                  NBE = IY(KEL)
                  KEL = KEL + 1
               ENDIF
               MEL = IC1(K,IEL)
               IFOIS = 0
               ISAUV = IEL
 120           DO 130 LDB = 1,NBE
                  IC2(K,JEL) = IC1(K,IEL)
                  JEL = JEL + 1
                  IEL = IEL + 1
 130           CONTINUE
               IFOIS = IFOIS + 1
               IF(IFOIS.LT.(ISPLH-1)) THEN
                  IEL = ISAUV
                  GO TO 120
               ENDIF
               GO TO 110
            ENDIF
         ENDIF
 140  CONTINUE
      DO 152 K=1,ICAS
      DO 151 I=1,LX
      DO 150 J=1,6*(ISPLH-1)**2
         NIK(K,J,I) = 0
 150  CONTINUE
 151  CONTINUE
 152  CONTINUE
      ALLOCATE(IPO(LX))
      DO 200 K=1,ICAS
         DO 160 KK=1,LX
            IPO(KK) = 1
 160     CONTINUE
         IA  = 1
         IX  = 1
         ILI = 1
         ICOMPT = 1
 170     IEL = 1
         JCL = 1
         IVAL = IC2(K,IX)
 180     IF(MAT(IC2(K,IX)).EQ.0) THEN
            IX  = IX  + 1
            IEL = IEL + 1
            JCL = JCL + 1
            IF(JCL.GT.2) JCL = 1
         ELSE
            IF(ILI+ISAU.GT.20) CALL XABORT('TRINTR: NBL OVERFLOW.')
            IDEB = IPO(IC2(K,IX))
            IFIN = IPO(IC2(K,IX)) + NBL(ILI+ISAU,JCL) - 1
            DO 190 J=IDEB,IFIN
               NIK(K,J,IC2(K,IX)) = ICOMPT
               ICOMPT = ICOMPT + 1
 190        CONTINUE
            IPO(IC2(K,IX)) = J
            IX  = IX  + 1
            IEL = IEL + 1
            JCL = JCL + 1
            IF(JCL.GT.2) JCL = 1
         ENDIF
         IF(IEL.LE.IY(IA)) GO TO   180
         IF(IX.GT.2*LX*(ISPLH-1))   GO TO  200
         IF(IC2(K,IX).NE.IVAL)   IA  = IA  + 1
         ILI = ILI + 1
         IF(ILI.LE.(ISPLH-1)) GO TO   170
         IF((ILI.GT.(ISPLH-1).AND.ILI.LE.2*(ISPLH-1)).AND.
     >       (IVAL.EQ.IC2(K,IX)))  GO TO  170
         ILI = 1
         IF(IA.GT.(1+2*(NC-1))) ILI = ISPLH
         IF(IX.LE.2*LX*(ISPLH-1)) GO TO   170
 200  CONTINUE
      LI4 = ICOMPT - 1
      DEALLOCATE(IPO,IY,IW)
      ALLOCATE(IXN(LI4),IDX(LI4),IDY(LI4))
      KEL = 0
      ICR = ISAU*(1+2*(ISPLH-2))
      DO 220 I = 1, LX
         IF(MAT(I).EQ.0) GO TO 220
         DO 210 J=1,6*(ISPLH-1)**2
            KEL = KEL + 1
            IXN(KEL) = NIK(1,J,I)
            IDX(NIK(1,J,I))=NIK(2,IRO(J+ICR,1),I)
            IDY(NIK(1,J,I))=NIK(3,IRO(J+ICR,2),I)
 210     CONTINUE
 220  CONTINUE
*
      CALL LCMPUT(IPTRK,'IKN',LI4,1,IXN)
      CALL LCMPUT(IPTRK,'ILX',LI4,1,IDX)
      CALL LCMPUT(IPTRK,'ILY',LI4,1,IDY)
      DEALLOCATE(IDY,IDX,IXN)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(NIK,IC2,IC1)
      RETURN
      END
