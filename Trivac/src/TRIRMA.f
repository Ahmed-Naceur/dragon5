*DECK TRIRMA
      SUBROUTINE TRIRMA(ISPLH,R,Q,RH,QH,RT,QT,LL,LC,ISR,QTHP,QTHZ,RTHG,
     > HW,HX,HY,HZ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the unit matrices for a mesh corner finite difference
* discretization in hexagonal geometry.
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
* ISPLH   hexagonal mesh-splitting flag:
*         =1 for complete hexagons; >1 for triangular elements.
* R       unit matrix.
* Q       unit matrix.
* RH      unit matrix.
* QH      unit matrix.
* RT      unit matrix.
* QT      unit matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE PARAMETERS
*----
      INTEGER ISPLH
      REAL R(2,2),Q(2,2),RH(6,6),QH(6,6),RT(3,3),QT(3,3)
      DOUBLE PRECISION QTHP(14,14),QTHZ(14,14),RTHG(14,14),
     >           HW(14,14),HX(14,14),HY(14,14),HZ(14,14)
*----
*  LOCAL VARIABLES
*----
      INTEGER IJ1(14),IJ2(14),ISR(8,25),ILIEN(6,3),IJ27(14),ISRH(8,6),
     > ISRT(8,7),IJ16(12),IJ26(12),IJ17(14)
      REAL HL(2,2),RFAC(28,7),RH2(7,7),QH2(7,7),RF6(24,6),RF7(28,7)
      DATA HL / 1.0,2*0.0,1.0/
      DATA ILIEN/6*4,2,1,5,6,7,3,1,5,6,7,3,2/
      DATA IJ16,IJ26 /1,2,3,4,5,6,1,2,3,4,5,6,6*1,6*2/
      DATA IJ17,IJ27 /1,2,3,4,5,6,7,1,2,3,4,5,6,7,7*1,7*2/
      DATA ISRT/2,1,5,6,7,3,1,8,1,5,6,7,3,2,2,9,9,8,12,13,14,10,3,10,
     >          8,12,13,14,10,9,4,11,6*0,5,12,6*0,6,13,6*0,7,14/
      DATA ISRH/2,1,4,5,6,3,1,7,1,4,5,6,3,2,2,8,8,7,10,11,12,9,3,9,
     >          7,10,11,12,9,8,4,10,6*0,5,11,6*0,6,12/
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
*----
*  COMPUTE THE HEXAGONAL MASS (RH2) AND STIFFNESS (QH2) MATRICES
*----
      IF(ISPLH.EQ.1) THEN
         LC=6
         DO 11 I=1,8
           DO 10 J=1,6
             ISR(I,J)=ISRH(I,J)
   10      CONTINUE
   11    CONTINUE
         DO 20 I=1,2*LC
           IJ1(I)=IJ16(I)
           IJ2(I)=IJ26(I)
   20    CONTINUE
         DO 31 I=1,4*LC
           DO 30 J=1,LC
             RFAC(I,J)=RF6(I,J)
   30      CONTINUE
   31    CONTINUE
         DO 41 I=1,LC
           DO 40 J=1,LC
             RH2(I,J)=RH(I,J)
             QH2(I,J)=QH(I,J)
   40      CONTINUE
   41    CONTINUE
      ELSE
         LC=7
         DO 51 I=1,8
           DO 50 J=1,7
             ISR(I,J)=ISRT(I,J)
   50      CONTINUE
   51    CONTINUE
         DO 60 I=1,2*LC
           IJ1(I)=IJ17(I)
           IJ2(I)=IJ27(I)
   60    CONTINUE
         DO 71 I=1,4*LC
           DO 70 J=1,LC
             RFAC(I,J)=RF7(I,J)
   70      CONTINUE
   71    CONTINUE
         DO 76 I=1,LC
           DO 75 J=1,LC
             RH2(I,J)=0.0
             QH2(I,J)=0.0
   75      CONTINUE
   76    CONTINUE
         DO 82 K=1,6
           DO 81 I=1,3
             NUMI=ILIEN(K,I)
             DO 80 J=1,3
               NUMJ=ILIEN(K,J)
               RH2(NUMI,NUMJ)=RH2(NUMI,NUMJ)+RT(I,J)
               QH2(NUMI,NUMJ)=QH2(NUMI,NUMJ)+QT(I,J)
   80        CONTINUE
   81      CONTINUE
   82    CONTINUE
      ENDIF
      LL=2*LC
      DO 91 I=1,LL
        I1=IJ1(I)
        I2=IJ2(I)
        DO 90 J=1,LL
          J1=IJ1(J)
          J2=IJ2(J)
          HW(I,J)  =RFAC(I1     ,J1) * HL(I2,J2)
          HX(I,J)  =RFAC(I1+LC  ,J1) * HL(I2,J2)
          HY(I,J)  =RFAC(I1+2*LC,J1) * HL(I2,J2)
          HZ(I,J)  =RFAC(I1+3*LC,J1)
          RTHG(I,J)=RH2(I1,J1) * R(I2,J2)
          QTHP(I,J)=QH2(I1,J1) * R(I2,J2)
          QTHZ(I,J)=RH2(I1,J1) * Q(I2,J2)
   90   CONTINUE
   91 CONTINUE
      RETURN
      END
