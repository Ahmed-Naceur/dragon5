*DECK SYBRTK
      SUBROUTINE SYBRTK (NA,NX,NREG,A,B,RAYRE,ILIGN,INORM,IRECT,IQW,LR,
     1 Z,LI,IZ,PREC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the tracking information related to a square or rectangular
* heterogeneous cell.
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
* NA      number of angles in (0,$\\pi$/4).
* NX      number of tracks in each sub domain for a given angle.
* NREG    number of regions in the cell.
* A       Cartesian dimension of the cell along the X axis.
* B       Cartesian dimension of the cell along the Y axis.
* RAYRE   radius of each cylinder (RAYRE(1)=0.0).
* ILIGN   tracking print flag (=1 to print the tracking).
* INORM   track normalization flag (=1 to avoid track normalization).
* IRECT   rectangular flag (=1 to avoid considering symmetries of
*         square cells).
* IQW     equal weight quadrature flag (=1 to use equal weight
*         quadratures in angle and space).
*
*Parameters: output
* LR      exact size of array Z.
*         L.LE.4+2*NA*(9+2*(NREG+1)*NX*NREG) for a square cell;
*         L.LE.4+2*NA*(9+2*(2*NREG+1)*NX*NREG) for a rectangular cell.
* Z       real tracking information.
*         Z(1) to Z(4) contain the numerical orthonormalization
*         factors.
* LI      size of array IZ.
*         L.LE.NREG+4+2*NA*(2+(NREG+1)*(3+2*NREG)) for a square cell;
*         L.LE.NREG+4+2*NA*(2+(2*NREG+1)*(3+2*NREG)) for a rectangular
*         cell.
* IZ      integer tracking information.
*         IZ(1)=5 and IZ(2)=NREG+1 for a square or rectangular cell;
*         IZ(3)=1 if the cell is rectangular or if IRECT=1;
*         IZ(3)=2 if the cell is square.
* PREC    accuracy obtained if the non-normalized tracks are used
*         to integrate the volumes.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NA,NX,NREG,ILIGN,INORM,IRECT,IQW,LR,LI,IZ(*)
      REAL A,B,RAYRE(NREG),Z(*),PREC
*----
*  LOCAL VARIABLES
*----
      PARAMETER (PIO2=1.570796327,PI=3.14159265358979)
      REAL ZX(64),WX(64),ZA(64),WA(64),ZXJ(64),WXJ(64)
      REAL, ALLOCATABLE, DIMENSION(:) :: VAP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(VAP(NREG))
*
      NA2=2*NA
      IF(NX.GT.10) CALL XABORT('SYBRTK: NX IS GREATER THAN 10.')
      IF(NA.GT.64) CALL XABORT('SYBRTK: NA IS GREATER THAN 64.')
      IF(2.0*RAYRE(NREG).GT.MIN(A,B)) CALL XABORT('SYBRTK: A RADIUS IS'
     1 //' GREATER THAN HALF A SIDE.')
*----
*  COMPUTE VOLUMES
*----
      VOL=A*B
      DO 10 IR=NREG,1,-1
      R2=PI*RAYRE(IR)**2
      Z(IR)=VOL-R2
      VOL=R2
   10 CONTINUE
*     
      IF(IQW.EQ.0) THEN
*        GAUSS-LEGENDRE AND GAUSS-JACOBI INTEGRATION POINTS.
         CALL ALGPT(NX,-1.,1.,ZX,WX)
         CALL ALGJP(NX,ZXJ,WXJ)
         CALL ALGPT(NA,-1.,0.,ZA,WA)
         CALL ALGPT(NA,0.,1.,ZA(NA+1),WA(NA+1))
      ELSE
*        EQUAL WEIGHT INTEGRATION POINTS.
         DO 15 I=1,NX
         ZX(I)=(2.0*REAL(I)-1.0)/REAL(NX)-1.0
         WX(I)=2.0/REAL(NX)
         ZXJ(I)=0.5*(2.0*REAL(I)-1.0)/REAL(NX)
         WXJ(I)=ZXJ(I)/REAL(NX)
   15    CONTINUE
         DO 20 I=1,NA2
         ZA(I)=(2.0*REAL(I)-1.0)/REAL(NA2)-1.0
         WA(I)=2.0/REAL(NA2)
   20    CONTINUE
      ENDIF
      IZ(1)=5
      IZ(2)=NREG+1
      IZ(3)=2
      IF((A.NE.B).OR.(IRECT.EQ.1)) IZ(3)=1
      IZ(4)=NA2
      AI=1.0/A
      BI=1.0/B
      AO2=0.5*A
      BO2=0.5*B
      AB=A*B
      PREC=0.0
      LI=4
      LR=NREG+4
*----
*  INTEGRATION IN ANGLE FROM 0 TO PI/2
*----
      ZN1=0.0
      ZN2=0.0
      ZN3=0.0
      DO 350 IA=1,NA2
      PHI=0.5*PIO2*(ZA(IA)+1.0)
      SI=SIN(PHI)
      CO=COS(PHI)
      TA=SI/CO
      ZN1=ZN1+SI*WA(IA)
      ZN2=ZN2+SI*SI*WA(IA)
      ZN3=ZN3+SI*SI*SI*WA(IA)
      Z(LR+1)=SI
      Z(LR+2)=CO
      Z(LR+3)=SI
      Z(LR+4)=CO
      Z(LR+5)=CO
      Z(LR+6)=SI
      Z(LR+7)=CO
      Z(LR+8)=SI
      Z(LR+9)=WA(IA)
      LR=LR+9
*----
*  FIRST ANGULAR DOMAIN
*----
      L4=LI+1
      IZ(LI+1)=0
      IZ(LI+2)=0
      LI=LI+2
      IF((IZ(3).EQ.2).AND.(PHI.GT.0.5*PIO2)) GO TO 120
      X1=0.0
      XLIM=MIN(A,B/TA)
      DLIM=BO2*CO+(AO2-XLIM)*SI
      DO 100 K0=NREG,1,-1
      KMAX=NREG-K0+1
      X2=MIN(XLIM,XLIM-(RAYRE(K0)-DLIM)/SI)
      L3=LR+1
      L5=LI+1
      LI=LI+3
      CALL XDRSET(VAP,NREG,0.0)
      DO 50 IX=1,NX
      IF(K0.EQ.NREG) THEN
         S=0.5*(X2-X1)*SI*WX(IX)
         X=X1+0.5*(X2-X1)*(1.0+ZX(IX))
      ELSE
*        FLURIG CHANGE OF VARIABLE.
         S=2.0*(X2-X1)*SI*WXJ(IX)
         X=X1+(X2-X1)*ZXJ(IX)**2
      ENDIF
      Z(LR+1)=S*WA(IA)
      LR=LR+1
      C=BO2*SI-(AO2-X)*CO
      D=BO2*CO+(AO2-X)*SI
      D=D*D
      SUM=0.0
      CORDE=0.0
      DO 30 K=NREG,K0+1,-1
      RR=RAYRE(K)**2-D
      CORDE=SQRT(RR)
      DEL=C-CORDE
      SUM=SUM+DEL
      Z(LR+NREG-K+1)=DEL
      VAP(K)=VAP(K)+DEL*S
      C=CORDE
   30 CONTINUE
      IF(KMAX.NE.1) THEN
         DEL=2.0*CORDE
         SUM=SUM+DEL
         Z(LR+KMAX)=DEL
         VAP(K)=VAP(K)+DEL*S
         DO 40 I=1,KMAX-2
         DEL=Z(LR+KMAX-I)
         SUM=SUM+DEL
         Z(LR+KMAX+I)=DEL
         VAP(K+I)=VAP(K+I)+DEL*S
   40    CONTINUE
      ENDIF
      LR=LR+2*KMAX-1
      DEL=X/CO-SUM
      Z(LR)=DEL
      VAP(NREG)=VAP(NREG)+DEL*S
   50 CONTINUE
      DO 60 K=KMAX-1,1-KMAX,-1
      IZ(LI+K+KMAX)=3+ABS(K)+1+NREG-KMAX
   60 CONTINUE
      LI=LI+2*KMAX
      IZ(L5)=2*KMAX-1
      IZ(L5+1)=NX
      IZ(L5+2)=2 ! ISURF
      IZ(LI)=1   ! JSURF
*----
*  VOLUME NORMALIZATION
*----
      IF((INORM.EQ.0).AND.(K0.LT.NREG)) THEN
         DLIM1=BO2*CO+(AO2-X2)*SI
         DLIM2=BO2*CO+(AO2-X1)*SI
         VW1=0.0
         SUMVAP=0.0
         DO 70 I=K0,NREG-1
         SUMVAP=SUMVAP+VAP(I)
         RW=RAYRE(I+1)
         VEX1=RW*RW*ACOS(DLIM1/RW)-DLIM1*SQRT(RW*RW-DLIM1*DLIM1)
         IF(RW.GT.DLIM2)
     1   VEX1=VEX1-(RW*RW*ACOS(DLIM2/RW)-DLIM2*SQRT(RW*RW-DLIM2*DLIM2))
         SUM=(VEX1-VW1)/VAP(I)
         PREC=MAX(PREC,ABS(1.0-SUM)*(VEX1-VW1)/AB)
         VW1=VEX1
         VAP(I)=SUM
   70    CONTINUE
         VEX1=0.5*(B*SI-(A-X1-X2)*CO)*(X2-X1)*SI
         VEX2=0.5*TA*(X2*X2-X1*X1)-VEX1
         SUM=(VEX1-0.5*VW1)/(VEX1-0.5*SUMVAP)
         PREC=MAX(PREC,ABS(1.0-SUM)*(VEX1-0.5*VW1)/AB)
         VEX1=SUM
         SUM=(VEX2-0.5*VW1)/(VEX2-0.5*SUMVAP)
         PREC=MAX(PREC,ABS(1.0-SUM)*(VEX2-0.5*VW1)/AB)
         VEX2=SUM
         DO 90 IX=1,NX
         KMAX=(IZ(L5)+1)/2
         Z(L3+KMAX)=Z(L3+KMAX)*VAP(K0)
         DO 80 I=1,KMAX-2
         Z(L3+KMAX-I)=Z(L3+KMAX-I)*VAP(K0+I)
         Z(L3+KMAX+I)=Z(L3+KMAX+I)*VAP(K0+I)
   80    CONTINUE
         Z(L3+1)=Z(L3+1)*VEX1
         Z(L3+2*KMAX-1)=Z(L3+2*KMAX-1)*VEX2
         L3=L3+2*KMAX
   90    CONTINUE
      ENDIF
      IZ(L4)=IZ(L4)+1
      IF(X2.GE.XLIM) GO TO 120
      X1=X2
  100 CONTINUE
*----
*  SECOND ANGULAR DOMAIN
*----
  120 IF(PHI.LE.ATAN(B*AI)) GO TO 240
      X1=B/TA
      XLIM=0.5*(A+X1)
      DO 230 K0=NREG,1,-1
      KMAX=NREG-K0+1
      X2=XLIM-RAYRE(K0)/SI
      IF(X2.LE.X1) GO TO 230
      L3=LR+1
      L5=LI+1
      LI=LI+3
      CALL XDRSET(VAP,NREG,0.0)
      DO 150 IX=1,NX
      IF(K0.EQ.NREG) THEN
         S=0.5*(X2-X1)*SI*WX(IX)
         X=X1+0.5*(X2-X1)*(1.0+ZX(IX))
      ELSE
*        FLURIG CHANGE OF VARIABLE.
         S=2.0*(X2-X1)*SI*WXJ(IX)
         X=X1+(X2-X1)*ZXJ(IX)**2
      ENDIF
      Z(LR+1)=S*WA(IA)
      LR=LR+1
      C=BO2*SI-(AO2-X)*CO
      D=BO2*CO+(AO2-X)*SI
      D=D*D
      SUM=0.0
      CORDE=0.0
      DO 130 K=NREG,K0+1,-1
      RR=RAYRE(K)**2-D
      CORDE=SQRT(RR)
      DEL=C-CORDE
      SUM=SUM+DEL
      Z(LR+NREG-K+1)=DEL
      VAP(K)=VAP(K)+DEL*S
      C=CORDE
  130 CONTINUE
      IF(KMAX.NE.1) THEN
         DEL=2.0*CORDE
         SUM=SUM+DEL
         Z(LR+KMAX)=DEL
         VAP(K)=VAP(K)+DEL*S
         DO 140 I=1,KMAX-2
         DEL=Z(LR+KMAX-I)
         SUM=SUM+DEL
         Z(LR+KMAX+I)=DEL
         VAP(K+I)=VAP(K+I)+DEL*S
  140    CONTINUE
      ENDIF
      LR=LR+2*KMAX-1
      DEL=B/SI-SUM
      Z(LR)=DEL
      VAP(NREG)=VAP(NREG)+DEL*S
  150 CONTINUE
      DO 160 K=KMAX-1,1-KMAX,-1
      IZ(LI+K+KMAX)=3+ABS(K)+1+NREG-KMAX
  160 CONTINUE
      LI=LI+2*KMAX
      IZ(L5)=2*KMAX-1
      IZ(L5+1)=NX
      IZ(L5+2)=2 ! ISURF
      IZ(LI)=3   ! JSURF
*----
*  VOLUME NORMALIZATION
*----
      IF((INORM.EQ.0).AND.(K0.LT.NREG)) THEN
         DLIM1=BO2*CO+(AO2-X2)*SI
         DLIM2=BO2*CO+(AO2-X1)*SI
         VW1=0.0
         SUMVAP=0.0
         DO 200 I=K0,NREG-1
         SUMVAP=SUMVAP+VAP(I)
         RW=RAYRE(I+1)
         VEX1=RW*RW*ACOS(DLIM1/RW)-DLIM1*SQRT(RW*RW-DLIM1*DLIM1)
         IF(RW.GT.DLIM2)
     1   VEX1=VEX1-(RW*RW*ACOS(DLIM2/RW)-DLIM2*SQRT(RW*RW-DLIM2*DLIM2))
         SUM=(VEX1-VW1)/VAP(I)
         PREC=MAX(PREC,ABS(1.0-SUM)*(VEX1-VW1)/AB)
         VW1=VEX1
         VAP(I)=SUM
  200    CONTINUE
         VEX1=0.5*(B*SI-(A-X1-X2)*CO)*(X2-X1)*SI
         VEX2=(X2-X1)*B-VEX1
         SUM=(VEX1-0.5*VW1)/(VEX1-0.5*SUMVAP)
         PREC=MAX(PREC,ABS(1.0-SUM)*(VEX1-0.5*VW1)/AB)
         VEX1=SUM
         SUM=(VEX2-0.5*VW1)/(VEX2-0.5*SUMVAP)
         PREC=MAX(PREC,ABS(1.0-SUM)*(VEX2-0.5*VW1)/AB)
         VEX2=SUM
         DO 220 IX=1,NX
         KMAX=(IZ(L5)+1)/2
         Z(L3+KMAX)=Z(L3+KMAX)*VAP(K0)
         DO 210 I=1,KMAX-2
         Z(L3+KMAX-I)=Z(L3+KMAX-I)*VAP(K0+I)
         Z(L3+KMAX+I)=Z(L3+KMAX+I)*VAP(K0+I)
  210    CONTINUE
         Z(L3+1)=Z(L3+1)*VEX1
         Z(L3+2*KMAX-1)=Z(L3+2*KMAX-1)*VEX2
         L3=L3+2*KMAX
  220    CONTINUE
      ENDIF
      IZ(L4)=IZ(L4)+1
      X1=X2
  230 CONTINUE
*----
*  THIRD ANGULAR DOMAIN
*----
  240 IF(IZ(3).EQ.2) GO TO 350
      IF(PHI.LE.ATAN(A*BI)) GO TO 350
      X1=A/TA
      XLIM=0.5*(B+X1)
      DO 340 K0=NREG,1,-1
      KMAX=NREG-K0+1
      X2=XLIM-RAYRE(K0)/SI
      IF(X2.LE.X1) GO TO 340
      L3=LR+1
      L5=LI+1
      LI=LI+3
      CALL XDRSET(VAP,NREG,0.0)
      DO 270 IX=1,NX
      IF(K0.EQ.NREG) THEN
         S=0.5*(X2-X1)*SI*WX(IX)
         X=X1+0.5*(X2-X1)*(1.0+ZX(IX))
      ELSE
*        FLURIG CHANGE OF VARIABLE.
         S=2.0*(X2-X1)*SI*WXJ(IX)
         X=X1+(X2-X1)*ZXJ(IX)**2
      ENDIF
      Z(LR+1)=S*WA(IA)
      LR=LR+1
      C=AO2*SI-(BO2-X)*CO
      D=AO2*CO+(BO2-X)*SI
      D=D*D
      SUM=0.0
      CORDE=0.0
      DO 250 K=NREG,K0+1,-1
      RR=RAYRE(K)**2-D
      CORDE=SQRT(RR)
      DEL=C-CORDE
      SUM=SUM+DEL
      Z(LR+NREG-K+1)=DEL
      VAP(K)=VAP(K)+DEL*S
      C=CORDE
  250 CONTINUE
      IF(KMAX.NE.1) THEN
         DEL=2.0*CORDE
         SUM=SUM+DEL
         Z(LR+KMAX)=DEL
         VAP(K)=VAP(K)+DEL*S
         DO 260 I=1,KMAX-2
         DEL=Z(LR+KMAX-I)
         SUM=SUM+DEL
         Z(LR+KMAX+I)=DEL
         VAP(K+I)=VAP(K+I)+DEL*S
  260    CONTINUE
      ENDIF
      LR=LR+2*KMAX-1
      DEL=A/SI-SUM
      Z(LR)=DEL
      VAP(NREG)=VAP(NREG)+DEL*S
  270 CONTINUE
      DO 280 K=KMAX-1,1-KMAX,-1
      IZ(LI+K+KMAX)=3+ABS(K)+1+NREG-KMAX
  280 CONTINUE
      LI=LI+2*KMAX
      IZ(L5)=2*KMAX-1
      IZ(L5+1)=NX
      IZ(L5+2)=0 ! ISURF
      IZ(LI)=1   ! JSURF
*----
*  VOLUME NORMALIZATION
*----
      IF((INORM.EQ.0).AND.(K0.LT.NREG)) THEN
         DLIM1=AO2*CO+(BO2-X2)*SI
         DLIM2=AO2*CO+(BO2-X1)*SI
         VW1=0.0
         SUMVAP=0.0
         DO 310 I=K0,NREG-1
         SUMVAP=SUMVAP+VAP(I)
         RW=RAYRE(I+1)
         VEX1=RW*RW*ACOS(DLIM1/RW)-DLIM1*SQRT(RW*RW-DLIM1*DLIM1)
         IF(RW.GT.DLIM2)
     1   VEX1=VEX1-(RW*RW*ACOS(DLIM2/RW)-DLIM2*SQRT(RW*RW-DLIM2*DLIM2))
         SUM=(VEX1-VW1)/VAP(I)
         PREC=MAX(PREC,ABS(1.0-SUM)*(VEX1-VW1)/AB)
         VW1=VEX1
         VAP(I)=SUM
  310    CONTINUE
         VEX1=0.5*(A*SI-(B-X1-X2)*CO)*(X2-X1)*SI
         VEX2=(X2-X1)*A-VEX1
         SUM=(VEX1-0.5*VW1)/(VEX1-0.5*SUMVAP)
         PREC=MAX(PREC,ABS(1.0-SUM)*(VEX1-0.5*VW1)/AB)
         VEX1=SUM
         SUM=(VEX2-0.5*VW1)/(VEX2-0.5*SUMVAP)
         PREC=MAX(PREC,ABS(1.0-SUM)*(VEX2-0.5*VW1)/AB)
         VEX2=SUM
         DO 330 IX=1,NX
         KMAX=(IZ(L5)+1)/2
         Z(L3+KMAX)=Z(L3+KMAX)*VAP(K0)
         DO 320 I=1,KMAX-2
         Z(L3+KMAX-I)=Z(L3+KMAX-I)*VAP(K0+I)
         Z(L3+KMAX+I)=Z(L3+KMAX+I)*VAP(K0+I)
  320    CONTINUE
         Z(L3+1)=Z(L3+1)*VEX1
         Z(L3+2*KMAX-1)=Z(L3+2*KMAX-1)*VEX2
         L3=L3+2*KMAX
  330    CONTINUE
      ENDIF
      IZ(L4)=IZ(L4)+1
      X1=X2
  340 CONTINUE
  350 CONTINUE
      ZN1=0.5*ZN1*PIO2
      ZN2=0.5*ZN2*PIO2
      ZN3=0.5*ZN3*PIO2
      Z(NREG+1)=1.0/SQRT(ZN1)
      Z(NREG+2)=1.0/SQRT(0.75*ZN3-0.7205061948*ZN2*ZN2/ZN1)
      Z(NREG+3)=Z(NREG+2)*0.8488263632*ZN2/ZN1
      Z(NREG+4)=2.0/SQRT(3.0*(ZN1-ZN3))
*----
*  TRACKING INFORMATION OUTPUT
*----
      IF(ILIGN.EQ.1) THEN
         L1I=IZ(1)-1
         L1R=IZ(2)-1
         WRITE(6,500) (Z(L1R+I),I=1,4)
         L1R=L1R+4
         L2=0
         DO 380 IA=1,NA2
         MNT=IZ(L1I+1)
         L1I=L1I+2
         ZSIN=Z(L1R+1)
         ZCOS=Z(L1R+2)
         L1R=L1R+9
         DO 370 IMNT=1,MNT
         NH=IZ(L1I+1)
         NX=IZ(L1I+2)
         L1I=L1I+3
         DO 360 IX=1,NX
         L2=L2+1
         IF((IMNT.EQ.1).AND.(IX.EQ.1)) THEN
            WRITE(6,510) L2,ZSIN,ZCOS,Z(L1R+1),NH,(Z(L1R+I+1),I=1,NH)
         ELSE
            WRITE(6,520) L2,Z(L1R+1),NH,(Z(L1R+I+1),I=1,NH)
         ENDIF
         L1R=L1R+NH+1
  360    CONTINUE
         L1I=L1I+NH+1
  370    CONTINUE
  380    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(VAP)
      RETURN
*
  500 FORMAT (1H1//30H TRACKING INFORMATION LISTING.//12H NUMERICAL O,
     1 27HRTHONORMALIZATION FACTORS =,1P,4E12.4//6H TRACK)
  510 FORMAT (1X,I5,7H  SIN =,1P,E10.3,7H  COS =,E10.3,9H WEIGHT =,
     1 E10.3,6H  NH =,I3,12H  SEGMENTS =,5E10.3:/(80X,5E10.3))
  520 FORMAT (1X,I5,34X,9H WEIGHT =,1P,E10.3,6H  NH =,I3,10H  SEGMENTS,
     1 2H =,5E10.3:/(80X,5E10.3))
      END
