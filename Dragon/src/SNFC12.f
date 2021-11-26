*DECK SNFC12
      SUBROUTINE SNFC12(LX,LY,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,NCODE,ZCODE,
     1 QEXT,LFIXUP,DU,DE,W,MRM,MRMY,DB,DA,DAL,PL,FLUX,XNEI,XNEJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 2D R-Z
* geometry for the diamond differencing method. Albedo boundary
* conditions.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* LX      number of meshes along X axis.
* LY      number of meshes along Y axis.
* NMAT    number of material mixtures.
* NPQ     number of SN directions in four octants (including zero-weight
*         directions).
* NSCT    maximum number of spherical harmonics moments of the flux.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* TOTAL   macroscopic total cross sections.
* NCODE   boundary condition indices.
* ZCODE   albedos.
* QEXT    Legendre components of the fixed source.
* LFIXUP  flag to enable negative flux fixup.
* DU      first direction cosines ($\\mu$).
* DE      second direction cosines ($\\eta$).
* W       weights.
* MRM     quadrature index.
* MRMY    quadrature index.
* DB      diamond-scheme parameter.
* DA      diamond-scheme parameter.
* DAL     diamond-scheme parameter.
* PL      discrete values of the spherical harmonics corresponding
*         to the 2D SN quadrature.
*
*Parameters: input/output
* XNEI    X-directed SN boundary fluxes.
* XNEJ    Y-directed SN boundary fluxes.
*
*Parameters: output
* FLUX    Legendre components of the flux.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LX,LY,NMAT,NPQ,NSCT,MAT(LX,LY),NCODE(4),MRM(NPQ),MRMY(NPQ)
      REAL VOL(LX,LY),TOTAL(0:NMAT),ZCODE(4),QEXT(NSCT,LX,LY),DU(NPQ),
     1 DE(NPQ),W(NPQ),DB(LX,NPQ),DA(LX,LY,NPQ),DAL(LX,LY,NPQ),
     2 PL(NSCT,NPQ),FLUX(NSCT,LX,LY),XNEI(LY,NPQ),XNEJ(LX,NPQ)
      LOGICAL LFIXUP
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION QQ,C1,XNM,XNJ,Q2(1,2)
      PARAMETER(RLOG=1.0E-8,PI=3.141592654)
      REAL, ALLOCATABLE, DIMENSION(:,:) :: XARN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: XNI
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(XARN(LX,LY),XNI(LY))
*----
*  MAIN LOOP OVER SN ANGLES.
*----
      CALL XDRSET(FLUX,LX*LY*NSCT,0.0)
      CALL XDRSET(XARN,LX*LY,0.0)
      C1=0.0
      XNM=0.0
      DO 170 M=1,NPQ
      WEIGHT=W(M)
      VU=DU(M)
      VE=DE(M)
      IF(NCODE(1).NE.4) THEN
         M1=MRM(M)
         IF(WEIGHT.EQ.0.0) THEN
            DO 10 J=1,LY
            XNEI(J,M)=XNEI(J,M1)
   10       CONTINUE
         ELSE IF(VU.GT.0.0) THEN
            DO 20 J=1,LY
            E1=XNEI(J,M)
            XNEI(J,M)=XNEI(J,M1)
            XNEI(J,M1)=E1
   20       CONTINUE
         ENDIF
      ENDIF
      IF(NCODE(3).NE.4) THEN
         M1=MRMY(M)
         IF(VE.GT.0) THEN
            DO 40 I=1,LX
            E1=XNEJ(I,M)
            XNEJ(I,M)=XNEJ(I,M1)
            XNEJ(I,M1)=E1
   40       CONTINUE
         ENDIF
      ENDIF
      IF(VE.GT.0.0) GOTO 70
      IF(VU.GT.0.0) GOTO 60
      IND=3
      GOTO 90
   60 IND=4
      GOTO 90
   70 IF(VU.GT.0.0) GOTO 80
      IND=2
      GOTO 90
   80 IND=1
*----
*  LOOP OVER X- AND Y-DIRECTED AXES.
*----
   90 DO 155 II=1,LX
      I=II
      IF((IND.EQ.2).OR.(IND.EQ.3)) I=LX+1-II
      IF((IND.EQ.1).OR.(IND.EQ.2)) THEN
         XNJ=XNEJ(I,M)*ZCODE(3)
      ELSE
         XNJ=XNEJ(I,M)*ZCODE(4)
      ENDIF
      DO 140 JJ=1,LY
      J=JJ
      IF((IND.EQ.3).OR.(IND.EQ.4)) J=LY+1-JJ
      C1=DAL(I,J,M)
      IF(II.EQ.1) THEN
         IF((IND.EQ.1).OR.(IND.EQ.4)) THEN
            XNI(J)=XNEI(J,M)
         ELSE
            XNI(J)=XNEI(J,M)*ZCODE(2)
         ENDIF
      ENDIF
      IF(MAT(I,J).EQ.0) GO TO 140
      QQ=0.0D0
      DO 110 K=1,NSCT
      QQ=QQ+QEXT(K,I,J)*PL(K,M)/(4.0*PI)
  110 CONTINUE
      VT=VOL(I,J)*TOTAL(MAT(I,J))
      XNM=XARN(I,J)
      Q2(1,1)=C1+2.0D0*ABS(DA(I,J,M))+2.0D0*ABS(DB(I,M))+VT
      Q2(1,2)=C1*XNM+2.0D0*ABS(DA(I,J,M))*XNI(J)+2.0D0*ABS(DB(I,M))
     1          *XNJ+VOL(I,J)*QQ
      IF(Q2(1,1).EQ.0.0D0) CALL XABORT('SNFC12: SINGULAR MATRIX.')
      Q2(1,2)=Q2(1,2)/Q2(1,1)
      IF(LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
      XNI(J)=2.0D0*Q2(1,2)-XNI(J)
      XNJ=2.0D0*Q2(1,2)-XNJ
      XARN(I,J)=REAL(2.0D0*Q2(1,2)-XNM)
      IF(LFIXUP.AND.(XARN(I,J).LE.RLOG)) XARN(I,J)=0.0
      IF(W(M).LE.RLOG) XARN(I,J)=REAL(Q2(1,2))
      IF(LFIXUP.AND.(XNI(J).LE.RLOG)) XNI(J)=0.0
      IF(LFIXUP.AND.(XNJ.LE.RLOG)) XNJ=0.0
      DO 135 K=1,NSCT
      FLUX(K,I,J)=FLUX(K,I,J)+2.0*W(M)*REAL(Q2(1,2))*PL(K,M)
  135 CONTINUE
  140 CONTINUE
      XNEJ(I,M)=REAL(XNJ)
  155 CONTINUE
      DO 165 J=1,LY
      XNEI(J,M)=REAL(XNI(J))
  165 ENDDO
  170 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XNI,XARN)
      RETURN
      END
