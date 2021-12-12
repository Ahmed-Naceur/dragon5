*DECK SNTT2D
      SUBROUTINE SNTT2D (IGE,IMPX,LX,LY,SIDE,IELEM,NLF,NPQ,NSCT,IQUAD,
     1 NCODE,ZCODE,MAT,XXX,YYY,VOL,IDL,DU,DE,W,MRM,MRMY,DB,DA,DAL,PL,
     2 LL4,NUN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a 2-D Cartesian or R-Z geometry with
* discrete ordinates approximation of the flux.
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
* IGE     type of 2D geometry (=0 Cartesian; =1 R-Z; =2 Hexagonal).
* IMPX    print parameter.
* LX      number of elements along the X axis.
* LY      number of elements along the Y axis.
* SIDE    side of an hexagon.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant - only for HODD, classical diamond scheme 
*         (default for HODD);
*         =2 linear - default for DG;
*         =3 parabolic;
*         =4 cubic - only for DG.
* NLF     SN order for the flux (even number).
* NPQ     number of SN directions in four octants (including zero-weight
*         directions).
* NSCT    maximum number of spherical harmonics moments of the flux.
* IQUAD   type of SN quadrature (1 Level symmetric, type IQUAD;
*         4 Legendre-Chebyshev; 5 symmetric Legendre-Chebyshev;
*         6 quadruple range).
* NCODE   type of boundary condition applied on each side
*         (i=1 X-;  i=2 X+;  i=3 Y-;  i=4 Y+):
*         =1: VOID; =2: REFL; =4: TRAN.
* ZCODE   ZCODE(I) is the albedo corresponding to boundary condition
*         'VOID' on each side (ZCODE(I)=0.0 by default).
* MAT     mixture index assigned to each element.
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
*
*Parameters: output
* VOL     volume of each element.
* IDL     isotropic flux indices.
* DU      first direction cosines ($\\mu$).
* DE      second direction cosines ($\\eta$).
* W       weights.
* MRM     quadrature index.
* MRMY    quadrature index.
* DB      diamond-scheme parameter.
* DA      diamond-scheme parameter.
* DAL     diamond-scheme angular redistribution parameter.
* PL      discrete values of the spherical harmonics corresponding
*         to the 2D SN quadrature.
* LL4     number of unknowns being solved for, over the domain. This 
*         includes the various moments of the isotropic (and if present,
*         anisotropic) flux. 
* NUN     total number of unknowns stored in the FLUX vector per group.
*         This includes LL4 (see above) as well as any surface boundary
*         fluxes, if present.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IGE,IMPX,LX,LY,IELEM,NLF,NPQ,NSCT,IQUAD,NCODE(4),
     1 MAT(LX,LY),IDL(LX*LY),MRM(NPQ),MRMY(NPQ),LL4,NUN
      REAL ZCODE(4),VOL(LX,LY),XXX(LX+1),YYY(LY+1),DU(NPQ),DE(NPQ),
     1 W(NPQ),DB(LX,NPQ),DA(LX,LY,NPQ),DAL(LX,LY,NPQ),PL(NSCT,NPQ)
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131
      LOGICAL L1,L2,L3,L4
      PARAMETER(RLOG=1.0E-8,PI=3.141592654)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JOP
      REAL, ALLOCATABLE, DIMENSION(:) :: XX,YY,UU,WW,TPQ,UPQ,VPQ,WPQ
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(XX(LX),YY(LY))
*----
*  UNFOLD FOUR-OCTANT QUADRATURES.
*----
      IF(MOD(NLF,2).EQ.1) CALL XABORT('SNTT2D: EVEN NLF EXPECTED.')
      NPQ0=NLF*(NLF/2+1)/4
      ALLOCATE(JOP(NLF/2),UU(NLF/2),WW(NLF/2),TPQ(NPQ0),UPQ(NPQ0),
     1 VPQ(NPQ0),WPQ(NPQ0))
      IF(IQUAD.EQ.1) THEN
         CALL SNQU01(NLF,JOP,UU,WW,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.2) THEN
         CALL SNQU02(NLF,JOP,UU,WW,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.3) THEN
         CALL SNQU03(NLF,JOP,UU,WW,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.4) THEN
         CALL SNQU04(NLF,JOP,UU,WW,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.5) THEN
         CALL XDRSET(UU,NLF/2,0.0)
         CALL SNQU05(NLF,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.6) THEN
         CALL XDRSET(UU,NLF/2,0.0)
         CALL SNQU06(NLF,TPQ,UPQ,VPQ,WPQ)
      ELSE
         CALL XABORT('SNTT2D: UNKNOWN QUADRATURE TYPE.')
      ENDIF
      M=0
      IOF=0
      DO 30 I=1,NLF/2
      IOF=IOF+1
      JOF=IOF+NLF-2*I+2
      KOF=IOF+(NLF+4)*NLF/4
      MRM(IOF)=JOF
      MRMY(IOF)=KOF
      DU(IOF)=-SQRT(1.0-UU(I)*UU(I))
      DE(IOF)=-UU(I)
      W(IOF)=0.0
      DO 10 J=0,NLF/2-I
      IOF=IOF+1
      KOF=IOF+(NLF+4)*NLF/4
      MRM(IOF)=JOF
      MRMY(IOF)=KOF
      DU(IOF)=-UPQ(M+J+1)
      DE(IOF)=-VPQ(M+J+1)
      W(IOF)=WPQ(M+J+1)
      JOF=JOF-1
   10 CONTINUE
      DO 20 J=NLF/2-I,0,-1
      IOF=IOF+1
      KOF=IOF+(NLF+4)*NLF/4
      MRM(IOF)=JOF
      MRMY(IOF)=KOF
      DU(IOF)=UPQ(M+J+1)
      DE(IOF)=-VPQ(M+J+1)
      W(IOF)=WPQ(M+J+1)
      JOF=JOF-1
   20 CONTINUE
      M=M+NLF/2-I+1
   30 CONTINUE
      M=0
      DO 60 I=1,NLF/2
      IOF=IOF+1
      JOF=IOF+NLF-2*I+2
      KOF=IOF-(NLF+4)*NLF/4
      MRM(IOF)=JOF
      MRMY(IOF)=KOF
      DU(IOF)=-SQRT(1.0-UU(I)*UU(I))
      DE(IOF)=UU(I)
      W(IOF)=0.0
      DO 40 J=0,NLF/2-I
      IOF=IOF+1
      KOF=IOF-(NLF+4)*NLF/4
      MRM(IOF)=JOF
      MRMY(IOF)=KOF
      DU(IOF)=-UPQ(M+J+1)
      DE(IOF)=VPQ(M+J+1)
      W(IOF)=WPQ(M+J+1)
      JOF=JOF-1
   40 CONTINUE
      DO 50 J=NLF/2-I,0,-1
      IOF=IOF+1
      KOF=IOF-(NLF+4)*NLF/4
      MRM(IOF)=JOF
      MRMY(IOF)=KOF
      DU(IOF)=UPQ(M+J+1)
      DE(IOF)=VPQ(M+J+1)
      W(IOF)=WPQ(M+J+1)
      JOF=JOF-1
   50 CONTINUE
      M=M+NLF/2-I+1
   60 CONTINUE
      DEALLOCATE(WPQ,VPQ,UPQ,TPQ,WW,UU,JOP)
      IF(IMPX.GE.4) THEN
         WRITE(6,'(/41H SNTT2D: FOUR-OCTANT ANGULAR QUADRATURES:/26X,
     1   2HMU,9X,3HETA,10X,2HXI,6X,6HWEIGHT)')
         SUM=0.0
         DO 70 M=1,NPQ
         SUM=SUM+W(M)
         ZI=SQRT(ABS(1.0-DU(M)**2-DE(M)**2))
         IF(ZI.LT.1.0E-3) ZI=0.0
         WRITE(6,'(1X,3I5,1P,4E12.4)') M,MRM(M),MRMY(M),DU(M),DE(M),ZI,
     1   W(M)
   70    CONTINUE
         WRITE(6,'(54X,10(1H-)/52X,1P,E12.4)') SUM
      ENDIF
*----
*  IDENTIFICATION OF THE GEOMETRY.
*----
      IF(IGE.EQ.0) THEN
* ----------
*        2D CARTESIAN
* ----------
         DO 82 M=1,NPQ
         VU=DU(M)
         VE=DE(M)
         DO 81 I=1,LX
         XX(I)=XXX(I+1)-XXX(I)
         DB(I,M)=VE*XX(I)
         DO 80 J=1,LY
         YY(J)=YYY(J+1)-YYY(J)
         DA(I,J,M)=VU*YY(J)
         DAL(I,J,M)=0.0
   80    CONTINUE
   81    CONTINUE
   82    CONTINUE
         DO 91 I=1,LX
         DO 90 J=1,LY
         VOL(I,J)=XX(I)*YY(J)
   90    CONTINUE
   91    CONTINUE
      ELSEIF(IGE.EQ.1) THEN
* ----------
*        2D TUBE
* ----------
         DO 95 J=1,LY
         YY(J)=YYY(J+1)-YYY(J)
   95    CONTINUE
         DO 102 M=1,NPQ
         VU=DU(M)*PI
         DO 101 I=1,LX
         XX(I)=XXX(I+1)-XXX(I)
         VE=(XXX(I)+XXX(I+1))*VU
         DO 100 J=1,LY
         DA(I,J,M)=VE*YY(J)
  100    CONTINUE
  101    CONTINUE
  102    CONTINUE
         CALL XDRSET(DAL,LY*LX*NPQ,0.)
         CALL XDRSET(DB,LX*NPQ,0.)
         DO 135 J=1,LY
         DO 111 I=1,LX
         VE=2.0*PI*(XXX(I+1)-XXX(I))*YY(J)
         DO 110 M=2,NPQ
         DB(I,M)=DB(I,M-1)-W(M)*DU(M)*VE
  110    CONTINUE
  111    CONTINUE
         DO 130 M=2,NPQ
         VE=W(M)
         IF(VE.LE.RLOG) GOTO 130
         DO 120 I=1,LX
         DAL(I,J,M)=(DB(I,M)+DB(I,M-1))/VE
  120    CONTINUE
  130    CONTINUE
  135    CONTINUE
         DO 155 I=1,LX
         VE=PI*XX(I)*(XXX(I+1)+XXX(I))
         DO 140 M=1,NPQ
         DB(I,M)=VE*DE(M)
  140    CONTINUE
         DO 150 J=1,LY
         VOL(I,J)=YY(J)*VE
  150    CONTINUE
  155    CONTINUE
      ELSEIF(IGE.EQ.2) THEN
* ----------
*        2D HEXAGONAL
* ----------
         DET = SQRT(3.0)*(SIDE**2)/2.0 
         DO 162 M=1,NPQ
         VU=DU(M)
         VE=DE(M)
         DO 161 I=1,LX
         DB(I,M)=VE
         DO 160 J=1,LY
         DA(I,J,M)=VU
         VOL(I,J)=DET
  160    CONTINUE
  161    CONTINUE
  162    CONTINUE
      ENDIF
*----
*  GENERATE SPHERICAL HARMONICS FOR SCATTERING SOURCE.
*----
      IOF=0
      DO 215 IL=0,NSCT-1
      DO 210 IM=-IL,IL
      IF(MOD(IL+IM,2).EQ.1) GO TO 210
      IOF=IOF+1
      IF(IOF.GT.NSCT) GO TO 220
      DO 200 M=1,NPQ
      ZI=SQRT(ABS(1.0-DU(M)**2-DE(M)**2))
      IF(ZI.LT.1.0E-3) ZI=0.0
      PL(IOF,M)=PNSH(IL,IM,ZI,DU(M),DE(M))
  200 CONTINUE
  210 CONTINUE
  215 CONTINUE
*----
*  COMPUTE ISOTROPIC FLUX INDICES.
*----
  220 LL4=LX*LY*NSCT*IELEM*IELEM
      IF(IGE.LT.2)THEN
         NUN=LL4 + (LX+LY)*NPQ*IELEM
         DO I=1,LX*LY
            IDL(I)=(I-1)*NSCT*IELEM*IELEM +1
         ENDDO
      ELSEIF(IGE.EQ.2)THEN
         NUN=LL4
         DO I=1,LX
            IDL(I)=(I-1)*NSCT*IELEM*IELEM +1
         ENDDO
      ELSE
         CALL XABORT('SNTT2D: CHECK SPATIAL SCHEME DISCRETISATION '//
     1       'PARAMETER.')
      ENDIF
*----
*  SET BOUNDARY CONDITIONS.
*----
      DO 240 I=1,4
      IF(NCODE(I).NE.1) ZCODE(I)=1.0
      IF(NCODE(I).EQ.7) CALL XABORT('SNTT2D: ZERO FLUX BC NOT ALLOWED.')
  240 CONTINUE
*----
*  CHECK FOR INVALID VIRTUAL ELEMENTS.
*----
      DO 295 I=2,LX-1
      DO 290 J=2,LY-1
      IF(MAT(I,J).EQ.0) THEN
         L1=(NCODE(1).NE.1)
         DO 250 J1=1,J-1
         L1=L1.OR.(MAT(I,J1).NE.0)
  250    CONTINUE
         L2=(NCODE(2).NE.1)
         DO 260 J1=J+1,LY
         L2=L2.OR.(MAT(I,J1).NE.0)
  260    CONTINUE
         L3=(NCODE(3).NE.1)
         DO 270 I1=1,I-1
         L3=L3.OR.(MAT(I1,J).NE.0)
  270    CONTINUE
         L4=(NCODE(4).NE.1)
         DO 280 I1=I+1,LX
         L4=L4.OR.(MAT(I1,J).NE.0)
  280    CONTINUE
         IF(L1.AND.L2.AND.L3.AND.L4) THEN
            WRITE(HSMG,'(17HSNTT2D: ELEMENT (,I3,1H,,I3,11H) CANNOT BE,
     1      9H VIRTUAL.)') I,J
            CALL XABORT(HSMG)
         ENDIF
      ENDIF
  290 CONTINUE
  295 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(YY,XX)
      RETURN
      END
