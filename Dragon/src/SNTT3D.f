*DECK SNTT3D
      SUBROUTINE SNTT3D(IGE,IMPX,LX,LY,LZ,SIDE,IELEM,ISPLH,NLF,NPQ,NSCT,
     1 IQUAD,NCODE,ZCODE,MAT,XXX,YYY,ZZZ,VOL,IDL,DU,DE,DZ,W,MRMX,MRMY,
     2 MRMZ,DC,DB,DA,PL,LL4,NUN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a 3-D Cartesian with discrete ordinates
* approximation of the flux.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): N. Martin
*
*Parameters: input
* IGE     geometry type with (=1 for Cartesian, =2 for hexagonal).
* IMPX    print parameter.
* LX      number of elements along the X axis.
* LY      number of elements along the Y axis.
* LZ      number of elements along the Z axis.
* SIDE    side of hexagon.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant - only for HODD, classical diamond scheme 
*         (default for HODD);
*         =2 linear - default for DG;
*         =3 parabolic;
*         =4 cubic - only for DG.
* ISPLH   hexagone splitting option.
* NLF     SN order for the flux (even number).
* NPQ     number of SN directions in eight octants.
* NSCT    maximum number of spherical harmonics moments of the flux.
* IQUAD   type of SN quadrature (1 Level symmetric, type IQUAD;
*         4 Legendre-Chebyshev; 5 symmetric Legendre-Chebyshev;
*         6 quadruple range).
* NCODE   type of boundary condition applied on each side
*         (i=1 X-;  i=2 X+;  i=3 Y-;  i=4 Y+;  i=5 Z-;  i=6 Z+):
*         =1 VOID;   =2 REFL;  =4 TRAN.
* ZCODE   ZCODE(I) is the albedo corresponding to boundary condition
*         'VOID' on each side (ZCODE(I)=0.0 by default).
* MAT     mixture index assigned to each element.
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
* ZZZ     Cartesian coordinates along the Z axis.
*
*Parameters: output
* VOL     volume of each element.
* IDL     isotropic flux indices.
* DU      first direction cosines ($\\mu$).
* DE      second direction cosines ($\\eta$).
* DZ      third direction cosines ($\\xi$).
* W       weights.
* MRMX    quadrature index.
* MRMY    quadrature index.
* MRMZ    quadrature index.
* DC      diamond-scheme parameter.
* DB      diamond-scheme parameter.
* DA      diamond-scheme parameter.
* PL      discrete values of the spherical harmonics corresponding
*         to the 3D SN quadrature.
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
      INTEGER IMPX,LX,LY,LZ,IELEM,NLF,NPQ,NSCT,IQUAD,NCODE(6),
     1 MAT(LX,LY,LZ),IDL(LX*LY*LZ),MRMX(NPQ),MRMY(NPQ),MRMZ(NPQ)
      REAL ZCODE(6),VOL(LX,LY,LZ),XXX(LX+1),YYY(LY+1),ZZZ(LZ+1),
     1 DU(NPQ),DE(NPQ),DZ(NPQ),W(NPQ),DC(LX,LY,NPQ),DB(LX,LZ,NPQ),
     2 DA(LY,LZ,NPQ),PL(NSCT,NPQ)
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131
      LOGICAL L1,L2,L3,L4,L5,L6
      PARAMETER(RLOG=1.0E-8,PI=3.141592654)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JOP
      REAL, ALLOCATABLE, DIMENSION(:) :: XX,YY,ZZ,UU,WW,TPQ,UPQ,VPQ,WPQ
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(XX(LX),YY(LY),ZZ(LZ))
*----
*  UNFOLD HEIGHT-OCTANT QUADRATURES.
*----
      IF(MOD(NLF,2).EQ.1) CALL XABORT('SNTT3D: EVEN NLF EXPECTED.')
      NPQ0=NLF*(NLF/2+1)/4
      ALLOCATE(JOP(NLF/2),UU(NLF/2),WW(NLF/2),TPQ(NPQ0),UPQ(NPQ0),
     1 VPQ(NPQ0),WPQ(NPQ0))
      IF(IQUAD.EQ.1) THEN
!Level-symmetric quadrature of type 1
         CALL SNQU01(NLF,JOP,UU,WW,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.2) THEN
!Level-symmetric quadrature of type 2
         CALL SNQU02(NLF,JOP,UU,WW,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.3) THEN
!Snow Level-symmetric type quadrature 
         CALL SNQU03(NLF,JOP,UU,WW,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.4) THEN
!Legendre-Chebyshev quadrature 
         CALL SNQU04(NLF,JOP,UU,WW,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.5) THEN
!Symmetric Legendre-Chebyshev quadrature
         CALL SNQU05(NLF,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.6) THEN
!Quadruple Range quadrature
         CALL SNQU06(NLF,TPQ,UPQ,VPQ,WPQ)
      ELSE
         CALL XABORT('SNTT3D: UNKNOWN QUADRATURE TYPE.')
      ENDIF
      M=0
      IOF=0
      DO 320 I=1,NLF/2
         JOF=IOF+NLF-2*I+2
         DO 330 J=0,NLF/2-I
            IOF=IOF+1
            KOF=IOF+(NLF+2)*NLF/4
            LOF=IOF+(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=-UPQ(M+J+1)
            DE(IOF)=-VPQ(M+J+1)
            DZ(IOF)=-TPQ(M+J+1)
            W(IOF)=WPQ(M+J+1)
            JOF=JOF-1
 330     CONTINUE
         DO 340 J=NLF/2-I,0,-1
            IOF=IOF+1
            KOF=IOF+(NLF+2)*NLF/4
            LOF=IOF+(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=UPQ(M+J+1)
            DE(IOF)=-VPQ(M+J+1)
            DZ(IOF)=-TPQ(M+J+1)
            W(IOF)=WPQ(M+J+1)
            JOF=JOF-1
 340     CONTINUE
         M=M+NLF/2-I+1
 320  CONTINUE
      M=0
      DO 350 I=1,NLF/2
         JOF=IOF+NLF-2*I+2
         DO 360 J=0,NLF/2-I
            IOF=IOF+1
            KOF=IOF-(NLF+2)*NLF/4
            LOF=IOF+(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=-UPQ(M+J+1)
            DE(IOF)=VPQ(M+J+1)
            DZ(IOF)=-TPQ(M+J+1)
            W(IOF)=WPQ(M+J+1)
            JOF=JOF-1
 360     CONTINUE
         DO 370 J=NLF/2-I,0,-1
            IOF=IOF+1
            KOF=IOF-(NLF+2)*NLF/4
            LOF=IOF+(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=UPQ(M+J+1)
            DE(IOF)=VPQ(M+J+1)
            DZ(IOF)=-TPQ(M+J+1)
            W(IOF)=WPQ(M+J+1)
            JOF=JOF-1
 370     CONTINUE
         M=M+NLF/2-I+1
 350  CONTINUE
      M=0
      DO 380 I=1,NLF/2
         JOF=IOF+NLF-2*I+2
         DO 390 J=0,NLF/2-I
            IOF=IOF+1
            KOF=IOF+(NLF+2)*NLF/4
            LOF=IOF-(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=-UPQ(M+J+1)
            DE(IOF)=-VPQ(M+J+1)
            DZ(IOF)=TPQ(M+J+1)
            W(IOF)=WPQ(M+J+1)
            JOF=JOF-1
 390     CONTINUE
         DO 400 J=NLF/2-I,0,-1
            IOF=IOF+1
            KOF=IOF+(NLF+2)*NLF/4
            LOF=IOF-(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=UPQ(M+J+1)
            DE(IOF)=-VPQ(M+J+1)
            DZ(IOF)=TPQ(M+J+1)
            W(IOF)=WPQ(M+J+1)
            JOF=JOF-1
 400     CONTINUE
         M=M+NLF/2-I+1
 380  CONTINUE
      M=0
      DO 410 I=1,NLF/2
         JOF=IOF+NLF-2*I+2
         DO 420 J=0,NLF/2-I
            IOF=IOF+1
            KOF=IOF-(NLF+2)*NLF/4
            LOF=IOF-(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=-UPQ(M+J+1)
            DE(IOF)=VPQ(M+J+1)
            DZ(IOF)=TPQ(M+J+1)
            W(IOF)=WPQ(M+J+1)
            JOF=JOF-1
 420     CONTINUE
         DO 430 J=NLF/2-I,0,-1
            IOF=IOF+1
            KOF=IOF-(NLF+2)*NLF/4
            LOF=IOF-(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=UPQ(M+J+1)
            DE(IOF)=VPQ(M+J+1)
            DZ(IOF)=TPQ(M+J+1)
            W(IOF)=WPQ(M+J+1)
            JOF=JOF-1
 430     CONTINUE
         M=M+NLF/2-I+1
 410  CONTINUE
      DEALLOCATE(WPQ,VPQ,UPQ,TPQ,WW,UU,JOP)
     
      IF(IMPX.GE.4) THEN
         WRITE(6,'(/41H SNTT3D:HEIGHT-OCTANT ANGULAR QUADRATURES:/26X,
     1   2HMU,9X,3HETA,10X,2HXI,6X,6HWEIGHT)')
         SUM=0.0
         DO 70 M=1,NPQ
         SUM=SUM+W(M)
         WRITE(6,'(1X,4I5,1P,4E12.4)') M,MRMX(M),MRMY(M),MRMZ(M),DU(M),
     1   DE(M),DZ(M),W(M)
   70    CONTINUE
         WRITE(6,'(54X,10(1H-)/52X,1P,E12.4)') SUM
      ENDIF
*----
*  IDENTIFICATION OF THE GEOMETRY.
*----
      IF(IGE.EQ.0) THEN
* ----------
*        3D CARTESIAN
* ----------
         DO 83 M=1,NPQ
         VU=DU(M)
         VE=DE(M)
         VZ=DZ(M)
         DO 82 I=1,LX
         DO 81 J=1,LY
         DO 80 K=1,LZ
         XX(I)=XXX(I+1)-XXX(I)
         YY(J)=YYY(J+1)-YYY(J)
         ZZ(K)=ZZZ(K+1)-ZZZ(K)
         DA(J,K,M)=VU*YY(J)*ZZ(K)
         DB(I,K,M)=VE*XX(I)*ZZ(K)
         DC(I,J,M)=VZ*XX(I)*YY(J)
         VOL(I,J,K)=XX(I)*YY(J)*ZZ(K)
   80    CONTINUE
   81    CONTINUE
   82    CONTINUE
   83    CONTINUE
      ELSEIF(IGE.EQ.2) THEN
* ----------
*        3D HEXAGONAL
* ----------
         DET = SQRT(3.0)*(SIDE**2)/2.0 
         DO 93 M=1,NPQ
         VU=DU(M)
         VE=DE(M)
         VZ=DZ(M)
         DO 92 K=1,LZ
         DO 91 J=1,LY
         DO 90 I=1,LX
         ZZ(K)=ZZZ(K+1)-ZZZ(K)
         DA(J,K,M)=VU*ZZ(K)
         DB(I,K,M)=VE*ZZ(K)
         DC(I,J,M)=VZ*DET
         VOL(I,J,K)=DET*ZZ(K)
   90    CONTINUE
   91    CONTINUE
   92    CONTINUE
   93    CONTINUE
      ENDIF
*----
*  GENERATE SPHERICAL HARMONICS FOR SCATTERING SOURCE.
*----
      IOF=0
      DO 215 IL=0,NSCT-1
      DO 210 IM=-IL,IL   
      IOF=IOF+1
      IF(IOF.GT.NSCT) GO TO 220
      DO 200 M=1,NPQ
      PL(IOF,M)=PNSH(IL,IM,DU(M),DE(M),DZ(M))
  200 CONTINUE
  210 CONTINUE
  215 CONTINUE
*----
*  COMPUTE ISOTROPIC FLUX INDICES.
*----
  220 LL4=LX*LY*LZ*NSCT*IELEM*IELEM*IELEM
      IF(IGE.EQ.0)THEN
         NUN=LL4+(LX*(LY+LZ)+LY*LZ)*NPQ*IELEM*IELEM
      ELSEIF(IGE.EQ.2)THEN
         NUN=LL4+(LX*LY*LZ)*ISPLH*ISPLH*NPQ*IELEM*IELEM
      ELSE
         CALL XABORT('SNTT3D: CHECK SPATIAL SCHEME DISCRETISATION '//
     1       'PARAMETER.')
      ENDIF
      DO I=1,LX*LY*LZ
         IDL(I)=(I-1)*NSCT*IELEM*IELEM*IELEM+1
      ENDDO
*----
*  SET BOUNDARY CONDITIONS.
*----
      DO 240 I=1,6
      IF(NCODE(I).NE.1) ZCODE(I)=1.0
      IF(NCODE(I).EQ.7) CALL XABORT('SNTT3D: ZERO FLUX BC NOT ALLOWED.')
  240 CONTINUE
*----
*  CHECK FOR INVALID VIRTUAL ELEMENTS.
*----
      DO 292 I=2,LX-1
      DO 291 J=2,LY-1
      DO 290 K=2,LZ-1
      IF(MAT(I,J,K).EQ.0) THEN
         L1=(NCODE(1).NE.1)
         DO 251 J1=1,J-1
         DO 250 K1=1,K-1
         L1=L1.OR.(MAT(I,J1,K1).NE.0)
  250    CONTINUE
  251    CONTINUE
         L2=(NCODE(2).NE.1)
         DO 261 J1=J+1,LY
         DO 260 K1=K+1,LZ
         L2=L2.OR.(MAT(I,J1,K1).NE.0)
  260    CONTINUE
  261    CONTINUE
         L3=(NCODE(3).NE.1)
         DO 271 I1=1,I-1
         DO 270 K1=1,K-1
         L3=L3.OR.(MAT(I1,J,K1).NE.0)
  270    CONTINUE
  271    CONTINUE
         L4=(NCODE(4).NE.1)
         DO 281 I1=I+1,LX
         DO 280 K1=K+1,LZ
         L4=L4.OR.(MAT(I1,J,K1).NE.0)
  280    CONTINUE
  281    CONTINUE
         L5=(NCODE(5).NE.1)
         DO 301 I1=1,I-1
         DO 300 J1=1,J-1
         L5=L5.OR.(MAT(I1,J1,K).NE.0)
  300    CONTINUE
  301    CONTINUE
         L6=(NCODE(6).NE.0)
         DO 311 I1=I+1,LX
         DO 310 J1=I+1,LY
         L6=L6.OR.(MAT(I1,J1,K).NE.0)
  310    CONTINUE
  311    CONTINUE
         IF(L1.AND.L2.AND.L3.AND.L4.AND.L5.AND.L6) THEN
            WRITE(HSMG,'(17HSNTT3D: ELEMENT (,I3,1H,,I3,11H) CANNOT BE,
     1      9H VIRTUAL.)') I,J,K
            CALL XABORT(HSMG)
         ENDIF
      ENDIF
  290 CONTINUE
  291 CONTINUE
  292 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ZZ,YY,XX)
      RETURN
      END
