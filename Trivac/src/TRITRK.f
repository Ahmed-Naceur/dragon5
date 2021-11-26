*DECK TRITRK
      SUBROUTINE TRITRK (MAXPTS,IPTRK,IPGEOM,IMPX,IELEM,ICOL,ICHX,ISEG,
     1 IMPV,NLF,NVD,ISPN,ISCAT,NADI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover of the geometry and tracking for TRIVAC.
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
* MAXPTS  allocated storage for arrays of dimension NEL.
* IPTRK   L_TRACK pointer to the TRIVAC tracking information.
* IPGEOM  L_GEOM pointer to the geometry.
* IMPX    print flag.
* IELEM   degree of the Lagrangian finite elements:
*         =1: linear finite elements or finite differences;
*         =2: parabolic finite elements;
*         =3: cubic finite elements;
*         =4: quartic finite elements.
* ICOL    type of quadrature used to integrate the mass matrix:
*         =1: analytical integration;
*         =2: Gauss-Lobatto quadrature (collocation method);
*         =3: Gauss-Legendre quadrature (superconvergent)
*         IELEM=1 and ICOL=2 are finite difference approximations.
* ICHX    type of discretization method:
*         =1: variational collocation method (primal finite elements
*             with Gauss-Lobatto quadrature);
*         =2: dual finite element approximations;
*         =3: nodal collocation method with full tensorial products
*            (dual finite elements with Gauss-Lobatto quadrature).
* ISEG    number of elements in a vector register. Equal to zero for
*         operations in scalar mode.
* IMPV    print parameter for supervectorial operations.
* NLF     number of Legendre orders for the flux. Equal to zero for
*         diffusion theory.
* NVD     type of void boundary condition if NLF>0 and ICOL=3.
* ISPN    type of transport solution:
*         =0: complete PN method;
*         =1: simplified PN method.
* ISCAT   source anisotropy:
*         =1: isotropic sources in laboratory system;
*         =2: linearly anisotropic sources in laboratory system.
* NADI    number of ADI iterations at the inner iterative level.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPGEOM
      INTEGER MAXPTS,IMPX,IELEM,ICOL,ICHX,ISEG,IMPV,NLF,NVD,ISPN,ISCAT,
     1 NADI
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      LOGICAL ILK,CYLIND,CHEX
      CHARACTER HSMG*131
      INTEGER ISTATE(NSTATE),IGP(NSTATE),NCODE(6),ICODE(6)
      REAL ZCODE(6)
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT,IDL,IPERT,KN,IQFR,IDP,
     1 IMX,ISPLX,ISPLY,ISPLZ,MUW,MUX,MUY,MUZ,IPW,IPX,IPY,IPZ,ISET
      REAL, DIMENSION(:), ALLOCATABLE :: VOL,XXX,YYY,ZZZ,XX,YY,ZZ,DD,
     1 QFR,FRZ,RR0,XR0,ANG
      REAL, DIMENSION(:,:), ALLOCATABLE :: V,H
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CTRAN
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBLW,LBLW,MUVW,IPVW,NBLX,
     1 LBLX,MUVX,IPVX,NBLY,LBLY,MUVY,IPVY,NBLZ,LBLZ,MUVZ,IPVZ
      REAL, DIMENSION(:), ALLOCATABLE :: BBW,BBX,BBY,BBZ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPBBW,IPBBX,IPBBY,IPBBZ
*
******************* TRIVAC GEOMETRICAL STRUCTURE. **********************
*                                                                      *
*   ITYPE       : =2 : CARTESIAN 1-D GEOMETRY;                         *
*                 =3 : TUBE 1-D GEOMETRY;                              *
*                 =5 : CARTESIAN 2-D GEOMETRY;                         *
*                 =6 : TUBE 2-D GEOMETRY;                              *
*                 =7 : CARTESIAN 3-D GEOMETRY;                         *
*                 =8 : HEXAGONAL 2-D GEOMETRY;                         *
*                 =9 : HEXAGONAL 3-D GEOMETRY.                         *
*   IHEX        : TYPE OF HEXAGONAL SYMMETRY.                          *
*   IDIAG       : =0 NO DIAGONAL SYMMETRY; =1 DIAGONAL SYMMETRY.       *
*   IELEM       : DEGREE OF THE LAGRANGIAN FINITE ELEMENTS.            *
*                 =1: LINEAR FINITE ELEMENTS OR FINITE DIFFERENCES;    *
*                 =2: PARABOLIC FINITE ELEMENTS;                       *
*                 =3: CUBIC FINITE ELEMENTS;                           *
*                 =4: QUARTIC FINITE ELEMENTS.                         *
*   ICOL        : TYPE OF QUADRATURE USED TO INTEGRATE THE MASS MATRIX.*
*                 =1: ANALYTICAL INTEGRATION;                          *
*                 =2: GAUSS-LOBATTO QUADRATURE (COLLOCATION METHOD);   *
*                 =3: GAUSS-LEGENDRE QUADRATURE (SUPERCONVERGENT).     *
*                 IELEM=1 AND ICOL=2 ARE FINITE DIFFERENCE APPROX.     *
*   ICHX        : TYPE OF DISCRETIZATION METHOD.                       *
*                 =1: VARIATIONAL COLLOCATION METHOD (PRIMAL FINITE    *
*                     ELEMENTS WITH GAUSS-LOBATTO QUADRATURE);         *
*                 =2: DUAL FINITE ELEMENT APPROXIMATIONS;              *
*                 =3: NODAL COLLOCATION METHOD WITH FULL TENSORIAL     *
*                     PRODUCTS (DUAL FINITE ELEMENTS WITH GAUSS-       *
*                     LOBATTO QUADRATURE).                             *
*   SIDE        : SIDE OF THE HEXAGONS.                                *
*   LL4         : ORDER OF THE MATRICES PER GROUP IN TRIVAC.           *
*   NCODE       : TYPES OF BOUNDARY CONDITIONS. DIMENSION=6            *
*   ZCODE       : ALBEDOS. DIMENSION=6                                 *
*   LX,LY,LZ    : NUMBER OF ELEMENTS ALONG THE X, Y AND Z AXIS.        *
*   XX          : X-DIRECTED MESH SPACINGS. DIMENSION=LX*LY*LZ         *
*   YY          : Y-DIRECTED MESH SPACINGS. DIMENSION=LX*LY*LZ         *
*   ZZ          : Z-DIRECTED MESH SPACINGS. DIMENSION=LX*LY*LZ         *
*   DD          : USED WITH CYLINDRICAL GEOMETRIES. DIMENSION=LX*LY*LZ *
*   KN          : ELEMENT-ORDERED UNKNOWN LIST. DIMENSION LX*LY*LZ*ICO *
*                 WHERE ICO IS THE NUMBER OF UNKNOWN PER ELEMENT.      *
*   QFR         : ELEMENT-ORDERED BOUNDARY CONDITIONS.                 *
*                 DIMENSION 6*LX*LY*LZ OR 8*LX*LZ                      *
*   IQFR        : ELEMENT-ORDERED PHYSICAL ALBEDO INDICES.             *
*                 DIMENSION 6*LX*LY*LZ OR 8*LX*LZ                      *
*   MUW         : INDICES USED WITH W-DIRECTED COMPRESSED DIAGONAL     *
*                 STORAGE MODE MATRICES. DIMENSION LL4W                *
*   MUX         : INDICES USED WITH X-DIRECTED COMPRESSED DIAGONAL     *
*                 STORAGE MODE MATRICES. DIMENSION LL4X                *
*   MUY         : INDICES USED WITH Y-DIRECTED COMPRESSED DIAGONAL     *
*                 STORAGE MODE MATRICES. DIMENSION LL4Y                *
*   MUZ         : INDICES USED WITH Z-DIRECTED COMPRESSED DIAGONAL     *
*                 STORAGE MODE MATRICES. DIMENSION LL4Z                *
*   IPW         : W-DIRECTED PERMUTATION MATRIX. DIMENSION LL4         *
*   IPX         : X-DIRECTED PERMUTATION MATRIX. DIMENSION LL4         *
*   IPY         : Y-DIRECTED PERMUTATION MATRIX. DIMENSION LL4         *
*   IPZ         : Z-DIRECTED PERMUTATION MATRIX. DIMENSION LL4         *
*                                                                      *
* SUPERVECTORIAL OPERATION INFORMATION:                                *
*   ISEG        : NUMBER OF ELEMENTS IN A VECTOR REGISTER. EQUAL TO    *
*                 ZERO FOR OPERATIONS IN SCALAR MODE.                  *
*   IMPV        : PRINT PARAMETER FOR SUPERVECTORIAL OPERATIONS.       *
*   LTSW        : MAXIMUM BANDWIDTH. =2 FOR TRIDIAGONAL SYSTEMS.       *
*   LONW        : NUMBER OF GROUPS OF LINEAR SYSTEMS FOR W-MATRICES.   *
*   LONX        : NUMBER OF GROUPS OF LINEAR SYSTEMS FOR X-MATRICES.   *
*   LONY        : NUMBER OF GROUPS OF LINEAR SYSTEMS FOR Y-MATRICES.   *
*   LONZ        : NUMBER OF GROUPS OF LINEAR SYSTEMS FOR Z-MATRICES.   *
*   NBLW        : NUMBER OF LINEAR SYSTEMS PER W-GROUP. DIMENSION LONW *
*   NBLX        : NUMBER OF LINEAR SYSTEMS PER X-GROUP. DIMENSION LONX *
*   NBLY        : NUMBER OF LINEAR SYSTEMS PER Y-GROUP. DIMENSION LONY *
*   NBLZ        : NUMBER OF LINEAR SYSTEMS PER Z-GROUP. DIMENSION LONZ *
*   LBLW        : NUMBER OF UNKNOWNS PER W-GROUP. DIMENSION LONW       *
*   LBLX        : NUMBER OF UNKNOWNS PER X-GROUP. DIMENSION LONX       *
*   LBLY        : NUMBER OF UNKNOWNS PER Y-GROUP. DIMENSION LONY       *
*   LBLZ        : NUMBER OF UNKNOWNS PER Z-GROUP. DIMENSION LONZ       *
*   MUVW        : INDICES USED WITH W-DIRECTED COMPRESSED DIAGONAL     *
*                 STORAGE MODE MATRICES IN VECTOR MODE. DIMENSION LL4W *
*   MUVX        : INDICES USED WITH X-DIRECTED COMPRESSED DIAGONAL     *
*                 STORAGE MODE MATRICES IN VECTOR MODE. DIMENSION LL4X *
*   MUVY        : INDICES USED WITH Y-DIRECTED COMPRESSED DIAGONAL     *
*                 STORAGE MODE MATRICES IN VECTOR MODE. DIMENSION LL4Y *
*   MUVZ        : INDICES USED WITH Z-DIRECTED COMPRESSED DIAGONAL     *
*                 STORAGE MODE MATRICES IN VECTOR MODE. DIMENSION LL4Z *
*   IPVW        : W-DIRECTED VECTOR PERMUTATION MATRIX. DIMENSION LL4  *
*   IPVX        : X-DIRECTED VECTOR PERMUTATION MATRIX. DIMENSION LL4  *
*   IPVY        : Y-DIRECTED VECTOR PERMUTATION MATRIX. DIMENSION LL4  *
*   IPVZ        : Z-DIRECTED VECTOR PERMUTATION MATRIX. DIMENSION LL4  *
*                                                                      *
* INFORMATION RELATED TO CYLINDRICAL CORRECTIONS IN CARTESIAN GEOMETRY *
*   NR0         : NUMBER OF RADII.                                     *
*   RR0         : RADII. DIMENSION NR0                                 *
*   XR0         : COORDINATES ON PRINCIPAL AXIS. DIMENSION NR0         *
*   ANG         : ANGLES FOR APPLYING CIRCULAR CORRECTION.             *
*                 DIMENSION NR0                                        *
*                                                                      *
************************************************************************
*
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MAT(MAXPTS),IDL(MAXPTS),VOL(MAXPTS))
*
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      ITYPE=ISTATE(1)
*
      IF(IMPX.GE.1) WRITE (6,'(/35H TRITRK: DEGREE OF FINITE ELEMENT I,
     1 6HELEM =,I3/9X,25HTYPE OF QUADRATURE ICOL =,I3/9X,10HTYPE OF DI,
     2 19HSCRETIZATION ICHX =,I3/)') IELEM,ICOL,ICHX
      IF((IMPX.GE.1).AND.(ISEG.GT.0)) WRITE (6,'(18H TRITRK: SUPERVECT,
     1 27HORIZATION OPTION ON. ISEG =,I4,8H  IMPV =,I3/)') ISEG,IMPV
      IF(ISTATE(9).EQ.0) THEN
         IF((ITYPE.NE.1).AND.(ITYPE.NE.2).AND.(ITYPE.NE.3).AND.
     1      (ITYPE.NE.5).AND.(ITYPE.NE.6).AND.(ITYPE.NE.7).AND.
     2      (ITYPE.NE.8).AND.(ITYPE.NE.9)) THEN
            CALL XABORT('TRITRK: DISCRETIZATION NOT AVAILABLE.')
         ENDIF
         ALLOCATE(XXX(MAXPTS+1),YYY(MAXPTS+1),ZZZ(MAXPTS+1))
*
         ALLOCATE(ISPLX(MAXPTS),ISPLY(MAXPTS),ISPLZ(MAXPTS))
         CALL READ3D(MAXPTS,MAXPTS,MAXPTS,MAXPTS,IPGEOM,IHEX,IR,ILK,
     1   SIDE,XXX,YYY,ZZZ,IMPX,LX,LY,LZ,MAT,NEL,NCODE,ICODE,ZCODE,
     2   ISPLX,ISPLY,ISPLZ,ISPLH,ISPLL)
         DEALLOCATE(ISPLX,ISPLY,ISPLZ)
         IF((ITYPE.GE.8).AND.(ICHX.EQ.2)) THEN
           IF(ISPLL.EQ.0) THEN
             CALL XABORT('TRITRK: SPLITL KEYWORD MISSING IN GEOMETRY.')
           ENDIF
           ISPLH=ISPLL
         ELSE IF(ITYPE.GE.8) THEN
           ISPLH=ISPLH+1  
         ENDIF
      ELSE
         CALL XABORT('TRITRK: DISCRETIZATION NOT AVAILABLE.')
      ENDIF
*----
*  UNFOLD HEXAGONAL GEOMETRY CASES.
*----
      CHEX=(ITYPE.EQ.8).OR.(ITYPE.EQ.9)
      IF(CHEX.AND.(IHEX.NE.9)) THEN
         ALLOCATE(IDP(MAXPTS),IMX(NEL))
         DO 30 I=1,NEL
         IMX(I)=MAT(I)
   30    CONTINUE
         LXOLD=LX
         CALL BIVALL(MAXPTS,IHEX,LXOLD,LX,IDP)
         DO 41 KZ=1,LZ
         DO 40 KX=1,LX
         KEL=IDP(KX)+(KZ-1)*LXOLD
         MAT(KX+(KZ-1)*LX)=IMX(KEL)
   40    CONTINUE
   41    CONTINUE
         DEALLOCATE(IMX,IDP)
         NEL=LX*LZ
      ENDIF
*----
*  PROCESS INFORMATION RELATED TO CYLINDRICAL CORRECTION IN CARTESIAN
*  GEOMETRIES.
*----
      CALL LCMLEN(IPGEOM,'RR0',NR0,ITYLCM)
      IF(NR0.GT.0) THEN
         IF((ITYPE.NE.5).AND.(ITYPE.NE.7)) CALL XABORT('TRITRK: CYLIND'
     1   //'RICAL CORRECTIONS ARE LIMITED TO CARTESIAN GEOMETRIES.')
         IF(IMPX.GT.0) WRITE(6,'(/33H TRITRK: PERFORM A CYLINDRICAL CO,
     2   35HRRECTION ON THE CARTESIAN BOUNDARY.)')
         ALLOCATE(RR0(NR0),XR0(NR0),ANG(NR0))
         CALL LCMGET(IPGEOM,'RR0',RR0)
         CALL LCMGET(IPGEOM,'XR0',XR0)
         CALL LCMGET(IPGEOM,'ANG',ANG)
         CALL LCMPUT(IPTRK,'RR0',NR0,2,RR0)
         CALL LCMPUT(IPTRK,'XR0',NR0,2,XR0)
         CALL LCMPUT(IPTRK,'ANG',NR0,2,ANG)
         DEALLOCATE(ANG,XR0,RR0)
      ENDIF
*
      IF(LX*LY*LZ.GT.MAXPTS) THEN
         WRITE (HSMG,'(39HTRITRK: MAXPTS SHOULD BE INCREASED FROM,I8,
     1   3H TO,I8)') MAXPTS,LX*LY*LZ
         CALL XABORT(HSMG)
      ENDIF
*----
*  1-D AND 2-D CASES.
*----
      IDIM=1
      IF((ITYPE.EQ.5).OR.(ITYPE.EQ.6).OR.(ITYPE.EQ.8)) IDIM=2
      IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) IDIM=3
      IF((NCODE(3).EQ.0).AND.(NCODE(4).EQ.0).AND.(.NOT.CHEX)) THEN
         IF((IDIM.NE.1).OR.(LY.NE.1)) CALL XABORT('TRITRK: INVALID 1D '
     1   //'GEOMETRY.')
         NCODE(3)=2
         NCODE(4)=5
         ZCODE(3)=1.0
         ZCODE(4)=1.0
         YYY(1)=0.0
         YYY(2)=2.0
      ENDIF
      IF((NCODE(5).EQ.0).AND.(NCODE(6).EQ.0)) THEN
         IF((IDIM.EQ.3).OR.(LZ.NE.1)) CALL XABORT('TRITRK: INVALID 1D '
     1   //'OR 2D GEOMETRY.')
         NCODE(5)=2
         NCODE(6)=5
         ZCODE(5)=1.0
         ZCODE(6)=1.0
         ZZZ(1)=0.0
         ZZZ(2)=2.0
      ENDIF
*----
*  2-D CYLINDRICAL CASES.
*----
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
      IF(ITYPE.EQ.6) THEN
         LY=LZ
         DO 45 I=1,LZ+1
         YYY(I)=ZZZ(I)
   45    CONTINUE
         NCODE(3)=NCODE(5)
         NCODE(4)=NCODE(6)
         ICODE(3)=ICODE(5)
         ICODE(4)=ICODE(6)
         ZCODE(3)=ZCODE(5)
         ZCODE(4)=ZCODE(6)
         NCODE(5)=0
         NCODE(6)=0
         ZCODE(5)=0.0
         ZCODE(6)=0.0
      ENDIF
*----
*  UNFOLD THE DOMAIN IN DIAGONAL SYMMETRY CASES.
*----
      IDIAG=0
      IF((NCODE(2).EQ.3).AND.(NCODE(3).EQ.3)) THEN
         IDIAG=1
         NCODE(3)=NCODE(1)
         NCODE(2)=NCODE(4)
         ICODE(3)=ICODE(1)
         ICODE(2)=ICODE(4)
         ZCODE(3)=ZCODE(1)
         ZCODE(2)=ZCODE(4)
         K=NEL
         DO 82 IZ=LZ,1,-1
         IOFF=(IZ-1)*LX*LY
         DO 81 IY=LY,1,-1
         DO 70 IX=LX,IY+1,-1
         MAT(IOFF+(IY-1)*LX+IX)=MAT(IOFF+(IX-1)*LY+IY)
   70    CONTINUE
         DO 80 IX=IY,1,-1
         MAT(IOFF+(IY-1)*LX+IX)=MAT(K)
         K=K-1
   80    CONTINUE
   81    CONTINUE
   82    CONTINUE
         NEL=LX*LY*LZ
         IF(K.NE.0) THEN
            CALL XABORT('TRITRK: UNABLE TO UNFOLD THE DOMAIN(1).')
         ENDIF
      ELSE IF((NCODE(1).EQ.3).AND.(NCODE(4).EQ.3)) THEN
         IDIAG=1
         NCODE(1)=NCODE(3)
         NCODE(4)=NCODE(2)
         ICODE(1)=ICODE(3)
         ICODE(4)=ICODE(2)
         ZCODE(1)=ZCODE(3)
         ZCODE(4)=ZCODE(2)
         K=NEL
         DO 92 IZ=LZ,1,-1
         IOFF=(IZ-1)*LX*LY
         DO 91 IY=LY,1,-1
         DO 90 IX=LX,IY,-1
         MAT(IOFF+(IY-1)*LX+IX)=MAT(K)
         K=K-1
   90    CONTINUE
   91    CONTINUE
   92    CONTINUE
         DO 102 IZ=1,LZ
         IOFF=(IZ-1)*LX*LY
         DO 101 IY=1,LY
         DO 100 IX=1,IY-1
         MAT(IOFF+(IY-1)*LX+IX)=MAT(IOFF+(IX-1)*LY+IY)
  100    CONTINUE
  101    CONTINUE
  102    CONTINUE
         NEL=LX*LY*LZ
         IF(K.NE.0) THEN
            CALL XABORT('TRITRK: UNABLE TO UNFOLD THE DOMAIN(2).')
         ENDIF
      ENDIF
      IF(IMPX.GT.5) THEN
         WRITE(6,600) 'NCODE',(NCODE(I),I=1,6)
         WRITE(6,600) 'MAT',(MAT(I),I=1,LX*LY*LZ)
      ENDIF
*
      CALL KDRCPU(TK1)
      MAXQF=6*NEL
      IF(CHEX) MAXQF=8*NEL
      IF((ICHX.EQ.1).AND.(.NOT.CHEX)) THEN
         MAXKN=NEL*(IELEM+1)**3
      ELSE IF((ICHX.EQ.2).AND.(.NOT.CHEX)) THEN
         MAXKN=NEL*(1+6*IELEM**2)
      ELSE IF((ICHX.EQ.3).AND.(.NOT.CHEX)) THEN
         MAXKN=6*NEL
      ELSE IF((ICHX.EQ.1).AND.CHEX) THEN
         IF(ISPLH.EQ.1) THEN
            MAXKN=12*NEL
         ELSE
            MAXKN=2*(1+ISPLH*(ISPLH-1)*3)*NEL
         ENDIF
      ELSE IF((ICHX.EQ.2).AND.CHEX) THEN
         MAXKN=(NEL*ISPLH**2)*(3+6*IELEM*IELEM*(IELEM+2))
         MAXQF=(NEL*ISPLH**2)*8
      ELSE IF((ICHX.EQ.3).AND.CHEX) THEN
         IF(ISPLH.EQ.1) THEN
            MAXKN=8*NEL
         ELSE
            MAXKN=(18*(ISPLH-1)**2+8)*NEL
         ENDIF
      ELSE
         CALL XABORT('TRITRK: INVALID TYPE OF DISCRETIZATION.')
      ENDIF
      IF(CYLIND) THEN
         MAXDD=NEL
      ELSE
         MAXDD=1
      ENDIF
      IF((ICHX.NE.2).AND.CHEX.AND.(IELEM.NE.1)) CALL XABORT('TRITRK: T'
     1 //'HIS HEXAGONAL DISCRETIZATIONS IS LIMITED TO LINEAR ORDER.')
      IF(CHEX.AND.(NCODE(1).EQ.5)) CALL XABORT('TRITRK: SYME BOUNDARY '
     1 //'CONDITION IS NOT AVAILABLE AROUND THE HEXAGONAL PLANE.')
      ALLOCATE(XX(NEL),YY(NEL),ZZ(NEL),DD(MAXDD),KN(MAXKN),QFR(MAXQF),
     1 IQFR(MAXQF))
      CALL XDISET(KN,MAXKN,0)
      CALL XDRSET(QFR,MAXQF,0.0)
      CALL XDISET(IQFR,MAXQF,0)
      LL4=0
      IF((ICHX.EQ.1).AND.(.NOT.CHEX)) THEN
         CALL TRIPKN(IELEM,LX,LY,LZ,LL4,CYLIND,XXX,YYY,ZZZ,XX,YY,ZZ,DD,
     1   KN,QFR,IQFR,VOL,MAT,NCODE,ICODE,ZCODE,IMPX)
         IF((IMPX.GT.0).AND.(IELEM.EQ.1)) THEN
            WRITE (6,'(/40H TRITRK: MESH CORNER FINITE DIFFERENCES.)')
         ENDIF
         LL4W=0
         LL4X=LL4
         LL4Y=LL4
         LL4Z=LL4
      ELSE IF((ICHX.EQ.2).AND.(.NOT.CHEX)) THEN
         CALL TRIDKN(IMPX,LX,LY,LZ,CYLIND,IELEM,LL4,LL4F,LL4X,LL4Y,
     1   LL4Z,NCODE,ICODE,ZCODE,MAT,VOL,XXX,YYY,ZZZ,XX,YY,ZZ,DD,KN,
     2   QFR,IQFR,IDL)
         MAXIP=LX*LY*LZ
         NUN=LL4
      ELSE IF((ICHX.EQ.3).AND.(.NOT.CHEX)) THEN
         MAXIP=LX*LY*LZ
         CALL TRIDFC(IMPX,LX,LY,LZ,CYLIND,NCODE,ICODE,ZCODE,MAT,XXX,
     1   YYY,ZZZ,LL0,VOL,XX,YY,ZZ,DD,KN,QFR,IQFR)
         IF(IELEM.EQ.1) THEN
            LL4=LL0
            IF(IMPX.GT.0) WRITE (6,'(/29H TRITRK: MESH CENTERED FINITE,
     1      13H DIFFERENCES.)')
         ELSE IF((IELEM.GT.1).AND.(ICHX.EQ.3)) THEN
            LL4=LL0*IELEM**IDIM
            IF(IMPX.GT.0) WRITE (6,'(/29H TRITRK: NODAL COLLOCATION ME,
     1      13HTHOD OF ORDER,I3,1H.)') IELEM
         ENDIF
*        COMPUTE INDICES IDL.
         IF(ICHX.EQ.3) THEN
*           NODAL COLLOCATION METHOD.
            NUN=0
            DO 110 K=1,NEL
            IDL(K)=0
            IF(MAT(K).EQ.0) GO TO 110
            NUN=NUN+1
            IDL(K)=1+IELEM*(NUN-1)
  110       CONTINUE
            NUN=LL4
         ENDIF
         LL4W=0
         LL4X=LL4
         LL4Y=LL4
         LL4Z=LL4
      ELSE IF((ICHX.EQ.1).AND.CHEX) THEN
         MAXIP=1
         IF(IELEM.NE.1) CALL XABORT('TRITRK: INVALID DISCRETIZATION.')
         CALL TRIPRH(ISPLH,IPTRK,LX,LZ,LL4,SIDE,ZZZ,ZZ,KN,QFR,IQFR,VOL,
     1   MAT,NCODE,ICODE,ZCODE,IMPX)
         IF(IMPX.GT.0) WRITE (6,'(/32H TRITRK: MESH CORNER FINITE DIFF,
     1   39HERENCES FOR HEXAGONAL GEOMETRY. ISPLH =,I3,1H.)') ISPLH
         LL4W=LL4
         LL4X=LL4
         LL4Y=LL4
         LL4Z=LL4
      ELSE IF((ICHX.EQ.2).AND.CHEX) THEN
         NEL=LX*LZ
         LXH=LX/(3*ISPLH**2)
         NBLOS=LXH*LZ*ISPLH**2
         NBC=INT((SQRT(REAL((4*LXH-1)/3))+1.)/2.)
         MAXIP=3*(2*LXH*ISPLH*IELEM+2*NBC-1)*ISPLH*LZ*IELEM**2
     1        +3*LXH*(LZ+1)*(ISPLH**2)*IELEM**2
         ALLOCATE(IPERT(NBLOS),FRZ(NBLOS))
         CALL TRISFH(IMPX,MAXKN,MAXIP,NBLOS,ISPLH,IELEM,LXH,LZ,MAT,SIDE,
     1   ZZZ,NCODE,ICODE,ZCODE,LL4,LL4F,LL4W,LL4X,LL4Y,LL4Z,VOL,IDL,
     2   IPERT,ZZ,FRZ,KN,QFR,IQFR)
         CALL LCMPUT(IPTRK,'IPERT',NBLOS,1,IPERT)
         CALL LCMPUT(IPTRK,'FRZ',NBLOS,2,FRZ)
         DEALLOCATE(FRZ,IPERT)
         NUN=LL4
         IF(IMPX.GT.0) WRITE (6,'(/32H TRITRK: THOMAS-RAVIART-SCHNEIDE,
     1   49HR FINITE ELEMENTS FOR HEXAGONAL GEOMETRY. ISPLH =,I3,1H.)')
     2   ISPLH
      ELSE IF((ICHX.EQ.3).AND.CHEX) THEN
         MAXIP=LX*LZ
         IF(IELEM.NE.1) CALL XABORT('TRITRK: INVALID DISCRETIZATION.')
         CALL TRIDFH(ISPLH,IPTRK,IDIM,LX,LZ,LL4,NUN,SIDE,ZZZ,ZZ,KN,QFR,
     1   IQFR,VOL,MAT,IDL,NCODE,ICODE,ZCODE,IMPX)
         IF(IMPX.GT.0) WRITE (6,'(/32H TRITRK: MESH CENTERED FINITE DI,
     1   41HFFERENCES FOR HEXAGONAL GEOMETRY. ISPLH =,I3,1H.)') ISPLH
         LL4W=LL4
         LL4X=LL4
         LL4Y=LL4
         LL4Z=LL4
      ENDIF
*----
*  APPEND THE PN FLUXES AT THE END OF UNKNOWN VECTOR.
*----
      IF(NLF.GE.2) THEN
         IF((ITYPE.EQ.2).OR.((ITYPE.EQ.5).AND.(ISPN.EQ.1)).OR.
     1                      ((ITYPE.EQ.7).AND.(ISPN.EQ.1))) THEN
            NUN=LL4+LL4*(NLF-2)/2
         ELSE IF((ITYPE.EQ.8).AND.(ISPN.EQ.1)) THEN
            NUN=NUN+NUN*(NLF-2)/2
         ELSE IF((ITYPE.EQ.9).AND.(ISPN.EQ.1)) THEN
            NUN=NUN+NUN*(NLF-2)/2
         ELSE
            CALL XABORT('TRITRK: GEOMETRY NOT SUPPORTED WITH PN.')
         ENDIF
      ENDIF
*----
*  COMPUTE INDICES IDL FOR PRIMAL FINITE ELEMENTS.
*----
      IF(ICHX.EQ.1) THEN
         NUN=LL4
         DO 130 K=1,NEL
         IF(MAT(K).EQ.0) THEN
            IDL(K)=0
         ELSE
            NUN=NUN+1
            IDL(K)=NUN
         ENDIF
  130    CONTINUE
      ENDIF
*
      IF(IMPX.GT.0) WRITE (6,'(/34H TRITRK: ORDER OF LINEAR SYSTEMS =,
     1 I8/9X,37HNUMBER OF UNKNOWNS PER ENERGY GROUP =,I8)') LL4,NUN
      DEALLOCATE(ZZZ,YYY,XXX)
      CALL KDRCPU(TK2)
      IF(IMPX.GE.2) WRITE(6,'(/37H TRITRK: CPU TIME FOR FINITE ELEMENT ,
     1 11HNUMBERING =,F7.2,2H S)') TK2-TK1
*----
*  COMPUTE INDICES MUW, MUX, MUY, MUZ, IPW, IPX, IPY AND IPZ.
*----
      CALL KDRCPU(TK1)
      IF(CHEX) ALLOCATE(MUW(LL4))
      ALLOCATE(MUX(LL4),MUY(LL4),MUZ(LL4))
      IF(CHEX) ALLOCATE(IPW(LL4))
      IF(ICHX.NE.2) THEN
         ALLOCATE(IPX(LL4),IPY(LL4),IPZ(LL4))
         DO 140 I=1,LL4
         IPX(I)=I
  140    CONTINUE
      ENDIF
*
      IF((ICHX.EQ.1).AND.(.NOT.CHEX)) THEN
         CALL BIVCOL(IPTRK,IMPX,IELEM,2)
         CALL TRICHP(IELEM,LX,LY,LZ,LL4,MAT,KN,MUX,MUY,MUZ,IPY,IPZ,IMPX)
      ELSE IF((ICHX.EQ.2).AND.(.NOT.CHEX)) THEN
         LL4W=0
         CALL BIVCOL(IPTRK,IMPX,IELEM,ICOL)
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         ALLOCATE(V((IELEM+1),IELEM))
         CALL LCMGET(IPTRK,'V',V)
         CALL LCMSIX(IPTRK,' ',2)
         ALLOCATE(IPBBX(2*IELEM*LL4X),IPBBY(2*IELEM*LL4Y),
     1   IPBBZ(2*IELEM*LL4Z))
         ALLOCATE(BBX(2*IELEM*LL4X),BBY(2*IELEM*LL4Y),BBZ(2*IELEM*LL4Z))
         CALL TRICHD(IMPX,LX,LY,LZ,CYLIND,IELEM,LL4,LL4F,LL4X,LL4Y,LL4Z,
     1   MAT,VOL,XX,YY,ZZ,DD,KN,V,MUX,MUY,MUZ,IPBBX,IPBBY,IPBBZ,BBX,BBY,
     2   BBZ)
         IF(LL4X.GT.0) THEN
            CALL LCMPUT(IPTRK,'IPBBX',2*IELEM*LL4X,1,IPBBX)
            CALL LCMPUT(IPTRK,'XB',2*IELEM*LL4X,2,BBX)
         ENDIF
         IF(LL4Y.GT.0) THEN
            CALL LCMPUT(IPTRK,'IPBBY',2*IELEM*LL4Y,1,IPBBY)
            CALL LCMPUT(IPTRK,'YB',2*IELEM*LL4Y,2,BBY)
         ENDIF
         IF(LL4Z.GT.0) THEN
            CALL LCMPUT(IPTRK,'IPBBZ',2*IELEM*LL4Z,1,IPBBZ)
            CALL LCMPUT(IPTRK,'ZB',2*IELEM*LL4Z,2,BBZ)
         ENDIF
         DEALLOCATE(BBZ,BBY,BBX,IPBBZ,IPBBY,IPBBX)
         DEALLOCATE(V)
      ELSE IF((ICHX.EQ.3).AND.(.NOT.CHEX)) THEN
         CALL TRICH1(IELEM,IDIM,LX,LY,LZ,LL4,MAT,KN,MUX,MUY,MUZ,IPY,
     1   IPZ,IMPX)
      ELSE IF((ICHX.EQ.1).AND.CHEX) THEN
         CALL BIVCOL(IPTRK,IMPX,IELEM,2)
         CALL TRICH3(ISPLH,IPTRK,LX,LZ,LL4,MAT,KN,MUW,MUX,MUY,MUZ,IPW,
     1   IPX,IPY,IPZ,IMPX)
      ELSE IF((ICHX.EQ.2).AND.CHEX) THEN
         LXH=LX/(3*ISPLH**2)
         NBLOS=LXH*LZ*ISPLH**2
         ALLOCATE(IPERT(NBLOS),FRZ(NBLOS))
         CALL LCMGET(IPTRK,'IPERT',IPERT)
         CALL LCMGET(IPTRK,'FRZ',FRZ)
         CALL BIVCOL(IPTRK,IMPX,IELEM,ICOL)
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         ALLOCATE(V((IELEM+1),IELEM),H((IELEM+1),IELEM))
         CALL LCMGET(IPTRK,'V',V)
         CALL LCMGET(IPTRK,'H',H)
         CALL LCMSIX(IPTRK,' ',2)
         ALLOCATE(IPBBW(2*IELEM*LL4W),IPBBX(2*IELEM*LL4X),
     1   IPBBY(2*IELEM*LL4Y),IPBBZ(2*IELEM*LL4Z))
         ALLOCATE(BBW(2*IELEM*LL4W),BBX(2*IELEM*LL4X),
     1   BBY(2*IELEM*LL4Y),BBZ(2*IELEM*LL4Z))
         ALLOCATE(CTRAN(((IELEM+1)*IELEM)**2))
         CALL TRICHH(IMPX,MAXKN,NBLOS,LXH,LZ,IELEM,ISPLH,LL4,LL4F,LL4W,
     1   LL4X,LL4Y,LL4Z,SIDE,ZZ,FRZ,IPERT,KN,V,H,MUW,MUX,MUY,MUZ,IPBBW,
     2   IPBBX,IPBBY,IPBBZ,BBW,BBX,BBY,BBZ,CTRAN)
         CALL LCMPUT(IPTRK,'CTRAN',((IELEM+1)*IELEM)**2,4,CTRAN)
         CALL LCMPUT(IPTRK,'IPBBW',2*IELEM*LL4W,1,IPBBW)
         CALL LCMPUT(IPTRK,'WB',2*IELEM*LL4W,2,BBW)
         CALL LCMPUT(IPTRK,'IPBBX',2*IELEM*LL4X,1,IPBBX)
         CALL LCMPUT(IPTRK,'XB',2*IELEM*LL4X,2,BBX)
         CALL LCMPUT(IPTRK,'IPBBY',2*IELEM*LL4Y,1,IPBBY)
         CALL LCMPUT(IPTRK,'YB',2*IELEM*LL4Y,2,BBY)
         IF(LL4Z.GT.0) THEN
            CALL LCMPUT(IPTRK,'IPBBZ',2*IELEM*LL4Z,1,IPBBZ)
            CALL LCMPUT(IPTRK,'ZB',2*IELEM*LL4Z,2,BBZ)
         ENDIF
         DEALLOCATE(BBZ,BBY,BBX,BBW,IPBBZ,IPBBY,IPBBX,IPBBW)
         DEALLOCATE(H,V,CTRAN,FRZ,IPERT)
      ELSE IF((ICHX.EQ.3).AND.CHEX) THEN
         CALL TRICH4(ISPLH,IPTRK,IDIM,LX,LZ,LL4,MAT,KN,MUW,MUX,MUY,MUZ,
     1   IPW,IPX,IPY,IPZ,IMPX)
      ENDIF
      CALL KDRCPU(TK2)
      IF(IMPX.GE.2) WRITE(6,'(/36H TRITRK: CPU TIME FOR ADI SPLITTING ,
     1 11HNUMBERING =,F7.2,2H S)') TK2-TK1
      IF(IMPX.GT.5) THEN
         I1=1
         DO 150 I=1,(NEL-1)/8+1
         I2=I1+7
         IF(I2.GT.NEL) I2=NEL
         WRITE (6,620) (J,J=I1,I2)
         WRITE (6,630) (MAT(J),J=I1,I2)
         WRITE (6,640) (IDL(J),J=I1,I2)
         WRITE (6,650) (VOL(J),J=I1,I2)
         I1=I1+8
  150    CONTINUE
      ENDIF
*----
*  SUPERVECTORIZATION CONTROL.
*----
      LTSW=0
      IF(ISEG.GT.0) THEN
         CALL KDRCPU(TK1)
         ALLOCATE(ISET(LL4))
         IF(CHEX) THEN
            ISET(1)=0
            K1=MUW(1)+1
            DO 160 I=2,LL4W
            ISET(I)=0
            K2=MUW(I)
            DO 155 J=I-K2+K1,I-1
            ISET(J)=1
  155       CONTINUE
            K1=K2+1
  160       CONTINUE
            NSYS=0
            DO 165 I=1,LL4W
            IF(ISET(I).EQ.0) NSYS=NSYS+1
  165       CONTINUE
            LONW=1+(NSYS-1)/ISEG
            ALLOCATE(NBLW(LONW),LBLW(LONW),MUVW(LONW),IPVW(LONW))
            CALL VECPER('W',IMPV,ISEG,LL4W,MUW,LONW,LTSW2,NBLW,LBLW,
     1      MUVW,IPVW)
            IMU=0
            DO 166 I=1,LONW
            IMU=IMU+LBLW(I)
  166       CONTINUE
            LTSW=MAX(LTSW,LTSW2)
            CALL LCMPUT(IPTRK,'NBLW',LONW,1,NBLW)
            CALL LCMPUT(IPTRK,'LBLW',LONW,1,LBLW)
            CALL LCMPUT(IPTRK,'MUVW',IMU,1,MUVW)
            CALL LCMPUT(IPTRK,'IPVW',LL4W,1,IPVW)
            DEALLOCATE(IPVW,MUVW,LBLW,NBLW)
            IMU=IMU*ISEG
            CALL LCMPUT(IPTRK,'LL4VW',1,1,IMU)
         ENDIF
         IF(IDIAG.EQ.0) THEN
            ISET(1)=0
            K1=MUX(1)+1
            DO 175 I=2,LL4X
            ISET(I)=0
            K2=MUX(I)
            DO 170 J=I-K2+K1,I-1
            ISET(J)=1
  170       CONTINUE
            K1=K2+1
  175       CONTINUE
            NSYS=0
            DO 180 I=1,LL4X
            IF(ISET(I).EQ.0) NSYS=NSYS+1
  180       CONTINUE
            LONX=1+(NSYS-1)/ISEG
            ALLOCATE(NBLX(LONX),LBLX(LONX),MUVX(LONX),IPVX(LONX))
            CALL VECPER('X',IMPV,ISEG,LL4X,MUX,LONX,LTSW2,NBLX,LBLX,
     1      MUVX,IPVX)
            IMU=0
            DO 185 I=1,LONX
            IMU=IMU+LBLX(I)
  185       CONTINUE
            LTSW=MAX(LTSW,LTSW2)
            CALL LCMPUT(IPTRK,'NBLX',LONX,1,NBLX)
            CALL LCMPUT(IPTRK,'LBLX',LONX,1,LBLX)
            CALL LCMPUT(IPTRK,'MUVX',IMU,1,MUVX)
            CALL LCMPUT(IPTRK,'IPVX',LL4X,1,IPVX)
            DEALLOCATE(IPVX,MUVX,LBLX,NBLX)
            IMU=IMU*ISEG
            CALL LCMPUT(IPTRK,'LL4VX',1,1,IMU)
         ENDIF
         IF(IDIM.GE.2) THEN
            ISET(1)=0
            K1=MUY(1)+1
            DO 200 I=2,LL4Y
            ISET(I)=0
            K2=MUY(I)
            DO 190 J=I-K2+K1,I-1
            ISET(J)=1
  190       CONTINUE
            K1=K2+1
  200       CONTINUE
            NSYS=0
            DO 210 I=1,LL4Y
            IF(ISET(I).EQ.0) NSYS=NSYS+1
  210       CONTINUE
            LONY=1+(NSYS-1)/ISEG
            ALLOCATE(NBLY(LONY),LBLY(LONY),MUVY(LONY),IPVY(LONY))
            CALL VECPER('Y',IMPV,ISEG,LL4Y,MUY,LONY,LTSW2,NBLY,LBLY,
     1      MUVY,IPVY)
            IMU=0
            DO 215 I=1,LONY
            IMU=IMU+LBLY(I)
  215       CONTINUE
            LTSW=MAX(LTSW,LTSW2)
            CALL LCMPUT(IPTRK,'NBLY',LONY,1,NBLY)
            CALL LCMPUT(IPTRK,'LBLY',LONY,1,LBLY)
            CALL LCMPUT(IPTRK,'MUVY',IMU,1,MUVY)
            CALL LCMPUT(IPTRK,'IPVY',LL4Y,1,IPVY)
            DEALLOCATE(IPVY,MUVY,LBLY,NBLY)
            IMU=IMU*ISEG
            CALL LCMPUT(IPTRK,'LL4VY',1,1,IMU)
         ENDIF
         IF(IDIM.EQ.3) THEN
            ISET(1)=0
            K1=MUZ(1)+1
            DO 230 I=2,LL4Z
            ISET(I)=0
            K2=MUZ(I)
            DO 220 J=I-K2+K1,I-1
            ISET(J)=1
  220       CONTINUE
            K1=K2+1
  230       CONTINUE
            NSYS=0
            DO 240 I=1,LL4Z
            IF(ISET(I).EQ.0) NSYS=NSYS+1
  240       CONTINUE
            LONZ=1+(NSYS-1)/ISEG
            ALLOCATE(NBLZ(LONZ),LBLZ(LONZ),MUVZ(LONZ),IPVZ(LONZ))
            CALL VECPER('Z',IMPV,ISEG,LL4Z,MUZ,LONZ,LTSW2,NBLZ,LBLZ,
     1      MUVZ,IPVZ)
            IMU=0
            DO 250 I=1,LONZ
            IMU=IMU+LBLZ(I)
  250       CONTINUE
            LTSW=MAX(LTSW,LTSW2)
            CALL LCMPUT(IPTRK,'NBLZ',LONZ,1,NBLZ)
            CALL LCMPUT(IPTRK,'LBLZ',LONZ,1,LBLZ)
            CALL LCMPUT(IPTRK,'MUVZ',IMU,1,MUVZ)
            CALL LCMPUT(IPTRK,'IPVZ',LL4Z,1,IPVZ)
            DEALLOCATE(IPVZ,MUVZ,LBLZ,NBLZ)
            IMU=IMU*ISEG
            CALL LCMPUT(IPTRK,'LL4VZ',1,1,IMU)
         ENDIF
         DEALLOCATE(ISET)
         CALL KDRCPU(TK2)
         IF(IMPX.GE.2) WRITE(6,'(/33H TRITRK: CPU TIME FOR SUPERVECTOR,
     1   19HIZATION NUMBERING =,F7.2,2H S)') TK2-TK1
      ENDIF
*----
*  SAVE STATE-VECTOR AND TRACKING INFORMATION.
*----
      CALL XDISET(IGP,NSTATE,0)
      IGP(1)=NEL
      IGP(2)=NUN
      IF(ILK) THEN
         IGP(3)=0
      ELSE
         IGP(3)=1
      ENDIF
      IGP(4)=ISTATE(7)
      IGP(5)=0
      IGP(6)=ITYPE
      IGP(7)=IHEX
      IGP(8)=IDIAG
      IGP(9)=IELEM
      IGP(10)=ICOL
      IGP(11)=LL4
      IGP(12)=ICHX
      IGP(13)=ISPLH
      IGP(14)=LX
      IGP(15)=LY
      IGP(16)=LZ
      IGP(17)=ISEG
      IF(ISEG.NE.0) THEN
         IGP(18)=IMPV
         IGP(19)=LTSW
         IGP(20)=LONW
         IGP(21)=LONX
         IGP(22)=LONY
         IGP(23)=LONZ
      ENDIF
      IGP(24)=NR0
      IF(ICHX.EQ.2) THEN
         IGP(25)=LL4F
         IGP(26)=LL4W
         IGP(27)=LL4X
         IGP(28)=LL4Y
         IGP(29)=LL4Z
      ENDIF
      IGP(30)=NLF
      IGP(31)=ISPN
      IGP(32)=ISCAT
      IGP(33)=NADI
      IGP(34)=NVD
      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,IGP)
      CALL LCMPUT(IPTRK,'MATCOD',NEL,1,MAT)
      CALL LCMPUT(IPTRK,'VOLUME',NEL,2,VOL)
      CALL LCMPUT(IPTRK,'KEYFLX',NEL,1,IDL)
      CALL LCMPUT(IPTRK,'NCODE',6,1,NCODE)
      CALL LCMPUT(IPTRK,'ZCODE',6,2,ZCODE)
      CALL LCMPUT(IPTRK,'ICODE',6,1,ICODE)
      CALL LCMPUT(IPTRK,'ZZ',NEL,2,ZZ)
      CALL LCMPUT(IPTRK,'KN',MAXKN,1,KN)
      CALL LCMPUT(IPTRK,'QFR',MAXQF,2,QFR)
      CALL LCMPUT(IPTRK,'IQFR',MAXQF,1,IQFR)
      IF(ICHX.NE.2) THEN
         CALL LCMPUT(IPTRK,'IPX',LL4,1,IPX)
         DEALLOCATE(IPX)
      ENDIF
      IF(CHEX) THEN
         CALL LCMPUT(IPTRK,'SIDE',1,2,SIDE)
         CALL LCMPUT(IPTRK,'MUW',LL4W,1,MUW)
         IF(ICHX.NE.2) THEN
            CALL LCMPUT(IPTRK,'IPW',LL4,1,IPW)
            DEALLOCATE(IPW)
         ENDIF
         DEALLOCATE(MUW)
      ELSE
         CALL LCMPUT(IPTRK,'XX',NEL,2,XX)
         CALL LCMPUT(IPTRK,'YY',NEL,2,YY)
         IF(.NOT.CYLIND) DD=0.0
         CALL LCMPUT(IPTRK,'DD',MAXDD,2,DD)
      ENDIF
      DEALLOCATE(XX,YY,ZZ,DD,KN,QFR,IQFR)
      IF((IDIAG.EQ.0).AND.(LL4X.GT.0)) THEN
         CALL LCMPUT(IPTRK,'MUX',LL4X,1,MUX)
      ENDIF
      IF((IDIM.GE.2).AND.(LL4Y.GT.0)) THEN
         CALL LCMPUT(IPTRK,'MUY',LL4Y,1,MUY)
         IF(ICHX.NE.2) THEN
            CALL LCMPUT(IPTRK,'IPY',LL4,1,IPY)
            DEALLOCATE(IPY)
         ENDIF
      ELSE
         IF(ICHX.NE.2) DEALLOCATE(IPY)
      ENDIF
      IF((IDIM.EQ.3).AND.(LL4Z.GT.0)) THEN
         CALL LCMPUT(IPTRK,'MUZ',LL4Z,1,MUZ)
         IF(ICHX.NE.2) THEN
            CALL LCMPUT(IPTRK,'IPZ',LL4,1,IPZ)
            DEALLOCATE(IPZ)
         ENDIF
      ELSE
         IF(ICHX.NE.2) DEALLOCATE(IPZ)
      ENDIF
      DEALLOCATE(MUZ,MUY,MUX)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(MAT,IDL,VOL)
      RETURN
*
  600 FORMAT(/26H TRITRK: VALUES OF VECTOR ,A6,4H ARE/(1X,1P,20I6))
  620 FORMAT (///11H REGION    ,8(I8,6X,1HI))
  630 FORMAT (   11H MIXTURE   ,8(I8,6X,1HI))
  640 FORMAT (   11H POINTER   ,8(I8,6X,1HI))
  650 FORMAT (   11H VOLUME    ,8(1P,E13.6,2H I))
      END
