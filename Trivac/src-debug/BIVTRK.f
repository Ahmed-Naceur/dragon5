*DECK BIVTRK
      SUBROUTINE BIVTRK (MAXPTS,IPTRK,IPGEOM,IMPX,IELEM,ICOL,NLF,NVD,
     1 ISPN,ISCAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover of the geometry and tracking for BIVAC.
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
* IPTRK   L_TRACK pointer to the tracking information.
* IPGEOM  L_GEOM pointer to the geometry.
* IMPX    print flag.
* IELEM   degree of the Lagrangian finite elements:
*         <0: order -IELEM primal finite elements;
*         >0: order IELEM dual finite elements.
* ICOL    type of quadrature used to integrate the mass matrix:
*         =1: analytical integration;
*         =2: Gauss-Lobatto quadrature (collocation method);
*         =3: Gauss Legendre quadrature (superconvergent).
*         =4: mesh centered finite differences in hexagonal geometry.
*         IELEM=-1 and ICOL=2 : mesh corner finite differences;
*         IELEM=1 and ICOL=2 : mesh centered finite differences.
* NLF     number of Legendre orders for the flux. Equal to zero for
*         diffusion theory.
* NVD     type of void boundary condition if NLF>0 and ICOL=3.
* ISPN    type of transport solution:
*         =0: complete PN method;
*         =1: simplified PN method.
* ISCAT   source anisotropy:
*         =1: isotropic sources in laboratory system;
*         =2: linearly anisotropic sources in laboratory system.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPGEOM
      INTEGER MAXPTS,IMPX,IELEM,ICOL,NLF,NVD,ISPN,ISCAT
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      LOGICAL ILK,CYLIND
      CHARACTER HSMG*131
      INTEGER ISTATE(NSTATE),IGP(NSTATE),NCODE(6),ICODE(6)
      REAL ZCODE(6)
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT,IDL,IPERT,KN,IQFR,MU
      REAL, DIMENSION(:), ALLOCATABLE :: VOL,XXX,YYY,ZZZ,XX,YY,DD,QFR,
     1 BFR,ISPLX,ISPLY,ISPLZ
*
******************* BIVAC GEOMETRICAL STRUCTURE. ***********************
*                                                                      *
*   ITYPE       : =2  : CARTESIAN 1-D GEOMETRY;                        *
*                 =3  : TUBE 1-D GEOMETRY;                             *
*                 =4  : SPHERICAL 1-D GEOMETRY;                        *
*                 =5  : CARTESIAN 2-D GEOMETRY;                        *
*                 =6  : TUBE 2-D GEOMETRY;                             *
*                 =8  : HEXAGONAL 2-D GEOMETRY.                        *
*   IHEX        : TYPE OF HEXAGONAL SYMMETRY.                          *
*   IELEM       : .LT.0 : ORDER -IELEM PRIMAL FINITE ELEMENTS;         *
*                 .GT.0 : ORDER IELEM DUAL FINITE ELEMENTS.            *
*   ICOL        : TYPE OF QUADRATURE USED TO INTEGRATE THE MASS MATRIX.*
*                 =1  : ANALYTICAL INTEGRATION;                        *
*                 =2  : GAUSS-LOBATTO QUADRATURE (COLLOCATION METHOD); *
*                 =3  : GAUSS LEGENDRE QUADRATURE (SUPERCONVERGENT).   *
*                IELEM=-1 AND ICOL=2 : MESH CORNER FINITE DIFFERENCES. *
*                IELEM=1 AND ICOL=2 : MESH CENTERED FINITE DIFFERENCES.*
*   ISPLH       : TYPE OF HEXAGONAL MESH-SPLITTING.                    *
*                 =1 : NO MESH SPLITTING (COMPLETE HEXAGONS);          *
*                 =K : 6*(K-1)*(K-1) TRIANGLES PER HEXAGON.            *
*   SIDE        : SIDE OF THE HEXAGONS.                                *
*   LL4         : ORDER OF THE MATRICES PER GROUP IN BIVAC.            *
*   NCODE       : TYPES OF BOUNDARY CONDITIONS. DIMENSION=6            *
*   ZCODE       : ALBEDOS. DIMENSION=6                                 *
*   LX          : NUMBER OF ELEMENTS ALONG THE X AXIS.                 *
*   LY          : NUMBER OF ELEMENTS ALONG THE Y AXIS.                 *
*   XX          : X-DIRECTED MESH SPACINGS. DIMENSION=LX*LY            *
*   YY          : Y-DIRECTED MESH SPACINGS. DIMENSION=LX*LY            *
*   DD          : USED WITH CYLINDRICAL GEOMETRIES. DIMENSION=LX*LY    *
*   KN          : ELEMENT-ORDERED UNKNOWN LIST. DIMENSION LX*LY*ICOEF  *
*                 WHERE ICOEF IS THE NUMBER OF UNKNOWN PER ELEMENT.    *
*   QFR         : ELEMENT-ORDERED BOUNDARY CONDITIONS.                 *
*                 DIMENSION 4*LX*LY                                    *
*   IQFR        : ELEMENT-ORDERED PHYSICAL ALBEDO INDICES.             *
*                 DIMENSION 4*LX*LY                                    *
*   BFR         : ELEMENT-ORDERED SURFACE FRACTIONS.                   *
*                 DIMENSION 4*LX*LY                                    *
*   MU          : INDICES USED WITH COMPRESSED DIAGONAL STORAGE MODE   *
*                 MATRICES. DIMENSION MAXEV                            *
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
      IF(ISTATE(9).EQ.0) THEN
         IF((ITYPE.NE.1).AND.(ITYPE.NE.2).AND.(ITYPE.NE.3).AND.
     1      (ITYPE.NE.4).AND.(ITYPE.NE.5).AND.(ITYPE.NE.6).AND.
     2      (ITYPE.NE.8)) THEN
            CALL XABORT('BIVTRK: DISCRETIZATION NOT AVAILABLE.')
         ENDIF
         ALLOCATE(XXX(MAXPTS+1),YYY(MAXPTS+1),ZZZ(MAXPTS+1))
*
         ALLOCATE(ISPLX(MAXPTS),ISPLY(MAXPTS),ISPLZ(MAXPTS))
         CALL READ3D(MAXPTS,MAXPTS,MAXPTS,MAXPTS,IPGEOM,IHEX,IR,ILK,
     1   SIDE,XXX,YYY,ZZZ,IMPX,LX,LY,LZ,MAT,NEL,NCODE,ICODE,ZCODE,
     2   ISPLX,ISPLY,ISPLZ,ISPLH,ISPLL)
         DEALLOCATE(ISPLX,ISPLY,ISPLZ)
         IF((ITYPE.EQ.8).AND.(IELEM.GT.0).AND.(ICOL.LE.3)) THEN
           IF(ISPLL.EQ.0) THEN
             CALL XABORT('BIVTRK: SPLITL KEYWORD MISSING IN GEOMETRY.')
           ENDIF
           ISPLH=ISPLL
         ELSE IF(ITYPE.EQ.8) THEN
           ISPLH=ISPLH+1  
         ENDIF
      ELSE
         CALL XABORT('BIVTRK: DISCRETIZATION NOT AVAILABLE.')
      ENDIF
      IF((IMPX.GE.1).AND.(ITYPE.NE.8)) THEN
         WRITE (6,'(/39H BIVTRK: TYPE OF FINITE ELEMENT IELEM =,I3,
     1   8H  ICOL =,I3/)') IELEM,ICOL
      ELSE IF(IMPX.GE.1) THEN
         WRITE (6,'(/39H BIVTRK: TYPE OF FINITE ELEMENT IELEM =,I3,
     1   8H  ICOL =,I3,9H  ISPLH =,I3/)') IELEM,ICOL,ISPLH
      ENDIF
*
      IF(LX*LY*LZ.GT.MAXPTS) THEN
         WRITE (HSMG,'(39HBIVTRK: MAXPTS SHOULD BE INCREASED FROM,I7,
     1   3H TO,I7)') MAXPTS,LX*LY*LZ
         CALL XABORT(HSMG)
      ENDIF
*----
*  1-D AND 2-D CYLINDRICAL CASES.
*----
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.4).OR.(ITYPE.EQ.6)
      IF((ITYPE.EQ.2).OR.(ITYPE.EQ.3)) THEN
         NCODE(3)=2
         NCODE(4)=5
         ICODE(3)=0
         ICODE(4)=0
         ZCODE(3)=1.0
         ZCODE(4)=1.0
         YYY(1)=0.0
         YYY(2)=2.0
      ELSE IF(ITYPE.EQ.6) THEN
         LY=LZ
         DO 10 I=1,LZ+1
         YYY(I)=ZZZ(I)
   10    CONTINUE
         NCODE(3)=NCODE(5)
         NCODE(4)=NCODE(6)
         ICODE(3)=ICODE(5)
         ICODE(4)=ICODE(6)
         ZCODE(3)=ZCODE(5)
         ZCODE(4)=ZCODE(6)
      ENDIF
*----
*  UNFOLD THE DOMAIN IN DIAGONAL SYMMETRY CASES.
*----
      IF((NCODE(2).EQ.3).AND.(NCODE(3).EQ.3)) THEN
         NCODE(3)=NCODE(1)
         NCODE(2)=NCODE(4)
         ICODE(3)=ICODE(1)
         ICODE(2)=ICODE(4)
         ZCODE(3)=ZCODE(1)
         ZCODE(2)=ZCODE(4)
         K=LX*(LX+1)/2
         DO 35 IY=LY,1,-1
         DO 20 IX=LX,IY+1,-1
         MAT((IY-1)*LX+IX)=MAT((IX-1)*LY+IY)
   20    CONTINUE
         DO 30 IX=IY,1,-1
         MAT((IY-1)*LX+IX)=MAT(K)
         K=K-1
   30    CONTINUE
   35    CONTINUE
         NEL=LX*LY
         IF(K.NE.0) THEN
            CALL XABORT('BIVTRK: UNABLE TO UNFOLD THE DOMAIN.')
         ENDIF
      ELSE IF((NCODE(1).EQ.3).AND.(NCODE(4).EQ.3)) THEN
         NCODE(1)=NCODE(3)
         NCODE(4)=NCODE(2)
         ICODE(1)=ICODE(3)
         ICODE(4)=ICODE(2)
         ZCODE(1)=ZCODE(3)
         ZCODE(4)=ZCODE(2)
         K=LX*(LX+1)/2
         DO 45 IY=LY,1,-1
         DO 40 IX=LX,IY,-1
         MAT((IY-1)*LX+IX)=MAT(K)
         K=K-1
   40    CONTINUE
   45    CONTINUE
         DO 55 IY=1,LY
         DO 50 IX=1,IY-1
         MAT((IY-1)*LX+IX)=MAT((IX-1)*LY+IY)
   50    CONTINUE
   55    CONTINUE
         NEL=LX*LY
         IF(K.NE.0) THEN
            CALL XABORT('BIVTRK: UNABLE TO UNFOLD THE DOMAIN.')
         ENDIF
      ENDIF
      IF(IMPX.GT.5) THEN
         WRITE(6,600) 'NCODE',(NCODE(I),I=1,4)
         WRITE(6,600) 'MAT',(MAT(I),I=1,LX*LY)
      ENDIF
*
      IF((IELEM.LT.0).AND.(ITYPE.NE.8)) THEN
         IEL=-IELEM
         MAXEV=(IEL*LX+1)*(IEL*LY+1)
         MAXKN=(IEL+1)*(IEL+1)*NEL
         MAXQF=4*NEL
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.3).AND.(NLF.NE.0)) THEN
*        PN METHOD / 1D CYLINDRICAL GEOMETRY.
         MAXEV=(2*LX+1)*(NLF/2)*(NLF/2+1)/2
         MAXKN=3*NEL*(NLF/2)*(NLF/2)
         MAXQF=2*NEL
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.4).AND.(NLF.NE.0)) THEN
*        PN METHOD / 1D SPHERICAL GEOMETRY.
         MAXEV=(2*LX+1)*(NLF/2)
         MAXKN=3*NEL*(NLF/2)
         MAXQF=2*NEL
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.5).AND.(NLF.NE.0).AND.
     1 (ISPN.EQ.0)) THEN
*        PN METHOD / 2D CARTESIAN GEOMETRY.
         MAXEV=0
         DO 60 IL=1,NLF-1,2
         MAXEV=MAXEV+(IL*LX+(IL+1)*(LX+1))*LY+(IL+1)*(LX+1)
   60    CONTINUE
         MAXKN=5*NEL*NLF*(NLF/2)
         MAXQF=4*NEL
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.5).AND.(NLF.NE.0).AND.
     1 (ISPN.EQ.1)) THEN
*        SPN METHOD / 2D CARTESIAN GEOMETRY.
         MAXEV=(LX+1)*LY*IELEM+LX*(LY+1)*IELEM+LX*LY*IELEM*IELEM
         MAXEV=MAXEV*NLF/2
         MAXKN=5*NEL
         MAXQF=4*NEL
      ELSE IF((IELEM.GT.0).AND.(ITYPE.NE.8)) THEN
         MAXEV=(LX+1)*LY*IELEM+LX*(LY+1)*IELEM+LX*LY*IELEM*IELEM
         MAXKN=5*NEL
         MAXQF=4*NEL
      ELSE IF((IELEM.LT.0).AND.(ITYPE.EQ.8)) THEN
         IEL=-IELEM
         NEL=LX
         IF(ISPLH.EQ.1) THEN
            MAXEV=6*NEL
            MAXKN=7*NEL
         ELSE
            MAXEV=(1+ISPLH*(ISPLH-1)*3)*NEL
            MAXKN=(6*(ISPLH-1)**2)*NEL*4
         ENDIF
         MAXQF=MAXKN
      ELSE IF((ICOL.EQ.4).AND.(ITYPE.EQ.8)) THEN
         NEL=LX
         IF(ISPLH.EQ.1) THEN
            MAXEV=NEL
            MAXKN=7*NEL
         ELSE
            MAXEV=(6*(ISPLH-1)**2)*NEL
            MAXKN=(6*(ISPLH-1)**2)*NEL*4
         ENDIF
         MAXQF=MAXKN
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.8)) THEN
         NEL=LX
         LXH=LX/(3*ISPLH**2)
         NBLOS=LXH*ISPLH**2
         NBC=INT((SQRT(REAL((4*LXH-1)/3))+1.)/2.)
         MAXEV=3*(2*NBLOS*IELEM+(2*NBC-1)*ISPLH)*IELEM+3*NBLOS*IELEM**2
         MAXKN=(LXH*ISPLH**2)*(4+6*IELEM*(IELEM+1))
         MAXQF=(LXH*ISPLH**2)*6
      ELSE
         CALL XABORT('BIVTRK: INVALID TYPE OF DISCRETIZATION.')
      ENDIF
      IF(CYLIND) THEN
         MAXDD=NEL
      ELSE
         MAXDD=1
      ENDIF
      IF((ICOL.EQ.4).AND.(ITYPE.EQ.8).AND.(IELEM.NE.1)) THEN
         CALL XABORT('BIVTRK: THIS HEXAGONAL MCFD DISCRETIZATIONS IS L'
     1   //'IMITED TO LINEAR ORDER.')
      ELSE IF((IELEM.LT.0).AND.(ITYPE.EQ.8).AND.(IELEM.NE.-1)) THEN
         CALL XABORT('BIVTRK: THIS HEXAGONAL PRIM DISCRETIZATIONS IS L'
     1   //'IMITED TO LINEAR ORDER.')
      ENDIF
      IF(ICOL.LE.3) CALL BIVCOL(IPTRK,IMPX,ABS(IELEM),ICOL)
      ALLOCATE(XX(NEL),YY(NEL),DD(MAXDD),KN(MAXKN),QFR(MAXQF),
     1 IQFR(MAXQF),BFR(MAXQF),MU(MAXEV))
      CALL XDISET(KN,MAXKN,0)
      CALL XDRSET(QFR,MAXQF,0.0)
      CALL XDISET(IQFR,MAXQF,0)
      CALL XDRSET(BFR,MAXQF,0.0)
      IF((IELEM.LT.0).AND.(ITYPE.NE.8)) THEN
         IEL=-IELEM
         CALL BIVPKN(MAXEV,IMPX,LX,LY,CYLIND,IEL,LL4,NCODE,ICODE,ZCODE,
     1   MAT,VOL,XXX,YYY,XX,YY,DD,KN,QFR,IQFR,BFR,MU)
      ELSE IF(((ITYPE.EQ.2).OR.((ITYPE.EQ.5).AND.(ISPN.EQ.1))).AND.
     1 (IELEM.GT.0).AND.(NLF.NE.0)) THEN
*        MIXED-DUAL SPN APPROXIMATION IN 1D OR 2D CARTESIAN GEOMETRY.
         CALL BIVDKN(MAXEV,IMPX,LX,LY,CYLIND,IELEM,ICOL,LL4,NCODE,
     1   ICODE,ZCODE,MAT,VOL,XXX,YYY,XX,YY,DD,KN,QFR,IQFR,BFR,IDL,MU)
         NUN=LL4*NLF/2
      ELSE IF((IELEM.GT.0).AND.(ITYPE.NE.8)) THEN
         CALL BIVDKN(MAXEV,IMPX,LX,LY,CYLIND,IELEM,ICOL,LL4,NCODE,
     1   ICODE,ZCODE,MAT,VOL,XXX,YYY,XX,YY,DD,KN,QFR,IQFR,BFR,IDL,MU)
         NUN=LL4
      ELSE IF((IELEM.LT.0).AND.(ITYPE.EQ.8)) THEN
*        HEXAGONAL GEOMETRY MESH CORNER FINITE DIFFERENCES.
         CALL BIVPRH(MAXEV,MAXKN,IMPX,ISPLH,LX,IHEX,NCODE,ICODE,ZCODE,
     1   MAT,SIDE,LL4,NELEM,VOL,KN,QFR,IQFR,BFR,MU)
         IF(ISPLH.EQ.1) THEN
            MAXKN=7*NELEM
            MAXQF=7*NELEM
         ELSE
            MAXKN=4*NELEM
            MAXQF=4*NELEM
         ENDIF
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.8).AND.(ICOL.EQ.4)) THEN
*        HEXAGONAL GEOMETRY MESH CENTERED FINITE DIFFERENCES.
         CALL BIVDFH(MAXEV,MAXKN,IMPX,ISPLH,LX,SIDE,LL4,NUN,IHEX,
     1   NCODE,ICODE,ZCODE,MAT,VOL,IDL,KN,QFR,IQFR,BFR,MU)
         IF(ISPLH.EQ.1) THEN
            MAXKN=7*LL4
            MAXQF=7*LL4
         ELSE
            MAXKN=4*LL4
            MAXQF=4*LL4
         ENDIF
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.8)) THEN
*        HEXAGONAL GEOMETRY THOMAS-RAVIART-SCHNEIDER FINITE ELEMENTS.
         NBLOS=LXH*ISPLH**2
         ALLOCATE(IPERT(NBLOS))
         CALL BIVSFH(MAXEV,NBLOS,IMPX,ISPLH,IELEM,LXH,MAT,SIDE,NCODE,
     1   ICODE,ZCODE,LL4,VOL,IDL,IPERT,KN,QFR,IQFR,BFR,MU)
         CALL LCMPUT(IPTRK,'IPERT',NBLOS,1,IPERT)
         DEALLOCATE(IPERT)
         NUN=LL4
      ENDIF
      DEALLOCATE(YYY,ZZZ)
*----
*  APPEND THE PN FLUXES AT THE END OF UNKNOWN VECTOR.
*----
      IF(NLF.GE.2) THEN
         IF((ITYPE.EQ.2).OR.((ITYPE.EQ.5).AND.(ISPN.EQ.1))) THEN
            NUN=LL4+LL4*(NLF-2)/2
         ELSE IF((ITYPE.EQ.8).AND.(ISPN.EQ.1)) THEN
            NUN=NUN+NUN*(NLF-2)/2
         ELSE IF((ITYPE.NE.2).AND.(ITYPE.NE.5).AND.(ITYPE.NE.8)) THEN
            CALL XABORT('BIVTRK: GEOMETRY NOT SUPPORTED WITH PN.')
         ENDIF
      ENDIF
*----
*  APPEND THE AVERAGED FLUXES AT THE END OF UNKNOWN VECTOR.
*----
      IF(IELEM.LT.0) THEN
         NUN=LL4
         DO 190 I=1,NEL
         IF(MAT(I).EQ.0) THEN
            IDL(I)=0
         ELSE
            NUN=NUN+1
            IDL(I)=NUN
         ENDIF
  190    CONTINUE
      ENDIF
*----
*  RESERVE A COMPONENT TO STORE THE SURFACE-AVERAGED FLUX.
*----
      IF(NLF.EQ.0) NUN=NUN+1
      IF(IMPX.GT.0) WRITE (6,'(/34H BIVTRK: ORDER OF LINEAR SYSTEMS =,
     1 I7/9X,37HNUMBER OF UNKNOWNS PER ENERGY GROUP =,I7)') LL4,NUN
*
      IF(IMPX.GT.5) THEN
         I1=1
         DO 200 I=1,(NEL-1)/8+1
         I2=I1+7
         IF(I2.GT.NEL) I2=NEL
         WRITE (6,620) (J,J=I1,I2)
         WRITE (6,630) (MAT(J),J=I1,I2)
         WRITE (6,640) (IDL(J),J=I1,I2)
         WRITE (6,650) (VOL(J),J=I1,I2)
         I1=I1+8
  200    CONTINUE
      ENDIF
*----
*  SAVE GENERAL AND BIVAC-SPECIFIC TRACKING INFORMATION.
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
      IGP(5)=1
      IGP(6)=ITYPE
      IGP(7)=IHEX
      IGP(8)=IELEM
      IGP(9)=ICOL
      IGP(10)=ISPLH
      IGP(11)=LL4
      IGP(12)=LX
      IGP(13)=LY
      IGP(14)=NLF
      IGP(15)=ISPN
      IGP(16)=ISCAT
      IGP(17)=NVD
      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,IGP)
      CALL LCMPUT(IPTRK,'MATCOD',NEL,1,MAT)
      CALL LCMPUT(IPTRK,'VOLUME',NEL,2,VOL)
      CALL LCMPUT(IPTRK,'KEYFLX',NEL,1,IDL)
      CALL LCMPUT(IPTRK,'NCODE',6,1,NCODE)
      CALL LCMPUT(IPTRK,'ZCODE',6,2,ZCODE)
      CALL LCMPUT(IPTRK,'ICODE',6,1,ICODE)
      CALL LCMPUT(IPTRK,'BC-REFL+TRAN',1,1,NUN)
      IF(ITYPE.EQ.4) CALL LCMPUT(IPTRK,'XXX',LX+1,2,XXX)
      DEALLOCATE(XXX)
      IF(ITYPE.EQ.8) THEN
         CALL LCMPUT(IPTRK,'SIDE',1,2,SIDE)
      ELSE
         CALL LCMPUT(IPTRK,'XX',LX*LY,2,XX)
         CALL LCMPUT(IPTRK,'YY',LX*LY,2,YY)
         IF(.NOT.CYLIND) DD(1)=0.0
         CALL LCMPUT(IPTRK,'DD',MAXDD,2,DD)
      ENDIF
      DEALLOCATE(XX,YY,DD)
      CALL LCMPUT(IPTRK,'KN',MAXKN,1,KN)
      DEALLOCATE(KN)
      CALL LCMPUT(IPTRK,'QFR',MAXQF,2,QFR)
      DEALLOCATE(QFR)
      CALL LCMPUT(IPTRK,'IQFR',MAXQF,1,IQFR)
      DEALLOCATE(IQFR)
      CALL LCMPUT(IPTRK,'BFR',MAXQF,2,BFR)
      DEALLOCATE(BFR)
      CALL LCMPUT(IPTRK,'MU',LL4,1,MU)
      DEALLOCATE(MU)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(MAT,IDL,VOL)
      RETURN
*
  600 FORMAT(/26H BIVTRK: VALUES OF VECTOR ,A6,4H ARE/(1X,1P,20I6))
  620 FORMAT (///11H REGION    ,8(I8,6X,1HI))
  630 FORMAT (   11H MIXTURE   ,8(I8,6X,1HI))
  640 FORMAT (   11H POINTER   ,8(I8,6X,1HI))
  650 FORMAT (   11H VOLUME    ,8(1P,E13.6,2H I))
      END
