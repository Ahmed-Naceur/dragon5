*DECK NSST
      SUBROUTINE NSST(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Nodal expansion method (NEM) tracking operator.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): create or modification type(L_TRACK);
*         HENTRY(2): read-only type(L_GEOM).
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      CHARACTER HENTRY(NENTRY)*12
      TYPE(C_PTR) KENTRY(NENTRY)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      CHARACTER TEXT4*4,TEXT12*12,TITLE*72,HSIGN*12
      DOUBLE PRECISION DFLOTT
      LOGICAL CYLIND,ILK,LUMP
      INTEGER ISTATE(NSTATE),NCODE(6),ICODE(6)
      REAL ZCODE(6)
      TYPE(C_PTR) IPGEO,IPTRK
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAT,IDL,ISPLX,ISPLY,ISPLZ
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: KN,IQFR
      REAL, ALLOCATABLE, DIMENSION(:) :: XXX,YYY,ZZZ,XX,YY,ZZ,DD,VOL
      REAL, ALLOCATABLE, DIMENSION(:,:) :: QFR
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.NE.2) CALL XABORT('NSST: TWO PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('NSST: L'
     1 //'CM OBJECT EXPECTED AT LHS.')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) CALL XABORT('NSST: E'
     1 //'NTRY IN CREATE OR MODIFICATION MODE EXPECTED.')
      IF((JENTRY(2).NE.2).OR.((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2)))
     1 CALL XABORT('NSST: LCM OBJECT IN READ-ONLY MODE EXPECTED AT R'
     2 //'HS.')
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_GEOM') THEN
         TEXT12=HENTRY(2)
         CALL XABORT('NSST: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_GEOM EXPECTED.')
      ENDIF
      IPTRK=KENTRY(1)
      IPGEO=KENTRY(2)
      HSIGN='L_TRACK'
      CALL LCMPTC(IPTRK,'SIGNATURE',12,1,HSIGN)
      HSIGN='TRIVAC'
      CALL LCMPTC(IPTRK,'TRACK-TYPE',12,1,HSIGN)
      CALL LCMGET(IPGEO,'STATE-VECTOR',ISTATE)
      IDIM=0
      ITYPE=ISTATE(1)
      IF(ITYPE.EQ.2) THEN
        IDIM=1
      ELSE IF(ITYPE.EQ.5) THEN
        IDIM=2
      ELSE IF(ITYPE.EQ.7) THEN
        IDIM=3
      ELSE
        CALL XABORT('NSST: 1D, 2D OR 3D CARTESIAN GEOMETRY EXPECTED.')
      ENDIF
      LX=ISTATE(3)
      LY=ISTATE(4)
      LZ=ISTATE(5)
      CALL LCMLEN(IPGEO,'BIHET',ILONG,ITYLCM)
      IF(ILONG.NE.0) CALL XABORT('NSST: DOUBLE-HETEROGENEITY NOT SUPPO'
     1 //'RTED.')
*
      IMPX=1
      TITLE=' '
      IGMAX=0
      LUMP=.FALSE.
      IF(JENTRY(1).EQ.1) THEN
         CALL LCMGTC(IPTRK,'SIGNATURE',12,1,HSIGN)
         IF(HSIGN.NE.'L_TRACK') THEN
            TEXT12=HENTRY(1)
            CALL XABORT('NSST: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_TRACK EXPECTED.')
         ENDIF
         CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,HSIGN)
         IF(HSIGN.NE.'TRIVAC') THEN
            TEXT12=HENTRY(3)
            CALL XABORT('NSST: TRACK-TYPE OF '//TEXT12//' IS '//HSIGN
     1      //'. TRIVAC EXPECTED.')
         ENDIF
         CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
         IGMAX=ISTATE(39)
         CALL LCMLEN(IPTRK,'TITLE',LENGT,ITYLCM)
         IF(LENGT.GT.0) CALL LCMGTC(IPTRK,'TITLE',72,1,TITLE)
      ENDIF
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.EQ.10) GO TO 30
      IF(INDIC.NE.3) CALL XABORT('NSST: CHARACTER DATA EXPECTED.')
      IF(TEXT4.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('NSST: INTEGER DATA EXPECTED(1).')
      ELSE IF(TEXT4.EQ.'TITL') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TITLE,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('NSST: TITLE EXPECTED.')
      ELSE IF(TEXT4.EQ.'LUMP') THEN
         LUMP=.TRUE.
      ELSE IF(TEXT4.EQ.'HYPE') THEN
        CALL REDGET(INDIC,IGMAX,FLOTT,TEXT4,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('NSST: INTEGER DATA EXPECTED(1).')
      ELSE IF(TEXT4.EQ.';') THEN
         GO TO 30
      ELSE
         CALL XABORT('NSST: '//TEXT4//' IS AN INVALID KEYWORD.')
      ENDIF
      GO TO 10
*----
*  SCRATCH STORAGE ALLOCATION
*----
   30 IF(IDIM.EQ.1) THEN
        MAXPTS=LX
      ELSE IF(IDIM.EQ.2) THEN
        MAXPTS=LX*LY
      ELSE
        MAXPTS=LX*LY*LZ
      ENDIF
      ALLOCATE(XXX(LX+1),YYY(LY+1),ZZZ(LZ+1),MAT(MAXPTS),IDL(MAXPTS),
     1 VOL(MAXPTS),XX(MAXPTS),YY(MAXPTS),ZZ(MAXPTS),DD(MAXPTS),
     2 KN(6,MAXPTS),QFR(6,MAXPTS),IQFR(6,MAXPTS))
*----
*  RECOVER TRACKING INFORMATION
*----
      ALLOCATE(ISPLX(LX),ISPLY(LY),ISPLZ(LZ))
      CALL READ3D(LX,LY,LZ,MAXPTS,IPGEO,IHEX,IR,ILK,SIDE,XXX,YYY,ZZZ,
     1 IMPX,LX,LY,LZ,MAT,NEL,NCODE,ICODE,ZCODE,ISPLX,ISPLY,ISPLZ,ISPLH,
     2 ISPLL)
      DEALLOCATE(ISPLX,ISPLY,ISPLZ)
      IF(IDIM.EQ.1) THEN
*       1D GEOMETRY
        LY=1
        NCODE(3)=2
        NCODE(4)=5
        ZCODE(3)=1.0
        ZCODE(4)=1.0
        YYY(1)=0.0
        YYY(2)=2.0
      ENDIF
      IF(IDIM.LE.2) THEN
*       1D OR 2D GEOMETRY
        LZ=1
        NCODE(5)=2
        NCODE(6)=5
        ZCODE(5)=1.0
        ZCODE(6)=1.0
        ZZZ(1)=0.0
        ZZZ(2)=2.0
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
         DO IZ=LZ,1,-1
           IOFF=(IZ-1)*LX*LY
           DO IY=LY,1,-1
             DO IX=LX,IY+1,-1
               MAT(IOFF+(IY-1)*LX+IX)=MAT(IOFF+(IX-1)*LY+IY)
             ENDDO
             DO IX=IY,1,-1
               MAT(IOFF+(IY-1)*LX+IX)=MAT(K)
               K=K-1
             ENDDO
           ENDDO
         ENDDO
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
         DO IZ=LZ,1,-1
           IOFF=(IZ-1)*LX*LY
           DO IY=LY,1,-1
             DO IX=LX,IY,-1
               MAT(IOFF+(IY-1)*LX+IX)=MAT(K)
               K=K-1
             ENDDO
           ENDDO
         ENDDO
         DO IZ=1,LZ
           IOFF=(IZ-1)*LX*LY
           DO IY=1,LY
             DO IX=1,IY-1
               MAT(IOFF+(IY-1)*LX+IX)=MAT(IOFF+(IX-1)*LY+IY)
             ENDDO
           ENDDO
         ENDDO
         NEL=LX*LY*LZ
         IF(K.NE.0) THEN
            CALL XABORT('TRITRK: UNABLE TO UNFOLD THE DOMAIN(2).')
         ENDIF
      ENDIF
      IF(IMPX.GT.5) THEN
         WRITE(6,120) 'NCODE',(NCODE(I),I=1,6)
         WRITE(6,120) 'MAT',(MAT(I),I=1,LX*LY*LZ)
      ENDIF
*----
*  SET TRACKING INFORMATION
*----
      CYLIND=.FALSE.
      CALL TRIDFC(IMPX,LX,LY,LZ,CYLIND,NCODE,ICODE,ZCODE,MAT,XXX,
     1 YYY,ZZZ,LL0,VOL,XX,YY,ZZ,DD,KN,QFR,IQFR)
      INUM=0
      DO KEL=1,NEL
        IF(MAT(KEL).EQ.0) CYCLE
        INUM=INUM+1
        IF(LUMP) THEN
          IDL(KEL)=INUM
        ELSE
          IDL(KEL)=(INUM-1)*IDIM*5+1
        ENDIF
      ENDDO
      IF(LUMP) THEN
        NUN=INUM
      ELSE
        NUN=INUM*IDIM*5
      ENDIF
      ISTATE(:)=0
      ISTATE(1)=NEL
      ISTATE(2)=NUN
      ISTATE(4)=MAXVAL(MAT(:NEL))
      ISTATE(6)=ITYPE ! Geometry type
      ISTATE(12)=4 ! Nodal expansion method
      ISTATE(14)=LX
      IF(IDIM.GE.2) ISTATE(15)=LY
      IF(IDIM.EQ.3) ISTATE(16)=LZ
      IF(LUMP) ISTATE(38)=1
      ISTATE(39)=IGMAX
      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPTRK,'NCODE',6,1,NCODE)
      CALL LCMPUT(IPTRK,'ZCODE',6,2,ZCODE)
      CALL LCMPUT(IPTRK,'ICODE',6,1,ICODE)
      CALL LCMPUT(IPTRK,'MATCOD',NEL,1,MAT)
      CALL LCMPUT(IPTRK,'VOLUME',NEL,2,VOL)
      CALL LCMPUT(IPTRK,'KEYFLX',NEL,1,IDL)
      CALL LCMPUT(IPTRK,'XX',LX,2,XX)
      IF(IDIM.GE.2) CALL LCMPUT(IPTRK,'YY',LY,2,YY)
      IF(IDIM.EQ.3) CALL LCMPUT(IPTRK,'ZZ',LY,2,ZZ)
      CALL LCMPUT(IPTRK,'KN',6*NEL,1,KN)
      CALL LCMPUT(IPTRK,'QFR',6*NEL,2,QFR)
      CALL LCMPUT(IPTRK,'IQFR',6*NEL,1,IQFR)
      IF(TITLE.NE.' ') CALL LCMPTC(IPTRK,'TITLE',72,1,TITLE)
      TEXT12=HENTRY(2)
      CALL LCMPTC(IPTRK,'LINK.GEOM',12,1,TEXT12)
      IF(IMPX.GT.1) THEN
         WRITE(6,100) TITLE
         CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
         WRITE(6,110) ISTATE(1:2),ISTATE(4),ISTATE(6),ISTATE(14:16),
     1   ISTATE(38),ISTATE(39)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DD,XX,ZZ,YY,ZZZ,YYY,XXX)
      DEALLOCATE(IQFR,QFR,KN,VOL,IDL,MAT)
      RETURN
*
  100 FORMAT(1H1,24HNN     NN  SSSSS   SSSSS,
     1  97(1H*)/26H NNN    NN SSSSSSS SSSSSSS,
     2 58(1H*),38H MULTIGROUP VERSION.  A. HEBERT (2021)/
     3 26H NNNN   NN SS   SS SS   SS/26H NN NN  NN  SSS     SSS   /
     4 26H NN  NN NN    SSS     SSS /26H NN   NNNN SS   SS SS   SS/
     5 26H NN    NNN SSSSSSS SSSSSSS/26H NN     NN  SSSSS   SSSSS //
     6 1X,A72//)
  110 FORMAT(/14H STATE VECTOR:/
     1 7H NREG  ,I8,22H   (NUMBER OF REGIONS)/
     2 7H NUN   ,I8,23H   (NUMBER OF UNKNOWNS)/
     3 7H NMIX  ,I8,23H   (NUMBER OF MIXTURES)/
     4 7H ITYPE ,I8,41H   (TYPE OF GEOMETRY -- 2:1D; 5:2D; 7:3D)/
     5 7H LX    ,I8,40H   (NUMBER OF ELEMENTS ALONG THE X AXIS)/
     6 7H LY    ,I8,40H   (NUMBER OF ELEMENTS ALONG THE Y AXIS)/
     7 7H LZ    ,I8,40H   (NUMBER OF ELEMENTS ALONG THE Z AXIS)/
     8 7H LUMP  ,I8,35H   (UNKNOWN LUMP OPTION 0/1:OFF/ON)/
     9 7H IGMAX ,I8,47H   (ENERGY GROUP LIMIT WITH HYPERBOLIC TRIAL FU,
     1 8HNCTIONS))
  120 FORMAT(/24H NSST: VALUES OF VECTOR ,A6,4H ARE/(1X,1P,20I6))
      END
