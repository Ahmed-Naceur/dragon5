*DECK READ3D
      SUBROUTINE READ3D (MAXX,MAXY,MAXZ,MAXPTS,IPGEOM,IHEX,IR,ILK,SIDE,
     1 XXX,YYY,ZZZ,IMPX,LX,LY,LZ,MAT,NMBLK,NCODE,ICODE,ZCODE,ISPLTX,
     2 ISPLTY,ISPLTZ,ISPLTH,ISPLTL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the input data for the description of a 1-D, 2-D or 3-D
* Cartesian, cylindrical, spherical or hexagonal domain.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):A. Hebert
*
*Parameters: input/output
* MAXX    allocated storage for arrays of dimension LX.
* MAXY    allocated storage for arrays of dimension LY.
* MAXZ    allocated storage for arrays of dimension LZ.
* MAXPTS  allocated storage for arrays of dimension NMBLK.
* IPGEOM  L_GEOM pointer to the geometry.
* IHEX    type of hexagonal geometry (=0 for non-hexagonal geometry).
* IR      number of mixtures.
* ILK     (ILK=.true. if neutron leakage through external boundary
*         is present).
* SIDE    side of the hexagons. XXX and YYY arrays are not used with
*         hexagonal geometry.
* XXX     Cartesian coordinates of the domain along the X-axis.
* YYY     Cartesian coordinates of the domain along the Y-axis.
* ZZZ     Cartesian coordinates of the domain along the Z-axis.
* IMPX    print flag. Minimum printing if IMPX=0.
* LX      number of elements along the X-axis after mesh-splitting
*         or number of hexagons in one axial plane.
* LY      number of elements along the Y-axis.
* LZ      number of elements along the Z-axis.
* MAT     index-number of the mixture type assigned to each volume
*         after mesh-splitting.
* NMBLK   number of elements in the domain.
* NCODE   boundary condition relative to each side of the domain:
*         =1: VOID ; =2: REFL ; =3: DIAG ; =4: TRAN ; =5: SYME
*         =6: ALBE ; =7: ZERO ; =20: CYLI.
* ICODE   physical albedo index on each side of the domain.
* ZCODE   albedo relative to each side of the domain.
* ISPLTX  mesh-splitting data for parallelepipeds along the X-axis
*         negative value is used for equal-volume splitting of tubes.
* ISPLTY  mesh-splitting data for parallelepipeds along the Y-axis.
* ISPLTZ  mesh-splitting data for parallelepipeds along the Z-axis.
* ISPLTH  mesh-splitting index for hexagons into triangles.
* ISPLTL  mesh-splitting index for hexagons into lozenges.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPGEOM
      INTEGER MAXX,MAXY,MAXZ,MAXPTS,IHEX,IR,IMPX,LX,LY,LZ,MAT(MAXPTS),
     1 NMBLK,NCODE(6),ICODE(6),ISPLTX(MAXX),ISPLTY(MAXY),ISPLTZ(MAXZ),
     2 ISPLTH,ISPLTL
      REAL SIDE,XXX(MAXX+1),YYY(MAXY+1),ZZZ(MAXZ+1),ZCODE(6)
      LOGICAL ILK
      INTEGER, ALLOCATABLE, DIMENSION(:) :: DPP,MX
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      LOGICAL LL1,LL2,LCYL,SWCEN,EMPTY,LCM
      CHARACTER HSMG*131,GEONAM*12,TEXT12*12
      INTEGER ISTATE(NSTATE)
      EQUIVALENCE (ITYPE,ISTATE(1)),(LR1,ISTATE(2)),(LX1,ISTATE(3)),
     1 (LY1,ISTATE(4)),(LZ1,ISTATE(5))
*
      IHEX=0
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      IF((ITYPE.EQ.8).OR.(ITYPE.EQ.9)) CALL LCMGET(IPGEOM,'IHEX',IHEX)
      IF((ISTATE(8).NE.0).OR.(ISTATE(9).NE.0).OR.(ISTATE(10).NE.0).OR.
     1 (ISTATE(13).NE.0)) CALL XABORT('READ3D: UNABLE TO PROCESS THE G'
     2 //'EOMETRY.')
      LCYL=(ITYPE.EQ.3).OR.(ITYPE.EQ.4).OR.(ITYPE.EQ.6)
      IDIM=1
      IF((ITYPE.EQ.5).OR.(ITYPE.EQ.6).OR.(ITYPE.EQ.8)) IDIM=2
      IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) IDIM=3
*----
*  RECOVER THE BOUNDARY CONDITIONS
*----
      CALL LCMGET(IPGEOM,'NCODE',NCODE)
      CALL LCMGET(IPGEOM,'ZCODE',ZCODE)
      CALL LCMGET(IPGEOM,'ICODE',ICODE)
      DO 10 I=1,6
      IF(NCODE(I).EQ.10) NCODE(I)=2
      IF(NCODE(I).EQ.2) ZCODE(I)=1.0
      IF(NCODE(I).EQ.6) NCODE(I)=1
      IF((NCODE(I).EQ.20).AND.(ITYPE.NE.5).AND.(ITYPE.NE.7)) CALL
     1 XABORT('READ3D: CYLINDRICAL CORRECTION IS LIMITED TO CARTESIAN '
     2 //'GEOMETRIES.')
      IF((NCODE(I).GE.8).AND.(NCODE(I).NE.20)) THEN
         CALL XABORT('READ3D: INVALID TYPE OF B.C.')
      ENDIF
   10 CONTINUE
*----
*  CHECK COHERENCE OF THE CYLINDRICAL EXTERNAL B.C.
*----
      SWCEN=.FALSE.
      ALBMAX=-1.0E35
      ALBMIN=+1.0E35
      DO 15 IC=1,6
      IF(NCODE(IC).NE.20) GO TO 15
      SWCEN=.TRUE.
      IF(ZCODE(IC).LT.ALBMIN) ALBMIN=ZCODE(IC)
      IF(ZCODE(IC).GT.ALBMAX) ALBMAX=ZCODE(IC)
   15 CONTINUE
      IF(SWCEN.AND.(ALBMIN.NE.ALBMAX)) CALL XABORT('READ3D: CYLINDRICA'
     1 //'L IMBEDDED EXTERNAL GEOMETRY: ALBEDOS ARE INCONSISTENT.')
*
      IF(ITYPE.GE.8) THEN
         IF((NCODE(2).NE.0).OR.(NCODE(3).NE.0).OR.(NCODE(4).NE.0))
     1   CALL XABORT('READ3D: INVALID TYPE OF HEXAGONAL B.C.')
         IF(NCODE(1).EQ.5) THEN
            IF(IHEX.EQ.1) THEN
               IHEX=10
            ELSE IF(IHEX.EQ.2) THEN
               IHEX=11
            ELSE
               CALL XABORT('READ3D: BOUNDARY CONDITION HBC WITH OPTION'
     1         //' SYME IS ONLY PERMITTED WITH S30 OR SA60 SYMMETRY.')
            ENDIF
         ELSE IF((NCODE(1).GT.2).AND.(NCODE(1).NE.7)) THEN
            CALL XABORT('READ3D: BOUNDARY CONDITION HBC CAN ONLY BE US'
     1      //'ED WITH OPTIONS VOID, REFL, SYME, ALBE OR ZERO.')
         ENDIF
      ENDIF
*----
*  RECOVER THE MIXTURE NUMBERS
*----
      IF(ISTATE(6).GT.MAXPTS) THEN
         WRITE (HSMG,690) 'NMBLK',ISTATE(6),'MAXPTS',MAXPTS
         CALL XABORT(HSMG)
      ENDIF
      CALL LCMGET(IPGEOM,'MIX',MAT)
      IR=0
      DO 20 I=1,ISTATE(6)
      IR=MAX(IR,MAT(I))
   20 CONTINUE
*----
*  RECOVER THE MESH COORDINATES.
*----
      IF(LCYL.AND.(LR1.GT.MAXX)) THEN
         WRITE (HSMG,690) 'LX',LR1,'MAXX',MAXX
         CALL XABORT(HSMG)
      ELSE IF(LX1.GT.MAXX) THEN
         WRITE (HSMG,690) 'LX',LX1,'MAXX',MAXX
         CALL XABORT(HSMG)
      ENDIF
      IF(LY1.GT.MAXY) THEN
         WRITE (HSMG,690) 'LY',LY1,'MAXY',MAXY
         CALL XABORT(HSMG)
      ENDIF
      IF(LZ1.GT.MAXZ) THEN
         WRITE (HSMG,690) 'LZ',LZ1,'MAXZ',MAXZ
         CALL XABORT(HSMG)
      ENDIF
      LL1=.FALSE.
      LL2=.FALSE.
      LY=1
      YYY(1)=0.0
      YYY(2)=1.0
      LZ=1
      ZZZ(1)=0.0
      ZZZ(2)=1.0
      IF(ITYPE.EQ.2) THEN
*        1-D CARTESIAN GEOMETRY.
         LX=LX1
         NMBLK=LX
         IF((NCODE(1).EQ.0).OR.(NCODE(2).EQ.0)) GO TO 610
         CALL LCMGET(IPGEOM,'MESHX',XXX)
      ELSE IF((ITYPE.EQ.3).OR.(ITYPE.EQ.4)) THEN
*        1-D CYLINDRICAL/SPHERICAL GEOMETRY.
         LX=LR1
         NMBLK=LX
         IF(NCODE(1).NE.0) GO TO 640
         IF(NCODE(2).EQ.0) GO TO 610
         NCODE(1)=2
         CALL LCMGET(IPGEOM,'RADIUS',XXX)
      ELSE IF(ITYPE.EQ.5) THEN
*        2-D CARTESIAN GEOMETRY.
         LX=LX1
         LY=LY1
         NMBLK=LX*LY
         I2=0
         DO 30 IC=1,4
         IF(NCODE(IC).EQ.0) GO TO 610
         IF(NCODE(IC).EQ.3) I2=I2+1
   30    CONTINUE
         IF(I2.NE.0) THEN
            IF((I2.NE.2).OR.(LX.NE.LY)) GO TO 630
            NMBLK=(LX+1)*LX/2
            LL1=(NCODE(2).EQ.3).AND.(NCODE(3).EQ.3)
            LL2=(NCODE(1).EQ.3).AND.(NCODE(4).EQ.3)
            IF((.NOT.LL1).AND.(.NOT.LL2)) GO TO 620
         ENDIF
         CALL LCMGET(IPGEOM,'MESHX',XXX)
         IF(LL1.OR.LL2) THEN
            CALL LCMGET(IPGEOM,'MESHX',YYY)
         ELSE
            CALL LCMGET(IPGEOM,'MESHY',YYY)
         ENDIF
      ELSE IF(ITYPE.EQ.6) THEN
*        2-D CYLINDRICAL GEOMETRY.
         LX=LR1
         LZ=LZ1
         NMBLK=LX*LZ
         IF(NCODE(1).NE.0) GO TO 650
         IF((NCODE(2).EQ.3).OR.(NCODE(3).EQ.3).OR.(NCODE(4).EQ.3))
     1   GO TO 660
         IF((NCODE(2).EQ.0).OR.(NCODE(5).EQ.0).OR.(NCODE(6).EQ.0))
     1   GO TO 610
         NCODE(1)=2
         CALL LCMGET(IPGEOM,'RADIUS',XXX)
         CALL LCMGET(IPGEOM,'MESHZ',ZZZ)
      ELSE IF(ITYPE.EQ.7) THEN
*        3-D CARTESIAN GEOMETRY.
         LX=LX1
         LY=LY1
         LZ=LZ1
         NMBLK=LX*LY*LZ
         I2=0
         DO 40 IC=1,4
         IF(NCODE(IC).EQ.0) GO TO 610
         IF(NCODE(IC).EQ.3) I2=I2+1
   40    CONTINUE
         IF(I2.NE.0) THEN
            IF((I2.NE.2).OR.(LX.NE.LY)) GO TO 630
            NMBLK=((LX+1)*LX/2)*LZ
            LL1=(NCODE(2).EQ.3).AND.(NCODE(3).EQ.3)
            LL2=(NCODE(1).EQ.3).AND.(NCODE(4).EQ.3)
            IF((.NOT.LL1).AND.(.NOT.LL2)) GO TO 620
         ENDIF
         CALL LCMGET(IPGEOM,'MESHX',XXX)
         IF(LL1.OR.LL2) THEN
            CALL LCMGET(IPGEOM,'MESHX',YYY)
         ELSE
            CALL LCMGET(IPGEOM,'MESHY',YYY)
         ENDIF
         CALL LCMGET(IPGEOM,'MESHZ',ZZZ)
      ELSE IF(ITYPE.EQ.8) THEN
*        2-D HEXAGONAL GEOMETRY.
         LX=LX1
         NMBLK=LX
         CALL LCMGET(IPGEOM,'SIDE',SIDE)
      ELSE IF(ITYPE.EQ.9) THEN
*        3-D HEXAGONAL GEOMETRY.
         LX=LX1
         LZ=LZ1
         NMBLK=LX*LZ
         CALL LCMGET(IPGEOM,'SIDE',SIDE)
         CALL LCMGET(IPGEOM,'MESHZ',ZZZ)
      ELSE
         CALL XABORT('READ3D: INVALID TYPE OF GEOMETRY.')
      ENDIF
      IF(NMBLK.NE.ISTATE(6)) THEN
         WRITE(HSMG,'(45HREAD3D: INVALID NUMBER OF REGIONS. NUMBER OF ,
     1   13HMIX ENTRIES =,I7,20H NUMBER OF REGIONS =,I7)') ISTATE(6),
     2   NMBLK
         CALL XABORT(HSMG)
      ENDIF
      DO 50 IC=1,6,2
      IF((NCODE(IC).EQ.4).AND.(NCODE(IC+1).NE.4)) GO TO 670
   50 CONTINUE
      IF((NCODE(2).EQ.3).AND.(NCODE(3).EQ.3)) THEN
         ZCODE(3)=ZCODE(1)
         ZCODE(2)=ZCODE(4)
      ELSE IF((NCODE(1).EQ.3).AND.(NCODE(4).EQ.3)) THEN
         ZCODE(1)=ZCODE(3)
         ZCODE(4)=ZCODE(2)
      ENDIF
*----
*  UNFOLD GEOMETRY IF HEXAGONAL IN LOZENGES
*----
      ISPLTL=0
      ISPLTH=0
      CALL LCMLEN(IPGEOM,'SPLITL',ILEN,ITYLCM)
      IF(ILEN.GT.0) CALL LCMGET(IPGEOM,'SPLITL',ISPLTL)
      CALL LCMLEN(IPGEOM,'SPLITH',ILEN,ITYLCM)
      IF(ILEN.GT.0) CALL LCMGET(IPGEOM,'SPLITH',ISPLTH)
      IF((ISPLTL.GT.0).AND.(IHEX.NE.9)) THEN
         ALLOCATE(DPP(MAXPTS),MX(LX*LZ))
         DO 150 I=1,LX*LZ
         MX(I)=MAT(I)
  150    CONTINUE
         LXOLD=LX
         CALL BIVALL(MAXPTS,IHEX,LXOLD,LX,DPP)
         DO 165 KZ=1,LZ
         DO 160 KX=1,LX
         KEL=DPP(KX)+(KZ-1)*LXOLD
         MAT(KX+(KZ-1)*LX)=MX(KEL)
  160    CONTINUE
  165    CONTINUE
         DEALLOCATE(DPP,MX)
         IHEX=9
      ENDIF
*----
*  MESH-SPLITTING
*----
      IF(ISTATE(11).NE.0) THEN
         CALL LCMLEN(IPGEOM,'SPLITR',ILEN1,ITYLCM)
         CALL LCMLEN(IPGEOM,'SPLITX',ILEN2,ITYLCM)
         IF(LCYL.AND.(ILEN1.GT.0)) THEN
            CALL LCMGET(IPGEOM,'SPLITR',ISPLTX)
         ELSE IF(ILEN2.GT.0) THEN
            CALL LCMGET(IPGEOM,'SPLITX',ISPLTX)
         ELSE
            DO 60 I=1,LX
            ISPLTX(I)=1
   60       CONTINUE
         ENDIF
         CALL LCMLEN(IPGEOM,'SPLITY',ILEN,ITYLCM)
         IF(ILEN.GT.0) THEN
            CALL LCMGET(IPGEOM,'SPLITY',ISPLTY)
         ELSE IF(LL1.OR.LL2) THEN
            DO 65 I=1,LX
            ISPLTY(I)=ISPLTX(I)
   65       CONTINUE
         ELSE
            DO 70 I=1,LY
            ISPLTY(I)=1
   70       CONTINUE
         ENDIF
         CALL LCMLEN(IPGEOM,'SPLITZ',ILEN,ITYLCM)
         IF(ILEN.GT.0) THEN
            CALL LCMGET(IPGEOM,'SPLITZ',ISPLTZ)
         ELSE
            DO 80 I=1,LZ
            ISPLTZ(I)=1
   80       CONTINUE
         ENDIF
         IF((ISPLTH.GT.0).AND.(ISPLTL.GT.0)) THEN
            CALL XABORT('READ3D: SPLITH AND SPLITL KEYWORDS ARE EXCLUS'
     1      //'IVE.')
         ENDIF
         CALL SPLIT0(MAXPTS,ITYPE,NCODE,LXOLD,LYOLD,LZOLD,ISPLTX,ISPLTY,
     1   ISPLTZ,0,ISPLTL,NMBLK,LX,LY,LZ,SIDE,XXX,YYY,ZZZ,MAT,.TRUE.,
     2   IMPX)
         IF(NMBLK.GT.MAXPTS) THEN
            WRITE (HSMG,690) 'NMBLK',NMBLK,'MAXPTS',MAXPTS
            CALL XABORT(HSMG)
         ENDIF
      ENDIF
*
      ILK=((NCODE(1).EQ.1).AND.(ZCODE(1).NE.1.0)).OR.(NCODE(1).EQ.7).OR.
     1    ((NCODE(2).EQ.1).AND.(ZCODE(2).NE.1.0)).OR.(NCODE(2).EQ.7).OR.
     2    ((NCODE(3).EQ.1).AND.(ZCODE(3).NE.1.0)).OR.(NCODE(3).EQ.7).OR.
     3    ((NCODE(4).EQ.1).AND.(ZCODE(4).NE.1.0)).OR.(NCODE(4).EQ.7).OR.
     4    ((NCODE(5).EQ.1).AND.(ZCODE(5).NE.1.0)).OR.(NCODE(5).EQ.7).OR.
     5    ((NCODE(6).EQ.1).AND.(ZCODE(6).NE.1.0)).OR.(NCODE(6).EQ.7).OR.
     6    ((NCODE(1).EQ.8).AND.(ZCODE(1).NE.1.0)).OR.
     7    ((NCODE(2).EQ.8).AND.(ZCODE(2).NE.1.0)).OR.
     8    ((NCODE(3).EQ.8).AND.(ZCODE(3).NE.1.0)).OR.
     9    ((NCODE(4).EQ.8).AND.(ZCODE(4).NE.1.0)).OR.
     1    ((NCODE(5).EQ.8).AND.(ZCODE(5).NE.1.0)).OR.
     2    ((NCODE(6).EQ.8).AND.(ZCODE(6).NE.1.0))
      IF(IMPX.GT.0) THEN
         IF(ITYPE.EQ.2) THEN
            WRITE (6,'(/19H 1-D SLAB GEOMETRY.)')
         ELSE IF(ITYPE.EQ.3) THEN
            WRITE (6,'(/26H 1-D CYLINDRICAL GEOMETRY.)')
         ELSE IF(ITYPE.EQ.4) THEN
            WRITE (6,'(/24H 1-D SPHERICAL GEOMETRY.)')
         ELSE IF(ITYPE.EQ.5) THEN
            WRITE (6,'(/24H 2-D CARTESIAN GEOMETRY.)')
         ELSE IF(ITYPE.EQ.6) THEN
            WRITE (6,'(/18H 2-D R-Z GEOMETRY.)')
         ELSE IF(ITYPE.EQ.7) THEN
            WRITE (6,'(/24H 3-D CARTESIAN GEOMETRY.)')
         ELSE IF(ITYPE.EQ.8) THEN
            WRITE (6,'(/24H 2-D HEXAGONAL GEOMETRY.)')
         ELSE IF(ITYPE.EQ.9) THEN
            WRITE (6,'(/24H 3-D HEXAGONAL GEOMETRY.)')
         ENDIF
         CALL LCMINF(IPGEOM,GEONAM,TEXT12,EMPTY,ILONG,LCM)
         WRITE (6,'(1H+,26X,18HBASED ON GEOMETRY ,A12,1H./)') GEONAM
         WRITE (6,770) LX,MAXX,LY,MAXY,LZ,MAXZ,IR
         IF(.NOT.ILK) WRITE (6,'(17H INFINITE DOMAIN./)')
      ENDIF
      RETURN
*
  610 CALL XABORT('READ3D: A BOUNDARY CONDITION IS MISSING.')
  620 CALL XABORT('READ3D: THE DIAGONAL CONDITIONS X+: DIAG Y-: DIAG A'
     1 //'ND X-: DIAG Y+: DIAG ARE THE ONLY PERMITTED.')
  630 CALL XABORT('READ3D: LX=LY WITH A DIAGONAL SYMMETRY.')
  640 CALL XABORT('READ3D: CYLINDRICAL GEOMETRY - ONLY THE R+: BOUNDAR'
     1 //'Y CONDITION IS REQUIRED.')
  650 CALL XABORT('READ3D: CYLINDRICAL GEOMETRY - ONLY THE R+:, Z-: AN'
     1 //'D Z+: BOUNDARY CONDITIONS ARE REQUIRED.')
  660 CALL XABORT('READ3D: CYLINDRICAL GEOMETRY : THE DIAG BOUNDARY CO'
     1 //'NDITION CANNOT BE USED.')
  670 CALL XABORT('READ3D: THE TRANSLATION CONDITIONS X-: TRAN X+: TRA'
     1 //'N, Y-: TRAN Y+: TRAN AND Z-: TRAN Z+: TRAN ARE THE ONLY PERM'
     1 //'ITTED.')
*
  690 FORMAT (29HREAD3D: INSUFFICIENT STORAGE.,5X,A6,1H=,I7,8H ; AVAIL,
     1 13HABLE STORAGE ,A6,1H=,I7)
  770 FORMAT (/44H NUMBER OF MESH INTERVALS ALONG THE X AXIS =,I5,5X,
     1 24HAVAILABLE STORAGE MAXX =,I7/26X,18HALONG THE Y AXIS =,I5,5X,
     2 24HAVAILABLE STORAGE MAXY =,I7/26X,18HALONG THE Z AXIS =,I5,5X,
     3 24HAVAILABLE STORAGE MAXZ =,I7/28H NUMBER OF DISTINCT MIXTURES,
     4 2H =,I7/)
      END
