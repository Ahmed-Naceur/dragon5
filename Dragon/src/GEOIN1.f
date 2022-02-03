*DECK GEOIN1
      RECURSIVE SUBROUTINE GEOIN1 (GEONAM,IPLIST,LEVEL,IMPX,MAXMIX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read and/or modify an object oriented geometry.
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
* GEONAM  name of the directory where the geometry is stored.
* IPLIST  pointer to the geometry LCM object (L_GEOM signature).
* LEVEL   hierarchical level of the geometry.
* IMPX    print flag (IMPX=0 for no print).
*
*Parameters: output
* MAXMIX  maximum number of mixtures, considering all sub-geometries.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIST
      INTEGER LEVEL,IMPX,MAXMIX
      CHARACTER GEONAM*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXCOD=21,MAXHEX=9,MAXTEX=4,MAXTUR=12,MAXTYP=30,
     1 MXCL=500,NSTATE=40,IOUT=6)
      LOGICAL LHEX,LTRI,EMPTY,LCM,SWANG
      LOGICAL LTOT,LHOT,LCOUR
      CHARACTER NAMT*12,COND(MAXCOD)*4,CHEX(MAXHEX)*8,CHET(MAXTEX)*8,
     1 CTUR(MAXTUR)*1,TYPE(0:MAXTYP)*16,TEXT4*4,CARLIR*12,TEXT12*12,
     2 DIR*1
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISTATE,JSTATE,NCODE,ICODE
      REAL, ALLOCATABLE, DIMENSION(:) :: ZCODE
      DOUBLE PRECISION DBLLIR,DREALIR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: CELL,MIX,ISECT,IMESH,MERGE,
     1 ITURN,NS,MILIE,MIXDL,MIXGR
      REAL, ALLOCATABLE, DIMENSION(:) :: MESH,CYL,CENT,XR0,RR0,ANG,
     1 ARPIN,RS,FRACT,POURC,PROCE
      INTEGER IMIXHT,NRINGH,NZ,NAP1,NSETM,KREG,JREG,ITRI,IAN
*----
*  Data
*----
      SAVE  COND,CHEX,CTUR,TYPE
      DATA  COND
     >      /'VOID','REFL','DIAG','TRAN','SYME',
     >       'ALBE','ZERO','PI/2','PI'  ,'SSYM',
     >        9*' ','CYLI','ACYL'/
      DATA  CHEX
     >      /'S30     ','SA60    ','SB60    ','S90     ','R120    ',
     >       'R180    ','SA180   ','SB180   ','COMPLETE'/
      DATA  CTUR
     >      /'A','B','C','D','E','F','G','H','I','J','K','L'/
      DATA  TYPE
     >      /'VIRTUAL         ','HOMOGENEOUS     ','CARTESIAN 1-D   ',
     >       'TUBE 1-D        ','SPHERE 1-D      ','CARTESIAN 2-D   ',
     >       'TUBE 2-D (Z)    ','CARTESIAN 3-D   ','HEXAGONAL 2-D   ',
     >       'HEXAGONE 3-D (Z)','TUBE 2-D (X)    ','TUBE 2-D (Y)    ',
     >       'HEX/TRIANGLE 2D ','HEX/TRIANGLE 3D ','                ',
     >       'R-THETA         ','TRIANGULAR 2-D  ','TRIANGULAR 3-D  ',
     >       '                ','                ','2-D RECT. CELL  ',
     >       '3-D RECT. CELL X','3-D RECT. CELL Y','3-D RECT. CELL Z',
     >       '2-D HEX. CELL   ','3-D HEX. CELL Z ','2-D HEXT CELL   ',
     >       '3-D HEXT CELL Z ','                ','                ',
     >       'DO-IT-YOURSELF  '/
*
      ALLOCATE(ISTATE(NSTATE),JSTATE(NSTATE),NCODE(6),ICODE(6))
      ALLOCATE(ZCODE(6))
      IMIXHT=0
      MINMIX=0
      MINICO=1
      NPIN=0
      IRLYZ=0
      LR=0
      LX=0
      LY=0
      LZ=0
      LREG=0
      CALL LCMLEN(IPLIST,'SIGNATURE',ILONG,ITYX)
      IF(ILONG.EQ.0) THEN
*       INPUT A NEW GEOMETRY.
        CALL XDISET(ISTATE,NSTATE,0)
        LHEX=.FALSE.
        LTRI=.FALSE.
        LCOUR=.FALSE.
        DO 20 I=1,6
        NCODE(I)=0
        ZCODE(I)=0.0
        ICODE(I)=0
   20   CONTINUE
      ELSE
*       MODIFY AN EXISTING GEOMETRY.
        CALL LCMGTC(IPLIST,'SIGNATURE',12,1,CARLIR)
        IF(CARLIR.NE.'L_GEOM') THEN
          NAMT=GEONAM
          CALL XABORT('GEOIN1: SIGNATURE OF '//NAMT//' IS '//CARLIR
     1    //'. L_GEOM EXPECTED.')
        ENDIF
        CALL LCMGET(IPLIST,'STATE-VECTOR',ISTATE)
        LR=ISTATE(2)
        LX=ISTATE(3)
        LY=ISTATE(4)
        LZ=ISTATE(5)
        LREG=ISTATE(6)
        LHEX=(ISTATE(1).EQ. 8) .OR. (ISTATE(1) .EQ. 9) .OR.
     1       (ISTATE(1).EQ.12) .OR. (ISTATE(1) .EQ.13) .OR.
     2       (ISTATE(1).EQ.24) .OR. (ISTATE(1) .EQ.25) .OR.
     3       (ISTATE(1).EQ.26) .OR. (ISTATE(1) .EQ.27)
        LTRI=(ISTATE(1).EQ.16).OR.(ISTATE(1).EQ.17)
        LCOUR=.FALSE.
        IF(LHEX) THEN
          CALL LCMGET(IPLIST,'IHEX',IHEX)
          LCOUR=IHEX.EQ.9
        ENDIF
        CALL LCMGET(IPLIST,'NCODE',NCODE)
        CALL LCMGET(IPLIST,'ZCODE',ZCODE)
        CALL LCMGET(IPLIST,'ICODE',ICODE)
        IF((LEVEL.EQ.1).AND.(ISTATE(1).EQ.0)) THEN
          GO TO 30
        ELSE IF(LEVEL.EQ.1) THEN
          GO TO 50
        ENDIF
      ENDIF
*
   30 CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      IF(ITYPLU.NE.3) CALL XABORT('GEOIN1: CHARACTER DATA EXPECTED(1).')
      IF(CARLIR.EQ.'VIRTUAL') THEN
        ISTATE(1)=0
      ELSE IF(CARLIR.EQ.'HOMOGE') THEN
        ISTATE(1)=1
        LREG=1
      ELSE IF(CARLIR.EQ.'CAR1D') THEN
        ISTATE(1)=2
        CALL REDGET(ITYPLU,LX,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        LREG=LX
      ELSE IF(CARLIR.EQ.'SPHERE') THEN
        ISTATE(1)=4
        CALL REDGET(ITYPLU,LR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        LREG=LR
      ELSE IF(CARLIR.EQ.'CAR2D') THEN
        ISTATE(1)=5
        CALL REDGET(ITYPLU,LX,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        CALL REDGET(ITYPLU,LY,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        LREG=LX*LY
      ELSE IF(CARLIR.EQ.'CAR3D') THEN
        ISTATE(1)=7
        CALL REDGET(ITYPLU,LX,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        CALL REDGET(ITYPLU,LY,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        CALL REDGET(ITYPLU,LZ,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        LREG=LX*LY*LZ
      ELSE IF(CARLIR.EQ.'HEX') THEN
        ISTATE(1)=8
        LHEX=.TRUE.
        CALL REDGET(ITYPLU,LX,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        LREG=LX
      ELSE IF(CARLIR.EQ.'HEXZ') THEN
        ISTATE(1)=9
        LHEX=.TRUE.
        CALL REDGET(ITYPLU,LX,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        CALL REDGET(ITYPLU,LZ,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        LREG=LX*LZ
      ELSE IF(CARLIR.EQ.'HEXT') THEN
        IMIXHT=1
        ISTATE(1)=12
        LHEX=.TRUE.
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        NRINGH=MAX(1,INTLIR)
        LX=NRINGH
        LREG=6*NRINGH*NRINGH
      ELSE IF(CARLIR.EQ.'HEXTZ') THEN
        IMIXHT=2
        ISTATE(1)=13
        LHEX=.TRUE.
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        NRINGH=MAX(1,INTLIR)
        LX=NRINGH
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        LZ=INTLIR
        LREG=6*NRINGH*NRINGH*LZ
      ELSE IF(CARLIR.EQ.'RTHETA') THEN
        ISTATE(1)=15
        CALL REDGET(ITYPLU,LR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        CALL REDGET(ITYPLU,LZ,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        LREG=LR*LZ
      ELSE IF(CARLIR.EQ.'TRI') THEN
        ISTATE(1)=16
        LTRI=.TRUE.
        CALL REDGET(ITYPLU,LX,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        LREG=LX
      ELSE IF(CARLIR.EQ.'TRIZ') THEN
        ISTATE(1)=17
        LTRI=.TRUE.
        CALL REDGET(ITYPLU,LX,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        CALL REDGET(ITYPLU,LZ,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        LREG=LX*LZ
      ELSE IF(CARLIR(1:4).EQ.'TUBE') THEN
        DIR=CARLIR(5:5)
        CALL REDGET(ITYPLU,LR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        IF(DIR.EQ.' ') THEN
          ISTATE(1)=3
          LX=1
          LY=1
          IRLXY=1
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.EQ.3) THEN
            IRLXY=0
          ELSE IF(ITYPLU.EQ.2) THEN
            CALL XABORT('GEOIN1: INVALID REAL DATA.')
          ELSE
            LX=INTLIR
            CALL REDGET(ITYPLU,LY,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTE'
     1      //'D.')
          ENDIF
          LREG=LR*LY*LX
          IF(IRLXY.EQ.0) GO TO 60
        ELSE
          LX=1
          LY=1
          LZ=1
          CALL REDGET(ITYPLU,LX,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.EQ.3) THEN
            IRLYZ=0
          ELSE IF(ITYPLU.EQ.2) THEN
            CALL XABORT('GEOIN1: REAL DATA NOT EXPECTED.')
          ELSE
            LY=INTLIR
            IRLYZ=1
            CALL REDGET(ITYPLU,LZ,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTE'
     1      //'D.')
          ENDIF
          LREG=LR*LY*LZ*LX
          IF(DIR.EQ.'X') THEN
            ISTATE(1)=10
            IF(IRLYZ.EQ.0) GO TO 60
          ELSE IF(DIR.EQ.'Y') THEN
            ISTATE(1)=11
            IF(IRLYZ.EQ.0) THEN
              LY=LX
              LX=1
              GO TO 60
            ENDIF
          ELSE IF(DIR.EQ.'Z') THEN
            ISTATE(1)=6
            IF(IRLYZ.EQ.0) THEN
              LZ=LX
              LX=1
              GO TO 60
            ENDIF
          ELSE
            CALL XABORT('GEOIN1: INVALID DATA IN TUBE CONSTRUCT.')
          ENDIF
        ENDIF
      ELSE IF(CARLIR(1:6).EQ.'CARCEL') THEN
        CALL REDGET(ITYPLU,LR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        DIR=CARLIR(7:7)
        IF(DIR.EQ.' ') THEN
          ISTATE(1)=20
          LX=1
          LY=1
          IRLXY=1
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.EQ.3) THEN
            IRLXY=0
          ELSE IF(ITYPLU.EQ.2) THEN
            CALL XABORT('GEOIN1: INVALID REAL DATA.')
          ELSE
            LX=INTLIR
            CALL REDGET(ITYPLU,LY,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTE'
     1      //'D.')
          ENDIF
          LREG=(LR+1)*LY*LX
          IF(IRLXY.EQ.0) GO TO 60
        ELSE
          LX=1
          LY=1
          LZ=1
          CALL REDGET(ITYPLU,LX,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.EQ.3) THEN
            IRLYZ=0
          ELSE IF(ITYPLU.EQ.2) THEN
            CALL XABORT('GEOIN1: INVALID REAL DATA.')
          ELSE
            LY=INTLIR
            IRLYZ=1
            CALL REDGET(ITYPLU,LZ,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTE'
     1      //'D.')
          ENDIF
          LREG=(LR+1)*LY*LZ*LX
          IF(DIR.EQ.'X') THEN
            ISTATE(1)=21
          ELSE IF(DIR.EQ.'Y') THEN
            ISTATE(1)=22
            IF(IRLYZ.EQ.0) THEN
              LY=LX
              LX=1
              GO TO 60
            ENDIF
          ELSE IF(DIR.EQ.'Z') THEN
            ISTATE(1)=23
            IF(IRLYZ.EQ.0) THEN
              LZ=LX
              LX=1
              GO TO 60
            ENDIF
          ELSE
            CALL XABORT('GEOIN1: INVALID DATA.')
          ENDIF
        ENDIF
      ELSE IF(CARLIR(1:6).EQ.'HEXCEL') THEN
        LHEX=.TRUE.
        CALL REDGET(ITYPLU,LR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        LX=1
        IF(CARLIR(7:7).EQ.' ') THEN
          ISTATE(1)=24
          LREG=LR+1
        ELSE IF(CARLIR(7:7).EQ.'Z') THEN
          ISTATE(1)=25
          CALL REDGET(ITYPLU,LZ,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
          LREG=(LR+1)*LZ
        ELSE
          CALL XABORT('GEOIN1: INVALID SUFFIX FOR HEXCEL.')
        ENDIF
      ELSE IF(CARLIR(1:7).EQ.'HEXTCEL') THEN
        IMIXHT=1
        LHEX=.TRUE.
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1'//
     >  ': Number of annular regions missing.')
        LR=INTLIR
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1'//
     >  ': Number of triangular crowns missing.')
        NRINGH=MAX(1,INTLIR)
        LX=NRINGH
        LREG=6*NRINGH*NRINGH*(LR+1)
        IF(CARLIR(8:8).EQ.' ') THEN
          ISTATE(1)=26
        ELSE IF(CARLIR(8:8).EQ.'Z') THEN
          IMIXHT=2
          ISTATE(1)=27
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.1) CALL XABORT('GEOIN1'//
     >    ': Number of z planes missing.')
          LZ=INTLIR
          LREG=LREG*LZ
        ELSE
          CALL XABORT('GEOIN1: INVALID SUFFIX FOR HEXTCEL.')
        ENDIF
      ELSE IF(CARLIR.EQ.'GROUP') THEN
*       DO-IT-YOURSELF OPTION.
        ISTATE(1)=30
        CALL REDGET(ITYPLU,LX,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        LREG=LX
      ELSE IF(CARLIR.EQ.':::') THEN
        GO TO 60
      ELSE IF(CARLIR.NE.GEONAM) THEN
*       COPY ATTRIBUTES FROM AN EXISTING GEOMETRY LOCATED ON A PARALLEL
*       DIRECTORY OF THE LCM OBJECT POINTED BY IPLIST.
        IF(LEVEL.EQ.1) CALL XABORT('GEOIN1: THE GEOMETRY NAME SHOULD A'
     1  //'PPEAR BEFORE THE ::.')
        CALL LCMSIX(IPLIST,' ',2)
        CALL LCMLEN(IPLIST,CARLIR,ILONG,ITYX)
        IF(ILONG.EQ.0) CALL XABORT('GEOIN1: UNKNOWN GEOMETRY.')
        CALL LCMSIX(IPLIST,CARLIR,1)
        IFILE=KDROPN('DUMMYSQ',0,2,0)
        IF(IFILE.LE.0) CALL XABORT('GEOIN1: KDROPN FAILURE.')
        CALL LCMEXP(IPLIST,0,IFILE,1,1)
        REWIND(IFILE)
        CALL LCMSIX(IPLIST,' ',2)
        CALL LCMSIX(IPLIST,GEONAM,1)
        CALL LCMEXP(IPLIST,0,IFILE,1,2)
        IRC=KDRCLS(IFILE,2)
        IF(IRC.LT.0) CALL XABORT('GEOIN1: KDRCLS FAILURE.')
        CALL LCMGET(IPLIST,'STATE-VECTOR',ISTATE)
        LR=ISTATE(2)
        LX=ISTATE(3)
        LY=ISTATE(4)
        LZ=ISTATE(5)
        LREG=ISTATE(6)
        LHEX=(ISTATE(1).EQ. 8) .OR. (ISTATE(1) .EQ. 9) .OR.
     1       (ISTATE(1).EQ.12) .OR. (ISTATE(1) .EQ.13) .OR.
     2       (ISTATE(1).EQ.24) .OR. (ISTATE(1) .EQ.25) .OR.
     3       (ISTATE(1).EQ.26) .OR. (ISTATE(1) .EQ.27)
        LTRI=(ISTATE(1).EQ.16).OR.(ISTATE(1).EQ.17)
        CALL LCMGET(IPLIST,'NCODE',NCODE)
        CALL LCMGET(IPLIST,'ZCODE',ZCODE)
        CALL LCMGET(IPLIST,'ICODE',ICODE)
      ENDIF
*
   50 CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      IF(ITYPLU.NE.3) CALL XABORT('GEOIN1: CHARACTER DATA EXPECTED(2).')
   60 IF(CARLIR.EQ.'EDIT') THEN
        CALL REDGET(ITYPLU,IMPX,REALIR,CARLIR,DREALIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
      ELSE IF((CARLIR.EQ.'MIX').OR.(CARLIR.EQ.'CELL')) THEN
*       INPUT MIXTURE NUMBERS OR FORCE SUB GEOMETRIES AT SPECIFIC
*       LOCATIONS.
        ALLOCATE(CELL(3*LREG),MIX(LREG))
        CALL XDISET(MIX,LREG,0)
        LTOT=.TRUE.
        I=0
        IKG=0
   70   I=I+1
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.EQ.3) THEN
          IF(CARLIR.EQ.'PLANE') THEN
            IF(I.EQ.1) THEN
              IF(ISTATE(1).EQ.7.OR.ISTATE(1).EQ.9) THEN
                IF(ISTATE(1).EQ.9) LY=1
                CALL GEOMIX(LX,LY,LZ,LCOUR,MIX,MINMIX,ISTATE(7))
                LTOT=.FALSE.
                GO TO 70
              ELSE
                CALL XABORT('GEOIN1: INVALID KEY WORD PLANE FOR NON '
     1          //' 3-D GEOMETRY')
              ENDIF
            ELSE
              CALL XABORT('GEOIN1: WRONG USE OF KEYWORD PLANE.')
            ENDIF
          ENDIF
          IF((CARLIR(2:2).EQ.'-').OR.(CARLIR(2:2).EQ.'+').OR.
     1    (CARLIR.EQ.'HBC').OR.(CARLIR(1:4).EQ.'MESH').OR.
     2    (CARLIR(1:5).EQ.'SPLIT').OR.(CARLIR.EQ.'SIDE').OR.
     3    (CARLIR(:3).EQ.'MIX').OR.(CARLIR.EQ.'MERGE').OR.
     4    (CARLIR.EQ.'TURN').OR.(CARLIR.EQ.'CLUSTER').OR.
     5    (CARLIR(2:4).EQ.'PIN').OR.(CARLIR.EQ.'BIHET').OR.
     6    (CARLIR.EQ.'POURCE').OR.(CARLIR.EQ.'PROCEL').OR.
     7    (CARLIR.EQ.'SECT').OR.(CARLIR.EQ.'RADIUS').OR.
     8    (CARLIR.EQ.'HMIX').OR.(CARLIR.EQ.';').OR.
     9    (CARLIR.EQ.':::')) GO TO 90
           IF(I.GT.LREG) CALL XABORT('GEOIN1: MIX/CELL INDEX OVERFLO'
     1     //'W.')
           DO 80 J=1,I-1
           JKG=-MIX(J)
           WRITE (TEXT12(:4),'(A4)') CELL(3*(JKG-1)+1)
           WRITE (TEXT12(5:8),'(A4)') CELL(3*(JKG-1)+2)
           WRITE (TEXT12(9:),'(A4)') CELL(3*(JKG-1)+3)
           IF(CARLIR.EQ.TEXT12) THEN
             MIX(I)=-JKG
             GO TO 70
           ENDIF
   80      CONTINUE
           IKG=IKG+1
           ISTATE(8)=1
           MIX(I)=-IKG
           READ (CARLIR(:4),'(A4)') CELL(3*(IKG-1)+1)
           READ (CARLIR(5:8),'(A4)') CELL(3*(IKG-1)+2)
           READ (CARLIR(9:),'(A4)') CELL(3*(IKG-1)+3)
        ELSE IF(ITYPLU.EQ.1) THEN
           IF(I.GT.LREG) CALL XABORT('GEOIN1: MIX INDEX OVERFLOW.')
           MIX(I)=INTLIR
           ISTATE(7)=MAX(ISTATE(7),MIX(I))
           MINMIX=MIN(MINMIX,MIX(I))
        ELSE
           CALL XABORT('GEOIN1: INTEGER OR CHARACTER DATA EXPECTED.')
        ENDIF
        GO TO 70
   90   CONTINUE
        IF(CARLIR.EQ.'REPEAT') THEN
           NBR=LREG/(I-1)
           NBRR=NBR*(I-1)
           IF(NBRR.NE.LREG) THEN
             WRITE(IOUT,530) I-1,LREG
             CALL XABORT('GEOIN1: IMPOSSIBLE TO REPEAT AN INTEGER NUMB'
     1       //'ER OF TIMES.')
           ENDIF
           JREP=I-1
           DO IREP=1,NBR-1
             DO II=1,I-1
               JREP=JREP+1
               MIX(JREP)=MIX(II)
             ENDDO
           ENDDO
           I=JREP+1
           CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
           IF(ITYPLU.NE.3) CALL XABORT('GEOIN1: CHARACTER DATA EXPECTE'
     1     //'D.')
        ENDIF
*-- Begin symmetric mixtures for HEXT, HEXTZ, HEXTCEL and HEXTCELZ
        IF(IMIXHT .GT. 0) THEN
          IF(I-1 .NE. LREG) THEN
            LTOT=.FALSE.
            NRINGH=LX
            NZ=MAX(LZ,1)
            NAP1=LR+1
            NSETM=NRINGH*NZ*NAP1
            IF(I-1 .EQ. NSETM) THEN
*----
*  Mixture given per hexagonal rings
*  create compatible complete mix array 
*----
              KREG=LREG
              DO IZ=NZ,1,-1
                DO IS=6,1,-1
                  DO IR=NRINGH,1,-1
                    JREG=(IZ-1)*NRINGH+IR-1
                    DO ITRI=2*IR-1,1,-1
                      DO IAN=NAP1,1,-1
                        MREG=JREG*NAP1+IAN
                        MIX(KREG)=MIX(MREG)
                        KREG=KREG-1
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
              IF(KREG.NE.0) CALL XABORT('GEOIN1'//
     >         ': Problem with mixture reformatting.')
            ELSE IF(I-1 .EQ. 6*NSETM) THEN
*----
*  Mixture given per hexagonal rings and per sector
*  Create compatible complete mix array 
*----
              KREG=LREG
              DO IZ=NZ,1,-1
                DO IS=6,1,-1
                  DO IR=NRINGH,1,-1
                    JREG=((IZ-1)*6+(IS-1))*NRINGH+IR-1
                    DO ITRI=2*IR-1,1,-1
                      DO IAN=NAP1,1,-1
                        MREG=JREG*NAP1+IAN
                        MIX(KREG)=MIX(MREG)
                        KREG=KREG-1
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
              IF(KREG.NE.0) CALL XABORT('GEOIN1'//
     >        ': Problem with mixture reformatting.')
            ELSE
              CALL XABORT('GEOIN1'//
     >        ': Number of mixtures provided is invalid.')
            ENDIF
          ENDIF
        ELSE
          ISTATE(6)=I-1
        ENDIF
*-- End symmetric mixtures for HEXT, HEXTZ, HEXTCEL and HEXTCELZ
        IF(LTOT) LREG=I-1
        IF(IKG.GT.0) CALL LCMPUT(IPLIST,'CELL',3*IKG,3,CELL)
        CALL LCMPUT(IPLIST,'MIX',LREG,1,MIX)
        DEALLOCATE(MIX,CELL)
        GO TO 60
      ELSE IF(CARLIR.EQ.'HMIX') THEN
*       INPUT MERGED MIXTURE NUMBERS for homogenization by geometry.
        ALLOCATE(MIX(LREG))
        CALL XDISET(MIX,LREG,0)
        LHOT=.TRUE.
        I=0
  540   I=I+1
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.EQ.3) GO TO 550
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        IF(I.GT.LREG) GO TO 500
        MIX(I)=INTLIR
        GO TO 540
  550   CONTINUE
        IF(CARLIR .EQ. 'REPEAT') THEN
          NBR=LREG/(I-1)
          NBRR=NBR*(I-1)
          IF(NBRR .NE. LREG ) THEN
            WRITE(IOUT,530) I-1,LREG
            CALL XABORT('GEOIN1: Impossible to repeat an integer numbe'
     1      //'r of times')
          ENDIF
          JREP=I-1
          DO IREP=1,NBR-1
            DO II=1,I-1
              JREP=JREP+1
              MIX(JREP)=MIX(II)
            ENDDO
          ENDDO
          I=JREP+1
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.3) CALL XABORT('GEOIN1: CHARACTER DATA EXPECTED')
        ENDIF
*-- Begin symmetric H-mixtures for HEXT, HEXTZ, HEXTCEL and HEXTCELZ
        IF(IMIXHT .GT. 0) THEN
          IF(I-1 .NE. LREG) THEN
            LHOT=.FALSE.
            NRINGH=LX
            NZ=MAX(LZ,1)
            NAP1=LR+1
            NSETM=NRINGH*NZ*NAP1
            IF(I-1 .EQ. NSETM) THEN
*----
*  Mixture given per hexagonal rings
*  create compatible complete mix array 
*----
              KREG=LREG
              DO IZ=NZ,1,-1
                DO IS=6,1,-1
                  DO IR=NRINGH,1,-1
                    JREG=(IZ-1)*NRINGH+IR-1
                    DO ITRI=2*IR-1,1,-1
                      DO IAN=NAP1,1,-1
                        MREG=JREG*NAP1+IAN
                        MIX(KREG)=MIX(MREG)
                        KREG=KREG-1
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
              IF(KREG.NE.0) CALL XABORT('GEOIN1'//
     >         ': Problem with mixture reformatting.')
            ELSE IF(I-1 .EQ. 6*NSETM) THEN
*----
*  Mixture given per hexagonal rings and per sector
*  Create compatible complete mix array 
*----
              KREG=LREG
              DO IZ=NZ,1,-1
                DO IS=6,1,-1
                  DO IR=NRINGH,1,-1
                    JREG=((IZ-1)*6+(IS-1))*NRINGH+IR-1
                    DO ITRI=2*IR-1,1,-1
                      DO IAN=NAP1,1,-1
                        MREG=JREG*NAP1+IAN
                        MIX(KREG)=MIX(MREG)
                        KREG=KREG-1
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
              IF(KREG.NE.0) CALL XABORT('GEOIN1'//
     >        ': Problem with mixture reformatting.')
            ELSE
              CALL XABORT('GEOIN1'//
     >        ': Number of mixtures provided is invalid.')
            ENDIF
          ENDIF
        ENDIF
*-- End symmetric H-mixtures for HEXT, HEXTZ, HEXTCEL and HEXTCELZ
        IF(LHOT) LREG=I-1
        CALL LCMPUT(IPLIST,'HMIX',LREG,1,MIX)
        DEALLOCATE(MIX)
        GO TO 60
      ELSE IF(CARLIR(1:4).EQ.'MESH') THEN
*       INPUT CARTESIAN COORDINATES.
        IF(CARLIR(5:5).EQ.'X') THEN
          IF(LX.EQ.0) CALL XABORT('GEOIN1: MESHX - LX=0.')
          LMESH=LX+1
        ELSE IF(CARLIR(5:5).EQ.'Y') THEN
          IF(LY.EQ.0) CALL XABORT('GEOIN1: MESHY - LY=0.')
          LMESH=LY+1
        ELSE IF(CARLIR(5:5).EQ.'Z') THEN
          IF(LZ.EQ.0) CALL XABORT('GEOIN1: MESHZ - LZ=0.')
          LMESH=LZ+1
        ELSE
          CALL XABORT('GEOIN1: INVALID MESH SUFFIX.')
        ENDIF
        ALLOCATE(MESH(LMESH))
        DO 100 I=1,LMESH
        CALL REDGET(ITYPLU,INTLIR,MESH(I),TEXT12,DBLLIR)
        IF(ITYPLU.NE.2) CALL XABORT('GEOIN1: REAL DATA EXPECTED.(1)')
        IF(I.GT.1) THEN
          IF(MESH(I).LE.MESH(I-1)) THEN
            CALL XABORT('GEOIN1: NON INCREASING MESHES.')
          ENDIF
        ENDIF
  100   CONTINUE
        CALL LCMPUT(IPLIST,CARLIR,LMESH,2,MESH)
        DEALLOCATE(MESH)
      ELSE IF(CARLIR.EQ.'RADIUS') THEN
*       INPUT TUBE RADIUS.
        IF(LR.EQ.0) CALL XABORT('GEOIN1: RADIUS WITH LR=0.')
        LCYL=LR+1
        ALLOCATE(CYL(LCYL))
        DO 110 I=1,LCYL
        CALL REDGET(ITYPLU,INTLIR,CYL(I),CARLIR,DBLLIR)
        IF(ITYPLU.NE.2) CALL XABORT('GEOIN1: REAL DATA EXPECTED.(2)')
        IF(I.GT.1) THEN
          IF(CYL(I).LE.CYL(I-1)) THEN
            CALL XABORT('GEOIN1: NON INCREASING RADII.')
          ENDIF
        ENDIF
  110   CONTINUE
        IF(CYL(1).NE.0.0) CALL XABORT('GEOIN1: INVALID FIRST RADI'
     1  //'US.')
        CALL LCMPUT(IPLIST,'RADIUS',LCYL,2,CYL)
        DEALLOCATE(CYL)
      ELSE IF(CARLIR.EQ.'OFFCENTER') THEN
*       INPUT TUBE CENTER LOCATION (USE FOR CARCEL* ONLY).
        IF(LR.EQ.0) CALL XABORT('GEOIN1: OFFCENTER WITH LR=0.')
        ALLOCATE(CENT(3))
        CALL XDRSET(CENT,3,0.0)
        DO 120 I=1,3
         CALL REDGET(ITYPLU,INTLIR,CENT(I),CARLIR,DBLLIR)
         IF(ITYPLU.NE.2) GO TO 130
  120   CONTINUE
        CALL REDGET(ITYPLU,INTLIR,CENT(I),CARLIR,DBLLIR)
  130   CALL LCMPUT(IPLIST,'OFFCENTER',3,2,CENT)
        DEALLOCATE(CENT)
        GO TO 60
      ELSE IF(CARLIR.EQ.'SIDE') THEN
*        INPUT THE SIDE LENGTH IN TRIANGULAR OR HEXAGONAL GEOMETRY.
         IF((.NOT.LHEX).AND.(.NOT.LTRI)) CALL XABORT('GEOIN1: SIDE PRO'
     1   //'HIBITED.')
         CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
         IF(ITYPLU.NE.2) CALL XABORT('GEOIN1: REAL DATA EXPECTED.(3)')
         SIDE=REALIR
         CALL LCMPUT(IPLIST,'SIDE',1,2,SIDE)
         IF((ISTATE(1).EQ.12).OR.(ISTATE(1).EQ.13).OR.
     1      (ISTATE(1).EQ.26).OR.(ISTATE(1).EQ.27) ) THEN
*        Hexagonal mesh for HEXT and HEXTZ
           SIDET=SIDE/FLOAT(LX)
           CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
           IF(ITYPLU.EQ.1) THEN
             CALL XABORT('GEOIN1: INVALID INTEGER DATA .')
           ELSE IF(ITYPLU.EQ.2) THEN
             SIDET=REALIR
             HEXMMX=SIDET*FLOAT(LX)-SIDE
             IF((HEXMMX.LT.0.0).OR.(HEXMMX.GT.SIDET)) THEN
               CALL XABORT('GEOIN1: Invalid sidet, nringh or side ->'//
     1         ' 0 <= sidet*nringh-side <= sidet required')
             ENDIF
           ENDIF
           CALL LCMPUT(IPLIST,'SIDET',1,2,SIDET)
           IF(ITYPLU.EQ.3) GO TO 60
         ENDIF
      ELSE IF(CARLIR.EQ.'SECT') THEN
*       INPUT THE TYPE OF SECTORIZATION.
        IF(ISTATE(1).LT.20) CALL XABORT('GEOIN1: SECT PROHIBITED.')
        CALL REDGET(ITYPLU,ISTATE(14),REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) THEN
          ISTATE(15)=0
        ELSE
          IF((INTLIR.LT.0).OR.(INTLIR.GE.LREG)) CALL XABORT('GEOIN1: I'
     >    //'NVALID VALUE FOR jsect.')
          ISTATE(15)=INTLIR
        ENDIF
        IF((.NOT.LHEX).AND.(ISTATE(14).EQ.-1)) THEN
*         X-TYPE SECTORIZATION IN CARTESIAN CELL.
          LREG=4*LREG-3*ISTATE(15)
        ELSE IF((.NOT.LHEX).AND.(ISTATE(14).EQ.0  )) THEN
*         NO SECTORIZATION IN CARTESIAN CELL.
          LREG=LREG
        ELSE IF((.NOT.LHEX).AND.(ISTATE(14).EQ.1)) THEN
*         +-TYPE SECTORIZATION IN CARTESIAN CELL.
          LREG=4*LREG-3*ISTATE(15)
        ELSE IF((.NOT.LHEX).AND.(ISTATE(14).EQ.2)) THEN
*         +-TYPE SECTORIZATION IN CARTESIAN CELL.
          LREG=8*LREG-7*ISTATE(15)
        ELSE IF((.NOT.LHEX).AND.(ISTATE(14).EQ.3)) THEN
*         SHIFTED + AND X-TYPE SECTORIZATION IN CARTESIAN CELL.
          LREG=8*LREG-7*ISTATE(15)
        ELSE IF((.NOT.LHEX).AND.(ISTATE(14).EQ.4)) THEN
*         FULL WINDMILL.
          LREG=4+8*LREG-7*ISTATE(15)
        ELSE IF(LHEX.AND.(ISTATE(14).EQ.-1)) THEN
*         X-TYPE SECTORIZATION IN HEXAGONAL CELL.
          ALLOCATE(ISECT(LREG))
          CALL XDISET(ISECT,LREG,2)
          CALL LCMPUT(IPLIST,'SECTOR',LREG,1,ISECT)
          DEALLOCATE(ISECT)
          LREG=6*LREG-5*ISTATE(15)
        ELSE IF(ISTATE(14).NE.-999) THEN
          CALL XABORT('GEOIN1: INVALID TYPE OF SECTORIZATION.')
        ENDIF
        IF(ITYPLU.NE.1) GO TO 60
      ELSE IF (CARLIR.EQ.'RADS') THEN
*       OPTIONS FOR CYLINDRICAL CORRECTION IN CARTESIAN GEOMETRY.
        IF((ISTATE(1).NE.5).AND.(ISTATE(1).NE.7)) CALL XABORT('GEO'
     1  //'IN1: OPTION RADS IS LIMITED TO CARTESIAN GEOMETRIES.')
        CALL REDGET(INDIC,NR0,REALIR,TEXT4,DREALIR)
        SWANG=TEXT4.EQ.'ANG'
        IF(SWANG) CALL REDGET(INDIC,NR0,REALIR,TEXT4,DREALIR)
        IF(INDIC.NE.1) CALL XABORT('GEO: INTEGER DATA EXPECTED.')
        IF(NR0.EQ.0) CALL XABORT('GEOIN1: NON-ZERO INTEGER EXPECTED.')
        ALLOCATE(XR0(NR0),RR0(NR0),ANG(NR0))
        DO 135 I=1,NR0
        CALL REDGET(INDIC,INTLIR,XR0(I),TEXT4,DREALIR)
        IF(INDIC.NE.2) CALL XABORT('GEOIN1: REAL DATA EXPECTED.(4)')
        CALL REDGET(INDIC,INTLIR,RR0(I),TEXT4,DREALIR)
        IF(INDIC.NE.2) CALL XABORT('GEOIN1: REAL DATA EXPECTED.(5)')
        IF(SWANG) THEN
          CALL REDGET(INDIC,INTLIR,ANG(I),TEXT4,DREALIR)
          IF(INDIC.NE.2) CALL XABORT('GEOIN1: REAL DATA EXPECTED.(6)')
        ELSE
*         USE PI/2 + 0.1
          ANG(I)=1.670796327
        ENDIF
  135   CONTINUE
        CALL LCMPUT(IPLIST,'XR0',NR0,2,XR0)
        CALL LCMPUT(IPLIST,'RR0',NR0,2,RR0)
        CALL LCMPUT(IPLIST,'ANG',NR0,2,ANG)
        DEALLOCATE(ANG,RR0,XR0)
      ELSE IF(CARLIR(1:5).EQ.'SPLIT') THEN
*       INPUT MESH SPLITTING FACTORS.
        ISTATE(11)=1
        IF(CARLIR(6:6).EQ.'X') THEN
          IF(LX.EQ.0) CALL XABORT('GEOIN1: SPLITX - LX=0.')
          LMESH=LX
        ELSE IF(CARLIR(6:6).EQ.'Y') THEN
          IF(LY.EQ.0) CALL XABORT('GEOIN1: SPLITY - LY=0.')
          LMESH=LY
        ELSE IF(CARLIR(6:6).EQ.'Z') THEN
          IF(LZ.EQ.0) CALL XABORT('GEOIN1: SPLITZ - LZ=0.')
          LMESH=LZ
        ELSE IF(CARLIR(6:6).EQ.'R') THEN
          IF(LR.EQ.0) CALL XABORT('GEOIN1: SPLITR - LR=0.')
          LMESH=LR
        ELSE IF(CARLIR(6:6).EQ.'H') THEN
          IF(LX.EQ.0) CALL XABORT('GEOIN1: SPLITH - LX=0.')
          LMESH=1
        ELSE IF(CARLIR(6:6).EQ.'L') THEN
          IF(LX.EQ.0) CALL XABORT('GEOIN1: SPLITL - LX=0.')
          LMESH=1
        ELSE
          CALL XABORT('GEOIN1: INVALID SPLIT SUFFIX.')
        ENDIF
        ALLOCATE(IMESH(LMESH))
        DO 140 I=1,LMESH
        CALL REDGET(ITYPLU,IMESH(I),REALIR,TEXT12,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        IF(CARLIR.EQ.'SPLITR') THEN
          IF(IMESH(I).EQ.0) THEN
            CALL XABORT('GEOIN1: INVALID MESH-SPLITTING INDEX(1).')
          ENDIF
        ELSE IF((CARLIR.EQ.'SPLITH').OR.(CARLIR.EQ.'SPLITL')) THEN
          IF(IMESH(I).LT.0) THEN
            CALL XABORT('GEOIN1: INVALID MESH-SPLITTING INDEX(2).')
          ENDIF
        ELSE
          IF(IMESH(I).LE.0) THEN
            CALL XABORT('GEOIN1: INVALID MESH-SPLITTING INDEX(3).')
          ENDIF
        ENDIF
  140   CONTINUE
        CALL LCMPUT(IPLIST,CARLIR,LMESH,1,IMESH)
        DEALLOCATE(IMESH)
      ELSE IF(CARLIR.EQ.'MERGE') THEN
*       INPUT CELL-MERGING ITYPLUES.
        ISTATE(10)=1
        ALLOCATE(MERGE(LREG))
        I=0
  150   I=I+1
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.EQ.3) GO TO 160
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        IF(I.GT.LREG) CALL XABORT('GEOIN1: MERGE INDEX OVERFLOW.')
        MERGE(I)=INTLIR
        GO TO 150
  160   LREG=I-1
        CALL LCMPUT(IPLIST,'MERGE',LREG,1,MERGE)
        DEALLOCATE(MERGE)
        GO TO 60
      ELSE IF(CARLIR.EQ.'TURN') THEN
*       INPUT ORIENTATION INFORMATION.
        ALLOCATE(ITURN(LREG))
        I=0
  170   I=I+1
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('GEOIN1: CHARACTER DATA EXPECTED.')
        DO 180 J=1,MAXTUR
        IF(CARLIR.EQ.CTUR(J)) THEN
          IF(I.GT.LREG) CALL XABORT('GEOIN1: TURN INDEX OVERFLOW(1).')
          ITURN(I)=J
          GO TO 170
        ELSE IF(CARLIR.EQ.'-'//CTUR(J)) THEN
          IF(I.GT.LREG) CALL XABORT('GEOIN1: TURN INDEX OVERFLOW(2).')
          ITURN(I)=MAXTUR+J
          GO TO 170
        ENDIF
  180   CONTINUE
        LREG=I-1
        CALL LCMPUT(IPLIST,'TURN',LREG,1,ITURN)
        DEALLOCATE(ITURN)
        GO TO 60
      ELSE IF(CARLIR.EQ.'CLUSTER') THEN
*       DEFINE CLUSTER SUB GEOMETRIES.
        ALLOCATE(CELL(3*MXCL))
        I=0
  190   I=I+1
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('GEOIN1: CHARACTER DATA EXPECTED.')
        IF((CARLIR(2:2).EQ.'-').OR.(CARLIR(2:2).EQ.'+').OR.
     1  (CARLIR.EQ.'HBC').OR.(CARLIR(1:4).EQ.'MESH').OR.(CARLIR(1:5)
     2  .EQ.'SPLIT').OR.(CARLIR.EQ.'SIDE').OR.(CARLIR(:3).EQ.'MIX').OR.
     3  (CARLIR.EQ.'CELL').OR.(CARLIR.EQ.'MERGE').OR.(CARLIR.EQ.'TURN')
     4  .OR.(CARLIR(2:4).EQ.'PIN').OR.(CARLIR.EQ.'BIHET').OR.
     5  (CARLIR.EQ.'POURCE').OR.(CARLIR.EQ.'PROCEL').OR.
     6  (CARLIR.EQ.'SECT').OR.(CARLIR.EQ.'RADIUS').OR.
     7  (CARLIR.EQ.'HMIX').OR.(CARLIR.EQ.';').OR.(CARLIR.EQ.':::'))
     8  GO TO 200
        IF(I.GT.MXCL) CALL XABORT('GEOIN1: CLUSTER INDEX OVERFLOW.')
        READ (CARLIR(:4),'(A4)') CELL(3*(I-1)+1)
        READ (CARLIR(5:8),'(A4)') CELL(3*(I-1)+2)
        READ (CARLIR(9:),'(A4)') CELL(3*(I-1)+3)
        GO TO 190
  200   CALL LCMPUT(IPLIST,'CLUSTER',3*(I-1),3,CELL)
        ISTATE(13)=I-1
        DEALLOCATE(CELL)
        GO TO 60
      ELSE IF(CARLIR(2:4).EQ.'PIN') THEN
        IF(ISTATE(1) .NE.  3 .AND. ISTATE(1) .NE.  6 .AND.
     1     ISTATE(1) .NE. 10 .AND. ISTATE(1) .NE. 11 .AND.
     2     ISTATE(1) .NE.  4) GO TO 500
        IF(CARLIR.EQ.'NPIN') THEN
          IF(NPIN.EQ.-1) CALL XABORT('GEOIN1: NPIN and DPIN cannot be '
     1    //'used simultneously')
*         INPUT NUMBER OF PINS IN CLUSTER RING.
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
          IF(INTLIR.LT.1) CALL XABORT('GEOIN1: NPIN > 0 required.')
          NPIN=INTLIR
          CALL LCMPUT(IPLIST,'NPIN',1,1,NPIN)
        ELSE IF(NAMT.EQ.'DPIN') THEN
          IF(NPIN.GE.1) CALL XABORT('GEOIN1: NPIN and DPIN cannot be u'
     1    //'sed simultneously')
*         INPUT DENSITY OF PIN IN CLUSTER.
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.2) CALL XABORT('GEOIN1: REAL DATA EXPECTED FOR '
     1    //NAMT//' KEYWORD.')
          NPIN=-1
          DPIN=REALIR
          ISTATE(18)=-1
          CALL LCMPUT(IPLIST,'NPIN',1,1,NPIN)
          CALL LCMPUT(IPLIST,'DPIN',1,2,DPIN)
        ELSE IF((CARLIR.EQ.'RPIN').OR.(CARLIR.EQ.'APIN')) THEN
          IF(NPIN.EQ.-1) THEN
            CALL XABORT('GEOIN1: RPIN and APIN not compatible with DPI'
     1      //'N')
          ELSE IF(NPIN.EQ.0) THEN
            CALL XABORT('GEOIN1: NPIN required before RPIN and APIN ar'
     1      //'e defined')
          ENDIF
          IF(ISTATE(18) .EQ. 2) CALL XABORT('GEOIN1: CPIN* cannot be m'
     1    //'ixed with RPIN and APIN')
          ISTATE(18)=1
          NAMT=CARLIR
*----
*  Allocate memory for APIN or RPIN
*----
          ALLOCATE(ARPIN(NPIN))
*         INPUT RADIUS/ANGLE OF CLUSTER RING.
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.2) CALL XABORT('GEOIN1: REAL DATA EXPECTED.(7)')
          ARPIN(1)=REALIR
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.2) THEN
            CALL LCMPUT(IPLIST,NAMT,1,2,ARPIN)
            DEALLOCATE(ARPIN)
            GO TO 60
          ENDIF
          IF(NPIN.EQ.1) CALL XABORT('GEOIN1: Only one APIN or RPIN per'
     1    //'mitted.')
          ARPIN(2)=REALIR
           DO IPIN=2,NPIN-1
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.NE.2) CALL XABORT('GEOIN1: REAL DATA EXPECTED.'
     1        //'(8)')
            ARPIN(IPIN+1)=REALIR
          ENDDO
          CALL LCMPUT(IPLIST,NAMT,NPIN,2,ARPIN)
          DEALLOCATE(ARPIN)
        ELSE IF(CARLIR(1:4).EQ.'CPIN') THEN
          IF(NPIN.EQ.-1) THEN
            CALL XABORT('GEOIN1: CPIN* not compatible with DPIN')
          ELSE IF(NPIN.EQ.0) THEN
            CALL XABORT('GEOIN1: NPIN required before CPIN* is defined')
          ENDIF
          IF(ISTATE(18).EQ.1) CALL XABORT('GEOIN1: CPIN* cannot be mix'
     1    //'ed with RPIN and APIN')
          ISTATE(18)=2
          IF(CARLIR(5:5). NE. 'X' .AND.
     1       CARLIR(5:5). NE. 'Y' .AND.
     2       CARLIR(5:5). NE. 'Z' ) THEN
               CALL XABORT('GEOIN1: Only CPINX, CPINY and CPINZ permit'
     1         //'ted -- '//CARLIR(1:5)//' provided')
             ENDIF
          NAMT=CARLIR
*----
*  Allocate memory for CPIN
*----
          ALLOCATE(ARPIN(NPIN))
*         INPUT Cartesian positions of pins.
          DO IPIN=0,NPIN-1
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.NE.2) CALL XABORT('GEOIN1: REAL CPIN DATA EXPECT'
     1      //'ED.')
            ARPIN(IPIN+1)=REALIR
          ENDDO
          CALL LCMPUT(IPLIST,NAMT,NPIN,2,ARPIN)
          DEALLOCATE(ARPIN)
        ENDIF
      ELSE IF(CARLIR.EQ.'BIHET') THEN
*       DOUBLE HETEROGENEITY OPTION.
        ISTATE(12)=1
        IF(LEVEL.NE.1) CALL XABORT('GEOIN1: BIHET DATA SHOULD BE WRI'
     1  //'TTEN ON FIRST DIRECTORY LEVEL.')
        CALL LCMSIX(IPLIST,'BIHET',1)
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('GEOIN1: CHARACTER DATA EXPECTED.')
        MICRO=0
        IF(CARLIR.EQ.'TUBE') THEN
          MICRO=3
        ELSE IF(CARLIR.EQ.'SPHE') THEN
          MICRO=4
        ELSE
          CALL XABORT('GEOIN1: PROHIBITED TYPE OF MICRO GEOMETRY.')
        ENDIF
        CALL REDGET(ITYPLU,NG,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        CALL REDGET(ITYPLU,NMILG,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
*
        ALLOCATE(NS(NG))
        NSMAX=0
        DO 210 I=1,NG
        CALL REDGET(ITYPLU,NS(I),REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
        NSMAX=MAX(NSMAX,NS(I))
  210   CONTINUE
        CALL LCMPUT(IPLIST,'NS',NG,1,NS)
*
        ALLOCATE(RS((NSMAX+1)*NG))
        DO 220 IOFJ=1,(NSMAX+1)*NG
        RS(IOFJ)=0.0
  220   CONTINUE
        DO 235 I=1,NG
        DO 230 J=1,NS(I)+1
        IOFJ=(I-1)*(NSMAX+1)+J
        CALL REDGET(ITYPLU,INTLIR,RS(IOFJ),CARLIR,DBLLIR)
        IF(ITYPLU.NE.2) CALL XABORT('GEOIN1: REAL DATA EXPECTED.(9)')
  230   CONTINUE
  235   CONTINUE
        CALL LCMPUT(IPLIST,'RS',(NSMAX+1)*NG,2,RS)
        DEALLOCATE(RS)
*
        ALLOCATE(MILIE(NMILG))
        DO 240 I=1,NMILG
        CALL REDGET(ITYPLU,MILIE(I),REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
  240   CONTINUE
        CALL LCMPUT(IPLIST,'MILIE',NMILG,1,MILIE)
        DEALLOCATE(MILIE)
*
        ALLOCATE(MIXDL(NMILG))
        DO 250 I=1,NMILG
        CALL REDGET(ITYPLU,MIXDL(I),REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
  250   CONTINUE
        CALL LCMPUT(IPLIST,'MIXDIL',NMILG,1,MIXDL)
        DEALLOCATE(MIXDL)
*
        ALLOCATE(MIXGR(NSMAX*NG*NMILG),FRACT(NG*NMILG))
        DO 260 IOFK=1,NSMAX*NG*NMILG
        MIXGR(IOFK)=0
  260   CONTINUE
        DO 300 I=1,NMILG
        DO 270 J=1,NG
        IOFJ=(I-1)*NG+J
        CALL REDGET(ITYPLU,INTLIR,FRACT(IOFJ),CARLIR,DBLLIR)
        IF(ITYPLU.NE.2) CALL XABORT('GEOIN1: REAL DATA EXPECTED.(10)')
  270   CONTINUE
        DO 290 J=1,NG
        IOFJ=(I-1)*NG+J
        IF(FRACT(IOFJ).GT.0.0) THEN
          DO 280 K=1,NS(J)
          IOFK=((I-1)*NG+(J-1))*NSMAX+K
          CALL REDGET(ITYPLU,MIXGR(IOFK),REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.1) CALL XABORT('GEOIN1: INTEGER DATA EXPECTED.')
  280     CONTINUE
        ENDIF
  290   CONTINUE
  300   CONTINUE
        CALL LCMPUT(IPLIST,'FRACT',NG*NMILG,2,FRACT)
        CALL LCMPUT(IPLIST,'MIXGR',NSMAX*NG*NMILG,1,MIXGR)
        DEALLOCATE(FRACT,MIXGR)
*
        DEALLOCATE(NS)
        CALL XDISET(JSTATE,NSTATE,0)
        JSTATE(1)=NG
        JSTATE(2)=NSMAX+1
        JSTATE(3)=NMILG
        JSTATE(4)=NSMAX*NG
        JSTATE(5)=MICRO
        CALL LCMPUT(IPLIST,'STATE-VECTOR',NSTATE,1,JSTATE)
        CALL LCMSIX(IPLIST,' ',2)
      ELSE IF(CARLIR.EQ.'POURCE') THEN
*       CELL PROPORTIONS FOR DO-IT-YOURSELF OPTION.
        IF(ISTATE(1).NE.30) CALL XABORT('GEOIN1: POURCE - KEY WORD LI'
     1  //'MITED TO DO-IT-YOURSELF GEOMETRY.')
        ALLOCATE(POURC(LX))
        DO 310 I=1,LX
        CALL REDGET(ITYPLU,INTLIR,POURC(I),CARLIR,DBLLIR)
        IF(ITYPLU.NE.2) CALL XABORT('GEOIN1: REAL DATA EXPECTED.(11)')
  310   CONTINUE
        CALL LCMPUT(IPLIST,'POURCE',LX,2,POURC)
        DEALLOCATE(POURC)
      ELSE IF(CARLIR.EQ.'PROCEL') THEN
*       CELL PROBABILITIES FOR DO-IT-YOURSELF OPTION.
        IF(ISTATE(1).NE.30) CALL XABORT('GEOIN1: PROCEL - KEY WORD LI'
     1  //'MITED TO DO-IT-YOURSELF GEOMETRY.')
        ALLOCATE(PROCE(LX*LX))
        DO 325 I=1,LX
        DO 320 J=1,LX
        IOFJ=(J-1)*LX+I
        CALL REDGET(ITYPLU,INTLIR,PROCE(IOFJ),CARLIR,DBLLIR)
        IF(ITYPLU.NE.2) CALL XABORT('GEOIN1: REAL DATA EXPECTED.(12)')
  320   CONTINUE
  325   CONTINUE
        CALL LCMPUT(IPLIST,'PROCEL',LX*LX,2,PROCE)
        DEALLOCATE(PROCE)
      ELSE IF((CARLIR(2:2).EQ.'+').OR.(CARLIR(2:2).EQ.'-').OR.
     1  (CARLIR.EQ.'HBC')) THEN
*       INPUT BOUNDARY CONDITIONS.
        ISURF=0
        IF(CARLIR.EQ.'X-') THEN
          ISURF=1
          IF(LX.EQ.0) CALL XABORT('GEOIN1: HBC X- -> LX=0.')
        ELSE IF(CARLIR.EQ.'X+') THEN
          ISURF=2
          IF(LX.EQ.0) CALL XABORT('GEOIN1: HBC X+ -> LX=0.')
        ELSE IF(CARLIR.EQ.'R+') THEN
          ISURF=2
          IF(ISTATE(1).EQ.10) THEN
            ISURF=4
          ELSE IF(ISTATE(1).EQ.11) THEN
            ISURF=6
          ENDIF
          IF(LR.EQ.0) CALL XABORT('GEOIN1: HBC R+ -> LR=0.')
        ELSE IF(CARLIR.EQ.'Y-') THEN
          ISURF=3
          IF(LY.EQ.0) CALL XABORT('GEOIN1: HBC Y- -> LY=0.')
        ELSE IF(CARLIR.EQ.'Y+') THEN
          ISURF=4
          IF(LY.EQ.0) CALL XABORT('GEOIN1: HBC Y+ -> LY=0.')
        ELSE IF(CARLIR.EQ.'Z-') THEN
          ISURF=5
          IF(LZ.EQ.0) CALL XABORT('GEOIN1: HBC Z- -> LZ=0.')
        ELSE IF(CARLIR.EQ.'Z+') THEN
          ISURF=6
          IF(LZ.EQ.0) CALL XABORT('GEOIN1: HBC Z+ -> LZ=0.')
        ELSE IF(CARLIR.EQ.'HBC') THEN
          ISURF=1
          IF(.NOT.LHEX) CALL XABORT('GEOIN1: HBC PROHIBITED.')
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.3) CALL XABORT('GEOIN1: CHARACTER DATA EXPECTE'
     1    //'D.')
          DO 330 I=1,MAXHEX
          IF(CARLIR.EQ.CHEX(I)) THEN
            IHEX=I
            GO TO 340
          ENDIF
  330     CONTINUE
          CALL XABORT('GEOIN1: INVALID TYPE OF HEXAGONAL SYMMETRY.')
  340     CALL LCMPUT(IPLIST,'IHEX',1,1,IHEX)
          LCOUR=IHEX.EQ.9
        ELSE IF(CARLIR.EQ.'TBC') THEN
          ISURF=1
          IF(.NOT.LTRI) CALL XABORT('GEOIN1: TBC PROHIBITED.')
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.3) CALL XABORT('GEOIN1: CHARACTER DATA EXPECTE'
     1    //'D.')
          DO 350 I=1,MAXTEX
          IF(CARLIR.EQ.CHET(I)) THEN
            ITRI=I
            GO TO 360
          ENDIF
  350     CONTINUE
          CALL XABORT('GEOIN1: INVALID TYPE OF TRIANGULAR SYMMETRY.')
  360     CALL LCMPUT(IPLIST,'ITRI',1,1,ITRI)
        ELSE
          CALL XABORT('GEOIN1: INVALID KEY WORD '//CARLIR//'.')
        ENDIF
        CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT4,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('GEOIN1: CHARACTER DATA EXPECTED.')
        DO 370 I=1,MAXCOD
        IF(TEXT4.EQ.COND(I)) THEN
          NCODE(ISURF)=I
          IF(TEXT4.EQ.'ACYL') NCODE(ISURF)=I-1
          IF(TEXT4.NE.'ALBE') ICODE(ISURF)=0
          GO TO 380
        ENDIF
  370   CONTINUE
        CALL XABORT('GEOIN1: INVALID TYPE OF BOUNDARY CONDITION.')
  380   IF(TEXT4.EQ.'ALBE') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT4,DBLLIR)
          IF(ITYPLU.EQ.1) THEN
            ICODE(ISURF)=INTLIR
            MINICO=MIN(MINICO,INTLIR)
          ELSE IF(ITYPLU.EQ.2) THEN
            ZCODE(ISURF)=REALIR
          ELSE
            CALL XABORT('GEOIN1: INTEGER OR REAL DATA EXPECTED.')
          ENDIF
        ELSE IF(TEXT4.EQ.'ACYL') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT4,DBLLIR)
          IF(ITYPLU.EQ.1) THEN
            ICODE(ISURF)=INTLIR
            MINICO=MIN(MINICO,INTLIR)
          ELSE IF(ITYPLU.EQ.2) THEN
            ZCODE(ISURF)=REALIR
          ELSE
            CALL XABORT('GEOIN1: INTEGER OR REAL DATA EXPECTED '
     1      //'AFTER ACYL.')
          ENDIF
        ELSE IF(TEXT4.EQ.'REFL') THEN
          ZCODE(ISURF)=1.0
        ELSE IF(TEXT4.EQ.'VOID') THEN
          ZCODE(ISURF)=0.0
        ENDIF
      ELSE IF(CARLIR.EQ.';') THEN
*       END-OF-GEOMETRY.
        GO TO 410
      ELSE IF(CARLIR.EQ.':::') THEN
*       INPUT A SUB GEOMETRY.
        IMPX2=IMPX
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('GEOIN1: CHARACTER DATA EXPECTED.')
        CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT4,DBLLIR)
        IF((ITYPLU.NE.3).OR.(TEXT4.NE.':=')) CALL XABORT('GEOIN1: := TO'
     1  //'KEN EXPECTED.')
        CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT4,DBLLIR)
        IF((ITYPLU.NE.3).OR.(TEXT4.NE.'GEO:')) CALL XABORT('GEOIN1: GEO'
     1  //': TOKEN EXPECTED.')
        CALL LCMLEN(IPLIST,CARLIR,ILONG,ITYX)
        IF(ILONG.NE.0) THEN
          IF(ITYX.NE.0) CALL XABORT('GEOIN1: INVALID GEOMETRY NAME.')
        ELSE
          ISTATE(9)=ISTATE(9)+1
        ENDIF
        CALL LCMSIX(IPLIST,CARLIR,1)
        CALL GEOIN1(CARLIR,IPLIST,LEVEL+1,IMPX2,MAXMI2)
        CALL LCMSIX(IPLIST,' ',2)
        ISTATE(7)=MAX(ISTATE(7),MAXMI2)
      ELSE IF(CARLIR.EQ.'MIX-NAMES') THEN
*       DEFINE MIXTURE CHARACTER NAMES.
        IF(LEVEL.NE.1) CALL XABORT('GEOIN1: MIX-NAMES DATA SHOULD BE '
     1  //'WRITTEN ON FIRST DIRECTORY LEVEL.')
        ALLOCATE(CELL(3*LREG))
        I=0
  390   I=I+1
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('GEOIN1: CHARACTER DATA EXPECTED.')
        IF((CARLIR(2:2).EQ.'-').OR.(CARLIR(2:2).EQ.'+').OR.
     1  (CARLIR.EQ.'HBC').OR.(CARLIR(1:4).EQ.'MESH').OR.(CARLIR(1:5)
     2  .EQ.'SPLIT').OR.(CARLIR.EQ.'SIDE').OR.(CARLIR(:3).EQ.'MIX').OR.
     3  (CARLIR.EQ.'CELL').OR.(CARLIR.EQ.'MERGE').OR.(CARLIR.EQ.'TURN')
     4  .OR.(CARLIR(2:4).EQ.'PIN').OR.(CARLIR.EQ.'BIHET').OR.
     5  (CARLIR.EQ.'POURCE').OR.(CARLIR.EQ.'PROCEL').OR.
     6  (CARLIR.EQ.'SECT').OR.(CARLIR.EQ.'RADIUS').OR.
     7  (CARLIR.EQ.'HMIX').OR.(CARLIR.EQ.';').OR. (CARLIR.EQ.':::'))
     8  GO TO 400
        IF(I.GT.LREG) CALL XABORT('GEOIN1: MIX-NAMES INDEX OVERFLOW.')
        READ (CARLIR(:4),'(A4)') CELL(3*(I-1)+1)
        READ (CARLIR(5:8),'(A4)') CELL(3*(I-1)+2)
        READ (CARLIR(9:),'(A4)') CELL(3*(I-1)+3)
        GO TO 390
  400   CALL LCMPUT(IPLIST,'MIX-NAMES',3*(I-1),3,CELL)
        ISTATE(13)=I-1
        DEALLOCATE(CELL)
        GO TO 60
      ELSE
        CALL XABORT('GEOIN1: '//CARLIR//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 50
*
  410 CARLIR='L_GEOM'
      CALL LCMPTC(IPLIST,'SIGNATURE',12,1,CARLIR)
      ISTATE(2)=LR
      ISTATE(3)=LX
      ISTATE(4)=LY
      ISTATE(5)=LZ
      ISTATE(6)=LREG
      CALL LCMPUT(IPLIST,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPLIST,'NCODE',6,1,NCODE)
      CALL LCMPUT(IPLIST,'ZCODE',6,2,ZCODE)
      CALL LCMPUT(IPLIST,'ICODE',6,1,ICODE)
      IF(MINMIX.LT.0)
     >  CALL XABORT('GEOIN1: NEGATIVE MIXTURE NUMBERS INVALID')
      IF(MINICO.LT.1)
     >  CALL XABORT('GEOIN1: ALBEDO NUMBER MUST BE GREATER THAN 0')
      MAXMIX=ISTATE(7)
      IF(IMPX.GT.0) THEN
        CALL LCMINF(IPLIST,CARLIR,TEXT12,EMPTY,ILONG,LCM)
        WRITE (IOUT,510) LEVEL,GEONAM,CARLIR,TYPE(ISTATE(1))
      ENDIF
      IF(IMPX.GT.1) THEN
        WRITE (IOUT,520) ISTATE(1),TYPE(ISTATE(1)),(ISTATE(I),I=2,12)
        WRITE (IOUT,525) (ISTATE(I),I=13,15),ISTATE(18)
      ENDIF
      IF((ISTATE(8).EQ.1).AND.(ISTATE(9).EQ.0)) CALL XABORT('GEOIN1: '
     1 //'CELL OPTION ACTIVATED WITHOUT SUB-GEOMETRIES.')
      DEALLOCATE(ZCODE)
      DEALLOCATE(ICODE,NCODE,JSTATE,ISTATE)
      RETURN
*
  500 CALL XABORT('GEOIN1: INVALID DATA.')
  510 FORMAT(/20H CREATION OF A LEVEL,I3,27H GEOMETRY ON THE DIRECTORY ,
     1 7HNAMED ',A12,21H' OF THE LCM OBJECT ',A12,12H' WITH TYPE ,A16,
     2 1H.)
  520 FORMAT(/14H STATE VECTOR:/
     1 7H ITYPE ,I6, 4H   (,A16,1H)/
     2 7H LR    ,I6,20H   (NUMBER OF TUBES)/
     3 7H LX    ,I6,22H   (X-DIMENSION INDEX)/
     4 7H LY    ,I6,22H   (Y-DIMENSION INDEX)/
     5 7H LZ    ,I6,22H   (Z-DIMENSION INDEX)/
     6 7H LREG  ,I6,22H   (NUMBER OF REGIONS)/
     7 7H MAXMIX,I6,48H   (MAX. NB. OF MIXTURES/0=TRANSPARENT GEOMETRY)/
     8 7H ISUB1 ,I6,34H   (1=COMMAND CELL IS USED/0=ELSE)/
     9 7H ISUB2 ,I6,29H   (NUMBER OF SUB GEOMETRIES)/
     1 7H IMERGE,I6,26H   (1=CELL-MERGING/0=ELSE)/
     2 7H ISPLIT,I6,28H   (1=MESH-SPLITTING/0=ELSE)/
     3 7H IBIHET,I6,34H   (1=DOUBLE HETEROGENEITY/0=ELSE))
  525 FORMAT(
     1 7H ICLUST,I6,28H   (NUMBER OF CLUSTER RINGS)/
     2 7H ISECT ,I6,26H   (TYPE OF SECTORIZATION)/
     3 7H JSECT ,I6,37H   (NUMBER OF NON-SECTORIZED ANNULII)/
     4 7H IPIN  ,I6,24H   (PIN LOCATION OPTION))
  530 FORMAT(' ***** Error in GEOIN1 *****'/
     1       ' Initial number of mixtures ',I10/
     2       ' Cannot be repeated an integer number of times',
     3       ' to fill ',I10,' mixtures')
      END
