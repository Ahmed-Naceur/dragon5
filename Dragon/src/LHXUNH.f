*DECK LHXUNH
      SUBROUTINE LHXUNH(IPTRK,IPGEOM,GEONAM,MESH,NCELA,IPLANZ,NCPHY,
     +              ICODE,ZCODE,MVOSU,NREGIO,ISURF,SIDE,ISTATE,NSMIN,
     +              NSMAX,MVOLUM,IHEX,LX,MCODE,IPLANI,VLAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read and analyse hexagonal geometry.
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): M. Ouisloumen
*
*Parameters: input
* IPTRK   pointer to tracking.
* IPGEOM  pointer to geometry.
* GEONAM  name of geometry
*
*Parameters: output
* ISURF   number of surfaces.
* SIDE    side of hexagone.
* ISTATE  state vector.
* NSMIN   surface minimum index.
* NSMAX   surface maximum index.
* MVOLUM  volume maximum index.
* IHEX    flag for hexagonal geometry (=1 if geometry hexagonal).
* LX      number of hexagones.
* MESH    dimension of array REMESH for mesh storage.
* NCELA   number of cells in assembly after unfolding.
* IPLANZ  number of Z planes.
* NCPHY   number of physical cells.
* ICODE   albedo indices.
* ZCODE   geometric albedos.
* NREGIO  number of physical regions.
* MVOSU   number of volumes ans surfaces.
* MCODE   =1 if NCODE(5)=5 ,=2 if NCODE(6)=5 ,=0 otherwise.
* IPLANI  plane identifier.
* VLAT    lattice indices for surface and volumes. 
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      PARAMETER(NLCM=25,NSTATE=40,NIXS=9,NIST=2)
      TYPE(C_PTR) IPTRK,IPGEOM
      INTEGER  LNLCM(NLCM),INVLCM(NIXS),INVSTA(NIST),ISTATE(NSTATE),
     +         NCODE(6),ICODE(6),JCODE(6),ISTOR(8),VLAT(*)
      REAL     ZCODE(6)
      LOGICAL  L1CELL,LGMERG,LGTURN,LSPLIT,LGCELL,LGSYM,LGSIDE,LGPASS
      LOGICAL  LMERG1,LTURN1
      CHARACTER LCMNM(NLCM)*12,GEONAM*12,TEXT12*12
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IGG,IAA,SURL,JNR1,JCEL,
     + ISECT,LSECT,ISP,ISZ,IBB,ICC,JSS,ISEC1,IXXX,ISS,ISSS,ICELL,IDD,
     + IAD
      REAL, ALLOCATABLE, DIMENSION(:) :: ZZZ,ZMZZ,ZZS,AQ,RR
      TYPE PP
        INTEGER, POINTER, DIMENSION(:) :: ILOCR
      END TYPE PP
      TYPE(PP), ALLOCATABLE, DIMENSION(:) :: JNR2,ISEC2
*----
*  DATA STATEMENTS
*----
      DATA INVLCM,INVSTA /2,3,7,8,12,15,16,17,18,4,13/
      DATA ISTOR /12,6,6,4,3,2,2,2/
      DATA LCMNM /   'MIX',  'MESHX',  'MESHY',  'MESHZ', 'RADIUS',
     >              'SIDE', 'SPLITX', 'SPLITY', 'SPLITZ', 'SPLITR',
     >              'CELL',  'COORD',  'MERGE',   'TURN', 'CLUSTER',
     >              'NPIN',   'RPIN',   'APIN',  'BIHET', 'POURCE',
     >            'PROCEL',   'IHEX',  'NCODE',  'ZCODE', 'ICODE'/
*
      IFCOUR(N)=NINT( (4.+SQRT(1.+4.*FLOAT(N-1)/3.)
     +                 +SQRT(1.+4.*FLOAT(N-2)/3.))*.25)
      IFONC(N,L)= 2+(N-1)*(L+3*(N-2))
*
      LGSYM=.FALSE.
      LGMERG=.FALSE.
      LGTURN=.FALSE.
      L1CELL=.FALSE.
      LSPLIT=.FALSE.
      LGCELL=.FALSE.
      MCODE=0
      NRAY=0
      KRAY=0
      KSECT=0
      IPLANI=0
*
*         LECTURE DES DIFFERENTS BLOCS SUR LCM
*
      DO 10 I=1,NLCM
         CALL LCMLEN(IPGEOM,LCMNM(I),LNLCM(I),ITP)
   10 CONTINUE
*
*        ELEMINATES OPTIONS NOT CHECKED BY THIS ROUTINE
*
      DO 20 I=1,NIXS
         IF(LNLCM(INVLCM(I)).NE.0)
     +           CALL XABORT('LHXUNH : '//GEONAM//' WAS NOT UNFOLDED ')
   20 CONTINUE
*
      CALL LCMLEN(IPGEOM,'STATE-VECTOR',ILEN,ITP)
      IF(ITP.NE.1.OR.ILEN.NE.NSTATE)
     +  CALL XABORT('LHXUNH : INVALID STATE VECTOR ')
      CALL XDRSET(ISTATE,NSTATE,0)
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
*
*         ELEMINATES THE INVALID OPTIONS
*
      DO 30 IST=1, NIST
         IF(ISTATE(INVSTA(IST)).NE.0)
     +        CALL XABORT('LHXUNH : '//GEONAM//'WAS NOT UNFOLDED')
   30 CONTINUE
      ITYPE=ISTATE(1)
      LX=ISTATE(3)
      LZ=ISTATE(5)
      LZZ=LZ
      LREG=ISTATE(6)
      ISTAT9=ISTATE(9)
      IF(ISTATE(8).EQ.1)LGCELL=.TRUE.
      IF(ISTATE(10).EQ.1)LGMERG=.TRUE.
      IF(ISTATE(11).EQ.1)LSPLIT=.TRUE.
      NDIM=2
      NCEL=LX
      IPLANZ=1
      IF(LZ.NE.0) THEN
        NCEL=LX*LZ
        NDIM=3
        NCELP=LX
        IPLANZ=LZ
      ENDIF
      IF(LX.EQ.1)L1CELL=.TRUE.
      IF(LGCELL)THEN
        IF(L1CELL) THEN
        IF(ITYPE.NE.8.AND.ITYPE.NE.9.AND.ITYPE.NE.24.AND.ITYPE.NE.25)
     +     CALL XABORT('LHXUNH : INVALID TYPE OF GEOMETRY')
        ELSE
        IF(ITYPE.NE.8.AND.ITYPE.NE.9)
     +     CALL XABORT('LHXUNH : INVALID TYPE OF GEOMETRY')
        ENDIF
      ELSE
        IF(ITYPE.NE.8.AND.ITYPE.NE.9.AND.ITYPE.NE.24.AND.ITYPE.NE.25)
     +     CALL XABORT('LHXUNH : INVALID TYPE OF GEOMETRY')
      ENDIF
*
*         RECOVERS BOUDARY CONDITIONS
*
      CALL LCMGET(IPGEOM,'NCODE',NCODE)
      CALL LCMGET(IPGEOM,'ZCODE',ZCODE)
      CALL LCMGET(IPGEOM,'ICODE',JCODE)
      DO 35 I=1,6
        IF( JCODE(I).GT.0 )THEN
           ICODE(I)= JCODE(I)
        ELSE
           ICODE(I)= -I
        ENDIF
        IF( NCODE(I).EQ.0 )THEN
           IF( NDIM.EQ.3.AND.(I.EQ.5.OR.I.EQ.6) )THEN
              CALL XABORT('LHXUNH : A BOUNDARY CONDITION IS MISSING')
           ENDIF
           IF( I.EQ.1 )THEN
              CALL XABORT('LHXUNH : A BOUNDARY CONDITION IS MISSING')
           ENDIF
           ICODE(I)= 0
        ENDIF
   35 CONTINUE
*
*       ELEMINATIONS DES OPTIONS IMCOMPATIBLES
*
      DO 40 I=1,6
         IF(NCODE(I).EQ.4.OR.NCODE(I).EQ.7.OR.NCODE(I).EQ.3 .OR.
     +      (NDIM.EQ.2.AND.NCODE(I).EQ.5) .OR.
     +      (NDIM.EQ.2.AND.NCODE(I).EQ.10) )
     +      CALL XABORT('LHXUNH : INVALID BOUNDARY CONDITION ')
         IF((NCODE(5) .EQ. 5 .OR. NCODE(5) .EQ. 10) .AND.
     +      (NCODE(6) .EQ. 5 .OR. NCODE(6) .EQ. 10)) CALL XABORT(
     +      'LHXUNH: UNE SEULE SYMETRIE SELON Z EST VALABLE')
         IF(NCODE(5) .EQ. 5 .OR. NCODE(5) .EQ. 10 ) THEN
           ZCODE(5)=ZCODE(6)
           ICODE(5)=ICODE(6)
         ELSEIF(NCODE(6) .EQ.5 .OR. NCODE(6) .EQ. 10) THEN
           ZCODE(6)=ZCODE(5)
           ICODE(6)=ICODE(5)
         ENDIF
   40 CONTINUE
      IF(L1CELL) THEN
*
*           CAS D'UNE SEULE CELLULE
*
         IF(NCODE(5) .EQ. 5 .OR. NCODE(5) .EQ. 10 .OR.
     >      NCODE(6) .EQ. 5 .OR. NCODE(6) .EQ. 10 ) CALL XABORT(
     +    'LHXUNH: THE B.C. SYME FOR ONE CELL IS NOT PERMITTED ')
        IF(LGCELL) THEN
          CALL LCMLEN(IPGEOM,'CELL',ILEN,ITP)
          IF(ILEN.NE.3*ISTAT9.OR.ITP.NE.3) CALL XABORT('LHXUNH: '
     +                 //'MISSING DIMENSION OR TYPE OF CELL ')
          ALLOCATE(JCEL(3*ISTAT9))
          CALL LCMGET(IPGEOM,'CELL',JCEL)
          WRITE(TEXT12,'(3A4)') (JCEL(ITC),ITC=1,3)
          CALL LCMLEN(IPGEOM,TEXT12,ILEN,ITP)
          IF(ILEN.NE.-1.OR.ITP.NE.0) CALL XABORT('LHXUNH: '
     +                 //'INVALID CELL DATA ')
          CALL LCMSIX(IPGEOM,TEXT12,1)
        ENDIF
        CALL LCMLEN(IPGEOM,'SIDE',ISIDE,ITS)
        IF(ISIDE.EQ.0) CALL XABORT('LHXUNH: SIDE NOT FOUND')
        IF(ITS.NE.2) CALL XABORT('LHXUNH: SIDE MUST BE REAL')
        CALL LCMGET(IPGEOM,'SIDE',SIDE)
        CALL LCMLEN(IPGEOM,'RADIUS',ILENR,ITPR)
        CALL LCMLEN(IPGEOM,'SECTOR',ILENS,ITPS)
        IXX=1
        IF(ILENR.GT.0) IXX=ILENR
        ALLOCATE(ISECT(IXX))
        IF(ILENS.GT.0) THEN
          IF(ILENS.NE.IXX)
     +           CALL XABORT('LHXUNH: INVALID SECTOR DIMENSION(1)')
          CALL LCMGET(IPGEOM,'SECTOR',ISECT)
        ELSE
          CALL XDISET(ISECT,IXX,1)
        ENDIF
        IF(ILENR.GT.0) THEN
          IF(ITPR.NE.2)CALL XABORT('LHXUNH: RADIUS MUST BE REAL')
          CALL LCMLEN(IPGEOM,'SPLITR',ILENSP,ITPS)
          NRAY=ILENR-1
          IF(ILENSP.GT.0) THEN
            IF(ILENSP.NE.NRAY)
     +        CALL XABORT('LHXUNH: INVALID SPLITR DIMENSION')
            IF(ITPS.NE.1) CALL XABORT('LHXUNH: SPLITR MUST BE INTEGER')
            ALLOCATE(ISP(ILENSP))
            CALL LCMGET(IPGEOM,'SPLITR',ISP)
            NRAY=0
            DO 36 J=1,ILENSP
               NRAY=NRAY+ABS(ISP(J))
   36       CONTINUE
            IF(ILENS.GT.0) THEN
               ALLOCATE(LSECT(NRAY+1))
               MADD=-1
               DO 336 J=1,ILENR-1
               KSECT=ISECT(J)
               DO 335 ISPS=1,ISP(J)
                 MADD=MADD+1
                 LSECT(MADD+1)=KSECT
  335          CONTINUE
  336          CONTINUE
               LSECT(NRAY+1)=ISECT(ILENR)
               DEALLOCATE(ISECT)
               ALLOCATE(ISECT(NRAY+1))
               DO 337 J=1,NRAY+1
                 ISECT(J)=LSECT(J)
  337          CONTINUE
               DEALLOCATE(LSECT)
            ENDIF
            DEALLOCATE(ISP)
          ENDIF
        ENDIF
*
        ISURF=6
        MVOSU=0
        DO 39 K=1,NRAY+1
           KSECT=ISECT(K)
           IF(KSECT.GT.1) THEN
             MVOSU=MVOSU+6*(KSECT-1)
           ELSE
             MVOSU=MVOSU+1
           ENDIF
   39   CONTINUE
        DEALLOCATE(ISECT)
        IF(KSECT.GT.1) ISURF=6*(KSECT-1)
        MVOLUM=1+NRAY
        IF(NDIM.EQ.2) THEN
          NCELA=1
          NCPHY=1
          MESH=2+NRAY
          NREGIO=MVOSU
        ELSEIF(NDIM.EQ.3) THEN
          CALL LCMLEN(IPGEOM,'MESHZ',LZ,ITZ)
          IF(LZ.EQ.0) CALL XABORT('LHXUNH: MESHZ NOT FOUND')
          IF(LZ.NE.2) CALL XABORT('LHXUNH: MISSING DIMENSION OF MESHZ')
          IF(ITZ.NE.2)CALL XABORT('LHXUNH: MESHZ MUST BE REAL')
          ALLOCATE(ZZZ(LZ))
          CALL LCMGET(IPGEOM,'MESHZ',ZZZ)
          IF(ZZZ(1).NE.0.) CALL XABORT('LHXUNH: FIRST MESHZ MUST BE 0')
          CALL LCMLEN(IPGEOM,'SPLITZ',LSZ,ITSZ)
          IF(LSZ.GT.0) THEN
            IF(ITSZ.NE.1)CALL XABORT('LHXUNH: SPLITZ MUST BE INTEGER')
            IF(LSZ.NE.LZ-1)CALL XABORT('LHXUNH: WRONG SPLITZ DIMENSION')
            ALLOCATE(ISZ(LSZ))
            CALL LCMGET(IPGEOM,'SPLITZ',ISZ)
            CALL LCMSIX(IPTRK,'DATA_DUP',1)
            CALL LCMPUT(IPTRK,'SPLITZD',LSZ,1,ISZ)
            CALL LCMSIX(IPTRK,' ',2)
            MZ=ISZ(1)+1
            ALLOCATE(ZMZZ(MZ))
            KSP=ISZ(1)
            PAS=(ZZZ(2)-ZZZ(1))/FLOAT(KSP)
            ZMZZ(1)=0.0
            DO 37 K=1,KSP
              ZMZZ(K+1)=PAS*FLOAT(K)+ZZZ(1)
   37       CONTINUE
            DEALLOCATE(ISZ)
          ELSE
            MZ=LZ
            ALLOCATE(ZMZZ(MZ))
            DO 38 J=1,MZ
              ZMZZ(J)=ZZZ(J)
   38       CONTINUE
          ENDIF
          DEALLOCATE(ZZZ)
          CALL LCMSIX(IPTRK,'DATA_DUP',1)
          CALL LCMPUT(IPTRK,'MESHZ',MZ,2,ZMZZ)
          CALL LCMSIX(IPTRK,' ',2)
          DEALLOCATE(ZMZZ)
          IPLANZ=MZ-1
          ISTATE(5)=IPLANZ
          NCELA=MZ-1
          ISURF=NCELA*ISURF
          NCPHY=MZ-1
          MESH=NCELA*(3+NRAY)
          NSMIN=MVOSU
          NSMAX=MVOSU
          MVOSU=NCELA*MVOSU
          MVOLUM=NCELA*MVOLUM
          NREGIO=MVOSU
        ENDIF
        NCOUR=1
        ISURF=ISURF+NSMIN+NSMAX
        RETURN
      ENDIF
*
*        TYPE DE LA SYMETRIE CONSIDEREE
*
      CALL LCMLEN(IPGEOM,'IHEX',ILEN,ITP)
      IF(ITP.NE.1.OR.ILEN.NE.1)CALL XABORT('LHXUNH: TYPE OF '//
     +   'SYMETRIE MUST BE INTEGER')
      CALL LCMGET(IPGEOM,'IHEX',IHEX)
      IF(IHEX.LT.1.OR.IHEX.GT.9) CALL XABORT('LHXUNH: INVALID TYPE OF '
     +   //'SYMETRIE')
      IF(NDIM.EQ.2) THEN
*
*           TRAITEMENT DE LA GEOMETRIE EN 2D
*
        NCELA=NCEL
        NCELAP=NCEL
*
*             DUPLICATION DE LA SYMETRIE CONSIDEREE
*
        IF(IHEX.LT.9) THEN
          ALLOCATE(IBB(ISTOR(IHEX)*NCEL))
          CALL DEPLIT(IHEX,NCEL,NCELA,IBB)
          NCELAP=NCELA
        ELSE
          ALLOCATE(IBB(NCEL))
          DO 137 I=1,NCEL
            IBB(I)=I
  137     CONTINUE
        ENDIF
          CALL LCMLEN(IPGEOM,'MIX',ILENX,ITPX)
          IF(ILENX.NE.NCEL.OR.ITPX.NE.1) THEN
            CALL XABORT('LHXUNH: MISSING TYPE OR DIMENSION OF MIX')
          ENDIF
          NBASE=NCEL+NCELA
          CALL LCMLEN(IPGEOM,'TURN',ILENT,ITP)
          IF(ILENT.GT.0) THEN
            IF(ILENT.NE.NCEL)CALL XABORT('LHXUNH: EXPECTED DIMENSION'//
     +      'OF TURN')
            IF(ITP.NE.1)CALL XABORT('LHXUNH: EXPECTED VALUE IN TURN ')
            LGTURN=.TRUE.
          ENDIF
          NBASE=NBASE*2
          NBASE=NBASE+NCEL+NCELA
          ALLOCATE(ICC(NBASE))
          NBB=2*NCEL
          IF(LGCELL) THEN
            CALL LCMGET(IPGEOM,'MIX',ICC)
            ICMAX=-1
            DO 147 IGX=1,NCEL
              ICC(IGX)=-ICC(IGX)
              ICMAX=MAX(ICMAX,ICC(IGX))
  147       CONTINUE
            IF(ISTAT9.GT.ICMAX)CALL XABORT('LHXUNH: THERE ARE DEFINED'
     +          //' CELLS NOT USED IN THE ASSEMBLY')
          ELSE
            CALL LCMGET(IPGEOM,'MIX',ICC)
          ENDIF
          IF(LGMERG) THEN
            CALL LCMLEN(IPGEOM,'MERGE',IAUXN,ITAUX)
            IF(IAUXN.NE.NCEL.OR.ITAUX.NE.1)CALL XABORT('LHXUNH: DIMEN'
     +         //'SION OR TYPE OF MERGE INVALID')
            CALL LCMGET(IPGEOM,'MERGE',ICC(NCEL+1))
            NCPHY=0
            DO 42 I=1,NCEL
              NCPHY=MAX(NCPHY,ICC(NCEL+I))
   42       CONTINUE
          ELSE
             DO 146 I=1,NCEL
               ICC(NCEL+I)=ICC(I)
  146        CONTINUE
            LGMERG=.TRUE.
            NCPHY=NCEL
          ENDIF
          IF(LGTURN) THEN
            CALL LCMLEN(IPGEOM,'TURN',IAUXN,ITAUX)
            IF(IAUXN.NE.NCEL.OR.ITAUX.NE.1)CALL XABORT('LHXUNH: DIMEN'
     +         //'SION OR TYPE OF TURN INVALID')
            CALL LCMGET(IPGEOM,'TURN',ICC(NBB+1))
          ELSE
            LGTURN=.TRUE.
            DO 777 ITT=1,NCEL
              ICC(NBB+ITT)=1
  777       CONTINUE
          ENDIF
          NBB=NBB+NCEL
*
*            DUPLICATION DE MERGE ET MIX
*
           DO 47 I=1,NCELA
            ICC(NBB+I)=ICC(IBB(I))
            ICC(NBB+NCELA+I)=ICC(IBB(I)+NCEL)
   47      CONTINUE
*
*            CAS OU SIDE EST ENTREE UNE SEULE FOIS
*
          LGSIDE=.TRUE.
          CALL LCMLEN(IPGEOM,'SIDE',ISIDE,IT)
          IF(ISIDE.GT.0) THEN
            CALL LCMGET(IPGEOM,'SIDE',SIDE)
            LGSIDE=.FALSE.
          ELSEIF(.NOT.LGCELL) THEN
            CALL XABORT('LHXUNH: SIDE NOT FOUND ')
          ENDIF
*
*            EVALUATION DE MESH
*
          MESH=2*NCELA
          NCOUR=1
          IF(NCELAP.GT.1)NCOUR=IFCOUR(NCELAP)
          ICELC0= IFONC(NCOUR,0)
          ICELC1= IFONC(NCOUR,1)-ICELC0
          IF(LGCELL) THEN
            ALLOCATE(JNR1(ISTAT9),JNR2(ISTAT9),JCEL(3*ISTAT9))
            CALL LCMGET(IPGEOM,'CELL',JCEL)
            NREGIO=0
            DO 48 I=1,ISTAT9
              IRTC=3*I-2
              WRITE(TEXT12,'(3A4)') (JCEL(ITC),ITC=IRTC,IRTC+2)
               CALL LCMSIX(IPGEOM,TEXT12,1)
               CALL LCMLEN(IPGEOM,'RADIUS',ILENR,ITP)
               CALL LCMLEN(IPGEOM,'SPLITR',ILENSP,ITP)
               IF(LGSIDE) THEN
                 CALL LCMLEN(IPGEOM,'SIDE',ISIDE,IT)
                 IF(ISIDE.EQ.0)CALL XABORT('LHXUNH: SIDE NOT FOUND')
                 CALL LCMGET(IPGEOM,'SIDE',SSIDE)
                 IF(I.GT.1) THEN
               IF(SSIDE.NE.SIDE)CALL XABORT('LHXUNH: INCOMPATIBLE SIDE')
                 ENDIF
                 SIDE=SSIDE
               ENDIF
               NRAY=ILENR
               IF(ILENR.GT.0)NRAY=ILENR-1
               IF(ILENSP.GT.0) THEN
                 ALLOCATE(JSS(ILENSP))
                 CALL LCMGET(IPGEOM,'SPLITR',JSS)
                 NRAY=0
                 DO 46 J=1,ILENSP
                   NRAY=NRAY+ABS(JSS(J))
   46            CONTINUE
                 DEALLOCATE(JSS)
               ENDIF
               JNR1(I)=NRAY
               ALLOCATE(JNR2(I)%ILOCR(NRAY+1))
               CALL LCMLEN(IPGEOM,'SECTOR',ILENS,ITPS)
               IF(ILENS.GT.0) THEN
                 IF(ILENS.NE.NRAY+1)
     +                CALL XABORT('LHXUNH: INVALID SECTOR DIMENSION(2)')
                 CALL LCMGET(IPGEOM,'SECTOR',JNR2(I)%ILOCR)
               ELSE
                 CALL XDISET(JNR2(I)%ILOCR,NRAY+1,1)
                 NREGIO=NREGIO+1+NRAY
                 GOTO 430
               ENDIF
               DO 43 K=1,NRAY+1
                 KSECT=JNR2(I)%ILOCR(K)
                 IF(KSECT.GT.1) THEN
                   NREGIO=NREGIO+6*(KSECT-1)
                 ELSE
                   NREGIO=NREGIO+1
                 ENDIF
   43          CONTINUE
  430          CONTINUE
               CALL LCMSIX(IPGEOM,' ',2)
   48       CONTINUE
            DEALLOCATE(JCEL)
            NRAY=0
            MVOSU=0
            MVOLUM=0
            ISURF=0
            DO 49 I=1,NCELA
               KRAY=JNR1(ICC(NBB+I))
               NRAY=NRAY+KRAY
               DO 490 JR=1,KRAY+1
                 KSECT=JNR2(ICC(NBB+I))%ILOCR(JR)
                 IF(KSECT.GT.1) THEN
                   MVOSU=MVOSU+6*(KSECT-1)
                 ELSE
                   MVOSU=MVOSU+1
                 ENDIF
  490          CONTINUE
               IF(I.GE.ICELC0) THEN
                 IF(KSECT.GT.1) THEN
                   ISURF=ISURF+2*(KSECT-1)
                   IF(MOD(I-ICELC0,ICELC1).EQ.0) ISURF=ISURF+KSECT-1
                 ELSE
                   ISURF=ISURF+2
                   IF(MOD(I-ICELC0,ICELC1).EQ.0) ISURF=ISURF+1
                 ENDIF
               ENDIF
               MVOLUM=MVOLUM+KRAY+1
   49       CONTINUE
            MESH=NRAY+MESH
            DO 495 I=1,ISTAT9
               DEALLOCATE(JNR2(I)%ILOCR)
  495       CONTINUE
            DEALLOCATE(JNR2,JNR1)
          ELSE
            CALL LCMLEN(IPGEOM,'SECTOR',ILENS,ITP)
            NREGIO=NCPHY
            MVOSU=NCELA
            MVOLUM=NCELA
            ALLOCATE(IXXX(NCEL))
            IF(ILENS.GT.0) THEN
              IF(ILENS.NE.NCEL)
     +               CALL XABORT('LHXUNH: INVALID SECTOR DIMENSION(3)')
              CALL LCMGET(IPGEOM,'SECTOR',IXXX)
            ELSE
              CALL XDISET(IXXX,NCEL,1)
            ENDIF
            ILESS=0
            NREGIO=0
            DO 50 L=1,ILENS
               MCEL=ICC(L)
               IF(MCEL.GT.ILESS) THEN
                 ILESS=MCEL
                 KSECT=IXXX(L)
                 IF(KSECT.GT.1) THEN
                   NREGIO=NREGIO+6*(KSECT-1)
                 ELSE
                   NREGIO=NREGIO+1
                 ENDIF
               ENDIF
   50       CONTINUE
            MVOSU=0
            ISURF=0
            DO 51 I=1,NCELA
               KSECT=IXXX(ICC(NBB+I))
               IF(KSECT.GT.1) THEN
                 MVOSU=MVOSU+6*(KSECT-1)
                 IF(I.GE.ICELC0) THEN
                   ISURF=ISURF+2*(KSECT-1)
                   IF(MOD(I-ICELC0,ICELC1).EQ.0) ISURF=ISURF+KSECT-1
                 ENDIF
               ELSE
                 MVOSU=MVOSU+1
                 IF(I.GE.ICELC0) THEN
                   ISURF=ISURF+2
                   IF(MOD(I-ICELC0,ICELC1).EQ.0) ISURF=ISURF+1
                 ENDIF
               ENDIF
   51       CONTINUE
            DEALLOCATE(IXXX)
          ENDIF
*
*             DUPLICATION DE TURN
*
          MBB=2*NCELA
          IF(LGTURN) THEN
           IF (IHEX.LT.9) THEN
            CALL DUTURN(IHEX,ICC(NBB-NCEL+1),NCEL,ICC(NBB+MBB+1),
     +                  NCELA,IBB)
           ELSE
            DO 187 I=1,NCELA
               ICC(NBB+MBB+I)=ICC(NBB-NCEL+I)
  187       CONTINUE
           ENDIF
          ENDIF
*
*             RESTORAGE DES DONNEES
*
          CALL LCMSIX(IPTRK,'DATA_DUP',1)
*
*--- CE VECTEUR EST UTILE DANS LA ROUTINE MESHST
          CALL LCMPUT(IPTRK,'GENER0',NCELA,1,IBB)
          IF(LGCELL) THEN
            CALL LCMPUT(IPTRK,'GENERATING',NCELA,1,ICC(NBB+1))
          ELSE
            CALL LCMPUT(IPTRK,'MIX',NCELA,1,ICC(NBB+1))
          ENDIF
          IF(LGMERG)
     +   CALL LCMPUT(IPTRK,'MERGE',NCELA,1,ICC(NCELA+NBB+1))
         IF(LGTURN)
     +   CALL LCMPUT(IPTRK,'TURN',NCELA,1,ICC(NBB+MBB+1))
          CALL LCMSIX(IPTRK,' ',2)
          DEALLOCATE(IBB,ICC)
*
      ELSEIF(NDIM.EQ.3) THEN
*
*         TRAITEMENT DU CAS 3D
*
         NCELA=NCEL
         CALL LCMLEN(IPGEOM,'MIX',ILENX,ITPX)
         IF(ILENX.NE.NCEL.OR.ITPX.NE.1) CALL XABORT('LHXUNH: MISSING'
     +                                              //'VECTOR MIX')
         NBASE=NCEL
         CALL LCMLEN(IPGEOM,'TURN',ILENT,ITPT)
         IF(ILENT.GT.0)THEN
           IF(ILENT.NE.NCEL.OR.ITPT.NE.1) CALL XABORT('LHXUNH: MISSING'
     +                                                //'VECTOR TURN')
           LGTURN=.TRUE.
         ENDIF
         NBASE=NCEL+NCEL
         NBASE=NCEL+NBASE
         ALLOCATE(ICC(NBASE))
         CALL LCMGET(IPGEOM,'MIX',ICC)
         IF(LGCELL) THEN
           DO 144 IGX=1,NCEL
             ICC(IGX)=-ICC(IGX)
  144      CONTINUE
         ENDIF
         NBB=NCEL+NCEL
         LMERG1=.FALSE.
         LTURN1=.FALSE.
         IF(LGMERG) THEN
           CALL LCMGET(IPGEOM,'MERGE',ICC(NCEL+1))
           NCPHY=0
           DO 52 I=1,NCEL
              NCPHY=MAX(NCPHY,ICC(NCEL+I))
   52      CONTINUE
           LMERG1=.TRUE.
         ELSE
           DO 53 I=1,NCEL
              ICC(I+NCEL)=ICC(I)
   53      CONTINUE
           LGMERG=.TRUE.
           NCPHY=NCEL
         ENDIF
         IF(LGTURN) THEN
           CALL LCMGET(IPGEOM,'TURN',ICC(NBB+1))
           LTURN1=.TRUE.
         ELSE
           LGTURN=.TRUE.
           DO 778 ITT=1,NCEL
             ICC(NBB+ITT)=1
  778      CONTINUE
         ENDIF
         NBB=NBB+NCEL
*
*           SPLITING DE LA DIRECTION Z
*
         CALL LCMLEN(IPGEOM,'MESHZ',ILENZ,ITPZ)
         IF(LGCELL.OR.IPLANZ.GT.1) THEN
         IF(ILENZ.NE.0)CALL XABORT('LHXUNH: INVALID POSITION OF MESHZ')
         ENDIF
         IF(ILENZ.GT.0) THEN
           IF(ILENZ.NE.LZ+1.OR.ITPZ.NE.2)CALL XABORT('LHXUNH: MISSING'
     +                                 //'DIMENSION OR VALUE IN MESHZ')
           ALLOCATE(ZZZ(ILENZ))
           CALL LCMGET(IPGEOM,'MESHZ',ZZZ)
           MZ=ILENZ
         ENDIF
*
*            EVALUATION DE NREGIO
*
         NREGIO=NCPHY
         CALL LCMLEN(IPGEOM,'SECTOR',ILENS,ITP)
         IF(LGCELL.OR.IPLANZ.GT.1) THEN
         IF(ILENS.NE.0)CALL XABORT('LHXUNH: INVALID POSITION OF SECTOR')
         ENDIF
         IF(ILENS.GT.0) THEN
           IF(ILENS.NE.NCEL)CALL XABORT('LHXUNH: INVALID DIMENSION OF'
     +                                  //'SECTOR')
           IF(ITP.NE.1)CALL XABORT('LHXUNH: SECTOR MUST BE INTEGER')
           ALLOCATE(ISEC1(ILENS))
           ILESS=0
           NREGIO=0
           CALL LCMGET(IPGEOM,'SECTOR',ISEC1)
           DO 54 L=1,ILENS
              MCEL=ICC(L)
              IF(MCEL.GT.ILESS) THEN
                ILESS=MCEL
                KSECT=ISEC1(L)
                IF(KSECT.GT.1) THEN
                   NREGIO=NREGIO+6*(KSECT-1)
                ELSE
                   NREGIO=NREGIO+1
                ENDIF
              ENDIF
   54      CONTINUE
         ENDIF
*
         CALL LCMLEN(IPGEOM,'SPLITZ',ILEN,ITP)
         IF(ILEN.NE.0) THEN
           IF(ILEN.NE.LZ)CALL XABORT('LHXUNH: INVALID DIMENSION OF'
     +                               //'SPLITZ')
           IF(ILENZ.LT.0)CALL XABORT('LHXUNH: MESHZ MUST BE DEFINED'
     +                               //'LIKE SPLITZ')
           ALLOCATE(ISS(LZ))
           CALL LCMGET(IPGEOM,'SPLITZ',ISS)
           CALL LCMSIX(IPTRK,'DATA_DUP',1)
           CALL LCMPUT(IPTRK,'SPLITZD',LZ,1,ISS)
           CALL LCMSIX(IPTRK,' ',2)
           MZ=0
           DO 55 K=1,LZ
              MZ=MZ+ISS(K)
   55      CONTINUE
           LSPLIT=.TRUE.
         ELSEIF(LGCELL) THEN
           CALL LCMLEN(IPGEOM,'CELL',ILEN,ITP)
           IF(ILEN.NE.3*ISTAT9.OR.ITP.NE.3)CALL XABORT('LHXUNH: '
     +                 //'MISSING DIMENSION OR TYPE OF CELL ')
           ALLOCATE(ICELL(3*ISTAT9),RR(3*ISTAT9),ISSS(ISTAT9))
           CALL LCMGET(IPGEOM,'CELL',ICELL)
           ALLOCATE(ISEC1(ISTAT9),ISEC2(ISTAT9))
           IAUX=0
           NREGIO=0
           DO 70 I=1,ISTAT9
              WRITE(TEXT12(1:4),'(A4)')ICELL(3*I-2)
              WRITE(TEXT12(5:8),'(A4)')ICELL(3*I-1)
              WRITE(TEXT12(9:12),'(A4)')ICELL(3*I)
              CALL LCMSIX(IPGEOM,TEXT12,1)
              CALL LCMGET(IPGEOM,'MESHZ',RR(IAUX+1))
              CALL LCMGET(IPGEOM,'SIDE',RR(2*ISTAT9+I))
              IAUX=IAUX+2
              CALL LCMLEN(IPGEOM,'SPLITZ',ILEN,ITP)
              ISSS(I)=1
              IF(ILEN.NE.0) THEN
               IF(ITP.NE.1)CALL XABORT('LHXUNH: INVALID TYPE OF SPLITZ')
                CALL LCMGET(IPGEOM,'SPLITZ',ISSS(I))
                LSPLIT=.TRUE.
              ENDIF
              CALL LCMLEN(IPGEOM,'RADIUS',ILENR,ITP)
              CALL LCMLEN(IPGEOM,'SPLITR',ILENSP,LT)
              NRAY=0
              NSECT=1
              IF(ILENR.GT.0) THEN
                IF(ITP.NE.2)CALL XABORT('LHXUNH: RADIUS MUST BE REAL ')
                NRAY=ILENR-1
                NSECT=ILENR
              ENDIF
              CALL LCMLEN(IPGEOM,'SECTOR',ILENS,ITP)
              ALLOCATE(ISECT(NSECT))
              IF(ILENS.GT.0) THEN
                 IF(ILENS.NE.NSECT)
     +                CALL XABORT('LHXUNH: INVALID SECTOR DIMENSION(4)')
                 CALL LCMGET(IPGEOM,'SECTOR',ISECT)
              ELSE
                 CALL XDISET(ISECT,NSECT,1)
              ENDIF
              IF(ILENSP.GT.0) THEN
                IF(LT.NE.1)CALLXABORT('LHXUNH: SPLITR MUST BE INTEGER')
                ALLOCATE(JSS(ILENSP))
                CALL LCMGET(IPGEOM,'SPLITR',JSS)
                NRAY=0
                DO 246 JP=1,ILENSP
                  NRAY=NRAY+ABS(JSS(JP))
  246           CONTINUE
                ALLOCATE(ISEC2(I)%ILOCR(NRAY+1))
                MAD=-1
                DO 248 JP=1,ILENS
                   KSECT=ISECT(JP)
                   DO 247 JT=1,JSS(JP)
                     MAD=MAD+1
                     ISEC2(I)%ILOCR(MAD+1)=KSECT
  247              CONTINUE
  248           CONTINUE
                ISEC2(I)%ILOCR(MAD+2)=ISECT(ILENS+1)
                DEALLOCATE(JSS)
              ELSE
                ALLOCATE(ISEC2(I)%ILOCR(NRAY+1))
                DO 249 JP=1,NSECT
                   ISEC2(I)%ILOCR(JP)=ISECT(JP)
  249           CONTINUE
              ENDIF
              DEALLOCATE(ISECT)
              ISEC1(I)=NRAY
              DO 255 K=1,NRAY+1
                 KSECT=ISEC2(I)%ILOCR(K)
                 IF(KSECT.GT.1) THEN
                   NREGIO=NREGIO+6*(KSECT-1)
                 ELSE
                   NREGIO=NREGIO+1
                 ENDIF
  255         CONTINUE
              CALL LCMSIX(IPGEOM,' ',2)
   70      CONTINUE
           DEALLOCATE(ICELL)
*
*                 COMPATIBILITE DES DONNEES DES CELLULES D'UN MEME PLAN
*
           LZZ=LZ
           IAUX=0
           LISP=0
           ALLOCATE(ISS(ISTAT9),ZZZ(LZ+1))
           MZZ=0
           DO 80 I=1,LZ
              SIDE=RR(2*ISTAT9+ICC(IAUX+1))
              LCC=2*(ICC(IAUX+1)-1)
              Z1=RR(LCC+1)
              Z2=RR(LCC+2)
              ZZZ(I)=Z1
              ZZZ(I+1)=Z2
              IF(LSPLIT) THEN
                 LISP=ISSS(ICC(IAUX+1))
                 ISS(I)=LISP
                 MZZ=MZZ+LISP
                 LZZ=MZZ
              ENDIF
              DO 75 K=2,NCELP
                 IAUX=IAUX+1
                 IF(SIDE.NE.RR(2*ISTAT9+ICC(IAUX+1)))
     +                       CALL XABORT('LHXUNH: INCOMPATIBLE SIDE ')
                 LCC=2*(ICC(IAUX+1)-1)
                 IF((Z1.NE.RR(LCC+1)).OR.(Z2.NE.RR(LCC+2)))
     +                    CALL XABORT('LHXUNH: INCOMPATIBLE MESHZ ')
                 IF(LSPLIT) THEN
                   IF(LISP.NE.ISSS(ICC(IAUX+1)))
     +                       CALL XABORT('LHXUNH: INCOMPATIBLE SPLITZ')
                 ENDIF
   75         CONTINUE
              IAUX=IAUX+1
   80      CONTINUE
           DEALLOCATE(RR,ISSS)
           MZ=LZ+1
           IF(LSPLIT)MZ=MZZ+1
         ENDIF
         ALLOCATE(ZZS(MZ))
         IF(LSPLIT) THEN
           LAUX=0
           Z1=ZZZ(1)
           ZZS(1)=Z1
           DO 83 K=1,LZ
              MSP=ISS(K)
              Z2=ZZZ(K+1)
              PAS=(Z2-Z1)/REAL(MSP)
              DO 81 L=1,MSP-1
                 LAUX=LAUX+1
                 ZZS(LAUX+1)=Z1+PAS*REAL(L)
   81         CONTINUE
              LAUX=LAUX+1
              ZZS(LAUX+1)=Z2
              Z1=Z2
   83      CONTINUE
         ELSE
           DO 88 J=1,MZ
             ZZS(J)=ZZZ(J)
   88      CONTINUE
         ENDIF
         DEALLOCATE(ZZZ)
*
*               STORAGE DES VALEURS DE Z
*
         CALL LCMSIX(IPTRK,'DATA_DUP',1)
         LSPLZM=1
         LSPLZP=1
         IF(LSPLIT) THEN
           LSPLZM=ISS(1)
           LSPLZP=ISS(LZ)
           CALL LCMPUT(IPTRK,'SPLITZD',LZ,1,ISS)
         ENDIF
         IF(NCODE(5) .EQ. 5 .OR. NCODE(5) .EQ. 10) THEN
           IF(NCODE(5) .EQ. 5) THEN 
             MDZ=2*MZ-1-LSPLZM
             IOFZD=2+LSPLZM
           ELSE
             MDZ=2*MZ-1
             IOFZD=2
           ENDIF
           ALLOCATE(AQ(MDZ))
           IQ=0
           DO 583 IW=MZ,IOFZD,-1
              AQ(IQ+1)=-ZZS(IW)+ZZS(1)+ZZS(2)
              IQ=IQ+1
  583      CONTINUE
           DO 584 IW=1,MZ
              AQ(IQ+1)=ZZS(IW)
              IQ=IQ+1
  584      CONTINUE
           CALL LCMPUT(IPTRK,'MESHZ',MDZ,2,AQ)
           DEALLOCATE(AQ)
         ELSEIF(NCODE(6) .EQ. 5 .OR. NCODE(6) .EQ. 10) THEN
           IF(NCODE(6) .EQ. 5) THEN 
             MDZ=2*MZ-LSPLZP-1
             IOFZD=LSPLZP+1
           ELSE
             MDZ=2*MZ-1
             IOFZD=1
           ENDIF
           ALLOCATE(AQ(MDZ))
           IQ=0
           DO 585 IW=1,MZ
              AQ(IQ+1)=ZZS(IW)
              IQ=IQ+1
  585      CONTINUE
           DO 586 IW=MZ-IOFZD,1,-1
             AQ(IQ+1)=AQ(IQ)+ZZS(IW+1)-ZZS(IW)
             IQ=IQ+1
  586      CONTINUE
           CALL LCMPUT(IPTRK,'MESHZ',MDZ,2,AQ)
           DEALLOCATE(AQ)
         ELSE
           CALL LCMPUT(IPTRK,'MESHZ',MZ,2,ZZS)
         ENDIF
         CALL LCMSIX(IPTRK,' ',2)
         DEALLOCATE(ZZS)
*
         IF(IHEX.LT.9) THEN
*
*             DUPLICATION DE LA GEOMETRIE CONSIDEREE
*
           ALLOCATE(IBB(ISTOR(IHEX)*NCEL))
           CALL DEPLIT(IHEX,NCELP,NCELAP,IBB)
           IAUX=NCELP
           JAUX=NCELAP
           DO 65 K=2,LZ
              DO 60 L=1,NCELAP
                 IBB(JAUX+L)=IBB(L)+IAUX
   60         CONTINUE
              IAUX=IAUX+NCELP
              JAUX=JAUX+NCELAP
   65      CONTINUE
         ELSEIF(IHEX.EQ.9) THEN
           NCELAP=NCELA/LZ
         ENDIF
         IF(LSPLIT)THEN
            NCELA=NCELAP*LZZ
         ELSE
            NCELA=NCELAP*LZ
         ENDIF
*
*        DUPLICATION DE MERGE ET MIX
*
         NBASE=NCELA+NCELA
         IF(LGTURN)NBASE=NBASE+NCELA
         ALLOCATE(IDD(NBASE),IAD(NCELA))
         IAUX=-1
         JAUX=0
         ICELSP=0
         MAXSP=0
         MINSP=99999
         DO 85 I=1,LZ
           IF(IHEX.LT.9) THEN
              MAXSP=IBB(JAUX+1)
              MINSP=MAXSP
           ENDIF
           DO 84 J=1,NCELAP
             IAUX=IAUX+1
             IF(IHEX.LT.9) THEN
                IAD(IAUX+1)=IBB(JAUX+J)+ICELSP
                MAXSP=MAX(MAXSP,IBB(JAUX+J))
                MINSP=MIN(MINSP,IBB(JAUX+J))
                IDD(IAUX+1)=ICC(IBB(JAUX+J))
                IDD(NCELA+IAUX+1)=ICC(NCEL+IBB(JAUX+J))
             ELSE
                IAD(IAUX+1)=IAUX+1
                IDD(IAUX+1)=ICC(JAUX+J)
                IDD(NCELA+IAUX+1)=ICC(NCEL+JAUX+J)
             ENDIF
   84      CONTINUE
           IF(LSPLIT) THEN
             LISP=ISS(I)
             DO 86 K=2,LISP
               IF(IHEX.LT.9) ICELSP=ICELSP+MAXSP-MINSP+1
               DO 82 J=1,NCELAP
                 IAUX=IAUX+1
                 IF(IHEX.LT.9) THEN
                    IAD(IAUX+1)=IBB(JAUX+J)+ICELSP
                    IDD(IAUX+1)=ICC(IBB(JAUX+J))
                    IDD(NCELA+IAUX+1)=ISS(NCEL+IBB(JAUX+J))
                 ELSE
                    IAD(IAUX+1)=IAUX+1
                    IDD(IAUX+1)=ICC(JAUX+J)
                    IDD(NCELA+IAUX+1)=ICC(NCEL+JAUX+J)
                 ENDIF
   82          CONTINUE
   86        CONTINUE
           ENDIF
           JAUX=JAUX+NCELAP
   85    CONTINUE
*
*          DUPLICATION DE TURN
*
         MBB=NCELA+NCELA
         IF(LGTURN) THEN
           IAUX=2*NCEL
           JAUX=0
           DO 95 I=1,LZ
           IF(IHEX.LT.9) THEN
             CALL DUTURN(IHEX,ICC(IAUX+1),NCELP,IDD(MBB+JAUX+1),
     +              NCELAP,IBB)
           ELSE
             DO 87 IV=1,NCELAP
               IDD(MBB+JAUX+IV)=ICC(IAUX+IV)
   87        CONTINUE
           ENDIF
             IAUX=IAUX+NCELP
             KAUX=JAUX-1
             JAUX=JAUX+NCELAP
             IF(LSPLIT) THEN
               LISP=ISS(I)
               DO 92 J=2,LISP
                  DO 90 K=1,NCELAP
                    IDD(MBB+JAUX+1)=IDD(MBB+KAUX+K+1)
                    JAUX=JAUX+1
   90             CONTINUE
   92          CONTINUE
             ENDIF
   95      CONTINUE
           DEALLOCATE(ISS)
         ENDIF
         IF(IHEX.NE.9) DEALLOCATE(IBB)
         DEALLOCATE(ICC)
         IF(NCODE(5) .EQ. 5 .OR. NCODE(5) .EQ. 10) THEN
*
*             DUPLICATION DE LA SYMETRIE SELON L'AXE Z-
*
           IF(NCODE(6) .EQ. 5 .OR. NCODE(6) .EQ. 10) CALL XABORT(
     >     'LHXUNH: ONLY ONE Z SYMETRY IS PERMITED')
           IF(NCODE(6) .EQ. 5) THEN 
             KBB=NCELA+NCELA-NCELAP*LSPLZM
             IOFZD=LSPLZM
           ELSE
             KBB=NCELA+NCELA
             IOFZD=0
           ENDIF
           MCODE=1
           KK=KBB
           KBB=2*KBB
           IF(LGTURN)KBB=KBB+KK
           ALLOCATE(IGG(KBB),IAA(KK))
           IAUX=0
           JAUX=NCELA-NCELAP
           LAUX=2*NCELA-NCELAP
*
*           DUPLICATION DE MIX ET MERGE
*
           NCOUR=IFCOUR(NCELAP)
           DO 110 I=IOFZD,LZZ-1
             IPOS=0
             IAA(IAUX+1)=IAD(JAUX+1)
             IGG(IAUX+1)=IDD(JAUX+1)
             IGG(KK+IAUX+1)=IDD(LAUX+1)
             DO 100 J1=2,NCOUR
             DO 99 J=IFONC(J1,5)-1,IFONC(J1,0)-1,-1
               IPOS=IPOS+1
               IAA(IAUX+IPOS+1)=IAD(JAUX+J+1)
               IGG(IAUX+IPOS+1)=IDD(JAUX+J+1)
               IGG(KK+IAUX+IPOS+1)=IDD(LAUX+J+1)
   99        CONTINUE
             NCLCOR=1+3*J1*(J1-1)
             DO 105 J=NCLCOR-1,IFONC(J1,5),-1
               IPOS=IPOS+1
               IAA(IAUX+IPOS+1)=IAD(JAUX+J+1)
               IGG(IAUX+IPOS+1)=IDD(JAUX+J+1)
               IGG(KK+IAUX+IPOS+1)=IDD(LAUX+J+1)
  105        CONTINUE
  100        CONTINUE
             IAUX=IAUX+NCELAP
             JAUX=JAUX-NCELAP
             LAUX=LAUX-NCELAP
  110      CONTINUE
           DO 115 K=1,NCELA
             IAA(IAUX+K)=IAD(K)
             IGG(IAUX+K)=IDD(K)
             IGG(KK+IAUX+K)=IDD(NCELA+K)
  115      CONTINUE
*
*             DUPLICATION DE TURN
*
           IF(LGTURN) THEN
             IAUX=0
             LAUX=3*NCELA-NCELAP*IOFZD
             DO 130 I=IOFZD,LZZ-1
               J=-1
               JAUX2=0
               JAUX1=0
               DO 120 J1=1,NCOUR
               LGPASS=.TRUE.
  118          CONTINUE
               DO 119 J2=JAUX2,JAUX1,-1
                 J=J+1
                 ITURN=IDD(LAUX+J2+1)
                 IF(ITURN.LE.6) THEN
                   IF(ITURN.EQ.1) THEN
                     IGG(2*KK+IAUX+J+1)=12
                   ELSEIF(ITURN.EQ.2) THEN
                     IGG(2*KK+IAUX+J+1)=11
                   ELSEIF(ITURN.EQ.3) THEN
                     IGG(2*KK+IAUX+J+1)=10
                   ELSEIF(ITURN.EQ.4) THEN
                     IGG(2*KK+IAUX+J+1)=9
                   ELSEIF(ITURN.EQ.5) THEN
                     IGG(2*KK+IAUX+J+1)=8
                   ELSEIF(ITURN.EQ.6) THEN
                     IGG(2*KK+IAUX+J+1)=7
                   ENDIF
                 ELSEIF(ITURN.LE.12) THEN
                   IF(ITURN.EQ.12) THEN
                     IGG(2*KK+IAUX+J+1)=1
                   ELSEIF(ITURN.EQ.11) THEN
                     IGG(2*KK+IAUX+J+1)=2
                   ELSEIF(ITURN.EQ.10) THEN
                     IGG(2*KK+IAUX+J+1)=3
                   ELSEIF(ITURN.EQ.9) THEN
                     IGG(2*KK+IAUX+J+1)=4
                   ELSEIF(ITURN.EQ.8) THEN
                     IGG(2*KK+IAUX+J+1)=5
                   ELSEIF(ITURN.EQ.7) THEN
                     IGG(2*KK+IAUX+J+1)=6
                   ENDIF
                 ELSE
                   CALL XABORT('LHXUNH: EXPECTED VALUE OF TURN ')
                 ENDIF
  119          CONTINUE
               IF(LGPASS) THEN
                 JAUX2=1+3*J1*(J1-1)-1
                 JAUX1=IFONC(J1,5)
                 LGPASS=.FALSE.
                 IF(J1.GT.1) GOTO 118
               ENDIF
               JAUX2=IFONC(J1+1,5)-1
               JAUX1=IFONC(J1+1,0)-1
  120          CONTINUE
               IAUX=IAUX+NCELAP
               LAUX=LAUX-NCELAP
  130        CONTINUE
           ENDIF
           DO 135 K=1,NCELA
             IGG(2*KK+IAUX+K)=IDD(2*NCELA+K)
  135      CONTINUE
           IPLANI=LZZ
           LZZ=2*LZZ-IOFZD
         ELSEIF(NCODE(6) .EQ. 5 .OR. NCODE(6) .EQ. 10) THEN
*
*           DULPLICATION DE LA GEOMETRIE SELON Z+
*
           IF(NCODE(6) .EQ. 5) THEN 
             KBB=NCELA+NCELA-NCELAP*LSPLZP
             IOFZD=LSPLZP
           ELSE
             KBB=NCELA+NCELA
             IOFZD=0
           ENDIF
           MCODE=2
           KK=KBB
           KBB=2*KBB
           IF(LGTURN)KBB=KBB+KK
           ALLOCATE(IGG(KBB),IAA(KK))
*
*           DUPLICATION DE MIX ET MERGE
*
           DO 140 I=1,NCELA
             IAA(I)=IAD(I)
             IGG(I)=IDD(I)
             IGG(KK+I)=IDD(NCELA+I)
  140      CONTINUE
           NCOUR=IFCOUR(NCELAP)
           LFIN=-NCELAP*IOFZD+NCELA
           IPOS=-1
           DO 145 I=IOFZD,LZZ-1
             LFIN=LFIN-NCELAP
             IPOS=IPOS+1
             IAA(NCELA+IPOS+1)=IAD(LFIN+1)
             IGG(NCELA+IPOS+1)=IDD(LFIN+1)
             IGG(KK+NCELA+IPOS+1)=IDD(NCELA+LFIN+1)
             DO 143 J1=2,NCOUR
             DO 141 J=IFONC(J1,5)-1,IFONC(J1,0)-1,-1
               IPOS=IPOS+1
               IAA(NCELA+IPOS+1)=IAD(LFIN+J+1)
               IGG(NCELA+IPOS+1)=IDD(LFIN+J+1)
               IGG(KK+NCELA+IPOS+1)=IDD(NCELA+LFIN+J+1)
  141        CONTINUE
             NCLCOR=1+3*J1*(J1-1)
             DO 142 J=NCLCOR-1,IFONC(J1,5),-1
               IPOS=IPOS+1
               IAA(NCELA+IPOS+1)=IAD(LFIN+J+1)
               IGG(NCELA+IPOS+1)=IDD(LFIN+J+1)
               IGG(KK+NCELA+IPOS+1)=IDD(NCELA+LFIN+J+1)
  142        CONTINUE
  143        CONTINUE
  145      CONTINUE
*
*          DUPLICATION DE TURN
*
           IF(LGTURN) THEN
             DO 150 I=1,NCELA
               IGG(2*KK+I)=IDD(2*NCELA+I)
  150        CONTINUE
             LFIN=-NCELAP*IOFZD+NCELA
             J=-1+NCELA
             DO 155 JP=IOFZD,LZZ-1
               LFIN=LFIN-NCELAP
               JAUX2=0
               JAUX1=0
               DO 154 J1=1,NCOUR
               LGPASS=.TRUE.
  152          CONTINUE
               DO 153 J2=JAUX2,JAUX1,-1
                 J=J+1
                 ITURN=IDD(2*NCELA+LFIN+J2+1)
                 IF(ITURN.LE.6) THEN
                   IF(ITURN.EQ.1) THEN
                     IGG(2*KK+J+1)=12
                   ELSEIF(ITURN.EQ.2) THEN
                     IGG(2*KK+J+1)=11
                   ELSEIF(ITURN.EQ.3) THEN
                     IGG(2*KK+J+1)=10
                   ELSEIF(ITURN.EQ.4) THEN
                     IGG(2*KK+J+1)=9
                   ELSEIF(ITURN.EQ.5) THEN
                     IGG(2*KK+J+1)=8
                   ELSEIF(ITURN.EQ.6) THEN
                     IGG(2*KK+J+1)=7
                   ENDIF
                 ELSEIF(ITURN.LE.12) THEN
                   IF(ITURN.EQ.12) THEN
                     IGG(2*KK+J+1)=1
                   ELSEIF(ITURN.EQ.11) THEN
                     IGG(2*KK+J+1)=2
                   ELSEIF(ITURN.EQ.10) THEN
                     IGG(2*KK+J+1)=3
                   ELSEIF(ITURN.EQ.9) THEN
                     IGG(2*KK+J+1)=4
                   ELSEIF(ITURN.EQ.8) THEN
                     IGG(2*KK+J+1)=5
                   ELSEIF(ITURN.EQ.7) THEN
                     IGG(2*KK+J+1)=6
                   ENDIF
                 ELSE
                   CALL XABORT('LHXUNH: EXPECTED VALUE OF TURN ')
                 ENDIF
  153          CONTINUE
               IF(LGPASS) THEN
                 JAUX2=1+3*J1*(J1-1)-1
                 JAUX1=IFONC(J1,5)
                 LGPASS=.FALSE.
                 IF(J1.GT.1) GOTO 152
               ENDIF
               JAUX2=IFONC(J1+1,5)-1
               JAUX1=IFONC(J1+1,0)-1
  154          CONTINUE
  155        CONTINUE
           ENDIF
           IPLANI=LZZ
           LZZ=2*LZZ-IOFZD
         ELSE
           ALLOCATE(IGG(NBASE),IAA(NCELA))
           DO 156 J=1,NBASE
             IGG(J)=IDD(J)
  156      CONTINUE
           DO 157 J=1,NCELA
             IAA(J)=IAD(J)
  157      CONTINUE
           KK=NCELA
         ENDIF
         DEALLOCATE(IDD,IAD)
*
*             RESTORAGE DES DONNEES
*
         NCELA=KK
         MVOSU=0
         MVOLUM=0
         MESH=3*KK
         NSMIN=0
         NSMAX=0
         ISURF=0
         NCOUR=1
         IF(NCELAP.GT.1)NCOUR=IFCOUR(NCELAP)
         ICELC0= IFONC(NCOUR,0)
         ICELC1= IFONC(NCOUR,1)-ICELC0
         MCOU=1
         LCOU=1
         ALLOCATE(SURL(2*LZZ))
         CALL XDISET(SURL,2*LZZ,0)
         IF(LGCELL) THEN
           CALL LCMSIX(IPTRK,'DATA_DUP',1)
           CALL LCMPUT(IPTRK,'GENER0',KK,1,IAA)
           CALL LCMPUT(IPTRK,'GENERATING',KK,1,IGG)
           CALL LCMSIX(IPTRK,' ',2)
           DEALLOCATE(IAA)
           NRAY=0
           NCDP=NCELA-NCELAP
           ICELL1=IFONC(NCOUR,5)
           IP=1
           ISAUX=0
           DO 200 I=0,KK-1
              KRAY=ISEC1(IGG(I+1))
              NRAY=NRAY+KRAY
              I1=I+1
              DO 205 JR=1,KRAY+1
                 KSECT=ISEC2(IGG(I+1))%ILOCR(JR)
                 IF(KSECT.GT.1) THEN
                   MSEC=6*(KSECT-1)
                   MVOSU=MVOSU+MSEC
                   IF(I1.LE.NCELAP)THEN
                     NSMIN=NSMIN+MSEC
                   ELSEIF(I1.GT.NCDP) THEN
                     NSMAX=NSMAX+MSEC
                   ENDIF
                 ELSE
                   MVOSU=MVOSU+1
                   IF(I1.LE.NCELAP)THEN
                     NSMIN=NSMIN+1
                   ELSEIF(I1.GT.NCDP) THEN
                     NSMAX=NSMAX+1
                   ENDIF
                 ENDIF
  205         CONTINUE
              IF(IHEX.EQ.9) THEN
                IF(MCODE.EQ.1) THEN
                  IF(I1.LE.NCELAP) THEN
                    IF(I1.LE.LCOU) THEN
                      VLAT(MCOU)=NSMIN
                      IF(I1.EQ.IFONC(MCOU,5)) VLAT(NCOUR+MCOU)=NSMIN
                    ELSE
                      LCOU=LCOU+6*MCOU
                      MCOU=MCOU+1
                    ENDIF
                  ENDIF
                ELSEIF(MCODE.EQ.2) THEN
                  IF(I1.GT.NCDP) THEN
                    I2=I1-NCDP
                    IF(I2.LE.LCOU) THEN
                      VLAT(MCOU)=NSMAX
                      IF(I2.EQ.IFONC(MCOU,5)) VLAT(NCOUR+MCOU)=NSMAX
                    ELSE
                      LCOU=LCOU+6*MCOU
                      MCOU=MCOU+1
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
              ICOX=INT(AINT(REAL(I1)/REAL(NCELAP)))
              IF(MOD(I1,NCELAP).NE.0)ICOX=ICOX+1
              ICOX=ICOX-1
              ICX=NCELAP*ICOX
              ICELC=ICELC0+ICX
              ICLIM=NCELAP+ICX
              IF(I1.GE.ICELC.AND.I1.LE.ICLIM) THEN
                IF(KSECT.GT.1) THEN
                  ISURF=ISURF+2*(KSECT-1)
                  ISURF0=ISURF
                  IF(MOD(I1-ICELC,ICELC1).EQ.0)ISURF=ISURF+KSECT-1
                  IF(I1.LE.NCELAP*IP) THEN
                    NXX=NCELAP*(IP-1)
                    IF(I1.EQ.ICELC0+NXX) ISAUX=KSECT-1
                    IF(I1.EQ.ICELL1+NXX) SURL(LZZ+IP)=ISURF0-ISAUX
                    SURL(IP)=ISURF
                  ELSE
                    ISAUX=0
                    IP=IP+1
                    NXX=NCELAP*(IP-1)
                    IF(I1.EQ.ICELC0+NXX) ISAUX=KSECT-1
                  ENDIF
                ELSE
                  ISURF=ISURF+2
                  ISURF0=ISURF
                  IF(MOD(I1-ICELC,ICELC1).EQ.0)ISURF=ISURF+1
                  IF(I1.LE.NCELAP*IP) THEN
                    NXX=NCELAP*(IP-1)
                    IF(I1.EQ.ICELC0+NXX) ISAUX=1
                    IF(I1.EQ.ICELL1+NXX) SURL(LZZ+IP)=ISURF0-ISAUX
                    SURL(IP)=ISURF
                  ELSE
                    ISAUX=0
                    IP=IP+1
                    NXX=NCELAP*(IP-1)
                    IF(I1.EQ.ICELC0+NXX) ISAUX=1
                  ENDIF
                ENDIF
              ENDIF
              MVOLUM=MVOLUM+KRAY+1
  200      CONTINUE
           DO 206 I=1,ISTAT9
             DEALLOCATE(ISEC2(I)%ILOCR)
  206      CONTINUE
           MESH=MESH+NRAY
           IF(IPLANZ.EQ.1) NSMAX=NSMIN
         ELSE
           CALL LCMSIX(IPTRK,'DATA_DUP',1)
           CALL LCMPUT(IPTRK,'MIX',KK,1,IGG)
           CALL LCMSIX(IPTRK,' ',2)
           NCDP=NCELA-NCELAP
           MVOLUM=KK
           NCOUR=IFCOUR(NCELAP)
           ICELL1=IFONC(NCOUR,5)
           IP=1
           ISAUX=0
           DO 210 I=0,KK-1
              I1=I+1
              KSECT=ISEC1(IGG(I+1))
              ICOX=INT(AINT(REAL(I1)/REAL(NCELAP)))
              IF(MOD(I1,NCELAP).NE.0)ICOX=ICOX+1
              ICOX=ICOX-1
              ICXX=NCELAP*ICOX
              ICELC=ICELC0+ICXX
              ICLIM=NCELAP+ICXX
              IF(KSECT.GT.1) THEN
                MSEC=6*(KSECT-1)
                MVOSU=MVOSU+MSEC
                IF(I1.LE.NCELAP)THEN
                   NSMIN=NSMIN+MSEC
                ELSEIF(I1.GT.NCDP) THEN
                   NSMAX=NSMAX+MSEC
                ENDIF
                IF(I1.GE.ICELC.AND.I1.LE.ICLIM) THEN
                  ISURF=ISURF+2*(KSECT-1)
                  ISURF0=ISURF
                  IF(MOD(I1-ICELC,ICELC1).EQ.0)ISURF=ISURF+KSECT-1
                  IF(I1.LE.NCELAP*IP) THEN
                    NXX=NCELAP*(IP-1)
                    IF(I1.EQ.ICELC0+NXX) ISAUX=KSECT-1
                    IF(I1.EQ.ICELL1+NXX) SURL(LZZ+IP)=ISURF0-ISAUX
                    SURL(IP)=ISURF
                  ELSE
                    ISAUX=0
                    IP=IP+1
                    NXX=NCELAP*(IP-1)
                    IF(I1.EQ.ICELC0+NXX) ISAUX=KSECT-1
                  ENDIF
                ENDIF
              ELSE
                MVOSU=MVOSU+1
                IF(I1.LE.NCELAP)THEN
                   NSMIN=NSMIN+1
                ELSEIF(I1.GT.NCDP) THEN
                   NSMAX=NSMAX+1
                ENDIF
                IF(I1.GE.ICELC.AND.I1.LE.ICLIM) THEN
                  ISURF=ISURF+2
                  ISURF0=ISURF
                  IF(MOD(I1-ICELC,ICELC1).EQ.0)ISURF=ISURF+1
                  IF(I1.LE.NCELAP*IP) THEN
                    NXX=NCELAP*(IP-1)
                    IF(I1.EQ.ICELC0+NXX) ISAUX=1
                    IF(I1.EQ.ICELL1+NXX) SURL(LZZ+IP)=ISURF0-ISAUX
                    SURL(IP)=ISURF
                  ELSE
                    ISAUX=0
                    IP=IP+1
                    NXX=NCELAP*(IP-1)
                    IF(I1.EQ.ICELC0+NXX) ISAUX=1
                  ENDIF
                ENDIF
              ENDIF
              IF(IHEX.EQ.9) THEN
                IF(MCODE.EQ.1) THEN
                  IF(I1.LE.NCELAP) THEN
                    IF(I1.LE.LCOU) THEN
                      VLAT(MCOU)=NSMIN
                      IF(I1.EQ.IFONC(MCOU,5)) VLAT(NCOUR+MCOU)=NSMIN
                    ELSE
                      LCOU=LCOU+6*MCOU
                      MCOU=MCOU+1
                    ENDIF
                  ENDIF
                ELSEIF(MCODE.EQ.2) THEN
                  IF(I1.GT.NCDP) THEN
                    I2=I1-NCDP
                    IF(I2.LE.LCOU) THEN
                      VLAT(MCOU)=NSMAX
                      IF(I2.EQ.IFONC(MCOU,5)) VLAT(NCOUR+MCOU)=NSMAX
                    ELSE
                      LCOU=LCOU+6*MCOU
                      MCOU=MCOU+1
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
  210      CONTINUE
           IF(IPLANZ.EQ.1) NSMAX=NSMIN
         ENDIF
         DEALLOCATE(ISEC2,ISEC1)
         CALL LCMSIX(IPTRK,'DATA_DUP',1)
         CALL LCMPUT(IPTRK,'SURL_HEX',2*LZZ,1,SURL)
         CALL LCMPUT(IPTRK,'MERGE',KK,1,IGG(KK+1))
         IF(LGTURN)CALL LCMPUT(IPTRK,'TURN',KK,1,IGG(2*KK+1))
         CALL LCMSIX(IPTRK,' ',2)
         DEALLOCATE(IGG,SURL)
*
*              MODIFICATION DU VECTEUR STATE
*
         ISTATE(5)=LZZ
         IPLANZ=LZZ
      ENDIF
      ISTATE(3)=NCELAP
      ISURF=ISURF+NSMIN+NSMAX
      RETURN
      END
