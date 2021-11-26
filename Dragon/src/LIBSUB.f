*DECK LIBSUB
      SUBROUTINE LIBSUB (MAXISO,IPLIB,IPROC,NGRO,NBISO,NLIB,ISONAM,TN,
     1 MASKI,IPRECI,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Production of an internal library with subgroups.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert 
*
*Parameters: input
* MAXISO  maximum number of isotopes permitted.
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* IPROC   type of microlib processing:
*         =2 perform temperature interpolation and build a
*             temperature-independent draglib;
*         =3 perform temperature interpolation and compute calendf-
*             type mathematical probability tables based on bin-type
*             cross-sections;
*         =4 perform temperature interpolation and compute physical
*             probability tables or slowing-down correlated probability
*             tables.
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NLIB    number of independent libraries.
* ISONAM  alias name of isotopes.
* TN      temperature of each isotope.
* MASKI   isotope masks (isotope with index I is process if
*         MASKI(I)=.true.).
* IPRECI  accuracy index for probability tables in CALENDF.
* IMPX    print flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER MAXISO,IPROC,NGRO,NBISO,NLIB,ISONAM(3,NBISO),IPRECI,IMPX
      LOGICAL MASKI(NBISO)
      REAL TN(NBISO)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXDIL=65,MAXED=50,NSTATE=40,MAXESP=4)
      TYPE(C_PTR) JPLIB,KPLIB,IPTMP,JPTMP,KPTMP,IPDRL
      CHARACTER NAMLBT*8,NAMFIL*64,HNISOR*12,HSMG*131,TEXT12*12,
     1 HSHI*12,HVECT(MAXED)*8,HNAMIS*12,HNAMIS2*12,TEXT4*4
      LOGICAL LLENG,LLSHI,LTRANC,LINDEX,MASK2(MAXDIL)
      INTEGER ISOR(3),IPAR(NSTATE),IPAR2(NSTATE),IESP(MAXESP+1)
      REAL DILUT(MAXDIL),EESP(MAXESP+1)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LSHI,NFS,ILLIB,NTFG,NIR,
     1 INAME,IJCEDM,KISONA,KISONR,KTYPE,KNAME,KCOH,KINC,KRSK,KNTFG,
     2 KNIR,KSHIN
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISONRF,ISHINA
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IHLIB
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASKJ
      REAL, ALLOCATABLE, DIMENSION(:) :: GIR,KGIR,KSN,KTN,ENER,EBIN
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SN
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(LSHI(NBISO),NFS(NGRO),ISONRF(3,MAXISO),ISHINA(3,MAXISO),
     2 IHLIB(2,MAXISO,4),ILLIB(MAXISO),NTFG(MAXISO),NIR(MAXISO))
      ALLOCATE(MASKJ(NBISO))
      ALLOCATE(SN(NGRO,NBISO),GIR(MAXISO))
*
      TKSUB=0.0
      TKTAB=0.0
*----
*  CHECK FOR DUPLICATE ISOTOPE NAMES. MASKI(I) NUST BE SET IN SUCH A WAY
*  THAT TWO IDENTICAL ISOTOPES ARE NEVER PROCESSED.
*----
      CALL LCMGET(IPLIB,'ILIBRARYINDX',ILLIB)
      DO 20 I=1,NBISO
      IF(MASKI(I).AND.(ILLIB(I).NE.0)) THEN
         DO 10 J=I+1,NBISO
         IF(MASKI(J).AND.(ISONAM(1,I).EQ.ISONAM(1,J)).AND.(ISONAM(2,I)
     1   .EQ.ISONAM(2,J)).AND.(ISONAM(3,I).EQ.ISONAM(3,J))) THEN
            WRITE (HSMG,300) (ISONAM(I0,I),I0=1,3)
            CALL XABORT(HSMG)
         ENDIF
   10    CONTINUE
      ENDIF
   20 CONTINUE
*----
*  PROCESS THE NON-RESONANT ISOTOPES.
*----
      CALL LCMGET(IPLIB,'ISOTOPESSHI',LSHI)
      DO 35 ISOT=1,NBISO
      MASKJ(ISOT)=MASKI(ISOT).AND.(LSHI(ISOT).EQ.0)
      DO 30 I=1,NGRO
      SN(I,ISOT)=1.0E10
   30 CONTINUE
   35 CONTINUE
      CALL KDRCPU(TK1)
*     -----------------------------------
      CALL LIBLIB(IPLIB,NBISO,MASKJ,IMPX)
*     -----------------------------------
      CALL KDRCPU(TK2)
      TKSUB=TKSUB+(TK2-TK1)
      CALL LCMLEN(IPLIB,'ENERGY',ILENG,ITYLCM)
      LLENG=(ILENG.EQ.NGRO+1)
      CALL LCMLEN(IPLIB,'INDEX',ILENG,ITYLCM)
      LINDEX=(ILENG.NE.0)
*----
*  RECOVER SOME LIBRARY PARAMETERS.
*----
      CALL LCMGET(IPLIB,'STATE-VECTOR',IPAR)
      NL=IPAR(4)
      NED=IPAR(13)
      IF(NED.GT.0) THEN
         IF(NED.GT.MAXED) CALL XABORT('LIBSUB: MAXED OVERFLOW.')
         CALL LCMGTC(IPLIB,'ADDXSNAME-P0',8,NED,HVECT)
      ENDIF
*----
*  RECOVER INFORMATION FROM THE /MICROLIB/ DIRECTORY.
*----
      CALL LCMGET(IPLIB,'ISOTOPERNAME',ISONRF)
      CALL LCMGET(IPLIB,'ILIBRARYTYPE',IHLIB(1,1,1))
      CALL LCMLEN(IPLIB,'ISOTOPESNTFG',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
         CALL LCMGET(IPLIB,'ISOTOPESNTFG',NTFG)
         CALL LCMGET(IPLIB,'ISOTOPESCOH',IHLIB(1,1,2))
         CALL LCMGET(IPLIB,'ISOTOPESINC',IHLIB(1,1,3))
         CALL LCMGET(IPLIB,'ISOTOPESRESK',IHLIB(1,1,4))
      ELSE
         CALL XDISET(NTFG,NBISO,0)
      ENDIF
      CALL LCMLEN(IPLIB,'ISOTOPESHIN',ILENG,ITYLCM)
      LLSHI=(ILENG.GT.0)
      IF(LLSHI) CALL LCMGET(IPLIB,'ISOTOPESHIN',ISHINA)
      CALL LCMLEN(IPLIB,'ISOTOPESNIR',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
         CALL LCMGET(IPLIB,'ISOTOPESNIR',NIR)
         CALL LCMGET(IPLIB,'ISOTOPESGIR',GIR)
      ELSE
         CALL XDISET(NIR,NBISO,0)
         CALL XDRSET(GIR,NBISO,1.0)
      ENDIF
*----
*  PROCESS THE RESONANT ISOTOPES.
*----
      JPLIB=LCMLID(IPLIB,'ISOTOPESLIST',NBISO)
      IMPY=MAX(0,IMPX-5)
      LTRANC=.FALSE.
      DO 200 ISOT=1,NBISO
      IF(MASKI(ISOT).AND.(LSHI(ISOT).NE.0)) THEN
         WRITE(HNAMIS,'(3A4)') (ISONAM(I0,ISOT),I0=1,3)
         IF(IMPX.GT.0) WRITE (6,'(/33H LIBSUB: PROCESSING ISOTOPE/MATER,
     1   5HIAL '',A12,2H''.)') HNAMIS
*
*        RECOVER MULTI-DILUTION INFORMATION.
*
*        FIND THE DILUTION VALUES.
         NDIL=0
         CALL LCMOP(IPTMP,'*TEMPORARY*',0,1,0)
         WRITE(HNISOR,'(3A4)') (ISONRF(I0,ISOT),I0=1,3)
         WRITE(NAMLBT,'(2A4)') IHLIB(1,ISOT,1),IHLIB(2,ISOT,1)
         ALLOCATE(INAME(16*NLIB))
         CALL LCMGET(IPLIB,'ILIBRARYNAME',INAME)
         ILIB=ILLIB(ISOT)
         WRITE(NAMFIL,'(16A4)') (INAME(16*(ILIB-1)+I),I=1,16)
         DEALLOCATE(INAME)
         IF(NAMLBT.EQ.'DRAGON') THEN
            CALL LCMOP(IPDRL,NAMFIL(:12),2,2,0)
            CALL LIBDI1(MAXDIL,IPDRL,HNISOR,NDIL,DILUT)
            CALL LCMCL(IPDRL,1)
         ELSE IF(NAMLBT.EQ.'MATXS') THEN
            CALL LIBDI2(MAXDIL,NAMFIL,HNISOR,NDIL,DILUT)
         ELSE IF(NAMLBT.EQ.'MATXS2') THEN
            CALL LIBDI3(MAXDIL,NAMFIL,HNISOR,NDIL,DILUT)
         ELSE IF(NAMLBT.EQ.'APLIB1') THEN
            WRITE(HSHI,'(3A4)') (ISHINA(I0,ISOT),I0=1,3)
            CALL LIBDI4(MAXDIL,NAMFIL,HSHI,NDIL,DILUT)
         ELSE IF(NAMLBT.EQ.'APLIB2') THEN
            WRITE(HSHI,'(3A4)') (ISHINA(I0,ISOT),I0=1,3)
            IF(HSHI.EQ.' ') THEN
               WRITE (HSMG,'(35HLIBSUB: SELF-SHIELDING ISOTOPE NOT ,
     1         25HDEFINED FOR MAIN ISOTOPE ,A12,1H.)') HNISOR
               CALL XABORT(HSMG)
            ENDIF
            CALL LIBDI5(MAXDIL,NAMFIL,HSHI,NDIL,DILUT)
         ELSE IF(NAMLBT.EQ.'APXSM') THEN
            WRITE(HSHI,'(3A4)') (ISHINA(I0,ISOT),I0=1,3)
            IF(HSHI.EQ.' ') THEN
               WRITE (HSMG,'(35HLIBSUB: SELF-SHIELDING ISOTOPE NOT ,
     1         25HDEFINED FOR MAIN ISOTOPE ,A12,1H.)') HNISOR
               CALL XABORT(HSMG)
            ENDIF
            CALL LIBXS6(MAXDIL,NAMFIL,HSHI,NDIL,DILUT)
         ELSE IF(NAMLBT.EQ.'WIMSAECL') THEN
            WRITE(HSHI,'(3A4)') (ISHINA(I0,ISOT),I0=1,3)
            CALL LIBDI6(MAXDIL,NGRO,NAMFIL,HNISOR,HSHI,NDIL,DILUT)
         ELSE IF(NAMLBT.EQ.'NDAS') THEN
            CALL LIBND7(MAXDIL,NGRO,NAMFIL,HNISOR,NDIL,DILUT)
         ELSE IF(NAMLBT.EQ.'WIMSD4') THEN
            WRITE(HSHI,'(3A4)') (ISHINA(I0,ISOT),I0=1,3)
            CALL LIBDI8(MAXDIL,NGRO,NAMFIL,HNISOR,HSHI,NDIL,DILUT)
         ELSE IF(NAMLBT.EQ.'WIMSE') THEN
            WRITE(HSHI,'(3A4)') (ISHINA(I0,ISOT),I0=1,3)
            CALL LIBDI9(MAXDIL,NGRO,NAMFIL,HNISOR,HSHI,NDIL,DILUT)
         ELSE
            CALL XABORT('LIBSUB: '//NAMLBT//' LIBRARY TREATMENT IS '
     1      //'NOT IMPLEMENTED.')
         ENDIF
         IF(NDIL.EQ.0) GO TO 70
   50    IF(DILUT(1).LT.1.0) THEN
            DO 60 I=2,NDIL+1
            DILUT(I-1)=DILUT(I)
   60       CONTINUE
            NDIL=NDIL-1
            GO TO 50
         ENDIF
         IF(IMPX.GT.4) THEN
            WRITE(6,'(/32H LIBSUB: DILUTIONS FOR ISOTOPE '',A12,
     1      2H'':/(1X,1P,10E12.4))') HNAMIS,(DILUT(I),I=1,NDIL+1)
         ENDIF
   70    IF(NDIL.EQ.0) THEN
            WRITE(HSMG,'(41HLIBSUB: NOT ENOUGH DILUTIONS FOR ISOTOPE ,
     1      A,1H.)') HNAMIS
            CALL XABORT(HSMG)
         ENDIF
*
*        PROCESS THE ISOTOPE FOR EACH DILUTION.
         TEXT12='L_LIBRARY'
         READ(TEXT12,'(3A4)') (ISOR(I),I=1,3)
         CALL LCMPUT(IPTMP,'SIGNATURE',3,3,ISOR)
         DO 80 I=1,NSTATE
         IPAR2(I)=IPAR(I)
   80    CONTINUE
         IPAR2(2)=NDIL+1
         CALL LCMPUT(IPTMP,'STATE-VECTOR',NSTATE,1,IPAR2)
         IF(NED.GT.0) THEN
            ALLOCATE(IJCEDM(2*NED))
            CALL LCMGET(IPLIB,'ADDXSNAME-P0',IJCEDM)
            CALL LCMPUT(IPTMP,'ADDXSNAME-P0',2*NED,3,IJCEDM)
            DEALLOCATE(IJCEDM)
         ENDIF
         IF(LINDEX) THEN
            CALL LCMSIX(IPLIB,'INDEX',1)
            CALL LCMSIX(IPTMP,'INDEX',1)
            CALL LCMEQU(IPLIB,IPTMP)
            CALL LCMSIX(IPTMP,' ',2)
            CALL LCMSIX(IPLIB,' ',2)
         ENDIF
*
*        BUILD A MICROLIB WITH NDIL+1 ISOTOPES.
         ALLOCATE(KISONA(3*(NDIL+1)),KISONR(3*(NDIL+1)),
     1   KTYPE(2*(NDIL+1)),KNAME(NDIL+1))
         IF(NTFG(ISOT).GT.0) THEN
            ALLOCATE(KCOH(2*(NDIL+1)),KINC(2*(NDIL+1)),KRSK(2*(NDIL+1)),
     1      KNTFG(NDIL+1))
         ENDIF
         IF(NIR(ISOT).NE.0) THEN
            ALLOCATE(KGIR(NDIL+1),KNIR(NDIL+1))
         ENDIF
         IF(LLSHI) ALLOCATE(KSHIN(3*(NDIL+1)))
         ALLOCATE(KSN(NGRO*(NDIL+1)),KTN(NDIL+1))
         DO 100 IDIL=1,NDIL+1
         MASK2(IDIL)=.TRUE.
         KISONA(3*(IDIL-1)+1)=ISONAM(1,ISOT)
         KISONA(3*(IDIL-1)+2)=ISONAM(2,ISOT)
         KISONR(3*(IDIL-1)+1)=ISONRF(1,ISOT)
         KISONR(3*(IDIL-1)+2)=ISONRF(2,ISOT)
         KISONR(3*(IDIL-1)+3)=ISONRF(3,ISOT)
         WRITE(TEXT4,'(I4.4)') IDIL
         READ(TEXT4,'(A4)') KISONA(3*(IDIL-1)+3)
         IF(NIR(ISOT).NE.0) THEN
            KGIR(IDIL)=GIR(ISOT)
            KNIR(IDIL)=NIR(ISOT)
         ENDIF
         KTYPE(2*(IDIL-1)+1)=IHLIB(1,ISOT,1)
         KTYPE(2*(IDIL-1)+2)=IHLIB(2,ISOT,1)
         KNAME(IDIL)=ILLIB(ISOT)
         IF(NTFG(ISOT).GT.0) THEN
            KCOH(2*(IDIL-1)+1)=IHLIB(1,ISOT,2)
            KCOH(2*(IDIL-1)+2)=IHLIB(2,ISOT,2)
            KINC(2*(IDIL-1)+1)=IHLIB(1,ISOT,3)
            KINC(2*(IDIL-1)+2)=IHLIB(2,ISOT,3)
            KRSK(2*(IDIL-1)+1)=IHLIB(1,ISOT,4)
            KRSK(2*(IDIL-1)+2)=IHLIB(2,ISOT,4)
            KNTFG(IDIL)=NTFG(ISOT)
         ENDIF
         IF(LLSHI) THEN
            KSHIN(3*(IDIL-1)+1)=ISHINA(1,ISOT)
            KSHIN(3*(IDIL-1)+2)=ISHINA(2,ISOT)
            KSHIN(3*(IDIL-1)+3)=ISHINA(3,ISOT)
         ENDIF
         DO 90 I=1,NGRO
         KSN((IDIL-1)*NGRO+I)=DILUT(IDIL)
   90    CONTINUE
         KTN(IDIL)=TN(ISOT)
  100    CONTINUE
         ALLOCATE(INAME(16*NLIB))
         CALL LCMGET(IPLIB,'ILIBRARYNAME',INAME)
         CALL LCMPUT(IPTMP,'ILIBRARYNAME',16*NLIB,3,INAME)
         DEALLOCATE(INAME)
         CALL LCMPUT(IPTMP,'ISOTOPESUSED',3*(NDIL+1),3,KISONA)
         CALL LCMPUT(IPTMP,'ISOTOPERNAME',3*(NDIL+1),3,KISONR)
         DEALLOCATE(KISONR,KISONA)
         IF(NIR(ISOT).NE.0) THEN
            CALL LCMPUT(IPTMP,'ISOTOPESGIR',NDIL+1,2,KGIR)
            CALL LCMPUT(IPTMP,'ISOTOPESNIR',NDIL+1,1,KNIR)
            DEALLOCATE(KNIR,KGIR)
         ENDIF
         CALL LCMPUT(IPTMP,'ILIBRARYTYPE',2*(NDIL+1),3,KTYPE)
         CALL LCMPUT(IPTMP,'ILIBRARYINDX',NDIL+1,1,KNAME)
         DEALLOCATE(KNAME,KTYPE)
         IF(NTFG(ISOT).GT.0) THEN
            CALL LCMPUT(IPTMP,'ISOTOPESCOH',2*(NDIL+1),3,KCOH)
            CALL LCMPUT(IPTMP,'ISOTOPESINC',2*(NDIL+1),3,KINC)
            CALL LCMPUT(IPTMP,'ISOTOPESRESK',2*(NDIL+1),3,KRSK)
            CALL LCMPUT(IPTMP,'ISOTOPESNTFG',NDIL+1,1,KNTFG)
            DEALLOCATE(KNTFG,KRSK,KINC,KCOH)
         ENDIF
         IF(LLSHI) THEN
            CALL LCMPUT(IPTMP,'ISOTOPESHIN',3*(NDIL+1),3,KSHIN)
            DEALLOCATE(KSHIN)
         ENDIF
         CALL LCMPUT(IPTMP,'ISOTOPESDSN',NGRO*(NDIL+1),2,KSN)
         CALL LCMPUT(IPTMP,'ISOTOPESDSB',NGRO*(NDIL+1),2,KSN)
         CALL LCMPUT(IPTMP,'ISOTOPESTEMP',NDIL+1,2,KTN)
         IF(NED.GT.0) CALL LCMPTC(IPTMP,'ADDXSNAME-P0',8,NED,HVECT)
         DEALLOCATE(KTN,KSN)
*
         CALL KDRCPU(TK1)
*        ------------------------------------
         CALL LIBLIB(IPTMP,NDIL+1,MASK2,IMPY)
*        ------------------------------------
         CALL KDRCPU(TK2)
         TKSUB=TKSUB+(TK2-TK1)
*
*        RECOVER THE SELF-SHIELDING GROUP LIMITS.
         CALL LCMGET(IPTMP,'STATE-VECTOR',IPAR2)
         LTRANC=LTRANC.OR.(IPAR2(5).NE.0)
         IPAR(9)=MIN(IPAR(9),IPAR2(9))
         IPAR(10)=MAX(IPAR(10),IPAR2(10))
         IPAR(16)=MAX(IPAR(16),IPAR2(16))
         IPAR(19)=MAX(IPAR(19),IPAR2(19))
*
*        RECOVER GROUP STRUCTURE
         IF(.NOT.LLENG) THEN
            ALLOCATE(ENER(NGRO+1))
            CALL LCMGET(IPTMP,'ENERGY',ENER)
            CALL LCMPUT(IPLIB,'ENERGY',NGRO+1,2,ENER)
            CALL LCMGET(IPTMP,'DELTAU',ENER)
            CALL LCMPUT(IPLIB,'DELTAU',NGRO,2,ENER)
            DEALLOCATE(ENER)
            LLENG=.TRUE.
         ENDIF
*
*        RECOVER ENERGY-DEPENDENT FISSION SPECTRA
         CALL LCMLEN(IPTMP,'CHI-LIMITS',NBESP,ITYLCM)
         IF(NBESP.GT.0) THEN
            NBESP=NBESP-1
            IF(NBESP.GT.MAXESP) CALL XABORT('LIBSUB: MAXESP OVERFLOW.')
            CALL LCMGET(IPTMP,'CHI-LIMITS',IESP)
            CALL LCMPUT(IPLIB,'CHI-LIMITS',NBESP+1,1,IESP)
            CALL LCMGET(IPTMP,'CHI-ENERGY',EESP)
            CALL LCMPUT(IPLIB,'CHI-ENERGY',NBESP+1,2,EESP)
         ENDIF
*
*        RECOVER BIN TYPE INFORMATION (IF AVAILABLE).
         JPTMP=LCMGID(IPTMP,'ISOTOPESLIST')
         KPLIB=LCMDIL(JPLIB,ISOT) ! set ISOT-th isotope
         CALL LCMLEL(JPTMP,NDIL+1,ILENG,ITYLCM)
         IF(ILENG.EQ.0) THEN
            TEXT12=HNAMIS(1:8)
            WRITE(TEXT12(9:12),'(I4.4)') NDIL+1
            CALL XABORT('LIBSUB: MISSING LIST ITEM FOR '//TEXT12)
         ENDIF
         KPTMP=LCMGIL(JPTMP,NDIL+1) ! set (NDIL+1)-th isotope
         CALL LCMGET(KPTMP,'AWR',AWR)
         CALL LCMLEN(KPTMP,'BIN-NFS',ILENG,ITYLCM)
         IF(ILENG.GT.0) THEN
            CALL LCMGET(KPTMP,'BIN-NFS',NFS)
            CALL LCMPUT(KPLIB,'BIN-NFS',NGRO,1,NFS)
            LBIN=0
            DO 130 I=1,NGRO
            LBIN=LBIN+NFS(I)
  130       CONTINUE
            ALLOCATE(EBIN(LBIN+1))
            CALL LCMLEN(KPTMP,'BIN-ENERGY',ILONG,ITYLCM)
            IF(ILONG.GT.LBIN+1) CALL XABORT('LIBSUB: NFS OVERFLOW.')
            CALL LCMGET(KPTMP,'BIN-ENERGY',EBIN)
            CALL LCMPUT(KPLIB,'BIN-ENERGY',LBIN+1,2,EBIN)
            CALL LCMGET(KPTMP,'BIN-NTOT0',EBIN)
            CALL LCMPUT(KPLIB,'BIN-NTOT0',LBIN,2,EBIN)
            CALL LCMGET(KPTMP,'BIN-SIGS00',EBIN)
            CALL LCMPUT(KPLIB,'BIN-SIGS00',LBIN,2,EBIN)
            CALL LCMLEN(KPTMP,'BIN-NUSIGF',ILENG,ITYLCM)
            IF(ILENG.GT.0) THEN
               CALL LCMGET(KPTMP,'BIN-NUSIGF',EBIN)
               CALL LCMPUT(KPLIB,'BIN-NUSIGF',LBIN,2,EBIN)
            ENDIF
            DEALLOCATE(EBIN)
         ENDIF
*
*        RESET CALENDF MAXIMUM ACCURACY FOR INTERMEDIATE ISOTOPES.
         IPRECJ=IPRECI
         IF((AWR.LT.220.0).AND.(IPRECI.GT.3)) IPRECJ=3
*
         CALL KDRCPU(TK1)
         CALL LIBPTW (KPLIB,IPTMP,IPROC,NGRO,NL,HNAMIS,NED,HVECT,
     1   NDIL,DILUT,AWR,IPRECJ,IMPX)
         CALL KDRCPU(TK2)
         TKTAB=TKTAB+(TK2-TK1)
      ENDIF
  200 CONTINUE
*----
*  COMPUTE CORRELATION INFORMATION BETWEEN PAIRS OF RESONANT ISOTOPES.
*----
      IF((IPROC.EQ.3).OR.(IPROC.EQ.4).OR.(IPROC.EQ.5)) THEN
         DO 220 ISOT=1,NBISO
         WRITE(HNAMIS,'(3A4)') (ISONAM(I0,ISOT),I0=1,3)
         IF(MASKI(ISOT).AND.(LSHI(ISOT).LT.0)) THEN
            DO 210 JSOT=1,ISOT-1
            IF(MASKI(JSOT).AND.(LSHI(ISOT).EQ.LSHI(JSOT))) THEN
               WRITE(HNAMIS2,'(3A4)') (ISONAM(I0,JSOT),I0=1,3)
               IF(IMPX.GT.0) WRITE (6,'(/26H LIBSUB: COMPUTING CORRELA,
     1         41HTION EFFECTS BETWEEN ISOTOPES/MATERIALS '',A12,
     2         7H'' AND '',A12,2H''.)') HNAMIS,HNAMIS2
               CALL LIBCOR(IPLIB,NGRO,ISOT,JSOT,HNAMIS,HNAMIS2)
            ENDIF
  210       CONTINUE
         ENDIF
  220    CONTINUE
      ENDIF
*----
*  RESET IPAR(5) FOR TRANSPORT CORRECTION AND SAVE STATE-VECTOR.
*----
      IF((IPROC.LE.2).AND.LTRANC) IPAR(5)=2
      CALL LCMPUT(IPLIB,'STATE-VECTOR',NSTATE,1,IPAR)
*
      IF(IMPX.GT.0) WRITE(6,'(/28H LIBSUB: CPU TIME IN LIBLIB=,F10.2,
     1 9H SECONDS./9X,19HCPU TIME IN LIBPTW=,F10.2,9H SECONDS.)') TKSUB,
     2 TKTAB
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GIR,SN)
      DEALLOCATE(MASKJ)
      DEALLOCATE(NIR,NTFG,ILLIB,IHLIB,ISHINA,ISONRF,NFS,LSHI)
      RETURN
*
  300 FORMAT(8HLIBSUB: ,3A4,34H IS A DUPLICATE ISOTOPE/MATERIAL N,
     1 4HAME.)
      END
