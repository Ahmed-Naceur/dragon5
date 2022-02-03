*DECK LIBINP
      SUBROUTINE LIBINP (MAXMIX,MAXED,MAXISO,IPLIB,INDREC,IMPX,NBISO,
     1 NGRO,NGT,NL,ITRANC,IPROB,ITIME,NLIB,NGF,IGRMAX,NDEPL,NCOMB,
     2 NEDMAC,NBMIX,NRES,IPROC,IMAC,NDEL,ISOADD,MAXISM,HVECT,IPRECI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read the information related to microscopic cross section libraries.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input/output
* MAXMIX  maximum value of NBMIX.
* MAXED   maximum value of NEDMAC.
* MAXISO  maximum number of isotopes permitted
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature)
* INDREC  type of action:
*         =1 a new library is created; =2 the library is updated;
*         =3 a read-only macrolib is copied in a new library.
* IMPX    print flag.
* NBISO   number of isotopes present in the calculation domain.
* NGRO    number of energy groups.
* NGT     number of energy groups to test.
* NL      number of Legendre orders required in the calculation.
*         NL=1 (for isotropic scattering) or higher.
* ITRANC  type of transport correction: =0 no transport correction
*         =1 Apollo type transport correction; =2 recover from
*         library; =3 WIMS-D type; =4 leakage correction alone.
* IPROB   adjoint macrolib flag:
*         =0 direct problem; =1 adjoint problem.
* ITIME   MATXS type of fission spectrum:
*         =1 steady-state; =2 prompt.
* NLIB    number of independent libraries.
* NGF     number of fast groups without self-shielding.
* IGRMAX  maximum group index with self-shielding.
* NDEPL   number of depleting isotopes (used by EVO:).
* NCOMB   number of depleting mixtures (used by EVO:).
* NEDMAC  number of extra vector edits from matxs.
* NBMIX   number of mixtures defined in the library.
* NRES    number of resonant mixtures (used by SHI:, TONE: or USS:).
* IPROC   type of microlib processing:
*         =0: perform temperature/dilution interpolation;
*         =1: perform temperature interpolation and compute physical
*             probability tables;
*         =2: perform temperature interpolation and build a
*             temperature-independent draglib;
*         =3: perform temperature interpolation and compute calendf
*             type mathematical probability tables based on bin-type
*             cross-sections;
*         =4: compute slowing-down correlated probability tables.
* IMAC    macrolib construction flag:
*         =0 do not compute an embedded macrolib;
*         =1 compute an embedded macrolib.
* NDEL    number of precursor groups for delayed neutrons.
* ISOADD  flag to complete the depletion chain:
*         =0 complete; =1 do not complete.
* MAXISM  maximum number of isotopes per mixture.
* HVECT   matxs names of the extra vector edits.
* IPRECI  accuracy index for probability tables in CALENDF.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER MAXMIX,MAXED,MAXISO,INDREC,IMPX,NBISO,NGRO,NGT,NL,ITRANC,
     > IPROB,ITIME,NLIB,NGF,IGRMAX,NDEPL,NCOMB,NEDMAC,NBMIX,NRES,IPROC,
     > IMAC,NDEL,ISOADD,MAXISM,IPRECI
      CHARACTER*(*) HVECT(MAXED)
*----
*  LOCAL PARAMETERS
*----
      PARAMETER (IOUT=6,NHOBL=16,MAXPAR=5,NSTATE=40)
      TYPE(C_PTR) JPLIB
      DOUBLE PRECISION DBLINP
      CHARACTER TEXT4*4,HOBL(NHOBL)*8,TEXT12*12,HSMG*131,NAMFIL*64,
     >          NAMFIL2*64,NAMLBT*8,NAMLCM*12,NAMMY*12
      LOGICAL LNEW,EMPTY,LCM,LSET
      INTEGER KCHAR(2),ISTATE(NSTATE)
      REAL TMPDAY(3)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISOMIX,NTFG,LSHI,NIR,ILLIB,
     > IEVOL,ITYP,INAME,KGAS
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISONAM,ISONRF,ISHINA
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IHLIB
      REAL, ALLOCATABLE, DIMENSION(:) :: DENISO,DENMIX,TMPISO,SNISO,
     > SBISO,GIR,TMPMIX,GSN,GSB
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASK,MASKI,MASKL
*----
*  DATA STATEMENTS
*----
      SAVE HOBL
      DATA HOBL  /'NFTOT   ','NG      ','N2N     ','N3N     ',
     >            'N4N     ','NA      ','NP      ','N2A     ',
     >            'NNP     ','ND      ','NT      ','TRANC   ',
     >            'BSTC    ','BSTR    ','CSTC    ','CSTR    '/
*----
*  SCRATCH STORAGE ALLOCATION
*   ISONAM  alias name of each isotope
*   ISONRF  library name of each isotope
*   ISOMIX  mix number of each isotope. can be zero
*   DENISO  density of each isotope
*   DENMIX  mixture density
*   MASK    mixture masks
*   TMPISO  temperature of each isotope
*   ISHINA  self-shielding name of each isotope
*   SNISO   dilution cross section of each isotope. a value of 1.0E10
*           is used for infinite dilution
*   SBISO   dilution cross section of each isotope used with Livolant-
*           Jeanpierre normalization
*   NTFG    number of thermal groups where the thermal inelastic
*           correction is applied
*   LSHI    resonant region number associated with i-th isotope.
*           Infinite dilution will be assumed if LSHI(I)=0. A negative
*           value is indicating correlation of cross sections with all
*           isotopes sharing the same LSHI value
*   GIR     Goldstein-Cohen IR parameter of each isotope
*   NIR     Goldstein-Cohen IR cutoff energy index. Use IR approximation
*           for groups with index.ge.nir; Use library value if NIR=0
*   MASKI   isotope masks
*   TMPMIX  mixture temperature
*   IHLIB   isotope options
*   ILLIB   xs library index for each isotope (.le.NLIB)
*   IEVOL   flag making an isotope non-depleting:
*           =1 to force an isotope to be non-depleting;
*           =2 to force an isotope to be depleting;
*           =3 to force an isotope to be at saturation
*   ITYP    isotopic type:
*           =1: the isotope is not fissile and not a fission product;
*           =2: the isotope is fissile; =3: is a fission product
*   INAME   library name
*   KGAS    state of mixture (used for stopping power correction):
*           =0: solid or liquid;
*           =1: gas
*----
      ALLOCATE(ISONAM(3,MAXISO),ISONRF(3,MAXISO),ISOMIX(MAXISO),
     > ISHINA(3,MAXISO),NTFG(MAXISO),LSHI(MAXISO),NIR(MAXISO),
     > IHLIB(2,MAXISO,4),ILLIB(MAXISO),IEVOL(MAXISO),ITYP(MAXISO),
     > KGAS(MAXMIX))
      ALLOCATE(DENISO(MAXISO),DENMIX(MAXMIX),TMPISO(MAXISO),
     > SNISO(MAXISO),SBISO(MAXISO),GIR(MAXISO),TMPMIX(MAXMIX))
      ALLOCATE(MASK(MAXMIX),MASKI(MAXISO))
*----
*  INITIALIZATIONS.
*----
      KEVOL=0
      IF(NGT .NE. NGRO) THEN
        WRITE(IOUT,400) NGT,NGRO
      ENDIF
      IF((INDREC.EQ.2).AND.(NBISO.GT.0)) THEN
*        THE LIBRARY IS UPDATED. READ OLD LIBRARY INFORMATION.
         IF(NBISO.GT.MAXISO) CALL XABORT('LIBINP: MAXISO OVERFLOW.')
         IF(NBMIX.GT.MAXMIX) CALL XABORT('LIBINP: MAXMIX OVERFLOW.')
         NNMIX=NBMIX
         CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
         CALL LCMGET(IPLIB,'ISOTOPESUSED',ISONAM)
         CALL LCMLEN(IPLIB,'ISOTOPERNAME',ILENG,ITYLCM)
         IF(ILENG.GT.0) THEN
            CALL LCMGET(IPLIB,'ISOTOPERNAME',ISONRF)
         ELSE
            CALL LCMGET(IPLIB,'ISOTOPESUSED',ISONRF)
         ENDIF
         TEXT4=' '
         READ(TEXT4,'(A4)') IVOID
         CALL XDISET(IHLIB,8*MAXISO,IVOID)
         CALL XDISET(ILLIB,MAXISO,0)
         CALL XDISET(ISHINA,3*MAXISO,IVOID)
         CALL LCMLEN(IPLIB,'ILIBRARYTYPE',ILENG,ITYLCM)
         IF(ILENG.GT.0) THEN
            CALL LCMGET(IPLIB,'ILIBRARYTYPE',IHLIB(1,1,1))
            CALL LCMGET(IPLIB,'ILIBRARYINDX',ILLIB)
         ENDIF
         CALL LCMLEN(IPLIB,'ISOTOPESNTFG',ILENG,ITYLCM)
         IF(ILENG.GT.0) THEN
            CALL LCMGET(IPLIB,'ISOTOPESNTFG',NTFG)
            CALL LCMGET(IPLIB,'ISOTOPESCOH',IHLIB(1,1,2))
            CALL LCMGET(IPLIB,'ISOTOPESINC',IHLIB(1,1,3))
         ELSE
            CALL XDISET(NTFG,MAXISO,0)
         ENDIF
        CALL LCMLEN(IPLIB,'ISOTOPESRESK',ILENG,ITYLCM)
         IF(ILENG.GT.0) THEN
            CALL LCMGET(IPLIB,'ISOTOPESRESK',IHLIB(1,1,4))
         ELSE
            NAMLBT=','
            DO ISOT=1,NBISO
               READ(NAMLBT,'(2A4)') IHLIB(1,ISOT,4),IHLIB(2,ISOT,4)
            ENDDO
         ENDIF
         CALL LCMLEN(IPLIB,'ISOTOPESHIN',ILENG,ITYLCM)
         IF(ILENG.GT.0) CALL LCMGET(IPLIB,'ISOTOPESHIN',ISHINA)
         CALL LCMLEN(IPLIB,'ISOTOPESSHI',ILENG,ITYLCM)
         IF(ILENG.GT.0) THEN
            CALL LCMGET(IPLIB,'ISOTOPESSHI',LSHI)
         ELSE
            CALL XDISET(LSHI,MAXISO,0)
         ENDIF
         CALL LCMLEN(IPLIB,'ISOTOPESNIR',ILENG,ITYLCM)
         IF(ILENG.GT.0) THEN
            CALL LCMGET(IPLIB,'ISOTOPESGIR',GIR)
            CALL LCMGET(IPLIB,'ISOTOPESNIR',NIR)
         ELSE
            CALL XDRSET(GIR,MAXISO,1.0)
            CALL XDISET(NIR,MAXISO,0)
         ENDIF
         CALL LCMGET(IPLIB,'ISOTOPESDENS',DENISO)
         CALL LCMGET(IPLIB,'ISOTOPESMIX',ISOMIX)
         CALL LCMGET(IPLIB,'ISOTOPESTEMP',TMPISO)
         CALL LCMLEN(IPLIB,'ISOTOPESTODO',ILENG,ITYLCM)
         IF(ILENG.GT.0) THEN
            CALL LCMGET(IPLIB,'ISOTOPESTODO',IEVOL)
         ELSE
            CALL XDISET(IEVOL,MAXISO,0)
         ENDIF
         CALL LCMLEN(IPLIB,'ISOTOPESTYPE',ILENG,ITYLCM)
         IF(ILENG.GT.0) THEN
            CALL LCMGET(IPLIB,'ISOTOPESTYPE',ITYP)
         ELSE
            CALL XDISET(ITYP,MAXISO,1)
         ENDIF
         CALL LCMLEN(IPLIB,'MIXTUREGAS',ILENG,ITYLCM)
         IF(ILENG.EQ.NBMIX) THEN
            CALL LCMGET(IPLIB,'MIXTUREGAS',KGAS)
         ELSE
            CALL XDISET(KGAS,NBMIX,0)
         ENDIF
         DO 40 IIIMIX=1,MAXMIX
         DENMIX(IIIMIX)=-1.0
         DO 30 IIISO=1,NBISO
         IF(ISOMIX(IIISO).EQ.IIIMIX) THEN
           TMPMIX(IIIMIX)=TMPISO(IIISO)
           GO TO 40
         ENDIF
   30    CONTINUE
         TMPMIX(IIIMIX)=-1.0
   40    CONTINUE
         IF(ISTATE(17).EQ.-1) THEN
            DO 44 IISO=1,NBISO
            MASKI(IISO)=.TRUE.
   44       CONTINUE
         ELSE
            DO 45 IISO=1,NBISO
            MASKI(IISO)=.FALSE.
   45       CONTINUE
         ENDIF
         CALL XDRSET(SNISO,NBISO,0.0)
         CALL XDRSET(SBISO,NBISO,0.0)
      ELSE
         NELSN=0
         NNMIX=0
         DO 50 IIIMIX=1,MAXMIX
         DENMIX(IIIMIX)=-1.0
         TMPMIX(IIIMIX)=-1.0
         KGAS(IIIMIX)=0
   50    CONTINUE
      ENDIF
*----
*  READ THE SPECIFICATION FOR EACH ISOTOPE.
*----
      TEXT12='MIXS'
      JLIB=0
      LSET=.TRUE.
   60 IF(TEXT12.EQ.'MIXS') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DBLINP)
         IF((INDIC.EQ.3).AND.(TEXT4.EQ.';')) THEN
            DO 65 IISO=1,NBISO
            MASKI(IISO)=.TRUE.
   65       CONTINUE
            GO TO 100
         ENDIF
         IF((INDIC.NE.3).OR.(TEXT4.NE.'LIB:'))
     >     CALL XABORT('LIBINP: KEYWORD LIB: EXPECTED')
         CALL REDGET(INDIC,NITMA,FLOTT,NAMLBT,DBLINP)
         IF(INDIC.NE.3)
     >     CALL XABORT('LIBINP: CHARACTER LIBRARY NAME REQUIRED.')
         IF( (NAMLBT.NE.'MATXS' ).AND.(NAMLBT.NE.'MATXS2').AND.
     >       (NAMLBT.NE.'APLIB1').AND.(NAMLBT.NE.'APLIB2').AND.
     >       (NAMLBT.NE.'DRAGON').AND.(NAMLBT.NE.'WIMSAECL').AND.
     >       (NAMLBT.NE.'WIMSD4').AND.(NAMLBT.NE.'WIMSE' ).AND.
     >       (NAMLBT.NE.'NDAS'  ).AND.(NAMLBT.NE.'APXSM' ).AND.
     >       (NAMLBT.NE.'MICROLIB')) THEN
           WRITE(HSMG,'(29HLIBINP: INVALID LIBRARY TYPE ,A8)') NAMLBT
           CALL XABORT(HSMG)
         ENDIF
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
         IF((INDIC.NE.3).OR.(TEXT12.NE.'FIL:'))
     >     CALL XABORT('LIBINP: FIL: EXPECTED.')
         NAMFIL=' '
         CALL REDGET(INDIC,NITMA,FLOTT,NAMFIL,DBLINP)
         IF(INDIC.NE.3) CALL XABORT('LIBINP: CHARACTER DATA EXPECTED'//
     >   '(1).')
         CALL LIBNRG(IPLIB,NAMLBT,NAMFIL,NGRO,NGT)
         ALLOCATE(INAME(16*(NLIB+1)))
         IF(NLIB.GT.0) CALL LCMGET(IPLIB,'ILIBRARYNAME',INAME)
         DO 66 ILIB=1,NLIB
           WRITE(NAMFIL2,'(16A4)') (INAME(16*(ILIB-1)+I),I=1,16)
           IF(NAMFIL2.EQ.NAMFIL) THEN
             JLIB=ILIB
             DEALLOCATE(INAME)
             GO TO 67
           ENDIF
   66    CONTINUE
         NLIB=NLIB+1
         READ(NAMFIL,'(16A4)') (INAME(16*(NLIB-1)+I),I=1,16)
         CALL LCMPUT(IPLIB,'ILIBRARYNAME',16*NLIB,3,INAME)
         DEALLOCATE(INAME)
         JLIB=NLIB
   67    CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
         IF(INDIC.NE.3) CALL XABORT('LIBINP: CHARACTER DATA EXPECTED'//
     >   '(2).')
         GO TO 60
      ELSE IF(TEXT12.EQ.';') THEN
         GO TO 100
      ELSE IF(TEXT12.EQ.'MIX') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
         IF(INDIC.EQ.1) THEN
           NNMIX=NITMA
           CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
         ELSE
           NNMIX=NNMIX+1
         ENDIF
         IF(NNMIX.GT.MAXMIX) THEN
           CALL XABORT('LIBINP: MIX NUMBER LARGER THAN MAXMIX.')
         ELSE IF(NNMIX.LE.0) THEN
           CALL XABORT('LIBINP: MIX NUMBER .LE. 0.')
         ENDIF
         NBMIX=MAX(NNMIX,NBMIX)
         IF(INDIC.EQ.3) THEN
*----
*  THIS MIXTURE IS A COMBINATION OF OTHER MIXTURES
*----
           IF(TEXT12.EQ.'COMB') THEN
              CALL LCMLEN(IPLIB,'MACROLIB',ILONG,ITYLCM)
              IF((ILONG.NE.0).AND.LSET) THEN
*               perform a reset of the macrolib to be safe
                CALL LCMSIX(IPLIB,'MACROLIB',1)
                CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
                ISTATE(4)=0
                CALL LCMPUT(IPLIB,'STATE-VECTOR',NSTATE,1,ISTATE)
                CALL LCMSIX(IPLIB,' ',2)
                CALL XDLSET(MASKI,NBISO,.TRUE.)
                LSET=.FALSE.
              ENDIF
              CALL LCMINF(IPLIB,NAMLCM,NAMMY,EMPTY,ILONG,LCM)
              VOLTOT=0.0
  70          VOLFRA=0.0
              MIXCMB=0
              CALL REDGET(INDIC,MIXCMB,FLOTT,TEXT12,DBLINP)
              IF(INDIC.EQ.3) THEN
                 IF(VOLTOT.EQ.0.0) CALL XABORT('LIBINP: TOTAL VOLUME F'
     >           //'RACTION OF 0.0 IS ILLEGAL.')
                 GO TO 60
              ENDIF
              IF(INDIC.EQ.2) CALL XABORT('LIBINP: MIXTURE NUMBER MISSI'
     >        //'NG FOR COMBINATION.')
              CALL REDGET(INDIC,NITMA,VOLFRA,TEXT12,DBLINP)
              IF((INDIC.EQ.1).OR.(INDIC.EQ.3)) CALL XABORT('LIBINP: VO'
     >        //'LUME FRACTION MISSING FOR COMBINATION.')
              IF(VOLFRA.EQ.0.0) CALL XABORT('LIBINP: INDIVIDUAL VOLUME'
     >        //' FRACTION OF 0.0 IS ILLEGAL.')
              CALL LIBCMB(MAXMIX,MAXISO,NBISO,NEWISO,NNMIX,MIXCMB,
     1        VOLTOT,VOLFRA,DENMIX,ISONAM,ISONRF,ISHINA,ISOMIX,IHLIB,
     2        ILLIB,DENISO,TMPISO,LSHI,SNISO,SBISO,NTFG,NIR,GIR,MASKI,
     3        IEVOL,ITYP)
              GO TO 70
           ELSE
              CALL XABORT('LIBINP: ONLY COMB KEYWORD CAN FOLLOW MIXTUR'
     1       //'E NUMBER.')
           ENDIF
         ELSE
           IF(INDIC.NE.2) CALL XABORT('LIBINP: REAL NUMBER EXPECTED.')
           TMPMIX(NNMIX)=FLOTT
           CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
           IF(INDIC.EQ.2) THEN
             DENMIX(NNMIX)=FLOTT
             CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
           ENDIF
           IF(INDIC.NE.3) CALL XABORT('LIBINP: CHARACTER DATA EXPECT'//
     >     'ED(3).')
           IF(TEXT12.EQ.'NOEV') THEN
             CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
             IF(INDIC.NE.3) CALL XABORT('LIBINP: CHARACTER DATA EXPE'//
     >       'CTED(4).')
             KEVOL=1
           ELSE IF(TEXT12.EQ.'EVOL') THEN
             CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
             IF(INDIC.NE.3) CALL XABORT('LIBINP: CHARACTER DATA EXPE'//
     >       'CTED(5).')
             KEVOL=2
           ELSE
             KEVOL=0
           ENDIF
           IF(TEXT12.EQ.'NOGAS') THEN
             CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
             IF(INDIC.NE.3) CALL XABORT('LIBINP: CHARACTER DATA EXPE'//
     >       'CTED(6).')
             KGAS(NNMIX)=0
           ENDIF
           IF(TEXT12.EQ.'GAS') THEN
             CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
             IF(INDIC.NE.3) CALL XABORT('LIBINP: CHARACTER DATA EXPE'//
     >       'CTED(7).')
             KGAS(NNMIX)=1
           ENDIF
           IF((TEXT12.EQ.'MIX').OR.(TEXT12.EQ.'MIXS').OR.
     >        (TEXT12.EQ.';')) THEN
              DO 80 IISO=1,NBISO
              IF(ISOMIX(IISO).EQ.NNMIX) THEN
                 TMPISO(IISO)=TMPMIX(NNMIX)
                 MASKI(IISO)=.TRUE.
              ENDIF
   80         CONTINUE
           ENDIF
         ENDIF
         GO TO 60
      ENDIF
      READ(TEXT12,'(2A4)') KCHAR(1),KCHAR(2)
      DO 81 I=1,NBISO
      IF((KCHAR(1).EQ.ISONAM(1,I)).AND.(KCHAR(2).EQ.ISONAM(2,I)).AND.
     > (NNMIX.EQ.ISOMIX(I))) THEN
*        UPDATE AN EXISTING ISOTOPE.
         NEWISO=I
         LNEW=.FALSE.
         GO TO 82
      ENDIF
   81 CONTINUE
      LNEW=.TRUE.
      NBISO=NBISO+1
      NEWISO=NBISO
      IF(NBISO.GT.MAXISO) CALL XABORT('LIBINP: MAXISO TOO SMALL.')
      READ(TEXT12,'(3A4)') (ISONAM(I0,NBISO),I0=1,3)
      READ(TEXT12,'(3A4)') (ISONRF(I0,NBISO),I0=1,3)
      TEXT12=' '
      READ(TEXT12,'(3A4)') (ISHINA(I0,NBISO),I0=1,3)
      READ(TEXT12,'(2A4)') IHLIB(1,NBISO,2),IHLIB(2,NBISO,2)
      READ(TEXT12,'(2A4)') IHLIB(1,NBISO,3),IHLIB(2,NBISO,3)
      READ(TEXT12,'(2A4)') IHLIB(1,NBISO,4),IHLIB(2,NBISO,4)
      NTFG(NBISO)=0
      LSHI(NBISO)=0
      GIR(NBISO)=1.0
      NIR(NBISO)=0
      ISOMIX(NBISO)=NNMIX
      DENISO(NBISO)=0.0
      SNISO(NBISO)=1.0E10
      SBISO(NBISO)=1.0E10
      IEVOL(NBISO)=KEVOL
      ITYP(NBISO)=1
*
   82 MASKI(NEWISO)=.TRUE.
      READ(NAMLBT,'(2A4)') IHLIB(1,NEWISO,1),IHLIB(2,NEWISO,1)
      ILLIB(NEWISO)=JLIB
      TMPISO(NEWISO)=TMPMIX(NNMIX)
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
      IF(INDIC.EQ.3.AND.TEXT12.EQ.'=') THEN
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
        IF(INDIC.EQ.3) THEN
           READ(TEXT12,'(3A4)') (ISONRF(I0,NEWISO),I0=1,3)
        ELSE
           CALL XABORT('LIBINP: LIBRARY ISOTOPE NAME MISSING AFTER =')
        ENDIF
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
      ENDIF
      IF(INDIC.EQ.3) THEN
*        USE THE VALUES ALREADY THERE.
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
         IF(INDIC.NE.3) CALL XABORT('LIBINP: CHARACTER DATA EXPECTED'//
     >   '(8).')
         GO TO 60
      ELSE IF(INDIC.NE.2) THEN
        CALL XABORT('LIBINP: ISOTOPIC DENSITY OR WEIGHT PERCENT EXPECT'
     >  //'ED.')
      ENDIF
      IF((.NOT.LNEW).AND.(DENMIX(NNMIX).NE.-1.0).AND.(ABS(DENISO(NEWISO)
     1 -FLOTT).GT.1.0E-4)) THEN
        CALL XABORT('LIBINP: PERTURBATION OF THE WEIGHT PERCENTS IS FOR'
     1  //'BIDDEN.')
      ENDIF
      DENISO(NEWISO)=FLOTT
   90 CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
      IF(INDIC.EQ.2) THEN
         SNISO(NEWISO)=FLOTT
         SBISO(NEWISO)=FLOTT
      ELSE IF(INDIC.EQ.1) THEN
         LSHI(NEWISO)=NITMA
         NRES=MAX(NRES,NITMA)
         IF(IPROC.EQ.3) THEN
            NIR(NEWISO)=1
            GIR(NEWISO)=-998.0
         ELSE IF(IPROC.EQ.4) THEN
            NIR(NEWISO)=1
            GIR(NEWISO)=-999.0
         ELSE IF(IPROC.EQ.5) THEN
            NIR(NEWISO)=1
            GIR(NEWISO)=-1000.0
         ENDIF
      ELSE IF(TEXT12.EQ.'INF') THEN
         SNISO(NEWISO)=1.0E10
         SBISO(NEWISO)=1.0E10
      ELSE IF(TEXT12.EQ.'SHIB') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
         IF(INDIC.NE.3) CALL XABORT('LIBINP: CHARACTER DATA EXPECTED'//
     >   '(9).')
         READ(TEXT12,'(3A4)') (ISHINA(I0,NEWISO),I0=1,3)
      ELSE IF(TEXT12.EQ.'THER') THEN
         CALL REDGET(INDIC,NTFG(NEWISO),FLOTT,TEXT12,DBLINP)
         IF(INDIC.NE.1) CALL XABORT('LIBINP: NUMBER OF THERMALIZED '//
     >   'GROUPS REQUIRED.')
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
         IF(INDIC.NE.3) CALL XABORT('LIBINP: CHARACTER DATA EXPECTED'//
     >   '(10).')
         READ(TEXT12,'(2A4)') IHLIB(1,NEWISO,3),IHLIB(2,NEWISO,3)
      ELSE IF(TEXT12.EQ.'TCOH') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
         IF(INDIC.NE.3) CALL XABORT('LIBINP: CHARACTER DATA EXPECTED'//
     >   '(11).')
         READ(TEXT12,'(2A4)') IHLIB(1,NEWISO,2),IHLIB(2,NEWISO,2)
      ELSE IF(TEXT12.EQ.'RESK') THEN
         TEXT12='RESK'
         READ(TEXT12,'(2A4)') IHLIB(1,NEWISO,4),IHLIB(2,NEWISO,4)
      ELSE IF(TEXT12.EQ.'DBYE') THEN
         CALL REDGET(INDIC,NITMA,TMPISO(NEWISO),TEXT12,DBLINP)
         IF(INDIC.NE.2) CALL XABORT('LIBINP: REAL DATA EXPECTED.')
      ELSE IF(TEXT12.EQ.'CORR') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
         IF(INDIC.NE.1) CALL XABORT('LIBINP: INTEGER DATA EXPECTED.')
         LSHI(NEWISO)=-NITMA
         NRES=MAX(NRES,NITMA)
      ELSE IF(TEXT12.EQ.'IRSET') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
         IF(INDIC.EQ.2) THEN
            GIR(NEWISO)=FLOTT
            IF((IPROC.EQ.3).AND.(FLOTT.NE.1.0)) CALL XABORT('LIBINP: P'
     >      //'T MAIN OPTION NOT EXPECTED.')
            IF((IPROC.EQ.4).AND.(FLOTT.NE.1.0)) CALL XABORT('LIBINP: P'
     >      //'TSL MAIN OPTION NOT EXPECTED.')
         ELSE IF(INDIC.EQ.3) THEN
            IF(TEXT12.EQ.'PT') THEN
               IF(IPROC.NE.3) CALL XABORT('LIBINP: PT MAIN OPTION NOT '
     >         //'SET.')
               GIR(NEWISO)=-998.0
            ELSE IF(TEXT12.EQ.'PTSL') THEN
               IF(IPROC.NE.4) CALL XABORT('LIBINP: PTSL MAIN OPTION NO'
     >         //'T SET.')
               GIR(NEWISO)=-999.0
            ELSE IF(TEXT12.EQ.'PTMC') THEN
               IF(IPROC.NE.5) CALL XABORT('LIBINP: PTMC MAIN OPTION NO'
     >         //'T SET.')
               GIR(NEWISO)=-1000.0
            ELSE
               CALL XABORT('LIBINP: PT, PTSL OR PTMC EXPECTED.')
            ENDIF
         ELSE
            CALL XABORT('LIBINP: REAL OR CHARACTER DATA EXPECTED.')
         ENDIF
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
         IF(INDIC.EQ.1) THEN
            IF((NITMA.LT.0).OR.(NITMA.GT.NGRO)) CALL
     >      XABORT('LIBINP: INVALID VALUE OF NIR.')
            NIR(NEWISO)=NITMA
         ELSE IF((INDIC.EQ.3).AND.(TEXT12.EQ.'NONE')) THEN
            NIR(NEWISO)=NGRO+1
         ELSE
            CALL XABORT('LIBINP: NONE OR INTEGER DATA EXPECTED.')
         ENDIF
      ELSE IF(TEXT12.EQ.'NOEV') THEN
         IEVOL(NEWISO)=1
      ELSE IF(TEXT12.EQ.'EVOL') THEN
         IEVOL(NEWISO)=2
      ELSE IF(TEXT12.EQ.'SAT') THEN
         IEVOL(NEWISO)=3
      ELSE
         IF(INDIC.NE.3) CALL XABORT('LIBINP: CHARACTER DATA EXPECTED'//
     >   '(12).')
         GO TO 60
      ENDIF
      GO TO 90
*----
*  INCLUDE SOME DEFAULT EXTRA EDITS.
*----
  100 IF((NGRO.EQ.0).OR.(NGT.EQ.0)) CALL XABORT('LIBINP: NUMBER OF GRO'
     > //'UPS REQUIRED.')
      DO 120 I=1,NHOBL
      DO 110 IED=1,NEDMAC
      IF(HVECT(IED).EQ.HOBL(I)) GO TO 120
  110 CONTINUE
      NEDMAC=NEDMAC+1
      IF(NEDMAC.GT.MAXED) CALL XABORT('LIBINP: TOO MANY EXTRA EDITS R'
     > //'EQUESTED.')
      HVECT(NEDMAC)=HOBL(I)
  120 CONTINUE
*----
*  ADD THE MISSING ISOTOPES FROM THE DEPLETION CHAIN.
*----
      IF((NDEPL.NE.0).AND.(ISOADD.EQ.0)) THEN
         NBISOL=NBISO
         CALL LCMSIX(IPLIB,'DEPL-CHAIN',1)
         CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
         IF(ISTATE(1).NE.NDEPL) CALL XABORT('LIBINP: INVALID NUMBER OF'
     >   //' DEPLETING ISOTOPES.')
         NFISS=ISTATE(2)
         NSUPF=ISTATE(5)
         NSUPS=ISTATE(7)
         NREAC=ISTATE(8)
         NPAR=ISTATE(9)
         CALL LIBEAD(IPLIB,MAXISO,MAXMIX,IMPX,NDEPL,NFISS,NSUPS,
     1   NREAC,NPAR,NBISO,ISONAM,ISONRF,IHLIB,ILLIB,ISOMIX,TMPISO,
     2   IEVOL,ITYP,NCOMB)
         CALL LCMSIX(IPLIB,' ',2)
*
         DO 140 ISOT=NBISOL+1,NBISO
         SNISO(ISOT)=1.0E10
         SBISO(ISOT)=1.0E10
         DENISO(ISOT)=0.0
         NTFG(ISOT)=0
         TEXT12=' '
         READ(TEXT12,'(3A4)') (ISHINA(I0,ISOT),I0=1,3)
         READ(TEXT12,'(2A4)') IHLIB(1,ISOT,2),IHLIB(2,ISOT,2)
         READ(TEXT12,'(2A4)') IHLIB(1,ISOT,3),IHLIB(2,ISOT,3)
         READ(TEXT12,'(2A4)') IHLIB(1,ISOT,4),IHLIB(2,ISOT,4)
         LSHI(ISOT)=0
         GIR(ISOT)=1.0
         NIR(ISOT)=0
         MASKI(ISOT)=.TRUE.
  140    CONTINUE
      ENDIF
*----
*  SET THE MIXTURE MASKS.
*----
      DO 170 I=1,NBMIX
        MASK(I)=.FALSE.
        DO 150 JJ=1,NBISO
          IF((ISOMIX(JJ).EQ.I).AND.MASKI(JJ)) THEN
            MASK(I)=.TRUE.
            GO TO 160
          ENDIF
  150   CONTINUE
  160   CONTINUE
  170 CONTINUE
*----
*  FIND AND NAME DISTINCT ISOTOPES.
*----
      DO 200 I=1,NBISO
      IF(MASKI(I).AND.(ILLIB(I).NE.0)) THEN
*       CATENATE THE 4-DIGIT MIXTURE SUFFIX.
        DO 190 J=1,I-1
        IF((ISONAM(1,I).NE.ISONAM(1,J)).OR.(ISONAM(2,I).NE.ISONAM(2,J)))
     >  GO TO 190
        IF((ISONRF(1,I).NE.ISONRF(1,J)).OR.(ISONRF(2,I).NE.ISONRF(2,J))
     >  .OR.(ISONRF(3,I).NE.ISONRF(3,J))) GO TO 190
        IF((ISHINA(1,I).NE.ISHINA(1,J)).OR.(ISHINA(2,I).NE.ISHINA(2,J))
     >  .OR.(ISHINA(3,I).NE.ISHINA(3,J))) GO TO 190
        IF((LSHI(I).NE.0).AND.(LSHI(J).NE.0).AND.(DENISO(I).EQ.0.0)
     >  .AND.(DENISO(J).NE.0.0)) GO TO 190
        IF((LSHI(I).NE.0).AND.(LSHI(J).NE.0).AND.(DENISO(I).NE.0.0)
     >  .AND.(DENISO(J).EQ.0.0)) GO TO 190
        DO 186 IOP=1,4
        DO 185 I0=1,2
          IF(IHLIB(I0,I,IOP).NE.IHLIB(I0,J,IOP)) GO TO 190
  185   CONTINUE
  186   CONTINUE
        IF(ILLIB(I).NE.ILLIB(J)) GO TO 190
        IF((NTFG(I).NE.NTFG(J)).OR.(GIR(I).NE.GIR(J)).OR.
     >  (NIR(I).NE.NIR(J)).OR.(TMPISO(I).NE.TMPISO(J))) GO TO 190
        IF(((LSHI(I).EQ.0).AND.(LSHI(J).EQ.0))
     >    .OR.((IPROC.NE.0).AND.(LSHI(I).EQ.LSHI(J)))) THEN
           MASKI(I)=.FALSE.
           WRITE(TEXT4,'(I4.4)') ISOMIX(J)
           GO TO 195
        ENDIF
  190   CONTINUE
        WRITE(TEXT4,'(I4.4)') ISOMIX(I)
  195   READ(TEXT4,'(A4)') ISONAM(3,I)
      ENDIF
  200 CONTINUE
*
      IF(IMPX.GT.1) THEN
         WRITE (IOUT,320)
         DO 210 I=1,NBISO
           IF(ISOMIX(I).EQ.0) GO TO 210
           DZN=DENMIX(ISOMIX(I))
           IF(DZN.EQ.-1.0) THEN
             WRITE (IOUT,330) I,(ISONAM(I0,I),I0=1,3),(ISONRF(I0,I),
     >       I0=1,3),IHLIB(1,I,1),IHLIB(2,I,1),ILLIB(I),ISOMIX(I),
     >       DENISO(I),TMPISO(I),SNISO(I),LSHI(I),ISHINA(1,I),
     >       ISHINA(2,I),NTFG(I),IHLIB(1,I,3),IHLIB(2,I,3),
     >       IHLIB(1,I,4),IHLIB(2,I,4),IHLIB(1,I,2),IHLIB(2,I,2)
           ELSE
             WRITE (IOUT,340) I,(ISONAM(I0,I),I0=1,3),(ISONRF(I0,I),
     >       I0=1,3),IHLIB(1,I,1),IHLIB(2,I,1),ILLIB(I),ISOMIX(I),
     >       DZN,DENISO(I),TMPISO(I),SNISO(I),LSHI(I),ISHINA(1,I),
     >       ISHINA(2,I),NTFG(I),IHLIB(1,I,3),IHLIB(2,I,3),IHLIB(1,I,4),
     >       IHLIB(2,I,4),IHLIB(1,I,2),IHLIB(2,I,2)
           ENDIF
  210    CONTINUE
      ENDIF
*----
*  READ OLD DILUTIONS IF PRESENT.
*----
      NGIS=NGRO*NBISO
      ALLOCATE(GSN(NGIS),GSB(NGIS))
      CALL XDRSET(GSN,NGIS,1.0E10)
      CALL XDRSET(GSB,NGIS,1.0E10)
      CALL LCMLEN(IPLIB,'ISOTOPESDSN',NELSN,ITYLCM)
      IF(NELSN.GT.0) THEN
        CALL LCMGET(IPLIB,'ISOTOPESDSN',GSN)
        CALL LCMGET(IPLIB,'ISOTOPESDSB',GSB)
      ENDIF
      ILOCSN=0
      ILOCSB=0
      DO 215 ISO=1,NBISO
        IF(SNISO(ISO).GT.0.0) THEN
          CALL XDRSET(GSN(ILOCSN+1),NGRO,SNISO(ISO))
          CALL XDRSET(GSB(ILOCSB+1),NGRO,SBISO(ISO))
        ENDIF
        ILOCSN=ILOCSN+NGRO
        ILOCSB=ILOCSB+NGRO
  215 CONTINUE
*----
*  SAVE THE LIBRARY SPECIFIC INFORMATION.
*----
      NBMIX=0
      DO 220 I=1,NBISO
      NBMIX=MAX(NBMIX,ISOMIX(I))
  220 CONTINUE
      IF(NBMIX.GT.MAXMIX) CALL XABORT('LIBINP: MAXMIX TOO SMALL.')
      TEXT12='L_LIBRARY'
      CALL LCMPTC(IPLIB,'SIGNATURE',12,1,TEXT12)
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=MAXMIX
      ISTATE(2)=NBISO
      ISTATE(3)=NGRO
      ISTATE(4)=NL
      ISTATE(5)=ITRANC
      ISTATE(6)=IPROB
      ISTATE(7)=ITIME
      ISTATE(8)=NLIB
      ISTATE(9)=MIN(NGF,NGRO+1)
      ISTATE(10)=IGRMAX
      ISTATE(11)=NDEPL
      ISTATE(12)=NCOMB
      ISTATE(13)=NEDMAC
      ISTATE(14)=NBMIX
      ISTATE(15)=NRES
      ISTATE(17)=IPROC
      ISTATE(18)=IMAC
      ISTATE(19)=NDEL
      ISTATE(20)=0
      ISTATE(21)=ISOADD
      ISTATE(22)=MAXISM
      ISTATE(23)=IPRECI
      CALL LCMPUT(IPLIB,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPLIB,'ISOTOPESUSED',3*NBISO,3,ISONAM)
      CALL LCMPUT(IPLIB,'ISOTOPERNAME',3*NBISO,3,ISONRF)
      CALL LCMPUT(IPLIB,'ISOTOPESMIX',NBISO,1,ISOMIX)
      CALL LCMPUT(IPLIB,'ISOTOPESTODO',NBISO,1,IEVOL)
      CALL LCMPUT(IPLIB,'ISOTOPESTYPE',NBISO,1,ITYP)
      IF(NLIB.GT.0) THEN
        CALL LCMPUT(IPLIB,'ILIBRARYTYPE',2*NBISO,3,IHLIB(1,1,1))
        CALL LCMPUT(IPLIB,'ILIBRARYINDX',NBISO,1,ILLIB)
        CALL LCMPUT(IPLIB,'ISOTOPESCOH',2*NBISO,3,IHLIB(1,1,2))
        CALL LCMPUT(IPLIB,'ISOTOPESINC',2*NBISO,3,IHLIB(1,1,3))
        CALL LCMPUT(IPLIB,'ISOTOPESRESK',2*NBISO,3,IHLIB(1,1,4))
        CALL LCMPUT(IPLIB,'ISOTOPESNTFG',NBISO,1,NTFG)
        CALL LCMPUT(IPLIB,'ISOTOPESHIN',3*NBISO,3,ISHINA)
        CALL LCMPUT(IPLIB,'ISOTOPESSHI',NBISO,1,LSHI)
        CALL LCMPUT(IPLIB,'ISOTOPESDSN',NGIS,2,GSN)
        CALL LCMPUT(IPLIB,'ISOTOPESDSB',NGIS,2,GSB)
        CALL LCMPUT(IPLIB,'ISOTOPESGIR',NBISO,2,GIR)
        CALL LCMPUT(IPLIB,'ISOTOPESNIR',NBISO,1,NIR)
      ENDIF
      CALL LCMPUT(IPLIB,'ISOTOPESTEMP',NBISO,2,TMPISO)
      IF(NEDMAC.GT.0) THEN
         CALL LCMPTC(IPLIB,'ADDXSNAME-P0',8,NEDMAC,HVECT)
      ENDIF
      CALL LCMPUT(IPLIB,'MIXTUREGAS',NBMIX,1,KGAS)
      DEALLOCATE(GSB,GSN)
*----
*  CHECK FOR DUPLICATE ALIAS.
*----
      DO 255 I=1,NBISO
      DO 250 J=I+1,NBISO
      IF((ISOMIX(I).EQ.ISOMIX(J)).AND.(ISONRF(1,I).EQ.ISONRF(1,J))
     1   .AND.(ISONRF(2,I).EQ.ISONRF(2,J))
     2   .AND.(ISONRF(3,I).EQ.ISONRF(3,J)).AND.(LSHI(I).NE.0)) THEN
         WRITE(HSMG,390) (ISONAM(I1,I),I1=1,3),(ISONAM(I1,J),I1=1,3),
     >   (ISONRF(I1,I),I1=1,3),ISOMIX(I)
         CALL XABORT(HSMG)
      ENDIF
  250 CONTINUE
  255 CONTINUE
*----
*  READ AND INTERPOLATE IN THE MICROSCOPIC X-SECTIONS LIBRARIES.
*----
      IF(NGRO.EQ.0) CALL XABORT('LIBINP: NUMBER OF GROUPS NOT DEFINED.')
      IF((IPROC.EQ.0).AND.(NLIB.GT.0)) THEN
*        ------------------------------------
         CALL LIBLIB (IPLIB,NBISO,MASKI,IMPX)
*        ------------------------------------
      ELSE IF((IPROC.GT.0).AND.(NLIB.GT.0)) THEN
         CALL LIBSUB (MAXISO,IPLIB,IPROC,NGRO,NBISO,NLIB,ISONAM,
     1   TMPISO,MASKI,IPRECI,IMPX)
      ELSE IF((IPROC.EQ.-1).AND.(NLIB.GT.0)) THEN
         JPLIB=LCMLID(IPLIB,'ISOTOPESLIST',NBISO)
      ENDIF
      CALL LCMVAL(IPLIB,' ')
*
      IF(IMPX.GT.0) THEN
         CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
         ITRANC=ISTATE(5)
         NGF=ISTATE(9)
         IGRMAX=ISTATE(10)
         NDEPL=ISTATE(11)
         NBESP=ISTATE(16)
         NDEL=ISTATE(19)
         NFISS=ISTATE(20)
         NPART=ISTATE(26)
         WRITE (IOUT,300) IMPX,IPROB,ITIME,NLIB,NGF,IGRMAX,NBISO,NBMIX,
     1   NRES,NCOMB,NEDMAC,NGRO,NL
         WRITE (IOUT,305) ITRANC,NBESP,IPROC,IMAC,NDEL,NDEPL,NFISS,
     1   NPART,ISOADD,MAXISM,IPRECI
         IF(NEDMAC.GT.0) WRITE (IOUT,310) (I,HVECT(I),I=1,NEDMAC)
         IF(NLIB.GT.0) THEN
            ALLOCATE(INAME(16*NLIB))
            CALL LCMGET(IPLIB,'ILIBRARYNAME',INAME)
            WRITE(IOUT,315)
            DO 260 ILIB=1,NLIB
            WRITE(NAMFIL,'(16A4)') (INAME(16*(ILIB-1)+I),I=1,16)
            WRITE(IOUT,'(1X,I4,4H -- ,A)') ILIB,NAMFIL
  260       CONTINUE
            DEALLOCATE(INAME)
         ENDIF
      ENDIF
*----
*  COMPUTE AND STORE THE EFFECTIVE DENSITY FROM AWR AND MATERIAL DENSITY
*----
      DO 270 IMX=1,NBMIX
      IF(MASK(IMX).AND.(DENMIX(IMX).GE.0.0)) THEN
         CALL LIBCON(IPLIB,IMX,NBISO,ISOMIX,DENISO,DENMIX(IMX),1)
      ENDIF
  270 CONTINUE
      IF(IMPX.GT.0) THEN
         WRITE (IOUT,370)
         DO 280 I=1,NBISO
           IF(ISOMIX(I).EQ.0) GO TO 280
           IF(MASK(ISOMIX(I))) THEN
              WRITE (IOUT,380) I,(ISONAM(I0,I),I0=1,3),(ISONRF(I0,I),
     >        I0=1,3),IHLIB(1,I,1),IHLIB(2,I,1),ILLIB(I),ISOMIX(I),
     >        DENISO(I),TMPISO(I),SNISO(I),LSHI(I),ISHINA(1,I),
     >        ISHINA(2,I),NTFG(I),IHLIB(1,I,3),IHLIB(2,I,3),
     >        IHLIB(1,I,4),IHLIB(2,I,4),IHLIB(1,I,2),IHLIB(2,I,2)
           ENDIF
  280    CONTINUE
      ENDIF
      CALL LCMPUT(IPLIB,'ISOTOPESDENS',NBISO,2,DENISO)
*----
*  STORE MIXTURES DENSITIES
*----
      CALL LCMPUT(IPLIB,'MIXTURESDENS',NBMIX,2,DENMIX)
*----
*  COMPUTE THE MACROSCOPIC X-SECTIONS.
*----
      IF(IMAC.EQ.1) THEN
         ALLOCATE(MASKL(NGRO))
         CALL XDLSET(MASKL,NGRO,.TRUE.)
         ITSTMP=0
         TMPDAY(1)=0.0
         TMPDAY(2)=0.0
         TMPDAY(3)=0.0
         CALL LIBMIX(IPLIB,NBMIX,NGRO,NBISO,ISONAM,ISOMIX,DENISO,MASK,
     >   MASKL,ITSTMP,TMPDAY)
         DEALLOCATE(MASKL)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(MASKI,MASK)
      DEALLOCATE(TMPMIX,GIR,SBISO,SNISO,TMPISO,DENMIX,DENISO)
      DEALLOCATE(KGAS,ITYP,IEVOL,ILLIB,IHLIB,NIR,LSHI,NTFG,ISHINA,
     > ISOMIX,ISONRF,ISONAM)
      RETURN
*
  300 FORMAT(/8H OPTIONS/8H -------/
     1 7H IMPX  ,I6,30H   (0=NO PRINT/1=SHORT/2=MORE)/
     2 7H IPROB ,I6,23H   (0=DIRECT/1=ADJOINT)/
     3 7H ITIME ,I6,28H   (1=STEADY-STATE/2=PROMPT)/
     4 7H NLIB  ,I6,32H   (NUMBER OF SETS OF LIBRARIES)/
     5 7H NGF   ,I6,48H   (NUMBER OF FAST GROUP WITHOUT SELF-SHIELDING)/
     6 7H IGRMAX,I6,41H   (LAST GROUP INDEX WITH SELF-SHIELDING)/
     7 7H NBISO ,I6,36H   (NUMBER OF ISOTOPES OR MATERIALS)/
     8 7H NBMIX ,I6,23H   (NUMBER OF MIXTURES)/
     9 7H NRES  ,I6,40H   (NUMBER OF SETS OF RESONANT MIXTURES)/
     1 7H NCOMB ,I6,33H   (NUMBER OF DEPLETING MIXTURES)/
     2 7H NEDMAC,I6,34H   (NUMBER OF CROSS SECTION EDITS)/
     3 7H NGRO  ,I6,28H   (NUMBER OF ENERGY GROUPS)/
     4 7H NL    ,I6,30H   (NUMBER OF LEGENDRE ORDERS))
  305 FORMAT(
     1 7H ITRANC,I6,45H   (0=NO TRANSPORT CORRECTION/1=APOLLO TYPE/2,
     2 57H=RECOVER FROM LIBRARY/3=WIMS-D TYPE/4=LEAKAGE CORRECTION)/
     3 7H NBESP ,I6,47H   (NUMBER OF ENERGY-DEPENDENT FISSION SPECTRA)/
     4 7H IPROC ,I6,48H   (-1=SKIP LIBRARY PROCESSING/0=DILUTION INTERP,
     5 48HOLATION/1=USE PHYSICAL TABLES/2=BUILD A DRAGLIB/,
     6 61H3=COMPUTE CALENDF TABLES/4=SLOWING-DOWN TABLES/5=ALL CALENDF)/
     7 7H IMAC  ,I6,45H   (0=DO NOT/1=DO BUILD AN EMBEDDED MACROLIB)/
     8 7H NDEL  ,I6,31H   (NUMBER OF PRECURSOR GROUPS)/
     9 7H NDEPL ,I6,33H   (NUMBER OF DEPLETING ISOTOPES)/
     1 7H NFISS ,I6,48H   (NUMBER OF FISSILE ISOTOPES WITH PYIELD DATA)/
     2 7H NPART ,I6,34H   (NUMBER OF COMPANION PARTICLES)/
     3 7H ISOADD,I6,37H   (0=COMPLETE BURNUP CHAIN/1=DO NOT)/
     4 7H MAXISM,I6,40H   (MAX. NUMBER OF ISOTOPES PER MIXTURE)/
     5 7H IPRECI,I6,34H   (CALENDF ACCURACY FLAG:1/2/3/4))
  310 FORMAT(/45H CROSS SECTION EDIT NAME (LCM DIRECTORY NAME)/1X,
     1 44(1H-)/(1X,I3,2X,A6,5X,I3,2X,A6,5X,I3,2X,A6))
  315 FORMAT(/35H AVAILABLE CROSS-SECTION LIBRARIES:)
  320 FORMAT(/' SPEC   LOCAL NAME    ISOTOPE       FROM LIBRARY  MI',
     1 'X   DENSITY     WEIGHT%     TEMP(K)    SIGZERO   SELF-SHIEL  T',
     2 'HERMAL CORRECTION'/' -----  ------------  ------------  ------',
     3 '------  ----  ----------  ----------  ---------  --------  ',
     4 '----------  ------------------')
  330 FORMAT(1X,I5,2X,3A4,2X,3A4,2X,2A4,I4,2X,I4,1P,E12.4,12X,E11.3,
     1 E10.2,I4,2X,2A4,I4,1X,6A4)
  340 FORMAT(1X,I5,2X,3A4,2X,3A4,2X,2A4,I4,2X,I4,1P,2E12.4,E11.3,E10.2,
     1 I4,2X,2A4,I4,1X,6A4)
  370 FORMAT(/56X,'NUMBER'/' SPEC   LOCAL NAME    ISOTOPE       FROM L',
     1 'IBRARY  MIX   DENSITY     TEMP(K)    SIGZERO    SELF-SHIEL  ',
     2 'THERMAL CORRECTION'/' -----  ------------  ------------  -----',
     3 '-------  ----  ----------  ---------  ---------  ---------- ',
     4 ' ------------------')
  380 FORMAT(1X,I5,2X,3A4,2X,3A4,2X,2A4,I4,2X,I4,1P,E12.4,2E11.3,I4,2X,
     1 2A4,I4,1X,6A4)
  390 FORMAT(9HLIBINP: ',3A4,7H' AND ',3A4,24H' ARE BOTH ALIAS FOR THE,
     1 23H SAME LIBRARY ISOTOPE ',3A4,12H' IN MIXTURE,I5,1H.)
  400 FORMAT(' Error in LIBINP : Invalid group structure',2I10)
      END
