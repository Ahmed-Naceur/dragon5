*DECK EDIGET
      SUBROUTINE EDIGET(IPEDIT,IFGEO,NGROUP,NGCOND,NREG,NBMIX,MATCOD,
     >                  ITMERG,NMERGE,IHF,IFFAC,ILUPS,NSAVES,NSTATS,
     >                  IGCR,EGCR,IMERGE,CURNAM,OLDNAM,IADF,NW,ICURR,
     >                  NBMICR,CARISO,NACTI,IACTI,IPRINT,MAXPTS,ICALL,
     >                  ISOTXS,LISO,IADJ,MACGEO,IEUR,NOUT,HVOUT,BB2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read edition option parameters.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPEDIT  pointer to the edition LCM object.
* IFGEO   unit file number of the surfacic file.
* NGROUP  number of groups.
* NREG    number of regions.
* NBMIX   maximum number of mixtures.
* MATCOD  mixture index in region.
*
*Parameters: output
* NGCOND  number of groups condensed.
* ITMERG  type of technique to compute merge indices:
*         = 0 no merge;
*         =-1 merge by geometry (equigeom);
*         =-2 merge by cell;
*         =-3 merge by HMIX defined in GEO:;
*         =-4 merge using IMERGE array directly.
* NMERGE  number of merged indices in array IMERGE.
* IHF     H-factor calculation (= 0 no; =1 yes).
* IFFAC   four factor calculation flag:
*         = 0 no four factors (default);
*         = 1 four factor evaluation.
* ILUPS   remove up-scattering from output.
* NSAVES  homogenized x-s computation+save:
*         = 0 no compute no save;
*         = 1 compute, no save;
*         = 2 compute and save.
* NSTATS  statistics level:
*         = 0 no statistics;
*         = 1 statistics on fluxes;
*         = 2 statistics on reaction rates;
*         = 3 statistics on fluxes and reaction rates;
*         =-1 delta sigma (MERG COMP only).
* IGCR    condensed group limits.
* EGCR    condensed energy limits.
* IMERGE  merged region positions.
* CURNAM  name of LCM directory where the current rates are to be
*         stored.
* OLDNAM  name of LCM directory where old reaction rates were stored.
* IADF    flag for computing boundary or ADF information:
*         = 0 do not compute them;
*         = 1 compute boundary currents using ALBS information;
*         = 2 recover averaged fluxes in boundary regions;
*         = -2 compute ADF using averaged fluxes in boundary regions;
*         = 3 compute boundary information using SYBIL/ARM or MOC
*         interface currents;
*         = 4 recover ADF information from input macrolib.
* NW      type of weighting for P1 cross section info:
*         =0 use flux to merge/condense P1 matrices;
*         =1 use current to merge/condense P1 matrices.
* ICURR   type of current approximation if NW=1:
*         =1: heterogeneous leakage;
*         =2: Todorova outscatter approximation;
*         =4: use spherical harmonic moments of the flux.
* NBMICR  type of microlib edition:
*         =-2: process only macroscopic residue;
*         =-1: process each isotope;
*         =0: process no isotope;
*         >0 number of isotopes to process.
* CARISO  names of the isotopes to process.
* NACTI   number of activation edit.
* IACTI   activation mixtures.
* IPRINT  print index.
* MAXPTS  maximum number of macro-regions.
* ICALL   maximum directory index in IPEDIT.
* ISOTXS  ISOTX file enabling flag (0: off; 1: binary; 2: ascii).
* LISO    =.TRUE. if we want to keep all the isotopes after 
*         homogeneization.
* IADJ    type of flux weighting:
*         =0: direct flux weighting;
*         =1: direct-adjoint flux weighting.
* MACGEO  name of the macro-geometry.
* IEUR    type of tracking tone on the macro-geometry:
*         =1: SYBIL or EXCELL;
*         =2: NXT;
*         =3: else.
* NOUT    number of output cross section types (set to zero to recover
*         all cross section types).
* HVOUT   MATXS names of the output cross section types.
* BB2     imposed leakage used in non-regression tests.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE EDIG2S_MOD
*----
*  SUBROUTINE ARGUMENTS
*----
      PARAMETER    (MAXED=100,MAXOUT=100)
      TYPE(C_PTR)   IPEDIT
      INTEGER       IFGEO,NGROUP,NGCOND,NREG,NBMIX,MATCOD(NREG),ITMERG,
     >              NMERGE,IHF,IFFAC,ILUPS,NSAVES,NSTATS,IGCR(NGROUP),
     >              IMERGE(NREG),IADF,NW,ICURR,NBMICR,NACTI,
     >              IACTI(NBMIX),IPRINT,MAXPTS,ICALL,ISOTXS,IADJ,
     >              IEUR,NOUT
      REAL          EGCR(NGROUP),BB2
      LOGICAL       LISO
      CHARACTER     CURNAM*12,OLDNAM*12,CARISO(MAXED)*12,MACGEO*12,
     >              HVOUT(MAXOUT)*8,HSMG*131
*----
*  LOCAL VARIABLES
*----
      CHARACTER     CARLIR*8,HTYPE*8
      REAL          REALIR
      DOUBLE PRECISION DBLLIR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIXMER,INADF,IOFGAP,IREMIX
      CHARACTER*8, ALLOCATABLE, DIMENSION(:) :: HADF
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MIXMER(0:NBMIX))
*----
*  INITIALIZE MIXMER
*----
      DO 10 IMATER=0,NBMIX
        MIXMER(IMATER)=IMATER
   10 CONTINUE
*----
*  READ OPTION NAME
*----
      ISOTXS=0
   20 CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
   30 IF(ITYPLU.EQ.10) GO TO 250
      IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTER VARI'
     > //'ABLE EXPECTED')
   40 IF(CARLIR.EQ.';') THEN
        GO TO 250
      ELSE IF(CARLIR.EQ.'EDIT') THEN
        CALL REDGET(ITYPLU,IPRINT,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('EDIGET: READ ERROR - INTEGER VARI'
     >  //'ABLE EXPECTED')
      ELSE IF(CARLIR.EQ.'NADF') THEN
        IADF=0
      ELSE IF(CARLIR.EQ.'ALBS') THEN
        IADF=1
      ELSE IF(CARLIR.EQ.'ADF') THEN
        IADF=2
        CALL REDGET(ITYPLU,INTLIR,REALIR,HTYPE,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTER*8 '
     >  //'TYPE EXPECTED(1)')
        IF(HTYPE.EQ.'*') THEN
          IADF=-2
          CALL REDGET(ITYPLU,INTLIR,REALIR,HTYPE,DBLLIR)
          IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTER*'
     >    //'8 TYPE EXPECTED(2)')
        ENDIF
        CALL LCMSIX(IPEDIT,'REF:ADF',1)
        CALL LCMLEN(IPEDIT,'NTYPE',ILONG,ITYLCM)
        IF(ILONG.EQ.0) THEN
          NTYPE=0
        ELSE
          CALL LCMGET(IPEDIT,'NTYPE',NTYPE)
        ENDIF
        ALLOCATE(INADF(NTYPE+1),HADF(NTYPE+1),IOFGAP(NREG))
        IF(NTYPE.GT.0) THEN
          CALL LCMGET(IPEDIT,'NADF',INADF)
          CALL LCMGTC(IPEDIT,'HADF',8,NTYPE,HADF)
        ENDIF
        CALL XDISET(IOFGAP,NREG,0)
        IGAP=0
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTER*4 '
     >  //'TYPE EXPECTED')
        IF(CARLIR(:4).EQ.'REGI') THEN
   50     CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.EQ.1) THEN
            IGAP=IGAP+1
            IF(IGAP.GT.NREG) THEN
              CALL XABORT('EDIGET: BOUNDARY REGI OVERFLOW(1)')
            ELSE IF(INTLIR.GT.NREG) THEN
              CALL XABORT('EDIGET: BOUNDARY REGO OVERFLOW(2)')
            ELSE IF(IOFGAP(IGAP).NE.0) THEN
              CALL XABORT('EDIGET: REGI ALREADY DEFINED')
            ENDIF
            IOFGAP(IGAP)=INTLIR
          ELSE IF((ITYPLU.EQ.3).AND.(CARLIR.EQ.'ENDR')) THEN
            GO TO 80
          ELSE
            CALL XABORT('EDIGET: INTEGER OR ENDR KEYWORD EXPECTED')
          ENDIF
          GO TO 50
        ELSE IF(CARLIR.EQ.'MIX') THEN
   60     CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.EQ.1) THEN
            DO 70 IREG=1,NREG
            IF(MATCOD(IREG).EQ.INTLIR) THEN
              IGAP=IGAP+1
              IF(IGAP.GT.NREG) THEN
                CALL XABORT('EDIGET: BOUNDARY MIX OVERFLOW(1)')
              ELSE IF(INTLIR.GT.NBMIX) THEN
                CALL XABORT('EDIGET: BOUNDARY MIX OVERFLOW(2)')
              ELSE IF(IOFGAP(IGAP).NE.0) THEN
                CALL XABORT('EDIGET: MIX ALREADY DEFINED')
              ENDIF
              IOFGAP(IGAP)=IREG
            ENDIF
   70       CONTINUE
          ELSE IF((ITYPLU.EQ.3).AND.(CARLIR.EQ.'ENDM')) THEN
            GO TO 80
          ELSE
            CALL XABORT('EDIGET: INTEGER OR ENDM KEYWORD EXPECTED')
          ENDIF
          GO TO 60
        ELSE
          CALL XABORT('EDIGET: REGI OR MIX KEYWORD EXPECTED(1)')
        ENDIF
   80   NTYPE=NTYPE+1
        INADF(NTYPE)=IGAP
        HADF(NTYPE)=HTYPE
*
        CALL LCMPUT(IPEDIT,'NTYPE',1,1,NTYPE)
        CALL LCMPUT(IPEDIT,'NADF',NTYPE,1,INADF)
        CALL LCMPTC(IPEDIT,'HADF',8,NTYPE,HADF)
        CALL LCMPUT(IPEDIT,HTYPE,IGAP,1,IOFGAP)
        CALL LCMSIX(IPEDIT,' ',2)
*
        DEALLOCATE(IOFGAP,HADF,INADF)
      ELSE IF(CARLIR.EQ.'JOUT') THEN
        IADF=3
      ELSE IF(CARLIR.EQ.'ADFM') THEN
        IADF=4
      ELSE IF(CARLIR(:4).EQ.'MGEO') THEN
         CALL REDGET(ITYPLU,INTLIR,REALIR,MACGEO,DBLLIR)
         IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTER'
     >   //' VARIABLE EXPECTED')
      ELSE IF(CARLIR.EQ.'UPS') THEN
        ILUPS=1
      ELSE IF(CARLIR.EQ.'P0W') THEN
*       FLUX WEIGHTING OF THE PN MATRICES.
        NW=0
        ICURR=0
      ELSE IF(CARLIR.EQ.'P1W_L') THEN
*       FUNDAMENTAL CURRENT WEIGHTING OF THE PN MATRICES.
        NW=1
        ICURR=1
      ELSE IF(CARLIR.EQ.'P1W_TO') THEN
*       TODOROVA OUTSCATTER CURRENT WEIGHTING OF THE PN MATRICES.
        NW=1
        ICURR=2
      ELSE IF(CARLIR.EQ.'PNW_SP') THEN
*       SPHERICAL HARMONICS WEIGHTING OF THE PN MATRICES.
        NW=1
        ICURR=4
      ELSE IF(CARLIR(:4).EQ.'MICR') THEN
        CALL REDGET(ITYPLU,NBMICR,REALIR,CARLIR,DBLLIR)
        IF((ITYPLU.EQ.3).AND.(CARLIR(:4).EQ.'ALLX')) THEN
*         TO REGISTER ALL ISOTOPES CROSS SECTION IN THE MERGED REGIONS
          LISO=.TRUE.
          CALL REDGET(ITYPLU,NBMICR,REALIR,CARLIR,DBLLIR)
        ENDIF
        IF((ITYPLU.EQ.3).AND.(CARLIR(:4).EQ.'ISOT')) THEN
          ISOTXS=1
          CALL REDGET(ITYPLU,NBMICR,REALIR,CARLIR,DBLLIR)
          IF((ITYPLU.EQ.3).AND.(CARLIR(:4).EQ.'ASCI')) THEN
            ISOTXS=2
            CALL REDGET(ITYPLU,NBMICR,REALIR,CARLIR,DBLLIR)
          ENDIF
        ENDIF
        IF((ITYPLU.EQ.3).AND.(CARLIR.EQ.'RES')) THEN
          NBMICR=-2
        ELSE IF((ITYPLU.EQ.3).AND.(CARLIR.EQ.'ALL')) THEN
          NBMICR=-1
        ELSE IF(ITYPLU.EQ.1) THEN
          IF(NBMICR.GT.MAXED) CALL XABORT('EDIGET: TOO MANY MICR')
          DO 90 IIII=1,NBMICR
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARISO(IIII),DBLLIR)
            IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTE'
     >      //'R VARIABLE EXPECTED')
   90     CONTINUE
        ELSE
          CALL XABORT('EDIGET: READ ERROR - KEY ISOTXS, ALL, NONE OR I'
     >    //'NTEGER VARIABLE EXPECTED AFTER MICR')
        ENDIF
      ELSE IF(CARLIR(:4).EQ.'REAC') THEN
        CALL REDGET(ITYPLU,NOUT,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('EDIGET: READ ERROR - INTEGER'
     >  //' VARIABLE EXPECTED')
        IF(NOUT.GT.MAXOUT) CALL XABORT('EDIGET: MAXOUT OVERFLOW')
        DO 100 IOT=1,NOUT
        CALL REDGET(ITYPLU,INTLIR,REALIR,HVOUT(IOT),DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTER'
     >  //' VARIABLE EXPECTED')
  100   CONTINUE
      ELSE IF(CARLIR(:4).EQ.'ACTI') THEN
        IF((ITYPLU.EQ.3).AND.(CARLIR(:4).EQ.'ISOT')) THEN
          ISOTXS=1
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF((ITYPLU.EQ.3).AND.(CARLIR(:4).EQ.'ASCI')) THEN
            ISOTXS=2
            CALL REDGET(ITYPLU,NBMICR,REALIR,CARLIR,DBLLIR)
          ENDIF
        ENDIF
        IF((ITYPLU.EQ.3).AND.(CARLIR.EQ.'NONE')) THEN
          NACTI=0
        ELSE
          DO 211 IREG=1,NBMIX
            IF(ITYPLU.EQ.1) THEN
              IF(INTLIR.GT.NBMIX) CALL XABORT('EDIGET: INVALID ACTIVAT'
     >        //'ION INDEX')
              NACTI=NACTI+1
              IACTI(NACTI)=INTLIR
            ELSE
              GO TO 30
            ENDIF
            IF(IREG.LT.NBMIX) THEN
              CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            ENDIF
  211     CONTINUE
        ENDIF
      ELSE IF(CARLIR(:4).EQ.'COND') THEN
*----
*  GROUP CONDENSATION DIRECTIVE ANALYSIS
*----
        DO 108 IGROUP=1,NGROUP+1
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.EQ.3) THEN
            IF(IGROUP.EQ.1) THEN
              IF(CARLIR.EQ.'NONE') THEN
                NGCOND=NGROUP
                DO 107 JGROUP=1,NGROUP
                  IGCR(JGROUP)=JGROUP
  107           CONTINUE
                CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
                GO TO 30
              ELSE
                NGCOND=1
                IGCR(NGCOND)=NGROUP
              ENDIF
            ENDIF
            GO TO 30
          ELSE IF(ITYPLU.EQ.1) THEN
            IF(INTLIR.GT.NGROUP) INTLIR=NGROUP
            IF(NGCOND.GT.0) THEN
              IF(INTLIR.GT.IGCR(NGCOND)) THEN
                NGCOND=NGCOND+1
                IGCR(NGCOND)=INTLIR
              ENDIF
            ELSE
              NGCOND=NGCOND+1
              IGCR(NGCOND)=INTLIR
            ENDIF
          ELSE
            IF(NGCOND.GT.0) THEN
              IF(REALIR.LT.EGCR(NGCOND)) THEN
                NGCOND=NGCOND+1
                EGCR(NGCOND)=REALIR
              ENDIF
            ELSE
              NGCOND=NGCOND+1
              EGCR(NGCOND)=REALIR
            ENDIF
          ENDIF
  108   CONTINUE
      ELSE IF(CARLIR(:4).EQ.'MERG') THEN
*----
*  MERGING DIRECTIVE ANALYSIS
*----
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTER VA'
     >  //'RIABLE EXPECTED')
        IF(CARLIR.EQ.'COMP') THEN
*----
*  COMPLETE MERGE
*----
          CALL XDISET(IMERGE,NREG,1)
          ITMERG=-4
          NMERGE=1
          GO TO 20
        ELSE IF(CARLIR.EQ.'GEO') THEN
*----
*  MERGE BY GEOMETRY
*----
          ITMERG=-1
          NMERGE=0
          GO TO 20
        ELSE IF(CARLIR.EQ.'CELL') THEN
*----
*  CELL-BY-CELL MERGE
*----
          ITMERG=-2
          NMERGE=0
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTER '
     >    //'VARIABLE EXPECTED')
          IF((CARLIR.EQ.'SYBIL').OR.(CARLIR.EQ.'EXCELL')) THEN
            IEUR=1
          ELSE IF(CARLIR.EQ.'NXT') THEN
            IEUR=2
          ELSE IF(CARLIR.EQ.'DEFAULT') THEN
            IEUR=3
          ELSE IF(CARLIR.EQ.'UNFOLD') THEN
            IEUR=4
          ELSE IF(CARLIR.EQ.'REMIX') THEN
            GO TO 105
          ELSE
            IEUR=3
            GO TO 40
          ENDIF
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTER '
     >    //'VARIABLE EXPECTED')
  105     IF(CARLIR.EQ.'REMIX') THEN
*           Data to further homogenize a cell-by-cell homogenization.
  110       CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.EQ.1) THEN
              NMERGE=NMERGE+1
              IF(NMERGE.GT.NREG) CALL XABORT('EDIGET: IMERGE(NREG) OVE'
     >        //'RFLOW')
              IMERGE(NMERGE)=INTLIR
              GO TO 110
            ELSE
              GO TO 40
            ENDIF
          ENDIF
          GO TO 40
        ELSE IF(CARLIR.EQ.'HMIX') THEN
*----
*  MERGE BY HOMOGENIZATION MIXTURES
*----
          ITMERG=-3
          NMERGE=0
          GO TO 20
        ELSE IF(CARLIR.EQ.'MIX') THEN
*----
*  MERGE BY MIXTURES
*----
          ITMERG=-4
          NMIXME=0
          DO 114 IREG=1,NREG
            IBM=MATCOD(IREG)
            IF(IBM.GT.NBMIX) CALL XABORT('EDIGET: NBMIX OVERFLOW.')
            NMIXME=MAX(NMIXME,IBM)
            IMERGE(IREG)=MIXMER(IBM)
  114     CONTINUE
          NMERGE=NMIXME
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.EQ.1) THEN
*----
*  SPECIFY MIXTURES TO BE MERGED
*----
            NMERGE=MAX(0,INTLIR)
            MIXMER(1)=INTLIR
            DO 115 IMATER=2,NMIXME
              CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
              IF(ITYPLU.NE.1) CALL XABORT('EDIGET: READ ERROR - INTEGE'
     >        //'R VARIABLE EXPECTED')
              NMERGE=MAX(NMERGE,INTLIR)
              MIXMER(IMATER)=INTLIR
  115       CONTINUE
            DO 116 IREG=1,NREG
              IMERGE(IREG)=MIXMER(MATCOD(IREG))
  116       CONTINUE
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.NE.3) THEN
              WRITE(HSMG,'(40HEDIGET: READ ERROR - CHARACTER VARIABLE ,
     >        10H EXPECTED.,I5,26H MIXTURE INDICES EXPECTED.)') NMIXME
              CALL XABORT(HSMG)
            ENDIF
            GO TO 40
          ELSE IF(ITYPLU.EQ.3) THEN
*----
*  ASSOCIATE ONE REGION BY MIXTURE
*----
            GO TO 40
          ELSE
            CALL XABORT('EDIGET: READ ERROR - INVALID TYPE READ')
          ENDIF
        ELSE IF(CARLIR(:4).EQ.'REGI') THEN
*----
*  MERGE BY REGIONS
*----
          ITMERG=-4
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.1) CALL XABORT('EDIGET: READ ERROR - INTEGE'
     >    //'R VARIABLE EXPECTED')
          NMERGE=MAX(0,INTLIR)
          IMERGE(1)=INTLIR
          DO 118 IREG=2,NREG
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.NE.1) THEN
               WRITE(CARLIR,'(I4)') NREG
               CALL XABORT('EDIGET: READ ERROR - INTEGER VARIABLE EXPE'
     >         //'CTED NREG='//CARLIR)
            ENDIF
            NMERGE=MAX(NMERGE,INTLIR)
            IMERGE(IREG)=INTLIR
  118     CONTINUE
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.3) THEN
            WRITE(HSMG,'(40HEDIGET: READ ERROR - CHARACTER VARIABLE ,
     >      10H EXPECTED.,I5,25H REGION INDICES EXPECTED.)') NREG
            CALL XABORT(HSMG)
          ENDIF
          GO TO 40
        ELSE IF(CARLIR.EQ.'G2S') THEN
          CALL EDIG2S(IPRINT,IFGEO,NREG,NMERGE,IMERGE)
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTER '
     >    //'VARIABLE EXPECTED')
          IF(CARLIR.EQ.'REMIX') THEN
*           REMIX option.
            NMEOLD=NMERGE
            NMERGE=0
            ALLOCATE(IREMIX(NMEOLD))
            DO II=1,NMEOLD
              CALL REDGET(ITYPLU,IREMIX(II),REALIR,CARLIR,DBLLIR)
              IF(ITYPLU.NE.1) CALL XABORT('EDIGET: READ ERROR - INTEGE'
     >        //'R VARIABLE EXPECTED')
            ENDDO
            DO IREG=1,NREG
              IM=IMERGE(IREG)
              IF(IM.GT.0) THEN
                IF(IM.GT.NMEOLD) CALL XABORT('EDIGET: IMERGE OVERFLOW')
                IMERGE(IREG)=IREMIX(IM)
                NMERGE=MAX(NMERGE,IMERGE(IREG))
              ENDIF
            ENDDO
            DEALLOCATE(IREMIX)
          ELSE
            GO TO 40
          ENDIF
        ELSE IF(CARLIR.EQ.'NONE') THEN
*----
*  NO MERGING
*----
          ITMERG=-4
          NMERGE=NREG
          DO 106 IREG=1,NREG
            IMERGE(IREG)=IREG
  106     CONTINUE
        ELSE
          CALL XABORT('EDIGET: READ ERROR - ILLEGAL KEYWORD '//
     >    'FOLLOWING MERG -- ALLOWED : COMP, MIX REGI, READ : '//
     >     CARLIR)
        ENDIF
      ELSE IF(CARLIR.EQ.'TAKE') THEN
*----
*  TAKE DIRECTIVE ANALYSIS
*----
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTER VA'
     >  //'RIABLE EXPECTED')
        IF(CARLIR.EQ.'MIX') THEN
*----
*  TAKE PER MIXTURE
*----
          NMIXME=0
          DO 120 IREG=1,NREG
            NMIXME=MAX(NMIXME,MATCOD(IREG))
  120     CONTINUE
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.EQ.1) THEN
            CALL XDISET(MIXMER,NMIXME,0)
*----
*  SPECIFY MIXTURES TO BE SELECTED
*----
            IF(INTLIR.LE.NMIXME.AND.INTLIR.GT.0) MIXMER(INTLIR)=1
            NMERGE=1
            DO 122 IMATER=2,NBMIX
              CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
              IF(ITYPLU.NE.1) GO TO 123
              IF(INTLIR.LE.NMIXME.AND.INTLIR.GT.0) MIXMER(INTLIR)=IMATER
              NMERGE=NMERGE+1
  122       CONTINUE
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          ENDIF
  123     CONTINUE
          WRITE(6,'(1X,A6,2X,2I10)') 'MIXMER',NMIXME,NMERGE
          WRITE(6,'(5I10)') (MIXMER(JJJ),JJJ=1,NMIXME)
          DO 124 IREG=1,NREG
            IMERGE(IREG)=MIXMER(MATCOD(IREG))
  124     CONTINUE
          GO TO 30
        ELSE IF(CARLIR(:4).EQ.'REGI') THEN
*----
*  TAKE PER REGIONS
*----
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE.1) CALL XABORT('EDIGET: AT LEAST ONE REGION'
     >      //' MUST BE SELECTED')
          DO 125 IREG=1,NREG
            IMERGE(IREG)=0
  125     CONTINUE
          NMERGE=1
          IMERGE(INTLIR)=1
          DO 126 IREG=2,NREG
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.NE.1) GO TO 30
            NMERGE=NMERGE+1
            IMERGE(INTLIR)=IREG
  126     CONTINUE
        ELSE
          CALL XABORT('EDIGET: READ ERROR - ILLEGAL KEYWORD '//
     >    'FOLLOWING TAKE -- ALLOWED : MIX REGI, READ : '// CARLIR)
        ENDIF
      ELSE IF(CARLIR.EQ.'SAVE') THEN
*----
*  SAVE DIRECTIVE ANALYSIS
*----
        NSAVES=2
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTER VA'
     >  //'RIABLE EXPECTED')
        IF(CARLIR.EQ.'ON') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CURNAM,DBLLIR)
          IF(ITYPLU.EQ.2) CALL XABORT('EDIGET: READ ERROR - REAL VARIA'
     >    //'BLE FORBIDDEN')
          IF(ITYPLU.EQ.1) THEN
            WRITE(CURNAM,'(8HREF-CASE,I4.4)') INTLIR
            ICALL=MAX(ICALL,INTLIR)
          ENDIF
        ELSE
          GO TO 40
        ENDIF
      ELSE IF(CARLIR.EQ.'STAT') THEN
*----
*  STAT DIRECTIVE ANALYSIS
*----
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTER VA'
     >  //'RIABLE EXPECTED')
        IF(CARLIR.EQ.'FLUX') THEN
          NSTATS=1
        ELSE IF(CARLIR.EQ.'RATE') THEN
          NSTATS=2
        ELSE IF(CARLIR.EQ.'ALL ') THEN
          NSTATS=3
        ELSE IF(CARLIR.EQ.'DELS') THEN
          NSTATS=-1
        ELSE
          CALL XABORT('EDIGET: READ ERROR - ILLEGAL KEYWORD '//
     >                 CARLIR)
        ENDIF
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('EDIGET: READ ERROR - CHARACTER VA'
     >  //'RIABLE EXPECTED')
        IF(CARLIR(:4).EQ.'REFE') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,OLDNAM,DBLLIR)
          IF(ITYPLU.EQ.2) CALL XABORT('EDIGET: READ ERROR - REAL VARIA'
     >    //'BLE FORBIDDEN')
          IF(ITYPLU.EQ.1) WRITE(OLDNAM,'(8HREF-CASE,I4.4)') INTLIR
        ELSE
          GO TO 40
        ENDIF
      ELSE IF(CARLIR.EQ.'NOHF') THEN
        IHF=0
      ELSE IF(CARLIR.EQ.'NBAL') THEN
        IFFAC=1000
      ELSE IF(CARLIR.EQ.'MAXR') THEN
        CALL REDGET(ITYPLU,MAXPTS,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('EDIGET: READ ERROR - INTEGER VARI'
     >  //'ABLE EXPECTED')
      ELSE IF(CARLIR(:4).EQ.'DIRE') THEN
        IADJ=0
      ELSE IF(CARLIR(:4).EQ.'PROD') THEN
        IADJ=1
      ELSE IF(CARLIR(:4).EQ.'LEAK') THEN
        CALL REDGET(ITYPLU,INTLIR,BB2,CARLIR,DBLLIR)
        IF(ITYPLU.NE.2) CALL XABORT('EDIGET: REAL DATA EXPECTED.')
      ELSE
        CALL XABORT('EDIGET:ILLEGAL KEYWORD '//CARLIR)
      ENDIF
      GO TO 20
*----
*  RETURN
*----
  250 IF(IPRINT.GE.2) NSAVES=MAX(1,NSAVES)
      IF((NSAVES.EQ.0).AND.((NSTATS.NE.0).OR.(IFFAC.NE.0))) NSAVES=1
      IF((NSAVES.GE.2).AND.(CURNAM.EQ.' ')) THEN
        ICALL=ICALL+1
        WRITE(CURNAM,'(8HREF-CASE,I4.4)') ICALL
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(MIXMER)
      RETURN
      END
