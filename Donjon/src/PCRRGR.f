*DECK PCRRGR
      SUBROUTINE PCRRGR(IPMAP,LCUBIC,NMIX,IMPX,NCAL,NCH,NB,NFUEL,
     1 NPARM,ITER,MAXNIS,TERP,NISO,HISO,CONC,LMIXC,XS_CALC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute TERP factors for PMAXS file interpolation. Use global and
* local parameters from a fuel-map object and optional user-defined
* values.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPMAP   address of the fuel-map object.
* LCUBIC  =.TRUE.: cubic Ceschino interpolation; =.FALSE: linear
*         Lagrange interpolation.
* NMIX    number of material mixtures in the fuel-map macrolib.
* IMPX    printing index (=0 for no print).
* NCAL    number of elementary calculations in the PMAXS file.
* NCH     number of reactor channels.
* NB      number of fuel bundles per channel.
* NFUEL   number of fuel types.
* NPARM   number of additional parameters (other than burnup) defined
*         in FMAP object
*
*Parameters: output
* ITER    completion flag (=0: all over; =1: use another PMAXS file;
*         =2 use another L_MAP + PMAXS file).
* MAXNIS  maximum value of NISO(I) in user data.
* TERP    interpolation factors.
* NISO    number of user-selected isotopes.
* HISO    name of the user-selected isotopes.
* CONC    user-defined number density of the user-selected isotopes.
* LMIXC   flag set to .true. for fuel-map mixtures to process.
* XS_CALC pointers towards PMAXS elementary calculations.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE PCRDATA
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER, PARAMETER::MAXISD=400
      TYPE(C_PTR) IPMAP
      INTEGER NMIX,IMPX,NCAL,NFUEL,NCH,NB,ITER,MAXNIS,NPARM,
     1 HISO(2,NMIX,MAXISD),NISO(NMIX)
      REAL TERP(NCAL,NMIX),CONC(NMIX,MAXISD)
      LOGICAL LCUBIC,LMIXC(NMIX)
      TYPE(XSBLOCK_ITEM) XS_CALC(NCAL)
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::MAXADD=10
      INTEGER, PARAMETER::MAXLIN=50
      INTEGER, PARAMETER::MAXPAR=50
      INTEGER, PARAMETER::MAXVAL=200
      REAL, PARAMETER::REPS=1.0E-4
      REAL BURN0, BURN1, FLOTT, SUM, VALR1, VALR2, VARVAL
      INTEGER I0, IBM, IBTYP, IB, ICAL, ICH, IFUEL, ILONG, IMIX,
     & IMPY, INDIC, IPAR, ISO, ITYPE, ITYP, IVARTY, I, JBM, JB,
     & JCAL, JPARM, JPAR, J, NCOMLI, NISOMI, NITMA, NPARMP, NPAR,
     & NTOT, N, IBRA, IBSET, NBURN, IND, II, INDELT
      CHARACTER TEXT12*12,PARKEY(MAXPAR)*12,HSMG*131,RECNAM*12,
     1 COMMEN(MAXLIN)*80,PARNAM*12,HCUBIC*12
      INTEGER NVALUE(MAXPAR),MUPLET(MAXPAR),MUTYPE(MAXPAR),
     1 MAPLET(MAXPAR,MAXADD),MATYPE(MAXPAR,MAXADD),IDLTA(MAXPAR,MAXADD),
     2 NDLTA(MAXPAR),IDLTA1,MUPLT2(MAXPAR),MUTYP2(MAXPAR),
     3 HISOMI(2,MAXISD)
      DOUBLE PRECISION DFLOTT
      REAL VALR(MAXPAR,2),VREAL(MAXVAL,MAXPAR),CONCMI(MAXISD),
     1 VALRA(MAXPAR,2,MAXADD)
      LOGICAL LDELT(MAXPAR),LDELT1,LSET(MAXPAR),LADD(MAXPAR),
     1 LSET1,LADD1,LDMAP(MAXPAR,2),LAMAP(MAXPAR,2,MAXADD),
     2 LCUB2(MAXPAR),LTST
      TYPE(C_PTR) JPMAP,KPMAP
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: FMIX,ZONEC
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ZONEDP,MUBASE
      REAL, ALLOCATABLE, DIMENSION(:) :: BRN0,BRN1,VARC,TERPA
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WPAR
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LPARM,LDELTA
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: HPAR
*----
*  SCRATCH STORAGE ALLOCATION
*   FMIX    fuel mixture indices per fuel bundle.
*   BRN0    contains either low burnup integration limits or
*           instantaneous burnups per fuel bundle.
*   BRN1    upper burnup integration limits per fuel bundle.
*   WPAR    other parameter distributions.
*   HPAR    'PARKEY' name of the other parameters.
*----
      ALLOCATE(LPARM(NPARM+1),FMIX(NCH*NB),ZONEDP(NCH,NB),
     1 ZONEC(NCH),BRN0(NCH*NB),BRN1(NCH*NB),WPAR(NCH*NB,NPARM),
     2 LDELTA(NMIX),HPAR(NPARM+1))
*----
*  RECOVER TABLE-OF-CONTENT INFORMATION FOR THE PMAXS FILE. THE I-TH
*  PMAXS FILE INFORMATION CORRESPONDS TO POINTERS bran_i and PMAX.
*----
      NPAR=bran_i%Nstat_var
      NVALUE(:NPAR)=0
      DO IPAR=1,bran_i%Nstat_var
        PARKEY(IPAR)=bran_i%var_nam(IPAR)
      ENDDO
      IF(PMAX%NBset.GT.0) THEN
        NPAR=NPAR+1
        PARKEY(NPAR)='B'
        NVALUE(NPAR)=PMAX%Bset(1)%NBURN
        VREAL(:,NPAR)=REAL(PMAX%Bset(1)%burns(:))
        VREAL(:,NPAR)=REAL(PMAX%Bset(1)%burns(:))*1000.0
      ENDIF
      IF(NPAR.GT.MAXPAR) CALL XABORT('PCRRGR: MAXPAR OVERFLOW.')
      IF(NHST.NE.1) CALL XABORT('PCRRGR: MULTIPLE HISTORY CASE NOT IMP'
     1 //'LEMENTED.')
      NCOMLI=6
      COMMEN(:6)=hcomment(:6)
      DO IBRA=1,NBRA
        DO IPAR=1,bran_i%Nstat_var
          FLOTT=REAL(bran_i%state(IPAR,IBRA))
          IF(PARKEY(IPAR).EQ.'TF') FLOTT=(FLOTT**2)-273.15
          IF(NVALUE(IPAR).EQ.0) THEN
            NVALUE(IPAR)=1
            VREAL(1,IPAR)=FLOTT
          ELSE
            DO I=1,NVALUE(IPAR)
              IF(FLOTT.EQ.VREAL(I,IPAR)) THEN
                GO TO 10
              ELSE IF(FLOTT.LT.VREAL(I,IPAR)) THEN
                DO J=NVALUE(IPAR),I,-1
                  VREAL(J+1,IPAR)=VREAL(J,IPAR)
                ENDDO
                VREAL(I,IPAR)=FLOTT
                NVALUE(IPAR)=NVALUE(IPAR)+1
                GO TO 10
              ENDIF
            ENDDO
            IF(FLOTT.GT.VREAL(NVALUE(IPAR),IPAR)) THEN
              NVALUE(IPAR)=NVALUE(IPAR)+1
              VREAL(NVALUE(IPAR),IPAR)=FLOTT
            ENDIF
          ENDIF
   10     CONTINUE
        ENDDO
      ENDDO
      IF((IMPX.GT.0).AND.(bran_i%Nstat_var.GT.0))THEN
        DO IPAR=1,NPAR
          WRITE(RECNAM,'(''pval'',I8.8)') IPAR
          WRITE(IOUT,'(13H PCRRGR: KEY=,A,18H TABULATED POINTS=,
     1    1P,6E12.4/(43X,6E12.4))') PARKEY(IPAR),(VREAL(I,IPAR),I=1,
     2    NVALUE(IPAR))
        ENDDO
      ENDIF
*----
*  PRINT PMAXS FILE AND FUELMAP STATISTICS
*----
      IF(IMPX.GT.0) THEN
        WRITE(IOUT,'(43H PCRRGR: NUMBER OF CALCULATIONS IN PMAXS FI,
     1  3HLE=,I6)') NCAL
        WRITE(IOUT,'(43H PCRRGR: NUMBER OF MATERIAL MIXTURES IN FUE,
     1  6HL MAP=,I6)') NMIX
        WRITE(IOUT,'(43H PCRRGR: NUMBER OF LOCAL VARIABLES INCLUDIN,
     1  9HG BURNUP=,I6)') NPAR
        WRITE(IOUT,'(28H PCRRGR: PMAXS FILE COMMENTS,60(1H-))')
        WRITE(IOUT,'(1X,A)') (COMMEN(I),I=1,NCOMLI)
        WRITE(IOUT,'(9H PCRRGR: ,79(1H-))')
      ENDIF
*----
*  SCAN THE PMAXS FILE INFORMATION TO RECOVER THE MUPLET DATABASE
*----
      IF(IMPX.GT.-5) THEN
        WRITE(IOUT,'(24H PCRRGR: MUPLET DATABASE/12H CALCULATION,4X,
     1  6HMUPLET)')
        WRITE(IOUT,'(16X,20A4)') PARKEY(:NPAR)
      ENDIF
      ALLOCATE(MUBASE(NPAR,NCAL))
      ICAL=0
      DO IBRA=1,NBRA
        INDELT=0
        DO IPAR=1,NPAR
          IF(bran_i%state_nam(IBRA).EQ.PARKEY(IPAR)) THEN
            INDELT=IPAR
            CYCLE
          ENDIF
        ENDDO
        IBSET=PMAX%BRANCH(IBRA,1)%IBSET
        NBURN=PMAX%Bset(IBSET)%NBURN
        DO IPAR=1,bran_i%Nstat_var
          FLOTT=REAL(bran_i%state(IPAR,IBRA))
          IF(PARKEY(IPAR).EQ.'TF') FLOTT=(FLOTT**2)-273.15
          IND=0
          DO I=1,NVALUE(IPAR)
            IF(FLOTT.EQ.VREAL(I,IPAR)) THEN
               IND=I
               EXIT
            ENDIF
          ENDDO
          IF(IND.EQ.0) THEN
            CALL XABORT('PCRRGR: MUPLET ALGORITHM FAILURE.')
          ELSE
            MUPLET(IPAR)=IND
          ENDIF
        ENDDO
        IF((NBURN.EQ.PMAX%Bset(1)%NBURN).OR.(NBURN.EQ.1)) THEN
          DO I=1,NBURN
            MUPLET(bran_i%Nstat_var+1)=I
            II=ICAL+I
            MUBASE(:bran_i%Nstat_var+1,II)=MUPLET(:bran_i%Nstat_var+1)
            XS_CALC(ICAL+I)%IBURN=I
            XS_CALC(ICAL+I)%XS=>PMAX%BRANCH(IBRA,1)%XS(I)
            XS_CALC(ICAL+I)%TIV=>PMAX%TIVB(1)%TIV(I)
            IF(INDELT.GT.0) THEN
              XS_CALC(ICAL+I)%DELTA=bran_i%state(INDELT,IBRA)-
     1        bran_i%state(INDELT,1)
            ELSE
              XS_CALC(ICAL+I)%DELTA=0.0
            ENDIF
          ENDDO
        ELSE
          CALL XABORT('PCRRGR: INVALID VALUE OF NBURN.')
        ENDIF
        IF(IMPX.GT.-5) THEN
          DO I=ICAL+1,ICAL+NBURN
            WRITE(IOUT,'(I8,2X,A2,2X,20I4/(14X,20I4))') I,
     1      bran_i%state_nam(IBRA),MUBASE(:NPAR,I)
          ENDDO
        ENDIF
        ICAL=ICAL+NBURN
      ENDDO !IBRA
      IF(ICAL.NE.NCAL) CALL XABORT('PCRRGR: MUPLET ALGORITHM FAILURE.')
*----
*  READ (INTERP_DATA) AND SET VALR PARAMETERS CORRESPONDING TO THE
*  INTERPOLATION POINT. FILL MUPLET FOR PARAMETERS SET WITHOUT
*  INTERPOLATION.
*----
      IBM=0
      MAXNIS=0
      NISOMI=0
      LDELT1=.FALSE.
      LADD1=.FALSE.
      NISO(:NMIX)=0
      LDELTA(:NMIX)=.FALSE.
      IDLTA1=0
      DO I=1,MAXPAR
         LSET(I)=.FALSE.
         LDELT(I)=.FALSE.
         LADD(I)=.FALSE.
         LDMAP(I,:2)=.FALSE.
         LAMAP(I,:2,:MAXADD)=.FALSE.
         NDLTA(I)=0
      ENDDO
      TERP(:NCAL,:NMIX)=0.0
      LMIXC(:NMIX)=.FALSE.
*----
*  ADD THE PARKEY NAME OF THE BURNUP FOR THIS PMAX FILE.
*----
      NPARMP=NPARM+1
      HPAR(NPARMP)='B'
*----
*  MAIN LOOP OF THE SUBROUTINE (UNTIL THE END)
*----
   20 CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF(INDIC.NE.3)CALL XABORT('PCRRGR: CHARACTER DATA EXPECTED(2).')
   30 IF(TEXT12.EQ.'MIX')THEN
        NISOMI=0
        IVARTY=0
        IBTYP=0
        MUPLET(:NPAR)=0
        MUTYPE(:NPAR)=0
        VALR(:NPAR,1)=0.0
        VALR(:NPAR,2)=0.0
        DO 35 I=1,MAXADD
          MAPLET(:NPAR,I)=0
          MATYPE(:NPAR,I)=0
          VALRA(:NPAR,1,I)=0.0
          VALRA(:NPAR,2,I)=0.0
   35   CONTINUE
        DO I=1,MAXPAR
          LSET(I)=.FALSE.
          LDELT(I)=.FALSE.
          LADD(I)=.FALSE.
          LDMAP(I,:2)=.FALSE.
          LAMAP(I,:2,:MAXADD)=.FALSE.
        ENDDO
        LCUB2(:NPAR)=LCUBIC
        CALL REDGET(INDIC,IBM,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.1)CALL XABORT('PCRRGR: INTEGER DATA EXPECTED.')
*       CHECK FUEL MIXTURE
        JPMAP=LCMGID(IPMAP,'FUEL')
        DO IFUEL=1,NFUEL
          KPMAP=LCMGIL(JPMAP,IFUEL)
          CALL LCMGET(KPMAP,'MIX',IMIX)
          IF(IMIX.EQ.IBM)GOTO 50
        ENDDO
        WRITE(IOUT,*)'PCRRGR: UNABLE TO FIND FUEL MIXTURE ',IBM
        CALL XABORT('PCRRGR: WRONG MIXTURE NUMBER.')
   50   CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.3)CALL XABORT('PCRRGR: CHARACTER DATA EXPECTED(3).')
        GOTO 30
      ELSEIF(TEXT12.EQ.'MICRO')THEN
        IF(IBM.EQ.0) CALL XABORT('PCRRGR: MIX NOT SET (1).')
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.3)CALL XABORT('PCRRGR: CHARACTER DATA EXPECTED(5).')
   60   IF(TEXT12.EQ.'ENDMIX')THEN
          GOTO 30
        ELSE
          NISOMI=NISOMI+1
          IF(NISOMI.GT.MAXISD) CALL XABORT('PCRRGR: MAXISD OVERFLOW.')
          MAXNIS=MAX(MAXNIS,NISOMI)
          READ(TEXT12,'(2A4)') (HISOMI(I0,NISOMI),I0=1,2)
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          IF(INDIC.EQ.2)THEN
            CONCMI(NISOMI)=FLOTT
          ELSE
            CALL XABORT('PCRRGR: INVALID HISO DATA.')
          ENDIF
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          IF(INDIC.NE.3)CALL XABORT('PCRRGR: CHARACTER DATA EXPECTED.')
          GOTO 60
        ENDIF
      ELSEIF((TEXT12.EQ.'SET').OR.(TEXT12.EQ.'DELTA').OR.
     1       (TEXT12.EQ.'ADD'))THEN
        IF(IBM.EQ.0) CALL XABORT('PCRRGR: MIX NOT SET (2).')
        ITYPE=0
        LSET1=.FALSE.
        LDELT1=.FALSE.
        LADD1=.FALSE.
        IF(TEXT12.EQ.'SET')THEN
          ITYPE=1
          LSET1=.TRUE.
        ELSEIF(TEXT12.EQ.'DELTA')THEN
          ITYPE=2
          LDELT1=.TRUE.
        ELSEIF(TEXT12.EQ.'ADD')THEN
          ITYPE=2
          LADD1=.TRUE.
          IDLTA1=IDLTA1+1
          DO 65 JPAR=1,NPAR
            MAPLET(JPAR,IDLTA1)=MUPLET(JPAR)
            MATYPE(JPAR,IDLTA1)=MUTYPE(JPAR)
   65     CONTINUE
        ENDIF
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.3)CALL XABORT('PCRRGR: CHARACTER DATA EXPECTED(7).')
        IF((TEXT12.EQ.'LINEAR').OR.(TEXT12.EQ.'CUBIC')) THEN
          HCUBIC=TEXT12
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        ELSE
          HCUBIC=' '
        ENDIF
        IF(INDIC.NE.3)CALL XABORT('PCRRGR: CHARACTER DATA EXPECTED(8).')
        IPAR=-99
        DO I=1,NPAR
          IF(TEXT12.EQ.PARKEY(I))THEN
            IPAR=I
            PARNAM=TEXT12
            GOTO 70
          ENDIF
        ENDDO
        WRITE(HSMG,'(18HPCRRGR: PARAMETER ,A,14H NOT FOUND(1).)') TEXT12
        CALL XABORT(HSMG)
*
   70   IF(HCUBIC.EQ.'LINEAR') THEN
          LCUB2(IPAR)=.FALSE.
        ELSE IF(HCUBIC.EQ.'CUBIC') THEN
          LCUB2(IPAR)=.TRUE.
        ENDIF
        IF((IPAR.GT.NPAR).OR.(IPAR.LE.NPAR))THEN
          CALL REDGET(INDIC,NITMA,VALR1,TEXT12,DFLOTT)
          IF(INDIC.EQ.2)THEN
            VALR2=VALR1
            IF(LSET1) THEN
              LSET(IPAR)=.TRUE.
              VALR(IPAR,1)=VALR1
              VALR(IPAR,2)=VALR1
            ENDIF
            IF(LDELT1) THEN
              LDELT(IPAR)=.TRUE.
              VALR(IPAR,1)=VALR1
              VALR(IPAR,2)=VALR1
            ELSEIF(LADD1) THEN
              LADD(IPAR)=.TRUE.
              VALRA(IPAR,1,IDLTA1)=VALR1
              VALRA(IPAR,2,IDLTA1)=VALR1
              NDLTA(IPAR)=NDLTA(IPAR)+1
              IF(NDLTA(IPAR).GT.MAXADD) CALL XABORT('PCRRGR: MAXADD OV'
     1        //'ERFLOW.')
              IDLTA(IPAR,NDLTA(IPAR))=IDLTA1
            ENDIF
          ELSEIF(TEXT12.EQ.'MAP')THEN
            IF(LDELT1)THEN
              LDELT(IPAR)=.TRUE.
              LDMAP(IPAR,1)=.TRUE.
            ELSEIF(LADD1)THEN
              LADD(IPAR)=.TRUE.
              NDLTA(IPAR)=NDLTA(IPAR)+1
              IF(NDLTA(IPAR).GT.MAXADD) CALL XABORT('PCRRGR: MAXADD OV'
     1        //'ERFLOW.')
              LAMAP(IPAR,1,NDLTA(IPAR))=.TRUE.
              IDLTA(IPAR,NDLTA(IPAR))=IDLTA1
            ENDIF
            IF(LSET1.AND.(.NOT.LSET(IPAR))) GO TO 20
          ELSE
            CALL XABORT('PCRRGR: real value or "MAP" expected(1).')
          ENDIF
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          IF(ITYPE.GE.2)THEN
            IF(INDIC.EQ.2)THEN
              VALR2=FLOTT
              IF(LDELT1)THEN
                VALR(IPAR,2)=VALR2
              ELSEIF(LADD1)THEN
                VALRA(IPAR,2,IDLTA1)=VALR2
              ENDIF
            ELSEIF(TEXT12.EQ.'MAP')THEN
              IF(LDELT1)THEN
                LDMAP(IPAR,2)=.TRUE.
              ELSEIF(LADD1)THEN
                LAMAP(IPAR,2,IDLTA1)=.TRUE.
              ENDIF
            ELSE
              CALL XABORT('PCRRGR: real value or "MAP" expected(2).')
            ENDIF
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          ENDIF
          LTST=.FALSE.
          IF(.NOT.LADD1)THEN
            IF(VALR(IPAR,1).EQ.VALR(IPAR,2)) LTST=.TRUE.
            MUPLET(IPAR)=-1
            MUTYPE(IPAR)=ITYPE
          ELSE
            MAPLET(IPAR,IDLTA1)=-1
            MATYPE(IPAR,IDLTA1)=2
          ENDIF
          IF((LTST).AND.(ITYPE.EQ.1))THEN
            DO J=1,NVALUE(IPAR)
              IF(ABS(VALR(IPAR,1)-VREAL(J,IPAR)).LE.REPS*
     1        ABS(VREAL(J,IPAR)))THEN
                MUPLET(IPAR)=J
                GOTO 30
              ENDIF
            ENDDO
          ENDIF
*----
*  ERRORS HANDLING
*----
          IF(VALR1.LT.VREAL(1,IPAR))THEN
*           OUTSIDE OF THE DOMAIN (1)
            WRITE(HSMG,'(23HPCRRGR: REAL PARAMETER ,A,10H WITH VALU,
     1       1HE,1P,E12.4,26H IS OUTSIDE THE DOMAIN.(1))') PARNAM,VALR1
             WRITE(6,*)'Domain:',VREAL(1,IPAR),' <-> ',
     1       VREAL(NVALUE(IPAR),IPAR)
            CALL XABORT(HSMG)
          ELSEIF(VALR2.GT.VREAL(NVALUE(IPAR),IPAR))THEN
*           OUTSIDE OF THE DOMAIN (2)
            WRITE(HSMG,'(23HPCRRGR: REAL PARAMETER ,A,10H WITH VALU,
     1       1HE,1P,E12.4,26H IS OUTSIDE THE DOMAIN.(2))') PARNAM,VALR2
             WRITE(6,*)'Domain:',VREAL(1,IPAR),' <-> ',
     1       VREAL(NVALUE(IPAR),IPAR)
            CALL XABORT(HSMG)
          ELSEIF((VALR1.GT.VALR2).AND.(ITYPE.EQ.1))THEN
*           ITYPE=1 correspond to an integral between VALR1 and VALR2
*           otherwise it is a simple difference
            WRITE(HSMG,'(23HPCRRGR: REAL PARAMETER ,A,9H IS DEFIN,
     1       7HED WITH,1P,E12.4,2H >,E12.4,4H.(1))') PARNAM,
     2       VALR1,VALR2
            CALL XABORT(HSMG)
          ENDIF
          IF((LADD1).AND.(TEXT12.EQ.'REF'))THEN
  120       IPAR=-99
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(TEXT12.EQ.'ENDREF') GOTO 140
            DO I=1,NPAR
              IF(TEXT12.EQ.PARKEY(I))THEN
                IPAR=I
                GOTO 130
              ENDIF
            ENDDO
            CALL XABORT('PCRRGR: PARAMETER '//TEXT12//' NOT FOUND(2).')
  130       CONTINUE
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.EQ.2)THEN
              VALRA(IPAR,1,IDLTA1)=FLOTT
              VALRA(IPAR,2,IDLTA1)=FLOTT
              MAPLET(IPAR,IDLTA1)=-1
              MATYPE(IPAR,IDLTA1)=1
              DO J=1,NVALUE(IPAR)
                IF(ABS(VALRA(IPAR,1,IDLTA1)-VREAL(J,IPAR)).LE.
     1             REPS*ABS(VREAL(J,IPAR)))THEN
                  MAPLET(IPAR,IDLTA1)=J
                  GOTO 120
                ENDIF
              ENDDO
            ELSEIF(TEXT12.EQ.'SAMEASREF')THEN
              MAPLET(IPAR,IDLTA1)=-1
              MATYPE(IPAR,IDLTA1)=-1
            ELSE
              CALL XABORT('PCRRGR: REAL or "SAMEASREF" expected')
            ENDIF
            GOTO 120
  140       CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          ELSE IF((LDELT1).AND.(TEXT12.EQ.'REF'))THEN
  150       IPAR=-99
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(TEXT12.EQ.'ENDREF') GOTO 170
            DO I=1,NPAR
              IF(TEXT12.EQ.PARKEY(I))THEN
                IPAR=I
                GOTO 160
              ENDIF
            ENDDO
            CALL XABORT('PCRRGR: PARAMETER '//TEXT12//' NOT FOUND(3).')
  160       CONTINUE
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.EQ.2)THEN
              VALR(IPAR,1)=FLOTT
              VALR(IPAR,2)=FLOTT
              MUPLET(IPAR)=-1
              MUTYPE(IPAR)=1
              DO J=1,NVALUE(IPAR)
                IF(ABS(VALR(IPAR,1)-VREAL(J,IPAR)).LE.REPS*
     1          ABS(VREAL(J,IPAR)))THEN
                  MUPLET(IPAR)=J
                  GOTO 150
                ENDIF
              ENDDO
            ELSEIF(TEXT12.EQ.'SAMEASREF')THEN
              MUPLET(IPAR)=-1
              MUTYPE(IPAR)=-1
            ELSE
              CALL XABORT('PCRRGR: REAL or "SAMEASREF" expected')
            ENDIF
            GOTO 150
  170       CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          ENDIF
          GOTO 30
        ENDIF
      ELSEIF(TEXT12.EQ.'TIMAV-BURN')THEN
        IF(IBM.EQ.0) CALL XABORT('PCRRGR: MIX NOT SET (3).')
        IBTYP=1
      ELSEIF(TEXT12.EQ.'INST-BURN')THEN
        IF(IBM.EQ.0) CALL XABORT('PCRRGR: MIX NOT SET (4).')
        IBTYP=2
      ELSEIF(TEXT12.EQ.'AVG-EX-BURN')THEN
        IF(IBM.EQ.0) CALL XABORT('PCRRGR: MIX NOT SET (5).')
        IBTYP=3
        CALL REDGET(INDIC,IVARTY,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.1)CALL XABORT('PCRRGR: INTEGER DATA EXPECTED.')
      ELSEIF(TEXT12.EQ.'ENDMIX')THEN
*----
*  RECOVER FUEL-MAP INFORMATION.
*----
        IF(IMPX.GT.0) THEN
          DO IPAR=1,NPAR
            IF(LCUB2(IPAR)) THEN
              WRITE(IOUT,'(26H PCRRGR: GLOBAL PARAMETER:,A12,5H ->CU,
     1        18HBIC INTERPOLATION.)') PARKEY(IPAR)
            ELSE
              WRITE(IOUT,'(26H PCRRGR: GLOBAL PARAMETER:,A12,5H ->LI,
     1        19HNEAR INTERPOLATION.)') PARKEY(IPAR)
            ENDIF
          ENDDO
        ENDIF
        FMIX(:NCH*NB)=0
        CALL LCMGET(IPMAP,'FLMIX',FMIX)
        CALL NCRMAP(IPMAP,NPARM,HPAR,NCH,NB,IBTYP,IMPX,BRN0,BRN1,WPAR,
     1              LPARM)
        IF(IBTYP.EQ.3) THEN
          IF(IVARTY.EQ.0) CALL XABORT('PCRRGR: IVARTY NOT SET.')
          CALL LCMGET(IPMAP,'B-ZONE',ZONEC)
          DO ICH=1,NCH
            DO J=1,NB
              IF(ZONEC(ICH).EQ.IVARTY) THEN
                ZONEDP(ICH,J)=1
              ELSE
                ZONEDP(ICH,J)=0
              ENDIF
            ENDDO
          ENDDO
          CALL LCMLEN(IPMAP,'B-VALUE',ILONG,ITYP)
          IF (ILONG.EQ.0) CALL XABORT('PCRRGR: NO SAVED VALUES FOR '
     1    //'THIS TYPE OF VARIABLE IN L_MAP')
          ALLOCATE(VARC(ILONG))
          CALL LCMGET(IPMAP,'B-VALUE',VARC)
          VARVAL=VARC(IVARTY)
          DEALLOCATE(VARC)
        ENDIF
*----
*  PERFORM INTERPOLATION OVER THE FUEL MAP.
*----
        DO 185 JPARM=1,NPARMP
          IPAR=-99
          DO I=1,NPAR
            IF(HPAR(JPARM).EQ.PARKEY(I))THEN
              IPAR=I
              IF(LSET(IPAR)) THEN
                WRITE(6,*) 'L_MAP values overwritten by the SET option'
     1          // ' for parameter '//HPAR(JPARM)
                IF(.NOT.LADD(IPAR)) LPARM(JPARM)=.FALSE.
              ENDIF
              GOTO 185
            ENDIF
          ENDDO
          LPARM(JPARM)=.FALSE.
  185   CONTINUE
*----
*  COMPUTE ALL THE MUPLETS FOR EACH BUNDLE
*----
        IMPY=MAX(0,IMPX-1)
        NTOT=0
        DO 285 JB=1,NB
        DO 280 ICH=1,NCH
        IB=(JB-1)*NCH+ICH
        IF(FMIX(IB).EQ.0) GO TO 280
        NTOT=NTOT+1
        IPAR=-99
        IF(FMIX(IB).EQ.IBM)THEN
          IF(NTOT.GT.NMIX) CALL XABORT('PCRRGR: NMIX OVERFLOW.')
          DO 260 JPARM=1,NPARMP
          IF(.NOT.LPARM(JPARM))GOTO 260
        DO I=1,NPAR
          IF(HPAR(JPARM).EQ.PARKEY(I))THEN
            IPAR=I
            PARNAM=HPAR(JPARM)
            GOTO 190
          ENDIF
        ENDDO
        WRITE(HSMG,'(18HPCRRGR: PARAMETER ,A,14H NOT FOUND(4).)')
     1  HPAR(JPARM)
        CALL XABORT(HSMG)
  190   CONTINUE
        ITYPE=0
        IF((JPARM.EQ.NPARMP).AND.(NPARMP.EQ.NPARM+1))THEN
*         parameter JPARAM is burnup
          IF(.NOT.LSET(IPAR))THEN
            MUTYPE(IPAR)=1
            MUPLET(IPAR)=-1
            BURN0=0.0
            BURN1=0.0
            IF(IBTYP.EQ.1)THEN
*             TIME-AVERAGE
              BURN0=BRN0(IB)
              BURN1=BRN1(IB)
            ELSEIF(IBTYP.EQ.2)THEN
*             INSTANTANEOUS
              BURN0=BRN0(IB)
              BURN1=BURN0
            ELSEIF(IBTYP.EQ.3)THEN
*             DIFFERENCIATION RELATIVE TO EXIT BURNUP
              ITYPE=3
              BURN0=BRN0(IB)
              BURN1=BRN1(IB)
            ENDIF
            VALR(IPAR,1)=BURN0
            VALR(IPAR,2)=BURN1
            VALR1=VALR(IPAR,1)
            VALR2=VALR(IPAR,2)
            ITYPE=1
          ENDIF
        ELSE
          IF(.NOT.LSET(IPAR))THEN
            VALR(IPAR,1)=WPAR(IB,JPARM)
            VALR(IPAR,2)=WPAR(IB,JPARM)
            MUPLET(IPAR)=-1
            MUTYPE(IPAR)=1
            VALR1=VALR(IPAR,1)
            VALR2=VALR(IPAR,2)
            ITYPE=1
          ENDIF
          IF(LDMAP(IPAR,1).OR.LDMAP(IPAR,2))THEN
            IF(LDMAP(IPAR,1)) VALR(IPAR,1)=WPAR(IB,JPARM)
            IF(LDMAP(IPAR,2)) VALR(IPAR,2)=WPAR(IB,JPARM)
            MUPLET(IPAR)=-1
            MUTYPE(IPAR)=2
            VALR1=VALR(IPAR,1)
            VALR2=VALR(IPAR,2)
            ITYPE=2
          ELSE IF(LADD(IPAR))THEN
            DO N=1,NDLTA(IPAR)
              IDLTA1=IDLTA(IPAR,N)
              IF(LAMAP(IPAR,1,IDLTA1)) THEN
                VALRA(IPAR,1,IDLTA1)=WPAR(IB,JPARM)
                MAPLET(IPAR,IDLTA1)=-1
                MATYPE(IPAR,IDLTA1)=2
              ENDIF
              IF(LAMAP(IPAR,2,IDLTA1)) THEN
                VALRA(IPAR,2,IDLTA1)=WPAR(IB,JPARM)
                MAPLET(IPAR,IDLTA1)=-1
                MATYPE(IPAR,IDLTA1)=2
              ENDIF
            ENDDO
            VALR1=VALRA(IPAR,1,IDLTA(IPAR,1))
            VALR2=VALRA(IPAR,2,IDLTA(IPAR,1))
            ITYPE=2
          ENDIF
        ENDIF
        IF(ITYPE.EQ.1)THEN
          IF(VALR1.EQ.VALR2)THEN
            DO J=1,NVALUE(IPAR)
             IF(ABS(VALR1-VREAL(J,IPAR)).LE.REPS*ABS(VREAL(J,IPAR)))THEN
               MUPLET(IPAR)=J
               MUTYPE(IPAR)=ITYPE
               GOTO 260
             ENDIF
            ENDDO
          ENDIF
        ENDIF
*----
*  ERRORS HANDLING
*----
        IF(VALR1.LT.VREAL(1,IPAR))THEN
*         OUTSIDE OF THE DOMAIN (1)
          WRITE(HSMG,'(23HPCRRGR: REAL PARAMETER ,A,10H WITH VALU,
     1    1HE,1P,E12.4,26H IS OUTSIDE THE DOMAIN(3).)') PARNAM,VALR1
          WRITE(6,*)'Domain:',VREAL(1,IPAR),' <-> ',
     1    VREAL(NVALUE(IPAR),IPAR)
          CALL XABORT(HSMG)
        ELSEIF(VALR2.GT.VREAL(NVALUE(IPAR),IPAR))THEN
*         OUTSIDE OF THE DOMAIN (2)
          WRITE(HSMG,'(23HPCRRGR: REAL PARAMETER ,A,10H WITH VALU,
     1    1HE,1P,E12.4,26H IS OUTSIDE THE DOMAIN(4).)') PARNAM,VALR2
          WRITE(6,*)'Domain:',VREAL(1,IPAR),' <-> ',
     1    VREAL(NVALUE(IPAR),IPAR)
          CALL XABORT(HSMG)
        ELSEIF((ITYPE.EQ.1).AND.(VALR1.GT.VALR2))THEN
*         VALR1 > VALR2
          WRITE(HSMG,'(23HPCRRGR: REAL PARAMETER ,A,9H IS DEFIN,
     1    7HED WITH,1P,E12.4,2H >,E12.4,4H.(2))') PARNAM,
     2    VALR1,VALR2
          CALL XABORT(HSMG)
        ENDIF
*----
*  COMPUTE THE TERP FACTORS USING TABLE-OF-CONTENT INFORMATION.
*----
  260     CONTINUE
          LMIXC(NTOT)=.TRUE.
          IF(IMPY.GT.2) WRITE(6,'(32H PCRRGR: COMPUTE TERP FACTORS IN,
     1    17H FUEL-MAP MIXTURE,I5,1H.)') NTOT
          NISO(NTOT)=NISOMI
          LDELTA(NTOT)=LDELT1
          DO ISO=1,NISOMI
            HISO(1,NTOT,ISO)=HISOMI(1,ISO)
            HISO(2,NTOT,ISO)=HISOMI(2,ISO)
            CONC(NTOT,ISO)=CONCMI(ISO)
          ENDDO
          DO JPAR=1,NPAR
            MUPLT2(JPAR)=MUPLET(JPAR)
          ENDDO
          IF(IBTYP.EQ.3)THEN
             IF(ZONEDP(ICH,JB).NE.0) THEN
                CALL PCRTRP(LCUB2,IMPY,NPAR,NCAL,NVALUE,MUPLT2,
     1                      MUTYPE,VALR,VARVAL,MUBASE,VREAL,
     2                      TERP(1,NTOT))
             ELSE
                TERP(:NCAL,NTOT)=0.0
             ENDIF
          ELSE
             CALL PCRTRP(LCUB2,IMPY,NPAR,NCAL,NVALUE,MUPLT2,
     1                   MUTYPE,VALR,VARVAL,MUBASE,VREAL,
     2                   TERP(1,NTOT))
          ENDIF
*         DELTA-ADD
          DO 270 IPAR=1,NPAR
            IF(LADD(IPAR))THEN
              DO N=1,NDLTA(IPAR)
                IDLTA1=IDLTA(IPAR,N)
                DO JPAR=1,NPAR
                  MUPLT2(JPAR)=MAPLET(JPAR,IDLTA1)
                  MUTYP2(JPAR)=MATYPE(JPAR,IDLTA1)
                ENDDO
                DO JPAR=1,NPAR
                  IF(MUTYP2(JPAR).LT.0)THEN
                    MUPLT2(JPAR)=MUPLET(JPAR)
                    MUTYP2(JPAR)=MUTYPE(JPAR)
                    VALRA(JPAR,1,IDLTA1)=VALR(JPAR,1)
                    VALRA(JPAR,2,IDLTA1)=VALR(JPAR,2)
                  ENDIF
                ENDDO
                ALLOCATE(TERPA(NCAL))
                CALL PCRTRP(LCUB2,IMPY,NPAR,NCAL,NVALUE,MUPLT2,
     1          MUTYP2,VALRA(1,1,IDLTA1),VARVAL,MUBASE,VREAL,
     2          TERPA)
                DO 275 JCAL=1,NCAL
                TERP(JCAL,NTOT)=TERP(JCAL,NTOT)+TERPA(JCAL)
  275           CONTINUE
                DEALLOCATE(TERPA)
              ENDDO
            ENDIF
  270     CONTINUE
        ENDIF
  280   CONTINUE
  285   CONTINUE
        IF(NTOT.GT.NMIX) CALL XABORT('PCRRGR: ALGORITHM FAILURE.')
        IBM=0
      ELSEIF((TEXT12.EQ.'PMAXS').OR.(TEXT12.EQ.'TABLE').OR.
     1  (TEXT12.EQ.'CHAIN').OR.(TEXT12.EQ.';')) THEN
*----
*  CHECK TERP FACTORS AND RETURN
*----
        IF(TEXT12.EQ.';') ITER=0
        IF(TEXT12.EQ.'PMAXS') ITER=1
        IF(TEXT12.EQ.'TABLE') ITER=2
        IF(TEXT12.EQ.'CHAIN') ITER=3
        DO 300 IBM=1,NMIX
        IF(.NOT.LMIXC(IBM)) GO TO 300
        IF(NISO(IBM).GT.MAXNIS) CALL XABORT('PCRRGR: MAXNIS OVERFLOW.')
        IF(LDELTA(IBM)) THEN
          SUM=0.0
        ELSE
          SUM=1.0
        ENDIF
        DO 290 ICAL=1,NCAL
        SUM=SUM-TERP(ICAL,IBM)
  290   CONTINUE
        IF(ABS(SUM).GT.1.0E-4) THEN
           WRITE(HSMG,'(43HPCRRGR: INVALID INTERPOLATION FACTORS IN MI,
     1     5HXTURE,I4,1H.)') IBM
           CALL XABORT(HSMG)
        ENDIF
  300   CONTINUE
*----
*  EXIT MAIN LOOP OF THE SUBROUTINE
*----
        GO TO 310
      ELSE
        CALL XABORT('PCRRGR: '//TEXT12//' IS AN INVALID KEYWORD.')
      ENDIF
      GOTO 20
*----
*  PRINT INTERPOLATION (TERP) FACTORS
*----
  310 IF(IMPX.GT.2) THEN
        WRITE(IOUT,'(/30H PCRRGR: INTERPOLATION FACTORS)')
        DO ICAL=1,NCAL
          DO IBM=1,NMIX
            IF(TERP(ICAL,IBM).NE.0.0) THEN
              WRITE(IOUT,320) ICAL,(TERP(ICAL,JBM),JBM=1,NMIX)
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(MUBASE)
      DEALLOCATE(HPAR,LDELTA,WPAR,BRN1,BRN0,ZONEC,ZONEDP,FMIX,LPARM)
      RETURN
  320 FORMAT(6H CALC=,I8,6H TERP=,1P,8E13.5/(20X,8E13.5))
      END
