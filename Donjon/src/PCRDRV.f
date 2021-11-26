*DECK PCRDRV
      SUBROUTINE PCRDRV(LCUBIC,NMIX,IMPX,NCAL,ITER,MAXNIS,TERP,NISO,
     1 HISO,CONC,LMIXC,XS_CALC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute TERP factors for PMAXS file interpolation. Use user-defined
* global and local parameters.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* LCUBIC  =.TRUE.: cubic Ceschino interpolation; =.FALSE: linear
*         Lagrange interpolation.
* NMIX    maximum number of material mixtures in the microlib.
* IMPX    print parameter (equal to zero for no print).
* NCAL    number of elementary calculations in the PMAXS file.
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
      INTEGER, PARAMETER::MAXISD=200
      INTEGER NMIX,IMPX,NCAL,ITER,MAXNIS,NISO(NMIX),HISO(2,NMIX,MAXISD)
      REAL TERP(NCAL,NMIX),CONC(NMIX,MAXISD)
      LOGICAL LCUBIC,LMIXC(NMIX)
      TYPE(XSBLOCK_ITEM) XS_CALC(NCAL)
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::MAXLIN=50
      INTEGER, PARAMETER::MAXPAR=50
      INTEGER, PARAMETER::MAXVAL=200
      REAL, PARAMETER::REPS=1.0E-4
      REAL FLOTT, SUM
      INTEGER I0, IBM, ICAL, INDIC, IPAR, ITYPE, I, JBM, J, NCOMLI,
     & NITMA, NPAR, IBRA, IBSET, II, IND, INDELT, NBURN
      CHARACTER TEXT12*12,PARKEY(MAXPAR)*12,HSMG*131,COMMEN(MAXLIN)*80,
     1 RECNAM*12,HCUBIC*12
      INTEGER NVALUE(MAXPAR),MUPLET(MAXPAR),MUTYPE(MAXPAR)
      DOUBLE PRECISION DFLOTT
      REAL VALR(MAXPAR,2),VREAL(MAXVAL,MAXPAR)
      LOGICAL LCUB2(MAXPAR)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: MUBASE
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LDELTA
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(LDELTA(NMIX))
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
      ENDIF
      IF(NPAR.GT.MAXPAR) CALL XABORT('PCRDRV: MAXPAR OVERFLOW.')
      IF(NHST.NE.1) CALL XABORT('PCRDRV: MULTIPLE HISTORY CASE NOT IMP'
     1 //'LEMENTED.')
      NCOMLI=6
      COMMEN(:6)=hcomment(:6)
      DO IBRA=1,NBRA
        DO IPAR=1,bran_i%Nstat_var
          FLOTT=REAL(bran_i%state(IPAR,IBRA))
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
          WRITE(IOUT,'(13H PCRDRV: KEY=,A,18H TABULATED POINTS=,
     1    1P,6E12.4/(43X,6E12.4))') PARKEY(IPAR),(VREAL(I,IPAR),I=1,
     2    NVALUE(IPAR))
        ENDDO
      ENDIF
*----
*  PRINT PMAXS FILE AND FUELMAP STATISTICS
*----
      IF(IMPX.GT.0) THEN
        WRITE(IOUT,'(43H PCRDRV: NUMBER OF CALCULATIONS IN PMAXS FI,
     1  3HLE=,I6)') NCAL
        WRITE(IOUT,'(43H PCRDRV: NUMBER OF MATERIAL MIXTURES IN FUE,
     1  6HL MAP=,I6)') NMIX
        WRITE(IOUT,'(43H PCRDRV: NUMBER OF LOCAL VARIABLES INCLUDIN,
     1  9HG BURNUP=,I6)') NPAR
        WRITE(IOUT,'(28H PCRDRV: PMAXS FILE COMMENTS,60(1H-))')
        WRITE(IOUT,'(1X,A)') (COMMEN(I),I=1,NCOMLI)
        WRITE(IOUT,'(9H PCRDRV: ,79(1H-))')
      ENDIF
*----
*  SCAN THE PMAXS FILE INFORMATION TO RECOVER THE MUPLET DATABASE
*----
      IF(IMPX.GT.-5) THEN
        WRITE(IOUT,'(24H PCRDRV: MUPLET DATABASE/12H CALCULATION,4X,
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
          IND=0
          DO I=1,NVALUE(IPAR)
            IF(FLOTT.EQ.VREAL(I,IPAR)) THEN
               IND=I
               EXIT
            ENDIF
          ENDDO
          IF(IND.EQ.0) THEN
            CALL XABORT('PCRDRV: MUPLET ALGORITHM FAILURE.')
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
          CALL XABORT('PCRDRV: INVALID VALUE OF NBURN.')
        ENDIF
        IF(IMPX.GT.-5) THEN
          DO I=ICAL+1,ICAL+NBURN
            WRITE(IOUT,'(I8,2X,A2,2X,20I4/(14X,20I4))') I,
     1      bran_i%state_nam(IBRA),MUBASE(:NPAR,I)
          ENDDO
        ENDIF
        ICAL=ICAL+NBURN
      ENDDO !IBRA
      IF(ICAL.NE.NCAL) CALL XABORT('PCRDRV: MUPLET ALGORITHM FAILURE.')
*----
*  READ (INTERP_DATA) AND SET VALR PARAMETERS CORRESPONDING TO THE
*  INTERPOLATION POINT. FILL MUPLET FOR PARAMETERS SET WITHOUT
*  INTERPOLATION.
*----
      NISO(:NMIX)=0
      TERP(:NCAL,:NMIX)=0.0
      LMIXC(:NMIX)=.FALSE.
*----
*  READ (INTERP_DATA) AND SET VALR PARAMETERS CORRESPONDING TO THE
*  INTERPOLATION POINT. FILL MUPLET FOR PARAMETERS SET WITHOUT
*  INTERPOLATION.
*----
      IBM=0
      MAXNIS=0
      NISO(:NMIX)=0
      LDELTA(:NMIX)=.FALSE.
   20 CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('PCRDRV: CHARACTER DATA EXPECTED.')
   30 IF(TEXT12.EQ.'MIX') THEN
         MUPLET(:NPAR)=0
         MUTYPE(:NPAR)=0
         VALR(:NPAR,1)=0.0
         VALR(:NPAR,2)=0.0
         LCUB2(:NPAR)=LCUBIC
         CALL REDGET(INDIC,IBM,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('PCRDRV: INTEGER DATA EXPECTED.')
         IF(IBM.GT.NMIX) CALL XABORT('PCRDRV: NMIX OVERFLOW.')
         LMIXC(IBM)=.TRUE.
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('PCRDRV: CHARACTER DATA EXPECTED.')
         GO TO 30
      ELSE IF(TEXT12.EQ.'MICRO') THEN
         IF(IBM.EQ.0) CALL XABORT('PCRDRV: MIX NOT SET (1).')
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('PCRDRV: CHARACTER DATA EXPECTED.')
   50    IF(TEXT12.EQ.'ENDMIX') THEN
            GO TO 30
         ELSE
            NISO(IBM)=NISO(IBM)+1
            IF(NISO(IBM).GT.MAXISD) CALL XABORT('PCRDRV: MAXISD OVERFL'
     1      //'OW.')
            MAXNIS=MAX(MAXNIS,NISO(IBM))
            READ(TEXT12,'(2A4)') (HISO(I0,IBM,NISO(IBM)),I0=1,2)
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.EQ.2) THEN
               CONC(IBM,NISO(IBM))=FLOTT
            ELSE
               CALL XABORT('PCRDRV: INVALID HISO DATA.')
            ENDIF
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.NE.3) CALL XABORT('PCRDRV: CHARACTER DATA EXPECTE'
     1      //'D.')
            GO TO 50
         ENDIF
      ELSE IF((TEXT12.EQ.'SET').OR.(TEXT12.EQ.'DELTA')) THEN
         IF(IBM.EQ.0) CALL XABORT('PCRDRV: MIX NOT SET (2).')
         ITYPE=0
         IF(TEXT12.EQ.'SET') THEN
            ITYPE=1
         ELSE IF(TEXT12.EQ.'DELTA') THEN
            ITYPE=2
            LDELTA(IBM)=.TRUE.
         ENDIF
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('PCRDRV: CHARACTER DATA EXPECTED.')
         IF((TEXT12.EQ.'LINEAR').OR.(TEXT12.EQ.'CUBIC')) THEN
            HCUBIC=TEXT12
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         ELSE
            HCUBIC=' '
         ENDIF
         IF(INDIC.NE.3) CALL XABORT('PCRDRV: CHARACTER DATA EXPECTED.')
         DO 60 I=1,NPAR
         IF(TEXT12.EQ.PARKEY(I)) THEN
            IPAR=I
            GO TO 70
         ENDIF
   60    CONTINUE
         CALL XABORT('PCRDRV: PARAMETER '//TEXT12//' NOT FOUND.')
   70    IF(HCUBIC.EQ.'LINEAR') THEN
            LCUB2(IPAR)=.FALSE.
         ELSE IF(HCUBIC.EQ.'CUBIC') THEN
            LCUB2(IPAR)=.TRUE.
         ENDIF
         CALL REDGET(INDIC,NITMA,VALR(IPAR,1),TEXT12,DFLOTT)
         IF(INDIC.NE.2) CALL XABORT('PCRDRV: REAL DATA EXPECTED.')
         VALR(IPAR,2)=VALR(IPAR,1)
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.EQ.2) THEN
            VALR(IPAR,2)=FLOTT
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         ENDIF
         IF(VALR(IPAR,1).EQ.VALR(IPAR,2)) THEN
            DO 80 J=1,NVALUE(IPAR)
            IF(ABS(VALR(IPAR,1)-VREAL(J,IPAR)).LE.REPS*
     1      ABS(VREAL(J,IPAR)))THEN
               MUPLET(IPAR)=J
               IF(ITYPE.NE.1) MUPLET(IPAR)=-1
               MUTYPE(IPAR)=ITYPE
               GO TO 30
            ENDIF
   80       CONTINUE
         ENDIF
         IF(VALR(IPAR,1).LT.VREAL(1,IPAR)) THEN
            WRITE(HSMG,'(23HPCRDRV: REAL PARAMETER ,A,10H WITH VALU,
     1      1HE,1P,E12.4,25H IS OUTSIDE THE DOMAIN (<,E12.4,1H))')
     2      PARKEY(IPAR),VALR(IPAR,1),VREAL(1,IPAR)
            CALL XABORT(HSMG)
         ELSE IF(VALR(IPAR,2).GT.VREAL(NVALUE(IPAR),IPAR)) THEN
            WRITE(HSMG,'(23HPCRDRV: REAL PARAMETER ,A,10H WITH VALU,
     1      1HE,1P,E12.4,25H IS OUTSIDE THE DOMAIN (>,E12.4,1H))')
     2      PARKEY(IPAR),VALR(IPAR,1),VREAL(NVALUE(IPAR),IPAR)
            CALL XABORT(HSMG)
         ELSE IF(VALR(IPAR,1).GT.VALR(IPAR,2)) THEN
            WRITE(HSMG,'(23HPCRDRV: REAL PARAMETER ,A,9H IS DEFIN,
     1      7HED WITH,1P,E12.4,2H >,E12.4,1H.)') PARKEY(IPAR),
     2      VALR(IPAR,1),VALR(IPAR,2)
            CALL XABORT(HSMG)
         ENDIF
         MUPLET(IPAR)=-1
         MUTYPE(IPAR)=ITYPE
         GO TO 30
      ELSE IF(TEXT12.EQ.'ENDMIX') THEN
*----
*  COMPUTE THE TERP FACTORS USING TABLE-OF-CONTENT INFORMATION.
*----
         IF(IMPX.GT.0) THEN
           DO IPAR=1,NPAR
             IF(LCUB2(IPAR)) THEN
               WRITE(IOUT,'(26H PCRDRV: GLOBAL PARAMETER:,A12,5H ->CU,
     1         18HBIC INTERPOLATION.)') PARKEY(IPAR)
             ELSE
               WRITE(IOUT,'(26H PCRDRV: GLOBAL PARAMETER:,A12,5H ->LI,
     1         19HNEAR INTERPOLATION.)') PARKEY(IPAR)
             ENDIF
           ENDDO
         ENDIF
         IF(IBM.GT.NMIX) CALL XABORT('PCRDRV: MIX OVERFLOW (MICROLIB).')
         IF(NCAL.EQ.1) THEN
           TERP(1,IBM)=1.0
         ELSE
           CALL PCRTRP(LCUB2,IMPX,NPAR,NCAL,NVALUE,MUPLET,MUTYPE,VALR,
     1     0.0,MUBASE,VREAL,TERP(1,IBM))
         ENDIF
         IBM=0
      ELSE IF((TEXT12.EQ.'PMAXS').OR.(TEXT12.EQ.'TABLE').OR.
     1   (TEXT12.EQ.';')) THEN
*----
*  CHECK TERP FACTORS AND RETURN
*----
         IF(TEXT12.EQ.';') ITER=0
         IF(TEXT12.EQ.'PMAXS') ITER=1
         IF(TEXT12.EQ.'TABLE') ITER=2
         DO 150 IBM=1,NMIX
         IF(.NOT.LMIXC(IBM)) GO TO 150
         IF(NISO(IBM).GT.MAXNIS) CALL XABORT('PCRDRV: MAXNIS OVERFLOW.')
         IF(LDELTA(IBM)) THEN
            SUM=0.0
         ELSE
            SUM=1.0
         ENDIF
         DO 140 ICAL=1,NCAL
         SUM=SUM-TERP(ICAL,IBM)
  140    CONTINUE
         IF(ABS(SUM).GT.1.0E-4) THEN
            WRITE(HSMG,'(43HPCRDRV: INVALID INTERPOLATION FACTORS IN MI,
     1      5HXTURE,I4,1H.)') IBM
            CALL XABORT(HSMG)
         ENDIF
  150    CONTINUE
         GO TO 160
      ELSE
         CALL XABORT('PCRDRV: '//TEXT12//' IS AN INVALID KEYWORD.')
      ENDIF
      GO TO 20
*----
*  PRINT INTERPOLATION (TERP) FACTORS
*----
  160 IF(IMPX.GT.2) THEN
        WRITE(IOUT,'(/30H PCRDRV: INTERPOLATION FACTORS)')
        DO ICAL=1,NCAL
          DO IBM=1,NMIX
            IF(TERP(ICAL,IBM).NE.0.0) THEN
              WRITE(IOUT,170) ICAL,(TERP(ICAL,JBM),JBM=1,NMIX)
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(LDELTA)
      RETURN
  170 FORMAT(6H CALC=,I8,6H TERP=,1P,8E13.5/(20X,8E13.5))
      END
