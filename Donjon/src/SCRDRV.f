*DECK SCRDRV
      SUBROUTINE SCRDRV(IPSAP,LCUBIC,NMIX,IMPX,NMIL,NCAL,MD2,ITER,
     1 MAXNIS,MIXC,TERP,NISO,LISO,HISO,CONC,ITODO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute TERP factors for Saphyb interpolation. Use user-defined
* global parameters.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPSAP   address of the Saphyb object.
* LCUBIC  =.TRUE.: cubic Ceschino interpolation; =.FALSE: linear
*         Lagrange interpolation.
* NMIX    maximum number of material mixtures in the microlib.
* IMPX    print parameter (equal to zero for no print).
* NMIL    number of material mixtures in the Saphyb.
* NCAL    number of elementary calculations in the Saphyb.
* MD2     number of particularized and macro isotopes in the Saphyb.
*
*Parameters: output
* ITER    completion flag (=0: all over; =1: use another Saphyb;
*         =2 use another L_MAP + Saphyb).
* MAXNIS  maximum value of NISO(I) in user data.
* MIXC    mixture index in the Saphyb corresponding to each microlib
*         mixture.
* TERP    interpolation factors.
* NISO    number of user-selected isotopes.
* LISO    type of treatment (=.true.: ALL; =.false.: ONLY).
* HISO    name of the user-selected isotopes.
* CONC    user-defined number density of the user-selected isotopes. A
*         value of -99.99 is set to indicate that the compo value is
*         used.
* ITODO   non-depletion mask (=1 to force a user-selected isotope to be
*         non-depleting)
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSAP
      INTEGER NMIX,IMPX,NMIL,NCAL,MD2,ITER,MAXNIS,MIXC(NMIX),
     1 HISO(2,NMIX,MD2),NISO(NMIX),ITODO(NMIX,MD2)
      REAL TERP(NCAL,NMIX),CONC(NMIX,MD2)
      LOGICAL LCUBIC,LISO(NMIX)
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::MAXLIN=50
      INTEGER, PARAMETER::MAXPAR=50
      INTEGER, PARAMETER::MAXVAL=200
      REAL, PARAMETER::REPS=1.0E-4
      INTEGER I, I0, IBM, IBMOLD, ICAL, INDIC, IPAR, ITYLCM, ITYPE, J
     &, JBM, LENGTH, NCOMLI, NITMA, NPAR, NVP
      REAL SUM, FLOTT
      CHARACTER TEXT12*12,PARKEY(MAXPAR)*4,PARFMT(MAXPAR)*8,HSMG*131,
     1 COMMEN(MAXLIN)*80,VALH(MAXPAR)*12,VCHAR(MAXVAL)*12,RECNAM*12,
     2 HCUBIC*12
      INTEGER DIMSAP(50),VALI(MAXPAR),NVALUE(MAXPAR),VINTE(MAXVAL),
     1 MUPLET(2*MAXPAR),MUTYPE(2*MAXPAR)
      DOUBLE PRECISION DFLOTT
      REAL VALR(2*MAXPAR,2),VREAL(MAXVAL)
      LOGICAL LCUB2(MAXPAR)
      TYPE(C_PTR) LPSAP
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LDELTA
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(LDELTA(NMIX))
*----
*  RECOVER TABLE-OF-CONTENT INFORMATION FOR THE SAPHYB.
*----
      CALL LCMGET(IPSAP,'DIMSAP',DIMSAP)
      NCOMLI=DIMSAP(1)
      NPAR=DIMSAP(8)
      NVP=DIMSAP(17)
      IF(NCOMLI.GT.MAXLIN) CALL XABORT('SCRDRV: MAXLIN OVERFLOW.')
      IF(NPAR.GT.MAXPAR) CALL XABORT('SCRDRV: MAXPAR OVERFLOW.')
      CALL LCMGTC(IPSAP,'COMMEN',80,NCOMLI,COMMEN)
      IF(NPAR.GT.0) THEN
         CALL LCMSIX(IPSAP,'paramdescrip',1)
         CALL LCMGTC(IPSAP,'PARKEY',4,NPAR,PARKEY)
         CALL LCMGTC(IPSAP,'PARFMT',8,NPAR,PARFMT)
         CALL LCMSIX(IPSAP,' ',2)
      ENDIF
      IF(IMPX.GT.0) WRITE(IOUT,'(1X,A)') (COMMEN(I),I=1,NCOMLI)
      TERP(:NCAL,:NMIX)=0.0
      MIXC(:NMIX)=0
*----
*  READ (INTERP_DATA) AND SET VALI, VALR AND VALH PARAMETERS
*  CORRESPONDING TO THE INTERPOLATION POINT. FILL MUPLET FOR
*  PARAMETERS SET WITHOUT INTERPOLATION.
*----
      IBM=0
      MAXNIS=0
      NISO(:NMIX)=0
      LISO(:NMIX)=.TRUE.
      LDELTA(:NMIX)=.FALSE.
      ITODO(:NMIX,:MD2)=0
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('SCRDRV: CHARACTER DATA EXPECTED.')
   20 IF(TEXT12.EQ.'MIX') THEN
         MUPLET(:NPAR)=0
         MUTYPE(:NPAR)=0
         VALI(:NPAR)=0
         VALR(:NPAR,1)=0.0
         VALR(:NPAR,2)=0.0
         DO 30 I=1,NPAR
         VALH(I)=' '
   30    CONTINUE
         LCUB2(:NPAR)=LCUBIC
         CALL REDGET(INDIC,IBM,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('SCRDRV: INTEGER DATA EXPECTED.')
         IF(IBM.GT.NMIX) THEN
            WRITE(HSMG,'(27HSCRDRV: NMIX OVERFLOW (IBM=,I8,6H NMIX=,I8,
     1      2H).)') IBM,NMIX
            CALL XABORT(HSMG)
         ENDIF
         IBMOLD=1
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('SCRDRV: CHARACTER DATA EXPECTED.')
         IF(TEXT12.EQ.'FROM') THEN
            CALL REDGET(INDIC,IBMOLD,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.NE.1) CALL XABORT('SCRDRV: INTEGER DATA EXPECTED.')
            MIXC(IBM)=IBMOLD
            GO TO 10
         ELSE IF(TEXT12.EQ.'USE') THEN
            MIXC(IBM)=IBM
            GO TO 10
         ENDIF
         MIXC(IBM)=IBMOLD
         GO TO 20
      ELSE IF(TEXT12.EQ.'MICRO') THEN
         IF(IBM.EQ.0) CALL XABORT('SCRDRV: MIX NOT SET (1).')
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('SCRDRV: CHARACTER DATA EXPECTED.')
         IF(TEXT12.EQ.'ALL') THEN
            LISO(IBM)=.TRUE.
         ELSE IF(TEXT12.EQ.'ONLY') THEN
            LISO(IBM)=.FALSE.
         ENDIF
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('SCRDRV: CHARACTER DATA EXPECTED.')
   40    IF(TEXT12.EQ.'ENDMIX') THEN
            GO TO 20
         ELSE IF(TEXT12.EQ.'NOEV') THEN
            IF(NISO(IBM).EQ.0) CALL XABORT('SCRDRV: MISPLACED NOEV.')
            ITODO(IBM,NISO(IBM))=1
         ELSE
            NISO(IBM)=NISO(IBM)+1
            IF(NISO(IBM).GT.MD2) CALL XABORT('SCRDRV: MD2 OVERFLOW.')
            MAXNIS=MAX(MAXNIS,NISO(IBM))
            READ(TEXT12,'(2A4)') (HISO(I0,IBM,NISO(IBM)),I0=1,2)
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.EQ.2) THEN
               CONC(IBM,NISO(IBM))=FLOTT
            ELSE IF((INDIC.EQ.3).AND.(TEXT12.EQ.'*')) THEN
               CONC(IBM,NISO(IBM))=-99.99
            ELSE
               CALL XABORT('SCRDRV: INVALID HISO DATA.')
            ENDIF
         ENDIF
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('SCRDRV: CHARACTER DATA EXPECTED.')
         GO TO 40
      ELSE IF((TEXT12.EQ.'SET').OR.(TEXT12.EQ.'DELTA')) THEN
         IF(IBM.EQ.0) CALL XABORT('SCRDRV: MIX NOT SET (2).')
         ITYPE=0
         IF(TEXT12.EQ.'SET') THEN
            ITYPE=1
         ELSE IF(TEXT12.EQ.'DELTA') THEN
            ITYPE=2
            LDELTA(IBM)=.TRUE.
         ENDIF
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('SCRDRV: CHARACTER DATA EXPECTED.')
         IF((TEXT12.EQ.'LINEAR').OR.(TEXT12.EQ.'CUBIC')) THEN
            HCUBIC=TEXT12
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         ELSE
            HCUBIC=' '
         ENDIF
         IF(INDIC.NE.3) CALL XABORT('SCRDRV: CHARACTER DATA EXPECTED.')
         DO 50 I=1,NPAR
         IF(TEXT12.EQ.PARKEY(I)) THEN
            IPAR=I
            GO TO 60
         ENDIF
   50    CONTINUE
         CALL XABORT('SCRDRV: PARAMETER '//TEXT12//' NOT FOUND.')
   60    IF(HCUBIC.EQ.'LINEAR') THEN
            LCUB2(IPAR)=.FALSE.
         ELSE IF(HCUBIC.EQ.'CUBIC') THEN
            LCUB2(IPAR)=.TRUE.
         ENDIF
         LPSAP=LCMGID(IPSAP,'paramdescrip')
         CALL LCMGET(LPSAP,'NVALUE',NVALUE)
         IF(NVALUE(IPAR).GT.MAXVAL) CALL XABORT('SCRDRV: MAXVAL OVERFL'
     1   //'OW.')
         WRITE(RECNAM,'(''pval'',I8)') IPAR
         LPSAP=LCMGID(IPSAP,'paramvaleurs')
         CALL LCMLEN(LPSAP,RECNAM,LENGTH,ITYLCM)
         IF(LENGTH.EQ.0) THEN
            WRITE(HSMG,'(25HSCRDRV: GLOBAL PARAMETER ,A,9H NOT SET.)')
     1      PARKEY(IPAR)
            CALL XABORT(HSMG)
         ENDIF
         IF(PARFMT(IPAR).EQ.'ENTIER') THEN
            IF(ITYPE.NE.1) CALL XABORT('SCRDRV: SET MANDATORY WITH INT'
     1      //'EGER PARAMETERS.')
            CALL REDGET(INDIC,VALI(IPAR),FLOTT,TEXT12,DFLOTT)
            IF(INDIC.NE.1) CALL XABORT('SCRDRV: INTEGER DATA EXPECTED.')
            CALL LCMGET(LPSAP,RECNAM,VINTE)
            DO 70 J=1,NVALUE(IPAR)
            IF(VALI(IPAR).EQ.VINTE(J)) THEN
               MUPLET(IPAR)=J
               MUTYPE(IPAR)=ITYPE
               GO TO 10
            ENDIF
   70       CONTINUE
            WRITE(HSMG,'(26HSCRDRV: INTEGER PARAMETER ,A,9H WITH VAL,
     1      2HUE,I5,30H NOT FOUND IN SAPHYB DATABASE.)') PARKEY(IPAR),
     2      VALI(IPAR)
            CALL XABORT(HSMG)
         ELSE IF(PARFMT(IPAR).EQ.'FLOTTANT') THEN
            CALL REDGET(INDIC,NITMA,VALR(IPAR,1),TEXT12,DFLOTT)
            IF(INDIC.NE.2) CALL XABORT('SCRDRV: REAL DATA EXPECTED.')
            VALR(IPAR,2)=VALR(IPAR,1)
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.EQ.2) THEN
               VALR(IPAR,2)=FLOTT
               CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            ENDIF
            CALL LCMGET(LPSAP,RECNAM,VREAL)
            IF(VALR(IPAR,1).EQ.VALR(IPAR,2)) THEN
               DO 80 J=1,NVALUE(IPAR)
               IF(ABS(VALR(IPAR,1)-VREAL(J)).LE.REPS*ABS(VREAL(J))) THEN
                  MUPLET(IPAR)=J
                  IF(ITYPE.NE.1) MUPLET(IPAR)=-1
                  MUTYPE(IPAR)=ITYPE
                  GO TO 20
               ENDIF
   80          CONTINUE
            ENDIF
            IF(VALR(IPAR,1).LT.VREAL(1)) THEN
               WRITE(HSMG,'(23HSCRDRV: REAL PARAMETER ,A,10H WITH VALU,
     1         1HE,1P,E12.4,25H IS OUTSIDE THE DOMAIN (<,E12.4,1H))')
     2         PARKEY(IPAR),VALR(IPAR,1),VREAL(1)
               CALL XABORT(HSMG)
            ELSE IF(VALR(IPAR,2).GT.VREAL(NVALUE(IPAR))) THEN
               WRITE(HSMG,'(23HSCRDRV: REAL PARAMETER ,A,10H WITH VALU,
     1         1HE,1P,E12.4,25H IS OUTSIDE THE DOMAIN (>,E12.4,1H))')
     2         PARKEY(IPAR),VALR(IPAR,1),VREAL(NVALUE(IPAR))
               CALL XABORT(HSMG)
            ELSE IF(VALR(IPAR,1).GT.VALR(IPAR,2)) THEN
               WRITE(HSMG,'(23HSCRDRV: REAL PARAMETER ,A,9H IS DEFIN,
     1         7HED WITH,1P,E12.4,2H >,E12.4,1H.)') PARKEY(IPAR),
     2         VALR(IPAR,1),VALR(IPAR,2)
               CALL XABORT(HSMG)
            ENDIF
            MUPLET(IPAR)=-1
            MUTYPE(IPAR)=ITYPE
            GO TO 20
         ELSE IF(PARFMT(IPAR).EQ.'CHAINE') THEN
            IF(ITYPE.NE.1) CALL XABORT('SCRDRV: SET MANDATORY WITH STR'
     1      //'ING PARAMETERS.')
            CALL REDGET(INDIC,NITMA,FLOTT,VALH(IPAR),DFLOTT)
            IF(INDIC.NE.3) CALL XABORT('SCRDRV: STRING DATA EXPECTED.')
            CALL LCMGTC(LPSAP,RECNAM,12,NVALUE(IPAR),VCHAR)
            DO 90 J=1,NVALUE(IPAR)
            IF(VALH(IPAR).EQ.VCHAR(J)) THEN
               MUPLET(IPAR)=J
               MUTYPE(IPAR)=ITYPE
               GO TO 10
            ENDIF
   90       CONTINUE
            WRITE(HSMG,'(25HSCRDRV: STRING PARAMETER ,A,10H WITH VALU,
     1      2HE ,A12,30H NOT FOUND IN SAPHYB DATABASE.)') PARKEY(IPAR),
     2      VALH(IPAR)
            CALL XABORT(HSMG)
         ELSE
            CALL XABORT('SCRDRV: INVALID FORMAT='//PARFMT(IPAR))
         ENDIF
      ELSE IF(TEXT12.EQ.'ENDMIX') THEN
*----
*  COMPUTE THE TERP FACTORS USING TABLE-OF-CONTENT INFORMATION.
*----
         IF(IMPX.GT.0) THEN
           DO IPAR=1,NPAR
             IF(PARFMT(IPAR).EQ.'FLOTTANT')THEN
               IF(LCUB2(IPAR)) THEN
                 WRITE(IOUT,'(26H SCRDRV: GLOBAL PARAMETER:,A12,5H ->CU,
     1           18HBIC INTERPOLATION.)') PARKEY(IPAR)
               ELSE
                 WRITE(IOUT,'(26H SCRDRV: GLOBAL PARAMETER:,A12,5H ->LI,
     1           19HNEAR INTERPOLATION.)') PARKEY(IPAR)
               ENDIF
             ENDIF
           ENDDO
         ENDIF
         IF(IBMOLD.GT.NMIL)CALL XABORT('SCRDRV: MIX OVERFLOW (SAPHYB).')
         IF(IBM.GT.NMIX)CALL XABORT('SCRDRV: MIX OVERFLOW (MICROLIB).')
         IF(NCAL.EQ.1) THEN
           TERP(1,IBM)=1.0
         ELSE
           CALL SCRTRP(IPSAP,LCUB2,IMPX,NVP,NPAR,NCAL,MUPLET,MUTYPE,
     1     VALR,0.0,TERP(1,IBM))
         ENDIF
         IBM=0
      ELSE IF((TEXT12.EQ.'SAPHYB').OR.(TEXT12.EQ.'TABLE').OR.
     1   (TEXT12.EQ.'CHAIN').OR.(TEXT12.EQ.';')) THEN
*----
*  CHECK TERP FACTORS AND RETURN
*----
         IF(TEXT12.EQ.';') ITER=0
         IF(TEXT12.EQ.'SAPHYB') ITER=1
         IF(TEXT12.EQ.'TABLE') ITER=2
         IF(TEXT12.EQ.'CHAIN') ITER=3
         DO 150 IBM=1,NMIX
         IBMOLD=MIXC(IBM)
         IF(IBMOLD.EQ.0) GO TO 150
         IF(NISO(IBM).GT.MAXNIS) CALL XABORT('SCRDRV: MAXNIS OVERFLOW.')
         IF(LDELTA(IBM)) THEN
            SUM=0.0
         ELSE
            SUM=1.0
         ENDIF
         DO 140 ICAL=1,NCAL
         SUM=SUM-TERP(ICAL,IBM)
  140    CONTINUE
         IF(ABS(SUM).GT.1.0E-4) THEN
            WRITE(HSMG,'(43HSCRDRV: INVALID INTERPOLATION FACTORS IN MI,
     1      5HXTURE,I4,1H.)') IBM
            CALL XABORT(HSMG)
         ENDIF
  150    CONTINUE
         GO TO 160
      ELSE
         CALL XABORT('SCRDRV: '//TEXT12//' IS AN INVALID KEYWORD.')
      ENDIF
      GO TO 10
*----
*  PRINT INTERPOLATION (TERP) FACTORS
*----
  160 IF(IMPX.GT.2) THEN
        WRITE(IOUT,'(/30H SCRDRV: INTERPOLATION FACTORS)')
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
