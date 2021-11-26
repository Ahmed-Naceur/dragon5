*DECK NCRRGR
      SUBROUTINE NCRRGR(IPCPO,IPMAP,LCUBIC,NMIX,IMPX,NMIL,NCAL,NCH,NB,
     1 NFUEL,NPARM,ITER,MAXNIS,MIXC,TERP,NISO,LISO,HISO,CONC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute TERP factors for multicompo interpolation. Use global and
* local parameters from a fuel-map object and optional user-defined
* values.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert, D. Sekki, R. Chambon
*
*Parameters: input
* IPCPO   address of the multicompo object.
* IPMAP   address of the fuel-map object.
* LCUBIC  =.TRUE.: cubic Ceschino interpolation; =.FALSE: linear
*         Lagrange interpolation.
* NMIX    number of material mixtures in the fuel-map macrolib.
* IMPX    printing index (=0 for no print).
* NMIL    number of material mixtures in the multicompo.
* NCAL    number of elementary calculations in the multicompo.
* NCH     number of reactor channels.
* NB      number of fuel bundles per channel.
* NFUEL   number of fuel types.
* NPARM   number of additional parameters (other than burnup) defined
*         in FMAP object
*
*Parameters: output
* ITER    completion flag (=0: all over; =1: use another multicompo;
*         =2 use another L_MAP + multicompo).
* MAXNIS  maximum value of NISO(I) in user data.
* MIXC    mixture index in the multicompo corresponding to each microlib
*         mixture.
* TERP    interpolation factors.
* NISO    number of user-selected isotopes.
* LISO    type of treatment (=.true.: ALL; =.false.: ONLY).
* HISO    name of the user-selected isotopes.
* CONC    user-defined number density of the user-selected isotopes. A
*         value of -99.99 is set to indicate that the multicompo value
*         is used.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER, PARAMETER::MAXISD=200
      TYPE(C_PTR) IPCPO,IPMAP
      INTEGER NMIX,IMPX,NMIL,NCAL,NFUEL,NCH,NB,ITER,MAXNIS,
     1        MIXC(NMIX),NPARM,HISO(2,NMIX,MAXISD),NISO(NMIX)
      REAL TERP(NCAL,NMIX),CONC(NMIX,MAXISD)
      LOGICAL LCUBIC,LISO(NMIX)
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::MAXADD=10
      INTEGER, PARAMETER::MAXLIN=50
      INTEGER, PARAMETER::MAXPAR=50
      INTEGER, PARAMETER::MAXVAL=200
      INTEGER, PARAMETER::NSTATE=40
      REAL, PARAMETER::REPS=1.0E-4
      REAL BURN0, BURN1, FLOTT, SUM, VALR1, VALR2, VARVAL
      INTEGER I0, IBMB, IBME, IBMOLD, IBM, IBTYP, IB, ICAL, ICH, IFUEL,
     & ILONG, IMIX, IMPY, INDIC, IPARTM, IPAR, ISO, ITYLCM, ITYPE, ITYP,
     & IVARTY, I, JBM, JB, JCAL, JPARM, JPAR, J, LENGTH, NCOMLI, NISOMI,
     & NITMA, NLOC, NMIXA, NPARMP, NPAR, NTOT, N
      CHARACTER TEXT12*12,PARKEY(MAXPAR)*12,PARFMT(MAXPAR)*8,
     1 PARKEL(MAXPAR)*12,HSMG*131,COMMEN(MAXLIN)*80,VALH(MAXPAR)*12,
     2 RECNAM*12,VCHAR(MAXVAL)*12,PARNAM*12,HCUBIC*12
      INTEGER ISTATE(NSTATE),VALI(MAXPAR),NVALUE(MAXPAR),VINTE(MAXVAL),
     1 MUPLET(2*MAXPAR),MUTYPE(2*MAXPAR),MAPLET(2*MAXPAR,MAXADD),
     2 MATYPE(2*MAXPAR,MAXADD),IDLTA(2*MAXPAR,MAXADD),NDLTA(2*MAXPAR),
     3 IDLTA1,MUPLT2(2*MAXPAR),MUTYP2(2*MAXPAR),HISOMI(2,MAXISD)
      DOUBLE PRECISION DFLOTT
      REAL VALR(2*MAXPAR,2),VREAL(MAXVAL),VALRA(2*MAXPAR,2,MAXADD),
     1 CONCMI(MAXISD)
      LOGICAL LDELT(2*MAXPAR),LDELT1,LSET(2*MAXPAR),LADD(2*MAXPAR),
     1 LSET1,LADD1,LDMAP(2*MAXPAR,2),LAMAP(2*MAXPAR,2,MAXADD),
     2 LCUB2(2*MAXPAR),LTST,LISOMI,LASBLY
      TYPE(C_PTR) JPMAP,KPMAP,JPCPO,KPCPO,LPCPO
      INTEGER, ALLOCATABLE, DIMENSION(:) :: FMIX,ZONEC,MIXA
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ZONEDP
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
      ALLOCATE(LPARM(NPARM+1),FMIX(NCH*NB),ZONEDP(NCH,NB),ZONEC(NCH),
     1 BRN0(NCH*NB),BRN1(NCH*NB),WPAR(NCH*NB,NPARM),LDELTA(NMIX),
     2 HPAR(NPARM+1))
*----
*  RECOVER TABLE-OF-CONTENT INFORMATION FOR THE MULTICOMPO.
*----
      CALL LCMGET(IPCPO,'STATE-VECTOR',ISTATE)
      NPAR=ISTATE(5)
      NLOC=ISTATE(6)
      NCOMLI=ISTATE(10)
* ASBLY : 
      LASBLY=.FALSE.
      CALL LCMGTC(IPCPO,'COMMENT',80,NCOMLI,COMMEN)
      IF(NPAR.GT.0)THEN
        CALL LCMSIX(IPCPO,'GLOBAL',1)
        CALL LCMGTC(IPCPO,'PARKEY',12,NPAR,PARKEY)
        CALL LCMGTC(IPCPO,'PARFMT',8,NPAR,PARFMT)
        CALL LCMGET(IPCPO,'NVALUE',NVALUE)
        IF(IMPX.GT.0)THEN
          DO IPAR=1,NPAR
            WRITE(RECNAM,'(''pval'',I8.8)') IPAR
            IF(PARFMT(IPAR).EQ.'INTEGER') THEN
              CALL LCMGET(IPCPO,RECNAM,VINTE)
              WRITE(IOUT,'(13H NCRRGR: KEY=,A,18H TABULATED POINTS=,
     1        1P,6I12/(43X,6I12))') PARKEY(IPAR),(VINTE(I),I=1,
     2        NVALUE(IPAR))
            ELSE IF(PARFMT(IPAR).EQ.'REAL') THEN
              CALL LCMGET(IPCPO,RECNAM,VREAL)
              WRITE(IOUT,'(13H NCRRGR: KEY=,A,18H TABULATED POINTS=,
     1        1P,6E12.4/(43X,6E12.4))') PARKEY(IPAR),(VREAL(I),I=1,
     2        NVALUE(IPAR))
            ELSE IF(PARFMT(IPAR).EQ.'STRING') THEN
              CALL LCMGTC(IPCPO,RECNAM,12,NVALUE(IPAR),VCHAR)
              WRITE(IOUT,'(13H NCRRGR: KEY=,A,18H TABULATED POINTS=,
     1        1P,6A12/(43X,6A12))') PARKEY(IPAR),(VCHAR(I),I=1,
     2        NVALUE(IPAR))
            ENDIF
          ENDDO
        ENDIF
        CALL LCMSIX(IPCPO,' ',2)
      ENDIF
      IF(NLOC.GT.0)THEN
        CALL LCMSIX(IPCPO,'LOCAL',1)
        CALL LCMGTC(IPCPO,'PARKEY',12,NLOC,PARKEL)
        CALL LCMSIX(IPCPO,' ',2)
        JPCPO=LCMGID(IPCPO,'MIXTURES')
        DO IBMOLD=1,NMIL
          KPCPO=LCMGIL(JPCPO,IBMOLD)
          LPCPO=LCMGID(KPCPO,'TREE')
          CALL LCMGET(LPCPO,'NVALUE',NVALUE)
          IF(IMPX.GT.0)THEN
            WRITE(IOUT,'(17H NCRRGR: MIXTURE=,I6)') IBMOLD
            DO IPAR=1,NLOC
              WRITE(RECNAM,'(''pval'',I8.8)') IPAR
              CALL LCMGET(LPCPO,RECNAM,VREAL)
              WRITE(IOUT,'(13H NCRRGR: KEY=,A,18H TABULATED POINTS=,
     1        1P,6E12.4/(43X,6E12.4))') PARKEL(IPAR),(VREAL(I),I=1,
     2        NVALUE(IPAR))
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      IF(IMPX.GT.0) THEN
        WRITE(IOUT,'(43H NCRRGR: NUMBER OF CALCULATIONS IN MULTICOM,
     1  3HPO=,I5)') NCAL
        WRITE(IOUT,'(43H NCRRGR: NUMBER OF MATERIAL MIXTURES IN MUL,
     1  8HTICOMPO=,I5)') NMIL
        WRITE(IOUT,'(43H NCRRGR: NUMBER OF MATERIAL MIXTURES IN FUE,
     1  6HL MAP=,I6)') NMIX
        WRITE(IOUT,'(1X,A)') (COMMEN(I),I=1,NCOMLI)
      ENDIF
      TERP(:NCAL,:NMIX)=0.0
      MIXC(:NMIX)=0
*----
*  READ (INTERP_DATA) AND SET VALI, VALR AND VALH PARAMETERS
*  CORRESPONDING TO THE INTERPOLATION POINT. FILL MUPLET FOR
*  PARAMETERS SET WITHOUT INTERPOLATION.
*----
      IBM=0
      IBMB=0
      IBME=0
      MAXNIS=0
      NISOMI=0
      LISOMI=.TRUE.
      LDELT1=.FALSE.
      LADD1=.FALSE.
      NISO(:NMIX)=0
      LISO(:NMIX)=.TRUE.
      LDELTA(:NMIX)=.FALSE.
      IDLTA1=0
      DO I=1,2*MAXPAR
         LSET(I)=.FALSE.
         LDELT(I)=.FALSE.
         LADD(I)=.FALSE.
         LDMAP(I,:2)=.FALSE.
         LAMAP(I,:2,:MAXADD)=.FALSE.
         NDLTA(I)=0
      ENDDO
*----
*  READ THE PARKEY NAME OF THE BURNUP FOR THIS MULTICOMPO.
*----
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF(INDIC.NE.3)CALL XABORT('NCRRGR: CHARACTER DATA EXPECTED(1).')
      IF((TEXT12.EQ.'MIX').OR.(TEXT12.EQ.';')) THEN
        NPARMP=NPARM
        GO TO 30
      ELSE
*       add burnup to parameters
        NPARMP=NPARM+1
        HPAR(NPARMP)=TEXT12
      ENDIF
*----
*  MAIN LOOP OF THE SUBROUTINE (UNTIL THE END)
*----
   20 CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF(INDIC.NE.3)CALL XABORT('NCRRGR: CHARACTER DATA EXPECTED(2).')
   30 IF(TEXT12.EQ.'MIX')THEN
        NISOMI=0
        LISOMI=.TRUE.
        IVARTY=0
        IBTYP=0
        MUPLET(:NPAR+NLOC)=0
        MUTYPE(:NPAR+NLOC)=0
        VALI(:NPAR)=0
        VALR(:NPAR+NLOC,1)=0.0
        VALR(:NPAR+NLOC,2)=0.0
        DO 35 I=1,MAXADD
          MAPLET(:NPAR+NLOC,I)=0
          MATYPE(:NPAR+NLOC,I)=0
          VALRA(:NPAR+NLOC,1,I)=0.0
          VALRA(:NPAR+NLOC,2,I)=0.0
   35   CONTINUE
        DO I=1,2*MAXPAR
          LSET(I)=.FALSE.
          LDELT(I)=.FALSE.
          LADD(I)=.FALSE.
          LDMAP(I,:2)=.FALSE.
          LAMAP(I,:2,:MAXADD)=.FALSE.
        ENDDO
        DO 40 I=1,NPAR
        VALH(I)=' '
   40   CONTINUE
        LCUB2(:NPAR+NLOC)=LCUBIC
        CALL REDGET(INDIC,IBM,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.1)CALL XABORT('NCRRGR: INTEGER DATA EXPECTED.')
*       CHECK FUEL MIXTURE
        JPMAP=LCMGID(IPMAP,'FUEL')
        DO IFUEL=1,NFUEL
          KPMAP=LCMGIL(JPMAP,IFUEL)
          CALL LCMGET(KPMAP,'MIX',IMIX)
          IF(IMIX.EQ.IBM)GOTO 50
        ENDDO
        WRITE(IOUT,*)'NCRRGR: UNABLE TO FIND FUEL MIXTURE ',IBM
        CALL XABORT('NCRRGR: WRONG MIXTURE NUMBER.')
   50   IBMOLD=1
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.3)CALL XABORT('NCRRGR: CHARACTER DATA EXPECTED(3).')
        IF(TEXT12.EQ.'FROM')THEN
          CALL REDGET(INDIC,IBMOLD,FLOTT,TEXT12,DFLOTT)
          IF(INDIC.NE.1)CALL XABORT('NCRRGR: INTEGER DATA EXPECTED.')
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          IF(INDIC.NE.3)CALL XABORT('NCRRGR: CHARACTER DATA EXPECTE'
     1     //'D.')
        ELSE IF(TEXT12.EQ.'USE') THEN
          IBMOLD=IBM
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          IF(INDIC.NE.3)CALL XABORT('NCRRGR: CHARACTER DATA EXPECTE'
     1     //'D.')
* ASBLY: automatically assembly-wise unfolded geometry
        ELSE IF(TEXT12.EQ.'ASBLY') THEN
          IF(LASBLY) DEALLOCATE(MIXA)
          IBMOLD=1
          LASBLY=.TRUE.
          JPMAP=LCMGID(IPMAP,'GEOMAP')
          CALL LCMGET(JPMAP,'STATE-VECTOR',ISTATE)
          CALL LCMLEN(JPMAP,'MIX-ASBLY',NITMA,INDIC)
          IF(NITMA.EQ.0)CALL XABORT('NCRRGR: No assembly defined')
          NMIXA=NITMA/2
*          NMIXA=ISTATE(39)
          ALLOCATE(MIXA(2*NMIXA))
          CALL LCMGET(JPMAP,'MIX-ASBLY',MIXA)
          DO I=1,NMIXA
            IF(IBM.EQ.MIXA(I)) THEN
              IBMB=MIXA(I+NMIXA)
              IBME=IBMB+NMIL-1
              CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
              IF(INDIC.NE.3)CALL XABORT('NCRRGR: CHARACTER DATA EXPECTE'
     1     //'D.')
              GOTO 30
            ENDIF
          ENDDO
          CALL XABORT('NCRRGR: WRONG ASSEMBLY MIXTURE.')
        ENDIF
* ASBLY: automatically assembly-wise unfolded geometry
        IBMB=IBM
        IBME=IBM
        GOTO 30
      ELSEIF(TEXT12.EQ.'MICRO')THEN
        IF(IBM.EQ.0) CALL XABORT('NCRRGR: MIX NOT SET (1).')
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.3)CALL XABORT('NCRRGR: CHARACTER DATA EXPECTED(4).')
        IF(TEXT12.EQ.'ALL')THEN
          LISOMI=.TRUE.
        ELSEIF(TEXT12.EQ.'ONLY')THEN
          LISOMI=.FALSE.
        ENDIF
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.3)CALL XABORT('NCRRGR: CHARACTER DATA EXPECTED(5).')
   60   IF(TEXT12.EQ.'ENDMIX')THEN
          GOTO 30
        ELSE
          NISOMI=NISOMI+1
          IF(NISOMI.GT.MAXISD) CALL XABORT('NCRRGR: MAXISD OVERFLOW.')
          MAXNIS=MAX(MAXNIS,NISOMI)
          READ(TEXT12,'(2A4)') (HISOMI(I0,NISOMI),I0=1,2)
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          IF(INDIC.EQ.2)THEN
            CONCMI(NISOMI)=FLOTT
          ELSEIF((INDIC.EQ.3).AND.(TEXT12.EQ.'*'))THEN
            CONCMI(NISOMI)=-99.99
          ELSE
            CALL XABORT('NCRRGR: INVALID HISO DATA.')
          ENDIF
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          IF(INDIC.NE.3)CALL XABORT('NCRRGR: CHARACTER DATA EXPECTED.')
          GOTO 60
        ENDIF
      ELSEIF((TEXT12.EQ.'SET').OR.(TEXT12.EQ.'DELTA').OR.
     1       (TEXT12.EQ.'ADD'))THEN
        IF(IBM.EQ.0) CALL XABORT('NCRRGR: MIX NOT SET (2).')
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
          DO 65 JPAR=1,NPAR+NLOC
            MAPLET(JPAR,IDLTA1)=MUPLET(JPAR)
            MATYPE(JPAR,IDLTA1)=MUTYPE(JPAR)
   65     CONTINUE
        ENDIF
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.3)CALL XABORT('NCRRGR: CHARACTER DATA EXPECTED(7).')
        IF((TEXT12.EQ.'LINEAR').OR.(TEXT12.EQ.'CUBIC')) THEN
          HCUBIC=TEXT12
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        ELSE
          HCUBIC=' '
        ENDIF
        IF(INDIC.NE.3)CALL XABORT('NCRRGR: CHARACTER DATA EXPECTED(8).')
* check if parameter is global
        IPAR=-99
        DO I=1,NPAR
          IF(TEXT12.EQ.PARKEY(I))THEN
            IPAR=I
            LPCPO=LCMGID(IPCPO,'GLOBAL')
            IPARTM=IPAR
            PARNAM=TEXT12
            GOTO 70
          ENDIF
        ENDDO
* check if parameter is local
        DO I=1,NLOC
          IF(TEXT12.EQ.PARKEL(I))THEN
            IPAR=NPAR+I
            JPCPO=LCMGID(IPCPO,'MIXTURES')
            KPCPO=LCMGIL(JPCPO,IBMOLD)
            LPCPO=LCMGID(KPCPO,'TREE')
            IPARTM=IPAR-NPAR
            PARNAM=TEXT12
            GOTO 70
          ENDIF
        ENDDO
        WRITE(HSMG,'(18HNCRRGR: PARAMETER ,A,14H NOT FOUND(1).)') TEXT12
        CALL XABORT(HSMG)
*
   70   IF(HCUBIC.EQ.'LINEAR') THEN
          LCUB2(IPAR)=.FALSE.
        ELSE IF(HCUBIC.EQ.'CUBIC') THEN
          LCUB2(IPAR)=.TRUE.
        ENDIF
        WRITE(RECNAM,'(''pval'',I8.8)') IPARTM
        CALL LCMGET(LPCPO,'NVALUE',NVALUE)
        IF(NVALUE(IPARTM).GT.MAXVAL)CALL XABORT('NCRRGR: MAXVAL OVERFL'
     1   //'OW.')
        CALL LCMLEN(LPCPO,RECNAM,LENGTH,ITYLCM)
        IF(LENGTH.EQ.0)THEN
          WRITE(HSMG,'(25HNCRRGR: GLOBAL PARAMETER ,A,9H NOT SET.)')
     1     PARNAM
          CALL XABORT(HSMG)
        ENDIF
        IF((IPAR.GT.NPAR).OR.
     1   ((IPAR.LE.NPAR).AND.(PARFMT(IPAR).EQ.'REAL')))THEN
          VALR1=VREAL(1)
          VALR2=VREAL(NVALUE(IPAR))
          CALL LCMGET(LPCPO,RECNAM,VREAL)
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
              IF(NDLTA(IPAR).GT.MAXADD) CALL XABORT('NCRRGR: MAXADD OV'
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
              IF(NDLTA(IPAR).GT.MAXADD) CALL XABORT('NCRRGR: MAXADD OV'
     1        //'ERFLOW.')
              LAMAP(IPAR,1,NDLTA(IPAR))=.TRUE.
              IDLTA(IPAR,NDLTA(IPAR))=IDLTA1
            ENDIF
            IF(LSET1.AND.(.NOT.LSET(IPAR))) GO TO 20
          ELSE
            CALL XABORT('NCRRGR: real value or "MAP" expected(1).')
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
              CALL XABORT('NCRRGR: real value or "MAP" expected(2).')
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
            DO J=1,NVALUE(IPARTM)
              IF(ABS(VALR(IPAR,1)-VREAL(J)).LE.REPS*ABS(VREAL(J)))THEN
                MUPLET(IPAR)=J
                GOTO 30
              ENDIF
            ENDDO
          ENDIF
*----
*  ERRORS HANDLING
*----
          IF(VALR1.LT.VREAL(1))THEN
*           OUTSIDE OF THE DOMAIN (1)
            WRITE(HSMG,'(23HNCRRGR: REAL PARAMETER ,A,10H WITH VALU,
     1       1HE,1P,E12.4,26H IS OUTSIDE THE DOMAIN.(1))') PARNAM,VALR1
             WRITE(6,*)'Domain:',VREAL(1),' <-> ',VREAL(NVALUE(IPAR))
            CALL XABORT(HSMG)
          ELSEIF(VALR2.GT.VREAL(NVALUE(IPARTM)))THEN
*           OUTSIDE OF THE DOMAIN (2)
            WRITE(HSMG,'(23HNCRRGR: REAL PARAMETER ,A,10H WITH VALU,
     1       1HE,1P,E12.4,26H IS OUTSIDE THE DOMAIN.(2))') PARNAM,VALR2
             WRITE(6,*)'Domain:',VREAL(1),' <-> ',VREAL(NVALUE(IPAR))
            CALL XABORT(HSMG)
          ELSEIF((VALR1.GT.VALR2).AND.(ITYPE.EQ.1))THEN
*           ITYPE=1 correspond to an integral between VALR1 and VALR2
*           otherwise it is a simple difference
            WRITE(HSMG,'(23HNCRRGR: REAL PARAMETER ,A,9H IS DEFIN,
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
            DO I=1,NLOC
              IF(TEXT12.EQ.PARKEL(I))THEN
                IPAR=NPAR+I
                GOTO 130
              ENDIF
            ENDDO
            CALL XABORT('NCRRGR: PARAMETER '//TEXT12//' NOT FOUND(2).')
  130       CONTINUE
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.EQ.2)THEN
              VALRA(IPAR,1,IDLTA1)=FLOTT
              VALRA(IPAR,2,IDLTA1)=FLOTT
              IF(IPAR.LE.NPAR)THEN
                LPCPO=LCMGID(IPCPO,'GLOBAL')
                IPARTM=IPAR
              ELSE
                JPCPO=LCMGID(IPCPO,'MIXTURES')
                KPCPO=LCMGIL(JPCPO,IBMOLD)
                LPCPO=LCMGID(KPCPO,'TREE')
                IPARTM=IPAR-NPAR
              ENDIF
              WRITE(RECNAM,'(''pval'',I8.8)') IPARTM
              CALL LCMGET(LPCPO,RECNAM,VREAL)
              CALL LCMGET(LPCPO,'NVALUE',NVALUE)
              MAPLET(IPAR,IDLTA1)=-1
              MATYPE(IPAR,IDLTA1)=1
              DO J=1,NVALUE(IPARTM)
                IF(ABS(VALRA(IPAR,1,IDLTA1)-VREAL(J)).LE.
     1             REPS*ABS(VREAL(J)))THEN
                  MAPLET(IPAR,IDLTA1)=J
                  GOTO 120
                ENDIF
              ENDDO
            ELSEIF(TEXT12.EQ.'SAMEASREF')THEN
              MAPLET(IPAR,IDLTA1)=-1
              MATYPE(IPAR,IDLTA1)=-1
            ELSE
              CALL XABORT('NCRRGR: REAL or "SAMEASREF" expected')
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
            DO I=1,NLOC
              IF(TEXT12.EQ.PARKEL(I))THEN
                IPAR=NPAR+I
                GOTO 160
              ENDIF
            ENDDO
            CALL XABORT('NCRRGR: PARAMETER '//TEXT12//' NOT FOUND(3).')
  160       CONTINUE
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.EQ.2)THEN
              VALR(IPAR,1)=FLOTT
              VALR(IPAR,2)=FLOTT
              IF(IPAR.LE.NPAR)THEN
                LPCPO=LCMGID(IPCPO,'GLOBAL')
                IPARTM=IPAR
              ELSE
                JPCPO=LCMGID(IPCPO,'MIXTURES')
                KPCPO=LCMGIL(JPCPO,IBMOLD)
                LPCPO=LCMGID(KPCPO,'TREE')
                IPARTM=IPAR-NPAR
              ENDIF
              WRITE(RECNAM,'(''pval'',I8.8)') IPARTM
              CALL LCMGET(LPCPO,RECNAM,VREAL)
              CALL LCMGET(LPCPO,'NVALUE',NVALUE)
              MUPLET(IPAR)=-1
              MUTYPE(IPAR)=1
              DO J=1,NVALUE(IPARTM)
                IF(ABS(VALR(IPAR,1)-VREAL(J)).LE.REPS*ABS(VREAL(J)))THEN
                  MUPLET(IPAR)=J
                  GOTO 150
                ENDIF
              ENDDO
            ELSEIF(TEXT12.EQ.'SAMEASREF')THEN
              MUPLET(IPAR)=-1
              MUTYPE(IPAR)=-1
            ELSE
              CALL XABORT('NCRRGR: REAL or "SAMEASREF" expected')
            ENDIF
            GOTO 150
  170       CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          ENDIF
          GOTO 30
        ELSEIF(PARFMT(IPAR).EQ.'INTEGER')THEN
          IF(ITYPE.NE.1)CALL XABORT('NCRRGR: SET MANDATORY WITH INT'
     1     //'EGER PARAMETERS.')
          CALL REDGET(INDIC,VALI(IPAR),FLOTT,TEXT12,DFLOTT)
          IF(INDIC.NE.1)CALL XABORT('NCRRGR: INTEGER DATA EXPECTED.')
          CALL LCMGET(LPCPO,RECNAM,VINTE)
          DO 175 J=1,NVALUE(IPAR)
            IF(VALI(IPAR).EQ.VINTE(J))THEN
              MUPLET(IPAR)=J
              MUTYPE(IPAR)=ITYPE
              GOTO 20
            ENDIF
  175     CONTINUE
          WRITE(HSMG,'(26HNCRRGR: INTEGER PARAMETER ,A,9H WITH VAL,
     1     2HUE,I5,34H NOT FOUND IN MULTICOMPO DATABASE.)')
     2     PARKEY(IPAR), VALI(IPAR)
           CALL XABORT(HSMG)
        ELSEIF(PARFMT(IPAR).EQ.'STRING')THEN
          IF(ITYPE.NE.1)CALL XABORT('NCRRGR: SET MANDATORY WITH STR'
     1     //'ING PARAMETERS.')
          CALL REDGET(INDIC,NITMA,FLOTT,VALH(IPAR),DFLOTT)
          IF(INDIC.NE.3)CALL XABORT('NCRRGR: STRING DATA EXPECTED.')
          CALL LCMGTC(LPCPO,RECNAM,12,NVALUE(IPAR),VCHAR)
          DO 180 J=1,NVALUE(IPAR)
            IF(VALH(IPAR).EQ.VCHAR(J))THEN
              MUPLET(IPAR)=J
              MUTYPE(IPAR)=ITYPE
              GOTO 20
            ENDIF
  180     CONTINUE
          WRITE(HSMG,'(25HNCRRGR: STRING PARAMETER ,A,10H WITH VALU,
     1     1HE,A12,34H NOT FOUND IN MULTICOMPO DATABASE.)')
     2     PARKEY(IPAR), VALH(IPAR)
          CALL XABORT(HSMG)
        ENDIF
      ELSEIF(TEXT12.EQ.'TIMAV-BURN')THEN
        IF(IBM.EQ.0) CALL XABORT('NCRRGR: MIX NOT SET (3).')
        IBTYP=1
      ELSEIF(TEXT12.EQ.'INST-BURN')THEN
        IF(IBM.EQ.0) CALL XABORT('NCRRGR: MIX NOT SET (4).')
        IBTYP=2
      ELSEIF(TEXT12.EQ.'AVG-EX-BURN')THEN
        IF(IBM.EQ.0) CALL XABORT('NCRRGR: MIX NOT SET (5).')
        IBTYP=3
        CALL REDGET(INDIC,IVARTY,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.1)CALL XABORT('NCRRGR: INTEGER DATA EXPECTED.')
      ELSEIF(TEXT12.EQ.'ENDMIX')THEN
*----
*  RECOVER FUEL-MAP INFORMATION.
*----
        IF(IMPX.GT.0) THEN
          DO IPAR=1,NPAR
            IF(PARFMT(IPAR).EQ.'REAL')THEN
              IF(LCUB2(IPAR)) THEN
                WRITE(IOUT,'(26H NCRRGR: GLOBAL PARAMETER:,A12,5H ->CU,
     1          18HBIC INTERPOLATION.)') PARKEY(IPAR)
              ELSE
                WRITE(IOUT,'(26H NCRRGR: GLOBAL PARAMETER:,A12,5H ->LI,
     1          19HNEAR INTERPOLATION.)') PARKEY(IPAR)
              ENDIF
            ENDIF
          ENDDO
          DO IPAR=1,NLOC
            IF(LCUB2(NPAR+IPAR)) THEN
              WRITE(IOUT,'(25H NCRRGR: LOCAL PARAMETER:,A12,8H ->CUBIC,
     1        14HINTERPOLATION.)') PARKEL(IPAR)
            ELSE
              WRITE(IOUT,'(25H NCRRGR: LOCAL PARAMETER:,A12,8H ->LINEA,
     1        16HR INTERPOLATION.)') PARKEL(IPAR)
            ENDIF
          ENDDO
        ENDIF
        FMIX(:NCH*NB)=0
        CALL LCMGET(IPMAP,'FLMIX',FMIX)
        CALL NCRMAP(IPMAP,NPARM,HPAR,NCH,NB,IBTYP,IMPX,BRN0,BRN1,WPAR,
     1              LPARM)
        IF(IBTYP.EQ.3) THEN
          IF(IVARTY.EQ.0) CALL XABORT('NCRRGR: IVARTY NOT SET.')
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
          IF (ILONG.EQ.0) CALL XABORT('NCRRGR: NO SAVED VALUES FOR '
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
                IF(IMPX.GT.0) WRITE(6,*) 'L_MAP values overwritten by '
     1          // 'the SET option for parameter '//HPAR(JPARM)
                IF(.NOT.LADD(IPAR)) LPARM(JPARM)=.FALSE.
              ENDIF
              GOTO 185
            ENDIF
          ENDDO
          DO I=1,NLOC
            IF(HPAR(JPARM).EQ.PARKEL(I))THEN
              IPAR=NPAR+I
              IF(LSET(IPAR)) THEN
                IF(IMPX.GT.0) WRITE(6,*) 'L_MAP values overwritten by '
     1          // 'the SET option for parameter '//HPAR(JPARM)
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
        DO 281 JB=1,NB
        DO 280 ICH=1,NCH
        IB=(JB-1)*NCH+ICH
        IF(FMIX(IB).EQ.0) GO TO 280
        NTOT=NTOT+1
* ASBLY: loop on multicompo mixtures        
        DO 285 IBM=IBMB,IBME
        IF(LASBLY) IBMOLD=IBM-IBMB+1      
* ASBLY: end        
        IPAR=-99
        IF(FMIX(IB).EQ.IBM)THEN
          IF(NTOT.GT.NMIX) CALL XABORT('NCRRGR: NMIX OVERFLOW.')
          DO 260 JPARM=1,NPARMP
          IF(.NOT.LPARM(JPARM))GOTO 260
* check if parameter is global
        DO I=1,NPAR
          IF(HPAR(JPARM).EQ.PARKEY(I))THEN
            IPAR=I
            LPCPO=LCMGID(IPCPO,'GLOBAL')
            IPARTM=IPAR
            PARNAM=HPAR(JPARM)
            GOTO 190
          ENDIF
        ENDDO
* check if parameter is local
        DO I=1,NLOC
          IF(HPAR(JPARM).EQ.PARKEL(I))THEN
            IPAR=NPAR+I
            JPCPO=LCMGID(IPCPO,'MIXTURES')
            KPCPO=LCMGIL(JPCPO,IBMOLD)
            LPCPO=LCMGID(KPCPO,'TREE')
            IPARTM=IPAR-NPAR
            PARNAM=HPAR(JPARM)
            GOTO 190
          ENDIF
        ENDDO
        WRITE(HSMG,'(18HNCRRGR: PARAMETER ,A,14H NOT FOUND(4).)')
     1  HPAR(JPARM)
        CALL XABORT(HSMG)
  190   CONTINUE
        WRITE(RECNAM,'(''pval'',I8.8)') IPARTM
        CALL LCMGET(LPCPO,'NVALUE',NVALUE)
        IF(NVALUE(IPARTM).GT.MAXVAL)CALL XABORT('NCRRGR: MAXVAL OVERFL'
     1   //'OW.')
        CALL LCMLEN(LPCPO,RECNAM,LENGTH,ITYLCM)
        IF(LENGTH.EQ.0)THEN
          WRITE(HSMG,'(25HNCRRGR: GLOBAL PARAMETER ,A,9H NOT SET.)')
     1    PARNAM
          CALL XABORT(HSMG)
        ENDIF
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
        LPCPO=LCMGID(IPCPO,'GLOBAL')
        CALL LCMGET(LPCPO,'NVALUE',NVALUE)
        IF(NVALUE(IPARTM).GT.MAXVAL) CALL XABORT('NCRRGR: MAXVAL OVE'
     1  //'RFLOW.')
        WRITE(RECNAM,'(''pval'',I8.8)') IPARTM
        CALL LCMLEN(LPCPO,RECNAM,LENGTH,ITYLCM)
        IF(LENGTH.EQ.0)THEN
          WRITE(HSMG,'(25HNCRRGR: GLOBAL PARAMETER ,A,9H NOT SET.)')
     1    PARNAM
          CALL XABORT(HSMG)
        ENDIF
        IF(LENGTH.GT.MAXVAL) CALL XABORT('NCRRGR: MAXVAL OVERFLOW.')
        CALL LCMGET(LPCPO,RECNAM,VREAL)
        IF(ITYPE.EQ.1)THEN
          IF(VALR1.EQ.VALR2)THEN
            DO J=1,NVALUE(IPARTM)
              IF(ABS(VALR1-VREAL(J)).LE.REPS*ABS(VREAL(J)))THEN
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
        IF(VALR1.LT.VREAL(1))THEN
*         OUTSIDE OF THE DOMAIN (1)
          WRITE(HSMG,'(23HNCRRGR: REAL PARAMETER ,A,10H WITH VALU,
     1    1HE,1P,E12.4,26H IS OUTSIDE THE DOMAIN(3).)') PARNAM,VALR1
          WRITE(6,*)'Domain:',VREAL(1),' <-> ',VREAL(NVALUE(IPAR))
          CALL XABORT(HSMG)
        ELSEIF(VALR2.GT.VREAL(NVALUE(IPAR)))THEN
*         OUTSIDE OF THE DOMAIN (2)
          WRITE(HSMG,'(23HNCRRGR: REAL PARAMETER ,A,10H WITH VALU,
     1    1HE,1P,E12.4,26H IS OUTSIDE THE DOMAIN(4).)') PARNAM,VALR2
          WRITE(6,*)'Domain:',VREAL(1),' <-> ',VREAL(NVALUE(IPAR))
          CALL XABORT(HSMG)
        ELSEIF((ITYPE.EQ.1).AND.(VALR1.GT.VALR2))THEN
*         VALR1 > VALR2
          WRITE(HSMG,'(23HNCRRGR: REAL PARAMETER ,A,9H IS DEFIN,
     1    7HED WITH,1P,E12.4,2H >,E12.4,4H.(2))') PARNAM,
     2    VALR1,VALR2
          CALL XABORT(HSMG)
        ENDIF
*----
*  COMPUTE THE TERP FACTORS USING TABLE-OF-CONTENT INFORMATION.
*----
  260     CONTINUE
          MIXC(NTOT)=IBMOLD
          IF(IBMOLD.GT.NMIL)
     1       CALL XABORT('NCRRGR: MIX OVERFLOW (COMPO).')
          IF(IMPY.GT.2) WRITE(6,'(32H NCRRGR: COMPUTE TERP FACTORS IN,
     1    12H NEW MIXTURE,I5,1H.)') NTOT
          NISO(NTOT)=NISOMI
          LISO(NTOT)=LISOMI
          LDELTA(NTOT)=LDELT1
          DO ISO=1,NISOMI
            HISO(1,NTOT,ISO)=HISOMI(1,ISO)
            HISO(2,NTOT,ISO)=HISOMI(2,ISO)
            CONC(NTOT,ISO)=CONCMI(ISO)
          ENDDO
          DO JPAR=1,NPAR+NLOC
            MUPLT2(JPAR)=MUPLET(JPAR)
          ENDDO
          IF(IBTYP.EQ.3)THEN
             IF(ZONEDP(ICH,JB).NE.0) THEN
                CALL NCRTRP(IPCPO,LCUB2,IMPY,IBMOLD,NPAR,NLOC,NCAL,
     1                      MUPLT2,MUTYPE,VALR,VARVAL,TERP(1,NTOT))
             ELSE
                TERP(:NCAL,NTOT)=0.0
             ENDIF
          ELSE
             CALL NCRTRP(IPCPO,LCUB2,IMPY,IBMOLD,NPAR,NLOC,NCAL,
     1                   MUPLT2,MUTYPE,VALR,VARVAL,TERP(1,NTOT))
          ENDIF
*         DELTA-ADD
          DO 270 IPAR=1,NPAR+NLOC
            IF(LADD(IPAR))THEN
              DO N=1,NDLTA(IPAR)
                IDLTA1=IDLTA(IPAR,N)
                DO JPAR=1,NPAR+NLOC
                  MUPLT2(JPAR)=MAPLET(JPAR,IDLTA1)
                  MUTYP2(JPAR)=MATYPE(JPAR,IDLTA1)
                ENDDO
                DO JPAR=1,NPAR+NLOC
                  IF(MUTYP2(JPAR).LT.0)THEN
                    MUPLT2(JPAR)=MUPLET(JPAR)
                    MUTYP2(JPAR)=MUTYPE(JPAR)
                    VALRA(JPAR,1,IDLTA1)=VALR(JPAR,1)
                    VALRA(JPAR,2,IDLTA1)=VALR(JPAR,2)
                  ENDIF
                ENDDO
                ALLOCATE(TERPA(NCAL))
                CALL NCRTRP(IPCPO,LCUB2,IMPY,IBMOLD,NPAR,NLOC,NCAL,
     1          MUPLT2,MUTYP2,VALRA(1,1,IDLTA1),VARVAL,TERPA)
                DO 275 JCAL=1,NCAL
                TERP(JCAL,NTOT)=TERP(JCAL,NTOT)+TERPA(JCAL)
  275           CONTINUE
                DEALLOCATE(TERPA)
              ENDDO
            ENDIF
  270     CONTINUE
        ENDIF
* ASBLY: next mixture        
  285   CONTINUE      
* ASBLY: end        
  280   CONTINUE
  281   CONTINUE
        IF(NTOT.GT.NMIX) CALL XABORT('NCRRGR: ALGORITHM FAILURE.')
        IBM=0
        IBMB=0
        IBME=0
      ELSEIF((TEXT12.EQ.'COMPO').OR.(TEXT12.EQ.'TABLE').OR.
     1  (TEXT12.EQ.';')) THEN
*----
*  CHECK TERP FACTORS AND RETURN
*----
        IF(TEXT12.EQ.';') ITER=0
        IF(TEXT12.EQ.'COMPO') ITER=1
        IF(TEXT12.EQ.'TABLE') ITER=2
        DO 300 IBM=1,NMIX
        IBMOLD=MIXC(IBM)
        IF(IBMOLD.EQ.0) GO TO 300
        IF(NISO(IBM).GT.MAXNIS) CALL XABORT('NCRRGR: MAXNIS OVERFLOW.')
        IF(LDELTA(IBM)) THEN
          SUM=0.0
        ELSE
          SUM=1.0
        ENDIF
        DO 290 ICAL=1,NCAL
        SUM=SUM-TERP(ICAL,IBM)
  290   CONTINUE
        IF(ABS(SUM).GT.1.0E-4) THEN
           WRITE(HSMG,'(43HNCRRGR: INVALID INTERPOLATION FACTORS IN MI,
     1     5HXTURE,I4,1H.)') IBM
           CALL XABORT(HSMG)
        ENDIF
  300   CONTINUE
*----
*  EXIT MAIN LOOP OF THE SUBROUTINE
*----
        GO TO 310
      ELSE
        CALL XABORT('NCRRGR: '//TEXT12//' IS AN INVALID KEYWORD.')
      ENDIF
      GOTO 20
*----
*  PRINT INTERPOLATION (TERP) FACTORS
*----
  310 IF(IMPX.GT.2) THEN
        WRITE(IOUT,'(/30H NCRRGR: INTERPOLATION FACTORS)')
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
      DEALLOCATE(HPAR,LDELTA,WPAR,BRN1,BRN0,ZONEC,ZONEDP,FMIX,LPARM)
      IF(LASBLY) DEALLOCATE(MIXA)
      RETURN
  320 FORMAT(6H CALC=,I8,6H TERP=,1P,8E13.5/(20X,8E13.5))
      END
