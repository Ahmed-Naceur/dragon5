*DECK SHISN3
      SUBROUTINE SHISN3 (IPLIB,IPTRK,IFTRAK,LEVEL,NGRO,NBISO,NBM,NREG,
     1 NUN,CDOOR,INRS,NBNRS,IMPX,ISONAM,MIX,DEN,SN,SB,LSHI,IPHASE,MAT,
     2 VOL,KEYFLX,LEAKSW,TITR,START,SIGT,SIGT3,NOCONV,BIEFF,LGC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one multidimensional self-shielding iteration using the
* generalized Stamm'ler algorithm with Nordheim (PIC) approximation.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIB   pointer to the internal microscopic cross section library
*         (L_LIBRARY signature).
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  unit number of the sequential binary tracking file.
* LEVEL   type of self-shielding model (=1 original Stamm'ler model
*         with Nordheim approximation; =2 Stamm'ler model with Nordheim
*         approximation and Riemann integration method).
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NBM     number of mixtures in the macrolib.
* NREG    number of volumes.
* NUN     number of unknowns in the flux or source vector in one
*         energy group.
* CDOOR   name of the geometry/solution module.
* INRS    index of the resonant isotope under consideration.
* NBNRS   number of totaly correlated resonant regions.
* IMPX    print flag.
* ISONAM  alias name of isotopes.
* MIX     mix number of each isotope (can be zero).
* DEN     density of each isotope.
* LSHI    resonant region number associated with each isotope.
*         Infinite dilution will be assumed if LSHI(i)=0.
* IPHASE  type of flux solution (=1 use a native flux solution door;
*         =2 use collision probabilities).
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  pointers of fluxes in unknown vector.
* LEAKSW  leakage flag (.TRUE. only if leakage is present on the outer
*         surface).
* TITR    title.
* START   beginning-of-iteration flag (.TRUE. if SHISN3 is called
*         for the first time).
* SIGT3   transport correction.
* NOCONV  mixture convergence flag. (.TRUE. if mixture IBM
*         is not converged in group L).
* BIEFF   Livolant-Jeanpierre normalization flag (.TRUE. to
*         activate).
* LGC     Goldstein-Cohen approximation flag (.TRUE. to activate).
*
*Parameters: input/output
* SN      on input, estimate of the dilution cross section in each 
*         energy group of each isotope. A value of 1.0e10 is used 
*         for infinitedilution.
*         On output, computed dilution cross section in each energy 
*         group of each isotope.
* SIGT    total macroscopic cross sections as modified by Shiba.
*
*Parameters: output
* SB      dilution cross section as used in Livolant-Jeanpierre
*         normalization.
*
*Reference:
* A. Hebert, Revisiting the Stamm'ler Self-Shielding Method, Presented
* at the 25th CNS annnual conference, June 6-9, Toronto, 2004.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      PARAMETER (NALPHA=9,NRAT=(NALPHA+1)/2)
      TYPE(C_PTR) IPLIB,IPTRK
      INTEGER IFTRAK,LEVEL,NGRO,NBISO,NBM,NREG,NUN,INRS,NBNRS,IMPX,
     1 ISONAM(3,NBISO),MIX(NBISO),LSHI(NBISO),IPHASE,MAT(NREG),
     2 KEYFLX(NREG)
      REAL DEN(NBISO),SN(NGRO,NBISO),SB(NGRO,NBISO),VOL(NREG),
     1 SIGT(NBM,NGRO),SIGT3(NBM,NGRO)
      CHARACTER CDOOR*12,TITR*72
      LOGICAL LEAKSW,START,NOCONV(NBM,NGRO),BIEFF,LGC
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) KPLIB
      CHARACTER HSMG*131,CGRPNM*12
      LOGICAL LOGDO
      REAL FACT(NALPHA),SIGX(NALPHA)
      COMPLEX SUM
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IRES,ISONR,NPSYS
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR,GAS,SIGE,VST,DIST,FUN,DILG
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIG0,SIG3,TOTAL,DILUT,GC
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: PICX
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASKI
      COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: COEF,DENOM
      COMPLEX*16, ALLOCATABLE, DIMENSION(:,:,:) :: XCOEF,XDENO
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  DATA STATEMENTS
*----
      DATA FACT/0.01,0.03162278,0.1,0.3162278,1.0,3.162278,10.0,
     1 31.62278,100.0/
*----
*  SCRATCH STORAGE ALLOCATION
*   SIG0    macroscopic xs of the resonant isotopes as interpolated.
*   SIG3    macroscopic transport correction.
*----
      ALLOCATE(IRES(NBM),ISONR(NBNRS),NPSYS(NGRO))
      ALLOCATE(SIG0(NBM,NGRO),SIG3(NBM,NGRO),TOTAL(NGRO,NBNRS),
     1 GAR(NGRO),GAS(NGRO),SIGE(NGRO),DILUT(NALPHA,NGRO),
     2 GC(NGRO,NBNRS),VST(NBNRS),DIST(NBNRS))
      ALLOCATE(PICX(NALPHA,NBNRS,NGRO))
      ALLOCATE(MASKI(NBISO))
      ALLOCATE(COEF(NRAT,NGRO),DENOM(NRAT,NGRO))
      ALLOCATE(XCOEF(NRAT,NBNRS,NGRO),XDENO(NRAT,NBNRS,NGRO))
      ALLOCATE(IPISO(NBISO))
*----
*  FIND THE RESONANT MIXTURE NUMBERS AND THE CORRELATED ISOTOPES
*  ASSOCIATED WITH REGION INRS
*----
      IRS=0
      DO 30 IBM=1,NBM
      LOGDO=.FALSE.
      DO 10 I=1,NREG
      LOGDO=LOGDO.OR.(MAT(I).EQ.IBM)
   10 CONTINUE
      IF(.NOT.LOGDO) GO TO 30
      DO 20 ISO=1,NBISO
      IF((MIX(ISO).EQ.IBM).AND.(LSHI(ISO).EQ.INRS)) THEN
         IRS=IRS+1
         ISONR(IRS)=ISO
         GO TO 30
      ENDIF
   20 CONTINUE
   30 CONTINUE
      IF(IRS.NE.NBNRS) CALL XABORT('SHISN3: INVALID VALUE OF NBNRS.')
      CALL XDISET(IRES,NBM,0)
      DO 40 IRS=1,NBNRS
      ISO=ISONR(IRS)
      IRES(MIX(ISO))=IRS
   40 CONTINUE
*----
*  SET THE LCM MICROLIB ISOTOPEWISE DIRECTORIES.
*----
      CALL LIBIPS(IPLIB,NBISO,IPISO)
*----
*  UNLOAD MICROSCOPIC X-S FROM LCM TO SCRATCH STORAGE. SET THE
*  GOLDSTEIN-COHEN TO ONE IN LEVEL 2 CALCULATIONS.
*----
      DO 50 IRS=1,NBNRS
      ISO=ISONR(IRS)
      KPLIB=IPISO(ISO) ! set ISO-th isotope
      CALL LCMGET(KPLIB,'NTOT0',TOTAL(1,IRS))
      CALL XDRSET(GC(1,IRS),NGRO,1.0)
   50 CONTINUE
*
      CALL XDRSET(VST,NBNRS,0.0)
      DO 60 I=1,NREG
      IF(MAT(I).EQ.0) GO TO 60
      IND=IRES(MAT(I))
      IF(IND.GT.0) VST(IND)=VST(IND)+VOL(I)
   60 CONTINUE
*
      CALL XDISET(NPSYS,NGRO,0)
      DO 110 LLL=1,NGRO
      LOGDO=.FALSE.
      DO 70 IRS=1,NBNRS
      LOGDO=LOGDO.OR.NOCONV(MIX(ISONR(IRS)),LLL)
   70 CONTINUE
      IF(LOGDO) THEN
         NPSYS(LLL)=LLL
*
*        COMPUTE THE LIGHT AND RESONANT COMPONENTS OF THE MACROSCOPIC
*        CROSS SECTIONS IN EACH RESONANT MIXTURE.
         DO 80 IRS=1,NBNRS
         ISO=ISONR(IRS)
         IBM=MIX(ISO)
         SIGT(IBM,LLL)=SIGT(IBM,LLL)-TOTAL(LLL,IRS)*DEN(ISO)
   80    CONTINUE
         DO 90 IBM=1,NBM
         SIG0(IBM,LLL)=0.0
         SIG3(IBM,LLL)=SIGT3(IBM,LLL)
   90    CONTINUE
         DO 100 IRS=1,NBNRS
         ISO=ISONR(IRS)
         SIG0(MIX(ISO),LLL)=TOTAL(LLL,IRS)*DEN(ISO)
         SIG3(MIX(ISO),LLL)=0.0
  100    CONTINUE
         IF(IMPX.GE.10) THEN
            WRITE (6,400) LLL,(SIG0(I,LLL),I=1,NBM)
            WRITE (6,410) LLL,(SIGT(I,LLL),I=1,NBM)
            WRITE (6,420) LLL,(SIGT3(I,LLL),I=1,NBM)
         ENDIF
      ENDIF
  110 CONTINUE
*----
*  SET UP VECTORS SIGE AND SB.
*----
      CALL LCMSIX(IPLIB,'SHIBA',1)
*
      CALL XDRSET(SIGE,NGRO,0.0)
      ALLOCATE(FUN(NUN*NGRO))
      CALL LCMSIX(IPLIB,'--AVERAGE--',1)
      CALL SHIDST(IPLIB,NPSYS,IPTRK,IFTRAK,CDOOR,IMPX,NBM,NREG,NUN,
     1 NGRO,IPHASE,MAT,VOL,KEYFLX,LEAKSW,IRES,SIG0,SIGT,SIGT3,TITR,
     2 FUN,SIGE)
      CALL LCMSIX(IPLIB,' ',2)
      DO 130 LLL=1,NGRO
      IF(NPSYS(LLL).NE.0) THEN
         DO 120 I=1,NREG
         IBM=MAT(I)
         IF(IBM.EQ.0) GO TO 120
         IND=IRES(IBM)
         IF(IND.GT.0) THEN
            IOF=(LLL-1)*NUN+KEYFLX(I)
            ISO=ISONR(IND)
            IF(NOCONV(IBM,LLL)) SB(LLL,ISO)=FUN(IOF)/
     1      SIG0(IBM,LLL)
         ENDIF
  120    CONTINUE
      ENDIF
  130 CONTINUE
      DEALLOCATE(FUN)
*----
*  SET UP VECTORS DILUT AND SIGX.
*----
      IF(START) THEN
*        USE A VERY CHEAP APPROXIMATION TO START ITERATIONS.
         DO 145 LLL=1,NGRO
         DO 140 IALP=1,NALPHA
         SIGX(IALP)=0.0
         DILUT(IALP,LLL)=SIGE(LLL)
  140    CONTINUE
  145    CONTINUE
      ELSE
         AVDEN=0.0
         VOLTOT=0.0
         DO 150 IRS=1,NBNRS
         AVDEN=AVDEN+DEN(ISONR(IRS))*VST(IRS)
         VOLTOT=VOLTOT+VST(IRS)
  150    CONTINUE
         AVDEN=AVDEN/VOLTOT
         DO 160 IRS=1,NBNRS
         ISO=ISONR(IRS)
         DIST(IRS)=DEN(ISO)/AVDEN
  160    CONTINUE
         DO 220 IALP=1,NALPHA
         DO 190 LLL=1,NGRO
         IF(NPSYS(LLL).NE.0) THEN
            SIGX(IALP)=FACT(IALP)*SIGE(LLL)
            DO 170 IBM=1,NBM
            SIG0(IBM,LLL)=0.0
            SIG3(IBM,LLL)=SIGT3(IBM,LLL)
  170       CONTINUE
            DO 180 IRS=1,NBNRS
            ISO=ISONR(IRS)
            SIG0(MIX(ISO),LLL)=SIGX(IALP)*DIST(IRS)
            SIG3(MIX(ISO),LLL)=0.0
  180       CONTINUE
         ENDIF
  190    CONTINUE
         ALLOCATE(DILG(NGRO),FUN(NUN*NGRO))
         WRITE(CGRPNM,'(8H--BAND--,I4.4)') IALP
         CALL LCMSIX(IPLIB,CGRPNM,1)
         CALL SHIDST(IPLIB,NPSYS,IPTRK,IFTRAK,CDOOR,IMPX,NBM,NREG,NUN,
     1   NGRO,IPHASE,MAT,VOL,KEYFLX,LEAKSW,IRES,SIG0,SIGT,SIG3,TITR,
     2   FUN,DILG)
         CALL LCMSIX(IPLIB,' ',2)
         DO 210 LLL=1,NGRO
         IF(NPSYS(LLL).NE.0) THEN
            DILUT(IALP,LLL)=DILG(LLL)
            DO 200 I=1,NREG
            IBM=MAT(I)
            IF(IBM.EQ.0) GO TO 200
            IND=IRES(IBM)
            IF(IND.GT.0) THEN
               IOF=(LLL-1)*NUN+KEYFLX(I)
               PICX(IALP,IND,LLL)=FUN(IOF)/SIG0(IBM,LLL)
            ENDIF
  200       CONTINUE
         ENDIF
  210    CONTINUE
         DEALLOCATE(FUN,DILG)
  220    CONTINUE
      ENDIF
      CALL LCMSIX(IPLIB,' ',2)
*----
*  COMPUTE AVERAGE MACROSCOPIC DILUTION X-S USING A N-TERM RATIONAL
*  APPROXIMATION.
*----
      DO 260 LLL=1,NGRO
      IF(NPSYS(LLL).NE.0) THEN
         DO 230 IALP=1,NALPHA
         SIGX(IALP)=FACT(IALP)*SIGE(LLL)
  230    CONTINUE
*        **********************************************************
         CALL SHIRAT(IMPX,NRAT,SIGX,DILUT(1,LLL),LLL,A,COEF(1,LLL),
     1   DENOM(1,LLL))
*        **********************************************************
*
*        COMPUTE THE PIC BASE POINTS FOR A N-TERM RATIONAL APPROXIMATION
         IF(START) THEN
            DO 245 IRS=1,NBNRS
            DO 240 I=1,NRAT
            XCOEF(I,IRS,LLL)=COEF(I,LLL)
            XDENO(I,IRS,LLL)=DENOM(I,LLL)
  240       CONTINUE
  245       CONTINUE
         ELSE
            CALL SHIDIL(NRAT,NALPHA,NBNRS,COEF(1,LLL),DENOM(1,LLL),
     1      DILUT(1,LLL),PICX(1,1,LLL),SIGX,DIST,VST,IMPX,LLL,
     2      XCOEF(1,1,LLL),XDENO(1,1,LLL))
         ENDIF
         IF(.NOT.START.AND.BIEFF) THEN
            DO 255 IRS=1,NBNRS
            ISO=ISONR(IRS)
            IF(NOCONV(MIX(ISO),LLL)) THEN
               SIGRES=TOTAL(LLL,IRS)*DEN(ISO)
               IF(NBNRS.EQ.1) THEN
                  SUM=0.0
                  DO 250 I=1,NRAT
                  SUM=SUM+COEF(I,LLL)/(SIGRES+DENOM(I,LLL))
  250             CONTINUE
                  PX0=REAL(SUM)
               ELSE
                  PX0=SB(LLL,ISO)
               ENDIF
               SB(LLL,ISO)=(1.0/PX0-SIGRES)/DEN(ISO)
               IF(SB(LLL,ISO).LT.0.0) THEN
                  WRITE (HSMG,515) (ISONAM(I0,ISO),I0=1,3),SB(LLL,ISO),
     1            LLL
                  CALL XABORT(HSMG)
               ENDIF
            ENDIF
  255       CONTINUE
         ENDIF
      ENDIF
  260 CONTINUE
*----
*  APPLY A GOLDSTEIN-COHEN CORRECTION SIMILAR TO THE CORRECTION USED
*  IN SHISN2.
*----
      IF((.NOT.START).AND.LGC) THEN
         DO 295 IRS=1,NBNRS
         ISO=ISONR(IRS)
         KPLIB=IPISO(ISO) ! set ISO-th isotope
         CALL LCMLEN(KPLIB,'NGOLD',LENGT,ITYLCM)
         IF(LENGT.EQ.NGRO) THEN
            IF(IMPX.GE.5) WRITE (6,390) (ISONAM(I0,ISO),I0=1,3)
            CALL LCMGET(KPLIB,'NGOLD',GC(1,IRS))
         ENDIF
         IF(LEVEL.EQ.2) GO TO 295
         DO 290 JSO=1,NBISO
         IF((MIX(JSO).EQ.MIX(ISO)).AND.(JSO.NE.ISO).AND.
     1   (LSHI(JSO).NE.0)) THEN
            KPLIB=IPISO(JSO) ! set JSO-th isotope
            CALL LCMLEN(KPLIB,'NGOLD',LENGT,ITYLCM)
            IF(LENGT.EQ.NGRO) THEN
               CALL LCMGET(KPLIB,'SIGS00',GAS)
               CALL LCMGET(KPLIB,'NGOLD',GAR)
               DO 280 LLL=1,NGRO
               IF((NOCONV(MIX(JSO),LLL)).AND.(GAR(LLL).NE.1.0)) THEN
                  DDD=(1.0-GAR(LLL))*GAS(LLL)*DEN(JSO)
                  DO 270 I=1,NRAT
                  XDENO(I,IRS,LLL)=XDENO(I,IRS,LLL)-DDD
  270             CONTINUE
               ENDIF
  280          CONTINUE
            ENDIF
         ENDIF
  290    CONTINUE
  295    CONTINUE
      ENDIF
*----
*  COMPUTE THE DILUTION PARAMETERS (BARN) FOR EACH RESONANT ISOTOPE IN
*  RESONANT MIXTURE INRS
*----
      CALL SHIEQU(IPLIB,LEVEL,NGRO,NBISO,NBM,NBNRS,NRAT,MIX,ISONAM,
     1 NOCONV,ISONR,GC,COEF,DENOM,XCOEF,XDENO,DEN,IMPX,SN)
*
      DO 320 LLL=1,NGRO
      LOGDO=.FALSE.
      DO 300 IRS=1,NBNRS
      LOGDO=LOGDO.OR.NOCONV(MIX(ISONR(IRS)),LLL)
  300 CONTINUE
      IF(LOGDO) THEN
         DO 310 IRS=1,NBNRS
         ISO=ISONR(IRS)
         IBM=MIX(ISO)
         IF(START.OR.(.NOT.BIEFF)) THEN
            SB(LLL,ISO)=SN(LLL,ISO)
         ELSE IF(SB(LLL,ISO).LT.0.97*SN(LLL,ISO)) THEN
            WRITE (6,520)(ISONAM(I0,ISO),I0=1,3),SB(LLL,ISO)/
     1      SN(LLL,ISO),0.97,LLL
            SB(LLL,ISO)=0.97*SN(LLL,ISO)
         ENDIF
         SIGT(IBM,LLL)=SIGT(IBM,LLL)+TOTAL(LLL,IRS)*DEN(ISO)
  310    CONTINUE
      ENDIF
  320 CONTINUE
*----
*  SAVE THE GROUP- AND ISOTOPE-DEPENDENT DILUTIONS
*----
      CALL LCMPUT(IPLIB,'ISOTOPESDSB',NBISO*NGRO,2,SB)
      CALL LCMPUT(IPLIB,'ISOTOPESDSN',NBISO*NGRO,2,SN)
*----
*  COMPUTE THE SELF-SHIELDED MICROSCOPIC CROSS SECTIONS AND UPDATE
*  VECTOR SIGT
*----
      DO 330 ISO=1,NBISO
      LOGDO=START.OR.(DEN(ISO).NE.0.)
      MASKI(ISO)=(LSHI(ISO).EQ.INRS).AND.LOGDO
  330 CONTINUE
      IMPX2=MAX(0,IMPX-1)
      CALL LIBLIB (IPLIB,NBISO,MASKI,IMPX2)
      DO 345 IRS=1,NBNRS
      ISO=ISONR(IRS)
      IBM=MIX(ISO)
      KPLIB=IPISO(ISO) ! set ISO-th isotope
      CALL LCMGET(KPLIB,'NTOT0',GAR)
      DO 340 LLL=1,NGRO
      TOTAL(LLL,IRS)=TOTAL(LLL,IRS)-GAR(LLL)
      SIGT(IBM,LLL)=SIGT(IBM,LLL)-DEN(ISO)*TOTAL(LLL,IRS)
  340 CONTINUE
  345 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IPISO)
      DEALLOCATE(XDENO,XCOEF)
      DEALLOCATE(DENOM,COEF)
      DEALLOCATE(MASKI)
      DEALLOCATE(PICX)
      DEALLOCATE(DIST,VST,GC,DILUT,SIGE,GAS,GAR,TOTAL,SIG3,SIG0)
      DEALLOCATE(NPSYS,ISONR,IRES)
      RETURN
*
  390 FORMAT(/53H SHISN3: GOLDSTEIN AND COHEN APPROXIMATION USED FOR I,
     1 8HSOTOPE ',3A4,2H'.)
  400 FORMAT(//51H TOTAL MACROSCOPIC CROSS SECTIONS OF THE RESONANT M,
     1 31HATERIALS IN EACH MIXTURE (GROUP,I5,2H):/(1X,1P,11E11.3))
  410 FORMAT(//51H TOTAL MACROSCOPIC CROSS SECTIONS OF THE OTHER MATE,
     1 28HRIALS IN EACH MIXTURE (GROUP,I5,2H):/(1X,1P,11E11.3))
  420 FORMAT(//1X,'TRANSPORT CORRECTION CROSS SECTIONS OF THE OTHER ',
     1'MATERIALS IN EACH MIXTURE (GROUP',I5,'):'/(1X,1P,11E11.3))
  515 FORMAT(30HSHISN3: THE RESONANT ISOTOPE ',3A4,14H' HAS A NEGATI,
     1 22HVE L-J CROSS-SECTION (,1P,E14.4,0P,10H) IN GROUP,I4,1H.)
  520 FORMAT(54H SHISN3: THE L-J EQUIVALENCE FACTOR OF RESONANT ISOTOP,
     1 3HE ',3A4,18H' WAS CHANGED FROM,F6.3,3H TO,F5.2,9H IN GROUP,I4,
     2 1H.)
      END
