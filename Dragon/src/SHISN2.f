*DECK SHISN2
      SUBROUTINE SHISN2 (IPLIB,IPTRK,IFTRAK,NGRO,NBISO,NBM,NREG,NUN,
     1 CDOOR,NRES,NBM2,IMPX,ISONAM,MIX,DEN,SN,SB,LSHI,IPHASE,MAT,VOL,
     2 KEYFLX,LEAKSW,TITR,START,SIGT,SIGT3,NOCONV,BIEFF,LGC,SIGE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one multidimensional self-shielding iteration using the
* generalized Stamm'ler algorithm without Nordheim (PIC) approximation.
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
* IPLIB   pointer to the internal microscopic cross section library
*         (L_LIBRARY signature).
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  unit number of the sequential binary tracking file.
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NBM     number of mixtures in the macrolib.
* NREG    number of volumes.
* NUN     number of unknowns in the flux or source vector in one
*         energy group.
* CDOOR   name of the geometry/solution module.
* NRES    number of resonant mixtures.
* NBM2    number of resonant isotopes.
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
* START   beginning-of-iteration flag (.TRUE. if SHISN2 is called
*         for the first time).
* SIGT3   transport correction.
* NOCONV  mixture convergence flag (.TRUE. if mixture IBM
*         is not converged in group L).
* BIEFF   Livolant-Jeanpierre normalization flag (.TRUE. to
*         activate).
* LGC     Goldstein-Cohen approximation flag (.TRUE. to activate).
*
*Parameters: output
* SN      on input, estimate of the dilution cross section in each 
*         energy group of each isotope. A value of 1.0e10 is used 
*         for infinite dilution.
*         On output, computed dilution cross section in each energy 
*         group of each isotope.
* SIGT    total macroscopic cross sections as modified by Shiba.
*
*Parameters: output
* SB      dilution cross section as used in Livolant-Jeanpierre
*         normalization.
* SIGE    computed macroscopic dilution cross section in each resonant
*         mixture and each energy group.
*
*Reference:
* A. Hebert and G. Marleau, Generalization of the Stamm'ler Method
* for the Self-Shielding of Resonant isotopes in Arbitrary Geometries,
* Nucl. Sci. Eng. 108, 230 (1991).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      PARAMETER (NALPHA=5)
      TYPE(C_PTR) IPLIB,IPTRK
      INTEGER IFTRAK,NGRO,NBISO,NBM,NREG,NUN,NRES,NBM2,IMPX,
     1 ISONAM(3,NBISO),MIX(NBISO),LSHI(NBISO),IPHASE,MAT(NREG),
     2 KEYFLX(NREG)
      REAL DEN(NBISO),SN(NGRO,NBISO),SB(NGRO,NBISO),VOL(NREG),
     1 SIGT(NBM,NGRO),SIGT3(NBM,NGRO),SIGE(NRES,NGRO)
      CHARACTER CDOOR*12,TITR*72,CGRPNM*12
      LOGICAL LEAKSW,START,NOCONV(NBM,NGRO),BIEFF,LGC
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131
      LOGICAL LOGDO
      COMPLEX COEF(3),DENOM(3),EAV
      PARAMETER (NRAT=(NALPHA+1)/2)
      TYPE(C_PTR) KPLIB
      REAL FACT(NALPHA),SIGX(NALPHA)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IRES,MIX2,IRNBM,NPSYS
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR,SIGRES,DILAV,FUN
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIG0,SIG1,SIG3,TOTAL,SIGOLD,
     1 DILUT
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASKI
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  DATA STATEMENTS
*----
      DATA FACT/0.01,0.1,1.0,10.0,100.0/
*----
*  SCRATCH STORAGE ALLOCATION
*   SIG0    macroscopic xs of the resonant isotopes as interpolated.
*   SIG1    macroscopic xs of the resonant isotopes at various SIGX.
*   SIG3    macroscopic transport correction.
*----
      ALLOCATE(IRES(NBM),MIX2(NBISO),IRNBM(NBM),NPSYS(NGRO))
      ALLOCATE(SIG0(NBM,NGRO),SIG1(NBM,NGRO),SIG3(NBM,NGRO),
     1 TOTAL(NGRO,NBM2),SIGOLD(NGRO,NBM2),GAR(NGRO),SIGRES(NBM),
     2 DILAV(NGRO),DILUT(NALPHA,NGRO))
      ALLOCATE(MASKI(NBISO))
      ALLOCATE(IPISO(NBISO))
*----
*  SET THE LCM MICROLIB ISOTOPEWISE DIRECTORIES.
*----
      CALL LIBIPS(IPLIB,NBISO,IPISO)
*----
*  UNLOAD MICROSCOPIC X-S FROM LCM TO SCRATCH STORAGE
*----
      IBM=0
      DO 20 ISO=1,NBISO
      MIX2(ISO)=0
      IF(LSHI(ISO).GT.0) THEN
         IBM=IBM+1
         MIX2(ISO)=IBM
         KPLIB=IPISO(ISO) ! set ISO-th isotope
         CALL LCMGET(KPLIB,'NTOT0',TOTAL(1,IBM))
         CALL LCMLEN(KPLIB,'NGOLD',LENGT,ITYLCM)
         IF((LENGT.EQ.NGRO).AND.(.NOT.START).AND.LGC) THEN
            IF(IMPX.GE.5) WRITE (6,390) (ISONAM(I0,ISO),I0=1,3)
            CALL LCMGET(KPLIB,'SIGS00',SIGOLD(1,IBM))
            CALL LCMGET(KPLIB,'NGOLD',GAR)
            DO 10 LLL=1,NGRO
            SIGOLD(LLL,IBM)=(1.0-GAR(LLL))*SIGOLD(LLL,IBM)
   10       CONTINUE
         ELSE
            CALL XDRSET(SIGOLD(1,IBM),NGRO,0.0)
         ENDIF
      ENDIF
   20 CONTINUE
*----
*  LOOP OVER RESONANT REGIONS. THE CP ARE STORED ON DIRECTORY SHIBA
*----
      CALL LCMSIX(IPLIB,'SHIBA',1)
      DO 260 INRS=1,NRES
*----
*  FIND THE RESONANT MIXTURE NUMBERS (IRNBM) ASSOCIATED WITH REGION INRS
*----
      NBNRS=0
      DO 50 IBM=1,NBM
      IRES(IBM)=0
      DO 40 ISO=1,NBISO
      IF((MIX(ISO).EQ.IBM).AND.(LSHI(ISO).EQ.INRS)) THEN
         NBNRS=NBNRS+1
         IRNBM(NBNRS)=IBM
         IRES(IBM)=1
         GO TO 50
      ENDIF
   40 CONTINUE
   50 CONTINUE
      IF(NBNRS.EQ.0) THEN
         IF(START.AND.(IMPX.GE.1)) WRITE(6,385) 'SHISN2',INRS
         GO TO 260
      ELSE IF(START.AND.(NBNRS.GT.1).AND.(IMPX.GE.5)) THEN
         WRITE (6,380) NBNRS,INRS
      ENDIF
*
      CALL XDISET(NPSYS,NGRO,0)
      DO 120 LLL=1,NGRO
      LOGDO=.FALSE.
      DO 60 I=1,NBNRS
      LOGDO=LOGDO.OR.NOCONV(IRNBM(I),LLL)
   60 CONTINUE
      IF(LOGDO) THEN
         NPSYS(LLL)=LLL
*
*        COMPUTE THE LIGHT AND RESONANT COMPONENTS OF THE MACROSCOPIC
*        CROSS SECTIONS IN EACH RESONANT MIXTURE.
         DO 80 I=1,NBNRS
         SIGRES(I)=0.0
         DO 70 ISO=1,NBISO
         IF((MIX(ISO).EQ.IRNBM(I)).AND.(LSHI(ISO).EQ.INRS)) THEN
            SIGRES(I)=SIGRES(I)+TOTAL(LLL,MIX2(ISO))*DEN(ISO)
         ENDIF
   70    CONTINUE
         SIGT(IRNBM(I),LLL)=SIGT(IRNBM(I),LLL)-SIGRES(I)
   80    CONTINUE
         DO 90 IBM=1,NBM
         SIG0(IBM,LLL)=0.0
         SIG1(IBM,LLL)=0.0
         SIG3(IBM,LLL)=SIGT3(IBM,LLL)
   90    CONTINUE
         DO 110 I=1,NBNRS
         SIG0(IRNBM(I),LLL)=SIGRES(I)
         SIG3(IRNBM(I),LLL)=0.0
  110    CONTINUE
         IF(IMPX.GE.10) THEN
            WRITE (6,400) LLL,(SIG0(I,LLL),I=1,NBM)
            WRITE (6,410) LLL,(SIGT(I,LLL),I=1,NBM)
            WRITE (6,420) LLL,(SIGT3(I,LLL),I=1,NBM)
         ENDIF
      ENDIF
  120 CONTINUE
*----
*  SET UP VECTORS DILUT AND SIGX.
*----
      CALL XDRSET (DILAV,NGRO,0.0)
      IF(START) THEN
*        USE A VERY CHEAP APPROXIMATION TO START ITERATIONS.
         ALLOCATE(FUN(NUN*NGRO))
         CALL LCMSIX(IPLIB,'--AVERAGE--',1)
         CALL SHIDST (IPLIB,NPSYS,IPTRK,IFTRAK,CDOOR,IMPX,NBM,NREG,
     1   NUN,NGRO,IPHASE,MAT,VOL,KEYFLX,LEAKSW,IRES,SIG0,SIGT,SIGT3,
     2   TITR,FUN,DILAV)
         CALL LCMSIX(IPLIB,' ',2)
         DEALLOCATE(FUN)
         DO 135 LLL=1,NGRO
         DO 130 IALP=1,NALPHA
         DILUT(IALP,LLL)=DILAV(LLL)
  130    CONTINUE
  135    CONTINUE
      ELSE
         DO 165 IALP=1,NALPHA
         DO 150 LLL=1,NGRO
         IF(NPSYS(LLL).NE.0) THEN
            DO 140 I=1,NBNRS
            SIG1(IRNBM(I),LLL)=FACT(IALP)*SIGE(INRS,LLL)
  140       CONTINUE
         ENDIF
  150    CONTINUE
         ALLOCATE(FUN(NUN*NGRO))
         WRITE(CGRPNM,'(8H--BAND--,I4.4)') IALP
         CALL LCMSIX(IPLIB,CGRPNM,1)
         CALL SHIDST (IPLIB,NPSYS,IPTRK,IFTRAK,CDOOR,IMPX,NBM,NREG,
     1   NUN,NGRO,IPHASE,MAT,VOL,KEYFLX,LEAKSW,IRES,SIG1,SIGT,SIG3,
     2   TITR,FUN,DILAV)
         CALL LCMSIX(IPLIB,' ',2)
         DEALLOCATE(FUN)
         DO 160 LLL=1,NGRO
         DILUT(IALP,LLL)=DILAV(LLL)
  160    CONTINUE
  165    CONTINUE
      ENDIF
*----
*  COMPUTE AVERAGE MACROSCOPIC DILUTION X-S (SIGE) USING A THREE-TERM
*  RATIONAL APPROXIMATION.
*----
      DO 200 LLL=1,NGRO
      IF(NPSYS(LLL).NE.0) THEN
         DO 170 IALP=1,NALPHA
         SIGX(IALP)=FACT(IALP)*SIGE(INRS,LLL)
  170    CONTINUE
         IMPX2=IMPX
         IF(START) IMPX2=MAX(0,IMPX-10)
*        **********************************************************
         CALL SHIRAT(IMPX2,NRAT,SIGX,DILUT(1,LLL),LLL,A,COEF,DENOM)
*        **********************************************************
         EAV=(COEF(1)*SQRT(DENOM(1))+COEF(2)*SQRT(DENOM(2))+
     1   COEF(3)*SQRT(DENOM(3)))**2
         SIGE(INRS,LLL)=REAL(EAV)
         IF((.NOT.START).AND.(BIEFF).AND.(NBNRS.EQ.1)) THEN
*           COMPUTE DILAV FOR THE L-J NORMALIZATION.
            SIGXX=SIG0(IRNBM(1),LLL)
            PXX=REAL(COEF(1)/(SIGXX+DENOM(1))+COEF(2)/(SIGXX+DENOM(2))
     1      +COEF(3)/(SIGXX+DENOM(3)))
            DILAV(LLL)=1.0/PXX-SIGXX
         ENDIF
*----
*  COMPUTE THE ISOTOPE DILUTION MICROSCOPIC CROSS SECTIONS (SN) USED
*  FOR LIBRARY INTERPOLATION.
*----
         DO 190 ISO=1,NBISO
         IF((LSHI(ISO).EQ.INRS).AND.(IRES(MIX(ISO)).EQ.1).AND.
     1   (DEN(ISO).NE.0.)) THEN
            SUM=0.0
            DO 180 JSO=1,NBISO
            IBM=MIX(JSO)
            IF((LSHI(JSO).EQ.INRS).AND.(IBM.EQ.MIX(ISO)).AND.
     1      (ISO.NE.JSO)) SUM=SUM+(TOTAL(LLL,MIX2(JSO))-
     2      SIGOLD(LLL,MIX2(JSO)))*DEN(JSO)
  180       CONTINUE
            SN(LLL,ISO)=REAL(COEF(1)*SQRT(DENOM(1)+SUM)+COEF(2)*
     1      SQRT(DENOM(2)+SUM)+COEF(3)*SQRT(DENOM(3)+SUM))**2/DEN(ISO)
            IF(SN(LLL,ISO).LE.0.0) THEN
               WRITE (HSMG,510) (ISONAM(I0,ISO),I0=1,3),SN(LLL,ISO),LLL
               CALL XABORT(HSMG)
            ENDIF
         ELSE IF((LSHI(ISO).EQ.INRS).AND.(IRES(MIX(ISO)).EQ.1).AND.
     1   (DEN(ISO).EQ.0.)) THEN
            SN(LLL,ISO)=1.0E10
         ENDIF
  190    CONTINUE
         IF((.NOT.START).AND.(IMPX.GE.10)) THEN
            DO 195 I=1,NBNRS
            PP=A-SIGT(IRNBM(I),LLL)
            QQ=SIGE(INRS,LLL)-SIGT(IRNBM(I),LLL)
            IF(ABS(PP).GT.1.0E-4*SIGT(IRNBM(I),LLL)) THEN
               BEL=QQ/PP
            ELSE
               BEL=0.0
            ENDIF
            WRITE (6,610) I,SIGE(INRS,LLL),BEL
  195       CONTINUE
         ENDIF
      ENDIF
  200 CONTINUE
*----
*  COMPUTE THE ISOTOPE DILUTION MICROSCOPIC CROSS SECTIONS (SB) USED
*  IN L-J NORMALIZATION.
*----
      IF((.NOT.START).AND.(BIEFF).AND.(NBNRS.GT.1)) THEN
*        COMPUTE DILAV FOR THE L-J NORMALIZATION.
         ALLOCATE(FUN(NUN*NGRO))
         CALL LCMSIX(IPLIB,'--AVERAGE--',1)
         CALL SHIDST (IPLIB,NPSYS,IPTRK,IFTRAK,CDOOR,IMPX,NBM,NREG,
     1   NUN,NGRO,IPHASE,MAT,VOL,KEYFLX,LEAKSW,IRES,SIG0,SIGT,SIGT3,
     2   TITR,FUN,DILAV)
         CALL LCMSIX(IPLIB,' ',2)
         DEALLOCATE(FUN)
      ENDIF
      DO 250 LLL=1,NGRO
      IF(NPSYS(LLL).NE.0) THEN
         DO 220 ISO=1,NBISO
         IF((LSHI(ISO).EQ.INRS).AND.(IRES(MIX(ISO)).EQ.1).AND.
     1   (DEN(ISO).NE.0.)) THEN
            SUM=0.0
            DO 210 JSO=1,NBISO
            IBM=MIX(JSO)
            IF((LSHI(JSO).EQ.INRS).AND.(IBM.EQ.MIX(ISO)).AND.
     1      (ISO.NE.JSO)) SUM=SUM+TOTAL(LLL,MIX2(JSO))*DEN(JSO)
  210       CONTINUE
            IF(START.OR.(.NOT.BIEFF)) THEN
               SB(LLL,ISO)=SN(LLL,ISO)
            ELSE
               SB(LLL,ISO)=(DILAV(LLL)+SUM)/DEN(ISO)
               IF(SB(LLL,ISO).LT.0.0) THEN
                  WRITE (HSMG,515) (ISONAM(I0,ISO),I0=1,3),SB(LLL,ISO),
     1            LLL
                  CALL XABORT(HSMG)
               ELSE IF(SB(LLL,ISO).LT.SN(LLL,ISO)) THEN
                  IF(SB(LLL,ISO).LT.0.99*SN(LLL,ISO)) WRITE (6,520)
     1            (ISONAM(I0,ISO),I0=1,3),SB(LLL,ISO)/SN(LLL,ISO),LLL
                  SB(LLL,ISO)=SN(LLL,ISO)
               ENDIF
            ENDIF
         ELSE IF((LSHI(ISO).EQ.INRS).AND.(IRES(MIX(ISO)).EQ.1).AND.
     1   (DEN(ISO).EQ.0.)) THEN
            SB(LLL,ISO)=1.0E10
         ENDIF
  220    CONTINUE
*
*        RESTORE SIGT ARRAY.
         DO 240 I=1,NBNRS
         SIGRES(I)=0.0
         DO 230 ISO=1,NBISO
         IF((MIX(ISO).EQ.IRNBM(I)).AND.(LSHI(ISO).EQ.INRS)) THEN
            SIGRES(I)=SIGRES(I)+TOTAL(LLL,MIX2(ISO))*DEN(ISO)
         ENDIF
  230    CONTINUE
         SIGT(IRNBM(I),LLL)=SIGT(IRNBM(I),LLL)+SIGRES(I)
  240    CONTINUE
      ENDIF
  250 CONTINUE
  260 CONTINUE
      CALL LCMSIX(IPLIB,' ',2)
*----
*  SAVE THE GROUP- AND ISOTOPE-DEPENDENT DILUTIONS
*----
      CALL LCMPUT(IPLIB,'ISOTOPESDSB',NBISO*NGRO,2,SB)
      CALL LCMPUT(IPLIB,'ISOTOPESDSN',NBISO*NGRO,2,SN)
*----
*  COMPUTE THE SELF-SHIELDED MICROSCOPIC CROSS SECTIONS AND UPDATE
*  VECTOR SIGT
*----
      DO 290 ISO=1,NBISO
      LOGDO=START.OR.(DEN(ISO).NE.0.)
      MASKI(ISO)=(LSHI(ISO).GT.0).AND.LOGDO
  290 CONTINUE
      IMPX2=MAX(0,IMPX-1)
      CALL LIBLIB (IPLIB,NBISO,MASKI,IMPX2)
      DO 320 ISO=1,NBISO
      IBM=MIX(ISO)
      IF((LSHI(ISO).GT.0).AND.(IBM.GT.0).AND.(DEN(ISO).NE.0.)) THEN
         KPLIB=IPISO(ISO) ! set ISO-th isotope
         CALL LCMGET(KPLIB,'NTOT0',GAR)
         DO 300 LLL=1,NGRO
         TOTAL(LLL,MIX2(ISO))=TOTAL(LLL,MIX2(ISO))-GAR(LLL)
  300    CONTINUE
         DO 310 LLL=1,NGRO
         IF(NOCONV(IBM,LLL)) SIGT(IBM,LLL)=SIGT(IBM,LLL)-DEN(ISO)*
     1   TOTAL(LLL,MIX2(ISO))
  310    CONTINUE
      ENDIF
  320 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IPISO)
      DEALLOCATE(MASKI)
      DEALLOCATE(DILUT,DILAV,SIGRES,GAR,SIGOLD,TOTAL,SIG3,SIG1,SIG0)
      DEALLOCATE(NPSYS,IRNBM,MIX2,IRES)
      RETURN
*
  380 FORMAT(/16H SHISN2: MERGING,I3,30H RESONANT MIXTURES IN RESONANT,
     1 14H REGION NUMBER,I3,1H.)
  385 FORMAT(A6,1X,': RESONANT REGION =',I10,1X,'NOT USED.')
  390 FORMAT(/53H SHISN2: GOLDSTEIN AND COHEN APPROXIMATION USED FOR I,
     1 8HSOTOPE ',3A4,2H'.)
  400 FORMAT(1X,'TOTAL MACROSCOPIC CROSS SECTIONS OF THE RESONANT ',
     1'MATERIALS IN EACH MIXTURE (GROUP',I5,'):'/(1X,1P,11E11.3))
  410 FORMAT(1X,'TOTAL MACROSCOPIC CROSS SECTIONS OF THE OTHER ',
     1'MATERIALS IN EACH MIXTURE (GROUP',I5,'):'/(1X,1P,11E11.3))
  420 FORMAT(//1X,'TRANSPORT CORRECTION CROSS SECTIONS OF THE OTHER ',
     1'MATERIALS IN EACH MIXTURE (GROUP',I5,'):'/(1X,1P,11E11.3))
  510 FORMAT(30HSHISN2: THE RESONANT ISOTOPE ',3A4,14H' HAS A NEGATI,
     1 27HVE DILUTION CROSS-SECTION (,1P,E14.4,0P,10H) IN GROUP,I4,1H.)
  515 FORMAT(30HSHISN2: THE RESONANT ISOTOPE ',3A4,14H' HAS A NEGATI,
     1 22HVE L-J CROSS-SECTION (,1P,E14.4,0P,10H) IN GROUP,I4,1H.)
  520 FORMAT(54H SHISN2: THE L-J EQUIVALENCE FACTOR OF RESONANT ISOTOP,
     1 3HE ',3A4,18H' WAS CHANGED FROM,F6.3,16H TO 1.0 IN GROUP,I4,1H.)
  610 FORMAT(8X,8HAVERAGE(,I2,1H),1P,E13.5/8X,11HBELL FACTOR,E13.5)
      END
