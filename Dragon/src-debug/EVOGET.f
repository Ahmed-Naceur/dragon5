*DECK EVOGET
      SUBROUTINE EVOGET(IPRINT,ITYPE ,ITIXS ,IEXTR ,IGLOB ,ISAT  ,
     >                  IDIRAC,ISAVE ,ISET  ,INR   ,IDEPL ,IFLMAC,
     >                  IYLMIX,RPAR  ,XT    ,NBMIX ,IPICK ,MIXBRN,
     >                  MIXPWR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read burnup input options.
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
* IPRINT  print flag.
* ITYPE   solution method:                         
*            ITYPE=1 for fifth order Cash-Karp algorithm;
*            ITYPE=2 for fourth order Kaps-Renthrop algorithm.
* ITIXS   flag for time dependent xs:
*            ITIXS= 0 for flag ON;
*            ITIXS= 1 for flag OFF.
* IEXTR   flag for power normalization: 
*            IEXTR= 0 for flag OFF;
*            IEXTR= 1 for flag ON.
* IGLOB   flag for power computation option: 
*            IGLOB= 0 for flag OFF;
*            IGLOB= 1 for flag ON.
* ISAT    flag for saturaton: 
*            ISAT= 0 for flag OFF;
*            ISAT= 1 for flag ON.
* IDIRAC  flag for delta Dirac saturation: 
*            IDIRAC= 0 for flag OFF;
*            IDIRAC= 1 for flag ON.
* ISAVE   flag for SAVE:
*            ISAVE=-1 for no SAVE. 
*            ISAVE= 0 for automatic SAVE.
*            ISAVE= 1 for manual SAVE.
* ISET    flag for SET:
*            ISET= 0 for automatic SET;
*            ISET= 1 for manual SET.
* INR     burnup option:
*            INR= 0 for out-of-core depletion;
*            INR= 1 for constant flux depletion;
*            INR= 2 for constant power per kg;
*            INR= 2 for constant power per volume (cc).
* IFLMAC  flag to recover neutron flux:
*            IFLMAC= 0 recover from L_FLUX object;
*            IFLMAC= 1 recover from embedded macrolib 
*            in L_LIBRARY object;
*            IFLMAC= 2 recover from 'FLUX-BUND' record 
*            in L_POWER object.
* IYLMIX  flag to recover fission yield data:
*            IYLMIX= 0 recover from DEPL-CHAIN data;
*            IYLMIX= 1 recover from isotopic PYIELD data.
* IDEPL   flag for depletion:
*            IDEPL= 0 for no depletion;
*            IDEPL= 1 for depletion.
* RPAR    burnup parameters:
*            RPAR(1) = EPS1   accuracy of ODE solver;
*            RPAR(2) = EPS2   accuracy for constant
*                             power iteration;
*            RPAR(3) = EXPMAX isotope saturation flag;
*            RPAR(4) = H1     guessed first time step;
*            RPAR(5) = FIT    flux (N/CM**2/S) OR
*                             power (kW/kG INITIAL HEAVY)
*                             W/CC  (W/CC)
*                             normalization factor.
* XT      time control table:
*            XT(1)   = initial time for depletion;
*            XT(2)   = final time for depletion;
*            XT(3)   = time for save;
*            XT(4)   = time for last set;
*            XT(5)   = time for current set.
* NBMIX   number of mixtures.
* IPICK   flag for burnup value recovery in a CLE-2000 variable:
*            IPICK= 0 for no recovery;
*            IPICK= 1 for recovery.
* MIXBRN  flags for mixtures to burn.
* MIXPWR  flags for mixtures to include in power normalization.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT         NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER          IPRINT,ITYPE,ITIXS,IEXTR,IGLOB,ISAT,
     >                 IDIRAC,ISAVE,ISET,INR,IDEPL,IFLMAC,IYLMIX,NBMIX,
     >                 IPICK,MIXBRN(NBMIX),MIXPWR(NBMIX)
      REAL             RPAR(5),XT(5)
*----
*  LOCAL VARIABLES
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      REAL             CSEC,CDAY,CYEAR
      PARAMETER       (IOUT=6,CSEC=1.0E-8,CDAY=8.64E-4,
     >                 CYEAR=3.1536E-1,NAMSBR='EVOGET')
*----
*  INPUT FILE PARAMETERS
*----
      CHARACTER        CARLIR*12
      INTEGER          ITYPLU,INTLIR
      REAL             REALIR
      DOUBLE PRECISION DBLLIR
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IRDT,KMIXB,KMIXP,IMIX
*----
*  READ THE BURNUP INPUT DATA.
*----
      ISAVE=0
      ISET=0
      KMIXB=0
      KMIXP=0
 100  CONTINUE
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
 101  CONTINUE
      IF(ITYPLU .EQ. 10) GO TO 105
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': CHARACTER KEYWORD EXPECTED.')
      IF(CARLIR(1:4) .EQ. ';') THEN
        IPICK=0
        GO TO 105
      ELSE IF(CARLIR(1:4) .EQ. 'PICK') THEN
        IPICK=1
        GO TO 105
      ELSE IF(CARLIR(1:4) .EQ. 'EDIT') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': INTEGER IPRINT EXPECTED')
        IPRINT=INTLIR
      ELSE IF(CARLIR(1:4) .EQ. 'EPS1') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': REAL EPS1 EXPECTED')
        RPAR(1)=REALIR
      ELSE IF(CARLIR(1:4) .EQ. 'EPS2') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': REAL EPS2 EXPECTED')
        RPAR(2)=REALIR
      ELSE IF(CARLIR(1:4) .EQ. 'EXPM') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': REAL EXPM EXPECTED')
        RPAR(3)=REALIR
      ELSE IF(CARLIR(1:6) .EQ. 'SATOFF') THEN
        RPAR(3)=0.0
      ELSE IF(CARLIR(1:9) .EQ. 'FLUX_FLUX') THEN
        IFLMAC=0
      ELSE IF(CARLIR(1:8) .EQ. 'FLUX_MAC') THEN
        IFLMAC=1
      ELSE IF(CARLIR(1:8) .EQ. 'FLUX_POW') THEN
        IFLMAC=2
      ELSE IF(CARLIR(1:8) .EQ. 'CHAIN') THEN
        IYLMIX=0
      ELSE IF(CARLIR(1:8) .EQ. 'PIFI') THEN
        IYLMIX=1
      ELSE IF(CARLIR(1:4) .EQ. 'H1') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  'REAL H1 EXPECTED')
        RPAR(4)=REALIR
      ELSE IF(CARLIR(1:4) .EQ. 'RUNG') THEN
        ITYPE=1
      ELSE IF(CARLIR(1:4) .EQ. 'KAPS') THEN
        ITYPE=2
      ELSE IF(CARLIR(1:4) .EQ. 'NOEX') THEN
        IEXTR=0
      ELSE IF(CARLIR(1:4) .EQ. 'EXTR') THEN
        IEXTR=1
      ELSE IF(CARLIR(1:4) .EQ. 'NOGL') THEN
        IGLOB=0
      ELSE IF(CARLIR(1:4) .EQ. 'GLOB') THEN
        IGLOB=1
      ELSE IF(CARLIR(1:4) .EQ. 'NSAT') THEN
        ISAT=0
      ELSE IF(CARLIR(1:4) .EQ. 'SAT ') THEN
        ISAT=1
      ELSE IF(CARLIR(1:4) .EQ. 'NODI') THEN
        IDIRAC=0
      ELSE IF(CARLIR(1:4) .EQ. 'DIRA') THEN
        IDIRAC=1
      ELSE IF(CARLIR(1:4) .EQ. 'TIXS') THEN
        ITIXS=1
      ELSE IF(CARLIR(1:4) .EQ. 'TDXS') THEN
        ITIXS=0
      ELSE IF(CARLIR(1:4) .EQ. 'MIXB') THEN
        KMIXB=1
        CALL XDISET(MIXBRN,NBMIX,0)
        DO IMIX=1,NBMIX
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) GO TO 101
          IF(INTLIR .GE. 1 .AND. INTLIR .LE. NBMIX) THEN
            MIXBRN(INTLIR)=1
          ENDIF
        ENDDO
      ELSE IF(CARLIR(1:4) .EQ. 'MIXP') THEN
        KMIXP=1
        CALL XDISET(MIXPWR,NBMIX,0)
        DO IMIX=1,NBMIX
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) GO TO 101
          IF(INTLIR .GE. 1 .AND. INTLIR .LE. NBMIX) THEN
            MIXPWR(INTLIR)=1
          ENDIF
        ENDDO
      ELSE IF(CARLIR(1:4) .EQ. 'NOSA') THEN
        ISAVE=-1
      ELSE IF(CARLIR(1:4) .EQ. 'SAVE') THEN
        ISAVE=1
*       SAVE THE LAST FLUX CALCULATION IN THE DEPLETION TABLE.
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': SAVE TIME EXPECTED')
        XT(3)=REALIR
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': TIME UNITS EXPECTED')
        IF(CARLIR(1:4) .EQ. 'S') THEN
          XT(3)=XT(3)*CSEC
        ELSE IF(CARLIR(1:4) .EQ. 'DAY') THEN
          XT(3)=XT(3)*CDAY
        ELSE IF(CARLIR(1:4) .EQ. 'YEAR') THEN
          XT(3)=XT(3)*CYEAR
        ELSE
          CALL XABORT(NAMSBR//': INVALID TIME UNITS')
        ENDIF
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': SAVE NORMALIZATION EXPECTED')
        IF(CARLIR(1:4) .EQ. 'FLUX') THEN
          INR=1
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >    ': BURNUP FLUX LEVEL EXPECTED')
          RPAR(5)=REALIR
          IF(RPAR(5) .EQ. 0.0) THEN
            INR=0
          ENDIF
        ELSE IF(CARLIR(1:4) .EQ. 'POWR') THEN
          INR=2
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >    ': BURNUP POWER EXPECTED')
          RPAR(5)=REALIR
          IF(RPAR(5) .EQ. 0.0) THEN
            INR=0
          ENDIF
        ELSE IF(CARLIR(1:4) .EQ. 'W/CC') THEN
          INR=3
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >    ': BURNUP POWER EXPECTED')
          RPAR(5)=REALIR
          IF(RPAR(5) .EQ. 0.0) THEN
            INR=0
          ENDIF
        ELSE IF(CARLIR(1:4) .EQ. 'KEEP') THEN
          INR=4
        ELSE
          GO TO 101
        ENDIF
      ELSE IF(CARLIR(1:4) .EQ. 'DEPL') THEN
        IDEPL=1
        IRDT=2
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': INITIAL OR INCREMENTAL DEPLETE TIME EXPECTED')
        XT(1)=REALIR
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .EQ. 1) THEN
          CALL XABORT(NAMSBR//': UNITS OR FINAL TIME EXPECTED')
        ELSE IF(ITYPLU .EQ. 2) THEN
          XT(2)=REALIR
          IRDT=1
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        ENDIF
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': TIME UNITS EXPECTED')
        IF(CARLIR(1:4) .EQ. 'S') THEN
          XT(1)=XT(1)*CSEC
          XT(2)=XT(2)*CSEC
        ELSE IF(CARLIR(1:4) .EQ. 'DAY') THEN
          XT(1)=XT(1)*CDAY
          XT(2)=XT(2)*CDAY
        ELSE IF(CARLIR(1:4) .EQ. 'YEAR') THEN
          XT(1)=XT(1)*CYEAR
          XT(2)=XT(2)*CYEAR
        ELSE
          CALL XABORT(NAMSBR//': INVALID TIME UNITS')
        ENDIF
        IF(IRDT .EQ. 2) THEN
          XT(2)=XT(4)+XT(1)
          XT(1)=XT(4)
        ENDIF
        IF(XT(2) .LE. XT(1)) CALL XABORT(NAMSBR//
     >  ': FINAL TIME IS LESS OR EQUAL TO THE INITIAL TIME.')
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': DEPLETION NORMALIZATION LEVEL EXPECTED')
        IF(CARLIR(1:4) .EQ. 'COOL') THEN
          INR=0
          RPAR(5)=0.0
        ELSE IF(CARLIR(1:4) .EQ. 'FLUX') THEN
          INR=1
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >    ': FLUX LEVEL EXPECTED')
          RPAR(5)=REALIR
          IF(RPAR(5) .EQ. 0.0) THEN
            INR=0
          ENDIF
        ELSE IF(CARLIR(1:4) .EQ. 'POWR') THEN
          INR=2
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >    ': POWER LEVEL EXPECTED')
          RPAR(5)=REALIR
          IF(RPAR(5) .EQ. 0.0) THEN
            INR=0
          ENDIF
        ELSE IF(CARLIR(1:4) .EQ. 'W/CC') THEN
          INR=3
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >    ': POWER LEVEL EXPECTED')
          RPAR(5)=REALIR
          IF(RPAR(5) .EQ. 0.0) THEN
            INR=0
          ENDIF
        ELSE IF(CARLIR(1:4) .EQ. 'KEEP') THEN
          INR=4
        ELSE
          GO TO 101
        ENDIF
      ELSE IF(CARLIR(1:4) .EQ. 'SET') THEN
        ISET=1
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': SET TIME EXPECTED')
        XT(5)=REALIR
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': TIME UNITS EXPECTED')
        IF(CARLIR(1:4) .EQ. 'S') THEN
          XT(5)=XT(5)*CSEC
        ELSE IF(CARLIR(1:4) .EQ. 'DAY') THEN
          XT(5)=XT(5)*CDAY
        ELSE IF(CARLIR(1:4) .EQ. 'YEAR') THEN
          XT(5)=XT(5)*CYEAR
        ELSE
          CALL XABORT(NAMSBR//': INVALID TIME UNITS')
        ENDIF
      ELSE
        CALL XABORT(NAMSBR//': '//CARLIR//' IS AN INVALID KEYWORD')
      ENDIF
      GO TO 100
*----
*  INPUT DATA READ COMPLETE
*----
 105  CONTINUE
      IF((ISET .EQ. 0) .AND. (IDEPL .EQ. 0)) ISET=-1
      XT(4)=0.0
      IF(ISAVE .EQ. -1) THEN
        XT(4)=-1.0
      ELSE IF(ISAVE .EQ. 0) THEN
        XT(3)=XT(1)
      ENDIF
      IF(ISET .EQ. 0) THEN
        XT(5)=XT(2)
      ENDIF
      IF(KMIXB .EQ. 1) THEN
        IF(KMIXP .EQ. 0) THEN
          DO IMIX=1,NBMIX
            MIXPWR(IMIX)=MIXBRN(IMIX)
          ENDDO
        ENDIF
      ENDIF
      RETURN
      END
