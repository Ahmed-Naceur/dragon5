*DECK TONSN3
      SUBROUTINE TONSN3 (IPLIB,IPTRK,IFTRAK,NGRO,NBISO,NBM,NREG,NUN,
     1 CDOOR,INRS,NBNRS,IMPX,ISONAM,MIX,DEN,SN,LSHI,IPHASE,MAT,VOL,
     2 KEYFLX,LEAKSW,TITR,START,SIGT,SIGT3,NOCONV,ICPIJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one multidimensional self-shielding iteration using the
* Tone's method with Nordheim (PIC) approximation.
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal
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
* IPTRK   pointer to the tracking. (L_TRACK signature).
* IFTRAK  unit number of the sequential binary tracking file.
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
* LEAKSW  leakage flag (=.TRUE. if leakage is present on the outer
*         surface).
* TITR    title.
* START   beginning-of-iteration flag (=.TRUE. if TONSN3 is called
*         for the first time).
* SIGT3   transport correction.
* NOCONV  mixture convergence flag. (NOCONV(IBM,L)=.TRUE. if mixture IBM
*         is not converged in group L).
*
*Parameters: input/output
* SN      estimate of the dilution cross section in each energy group
*         of each isotope on input and computed dilution cross section
*         in each energy group of each isotope at output.
* SIGT    total macroscopic cross sections on ipput and total 
*         macroscopic cross sections as modified by Tone's method
*         at output.
*
*Parameters: output
* ICPIJ   number of flux solution door calls.
*
*Reference:
* A. Hebert, 'Revisiting the Stamm'ler Self-Shielding Method', Presented
* at the 25th CNS annnual conference, June 6-9, Toronto, 2004.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPTRK
      INTEGER IFTRAK,NGRO,NBISO,NBM,NREG,NUN,INRS,NBNRS,IMPX,
     1 ISONAM(3,NBISO),MIX(NBISO),LSHI(NBISO),IPHASE,MAT(NREG),
     2 KEYFLX(NREG),ICPIJ
      REAL DEN(NBISO),SN(NGRO,NBISO),VOL(NREG),SIGT(NBM,NGRO),
     1 SIGT3(NBM,NGRO)
      CHARACTER CDOOR*12,TITR*72
      LOGICAL LEAKSW,START,NOCONV(NBM,NGRO)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) KPLIB
      CHARACTER TEXT12*12,HNAMIS*12
      LOGICAL LOGDO
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IRES,ISONR,NPSYS
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR,GAS,VST,DENM
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIG0,TOTAL,SIGE
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASKI
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  SCRATCH STORAGE ALLOCATION
*   SIG0    macroscopic xs of the resonant isotopes as interpolated.
*----
      ALLOCATE(IRES(NBM),ISONR(NBISO),NPSYS(NGRO))
      ALLOCATE(SIG0(NBM,NGRO),TOTAL(NGRO,NBNRS),
     1 DENM(NBM),GAR(NGRO),GAS(NGRO),SIGE(NBNRS,NGRO),VST(NBNRS))
      ALLOCATE(MASKI(NBISO))
      ALLOCATE(IPISO(NBISO))
*----
*  FIND THE RESONANT MIXTURE NUMBERS AND THE CORRELATED ISOTOPES
*  ASSOCIATED WITH REGION INRS
*----
      CALL XDISET(IRES,NBM,0)
      CALL XDISET(ISONR,NBISO,0)
      CALL XDLSET(MASKI,NBISO,.FALSE.)
      IRS=0
      TEXT12=' '
      DO 30 IBM=1,NBM
      LOGDO=.FALSE.
      DO 10 I=1,NREG
      LOGDO=LOGDO.OR.(MAT(I).EQ.IBM)
   10 CONTINUE
      IF(.NOT.LOGDO) GO TO 30
      DO 20 ISO=1,NBISO
      LOGDO=START.OR.(DEN(ISO).NE.0.)
      IF((MIX(ISO).EQ.IBM).AND.(LSHI(ISO).EQ.INRS)) THEN
         WRITE(HNAMIS,'(3A4)') (ISONAM(I0,ISO),I0=1,3)
         IF(HNAMIS.NE.TEXT12) THEN
           IRS=IRS+1
           TEXT12=HNAMIS
           IF(LOGDO) MASKI(ISO)=.TRUE.
         ENDIF
         ISONR(ISO)=IRS
         IRES(IBM)=IRS
      ENDIF
   20 CONTINUE
   30 CONTINUE
      IF(IRS.NE.NBNRS) CALL XABORT('TONSN3: INVALID VALUE OF NBNRS.')
*----
*  SET THE LCM MICROLIB ISOTOPEWISE DIRECTORIES.
*----
      CALL LIBIPS(IPLIB,NBISO,IPISO)
*----
*  UNLOAD MICROSCOPIC X-S FROM LCM TO SCRATCH STORAGE.
*----
      DO 40 ISO=1,NBISO
      IRS=ISONR(ISO)
      IF(IRS.GT.0) THEN
        KPLIB=IPISO(ISO) ! set ISO-th isotope
        CALL LCMGET(KPLIB,'NTOT0',TOTAL(1,IRS))
      ENDIF
   40 CONTINUE
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
      DO 70 IBM=1,NBM
      IRS=IRES(IBM)
      IF(IRS.GT.0) LOGDO=LOGDO.OR.NOCONV(IBM,LLL)
   70 CONTINUE
      IF(LOGDO) THEN
         NPSYS(LLL)=LLL
*
*        COMPUTE THE LIGHT AND RESONANT COMPONENTS OF THE MACROSCOPIC
*        CROSS SECTIONS IN EACH RESONANT MIXTURE.
         DO 80 IBM=1,NBM
         DENM(IBM)=0.0
         SIG0(IBM,LLL)=0.0
   80    CONTINUE
         DO 90 ISO=1,NBISO
         IRS=ISONR(ISO)
         IF(IRS.GT.0) THEN
           IBM=MIX(ISO)
           DENM(IBM)=DEN(ISO)
           SIGT(IBM,LLL)=SIGT(IBM,LLL)-TOTAL(LLL,IRS)*DEN(ISO)
           SIG0(IBM,LLL)=TOTAL(LLL,IRS)*DEN(ISO)
         ENDIF
   90    CONTINUE
         IF(IMPX.GE.10) THEN
            WRITE (6,400) LLL,(SIG0(I,LLL),I=1,NBM)
            WRITE (6,410) LLL,(SIGT(I,LLL),I=1,NBM)
         ENDIF
      ENDIF
  110 CONTINUE
*----
*  SET UP VECTOR SIGE.
*----
      CALL LCMSIX(IPLIB,'SHIBA',1)
*
      CALL XDRSET(SIGE,NBNRS*NGRO,0.0)
      CALL LCMSIX(IPLIB,'--AVERAGE--',1)
      CALL TONDST(IPLIB,NPSYS,IPTRK,IFTRAK,CDOOR,IMPX,NBM,NBNRS,NREG,
     1 NUN,NGRO,IPHASE,MAT,VOL,KEYFLX,LEAKSW,IRES,DENM,SIG0,SIGT,SIGT3,
     2 TITR,SIGE)
      CALL LCMSIX(IPLIB,' ',2)
      DO 130 LLL=1,NGRO
      IF(NPSYS(LLL).NE.0) THEN
         ICPIJ=ICPIJ+2
         DO 120 ISO=1,NBISO
         IRS=ISONR(ISO)
         IF((LSHI(ISO).EQ.INRS).AND.(IRS.GT.0).AND.
     1   (DEN(ISO).NE.0.0)) THEN
            SN(LLL,ISO)=MAX(1.0,SIGE(IRS,LLL))
         ELSE IF((LSHI(ISO).EQ.INRS).AND.(IRS.GT.0)) THEN
            SN(LLL,ISO)=1.0E10
         ENDIF
  120    CONTINUE
      ENDIF
  130 CONTINUE
      CALL LCMSIX(IPLIB,' ',2)
*
      DO 320 LLL=1,NGRO
      LOGDO=.FALSE.
      DO 300 IBM=1,NBM
      IRS=IRES(IBM)
      IF(IRS.GT.0) LOGDO=LOGDO.OR.NOCONV(IBM,LLL)
  300 CONTINUE
      IF(LOGDO) THEN
         DO 310 ISO=1,NBISO
         IRS=ISONR(ISO)
         IF(IRS.GT.0) THEN
           IBM=MIX(ISO)
           SIGT(IBM,LLL)=SIGT(IBM,LLL)+TOTAL(LLL,IRS)*DEN(ISO)
         ENDIF
  310    CONTINUE
      ENDIF
  320 CONTINUE
*----
*  SAVE THE GROUP- AND ISOTOPE-DEPENDENT DILUTIONS
*----
      CALL LCMPUT(IPLIB,'ISOTOPESDSB',NBISO*NGRO,2,SN)
      CALL LCMPUT(IPLIB,'ISOTOPESDSN',NBISO*NGRO,2,SN)
*----
*  COMPUTE THE SELF-SHIELDED MICROSCOPIC CROSS SECTIONS AND UPDATE
*  VECTOR SIGT
*----
      IMPX2=MAX(0,IMPX-1)
      CALL LIBLIB (IPLIB,NBISO,MASKI,IMPX2)
      DO 350 ISO=1,NBISO
      IRS=ISONR(ISO)
      IF(IRS.GT.0) THEN
        IBM=MIX(ISO)
        KPLIB=IPISO(ISO) ! set ISO-th isotope
        CALL LCMGET(KPLIB,'NTOT0',GAR)
        DO 340 LLL=1,NGRO
        TOTAL(LLL,IRS)=TOTAL(LLL,IRS)-GAR(LLL)
        SIGT(IBM,LLL)=SIGT(IBM,LLL)-DEN(ISO)*TOTAL(LLL,IRS)
  340   CONTINUE
      ENDIF
  350 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IPISO)
      DEALLOCATE(MASKI)
      DEALLOCATE(VST,SIGE,GAS,GAR,DENM,TOTAL,SIG0)
      DEALLOCATE(NPSYS,ISONR,IRES)
      RETURN
*
  400 FORMAT(//51H TOTAL MACROSCOPIC CROSS SECTIONS OF THE RESONANT M,
     1 31HATERIALS IN EACH MIXTURE (GROUP,I5,2H):/(1X,1P,11E11.3))
  410 FORMAT(//51H TOTAL MACROSCOPIC CROSS SECTIONS OF THE OTHER MATE,
     1 28HRIALS IN EACH MIXTURE (GROUP,I5,2H):/(1X,1P,11E11.3))
      END
