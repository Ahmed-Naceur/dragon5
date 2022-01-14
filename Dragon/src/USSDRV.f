*DECK USSDRV
      SUBROUTINE USSDRV(IPLI0,IPTRK,IPLIB,IFTRAK,INDREC,CDOOR,IMPX,
     1 IGRMIN,IGRMAX,NGRP,NBMIX,NREG,NUN,NBISO,NL,NED,NDEL,LEAKSW,
     2 ITRANC,IPHASE,TITR,KSPH,NRES,NPASS,ICALC,ICORR,ISUBG,MAXST,
     3 LFLAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for a resonance self-shielding calculation.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLI0   pointer to the internal microscopic cross section library
*         builded by the self-shielding module (L_LIBRARY signature).
* IPTRK   pointer to the tracking (L_TRACK signature).
* IPLIB   pointer to the internal microscopic cross section library
*         with subgroups (L_LIBRARY signature).
* IFTRAK  unit number of the sequential binary tracking file.
* INDREC  access flag for the internal microscopic cross section library
*         builded by the self-shielding module (=1 IPLI0 access in
*         creation mode; =2 in modification mode).
* CDOOR   name of the geometry/solution operator.
* IMPX    print flag.
* IGRMIN  first group where the self-shielding is applied.
* IGRMAX  most thermal group where the self-shielding is applied.
* NGRP    number of energy groups.
* NBMIX   number of mixtures in the internal library.
* NREG    number of regions.
* NUN     number of unknowns per energy group.
* NBISO   number of isotopes specifications in the internal library.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NED     number of extra vector edits.
* NDEL    number of delayed neutron precursor groups.
* LEAKSW  leakage flag (LEAKSW=.TRUE. if neutron leakage through
*         external boundary is present).
* ITRANC  type of transport correction.
* IPHASE  type of flux solution (=1 use a native flux solution door;
*         =2 use collision probabilities).
* TITR    title.
* KSPH    SPH equivalence flag (=0 no SPH correction; =1 SPH correction
*         in the fuel).
* NRES    number of self-shielding zones, as given by LIB:.
* NPASS   number of outer iterations.
* ICALC   simplified self-shielding flag (=1 IPLI0 is containing ICALC
*         data. =0 no ICALC data).
* ICORR   mutual resonance shielding flag (=1 to suppress the model
*         in cases it is required in LIB operator).
* ISUBG   type of self-shielding model (=1 use physical probability
*         tables; =3 use original Ribon method; =4 use Ribon extended
*         method).
* MAXST   maximum number of fixed point iterations for the ST scattering
*         source.
* LFLAT   force the initial subgroup flux to be flat if IPLI0 is open
*         in modification mode.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLI0,IPTRK,IPLIB
      INTEGER IFTRAK,INDREC,IMPX,IGRMIN,IGRMAX,NGRP,NBMIX,NREG,NUN,
     1 NBISO,NL,NED,NDEL,ITRANC,IPHASE,KSPH,NRES,NPASS,ICALC,ICORR,
     2 ISUBG,MAXST
      CHARACTER CDOOR*12,TITR*72
      LOGICAL LEAKSW,LFLAT
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,MAXRSS=300,MAXESP=4)
      TYPE(C_PTR) JPLI0,KPLI0,JPLIB,KPLIB
      CHARACTER HSMG*131,HCAL*12,TEXT4*4,NAM1*4,FNAM1*4,NAM2*12,
     1 FNAM2*12,CBDPNM*12,TEXT8*8
      INTEGER IPAR(NSTATE),IRSS(MAXRSS),IESP(MAXESP+1)
      REAL TMPDAY(3),EESP(MAXESP+1)
      LOGICAL LTEST
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAT,KEYFLX,MIX,IEVOL,ITYPE,
     1 LSHI,IAPT,IHSUF,IREX,ILLIB,JCEDM,LSHI2
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISONAM,ISONRF,IHLIB
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL,TN,DEN,ENER,GS
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASK,MASKL
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MAT(NREG),KEYFLX(NREG),ISONAM(3,NBISO),ISONRF(3,NBISO),
     3 MIX(NBISO),IEVOL(NBISO),ITYPE(NBISO),LSHI(NBISO),IAPT(NBISO),
     4 IHSUF(NBISO),IREX(NBMIX),IHLIB(2,NBISO),ILLIB(NBISO))
      ALLOCATE(VOL(NREG),TN(NBISO),DEN(NBISO))
*----
*  RECOVER USEFUL INFORMATION FROM TRACKING OBJECT.
*----
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMGET(IPTRK,'KEYFLX',KEYFLX)
*----
*  RECOVER USEFUL INFORMATION FROM LIBRARY OBJECTS.
*----
      CALL LCMGET(IPLIB,'ISOTOPESUSED',ISONAM)
      CALL LCMGET(IPLIB,'ISOTOPERNAME',ISONRF)
      CALL LCMGET(IPLIB,'ISOTOPESMIX',MIX)
      CALL LCMGET(IPLIB,'ISOTOPESTODO',IEVOL)
      CALL LCMGET(IPLIB,'ISOTOPESTYPE',ITYPE)
      CALL LCMGET(IPLIB,'ISOTOPESTEMP',TN)
*
      CALL LCMPUT(IPLI0,'ISOTOPESMIX',NBISO,1,MIX)
      CALL LCMPUT(IPLI0,'ISOTOPESTODO',NBISO,1,IEVOL)
      CALL LCMPUT(IPLI0,'ISOTOPESTYPE',NBISO,1,ITYPE)
      CALL LCMPUT(IPLI0,'ISOTOPESTEMP',NBISO,2,TN)
      IF(INDREC.EQ.1) THEN
         CALL LCMGET(IPLIB,'ISOTOPESDENS',DEN)
         CALL LCMPUT(IPLI0,'ISOTOPESDENS',NBISO,2,DEN)
      ELSE IF(INDREC.EQ.2) THEN
         CALL LCMGET(IPLI0,'ISOTOPESDENS',DEN)
      ENDIF
      CALL LCMGET(IPLIB,'ISOTOPESSHI',LSHI)
      CALL LCMLEN(IPLIB,'ISOTOPESDSN',NELSN,ITYLCM)
      IF(NELSN.GT.0) THEN
        NGIS=NGRP*NBISO
        ALLOCATE(GS(NGIS))
        CALL LCMGET(IPLIB,'ISOTOPESDSN',GS)
        CALL LCMPUT(IPLI0,'ISOTOPESDSN',NGIS,2,GS)
        CALL LCMGET(IPLIB,'ISOTOPESDSB',GS)
        CALL LCMPUT(IPLI0,'ISOTOPESDSB',NGIS,2,GS)
        DEALLOCATE(GS)
      ENDIF
      ALLOCATE(ENER(NGRP+1))
      CALL LCMGET(IPLIB,'ENERGY',ENER)
      CALL LCMPUT(IPLI0,'ENERGY',NGRP+1,2,ENER)
      CALL LCMGET(IPLIB,'DELTAU',ENER)
      CALL LCMPUT(IPLI0,'DELTAU',NGRP,2,ENER)
      DEALLOCATE(ENER)
      CALL LCMLEN(IPLIB,'CHI-LIMITS',NBESP,ITYLCM)
      IF(NBESP.GT.0) THEN
         NBESP=NBESP-1
         IF(NBESP.GT.MAXESP) CALL XABORT('USSDRV: MAXESP OVERFLOW.')
         CALL LCMGET(IPLIB,'CHI-LIMITS',IESP)
         CALL LCMPUT(IPLI0,'CHI-LIMITS',NBESP+1,1,IESP)
         CALL LCMGET(IPLIB,'CHI-ENERGY',EESP)
         CALL LCMPUT(IPLI0,'CHI-ENERGY',NBESP+1,2,EESP)
      ENDIF
      DO 10 ISO=1,NBISO
      DO 5 I=1,NREG
      IF(MAT(I).EQ.MIX(ISO)) GO TO 10
    5 CONTINUE
      LSHI(ISO)=0
   10 CONTINUE
*
      CALL LCMLEN(IPLIB,'ILIBRARYTYPE',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
         CALL LCMGET(IPLIB,'ILIBRARYTYPE',IHLIB(1,1))
         CALL LCMGET(IPLIB,'ILIBRARYINDX',ILLIB)
         DO 15 ISO=1,NBISO
         IF(LSHI(ISO).NE.0) THEN
            TEXT8=' '
            READ(TEXT8,'(2A4)') IHLIB(1,ISO),IHLIB(2,ISO)
            ILLIB(ISO)=0
         ENDIF
   15    CONTINUE
         CALL LCMPUT(IPLI0,'ILIBRARYTYPE',2*NBISO,1,IHLIB(1,1))
         CALL LCMPUT(IPLI0,'ILIBRARYINDX',NBISO,1,ILLIB)
      ENDIF
*
      JPLIB=LCMGID(IPLIB,'ISOTOPESLIST')
      JPLI0=LCMLID(IPLI0,'ISOTOPESLIST',NBISO)
      IF(INDREC.EQ.1) THEN
*        COPY THE NON RESONANT ISOTOPES.
         CALL KDRCPU(TK1)
         DO 20 ISO=1,NBISO
         IF((LSHI(ISO).EQ.0).OR.(DEN(ISO).EQ.0.0)) THEN
            CALL LCMLEL(JPLIB,ISO,ILEN,ITYLCM)
            IF(ILEN.EQ.0) THEN
              DO JSO=1,ISO-1
                CALL LCMLEL(JPLIB,JSO,ILEN,ITYLCM)
                IF(ILEN.EQ.0) CYCLE
                IF((ISONAM(1,ISO).EQ.ISONAM(1,JSO)).AND.(ISONAM(2,ISO)
     1          .EQ.ISONAM(2,JSO)).AND.(ISONAM(3,ISO).EQ.ISONAM(3,JSO)))
     2          THEN
                  IF(LSHI(JSO).GT.0) THEN
                    KPLIB=LCMGIL(JPLIB,JSO) ! set JSO-th isotope
                    GO TO 16
                  ELSE
                    GO TO 20
                  ENDIF
                ENDIF
              ENDDO
            ELSE
              KPLIB=LCMGIL(JPLIB,ISO) ! set ISO-th isotope
              GO TO 16
            ENDIF
            GO TO 20
   16       CALL LCMLEL(JPLI0,ISO,ILEN,ITYLCM)
            IF(ILEN.NE.0) GO TO 20
            KPLI0=LCMDIL(JPLI0,ISO) ! set ISO-th isotope
            CALL LCMEQU(KPLIB,KPLI0)
         ENDIF
   20    CONTINUE
         CALL KDRCPU(TK2)
         IF(IMPX.GT.1) WRITE(6,'(/33H USSDRV: CPU TIME SPENT TO COPY T,
     1   26HHE NON-RESONANT ISOTOPES =,F8.1,8H SECOND.)') TK2-TK1
*
*        WRITE THE OUTPUT INTERNAL LIBRARY PARAMETERS.
         CALL LCMGET(IPLIB,'STATE-VECTOR',IPAR)
         IPAR(8)=0
         IPAR(17)=0
         CALL LCMPUT(IPLI0,'STATE-VECTOR',NSTATE,1,IPAR)
         IF(NED.GT.0) THEN
            ALLOCATE(JCEDM(2*NED))
            CALL LCMGET(IPLIB,'ADDXSNAME-P0',JCEDM)
            CALL LCMPUT(IPLI0,'ADDXSNAME-P0',2*NED,3,JCEDM)
            DEALLOCATE(JCEDM)
         ENDIF
         CALL LCMLEN(IPLIB,'DEPL-CHAIN',ILENG,ITYLCM)
         IF(ILENG.NE.0) THEN
            CALL LCMSIX(IPLIB,'DEPL-CHAIN',1)
            CALL LCMSIX(IPLI0,'DEPL-CHAIN',1)
            CALL LCMEQU(IPLIB,IPLI0)
            CALL LCMSIX(IPLI0,' ',2)
            CALL LCMSIX(IPLIB,' ',2)
         ENDIF
      ENDIF
*----
*  RECOMPUTE THE VECTOR LSHI.
*----
      ALLOCATE(LSHI2(NBISO))
      NRES2=0
      NRES3=0
      DO 30 ISO=1,NBISO
      IF(LSHI(ISO).NE.0) NRES3=NRES3+1
      LSHI2(ISO)=0
   30 CONTINUE
      DO 80 INRS=1,NRES
   40 DENMAX=0.0
      KSOT=0
      DO 60 ISO=1,NBISO
      IF(LSHI2(ISO).EQ.0) THEN
         VOLISO=0.0
         DO 50 I=1,NREG
         IF(MAT(I).EQ.MIX(ISO)) VOLISO=VOLISO+VOL(I)
   50    CONTINUE
         IF((ABS(LSHI(ISO)).EQ.INRS).AND.(DEN(ISO)*VOLISO.GT.DENMAX))
     1   THEN
            KSOT=ISO
            DENMAX=DEN(ISO)*VOLISO
         ENDIF
      ENDIF
   60 CONTINUE
      IF(KSOT.GT.0) THEN
         NRES2=NRES2+1
         DO 70 ISO=1,NBISO
         LTEST=((ISONRF(1,ISO).EQ.ISONRF(1,KSOT)).AND.
     1      (ISONRF(2,ISO).EQ.ISONRF(2,KSOT)).AND.
     2      (ISONRF(3,ISO).EQ.ISONRF(3,KSOT)).AND.
     3      (ABS(LSHI(ISO)).EQ.INRS))
         LTEST=LTEST.OR.((ISONAM(1,ISO).EQ.ISONAM(1,KSOT)).AND.
     1      (ISONAM(2,ISO).EQ.ISONAM(2,KSOT)).AND.
     2      (ABS(LSHI(ISO)).EQ.INRS))
         IF(LTEST) LSHI2(ISO)=NRES2
         IF(LTEST.AND.(LSHI(ISO).EQ.-INRS)) THEN
            DO 65 JSO=1,NBISO
            IF(LSHI(JSO).EQ.LSHI(ISO)) LSHI2(JSO)=NRES2
   65       CONTINUE
         ENDIF
   70    CONTINUE
         GO TO 40
      ENDIF
   80 CONTINUE
*----
* FIND THE ISOTOPE-NAME SUFFIX VALUES.
*----
      TEXT4=' '
      READ(TEXT4,'(A4)') IHBLK
      DO 90 ISO=1,NBISO
      IF((LSHI2(ISO).NE.0).AND.(DEN(ISO).NE.0.0)) THEN
         WRITE(TEXT4,'(I4.4)') MIX(ISO)
         READ(TEXT4,'(A4)') IHSUF(ISO)
      ELSE
         IHSUF(ISO)=IHBLK
      ENDIF
   90 CONTINUE
      IF(ICALC.EQ.1) THEN
         CALL LCMSIX(IPLI0,'SHIBA_SG',1)
         CALL LCMSIX(IPLI0,'-DATA-CALC-',1)
         NAM1=' '
         CALL LCMNXT(IPLI0,NAM1)
         FNAM1=NAM1
  100    CALL LCMSIX(IPLI0,NAM1,1)
            NAM2=' '
            CALL LCMNXT(IPLI0,NAM2)
            FNAM2=NAM2
  110       CALL LCMLEN(IPLI0,NAM2,NRSS,ITYLCM)
            CALL LCMGET(IPLI0,NAM2,IRSS)
            READ(NAM2,'(2A4)') IN1,IN2
            DO 130 ISO=1,NBISO
            IF((ISONAM(1,ISO).EQ.IN1).AND.(ISONAM(2,ISO).EQ.IN2).AND.
     1      (LSHI2(ISO).NE.0)) THEN
               IF((NRSS.EQ.1).AND.(IRSS(1).EQ.-999)) THEN
                  READ(NAM1,'(A4)') IHSUF(ISO)
               ELSE
                  DO 120 I=1,NRSS
                  IF(IRSS(I).EQ.MIX(ISO)) READ(NAM1,'(A4)') IHSUF(ISO)
  120             CONTINUE
               ENDIF
            ENDIF
  130       CONTINUE
            CALL LCMNXT(IPLI0,NAM2)
            IF(NAM2.EQ.FNAM2) GO TO 140
            GO TO 110
  140    CALL LCMSIX(IPLI0,' ',2)
         CALL LCMNXT(IPLI0,NAM1)
         IF(NAM1.EQ.FNAM1) THEN
            CALL LCMSIX(IPLI0,' ',2)
            CALL LCMSIX(IPLI0,' ',2)
            GO TO 150
         ENDIF
         GO TO 100
      ENDIF
*
  150 NPASS2=NPASS
      IF(NRES3.EQ.1) NPASS2=1
      DO 265 IPASS=1,NPASS2
      IF((IMPX.GT.0).AND.(NPASS2.GT.1)) WRITE (6,'(/15H USSDRV: SELF S,
     1 25HHIELDING ITERATION NUMBER,I4,1H.)') IPASS
      DO 260 INRS=1,NRES2
*----
*  COMPUTE THE NUMBER OF RESONANT ISOTOPES IN REGION INRS AND THE
*  RESONANT ISOTOPE INDEX ASSOCIATED TO EACH ISOTOPE SPECIFICATION.
*----
      NIRES=0
      DO 200 ISO=1,NBISO
      IAPT(ISO)=0
      IF((LSHI2(ISO).EQ.INRS).AND.(DEN(ISO).NE.0.0)) THEN
         DO 170 I=1,NREG
         IF(MAT(I).EQ.MIX(ISO)) GO TO 180
  170    CONTINUE
         GO TO 200
  180    DO 190 JSO=1,ISO-1
         IF((ISONAM(1,ISO).EQ.ISONAM(1,JSO)).AND.
     1      (ISONAM(2,ISO).EQ.ISONAM(2,JSO)).AND.
     2      (ISONAM(3,ISO).EQ.ISONAM(3,JSO)).AND.
     3      (LSHI2(JSO).EQ.INRS).AND.
     4      (DEN(JSO).NE.0.0).AND.(IAPT(JSO).NE.0)) THEN
              IAPT(ISO)=IAPT(JSO)
              GO TO 200
         ENDIF
  190    CONTINUE
         IIII=ISO
         NIRES=NIRES+1
         IAPT(ISO)=NIRES
      ENDIF
  200 CONTINUE
      WRITE(HCAL,'(1HC,I5.5,1H/,I5.5)') IIII,NBISO
      IF(NIRES.EQ.0) THEN
         WRITE(HSMG,'(45HUSSDRV: NO RESONANT ISOTOPES IN RESONANT REGI,
     1   9HON NUMBER,I4,7H (HCAL=,A12,2H).)') INRS,HCAL
         CALL XABORT(HSMG)
      ENDIF
      IF(IMPX.GT.0) WRITE (6,'(/35H USSDRV: PERFORMING SELF-SHIELDING ,
     1 18HCALCULATION NAMED ,A12,1H.)') HCAL
*----
*  FIND THE NUMBER OF FUEL REGIONS AND THE FUEL REGION INDICES ASSIGNED
*  TO EACH RESONANT MIXTURE.
*----
      NBNRS=0
      DO 210 IBM=1,NBMIX
      IREX(IBM)=0
  210 CONTINUE
      DO 230 ISO=1,NBISO
      IBM=MIX(ISO)
      IF((IAPT(ISO).GT.0).AND.(IREX(IBM).EQ.0)) THEN
         DO 220 JSO=1,ISO-1
         IF((IHSUF(JSO).EQ.IHSUF(ISO)).AND.(IAPT(JSO).EQ.IAPT(ISO)))
     1   THEN
            IREX(IBM)=IREX(MIX(JSO))
            GO TO 230
         ENDIF
  220    CONTINUE
         IF(IMPX.GT.0) WRITE(6,'(9X,3H-->,3A4)') (ISONAM(J,ISO),J=1,3)
         NBNRS=NBNRS+1
         IREX(IBM)=NBNRS
      ELSE IF(IAPT(ISO).GT.0) THEN
         IF(IMPX.GT.0) WRITE(6,'(9X,3H-->,3A4)') (ISONAM(J,ISO),J=1,3)
      ENDIF
  230 CONTINUE
      IF(NBNRS.EQ.0) THEN
         WRITE (HSMG,'(33HUSSDRV: INVALID RESONANT REGION =,I10)') INRS
         CALL XABORT(HSMG)
      ENDIF
      IF(IMPX.GE.1) WRITE(6,410) NIRES,NBNRS,INRS
*----
*  DETERMINE WHICH MODERATOR ISOTOPES ARE MIXED WITH RESONANT ONES.
*----
      DO 250 ISO=1,NBISO
      IF((IAPT(ISO).EQ.0).AND.(IREX(MIX(ISO)).GT.0)) IAPT(ISO)=NIRES+1
  250 CONTINUE
*----
*  ERASE OLD GROUP_INFO, ASSEMB_PHYS AND ASSEMB_RIBON DIRECTORIES.
*----
      IF(LFLAT.AND.(IPASS.EQ.1).AND.(INDREC.EQ.2)) THEN
         CALL LCMSIX(IPLI0,'SHIBA_SG',1)
         CALL LCMSIX(IPLI0,HCAL,1)
         DO IRES=1,NIRES
            WRITE(CBDPNM,'(3HCOR,I4.4,1H/,I4.4)') IRES,NIRES
            CALL LCMSIX(IPLI0,CBDPNM,1)
            CALL LCMLEN(IPLI0,'GROUP_INFO',ILONG,ITYLCM)
            IF(ILONG.GT.0) CALL LCMDEL(IPLI0,'GROUP_INFO')
            CALL LCMLEN(IPLI0,'ASSEMB_PHYS',ILONG,ITYLCM)
            IF(ILONG.GT.0) CALL LCMDEL(IPLI0,'ASSEMB_PHYS')
            CALL LCMLEN(IPLI0,'ASSEMB_RIBON',ILONG,ITYLCM)
            IF(ILONG.GT.0) CALL LCMDEL(IPLI0,'ASSEMB_RIBON')
            CALL LCMSIX(IPLI0,' ',2)
         ENDDO
         CALL LCMSIX(IPLI0,' ',2)
         CALL LCMSIX(IPLI0,' ',2)
      ENDIF
*----
*  PERFORM A SELF-SHIELDING CALCULATION NAMED HCAL.
*----
      CALL USSONE(IPLI0,IPTRK,IPLIB,IFTRAK,CDOOR,IMPX,IGRMIN,IGRMAX,
     1 NIRES,NBNRS,IREX,NGRP,NBMIX,NREG,NUN,NBISO,NL,NED,NDEL,ISONAM,
     2 IHSUF,HCAL,DEN,MIX,IAPT,MAT,VOL,KEYFLX,LEAKSW,ITRANC,IPHASE,
     3 TITR,KSPH,ICORR,ISUBG,MAXST)
  260 CONTINUE
  265 CONTINUE
      DEALLOCATE(LSHI2)
      IF(IMPX.GE.3) CALL LCMLIB(IPLI0)
*----
*  BUILD THE MACROLIB IN THE OUTPUT INTERNAL LIBRARY.
*----
      ALLOCATE(MASK(NBMIX))
      DO 280 IBM=1,NBMIX
      MASK(IBM)=.TRUE.
      DO 270 I=1,NREG
      IF(MAT(I).EQ.IBM) GO TO 280
  270 CONTINUE
      MASK(IBM)=.FALSE.
  280 CONTINUE
      ALLOCATE(MASKL(NGRP))
      DO 290 I=1,NGRP
      MASKL(I)=.TRUE.
  290 CONTINUE
*
      ITSTMP=0
      TMPDAY(1)=0.0
      TMPDAY(2)=0.0
      TMPDAY(3)=0.0
      CALL KDRCPU(TK1)
      CALL LCMGET(IPLI0,'ISOTOPESUSED',ISONAM)
      CALL LIBMIX(IPLI0,NBMIX,NGRP,NBISO,ISONAM,MIX,DEN,MASK,MASKL,
     1 ITSTMP,TMPDAY)
      CALL KDRCPU(TK2)
      IF(IMPX.GT.1) WRITE(6,'(/37H USSDRV: CPU TIME SPENT TO BUILD THE ,
     1 19HEMBEDDED MACROLIB =,F8.1,8H SECOND.)') TK2-TK1
      DEALLOCATE(MASKL,MASK)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DEN,TN,VOL)
      DEALLOCATE(ILLIB,IHLIB,IREX,IHSUF,IAPT,LSHI,ITYPE,IEVOL,MIX,
     1 ISONRF,ISONAM,KEYFLX,MAT)
      RETURN
*
  410 FORMAT(/48H USSDRV: NUMBER OF CORRELATED RESONANT ISOTOPES=,I4/9X,
     1 35HNUMBER OF CORRELATED FUEL MIXTURES=,I4,19H IN RESONANT REGION,
     2 I3)
      END
