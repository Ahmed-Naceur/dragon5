*DECK CREMAC
      SUBROUTINE CREMAC(IPCPO,NISO,NGRP,NL,IMPX,HISO,DENSIT,ILEAK,TOTAL,
     1           ZNUG,SNUGF,CHI,OVERV,DIFFX,DIFFY,DIFFZ,H,SCAT,FLUX,UPS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Add the microscopic x-sections of the extracted isotopes to the
* macroscopic residual.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert
*
*Update(s):
* E. Varin (2010/01/26)
*
*Parameters: input
* IPCPO   pointer to l_compo information.
* NISO    1+number of extracted isotopes.
* NGRP    number of energy groups.
* NL      number of legendre orders (=1 for isotropic scattering).
* IMPX    print parameter (=0 for no print).
* HISO    hollerith name information for extracted isotopes.
* DENSIT  number densities.
* UPS     =.true.: no upscatering cross sections will be stored.
*
*Parameters: output
* ILEAK   diffusion coefficient flag (=1: isotropic; =2: anisotropic).
* TOTAL   total macroscopic x-sections.
* ZNUG    nu*fission macroscopic x-sections.
* SNUGF   fission macroscopic x-sections.
* CHI     fission spectrum.
* OVERV   reciprocal neutron velocities.
* DIFFX   x-directed diffusion coefficients.
* DIFFY   y-directed diffusion coefficients.
* DIFFZ   z-directed diffusion coefficients.
* H       h-factors (kappa*fission macroscopic x-sections).
* SCAT    scattering macroscopic x-sections.
* FLUX    integrated fluxes.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPCPO
      INTEGER NISO,NGRP,NL,IMPX,ILEAK,HISO(3*NISO)
      REAL    DENSIT(NISO),TOTAL(NGRP),ZNUG(NGRP),SNUGF(NGRP),CHI(NGRP),
     1        OVERV(NGRP),DIFFX(NGRP),DIFFY(NGRP),DIFFZ(NGRP),H(NGRP),
     2        SCAT(NL,NGRP,NGRP),FLUX(NGRP)
      LOGICAL UPS
*----
*  LOCAL VARIABLES
*----
      CHARACTER HMICRO*12,CM*2
      LOGICAL   LFISS
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,INDXS
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK2,ENGFIS
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WORK1
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(NGRP),NJJ(NGRP),WORK1(NGRP,3),WORK2(NGRP*NGRP),
     1 INDXS(21+NL),ENGFIS(NISO))
*----
*  RECOVER MACROSCOPIC RESIDUAL OF VECTORIAL X-SECTIONS
*----
      DO 10 IGR=1,NGRP
      TOTAL(IGR)=0.0
      DIFFX(IGR)=0.0
      DIFFY(IGR)=0.0
      DIFFZ(IGR)=0.0
      ZNUG(IGR)=0.0
      SNUGF(IGR)=0.0
      CHI(IGR)=0.0
   10 CONTINUE
      CALL LCMGET(IPCPO,'FLUX-INTG',FLUX)
      CALL LCMGET(IPCPO,'OVERV',OVERV)
      CALL LCMGET(IPCPO,'ISOTOPES-EFJ',ENGFIS)
      CALL LCMSIX(IPCPO,'MACR',1)
      CALL LCMGET(IPCPO,'XS-SAVED',INDXS)
      IF(INDXS(1).EQ.1)CALL LCMGET(IPCPO,'TOTAL',TOTAL)
      ILEAK=0
      IF(INDXS(17).EQ.1)THEN
         ILEAK=1
         CALL LCMGET(IPCPO,'STRD',DIFFX)
      ELSE IF(INDXS(18).EQ.1)THEN
         ILEAK=2
         CALL LCMGET(IPCPO,'STRD X',DIFFX)
         CALL LCMGET(IPCPO,'STRD Y',DIFFY)
         CALL LCMGET(IPCPO,'STRD Z',DIFFZ)
      ENDIF
      IF(INDXS(3).EQ.1)THEN
        CALL LCMGET(IPCPO,'NUSIGF',ZNUG)
        CALL LCMGET(IPCPO,'NFTOT',SNUGF)
        CALL LCMGET(IPCPO,'CHI',CHI)
      ENDIF
      DO 11 IGR=1,NGRP
      H(IGR)=ENGFIS(1)*SNUGF(IGR)
   11 CONTINUE
      CALL LCMSIX(IPCPO,' ',2)
*----
*  RECOVER MICROSCOPIC CONTRIBUTIONS OF VECTORIAL X-SECTIONS
*----
      LFISS=.FALSE.
      DO 40 ISO=2,NISO
      IF(DENSIT(ISO).EQ.0.)GOTO 40
      WRITE(HMICRO,'(3A4)') (HISO((ISO-1)*3+I),I=1,3)
      CALL LCMLEN(IPCPO,HMICRO,ILENG,ITYLCM)
      IF(ILENG.EQ.0)GOTO 40
      IF(IMPX.GT.1)WRITE(6,'(/29H CREMAC: PROCESSING ISOTOPE '',A12,
     1   16H'' WITH DENSITY =,1P,E13.5,2H .)') HMICRO,DENSIT(ISO)
      CALL LCMSIX(IPCPO,HMICRO,1)
      CALL LCMGET(IPCPO,'XS-SAVED',INDXS)
      IF(INDXS(1).EQ.1)THEN
        CALL LCMGET(IPCPO,'TOTAL',WORK1(1,1))
        DO 20 IGR=1,NGRP
        TOTAL(IGR)=TOTAL(IGR)+DENSIT(ISO)*WORK1(IGR,1)
   20   CONTINUE
      ENDIF
      IF(INDXS(17).EQ.1)THEN
        CALL LCMGET(IPCPO,'STRD',WORK1(1,1))
        DO 21 IGR=1,NGRP
        DIFFX(IGR)=DIFFX(IGR)+DENSIT(ISO)*WORK1(IGR,1)
   21   CONTINUE
      ELSE IF(INDXS(18).EQ.1)THEN
        CALL LCMGET(IPCPO,'STRD X',WORK1(1,1))
        CALL LCMGET(IPCPO,'STRD Y',WORK1(1,2))
        CALL LCMGET(IPCPO,'STRD Z',WORK1(1,3))
        DO 22 IGR=1,NGRP
        DIFFX(IGR)=DIFFX(IGR)+DENSIT(ISO)*WORK1(IGR,1)
        DIFFY(IGR)=DIFFY(IGR)+DENSIT(ISO)*WORK1(IGR,2)
        DIFFZ(IGR)=DIFFZ(IGR)+DENSIT(ISO)*WORK1(IGR,3)
   22   CONTINUE
      ENDIF
      IF(INDXS(3).EQ.1)THEN
        CALL LCMGET(IPCPO,'NUSIGF',WORK1(1,1))
        CALL LCMGET(IPCPO,'NFTOT',WORK1(1,2))
        CALL LCMGET(IPCPO,'CHI',WORK1(1,3))
        DO 30 IGR=1,NGRP
        LFISS=LFISS.OR.(CHI(IGR).NE.WORK1(IGR,3))
        ZNUG(IGR)=ZNUG(IGR)+DENSIT(ISO)*WORK1(IGR,1)
        SNUGF(IGR)=SNUGF(IGR)+DENSIT(ISO)*WORK1(IGR,2)
        H(IGR)=H(IGR)+DENSIT(ISO)*WORK1(IGR,2)*ENGFIS(ISO)
   30   CONTINUE
      ENDIF
      CALL LCMSIX(IPCPO,' ',2)
   40 CONTINUE
*----
*  COMPUTE AN AVERAGE FISSION SPECTRUM
*----
      IF(LFISS)THEN
        CALL LCMGET(IPCPO,'FLUX-INTG',WORK1(1,1))
        CALL LCMSIX(IPCPO,'MACR',1)
        CALL LCMGET(IPCPO,'XS-SAVED',INDXS)
        IF(INDXS(3).EQ.1)THEN
          CALL LCMGET(IPCPO,'NUSIGF',WORK1(1,2))
          CALL LCMGET(IPCPO,'CHI',WORK1(1,3))
          DO 55 JGR=1,NGRP
          DO 50 IGR=1,NGRP
          SCAT(1,IGR,JGR)=WORK1(IGR,1)*WORK1(IGR,2)*WORK1(JGR,3)
   50     CONTINUE
   55     CONTINUE
        ELSE
          DO 65 JGR=1,NGRP
          DO 60 IGR=1,NGRP
          SCAT(1,IGR,JGR)=0.
   60     CONTINUE
   65     CONTINUE
        ENDIF
        CALL LCMSIX(IPCPO,' ',2)
        DO 80 ISO=2,NISO
        IF(DENSIT(ISO).EQ.0.)GOTO 80
        WRITE(HMICRO,'(3A4)') (HISO((ISO-1)*3+I),I=1,3)
        CALL LCMLEN(IPCPO,HMICRO,ILENG,ITYLCM)
        IF(ILENG.EQ.0)GOTO 80
        CALL LCMSIX(IPCPO,HMICRO,1)
        CALL LCMGET(IPCPO,'XS-SAVED',INDXS)
        IF(INDXS(3).EQ.1)THEN
          CALL LCMGET(IPCPO,'NUSIGF',WORK1(1,2))
          CALL LCMGET(IPCPO,'CHI',WORK1(1,3))
          DO 75 JGR=1,NGRP
          DO 70 IGR=1,NGRP
          SCAT(1,IGR,JGR)=SCAT(1,IGR,JGR)+DENSIT(ISO)*
     1         WORK1(IGR,1)*WORK1(IGR,2)*WORK1(JGR,3)
   70     CONTINUE
   75     CONTINUE
        ENDIF
        CALL LCMSIX(IPCPO,' ',2)
   80   CONTINUE
        SSUM=0.
        DO 95 JGR=1,NGRP
        CHI(JGR)=0.
        DO 90 IGR=1,NGRP
          SSUM=SSUM+SCAT(1,IGR,JGR)
          CHI(JGR)=CHI(JGR)+SCAT(1,IGR,JGR)
   90   CONTINUE
   95   CONTINUE
        DO 100 JGR=1,NGRP
        CHI(JGR)=CHI(JGR)/SSUM
  100   CONTINUE
      ENDIF
*----
*  RECOVER MACROSCOPIC RESIDUAL OF SCATTERING X-SECTIONS
*----
      CALL LCMSIX(IPCPO,'MACR',1)
      CALL LCMLEN(IPCPO,'SCAT-SAVED',ILONG,ITYP)
      IF(ILONG.EQ.0)THEN
        CALL LCMGET(IPCPO,'XS-SAVED',INDXS)
      ELSE
        CALL LCMGET(IPCPO,'SCAT-SAVED',INDXS(21))
      ENDIF
      DO 130 IL=1,NL
      DO 115 JGR=1,NGRP
      DO 110 IGR=1,NGRP
      SCAT(IL,IGR,JGR)=0.
  110 CONTINUE
  115 CONTINUE
      WRITE (CM,'(I2.2)') IL-1
      IF(INDXS(20+IL).EQ.1)THEN
*     OLD COMPO DEFINITION  
        CALL LCMLEN(IPCPO,'SCAT'//CM,ILONG,ITYP)
        IF(ILONG.EQ.0)THEN
          WRITE (CM,'(I2)') IL-1
          CALL LCMGET(IPCPO,'SCAT'//CM,WORK2)
          CALL LCMGET(IPCPO,'NJJ '//CM,NJJ)
          CALL LCMGET(IPCPO,'IJJ '//CM,IJJ)
        ELSE
          CALL LCMGET(IPCPO,'SCAT'//CM,WORK2)
          CALL LCMGET(IPCPO,'NJJS'//CM,NJJ)
          CALL LCMGET(IPCPO,'IJJS'//CM,IJJ)
        ENDIF
        IGAR=0
        DO 125 JGR=1,NGRP
        DO 120 IGR=IJJ(JGR),IJJ(JGR)-NJJ(JGR)+1,-1
          IGAR=IGAR+1
          SCAT(IL,IGR,JGR)=WORK2(IGAR)
  120   CONTINUE
  125   CONTINUE
      ENDIF
  130 CONTINUE
      CALL LCMSIX(IPCPO,' ',2)
*----
*  RECOVER MICROSCOPIC CONTRIBUTIONS OF SCATTERING X-SECTIONS
*----
      DO 160 ISO=2,NISO
      IF(DENSIT(ISO).EQ.0.)GOTO 160
        WRITE(HMICRO,'(3A4)') (HISO((ISO-1)*3+I),I=1,3)
        CALL LCMLEN(IPCPO,HMICRO,ILENG,ITYLCM)
        IF(ILENG.EQ.0)GOTO 160
        CALL LCMSIX(IPCPO,HMICRO,1)
        CALL LCMLEN(IPCPO,'SCAT-SAVED',ILONG,ITYP)
*EV
        IF(ILONG.EQ.0)THEN
          CALL LCMGET(IPCPO,'XS-SAVED',INDXS)
        ELSE
          CALL LCMGET(IPCPO,'SCAT-SAVED',INDXS(21))
        ENDIF
*EV
        DO 150 IL=1,NL
        WRITE (CM,'(I2.2)') IL-1
        IF(INDXS(20+IL).EQ.1)THEN
*       OLD COMPO DEFINITION  
          CALL LCMLEN(IPCPO,'SCAT'//CM,ILONG,ITYP)
          IF(ILONG.EQ.0)THEN
            WRITE (CM,'(I2)') IL-1
            CALL LCMGET(IPCPO,'SCAT'//CM,WORK2)
            CALL LCMGET(IPCPO,'NJJ '//CM,NJJ)
            CALL LCMGET(IPCPO,'IJJ '//CM,IJJ)
          ELSE
            CALL LCMGET(IPCPO,'SCAT'//CM,WORK2)
            CALL LCMGET(IPCPO,'NJJS'//CM,NJJ)
            CALL LCMGET(IPCPO,'IJJS'//CM,IJJ)
          ENDIF
          IGAR=0
          DO 145 JGR=1,NGRP
          DO 140 IGR=IJJ(JGR),IJJ(JGR)-NJJ(JGR)+1,-1
          IGAR=IGAR+1
          SCAT(IL,IGR,JGR)=SCAT(IL,IGR,JGR)
     1                    +DENSIT(ISO)*WORK2(IGAR)
  140     CONTINUE
  145     CONTINUE
        ENDIF
  150   CONTINUE
        CALL LCMSIX(IPCPO,' ',2)
  160 CONTINUE
*----
*  COMPUTE DIFFUSION COEFFICIENTS FROM STRD X-SECTIONS
*----
      CALL LCMSIX(IPCPO,'MACR',1)
      CALL LCMGET(IPCPO,'XS-SAVED',INDXS)
      CALL LCMSIX(IPCPO,' ',2)
      IF(INDXS(17).EQ.1)THEN
        DO 170 IGR=1,NGRP
        DIFFX(IGR)=1.0/(3.0*DIFFX(IGR))
  170   CONTINUE
      ELSE IF(INDXS(18).EQ.1)THEN
        DO 180 IGR=1,NGRP
        DIFFX(IGR)=1.0/(3.0*DIFFX(IGR))
        DIFFY(IGR)=1.0/(3.0*DIFFY(IGR))
        DIFFZ(IGR)=1.0/(3.0*DIFFZ(IGR))
  180   CONTINUE
      ENDIF
*----
*  COMPUTE TOTAL CROSS SECTION FOR UPSCATERING CORRECTION
*----
      IF((UPS).AND.(NGRP.EQ.2))THEN
        DO 200 IL=1,NL
        TOTAL(2)=TOTAL(2)-SCAT(IL,2,1)
        SCAT(IL,2,1)=0.
  200   CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ENGFIS,INDXS,WORK2,WORK1,NJJ,IJJ)
      RETURN
      END
