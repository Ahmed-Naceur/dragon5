*DECK KINXSD
      SUBROUTINE KINXSD(IPMAC,NGR,NBM,NBFIS,NDG,EVL,DT,DNF,DNS,LNUD,
     1 LCHD,OVR,CHI,CHD,SGF,SGD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the 1/v and fission properties from L_MACROLIB which will be
* used for assembling source and kinetics matrix systems.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMAC  pointer to L_MACROLIB object.
* NGR    number of energy groups.
* NBM    number of material mixtures.
* NBFIS  number of fissile isotopes.
* NDG    number of delayed-neutron groups.
* EVL    steady-state eigenvalue.
* DNF    delayed neutron fractions (from module input).
* DNS    delayed neutron spectrum (from module input).
* LNUD   flag: =.true. if DNF provided from module input.
* LCHD   flag: =.true. if DNS provided from module input.
*
*Parameters: output
* OVR    reciprocal neutron velocities/DT.
* CHI    steady-state fission spectrum.
* CHD    delayed fission spectrum
* SGF    nu*fission macroscopic x-sections/keff.
* SGD    delayed nu*fission macroscopic x-sections/keff.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC
      INTEGER NGR,NBM,NBFIS,NDG
      REAL EVL,DT,DNF(NDG),DNS(NDG,NGR),OVR(NBM,NGR),CHI(NBM,NBFIS,NGR),
     1 CHD(NBM,NBFIS,NGR,NDG),SGF(NBM,NBFIS,NGR),SGD(NBM,NBFIS,NGR,NDG)
      LOGICAL LNUD,LCHD
*----
*  LOCAL VARIABLES (AUTOMATIC ALLOCATION)
*----
      LOGICAL LFIS,LFISD
      CHARACTER TEXT12*12
      TYPE(C_PTR) JPMAC,KPMAC
*----
*  PROCESS FISSION SPECTRUM TERMS.
*----
      CALL XDRSET(CHI,NBM*NBFIS*NGR,0.)
      CALL XDRSET(CHD,NBM*NBFIS*NGR*NDG,0.)
      CALL XDRSET(SGF,NBM*NBFIS*NGR,0.)
      CALL XDRSET(SGD,NBM*NBFIS*NGR*NDG,0.)
      JPMAC=LCMGID(IPMAC,'GROUP')
      KPMAC=LCMGIL(JPMAC,1)
      CALL LCMLEN(KPMAC,'CHI',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
        IF(LENGT.NE.NBM*NBFIS) CALL XABORT('@KINXSD: INVALID LENGTH FO'
     1  //'R CHI INFORMATION.')
        DO 10 IGR=1,NGR
        KPMAC=LCMGIL(JPMAC,IGR)
        CALL LCMGET(KPMAC,'CHI',CHI(1,1,IGR))
   10   CONTINUE
      ELSE
        DO 22 IBM=1,NBM
        DO 21 IFIS=1,NBFIS
        CHI(IBM,IFIS,1)=1.0
        DO 20 IGR=2,NGR
        CHI(IBM,IFIS,IGR)=0.0
   20   CONTINUE
   21   CONTINUE
   22   CONTINUE
      ENDIF
      IF(LCHD) THEN
        DO 33 IDEL=1,NDG
        DO 32 IGR=1,NGR
        DO 31 IFIS=1,NBFIS
        DO 30 IBM=1,NBM
        CHD(IBM,IFIS,IGR,IDEL)=DNS(IDEL,IGR)
   30   CONTINUE
   31   CONTINUE
   32   CONTINUE
   33   CONTINUE
      ELSE
        KPMAC=LCMGIL(JPMAC,1)
        CALL LCMLEN(KPMAC,'CHI01',LENGT,ITYLCM)
        IF(LENGT.GT.0) THEN
          IF(LENGT.NE.NBM*NBFIS) CALL XABORT('@KINXSD: INVALID LENGTH '
     1    //'FOR DELAYED CHI INFORMATION.')
          DO 42 IDEL=1,NDG
          WRITE(TEXT12,'(3HCHI,I2.2)') IDEL
          DO 40 IGR=1,NGR
          KPMAC=LCMGIL(JPMAC,IGR)
          CALL LCMGET(KPMAC,TEXT12,CHD(1,1,IGR,IDEL))
   40     CONTINUE
   42     CONTINUE
        ELSE
          CALL XDRSET(CHD,NBM*NBFIS*NGR*NDG,0.)
        ENDIF
      ENDIF
      LFIS=.FALSE.
      LFISD=.FALSE.
      DO 52 IGR=1,NGR
      DO 51 IFIS=1,NBFIS
      DO 50 IBM=1,NBM
      LFIS=LFIS.OR.(CHI(IBM,IFIS,IGR).NE.0.0)
      LFISD=LFISD.OR.(CHD(IBM,IFIS,IGR,1).NE.0.0)
   50 CONTINUE
   51 CONTINUE
   52 CONTINUE
*
      DO 85 IGR=1,NGR
      KPMAC=LCMGIL(JPMAC,IGR)
*----
*  PROCESS FISSION NUSIGF TERMS.
*----
      IF(LFIS) THEN
        CALL LCMLEN(KPMAC,'NUSIGF',LENGT,ITYLCM)
        IF(LENGT.NE.NBM*NBFIS) CALL XABORT('@KINXSD: INVALID LENGTH FO'
     1  //'R NUSIGF INFORMATION.')
        IF(LENGT.GT.0) CALL LCMGET(KPMAC,'NUSIGF',SGF(1,1,IGR))
      ENDIF
      IF(LNUD) THEN
        DO 62 IDEL=1,NDG
        DO 61 IFIS=1,NBFIS
        DO 60 IBM=1,NBM
        SGD(IBM,IFIS,IGR,IDEL)=SGF(IBM,IFIS,IGR)*DNF(IDEL)
   60   CONTINUE
   61   CONTINUE
   62   CONTINUE
      ELSE IF(LFISD) THEN
        DO 70 IDEL=1,NDG
        WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL
        CALL LCMLEN(KPMAC,TEXT12,LENGT,ITYLCM)
        IF(LENGT.NE.NBM*NBFIS) CALL XABORT('@KINXSD: INVALID LENGTH FO'
     1  //'R DELAYED NUSIGF INFORMATION.')
        IF(LENGT.GT.0) CALL LCMGET(KPMAC,TEXT12,SGD(1,1,IGR,IDEL))
   70   CONTINUE
      ENDIF
*----
*  PROCESS 1/V TERMS.
*----
      CALL LCMLEN(KPMAC,'OVERV',LENGT,ITYLCM)
      IF(LENGT.EQ.NBM)THEN
        CALL LCMGET(KPMAC,'OVERV',OVR(1,IGR))
      ELSEIF(LENGT.EQ.0)THEN
        CALL XABORT('@KINXSD: MISSING OVERV DATA.')
      ELSE
        CALL XABORT('@KINXSD: INVALID OVERV DATA.')
      ENDIF
      DO 80 IBM=1,NBM
      OVR(IBM,IGR)=OVR(IBM,IGR)/DT
   80 CONTINUE
   85 CONTINUE
*
      DO 93 IGR=1,NGR
      DO 92 IFIS=1,NBFIS
      DO 91 IBM=1,NBM
      SGF(IBM,IFIS,IGR)=SGF(IBM,IFIS,IGR)/EVL
      DO 90 IDEL=1,NDG
      SGD(IBM,IFIS,IGR,IDEL)=SGD(IBM,IFIS,IGR,IDEL)/EVL
   90 CONTINUE
   91 CONTINUE
   92 CONTINUE
   93 CONTINUE
      RETURN
      END
