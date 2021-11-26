*DECK FLPOWB
      SUBROUTINE FLPOWB(IPPOW,IPMAP,IPMTX,NMIX,NMAT,NGRP,NCH,NB,NEL,MAT,
     1 VOL,HFAC,FLUX,POWB,POWC,IMPX,PTOT,FSTH,LFSTH,FMIX,FLUB,IGEO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the channel and bundle powers over the fuel lattice.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki, M. Guyot, V. Descotes
*
*Parameters: input
* IPPOW  pointer to power information.
* IPMAP  pointer to fuel-map information.
* IPMTX  pointer to matex information.
* NMIX   maximum number of material mixtures.
* NMAT   total number of mixtures (includes virtual regions).
* NGRP   number of energy groups.
* NCH    number of reactor channels.
* NB     number of fuel bundles per channel.
* NEL    total number of finite elements.
* MAT    index-number of mixture assigned to each volume.
* VOL    element-ordered mesh-splitted volumes.
* HFAC   h-factors over the reactor core.
* FLUX   normalized average fluxes associated with each volume.
* IMPX   printing index (=0 for no print).
* PTOT   total power in MW
* FSTH   thermal to fission ratio power
* LFSTH  boolean =.true. if FSTH is specified
* FMIX   fuel bundle indices.
* FLUB   normalized average fluxes associated with each bundle
* IGEO   type of the geometry (=7 or =9)
*
*Parameters: output
* POWB   bundle powers in kW.
* POWC   channel powers in kW.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPPOW,IPMAP,IPMTX
      INTEGER NMIX,NMAT,NCH,NB,NGRP,NEL,MAT(NEL),IMPX,NX,NY,NZ,IGEO,
     1        FMIX(NCH*NB)
      REAL HFAC(NMIX,NGRP),VOL(NEL),FLUX(NEL,NGRP),BFACT1,CFACT1,
     1     POWB(NCH,NB),POWC(NCH),FSTH,FLUB(NCH,NB,NGRP)
      LOGICAL LFSTH
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40,IOUT=6)
      INTEGER IGST(NSTATE),FMAT(NMAT),IB,ICH,IEL,ICMX,IBMX,IPCH
      DOUBLE PRECISION POWER,BAVG,CAVG,BFACT
      REAL PBMX,VTOT,VOLB(NCH,NB)
      CHARACTER TEXT*12
      TYPE(C_PTR) JPMAP
*----
*  BUNDLE POWERS
*----
      PBMX=0.
      BAVG=0.0D0
      POWER=0.0D0
      IBMX=0
      ICMX=0
      CALL XDISET(FMAT,NMAT,0)
      CALL LCMGET(IPMTX,'MAT',FMAT)
      CALL XDRSET(POWB,NCH*NB,0.)
      IF(IMPX.GT.0)WRITE(IOUT,1004)
*
      NTOT=0
      DO 35 IB=1,NB
      DO 30 ICH=1,NCH
      POWB(ICH,IB)=0.0
      VOLB(ICH,IB)=0.0
      NUM=(IB-1)*NCH+ICH
      IF(FMIX(NUM).EQ.0) GO TO 30
      NTOT=NTOT+1
      DO 20 IEL=1,NEL
      IF((FMAT(IEL).EQ.-NTOT).AND.(MAT(IEL).GT.0)) THEN
         DO 10 JGR=1,NGRP
         POWB(ICH,IB)=POWB(ICH,IB)+
     1     FLUX(IEL,JGR)*HFAC(MAT(IEL),JGR)*VOL(IEL)
   10    CONTINUE
         VOLB(ICH,IB)=VOLB(ICH,IB)+VOL(IEL)
      ENDIF
   20 CONTINUE
      POWER=POWER+DBLE(POWB(ICH,IB))
   30 CONTINUE
   35 CONTINUE
      POWER=POWER/(10**6)
      VTOT=0.0
      DO 45 IB=1,NB
      DO 40 ICH=1,NCH
      POWB(ICH,IB)=POWB(ICH,IB)/1000.
      IF(POWB(ICH,IB).GT.PBMX)THEN
        PBMX=POWB(ICH,IB)
        ICMX=ICH
        IBMX=IB
      ENDIF
      BAVG=BAVG+DBLE(POWB(ICH,IB)*VOLB(ICH,IB))
      VTOT=VTOT+VOLB(ICH,IB)
   40 CONTINUE
   45 CONTINUE
      BAVG=BAVG/VTOT
      BFACT=BAVG/PBMX

*     CHECK TOTAL POWER
      IF(IMPX.EQ.99)WRITE(IOUT,1000)POWER
      IF((IMPX.EQ.0).OR.(IMPX.GT.1))GOTO 50
      WRITE(TEXT,'(A9,I3.3)')'CHANNEL #',ICMX
      IF(PBMX.LT.1000.)THEN
        WRITE(IOUT,1001)PBMX,TEXT,IBMX
      ELSE
        WRITE(IOUT,1011)PBMX,TEXT,IBMX
      ENDIF
      IF(BAVG.LT.1000.)THEN
        WRITE(IOUT,1007)BAVG
      ELSE
        WRITE(IOUT,1012)BAVG
      ENDIF
      FACT=1./REAL(BFACT)
      WRITE(IOUT,1009)BFACT,FACT
*----
*  CHANNEL POWERS
*----
   50 PCMX=0.
      CAVG=0.0D0
      POWER=0.0D0
      CALL XDRSET(POWC,NCH,0.)
      DO 70 ICH=1,NCH
      POWC(ICH)=0.
      VOLCH=0.0
      DO 60 IB=1,NB
      POWC(ICH)=POWC(ICH)+POWB(ICH,IB)
      VOLCH=VOLCH+VOLB(ICH,IB)
   60 CONTINUE
      POWER=POWER+DBLE(POWC(ICH))
      IF(POWC(ICH).GT.PCMX)THEN
        PCMX=POWC(ICH)
        IPCH=ICH
      ENDIF
      CAVG=CAVG+DBLE(POWC(ICH)*VOLCH)
   70 CONTINUE
      POWER=POWER/(10**3)
      CAVG=CAVG/VTOT
      CFACT=REAL(CAVG)/PCMX
*----
*  THERMAL TO FISSION RATIO POWER
*----
      IF(LFSTH) THEN
        CALL FLFSTH(PTOT,POWER,POWC,POWB,FLUX,NGRP,NCH,
     +              NB,NEL,FSTH,FLUB)
      ENDIF

      IF(IMPX.EQ.0)GOTO 90
*     CHECK TOTAL POWER
      IF(IMPX.EQ.99)WRITE(IOUT,1002)POWER
      IF(IMPX.GT.1)GOTO 80
      WRITE(TEXT,'(A9,I3.3)')'CHANNEL #',IPCH
      IF(PCMX.LT.10000.)THEN
        WRITE(IOUT,1003)PCMX,TEXT
      ELSE
        WRITE(IOUT,1013)PCMX,TEXT
      ENDIF
      IF(CAVG.LT.10000.)THEN
        WRITE(IOUT,1008)CAVG
      ELSE
        WRITE(IOUT,1014)CAVG
      ENDIF
      FACT=1./CFACT
      WRITE(IOUT,1010)CFACT,FACT
      GOTO 90
*----
*  PRINTING
*----
   80 JPMAP=LCMGID(IPMAP,'GEOMAP')
      CALL LCMGET(JPMAP,'STATE-VECTOR',IGST)
      NX=IGST(3)
      NY=IGST(4)
      NZ=IGST(5)
      IF(IGEO.NE.IGST(1)) CALL XABORT('@FLPOWB: WRONG GEOMETRY '
     1 // 'EMBEDDED IN THE FUEL MAP')
      IF(IGEO.EQ.7) THEN
        CALL FLPRNT(IPMAP,NCH,NB,NX,NY,NZ,POWB,PBMX,ICMX,
     1   IBMX,POWC,PCMX,IPCH,BAVG,BFACT,CAVG,CFACT,IMPX)
      ELSEIF(IGEO.EQ.9) THEN
        CALL FLPHPR(IPMAP,NCH,NB,NX,NZ,POWB,PBMX,ICMX,
     1   IBMX,POWC,PCMX,BAVG,BFACT,CAVG,CFACT,IMPX)
      ENDIF
   90 BFACT1=1./REAL(BFACT)
      CALL LCMPUT(IPPOW,'PMAX-BUND',1,2,PBMX)
      CALL LCMPUT(IPPOW,'FORM-BUND',1,2,BFACT1)
      CFACT1=1./CFACT
      CALL LCMPUT(IPPOW,'PMAX-CHAN',1,2,PCMX)
      CALL LCMPUT(IPPOW,'FORM-CHAN',1,2,CFACT1)
      RETURN
*
 1000 FORMAT(1X,'COMPUTED TOTAL POWER OVER ',
     1  'ALL BUNDLES =>',1P,E13.6,1X,'MW')
 1001 FORMAT(1X,'MAXIMUM BUNDLE POWER =',1X,F9.1,
     1  1X,'kW',2X,'=>',2X,A12,2X,'BUNDLE #',I2.2)
 1002 FORMAT(1X,'COMPUTED TOTAL POWER OVER',
     1  'ALL CHANNELS =>',1P,E13.6,1X,'MW')
 1003 FORMAT(1X,'MAXIMUM CHANNEL POWER =',1X,F9.1,
     1  1X,'kW',2X,'=>',2X,A12)
 1004 FORMAT(/1X,'** COMPUTING CHANNEL AND',
     1 1X,'BUNDLE POWERS **'/)
 1007 FORMAT(1X,'AVERAGE POWER OVER ALL BUNDLES',
     1 1X,'=',1X,F9.1,1X,'kW')
 1008 FORMAT(1X,'AVERAGE POWER OVER ALL CHANNELS',
     1 1X,'=',1X,F9.1,1X,'kW')
 1009 FORMAT(1X,'BUNDLE-POWER FORM FACTOR',2X,'=>',2X,
     1 'AVG/MAX =',1X,F8.4,3X,'(MAX/AVG = ',F8.4,')'/)
 1010 FORMAT(1X,'CHANNEL-POWER FORM FACTOR',2X,'=>',2X,
     1 'AVG/MAX =',1X,F8.4,3X,'(MAX/AVG = ',F8.4,')'/)
 1011 FORMAT(1X,'MAXIMUM BUNDLE POWER =',1X,F9.1,
     1  1X,'kW',2X,'=>',2X,A12,2X,'BUNDLE #',I2.2)
 1012 FORMAT(1X,'AVERAGE POWER OVER ALL BUNDLES',
     1 1X,'=',1X,F9.1,1X,'kW')
 1013 FORMAT(1X,'MAXIMUM CHANNEL POWER =',1X,F9.1,
     1  1X,'kW',2X,'=>',2X,A12)
 1014 FORMAT(1X,'AVERAGE POWER OVER ALL CHANNELS',
     1 1X,'=',1X,F9.1,1X,'kW')
      END
