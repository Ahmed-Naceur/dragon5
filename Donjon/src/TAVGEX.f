*DECK TAVGEX
      SUBROUTINE TAVGEX(IPMAP,IPPOW,NCH,NCOMB,NX,NY,NZ,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the core-average exit burnup and channel refuelling rates.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input/output
* IPMAP  pointer to fuel-map information.
* IPPOW  pointer to power information.
* NCH    number of reactor channels.
* NCOMB  number of combustion zones.
* NX     number of elements along x-axis in fuel map.
* NY     number of elements along y-axis in fuel map.
* NZ     number of elements along z-axis in fuel map.
* IMPX   printing index (=0 for no print).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP,IPPOW
      INTEGER NCH,NCOMB,NX,NY,NZ,IMPX
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      INTEGER MIX(NX*NY*NZ),BZONE(NCH),NSCH(NCH),NAMX(NX),NAMY(NY)
      REAL BVAL(NCOMB),RATE(NCH),POWC(NCH)
      DOUBLE PRECISION SUMR,SUMB
      CHARACTER TEXT*12,CHANX*2,CHANY*2
*----
*  RECOVER INFORMATION
*----
      CALL XDISET(MIX,NX*NY*NZ,0)
      CALL LCMGET(IPMAP,'BMIX',MIX)
*     CHANNEL POWERS
      CALL XDRSET(POWC,NCH,0.)
      CALL LCMGET(IPPOW,'POWER-CHAN',POWC)
*     REFUELLING SCHEME
      CALL XDISET(NSCH,NCH,0)
      CALL LCMGET(IPMAP,'REF-SCHEME',NSCH)
*     AVERAGE EXIT BURNUPS
      CALL XDRSET(BVAL,NCOMB,0.)
      CALL LCMGET(IPMAP,'BURN-AVG',BVAL)
*     COMBUSTION-ZONE INDEX
      CALL XDISET(BZONE,NCH,0)
      CALL LCMGET(IPMAP,'B-ZONE',BZONE)
*     CHANNEL NAMES
      CALL XDISET(NAMX,NX,0)
      CALL LCMGET(IPMAP,'XNAME',NAMX)
      CALL XDISET(NAMY,NY,0)
      CALL LCMGET(IPMAP,'YNAME',NAMY)
*----
*  CALCULATION OVER EACH CHANNEL
*----
      IF(IMPX.GT.0)WRITE(IOUT,1000)
      CALL XDRSET(RATE,NCH,0.)
      IEL=0
      ICH=0
      SUMR=0.0D0
      SUMB=0.0D0
      DO 15 J=1,NY
      DO 10 I=1,NX
      IEL=IEL+1
      IF(MIX(IEL).EQ.0)GOTO 10
      ICH=ICH+1
*     REFUELLING RATE
      RATE(ICH)=POWC(ICH)/BVAL(BZONE(ICH))
      SUMR=SUMR+RATE(ICH)
      SUMB=SUMB+BVAL(BZONE(ICH))*RATE(ICH)
      IF(IMPX.LT.4)GOTO 10
*     PRINT RATE
      WRITE(TEXT,'(A9,I3.3)')'CHANNEL #',ICH
      WRITE(CHANX,'(A2)') (NAMX(I))
      WRITE(CHANY,'(A2)') (NAMY(J))
      WRITE(IOUT,1001)TEXT,CHANY,CHANX,NSCH(ICH),RATE(ICH)
   10 CONTINUE
   15 CONTINUE
*     EXIT BURNUP
      BEXIT=REAL(SUMB/SUMR)
      IF(IMPX.EQ.0)GOTO 20
      IF(BEXIT.LT.10000.)THEN
        WRITE(IOUT,1002)BEXIT
      ELSE
        WRITE(IOUT,1003)BEXIT
      ENDIF
   20 CALL LCMPUT(IPMAP,'B-EXIT',1,2,BEXIT)
      CALL LCMPUT(IPMAP,'REF-RATE',NCH,2,RATE)
      RETURN
*
 1000 FORMAT(/1X,'**',1X,'COMPUTING CHANNEL',
     1 1X,'REFUELLING',1X,'RATES',1X,'**'/)
 1001 FORMAT(5X,A12,5X,'NAME:',1X,A2,A2,5X,'RE',
     1 'F-SCHEME:',1X,I2,5X,'REF-RATE: ',F6.4/)
 1002 FORMAT(/1X,'CORE-AVERAGE EXIT BURNUP',
     1  1X,'=',1X,F7.2,1X,'MW*DAY/T'/)
 1003 FORMAT(/1X,'CORE-AVERAGE EXIT BURNUP',
     1  1X,'=',1X,F8.2,1X,'MW*DAY/T'/)
      END
