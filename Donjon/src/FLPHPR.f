*DECK FLPHPR
      SUBROUTINE FLPHPR(IPMAP,NCH,NB,NX,NZ,POWB,PBNM,ICHM,IBNM,POWC,
     1 PCHM,BAVG,BFACT,CAVG,CFACT,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print the bundle and channel powers over the fuel lattice. Adapted
* from FLPRNT.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): 
* V. Descotes
*
*Parameters: input/output
* IPMAP  pointer to fuel-map information.
* NCH    number of reactor channels.
* NB     number of fuel bundles per channel.
* NX     number of elements along x-axis.
* NZ     number of elements along z-axis.
* POWB   bundle powers in kW.
* PBNM   maximum bundle power.
* ICHM   maximum-power channel number.
* IBNM   maximum-power bundle number.
* POWC   channel powers in kW.
* PCHM   maximum channel power.
* BAVG   average bundle power.
* BFACT  bundle power-form factor.
* CAVG   average channel power.
* CFACT  channel power-form factor.
* IMPX   printing index:  0 = no print
*                         1 = minimal printing
*                         2 = channel power only
*                         3 = bundle power by plane only
*                         10 = bundle power by channel
*                         any added values of 2, 3 and 10: 5,12,13,15
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER NCH,NB,NX,NZ,ICHM,IBNM,IMPX
      REAL POWB(NCH,NB),POWC(NCH),PBNM,PCHM
      DOUBLE PRECISION BAVG,CAVG,BFACT
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      REAL RADB(NX,NB),RADC(NX)
      INTEGER MIX(NX,NZ)
      CHARACTER TEXT*12,
     1          TEXT1A*17,TEXT2A*17,TEXT1B*17,TEXT2B*17
*
      CALL XDISET(MIX,NX*NZ,0)
      CALL LCMGET(IPMAP,'BMIX',MIX)
*----
*  BUNDLE POWERS OVER EACH CHANNEL
*----
      WRITE(IOUT,1009)
      IEL=0
      ICH=0
      DO 10 I=1,NX
      IEL=IEL+1
      DO 5 K=1,NZ
      IF(MIX(IEL,K).NE.0)GOTO 6
    5 CONTINUE
      GO TO 10
    6 ICH=ICH+1
      IF(IMPX.GE.10) THEN
        WRITE(TEXT,'(A9,I3.3)')'CHANNEL #',ICH
        IF(PCHM.LT.10000.)THEN
          WRITE(IOUT,1000)TEXT,POWC(ICH)
        ELSE
          WRITE(IOUT,1014)TEXT,POWC(ICH)
        ENDIF
        IF(PBNM.LT.1000.)THEN
          WRITE(IOUT,1001)(POWB(ICH,IB),IB=1,NB)
        ELSE
          WRITE(IOUT,1015)(POWB(ICH,IB),IB=1,NB)
        ENDIF
      ENDIF
   10 CONTINUE
*
      WRITE(TEXT1A,'(A6,I2,A8)')'(A5,',(NX/2),'(A5,3X))'
      WRITE(TEXT2A,'(A4,I2,A10)')'(A2,',(NX/2),'(F7.1,1X))'
      WRITE(TEXT1B,'(A6,I2,A8)')'(A5,',NX-(NX/2),'(A5,3X))'
      WRITE(TEXT2B,'(A4,I2,A10)')'(A2,',NX-(NX/2),'(F7.1,1X))'
      IF((IMPX.LT.3).OR.((IMPX.GE.10).AND.(IMPX.LT.13)))GOTO 50
*----
*  BUNDLE POWERS PER RADIAL PLANE
*----
      CALL XDRSET(RADB,NX*NB,0.)
      WRITE(IOUT,1010)
      DO IB=1,NB
        IEL=0
        ICH=0
        DO 20 I=1,NX
        IEL=IEL+1
        DO 15 K=1,NZ
        IF(MIX(IEL,K).NE.0)GOTO 16
   15   CONTINUE
        GO TO 20
   16   ICH=ICH+1
        RADB(I,IB)=POWB(ICH,IB)
   20   CONTINUE
      ENDDO
      DO IB=1,NB
        WRITE(IOUT,1011)IB

      ENDDO
   50 IF((IMPX.EQ.0).OR.(IMPX.EQ.1).OR.(IMPX.EQ.3).OR.(IMPX.EQ.4)
     1 .OR.(IMPX.EQ.10).OR.(IMPX.EQ.11).OR.(IMPX.EQ.13).OR.(IMPX.EQ.14))
     2  GOTO 90
*----
*  CHANNEL POWERS IN RADIAL PLANE
*----
      CALL XDRSET(RADC,NX,0.)
      WRITE(IOUT,1013)
      IEL=0
      ICH=0
      DO 60 I=1,NX
      IEL=IEL+1
      DO 55 K=1,NZ
      IF(MIX(IEL,K).NE.0)GOTO 56
   55 CONTINUE
      GO TO 60
   56 ICH=ICH+1
      RADC(I)=POWC(ICH)
   60 CONTINUE
      WRITE(IOUT,*)' '
*----
*  FINAL INFORMATION
*----
   90 WRITE(IOUT,1002)      
      IF(PBNM.LT.1000.)THEN
        WRITE(IOUT,1003)PBNM,IBNM
      ELSE
        WRITE(IOUT,1016)PBNM,ICHM,IBNM
      ENDIF
      IF(BAVG.LT.1000.)THEN
        WRITE(IOUT,1005)BAVG
      ELSE
        WRITE(IOUT,1017)BAVG
      ENDIF
      FACT=1./REAL(BFACT)
      WRITE(IOUT,1006)BFACT,FACT
      IF(PCHM.LT.10000.)THEN
        WRITE(IOUT,1004)PCHM,ICHM
      ELSE
        WRITE(IOUT,1018)PCHM,ICHM
      ENDIF
      IF(CAVG.LT.10000.)THEN
        WRITE(IOUT,1007)CAVG
      ELSE
        WRITE(IOUT,1019)CAVG
      ENDIF
      FACT=1./CFACT
      WRITE(IOUT,1008)CFACT,FACT
      RETURN
*
 1000 FORMAT(/5X,A12,10X,'CHANNEL POWER =',1X,F9.1,1X,'kW')
 1001 FORMAT(1X,12F8.3)
 1002 FORMAT(/5X,5('--o--',6X)/)
 1003 FORMAT(/1X,'MAXIMUM BUNDLE POWER  =',1X,F9.1,1X,'kW',
     1 3X,'=>',3X,'BUNDLE #',I2.2)
 1004 FORMAT(/1X,'MAXIMUM CHANNEL POWER =',1X,F9.1,1X,'kW',
     1 3X,'=>',3X,'CHANNEL #',I2,3X)
 1005 FORMAT(1X,'AVERAGE POWER OVER ALL BUNDLES',
     1 1X,'=',1X,F9.1,1X,'kW')
 1006 FORMAT(1X,'BUNDLE-POWER FORM FACTOR',2X,'=>',2X,
     1 'AVG/MAX =',1X,F8.4,3X,'(MAX/AVG = ',F8.4,')')
 1007 FORMAT(1X,'AVERAGE POWER OVER ALL CHANNELS',
     1 1X,'=',1X,F9.1,1X,'kW')
 1008 FORMAT(1X,'CHANNEL-POWER FORM FACTOR',2X,'=>',2X,
     1 'AVG/MAX =',1X,F8.4,2X,'(MAX/AVG = ',F8.4,')'/)
 1009 FORMAT(/20X,'** BUNDLE POWERS OVER EACH',
     1 1X,'CHANNEL (kW) **'/)
 1010 FORMAT(//20X,'** BUNDLE POWERS PER RADIAL',
     1 1X,'PLANE **'/)
 1011 FORMAT(//1X,'BUNDLE POWERS',1X,'(kW)',1X,
     1 '=>',1X,'RADIAL PLANE',1X,'#',I2.2/)
 1013 FORMAT(//20X,'** CHANNEL POWERS IN RADIAL',
     1 1X,'PLANE (kW) **'/)
 1014 FORMAT(/5X,A12,5X,
     1 'CHANNEL POWER =',1X,F9.1,1X,'kW')
 1015 FORMAT(1X,12(F9.1,1X))
 1016 FORMAT(/1X,'MAXIMUM BUNDLE POWER  =',1X,F9.1,1X,'kW',
     1 3X,'=>',3X,'CHANNEL #',1X,I2,3X,'BUNDLE #',I2.2)
 1017 FORMAT(1X,'AVERAGE POWER OVER ALL BUNDLES',
     1 1X,'=',1X,F9.1,1X,'kW')
 1018 FORMAT(/1X,'MAXIMUM CHANNEL POWER =',1X,F9.1,1X,'kW',
     1 3X,'=>',3X,'CHANNEL #',1X,I2)
 1019 FORMAT(1X,'AVERAGE POWER OVER ALL CHANNELS',
     1 1X,'=',1X,F9.1,1X,'kW')
      END
