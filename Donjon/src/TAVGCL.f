*DECK TAVGCL
      SUBROUTINE TAVGCL(IPMAP,IPPOW,NCH,NB,NCOMB,NX,NY,NZ,ARP,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute burnup limits over the fuel lattice for the time-average
* integration, based on the axial power shape over each channel.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
*
*Author(s): 
* D. Sekki, R. Chambon
*
*Parameters: input/output
* IPMAP   pointer to fuel-map information.
* IPPOW   pointer to power information.
* NCH     number of reactor channels.
* NB      number of fuel bundles per channel.
* NCOMB   number of combustion zones.
* NX      number of elements along x-axis in fuel map.
* NY      number of elements along y-axis in fuel map.
* NZ      number of elements along z-axis in fuel map.
* ARP     relaxation parameter for shape convergence.
* IMPX    printing index (=0 for no print).
*
*Parameters: scratch
* BURN0   low burnup integration limits.
* BURN1   upper burnup integration limits.
*
*----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP,IPPOW
      INTEGER NCH,NB,NCOMB,NX,NY,NZ,IMPX
      REAL ARP
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      INTEGER MIX(NX*NY*NZ),NAMX(NX),NAMY(NY),
     1 IVECT(NCOMB,NB),NSCH(NCH),BZONE(NCH),IGAR(NB)
      REAL POWB(NCH,NB),POWC(NCH),PSI(NB),BVAL(NCOMB),SOLD(NCH,NB),
     1 BURN0(NCH,NB),BURN1(NCH,NB),B0(NB),B1(NB),SNEW(NCH,NB)
      CHARACTER TEXT*12,CHANX*2,CHANY*2
      DOUBLE PRECISION PNUM,PDEN
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ICHMAP
*----
*  RECOVER INFORMATION
*----
      CALL XDISET(MIX,NX*NY*NZ,0)
      CALL LCMGET(IPMAP,'BMIX',MIX)
*     CHANNEL NAMES
      CALL XDISET(NAMX,NX,0)
      CALL LCMGET(IPMAP,'XNAME',NAMX)
      CALL XDISET(NAMY,NY,0)
      CALL LCMGET(IPMAP,'YNAME',NAMY)
*     COMBUSTION-ZONE INDEX
      CALL XDISET(BZONE,NCH,0)
      CALL LCMGET(IPMAP,'B-ZONE',BZONE)
*     AVERAGE EXIT BURNUPS
      CALL XDRSET(BVAL,NCOMB,0.)
      CALL LCMGET(IPMAP,'BURN-AVG',BVAL)
*     REFUELLING SCHEME
      CALL XDISET(NSCH,NCH,0)
      CALL LCMGET(IPMAP,'REF-SCHEME',NSCH)
*     REFUELLING VECTOR
      CALL XDISET(IVECT,NCOMB*NB,0)
      CALL LCMGET(IPMAP,'REF-VECTOR',IVECT)
*     PREVIOUS AXIAL SHAPE
      CALL XDRSET(SOLD,NCH*NB,0.)
      CALL LCMGET(IPMAP,'AX-SHAPE',SOLD)
*     CHANNEL POWERS
      CALL XDRSET(POWC,NCH,0.)
      CALL LCMGET(IPPOW,'POWER-CHAN',POWC)
*     BUNDLE POWERS
      CALL XDRSET(POWB,NCH*NB,0.)
      CALL LCMGET(IPPOW,'POWER-BUND',POWB)
*----
*  SET THE CHANNEL INDEX MAP
*----
      ALLOCATE(ICHMAP(NX,NY))
      CALL XDISET(ICHMAP,NX*NY,0)
      ICH=0
      DO 35 J=1,NY
      DO 30 I=1,NX
      IEL=(J-1)*NX+I
      DO 10 IZ=1,NZ
      IF(MIX((IZ-1)*NX*NY+IEL).NE.0) GO TO 20
  10  CONTINUE
      GO TO 30
  20  ICH=ICH+1
      ICHMAP(I,J)=ICH
  30  CONTINUE
  35  CONTINUE
      IF(ICH.NE.NCH) CALL XABORT('@TAVGCL: INVALID NUMBER OF CHANNELS')
*----
*  CALCULATION OVER EACH CHANNEL
*----
      IF(IMPX.GT.0)WRITE(IOUT,1005)
      CALL XDRSET(BURN0,NCH*NB,0.)
      CALL XDRSET(BURN1,NCH*NB,0.)
      ICH=0
      PNUM=0.0D0
      PDEN=0.0D0
      DO 45 J=1,NY
      DO 40 I=1,NX
      IF(ICHMAP(I,J).EQ.0)GOTO 40
      ICH=ICH+1
*     POWER-SHAPE
      DO IB=1,NB
        IF(POWC(ICH).EQ.0.0) CALL XABORT('TAVGCL: ZERO CHANNEL POWER.')
        PSI(IB)=ARP*(POWB(ICH,IB)/POWC(ICH))+(1.-ARP)*SOLD(ICH,IB)
        SNEW(ICH,IB)=PSI(IB)
        PNUM=PNUM+(SNEW(ICH,IB)-SOLD(ICH,IB))**2
        PDEN=PDEN+SNEW(ICH,IB)**2
        IGAR(IB)=IVECT(BZONE(ICH),IB)
      ENDDO
      IBSH=ABS(NSCH(ICH))
*     INTEGRATION LIMITS
      CALL TAVGLM(NB,IBSH,BVAL(BZONE(ICH)),PSI,B0,B1,IGAR,NSCH(ICH))
      DO IB=1,NB
        BURN0(ICH,IB)=B0(IB)
        BURN1(ICH,IB)=B1(IB)
      ENDDO
      IF(IMPX.GE.3) THEN
*       PRINT BURNUP LIMITS
        WRITE(TEXT,'(A9,I3.3)')'CHANNEL #',ICH
        WRITE(CHANX,'(A2)') (NAMX(I))
        WRITE(CHANY,'(A2)') (NAMY(J))
        WRITE(IOUT,1000)TEXT,CHANY,CHANX,NSCH(ICH)
        WRITE(IOUT,1001)'B0',(B0(IB),IB=1,NB)
        WRITE(IOUT,1001)'B1',(B1(IB),IB=1,NB)
      ENDIF
   40 CONTINUE
   45 CONTINUE
*     AXIAL-SHAPE ERROR
      EPS=REAL(SQRT(PNUM/PDEN))
*----
*  PRINT INFORMATION
*----
      IF(IMPX.GT.0)WRITE(IOUT,1002)EPS,ARP
      IF(IMPX.GE.3) THEN
*       PRINT SHAPE
        WRITE(IOUT,1003)
        DO ICH=1,NCH
          WRITE(TEXT,'(A6,I3.3)')'CHAN #',ICH
          WRITE(IOUT,1004)TEXT,(SNEW(ICH,IB),IB=1,NB)
        ENDDO
      ENDIF
*----
*  STORE INFORMATION
*----
      CALL LCMPUT(IPMAP,'BURN-BEG',NCH*NB,2,BURN0)
      CALL LCMPUT(IPMAP,'BURN-END',NCH*NB,2,BURN1)
      CALL LCMPUT(IPMAP,'EPS-AX',1,2,EPS)
      CALL LCMPUT(IPMAP,'AX-SHAPE',NCH*NB,2,SNEW)
      DEALLOCATE(ICHMAP)
      RETURN
*
 1000 FORMAT(/5X,A12,5X,'NAME:',1X,A2,A2,
     1  5X,'REFUELLING SCHEME:',1X,I2)
 1001 FORMAT(A3,12(F8.1,1X))
 1002 FORMAT(1X,'AXIAL-SHAPE ERROR =>',1P,E13.6,5X,
     1 'RELAXATION PARAMETER =>',E13.6/)
 1003 FORMAT(/20X,'** AXIAL SHAPE OVER EACH',
     1  1X,'CHANNEL **'/)
 1004 FORMAT(1X,A10,(2X,12(F6.4,1X)))
 1005 FORMAT(/1X,'** COMPUTING BURNUP INTEG',
     1 'RATION',1X,'LIMITS **'/)
      END
