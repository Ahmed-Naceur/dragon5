*DECK FLPOWR
      SUBROUTINE FLPOWR(NMIX,NGRP,NEL,LX,LY,LZ,MAT,VOL,FLUX,HFAC,PXYZ,
     1 VTOT,IMPX,LPOW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute and print a power distribution over the whole reactor core.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input
* NMIX   maximum number of material mixtures.
* NGRP   total number of energy groups.
* NEL    total number of finite elements.
* LX     number of elements along x-axis.
* LY     number of elements along y-axis.
* LZ     number of elements along z-axis.
* FLUX   normalized fluxes associated with each volume.
* MAT    index-number of mixture assigned to each volume.
* VOL    element-ordered mesh-splitted volumes.
* HFAC   h-factors over the reactor core.
* VTOT   total reactor core volume.
* IMPX   printing index (=0 for no print).
* LPOW   file printing flag: =.true. print on file.
*
*Parameters: output
* PXYZ   power distribution over the reactor core.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGRP,NMIX,NEL,LX,LY,LZ,MAT(NEL),IMPX
      REAL FLUX(NEL,NGRP),VOL(NEL),HFAC(NMIX,NGRP),PXYZ(LX,LY,LZ)
      DOUBLE PRECISION VTOT
      LOGICAL LPOW
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6,INIT=1)
      DOUBLE PRECISION PTOT
      CHARACTER TEXT*12
*----
*  CHECK TOTAL POWER
*----
      IF(IMPX.GT.0)WRITE(IOUT,1005)
      PTOT=0.0D0
      DO 20 IEL=1,NEL
      IF(MAT(IEL).EQ.0)GOTO 20
      DO 10 JGR=1,NGRP
      PTOT=PTOT+FLUX(IEL,JGR)*VOL(IEL)*HFAC(MAT(IEL),JGR)
   10 CONTINUE
   20 CONTINUE
      PAVG=REAL(PTOT/VTOT)
      IF(IMPX.GT.0)WRITE(IOUT,1001)PTOT,PAVG
*----
*  PERFORM CALCULATION
*----
      CALL XDRSET(PXYZ,NEL,0.)
      IEL=0
      PMAX=0.
      DO 52 K=1,LZ
      DO 51 J=1,LY
      DO 50 I=1,LX
      IEL=IEL+1
      IF(MAT(IEL).EQ.0)GOTO 50
      DO 40 JGR=1,NGRP
      PXYZ(I,J,K)=PXYZ(I,J,K)+
     1            HFAC(MAT(IEL),JGR)*FLUX(IEL,JGR)*VOL(IEL)
   40 CONTINUE
      IF(PXYZ(I,J,K).GT.PMAX)THEN
        PMAX=PXYZ(I,J,K)
        IMX=I
        JMX=J
        KMX=K
      ENDIF
   50 CONTINUE
   51 CONTINUE
   52 CONTINUE
      IF(IMPX.GT.0)WRITE(IOUT,1000)PMAX,IMX,JMX,KMX
      IF(.NOT.LPOW)GOTO 70
*----
*  PRINTING
*----
      TEXT='Pdistr.res'
      OPEN(UNIT=INIT,FILE=TEXT,STATUS='UNKNOWN')
      WRITE(INIT,1008)LX,LY,LZ
      DO 65 K=1,LZ
      DO 60 J=1,LY
      WRITE(INIT,1007)J,K
      WRITE(INIT,1002) (PXYZ(I,J,K),I=1,LX)
   60 CONTINUE
   65 CONTINUE
      CLOSE(UNIT=INIT)
      IF(IMPX.GT.0)WRITE(IOUT,1006)TEXT
   70 RETURN
*
 1000 FORMAT(/1X,'MAX POWER =',1P,E13.6,1X,'WATTS',4X,
     1 'AT COORD :',1X,'I =',I3,2X,'J =',I3,2X,'K =',I3/)
 1001 FORMAT(1X,'COMPUTED TOTAL POWER :',1P,E15.8,1X,'WATTS'/
     1 1X,'MEAN POWER DENSITY',3X,':',1P,E15.8,1X,'WATTS/CM3')
 1002 FORMAT(6(1P,E15.8))
 1005 FORMAT(/1X,'** COMPUTING POWER-DISTRIBUTION OVER',
     1 1X,'THE REACTOR CORE **'/)
 1006 FORMAT(/1X,'PRINTING POWER-DISTRIBUTION ON FILE:',
     1 1X,'<',A10,'>',3X,'=>',2X,'DONE.')
 1007 FORMAT(//3X,'PLANE-Y #',I2.2,5X,'PLANE-Z #',I2.2/)
 1008 FORMAT(/10X,5('*'),3X,'POWER-DISTRIBUTION OVER THE',
     1 1X,'REACTOR CORE',3X,5('*')//25X,'NX=',I2,',',2X,
     2 'NY=',I2,',',2X,'NZ=',I2)
      END
