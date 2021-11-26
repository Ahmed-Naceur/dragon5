*DECK FLPFLX
      SUBROUTINE FLPFLX(NGRP,NEL,LX,LY,LZ,MAT,VOL,FLUX,FXYZ,RATIO,VTOT,
     1 IMPX,LFLX,LRAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the flux distributions and ratios over the whole reactor core;
* print the normalized fluxes on files.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input
* NGRP   total number of energy groups.
* NEL    total number of finite elements.
* LX     number of elements along x-axis.
* LY     number of elements along y-axis.
* LZ     number of elements along z-axis.
* MAT    index-number of mixture assigned to each volume.
* VOL    element-ordered mesh-splitted volumes.
* FLUX   normalized fluxes associated with each volume.
* IMPX   screen printing index (=0 for no print).
* LFLX   fluxes printing flag: =.true. print on files.
* LRAT   ratios printing flag: =.true. print on files.
*
*Parameters: output
* FXYZ   mesh-ordered fluxes.
* RATIO  fluxes ratios with respect to thermal fluxes.
* VTOT   total reactor-core volume.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGRP,NEL,LX,LY,LZ,MAT(NEL),IMPX
      REAL FXYZ(LX,LY,LZ,NGRP),FLUX(NEL,NGRP),
     1     RATIO(LX,LY,LZ,NGRP-1),VOL(NEL)
      DOUBLE PRECISION VTOT
      LOGICAL LFLX,LRAT
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6,INIT=1)
      CHARACTER TEXT*12
      DOUBLE PRECISION FAVG(NGRP)
*----
*  PERFORM CALCULATION
*----
      IF(IMPX.GT.0)WRITE(IOUT,1010)
      FMAX=0.
      IM=0
      KM=0
      MGR=0
      VTOT=0.0D0
      CALL XDDSET(FAVG,NGRP,0.0D0)
      CALL XDRSET(FXYZ,NEL*NGRP,0.)
      IEL=0
      DO 12 K=1,LZ
      DO 11 J=1,LY
      DO 10 I=1,LX
      IEL=IEL+1
      IF(MAT(IEL).EQ.0)GOTO 10
      VTOT=VTOT+VOL(IEL)
      DO JGR=1,NGRP
        FXYZ(I,J,K,JGR)=FLUX(IEL,JGR)
        IF(ABS(FXYZ(I,J,K,JGR)).GT.FMAX)THEN
          FMAX=FXYZ(I,J,K,JGR)
          IM=I
          JM=J
          KM=K
          MGR=JGR
        ENDIF
        FAVG(JGR)=FAVG(JGR)+FXYZ(I,J,K,JGR)*VOL(IEL)
      ENDDO
   10 CONTINUE
   11 CONTINUE
   12 CONTINUE
*     MAX AND CORE-AVERAGE FLUXES
      IF(IMPX.GT.0)WRITE(IOUT,1000)FMAX,IM,JM,KM,MGR
      DO JGR=1,NGRP
        FAVG(JGR)=FAVG(JGR)/VTOT
        IF(IMPX.GT.0)WRITE(IOUT,1001)FAVG(JGR),JGR
      ENDDO
      IF(MGR.EQ.0) CALL XABORT('FLPFLX: FLUX NORMALIZATION FAILURE.')
      FACT=REAL(FAVG(MGR))/FMAX
      FACT2=1./FACT 
      IF(IMPX.GT.0)WRITE(IOUT,1002)MGR,FACT,FACT2,VTOT
*     FLUXES RATIOS
      CALL XDRSET(RATIO,NEL*(NGRP-1),0.)
      DO 32 K=1,LZ
      DO 31 J=1,LY
      DO 30 I=1,LX
      IF(FXYZ(I,J,K,NGRP).EQ.0.)GOTO 30
      DO 20 JGR=1,NGRP-1
      RATIO(I,J,K,JGR)=FXYZ(I,J,K,JGR)/FXYZ(I,J,K,NGRP)
   20 CONTINUE
   30 CONTINUE
   31 CONTINUE
   32 CONTINUE
      IF(.NOT.LFLX)GOTO 60
*----
*  PRINTING
*----
      IF(IMPX.GT.0)WRITE(IOUT,1006)
*     FLUXES
      DO 50 JGR=1,NGRP
        WRITE(TEXT,'(A4,I2.2,A4)')'Flux',JGR,'.res'
        OPEN(UNIT=INIT,FILE=TEXT,STATUS='UNKNOWN')
        WRITE(INIT,1011)LX,LY,LZ,JGR
        DO 40 K=1,LZ
        DO J=1,LY
        WRITE(INIT,1009)J,K
        WRITE(INIT,1005) (FXYZ(I,J,K,JGR),I=1,LX)
        ENDDO
   40   CONTINUE
        CLOSE(UNIT=INIT)
        IF(IMPX.GT.0)WRITE(IOUT,1003)JGR,TEXT
   50 CONTINUE
*
   60 IF(.NOT.LRAT)GOTO 90
      IF(IMPX.GT.0)WRITE(IOUT,1007)
*     RATIOS
      DO 80 JGR=1,NGRP-1
        WRITE(TEXT,'(A4,I2.2,A4)')'Rati',JGR,'.res'
        OPEN(UNIT=INIT,FILE=TEXT,STATUS='UNKNOWN')
        WRITE(INIT,1012)JGR,NGRP,LX,LY,LZ
        DO 70 K=1,LZ
        DO J=1,LY
        WRITE(INIT,1009)J,K
        WRITE(INIT,1008) (RATIO(I,J,K,JGR),I=1,LX)
        ENDDO
   70   CONTINUE
        CLOSE(UNIT=INIT)
        IF(IMPX.GT.0)WRITE(IOUT,1004)JGR,NGRP,TEXT
   80 CONTINUE
   90 RETURN
*
 1000 FORMAT(1X,'MAX FLUX =',1P,E13.6,4X,'AT COORD :',1X,
     1 'I =',I3,2X,'J =',I3,2X,'K =',I3,2X,'GROUP #',I2.2/)
 1001 FORMAT(1X,'CORE-AVERAGE FLUX =',1P,E13.6,
     1 2X,'=>',2X,'GROUP #',I2.2)
 1002 FORMAT(/1X,'OVERALL FLUX-FORM FACTOR FOR GROUP #',I2.2,
     1 2X,'=>',2X,'AVG/MAX =',1X,F8.4,2X,'(MAX/AVG = ',F8.4,
     2 ')'/1X,'TOTAL CORE VOLUME =',1P,E13.6,1X,'CM3'/)
 1003 FORMAT(1X,'FLUXES',2X,'=>',2X,'GROUP #',I2.2,2X,
     1 '=>',2X,'FILE NAME: <',A10,'>',2X,'=>',2X,'DONE.')
 1004 FORMAT(1X,'FLUX RATIOS',2X,'=>',2X,'GR.#',I2.2,
     1 '/GR.#',I2.2,2X,'=>',2X,'FILE NAME: <',A10,'>',
     1 2X,'=>',2X,'DONE.')
 1005 FORMAT(1X,1P,6E16.8)
 1006 FORMAT(/15X,'** PRINTING OF FLUXES ON FILES **'/)
 1007 FORMAT(/15X,'** PRINTING OF RATIOS ON FILES **'/)
 1008 FORMAT(1X,1P,6E14.6)
 1009 FORMAT(//3X,'PLANE-Y #',I2.2,5X,'PLANE-Z #',I2.2/)
 1010 FORMAT(/1X,'** COMPUTING FLUX-DISTRIBUTION',
     1 1X,'OVER THE REACTOR CORE **'/)
 1011 FORMAT(/10X,5('*'),3X,'FLUX-DISTRIBUTION OVER THE',
     1 1X,'REACTOR CORE',3X,5('*')//21X,'NX=',I2,',',2X,
     2 'NY=',I2,',',2X,'NZ=',I2,',',2X,'GROUP #',I2.2)
 1012 FORMAT(/10X,5('*'),3X,'FLUXES RATIO',1X,'#',I2.2,
     1 '/#',I2.2,1X,'OVER THE REACTOR CORE',3X,5('*')//
     2 25X,'NX=',I2,',',2X,'NY=',I2,',',2X,'NZ=',I2)
      END
