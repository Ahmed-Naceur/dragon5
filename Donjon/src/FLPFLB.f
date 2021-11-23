*DECK FLPFLB
      SUBROUTINE FLPFLB(IPMTX,NMAT,NGRP,NEL,NCH,NB,FLUX,VOL,FMIX,VOLB,
     1 FLXB,IMPX,LMAP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute average fluxes per fuel bundle and other related quantities;
* print bundle fluxes on file.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input
* IPMTX  pointer to matex information.
* NMAT   total number of mixtures (includes virtual regions).
* NGRP   number of energy groups.
* NEL    total number of elements.
* NCH    number of reactor channels.
* NB     number of fuel bundles per channel.
* FLUX   normalized fluxes associated with each volume.
* VOL    element-ordered mesh-splitted volumes.
* IMPX   printing index (=0 for no print).
* LMAP   flux printing flag (=.true. print on file).
* FMIX   fuel bundle indices.
*
*Parameters: output
* VOLB   bundle volumes.
* FLXB   bundle fluxes.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMTX
      INTEGER NMAT,NEL,NGRP,NCH,NB,IMPX,FMIX(NCH*NB)
      REAL FLUX(NEL,NGRP),VOL(NEL),VOLB(NCH,NB),FLXB(NCH,NB,NGRP)
      LOGICAL LMAP
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6,INIT=1)
      INTEGER FMAT(NMAT)
      CHARACTER TEXT*12,FNAM*12
      REAL RATIO(NCH,NB,NGRP-1)
      DOUBLE PRECISION FAVG(NGRP)
*----
*  PERFORM CALCULATION
*----
      FMAX=0.
      ICM=0
      IBM=0
      MGR=0
      CALL XDDSET(FAVG,NGRP,0.0D0)
      CALL XDISET(FMAT,NMAT,0)
      CALL LCMGET(IPMTX,'MAT',FMAT)
      CALL XDRSET(FLXB,NCH*NB*NGRP,0.)
      IF(IMPX.GT.0)WRITE(IOUT,1004)
      NTOT=0
      VTOT=0.0
      DO 45 IB=1,NB
      DO 40 ICH=1,NCH
      NUM=(IB-1)*NCH+ICH
      VOLB(ICH,IB)=0.0
      IF(FMIX(NUM).EQ.0) GO TO 40
      NTOT=NTOT+1
      DO 20 IEL=1,NEL
      IF(FMAT(IEL).NE.-NTOT)GOTO 20
      DO 10 JGR=1,NGRP
      FLXB(ICH,IB,JGR)=FLXB(ICH,IB,JGR)+FLUX(IEL,JGR)*VOL(IEL)
   10 CONTINUE
      VOLB(ICH,IB)=VOLB(ICH,IB)+VOL(IEL)
   20 CONTINUE
      DO JGR=1,NGRP
        FLXB(ICH,IB,JGR)=FLXB(ICH,IB,JGR)/VOLB(ICH,IB)
        IF(ABS(FLXB(ICH,IB,JGR)).GT.FMAX)THEN
          FMAX=FLXB(ICH,IB,JGR)
          ICM=ICH
          IBM=IB
          MGR=JGR
        ENDIF
        FAVG(JGR)=FAVG(JGR)+FLXB(ICH,IB,JGR)*VOLB(ICH,IB)
      ENDDO
      VTOT=VTOT+VOLB(ICH,IB)
   40 CONTINUE
   45 CONTINUE
*     MAX AND CORE-AVERAGE FLUXES
      IF(IMPX.GT.0)WRITE(IOUT,1007)FMAX,ICM,IBM,MGR
      DO JGR=1,NGRP
        FAVG(JGR)=FAVG(JGR)/VTOT
        IF(IMPX.GT.0)WRITE(IOUT,1008)FAVG(JGR),JGR
      ENDDO
*     FORM FACTOR
      IF(MGR.EQ.0) CALL XABORT('FLPFLB: FLUX NORMALIZATION FAILURE.')
      FACT=REAL(FAVG(MGR))/FMAX
      FACT2=1./FACT
      IF(IMPX.GT.0)WRITE(IOUT,1009)MGR,FACT,FACT2,VTOT
*     FLUXES RATIOS
      CALL XDRSET(RATIO,NCH*NB*(NGRP-1),0.)
      DO 52 IB=1,NB
      DO 51 ICH=1,NCH
      DO 50 JGR=1,NGRP-1
      RATIO(ICH,IB,JGR)=FLXB(ICH,IB,JGR)/FLXB(ICH,IB,NGRP)
   50 CONTINUE
   51 CONTINUE
   52 CONTINUE
      IF(.NOT.LMAP)GOTO 80
*----
*  PRINTING
*----
      FNAM='FluxMAP.res'
      OPEN(UNIT=INIT,FILE=FNAM,STATUS='UNKNOWN')
      WRITE(INIT,1000)NCH,NB,NGRP
      DO 65 JGR=1,NGRP
      WRITE(INIT,1001)JGR
      DO 60 ICH=1,NCH
      WRITE(TEXT,'(A9,I3.3)')'CHANNEL #',ICH
      WRITE(INIT,1002)TEXT
      WRITE(INIT,1003)(FLXB(ICH,IB,JGR),IB=1,NB)
   60 CONTINUE
   65 CONTINUE
      WRITE(INIT,1010)
      DO 75 JGR=1,NGRP-1
      WRITE(INIT,1011)JGR,NGRP
      DO 70 ICH=1,NCH
      WRITE(TEXT,'(A9,I3.3)')'CHANNEL #',ICH
      WRITE(INIT,1002)TEXT
      WRITE(INIT,1012)(RATIO(ICH,IB,JGR),IB=1,NB)
   70 CONTINUE
   75 CONTINUE
      CLOSE(UNIT=INIT)
      IF(IMPX.GT.0)WRITE(IOUT,1006)FNAM
   80 RETURN
*
 1000 FORMAT(/20X,5('*'),3X,'AVERAGE FUEL-BUNDLES ',
     1 'FLUXES',3X,5('*')//5X,'NUMBER OF CHANNELS:',
     2 1X,I3,4X,'NUMBER OF BUNDLES:',1X,I2,4X,
     3 'NUMBER OF GROUPS:',I2)
 1001 FORMAT(//18X,'ENERGY GROUP =>',1X,I2.2)
 1002 FORMAT(/1X,A12)
 1003 FORMAT(6(1P,E15.8))
 1004 FORMAT(/1X,'** COMPUTING AVERAGE',1X,'BUNDLE FLUXES **'/)
 1006 FORMAT(/1X,'PRINTING BUNDLE FLUXES ON FILE:',
     1 1X,'<',A11,'>',3X,'=>',2X,'DONE.')
 1007 FORMAT(1X,'MAX FLUX =',1P,E13.6,2X,'=>',
     1 2X,'CHANNEL #',I3.3,2X,'BUNDLE #',I2.2,
     2 2X,'GROUP #',I2.2/)
 1008 FORMAT(1X,'FUEL-ZONE AVERAGE FLUX =',
     1 1P,E13.6,3X,'=>',2X,'GROUP #',I2.2)
 1009 FORMAT(/1X,'FLUX-FORM FACTOR FOR GROUP #',I2.2,
     1 2X,'=>',2X,'AVG/MAX = ',F8.4,2X,'(MAX/AVG = ',
     2 F8.4,')'/' FUEL-ZONE VOLUME =',1P,E13.6,' CM3'/)
 1010 FORMAT(//16X,5('*'),3X,'FUEL-BUNDLES',
     1 1X,'FLUXES RATIOS',3X,5('*')/)
 1011 FORMAT(/18X,'FLUX RATIO: GROUP #',I2.2,
     1 1X,'=>',1X,'GROUP #',I2.2)
 1012 FORMAT(6(1P,E13.6))
      END
