*DECK DOORAB
      SUBROUTINE DOORAB (CDOOR,JPSYS,NPSYS,IPTRK,IMPX,NGRP,NREG,NBMIX,
     1 NANI,MAT,VOL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Double heterogeneity treatment (part 1). Vectorial version.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* CDOOR   name of the geometry/solution operator.
* JPSYS   pointer to the PIJ LCM object (L_PIJ signature). JPSYS is a
*         list of NGRP directories.
* NPSYS   index array pointing to the JPSYS list component corresponding
*         to each energy group. Set to zero if a group is not to be
*         processed. Usually, NPSYS(I)=I.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IMPX    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* NBMIX   number of mixtures (NBMIX=max(MAT(i))).
* NANI    number of Legendre orders (usually equal to one).
*
*Parameters: input/output
* NREG    total number of regions for which specific values of the
*         neutron flux and reaction rates are required on input and
*         number of regions in the macro-geometry at output.
* MAT     index-number of the mixture type assigned to each volume
*         on input and index-number of the mixture type assigned 
*         to each macro-volume at output.
* VOL     volume on input and macro-volumes at output.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER CDOOR*12
      TYPE(C_PTR) JPSYS,IPTRK
      INTEGER NPSYS(NGRP)
      INTEGER IMPX,NGRP,NREG,NBMIX,NANI,MAT(NREG)
      REAL VOL(NREG)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      INTEGER IPAR(8)
      TYPE(C_PTR) KPSYS,KPBIH
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NS,IDIL,MIXGR,IBI,NCO
      REAL, ALLOCATABLE, DIMENSION(:) :: RS,FRACT,VOLK,VOL2,SGAR,SGAS,
     > RRRR,QKDEL,QKOLD,PKL,P1I,P1DI,P1KI,SIGA1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: COEF
*---
*  RECOVER DOUBLE HETEROGENEITY DATA
*----
      CALL LCMSIX(IPTRK,'BIHET',1)
      CALL LCMGET(IPTRK,'PARAM',IPAR)
      IR1=IPAR(1)
      IR2=IPAR(2)
      NREG2=IPAR(3)
      NG=IPAR(4)
      NSMAX=IPAR(5)
      IBIHET=IPAR(6)
      MICRO=IPAR(7)
      IQUA10=IPAR(8)
      IF(IR1.NE.NBMIX) CALL XABORT('DOORAB: INVALID DATA IN TRACKING.')
      IF(IBIHET.EQ.0) CALL XABORT('DOORAB: BIHET MODEL NOT SET.')
      NMILG=IR2-IR1
      ALLOCATE(NS(NG),IDIL(NMILG),MIXGR(NSMAX*NG*NMILG),IBI(NREG2))
      ALLOCATE(RS(NG*(1+NSMAX)),FRACT(NG*IR2),VOLK(NG*NSMAX),
     1 VOL2(NREG2))
      CALL LCMGET(IPTRK,'NS',NS)
      CALL LCMGET(IPTRK,'RS',RS)
      CALL LCMGET(IPTRK,'FRACT',FRACT)
      CALL LCMGET(IPTRK,'VOLK',VOLK)
      CALL LCMGET(IPTRK,'IDIL',IDIL)
      CALL LCMGET(IPTRK,'MIXGR',MIXGR)
      CALL LCMGET(IPTRK,'VOLUME',VOL2)
      CALL LCMGET(IPTRK,'IBI',IBI)
      CALL LCMGET(IPTRK,'FRTM',FRTM)
      CALL LCMSIX(IPTRK,' ',2)
      IF(IMPX.GT.1) THEN
        WRITE(IOUT,'(/43H DOORAB: RECOVER DOUBLE-HETEROGENEITY DATA.)')
        WRITE(IOUT,'(35H   NUMBER OF ORDINARY MIXTURES     ,I4)') IR1
        WRITE(IOUT,'(35H   NUMBER OF COMPOSITE MIXTURES    ,I4)') NMILG
        WRITE(IOUT,'(35H   NUMBER OF KIND OF GRAINS        ,I4)') NG
        WRITE(IOUT,'(35H   NUMBER OF ORDINARY VOLUMES      ,I4)') NREG2
        WRITE(IOUT,'(35H   NUMBER OF COMPOSITE VOLUMES     ,I4)') NREG
        WRITE(IOUT,'(35H   TYPE OF MICRO VOLUMES (3 OR 4)  ,I4)') MICRO
        WRITE(IOUT,'(35H   MAX. NUMBER OF VOLUMES PER GRAIN,I4)') NSMAX
        WRITE(IOUT,'(35H   QUADRATURE PARAMETER FOR GRAINS ,I4)') IQUA10
        WRITE(IOUT,'(35H   DOUBLE HETEROGENEITY MODEL      ,I4)') IBIHET
        IF(IBIHET.EQ.3.OR.IBIHET.EQ.4) THEN
          WRITE(IOUT,'(35H   MINIMUM GRAINS FRACTION         ,F8.4)') 
     >     FRTM
        ENDIF
      ENDIF
*----
*  COMPUTE THE EQUIVALENT CROSS SECTIONS IN COMPOSITE REGIONS
*----
      NB1=NBMIX+1
      ALLOCATE(SGAR((NB1+NMILG)),SGAS((NB1+NMILG)*NANI))
      CALL XDRSET(SGAS,(NB1+NMILG)*NANI,0.0)
      DO 100 IGR=1,NGRP
      IOFSET=NPSYS(IGR)
      IF(IOFSET.NE.0) THEN
         IF(IMPX.GT.10) WRITE(IOUT,'(/25H DOORAB: PROCESSING GROUP,I5,
     >   6H WITH ,A,1H.)') IGR,CDOOR
         KPSYS=LCMGIL(JPSYS,IOFSET)
         CALL LCMGET(KPSYS,'DRAGON-TXSC',SGAR)
         CALL LCMGET(KPSYS,'DRAGON-S0XSC',SGAS)
*----
*   MORE MEMORY ALLOCATION
*----
         IF((IBIHET.EQ.1).OR.(IBIHET.EQ.2)) THEN
           ALLOCATE(NCO(NMILG))
           ALLOCATE(RRRR(NMILG),QKDEL(NG*NSMAX*NMILG),
     >     QKOLD(NG*NSMAX*NMILG),PKL(NG*NSMAX*NSMAX*NMILG))
           ALLOCATE(COEF(NMILG*(1+NG*NSMAX)**2))
         ELSEIF ((IBIHET.EQ.3).OR.(IBIHET.EQ.4)) THEN
           ALLOCATE(P1I(NG*NMILG),P1KI(NSMAX*NG*NMILG),
     >     P1DI(NG*NMILG),SIGA1(NG*NMILG))
         ENDIF
*----
*   DOUBLE HETEROGENEITY TREATMENT -- PART 1
*----
         IF(IBIHET.EQ.1) THEN
*           SANCHEZ-POMRANING MODEL.
            CALL XDRH11(NBMIX,NMILG,NG,NSMAX,MICRO,IQUA10,NS,IDIL,
     >      MIXGR,RS,FRACT,VOLK,SGAR,SGAS,NCO,RRRR,QKOLD,QKDEL,PKL,
     >      COEF)
         ELSEIF(IBIHET.EQ.2) THEN
*           HEBERT MODEL.
            CALL XDRH12(NBMIX,NMILG,NG,NSMAX,MICRO,IQUA10,NS,IDIL,
     >      MIXGR,RS,FRACT,VOLK,SGAR,SGAS,NCO,RRRR,QKDEL,PKL,COEF)
         ELSEIF((IBIHET.EQ.3).OR.(IBIHET.EQ.4)) THEN
*           SHE-LIU-SHI MODEL.
            CALL XDRH13(NBMIX,NMILG,NG,NSMAX,IQUA10,FRTM,NS,IDIL,
     >      MIXGR,RS,FRACT,SGAR,SGAS,P1I,P1DI,P1KI,SIGA1)
         ELSE
            CALL XABORT('DOORAB: INVALID DOUBLE HETEROGENEITY MODEL.')
         ENDIF
*
         KPBIH=LCMDID(KPSYS,'BIHET')
         CALL LCMPUT(KPBIH,'DRAGON-TXSC',1+IR2,2,SGAR)
         CALL LCMPUT(KPBIH,'DRAGON-S0XSC',(1+IR2)*NANI,2,SGAS)
*
         IF((IBIHET.EQ.1).OR.(IBIHET.EQ.2)) THEN
           CALL LCMPUT(KPSYS,'NCO',NMILG,1,NCO)
           CALL LCMPUT(KPSYS,'RRRR',NMILG,2,RRRR)
           CALL LCMPUT(KPSYS,'QKOLD',NG*NSMAX*NMILG,2,QKOLD)
           CALL LCMPUT(KPSYS,'QKDEL',NG*NSMAX*NMILG,2,QKDEL)
           CALL LCMPUT(KPSYS,'PKL',NG*NSMAX*NSMAX*NMILG,2,PKL)
           CALL LCMPUT(KPSYS,'COEF',NMILG*(1+NG*NSMAX)**2,4,COEF)
           DEALLOCATE(NCO)
           DEALLOCATE(COEF)
           DEALLOCATE(PKL,QKOLD,QKDEL,RRRR)
         ELSEIF((IBIHET.EQ.3).OR.(IBIHET.EQ.4)) THEN
           CALL LCMPUT(KPSYS,'P1I',NG*NMILG,2,P1I)
           CALL LCMPUT(KPSYS,'P1DI',NG*NMILG,2,P1DI)
           CALL LCMPUT(KPSYS,'P1KI',NG*NSMAX*NMILG,2,P1KI)
           CALL LCMPUT(KPSYS,'SIGA1',NG*NMILG,2,SIGA1)
           DEALLOCATE(P1I,P1DI,P1KI,SIGA1)
         ENDIF
      ENDIF
  100 CONTINUE
      DEALLOCATE(SGAS,SGAR)
*----
*  SET MACRO-GEOMETRY.
*----
      NREG=NREG2
      NBMIX=NBMIX+NMILG
      DO I=1,NREG2
         MAT(I)=IBI(I)
         VOL(I)=VOL2(I)
      ENDDO
      DEALLOCATE(RS,FRACT,VOLK,VOL2)
      DEALLOCATE(NS,IDIL,MIXGR,IBI)
      RETURN
      END
