*DECK M2TDRV
      SUBROUTINE M2TDRV(IMPX,LOUT,IPMAC,NGRP,NBMIX,MAXMIX,NL,NBFIS,ICTR,
     1 IGMAIL,BUP,TEMP,HBM,NBM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Build an Apotrim interface file.
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
* IMPX    print index.
* LOUT    Apotrim file unit number.
* IPMAC   LCM pointer to the Macrolib.
* NGRP    number of energy groups.
* NBMIX   number of material mixtures in the Apotrim file.
* MAXMIX  number of material mixtures in the Macrolib.
* NL      maximum anisotropy level in the Apotrim file (=1 for
*         isotropic collision in LAB).
* NBFIS   maximum number of fissile isotopes in a mixture.
* ICTR    flag set to 1 if the Apotrim xs are transport corrected.
* IGMAIL  flag set to 1 to avoid writing the energy mesh on file.
* BUP     burnup of each Apotrim mixture.
* TEMP    temperature of each Apotrim mixture in Celsius.
* HBM     name of material mixtures in the Apotrim file.
* NBM     corresponding material mixtures indices in the Macrolib.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC
      INTEGER IMPX,LOUT,NGRP,NBMIX,MAXMIX,NL,NBFIS,ICTR,IGMAIL,
     1 HBM(5,NBMIX),NBM(NBMIX)
      REAL BUP(NBMIX),TEMP(NBMIX)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAC,KPMAC
      CHARACTER TEXT20*20,FMTOUT*80,CM*2
      PARAMETER(FMTOUT='(1P,6E13.5)',IOUT=6)
      INTEGER FFAGGM,LLAGGM,FFDGGM,WWGALM,FFAGM,LLAGM,NNPSNM
*----
*  ALLOCATABLE STATEMENTS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IFDG,IADR,IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR1,XTRAN,SIG,WORK,TRAN
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GAR2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IFDG(NGRP),IADR(NGRP+1),IJJ(MAXMIX),NJJ(MAXMIX),
     1 IPOS(MAXMIX))
      ALLOCATE(GAR1(NGRP),GAR2(5,NGRP),XTRAN(NGRP*NGRP),SIG(MAXMIX),
     1 WORK(NGRP*MAXMIX),TRAN(NGRP))
*----
*  RECOVER THE ENERGY MESH
*----
      IF(IGMAIL.EQ.0) THEN
        CALL LCMGET(IPMAC,'ENERGY',GAR1)
        DO 10 I=1,NGRP+1
        GAR1(I)=1.0E-6*GAR1(I)
   10   CONTINUE
        WRITE(LOUT,'(2I8)') NBMIX,NGRP
        WRITE(LOUT,FMTOUT) (GAR1(I),I=1,NGRP+1)
        IF(IMPX.GE.1) THEN
          WRITE(IOUT,4000) NBMIX,NGRP
          WRITE(IOUT,4100) (GAR1(I),I=1,NGRP+1)
        ENDIF
      ENDIF
*----
*  MIXTURE LOOP
*----
      DO 100 IMED=1,NBMIX
      WRITE(TEXT20,'(5A4)') (HBM(I0,IMED),I0=1,5)
      IF(IMPX.GT.0) WRITE(IOUT,'(/25H M2TDRV: PROCESS MIXTURE ,A20)')
     1 TEXT20
      IBM=NBM(IMED)
      JPMAC=LCMGID(IPMAC,'GROUP')
*----
*  RECOVER FISSION INFORMATION
*----
      LFIS=0
      IF(NBFIS.EQ.1) THEN
        DO 20 IGR=1,NGRP
        KPMAC=LCMGIL(JPMAC,IGR)
        CALL LCMLEN(KPMAC,'CHI',ILONG,ITYLCM)
        IF(ILONG.GT.0) THEN
          CALL LCMGET(KPMAC,'CHI',SIG)
          GAR1(IGR)=SIG(IBM)
          IF(GAR1(IGR).NE.0.0) LFIS=1
        ELSE
          GAR1(IGR)=0.0
        ENDIF
   20   CONTINUE
        IF((LFIS.EQ.1).AND.(IMPX.GE.1)) THEN
          WRITE(IOUT,1110)
          WRITE(IOUT,4100) (GAR1(IGR),IGR=1,NGRP)
        ENDIF
      ENDIF
      WRITE(LOUT,'(A20,2I5,3I3,2I10)') TEXT20,IMED,NGRP,LFIS,ICTR,NL-1,
     1 NINT(TEMP(IMED)),NINT(BUP(IMED))
      IF(LFIS.EQ.1) WRITE(LOUT,FMTOUT) (GAR1(IGR),IGR=1,NGRP)
*----
*  RECOVER TRANSPORT CORRECTION
*----
      IF(ICTR.GT.0) THEN
        DO 25 IGR=1,NGRP
        KPMAC=LCMGIL(JPMAC,IGR)
        CALL LCMLEN(KPMAC,'TRANC',ILONG1,ITYLCM)
        CALL LCMLEN(KPMAC,'SIGS01',ILONG2,ITYLCM)
        IF(ILONG1.GT.0) THEN
          CALL LCMGET(KPMAC,'TRANC',SIG)
          TRAN(IGR)=SIG(IBM)
        ELSE IF(ILONG2.GT.0) THEN
          CALL LCMGET(KPMAC,'SIGS01',SIG)
          TRAN(IGR)=SIG(IBM)
        ELSE
          TRAN(IGR)=0.0
        ENDIF
   25   CONTINUE
      ENDIF
*----
*  RECOVER REMAINING VECTOR XS INFORMATION
*----
      IF(ICTR.EQ.0) THEN
        IOF=0
        NXS=4
      ELSE
        IOF=1
        NXS=5
      ENDIF
      DO 30 IGR=1,NGRP
      KPMAC=LCMGIL(JPMAC,IGR)
      CALL LCMGET(KPMAC,'NTOT0',SIG)
      GAR2(IOF+1,IGR)=SIG(IBM)
      IF(ICTR.GT.0) THEN
        GAR2(1,IGR)=TRAN(IGR)
        GAR2(IOF+1,IGR)=GAR2(IOF+1,IGR)-TRAN(IGR)
      ENDIF
      GAR2(IOF+2,IGR)=SIG(IBM)
      CALL LCMLEN(KPMAC,'SIGS00',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
        CALL LCMGET(KPMAC,'SIGS00',SIG)
        GAR2(IOF+2,IGR)=GAR2(IOF+2,IGR)-SIG(IBM)
      ENDIF
      CALL LCMLEN(KPMAC,'N2N',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
        CALL LCMGET(KPMAC,'N2N',SIG)
        GAR2(IOF+2,IGR)=GAR2(IOF+2,IGR)+SIG(IBM)
        GAR2(IOF+4,IGR)=SIG(IBM)
      ELSE
        GAR2(IOF+4,IGR)=0.0
      ENDIF
      CALL LCMLEN(KPMAC,'N3N',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
        CALL LCMGET(KPMAC,'N3N',SIG)
        GAR2(IOF+2,IGR)=GAR2(IOF+2,IGR)+2.0*SIG(IBM)
        GAR2(IOF+4,IGR)=GAR2(IOF+4,IGR)+2.0*SIG(IBM)
      ENDIF
      CALL LCMLEN(KPMAC,'N4N',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
        CALL LCMGET(KPMAC,'N4N',SIG)
        GAR2(IOF+2,IGR)=GAR2(IOF+2,IGR)+3.0*SIG(IBM)
        GAR2(IOF+4,IGR)=GAR2(IOF+4,IGR)+3.0*SIG(IBM)
      ENDIF
      CALL LCMLEN(KPMAC,'NUSIGF',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
        CALL LCMGET(KPMAC,'NUSIGF',SIG)
        GAR2(IOF+3,IGR)=SIG(IBM)
      ELSE
        GAR2(IOF+3,IGR)=0.0
      ENDIF
   30 CONTINUE
      WRITE(IOUT,4300) NL-1
      DO 40 IGR=1,NGRP
      WRITE(LOUT,FMTOUT) (GAR2(II,IGR),II=1,NXS)
   40 CONTINUE
      IF(IMPX.GE.1) THEN
        WRITE(IOUT,1000)
        DO 50 IGR=1,NGRP
        WRITE(IOUT,'(8X,I7,1P,6E15.6)') IGR,(GAR2(II,IGR),II=1,NXS)
   50   CONTINUE
      ENDIF
*----
*  RECOVER TRANSFER XS INFORMATION
*----
      DO 90 INL=1,NL
      WRITE (CM,'(I2.2)') INL-1
      IADR(1)=1
      NNPSNM=0
      FFAGGM=NGRP+1
      LLAGGM=0
      FFDGGM=NGRP+1
      WWGALM=0
      FFAGM=1
      LLAGM=NGRP
      DO 70 IGR=1,NGRP
      KPMAC=LCMGIL(JPMAC,IGR)
      IFDG(IGR)=NGRP+1
      CALL LCMGET(KPMAC,'IJJS'//CM,IJJ)
      CALL LCMGET(KPMAC,'NJJS'//CM,NJJ)
      CALL LCMGET(KPMAC,'IPOS'//CM,IPOS)
      CALL LCMGET(KPMAC,'SCAT'//CM,WORK)
      IF(ICTR.GT.0) THEN
        IOF=IPOS(IBM)-IGR+IJJ(IBM)
        WORK(IOF)=WORK(IOF)-TRAN(IGR)
      ENDIF
      IFDG(IGR)=MIN(IFDG(IGR),IJJ(IBM)-NJJ(IBM)+1)
      IPO=IPOS(IBM)+NJJ(IBM)
      DO 60 IB=1,NJJ(IBM)
      NNPSNM=NNPSNM+1
      XTRAN(NNPSNM)=WORK(IPO-IB)*REAL(2*INL-1)
   60 CONTINUE
      IADR(IGR+1)=IADR(IGR)+(IJJ(IBM)-IFDG(IGR)+1)
   70 CONTINUE
      WRITE(LOUT,'(A20,2I5,3I3,2I10)') TEXT20,IMED,NGRP,LFIS,ICTR,INL-1,
     1 NINT(TEMP(IMED)),NINT(BUP(IMED))
      WRITE(LOUT,'(8I8)') FFAGGM,LLAGGM,FFDGGM,WWGALM,FFAGM,LLAGM
      WRITE(LOUT,'(8I8)') (IFDG(IGR),IGR=1,NGRP)
      WRITE(LOUT,'(8I8)') (IADR(IGR),IGR=1,NGRP+1)
      WRITE(LOUT,'(I10)') NNPSNM
      WRITE(LOUT,FMTOUT) (XTRAN(II),II=1,NNPSNM)
      IF(IMPX.GE.2) THEN
        WRITE(IOUT,3000) INL-1
        WRITE(IOUT,3050) FFAGGM,LLAGGM,FFDGGM,WWGALM,FFAGM,LLAGM,NNPSNM
        WRITE(IOUT,3100)
        WRITE(IOUT,4200) (IFDG(IGR),IGR=1,NGRP)
        WRITE(IOUT,3200)
        WRITE(IOUT,4200) (IADR(IGR),IGR=1,NGRP+1)
      ENDIF
*     PRINT TRANSFERT MATRICES ON LISTING, WIDLY AS THEY ARE CODED
*     IN MACROLIB FOR IMPX.EQ.2, EXPLICITLY FOR IMPX.EQ.3
      IF(IMPX.EQ.2) THEN
        WRITE(IOUT,3300)
        WRITE(IOUT,4100) (XTRAN(II),II=1,NNPSNM)
      ENDIF
      IF(IMPX.EQ.3) THEN
        WRITE(IOUT,3300)
        DO 85 IG=1,NGRP
        DO 80 IGP=1,NGRP
        SECT=0.0
        IF((IG.GE.FFAGGM).AND.(IG.LE.LLAGGM).AND.
     1    (IGP.GE.FFDGGM).AND.(IGP.LE.(FFDGGM+WWGALM-1))) THEN
             SECT=XTRAN((IG-FFAGGM)*WWGALM+IGP-FFDGGM+1)
             WRITE(IOUT,3060) IGP,IG,SECT
        ELSE IF((IGP.GE.IFDG(IG)).AND.
     1    (IGP.LE.(IADR(IG+1)-IADR(IG)+IFDG(IG)-1))
     2    .AND.(IG.GE.FFAGM).AND.(IG.LE.LLAGM)) THEN
             SECT=XTRAN(IADR(IG)+IGP-IFDG(IG))
             WRITE(IOUT,3060) IGP,IG,SECT
        ENDIF
   80   CONTINUE
   85   CONTINUE
      ENDIF
   90 CONTINUE
  100 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(TRAN,WORK,SIG,XTRAN,GAR2,GAR1)
      DEALLOCATE(IPOS,NJJ,IJJ,IADR,IFDG)
      RETURN
*
 1000 FORMAT (//29X,10(2H**)/29X,'** CROSS SECTIONS **'/29X,
     1 10(2H**)//)
 1110 FORMAT (//31X,10(2H**)/31X,'* FISSION SPECTRUM *',
     1 /31X,10(2H**))
 3000 FORMAT (//26X,15(2H**)/26X,'* P',I1,' TRANSFER CROSS SECTIONS *'/
     1 26X,15(2H**)/)
 3050 FORMAT (//10X,'FAGGM = ',I6,10X,'LAGGM = ',I6,10X,'FDGGM = ',I6
     1 /10X,'WGALM = ',I6,10X,'FAGM =  ',I6,10X,'LAGM =  ',I6
     2 /10X,'NPSNM = ',I10)
 3060 FORMAT (1X,I3,' ==>',I3,1P,E13.5)
 3100 FORMAT (//26X,6(2H**)/26X,'*   FDGM   *'/26X,6(2H**)/)
 3200 FORMAT (//26X,6(2H**)/26X,'*   IADM   *'/26X,6(2H**)/)
 3300 FORMAT (//26X,6(2H**)/26X,'*   XTRAN  *'/26X,6(2H**)/)
 4000 FORMAT (//25X,11(3H***)/25X,'*   NUMBER OF MIXTURES : ',I5,
     1 '  *'/25X,'*   ',I5,'-GROUP ENERGY MESH     *'/25X,11(3H***))
 4100 FORMAT (2X,1P,5E15.6)
 4200 FORMAT (3X,5I10)
 4300 FORMAT (//28X,13(2H**)/28X,'*  ANISOTROPY LEVEL : P',I1,' *'/
     1  28X,13(2H**))
      END
