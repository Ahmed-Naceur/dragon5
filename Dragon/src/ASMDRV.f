*DECK ASMDRV
      SUBROUTINE ASMDRV(IPSYS,IPTRK,IPMACR,IFTRAK,CDOOR,IPRNTP,NGROUP,
     >                  NBMIX,NREGIO,NANI,NANIST,NW,MATCOD,VOLUME,
     >                  LEAKSW,ITRANC,LDIFF,IBFP,TITRE,ITPIJ,LNORM,
     >                  IPHASE,ISTRM,KNORM,NALBP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Dragon assembly and pij phases.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPSYS   pointer to the pij LCM object.
* IPTRK   pointer to the tracking LCM object.
* IPMACR  pointer to the macrolib LCM object.
* IFTRAK  file unit number for tracks.
* CDOOR   name of the pij calculation door.
* IPRNTP  print option for pij calculations.
* NGROUP  number of groups treated.
* NBMIX   number of mixtures considered.
* NREGIO  number of regions considered .
* NANI    number of Legendre orders for scattering cross sections.
* NANIST  number of Legendre orders for scattering cross sections
*         if streaming leakage is present.
* NW      type of weighting for P1 cross section info (=0: P0; =1: P1).
* MATCOD  mixture code in each region.
* VOLUME  volume of each region.
* LEAKSW  leakage switch.
* ITRANC  type of transport correction.
* LDIFF   diffusion coefficient switch.
* IBFP    Fokker-Planck solution (=0: off; =1/2: on).
* TITRE   execution title.
* ITPIJ   type of collision probability available:
*         =1 scatt mod pij (wij);
*         =2 stand. pij;
*         =3 scatt mod pij+pijk (wij,wijk);
*         =4 stand. pij+pijk.
* LNORM   switch for removing leakage from collision probabilities and
*         keeping the pis information.
* IPHASE  type of assembly (=1 for ass and 2 for pij).
* ISTRM   type of streaming effect:
*         =1 no streaming effect;
*         =2 isotropic streaming effect;
*         =3 anisotropic streaming effect.
* KNORM   type of pij normalization.
* NALBP   number of physical albedos.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER   CDOOR*12,TITRE*72,TEXT12*12
      LOGICAL     LEAKSW,LDIFF,LNORM
      TYPE(C_PTR) IPSYS,IPTRK,IPMACR
      INTEGER     IFTRAK,IPRNTP,NGROUP,NBMIX,NREGIO,NANI,NANIST,NW,
     >            MATCOD(NREGIO),ITRANC,IBFP,ITPIJ,IPHASE,ISTRM,KNORM,
     >            NALBP
      REAL        VOLUME(NREGIO)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,ILCMUP=1,ILCMDN=2)
      LOGICAL     LTRANC
      CHARACTER   HSMG*130,CM*2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NPSYS
      TYPE(C_PTR) JPSYS,KPSYS,JPMACR,KPMACR
      REAL, ALLOCATABLE, DIMENSION(:) :: TEMP,ENERGY
      REAL, ALLOCATABLE, DIMENSION(:,:) :: XSSCOR,XSDIFF,ALBP
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XSSIGT,XSSIGW,ESTOPW
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NPSYS(NGROUP))
      ALLOCATE(XSSIGT(0:NBMIX,NGROUP,NW+1),XSSCOR(NBMIX,NGROUP),
     >         XSSIGW(0:NBMIX,NANIST,NGROUP),XSDIFF(0:NBMIX,NGROUP),
     >         ESTOPW(0:NBMIX,2,NGROUP),ENERGY(NGROUP+1))
      ALLOCATE(ALBP(NALBP,NGROUP))
*----
*  CHECK MIXTURE INDICES.
*----
      DO 10 I=1,NREGIO
      IF(MATCOD(I).GT.NBMIX) THEN
         WRITE (HSMG,5090) NBMIX
         CALL XABORT(HSMG)
      ENDIF
   10 CONTINUE
*----
*  RECOVER PHYSICAL ALBEDOS.
*----
      IF(NALBP.GT.0) THEN
        CALL LCMLEN(IPMACR,'ALBEDO',ILCMLN,ITYLCM)
        IF(ILCMLN.EQ.NALBP*NGROUP) THEN
          CALL LCMGET(IPMACR,'ALBEDO',ALBP)
        ELSE
          CALL LCMLIB(IPMACR)
          CALL XABORT('ASMDRV: READ ERROR ON LCM RECORD= ALBEDO')
        ENDIF
      ENDIF
*----
*  RECOVER ENERGY MESH VALUES.
*----
      IF(IBFP.GT.0) CALL LCMGET(IPMACR,'ENERGY',ENERGY)
*----
*  READ X-SECTIONS AND COMPUTE TRANSPORT CORRECTED X-SECTIONS.
*----
      IF(IPRNTP.GE.1) THEN
        IF(IPHASE.EQ.1) THEN
          WRITE(IUNOUT,6200) CDOOR
        ELSE
          WRITE(IUNOUT,6201) CDOOR
        ENDIF
        IF(ITRANC.NE.0) WRITE(IUNOUT,6101) ITRANC
      ENDIF
      IF(IPRNTP.GE.2) THEN
        WRITE(IUNOUT,6000)
        WRITE(IUNOUT,6001) (IREGIO,VOLUME(IREGIO),MATCOD(IREGIO),
     >                      IREGIO=1,NREGIO)
      ENDIF
      CALL LCMLEN(IPMACR,'GROUP',ILON,ITYLCM)
      IF(ILON.NE.NGROUP) CALL XABORT('ASMDRV: INVALID MACROLIB.')
      JPMACR=LCMGID(IPMACR,'GROUP')
      JPSYS=LCMLID(IPSYS,'GROUP',NGROUP)
      ITRAN2=0
      DO 60 IGR=1,NGROUP
        KPMACR=LCMGIL(JPMACR,IGR)
        DO 20 IW=1,MIN(NW+1,10)
        WRITE(TEXT12,'(4HNTOT,I1)') IW-1
        CALL LCMLEN(KPMACR,TEXT12,ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          CALL LCMGET(KPMACR,TEXT12,XSSIGT(1,IGR,IW))
          XSSIGT(0,IGR,IW)=0.0
        ELSE IF(IW.EQ.1) THEN
          CALL LCMLIB(KPMACR)
          CALL XABORT('ASMDRV: READ ERROR ON LCM RECORD= TOTAL')
        ELSE
          CALL LCMGET(KPMACR,'NTOT0',XSSIGT(1,IGR,IW))
          XSSIGT(0,IGR,IW)=0.0
        ENDIF
   20   CONTINUE
        DO 30 IL=1,NANIST
        WRITE(CM,'(I2.2)') IL-1
        CALL LCMLEN(KPMACR,'SIGW'//CM,ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          CALL LCMGET(KPMACR,'SIGW'//CM,XSSIGW(1,IL,IGR))
          XSSIGW(0,IL,IGR)=0.0
        ELSE IF((MOD(ITPIJ,2).EQ.1).AND.(IL.EQ.1)) THEN
          CALL LCMLIB(KPMACR)
          CALL XABORT('ASMDRV: READ ERROR ON LCM RECORD= SIGW'//CM)
        ELSE
          CALL XDRSET(XSSIGW(0,IL,IGR),NBMIX+1,0.0)
        ENDIF
*       Set XSSIGW to zero if the PIJ matrices are not to be scattering-
*       reduced.
        IF(MOD(ITPIJ,2).EQ.0) CALL XDRSET(XSSIGW(0,IL,IGR),NBMIX+1,0.0)
   30   CONTINUE
        CALL LCMLEN(KPMACR,'TRANC',ILCMLN,ITYLCM)
        LTRANC=ILCMLN.GT.0
        IF((ITRANC.NE.0).AND.LTRANC) THEN
*         TRANSPORT CORRECTION (INCLUDE THE LEAKAGE CORRECTION).
          ITRAN2=ITRANC
          CALL LCMGET(KPMACR,'TRANC',XSSCOR(1,IGR))
        ELSE
          ITRAN2=0
          CALL XDRSET(XSSCOR(1,IGR),NBMIX,0.0)
        ENDIF
        IF(ITRAN2.NE.0) THEN
*         INCLUDE TRANSPORT CORRECTION.
          DO 40 IMAT=1,NBMIX
            XSSIGT(IMAT,IGR,1)=XSSIGT(IMAT,IGR,1)-XSSCOR(IMAT,IGR)
            IF(MOD(ITPIJ,2).NE.0) XSSIGW(IMAT,1,IGR)=XSSIGW(IMAT,1,IGR)-
     >      XSSCOR(IMAT,IGR)
*           Tibere is using transport-corrected XS for the second
*           equation. Scattering reduction must be performed with
*           transport-corrected SIGS1 values.
            IF((ITPIJ.EQ.3).OR.(ISTRM.EQ.3)) THEN
              XSSIGW(IMAT,2,IGR)=XSSIGW(IMAT,2,IGR)-XSSCOR(IMAT,IGR)
            ENDIF
   40     CONTINUE
        ENDIF
*
        IF(NW.GT.0) THEN
*         PERFORM A P0_TOTAL LEAKAGE CORRECTION.
          DO 55 IW=2,MIN(NANIST,NW+1)
          DO 50 IMAT=1,NBMIX
          DELTA=XSSIGT(IMAT,IGR,1)-XSSIGT(IMAT,IGR,IW)
          IF((ITRAN2.NE.0).AND.(DELTA.NE.0.0)) THEN
            CALL XABORT('ASMDRV: CANNOT PERFORM BOTH TRANSPORT AND LEA'
     >      //'KAGE CORRECTIONS.')
          ENDIF
          XSSIGT(IMAT,IGR,IW)=XSSIGT(IMAT,IGR,1)
          XSSIGW(IMAT,IW,IGR)=XSSIGW(IMAT,IW,IGR)+DELTA
   50     CONTINUE
   55     CONTINUE
        ENDIF
*
        IF(IPRNTP.GE.3) THEN
          WRITE(IUNOUT,6002)  IGR
          WRITE(IUNOUT,6003) (IMIX,XSSIGT(IMIX,IGR,1),XSSIGW(IMIX,1,IGR)
     >    ,IMIX=1,NBMIX)
        ENDIF
        IF(LDIFF) THEN
*         INCLUDE DIFFUSION COEFFICIENTS.
          CALL LCMGET(KPMACR,'DIFF',XSDIFF(1,IGR))
          XSDIFF(0,IGR)=1.0E10
        ENDIF
        IF(IBFP.GT.0) THEN
*         INCLUDE RESTRICTED STOPPING POWER.
          ALLOCATE(TEMP(2*NBMIX))
          CALL LCMGET(KPMACR,'ESTOPW',TEMP)
          ESTOPW(0,:2,IGR)=0.0
          ESTOPW(1:NBMIX,1,IGR)=TEMP(:NBMIX)
          ESTOPW(1:NBMIX,2,IGR)=TEMP(NBMIX+1:)
          DEALLOCATE(TEMP)
        ENDIF
   60 CONTINUE
*----
*  COMPUTE THE CP OR RESPONSE MATRIX INFORMATION FOR THE SOLUTION OF
*  THE BALANCE EQUATION.
*----
      IPIJK=1
      IF(ISTRM.EQ.3) IPIJK=4
      DO 70 IGR=1,NGROUP
      NPSYS(IGR)=IGR
      KPSYS=LCMDIL(JPSYS,IGR)
      CALL LCMPUT(KPSYS,'DRAGON-TXSC',NBMIX+1,2,XSSIGT(0,IGR,1))
      IF(NW.GT.0) THEN
         CALL LCMPUT(KPSYS,'DRAGON-T1XSC',NBMIX+1,2,XSSIGT(0,IGR,2))
      ENDIF
      CALL LCMPUT(KPSYS,'DRAGON-S0XSC',(NBMIX+1)*NANI,2,XSSIGW(0,1,IGR))
      IF(LDIFF) CALL LCMPUT(KPSYS,'DRAGON-DIFF',NBMIX+1,2,XSDIFF(0,IGR))
      IF(IBFP.GT.0) THEN
        CALL LCMPUT(KPSYS,'DRAGON-ESTOP',(NBMIX+1)*2,2,ESTOPW(0,1,IGR))
        DELTAE=(ENERGY(IGR)-ENERGY(IGR+1))/1.0E6
        CALL LCMPUT(KPSYS,'DRAGON-DELTE',1,2,DELTAE)
        IF(IGR.EQ.NGROUP) THEN
          CALL LCMPUT(KPSYS,'DRAGON-ISLG',1,1,1)
        ELSE
          CALL LCMPUT(KPSYS,'DRAGON-ISLG',1,1,0)
        ENDIF
      ENDIF
      IF(NALBP.GT.0) CALL LCMPUT(KPSYS,'ALBEDO',NALBP,2,ALBP(1,IGR))
   70 CONTINUE
      IF(IBFP.GT.0) THEN
        CALL LCMPUT(IPSYS,'ECUTOFF',1,2,ENERGY(NGROUP+1)/1.0E6)
      ENDIF
      IF(IPHASE.EQ.2) THEN
         CALL DOORPV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRNTP,NGROUP,
     >   NREGIO,NBMIX,NANI,MATCOD,VOLUME,KNORM,IPIJK,LEAKSW,LNORM,
     >   TITRE,NALBP)
      ELSE
         CALL DOORAV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRNTP,NGROUP,
     >   NREGIO,NBMIX,NANI,NW,MATCOD,VOLUME,KNORM,LEAKSW,TITRE,NALBP,
     >   ISTRM)
      ENDIF
*----
*  COMPUTE THE P1 CP OR RESPONSE MATRIX INFORMATION FOR THE ECCO
*  ISOTROPIC STREAMING MODEL.
*----
      IF(ISTRM.EQ.2) THEN
         CALL LCMSIX(IPSYS,'STREAMING',ILCMUP)
         JPSYS=LCMLID(IPSYS,'GROUP',NGROUP)
         IF(ITRAN2.NE.0) THEN
*          REMOVE TRANSPORT CORRECTION.
           DO 85 IGR=1,NGROUP
           DO 80 IMAT=1,NBMIX
             XSSIGT(IMAT,IGR,1)=XSSIGT(IMAT,IGR,1)+XSSCOR(IMAT,IGR)
             XSSIGW(IMAT,1,IGR)=XSSIGW(IMAT,1,IGR)+XSSCOR(IMAT,IGR)
   80      CONTINUE
   85      CONTINUE
         ENDIF
         IF(NANIST.LE.1) CALL XABORT('ASMDRV: MISSING P1 XS INFO.')
         DO 90 IGR=1,NGROUP
         NPSYS(IGR)=IGR
         KPSYS=LCMDIL(JPSYS,IGR)
         CALL LCMPUT(KPSYS,'DRAGON-TXSC',NBMIX+1,2,XSSIGT(0,IGR,1))
         CALL LCMPUT(KPSYS,'DRAGON-S0XSC',(NBMIX+1)*(NANIST-1),2,
     >   XSSIGW(0,2,IGR))
         IF(LDIFF) CALL LCMPUT(KPSYS,'DRAGON-DIFF',NBMIX+1,2,
     >   XSDIFF(0,IGR))
   90    CONTINUE
         IPIJK=1
         IF(IPHASE.EQ.2) THEN
            CALL DOORPV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRNTP,NGROUP,
     >      NREGIO,NBMIX,NANIST-1,MATCOD,VOLUME,KNORM,IPIJK,LEAKSW,
     >      LNORM,TITRE,NALBP)
         ELSE
            CALL DOORAV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRNTP,NGROUP,
     >      NREGIO,NBMIX,NANIST-1,NW,MATCOD,VOLUME,KNORM,LEAKSW,TITRE,
     >      NALBP,ISTRM)
         ENDIF
         CALL LCMSIX(IPSYS,' ',ILCMDN)
      ENDIF
*----
*  COMPUTE RESPONSE MATRIX INFORMATION FOR THE TIBERE
*  ANISOTROPIC STREAMING MODEL.
*----
      IF((ISTRM.EQ.3).AND.(IPHASE.EQ.1)) THEN
         CALL LCMSIX(IPSYS,'STREAMING',ILCMUP)
         JPSYS=LCMLID(IPSYS,'GROUP',NGROUP)
         IF(NANIST.LE.1) CALL XABORT('ASMDRV: MISSING P1 XS INFO.')
         DO 100 IGR=1,NGROUP
         NPSYS(IGR)=IGR
         KPSYS=LCMDIL(JPSYS,IGR)
         CALL LCMPUT(KPSYS,'DRAGON-TXSC',NBMIX+1,2,XSSIGT(0,IGR,1))
         CALL LCMPUT(KPSYS,'DRAGON-S0XSC',(NBMIX+1)*(NANIST-1),2,
     >   XSSIGW(0,2,IGR))
         IF(LDIFF) CALL LCMPUT(KPSYS,'DRAGON-DIFF',NBMIX+1,2,
     >   XSDIFF(0,IGR))
  100    CONTINUE
         CALL DOORAV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRNTP,NGROUP,
     >   NREGIO,NBMIX,NANIST-1,NW,MATCOD,VOLUME,KNORM,LEAKSW,TITRE,
     >   NALBP,ISTRM)
         CALL LCMSIX(IPSYS,' ',ILCMDN)
      ENDIF
      IF(LNORM) LEAKSW=.FALSE.
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ALBP)
      DEALLOCATE(ENERGY,ESTOPW,XSDIFF,XSSCOR,XSSIGW,XSSIGT)
      DEALLOCATE(NPSYS)
      RETURN
*----
*  FORMATS
*----
 5090 FORMAT(32HASMDRV: INVALID VALUE OF NBMIX (,I5,2H).)
 6000 FORMAT(//30X,' EDITION REGION/VOLUME/MIXTURE '//
     >3(5X,'REGION',5X,'VOLUME  ',5X,'MIXTURE')/)
 6001 FORMAT(1P,3(5X,I4,4X,E12.5,4X,I4,4X))
 6002 FORMAT(//30X,' G R O U P : ',I5//31X,
     >'TOTAL AND WITHIN-GROUP MACROSCOPIC CROSS SECTIONS PER MIXTURE '/)
 6003 FORMAT(3(1X,'MIXTURE',4X,'NTOT0',11X,'SIGW',3X)/
     >1P,3(1X,I4,3X,E12.5,3X,E12.5))
 6101 FORMAT(//' USE TRANSPORT CORRECTED CROSS-SECTIONS (ITRANC=',I4,
     >' )')
 6200 FORMAT(//' COMPUTATION OF DRAGON RESPONSE MATRICES BY DOOR =',
     >3X,A12)
 6201 FORMAT(//' COMPUTATION OF DRAGON COMPLETE CP BY DOOR =',
     >3X,A12)
      END
