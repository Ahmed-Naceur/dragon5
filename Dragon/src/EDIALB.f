*DECK EDIALB
      SUBROUTINE EDIALB(IPMAC2,IPFLUX,IPMACR,IPSYS,IPRINT,NBMIX,
     1 NW,B2,NGROUP,NIFISS,NGCOND,ITRANC,ILEAKS,NREGIO,MATCOD,
     2 VOLUME,KEYFLX,IGCOND,FLUXES)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute boundary current from ALBS information for use with SPH
* equivalence techniques.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMAC2  pointer to condensed macrolib information (L_MACROLIB
*         signature) built by EDI:.
* IPFLUX  pointer to the reference solution (L_FLUX signature).
* IPMACR  pointer to the reference macrolib (L_MACROLIB signature).
* IPSYS   pointer to the reference pij LCM object (L_PIJ signature).
* IPRINT  print index.
* NBMIX   number of mixtures in the reference geometry.
* NW      type of weighting for P1 cross section information:
*         = 0 P0; = 1 P1.
* B2      square buckling array.
*         For ILEAKS = 1 or 2, B2(4) is the homogeneous square buckling;
*         for ILEAKS = 3, B2(1),B2(2),B2(3) are the directional
*         heterogeneous and B2(4) is the homogeneous square buckling.
* NGROUP  number of energy groups in the reference calculation.
* NIFISS  number of fissile isotopes.
* NGCOND  number of condensed groups.
* ITRANC  type of transport correction.
* ILEAKS  type of leakage calculation: =0: no leakage; =1: homogeneous
*         leakage (Diffon); =2: isotropic streaming (Ecco);
*         =3: anisotropic streaming (Tibere).
* NREGIO  number of regions in the reference geometry.
* MATCOD  mixture index in region.
* VOLUME  volume of region.
* KEYFLX  position of average fluxes.
* IGCOND  limit of condensed groups.
* FLUXES  fluxes.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC2,IPFLUX,IPMACR,IPSYS
      INTEGER IPRINT,NBMIX,NW,NGROUP,NIFISS,NGCOND,ITRANC,ILEAKS,
     1 NREGIO,MATCOD(NREGIO),KEYFLX(NREGIO),IGCOND(NGCOND)
      REAL B2(4),VOLUME(NREGIO),FLUXES(NREGIO,NGROUP,NW+1)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPSYS,KPSYS,IPSYS2,JPMACR,KPMACR,JPFLUX
      CHARACTER TEXT5*5,SUFF(2)*2
      DOUBLE PRECISION SUM,SUD
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NJJ,IJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: COURIN,SIGMA,XSCAT,DIFF,GAMMA,
     1 SIG1,WORKD
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WORK
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: PRODUC
      DATA SUFF/'00','01'/
*----
*  SCRATCH STORAGE ALLOCATION
*   COURIN  ingoing currents (4*J-/S). PIS information must be available
*           on LCM.
*----
      ALLOCATE(NJJ(NBMIX),IJJ(NBMIX),IPOS(NBMIX))
      ALLOCATE(COURIN(NGCOND),PRODUC(NREGIO,NGCOND,NIFISS),
     1 WORK(NREGIO,2),SIGMA(0:NBMIX*NIFISS),XSCAT(NBMIX*NGROUP),
     2 DIFF(NGROUP),GAMMA(NGROUP))
*----
*  CONSISTENCY TESTS
*----
      IF(.NOT.C_ASSOCIATED(IPSYS)) THEN
         CALL XABORT('EDIALB: THE L_PIJ INFO IS NOT AVAILABLE.')
      ENDIF
      JPSYS=LCMGID(IPSYS,'GROUP')
      KPSYS=LCMGIL(JPSYS,1)
      CALL LCMLEN(KPSYS,'DRAGON-WIS',IXSLEN,ITYLCM)
      IF(IXSLEN.NE.NREGIO) THEN
        CALL LCMLIB(KPSYS)
        WRITE(TEXT5,'(I5)') NREGIO
        CALL XABORT('EDIALB: THE ALBS OPTION OF THE ASM: MODULE HAS NO'
     >  //'T BEEN ACTIVATED. NREGIO='//TEXT5)
      ENDIF
*----
*  COMPUTE THE FISSION RATE INFORMATION
*----
      SIGMA(0)=0.0
      CALL XDRSET(PRODUC,NIFISS*NGCOND*NREGIO,0.0)
      IGRFIN=0
      JPMACR=LCMGID(IPMACR,'GROUP')
      DO 45 IGRCD=1,NGCOND
      IGRDEB=IGRFIN+1
      IGRFIN=IGCOND(IGRCD)
      DO 40 IGR=IGRDEB,IGRFIN
      KPMACR=LCMGIL(JPMACR,IGR)
      CALL LCMLEN(KPMACR,'NUSIGF',ILCMLN,ITYLCM)
      IF(ILCMLN.GT.0) THEN
         CALL LCMGET(KPMACR,'NUSIGF',SIGMA(1))
         DO 35 IFIS=1,NIFISS
         DO 30 IREG=1,NREGIO
         IBM=MATCOD(IREG)
         IF(IBM.GT.0) THEN
            SS=FLUXES(IREG,IGR,1)*SIGMA((IFIS-1)*NBMIX+IBM)
            PRODUC(IREG,IGRCD,IFIS)=PRODUC(IREG,IGRCD,IFIS)+SS
         ENDIF
  30     CONTINUE
  35     CONTINUE
      ENDIF
  40  CONTINUE
  45  CONTINUE
      CALL LCMLEN(IPFLUX,'K-EFFECTIVE',ILCMLN,ITYLCM)
      IF(ILCMLN.EQ.1) THEN
         CALL LCMGET(IPFLUX,'K-EFFECTIVE',EIGENK)
      ELSE
         EIGENK=1.0
      ENDIF
      IF(IPRINT.GT.5) WRITE(6,'(/16H EDIALB: EIGENK=,1P,E12.4)') EIGENK
*----
*  COMPUTE MERGED/CONDENSED CROSS SECTIONS
*----
      IF(ILEAKS.EQ.1) CALL LCMGET(IPFLUX,'DIFFB1HOM',DIFF)
      IF(NW.EQ.1) CALL LCMGET(IPFLUX,'GAMMA',GAMMA)
      CALL LCMSIX(IPMAC2,'ADF',1)
      DO 180 INL=1,NW+1
      IGRFIN=0
      DO 175 IGRCD=1,NGCOND
      COURIN(IGRCD)=0.0
      IGRDEB=IGRFIN+1
      IGRFIN=IGCOND(IGRCD)
      IF((ILEAKS.EQ.2).OR.(ILEAKS.EQ.3)) THEN
         CALL LCMLEN(IPFLUX,'FLUX',ILON,ITYLCM)
         IF(ILON.EQ.0) CALL XABORT('EDIALB: MISSING FLUX INFO.')
         JPFLUX=LCMGID(IPFLUX,'FLUX')
      ENDIF
      DO 170 IGR=IGRDEB,IGRFIN
      KPMACR=LCMGIL(JPMACR,IGR)
      CALL LCMLEN(KPMACR,'NTOT0',ILCMLN,ITYLCM)
      IF(ILCMLN.GT.0) THEN
         CALL LCMGET(KPMACR,'NTOT0',SIGMA(1))
      ELSE
         CALL XABORT('EDIALB: READ ERROR ON LCM RECORD= TOTAL.')
      ENDIF
      IF((ITRANC.NE.0).AND.(INL.EQ.1)) THEN
*        TRANSPORT CORRECTION.
         ALLOCATE(SIG1(NBMIX))
         CALL LCMGET(KPMACR,'TRANC',SIG1)
         DO 50 IMAT=1,NBMIX
         SIGMA(IMAT)=SIGMA(IMAT)-SIG1(IMAT)
   50    CONTINUE
         DEALLOCATE(SIG1)
      ENDIF
      IF((ILEAKS.EQ.2).OR.(ILEAKS.EQ.3)) THEN
         CALL LCMLEL(JPFLUX,IGR,ILCMLN,ITYLCM)
         IF(ILCMLN.EQ.0) CALL XABORT('EDIALB: MISSING FLUX INFO.')
         ALLOCATE(WORKD(ILCMLN))
         CALL LCMGDL(JPFLUX,IGR,WORKD)
      ENDIF
      ZNUM=0.0
      IF((NW.EQ.1).AND.(INL.EQ.2)) THEN
*        USE WITH THE FUNDAMENTAL CURRENT EQUATION OF THE ECCO MODEL.
         ZDEN=0.0
         DO 55 IREG=1,NREGIO
         IBM=MATCOD(IREG)
         ZNUM=ZNUM+SIGMA(IBM)*FLUXES(IREG,IGR,1)*VOLUME(IREG)
         ZDEN=ZDEN+FLUXES(IREG,IGR,1)*VOLUME(IREG)
   55    CONTINUE
         ZNUM=ZNUM/ZDEN
      ENDIF
      DO 60 IREG=1,NREGIO
      IBM=MATCOD(IREG)
      IF((NW.EQ.1).AND.(INL.EQ.1)) THEN
        ZLEAK=B2(4)*FLUXES(IREG,IGR,2)
      ELSE IF((NW.EQ.1).AND.(INL.EQ.2)) THEN
        ZLEAK=(-(1.0-GAMMA(IGR))*(ZNUM-SIGMA(IBM))*FLUXES(IREG,IGR,2)
     >        -FLUXES(IREG,IGR,1)/3.0)/GAMMA(IGR)
      ELSE IF(ILEAKS.EQ.1) THEN
        ZLEAK=DIFF(IGR)*B2(4)*FLUXES(IREG,IGR,1)
      ELSE IF(ILEAKS.EQ.2) THEN
           ZLEAK=B2(4)*WORKD(KEYFLX(IREG)+ILCMLN/2)
      ELSE IF(ILEAKS.EQ.3) THEN
        ZLEAK=B2(1)*WORKD(KEYFLX(IREG)+ILCMLN/4)+
     >        B2(2)*WORKD(KEYFLX(IREG)+ILCMLN/2)+
     >        B2(3)*WORKD(KEYFLX(IREG)+3*ILCMLN/4)
      ELSE
        ZLEAK=0.0
      ENDIF
      WORK(IREG,1)=-ZLEAK
   60 CONTINUE
      IF((ILEAKS.EQ.2).OR.(ILEAKS.EQ.3)) DEALLOCATE(WORKD)
*
      CALL LCMLEN(KPMACR,'CHI',ILCMLN,ITYLCM)
      IF((ILCMLN.GT.0).AND.(INL.EQ.1)) THEN
         DO 85 IFIS=1,NIFISS
         CALL LCMGET(KPMACR,'CHI',SIGMA(1))
         DO 80 IREG=1,NREGIO
         IBM=MATCOD(IREG)
         IF(IBM.GT.0) THEN
            DO 70 JGRCD=1,NGCOND
            SS=SIGMA((IFIS-1)*NBMIX+IBM)*PRODUC(IREG,JGRCD,IFIS)/EIGENK
            WORK(IREG,1)=WORK(IREG,1)+SS
   70       CONTINUE
         ENDIF
   80    CONTINUE
   85    CONTINUE
      ENDIF
*
      CALL LCMLEN(KPMACR,'SIGW'//SUFF(INL),ILCMLN,ITYLCM)
      IF(ILCMLN.GT.0) THEN
         CALL LCMGET(KPMACR,'SIGW'//SUFF(INL),SIGMA(1))
      ELSE
         CALL XDRSET(SIGMA,NBMIX,0.0)
      ENDIF
      IF((ITRANC.NE.0).AND.(INL.EQ.1)) THEN
*        TRANSPORT CORRECTION.
         ALLOCATE(SIG1(NBMIX))
         CALL LCMGET(KPMACR,'TRANC',SIG1)
         DO 120 IMAT=1,NBMIX
         SIGMA(IMAT)=SIGMA(IMAT)-SIG1(IMAT)
  120    CONTINUE
         DEALLOCATE(SIG1)
      ENDIF
      CALL LCMLEN(KPMACR,'NJJS'//SUFF(INL),ILCMLN,ITYLCM)
      IF(ILCMLN.GT.0) THEN
         CALL LCMGET(KPMACR,'NJJS'//SUFF(INL),NJJ)
         CALL LCMGET(KPMACR,'IJJS'//SUFF(INL),IJJ)
         CALL LCMGET(KPMACR,'IPOS'//SUFF(INL),IPOS)
         CALL LCMGET(KPMACR,'SCAT'//SUFF(INL),XSCAT)
         DO 150 IREG=1,NREGIO
         IBM=MATCOD(IREG)
         IF(IBM.GT.0) THEN
            JGRFIN=0
            DO 140 JGRCD=1,NGCOND
            SS=0.0D0
            JGRDEB=JGRFIN+1
            JGRFIN=IGCOND(JGRCD)
            J2=MIN(JGRFIN,IJJ(IBM))
            J1=MAX(JGRDEB,IJJ(IBM)-NJJ(IBM)+1)
            IPO=IPOS(IBM)+IJJ(IBM)-J2
            DO 130 JGR=J2,J1,-1
            IF(IGR.EQ.JGR) THEN
               SS=SS+SIGMA(IBM)*FLUXES(IREG,JGR,INL)
            ELSE
               SS=SS+XSCAT(IPO)*FLUXES(IREG,JGR,INL)
            ENDIF
            IPO=IPO+1
  130       CONTINUE
            IF(INL.EQ.2) SS=SS/GAMMA(IGR)
            WORK(IREG,1)=WORK(IREG,1)+SS
  140       CONTINUE
         ENDIF
  150    CONTINUE
      ENDIF
*----
*  COMPUTE BOUNDARY CURRENTS
*----
      IF(INL.EQ.1) THEN
         JPSYS=LCMGID(IPSYS,'GROUP')
      ELSE IF(INL.EQ.2) THEN
         IPSYS2=LCMGID(IPSYS,'STREAMING')
         JPSYS=LCMGID(IPSYS2,'GROUP')
      ENDIF
      KPSYS=LCMGIL(JPSYS,IGR)
      CALL LCMLEN(KPSYS,'DRAGON-WIS',ILCMLN,ITYLCM)
      IF(ILCMLN.EQ.NREGIO) THEN
         CALL LCMGET(KPSYS,'DRAGON-WIS',WORK(1,2))
         CALL LCMGET(KPSYS,'DRAGON-TXSC',SIGMA(0))
         SUM=0.0D0
         SUD=0.0D0
         DO 160 IREG=1,NREGIO
         FACTOR=VOLUME(IREG)
         IBM=MATCOD(IREG)
         SUM=SUM+SIGMA(IBM)*FACTOR*FLUXES(IREG,IGR,1)
     >          -WORK(IREG,1)*FACTOR*(1.0-WORK(IREG,2))
         SUD=SUD+SIGMA(IBM)*FACTOR*WORK(IREG,2)
  160    CONTINUE
         COURIN(IGRCD)=COURIN(IGRCD)+REAL(SUM/SUD)
      ENDIF
  170 CONTINUE
  175 CONTINUE
      CALL LCMPUT(IPMAC2,'ALBS'//SUFF(INL),NGCOND,2,COURIN)
      IF(IPRINT.GT.3) THEN
         WRITE(6,900) SUFF(INL),(COURIN(IGR),IGR=1,NGCOND)
         WRITE(6,'(/)')
      ENDIF
  180 CONTINUE
      CALL LCMSIX(IPMAC2,' ',2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAMMA,DIFF,XSCAT,SIGMA,WORK,PRODUC,COURIN)
      DEALLOCATE(IPOS,IJJ,NJJ)
      RETURN
*
  900 FORMAT(/10H EDIALB: P,A2,36H IN-CURRENTS (4J-/S) PER MACRO-GROUP,
     > 5HS ARE/(1X,1P,10E13.5))
      END
