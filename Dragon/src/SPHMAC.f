*DECK SPHMAC
      SUBROUTINE SPHMAC(IPMACR,IPRINT,NMERGE,NALBP,NGCOND,ISCAT,NW,
     1 NIFISS,ILEAKS,VOLMER,FLXMER,SUNMER,SIGT,SIGW,DIFF,ZLEAK,OUTG,
     2 ALB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recovery of the reference merged/condensed set of cross sections
* to be used by an SPH homogenization algorithm.
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
* IPMACR  pointer to the condensed macrolib (L_MACROLIB signature).
* IPRINT  print flag (equal to 0 for no print).
* NMERGE  number of merged regions.
* NALBP   number of physical albedos.
* NGCOND  number of condensed groups.
* ISCAT   scattering anisotropy in the reference set of cross sections
*         (=1 isotropic in LAB; =2 linearly-anisotropic in LAB).
* NW      type of weighting for PN cross section info (=0 P0; =1 P1).
* NIFISS  number of fissile isotopes.
* ILEAKS  type of leakage calculation: =0 no leakage; =1 homogeneous
*         leakage (Diffon); =2 isotropic streaming (Ecco);
*         =3 anisotropic streaming (Tibere).
*
*Parameters: output
* VOLMER  merged volumes.
* FLXMER  merged/condensed averaged fluxes.
* SUNMER  merged/condensed production (fission + scattering) cross
*         sections. The third dimension is for secondary neutrons.
* SIGT    merged/condensed total P0 and P1 cross sections.
* SIGW    merged/condensed within-group scattering cross sections.
* DIFF    merged/condensed diffusion coefficients.
* ZLEAK   merged/condensed DB2 leakage rates.
* OUTG    merged/condensed leakage rates.
* ALB     physical albedos.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMACR
      INTEGER IPRINT,NMERGE,NALBP,NGCOND,ISCAT,NW,NIFISS,ILEAKS
      REAL VOLMER(NMERGE),FLXMER(NMERGE,NGCOND),
     1 SUNMER(NMERGE,NGCOND,NGCOND,ISCAT),SIGT(NMERGE,NGCOND,NW+1),
     2 SIGW(NMERGE,NGCOND,ISCAT+1),DIFF(NMERGE,NGCOND),
     3 ZLEAK(NMERGE,NGCOND),OUTG(NGCOND),ALB(NALBP,NGCOND)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE)
      DOUBLE PRECISION ZNUMER,ZLEAKA,ZDENUM
      CHARACTER HSIGN*12,SUFF*2,TEXT12*12
      TYPE(C_PTR) JPMACR,KPMACR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NJJ,IJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGMA,XSCAT,DIFHOM,SIG1
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: PRODUC
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NJJ(NMERGE),IJJ(NMERGE),IPOS(NMERGE))
      ALLOCATE(PRODUC(NMERGE,NGCOND,NIFISS),SIGMA(NMERGE*NIFISS),
     1 XSCAT(NMERGE*NGCOND),DIFHOM(NGCOND))
*----
*  RECOVER MACROLIB INFORMATION
*----
      CALL LCMGTC(IPMACR,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_MACROLIB') CALL XABORT('SPHMAC: MACROLIB EXPECTED')
      CALL LCMGET(IPMACR,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.NGCOND) CALL XABORT('SPHMAC: INVALID NGCOND')
      IF(ISTATE(2).NE.NMERGE) CALL XABORT('SPHMAC: INVALID NMERGE')
      NL=ISTATE(3)
      IF(ISTATE(4).NE.NIFISS) CALL XABORT('SPHMAC: INVALID NIFISS')
      ITRANC=ISTATE(6)
      ILEAK=ISTATE(9)
      IF(MAX(1,ISTATE(10)).NE.NW) CALL XABORT('SPHMAC: INVALID NW')
*----
*  SET OUTPUT INFORMATION TO ZERO
*----
      CALL XDRSET(PRODUC,NMERGE*NGCOND*NIFISS,0.0)
      CALL XDRSET(FLXMER,NMERGE*NGCOND,0.0)
      CALL XDRSET(SIGT,NMERGE*NGCOND*(NW+1),0.0)
      CALL XDRSET(SIGW,NMERGE*NGCOND*(ISCAT+1),0.0)
      CALL XDRSET(SUNMER,NMERGE*NGCOND*NGCOND*ISCAT,0.0)
      CALL XDRSET(ZLEAK,NMERGE*NGCOND,0.0)
*----
*  RECOVER FLUX AND COMPUTE THE FISSION RATE INFORMATION
*----
      ZNUMER=0.0D0
      ZLEAKA=0.0D0
      ZDENUM=0.0D0
      CALL LCMGET(IPMACR,'VOLUME',VOLMER)
      IF(NALBP.GT.0) CALL LCMGET(IPMACR,'ALBEDO',ALB)
      JPMACR=LCMGID(IPMACR,'GROUP')
      DO 40 IGR=1,NGCOND
      KPMACR=LCMGIL(JPMACR,IGR)
      CALL LCMLEN(KPMACR,'FLUX-INTG',ILCMLN,ITYLCM)
      IF(ILCMLN.EQ.0) CALL XABORT('SPHMAC: MISSING FLUX-INTG INFO')
      CALL LCMGET(KPMACR,'FLUX-INTG',FLXMER(1,IGR))
      CALL LCMLEN(KPMACR,'NUSIGF',ILCMLN,ITYLCM)
      IF(ILCMLN.GT.0) THEN
         CALL LCMGET(KPMACR,'NUSIGF',SIGMA(1))
         DO 35 IFIS=1,NIFISS
         DO 30 IBM=1,NMERGE
         SS=FLXMER(IBM,IGR)*SIGMA((IFIS-1)*NMERGE+IBM)
         PRODUC(IBM,IGR,IFIS)=PRODUC(IBM,IGR,IFIS)+SS
         ZNUMER=ZNUMER+SS
   30    CONTINUE
   35    CONTINUE
      ENDIF
   40 CONTINUE
*----
*  RECOVER EIGENVALUES
*----
      CALL LCMLEN(IPMACR,'K-EFFECTIVE',ILCMLN,ITYLCM)
      IF(ILCMLN.EQ.1) THEN
         CALL LCMGET(IPMACR,'K-EFFECTIVE',EIGENK)
      ELSE
         EIGENK=1.0
      ENDIF
      IF(IPRINT.GT.5) WRITE(6,'(/16H SPHMAC: EIGENK=,1P,E12.4)') EIGENK
      CALL LCMLEN(IPMACR,'B2  B1HOM',ILCMLN,ITYLCM)
      IF(ILCMLN.EQ.1) THEN
         CALL LCMGET(IPMACR,'B2  B1HOM',B2)
      ELSE
         B2=0.0
      ENDIF
      IF(IPRINT.GT.5) WRITE(6,'(/12H SPHMAC: B2=,1P,E12.4)') B2
*----
*  RECOVER MERGED/CONDENSED CROSS SECTIONS
*----
      DO 175 IGR=1,NGCOND
      KPMACR=LCMGIL(JPMACR,IGR)
      CALL LCMLEN(KPMACR,'NTOT0',ILCMLN,ITYLCM)
      IF(ILCMLN.GT.0) THEN
         CALL LCMGET(KPMACR,'NTOT0',SIGMA(1))
      ELSE
         CALL XABORT('SPHMAC: MISSING NTOT0 INFO')
      ENDIF
      DO 45 IBM=1,NMERGE
      ZDENUM=ZDENUM+SIGMA(IBM)*FLXMER(IBM,IGR)
   45 CONTINUE
      IF(ITRANC.NE.0) THEN
*        TRANSPORT CORRECTION.
         ALLOCATE(SIG1(NMERGE))
         CALL LCMGET(KPMACR,'TRANC',SIG1)
         DO 50 IBM=1,NMERGE
         SIGMA(IBM)=SIGMA(IBM)-SIG1(IBM)
   50    CONTINUE
         DEALLOCATE(SIG1)
      ENDIF
      IF(ILEAKS.EQ.1) THEN
         CALL LCMLEN(KPMACR,'DIFF',ILCMLN,ITYLCM)
         IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPMACR,'DIFF',DIFF(1,IGR))
         ELSE
            CALL LCMGET(IPMACR,'DIFHOMB1HOM',DIFHOM)
            DO 55 I=1,NMERGE
            DIFF(I,IGR)=DIFHOM(IGR)
   55       CONTINUE
         ENDIF
      ENDIF
      DO 60 IBM=1,NMERGE
      IF(ILEAKS.EQ.1) THEN
         ZLEAK(IBM,IGR)=DIFF(IBM,IGR)*B2*FLXMER(IBM,IGR)
      ELSE IF(ILEAKS.GT.1) THEN
         CALL XABORT('SPHMAC: LEAKAGE MODEL NOT IMPLEMENTED')
      ELSE
         ZLEAK(IBM,IGR)=0.0
      ENDIF
      SIGT(IBM,IGR,1)=SIGMA(IBM)
      ZLEAKA=ZLEAKA+ZLEAK(IBM,IGR)
   60 CONTINUE
*
      CALL LCMLEN(KPMACR,'CHI',ILCMLN,ITYLCM)
      IF(ILCMLN.GT.0) THEN
         DO 75 IFIS=1,NIFISS
         CALL LCMGET(KPMACR,'CHI',SIGMA(1))
         DO 70 IBM=1,NMERGE
         DO 65 JGR=1,NGCOND
         SS=SIGMA((IFIS-1)*NMERGE+IBM)*PRODUC(IBM,JGR,IFIS)/EIGENK
         SUNMER(IBM,JGR,IGR,1)=SUNMER(IBM,JGR,IGR,1)+SS
   65    CONTINUE
   70    CONTINUE
   75    CONTINUE
      ENDIF
      DO 90 IW=2,MIN(NW+1,10)
      WRITE(TEXT12,'(4HNTOT,I1)') IW-1
      CALL LCMLEN(KPMACR,TEXT12,ILCMLN,ITYLCM)
      IF(ILCMLN.GT.0) THEN
         CALL LCMGET(KPMACR,TEXT12,SIGMA(1))
         DO 80 IBM=1,NMERGE
         SIGT(IBM,IGR,IW)=SIGMA(IBM)
   80    CONTINUE
      ELSE
         DO 85 IBM=1,NMERGE
         SIGT(IBM,IGR,IW)=SIGT(IBM,IGR,1)
   85    CONTINUE
      ENDIF
   90 CONTINUE
*----
*  PROCESS SCATTERING INFORMATION
*----
      DO 170 INL=1,ISCAT
      WRITE(SUFF,'(I2.2)') INL-1
      CALL LCMLEN(KPMACR,'SIGW'//SUFF,ILCMLN,ITYLCM)
      IF(ILCMLN.GT.0) THEN
         CALL LCMGET(KPMACR,'SIGW'//SUFF,SIGMA(1))
      ELSE
         CALL XDRSET(SIGMA,NMERGE,0.0)
      ENDIF
      IF((ITRANC.NE.0).AND.(INL.EQ.1)) THEN
*        TRANSPORT CORRECTION.
         ALLOCATE(SIG1(NMERGE))
         CALL LCMGET(KPMACR,'TRANC',SIG1)
         DO 120 IBM=1,NMERGE
         SIGMA(IBM)=SIGMA(IBM)-SIG1(IBM)
  120    CONTINUE
         DEALLOCATE(SIG1)
      ENDIF
      CALL LCMLEN(KPMACR,'NJJS'//SUFF,ILCMLN,ITYLCM)
      IF(ILCMLN.GT.0) THEN
         CALL LCMGET(KPMACR,'NJJS'//SUFF,NJJ)
         CALL LCMGET(KPMACR,'IJJS'//SUFF,IJJ)
         CALL LCMGET(KPMACR,'IPOS'//SUFF,IPOS)
         CALL LCMGET(KPMACR,'SCAT'//SUFF,XSCAT)
         DO 150 IBM=1,NMERGE
         IPO=IPOS(IBM)
         DO 130 JGR=IJJ(IBM),IJJ(IBM)-NJJ(IBM)+1,-1
         IF(IGR.EQ.JGR) THEN
            SS=SIGMA(IBM)*FLXMER(IBM,JGR)
         ELSE
            SS=XSCAT(IPO)*FLXMER(IBM,JGR)
         ENDIF
         IF(IGR.EQ.JGR) SIGW(IBM,IGR,INL)=SIGW(IBM,IGR,INL)+SS
         SUNMER(IBM,JGR,IGR,INL)=SUNMER(IBM,JGR,IGR,INL)+SS
         IF(INL.EQ.1) ZDENUM=ZDENUM-SS
         IPO=IPO+1
  130    CONTINUE
  150    CONTINUE
      ENDIF
  170 CONTINUE
  175 CONTINUE
*----
*  COMPUTE REFERENCE LEAKAGE RATES
*----
      DO 182 IGR=1,NGCOND
      OUTG(IGR)=0.0
      DO 181 IBM=1,NMERGE
      OUTG(IGR)=OUTG(IGR)-SIGT(IBM,IGR,1)*FLXMER(IBM,IGR)-
     1 ZLEAK(IBM,IGR)
      DO 180 JGR=1,NGCOND
      OUTG(IGR)=OUTG(IGR)+SUNMER(IBM,JGR,IGR,1)
  180 CONTINUE
  181 CONTINUE
  182 CONTINUE
*
      DO 202 INL=1,ISCAT
      DO 201 IGR=1,NGCOND
      DO 200 IBM=1,NMERGE
      IF(VOLMER(IBM).NE.0.0) THEN
         SIGW(IBM,IGR,INL)=SIGW(IBM,IGR,INL)/FLXMER(IBM,IGR)
         DO 190 JGR=1,NGCOND
         SUNMER(IBM,JGR,IGR,INL)=SUNMER(IBM,JGR,IGR,INL)/FLXMER(IBM,JGR)
  190    CONTINUE
      ENDIF
  200 CONTINUE
  201 CONTINUE
  202 CONTINUE
      DO 215 IGR=1,NGCOND
      DO 210 IBM=1,NMERGE
      IF(VOLMER(IBM).NE.0.0) THEN
         ZLEAK(IBM,IGR)=ZLEAK(IBM,IGR)/FLXMER(IBM,IGR)
         FLXMER(IBM,IGR)=FLXMER(IBM,IGR)/VOLMER(IBM)
      ENDIF
  210 CONTINUE
  215 CONTINUE
*----
*  PRINT INFORMATION
*----
      IF(IPRINT.GT.4) THEN
         WRITE(6,'(/33H SPHMAC: type of PN weighting NW=,I2)') NW
         WRITE(6,240) ZNUMER/ZDENUM,ZNUMER/(ZDENUM+ZLEAKA)
         WRITE(6,250) 'VOLMER',(VOLMER(IKK),IKK=1,NMERGE)
         DO 220 IW=1,NW+1
         WRITE(TEXT12,'(4HNTOT,I1)') IW-1
         WRITE(6,250) TEXT12,((SIGT(IKK,IGR,IW),IKK=1,NMERGE),
     >   IGR=1,NGCOND)
  220    CONTINUE
         WRITE(6,250) 'FLXMER',((FLXMER(IKK,IGR),IKK=1,NMERGE),
     >   IGR=1,NGCOND)
         WRITE(6,250) 'ZLEAK',((ZLEAK(IKK,IGR),IKK=1,NMERGE),
     >   IGR=1,NGCOND)
         IF(NALBP.GT.0) THEN
           WRITE(6,250) 'ALBEDO',((ALB(IAL,IGR),IAL=1,NALBP),
     >     IGR=1,NGCOND)
           WRITE(6,250) 'OUTG',(OUTG(IGR),IGR=1,NGCOND)
         ENDIF
         DO 230 INL=1,ISCAT
         WRITE(SUFF,'(I2.2)') INL-1
         WRITE(6,250) 'SIGW'//SUFF,((SIGW(IKK,IGR,INL),IKK=1,NMERGE)
     >   ,IGR=1,NGCOND)
         WRITE(6,250) 'SUNMER'//SUFF,(((SUNMER(IKK,IGR,JGR,INL),
     >   IKK=1,NMERGE),IGR=1,NGCOND),JGR=1,NGCOND)
  230    CONTINUE
         IF(ILEAKS.EQ.1) THEN
            WRITE(6,250) 'DIFF',((DIFF(IKK,IGR),IKK=1,NMERGE),IGR=1,
     >      NGCOND)
         ENDIF
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DIFHOM,XSCAT,SIGMA,PRODUC)
      DEALLOCATE(IPOS,IJJ,NJJ)
      RETURN
*
  240 FORMAT(/20H SPHMAC: K-INFINITY=,1P,D13.6/8X,12HK-EFFECTIVE=,D13.6,
     > 25H (FUNDAMENTAL MODE VALUE))
  250 FORMAT(/26H SPHMAC: VALUES OF VECTOR ,A,4H ARE/(1X,1P,10E13.5))
      END
