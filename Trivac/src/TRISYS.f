*DECK TRISYS
      SUBROUTINE TRISYS(IPTRK,IPMACR,IPMACP,IPSYS,IMPX,NGRP,NEL,NBFIS,
     1 NALBP,IPR,MAT,VOL,NBMIX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the diffusion coefficient and cross-section data in the LCM
* object with pointer IPMACR, compute and store the corresponding 
* Trivac system matrices (or a perturbation to the system matrices).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPTRK   L_TRACK pointer to the TRIVAC tracking information.
* IPMACR  L_MACROLIB pointer to the unperturbed cross sections.
* IPMACP  L_MACROLIB pointer to the perturbed cross sections if
*         IPR.gt.0. Equal to IPMACR if IPR=0.
* IPSYS   L_SYSTEM pointer to system matrices.
* IMPX    print parameter (equal to zero for no print).
* NGRP    number of energy groups.
* NEL     total number of finite elements.
* NBFIS   number of fissionable isotopes.
* NALBP   number of physical albedos per energy group.
* IPR     type of assembly:
*         =0: calculation of the system matrices;
*         =1: calculation of the derivative of these matrices;
*         =2: calculation of the first variation of these matrices;
*         =3: identical to IPR=2, but these variation are added to
*         unperturbed system matrices.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* NBMIX   total number of material mixtures in the macrolib.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPMACR,IPMACP,IPSYS
      INTEGER IMPX,NGRP,NEL,NBFIS,NALBP,IPR,MAT(NEL),NBMIX
      REAL VOL(NEL)
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXDIG*12,HSMG*131
      LOGICAL LFIS
      TYPE(C_PTR) JPMACR,KPMACR,JPMACP,KPMACP
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJJ,NJJ,IPOS
      REAL, DIMENSION(:), ALLOCATABLE :: WORK
      REAL, DIMENSION(:,:), ALLOCATABLE :: GAMMA,SGD,DSGD,ZUFIS
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: CHI
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(NBMIX),NJJ(NBMIX),IPOS(NBMIX))
      ALLOCATE(GAMMA(NALBP,NGRP),SGD(NBMIX,4),DSGD(NBMIX,4),
     1 WORK(NBMIX*NGRP),CHI(NBMIX,NBFIS,NGRP),ZUFIS(NBMIX,NBFIS))
*----
*  PROCESS PHYSICAL ALBEDOS.
*----
      IF(NALBP.GT.0) THEN
         CALL TRIALB(IPTRK,IPMACR,IPMACP,IPSYS,NGRP,NALBP,IPR,GAMMA)
      ENDIF
*----
*  LOOP OVER ENERGY GROUPS
*----
      JPMACR=LCMGID(IPMACR,'GROUP')
      JPMACP=LCMGID(IPMACP,'GROUP')
      DO 110 IGR=1,NGRP
*     PROCESS SECONDARY GROUP IGR.
      KPMACR=LCMGIL(JPMACR,IGR)
      KPMACP=LCMGIL(JPMACP,IGR)
*----
*  PROCESS LEAKAGE AND REMOVAL TERMS
*----
      CALL LCMLEN(KPMACR,'NTOT0',LENGT,ITYLCM)
      IF(LENGT.EQ.0) THEN
         CALL XABORT('TRISYS: NO TOTAL CROSS SECTIONS.')
      ELSE IF(LENGT.GT.NBMIX) THEN
         CALL XABORT('TRISYS: INVALID LENGTH FOR TOTAL CROSS SECTIONS.')
      ENDIF
      CALL LCMGET(KPMACR,'NTOT0',SGD(1,4))
      CALL LCMLEN(KPMACR,'SIGW00',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
         IF(LENGT.GT.NBMIX) CALL XABORT('TRISYS: INVALID LENGTH FOR '
     1   //'''SIGW00'' CROSS SECTIONS.')
         CALL LCMGET(KPMACR,'SIGW00',SGD(1,1))
         DO 10 IBM=1,LENGT
         SGD(IBM,4)=SGD(IBM,4)-SGD(IBM,1)
   10    CONTINUE
      ENDIF
      CALL LCMLEN(KPMACR,'DIFF',LENGT1,ITYLCM)
      IF(LENGT1.GT.0) THEN
         IF(LENGT1.GT.NBMIX) CALL XABORT('TRISYS: INVALID LENGTH FOR'
     1   //' DIFF (ISOTROPIC DIFFUSION COEFFICIENT).')
         CALL LCMGET(KPMACR,'DIFF',SGD(1,1))
         DO 20 IBM=1,LENGT1
         SGD(IBM,2)=SGD(IBM,1)
         SGD(IBM,3)=SGD(IBM,1)
   20    CONTINUE
      ENDIF
      CALL LCMLEN(KPMACR,'DIFFX',LENGT2,ITYLCM)
      IF(LENGT2.GT.0) THEN
         IF(LENGT2.GT.NBMIX) CALL XABORT('TRISYS: INVALID LENGTH FOR'
     1   //' DIFFX (ANISOTROPIC DIFFUSION COEFFICIENT).')
         CALL LCMGET(KPMACR,'DIFFX',SGD(1,1))
         DO 30 IBM=1,LENGT2
         SGD(IBM,2)=SGD(IBM,1)
         SGD(IBM,3)=SGD(IBM,1)
   30    CONTINUE
      ENDIF
      CALL LCMLEN(KPMACR,'DIFFY',LENGT3,ITYLCM)
      IF(LENGT3.GT.0) THEN
         IF(LENGT3.GT.NBMIX) CALL XABORT('TRISYS: INVALID LENGTH FOR'
     1   //' DIFFY (ANISOTROPIC DIFFUSION COEFFICIENT).')
         CALL LCMGET(KPMACR,'DIFFY',SGD(1,2))
      ENDIF
      CALL LCMLEN(KPMACR,'DIFFZ',LENGT3,ITYLCM)
      IF(LENGT3.GT.0) THEN
         IF(LENGT3.GT.NBMIX) CALL XABORT('TRISYS: INVALID LENGTH FOR'
     1   //' DIFFZ (ANISOTROPIC DIFFUSION COEFFICIENT).')
         CALL LCMGET(KPMACR,'DIFFZ',SGD(1,3))
      ENDIF
      IF((LENGT1.EQ.0).AND.(LENGT2.EQ.0)) THEN
         CALL XABORT('TRISYS: NO DIFFUSION COEFFICIENTS.')
      ENDIF
      WRITE(TEXDIG,'(1HA,2I3.3)') IGR,IGR
      IF(IPR.EQ.0) THEN
*        COMPUTE UNPERTURBED SYSTEM MATRICES.
         DO 35 IBM=1,NBMIX
         IF((SGD(IBM,1).LT.0.0).OR.(SGD(IBM,4).LT.0.0)) THEN
            WRITE(HSMG,'(28HTRISYS: NEGATIVE XS IN GROUP,I5)') IGR
            CALL XABORT(HSMG)
         ENDIF
   35    CONTINUE
         CALL TRIASM(TEXDIG,IPTRK,IPSYS,IMPX,NBMIX,NEL,NALBP,0,MAT,VOL,
     1   GAMMA(1,IGR),SGD,SGD)
      ELSE
*        COMPUTE A PERTURBATION TO THE SYSTEM MATRICES
         DO 45 J=1,4
         DO 40 IBM=1,NBMIX
         DSGD(IBM,J)=0.0
   40    CONTINUE
   45    CONTINUE
         CALL LCMLEN(KPMACP,'NTOT0',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            IF(LENGT.GT.NBMIX) CALL XABORT('TRISYS: INVALID LENGTH FOR'
     1      //' DELTA TOTAL CROSS SECTIONS.')
            CALL LCMGET(KPMACP,'NTOT0',DSGD(1,4))
         ENDIF
         CALL LCMLEN(KPMACP,'SIGW00',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            IF(LENGT.GT.NBMIX) CALL XABORT('TRISYS: INVALID LENGTH FOR'
     1      //' DELTA ''SIGW00'' CROSS SECTIONS.')
            CALL LCMGET(KPMACP,'SIGW00',DSGD(1,1))
            DO 50 IBM=1,LENGT
            DSGD(IBM,4)=DSGD(IBM,4)-DSGD(IBM,1)
            DSGD(IBM,1)=0.0
   50       CONTINUE
         ENDIF
         CALL LCMLEN(KPMACP,'DIFF',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            IF(LENGT.GT.NBMIX) CALL XABORT('TRISYS: INVALID LENGTH FOR'
     1      //' DELTA DIFF (ISOTROPIC DIFFUSION COEFFICIENT).')
            CALL LCMGET(KPMACP,'DIFF',DSGD(1,1))
            DO 60 IBM=1,LENGT
            DSGD(IBM,2)=DSGD(IBM,1)
            DSGD(IBM,3)=DSGD(IBM,1)
   60       CONTINUE
         ENDIF
         CALL LCMLEN(KPMACP,'DIFFX',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            IF(LENGT.GT.NBMIX) CALL XABORT('TRISYS: INVALID LENGTH FOR'
     1      //' DELTA DIFFX (ANISOTROPIC DIFFUSION COEFFICIENT).')
            CALL LCMGET(KPMACP,'DIFFX',DSGD(1,1))
            CALL LCMGET(KPMACP,'DIFFY',DSGD(1,2))
            CALL LCMGET(KPMACP,'DIFFZ',DSGD(1,3))
         ENDIF
         CALL TRIASM(TEXDIG,IPTRK,IPSYS,IMPX,NBMIX,NEL,NALBP,IPR,MAT,
     1   VOL,GAMMA(1,IGR),SGD,DSGD)
      ENDIF
*----
*  PROCESS SCATTERING TERMS
*----
      CALL LCMLEN(KPMACP,'NJJS00',LENGT,ITYLCM)
      IF(LENGT.GT.NBMIX) CALL XABORT('TRISYS: INVALID LENGTH FOR ''N'
     1 //'JJS00'' INFORMATION.')
      IF(LENGT.GT.0) THEN
         CALL LCMGET(KPMACP,'NJJS00',NJJ)
         CALL LCMGET(KPMACP,'IJJS00',IJJ)
         JGRMIN=IGR
         JGRMAX=IGR
         DO 80 IBM=1,LENGT
         JGRMIN=MIN(JGRMIN,IJJ(IBM)-NJJ(IBM)+1)
         JGRMAX=MAX(JGRMAX,IJJ(IBM))
   80    CONTINUE
         CALL LCMGET(KPMACP,'IPOS00',IPOS)
         CALL LCMGET(KPMACP,'SCAT00',WORK)
         DO 100 JGR=JGRMAX,JGRMIN,-1
         IF(JGR.EQ.IGR) GO TO 100
         DO 90 IBM=1,LENGT
         IF((JGR.GT.IJJ(IBM)-NJJ(IBM)).AND.(JGR.LE.IJJ(IBM))) THEN
            SGD(IBM,1)=WORK(IPOS(IBM)+IJJ(IBM)-JGR)
         ELSE
            SGD(IBM,1)=0.0
         ENDIF
   90    CONTINUE
         WRITE (TEXDIG,'(1HA,2I3.3)') IGR,JGR
         CALL TRIDIG(TEXDIG,IPTRK,IPSYS,IMPX,NBMIX,NEL,IPR,MAT,
     1   VOL,SGD)
  100    CONTINUE
      ENDIF
  110 CONTINUE
*----
*  PROCESS FISSION SPECTRUM TERMS
*----
      KPMACP=LCMGIL(JPMACP,1)
      CALL LCMLEN(KPMACP,'CHI',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
         IF(LENGT.NE.NBMIX*NBFIS) CALL XABORT('TRISYS: INVALID LENGTH '
     1   //'FOR CHI INFORMATION.')
         DO 120 IGR=1,NGRP
         KPMACP=LCMGIL(JPMACP,IGR)
         CALL LCMGET(KPMACP,'CHI',CHI(1,1,IGR))
  120    CONTINUE
      ELSE
         DO 132 IBM=1,NBMIX
         DO 131 IFISS=1,NBFIS
         CHI(IBM,IFISS,1)=1.0
         DO 130 IGR=2,NGRP
         CHI(IBM,IFISS,IGR)=0.0
  130    CONTINUE
  131    CONTINUE
  132    CONTINUE
      ENDIF
*----
*  PROCESS FISSION NUSIGF TERMS
*----
      DO 170 IGR=1,NGRP
*     PROCESS SECONDARY GROUP IGR.
      LFIS=.FALSE.
      DO 141 IBM=1,NBMIX
      DO 140 IFISS=1,NBFIS
      LFIS=LFIS.OR.(CHI(IBM,IFISS,IGR).NE.0.0)
  140 CONTINUE
  141 CONTINUE
      IF(LFIS) THEN
         DO 160 JGR=1,NGRP
         KPMACP=LCMGIL(JPMACP,JGR)
         CALL LCMLEN(KPMACP,'NUSIGF',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            IF(LENGT.NE.NBMIX*NBFIS) CALL XABORT('TRISYS: INVALID LENG'
     1      //'TH FOR NUSIGF INFORMATION.')
            CALL LCMGET(KPMACP,'NUSIGF',ZUFIS)
            CALL XDRSET(SGD(1,1),NBMIX,0.0)
            DO 151 IBM=1,NBMIX
            DO 150 IFISS=1,NBFIS
            SGD(IBM,1)=SGD(IBM,1)+CHI(IBM,IFISS,IGR)*ZUFIS(IBM,IFISS)
  150       CONTINUE
  151       CONTINUE
            WRITE(TEXDIG,'(4HFISS,2I3.3)') IGR,JGR
            CALL LCMPUT(IPSYS,TEXDIG,NBMIX,2,SGD(1,1))
            WRITE (TEXDIG,'(1HB,2I3.3)') IGR,JGR
            CALL TRIDIG(TEXDIG,IPTRK,IPSYS,IMPX,NBMIX,NEL,IPR,MAT,VOL,
     1      SGD)
         ENDIF
  160    CONTINUE
      ENDIF
  170 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAMMA,SGD,DSGD,WORK,CHI,ZUFIS)
      DEALLOCATE(IJJ,NJJ,IPOS)
      RETURN
      END
