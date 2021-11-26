*DECK TRISPS
      SUBROUTINE TRISPS(IPTRK,IPMACR,IPMACP,IPSYS,IMPX,NGRP,NEL,NLF,
     1 NANI,NBFIS,NALBP,LDIFF,IPR,MAT,VOL,NBMIX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the cross-section data in LCM object with pointer IPMACR,
* compute and store the corresponding Trivac system matrices for a
* simplified PN approximation (or a perturbation to the system
* matrices).
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
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
* NLF     number of Legendre orders for the flux (even number).
* NANI    number of Legendre orders for the scattering cross sections.
* NBFIS   number of fissionable isotopes.
* NALBP   number of physical albedos per energy group.
* LDIFF   flag set to .true. to use 1/3D as 'NTOT1' cross sections.
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
      INTEGER IMPX,NGRP,NEL,NLF,NANI,NBFIS,NALBP,IPR,MAT(NEL),NBMIX
      REAL VOL(NEL)
      LOGICAL LDIFF
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXDIG*12,TEXT12*12,CM*2
      LOGICAL LFIS
      TYPE(C_PTR) JPMACP,KPMACP
      REAL, DIMENSION(:), ALLOCATABLE :: WORK
      REAL, DIMENSION(:,:), ALLOCATABLE :: GAMMA,SGD,ZUFIS
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: CHI
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GAR
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: RCAT,RCATI,
     1 RCAT2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(GAMMA(NALBP,NGRP),SGD(NBMIX,2*NLF),WORK(NBMIX*NGRP),
     1 CHI(NBMIX,NBFIS,NGRP),ZUFIS(NBMIX,NBFIS))
      ALLOCATE(RCAT(NGRP,NGRP,NBMIX),RCATI(NGRP,NGRP,NBMIX))
*----
*  PROCESS PHYSICAL ALBEDOS.
*----
      IF(NALBP.GT.0) THEN
         CALL TRIALB(IPTRK,IPMACR,IPMACP,IPSYS,NGRP,NALBP,IPR,GAMMA)
      ENDIF
*----
*  PROCESS MACROLIB INFORMATION FOR VARIOUS LEGENDRE ORDERS AND
*  INVERSION OF THE REMOVAL MATRIX.
*----
      IF(NLF.EQ.0) CALL XABORT('TRISPS: SPN APPROXIMATION REQUESTED.')
      DO 142 IL=1,NLF
      WRITE(CM,'(I2.2)') IL-1
      CALL TRIRCA(IPMACR,IPMACR,NGRP,NBMIX,NANI,LDIFF,IL,0,RCAT)
      IF(IPR.EQ.0) THEN
         DO 20 IBM=1,NBMIX
         DO 15 JGR=1,NGRP
         DO 10 IGR=1,NGRP
         RCATI(IGR,JGR,IBM)=RCAT(IGR,JGR,IBM)
   10    CONTINUE
   15    CONTINUE
         CALL ALINVD(NGRP,RCATI(1,1,IBM),NGRP,IER)
         IF(IER.NE.0) CALL XABORT('TRISPS: SINGULAR MATRIX(1).')
   20    CONTINUE
      ELSE
         ALLOCATE(RCAT2(NGRP,NGRP,NBMIX),GAR(NGRP))
         CALL TRIRCA(IPMACR,IPMACP,NGRP,NBMIX,NANI,LDIFF,IL,IPR,RCAT2)
         IF(IPR.EQ.1) THEN
            DO 62 IBM=1,NBMIX
            DO 31 JGR=1,NGRP
            DO 30 IGR=1,NGRP
            RCATI(IGR,JGR,IBM)=RCAT(IGR,JGR,IBM)
            RCAT(IGR,JGR,IBM)=RCAT2(IGR,JGR,IBM)
   30       CONTINUE
   31       CONTINUE
            CALL ALINVD(NGRP,RCATI(1,1,IBM),NGRP,IER)
            IF(IER.NE.0) CALL XABORT('TRISPS: SINGULAR MATRIX(2).')
            DO 42 JGR=1,NGRP
            CALL XDDSET(RCAT2(1,JGR,IBM),NGRP,0.0D0)
            DO 41 IGR=1,NGRP
            DO 40 KGR=1,NGRP
            RCAT2(IGR,JGR,IBM)=RCAT2(IGR,JGR,IBM)+RCATI(IGR,KGR,IBM)*
     1      RCAT(KGR,JGR,IBM)
   40       CONTINUE
   41       CONTINUE
   42       CONTINUE
            DO 61 JGR=1,NGRP
            CALL XDDSET(GAR,NGRP,0.0D0)
            DO 51 IGR=1,NGRP
            DO 50 KGR=1,NGRP
            GAR(IGR)=GAR(IGR)+RCAT2(IGR,KGR,IBM)*RCATI(KGR,JGR,IBM)
   50       CONTINUE
   51       CONTINUE
            DO 60 KGR=1,NGRP
            RCATI(KGR,JGR,IBM)=-GAR(KGR)
   60       CONTINUE
   61       CONTINUE
   62       CONTINUE
         ELSE IF(IPR.EQ.2) THEN
            DO 82 IBM=1,NBMIX
            DO 71 JGR=1,NGRP
            DO 70 IGR=1,NGRP
            RCATI(IGR,JGR,IBM)=RCAT(IGR,JGR,IBM)
            RCAT(IGR,JGR,IBM)=RCAT2(IGR,JGR,IBM)
            RCAT2(IGR,JGR,IBM)=RCAT(IGR,JGR,IBM)+RCATI(IGR,JGR,IBM)
   70       CONTINUE
   71       CONTINUE
            CALL ALINVD(NGRP,RCATI(1,1,IBM),NGRP,IER)
            IF(IER.NE.0) CALL XABORT('TRISPS: SINGULAR MATRIX(3).')
            CALL ALINVD(NGRP,RCAT2(1,1,IBM),NGRP,IER)
            IF(IER.NE.0) CALL XABORT('TRISPS: SINGULAR MATRIX(4).')
            DO 81 JGR=1,NGRP
            DO 80 IGR=1,NGRP
            RCATI(IGR,JGR,IBM)=RCAT2(IGR,JGR,IBM)-RCATI(IGR,JGR,IBM)
   80       CONTINUE
   81       CONTINUE
   82       CONTINUE
         ELSE IF(IPR.EQ.3) THEN
            DO 100 IBM=1,NBMIX
            DO 91 JGR=1,NGRP
            DO 90 IGR=1,NGRP
            RCAT(IGR,JGR,IBM)=RCAT(IGR,JGR,IBM)+RCAT2(IGR,JGR,IBM)
            RCATI(IGR,JGR,IBM)=RCAT(IGR,JGR,IBM)
   90       CONTINUE
   91       CONTINUE
            CALL ALINVD(NGRP,RCATI(1,1,IBM),NGRP,IER)
            IF(IER.NE.0) CALL XABORT('TRISPS: SINGULAR MATRIX(5).')
  100       CONTINUE
         ENDIF
         DEALLOCATE(GAR,RCAT2)
      ENDIF
*
      DO 141 IGR=1,NGRP
      IGMIN=IGR
      IGMAX=IGR
      DO 111 IBM=1,NBMIX
      DO 110 JGR=1,NGRP
      IF((RCAT(IGR,JGR,IBM).NE.0.0).OR.(RCATI(IGR,JGR,IBM).NE.0.0)) THEN
         IGMIN=MIN(IGMIN,JGR)
         IGMAX=MAX(IGMAX,JGR)
      ENDIF
  110 CONTINUE
  111 CONTINUE
      DO 140 JGR=IGMIN,IGMAX
      DO 120 IBM=1,NBMIX
      WORK(IBM)=REAL(RCAT(IGR,JGR,IBM))
  120 CONTINUE
      WRITE(TEXT12,'(4HSCAR,A2,2I3.3)') CM,IGR,JGR
      CALL LCMPUT(IPSYS,TEXT12,NBMIX,2,WORK)
      DO 130 IBM=1,NBMIX
      WORK(IBM)=REAL(RCATI(IGR,JGR,IBM))
  130 CONTINUE
      WRITE(TEXT12,'(4HSCAI,A2,2I3.3)') CM,IGR,JGR
      CALL LCMPUT(IPSYS,TEXT12,NBMIX,2,WORK)
  140 CONTINUE
  141 CONTINUE
  142 CONTINUE
*----
*  COMPUTE AND FACTORIZE THE DIAGONAL SYSTEM MATRICES.
*----
      DO 162 IGR=1,NGRP
      DO 150 IL=1,NLF
      WRITE(TEXT12,'(4HSCAR,I2.2,2I3.3)') IL-1,IGR,IGR
      CALL LCMGET(IPSYS,TEXT12,SGD(1,IL))
      WRITE(TEXT12,'(4HSCAI,I2.2,2I3.3)') IL-1,IGR,IGR
      CALL LCMGET(IPSYS,TEXT12,SGD(1,NLF+IL))
  150 CONTINUE
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
      CALL TRIASN(TEXT12,IPTRK,IPSYS,IMPX,NBMIX,NEL,NLF,NALBP,IPR,MAT,
     1 VOL,GAMMA(1,IGR),SGD(1,1),SGD(1,1+NLF))
*----
*  PUT A FLAG IN IPSYS TO IDENTIFY NON-ZERO SCATTERING TERMS.
*----
      DO 161 IL=1,NLF
      DO 160 JGR=1,NGRP
      WRITE(TEXT12,'(4HSCAR,I2.2,2I3.3)') IL-1,IGR,JGR
      CALL LCMLEN(IPSYS,TEXT12,LENGT,ITYLCM)
      IF(LENGT.EQ.NBMIX) THEN
         WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
         CALL LCMPUT(IPSYS,TEXT12,1,2,0.0)
      ENDIF
  160 CONTINUE
  161 CONTINUE
  162 CONTINUE
*----
*  PROCESS FISSION SPECTRUM TERMS
*----
      JPMACP=LCMGID(IPMACP,'GROUP')
      KPMACP=LCMGIL(JPMACP,1)
      CALL LCMLEN(KPMACP,'CHI',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
         IF(LENGT.NE.NBMIX*NBFIS) CALL XABORT('TRISPS: INVALID LENGTH '
     1   //'FOR CHI INFORMATION.')
         DO 180 IGR=1,NGRP
         KPMACP=LCMGIL(JPMACP,IGR)
         CALL LCMGET(KPMACP,'CHI',CHI(1,1,IGR))
  180    CONTINUE
      ELSE
         DO 192 IBM=1,NBMIX
         DO 191 IFISS=1,NBFIS
         CHI(IBM,IFISS,1)=1.0
         DO 190 IGR=2,NGRP
         CHI(IBM,IFISS,IGR)=0.0
  190    CONTINUE
  191    CONTINUE
  192    CONTINUE
      ENDIF
*----
*  PROCESS FISSION NUSIGF TERMS
*----
      DO 230 IGR=1,NGRP
*     PROCESS SECONDARY GROUP IGR.
      LFIS=.FALSE.
      DO 201 IBM=1,NBMIX
      DO 200 IFISS=1,NBFIS
      LFIS=LFIS.OR.(CHI(IBM,IFISS,IGR).NE.0.0)
  200 CONTINUE
  201 CONTINUE
      IF(LFIS) THEN
         DO 220 JGR=1,NGRP
         KPMACP=LCMGIL(JPMACP,JGR)
         CALL LCMLEN(KPMACP,'NUSIGF',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            IF(LENGT.NE.NBMIX*NBFIS) CALL XABORT('TRISPS: INVALID LENG'
     1      //'TH FOR NUSIGF INFORMATION.')
            CALL LCMGET(KPMACP,'NUSIGF',ZUFIS)
            CALL XDRSET(SGD(1,1),NBMIX,0.0)
            DO 211 IBM=1,NBMIX
            DO 210 IFISS=1,NBFIS
            SGD(IBM,1)=SGD(IBM,1)+CHI(IBM,IFISS,IGR)*ZUFIS(IBM,IFISS)
  210       CONTINUE
  211       CONTINUE
            WRITE(TEXDIG,'(4HFISS,2I3.3)') IGR,JGR
            CALL LCMPUT(IPSYS,TEXDIG,NBMIX,2,SGD(1,1))
            WRITE (TEXDIG,'(1HB,2I3.3)') IGR,JGR
            CALL TRIDIG(TEXDIG,IPTRK,IPSYS,IMPX,NBMIX,NEL,IPR,MAT,VOL,
     1      SGD)
         ENDIF
  220    CONTINUE
      ENDIF
  230 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(RCAT,RCATI)
      DEALLOCATE(GAMMA,SGD,WORK,CHI,ZUFIS)
      RETURN
      END
