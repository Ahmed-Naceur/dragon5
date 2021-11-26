*DECK TRIRCA
      SUBROUTINE TRIRCA(IPMACR,IPMACP,NGRP,NBMIX,NANI,LDIFF,IL,IPR,RCAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the RCAT removal matrix in SPN cases.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMACR  L_MACROLIB pointer to the unperturbed cross sections.
* IPMACP  L_MACROLIB pointer to the perturbed cross sections if
*         IPR.gt.0. Equal to IPMACR if IPR=0.
* NGRP    number of energy groups.
* NBMIX   total number of material mixtures in the macrolib.
* NANI    maximum scattering order recovered from tracking and macrolib.
* LDIFF   flag set to .true. to use 1/3D as 'NTOT1' cross sections.
* IL      scattering Legendre order.
* IPR     type of assembly:
*         =0: calculation of the system matrices;
*         =1: calculation of the derivative of these matrices;
*         =2: calculation of the first variation of these matrices;
*         =3: identical to IPR=2, but these variation are added to
*         unperturbed system matrices.
*
*Parameters: output
* RCAT    removal matrix.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMACR,IPMACP
      INTEGER NGRP,NBMIX,NANI,IL,IPR
      LOGICAL LDIFF
      DOUBLE PRECISION RCAT(NGRP,NGRP,NBMIX)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMACR,JPMACP,KPMACR,KPMACP
      CHARACTER TEXT12*12,CM*2,HSMG*131
      DOUBLE PRECISION OTH
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SGD
      PARAMETER(OTH=1.0D0/3.0D0)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(NBMIX),NJJ(NBMIX),IPOS(NBMIX))
      ALLOCATE(SGD(NBMIX,3),WORK(NBMIX*NGRP))
*      
      JPMACR=LCMGID(IPMACR,'GROUP')
      JPMACP=LCMGID(IPMACP,'GROUP')
      WRITE(CM,'(I2.2)') IL-1
      CALL XDDSET(RCAT,NGRP*NGRP*NBMIX,0.0D0)
      DO 100 IGR=1,NGRP
*     PROCESS SECONDARY GROUP IGR.
      KPMACP=LCMGIL(JPMACP,IGR)
      CALL XDRSET(SGD(1,1),NBMIX,0.0)
      CALL LCMLEN(KPMACP,'SIGW'//CM,LENGT,ITYLCM)
      IF((LENGT.GT.0).AND.(IL.LE.NANI)) THEN
         IF(LENGT.GT.NBMIX) CALL XABORT('TRIRCA: INVALID LENGTH FOR'
     1   //' SIGW'//CM//' CROSS SECTIONS.')
         CALL LCMGET(KPMACP,'SIGW'//CM,SGD(1,1))
      ENDIF
      WRITE(TEXT12,'(4HNTOT,I1)') MIN(IL-1,9)
      CALL LCMLEN(KPMACP,TEXT12,LENGT,ITYLCM)
      CALL LCMLEN(KPMACP,'NTOT1',LENGT1,ITYLCM)
      IF((IL.EQ.1).AND.(LENGT.NE.NBMIX)) CALL XABORT('TRIRCA: NO NTOT0'
     1 //' CROSS SECTIONS.')
      IF(MOD(IL-1,2).EQ.0) THEN
*        macroscopic total cross section in even-parity equations.
         IF(LENGT.EQ.NBMIX) THEN
            CALL LCMGET(KPMACP,TEXT12,SGD(1,2))
         ELSE
            CALL LCMGET(KPMACP,'NTOT0',SGD(1,2))
         ENDIF
         DO 10 IBM=1,NBMIX
         IF((SGD(IBM,2)-SGD(IBM,1).LT.0.0).AND.(IPR.EQ.0)) THEN
            WRITE(HSMG,'(28HTRIRCA: NEGATIVE XS IN GROUP,I5)') IGR
            CALL XABORT(HSMG)
         ENDIF
         RCAT(IGR,IGR,IBM)=SGD(IBM,2)-SGD(IBM,1)
   10    CONTINUE
      ELSE
*        macroscopic total cross section in odd-parity equations.
         IF(LDIFF) THEN
            CALL LCMLEN(KPMACP,'DIFF',LENGT,ITYLCM)
            IF(LENGT.EQ.0) CALL XABORT('TRIRCA: DIFFUSION COEFFICIENTS'
     1      //' EXPECTED IN THE MACROLIB.')
            IF(LENGT.GT.NBMIX) CALL XABORT('TRIRCA: INVALID LENGTH FOR'
     1      //' DIFFUSION COEFFICIENTS.')
            CALL LCMGET(KPMACP,'DIFF',SGD(1,2))
            IF(IPR.EQ.0) THEN
              DO 20 IBM=1,NBMIX
              RCAT(IGR,IGR,IBM)=OTH/SGD(IBM,2)-SGD(IBM,1)
   20         CONTINUE
            ELSE IF(IPR.EQ.1) THEN
              KPMACR=LCMGIL(JPMACR,IGR)
              CALL LCMGET(KPMACR,'DIFF',SGD(1,3))
              DO 30 IBM=1,NBMIX
              RCAT(IGR,IGR,IBM)=-OTH*SGD(IBM,2)/SGD(IBM,3)**2-SGD(IBM,1)
   30         CONTINUE
            ELSE IF(IPR.EQ.2) THEN
              KPMACR=LCMGIL(JPMACR,IGR)
              CALL LCMGET(KPMACR,'DIFF',SGD(1,3))
              DO 40 IBM=1,NBMIX
              RCAT(IGR,IGR,IBM)=OTH/(SGD(IBM,2)+SGD(IBM,3))
     1                         -OTH/SGD(IBM,3)-SGD(IBM,1)
   40         CONTINUE
            ELSE IF(IPR.EQ.3) THEN
              KPMACR=LCMGIL(JPMACR,IGR)
              CALL LCMGET(KPMACR,'DIFF',SGD(1,3))
              DO 50 IBM=1,NBMIX
              RCAT(IGR,IGR,IBM)=OTH/(SGD(IBM,2)+SGD(IBM,3))-SGD(IBM,1)
   50         CONTINUE
            ENDIF
         ELSE
            IF(LENGT.EQ.NBMIX) THEN
               CALL LCMGET(KPMACP,TEXT12,SGD(1,2))
            ELSE IF(LENGT1.EQ.NBMIX) THEN
               CALL LCMGET(KPMACP,'NTOT1',SGD(1,2))
            ELSE
               CALL LCMGET(KPMACP,'NTOT0',SGD(1,2))
            ENDIF
            DO 60 IBM=1,NBMIX
            RCAT(IGR,IGR,IBM)=SGD(IBM,2)-SGD(IBM,1)
   60       CONTINUE
         ENDIF
         IF(IPR.EQ.0) THEN
            DO 65 IBM=1,NBMIX
            IF(RCAT(IGR,IGR,IBM).LT.0.0) THEN
               WRITE(HSMG,'(39HTRIRCA: INVALID CROSS-SECTION DATA (IL=,
     1         I3,2H).)') IL
               CALL XABORT(HSMG)
            ENDIF
   65       CONTINUE
         ENDIF
      ENDIF
      CALL LCMLEN(KPMACP,'NJJS'//CM,LENGT,ITYLCM)
      IF(LENGT.GT.NBMIX) CALL XABORT('TRIRCA: INVALID LENGTH FOR NJJS'
     1 //CM//' INFORMATION.')
      IF((LENGT.GT.0).AND.(IL.LE.NANI)) THEN
         CALL LCMGET(KPMACP,'NJJS'//CM,NJJ)
         CALL LCMGET(KPMACP,'IJJS'//CM,IJJ)
         IGMIN=IGR
         IGMAX=IGR
         DO 70 IBM=1,NBMIX
         IGMIN=MIN(IGMIN,IJJ(IBM)-NJJ(IBM)+1)
         IGMAX=MAX(IGMAX,IJJ(IBM))
   70    CONTINUE
         CALL LCMGET(KPMACP,'IPOS'//CM,IPOS)
         CALL LCMGET(KPMACP,'SCAT'//CM,WORK)
         DO 90 JGR=IGMAX,IGMIN,-1
         IF(JGR.EQ.IGR) GO TO 90
         DO 80 IBM=1,NBMIX
         IF((JGR.GT.IJJ(IBM)-NJJ(IBM)).AND.(JGR.LE.IJJ(IBM))) THEN
            RCAT(IGR,JGR,IBM)=-WORK(IPOS(IBM)+IJJ(IBM)-JGR)
         ENDIF
   80    CONTINUE
   90    CONTINUE
      ENDIF
  100 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK,SGD)
      DEALLOCATE(IPOS,NJJ,IJJ)
      RETURN
      END
