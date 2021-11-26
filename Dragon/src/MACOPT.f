*DECK MACOPT
      SUBROUTINE MACOPT(IPMAC,IPOPT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Update a Macrolib using control variables from a L_OPTIMIZE object.
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
* IPMAC   pointer to the macrolib to be updated.
* IPOPT   pointer to the L_OPTIMIZE object open in read-only mode.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC,IPOPT
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      TYPE(C_PTR) JPMAC,KPMAC
      INTEGER ISTATE(NSTATE)
*----
* ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SPH,ALB
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: VARV
*----
*  GET L_OPTIMIZE INFORMATION
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPOPT,'DEL-STATE',ISTATE)
      NGRP=ISTATE(1)
      NMIX=ISTATE(2)
      ITYPE=ISTATE(3)
      IDELTA=ISTATE(4)
      NGR1=ISTATE(5)
      NGR2=ISTATE(6)
      IBM1=ISTATE(7)
      IBM2=ISTATE(8)
*----
*  GET MACROLIB INFORMATION
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.NGRP) CALL XABORT('MACOPT: INVALID NGRP.')
      IF(ISTATE(2).NE.NMIX) CALL XABORT('MACOPT: INVALID NMIX.')
      NIFISS=ISTATE(4)
      NED=ISTATE(5)
      NALBP=ISTATE(8)
      ILEAK=ISTATE(9)
      IF(ITYPE.EQ.2) THEN
         ISTATE(10)=MAX(1,ISTATE(10))
         CALL LCMPUT(IPMAC,'STATE-VECTOR',NSTATE,1,ISTATE)
      ENDIF
      NPERT=(NGR2-NGR1+1)*(NALBP+IBM2-IBM1+1)
      IF(IDELTA.EQ.5) NPERT=(NGR2-NGR1+1)*NALBP
*----
*  CORRECT MACROLIB
*----
      ALLOCATE(VARV(NPERT))
      CALL LCMGET(IPOPT,'VAR-VALUE',VARV)
      IF(IDELTA.LE.2) THEN
*----
*  UPDATE ONLY LEAKAGE INFORMATION IN MACROLIB
*----
         ALLOCATE(GAR(3*NMIX))
         JPMAC=LCMGID(IPMAC,'GROUP')
         IPERT=0
         DO 70 IGR=NGR1,NGR2
         KPMAC=LCMDIL(JPMAC,IGR)
         IF((IDELTA.EQ.1).AND.(ITYPE.EQ.1).AND.(ILEAK.EQ.1)) THEN
            CALL LCMGET(KPMAC,'DIFF',GAR)
            DO 10 IBM=IBM1,IBM2
            IPERT=IPERT+1
            GAR(IBM)=REAL(VARV(IPERT))
   10       CONTINUE
            CALL LCMPUT(KPMAC,'DIFF',NMIX,2,GAR)
         ELSE IF((IDELTA.EQ.1).AND.(ITYPE.EQ.1).AND.(ILEAK.EQ.2)) THEN
            CALL LCMGET(KPMAC,'DIFFX',GAR)
            CALL LCMGET(KPMAC,'DIFFY',GAR(NMIX+1))
            CALL LCMGET(KPMAC,'DIFFZ',GAR(2*NMIX+1))
            DO 20 IBM=IBM1,IBM2
            IPERT=IPERT+1
            GAR(IBM)=REAL(VARV(IPERT))
            GAR(NMIX+IBM)=REAL(VARV(IPERT))
            GAR(2*NMIX+IBM)=REAL(VARV(IPERT))
   20       CONTINUE
            CALL LCMPUT(KPMAC,'DIFFX',NMIX,2,GAR)
            CALL LCMPUT(KPMAC,'DIFFY',NMIX,2,GAR(NMIX+1))
            CALL LCMPUT(KPMAC,'DIFFZ',NMIX,2,GAR(2*NMIX+1))
         ELSE IF((IDELTA.EQ.1).AND.(ITYPE.EQ.2)) THEN
            CALL LCMLEN(KPMAC,'NTOT1',ILONG,ITYLCM)
            IF(ILONG.NE.0.0) THEN
               CALL LCMGET(KPMAC,'NTOT1',GAR)
            ELSE
               CALL LCMGET(KPMAC,'NTOT0',GAR)
            ENDIF
            DO 30 IBM=IBM1,IBM2
            IPERT=IPERT+1
            GAR(IBM)=REAL(VARV(IPERT))
   30       CONTINUE
            CALL LCMPUT(KPMAC,'NTOT1',NMIX,2,GAR)
         ELSE IF((IDELTA.EQ.2).AND.(ITYPE.EQ.1).AND.(ILEAK.EQ.1)) THEN
            CALL LCMGET(KPMAC,'DIFF',GAR)
            DO 40 IBM=IBM1,IBM2
            IPERT=IPERT+1
            FACT=REAL(VARV(IPERT))
            GAR(IBM)=GAR(IBM)*FACT
   40       CONTINUE
            CALL LCMPUT(KPMAC,'DIFF',NMIX,2,GAR)
         ELSE IF((IDELTA.EQ.2).AND.(ITYPE.EQ.1).AND.(ILEAK.EQ.2)) THEN
            CALL LCMGET(KPMAC,'DIFFX',GAR)
            CALL LCMGET(KPMAC,'DIFFY',GAR(NMIX+1))
            CALL LCMGET(KPMAC,'DIFFZ',GAR(2*NMIX+1))
            DO 50 IBM=IBM1,IBM2
            IPERT=IPERT+1
            FACT=REAL(VARV(IPERT))
            GAR(IBM)=GAR(IBM)*FACT
            GAR(NMIX+IBM)=GAR(NMIX+IBM)*FACT
            GAR(2*NMIX+IBM)=GAR(2*NMIX+IBM)*FACT
   50       CONTINUE
            CALL LCMPUT(KPMAC,'DIFFX',NMIX,2,GAR)
            CALL LCMPUT(KPMAC,'DIFFY',NMIX,2,GAR(NMIX+1))
            CALL LCMPUT(KPMAC,'DIFFZ',NMIX,2,GAR(2*NMIX+1))
         ELSE IF((IDELTA.EQ.2).AND.(ITYPE.EQ.2)) THEN
            CALL LCMLEN(KPMAC,'NTOT1',ILONG,ITYLCM)
            IF(ILONG.NE.0.0) THEN
               CALL LCMGET(KPMAC,'NTOT1',GAR)
            ELSE
               CALL LCMGET(KPMAC,'NTOT0',GAR)
            ENDIF
            DO 60 IBM=IBM1,IBM2
            IPERT=IPERT+1
            FACT=REAL(VARV(IPERT))
            GAR(IBM)=GAR(IBM)*FACT
   60       CONTINUE
            CALL LCMPUT(KPMAC,'NTOT1',NMIX,2,GAR)
         ENDIF
   70    CONTINUE
         DEALLOCATE(GAR)
      ELSE IF(IDELTA.EQ.5) THEN
*----
*  CORRECT ONLY THE ALBEDO
*----
        ALLOCATE(ALB(NALBP,NGRP))
        CALL LCMGET(IPMAC,'ALBEDO',ALB)
        IPERT=0
        DO 90 IGR=NGR1,NGR2
        DO 80 IAL=1,NALBP
        IPERT=IPERT+1
        FACT=0.5*(1.0-ALB(IAL,IGR))/(1.0+ALB(IAL,IGR))*REAL(VARV(IPERT))
        ALB(IAL,IGR)=(1.0-2.0*FACT)/(1.0+2.0*FACT)
   80   CONTINUE
   90   CONTINUE
        CALL LCMPUT(IPMAC,'ALBEDO',NGRP*NALBP,2,ALB)
        DEALLOCATE(ALB)
      ELSE
*----
*  APPLY A FULL SPH CORRECTION
*----
        IPRINT=0
        IMC=IDELTA-2
        ALLOCATE(SPH(NMIX+NALBP,NGRP))
        CALL XDRSET(SPH,(NMIX+NALBP)*NGRP,1.0)
        IPERT=0
        DO 120 IGR=NGR1,NGR2
        DO 100 IBM=IBM1,IBM2
        IPERT=IPERT+1
        SPH(IBM,IGR)=REAL(VARV(IPERT))
  100   CONTINUE
        DO 110 IAL=1,NALBP
        IPERT=IPERT+1
        SPH(NMIX+IAL,IGR)=REAL(VARV(IPERT))
  110   CONTINUE
  120   CONTINUE
        CALL SPHCMA(IPMAC,IPRINT,IMC,NMIX,NGRP,NIFISS,NED,NALBP,SPH)
        DEALLOCATE(SPH)
      ENDIF
      DEALLOCATE(VARV)
      IF(IPERT.NE.NPERT) CALL XABORT('MACOPT: UPDATE FAILURE.')
      RETURN
      END
