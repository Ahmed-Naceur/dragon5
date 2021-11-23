*DECK COMDEP
      SUBROUTINE COMDEP(IPRINT,IPEDIT,IPWORK,ITRES,NISOP,NOMEVO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Creation of a lumped depletion chain in the multicompo.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input/output
* IPRINT  print parameter.
* IPEDIT  pointer to the edition object (L_EDIT signature).
* IPWORK  pointer to the LCM object where the lumped depletion chain is
*         written.
* ITRES   creation index for the macroscopic residual (=0: not created;
*         =1: not a FP precursor; =2: is a FP precursor).
* NISOP   number of user-requested particularized isotopes. Equal to
*         zero if all EDI: isotopes are particularized.
* NOMEVO  library names of user-requested particularized isotopes.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPRINT,ITRES,NISOP
      TYPE(C_PTR) IPEDIT,IPWORK
      CHARACTER NOMEVO(NISOP)*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,MAXBCH=500)
      INTEGER ISTATE(NSTATE),IHICH(3,MAXBCH)
      LOGICAL LISO
      CHARACTER TEXT12*12
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MYLIS,IHREAC
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IHISO,IDREA,IPREA
      REAL, ALLOCATABLE, DIMENSION(:) :: DDECA
      REAL, ALLOCATABLE, DIMENSION(:,:) :: DENER,PRATE,YIELD
*----
*  RECOVER DEPLETION INFORMATION FROM EDITION OBJECT
*----
      CALL LCMSIX(IPEDIT,'DEPL-CHAIN',1)
      IF(NISOP.GT.0) THEN
         CALL LCMGET(IPEDIT,'STATE-VECTOR',ISTATE)
         NBISO=ISTATE(1)
         IF(ITRES.EQ.2) NBISO=NBISO+1
         NBFISS=ISTATE(2)
         NBDPF=ISTATE(3)
         NSUPS=ISTATE(7)
         NREAC=ISTATE(8)
         NFATH=ISTATE(9)
         MAXFP=NBDPF+30 ! reserve 30 location for lumped fp daughters
         ALLOCATE(IHISO(3,NBISO),MYLIS(NBISO),IHREAC(2*NREAC),
     1   IDREA(NREAC,NBISO),DENER(NREAC,NBISO),DDECA(NBISO),
     2   IPREA(NFATH,NBISO),PRATE(NFATH,NBISO),YIELD(NBFISS,MAXFP))
         CALL LCMGET(IPEDIT,'ISOTOPESDEPL',IHISO)
         CALL LCMGET(IPEDIT,'CHARGEWEIGHT',MYLIS)
         CALL LCMGET(IPEDIT,'DEPLETE-IDEN',IHREAC)
         CALL LCMGET(IPEDIT,'DEPLETE-REAC',IDREA)
         CALL LCMGET(IPEDIT,'DEPLETE-ENER',DENER)
         CALL LCMGET(IPEDIT,'DEPLETE-DECA',DDECA)
         CALL LCMGET(IPEDIT,'PRODUCE-REAC',IPREA)
         CALL LCMGET(IPEDIT,'PRODUCE-RATE',PRATE)
         IF(NBFISS*NBDPF.GT.0) THEN
            CALL LCMGET(IPEDIT,'FISSIONYIELD',YIELD)
         ENDIF
*----
*  DESCRIBE FISSILE ISOTOPE *MAC*RES
*----
         IF(ITRES.EQ.2) THEN
           IF(IPRINT.GT.1) THEN
              WRITE(6,'(/42H COMDEP: ADD *MAC*RES RESIDUAL ISOTOPE TO ,
     1        17HDEPLETION CHAINS.)')
            ENDIF
            TEXT12='*MAC*RES'
            READ(TEXT12,'(3A4)') (IHISO(I0,NBISO),I0=1,3)
            MYLIS(NBISO)=0
            IDREA(:,NBISO)=0
            DENER(:,NBISO)=0.0
            IDREA(1,NBISO)=4
            DDECA(NBISO)=0.0
            IPREA(:,NBISO)=0
            PRATE(:,NBISO)=0.0
         ENDIF
*----
*  CREATE LUMPED DEPLETION CHAIN
*----
         CALL LCMSIX(IPWORK,'DEPL-CHAIN',1)
         LISO=.FALSE.
         NBCH=0
         DO 20 ISO=1,NBISO
         WRITE(TEXT12,'(3A4)') (IHISO(I0,ISO),I0=1,3)
         DO JSO=1,NISOP
           IF((TEXT12.EQ.NOMEVO(JSO)).AND.(TEXT12.NE.'*MAC*RES')) THEN
             NBCH=NBCH+1
             IF(NBCH.GT.MAXBCH) CALL XABORT('COMDEP: MAXBCH OVERFLOW.')
             READ(TEXT12,'(3A4)') (IHICH(I0,NBCH),I0=1,3)
             GO TO 20
           ENDIF
         ENDDO
         IF((TEXT12.EQ.'*MAC*RES').AND.(ITRES.EQ.2)) THEN
           NBCH=NBCH+1
           IF(NBCH.GT.MAXBCH) CALL XABORT('COMDEP: MAXBCH OVERFLOW.')
           READ(TEXT12,'(3A4)') (IHICH(I0,NBCH),I0=1,3)
         ENDIF
   20    CONTINUE
         CALL EDILUM(IPRINT,IPWORK,MAXFP,NBISO,NBFISS,NBDPF,NSUPS,
     1   NREAC,NFATH,NBCH,IHICH,IHISO,MYLIS,IHREAC,IDREA,DENER,DDECA,
     2   IPREA,PRATE,YIELD,LISO,NBFISS,NBCH)
         DEALLOCATE(YIELD,PRATE,IPREA,DDECA,DENER,IDREA,IHREAC,MYLIS,
     1   IHISO)
      ELSE
*----
*  RECOVER THE DEPLETION CHAIN WITHOUT LUMPING
*----
         CALL LCMSIX(IPWORK,'DEPL-CHAIN',1)
         CALL LCMEQU(IPEDIT,IPWORK)
      ENDIF
      CALL LCMSIX(IPWORK,' ',2)
      CALL LCMSIX(IPEDIT,' ',2)
      RETURN
      END
