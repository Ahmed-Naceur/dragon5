*DECK EDIDEP
      SUBROUTINE EDIDEP(IPRINT,IPLIB,IPEDIT,NBNISO,HNNRF,ILNRF,IEVOL,
     1           LISO,NBCH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Create the 'DEPL-CHAIN' directory on the edition LCM object.
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
* IPRINT  print parameter.
* IPLIB   pointer to the internal library LCM object.
* IPEDIT  pointer to the edition LCM object.
* NBNISO  number of available isotopes in the edition LCM object.
* HNNRF   reference names of the available isotopes in the edition
*         LCM object.
* ILNRF   selection flag of the available isotopes in the edition
*         LCM object (=1 if selected).
* IEVOL   flag making an isotope non-depleting:
*         =1 to force an isotope to be non-depleting.
* LISO    =.true. if we want to register each isotope after merging.
*
*Parameters: output
* NBCH    number of depleting nuclides after lumping
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPEDIT
      INTEGER IPRINT,NBNISO,HNNRF(3,NBNISO),ILNRF(NBNISO),IEVOL(NBNISO),
     & NBCH
      LOGICAL LISO
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,MAXBCH=500)
      INTEGER ISTATE(NSTATE),HICH(3,MAXBCH)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MYLIS,IHREAC,IDREA,IPREA
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IHISO
      REAL, ALLOCATABLE, DIMENSION(:) :: DENER,DDECA,PRATE,YIELD
*----
*  FIND THE DEPLETING ISOTOPES IN THE EDITION MICROLIB
*----
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      NBISO=ISTATE(1)
      NBFISS=ISTATE(2)
      NBDPF=ISTATE(3)
      NSUPS=ISTATE(7)
      NREAC=ISTATE(8)
      NFATH=ISTATE(9)
      ALLOCATE(IHISO(3,NBISO))
      CALL LCMGET(IPLIB,'ISOTOPESDEPL',IHISO)
*     WE HAVE TO REGISTER SEVERAL TIMES THE SAME ISOTOPE IN THE NEW
*     DEPL-CHAIN IF WE WANT IT TO DEPLETE
      NBCH=0
      DO 20 ISO=1,NBISO
      DO JSO=1,NBNISO
        IF((ILNRF(JSO).EQ.0).OR.(IEVOL(JSO).EQ.1)) CYCLE
        IF((IHISO(1,ISO).EQ.HNNRF(1,JSO)).AND.
     &     (IHISO(2,ISO).EQ.HNNRF(2,JSO))) THEN
          IF(LISO) THEN
            NBCH=NBCH+1
            IF(NBCH.GT.MAXBCH) CALL XABORT('EDIDEP: MAXBCH OVERFLOW(1)')
            HICH(1,NBCH)=IHISO(1,ISO)
            HICH(2,NBCH)=IHISO(2,ISO)
          ELSE
            GO TO 10
          ENDIF
        ENDIF
      ENDDO
      GO TO 20
   10 IF(.NOT.LISO) THEN
        NBCH=NBCH+1
        IF(NBCH.GT.MAXBCH) CALL XABORT('EDIDEP: MAXBCH OVERFLOW(2)')
        HICH(1,NBCH)=IHISO(1,ISO)
        HICH(2,NBCH)=IHISO(2,ISO)
      ENDIF
   20 CONTINUE
*----
*  GENERATE THE DEPLETION INFORMATION CORRESPONDING TO THE AVAILABLE
*  ISOTOPES
*----
      IF(NBCH.GT.0) THEN
         MAXFP=NBDPF+30 ! reserve 30 location for lumped fp daughters
         NBFPCH=NBCH
         ALLOCATE(MYLIS(NBISO),IHREAC(2*NREAC),IDREA(NREAC*NBISO),
     1   DENER(NREAC*NBISO),DDECA(NBISO),IPREA(NFATH*NBISO),
     2   PRATE(NFATH*NBISO),YIELD(NBFISS*MAXFP))
         CALL LCMGET(IPLIB,'CHARGEWEIGHT',MYLIS)
         CALL LCMGET(IPLIB,'DEPLETE-IDEN',IHREAC)
         CALL LCMGET(IPLIB,'DEPLETE-REAC',IDREA)
         CALL LCMGET(IPLIB,'DEPLETE-ENER',DENER)
         CALL LCMGET(IPLIB,'DEPLETE-DECA',DDECA)
         CALL LCMGET(IPLIB,'PRODUCE-REAC',IPREA)
         CALL LCMGET(IPLIB,'PRODUCE-RATE',PRATE)
         IF(NBFISS*NBDPF.GT.0) THEN
            CALL LCMGET(IPLIB,'FISSIONYIELD',YIELD)
         ENDIF
*
         CALL LCMSIX(IPEDIT,'DEPL-CHAIN',1)
         IF(LISO) THEN
           NBFISS2=NBFPCH
           NBFPCH2=NBFPCH
         ELSE
           NBFISS2=NBFISS
           NBFPCH2=NBFPCH
         ENDIF
         CALL EDILUM(IPRINT,IPEDIT,MAXFP,NBISO,NBFISS,NBDPF,NSUPS,
     &   NREAC,NFATH,NBCH,HICH,IHISO,MYLIS,IHREAC,IDREA,DENER,DDECA,
     &   IPREA,PRATE,YIELD,LISO,NBFISS2,NBFPCH2)
         CALL LCMSIX(IPEDIT,' ',2)
*
         DEALLOCATE(YIELD,PRATE,IPREA,DDECA,DENER,IDREA,IHREAC,MYLIS)
      ENDIF
      DEALLOCATE(IHISO)
      RETURN
      END
