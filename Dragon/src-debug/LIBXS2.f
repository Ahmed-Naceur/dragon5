*DECK LIBXS2
      SUBROUTINE LIBXS2(CFILNA,MAXR,NEL,NMDEPL,ITNAM,ITZEA,KPAX,BPAX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read depletion data on an APOLIB-XSM formatted library.
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input
* CFILNA  APOLIB-XSM file name.
* MAXR    number of reaction types.
* NEL     number of isotopes on library.
* NMDEPL  names of reactions:
*           NMDEPL(1)='DECAY'; NMDEPL(2)='NFTOT';
*           NMDEPL(3)='NG'   ; NMDEPL(4)='N2N';
*           etc.
*
*Parameters: output
* ITNAM   reactive isotope names in chain.
* ITZEA   6-digit nuclide identifier:
*         atomic number z*10000 (digits) + mass number a*10 +
*         energy state (0 = ground state, 1 = first state, etc.).
* KPAX    complete reaction type matrix.
* BPAX    complete branching ratio matrix.
*
*Comments:
*  INPUT FORMAT
*    LIB: APLIB2 FIL: CFILNA CHAIN
*    [[ hnamson
*    [ FROM  [[ { DECAY | reaction } yield hnampar ]] ]
*    ]]
*    ENDCHAIN
*-----------------------------------------------------------------------
*
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER CFILNA*(*),NMDEPL(MAXR)*8
      INTEGER MAXR,NEL,ITNAM(3,NEL),ITZEA(NEL),KPAX(NEL+MAXR,NEL)
      REAL BPAX(NEL+MAXR,NEL)
*
      TYPE(C_PTR) IPAP
      PARAMETER (IOUT=6)
      CHARACTER TEXT20*20,TEXT12*12,HNISOR*20,HITNAM*20,HSMG*131
      DOUBLE PRECISION DBLINP
      REAL E458(9)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NOM,IA,IZ,NFG,IKEEP
      REAL, ALLOCATABLE, DIMENSION(:) :: GAMMA,RTSEGM
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HSECTT
*----
*  OPEN APOLIB FILE
*----
      CALL LCMOP(IPAP,CFILNA,2,2,0)
*----
*  RECOVER INFORMATION FROM PHEAD DIRECTORY
*----
      CALL LCMSIX(IPAP,'PHEAD',1)
      CALL LCMLEN(IPAP,'NOM',NV,ITYLCM)
      NISOT=NV/5
      ALLOCATE(NOM(5*NISOT))
      CALL LCMGET(IPAP,'NOM',NOM)
      DO 20 ISO=1,NISOT
      WRITE(HNISOR,'(5A4)') (NOM((ISO-1)*5+II),II=1,5)
      READ(HNISOR,'(3A4)') (ITNAM(II,ISO),II=1,3)
   20 CONTINUE
      CALL LCMSIX(IPAP,' ',2)
*----
*  RECOVER INFORMATION FROM PCONST DIRECTORY
*----
      CALL LCMSIX(IPAP,'PCONST',1)
      CALL LCMLIB(IPAP)
      CALL LCMLEN(IPAP,'A',NV,ITYLCM)
      IF(NV.NE.NISOT) CALL XABORT('LIBXS2: IA OVERFLOW')
      ALLOCATE(IA(NISOT),IZ(NISOT),NFG(NISOT))
      CALL LCMGET(IPAP,'A',IA)
      CALL LCMGET(IPAP,'Z',IZ)
      CALL LCMGET(IPAP,'NFG',NFG)
      CALL LCMSIX(IPAP,' ',2)
*----
*  RECOVER INFORMATION FROM PNUMF DIRECTORY
*----
      CALL LCMSIX(IPAP,'PNUMF',1)
      CALL LCMLEN(IPAP,'GAMMA',NGAMMA,ITYLCM)
      CALL LCMLEN(IPAP,'NOMFIS',NBFISS,ITYLCM)
      CALL LCMLEN(IPAP,'NOMPF',NBPF,ITYLCM)
      NBFISS=NBFISS/2
      NBPF=NBPF/2
      ALLOCATE(GAMMA(NGAMMA))
      CALL LCMGET(IPAP,'GAMMA',GAMMA)
      NMGY=NGAMMA/(NBFISS*NBPF)
      CALL LCMSIX(IPAP,' ',2)
*----
*  LOOP OVER ISOTOPES
*----
      CALL LCMSIX(IPAP,'QFIX',1)
      DO 260 ISO=1,NISOT
      WRITE(HNISOR,'(5A4)') (NOM((ISO-1)*5+II),II=1,5)
      WRITE(TEXT12,'(4HISOT,I8.8)') ISO
      CALL LCMSIX(IPAP,TEXT12,1)
      CALL LCMSIX(IPAP,'ISOTOP',1)
*     NG ENERGY.
      CALL LCMLEN(IPAP,'EGAMM',NV,ITYLCM)
      IF(NV.NE.0) THEN
         KPAX(NEL+3,ISO)=1
         CALL LCMGET(IPAP,'EGAMM',BPAX(NEL+3,ISO))
      ENDIF
*     FISSION ENERGIES.
      CALL LCMLEN(IPAP,'EF',NV,ITYLCM)
      IF(NV.NE.0) THEN
         KPAX(NEL+2,ISO)=1
         CALL LCMGET(IPAP,'EF',BPAX(NEL+2,ISO))
      ENDIF
      CALL LCMLEN(IPAP,'ENER_458',NV,ITYLCM)
      IF(NV.NE.0) THEN
         KPAX(NEL+2,ISO)=1
         CALL LCMGET(IPAP,'ENER_458',E458)
         BPAX(NEL+2,ISO)=E458(8)
      ENDIF
*     RADIOACTIVE DECAY CONSTANTS.
      CALL LCMLEN(IPAP,'LAMBD0',NCHANN,ITYLCM)
      IF(NCHANN.GT.0) THEN
        ALLOCATE(RTSEGM(NCHANN))
        CALL LCMGET(IPAP,'LAMBD0',RTSEGM)
        SUM=0.0
        DO 140 I=1,NCHANN
        SUM=SUM+RTSEGM(I)
  140   CONTINUE
        DEALLOCATE(RTSEGM)
        IF(SUM.NE.0.0) BPAX(NEL+1,ISO)=SUM*1.0E8
      ENDIF
*     X-S NAMES.
      CALL LCMLEN(IPAP,'TYSECT',NV,ITYLCM)
      NSECTT=NV/2
      ALLOCATE(HSECTT(NSECTT))
      CALL LCMGTC(IPAP,'TYSECT',8,NSECTT,HSECTT)
      DO 150 IS=1,NSECTT
      IF(HSECTT(IS).EQ.'SIGA') THEN
        KPAX(NEL+3,ISO)=1
      ELSE IF(HSECTT(IS).EQ.'NEXCESS') THEN
        KPAX(NEL+4,ISO)=1
      ELSE IF(HSECTT(IS).EQ.'SIGF') THEN
        KPAX(NEL+2,ISO)=1
      ELSE IF(HSECTT(IS).EQ.'CREA-A') THEN
        KPAX(NEL+7,ISO)=1
      ELSE IF(HSECTT(IS).EQ.'CREA-P') THEN
        KPAX(NEL+8,ISO)=1
      ELSE IF(HSECTT(IS).EQ.'CREA-H2') THEN
        KPAX(NEL+11,ISO)=1
      ELSE IF(HSECTT(IS).EQ.'CREA-H3') THEN
        KPAX(NEL+12,ISO)=1
      ENDIF
  150 CONTINUE
      DEALLOCATE(HSECTT)
*----
*  SET OTHER INFORMATION.
*----
      ITZEA(ISO)=IZ(ISO)*10000+IA(ISO)*10
      IPF=NFG(ISO)
      IF(IPF.LT.0) THEN
         KPAX(NEL+2,ISO)=-1
         DO 250 JSO=1,NISOT
         IFI=NFG(JSO)
         IF(IFI.GT.0) THEN
            IOFSET=((-IPF-1)*NBFISS+(IFI-1))*NMGY+NMGY
            IF(IOFSET.GT.NGAMMA) CALL XABORT('LIBXS2: GAMMA OVERFLOW.')
            BPAX(ISO,JSO)=GAMMA(IOFSET)
            IF(BPAX(ISO,JSO).NE.0.0) KPAX(ISO,JSO)=2
         ENDIF
  250    CONTINUE
      ENDIF
      CALL LCMSIX(IPAP,' ',2)
      CALL LCMSIX(IPAP,' ',2)
  260 CONTINUE
      CALL LCMSIX(IPAP,' ',2)
*
      DEALLOCATE(GAMMA,NFG,IZ,IA,NOM)
      CALL LCMCL(IPAP,1)
*----
*  RECOVER INFORMATION FROM INPUT DATA STREAM.
*----
      ALLOCATE(IKEEP(NEL))
      CALL XDISET(IKEEP,NEL,0)
      TEXT12=' '
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
      IF(INDIC.NE.3.OR.TEXT12.NE.'CHAIN')
     >  CALL XABORT('LIBXS2: KEYWORD CHAIN MISSING')
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
      DO 340 IEL=1,NEL
      IF(TEXT12.EQ.'ENDCHAIN') GO TO 350
      IF(INDIC.NE.3) CALL XABORT('LIBXS2: ISOTOPE NAME hnamson MISSING')
      I1=INDEX(TEXT12,'_')
      HNISOR=' '
      IF(I1.EQ.0) THEN
         HNISOR(:12)=TEXT12
      ELSE
         HNISOR(:I1-1)=TEXT12(:I1-1)
      ENDIF
      IDEPL=0
      DO 270 JEL=1,NEL
      WRITE(TEXT12,'(3A4)') (ITNAM(II,JEL),II=1,3)
      I1=INDEX(TEXT12,'_')
      HITNAM=' '
      IF(I1.EQ.0) THEN
         HITNAM(:12)=TEXT12
       ELSE
         HITNAM(:I1-1)=TEXT12(:I1-1)
      ENDIF
      IF(HNISOR.EQ.HITNAM) THEN
         IDEPL=JEL
         GO TO 280
      ENDIF
  270 CONTINUE
      WRITE(HSMG,'(25HLIBXS2: MISSING ISOTOPE '',A12,5H''(1).)')
     > HNISOR
      CALL XABORT(HSMG)
  280 IKEEP(IDEPL)=1
      IF(BPAX(NEL+1,IDEPL).NE.0.0) KPAX(NEL+1,IDEPL)=1
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
      IF(INDIC.NE.3) CALL XABORT('LIBXS2: REACTION TYPE EXPECTED')
      IF(TEXT12.EQ.'FROM') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
  290    IF(INDIC.NE.3) CALL XABORT('LIBXS2: REACTION TYPE EXPECTED')
         DO 330 IREAC=1,MAXR
         RRAT=1.0
         IF(TEXT12.EQ.NMDEPL(IREAC)) THEN
            DO 320 JEL=1,NEL
            CALL REDGET(INDIC,ISOT,RRAT,TEXT12,DBLINP)
            IF(INDIC.NE.2) GO TO 290
            CALL REDGET(INDIC,ISOT,FLOTT,TEXT12,DBLINP)
            IF(INDIC.NE.3) CALL XABORT('LIBXS2: ISOTOPE NAME HNAMPAR '
     >      //'MISSING')
            I1=INDEX(TEXT12,'_')
            TEXT20=' '
            IF(I1.EQ.0) THEN
               TEXT20(:12)=TEXT12
            ELSE
               TEXT20(:I1-1)=TEXT12(:I1-1)
            ENDIF
            JDEPL=0
            DO 300 JREL=1,NEL
            WRITE(TEXT12,'(3A4)') (ITNAM(II,JREL),II=1,3)
            I1=INDEX(TEXT12,'_')
            HITNAM=' '
            IF(I1.EQ.0) THEN
               HITNAM(:12)=TEXT12
            ELSE
               HITNAM(:I1-1)=TEXT12(:I1-1)
            ENDIF
            IF(TEXT20.EQ.HITNAM) THEN
               JDEPL=JREL
               GO TO 310
            ENDIF
  300       CONTINUE
            WRITE(HSMG,'(25HLIBXS2: MISSING ISOTOPE '',A12,5H''(2).)')
     >      TEXT20
            CALL XABORT(HSMG)
  310       KPAX(IDEPL,JDEPL)=IREAC
            BPAX(IDEPL,JDEPL)=RRAT
  320       CONTINUE
            CALL XABORT('LIBXS2: TO MANY PARENT ISOTOPES')
         ENDIF
  330    CONTINUE
      ENDIF
  340 CONTINUE
      IF(INDIC.NE.3.OR.TEXT12.NE.'ENDCHAIN')
     >  CALL XABORT('LIBXS2: KEYWORD ENDCHAIN MISSING')
  350 DO 380 JEL=1,NEL
      IF(IKEEP(JEL).EQ.0) THEN
         DO 360 IREAC=1,NEL+MAXR
         KPAX(IREAC,JEL)=0
  360    CONTINUE
         DO 370 IEL=1,NEL
         KPAX(JEL,IEL)=0
  370    CONTINUE
      ENDIF
  380 CONTINUE
      DEALLOCATE(IKEEP)
      RETURN
      END
