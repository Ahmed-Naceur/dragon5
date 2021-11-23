*DECK EDILUM
      SUBROUTINE EDILUM(IPRINT,IPEDIT,MAXFP,NBISO,NBFISS,NBDPF,NSUPS,
     & NREAC,NFATH,NBCH,HICH,HISO,MYLIST,HREAC,IDREAC,DENER,DDECA,
     & IPREAC,PRATE,YIELD,LISO,NBFISS2,NBFPCH2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Complete and lump the burnup chain from NBISO to NBCH isotopes.
* Write the lumped chain on the LCM object. Based on subroutine
* dralum.f in dragr module.
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
* IPEDIT  pointer to the edition LCM object.
* MAXFP   second dimension of array 'YIELD'.
* NBISO   number of depleting nuclides before lumping.
* NBFISS  number of fissile isotopes with fission yields.
* NBDPF   number of fission products before lumping.
* NSUPS   number of stable isotopes producing energy.
* NREAC   number of depleting reactions including radioactive decay
* NFATH   maximum number of parent isotopes leading to the production of
*         an isotope in the depletion chain before lumping.
* NBCH    number of depleting nuclides after lumping.
* HICH    names of remaining isotopes after lumping.
* HISO    'ISOTOPESDEPL' record before lumping.
* MYLIST  'CHARGEWEIGHT' record before lumping.
* HREAC   'DEPLETE-IDEN' record before lumping.
* IDREAC  'DEPLETE-REAC' record before lumping.
* DENER   'DEPLETE-ENER' record before lumping.
* DDECA   'DEPLETE-DECA' record before lumping.
* IPREAC  'PRODUCE-REAC' record before lumping.
* PRATE   'PRODUCE-RATE' record before lumping.
* YIELD   'FISSIONYIELD' record before lumping.
* LISO    =.true. when all isotopes are kept separately during 
*         the merging.
* NBFISS2 new number of fissile isotopes with fission yields.
* NBFPCH2 new maximum number of fission products after lumping.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPEDIT
      INTEGER IPRINT,MAXFP,NBISO,NBFISS,NBDPF,NSUPS,NREAC,NFATH,
     & MYLIST(NBISO),HICH(3,NBCH),HISO(3,NBISO),NBCH,HREAC(2,NREAC),
     & IDREAC(NREAC,NBISO),IPREAC(NFATH,NBISO),NBFISS2,NBFPCH2
      REAL DENER(NREAC,NBISO),DDECA(NBISO),PRATE(NFATH,NBISO),
     & YIELD(NBFISS,MAXFP)
      LOGICAL LISO
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,MAXIT=20,NSYSO=6)
      CHARACTER HNAME*8,TEXT4*4
      INTEGER ISTATE(NSTATE),MIX(NBCH),IMIX
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: JPREAC,JDREAC,IPOS,HHHH
      REAL, ALLOCATABLE, DIMENSION(:) :: EYIEL,DDDD
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RRATE,EENER
*----
*  SCRATCH STORAGE ALLOCATION
*----
      PARAMETER (MAXFAT=25)
      ALLOCATE(JPREAC(MAXFAT,NBCH),JDREAC(NREAC,NBCH),IPOS(NBCH,3),
     & HHHH(3,NBCH))
      ALLOCATE(RRATE(MAXFAT,NBCH),EENER(NREAC,NBCH),
     & EYIEL(NBFISS2*NBFPCH2),DDDD(NBCH))
*----
*  FIND THE POSITION OF THE LUMPED ISOTOPES IN THE COMPLETE CHAIN
*----
      J0OLD=0
      J0=0
      DO ISO=1,NBCH
        DO JSO=1,NBISO
          J0=JSO
          IF((HICH(1,ISO).EQ.HISO(1,JSO)).AND.
     &       (HICH(2,ISO).EQ.HISO(2,JSO))) GO TO 10
        ENDDO
        WRITE(NSYSO,'(/35H EDILUM: LIST OF AVAILABLE ISOTOPES)')
        DO JSO=1,NBISO
          WRITE(NSYSO,'(2X,I5,5H --> ,2A4)') JSO,HISO(1,JSO),HISO(2,JSO)
        ENDDO
        WRITE(HNAME,'(2A4)') HICH(1,ISO),HICH(2,ISO)
        CALL XABORT('EDILUM: UNABLE TO FIND '//HNAME)
 10     IF(.NOT.LISO) THEN
          IF(J0.LE.J0OLD) THEN
             WRITE(HNAME,'(2A4)') HISO(1,J0),HISO(2,J0)
             CALL XABORT('EDILUM: ISOTOPIC DATA NOT SORTED:'//HNAME)
          ENDIF
        ENDIF
        J0OLD=J0
        IPOS(ISO,1)=J0
      ENDDO
      DO I=1,NBCH
        IMIX=1
        IF(I.EQ.1) GOTO 15
        DO J=1,I-1
          IF(IPOS(I,1).EQ.IPOS(J,1)) IMIX=IMIX+1
        ENDDO
 15     MIX(I)=IMIX
      ENDDO
*----
*  REMOVE THE LUMPED FATHERS THAT ARE NEUTRON INDUCED REACTIONS.
*----
      DO ISO=1,NBISO
        DO IFATH=1,NFATH
          IF(IPREAC(IFATH,ISO).EQ.0) GO TO 20
          IF(MOD(IPREAC(IFATH,ISO),100).EQ.1) GO TO 20
          JND=IPREAC(IFATH,ISO)/100
          DO J=1,NBCH
            IF(IPOS(J,1).EQ.JND) GO TO 20
          ENDDO
          IPREAC(IFATH,ISO)=0
          PRATE(IFATH,ISO)=0.0
   20     CONTINUE
        ENDDO
      ENDDO
*----
*  LUMP IDREAC, DENER, IPREAC AND PRATE.
*----
      ITER=0
   40 ITER=ITER+1
      IF(ITER.GT.MAXIT) CALL XABORT('EDILUM: TOO MANY ITERATIONS.')
      NLUMP=0
      DO ISO=1,NBCH
        IND=IPOS(ISO,1)
        DO IFATH=1,NFATH
          IF(IPREAC(IFATH,IND).EQ.0) GO TO 50
          IF(MOD(IPREAC(IFATH,IND),100).NE.1) GO TO 50
          JND=IPREAC(IFATH,IND)/100
          IF(MYLIST(JND).EQ.0) GO TO 50
          DO J=1,NBCH
            IF(IPOS(J,1).EQ.JND) GO TO 50
          ENDDO
          NLUMP=NLUMP+1 ! ISOTOPE JND IS LUMPED
          DO IDA=1,NBISO
            DO IFA=1,NFATH
              IPGAR=IPREAC(IFA,IDA)
              IF((IPGAR/100.EQ.JND).AND.(MOD(IPGAR,100).EQ.1)) THEN
                IF(MYLIST(IDA).EQ.0) GO TO 50
              ENDIF
            ENDDO
          ENDDO
          DO IDA=1,NBISO
            DO IFA=1,NFATH
              IPGAR=IPREAC(IFA,IDA)
              IF((IPGAR/100.EQ.JND).AND.(MOD(IPGAR,100).EQ.1)) THEN
                IF(IDA.EQ.JND) CALL XABORT('EDILUM: BUG.')
                PRGAR=PRATE(IFA,IDA)
                DO IM=IFA,NFATH-1
                  IPREAC(IM,IDA)=IPREAC(IM+1,IDA)
                  PRATE(IM,IDA)=PRATE(IM+1,IDA)
                ENDDO
                IPREAC(NFATH,IDA)=0
                PRATE(NFATH,IDA)=0.0
                DO JFATH=1,NFATH
                  IF(IPREAC(JFATH,JND).EQ.0) GO TO 45
                  IM=NFATH+1
                  DO K=NFATH,1,-1
                    IF(IPREAC(K,IDA).EQ.IPREAC(JFATH,JND)) THEN
                      PRATE(K,IDA)=PRATE(K,IDA)+PRGAR*PRATE(JFATH,JND)
                      GO TO 44
                    ENDIF
                    IF(IPREAC(K,IDA).EQ.0) IM=K
                  ENDDO
                  IF(IM.GT.NFATH) CALL XABORT('EDILUM: NFATH OVERFLOW.')
                  IPREAC(IM,IDA)=IPREAC(JFATH,JND)
                  PRATE(IM,IDA)=PRGAR*PRATE(JFATH,JND)
   44             CONTINUE
                ENDDO
   45           IF(MOD(IDREAC(2,JND),100).EQ.5) THEN
                  JFP=IDREAC(2,JND)/100
                  IF(MOD(IDREAC(2,IDA),100).EQ.5) THEN
                    IFP=IDREAC(2,IDA)/100
                  ELSE
                    NBDPF=NBDPF+1
                    IF(NBDPF.GT.MAXFP) THEN
                      CALL XABORT('EDILUM: MAXFP OVERFLOW.')
                    ENDIF
                    IFP=NBDPF
                    CALL XDRSET(YIELD(1,IFP),NBFISS,0.0)
                  ENDIF
                  DO IFI=1,NBFISS
                    YIELD(IFI,IFP)=YIELD(IFI,IFP)+YIELD(IFI,JFP)*PRGAR
                  ENDDO
                  IDREAC(2,IDA)=IFP*100+5
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          DO JFATH=1,NFATH
            IF(IPREAC(JFATH,JND).GT.0) THEN
              KT=MOD(IPREAC(JFATH,JND),100)
              KND=IPREAC(JFATH,JND)/100
              DENER(KT,KND)=DENER(KT,KND)+PRATE(JFATH,JND)*DENER(1,JND)
            ENDIF
            IPREAC(JFATH,JND)=0
            PRATE(JFATH,JND)=0.0
          ENDDO
          YMAX=0.0
          IF(MOD(IDREAC(2,JND),100).EQ.5) THEN
            JFP=IDREAC(2,JND)/100
            DO KSO=1,NBISO
              IF(MOD(IDREAC(2,KSO),100).EQ.4) THEN
                IFI=IDREAC(2,KSO)/100
                DENER(2,KSO)=DENER(2,KSO)+YIELD(IFI,JFP)*DENER(1,JND)
              ENDIF
            ENDDO
            DO IFI=1,NBFISS
              YMAX=MAX(YMAX,ABS(YIELD(IFI,JFP)))
              YIELD(IFI,JFP)=0.0
            ENDDO
            DENER(2,JND)=0.0
            IDREAC(2,JND)=0
          ENDIF
          DENER(1,JND)=0.0
          IDREAC(1,JND)=0
          HALF=1.0E8*LOG(2.0)/DDECA(JND)/86400.0
          IF(DDECA(JND).EQ.0.0) THEN
            WRITE(HNAME,'(2A4)') HISO(1,JND),HISO(2,JND)
            IF(IPRINT.GT.2) THEN
              WRITE(NSYSO,'('' WARNING: ISOTOPE '',A8,'' IS LUMPED A'',
     &        ''ND IS STABLE. MAX FISSION YIELD='',1P,E8.1,''%'')') 
     &        HNAME,YMAX*100.0
            ENDIF
            IF(YMAX.GT.1.0E-2) CALL XABORT('EDILUM: ISOTOPE '//HNAME//
     &      ' SHOULD NOT BE LUMPED.')
          ELSE IF((HALF.GT.30.0).AND.(HALF.LT.999999.99)) THEN
            WRITE(HNAME,'(2A4)') HISO(1,JND),HISO(2,JND)
            IF(IPRINT.GT.2) THEN
              WRITE(NSYSO,'('' WARNING: ISOTOPE '',A8,'' IS LUMPED A'',
     &        ''ND HAS A HALF-LIFE OF'',F10.2,'' DAYS. MAX FISSION Y'',
     &        ''IELD='',1P,E8.1,''%'')') HNAME,HALF,YMAX*100.0
            ENDIF
            IF(YMAX.GT.1.0E-2) CALL XABORT('EDILUM: ISOTOPE '//HNAME//
     &      ' SHOULD NOT BE LUMPED.')
          ELSE IF(HALF.GT.30.0) THEN
            WRITE(HNAME,'(2A4)') HISO(1,JND),HISO(2,JND)
            IF(IPRINT.GT.2) THEN
              WRITE(NSYSO,'('' WARNING: ISOTOPE '',A8,'' IS LUMPED A'',
     &        ''ND HAS A HALF-LIFE OF'',1P,E10.3,'' DAYS. MAX FISSIO'',
     &        ''N YIELD='',E8.1,''%'')') HNAME,HALF,YMAX*100.0
            ENDIF
            IF(YMAX.GT.1.0E-2) CALL XABORT('EDILUM: ISOTOPE '//HNAME//
     &      ' SHOULD NOT BE LUMPED.')
          ENDIF
          DDECA(JND)=0.0
          MYLIST(JND)=0
   50     CONTINUE
        ENDDO
      ENDDO
      IF(IPRINT.GT.2) WRITE(NSYSO,'('' ......... NLUMP='',I5)') NLUMP
      IF(NLUMP.GT.0) GO TO 40
*----
*  WRITE VECTORS 'PRODUCE-REAC' AND 'PRODUCE-RATE' TO THE LCM OBJECT
*----
      DO ISO=1,NBCH
        DO IFATH=1,MAXFAT
          JPREAC(IFATH,ISO)=0
          RRATE(IFATH,ISO)=0.0
        ENDDO
        IND=IPOS(ISO,1)
        NN=0
        DO IFATH=1,NFATH
          IF(IPREAC(IFATH,IND).NE.0) THEN
            DO J=1,IFATH-1
              IF(IPREAC(IFATH,IND).EQ.IPREAC(J,IND)) THEN
                JND1=IPREAC(IFATH,IND)/100
                JND2=IPREAC(J,IND)/100
                WRITE(NSYSO,'(/27H EDILUM: DUPLICATE FATHERS:,2A4,
     &          1X,2A4)') HISO(1,JND1),HISO(2,JND1),HISO(1,JND2),
     &          HISO(2,JND2)
                WRITE(HNAME,'(2A4)') HISO(1,IND),HISO(2,IND)
                CALL XABORT('EDILUM: DUPLICATE FATHERS FOR '//HNAME)
              ENDIF
            ENDDO
            DO I=1,NBCH
              JSO=I
              IF((IPOS(I,1).EQ.IPREAC(IFATH,IND)/100).AND.
     &            (MIX(I).EQ.MIX(ISO)))GO TO 70
            ENDDO
            JSO=-1
   70       IF(JSO.EQ.-1) THEN
              JND=IPREAC(IFATH,IND)/100
              IF(IPRINT.GT.2) THEN
                WRITE(NSYSO,'(/24H EDILUM: UNKNOWN FATHER ,2A4,5H FOR ,
     &          2A4)') HISO(1,JND),HISO(2,JND),HISO(1,IND),HISO(2,IND)
              ENDIF
            ELSE
              NN=NN+1
              IF(NN.GT.MAXFAT) THEN
                WRITE(TEXT4,'(I4)') NN
                CALL XABORT('EDILUM: MAXFAT OVERFLOW NN='//TEXT4)
              ENDIF
              JPREAC(NN,ISO)=100*JSO+MOD(IPREAC(IFATH,IND),100)
              RRATE(NN,ISO)=PRATE(IFATH,IND)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      CALL LCMPUT(IPEDIT,'PRODUCE-REAC',MAXFAT*NBCH,1,JPREAC)
      CALL LCMPUT(IPEDIT,'PRODUCE-RATE',MAXFAT*NBCH,2,RRATE)
*----
*  WRITE THE ISOTOPE ASCII NAMES ON LCM OBJECT
*----
      CALL LCMPUT(IPEDIT,'DEPLETE-IDEN',2*NREAC,3,HREAC(1,1))
*----
*  WRITE THE LUMPED FISSION YIELD MATRIX TO THE LCM OBJECT
*----
      IBFI=0
      IBFP=0
      DO ISO=1,NBCH
        IND=IPOS(ISO,1)
        IPOS(ISO,2)=0
        IPOS(ISO,3)=0
        IF(MOD(IDREAC(2,IND),100).EQ.4) THEN
          IBFI=IBFI+1
          IF(IBFI.GT.NBFISS2) CALL XABORT('EDILUM: NBFISS OVERFLOW.')
          IPOS(ISO,2)=IBFI
        ELSE IF(MOD(IDREAC(2,IND),100).EQ.5) THEN
          IBFP=IBFP+1
          IF(IBFP.GT.NBFPCH2) CALL XABORT('EDILUM: NBFPCH2 OVERFLOW.')
          IPOS(ISO,3)=IBFP
        ENDIF
      ENDDO
      DO ISO=1,NBCH
        IFI=IPOS(ISO,2)
        IF(IFI.GT.0) THEN
          DO JSO=1,NBCH
            IFP=IPOS(JSO,3)
            IF(IFP.GT.0) THEN
              IND=IPOS(ISO,1)
              JND=IPOS(JSO,1)
              IF(MIX(ISO).NE.MIX(JSO)) THEN
                EYIEL((IFP-1)*IBFI+IFI)=0.0
              ELSE
                EYIEL((IFP-1)*IBFI+IFI)=YIELD(IDREAC(2,IND)/100,
     &          IDREAC(2,JND)/100)
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      IF(IBFI*IBFP.GT.0) CALL LCMPUT(IPEDIT,'FISSIONYIELD',IBFI*IBFP,2,
     & EYIEL)
*----
*  WRITE VECTORS 'DEPLETE-REAC' AND 'DEPLETE-ENER' TO THE LCM OBJECT
*----
      DO ISO=1,NBCH
        IND=IPOS(ISO,1)
        DO I=1,NREAC
          IF(IDREAC(I,IND)/100.GT.0) THEN
            KREAC=MOD(IDREAC(I,IND),100)
            IF((KREAC.LE.0).OR.(KREAC.GT.5)) THEN
              CALL XABORT('EDILUM: INVALID REACTION.')
            ENDIF
          ENDIF
          IF((I.EQ.2).AND.(MOD(IDREAC(I,IND),100).EQ.4)) THEN
            IF(IPOS(ISO,2).EQ.0) CALL XABORT('EDILUM: FAILURE 1.')
            JDREAC(I,ISO)=IPOS(ISO,2)*100+4
            IF(IPRINT.GT.2) THEN
              WRITE(NSYSO,6010) (HISO(I0,IPOS(ISO,1)),I0=1,3)
            ENDIF
          ELSE IF((I.EQ.2).AND.(MOD(IDREAC(I,IND),100).EQ.5)) THEN
            IF(IPOS(ISO,3).EQ.0) CALL XABORT('EDILUM: FAILURE 2.')
            JDREAC(I,ISO)=IPOS(ISO,3)*100+5
          ELSE
            JDREAC(I,ISO)=IDREAC(I,IND)
          ENDIF
          EENER(I,ISO)=DENER(I,IND)
        ENDDO
      ENDDO
      CALL LCMPUT(IPEDIT,'DEPLETE-REAC',NREAC*NBCH,1,JDREAC)
      CALL LCMPUT(IPEDIT,'DEPLETE-ENER',NREAC*NBCH,2,EENER)
*----
*  WRITE VECTORS 'CHARGEWEIGHT', 'DEPLETE-DECA', 'ISOTOPESDEPL' AND
*  'STATE-VECTOR' TO THE LCM OBJECT
*----
      NBHEAV=0
      NSUPS2=0
      DO ISO=1,NBCH
        IF(IPOS(ISO,1).GT.NBISO-NSUPS) NSUPS2=NSUPS2+1
        HHHH(1,ISO)=HISO(1,IPOS(ISO,1))
        HHHH(2,ISO)=HISO(2,IPOS(ISO,1))
        HHHH(3,ISO)=HISO(3,IPOS(ISO,1))
        DDDD(ISO)=DDECA(IPOS(ISO,1))
        IPOS(ISO,1)=MYLIST(IPOS(ISO,1))
      ENDDO
      CALL LCMPUT(IPEDIT,'ISOTOPESDEPL',3*NBCH,3,HHHH)
      CALL LCMPUT(IPEDIT,'CHARGEWEIGHT',NBCH,1,IPOS(1,1))
      CALL LCMPUT(IPEDIT,'DEPLETE-DECA',NBCH,2,DDDD)
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NBCH
      ISTATE(2)=IBFI
      ISTATE(3)=IBFP
      ISTATE(4)=NBHEAV
      ISTATE(5)=NBCH-NBHEAV
      ISTATE(7)=NSUPS2
      ISTATE(8)=NREAC
      ISTATE(9)=MAXFAT
      CALL LCMPUT(IPEDIT,'STATE-VECTOR',NSTATE,1,ISTATE)
      IF(IPRINT.GT.2) WRITE(NSYSO,6000) (ISTATE(ISTA),ISTA=1,9)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DDDD,EYIEL,EENER,RRATE)
      DEALLOCATE(HHHH,IPOS,JDREAC,JPREAC)
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(/' STATE-VECTOR FOR LUMPED DEPLETION CHAIN'/' -------'/
     >  ' NDEPL  ',I6,' (NUMBER OF DEPLETING ISOTOPES)'/
     >  ' NDFI   ',I6,' (NUMBER OF DIRECT FISSILE ISOTOPES)'/
     >  ' NDFP   ',I6,' (NUMBER OF DIRECT FISSION PRODUCT)'/
     >  ' NHEAVY ',I6,' (NUMBER OF HEAVY ISOTOPES)'/
     >  ' NLIGHT ',I6,' (NUMBER OF FISSION PRODUCTS)'/
     >  ' NOTHER ',I6,' (NUMBER OF OTHER ISOTOPES)'/
     >  ' NSTABL ',I6,' (NUMBER OF STABLE ISOTOPES PRODUCING ENERGY)'/
     >  ' NREAC  ',I6,' (MAXIMUM NUMBER OF DEPLETION REACTIONS)'/
     >  ' NPAR   ',I6,' (MAXIMUM NUMBER OF PARENT REACTIONS)'/)
 6010 FORMAT(45H EDILUM: FISSILE ISOTOPE WITH FISSION YIELD: ,3A4)
      END
