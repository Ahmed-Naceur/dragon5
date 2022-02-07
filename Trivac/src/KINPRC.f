*DECK KINPRC
      SUBROUTINE KINPRC(IPTRK,IPSYS,CMOD,NGR,NBM,NBFIS,NDG,NEL,LL4,NUN,
     1 NUP,MAT,VOL,IDLPC,FLN,FLO,CHD,CHO,SGD,SGO,PDC,DT,ADJ,TTP,PC,IPR,
     2 IEXP,OMEGA,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the precursors unknowns for the current time step according
* to the pre-defined temporal integration scheme.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal.
*
*Author(s): D. Sekki
*
*Parameters: input/output
* IPTRK  pointer to L_TRACK object.
* IPSYS  pointer to L_SYSTEM object.
* CMOD   name of the assembly door (BIVAC or TRIVAC).
* NGR    number of energy groups.
* NBM    number of material mixtures.
* NBFIS  number of fissile isotopes.
* NDG    number of delayed-neutron groups.
* NEL    total number of finite elements.
* LL4    number of flux unknowns per energy group.
* NUN    total number of unknowns per energy group.
* NUP    number of precursor unknowns per delayed group.
* MAT    mixture index assigned to each volume.
* VOL    volume of each element.
* IDLPC  position of averaged precursor values in unknown vector.
* FLN    unknown flux vector at current time step.
* FLO    unknown flux vector at previous time step.
* CHD    current delayed fission spectrum.
* CHO    previous delayed fission spectrum.
* SGD    current delayed nu*fission macroscopic x-sections/keff.
* SGO    previous delayed nu*fission macroscopic x-sections/keff.
* PDC    precursor decay constants.
* DT     current time increment.
* ADJ    flag for adjoint space-time kinetics calculation.
* TTP    value of theta-parameter for precursors.
* PC     unknown vector for precursors.
* IPR    integration scheme for precursors: =1 implicit;
*        =2 Crank-Nicholson; =3 theta; =4 exponential.
* IEXP   exponential transformation flag (=1 to activate).
* OMEGA  exponential transformation parameter.
* IMPX   printing parameter (=0 for no print).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS
      INTEGER NGR,NBM,NBFIS,NDG,NEL,LL4,NUN,NUP,MAT(NEL),IDLPC(NEL),
     1 IPR,IEXP,IMPX
      REAL VOL(NEL),PDC(NDG),DT,TTP,PC(NUP,NDG,NBFIS),FLN(NUN,NGR),
     1 FLO(NUN,NGR),CHD(NBM,NBFIS,NGR,NDG),CHO(NBM,NBFIS,NGR,NDG),
     2 SGD(NBM,NBFIS,NGR,NDG),SGO(NBM,NBFIS,NGR,NDG),OMEGA(NBM,NGR)
      CHARACTER CMOD*12
      LOGICAL ADJ
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOS=6)
      DOUBLE PRECISION DK,DTP,TK1(NDG),TK2(NDG),TK3(NDG)
      REAL, DIMENSION(:), ALLOCATABLE :: GAR1,GAR2
      REAL, DIMENSION(:,:), ALLOCATABLE :: XSEXP
      REAL, DIMENSION(:), POINTER :: RM
      TYPE(C_PTR) RM_PTR
*----
*  COMPUTE THE KINETICS FACTORS
*----
      CALL XDDSET(TK1,NDG,0.0D0)
      CALL XDDSET(TK2,NDG,0.0D0)
      CALL XDDSET(TK3,NDG,0.0D0)
      DTP=9999.0D0
      IF(IPR.EQ.2)THEN
*     CRANK-NICHOLSON
        DTP=0.5D0
      ELSEIF(IPR.EQ.3)THEN
*     THETA
        DTP=DBLE(TTP)
      ENDIF
      DO 10 L=1,NDG
      DK=PDC(L)*DT
      IF(IPR.EQ.1)THEN
*     IMPLICIT
        TK1(L)=1.0D0/(1.0D0+DK)
        TK2(L)=DT/(1.0D0+DK)
      ELSEIF(IPR.EQ.4)THEN
*     EXPONENTIAL
        TK1(L)=DEXP(-DK)
        TK2(L)=(1.0D0-(1.0D0-TK1(L))/DK)/PDC(L)
        TK3(L)=((1.0D0-TK1(L))/DK-TK1(L))/PDC(L)
      ELSE
*     GENERAL
        TK1(L)=(1.0D0-(1.0D0-DTP)*DK)/(1.0D0+DTP*DK)
        TK2(L)=DTP*DT/(1.0D0+DTP*DK)
        TK3(L)=(1.0D0-DTP)*DT/(1.0D0+DTP*DK)
      ENDIF
   10 CONTINUE
*----
*  COMPUTE THE PRECURSOR UNKNOWN VECTOR
*----
      IF(IMPX.GT.0)WRITE(IOS,1001)CMOD
      ALLOCATE(GAR1(NUP),GAR2(NUP),XSEXP(NBM,NGR))
      DO 220 IFIS=1,NBFIS
      DO 210 IDG=1,NDG
      DO 20 I=1,NUP
      PC(I,IDG,IFIS)=REAL(TK1(IDG))*PC(I,IDG,IFIS)
   20 CONTINUE
      CALL XDRSET(GAR2,NUP,0.0)
      IF(.NOT.ADJ) THEN
        DO 40 IGR=1,NGR
        DO 30 IBM=1,NBM
        IF(IEXP.EQ.0) THEN
          XSEXP(IBM,IGR)=SGD(IBM,IFIS,IGR,IDG)
        ELSE
*         exponential transformation
          XSEXP(IBM,IGR)=SGD(IBM,IFIS,IGR,IDG)*EXP(OMEGA(IBM,IGR)*DT)
        ENDIF
   30   CONTINUE
   40   CONTINUE
      ELSE
        DO 60 IGR=1,NGR
        DO 50 IBM=1,NBM
        IF(IEXP.EQ.0) THEN
          XSEXP(IBM,IGR)=PDC(IDG)*CHD(IBM,IFIS,IGR,IDG)
        ELSE
*         exponential transformation
          XSEXP(IBM,IGR)=PDC(IDG)*CHD(IBM,IFIS,IGR,IDG)*
     1    EXP(OMEGA(IBM,IGR)*DT)
        ENDIF
   50   CONTINUE
   60   CONTINUE
      ENDIF
      IF(CMOD.EQ.'BIVAC')THEN
        ITY=1
        DO 80 IGR=1,NGR
        CALL KINBLM(IPTRK,NBM,NUP,XSEXP(1,IGR),FLN(1,IGR),GAR1)
        DO 70 IND=1,NUP
        GAR2(IND)=GAR2(IND)+GAR1(IND)
   70   CONTINUE
   80   CONTINUE
        CALL MTLDLS('RM',IPTRK,IPSYS,LL4,1,GAR2)
      ELSEIF(CMOD.EQ.'TRIVAC')THEN
        DO 100 IGR=1,NGR
        CALL KINTLM(IPTRK,NBM,NUP,XSEXP(1,IGR),FLN(1,IGR),GAR1)
        DO 90 IND=1,NUP
        GAR2(IND)=GAR2(IND)+GAR1(IND)
   90   CONTINUE
  100   CONTINUE
        CALL LCMLEN(IPSYS,'RM',ILONG,ITYLCM)
        CALL LCMGPD(IPSYS,'RM',RM_PTR)
        CALL C_F_POINTER(RM_PTR,RM,(/ ILONG /))
        DO 110 IND=1,ILONG
        GAR2(IND)=GAR2(IND)/RM(IND)
  110   CONTINUE
      ENDIF
      DO 120 IND=1,NUP
      PC(IND,IDG,IFIS)=PC(IND,IDG,IFIS)+REAL(TK2(IDG))*GAR2(IND)
  120 CONTINUE
      IF(IPR.GT.1) THEN
        CALL XDRSET(GAR2,NUP,0.0)
        IF(CMOD.EQ.'BIVAC')THEN
          IF(ADJ) CALL XABORT('KINPRC: ADJOINT CALCULATION NOT IMPLEME'
     1    //'NTED WITH BIVAC.')
          ITY=1
          DO 140 IGR=1,NGR
          CALL KINBLM(IPTRK,NBM,NUP,SGO(1,IFIS,IGR,IDG),FLO(1,IGR),
     1    GAR1)
          DO 130 IND=1,NUP
          GAR2(IND)=GAR2(IND)+GAR1(IND)
  130     CONTINUE
  140     CONTINUE
          CALL MTLDLS('RM',IPTRK,IPSYS,LL4,1,GAR2)
          CALL FLDBIV(IPTRK,NEL,NUP,GAR2,MAT,VOL,IDLPC)
        ELSEIF(CMOD.EQ.'TRIVAC')THEN
          IF(.NOT.ADJ) THEN
            DO 160 IGR=1,NGR
            CALL KINTLM(IPTRK,NBM,NUP,SGO(1,IFIS,IGR,IDG),FLO(1,IGR),
     1      GAR1)
            DO 150 IND=1,NUP
            GAR2(IND)=GAR2(IND)+GAR1(IND)
  150       CONTINUE
  160       CONTINUE
          ELSE
            DO 180 IGR=1,NGR
            CALL KINTLM(IPTRK,NBM,NUP,CHO(1,IFIS,IGR,IDG),FLO(1,IGR),
     1      GAR1)
            DO 170 IND=1,NUP
            GAR2(IND)=GAR2(IND)+PDC(IDG)*GAR1(IND)
  170       CONTINUE
  180       CONTINUE
          ENDIF
          CALL LCMLEN(IPSYS,'RM',ILONG,ITYLCM)
          CALL LCMGPD(IPSYS,'RM',RM_PTR)
          CALL C_F_POINTER(RM_PTR,RM,(/ ILONG /))
          DO 190 IND=1,ILONG
          GAR2(IND)=GAR2(IND)/RM(IND)
  190     CONTINUE
        ENDIF
        DO 200 IND=1,NUP
        PC(IND,IDG,IFIS)=PC(IND,IDG,IFIS)+REAL(TK3(IDG))*GAR2(IND)
  200   CONTINUE
      ENDIF
      IF(CMOD.EQ.'BIVAC')THEN
        CALL FLDBIV(IPTRK,NEL,NUP,PC(1,IDG,IFIS),MAT,VOL,IDLPC)
      ELSEIF(CMOD.EQ.'TRIVAC')THEN
        CALL FLDTRI(IPTRK,NEL,NUP,PC(1,IDG,IFIS),MAT,VOL,IDLPC)
      ENDIF
  210 CONTINUE
  220 CONTINUE
      DEALLOCATE(XSEXP,GAR1,GAR2)
*----
*  EDITION
*----
      IF(IMPX.GT.5) THEN
        WRITE(IOS,1002)
        DO 240 IFIS=1,NBFIS
        DO 230 IDG=1,NDG
        WRITE(IOS,1003) IDG,IFIS,(PC(IND,IDG,IFIS),IND=1,LL4)
  230   CONTINUE
  240   CONTINUE
      ENDIF
      IF(IMPX.GT.2) THEN
        DO 260 IFIS=1,NBFIS
        WRITE(IOS,1004) IFIS,(IDG,IDG=1,NDG)
        DO 250 IEL=1,NEL
        IND=IDLPC(IEL)
        IF(IND.EQ.0) GO TO 250
        WRITE(IOS,1005) IEL,(PC(IND,IDG,IFIS),IDG=1,NDG)
  250   CONTINUE
        WRITE(IOS,'(/)')
  260   CONTINUE
      ENDIF
      RETURN
*
 1001 FORMAT(/1X,'COMPUTING THE PRECURSOR UNKNOWN VECTOR',
     1 1X,'ACCORDING TO THE TRACKING TYPE: ',A6/)
 1002 FORMAT(/1X,'=> COMPUTED PRECURSOR UNKNOWN VECTOR')
 1003 FORMAT(/17H PRECURSOR GROUP=,I5,18H  FISSILE ISOTOPE=,I5/
     1 (1P,8E14.5))
 1004 FORMAT(/51H KINPRC: PRECURSOR UNKNOWN VECTOR (FISSILE ISOTOPE=,
     1 I5,1H)/(9X,6I13,:))
 1005 FORMAT(1X,I6,2X,1P,6E13.5,:/(9X,6E13.5,:))
      END
