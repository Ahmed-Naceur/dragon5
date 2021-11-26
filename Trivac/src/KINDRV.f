*DECK KINDRV
      SUBROUTINE KINDRV(NEN,KEN,CMOD,NGR,NBM,NBFIS,NDG,NLF,ITY,NEL,
     1 LL4,NUN,NUP,TTF,TTP,DT,IMPH,ICL1,ICL2,NADI,ADJ,MAXOUT,EPSOUT,
     2 MAXINR,EPSINR,IFL,IPR,IEXP,INORM,IMPX,POWTOT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver to perform the space-time kinetics calculations.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal.
*
*Author(s): D. Sekki
*
*Parameters: input
* NEN     number of LCM objects used in the module.
* KEN     addresses of LCM objects: (1) L_KINET; (2) L_MACROLIB;
*         (3) L_TRACK; (4) L_SYSTEM; (5) L_MACROLIB.
* CMOD    name of the assembly door (BIVAC or TRIVAC).
* NGR     number of energy groups.
* NBM     number of material mixtures.
* NBFIS   number of fissile isotopes.
* NDG     number of delayed-neutron groups.
* NLF     number of Legendre orders for fluxes.
* ITY     type of finite elements and tracking.
* NEL     total number of finite elements.
* LL4     number of flux unknowns per energy group.
* NUN     total number of unknowns per energy group.
* NUP     number of precursor unknowns per delayed group.
* TTF     value of theta-parameter for fluxes.
* TTP     value of theta-parameter for precursors.
* DT      current time increment.
* IMPH    management of convergence histogram.
* ICL1    number of free iterations in one cycle of the inverse power
*         method
* ICL2    number of accelerated iterations in one cycle
* NADI    number of inner adi iterations per outer iteration
* ADJ     flag for adjoint space-time kinetics calculation
* MAXOUT  maximum number of outer iterations
* EPSOUT  convergence criteria for the flux
* MAXINR  maximum number of thermal iterations.
* EPSINR  thermal iteration epsilon.
* IFL     temporal integration scheme for fluxes.
* IPR     temporal integration scheme for precursors.
* IEXP    exponential transformation flag (=1 to activate).
* INORM   type of flux normalization (=0: no normalization; =1: imposed
*         factor; =2: maximum flux; =3 initial power).
* IMPX    printing parameter (=0 for no print).
*
*Parameter: output
* POWTOT  power.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NEN,NGR,NBM,NBFIS,NDG,NLF,ITY,NEL,LL4,NUN,NUP,IMPH,ICL1,
     1 ICL2,NADI,MAXOUT,MAXINR,IFL,IPR,IEXP,INORM,IMPX
      TYPE(C_PTR) KEN(NEN)
      REAL TTF,TTP,DT,EPSOUT,EPSINR,POWTOT
      CHARACTER CMOD*12
      LOGICAL ADJ
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOS=6)
      INTEGER MAT(NEL),IDL(NEL),IDLPC(NEL)
      REAL VOL(NEL),PDC(NDG),PMAX(NDG,NBFIS)
      LOGICAL LNUD,LCHD
      TYPE(C_PTR) IPMAC,IPSYS
      REAL, DIMENSION(:), ALLOCATABLE :: DNF,AVG1,AVG2,WORK1,RM
      REAL, DIMENSION(:,:), ALLOCATABLE :: EVECT,DNS,PHO,OVR,OMEGA
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: PC,CHI,SGF
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: SGD,CHD,SGO
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: SRC
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(EVECT(NUN,NGR),PC(NUP,NDG,NBFIS),SGD(NBM,NBFIS,NGR,NDG),
     1 OMEGA(NBM,NGR))
*----
*  RECOVER INFORMATION
*----
      CALL KDRCPU(TK1)
      TA1=TK1
      IF(IMPX.GT.0) WRITE(IOS,1001)
      CALL LCMGET(KEN(1),'E-KEFF',EVL)
      CALL LCMGET(KEN(1),'LAMBDA-D',PDC)
      CALL LCMGET(KEN(1),'E-IDLPC',IDLPC)
      CALL LCMLEN(KEN(1),'OMEGA',ILONG,ITYLCM)
      IF((IEXP.EQ.0).OR.(ILONG.EQ.0)) THEN
        CALL XDRSET(OMEGA,NBM*NGR,0.0)
      ELSE
        CALL LCMGET(KEN(1),'OMEGA',OMEGA)
      ENDIF
      CALL LCMGET(KEN(3),'VOLUME',VOL)
      CALL LCMGET(KEN(3),'MATCOD',MAT)
      CALL LCMGET(KEN(3),'KEYFLX',IDL)
*----
*  RECOVER CROSS SECTIONS (BEGINNING-OF-STEP)
*----
      ALLOCATE(DNF(NDG),DNS(NGR,NDG))
      CALL LCMLEN(KEN(1),'BETA-D',LEN,ITYL)
      LNUD=(LEN.EQ.NDG)
      IF(LNUD) CALL LCMGET(KEN(1),'BETA-D',DNF)
      CALL LCMLEN(KEN(1),'CHI-D',LEN,ITYL)
      LCHD=(LEN.EQ.NGR*NDG)
      IF(LCHD) CALL LCMGET(KEN(1),'CHI-D',DNS)
      ALLOCATE(OVR(NBM,NGR),CHI(NBM,NBFIS,NGR),CHD(NBM,NBFIS,NGR,NDG),
     1 SGF(NBM,NBFIS,NGR),SGO(NBM,NBFIS,NGR,NDG))
      IF(NEN.EQ.4) THEN
        IPMAC=KEN(2)
        IPSYS=KEN(4)
      ELSE IF(NEN.EQ.6) THEN
        IPMAC=KEN(5)
        IPSYS=KEN(6)
      ENDIF
      CALL KINXSD(IPMAC,NGR,NBM,NBFIS,NDG,EVL,DT,DNF,DNS,LNUD,LCHD,OVR,
     1 CHI,CHD,SGF,SGO)
*----
*  COMPUTE THE SOURCE TERM
*----
      LL4B=LL4
      IF((ITY.EQ.11).OR.(ITY.EQ.13)) LL4B=LL4*NLF/2
      ALLOCATE(PHO(NUN,NGR),SRC(NUN,NGR))
      CALL LCMGET(KEN(1),'E-PREC',PC)
      CALL LCMGET(KEN(1),'E-VECTOR',PHO)
      CALL KINSRC(KEN(3),IPSYS,CMOD,IMPX,IFL,IPR,IEXP,NGR,NBM,NBFIS,NDG,
     1 ITY,LL4B,NUN,NUP,PDC,TTF,TTP,DT,ADJ,OVR,CHI,CHD,SGF,SGO,OMEGA,
     2 PHO,PC,SRC)
*----
*  RECOVER CROSS SECTIONS (END-OF-STEP)
*----
      CALL KINXSD(KEN(2),NGR,NBM,NBFIS,NDG,EVL,DT,DNF,DNS,LNUD,LCHD,
     1 OVR,CHI,CHD,SGF,SGD)
      DEALLOCATE(DNS,DNF)
*----
*  RECOVER THE BEGINNING-OF-STEP FLUX
*----
      IF(IMPX.GT.0)THEN
        CALL KDRCPU(TA2)
        WRITE(IOS,1002) TA2-TA1
        WRITE(IOS,1003)
      ENDIF
      DO 15 IGR=1,NGR
      DO 10 IND=1,NUN
      EVECT(IND,IGR)=PHO(IND,IGR)
   10 CONTINUE
   15 CONTINUE
*----
*  COMPUTE THE FLUX SOLUTION
*----
      IF(CMOD.EQ.'BIVAC')THEN
        IF(ADJ) CALL XABORT('KINDRV: ADJOINT CALCULATION NOT IMPLEMENT'
     1  //'ED WITH BIVAC.')
        CALL KINSLB(KEN(3),KEN(4),KEN(1),LL4B,ITY,NUN,NGR,IFL,IPR,
     1  IEXP,NBM,NBFIS,NDG,ICL1,ICL2,IMPX,IMPH,TITR,EPSOUT,MAXINR,
     2  EPSINR,MAXOUT,PDC,TTF,TTP,DT,OVR,CHI,CHD,SGF,SGD,OMEGA,EVECT,
     3  SRC)
      ELSEIF(CMOD.EQ.'TRIVAC')THEN
        CALL KINSLT(KEN(3),KEN(4),KEN(1),LL4B,ITY,NUN,NGR,IFL,IPR,
     1  IEXP,NBM,NBFIS,NDG,ICL1,ICL2,IMPX,IMPH,TITR,EPSOUT,MAXINR,
     2  EPSINR,NADI,ADJ,MAXOUT,PDC,TTF,TTP,DT,OVR,CHI,CHD,SGF,SGD,
     3  OMEGA,EVECT,SRC)
      ENDIF
      DEALLOCATE(SRC)
*----
*  COMPUTE THE PRECURSOR SOLUTION
*----
      CALL KINPRC(KEN(3),KEN(4),CMOD,NGR,NBM,NBFIS,NDG,NEL,LL4,NUN,NUP,
     1 MAT,VOL,IDLPC,EVECT,PHO,CHD,CHO,SGD,SGO,PDC,DT,ADJ,TTP,PC,IPR,
     2 IEXP,OMEGA,IMPX)
      CALL LCMPUT(KEN(1),'E-PREC',NDG*NUP*NBFIS,2,PC)
*----
*  COMPUTE THE EXPONENTIAL TRANSFORMATION FACTORS
*----
      IF(IEXP.EQ.1) THEN
        ALLOCATE(WORK1(LL4),RM(LL4),AVG1(NBM),AVG2(NBM))
        DO 35 IGR=1,NGR
        DO 20 IBM=1,NBM
        AVG1(IBM)=EXP(OMEGA(IBM,IGR)*DT)
   20   CONTINUE
        IF(CMOD.EQ.'BIVAC')THEN
          CALL KINBLM(KEN(3),NBM,LL4,AVG1,EVECT(1,IGR),WORK1)
          CALL MTLDLS('RM',KEN(3),KEN(4),LL4,1,WORK1)
        ELSEIF(CMOD.EQ.'TRIVAC')THEN
          CALL KINTLM(KEN(3),NBM,LL4,AVG1,EVECT(1,IGR),WORK1)
          CALL LCMLEN(KEN(4),'RM',ILONG,ITYLCM)
          CALL LCMGET(KEN(4),'RM',RM)
          DO 25 IND=1,ILONG
          FACT=RM(IND)
          IF(FACT.EQ.0.0) CALL XABORT('KINDRV: SINGULAR RM.')
          WORK1(IND)=WORK1(IND)/FACT
   25     CONTINUE
        ENDIF
        DO 30 IND=1,LL4
        EVECT(IND,IGR)=WORK1(IND)
   30   CONTINUE
   35   CONTINUE
        CALL LCMPUT(KEN(1),'E-VECTOR',NGR*NUN,2,EVECT)
*
        DO 60 IGR=1,NGR
        CALL XDRSET(AVG1,NBM,0.0)
        CALL XDRSET(AVG2,NBM,0.0)
        DO 40 IEL=1,NEL
        IBM=MAT(IEL)
        IF(IBM.GT.0) THEN
          AVG1(IBM)=AVG1(IBM)+VOL(IEL)*PHO(IDL(IEL),IGR)
          AVG2(IBM)=AVG2(IBM)+VOL(IEL)*EVECT(IDL(IEL),IGR)
        ENDIF
   40   CONTINUE
        DO 50 IBM=1,NBM
        RATIO=MIN(10.0,ABS(AVG2(IBM)/AVG1(IBM)))
        OMEGA(IBM,IGR)=LOG(RATIO)/DT
   50   CONTINUE
        IF(IMPX.GT.1) THEN
          WRITE(IOS,1006) (OMEGA(IBM,IGR),IBM=1,NBM)
        ENDIF
   60   CONTINUE
        CALL LCMPUT(KEN(1),'OMEGA',NBM*NGR,2,OMEGA)
        DEALLOCATE(AVG2,AVG1,RM,WORK1)
      ENDIF
*----
*  COMPUTE AVERAGED FLUX VALUES.
*----
      DO 70 IGR=1,NGR
      IF(CMOD.EQ.'BIVAC')THEN
        CALL FLDBIV(KEN(3),NEL,NUP,EVECT(1,IGR),MAT,VOL,IDL)
      ELSEIF(CMOD.EQ.'TRIVAC')THEN
        CALL FLDTRI(KEN(3),NEL,NUP,EVECT(1,IGR),MAT,VOL,IDL)
      ENDIF
   70 CONTINUE
      CALL LCMPUT(KEN(1),'E-VECTOR',NGR*NUN,2,EVECT)
*----
*  FIND THE MAXIMUM FLUX VALUE
*----
      FMAX=0.0
      IDMX=0
      DO 85 IGR=1,NGR
      DO 80 IEL=1,NEL
      IND=IDL(IEL)
      IF(IND.EQ.0) GO TO 80
      IF(ABS(EVECT(IND,IGR)).GT.FMAX) THEN
        FMAX=EVECT(IND,IGR)
        IDMX=IEL
        IGMX=IGR
      ENDIF
   80 CONTINUE
   85 CONTINUE
      IF(IDMX.EQ.0) CALL XABORT('KINDRV: UNABLE TO SET FMAX.')
      IND=IDLPC(IDMX)
      IF(IND.EQ.0) CALL XABORT('KINDRV: UNABLE TO SET PMAX.')
      DO 95 IFIS=1,NBFIS
      DO 90 IDG=1,NDG
      PMAX(IDG,IFIS)=PC(IND,IDG,IFIS)
   90 CONTINUE
   95 CONTINUE
      IF(IMPX.GT.0) THEN
        WRITE(IOS,1004) FMAX,IDMX,IGMX
        CALL KDRCPU(TK2)
        WRITE(IOS,1005)TK2-TK1
      ENDIF
      CALL LCMPUT(KEN(1),'CTRL-FLUX',1,2,FMAX)
      CALL LCMPUT(KEN(1),'CTRL-PREC',NDG*NBFIS,2,PMAX)
      CALL LCMPUT(KEN(1),'CTRL-IDL',1,1,IDMX)
      CALL LCMPUT(KEN(1),'CTRL-IGR',1,1,IGMX)
*----
*  COMPUTE REACTOR POWER
*----
      IF(INORM.EQ.3) THEN
        CALL KINPOW(KEN(2),NGR,NBM,NUN,NEL,MAT,VOL,IDL,EVECT,POWTOT)
        CALL LCMPUT(KEN(1),'E-POW',1,2,POWTOT)
        IF(IMPX.GT.0) WRITE(6,*) 'REACTOR POWER (MW) =',POWTOT
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(PHO,SGO,SGF,CHD,CHI,OVR)
      DEALLOCATE(OMEGA,SGD,PC,EVECT)
      RETURN
*
 1001 FORMAT(/1X,'=> ASSEMBLY OF THE SYSTEM MATRICES'/)
 1002 FORMAT(/1X,'TOTAL CPU TIME USED FOR THE ASSEMBLING',
     1 1X,'OF ALL SYSTEM MATRICES =',F6.3/)
 1003 FORMAT(/1X,'=> COMPUTING THE KINETICS SOLUTION'/)
 1004 FORMAT(/1X,'CONTROLLING PARAMETERS:',2X,'MAX-VA',
     1 'L',1X,1PE12.5,3X,'IDL #',I5.5,3X,'IGR #',I2.2/)
 1005 FORMAT(/1X,'TOTAL CPU TIME USED FOR KINETICS CALC',
     1 'ULATIONS =',F10.3//1X,'=> SPACE-TIME',1X,
     2 'KINETICS CALCULATION IS DONE.')
 1006 FORMAT(39H KINDRV: MIXTURE-ORDERED OMEGA FACTORS:/(1P,10E14.6))
      END
