*DECK KINSLB
      SUBROUTINE KINSLB (IPTRK,IPSYS,IPKIN,LL4,ITY,NUN,NGR,IFL,IPR,IEXP,
     1 NBM,NBFIS,NDG,ICL1,ICL2,IMPX,IMPH,TITR,EPS2,MAXINR,EPSINR,MAXX0,
     2 PDC,TTF,TTP,DT,OVR,CHI,CHD,SGF,SGD,OMEGA,EVECT,SRC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of the kinetics multigroup linear systems for the transient
* neutron fluxes in Bivac. Use the inverse power method with a
* two-parameter SVAT acceleration technique.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPTRK   L_TRACK pointer to the tracking information.
* IPSYS   L_SYSTEM pointer to system matrices.
* IPKIN   L_KINET pointer to the KINET object.
* LL4     order of the system matrices.
* ITY     type of solution (1: classical Bivac/diffusion;
*         11: Bivac/SPN).
* NUN     number of unknowns in each energy group.
* NGR     number of energy groups.
* IFL     integration scheme for fluxes: =1 implicit;
*         =2 Crank-Nicholson; =3 theta.
* IPR     integration scheme for precursors: =1 implicit;
*         =2 Crank-Nicholson; =3 theta; =4 exponential.
* IEXP    exponential transformation flag (=1 to activate).
* NBM     number of material mixtures.
* NBFIS   number of fissile isotopes.
* NDG     number of delayed-neutron groups.
* ICL1    number of free iterations in one cycle of the inverse power
*         method
* ICL2    number of accelerated iterations in one cycle
* IMPX    print parameter. =0: no print ; =1: minimum printing ;
*         =2: iteration history is printed. =3: solution is printed
* IMPH    =0: no action is taken
*         =1: the flux is compared to a reference flux stored on lcm
*         =2: the convergence histogram is printed
*         =3: the convergence histogram is printed with axis and
*            titles. The plotting file is completed
*         =4: the convergence histogram is printed with axis, acce-
*            leration factors and titles. The plotting file is
*            completed
* TITR    character*72 title
* EPS2    convergence criteria for the flux
* MAXINR  maximum number of thermal iterations.
* EPSINR  thermal iteration epsilon.
* MAXX0   maximum number of outer iterations
* PDC     precursor decay constants.
* TTF     value of theta-parameter for fluxes.
* TTP     value of theta-parameter for precursors.
* DT      current time increment.
* OVR     reciprocal neutron velocities/DT.
* CHI     steady-state fission spectrum.
* CHD     delayed fission spectrum
* SGF     nu*fission macroscopic x-sections/keff.
* SGD     delayed nu*fission macroscopic x-sections/keff.
* OMEGA   exponential transformation parameter.
* SRC     fixed source
*
*Parameters: output
* EVECT    converged solution
*
*References:
* A. H\'ebert, 'Preconditioning the power method for reactor
* calculations', Nucl. Sci. Eng., 94, 1 (1986).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER TITR*72
      TYPE(C_PTR) IPTRK,IPSYS,IPKIN
      INTEGER LL4,ITY,NUN,NGR,IFL,IPR,IEXP,NBM,NBFIS,NDG,ICL1,ICL2,IMPX,
     1 IMPH,MAXINR,MAXX0
      REAL EPS2,EPSINR,PDC(NDG),TTF,TTP,DT,OVR(NBM,NGR),
     1 CHI(NBM,NBFIS,NGR),CHD(NBM,NBFIS,NGR,NDG),SGF(NBM,NBFIS,NGR),
     2 SGD(NBM,NBFIS,NGR,NDG),OMEGA(NBM,NGR),EVECT(NUN,NGR)
      DOUBLE PRECISION SRC(NUN,NGR)
*----
*  LOCAL VARIABLES
*----
      CHARACTER*12 TEXT12
      LOGICAL LOGTES,LMPH
      DOUBLE PRECISION D2F(2,3),ALP,BET,DTF,DTP,DARG,DK
      REAL ERR(250),ALPH(250),BETA(250),TKT,TKB
      INTEGER  ITITR(18)
      REAL, DIMENSION(:,:), ALLOCATABLE :: GRAD1,GRAD2
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: GAR1,GAR2,GAR3
      REAL, DIMENSION(:), ALLOCATABLE :: WORK1,WORK2,WORK3,WORK4
      DATA EPS1,MMAXX/1.0E-4,250/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(GRAD1(NUN,NGR),GRAD2(NUN,NGR),GAR1(NUN,NGR),
     1 GAR2(NUN,NGR),GAR3(NUN,NGR),WORK1(LL4),WORK2(LL4),WORK3(NBM))
*
      CALL MTOPEN(IMPX,IPTRK,LL4)
      IF(LL4.GT.NUN) CALL XABORT('KINSLB: INVALID NUMBER OF UNKNOWNS.')
*----
*  INVERSE POWER METHOD.
*----
      DTF=9999.0D0
      DTP=9999.0D0
      TEST=0.0
      IF(IFL.EQ.1)THEN
        DTF=1.0D0
      ELSEIF(IFL.EQ.2)THEN
        DTF=0.5D0
      ELSEIF(IFL.EQ.3)THEN
        DTF=DBLE(TTF)
      ENDIF
      IF(IPR.EQ.2)THEN
        DTP=0.5D0
      ELSEIF(IPR.EQ.3)THEN
        DTP=DBLE(TTP)
      ENDIF
      DCRIT=MINVAL(DT*PDC(:))
*
      ISTART=1
      IF(IMPX.GE.1) WRITE (6,600)
      IF(IMPX.GE.2) WRITE (6,610)
      M=0
   10 M=M+1
*
      DO 84 IGR=1,NGR
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
      CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,EVECT(1,IGR),WORK1)
      DO 15 IND=1,LL4
      GAR1(IND,IGR)=DTF*WORK1(IND)
   15 CONTINUE
      IF(IEXP.EQ.0) THEN
        DO 16 IBM=1,NBM
        WORK3(IBM)=OVR(IBM,IGR)
   16   CONTINUE
      ELSE
        DO 17 IBM=1,NBM
        WORK3(IBM)=OVR(IBM,IGR)*(1.0+OMEGA(IBM,IGR)*DT)
   17   CONTINUE
      ENDIF
      CALL KINBLM(IPTRK,NBM,LL4,WORK3,EVECT(1,IGR),WORK1)
      DO 20 IND=1,LL4
      GAR1(IND,IGR)=GAR1(IND,IGR)+WORK1(IND)
   20 CONTINUE
      DO 83 JGR=1,NGR
      IF(JGR.EQ.IGR) GO TO 40
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
      CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 40
      CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,EVECT(1,JGR),WORK1)
      DO 30 IND=1,LL4
      GAR1(IND,IGR)=GAR1(IND,IGR)-DTF*WORK1(IND)
   30 CONTINUE
   40 DO 82 IFIS=1,NBFIS
      DO 50 IBM=1,NBM
      WORK3(IBM)=CHI(IBM,IFIS,IGR)*SGF(IBM,IFIS,JGR)
   50 CONTINUE
      CALL KINBLM(IPTRK,NBM,LL4,WORK3,EVECT(1,JGR),WORK1)
      DO 60 IND=1,LL4
      GAR1(IND,IGR)=GAR1(IND,IGR)-DTF*WORK1(IND)
   60 CONTINUE
      DO 81 IDG=1,NDG
      DARG=PDC(IDG)*DT
      IF(IPR.EQ.1)THEN
        DK=1.0D0/(1.0D0+DARG)
      ELSEIF(IPR.EQ.4)THEN
        DK=(1.0D0-DEXP(-DARG))/DARG
      ELSE
        DK=1.0D0/(1.0D0+DTP*DARG)
      ENDIF
      DO 70 IBM=1,NBM
      WORK3(IBM)=CHD(IBM,IFIS,IGR,IDG)*SGD(IBM,IFIS,JGR,IDG)
   70 CONTINUE
      CALL KINBLM(IPTRK,NBM,LL4,WORK3,EVECT(1,JGR),WORK1)
      DO 80 IND=1,LL4
      GAR1(IND,IGR)=GAR1(IND,IGR)+DTF*DK*WORK1(IND)
   80 CONTINUE
   81 CONTINUE
   82 CONTINUE
   83 CONTINUE
   84 CONTINUE
*----
*  DIRECTION EVALUATION.
*----
      DO 120 IGR=1,NGR
      DO 90 IND=1,LL4
      GRAD1(IND,IGR)=REAL(SRC(IND,IGR)-GAR1(IND,IGR))
   90 CONTINUE
      DO 110 JGR=1,IGR-1
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
      CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 110
      CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD1(1,JGR),WORK1)
      DO 100 IND=1,LL4
      GRAD1(IND,IGR)=GRAD1(IND,IGR)+REAL(DTF)*WORK1(IND)
  100 CONTINUE
  110 CONTINUE
      CALL KDRCPU(TK2)
      TKB=TKB+(TK2-TK1)
*
      CALL KDRCPU(TK1)
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
      CALL MTLDLS(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD1(1,IGR))
      CALL KDRCPU(TK2)
      DO 115 IND=1,LL4
      GRAD1(IND,IGR)=GRAD1(IND,IGR)/REAL(DTF)
  115 CONTINUE
      TKT=TKT+(TK2-TK1)
  120 CONTINUE
*----
*  PERFORM THERMAL (UP-SCATTERING) ITERATIONS
*----
      KTER=0
      NADI=5 ! used with SPN approximations
      IF(MAXINR.GT.1) THEN
         CALL FLDBHR(IPTRK,IPSYS,.FALSE.,LL4,ITY,NUN,NGR,ICL1,ICL2,IMPX,
     1   NADI,MAXINR,EPSINR,KTER,TKT,TKB,GRAD1)
      ENDIF
*----
*  EVALUATION OF THE DISPLACEMENT AND OF THE TWO ACCELERATION PARAMETERS
*  ALP AND BET.
*----
      DO 204 IGR=1,NGR
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
      CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD1(1,IGR),WORK1)
      DO 130 IND=1,LL4
      GAR2(IND,IGR)=DTF*WORK1(IND)
  130 CONTINUE
      IF(IEXP.EQ.0) THEN
        DO 135 IBM=1,NBM
        WORK3(IBM)=OVR(IBM,IGR)
  135   CONTINUE
      ELSE
        DO 136 IBM=1,NBM
        WORK3(IBM)=OVR(IBM,IGR)*(1.0+OMEGA(IBM,IGR)*DT)
  136   CONTINUE
      ENDIF
      CALL KINBLM(IPTRK,NBM,LL4,WORK3,GRAD1(1,IGR),WORK1)
      DO 140 IND=1,LL4
      GAR2(IND,IGR)=GAR2(IND,IGR)+WORK1(IND)
  140 CONTINUE
      DO 203 JGR=1,NGR
      IF(JGR.EQ.IGR) GO TO 160
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
      CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 160
      CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD1(1,JGR),WORK1)
      DO 150 IND=1,LL4
      GAR2(IND,IGR)=GAR2(IND,IGR)-DTF*WORK1(IND)
  150 CONTINUE
  160 DO 202 IFIS=1,NBFIS
      DO 170 IBM=1,NBM
      WORK3(IBM)=CHI(IBM,IFIS,IGR)*SGF(IBM,IFIS,JGR)
  170 CONTINUE
      CALL KINBLM(IPTRK,NBM,LL4,WORK3,GRAD1(1,JGR),WORK1)
      DO 180 IND=1,LL4
      GAR2(IND,IGR)=GAR2(IND,IGR)-DTF*WORK1(IND)
  180 CONTINUE
      DO 201 IDG=1,NDG
      DARG=PDC(IDG)*DT
      IF(IPR.EQ.1)THEN
        DK=1.0D0/(1.0D0+DARG)
      ELSEIF(IPR.EQ.4)THEN
        DK=(1.0D0-DEXP(-DARG))/DARG
      ELSE
        DK=1.0D0/(1.0D0+DTP*DARG)
      ENDIF
      DO 190 IBM=1,NBM
      WORK3(IBM)=CHD(IBM,IFIS,IGR,IDG)*SGD(IBM,IFIS,JGR,IDG)
  190 CONTINUE
      CALL KINBLM(IPTRK,NBM,LL4,WORK3,GRAD1(1,JGR),WORK1)
      DO 200 IND=1,LL4
      GAR2(IND,IGR)=GAR2(IND,IGR)+DTF*DK*WORK1(IND)
  200 CONTINUE
  201 CONTINUE
  202 CONTINUE
  203 CONTINUE
  204 CONTINUE
*
  270 ALP=1.0D0
      BET=0.0D0
      CALL XDDSET(D2F,6,0.0D0)
      IF(1+MOD(M-ISTART,ICL1+ICL2).GT.ICL1) THEN
         IF(DCRIT.GT.1.0E-6) THEN
*           TWO-PARAMETER ACCELERATION. SOLUTION OF A LINEAR SYSTEM.
            DO 285 IGR=1,NGR
            DO 280 I=1,LL4
            D2F(1,1)=D2F(1,1)+GAR2(I,IGR)**2
            D2F(1,2)=D2F(1,2)+GAR2(I,IGR)*GAR3(I,IGR)
            D2F(2,2)=D2F(2,2)+GAR3(I,IGR)**2
            D2F(1,3)=D2F(1,3)-(GAR1(I,IGR)-SRC(I,IGR))*GAR2(I,IGR)
            D2F(2,3)=D2F(2,3)-(GAR1(I,IGR)-SRC(I,IGR))*GAR3(I,IGR)
  280       CONTINUE
  285       CONTINUE
            D2F(2,1)=D2F(1,2)
            CALL ALSBD(2,1,D2F,IER,2)
            IF(IER.NE.0) THEN
               DCRIT=1.0E-6
               GO TO 270
            ENDIF
            ALP=D2F(1,3)
            BET=D2F(2,3)/ALP
            IF((ALP.LT.1.0D0).AND.(ALP.GT.0.0D0)) THEN
               ALP=1.0D0
               BET=0.0D0
            ELSE IF(ALP.LE.0.0D0) THEN
               ISTART=M+1
               ALP=1.0D0
               BET=0.0D0
            ENDIF
         ELSE
*           ONE-PARAMETER ACCELERATION.
            DO 295 IGR=1,NGR
            DO 290 I=1,LL4
            D2F(1,1)=D2F(1,1)+GAR2(I,IGR)**2
            D2F(1,3)=D2F(1,3)-(GAR1(I,IGR)-SRC(I,IGR))*GAR2(I,IGR)
  290       CONTINUE
  295       CONTINUE
            IF(D2F(1,1).NE.0.0D0) THEN
               ALP=D2F(1,3)/D2F(1,1)
            ELSE
               ISTART=M+1
            ENDIF
         ENDIF
         DO 305 IGR=1,NGR
         DO 300 I=1,LL4
         GRAD1(I,IGR)=REAL(ALP)*(GRAD1(I,IGR)+REAL(BET)*GRAD2(I,IGR))
         GAR2(I,IGR)=ALP*(GAR2(I,IGR)+BET*GAR3(I,IGR))
  300    CONTINUE
  305    CONTINUE
      ENDIF
*
      LOGTES=(M.LT.ICL1).OR.(MOD(M-ISTART,ICL1+ICL2).EQ.ICL1-1)
      IF(LOGTES) THEN
        ALLOCATE(WORK4(LL4))
        DELT=0.0
        DO 350 IGR=1,NGR
        CALL XDRSET(WORK1,LL4,0.0)
        CALL XDRSET(WORK2,LL4,0.0)
        DO 320 JGR=1,NGR
        WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
        CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
        IF(ILONG.EQ.0) GO TO 320
        CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,EVECT(1,JGR),WORK4)
        DO 310 I=1,LL4
        WORK1(I)=WORK1(I)+WORK4(I)
  310   CONTINUE
        CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD1(1,JGR),WORK4)
        DO 315 I=1,LL4
        WORK2(I)=WORK2(I)+WORK4(I)
  315   CONTINUE
  320   CONTINUE
        DELN=0.0
        DELD=0.0
        DO 340 I=1,LL4
        EVECT(I,IGR)=EVECT(I,IGR)+GRAD1(I,IGR)
        GAR1(I,IGR)=GAR1(I,IGR)+GAR2(I,IGR)
        GRAD2(I,IGR)=GRAD1(I,IGR)
        GAR3(I,IGR)=GAR2(I,IGR)
        DELN=MAX(DELN,ABS(WORK2(I)))
        DELD=MAX(DELD,ABS(WORK1(I)))
  340   CONTINUE
        IF(DELD.NE.0.0) DELT=MAX(DELT,DELN/DELD)
  350   CONTINUE
        DEALLOCATE(WORK4)
        IF(IMPX.GE.2) WRITE (6,620) M,ALP,BET,DELT
*       COMPUTE THE CONVERGENCE HISTOGRAM.
        IF((IMPH.GE.1).AND.(M.LE.250)) THEN
          LMPH=IMPH.GE.1
          CALL FLDXCO(IPKIN,LL4,NUN,EVECT(1,NGR),LMPH,ERR(M))
          ALPH(M)=REAL(ALP)
          BETA(M)=REAL(BET)
        ENDIF
        IF(DELT.LT.EPS2) GO TO 370
      ELSE
        DO 365 IGR=1,NGR
        DO 360 I=1,LL4
        EVECT(I,IGR)=EVECT(I,IGR)+GRAD1(I,IGR)
        GAR1(I,IGR)=GAR1(I,IGR)+GAR2(I,IGR)
        GRAD2(I,IGR)=GRAD1(I,IGR)
        GAR3(I,IGR)=GAR2(I,IGR)
  360   CONTINUE
  365   CONTINUE
        IF(IMPX.GE.2) WRITE (6,620) M,ALP,BET
*       COMPUTE THE CONVERGENCE HISTOGRAM.
        IF((IMPH.GE.1).AND.(M.LE.250)) THEN
          LMPH=IMPH.GE.1
          CALL FLDXCO(IPKIN,LL4,NUN,EVECT(1,NGR),LMPH,ERR(M))
          ALPH(M)=REAL(ALP)
          BETA(M)=REAL(BET)
        ENDIF
      ENDIF
      IF(M.EQ.1) TEST=DELT
      IF((M.GT.30).AND.(DELT.GT.TEST)) CALL XABORT('KINSLB: CONVERGENC'
     1 //'E FAILURE.')
      IF(M.GE.MIN(MAXX0,MMAXX)) THEN
         WRITE (6,710)
         GO TO 370
      ENDIF
      GO TO 10
*----
*  SOLUTION EDITION.
*----
  370 IF(IMPX.EQ.1) WRITE (6,640) M
      IF(IMPX.GE.3) THEN
         DO 380 IGR=1,NGR
         WRITE (6,690) IGR,(EVECT(I,IGR),I=1,LL4)
  380    CONTINUE
      ENDIF
      IF(IMPH.GE.2) THEN
         IGRAPH=0
  390    IGRAPH=IGRAPH+1
         WRITE (TEXT12,'(5HHISTO,I3)') IGRAPH
         CALL LCMLEN (IPKIN,TEXT12,ILENG,ITYLCM)
         IF(ILENG.EQ.0) THEN
            MDIM=MIN(250,M)
            READ (TITR,'(18A4)') ITITR
            CALL LCMSIX (IPKIN,TEXT12,1)
            CALL LCMPUT (IPKIN,'HTITLE',18,3,ITITR)
            CALL LCMPUT (IPKIN,'ALPHA',MDIM,2,ALPH)
            CALL LCMPUT (IPKIN,'BETA',MDIM,2,BETA)
            CALL LCMPUT (IPKIN,'ERROR',MDIM,2,ERR)
            CALL LCMPUT (IPKIN,'IMPH',1,1,IMPH)
            CALL LCMSIX (IPKIN,' ',2)
         ELSE
            GO TO 390
         ENDIF
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GRAD1,GRAD2,GAR1,GAR2,GAR3,WORK1,WORK2,WORK3)
      RETURN
*
  600 FORMAT(1H1/50H KINSLB: ITERATIVE PROCEDURE BASED ON INVERSE POWE,
     1 8HR METHOD/9X,30HSPACE-TIME KINETICS EQUATIONS.)
  610 FORMAT(/11X,5HALPHA,3X,4HBETA,6X,8HACCURACY,12(1H.))
  620 FORMAT(1X,I3,4X,2F8.3,1PE13.2)
  640 FORMAT(/23H KINSLB: CONVERGENCE IN,I4,12H ITERATIONS.)
  690 FORMAT(//52H KINSLB: SPACE-TIME KINETICS SOLUTION CORRESPONDING ,
     1 12HTO THE GROUP,I4//(5X,1P,8E14.5))
  710 FORMAT(/53H KINSLB: ***WARNING*** THE MAXIMUM NUMBER OF OUTER IT,
     1 20HERATIONS IS REACHED.)
      END
