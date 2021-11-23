*DECK FLDDIR
      SUBROUTINE FLDDIR(IPTRK,IPSYS,IPFLUX,LL4,ITY,NUN,NGRP,ICL1,ICL2,
     1 IMPX,IMPH,TITR,EPS2,NADI,MAXOUT,MAXINR,EPSINR,EVECT,FKEFF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of a multigroup eigenvalue system for the calculation of the
* direct neutron flux in Trivac. Use the preconditioned power method
* with a two-parameter SVAT acceleration technique.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
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
* IPFLUX  L_FLUX pointer to the solution.
* LL4     order of the system matrices.
* ITY     type of solution (2: classical Trivac; 3: Thomas-Raviart).
* NUN     number of unknowns in each energy group.
* NGRP    number of energy groups.
* ICL1    number of free iterations in one cycle of the inverse power
*         method.
* ICL2    number of accelerated iterations in one cycle.
* IMPX    print parameter: =0: no print; =1: minimum printing;
*         =2: iteration history is printed; =3: solution is printed.
* IMPH    type of histogram processing:
*         =0: no action is taken;
*         =1: the flux is compared to a reference flux stored on LCM;
*         =2: the convergence histogram is printed;
*         =3: the convergence histogram is printed with axis and
*            titles. The plotting file is completed;
*         =4: the convergence histogram is printed with axis, acce-
*            leration factors and titles. The plotting file is
*            completed.
* TITR    title.
* EPS2    convergence criteria for the flux.
* NADI    number of inner ADI iterations per outer iteration.
* MAXOUT  maximum number of outer iterations.
* MAXINR  maximum number of thermal iterations.
* EPSINR  thermal iteration epsilon.
* EVECT   initial estimate of the unknown vector.
*
*Parameters: output
* FKEFF   effective multiplication factor.
* EVECT   converged unknown vector.
*
*Reference:
* A. H\'ebert, 'Preconditioning the power method for reactor
* calculations', Nucl. Sci. Eng., 94, 1 (1986).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS,IPFLUX
      CHARACTER TITR*72
      INTEGER LL4,ITY,NUN,NGRP,ICL1,ICL2,IMPX,IMPH,NADI,MAXOUT,MAXINR
      REAL FKEFF,EPS2,EPSINR,EVECT(NUN,NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MMAXX=250,EPS1=1.0E-5)
      CHARACTER TEXT12*12
      LOGICAL LOGTES
      DOUBLE PRECISION AEAE,AEAG,AEAH,AGAG,AGAH,AHAH,BEBE,BEBG,BEBH,
     1 BGBG,BGBH,BHBH,AEBE,AEBG,AEBH,AGBE,AGBG,AGBH,AHBE,AHBG,AHBH,
     2 X,DXDA,DXDB,Y,DYDA,DYDB,Z,DZDA,DZDB,F,D2F(2,3),EVAL,ALP,BET,
     3 FMIN
      REAL ERR(MMAXX),ALPH(MMAXX),BETA(MMAXX)
      DOUBLE PRECISION, PARAMETER :: ALP_TAB(24) = (/ 0.2, 0.4, 0.6,
     1  0.8, 1.0, 1.2, 1.5, 2.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0,
     2  40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0 /)
      DOUBLE PRECISION, PARAMETER :: BET_TAB(11) = (/ -1.0, -0.8, -0.6,
     1 -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 /)
      REAL, DIMENSION(:,:), ALLOCATABLE :: GRAD1,GRAD2,GAR1,GAR2,GAR3
      REAL, DIMENSION(:), ALLOCATABLE :: GAF1,GAF2,GAF3
      REAL, DIMENSION(:), POINTER :: AGAR
      TYPE(C_PTR) AGAR_PTR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(GRAD1(NUN,NGRP),GRAD2(NUN,NGRP),GAR1(NUN,NGRP),
     1 GAR2(NUN,NGRP),GAR3(NUN,NGRP),GAF1(NUN),GAF2(NUN),GAF3(NUN))
*
*     TKT : CPU TIME FOR THE SOLUTION OF LINEAR SYSTEMS.
*     TKB : CPU TIME FOR BILINEAR PRODUCT EVALUATIONS.
      TKT=0.0
      TKB=0.0
      CALL KDRCPU(TK1)
      CALL MTOPEN(IMPX,IPTRK,LL4)
      IF(LL4.GT.NUN) CALL XABORT('FLDDIR: INVALID NUMBER OF UNKNOWNS.')
*----
*  PRECONDITIONED POWER METHOD
*----
      EVAL=1.0D0
      VVV=0.0
      ISTART=1
      NNADI=NADI
      TEST=0.0
      IF(IMPX.GE.1) WRITE (6,600) NADI
      IF(IMPX.GE.2) WRITE (6,610)
      DO 35 IGR=1,NGRP
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
      CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,EVECT(1,IGR),GAR1(1,IGR))
      DO 30 JGR=1,NGRP
      IF(JGR.EQ.IGR) GO TO 30
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
      CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 30
      IF(ITY.EQ.13) THEN
         CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,EVECT(1,JGR),GAF1)
         DO 10 I=1,LL4
         GAR1(I,IGR)=GAR1(I,IGR)-GAF1(I)
   10    CONTINUE
      ELSE
         CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ ILONG /))
         DO 20 I=1,ILONG
         GAR1(I,IGR)=GAR1(I,IGR)-AGAR(I)*EVECT(I,JGR)
   20    CONTINUE
      ENDIF
   30 CONTINUE
   35 CONTINUE
      CALL KDRCPU(TK2)
      TKB=TKB+(TK2-TK1)
*
      M=0
   40 M=M+1
*----
*  EIGENVALUE EVALUATION
*----
      CALL KDRCPU(TK1)
      AEBE=0.0D0
      BEBE=0.0D0
      DO 95 IGR=1,NGRP
      DO 50 I=1,LL4
      GAF1(I)=0.0
   50 CONTINUE
      DO 80 JGR=1,NGRP
      WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
      CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 80
      CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
      CALL C_F_POINTER(AGAR_PTR,AGAR,(/ ILONG /))
      DO 60 I=1,ILONG
      GAF1(I)=GAF1(I)+AGAR(I)*EVECT(I,JGR)
   60 CONTINUE
   80 CONTINUE
      DO 90 I=1,LL4
      AEBE=AEBE+GAR1(I,IGR)*GAF1(I)
      BEBE=BEBE+GAF1(I)**2
      GRAD1(I,IGR)=GAF1(I)
   90 CONTINUE
   95 CONTINUE
      EVAL=AEBE/BEBE
      CALL KDRCPU(TK2)
      TKB=TKB+(TK2-TK1)
*----
*  DIRECTION EVALUATION
*----
      DO 140 IGR=1,NGRP
      CALL KDRCPU(TK1)
      DO 100 I=1,LL4
      GRAD1(I,IGR)=REAL(EVAL)*GRAD1(I,IGR)-GAR1(I,IGR)
  100 CONTINUE
      DO 130 JGR=1,IGR-1
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
      CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 130
      IF(ITY.EQ.13) THEN
         CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD1(1,JGR),GAF1)
         DO 110 I=1,LL4
         GRAD1(I,IGR)=GRAD1(I,IGR)+GAF1(I)
  110    CONTINUE
      ELSE
         CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ ILONG /))
         DO 120 I=1,ILONG
         GRAD1(I,IGR)=GRAD1(I,IGR)+AGAR(I)*GRAD1(I,JGR)
  120    CONTINUE
      ENDIF
  130 CONTINUE
      CALL KDRCPU(TK2)
      TKB=TKB+(TK2-TK1)
*
      CALL KDRCPU(TK1)
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
      CALL FLDADI(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD1(1,IGR),NNADI)
      CALL KDRCPU(TK2)
      TKT=TKT+(TK2-TK1)
  140 CONTINUE
*----
*  PERFORM THERMAL (UP-SCATTERING) ITERATIONS
*----
      IF(MAXINR.GT.1) THEN
         CALL FLDTHR(IPTRK,IPSYS,IPFLUX,.FALSE.,LL4,ITY,NUN,NGRP,ICL1,
     1   ICL2,IMPX,NNADI,0,MAXINR,EPSINR,ITER,TKT,TKB,GRAD1)
      ENDIF
*----
*  DISPLACEMENT EVALUATION
*----
      F=0.0D0
      DELS=ABS(REAL((EVAL-VVV)/EVAL))
      VVV=REAL(EVAL)
      CALL KDRCPU(TK1)
*----
*  EVALUATION OF THE TWO ACCELERATION PARAMETERS ALP AND BET
*----
      ALP=1.0D0
      BET=0.0D0
      N=0
      AEAE=0.0D0
      AEAG=0.0D0
      AEAH=0.0D0
      AGAG=0.0D0
      AGAH=0.0D0
      AHAH=0.0D0
      BEBG=0.0D0
      BEBH=0.0D0
      BGBG=0.0D0
      BGBH=0.0D0
      BHBH=0.0D0
      AEBG=0.0D0
      AEBH=0.0D0
      AGBE=0.0D0
      AGBG=0.0D0
      AGBH=0.0D0
      AHBE=0.0D0
      AHBG=0.0D0
      AHBH=0.0D0
      DO 175 IGR=1,NGRP
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
      CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD1(1,IGR),GAR2(1,IGR))
      DO 170 JGR=1,NGRP
      IF(JGR.EQ.IGR) GO TO 170
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
      CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 170
      IF(ITY.EQ.13) THEN
         CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD1(1,JGR),GAF1)
         DO 150 I=1,LL4
         GAR2(I,IGR)=GAR2(I,IGR)-GAF1(I)
  150    CONTINUE
      ELSE
         CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ ILONG /))
         DO 160 I=1,ILONG
         GAR2(I,IGR)=GAR2(I,IGR)-AGAR(I)*GRAD1(I,JGR)
  160    CONTINUE
      ENDIF
  170 CONTINUE
  175 CONTINUE
      IF(1+MOD(M-ISTART,ICL1+ICL2).GT.ICL1) THEN
         DO 205 IGR=1,NGRP
         CALL XDRSET(GAF1,LL4,0.0)
         CALL XDRSET(GAF2,LL4,0.0)
         CALL XDRSET(GAF3,LL4,0.0)
         DO 190 JGR=1,NGRP
         WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 190
         CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ ILONG /))
         DO 180 I=1,ILONG
         GAF1(I)=GAF1(I)+AGAR(I)*EVECT(I,JGR)
         GAF2(I)=GAF2(I)+AGAR(I)*GRAD1(I,JGR)
         GAF3(I)=GAF3(I)+AGAR(I)*GRAD2(I,JGR)
  180    CONTINUE
  190    CONTINUE
         DO 200 I=1,LL4
*        COMPUTE (A ,A )
         AEAE=AEAE+GAR1(I,IGR)**2
         AEAG=AEAG+GAR1(I,IGR)*GAR2(I,IGR)
         AEAH=AEAH+GAR1(I,IGR)*GAR3(I,IGR)
         AGAG=AGAG+GAR2(I,IGR)**2
         AGAH=AGAH+GAR2(I,IGR)*GAR3(I,IGR)
         AHAH=AHAH+GAR3(I,IGR)**2
*        COMPUTE (B ,B )
         BEBG=BEBG+GAF1(I)*GAF2(I)
         BEBH=BEBH+GAF1(I)*GAF3(I)
         BGBG=BGBG+GAF2(I)**2
         BGBH=BGBH+GAF2(I)*GAF3(I)
         BHBH=BHBH+GAF3(I)**2
*        COMPUTE (A ,B )
         AEBG=AEBG+GAR1(I,IGR)*GAF2(I)
         AEBH=AEBH+GAR1(I,IGR)*GAF3(I)
         AGBE=AGBE+GAR2(I,IGR)*GAF1(I)
         AGBG=AGBG+GAR2(I,IGR)*GAF2(I)
         AGBH=AGBH+GAR2(I,IGR)*GAF3(I)
         AHBE=AHBE+GAR3(I,IGR)*GAF1(I)
         AHBG=AHBG+GAR3(I,IGR)*GAF2(I)
         AHBH=AHBH+GAR3(I,IGR)*GAF3(I)
  200    CONTINUE
  205    CONTINUE
*
  210    N=N+1
         IF(N.GT.10) THEN
            IF(IMPX.GT.0) WRITE(6,'(/30H FLDDIR: FAILURE OF THE NEWTON,
     1      55H-RAPHTON ALGORIHTHM FOR COMPUTING THE OVERRELAXATION PA,
     2      12HRAMETERS(1).)')
            GO TO 215
         ENDIF
*        COMPUTE X(M+1)
         X=BEBE+ALP*ALP*BGBG+BET*BET*BHBH+2.0D0*(ALP*BEBG+BET*BEBH
     1   +ALP*BET*BGBH)
         DXDA=2.0D0*(BEBG+ALP*BGBG+BET*BGBH)
         DXDB=2.0D0*(BEBH+ALP*BGBH+BET*BHBH)
*        COMPUTE Y(M+1)
         Y=AEAE+ALP*ALP*AGAG+BET*BET*AHAH+2.0D0*(ALP*AEAG+BET*AEAH
     1   +ALP*BET*AGAH)
         DYDA=2.0D0*(AEAG+ALP*AGAG+BET*AGAH)
         DYDB=2.0D0*(AEAH+ALP*AGAH+BET*AHAH)
*        COMPUTE Z(M+1)
         Z=AEBE+ALP*ALP*AGBG+BET*BET*AHBH+ALP*(AEBG+AGBE)
     1   +BET*(AEBH+AHBE)+ALP*BET*(AGBH+AHBG)
         DZDA=AEBG+AGBE+2.0D0*ALP*AGBG+BET*(AGBH+AHBG)
         DZDB=AEBH+AHBE+ALP*(AGBH+AHBG)+2.0D0*BET*AHBH
*        COMPUTE F(M+1)
         F=X*Y-Z*Z
         D2F(1,1)=2.0D0*(BGBG*Y+DXDA*DYDA+X*AGAG-DZDA**2-2.0D0*Z*AGBG)
         D2F(1,2)=2.0D0*BGBH*Y+DXDA*DYDB+DXDB*DYDA+2.0D0*X*AGAH
     1   -2.0D0*DZDA*DZDB-2.0D0*Z*(AGBH+AHBG)
         D2F(2,2)=2.0D0*(BHBH*Y+DXDB*DYDB+X*AHAH-DZDB**2-2.0D0*Z*AHBH)
         D2F(2,1)=D2F(1,2)
         D2F(1,3)=DXDA*Y+X*DYDA-2.0D0*Z*DZDA
         D2F(2,3)=DXDB*Y+X*DYDB-2.0D0*Z*DZDB
*        SOLUTION OF A LINEAR SYSTEM.
         CALL ALSBD(2,1,D2F,IER,2)
         IF(IER.NE.0) CALL XABORT('FLDDIR: SINGULAR MATRIX.')
         ALP=ALP-D2F(1,3)
         BET=BET-D2F(2,3)
         IF(ALP.GT.100.0) THEN
            IF(IMPX.GT.0) WRITE(6,'(/30H FLDDIR: FAILURE OF THE NEWTON,
     1      55H-RAPHTON ALGORIHTHM FOR COMPUTING THE OVERRELAXATION PA,
     2      12HRAMETERS(2).)')
            GO TO 215
         ENDIF
         IF((ABS(D2F(1,3)).LE.1.0D-4).AND.(ABS(D2F(2,3)).LE.1.0D-4))
     1   GO TO 220
         GO TO 210
*
*        alternative algorithm in case of Newton-Raphton failure
  215    IAMIN=999
         IBMIN=999
         FMIN=HUGE(FMIN)
         DO IA=1,SIZE(ALP_TAB)
           ALP=ALP_TAB(IA)
           DO IB=1,SIZE(BET_TAB)
             BET=BET_TAB(IB)
*            COMPUTE X
             X=BEBE+ALP*ALP*BGBG+BET*BET*BHBH+2.0D0*(ALP*BEBG+BET*BEBH
     1       +ALP*BET*BGBH)
*            COMPUTE Y
             Y=AEAE+ALP*ALP*AGAG+BET*BET*AHAH+2.0D0*(ALP*AEAG+BET*AEAH
     1       +ALP*BET*AGAH)
*            COMPUTE Z
             Z=AEBE+ALP*ALP*AGBG+BET*BET*AHBH+ALP*(AEBG+AGBE)
     1       +BET*(AEBH+AHBE)+ALP*BET*(AGBH+AHBG)
*            COMPUTE F
             F=X*Y-Z*Z
             IF(F.LT.FMIN) THEN
               IAMIN=IA
               IBMIN=IB
               FMIN=F
             ENDIF
           ENDDO
         ENDDO
         ALP=ALP_TAB(IAMIN)
         BET=BET_TAB(IBMIN)
  220    BET=BET/ALP
         IF((ALP.LT.1.0D0).AND.(ALP.GT.0.0D0)) THEN
            ALP=1.0D0
            BET=0.0D0
         ELSE IF(ALP.LE.0.0D0) THEN
            ISTART=M+1
            ALP=1.0D0
            BET=0.0D0
         ENDIF
         DO 235 IGR=1,NGRP
         DO 230 I=1,LL4
         GRAD1(I,IGR)=REAL(ALP)*(GRAD1(I,IGR)+REAL(BET)*GRAD2(I,IGR))
         GAR2(I,IGR)=REAL(ALP)*(GAR2(I,IGR)+REAL(BET)*GAR3(I,IGR))
  230    CONTINUE
  235    CONTINUE
      ENDIF
      CALL KDRCPU(TK2)
      TKB=TKB+(TK2-TK1)
*
      LOGTES=(M.LT.ICL1).OR.(MOD(M-ISTART,ICL1+ICL2).EQ.ICL1-1)
      IF(LOGTES.AND.(DELS.LE.EPS1)) THEN
         DELT=0.0
         DO 290 IGR=1,NGRP
         CALL XDRSET(GAF1,LL4,0.0)
         CALL XDRSET(GAF2,LL4,0.0)
         DO 250 JGR=1,NGRP
         WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 250
         CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ ILONG /))
         DO 240 I=1,ILONG
         GAF1(I)=GAF1(I)+AGAR(I)*EVECT(I,JGR)
         GAF2(I)=GAF2(I)+AGAR(I)*GRAD1(I,JGR)
  240    CONTINUE
  250    CONTINUE
         DELN=0.0
         DELD=0.0
         DO 280 I=1,LL4
         EVECT(I,IGR)=EVECT(I,IGR)+GRAD1(I,IGR)
         GAR1(I,IGR)=GAR1(I,IGR)+GAR2(I,IGR)
         GRAD2(I,IGR)=GRAD1(I,IGR)
         GAR3(I,IGR)=GAR2(I,IGR)
         DELN=MAX(DELN,ABS(GAF2(I)))
         DELD=MAX(DELD,ABS(GAF1(I)))
  280    CONTINUE
         IF(DELD.NE.0.0) DELT=MAX(DELT,DELN/DELD)
  290    CONTINUE
         IF(IMPX.GE.2) WRITE (6,615) M,AEAE,AEAG,AEAH,AGAG,AGAH,AHAH,
     1   BEBE,ALP,REAL(BET),EVAL,F,DELS,DELT,N,BEBG,BEBH,BGBG,BGBH,BHBH,
     2   AEBE,AEBG,AEBH,AGBE,AGBG,AGBH,AHBE,AHBG,AHBH
*        COMPUTE THE CONVERGENCE HISTOGRAM.
         IF((IMPH.GE.1).AND.(M.LE.MMAXX)) THEN
            CALL FLDXCO(IPFLUX,LL4,NUN,EVECT(1,NGRP),.TRUE.,ERR(M))
            ALPH(M)=REAL(ALP)
            BETA(M)=REAL(BET)
         ENDIF
         IF(DELT.LE.EPS2) GO TO 310
      ELSE
         DO 305 IGR=1,NGRP
         DO 300 I=1,LL4
         EVECT(I,IGR)=EVECT(I,IGR)+GRAD1(I,IGR)
         GAR1(I,IGR)=GAR1(I,IGR)+GAR2(I,IGR)
         GRAD2(I,IGR)=GRAD1(I,IGR)
         GAR3(I,IGR)=GAR2(I,IGR)
  300    CONTINUE
  305    CONTINUE
         IF(IMPX.GE.2) WRITE (6,620) M,AEAE,AEAG,AEAH,AGAG,AGAH,AHAH,
     1   BEBE,ALP,BET,EVAL,F,DELS,N,BEBG,BEBH,BGBG,BGBH,BHBH,AEBE,
     2   AEBG,AEBH,AGBE,AGBG,AGBH,AHBE,AHBG,AHBH
*        COMPUTE THE CONVERGENCE HISTOGRAM.
         IF((IMPH.GE.1).AND.(M.LE.MMAXX)) THEN
            CALL FLDXCO(IPFLUX,LL4,NUN,EVECT(1,NGRP),.TRUE.,ERR(M))
            ALPH(M)=REAL(ALP)
            BETA(M)=REAL(BET)
         ENDIF
      ENDIF
*
      IF(M.EQ.1) TEST=DELS
      IF((M.GT.5).AND.(DELS.GT.TEST)) CALL XABORT('FLDDIR: CONVERGENCE'
     1 //' FAILURE.')
      IF(M.GE.MAXOUT) THEN
         WRITE (6,690)
         GO TO 310
      ENDIF
      IF(MOD(M,36).EQ.0) THEN
         ISTART=M+1
         NNADI=NNADI+1
         IF(IMPX.GE.1) WRITE (6,700) NNADI
      ENDIF
      GO TO 40
*----
*  SOLUTION EDITION
*----
  310 FKEFF=REAL(1.0D0/EVAL)
      IF(IMPX.EQ.1) WRITE (6,640) M
      IF(IMPX.GE.1) THEN
         WRITE (6,650) TKT,TKB,TKT+TKB
         WRITE (6,670) FKEFF
      ENDIF
      IF(IMPX.EQ.3) THEN
         DO 320 IGR=1,NGRP
         WRITE (6,680) IGR,(EVECT(I,IGR),I=1,LL4)
  320    CONTINUE
      ENDIF
      IF(IMPH.EQ.1) THEN
         CALL LCMLEN(IPFLUX,'REF',ILONG,ITYLCM)
         IF(ILONG.EQ.0) THEN
            WRITE(6,'(40H FLDDIR: STORE A REFERENCE THERMAL FLUX.)')
            CALL LCMPUT(IPFLUX,'REF',NUN,2,EVECT(1,NGRP))
         ENDIF
      ELSE IF(IMPH.GE.2) THEN
         IGRAPH=0
  330    IGRAPH=IGRAPH+1
         WRITE (TEXT12,'(5HHISTO,I3)') IGRAPH
         CALL LCMLEN (IPFLUX,TEXT12,ILENG,ITYLCM)
         IF(ILENG.EQ.0) THEN
            MM=MIN(M,MMAXX)
            CALL LCMSIX (IPFLUX,TEXT12,1)
            CALL LCMPTC (IPFLUX,'HTITLE',72,1,TITR)
            CALL LCMPUT (IPFLUX,'ALPHA',MM,2,ALPH)
            CALL LCMPUT (IPFLUX,'BETA',MM,2,BETA)
            CALL LCMPUT (IPFLUX,'ERROR',MM,2,ERR)
            CALL LCMPUT (IPFLUX,'IMPH',1,1,IMPH)
            CALL LCMSIX (IPFLUX,' ',2)
         ELSE
            GO TO 330
         ENDIF
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GRAD1,GRAD2,GAR1,GAR2,GAR3,GAF1,GAF2,GAF3)
      RETURN
*
  600 FORMAT(1H1/50H FLDDIR: ITERATIVE PROCEDURE BASED ON PRECONDITION,
     1 17HED POWER METHOD (,I2,37H ADI ITERATIONS PER OUTER ITERATION)./
     2 9X,16HDIRECT EQUATION.)
  610 FORMAT(//5X,17HBILINEAR PRODUCTS,48X,5HALPHA,3X,4HBETA,3X,
     1 12HEIGENVALUE..,12X,8HACCURACY,11(1H.),2X,1HN)
  615 FORMAT(1X,I3,1P,7E9.1,0P,2F8.3,E14.6,3E10.2,I4/
     1 (4X,1P,7E9.1))
  620 FORMAT(1X,I3,1P,7E9.1,0P,2F8.3,E14.6,2E10.2,10X,I4/
     1 (4X,1P,7E9.1))
  640 FORMAT(/23H FLDDIR: CONVERGENCE IN,I4,12H ITERATIONS.)
  650 FORMAT(/53H FLDDIR: CPU TIME USED TO SOLVE THE TRIANGULAR LINEAR,
     1 10H SYSTEMS =,F10.3/23X,34HTO COMPUTE THE BILINEAR PRODUCTS =,
     2 F10.3,20X,16HTOTAL CPU TIME =,F10.3)
  670 FORMAT(//42H FLDDIR: EFFECTIVE MULTIPLICATION FACTOR =,1P,E17.10/)
  680 FORMAT(//47H FLDDIR: EIGENVECTOR CORRESPONDING TO THE GROUP,I4
     1 //(5X,1P,8E14.5))
  690 FORMAT(/53H FLDDIR: ***WARNING*** THE MAXIMUM NUMBER OF OUTER IT,
     1 20HERATIONS IS REACHED.)
  700 FORMAT(/53H FLDDIR: INCREASING THE NUMBER OF INNER ITERATIONS TO,
     1 I3,36H ADI ITERATIONS PER OUTER ITERATION./)
      END
