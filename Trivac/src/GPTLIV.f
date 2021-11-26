*DECK GPTLIV
      SUBROUTINE GPTLIV(IPTRK,IPSYS,IPFLUP,LADJ,LL4,ITY,NUN,NGRP,ICL1,
     1 ICL2,IMPX,IMPH,TITR,NADI,MAXINR,MAXX0,EPS2,EPSINR,EVAL,EVECT,
     2 ADECT,EASS,SOUR,TKT,TKB,ZNORM,M)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of a multigroup fixed source eigenvalue problem for the
* calculation of a gpt solution in Trivac. Use the preconditioned power
* method with two parameter SVAT  acceleration.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
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
* IPFLUP  L_FLUX pointer to the gpt solution
* LADJ    flag set to .TRUE. for adjoint solution acceleration.
* LL4     order of the system matrices.
* ITY     type of solution (2: classical Trivac; 3: Thomas-Raviart).
* NUN     number of unknowns in each energy group.
* NGRP    number of energy groups.
* ICL1    number of free up-scattering iterations in one cycle of the
*         inverse power method.
* ICL2    number of accelerated up-scattering iterations in one cycle.
* IMPX    print parameter. =0: no print; =1: minimum printing;
*         =2: iteration history is printed; =3: solution is printed.
* IMPH    =0: no action is taken
*         =1: the flux is compared to a reference flux stored on lcm
*         =2: the convergence histogram is printed
*         =3: the convergence histogram is printed with axis and
*            titles. the plotting file is completed
*         =4: the convergence histogram is printed with axis, acce-
*            leration factors and titles. the plotting file is
*            completed.
* TITR    character*72 title
* NADI    initial number of inner ADI iterations per outer iteration.
* MAXINR  maximum number of thermal iterations.
* EPSINR  thermal iteration epsilon.
* EVAL    eigenvalue.
* EVECT   unknown vector for the non perturbed direct flux
* ADECT   unknown vector for the non perturbed adjoint flux
* EASS    solution of the fixed source eigenvalue problem
* SOUR    fixed source
*
*Parameters: input/output
* TKT     CPU time spent to compute the solution of linear systems.
* TKB     CPU time spent to compute the bilinear products.
* ZNORM   Hotelling deflation accuracy.
* M       number of iterations.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS,IPFLUP
      CHARACTER TITR*72
      LOGICAL LADJ
      INTEGER LL4,ITY,NUN,NGRP,ICL1,ICL2,IMPX,IMPH,NNADI,MAXINR,MAXX0,M
      REAL EPS2,EPSINR,EVECT(NUN,NGRP),ADECT(NUN,NGRP),EASS(NUN,NGRP),
     1 SOUR(NUN,NGRP),TKT,TKB
      DOUBLE PRECISION EVAL,ZNORM
*----
*  LOCAL VARIABLES
*----
      CHARACTER*12 TEXT12
      LOGICAL LGAR1,LOGTES,LMPH
      DOUBLE PRECISION D2F(2,3),ALP,BET
      REAL ERR(250),ALPH(250),BETA(250)
      REAL, DIMENSION(:,:), ALLOCATABLE :: GRAD1,GRAD2,GAR1,GAR2,GAR3
      REAL, DIMENSION(:), ALLOCATABLE :: WORK1,WORK2
      REAL, DIMENSION(:), POINTER :: AGAR
      TYPE(C_PTR) AGAR_PTR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(GRAD1(NUN,NGRP),GRAD2(NUN,NGRP),GAR1(NUN,NGRP),
     1 GAR2(NUN,NGRP),GAR3(NUN,NGRP),WORK1(NUN),WORK2(NUN))
*
      TEST=0.0
      ISTART=1
      NNADI=NADI
      IF(IMPX.GE.2) WRITE(6,500)
      M=0
  100 M=M+1
*
      LGAR1=(MOD(M-ISTART+1,ICL1+ICL2).EQ.1).OR.(M.EQ.1)
      CALL GPTGRA(IPTRK,IPSYS,IPFLUP,LADJ,LGAR1,LL4,ITY,NUN,NGRP,ICL1,
     1 ICL2,IMPX,NNADI,MAXINR,EPSINR,EVAL,EVECT,ADECT,EASS,SOUR,GAR1,
     2 ITER,TKT,TKB,ZNORM,GRAD1)
*----
*  EVALUATION OF THE DISPLACEMENT AND OF THE TWO ACCELERATION PARAMETERS
*  ALP AND BET.
*----
      ALP=1.0D0
      BET=0.0D0
      DO 240 I=1,2
      DO 230 J=1,3
      D2F(I,J)=0.0D0
  230 CONTINUE
  240 CONTINUE
      DO 285 IGR=1,NGRP
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
      CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD1(1,IGR),GAR2(1,IGR))
      DO 280 JGR=1,NGRP
      IF(JGR.EQ.IGR) GO TO 260
      IF(LADJ) THEN
        WRITE(TEXT12,'(1HA,2I3.3)') JGR,IGR
      ELSE
        WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
      ENDIF
      CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 260
      IF(ITY.EQ.13) THEN
         CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD1(1,JGR),WORK1)
         DO 245 I=1,LL4
         GAR2(I,IGR)=GAR2(I,IGR)-WORK1(I)
  245    CONTINUE
      ELSE
         CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
         DO 250 I=1,ILONG
         GAR2(I,IGR)=GAR2(I,IGR)-AGAR(I)*GRAD1(I,JGR)
  250    CONTINUE
      ENDIF
  260 IF(LADJ) THEN
        WRITE(TEXT12,'(1HB,2I3.3)') JGR,IGR
      ELSE
        WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
      ENDIF
      CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 280
      CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
      CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
      DO 270 I=1,ILONG
      GAR2(I,IGR)=GAR2(I,IGR)-REAL(EVAL)*AGAR(I)*GRAD1(I,JGR)
  270 CONTINUE
  280 CONTINUE
  285 CONTINUE
      IF(1+MOD(M-ISTART,ICL1+ICL2).GT.ICL1) THEN
         DO 295 IGR=1,NGRP
         DO 290 I=1,LL4
         D2F(1,1)=D2F(1,1)+GAR2(I,IGR)**2
         D2F(1,2)=D2F(1,2)+GAR2(I,IGR)*GAR3(I,IGR)
         D2F(2,2)=D2F(2,2)+GAR3(I,IGR)**2
         D2F(1,3)=D2F(1,3)-(GAR1(I,IGR)+SOUR(I,IGR))*GAR2(I,IGR)
         D2F(2,3)=D2F(2,3)-(GAR1(I,IGR)+SOUR(I,IGR))*GAR3(I,IGR)
  290    CONTINUE
  295    CONTINUE
         D2F(2,1)=D2F(1,2)
*        SOLUTION OF A LINEAR SYSTEM.
         CALL ALSBD(2,1,D2F,IER,2)
         IF(IER.NE.0) CALL XABORT('GPTLIV: SINGULAR MATRIX.')
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
         DO 305 IGR=1,NGRP
         DO 300 I=1,LL4
         GRAD1(I,IGR)=REAL(ALP)*(GRAD1(I,IGR)+REAL(BET)*GRAD2(I,IGR))
         GAR2(I,IGR)=REAL(ALP)*(GAR2(I,IGR)+REAL(BET)*GAR3(I,IGR))
  300    CONTINUE
  305    CONTINUE
      ENDIF
*
      LOGTES=(M.LT.ICL1).OR.(MOD(M-ISTART,ICL1+ICL2).EQ.ICL1-1)
      IF(LOGTES) THEN
         DELT=0.0
         DO 350 IGR=1,NGRP
         CALL XDRSET(WORK1,LL4,0.0)
         CALL XDRSET(WORK2,LL4,0.0)
         DO 320 JGR=1,NGRP
         IF(LADJ) THEN
           WRITE(TEXT12,'(1HB,2I3.3)') JGR,IGR
         ELSE
           WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
         ENDIF
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 320
         CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
         DO 310 I=1,ILONG
         WORK1(I)=WORK1(I)+AGAR(I)*EASS(I,JGR)
         WORK2(I)=WORK2(I)+AGAR(I)*GRAD1(I,JGR)
  310    CONTINUE
  320    CONTINUE
         DELN=0.0
         DELD=0.0
         DO 340 I=1,LL4
         EASS(I,IGR)=EASS(I,IGR)+GRAD1(I,IGR)
         GAR1(I,IGR)=GAR1(I,IGR)+GAR2(I,IGR)
         GRAD2(I,IGR)=GRAD1(I,IGR)
         GAR3(I,IGR)=GAR2(I,IGR)
         DELN=MAX(DELN,ABS(WORK2(I)))
         DELD=MAX(DELD,ABS(WORK1(I)))
  340    CONTINUE
         IF(DELD.NE.0.0) DELT=MAX(DELT,DELN/DELD)
  350    CONTINUE
         IF(IMPX.GE.2) WRITE(6,510) M,ALP,BET,ZNORM,DELT,ITER
*        COMPUTE THE CONVERGENCE HISTOGRAM.
         IF(IMPH.GE.1) THEN
            LMPH=IMPH.GE.1
            CALL FLDXCO(IPFLUP,LL4,NUN,EASS(1,NGRP),LMPH,ERR(M))
            ALPH(M)=REAL(ALP)
            BETA(M)=REAL(BET)
         ENDIF
         IF(DELT.LT.EPS2) GO TO 370
      ELSE
         DO 365 IGR=1,NGRP
         DO 360 I=1,LL4
         EASS(I,IGR)=EASS(I,IGR)+GRAD1(I,IGR)
         GAR1(I,IGR)=GAR1(I,IGR)+GAR2(I,IGR)
         GRAD2(I,IGR)=GRAD1(I,IGR)
         GAR3(I,IGR)=GAR2(I,IGR)
  360    CONTINUE
  365    CONTINUE
         IF(IMPX.GE.2) WRITE(6,510) M,ALP,BET,ZNORM,0.0,ITER
*        COMPUTE THE CONVERGENCE HISTOGRAM.
         IF(IMPH.GE.1) THEN
            LMPH=IMPH.GE.1
            CALL FLDXCO(IPFLUP,LL4,NUN,EASS(1,NGRP),LMPH,ERR(M))
            ALPH(M)=REAL(ALP)
            BETA(M)=REAL(BET)
         ENDIF
      ENDIF
      IF(M.EQ.1) TEST=DELT
      IF((M.GT.20).AND.(DELT.GT.TEST)) CALL XABORT('GPTLIV: CONVERGENC'
     1 //'E FAILURE.')
      IF(M.GE.MAXX0) THEN
         WRITE(6,520)
         GO TO 370
      ENDIF
      IF(MOD(M,36).EQ.0) THEN
         ISTART=M+1
         NNADI=NNADI+1
         IF(IMPX.NE.0) WRITE(6,530) NNADI
      ENDIF
      GO TO 100
*----
*  SAVE THE CONVERGENCE HISTOGRAM ON LCM.
*----
  370 IF(IMPH.GE.2) THEN
         IGRAPH=0
  390    IGRAPH=IGRAPH+1
         WRITE(TEXT12,'(5HHISTO,I3)') IGRAPH
         CALL LCMLEN (IPFLUP,TEXT12,ILENG,ITYLCM)
         IF(ILENG.EQ.0) THEN
            CALL LCMSIX (IPFLUP,TEXT12,1)
            CALL LCMPTC (IPFLUP,'HTITLE',72,1,TITR)
            CALL LCMPUT (IPFLUP,'ALPHA',M,2,ALPH)
            CALL LCMPUT (IPFLUP,'BETA',M,2,BETA)
            CALL LCMPUT (IPFLUP,'ERROR',M,2,ERR)
            CALL LCMPUT (IPFLUP,'IMPH',1,1,IMPH)
            CALL LCMSIX (IPFLUP,' ',2)
         ELSE
            GO TO 390
         ENDIF
      ENDIF
      DEALLOCATE(WORK2,WORK1,GAR3,GAR2,GAR1,GRAD2,GRAD1)
      RETURN
*
  500 FORMAT (/29X,15HORTHONORMALIZA-/11X,5HALPHA,3X,4HBETA,6X,
     1 11HTION FACTOR,6X,8HACCURACY,5X,7HTHERMAL)
  510 FORMAT (1X,I3,4X,2F8.3,1P,E14.2,6X,E10.2,5X,1H(,I4,1H))
  520 FORMAT(/53H GPTLIV: ***WARNING*** THE MAXIMUM NUMBER OF OUTER IT,
     1 20HERATIONS IS REACHED.)
  530 FORMAT(/53H GPTLIV: INCREASING THE NUMBER OF INNER ITERATIONS TO,
     1 I3,36H ADI ITERATIONS PER OUTER ITERATION./)
      END
