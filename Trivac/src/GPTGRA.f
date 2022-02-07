*DECK GPTGRA
      SUBROUTINE GPTGRA(IPTRK,IPSYS,IPFLUP,LADJ,LGAR1,LL4,ITY,NUN,NGRP,
     1 ICL1,ICL2,IMPX,NNADI,MAXINR,EPSINR,EVAL,EVECT,ADECT,EASS,SOUR,
     2 GAR1,ITER,TKT,TKB,ZNORM,GRAD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute multigroup delta flux in a fixed source eigenvalue iteration.
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
* LGAR1   flag set to .TRUE. for recomputing GAR1.
* LL4     order of the system matrices.
* ITY     type of solution (2: classical Trivac; 3: Thomas-Raviart).
* NUN     number of unknowns in each energy group.
* NGRP    number of energy groups.
* ICL1    number of free up-scattering iterations in one cycle of the
*         inverse power method.
* ICL2    number of accelerated up-scattering iterations in one cycle.
* IMPX    print parameter (set to 0 for no printing).
* NNADI   number of inner ADI iterations per outer iteration.
* MAXINR  maximum number of thermal iterations.
* EPSINR  thermal iteration epsilon.
* EVAL    eigenvalue.
* EVECT   unknown vector for the non perturbed direct flux
* ADECT   unknown vector for the non perturbed adjoint flux
* EASS    solution of the fixed source eigenvalue problem
* SOUR    fixed source
* GAR1    delta flux for this iteration before Hotelling deflation.
*
*Parameters: input/output
* ITER    actual number of thermal iterations.
* TKT     CPU time spent to compute the solution of linear systems.
* TKB     CPU time spent to compute the bilinear products.
* ZNORM   Hotelling deflation accuracy.
* GRAD    delta flux for this iteration.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS,IPFLUP
      LOGICAL LADJ,LGAR1
      INTEGER LL4,ITY,NUN,NGRP,ICL1,ICL2,IMPX,NNADI,MAXINR,ITER
      REAL EPSINR,EVECT(NUN,NGRP),ADECT(NUN,NGRP),EASS(NUN,NGRP),
     1 SOUR(NUN,NGRP),GAR1(NUN,NGRP),TKT,TKB,GRAD(NUN,NGRP)
      DOUBLE PRECISION EVAL,ZNORM
*----
*  LOCAL VARIABLES
*----
      CHARACTER*12 TEXT12
      DOUBLE PRECISION DDELN1,DDELD1
      REAL, DIMENSION(:), ALLOCATABLE :: WORK1,WORK3
      REAL, DIMENSION(:), POINTER :: AGAR
      TYPE(C_PTR) AGAR_PTR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK1(NUN))
*
      IF(LADJ) THEN
         CALL KDRCPU(TK1)
*        ADJOINT SOLUTION
         IF(LGAR1) THEN
            DO 55 IGR=1,NGRP
            WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
            CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,EASS(1,IGR),
     1      GAR1(1,IGR))
            DO 50 JGR=1,NGRP
            IF(JGR.EQ.IGR) GO TO 30
            WRITE(TEXT12,'(1HA,2I3.3)') JGR,IGR
            CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
            IF(ILONG.EQ.0) GO TO 30
            IF(ITY.EQ.13) THEN
               ALLOCATE(WORK3(LL4))
               CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,EASS(1,JGR),WORK3)
               DO 10 I=1,LL4
               GAR1(I,IGR)=GAR1(I,IGR)-WORK3(I)
   10          CONTINUE
               DEALLOCATE(WORK3)
            ELSE
               CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
               CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
               DO 20 I=1,ILONG
               GAR1(I,IGR)=GAR1(I,IGR)-AGAR(I)*EASS(I,JGR)
   20          CONTINUE
            ENDIF
   30       WRITE(TEXT12,'(1HB,2I3.3)') JGR,IGR
            CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
            IF(ILONG.EQ.0) GO TO 50
            CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
            CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
            DO 40 I=1,ILONG
            GAR1(I,IGR)=GAR1(I,IGR)-REAL(EVAL)*AGAR(I)*EASS(I,JGR)
   40       CONTINUE
   50       CONTINUE
   55       CONTINUE
         ENDIF
*----
*  DIRECTION EVALUATION.
*----
         DO 100 IGR=NGRP,1,-1
         DO 60 I=1,LL4
         GRAD(I,IGR)=-SOUR(I,IGR)-GAR1(I,IGR)
   60    CONTINUE
         DO 90 JGR=NGRP,IGR+1,-1
         WRITE(TEXT12,'(1HA,2I3.3)') JGR,IGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 90
         IF(ITY.EQ.13) THEN
            CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD(1,JGR),WORK1)
            DO 70 I=1,LL4
            GRAD(I,IGR)=GRAD(I,IGR)+WORK1(I)
   70       CONTINUE
         ELSE
            CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
            CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
            DO 80 I=1,ILONG
            GRAD(I,IGR)=GRAD(I,IGR)+AGAR(I)*GRAD(I,JGR)
   80       CONTINUE
         ENDIF
   90    CONTINUE
*
         WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
         CALL FLDADI(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD(1,IGR),NNADI)
  100    CONTINUE
         CALL KDRCPU(TK2)
         TKB=TKB+(TK2-TK1)
*----
*  PERFORM THERMAL (UP-SCATTERING) ITERATIONS
*----
         ITER=1
         IF(MAXINR.GT.1) THEN
            CALL FLDTHR(IPTRK,IPSYS,IPFLUP,.TRUE.,LL4,ITY,NUN,NGRP,
     1      ICL1,ICL2,IMPX,NNADI,0,MAXINR,EPSINR,ITER,TKT,TKB,GRAD)
         ENDIF
*----
*  HOTELLING DEFLATION.
*----
         CALL KDRCPU(TK1)
         DDELN1=0.0D0
         DDELD1=0.0D0
         DO 135 IGR=1,NGRP
         CALL XDRSET(WORK1,LL4,0.0)
         DO 120 JGR=1,NGRP
         WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 120
         CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
         DO 110 I=1,ILONG
         WORK1(I)=WORK1(I)+AGAR(I)*EVECT(I,JGR)
  110    CONTINUE
  120    CONTINUE
         DO 130 I=1,LL4
         DDELN1=DDELN1+WORK1(I)*EASS(I,IGR)
         DDELD1=DDELD1+WORK1(I)*ADECT(I,IGR)
  130    CONTINUE
  135    CONTINUE
         ZNORM=DDELN1/DDELD1
         DO 145 IGR=1,NGRP
         DO 140 I=1,LL4
         GRAD(I,IGR)=GRAD(I,IGR)-REAL(ZNORM)*ADECT(I,IGR)
  140    CONTINUE
  145    CONTINUE
         CALL KDRCPU(TK2)
         TKB=TKB+(TK2-TK1)
      ELSE
         CALL KDRCPU(TK1)
*        DIRECT SOLUTION
         IF(LGAR1) THEN
            DO 195 IGR=1,NGRP
            WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
            CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,EASS(1,IGR),
     1      GAR1(1,IGR))
            DO 190 JGR=1,NGRP
            IF(JGR.EQ.IGR) GO TO 170
            WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
            CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
            IF(ILONG.EQ.0) GO TO 170
            IF(ITY.EQ.13) THEN
               ALLOCATE(WORK3(LL4))
               CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,EASS(1,JGR),WORK3)
               DO 150 I=1,LL4
               GAR1(I,IGR)=GAR1(I,IGR)-WORK3(I)
  150          CONTINUE
               DEALLOCATE(WORK3)
            ELSE
               CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
               CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
               DO 160 I=1,ILONG
               GAR1(I,IGR)=GAR1(I,IGR)-AGAR(I)*EASS(I,JGR)
  160          CONTINUE
            ENDIF
  170       WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
            CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
            IF(ILONG.EQ.0) GO TO 190
            CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
            CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
            DO 180 I=1,ILONG
            GAR1(I,IGR)=GAR1(I,IGR)-REAL(EVAL)*AGAR(I)*EASS(I,JGR)
  180       CONTINUE
  190       CONTINUE
  195       CONTINUE
         ENDIF
*----
*  DIRECTION EVALUATION.
*----
         DO 240 IGR=1,NGRP
         DO 200 I=1,LL4
         GRAD(I,IGR)=-SOUR(I,IGR)-GAR1(I,IGR)
  200    CONTINUE
         DO 230 JGR=1,IGR-1
         WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 230
         IF(ITY.EQ.13) THEN
            CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD(1,JGR),WORK1)
            DO 210 I=1,LL4
            GRAD(I,IGR)=GRAD(I,IGR)+WORK1(I)
  210       CONTINUE
         ELSE
            CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
            CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
            DO 220 I=1,ILONG
            GRAD(I,IGR)=GRAD(I,IGR)+AGAR(I)*GRAD(I,JGR)
  220       CONTINUE
         ENDIF
  230    CONTINUE
*
         WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
         CALL FLDADI(TEXT12,IPTRK,IPSYS,LL4,ITY,GRAD(1,IGR),NNADI)
  240    CONTINUE
         CALL KDRCPU(TK2)
         TKB=TKB+(TK2-TK1)
*----
*  PERFORM THERMAL (UP-SCATTERING) ITERATIONS
*----
         ITER=1
         IF(MAXINR.GT.1) THEN
            CALL FLDTHR(IPTRK,IPSYS,IPFLUP,.FALSE.,LL4,ITY,NUN,NGRP,
     1      ICL1,ICL2,IMPX,NNADI,0,MAXINR,EPSINR,ITER,TKT,TKB,GRAD)
         ENDIF
*----
*  HOTELLING DEFLATION.
*----
         CALL KDRCPU(TK1)
         DDELN1=0.0D0
         DDELD1=0.0D0
         DO 275 IGR=1,NGRP
         CALL XDRSET(WORK1,LL4,0.0)
         DO 260 JGR=1,NGRP
         WRITE(TEXT12,'(1HB,2I3.3)') JGR,IGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 260
         CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
         DO 250 I=1,ILONG
         WORK1(I)=WORK1(I)+AGAR(I)*ADECT(I,JGR)
  250    CONTINUE
  260    CONTINUE
         DO 270 I=1,LL4
         DDELN1=DDELN1+WORK1(I)*EASS(I,IGR)
         DDELD1=DDELD1+WORK1(I)*EVECT(I,IGR)
  270    CONTINUE
  275    CONTINUE
         ZNORM=DDELN1/DDELD1
         DO 285 IGR=1,NGRP
         DO 280 I=1,LL4
         GRAD(I,IGR)=GRAD(I,IGR)-REAL(ZNORM)*EVECT(I,IGR)
  280    CONTINUE
  285    CONTINUE
         CALL KDRCPU(TK2)
         TKB=TKB+(TK2-TK1)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK1)
      RETURN
      END
