*DECK FLDBHR
      SUBROUTINE FLDBHR(IPTRK,IPSYS,LADJ,LL4,ITY,NUN,NGRP,ICL1,ICL2,
     1 IMPX,NADI,MAXINR,EPSINR,ITER,TKT,TKB,GRAD1)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform thermal (up-scattering) iterations in Bivac.
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
* LADJ    flag set to .TRUE. for adjoint solution acceleration.
* LL4     order of the system matrices.
* ITY     type of solution (2: classical Bivac; 3: Thomas-Raviart).
* NUN     number of unknowns in each energy group.
* NGRP    number of energy groups.
* ICL1    number of free iretations in one cycle of the up-scattering
*         iterations.
* ICL2    number of accelerated up-scattering iterations in one cycle.
* IMPX    print parameter (set to 0 for no printing).
* NADI    number of inner ADI iterations per outer iteration (used with
*         SPN approximations).
* MAXINR  maximum number of thermal iterations.
* EPSINR  thermal iteration epsilon.
*
*Parameters: input/output
* ITER    actual number of thermal iterations.
* TKT     CPU time spent to compute the solution of linear systems.
* TKB     CPU time spent to compute the bilinear products.
* GRAD1   delta flux for this outer iteration.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS
      INTEGER LL4,ITY,NUN,NGRP,ICL1,ICL2,IMPX,NADI,MAXINR,ITER
      REAL EPSINR,TKT,TKB,GRAD1(NUN,NGRP)
      LOGICAL LADJ
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE)
      CHARACTER TEXT12*12,TEXT3*3
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, DIMENSION(:), ALLOCATABLE :: WORK2
      REAL, DIMENSION(:,:), ALLOCATABLE :: GAR2
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: WORK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      IF(MAXINR.EQ.0) RETURN
      ALLOCATE(GAR2(NUN,NGRP),WORK(LL4,NGRP,3),WORK2(LL4))
*
      CALL LCMGET(IPSYS,'STATE-VECTOR',ISTATE)
      NAN=ISTATE(8)
      NCTOT=ICL1+ICL2
      IF(ICL2.EQ.0) THEN
         NCPTM=NCTOT+1
      ELSE
         NCPTM=ICL1
      ENDIF
      DO 15 IGR=1,NGRP
      DO 10 I=1,LL4
      WORK(I,IGR,1)=0.0
      WORK(I,IGR,2)=0.0
      WORK(I,IGR,3)=GRAD1(I,IGR)
   10 CONTINUE
   15 CONTINUE
      IGDEB=1
*----
*  PERFORM THERMAL (UP-SCATTERING) ITERATIONS
*----
      TEXT3='NO '
      ITER=2
      DO
        CALL KDRCPU(TK1)
        IF(LADJ) THEN
*         ADJOINT SOLUTION
          DO 35 IGR=IGDEB,NGRP
          WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
          CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,WORK(1,IGR,3),
     1    GAR2(1,IGR))
          DO 30 JGR=1,NGRP
          IF(JGR.EQ.IGR) GO TO 30
          WRITE(TEXT12,'(1HA,2I3.3)') JGR,IGR
          CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
          IF(ILONG.EQ.0) GO TO 30
          CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,WORK(1,JGR,3),WORK2)
          DO 20 I=1,LL4
          GAR2(I,IGR)=GAR2(I,IGR)-WORK2(I)
   20     CONTINUE
   30     CONTINUE
   35     CONTINUE
          DO 65 IGR=NGRP,IGDEB,-1
          DO 50 JGR=NGRP,IGR+1,-1
          WRITE(TEXT12,'(1HA,2I3.3)') JGR,IGR
          CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
          IF(ILONG.EQ.0) GO TO 50
          CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GAR2(1,JGR),WORK2)
          DO 40 I=1,LL4
          GAR2(I,IGR)=GAR2(I,IGR)+WORK2(I)
   40     CONTINUE
   50     CONTINUE
          CALL KDRCPU(TK2)
          TKB=TKB+(TK2-TK1)
*
          CALL KDRCPU(TK1)
          WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
          IF(ITY.EQ.11) THEN
*           SIMPLIFIED PN BIVAC TRACKING.
            IF(NAN.EQ.0) CALL XABORT('FLDBHR: SPN-ONLY ALGORITHM.')
            CALL FLDBSS(TEXT12,IPTRK,IPSYS,LL4,NBMIX,NAN,GAR2(1,IGR),
     1      NADI)
          ELSE
            CALL MTLDLS(TEXT12,IPTRK,IPSYS,LL4,ITY,GAR2(1,IGR))
          ENDIF
          DO 60 I=1,LL4
          WORK(I,IGR,1)=WORK(I,IGR,2)
          WORK(I,IGR,2)=WORK(I,IGR,3)
          WORK(I,IGR,3)=GRAD1(I,IGR)+(WORK(I,IGR,2)-GAR2(I,IGR))
   60     CONTINUE
   65     CONTINUE
        ELSE
*         DIRECT SOLUTION
          DO 85 IGR=IGDEB,NGRP
          WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
          CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,WORK(1,IGR,3),
     1    GAR2(1,IGR))
          DO 80 JGR=1,NGRP
          IF(JGR.EQ.IGR) GO TO 80
          WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
          CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
          IF(ILONG.EQ.0) GO TO 80
          CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,WORK(1,JGR,3),WORK2)
          DO 70 I=1,LL4
          GAR2(I,IGR)=GAR2(I,IGR)-WORK2(I)
   70     CONTINUE
   80     CONTINUE
   85     CONTINUE
          DO 115 IGR=IGDEB,NGRP
          DO 100 JGR=1,IGR-1
          WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
          CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
          IF(ILONG.EQ.0) GO TO 100
          CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GAR2(1,JGR),WORK2)
          DO 90 I=1,LL4
          GAR2(I,IGR)=GAR2(I,IGR)+WORK2(I)
   90     CONTINUE
  100     CONTINUE
          CALL KDRCPU(TK2)
          TKB=TKB+(TK2-TK1)
*
          CALL KDRCPU(TK1)
          WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
          IF(ITY.EQ.11) THEN
*           SIMPLIFIED PN BIVAC TRACKING.
            IF(NAN.EQ.0) CALL XABORT('FLDBHR: SPN-ONLY ALGORITHM.')
            CALL FLDBSS(TEXT12,IPTRK,IPSYS,LL4,NBMIX,NAN,GAR2(1,IGR),
     1      NADI)
          ELSE
            CALL MTLDLS(TEXT12,IPTRK,IPSYS,LL4,ITY,GAR2(1,IGR))
          ENDIF
          DO 110 I=1,LL4
          WORK(I,IGR,1)=WORK(I,IGR,2)
          WORK(I,IGR,2)=WORK(I,IGR,3)
          WORK(I,IGR,3)=GRAD1(I,IGR)+(WORK(I,IGR,2)-GAR2(I,IGR))
  110     CONTINUE
  115     CONTINUE
        ENDIF
        IF(MOD(ITER-2,NCTOT).GE.NCPTM) THEN
          CALL FLD2AC(NGRP,LL4,IGDEB,WORK,ZMU)
        ELSE
          ZMU=1.0
        ENDIF
        IGDEBO=IGDEB
        DO 130 IGR=IGDEBO,NGRP
        GINN=0.0
        FINN=0.0
        DO 120 I=1,LL4
        GINN=MAX(GINN,ABS(WORK(I,IGR,2)-WORK(I,IGR,3)))
        FINN=MAX(FINN,ABS(WORK(I,IGR,3)))
  120   CONTINUE
        GINN=GINN/FINN
        IF((GINN.LT.EPSINR).AND.(IGDEB.EQ.IGR)) IGDEB=IGDEB+1
  130   CONTINUE
        CALL KDRCPU(TK2)
        TKT=TKT+(TK2-TK1)
        IF(GINN.LT.EPSINR) TEXT3='YES'
        IF(IMPX.GT.2) WRITE(6,1000) ITER,GINN,EPSINR,IGDEB,ZMU,TEXT3
        IF((GINN.LT.EPSINR).OR.(ITER.EQ.MAXINR)) GO TO 160
        ITER=ITER+1
      ENDDO
*----
*  END OF THERMAL ITERATIONS
*----
  160 DO 175 I=1,LL4
      DO 170 IGR=1,NGRP
      GRAD1(I,IGR)=WORK(I,IGR,3)
  170 CONTINUE
  175 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAR2,WORK,WORK2)
      RETURN
*
 1000 FORMAT (10X,3HIN(,I3,6H) FLX:,5H PRC=,1P,E9.2,5H TAR=,E9.2,
     1 7H IGDEB=, I13,6H ACCE=,0P,F12.5,12H  CONVERGED=,A3)
      END
