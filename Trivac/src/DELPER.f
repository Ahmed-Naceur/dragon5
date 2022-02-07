*DECK DELPER
      SUBROUTINE DELPER (IPTRK,IPSYS0,IPSYSP,ADJ,NUN,NGRP,FKEFF,IMPX,
     1 EVECT,ADECT,DELKEF,SOUR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the source term for a direct or adjoint fixed source
* eigenvalue problem.
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
* IPSYS0  L_SYSTEM pointer to unperturbed system matrices.
* IPSYSP  L_SYSTEM pointer to delta system matrices.
* ADJ     adjoint flag. If ADJ=.true., we compute the source term for an
*         adjoint fixed source eigenvalue problem.
* NUN     total number of unknowns per energy group.
* NGRP    number of energy groups.
* FKEFF   reference k-effective.
* IMPX    delta k-effective is printed if impx.ge.1.
* EVECT   reference solution of the associated direct eigenvalue
*         problem.
* ADECT   reference solution of the associated adjoint eigenvalue
*         problem.
*
*Parameters: output
* DELKEF  delta k-effective.
* SOUR    fixed source term.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS0,IPSYSP
      INTEGER NUN,NGRP,IMPX
      LOGICAL ADJ
      REAL FKEFF,EVECT(NUN,NGRP),ADECT(NUN,NGRP),DELKEF,SOUR(NUN,NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      PARAMETER (EPS1=1.0E-4)
      INTEGER ISTATE(NSTATE)
      CHARACTER*12 TEXT12
      DOUBLE PRECISION AIL,BIL,EVAL,DEVAL
      REAL, DIMENSION(:), ALLOCATABLE :: WORK,WORK1
      REAL, DIMENSION(:), POINTER :: AGAR
      TYPE(C_PTR) AGAR_PTR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK(NUN))
*
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      LL4=ISTATE(11)
      NLF=ISTATE(30)
      IF(NLF.GT.0) LL4=LL4*NLF/2
      ITY=2
      IF(ISTATE(12).EQ.2) ITY=3
      IF((NLF.GT.0).AND.(ITY.GE.3)) ITY=10+ITY
      CALL MTOPEN(IMPX,IPTRK,LL4)
      IF(LL4.GT.NUN) CALL XABORT('DELPER: INVALID NUMBER OF UNKNOWNS.')
*----
*  COMPUTE THE NON-PERTURBED K-EFFECTIVE.
*----
      AIL=0.0D0
      BIL=0.0D0
      DO 85 IGR=1,NGRP
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
      CALL MTLDLM(TEXT12,IPTRK,IPSYS0,LL4,ITY,EVECT(1,IGR),SOUR(1,IGR))
      DO 10 I=1,LL4
      WORK(I)=0.0
   10 CONTINUE
      DO 70 JGR=1,NGRP
      IF(JGR.EQ.IGR) GO TO 40
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
      CALL LCMLEN(IPSYS0,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 40
      IF(ITY.EQ.13) THEN
         ALLOCATE(WORK1(LL4))
         CALL MTLDLM(TEXT12,IPTRK,IPSYS0,LL4,ITY,EVECT(1,JGR),WORK1)
         DO 20 I=1,LL4
         SOUR(I,IGR)=SOUR(I,IGR)-WORK1(I)
   20    CONTINUE
         DEALLOCATE(WORK1)
      ELSE
         CALL LCMGPD(IPSYS0,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
         DO 30 I=1,ILONG
         SOUR(I,IGR)=SOUR(I,IGR)-AGAR(I)*EVECT(I,JGR)
   30    CONTINUE
      ENDIF
   40 WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
      CALL LCMLEN(IPSYS0,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 70
      CALL LCMGPD(IPSYS0,TEXT12,AGAR_PTR)
      CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
      DO 50 I=1,ILONG
      WORK(I)=WORK(I)+AGAR(I)*EVECT(I,JGR)
   50 CONTINUE
   70 CONTINUE
      DO 80 I=1,LL4
      AIL=AIL+ADECT(I,IGR)*SOUR(I,IGR)
      BIL=BIL+ADECT(I,IGR)*WORK(I)
   80 CONTINUE
   85 CONTINUE
      EVAL=AIL/BIL
      IF(ABS(FKEFF-1.0/EVAL).GT.EPS1) CALL XABORT('DELPER: INCOMPATIBIL'
     1 //'ITY BETWEEN THE PROVIDED AND CALCULATED KEFF.')
*----
*  COMPUTE THE DIRECT OR ADJOINT SOURCE TERM.
*----
      IF(ADJ) THEN
         DO 155 IGR=1,NGRP
         WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
         CALL MTLDLM(TEXT12,IPTRK,IPSYSP,LL4,ITY,ADECT(1,IGR),
     1   SOUR(1,IGR))
         DO 150 JGR=1,NGRP
         IF(JGR.EQ.IGR) GO TO 120
         WRITE(TEXT12,'(1HA,2I3.3)') JGR,IGR
         CALL LCMLEN(IPSYSP,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 120
         IF(ITY.EQ.13) THEN
            ALLOCATE(WORK1(LL4))
            CALL MTLDLM(TEXT12,IPTRK,IPSYSP,LL4,ITY,ADECT(1,JGR),WORK1)
            DO 100 I=1,LL4
            SOUR(I,IGR)=SOUR(I,IGR)-WORK1(I)
  100       CONTINUE
            DEALLOCATE(WORK1)
         ELSE
            CALL LCMGPD(IPSYSP,TEXT12,AGAR_PTR)
            CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
            DO 110 I=1,ILONG
            SOUR(I,IGR)=SOUR(I,IGR)-AGAR(I)*ADECT(I,JGR)
  110       CONTINUE
         ENDIF
  120    WRITE(TEXT12,'(1HB,2I3.3)') JGR,IGR
         CALL LCMLEN(IPSYSP,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 150
         CALL LCMGPD(IPSYSP,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
         DO 130 I=1,ILONG
         SOUR(I,IGR)=SOUR(I,IGR)-REAL(EVAL)*AGAR(I)*ADECT(I,JGR)
  130    CONTINUE
  150    CONTINUE
  155    CONTINUE
         AIL=0.0D0
         DO 165 IGR=1,NGRP
         DO 160 I=1,LL4
         AIL=AIL+SOUR(I,IGR)*EVECT(I,IGR)
  160    CONTINUE
  165    CONTINUE
         DEVAL=AIL/BIL
         DO 215 IGR=1,NGRP
         DO 170 I=1,LL4
         WORK(I)=0.0
  170    CONTINUE
         DO 200 JGR=1,NGRP
         WRITE(TEXT12,'(1HB,2I3.3)') JGR,IGR
         CALL LCMLEN(IPSYS0,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 200
         CALL LCMGPD(IPSYS0,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
         DO 180 I=1,ILONG
         WORK(I)=WORK(I)+AGAR(I)*ADECT(I,JGR)
  180    CONTINUE
  200    CONTINUE
         DO 210 I=1,LL4
         SOUR(I,IGR)=SOUR(I,IGR)-REAL(DEVAL)*WORK(I)
  210    CONTINUE
  215    CONTINUE
      ELSE
         DO 285 IGR=1,NGRP
         WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
         CALL MTLDLM(TEXT12,IPTRK,IPSYSP,LL4,ITY,EVECT(1,IGR),
     1   SOUR(1,IGR))
         DO 280 JGR=1,NGRP
         IF(JGR.EQ.IGR) GO TO 250
         WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
         CALL LCMLEN(IPSYSP,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 250
         IF(ITY.EQ.13) THEN
            ALLOCATE(WORK1(LL4))
            CALL MTLDLM(TEXT12,IPTRK,IPSYSP,LL4,ITY,EVECT(1,JGR),WORK1)
            DO 220 I=1,LL4
            SOUR(I,IGR)=SOUR(I,IGR)-WORK1(I)
  220       CONTINUE
            DEALLOCATE(WORK1)
         ELSE
            CALL LCMGPD(IPSYSP,TEXT12,AGAR_PTR)
            CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
            DO 230 I=1,ILONG
            SOUR(I,IGR)=SOUR(I,IGR)-AGAR(I)*EVECT(I,JGR)
  230       CONTINUE
         ENDIF
  250    WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
         CALL LCMLEN(IPSYSP,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 280
         CALL LCMGPD(IPSYSP,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
         DO 260 I=1,ILONG
         SOUR(I,IGR)=SOUR(I,IGR)-REAL(EVAL)*AGAR(I)*EVECT(I,JGR)
  260    CONTINUE
  280    CONTINUE
  285    CONTINUE
         AIL=0.0D0
         DO 295 IGR=1,NGRP
         DO 290 I=1,LL4
         AIL=AIL+ADECT(I,IGR)*SOUR(I,IGR)
  290    CONTINUE
  295    CONTINUE
         DEVAL=AIL/BIL
         DO 345 IGR=1,NGRP
         DO 300 I=1,LL4
         WORK(I)=0.0
  300    CONTINUE
         DO 330 JGR=1,NGRP
         WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
         CALL LCMLEN(IPSYS0,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 330
         CALL LCMGPD(IPSYS0,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
         DO 310 I=1,ILONG
         WORK(I)=WORK(I)+AGAR(I)*EVECT(I,JGR)
  310    CONTINUE
  330    CONTINUE
         DO 340 I=1,LL4
         SOUR(I,IGR)=SOUR(I,IGR)-REAL(DEVAL)*WORK(I)
  340    CONTINUE
  345    CONTINUE
      ENDIF
      DELKEF=-REAL(DEVAL/(EVAL*EVAL))
      IF(IMPX.GE.1) WRITE (6,'(/21H DELPER: DELTA KEFF =,1P,E17.9/)')
     1 DELKEF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK)
      RETURN
      END
