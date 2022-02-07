*DECK FLDDEF
      SUBROUTINE FLDDEF (MAX,IPTRK,IPSYS,LL4,ITY,NGRP,IMOD,LMOD,EVECT,
     1 ADECT,VEC,IADJ,VEA,VEB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multigroup Hotelling deflation procedure (1- multiplication of the
* 'A' matrix by a vector; 2- multiplication of a deflated 'B' matrix
* by the same vector).
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
* MAX     first dimension of arrays EVECT, ADECT, VEC, VEA and VEB.
* IPTRK   L_TRACK pointer to the tracking information.
* IPSYS   L_SYSTEM pointer to system matrices.
* LL4     order of the system matrices.
* ITY     type of solution (2: classical Trivac; 3: Thomas-Raviart).
* NGRP    number of energy groups.
* IMOD    number of the harmonic to be deflated.
* LMOD    total number of harmonics.
* EVECT   direct eigenvector.
* ADECT   adjoint eigenvector.
* VEC     vector to be multiplied.
* IADJ    type of deflation:
*         =1 for a direct deflation; =2 for an adjoint deflation.
*
*Parameters: output
* VEA     result of the multiplication to the 'A' matrix.
* VEB     result of the multiplication to the 'B' matrix.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS
      INTEGER MAX,LL4,ITY,NGRP,IMOD,LMOD,IADJ
      REAL EVECT(MAX,NGRP,LMOD),ADECT(MAX,NGRP,LMOD),VEC(MAX,NGRP),
     1 VEA(MAX,NGRP),VEB(MAX,NGRP)
*----
*  LOCAL VARIABLES
*----
      CHARACTER*12 TEXT12
      DOUBLE PRECISION DDELN1,DDELD1
      REAL, DIMENSION(:), ALLOCATABLE :: GAF,W
      REAL, DIMENSION(:), POINTER :: AGARM
      TYPE(C_PTR) AGARM_PTR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(GAF(LL4))
*
      IF(IADJ.EQ.1) THEN
*        DIRECT CASE.
         DO 45 IGR=1,NGRP
         CALL XDRSET(VEB(1,IGR),LL4,0.0)
         DO 40 JGR=1,NGRP
         WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 40
         CALL LCMGPD(IPSYS,TEXT12,AGARM_PTR)
         CALL C_F_POINTER(AGARM_PTR,AGARM,(/ ILONG /))
         DO 20 I=1,ILONG
         VEB(I,IGR)=VEB(I,IGR)+AGARM(I)*VEC(I,JGR)
   20    CONTINUE
   40    CONTINUE
   45    CONTINUE
         DO 132 JMOD=1,IMOD-1
         DDELN1=0.0D0
         DDELD1=0.0D0
         DO 125 IGR=1,NGRP
         WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
         CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,EVECT(1,IGR,JMOD),
     1   VEA(1,IGR))
         CALL XDRSET(GAF,LL4,0.0)
         DO 110 JGR=1,NGRP
         IF(JGR.EQ.IGR) GO TO 80
         WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 80
         IF(ITY.EQ.13) THEN
            ALLOCATE(W(LL4))
            CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,EVECT(1,JGR,JMOD),W)
            DO 50 I=1,LL4
            VEA(I,IGR)=VEA(I,IGR)-W(I)
   50       CONTINUE
            DEALLOCATE(W)
         ELSE
            CALL LCMGPD(IPSYS,TEXT12,AGARM_PTR)
            CALL C_F_POINTER(AGARM_PTR,AGARM,(/ ILONG /))
            DO 60 I=1,ILONG
            VEA(I,IGR)=VEA(I,IGR)-AGARM(I)*EVECT(I,JGR,JMOD)
   60       CONTINUE
         ENDIF
   80    WRITE(TEXT12,'(1HB,2I3.3)') JGR,IGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 110
         CALL LCMGPD(IPSYS,TEXT12,AGARM_PTR)
         CALL C_F_POINTER(AGARM_PTR,AGARM,(/ ILONG /))
         DO 90 I=1,ILONG
         GAF(I)=GAF(I)+AGARM(I)*ADECT(I,JGR,JMOD)
   90    CONTINUE
  110    CONTINUE
         DO 120 I=1,LL4
         DDELN1=DDELN1+GAF(I)*VEC(I,IGR)
         DDELD1=DDELD1+ADECT(I,IGR,JMOD)*VEA(I,IGR)
  120    CONTINUE
  125    CONTINUE
         DDELN1=DDELN1/DDELD1
         DO 131 IGR=1,NGRP
         DO 130 I=1,LL4
         VEB(I,IGR)=VEB(I,IGR)-VEA(I,IGR)*REAL(DDELN1)
  130    CONTINUE
  131    CONTINUE
  132    CONTINUE
         DO 165 IGR=1,NGRP
         WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
         CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,VEC(1,IGR),VEA(1,IGR))
         DO 160 JGR=1,NGRP
         IF(JGR.EQ.IGR) GO TO 160
         WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 160
         IF(ITY.EQ.13) THEN
            ALLOCATE(W(LL4))
            CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,VEC(1,JGR),W)
            DO 135 I=1,LL4
            VEA(I,IGR)=VEA(I,IGR)-W(I)
  135       CONTINUE
            DEALLOCATE(W)
         ELSE
            CALL LCMGPD(IPSYS,TEXT12,AGARM_PTR)
            CALL C_F_POINTER(AGARM_PTR,AGARM,(/ ILONG /))
            DO 140 I=1,ILONG
            VEA(I,IGR)=VEA(I,IGR)-AGARM(I)*VEC(I,JGR)
  140       CONTINUE
         ENDIF
  160    CONTINUE
  165    CONTINUE
      ELSE IF(IADJ.EQ.2) THEN
*        ADJOINT CASE.
         DO 205 IGR=1,NGRP
         CALL XDRSET(VEB(1,IGR),LL4,0.0)
         DO 200 JGR=1,NGRP
         WRITE(TEXT12,'(1HB,2I3.3)') JGR,IGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 200
         CALL LCMGPD(IPSYS,TEXT12,AGARM_PTR)
         CALL C_F_POINTER(AGARM_PTR,AGARM,(/ ILONG /))
         DO 180 I=1,ILONG
         VEB(I,IGR)=VEB(I,IGR)+AGARM(I)*VEC(I,JGR)
  180    CONTINUE
  200    CONTINUE
  205    CONTINUE
         DO 292 JMOD=1,IMOD-1
         DDELN1=0.0D0
         DDELD1=0.0D0
         DO 285 IGR=1,NGRP
         WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
         CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,ADECT(1,IGR,JMOD),
     1   VEA(1,IGR))
         CALL XDRSET(GAF,LL4,0.0)
         DO 270 JGR=1,NGRP
         IF(JGR.EQ.IGR) GO TO 240
         WRITE(TEXT12,'(1HA,2I3.3)') JGR,IGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 240
         IF(ITY.EQ.13) THEN
            ALLOCATE(W(LL4))
            CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,ADECT(1,JGR,JMOD),W)
            DO 210 I=1,LL4
            VEA(I,IGR)=VEA(I,IGR)-W(I)
  210       CONTINUE
            DEALLOCATE(W)
         ELSE
            CALL LCMGPD(IPSYS,TEXT12,AGARM_PTR)
            CALL C_F_POINTER(AGARM_PTR,AGARM,(/ ILONG /))
            DO 220 I=1,ILONG
            VEA(I,IGR)=VEA(I,IGR)-AGARM(I)*ADECT(I,JGR,JMOD)
  220       CONTINUE
         ENDIF
  240    WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 270
         CALL LCMGPD(IPSYS,TEXT12,AGARM_PTR)
         CALL C_F_POINTER(AGARM_PTR,AGARM,(/ ILONG /))
         DO 250 I=1,ILONG
         GAF(I)=GAF(I)+AGARM(I)*EVECT(I,JGR,JMOD)
  250    CONTINUE
  270    CONTINUE
         DO 280 I=1,LL4
         DDELN1=DDELN1+GAF(I)*VEC(I,IGR)
         DDELD1=DDELD1+EVECT(I,IGR,JMOD)*VEA(I,IGR)
  280    CONTINUE
  285    CONTINUE
         DDELN1=DDELN1/DDELD1
         DO 291 IGR=1,NGRP
         DO 290 I=1,LL4
         VEB(I,IGR)=VEB(I,IGR)-VEA(I,IGR)*REAL(DDELN1)
  290    CONTINUE
  291    CONTINUE
  292    CONTINUE
         DO 325 IGR=1,NGRP
         WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
         CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,VEC(1,IGR),VEA(1,IGR))
         DO 320 JGR=1,NGRP
         IF(JGR.EQ.IGR) GO TO 320
         WRITE(TEXT12,'(1HA,2I3.3)') JGR,IGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 320
         IF(ITY.EQ.13) THEN
            ALLOCATE(W(LL4))
            CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,VEC(1,JGR),W)
            DO 295 I=1,LL4
            VEA(I,IGR)=VEA(I,IGR)-W(I)
  295       CONTINUE
            DEALLOCATE(W)
         ELSE
            CALL LCMGPD(IPSYS,TEXT12,AGARM_PTR)
            CALL C_F_POINTER(AGARM_PTR,AGARM,(/ ILONG /))
            DO 300 I=1,ILONG
            VEA(I,IGR)=VEA(I,IGR)-AGARM(I)*VEC(I,JGR)
  300       CONTINUE
         ENDIF
  320    CONTINUE
  325    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAF)
      RETURN
      END
