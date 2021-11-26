*DECK FLDTHR
      SUBROUTINE FLDTHR(IPTRK,IPSYS,IPFLUX,LADJ,LL4,ITY,NUN,NGRP,ICL1,
     1 ICL2,IMPX,NADI,NSTARD,MAXINR,EPSINR,ITER,TKT,TKB,GRAD1)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform thermal (up-scattering) iterations in Trivac.
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
* LADJ    flag set to .TRUE. for adjoint solution acceleration.
* LL4     order of the system matrices.
* ITY     type of solution (2: classical Trivac; 3: Thomas-Raviart).
* NUN     number of unknowns in each energy group.
* NGRP    number of energy groups.
* ICL1    number of free iterations in one cycle of the inverse power
*         method.
* ICL2    number of accelerated iterations in one cycle.
* IMPX    print parameter (set to 0 for no printing).
* NADI    number of inner ADI iterations per outer iteration.
* NSTARD  number of restarting iterations with GMRES.
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
      TYPE(C_PTR) IPTRK,IPSYS,IPFLUX
      INTEGER LL4,ITY,NUN,NGRP,ICL1,ICL2,IMPX,NADI,NSTARD,MAXINR,ITER
      REAL EPSINR,TKT,TKB,GRAD1(NUN,NGRP)
      LOGICAL LADJ
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE)
      REAL(KIND=8) DERTOL
      CHARACTER TEXT12*12,TEXT3*3
      INTERFACE
        FUNCTION FLDONE(X,B,N,IPTRK,IPSYS,IPFLUX) RESULT(Y)
          USE GANLIB
          INTEGER, INTENT(IN) :: N
          REAL(KIND=8), DIMENSION(N), INTENT(IN) :: X, B
          REAL(KIND=8), DIMENSION(N) :: Y
          TYPE(C_PTR) IPTRK,IPSYS,IPFLUX
        END FUNCTION FLDONE
      END INTERFACE
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, DIMENSION(:), ALLOCATABLE :: W
      REAL, DIMENSION(:,:), ALLOCATABLE :: GAR2
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: WORK
      REAL, DIMENSION(:), POINTER :: AGAR
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: DWORK1,DWORK2
      TYPE(C_PTR) AGAR_PTR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      IF(MAXINR.EQ.0) RETURN
      ALLOCATE(GAR2(NUN,NGRP),WORK(LL4,NGRP,3))
*
      IF(NSTARD.GT.0) CALL LCMGET(IPFLUX,'STATE-VECTOR',ISTATE)
      NCTOT=ICL1+ICL2
      IF(ICL2.EQ.0) THEN
         NCPTM=NCTOT+1
      ELSE
         NCPTM=ICL1
      ENDIF
      DO 11 IGR=1,NGRP
      DO 10 I=1,LL4
      WORK(I,IGR,1)=0.0
      WORK(I,IGR,2)=0.0
      WORK(I,IGR,3)=GRAD1(I,IGR)
   10 CONTINUE
   11 CONTINUE
      IGDEB=1
*----
*  PERFORM THERMAL (UP-SCATTERING) ITERATIONS
*----
      TEXT3='NO '
      ITER=2
      DO
         CALL KDRCPU(TK1)
         IF(LADJ) THEN
*           ADJOINT SOLUTION
            DO 31 IGR=IGDEB,NGRP
            WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
            CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,WORK(1,IGR,3),
     1      GAR2(1,IGR))
            DO 30 JGR=1,NGRP
            IF(JGR.EQ.IGR) GO TO 30
            WRITE(TEXT12,'(1HA,2I3.3)') JGR,IGR
            CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
            IF(ILONG.EQ.0) GO TO 30
            IF(ITY.EQ.13) THEN
               ALLOCATE(W(LL4))
               CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,WORK(1,JGR,3),W)
               DO 15 I=1,LL4
               GAR2(I,IGR)=GAR2(I,IGR)-W(I)
   15          CONTINUE
               DEALLOCATE(W)
            ELSE
               CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
               CALL C_F_POINTER(AGAR_PTR,AGAR,(/ ILONG /))
               DO 20 I=1,ILONG
               GAR2(I,IGR)=GAR2(I,IGR)-AGAR(I)*WORK(I,JGR,3)
   20          CONTINUE
            ENDIF
   30       CONTINUE
   31       CONTINUE
            DO 61 IGR=NGRP,IGDEB,-1
            DO 50 JGR=NGRP,IGR+1,-1
            WRITE(TEXT12,'(1HA,2I3.3)') JGR,IGR
            CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
            IF(ILONG.EQ.0) GO TO 50
            IF(ITY.EQ.13) THEN
               ALLOCATE(W(LL4))
               CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GAR2(1,JGR),W)
               DO 35 I=1,LL4
               GAR2(I,IGR)=GAR2(I,IGR)+W(I)
   35          CONTINUE
               DEALLOCATE(W)
            ELSE
               CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
               CALL C_F_POINTER(AGAR_PTR,AGAR,(/ ILONG /))
               DO 40 I=1,ILONG
               GAR2(I,IGR)=GAR2(I,IGR)+AGAR(I)*GAR2(I,JGR)
   40          CONTINUE
            ENDIF
   50       CONTINUE
            CALL KDRCPU(TK2)
            TKB=TKB+(TK2-TK1)
*
            CALL KDRCPU(TK1)
            IF(NSTARD.EQ.0) THEN
               WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
               CALL FLDADI(TEXT12,IPTRK,IPSYS,LL4,ITY,GAR2(1,IGR),NADI)
               JTER=NADI
            ELSE
*              use a GMRES solution of the linear system
               DERTOL=EPSINR
               ISTATE(39)=IGR
               CALL LCMPUT(IPFLUX,'STATE-VECTOR',NSTATE,1,ISTATE)
               ALLOCATE(DWORK1(LL4),DWORK2(LL4))
               DWORK1(:LL4)=GAR2(:LL4,IGR)   ! source
               DWORK2(:LL4)=WORK(:LL4,IGR,3) ! estimate of the flux
               CALL FLDMRA(DWORK1,FLDONE,LL4,DERTOL,NSTARD,NADI,IMPX,
     1         IPTRK,IPSYS,IPFLUX,DWORK2,JTER)
               GAR2(:LL4,IGR)=REAL(DWORK2(:LL4))
               DEALLOCATE(DWORK2,DWORK1)
            ENDIF
            DO 60 I=1,LL4
            WORK(I,IGR,1)=WORK(I,IGR,2)
            WORK(I,IGR,2)=WORK(I,IGR,3)
            WORK(I,IGR,3)=GRAD1(I,IGR)+(WORK(I,IGR,2)-GAR2(I,IGR))
   60       CONTINUE
   61       CONTINUE
         ELSE
*           DIRECT SOLUTION
            DO 81 IGR=IGDEB,NGRP
            WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
            CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,WORK(1,IGR,3),
     1      GAR2(1,IGR))
            DO 80 JGR=1,NGRP
            IF(JGR.EQ.IGR) GO TO 80
            WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
            CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
            IF(ILONG.EQ.0) GO TO 80
            IF(ITY.EQ.13) THEN
               ALLOCATE(W(LL4))
               CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,WORK(1,JGR,3),W)
               DO 65 I=1,LL4
               GAR2(I,IGR)=GAR2(I,IGR)-W(I)
   65          CONTINUE
               DEALLOCATE(W)
            ELSE
               CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
               CALL C_F_POINTER(AGAR_PTR,AGAR,(/ ILONG /))
               DO 70 I=1,ILONG
               GAR2(I,IGR)=GAR2(I,IGR)-AGAR(I)*WORK(I,JGR,3)
   70          CONTINUE
            ENDIF
   80       CONTINUE
   81       CONTINUE
            DO 115 IGR=IGDEB,NGRP
            DO 100 JGR=1,IGR-1
            WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
            CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
            IF(ILONG.EQ.0) GO TO 100
            IF(ITY.EQ.13) THEN
               ALLOCATE(W(LL4))
               CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,GAR2(1,JGR),W)
               DO 85 I=1,LL4
               GAR2(I,IGR)=GAR2(I,IGR)+W(I)
   85          CONTINUE
               DEALLOCATE(W)
            ELSE
               CALL LCMGPD(IPSYS,TEXT12,AGAR_PTR)
               CALL C_F_POINTER(AGAR_PTR,AGAR,(/ ILONG /))
               DO 90 I=1,ILONG
               GAR2(I,IGR)=GAR2(I,IGR)+AGAR(I)*GAR2(I,JGR)
   90          CONTINUE
            ENDIF
  100       CONTINUE
            CALL KDRCPU(TK2)
            TKB=TKB+(TK2-TK1)
*
            CALL KDRCPU(TK1)
            WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
            IF(NSTARD.EQ.0) THEN
               WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
               CALL FLDADI(TEXT12,IPTRK,IPSYS,LL4,ITY,GAR2(1,IGR),NADI)
               JTER=NADI
            ELSE
*              use a GMRES solution of the linear system
               DERTOL=EPSINR
               ISTATE(39)=IGR
               CALL LCMPUT(IPFLUX,'STATE-VECTOR',NSTATE,1,ISTATE)
               ALLOCATE(DWORK1(LL4),DWORK2(LL4))
               DWORK1(:LL4)=GAR2(:LL4,IGR)   ! source
               DWORK2(:LL4)=WORK(:LL4,IGR,3) ! estimate of the flux
               CALL FLDMRA(DWORK1,FLDONE,LL4,DERTOL,NSTARD,NADI,IMPX,
     1         IPTRK,IPSYS,IPFLUX,DWORK2,JTER)
               GAR2(:LL4,IGR)=REAL(DWORK2(:LL4))
               DEALLOCATE(DWORK2,DWORK1)
            ENDIF
            DO 110 I=1,LL4
            WORK(I,IGR,1)=WORK(I,IGR,2)
            WORK(I,IGR,2)=WORK(I,IGR,3)
            WORK(I,IGR,3)=GRAD1(I,IGR)+(WORK(I,IGR,2)-GAR2(I,IGR))
  110       CONTINUE
  115       CONTINUE
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
  120    CONTINUE
         GINN=GINN/FINN
         IF((GINN.LT.EPSINR).AND.(IGDEB.EQ.IGR)) IGDEB=IGDEB+1
  130    CONTINUE
         CALL KDRCPU(TK2)
         TKT=TKT+(TK2-TK1)
         IF(GINN.LT.EPSINR) TEXT3='YES'
         IF(IMPX.GT.2) WRITE(6,1000) ITER,GINN,EPSINR,IGDEB,ZMU,TEXT3,
     1   JTER
         IF((GINN.LT.EPSINR).OR.(ITER.EQ.MAXINR)) EXIT
         ITER=ITER+1
      ENDDO
*----
*  END OF THERMAL ITERATIONS
*----
      DO 175 I=1,LL4
      DO 170 IGR=1,NGRP
      GRAD1(I,IGR)=WORK(I,IGR,3)
  170 CONTINUE
  175 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAR2,WORK)
      RETURN
*
 1000 FORMAT (10X,3HIN(,I3,6H) FLX:,5H PRC=,1P,E9.2,5H TAR=,E9.2,
     1 7H IGDEB=, I13,6H ACCE=,0P,F12.5,12H  CONVERGED=,A3,6H JTER=,
     2 I4)
      END
