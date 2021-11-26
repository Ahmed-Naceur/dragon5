*DECK GPTAFL
      SUBROUTINE GPTAFL (IPTRK,IPSYS0,IPFLUP,LL4,ITY,NUN,NGRP,ICL1,ICL2,
     1 NSTART,IMPX,IMPH,TITR,EPS2,MAXINR,EPSINR,NADI,MAXX0,FKEFF,EVECT,
     2 ADECT,FKEFF2,EASS,SOUR)
*
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of a multigroup fixed source eigenvalue problem for the
* calculation of an adjoint gpt solution in Trivac. use the precondi-
* tioned power method.
*
*Copyright:
* Copyright (C) 1987 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPTRK   L_TRACK pointer to the tracking information
* IPSYS0  L_SYSTEM pointer to unperturbed system matrices
* IPFLUP  L_FLUX pointer to the gpt solution
* LL4     order of the system matrices.
* ITY     type of solution (2: classical Trivac; 3: Thomas-Raviart).
* NUN     number of unknowns in each energy group.
* NGRP    number of energy groups.
* ICL1    number of free iterations in one cycle of the inverse power
*         method
* ICL2    number of accelerated iterations in one cycle
* NSTART  GMRES method flag. =0: use Livolant acceleration; 
*         >0: restarts the GMRES method every NSTART iterations.
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
* EPS2    convergence criteria for the flux
* MAXINR  maximum number of thermal iterations.
* EPSINR  thermal iteration epsilon.
* NADI    initial number of inner adi iterations per outer iteration
* MAXX0   maximum number of outer iterations
* FKEFF   effective multiplication factor
* EVECT   unknown vector for the non perturbed direct flux
* ADECT   unknown vector for the non perturbed adjoint flux
* SOUR    fixed source
*
*Parameters: output
* FKEFF2  perturbed effective multiplication factor
* EASS    converged generalized adjoint
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
      TYPE(C_PTR) IPTRK,IPSYS0,IPFLUP
      CHARACTER TITR*72,HSMG*131
      INTEGER LL4,ITY,NUN,NGRP,ICL1,ICL2,NSTART,IMPX,IMPH,MAXINR,NADI,
     1 MAXX0
      REAL EPS2,EPSINR,FKEFF,EVECT(NUN,NGRP),ADECT(NUN,NGRP),FKEFF2,
     1 EASS(NUN,NGRP),SOUR(NUN,NGRP)
*----
*  LOCAL VARIABLES
*----
      CHARACTER*12 TEXT12
      DOUBLE PRECISION AIL,BIL,EVAL,ZNORM,GAZ,DAZ
      REAL TKT,TKB
      REAL, DIMENSION(:), ALLOCATABLE :: WORK1,WORK3
      REAL, DIMENSION(:,:), ALLOCATABLE :: GRAD1,GAR1
      REAL, DIMENSION(:), POINTER :: AGAR
      TYPE(C_PTR) AGAR_PTR
      DATA EPS1/1.0E-4/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(GRAD1(NUN,NGRP),GAR1(NUN,NGRP),WORK1(NUN))
*
      CALL MTOPEN(IMPX,IPTRK,LL4)
      IF(LL4.GT.NUN) CALL XABORT('DELDFL: INVALID NUMBER OF UNKNOWNS.')
*----
*  UNPERTURBED EIGENVALUE CALCULATION.
*----
      AIL=0.0D0
      BIL=0.0D0
      DO 85 IGR=1,NGRP
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
      CALL MTLDLM(TEXT12,IPTRK,IPSYS0,LL4,ITY,EVECT(1,IGR),GRAD1(1,IGR))
      CALL XDRSET(WORK1,LL4,0.0)
      DO 70 JGR=1,NGRP
      IF(JGR.EQ.IGR) GO TO 40
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,JGR
      CALL LCMLEN(IPSYS0,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 40
      IF(ITY.EQ.13) THEN
         ALLOCATE(WORK3(LL4))
         CALL MTLDLM(TEXT12,IPTRK,IPSYS0,LL4,ITY,EVECT(1,JGR),WORK3)
         DO 20 I=1,LL4
         GRAD1(I,IGR)=GRAD1(I,IGR)-WORK3(I)
   20    CONTINUE
         DEALLOCATE(WORK3)
      ELSE
         CALL LCMGPD(IPSYS0,TEXT12,AGAR_PTR)
         CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
         DO 30 I=1,ILONG
         GRAD1(I,IGR)=GRAD1(I,IGR)-AGAR(I)*EVECT(I,JGR)
   30    CONTINUE
      ENDIF
*
   40 WRITE(TEXT12,'(1HB,2I3.3)') IGR,JGR
      CALL LCMLEN(IPSYS0,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 70
      CALL LCMGPD(IPSYS0,TEXT12,AGAR_PTR)
      CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
      DO 60 I=1,ILONG
      WORK1(I)=WORK1(I)+AGAR(I)*EVECT(I,JGR)
   60 CONTINUE
*
   70 CONTINUE
      DO 80 I=1,LL4
      AIL=AIL+ADECT(I,IGR)*GRAD1(I,IGR)
      BIL=BIL+ADECT(I,IGR)*WORK1(I)
   80 CONTINUE
   85 CONTINUE
      EVAL=AIL/BIL
      FKEFF2=REAL(1.0D0/EVAL)
      IF(ABS(FKEFF-1.0/EVAL).GT.EPS1) CALL XABORT('GPTAFL: THE COMPUTE'
     1 //'D AND PROVIDED K-EFFECTIVES ARE INCONSISTENTS.')
*----
*  VALIDATION OF THE FIXED SOURCE TERM.
*----
      AIL=0.0D0
      BIL=0.0D0
      DO 95 IGR=1,NGRP
      DO 90 I=1,LL4
      GAZ=EVECT(I,IGR)*SOUR(I,IGR)
      DAZ=EVECT(I,IGR)**2
      AIL=AIL+GAZ
      BIL=BIL+DAZ
   90 CONTINUE
   95 CONTINUE
      GAZ=ABS(AIL)/ABS(BIL)/REAL(LL4)
      IF(AIL.EQ.0.0) THEN
        CALL XDRSET(EASS,NUN*NGRP,0.0)
        FKEFF2=0.0
        DEALLOCATE(GRAD1,GAR1,WORK1)
        RETURN
      ENDIF
      IF(IMPX.GE.1) THEN
        WRITE(6,'(/28H GPTAFL: ORTHONORMALIZATION=,1P,E11.4)') GAZ
      ENDIF
      IF(GAZ.GT.EPS2) THEN
        WRITE(HSMG,'(46HGPTAFL: THE SOURCE TERM IS NOT ORTHOGONAL TO T,
     1  26HHE DIRECT REFERENCE FLUX (,1P,E11.4,2H).)') GAZ
        CALL XABORT(HSMG)
      ENDIF
*----
*  ORTHONORMALIZATION OF THE SOURCE TERM.
*----
      AIL=0.0D0
      BIL=0.0D0
      CALL XDRSET(GAR1,NUN*NGRP,0.0)
      DO 110 IGR=1,NGRP
      DO 100 JGR=1,NGRP
      WRITE(TEXT12,'(1HB,2I3.3)') JGR,IGR
      CALL LCMLEN(IPSYS0,TEXT12,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 100
      CALL LCMGPD(IPSYS0,TEXT12,AGAR_PTR)
      CALL C_F_POINTER(AGAR_PTR,AGAR,(/ NUN /))
      DO I=1,ILONG
        GAR1(I,IGR)=GAR1(I,IGR)+AGAR(I)*ADECT(I,JGR)
      ENDDO
  100 CONTINUE
      DO I=1,LL4
        AIL=AIL+EVECT(I,IGR)*SOUR(I,IGR)
        BIL=BIL+EVECT(I,IGR)*GAR1(I,IGR)
      ENDDO
  110 CONTINUE
      DO 125 IGR=1,NGRP
      DO 120 I=1,LL4
      SOUR(I,IGR)=SOUR(I,IGR)-REAL(AIL/BIL)*GAR1(I,IGR)
  120 CONTINUE
  125 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION.
*----
      DEALLOCATE(GRAD1,GAR1,WORK1)
*----
*  LIVOLANT ACCELERATION.
*----
      IF(IMPX.GE.1) WRITE (6,600) NADI
      IF(NSTART.EQ.0) THEN
        CALL GPTLIV(IPTRK,IPSYS0,IPFLUP,.TRUE.,LL4,ITY,NUN,NGRP,ICL1,
     1  ICL2,IMPX,IMPH,TITR,NADI,MAXINR,MAXX0,EPS2,EPSINR,EVAL,EVECT,
     2  ADECT,EASS,SOUR,TKT,TKB,ZNORM,M)
*----
*  GMRES.
*----
      ELSE IF(NSTART.GT.0) THEN
        CALL GPTMRA(IPTRK,IPSYS0,IPFLUP,.TRUE.,LL4,ITY,NUN,NGRP,ICL1,
     1  ICL2,IMPX,NADI,MAXINR,NSTART,MAXX0,EPS2,EPSINR,EVAL,EVECT,ADECT,
     2  EASS,SOUR,TKT,TKB,ZNORM,M)
      ENDIF
*----
*  SOLUTION EDITION.
*----
      IF(IMPX.GE.1) WRITE (6,610) M
      IF(IMPX.GE.3) THEN
         DO 130 IGR=1,NGRP
         WRITE (6,620) IGR,(EASS(I,IGR),I=1,LL4)
  130    CONTINUE
      ENDIF
      RETURN
*
  600 FORMAT(1H1/50H GPTAFL: ITERATIVE PROCEDURE BASED ON PRECONDITION,
     1 17HED POWER METHOD (,I2,37H ADI ITERATIONS PER OUTER ITERATION)./
     2 9X,40HADJOINT FIXED SOURCE EIGENVALUE PROBLEM.)
  610 FORMAT(/23H GPTAFL: CONVERGENCE IN,I4,12H ITERATIONS.)
  620 FORMAT(//52H GPTAFL: GENERALIZED ADJOINT CORRESPONDING TO THE GR,
     1 3HOUP,I4//(5X,1P,8E14.5))
      END
