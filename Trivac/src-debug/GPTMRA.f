*DECK GPTMRA
      SUBROUTINE GPTMRA(IPTRK,IPSYS,IPFLUP,LADJ,LL4,ITY,NUN,NGRP,ICL1,
     1 ICL2,IMPX,NADI,MAXINR,NSTART,MAXX0,EPS2,EPSINR,EVAL,EVECT,ADECT,
     2 EASS,SOUR,TKT,TKB,ZNORM,ITER)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of a multigroup fixed source eigenvalue problem for the
* calculation of a gpt solution in Trivac. Use the preconditioned power
* method with GMRES(m) acceleration.
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
* IMPX    print parameter (set to 0 for no printing).
* NADI    initial number of inner ADI iterations per outer iteration.
* MAXINR  maximum number of thermal iterations.
* NSTART  restarts the GMRES method every NSTART iterations.
* MAXX0   maximum number of outer iterations
* EPS2    outer iteration convergence criterion
* EPSINR  thermal iteration convergence criterion
* EVAL    eigenvalue.
* EVECT   unknown vector for the non perturbed direct flux
* ADECT   unknown vector for the non perturbed adjoint flux
* SOUR    fixed source
*
*Parameters: input/output
* EASS    solution of the fixed source eigenvalue problem
* TKT     CPU time spent to compute the solution of linear systems.
* TKB     CPU time spent to compute the bilinear products.
* ZNORM   Hotelling deflation accuracy.
* ITER    number of iterations.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS,IPFLUP
      LOGICAL LADJ
      INTEGER LL4,ITY,NUN,NGRP,ICL1,ICL2,IMPX,NADI,MAXINR,NSTART,MAXX0,
     1 ITER
      REAL EPS2,EPSINR,EVECT(NUN,NGRP),ADECT(NUN,NGRP),EASS(NUN*NGRP),
     1 SOUR(NUN,NGRP),TKT,TKB,SDOT
      DOUBLE PRECISION EVAL,ZNORM
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6)
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: RR,QQ,VV,GAR1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: V,H
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: G,C,S,X
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(V(NUN*NGRP,NSTART+1),G(NSTART+1),H(NSTART+1,NSTART+1),
     1 C(NSTART+1),S(NSTART+1),X(NUN*NGRP),GAR1(NUN*NGRP))
*----
*  GLOBAL GMRES ITERATION.
*----
      ALLOCATE(RR(NUN*NGRP),QQ(NUN*NGRP),VV(NUN*NGRP))
      
      EPS1=EPS2*SQRT(SDOT(NUN*NGRP,SOUR,1,SOUR,1))
      RHO=1.0E10
      ITER=0
      NITER=1
      NNADI=NADI
      DO WHILE((RHO.GT.EPS1).AND.(ITER.LT.MAXX0))
        CALL GPTGRA(IPTRK,IPSYS,IPFLUP,LADJ,.TRUE.,LL4,ITY,NUN,NGRP,
     1  ICL1,ICL2,IMPX,NNADI,MAXINR,EPSINR,EVAL,EVECT,ADECT,EASS,SOUR,
     2  GAR1,JTER0,TKT,TKB,ZNORM,RR)
        NITER=NITER+1
        DO I=1,NUN*NGRP
          X(I)=RR(I)
        ENDDO
        RHO=SQRT(DDOT(NUN*NGRP,X,1,X,1))
*----
*  TEST FOR TERMINATION ON ENTRY
*----
        IF(RHO.LT.EPS1) THEN
           DEALLOCATE(VV,QQ,RR)
           GO TO 100
        ENDIF
*
        CALL XDDSET(G,NSTART+1,0.0D0)
        CALL XDDSET(H,(NSTART+1)*NSTART,0.0D0)
        CALL XDDSET(V,NUN*NGRP*(NSTART+1),0.0D0)
        CALL XDDSET(C,NSTART+1,0.0D0)
        CALL XDDSET(S,NSTART+1,0.0D0)
        G(1)=RHO
        DO I=1,NUN*NGRP
          V(I,1)=X(I)/RHO
        ENDDO
*----
*  GMRES(1) ITERATION
*----
        K=0
        DO WHILE((RHO.GT.EPS1).AND.(K.LT.NSTART).AND.(ITER.LT.MAXX0))
          K=K+1
          ITER=ITER+1
          IF(IMPX.GT.1) WRITE(IUNOUT,300) ITER,RHO,JTER0
          DO I=1,NUN*NGRP
            VV(I)=REAL(V(I,K))
            QQ(I)=0.0
          ENDDO
          CALL GPTGRA(IPTRK,IPSYS,IPFLUP,LADJ,.TRUE.,LL4,ITY,NUN,NGRP,
     1    ICL1,ICL2,IMPX,NNADI,MAXINR,EPSINR,EVAL,EVECT,ADECT,VV,QQ,
     2    GAR1,JTER,TKT,TKB,ZNORM,RR)
          IF(JTER.NE.JTER0) CALL XABORT('GPTMRA: INCONSISTENT PRECONDIT'
     1    //'IONING IN GMRES.')
          NITER=NITER+1
          DO I=1,NUN*NGRP
            V(I,K+1)=-RR(I)
          ENDDO
*----
*  MODIFIED GRAM-SCHMIDT
*----
          DO J=1,K
            HR=DDOT(NUN*NGRP,V(1,J),1,V(1,K+1),1)
            H(J,K)=HR
            DO I=1,NUN*NGRP
              V(I,K+1)=V(I,K+1)-HR*V(I,J)
            ENDDO
          ENDDO
          H(K+1,K)=SQRT(DDOT(NUN*NGRP,V(1,K+1),1,V(1,K+1),1))
*----
*  REORTHOGONALIZE
*----
          DO J=1,K
            HR=DDOT(NUN*NGRP,V(1,J),1,V(1,K+1),1)
            H(J,K)=H(J,K)+HR
            DO I=1,NUN*NGRP
              V(I,K+1)=V(I,K+1)-HR*V(I,J)
            ENDDO
          ENDDO
          H(K+1,K)=SQRT(DDOT(NUN*NGRP,V(1,K+1),1,V(1,K+1),1))
*----
*  WATCH OUT FOR HAPPY BREAKDOWN 
*----
          IF(H(K+1,K).NE.0.0) THEN
            DO I=1,NUN*NGRP
              V(I,K+1)=V(I,K+1)/H(K+1,K)
            ENDDO
          ENDIF
*----
*  FORM AND STORE THE INFORMATION FOR THE NEW GIVENS ROTATION
*----
          DO I=1,K-1
            W1=C(I)*H(I,K)-S(I)*H(I+1,K)
            W2=S(I)*H(I,K)+C(I)*H(I+1,K)
            H(I,K)=W1
            H(I+1,K)=W2
          ENDDO
          ZNU=SQRT(H(K,K)**2+H(K+1,K)**2)
          IF(ZNU.NE.0.0) THEN
            C(K)=H(K,K)/ZNU
            S(K)=-H(K+1,K)/ZNU
            H(K,K)=C(K)*H(K,K)-S(K)*H(K+1,K)
            H(K+1,K)=0.0D0
            W1=C(K)*G(K)-S(K)*G(K+1)
            W2=S(K)*G(K)+C(K)*G(K+1)
            G(K)=W1
            G(K+1)=W2
          ENDIF
*----
*  UPDATE THE RESIDUAL NORM
*----
          RHO=ABS(G(K+1))
        ENDDO
*----
*  AT THIS POINT EITHER K > NSTART OR RHO < EPS1.
*  IT'S TIME TO COMPUTE X AND CYCLE.
*----
        DO J=1,K
          H(J,K+1)=G(J)
        ENDDO
        CALL ALSBD(K,1,H,IER,NSTART+1)
        IF(IER.NE.0) CALL XABORT('GPTMRA: SINGULAR MATRIX.')
        DO I=1,NUN*NGRP
          EASS(I)=EASS(I)+REAL(DDOT(K,V(I,1),NUN*NGRP,H(1,K+1),1))
        ENDDO
        IF(K.EQ.NSTART) THEN
          NNADI=NNADI+1
          IF(IMPX.NE.0) WRITE (6,310) NNADI
        ENDIF
      ENDDO
      DEALLOCATE(VV,QQ,RR)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  100 DEALLOCATE(GAR1,X,S,C,H,G,V)
      RETURN
*
  300 FORMAT(24H GPTMRA: OUTER ITERATION,I4,10H  L2 NORM=,1P,E11.4,
     1 28H (NB. OF THERMAL ITERATIONS=,I4,1H))
  310 FORMAT(/53H GPTMRA: INCREASING THE NUMBER OF INNER ITERATIONS TO,
     1 I3,36H ADI ITERATIONS PER OUTER ITERATION./)
      END
