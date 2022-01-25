*DECK PKIDRV
      SUBROUTINE PKIDRV(IPMAP,NALPHA,NGROUP,LAMBDA,EPSILON,BETAI,
     1 LAMBDAI,DT,PARAMI,PARAMB,T,Y)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of the point kinetic equations using the Runge-Kutta method.
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPMAP   pointer to the point kinetic directory
* NALPHA  number of feedback parameters
* NGROUP  number of delayed precursor groups
* LAMBDA  prompt neutron generation time
* EPSILON Runge-Kutta epsilon
* BETAI   delayed neutron fraction vector
* LAMBDAI delayed neutron time constant vector
* DT      stage duration (double precision value)
* PARAMI  initial values of the global parameters corresponding to
*         RHO=0
* PARAMB  values of global parameters at beginning of stage
* T       time at beggining of stage (double precision value)
* Y       solution of the point kinetic equations at beginning of stage
*
*Parameters: ouput
* PARAMB  values of global parameters at end of stage
* T       time at end of stage
* Y       solution of the point kinetic equations at end of stage
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  Subroutine arguments
*----
      TYPE(C_PTR) IPMAP
      INTEGER NALPHA,NGROUP
      REAL LAMBDA,EPSILON,BETAI(NGROUP),LAMBDAI(NGROUP),PARAMI(NALPHA),
     1 PARAMB(NALPHA)
      DOUBLE PRECISION DT,T,Y(NGROUP+1)
*----
*  Local variables
*----
      PARAMETER(NRKMIN=100,NRKMAX=100000)
      DOUBLE PRECISION DH,DPP,P0,T1,BETA,MAXI,RHO(3)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: YSAV,YSUM,Y1,Y2,
     1 Y3,Y4
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A
*----
*  Scratch storage allocation
*----
      ALLOCATE(YSAV(NGROUP+1),YSUM(NGROUP+1),Y1(NGROUP+1),Y2(NGROUP+1),
     1 Y3(NGROUP+1),Y4(NGROUP+1),A(NGROUP+1,NGROUP+1))
*----
*   Runge-Kutta and calcul parameters initialisation
*----
      DPP=1.D0
      NRK=NRKMIN
      P0=-1.D0
*----
*  Set the Runge-Kutta evolution matrix
*----
      BETA=0.D0
      DO I=1,NGROUP
         BETA=BETA+BETAI(I)
      ENDDO
      CALL XDDSET(A,(NGROUP+1)**2,0.0D0)
      DO I=2,NGROUP+1
         A(I,1)=BETAI(I-1)/LAMBDA
      ENDDO
      DO I=2,NGROUP+1
         A(1,I)=LAMBDAI(I-1)
      ENDDO
      DO I=2,NGROUP+1
         A(I,I)=-LAMBDAI(I-1)
      ENDDO
*----
*  Runge-Kutta convergence loop
*----
      CALL XDDSET(RHO,3,0.0D0)
      DO WHILE ((DPP.GE.EPSILON).AND.(NRK.LE.NRKMAX))
*       time and time-step initialisation
        DH=DT/REAL(NRK)
        T1=T
*
*       save of the working vector
        DO I=1,NGROUP+1
          YSAV(I)=Y(I)
        ENDDO       
*
*       Runge-Kutta iteration loop
        DO I=1,NRK
*         total reactivity calculation with feedback
          IF(NALPHA.GT.0) CALL PKIRHO(IPMAP,NALPHA,T,DH,PARAMI,PARAMB,
     1    RHO)
*
*         Runge-Kutta procedure
          A(1,1)=(RHO(1)-BETA)/LAMBDA
          DO J=1,NGROUP+1
             Y1(J)=0.0D0
             DO K=1,NGROUP+1
                Y1(J)=Y1(J)+A(J,K)*Y(K)
             ENDDO
          ENDDO
          DO J=1,NGROUP+1
             YSUM(J)=Y(J)+(DH/6.D0)*Y1(J)
             Y1(J)=Y(J)+DH/2.D0*Y1(J)
          ENDDO
          A(1,1)=(RHO(2)-BETA)/LAMBDA
          DO J=1,NGROUP+1
             Y2(J)=0.0D0
             DO K=1,NGROUP+1
                Y2(J)=Y2(J)+A(J,K)*Y1(K)
             ENDDO
          ENDDO
          DO J=1,NGROUP+1
             YSUM(J)=YSUM(J)+(DH/3.D0)*Y2(J)
             Y2(J)=Y(J)+DH/2.D0*Y2(J)
          ENDDO
          A(1,1)=(RHO(2)-BETA)/LAMBDA
          DO J=1,NGROUP+1
             Y3(J)=0.0D0
             DO K=1,NGROUP+1
                Y3(J)=Y3(J)+A(J,K)*Y2(K)
             ENDDO
          ENDDO
          DO J=1,NGROUP+1
             YSUM(J)=YSUM(J)+(DH/3.D0)*Y3(J)
             Y3(J)=Y(J)+DH*Y3(J)
          ENDDO
          A(1,1)=(RHO(3)-BETA)/LAMBDA
          DO J=1,NGROUP+1
             Y4(J)=0.0D0
             DO K=1,NGROUP+1
                Y4(J)=Y4(J)+A(J,K)*Y3(K)
             ENDDO
          ENDDO
          DO J=1,NGROUP+1
             YSUM(J)=YSUM(J)+(DH/6.D0)*Y4(J)
             Y(J)=YSUM(J)
          ENDDO
          T=T+DH
*
*         convergence test initialisation
          MAXI=0.D0
          DO J=1,NGROUP+1
             MAXI=MAX(ABS(Y(J)),MAXI)
          ENDDO
          IF(MAXI.GT.1.0D30) GOTO 100
        ENDDO
* 
*       convergence test
 100    IF(P0.NE.-1.D0) DPP=ABS(Y(1)-P0)/ABS(P0)
        P0=Y(1)
*
*       reinitialisation of the number of Runge-Kutta time-steps
        NRK=2*NRK
*
*       reinitialisation of the working vector if not converged
        IF((DPP.GE.EPSILON).AND.(NRK.LE.NRKMAX)) THEN
          DO I=1,NGROUP+1
             Y(I)=YSAV(I)
          ENDDO
          T=T1
        ENDIF
      ENDDO
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(A,Y4,Y3,Y2,Y1,YSUM,YSAV)
      RETURN
      END
