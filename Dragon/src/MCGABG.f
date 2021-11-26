*DECK MCGABG
      SUBROUTINE MCGABG(IPRINT,LFORW,PACA,N,LC,EPSM,MAXM,IM,MCU,JU,
     1                  DIAGF,CF,ILUDF,ILUCF,RHS,F,FAC,LC0,IM0,MCU0)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve the ACA corrective system using BICGSTAB.
*
*Reference: (p382)
* MEURANT, G. 1999. "Computer Solution of Large Linear Systems".
* Studies in Mathematics and its Applications vol.28. North Holland.
* 776p.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* IPRINT  print parameter.
* LFORW   flag set to .false. to transpose the coefficient matrix.
* PACA    type of preconditioner to solve the ACA corrective system.
* N       dimension of the corrective system.
* LC      dimension of profiled matrices MCU and CQ.
* IM      connection matrix.
* MCU     connection matrix.
* DIAGF   diagonal elements of the matrix to inverse.
* CF      non-diagonal elements of the matrix to inverse.
* RHS     right hand-side of the corrective system (already
*         preconditioned).
* FAC     scaling factor for precision.
* LC0     used in ILU0-ACA acceleration.
* IM0     used in ILU0-ACA acceleration.
* MCU0    used in ILU0-ACA acceleration.
* EPSM    stopping criterion.
* MAXM    maximum number of iterations allowed.
*
*Parameters: output
* F       corrective fluxes and currents.
*
*Parameters: scratch
* JU      undefined.
* ILUDF   undefined.
* ILUCF   undefined.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPRINT,PACA,N,LC,IM(N+1),MCU(LC),JU(N),MAXM,LC0,IM0(*),
     1 MCU0(*)
      REAL DIAGF(N),CF(LC),ILUDF(N),ILUCF(LC),EPSM,FAC
      DOUBLE PRECISION RHS(N),F(N)
      LOGICAL LFORW
*----
* LOCAL VARIABLE
*----
      REAL EPSMAX,EPSINF,EPS2
      PARAMETER (EPSMAX=1E-7)
      INTEGER I,J,ITER
      DOUBLE PRECISION R,BI,WI,RT1
      DOUBLE PRECISION DDOT,AUX(2),EPS,FNORM,RHSN,ASIN,ASIN2,SIN,CN,SQ2
      LOGICAL DEBUG
      INTRINSIC SQRT,ABS,SIGN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PI,RI,SI,ROT,API,
     1 ASI
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PI(N),RI(N),SI(N),ROT(N),API(N),ASI(N))
*
      SQ2=1.D0/SQRT(2.D0)
*----
      DEBUG=.FALSE.
      EPSINF=EPSMAX*FAC
      ITER=0
*
      RHSN=0.0
      DO I=1,N
         RHSN=MAX(RHSN,ABS(RHS(I)))
      ENDDO
      IF(RHSN.LT.EPSINF) THEN
         DO I=1,N
            F(I)=0.0
         ENDDO
         IF(DEBUG) WRITE(6,200) RHSN,EPSINF
         GO TO 40
      ENDIF
      EPS2=EPSMAX*REAL(RHSN)
      EPS2=EPS2*EPS2
*---
*     initial corrective flux is set to rhs
*     calculate (P times (D times RHS)) -> RI
      CALL MCGPRA(LFORW,3,PACA,.FALSE.,N,LC,IM,MCU,JU,DIAGF,CF,ILUDF,
     1     ILUCF,DIAGF,RHS,RI,LC0,IM0,MCU0,CF)
      DO I=1,N
         F(I)=RHS(I)
         RI(I)=RHS(I)-RI(I)
         PI(I)=RI(I)
         ROT(I)=RI(I)
      ENDDO
      R=DDOT(N,RI,1,RI,1)
      FNORM=DDOT(N,F,1,F,1)
      EPS=SQRT(R/FNORM)
      IF(DEBUG) WRITE(6,100) ITER,EPS,EPSM
      IF(EPS.LE.EPSM) THEN
         IF(IPRINT.GT.2) WRITE(6,100) ITER,EPS,EPSM
         GO TO 40
      ENDIF
      AUX(1)=R
*     
      DO WHILE (ITER.LT.MAXM)
*     BiCGSTAB iterations
         ITER=ITER+1
*        calculate (P times (D times PI)) -> API
         CALL MCGPRA(LFORW,3,PACA,.FALSE.,N,LC,IM,MCU,JU,DIAGF,CF,
     1        ILUDF,ILUCF,DIAGF,PI,API,LC0,IM0,MCU0,CF)
*
         AUX(2)=AUX(1)/DDOT(N,API,1,ROT,1)
         DO J=1,N
            SI(J)=RI(J)-AUX(2)*API(J)
         ENDDO
         ITER=ITER+1
*        calculate (P times (D times SI)) -> ASI
         CALL MCGPRA(LFORW,3,PACA,.FALSE.,N,LC,IM,MCU,JU,DIAGF,CF,
     1        ILUDF,ILUCF,DIAGF,SI,ASI,LC0,IM0,MCU0,CF)
*
!!!!         ASIN=DDOT(N,ASI,1,ASI,1)
!!!!         ASIN2=DDOT(N,ASI,1,SI,1)
!!!!         IF(ASIN.GT.EPSMAX*ASIN2) THEN
!!!!            WI=ASIN2/ASIN
!!!!         ELSE
!!!!*        assuming lucky breakdown
!!!!            WI=1.0
!!!!         ENDIF
*        Modification proposed by Sleijpen and Van der Vorst (Numerical Algorithms, 10:203-223, 1995)
         ASIN2=DDOT(N,ASI,1,SI,1)
         ASIN=SQRT(DDOT(N,ASI,1,ASI,1))
         SIN=SQRT(DDOT(N,SI,1,SI,1))
         CN=ASIN*SIN
         IF(CN.GT.EPSMAX*ASIN2) THEN
            CN=ASIN2/CN
            WI=MAX(ABS(CN),SQ2)*SIN/ASIN
            WI=SIGN(WI,CN)
         ELSE
*        assuming lucky breakdown
            WI=1.0
         ENDIF           
*        calculate new iterate
         DO J=1,N
            F(J)=F(J)+AUX(2)*PI(J)+WI*SI(J)
            RI(J)=SI(J)-WI*ASI(J)
         ENDDO
         R=DDOT(N,RI,1,RI,1)
         FNORM=DDOT(N,F,1,F,1)
         IF(FNORM.LT.EPS2) GOTO 30
         EPS=SQRT(R/FNORM)
         IF(DEBUG) WRITE(6,100) ITER,EPS,EPSM
         IF(EPS.LE.EPSM) GO TO 20
         RT1=AUX(1)
         AUX(1)=DDOT(N,RI,1,ROT,1)
         BI=AUX(1)/RT1*AUX(2)/WI
         DO J=1,N
            PI(J)=RI(J)+BI*(PI(J)-WI*API(J))
         ENDDO
      ENDDO    
 20   CONTINUE
*     determine final residual norm
      ITER=ITER+1
*     calculate (P times (D times F)) -> RI
      CALL MCGPRA(LFORW,3,PACA,.FALSE.,N,LC,IM,MCU,JU,DIAGF,CF,ILUDF,
     1     ILUCF,DIAGF,F,RI,LC0,IM0,MCU0,CF)
      DO I=1,N
         RI(I)=RHS(I)-RI(I)
      ENDDO
*
      R=DDOT(N,RI,1,RI,1)
      FNORM=DDOT(N,F,1,F,1)
      IF(FNORM.LT.EPS2) GOTO 30
      EPS=SQRT(R/FNORM)
      IF(IPRINT.GT.2) WRITE(6,100) ITER,EPS,EPSM
!!!!      IF(EPS.GT.EPSM) THEN
!!!!         DO I=1,N
!!!!            PI(I)=RI(I)
!!!!            ROT(I)=RI(I)
!!!!         ENDDO 
!!!!         GO TO 10
!!!!      ENDIF
      GO TO 40
*
 30   IF(DEBUG) WRITE(6,300) ITER,FNORM,EPS2
      CALL XDRSET(F,N,0.0)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
 40   DEALLOCATE(ASI,API,ROT,SI,RI,PI)
      RETURN
*
 100  FORMAT(12X,14H MCGABG: ITER=,I3,5H EPS=,E9.2,5H TAR=,E9.2)
 200  FORMAT(12X,27H MCGABG: RHS INFINITE NORM=,E9.2,5H LIM=,E9.2/
     1       12X,33H -> ACA CORRECTION IS SET TO ZERO)
 300  FORMAT(12X,14H MCGABG: ITER=,I3,7H FNORM=,E9.2,5H LIM=,E9.2)
      END
