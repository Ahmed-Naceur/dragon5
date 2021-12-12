*DECK MCGABGR
      SUBROUTINE MCGABGR(IPRINT,LFORW,PACA,N,NG,NFIRST,NGEFF,M,LC,NGIND,
     1                   NGINDV,NCONV,KPSYS,JPMACR,NZON,IPERM,IM,MCU,JU,
     2                   EPSM,MAXM,RHS,F,FAC,LC0,IM0,MCU0)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve the ACA corrective system (in a rebalancing form) using 
* BICGSTAB.
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
* N       number of unknowns per group.
* NG      total number of groups.
* NFIRST  first group to proceed.
* NGEFF   number of unconverged groups.
* M       number of material mixtures.
* LC      dimension of profiled matrices MCU and CQ.
* NGIND   index of the groups to process.
* NGINDV  index to pass from "NGEFF format" to "NG format".
* NCONV   logical array of convergence status for each group (.TRUE.:
*         not converged).
* KPSYS   pointer array for each group properties.
* JPMACR  pointer to the macrolib LCM object ('GROUP' directory).
* NZON    index-number of the mixture type assigned to each volume.
* IPERM   permutation array for ACA.
* IM      connection matrix.
* MCU     connection matrix.
* JU      used for ilu0 preconditioner.
* EPSM    stopping criterion.
* MAXM    maximum number of iterations allowed.
* RHS     right hand-side of the corrective system (already
*         preconditioned).
* FAC     scaling factor for precision.
* LC0     used in ILU0-ACA acceleration.
* IM0     used in ILU0-ACA acceleration.
* MCU0    used in ILU0-ACA acceleration.
*
*Parameters: output
* F       corrective fluxes and currents.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPRINT,PACA,N,NG,NFIRST,NGEFF,M,LC,NGIND(NGEFF),
     1 NGINDV(NG),KPSYS(NGEFF),JPMACR,NZON(N),IPERM(N),IM(N+1),MCU(LC),
     2 JU(N),MAXM,LC0,IM0(*),MCU0(*)
      REAL EPSM,FAC
      DOUBLE PRECISION RHS(N,NGEFF),F(N,NGEFF)
      LOGICAL LFORW,NCONV(NGEFF)
*----
* LOCAL VARIABLE
*----
      REAL EPSMAX,EPSINF,EPS2
      PARAMETER (EPSMAX=1E-7)
      INTEGER I,II,J,ITER
      DOUBLE PRECISION R,BI,WI,RT1,ASIN,ASIN2,SQ2
      DOUBLE PRECISION DDOT,AUX(2),EPS,FNORM,RHSN
      LOGICAL DEBUG
      INTRINSIC SQRT,ABS
      REAL, ALLOCATABLE, DIMENSION(:) :: TEMP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: PI,RI,SI,ROT,API,
     1 ASI
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PI(N,NGEFF),RI(N,NGEFF),SI(N,NGEFF),ROT(N,NGEFF),
     1 API(N,NGEFF),ASI(N,NGEFF))
*
      SQ2=1.D0/SQRT(2.D0)
*---
      DEBUG=.FALSE.
      EPSINF=EPSMAX*FAC
      ITER=0
      ALLOCATE(TEMP(2*N))
*
      RHSN=0.0
      DO II=NFIRST,NGEFF
         IF (NCONV(II)) THEN
            DO I=1,N
               RHSN=MAX(RHSN,ABS(RHS(I,II)))
            ENDDO
         ENDIF
      ENDDO
      IF (RHSN.LT.EPSINF) THEN
         DO II=NFIRST,NGEFF
            IF (NCONV(II)) THEN
               DO I=1,N
                  F(I,II)=0.0
               ENDDO
            ENDIF
         ENDDO
         IF (DEBUG) WRITE(6,200) RHSN,EPSINF
         DEALLOCATE(TEMP)
         GO TO 40
      ENDIF
      EPS2=EPSMAX*REAL(RHSN)
      EPS2=EPS2*EPS2
*---
*     initial corrective flux is set to rhs
*     calculate (P times (D times RHS)) -> RI
      CALL MCGACA(LFORW,PACA,N,NG,NFIRST,NGEFF,M,LC,NGIND,NGINDV,NCONV,
     1     KPSYS,JPMACR,NZON,IPERM,IM,MCU,JU,RHS,RI,TEMP,LC0,IM0,MCU0)
      R=0.0
      FNORM=0.0
      DO II=NFIRST,NGEFF
         IF (NCONV(II)) THEN
            DO I=1,N
               F(I,II)=REAL(RHS(I,II))
               RI(I,II)=REAL(RHS(I,II))-RI(I,II)
               PI(I,II)=RI(I,II)
               ROT(I,II)=RI(I,II)
            ENDDO
            R=R+DDOT(N,RI(1,II),1,RI(1,II),1)
            FNORM=FNORM+DDOT(N,F(1,II),1,F(1,II),1)
         ENDIF
      ENDDO
      EPS=SQRT(R/FNORM)
      IF (DEBUG) WRITE(6,100) ITER,EPS,EPSM
      IF (EPS.LE.EPSM) THEN
         DEALLOCATE(TEMP)
         GO TO 40
      ENDIF
      AUX(1)=R !!DDOT(N,RI,1,ROT,1)
*     
      DO WHILE (ITER.LT.MAXM)
*     BiCGSTAB iterations
         ITER=ITER+1
*        calculate (P times (D times PI)) -> API
         CALL MCGACA(LFORW,PACA,N,NG,NFIRST,NGEFF,M,LC,NGIND,NGINDV,
     1        NCONV,KPSYS,JPMACR,NZON,IPERM,IM,MCU,JU,PI,API,TEMP,LC0,
     2        IM0,MCU0)
*
         AUX(2)=0.0
         DO II=NFIRST,NGEFF
            IF (NCONV(II)) THEN
               AUX(2)=AUX(2)+DDOT(N,API(1,II),1,ROT(1,II),1)
            ENDIF
         ENDDO
         AUX(2)=AUX(1)/AUX(2) 
         DO II=NFIRST,NGEFF
            IF (NCONV(II)) THEN
               DO J=1,N
                  SI(J,II)=RI(J,II)-AUX(2)*API(J,II)
               ENDDO
            ENDIF
         ENDDO
         ITER=ITER+1
*        calculate (P times (D times SI)) -> ASI
         CALL MCGACA(LFORW,PACA,N,NG,NFIRST,NGEFF,M,LC,NGIND,NGINDV,
     1        NCONV,KPSYS,JPMACR,NZON,IPERM,IM,MCU,JU,SI,ASI,TEMP,LC0,
     2        IM0,MCU0)
*
         ASIN2=0.0
         ASIN=0.0
         DO II=NFIRST,NGEFF
            IF (NCONV(II)) THEN
               ASIN2=ASIN2+DDOT(N,ASI(1,II),1,SI(1,II),1)
               ASIN=ASIN+DDOT(N,ASI(1,II),1,ASI(1,II),1)
            ENDIF
         ENDDO
         IF (ASIN.GT.EPSMAX*ASIN2) THEN
            WI=ASIN2/ASIN
         ELSE
*        assuming lucky breakdown
            WI=1.0
         ENDIF
!!!!*        Modification proposed by Sleijpen and Van der Vorst (Numerical Algorithms, 10:203-223, 1995)
!!!!         ASIN2=0.0
!!!!         ASIN=0.0
!!!!         SIN=0.0
!!!!         DO II=NFIRST,NGEFF
!!!!            IF (NCONV(II)) THEN
!!!!               ASIN2=ASIN2+DDOT(N,ASI(1,II),1,SI(1,II),1)
!!!!               ASIN=ASIN+DDOT(N,ASI(1,II),1,ASI(1,II),1)
!!!!               SIN=SIN+DDOT(N,SI(1,II),1,SI(1,II),1)
!!!!            ENDIF
!!!!         ENDDO
!!!!         ASIN=SQRT(ASIN)
!!!!         SIN=SQRT(SIN)
!!!!         CN=ASIN*SIN
!!!!         IF (CN.GT.EPSMAX*ASIN2) THEN
!!!!            CN=ASIN2/CN
!!!!            WI=MAX(ABS(CN),SQ2)*SIN/ASIN
!!!!            WI=SIGN(WI,CN)
!!!!         ELSE
!!!!*        assuming lucky breakdown
!!!!            WI=1.0
!!!!         ENDIF
*        calculate new iterate
         R=0.0
         FNORM=0.0
         DO II=NFIRST,NGEFF
            IF (NCONV(II)) THEN
               DO J=1,N
                  F(J,II)=F(J,II)+AUX(2)*PI(J,II)+WI*SI(J,II)
                  RI(J,II)=SI(J,II)-WI*ASI(J,II)
               ENDDO
               R=R+DDOT(N,RI(1,II),1,RI(1,II),1)
               FNORM=FNORM+DDOT(N,F(1,II),1,F(1,II),1)
            ENDIF
         ENDDO
         IF (FNORM.LT.EPS2) GOTO 30
         EPS=SQRT(R/FNORM)
         IF (DEBUG) WRITE(6,100) ITER,EPS,EPSM
         IF (EPS.LE.EPSM) GO TO 20
         RT1=AUX(1)
         AUX(1)=0.0
         DO II=NFIRST,NGEFF
            IF (NCONV(II)) THEN         
               AUX(1)=AUX(1)+DDOT(N,RI(1,II),1,ROT(1,II),1)
            ENDIF
         ENDDO
         BI=AUX(1)/RT1*AUX(2)/WI
         DO II=NFIRST,NGEFF
            IF (NCONV(II)) THEN    
               DO J=1,N  
                  PI(J,II)=RI(J,II)+BI*(PI(J,II)-WI*API(J,II))
               ENDDO
            ENDIF
         ENDDO
      ENDDO    
 20   CONTINUE
*     determine final residual norm
      ITER=ITER+1
*     calculate (P times (D times F)) -> RI
      CALL MCGACA(LFORW,PACA,N,NG,NFIRST,NGEFF,M,LC,NGIND,NGINDV,NCONV,
     1     KPSYS,JPMACR,NZON,IPERM,IM,MCU,JU,F,RI,TEMP,LC0,IM0,MCU0)
*
      R=0.0
      FNORM=0.0 
      DO II=NFIRST,NGEFF
         IF (NCONV(II)) THEN    
            DO I=1,N
               RI(I,II)=RHS(I,II)-RI(I,II)
            ENDDO
            R=R+DDOT(N,RI(1,II),1,RI(1,II),1)
            FNORM=FNORM+DDOT(N,F(1,II),1,F(1,II),1)
         ENDIF
      ENDDO
      IF (FNORM.LT.EPS2) GOTO 30
      EPS=SQRT(R/FNORM)
      IF (IPRINT.GT.0) WRITE(6,400) EPS,ITER
!!!!      R=0.0
!!!!      FNORM=0.0
!!!!      EPS=0.0
!!!!      DO II=NFIRST,NGEFF
!!!!         IF (NCONV(II)) THEN
!!!!            DO I=1,N
!!!!               R=MAX(R,ABS(RI(I,II)))
!!!!               FNORM=MAX(FNORM,ABS(F(I,II)))
!!!!            ENDDO
!!!!            EPS=MAX(EPS,R/FNORM)
!!!!         ENDIF
!!!!      ENDDO
!!!!      WRITE(*,*) '              PRC=',EPS
      DEALLOCATE(TEMP)
      GO TO 40
*
 30   IF (DEBUG) WRITE(6,300) ITER,FNORM,EPS2 
      DO II=NFIRST,NGEFF
         IF (NCONV(II)) THEN    
            DO I=1,N
               F(I,II)=0.0
            ENDDO
         ENDIF
      ENDDO
      DEALLOCATE(TEMP)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  40  DEALLOCATE(ASI,API,ROT,SI,RI,PI)
      RETURN
*
 100  FORMAT(9X,14H MCGABGR:ITER=,I3,5H EPS=,E9.2,5H TAR=,E9.2)
 200  FORMAT(9X,27H MCGABGR:RHS INFINITE NORM=,E9.2,5H LIM=,E9.2/
     1       9X,33H -> ACA CORRECTION IS SET TO ZERO)
 300  FORMAT(9X,14H MCGABGR:ITER=,I3,7H FNORM=,E9.2,5H LIM=,E9.2)
 400  FORMAT(10X,48HACA: UP-SCATTE. GROUPS: MULTIGROUP BICGSTAB: PRC:,
     1       E9.2,2H (,I4,12H ITERATIONS))
      END
