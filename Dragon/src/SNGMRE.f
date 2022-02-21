*DECK SNGMRE
      SUBROUTINE SNGMRE (KPSYS,NGIND,IPTRK,IMPX,NGEFF,NREG,NBMIX,
     1 NUN,NSTART,MAXIT,EPSINR,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,NBS,KPSOU1,
     2 KPSOU2,FLUXC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve N-group transport equation for fluxes using the discrete
* ordinates (SN) method with a GMRES(m) acceleration.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* KPSYS   pointer to the assembly matrices. KPSYS is an array of
*         directories.
* NGIND   energy group indices assign to the NGEFF set.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IMPX    print flag (equal to zero for no print).
* NGEFF   number of energy groups processed in parallel.
* NREG    total number of regions for which specific values of the
*         neutron flux and reactions rates are required.
* NBMIX   number of mixtures.
* NUN     total number of unknowns in vectors SUNKNO and FUNKNO.
* NSTART  restarts the GMRES method every NSTART iterations.
* MAXIT   maximum number of inner iterations.
* EPSINR  convergence criterion.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  position of averaged flux elements in FUNKNO vector.
* SUNKNO  input source vector.
*
*Parameters: input/output
* FUNKNO  unknown vector.
* FLUXC   flux at the cutoff energy.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER     NGEFF,NGIND(NGEFF),IMPX,NREG,NBMIX,NUN,NSTART,MAXIT,
     1            MAT(NREG),KEYFLX(NREG),NBS(NGEFF)
      TYPE(C_PTR) KPSYS(NGEFF),IPTRK,KPSOU1(NGEFF),KPSOU2(NGEFF)
      REAL        EPSINR,VOL(NREG),FUNKNO(NUN,NGEFF),SUNKNO(NUN,NGEFF)
      REAL,OPTIONAL :: FLUXC(NREG)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6)
      REAL        SDOT
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KMAX
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RR,QQ
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: EPS1,RHO
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: G,C,S,X
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: V,H
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: INCONV
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(G(NSTART+1,NGEFF),C(NSTART+1,NGEFF),S(NSTART+1,NGEFF),
     1 X(NUN,NGEFF),V(NUN,NSTART+1,NGEFF),H(NSTART+1,NSTART,NGEFF),
     1 EPS1(NGEFF),INCONV(NGEFF),KMAX(NGEFF),RHO(NGEFF))
      ALLOCATE(RR(NUN,NGEFF),QQ(NUN,NGEFF))
*----
*  GLOBAL GMRES ITERATION.
*----
      DO II=1,NGEFF
        EPS1(II)=EPSINR*SQRT(SDOT(NUN,SUNKNO(1,II),1,SUNKNO(1,II),1))
        RHO(II)=1.0E10
      ENDDO
      LNCONV=NGEFF
      CALL XDLSET(INCONV,NGEFF,.TRUE.)
      ITER=0
      NITER=1
      DO WHILE((LNCONV.GT.0).AND.(ITER.LT.MAXIT))
        RR(:NUN,:NGEFF)=FUNKNO(:NUN,:NGEFF)
        CALL SNFLUX(KPSYS,INCONV,NGIND,IPTRK,IMPX,NGEFF,NREG,
     1  NBMIX,NUN,MAT,VOL,KEYFLX,RR,SUNKNO,NITER,NBS,KPSOU1,KPSOU2,
     2  FLUXC)
        NITER=NITER+1
        DO II=1,NGEFF
          IF(.NOT.INCONV(II)) CYCLE
          X(:NUN,II)=RR(:NUN,II)-FUNKNO(:NUN,II)
          RHO(II)=SQRT(DOT_PRODUCT(X(:NUN,II),X(:NUN,II)))
          IF(RHO(II).LT.EPS1(II)) THEN
            LNCONV=LNCONV-1
            INCONV(II)=.FALSE.
          ENDIF
        ENDDO
*----
*  TEST FOR TERMINATION ON ENTRY
*----
        IF(LNCONV.EQ.0) GO TO 100
*
        CALL XDDSET(G,(NSTART+1)*NGEFF,0.0D0)
        CALL XDDSET(H,(NSTART+1)*NSTART*NGEFF,0.0D0)
        CALL XDDSET(V,NUN*(NSTART+1)*NGEFF,0.0D0)
        CALL XDDSET(C,(NSTART+1)*NGEFF,0.0D0)
        CALL XDDSET(S,(NSTART+1)*NGEFF,0.0D0)
        CALL XDISET(KMAX,NGEFF,0)
        DO II=1,NGEFF
          IF(.NOT.INCONV(II)) CYCLE
          G(1,II)=RHO(II)
          DO I=1,NUN
            V(I,1,II)=X(I,II)/RHO(II)
            X(I,II)=FUNKNO(I,II)
          ENDDO
        ENDDO
*----
*  GMRES(1) ITERATION
*----
        K=0
        DO WHILE((LNCONV.GT.0).AND.(K.LT.NSTART).AND.(ITER.LT.MAXIT))
          K=K+1
          ITER=ITER+1
          IF(IMPX.GT.2) WRITE(IUNOUT,300) ITER,MAXVAL(RHO(:NGEFF)),
     1    LNCONV
          RR(:NUN,:NGEFF)=0.0
          QQ(:NUN,:NGEFF)=0.0
          DO II=1,NGEFF
            IF(.NOT.INCONV(II)) CYCLE
            RR(:NUN,II)=REAL(V(:NUN,K,II))
          ENDDO
          CALL SNFLUX(KPSYS,INCONV,NGIND,IPTRK,IMPX,NGEFF,NREG,
     1    NBMIX,NUN,MAT,VOL,KEYFLX,RR,QQ,NITER,NBS,KPSOU1,KPSOU2,FLUXC)
          NITER=NITER+1
          DO II=1,NGEFF
            IF(.NOT.INCONV(II)) CYCLE
            V(:NUN,K+1,II)=V(:NUN,K,II)-RR(:NUN,II)
            KMAX(II)=K
*----
*  MODIFIED GRAM-SCHMIDT
*----
            DO J=1,K
              HR=DOT_PRODUCT(V(:NUN,J,II),V(:NUN,K+1,II))
              H(J,K,II)=HR
              V(:NUN,K+1,II)=V(:NUN,K+1,II)-HR*V(:NUN,J,II)
            ENDDO
            H(K+1,K,II)=SQRT(DOT_PRODUCT(V(:NUN,K+1,II),V(:NUN,K+1,II)))
*----
*  REORTHOGONALIZE
*----
            DO J=1,K
              HR=DOT_PRODUCT(V(:NUN,J,II),V(:NUN,K+1,II))
              H(J,K,II)=H(J,K,II)+HR
              V(:NUN,K+1,II)=V(:NUN,K+1,II)-HR*V(:NUN,J,II)
            ENDDO
            H(K+1,K,II)=SQRT(DOT_PRODUCT(V(:NUN,K+1,II),V(:NUN,K+1,II)))
*----
*  WATCH OUT FOR HAPPY BREAKDOWN 
*----
            IF(H(K+1,K,II).NE.0.0) THEN
              V(:NUN,K+1,II)=V(:NUN,K+1,II)/H(K+1,K,II)
            ENDIF
*----
*  FORM AND STORE THE INFORMATION FOR THE NEW GIVENS ROTATION
*----
            DO I=1,K-1
              W1=C(I,II)*H(I,K,II)-S(I,II)*H(I+1,K,II)
              W2=S(I,II)*H(I,K,II)+C(I,II)*H(I+1,K,II)
              H(I,K,II)=W1
              H(I+1,K,II)=W2
            ENDDO
            ZNU=SQRT(H(K,K,II)**2+H(K+1,K,II)**2)
            IF(ZNU.NE.0.0) THEN
              C(K,II)=H(K,K,II)/ZNU
              S(K,II)=-H(K+1,K,II)/ZNU
              H(K,K,II)=C(K,II)*H(K,K,II)-S(K,II)*H(K+1,K,II)
              H(K+1,K,II)=0.0D0
              W1=C(K,II)*G(K,II)-S(K,II)*G(K+1,II)
              W2=S(K,II)*G(K,II)+C(K,II)*G(K+1,II)
              G(K,II)=W1
              G(K+1,II)=W2
            ENDIF
*----
*  UPDATE THE RESIDUAL NORM
*----
            RHO(II)=ABS(G(K+1,II))
            IF(RHO(II).LE.EPS1(II)) THEN
              INCONV(II)=.FALSE.
              LNCONV=LNCONV-1
            ENDIF
          ENDDO
        ENDDO
*----
*  AT THIS POINT EITHER K > NSTART OR RHO < EPS1.
*  IT'S TIME TO COMPUTE X AND CYCLE.
*----
        DO II=1,NGEFF
          K=KMAX(II)
          IF(K.EQ.0) CYCLE
          G(K,II)=G(K,II)/H(K,K,II)
          DO L=K-1,1,-1
            W1=G(L,II)-DOT_PRODUCT(H(L,L+1:K,II),G(L+1:K,II))
            G(L,II)=W1/H(L,L,II)
          ENDDO
          DO J=1,K
            X(:,II)=X(:,II)+G(J,II)*V(:,J,II)
          ENDDO
          FUNKNO(:,II)=REAL(X(:,II))
        ENDDO
      ENDDO
*
      IF(IMPX.GT.1) WRITE(IUNOUT,'(32H SNGMRE: NUMBER OF ONE-SPEED ITE,
     1 8HRATIONS=,I5,1H.)') ITER
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  100 DEALLOCATE(QQ,RR)
      DEALLOCATE(RHO,KMAX,INCONV,EPS1,H,V,X,S,C,G)
      RETURN
*
  300 FORMAT(28H SNGMRE: ONE-SPEED ITERATION,I4,10H  L2 NORM=,1P,E11.4,
     1 23H  NON-CONVERGED GROUPS=,I5)
      END
