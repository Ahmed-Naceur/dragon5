*DECK MCGMRE
      SUBROUTINE MCGMRE(SUBFFI,SUBFFA,SUBLDC,SUBSCH,CYCLIC,KPSYS,
     1           IPRINT,IPTRK,IFTRAK,IPMACR,NDIM,KV,NUN,NLONG,PHIOUT,
     2           NZON,MATALB,M,NANI,NMU,NMAX,NANGL,NREG,NSOUT,SOUR,IAAC,
     3           ISCR,LC,LFORW,PACA,ITST,MAXI,QFR,PHIIN,CAZ0,CAZ1,CAZ2,
     4           CPO,ZMU,WZMU,VOL,EPS,ERRTOL,REPSI,NSTART,SIGAL,LPS,
     5           NGROUP,NGEFF,NGIND,NCONV,LNCONV,NLIN,NFUNL,KEYFLX,
     6           KEYCUR,STIS,NPJJM,REBFLG,LPRISM,N2REG,N2SOU,NZP,DELU,
     7           FACSYM,IDIR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve the linear system obtained by the characteristics formalism
* with GMRES iterative approach.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier and A. Hebert
*
*Parameters:
* SUBFFI  flux integration subroutine with isotropic source.
* SUBFFA  flux integration subroutine with anisotropic source.
* SUBLDC  flux integration subroutine with linear-discontinuous source.
* SUBSCH  track coefficients calculation subroutine.
* CYCLIC  cyclic tracking flag.
* KPSYS   pointer array for each group properties.
* IPRINT  print parameter (equal to zero for no print).
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  tracking file unit number.
* IPMACR  pointer to the macrolib LCM object.
* NDIM    number of dimensions for the geometry.
* KV      total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NUN     total number of unknowns in vectors QFR and PHIIN.
* NLONG   number of spatial unknowns.
* PHIOUT  output flux vector.
* NZON    mixture-albedo index array in MCCG format.
* MATALB  albedo-mixture index array in MOCC format.
* M       number of material mixtures.
* NANI    scattering anisotropy (=1 for isotropic scattering).
* NMU     order of the polar quadrature set.
* NMAX    maximum number of elements in a track.
* NANGL   number of tracking angles in the plane.
* NREG    number of regions (volumes).
* NSOUT   number of outer surfaces.
* SOUR    undefined.
* IAAC    no acceleration / CDD acceleration of inner iterations (0/1).
* ISCR    no acceleration / SCR acceleration of inner iterations (0/1).
* LC      dimension of profiled matrices MCU and CQ.
* LFORW   flag set to .false. to transpose the coefficient matrix.
* PACA    type of preconditioner to solve the ACA corrective system.
* ITST    output: number of inner iterations.
* MAXI    maximum number of inner iterations allowed.
* QFR     input source vector.
* PHIIN   input flux vector.
* CAZ0    cosines of the tracking polar angles in 3D.
* CAZ1    first cosines of the different tracking azimuthal angles.
* CAZ2    second cosines of the different tracking azimuthal angles.
* CPO     cosines of the different tracking polar angles in 2D.
* ZMU     polar quadrature set.
* WZMU    polar quadrature set.
* VOL     volumes.
* EPS     precision reached after min(MAXI,ITST) iterations.
* ERRTOL  tolerance for stopping criterion. process is stopped
*         as soon as ||Phi(n+1)-Phi(n)||/||Phi(n)|| <= EPS
*         with ||.|| the euclidean norm.
* REPSI   array containing precision of each iteration.
* NSTART  undefined.
* SIGAL   total cross-section and albedo array.
* LPS     used in scr acceleration.
* NGROUP  number of groups.
* NGEFF   number of groups to process.
* NGIND   index of the groups to process.
* NCONV   array of convergence flag for each group.
* LNCONV  number of unconverged groups.
* NLIN    number of polynomial components in flux spatial expansion.
* NFUNL   number of moments of the flux (in 2D: NFUNL=NANI*(NANI+1)/2).
* KEYFLX  position of flux elements in PHIOUT vector.
* KEYCUR  position of current elements in PHIOUT vector.
* STIS    Source term isolation option for flux integration.
* NPJJM   number of pjj modes to store for STIS option.
* REBFLG  ACA or SCR rebalancing flag.
* LPRISM  3D prismatic extended tracking flag.
* N2REG   number of regions in the 2D tracking if LPRISM.
* N2SOU   number of external surfaces in the 2D tracking if LPRISM.
* NZP     number of z-plans if LPRISM.
* DELU    input track spacing for 3D track reconstruction if LPRISM.
* FACSYM  tracking symmetry factor for maximum track length if LPRISM.
* IDIR    direction of fundamental current for TIBERE with MoC 
*         (=0,1,2,3). 
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*----
*  Subroutine arguments
*----
      TYPE(C_PTR) KPSYS(NGEFF),IPTRK,IPMACR
      INTEGER NGEFF,IPRINT,IFTRAK,NDIM,KV,NUN,NLONG,NGROUP,NZON(NLONG),
     1 M,NANI,NMU,NMAX,IAAC,LC,PACA,ITST(NGEFF),NSTART,MAXI,NANGL,NREG,
     2 NSOUT,ISCR,LPS,NGIND(NGEFF),LNCONV,NLIN,NFUNL,KEYFLX(NREG,NLIN,
     3 NFUNL),KEYCUR(NLONG-NREG),STIS,NPJJM,MATALB(-NSOUT:NREG),N2REG,
     4 N2SOU,NZP,IDIR
      REAL QFR(NUN,NGROUP),PHIIN(NUN,NGEFF),CPO(NMU),ZMU(NMU),
     1 WZMU(NMU),VOL(NLONG),EPS(NGEFF),ERRTOL,REPSI(MAXI,NGEFF),
     2 SIGAL(-6:M,NGEFF),DELU,FACSYM
      DOUBLE PRECISION PHIOUT(NUN,NGEFF),CAZ0(NANGL),CAZ1(NANGL),
     1 CAZ2(NANGL),SOUR(NUN,NGEFF)
      LOGICAL LFORW,CYCLIC,NCONV(NGEFF),REBFLG,LPRISM
      EXTERNAL SUBFFI,SUBFFA,SUBLDC,SUBSCH
*---
*   Local variables
*---
      REAL EPSINTO,FAC
      LOGICAL RHSFLG
      PARAMETER (FAC=100.0,MAXINT=200)
*---
*   Allocatable arrays
*---
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KMAX
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RHS, GAR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DENOM, RHO
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: R, G, C, S, FLOUT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: V, H
*---
*   Scratch storage allocation
*---
      ALLOCATE(DENOM(NGEFF),RHO(NGEFF),KMAX(NGEFF),FLOUT(NUN,NGEFF))
      ALLOCATE(RHS(NUN,NGEFF),GAR(NUN,NGEFF),V(NUN,NGEFF,NSTART+1),
     1 G(NSTART+1,NGEFF),H(NSTART+1,NSTART,NGEFF),C(NSTART+1,NGEFF),
     2 S(NSTART+1,NGEFF),R(NUN,NGEFF))
*
      IF(MAXI.LT.3) CALL XABORT('MCGMRE: MAXI MUST BE >= 3.')
      RHSFLG=.TRUE.
      EPSINTO=ERRTOL/FAC
      MAXIT=MAXI-1
      NCONV(:)=.FALSE.
      RHO(:)=0.0D0
      LNCONV=0
      DO II=1,NGEFF
        IG=NGIND(II)
        DENOM(II)=SQRT(DOT_PRODUCT(QFR(:,IG),QFR(:,IG)))
        RHO(II)=1.0D20
        NCONV(II)=(DENOM(II) /= 0.0D0)
        IF(NCONV(II)) LNCONV=LNCONV+1
        ITST(II)=0
        EPS(II)=0.0
      ENDDO
*---
*   Global GMRES(M) iteration
*---
      ITER=0
      DO WHILE((LNCONV /= 0) .AND. (ITER < MAXIT))
        ITER=ITER+1
        CALL MCGFL1(SUBFFI,SUBFFA,SUBLDC,SUBSCH,CYCLIC,KPSYS,IPRINT,
     1       IPTRK,IFTRAK,IPMACR,NDIM,KV,NUN,NLONG,PHIOUT,NZON,MATALB,
     2       M,NANI,NMU,NMAX,NANGL,NREG,NSOUT,NGROUP,NGEFF,NGIND,SOUR,
     3       IAAC,ISCR,LC,LFORW,PACA,EPSINTO,MAXINT,NLIN,NFUNL,KEYFLX,
     4       KEYCUR,QFR,PHIIN,CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,VOL,SIGAL,
     5       LPS,NCONV,.TRUE.,.FALSE.,STIS,NPJJM,REBFLG,LPRISM,N2REG,
     6       N2SOU,NZP,DELU,FACSYM,IDIR)
        ERROR=0.0D0
        DO II=1,NGEFF
          REPSI(ITER,II)=0.0
          IF(.NOT.NCONV(II)) CYCLE
          IG=NGIND(II)
          R(:,II)=PHIOUT(:,II)-PHIIN(:,II)
          RHO(II)=SQRT(DOT_PRODUCT(R(:,II),R(:,II)))
          REPSI(ITER,II)=REAL(RHO(II)/DENOM(II))
          IF(IPRINT.GT.4) WRITE(6,200) ITER,II,REPSI(ITER,II)
          EPS(II)=REPSI(ITER,II)
          ITST(II)=ITER
          ERROR=MAX(ERROR,RHO(II)/DENOM(II))
          IF(RHO(II) < ERRTOL*DENOM(II)) THEN
            NCONV(II)=.FALSE.
            LNCONV=LNCONV-1
          ENDIF
        ENDDO
*
*       Test do termination on entry
        IF(LNCONV == 0) EXIT
*
        H(:,:,:)=0.0D0
        V(:,:,:)=0.0D0
        C(:,:)=0.0D0
        S(:,:)=0.0D0
        G(:,:)=0.0D0
        KMAX(:)=0
        DO II=1,NGEFF
          IF(.NOT.NCONV(II)) CYCLE
          G(1,II)=RHO(II)
          V(:,II,1)=R(:,II)/RHO(II)
        ENDDO
*---
*   Evaluate RHS of the linear system
*---
        IF(RHSFLG) THEN
          IF(IPRINT > 3) WRITE(6,100) IAAC,ISCR
          RHS(:,:)=0.0
          CALL MCGFL1(SUBFFI,SUBFFA,SUBLDC,SUBSCH,CYCLIC,KPSYS,IPRINT,
     1         IPTRK,IFTRAK,IPMACR,NDIM,KV,NUN,NLONG,FLOUT,NZON,MATALB,
     2         M,NANI,NMU,NMAX,NANGL,NREG,NSOUT,NGROUP,NGEFF,NGIND,SOUR,
     3         IAAC,ISCR,LC,LFORW,PACA,EPSINTO,MAXINT,NLIN,NFUNL,KEYFLX,
     4         KEYCUR,QFR,RHS,CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,VOL,SIGAL,
     5         LPS,NCONV,.TRUE.,.FALSE.,STIS,NPJJM,REBFLG,LPRISM,N2REG,
     6         N2SOU,NZP,DELU,FACSYM,IDIR)
          DO II=1,NGEFF
            IF(NCONV(II)) RHS(:,II)=REAL(FLOUT(:,II))
          ENDDO
          RHSFLG=.FALSE.
        ENDIF
*---
*   GMRES(1) iteration
*---
        K=0
        DO WHILE((LNCONV /= 0) .AND. (K < NSTART) .AND. (ITER < MAXIT))
          K=K+1
          ITER=ITER+1
          GAR(:,:)=REAL(V(:,:,K))
          CALL MCGFL1(SUBFFI,SUBFFA,SUBLDC,SUBSCH,CYCLIC,KPSYS,IPRINT,
     1         IPTRK,IFTRAK,IPMACR,NDIM,KV,NUN,NLONG,FLOUT,NZON,MATALB,
     2         M,NANI,NMU,NMAX,NANGL,NREG,NSOUT,NGROUP,NGEFF,NGIND,SOUR,
     3         IAAC,ISCR,LC,LFORW,PACA,EPSINTO,MAXINT,NLIN,NFUNL,KEYFLX,
     4         KEYCUR,QFR,GAR,CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,VOL,SIGAL,
     5         LPS,NCONV,.TRUE.,.FALSE.,STIS,NPJJM,REBFLG,LPRISM,N2REG,
     6         N2SOU,NZP,DELU,FACSYM,IDIR)
          V(:,:,K+1)=V(:,:,K)-FLOUT(:,:)+RHS(:,:)
          ERROR=0.0D0
          DO II=1,NGEFF
            REPSI(ITER,II)=0.0
            IF(.NOT.NCONV(II)) CYCLE
            IG=NGIND(II)
            KMAX(II)=K
*
*           Modified Gram-Schmidt
            DO J=1,K
              H(J,K,II)=DOT_PRODUCT(V(:,II,J),V(:,II,K+1))
              V(:,II,K+1)=V(:,II,K+1)-H(J,K,II)*V(:,II,J)
            ENDDO
            H(K+1,K,II)=SQRT(DOT_PRODUCT(V(:,II,K+1),V(:,II,K+1)))
*
*           Reorthogonalize
            DO J=1,K
              HR=DOT_PRODUCT(V(:,II,J),V(:,II,K+1))
              H(J,K,II)=H(J,K,II)+HR
              V(:,II,K+1)=V(:,II,K+1)-HR*V(:,II,J)
            ENDDO
            H(K+1,K,II)=SQRT(DOT_PRODUCT(V(:,II,K+1),V(:,II,K+1)))
*
            ! Watch out do happy breakdown 
            IF(H(K+1,K,II) /= 0.0D0) V(:,II,K+1)=V(:,II,K+1)/H(K+1,K,II)
*
*           Form and store the information for the new Givens rotation
            DO I=1,K-1
              W1=C(I,II)*H(I,K,II)-S(I,II)*H(I+1,K,II)
              W2=S(I,II)*H(I,K,II)+C(I,II)*H(I+1,K,II)
              H(I,K,II)=W1
              H(I+1,K,II)=W2
            ENDDO
            ZNU=SQRT(H(K,K,II)**2+H(K+1,K,II)**2)
            IF(ZNU /= 0.0D0) THEN
              C(K,II)=H(K,K,II)/ZNU
              S(K,II)=-H(K+1,K,II)/ZNU
              H(K,K,II)=C(K,II)*H(K,K,II)-S(K,II)*H(K+1,K,II)
              H(K+1,K,II)=0.0D0
              W1=C(K,II)*G(K,II)-S(K,II)*G(K+1,II)
              W2=S(K,II)*G(K,II)+C(K,II)*G(K+1,II)
              G(K,II)=W1
              G(K+1,II)=W2
            ENDIF
*
*           Update the residual norm
            RHO(II)=ABS(G(K+1,II))
            REPSI(ITER,II)=REAL(RHO(II)/DENOM(II))
            IF(IPRINT.GT.4) WRITE(6,200) ITER,II,REPSI(ITER,II)
            EPS(II)=REPSI(ITER,II)
            ITST(II)=ITER
            IF(RHO(II) < ERRTOL*DENOM(II)) THEN
              NCONV(II)=.FALSE.
              LNCONV=LNCONV-1
            ENDIF
            ERROR=MAX(ERROR,RHO(II)/DENOM(II))
          ENDDO
        ENDDO
*---
*   At this point either K > NSTART or RHOGRP < ERRTOL. It's time to
*   compute PHIIN and cycle.
*---
        DO II=1,NGEFF
          K=KMAX(II)
          IF(K == 0) CYCLE
          IG=NGIND(II)
          G(K,II)=G(K,II)/H(K,K,II)
          DO L=K-1,1,-1
            W1=G(L,II)-DOT_PRODUCT(H(L,L+1:K,II),G(L+1:K,II))
            G(L,II)=W1/H(L,L,II)
          ENDDO
          DO J=1,K
            PHIIN(:,II)=PHIIN(:,II)+REAL(G(J,II)*V(:,II,J))
          ENDDO
        ENDDO
      ENDDO
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(R, S, C, H, G, V, GAR, RHS)
      DEALLOCATE(FLOUT, KMAX, RHO, DENOM)
      RETURN
*
  100 FORMAT(' MCGMRE: RHS CALCULATED WITH AAC-SCR : ',I2,1H-,I1)
  200 FORMAT(' MCGMRE: EPS Iteration ',I5,' Group ',I5,2X,1P,E16.7)
      END
