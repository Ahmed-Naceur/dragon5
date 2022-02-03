*DECK MCGBIC
      SUBROUTINE  MCGBIC(SUBFFI,SUBFFA,SUBDLC,SUBSCH,CYCLIC,KPSYS,
     1            IPRINT,IPTRK,IFTRAK,IPMACR,NDIM,K,N,NLONG,PHIOUT,
     2            NZON,MATALB,M,NANI,NMU,NMAX,NANGL,NREG,NSOUT,SOUR,
     3            IAAC,ISCR,LC,LFORW,PACA,ITST,MAXI,QFR,PHIIN,CAZ0,
     4            CAZ1,CAZ2,CPO,ZMU,WZMU,V,EPS,EPSI,REPSI,SIGAL,LPS,
     5            NG,NGEFF,NGIND,NCONV,LNCONV,NLIN,NFUNL,KEYFLX,KEYCUR,
     6            STIS,NPJJM,REBFLG,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,
     7            IDIR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve the linear system obtained by the characteristics formalism
* with BiCGSTAB iterative approach.
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
*Parameters:
* SUBFFI  flux integration subroutine with isotropic source.
* SUBFFA  flux integration subroutine with anisotropic source.
* SUBDLC  flux integration subroutine with linear-discontinuous source.
* SUBSCH  track coefficients calculation subroutine.
* CYCLIC  cyclic tracking flag.
* KPSYS   pointer array for each group properties.
* IPRINT  print parameter (equal to zero for no print).
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  tracking file unit number.
* IPMACR  pointer to the macrolib LCM object.
* NDIM    number of dimensions for the geometry.
* K       total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* N       total number of unknowns in vectors SUNKNO and FUNKNO.
* NLONG   number of spatial unknowns.
* PHIOUT  output flux vector.
* NZON    mixture-albedo index array in MCCG format.
* MATALB  albedo-mixture index array in MOCC format.
* M       number of material mixtures.
* NANI    scattering anisotropy (=1 for isotropic scattering).
* NMU     order of the polar quadrature set.
* NMAX    maximum number of elements in a track.
* NANGL   number of tracking angles in the plan.
* NREG    number of regions (volumes).
* NSOUT   number of outer surfaces.
* SOUR    scratch.
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
* V       volumes.
* EPS     precision reached after min(MAXI,ITST) iterations.
* EPSI    tolerance for stopping criterion. Process is stopped
*         as soon as: ||Phi(n+1)-Phi(n)||/||Phi(n)|| <= eps
*         with ||.|| the euclidean norm.
* REPSI   array containing precision of each iteration.
* SIGAL   total cross-section and albedo array.
* LPS     used in scr acceleration.
* NG      number of groups.
* NGEFF   number of groups to process.
* NGIND   index of the groups to process.
* NCONV   array of convergence flag for each group.
* LNCONV  number of unconverged groups.
* NLIN    number of polynomial components in flux spatial expansion.
* NFUNL   number of moments of the flux (in 2D: NFUNL=NANI*(NANI+1)/2).
* KEYFLX  position of flux elements in FUNKNO vector.
* KEYCUR  position of current elements in FUNKNO vector.
* STIS    'Source term isolation' option for flux integration.
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
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KPSYS(NGEFF),IPTRK,IPMACR
      INTEGER NGEFF,IPRINT,IFTRAK,NDIM,K,N,NLONG,NG,NZON(NLONG),M,
     1 NANI,NMU,NMAX,IAAC,LC,PACA,ITST(NGEFF),MAXI,NANGL,NREG,NSOUT,
     2 ISCR,LPS,NGIND(NGEFF),LNCONV,NLIN,NFUNL,KEYFLX(NREG,NLIN,NFUNL),
     3 KEYCUR(NLONG-NREG),STIS,NPJJM,MATALB(-NSOUT:NREG),N2REG,N2SOU,
     4 NZP,IDIR
      REAL QFR(N,NG),PHIIN(N,NG),CPO(NMU),ZMU(NMU),WZMU(NMU),V(NLONG),
     1 EPS(NGEFF),EPSI,REPSI(MAXI,NGEFF),SIGAL(-6:M,NGEFF),DELU,FACSYM
      DOUBLE PRECISION CAZ0(NANGL),CAZ1(NANGL),CAZ2(NANGL),
     1 PHIOUT(N,NGEFF),SOUR(N,NGEFF)
      LOGICAL LFORW,CYCLIC,NCONV(NGEFF),REBFLG,LPRISM
      EXTERNAL SUBFFI,SUBFFA,SUBDLC,SUBSCH
*---
*   LOCAL VARIABLES
*---
      INTEGER J,II,IG,ITER,MAXINT
      REAL R,BI,WI,RT1,EPSINT,REPSMAX
      REAL SDOT
      DOUBLE PRECISION DDOT
      LOGICAL RHSFLG
      INTRINSIC SQRT,ABS
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RHS,PI,RI,SI,ROT,API,AUX
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(RHS(N,NGEFF),PI(N,NGEFF),RI(N,NGEFF),SI(N,NGEFF),
     1 ROT(N,NGEFF),API(N,NGEFF),AUX(2,NGEFF))
*---
      IF(MAXI.LT.4) CALL XABORT('MCGBIC: MAXI MUST BE >= 4.')
      MAXINT=MAXI-1
      EPSINT=0.01*EPSI
      RHSFLG=.TRUE.
*---
      ITER=1
*     a first iteration
      CALL MCGFL1(SUBFFI,SUBFFA,SUBDLC,SUBSCH,CYCLIC,KPSYS,IPRINT,IPTRK,
     1     IFTRAK,IPMACR,NDIM,K,N,NLONG,PHIOUT,NZON,MATALB,M,NANI,NMU,
     2     NMAX,NANGL,NREG,NSOUT,NG,NGEFF,NGIND,SOUR,IAAC,ISCR,LC,LFORW,
     3     PACA,EPSINT,MAXINT,NLIN,NFUNL,KEYFLX,KEYCUR,QFR,PHIIN,CAZ0,
     4     CAZ1,CAZ2,CPO,ZMU,WZMU,V,SIGAL,LPS,NCONV,.TRUE.,.FALSE.,
     5     STIS,NPJJM,REBFLG,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
      DO II=1,NGEFF
         IG=NGIND(II)
         IF (NCONV(II)) THEN
            DO J=1,N
               RI(J,II)=REAL(PHIOUT(J,II))-PHIIN(J,IG)
               PHIIN(J,IG)=REAL(PHIOUT(J,II))
            ENDDO
            R=SDOT(N,RI(1,II),1,RI(1,II),1)
            REPSI(ITER,II)=REAL(SQRT(R/DDOT(N,PHIOUT(1,II),1,
     1                     PHIOUT(1,II),1)))
            IF (REPSI(ITER,II).LE.EPSI) THEN 
               NCONV(II)=.FALSE.
               ITST(II)=ITER
               EPS(II)=REPSI(ITER,II)
               LNCONV=LNCONV-1
            ENDIF
            IF (LNCONV.EQ.0) RETURN
         ENDIF
      ENDDO
*---
 10   ITER=ITER+1
*     compute initial residual vector
      CALL MCGFL1(SUBFFI,SUBFFA,SUBDLC,SUBSCH,CYCLIC,KPSYS,IPRINT,IPTRK,
     1     IFTRAK,IPMACR,NDIM,K,N,NLONG,PHIOUT,NZON,MATALB,M,NANI,NMU,
     2     NMAX,NANGL,NREG,NSOUT,NG,NGEFF,NGIND,SOUR,IAAC,ISCR,LC,LFORW,
     3     PACA,EPSINT,MAXINT,NLIN,NFUNL,KEYFLX,KEYCUR,QFR,PHIIN,CAZ0,
     4     CAZ1,CAZ2,CPO,ZMU,WZMU,V,SIGAL,LPS,NCONV,.TRUE.,.FALSE.,
     5     STIS,NPJJM,REBFLG,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
      DO II=1,NGEFF
         IG=NGIND(II)
         IF (NCONV(II)) THEN
            DO J=1,N
               RI(J,II)=REAL(PHIOUT(J,II))-PHIIN(J,IG)
            ENDDO
            R=SDOT(N,RI(1,II),1,RI(1,II),1)
            REPSI(ITER,II)=REAL(SQRT(R/DDOT(N,PHIOUT(1,II),1,
     1                     PHIOUT(1,II),1)))
            IF (REPSI(ITER,II).LE.EPSI) THEN 
               NCONV(II)=.FALSE.
               ITST(II)=ITER
               EPS(II)=REPSI(ITER,II)
               LNCONV=LNCONV-1
               DO J=1,N
                  PHIIN(J,IG)=REAL(PHIOUT(J,II))
               ENDDO
            ENDIF
            IF (LNCONV.EQ.0) RETURN
            DO J=1,N
               PI(J,II)=RI(J,II)
               ROT(J,II)=RI(J,II)
            ENDDO
*            RT2=R !!SDOT(N,RI,1,ROT,1)
            AUX(1,II)=R !!SDOT(N,RI(1,II),1,ROT(1,II),1)
         ENDIF
      ENDDO
*---
      IF (RHSFLG) THEN
*     evaluate RHS of the linear system
         IF (IPRINT.GT.3) THEN
            WRITE(6,100) IAAC,ISCR
         ENDIF
         CALL XDRSET(RHS,(N*NGEFF),0.0)
         CALL MCGFL1(SUBFFI,SUBFFA,SUBDLC,SUBSCH,CYCLIC,KPSYS,IPRINT,
     1        IPTRK,IFTRAK,IPMACR,NDIM,K,N,NLONG,PHIOUT,NZON,MATALB,
     2        M,NANI,NMU,NMAX,NANGL,NREG,NSOUT,NG,NGEFF,NGIND,SOUR,IAAC,
     2        ISCR,LC,LFORW,PACA,EPSINT,MAXINT,NLIN,NFUNL,KEYFLX,KEYCUR,
     3        QFR,RHS,CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,V,SIGAL,LPS,NCONV,
     4        .TRUE.,.FALSE.,STIS,NPJJM,REBFLG,LPRISM,N2REG,N2SOU,NZP,
     5        DELU,FACSYM,IDIR)
         DO II=1,NGEFF
            IF (NCONV(II)) THEN
               DO J=1,N
                  RHS(J,II)=REAL(PHIOUT(J,II))
               ENDDO
            ENDIF
         ENDDO
         RHSFLG=.FALSE.
      ENDIF
*
      CALL MCGFL1(SUBFFI,SUBFFA,SUBDLC,SUBSCH,CYCLIC,KPSYS,IPRINT,IPTRK,
     1     IFTRAK,IPMACR,NDIM,K,N,NLONG,PHIOUT,NZON,MATALB,M,NANI,NMU,
     2     NMAX,NANGL,NREG,NSOUT,NG,NGEFF,NGIND,SOUR,IAAC,ISCR,LC,LFORW,
     3     PACA,EPSINT,MAXINT,NLIN,NFUNL,KEYFLX,KEYCUR,QFR,PI,CAZ0,CAZ1,
     4     CAZ2,CPO,ZMU,WZMU,V,SIGAL,LPS,NCONV,.TRUE.,.FALSE.,STIS,
     5     NPJJM,REBFLG,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
      ITER=ITER+1
      DO II=1,NGEFF
         IF (NCONV(II)) THEN
            DO J=1,N
               API(J,II)=PI(J,II)-REAL(PHIOUT(J,II))+RHS(J,II)
            ENDDO
            REPSI(ITER,II)=REPSI((ITER-1),II)
         ENDIF
      ENDDO

*     
      DO WHILE (ITER.LT.(MAXI-1))
*     BiCGSTAB iterations
         ITER=ITER+1
*
         DO II=1,NGEFF
            IF (NCONV(II)) THEN
               AUX(2,II)=AUX(1,II)/SDOT(N,API(1,II),1,ROT(1,II),1)
               DO J=1,N
                  SI(J,II)=RI(J,II)-AUX(2,II)*API(J,II)
               ENDDO
            ENDIF
         ENDDO
*
         CALL MCGFL1(SUBFFI,SUBFFA,SUBDLC,SUBSCH,CYCLIC,KPSYS,IPRINT,
     1        IPTRK,IFTRAK,IPMACR,NDIM,K,N,NLONG,PHIOUT,NZON,MATALB,
     2        M,NANI,NMU,NMAX,NANGL,NREG,NSOUT,NG,NGEFF,NGIND,SOUR,IAAC,
     3        ISCR,LC,LFORW,PACA,EPSINT,MAXINT,NLIN,NFUNL,KEYFLX,KEYCUR,
     4        QFR,SI,CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,V,SIGAL,LPS,NCONV,
     5        .TRUE.,.FALSE.,STIS,NPJJM,REBFLG,LPRISM,N2REG,N2SOU,NZP,
     6        DELU,FACSYM,IDIR)
         REPSMAX=0.0
         DO II=1,NGEFF
            IG=NGIND(II)
            IF (NCONV(II)) THEN
               DO J=1,N
                  RI(J,II)=SI(J,II)-REAL(PHIOUT(J,II))+RHS(J,II)
               ENDDO
               WI=SDOT(N,RI(1,II),1,SI(1,II),1)/
     1            SDOT(N,RI(1,II),1,RI(1,II),1)
               DO J=1,N
                  PHIIN(J,IG)=PHIIN(J,IG)+AUX(2,II)*PI(J,II)+WI*SI(J,II)
                  RI(J,II)=SI(J,II)-WI*RI(J,II)
               ENDDO
               R=SDOT(N,RI(1,II),1,RI(1,II),1)
               REPSI(ITER,II)=SQRT(R/SDOT(N,PHIIN(1,IG),1,
     1                                      PHIIN(1,IG),1))
               REPSMAX=MAX(REPSMAX,REPSI(ITER,II))
               IF (REPSI(ITER,II).LE.EPSI) THEN 
                  NCONV(II)=.FALSE.
                  ITST(II)=ITER
                  EPS(II)=REPSI(ITER,II)
                  LNCONV=LNCONV-1
               ENDIF
               IF (LNCONV.EQ.0) GO TO 20
               RT1=AUX(1,II)
               AUX(1,II)=SDOT(N,RI(1,II),1,ROT(1,II),1)
               BI=AUX(1,II)/RT1*AUX(2,II)/WI
               DO J=1,N
                  PI(J,II)=RI(J,II)+BI*(PI(J,II)-WI*API(J,II))
               ENDDO
            ENDIF
         ENDDO
*
         CALL MCGFL1(SUBFFI,SUBFFA,SUBDLC,SUBSCH,CYCLIC,KPSYS,IPRINT,
     1        IPTRK,IFTRAK,IPMACR,NDIM,K,N,NLONG,PHIOUT,NZON,MATALB,
     2        M,NANI,NMU,NMAX,NANGL,NREG,NSOUT,NG,NGEFF,NGIND,SOUR,IAAC,
     3        ISCR,LC,LFORW,PACA,EPSINT,MAXINT,NLIN,NFUNL,KEYFLX,KEYCUR,
     4        QFR,PI,CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,V,SIGAL,LPS,NCONV,
     5        .TRUE.,.FALSE.,STIS,NPJJM,REBFLG,LPRISM,N2REG,N2SOU,NZP,
     6        DELU,FACSYM,IDIR)
         ITER=ITER+1
         DO II=1,NGEFF
            IF (NCONV(II)) THEN
               DO J=1,N
                  API(J,II)=PI(J,II)-REAL(PHIOUT(J,II))+RHS(J,II)
               ENDDO
               REPSI(ITER,II)=REPSI((ITER-1),II)
            ENDIF
         ENDDO
      ENDDO
*     
 20   CONTINUE
*     determine final residual norm
      ITER=ITER+1
      DO II=1,NGEFF
         IF (NCONV(II)) THEN
            ITST(II)=ITER
         ELSE
            IF (ITST(II).NE.1) THEN
               NCONV(II)=.TRUE.
               ITST(II)=ITST(II)+1
            ENDIF
         ENDIF
      ENDDO
      CALL MCGFL1(SUBFFI,SUBFFA,SUBDLC,SUBSCH,CYCLIC,KPSYS,IPRINT,IPTRK,
     1     IFTRAK,IPMACR,NDIM,K,N,NLONG,PHIOUT,NZON,MATALB,M,NANI,NMU,
     2     NMAX,NANGL,NREG,NSOUT,NG,NGEFF,NGIND,SOUR,0,0,LC,LFORW,PACA,
     3     EPSINT,MAXINT,NLIN,NFUNL,KEYFLX,KEYCUR,QFR,PHIIN,CAZ0,CAZ1,
     4     CAZ2,CPO,ZMU,WZMU,V,SIGAL,LPS,NCONV,.TRUE.,.FALSE.,STIS,
     5     NPJJM,REBFLG,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
      DO II=1,NGEFF
         IG=NGIND(II)
         IF (NCONV(II)) THEN
            DO J=1,N
               RI(J,II)=REAL(PHIOUT(J,II))-PHIIN(J,IG)
            ENDDO
            R=SDOT(N,RI(1,II),1,RI(1,II),1)
            REPSI(ITST(II),II)=SQRT(R/SDOT(N,PHIIN(1,IG),1,
     1                                       PHIIN(1,IG),1))
            DO J=1,N
               PHIIN(J,IG)=REAL(PHIOUT(J,II))
            ENDDO
            EPS(II)=REPSI(ITST(II),II)
         ENDIF
      ENDDO
      LNCONV=0
      IF (ITER.LT.MAXI) THEN
         DO II=1,NGEFF
            IF (EPS(II).GT.EPSI) THEN
               IF ((IAAC.GT.0).OR.(ISCR.GT.0)) THEN
                  IAAC=0
                  ISCR=0
                  RHSFLG=.TRUE.
               ENDIF
               IF (IPRINT.GT.2) WRITE(6,200) ITER
               NCONV(II)=.TRUE.
               LNCONV=LNCONV+1
            ELSE
               NCONV(II)=.FALSE.
            ENDIF
         ENDDO
         IF (LNCONV.GT.0) GO TO 10
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(AUX,API,ROT,SI,RI,PI,RHS)
      RETURN
*
 100  FORMAT(31H RHS CALCULATED WITH AAC-SCR : ,I1,1H-,I1)
 200  FORMAT(37H WARNING : BAD PREVISION, RESTART AT ,I4,10H ITERATION)
      END
