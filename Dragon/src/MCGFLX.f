*DECK MCGFLX
      SUBROUTINE MCGFLX(SUBFFI,SUBFFA,SUBSCH,SUBLDC,CYCLIC,KPSYS,IPRINT,
     1           IPTRK,IFTRAK,IPMACR,NDIM,K,KPN,NLONG,NREG,NSOUT,NG,
     2           NGEFF,NGIND,NZON,MATALB,V,FIMEM,QFR,M,NANI,MAXI,IAAC,
     3           KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,CAZ0,CAZ1,CAZ2,CPO,
     4           ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEYFLX,KEYCUR,SIGAL,LPS,
     5           REPS,EPS,ITST,NCONV,LNCONV,REBFLG,STIS,NPJJM,LPRISM,
     6           N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* MOC solution of the transport equation in 2D,3D-irregular geometry.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and R. Le Tellier
*
*Parameters: input/output
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
* K       total number of volumes-surfaces for which specific values
*         of the neutron flux and reactions rates are required.
* KPN     total number of unknowns in vectors SUNKNO and FUNKNO.
* NLONG   number of spatial unknowns.
* NREG    number of regions (volumes).
* NSOUT   number of outer surfaces.
* NG      number of groups.
* NGEFF   number of groups to process.
* NGIND   index of the groups to process.
* NZON    mixture-albedo index array in MCCG format.
* MATALB  albedo-mixture index array in MOCC format.
* V       volumes and surfaces.
* FIMEM   unknown vector.
* QFR     input source vector.
* M       number of material mixtures.
* NANI    scattering anisotropy (=1 for isotropic scattering).
* MAXI    maximum number of inner iterations.
* IAAC    no acceleration / CDD acceleration of inner iterations (0/1).
* KRYL    Bi-CGSTAB scheme used / GMRES scheme not used / GMRES scheme
*         used (Krylov subspace dimension = KRYL) (<0 / 0 / >0).
* ISCR    no acceleration / SCR acceleration of inner iterations (0/1).
* NMU     order of the polar quadrature in 2D / 1 in 3D.
* NANGL   number of tracking angles in the plane.
* NMAX    maximum number of elements in a track.
* LC      dimension of MCU vector.
* EPSI    tolerance of inner iterations.
* CAZ0    cosines of the tracking polar angles in 3D.
* CAZ1    first cosines of the different tracking azimuthal angles.
* CAZ2    second cosines of the different tracking azimuthal angles.
* CPO     cosines of the different tracking polar angles in 2D.
* ZMU     polar quadrature set in 2D.
* WZMU    polar quadrature set in 2D.
* LFORW   flag set to .false. to transpose the coefficient matrix.
* PACA    type of preconditioner to solve the ACA corrective system.
* NLIN    number of polynomial components in flux spatial expansion.
* NFUNL   number of moments of the flux (in 2D: NFUNL=NANI*(NANI+1)/2).
* KEYFLX  position of flux elements in FUNKNO vector.
* KEYCUR  position of current elements in FUNKNO vector.
* SIGAL   total cross-section and albedo array.
* LPS     dimension of PSJ vector for SCR acceleration.
* EPS     array of the precision reached for each group.
* REPS    array of the precision for each iteration and each group.
* ITST    array of the number of iterations for each group.
* NCONV   array of convergence flag for each group.
* LNCONV  number of unconverged groups.
* REBFLG  ACA or SCR rebalancing flag.
* STIS    'Source term isolation' option for flux integration.
* NPJJM   number of pjj modes to store for STIS option.
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
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPMACR,KPSYS(NGEFF)
      INTEGER NGEFF,IPRINT,IFTRAK,NDIM,K,KPN,NLONG,NREG,NSOUT,NG,
     1 NGIND(NGEFF),NZON(NLONG),MATALB(-NSOUT:NREG),M,NANI,MAXI,
     2 IAAC,KRYL,NMU,NANGL,NMAX,LC,PACA,ISCR,LPS,ITST(NGEFF),LNCONV,
     3 NLIN,NFUNL,KEYFLX(NREG,NLIN,NFUNL),KEYCUR(NLONG-NREG),STIS,
     4 NPJJM,N2REG,N2SOU,NZP,IDIR
      REAL V(NLONG),FIMEM(KPN,NGEFF),QFR(KPN,NG),EPSI,CPO(NMU),ZMU(NMU),
     1 WZMU(NMU),SIGAL(-6:M,NGEFF),REPS(MAXI,NGEFF),EPS(NGEFF),DELU,
     2 FACSYM
      DOUBLE PRECISION CAZ0(NANGL),CAZ1(NANGL),CAZ2(NANGL)
      LOGICAL LFORW,CYCLIC,NCONV(NGEFF),REBFLG,LPRISM
      EXTERNAL SUBFFI,SUBFFA,SUBLDC,SUBSCH
*----
*  ALLOCATABLE ARRAYS
*----
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: SOUR,FLUX
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IUNOUT=6)
      REAL MAXDIF,MAXFL
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SOUR(KPN,NGEFF),FLUX(KPN,NGEFF))
      CALL XDDSET(SOUR,KPN*NGEFF,0.0D0)
*---
      IF(IPRINT.GT.5) THEN
         DO II=1,NGEFF
            WRITE(6,*) 'GROUP(',NGIND(II),')'
            CALL PRINAM('FI-0  ',FIMEM(1,II),KPN)
         ENDDO
      ENDIF
*----
*  INNER ITERATIONS FOR THE TRANSPORT SOLUTION
*----
      IF (KRYL.EQ.0) THEN
*     ---------------------------
*     Richardson Iterative Scheme
*     ---------------------------
         EPSINT=0.1*EPSI
         MAXINT=200
         ITER=0
         DO WHILE ((LNCONV.GT.0).AND.(ITER.LT.MAXI))
            ITER=ITER+1
            CALL MCGFL1(SUBFFI,SUBFFA,SUBLDC,SUBSCH,CYCLIC,KPSYS,IPRINT,
     1           IPTRK,IFTRAK,IPMACR,NDIM,K,KPN,NLONG,FLUX,NZON,MATALB,
     2           M,NANI,NMU,NMAX,NANGL,NREG,NSOUT,NG,NGEFF,NGIND,SOUR,
     3           IAAC,ISCR,LC,LFORW,PACA,EPSINT,MAXINT,NLIN,NFUNL,
     4           KEYFLX,KEYCUR,QFR,FIMEM,CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,V,
     5           SIGAL,LPS,NCONV,.TRUE.,(MAXI.EQ.1),STIS,NPJJM,REBFLG,
     6           LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
*           residual calculation and update NCONV
            DO II=1,NGEFF
               IF (NCONV(II)) THEN
                  IG=NGIND(II)
                  IF (MAXI.GT.1) THEN
                     MAXFL=0.0
                     MAXDIF=0.0
                     DO I=1,KPN
                        TEMP=REAL(ABS(FLUX(I,II)))
                        MAXFL=MAX(TEMP,MAXFL)
                        FIMEM(I,II)=FIMEM(I,II)-REAL(FLUX(I,II))
                        MAXDIF=MAX(ABS(FIMEM(I,II)),MAXDIF)
                     ENDDO
                     IF (MAXFL.EQ.0.0) MAXFL=1.0
                     REPS(ITER,II)=MAXDIF/MAXFL
                     IF (ITER.GT.2) THEN
                     IF (REPS(ITER,II).GT.REPS(ITER-1,II)) THEN
*                    preconditioning cutoff
                        IAAC=0
                        ISCR=0
                     ENDIF
                     ENDIF
                     IF ((REPS(ITER,II).LT.EPSI).OR.(ITER.EQ.MAXI)) THEN
                        NCONV(II)=.FALSE.
                        ITST(II)=ITER
                        EPS(II)=REPS(ITER,II)
                        LNCONV=LNCONV-1
                     ENDIF
                  ENDIF
                  DO I=1,KPN
                     FIMEM(I,II)=REAL(FLUX(I,II))
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ELSEIF(KRYL.GT.0) THEN
*     ----------------------
*     GMRES Iterative Scheme
*     ----------------------
         CALL MCGMRE(SUBFFI,SUBFFA,SUBLDC,SUBSCH,CYCLIC,KPSYS,IPRINT,
     1        IPTRK,IFTRAK,IPMACR,NDIM,K,KPN,NLONG,FLUX,NZON,MATALB,
     2        M,NANI,NMU,NMAX,NANGL,NREG,NSOUT,SOUR,IAAC,ISCR,LC,LFORW,
     3        PACA,ITST,MAXI,QFR,FIMEM,CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,V,
     4        EPS,EPSI,REPS,KRYL,SIGAL,LPS,NG,NGEFF,NGIND,NCONV,LNCONV,
     5        NLIN,NFUNL,KEYFLX,KEYCUR,STIS,NPJJM,REBFLG,LPRISM,N2REG,
     6        N2SOU,NZP,DELU,FACSYM,IDIR)
      ELSE
*     --------------------------
*     Bi-CGSTAB Iterative Scheme
*     --------------------------
         CALL MCGBIC(SUBFFI,SUBFFA,SUBLDC,SUBSCH,CYCLIC,KPSYS,IPRINT,
     1        IPTRK,IFTRAK,IPMACR,NDIM,K,KPN,NLONG,FLUX,NZON,MATALB,
     2        M,NANI,NMU,NMAX,NANGL,NREG,NSOUT,SOUR,IAAC,ISCR,LC,LFORW,
     3        PACA,ITST,MAXI,QFR,FIMEM,CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,V,
     4        EPS,EPSI,REPS,SIGAL,LPS,NG,NGEFF,NGIND,NCONV,LNCONV,NLIN,
     5        NFUNL,KEYFLX,KEYCUR,STIS,NPJJM,REBFLG,LPRISM,N2REG,N2SOU,
     6        NZP,DELU,FACSYM,IDIR)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUX,SOUR)
      RETURN
      END
