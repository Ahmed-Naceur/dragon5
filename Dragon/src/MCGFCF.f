*DECK MCGFCF
      SUBROUTINE MCGFCF(SUBFFI,SUBFFA,SUBLDC,SUBSCH,IFTRAK,NBTR,NMAX,
     1                  NDIM,KPN,K,NREG,M,NGEFF,NANGL,NMU,NLF,NFUNL,
     2                  NMOD,NLFX,NLIN,NFUNLX,KEYFLX,KEYCUR,NZON,NCONV,
     3                  CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,PHI,S,SIGAL,ISGNR,
     4                  IDIR,NSOUT,XSI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Flux integration upon the non-cyclic tracking.
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
* SUBFFI  flux integration subroutine with isotropic source.
* SUBFFA  flux integration subroutine with anisotropic source.
* SUBLDC  flux integration subroutine with linear-discontinuous source.
* SUBSCH  track coefficients calculation subroutine.
* IFTRAK  tracking file unit number.
* NBTR    total number of tracking lines.
* NMAX    maximum number of elements in a track.
* NDIM    number of dimensions for the geometry.
* KPN     total number of unknowns in vectors PHI.
* K       total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NREG    number of volumes.
* M       number of material mixtures.
* NGEFF   number of groups to process.
* NANGL   number of tracking angles in the tracking file.
* NMU     order of the polar quadrature in 2D / 1 in 3D.
* NLF     number of Legendre orders for the flux.
* NFUNL   number of moments of the flux (in 2D: NFUNL=NLF*(NLF+1)/2).
* NMOD    first dimension of ISGNR.
* NLFX    scattering anisotropy used to compute spherical harmonics.
* NLIN    linear discontinuous flag (=1 SC/DD0; =3 LDC/DD1).
* NFUNLX  number of spherical harmonics components.
* KEYFLX  position of flux elements in PHI vector.
* KEYCUR  position of current elements in PHI vector.
* NZON    index-number of the mixture type assigned to each volume.
* NCONV   logical array of convergence status for each group (.TRUE.
*         not converged).
* CAZ0    cosines of the tracking polar angles in 3D.
* CAZ1    first cosines of the different tracking azimuthal angles.
* CAZ2    second cosines of the different tracking azimuthal angles.
* CPO     cosines of the different tracking polar angles in 2D.
* ZMU     polar quadrature set in 2D.
* WZMU    polar quadrature set in 2D.
* S       total source vector components.
* SIGAL   total cross-section and albedo array.
* ISGNR   sign of correction.
* IDIR    direction of fundamental current for TIBERE with MoC 
*         (=0,1,2,3). 
* NSOUT   number of outer surfaces.
* XSI     x,y and z component of the shape parameter for TIBERE. 
*
*Parameters: input/output
* PHI     vector containing the zonal scalar flux.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGEFF,K,KPN,M,NMAX,NDIM,NMU,NZON(K),NLF,NFUNL,NMOD,
     1 NLFX,NLIN,NFUNLX,NREG,KEYFLX(NREG,NLIN,NFUNL),KEYCUR(K-NREG),
     2 IFTRAK,NBTR,NANGL,ISGNR(NMOD,NFUNLX),IDIR,NSOUT
      REAL CPO(NMU),ZMU(NMU),WZMU(NMU),SIGAL(-6:M,NGEFF)
      DOUBLE PRECISION CAZ0(NANGL),CAZ1(NANGL),CAZ2(NANGL),
     1 PHI(KPN,NGEFF),S(KPN,NGEFF),XSI(NSOUT)
      LOGICAL NCONV(NGEFF)
      EXTERNAL SUBFFI,SUBFFA,SUBLDC,SUBSCH
*----
*  LOCAL VARIABLES
*----
      INTEGER II,ILINE,IANG,NSEG,I,IMU,IANG0,NOMP,INDP,NOMM,INDM,
     1 NOMI,JF,IND,NSUB,INDX,INDY,IREG
      REAL XMUANG(1)
      DOUBLE PRECISION WEITF,WEIGHT,Q0,Q1,Q0X,Q1X,Q0Y,Q1Y,ZMUI,
     1 OMEGA2(3),ZZZ
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NOM
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: HTF,T2D
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RHARM
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: TRHAR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: B
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: COEFI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: PHIV,DPHIV
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: STOT,DSTOT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NOM(NMAX),HTF(NMAX))
*---
* Compute flux and currents for this tracking line
*---
      IF((NLF.EQ.1).AND.(NLIN.EQ.1)) THEN
*     --------------------
*     Isotropic Scattering
*     --------------------
         ALLOCATE(B(2*NMAX))
         IF(NDIM.EQ.3) THEN
*        ---
*        3D calculation -> no loop over the polar angle
*        ---
         DO ILINE=1,NBTR
            READ(IFTRAK) NSUB,NSEG,WEITF,IANG,(NOM(I),I=1,NSEG),
     1                                        (HTF(I),I=1,NSEG)
            IF(NSUB.NE.1) CALL XABORT('MCGFCF: NSUB.NE.1.')
            DO II=1,NGEFF
            IF(NCONV(II)) THEN
*                   MCGFFI:  'Source Term Isolation' Strategy turned on
*                   MCGFFIS: 'Source Term Isolation' Strategy turned off
*                   MCGFFIT: 'MOCC/MCI' Iterative Strategy
               OMEGA2(3)=CAZ0(IANG)*CAZ0(IANG)
               ZZZ=1.0D0/SQRT(1.0D0-OMEGA2(3))
               OMEGA2(1)=3.0D0*(CAZ1(IANG)/ZZZ)**2
               OMEGA2(2)=3.0D0*(CAZ2(IANG)/ZZZ)**2
               OMEGA2(3)=3.0D0*OMEGA2(3)
               CALL SUBFFI(SUBSCH,K,KPN,M,NSEG,HTF,NOM,NZON,SIGAL(0,II),
     1              S(1,II),NREG,KEYFLX,KEYCUR,PHI(1,II),B,WEITF,
     2              OMEGA2,IDIR,NSOUT,XSI)
            ENDIF
            ENDDO
         ENDDO
         ELSE
*        ---
*        2D calculation -> loop over the polar angle
*        ---
         ALLOCATE(T2D(NMAX))
         DO ILINE=1,NBTR
            READ(IFTRAK) NSUB,NSEG,WEITF,IANG,(NOM(I),I=1,NSEG),
     1                                        (HTF(I),I=1,NSEG)
            IF(NSUB.NE.1) CALL XABORT('MCGFCF: NSUB.NE.1.')
            DO IMU=1,NMU
               OMEGA2(1)=3.0D0*(CAZ1(IANG)/ZMU(IMU))**2
               OMEGA2(2)=3.0D0*(CAZ2(IANG)/ZMU(IMU))**2
               OMEGA2(3)=3.0D0*(1.0-1.0/ZMU(IMU)**2)
               ZMUI=ZMU(IMU)
               WEIGHT=WEITF*DBLE(WZMU(IMU))
               DO I=2,NSEG-1
                  T2D(I)=HTF(I)*ZMUI
               ENDDO
               DO II=1,NGEFF            
               IF(NCONV(II)) THEN
                  CALL SUBFFI(SUBSCH,K,KPN,M,NSEG,T2D,NOM,NZON,
     1                 SIGAL(0,II),S(1,II),NREG,KEYFLX,KEYCUR,PHI(1,II),
     2                 B,WEIGHT,OMEGA2,IDIR,NSOUT,XSI)
               ENDIF
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE(T2D)
         ENDIF
         DEALLOCATE(B)
      ELSE IF(NLIN.EQ.1) THEN
*     ----------------------
*     Anisotropic Scattering
*     ----------------------
         ALLOCATE(STOT(NMAX,NMU,NGEFF,2),B(2*NMAX))
         CALL XDDSET(STOT,2*NMU*NMAX*NGEFF,0.0D0)
         ALLOCATE(RHARM(NMU,NFUNL),TRHAR(NMU,NFUNL,NMOD))
         IANG0=0
         IF(NDIM.EQ.3) THEN
*        ---
*        3D calculation -> no loop over the polar angle
*        ---
         DO ILINE=1,NBTR
            READ(IFTRAK) NSUB,NSEG,WEITF,IANG,(NOM(I),I=1,NSEG),
     1                                        (HTF(I),I=1,NSEG)
            IF(NSUB.NE.1) CALL XABORT('MCGFCF: NSUB.NE.1.')
            IF(IANG.NE.IANG0) THEN
               IANG0=IANG
               XMUANG(1)=REAL(CAZ0(IANG))
               CALL MOCCHR(3,NLF-1,NFUNL,1,XMUANG,CAZ1(IANG),
     1                     CAZ2(IANG),RHARM)
               DO 10 JF=1,NFUNL
                  TRHAR(1,JF,1)=ISGNR(1,JF)*RHARM(1,JF)
                  TRHAR(1,JF,NMOD)=ISGNR(NMOD,JF)*RHARM(1,JF)
 10            CONTINUE
            ENDIF
            DO II=1,NGEFF            
            IF(NCONV(II)) THEN
*              incoming flux in + direction
               NOMP=NOM(1)
               INDP=KEYCUR(-NOMP)
*              incoming flux in - direction
               NOMM=NOM(NSEG)
               INDM=KEYCUR(-NOMM)
               STOT(1,1,II,1)=WEITF*S(INDP,II)
               STOT(NSEG,1,II,2)=WEITF*S(INDM,II)
*              regional sources   
               DO I=2,NSEG-1
                  NOMI=NOM(I)
                  Q0=0.0D0
                  Q1=0.0D0
                  DO JF=1,NFUNL
                     IND=KEYFLX(NOMI,1,JF)         
                     Q0=Q0+S(IND,II)*TRHAR(1,JF,1)
                     Q1=Q1+S(IND,II)*TRHAR(1,JF,NMOD)
                  ENDDO                       
                  STOT(I,1,II,1)=WEITF*Q0
                  STOT(I,1,II,2)=WEITF*Q1
               ENDDO
            ENDIF
            ENDDO
            DO II=1,NGEFF
            IF(NCONV(II)) THEN
*                   MCGFFA:  'Source Term Isolation' Strategy turned on
*                   MCGFFAS: 'Source Term Isolation' Strategy turned off
*                   MCGFFAT: 'MOCC/MCI' Iterative Strategy
               CALL SUBFFA(SUBSCH,K,KPN,M,NSEG,HTF,NOM,NZON,SIGAL(0,II),
     1              STOT(1,1,II,1),STOT(1,1,II,2),NREG,1,NLF,NFUNL,
     2              NMOD,TRHAR,KEYFLX,KEYCUR,1,PHI(1,II),B,1,NMOD)
            ENDIF
            ENDDO
         ENDDO
*        ---
         ELSE
*        ---
*        2D calculation -> loop over the polar angle
*        ---
         ALLOCATE(T2D(NMAX))
         DO ILINE=1,NBTR
            READ(IFTRAK) NSUB,NSEG,WEITF,IANG,(NOM(I),I=1,NSEG),
     1                                        (HTF(I),I=1,NSEG)
            IF(NSUB.NE.1) CALL XABORT('MCGFCF: NSUB.NE.1.')
            IF(IANG.NE.IANG0) THEN
               IANG0=IANG
               CALL MOCCHR(2,NLF-1,NFUNL,NMU,CPO,CAZ1(IANG),CAZ2(IANG),
     1                     RHARM)
               DO 25 JF=1,NFUNL
               DO 20 IMU=1,NMU
                  TRHAR(IMU,JF,1)=ISGNR(1,JF)*RHARM(IMU,JF)
                  TRHAR(IMU,JF,NMOD)=ISGNR(NMOD,JF)*RHARM(IMU,JF)
 20            CONTINUE
 25            CONTINUE
            ENDIF
            DO II=1,NGEFF            
            IF(NCONV(II)) THEN
*              incoming flux in + direction
               NOMP=NOM(1)
               INDP=KEYCUR(-NOMP)
*              incoming flux in - direction
               NOMM=NOM(NSEG)
               INDM=KEYCUR(-NOMM)
               DO IMU=1,NMU
                  WEIGHT=WEITF*DBLE(WZMU(IMU))
                  STOT(1,IMU,II,1)=WEIGHT*S(INDP,II)
                  STOT(NSEG,IMU,II,2)=WEIGHT*S(INDM,II)
               ENDDO
*              regional sources               
               DO I=2,NSEG-1
                  NOMI=NOM(I)
                  DO IMU=1,NMU
                     Q0=0.0D0
                     Q1=0.0D0
                     WEIGHT=WEITF*DBLE(WZMU(IMU))
                     DO JF=1,NFUNL
                        IND=KEYFLX(NOMI,1,JF)         
                        Q0=Q0+S(IND,II)*TRHAR(IMU,JF,1)
                        Q1=Q1+S(IND,II)*TRHAR(IMU,JF,NMOD)
                     ENDDO                       
                     STOT(I,IMU,II,1)=WEIGHT*Q0
                     STOT(I,IMU,II,2)=WEIGHT*Q1
                  ENDDO
               ENDDO
            ENDIF
            ENDDO
            DO IMU=1,NMU
               ZMUI=ZMU(IMU)
               WEIGHT=WEITF*DBLE(WZMU(IMU))
               DO I=2,NSEG-1
                  T2D(I)=HTF(I)*ZMUI
               ENDDO
               DO II=1,NGEFF            
               IF(NCONV(II)) THEN
                  CALL SUBFFA(SUBSCH,K,KPN,M,NSEG,T2D,NOM,NZON,
     1                 SIGAL(0,II),STOT(1,IMU,II,1),STOT(1,IMU,II,2),
     2                 NREG,NMU,NLF,NFUNL,NMOD,TRHAR,KEYFLX,KEYCUR,
     3                 IMU,PHI(1,II),B,1,NMOD)
               ENDIF
               ENDDO
            ENDDO
         ENDDO
*        ---
         DEALLOCATE(T2D)
         ENDIF
         DEALLOCATE(TRHAR,RHARM,B,STOT)
      ELSE IF(NLIN.EQ.3) THEN
*     -----------------------------------------
*     Linear discontinuous source approximation
*     -----------------------------------------
         ALLOCATE(RHARM(NMU,NFUNLX),TRHAR(NMU,NFUNLX,NMOD),B(6*NMAX))
         ALLOCATE(PHIV(NFUNLX,NREG,NGEFF),DPHIV(NDIM*NFUNLX,NREG,NGEFF))
         ALLOCATE(STOT(NMAX,NMU,NGEFF,2),DSTOT(NMAX,NMU,NGEFF,2))
         CALL XDDSET(STOT,2*NMU*NMAX*NGEFF,0.0D0)
         CALL XDDSET(DSTOT,2*NMU*NMAX*NGEFF,0.0D0)
         DO II=1,NGEFF            
           IF(NCONV(II)) THEN
             CALL XDDSET(PHIV(1,1,II),NFUNLX*NREG,0.0D0)
             CALL XDDSET(DPHIV(1,1,II),NDIM*NFUNLX*NREG,0.0D0)
           ENDIF
         ENDDO
         IF(NDIM.EQ.3) THEN
            CALL XABORT('MCGFCF: 3D LDC APPROXIMATION NOT IMPLEMENTED')
         ELSE
*        ---
*        2D calculation -> loop over the polar angle
*        ---
         ALLOCATE(T2D(NMAX))
         IANG0=0
         DO ILINE=1,NBTR
            READ(IFTRAK) NSUB,NSEG,WEITF,IANG,(NOM(I),I=1,NSEG),
     1                                        (HTF(I),I=1,NSEG)
            IF(NSUB.NE.1) CALL XABORT('MCGFCF: NSUB.NE.1.')
            IF(IANG.NE.IANG0) THEN
               IANG0=IANG
               CALL MOCCHR(2,NLFX-1,NFUNLX,NMU,CPO,CAZ1(IANG),
     1                     CAZ2(IANG),RHARM)
               DO 35 JF=1,NFUNLX
               DO 30 IMU=1,NMU
                  TRHAR(IMU,JF,1)=ISGNR(1,JF)*RHARM(IMU,JF)
                  TRHAR(IMU,JF,NMOD)=ISGNR(NMOD,JF)*RHARM(IMU,JF)
 30            CONTINUE
 35            CONTINUE
            ENDIF
            DO II=1,NGEFF            
            IF(NCONV(II)) THEN
*              incoming flux in + direction
               NOMP=NOM(1)
               INDP=KEYCUR(-NOMP)
*              incoming flux in - direction
               NOMM=NOM(NSEG)
               INDM=KEYCUR(-NOMM)
               DO IMU=1,NMU
                  STOT(1,IMU,II,1)=S(INDP,II)
                  STOT(NSEG,IMU,II,2)=S(INDM,II)
               ENDDO
*              regional sources               
               DO I=2,NSEG-1
                  NOMI=NOM(I)
                  DO IMU=1,NMU
                     Q0=0.0D0
                     Q1=0.0D0
                     Q0X=0.0D0
                     Q1X=0.0D0
                     Q0Y=0.0D0
                     Q1Y=0.0D0
                     DO JF=1,NFUNL
                        IND=KEYFLX(NOMI,1,JF)         
                        INDX=KEYFLX(NOMI,2,JF)         
                        INDY=KEYFLX(NOMI,3,JF)         
                        Q0=Q0+S(IND,II)*TRHAR(IMU,JF,1)
                        Q1=Q1+S(IND,II)*TRHAR(IMU,JF,NMOD)
                        Q0X=Q0X+S(INDX,II)*TRHAR(IMU,JF,1)
                        Q1X=Q1X+S(INDX,II)*TRHAR(IMU,JF,NMOD)
                        Q0Y=Q0Y+S(INDY,II)*TRHAR(IMU,JF,1)
                        Q1Y=Q1Y+S(INDY,II)*TRHAR(IMU,JF,NMOD)
                     ENDDO                       
                     STOT(I,IMU,II,1)=Q0
                     STOT(I,IMU,II,2)=Q1
                     DSTOT(I,IMU,II,1)=Q0X*CAZ1(IANG)+Q0Y*CAZ2(IANG)
                     DSTOT(I,IMU,II,2)=-Q1X*CAZ1(IANG)-Q1Y*CAZ2(IANG)
                  ENDDO
               ENDDO
            ENDIF
            ENDDO
            DO IMU=1,NMU
               ZMUI=ZMU(IMU)
               WEIGHT=WEITF*DBLE(WZMU(IMU))
               DO I=2,NSEG-1
                  T2D(I)=HTF(I)*ZMUI
               ENDDO
               DO II=1,NGEFF            
               IF(NCONV(II)) THEN
*                   MCGFFAL: 'Source Term Isolation' Strategy turned off
                  CALL SUBLDC(SUBSCH,K,KPN,M,NSEG,T2D,NOM,NZON,WEIGHT,
     1                 SIGAL(0,II),STOT(1,IMU,II,1),STOT(1,IMU,II,2),
     2                 DSTOT(1,IMU,II,1),DSTOT(1,IMU,II,2),NREG,NMU,
     3                 NLF,NFUNLX,NMOD,TRHAR,KEYCUR,IMU,B,1,NMOD,
     4                 PHI(1,II),PHIV(1,1,II),DPHIV(1,1,II))
               ENDIF
               ENDDO
            ENDDO
         ENDDO
         ALLOCATE(COEFI(2*NFUNLX,2*NFUNLX))
         CALL MCGCOEF(NFUNLX,NMU,ZMU,WZMU,NANGL,CAZ1,CAZ2,COEFI)
         DO II=1,NGEFF            
         IF(NCONV(II)) THEN
           DO IREG=1,NREG
             DPHIV(:,IREG,II)=MATMUL(COEFI,DPHIV(:,IREG,II))
             DO JF=1,NFUNL
               PHI(KEYFLX(IREG,1,JF),II)=PHIV(JF,IREG,II)
               PHI(KEYFLX(IREG,2,JF),II)=DPHIV(JF,IREG,II)
               PHI(KEYFLX(IREG,3,JF),II)=DPHIV(NFUNLX+JF,IREG,II)
             ENDDO
           ENDDO
         ENDIF
         ENDDO
         DEALLOCATE(COEFI)
*        ---
         DEALLOCATE(T2D)
         ENDIF
         DEALLOCATE(DSTOT,STOT,DPHIV,PHIV)
         DEALLOCATE(B,TRHAR,RHARM)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(NOM,HTF)
      RETURN
      END
