*DECK MOCFCF
      SUBROUTINE MOCFCF(SUBFFI,SUBFFA,SUBLDC,SUBSCH,IFTRAK,NBTR,MXSUB,
     1                  MXSEG,NUNKNO,NREG,NSOUT,NMAT,NALB,NGEFF,NPHI,
     2                  NGSS,NLF,NFUNL,NMOD,NLFX,NLIN,NFUNLX,KEYFLX,
     3                  MATALB,NCONV,SIGANG,CAZ1,CAZ2,XGSS,YGSS,WGSS,
     4                  PHIOUT,SOUR,ISGNR,IDIR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Flux integration upon the cyclic tracking.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy and R. Le Tellier
*
*Parameters: input
* SUBFFI  Isotropic flux integration subroutine.
* SUBFFA  Anisotropic flux integration subroutine.
* SUBLDC  flux integration subroutine with linear-discontinuous source.
* SUBSCH  Track coefficients calculation subroutine.
* IFTRAK  tracking file unit number.
* NGEFF   number of groups to process.
* NMAT    number of mixtures.
* NLF     number of Legendre orders for the flux.
* NREG    number of regions.
* NUNKNO  number of unknowns per energy group including spherical
*         harmonic terms and fundamental currents.
* NGSS    number of polar angles.
* NSOUT   number of surfaces.
* MXSUB   maximun number of subtracks in a track.
* MXSEG   maximun number of segments in a track.
* NBTR    number of tracks.
* NPHI    number of angles in the plane.
* NFUNL   number of moments of the flux (in 2D: NFUNL=NLF*(NLF+1)/2).
* NMOD    first dimension of ISGNR.
* NLFX    scattering anisotropy used to compute spherical harmonics.
* NLIN    linear discontinuous flag (=1 SC/DD0; =3 LDC/DD1).
* NFUNLX  number of spherical harmonics components.
* ISGNR   array of spherical harmonics signs.
* KEYFLX  position of flux elements in PHIIN vector.
* MATALB  mixture and albedo indices.
* WGSS    polar weights.
* XGSS    polar angle cosines.
* YGSS    polar angle sines.
* NALB    number of albedos.
* SIGANG  arrays of total cross-sections and albedos.
* CAZ1    first cosines of the different tracking azimuthal angles.
* CAZ2    second cosines of the different tracking azimuthal angles.
* SOUR    total source vector components.
* NCONV   logical array of convergence status for each group (.TRUE. for
*         not converged).
* IDIR    direction of fundamental current for TIBERE with MoC 
*         (=0,1,2,3). 
*
*Parameters: output
* PHIOUT  vector containing the zonal flux moments.
*
*-----------------------------------------------------------------------
*
      IMPLICIT  NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IFTRAK,NGEFF,NMAT,NLF,NREG,NUNKNO,NGSS,NSOUT,MXSUB,MXSEG,
     1 NBTR,NPHI,NFUNL,NMOD,NLFX,NLIN,NFUNLX,ISGNR(NMOD,NFUNLX),
     2 KEYFLX(NREG,NLIN,NFUNL),MATALB(-NSOUT:NREG),NALB,IDIR
      REAL WGSS(NGSS),XGSS(NGSS),YGSS(NGSS),
     1 SIGANG(-NALB:NMAT,NGEFF)
      DOUBLE PRECISION CAZ1(NPHI),CAZ2(NPHI),SOUR(NUNKNO,NGEFF),
     1 PHIOUT(NUNKNO,NGEFF)
      LOGICAL NCONV(NGEFF)
      EXTERNAL SUBFFI,SUBFFA,SUBLDC,SUBSCH
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NRSEG,KANGL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SEGLEN,DSIG
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RHARM
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRHAR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: COEFI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: OMEGAX,OMEGAY
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: OMG2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: PHIV,DPHIV
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: EXPT,EXP2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: FLM,FLP,CYM,CYP,
     1 DFLM,DFLP,CYM2,CYP2
*----
*  LOCAL VARIABLES
*----
      INTEGER MXE
      PARAMETER (MXE=64)
      INTEGER ILINE,IANG,ISUB,NSEG,I,IE,II,NOMI,NZI,IND,JF,NSUB,INDX,
     1 INDY,IREG
      DOUBLE PRECISION WEIGHT
      DOUBLE PRECISION DWEIG(MXE),Q0,Q1,Q0X,Q1X,Q0Y,Q1Y
      LOGICAL LNEW
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NRSEG(MXSEG),SEGLEN(MXSEG),KANGL(MXSUB),FLM(NGSS,MXSEG),
     1 FLP(NGSS,MXSEG),CYM(NGSS,MXSEG),CYP(NGSS,MXSEG),EXPT(NGSS*MXSEG))
      ALLOCATE(OMEGAX(NPHI,NGSS),OMEGAY(NPHI,NGSS),OMG2(NGSS,3))
*----
*  ISOTROPIC SCATTERING
*----
      CALL XDDSET(OMG2,3*NGSS,0.D0)
      IF((NLF.EQ.1).AND.(NLIN.EQ.1)) THEN
      ALLOCATE(EXP2(2*NGSS*MXSEG))
      DO IE=1,NGSS
        DO IANG=1,NPHI
          OMEGAX(IANG,IE)=CAZ1(IANG)/YGSS(IE)
          OMEGAY(IANG,IE)=CAZ2(IANG)/YGSS(IE)
        ENDDO   
        OMG2(IE,3)=3.0D0*(1.0-1.0/YGSS(IE)**2)
      ENDDO   
      DO ILINE=1,NBTR
         READ(IFTRAK) NSUB,NSEG,WEIGHT,(KANGL(I),I=1,NSUB),
     1   (NRSEG(I),I=1,NSEG),(SEGLEN(I),I=1,NSEG)
         IF(NSUB.GT.MXSUB) CALL XABORT('MOCFCF: MXSUB OVERFLOW.')
         DO IE=1,NGSS
            DWEIG(IE)=WEIGHT*WGSS(IE)
         ENDDO
         DO II=1,NGEFF
         IF(NCONV(II)) THEN
            ISUB=0
            LNEW=.TRUE.
            DO I=1,NSEG
               NOMI=NRSEG(I)
               NZI=MATALB(NOMI)
               IF(NZI.LE.0) THEN
                  LNEW=.TRUE.
                  DO IE=1,NGSS
                     FLM(IE,I)=0.D0
                     FLP(IE,I)=0.D0
                  ENDDO
               ELSE
                  IF(LNEW) THEN
                    ISUB=ISUB+1
                    LNEW=.FALSE.
                  ENDIF
                  IND=KEYFLX(NOMI,1,1)
                  IF(IDIR.EQ.0) THEN
                   DO IE=1,NGSS
                     FLM(IE,I)=DWEIG(IE)*SOUR(IND,II)
                     FLP(IE,I)=FLM(IE,I)
                   ENDDO
                  ELSEIF(IDIR.EQ.1) THEN
                   DO IE=1,NGSS
                     OMG2(IE,1)=3.0D0*OMEGAX(KANGL(ISUB),IE)**2
                     FLM(IE,I)=DWEIG(IE)*SOUR(IND,II)*OMG2(IE,1)
                     FLP(IE,I)=FLM(IE,I)
                   ENDDO
                  ELSEIF(IDIR.EQ.2) THEN
                   DO IE=1,NGSS
                     OMG2(IE,2)=3.0D0*OMEGAY(KANGL(ISUB),IE)**2
                     FLM(IE,I)=DWEIG(IE)*SOUR(IND,II)*OMG2(IE,2)
                     FLP(IE,I)=FLM(IE,I)
                   ENDDO
                  ELSEIF(IDIR.EQ.3) THEN
                   DO IE=1,NGSS
                     FLM(IE,I)=DWEIG(IE)*SOUR(IND,II)*OMG2(IE,3)
                     FLP(IE,I)=FLM(IE,I)
                   ENDDO
                  ENDIF
               ENDIF
            ENDDO
*                MOCFFI:  'Source Term Isolation' Strategy turned on
*                MOCFFIS: 'Source Term Isolation' Strategy turned off
*                MOCFFIT: 'MOCC/MCI' Iterative Strategy
            CALL SUBFFI(SUBSCH,NREG,NSOUT,NUNKNO,NMAT,NSEG,SEGLEN,NRSEG,
     1           NGSS,MATALB,SIGANG(-NALB,II),KEYFLX,YGSS,PHIOUT(1,II),
     2           EXPT,EXP2,FLM,FLP,CYM,CYP,IDIR,OMG2)
         ENDIF
         ENDDO
      ENDDO
      DEALLOCATE(OMG2,OMEGAY,OMEGAX)
      DEALLOCATE(EXP2)
*----
*  ANISOTROPIC SCATTERING
*----
      ELSE IF(NLIN.EQ.1) THEN
      ALLOCATE(EXP2(2*NGSS*MXSEG))
      ALLOCATE(RHARM(NGSS,NFUNL),TRHAR(NGSS,NFUNL,NPHI,2))
      DO IANG=1,NPHI
         CALL MOCCHR(2,NLF-1,NFUNL,NGSS,XGSS,CAZ1(IANG),CAZ2(IANG),
     1   RHARM)
         DO 15 JF=1,NFUNL
         DO 10 IE=1,NGSS
            TRHAR(IE,JF,IANG,1)=ISGNR(1,JF)*RHARM(IE,JF)
            TRHAR(IE,JF,IANG,2)=ISGNR(NMOD,JF)*RHARM(IE,JF)
 10      CONTINUE
 15      CONTINUE
      ENDDO
      DO ILINE=1,NBTR
         READ(IFTRAK) NSUB,NSEG,WEIGHT,(KANGL(I),I=1,NSUB),
     1   (NRSEG(I),I=1,NSEG),(SEGLEN(I),I=1,NSEG)
         IF(NSUB.GT.MXSUB) CALL XABORT('MOCFCF: MXSUB OVERFLOW.')
         DO IE=1,NGSS
            DWEIG(IE)=WEIGHT*WGSS(IE)
         ENDDO
         DO II=1,NGEFF
         IF(NCONV(II)) THEN
            CALL XDDSET(FLM,NGSS*MXSEG,0.D0)
            CALL XDDSET(FLP,NGSS*MXSEG,0.D0)
            ISUB=0
            LNEW=.TRUE.
            DO I=1,NSEG
               NOMI=NRSEG(I)
               NZI=MATALB(NOMI)
               IF(NZI.LE.0) THEN
                  LNEW=.TRUE.
               ELSE
                  IF(LNEW) THEN
                     ISUB=ISUB+1
                     LNEW=.FALSE.
                  ENDIF
                  DO IE=1,NGSS
                     Q0=0.D0
                     Q1=0.D0
                     DO JF=1,NFUNL
                        IND=KEYFLX(NOMI,1,JF)
                        Q0=Q0+SOUR(IND,II)*TRHAR(IE,JF,KANGL(ISUB),2)
                        Q1=Q1+SOUR(IND,II)*TRHAR(IE,JF,KANGL(ISUB),1)
                     ENDDO
                     FLM(IE,I)=DWEIG(IE)*Q0
                     FLP(IE,I)=DWEIG(IE)*Q1
                  ENDDO
               ENDIF
            ENDDO
            IF(ISUB.NE.NSUB) CALL XABORT('MOCFCF: NSUB INCONSISTENCY')
*                MOCFFA:  'Source Term Isolation' Strategy turned on
*                MOCFFAS: 'Source Term Isolation' Strategy turned off
*                MOCFFAT: 'MOCC/MCI' Iterative Strategy
            CALL SUBFFA(SUBSCH,NREG,NSOUT,NUNKNO,NMAT,NSEG,SEGLEN,
     1           NRSEG,NGSS,NFUNL,MATALB,SIGANG(-NALB,II),KEYFLX,YGSS,
     2           PHIOUT(1,II),EXPT,EXP2,FLM,FLP,CYM,CYP,NPHI,NSUB,KANGL,
     3           TRHAR)
         ENDIF
         ENDDO
      ENDDO
      DEALLOCATE(EXP2)
      DEALLOCATE(TRHAR,RHARM)
*----
*  LINEAR DISCONTINUOUS SOURCE APPROXIMATION
*----
      ELSE IF(NLIN.EQ.3) THEN
      ALLOCATE(PHIV(NFUNLX,NREG,NGEFF),DPHIV(2*NFUNLX,NREG,NGEFF))
      ALLOCATE(EXP2(5*NGSS*MXSEG))
      ALLOCATE(RHARM(NGSS,NFUNLX),TRHAR(NGSS,NFUNLX,NPHI,2))
      DO IANG=1,NPHI
         CALL MOCCHR(2,NLFX-1,NFUNLX,NGSS,XGSS,CAZ1(IANG),CAZ2(IANG),
     1               RHARM)
         DO 25 JF=1,NFUNLX
         DO 20 IE=1,NGSS
            TRHAR(IE,JF,IANG,1)=ISGNR(1,JF)*RHARM(IE,JF)
            TRHAR(IE,JF,IANG,2)=ISGNR(NMOD,JF)*RHARM(IE,JF)
 20      CONTINUE
 25      CONTINUE
      ENDDO
      DO II=1,NGEFF            
        IF(NCONV(II)) THEN
          CALL XDDSET(PHIV(1,1,II),NFUNLX*NREG,0.0D0)
          CALL XDDSET(DPHIV(1,1,II),2*NFUNLX*NREG,0.0D0)
        ENDIF
      ENDDO
      ALLOCATE(DFLM(NGSS,MXSEG),DFLP(NGSS,MXSEG),DSIG(MXSEG),
     1         CYM2(NGSS,MXSEG),CYP2(NGSS,MXSEG))
      DO ILINE=1,NBTR
         READ(IFTRAK) NSUB,NSEG,WEIGHT,(KANGL(I),I=1,NSUB),
     1   (NRSEG(I),I=1,NSEG),(SEGLEN(I),I=1,NSEG)
         IF(NSUB.GT.MXSUB) CALL XABORT('MOCFCF: MXSUB OVERFLOW.')
         DO IE=1,NGSS
            DWEIG(IE)=WEIGHT*WGSS(IE)
         ENDDO
         DO II=1,NGEFF
         IF(NCONV(II)) THEN
            CALL XDDSET(FLM,NGSS*MXSEG,0.D0)
            CALL XDDSET(FLP,NGSS*MXSEG,0.D0)
            CALL XDDSET(DFLM,NGSS*MXSEG,0.D0)
            CALL XDDSET(DFLP,NGSS*MXSEG,0.D0)
            ISUB=0
            IANG=0
            LNEW=.TRUE.
            DO I=1,NSEG
               NOMI=NRSEG(I)
               NZI=MATALB(NOMI)
               IF(NZI.LE.0) THEN
                  LNEW=.TRUE.
               ELSE
                  IF(LNEW) THEN
                     ISUB=ISUB+1
                     IANG=KANGL(ISUB)
                     LNEW=.FALSE.
                  ENDIF
                  DO IE=1,NGSS
                     Q0=0.D0
                     Q1=0.D0
                     Q0X=0.0D0
                     Q1X=0.0D0
                     Q0Y=0.0D0
                     Q1Y=0.0D0
                     DO JF=1,NFUNL
                        IND=KEYFLX(NOMI,1,JF)
                        INDX=KEYFLX(NOMI,2,JF)         
                        INDY=KEYFLX(NOMI,3,JF)         
                        Q0=Q0+SOUR(IND,II)*TRHAR(IE,JF,IANG,2)
                        Q1=Q1+SOUR(IND,II)*TRHAR(IE,JF,IANG,1)
                        Q0X=Q0X+SOUR(INDX,II)*TRHAR(IE,JF,IANG,2)
                        Q1X=Q1X+SOUR(INDX,II)*TRHAR(IE,JF,IANG,1)
                        Q0Y=Q0Y+SOUR(INDY,II)*TRHAR(IE,JF,IANG,2)
                        Q1Y=Q1Y+SOUR(INDY,II)*TRHAR(IE,JF,IANG,1)
                     ENDDO
                     FLM(IE,I)=Q0
                     FLP(IE,I)=Q1
                     DFLM(IE,I)=-Q0X*CAZ1(IANG)-Q0Y*CAZ2(IANG)
                     DFLP(IE,I)=Q1X*CAZ1(IANG)+Q1Y*CAZ2(IANG)
                  ENDDO
               ENDIF
            ENDDO
            IF(ISUB.NE.NSUB) CALL XABORT('MOCFCF: NSUB INCONSISTENCY')
*                MOCFFAL: 'Source Term Isolation' Strategy turned off
            CALL SUBLDC(SUBSCH,NREG,NSOUT,NMAT,NSEG,SEGLEN,NRSEG,NGSS,
     1           NFUNLX,MATALB,DWEIG,SIGANG(-NALB,II),YGSS,FLM,FLP,
     2           DFLM,DFLP,NPHI,NSUB,KANGL,TRHAR,PHIV(1,1,II),
     3           DPHIV(1,1,II),DSIG,EXPT,EXP2,CYM,CYP,CYM2,CYP2)
         ENDIF
         ENDDO
      ENDDO
      DEALLOCATE(CYP2,CYM2,DSIG,DFLP,DFLM)
      DEALLOCATE(EXP2)
      DEALLOCATE(TRHAR,RHARM)
      ALLOCATE(COEFI(2*NFUNLX,2*NFUNLX))
      CALL MCGCOEF(NFUNLX,NGSS,YGSS,WGSS,NPHI,CAZ1,CAZ2,COEFI)
      DO II=1,NGEFF            
      IF(NCONV(II)) THEN
        DO IREG=1,NREG
          DPHIV(:,IREG,II)=MATMUL(COEFI,DPHIV(:,IREG,II))
          DO JF=1,NFUNL
            PHIOUT(KEYFLX(IREG,1,JF),II)=PHIV(JF,IREG,II)
            PHIOUT(KEYFLX(IREG,2,JF),II)=DPHIV(JF,IREG,II)
            PHIOUT(KEYFLX(IREG,3,JF),II)=DPHIV(NFUNLX+JF,IREG,II)
          ENDDO
        ENDDO
      ENDIF
      ENDDO
      DEALLOCATE(COEFI,DPHIV,PHIV)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(EXPT,CYP,CYM,FLP,FLM,KANGL,SEGLEN,NRSEG)
*
      RETURN
      END
