*DECK MCGPTF
      SUBROUTINE MCGPTF(SUBFFI,SUBFFA,SUBSCH,IFTRAK,N2BTR,N2MAX,
     1                  KPN,K,NREG,M,NGEFF,NANGL,NMU,NANI,NFUNL,NMOD,
     2                  KEYFLX,KEYCUR,NZON,NCONV,CAZ1,CAZ2,XMU,WZMU,PHI,
     3                  S,SIGAL,ISGNR,NMAX,NZP,N2REG,N2SOU,DELU,INDREG,
     4                  Z,VNORF,CMU,CMUI,SMU,SMUI,TMU,TMUI,SSYM,IDIR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Flux integration upon the tracking (3D prismatic extended tracking).
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
* SUBFFI  isotropic flux integration subroutine.
* SUBFFA  anisotropic flux integration subroutine.
* SUBSCH  track coefficients calculation subroutine.
* IFTRAK  tracking file unit number.
* N2BTR   total number of 2D tracking lines.
* N2MAX   maximum number of elements in a 2D track.
* KPN     total number of unknowns in vectors PHI.
* K       total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NREG    number of volumes.
* M       number of material mixtures.
* NGEFF   number of groups to process.
* NANGL   number of tracking angles in the tracking file.
* NMU     order of the polar quadrature in 2D / 1 in 3D.
* NANI    scattering anisotropy (=1 for isotropic scattering).
* NFUNL   number of moments of the flux (in 2D: NFUNL=NANI*(NANI+1)/2).
* NMOD    first dimension of ISGNR.
* KEYFLX  position of flux elements in PHI vector.
* KEYCUR  position of current elements in PHI vector.
* NZON    index-number of the mixture type assigned to each volume.
* NCONV   logical array of convergence status for each group (.TRUE.
*         not converged).
* CAZ1    first cosines of the different tracking azimuthal angles.
* CAZ2    second cosines of the different tracking azimuthal angles.
* XMU     cosines of the different tracking polar angles.
*         (polar quadrature in 2D / tracking angles in 3D).
* WZMU    polar quadrature set in 2D.
* S       total source vector components.
* SIGAL   total cross-section and albedo array.
* ISGNR   spherical harmonic signs.
* NMAX    maximum number of segments for the 3D tracks.
* NZP     number of z-plans.
* N2SOU   number of external surfaces in the 2D tracking.
* N2REG   number of regions in the 2D tracking.
* DELU    input track spacing for 3D track reconstruction.
* INDREG  region/surface index to go from the 2D to the 3D geometry.
* Z       z-plan coordinates.
* VNORF   normalization factors per angle.
* CMU     polar angle cosines.
* CMUI    inverse of polar angle cosines.
* SMU     polar angle sines.
* SMUI    inverse of polar angle sines.
* TMU     polar angle tangents.
* TMUI    inverse of polar angle tangents.
* SSYM    symmetry flag.
* IDIR    direction of fundamental current for TIBERE with MoC 
*         (=0,1,2,3).
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
      INTEGER NGEFF,K,KPN,M,N2MAX,NMU,NZON(K),NANI,NFUNL,NMOD,
     1 NREG,KEYFLX(NREG,NFUNL),KEYCUR(K-NREG),IFTRAK,N2BTR,NANGL,
     2 ISGNR(NMOD,NFUNL),NMAX,NZP,N2REG,N2SOU,
     3 INDREG(-N2SOU:N2REG,0:NZP+1),SSYM,IDIR
      REAL WZMU(NMU),SIGAL(-6:M,NGEFF),XMU(NMU),DELU,Z(0:NZP)
      DOUBLE PRECISION CAZ1(NANGL),CAZ2(NANGL),PHI(KPN,NGEFF),
     1 S(KPN,NGEFF),VNORF(NREG,NANGL,NMU,2),CMU(NMU),CMUI(NMU),
     2 SMU(NMU),SMUI(NMU),TMU(NMU),TMUI(NMU)
      LOGICAL NCONV(NGEFF)
      EXTERNAL SUBFFI,SUBFFA,SUBSCH
*----
*  LOCAL VARIABLES
*----
      INTEGER MODUR,MODDL,MODDR,MODUL
      PARAMETER(MODUR=1,MODDL=8,MODDR=5,MODUL=4)
      INTEGER II,I2LIN,IANG,N2SEG,NR2D,NBTR,KST,IST,ILINE,N3D,I,I1,I2,
     1 IMU,IANG0,NOMP,INDP,NOMM,INDM,NOMI,JF,JM,IND,TIN,N3DP,NSUB
      DOUBLE PRECISION W2D,Q0,Q1,CPO,CPOI,SPO,SPOI,TPO,TPOI,LTOT,DELTE,
     1 DELZE,T,Z1,Z2,TP,Z1P,W3DPO,WPO,W3D,OMEGAX,OMEGAY,OMEGAZ
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NOM2D,NOM3D
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RHARM
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: TRHAR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: H2D,H3D,T2D,B
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: STOT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NOM2D(N2MAX),H2D(N2MAX),NOM3D(NMAX),H3D(NMAX),B(4*NMAX),
     1 T2D(0:N2MAX-1))
*
      IF(NANI.EQ.1) THEN
*---
*  ISOTROPIC SCATTERING
*---
      DO I2LIN=1,N2BTR
         READ(IFTRAK) NSUB,N2SEG,W2D,IANG,(NOM2D(I),I=1,N2SEG),
     1                                    (H2D(I),I=1,N2SEG)
         IF(NSUB.NE.1) CALL XABORT('MCGPTF: NSUB.NE.1.')
         NR2D=N2SEG-2
         T2D(0)=0.0
         DO II=1,NR2D
            T2D(II)=T2D(II-1)+H2D(II+1)
         ENDDO
         DO IMU=1,NMU
         CPO=CMU(IMU)
         CPOI=CMUI(IMU)
         SPO=SMU(IMU)
         SPOI=SMUI(IMU)
         TPO=TMU(IMU)
         TPOI=TMUI(IMU)
         WPO=WZMU(IMU)
         IF(SSYM.EQ.2) GOTO 15
*---
*  CONSTRUCT THE 3D TRACKS WHICH ENTER THE GEOMETRY THROUGH A BOTTOM/TOP
*  SURFACE
*---
*        length of the spatial integration interval 
         LTOT=T2D(NR2D)*CPO
*        number of 3D tracks generated for this x-y track and this polar
*        direction
         NBTR=INT(LTOT/DELU)+1
*        effective track spacing in T
         DELTE=T2D(NR2D)/DBLE(NBTR)
         W3DPO=W2D*DELTE*CPO
         W3D=WPO*W3DPO
         T=-0.5D0*DELTE
         KST=1
         DO 10 ILINE=1,NBTR
            T=T+DELTE
            TP=T
            DO WHILE (T2D(KST).LT.T)
               KST=KST+1
            ENDDO
            K=KST
*           ---
*           positive polar sine track
*           ---
            I1=1
            Z1=Z(I1-1)
            TIN=0
            N3D=1
            NOM3D(N3D)=INDREG(NOM2D(K+1),0)
            H3D(N3D)=0.5D0
            CALL MCGPT1(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I1,K,Z1,
     1           T,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
            DO II=2,N3D-1
               H3D(II)=H3D(II)*VNORF(NOM3D(II),IANG,IMU,1)
            ENDDO
            IF(SSYM.EQ.1) THEN
*           the top boundary condition is a surface symmetry
            IF(TIN.EQ.0) THEN
*              this track has encountered the top boundary -> it is reflected
               N3DP=N3D
               N3D=N3D-1
               I1=I1-1
               CALL MCGPT2(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I1,K,
     1              Z1,T,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
               DO II=N3DP,N3D-1
                  H3D(II)=H3D(II)*VNORF(NOM3D(II),IANG,IMU,2)
               ENDDO
            ENDIF
            ENDIF
            DO II=1,NGEFF
            IF(NCONV(II)) THEN
*                   MCGFFI:  'Source Term Isolation' Strategy turned on
*                   MCGFFIS: 'Source Term Isolation' Strategy turned off
*                   MCGFFIT: 'MOCC/MCI' Iterative Strategy
               OMEGAX=0.0D0
               OMEGAY=0.0D0
               OMEGAZ=0.0D0
               IDIR=0
               CALL SUBFFI(SUBSCH,K,KPN,M,N3D,H3D,NOM3D,NZON,
     1              SIGAL(0,II),S(1,II),NREG,KEYFLX,KEYCUR,PHI(1,II),
     2              B,W3D,OMEGAX,OMEGAY,OMEGAZ,IDIR)
            ENDIF
            ENDDO
            T=TP
            IF(SSYM.EQ.1) GOTO 10
            K=KST
*           ---
*           negative polar sine track
*           ---
            I2=NZP
            Z2=Z(I2)
            TIN=0
            N3D=1
            NOM3D(N3D)=INDREG(NOM2D(K+1),NZP+1)
            H3D(N3D)=0.5D0
            CALL MCGPT2(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I2,K,Z2,
     1           T,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
            DO II=2,N3D-1
               H3D(II)=H3D(II)*VNORF(NOM3D(II),IANG,IMU,2)
            ENDDO
            DO II=1,NGEFF
            IF(NCONV(II)) THEN
*                   MCGFFI:  'Source Term Isolation' Strategy turned on
*                   MCGFFIS: 'Source Term Isolation' Strategy turned off
*                   MCGFFIT: 'MOCC/MCI' Iterative Strategy
               OMEGAX=0.0D0
               OMEGAY=0.0D0
               OMEGAZ=0.0D0
               IDIR=0
               CALL SUBFFI(SUBSCH,K,KPN,M,N3D,H3D,NOM3D,NZON,
     1              SIGAL(0,II),S(1,II),NREG,KEYFLX,KEYCUR,PHI(1,II),
     2              B,W3D,OMEGAX,OMEGAY,OMEGAZ,IDIR)
            ENDIF
            ENDDO
*           ---
            T=TP
 10      CONTINUE
*---
*  CONSTRUCT THE 3D TRACKS WHICH ENTER THE GEOMETRY THROUGH A LATERAL
*  SURFACE
*---
*        length of the spatial integration interval
 15      LTOT=Z(NZP)*SPO
!         LTOT=(Z(NZP)-Z(0))*SPO with Z(0)=0.0
*        number of 3D tracks generated for this x-y track and this polar
*        direction
         NBTR=INT(LTOT/DELU)+1
*        effective track spacing in Z
         DELZE=Z(NZP)/DBLE(NBTR)
!         DELZE=(Z(NZP)-Z(0))/DBLE(NBTR) with Z(0)=0.0
         W3DPO=W2D*DELZE*SPO
         W3D=WPO*W3DPO
         Z1=-0.5D0*DELZE
!         Z1=Z(0)-0.5D0*DELZE with Z(0)=0.0
         IST=1
         DO 20 ILINE=1,NBTR
            Z1=Z1+DELZE
            Z1P=Z1
            DO WHILE (Z(IST).LT.Z1)
               IST=IST+1
            ENDDO
            I=IST
*           ---
*           positive polar sine track
*           ---
            K=1
            T=T2D(K-1)
            TIN=1
            N3D=1
            N3DP=2
            NOM3D(N3D)=INDREG(NOM2D(1),IST)
            H3D(N3D)=0.5D0
 21         CONTINUE
            CALL MCGPT1(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I,K,Z1,
     1           T,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
            DO II=N3DP,N3D-1
               H3D(II)=H3D(II)*VNORF(NOM3D(II),IANG,IMU,1)
            ENDDO
            IF(SSYM.GT.0) THEN
*           the top boundary condition is a surface symmetry
            IF(TIN.EQ.0) THEN
*           this track has encountered the top boundary -> it is reflected
               N3DP=N3D
               N3D=N3D-1
               I=I-1
               CALL MCGPT2(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I,K,
     1              Z1,T,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
               DO II=N3DP,N3D-1
                  H3D(II)=H3D(II)*VNORF(NOM3D(II),IANG,IMU,2)
               ENDDO
               IF((SSYM.EQ.2).AND.(TIN.EQ.0)) THEN
*              the bottom boundary is a surface symmetry 
*              this track has encountered the bottom boundary -> it is
*              reflected
               N3DP=N3D
               N3D=N3D-1
               I=I+1
               GOTO 21
            ENDIF
            ENDIF
            ENDIF
            DO II=1,NGEFF
            IF(NCONV(II)) THEN
*                   MCGFFI:  'Source Term Isolation' Strategy turned on
*                   MCGFFIS: 'Source Term Isolation' Strategy turned off
*                   MCGFFIT: 'MOCC/MCI' Iterative Strategy
               OMEGAX=0.0D0
               OMEGAY=0.0D0
               OMEGAZ=0.0D0
               IDIR=0
               CALL SUBFFI(SUBSCH,K,KPN,M,N3D,H3D,NOM3D,NZON,
     1              SIGAL(0,II),S(1,II),NREG,KEYFLX,KEYCUR,PHI(1,II),
     2              B,W3D,OMEGAX,OMEGAY,OMEGAZ,IDIR)
            ENDIF
            ENDDO
            Z1=Z1P
            I=IST
*           ---
*           negative polar sine track
*           ---
            K=1
            T=T2D(K-1)
            TIN=1
            N3D=1
            N3DP=2
            NOM3D(N3D)=INDREG(NOM2D(1),IST)
            H3D(N3D)=0.5D0
 22         CONTINUE
            CALL MCGPT2(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I,K,Z1,
     1           T,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
            DO II=N3DP,N3D-1
               H3D(II)=H3D(II)*VNORF(NOM3D(II),IANG,IMU,2)
            ENDDO
            IF(SSYM.EQ.2) THEN
*           the bottom boundary is a surface symmetry 
            IF(TIN.EQ.0) THEN
*           this track has encountered the bottom boundary -> it is
*           reflected
               N3DP=N3D
               N3D=N3D-1
               I=I+1
               CALL MCGPT1(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I,K,
     1              Z1,T,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
               DO II=N3DP,N3D-1
                  H3D(II)=H3D(II)*VNORF(NOM3D(II),IANG,IMU,1)
               ENDDO
               IF(TIN.EQ.0) THEN
*              the top boundary is a surface symmetry 
*              this track has encountered the top boundary -> it is
*              reflected
                  N3DP=N3D
                  N3D=N3D-1
                  I=I-1
                  GOTO 22
               ENDIF
            ENDIF
            ENDIF
            DO II=1,NGEFF
            IF(NCONV(II)) THEN
*                   MCGFFI:  'Source Term Isolation' Strategy turned on
*                   MCGFFIS: 'Source Term Isolation' Strategy turned off
*                   MCGFFIT: 'MOCC/MCI' Iterative Strategy
               OMEGAX=0.0D0
               OMEGAY=0.0D0
               OMEGAZ=0.0D0
               IDIR=0
               CALL SUBFFI(SUBSCH,K,KPN,M,N3D,H3D,NOM3D,NZON,
     1              SIGAL(0,II),S(1,II),NREG,KEYFLX,KEYCUR,PHI(1,II),
     2              B,W3D,OMEGAX,OMEGAY,OMEGAZ,IDIR)
            ENDIF
            ENDDO
*           ---
            Z1=Z1P
 20      CONTINUE
         ENDDO
      ENDDO
      ELSE
*---
*  ANISOTROPIC SCATTERING
*---
      ALLOCATE(STOT(NMAX,NMU,NGEFF,2))
      CALL XDDSET(STOT,2*NMU*NMAX*NGEFF,0.0D0)
      ALLOCATE(RHARM(NMU,NFUNL),TRHAR(NMU,NFUNL,NMOD))
      IANG0=0
      DO I2LIN=1,N2BTR
         READ(IFTRAK) NSUB,N2SEG,W2D,IANG,(NOM2D(I),I=1,N2SEG),
     1                                    (H2D(I),I=1,N2SEG)
         IF(NSUB.NE.1) CALL XABORT('MCGPTF: NSUB.NE.1.')
         NR2D=N2SEG-2
         T2D(0)=0.0
         DO II=1,NR2D
            T2D(II)=T2D(II-1)+H2D(II+1)
         ENDDO
         IF(IANG.NE.IANG0) THEN
            IANG0=IANG
            CALL MOCCHR(3,NANI-1,NFUNL,NMU,XMU,CAZ1(IANG),CAZ2(IANG),
     1                  RHARM)
            DO 27 JM=1,NMOD
            DO 26 JF=1,NFUNL
            DO 25 IMU=1,NMU
               TRHAR(IMU,JF,JM)=ISGNR(JM,JF)*RHARM(IMU,JF)
 25         CONTINUE
 26         CONTINUE
 27         CONTINUE
         ENDIF
         DO IMU=1,NMU
         CPO=CMU(IMU)
         CPOI=CMUI(IMU)
         SPO=SMU(IMU)
         SPOI=SMUI(IMU)
         TPO=TMU(IMU)
         TPOI=TMUI(IMU)
         WPO=WZMU(IMU)
*---
*  CONSTRUCT THE 3D TRACKS WHICH ENTER THE GEOMETRY THROUGH A BOTTOM/TOP
*  SURFACE
*---
*        length of the spatial integration interval 
         LTOT=T2D(NR2D)*CPO
*        number of 3D tracks generated for this x-y track and this polar
*        direction
         NBTR=INT(LTOT/DELU)+1
*        effective track spacing in T
         DELTE=T2D(NR2D)/DBLE(NBTR)
         W3DPO=W2D*DELTE*CPO
         W3D=WPO*W3DPO
         T=-0.5D0*DELTE
         KST=1
         DO 30 ILINE=1,NBTR
            T=T+DELTE
            TP=T
            DO WHILE (T2D(KST).LT.T)
               KST=KST+1
            ENDDO
            K=KST
*           ---
*           positive polar sine track
*           ---
            I1=1
            Z1=Z(I1-1)
            TIN=0
            N3D=1
            NOM3D(N3D)=INDREG(NOM2D(K+1),0)
            H3D(N3D)=0.5D0
            CALL MCGPT1(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I1,K,Z1,
     1           T,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
            DO II=2,N3D-1
               H3D(II)=H3D(II)*VNORF(NOM3D(II),IANG,IMU,1)
            ENDDO
            DO II=1,NGEFF
            IF(NCONV(II)) THEN
*              incoming flux in + direction
               NOMP=NOM3D(1)
               INDP=KEYCUR(-NOMP)
*              incoming flux in - direction
               NOMM=NOM3D(N3D)
               INDM=KEYCUR(-NOMM)
               STOT(1,IMU,II,1)=W3D*S(INDP,II)
               STOT(N3D,IMU,II,2)=W3D*S(INDM,II)
*              regional sources   
               DO I=2,N3D-1
                  NOMI=NOM3D(I)
                  Q0=0.0D0
                  Q1=0.0D0
                  DO JF=1,NFUNL
                     IND=KEYFLX(NOMI,JF)         
                     Q0=Q0+S(IND,II)*TRHAR(IMU,JF,1)
                     Q1=Q1+S(IND,II)*TRHAR(IMU,JF,NMOD)
                  ENDDO                       
                  STOT(I,IMU,II,1)=W3D*Q0
                  STOT(I,IMU,II,2)=W3D*Q1
               ENDDO
*              MCGFFA:  'Source Term Isolation' Strategy turned on
*              MCGFFAS: 'Source Term Isolation' Strategy turned off
*              MCGFFAT: 'MOCC/MCI' Iterative Strategy
               CALL SUBFFA(SUBSCH,K,KPN,M,N3D,H3D,NOM3D,NZON,
     1              SIGAL(0,II),STOT(1,IMU,II,1),STOT(1,IMU,II,2),
     2              NREG,NMU,NANI,NFUNL,NMOD,TRHAR,KEYFLX,KEYCUR,IMU,
     3              PHI(1,II),B,MODUR,MODDL)
            ENDIF
            ENDDO
            T=TP
            K=KST
*           ---
*           negative polar sine track
*           ---
            I2=NZP
            Z2=Z(I2)
            TIN=0
            N3D=1
            NOM3D(N3D)=INDREG(NOM2D(K+1),NZP+1)
            H3D(N3D)=0.5D0
            CALL MCGPT2(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I2,K,Z2,
     1           T,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
            DO II=2,N3D-1
               H3D(II)=H3D(II)*VNORF(NOM3D(II),IANG,IMU,2)
            ENDDO 
            DO II=1,NGEFF
            IF(NCONV(II)) THEN
*              incoming flux in + direction
               NOMP=NOM3D(1)
               INDP=KEYCUR(-NOMP)
*              incoming flux in - direction
               NOMM=NOM3D(N3D)
               INDM=KEYCUR(-NOMM)
               STOT(1,IMU,II,1)=W3D*S(INDP,II)
               STOT(N3D,IMU,II,2)=W3D*S(INDM,II)
*              regional sources   
               DO I=2,N3D-1
                  NOMI=NOM3D(I)
                  Q0=0.0D0
                  Q1=0.0D0
                  DO JF=1,NFUNL
                     IND=KEYFLX(NOMI,JF)         
                     Q0=Q0+S(IND,II)*TRHAR(IMU,JF,1)
                     Q1=Q1+S(IND,II)*TRHAR(IMU,JF,NMOD)
                  ENDDO                       
                  STOT(I,IMU,II,1)=W3D*Q0
                  STOT(I,IMU,II,2)=W3D*Q1
               ENDDO
*              MCGFFA:  'Source Term Isolation' Strategy turned on
*              MCGFFAS: 'Source Term Isolation' Strategy turned off
*              MCGFFAT: 'MOCC/MCI' Iterative Strategy
               CALL SUBFFA(SUBSCH,K,KPN,M,N3D,H3D,NOM3D,NZON,
     1              SIGAL(0,II),STOT(1,IMU,II,1),STOT(1,IMU,II,2),
     2              NREG,NMU,NANI,NFUNL,NMOD,TRHAR,KEYFLX,KEYCUR,IMU,
     3              PHI(1,II),B,MODDR,MODUL)
            ENDIF
            ENDDO
*           ---
            T=TP
 30      CONTINUE
*---
*  CONSTRUCT THE 3D TRACKS WHICH ENTER THE GEOMETRY THROUGH A LATERAL
*  SURFACE
*---
*        length of the spatial integration interval 
         LTOT=Z(NZP)*SPO
!         LTOT=(Z(NZP)-Z(0))*SPO with Z(0)=0.0
*        number of 3D tracks generated for this x-y track and this polar
*        direction
         NBTR=INT(LTOT/DELU)+1
*        effective track spacing in Z
         DELZE=Z(NZP)/DBLE(NBTR)
!         DELZE=(Z(NZP)-Z(0))/DBLE(NBTR) with Z(0)=0.0
         W3DPO=W2D*DELZE*SPO
         W3D=WPO*W3DPO
         Z1=-0.5D0*DELZE
!         Z1=Z(0)-0.5D0*DELZE with Z(0)=0.0
         IST=1
         DO 40 ILINE=1,NBTR
            Z1=Z1+DELZE
            Z1P=Z1
            DO WHILE (Z(IST).LT.Z1)
               IST=IST+1
            ENDDO
            I=IST
*           ---
*           positive polar sine track
*           ---
            K=1
            T=T2D(K-1)
            TIN=1
            N3D=1
            NOM3D(N3D)=INDREG(NOM2D(1),IST)
            H3D(N3D)=0.5D0
            CALL MCGPT1(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I,K,Z1,
     1           T,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
            DO II=2,N3D-1
               H3D(II)=H3D(II)*VNORF(NOM3D(II),IANG,IMU,1)
            ENDDO
            DO II=1,NGEFF
            IF(NCONV(II)) THEN
*              incoming flux in + direction
               NOMP=NOM3D(1)
               INDP=KEYCUR(-NOMP)
*              incoming flux in - direction
               NOMM=NOM3D(N3D)
               INDM=KEYCUR(-NOMM)
               STOT(1,IMU,II,1)=W3D*S(INDP,II)
               STOT(N3D,IMU,II,2)=W3D*S(INDM,II)
*              regional sources   
               DO I=2,N3D-1
                  NOMI=NOM3D(I)
                  Q0=0.0D0
                  Q1=0.0D0
                  DO JF=1,NFUNL
                     IND=KEYFLX(NOMI,JF)
                     Q0=Q0+S(IND,II)*TRHAR(IMU,JF,1)
                     Q1=Q1+S(IND,II)*TRHAR(IMU,JF,NMOD)
                  ENDDO                       
                  STOT(I,IMU,II,1)=W3D*Q0
                  STOT(I,IMU,II,2)=W3D*Q1
               ENDDO
*              MCGFFA:  'Source Term Isolation' Strategy turned on
*              MCGFFAS: 'Source Term Isolation' Strategy turned off
*              MCGFFAT: 'MOCC/MCI' Iterative Strategy
               CALL SUBFFA(SUBSCH,K,KPN,M,N3D,H3D,NOM3D,NZON,
     1              SIGAL(0,II),STOT(1,IMU,II,1),STOT(1,IMU,II,2),
     2              NREG,NMU,NANI,NFUNL,NMOD,TRHAR,KEYFLX,KEYCUR,IMU,
     3              PHI(1,II),B,MODUR,MODDL)
            ENDIF
            ENDDO
            Z1=Z1P
            I=IST
*           ---
*           negative polar sine track
*           ---
            K=1
            T=T2D(K-1)
            TIN=1
            N3D=1
            NOM3D(N3D)=INDREG(NOM2D(1),IST)
            H3D(N3D)=0.5D0
            CALL MCGPT2(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I,K,Z1,
     1           T,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
            DO II=2,N3D-1
               H3D(II)=H3D(II)*VNORF(NOM3D(II),IANG,IMU,2)
            ENDDO
            DO II=1,NGEFF
            IF(NCONV(II)) THEN
*              incoming flux in + direction
               NOMP=NOM3D(1)
               INDP=KEYCUR(-NOMP)
*              incoming flux in - direction
               NOMM=NOM3D(N3D)
               INDM=KEYCUR(-NOMM)
               STOT(1,IMU,II,1)=W3D*S(INDP,II)
               STOT(N3D,IMU,II,2)=W3D*S(INDM,II)
*              regional sources   
               DO I=2,N3D-1
                  NOMI=NOM3D(I)
                  Q0=0.0D0
                  Q1=0.0D0
                  DO JF=1,NFUNL
                     IND=KEYFLX(NOMI,JF)
                     Q0=Q0+S(IND,II)*TRHAR(IMU,JF,1)
                     Q1=Q1+S(IND,II)*TRHAR(IMU,JF,NMOD)
                  ENDDO                       
                  STOT(I,IMU,II,1)=W3D*Q0
                  STOT(I,IMU,II,2)=W3D*Q1
               ENDDO
*              MCGFFA:  'Source Term Isolation' Strategy turned on
*              MCGFFAS: 'Source Term Isolation' Strategy turned off
*              MCGFFAT: 'MOCC/MCI' Iterative Strategy
               CALL SUBFFA(SUBSCH,K,KPN,M,N3D,H3D,NOM3D,NZON,
     1              SIGAL(0,II),STOT(1,IMU,II,1),STOT(1,IMU,II,2),
     2              NREG,NMU,NANI,NFUNL,NMOD,TRHAR,KEYFLX,KEYCUR,IMU,
     3              PHI(1,II),B,MODDR,MODUL)
            ENDIF
            ENDDO
*           ---
            Z1=Z1P
 40      CONTINUE
         ENDDO
      ENDDO
      DEALLOCATE(TRHAR,RHARM,STOT)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(B,H3D,NOM3D,H2D,NOM2D,T2D)
*
      RETURN
      END
