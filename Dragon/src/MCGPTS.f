*DECK MCGPTS
      SUBROUTINE MCGPTS(SUBPJJ,NFI,NREG,M,NANI,NFUNL,NANGL,NMU,NMOD,LPS,
     1                  NPJJM,NGEFF,IANGL,IANGL0,NSEG,ISGNR,NZON,NOM2D,
     2                  IS,JS,PJJIND,W2D,XMU,CAZ1,CAZ2,RHARM,ZMU,WZMU,
     3                  TRHAR,SIGAL,T2D,PSJ,PJJD,LPJJAN,NR2D,NMAX,NZP,
     4                  N2REG,N2SOU,DELU,INDREG,NOM3D,H3D,Z,VNORF,CMU,
     5                  CMUI,SMU,SMUI,TMU,TMUI,SSYM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the PJJ and PSJ (3D prismatic extended tracking).
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
* SUBPJJ  PJJ calculation subroutine.
* NFI     total number of volumes and surfaces for which specific values
*         of the neutron flux and reactions rates are required.
* NREG    number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* M       number of material mixtures.
* NANI    number of Legendre orders.
* NFUNL   number of moments of the flux (in 2D: NFUNL=NANI*(NANI+1)/2).
* NANGL   number of tracking angles in the plane.
* NMU     order of the polar quadrature in 2D / 1 in 3D.
* NMOD    first dimension of ISGNR.
* LPS     dimension of JS.
* NPJJM   number of pjj modes to store for STIS option.
* NGEFF   number of groups to process.
* IANGL   direction index of this track.
* IANGL0  direction index of the previous track.
* NSEG    number of elements in the current track.
* ISGNR   spherical harmonic signs.
* NZON    index-number of the mixture type assigned to each volume.
* NOM2D   vector containing the region number of the different segments
*         of this 2D track.
* IS      arrays for surfaces neighbors.
* JS      JS(IS(ISOUT)+1:IS(ISOUT+1)) give the neighboring regions to
*         surface ISOUT.
* PJJIND  index of the modes for LPJJAN option.
* W2D     track weight.
* XMU     polar angle cosines.
* CAZ1    first cosines of the different tracking azimuthal angles.
* CAZ2    second cosines of the different tracking azimuthal angles.
* ZMU     polar quadrature set.
* WZMU    polar quadrature set.
* SIGAL   albedos and total cross sections array.
* T2D     vector containing the local coordinates of the segments
*         boundaries for this 2D track.
* LPJJAN  anisotropic scattering flag.
* NR2D    number of segments corresponding to regions for this 2D track.
* NMAX    maximum number of segments for the 3D tracks.
* NZP     number of z-planes.
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
*
*Parameters: input/output
* PJJD    collision probabilities.
* PSJ     leakage probabilities.
*
*Parameters: scratch
* RHARM   undefined.
* TRHAR   undefined.
* NOM3D   undefined.
* H3D     undefined.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NFI,NREG,M,NANI,NFUNL,NANGL,NMU,NMOD,LPS,NPJJM,NGEFF,
     1 IANGL,IANGL0,NSEG,ISGNR(NMOD,NFUNL),NZON(NFI),NOM2D(NSEG),
     2 IS(NFI-NREG+1),JS(LPS),PJJIND(NPJJM,2),NR2D,NMAX,NZP,N2REG,
     3 N2SOU,INDREG(-N2SOU:N2REG,0:NZP+1),NOM3D(NMAX),SSYM
      REAL XMU(NMU),RHARM(NMU,NFUNL),ZMU(NMU),WZMU(NMU),
     1 TRHAR(NMU,NFUNL,NMOD),SIGAL(-6:M,NGEFF),PSJ(LPS,NGEFF),
     2 DELU,Z(0:NZP)
      DOUBLE PRECISION PJJD(NREG,NPJJM,NGEFF),VNORF(NREG,NANGL,NMU,2),
     1 CMU(NMU),CMUI(NMU),SMU(NMU),SMUI(NMU),TMU(NMU),TMUI(NMU),
     2 W2D,T2D(0:NR2D),H3D(NMAX),CAZ1(NANGL),CAZ2(NANGL)
      LOGICAL LPJJAN
      EXTERNAL SUBPJJ
*----
*  LOCAL VARIABLES
*----
      INTEGER MODUR,MODDR
      PARAMETER(MODUR=1,MODDR=5)
      INTEGER JF,IE,IMU,NBTR,KST,IST,ILINE,I,I1,I2,K,N3D,II,TIN,N3DP
      DOUBLE PRECISION CPO,CPOI,SPO,SPOI,TPO,TPOI,LTOT,DELTE,DELZE,T,
     1 Z1,Z2,TP,Z1P,WPO,W3D,W3DPO
*
      IF (LPJJAN) THEN
         IF (IANGL.NE.IANGL0) THEN
            IF(MODDR.GT.NMOD) CALL XABORT('MCGPTS: NMOD OVERFLOW')
            IANGL0=IANGL
            CALL MOCCHR(3,NANI-1,NFUNL,NMU,XMU,CAZ1(IANGL),CAZ2(IANGL),
     1                  RHARM)
            DO 6 JF=1,NFUNL
            DO 5 IE=1,NMU
               TRHAR(IE,JF,MODUR)=ISGNR(MODUR,JF)*RHARM(IE,JF)
               TRHAR(IE,JF,MODDR)=ISGNR(MODDR,JF)*RHARM(IE,JF)
  5         CONTINUE
  6         CONTINUE
         ENDIF
      ENDIF
*
      DO IMU=1,NMU
      CPO=CMU(IMU)
      CPOI=CMUI(IMU)
      SPO=SMU(IMU)
      SPOI=SMUI(IMU)
      TPO=TMU(IMU)
      TPOI=TMUI(IMU)
      WPO=WZMU(IMU)
      IF (SSYM.EQ.2) GOTO 15
*---
* CONSTRUCT THE 3D TRACKS WHICH ENTER THE GEOMETRY THROUGH A BOTTOM/TOP SURFACE
*---
*     length of the spatial integration interval 
      LTOT=T2D(NR2D)*CPO
*     number of 3D tracks generated for this x-y track and this polar direction
      NBTR=INT(LTOT/DELU)+1
*     effective track spacing in T
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
*        ---
*        positive polar sine track
*        ---
         I1=1
         Z1=Z(I1-1)
         TIN=0
         N3D=1
         NOM3D(N3D)=NREG-INDREG(NOM2D(K+1),0)
         H3D(N3D)=0.5D0
         CALL MCGPT1(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I1,K,Z1,T,
     1        TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
         DO II=2,N3D-1
            H3D(II)=H3D(II)*VNORF(NOM3D(II),IANGL,IMU,1)
         ENDDO
         IF (SSYM.EQ.1) THEN
*        the top boundary condition is a surface symmetry
         IF (TIN.EQ.0) THEN
*        this track has encountered the top boundary -> it is reflected
            N3DP=N3D
            N3D=N3D-1
            I1=I1-1
            CALL MCGPT2(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I1,K,Z1,
     1           T,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
            DO II=N3DP,N3D-1
               H3D(II)=H3D(II)*VNORF(NOM3D(II),IANGL,IMU,2)
            ENDDO
         ENDIF
         ENDIF
         NOM3D(N3D)=NREG-NOM3D(N3D)
         CALL MCGDS4(SUBPJJ,N3D,NMU,LPS,NFUNL,NMOD,NGEFF,W3D,TRHAR,
     1        H3D,ZMU,WZMU,NOM3D,NZON,NFI,NREG,3,M,IS,JS,PJJD,PSJ,
     2        LPJJAN,NPJJM,PJJIND,SIGAL,IMU,MODUR)
         T=TP
         IF (SSYM.EQ.1) GOTO 10
         K=KST
*        ---
*        negative polar sine track
*        ---
         I2=NZP
         Z2=Z(I2)
         TIN=0
         N3D=1
         NOM3D(N3D)=NREG-INDREG(NOM2D(K+1),NZP+1)
         H3D(N3D)=0.5D0
         CALL MCGPT2(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I2,K,Z2,T,
     1        TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
         NOM3D(N3D)=NREG-NOM3D(N3D)
         DO II=2,N3D-1
            H3D(II)=H3D(II)*VNORF(NOM3D(II),IANGL,IMU,2)
         ENDDO
         CALL MCGDS4(SUBPJJ,N3D,NMU,LPS,NFUNL,NMOD,NGEFF,W3D,TRHAR,
     1        H3D,ZMU,WZMU,NOM3D,NZON,NFI,NREG,3,M,IS,JS,PJJD,PSJ,
     2        LPJJAN,NPJJM,PJJIND,SIGAL,IMU,MODDR)
*        ---
         T=TP
 10   CONTINUE
*---
* CONSTRUCT THE 3D TRACKS WHICH ENTER THE GEOMETRY THROUGH A LATERAL SURFACE
*---
*     length of the spatial integration interval 
 15   LTOT=Z(NZP)*SPO
!      LTOT=(Z(NZP)-Z(0))*SPO with Z(0)=0.0
*     number of 3D tracks generated for this x-y track and this polar direction
      NBTR=INT(LTOT/DELU)+1
*     effective track spacing in Z
      DELZE=Z(NZP)/DBLE(NBTR)
!      DELZE=(Z(NZP)-Z(0))/DBLE(NBTR) with Z(0)=0.0
      W3DPO=W2D*DELZE*SPO
      W3D=WPO*W3DPO
      Z1=-0.5D0*DELZE
!      Z1=Z(0)-0.5D0*DELZE with Z(0)=0.0
      IST=1
      DO 20 ILINE=1,NBTR
         Z1=Z1+DELZE
         Z1P=Z1
         DO WHILE (Z(IST).LT.Z1)
            IST=IST+1
         ENDDO
         I=IST
*        ---
*        positive polar sine track
*        ---
         K=1
         T=T2D(K-1)
         TIN=1
         N3D=1
         N3DP=2
         NOM3D(N3D)=NREG-INDREG(NOM2D(1),IST)
         H3D(N3D)=0.5D0
 21      CONTINUE
         CALL MCGPT1(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I,K,Z1,T,
     1      TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
         DO II=N3DP,N3D-1
            H3D(II)=H3D(II)*VNORF(NOM3D(II),IANGL,IMU,1)
         ENDDO
         IF (SSYM.GT.0) THEN
*        the top boundary condition is a surface symmetry
         IF (TIN.EQ.0) THEN
*        this track has encountered the top boundary -> it is reflected
            N3DP=N3D
            N3D=N3D-1
            I=I-1
            CALL MCGPT2(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I,K,Z1,
     1           T,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
            DO II=N3DP,N3D-1
               H3D(II)=H3D(II)*VNORF(NOM3D(II),IANGL,IMU,2)
            ENDDO
            IF ((SSYM.EQ.2).AND.(TIN.EQ.0)) THEN
*           the bottom boundary is a surface symmetry 
*           this track has encountered the bottom boundary -> it is reflected
               N3DP=N3D
               N3D=N3D-1
               I=I+1
               GOTO 21
            ENDIF
         ENDIF
         ENDIF
         NOM3D(N3D)=NREG-NOM3D(N3D)
         CALL MCGDS4(SUBPJJ,N3D,NMU,LPS,NFUNL,NMOD,NGEFF,W3D,TRHAR,
     1        H3D,ZMU,WZMU,NOM3D,NZON,NFI,NREG,3,M,IS,JS,PJJD,PSJ,
     2        LPJJAN,NPJJM,PJJIND,SIGAL,IMU,MODUR)
         Z1=Z1P
         I=IST
*        ---
*        negative polar sine track
*        ---
         K=1
         T=T2D(K-1)
         TIN=1
         N3D=1
         N3DP=2
         NOM3D(N3D)=NREG-INDREG(NOM2D(1),IST)
         H3D(N3D)=0.5D0
 22      CONTINUE
         CALL MCGPT2(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I,K,Z1,T,
     1        TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
         DO II=N3DP,N3D-1
            H3D(II)=H3D(II)*VNORF(NOM3D(II),IANGL,IMU,2)
         ENDDO
         IF (SSYM.EQ.2) THEN
*        the bottom boundary is a surface symmetry 
         IF (TIN.EQ.0) THEN
*        this track has encountered the bottom boundary -> it is reflected
            N3DP=N3D
            N3D=N3D-1
            I=I+1
            CALL MCGPT1(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I,K,Z1,
     1           T,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
            DO II=N3DP,N3D-1
               H3D(II)=H3D(II)*VNORF(NOM3D(II),IANGL,IMU,1)
            ENDDO
            IF (TIN.EQ.0) THEN
*           the top boundary is a surface symmetry 
*           this track has encountered the top boundary -> it is reflected
               N3DP=N3D
               N3D=N3D-1
               I=I-1
               GOTO 22
            ENDIF
         ENDIF
         ENDIF
         NOM3D(N3D)=NREG-NOM3D(N3D)
         CALL MCGDS4(SUBPJJ,N3D,NMU,LPS,NFUNL,NMOD,NGEFF,W3D,TRHAR,
     1        H3D,ZMU,WZMU,NOM3D,NZON,NFI,NREG,3,M,IS,JS,PJJD,PSJ,
     2        LPJJAN,NPJJM,PJJIND,SIGAL,IMU,MODDR)
*        ---
         Z1=Z1P
 20   CONTINUE
*
      ENDDO
      RETURN
      END
