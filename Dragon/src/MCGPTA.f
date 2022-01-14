*DECK MCGPTA
      SUBROUTINE MCGPTA(NFI,NREG,NLONG,M,NANGL,NMU,LC,NGEFF,
     1                  IANGL,NSEG,NOM2D,NZONA,IPERM,KM,IM,MCU,PREV,
     2                  NEXT,W2D,ZMU,WZMU,SIGAL,XSW,T2D,DIAGQ,CQ,
     3                  DIAGF,CF,WORK,LTMT,SUBDS2,SUBDSP,SUBDSC,NR2D,
     4                  NMAX,NZP,N2REG,N2SOU,DELU,INDREG,NOM3D,NOM3D0,
     5                  H3D,H3D0,Z,VNORF,CMU,CMUI,SMU,SMUI,TMU,TMUI,
     6                  N3TR,N3TRTMT,N3SE,N3SETMT,N2TPROC,SSYM)
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
* NFI     total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NREG    number of volumes.
* NLONG   order of the corrective system.
* M       number of material mixtures.
* NANGL   number of tracking angles in the plane.
* NMU     order of the polar quadrature.
* LC      dimension of vector MCU.
* NGEFF   number of energy groups to process.
* IANGL   direction index for the current 2D track.
* NSEG    number of segments for the current 2D track.
* NOM2D   vector containing the region number of the different 
*         segments of this 2D track.
* NZONA   index-number of the mixture type assigned to each volume
*         for ACA.
* IPERM   permutation array.
* KM      used in CDD acceleration.
* IM      used in CDD acceleration.
* MCU     used in CDD acceleration.
* W2D     2D track weight.
* ZMU     polar quadrature set in 2D.
* WZMU    polar quadrature set in 2D.
* SIGAL   total cross-section and albedo array.
* XSW     scattering cross sections array. 
* T2D     vector containing the local coordinates of the segments 
*         boundaries for this 2D track.
* LTMT    track merging flag.
* SUBDS2  ACA coefficients summation subroutine.
* SUBDSP  ACA coefficients position subroutine.
* SUBDSC  ACA coefficients calculation subroutine for this 2D track.
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
* N2TPROC number of 2D tracks corresponding to this merged track (if LTMT).
* SSYM    symmetry flag.
*
*Parameters: input/output
* CQ      undefined.
* CF      undefined.
* DIAGQ   undefined.
* DIAGF   undefined.
* N3TR    total number of 3D tracks generated.
* N3TRTMT total number of 3D merged tracks.
* N3SE    total number of segments on the 3D tracks generated.
* N3SETMT total number of segments on the 3D merged tracks.
*
*Parameters: scratch
* PREV    undefined.
* NEXT    undefined.
* WORK    undefined.
* NOM3D   undefined.
* NOM3D0  undefined.
* H3D     undefined.
* H3D0    undefined.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NFI,NREG,NLONG,M,NANGL,NMU,LC,NGEFF,IANGL,NSEG,
     1 NOM2D(NSEG),NZONA(NFI),IPERM(NFI),KM(NLONG),IM(NLONG),MCU(LC),
     2 NMAX,PREV(NMAX),NEXT(NMAX),NR2D,NZP,N2REG,N2SOU,
     3 INDREG(-N2SOU:N2REG,0:NZP+1),NOM3D(NMAX),NOM3D0(NMAX,2),N3TR,
     4 N3TRTMT,N3SE,N3SETMT,N2TPROC,SSYM
      REAL ZMU(NMU),WZMU(NMU),SIGAL(-6:M,NGEFF),XSW(0:M,NGEFF),
     1 DIAGQ(NLONG,NGEFF),CQ(LC,NGEFF),DELU,Z(0:NZP)
      DOUBLE PRECISION W2D,DIAGF(NLONG,NGEFF),CF(LC,NGEFF),WORK(NMAX,3),
     1 VNORF(NREG,NANGL,NMU,2),CMU(NMU),CMUI(NMU),SMU(NMU),SMUI(NMU),
     2 TMU(NMU),TMUI(NMU),H3D(NMAX),H3D0(NMAX,2),T2D(0:NR2D)
      LOGICAL LTMT
      EXTERNAL SUBDS2,SUBDSP,SUBDSC
*----
*  LOCAL VARIABLES
*----
      INTEGER IMU,NBTR,KST,IST,ILINE,I,I1,I2,K,N3D,II,NMERG1,NMERG2,
     1 N3D01,N3D02,TIN,NTR,NSE,NTPROC,N3DP
      DOUBLE PRECISION CPO,CPOI,SPO,SPOI,TPO,TPOI,LTOT,DELTE,DELZE,T,
     1 Z1,Z2,TP,Z1P,W3D,W3D01,W3D02,W3DPO,W3DS,WPO
      LOGICAL LFORC
*
      NTR=0
      NSE=0
      DO IMU=1,NMU
*     ------------polar angle loop
      NMERG1=0
      NMERG2=0
      LFORC=.FALSE.
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
      W3DS=WPO*W3DPO
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
         W3D=W3DS
         N3D=1
         NOM3D(N3D)=INDREG(NOM2D(K+1),0)
         H3D(N3D)=0.5
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
         IF (LTMT) THEN
            NTR=NTR+1
            NSE=NSE+N3D
            CALL MCGTMT(NMERG1,N3TRTMT,N3SETMT,N3D,N3D01,NOM3D,
     1           NOM3D0(1,1),W3D,W3D01,H3D,H3D0(1,1),LFORC,NTPROC)
            IF (NTPROC.EQ.0) GOTO 31
         ENDIF
         NOM3D(1)=NREG-NOM3D(1)
         NOM3D(N3D)=NREG-NOM3D(N3D)
         DO II=1,N3D
            NOM3D(II)=IPERM(NOM3D(II))
         ENDDO
         CALL MCGDS1(SUBDS2,SUBDSP,SUBDSC,N3D,NMU,NGEFF,W3D,
     1        H3D,ZMU,WZMU,NOM3D,NZONA,NLONG,NFI,3,LC,M,KM,IM,
     2        MCU,DIAGF,DIAGQ,CF,CQ,PREV,NEXT,SIGAL,XSW,WORK)
 31      T=TP
         IF (SSYM.EQ.1) GOTO 10
         K=KST
*        ---
*        negative polar sine track
*        ---
         I2=NZP
         Z2=Z(I2)
         TIN=0
         W3D=W3DS
         N3D=1
         NOM3D(N3D)=INDREG(NOM2D(K+1),NZP+1)
         H3D(N3D)=0.5
         CALL MCGPT2(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I2,K,Z2,T,
     1        TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
         DO II=2,N3D-1
            H3D(II)=H3D(II)*VNORF(NOM3D(II),IANGL,IMU,2)
         ENDDO
         IF (LTMT) THEN
            NTR=NTR+1
            NSE=NSE+N3D
            CALL MCGTMT(NMERG2,N3TRTMT,N3SETMT,N3D,N3D02,NOM3D,
     1           NOM3D0(1,2),W3D,W3D02,H3D,H3D0(1,2),LFORC,NTPROC)
            IF (NTPROC.EQ.0) GOTO 32
         ENDIF
         NOM3D(1)=NREG-NOM3D(1)
         NOM3D(N3D)=NREG-NOM3D(N3D)
         DO II=1,N3D
            NOM3D(II)=IPERM(NOM3D(II))
         ENDDO
         CALL MCGDS1(SUBDS2,SUBDSP,SUBDSC,N3D,NMU,NGEFF,W3D,
     1        H3D,ZMU,WZMU,NOM3D,NZONA,NLONG,NFI,3,LC,M,KM,IM,
     2        MCU,DIAGF,DIAGQ,CF,CQ,PREV,NEXT,SIGAL,XSW,WORK)        
*        ---
 32      T=TP
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
      W3DS=WPO*W3DPO
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
         W3D=W3DS
         N3D=1
         N3DP=2
         NOM3D(N3D)=INDREG(NOM2D(1),IST)
         H3D(N3D)=0.5
 21      CONTINUE
         CALL MCGPT1(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I,K,Z1,T,
     1        TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
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
         IF (LTMT) THEN
            NTR=NTR+1
            NSE=NSE+N3D
            IF (ILINE.EQ.NBTR) LFORC=.TRUE.
            CALL MCGTMT(NMERG1,N3TRTMT,N3SETMT,N3D,N3D01,NOM3D,
     1           NOM3D0(1,1),W3D,W3D01,H3D,H3D0(1,1),LFORC,NTPROC)
            IF (NTPROC.EQ.0) GOTO 41
         ENDIF
         NOM3D(1)=NREG-NOM3D(1)
         NOM3D(N3D)=NREG-NOM3D(N3D)
         DO II=1,N3D
            NOM3D(II)=IPERM(NOM3D(II))
         ENDDO
         CALL MCGDS1(SUBDS2,SUBDSP,SUBDSC,N3D,NMU,NGEFF,W3D,
     1        H3D,ZMU,WZMU,NOM3D,NZONA,NLONG,NFI,3,LC,M,KM,IM,
     2        MCU,DIAGF,DIAGQ,CF,CQ,PREV,NEXT,SIGAL,XSW,WORK)
 41      Z1=Z1P
         I=IST
*        ---
*        negative polar sine track
*        ---
         K=1
         T=T2D(K-1)
         TIN=1
         W3D=W3DS
         N3D=1
         N3DP=2
         NOM3D(N3D)=INDREG(NOM2D(1),IST)
         H3D(N3D)=0.5
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
         IF (LTMT) THEN
            NTR=NTR+1
            NSE=NSE+N3D
            CALL MCGTMT(NMERG2,N3TRTMT,N3SETMT,N3D,N3D02,NOM3D,
     1           NOM3D0(1,2),W3D,W3D02,H3D,H3D0(1,2),LFORC,NTPROC)
            IF (NTPROC.EQ.0) GOTO 42
         ENDIF
         NOM3D(1)=NREG-NOM3D(1)
         NOM3D(N3D)=NREG-NOM3D(N3D)
         DO II=1,N3D
            NOM3D(II)=IPERM(NOM3D(II))
         ENDDO
         CALL MCGDS1(SUBDS2,SUBDSP,SUBDSC,N3D,NMU,NGEFF,W3D,
     1        H3D,ZMU,WZMU,NOM3D,NZONA,NLONG,NFI,3,LC,M,KM,IM,
     2        MCU,DIAGF,DIAGQ,CF,CQ,PREV,NEXT,SIGAL,XSW,WORK)
*        ---
 42      Z1=Z1P
 20   CONTINUE 
      IF (LTMT) THEN
*     process last positive polar sine track
         CALL MCGTMT(NMERG1,N3TRTMT,N3SETMT,N3D,N3D01,NOM3D,
     1        NOM3D0(1,1),W3D,W3D01,H3D,H3D0(1,1),LFORC,NTPROC)
         NOM3D(1)=NREG-NOM3D(1)
         NOM3D(N3D)=NREG-NOM3D(N3D)
         DO II=1,N3D
            NOM3D(II)=IPERM(NOM3D(II))
         ENDDO
         CALL MCGDS1(SUBDS2,SUBDSP,SUBDSC,N3D,NMU,NGEFF,W3D,
     1        H3D,ZMU,WZMU,NOM3D,NZONA,NLONG,NFI,3,LC,M,KM,IM,
     2        MCU,DIAGF,DIAGQ,CF,CQ,PREV,NEXT,SIGAL,XSW,WORK)
*     process last negative polar sine track
         CALL MCGTMT(NMERG2,N3TRTMT,N3SETMT,N3D,N3D02,NOM3D,
     1        NOM3D0(1,2),W3D,W3D02,H3D,H3D0(1,2),LFORC,NTPROC)
         NOM3D(1)=NREG-NOM3D(1)
         NOM3D(N3D)=NREG-NOM3D(N3D)
         DO II=1,N3D
            NOM3D(II)=IPERM(NOM3D(II))
         ENDDO      
         CALL MCGDS1(SUBDS2,SUBDSP,SUBDSC,N3D,NMU,NGEFF,W3D,
     1        H3D,ZMU,WZMU,NOM3D,NZONA,NLONG,NFI,3,LC,M,KM,IM,
     2        MCU,DIAGF,DIAGQ,CF,CQ,PREV,NEXT,SIGAL,XSW,WORK)
      ENDIF 
*     ------------polar angle loop   
      ENDDO
*
      N3TR=N3TR+N2TPROC*NTR
      N3SE=N3SE+N2TPROC*NSE     
      RETURN
*
      END
