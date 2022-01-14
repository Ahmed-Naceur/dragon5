*DECK MCGPTV
      SUBROUTINE MCGPTV(N2SOU,N2REG,NZP,SSYM,N3REG,N3SOU,N2D,NR2D,NANGL,
     1                  NMU,LMCU,LMXMCU,IANGL,INDREG,NOM2D,MCUW,MCUI,Z,
     2                  T2D,W2D,CMU,CMUI,SMU,SMUI,TMU,TMUI,WZMU,DELU,
     3                  NOM3D,H3D,SURF,VNUM,ACFLAG)
*     
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the contribution of reconstructed tracks to the numerical 
* surfaces/volumes and connection matrices for a 3D prismatic extended
* tracking (from a 2D EXCELT tracking).
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
* N2SOU   number of external surfaces in the 2D tracking.
* N2REG   number of regions in the 2D tracking.
* NZP     number of z-planes.
* SSYM    symmetry flag (0=none, 1=top, 2=top and bottom).
* N3REG   number of regions in the 3D geometry.
* N3SOU   number of external surfaces in the 3D geometry.
* N2D     number of segments for this 2D track.
* NR2D    number of segments corresponding to regions for this 2D track.
* NANGL   number of plan tracking angles.
* NMU     number of polar angles.
* LMXMCU  maximum dimension for the connection matrix.
* IANGL   index of the tracking angle considered.
* INDREG  region/surface index to go from the 2D to the 3D geometry.
* NOM2D   vector containing the region number of the different segments 
*         of this 2D track.
* Z       z-plan coordinates.
* T2D     vector containing the local coordinates of the segments 
*         boundaries for this 2D track.
* W2D     weight for this 2D track.
* WZMU    polar quadrature weight.
* DELU    input track spacing for 3D track reconstruction.
* ACFLAG  preconditioning flag.
*
*Parameters: input/output
* LMCU    number of elements in the connection matrix.
* MCUW    temporary connection matrix.
* MCUI    temporary connection matrix.
* SURF    numerical surfaces.
* VNUM    numerical volumes.
* 
*Parameters: 
* CMU     undefined.
* CMUI    undefined.
* SMU     undefined.
* SMUI    undefined.
* TMU     undefined.
* TMUI    undefined.
*
*Parameters: scratch
* NOM3D   undefined.
* H3D     undefined.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N2SOU,N2REG,NZP,SSYM,N3REG,N3SOU,N2D,NR2D,NANGL,NMU,
     1 LMCU,LMXMCU,IANGL,INDREG(-N2SOU:N2REG,0:NZP+1),NOM2D(N2D),
     2 MCUW(LMXMCU),MCUI(LMXMCU),NOM3D(*)
      REAL Z(0:NZP),WZMU(NMU),DELU
      DOUBLE PRECISION W2D,T2D(0:NR2D),H3D(*),SURF(N3SOU),CMU(NMU),
     1 CMUI(NMU),SMU(NMU),SMUI(NMU),TMU(NMU),TMUI(NMU),
     2 VNUM(N3REG,NANGL,NMU,2)
      LOGICAL ACFLAG
*---
* LOCAL VARIABLES
*---
      INTEGER IMU,NBTR,KST,IST,ILINE,I,I1,I2,K,N3D,II,ITEMP,TIN,N3DP
      DOUBLE PRECISION CPO,CPOI,SPO,SPOI,TPO,TPOI,LTOT,DELTE,DELZE,T,Z1,
     1 Z2,TP,Z1P,WPO,W3D,W3DPO
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
*  CONSTRUCT THE 3D TRACKS WHICH ENTER THE GEOMETRY THROUGH A BOTTOM/TOP
*  SURFACE
*---
*     length of the spatial integration interval 
      LTOT=T2D(NR2D)*CPO
*     number of 3D tracks generated for this x-y track and this polar
*     direction
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
         NOM3D(N3D)=INDREG(NOM2D(K+1),0)
         H3D(N3D)=0.5D0
         CALL MCGPT1(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I1,K,Z1,T,
     1        TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
         SURF(-NOM3D(1))=SURF(-NOM3D(1))+W3D
         NOM3D(1)=N3REG-NOM3D(1)
         DO II=2,N3D-1
            VNUM(NOM3D(II),IANGL,IMU,1)=VNUM(NOM3D(II),IANGL,IMU,1)
     1                                 +H3D(II)*W3DPO
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
               VNUM(NOM3D(II),IANGL,IMU,2)=VNUM(NOM3D(II),IANGL,IMU,2)
     1                                    +H3D(II)*W3DPO
            ENDDO
         ENDIF
         ENDIF
         SURF(-NOM3D(N3D))=SURF(-NOM3D(N3D))+W3D
         NOM3D(N3D)=N3REG-NOM3D(N3D)
         IF(ACFLAG) THEN
            CALL MCGCAL(N3D,NOM3D,N3REG,MCUW,MCUI,LMCU,LMXMCU)
            DO II=1,N3D/2
               ITEMP=NOM3D(II)
               NOM3D(II)=NOM3D(N3D+1-II)
               NOM3D(N3D+1-II)=ITEMP
            ENDDO
            CALL MCGCAL(N3D,NOM3D,N3REG,MCUW,MCUI,LMCU,LMXMCU)
         ENDIF
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
         NOM3D(N3D)=INDREG(NOM2D(K+1),NZP+1)
         H3D(N3D)=0.5D0
         CALL MCGPT2(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I2,K,Z2,T,
     1        TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
         SURF(-NOM3D(1))=SURF(-NOM3D(1))+W3D
         NOM3D(1)=N3REG-NOM3D(1)
         DO II=2,N3D-1
            VNUM(NOM3D(II),IANGL,IMU,2)=VNUM(NOM3D(II),IANGL,IMU,2)
     1                               +H3D(II)*W3DPO
         ENDDO
         SURF(-NOM3D(N3D))=SURF(-NOM3D(N3D))+W3D
         NOM3D(N3D)=N3REG-NOM3D(N3D)
         IF(ACFLAG) THEN
            CALL MCGCAL(N3D,NOM3D,N3REG,MCUW,MCUI,LMCU,LMXMCU)
            DO II=1,N3D/2
               ITEMP=NOM3D(II)
               NOM3D(II)=NOM3D(N3D+1-II)
               NOM3D(N3D+1-II)=ITEMP
            ENDDO
            CALL MCGCAL(N3D,NOM3D,N3REG,MCUW,MCUI,LMCU,LMXMCU)
         ENDIF
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
         NOM3D(N3D)=INDREG(NOM2D(1),IST)
         H3D(N3D)=0.5D0
         SURF(-NOM3D(N3D))=SURF(-NOM3D(N3D))+W3D
         NOM3D(N3D)=N3REG-NOM3D(N3D)
 21      CONTINUE
         CALL MCGPT1(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I,K,Z1,T,
     1        TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
         DO II=N3DP,N3D-1
            VNUM(NOM3D(II),IANGL,IMU,1)=VNUM(NOM3D(II),IANGL,IMU,1)
     1                               +H3D(II)*W3DPO
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
               VNUM(NOM3D(II),IANGL,IMU,2)=VNUM(NOM3D(II),IANGL,IMU,2)
     1                                    +H3D(II)*W3DPO
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
         SURF(-NOM3D(N3D))=SURF(-NOM3D(N3D))+W3D
         NOM3D(N3D)=N3REG-NOM3D(N3D)
         IF(ACFLAG) THEN
            CALL MCGCAL(N3D,NOM3D,N3REG,MCUW,MCUI,LMCU,LMXMCU)
            DO II=1,N3D/2
               ITEMP=NOM3D(II)
               NOM3D(II)=NOM3D(N3D+1-II)
               NOM3D(N3D+1-II)=ITEMP
            ENDDO
            CALL MCGCAL(N3D,NOM3D,N3REG,MCUW,MCUI,LMCU,LMXMCU)
         ENDIF
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
         NOM3D(N3D)=INDREG(NOM2D(1),IST)
         H3D(N3D)=0.5D0
         SURF(-NOM3D(N3D))=SURF(-NOM3D(N3D))+W3D
         NOM3D(N3D)=N3REG-NOM3D(N3D)
 22      CONTINUE
         CALL MCGPT2(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I,K,Z1,T,
     1        TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
         DO II=N3DP,N3D-1
            VNUM(NOM3D(II),IANGL,IMU,2)=VNUM(NOM3D(II),IANGL,IMU,2)
     1                               +H3D(II)*W3DPO
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
               VNUM(NOM3D(II),IANGL,IMU,1)=VNUM(NOM3D(II),IANGL,IMU,1)
     1                                    +H3D(II)*W3DPO
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
         SURF(-NOM3D(N3D))=SURF(-NOM3D(N3D))+W3D
         NOM3D(N3D)=N3REG-NOM3D(N3D)
         IF(ACFLAG) THEN
            CALL MCGCAL(N3D,NOM3D,N3REG,MCUW,MCUI,LMCU,LMXMCU)
            DO II=1,N3D/2
               ITEMP=NOM3D(II)
               NOM3D(II)=NOM3D(N3D+1-II)
               NOM3D(N3D+1-II)=ITEMP
            ENDDO
            CALL MCGCAL(N3D,NOM3D,N3REG,MCUW,MCUI,LMCU,LMXMCU)
         ENDIF
*        ---
         Z1=Z1P
 20   CONTINUE
      ENDDO
*
      RETURN
      END
