*DECK MCGPT2
      SUBROUTINE MCGPT2(N2SOU,N2REG,NZP,NR2D,INDREG,Z,NOM2D,T2D,I,K,ZI,
     1                  TK,TIN,CPOI,SPOI,TPO,TPOI,N3D,NOM3D,H3D)
*     
*-----------------------------------------------------------------------
*
*Purpose:
* Reconstruct 3D track for a 3D prismatic geometry from a 2D track.
* polar angle in [pi/2, pi] case.
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
* NR2D    number of segments corresponding to regions for this 2D track.
* INDREG  region/surface index to go from the 2D to the 3D geometry.
* Z       z-plan coordinates.
* NOM2D   vector containing the region number of the different segments of
*         this 2D track.
* T2D     vector containing the local coordinates of the segments
*         boundaries for this 2D track.
* CPOI    inverse of the polar cosine.
* SPOI    inverse of the polar sine.
* TPO     polar tangent.
* TPOI    polar cotangent.
*
*Parameters: input/output
* I       starting/ending z plan.
* K       starting/ending x-y tracking segment.
* ZI      strating/ending z coordinate.
* TK      starting/ending x-y tracking coordinate.
* TIN     orientation of the starting/surface.
*
*Parameters: output
* N3D     number of segments for this 3D track.
* NOM3D   vector containing the region number of the different segments
*         of this 3D track.
* H3D     vector containing the length of the different segments of this
*         3D track.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N2SOU,N2REG,NZP,NR2D,INDREG(-N2SOU:N2REG,0:NZP+1),
     1 NOM2D(NR2D+2),I,K,N3D,NOM3D(*),TIN
      REAL Z(0:NZP)
      DOUBLE PRECISION ZI,TK,CPOI,SPOI,TPO,TPOI,H3D(N3D),T2D(0:NR2D)
*---
* LOCAL VARIABLES
*---
!!      INTEGER N3DP
      DOUBLE PRECISION DELZ,DELT,X,L
!!      REAL Tstart,Zstart
*
!!!!      Tstart=TK
!!!!      Zstart=ZI
!!!!      N3DP=N3D
      DO WHILE ((I.GE.1).AND.(K.LE.NR2D))
         N3D=N3D+1
         NOM3D(N3D)=INDREG(NOM2D(K+1),I)
         IF (TIN.EQ.0) THEN
*        track enters the region through top boundary
            DELZ=Z(I)-Z(I-1)
            DELT=T2D(K)-TK
         ELSE
*        track enters the region through right boundary
            DELZ=ZI-Z(I-1)
            DELT=T2D(K)-T2D(K-1)
         ENDIF
         X=DELT/DELZ
         IF (TPO.LT.X) THEN
*        track leaves the region through the bottom boundary
            ZI=Z(I-1)
            I=I-1
            TK=TK+DELZ*TPO
            L=DELZ*CPOI
            TIN=0
         ELSE
*        track leaves the region through the right boundary
            TK=T2D(K)
            K=K+1
            ZI=ZI-DELT*TPOI
            L=DELT*SPOI
            TIN=1
         ENDIF
         IF (L.GT.0.0) THEN
            H3D(N3D)=DBLE(L)
         ELSE
            N3D=N3D-1
         ENDIF
      ENDDO
      N3D=N3D+1
      H3D(N3D)=0.5
      NOM3D(N3D)=INDREG(NOM2D(K+1),I)
!!!!      call MCGPTP(8,N3D-N3DP+1,NOM3D(N3DP),Tstart,ZStart,CPO,SPO,
!!!!     1     H3D(N3DP),2)
*
      RETURN
*
      END
