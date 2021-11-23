*DECK DETINT
      SUBROUTINE DETINT(NX,NY,NZ,NEL,NUN,LPARAB,MESHX,MESHY,MESHZ,
     +           KEYF,FLUX,NGRP,DEVPOS,RESP,IPRT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for the interpolation.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): 
* E. Varin, M. Guyot
*
*Parameters:
* NX     number of x mesh-splitted elements 
* NY     number of y mesh-splitted elements 
* NZ     number of z mesh-splitted elements
* NEL    number of finite elements
* NUN    number of unknowns
* LPARAB =.TRUE. if parabolic interpolation is performed
* MESHX  regions coordinates according to x
* MESHY  regions coordinates according to y
* MESHZ  regions coordinates according to z
* KEYF   keyflux recover from L_TRACk object
* FLUX   flux for each mesh-splitted elements
* NGRP   number of energy groups
* DEVPOS detector coordinates
* RESP   flux reads by the detector
* IPRT   printing index
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NX,NY,NZ,NEL,NUN,NGRP,IPRT,KEYF(NEL)
      REAL MESHX(NX+1),MESHY(NY+1),MESHZ(NZ+1),FLUX(NUN,NGRP),RESP,
     1 DEVPOS(6)
      LOGICAL LPARAB
*----
*  LOCAL VARIABLES
*----
      INTEGER NXP1,NYP1,NZP1,NDET,I,IM
      REAL COR(3)
      REAL, ALLOCATABLE, DIMENSION(:) :: XCT,YCT,ZCT
*----
*  SCRATCH STORAGE ALLOCATION
*   XCT    center coordinates of each mesh-splitted elements for x
*   YCT    center coordinates of each mesh-splitted elements for y
*   ZCT    center coordinates of each mesh-splitted elements for z
*   COR    center detector coordinates
*----
      ALLOCATE(XCT(NX),YCT(NY),ZCT(NZ))
*
      NXP1 = NX+1
      NYP1 = NY+1
      NZP1 = NZ+1

      IF(IPRT.GT.1) 
     +  WRITE(6,*) 'INTERPOLATION POLYNOMIALE DES LECTURES AUX VANADIUM'
      NDET = 1
*----
*  CENTER MESH CALCULATION
*----
      DO 10 I=1,NX
        XCT(I) = (MESHX(I+1) + MESHX(I)) /2.
   10 CONTINUE
      DO 11 I=1,NY
        YCT(I) = (MESHY(I+1) + MESHY(I)) /2.
   11 CONTINUE
      DO 12 I=1,NZ
        ZCT(I) = (MESHZ(I+1) + MESHZ(I)) /2.
   12 CONTINUE
*----
*  CENTER DETECTOR COORDINATE
*----
      DO 13 I=1,3
        COR(I) = (DEVPOS(2*I) + DEVPOS(2*I-1)) /2.
   13 CONTINUE
      IF(LPARAB) THEN
*----
*  POLYNOMIAL FLUX INTERPOLATION AT DETECTOR SITES
*----
         CALL DETCTL(NX,NY,NZ,NEL,FLUX(1,2),RESP,NDET,XCT,YCT,ZCT,COR,
     >   KEYF,IPRT)
      ELSE
         IM = MAX(NX,NY)
         IM = MAX(IM,NZ)
         CALL DETSPL(NX,NY,NZ,IM,FLUX(1,2),RESP,NDET,XCT,YCT,ZCT,COR,
     >   KEYF,IPRT)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ZCT,YCT,XCT)
      RETURN
      END
