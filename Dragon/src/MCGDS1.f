*DECK MCGDS1
      SUBROUTINE MCGDS1(SUBDS2,SUBDSP,SUBDSC,N,NMU,NGEFF,WEITF,HTF,ZMU,
     1                  WZMU,NOM,NZON,NLONG,NFI,NDIM,LC,M,KM,IM,MCU,
     2                  DIAGF,DIAGQ,CF,CQ,PREV,NEXT,SIGAL,XSW,WORK)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the contributions in preconditionning matrices 
* of a 2D-track.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): I. Suslov and R. Le Tellier
*
*Parameters: input
* SUBDS2  ACA coefficients summation subroutine.
* SUBDSP  ACA coefficients position subroutine.
* SUBDSC  ACA coefficients calculation subroutine.
* N       number of elements in the current track.
* NMU     order of the polar quadrature set.
* NGEFF   number of energy groups to process.
* NFI     total number of volumes and surfaces.
* NDIM    number of dimensions in the geometry.
* NLONG   total number of cells with unknowns quantities.
* M       number of material mixtures.
* LC      dimension of vector MCU.
* NZON    index-number of the mixture type assigned to each volume.
* WEITF   track weight.
* NOM     integer tracking elements.
* HTF     real tracking elements.
* ZMU     polar quadrature set.
* WZMU    polar quadrature set.
* KM      used in CDD acceleration.
* IM      used in CDD acceleration.
* MCU     used in CDD acceleration.
* SIGAL   albedos and total cross sections array.
* XSW     scattering cross sections array.
*
*Parameters: input/output
* CQ      undefined.
* CF      undefined.
* DIAGQ   undefined.
* DIAGF   undefined.
*
*Parameters: scratch
* PREV    undefined.
* NEXT    undefined.
* WORK    undefined.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NLONG,NFI,NDIM,LC,NGEFF,M,N,NMU,NOM(N),NZON(NFI),
     1 KM(NLONG),IM(NLONG),MCU(LC),PREV(N),NEXT(N)
      DOUBLE PRECISION WEITF,HTF(N)
      REAL ZMU(NMU),WZMU(NMU),DIAGQ(NLONG,NGEFF),CQ(LC,NGEFF),
     1 SIGAL(-6:M,NGEFF),XSW(0:M,NGEFF)
      DOUBLE PRECISION DIAGF(NLONG,NGEFF),CF(LC,NGEFF),WORK(N,3)
      EXTERNAL SUBDS2,SUBDSP,SUBDSC
*----
*  LOCAL VARIABLES
*----
      INTEGER IMU,I,II
      REAL ZMUI
      DOUBLE PRECISION W
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: HG
*----
*  POSITION OF COEFFICIENTS FOR THIS TRACK IN ACA MATRICES
*----
*          MCGDSP: non cyclic tracking
*          MOCDSP: cyclic tracking
      CALL SUBDSP(N,NFI,NLONG,LC,NZON,NOM,KM,MCU,IM,PREV,NEXT,HTF)
*----
*  CALCULATION OF COEFFICIENTS
*----
      IF (NDIM.EQ.3) THEN
*     3D calculation -> no loop over a polar angle   
         DO II=1,NGEFF
*                MCGDS2: non cyclic tracking
*                MOCDS2: cyclic tracking
            CALL SUBDS2(SUBDSC,LC,M,N,HTF,NOM,NZON,SIGAL(0,II),
     1           XSW(0,II),WEITF,NFI,DIAGF(1,II),DIAGQ(1,II),
     2           CF(1,II),CQ(1,II),PREV,NEXT,WORK(1,1),WORK(1,2),
     3           WORK(1,3))
         ENDDO
      ELSE
*     2D calculation -> loop over the polar angle
         ALLOCATE(HG(N))
         DO IMU=1,NMU
            ZMUI=ZMU(IMU)
            W=WEITF*WZMU(IMU)
            DO I=1,N
               IF(NZON(NOM(I)).GE.0) THEN
                  HG(I)=HTF(I)*ZMUI
               ENDIF  
            ENDDO             
            DO II=1,NGEFF
               CALL SUBDS2(SUBDSC,LC,M,N,HG,NOM,NZON,SIGAL(0,II),
     1              XSW(0,II),W,NFI,DIAGF(1,II),DIAGQ(1,II),CF(1,II),
     2              CQ(1,II),PREV,NEXT,WORK(1,1),WORK(1,2),WORK(1,3))
            ENDDO
         ENDDO
         DEALLOCATE(HG)
      ENDIF
*
      RETURN
      END
