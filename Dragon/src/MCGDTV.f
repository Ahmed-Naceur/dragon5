*DECK MCGDTV
      SUBROUTINE MCGDTV(NDIM,NFI,NREG,NSOU,NSEG,NMU,LMCU,LMXMCU,NZONA,
     1                  NRSEG,MCUW,MCUI,WEI2D,SEGLEN,WZMU,SURFD,CYCLIC,
     2                  ACFLAG,ZMU,XSIXYZ,CAZ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the contribution of a track to the numerical surfaces and 
* connection matrices for an EXCELT tracking.
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
* NDIM    number of dimensions for the geometry.
* NFI     total number of volumes and surfaces.
* NREG    number of regions.
* NSOU    number of external surfaces.
* NSEG    number of segments for this track.
* NMU     number of polar angles.
* LMXMCU  maximum dimension for the connection matrix.
* NZONA   index-number of the mixture/albedo type assigned to
*         each volume/surface.
* NRSEG   vector containing the region number of the different segments
*         of this track.
* WEI2D   weight for this track.
* SEGLEN  vector containing the length of the different segments of this
*         track.
* ZMU     polar quadrature points.
* WZMU    polar quadrature weights.
* CYCLIC  cyclic tracking flag.
* ACFLAG  preconditioning techniques flag.
* CAZ     directional cosines.
*
*Parameters: input/output
* LMCU    number of elements in the connection matrix.
* MCUW    temporary connection matrix.
* MCUI    temporary connection matrix.
* SURFD   numerical surfaces.
* XSIXYZ  XSI for B1 leakage.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NDIM,NFI,NREG,NSOU,NSEG,NMU,LMCU,LMXMCU,
     1 NZONA(NFI),NRSEG(NSEG),MCUW(LMXMCU),MCUI(LMXMCU)
      REAL WZMU(NMU)
      DOUBLE PRECISION WEI2D,SEGLEN(NSEG),SURFD(NSOU)
      LOGICAL CYCLIC,ACFLAG
      REAL ZMU(NMU)
      DOUBLE PRECISION CAZ(NDIM),XSIXYZ(NSOU,3)
      DOUBLE PRECISION OMEGA2(3)
      INTEGER IDIR
*---
* LOCAL VARIABLES
*---
      INTEGER II,NOMCEL,IMU,ITEMP
      DOUBLE PRECISION WEIGHT
      IF(.NOT.CYCLIC) THEN
*     non cyclic tracking: calculate numerical surfaces
         DO II=1,NSEG,NSEG-1
            NOMCEL=-NRSEG(II)
            IF (NOMCEL.GT.0) THEN
            IF (NDIM.EQ.2) THEN
               DO IMU=1,NMU
                  OMEGA2(1)=(CAZ(1)/DBLE(ZMU(IMU)))**2
                  OMEGA2(2)=(CAZ(2)/DBLE(ZMU(IMU)))**2
                  OMEGA2(3)=1.0D0-1.0D0/DBLE(ZMU(IMU)**2)
                  WEIGHT=WEI2D*DBLE(WZMU(IMU))
                  SURFD(NOMCEL)=SURFD(NOMCEL)+WEIGHT
                  DO IDIR=1,3
                    XSIXYZ(NOMCEL,IDIR)=XSIXYZ(NOMCEL,IDIR)+
     1                                  3.0D0*OMEGA2(IDIR)*WEIGHT
                  ENDDO
               ENDDO
            ELSE
               WEIGHT=WEI2D
               SURFD(NOMCEL)=SURFD(NOMCEL)+WEIGHT
               DO IDIR=1,3
                 XSIXYZ(NOMCEL,IDIR)=XSIXYZ(NOMCEL,IDIR)+
     1                               3.0D0*CAZ(IDIR)*CAZ(IDIR)*WEIGHT
               ENDDO
            ENDIF
            ENDIF
         ENDDO
      ENDIF
*
      IF(ACFLAG) THEN
*     SCR or ACA acceleration required
         DO II=1,NSEG
            IF(NRSEG(II).LT.0) THEN
               NRSEG(II)=NREG-NRSEG(II)
            ELSE IF(NRSEG(II).EQ.0) THEN
               NRSEG(II)=NREG+1
            ENDIF
         ENDDO
         IF (CYCLIC) THEN
*        cyclic tracking: "unfold" the tracking line
*                         calculate connection matrices
            CALL MCGTRK(NFI,NZONA,NSEG,NRSEG,SEGLEN)
            CALL MOCCAL(NSEG,NRSEG,NREG,MCUW,MCUI,LMCU,LMXMCU)
            DO II=1,NSEG/2
               ITEMP=NRSEG(II)
               NRSEG(II)=NRSEG(NSEG+1-II)
               NRSEG(NSEG+1-II)=ITEMP
            ENDDO
            CALL MOCCAL(NSEG,NRSEG,NREG,MCUW,MCUI,LMCU,LMXMCU)
         ELSE
*        non-cyclic tracking: calculate connection matrices
            CALL MCGCAL(NSEG,NRSEG,NREG,MCUW,MCUI,LMCU,LMXMCU)
            DO II=1,NSEG/2
               ITEMP=NRSEG(II)
               NRSEG(II)=NRSEG(NSEG+1-II)
               NRSEG(NSEG+1-II)=ITEMP
            ENDDO
            CALL MCGCAL(NSEG,NRSEG,NREG,MCUW,MCUI,LMCU,LMXMCU)
         ENDIF
      ENDIF
*     
      RETURN
      END
