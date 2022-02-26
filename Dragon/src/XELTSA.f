*DECK XELTSA
      SUBROUTINE XELTSA(NDIM,ABSC,INDC,DENS,ANGTSA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To compute the integration points and periodic density
* for cyclic tracking.
*
*Copyright:
* Copyright (C) 1994 Ecole Polytechnique de Montreal.
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* NDIM   number of dimensions for the problem.
* ABSC   multidimensional width of the cell.
* INDC   index of each coordinate of the angles.
*
*Parameters: output
* DENS   effective periodic density.
* ANGTSA tracking direction and its normal.
*
*Reference:
* R. Roy, G. Marleau, A. Hebert and D. Rozon,
* A Cyclic Tracking Procedure for Collision Probability Calculations
* in 2-D Lattice', Advances in Mathematics, Computations and 
* Reactor Physics, Pittsburgh, PA, April 28 - May 2 (1991).
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER          NDIM
      DOUBLE PRECISION ABSC(2)
      INTEGER          INDC(NDIM)
      DOUBLE PRECISION DENS
      DOUBLE PRECISION ANGTSA(NDIM,NDIM)
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='XELTSA')
*----
*  LOCAL VARIABLES
*----
      INTEGER          IDIM
      IF(NDIM .NE. 2) CALL XABORT(NAMSBR//
     >': Only 2-D problems permitted yet')
      DENS=0.0D0
      DO 10 IDIM= 1,NDIM
         ANGTSA(IDIM,1)= INDC(IDIM)/ABSC(NDIM+1-IDIM)
         DENS= DENS + ANGTSA(IDIM,1)*ANGTSA(IDIM,1)
   10 CONTINUE
      DENS= SQRT(DENS)
*----
*  ANGTSA(*,1) is the track direction 
*----
      DO 20 IDIM= 1,NDIM
         ANGTSA(IDIM,1)= ANGTSA(IDIM,1)/DENS
   20 CONTINUE
*----
*  ANGTSA(*,2) is a normal to track direction
*----
      ANGTSA(2,2)= -ANGTSA(1,1)
      ANGTSA(1,2)=  ANGTSA(2,1)
*----
*  Processing finished, return
*----
      RETURN
      END
