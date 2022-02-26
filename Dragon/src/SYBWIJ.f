*DECK SYBWIJ
      SUBROUTINE SYBWIJ (NREG,MAXPTS,SIGW,PIJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Scattering reduction for collision probabilities.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NREG    total number of regions.
* MAXPTS  first dimension of matrix PIJ.
* SIGW    P0 within-group scattering macroscopic cross sections
*         ordered by volume.
*
*Parameters: input/output
*  PIJ    reduced collision probability matrix on input and
*         scattering-reduced collision probability matrix at output.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NREG,MAXPTS
      REAL SIGW(NREG),PIJ(MAXPTS,NREG)
*----
*  LOCAL VARIABLES
*----
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WIJ
*
      ALLOCATE(WIJ(NREG,2*NREG))
      DO 20 I=1,NREG
      DO 10 J=1,NREG
      WIJ(I,NREG+J)=PIJ(I,J)
      WIJ(I,J)=-PIJ(I,J)*SIGW(J)
   10 CONTINUE
      WIJ(I,I)=1.0+WIJ(I,I)
   20 CONTINUE
      CALL ALSB(NREG,NREG,WIJ,IER,NREG)
      IF(IER.NE.0) CALL XABORT('SYBWIJ: SINGULAR MATRIX.')
      DO 40 J=1,NREG
      DO 30 I=1,NREG
      PIJ(I,J)=WIJ(I,NREG+J)
   30 CONTINUE
   40 CONTINUE
      DEALLOCATE(WIJ)
      RETURN
      END
