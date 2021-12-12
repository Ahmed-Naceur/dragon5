*DECK DUO005
      SUBROUTINE DUO005(NGRP,DLK,ALK,FLUX,AFLUX,V,W)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute a consistent direct-adjoint leakage multigroup term using
* the Villarino-Stamm'ler normalization method.
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NGRP    number of energy groups.
* DLK     leakage rates for the direct balance equation.
* ALK     leakage rates for the adjoint balance equation.
* FLUX    integrated direct flux.
* AFLUX   integrated adjoint flux.
*
*Parameters: output
* V       first Villarino-Stamm'ler normalization vector.
* W       second Villarino-Stamm'ler normalization vector.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGRP
      REAL DLK(NGRP),ALK(NGRP),FLUX(NGRP),AFLUX(NGRP),V(NGRP),W(NGRP)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION ZNUMD,ZNUMA,GAR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DFF,DFF2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(DFF(2*NGRP+1,2*NGRP),DFF2(2*NGRP,2*NGRP+1))
*
      ZNUMD=0.0D0
      ZNUMA=0.0D0
      DO IGR=1,NGRP
        ZNUMD=ZNUMD+FLUX(IGR)
        ZNUMA=ZNUMA+AFLUX(IGR)
      ENDDO
      CALL XDDSET(DFF,(2*NGRP+1)*(2*NGRP),0.0D0)
      DO IGR=1,NGRP
        DFF(IGR,IGR)=ZNUMD
        DFF(NGRP+IGR,NGRP+IGR)=ZNUMA
        DO JGR=1,NGRP
          DFF(IGR,NGRP+JGR)=FLUX(JGR)
          DFF(NGRP+IGR,JGR)=AFLUX(JGR)
        ENDDO
        DFF(2*NGRP+1,IGR)=AFLUX(IGR)
        DFF(2*NGRP+1,NGRP+IGR)=FLUX(IGR)
      ENDDO
      CALL XDDSET(DFF2,(2*NGRP)*(2*NGRP+1),0.0D0)
      DO IGR=1,2*NGRP
        GAR=0.0D0
        DO JGR=1,NGRP
          GAR=GAR+DFF(JGR,IGR)*DLK(JGR)+DFF(NGRP+JGR,IGR)*ALK(JGR)
        ENDDO
        DFF2(IGR,2*NGRP+1)=GAR
        DO JGR=1,2*NGRP
          DO KGR=1,2*NGRP+1
            DFF2(IGR,JGR)=DFF2(IGR,JGR)+DFF(KGR,IGR)*DFF(KGR,JGR)
          ENDDO
        ENDDO
      ENDDO
      CALL ALSBD(2*NGRP,1,DFF2,IER,2*NGRP)
      IF(IER.NE.0) CALL XABORT('DUO005: SINGULAR MATRIX.')
      DO IGR=1,NGRP
        V(IGR)=REAL(DFF2(IGR,2*NGRP+1))
        W(IGR)=REAL(DFF2(NGRP+IGR,2*NGRP+1))
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DFF2,DFF)
      RETURN
      END
