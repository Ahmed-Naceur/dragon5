*DECK NSS3TR
      SUBROUTINE NSS3TR(ITRIAL,NEL,NMIX,MAT,XX,DIFF,SIGR,SIGT,B11)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of fission system matrices for the nodal expansion method.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* ITRIAL  type of base (=1: polynomial; =2: hyperbolic)
* NEL     number of nodes
* NMIX    number of mixtures
* MAT     node mixtures
* XX      node widths
* DIFF    diffusion coefficients.
* SIGR    macroscopic removal cross section.
* SIGT    fission cross section.
*
*Parameters: output
* B11     assembly matrix.
*
*-----------------------------------------------------------------------
*
      INTEGER ITRIAL(NMIX),NEL,NMIX,MAT(NEL)
      REAL XX(NEL),DIFF(NMIX),SIGR(NMIX),SIGT(NMIX),B11(5*NEL,5*NEL)
*
      B11(:5*NEL,:5*NEL)=0.0
      NUM1=0
      DO KEL=1,NEL
        IBM=MAT(KEL)
        SIGG=SIGT(IBM)
        ETA=XX(KEL)*SQRT(SIGR(IBM)/DIFF(IBM))
        ! WEIGHT RESIDUAL EQUATIONS:
        B11(NUM1+1,NUM1+1)=SIGG
        B11(NUM1+2,NUM1+2)=SIGG/12.0
        B11(NUM1+3,NUM1+3)=SIGG/20.0
        IF(ITRIAL(IBM) == 1) THEN
          B11(NUM1+2,NUM1+4)=-SIGG/120.0
          B11(NUM1+3,NUM1+5)=-SIGG/700.0
        ELSE
          ALP1=ETA*COSH(ETA/2.0)-2.0*SINH(ETA/2.0)
          ALP2=((12.0+ETA**2)*SINH(ETA/2.0)-6.0*ETA*COSH(ETA/2.0))/ETA
          B11(NUM1+2,NUM1+4)=SIGG*ALP1/(ETA**2)
          B11(NUM1+3,NUM1+5)=SIGG*ALP2/(ETA**2)
        ENDIF
        NUM1=NUM1+5
      ENDDO
      RETURN
      END SUBROUTINE NSS3TR
