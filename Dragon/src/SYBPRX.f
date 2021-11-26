*DECK SYBPRX
      SUBROUTINE SYBPRX (NCOUR,IPAS,IKG,SIGT,P,PIS,PSS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print the cell-wise collision probabilities in SYBRX- modules.
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
* NCOUR   total number of surfaces.
* IPAS    total number of volumes.
* IKG     generating cell indices.
* SIGT    total macroscopic cross sections.
* P       reduced collision probabilities.
* PIS     volume to surface probabilities.
* PSS     surface to surface probabilities.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NCOUR,IPAS,IKG
      REAL SIGT(IPAS),P(IPAS,IPAS),PIS(NCOUR,IPAS),PSS(NCOUR,NCOUR)
*
      WRITE (6,100) IKG
      WRITE (6,110) (SIGT(I),I=1,IPAS)
      WRITE (6,'(//16H P(I,J) MATRIX :/)')
      DO 10 I=1,IPAS
      WRITE (6,110) (P(I,J),J=1,IPAS)
10    CONTINUE
      WRITE (6,'(//16H PIS(I) MATRIX :/)')
      DO 20 I=1,IPAS
      WRITE (6,110) (PIS(J,I),J=1,NCOUR)
20    CONTINUE
      WRITE (6,'(//13H PSS MATRIX :/)')
      DO 30 I=1,NCOUR
      WRITE (6,110) (PSS(I,J),J=1,NCOUR)
30    CONTINUE
      WRITE (6,'(//)')
      RETURN
100   FORMAT (1H1//19H GENERATING CELL NB,I4//21H TOTAL MACROSCOPIC CR,
     $ 14HOSS SECTIONS :/)
110   FORMAT (1X,1P,10E13.5)
      END
