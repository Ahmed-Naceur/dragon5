*DECK FMAC04 
      SUBROUTINE FMAC04(NGPRT,NGP,NPART,NK,H)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Reorganize boundary information from FMAC-M array.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): C. Bienvenue
*
*Parameters: input
* NGPRT   number of energy groups per particle.
* NGP     number of energy groups for all particles.
* NPART   number of particles in FMAC-M file.
* NK      number of mixtures
* H       boundary informations array from FMAC-M file.
*
*Parameters: output
* H       boundary informations array organized by particle.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGPRT(NPART),NGP,NPART
      REAL H(NK,NGP+NPART),HTEMP(NK,NGP+NPART)
*
      NGT1=0
      NGT2=0
      DO I=1,NPART
        DO J=1,NGPRT(I)+1
          DO K=1,NK
            IF(J.LE.NGPRT(I)) THEN
              HTEMP(k,NGT2+J)=H(k,NGT1+J)
            ELSE
              HTEMP(k,NGT2+J)=H(k,NGP+NPART+1-I)
            ENDIF
          ENDDO
        ENDDO
        NGT1=NGT1+NGPRT(I)
        NGT2=NGT2+NGPRT(I)+1
      ENDDO
      
      H=HTEMP

      RETURN
      END
