*DECK MOCDDFT
      SUBROUTINE MOCDDFT(N,NREG,NSOUT,M,NOM,NZON,H,SIGANG,EXPT,EXP2,
     1           NMU,ZMU)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculate coefficients of a track for the cyclic characteristics 
* integration: Diamond-Differencing scheme without fix-up and
* 'MOCC/MCI' integration strategy.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy and R. Le Tellier
*
*Parameters: input
* N       number of elements in the current track.
* NREG    number of volumes.
* NSOUT   number of surfaces.
* M       number of material mixtures.
* NOM     vector containing the region number of the different segments
*         of this track.
* NZON    index-number of the mixture type assigned to each volume.
* H       vector containing the lenght of the different segments of this
*         track.
* SIGANG  macroscopic total cross sections and albedos.
* NMU     order of the polar quadrature set.
* ZMU     inverse of polar quadrature cosines.
*
*Parameters: output
* EXPT    track coefficient.
* EXP2    quadratic expansion of (1-exp(-a*L))/L with small argument.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N,NREG,NSOUT,M,NOM(N),NZON(-NSOUT:NREG),NMU
      REAL SIGANG(-6:M),ZMU(NMU)
      DOUBLE PRECISION H(N),EXPT(NMU,N),EXP2(NMU,N)
*----
*  LOCAL VARIABLES
*----
      INTEGER I,NOMI,NUMOLD,NZI,IMU
      DOUBLE PRECISION TAUD
*
      NUMOLD=NOM(1)
      DO I=1,N
         NOMI=NOM(I)
         NZI=NZON(NOMI)
         IF (NZI.LT.0)THEN
            IF (NUMOLD.NE.NOMI)THEN
               DO IMU=1,NMU
                  EXP2(IMU,I)=1.D0-SIGANG(NZI)
               ENDDO
            ELSE
               DO IMU=1,NMU
                  EXP2(IMU,I)=0.D0
               ENDDO
            ENDIF
         ELSE
            DO IMU=1,NMU
               TAUD=SIGANG(NZI)*H(I)*ZMU(IMU)
               EXP2(IMU,I)=2.D0*TAUD/(2.D0+TAUD)
            ENDDO
         ENDIF
         NUMOLD=NOMI
      ENDDO
      CALL XDDSET(EXPT,NMU*N,0.0D0)
*
      RETURN
      END
