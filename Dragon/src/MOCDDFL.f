*DECK MOCDDFL
      SUBROUTINE MOCDDFL(N,NREG,NSOUT,M,NOM,NZON,H,SIGANG,DSIG,EXPT,
     1           EXP2,NMU,ZMU)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculate coefficients of a track for the cyclic characteristics 
* integration: DD1 Diamond-Differencing scheme without fix-up and
* 'source term isolation' option turned off.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
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
* DSIG    macroscopic total cross sections.
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
      DOUBLE PRECISION H(N),DSIG(N),EXPT(NMU,N),EXP2(5,NMU,N)
*----
*  LOCAL VARIABLES
*----
      INTEGER I,NOMI,NUMOLD,NZI,IMU
      DOUBLE PRECISION SQ3,TAUD,HID,H2,H3,DEN
*
      SQ3=SQRT(3.0D0)
      NUMOLD=NOM(1)
      DO I=1,N
         NOMI=NOM(I)
         NZI=NZON(NOMI)
         DSIG(I)=SIGANG(NZI)
         IF(NZI.LT.0) THEN
            DO IMU=1,NMU
               EXP2(2,IMU,I)=0.D0
               EXP2(3,IMU,I)=0.D0
               EXP2(5,IMU,I)=0.D0
            ENDDO
            IF(NUMOLD.NE.NOMI) THEN
               DO IMU=1,NMU
                  EXP2(1,IMU,I)=SIGANG(NZI)
                  EXP2(4,IMU,I)=EXP2(1,IMU,I)
                  EXPT(IMU,I)=EXP2(1,IMU,I)
               ENDDO
            ELSE
               DO IMU=1,NMU
                  EXP2(1,IMU,I)=1.D0
                  EXP2(4,IMU,I)=1.D0
                  EXPT(IMU,I)=EXP2(1,IMU,I)
               ENDDO
            ENDIF
         ELSE
            DO IMU=1,NMU
               HID=DBLE(H(I)*ZMU(IMU))
               TAUD=SIGANG(NZI)*HID
               H2=HID*HID
               H3=H2*HID
               DEN=(TAUD+6.0D0)*TAUD+12.0D0
               EXPT(IMU,I)=((TAUD-6.0D0)*TAUD+12.0D0)/DEN
               EXP2(1,IMU,I)=12.0D0*HID/DEN
               EXP2(2,IMU,I)=H2*(TAUD+6.0D0)/DEN
               EXP2(3,IMU,I)=-2.0D0*SQ3*H3/DEN
               EXP2(4,IMU,I)=-2.0D0*SQ3*TAUD/DEN
               EXP2(5,IMU,I)=H3/DEN
            ENDDO
         ENDIF
         NUMOLD=NOMI
      ENDDO
*
      RETURN
      END
