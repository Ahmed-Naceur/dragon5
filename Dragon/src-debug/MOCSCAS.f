*DECK MOCSCAS
      SUBROUTINE MOCSCAS(N,NREG,NSOUT,M,NOM,NZON,H,SIGANG,EXPT,EXP2,
     1           NMU,ZMU)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculate coefficients of a track for the cyclic characteristics 
* integration: Step-Characteristics scheme with tabulated exponential
* and 'source term isolation' option turned off.
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
      DOUBLE PRECISION H(N),EXPT(NMU,N),EXP2(2,NMU,N)
*----
*  LOCAL VARIABLES
*----
      INTEGER I,NOMI,NUMOLD,NZI,IMU
      REAL TAU
      DOUBLE PRECISION TAUDMIN,HID,TAUD,TEMP,HID2,TAUD3,TAUD4,TAUD5
*     tabulated exponential common block
      REAL             E0, E1, PAS1, DX1, XLIM1
      INTEGER          MEX1, LAU
      PARAMETER      ( MEX1=7936, TAUDMIN=2.D-2 )
      COMMON /EXP1/ E0(0:MEX1),E1(0:MEX1),PAS1,DX1,XLIM1
*
      NUMOLD=NOM(1)
      DO I=1,N
         NOMI=NOM(I)
         NZI=NZON(NOMI)
         IF(NZI.LT.0) THEN
            IF(NUMOLD.NE.NOMI) THEN
               DO IMU=1,NMU
                  EXP2(1,IMU,I)=SIGANG(NZI)
                  EXPT(IMU,I)=EXP2(1,IMU,I)
               ENDDO
            ELSE
               DO IMU=1,NMU
                  EXP2(1,IMU,I)=1.D0
                  EXPT(IMU,I)=EXP2(1,IMU,I)
               ENDDO
            ENDIF
         ELSE
            DO IMU=1,NMU
               HID=DBLE(H(I)*ZMU(IMU))
               TAUD=SIGANG(NZI)*HID
               TAU=REAL(TAUD)
               IF(TAU.GE.XLIM1) THEN
*              Out of the table range
                  EXPT(IMU,I)=0.D0
                  TEMP=1.D0/TAUD
                  EXP2(1,IMU,I)=TEMP*HID
                  EXP2(2,IMU,I)=HID*(1.D0-TEMP)/DBLE(SIGANG(NZI))
               ELSE
*              Linear interpolation in table of (1-exp(-x))/x
                  LAU=INT(TAU*PAS1)
                  TEMP=DBLE(E0(LAU)+E1(LAU)*TAU)
                  EXPT(IMU,I)=1.D0-TEMP*TAUD
                  EXP2(1,IMU,I)=TEMP*HID
                  IF(TAUD.LE.TAUDMIN) THEN 
*                 and expansion in Taylor serie in O(TAUD^3)
                     TAUD3=TAUD/3.D0
                     TAUD4=0.125D0*TAUD
                     TAUD5=0.2D0*TAUD
                     HID2=HID*HID
                     EXP2(2,IMU,I)=HID2*(0.5D0-TAUD3*(0.5D0-TAUD4
     1                                   *(1.D0-TAUD5)))
                  ELSE
                     EXP2(2,IMU,I)=HID*(1.D0-TEMP)/DBLE(SIGANG(NZI))
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
         NUMOLD=NOMI
      ENDDO
*
      RETURN
      END
