*DECK LIBSDC
      SUBROUTINE LIBSDC(NBMIX,NGROUP,NBISO,ISONRF,MIX,DEN,MASK,ENER,
     1 KGAS,DENMAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Apply Sternheimer density correction to the collision stopping power.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and A. Naceur
*
*Parameters: input
* NBMIX   number of mixtures present in the calculation domain.
* NGROUP  number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* ISONRF  character*12 reference names of isotopes.
* MIX     mixture number of each isotope (can be zero).
* DEN     density of each isotope.
* MASK    mixture mask (=.true. if a mixture is to be made).
* ENER    energy groups limits.
* KGAS    state of each mixture (=0: solid/liquid; =1: gas).
*
*Parameters: input/output
* DENMAT  Sterheimer density correction (delta).
*
*References:
* [1]. L. J. Lorence Jr., J. E. Morel, and G. D. Valdez, "Physics guide
*      to CEPXS: A multigroup coupled electron-photon cross section
*      generating code," Technical report SAND89-1685, Sandia National
*      Laboratories, Albuquerque, New Mexico 87185 and Livermore,
*      California 94550.
*
* [2]. Sternheimer, R. M. (1956). Density effect for the ionization
*      loss in various materials. Physical Review, 103(3), 511.
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER :: NBMIX,NGROUP,NBISO,MIX(NBISO),KGAS(NBMIX)
      CHARACTER(LEN=12) :: ISONRF(NBISO)
      REAL :: DEN(NBISO),ENER(NGROUP+1) !,ESTOP(NBMIX,NGROUP+1)
      REAL :: DENMAT(NBMIX,NGROUP+1)
      LOGICAL :: MASK(NBMIX)
*----
*  LOCAL VARIABLES
*----
      INTEGER :: I
      DOUBLE PRECISION :: PITT(12)
      CHARACTER(LEN=2) :: ELEMNT(100)
      CHARACTER(LEN=131) :: HSMG
      DOUBLE PRECISION, PARAMETER :: DM=3.0D0 ! Sternheimer exponent
      DOUBLE PRECISION, PARAMETER :: CEMASS=0.510976D0 !Electron mass MeV
      DOUBLE PRECISION, PARAMETER :: C2=0.249467D0 ! 3/8 Thomson xs
      DOUBLE PRECISION :: XDRCST
*----
*  DATA STATEMENTS
*----
      DATA (ELEMNT(I),I=1,100) /
     1 'h',  'he', 'li', 'be', 'b',  'c',  'n',
     2 'o',  'f',  'ne', 'na', 'mg', 'al', 'si',
     3 'p',  's',  'cl', 'ar', 'k',  'ca', 'sc',
     4 'ti', 'v',  'cr', 'mn', 'fe', 'co', 'ni',
     5 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br',
     6 'kr', 'rb', 'sr', 'y',  'zr', 'nb', 'mo',
     7 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in',
     8 'sn', 'sb', 'te', 'i',  'xe', 'cs', 'ba',
     9 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu',
     1 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb',
     2 'lu', 'hf', 'ta', 'w',  're', 'os', 'ir',
     3 'pt', 'au', 'hg', 'tl', 'pb', 'bi', 'po',
     4 'at', 'rn', 'fr', 'ra', 'ac', 'th', 'pa',
     5 'u',  'np', 'pu', 'am', 'cm', 'bk', 'cf',
     6 'es', 'fm' /
      DATA (PITT(I),I=1,12) /
     1  18.7D0, 42.0D0, 38.0D0, 60.0D0, 71.0D0, 78.0D0,
     2  85.0D0, 89.0D0, 92.0D0, 131.0D0, 146.0D0, 156.0D0 /
*----
*  MAIN LOOP OVER MIXTURES
*----
      AVCON=1.0D-24*XDRCST('Avogadro','N/moles')
      DO IBM=1,NBMIX
        IF(MASK(IBM)) THEN
*----
*  CALCULATE THE MEAN IONIZATION ENERGY, EION, IN EV
*----
          EION=0.0D0
          ZZA=0.0D0
          DO ISOT=1,NBISO
            IF((MIX(ISOT).NE.IBM).OR.(DEN(ISOT).EQ.0.0)) CYCLE
            IZ=0
            DO I=1,100
              IF(ISONRF(ISOT)(:2).EQ.ELEMNT(I)) THEN
                IZ=I
                EXIT
              ENDIF
            ENDDO
            IF(IZ.EQ.0) THEN
              WRITE(HSMG,'(40HLIBSDC: UNABLE TO ASSIGN AN ATOMIC NUMBE,
     1        5HR TO ,A,1H.)') ISONRF(ISOT)(:2)
              CALL XABORT(HSMG)
            ENDIF
            WAZ=DEN(ISOT)*REAL(IZ)/AVCON
            IF(IZ.GE.13) THEN
*             for Z > 13, use definition of mean ionization given by
*             Sternheimer
              PIT=(9.76D0+58.8D0*(REAL(IZ)**(-1.19D0)))*REAL(IZ)
            ELSE
*             obtain ionization energy from the data statement
              PIT=PITT(IZ)
            ENDIF
            EION=EION+WAZ*LOG(PIT)
            ZZA=ZZA+WAZ
          ENDDO
          EION = EXP(EION/ZZA)  
*----
*  EVALUATE PLANCK'S CONSTANT TIMES THE PLASMA FREQUENCY IN EV
*----
          HNUP = 28.8D0*SQRT(ZZA)
*----
*  EVALUATE PARAMETERS IN THE STERNHEIMER FORMALISM
*----
          C = - (2.0D0*LOG(EION/HNUP) + 1.0D0)
          IF(KGAS(IBM).EQ.0) THEN
*           The material is a solid/liquid
            IF(EION .GE. 100.D0) THEN
              X1 = 3.0D0
              IF(-C .GE. 5.215D0) THEN
                X0 = -0.326D0 * C - 1.5D0
              ELSE
                X0 = 0.2D0
              ENDIF
            ELSE
              X1 = 2.0D0
              IF(-C .GE. 3.681D0) THEN
                X0 = -0.326D0 * C - 1.0D0
              ELSE
                X0 = 0.2D0
              ENDIF
            ENDIF
          ELSE
*           The material is a gas
            IF(-C .LT. 12.25D0) THEN
              X1 = 4.0D0
              X0 = 2.0D0
              IF(-C .LT. 11.5D0) X0 = 1.9D0
              IF(-C .LT. 11.0D0) X0 = 1.8D0
              IF(-C .LT. 10.5D0) X0 = 1.7D0
              IF(-C .LT. 10.0D0) X0 = 1.6D0
            ELSE IF(-C .GE. 13.804D0) THEN
              X1 = 5.0D0
              X0 = -.326D0 * C - 2.5D0
            ELSE
              X1 = 5.0D0
              X0 = 2.0D0
            ENDIF
          ENDIF
*
          IF(X1.LT.X0) THEN
            WRITE(HSMG,*) 'LIBSDC: NEGATIVE REAL TO REAL POWER. HAVE ',
     1      'YOU NEGLECTED THE "GAS" KEYWORD FOR A GASEOUS MIXTURE?'
            CALL XABORT(HSMG)
          ENDIF
          B = (-C - 4.606D0*X0) / (X1 - X0)**DM
          CONV = 2.0D0*AVCON*C2*CEMASS*ZZA
*----
*  CALCULATE THE DENSITY CORRECTION FACTOR
*----
          DO LLL=1,NGROUP+1
            T = ENER(LLL)/CEMASS/1.0D6
            T1 = T + 1.0D0
            T2 = T + 2.0D0
            BSQ = T * T2 / T1**2
*----
*  EVALUATE THE ELECTRON'S MOMENTUM / (MASS OF THE ELECTRON * C)
*----
            PMC = SQRT(2.0D0*T + T*T)
*----
*  EVALUTE THE ENERGY PARAMETER IN THE STERNHEIMER FORMALISM
*----
            X = LOG10(PMC)
            IF(X.LE.X0) THEN
              DS = 0.0D0
            ELSE IF(X.LE.X1) THEN
              DS = 4.606D0*X+C+B*((X1-X)**DM)
            ELSE
              DS = 4.606D0*X+C
            ENDIF
            IF(DS.LT.0.0) DS=0.0D0
            DENMAT(IBM,LLL)= REAL(CONV*DS/BSQ)
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END
