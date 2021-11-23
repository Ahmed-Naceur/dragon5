*DECK THMCCD
      REAL FUNCTION THMCCD(TEMP,POROS,FRACPU)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the product of the heat capacity of fuel (in J/Kg/K) times
* its density (in Kg/m^3).
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal.
*
*Author(s): 
* P. Gallet
*
*Parameters: input
* TEMP    fuel temperature in Kelvin.
* POROS   oxyde porosity.
* FRACPU  plutonium mass fraction in fuel.
* 
*Parameters: output
* THMCCD  product of the heat capacity of fuel times its density
*         (in J/K/m^3).
* 
*Reference:
* J. J. Carbajo, G. L. Yoder, S. G. Popov and V. K. Ivanov, "A review of
* the thermophysical properties of MOX and UO2 fuels," J. of Nuclear
* Materials, 299, 181-198 (2001).
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL TEMP,POROS,FRACPU
*----
*  LOCAL VARIABLES
*  CP:    heat capacity in J/Kg/K
*  DENS0: fuel density with zero porosity
*  ROURA: uranium density in Kg/m^3
*  ROPLU: plutonium density in Kg/m^3
*----
      REAL CP,DENS0,RO,ROURA,ROPLU,A1,A2,A3,A4,CORR,T2,T1,C1U,C2U,C3U,
     > C4U,C5U,C6U,C1PU,C2PU,C3PU,C4PU,C5PU,C6PU,CPU,CPPU
      PARAMETER (ROURA=10970.0,ROPLU=11460.0,A1=0.99672,A2=1.179E-05,
     > A3=-2.429E-09,A4=1.219E-12,C1U=193.238,C2U=325.7294,
     > C3U=-312.0042,C4U=116.8224,C5U=-9.7535,C6U=-2.6441,C1PU=311.7866,
     > C2PU=39.258,C3PU=-2.256,C4PU=0.0,C5PU=0.0,C6PU=-7.0131)
*
      T2=MAX(0.0,TEMP)
      T1=T2/1000.0
*     temperature correction coefficient for density calculation
      CORR=1.0/(A1+A2*T2+A3*T2**2.0+A4*T2**3.0)**3.0
      IF(FRACPU.EQ.0.0) THEN
*        UOX
*        density of the UOX fuel
         RO=(1.0-POROS)*ROURA*CORR
*        heat capacity of the UOX fuel
         CPU=C1U+C2U*T1+C3U*T1**2.0+C4U*T1**3.0+C5U*T1**4.0+C6U
     >       /(T1**2.0)
         CPPU=0.00
         CP=CPU
      ELSE
*        MOX
*        density of the MOX fuel
         DENS0=100.0*CORR/((FRACPU/ROPLU)+((100.0-FRACPU)/ROURA))
         RO=(1.-POROS)*DENS0
*        heat capacity of the MOX fuel
         CPU=C1U+C2U*T1+C3U*T1**2.0+C4U*T1**3.0+C5U*T1**4.0+C6U
     >       /(T1**2.0)
         CPPU=C1PU+C2PU*T1+C3PU*T1**2.0+C6PU/(T1**2.0)
         CP=((100.0-FRACPU)*CPU+FRACPU*CPPU)/100.0
      ENDIF
*     total internal energy of the fuel
      THMCCD=RO*CP
      RETURN
      END
