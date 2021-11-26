*DECK THMCCD
      REAL FUNCTION THMGCD(TEMP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the product of the heat capacity of cladding (in J/Kg/K) times
* its density (in Kg/m^3).
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal.
*
*Author(s): 
* P. Gallet
*
*Parameters: input
* TEMP    cladding temperature in Kelvin.
* 
*Parameters: output
* THMGCD  product of the heat capacity of the cladding times its density
*         (in J/K/m^3).
*         
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL TEMP
*----
*  LOCAL VARIABLES
*  CP:   cladding heat capacity in J/Kg/K
*  RO:   cladding density with zero porosity in kg/m^3
*----
      REAL CP,RO,DKELV,T0,T1,T2
      PARAMETER (DKELV=273.15,T0=1090.0,T1=1169.0,T2=1243.0)
*
*     calculation of the density of the cladding with the value of the
*     temperature
      RO=6690.0-0.1855*TEMP
*     calculation of the heat capacity of the cladding in J/kg/K
      IF(TEMP.LE.T0) THEN
*        for : T<1090.0
         CP=226.7+0.2066*TEMP-0.6492E-04*TEMP**2.0
      ELSE IF(TEMP.LE.T1) THEN
*        for : 1090<=T<1169.0
         CP=6.94*TEMP-7189.0
      ELSE IF(TEMP.LE.T2) THEN
*        for : 1169<=T<1243.0
         CP=9312.9-7.177*TEMP
      ELSE
*        for T>=1243.0
         CP=356.0
      ENDIF
*     calculation of internal energy of the cladding in J/m^3/K
      THMGCD=RO*CP
      RETURN
      END
