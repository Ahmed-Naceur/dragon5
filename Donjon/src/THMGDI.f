*DECK THMGDI
      FUNCTION THMGDI(T2K,T1K,ICONDC,NCONDC,KCONDC,UCONDC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the average thermal conductivity of the cladding
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal.
*
*Author(s): 
* P. Gallet, V. Salino
*
*Parameters: input
* T2K     final temperature in Kelvin.
* T1K     initial temperature in Kelvin.
* ICONDC  clad conductivity flag (0=default/1=user-provided
*         polynomial).
* NCONDC  degree of user-provided clad conductivity polynomial.
* KCONDC  polynomial coefficients for clad conductivity in W/m/K^(k+1).
* UCONDC  required unit of temperature in polynomial for clad
*         conductivity (KELVIN or CELSIUS).
*
*Parameters: output
* THMGDI  thermal conductivity of the cladding in W/m/K.
*
*Reference:
* A. Poncot, "Assimilation de donnees pour la dynamique du xenon dans
* les coeurs de centrale nucleaire", Ph.D Thesis, Universite de
* Toulouse, France, 2008.
* 
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ICONDC,NCONDC
      REAL T1K,T2K,KCONDC(NCONDC+1),THMGDI
      CHARACTER UCONDC*12
*----
*  LOCAL VARIABLES
*----
      INTEGER K
      REAL T1,T2,TM,TMK,ZKELV
*
      PARAMETER ( ZKELV=273.15 )
*
      IF(MIN(T1K,T2K).LE.0.0) THEN
         CALL XABORT('@THMGDI: NEGATIVE TEMPERATURE.')
      ENDIF
      T1=T1K-ZKELV
      T2=T2K-ZKELV
*
      TM=(T1+T2)*0.5
      IF(ICONDC.EQ.1) THEN
*        User-given conductivity, as a polynomial of temperature
         THMGDI=0.0
         IF(UCONDC.EQ.'KELVIN') THEN
            TMK = TM + ZKELV
            DO K=1,NCONDC+1
               THMGDI = THMGDI + KCONDC(K)*TMK**(K-1)
            ENDDO
         ELSE
            DO K=1,NCONDC+1
               THMGDI = THMGDI + KCONDC(K)*TM**(K-1)
            ENDDO
         ENDIF
      ELSE
*        thermal conductivity of the cladding in W/m/K
         THMGDI=12.0+1.25E-2*TM
      ENDIF
      
      RETURN
      END
