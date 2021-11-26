*DECK THMCDI
      FUNCTION THMCDI(T2K,T1K,BURN,POROS,FRACPU,ICONDF,NCONDF,KCONDF,
     > UCONDF,IFRCDI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the thermal conductivity integral of UOX or MOX fuel.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert, V. Salino
*
*Parameters: input
* T2K     final temperature in Kelvin.
* T1K     initial temperature in Kelvin.
* BURN    fuel burnup in MWday/tonne.
* POROS   fuel porosity.
* FRACPU  plutonium mass fraction in fuel.
* ICONDF  fuel conductivity flag (0=Stora-Chenebault or COMETHE/
*         1=user-provided polynomial + inverse term).
* NCONDF  degree of user-provided fuel conductivity polynomial.
* KCONDF  polynomial coefficients for fuel conductivity in W/m/K^(k+1)
*         (except for the two last coefficients which belongs to the
*         inverse term).
* UCONDF  required unit of temperature in polynomial for fuel
*         conductivity (KELVIN or CELSIUS).
* IFRCDI  flag indicating if average approximation is forced during
*         fuel conductivity evaluation (0=default/1=average
*         approximation forced).
*
*Parameters: output
* THMCDI  thermal conductivity integral in Watt/m/K.
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
      INTEGER ICONDF,NCONDF,IFRCDI
      REAL T1K,T2K,BURN,POROS,FRACPU,KCONDF(NCONDF+3),THMCDI
      CHARACTER UCONDF*12
*----
*  LOCAL VARIABLES
* NPAS    number of rectangles in the quadrature
* DT      rectangle width
* T2T1    temperature difference
* DTMIN   cutoff criterion for selecting the approximation
* FPI     burnup correcting factor
* CIRRA   burnup correction constant
* HV*     coefficients of the Stora-Chenebault correlation
* HK*     coefficients of the Comethe correlation
*----
      INTEGER NPAS,I,K
      REAL T1,T2,DT,TM,DTMIN,T2T1,FPI,TT,TMK,TEMP,FTP,CINT,FP,TTK
      REAL HV1, HV2, HV3
      REAL HK1, HK2, HK4, HK5
      REAL ZKELV,CIRRA
*
      PARAMETER ( ZKELV=273.15 )
      PARAMETER ( HV1= 1.3324E-08 , HV2 = -4.3554E-05 ,
     & HV3 = 5.8915E-02 )
      PARAMETER ( HK1= 40.05 , HK2 = 129.4 , HK4 = 0.8 ,
     & HK5 = 0.6416E-12 )
      PARAMETER ( CIRRA= 0.124E-02 )
*
      REAL A
      DATA NPAS /10/
      DATA DTMIN /10./
*
      IF(MIN(T1K,T2K).LE.0.0) THEN
         CALL XABORT('@THMCDI: NEGATIVE TEMPERATURE.')
      ENDIF
      T1=T1K-ZKELV
      T2=T2K-ZKELV
*
      T2T1 = T2-T1
      DT = T2T1/NPAS
      TM = (T1+T2)/2.0
      IF(ICONDF.EQ.1) THEN
*        User-given conductivity, as a function of temperature
         IF((ABS(T2T1).LT.DTMIN).OR.(IFRCDI.EQ.1)) THEN
*           Use the average value approximation
            THMCDI=0.0
            IF(UCONDF.EQ.'KELVIN') THEN
               TMK = TM + ZKELV
               DO K=1,NCONDF+1
                  THMCDI=THMCDI + KCONDF(K)*TMK**(K-1)
               ENDDO
               THMCDI=THMCDI + KCONDF(NCONDF+2)/(TMK-KCONDF(NCONDF+3))
            ELSE
               DO K=1,NCONDF+1
                  THMCDI=THMCDI + KCONDF(K)*TM**(K-1)
               ENDDO
               THMCDI=THMCDI + KCONDF(NCONDF+2)/(TM-KCONDF(NCONDF+3))
            ENDIF
         ELSE
*           Use the rectangle quadrature approximation
            TT=T1-DT*0.5
            CINT=0.
            DO I=1,NPAS
               TT=TT+DT
               IF(UCONDF.EQ.'KELVIN') THEN
                  TTK = TT + ZKELV
                  DO K=1,NCONDF+1
                     CINT=CINT + KCONDF(K)*TTK**(K-1)
                  ENDDO
                  CINT=CINT + KCONDF(NCONDF+2)/(TTK-KCONDF(NCONDF+3))
               ELSE
                  DO K=1,NCONDF+1
                     CINT=CINT + KCONDF(K)*TT**(K-1)
                  ENDDO
                  CINT=CINT + KCONDF(NCONDF+2)/(TT-KCONDF(NCONDF+3))
               ENDIF
            ENDDO
            THMCDI=CINT/NPAS
         ENDIF
      ELSE IF(FRACPU.GT.0.) THEN
*        Use the Comethe correlation for MOX fuel
         FPI=CIRRA*BURN
         IF((ABS(T2T1).LT.DTMIN).OR.(IFRCDI.EQ.1)) THEN
*           Use the average value approximation
            IF(TM.GT.1000.0) THEN
               A=2.0
            ELSE
               A=2.58-0.58E-03*TM
            ENDIF
            FP=(1.0-A*POROS)/(1.0-A*0.05)
            TMK = TM + ZKELV
            TEMP = HK2 + (1.0 + HK4*FRACPU*1.E-02) * TMK
            FTP = FP * (HK1/TEMP + HK5*TMK*TMK*TMK) *100.0
            IF(TM.EQ.0.) THEN
               THMCDI=FTP
            ELSE
               THMCDI=1.0/(1.0/FTP+FPI/TM)
            ENDIF
         ELSE
*           Use the rectangle quadrature approximation
            TT=T1-DT*0.5
            CINT=0.
            DO I=1,NPAS
               TT=TT+DT
               IF(TT.GT.1000.0) THEN
                  A=2.0
               ELSE
                  A=2.58-0.58E-03*TT
               ENDIF
               FP=(1.0-A*POROS)/(1.0-A*0.05)
               TTK = TT + ZKELV
               TEMP = HK2 + (1.0 + HK4*FRACPU*1.E-02) * TTK
               FTP = FP * (HK1/TEMP  + HK5*TTK*TTK*TTK) *100.0
               IF(TT.EQ.0.0) THEN
                  CINT=CINT+FTP
               ELSE
                  CINT=CINT+1.0/(1.0/FTP+FPI/TT)
               ENDIF
            ENDDO
            THMCDI=CINT/NPAS
         ENDIF
      ELSE
*        Use the Stora-Chenebault correlation for UOX fuel
*        (also called the "HGAP Variable 88" correlation)
         FPI=CIRRA*BURN
         IF((ABS(T2T1).LT.DTMIN).OR.(IFRCDI.EQ.1)) THEN
*           Use the average value approximation
            IF(TM.GT.1000.) THEN
               A=2.0
            ELSE
               A=2.58-0.58E-03*TM
            ENDIF
            FP=(1.0-A*POROS)/(1.0-A*0.034)
            FTP=FP*(HV3+HV2*TM+HV1*TM*TM)*100.0
            IF(TM.EQ.0.0) THEN
               THMCDI=FP*HV3*100.0
            ELSE
               THMCDI=1.0/(1.0/FTP+FPI/TM)
            ENDIF
         ELSE
*           Use the rectangle quadrature approximation
            TT=T1-DT*0.5
            CINT=0.
            DO I=1,NPAS
               TT=TT+DT
               IF(TT.GT.1000.) THEN
                  A=2.0
               ELSE
                  A=2.58-0.58E-03*TT
               ENDIF
               FP=(1.0-A*POROS)/(1.0-A*0.034)
               FTP=FP*(HV3+HV2*TT+HV1*TT*TT)*100.0
               IF(TT.EQ.0.0) THEN
                  CINT=CINT+FP*HV3*100.0
               ELSE
                  CINT=CINT+1.0/(1.0/FTP+FPI/TT)
               ENDIF
            ENDDO
            THMCDI=CINT*DT/T2T1
         ENDIF
      ENDIF
      RETURN
      END
