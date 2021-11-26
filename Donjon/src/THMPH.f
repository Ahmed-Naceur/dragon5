*DECK THMPH
      SUBROUTINE THMPH(IFLUID,PP,HH,RHO,TEMP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Backwards inversion of steam tables to find water density and
* temperature.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert, P. Gallet
*
*Parameters: input
* IFLUID  type of fluid (0=H2O; 1=D2O).
* PP      pressure (Pa)
* HH      enthalpy (J/Kg)
*
*Parameters: output
* RHO     water density (Kg/m^3)
* TEMP    temperature (K)
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IFLUID
      REAL PP,HH,RHO,TEMP
*----
*  LOCAL VARIABLES
*----
      PARAMETER(DT=-0.01,ZKELV=273.15,S=1.0)
      REAL A(15),R1,R3,R4,R5,RV,HSAT,HV,XTH
      DATA A/
     > 0.2873E+03,-0.5098E+00,-0.3459E+00,0.1910E+00,-0.2840E-01,
     > 0.8266E+02,0.1141E+01,-0.2724E+01,0.1077E+00,-0.1144E+02,
     > 0.9500E+01,-0.2715E+01,-0.1290E+02,0.9148E+01,-0.8093E+01 /
*----
*  INITIAL APPROXIMATION OF T1
*----
      IF(IFLUID.EQ.0) THEN
         CALL THMSAT(PP,TSAT)
         CALL THMTX(TSAT,1.0,RV,HV,R3,R4,R5)
         CALL THMTX(TSAT,0.0,R1,HSAT,R3,R4,R5)
      ELSE IF(IFLUID.EQ.1) THEN
         CALL THMHST(PP,TSAT)
         CALL THMHTX(TSAT,1.0,RV,HV,R3,R4,R5)
         CALL THMHTX(TSAT,0.0,R1,HSAT,R3,R4,R5)
      ENDIF
      IF((ABS(HSAT-HH)/HSAT).LE.1.0E-5) THEN
         T1=TSAT
         GO TO 20
      ELSEIF(HH.LE.HSAT) THEN
         T=(HH-1270.0E3)/420.0E3
         P=(PP-140.0E5)/30.0E5
         H1=A(1)+P*(A(2)+P*(A(3)+P*(A(4)+P*A(5))))
         H2=A(6)+P*(A(7)+P*(A(8)+P*(A(9))))
         H3=A(10)+P*(A(11)+P*(A(12)))
         H4=A(13)+P*A(14)
         H5=A(15)
         T1=H1+T*(H2+T*(H3+T*(H4+T*(H5))))+ZKELV
* INLET TEMPERATURE WAS VERIFIED TO BE GREATER THAN 0 C. T1 INITIAL
* GUESS LOWER THAN THAT SHOULD BE INTERPRETED AS FLAWED (FAR FROM
* FITTING REGION). CORRECTING WITH AN ABOVE-0 C GUESS.
         IF(T1.LT.ZKELV) T1=10.0+ZKELV
      ELSEIF(HH.LE.HV) THEN
*        saturated steam
         TEMP=TSAT
         XTH=(HH-HSAT)/(HV-HSAT)
         RHO=1.0/(XTH/RV+(1.0-XTH)/R1)
         GO TO 30
      ELSE
*        superheated steam
         T1=TSAT
      ENDIF
*----
*  NEWTON ITERATIONS
*----
      ITER=0
   10 ITER=ITER+1
      IF(ITER.GT.30) CALL XABORT('THMPH: CONVERGENCE FAILURE.')
      IF(IFLUID.EQ.0) THEN
         CALL THMPT(PP,T1,R1,H1,R3,R4,R5)
         CALL THMPT(PP,T1+DT,R1,H1P,R3,R4,R5)
      ELSE IF(IFLUID.EQ.1) THEN
         CALL THMHPT(PP,T1,R1,H1,R3,R4,R5)
         CALL THMHPT(PP,T1+DT,R1,H1P,R3,R4,R5)
      ENDIF
      IF(ABS((HH-H1)/HH).LT.1.E-05) GO TO 20
      T1=T1+(HH-H1)*DT/(H1P-H1)
      IF((HH.LE.HSAT).AND.(T1.GE.TSAT)) T1=TSAT
      IF((HH.GE.HV).AND.(T1.LE.TSAT)) T1=TSAT
      GO TO 10
   20 RHO=R1
      TEMP=T1
   30 RETURN
      END
