*DECK THMH2O
      SUBROUTINE THMH2O(ITIME,I,J,K,K0,PINLET,MFLOW,HMAVG,ENT,HD,IFLUID,
     > IHCONV,KHCONV,ISUBM,RADCL,ZF,PHI,XFL,EPS,SLIP,ACOOL,PCH,DZ,TCALO,
     > RHO,RHOLAV,TSCLAD,KWA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Nucleate boiling correlations along a single coolant channel.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert, P. Gallet
*
*Parameters: input
* ITIME   type of calculation  (0=steady-state; 1=transient).
* I       position of channel alon X-axis
* J       position of channel alon Y-axis
* K       position of channel alon Z-axis
* K0      onser of nuclear boiling point
* PINLET  pressure in Pascal
* MFLOW   massic coolant flow rate in Kg/m^2/s
* HMAVG   averaged enthalpy
* ENT     four values of enthalpy in J/Kg to be used in Gaussian
*         integration
* HD      hydraulic diameter in m
* IFLUID  type of fluid (0=H2O; 1=D2O).
* IHCONV  flag indicating HCONV chosen (0=default/1=user-provided).
* KHCONV  fixed user-provided HCONV value in W/m^2/K.
* ISUBM   subcooling model (0: one-phase; 1: Jens-Lottes model;
*         2: Saha-Zuber model).
* RADCL   outer clad radius in m
* ZF      parameters used to compute heat flux on clad surface in
*         transient cases.
* PHI     heat flow exchanged between clad and fluid in W/m^2.
*         Given in steady-state cases.
* XFL     input coolant flow quality
* EPS     input coolant void fraction
* SLIP    input slip ratio of vapor phase speed to liquid phase speed.
* ACOOL   coolant cross section area in m^2.
* PCH     heating perimeter in m.
* DZ      axial mesh width in m.
*
*Parameters: output
* PHI     heat flow exchanged between clad and fluid in W/m^2.
*         Computed in transient cases.
* XFL     output coolant flow quality
* EPS     output coolant void fraction
* SLIP    output slip ratio of vapor phase speed to liquid phase speed.
* TCALO   coolant temperature in K
* RHO     coolant density in Kg/m^3
* RHOLAV  liquid density in kg/m^3
* TSCLAD  clad temperature in K
* KWA     flow regime (=0: single-phase; =1: subcooled; =2: nucleate
*         boiling; =3 superheated steam)
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER I,J,K,K0,IFLUID,IHCONV,ISUBM,KWA
      REAL PINLET,MFLOW,HMAVG,ENT(4),HD,KHCONV,RADCL,ZF(2),PHI,TCALO,
     > RHO,RHOLAV,TSCLAD,XFL,EPS,SLIP,ACOOL,PCH,DZ
*----
*  LOCAL VARIABLES
*----
      REAL W(4),HL(4),JL,JG
      CHARACTER HSMG*131
      LOGICAL LFIRST
*----
*  SAVE VARIABLES
*----
      SAVE DHSUB,DSAT,W
      DATA W /0.347855,0.652145,0.652145,0.347855/
*----
*  COMPUTE THE PROPERTIES OF THE SATURATED STEAM
*----
      IF(HMAVG.LT.0.0) CALL XABORT('THMH2O: NEGATIVE INPUT ENTHALPY.')
      IF(IFLUID.EQ.0) THEN
         CALL THMSAT(PINLET,TSAT)
         CALL THMTX(TSAT,0.0,RHOL,HLSAT,ZKL,ZMUL,CPL)
         CALL THMTX(TSAT,1.0,RHOG,HGSAT,ZKG,ZMUG,CPG)
      ELSE IF(IFLUID.EQ.1) THEN
         CALL THMHST(PINLET,TSAT)
         CALL THMHTX(TSAT,0.0,RHOL,HLSAT,ZKL,ZMUL,CPL)
         CALL THMHTX(TSAT,1.0,RHOG,HGSAT,ZKG,ZMUG,CPG)
      ENDIF
*----
*  COMPUTE THE DENSITY AND TEMPERATURE OF THE LIQUID
*----
      HL(1)=MIN1(ENT(1),HLSAT)
      HL(2)=MIN1(ENT(2),HLSAT)
      HL(3)=MIN1(ENT(3),HLSAT)
      HL(4)=MIN1(ENT(4),HLSAT)
      CALL THMPH(IFLUID,PINLET,HL(1),R11,TL1)
      CALL THMPH(IFLUID,PINLET,HL(2),R11,TL2)
      CALL THMPH(IFLUID,PINLET,HL(3),R11,TL3)
      CALL THMPH(IFLUID,PINLET,HL(4),R11,TL4)
      IF(IFLUID.EQ.0) THEN
        CALL THMPT(PINLET,TL1,RHO1,R2,R3,R4,CP1)
        CALL THMPT(PINLET,TL2,RHO2,R2,R3,R4,CP2)
        CALL THMPT(PINLET,TL3,RHO3,R2,R3,R4,CP3)
        CALL THMPT(PINLET,TL4,RHO4,R2,R3,R4,CP4)
        IF(ABS(TSAT-TL1).LT.0.1) CALL THMTX(TSAT,0.0,RHO1,R2,R3,R4,CP1)
        IF(ABS(TSAT-TL2).LT.0.1) CALL THMTX(TSAT,0.0,RHO2,R2,R3,R4,CP2)
        IF(ABS(TSAT-TL3).LT.0.1) CALL THMTX(TSAT,0.0,RHO3,R2,R3,R4,CP3)
        IF(ABS(TSAT-TL4).LT.0.1) CALL THMTX(TSAT,0.0,RHO4,R2,R3,R4,CP4)
      ELSE IF(IFLUID.EQ.1) THEN
        CALL THMHPT(PINLET,TL1,RHO1,R2,R3,R4,CP1)
        CALL THMHPT(PINLET,TL2,RHO2,R2,R3,R4,CP2)
        CALL THMHPT(PINLET,TL3,RHO3,R2,R3,R4,CP3)
        CALL THMHPT(PINLET,TL4,RHO4,R2,R3,R4,CP4)
        IF(ABS(TSAT-TL1).LT.0.1) CALL THMHTX(TSAT,0.0,RHO1,R2,R3,R4,CP1)
        IF(ABS(TSAT-TL2).LT.0.1) CALL THMHTX(TSAT,0.0,RHO2,R2,R3,R4,CP2)
        IF(ABS(TSAT-TL3).LT.0.1) CALL THMHTX(TSAT,0.0,RHO3,R2,R3,R4,CP3)
        IF(ABS(TSAT-TL4).LT.0.1) CALL THMHTX(TSAT,0.0,RHO4,R2,R3,R4,CP4)
      ENDIF
      TL=0.5*(W(1)*TL1+W(2)*TL2+W(3)*TL3+W(4)*TL4)
      RHOLAV=0.5*(W(1)*RHO1+W(2)*RHO2+W(3)*RHO3+W(4)*RHO4)
      CPLAV=0.5*(W(1)*CP1+W(2)*CP2+W(3)*CP3+W(4)*CP4)
*----
*  COMPUTE THE STEAM FLOW QUALITY AND LIQUID ENTHALPY
*  Reference: R. T. Lahey Jr. and F. J. Moody, "The thermal hydraulics
*  of a Boiling water nuclear reactor," American Nuclear Society, 1977.
*  Equation (5.177), page 224
*  F2: Thermodynamic quality
*----
      TSCLAD=600.0
      IF(K0.GT.0) TSCLAD=TSAT+DSAT
      XFL0=XFL
      EPS0=EPS
      SLIP0=SLIP
      LFIRST=.TRUE.
      HLAVG=HMAVG
      F2=0.0
      F3=0.0
      IF(K0.GT.0) THEN
        HLV=HGSAT-HLSAT
        IF((HLV.GT.0.0).AND.(DHSUB.GT.0.0)) THEN
          F2=(HMAVG-HLSAT)/HLV
          F3=(DHSUB/HLV)*EXP(-(HMAVG-HLSAT)/DHSUB-1.0)
        ENDIF
        IF(HMAVG.GE.HGSAT) THEN
          XFL=1.0
          EPS=1.0
          SLIP=1.0
          HLAVG=0.0
        ELSE
          IF(ISUBM.EQ.1) THEN
*           Use the Paul Gallet thesis model.
            PI=RHOLAV*CPLAV*(TSCLAD-TL)/(RHOG*HLV)
            XFL=XFL0+PCH*PHI*DZ/(MFLOW*ACOOL*HLV)/(1.0+PI)
          ELSE IF(ISUBM.EQ.2) THEN
*           Use a profile fit model.
            XFL=MAX(XFL0,(F2+F3)/(1.0+F3))
          ENDIF
          HLAVG=MIN(HLSAT,(HMAVG-XFL*HGSAT)/(1.0-XFL))
        ENDIF
*----
*  RECOMPUTE THE LIQUID PROPERTIES
*----
        IF(HLAVG.GT.0.0) THEN
          CALL THMPH(IFLUID,PINLET,HLAVG,RHOL,TL)
          IF(IFLUID.EQ.0) THEN
            CALL THMPT(PINLET,TL,R1,R2,ZKL,ZMUL,CPL)
          ELSE IF(IFLUID.EQ.1) THEN
            CALL THMHPT(PINLET,TL,R1,R2,ZKL,ZMUL,CPL)
          ENDIF
*----
*  COMPUTE THE COOLANT VOID FRACTION AND SLIP RATIO
*  A drift-flux model is proposed by means of the concentration
*  parameter CO and the drift velocity VGJ (their correspondent 
*  correlations are supposed to work fine under different flow regimes.
*----
          IF(HGSAT.GT.HLSAT) THEN
            CO=1.13
            PR=PINLET/10**6
            SIGM=-7.2391E-6*PR**3+2.8345E-4*PR**2-5.1566E-3*PR+4.2324E-2
            VGJ=1.18*((SIGM*9.81*(RHOL-RHOG))/RHOL**2)**0.25
            F4=CO*((XFL*RHOL)+((1.0-XFL)*RHOG))+(RHOL*RHOG*VGJ/MFLOW)
            EPS=(XFL*RHOL)/F4
            JL=(1.0-XFL)*MFLOW/RHOL
            JG=XFL*MFLOW/RHOG
            IF(EPS.NE.0) SLIP=JG*(1.0-EPS)/(JL*EPS)
          ENDIF
        ELSE
*         superheated steam
          CALL THMPH(IFLUID,PINLET,HMAVG,RHOG,TCALO)
          IF(IFLUID.EQ.0) THEN
            CALL THMPT(PINLET,TCALO,R1,R2,ZKG,ZMUG,CPG)
          ELSE IF(IFLUID.EQ.1) THEN
            CALL THMHPT(PINLET,TCALO,R1,R2,ZKG,ZMUG,CPG)
          ENDIF
        ENDIF
      ENDIF
*----
*  COMPUTE THE FLUID PROPERTIES
*  RHO: fluid density
*  REL: Reynolds number of liquid phase
*  PRL: Prandtl number of liquid phase
*----
      IF(XFL.EQ.0.0) THEN
*       One phase liquid
        TB=TSAT-0.1
        IF(TL.LT.TB) THEN
          TCALO=TL
        ELSE
          TCALO=TB
        ENDIF
        IF(IFLUID.EQ.0) THEN
          CALL THMPT(PINLET,TCALO,R1,R2,ZKONE,ZMUONE,CPONE)
        ELSE IF(IFLUID.EQ.1) THEN
          CALL THMHPT(PINLET,TCALO,R1,R2,ZKONE,ZMUONE,CPONE)
        ENDIF
        RHO=RHOLAV
        REL=MFLOW*HD/ZMUONE
        PRL=ZMUONE*CPONE/ZKONE
      ELSE IF(HMAVG.LT.HGSAT) THEN
*       Two-phase flow
        TCALO=EPS*TSAT+(1.0-EPS)*TL
        ZKONE=ZKL
        CPONE=CPL
        RHO=EPS*RHOG+(1.0-EPS)*RHOL
        REL=MFLOW*(1.0-XFL)*HD/ZMUL
        PRL=ZMUL*CPL/ZKL
      ELSE
*       superheated steam
        RHO=RHOG
        REL=MFLOW*HD/ZMUG
        PRL=ZMUG*CPG/ZKG
      ENDIF
*----
*  THERMAL EXCHANGE BETWEEN CLAD AND FLUID USING THE DITTUS AND BOELTER
*  CORRELATION (SINGLE PHASE) OR CHEN CORRELATION (SATURATED BOILING)
*----
      IF(IHCONV.EQ.0) THEN
        ITER=0
        KWA=99
        DO
          ITER=ITER+1
          IF(ITER.GT.50) THEN
            WRITE(HSMG,'(30HTHMH2O: HCONV FAILURE IN SLICE,I5,1H.)') K
            CALL XABORT(HSMG)
          ENDIF
          HA=0.023*(ZKONE/HD)*REL**0.8*PRL**0.4
          F=1.0
          S=1.0
          IF((XFL.EQ.XFL0).OR.(TSCLAD.LE.TSAT).OR.(KWA.EQ.0)) THEN
*           Single-phase convection. Use Dittus-Boelter correlation
            KWA=0
            HB=0.0
            K0=0
            XFL=XFL0
            EPS=EPS0
            SLIP=SLIP0
          ELSE IF(HMAVG.LT.HGSAT) THEN
*           Subcooled convection. Use Dittus-Boelter and Forster-Zuber
*           correlations
*           XM: Martinelli parameter
*           F: Reynolds number factor
*           S: nucleate boiling suppression factor
*           SIGM: surface tension in N/m
*           HA: Dittus-Boelter coefficient
*           HB: Forster-Zuber coefficient
*
            IF(HMAVG.LT.HLSAT) THEN
              KWA=1
            ELSE
              KWA=2
            ENDIF
            XM=(XFL/(1.0-XFL))**0.9*(RHOL/RHOG)**0.5*(ZMUG/ZMUL)**0.1
            IF(XM.LE.0.100207) THEN
              F=1.0
            ELSE
              F=2.35*(0.213+XM)**0.736
            ENDIF
            RE=REL*F**1.25
            S=1.0/(1.0+2.53E-6*RE**1.17)
            PR=PINLET/10**6
            SIGM=-7.2391E-6*PR**3+2.8345E-4*PR**2-5.1566E-3*PR+4.2324E-2
            HA=0.023*(ZKL/HD)*REL**0.8*PRL**0.4
            DTSAT=TSCLAD-TSAT
            IF(IFLUID.EQ.0) THEN
              CALL THMSAP(PW, TSCLAD)
            ELSE
              CALL THMHSP(PW, TSCLAD)
            ENDIF
            DP=PW-PINLET
*           Forster-Zuber equation
            HLV=HGSAT-HLSAT
            HB=0.00122*ZKL**0.79*CPL**0.45*RHOL**0.49/(ZMUL**0.29*
     >      SIGM**0.5*HLV**0.24*RHOG**0.24)*DTSAT**0.24*DP**0.75
          ELSE
*           Superheated steam. Use Mokry correlation
            KWA=3
            IF(IFLUID.EQ.0) THEN
              CALL THMPT(PINLET,TSCLAD,RHOW,R2,R3,R4,R5)
            ELSE IF(IFLUID.EQ.1) THEN
              CALL THMHPT(PINLET,TSCLAD,RHOW,R2,R3,R4,R5)
            ENDIF
            HA=0.0061*(ZKG/HD)*REL**0.904*PRL**0.684*(RHOW/RHO)**0.564
            HB=0.0
          ENDIF
*         Chen correlation
          HCONV=F*HA+S*HB
          IF(HCONV.LE.0.0) THEN
            WRITE(HSMG,'(34HTHMH2O: DRY OUT REACHED IN CHANNEL,3I5)')
     >      I,J,K
            CALL XABORT(HSMG)
          ENDIF
          IF(ITIME.EQ.0) THEN
            TWAL=(PHI+S*HB*TSAT+F*HA*TCALO)/(S*HB+F*HA)
          ELSE
            ZNUM=ZF(1)+RADCL*S*HB*TSAT+RADCL*F*HA*TCALO
            ZDEN=ZF(2)+RADCL*S*HB+RADCL*F*HA
            TWAL=MAX(273.15,ZNUM/ZDEN)
            PHI=MAX(0.0,(ZF(1)-TWAL*ZF(2))/RADCL)
          ENDIF
          IF(ABS(TSCLAD-TWAL).GT.1.0E-5*TSCLAD) THEN
            TSCLAD=TWAL
          ELSE
            EXIT
          ENDIF
        ENDDO
      ELSE IF(IHCONV.EQ.1) THEN
        IF(ITIME.EQ.0) THEN
          TSCLAD=TCALO+PHI/KHCONV
        ELSE
          RCHC=RADCL*KHCONV
          TSCLAD=MAX(273.15,(ZF(1)+RCHC*TCALO)/(ZF(2)+RCHC))
          PHI=(ZF(1)-TSCLAD*ZF(2))/RADCL
        ENDIF
      ENDIF
*----
*  COMPUTE INITIAL BULK LIQUID ENTHALPY SUBCOOLING DHSUB
*----
      IF((ISUBM.GT.0).AND.(K0.EQ.0).AND.LFIRST) THEN
        DTSUB=0.0
        IF(ISUBM.EQ.1) THEN
*         Bowring correlation
*         Reference: R. W. Bowring, "Physical Model, Based on Bubble
*         Detachment, and Calculation of Steam Voidage in the Subcooled
*         Region of a Heated Channel," OECD Report HPR-10, 1962.
*         Equation3 (3) and (17)
          VC=MFLOW/RHOL
          ETA=14.0+0.1*PINLET/1.01325E+05
          DTSUB=ETA*PHI/VC*1.0E-6
         ELSE IF(ISUBM.EQ.2) THEN
*         Saha-Zuber subcooling model
*         PE: Peclet number
          PE=MFLOW*CPL*HD/ZKL
          IF(PE.LE.70000.0) THEN
            DTSUB=PHI*HD/(455.0*ZKL)
          ELSE
*           reactor conditions
            DTSUB=154.0*PHI/(MFLOW*CPL)
          ENDIF
        ENDIF
        IF(TCALO.GE.TSAT-DTSUB) K0=K
        DSAT=TSCLAD-TCALO-DTSUB
        DHSUB=CPL*DTSUB
        LFIRST=.FALSE.
      ENDIF
      RETURN
      END
