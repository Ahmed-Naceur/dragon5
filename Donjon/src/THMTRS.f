*DECK THMTRS
      SUBROUTINE THMTRS(MPTHMI,MPTHM,IMPX,IX,IY,NZ,XBURN,VOLXY,HZ,DTIME,
     > CFLUX,POROS,FNFU,NFD,NDTOT,IFLUID,FCOOL,FFUEL,ACOOL,HD,PCH,
     > MAXITC,MAXIT1,MAXITL,ERMAXT,ERMAXC,SPDIN,TINLET,POULET,FRACPU,
     > ICONDF,NCONDF,KCONDF,UCONDF,ICONDC,NCONDC,KCONDC,UCONDC,IHGAP,
     > KHGAP,IHCONV,KHCONV,WTEFF,IFRCDI,ISUBM,FRO,POW,TCOMB,DCOOL,TCOOL,
     > TSURF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver of the transient thermal-hydraulics module for a single time
* iteration
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal.
*
*Author(s): 
* P. Gallet and A. Hebert
*
*Parameters: input
* MPTHMI  directory of the THM object containing steady-state
*         thermohydraulics data at t-1.
* MPTHM   directory of the THM object containing steady-state
*         thermohydraulics data at t.
* IMPX    printing index (=0 for no print).
* IX      position of mesh along X direction.
* IY      position of mesh along Y direction.
* NZ      number of meshes along Z direction (channel direction).
* XBURN   burnup distribution in MWday/tonne.
* VOLXY   mesh area in the radial plane.
* HZ      Z-directed mesh widths.
* DTIME   time step in s.
* CFLUX   critical heat flux in W/m^2.
* POROS   oxyde porosity.
* FNFU    number of active fuel rods in the fuel bundle.
* NFD     number of discretisation points in fuel regions.
* NDTOT   number of total discretization points in the the fuel
*         pellet and the cladding.
* IFLUID  type of fluid (0=H2O; 1=D2O).
* FCOOL   power density fraction in coolant.
* FFUEL   power density fraction in fuel.
* ACOOL   coolant cross section area in m^2.
* HD      hydraulic diameter of one assembly in m.
* PCH     heating perimeter in m.
* MAXITC  maximum number of flow iterations.
* MAXIT1  maximum number of conduction iterations.
* MAXITL  maximum number of center-pellet iterations.
* ERMAXT  convergence criterion for temperature in fuel pin in K.
* ERMAXC  convergence criterion for coolant parameters (relative error).
* SPDIN   inlet flow velocity at t in m/s.
* TINLET  inlet temperature at t in K.
* POULET  outlet pressure at t in Pa.
* FRACPU  plutonium fraction in fuel.
* ICONDF  fuel conductivity flag (0=Stora-Chenebault or COMETHE/
*         1=user-provided polynomial + inverse term).
* NCONDF  degree of user-provided fuel conductivity polynomial.
* KCONDF  polynomial coefficients for fuel conductivity in W/m/K^(k+1)
*         (except for the two last coefficients which belongs to the
*         inverse term).
* UCONDF  required unit of temperature in polynomial for fuel
*         conductivity (KELVIN or CELSIUS).
* ICONDC  clad conductivity flag (0=default/1=user-provided
*         polynomial).
* NCONDC  degree of user-provided clad conductivity polynomial.
* KCONDC  polynomial coefficients for clad conductivity in W/m/K^(k+1).
* UCONDC  required unit of temperature in polynomial for clad
*         conductivity (KELVIN or CELSIUS).
* IHGAP   flag indicating HGAP chosen (0=default/1=user-provided).
* KHGAP   fixed user-provided HGAP value in W/m^2/K.
* IHCONV  flag indicating HCONV chosen (0=default/1=user-provided).
* KHCONV  fixed user-provided HCONV value in W/m^2/K.
* WTEFF   surface temperature's weighting factor in effective fuel
*         temperature.
* IFRCDI  flag indicating if average approximation is forced during
*         fuel conductivity evaluation (0=default/1=average
*         approximation forced).
* ISUBM   subcooling model (0: one-phase; 1: Bowring model; 2: Saha-
*         Zuber model).
* FRO     radial power form factors.
* POW     power distribution at t in W.
*
*Parameters: output
* TCOMB   averaged fuel temperature distribution in K.
* DCOOL   averaged coolant density distribution in g/cc.
* TCOOL   averaged coolant temperature distribution in K.
* TSURF   surface fuel temperature distribution in K.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) MPTHMI,MPTHM
      INTEGER IMPX,IX,IY,NZ,NFD,NDTOT,IFLUID,MAXITC,MAXIT1,MAXITL,IHGAP
      REAL XBURN(NZ),VOLXY,HZ(NZ),DTIME,CFLUX,POROS,FNFU,FFUEL,ACOOL,
     > HD,PCH,ERMAXT,ERMAXC,FCOOL,SPDIN,TINLET,POULET,FRACPU,
     > KCONDF(NCONDF+3),KCONDC(NCONDC+1),KHGAP,KHCONV,WTEFF,FRO(NFD-1),
     > POW(NZ),TCOMB(NZ),DCOOL(NZ),TCOOL(NZ),TSURF(NZ)
      CHARACTER UCONDF*12,UCONDC*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER(KMAXO=100,MAXNPO=40,PES=9.81)
      REAL ENT(4),RHOINL,MFLXIN,RHOIN0,MFLXIN0,HINLET,HINLE0,MUIN,
     > DV(NZ),PARAM1,PARAM2,PARAM3,ERRG,ERRP,ERRH,ERR,DELTH,HMINF,
     > POWLIN(NZ),PHI(NZ),MUT(NZ),RESM(NZ),RESP(NZ),RESH(NZ),QFUEL(NZ),
     > QCOOL(NZ),TC1,AGM(NZ),PC(NZ),TSAT,PHIC(NZ),TP(NZ),TLC(NZ),
     > HZC(NZ),XFL(NZ),EPS(NZ),TB,HGSAT,TCLAD(NZ),MFLXT0(NZ),ENTH(NZ),
     > MFLXT(NZ),SLIP(NZ),K11
      INTEGER KWA(NZ)
      REAL TRE10(MAXNPO),TRE11(MAXNPO),RADD(MAXNPO),XX2(MAXNPO),
     > XX3(MAXNPO),ZF(2)
      CHARACTER HSMG*131
      REAL XS(4)
      DATA XS/-0.861136,-0.339981,0.339981,0.861136/
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: VELOT0,DCOOL0,PREST0,ENTHT0,
     > DLIQT0,VELOT,PREST,ENTHT,TCENTT,DLIQT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RAD,TEMPT0,TEMPT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(RAD(NDTOT-1,NZ),VELOT0(NZ),DCOOL0(NZ),PREST0(NZ),
     > ENTHT0(NZ),TEMPT0(NDTOT,NZ),DLIQT0(NZ),VELOT(NZ),PREST(NZ),
     > ENTHT(NZ),TEMPT(NDTOT,NZ),TCENTT(NZ),DLIQT(NZ))
*----
*  RECOVER DATA FROM FORMER TIME STEP OR STEADY-STATE CALCULATION IN THM
*----
      CALL LCMGET(MPTHMI,'DENSITY',DCOOL0)
      CALL LCMGET(MPTHMI,'PRESSURE',PREST0)
      CALL LCMGET(MPTHMI,'ENTHALPY',ENTHT0)
      CALL LCMGET(MPTHMI,'VELOCITIES',VELOT0)
      CALL LCMGET(MPTHMI,'TEMPERATURES',TEMPT0)
      CALL LCMGET(MPTHMI,'LIQUID-DENS',DLIQT0)
      CALL LCMGET(MPTHMI,'POULET',POUT0)
      CALL LCMGET(MPTHMI,'TINLET',TIN0)
      CALL LCMGET(MPTHMI,'RADII',RAD)
*----
*  CALCULATE THE INVERSE TIME STEP
*----
      IF(DTIME.EQ.0.0) THEN
        CALL XABORT('THMTRS: TIME STEP NOT DEFINED')
      ELSE
        DTINV=1.0/DTIME
      ENDIF
*----
*  COMPUTE THE INLET FLOW ENTHALPY AND MASS FLOW RATE
*----
      IF(IFLUID.EQ.0) THEN
         CALL THMSAT(POULET,TSAT)
      ELSE IF(IFLUID.EQ.1) THEN
         CALL THMHST(POULET,TSAT)
      ENDIF
      IF(TINLET.GT.TSAT) THEN
         WRITE(HSMG,'(28HTHMTRS: OUTLET TEMPERATURE (,1P,E12.4,
     1   40H K) GREATER THAN SATURATION TEMPERATURE.)') TINLET
         CALL XABORT(HSMG)
      ENDIF
      IF(IFLUID.EQ.0) THEN
        CALL THMPT(POUT0,TIN0,RHOIN0,HINLE0,R3,R4,R5)
      ELSE IF(IFLUID.EQ.1) THEN
        CALL THMHPT(POUT0,TIN0,RHOIN0,HINLE0,R3,R4,R5)
      ENDIF
      MFLXIN0=SPDIN*RHOIN0
      CALL THMPT(POULET,TINLET,RHOINL,HINLET,R3,MUIN,CPVIN)
      MFLXIN=SPDIN*RHOINL
      IF(NDTOT.GT.MAXNPO) CALL XABORT('THMTRS: MAXNPO OVERFLOW')
*----
*  MAIN LOOP ALONG THE 1D CHANNEL.
      DO K=1,NZ
*----
*  COMPUTE THE LINEAR POWER, THE VOLUMIC POWER, THE THERMAL EXCHANGE
*  COEFFICIENT OF THE GAP AND THE THERMAL HEAT FLUX ALONG THE CHANNEL
*----
         DV(K)=VOLXY*HZ(K)
*        linear power in W/m.
         POWLIN(K)=(POW(K)/DV(K))*VOLXY/FNFU
*        volumic power in W/m^3.
         QFUEL(K)=POW(K)*FFUEL/DV(K)
         QCOOL(K)=POW(K)*FCOOL/DV(K)
*----
*  INITIALIZATION OF THE THERMO-HYDRAULICAL PROPERTIES OF THE FLUID
*----
         DCOOL(K)=DCOOL0(K)
         MUT(K)=MUIN
         VELOT(K)=VELOT0(K)
         MFLXT0(K)=DCOOL0(K)*VELOT(K)
         MFLXT(K)=MFLXT0(K)
         PREST(K)=PREST0(K)
         ENTHT(K)=ENTHT0(K)
         DLIQT(K)=DLIQT0(K)
         DO L=1,NDTOT
           TEMPT(L,K)=TEMPT0(L,K)
         ENDDO
         RESM(K)=MFLXT(K)
         RESP(K)=PREST(K)
         RESH(K)=ENTHT(K)
      ENDDO
*----
*  ITERATIVE PROCEDURE FOR EACH CHANNEL
*----
      DO K=1,NZ
         XFL(K)=0.0
         EPS(K)=0.0
         XFL(K)=0.0
         MFLXT(K)=0.0
         SLIP(K)=1.0
         KWA(K)=0
      ENDDO
      KMIN=1
      DO K=1,NZ
        IF(POW(K).NE.0.0) THEN
          KMIN=K
          EXIT
        ENDIF
      ENDDO
      ITERC=0
   20 ITERC=ITERC+1
      IF(ITERC.GT.MAXITC) THEN
         CALL XABORT('THMTRS: CONVERGENCE FAILURE IN FLOW CALCULATION.')
      ENDIF
*----
*  MAIN LOOP ALONG THE 1D CHANNEL.
*----
      K0=0 ! onset of nuclear boiling point
      DO K=KMIN,NZ
        IF(POW(K).EQ.0.0) CYCLE
        IF(IMPX.GT.4) WRITE(6,190) K
*----
*  SOLVE THE CONDUCTION EQUATIONS INSIDE THE FUEL ROD
*----
        DO L=1,NDTOT-1
          TRE10(L)=TEMPT0(L,K)
          TRE11(L)=TEMPT(L,K)
          RADD(L)=RAD(L,K)
        ENDDO
        TSCLAD=TEMPT(NDTOT,K)
        CALL THMROD(IMPX,NFD,NDTOT-1,MAXIT1,MAXITL,ERMAXT,DTINV,
     1  RADD,TRE10,TRE11,QFUEL(K),FRO,TSCLAD,POWLIN(K),XBURN(K),
     2  POROS,FRACPU,ICONDF,NCONDF,KCONDF,UCONDF,ICONDC,NCONDC,
     3  KCONDC,UCONDC,IHGAP,KHGAP,IFRCDI,TC1,XX2,XX3,ZF)
*----
*  COMPUTE THE HEAT FLUX FROM CLAD TO COOLANT IN W/m^2
*----
        PHI(K)=(ZF(1)-TSCLAD*ZF(2))/RAD(NDTOT-1,K)
        IF(PHI(K).GT.CFLUX) THEN
          WRITE(HSMG,'(23HTHMTRS: THE HEAT FLUX (,1P,E12.4,5H) IS ,
     >    37HGREATER THAN THE CRITICAL HEAT FLUX (,E12.4,2H).)')
     >    PHI(K),CFLUX
          WRITE(6,'(/1X,A)') HSMG
        ENDIF
*----
*  FLOW RATE CALCULATION WITH MASS CONSERVATION EQUATION
*----
        PARAM1=0.5*(DCOOL0(K)-DCOOL(K))*DTINV*HZ(K)
        IF(K.EQ.KMIN) THEN
          PARAM1=PARAM1+0.5*(RHOIN0-RHOINL)*DTINV*HZ(K)
          MFLXT(K)=MFLXIN+PARAM1
        ELSE
          PARAM1=PARAM1+0.5*(DCOOL0(K-1)-DCOOL(K-1))*DTINV*HZ(K)
          MFLXT(K)=MFLXT(K-1)+PARAM1
        ENDIF
*----
*  ENTHALPY VECTOR CALCULATION WITH ENERGY CONSERVATION EQUATION
*----
        PARAM1=0.5*DCOOL(K)*DTINV*HZ(K)+MFLXT(K)
        PARAM2=0.5*DCOOL0(K)*ENTHT0(K)*DTINV*HZ(K)
        PARAM3=(QCOOL(K)+PHI(K)*PCH/ACOOL)*HZ(K)
        IF(K.EQ.KMIN) THEN
          PARAM2=PARAM2+0.5*(RHOIN0*HINLE0-RHOINL*HINLET)*DTINV*HZ(K)
          PARAM2=PARAM2+MFLXIN*HINLET
          HMINF=HINLET
        ELSE
          PARAM2=PARAM2+0.5*(DCOOL0(K-1)*ENTHT0(K-1)-
     1    DCOOL(K-1)*ENTHT(K-1))*DTINV*HZ(K)
          PARAM2=PARAM2+MFLXT(K-1)*ENTHT(K-1)
          HMINF=ENTHT(K-1)
        ENDIF
        ENTHT(K)=(PARAM2+PARAM3)/PARAM1
        DELTH=ENTHT(K)-HMINF
*----
*  COMPUTE THE COOLANT TEMPERATURE AND THE OUTER CLADDING TEMPERATURE
*----
        DO I1=1,4
          POINT=(1.0+XS(I1))/2.0
          ENT(I1)=HMINF+POINT*DELTH
        ENDDO
        IF(K.GT.1) THEN
          XFL(K)=XFL(K-1)
          EPS(K)=EPS(K-1)
          SLIP(K)=SLIP(K-1)
        ENDIF
        CALL THMH2O(1,IX,IY,K,K0,PREST(K),MFLXT(K),ENTHT(K),ENT,HD,
     >  IFLUID,IHCONV,KHCONV,ISUBM,RAD(NDTOT-1,K),ZF,PHI(K),XFL(K),
     >  EPS(K),SLIP(K),ACOOL,PCH,HZ(K),TCALO,DCOOL(K),DLIQT(K),
     >  TRE11(NDTOT),KWA(K))
*
        DO L=1,NDTOT-1
          TRE11(L)=XX2(L)+TRE11(NDTOT)*XX3(L)
          TEMPT(L,K)=TRE11(L)
        ENDDO
        TEMPT(NDTOT,K)=TRE11(NDTOT)
*----
*  RECOVER MESHWISE TEMPERATURES AND FLUID DENSITY. BY DEFAULT, USE THE
*  ROWLANDS FORMULA TO COMPUTE THE EFFECTIVE FUEL TEMPERATURE, OTHERWISE
*  USE USER-SPECIFIED WEIGHTING FACTOR.
*----
        TCOMB(K)=(1.0-WTEFF)*TC1+WTEFF*TRE11(NFD)
        TCOOL(K)=TCALO
        TCENTT(K)=TC1
        TSURF(K)=TRE11(NFD)
        TCLAD(K)=TRE11(NDTOT)
      ENDDO
*----
*  MOMENTUM VECTOR CALCULATION WITH MOMENTUM CONSERVATION EQUATION
*----
*      DO K=NZ,1,-1
*        IF(POW(K).EQ.0.0) CYCLE
*        RET=ABS(MFLXT(K))*(1.0-XFL(K))*HD/MUT(K)
*        PARAM1=0.5*(MFLXT(K)-MFLXT0(K))*DTINV*HZ(K)
*        PARAM2=MFLXT(K)**2.0/DCOOL(K)
*        CALL THMFRI(RET,F)
*        IF(XFL(K).GT.0.0) THEN
*          CALL THMPLO(PREST(K),XFL(K),PHIL0)
*        ELSE
*          PHIL0=1.0
*        ENDIF
*        PARAM31=DCOOL(K)*PES
*        PARAM32=0.5*F*MFLXT(K)**2.0/HD/DLIQT0(K)*PHIL0
*        PARAM3=(PARAM31+PARAM32)*HZ(K)
*        IF(K.EQ.1) THEN
*          PARAM1=PARAM1+0.5*(MFLXIN-MFLXIN0)*DTINV*HZ(1)
*          PARAM2=PARAM2-MFLXIN**2.0/RHOINL
*          PREST(1)=PREST(2)+PARAM1+PARAM2+PARAM3
*        ELSE IF(K.LT.NZ) THEN
*          PARAM1=PARAM1+0.5*(MFLXT(K-1)-MFLXT0(K-1))*DTINV*
*     1    HZ(K)
*          PARAM2=PARAM2-MFLXT(K-1)**2.0/DCOOL(K-1)
*          PREST(K)=PREST(K+1)+PARAM1+PARAM2+PARAM3
*        ELSE IF(K.EQ.NZ) THEN
*          PARAM1=PARAM1+0.5*(MFLXT(NZ-1)-MFLXT0(NZ-1))*DTINV*
*     1    HZ(NZ)
*          PARAM2=PARAM2-MFLXT(K-1)**2.0/DCOOL(K-1)
*          PREST(NZ)=POULET+PARAM1+PARAM2+PARAM3
*        ENDIF
*      ENDDO
      PINLET=PREST(KMIN)
*----
*  CALCULATE THE VOID FRACTION COEFFICIENT AND THE STEAM QUALITY
*----
      DO K=1,NZ
        HZC(K)=HZ(K)
        PHIC(K)=PHI(K)
        TP(K)=TCLAD(K)
        TLC(K)=TCOOL(K)
        ENTH(K)=ENTHT(K)    
        AGM(K)=MFLXT(K)
        PC(K)=PREST(K)
      ENDDO
*----
*  COMPUTE NEW VALUES OF DENSITIES AND VELOCITIES OVER CHANNEL
*----
       DO K=1,NZ
         IF(EPS(K).GT.0.0) THEN
           IF(IFLUID.EQ.0) THEN
             CALL THMSAT(PREST(K),TSAT)
             CALL THMTX(TSAT,1.0,RGSAT,HGSAT,R3,R4,R5)
           ELSE IF(IFLUID.EQ.1) THEN
             CALL THMHST(PREST(K),TSAT)
             CALL THMHTX(TSAT,1.0,RGSAT,HGSAT,R3,R4,R5)
           ENDIF
           DCOOL(K)=DLIQT(K)*(1.0-EPS(K))+EPS(K)*RGSAT
         ELSE
           DCOOL(K)=DLIQT(K)
         ENDIF
         VELOT(K)=MFLXT(K)/DCOOL(K)
       ENDDO
*----
*  CONVERGENCE TEST FOR THE ENTHALPY, PRESSURE DENSITY AND
*  MASS FLUX CALCULATION.
*----
      ERRG=0.0
      ERRP=0.0
      ERRH=0.0
      ERR=0.0
      ERX=0.0
      DO K=1,NZ
        IF(POW(K).EQ.0.0) CYCLE
        CALL THMSAT(PREST(K),TSAT)
        TB=TSAT-0.1
        IF(TCOOL(K).LT.TB) THEN
           CALL THMPT(PREST(K),TCOOL(K),R11,H11,K11,MUT(K),C11)
        ELSE
           CALL THMPT(PREST(K),TB,R11,H11,K11,MUT(K),C11)
        ENDIF
        ERRG=MAX(ERRG,ABS(MFLXT(K)-RESM(K))/MFLXT(K))
        ERRP=MAX(ERRP,ABS(PREST(K)-RESP(K))/PREST(K))
        ERRH=MAX(ERRH,ABS(ENTHT(K)-RESH(K))/ENTHT(K))
        RESM(K)=MFLXT(K)
        RESP(K)=PREST(K)
        RESH(K)=ENTHT(K)
      ENDDO
      ERR=MAX(ERRG,ERRP,ERRH)
      IF(IMPX.GT.1) WRITE(6,200) ITERC,ERRG,ERRP,ERRH
      CALL THMPT(PINLET,TINLET,RHOINL,HINLET,R3,MUIN,CPVIN)
      IF((ERR.LT.ERMAXC).AND.(ITERC.GT.1)) THEN
        GO TO 30
      ELSE
        GO TO 20
      ENDIF
*----
*  PRINT THE OUTLET THERMOHYDRAULICAL PARAMETERS
*----
   30 IF(IMPX.GT.3) THEN
         WRITE(6,'(/16H THMTRS: CHANNEL,2I6/1X,27(1H-))') IX,IY
         WRITE(6,210) ' ___________________________________________',
     >          '____________________________________________________',
     >          '____________________________________________________',
     >          '_______________________________'
         WRITE(6,210) '|     |   TFUEL    |   TSURF    |    MFLXT  ',
     >          '  |    DCOOL    |    TCOOL    |    PCOOL    |    HCO',
     >          'OL    |    QFUEL    |    QCOOL    |    VOID   |     ',
     >          'QUAL    |     SLIP    |  FLOW  |',
     >          '|     |     K      |     K      |   Kg/m2/s   |    K',
     >          'g/m3    |      K      |     Pa      |    J/Kg     | ',
     >          '   W/m3     |    W/m3     |           |             ',
     >          '|             | REGIME |'
         WRITE(6,210) '|_____|____________|____________|___________',
     >          '__|_____________|_____________|_____________|_______',
     >          '______|_____________|_____________|___________|_____',
     >          '________|_____________|________|'
         DO L=NZ,1,-1
           IF(L.EQ.1) THEN
             WRITE(6,220) '| BOT |',TCOMB(L),' |',TSURF(L),
     >            ' |',MFLXT(L),' |',DCOOL(L),' |',TCOOL(L),
     >            ' |',PREST(L),' |',ENTHT(L),' |',QFUEL(L),
     >            ' |',QCOOL(L),' |',EPS(L),' |',XFL(L),' |',SLIP(L),
     >            ' |',KWA(L),' |'
           ELSEIF(L.EQ.NZ) THEN
             WRITE(6,220) '| TOP |',TCOMB(L),' |',TSURF(L),
     >            ' |',MFLXT(L),' |',DCOOL(L),' |',TCOOL(L),
     >            ' |',PREST(L),' |',ENTHT(L),' |',QFUEL(L),
     >            ' |',QCOOL(L),' |',EPS(L),' |',XFL(L),' |',SLIP(L),
     >            ' |',KWA(L),' |'
           ELSE
             WRITE(6,225) '| ',L,' |',TCOMB(L),' |',TSURF(L),
     >            ' |',MFLXT(L),' |',DCOOL(L),' |',TCOOL(L),
     >            ' |',PREST(L),' |',ENTHT(L),' |',QFUEL(L),
     >            ' |',QCOOL(L),' |',EPS(L),' |',XFL(L),' |',SLIP(L),
     >            ' |',KWA(L),' |'
           ENDIF
         ENDDO
         WRITE(6,210) '|_____|____________|____________|___________',
     >          '__|_____________|_____________|_____________|_______',
     >          '______|_____________|_____________|___________|_____',
     >          '________|_____________|________|'

      ENDIF
*----
*  MODIFICATION OF THE VECTORS TO FIT THE GEOMETRY OF THE CHANNELS AND
*  THE BUNDLES AND WRITE THE DATA IN LCM OBJECT THM 
*----
      CALL LCMPUT(MPTHM,'PRESSURE',NZ,2,PREST)
      CALL LCMPUT(MPTHM,'DENSITY',NZ,2,DCOOL)
      CALL LCMPUT(MPTHM,'ENTHALPY',NZ,2,ENTHT)
      CALL LCMPUT(MPTHM,'VELOCITIES',NZ,2,VELOT)
      CALL LCMPUT(MPTHM,'CENTER-TEMPS',NZ,2,TCENTT)
      CALL LCMPUT(MPTHM,'COOLANT-TEMP',NZ,2,TCOOL)
      CALL LCMPUT(MPTHM,'LIQUID-DENS',NZ,2,DLIQT)
      CALL LCMPUT(MPTHM,'PINLET',1,2,PINLET)
      CALL LCMPUT(MPTHM,'TINLET',1,2,TINLET)
      CALL LCMPUT(MPTHM,'VINLET',1,2,SPEED)
      CALL LCMPUT(MPTHM,'POWER',NZ,2,POW)
      CALL LCMPUT(MPTHM,'POULET',1,2,POULET)
      CALL LCMPUT(MPTHM,'TEMPERATURES',NDTOT*NZ,2,TEMPT)
      CALL LCMPUT(MPTHM,'RADII',(NDTOT-1)*NZ,2,RAD)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DLIQT,TCENTT,TEMPT,ENTHT,PREST,VELOT,DLIQT0,TEMPT0,
     > ENTHT0,PREST0,DCOOL0,VELOT0,RAD)
      RETURN
*
  190 FORMAT(/21H THMTRS: AXIAL SLICE=,I5)
  200 FORMAT(/24H THMTRS: FLOW ITERATION=,I5,1P,8H  ERROR=,3E12.4)
  210 FORMAT(1X,A,A,A,A)
  220 FORMAT(1X,A,F11.2,A,F11.2,A,F12.4,A,F12.4,A,F12.2,A,3P,E12.4,
     >       A,1P,E12.4,A,1P,E12.4,A,1P,E12.4,A,0P,F10.4,A,E12.4,A,
     >       E12.4,A,I5,2X,A)
  225 FORMAT(1X,A,I3,A,F11.2,A,F11.2,A,F12.4,A,F12.4,A,F12.2,A,3P,
     >       E12.4,A,1P,E12.4,A,1P,E12.4,A,1P,E12.4,A,0P,F10.4,A,
     >       E12.4,A,E12.4,A,I5,2X,A)
      END
