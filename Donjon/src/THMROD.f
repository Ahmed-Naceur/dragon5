*DECK THMROD
      SUBROUTINE THMROD(IMPX,NFD,NDTOT,MAXIT1,MAXIT2,ERMAXT,DTINV,
     1 RAD,XX0,XX1,QFUEL,FRO,TSURF,POWLIN,BURN,POROS,FRACPU,ICONDF,
     2 NCONDF,KCONDF,UCONDF,ICONDC,NCONDC,KCONDC,UCONDC,IHGAP,KHGAP,
     3 IFRCDI,TC1,XX2,XX3,ZF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of the discretized thermal conduction equations in a single
* fuel rod in an axial slice of the fuel channel.
*
*Copyright:
* Copyright (C) 2018 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IMPX    print parameter.
* NFD     number of fuel discretized points in the cladded fuel rod.
*         The last point of the discretization (i=NFD) is taken at the
*         surface of the fuel pellet.
* NDTOT   total number of discretized points in the cladded fuel rod
*         with the radial zones in the cladding. The points which are
*         located at i=NFD+1 and i=NDTOT are respectively taken at the
*         inner surface of clad and in the center of external clad ring.
* MAXIT1  maximum number of conduction iterations.
* MAXIT2  maximum number of center-pellet iterations.
* ERMAXT  convergence criterion.
* DTINV   inverse time step. Equal to 1/DT in transient cases. Equal to
*         0 in steady-state cases.
* RAD     fuel and clad radii (m).
* XX0     temperatures at time n-1 (K).
* XX1     estimate of the temperatures at time n (K).
* QFUEL   volumic power in fuel at time n (W/m^3).
* FRO     radial power form factors. All components are set to 1.0 for
*         a constant power source in fuel.
* TSURF   estimate of the external clad surface temperature at
*         time n (K).
* POWLIN  estimate of the linear power at time n (W/m).
* BURN    fuel burnup in MWday/tonne.
* POROS   fuel porosity.
* FRACPU  plutonium percent fraction.
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
* IFRCDI  flag indicating if average approximation is forced during
*         fuel conductivity evaluation (0=default/1=average
*         approximation forced).
*
*Parameters: output
* TC1     estimate of center-pellet temperature at time n (K).
* XX1     estimate of the temperatures at time n (K).
* XX2     first component of temperatures at time n (K).
* XX3     second component of temperatures at time n. The actual
*         temperatures are given as XX2(:)+TSURF*XX3(:) where TSURF
*         is the temperature of the external clad surface.
* ZF      components of the linear power transmitted from clad to fluid.
*         The linear power (W/m) is given as 2*PI*(ZF(1)-TSURF*ZF(2)).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,NFD,NDTOT,MAXIT1,MAXIT2,ICONDF,NCONDF,ICONDC,NCONDC,
     1 IHGAP,IFRCDI
      REAL ERMAXT,DTINV,RAD(NDTOT),XX0(NDTOT),XX1(NDTOT),QFUEL,
     1 FRO(NFD-1),TSURF,POWLIN,BURN,POROS,FRACPU,KCONDF(NCONDF+3),
     2 KCONDC(NCONDC+1),KHGAP,TC1,XX2(NDTOT),XX3(NDTOT),ZF(2)
      CHARACTER UCONDF*12,UCONDC*12
*----
*  LOCAL VARIABLES
*----
      REAL COEF(3)
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: DAR,ZK,CONDXA
      REAL, ALLOCATABLE, DIMENSION(:,:) :: TRID
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(DAR(NDTOT),ZK(NDTOT),CONDXA(NDTOT),TRID(NDTOT,NDTOT+2))
*----
*  COMPUTE ARs AND VOLUMES
*----
      CALL XDRSET(DAR,NDTOT,0.0)
      ARF=0.5*RAD(NFD)**2    ! at fuel radius
      ARCI=0.5*RAD(NFD+1)**2 ! at internal clad radius
      ARCE=0.5*RAD(NDTOT)**2 ! at external clad radius
      DO I=1,NFD-1
        DAR(I)=0.5*(RAD(I+1)**2-RAD(I)**2)
      ENDDO
      DO I=NFD+2,NDTOT
        DAR(I)=0.5*(RAD(I)**2-RAD(I-1)**2)
      ENDDO
*----
*  COMPUTE THE THERMAL CONDUCTIVITY INTEGRALS AT TIME n-1
*----
      CALL XDRSET(ZK,NDTOT,0.0)
      DO I=1,NFD-1
        ZK(I)=THMCDI(XX0(I),XX0(I+1),BURN,POROS,FRACPU,ICONDF,NCONDF,
     >  KCONDF,UCONDF,IFRCDI)
      ENDDO
      DO I=NFD+1,NDTOT-1
        ZK(I)=THMGDI(XX0(I),XX0(I+1),ICONDC,NCONDC,KCONDC,UCONDC)
      ENDDO
      ZK(NDTOT)=THMGDI(XX0(NDTOT),TSURF,ICONDC,NCONDC,KCONDC,UCONDC)
*----
*  COMPUTE CONDXA
*----
      CALL XDRSET(CONDXA,NDTOT,0.0)
      COEF(1)=0.0
      COEF(2)=0.0
      COEF(3)=0.0
      DO I=1,NDTOT
        IF(I.LE.NFD-1) THEN
          CONDXA(I)=THMCCD(XX0(I),POROS,FRACPU)*DTINV*XX0(I)*DAR(I)
        ELSE IF(I.GE.NFD+2) THEN
          CONDXA(I)=THMGCD(XX0(I))*DTINV*XX0(I)*DAR(I)
        ENDIF
      ENDDO
*----
*  ITERATIVE PROCEDURE
*----
      ITERT=0
   10 ITERT=ITERT+1
      IF(ITERT.GT.MAXIT1) CALL XABORT('THMROD: CONVERGENCE FAILURE(1).')
*----
*  COMPUTE THE THERMAL CONDUCTIVITY INTEGRALS AT TIME n
*----
      CALL XDRSET(ZK,NDTOT,0.0)
      DO I=1,NFD-1
        ZK(I)=THMCDI(XX1(I),XX1(I+1),BURN,POROS,FRACPU,ICONDF,NCONDF,
     >  KCONDF,UCONDF,IFRCDI)
      ENDDO
      DO I=NFD+1,NDTOT-1
        ZK(I)=THMGDI(XX1(I),XX1(I+1),ICONDC,NCONDC,KCONDC,UCONDC)
      ENDDO
      ZK(NDTOT)=THMGDI(XX1(NDTOT),TSURF,ICONDC,NCONDC,KCONDC,UCONDC)
      IF(IHGAP.EQ.0) THEN
         CALL THMGAP(POWLIN,BURN,HGAP)
      ELSE IF(IHGAP.EQ.1) THEN
         HGAP=KHGAP
      ENDIF
*----
*  BUILD THE TRIDIAGONAL SYSTEM
*----
      CALL XDRSET(TRID,NDTOT*(NDTOT+2),0.0)
      COEF(1)=0.0
      COEF(2)=0.0
      COEF(3)=0.0
      DO I=1,NDTOT
        TRID(I,NDTOT+1)=CONDXA(I)
        IF(I.LE.NFD-2) THEN
          ARI=0.5*RAD(I+1)**2
          COEF(3)=4.0*ARI*ZK(I)/(DAR(I)+DAR(I+1))
          TRID(I,NDTOT+1)=TRID(I,NDTOT+1)+QFUEL*FRO(I)*DAR(I)
        ELSE IF(I.EQ.NFD-1) THEN
          ARI=0.5*RAD(I+1)**2
          COEF(3)=4.0*ARI*ZK(I)/DAR(I)
          TRID(I,NDTOT+1)=TRID(I,NDTOT+1)+QFUEL*FRO(I)*DAR(I)
        ELSE IF(I.EQ.NFD) THEN
          RAVG=2.0*RAD(NFD)*RAD(NFD+1)/(RAD(NFD)+RAD(NFD+1))
          COEF(3)=RAVG*HGAP
        ELSE IF(I.EQ.NFD+1) THEN
          COEF(3)=4.0*ARCI*ZK(I)/DAR(I+1)
        ELSE IF(I.LE.NDTOT-1) THEN
          ARI=0.5*RAD(I)**2
          COEF(3)=4.0*ARI*ZK(I)/(DAR(I)+DAR(I+1))
        ELSE IF(I.EQ.NDTOT) THEN
          COEF(3)=4.0*ARCE*ZK(I)/DAR(I)
          TRID(I,NDTOT+2)=TRID(I,NDTOT+2)+COEF(3)
        ENDIF
        COEF(2)=COEF(1)+COEF(3)
        IF(I.GT.1) TRID(I,I-1)=-COEF(1)
        IF(I.LE.NFD-1) THEN
          TRID(I,I)=THMCCD(XX1(I),POROS,FRACPU)*DTINV*DAR(I)
        ELSE IF(I.GE.NFD+2) THEN
          TRID(I,I)=THMGCD(XX1(I))*DTINV*DAR(I)
        ENDIF
        TRID(I,I)=TRID(I,I)+COEF(2)
        IF(I.LT.NDTOT) THEN
          TRID(I,I+1)=-COEF(3)
          COEF(1)=COEF(3)
        ENDIF
      ENDDO
      ZWORK=COEF(3)
*----
*  SOLVE LINEAR SYSTEM
*----
      CALL ALSB(NDTOT,2,TRID,IER,NDTOT)
      IF(IER.NE.0) CALL XABORT('THMROD: SINGULAR MATRIX')
*----
*  SET TEMPERATURE AT TIME n
*----
      ERR=0.0
      IMAX=0
      DO I=1,NDTOT
        TNEW=TRID(I,NDTOT+1)+TSURF*TRID(I,NDTOT+2)
        IF(ABS(XX1(I)-TNEW).GT.ERR) THEN
          ERR=ABS(XX1(I)-TNEW)
          IMAX=I
        ENDIF
        IF(ITERT.LE.20) THEN
          XX1(I)=TNEW
        ELSE
*         perform under-relaxation
          XX1(I)=0.5*(TNEW+XX1(I))
        ENDIF
      ENDDO
      ZF(1)=ZWORK*TRID(NDTOT,NDTOT+1)
      ZF(2)=ZWORK*(1.0-TRID(NDTOT,NDTOT+2))
      IF(IMPX.GT.4) WRITE(6,100) ITERT,ERR,ERMAXT,IMAX
      IF((ERR.LT.ERMAXT).AND.(ITERT.NE.1)) GO TO 20
      GO TO 10
*----
*  SCRATCH STORAGE DEALLOCATION
*----
   20 DO I=1,NDTOT
        XX2(I)=TRID(I,NDTOT+1)
        XX3(I)=TRID(I,NDTOT+2)
      ENDDO
      DEALLOCATE(TRID,CONDXA,ZK,DAR)
*----
*  COMPUTE THE CENTER-PELLET TEMPERATURE.
*----
      TC=0.5*(XX1(1)+XX1(2))
      ITERC=0
   30 ITERC=ITERC+1
      IF(ITERC.GT.MAXIT2) CALL XABORT('THMROD: CONVERGENCE FAILURE(2).')
      TCOLD=TC
      CC1=THMCDI(XX1(1),TC,BURN,POROS,FRACPU,ICONDF,NCONDF,KCONDF,
     > UCONDF,IFRCDI)
      CC2=THMCDI(TC,XX1(2),BURN,POROS,FRACPU,ICONDF,NCONDF,KCONDF,
     > UCONDF,IFRCDI)
      TC=(CC1*XX1(1)+CC2*XX1(2))/(CC1+CC2)
      IF(ITERT.GT.20) TC=0.5*(TC+TCOLD)
      DELTAA=ABS(TC-TCOLD)
      IF(IMPX.GT.4) WRITE(6,110) ITERC,DELTAA,ERMAXT
      IF((DELTAA.LT.ERMAXT).AND.(ITERC.NE.1)) GO TO 40
      GO TO 30
   40 TC1=2.0*XX1(1)-TC
      RETURN
  100 FORMAT(/15H THMROD: ITERT=,I5,1P,7H ERROR=,E12.4,5H EPS=,E12.4,
     > 5H POS=,I5)
  110 FORMAT(/15H THMROD: ITERC=,I5,1P,7H ERROR=,E12.4,5H EPS=,E12.4)
      END
