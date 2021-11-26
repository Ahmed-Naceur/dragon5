*DECK EVOSIG
      SUBROUTINE EVOSIG(IMPX,INR,IGLOB,NGROUP,NBMIX,NBISO,NCOMB,
     1 ISONAM,IPISO,DEN,FLUMIX,VX,MILVO,JM,NVAR,NSUPS,NREAC,HREAC,
     2 IDR,RER,RRD,FIT,FUELDN,NXSPER,DELTAT,MIXPWR,PFACT,SIG,VPHV)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute and normalize the microscopic depletion reaction rates.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IMPX    print flag (equal to zero for no print).
* INR     type of flux normalization:
*         =0: out-of-core depletion;
*         =1: constant flux depletion;
*         =2: constant fuel power depletion;
*         =3: constant assembly power depletion.
* IGLOB   out-of-fuel power in flux normalization:
*         =0: compute the burnup using the power released in the fuel;
*         =1: compute the burnup using the power released in the global
*         geometry.
* NGROUP  number of energy groups.
* NBMIX   number of mixtures.
* NBISO   number of isotopes/materials including non-depleting ones.
* NCOMB   number of depleting mixtures.
* ISONAM  alias name of isotopes.
* IPISO   pointer array towards microlib isotopes.
* DEN     density of each isotope.
* FLUMIX  average fluxes in mixtures.
* VX      volumes of the depleting mixtures.
* MILVO   mixture index corresponding to each depleting mixture.
* JM      position in isotope list of each nuclide of the depletion
*         chain. A negative value indicates a non-depleting isotope
*         producing energy.
* NVAR    number of depleting nuclides.
* NSUPS   number of non-depleting isotopes producing energy.
* NREAC   maximum number of depletion reactions.
* HREAC   names of used depletion reactions:
*         HREAC(1)='DECAY'; HREAC(2)='NFTOT';
*         HREAC(3)='NG'   ; HREAC(4)='N2N';  etc.
* IDR     identifier for each depleting reaction.
* RER     energy (Mev) per reaction. If RER(3,J)=0., the fission energy
*         is including radiative capture energy. Neutrino energy is
*         never included.
* RRD     sum of radioactive decay constants in 10**-8/s.
* FIT     flux normalization factor:
*         n/cm**2/s if INR=1;
*         MW/tonne of initial heavy elements if INR=2;
*         W/cc of assembly volume if INR=3.
* FUELDN  fuel initial density and mass.
* NXSPER  perturbation order for cross sections.
* DELTAT  perturbation coefficients for cross sections.
* MIXPWR  flags for mixtures to include in power normalization.
*
*Parameters: output
* PFACT   form factor for out-of-fuel power production.
* SIG     microscopic reaction rates for nuclide I in mixture IBM:
*         SIG(I,1,IBM) fission reaction rate;
*         SIG(I,2,IBM) gamma reaction rate;
*         SIG(I,3,IBM) N2N reaction rate;
*         ...;
*         SIG(I,NREAC,IBM) neutron-induced energy released;
*         SIG(I,NREAC+1,IBM) decay energy released (10**-8 MeV/s).
* VPHV    integrated fluxes in mixtures.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPISO(NBISO)
      INTEGER IMPX,INR,IGLOB,NGROUP,NBMIX,NBISO,NCOMB,ISONAM(3,NBISO),
     1 MILVO(NCOMB),JM(NBMIX,NVAR+NSUPS),NVAR,NSUPS,NREAC,
     2 HREAC(2,NREAC),IDR(NREAC,NVAR+NSUPS),NXSPER,MIXPWR(NBMIX)
      REAL DEN(NBISO),VX(NBMIX),RER(NREAC,NVAR+NSUPS),RRD(NVAR+NSUPS),
     1 FIT,FUELDN(3),DELTAT(2),PFACT,SIG(NVAR+1,NREAC+1,NBMIX),
     2 VPHV(NBMIX),FLUMIX(NGROUP,NBMIX)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6,MAXREA=20)
      TYPE(C_PTR) KPLIB
      CHARACTER HSMG*131,NAMDXS(MAXREA)*6
      DOUBLE PRECISION GAR,GAR1,GAR2,GARD,XDRCST,EVJ,FITD,PHI,FNORM,VPH
      INTEGER IPRLOC
      REAL, ALLOCATABLE, DIMENSION(:,:) :: XSREC
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(XSREC(NGROUP,NREAC-1))
*----
*  COMPUTE MICRO RATES
*----
      IPRLOC=0
      EVJ=XDRCST('eV','J')*1.0E22
      VPH=0.0
      CALL XDRSET(VPHV,NBMIX,0.0)
      DO 55 IU=1,NGROUP
      DO 40 IBM=1,NBMIX
      VPHV(IBM)=VPHV(IBM)+VX(IBM)*FLUMIX(IU,IBM)
   40 CONTINUE
      DO 50 ICMB=1,NCOMB
      IBM=MILVO(ICMB)
      IF(MIXPWR(IBM).EQ.1) VPH=VPH+VX(IBM)*FLUMIX(IU,IBM)
   50 CONTINUE
   55 CONTINUE
      DO 62 IBM=1,NBMIX
      DO 61 IQ=1,NREAC+1
      DO 60 IS=1,NVAR+1
      SIG(IS,IQ,IBM)=0.0
   60 CONTINUE
   61 CONTINUE
   62 CONTINUE
      IF(NREAC-1.GT.MAXREA) CALL XABORT('EVOSIG: MAXREA OVERFLOW.')
      DO 70 IREAC=2,NREAC
      WRITE(NAMDXS(IREAC-1),'(A4,A2)') HREAC(1,IREAC),HREAC(2,IREAC)
   70 CONTINUE
      DO 220 IBM=1,NBMIX
      IF(VX(IBM).EQ.0) GO TO 220
      DO 210 IST=1,NVAR+NSUPS
      K=JM(IBM,IST)
      IF(K.EQ.0) THEN
         GO TO 210
      ELSE IF(K.GT.0) THEN
*        DEPLETING ISOTOPE.
         IS=IST
         FACT=1.0
      ELSE
*        STABLE ISOTOPE PRODUCING ENERGY.
         K=-K
         IS=NVAR+1
         FACT=DEN(K)*VX(IBM)
      ENDIF
      SIG(IS,NREAC+1,IBM)=SIG(IS,NREAC+1,IBM)+FACT*RER(1,IST)*RRD(IST)
      IF(INR.EQ.0) GO TO 210
*----
*  RECOVER MULTIGROUP XS
*----
      KPLIB=IPISO(K) ! set K-th isotope
      IF(.NOT.C_ASSOCIATED(KPLIB)) THEN
        WRITE(HSMG,'(17HEVOSIG: ISOTOPE '',3A4,19H'' IS NOT AVAILABLE ,
     >  16HIN THE MICROLIB.)') (ISONAM(I0,K),I0=1,3)
        CALL XABORT(HSMG)
      ENDIF
      DO 150 IXSPER=1,NXSPER
      CALL XDRLXS(KPLIB,-1,IPRLOC,NREAC-1,NAMDXS,IXSPER,NGROUP,XSREC)
      DO 140 IREAC=2,NREAC
      CALL LCMLEN(KPLIB,NAMDXS(IREAC-1),LENGT,ITYLCM)
      IF((LENGT.NE.NGROUP).AND.(IDR(IREAC,IST).GT.0)) THEN
         IF((IREAC.EQ.2).AND.(MOD(IDR(2,IST),100).EQ.5)) GO TO 120
         IF(IMPX.GT.90) CALL LCMLIB(KPLIB)
         IF(IMPX.GT.3) THEN
           WRITE(HSMG,'(17HEVOSIG: REACTION ,A6,18H IS MISSING FOR IS,
     1     7HOTOPE '',3A4,2H''.)') NAMDXS(IREAC-1),(ISONAM(I0,K),I0=1,3)
           WRITE(IOUT,'(1X,A)') HSMG
         ENDIF
      ENDIF
  120 GAR=0.0D0
      DO 130 IU=1,NGROUP
      GAR=GAR+DBLE(XSREC(IU,IREAC-1)*FLUMIX(IU,IBM))
  130 CONTINUE
      SIG(IS,IREAC-1,IBM)=SIG(IS,IREAC-1,IBM)+1.0E-3*FACT*REAL(GAR)*
     1 DELTAT(IXSPER)
      SIG(IS,NREAC,IBM)=SIG(IS,NREAC,IBM)+1.0E-3*FACT*RER(IREAC,IST)*
     1 REAL(GAR)*DELTAT(IXSPER)
  140 CONTINUE
  150 CONTINUE
  210 CONTINUE
  220 CONTINUE
*----
*  CONSTANT FLUX OR CONSTANT POWER NORMALIZATION
*----
      PFACT=1.0
      PHI=0.0
      VTOT=0.0
      DO 230 ICMB=1,NCOMB
      IBM=MILVO(ICMB)
      IF(MIXPWR(IBM).EQ.1) VTOT=VTOT+VX(IBM)
  230 CONTINUE
      IF(INR.EQ.1) THEN
         PHI=FIT*1.E-13
      ELSE IF(INR.GE.2) THEN
         GAR=0.0D0
         GARD=0.0D0
         DO 245 ICMB=1,NCOMB
         IBM=MILVO(ICMB)
         IF(MIXPWR(IBM).EQ.1) THEN
            DO 240 IS=1,NVAR
            IF(JM(IBM,IS).GT.0) THEN
               GAR=GAR+VX(IBM)*DEN(JM(IBM,IS))*SIG(IS,NREAC,IBM)
               GARD=GARD+VX(IBM)*DEN(JM(IBM,IS))*SIG(IS,NREAC+1,IBM)
            ENDIF
  240       CONTINUE
         ENDIF
  245    CONTINUE
         GAR1=GAR
         DO 250 ICMB=1,NCOMB
         IBM=MILVO(ICMB)
         IF(MIXPWR(IBM).EQ.1) GAR1=GAR1+SIG(NVAR+1,NREAC,IBM)
  250    CONTINUE
         GAR2=GAR
         DO 260 IBM=1,NBMIX
         IF(MIXPWR(IBM).EQ.1) GAR2=GAR2+SIG(NVAR+1,NREAC,IBM)
  260    CONTINUE
         PFACT=REAL(GAR2/GAR1)
         IF((IGLOB.EQ.1).OR.(INR.EQ.3)) THEN
            GAR=GAR2
         ELSE IF(IGLOB.EQ.0) THEN
            GAR=GAR1
         ENDIF
         IF(GAR.EQ.0.0D0) CALL XABORT('EVOSIG: UNABLE TO NORMALIZE.')
         IF(INR.EQ.2) THEN
*           FITD IS THE DECAY POWER IN WATT PER GRAM.
            FITD=(EVJ*GARD)/(FUELDN(1)*VTOT)
            IF(FITD.GT.FIT) THEN
               WRITE(HSMG,'(35HEVOSIG: NEGATIVE FIT(1) FIT(DECAY)=,1P,
     1         E11.4,12H FIT(INPUT)=,E11.4,1H.)') FITD,FIT
               CALL XABORT(HSMG)
            ENDIF
            PHI=(FIT-FITD)*FUELDN(1)*VPH/(EVJ*GAR)
         ELSE IF(INR.EQ.3) THEN
*           FITD IS THE DECAY POWER IN WATT PER CUBIC CENTIMETER.
            FITD=(EVJ*GARD*FUELDN(3))/(FUELDN(1)*VTOT)
            IF(FITD.GT.FIT) THEN
               WRITE(HSMG,'(35HEVOSIG: NEGATIVE FIT(2) FIT(DECAY)=,1P,
     1         E11.4,12H FIT(INPUT)=,E11.4,1H.)') FITD,FIT
               CALL XABORT(HSMG)
            ENDIF
            PHI=(FIT-FITD)*FUELDN(1)*VPH/(EVJ*GAR*FUELDN(3))
         ENDIF
      ENDIF
      IF(IMPX.GT.0) WRITE(IOUT,6000) PHI*1.0E+13
      IF(INR.GT.0) THEN
         FNORM=PHI*VTOT/VPH
         DO 290 IBM=1,NBMIX
         VPHV(IBM)=VPHV(IBM)*REAL(FNORM)
         DO 280 IQ=1,NREAC
         DO 270 IS=1,NVAR+1
         SIG(IS,IQ,IBM)=SIG(IS,IQ,IBM)*REAL(FNORM)
  270    CONTINUE
  280    CONTINUE
  290    CONTINUE
      ELSE
         CALL XDRSET(VPHV,NBMIX,0.0)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XSREC)
      RETURN
*
 6000 FORMAT(/' EVOSIG: flux level  =',1P,E12.4,' n/cm^2/s.')
      END
