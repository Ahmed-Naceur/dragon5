*DECK EVOBLD
      SUBROUTINE EVOBLD(IMPX,INR,IGLOB,NBMIX,NBISO,NCOMB,ISONAM,YDPL,
     1 VX,MILVO,JM,NVAR,NDFP,NSUPS,NREAC,NPAR,NFISS,XT,EPS1,EPS2,EXPMAX,
     2 H1,ITYPE,IDIRAC,FIT,DELTA,ENERG,KPAR,BPAR,YIELD,IDR,RER,RRD,AWR,
     3 FUELDN,SIG,VPH,VPHV,MIXPWR,VTOTD,IEVOLB,KFISS,KPF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform flux normalization and call EVOSOL to solve the depletion
* system for each depleting mixture between times XT(1) and XT(2).
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
*Parameters: input/output
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
* NBMIX   number of mixtures.
* NBISO   number of isotopes/materials including non-depleting ones.
* NCOMB   number of depleting mixtures.
* ISONAM  alias name of isotopes.
* YDPL    initial/final number density of isotope in the depletion
*         chain. YDPL(NVAR+1,2,ICMB) is the stage burnup increment
*         in region ICMB.
* VX      volumes of the depleting mixtures.
* MILVO   mixture index corresponding to each depleting mixture.
* JM      position in isotope list of each nuclide of the depletion
*         chain.
* NVAR    number of depleting nuclides.
* NDFP    number of direct fission products (fission fragments).
* NSUPS   number of non-depleting isotopes producing energy.
* NREAC   maximum number of depletion reactions.
* NPAR    maximum number of parent nuclides in the depletion chain.
* NFISS   number of fissile isotopes producing fission products.
* XT      initial and final time (independent variable).
* EPS1    required accuracy for the ode solver.
* EPS2    required accuracy for constant power iterations.
* EXPMAX  saturation limit. A nuclide is saturating if
*         -ADPL(MU1(I))*(XT(2)-XT(1)).GT.EXPMAX. Suggested value:
*         EXPMAX=80.0.
* H1      guessed first stepsize.
* ITYPE   type of ODE solution:
*         =1 fifth-order Runge-Kutta method;
*         =2 fourth-order Kaps-Rentrop method.
* IDIRAC  saturation model flag (=1 to use Dirac function contributions
*         in the saturating nuclide number densities).
* FIT     flux normalization factor:
*         n/cm**2/s if INR=1;
*         MW/tonne of initial heavy elements if INR=2;
*         W/cc of assembly volume if INR=3.
* DELTA   burnup stage increments:
*         DELTA(1): increment in fuel burnup for this stage;
*         DELTA(2): increment in fuel neutron exposure for this stage;
*         DELTA(3): target increment in fuel burnup for this stage.
*         Cross section should be tabulated with respect to the sum
*         of the DELTA(1) of all the previous stages.
* ENERG   increment in fuel burnup for this stage in each mixture.
* KPAR    position in chain of the parent nuclide and type of
*         reaction.
* BPAR    branching ratio for neutron induced reactions.
* YIELD   mixture-dependent fission yields.
* IDR     identifier for each depleting reaction.
* RER     energy (Mev) per reaction. If RER(3,J)=0., the fission energy
*         includes radiative capture energy. Neutrino energy is
*         never included.
* RRD     sum of radioactive decay constants in 10**-8/s.
* AWR     mass of the nuclides in unit of neutron mass.
* FUELDN  fuel initial density and mass.
* SIG     initial/final microscopic depletion reaction rates for nuclide
*         I in mixture IBM:
*         SIG(I,1,IBM,:) fission reaction rate;
*         SIG(I,2,IBM,:) gamma reaction rate;
*         SIG(I,3,IBM,:) N2N reaction rate;
*         cont...;
*         SIG(I,NREAC,IBM,:) neutron-induced energy released;
*         SIG(I,NREAC+1,IBM,:) decay energy released (10**-8 MeV/s).
* VPH     initial/final integrated flux in fuel.
* VPHV    initial/final integrated flux in each mixture.
* MIXPWR  flags for mixtures to include in power normalization.
* VTOTD   total fuel volume.
* IEVOLB  flag making an isotope non-depleting:
*         =0 the isotope is depleting;
*         =1 to force an isotope to be non-depleting;
*         =2 to force an isotope to be depleting;
*         =3 to force an isotope to be at saturation
* KFISS   position in chain of the fissile isotopes.
* KPF     position in chain of the direct fission products (fission
*         fragments).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,INR,IGLOB,NBMIX,NBISO,NCOMB,ISONAM(3,NBISO),
     1 MILVO(NCOMB),JM(NBMIX,NVAR+NSUPS),NVAR,NDFP,NSUPS,NREAC,NPAR,
     2 NFISS,ITYPE,IDIRAC,KPAR(NPAR,NVAR),IDR(NREAC,NVAR+NSUPS),
     3 MIXPWR(NBMIX),IEVOLB(NVAR,NBMIX),KFISS(NFISS,NBMIX),
     4 KPF(NDFP,NBMIX)
      REAL YDPL(NVAR+1,2,NCOMB),VX(NBMIX),XT(2),EPS1,EPS2,EXPMAX,H1,FIT,
     1 DELTA(3),ENERG(NBMIX),BPAR(NPAR,NVAR),YIELD(NFISS,NDFP,NBMIX),
     2 RER(NREAC,NVAR+NSUPS),RRD(NVAR+NSUPS),AWR(NVAR),FUELDN(3),
     3 SIG(NVAR+1,NREAC+1,NBMIX,2),VPH(2),VPHV(NBMIX,0:2)
      DOUBLE PRECISION VTOTD
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXT8*8,HSMG*131
      DOUBLE PRECISION GAR,GARD,XDRCST,EVJ,FITD,PHI2
      LOGICAL LCOOL,LSIMPL
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MU1,IMA,LP,CHAIN
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MU1(NVAR+1),IMA(NVAR+1),LP(NVAR))
*----
*  CHECK IF ONLY THE HEAVY ISOTOPES ARE PRODUCING ENERGY. IN THIS CASE,
*  SOME SIMPLIFICATIONS ARE POSSIBLE
*----
      LSIMPL=.TRUE.
      DO 10 IS=1,NVAR
      LSIMPL=LSIMPL.AND.(RER(3,IS).EQ.0.0)
   10 CONTINUE
*
      EVJ=XDRCST('eV','J')*1.0E22
      LCOOL=(INR.EQ.0)
      IF(LCOOL) GO TO 410
      IF(IMPX.GT.1) WRITE (6,640) VPH(1)/VTOTD,VPH(2)/VTOTD
*----
*  CONVERGE ON A FIXED FINAL POWER. SOLVE THE DEPLETION CHAIN WITHOUT
*  THE FISSION PRODUCTS.
*----
      IF((INR.GE.2).AND.(EPS2.LT.10.0)) THEN
         ITER=0
  250    ITER=ITER+1
         IF(ITER.GT.20) CALL XABORT('EVOBLD: UNABLE TO CONVERGE.')
         DO 330 ICMB=1,NCOMB
*        DETERMINE THE ROW AND COLUMN PROFILE OF THE ADPL MATRIX AND
*        COMPUTE NVAR2, THE NUMBER OF DEPLETING NUCLIDES IN REGION ICMB.
         IBM=MILVO(ICMB)
         IF(IMPX.GT.3) WRITE(6,'(/34H EVOBLD: PROCESS DEPLETING MIXTURE,
     1   I5,15H (REAL MIXTURE=,I5,2H).)') ICMB,IBM
         NVAR2=0
         NSUPL2=0
         CALL XDISET(LP,NVAR,0)
         DO 270 IS=1,NVAR
         KDRI=IDR(2,IS)
         IF(KDRI.EQ.0) GO TO 270
         IF((MOD(KDRI,100).NE.3).AND.(MOD(KDRI,100).NE.4)) GO TO 270
         IF(JM(IBM,IS).GT.0) THEN
            NVAR2=NVAR2+1
            LP(IS)=NVAR2
         ENDIF
  270    CONTINUE
         NSUPL2=NVAR2
         DO 280 IS=1,NVAR
         IF(LSIMPL.AND.(AWR(IS).LE.210.0)) GO TO 280
         IF((JM(IBM,IS).GT.0).AND.(LP(IS).EQ.0)) THEN
            NVAR2=NVAR2+1
            LP(IS)=NVAR2
         ENDIF
  280    CONTINUE
         IF(NVAR2.EQ.0) GO TO 330
*        CHECK IF ONLY THE HEAVY ISOTOPES ARE PRODUCING ENERGY. IN
*        THIS CASE, IT IS POSSIBLE TO AVOID THE SOLUTION FOR FISSION
*        PRODUCTS.
         NSUPF2=NVAR2-NSUPL2
         IF(LSIMPL) NSUPF2=0
         CALL EVOMU1(IMPX,NVAR,NREAC,LP,XT,LCOOL,NPAR,KPAR,RRD,
     1   SIG(1,1,IBM,1),SIG(1,1,IBM,2),EXPMAX,IEVOLB(1,IBM),MU1,
     2   IMA,MAXA)
         MU1(NVAR2+1)=IMA(NVAR2)+NVAR2+1
         IMA(NVAR2+1)=IMA(NVAR2)+NVAR2+1
         MAXA=MAXA+10*(NVAR2+1)
         NFISS2=0
         DO 300 I=1,NFISS
         IF(KFISS(I,IBM).EQ.0) GO TO 300
         IF(LP(KFISS(I,IBM)).GT.0) NFISS2=NFISS2+1
  300    CONTINUE
         ALLOCATE(CHAIN(2*(NVAR2+1)))
         DO 310 IS=1,NVAR
         IF(LP(IS).GT.0) THEN
            K=JM(IBM,IS)
            CHAIN((LP(IS)-1)*2+1)=ISONAM(1,K)
            CHAIN((LP(IS)-1)*2+2)=ISONAM(2,K)
         ENDIF
  310    CONTINUE
         TEXT8='*POWER*'
         READ(TEXT8,'(2A4)') (CHAIN(NVAR2*2+I0),I0=1,2)
         CALL EVOSOL(IMPX,LCOOL,NVAR,NREAC,NDFP,NPAR,NFISS,XT,EPS1,
     1   EXPMAX,H1,ITYPE,IDIRAC,RRD,KPAR,BPAR,KFISS(1,IBM),KPF(1,IBM),
     2   YIELD(1,1,IBM),LP,IEVOLB(1,IBM),SIG(1,1,IBM,1),SIG(1,1,IBM,2),
     3   NVAR2,NFISS2,NSUPF2,MU1,IMA,MAXA,YDPL(1,1,ICMB),CHAIN)
*
         DEALLOCATE(CHAIN)
  330    CONTINUE
         GAR=0.0D0
         GARD=0.0D0
         DO 360 IS=1,NVAR
         DO 350 ICMB=1,NCOMB
         IBM=MILVO(ICMB)
         IF(MIXPWR(IBM).EQ.1) THEN
            GAR=GAR+VX(IBM)*YDPL(IS,2,ICMB)*SIG(IS,NREAC,IBM,2)
            GARD=GARD+VX(IBM)*YDPL(IS,2,ICMB)*SIG(IS,NREAC+1,IBM,2)
         ENDIF
  350    CONTINUE
  360    CONTINUE
         IF((IGLOB.EQ.1).OR.(INR.EQ.3)) THEN
            DO 370 IBM=1,NBMIX
            IF(MIXPWR(IBM).EQ.1) THEN
               GAR=GAR+SIG(NVAR+1,NREAC,IBM,2)
               GAR=GAR+SIG(NVAR+1,NREAC+1,IBM,2)
            ENDIF
  370       CONTINUE
         ELSE IF(IGLOB.EQ.0) THEN
            DO 380 ICMB=1,NCOMB
            IBM=MILVO(ICMB)
            IF(MIXPWR(IBM).EQ.1) THEN
               GAR=GAR+SIG(NVAR+1,NREAC,IBM,2)
               GAR=GAR+SIG(NVAR+1,NREAC+1,IBM,2)
            ENDIF
  380       CONTINUE
         ENDIF
         IF(GAR.EQ.0.0D0) CALL XABORT('EVOBLD: UNABLE TO NORMALIZE.')
         IF(INR.EQ.3) THEN
*           FITD IS THE DECAY POWER IN WATT PER CUBIC CENTIMETER.
            FITD=(EVJ*GARD*FUELDN(3))/(FUELDN(1)*VTOTD)
            IF(FITD.GT.FIT) THEN
               WRITE(HSMG,'(35HEVOBLD: NEGATIVE FIT(1) FIT(DECAY)=,1P,
     1         E11.4,12H FIT(INPUT)=,E11.4,1H.)') FITD,FIT
               CALL XABORT(HSMG)
            ENDIF
            PHI2=(FIT-FITD)*FUELDN(1)*VPH(2)/(EVJ*GAR*FUELDN(3))
         ELSE
*           FITD IS THE DECAY POWER IN WATT PER GRAM.
            FITD=(EVJ*GARD)/(FUELDN(1)*VTOTD)
            IF(FITD.GT.FIT) THEN
               WRITE(HSMG,'(35HEVOBLD: NEGATIVE FIT(2) FIT(DECAY)=,1P,
     1         E11.4,12H FIT(INPUT)=,E11.4,1H.)') FITD,FIT
               CALL XABORT(HSMG)
            ENDIF
            PHI2=(FIT-FITD)*FUELDN(1)*VPH(2)/(EVJ*GAR)
         ENDIF
         ERROR=REAL(ABS(PHI2-VPH(2)/VTOTD)/ABS(PHI2))
         DO 400 IBM=1,NBMIX
         VPHV(IBM,2)=VPHV(IBM,2)*REAL(PHI2*VTOTD)/VPH(2)
         DO 395 IQ=1,NREAC
         DO 390 IS=1,NVAR+1
         SIG(IS,IQ,IBM,2)=SIG(IS,IQ,IBM,2)*REAL(PHI2*VTOTD)/VPH(2)
  390    CONTINUE
  395    CONTINUE
  400    CONTINUE
         VPH(2)=REAL(PHI2*VTOTD)
         IF(IMPX.GT.3) THEN
            WRITE (6,650) ITER,ERROR,VPH(1)/VTOTD,VPH(2)/VTOTD
         ENDIF
         IF(ERROR.LT.EPS2) THEN
            IF(IMPX.GT.-1) WRITE(6,'(/29H EVOBLD: POWER CONVERGENCE IN,
     1      I3,19H ITERATIONS. ERROR=,1P,E9.2,1H.)') ITER,ERROR
            GO TO 410
         ELSE
            GO TO 250
         ENDIF
      ENDIF
*----
*  SOLVE THE COMPLETE DEPLETION CHAIN, INCLUDING FISSION PRODUCTS
*----
  410 DELTA(1)=0.0
      CALL XDRSET(ENERG,NBMIX,0.0)
      DO 500 ICMB=1,NCOMB
*     DETERMINE THE ROW AND COLUMN PROFILE OF THE ADPL MATRIX AND
*     COMPUTE NVAR2, THE NUMBER OF DEPLETING NUCLIDES IN REGION ICMB.
      IBM=MILVO(ICMB)
      IF(IMPX.GT.3) WRITE(6,'(/34H EVOBLD: PROCESS DEPLETING MIXTURE,
     1 I5,15H (REAL MIXTURE=,I5,2H).)') ICMB,IBM
      NVAR2=0
      NSUPL2=0
      CALL XDISET(LP,NVAR,0)
      DO 440 IS=1,NVAR
      KDRI=IDR(2,IS)
      IF(KDRI.EQ.0) GO TO 440
      IF((MOD(KDRI,100).NE.3).AND.(MOD(KDRI,100).NE.4)) GO TO 440
      IF(JM(IBM,IS).GT.0) THEN
         NVAR2=NVAR2+1
         LP(IS)=NVAR2
      ENDIF
  440 CONTINUE
      NSUPL2=NVAR2
      DO 450 IS=1,NVAR
      IF((JM(IBM,IS).GT.0).AND.(LP(IS).EQ.0)) THEN
         NVAR2=NVAR2+1
         LP(IS)=NVAR2
      ENDIF
  450 CONTINUE
      CALL EVOMU1(IMPX,NVAR,NREAC,LP,XT,LCOOL,NPAR,KPAR,RRD,
     1 SIG(1,1,IBM,1),SIG(1,1,IBM,2),EXPMAX,IEVOLB(1,IBM),MU1,
     2 IMA,MAXA)
      MU1(NVAR2+1)=IMA(NVAR2)+NVAR2+1
      IMA(NVAR2+1)=IMA(NVAR2)+NVAR2+1
      MAXA=MAXA+10*(NVAR2+1)
      NFISS2=0
      DO 460 I=1,NFISS
      IF(KFISS(I,IBM).EQ.0) GO TO 460
      IF(LP(KFISS(I,IBM)).GT.0) NFISS2=NFISS2+1
  460 CONTINUE
      NSUPF2=NVAR2-NSUPL2
      ALLOCATE(CHAIN(2*(NVAR2+1)))
      DO 470 IS=1,NVAR
      IF(LP(IS).GT.0) THEN
         K=JM(IBM,IS)
         CHAIN((LP(IS)-1)*2+1)=ISONAM(1,K)
         CHAIN((LP(IS)-1)*2+2)=ISONAM(2,K)
      ENDIF
  470 CONTINUE
      TEXT8='*POWER*'
      READ(TEXT8,'(2A4)') (CHAIN(NVAR2*2+I0),I0=1,2)
      CALL EVOSOL(IMPX,LCOOL,NVAR,NREAC,NDFP,NPAR,NFISS,XT,EPS1,
     1 EXPMAX,H1,ITYPE,IDIRAC,RRD,KPAR,BPAR,KFISS(1,IBM),KPF(1,IBM),
     2 YIELD(1,1,IBM),LP,IEVOLB(1,IBM),SIG(1,1,IBM,1),SIG(1,1,IBM,2),
     3 NVAR2,NFISS2,NSUPF2,MU1,IMA,MAXA,YDPL(1,1,ICMB),CHAIN)
*
      DEALLOCATE(CHAIN)
      IF(MIXPWR(IBM).EQ.1) THEN
         DELTA(1)=DELTA(1)+YDPL(NVAR+1,2,ICMB)*VX(IBM)
      ENDIF
  500 CONTINUE
*----
*  BURNUP CALCULATION. TAKE THE CONTRIBUTION OBTAINED FROM THE ODE
*  SOLVER AND ADD THE CONTRIBUTION FROM THE NON-DEPLETING ISOTOPES
*  PRODUCING ENERGY
*----
      DO 510 IBM=1,NBMIX
      IF(MIXPWR(IBM).EQ.1) THEN
         GAR=0.5D0*(SIG(NVAR+1,NREAC,IBM,1)+SIG(NVAR+1,NREAC,IBM,2)
     1             +SIG(NVAR+1,NREAC+1,IBM,1)+SIG(NVAR+1,NREAC+1,IBM,2))
         ENERG(IBM)=ENERG(IBM)+REAL(GAR*(XT(2)-XT(1))*EVJ)
      ENDIF
  510 CONTINUE
      DO 516 ICMB=1,NCOMB
      IBM=MILVO(ICMB)
      IF(MIXPWR(IBM).EQ.1) THEN
         DO 515 IS=1,NVAR
         GAR=0.5D0*(YDPL(IS,1,ICMB)*SIG(IS,NREAC,IBM,1)
     1             +YDPL(IS,2,ICMB)*SIG(IS,NREAC,IBM,2)
     1             +YDPL(IS,1,ICMB)*SIG(IS,NREAC+1,IBM,1)
     1             +YDPL(IS,2,ICMB)*SIG(IS,NREAC+1,IBM,2))
         ENERG(IBM)=ENERG(IBM)+REAL(GAR*(XT(2)-XT(1))*VX(IBM)*EVJ)
  515    CONTINUE
      ENDIF
  516 CONTINUE
      DELTA(3)=0.0
      IF(IGLOB.EQ.0) THEN
         DO 520 ICMB=1,NCOMB
         IBM=MILVO(ICMB)
         IF(MIXPWR(IBM).EQ.1) THEN
            GAR=0.5D0*(SIG(NVAR+1,NREAC,IBM,1)+SIG(NVAR+1,NREAC,IBM,2)
     1                +SIG(NVAR+1,NREAC+1,IBM,1)
     2                +SIG(NVAR+1,NREAC+1,IBM,2))
            DELTA(1)=DELTA(1)+REAL(GAR)*(XT(2)-XT(1))
            DELTA(3)=DELTA(3)+REAL(GAR)*(XT(2)-XT(1))
         ENDIF
  520    CONTINUE
      ELSE IF(IGLOB.EQ.1) THEN
         DO 530 IBM=1,NBMIX
         IF(MIXPWR(IBM).EQ.1) THEN
            GAR=0.5D0*(SIG(NVAR+1,NREAC,IBM,1)+SIG(NVAR+1,NREAC,IBM,2)
     1                +SIG(NVAR+1,NREAC+1,IBM,1)
     2                +SIG(NVAR+1,NREAC+1,IBM,2))
            DELTA(1)=DELTA(1)+REAL(GAR)*(XT(2)-XT(1))
            DELTA(3)=DELTA(3)+REAL(GAR)*(XT(2)-XT(1))
         ENDIF
  530    CONTINUE
      ENDIF
      DO 545 ICMB=1,NCOMB
      IBM=MILVO(ICMB)
      IF(MIXPWR(IBM).EQ.1) THEN
         DO 540 IS=1,NVAR
         DELTA(3)=DELTA(3)+0.5*(YDPL(IS,1,ICMB)*SIG(IS,NREAC,IBM,1)
     1                 +YDPL(IS,2,ICMB)*SIG(IS,NREAC,IBM,2))*VX(IBM)
     2                 *(XT(2)-XT(1))
         DELTA(3)=DELTA(3)+0.5*(YDPL(IS,1,ICMB)*SIG(IS,NREAC+1,IBM,1)
     1                 +YDPL(IS,2,ICMB)*SIG(IS,NREAC+1,IBM,2))*VX(IBM)
     2                 *(XT(2)-XT(1))
  540    CONTINUE
      ENDIF
  545 CONTINUE
      IF(FUELDN(2) .EQ. 0.0) THEN
        DELTA(1)=0.0
        DELTA(3)=0.0
      ELSE
        DELTA(1)=DELTA(1)*REAL(EVJ)/FUELDN(2)
        DELTA(3)=DELTA(3)*REAL(EVJ)/FUELDN(2)
      ENDIF
      DELTA(2)=0.5*(VPH(1)+VPH(2))*(XT(2)-XT(1))/REAL(VTOTD)
      IF((.NOT.LCOOL).AND.(IMPX.GT.0)) THEN
         IF(DELTA(1) .EQ. 0.0) THEN
           WRITE (6,661) DELTA(2)
         ELSE
           WRITE (6,660) DELTA(1),DELTA(1)/8.64E-4,DELTA(2)
           WRITE (6,665) DELTA(3),DELTA(3)/8.64E-4
         ENDIF
      ENDIF
*----
*  PRINT THE BEGINNING- AND END-OF-STAGE NORMALIZATION POWERS
*----
      IF(IMPX.GT.-1) THEN
         DO 580 IP=1,2
         DELTA1=0.0
         DELTA2=0.0
         IF(IGLOB.EQ.0) THEN
            DO 550 ICMB=1,NCOMB
            IBM=MILVO(ICMB)
            IF(MIXPWR(IBM).EQ.1) THEN
               DELTA1=DELTA1+SIG(NVAR+1,NREAC,IBM,IP)
               DELTA2=DELTA2+SIG(NVAR+1,NREAC+1,IBM,IP)
            ENDIF
  550       CONTINUE
         ELSE IF(IGLOB.EQ.1) THEN
            DO 560 IBM=1,NBMIX
            IF(MIXPWR(IBM).EQ.1) THEN
               DELTA1=DELTA1+SIG(NVAR+1,NREAC,IBM,IP)
               DELTA2=DELTA2+SIG(NVAR+1,NREAC+1,IBM,IP)
            ENDIF
  560       CONTINUE
         ENDIF
         DO 575 ICMB=1,NCOMB
         IBM=MILVO(ICMB)
         IF(MIXPWR(IBM).EQ.1) THEN
           DO 570 IS=1,NVAR
           DELTA1=DELTA1+YDPL(IS,IP,ICMB)*SIG(IS,NREAC,IBM,IP)*VX(IBM)
           DELTA2=DELTA2+YDPL(IS,IP,ICMB)*SIG(IS,NREAC+1,IBM,IP)*VX(IBM)
  570      CONTINUE
         ENDIF
  575    CONTINUE
         IF(FUELDN(2) .EQ. 0.0) THEN
           DELTA1=0.0
           DELTA2=0.0
         ELSE
           DELTA1=DELTA1*REAL(EVJ)/FUELDN(2)
           DELTA2=DELTA2*REAL(EVJ)/FUELDN(2)
           IF(IP.EQ.1) THEN
              WRITE(6,680) 'BEGINNING-OF-STAGE',DELTA1
              WRITE(6,690) 'BEGINNING-OF-STAGE',DELTA2
           ELSE IF(IP.EQ.2) THEN
              WRITE(6,680) 'END-OF-STAGE',DELTA1
              WRITE(6,690) 'END-OF-STAGE',DELTA2
           ENDIF
         ENDIF
  580    CONTINUE
         WRITE(6,'(/48H NOTE: POWER MAY EXIBITS VARIATIONS OUTSIDE THE ,
     1   52HBEGINNING- AND END-OF-STAGE VALUES DURING THE STAGE.)')
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(LP,IMA,MU1)
      RETURN
*
  640 FORMAT(/53H EVOBLD: STARTING VALUES OF INITIAL/FINAL AVERAGED FL,
     1 11HUX IN FUEL=,1P,2E12.4,15H E+13 N/CM**2/S)
  650 FORMAT(/37H EVOBLD: ITERATION ON FINAL POWER NB.,I3,3X,6HERROR=,
     1 1P,E12.4,3X,36HINITIAL/FINAL AVERAGED FLUX IN FUEL=,2E12.4,
     2 15H E+13 N/CM**2/S)
  660 FORMAT(/' EVOBLD: ',
     1 'FUEL BURNUP INCREMENT DURING THIS STAGE =',1P,
     1 E12.4,' MW*S**8/TONNE (',E12.4,' MW*DAY/TONNE).'/9X,
     2 'NEUTRON EXPOSURE (FLUENCE) INCREMENT =',E12.4,' N/KB.')
  661 FORMAT(/' EVOBLD: ',
     1 'NEUTRON EXPOSURE (FLUENCE) INCREMENT DURING THIS STAGE =',1P,
     2 E12.4,' N/KB.')
  665 FORMAT(' EVOBLD: ',
     1 'TARGET FUEL BURNUP INCREMENT DURING THIS STAGE =',1P,
     2 E12.4,' MW*S**8/TONNE (',E12.4,' MW*DAY/TONNE).')
  680 FORMAT(/34H EVOBLD: NEUTRON-INDUCED POWER AT ,A,2H =,1P,E12.4,
     > 10H MW/TONNE.)
  690 FORMAT(/24H EVOBLD: DECAY POWER AT ,A,2H =,1P,E12.4,10H MW/TONNE.)
      END
