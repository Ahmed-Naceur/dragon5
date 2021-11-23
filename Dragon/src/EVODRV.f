*DECK EVODRV
      SUBROUTINE EVODRV(IPDEPL,IPLIB,INDREC,IMPX,NBISO,NGROUP,NBMIX,
     1 ISONAM,ISONRF,MIX,DEN,IEVOL,ISTYP,VX,NDEPL,NSUPS,NREAC,NCOMB,
     2 EPS1,EPS2,EXPMAX,H1,ITYPE,INR,IEXTR,IGLOB,ISAT,IDIRAC,ITIXS,
     3 IFLMAC,IYLMIX,FIT,ISAVE,ISET,IDEPL,XTI,XTF,XT,LMACRO,FLUMIX,
     4 IPICK,MIXBRN,MIXPWR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Isotopic depletion calculation main driver.
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
* IPDEPL  pointer to the depletion history (L_BURNUP signature).
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* INDREC  beginning of depletion flag (=1: beginning; =2: else).
* IMPX    print flag (equal to zero for no print).
* NBISO   number of isotopes present in the calculation domain.
* NGROUP  number of energy groups.
* NBMIX   number of mixtures.
* ISONAM  alias name of isotopes.
* ISONRF  library name of isotopes.
* MIX     mix number of each isotope (can be zero).
* DEN     density of each isotope.
* IEVOL   non-depleting flag (=1 to force an isotope to be
*         non-depleting; =2 to force an isotope to at saturation).
* ISTYP   isotope type (=1 not fissile nor fission product; =2: fissile;
*         =3: fission product).
* VX      volume occupied by each mixture.
* NDEPL   number of isotopes in the depletion chain.
* NSUPS   number of non-depleting isotopes producing energy.
* NREAC   maximum number of depletion reactions.
* NCOMB   number of depleting mixtures.
* EPS1    required accuracy for the ODE solver.
* EPS2    required accuracy for constant power iterations.
* EXPMAX  saturation limit. A nuclide is saturating if
*         -ADPL(MU1(I))*(XT(2)-XT(1)).GT.EXPMAX. Suggested value:
*         EXPMAX=80.0.
* H1      guessed first stepsize.
* ITYPE   type of ODE solution:
*         =1 fifth-order Runge-Kutta method;
*         =2 fourth-order Kaps-Rentrop method.
* INR     type of flux normalization:
*         =0: out-of-core depletion;
*         =1: constant flux depletion;
*         =2: constant fuel power depletion;
*         =3: constant assembly power depletion.
* IEXTR   flux extrapolation flag (=0: no extrapolation; =1: perform
*         extrapolation).
* IGLOB   out-of-fuel power in flux normalization:
*         =0: compute the burnup using the power released in the fuel;
*         =1: compute the burnup using the power released in the global
*         geometry.
* ISAT    initial saturation flag (=1 to save initial saturated number
*         densities).
* IDIRAC  saturation model flag (=1 to use Dirac function contributions
*         in the saturating nuclide number densities).
* ITIXS   flag for time-dependent cross sections (=0/1: on/off).
* IFLMAC  0/1/2 flag to recover fluxes from L_FLUX/L_MACROLIB/L_POWER.
* IYLMIX  0/1 flag to recover fission yield data from DEPL-CHAIN/PYIELD
*         data.
* FIT     flux normalization factor:
*         n/cm**2/s if INR=1;
*         MW/tonne of initial heavy elements if INR=2;
*         W/cc of assembly volume if INR=3.
* ISAVE   save flag:
*         =-1: do not save the last flux calculation in the depletion
*         table;
*         .GE.0 save the last flux calculation in the depletion
*         table at time XTI.
* ISET    set flag:
*         =-1: do not set the number densities to a selected time;
*         .GE.0 set the number densities to time XTF of the depletion
*         table.
* IDEPL   depletion flag:
*         =0: do not perform a depletion calculation
*         =1: perform a depletion calculation.
* XTI     initial save time (save the last flux calculation in the
*         depletion table at time XTI).
* XTF     final set time (recover the number densities from the
*         depletion table at time XTF and modify the internal library).
* XT      time variable (independent variable) for the depletion
*         calculation.
*         XT(1) initial time
*         XT(2) final time
* LMACRO  macrolib building flag (=.true. to compute the embedded
*         macrolib).
* FLUMIX  average fluxes in mixtures.
* IPICK   burnup recovery flag:
*         =0: do not recover the burnup in a CLE-2000 variable
*         =1: recover the burnup in a CLE-2000 variable.
* MIXBRN  flags for mixtures to burn.
* MIXPWR  flags for mixtures to include in power normalization.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDEPL,IPLIB
      INTEGER INDREC,IMPX,NBISO,NGROUP,NBMIX,ISONAM(3,NBISO),
     1 ISONRF(3,NBISO),MIX(NBISO),IEVOL(NBISO),ISTYP(NBISO),NDEPL,
     2 NSUPS,NREAC,NCOMB,ITYPE,INR,IEXTR,IGLOB,ISAT,IDIRAC,ITIXS,
     3 IFLMAC,IYLMIX,ISAVE,ISET,IDEPL,IPICK,MIXBRN(NBMIX),
     4 MIXPWR(NBMIX)
      REAL DEN(NBISO),VX(NBMIX),EPS1,EPS2,EXPMAX,H1,FIT,XTI,XTF,
     1 XT(2),FLUMIX(NGROUP,NBMIX)
      LOGICAL LMACRO
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IUNOUT=6,NSTATE=40,MAXTIM=1000)
      TYPE(C_PTR) KPLIB
      CHARACTER TEXT12*12,HSMG*131
      LOGICAL LCOOL
      INTEGER IDIM(NSTATE),IPAR(NSTATE)
      REAL TIMES(MAXTIM),DELTA(3),RPAR(5),BRNWIR(2),TMPDAY(3),VPH(2),
     1 FUELDN(3),DELTAT(2,2),TIMEP(2,3)
      DOUBLE PRECISION XDRCST,AVCON,VTOTD,VPHINI,DBLLIR
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MILVO,ISOCMB,NFISS2,
     1 NDFP2,HREAC,IPIFI
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: JM,INADPL,IEVOLB,KFISS,
     1 KPAR,IDR,KPF
      REAL, ALLOCATABLE, DIMENSION(:) :: ENERG,RRD,AWR,PYIELD
      REAL, ALLOCATABLE, DIMENSION(:,:) :: BPAR,RER,YIELD2,VPHV
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: YDPL,YIELD
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: SIG
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASK,MASKL
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  SCRATCH STORAGE ALLOCATION
*   MILVO   mixture index corresponding to each depleting mixture.
*----
      ALLOCATE(JM(NBMIX,NDEPL),MILVO(NCOMB),ISOCMB(NBISO),
     1 INADPL(3,NDEPL),IEVOLB(NDEPL,NBMIX))
      ALLOCATE(SIG(NDEPL-NSUPS+1,NREAC+1,NBMIX,0:2),VPHV(NBMIX,0:2),
     1 ENERG(NBMIX),AWR(NDEPL),YDPL(NDEPL-NSUPS+1,2,NCOMB))
      ALLOCATE(MASK(NBMIX),MASKL(NGROUP))
      ALLOCATE(IPISO(NBISO))
*----
*  INITIALIZE DATA.
*----
      AVCON=XDRCST('Neutron mass','amu')/
     1      (1.0D-24*XDRCST('Avogadro','N/moles'))
*
      IF(INDREC.EQ.1) THEN
*        BEGINNING OF DEPLETION.
         NTIM=0
      ELSE IF(INDREC.EQ.2) THEN
         CALL LCMLEN(IPDEPL,'DEPL-TIMES',NTIM,ITYLCM)
      ENDIF
      IF(NTIM.EQ.0) THEN
        CALL XDRSET(TIMES,MAXTIM,0.0)
        CALL LCMPUT(IPDEPL,'DEPL-TIMES',1,2,TIMES)
      ENDIF
*----
*  RECOVER DEPLETION CHAIN INFO FROM LCM
*----
      NVAR=NDEPL-NSUPS
      CALL LCMLEN(IPLIB,'DEPL-CHAIN',LENGT,ITYLCM)
      IF (LENGT.EQ.0) CALL XABORT('EVODRV: DEPLETION CHAIN MISSING.')
      CALL LCMSIX(IPLIB,'DEPL-CHAIN',1)
      CALL LCMGET(IPLIB,'STATE-VECTOR',IDIM)
      IF((NDEPL.NE.IDIM(1)).OR.(NSUPS.NE.IDIM(7)))
     1 CALL XABORT('EVODRV: INCONSISTENT NDEPL OR NSUPS.')
      IF(NVAR.EQ.0) CALL XABORT('EVODRV: NO DEPLETING ISOTOPES')
      NFISS=IDIM(2)
      NDFP=IDIM(3)
      NPAR=IDIM(9)
      ALLOCATE(KPAR(NDEPL,NPAR),HREAC(2*NREAC),IDR(NREAC,NDEPL))
      ALLOCATE(BPAR(NDEPL,NPAR),YIELD2(NFISS,NDFP),RER(NREAC,NDEPL),
     1 RRD(NDEPL))
      CALL LCMGET(IPLIB,'ISOTOPESDEPL',INADPL)
      CALL LCMGET(IPLIB,'PRODUCE-REAC',KPAR)
      CALL LCMGET(IPLIB,'PRODUCE-RATE',BPAR)
      CALL LCMGET(IPLIB,'DEPLETE-IDEN',HREAC)
      CALL LCMGET(IPLIB,'DEPLETE-REAC',IDR)
      CALL LCMGET(IPLIB,'DEPLETE-ENER',RER)
      CALL LCMGET(IPLIB,'DEPLETE-DECA',RRD)
      IF(NFISS*NDFP.GT.0) CALL LCMGET(IPLIB,'FISSIONYIELD',YIELD2)
      CALL LCMSIX(IPLIB,' ',2)
*----
*  SET THE LCM MICROLIB ISOTOPEWISE DIRECTORIES.
*----
      CALL LIBIPS(IPLIB,NBISO,IPISO)
*----
*  DETECT THE DEPLETING ISOTOPES AND MIXTURES IN THE MICROLIB.
*----
      ICOMB=0
      CALL XDISET(JM,NBMIX*NDEPL,0)
      CALL XDISET(IEVOLB,NDEPL*NBMIX,1)
      ALLOCATE(NFISS2(NBMIX),NDFP2(NBMIX))
      CALL XDISET(NFISS2,NBMIX,0)
      CALL XDISET(NDFP2,NBMIX,0)
      CALL XDRSET(AWR,NDEPL,0.0)
      DO 30 ISOT=1,NBISO
      IBM=MIX(ISOT)
      IF(IBM.EQ.0) GO TO 30
      IF(MIXBRN(IBM).EQ.0) GO TO 30
      KPLIB=IPISO(ISOT) ! set ISOT-th isotope
      DO 20 INUCL=1,NDEPL
      IF((ISONRF(1,ISOT).EQ.INADPL(1,INUCL)).AND.(ISONRF(2,ISOT).EQ.
     1 INADPL(2,INUCL)).AND.(ISONRF(3,ISOT).EQ.INADPL(3,INUCL))) THEN
         IF(JM(IBM,INUCL).GT.0) GO TO 20
         IF(C_ASSOCIATED(KPLIB)) THEN
           CALL LCMLEN(KPLIB,'AWR',ILONG,ITYLCM)
           IF(ILONG.EQ.1) CALL LCMGET(KPLIB,'AWR',AWR(INUCL))
         ENDIF
         IF(INUCL.GT.NVAR) THEN
            IF(JM(IBM,INUCL).EQ.0) THEN
               IEVOLB(INUCL,IBM)=1
               JM(IBM,INUCL)=-ISOT
            ENDIF
         ELSE IF(IEVOL(ISOT).EQ.1) THEN
            IF(JM(IBM,INUCL).EQ.0) THEN
               IEVOLB(INUCL,IBM)=1
               JM(IBM,INUCL)=-ISOT
            ENDIF
            IF(ISTYP(ISOT).EQ.2) THEN
               NFISS2(IBM)=NFISS2(IBM)+1
               JM(IBM,INUCL)=ISOT
            ENDIF
         ELSE
            IF(ISTYP(ISOT).EQ.2) THEN
               NFISS2(IBM)=NFISS2(IBM)+1
            ELSE IF(ISTYP(ISOT).EQ.3) THEN
               NDFP2(IBM)=NDFP2(IBM)+1
            ENDIF
            IEVOLB(INUCL,IBM)=IEVOL(ISOT)
            JM(IBM,INUCL)=ISOT
            DO 10 J=1,ICOMB
            IF(IBM.EQ.MILVO(J)) GO TO 30
   10       CONTINUE
            ICOMB=ICOMB+1
            MILVO(ICOMB)=IBM
         ENDIF
         GO TO 30
      ENDIF
   20 CONTINUE
   30 CONTINUE
      IF(ICOMB.NE.NCOMB) CALL XABORT('EVODRV: INVALID VALUE OF NCOMB.')
      IF(IYLMIX.EQ.1) THEN
         NFISS=MAXVAL(NFISS2)
         NDFP=MAXVAL(NDFP2)
      ENDIF
      DEALLOCATE(NDFP2,NFISS2)
      IF(IMPX.GT.0) WRITE(IUNOUT,500) NFISS,NDFP
*----
*  SET KFISS, KPF AND YIELD
*----
      ALLOCATE(KFISS(NFISS,NBMIX),KPF(NDFP,NBMIX),
     1         YIELD(NFISS,NDFP,NBMIX))
      CALL XDISET(KFISS,NFISS*NBMIX,0)
      CALL XDISET(KPF,NDFP*NBMIX,0)
      IF(IYLMIX.EQ.0) THEN
*        Use fission yield data from 'DEPL-CHAIN'
         DO 40 IS=1,NVAR
         KDRI=IDR(2,IS)/100
         IF((KDRI.GT.0).AND.(MOD(IDR(2,IS),100).EQ.4)) THEN
            IF(KDRI.GT.NFISS) CALL XABORT('EVODRV: INVALID NFISS.')
            KFISS(KDRI,:NBMIX)=IS
         ELSE IF((KDRI.GT.0).AND.(MOD(IDR(2,IS),100).EQ.5)) THEN
            IF(KDRI.GT.NDFP) CALL XABORT('EVODRV: INVALID NDFP.')
            KPF(KDRI,:NBMIX)=IS
         ENDIF
   40    CONTINUE
         DO IDFP=1,NDFP
            DO IFISS=1,NFISS
               YIELD(IFISS,IDFP,:NBMIX)=YIELD2(IFISS,IDFP)
            ENDDO
         ENDDO
      ELSE IF(IYLMIX.EQ.1) THEN
*        Use isotopic PIFI/PYIELD fission yield data
         CALL XDRSET(YIELD,NFISS*NDFP*NBMIX,0.0)
         DO 65 IBM=1,NBMIX
         IFISS=0
         IDFP=0
         DO 50 IS=1,NVAR
         ISOT=ABS(JM(IBM,IS))
         IF(ISOT.EQ.0) GO TO 50
         IF(ISTYP(ISOT).EQ.2) THEN
            IFISS=IFISS+1
            IF(IFISS.GT.NFISS) CALL XABORT('EVODRV: NFISS OVERFLOW.')
            KFISS(IFISS,IBM)=IS
            IF(IDR(2,IS).EQ.0) IDR(2,IS)=4
         ELSE IF(ISTYP(ISOT).EQ.3) THEN
            IDFP=IDFP+1
            IF(IDFP.GT.NDFP) CALL XABORT('EVODRV: NDFP OVERFLOW.')
            KPF(IDFP,IBM)=IS
            IF(IDR(2,IS).EQ.0) IDR(2,IS)=5
         ENDIF
   50    CONTINUE
         DO 60 IS=1,NVAR
         ISOT=JM(IBM,IS)
         IF(ISOT.LE.0) GO TO 60
         KPLIB=IPISO(ISOT) ! set ISOT-th isotope
         IF(.NOT.C_ASSOCIATED(KPLIB)) THEN
           WRITE(HSMG,'(17HEVODRV: ISOTOPE '',3A4,16H'' IS NOT AVAILAB,
     >     22HLE IN THE MICROLIB(1).)') (ISONAM(I0,ISOT),I0=1,3)
           CALL XABORT(HSMG)
         ENDIF
         CALL LCMLEN(KPLIB,'PIFI',NDFI,ITYLCM)
         IF(NDFI.GT.0) THEN
            ALLOCATE(IPIFI(NDFI),PYIELD(NDFI))
            CALL LCMGET(KPLIB,'PIFI',IPIFI)
            CALL LCMGET(KPLIB,'PYIELD',PYIELD)
            IDFP=0
            DO I=1,NDFP
               IF(KPF(I,IBM).EQ.IS) THEN
                  IDFP=I
                  EXIT
               ENDIF
            ENDDO
            IF(IDFP.EQ.0) THEN
               WRITE(IUNOUT,510) 'FISSION PRODUCT',IS,(KPF(I,IBM),
     1         I=1,NDFP)
               WRITE(HSMG,'(39HEVODRV: UNABLE TO FIND FP INDEX FOR ISO,
     1         5HTOPE ,3A4,5H (1).)')(ISONAM(I0,ISOT),I0=1,3)
               CALL XABORT(HSMG)
            ENDIF
            DO 55 I=1,NDFI
               IF(IPIFI(I).EQ.0) GO TO 55
               DO JST=1,NVAR
                  IF(ABS(JM(IBM,JST)).EQ.IPIFI(I)) THEN
                     IFISS=0
                     DO J=1,NFISS
                        IF(KFISS(J,IBM).EQ.JST) THEN
                           IFISS=J
                           EXIT
                        ENDIF
                     ENDDO
                     IF(IFISS.EQ.0) THEN
                        WRITE(IUNOUT,510) 'FISSILE ISOTOPE',JST,
     1                  (KFISS(J,IBM),J=1,NFISS)
                        CALL XABORT('EVODRV: UNABLE TO FIND FISSILE I'
     1                  //'SOTOPE INDEX')
                     ENDIF
                     YIELD(IFISS,IDFP,IBM)=PYIELD(I)
                     GO TO 55
                  ENDIF
               ENDDO
               WRITE(HSMG,'(39HEVODRV: UNABLE TO FIND FP INDEX FOR ISO,
     1         5HTOPE ,3A4,5H (2).)')(ISONAM(I0,ISOT),I0=1,3)
               CALL XABORT(HSMG)
   55       CONTINUE
            DEALLOCATE(PYIELD,IPIFI)
         ENDIF
   60    CONTINUE
   65    CONTINUE
      ELSE
         CALL XABORT('EVODRV: INVALID VALUE OF FLAG IYLMIX.')
      ENDIF
      DEALLOCATE(YIELD2)
*----
*  COMPUTE THE INITIAL INTEGRATED FLUX
*----
      VPHINI=0.0D0
      DO 85 IU=1,NGROUP
      DO 80 ICMB=1,NCOMB
      IBM=MILVO(ICMB)
      VPHINI=VPHINI+FLUMIX(IU,IBM)*VX(IBM)
   80 CONTINUE
   85 CONTINUE
*----
*  CHECK IF PERTURBATION XS OR STANDARD XS.
*----
      NXSPER=1
      CALL LCMLEN(IPLIB,'TIMESPER',LENGTH,ITYLCM)
      IF((LENGTH.GE.2).AND.(LENGTH.LE.6)) THEN
         CALL LCMGET(IPLIB,'TIMESPER',TIMEP)
         DELTAT(1,1)=TIMEP(1,1)
         DELTAT(2,1)=TIMEP(2,1)
         TMPREF=DELTAT(1,1)
         NXSPER=2
         IF(ITIXS.EQ.0) THEN
            DO 90 IP=1,2
            DELTAT(1,IP)=1.0
            DELTAT(2,IP)=XT(IP)/8.64E-4-TMPREF
   90       CONTINUE
         ELSE
            XREF=XT(1)/8.64E-4
            DO 100 IP=1,2
            DELTAT(1,IP)=1.0
            DELTAT(2,IP)=XREF-TMPREF
  100       CONTINUE
         ENDIF
      ELSE
         DO 110 IP=1,2
         DELTAT(1,IP)=1.0
         DELTAT(2,IP)=0.0
  110    CONTINUE
      ENDIF
*----
*  COMPUTE AND SAVE THE INITIAL MASS OF HEAVY ELEMENTS
*----
      VTOTD=0.0D0
      DO ICMB=1,NCOMB
        IBM=MILVO(ICMB)
        IF(MIXPWR(IBM).GT.0) VTOTD=VTOTD+DBLE(VX(IBM))
      ENDDO
      IF(INDREC.EQ.1) THEN
*
*        COMPUTE THE GLOBAL HEAVY-ELEMENT MASS FOR ISOTOPES IN MIXPWR.
         FUELDN(1)=0.0
         IF(IGLOB.EQ.0) THEN
           DO 120 ICMB=1,NCOMB
           IBM=MILVO(ICMB)
           IF(MIXPWR(IBM).GT.0) THEN
             DO 115 IS=1,NBISO
             IF((MIX(IS).NE.IBM).OR.(IEVOL(IS).NE.0)) GO TO 115
             KPLIB=IPISO(IS) ! set IS-th isotope
             IF(C_ASSOCIATED(KPLIB)) THEN
               AWRGAR=0.0
               CALL LCMLEN(KPLIB,'AWR',ILONG,ITYLCM)
               IF(ILONG.EQ.1) CALL LCMGET(KPLIB,'AWR',AWRGAR)
               IF((AWRGAR.GT.210.0).OR.(ISTYP(IS).EQ.2)) THEN
                 FUELDN(1)=FUELDN(1)+AWRGAR*DEN(IS)*VX(IBM)
               ENDIF
             ENDIF
  115        CONTINUE
           ENDIF
  120      CONTINUE
         ELSE IF(IGLOB.EQ.1) THEN
           DO 125 IBM=1,NBMIX
           IF(MIXPWR(IBM).GT.0) THEN
             DO 124 IS=1,NBISO
             IF(MIX(IS).NE.IBM) GO TO 124
             KPLIB=IPISO(IS) ! set IS-th isotope
             IF(C_ASSOCIATED(KPLIB)) THEN
               AWRGAR=0.0
               CALL LCMLEN(KPLIB,'AWR',ILONG,ITYLCM)
               IF(ILONG.EQ.1) CALL LCMGET(KPLIB,'AWR',AWRGAR)
               IF((AWRGAR.GT.210.0).OR.(ISTYP(IS).EQ.2)) THEN
                 FUELDN(1)=FUELDN(1)+AWRGAR*DEN(IS)*VX(IBM)
               ENDIF
             ENDIF
  124        CONTINUE
           ENDIF
  125      CONTINUE
         ENDIF
         IF(FUELDN(1).EQ.0) THEN
           IF(INR.LE.1) THEN
             FUELDN(1)=0.0
           ELSE
             CALL XABORT('EVODRV: Burnup at fixed power without '//
     1       'heavy fissile isotopes is forbidden')
           ENDIF
         ELSE
           FUELDN(1)=FUELDN(1)*REAL(AVCON/VTOTD)
         ENDIF
         FUELDN(2)=FUELDN(1)*REAL(VTOTD)
         VASSMB=0.0
         DO 130 IBM=1,NBMIX
         IF(MIXPWR(IBM).EQ.1) THEN
            VASSMB=VASSMB+VX(IBM)
         ENDIF
  130    CONTINUE
         FUELDN(3)=FUELDN(2)/VASSMB
         CALL LCMPUT(IPDEPL,'FUELDEN-INIT',3,2,FUELDN)
*
*        COMPUTE THE HEAVY-ELEMENT MASS PER MIXTURE.
         DO 150 IBM=1,NBMIX
         ENERG(IBM)=0.0
         DO 140 IS=1,NBISO
         IF((MIX(IS).NE.IBM).OR.(IEVOL(IS).NE.0)) GO TO 140
         KPLIB=IPISO(IS) ! set IS-th isotope
         IF(C_ASSOCIATED(KPLIB)) THEN
           AWRGAR=0.0
           CALL LCMLEN(KPLIB,'AWR',ILONG,ITYLCM)
           IF(ILONG.EQ.1) CALL LCMGET(KPLIB,'AWR',AWRGAR)
           IF((AWRGAR.GT.210.0).OR.(ISTYP(IS).EQ.2)) THEN
             ENERG(IBM)=ENERG(IBM)+AWRGAR*DEN(IS)*VX(IBM)
           ENDIF
         ENDIF
  140    CONTINUE
         ENERG(IBM)=ENERG(IBM)*REAL(AVCON)
  150    CONTINUE
         CALL LCMPUT(IPDEPL,'FUELDEN-MIX',NBMIX,2,ENERG)
*
*        COMPUTE THE TOTAL MASS PER MIXTURE.
         CALL XDRSET(ENERG,NBMIX,0.0)
         DO 170 IS=1,NBISO
         IF(DEN(IS).EQ.0.0) GO TO 170
         KPLIB=IPISO(IS) ! set IS-th isotope
         AWRGAR=0.0
         IF(C_ASSOCIATED(KPLIB)) THEN
            CALL LCMLEN(KPLIB,'AWR',ILONG,ITYLCM)
            IF(ILONG.EQ.1) CALL LCMGET(KPLIB,'AWR',AWRGAR)
         ENDIF
         DO 160 IBM=1,NBMIX
         IF(MIX(IS).EQ.IBM) THEN
            ENERG(IBM)=ENERG(IBM)+REAL(AWRGAR*DEN(IS)*VX(IBM)*AVCON)
         ENDIF
  160    CONTINUE
  170    CONTINUE
         CALL LCMPUT(IPDEPL,'WEIGHT-MIX',NBMIX,2,ENERG)
      ELSE
         CALL LCMGET(IPDEPL,'FUELDEN-INIT',FUELDN)
      ENDIF
      IF(IMPX.GT.0) THEN
         WRITE (6,610) FUELDN(1),VTOTD,FUELDN(2),FUELDN(3)
      ENDIF
*----
*  SAVE THE LAST FLUX CALCULATION SET POINT IN THE DEPLETION TABLE.
*  CROSS-SECTION PERTURBATION IS ENABLED IF ISAVE=0 AND NXSPER=2.
*----
      IF(ISAVE.EQ.0) THEN
        DO 200 IP=1,NXSPER
        ITIM=0
        IF(NTIM.GT.0) THEN
          CALL XDRSET(TIMES,MAXTIM,0.0)
          CALL LCMGET(IPDEPL,'DEPL-TIMES',TIMES)
          DO 180 I=1,NTIM
            IF(ABS(TIMES(I)-XT(IP)).LE.1.0E-4*XT(IP)) ITIM=I
  180     CONTINUE
        ENDIF
        IF(ITIM.EQ.0) THEN
          IF(NTIM.GT.0) THEN
            IF(XT(IP).LT.TIMES(NTIM)) CALL XABORT('EVODRV: INVALID X1.')
          ENDIF
          NTIM=NTIM+1
          IF(NTIM.GT.MAXTIM) CALL XABORT('EVODRV: INVALID MAXTIM(1).')
          TIMES(NTIM)=XT(IP)
          CALL LCMPUT(IPDEPL,'DEPL-TIMES',NTIM,2,TIMES)
          ITIM=NTIM
        ENDIF
        WRITE(TEXT12,'(8HDEPL-DAT,I4.4)') ITIM
        IF(IMPX.GT.0) WRITE(IUNOUT,530) XT(IP),XT(IP)/8.64E-4,TEXT12
        CALL LCMSIX(IPDEPL,TEXT12,1)
        CALL LCMPUT(IPDEPL,'ISOTOPESDENS',NBISO,2,DEN)
*----
*  COMPUTE, NORMALIZE AND SAVE THE MICROSCOPIC REACTION RATES.
*----
        CALL EVOSIG(IMPX,INR,IGLOB,NGROUP,NBMIX,NBISO,NCOMB,ISONAM,
     1  IPISO,DEN,FLUMIX,VX,MILVO,JM,NVAR,NSUPS,NREAC,HREAC,IDR,
     2  RER,RRD,FIT,FUELDN,NXSPER,DELTAT(1,IP),MIXPWR,PFACT,
     3  SIG(1,1,1,IP),VPHV(1,IP))
        NLENGT=(NVAR+1)*(NREAC+1)*NBMIX
        CALL LCMPUT(IPDEPL,'MICRO-RATES',NLENGT,2,SIG(1,1,1,IP))
        CALL LCMPUT(IPDEPL,'INT-FLUX',NBMIX,2,VPHV(1,IP))
        VPHD=0.0
        DO 190 ICMB=1,NCOMB
        VPHD=VPHD+VPHV(MILVO(ICMB),IP)
  190   CONTINUE
        VPH(IP)=VPHD
        IF(INR.NE.0) THEN
           FNORM=VPH(IP)/REAL(VPHINI)
           IF(INR.EQ.3) CALL LCMPUT(IPDEPL,'FORM-POWER',1,2,PFACT)
        ELSE
           FNORM=0.0
        ENDIF
        CALL LCMPUT(IPDEPL,'FLUX-NORM',1,2,FNORM)
        IF((INDREC.EQ.1).AND.(IP.EQ.1)) THEN
          BRNWIR(1)=0.0
          BRNWIR(2)=0.0
          CALL LCMPUT(IPDEPL,'BURNUP-IRRAD',2,2,BRNWIR)
        ENDIF
        CALL LCMSIX(IPDEPL,' ',2)
  200   CONTINUE
      ELSE IF(ISAVE.EQ.1) THEN
        IP=1
        ITIM=0
        IF(NTIM.GT.0) THEN
          CALL XDRSET(TIMES,MAXTIM,0.0)
          CALL LCMGET(IPDEPL,'DEPL-TIMES',TIMES)
          DO 210 I=1,NTIM
            IF(ABS(TIMES(I)-XTI).LE.1.0E-4*XTI) ITIM=I
  210     CONTINUE
        ENDIF
        IF(ITIM.EQ.0) THEN
          IF(NTIM.GT.0) THEN
            IF(XTI.LT.TIMES(NTIM)) CALL XABORT('EVODRV: INVALID X1.')
          ENDIF
          NTIM=NTIM+1
          IF(NTIM.GT.MAXTIM) CALL XABORT('EVODRV: INVALID MAXTIM(2).')
          TIMES(NTIM)=XTI
          CALL LCMPUT(IPDEPL,'DEPL-TIMES',NTIM,2,TIMES)
          ITIM=NTIM
        ENDIF
        WRITE(TEXT12,'(8HDEPL-DAT,I4.4)') ITIM
        IF(IMPX.GT.0) WRITE(IUNOUT,530) XTI,XTI/8.64E-4,TEXT12
        CALL LCMSIX(IPDEPL,TEXT12,1)
        CALL LCMPUT(IPDEPL,'ISOTOPESDENS',NBISO,2,DEN)
*----
*  COMPUTE, NORMALIZE AND SAVE THE MICROSCOPIC REACTION RATES.
*----
        CALL EVOSIG(IMPX,INR,IGLOB,NGROUP,NBMIX,NBISO,NCOMB,ISONAM,
     1  IPISO,DEN,FLUMIX,VX,MILVO,JM,NVAR,NSUPS,NREAC,HREAC,IDR,
     2  RER,RRD,FIT,FUELDN,NXSPER,DELTAT(1,IP),MIXPWR,PFACT,
     3  SIG(1,1,1,IP),VPHV(1,IP))
        NLENGT=(NVAR+1)*(NREAC+1)*NBMIX
        CALL LCMPUT(IPDEPL,'MICRO-RATES',NLENGT,2,SIG(1,1,1,IP))
        CALL LCMPUT(IPDEPL,'INT-FLUX',NBMIX,2,VPHV(1,IP))
        VPHD=0.0
        DO 220 ICMB=1,NCOMB
        VPHD=VPHD+VPHV(MILVO(ICMB),IP)
  220   CONTINUE
        VPH(IP)=VPHD
        IF(INR.NE.0) THEN
           FNORM=VPH(IP)/REAL(VPHINI)
           IF(INR.EQ.3) CALL LCMPUT(IPDEPL,'FORM-POWER',1,2,PFACT)
        ELSE
           FNORM=0.0
        ENDIF
        CALL LCMPUT(IPDEPL,'FLUX-NORM',1,2,FNORM)
        IF((INDREC.EQ.1).AND.(IP.EQ.1)) THEN
          BRNWIR(1)=0.0
          BRNWIR(2)=0.0
          CALL LCMPUT(IPDEPL,'BURNUP-IRRAD',2,2,BRNWIR)
        ENDIF
        CALL LCMSIX(IPDEPL,' ',2)
      ENDIF
*
      IF(IDEPL.EQ.1) THEN
*----
*  PERFORM A DEPLETION CALCULATION BETWEEN TIMES XT(1) AND XT(2).
*----
        CALL XDRSET(TIMES,MAXTIM,0.0)
        CALL LCMGET(IPDEPL,'DEPL-TIMES',TIMES)
        IF(IMPX.GT.0) WRITE(IUNOUT,600) XT(1),XT(2)
        LCOOL=.TRUE.
        DO 300 IP=1,2
        ITIM=0
        DO 230 I=1,NTIM
        IF(ABS(TIMES(I)-XT(IP)).LE.1.0E-4*XT(IP)) ITIM=I
  230   CONTINUE
        IF(ITIM.GT.0) THEN
          WRITE(TEXT12,'(8HDEPL-DAT,I4.4)') ITIM
          IF(IMPX.GT.0) WRITE(IUNOUT,520) XT(IP),XT(IP)/8.64E-4,TEXT12
          CALL LCMSIX(IPDEPL,TEXT12,1)
          CALL LCMGET(IPDEPL,'ISOTOPESDENS',DEN)
          IF(IP.EQ.1) CALL LCMGET(IPDEPL,'BURNUP-IRRAD',BRNWIR)
          CALL LCMLEN(IPDEPL,'MICRO-RATES',LENGT,ITYLCM)
          IF((LENGT.GT.0).OR.(IP.EQ.1)) THEN
            CALL LCMGET(IPDEPL,'MICRO-RATES',SIG(1,1,1,IP))
            CALL LCMGET(IPDEPL,'INT-FLUX',VPHV(1,IP))
          ENDIF
          IF((IP.EQ.1).AND.(INR.EQ.3)) THEN
             CALL LCMGET(IPDEPL,'FORM-POWER',PFACT)
          ENDIF
          CALL LCMSIX(IPDEPL,' ',2)
        ELSE
          IF(IP.EQ.1) CALL XABORT('EVODRV: NO DEPLETION DATA STORED.')
          IF((IEXTR.EQ.1).AND.(NTIM.GE.2).AND.(INR.NE.0)) THEN
*            PERFORM MICRO REACTION RATE EXTRAPOLATION.
             ITIM=0
             DO 240 I=1,NTIM
             IF(ABS(TIMES(I)-XT(1)).LE.1.0E-4*XT(1)) ITIM=I
  240        CONTINUE
             IF(ITIM.EQ.0) THEN
                CALL XABORT('EVODRV: TABLE LOOKUP FAILURE.')
             ELSE IF(ITIM.GE.2) THEN
                NLENGT=(NVAR+1)*(NREAC+1)*NBMIX
                T1=TIMES(ITIM-1)
                T2=TIMES(ITIM)
                IF(T1.GE.T2) CALL XABORT ('EVODRV: ALGORITHM FAILURE.')
                FACT=(XT(2)-T1)/(T2-T1)
                WRITE(TEXT12,'(8HDEPL-DAT,I4.4)') ITIM-1
                CALL LCMSIX(IPDEPL,TEXT12,1)
                CALL LCMLEN(IPDEPL,'MICRO-RATES',LENGT,ITYLCM)
                IF(LENGT.EQ.NLENGT) THEN
                  CALL LCMGET(IPDEPL,'MICRO-RATES',SIG(1,1,1,0))
                  CALL LCMGET(IPDEPL,'INT-FLUX',VPHV(1,0))
                  DO 252 IBM=1,NBMIX
                  VPHV(IBM,2)=VPHV(IBM,1)+FACT*(VPHV(IBM,1)-VPHV(IBM,0))
                  DO 251 IQ=1,NREAC+1
                  DO 250 IS=1,NVAR+1
                  SIG(IS,IQ,IBM,2)=SIG(IS,IQ,IBM,1)+
     1               FACT*(SIG(IS,IQ,IBM,1)-SIG(IS,IQ,IBM,0))
  250             CONTINUE
  251             CONTINUE
  252             CONTINUE
                  IF(IMPX.GT.0) WRITE(IUNOUT,'(/18H EVODRV: USE EXTRA,
     1            45HPOLATED MICRO REACTION RATES AT END-OF-STAGE.)')
                ELSE
                  DO 262 IBM=1,NBMIX
                  VPHV(IBM,2)=VPHV(IBM,1)
                  DO 261 IQ=1,NREAC+1
                  DO 260 IS=1,NVAR+1
                  SIG(IS,IQ,IBM,2)=SIG(IS,IQ,IBM,1)
  260             CONTINUE
  261             CONTINUE
  262             CONTINUE
                  IF(IMPX.GT.0) WRITE(IUNOUT,'(/18H EVODRV: USE BEGIN,
     1            48HNNING-OF-STAGE MICRO REACTION RATES AT END-OF-ST,
     2            7HAGE(1).)')
                ENDIF
                CALL LCMSIX(IPDEPL,' ',2)
             ENDIF
          ELSE
             DO 267 IBM=1,NBMIX
             VPHV(IBM,2)=VPHV(IBM,1)
             DO 266 IQ=1,NREAC+1
             DO 265 IS=1,NVAR+1
             SIG(IS,IQ,IBM,2)=SIG(IS,IQ,IBM,1)
  265        CONTINUE
  266        CONTINUE
  267        CONTINUE
             IF(IMPX.GT.0) WRITE(IUNOUT,'(/23H EVODRV: USE BEGINNING-,
     1       46HOF-STAGE MICRO REACTION RATES AT END-OF-STAGE.)')
          ENDIF
        ENDIF
        VPHD=0.0
        DO 270 ICMB=1,NCOMB
        VPHD=VPHD+VPHV(MILVO(ICMB),IP)
  270   CONTINUE
        VPH(IP)=VPHD
*
        DO 285 ICMB=1,NCOMB
        IBM=MILVO(ICMB)
        DO 280 IS=1,NVAR
          IF(JM(IBM,IS).GT.0) THEN
             YDPL(IS,IP,ICMB)=DEN(JM(IBM,IS))
          ELSE
             YDPL(IS,IP,ICMB)=0.0
          ENDIF
  280   CONTINUE
  285   CONTINUE
*
        IF(INR.NE.0) THEN
*          CHECK FOR OUT-OF-CORE DEPLETION.
           DO 295 IBM=1,NBMIX
           DO 290 IU=1,NGROUP
           LCOOL=LCOOL.AND.(FLUMIX(IU,IBM).EQ.0.)
  290      CONTINUE
  295      CONTINUE
        ENDIF
  300   CONTINUE
*
        IF(LCOOL.AND.(FIT.NE.0.0)) CALL XABORT('EVODRV: NEUTRON FLUX I'
     1  //'S ZERO. UNABLE TO NORMALIZE.')
        IF(LCOOL.AND.(IMPX.GT.1)) THEN
           WRITE (IUNOUT,'(/31H EVODRV: OUT-OF-CORE DEPLETION.)')
        ELSE IF(IMPX.GT.1) THEN
           WRITE (IUNOUT,'(/27H EVODRV: IN-CORE DEPLETION.)')
           IF((FUELDN(3).GT.0.0).AND.(INR.EQ.3)) THEN
              WRITE(IUNOUT,'(/31H EVODRV: FUEL POWER NORMALISATI,
     1        10HON FACTOR=,1P,E12.5,29H MW/TONNE. OUT-OF-FUEL POWER ,
     2        12HFORM FACTOR=,E12.5)') FIT/FUELDN(3)/PFACT,PFACT
           ENDIF
        ENDIF
        INR2=INR
        IF(LCOOL) INR2=0
*----
*  PERFORM THE DEPLETION CALCULATION
*----
        CALL EVOBLD(IMPX,INR2,IGLOB,NBMIX,NBISO,NCOMB,ISONAM,YDPL,VX,
     1  MILVO,JM,NVAR,NDFP,NSUPS,NREAC,NPAR,NFISS,XT,EPS1,EPS2,EXPMAX,
     2  H1,ITYPE,IDIRAC,FIT,DELTA,ENERG,KPAR,BPAR,YIELD,IDR,RER,RRD,
     3  AWR,FUELDN,SIG(1,1,1,1),VPH,VPHV,MIXPWR,VTOTD,IEVOLB,KFISS,KPF)
*----
*  SAVE THE INITIAL SATURATED NUMBER DENSITIES IN THE DEPLETION TABLE
*----
        IF((ISAVE.GE.0).AND.(ISAT.EQ.1)) THEN
          DO 315 ICMB=1,NCOMB
            IBM=MILVO(ICMB)
            DO 310 IS=1,NVAR
              IF(JM(IBM,IS).GT.0) DEN(JM(IBM,IS))=YDPL(IS,1,ICMB)
  310       CONTINUE
  315     CONTINUE
          ITIM=0
          DO 320 I=1,NTIM
            IF(ABS(TIMES(I)-XT(1)).LE.1.0E-4*XT(1)) ITIM=I
  320     CONTINUE
          IF(ITIM.EQ.0) CALL XABORT('EVODRV: MISSING TIME ENTRY.')
          WRITE(TEXT12,'(8HDEPL-DAT,I4.4)') ITIM
          IF(IMPX.GT.0) WRITE(IUNOUT,530) XT(1),XT(1)/8.64E-4,TEXT12
          CALL LCMSIX(IPDEPL,TEXT12,1)
          CALL LCMPUT(IPDEPL,'ISOTOPESDENS',NBISO,2,DEN)
          CALL LCMSIX(IPDEPL,' ',2)
        ENDIF
*----
*  SAVE THE DEPLETION CALCULATION RESULT IN THE DEPLETION TABLE
*----
        DO 335 ICMB=1,NCOMB
          IBM=MILVO(ICMB)
          DO 330 IS=1,NVAR
            IF(JM(IBM,IS).GT.0) DEN(JM(IBM,IS))=YDPL(IS,2,ICMB)
  330     CONTINUE
  335   CONTINUE
        ITIM=0
        DO 340 I=1,NTIM
          IF(ABS(TIMES(I)-XT(2)).LE.1.0E-4*XT(2)) ITIM=I
  340   CONTINUE
        IF(ITIM.EQ.0) THEN
          IF(XT(2).LT.TIMES(NTIM)) CALL XABORT('EVODRV: INVALID X2')
          NTIM=NTIM+1
          IF(NTIM.GT.MAXTIM) CALL XABORT('EVODRV: INVALID MAXTIM(3).')
          TIMES(NTIM)=XT(2)
          CALL LCMPUT(IPDEPL,'DEPL-TIMES',NTIM,2,TIMES)
          ITIM=NTIM
        ENDIF
        WRITE(TEXT12,'(8HDEPL-DAT,I4.4)') ITIM
        IF(IMPX.GT.0) WRITE(IUNOUT,530) XT(2),XT(2)/8.64E-4,TEXT12
        CALL LCMSIX(IPDEPL,TEXT12,1)
        CALL LCMPUT(IPDEPL,'ISOTOPESDENS',NBISO,2,DEN)
        NLENGT=(NVAR+1)*(NREAC+1)*NBMIX
        CALL LCMPUT(IPDEPL,'MICRO-RATES',NLENGT,2,SIG(1,1,1,2))
        CALL LCMPUT(IPDEPL,'INT-FLUX',NBMIX,2,VPHV(1,2))
        CALL LCMPUT(IPDEPL,'ENERG-MIX',NBMIX,2,ENERG)
*       We use DELTA(3) instead of DELTA(1) in order to avoid different
*       base points in multi-D tables.
        BRNWIR(1)=BRNWIR(1)+DELTA(3)/8.64E-4
        BRNWIR(2)=BRNWIR(2)+DELTA(2)
        CALL LCMPUT(IPDEPL,'BURNUP-IRRAD',2,2,BRNWIR)
        IF(IMPX.GE.1) WRITE(IUNOUT,580) XT(2)/8.64E-4,
     >                BRNWIR(1),BRNWIR(2)
        CALL LCMSIX(IPDEPL,' ',2)
      ENDIF
*----
*  RELEASE THE ALLOCATED MEMORY
*----
      DEALLOCATE(IDR,HREAC,KPAR)
      DEALLOCATE(RRD,RER,YIELD,BPAR)
      DEALLOCATE(KPF,KFISS)
*----
*  USE THE RESULT OF A DEPLETION CALCULATION IN THE FOLLOWING RUN
*----
      IF(ISET.GE.0) THEN
         CALL XDRSET(TIMES,MAXTIM,0.0)
         CALL LCMGET(IPDEPL,'DEPL-TIMES',TIMES)
         ITIM=0
         DO 350 I=1,NTIM
           IF(ABS(TIMES(I)-XTF).LE.1.0E-4*XTF) ITIM=I
  350    CONTINUE
         IF(ITIM.EQ.0) CALL XABORT('EVODRV: NO DEPLETION DATA STORED.')
         WRITE(TEXT12,'(8HDEPL-DAT,I4.4)') ITIM
         IF(IMPX.GT.0) WRITE(IUNOUT,520) XTF,XTF/8.64E-4,TEXT12
         CALL LCMSIX(IPDEPL,TEXT12,1)
         CALL LCMGET(IPDEPL,'ISOTOPESDENS',DEN)
         CALL LCMGET(IPDEPL,'BURNUP-IRRAD',BRNWIR)
         CALL LCMSIX(IPDEPL,' ',2)
         IF(IMPX.GT.1) THEN
           WRITE(IUNOUT,550) XTF,XTF/8.64E-4
           DO 370 ICMB=1,NCOMB
             IMIXC=MILVO(ICMB)
             NISOCC=0
             DO 360 ISOT=1,NBISO
               IMIXI=MIX(ISOT)
               IF(IMIXI.EQ.IMIXC) THEN
                 NISOCC=NISOCC+1
                 ISOCMB(NISOCC)=ISOT
               ENDIF
 360         CONTINUE
             WRITE(IUNOUT,560) IMIXC
             WRITE(IUNOUT,570) ((ISONAM(I0,ISOCMB(I)),I0=1,2),
     1       DEN(ISOCMB(I)),I=1,NISOCC)
 370       CONTINUE
         ENDIF
      ENDIF
*----
*  RECOVER THE BURNUP AND SAVE IT IN A CLE-2000 VARIABLE
*----
      IF(IPICK.EQ.1) THEN
         CALL XDRSET(TIMES,MAXTIM,0.0)
         CALL LCMGET(IPDEPL,'DEPL-TIMES',TIMES)
         ITIM=0
         DO 375 I=1,NTIM
           IF(ABS(TIMES(I)-XTF).LE.1.0E-4*XTF) ITIM=I
  375    CONTINUE
         IF(ITIM.EQ.0) CALL XABORT('EVODRV: NO DEPLETION DATA STORED.')
         WRITE(TEXT12,'(8HDEPL-DAT,I4.4)') ITIM
         IF(IMPX.GT.0) WRITE(IUNOUT,520) XTF,XTF/8.64E-4,TEXT12
         CALL LCMSIX(IPDEPL,TEXT12,1)
         CALL LCMGET(IPDEPL,'BURNUP-IRRAD',BRNWIR)
         CALL LCMSIX(IPDEPL,' ',2)
         CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT12,DBLLIR)
         IF(ITYPLU.NE.-2) CALL XABORT('EVODRV: OUTPUT REAL EXPECTED.')
         ITYPLU=2
         REALIR=BRNWIR(1)
         IF(IMPX.GT.2) WRITE(IUNOUT,540) REALIR
         CALL REDPUT(ITYPLU,INTLIR,REALIR,TEXT12,DBLLIR)
         CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT12,DBLLIR)
         IF((ITYPLU.NE.3).OR.(TEXT12.NE.';')) THEN
           CALL XABORT('EVODRV: ; CHARACTER EXPECTED.')
         ENDIF      
      ENDIF      
*
      IF((IDEPL.EQ.1).OR.(ISET.EQ.1).AND.LMACRO) THEN
*       COMPUTE THE NEW DEPLETED MACROSCOPIC CROSS SECTIONS.
        CALL XDLSET(MASKL,NGROUP,.TRUE.)
        CALL XDLSET(MASK,NBMIX,.FALSE.)
        DO 380 ICOMB=1,NCOMB
        MASK(MILVO(ICOMB))=.TRUE.
 380    CONTINUE
*
        ITSTMP=2
        TMPDAY(1)=XT(2)/8.64E-4
        TMPDAY(2)=BRNWIR(1)
        TMPDAY(3)=BRNWIR(2)
        CALL LIBMIX(IPLIB,NBMIX,NGROUP,NBISO,ISONAM,MIX,DEN,MASK,MASKL,
     1  ITSTMP,TMPDAY)
      ENDIF
*----
*  STORE THE GENERAL DEPLETION RELATED PARAMETERS
*----
      CALL LCMPUT(IPDEPL,'VOLUME-MIX',NBMIX,2,VX)
      CALL LCMPUT(IPDEPL,'DEPLETE-MIX',NVAR*NBMIX,1,JM)
      CALL LCMPUT(IPDEPL,'MIXTURESBurn',NBMIX,1,MIXBRN)
      CALL LCMPUT(IPDEPL,'MIXTURESPowr',NBMIX,1,MIXPWR)
      CALL XDISET(IPAR,NSTATE,0)
      IPAR(1)=ITYPE
      IPAR(2)=INR
      IPAR(3)=NTIM
      IPAR(4)=NBISO
      IPAR(5)=NCOMB
      IPAR(6)=NREAC
      IPAR(7)=NVAR
      IPAR(8)=NBMIX
      IPAR(9)=IEXTR
      IPAR(10)=IGLOB
      IPAR(11)=ISAT
      IPAR(12)=IDIRAC
      IPAR(13)=ITIXS
      IPAR(14)=IFLMAC
      IPAR(15)=IYLMIX
      CALL LCMPUT(IPDEPL,'STATE-VECTOR',NSTATE,1,IPAR)
      RPAR(1)=EPS1
      RPAR(2)=EPS2
      RPAR(3)=EXPMAX
      RPAR(4)=H1
      RPAR(5)=FIT
      CALL LCMPUT(IPDEPL,'EVOLUTION-R',5,2,RPAR)
      IF((IMPX.GT.1).OR.((IMPX.GT.0).AND.(INDREC.EQ.1))) THEN
         WRITE(IUNOUT,590) IMPX,ITYPE,INR,NTIM,NBISO,NCOMB,NREAC,
     1   NVAR,NBMIX,IEXTR
         WRITE(IUNOUT,595) IGLOB,ISAT,IDIRAC,ITIXS,IFLMAC,IYLMIX
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IPISO)
      DEALLOCATE(MASKL,MASK)
      DEALLOCATE(YDPL,AWR,ENERG,VPHV,SIG)
      DEALLOCATE(IEVOLB,INADPL,ISOCMB,MILVO,JM)
      RETURN
*
  500 FORMAT(/36H EVODRV: NUMBER OF FISSILE ISOTOPES=,I4/9X,8HNUMBER O,
     1 19HF FISSION PRODUCTS=,I4)
  510 FORMAT(/24H EVODRV: UNABLE TO FIND ,A15,6H INDEX,I5,6H AMONG,10I5
     1 /(56X,10I5))
  520 FORMAT(/44H EVODRV: RECOVER INFORMATION RELATED TO TIME,1P,E12.4,
     1 8H E+8 S (,E12.4,32H DAY) FROM LCM DIRECTORY NAMED ',A12,2H'.)
  530 FORMAT(/41H EVODRV: SAVE INFORMATION RELATED TO TIME,1P,E12.4,
     1 8H E+8 S (,E12.4,30H DAY) ON LCM DIRECTORY NAMED ',A12,2H'.)
  540 FORMAT(/21H EVODRV: PICK BURNUP=,1P,E12.4,10H MWd/tonne)
  550 FORMAT(/' EVODRV: NUMBER DENSITIES PER ISOTOPE AT TIME',1P,
     1 E12.4,' E+8 S (',E12.4,' DAY)')
  560 FORMAT(/' ISOTOPIC DENSITIES AFTER BURNUP FOR MIXTURE = ',I5,
     1 22H (10**24 PARTICLES/CC))
  570 FORMAT(1P,5(4X,2A4,':',E12.4))
  580 FORMAT(' -> FINAL BURNUP AT TIME  = ',1P,E14.6,' DAYS'/
     >       '    FUEL BURNUP           = ',E14.6,' MW*D/TONNE'/
     >       '    NEUTRON EXPOSURE      = ',E14.6,' N/KB')
  590 FORMAT(/8H OPTIONS/8H -------/
     1 7H IMPX  ,I6,30H   (0=NO PRINT/1=SHORT/2=MORE)/
     2 7H ITYPE ,I6,31H   (1=CASH-KARP/2=KAPS-RENTROP)/
     3 7H INR   ,I6,47H   (0=OUT-OF-CORE/1=CONSTANT FLUX/2=CONSTANT PO,
     4 4HWER)/
     5 7H NTIM  ,I6,35H   (NUMBER OF DEPLETION SET POINTS)/
     6 7H NBISO ,I6,29H   (TOTAL NUMBER OF ISOTOPES)/
     7 7H NCOMB ,I6,33H   (NUMBER OF DEPLETING MIXTURES)/
     8 7H NREAC ,I6,34H   (NUMBER OF DEPLETING REACTIONS)/
     9 7H NVAR  ,I6,33H   (NUMBER OF DEPLETING ISOTOPES)/
     1 7H NBMIX ,I6,23H   (NUMBER OF MIXTURES)/
     2 7H IEXTR ,I6,47H   (0=NO FLUX EXTRAPOLATION/1=FLUX EXTRAPOLATIO,
     3 2HN))
  595 FORMAT(
     1 7H IGLOB ,I6,47H   (0=COMPUTE BURNUP IN FUEL/1=COMPUTE BURNUP I,
     2 14HN GLOBAL CELL)/
     3 7H ISAT  ,I6,47H   (0/1=DO NOT/DO SAVE SATURATED INITIAL NUMBER,
     4 11H DENSITIES)/
     5 7H IDIRAC,I6,47H   (0/1=DO NOT/DO USE DIRAC FUNCTION CONTRIBUTI,
     6 34HONS IN SATURATED NUMBER DENSITIES)/
     7 7H ITIXS ,I6,38H   (0/1=TIME-DEPENDENT XS FLAG ON/OFF)/
     8 7H IFLMAC,I6,47H   (0/1/2=RECOVER FLUX FROM L_FLUX/L_MACROLIB/L,
     9 7H_POWER)/
     1 7H IYLMIX,I6,47H   (0/1=RECOVER FISSION YIELD DATA FROM DEPL-CH,
     2 16HAIN/PYIELD DATA))
  600 FORMAT(/54H EVODRV: SOLUTION OF A DEPLETION SYSTEM BETWEEN TIMES ,
     1 1P,E12.4,4H AND,E12.4,6H E+8 S)
  610 FORMAT(/' EVODRV:  FUEL INITIAL DENSITY   = ',1P,E14.6,' G/CC'/
     >        '          FUEL TOTAL VOLUME      = ',E14.6,' CC'/
     >        '          FUEL INITIAL MASS      = ',E14.6,' G'/
     >        '          FUEL INITIAL MASS/CELL VOL = ',E14.6,' G/CC')
      END
