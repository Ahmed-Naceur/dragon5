*DECK XENCAL
      SUBROUTINE XENCAL(IPLIB,IPPOW,NB,NCH,NGRP,NMIX,NBISO,XEN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the Xenon distribution according to the bundle flux
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
*
*Author(s): 
* M. Guyot
*
*Parameters: input/output
* IPLIB   adress of the L_LIBRARY
* IPPOW   adress of the L_POWER
* NB      number of fuel bundles per channel
* NCH     number of channels
* NGRP    number of energy groups
* NMIX    number of mixtures present in the library
* NBISO   number of isotopes
* XEN     xenon concentrations in each bundle
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPPOW
      INTEGER NB,NCH,NGRP,NMIX
      REAL XEN(NMIX)
*----
*  LOCAL VARIABLES
*----
      INTEGER ICH,IB,IGRP,IBM
      REAL TAUF(NMIX,NGRP),TAUX(NMIX,NGRP),XLAMBDA,GAMMAI,GAMMAX,CF,TF,
     1 TX
      TYPE(C_PTR) JPLIB,KPLIB,LPLIB,MPLIB
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IMIX
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGX,SIGF
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: FLUB
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: HNAMIS
*----
*  SCRATCH STORAGE ALLOCATION
*  SIGX    microscopic capture cross section of Xe-135
*  SIGF    macroscopic fission cross section
*  FLUB    bundle fluxes
*----
      ALLOCATE(SIGX(NGRP),SIGF(NMIX),FLUB(NCH,NB,NGRP))
      ALLOCATE(HNAMIS(NBISO),IMIX(NBISO))
*----
* SET THE YIELD AND THE DECAY CONSTANTE FOR XENON
*----
      XLAMBDA = 2.09E-5
      GAMMAI = 0.0631
      GAMMAX = 0.0045
      CF=1.0E-24
*----
*  COMPUTE FISSION AND XENON REACTION RATES IN EACH BUNDLE
*---- 
      CALL XDRSET(FLUB,NCH*NB*NGRP,0.0)
      CALL LCMGET(IPPOW,'FLUX-BUND',FLUB)
      CALL XDRSET(TAUX,NMIX*NGRP,0.0)
      CALL XDRSET(TAUF,NMIX*NGRP,0.0)
      CALL LCMSIX(IPLIB,'MACROLIB',1)
      JPLIB=LCMGID(IPLIB,'GROUP')
      CALL LCMSIX(IPLIB,' ',2)
      CALL LCMGTC(IPLIB,'ISOTOPESUSED',12,NBISO,HNAMIS)
      CALL LCMGET(IPLIB,'ISOTOPESMIX',IMIX)
      LPLIB=LCMGID(IPLIB,'ISOTOPESLIST')
*
      DO 10 IGRP=1,NGRP 
        KPLIB=LCMGIL(JPLIB,IGRP)
        CALL LCMGET(KPLIB,'NFTOT',SIGF)
        DO 20 ICH=1,NCH
          DO 30 IB=1,NB
            IBM=NCH*(IB-1)+ICH
            ISO=0
            DO JSO=1,NBISO
              IF((HNAMIS(JSO).EQ.'Xe135').AND.(IMIX(JSO).EQ.IBM)) THEN
                ISO=JSO
                GO TO 35
              ENDIF
            ENDDO
            CALL XABORT('XENCAL: UNABLE TO FIND ISOTOPE=Xe135.')
   35       MPLIB=LCMGIL(LPLIB,ISO)
            CALL LCMGET(MPLIB,'NG',SIGX)
            TAUX(IBM,IGRP)=TAUX(IBM,IGRP)+FLUB(ICH,IB,IGRP)*
     +        SIGX(IGRP)
            TAUF(IBM,IGRP)=TAUF(IBM,IGRP)+FLUB(ICH,IB,IGRP)*
     +        SIGF(IBM)
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
*
      DO 40 IBM=1,NMIX
        TF=0.0
        TX=0.0
        DO 50 IGRP=1,NGRP
          TF=TF+TAUF(IBM,IGRP)
          TX=TX+TAUX(IBM,IGRP)
   50   CONTINUE
        XEN(IBM)=CF*(GAMMAX+GAMMAI)*TF/(XLAMBDA+TX*CF)
   40 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IMIX,HNAMIS)
      DEALLOCATE(FLUB,SIGF,SIGX)
      RETURN
      END
