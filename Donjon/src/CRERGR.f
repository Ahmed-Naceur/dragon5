*DECK CRERGR
      SUBROUTINE CRERGR(IPCPO,IPMAP,NISO,NGRP,NMIXT,NL,IBM,IMPX,IBTYP,
     1 DERIV,UPS,NBURN,BURNUP,ILEAK,TOTAL,ZNUG,SNUGF,CHI,OVERV,DIFFX,
     2 DIFFY,DIFFZ,H,SCAT,IJJ,NJJ,HISO,ITY,CONC,FMIX,BRN0,BRN1,NCH,NB,
     3 IVARTY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform interpolation of fuel properties over the fuel lattice.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert, D. Sekki
*
*Parameters: input
* IPCPO   pointer to L_COMPO information.
* IPMAP   pointer to L_MAP information.
* NISO    1+number of extracted isotopes.
* NGRP    number of energy groups.
* NMIXT   number of material mixtures in the fuel-map macrolib.
* NL      number of legendre orders (=1 for isotropic scattering).
* IBM     mixture number to be treat.
* IMPX    printing index (=0 for no print).
* IBTYP   type of interpolation: =1 time-average; =2 instantaneous;
*         derivative with respect to a single exit burnup.
* DERIV   =.true.: derivative of macrolib info is computed with
*         respect to burn1.
* UPS     =.true.: no upscatering cross sections will be stored.
* NBURN   number of tabulated burnup steps.
* BURNUP  burnup tabulated values from compo file.
* HISO    hollerith name information for extracted isotopes.
* ITY     =0: do not process the isotope; =1: use number density
*         stored in conc(i); =2: use number density stored in compo.
* CONC    user defined number density.
* NCH     number of reactor channels.
* NB      number of fuel bundles per channel.
* FMIX    fuel mixture indices per fuel bundle.
* BRN0    contains either low burnup integration limits or 
*         instantaneous burnups per fuel bundle.
* BRN1    upper burnup integration limits per fuel bundle.
* IVARTY  index of the exit burnup used to compute derivatives. Used
*         if IBTYP=3.
*
*Parameters: output
* ILEAK   diffusion coefficient flag (=1: isotropic; =2: anisotropic).
* TOTAL   total macroscopic x-sections.
* ZNUG    nu*fission macroscopic x-sections.
* SNUGF   fission macroscopic x-sections.
* CHI     fission spectrum.
* OVERV   reciprocal neutron velocities.
* DIFFX   x-directed diffusion coefficients.
* DIFFY   y-directed diffusion coefficients.
* DIFFZ   z-directed diffusion coefficients.
* H       h-factors (kappa*fission macroscopic x-sections).
* SCAT    scattering macroscopic x-sections.
*
*Parameters: 
* IJJ
* NJJ
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPCPO,IPMAP
      INTEGER NISO,IBTYP,IBM,NMIXT,NBURN,NGRP,NL,IMPX,NCH,NB,ILEAK,
     1     IJJ(NMIXT,NL,NGRP),NJJ(NMIXT,NL,NGRP),FMIX(NCH*NB),
     2     HISO(3*NISO),ITY(NISO),IVARTY
      REAL CONC(NISO),TOTAL(NMIXT,NGRP),BURNUP(NBURN),SNUGF(NMIXT,NGRP),
     1     CHI(NMIXT,NGRP),OVERV(NMIXT,NGRP),DIFFX(NMIXT,NGRP),
     2     DIFFY(NMIXT,NGRP),DIFFZ(NMIXT,NGRP),BRN0(NCH*NB),
     3     BRN1(NCH*NB),H(NMIXT,NGRP),SCAT(NMIXT,NL,NGRP,NGRP),
     4     ZNUG(NMIXT,NGRP)
      LOGICAL DERIV,UPS
*----
*  LOCAL VARIABLES
*----
      LOGICAL LCUBIC
      PARAMETER(LCUBIC=.TRUE.)
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ZONEDP
      REAL, ALLOCATABLE, DIMENSION(:) :: TERP,TERPW
      REAL, ALLOCATABLE, DIMENSION(:) :: YTOTAL,YZNUG,YNUGF,YCHI,YOVERV,
     1 YDIFX,YDIFY,YDIFZ,YH,YSCAT,YFLUX
      REAL, ALLOCATABLE, DIMENSION(:,:) :: ZTOTAL,ZZNUG,ZNUGF,ZCHI,
     1 ZOVERV,ZDIFX,ZDIFY,ZDIFZ,ZH,ZSCAT,ZFLUX
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(TERP(NBURN),ZONEDP(NCH,NB),TERPW(NBURN))
*
      CALL XDRSET(BURNUP,NBURN,0.)
      CALL LCMGET(IPCPO,'BURNUP',BURNUP)
*----
*  FUEL-MAP INFORMATION
*----
      CALL CREGET(IPMAP,NCH,NB,IBTYP,IMPX,BRN0,BRN1,FMIX,ZONEDP,
     1 IVARTY,VARVAL)
*----
*  CREATE BURNUP-DEPENDENT TABLE
*----
      ALLOCATE(YTOTAL(NGRP),YZNUG(NGRP),YNUGF(NGRP),YCHI(NGRP),
     1 YOVERV(NGRP),YDIFX(NGRP),YDIFY(NGRP),YDIFZ(NGRP),YH(NGRP),
     2 YSCAT(NL*NGRP*NGRP),YFLUX(NGRP))
*
      CALL XDRSET(YTOTAL,NGRP,0.)
      CALL XDRSET(YZNUG,NGRP,0.)
      CALL XDRSET(YNUGF,NGRP,0.)
      CALL XDRSET(YCHI,NGRP,0.)
      CALL XDRSET(YOVERV,NGRP,0.)
      CALL XDRSET(YDIFX,NGRP,0.)
      CALL XDRSET(YDIFY,NGRP,0.)
      CALL XDRSET(YDIFZ,NGRP,0.)
      CALL XDRSET(YH,NGRP,0.)
      CALL XDRSET(YSCAT,NL*NGRP*NGRP,0.)
      CALL XDRSET(YFLUX,NGRP,0.)
*
      ALLOCATE(ZTOTAL(NGRP,NBURN),ZZNUG(NGRP,NBURN),ZNUGF(NGRP,NBURN),
     1 ZCHI(NGRP,NBURN),ZOVERV(NGRP,NBURN),ZDIFX(NGRP,NBURN),
     2 ZDIFY(NGRP,NBURN),ZDIFZ(NGRP,NBURN),ZH(NGRP,NBURN),
     3 ZSCAT(NL*NGRP*NGRP,NBURN),ZFLUX(NGRP,NBURN))
*
      CALL XDRSET(ZTOTAL,NGRP*NBURN,0.)
      CALL XDRSET(ZZNUG,NGRP*NBURN,0.)
      CALL XDRSET(ZNUGF,NGRP*NBURN,0.)
      CALL XDRSET(ZCHI,NGRP*NBURN,0.)
      CALL XDRSET(ZOVERV,NGRP*NBURN,0.)
      CALL XDRSET(ZDIFX,NGRP*NBURN,0.)
      CALL XDRSET(ZDIFY,NGRP*NBURN,0.)
      CALL XDRSET(ZDIFZ,NGRP*NBURN,0.)
      CALL XDRSET(ZH,NGRP*NBURN,0.)
      CALL XDRSET(ZSCAT,NL*NGRP*NGRP*NBURN,0.)
      CALL XDRSET(ZFLUX,NGRP*NBURN,0.)
*
      CALL CRETAB(IPCPO,NISO,NGRP,NL,IMPX,HISO,NBURN,ITY,CONC,ILEAK,
     1 ZTOTAL,ZZNUG,ZNUGF,ZCHI,ZOVERV,ZDIFX,ZDIFY,ZDIFZ,ZH,ZSCAT,ZFLUX,
     2 UPS)
*----
*  PERFORM INTERPOLATION
*----
      DO 105 ICH=1,NCH
      DO 100 J=1,NB
      IB=(J-1)*NCH+ICH
      IF(FMIX(IB).EQ.IBM)THEN
        IF(IBTYP.EQ.1)THEN
*         TIME-AVERAGE
          BURN0=BRN0(IB)
          BURN1=BRN1(IB)
          IF(BURN0.GE.BURN1) CALL XABORT('@CRERGR: INVALID BURNUP LIMI'
     1    //'TS(1).')
          CALL ALTERI(LCUBIC,NBURN,BURNUP,BURN0,BURN1,TERP)
          DO 20 I=1,NBURN
          TERP(I)=TERP(I)/(BURN1-BURN0)
   20     CONTINUE
        ELSEIF(IBTYP.EQ.2)THEN
*         INSTANTANEOUS
          BURN0=BRN0(IB)
          BURN1=BURN0
          IF(NBURN.EQ.1) THEN
            TERP(1)=1.0
          ELSE
            CALL ALTERP(LCUBIC,NBURN,BURNUP,BURN0,DERIV,TERP)
          ENDIF
        ELSEIF(IBTYP.EQ.3)THEN
*         DERIVATIVE WITH RESPECT TO A SINGLE EXIT BURNUP. USE EQ.(3.3)
*         OF RICHARD CHAMBON'S THESIS.
          IF(ZONEDP(ICH,J).NE.0) THEN
            BURN0=BRN0(IB)
            BURN1=BRN1(IB)
            IF(BURN0.GE.BURN1) CALL XABORT('@CRERGR: INVALID BURNUP LI'
     1      //'MITS(2).')
            CALL ALTERI(LCUBIC,NBURN,BURNUP,BURN0,BURN1,TERPW)
            DO 30 I=1,NBURN
            TERP(I)=-TERPW(I)
   30       CONTINUE
            CALL ALTERP(LCUBIC,NBURN,BURNUP,BURN0,.FALSE.,TERPW)
            DO 40 I=1,NBURN
            TERP(I)=TERP(I)-TERPW(I)*BURN0
   40       CONTINUE
            CALL ALTERP(LCUBIC,NBURN,BURNUP,BURN1,.FALSE.,TERPW)
            DO 50 I=1,NBURN
            TERP(I)=(TERP(I)+TERPW(I)*BURN1)/(VARVAL*(BURN1-BURN0))
   50       CONTINUE
          ELSE
            CALL XDRSET(TERP,NBURN,0.0)
          ENDIF
        ENDIF
        IF(BURN1.GT.BURNUP(NBURN))THEN
          WRITE(*,*)'@CRERGR: BURN1 VALUE  :',BURN1
          WRITE(*,*)'@CRERGR: BURNUP LIMIT :',BURNUP(NBURN)
          CALL XABORT('@CRERGR: INTERPOLATION IS OUT OF BURNUP LIMIT.')
        ENDIF
*
        IF((IBTYP.EQ.3).AND.(ZONEDP(ICH,J).EQ.0)) THEN
           CALL XDRSET(YTOTAL,NGRP,0.)
           CALL XDRSET(YZNUG,NGRP,0.)
           CALL XDRSET(YNUGF,NGRP,0.)
           CALL XDRSET(YCHI,NGRP,0.)
           CALL XDRSET(YOVERV,NGRP,0.)
           CALL XDRSET(YDIFX,NGRP,0.)
           CALL XDRSET(YDIFY,NGRP,0.)
           CALL XDRSET(YDIFZ,NGRP,0.)
           CALL XDRSET(YH,NGRP,0.)
           CALL XDRSET(YSCAT,NL*NGRP*NGRP,0.)
           CALL XDRSET(YFLUX,NGRP,0.)
        ELSE
           CALL CREITP(NGRP,NL,NBURN,TERP,YTOTAL,YZNUG,YNUGF,YCHI,
     1     YOVERV,YDIFX,YDIFY,YDIFZ,YH,YSCAT,YFLUX,ZTOTAL,ZZNUG,ZNUGF,
     2     ZCHI,ZOVERV,ZDIFX,ZDIFY,ZDIFZ,ZH,ZSCAT,ZFLUX)
        ENDIF
*       DATA STORAGE
        DO 72 JGR=1,NGRP
        TOTAL(IB,JGR)=YTOTAL(JGR)
        ZNUG(IB,JGR)=YZNUG(JGR)
        SNUGF(IB,JGR)=YNUGF(JGR)
        CHI(IB,JGR)=YCHI(JGR)
        OVERV(IB,JGR)=YOVERV(JGR)
        DIFFX(IB,JGR)=YDIFX(JGR)
        DIFFY(IB,JGR)=YDIFY(JGR)
        DIFFZ(IB,JGR)=YDIFZ(JGR)
        H(IB,JGR)=YH(JGR)
        DO 71 IGR=1,NGRP
        DO 70 IL=1,NL
        SCAT(IB,IL,IGR,JGR)=YSCAT(NL*((JGR-1)*NGRP+IGR-1)+IL)
   70   CONTINUE
   71   CONTINUE
   72   CONTINUE
*       JGR IS THE SECONDARY GROUP.
        DO 85 JGR=1,NGRP
        DO 80 IL=1,NL
        IGMIN=JGR
        IGMAX=JGR
        DO IGR=NGRP,1,-1
          IF(SCAT(IB,IL,IGR,JGR).NE.0.)THEN
            IGMIN=MIN(IGMIN,IGR)
            IGMAX=MAX(IGMAX,IGR)
          ENDIF
        ENDDO
        IJJ(IB,IL,JGR)=IGMAX
        NJJ(IB,IL,JGR)=IGMAX-IGMIN+1
   80   CONTINUE
   85   CONTINUE
      ENDIF
  100 CONTINUE
  105 CONTINUE
*
      DEALLOCATE(YFLUX,YSCAT,YH,YDIFZ,YDIFY,YDIFX,YOVERV,YCHI,YNUGF,
     1 YZNUG,YTOTAL)
*
      DEALLOCATE(ZFLUX,ZSCAT,ZH,ZDIFZ,ZDIFY,ZDIFX,ZOVERV,ZCHI,ZNUGF,
     1 ZZNUG,ZTOTAL)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(TERPW,ZONEDP,TERP)
      RETURN
      END
