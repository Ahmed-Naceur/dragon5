*DECK CREXSI
      SUBROUTINE CREXSI(IPMAP,NENTRY,HENTRY,KENTRY,NMIX,NGRP,NL,ILEAK,
     1 IMPX,TOTAL,ZNUG,SNUGF,CHI,OVERV,DIFFX,DIFFY,DIFFZ,H,IJJ,NJJ,SCAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover and/or interpolate l_compo data.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert, D. Sekki
*
*Parameters: input
* IPMAP   pointer to the fuel-map information.
* NENTRY  number of lcm or xsm objects used by the module.
* HENTRY  character*12 name of each lcm or xsm objects.
* KENTRY  pointers to the lcm or xsm objects.
* NMIX    maximum number of material mixtures.
* NGRP    number of energy groups.
* NL      number of legendre orders (=1 for isotropic scattering).
* ILEAK   diffusion coefficient flag (=1: isotropic; =2: anisotropic).
* IMPX    printing index (=0 for no print).
*
*Parameters: output
* TOTAL   total macroscopic x-sections.
* ZNUG    nu*fission macroscopic x-sections.
* SNUGF   fission macroscopic x-sections.
* CHI     fission spectrum.
* OVERV   reciprocal neutron velocities.
* DIFFX   x-directed diffusion coefficients.
* DIFFY   y-directed diffusion coefficients.
* DIFFZ   z-directed diffusion coefficients.
* H       h-factors (kappa*fission macroscopic x-sections).
* IJJ     profile storage index.
* NJJ     profile storage width.
* SCAT    scattering macroscopic x-sections.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NENTRY,NMIX,NGRP,NL,ILEAK,IMPX,IJJ(NMIX,NL,NGRP),
     1 NJJ(NMIX,NL,NGRP)
      TYPE(C_PTR) IPMAP,KENTRY(NENTRY)
      CHARACTER HENTRY(NENTRY)*12
      REAL TOTAL(NMIX,NGRP),ZNUG(NMIX,NGRP),SNUGF(NMIX,NGRP),
     1     CHI(NMIX,NGRP),OVERV(NMIX,NGRP),DIFFX(NMIX,NGRP),
     2     DIFFY(NMIX,NGRP),DIFFZ(NMIX,NGRP),H(NMIX,NGRP),
     3     SCAT(NMIX,NL,NGRP,NGRP)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPCPO,JPCPO,JPMAP,KPMAP
      PARAMETER(NSTATE=40,IOUT=6)
      CHARACTER TEXT*12,NAMDIR*12,HCOMPO*12,HSMG*131
      INTEGER IPAR(NSTATE),IDATA(NSTATE)
      LOGICAL DERIV,UPS,LTAB
      DOUBLE PRECISION DFLOT
      REAL, ALLOCATABLE, DIMENSION(:) :: YTOTAL,YZNUG,YNUGF,YCHI,YOVERV,
     1 YDIFX,YDIFY,YDIFZ,YH,YSCAT,YFLUX
      INTEGER, ALLOCATABLE, DIMENSION(:) :: HISO,ITY,FMIX
      REAL, ALLOCATABLE, DIMENSION(:) :: CONC,BURNU,BRN0,BRN1
*
      IVARTY=0
      UPS=.FALSE.
      LTAB=.FALSE.
      DERIV=.FALSE.
      NFUEL=0
      MAXEN=NENTRY
      IF(C_ASSOCIATED(IPMAP))THEN
        CALL LCMGET(IPMAP,'STATE-VECTOR',IDATA)
        IF(IDATA(4).NE.NGRP)CALL XABORT('@CREXSI: DIFFERENT NUM'
     1   //'BER OF ENERGY GROUPS IN COMPO AND FUEL MAP.')
        NB=IDATA(1)
        NCH=IDATA(2)
        NFUEL=IDATA(7)
        MAXEN=MAXEN-1
        LTAB=.TRUE.
      ENDIF
*----
*  READ INTERPOLATION OPTION
*----
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      DO 200 IEN=2,MAXEN
*     KEYWORD COMPO OR TABLE
      IF(TEXT.EQ.'COMPO')THEN
        IF(C_ASSOCIATED(IPMAP))CALL XABORT('@CREXSI: ONLY USE '
     1   //'OF EITHER COMPO OR TABLE OPTION. BOTH OPTIONS ARE '
     2   //'NOT ALLOWED.')
      ELSEIF(TEXT.EQ.'TABLE')THEN
        IF(.NOT.C_ASSOCIATED(IPMAP))CALL XABORT('@CREXSI: MISS'
     1   //'ING FUEL MAP.')
      ELSE
        CALL XABORT('@CREXSI: KEYWORD COMPO OR TABLE EXPECTED.')
      ENDIF
*     COMPO NAME
      CALL REDGET(ITYP,NITMA,FLOT,HCOMPO,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@CREXSI: COMPO NAME EXPECTED.')
      DO JEN=2,MAXEN
        IF(HCOMPO.EQ.HENTRY(JEN))THEN
          IPCPO=KENTRY(JEN)
          IF(IMPX.GT.1)CALL LCMLIB(IPCPO)
          GOTO 10
        ENDIF
      ENDDO
      WRITE(HSMG,'(44HCREXSI: UNABLE TO FIND THE COMPO WITH NAME '',
     1      A12,2H''.)') TEXT
      CALL XABORT(HSMG)
*----
*  READ MIX INFO
*----
   10 CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(TEXT.NE.'MIX')CALL XABORT('@CREXSI: KEYWORD MIX EXPECTED.')
      CALL LCMGET(IPCPO,'STATE-VECTOR',IPAR)
      NGRP1=IPAR(2)
      NL1=IPAR(4)
      NISO=IPAR(3)
      IF(NGRP1.NE.NGRP)THEN
        WRITE(HSMG,'(43HCREXSI: INCONSISTENT NB OF GROUPS. IN MACRO,
     1        5HLIB =,I5,11H IN COMPO =,I5)') NGRP,NGRP1
        CALL XABORT(HSMG)
      ENDIF
      IF(NL1.LT.NL)THEN
        WRITE(HSMG,'(43HCREXSI: INCONSISTENT NB OF LEGENDRE ORDERS.,
     1        14H IN MACROLIB =,I5,11H IN COMPO =,I5)') NL,NL1
        CALL XABORT(HSMG)
      ENDIF
   20 ALLOCATE(HISO(3*NISO),ITY(NISO),CONC(NISO))
      CALL CREXSR(IPCPO,LTAB,HCOMPO,NMIX,IMPX,NISO,IBM,DERIV,UPS,
     1            NAMDIR,NISO1,HISO,ITY,CONC,NBURN,KBURN,IVARTY,
     2            IBTYP,BURN0,BURN1)
      JPCPO=LCMGID(IPCPO,NAMDIR)
*----
*  TABLE-OPTION INTERPOLATION
*----
      IF(LTAB)THEN
*       CHECK FUEL MIXTURE
        JPMAP=LCMGID(IPMAP,'FUEL')
        DO 30 IFUEL=1,NFUEL
        KPMAP=LCMGIL(JPMAP,IFUEL)
        CALL LCMGET(KPMAP,'MIX',IMIX)
        IF(IMIX.EQ.IBM)GOTO 40
        CALL LCMLEN(KPMAP,'MIX-VOID',LENGT,ITYP)
        IF(LENGT.EQ.0)GOTO 30
        CALL LCMGET(KPMAP,'MIX-VOID',IMIX)
        IF(IMIX.EQ.IBM)GOTO 40
   30   CONTINUE
        WRITE(IOUT,*)'@CREXSI: UNABLE TO FIND FUEL MIXTURE ',IBM
        CALL XABORT('@CREXSI: WRONG MIXTURE NUMBER.')
*
   40   ALLOCATE(BURNU(NBURN),BRN0(NCH*NB),BRN1(NCH*NB),FMIX(NCH*NB))
        CALL CRERGR(JPCPO,IPMAP,NISO1,NGRP,NMIX,NL,IBM,IMPX,IBTYP,DERIV,
     1  UPS,NBURN,BURNU,ILEAK,TOTAL,ZNUG,SNUGF,CHI,OVERV,DIFFX,DIFFY,
     2  DIFFZ,H,SCAT,IJJ,NJJ,HISO,ITY,CONC,FMIX,BRN0,BRN1,NCH,NB,IVARTY)
        DEALLOCATE(FMIX,BRN1,BRN0,BURNU)
        DEALLOCATE(CONC,ITY,HISO)
*----
*  COMPO-OPTION INTERPOLATION
*----
      ELSE
        ALLOCATE(YTOTAL(NGRP),YZNUG(NGRP),YNUGF(NGRP),YCHI(NGRP),
     1  YOVERV(NGRP),YDIFX(NGRP),YDIFY(NGRP),YDIFZ(NGRP),YH(NGRP),
     2  YSCAT(NL*NGRP*NGRP),YFLUX(NGRP))
        CALL CREINT(JPCPO,NISO1,DERIV,NBURN,KBURN,BURN0,BURN1,NGRP,
     1  NL,IMPX,HISO,ITY,CONC,ILEAK,YTOTAL,YZNUG,YNUGF,YCHI,YOVERV,
     2  YDIFX,YDIFY,YDIFZ,YH,YSCAT,YFLUX,UPS)
*       DATA STORAGE.
        DO 112 JGR=1,NGRP
          TOTAL(IBM,JGR)=YTOTAL(JGR)
          ZNUG(IBM,JGR)=YZNUG(JGR)
          SNUGF(IBM,JGR)=YNUGF(JGR)
          CHI(IBM,JGR)=YCHI(JGR)
          OVERV(IBM,JGR)=YOVERV(JGR)
          DIFFX(IBM,JGR)=YDIFX(JGR)
          DIFFY(IBM,JGR)=YDIFY(JGR)
          DIFFZ(IBM,JGR)=YDIFZ(JGR)
          H(IBM,JGR)=YH(JGR)
          DO 111 IGR=1,NGRP
          DO 110 IL=1,NL
            SCAT(IBM,IL,IGR,JGR)=YSCAT(NL*((JGR-1)*NGRP+IGR-1)+IL)
  110     CONTINUE
  111     CONTINUE
  112   CONTINUE
        DEALLOCATE(YFLUX,YSCAT,YH,YDIFZ,YDIFY,YDIFX,YOVERV,YCHI,YNUGF,
     1  YZNUG,YTOTAL)
        DEALLOCATE(CONC,ITY,HISO)
*       JGR IS THE SECONDARY GROUP.
        DO 135 JGR=1,NGRP
        DO 130 IL=1,NL
        IGMIN=JGR
        IGMAX=JGR
        DO IGR=NGRP,1,-1
          IF(SCAT(IBM,IL,IGR,JGR).NE.0.)THEN
            IGMIN=MIN(IGMIN,IGR)
            IGMAX=MAX(IGMAX,IGR)
          ENDIF
        ENDDO
        IJJ(IBM,IL,JGR)=IGMAX
        NJJ(IBM,IL,JGR)=IGMAX-IGMIN+1
  130   CONTINUE
  135   CONTINUE
      ENDIF
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(TEXT.EQ.'MIX')GOTO 20
  200 CONTINUE
      IF(TEXT.NE.';') CALL XABORT('@CREXSI: FINAL ; EXPECTED.')
      RETURN
      END
