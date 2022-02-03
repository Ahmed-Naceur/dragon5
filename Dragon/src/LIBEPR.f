*DECK LIBEPR
      SUBROUTINE LIBEPR(IMPX,NDEPL,NSTABL,NDFI,NDFP,NREAC,NPAR,HNADPL,
     >                  NMDEPL,IDR,RER,RRD,KPAR,BPAR,YIELD,IZAE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print isotopic depletion chain.
* 
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau
*
*Parameters: input
* IMPX    print parameter:
*         .ge.5 nice formatted printout;
*         .ge.10 dragon input formatted;
*         .ge.20 dragon/APOLIB-2 input formatted.
* NDEPL   number of depleting isotopes.
* NSTABL  number of non-depleting isotopes producing energy.
* NDFI    number of direct fissile isotopes.
* NDFP    number of direct fission product.
* NREAC   maximum number of depletion reaction in the depletion chain.
* NPAR    maximum number of parent isopopes from decay and 
*         neutron-induced reactions.
* HNADPL  reactive isotope names in chain.
* NMDEPL  names of used depletion reactions
*         (NMDEPL(1)='DECAY'; NMDEPL(2)='NFTOT';
*          NMDEPL(3)='NG'   ; NMDEPL(4)='N2N';
*          etc.).
* IDR     DEPLETE-REAC matrix (reaction identifiers).
* RER     DEPLETE-ENER matrix (MeV/reaction values).
* RRD     DEPLETE-DECA vector (decay constant values).
* KPAR    PRODUCE-REAC matrix (production identifiers).
* BPAR    PRODUCE-RATE matrix (branching ratios).
* YIELD   fission product yield matrix.
* IZAE    6-digit nuclide identifier:
*         atomic number z*10000 (digits) + mass number a*10 +
*         energy state (0 = ground state, 1 = first state, etc.).
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER    IMPX,NDEPL,NSTABL,NDFI,NDFP,NREAC,NPAR,HNADPL(3,NDEPL),
     >           NMDEPL(2,NREAC),IDR(NREAC,NDEPL),KPAR(NPAR,NDEPL),
     >           IZAE(NDEPL)
      REAL       RER(NREAC,NDEPL),RRD(NDEPL),BPAR(NPAR,NDEPL),
     >           YIELD(NDFI,NDFP)
*----
*  LOCAL VARIABLES
*----
      INTEGER    I,J,ISO,KRD,KRF,ISOF,IPF,KRC,IOF,KROTH,NBPROD,IPAR,
     >           JSO,KREAC,ISOFP,IOFS,IREAC,NBFIS,IOUT
      REAL       RATD1,RATD2,RATF,RATFS,RATC,RREAC
*----
*  INTERNAL PARAMETERS
*----
      REAL        CNVDAY
      PARAMETER  (IOUT=6,CNVDAY=1.0E+8/86400.0)
      CHARACTER   NONOFF(2)*3
      CHARACTER   LINE*130
      SAVE        NONOFF
      DATA        NONOFF/'NO ','YES'/
*----
*  PRINT THE DEPLETION CHAIN.
*----
      IF(IMPX.GE.5) THEN
        WRITE(IOUT,6000) ((NMDEPL(J,I),J=1,2),I=4,NREAC)
        WRITE(IOUT,6010)
        DO 110 ISO=1,NDEPL
          KRD=1
          RATD1=0.0
          RATD2=0.0
          IF(IDR(1,ISO).NE.0) THEN
            RATD1=CNVDAY/RRD(ISO)
            RATD2=RER(1,ISO)
            KRD=2
          ENDIF
          KRF=1
          RATF=0.0
          RATFS=0.0
          IF(MOD(IDR(2,ISO),100).EQ.3) THEN
            KRF=2
            RATF=RER(2,ISO)
          ELSE IF(MOD(IDR(2,ISO),100).EQ.4) THEN
            KRF=2
            RATF=RER(2,ISO)
            ISOF=IDR(2,ISO)/100
            IF(ISOF.NE.0) THEN
              DO 111 IPF=1,NDFP
                RATFS=RATFS+YIELD(ISOF,IPF)
 111          CONTINUE
            ENDIF
          ENDIF
          KRC=1
          RATC=0.0
          IF(IDR(3,ISO).NE.0) THEN
            RATC=RER(3,ISO)
            KRC=2
          ENDIF
          LINE=' '
*----
*  WRITE ISOTOPE PROPERTIES
*----
          WRITE(LINE(:18),'(I5,1X,2A4,1X,A3)') ISO,HNADPL(1,ISO),
     >    HNADPL(2,ISO),NONOFF(KRD)
          IF(KRD.EQ.2) WRITE(LINE(19:42),'(1P,2E12.4)') RATD1,RATD2
          WRITE(LINE(45:47),'(A3)') NONOFF(KRF)
          IF(KRF.EQ.2) WRITE(LINE(48:71),'(1P,2E12.4)') RATF,RATFS
          WRITE(LINE(74:76),'(A3)') NONOFF(KRC)
          IF(KRC.EQ.2) WRITE(LINE(77:88),'(1P,E12.4)') RATC
          IOF=91
          DO 112 I=4,NREAC
            KROTH=1
            IF(IDR(I,ISO).GT.0) KROTH=2
            IF(IOF+7.GT.130) THEN
              WRITE(IOUT,'(1X,A)') LINE
              IOF=91
              LINE=' '
            ENDIF
            WRITE(LINE(IOF:IOF+7),'(A3,5X)') NONOFF(KROTH)
            IOF=IOF+8
 112      CONTINUE
          WRITE(IOUT,'(1X,A)') LINE
 110    CONTINUE
*----
*  WRITE PARENTS FROM ALL REACTION EXCEPT FISSION
*----
        WRITE(IOUT,7000)
        DO 210 ISO=1,NDEPL-NSTABL
          NBPROD=0
          DO 120 IPAR=1,NPAR
            JSO=KPAR(IPAR,ISO)/100
            KREAC=MOD(KPAR(IPAR,ISO),100)
            RREAC=BPAR(IPAR,ISO)
            IF(JSO.GT.0) THEN
              IF((KREAC.LE.0).OR.(KREAC.GT.NREAC)) CALL XABORT('LIBEPR'
     >        //': INDALID REACTION INDEX')
              NBPROD=NBPROD+1
              IF(NBPROD.EQ.1) THEN
                WRITE(IOUT,6012) ISO,HNADPL(1,ISO),HNADPL(2,ISO),
     >          NMDEPL(1,KREAC),NMDEPL(2,KREAC),JSO,HNADPL(1,JSO),
     >          HNADPL(2,JSO),RREAC
              ELSE
                WRITE(IOUT,6013) NMDEPL(1,KREAC),NMDEPL(2,KREAC),
     >          JSO,HNADPL(1,JSO),HNADPL(2,JSO),RREAC
              ENDIF
            ENDIF
 120      CONTINUE
*----
*  WRITE PARENTS FROM FISSION IF REQUIRED
*----
          IF(MOD(IDR(2,ISO),100).EQ.2) THEN
            GO TO 210
          ELSE IF(MOD(IDR(2,ISO),100).EQ.5) THEN
            ISOFP=IDR(2,ISO)/100
            IF(ISOFP.GT.NDFP)
     >        CALL XABORT('LIBEPR: INVALID FISSION PRODUCT NUMBER')
            DO 130 JSO=1,NDEPL
              IF(MOD(IDR(2,JSO),100).EQ.4) THEN
                ISOF=IDR(2,JSO)/100
                IF(ISOF.GT.NDFI) THEN
                  CALL XABORT('LIBEPR: INVALID FISSILE NUMBER')
                ELSE IF(ISOF.GT.0) THEN
                  RREAC=YIELD(ISOF,ISOFP)
                  IF(RREAC.GT.0.0) THEN
                    NBPROD=NBPROD+1
                    IF(NBPROD.EQ.1) THEN
                      WRITE(IOUT,6012) ISO,HNADPL(1,ISO),HNADPL(2,ISO),
     >                NMDEPL(1,2),NMDEPL(2,2),JSO,HNADPL(1,JSO),
     >                HNADPL(2,JSO),RREAC
                    ELSE
                      WRITE(IOUT,6013) NMDEPL(1,2),NMDEPL(2,2),JSO,
     >                HNADPL(1,JSO),HNADPL(2,JSO),RREAC
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
 130        CONTINUE
          ENDIF
 210    CONTINUE
      ENDIF
*
      IF(IMPX.GE.10) THEN
        WRITE(IOUT,'(/1X,4HDEPL,I6,6H CHAIN)') NDEPL
        DO 330 ISO=1,NDEPL
        LINE=' '
        WRITE(LINE(:17),'(1H'',2A4,1H'',I7)') HNADPL(1,ISO),
     >  HNADPL(2,ISO),IZAE(ISO)
        IOFS=18
        IF(IDR(1,ISO).NE.0) THEN
          WRITE(LINE(IOFS:IOFS+20),'(1X,5HDECAY,1P,E14.6)') RRD(ISO)
          IOFS=IOFS+20
*          IF(RER(1,ISO).NE.0.0) THEN
*            WRITE(LINE(IOFS:IOFS+13),'(1X,1P,E12.5)') RER(1,ISO)
*            IOFS=IOFS+13
*          ENDIF
        ENDIF
        IF((MOD(IDR(2,ISO),100).EQ.3).OR.(MOD(IDR(2,ISO),100).EQ.4))
     >  THEN
          WRITE(LINE(IOFS:IOFS+15),'(1X,5HNFTOT,F9.4)') RER(2,ISO)
          IOFS=IOFS+15
        ENDIF
        IF(IDR(3,ISO).NE.0) THEN
          WRITE(LINE(IOFS:IOFS+17),'(1X,2HNG,1P,E14.6)') RER(3,ISO)
          IOFS=IOFS+17
        ENDIF
        DO 300 IREAC=4,NREAC
        IF(IDR(IREAC,ISO).NE.0) THEN
          IF(IOFS+9.GT.71) THEN
             WRITE(IOUT,'(1X,A)') LINE
             LINE=' '
             IOFS=18
          ENDIF
          WRITE(LINE(IOFS:IOFS+9),'(1X,2A4)') NMDEPL(1,IREAC),
     >    NMDEPL(2,IREAC)
          IOFS=IOFS+9
        ENDIF
 300    CONTINUE
        IF(ISO.GT.NDEPL-NSTABL) THEN
          WRITE(LINE(IOFS:IOFS+7),'(7H STABLE)')
          IOFS=IOFS+7
          IF(IOFS.GT.71) CALL XABORT('LIBEPR: LINE OVERFLOW.')
        ENDIF
        IF(LINE.NE.' ') WRITE(IOUT,'(1X,A)') LINE
        LINE=' '
        IOFS=7
        NBPROD=0
        DO 310 IPAR=1,NPAR
          JSO=KPAR(IPAR,ISO)/100
          KREAC=MOD(KPAR(IPAR,ISO),100)
          RREAC=BPAR(IPAR,ISO)
          IF(JSO.GT.0) THEN
            NBPROD=NBPROD+1
            IF(NBPROD.EQ.1) WRITE(LINE(2:6),'(5HFROM )')
            IF(IOFS+34.GE.71) THEN
               WRITE(IOUT,'(1X,A)') LINE
               LINE=' '
               IOFS=7
            ENDIF
            WRITE(LINE(IOFS:IOFS+34),6100)
     >      NMDEPL(1,KREAC),NMDEPL(2,KREAC),RREAC,HNADPL(1,JSO),
     >      HNADPL(2,JSO)
            IOFS=IOFS+34
          ENDIF
 310    CONTINUE
        IF(LINE.NE.' ') WRITE(IOUT,'(1X,A)') LINE
        LINE=' '
        IOFS=16
        NBFIS=0
        IF(MOD(IDR(2,ISO),100).EQ.2) THEN
          GO TO 330
        ELSE IF(MOD(IDR(2,ISO),100).EQ.5) THEN
          ISOFP=IDR(2,ISO)/100
          DO 320 JSO=1,NDEPL
            IF(MOD(IDR(2,JSO),100).EQ.4) THEN
              ISOF=IDR(2,JSO)/100
              IF(ISOF.GT.0) THEN
                RREAC=YIELD(ISOF,ISOFP)
                IF(RREAC.GT.0.0) THEN
                  NBPROD=NBPROD+1
                  IF(NBPROD.EQ.1) WRITE(LINE(2:6),'(5HFROM )')
                  NBFIS=NBFIS+1
                  IF(NBFIS.EQ.1) WRITE(LINE(7:14),'(A8)') 'NFTOT   '
                  WRITE(LINE(IOFS:IOFS+24),6101)
     >            RREAC,HNADPL(1,JSO),HNADPL(2,JSO)
                  IOFS=IOFS+24
                  IF(IOFS.GE.60) THEN
                     WRITE(IOUT,'(1X,A)') LINE
                     LINE=' '
                     IOFS=16
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
 320      CONTINUE
          IF(LINE.NE.' ') WRITE(IOUT,'(1X,A)') LINE
        ENDIF
 330    CONTINUE
        WRITE(IOUT,'(9H ENDCHAIN)')
      ENDIF
*
      IF(IMPX.GE.500) THEN
        WRITE(IOUT,'(/1X,33HDEPL LIB: APLIB2 FIL: XXXXX CHAIN)')
        DO 350 ISO=1,NDEPL-NSTABL
        LINE=' '
        WRITE(LINE(:8),'(2A4)') HNADPL(1,ISO),HNADPL(2,ISO)
        IOFS=14
        NBPROD=0
        DO 340 IPAR=1,NPAR
          JSO=KPAR(IPAR,ISO)/100
          KREAC=MOD(KPAR(IPAR,ISO),100)
          RREAC=BPAR(IPAR,ISO)
          IF(JSO.GT.0) THEN
            NBPROD=NBPROD+1
            IF(NBPROD.EQ.1) WRITE(LINE(10:13),'(4HFROM)')
            WRITE(LINE(IOFS:IOFS+29),'(1X,2A4,1P,E11.4,1X,2A4)')
     >      NMDEPL(1,KREAC),NMDEPL(2,KREAC),RREAC,HNADPL(1,JSO),
     >      HNADPL(2,JSO)
            IOFS=IOFS+29
            IF(IOFS.GE.71) THEN
               WRITE(IOUT,'(1X,A)') LINE
               LINE=' '
               IOFS=14
            ENDIF
          ENDIF
 340    CONTINUE
        IF(LINE.NE.' ') WRITE(IOUT,'(1X,A)') LINE
 350    CONTINUE
        WRITE(IOUT,'(9H ENDCHAIN)')
      ENDIF
*----
*  RETURN FROM LIBEPR
*----
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(/' DEPLETION CHAIN DATA'/' --------------------'//
     >       43X,'DEPLETION REACTIONS'/
     >       '                -----------------------------',
     >       '--------------------------------------------'/
     >       '   ISOTOPE      ......DECAY................  ..........',
     >       '.FISSION.........  ......NG.......  ',10A4:/(91X,10A4))
 6010 FORMAT('   NB. NAME         MLIFE(DAYS) ENERGY(MEV)      ENERGY',
     >       '(MEV) TOTAL YIELD      ENERGY(MEV)')
 6012 FORMAT(I5,1X,2A4,4X,2A4,3X,I5,1X,2A4,1X,1P,E14.6)
 6013 FORMAT(18X,2A4,3X,I5,1X,2A4,1X,1P,E14.6)
 6100 FORMAT(2A4,1X,1P,E13.6,2H ',2A4,2H' )
 6101 FORMAT(1P,E13.6,2H ',2A4,1H')
 7000 FORMAT(//27X,'PRODUCTION REACTIONS'/
     >       18X,'---------------------------------------'/
     >       '  ISOTOPE         REACTION     ISOTOPE',14X,'YIELD'/
     >       '  NB. NAME                     NB. NAME    ')
      END
