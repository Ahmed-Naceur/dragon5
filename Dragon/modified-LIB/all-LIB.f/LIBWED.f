*DECK LIBWED
      SUBROUTINE LIBWED(MAXR,NEL,NDEPL,NDFI,NDFP,NHEAVY,NLIGHT,NOTHER,
     >                  NREAC,NPAR,ITNAM,ITZEA,MATNO,KPAX,BPAX,HNADPL,
     >                  IZEA,IDR,RER,RRD,KPAR,BPAR,YIELD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Create /depletion/ records in /microlib/.
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
* MAXR    number of reaction types.
* NEL     number of isotopes on library.
* NDEPL   number of depleting isotopes.
* NDFI    number of direct fissile isotopes.
* NDFP    number of direct fission product.
* NHEAVY  number of heavy isotopes (fissile isotopes + decay +
*         capture isotopes).
* NLIGHT  number of light isotopes (fission product + decay +
*         capture isotopes).
* NOTHER  number of other isotopes (depleting isotopes + decay +
*         capture isotopes).
* NREAC   maximum number of depletion reaction in the depletion chain.
* NPAR    maximum number of parent isopopes from decay and capture.
* ITNAM   reactive isotope names in chain.
* ITZEA   6-digit nuclide identifier:
*         atomic number z*10000 (digits) + mass number a*10 +
*         energy state (0 = ground state, 1 = first state, etc.).
* MATNO   reaction material index.
* KPAX    complete reaction type matrix.
* BPAX    complete branching ratio matrix.
*
*Parameters: output
* HNADPL  reactive isotope names in chain.
* IZEA    6-digit nuclide identifier:
*         atomic number z*10000 (digits) + mass number a*10 +
*         energy state (0 = ground state, 1 = first state, etc.).
* IDR     DEPLETE-REAC matrix (reaction identifiers).
* RER     DEPLETE-ENER matrix (MeV/reaction values).
* RRD     DEPLETE-DECA vector (decay constant values).
* KPAR    PRODUCE-REAC matrix (production identifiers).
* BPAR    PRODUCE-RATE matrix (branching ratios).
* YIELD   fission product yield matrix.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXR,NEL,NDEPL,NDFI,NDFP,NHEAVY,NLIGHT,NOTHER,NREAC,
     > NPAR,ITNAM(3,NEL),ITZEA(NEL),MATNO(NEL),KPAX(NEL+MAXR,NEL),
     > HNADPL(3,NDEPL),IZEA(NDEPL),IDR(NREAC,NDEPL),KPAR(NPAR,NDEPL)
      REAL BPAX(NEL+MAXR,NEL),RER(NREAC,NDEPL),RRD(NDEPL),
     > BPAR(NPAR,NDEPL),YIELD(NDFI,NDFP)
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXT4*4
      INTEGER ITEXT4,ISOHEA,ISOLIG,ISOOTH,ISOSTA,ISO,ISONUM,IT,IEL,
     > IPAR,JEL,JSO,IFI,IFP
*----
*  INTERNAL PARAMETERS
*----
      INTEGER KDECAY,KFISSP,KFISSI,KHEAT
      PARAMETER (KDECAY=1,KFISSP=2,KFISSI=6,KHEAT=9)
*----
*  INITIALIZE DECAY CHAIN
*----
      TEXT4='    '
      READ(TEXT4,'(A4)') ITEXT4
      CALL XDISET(HNADPL,3*NDEPL,ITEXT4)
      CALL XDISET(IZEA,NDEPL,0)
      CALL XDISET(IDR,NREAC*NDEPL,0)
      CALL XDRSET(RER,NREAC*NDEPL,0.0)
      CALL XDRSET(RRD,NDEPL,0.0)
      CALL XDISET(KPAR,NPAR*NDEPL,0)
      CALL XDRSET(BPAR,NPAR*NDEPL,0.0)
      CALL XDRSET(YIELD,NDFI*NDFP,0.0)
*----
*  RENUMBER ISOTOPES AND SAVE IDR, RER AND RRD
*----
      ISOHEA=0
      ISOLIG=ISOHEA+NHEAVY
      ISOOTH=ISOLIG+NLIGHT
      ISOSTA=ISOOTH+NOTHER
      DO 100 ISO=NEL,1,-1
        ISONUM=0
        IF(MATNO(ISO).EQ.-KFISSI) THEN
          ISOHEA=ISOHEA+1
          ISONUM=ISOHEA
        ELSE IF(MATNO(ISO).EQ.-KFISSP) THEN
          ISOLIG=ISOLIG+1
          ISONUM=ISOLIG
        ELSE IF(MATNO(ISO).EQ.-KDECAY) THEN
          ISOOTH=ISOOTH+1
          ISONUM=ISOOTH
        ELSE IF(MATNO(ISO).EQ.-KHEAT) THEN
          ISOSTA=ISOSTA+1
          ISONUM=ISOSTA
        ENDIF
        IF(ISONUM.GT.0) THEN
          MATNO(ISO)=ISONUM
          HNADPL(1,ISONUM)=ITNAM(1,ISO)
          HNADPL(2,ISONUM)=ITNAM(2,ISO)
          HNADPL(3,ISONUM)=ITNAM(3,ISO)
          IZEA(ISONUM)=ITZEA(ISO)
          IDR(1,ISONUM)=KPAX(NEL+1,ISO)
          RRD(ISONUM)=BPAX(NEL+1,ISO)
          DO 101 IT=2,NREAC
            IDR(IT,ISONUM)=KPAX(NEL+IT,ISO)
            RER(IT,ISONUM)=BPAX(NEL+IT,ISO)
 101      CONTINUE
        ENDIF
 100  CONTINUE
*----
*  CREATE KPAR AND BPAR MATRIX
*----
      DO 110 IEL=1,NEL
        ISO=MATNO(IEL)
        IF(ISO.GT.0) THEN
          IPAR=0
          DO 111 JEL=1,NEL
            JSO=MATNO(JEL)
            IF(JSO.GT.0) THEN
              IF((KPAX(IEL,JEL).NE.0).AND.(KPAX(IEL,JEL).NE.KFISSP))
     >        THEN
                IPAR=IPAR+1
                IF(IPAR.GT.NPAR)
     >            CALL XABORT('LIBWED: TOO MANY DECAY PARENTS')
                KPAR(IPAR,ISO)=JSO*100+KPAX(IEL,JEL)
                BPAR(IPAR,ISO)=BPAX(IEL,JEL)
              ENDIF
            ENDIF
 111      CONTINUE
        ENDIF
 110  CONTINUE
*----
*  CREATE YIELD MATRIX
*----
      DO 120 IEL=1,NEL
        ISO=MATNO(IEL)
        IF(ISO.GT.0) THEN
          IF(MOD(IDR(KFISSP,ISO),100).EQ.4) THEN
            IFI=IDR(KFISSP,ISO)/100
            IF(IFI.EQ.0) GO TO 120
            IF(IFI.GT.NDFI)
     >        CALL XABORT('LIBWED: INVALID FISSILE ISOTOPE NUMBER')
            DO 121 JEL=1,NEL
              JSO=MATNO(JEL)
              IF(JSO.GT.0) THEN
                IF(MOD(IDR(KFISSP,JSO),100).EQ.2) THEN
                  GO TO 121
                ELSE IF(MOD(IDR(KFISSP,JSO),100).EQ.5) THEN
                  IFP=IDR(KFISSP,JSO)/100
                  IF(IFP.GT.NDFP) CALL XABORT('LIBWED: INVALID FI'//
     >               'SSION PRODUCT NUMBER')
                  IF(KPAX(JEL,IEL) .EQ. KFISSP) THEN
                    YIELD(IFI,IFP)=BPAX(JEL,IEL)
                  ENDIF
                ENDIF
              ENDIF
 121        CONTINUE
          ENDIF
        ENDIF
 120  CONTINUE
*----
*  RETURN
*----
      RETURN
      END
