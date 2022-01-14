*DECK LIBCMB
      SUBROUTINE LIBCMB(MAXMIX,MAXISO,NBISO,NEWISO,NNMIX,MIXCMB,VOLTOT,
     >                  VOLFRA,DENMIX,ISONAM,ISONRF,ISHINA,ISOMIX,IHLIB,
     >                  ILLIB,DENISO,TMPISO,LSHI,SNISO,SBISO,NTFG,NIR,
     >                  GIR,MASKI,IEVOL,ITYP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Combine mixtures by volume fraction.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input/output
* MAXMIX  maximum value of nbmix.
* MAXISO  maximum number of isotopes permitted.
* NBISO   number of isotopes before combination.
* NEWISO  number of isotopes after combination.
* NNMIX   new mixture to create or modify.
* MIXCMB  mixture to add.
* VOLTOT  total volume fraction to date.
* VOLFRA  volume fraction of current mixture.
* DENMIX  density of each mixture.
* ISONAM  name of isotopes.
* ISONRF  reference name of isotopes.
* ISHINA  self-shielding name of isotopes.
* ISOMIX  mix number of each isotope.
* IHLIB   isotope options.
* ILLIB   xs library index for each isotope.
* DENISO  density of isotopes.
* TMPISO  temperature of isotopes.
* LSHI    self-shielding flag.
* SNISO   dilution cross section.
* SBISO   dilution cross section used in Livolant-Jeanpierre
*         normalization.
* NTFG    number of thermal inelastic groups,
* NIR     Goldstein-Cohen flag:
*         use IR approximation for groups with index.ge.NIR;
*         use library value if NIR=0.
* GIR     Goldstein-Cohen IR parameter of each isotope.
* MASKI   treat isotope logical.
* IEVOL   depletion suppression flag (=1/2 to suppress/force depletion).
* ITYP    type of isotope.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      MAXMIX,MAXISO,NBISO,NEWISO,NNMIX,MIXCMB,
     >             ISONAM(3,MAXISO),ISONRF(3,MAXISO),ISHINA(3,MAXISO),
     >             ISOMIX(MAXISO),IHLIB(2,MAXISO,4),ILLIB(MAXISO),
     >             LSHI(MAXISO),NTFG(MAXISO),NIR(MAXISO),IEVOL(MAXISO),
     >             ITYP(MAXISO)
      LOGICAL      MASKI(MAXISO)
      REAL         VOLTOT,VOLFRA,DENMIX(MAXMIX),DENISO(MAXISO),
     >             TMPISO(MAXISO),SNISO(MAXISO),SBISO(MAXISO),
     >             GIR(MAXISO)
      DOUBLE PRECISION TOTWPC
*----
*  LOCAL PARAMETERS
*----
      CHARACTER    HSMG*131
*
      CMBVOL=VOLTOT+VOLFRA
      IF(MIXCMB.EQ.NNMIX) GO TO 150
      RMAS1=1.0
      RMAS2=1.0
      IF(MIXCMB.EQ.0) THEN
*----
*  MIXTURE TO ADD IS VOID
*----
        IF(DENMIX(NNMIX).EQ.-1.0) THEN
*----
*  REDUCE ATOMIC DENSITY
*----
          RMAS1=VOLTOT/CMBVOL
        ELSE
*----
*  REDUCE MIXTURE DENSITY BUT NOT WEIGHT PERCENT
*----
          DENMIX(NNMIX)=DENMIX(NNMIX)*VOLTOT/CMBVOL
        ENDIF
      ELSE
*----
*  MIXTURE TO ADD IS NOT VOID
*----
        IF(DENMIX(NNMIX).EQ.-1.0) THEN
          IF(DENMIX(MIXCMB).EQ.-1.0) THEN
*----
*  REDUCE ATOMIC DENSITY
*----
            RMAS1=VOLTOT/CMBVOL
            RMAS2=VOLFRA/CMBVOL
          ELSE
            IF(VOLTOT.GT.0.0)
     >        CALL XABORT('LIBCMB: CANNOT COMBINE MIXTURE WITH '//
     >                    ' WEIGHT PERCENT AND ATOM CONTENTS')
*----
*  TRANSFER MIXTURE DENSITY WITH INITIAL WEIGHT PERCENT TO NEWISO
*----
            DENMIX(NNMIX)=DENMIX(MIXCMB)
          ENDIF
        ELSE
          IF(DENMIX(MIXCMB).EQ.-1.0)
     >      CALL XABORT('LIBCMB: CANNOT COMBINE MIXTURE WITH '//
     >                  ' WEIGHT PERCENT AND ATOM CONTENTS')
*----
*  REDUCE MIXTURE DENSITY AND WEIGHT PERCENT FOR OLD ISO
*  TRANSFER MIXTURE DENSITY WITH REDUCED WEIGHT PERCENT TO NEWISO
*----
          RMAS1=VOLTOT*DENMIX(NNMIX)
          RMAS2=VOLFRA*DENMIX(MIXCMB)
          CMBMAS=RMAS1+RMAS2
          RMAS1=RMAS1/CMBMAS
          RMAS2=RMAS2/CMBMAS
          DENMIX(NNMIX)=CMBMAS/CMBVOL
        ENDIF
      ENDIF
      NEWISO=NBISO
*----
*  RESET OLD DENSITIES
*----
      IF(VOLTOT.EQ.0.0) THEN
        DO 90 ISO=1,NBISO
          IF(ISOMIX(ISO).EQ.NNMIX) THEN
            IF(MASKI(ISO)) THEN
              WRITE(HSMG,'(15HLIBCMB: MIXTURE,I6,18H IS ALREADY DEFINE,
     >        14HD FOR ISOTOPE ,3A4,1H.)') NNMIX,(ISONAM(I,ISO),I=1,3)
              CALL XABORT(HSMG)
            ENDIF
            ISOMIX(ISO)=0
          ENDIF
  90    CONTINUE
      ENDIF
      IF(DENMIX(MIXCMB).EQ.-1.0) THEN
        TOTWPC=1.0D0
      ELSE
        TOTWPC=0.0D0
        DO ISO=1,NBISO
          IF(ISOMIX(ISO).EQ.MIXCMB) THEN
            TOTWPC=TOTWPC+DBLE(DENISO(ISO))
          ENDIF
        ENDDO
        TOTWPC=1.0D0/TOTWPC
      ENDIF
      DO 100 ISO=1,NBISO
        IF(ISOMIX(ISO).EQ.NNMIX) THEN
          DENISO(ISO)=DENISO(ISO)*RMAS1
        ENDIF
 100  CONTINUE
      DO 110 ISO=1,NBISO
        IF(ISOMIX(ISO).EQ.MIXCMB) THEN
*----
*  SCAN ISO IN NNMIX TO IDENTIFY IDENTICAL ISOTOPES
*----
          DO 111 JSO=1,NBISO
            IF(ISOMIX(JSO).EQ.NNMIX) THEN
              IF(ISONRF(1,JSO).EQ.ISONRF(1,ISO).AND.
     >           ISONRF(2,JSO).EQ.ISONRF(2,ISO)) THEN
                IF(ISONAM(1,JSO).NE.ISONAM(1,ISO).OR.
     >             ISONAM(2,JSO).NE.ISONAM(2,ISO).OR.
     >             TMPISO(JSO)  .NE.TMPISO(ISO)  .OR.
     >             LSHI(JSO)    .NE.LSHI(ISO)    .OR.
     >             SNISO(JSO)   .NE.SNISO(ISO)   .OR.
     >             SBISO(JSO)   .NE.SBISO(ISO)        ) THEN
                  WRITE(HSMG,'(17HLIBCMB: ISOTOPES ,3A4,5H AND ,3A4,
     >            18H CANNOT BE MERGED.)') (ISONAM(I,ISO),I=1,3),
     >            (ISONAM(I,JSO),I=1,3)
                  CALL XABORT(HSMG)
                ENDIF
                DENISO(JSO)=DENISO(JSO)+REAL(TOTWPC)*DENISO(ISO)*RMAS2
                GO TO 115
              ENDIF
            ENDIF
 111      CONTINUE
          ISO2=0
          DO 112 JSO=1,NBISO
            IF(ISOMIX(JSO).EQ.0) THEN
              ISO2=JSO
              GO TO 113
            ENDIF
 112      CONTINUE
          NEWISO=NEWISO+1
          IF(NEWISO.GT.MAXISO) CALL XABORT('LIBCMB: MAXISO OVERFLOW.')
          ISO2=NEWISO
 113      ISOMIX(ISO2)=NNMIX
          DENISO(ISO2)=REAL(TOTWPC)*DENISO(ISO)*RMAS2
          TMPISO(ISO2)=TMPISO(ISO)
          NTFG(ISO2)=NTFG(ISO)
          NIR(ISO2)=NIR(ISO)
          GIR(ISO2)=GIR(ISO)
          SNISO(ISO2)=SNISO(ISO)
          SBISO(ISO2)=SBISO(ISO)
          LSHI(ISO2)=LSHI(ISO)
          MASKI(ISO2)=.TRUE.
          IEVOL(ISO2)=IEVOL(ISO)
          ITYP(ISO2)=ITYP(ISO)
          DO 120 ITC=1,3
            ISONAM(ITC,ISO2)=ISONAM(ITC,ISO)
            ISONRF(ITC,ISO2)=ISONRF(ITC,ISO)
            ISHINA(ITC,ISO2)=ISHINA(ITC,ISO)
 120      CONTINUE
          DO 140 ILC=1,4
            DO 141 ITC=1,2
              IHLIB(ITC,ISO2,ILC)=IHLIB(ITC,ISO,ILC)
 141        CONTINUE
 140      CONTINUE
          ILLIB(ISO2)=ILLIB(ISO)
        ENDIF
 115    CONTINUE
 110  CONTINUE
 150  NBISO=NEWISO
      VOLTOT=CMBVOL
      RETURN
      END
