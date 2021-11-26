*DECK HEAT
      SUBROUTINE HEAT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for energy and charge deposition calculation.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1) : create type(L_MACROLIB);
*         HENTRY(2) : read-only ascii file containing HEAT-M data.
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPDEP,IPMAC,JPMAC,KPMAC
      PARAMETER(NSTATE=40)
      CHARACTER TEXT12*12,HSIGN*12,HSMG*131
      INTEGER ISTATE(NSTATE)
      DOUBLE PRECISION DFLOTT,ZNORM,ESUM,CSUM
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL,SGD,FLIN
      REAL, DIMENSION(:,:), ALLOCATABLE :: ZUFIS
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.LE.1) CALL XABORT('HEAT: >1 PARAMETERS EXPECTED.')
      IPDEP=KENTRY(1)
      IPMAC=KENTRY(2)
      IF((IENTRY(1).LE.2).AND.(JENTRY(1).EQ.0)) THEN
         HSIGN='L_DEPOSITION'
         CALL LCMPTC(IPDEP,'SIGNATURE',12,1,HSIGN)
      ELSE IF(IENTRY(1).LE.2) THEN
         CALL LCMGTC(IPDEP,'SIGNATURE',12,1,HSIGN)
         IF(HSIGN.NE.'L_DEPOSITION') THEN
            TEXT12=HENTRY(1)
            CALL XABORT('HEAT: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_DEPOSITION EXPECTED.')
         ENDIF
      ELSE
         CALL XABORT('HEAT: L_DEPOSITION LCM OBJECT EXPECTED.')
      ENDIF
      NGRP=0
      NBMIX=0
      NBFIS=0
      IPICK=0
      DO I=2,NENTRY
        IF((IENTRY(I).NE.1).AND.(IENTRY(I).NE.2)) CALL XABORT('HEAT: L'
     1  //'CM OBJECT EXPECTED AT LHS.')
        IF(JENTRY(I).NE.2) CALL XABORT('HEAT: ENTRY IN READ-ONLY MODE '
     1  //'EXPECTED.')
        HSIGN='L_MACROLIB'
        CALL LCMGTC(KENTRY(I),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_MACROLIB') THEN
          TEXT12=HENTRY(I)
          CALL XABORT('HEAT: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     >    '. L_MACROLIB EXPECTED.')
        ENDIF
        CALL LCMGET(KENTRY(I),'STATE-VECTOR',ISTATE)
        IF(NGRP.EQ.0) THEN
          NGRP=ISTATE(1)
          NBMIX=ISTATE(2)
          NBFIS=ISTATE(4)
        ELSE
          IF(ISTATE(1).NE.NGRP) THEN
            WRITE(HSMG,'(39HHEAT: INVALID NUMBER OF ENERGY GROUPS (,I5,
     1      3H VS,I5,20H) IN MACROLIB NUMBER,I3,1H.)') ISTATE(1),NGRP,I
            CALL XABORT(HSMG)
          ELSE IF(ISTATE(2).NE.NBMIX) THEN
            WRITE(HSMG,'(34HHEAT: INVALID NUMBER OF MIXTURES (,I5,
     1      3H VS,I5,20H) IN MACROLIB NUMBER,I3,1H.)') ISTATE(2),NBMIX,I
            CALL XABORT(HSMG)
          ELSE IF(ISTATE(4).NE.NBFIS) THEN
            WRITE(HSMG,'(42HHEAT: INVALID NUMBER OF FISSILE ISOTOPES (,
     1      I5,3H VS,I5,20H) IN MACROLIB NUMBER,I3,1H.)') ISTATE(4),
     2      NBMIX,I
            CALL XABORT(HSMG)
          ENDIF
        ENDIF
      ENDDO
*----
*  READ INPUT DATA
*----
      IMPX=1
      ZNORM=1.0D0
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF(INDIC.EQ.10) GO TO 60
      IF(INDIC.NE.3) CALL XABORT('HEAT: CHARACTER DATA EXPECTED.')
*
      IF(TEXT12.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('HEAT: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT12.EQ.'POWR') THEN
*        NORMALIZATION TO A GIVEN FISSION POWER.
         CALL REDGET (INDIC,NITMA,POWER,TEXT12,DFLOTT)
         IF(INDIC.NE.2) CALL XABORT('HEAT: REAL DATA EXPECTED.')
*        NORMALIZATION FACTOR FOR THE DIRECT FLUX.
         ALLOCATE(SGD(NBMIX),FLIN(NBMIX))
         ZNORM=0.0D0
         JPMAC=LCMGID(IPMAC,'GROUP')
         DO 30 IGR=1,NGRP
         KPMAC=LCMGIL(JPMAC,IGR)
         CALL LCMGET(KPMAC,'FLUX-INTG',FLIN)
         CALL LCMLEN(KPMAC,'H-FACTOR',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            CALL LCMGET(KPMAC,'H-FACTOR',SGD)
         ELSE
            WRITE(6,'(/43H HEAT: *** WARNING *** NO H-FACTOR FOUND ,
     1      28HON LCM. USE NU*SIGF INSTEAD.)')
            ALLOCATE(ZUFIS(NBMIX,NBFIS))
            CALL XDRSET(SGD,NBMIX,0.0)
            CALL LCMGET(KPMAC,'NUSIGF',ZUFIS)
            DO IBM=1,NBMIX
              DO IFISS=1,NBFIS
                SGD(IBM)=SGD(IBM)+ZUFIS(IBM,IFISS)
              ENDDO
            ENDDO
            DEALLOCATE(ZUFIS)
         ENDIF
         DO 20 IBM=1,NBFIS
         ZNORM=ZNORM+FLIN(IBM)*SGD(IBM)
   20    CONTINUE
   30    CONTINUE
         ZNORM=POWER/ZNORM
         WRITE(6,300) ' DIRECT',ZNORM
         DEALLOCATE(FLIN,SGD)
      ELSE IF(TEXT12.EQ.'SOUR') THEN
*        NORMALIZATION TO A GIVEN SOURCE INTENSITY.
         CALL REDGET (INDIC,NITMA,SNUMB,TEXT12,DFLOTT)
         IF(INDIC.NE.2) CALL XABORT('HEAT: REAL DATA EXPECTED.')
*        NORMALIZATION FACTOR FOR THE DIRECT FLUX.
         ALLOCATE(VOL(NBMIX),SGD(NBMIX))
         CALL LCMGET(IPMAC,'VOLUME',VOL)
         ZNORM=0.0D0
         JPMAC=LCMGID(IPMAC,'GROUP')
         DO 50 IGR=1,NGRP
         KPMAC=LCMGIL(JPMAC,IGR)
         CALL LCMLEN(KPMAC,'FIXE',LENGT,ITYLCM)
         IF(LENGT.EQ.0) THEN
            CALL LCMLIB(KPMAC)
            CALL XABORT('HEAT: SOURCE RECORD MISSING IN MACROLIB.')
         ENDIF
         CALL LCMGET(KPMAC,'FIXE',SGD)
         DO 40 IBM=1,NBMIX
         ZNORM=ZNORM+VOL(IBM)*SGD(IBM)
   40    CONTINUE
   50    CONTINUE
         ZNORM=SNUMB/ZNORM
         WRITE(6,310) ' DIRECT',ZNORM
         DEALLOCATE(SGD,VOL)
      ELSE IF(TEXT12.EQ.';') THEN
         IPICK=0
         GO TO 60
      ELSE IF(TEXT12.EQ.'PICKE') THEN
         IPICK=1
         GO TO 60
      ELSE IF(TEXT12.EQ.'PICKC') THEN
         IPICK=2
         GO TO 60
      ELSE
         CALL XABORT('HEAT: '//TEXT12//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 10
*----
*  COMPUTE THE ENERGY AND CHARGE DEPOSITION
*----
   60 CALL HEADRV(IPDEP,NENTRY-1,KENTRY(2),NBMIX,NGRP,ZNORM,IMPX,ESUM,
     1 CSUM)
*----
*  RECOVER THE TOTAL ENERGY OR CHARGE DEPOSITION AND SAVE IT IN A
*  CLE-2000 VARIABLE
*----
      IF(IPICK.EQ.1) THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.-2) CALL XABORT('HEAT: OUTPUT REAL EXPECTED(1).')
         INDIC=2
         FLOTT=REAL(ESUM)
         CALL REDPUT(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF((INDIC.NE.3).OR.(TEXT12.NE.';')) THEN
           CALL XABORT('HEAT: ; CHARACTER EXPECTED(1).')
         ENDIF      
      ELSE IF(IPICK.EQ.2) THEN
         IF(JENTRY(1).NE.2) CALL XABORT('HEAT: SECOND ENTRY IN READ-O'
     1  //'NLY MODE EXPECTED.')
         CALL REDGET(INDIC,NITMA,FLOT,TEXT12,DFLOTT)
         IF(INDIC.NE.-2) CALL XABORT('HEAT: OUTPUT REAL EXPECTED(2).')
         INDIC=2
         FLOTT=REAL(CSUM)
         CALL REDPUT(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF((INDIC.NE.3).OR.(TEXT12.NE.';')) THEN
           CALL XABORT('HEAT: ; CHARACTER EXPECTED(2).')
         ENDIF      
      ENDIF      
      RETURN
*
  300 FORMAT(/7H HEAT: ,A7,28H FLUX NORMALIZATION FACTOR =,1P,E13.5)
  310 FORMAT(/7H HEAT: ,A7,30H SOURCE NORMALIZATION FACTOR =,1P,E13.5)
      END
