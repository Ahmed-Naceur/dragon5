*DECK MRGGET
      SUBROUTINE MRGGET(IPRINT,NSOUTO,NVOUTO,NSOUTN,NVOUTN,
     >                  IUPD,IMERGE,MIXN,ALBEDN)
*
*----------
*
*Purpose:
* Read merge options.
*
*Copyright:
* Copyright (C) 1997 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPRINT  print level.
* NSOUTO  old number of surfaces.
* NVOUTO  old number of regions.
*
*Parameters: output
* NSOUTN  new number of surfaces.
* NVOUTN  new number of regions.
* IUPD    type of merge required:
*         IUPD(1) for region merge and track splitting;
*         IUPD(2) for surface merge;
*         IUPD(3) for material modification;
*         IUPD(4) for albedo modification.
* IMERGE  merged position.
* MIXN    new material for old regions.
* ALBEDN  new surface albedo.
*
*Comments:
*  Input options: 
*     [ EDIT iprint ]
*     [ REGI (imerge(ii),ii=1,nvouto)]             -> IUPD(1) > 0
*     [ EXTR (imerge(ii),ii=1,-IUPD(1))]           -> IUPD(1) < 0
*     [ SURF (imerge(ii),ii=-1,-nsouto)]           -> IUPD(2) < 0
*     [ { OLDM (mixn(ii),ii=1,nvouto) |            -> IUPD(3) > 0
*         NEWM (mixn(ii),ii=1,nvouto) } ]          -> IUPD(3) < 0
*     [ ALBE (albedn(ii),ii=1,6)]                  -> IUPD(4) > 0
*
*----------
*
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='MRGGET')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          IPRINT,NSOUTO,NVOUTO,NSOUTN,NVOUTN,
     >                 IUPD(4),
     >                 IMERGE(-NSOUTO:NVOUTO),
     >                 MIXN(NVOUTO)
      REAL             ALBEDN(6)
*----
*  LOCAL VARIABLES
*----
      INTEGER          IVSN,IVSO,ITYPLU,INTLIR
      CHARACTER        CARLIR*12
      REAL             REALIR
      DOUBLE PRECISION DBLINP
*----
*  INITIALIZE IMERGE
*----
      CALL XDISET(IUPD,4,0)
      CALL XDISET(MIXN,NVOUTO,0)
      CALL XDRSET(ALBEDN,6,0.0)
      DO IVSO=-NSOUTO,NVOUTO
        IMERGE(IVSO)=IVSO
      ENDDO
      NSOUTN=0
      NVOUTN=0
      IPRINT = 1
*----
*  READ OPTION NAME
*----
 110  CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
 111  IF(ITYPLU.NE.3) CALL XABORT(NAMSBR//
     >  ': READ ERROR - CHARACTER VARIABLE EXPECTED')
      IF(CARLIR(1:1) .EQ. ';') THEN
        GO TO 115
      ELSE IF(CARLIR .EQ. 'EDIT') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
        IF(ITYPLU .NE. 1) GO TO 111
        IPRINT=INTLIR
      ELSE IF(CARLIR(1:4) .EQ. 'REGI') THEN
        IF(IUPD(1) .NE. 0 ) CALL XABORT(NAMSBR//
     >      ': A single REGI or EXTR permitted')
        DO IVSO=1,NVOUTO
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
          IF(ITYPLU .NE. 1) GO TO 111
          IF(INTLIR .GT. NVOUTO) CALL XABORT(NAMSBR//
     >      ': FINAL REGION NUMBER MUST BE SMALLER '//
     >      'THAN NUMBER OF ORIGINAL REGIONS')
          IF(INTLIR .LE. 0) CALL XABORT(NAMSBR//
     >      ': FINAL REGION NUMBER MUST LARGER THAN 0 ')
          IUPD(1)=IUPD(1)+1
          NVOUTN=MAX(NVOUTN,INTLIR)
          IMERGE(IVSO)=INTLIR
        ENDDO
      ELSE IF(CARLIR(1:4) .EQ. 'EXTR') THEN
        IF(IUPD(1) .NE. 0 ) CALL XABORT(NAMSBR//
     >      ': A single REGI or EXTR permitted')
        DO IVSO=1,NVOUTO
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
          IF(ITYPLU .NE. 1) GO TO 111
          IF(INTLIR .GT. NVOUTO) CALL XABORT(NAMSBR//
     >      ': FINAL REGION NUMBER MUST BE SMALLER '//
     >      'THAN NUMBER OF ORIGINAL REGIONS')
          IF(INTLIR .LE. 0) CALL XABORT(NAMSBR//
     >      ': FINAL REGION NUMBER MUST LARGER THAN 0 ')
          IUPD(1)=IUPD(1)-1
          NVOUTN=MAX(NVOUTN,INTLIR)
          IMERGE(IVSO)=INTLIR
        ENDDO
      ELSE IF(CARLIR(1:4).EQ.'SURF') THEN
        DO IVSO=1,NSOUTO
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
          IF(ITYPLU.NE.1) GO TO 111
          IF(INTLIR.GT.NSOUTO) CALL XABORT(NAMSBR//
     >      ': FINAL SURFACE NUMBER MUST BE SMALLER '//
     >      'THAN NUMBER OF ORIGINAL SURFACES')
          IF(INTLIR .LE. 0) CALL XABORT(NAMSBR//
     >      ': FINAL SURFACE NUMBER MUST LARGER THAN 0 ')
          IUPD(2)=IUPD(2)-1
          NSOUTN=MAX(NSOUTN,INTLIR)
          IMERGE(-IVSO)=-INTLIR
        ENDDO
      ELSE IF(CARLIR(1:4) .EQ. 'OLDM') THEN
        DO IVSO=1,NVOUTO
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
          IF(ITYPLU .NE. 1) GO TO 111
          IF(INTLIR .LT. 0) CALL XABORT(NAMSBR//
     >      ': FINAL MIXTURE NUMBER MUST LARGER OR EQUAL TO 0 ')
          IUPD(3)=IUPD(3)+1
          MIXN(IVSO)=INTLIR
        ENDDO
      ELSE IF(CARLIR(1:4) .EQ. 'NEWM') THEN
        DO IVSO=1,NVOUTO
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
          IF(ITYPLU .NE. 1) GO TO 111
          IF(INTLIR .LT. 0) CALL XABORT(NAMSBR//
     >      ': FINAL MIXTURE NUMBER MUST LARGER OR EQUAL TO 0 ')
          IUPD(3)=IUPD(3)-1
          MIXN(IVSO)=INTLIR
        ENDDO
      ELSE IF(CARLIR(1:4) .EQ. 'ALBE') THEN
        DO IVSO=1,6
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLINP)
          IF(ITYPLU .NE. 3) GO TO 111
          IF(REALIR .LT. 0.0) CALL XABORT(NAMSBR//
     >      ': FINAL ALBEDO MUST LARGER THAN 0.0 ')
          IUPD(4)=IUPD(4)+1
          ALBEDN(IVSO)=REALIR
        ENDDO
      ELSE
        CALL XABORT(NAMSBR//': LEGAL KEYWORD '//CARLIR)
      ENDIF
      GO TO 110
 115  CONTINUE
*----
*  CHECK IF ALL THE SUCCESSIVE REGIONS CONSIDERED
*----
      IF(IUPD(1) .EQ. 0) THEN
        NVOUTN=NVOUTO
      ENDIF
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6000) NAMSBR
        IF(IUPD(1) .GT. 0) THEN
          WRITE(IOUT,6010)
          WRITE(IOUT,6020) (IVSO,IMERGE(IVSO),IVSO=1,IUPD(1))
        ELSE IF(IUPD(1) .LT. 0) THEN
          WRITE(IOUT,6015)
          WRITE(IOUT,6025) (IMERGE(IVSO),IVSO=1,-IUPD(1))
        ENDIF
        IF(IUPD(2) .LT. 0) THEN
          WRITE(IOUT,6011)
          WRITE(IOUT,6020) (IVSO,IMERGE(IVSO),IVSO=-1,IUPD(2),-1)
        ENDIF
        IF(IUPD(3) .LT. 0) THEN
          WRITE(IOUT,6012)
          WRITE(IOUT,6020) (IVSO,MIXN(IVSO),IVSO=1,-IUPD(3))
        ELSE IF(IUPD(3) .GT. 0) THEN
          WRITE(IOUT,6013)
          WRITE(IOUT,6020) (IVSO,MIXN(IVSO),IVSO=1,IUPD(3))
        ENDIF
        IF(IUPD(4) .GT. 0) THEN
          WRITE(IOUT,6014)
          WRITE(IOUT,6021) (IVSO,ALBEDN(IVSO),IVSO=1,6)
        ENDIF
        WRITE(IOUT,6001)
      ENDIF
      IF(IUPD(1) .GT. 0) THEN
        DO IVSN=1,NVOUTN
          DO IVSO=1,NVOUTO
            IF(IMERGE(IVSO) .EQ. IVSN) GO TO 205
          ENDDO
          CALL XABORT(NAMSBR//
     >      ': NEW REGION NUMBERS NOT SUCCESSIVE')
 205      CONTINUE
        ENDDO
      ENDIF
*----
*  CHECK IF ALL THE SUCCESSIVE SURFACE CONSIDERED
*----
      IF(IUPD(2) .EQ. 0) THEN
        NSOUTN=NSOUTO
      ENDIF
      IF(IUPD(2) .GT. 0) THEN
        DO IVSN=-NSOUTN,-1
          DO IVSO=-NSOUTO,-1
            IF(IMERGE(IVSO) .EQ. IVSN) GO TO 215
          ENDDO
          CALL XABORT(NAMSBR//
     >      ': NEW SURFACE NUMBERS NOT SUCCESSIVE')
 215      CONTINUE
        ENDDO
      ENDIF
*----
*  RETURN
*----
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(' ------  OUTPUT FROM ROUTINE = ',A6)
 6001 FORMAT(' --------------------------------------')
 6010 FORMAT(' REGIONAL MERGE ',/
     >       3(' OLD NUNBER -> NEW NUMBER '))
 6011 FORMAT(' SURFACE  MERGE ',/
     >       3(' OLD NUNBER -> NEW NUMBER'))
 6012 FORMAT(' MIXTURE MODIFICATION ',/
     >       3(' OLD REGION -> MIXTURE   '))
 6013 FORMAT(' MIXTURE MODIFICATION ',/
     >       3(' NEW REGION -> MIXTURE   '))
 6014 FORMAT(' ALBEDO MODIFICATION  ',/
     >       3(' SURFACE    -> ALBEDO    '))
 6015 FORMAT(' REGION EXTRACTED FROM TRACK FILE ')
 6020 FORMAT(3(1X,I10,4X,I10))
 6021 FORMAT(3(1X,I10,4X,F10.7))
 6025 FORMAT(6(1X,I10))
      END
