*DECK EPC
      SUBROUTINE EPC(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Process error propagation parameters.
*
*Copyright:
* Copyright (C) 2009 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* NENTRY  number of data structures transfered to this module.
* HENTRY  name of the data structures.
* IENTRY  data structure type where:
*         =1 for LCM memory object;
*         =2 for XSM file;
*         =3 for sequential binary file;
*         =4 for sequential ASCII file.
* JENTRY  access permission for the data structure where:
*         =0 for a data structure in creation mode;
*         =1 for a data structure in modifications mode;
*         =2 for a data structure in read-only mode.
* KENTRY  data structure pointer.
*
*Comments:
* Instructions for the use of the EPC: module:
*   Param := EPC: [ Param ] [[ ParamDist ]]  :: (EPCpara) (EPCget) ;
* where
*     Param     : parameter data structure.
*     ParamDist : sequential binary/ASCII parameter distributions
*     (EPCpara) : PARA keyword processing options (routine EPCPAR)
*     (EPCget)  : GET keyword processing options (routine EPCGET)
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)      KENTRY(NENTRY)
      CHARACTER        HENTRY(NENTRY)*12
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='EPC   ')
      INTEGER          ILCMUP,ILCMDN,MXFIL,MXOPT
      PARAMETER       (ILCMUP=1,ILCMDN=2,MXFIL=20,MXOPT=20)
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
*----
*  Input and output parameters
*----
      INTEGER          ITYPLU,INTLIR
      CHARACTER        CARLIR*12
      REAL             REALIR
      DOUBLE PRECISION DBLLIR
*----
*  Local variables
*----
      CHARACTER*12     CARRET,SENTRY(MXFIL),HET
      INTEGER          IEN
      CHARACTER        HSIGN*12
      INTEGER          IPRINT,NOPT,IOPT(MXOPT)
      INTEGER          ISTATE(NSTATE)
*----
*  Validate entry parameters
*----
      IF(NENTRY .GT. MXFIL) CALL XABORT(NAMSBR//
     >  ': Too many files or data structures for this module.')
*----
*  Scan data structure to determine signature (input or update)
*----
      DO IEN=1,NENTRY
        SENTRY(IEN)='            '
        IF(IENTRY(IEN) .EQ. 1 .OR. IENTRY(IEN) .EQ. 2) THEN
          IF(JENTRY(IEN) .NE. 0) THEN
            CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
            SENTRY(IEN)=HSIGN
          ENDIF
        ENDIF
      ENDDO
*----
*  Read EDIT and main option
*  Default option is NDIST option
*----
      CARRET=';'
      NOPT=MXOPT
      CALL XDISET(IOPT,NOPT,0)
      IPRINT=1
 100  CONTINUE
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      IF(ITYPLU .EQ. 10) GO TO 105
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': Read error -- Character variable expected')
      CARRET=CARLIR
      IF(CARLIR .EQ. ';') THEN
        GO TO 105
      ELSE IF(CARLIR .EQ. 'EDIT') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- integer value for EDIT expected.')
        IPRINT=INTLIR
      ELSE IF(CARLIR .EQ. 'RNDPhysParam') THEN
        IOPT(1)=1
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': Read error -- Character variable expected')
        CARRET=CARLIR
        GO TO 105
      ELSE IF(CARLIR .EQ. 'GPTPhysParam') THEN
        IOPT(1)=2
        GO TO 105
      ELSE IF(CARLIR .EQ. 'RNDMicXS') THEN
        IOPT(1)=3
        GO TO 105
      ELSE IF(CARLIR .EQ. 'GPTMicXS') THEN
        IOPT(1)=4
        GO TO 105
      ELSE
        CALL XABORT(NAMSBR//': Read error -- keyword'//CARLIR//
     >'for processing option is invalid.')
      ENDIF
      GO TO 100
 105  CONTINUE
*----
*  Analyze structures
*  Find the first L_EPC structure or the first new structure
*----
      DO IEN=1,NENTRY
        HET=HENTRY(IEN)
        IF(IENTRY(IEN) .EQ. 1 .OR. IENTRY(IEN) .EQ. 2) THEN
          IF(JENTRY(IEN) .EQ. 1) THEN
            IF(SENTRY(IEN) .EQ. 'L_EPC       ') THEN
              CALL LCMGET(KENTRY(IEN),'STATE-VECTOR',ISTATE)
              IF(IOPT(1) .NE. ISTATE(1)) CALL XABORT(NAMSBR//
     >': Structure L_EPC :'//HET//
     >' not compatible with processing option.')
              IOPT(2)=IEN
              IOPT(4)=ISTATE(2)
              GO TO 120
            ENDIF
          ELSE
            IOPT(2)=IEN
            SENTRY(IEN)='L_EPC       '
            ISTATE(1)=IOPT(1)
            HSIGN=SENTRY(IEN)
            CALL LCMPTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
            CALL LCMPUT(KENTRY(IEN),'STATE-VECTOR',NSTATE,1,ISTATE)
            IF(ISTATE(1) .EQ. 0) THEN
              HSIGN='            '
              CALL LCMPTC(KENTRY(IEN),'ParametreNom',12,1,HSIGN)
              CALL LCMPUT(KENTRY(IEN),'ParametreNbr',1,1,0)
              CALL LCMPUT(KENTRY(IEN),'ParametreNxt',1,1,0)
              CALL LCMPUT(KENTRY(IEN),'ParametreRef',1,2,0.0)
            ENDIF
            GO TO 120
          ENDIF
        ENDIF
      ENDDO
      CALL XABORT(NAMSBR//': No structure found for L_EPC.')
 120  CONTINUE
*----
*  Find the first read_only SEQ_ASCII or BINARY file
*----
      DO IEN=1,NENTRY
        IF(JENTRY(IEN) .EQ. 2) THEN
          IF(IENTRY(IEN) .EQ. 3) THEN
            IOPT(3)=-IEN
            GO TO 130
          ELSE IF(IENTRY(IEN) .EQ. 4) THEN
            IOPT(3)=IEN
            GO TO 130
          ENDIF
        ENDIF
      ENDDO
 130  CONTINUE
*----
*  Process option
*----
      IF(IOPT(1) .EQ. 1) THEN
*----
*  Option RNDPhysParam
*----
        CALL EPCRPD(NENTRY,KENTRY,IPRINT,NOPT,IOPT,CARRET)
      ELSE IF(IOPT(1) .EQ. 2) THEN
*----
*  Option GPTPhysParam
*----
        CALL XABORT(NAMSBR//' Option GPTPhysParam not programmed yet')
*        CALL EPCGPD(NENTRY,KENTRY,IPRINT,NOPT,IOPT,CARRET)
      ELSE IF(IOPT(1) .EQ. 3) THEN
*----
*  Locate microlib
*----
        DO IEN=1,NENTRY
          IF(IENTRY(IEN) .EQ. 1 .OR. IENTRY(IEN) .EQ. 2) THEN
            IF(JENTRY(IEN) .EQ. 1) THEN
              CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
              IF(HSIGN .EQ. 'L_LIBRARY   ') THEN
                IOPT(4)=IEN
                GO TO 140
              ENDIF
            ENDIF
          ENDIF
        ENDDO
*        CALL XABORT(NAMSBR//' Option RNDMicXS requires a microlib')
 140    CONTINUE
*----
*  Process
*----
        CALL EPCRMD(NENTRY,KENTRY,IPRINT,NOPT,IOPT,CARRET)
      ELSE IF(IOPT(1) .EQ. 4) THEN
*----
*  Option GPTMicXS
*----
        CALL XABORT(NAMSBR//' Option GPTMicXS not programmed yet')
*        CALL EPCGMD(NENTRY,KENTRY,IPRINT,NOPT,IOPT,CARRET)
      ENDIF
*----
*  Processing finished, return
*----
      RETURN
*----
*  Warning formats
*----
      END
