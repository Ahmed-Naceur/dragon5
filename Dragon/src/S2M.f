*DECK S2M
      SUBROUTINE S2M(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
*  To recover information from a SERPENT output file and translate the
*  requested data towards a macrolib.
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal
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
*         HENTRY(1) create or modification ascii file containing
*         Apotrim data;
*         HENTRY(2) read-only type(L_MACROLIB).
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
      TYPE(C_PTR) IPMAC,JPMAC,KPMAC
      PARAMETER (NSTATE=40,IOUT=6)
      CHARACTER TEXT12*12,HSIGN*12,HLINE*512,CM*2
      LOGICAL LB1,LFIS
      DOUBLE PRECISION DFLOTT
      INTEGER ISTATE(NSTATE)
      PARAMETER(CM='00')
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: XS,FLUX,FISS,CHI,DIFF,GAR3
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SCAT
*
      IF(NENTRY.LE.1) CALL XABORT('S2M: MINIMUM OF 2 OBJECTS EXPECTED.')
      TEXT12=HENTRY(1)
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.4)) CALL XABORT('S2M: LCM '
     1 //'OBJECT EXPECTED AT LHS.')
      IF(JENTRY(1).EQ.2) CALL XABORT('S2M: LCM OBJECT IN CREATE OR MOD'
     1 //'IFICATION MODE EXPECTED AT LHS.')
      IPMAC=KENTRY(1)
      IF(IENTRY(2).NE.4) CALL XABORT('S2M: ASCII FILE NAMED '//TEXT12
     1 //' EXPECTED AT RHS.')
      IF(JENTRY(2).NE.2) CALL XABORT('S2M: ASCII FILE IN READ-ONLY MOD'
     1 //'E EXPECTED AT RHS.')
      IFIN=FILUNIT(KENTRY(2))
      HSIGN='L_MACROLIB'
      CALL LCMPTC(IPMAC,'SIGNATURE',12,1,HSIGN)
*----
*  READ THE INPUT DATA
*----
      IMPX=1
      IDX=1
      NGRP=0
      LB1=.FALSE.
   20 CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('S2M: CHARACTER DATA EXPECTED(1).')
      IF(TEXT12.EQ.'EDIT') THEN
*        READ THE PRINT INDEX.
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('S2M: INTEGER DATA EXPECTED(1).')
      ELSE IF(TEXT12.EQ.'IDX') THEN
*        READ THE INSTANCE INDEX.
         CALL REDGET(INDIC,IDX,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('S2M: INTEGER DATA EXPECTED(2).')
         IF(IDX.LE.0) CALL XABORT('S2M: INVALID VALUE OF IDX.')
      ELSE IF(TEXT12.EQ.'B1') THEN
         LB1=.TRUE.
      ELSE IF(TEXT12.EQ.';') THEN
         GO TO 30
      ELSE
         CALL XABORT('S2M: '//TEXT12//' IS AN INVALID KEYWORD.')
      ENDIF
      GO TO 20
*----
*  FIND THE NUMBER OF ENERGY GROUPS
*----
   30 DO
        READ(IFIN,'(A)',END=40) HLINE
        IND1=INDEX(HLINE,'GC_NE                     (idx, 1)')
        IF(IND1.GT.0) THEN
           IND2=INDEX(HLINE,';')
           READ(HLINE(45:IND2-1),'(I5)') NGRP
           GO TO 50
        ENDIF
      ENDDO
   40 CALL XABORT('S2M: UNABLE TO FING NUMBER OF ENERGY GROUPS.')
   50 IF(IMPX.GT.0) WRITE(IOUT,100) NGRP
*----
*  DETERMINE IF THE ISOTOPE IS FISSILE
*----
      LFIS=.FALSE.
      DO
        READ(IFIN,'(A)',END=40) HLINE
        IND1=INDEX(HLINE,'CHI')
        IF(IND1.GT.0) THEN
           LFIS=.TRUE.
           GO TO 60
        ENDIF
      ENDDO
   60 IF(IMPX.GT.0) WRITE(IOUT,110) LFIS
*----
*  RECOVER CROSS SECTIONS
*----
      ALLOCATE(XS(NGRP+1),FLUX(NGRP+1),CHI(NGRP),FISS(NGRP+1),
     > SCAT(NGRP,NGRP),DIFF(NGRP),GAR3(NGRP*NGRP))
      CALL S2MGET(IFIN,'GC_BOUNDS',IDX,.FALSE.,NGRP+1,XS)
      DO IGRP=1,NGRP+1
        XS(IGRP)=XS(IGRP)*1.0E6
      ENDDO
      CALL LCMPUT(IPMAC,'ENERGY',NGRP+1,2,XS)
      CALL S2MGET(IFIN,'TOTXS',IDX,.TRUE.,NGRP+1,XS)
      CALL S2MGET(IFIN,'FLUX',IDX,.TRUE.,NGRP+1,FLUX)
      CALL S2MGET(IFIN,'GTRANSFXS',IDX,.TRUE.,NGRP*NGRP,SCAT) ! I -> J
      IF(LFIS) THEN
        CALL S2MGET(IFIN,'CHI  ',IDX,.TRUE.,NGRP,CHI)
        CALL S2MGET(IFIN,'NSF  ',IDX,.TRUE.,NGRP+1,FISS)
      ENDIF
      IF(LB1) THEN
        CALL S2MGET(IFIN,'B1_DIFFCOEF',IDX,.TRUE.,NGRP,DIFF)
      ENDIF
      JPMAC=LCMLID(IPMAC,'GROUP',NGRP)
      DO IGRP=1,NGRP
        KPMAC=LCMDIL(JPMAC,IGRP)
        CALL LCMPUT(KPMAC,'NTOT0',1,2,XS(IGRP+1))
        CALL LCMPUT(KPMAC,'NWT0',1,2,FLUX(IGRP+1))
        IF(LFIS) THEN
          CALL LCMPUT(KPMAC,'CHI',1,2,CHI(IGRP))
          CALL LCMPUT(KPMAC,'NUSIGF',1,2,FISS(IGRP+1))
        ENDIF
        IF(LB1) THEN
          CALL LCMPUT(KPMAC,'DIFF',1,2,DIFF(IGRP))
        ENDIF
        IPOSDE=0
        IPOS=1
        IGMIN=IGRP
        IGMAX=IGRP
        DO JGRP=1,NGRP
          IF(SCAT(JGRP,IGRP).NE.0.0) THEN
            IGMIN=MIN(IGMIN,JGRP)
            IGMAX=MAX(IGMAX,JGRP)
          ENDIF
        ENDDO
        IJJ=IGMAX
        NJJ=IGMAX-IGMIN+1
        DO JGRP=IGMAX,IGMIN,-1
          IPOSDE=IPOSDE+1
          GAR3(IPOSDE)=SCAT(JGRP,IGRP)
        ENDDO
        CALL LCMPUT(KPMAC,'SCAT'//CM,IPOSDE,2,GAR3)
        CALL LCMPUT(KPMAC,'SIGW'//CM,1,2,SCAT(IGRP,IGRP))
        CALL LCMPUT(KPMAC,'NJJS'//CM,1,1,NJJ)
        CALL LCMPUT(KPMAC,'IJJS'//CM,1,1,IJJ)
        CALL LCMPUT(KPMAC,'IPOS'//CM,1,1,IPOS)
      ENDDO
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NGRP
      ISTATE(2)=1
      ISTATE(3)=1
      IF(LFIS) ISTATE(4)=1
      CALL LCMPUT(IPMAC,'STATE-VECTOR',NSTATE,1,ISTATE)
      DEALLOCATE(GAR3,DIFF,SCAT,FISS,CHI,FLUX,XS)
      RETURN
*
  100 FORMAT(/30H S2M: NUMBER OF ENERGY GROUPS=,I5)
  110 FORMAT(/19H S2M: FISSILE FLAG=,L1)
      END
