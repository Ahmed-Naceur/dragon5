*DECK IDET
      SUBROUTINE IDET(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Detector integrated response evaluation
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* NENTRY  number of data structures transfered to this module.
* HENTRY  name of the data structures.
* IENTRY  data structure type where:
*         IENTRY=1 for LCM memory object;
*         IENTRY=2 for XSM file;
*         IENTRY=3 for sequential binary file;
*         IENTRY=4 for sequential ASCII file.
* JENTRY  access permission for the data structure where:
*         JENTRY=0 for a data structure in creation mode;
*         JENTRY=1 for a data structure in modifications mode;
*         JENTRY=2 for a data structure in read-only mode.
* KENTRY  data structure pointer.
*
*Comments:
* The IDET: module specification is:
* IDETEC := IDET: [ IDETEC ] TRKNAM FLUNAM LIBNAM [ FMAP ] :: (descidet) ;
* where
*   IDETEC : name of a \emph{idetect} data structure, (L\_INTDETEC signature) 
*     that will be created or updated by the IDET: module.
*   TRKNAM : name of the read-only \emph{tracking} data structure 
*     (L\_TRACK signature) containing the finite-element tracking.
*   FLUNAM : name of the read-only \emph{fluxunk data structure 
*     (L\_FLUX signature) containing the finite-element solution.
*   LIBNAM : name of the read-only \emph{macrolib} data structure 
*     (L\_LIBRARY signature) that contains the interpolated microscopic
*     cross sections.
*   FMAP   : name of the read-only  \emph{fmap} data structure
*     (L\_MAP signature) containing renumbered mixture indices. This object 
*     is optionnal.
*   (descidet) : structure describing the input data to the IDET: module.
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
      INTEGER MAXCO
      PARAMETER (MAXCO=300,MAXNI=10,NSTATE=40)
      INTEGER INDIC,NITMA,ISTATE(NSTATE)
      DOUBLE PRECISION DFLOT
      CHARACTER CMODUL*12,HSIGN*12,TEXT12*12,DETNAM*12,REANAM*12
      REAL FLOT
      TYPE(C_PTR) IPIDET,IPTRK,IPFLU,IPLIB,IPMAP
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, DIMENSION(:), ALLOCATABLE :: NINX,NINY,NINZ
      REAL, DIMENSION(:), ALLOCATABLE :: DETECT
      REAL, DIMENSION(:,:), ALLOCATABLE :: COORD1,COORD2,COORD3
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(COORD1(MAXNI,MAXCO),COORD2(MAXNI,MAXCO),
     > COORD3(MAXNI,MAXCO),NINX(MAXCO),NINY(MAXCO),NINZ(MAXCO))
*----
*  PARAMETER VALIDATION
*----
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('IDET: LCM'
     > //' object expected at LHS.')
      IF(JENTRY(1).EQ.2) CALL XABORT('IDET: L_INTDETEC entry in create'
     > //' or modification mode expected.')
      IPIDET=KENTRY(1)
      IF(JENTRY(1).EQ.0) THEN
        HSIGN='L_INTDETEC'
        CALL LCMPTC(IPIDET,'SIGNATURE',12,1,HSIGN)
        DETNAM='U235'
        REANAM='NFTOT'
      ELSE IF(JENTRY(1).EQ.1) THEN
        CALL LCMGTC(IPIDET,'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_INTDETEC') THEN
          TEXT12=HENTRY(3)
          CALL XABORT('IDET: signature of '//TEXT12//' IS '//HSIGN//
     >    '. L_INTDETEC expected.')
        ENDIF
        CALL LCMGET(IPIDET,'STATE-VECTOR',ISTATE)
        IF(ISTATE(1).NE.MAXNI) CALL XABORT('IDET: invalid MAXNI.')
        NDETC=ISTATE(2)
        IF(NDETC.GT.MAXCO) CALL XABORT('IDET: MAXCO overflow.')
        CALL LCMGET(IPIDET,'NINX',NINX)
        CALL LCMGET(IPIDET,'NINY',NINY)
        CALL LCMGET(IPIDET,'NINZ',NINZ)
        CALL LCMGET(IPIDET,'COORD1',COORD1)
        CALL LCMGET(IPIDET,'COORD2',COORD2)
        CALL LCMGET(IPIDET,'COORD3',COORD3)
        CALL LCMGTC(IPIDET,'DETNAM',12,1,DETNAM)
        CALL LCMGTC(IPIDET,'REANAM',12,1,REANAM)
      ENDIF
      IPFLU=C_NULL_PTR
      IPTRK=C_NULL_PTR
      IPLIB=C_NULL_PTR
      IPMAP=C_NULL_PTR
      CMODUL=' '
      DO I=2,NENTRY
        IF(IENTRY(I).GT.2) CALL XABORT('IDET: LCM object expected.')
        IF(JENTRY(I).NE.2) CALL XABORT('IDET: LCM object in read-only '
     > //'MODE EXPECTED AT RHS.')
        CALL LCMGTC(KENTRY(I),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.EQ.'L_FLUX') THEN
          IPFLU=KENTRY(I)
        ELSEIF(HSIGN.EQ.'L_TRACK') THEN
          IPTRK=KENTRY(I)
          CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,CMODUL)
        ELSEIF(HSIGN.EQ.'L_LIBRARY') THEN
          IPLIB=KENTRY(I)
        ELSEIF(HSIGN.EQ.'L_MAP') THEN
          IPMAP=KENTRY(I)
        ELSE
          TEXT12=HENTRY(I)
          CALL XABORT('IDET: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     >    '. L_FLUX, L_TRACK or L_LIBRARY expected.')
        ENDIF
      ENDDO
      IF(CMODUL.NE.'TRIVAC') CALL XABORT('IDET: TRIVAC tracking expect'
     > //'ed.')
*----
*  READ INPUTS
*----
      IMPX=1
      NDETC=0
   10 CALL REDGET(INDIC,NITMA,FLOT,TEXT12,DFLOT)
      IF(INDIC.NE.3) CALL XABORT('IDET: character data expected.')
      IF(TEXT12.EQ.'EDIT') THEN
        CALL REDGET(INDIC,IMPX,FLOT,TEXT12,DFLOT)
        IF(INDIC.NE.1) CALL XABORT('IDET: integer data expected.')
      ELSE IF(TEXT12.EQ.'DETNAME') THEN
        CALL REDGET(INDIC,NITMA,FLOT,DETNAM,DFLOT)
        IF(INDIC.NE.3) CALL XABORT('IDET: character data expected(1).')
      ELSE IF(TEXT12.EQ.'REANAME') THEN
        CALL REDGET(INDIC,NITMA,FLOT,REANAM,DFLOT)
        IF(INDIC.NE.3) CALL XABORT('IDET: character data expected(2).')
      ELSE IF(TEXT12.EQ.'DETECTOR') THEN
   20   CALL REDGET(INDIC,NITMA,FLOT,TEXT12,DFLOT)
        IF(INDIC.NE.3) CALL XABORT('IDET: character data expected.')
   30   IF(TEXT12.EQ.'POSITION') THEN
*         Cartesian position of a single detector
          NDETC=NDETC+1
          IF(NDETC.GT.MAXCO) CALL XABORT('IDET: MAXCO overflow.')
          NINX(NDETC)=1
          NINY(NDETC)=1
          NINZ(NDETC)=1
          CALL REDGET(INDIC,NITMA,FLOT,TEXT12,DFLOT)
          IF(INDIC.EQ.2) THEN
            COORD1(1,NDETC)=FLOT
          ELSE IF((INDIC.EQ.3).AND.(TEXT12.EQ.'INTEG')) THEN
            NINX(NDETC)=MAXNI
            CALL REDGET(INDIC,NITMA,COO1,TEXT12,DFLOT)
            IF(INDIC.NE.2) CALL XABORT('IDET: COORD1 data1 expected.')
            CALL REDGET(INDIC,NITMA,COO2,TEXT12,DFLOT)
            IF(INDIC.NE.2) CALL XABORT('IDET: COORD1 data2 expected.')
            IF(COO2.LE.COO1) CALL XABORT('IDET: COORD1 data2<=data1.')
            DELTA=(COO2-COO1)/REAL(MAXNI-1)
            DO INX=1,MAXNI
              COORD1(INX,NDETC)=COO1+REAL(INX-1)*DELTA
            ENDDO
          ELSE
            CALL XABORT('IDET: COORD1 data or INTEG keyword expected.')
          ENDIF
          CALL REDGET(INDIC,NITMA,FLOT,TEXT12,DFLOT)
          IF(INDIC.EQ.2) THEN
            COORD2(1,NDETC)=FLOT
          ELSE IF((INDIC.EQ.3).AND.(TEXT12.EQ.'INTEG')) THEN
            NINY(NDETC)=MAXNI
            CALL REDGET(INDIC,NITMA,COO1,TEXT12,DFLOT)
            IF(INDIC.NE.2) CALL XABORT('IDET: COORD2 data1 expected.')
            CALL REDGET(INDIC,NITMA,COO2,TEXT12,DFLOT)
            IF(INDIC.NE.2) CALL XABORT('IDET: COORD2 data2 expected.')
            IF(COO2.LE.COO1) CALL XABORT('IDET: COORD2 data2<=data1.')
            DELTA=(COO2-COO1)/REAL(MAXNI-1)
            DO INY=1,MAXNI
              COORD2(INY,NDETC)=COO1+REAL(INY-1)*DELTA
            ENDDO
          ELSE
            CALL XABORT('IDET: COORD2 data or INTEG keyword expected.')
          ENDIF
          CALL REDGET(INDIC,NITMA,FLOT,TEXT12,DFLOT)
          IF(INDIC.EQ.2) THEN
            COORD3(1,NDETC)=FLOT
            GO TO 20
          ELSE IF(INDIC.EQ.3) THEN
            IF(TEXT12.EQ.'INTEG') THEN
              NINZ(NDETC)=MAXNI
              CALL REDGET(INDIC,NITMA,COO1,TEXT12,DFLOT)
              IF(INDIC.NE.2) CALL XABORT('IDET: COORD3 data1 expected.')
              CALL REDGET(INDIC,NITMA,COO2,TEXT12,DFLOT)
              IF(INDIC.NE.2) CALL XABORT('IDET: COORD3 data2 expected.')
              IF(COO2.LE.COO1) CALL XABORT('IDET: COORD3 data2<=data1.')
              DELTA=(COO2-COO1)/REAL(MAXNI-1)
              DO INZ=1,MAXNI
                COORD3(INZ,NDETC)=COO1+REAL(INZ-1)*DELTA
              ENDDO
              GO TO 20
            ELSE
              COORD3(1,NDETC)=1.0
              GO TO 30
            ENDIF
          ELSE
            CALL XABORT('IDET: real or character data expected.')
          ENDIF
        ELSE IF(TEXT12.EQ.'ENDD') THEN
          GO TO 10
        ELSE
          CALL XABORT('IDET: POSITION, MIXTURE or ENDP keyword expec'
     >    //'ted.')
        ENDIF
        GO TO 20
      ELSE IF(TEXT12.EQ.';') THEN
        GO TO 40
      ELSE
        CALL XABORT('IDET: unknownn keyword-->'//TEXT12)
      ENDIF
      GO TO 10
*----
*  PERFORM FLUX INTERPOLATION OVER DETECTOR LOCATIONS
*----
   40 IF(NDETC.EQ.0) CALL XABORT('IDET: no detector defined.')
      ALLOCATE(DETECT(NDETC))
      CALL IDET01(IPTRK,IPFLU,IPLIB,IPMAP,IMPX,NDETC,MAXNI,NINX,NINY,
     > NINZ,COORD1,COORD2,COORD3,DETNAM,REANAM,DETECT)
*----
*  PRINT DETECTOR RESPONSE
*----
      IF(IMPX.GT.0) THEN
        WRITE(6,'(/25H DET: DETECTOR READINGS (,2A12,1H))') DETNAM,
     >  REANAM
        WRITE(6,'(10X,8HDETECTOR,5X,7HREADING)')
        DO I=1,NDETC
          WRITE(6,'(8X,I10,1P,E16.5)') I,DETECT(I)
        ENDDO
      ENDIF
*----
*  SAVE DETECTOR INFORMATION ON LCM
*----
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=MAXNI
      ISTATE(2)=NDETC
      CALL LCMPUT(IPIDET,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPIDET,'NINX',NDETC,1,NINX)
      CALL LCMPUT(IPIDET,'NINY',NDETC,1,NINY)
      CALL LCMPUT(IPIDET,'NINZ',NDETC,1,NINZ)
      CALL LCMPUT(IPIDET,'COORD1',MAXNI*NDETC,2,COORD1)
      CALL LCMPUT(IPIDET,'COORD2',MAXNI*NDETC,2,COORD2)
      CALL LCMPUT(IPIDET,'COORD3',MAXNI*NDETC,2,COORD3)
      CALL LCMPTC(IPIDET,'DETNAM',12,1,DETNAM)
      CALL LCMPTC(IPIDET,'REANAM',12,1,DETNAM)
      CALL LCMPUT(IPIDET,'RESPON',NDETC,2,DETECT)
*----
*  RELEASE MEMORY
*----
      DEALLOCATE(DETECT,NINZ,NINY,NINX,COORD3,COORD2,COORD1)
      RETURN
      END
