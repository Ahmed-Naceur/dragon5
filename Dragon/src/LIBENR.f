*DECK LIBENR
      SUBROUTINE LIBENR(CFILNA,IVERW,MAXR,NEL,ITNAM,KPAX,BPAX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read depletion data on a WIMA-D4 or WIMSE formatted library.
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
* CFILNA  WIMS-D4 or WIMS-E file name.
* IVERW   type of file (=4: WIMS-D4; =5: WIMS-E).
* MAXR    number of reaction types.
* NEL     number of isotopes on library.
*
*Parameters: output
* ITNAM   reactive isotope names in chain.
* KPAX    complete reaction type matrix.
* BPAX    complete branching ratio matrix.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER CFILNA*8
      INTEGER IVERW,MAXR,NEL,ITNAM(3,NEL),KPAX(NEL+MAXR,NEL)
      REAL BPAX(NEL+MAXR,NEL)
*----
*  INTERNAL PARAMETERS
*   CONVE  : ENERGY CONVERSION FACTOR FROM JOULES/(MOLES*10**-24)
*            TO MEV/NUCLIDE = 1.03643526E+13
*   CONVD  : DECAY CONSTANT CONVERSION FACTOR FROM S**(-1) TO
*            10**(-8)*S**(-1) = 1.0+8
*----
      INTEGER      KCAPTU,KDECAY,KFISSP
      REAL         CONVE,CONVD
      PARAMETER   (KCAPTU=3,KDECAY=1,KFISSP=2,
     >             CONVE=1.03643526E+13,CONVD=1.0E+8)
      CHARACTER    TEXT8*8
*----
*  WIMS-D4 LIBRARY PARAMETERS
*   IUTYPE : TYPE OF FILE = 2 (BINARY)
*   LRIND  : LENGHT RECORD ON DA FILE = 0
*   IACTO  : OPEN ACTION = 2 (READ ONLY)
*   IACTC  : CLOSE ACTION = 2 (KEEP)
*   MAXISO : MAX. NB. ISOTOPE = 246
*   MLDEP  : MAXIMUM NUMBER OF REACTION PER
*            ISOTOPE IN WIMS-D4 = MAXISO+4
*   LPZ    : LENGTH OF WIMS PARAMETER ARRAY = 8
*   NPZ    : LIST OF MAIN PARAMETERS
*   IWISO  : ID OF ISOTOPE
*   IBURN  : INTEGER BURNUP DATA
*   RBURN  : REAL BURNUP DATA
*----
      INTEGER      IUTYPE,LRIND,IACTO,IACTC,MAXISO,MLDEP,LPZ
      PARAMETER   (IUTYPE=2,LRIND=0,IACTO=2,IACTC=1,MAXISO=246,
     >             MLDEP=MAXISO+4,LPZ=8)
      INTEGER      NPZ(LPZ),IWISO(MAXISO),IBURN(MLDEP)
      REAL         RBURN(MLDEP),RTEMP
*----
*  EXTERNAL FUNCTIONS
*----
      INTEGER      KDROPN,LIBWID,KDRCLS
*----
*  LOCAL VARIABLES
*----
      INTEGER      IUNIT,II,J,ISO,JC,JB,JSO,IT,IERR
*----
*  OPEN WIMS-D4 OR WIMSE LIBRARY
*  READ GENERAL DIMENSIONING
*  READ ISOTOPE ID NUMBER AND CREATE EQUIVALENT ISOTOPE NAME
*----
      IUNIT=KDROPN(CFILNA,IACTO,IUTYPE,LRIND)
      IF(IUNIT.LE.0) CALL XABORT('LIBENR: WIMS-D4 LIBRARY '//
     >    CFILNA//' CANNOT BE OPENED FOR DEPLETION')
      READ(IUNIT) (NPZ(II),II=1,LPZ)
      IF(NPZ(1).NE.NEL) CALL XABORT('LIBENR: TOO MANY ISOTOPES '//
     >    'ON WIMS-D4 LIBRARY'//CFILNA)
      READ(IUNIT) (IWISO(J),J=1,NEL)
      DO 10 ISO=1,NEL
        TEXT8='        '
        IF     (IWISO(ISO).LT.10) THEN
          WRITE(TEXT8,'(I1)') IWISO(ISO)
        ELSE IF(IWISO(ISO).LT.100) THEN
          WRITE(TEXT8,'(I2)') IWISO(ISO)
        ELSE IF(IWISO(ISO).LT.1000) THEN
          WRITE(TEXT8,'(I3)') IWISO(ISO)
        ELSE IF(IWISO(ISO).LT.10000) THEN
          WRITE(TEXT8,'(I4)') IWISO(ISO)
        ELSE IF(IWISO(ISO).LT.100000) THEN
          WRITE(TEXT8,'(I5)') IWISO(ISO)
        ELSE IF(IWISO(ISO).LT.1000000) THEN
          WRITE(TEXT8,'(I6)') IWISO(ISO)
        ELSE IF(IWISO(ISO).LT.10000000) THEN
          WRITE(TEXT8,'(I7)') IWISO(ISO)
        ELSE IF(IWISO(ISO).LT.100000000) THEN
          WRITE(TEXT8,'(I8)') IWISO(ISO)
        ENDIF
        READ(TEXT8,'(2A4)') ITNAM(1,ISO),ITNAM(2,ISO)
 10   CONTINUE
*---
*  READ TWO ADDITIONAL RECORDS BEFORE DEPLETION DATA
*----
      READ(IUNIT) (RTEMP,J=1,NPZ(2)+1)
      IF(IVERW.EQ.4) READ(IUNIT) (RTEMP,J=1,NPZ(3))
*----
*  READ DEPLETION CHAIN FOR EACH ISOTOPES
*----
      DO 100 ISO=1,NEL
        RBURN(1)=0.0
        READ(IUNIT) JC,IBURN(1),
     >    (RBURN(JB),IBURN(JB),JB=2,JC/2)
        IF(JC/2.GT.MLDEP) CALL XABORT('LIBENR: MLDEP OVERFLOW.')
*----
*  CAPTURE -> RBURN(2) > ALWAYS PRESENT
*   IF ISOTOPE RESULTING FROM CAPTURE IS KNOWN STORE IN ADEQUATE
*   POSITION ELSE STORE IN NEL+1
*  DECAY   -> RBURN(3) > 0.0
*   IF ISOTOPE RESULTING FROM DECAY IS KNOWN STORE IN ADEQUATE
*   POSITION ELSE STORE IN NEL+2
*  FISSILE -> IBURN(4) > 1
*   JC=8 -> ISOTOPE RESULTING FROM FISSION NOT KNOWN STORE IN NEL+3
*   JC>8 -> ISOTOPE RESULTING FROM FISSION KNOWN STORE IN ADEQUATE
*   POSITION
*----
        IF(JC.GE.8) THEN
*         radiative capture, always present
          JSO=LIBWID(NEL,IWISO,IBURN(2))
          IF(JSO.GT.0) THEN
            IF(KPAX(JSO,ISO) .EQ. 0) THEN
              KPAX(JSO,ISO)=KCAPTU
              BPAX(JSO,ISO)=RBURN(2)
              KPAX(NEL+KCAPTU,JSO)=1
            ENDIF
          ENDIF
          KPAX(NEL+KCAPTU,ISO)=1
*
*         radioactive decay, optionnal
          IF(RBURN(3).GT.0.0) THEN
            JSO=LIBWID(NEL,IWISO,IBURN(3))
            IF(JSO.GT.0) THEN
              IF(KPAX(JSO,ISO) .EQ. 0) THEN
                KPAX(JSO,ISO)=KDECAY
                BPAX(JSO,ISO)=1.0
                KPAX(NEL+KCAPTU,JSO)=1
              ENDIF
            ENDIF
            KPAX(NEL+KDECAY,ISO)=1
            BPAX(NEL+KDECAY,ISO)=RBURN(3)*CONVD
          ENDIF
*
*         fission energy, optionnal
          IF(IBURN(4).GT.1) THEN
            KPAX(NEL+KFISSP,ISO)=1
            BPAX(NEL+KFISSP,ISO)=RBURN(4)*CONVE
          ENDIF
*
*         fission yields and non-fission energy, optionnal
          DO 102 IT=5,JC/2
            IF(IBURN(IT).EQ.-1) THEN
*             radiative capture energy, extension to the WIMS-D4 and
*             WIMS-E specifications
              BPAX(NEL+KCAPTU,ISO)=RBURN(IT)*CONVE
            ELSE IF(IBURN(IT).EQ.-2) THEN
*             radioactive decay energy, extension to the WIMS-D4 and
*             WIMS-E specifications
              BPAX(NEL+KDECAY,ISO)=RBURN(IT)*CONVE
            ELSE IF(RBURN(IT).GT.0.0) THEN
*             fission yields
              JSO=LIBWID(NEL,IWISO,IBURN(IT))
              IF(JSO.GT.0) THEN
                IF(KPAX(JSO,ISO) .EQ. 0) THEN
                  KPAX(JSO,ISO)=KFISSP
                  BPAX(JSO,ISO)=RBURN(IT)
                  KPAX(NEL+KFISSP,JSO)=-1
                  KPAX(NEL+KCAPTU,JSO)=1
                ENDIF
              ENDIF
            ENDIF
 102      CONTINUE
        ENDIF
 100  CONTINUE
*----
*  CLOSE WIMS-D4 OR WIMSE LIBRARY
*----
      IERR=KDRCLS(IUNIT,IACTC)
      IF(IERR.LT.0)
     >  CALL XABORT('LIBENR: WIMS LIBRARY '//CFILNA//
     >    ' CANNOT BE CLOSED')
*----
*  RETURN
*----
      RETURN
      END
