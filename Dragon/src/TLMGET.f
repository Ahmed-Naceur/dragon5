*DECK TLMGET
      SUBROUTINE TLMGET(IPRINT,NPLOTS,NDIM  ,CARLST,
     >                  IPLP  ,DPLP  )
*
*-----------------------------------------------------------------------
*
*Purpose:
* To read from the input file the TLM: module processing options.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
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
* NPLOTS  number of plots generated.
* NDIM    number of dimensions.
*
*Parameters: input/output
* CARLST  last character string read.
*
*Parameters: output
* IPLP    integer plot parameters.
* DPLP    double precision plot parameters.
*
*Comments:
* Input data is of the form:
*    {
*      POINTS [NoPause] |
*      DIRECTIONS [NoPause] DIR idir [ PLAN iplan { U iu | V iv } ] |
*      PLANA [NoPause] A ra B rb [ C rc ] D rd |
*      PLANP [NoPause] DIR idir DIST dist [ PLAN iplan ]
*    }
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      CHARACTER        CARLST*72
      INTEGER          IPRINT,NPLOTS,NDIM
      INTEGER          IPLP(6,NPLOTS)
      DOUBLE PRECISION DPLP(4,NPLOTS)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='TLMGET')
*----
*  Variables for input via REDGET
*----
      INTEGER          ITYPLU,INTLIR
      CHARACTER        CARLIR*72
      REAL             REALIR
      DOUBLE PRECISION DBLLIR
*----
*  Local variables
*----
      INTEGER          IPLOT
*----
*  Processing starts:
*  print routine opening header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
*----
*  Get data from input file
*----
      CALL XDISET(IPLP,6*NPLOTS,0)
      CALL XDDSET(DPLP,4*NPLOTS,0.0)
      ITYPLU=3
      CARLIR=CARLST
      IPLOT=0
 100  CONTINUE
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': Read error -- Character variable expected')
      IF(CARLIR(1:4) .EQ. ';') GO TO 105
      IPLOT=IPLOT+1
      IF(IPLOT .GT. NPLOTS) THEN
        WRITE(IOUT,9000) NAMSBR,NPLOTS
 110    CONTINUE
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .EQ. 3) THEN
          IF(CARLIR(1:4) .EQ. ';') GO TO 105
        ENDIF
        GO TO 110
      ENDIF
      IF(CARLIR .EQ. 'POINTS') THEN
        IPLP(1,IPLOT)=1
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': Read error -- Character variable expected')
        IF(CARLIR .NE. 'NoPause') GO TO 100
        IPLP(1,IPLOT)=-1
      ELSE IF(CARLIR .EQ. 'DIRECTIONS') THEN
        IPLP(1,IPLOT)=2
 120    CONTINUE
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': Read error -- Character variable expected')
        IF(CARLIR .EQ. 'DIR') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- integer variable expected for DIR')
          IPLP(2,IPLOT)=INTLIR
          GO TO 120
        ELSE IF(CARLIR .EQ. 'PLAN') THEN
          IF(NDIM .EQ. 2) WRITE(IOUT,9001) NAMSBR,CARLIR(1:12)
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- integer variable expected for PLAN')
          IPLP(3,IPLOT)=INTLIR
          GO TO 120
        ELSE IF(CARLIR .EQ. 'U') THEN
          IF(NDIM .EQ. 2) WRITE(IOUT,9001) NAMSBR,CARLIR(1:12)
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- integer variable expected for U')
          IPLP(4,IPLOT)=INTLIR
          GO TO 120
        ELSE IF(CARLIR .EQ. 'V') THEN
          IF(NDIM .EQ. 2) WRITE(IOUT,9001) NAMSBR,CARLIR(1:12)
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- integer variable expected for V')
          IPLP(5,IPLOT)=INTLIR
          GO TO 120
        ELSE IF(CARLIR .EQ. 'NoPause') THEN
          IPLP(1,IPLOT)=-2
          GO TO 120
        ELSE IF(CARLIR .EQ. 'SPoints') THEN
          IPLP(6,IPLOT)=1
          GO TO 120
        ELSE
          GO TO 100
        ENDIF
      ELSE IF(CARLIR .EQ. 'PLANA') THEN
        IPLP(1,IPLOT)=3
 130    CONTINUE
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': Read error -- Character variable expected')
        IF(CARLIR .EQ. 'A') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .EQ. 2) THEN
            DPLP(1,IPLOT)=REALIR
          ELSE IF(ITYPLU .EQ. 4) THEN
            DPLP(1,IPLOT)=DBLLIR
          ELSE
            CALL XABORT(NAMSBR//
     >      ': Read error -- real variable expected for A')
          ENDIF
          GO TO 130
        ELSE IF(CARLIR .EQ. 'B') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .EQ. 2) THEN
            DPLP(2,IPLOT)=REALIR
          ELSE IF(ITYPLU .EQ. 4) THEN
            DPLP(2,IPLOT)=DBLLIR
          ELSE
            CALL XABORT(NAMSBR//
     >      ': Read error -- real variable expected for B')
          ENDIF
          GO TO 130
        ELSE IF(CARLIR .EQ. 'C') THEN
          IF(NDIM .EQ. 2) WRITE(IOUT,9001) NAMSBR,CARLIR(1:12)
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .EQ. 2) THEN
            DPLP(3,IPLOT)=REALIR
          ELSE IF(ITYPLU .EQ. 4) THEN
            DPLP(3,IPLOT)=DBLLIR
          ELSE
            CALL XABORT(NAMSBR//
     >      ': Read error -- real variable expected for C')
          ENDIF
          GO TO 130
        ELSE IF(CARLIR .EQ. 'D') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .EQ. 2) THEN
            DPLP(4,IPLOT)=REALIR
          ELSE IF(ITYPLU .EQ. 4) THEN
            DPLP(4,IPLOT)=DBLLIR
          ELSE
            CALL XABORT(NAMSBR//
     >      ': Read error -- real variable expected for D')
          ENDIF
          GO TO 130
        ELSE IF(CARLIR .EQ. 'NoPause') THEN
          IPLP(1,IPLOT)=-3
          GO TO 130
        ELSE
          GO TO 100
        ENDIF
      ELSE IF(CARLIR .EQ. 'PLANP') THEN
        IPLP(1,IPLOT)=4
 140    CONTINUE
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': Read error -- Character variable expected')
        IF(CARLIR .EQ. 'DIR') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- integer variable expected for DIR')
          IPLP(2,IPLOT)=INTLIR
          GO TO 140
        ELSE IF(CARLIR .EQ. 'DIST') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .EQ. 2) THEN
            DPLP(1,IPLOT)=REALIR
          ELSE IF(ITYPLU .EQ. 4) THEN
            DPLP(1,IPLOT)=DBLLIR
          ELSE
            CALL XABORT(NAMSBR//
     >      ': Read error -- real variable expected for DIST')
          ENDIF
          GO TO 140
        ELSE IF(CARLIR .EQ. 'PLAN') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(NDIM .EQ. 2) WRITE(IOUT,9001) NAMSBR,CARLIR(1:12)
          IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- integer variable expected for PLAN')
          IPLP(3,IPLOT)=INTLIR
          GO TO 140
        ELSE IF(CARLIR .EQ. 'NoPause') THEN
          IPLP(1,IPLOT)=-4
          GO TO 140
        ELSE
          GO TO 100
        ENDIF
      ELSE IF(CARLIR .EQ. 'REGIONS') THEN
        IPLP(1,IPLOT)=5
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .EQ. 1) THEN
          IPLP(2,IPLOT)=INTLIR
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          GO TO 100
        ENDIF
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': Read error -- Character variable expected')
        IF(CARLIR .NE. 'NoPause') GO TO 100
        IPLP(1,IPLOT)=-5
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .EQ. 1) THEN
          IPLP(2,IPLOT)=INTLIR
        ENDIF
      ELSE
        CALL XABORT(NAMSBR//': Keyword '//CARLIR(1:12)//' is invalid.')
      ENDIF
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      GO TO 100
 105  CONTINUE
      CARLST=CARLIR
*----
*  Processing finished, return
*----
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
*----
*  Warning formats
*----
 9000 FORMAT(1X,'Warning from ',A6,2X,'Number of plots exceeded'/
     >1X,'Only first ',I10,1X,'plots considered')
 9001 FORMAT(1X,'Warning from ',A6,2X,'Invalid keyword '/
     >1X,'Keyword : ',A12,1X,'Not used in 2-D')
      END
