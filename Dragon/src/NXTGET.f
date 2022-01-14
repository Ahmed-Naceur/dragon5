*DECK NXTGET
      SUBROUTINE NXTGET(NSTATE,IPRINT,TITLE,ISTATU,RSTATU,NBSLIN,IQUA10,
     >           IBIHET)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To read from the input file the NXT: module processing options.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* NSTATE  dimensions of tracking state vectors.
*
*Parameters: input/output
* IPRINT  print level.
* TITLE   execution title.
* ISTATU  integer parameters for tracking:
*         ISTATU( 1) is the number of regions;
*         ISTATU( 2) is the number of unknown;
*         ISTATU( 3) is the leakage flag;
*         ISTATU( 4) is the maximum number of mixture used;
*         ISTATU( 5) is the number of outer surfaces;
*         ISTATU( 6) is the flux anisotropy order;
*         ISTATU( 7) is the tracking option used;
*         ISTATU( 8) is the track normalization option;
*         ISTATU( 9) is the type of tracks considered;
*         ISTATU(10) is the CP calculation option;
*         ISTATU(11) is the azimuthal quadrature level;
*         ISTATU(12) is the symmetry option;
*         ISTATU(13) is the polar quadrature type;
*         ISTATU(14) is the polar quadrature level;
*         ISTATU(15) is the azimuthal quadrature type;
*         ISTATU(16) is the number of dimensions;
*         ISTATU(17) is the number of tracking points per line;
*         ISTATU(18) is the maximum length of a track;
*         ISTATU(19) is the total number of tracks;
*         ISTATU(20) is the number of tracks directions;
*         ISTATU(21) line format (by default a short 
*         format is considered but the complete format for TLM:
*         can be generated using the keyword LONG);
*         ISTATU(22) is the vectorization option;
*         ISTATU(23) is the tracking flag (-1 MC; 0 NOTR;
*               1 tracking available).
* RSTATU  real parameters for tracking:
*         RSTATU( 1) is the track length cutoff for
*               exponential functions;
*         RSTATU( 2) is the 1D line or 2D plane
*               quadrature line density;
*         RSTATU( 3) is the corner identification cutoff;
*         RSTATU( 4) is the axial quadrature line density;
*         RSTATU( 5) contains the linear track spacing
*         for general 2--D geometry and for 3--D Cartesian and
*         geometries;
*         RSTATU( 6) is the $X$ cell center;
*         RSTATU( 7) is the $y$ cell center;
*         RSTATU( 8) is the $Z$ cell center;
*         RSTATU(11) is the spatial cutoff factor for
*               tracking;
*         RSTATU(39) is the minimum volume fraction of the
*               grain in the representative volume for She-Liu-Shi 
*               model. 
* NBSLIN  maximum number of segments in a single tracking line
*         (computed by default in NXTTCG but limited to 100000
*         elements). This default value can be bypassed using
*         keyword NBSLIN.
* IQUA10  quadrature parameter for micro-structures in Bihet.
* IBIHET  type of double-heterogeneity method (=1 Sanchez-Pomraning
*         model; =2 Hebert model; =3 She-Liu-Shi model (no shadow);
*         =4 She-Liu-Shi model (with shadow)).
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*Comments:
* Input data is of the form:
*    [ EDIT iprint ]
*    [ TITLE trackt ]
*    [ NBSLIN nbslin ]
*    [ ANIS nanis ]
*    [ { RENO | REND | NORE } ]
*    [ { PISO | PSPC } ]
*    [ { PRIX | PRIY | PRIZ } denspr ]
*    [ { GAUS | CACA | CACB | LCMD | OPP1 | OGAU } npol ]
*    [ { TISO  [ { EQW | PNTN | SMS  | GAUS | LSN | QRN } ]
*                                        nangl dens [ densz ]  |
*        TSPC  [ EQW | MEDI | EQW2 ] nangl dens [ densz ] } ]
*    [ CORN cutofc ]
*    [ CUT  cutofx ]
*    [ { SYMM isymm | NOSY } ]
*    [ { NOTR | MC } ]
*    [ [ QUAB iqua10 ] [ { SAPO | HEBE | SLSI [frtm] | SLSS [frtm] } ] ]
*    with frtm minimum volume fraction of the grain in the  
*    representative volume for She-Liu-Shi model.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NSTATE
      INTEGER          IPRINT,IQUA10,IBIHET
      CHARACTER        TITLE*72
      INTEGER          ISTATU(NSTATE)
      REAL             RSTATU(NSTATE)
      INTEGER          NBSLIN
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTGET')
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
      INTEGER          IRT,IRMXR
*----
*  Initialize default values for IPRINT
*----
      IPRINT=1
      IRT=0
      IRMXR=0
      IBIHET=2
      IQUA10=5
*----
*  Get data from input file
*----
 100  CONTINUE
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
 101  CONTINUE
      IF(ITYPLU .EQ. 10) GO TO 105
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': Read error -- Character variable expected')
      IF(CARLIR(1:4) .EQ. ';') THEN
        GO TO 105
      ELSE IF(CARLIR(1:4) .EQ. 'EDIT') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- print level expected after EDIT.')
        IPRINT=INTLIR
      ELSE IF(CARLIR(1:4) .EQ. 'TITL') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': Read error -- title expected after TITL.')
        TITLE=CARLIR
      ELSE IF(CARLIR(1:4) .EQ. 'ANIS') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error --  anisotropy level expected after ANIS.')
        IF(INTLIR .LE. 0) WRITE(IOUT,9000) NAMSBR
        ISTATU(6)=MAX(ISTATU(6),INTLIR)
      ELSE IF(CARLIR(1:4) .EQ. 'RENO') THEN
        ISTATU(8)=0
      ELSE IF(CARLIR(1:4) .EQ. 'REND') THEN
        ISTATU(8)=-1
      ELSE IF(CARLIR(1:4) .EQ. 'NORE') THEN
        ISTATU(8)=1
      ELSE IF(CARLIR(1:4) .EQ. 'PISO') THEN
        ISTATU(10)=0
      ELSE IF(CARLIR(1:4) .EQ. 'PSPC') THEN
        ISTATU(10)=-1
      ELSE IF(CARLIR(1:3) .EQ. 'PRI')  THEN
         IF (CARLIR(4:4).EQ.'Z') THEN
            ISTATU(39)=3
         ELSEIF (CARLIR(4:4).EQ.'Y') THEN
            ISTATU(39)=2
         ELSEIF (CARLIR(4:4).EQ.'X') THEN
            ISTATU(39)=1
         ELSE
            CALL XABORT('NXTGET: INVALID PROJECTION AXIS FOR 3D PRISM.')
         ENDIF
         CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
         IF(ITYPLU.NE.2) THEN
            CALL XABORT('NXTGET: REAL DATA EXPECTED')
         ELSE
            RSTATU(40)=1.0/REALIR
            IF (RSTATU(40).LT.0.0)
     >        CALL XABORT('NXTGET: DELU > 0.0 EXPECTED')
         ENDIF
      ELSEIF(CARLIR(1:4) .EQ. 'GAUS') THEN
         ISTATU(13)=0
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF (ITYPLU.NE.1) GOTO 101
          ISTATU(14)=MAX(ISTATU(14),INTLIR)
      ELSEIF(CARLIR(1:4) .EQ. 'CACA') THEN
         ISTATU(13)=1
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF (ITYPLU.NE.1) GOTO 101
          ISTATU(14)=MAX(ISTATU(14),INTLIR)
      ELSEIF(CARLIR(1:4) .EQ. 'CACB') THEN
         ISTATU(13)=2
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF (ITYPLU.NE.1) GOTO 101
          ISTATU(14)=MAX(ISTATU(14),INTLIR)
      ELSEIF(CARLIR(1:4) .EQ. 'LCMD') THEN
         ISTATU(13)=3
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF (ITYPLU.NE.1) GOTO 101
          ISTATU(14)=MAX(ISTATU(14),INTLIR)
      ELSEIF(CARLIR(1:4) .EQ. 'OPP1') THEN
         ISTATU(13)=4
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF (ITYPLU.NE.1) GOTO 101
          ISTATU(14)=MAX(ISTATU(14),INTLIR)
      ELSEIF(CARLIR(1:4) .EQ. 'OGAU') THEN
         ISTATU(13)=5
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF (ITYPLU.NE.1) GOTO 101
          ISTATU(14)=MAX(ISTATU(14),INTLIR)
      ELSE IF(CARLIR(1:4) .EQ. 'TISO' .OR.
     >        CARLIR(1:4) .EQ. 'TSPC' ) THEN
        ISTATU(9)=0
        IF(CARLIR(1:4) .EQ. 'TSPC') THEN
          ISTATU(9)=1
          ISTATU(10)=-1
        ENDIF
*----
*  Azimuthal or 3-D quadrature type
*----
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .EQ. 3) THEN
          IF(CARLIR(1:4) .EQ. 'EQW') THEN
            ISTATU(15)=1
          ELSE IF(CARLIR(1:4) .EQ. 'GAUS') THEN
            ISTATU(15)=2
          ELSE IF(CARLIR(1:4) .EQ. 'MEDI') THEN
            ISTATU(15)=3
          ELSE IF(CARLIR(1:4) .EQ. 'PNTN') THEN
            ISTATU(15)=4
          ELSE IF(CARLIR(1:3) .EQ. 'SMS') THEN
            ISTATU(15)=5
          ELSE IF(CARLIR(1:3) .EQ. 'LSN') THEN
            ISTATU(15)=6
          ELSE IF(CARLIR(1:3) .EQ. 'QRN') THEN
            ISTATU(15)=7
          ELSE IF(CARLIR(1:4) .EQ. 'EQW2') THEN
            ISTATU(15)=8
          ELSE
            CALL XABORT(NAMSBR//':'//CARLIR(1:4)//
     >      ' is an invalid azimuthal or 3D quadrature type')
          ENDIF
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        ENDIF
        IF(ITYPLU .EQ. 1) THEN
          IF(INTLIR .LE. 0) WRITE(IOUT,9001) NAMSBR
          ISTATU(11)=MAX(ISTATU(11),INTLIR)
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        ENDIF
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Real density number expected')
        RSTATU(2)=REALIR
        IF(REALIR .LE. 0.0) THEN
          WRITE(IOUT,9010) NAMSBR
          RSTATU(2)=1.0
        ENDIF
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .EQ. 2) THEN
          RSTATU(4)=REALIR
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(REALIR .LE. 0.0) THEN
            WRITE(IOUT,9011) NAMSBR
            RSTATU(4)=1.0
          ENDIF
        ENDIF
        GO TO 101
      ELSE IF(CARLIR(1:4) .EQ. 'CORN') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Real value expected for CORN')
        RSTATU(3)=REALIR
        IF(REALIR .LT. 0.0) THEN
          WRITE(IOUT,9012) NAMSBR
          RSTATU(3)=0.0
        ENDIF
      ELSE IF(CARLIR(1:4) .EQ. 'CUT') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Real value expected for CUT')
        RSTATU(1)=REALIR
        IF(REALIR .LT. 0.0) THEN
          WRITE(IOUT,9013) NAMSBR
          RSTATU(1)=0.0
        ENDIF
      ELSE IF(CARLIR(1:4) .EQ. 'SYMM') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Integer value expected for SYMM')
        ISTATU(12)=INTLIR
      ELSE IF(CARLIR(1:4) .EQ. 'NOSY') THEN
        ISTATU(12)=0
      ELSE IF(CARLIR(1:4) .EQ. 'NOTR') THEN
        ISTATU(23)=0
      ELSE IF(CARLIR(1:2) .EQ. 'MC') THEN
        ISTATU(23)=-1
      ELSE IF(CARLIR(1:4) .EQ. 'TRAK') THEN
        IRT=1
      ELSE IF(CARLIR(1:4) .EQ. 'MAXR') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Integer value expected for MAXR')
        IRMXR=MAX(INTLIR,1)
      ELSE IF(CARLIR .EQ. 'NBSLIN') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- nbslin value expected.')
        NBSLIN=MAX(INTLIR,NBSLIN)
      ELSE IF(CARLIR(1:4) .EQ. 'SCFT') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Real value expected for SCFT')
        RSTATU(11)=REALIR
        IF(REALIR .LT. 0.0) THEN
          WRITE(IOUT,9012) NAMSBR
          RSTATU(11)=1.0
        ENDIF
      ELSE IF(CARLIR(1:4) .EQ. 'ONEG') THEN
        ISTATU(22)=0
      ELSE IF(CARLIR(1:4) .EQ. 'ALLG') THEN
        ISTATU(22)=1
      ELSE IF(CARLIR(1:4) .EQ. 'XCLL') THEN
        ISTATU(22)=2
      ELSE IF(CARLIR(1:4) .EQ. 'QUAB') THEN
        CALL REDGET(ITYPLU,IQUA10,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Integer value expected for QUAB')
      ELSE IF(CARLIR .EQ. 'LONG') THEN
        ISTATU(21)=1
      ELSE IF(CARLIR(1:4) .EQ. 'SAPO') THEN
        IBIHET=1
      ELSE IF(CARLIR(1:4) .EQ. 'HEBE') THEN
        IBIHET=2
      ELSE IF(CARLIR(1:4) .EQ. 'SLSI') THEN
        IBIHET=3
        RSTATU(39)=0.05
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF (ITYPLU.NE.2) GOTO 101
        RSTATU(39)=REALIR
      ELSE IF(CARLIR(1:4) .EQ. 'SLSS') THEN
        IBIHET=4
        RSTATU(39)=0.05
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF (ITYPLU.NE.2) GOTO 101
        RSTATU(39)=REALIR
      ELSE
        CALL XABORT(NAMSBR//': Keyword '//CARLIR(1:4)//' is invalid.')
      ENDIF
      GO TO 100
 105  CONTINUE
      IF( ISTATU(9) .EQ. 0) THEN
        ISTATU(13)=0
      ENDIF
*----
*  Processing finished, return
*----
      IF(IRT .GT. 0) WRITE(IOUT,9020) NAMSBR
      IF(IRMXR .GT. 0) WRITE(IOUT,9021) NAMSBR
      RETURN
*----
*  Warning formats
*----
 9000 FORMAT(1X,'Warning from ',A6,2X,'Invalid anisotropy level'/
     >1X,'Use default value : nanis=1')
 9001 FORMAT(1X,'Warning from ',A6,2X,'Invalid number of angles'/
     >1X,'Use default value : nangle=1')
 9010 FORMAT(1X,'Warning from ',A6,2X,'Invalid tracking density'/
     >1X,'Use default value : dens=1.0')
 9011 FORMAT(1X,'Warning from ',A6,2X,'Invalid axial tracking density'/
     >1X,'Use default value : densz=1.0')
 9012 FORMAT(1X,'Warning from ',A6,2X,'Invalid corner proximity'/
     >1X,'Use default value : pcorn=0.0')
 9013 FORMAT(1X,'Warning from ',A6,2X,'Invalid exponential cutoff'/
     >1X,'Use default value : cutofx=0.0')
 9020 FORMAT(1X,'Warning from ',A6,1X,'-- Keyword TRAK not used ',
     >'by module NXT: but kept for compatibility with module EXCELT:')
 9021 FORMAT(1X,'Warning from ',A6,1X,'-- Keyword MAXR not used ',
     >'by module NXT: but kept for compatibility with module EXCELT:')
      END
