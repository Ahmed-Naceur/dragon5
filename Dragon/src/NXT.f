*DECK NXT
      SUBROUTINE NXT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Module used to analyze and track a geometry data structure based 
* on the new EXCELL type procedure.
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
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*Comments:
* Instructions for the use of the NXT: module:
*   Option 1 -- Analyze and optionnally track a basic geometry
*   [ TRKFIL ] VOLTRK     := NXT: GEOMETRY :: (nxtget) ;
*   Option 2 -- Track a geometry already analyzed
*   TRKFIL     VOLTRK     := NXT: VOLTRK   :: (nxtget) ;
*   where
*     TRKFIL   : sequential binary tracking file to be created
*     VOLTRK   : tracking data structure
*                (signature L_TRACK)
*     GEOMETRY : geometry data structure
*                (signature L_GEOM)
*     (nxtget) : Processing options
*                (read from input using the NXTGET routine).
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
      PARAMETER       (IOUT=6,NAMSBR='NXT   ')
      INTEGER          ILCMUP,ILCMDN
      PARAMETER       (ILCMUP=1,ILCMDN=2)
      INTEGER          NSTATE,MAXENT
      PARAMETER       (NSTATE=40,MAXENT=2)
      INTEGER          IUTYPE
      PARAMETER       (IUTYPE=2)
*----
*  Local variables
*----
      TYPE(C_PTR)      IPGEO,IPTRK
      INTEGER          IMGEO,IMTRK,IFTRK,IMFTRK
      INTEGER          IGANA,IGTRK
      INTEGER          IEN,ITC
      INTEGER          IQUA10,IBIHET
      CHARACTER        HSIGN*12
      CHARACTER        TEXT12*12
      INTEGER          ISTATT(NSTATE)
      REAL             RSTATT(NSTATE)
      CHARACTER        TITLE*72
      INTEGER          IPRINT,ITITL(18)
      INTEGER          NBSLIN
      INTEGER          ILONG,ITYLCM
*----
*  Validate entry parameters
*----
      IF(NENTRY .LT. 2) CALL XABORT(NAMSBR//
     >  ': At least two data structures required')
      IF(NENTRY .GT. 3) CALL XABORT(NAMSBR//
     >  ': Maximum of three data structures permitted')
      IPGEO=C_NULL_PTR
      IMGEO=0
      IPTRK=C_NULL_PTR
      IMTRK=0
      IFTRK=0
      IMFTRK=0
      IGANA=0
      IGTRK=0
      NBSLIN=100000
*----
*  Scan data structure to determine type and mode
*----
      DO IEN=1,NENTRY
        IF(IENTRY(IEN) .EQ. 1 .OR. IENTRY(IEN) .EQ. 2) THEN
          IF(JENTRY(IEN) .EQ. 0) THEN
            IPTRK=KENTRY(IEN)
            IMTRK=2
            HSIGN='L_TRACK     '
            CALL LCMPTC(IPTRK,'SIGNATURE',12,1,HSIGN)
            HSIGN='EXCELL'
            CALL LCMPTC(IPTRK,'TRACK-TYPE',12,1,HSIGN)
          ELSE
            CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
            IF(HSIGN .EQ. 'L_GEOM') THEN
              IPGEO=KENTRY(IEN)
              IF(JENTRY(IEN) .NE. 2) CALL XABORT(NAMSBR//
     >        ': Geometry data structure not in read-only mode')
              TEXT12=HENTRY(IEN)
              CALL LCMPTC(IPTRK,'LINK.GEOM',12,1,TEXT12)
              IMGEO=-1
            ELSE IF(HSIGN .EQ. 'L_TRACK') THEN
              IPTRK=KENTRY(IEN)
              IF(JENTRY(IEN) .NE. 1) CALL XABORT(NAMSBR//
     >        ': Tracking data structure not in update mode')
              IMTRK=1
              CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,HSIGN)
              IF(HSIGN .NE. 'EXCELL') CALL XABORT(NAMSBR//
     >        ': Tracking data structure type is invalid')
            ELSE
              CALL XABORT(NAMSBR//
     >        ': Invalid signature for '//HENTRY(IEN))
            ENDIF
          ENDIF
        ELSE IF(IENTRY(IEN) .EQ. 3) THEN
          IF(JENTRY(IEN) .NE. 0) CALL XABORT(NAMSBR//
     >        ': Geometry data structure not in creation mode')
          IFTRK=FILUNIT(KENTRY(IEN))
          IMFTRK=2
        ELSE
          CALL XABORT(NAMSBR//
     >    ': Invalid data structure format for '//HENTRY(IEN))
        ENDIF
      ENDDO
*----
*  Select processing option from data structures provided
*----
      IF(IMGEO .EQ. -1) THEN
        IF(IMTRK .NE. 2) CALL XABORT(NAMSBR//
     >  ': Creation mode tracking data structure required')
        IGANA=1
        IF(IMFTRK .EQ. 2) IGTRK=1
      ELSE IF(IMTRK .EQ. 1) THEN
        IF(IMFTRK .NE. 2) CALL XABORT(NAMSBR//
     >  ': Creation mode tracking file required')
        IGTRK=1
      ELSE
        CALL XABORT(NAMSBR//': No processing option identified')
      ENDIF
*----
*  Initialize tracking parameters to 0
*----
      CALL XDISET(ISTATT,NSTATE,0)
      CALL XDRSET(RSTATT,NSTATE,0.0)
*----
*  Read state vectors available
*----
      IF(IMTRK .EQ.  1) THEN
        CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATT)
        CALL LCMGET(IPTRK,'EXCELTRACKOP',RSTATT)
        CALL LCMGET(IPTRK,'TITLE',ITITL)
        WRITE(TITLE,'(18A4)') (ITITL(ITC),ITC=1,18)
        IF(ISTATT(7) .NE. 4 ) CALL XABORT(NAMSBR//
     >  ': Tracking data structure incompatible with current module')
        ISTATT(23)=1
      ELSE
*----
*  Define default tracking options that are different from 0
*----
        ISTATT(6)=1
        ISTATT(7)=4
        ISTATT(11)=1
        ISTATT(12)=-1
        ISTATT(13)=1
        ISTATT(15)=1
        ISTATT(22)=0
        ISTATT(23)=1
        IF(IMFTRK .EQ. 0) ISTATT(22)=3
        IF(IMTRK .EQ. 2 .AND. IMGEO .EQ. -1) THEN
          CALL LCMLEN(IPGEO,'BIHET',ILONG,ITYLCM)
          IF(ILONG.NE.0) ISTATT(40)=1
        ENDIF
        RSTATT(11)=1.0
        TITLE=' '
        HSIGN='EXCELL'
      ENDIF
*----
*  Recover processing option
*----
      CALL NXTGET(NSTATE,IPRINT,TITLE ,ISTATT,RSTATT,NBSLIN,IQUA10,
     >            IBIHET)
*----
*  Compute Bickley functions
*----
      IF((ISTATT(23).GT.0).AND.(IMTRK.EQ.2)) THEN
        CALL XDRTA1(IPTRK,.TRUE.,.TRUE.)
      ENDIF
*----
*  Save updated STATE-VECTOR, TITLE and EXCELL track options
*  on tracking data structure
*----
      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,ISTATT)
      CALL LCMPUT(IPTRK,'EXCELTRACKOP',NSTATE,2,RSTATT)
      READ(TITLE,'(18A4)') (ITITL(ITC),ITC=1,18)
      CALL LCMPUT(IPTRK,'TITLE',18,3,ITITL)
*----
*  Analyse geometry if required
*----
      IF(IGANA .EQ. 1) THEN
        CALL NXTACG(IPGEO ,IPTRK ,IPRINT)
      ENDIF
*----
*  If a prismatic 3D tracking is requested,
*  create 2D projected geometry analysis 
*----
      IF(ISTATT(39) .NE. 0) THEN
         CALL NXTPR3(IPTRK)     
      ENDIF
*----
*  Track geometry if required
*----
      IF(ISTATT(9) .GE. 0 .AND. ISTATT(23) .EQ. 1) THEN
        IF(ISTATT(39) .NE. 0) CALL LCMSIX(IPTRK,'PROJECTION',1)
        CALL NXTTCG(IPTRK ,IFTRK ,IPRINT,IGTRK ,NBSLIN)
        IF(ISTATT(39) .NE. 0) CALL LCMSIX(IPTRK,' ',2)
      ENDIF
*----
*  Add useful information for the Monte-Carlo method
*----
      IF(ISTATT(23) .EQ. -1) THEN
         CALL NXTMCA(IPTRK)
      ENDIF
*----
*  Process double heterogeneity (BIHET) data (if available)
*----
      IF(ISTATT(40) .NE. 0) THEN
         CALL XDRTBH(IPGEO,IPTRK,IQUA10,IBIHET,IPRINT,RSTATT(39))
      ENDIF
*----
*  Processing finished, return
*----
      IF(IPRINT .GT. 1) THEN
         CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATT)
         WRITE(IOUT,100) (ISTATT(ITC),ITC=1,10)
         WRITE(IOUT,120) (ISTATT(ITC),ITC=11,22)
         WRITE(IOUT,130) ISTATT(23),ISTATT(40)
      ENDIF
      RETURN
*----
*  Warning formats
*----
  100 FORMAT(/
     1 14H STATE VECTOR:/
     2 7H NREG  ,I9,22H   (NUMBER OF REGIONS)/
     3 7H KPN   ,I9,23H   (NUMBER OF UNKNOWNS)/
     4 7H ILK   ,I9,39H   (0=LEAKAGE PRESENT/1=LEAKAGE ABSENT)/
     5 7H NBMIX ,I9,36H   (MAXIMUM NUMBER OF MIXTURES USED)/
     6 7H NSURF ,I9,29H   (NUMBER OF OUTER SURFACES)/
     7 7H NANI  ,I9,48H   (1=P0 CROSS SECTIONS/2=P1 CROSS SECTIONS/...)/
     8 7H GEOT  ,I9,21H   (TYPE OF GEOMETRY)/
     9 7H NORM  ,I9,48H   (NORMALIZATION OPTION 1=ABSENT/0=GLOBAL/-1=NO,
     1 21HRMALIZATION BY ANGLE)/
     2 7H TRKT  ,I9,36H   (TRACKING TYPE 0=FINITE/1=CYCLIC)/
     3 7H BOUND ,I9,48H   (BOUNDARY CONDITIONS TYPE 0=ISOTROPIC/1=SPECU,
     4 4HLAR))
  120 FORMAT(
     1 7H NANG  ,I9,30H   (NUMBER OF TRACKING ANGLES)/
     2 7H ASYM  ,I9,28H   (ANGULAR SYMMETRY FACTOR)/
     3 7H POLQUA,I9,32H   (POLAR ANGLE QUADRATURE TYPE)/
     4 7H POLOAQ,I9,33H   (POLAR ANGLE QUADRATURE ORDER)/
     5 7H AZMQUA,I9,47H   (AZIMUTHAL OR SOLID ANGULAR QUADRATURE TYPE)/
     6 7H NDIM  ,I9,25H   (NUMBER OF DIMENSIONS)/
     7 7H NPOINT,I9,40H   (NUMBER OF TRACKING POINTS ON A LINE)/
     8 7H MAXSGL,I9,30H   (MAXIMUM LENGTH OF A TRACK)/
     9 7H NTLINE,I9,37H   (TOTAL NUMBER OF TRACKS GENERATED)/
     1 7H NBTDIR,I9,47H   (TOTAL NUMBER OF TRACK DIRECTIONS PROCESSED)/
     2 7H NANGL ,I9,47H   (NUMBER OF TRACK DIRECTION ANGLES CONSIDERED,
     3 20H IN THE INTEGRATION)/
     4 7H INSB  ,I9,25H   (VECTORIZATION OPTION))
  130 FORMAT(
     1 7H ITRACK,I9,47H   (-1=MONTE-CARLO/0=DESACTIVATES TRACKING FILE,
     2 39H BUILD/1=ACTIVATES TRACKING FILE BUILD)/
     3 7H IBIHET,I9,46H   (0/1=DOUBLE HETEROGENEITY IS NOT/IS ACTIVE))
      END
