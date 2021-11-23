*DECK SALT
      SUBROUTINE SALT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To analyze and track a geometry data structure using the SALT
* algorithm.
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* A. Hebert and G. Marleau
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
*         =0 for a data structure in creation mode;
*         =1 for a data structure in modifications mode;
*         =2 for a data structure in read-only mode.
* KENTRY  data structure pointer.
*
*Comments:
* Instructions for the use of the SALT: module:
*   TRKFIL VOLTRK := SALT: SURFIL [ GEOMETRY ] :: (saltget) ;
*   where
*     TRKFIL   : sequential binary tracking file to be created
*     VOLTRK   : tracking data structure (signature L_TRACK)
*     SURFIL   : sequential ascii file used to store the surfacic
*                elements of the geometry.
*     GEOMETRY : optional geometry data structure used id BIHET is set
*                (signature L_GEOM)
*     (saltget): Processing options
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
      PARAMETER       (IOUT=6,NAMSBR='SALT  ')
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
      INTEGER          IMGEO,IFTRK,FGEO
      INTEGER          IGTRK
      INTEGER          IEN,ITC
      INTEGER          IQUA10,IBIHET
      CHARACTER        HSIGN*12
      INTEGER          ISTATT(NSTATE)
      REAL             RSTATT(NSTATE)
      CHARACTER        TITLE*72
      INTEGER          IPRINT
      INTEGER          NBSLIN
      INTEGER          ILONG,ITYLCM
      DOUBLE PRECISION RCUTOF
*----
*  Validate entry parameters
*----
      IF((NENTRY.LT.3).OR.(NENTRY.GT.4)) CALL XABORT(NAMSBR//
     >  ': Three or four data structures permitted')
      IPGEO=C_NULL_PTR
      FGEO=0
      IMGEO=0
      NBSLIN=100000
*----
*  Scan data structure to determine type and mode
*----
      DO IEN=1,2
        IF(JENTRY(IEN).NE.0) CALL XABORT(NAMSBR//
     >  ': Object in creation mode expected')
        IF((IENTRY(IEN).EQ.1).OR.(IENTRY(IEN).EQ.2)) THEN
          IPTRK=KENTRY(IEN)
          HSIGN='L_TRACK     '
          CALL LCMPTC(IPTRK,'SIGNATURE',12,1,HSIGN)
          HSIGN='EXCELL'
          CALL LCMPTC(IPTRK,'TRACK-TYPE',12,1,HSIGN)
        ELSE IF(IENTRY(IEN).EQ.3) THEN
          IFTRK=FILUNIT(KENTRY(IEN))
        ENDIF
      ENDDO
      DO IEN=3,NENTRY
        IF(JENTRY(IEN).NE.2) CALL XABORT(NAMSBR//
     >  ': Object in read-only mode expected')
        IF((IENTRY(IEN).EQ.1).OR.(IENTRY(IEN).EQ.2)) THEN
          CALL LCMGTC(IPTRK,'SIGNATURE',12,1,HSIGN)
          IF(HSIGN.NE.'L_GEOM') CALL XABORT(NAMSBR//
     >    ': L_GEOM signature expected for '//HENTRY(IEN))
          IPGEO=KENTRY(IEN)
          IMGEO=-1
        ELSE IF(IENTRY(IEN).EQ.4) THEN
          FGEO=FILUNIT(KENTRY(IEN))
        ENDIF
      ENDDO
      IF(FGEO.EQ.0) CALL XABORT(NAMSBR//
     >  ': The surfacic file is not defined')
*----
*  Initialize tracking parameters to 0
*----
      CALL XDISET(ISTATT,NSTATE,0)
      CALL XDRSET(RSTATT,NSTATE,0.0)
*----
*  Define default tracking options that are different from 0
*----
      ISTATT(6)=1
      ISTATT(7)=4
      ISTATT(11)=1
      ISTATT(12)=-1
      ISTATT(13)=1
      ISTATT(15)=1
      ISTATT(16)=2
      ISTATT(22)=0
      ISTATT(23)=1
      IF(IMGEO .EQ. -1) THEN
        CALL LCMLEN(IPGEO,'BIHET',ILONG,ITYLCM)
        IF(ILONG.NE.0) ISTATT(40)=1
      ENDIF
      RSTATT(11)=1.0
      TITLE=' '
*----
*  Recover processing option
*----
      CALL NXTGET(NSTATE,IPRINT,TITLE ,ISTATT,RSTATT,NBSLIN,IQUA10,
     >            IBIHET)
      IF((ISTATT(9).EQ.1).AND.(ISTATT(15).EQ.1)) THEN
        ISTATT(15)=8 ! replace EQW by EQW2
      ENDIF
*----
*  Compute Bickley functions
*----
      IF(ISTATT(23).GT.0) CALL XDRTA1(IPTRK,.TRUE.,.TRUE.)
*----
*  Save updated STATE-VECTOR, TITLE and EXCELL track options
*  on tracking data structure
*----
      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,ISTATT)
      CALL LCMPUT(IPTRK,'EXCELTRACKOP',NSTATE,2,RSTATT)
      CALL LCMPTC(IPTRK,'TITLE',72,1,TITLE)
*----
*  Analyse geometry if required
*----
      RCUTOF=DBLE(RSTATT(3))
      CALL SALACG(FGEO ,IPTRK, NBSLIN, RCUTOF, IPRINT)
*----
*  Track geometry if required
*----
      IF(ISTATT(9) .GE. 0 .AND. ISTATT(23) .EQ. 1) THEN
        IGTRK=1
        CALL SALTCG(IPTRK,IFTRK,IPRINT,IGTRK,NBSLIN)
      ENDIF
*----
*  Release allocated memory in SALT module
*----
      CALL SALEND()
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
*  Formats
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
