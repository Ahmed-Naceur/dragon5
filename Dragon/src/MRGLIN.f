*DECK MRGLIN
      SUBROUTINE MRGLIN(IPRINT,IFTRKO,NSOUTO,NVOUTO,IFTRKN,
     >                  IMERGE,MXSEG)
*
*----------
*
*Purpose:
* Merge volume surface information on track file.
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
* IFTRKO  old tracking file.
* NSOUTO  old number of surfaces.
* NVOUTO  old number of regions.
* IFTRKN  new tracking file.
* IMERGE  merged position.
* MXSEG   maximum number of segments.
*
*----------
*
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='MRGLIN')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          IPRINT,IFTRKO,NSOUTO,NVOUTO,IFTRKN,MXSEG
      INTEGER          IMERGE(-NSOUTO:NVOUTO)
*----
*  LOCAL VARIABLES
*----
      INTEGER          ITRAK,IANG,NLINEO,NLINEN,ILINE,
     >                 ISEG,IVSO
      REAL             WEIGHT
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NRSEG
      REAL, ALLOCATABLE, DIMENSION(:) :: PATH
*----
*  LOOP OVER TRACKS
*----
      ALLOCATE(NRSEG(MXSEG),PATH(MXSEG))
      ITRAK=0
 1000 CONTINUE
        READ (IFTRKO,END=1010) IANG,NLINEO,WEIGHT,
     >                         (NRSEG(ILINE),ILINE=1,NLINEO),
     >                         (PATH(ILINE),ILINE=1,NLINEO)
*----
*  SCAN NRSEG AND RESET TO NEW VOLUME AND SURFACE NUMBER
*----
        ITRAK=ITRAK+1
        IF(IPRINT.GE.1000) THEN
          WRITE(IOUT,6000) ITRAK,IANG,NLINEO,WEIGHT
          WRITE(IOUT,6010)
     >      (NRSEG(ILINE),PATH(ILINE),ILINE=1,NLINEO)
        ENDIF
        DO 100 ILINE=1,NLINEO
          DO 110 IVSO=-NSOUTO,NVOUTO
            IF(NRSEG(ILINE) .EQ. IVSO ) THEN
              NRSEG(ILINE) = IMERGE(IVSO)
              GO TO 115
            ENDIF
 110      CONTINUE
 115      CONTINUE
 100    CONTINUE
*----
*  COMPRESS REGION OF SUCCESSIVE IDENTICAL REGION
*  EXCEPT FOR SURFACES
*----
        NLINEN=1
        ISEG=NRSEG(NLINEN)
        DO 120 ILINE=2,NLINEO
          IF(NRSEG(ILINE) .EQ. ISEG .AND.
     >       ISEG .GT. 0                  ) THEN
            PATH(NLINEN)=PATH(NLINEN)+PATH(ILINE)
          ELSE
            NLINEN=NLINEN+1
            NRSEG(NLINEN)=NRSEG(ILINE)
            PATH(NLINEN)=PATH(ILINE)
            ISEG=NRSEG(NLINEN)
          ENDIF
 120    CONTINUE
        WRITE(IFTRKN) IANG,NLINEN,WEIGHT,
     >               (NRSEG(ILINE),ILINE=1,NLINEN),
     >               (PATH(ILINE),ILINE=1,NLINEN)
        IF(IPRINT.GE.100) THEN
          WRITE(IOUT,6001) ITRAK,IANG,NLINEN,WEIGHT
          WRITE(IOUT,6010)
     >      (NRSEG(ILINE),PATH(ILINE),ILINE=1,NLINEN)
        ENDIF
      GO TO 1000
 1010 CONTINUE
      DEALLOCATE(PATH,NRSEG)
*----
*  FORMAT
*----
 6000 FORMAT(' INITIAL LINE = ',I10/
     >       ' PARAMETERS = ',2I10,1P,E15.7)
 6001 FORMAT(' FINAL LINE = ',I10/
     >       ' PARAMETERS = ',2I10,1P,E15.7)
 6010 FORMAT(1P,5(I10,E15.7))
      RETURN
      END
