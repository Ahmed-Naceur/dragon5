*DECK MRGXTC
      SUBROUTINE MRGXTC(IFTRKO,IFTRKN,IFTRKE,IPRINT,IUPD,NDIM,
     >                  NALBG ,NANGL ,NSOUTO,NVOUTO,MXSEG,IMERGE)
*
*----------
*
*Purpose:
* Subdivide tracking file into 2 sets.
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
* IFTRKO  old tracking file.
* IFTRKN  old part r (residual) tracking file.
* IFTRKE  old part e (extracted) tracking file.
* IPRINT  print level.
* IUPD    type of merge required:
*         IUPD(1) for region merge;
*         IUPD(2) for surface merge;
*         IUPD(3) for material merge;
*         IUPD(4) for albedo merge.
* NDIM    number of dimensions.
* NALBG   number of albedos.
* NANGL   number of track directions.
* NSOUTO  old number of surfaces.
* NVOUTO  old number of regions.
* MXSEG   maximum number of segments.
* IMERGE  merged position.
*
*----------
*
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='MRGXTC')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          IFTRKO,IFTRKN,IFTRKE
      INTEGER          IPRINT,IUPD(4),NDIM,
     >                 NALBG,NANGL,NSOUTO,NVOUTO,MXSEG
      INTEGER          IMERGE(-NSOUTO:NVOUTO)
*----
*  LOCAL VARIABLES
*----
      INTEGER          IVSO,ITRAK,IANG,ILINE,NLINE,IEXT
      REAL             WEIGHT
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MATO,ICODE,NRSEG
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLO,ALBD
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ANGLES,DENSTY
      REAL, ALLOCATABLE, DIMENSION(:) :: PATH
*----
*  Processing starts:
*  print routine openning output header if required
*----
      ALLOCATE(MATO(-NSOUTO:NVOUTO),VOLO(-NSOUTO:NVOUTO))
      IF(IPRINT.GE.10) THEN
        WRITE(IOUT,6000)
      ENDIF
*----
*  Read old Volume/surface and albedo/material information
*  and print if required
*----
      READ(IFTRKO) (VOLO(IVSO),IVSO=-NSOUTO,NVOUTO)
      READ(IFTRKO) (MATO(IVSO),IVSO=-NSOUTO,NVOUTO)
      IF(IPRINT.GE.10) THEN
        WRITE(IOUT,6010)
        WRITE(IOUT,6011) (MATO(IVSO),VOLO(IVSO),IVSO=-NSOUTO,NVOUTO)
      ENDIF
*----
*  Save new Volume/surface and albedo/material information
*----
      WRITE(IFTRKN) (VOLO(IVSO),IVSO=-NSOUTO,NVOUTO)
      WRITE(IFTRKE) (VOLO(IVSO),IVSO=-NSOUTO,NVOUTO)
      WRITE(IFTRKN) (MATO(IVSO),IVSO=-NSOUTO,NVOUTO)
      WRITE(IFTRKE) (MATO(IVSO),IVSO=-NSOUTO,NVOUTO)
*----
*  Read and save BC and tracking directions and density
*----
      ALLOCATE(ICODE(NALBG),ALBD(NALBG),ANGLES(NDIM*NANGL),
     >         DENSTY(NANGL))
      READ(IFTRKO) (ICODE(IVSO),IVSO=1,NALBG)
      READ(IFTRKO) (ALBD(IVSO),IVSO=1,NALBG)
      READ(IFTRKO) (ANGLES(IVSO),IVSO=1,NDIM*NANGL)
      READ(IFTRKO) (DENSTY(IVSO),IVSO=1,NANGL)
      WRITE(IFTRKN) (ICODE(IVSO),IVSO=1,NALBG)
      WRITE(IFTRKN) (ALBD(IVSO),IVSO=1,NALBG)
      WRITE(IFTRKN) (ANGLES(IVSO),IVSO=1,NDIM*NANGL)
      WRITE(IFTRKN) (DENSTY(IVSO),IVSO=1,NANGL)
      WRITE(IFTRKE) (ICODE(IVSO),IVSO=1,NALBG)
      WRITE(IFTRKE) (ALBD(IVSO),IVSO=1,NALBG)
      WRITE(IFTRKE) (ANGLES(IVSO),IVSO=1,NDIM*NANGL)
      WRITE(IFTRKE) (DENSTY(IVSO),IVSO=1,NANGL)
      DEALLOCATE(DENSTY,ANGLES,ALBD,ICODE)
*----
*  select track for tracking files
      ALLOCATE(NRSEG(MXSEG),PATH(MXSEG))
      ITRAK=0
 1000 CONTINUE
        READ (IFTRKO,END=1010) IANG,NLINE,WEIGHT,
     >                         (NRSEG(ILINE),ILINE=1,NLINE),
     >                         (PATH(ILINE),ILINE=1,NLINE)
C----
C  SCAN NRSEG AND RESET TO NEW VOLUME AND SURFACE NUMBER
C----
        ITRAK=ITRAK+1
        IEXT=1
        DO ILINE=1,NLINE
          DO IVSO=1,-IUPD(1)
            IF(NRSEG(ILINE) .EQ. IMERGE(IVSO)) GO TO 1005
          ENDDO
        ENDDO
        IEXT=0
 1005   CONTINUE
        IF(IEXT .EQ. 1) THEN
          WRITE(IFTRKE) IANG,NLINE,WEIGHT,
     >                 (NRSEG(ILINE),ILINE=1,NLINE),
     >                 (PATH(ILINE),ILINE=1,NLINE)
          IF(IPRINT.GE.1000) THEN
            WRITE(IOUT,6020) ITRAK,IANG,NLINE,WEIGHT
            WRITE(IOUT,6025)
     >       (NRSEG(ILINE),PATH(ILINE),ILINE=1,NLINE)
          ENDIF
        ELSE
          WRITE(IFTRKN) IANG,NLINE,WEIGHT,
     >                 (NRSEG(ILINE),ILINE=1,NLINE),
     >                 (PATH(ILINE),ILINE=1,NLINE)
          IF(IPRINT.GE.1000) THEN
            WRITE(IOUT,6021) ITRAK,IANG,NLINE,WEIGHT
            WRITE(IOUT,6025)
     >       (NRSEG(ILINE),PATH(ILINE),ILINE=1,NLINE)
          ENDIF
        ENDIF
        GO TO 1000
 1010 CONTINUE
      DEALLOCATE(PATH,NRSEG)
      DEALLOCATE(VOLO,MATO)
*----
*  Print output header if required
*  and return
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  PRINT FORMATS
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(' Region description '/
     >       4(' Region  ->               Volume'))
 6011 FORMAT(4(1X,I6,5X,F20.8))
 6020 FORMAT(' Line = ',I10,' is extracted '/
     >       ' Parameter = ',2I10,1P,E15.7)
 6021 FORMAT(' Line = ',I10,' is kept '/
     >       ' Parameter = ',2I10,1P,E15.7)
 6025 FORMAT(1P,5(I10,E15.7))
      END
