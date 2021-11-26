*DECK MRGVST
      SUBROUTINE MRGVST(IFTRKO,IFTRKN,IPRINT,IUPD  ,NDIM  ,NALBG,NANGL,
     >                  NSOUTO,NVOUTO,NSOUTN,NVOUTN,IMERGE,MIXN)
*
*----------
*
*Purpose:
* Merge volume and surface on track file and save.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
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
* IFTRKN  old tracking file.
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
* NSOUTN  new number of surfaces.
* NVOUTN  new number of regions.
* IMERGE  merging index.
* MIXN    new albedos and material for old regions.
*
*----------
*
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='MRGVST')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          IFTRKO,IFTRKN
      INTEGER          IPRINT,IUPD(4),NDIM,NALBG,NANGL,
     >                 NSOUTO,NVOUTO,NSOUTN,NVOUTN
      INTEGER          IMERGE(-NSOUTO:NVOUTO),MIXN(NVOUTO)
*----
*  LOCAL VARIABLES
*----
      INTEGER          IVSN,IVSO,IVST,ITC
      DOUBLE PRECISION DVOL
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MATO,MATN,ICODE
      REAL, ALLOCATABLE, DIMENSION(:)  :: VOLO,VOLN
      REAL, ALLOCATABLE, DIMENSION(:) :: ALBD
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ANGLES,DENSTY
*----
*  Processing starts:
*  print routine openning output header if required
*----
      IF(IPRINT.GE.10) THEN
        WRITE(IOUT,6000)
      ENDIF
*----
*  Read old Volume/surface and albedo/material information
*  and print if required
*----
      ALLOCATE(MATO(-NSOUTO:NVOUTO),VOLO(-NSOUTO:NVOUTO))
      ALLOCATE(MATN(-NSOUTN:NVOUTN),VOLN(-NSOUTN:NVOUTN))
      READ(IFTRKO) (VOLO(ITC),ITC=-NSOUTO,NVOUTO)
      READ(IFTRKO) (MATO(ITC),ITC=-NSOUTO,NVOUTO)
      IF(IPRINT.GE.10) THEN
        WRITE(IOUT,6010)
        WRITE(IOUT,6020) (MATO(IVSO),VOLO(IVSO),IVSO=-1,-NSOUTO,-1)
      ENDIF
*----
*  VERIFY IF BOUNDARY CONDITIONS ARE ADEQUATE FOR MERGE VECTOR
*----
      IF(IUPD(2). LT. 0) THEN
        DO 100 IVSN=-NSOUTN,-1
          MATN(IVSN)=0
          DVOL=0.0D0
          DO 101 IVSO=-NSOUTO,-1
           IF(IMERGE(IVSO) .EQ. IVSN) THEN
             IF(MATN(IVSN) .EQ. 0) THEN
               MATN(IVSN)=MATO(IVSO)
             ELSE IF(MATN(IVSN) .NE. MATO(IVSO))THEN
               WRITE(IOUT,6100) NAMSBR,IVSN,MATN(IVSN),IVSO,MATO(IVSO)
               WRITE(IOUT,6101)
               WRITE(IOUT,6102) (IVST,MATO(IVST),IVST=-NSOUTO,-1)
               CALL XABORT(NAMSBR//
     >           ': BOUNDARY CONDITIONS INCOMPATIBLE FOR MERGE')
             ENDIF
             DVOL=DVOL+DBLE(VOLO(IVSO))
           ENDIF
 101      CONTINUE
          VOLN(IVSN)=REAL(DVOL)
 100    CONTINUE
      ELSE
        DO 110 IVSO=-NSOUTO,-1
          MATN(IVSO)=MATO(IVSO)
          VOLN(IVSO)=VOLO(IVSO)
 110    CONTINUE
      ENDIF
      IF(IPRINT.GE.10) THEN
        WRITE(IOUT,6011)
        WRITE(IOUT,6020) (MATN(IVSN),VOLN(IVSN),IVSN=-1,-NSOUTN,-1)
        WRITE(IOUT,6012)
        WRITE(IOUT,6020) (MATO(IVSO),VOLO(IVSO),IVSO=1,NVOUTO)
      ENDIF
*----
*  CHANGE ORIGINAL MATERIAL IF REQUESTED
*----
      IF(IUPD(3) .GT. 0) THEN
        DO 120 IVSO=1,IUPD(3)
          MATO(IVSO)=MIXN(IVSO)
 120    CONTINUE
      ENDIF
*----
*  VERIFY IF MATERIALS ARE ADEQUATE FOR MERGE VECTOR
*----
      IF(IUPD(1) .GT. 0) THEN
        MATN(0)=0
        VOLN(0)=0.0
        DO 130 IVSN=1,NVOUTN
          MATN(IVSN)=0
          DVOL=0.0D0
          DO 131 IVSO=1,NVOUTO
           IF(IMERGE(IVSO) .EQ. IVSN) THEN
             IF(MATN(IVSN) .EQ. 0) THEN
               MATN(IVSN)=MATO(IVSO)
             ELSE IF(MATN(IVSN) .NE. MATO(IVSO))THEN
               WRITE(IOUT,6200) NAMSBR,IVSN,MATN(IVSN),IVSO,MATO(IVSO)
               WRITE(IOUT,6201)
               WRITE(IOUT,6202) (IVST,MATO(IVST),IVST=1,NVOUTO)
               CALL XABORT(NAMSBR//
     >           ': MATERIALS INCOMPATIBLE FOR MERGE')
             ENDIF
             DVOL=DVOL+DBLE(VOLO(IVSO))
           ENDIF
 131      CONTINUE
          VOLN(IVSN)=REAL(DVOL)
 130    CONTINUE
      ELSE
        DO 140 IVSO=1,NVOUTO
          MATN(IVSO)=MATO(IVSO)
          VOLN(IVSO)=VOLO(IVSO)
 140    CONTINUE
      ENDIF
*----
*  CHANGE FINAL MATERIAL IF REQUESTED
*----
      IF(IUPD(3) .LT. 0) THEN
        DO 150 IVSN=1,-IUPD(3)
          MATN(IVSN)=MIXN(IVSN)
 150    CONTINUE
      ENDIF
      IF(IPRINT.GE.10) THEN
        WRITE(IOUT,6013)
        WRITE(IOUT,6020) (MATN(IVSN),VOLN(IVSN),IVSN=1,NVOUTN)
      ENDIF
*----
*  Save new Volume/surface and albedo/material information
*----
      WRITE(IFTRKN) (VOLN(ITC),ITC=-NSOUTN,NVOUTN)
      WRITE(IFTRKN) (MATN(ITC),ITC=-NSOUTN,NVOUTN)
*----
*  Read and save BC and tracking directions and density
*----
      ALLOCATE(ICODE(NALBG),ALBD(NALBG),ANGLES(NDIM*NANGL),
     >         DENSTY(NANGL))
      READ(IFTRKO) (ICODE(ITC),ITC=1,NALBG)
      READ(IFTRKO) (ALBD(ITC),ITC=1,NALBG)
      READ(IFTRKO) (ANGLES(ITC),ITC=1,NDIM*NANGL)
      READ(IFTRKO) (DENSTY(ITC),ITC=1,NANGL)
      WRITE(IFTRKN) (ICODE(ITC),ITC=1,NALBG)
      WRITE(IFTRKN) (ALBD(ITC),ITC=1,NALBG)
      WRITE(IFTRKN) (ANGLES(ITC),ITC=1,NDIM*NANGL)
      WRITE(IFTRKN) (DENSTY(ITC),ITC=1,NANGL)
      DEALLOCATE(DENSTY,ANGLES,ALBD,ICODE)
      DEALLOCATE(MATN,VOLN)
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
 6010 FORMAT(' Initial surfaces and albedos')
 6011 FORMAT(' Final surfaces and albedos')
 6012 FORMAT(' Initial volumes and materials')
 6013 FORMAT(' Final volumes and materials')
 6020 FORMAT(1P,5(I10,E15.7))
*----
*  ABORT FORMATS
*----
 6100 FORMAT(' ------ ABORT IN ROUTINE ',A6,'  ------'/
     >       ' BOUNDARY CONDITIONS INCOMPATIBLE FOR MERGE '/
     >       ' NEW REGION = ',I10,5X,'SURFACE =',I10/
     >       ' OLD REGION = ',I10,5X,'SURFACE =',I10/
     >       ' ----------------------------------------')
 6101 FORMAT(' SURFACE DESCRIPTION '/
     >       4(' SURFACE  ->  ALBEDO'))
 6102 FORMAT(4(1X,I7,6X,I6))
 6200 FORMAT(' ------ ABORT IN ROUTINE ',A6,'  ------'/
     >       ' MATERIAL INCOMPATIBLE FOR MERGE '/
     >       ' NEW REGION = ',I10,5X,'MATERIAL =',I10/
     >       ' OLD REGION = ',I10,5X,'MATERIAL =',I10/
     >       ' ----------------------------------------')
 6201 FORMAT(' REGION DESCRIPTION '/
     >       4(' REGION  -> MATERIAL'))
 6202 FORMAT(4(1X,I6,5X,I8))
      END
