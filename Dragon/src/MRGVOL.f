*DECK MRGVOL
      SUBROUTINE MRGVOL(IUPD  ,NSOUTO,NVOUTO,NSOUTN,NVOUTN,NUNN  ,
     >                  IMERGE,MIXN  ,MATO  ,VOLO  ,MATN  ,VOLN  ,
     >                  KEYN  ,MATRTO,MATRTN,MAXMN ,NETVOL,NETSUR,
     >                  MATRO ,KEYRO ,MATRN ,KEYRN )
*
*----------
*
*Purpose:
* Merge information on data structure.
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
* IUPD    type of merge required:
*         IUPD(1) for region merge;
*         IUPD(2) for surface merge;
*         IUPD(3) for material merge;
*         IUPD(4) for albedo merge.
* NSOUTO  old number of surfaces.
* NVOUTO  old number of regions.
* NSOUTN  new number of surfaces.
* NVOUTN  new number of regions.
* NUNN    new number of unknowns.
* IMERGE  merged position.
* MIXN    new material for old regions.
* MATO    old material per region.
* VOLO    old volumes.
* MATRTO  old B.C. conditions.
* NETVOL  number of original regions.
* NETSUR  number of original surfaces.
* MATRO   old regional MATALB.
* KEYRO   old regional KEYMRG.
*
*Parameters: output
* MATN    new material per region.
* VOLN    new volumes.
* KEYN    new keyflux.
* MATRTN  new B.C. conditions.
* MAXMN   new maximum number of mixture.
* MATRN   new regional MATALB.
* KEYRN   new regional KEYMRG.
*
*----------
*
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='MRGVOL')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          IUPD(4),NSOUTO,NVOUTO,NSOUTN,NVOUTN,NUNN,
     >                 MAXMN,NETVOL,NETSUR
      INTEGER          IMERGE(-NSOUTO:NVOUTO),MIXN(NVOUTO),
     >                 MATO(NVOUTO),MATN(NVOUTN),KEYN(NUNN),
     >                 MATRTO(NSOUTO),MATRTN(NSOUTN),
     >                 MATRO(-NETSUR:NETVOL),KEYRO(-NETSUR:NETVOL),
     >                 MATRN(-NSOUTN:NVOUTN),KEYRN(-NSOUTN:NVOUTN)
      REAL             VOLO(NVOUTO),VOLN(NVOUTN)
*----
*  LOCAL VARIABLES
*----
      INTEGER          IVSN,IVSO
      DOUBLE PRECISION DVOL
*----
*  TRANSFER OLD KEYMRG AND MATALB TO NEW VECTOR
*----
      DO 90 IVSN=-NETSUR,NETVOL
        KEYRN(IVSN)=KEYRO(IVSN)
        MATRN(IVSN)=MATRO(IVSN)
 90   CONTINUE
*----
*  CHANGE ORIGINAL MATERIAL IF REQUESTED
*----
      IF(IUPD(3) .GT. 0) THEN
        DO 100 IVSN=1,IUPD(3)
          MATO(IVSN)=MIXN(IVSN)
          DO 101 IVSO=1,NETVOL
            IF(KEYRO(IVSO) .EQ. IVSN) THEN
              MATRN(IVSO)=MATO(IVSN)
            ENDIF
 101      CONTINUE
 100    CONTINUE
      ENDIF
      IF(IUPD(1) .GT. 0) THEN
*----
*  MERGE MATERIAL VOLUME AND KEY
*----
        DO 110 IVSN=1,NVOUTN
          MATN(IVSN)=0
          DVOL=0.0D0
          DO 111 IVSO=1,NVOUTO
           IF(IMERGE(IVSO) .EQ. IVSN) THEN
             IF(MATN(IVSN) .EQ. 0) THEN
               MATN(IVSN)=MATO(IVSO)
             ELSE IF(MATN(IVSN) .NE. MATO(IVSO))THEN
               WRITE(IOUT,6000) NAMSBR,IVSN,MATN(IVSN),IVSO,MATO(IVSO)
               CALL XABORT(NAMSBR//
     >           ': MATERIAL INCOMPATIBLE FOR MERGE')
             ENDIF
             DVOL=DVOL+DBLE(VOLO(IVSO))
           ENDIF
 111      CONTINUE
          VOLN(IVSN)=REAL(DVOL)
          KEYN(IVSN)=IVSN
 110    CONTINUE
        DO 112 IVSN=NVOUTN+1,NUNN
          KEYN(IVSN)=0
 112    CONTINUE
        DO 113 IVSO=1,NVOUTO
          DO 114 IVSN=1,NETVOL
            IF(KEYRO(IVSN) .EQ. IVSO) THEN
              KEYRN(IVSN)=IMERGE(IVSO)
            ENDIF
 114      CONTINUE
 113    CONTINUE
      ELSE
*----
*  NO MERGE TRANSFER INFORMATION TO NEW VECTORS
*----
        DO 120 IVSO=1,NVOUTO
          MATN(IVSO)=MATO(IVSO)
          VOLN(IVSO)=VOLO(IVSO)
 120    CONTINUE
      ENDIF
*----
*  CHANGE FINAL MATERIAL IF REQUESTED
*----
      IF(IUPD(3) .LT. 0) THEN
        DO 130 IVSN=1,-IUPD(3)
          MATN(IVSN)=MIXN(IVSN)
          DO 131 IVSO=1,NETVOL
            IF(KEYRO(IVSO) .EQ. IVSN) THEN
              MATRN(IVSO)=MIXN(IVSN)
            ENDIF
 131      CONTINUE
 130    CONTINUE
      ENDIF
*----
*  FIND NEW MAXIMUM NUMBER OF MIXTURE
*----
      MAXMN=0
      DO 140 IVSN=1,NVOUTN
        MAXMN=MAX(MAXMN,MATN(IVSN))
 140  CONTINUE
*----
*  MERGE REFLECTION/TRANSMISSION MATRIX
*----
      IF(IUPD(2).EQ.0) THEN
        DO 150 IVSN=1,NSOUTO
          MATRTN(IVSN)=MATRTO(IVSN)
 150    CONTINUE
      ELSE
        DO 160 IVSN=-NSOUTN,-1,1
          DO 161 IVSO=-NSOUTO,-1,1
            IF(IMERGE(IVSO).EQ.IVSN) THEN
              MATRTN(-IVSN)=-IMERGE(-MATRTO(-IVSO))
              GO TO 165
            ENDIF
 161      CONTINUE
 165      CONTINUE
 160    CONTINUE
*----
*  TEST IF MATRTN IS COHERENT
*----
        DO 162 IVSN=1,NSOUTN
          IVSO=MATRTN(IVSN)
          IF(MATRTN(IVSO).NE.IVSN) THEN
            CALL XABORT(NAMSBR//
     >           ': SURFACES BC INCOMPATIBLE FOR MERGE')
          ENDIF
 162    CONTINUE
        DO 163 IVSO=-1,-NSOUTO,-1
          DO 164 IVSN=-1,-NETSUR,-1
            IF(KEYRO(IVSN) .EQ. IVSO) THEN
              KEYRN(IVSN)=IMERGE(IVSO)
            ENDIF
 164      CONTINUE
 163    CONTINUE
      ENDIF
      RETURN
*----
*  ABORT FORMATS
*----
 6000 FORMAT(' ------ ABORT IN ROUTINE ',A6,'  ------'/
     >       ' MATERIAL INCOMPATIBLE FOR MERGE '/
     >       ' NEW REGION = ',I10,5X,'MATERIAL =',I10/
     >       ' OLD REGION = ',I10,5X,'MATERIAL =',I10/
     >       ' ----------------------------------------')
      END
