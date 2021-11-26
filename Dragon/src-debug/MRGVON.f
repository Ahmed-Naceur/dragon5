*DECK MRGVON
      SUBROUTINE MRGVON(IUPD  ,NSOUTO,NVOUTO,NSOUTN,NVOUTN,
     >                  NETSUR,NETVOL,NUNN  ,MAXMN ,
     >                  IMERGE,MATO  ,VOLO  ,MATRTO,
     >                  MATN  ,VOLN  ,KEYN  ,MATRTN,
     >                  NEXMAT,NEXKEY)
*
*----------
*
*Purpose:
* Merge volume and surface for NXT geometry.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Harrisson
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
* NETVOL  number of original regions.
* NETSUR  number of original surfaces.
* NUNN    new number of unknowns.
* IMERGE  merged position.
* MATO    old material per region.
* VOLO    old volumes.
* MATRTO  old B.C. conditions.
*
*Parameters: input/output
* NEXMAT  old/new NXTRecord MATALB for albedo number modification.
* NEXKEY  old/new KEYMRG index for NXT.
*
*Parameters: output
* MAXMN   new maximum number of mixture.
* MATN    new material per region.
* VOLN    new volumes.
* KEYN    new keyflux.
* MATRTN  new B.C. conditions.
*
*----------
*
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='MRGVON')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          IUPD(4),NSOUTO,NVOUTO,NSOUTN,NVOUTN,
     >                 NETSUR,NETVOL,NUNN,MAXMN
      INTEGER          IMERGE(-NSOUTO:NVOUTO),
     >                 MATO(NVOUTO),MATRTO(NSOUTO),
     >                 MATN(NVOUTN),KEYN(NUNN),MATRTN(NSOUTN),
     >                 NEXMAT(-NETSUR:NETVOL),
     >                 NEXKEY(-NETSUR:NETVOL)
      REAL             VOLO(NVOUTO),VOLN(NVOUTN)
*----
*  LOCAL VARIABLES
*----
      INTEGER          IVSN,IVSO
      DOUBLE PRECISION DVOL
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDABL
*----
*  CHANGE ORIGINAL AND/OR FINAL MATERIAL AND ORIGINAL ALBEDO IF REQUESTED
*----
      IF(IUPD(3) .GT. 0) THEN
        WRITE(IOUT,6300)
      ELSE IF(IUPD(3) .LT. 0) THEN
        WRITE(IOUT,6400)
      ELSE IF(IUPD(4) .GT. 0) THEN
        WRITE(IOUT,6500)
      ENDIF
*----
*  MERGE MATERIAL VOLUME AND KEY
*----
      IF(IUPD(1) .GT. 0) THEN
        DO IVSN=1,NVOUTN
          MATN(IVSN)=0
          DVOL=0.0D0
          DO IVSO=1,NVOUTO
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
          ENDDO
          VOLN(IVSN)=REAL(DVOL)
          KEYN(IVSN)=IVSN
        ENDDO
        DO IVSN=NVOUTN+1,NUNN
          KEYN(IVSN)=0
        ENDDO
        DO  IVSN=1,NETVOL
          DO  IVSO=1,NVOUTO
            IF(NEXKEY(IVSN) .EQ. IVSO) THEN
              NEXKEY(IVSN)=IMERGE(IVSO)
              GO TO 100
            ENDIF
          ENDDO
 100      CONTINUE
        ENDDO
      ELSE
*----
*  NO MERGE TRANSFER INFORMATION TO NEW VECTORS
*----
        DO IVSO=1,NVOUTO
          MATN(IVSO)=MATO(IVSO)
          VOLN(IVSO)=VOLO(IVSO)
        ENDDO
      ENDIF
*----
*  FIND NEW MAXIMUM NUMBER OF MIXTURE
*----
      MAXMN=0
      DO IVSN=1,NVOUTN
        MAXMN=MAX(MAXMN,MATN(IVSN))
      ENDDO
*----
*  MERGE REFLECTION/TRANSMISSION MATRIX
*----
      IF(IUPD(2).EQ.0) THEN
        DO IVSN=1,NSOUTO
          MATRTN(IVSN)=MATRTO(IVSN)
        ENDDO
      ELSE
        DO IVSN=-NSOUTN,-1,1
          DO IVSO=-NSOUTO,-1,1
            IF(IMERGE(IVSO).EQ.IVSN) THEN
              MATRTN(-IVSN)=-IMERGE(-MATRTO(-IVSO))
              GO TO 110
            ENDIF
          ENDDO
 110      CONTINUE
        ENDDO
*----
*  TEST IF MATRTN IS COHERENT
*----
        DO IVSN=1,NSOUTN
          IVSO=MATRTN(IVSN)
          IF(MATRTN(IVSO).NE.IVSN) THEN
            CALL XABORT(NAMSBR//
     >           ': SURFACES BC INCOMPATIBLE FOR MERGE')
          ENDIF
        ENDDO
        DO  IVSN=-1,-NETSUR,-1
          DO  IVSO=-1,-NSOUTO,-1
            IF(NEXKEY(IVSN) .EQ. IVSO) THEN
              NEXKEY(IVSN)=IMERGE(IVSO)
              GO TO 120
            ENDIF
          ENDDO
 120      CONTINUE
        ENDDO
*----
*  MERGING SURFACES WITH DIFFERENT ALBEDO NUMBER
*  USEFUL TO ACHIEVE SYME SYMMETRY
*----
         ALLOCATE(IDABL(NSOUTN))
         CALL XDISET(IDABL,NSOUTN,0)
         DO IVSN=1,NSOUTN
           DO IVSO=1,NETSUR
             IF (IMERGE(-IVSO) .EQ. -IVSN) THEN
               IF (IDABL(IVSN) .EQ. 0) THEN
                 IDABL(IVSN)=NEXMAT(-IVSO)
               ELSE
                 NEXMAT(-IVSO)=IDABL(IVSN)
               ENDIF
             ENDIF
           ENDDO
         ENDDO
         DEALLOCATE(IDABL)
      ENDIF
      RETURN
*----
*  FORMATS
*----
 6000 FORMAT(' ------ ABORT IN ROUTINE ',A6,'  ------'/
     >       ' MATERIAL INCOMPATIBLE FOR MERGE '/
     >       ' NEW REGION = ',I10,5X,'MATERIAL =',I10/
     >       ' OLD REGION = ',I10,5X,'MATERIAL =',I10/
     >       ' ----------------------------------------')
 6300 FORMAT(' ***** WARNING:  OPTION OLDM IS INVALID FOR GEOMETRIES'/
     >       '       TRACKED WITH NXT:. ORIGINAL MIXTURES ARE USED')
 6400 FORMAT(' ***** WARNING:  OPTION NEWM IS INVALID FOR GEOMETRIES'/
     >       '       TRACKED WITH NXT:. ORIGINAL MIXTURES ARE USED')
 6500 FORMAT(' ***** WARNING:  OPTION ALBE IS INVALID FOR GEOMETRIES'/
     >       '       TRACKED WITH NXT:. ORIGINAL ALBEDO ARE USED')
      END
