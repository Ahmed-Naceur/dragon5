*DECK DMAGET
      SUBROUTINE DMAGET(IPDMA,NGRP,NFREG,NBMIX,MATCOD,IPRINT,NMERGE,
     1 NGCOND)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read from the input file the DMAC: module input options.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPDMA   pointer to the DMA data structure.
* NGRP    number of energy groups.
* NFREG   number of regions.
* NBMIX   maximum number of mixtures.
* MATCOD  region material.
*
*Parameters: input/output
* IPRINT  print parameter.
* NMERGE  number of merged regions.
* NGCOND  number of condensed energy groups.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDMA
      INTEGER NGRP,NFREG,NBMIX,MATCOD(NFREG),IPRINT,NMERGE,NGCOND
*----
*  LOCAL VARIABLES
*----
      INTEGER INDIC,NITMA,IREGIO,IMATER,IGROUP,JGROUP
      REAL FLOTT
      CHARACTER TEXT*12
      DOUBLE PRECISION DFLOTT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IMERGE,MIXMER,IGCR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IMERGE(NFREG),MIXMER(0:NBMIX),IGCR(NGRP))
*----
*  READ INPUT
*----
      IPRINT=1
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT,DFLOTT)
      IF(INDIC.EQ.10) GO TO 200
      IF(INDIC.NE.3) CALL XABORT('DMAGET: CHARACTER DATA EXPECTED(1)')
      IF(TEXT(1:4).EQ.'EDIT') THEN
        CALL REDGET(INDIC,IPRINT,FLOTT,TEXT,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('DMAGET: INTEGER DATA EXPECTED FOR'
     <  //' IPRINT')
      ELSE IF(TEXT(1:5).EQ.'RATE') THEN
*       DEFINE A TALLY
        NMERGE=0
        NGCOND=0
   20   CALL REDGET(INDIC,NITMA,FLOTT,TEXT,DFLOTT)
   30   IF(INDIC.NE.3) CALL XABORT('DMAGET: CHARACTER DATA EXPECTED(2)')
        IF(TEXT(:4).EQ.'MERG') THEN
*----
*  MERGING DIRECTIVE ANALYSIS
*----
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT,DFLOTT)
          IF(INDIC.NE.3) CALL XABORT('DMAGET: CHARACTER DATA EXPECTED'
     <    //'(3)')
          IF(TEXT.EQ.'COMP') THEN
*----
*  COMPLETE MERGE
*----
            CALL XDISET(IMERGE,NFREG,1)
            NMERGE=1
            GO TO 20
          ELSE IF(TEXT.EQ.'MIX') THEN
*----
*  MERGE BY MIXTURES
*----
            DO 40 IMATER=0,NBMIX
              MIXMER(IMATER)=IMATER
   40       CONTINUE
            DO 50 IREGIO=1,NFREG
              NMERGE=MAX(NMERGE,MATCOD(IREGIO))
              IMERGE(IREGIO)=MIXMER(MATCOD(IREGIO))
   50       CONTINUE
            NMERGE=NBMIX
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT,DFLOTT)
            IF(INDIC.EQ.1) THEN
*----
*  SPECIFY MIXTURES TO BE MERGED
*----
              NMERGE=MAX(0,NITMA)
              MIXMER(1)=NITMA
              DO 60 IMATER=2,NBMIX
                CALL REDGET(INDIC,NITMA,FLOTT,TEXT,DFLOTT)
                IF(INDIC.NE.1) CALL XABORT('DMAGET: INTEGER DATA EXPEC'
     <          //'TED FOR IMIXM')
                NMERGE=MAX(NMERGE,NITMA)
                MIXMER(IMATER)=NITMA
   60         CONTINUE
              DO 70 IREGIO=1,NFREG
                IMERGE(IREGIO)=MIXMER(MATCOD(IREGIO))
   70         CONTINUE
            ELSE IF(INDIC.EQ.3) THEN
*----
*  ASSOCIATE ONE REGION BY MIXTURE
*----
              GO TO 30
            ELSE
              CALL XABORT('DMAGET: READ ERROR - INVALID TYPE READ')
            ENDIF
          ELSE IF(TEXT.EQ.'REGI') THEN
*----
*  MERGE BY REGIONS
*----
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT,DFLOTT)
            IF(INDIC.NE.1) CALL XABORT('DMAGET: INTEGER DATA EXPEC'
     <          //'TED FOR IREGM')
            NMERGE=MAX(0,NITMA)
            IMERGE(1)=NITMA
            DO 80 IREGIO=2,NFREG
              CALL REDGET(INDIC,NITMA,FLOTT,TEXT,DFLOTT)
              IF(INDIC.NE.1) CALL XABORT('DMAGET: INTEGER DATA EXPECTE'
     <        //'D FOR IREGM')
              NMERGE=MAX(NMERGE,NITMA)
              IMERGE(IREGIO)=NITMA
   80       CONTINUE
          ELSE IF(TEXT.EQ.'NONE') THEN
*----
*  NO MERGING
*----
            NMERGE=NFREG
            DO 90 IREGIO=1,NFREG
              IMERGE(IREGIO)=IREGIO
   90       CONTINUE
          ELSE
            CALL XABORT('DMAGET: '//TEXT//' IS AN INVALID KEYWORD(1)')
          ENDIF
        ELSE IF(TEXT(:4).EQ.'COND') THEN
*----
*  GROUP CONDENSATION DIRECTIVE ANALYSIS
*----
          DO 110 IGROUP=1,NGRP+1
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT,DFLOTT)
            IF(INDIC.EQ.3) THEN
              IF(IGROUP.EQ.1) THEN
                IF(TEXT.EQ.'NONE') THEN
                  NGCOND=NGRP
                  DO 100 JGROUP=1,NGRP
                    IGCR(JGROUP)=JGROUP
  100             CONTINUE
                  GO TO 20
                ELSE
                  NGCOND=1
                  IGCR(NGCOND)=NGRP
                ENDIF
              ENDIF
              IF(IGCR(NGCOND).NE.NGRP) THEN
                 NGCOND=NGCOND+1
                 IGCR(NGCOND)=NGRP
              ENDIF
              GO TO 30
            ELSE IF(INDIC.EQ.1) THEN
              IF(NITMA.GT.NGRP) NITMA=NGRP
              IF(NGCOND.GT.0) THEN
                IF(NITMA.GT.IGCR(NGCOND)) THEN
                  NGCOND=NGCOND+1
                  IGCR(NGCOND)=NITMA
                ENDIF
              ELSE
                NGCOND=NGCOND+1
                IGCR(NGCOND)=NITMA
              ENDIF
            ENDIF
  110     CONTINUE
        ELSE IF(TEXT(:4).EQ.'ENDR') THEN
          GO TO 120
        ELSE 
          CALL XABORT('DMAGET: '//TEXT//' IS AN INVALID KEYWORD(2)')
        ENDIF
        GO TO 20   
  120   CALL LCMPUT(IPDMA,'REF:IMERGE',NFREG,1,IMERGE)
        CALL LCMPUT(IPDMA,'REF:IGCOND',NGCOND,1,IGCR)
      ELSE IF(TEXT(1:1).EQ.';') THEN
        GO TO 200
      ELSE 
        CALL XABORT('DMAGET: '//TEXT//' IS AN INVALID KEYWORD(3)')
      ENDIF
      GO TO 10
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  200 DEALLOCATE(IGCR,MIXMER,IMERGE)
      RETURN
      END
