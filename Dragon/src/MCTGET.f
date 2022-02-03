*DECK MCTGET
      SUBROUTINE MCTGET(IPOUT,NGRP,NFREG,NBMIX,MATCOD,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To read from the input file the MC: module input options.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): B. Arsenault
*
*Parameters: input
* IPOUT   pointer to the MC: data structure.
* NGRP    number of energy groups.
* NFREG   number of regions.
* NBMIX   maximum number of mixtures.
* MATCOD  region material.
*
*Parameters: input/output
* IPRINT  print parameter.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPOUT
      INTEGER NSTATE,NGRP,NFREG,NBMIX,MATCOD(NFREG),IPRINT
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE),INDIC,NITMA,NMERGE,IREGIO,IMATER,NGCOND,
     1 IGROUP,JGROUP
      REAL FLOTT
      CHARACTER TEXT*12
      DOUBLE PRECISION DFLOTT
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IMERGE,MIXMER,IGCR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IMERGE(NFREG),MIXMER(0:NBMIX),IGCR(NGRP))
*----
*  READ INPUT
*----
      IPRINT=1
      CALL XDISET(ISTATE,NSTATE,0)
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT,DFLOTT)
      IF(INDIC.EQ.10) GO TO 200
      IF(INDIC.NE.3) CALL XABORT('MCTGET: CHARACTER DATA EXPECTED(1)')
      IF(TEXT(1:4).EQ.'EDIT') THEN
        CALL REDGET(INDIC,IPRINT,FLOTT,TEXT,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('MCTGET: INTEGER DATA EXPECTED FOR'
     <  //' IPRINT')
      ELSE IF(TEXT(1:5).EQ.'KCODE') THEN
*       READ THE NSRCK PARAMETER
        CALL REDGET(INDIC,ISTATE(1),FLOTT,TEXT,DFLOTT)
        IF (INDIC.NE.1) CALL XABORT('MCTGET: INTEGER DATA EXPECTED FOR'
     <  //' NSRCK') 
*       READ THE IKZ PARAMETER
        CALL REDGET(INDIC,ISTATE(2),FLOTT,TEXT,DFLOTT)
        IF (INDIC.NE.1) CALL XABORT('MCTGET: INTEGER DATA EXPECTED FOR'
     <  //' IKZ') 
*       READ THE KCT PARAMETER
        CALL REDGET(INDIC,ISTATE(3),FLOTT,TEXT,DFLOTT)
        IF (INDIC.NE.1) CALL XABORT('MCTGET: INTEGER DATA EXPECTED FOR'
     <  //' KCT')
      ELSE IF(TEXT(1:4).EQ.'SEED') THEN
*       INPUT A SEED INTEGER
        CALL REDGET(INDIC,ISTATE(4),FLOTT,TEXT,DFLOTT)
        IF (INDIC.NE.1) CALL XABORT('MCTGET: INTEGER DATA EXPECTED FOR'
     <  //' SEED')
      ELSE IF(TEXT(1:3).EQ.'N2N') THEN
*       N2N CROSS SECTION RECOVERY FLAG
        ISTATE(5)=1
      ELSE IF(TEXT(1:5).EQ.'TALLY') THEN
*       DEFINE A TALLY
        IF(ISTATE(6).NE.0) CALL XABORT('MCTGET: TALLY EXISTS')
        ISTATE(6)=1
        NMERGE=0
        NGCOND=0
   20   CALL REDGET(INDIC,NITMA,FLOTT,TEXT,DFLOTT)
   30   IF(INDIC.NE.3) CALL XABORT('MCTGET: CHARACTER DATA EXPECTED(2)')
        IF(TEXT(:4).EQ.'MERG') THEN
*----
*  MERGING DIRECTIVE ANALYSIS
*----
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT,DFLOTT)
          IF(INDIC.NE.3) CALL XABORT('MCTGET: CHARACTER DATA EXPECTED'
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
                IF(INDIC.NE.1) CALL XABORT('MCTGET: INTEGER DATA EXPEC'
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
              CALL XABORT('MCTGET: READ ERROR - INVALID TYPE READ')
            ENDIF
          ELSE IF(TEXT.EQ.'REGI') THEN
*----
*  MERGE BY REGIONS
*----
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT,DFLOTT)
            IF(INDIC.NE.1) CALL XABORT('MCTGET: INTEGER DATA EXPEC'
     <          //'TED FOR IREGM')
            NMERGE=MAX(0,NITMA)
            IMERGE(1)=NITMA
            DO 80 IREGIO=2,NFREG
              CALL REDGET(INDIC,NITMA,FLOTT,TEXT,DFLOTT)
              IF(INDIC.NE.1) CALL XABORT('MCTGET: INTEGER DATA EXPECTE'
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
            CALL XABORT('MCTGET: '//TEXT//' IS AN INVALID KEYWORD(1)')
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
        ELSE IF(TEXT(:4).EQ.'ENDT') THEN
          GO TO 120
        ELSE 
          CALL XABORT('MCTGET: '//TEXT//' IS AN INVALID KEYWORD(2)')
        ENDIF
        GO TO 20   
  120   CALL LCMPUT(IPOUT,'REF:IMERGE',NFREG,1,IMERGE)
        CALL LCMPUT(IPOUT,'REF:IGCOND',NGCOND,1,IGCR)
        IF((NMERGE.GT.0).AND.(NGCOND.GT.0)) ISTATE(6)=2
        ISTATE(7)=NMERGE
        ISTATE(8)=NGCOND
      ELSE IF(TEXT(1:1).EQ.';') THEN
        GO TO 200
      ELSE 
        CALL XABORT('MCTGET: '//TEXT//' IS AN INVALID KEYWORD(3)')
      ENDIF
      GO TO 10
  200 ISTATE(9)=NFREG
      CALL LCMPUT(IPOUT,'STATE-VECTOR',NSTATE,1,ISTATE)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IGCR,MIXMER,IMERGE)
      RETURN
      END
