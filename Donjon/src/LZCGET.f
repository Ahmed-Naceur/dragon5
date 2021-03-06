*DECK LZCGET
      SUBROUTINE LZCGET(KPDEV,NTOT,NMIX,NTOT2,MIX,ID,LIMIT,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read the specification for a given liquid zone controller.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input/output
* KPDEV  pointer to DEV_LZC directory for lzc information.
* NTOT   old total number of all mixtures.
* NMIX   old maximum number of material mixtures.
* NTOT2  new total number of all mixtures.
* MIX    new mixture index of all mixtures.
* ID     current lzc identification number.
* LIMIT  core limiting coordinates.
* IMPX   printing index (=0 for no print).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KPDEV
      INTEGER NTOT,NMIX,NTOT2,MIX(NTOT2),ID,IMPX
      REAL LIMIT(6)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      INTEGER EMIX(2),FMIX(2)
      REAL MAXPOS(6),EMPTPOS(6),FULLPOS(6),LEVEL
      DOUBLE PRECISION DFLOT
      CHARACTER TEXT*12,AXIS
*----
*  REACTOR CORE LIMITS
*----
      XMIN=LIMIT(1)
      XMAX=LIMIT(2)
      YMIN=LIMIT(3)
      YMAX=LIMIT(4)
      ZMIN=LIMIT(5)
      ZMAX=LIMIT(6)
*----
*  WHOLE LZC POSITION
*----
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(TEXT.NE.'MAXPOS')CALL XABORT('@LZCGET: KEYWORD MAXPOS EXP'
     1 //'ECTED.')
      DO I=1,6
        CALL REDGET(ITYP,NITMA,MAXPOS(I),TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@LZCGET: REAL FOR MAXPOS EXPECTED.')
      ENDDO
*----
*  CHECK LZC POSITION
*----
      IF(MAXPOS(2).LT.MAXPOS(1))CALL XABORT('@LZCGET: WRONG X '
     1 //'LZC COORDINATES: X- > X+')
      IF(MAXPOS(1).LT.XMIN)CALL XABORT('@LZCGET: WRONG X- VALUE.')
      IF(MAXPOS(2).GT.XMAX)CALL XABORT('@LZCGET: WRONG X+ VALUE.')
*
      IF(MAXPOS(4).LT.MAXPOS(3))CALL XABORT('@LZCGET: WRONG Y '
     1 //'LZC COORDINATES: Y- > Y+')
      IF(MAXPOS(3).LT.YMIN)CALL XABORT('@LZCGET: WRONG Y- VALUE.')
      IF(MAXPOS(4).GT.YMAX)CALL XABORT('@LZCGET: WRONG Y+ VALUE.')
*
      IF(MAXPOS(6).LT.MAXPOS(5))CALL XABORT('@LZCGET: WRONG Z '
     1 //'LZC COORDINATES: Z- > Z+')
      IF(MAXPOS(5).LT.ZMIN)CALL XABORT('@LZCGET: WRONG Z- VALUE.')
      IF(MAXPOS(6).GT.ZMAX)CALL XABORT('@LZCGET: WRONG Z+ VALUE.')
*----
*  MAX-FULL COORDINATE
*----
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(TEXT.NE.'MAX-FULL')CALL XABORT('@LZCGET: KEYWORD MAX-FULL'
     1 //' EXPECTED.')
      CALL REDGET(ITYP,NITMA,FULMAX,TEXT,DFLOT)
      IF(ITYP.NE.2)CALL XABORT('@LZCGET: REAL FOR MAX-FULL COORDIN'
     1 //'ATE EXPECTED.')
*----
*  LZC FILLING AXIS
*----
      HEIGHT=0.
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(TEXT.NE.'AXIS')CALL XABORT('@LZCGET: KEYWORD AXIS EXPECTED.')
      CALL REDGET(ITYP,NITMA,FLOT,AXIS,DFLOT)
      IF(AXIS.NE.'X')THEN
      IF(AXIS.NE.'Y')THEN
      IF(AXIS.NE.'Z')THEN
       CALL XABORT('@LZCGET: X, Y OR Z EXPECTED FOR AXIS.')
      ELSE
        IAXIS=3
        IF(FULMAX.GT.MAXPOS(4))CALL XABORT('@LZCGET: WRONG MAX-FULL VA'
     1  //'LUE: > Z+.')
        IF(FULMAX.LT.MAXPOS(3))CALL XABORT('@LZCGET: WRONG MAX-FULL VA'
     1  //'LUE: < Z-.')
        HEIGHT=MAXPOS(6)-FULMAX
      ENDIF
      ELSE
        IAXIS=2
        IF(FULMAX.GT.MAXPOS(4))CALL XABORT('@LZCGET: WRONG MAX-FULL VA'
     1  //'LUE: > Y+.')
        IF(FULMAX.LT.MAXPOS(3))CALL XABORT('@LZCGET: WRONG MAX-FULL VA'
     1  //'LUE: < Y-.')
        HEIGHT=MAXPOS(4)-FULMAX
      ENDIF
      ELSE
        IAXIS=1
        IF(FULMAX.GT.MAXPOS(2))CALL XABORT('@LZCGET: WRONG MAX-FULL VA'
     1  //'LUE: > X+.')
        IF(FULMAX.LT.MAXPOS(1))CALL XABORT('@LZCGET: WRONG MAX-FULL VA'
     1  //'LUE: < X-.')
        HEIGHT=MAXPOS(2)-FULMAX
      ENDIF
      IF(HEIGHT.EQ.0.)CALL XABORT('@LZCGET: MAX-FULL WATER HEIGHT =0.')
*----
*  LZC FILLING LEVEL
*----
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(TEXT.NE.'LEVEL')CALL XABORT('@LZCGET: KEYWORD LEVEL EXPECTED.')
      CALL REDGET(ITYP,NITMA,LEVEL,TEXT,DFLOT)
      IF(ITYP.NE.2)CALL XABORT('@LZCGET: REAL FOR FILLING LEVEL EX'
     1 //'PECTED.')
      IF(LEVEL.GT.1.)CALL XABORT('@LZCGET: WRONG FILLING LEVEL: > 1.')
      IF(LEVEL.LT.0.)CALL XABORT('@LZCGET: WRONG FILLING LEVEL: < 0.')
*----
*  LZC FILLING RATE
*----
      RATE=0.
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@LZCGET: CHARACTER DATA EXPECTED(1).')
      IF(TEXT.NE.'RATE')GOTO 10
      CALL REDGET(ITYP,NITMA,RATE,TEXT,DFLOT)
      IF(ITYP.NE.2)CALL XABORT('@LZCGET: REAL FOR RATE EXPECTED.')
      IF(RATE.LT.0.)CALL XABORT('@DEVSET: WRONG RATE VALUE < 0.')
*----
*  LZC FILLING TIME
*----
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@LZCGET: CHARACTER DATA EXPECTED(2).')
   10 TIME=0.
      IF(TEXT.NE.'TIME')GOTO 20
      CALL REDGET(ITYP,NITMA,TIME,TEXT,DFLOT)
      IF(ITYP.NE.2)CALL XABORT('@LZCGET: REAL FOR TIME EXPECTED.')
      IF(TIME.LT.0.)CALL XABORT('@DEVSET: WRONG TIME VALUE < 0.')
*----
*  LZC MIXTURES
*----
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@LZCGET: CHARACTER DATA EXPECTED(3).')
*     EMPTY PART
   20 IF(TEXT.NE.'EMPTY-MIX')CALL XABORT('@LZCGET: KEYWORD EMPTY-MI'
     1 //'X EXPECTED.')
      DO I=1,2
        CALL REDGET(ITYP,EMIX(I),FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@LZCGET: INTEGER EMPTY-MIX NUMBER'
     1 //' EXPECTED.')
        MIX(NTOT+(ID-1)*4+I)=EMIX(I)
        EMIX(I)=NMIX+(ID-1)*4+I
      ENDDO
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@LZCGET: CHARACTER DATA EXPECTED(4).')
*     FULL PART
      IF(TEXT.NE.'FULL-MIX')CALL XABORT('@LZCGET: KEYWORD FULL-MIX '
     1 //'EXPECTED.')
      DO I=1,2
        CALL REDGET(ITYP,FMIX(I),FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@LZCGET: INTEGER FULL-MIX NUMBER '
     1 //'EXPECTED.')
        MIX(NTOT+(ID-1)*4+I+2)=FMIX(I)
        FMIX(I)=NMIX+(ID-1)*4+I+2
      ENDDO
*----
*  CURRENT LZC POSITION
*----
      DELH=LEVEL*HEIGHT
      DO I=1,6
        EMPTPOS(I)=MAXPOS(I)
        FULLPOS(I)=MAXPOS(I)
      ENDDO
      IF(IAXIS.EQ.1)THEN
        FULLPOS(1)=MAXPOS(2)-DELH
        EMPTPOS(2)=FULLPOS(1)
      ELSEIF(IAXIS.EQ.2)THEN
        FULLPOS(3)=MAXPOS(4)-DELH
        EMPTPOS(4)=FULLPOS(3)
      ELSEIF(IAXIS.EQ.3)THEN
        FULLPOS(5)=MAXPOS(6)-DELH
        EMPTPOS(6)=FULLPOS(5)
      ENDIF
*----
*  STORE LZC DATA
*----
      CALL LCMPUT(KPDEV,'LZC-ID',1,1,ID)
      CALL LCMPUT(KPDEV,'MAX-POS',6,2,MAXPOS)
      CALL LCMPUT(KPDEV,'AXIS',1,1,IAXIS)
      CALL LCMPUT(KPDEV,'HEIGHT',1,2,HEIGHT)
      CALL LCMPUT(KPDEV,'LEVEL',1,2,LEVEL)
      CALL LCMPUT(KPDEV,'EMPTY-POS',6,2,EMPTPOS)
      CALL LCMPUT(KPDEV,'FULL-POS',6,2,FULLPOS)
      CALL LCMPUT(KPDEV,'EMPTY-MIX',2,1,EMIX)
      CALL LCMPUT(KPDEV,'FULL-MIX',2,1,FMIX)
      CALL LCMPUT(KPDEV,'RATE',1,2,RATE)
      CALL LCMPUT(KPDEV,'TIME',1,2,TIME)
*
      IF(IMPX.GT.1)WRITE(IOUT,1000)MAXPOS(1),MAXPOS(3),MAXPOS(5),
     1 MAXPOS(2),MAXPOS(4),MAXPOS(6),HEIGHT,AXIS,LEVEL,EMPTPOS(1),
     2 EMPTPOS(3),EMPTPOS(5),EMPTPOS(2),EMPTPOS(4),EMPTPOS(6),
     3 FULLPOS(1),FULLPOS(3),FULLPOS(5),FULLPOS(2),FULLPOS(4),
     4 FULLPOS(6),RATE,TIME
      RETURN
*
 1000 FORMAT(/5X,'WHOLE POSITION :',
     1   4X,'X-',F10.4,5X,'Y-',F10.4,5X,'Z-',F10.4/
     2  37X,'X+',F10.4,5X,'Y+',F10.4,5X,'Z+',F10.4/
     3  /5X,'FIL-HEIGHT =',F9.4/
     4  /5X,'FIL-AXIS : ',A1,5X,'FIL-LEVEL =',F8.4/
     5  /5X,'EMPTY-PART POSITION :',
     6   5X,'X-',F10.4,5X,'Y-',F10.4,5X,'Z-',F10.4/
     7  32X,'X+',F10.4,5X,'Y+',F10.4,5X,'Z+',F10.4/
     8  /5X,'FULL-PART POSITION :',
     9   5X,'X-',F10.4,5X,'Y-',F10.4,5X,'Z-',F10.4/
     1  32X,'X+',F10.4,5X,'Y+',F10.4,5X,'Z+',F10.4/
     2  /5X,'FIL-RATE =',E11.4,5X,'FIL-TIME =',E11.4/)
      END
