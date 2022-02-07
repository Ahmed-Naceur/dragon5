*DECK GEODMI
      SUBROUTINE GEODMI(LX,LY,LZ,LCOUR,MIX,MINMIX,MAXMIX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Build array MIX from plane-defined information.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): E. Varin and R. Roy
*
*Parameters: input
* LX      number of meshes along X-axis.
* LY      number of meshes along Y-axis.
* LZ      number of meshes along Z-axis.
* LCOUR   flag indicating if 'CROWN' or 'UPTO' keywords are allowed.
* MIX     array of material mixtures.
*
*Parameters: output
* MINMIX  minimum number of mixtures, considering all sub-geometries.
* MAXMIX  maximum number of mixtures, considering all sub-geometries.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER    LX,LY,LZ,MIX(LX,LY,LZ),MINMIX,MAXMIX
      LOGICAL    LCOUR
*----
*  LOCAL VARIABLES
*----
      INTEGER    IZ1,IZ2,IZ3,IX,IY,NZ,NC,NCSAME,IC,INDIC,NITMA,IHEX
      CHARACTER  TEXT12*12
      REAL       FLOTT
      DOUBLE PRECISION DFLOTT
*----
*  READ AN OPTION KEYWORD
*----
      IHEX= 0
      NC = -1
      NZ =  0
  5   IZ2 = 0
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF(INDIC.EQ.1) THEN
         NZ = NZ + 1
         IZ1 = NITMA
         IF( IZ1.LT.1.OR.IZ1.GT.LZ )THEN
            CALL XABORT('GEODMI: INVALID PLANE NUMBER'//
     >                  '(GREATER THAN *LZ*).')
         ENDIF
      ELSE
         CALL XABORT('GEODMI: PLANE NUMBER MUST BE READ'//
     >               '(INTEGER EXPECTED).')
      ENDIF
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF( INDIC.EQ.3 )THEN
         IF(TEXT12.EQ.'SAME') THEN
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.EQ.1) THEN
               IZ2 = NITMA
               IF( IZ2.GT.IZ1 )THEN
                  CALL XABORT('GEODMI: INVALID PLANE NUMBER'//
     >                        '(GREATER THAN PREVIOUS).')
               ENDIF
               GOTO 20
            ELSE
               CALL XABORT('GEODMI: SAME AS WHICH PLANE? '//
     >               '(INTEGER EXPECTED).')
            ENDIF
         ELSEIF((TEXT12.EQ.'CROWN'.OR.TEXT12.EQ.'UPTO').AND.LCOUR) THEN
            NCSAME= 1
            IF( TEXT12.EQ.'UPTO' )THEN
               CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
               IF(INDIC.NE.1) CALL XABORT('GEODMI: INTEGER DATA'//
     >                       ' EXPECTED AFTER *UPTO* KEYWORD')
               NCSAME= NITMA
            ENDIF
            GO TO 30
         ELSEIF(.NOT.LCOUR.AND.(TEXT12.EQ.'CROWN'.OR.TEXT12.EQ.'UPTO'))
     >         THEN
            CALL XABORT('GEODMI: UNSUPPORTED KEYWORD *CROWN* OR *UPTO*'
     >                 //': HEX3D COMPLETE ONLY')
         ELSE
            CALL XABORT('GEODMI: INVALID CHARACTER VARIABLE '//TEXT12)
         ENDIF
      ELSEIF (INDIC.EQ.1) THEN
         GOTO 20
      ELSE
         CALL XABORT('GEODMI: INTEGER OR CHARACTER VARIABLE EXPECTED')
      ENDIF
*----
*  READ A CWOWN
*----
 30   CONTINUE
      NC= NC+1
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF( INDIC.EQ.3 )THEN
         IF(TEXT12.EQ.'SAME') THEN
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.EQ.1) THEN
               IZ3= NITMA
               IF( IZ3.GT.IZ1 )THEN
                  CALL XABORT('GEODMI: INVALID PLANE NUMBER'//
     >                        '(GREATER THAN PREVIOUS).')
               ENDIF
               DO 41 IC= 1, NCSAME
                  IF( NC.EQ.0 )THEN
                     MIX(IHEX+1,1,IZ1)= MIX(IHEX+1,1,IZ3)
                     IHEX= IHEX+1
                  ELSE
                     DO 31 IX= IHEX+1, IHEX+6*NC
                        MIX(IX,1,IZ1)= MIX(IX,1,IZ3)
 31                  CONTINUE
                     IHEX= IHEX+6*NC
                  ENDIF
                  NC= NC+1
 41            CONTINUE
               NC= NC -1
            ELSE
               CALL XABORT('GEODMI: SAME AS WHICH PLANE? '//
     >               '(INTEGER EXPECTED).')
            ENDIF
         ELSEIF(TEXT12.EQ.'ALL') THEN
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.EQ.1) THEN
               MAXMIX=MAX(MAXMIX,NITMA)
               MINMIX=MIN(MINMIX,NITMA)
               DO 42 IC= 1, NCSAME
                  IF( NC.EQ.0 )THEN
                     MIX(IHEX+1,1,IZ1)= NITMA
                     IHEX= IHEX+1
                  ELSE
                     DO 32 IX= IHEX+1, IHEX+6*NC
                        MIX(IX,1,IZ1)= NITMA
 32                  CONTINUE
                     IHEX= IHEX+6*NC
                  ENDIF
                  NC= NC+1
 42            CONTINUE
               NC= NC -1
            ELSE
               CALL XABORT('GEODMI: ALL OF WHICH MIX? '//
     >               '(INTEGER EXPECTED).')
            ENDIF
          ELSE
            CALL XABORT('GEODMI: KEYWORD *SAME* OR *ALL* '//
     >                  '(CHARACTER EXPECTED).')
          ENDIF
      ELSEIF( INDIC.EQ.1 )THEN
        IF( NCSAME.NE.1 )THEN
           CALL XABORT('GEODMI: INVALID INTEGER WITH *UPTO* ')
        ENDIF
        IF( NC.EQ.0 )THEN
           MIX(IHEX+1,1,IZ1)= NITMA
           IHEX= IHEX+1
        ELSE
         DO 33 IX= 1, 6*NC
         IF(.NOT.(IX.EQ.1) ) THEN
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          IF(INDIC.NE.1)THEN
           WRITE(6,*) 'NC=',NC,' IZ1=',IZ1,' NCSAME=',NCSAME
           WRITE(6,*) 'IHEX=',IHEX,' INDIC=',INDIC,' C=',TEXT12
           CALL XABORT('GEODMI: 1. INTEGER DATA EXPECTED')
          ENDIF
         ENDIF
         MIX(IHEX+IX,1,IZ1) = NITMA
         MAXMIX=MAX(MAXMIX,NITMA)
         MINMIX=MIN(MINMIX,NITMA)
 33      CONTINUE
         IHEX= IHEX+6*NC
        ENDIF
      ELSE
        CALL XABORT('GEODMI: MIXTURE # EXPECTED '//
     >              '(INTEGER EXPECTED).')
      ENDIF
      IF( IHEX.LT.LX )THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.3)THEN
            WRITE(6,*) ' TEST IZ1-2-3',IZ1,IZ2,IZ3,' IHEX NC',IHEX,NC
            CALL XABORT('GEODMI: KEYWORD *CROWN* OR *UPTO*'//
     >                  ' MUST BE READ.')
         ENDIF
         IF( TEXT12.EQ.'CROWN') THEN
            NCSAME= 1
         ELSEIF( TEXT12.EQ.'UPTO') THEN
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.NE.1) CALL XABORT('GEODMI: INTEGER DATA'//
     >                    ' EXPECTED AFTER *UPTO* KEYWORD')
            NCSAME= NITMA-NC-1
         ELSE
            CALL XABORT('GEODMI: KEYWORD *CROWN* OR *UPTO*'//
     >                  ' MUST BE READ.')
         ENDIF
         GO TO 30
      ELSEIF( IHEX.EQ.LX )THEN
         GO TO 25
      ELSE
         CALL XABORT('GEODMI: INVALID # OF MIX IN THIS PLANE.')
      ENDIF
*----
*  READ MIXTURE INDICES BY PLANE
*----
 20   IF (IZ2.EQ.0) THEN
       DO 22 IY=1,LY
        DO 21 IX=1,LX
         IF(.NOT.((IX.EQ.1).AND.(IY.EQ.1))) THEN
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          IF(INDIC.NE.1) CALL XABORT('GEODMI: 2. INTEGER DATA EXPECTED')
         ENDIF
         MIX(IX,IY,IZ1) = NITMA
         MAXMIX=MAX(MAXMIX,NITMA)
         MINMIX=MIN(MINMIX,NITMA)
 21     CONTINUE
 22    CONTINUE
      ELSE
       DO 24 IY=1,LY
        DO 23 IX=1,LX
         MIX(IX,IY,IZ1) = MIX(IX,IY,IZ2)
 23     CONTINUE
 24    CONTINUE
      ENDIF
*
 25   CONTINUE
      IF (NZ.LT.LZ) THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.3.OR.TEXT12.NE.'PLANE') THEN
            CALL XABORT('GEODMI: KEYWORD *PLANE* MUST BE READ.')
         ENDIF
         NC= -1
         IHEX= 0
         GO TO 5
      ENDIF
      IF (NZ.NE.LZ) CALL XABORT('GEODMI: WRONG NUMBER OF PLANES')
*
      RETURN
      END
