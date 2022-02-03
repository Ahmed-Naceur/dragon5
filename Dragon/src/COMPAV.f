*DECK COMPAV
      SUBROUTINE COMPAV(IPCPO,IPAR,NPAR,TYPE,RVAL,IVAL,CVAL,IV,LGNEW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Return the index of a global or local parameter value. Reorganize the
* 'GLOBAL' or 'LOCAL' directory if required.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPCPO   pointer to the 'GLOBAL' or 'LOCAL' directory of the
*         multicompo.
* IPAR    index of the parameter.
* NPAR    total number of parameters.
* TYPE    type of the parameter value.
* RVAL    parameter value if TYPE='REAL'.
* IVAL    parameter value if TYPE='INTEGER'.
* CVAL    parameter value if TYPE='STRING'.
*
*Parameters: output
* IV      index of the global or local parameter value.
* LGNEW   new parameter flag (=.true. if the parameter value is new).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPCPO
      INTEGER IPAR,NPAR,IV,IVAL
      REAL RVAL
      LOGICAL LGNEW
      CHARACTER TYPE*8,CVAL*(*)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (REPS=1.0E-4,MAXPAR=50,MAXVAL=1000)
      CHARACTER RECNAM*12,VCHAR(MAXVAL)*12
      INTEGER NVALUE(MAXPAR)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: VINTE
      REAL, ALLOCATABLE, DIMENSION(:) :: VREAL
*
      CALL LCMLEN(IPCPO,'NVALUE',ILONG,ITYLCM)
      IF(ILONG.EQ.NPAR) THEN
         CALL LCMGET(IPCPO,'NVALUE',NVALUE)
      ELSE
         CALL XDISET(NVALUE,NPAR,0)
      ENDIF
      WRITE(RECNAM,'(''pval'',I8.8)') IPAR
*
      LGNEW=.TRUE.
      IF(TYPE.EQ.'REAL') THEN
         ALLOCATE(VREAL(NVALUE(IPAR)+1))
         IF(NVALUE(IPAR).EQ.0) THEN
            IV=1
            VREAL(IV)=RVAL
            NVALUE(IPAR)=1
         ELSE
            CALL LCMLEN(IPCPO,RECNAM,ILONG,ITYLCM)
            IF(ILONG.GT.NVALUE(IPAR)) CALL XABORT('COMPAV: NVALUE OVER'
     1      //'FLOW(1).')
            CALL LCMGET(IPCPO,RECNAM,VREAL)
            DO 20 I=1,NVALUE(IPAR)
               IF(RVAL.LE.VREAL(I)*(1.+REPS))THEN
                  IV=I
                  LGNEW=RVAL.LT.VREAL(IV)*(1.-REPS)
                  GO TO 30
               ENDIF
   20       CONTINUE
            IV=NVALUE(IPAR)+1
   30       IF(LGNEW) THEN
               NVALUE(IPAR)=NVALUE(IPAR)+1
               DO 40 J=NVALUE(IPAR)-1,IV,-1
                  VREAL(J+1)=VREAL(J)
   40          CONTINUE
               VREAL(IV)=RVAL
            ENDIF
         ENDIF
         IF(LGNEW) CALL LCMPUT(IPCPO,RECNAM,NVALUE(IPAR),2,VREAL)
         DEALLOCATE(VREAL)
      ELSE IF(TYPE.EQ.'INTEGER') THEN
         ALLOCATE(VINTE(NVALUE(IPAR)+1))
         IF(NVALUE(IPAR).EQ.0) THEN
            IV=1
            VINTE(IV)=IVAL
            NVALUE(IPAR)=1
         ELSE
            CALL LCMLEN(IPCPO,RECNAM,ILONG,ITYLCM)
            IF(ILONG.GT.NVALUE(IPAR)) CALL XABORT('COMPAV: NVALUE OVER'
     1      //'FLOW(2).')
            CALL LCMGET(IPCPO,RECNAM,VINTE)
            DO 50 I=1,NVALUE(IPAR)
               IF(IVAL.LE.VINTE(I))THEN
                  IV=I
                  LGNEW=IVAL.LT.VINTE(IV)
                  GO TO 60
               ENDIF
   50       CONTINUE
            IV=NVALUE(IPAR)+1
   60       IF(LGNEW) THEN
               NVALUE(IPAR)=NVALUE(IPAR)+1
               DO 70 J=NVALUE(IPAR)-1,IV,-1
                  VINTE(J+1)=VINTE(J)
   70          CONTINUE
               VINTE(IV)=IVAL
            ENDIF
         ENDIF
         IF(LGNEW) CALL LCMPUT(IPCPO,RECNAM,NVALUE(IPAR),1,VINTE)
         DEALLOCATE(VINTE)
      ELSE IF(TYPE.EQ.'STRING') THEN
         IF(NVALUE(IPAR).EQ.0) THEN
            IV=1
            VCHAR(IV)=CVAL
            NVALUE(IPAR)=1
         ELSE
            CALL LCMGTC(IPCPO,RECNAM,12,NVALUE(IPAR),VCHAR)
            DO 80 I=1,NVALUE(IPAR)
               IF(CVAL.EQ.VCHAR(I))THEN
                  IV=I
                  LGNEW=.FALSE.
                  GO TO 90
               ENDIF
   80       CONTINUE
            IV=NVALUE(IPAR)+1
   90       IF(LGNEW) THEN
               NVALUE(IPAR)=NVALUE(IPAR)+1
               IF(NVALUE(IPAR).GT.MAXVAL) CALL XABORT('COMPAV: MAXVAL '
     1         //'OVERFLOW.')
               VCHAR(IV)=CVAL
            ENDIF
         ENDIF
         IF(LGNEW) CALL LCMPTC(IPCPO,RECNAM,12,NVALUE(IPAR),VCHAR)
      ENDIF
*
      IF(LGNEW) CALL LCMPUT(IPCPO,'NVALUE',NPAR,1,NVALUE)
      RETURN
      END
