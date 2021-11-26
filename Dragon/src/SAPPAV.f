*DECK SAPPAV
      SUBROUTINE SAPPAV(IPSAP,IPAR,NPAR,TYPE,RVAL,IVAL,CVAL,IV,LGNEW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To return the index of a global parameter value. Reorganize the
* 'paramvaleur' directory if required.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPSAP   pointer to the Saphyb.
* IPAR    index of the global parameter.
* NPAR    total number of global parameters.
* TYPE    type of the global parameter value.
* RVAL    global parameter value if TYPE='FLOTTANT'.
* IVAL    global parameter value if TYPE='ENTIER'.
* CVAL    global parameter value if TYPE='CHAINE'.
*
*Parameters: output
* IV      index of the global parameter value.
* LGNEW   new parameter flag (=.true. if the parameter value is new).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSAP
      INTEGER IPAR,NPAR,IV,IVAL
      REAL RVAL
      LOGICAL LGNEW
      CHARACTER TYPE*8,CVAL*(*)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (REPS=1.0E-5,MAXPAR=50,MAXVAL=1000)
      CHARACTER RECNAM*12,VCHAR(MAXVAL)*12
      INTEGER NVALUE(MAXPAR)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: VINTE
      REAL, ALLOCATABLE, DIMENSION(:) :: VREAL
*
      CALL LCMSIX(IPSAP,'paramdescrip',1)
      CALL LCMGET(IPSAP,'NVALUE',NVALUE)
      CALL LCMSIX(IPSAP,' ',2)
      CALL LCMSIX(IPSAP,'paramvaleurs',1)
      WRITE(RECNAM,'(''pval'',I8)') IPAR
*
      LGNEW=.TRUE.
      IF(TYPE.EQ.'FLOTTANT') THEN
         ALLOCATE(VREAL(NVALUE(IPAR)+1))
         IF(NVALUE(IPAR).EQ.0) THEN
            IV=1
            VREAL(IV)=RVAL
            NVALUE(IPAR)=1
         ELSE
            CALL LCMLEN(IPSAP,RECNAM,ILONG,ITYLCM)
            IF(ILONG.GT.NVALUE(IPAR)) CALL XABORT('SAPPAV: NVALUE OVER'
     1      //'FLOW(1).')
            CALL LCMGET(IPSAP,RECNAM,VREAL)
            DO 10 I=1,NVALUE(IPAR)
               IF(RVAL.LE.VREAL(I)*(1.+REPS))THEN
                  IV=I
                  LGNEW=RVAL.LT.VREAL(IV)*(1.-REPS)
                  GO TO 20
               ENDIF
   10       CONTINUE
            IV=NVALUE(IPAR)+1
   20       IF(LGNEW) THEN
               NVALUE(IPAR)=NVALUE(IPAR)+1
               DO 30 J=NVALUE(IPAR)-1,IV,-1
                  VREAL(J+1)=VREAL(J)
   30          CONTINUE
               VREAL(IV)=RVAL
            ENDIF
         ENDIF
         IF(LGNEW) CALL LCMPUT(IPSAP,RECNAM,NVALUE(IPAR),2,VREAL)
         DEALLOCATE(VREAL)
      ELSE IF(TYPE.EQ.'ENTIER') THEN
         ALLOCATE(VINTE(NVALUE(IPAR)+1))
         IF(NVALUE(IPAR).EQ.0) THEN
            IV=1
            VINTE(IV)=IVAL
            NVALUE(IPAR)=1
         ELSE
            CALL LCMLEN(IPSAP,RECNAM,ILONG,ITYLCM)
            IF(ILONG.GT.NVALUE(IPAR)) CALL XABORT('SAPPAV: NVALUE OVER'
     1      //'FLOW(2).')
            CALL LCMGET(IPSAP,RECNAM,VINTE)
            DO 40 I=1,NVALUE(IPAR)
               IF(IVAL.LE.VINTE(I))THEN
                  IV=I
                  LGNEW=IVAL.LT.VINTE(IV)
                  GO TO 50
               ENDIF
   40       CONTINUE
            IV=NVALUE(IPAR)+1
   50       IF(LGNEW) THEN
               NVALUE(IPAR)=NVALUE(IPAR)+1
               DO 60 J=NVALUE(IPAR)-1,IV,-1
                  VINTE(J+1)=VINTE(J)
   60          CONTINUE
               VINTE(IV)=IVAL
            ENDIF
         ENDIF
         IF(LGNEW) CALL LCMPUT(IPSAP,RECNAM,NVALUE(IPAR),1,VINTE)
         DEALLOCATE(VINTE)
      ELSE IF(TYPE.EQ.'CHAINE') THEN
         IF(NVALUE(IPAR).EQ.0) THEN
            IV=1
            VCHAR(IV)=CVAL
            NVALUE(IPAR)=1
         ELSE
            CALL LCMGTC(IPSAP,RECNAM,12,NVALUE(IPAR),VCHAR)
            DO 70 I=1,NVALUE(IPAR)
               IF(CVAL.EQ.VCHAR(I))THEN
                  IV=I
                  LGNEW=.FALSE.
                  GO TO 80
               ENDIF
   70       CONTINUE
            IV=NVALUE(IPAR)+1
   80       IF(LGNEW) THEN
               NVALUE(IPAR)=NVALUE(IPAR)+1
               IF(NVALUE(IPAR).GT.MAXVAL) CALL XABORT('SAPPAV: MAXVAL '
     1         //'OVERFLOW.')
               VCHAR(IV)=CVAL
            ENDIF
         ENDIF
         IF(LGNEW) CALL LCMPTC(IPSAP,RECNAM,12,NVALUE(IPAR),VCHAR)
      ENDIF
      CALL LCMSIX(IPSAP,' ',2)
*
      IF(LGNEW) THEN
         CALL LCMSIX(IPSAP,'paramdescrip',1)
         CALL LCMPUT(IPSAP,'NVALUE',NPAR,1,NVALUE)
         CALL LCMSIX(IPSAP,' ',2)
      ENDIF
      RETURN
      END
