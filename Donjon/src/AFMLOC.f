*DECK AFMLOC
      SUBROUTINE AFMLOC(NBURN,NTP,XBMAX,XBMIN,XBURN,MAX,MIN,COF,ILIN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Burnup localisation and interpolation
*
*Copyright:
* Copyright (C) 1996 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* M.T. Sissaoui
*
*Parameters: input
*  NBURN  total number of burnup steps.
*  XBURN  burnup steps dimemsion (NBURN).
*  XBMAX  higher burnup value.
*  XBMIN  lower burnup value.
*
*Parameters: output
*  MAX    maximum burnup number
*  MIN    minimum burnup number
*  COF    interpolation coefficient (Lagrange)
*
*Parameters: 
* NTP     
* MAX     
* MIN     
* ILIN    
*
*---------------------------------------------------------------*
*
      DIMENSION XBURN(NBURN),ELMT(3)
      DOUBLE PRECISION COF(3),XCOF(1)
      NTP=2
      COF(1)=0.0D0
      COF(2)=0.0D0
      COF(3)=0.0D0
      IF(XBMAX.EQ.XBMIN) NTP=1
      IF(XBMAX.GT.XBURN(NBURN)) THEN
         WRITE(6,100) XBMAX,XBURN(NBURN)
         CALL XABORT('AFMLOC: THE HIGHER BURNUP VALUE IS BEYOND'
     1   //' THE MAXIMUM BURNUP IN THE DATABASE')
      ELSE IF(NBURN.EQ.1.AND.NTP.EQ.2) THEN
         CALL XABORT('AFMLOC: TIME AVERAGE CALCULATION REQUIRE'
     1    //' AT LEAST TWO IRRADIATIONS STEPS')
      ELSE IF(NBURN.EQ.1.AND.NTP.EQ.1) THEN
         COF(1)=1.0D0
         MIN=1
         MAX=1
      ELSE IF(NBURN.EQ.2) THEN
         MIN=1
         MAX=2
         IF(NTP.EQ.1) THEN
           XIRAD=XBMIN
           IF(ILIN.EQ.1) THEN
             NTOX=-1
           ELSE
             NTOX=2
           ENDIF
           NELE=2
           ELMT(1)=XBURN(1)
           ELMT(2)=XBURN(2)
           CALL LIBLEX(NELE,XIRAD,ELMT,NTOX,XCOF)
         ENDIF
      ELSE IF(NBURN.GE.3) THEN
         DO 85 IV=1,NTP
           IF(IV.EQ.1) THEN
             XIRAD=XBMIN
           ELSE
             XIRAD=XBMAX
           ENDIF
*
           DO 80 I=2,NBURN
             IF(XIRAD.GE.XBURN(I-1).AND.XIRAD.LE.XBURN(I)) THEN
               IF(NTP.EQ.2) THEN
                 IF(IV.EQ.1) THEN
                   MIN=I-1
                 ELSE
                   IF(I+1.LE.NBURN) THEN
                     MAX=I+1
                   ELSE
                    MAX=I
                   ENDIF
                 ENDIF
               ELSE
                 IF(I+1.LE.NBURN) THEN
                   MIN=I-1
                   MAX=I+1
                 ELSE
                   MIN=I-2
                   MAX=I
                 ENDIF
               ENDIF
             ENDIF
 80        CONTINUE
 85      CONTINUE
         IF(NTP.EQ.1) THEN
           IF(ILIN.EQ.1) THEN
             NTOX=-1
           ELSE
             NTOX=3
           ENDIF
           NELE=3
           ELMT(1)=XBURN(MAX-2)
           ELMT(2)=XBURN(MAX-1)
           ELMT(3)=XBURN(MAX)
           CALL LIBLEX(NELE,XIRAD,ELMT,NTOX,COF)
         ENDIF
      ENDIF
      RETURN
*
  100 FORMAT(/30H AFMLOC: MAXIMUM BURNUP VALUE=,1P,E12.4/
     1  9X,25HMAXIMUM TABULATED BURNUP=,E12.4)
      END
