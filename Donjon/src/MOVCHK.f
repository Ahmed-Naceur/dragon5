*DECK MOVCHK
      SUBROUTINE MOVCHK(IMPX,IMODE,NPART,IAXIS,ITOP,DELH,LENG,RODPOS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute and set the new rod position and insertion level for the
* fading or moving rod.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki and A. Hebert
*
*Parameters: input
* IMPX    printing index (=0 for no print).
* IMODE   type of displacement: =1 for FADE; =2 for MOVE (DONJON3-type
*         movement).
* NPART   number of parts in the control rod.
* IAXIS   axis of rod movement: =1 for X; =2 for Y; =3 for Z.
* ITOP    rod insertion: = +1 from the top; = -1 from the bottom.
* DELH    rod displacement along the IAXIS of movement if FADE; position
*         of moving end in core if MOVE.
* LENG    fully-inserted complete rod position.
* RODPOS  fully inserted rod position.
*
*Parameters: output
* RODPOS  new rod position.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,IMODE,NPART,IAXIS,ITOP
      REAL DELH,RODPOS(6,NPART),LENG(2),LIMINF
*
      PARAMETER(IOUT=6)
*----
*  X-AXIS MOVEMENT
*----
      IF(IMPX.GT.1) WRITE(IOUT,1000) DELH
      SUP=DELH
      DO 10 IPART=1,NPART
      IF(IAXIS.EQ.1) THEN
        IF((ITOP.EQ.1).AND.(IMODE.EQ.1)) THEN
          RODPOS(1,IPART)=MAX(RODPOS(1,IPART),LENG(2)-DELH)
          RODPOS(1,IPART)=MIN(RODPOS(1,IPART),RODPOS(2,IPART))
        ELSE IF((ITOP.EQ.-1).AND.(IMODE.EQ.1)) THEN
          RODPOS(2,IPART)=MIN(RODPOS(2,IPART),LENG(1)+DELH)
          RODPOS(2,IPART)=MAX(RODPOS(1,IPART),RODPOS(2,IPART))
        ELSE IF((ITOP.EQ.1).AND.(IMODE.EQ.2)) THEN
          DELTA=RODPOS(2,IPART)-RODPOS(1,IPART)
          RODPOS(1,IPART)=SUP
          RODPOS(2,IPART)=SUP+DELTA
          SUP=RODPOS(2,IPART)
        ELSE IF((ITOP.EQ.-1).AND.(IMODE.EQ.2)) THEN
          DELTA=RODPOS(2,IPART)-RODPOS(1,IPART)
          RODPOS(2,IPART)=SUP
          RODPOS(1,IPART)=SUP-DELTA
          SUP=RODPOS(1,IPART)
        ENDIF
*----
*  Y-AXIS MOVEMENT
*----
      ELSE IF(IAXIS.EQ.2) THEN
        IF((ITOP.EQ.1).AND.(IMODE.EQ.1)) THEN
          RODPOS(3,IPART)=MAX(RODPOS(3,IPART),LENG(2)-DELH)
          RODPOS(3,IPART)=MIN(RODPOS(3,IPART),RODPOS(4,IPART))
        ELSE IF((ITOP.EQ.-1).AND.(IMODE.EQ.1)) THEN
          RODPOS(4,IPART)=MIN(RODPOS(4,IPART),LENG(1)+DELH)
          RODPOS(4,IPART)=MAX(RODPOS(3,IPART),RODPOS(4,IPART))
        ELSE IF((ITOP.EQ.1).AND.(IMODE.EQ.2)) THEN
          DELTA=RODPOS(4,IPART)-RODPOS(3,IPART)
          RODPOS(3,IPART)=SUP
          RODPOS(4,IPART)=SUP+DELTA
          SUP=RODPOS(4,IPART)
        ELSE IF((ITOP.EQ.-1).AND.(IMODE.EQ.2)) THEN
          DELTA=RODPOS(4,IPART)-RODPOS(3,IPART)
          RODPOS(4,IPART)=SUP
          RODPOS(3,IPART)=SUP-DELTA
          SUP=RODPOS(3,IPART)
        ENDIF
*----
*  Z-AXIS MOVEMENT
*----
      ELSE IF(IAXIS.EQ.3) THEN
        IF((ITOP.EQ.1).AND.(IMODE.EQ.1)) THEN
          RODPOS(5,IPART)=MAX(RODPOS(5,IPART),LENG(2)-DELH)
          RODPOS(5,IPART)=MIN(RODPOS(5,IPART),RODPOS(6,IPART))
        ELSE IF((ITOP.EQ.-1).AND.(IMODE.EQ.1)) THEN
          RODPOS(6,IPART)=MIN(RODPOS(6,IPART),LENG(1)+DELH)
          RODPOS(6,IPART)=MAX(RODPOS(5,IPART),RODPOS(6,IPART))
        ELSE IF((ITOP.EQ.1).AND.(IMODE.EQ.2)) THEN
          DELTA=RODPOS(6,IPART)-RODPOS(5,IPART)
          RODPOS(5,IPART)=SUP
          RODPOS(6,IPART)=SUP+DELTA
          SUP=RODPOS(6,IPART)
        ELSE IF((ITOP.EQ.-1).AND.(IMODE.EQ.2)) THEN
          DELTA=RODPOS(6,IPART)-RODPOS(5,IPART)
          RODPOS(6,IPART)=SUP
          RODPOS(5,IPART)=SUP-DELTA
          SUP=RODPOS(5,IPART)
        ENDIF
      ENDIF
*----
*  PRINT NEW POSITION
*----
      IF(IMPX.GT.1) THEN
        WRITE(IOUT,1001) IPART,RODPOS(1,IPART),RODPOS(3,IPART),
     1         RODPOS(5,IPART),RODPOS(2,IPART),RODPOS(4,IPART),
     2         RODPOS(6,IPART)
        ENDIF
   10 CONTINUE
*----
*  CONSISTENCY CHECK
*----
      LIMINF=0
      IF(IMODE.EQ.2) THEN
        IF(ITOP.EQ.-1) THEN
          LIMINF=DELH-(LENG(2)-LENG(1))
        ELSE IF(ITOP.EQ.1) THEN
          LIMINF=DELH+(LENG(2)-LENG(1))
        ENDIF
        IF(ABS(SUP-LIMINF).GT.1.E-3) CALL XABORT('@MOVCHK: WRONG LENG'
     1  //'TH OF ADJUSTER')
      ENDIF
      RETURN
*
 1000 FORMAT(/5X,'MOVCHK: MOVE A ROD BY',F10.4)
 1001 FORMAT(
     1 /5X,'MOVCHK: PART =',I5/
     2  5X,'NEW ROD POSITION :'/
     3  5X,'X-',F10.4,5X,'Y-',F10.4,5X,'Z-',F10.4/
     4  5X,'X+',F10.4,5X,'Y+',F10.4,5X,'Z+',F10.4)
      END
