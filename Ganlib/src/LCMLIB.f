*DECK LCMLIB
      SUBROUTINE LCMLIB(IPLIST)
*
*----------------------------------------------------------------------
*
*Purpose:
* List the LCM entries contained in a table or a XSM file.
*
*Copyright:
* Copyright (C) 1993 Ecole Polytechnique de Montreal
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIST  address of the table or handle to the XSM file.
*
*----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIST
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NTYPE=11)
      CHARACTER NAMT*12,NAMLCM*12,MYNAME*12,FIRST*12,CTYPE(NTYPE)*16,
     1 CMEDIU(2)*8
      LOGICAL EMPTY,LCM
      SAVE CTYPE,CMEDIU
      CHARACTER(LEN=12) HSIGN
      DATA (CTYPE(ITY),ITY=1,NTYPE)/'DIRECTORY','INTEGER','REAL',
     > 'CHARACTER','DOUBLE PRECISION','LOGICAL','COMPLEX','UNDEFINED',
     > ' ',' ','LIST'/
      DATA (CMEDIU(II),II=1,2)/'TABLE','XSM FILE'/
*
      CALL LCMINF(IPLIST,NAMLCM,MYNAME,EMPTY,ILONG,LCM)
      IMED=1
      IF(.NOT.LCM) IMED=2
      ITOT=0
      NAMT=' '
      IF(ILONG.EQ.-1) THEN
         IF(EMPTY) THEN
            WRITE (6,80) MYNAME,CMEDIU(IMED),NAMLCM
            RETURN
         ENDIF
         CALL LCMNXT(IPLIST,NAMT)
         FIRST=NAMT
         WRITE(6,100) MYNAME,CMEDIU(IMED),NAMLCM
         INMT=0
*
   10    INMT=INMT+1
         CALL LCMLEN(IPLIST,NAMT,ILONG,ITYLCM)
         IF((ITYLCM.EQ.0).OR.(ITYLCM.EQ.10)) THEN
            WRITE (6,120) INMT,NAMT,ILONG,CTYPE(ITYLCM+1)
         ELSE IF((ITYLCM.GE.1).AND.(ITYLCM.LE.6)) THEN
            IF((NAMT.EQ.'SIGNATURE').AND.(ITYLCM.EQ.3)) THEN
               CALL LCMGTC(IPLIST,NAMT,12,1,HSIGN)
               WRITE (6,110) INMT,NAMT,ILONG,CTYPE(ITYLCM+1),
     1         HSIGN
            ELSE
               WRITE (6,120) INMT,NAMT,ILONG,CTYPE(ITYLCM+1)
            ENDIF
            ITOT=ITOT+ILONG
         ELSE
            WRITE (6,120) INMT,NAMT,ILONG,CTYPE(8)
         ENDIF
         CALL LCMNXT(IPLIST,NAMT)
         IF(NAMT.EQ.FIRST) GO TO 20
         GO TO 10
*
   20    WRITE(6,130) MYNAME,ITOT
      ELSE
         IF(ILONG.EQ.0) THEN
            WRITE (6,90) MYNAME,CMEDIU(IMED),NAMLCM
            RETURN
         ENDIF
         WRITE(6,100) MYNAME,CMEDIU(IMED),NAMLCM
         DO 30 INMT=1,ILONG
         CALL LCMLEL(IPLIST,INMT,ILONG,ITYLCM)
         IF((ITYLCM.EQ.0).OR.(ITYLCM.EQ.10)) THEN
            WRITE (6,120) INMT,NAMT,ILONG,CTYPE(ITYLCM+1)
         ELSE IF((ITYLCM.GE.1).AND.(ITYLCM.LE.6)) THEN
            WRITE (6,120) INMT,NAMT,ILONG,CTYPE(ITYLCM+1)
            ITOT=ITOT+ILONG
         ELSE
            WRITE (6,120) INMT,NAMT,ILONG,CTYPE(8)
         ENDIF
   30    CONTINUE
         WRITE(6,140) MYNAME,ITOT
      ENDIF
      RETURN
*
   80 FORMAT (/10H LCMLIB: ',A12,31H' IS AN EMPTY DIRECTORY OF THE ,A8,
     1 2H ',A12,2H'.)
   90 FORMAT (/10H LCMLIB: ',A12,26H' IS AN EMPTY LIST OF THE ,A8,2H ',
     1 A12,2H'.)
  100 FORMAT (//38H LCMLIB: CONTENT OF ACTIVE DIRECTORY ',A12,
     1 9H' OF THE ,A8,2H ',A12,2H'://5X,10HBLOCK NAME,10(1H-),4X,
     2 6HLENGTH,4X,4HTYPE/)
  110 FORMAT (1X,I8,3H  ',A12,1H',I10,4X,A16,2H=',A12,1H')
  120 FORMAT (1X,I8,3H  ',A12,1H',I10,4X,A16)
  130 FORMAT (//37H TOTAL NUMBER OF WORDS ON DIRECTORY ',A12,3H' =,
     > I10/)
  140 FORMAT (//32H TOTAL NUMBER OF WORDS ON LIST ',A12,3H' =,I10/)
      END
