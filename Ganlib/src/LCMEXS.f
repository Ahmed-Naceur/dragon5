*DECK LCMEXS
      SUBROUTINE LCMEXS(IPLIST,IMPX,NUNIT,IMODE,IDIR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Export/import the content of a table or xsm file using the contour
* method. Export start from the active directory. This version is
* backward compatible with the Saphyr version of xsm file export
* format.
*
*Copyright:
* Copyright (C) 1993 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIS1  address of the table or handle to the XSM file.
* IMPX    equal to zero for no print.
* NUNIT   file unit number where the export/import is performed.
* IMODE   type of export/import file: =1 sequential unformatted;
*         =2 sequential formatted (ascii).
* IDIR    type of operation: =1 to export; =2 to import.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIST
      INTEGER IMPX,NUNIT,IMODE,IDIR
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NLEVEL=50)
      CHARACTER NAMT*12,MYNAME*12,PATH(NLEVEL)*12,FIRST(NLEVEL)*12,
     1 NAMLCM*12,HSMG*131,CMEDIU(2)*8
      LOGICAL EMPTY,LCM
      TYPE(C_PTR) PT_DATA
      DATA (CMEDIU(II),II=1,2)/'TABLE','XSM FILE'/
*
      CALL LCMINF(IPLIST,NAMLCM,MYNAME,EMPTY,ILONG,LCM)
      IMED=2
      IF(LCM) IMED=1
      IF(ILONG.GE.0) THEN
         WRITE(HSMG,'(46HLCMEXS: UNABLE TO IMPORT/EXPORT A LIST IN THE ,
     1   A8,8H NAMED '',A12,2H''.)') CMEDIU(IMED),NAMLCM
         CALL XABORT(HSMG)
      ENDIF
      IF((IMODE.LT.1).OR.(IMODE.GT.2)) THEN
         WRITE(HSMG,'(33HLCMEXS: INVALID FILE TYPE ON THE ,A8,
     1   8H NAMED '',A12,2H''.)') CMEDIU(IMED),NAMLCM
         CALL XABORT(HSMG)
      ENDIF
      ITOT=0
      ILEVEL=1
      IF(IDIR.EQ.1) THEN
         IF(IMPX.GT.0)THEN
            WRITE(6,300) 'EXPORT',CMEDIU(IMED),NAMLCM,MYNAME
         ENDIF
         CALL LCMVAL(IPLIST,' ')
         GO TO 10
      ELSE IF(IDIR.EQ.2) THEN
         IF(IMPX.GT.0)THEN
            WRITE(6,300) 'IMPORT',CMEDIU(IMED),NAMLCM,MYNAME
         ENDIF
         GO TO 50
      ELSE IF(EMPTY) THEN
         WRITE(HSMG,'(14HLCMEXS: EMPTY ,A8,8H NAMED '',A12,2H''.)')
     1   CMEDIU(IMED),NAMLCM
         CALL XABORT(HSMG)
      ELSE
         WRITE(HSMG,'(30HLCMEXS: INVALID ACTION ON THE ,A8,8H NAMED '',
     1   A12,2H''.)') CMEDIU(IMED),NAMLCM
         CALL XABORT(HSMG)
      ENDIF
*----
*  FILE EXPORT.
*----
   10 NAMT=' '
      LENNAM=12
      CALL LCMNXT(IPLIST,NAMT)
      IF(NAMT.EQ.' ') THEN
         IF(ILEVEL.EQ.1) RETURN
         NAMT=PATH(ILEVEL)
         ILEVEL=ILEVEL-1
         CALL LCMSIX(IPLIST,' ',2)
         IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
            WRITE(NUNIT) 0,0,0,0
         ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
            WRITE(NUNIT,310) 0,0,0,0
         ENDIF
         IF(IMPX.GT.0) WRITE(6,350) ILEVEL
         GO TO 30
      ENDIF
      FIRST(ILEVEL)=NAMT
*
   20 CALL LCMLEN(IPLIST,NAMT,ILONG,ITYLCM)
      IF(ITYLCM.EQ.0) ILONG=1
      IF(IMPX.GT.0) WRITE(6,320) ILEVEL,NAMT,ITYLCM,ILONG
      IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
         WRITE(NUNIT) ILEVEL,LENNAM,ITYLCM,ILONG
         WRITE(NUNIT) NAMT
      ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
         WRITE(NUNIT,310) ILEVEL,LENNAM,ITYLCM,ILONG
         WRITE(NUNIT,'(A12,68(1H ))') NAMT
      ENDIF
      IF(ITYLCM.EQ.0) THEN
*        DIRECTORY DATA.
         ILEVEL=ILEVEL+1
         IF(ILEVEL.GT.NLEVEL) CALL XABORT('LCMEXS: TOO MANY DIRECTORY '
     1   //'LEVELS.')
         CALL LCMSIX(IPLIST,NAMT,1)
         PATH(ILEVEL)=NAMT
         GO TO 10
      ELSE IF((ILONG.NE.0).AND.(ITYLCM.LE.6)) THEN
         ITOT=ITOT+ILONG
         IF(NUNIT.NE.0) THEN
            CALL LCMGPD(IPLIST,NAMT,PT_DATA)
*           ------------------ EXPORT A NODE -----------------
            CALL LCMNOS(NUNIT,IMODE,IDIR,ILONG,ITYLCM,PT_DATA)
*           --------------------------------------------------
         ENDIF
      ENDIF
   30 CALL LCMNXT(IPLIST,NAMT)
      IF(NAMT.EQ.FIRST(ILEVEL)) THEN
         IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
            WRITE(NUNIT) 0,0,0,0
         ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
            WRITE(NUNIT,310) 0,0,0,0
         ENDIF
         IF(IMPX.GT.0) WRITE(6,350) ILEVEL
         IF(ILEVEL.EQ.1) GO TO 40
         NAMT=PATH(ILEVEL)
         ILEVEL=ILEVEL-1
         CALL LCMSIX(IPLIST,' ',2)
         GO TO 30
      ENDIF
      GO TO 20
   40 IF(IMPX.GT.0) WRITE(6,330) 'EXPORTED',ITOT
      RETURN
*----
*  FILE IMPORT.
*----
   50 IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
         READ(NUNIT,END=80) JLEVEL,LENNAM,ITYLCM,ILONG
         IF(LENNAM.GT.12) THEN
            CALL XABORT('LCMEXS: A RECORD NAME IS GREATER THAN 12 CHAR'
     1      //'ACTERS(1).')
         ENDIF
         READ(NUNIT) NAMT
      ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
         READ(NUNIT,340,END=70) JLEVEL,LENNAM,ITYLCM,ILONG
         IF(LENNAM.GT.12) THEN
            CALL XABORT('LCMEXS: A RECORD NAME IS GREATER THAN 12 CHAR'
     1      //'ACTERS(2).')
         ENDIF
         READ(NUNIT,'(A12)') NAMT
      ENDIF
      IF(JLEVEL.NE.1) THEN
         WRITE(HSMG,'(29HLCMEXS: UNABLE TO IMPORT THE ,A8,9H LOCATED ,
     1   7HON UNIT,I3,1H.)') CMEDIU(IMED),NUNIT
         CALL XABORT(HSMG)
      ENDIF
*
   60 IF(ITYLCM.EQ.0) THEN
*        DIRECTORY DATA.
         IF(IMPX.GT.0) WRITE(6,320) JLEVEL,NAMT,ITYLCM
         ILEVEL=ILEVEL+1
         CALL LCMSIX(IPLIST,NAMT,1)
      ELSE
         IF(IMPX.GT.0) WRITE(6,320) JLEVEL,NAMT,ITYLCM,ILONG
         JLONG=ILONG
         IF((ITYLCM.EQ.4).OR.(ITYLCM.EQ.6)) JLONG=2*ILONG
         PT_DATA = LCMARA(JLONG)
*        ----------------- IMPORT A NODE ------------------
         CALL LCMNOS(NUNIT,IMODE,IDIR,ILONG,ITYLCM,PT_DATA)
*        --------------------------------------------------
         CALL LCMPPD(IPLIST,NAMT,ILONG,ITYLCM,PT_DATA)
         ITOT=ITOT+ILONG
      ENDIF
   70 IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
         READ(NUNIT,END=70) JLEVEL,LENNAM,ITYLCM,ILONG
      ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
         READ(NUNIT,340,END=70) JLEVEL,LENNAM,ITYLCM,ILONG
      ENDIF
      IF(JLEVEL.EQ.0) THEN
         IF(IMPX.GT.0) WRITE(6,350) ILEVEL
         ILEVEL=ILEVEL-1
         IF(ILEVEL.EQ.0) GO TO 80
         CALL LCMSIX(IPLIST,' ',2)
         GO TO 70
      ELSE
         IF(JLEVEL.NE.ILEVEL) THEN
            CALL XABORT('LCMEXS: IMPORT FAILURE.')
         ELSE IF(LENNAM.GT.12) THEN
            CALL XABORT('LCMEXS: A RECORD NAME IS GREATER THAN 12 CHAR'
     1      //'ACTERS(3).')
         ENDIF
         IF((NUNIT.NE.0).AND.(IMODE.EQ.1)) THEN
            READ(NUNIT) NAMT
         ELSE IF((NUNIT.NE.0).AND.(IMODE.EQ.2)) THEN
            READ(NUNIT,'(A12)') NAMT
         ENDIF
         GO TO 60
      ENDIF
*
   80 IF(IMPX.GT.0) WRITE(6,330) 'IMPORTED',ITOT
      RETURN
*
  300 FORMAT (//9H LCMEXS: ,A6,1H ,A8,8H NAMED ',A12,15H' FROM ACTIVE D,
     1 10HIRECTORY ',A12,3H' ://18H LEVEL  BLOCK NAME,4(1H-),4X,5HTYPE ,
     2 7H LENGTH/)
  310 FORMAT ('->',4I8,46(1H ))
  320 FORMAT ('&*',I5,'  ''',A12,'''',2I8)
  330 FORMAT (/23H TOTAL NUMBER OF WORDS ,A8,2H =,I10/)
  340 FORMAT (2X,4I8)
  350 FORMAT ('&*',I5,2X,14('-'))
      END
