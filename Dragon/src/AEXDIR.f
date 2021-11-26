*DECK AEXDIR
      SUBROUTINE AEXDIR (NFICH,LBLOC,DATA,IADRES,LGSEG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read infomation from a direct access file. Component of a FORTRAN-77
* emulator of the SAPHYR archive system.
*
*Copyright:
* Copyright (C) 1999 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NFICH   unit number of the direct access file.
* LBLOC   direct access buffer length.
* IADRES  offset, from start of file where data is extracted from
*         or where data is to be stored.
* LGSEG   number of words to read from or write into file.
*
*Parameters: output
* DATA    address in memory where data is to be moved or extracted.
*
*-----------------------------------------------------------------------
*
      IMPLICIT INTEGER(A-Z)
      INTEGER DATA(LGSEG),LNEWAD(2)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: WRK
*
      ALLOCATE(WRK(LBLOC))
      INDEX=IADRES
      ID=0
      NROLD=0
   10 NREC=1+INDEX/LBLOC
      N=MOD(INDEX,LBLOC)
      LMIN=1
   20 IF(NREC.NE.NROLD) THEN
*       --------------------------------------------------------
        READ(NFICH,REC=NREC,ERR=90,IOSTAT=IR) (WRK(I),I=1,LBLOC)
*       --------------------------------------------------------
        NROLD=NREC
      ENDIF
      NGRO=MIN(LBLOC+LMIN-N-1,2)
      DO 30 L=LMIN,NGRO
      N=N+1
      LNEWAD(L)=WRK(N)
   30 CONTINUE
      IF(NGRO.EQ.2) GO TO 40
      NREC=NREC+1
      N=0
      LMIN=NGRO+1
      GO TO 20
   40 LINFO=LNEWAD(2)
      IF(ID+LINFO.GT.LGSEG) CALL XABORT('AEXDIR: DIRECT ACCESS READ FA'
     1 //'ILURE(1).')
      NREC=1+(INDEX+2)/LBLOC
      N=MOD(INDEX+2,LBLOC)
      LMIN=1
   50 IF(NREC.NE.NROLD) THEN
*       --------------------------------------------------------
        READ(NFICH,REC=NREC,ERR=90,IOSTAT=IR) (WRK(I),I=1,LBLOC)
*       --------------------------------------------------------
        NROLD=NREC
      ENDIF
      NGRO=MIN(LBLOC+LMIN-N-1,LINFO)
      DO 60 L=LMIN,NGRO
      N=N+1
      DATA(ID+L)=WRK(N)
   60 CONTINUE
      IF(NGRO.EQ.LINFO) GO TO 70
      NREC=NREC+1
      N=0
      LMIN=NGRO+1
      GO TO 50
*
   70 INDEX=LNEWAD(1)
      ID=ID+LNEWAD(2)
      IF(ID.EQ.LGSEG) GO TO 80
      GO TO 10
   80 DEALLOCATE(WRK)
      IF(LNEWAD(1).NE.-1) CALL XABORT('AEXDIR: DIRECT ACCESS READ FAIL'
     1 //'URE(3).')
      RETURN
   90 CALL XABORT('AEXDIR: DIRECT ACCESS READ FAILURE(2).')
      END
