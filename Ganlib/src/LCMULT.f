*DECK LCMULT
      SUBROUTINE LCMULT(IPLIST,FLOTT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multiply the floating point information contained in the active
* directory of a table or XSM file by a real number.
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
* FLOTT   real number.
*
*Parameters: output
* IPLIS1  address of the table or handle to the XSM file.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIST
      REAL FLOTT
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXLEV=50)
      TYPE(C_PTR) KDATA(MAXLEV)
      CHARACTER NAMT*12,HSMG*131,MYNAME*12,NAMLCM*12,PATH(MAXLEV)*12,
     1 FIRST(MAXLEV)*12
      LOGICAL EMPTY,LCM
      INTEGER IVEC(MAXLEV),KJLON(MAXLEV),IGO(MAXLEV)
      TYPE(C_PTR) :: PT_DATA
      REAL, POINTER :: RRR(:)
      DOUBLE PRECISION, POINTER :: DDD(:)
      COMPLEX, POINTER :: CCC(:)
*
      CALL LCMVAL(IPLIST,' ')
      ILEV=1
      KDATA(1)=IPLIST
      KJLON(1)=-1
      IVEC(1)=1
      IGO(1)=5
*
* ASSOCIATIVE TABLE.
   10 CALL LCMINF(IPLIST,MYNAME,NAMLCM,EMPTY,ILONG,LCM)
      IF(EMPTY) GO TO (100,100,240,240,250),IGO(ILEV)
      NAMT=' '
      CALL LCMNXT(IPLIST,NAMT)
*
      FIRST(ILEV)=NAMT
   15 CALL LCMLEN(IPLIST,NAMT,ILON1,ITY1)
      IF(ITY1.EQ.0) THEN
*        ASSOCIATIVE TABLE DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LCMULT: TOO MANY DIRECTORY ',
     1      'LEVELS ON ''',MYNAME,'''(1).'
            CALL XABORT(HSMG)
         ENDIF
         KJLON(ILEV)=-1
         KDATA(ILEV)=LCMGID(IPLIST,NAMT)
         PATH(ILEV)=NAMT
         IPLIST=KDATA(ILEV)
         IVEC(ILEV)=1
         IGO(ILEV)=1
         GO TO 10
      ELSE IF(ITY1.EQ.10) THEN
*        LIST DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LCMULT: TOO MANY DIRECTORY ',
     1      'LEVELS ON ''',MYNAME,'''(2).'
            CALL XABORT(HSMG)
         ENDIF
         KJLON(ILEV)=ILON1
         KDATA(ILEV)=LCMGID(IPLIST,NAMT)
         PATH(ILEV)=NAMT
         IPLIST=KDATA(ILEV)
         IVEC(ILEV)=0
         IGO(ILEV)=2
         GO TO 190
      ELSE IF(ITY1.EQ.2) THEN
*        SINGLE PRECISION DATA.
         CALL LCMGPD(IPLIST,NAMT,PT_DATA)
         CALL C_F_POINTER(PT_DATA, RRR, (/ ILON1 /))
         DO 70 I=1,ILON1
         RRR(I)=FLOTT*RRR(I)
   70    CONTINUE
         CALL LCMPPD(IPLIST,NAMT,ILON1,ITY1,PT_DATA)
      ELSE IF(ITY1.EQ.4) THEN
*        DOUBLE PRECISION DATA.
         CALL LCMGPD(IPLIST,NAMT,PT_DATA)
         CALL C_F_POINTER(PT_DATA, DDD, (/ ILON1 /))
         DO 80 I=1,ILON1
         DDD(I)=FLOTT*DDD(I)
   80    CONTINUE
         CALL LCMPPD(IPLIST,NAMT,ILON1,ITY1,PT_DATA)
      ELSE IF(ITY1.EQ.6) THEN
*        COMPLEX DATA.
         CALL LCMGPD(IPLIST,NAMT,PT_DATA)
         CALL C_F_POINTER(PT_DATA, CCC, (/ ILON1 /))
         DO 90 I=1,ILON1
         CCC(I)=FLOTT*CCC(I)
   90    CONTINUE
         CALL LCMPPD(IPLIST,NAMT,ILON1,ITY1,PT_DATA)
      ENDIF
      CALL LCMNXT(IPLIST,NAMT)
      IF(NAMT.NE.FIRST(ILEV)) GO TO 15
      GO TO (100,100,240,240,250),IGO(ILEV)
*
  100 NAMT=PATH(ILEV)
      ILEV=ILEV-1
      IPLIST=KDATA(ILEV)
      CALL LCMNXT(IPLIST,NAMT)
      IF(NAMT.NE.FIRST(ILEV)) GO TO 15
      GO TO (100,100,240,240,250),IGO(ILEV)
*
* LIST.
  190 IVEC(ILEV)=IVEC(ILEV)+1
      IF(IVEC(ILEV).GT.KJLON(ILEV)) THEN
         GO TO (100,100,240,240,250),IGO(ILEV)
      ENDIF
      CALL LCMLEL(KDATA(ILEV),IVEC(ILEV),ILON1,ITY1)
      IF((ILON1.NE.0).AND.(ITY1.EQ.0)) THEN
*        ASSOCIATIVE TABLE DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LCMULT: TOO MANY DIRECTORY ',
     1      'LEVELS ON ''',MYNAME,'''(3).'
            CALL XABORT(HSMG)
         ENDIF
         KJLON(ILEV)=-1
         KDATA(ILEV)=LCMGIL(IPLIST,IVEC(ILEV-1))
         IPLIST=KDATA(ILEV)
         IVEC(ILEV)=1
         IGO(ILEV)=3
         GO TO 10
      ELSE IF((ILON1.NE.0).AND.(ITY1.EQ.10)) THEN
*        LIST DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LCMULT: TOO MANY DIRECTORY ',
     1      'LEVELS ON ''',MYNAME,'''(4).'
            CALL XABORT(HSMG)
         ENDIF
         KJLON(ILEV)=ILON1
         KDATA(ILEV)=LCMGIL(IPLIST,IVEC(ILEV-1))
         IPLIST=KDATA(ILEV)
         IVEC(ILEV)=0
         IGO(ILEV)=4
         GO TO 190
      ELSE IF((ILON1.NE.0).AND.(ITY1.LE.6)) THEN
         IF((ITY1.EQ.2).OR.(ITY1.EQ.6)) THEN
*           SINGLE PRECISION DATA.
            CALL LCMGPL(IPLIST,IVEC(ILEV),PT_DATA)
            CALL C_F_POINTER(PT_DATA, RRR, (/ ILON1 /))
            DO 210 I=1,ILON1
            RRR(I)=FLOTT*RRR(I)
  210       CONTINUE
            CALL LCMPPL(IPLIST,IVEC(ILEV),ILON1,ITY1,PT_DATA)
         ELSE IF(ITY1.EQ.4) THEN
*           DOUBLE PRECISION DATA.
            CALL LCMGPL(IPLIST,IVEC(ILEV),PT_DATA)
            CALL C_F_POINTER(PT_DATA, DDD, (/ ILON1 /))
            DO 220 I=1,ILON1
            DDD(I)=FLOTT*DDD(I)
  220       CONTINUE
            CALL LCMPPL(IPLIST,IVEC(ILEV),ILON1,ITY1,PT_DATA)
         ELSE IF(ITY1.EQ.6) THEN
*           COMPLEX DATA.
            CALL LCMGPL(IPLIST,IVEC(ILEV),PT_DATA)
            CALL C_F_POINTER(PT_DATA, CCC, (/ ILON1 /))
            DO 230 I=1,ILON1
            CCC(I)=FLOTT*CCC(I)
  230       CONTINUE
            CALL LCMPPL(IPLIST,IVEC(ILEV),ILON1,ITY1,PT_DATA)
         ENDIF
      ENDIF
      GO TO 190
*
  240 ILEV=ILEV-1
      IPLIST=KDATA(ILEV)
      GO TO 190
*
  250 RETURN
      END
