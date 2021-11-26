*DECK LCMADD
      SUBROUTINE LCMADD(IPLIS1,IPLIS2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Add the floating point information contained in the active directories
* of two tables or XSM files pointed by IPLIS1 and IPLIS2 and store the
* result in the table or XSM file pointed by IPLIS2.
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
* IPLIS2  address of the table or handle to the XSM file.
*
*Parameters: output
* IPLIS2  address of the table or handle to the XSM file.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIS1,IPLIS2
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXLEV=50)
      CHARACTER NAMT*12,HSMG*131,CTMP1*4,CTMP2*4,HNAME1*12,HNAME2*12,
     1 NAMMY*12,PATH(MAXLEV)*12,FIRST(MAXLEV)*12
      TYPE(C_PTR) KDATA1(MAXLEV),KDATA2(MAXLEV)
      INTEGER IVEC(MAXLEV),KJLON(MAXLEV),IGO(MAXLEV)
      LOGICAL EMPTY,LCM1,LCM2
      TYPE(C_PTR) :: PT_DATA1,PT_DATA2
      INTEGER, POINTER :: III1(:),III2(:)
      REAL, POINTER :: RRR1(:),RRR2(:)
      LOGICAL, POINTER :: LLL1(:),LLL2(:)
      DOUBLE PRECISION, POINTER :: DDD1(:),DDD2(:)
      COMPLEX, POINTER :: CCC1(:),CCC2(:)
*
      IF(C_ASSOCIATED(IPLIS1,IPLIS2)) THEN
         WRITE(HSMG,'(45HLCMADD: TWO TABLES OR XSM FILES HAVE THE SAME,
     1   8H HANDLE.)')
         CALL XABORT(HSMG)
      ENDIF
      CALL LCMVAL(IPLIS1,' ')
      CALL LCMVAL(IPLIS2,' ')
      ILEV=1
      KDATA1(1)=IPLIS1
      KDATA2(1)=IPLIS2
      KJLON(1)=-1
      IVEC(1)=1
      IGO(1)=5
*
* ASSOCIATIVE TABLE.
   10 CALL LCMINF(IPLIS1,HNAME1,NAMMY,EMPTY,ILONG,LCM1)
      CALL LCMINF(IPLIS2,HNAME2,NAMMY,EMPTY,ILONG,LCM2)
      IF(EMPTY) GO TO (150,150,270,270,380),IGO(ILEV)
      NAMT=' '
      CALL LCMNXT(IPLIS1,NAMT)
*
      FIRST(ILEV)=NAMT
   15 CALL LCMLEN(IPLIS1,NAMT,ILON1,ITY1)
      CALL LCMLEN(IPLIS2,NAMT,ILON2,ITY2)
      IF((ILON1.NE.ILON2).OR.(ITY1.NE.ITY2)) THEN
         WRITE(6,'(/21H LCMADD: TWO BLOCKS '',A12,6H'' OF '',A12,
     1   7H'' AND '',A12,23H'' ARE OF UNEQUAL TYPE (,2I4,8H) OR LEN,
     2   5HGTH (,2I7,2H).)') NAMT,HNAME1,HNAME2,ITY1,ITY2,ILON1,ILON2
         GO TO 10
      ENDIF
      IF(ITY1.EQ.0) THEN
*        ASSOCIATIVE TABLE DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LCMADD: TOO MANY DIRECTORY ',
     1      'LEVELS ON ''',HNAME2,'''(1).'
            CALL XABORT(HSMG)
         ENDIF
         KJLON(ILEV)=-1
         KDATA1(ILEV)=LCMGID(IPLIS1,NAMT)
         KDATA2(ILEV)=LCMGID(IPLIS2,NAMT)
         PATH(ILEV)=NAMT
         IPLIS1=KDATA1(ILEV)
         IPLIS2=KDATA2(ILEV)
         IVEC(ILEV)=1
         IGO(ILEV)=1
         GO TO 10
      ELSE IF(ITY1.EQ.10) THEN
*        LIST DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LCMADD: TOO MANY DIRECTORY ',
     1      'LEVELS ON ''',HNAME2,'''(2).'
            CALL XABORT(HSMG)
         ENDIF
         KJLON(ILEV)=ILON1
         KDATA1(ILEV)=LCMGID(IPLIS1,NAMT)
         KDATA2(ILEV)=LCMGID(IPLIS2,NAMT)
         PATH(ILEV)=NAMT
         IPLIS1=KDATA1(ILEV)
         IPLIS2=KDATA2(ILEV)
         IVEC(ILEV)=0
         IGO(ILEV)=2
         GO TO 190
      ELSE IF(ITY1.LE.6) THEN
         IF(ITY1.EQ.1) THEN
*           INTEGER DATA.
            CALL LCMGPD(IPLIS1,NAMT,PT_DATA1)
            CALL LCMGPD(IPLIS2,NAMT,PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, III1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, III2, (/ ILON2 /))
            DO 80 I=1,ILON1
            IF(III1(I).NE.III2(I)) THEN
               WRITE(HSMG,'(39HLCMADD: INCONSISTENT INTEGER DATA ON TH,
     1         27HE TWO DIRECTORIES. RECORD='',A12,1H'')') NAMT
               CALL XABORT(HSMG)
            ENDIF
   80       CONTINUE
         ELSE IF(ITY1.EQ.2) THEN
*           SINGLE PRECISION DATA.
            CALL LCMGPD(IPLIS1,NAMT,PT_DATA1)
            CALL LCMGPD(IPLIS2,NAMT,PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, RRR1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, RRR2, (/ ILON2 /))
            DO 90 I=1,ILON1
            RRR2(I)=RRR1(I)+RRR2(I)
   90       CONTINUE
            CALL LCMPPD(IPLIS2,NAMT,ILON2,ITY2,PT_DATA2)
         ELSE IF(ITY1.EQ.3) THEN
*           CHARACTER*4 DATA.
            CALL LCMGPD(IPLIS1,NAMT,PT_DATA1)
            CALL LCMGPD(IPLIS2,NAMT,PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, III1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, III2, (/ ILON2 /))
            DO 100 I=1,ILON1
              WRITE(CTMP1,'(A4)') III1(I)
              WRITE(CTMP2,'(A4)') III2(I)
              IF(CTMP1.NE.CTMP2) THEN
                 WRITE(HSMG,'(37HLCMADD: INCONSISTENT CHARACTER DATA O,
     1           31HN THE TWO DIRECTORIES. RECORD='',A12,1H'')') NAMT
                 CALL XABORT(HSMG)
              ENDIF
  100       CONTINUE
         ELSE IF(ITY1.EQ.4) THEN
*           DOUBLE PRECISION DATA.
            CALL LCMGPD(IPLIS1,NAMT,PT_DATA1)
            CALL LCMGPD(IPLIS2,NAMT,PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, DDD1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, DDD2, (/ ILON2 /))
            DO 110 I=1,ILON1
            DDD2(I)=DDD1(I)+DDD2(I)
  110       CONTINUE
            CALL LCMPPD(IPLIS2,NAMT,ILON2,ITY2,PT_DATA2)
         ELSE IF(ITY1.EQ.5) THEN
*           LOGICAL DATA.
            CALL LCMGPD(IPLIS1,NAMT,PT_DATA1)
            CALL LCMGPD(IPLIS2,NAMT,PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, LLL1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, LLL2, (/ ILON2 /))
            DO 120 I=1,ILON1
            IF(LLL1(I).NEQV.LLL2(I)) THEN
               WRITE(HSMG,'(39HLCMADD: INCONSISTENT LOGICAL DATA ON TH,
     1         27HE TWO DIRECTORIES. RECORD='',A12,1H'')') NAMT
               CALL XABORT(HSMG)
            ENDIF
  120       CONTINUE
         ELSE IF(ITY1.EQ.6) THEN
*           COMPLEX DATA.
            CALL LCMGPD(IPLIS1,NAMT,PT_DATA1)
            CALL LCMGPD(IPLIS2,NAMT,PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, CCC1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, CCC2, (/ ILON2 /))
            DO 130 I=1,ILON1
            CCC2(I)=CCC1(I)+CCC2(I)
  130       CONTINUE
            CALL LCMPPD(IPLIS2,NAMT,ILON2,ITY2,PT_DATA2)
         ELSE
            CALL XABORT('LCMADD: INVALID DATA TYPE(1).')
         ENDIF
      ENDIF
      CALL LCMNXT(IPLIS1,NAMT)
      IF(NAMT.NE.FIRST(ILEV)) GO TO 15
      GO TO (150,150,270,270,380),IGO(ILEV)
*
  150 NAMT=PATH(ILEV)
      ILEV=ILEV-1
      IPLIS1=KDATA1(ILEV)
      IPLIS2=KDATA2(ILEV)
      CALL LCMNXT(IPLIS1,NAMT)
      IF(NAMT.NE.FIRST(ILEV)) GO TO 15
      GO TO (150,150,270,270,380),IGO(ILEV)
*
* LIST.
  190 IVEC(ILEV)=IVEC(ILEV)+1
      IF(IVEC(ILEV).GT.KJLON(ILEV)) THEN
         GO TO (150,150,270,270,380),IGO(ILEV)
      ENDIF
      CALL LCMLEL(KDATA1(ILEV),IVEC(ILEV),ILON1,ITY1)
      CALL LCMLEL(KDATA2(ILEV),IVEC(ILEV),ILON2,ITY2)
      IF((ILON1.NE.ILON2).OR.(ITY1.NE.ITY2)) THEN
         WRITE(6,'(/24H LCMADD: TWO COMPONENTS ,I6,5H OF '',A12,
     1   7H'' AND '',A12,23H'' ARE OF UNEQUAL TYPE (,2I4,8H) OR LEN,
     2   5HGTH (,2I7,2H).)') IVEC(ILEV),HNAME1,HNAME2,ITY1,ITY2,ILON1,
     3   ILON2
         GO TO 190
      ENDIF
      IF((ILON1.NE.0).AND.(ITY1.EQ.0)) THEN
*        ASSOCIATIVE TABLE DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LCMADD: TOO MANY DIRECTORY ',
     1      'LEVELS ON ''',HNAME2,'''(3).'
            CALL XABORT(HSMG)
         ENDIF
         KJLON(ILEV)=-1
         KDATA1(ILEV)=LCMGIL(IPLIS1,IVEC(ILEV-1))
         KDATA2(ILEV)=LCMGIL(IPLIS2,IVEC(ILEV-1))
         IPLIS1=KDATA1(ILEV)
         IPLIS2=KDATA2(ILEV)
         IVEC(ILEV)=1
         IGO(ILEV)=3
         GO TO 10
      ELSE IF((ILON1.NE.0).AND.(ITY1.EQ.10)) THEN
*        LIST DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LCMADD: TOO MANY DIRECTORY ',
     1      'LEVELS ON ''',HNAME2,'''(4).'
            CALL XABORT(HSMG)
         ENDIF
         KJLON(ILEV)=ILON1
         KDATA1(ILEV)=LCMGIL(IPLIS1,IVEC(ILEV-1))
         KDATA2(ILEV)=LCMGIL(IPLIS2,IVEC(ILEV-1))
         IPLIS1=KDATA1(ILEV)
         IPLIS2=KDATA2(ILEV)
         IVEC(ILEV)=0
         IGO(ILEV)=4
         GO TO 190
      ELSE IF((ILON1.NE.0).AND.(ITY1.LE.6)) THEN
         IF(ITY1.EQ.1) THEN
*           INTEGER DATA.
            CALL LCMGPL(IPLIS1,IVEC(ILEV),PT_DATA1)
            CALL LCMGPL(IPLIS2,IVEC(ILEV),PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, III1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, III2, (/ ILON2 /))
            DO 220 I=1,ILON1
            IF(III1(I).NE.III2(I)) THEN
               WRITE(HSMG,'(39HLCMADD: INCONSISTENT INTEGER DATA ON TH,
     1         32HE TWO DIRECTORIES. LIST ELEMENT=,I5,1H.)') IVEC(ILEV)
               CALL XABORT(HSMG)
            ENDIF
  220       CONTINUE
          ELSE IF((ITY1.EQ.2).OR.(ITY1.EQ.6)) THEN
*           SINGLE PRECISION DATA.
            CALL LCMGPL(IPLIS1,IVEC(ILEV),PT_DATA1)
            CALL LCMGPL(IPLIS2,IVEC(ILEV),PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, RRR1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, RRR2, (/ ILON2 /))
            DO 230 I=1,ILON1
            RRR2(I)=RRR1(I)+RRR2(I)
  230       CONTINUE
            CALL LCMPPL(IPLIS2,IVEC(ILEV),ILON2,ITY2,PT_DATA2)
         ELSE IF(ITY1.EQ.3) THEN
*           CHARACTER*4 DATA.
            CALL LCMGPL(IPLIS1,IVEC(ILEV),PT_DATA1)
            CALL LCMGPL(IPLIS2,IVEC(ILEV),PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, III1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, III2, (/ ILON2 /))
            DO 240 I=1,ILON1
              WRITE(CTMP1,'(A4)') III1(I)
              WRITE(CTMP2,'(A4)') III2(I)
              IF(CTMP1.NE.CTMP2) THEN
                 WRITE(HSMG,'(38HLCMADD: INCONSISTENT CHARACTER DATA ON,
     1           35H THE TWO DIRECTORIES. LIST ELEMENT=,I5,1H.)')
     2           IVEC(ILEV)
                 CALL XABORT(HSMG)
              ENDIF
  240       CONTINUE
         ELSE IF(ITY1.EQ.4) THEN
*           DOUBLE PRECISION DATA.
            CALL LCMGPL(IPLIS1,IVEC(ILEV),PT_DATA1)
            CALL LCMGPL(IPLIS2,IVEC(ILEV),PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, DDD1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, DDD2, (/ ILON2 /))
            DO 250 I=1,ILON1
            DDD2(I)=DDD1(I)+DDD2(I)
  250       CONTINUE
            CALL LCMPPL(IPLIS2,IVEC(ILEV),ILON2,ITY2,PT_DATA2)
         ELSE IF(ITY1.EQ.5) THEN
*           LOGICAL DATA.
            CALL LCMGPL(IPLIS1,IVEC(ILEV),PT_DATA1)
            CALL LCMGPL(IPLIS2,IVEC(ILEV),PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, LLL1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, LLL2, (/ ILON2 /))
            DO 260 I=1,ILON1
            IF(LLL1(I).NEQV.LLL2(I)) THEN
               WRITE(HSMG,'(39HLCMADD: INCONSISTENT LOGICAL DATA ON TH,
     1         32HE TWO DIRECTORIES. LIST ELEMENT=,I5,1H.)') IVEC(ILEV)
               CALL XABORT(HSMG)
            ENDIF
  260       CONTINUE
          ELSE IF(ITY1.EQ.6) THEN
*           COMPLEX DATA.
            CALL LCMGPL(IPLIS1,IVEC(ILEV),PT_DATA1)
            CALL LCMGPL(IPLIS2,IVEC(ILEV),PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, CCC1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, CCC2, (/ ILON2 /))
            DO 265 I=1,ILON1
            CCC2(I)=CCC1(I)+CCC2(I)
  265       CONTINUE
            CALL LCMPPL(IPLIS2,IVEC(ILEV),ILON2,ITY2,PT_DATA2)
         ELSE
            CALL XABORT('LCMADD: INVALID DATA TYPE(2).')
         ENDIF
      ENDIF
      GO TO 190
*
  270 ILEV=ILEV-1
      IPLIS1=KDATA1(ILEV)
      IPLIS2=KDATA2(ILEV)
      GO TO 190
*
  380 RETURN
      END
