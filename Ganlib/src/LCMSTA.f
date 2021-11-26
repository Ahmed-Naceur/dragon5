*DECK LCMSTA
      SUBROUTINE LCMSTA(IPLIS1,IPLIS2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compare the floating point information contained in the active
* directories of two tables or XSM files pointed by IPLIS1 and IPLIS2.
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
      TYPE(C_PTR) KDATA1(MAXLEV),KDATA2(MAXLEV)
      CHARACTER NAMT*12,CTMP1*4,CTMP2*4,HNAME1*12,HNAME2*12,NAMMY*12,
     1 HSMG*131,PATH(MAXLEV)*12,FIRST(MAXLEV)*12,MYDIR(MAXLEV)*12
      INTEGER IVEC(MAXLEV),KJLON(MAXLEV),IGO(MAXLEV)
      LOGICAL EMPTY,LCM
      TYPE(C_PTR) :: PT_DATA1,PT_DATA2
      INTEGER, POINTER :: III1(:),III2(:)
      REAL, POINTER :: RRR1(:),RRR2(:)
      LOGICAL, POINTER :: LLL1(:),LLL2(:)
      DOUBLE PRECISION, POINTER :: DDD1(:),DDD2(:)
      COMPLEX, POINTER :: CCC1(:),CCC2(:)
*
      CALL LCMVAL(IPLIS1,' ')
      CALL LCMVAL(IPLIS2,' ')
      ILEV=1
      KDATA1(1)=IPLIS1
      KDATA2(1)=IPLIS2
      KJLON(1)=-1
      IVEC(1)=1
      IGO(1)=5
      WRITE(6,'(/39H LCMSTA: COMPARISON OF TWO LCM OBJECTS.)')
*
* ASSOCIATIVE TABLE.
   10 CALL LCMINF(IPLIS1,HNAME1,NAMMY,EMPTY,ILONG,LCM)
      CALL LCMINF(IPLIS2,HNAME2,NAMMY,EMPTY,ILONG,LCM)
      MYDIR(ILEV)=NAMMY
      IF(EMPTY) GO TO (185,185,370,370,380),IGO(ILEV)
      NAMT=' '
      CALL LCMNXT(IPLIS1,NAMT)
*
      FIRST(ILEV)=NAMT
   15 CALL LCMLEN(IPLIS1,NAMT,ILON1,ITY1)
      CALL LCMLEN(IPLIS2,NAMT,ILON2,ITY2)
      IF((ILON1.NE.ILON2).OR.(ITY1.NE.ITY2)) THEN
         WRITE(6,'(/13H TWO BLOCKS '',A12,6H'' OF '',A12,7H'' AND '',
     1   A12,23H'' ARE OF UNEQUAL TYPE (,2I4,13H) OR LENGTH (,2I7,
     2   14H). DIRECTORY='',A12,2H''.)') NAMT,HNAME1,HNAME2,ITY1,
     3   ITY2,ILON1,ILON2,MYDIR(ILEV)
         GO TO 180
      ENDIF
      IF(ITY1.EQ.0) THEN
*        ASSOCIATIVE TABLE DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LCMSTA: TOO MANY DIRECTORY ',
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
            WRITE(HSMG,'(2A,A12,A)') 'LCMSTA: TOO MANY DIRECTORY ',
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
               WRITE(6,'(/40H INCONSISTENT INTEGER DATA ON THE TWO DI,
     1         19HRECTORIES. RECORD='',A12,13H'' DIRECTORY='',A12,
     2         1H'')') NAMT,MYDIR(ILEV)
               GO TO 180
            ENDIF
   80       CONTINUE
         ELSE IF((ITY1.EQ.2).OR.(ITY1.EQ.6)) THEN
*           COMPARE THE TWO SINGLE PRECISION OR COMPLEX BLOCKS.
            CALL LCMGPD(IPLIS1,NAMT,PT_DATA1)
            CALL LCMGPD(IPLIS2,NAMT,PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, RRR1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, RRR2, (/ ILON2 /))
            EPSMAX=0.0
            EPSAVG=0.0
            INGRO=0
            WRITE(6,'(/32H COMPARE REAL OR COMPLEX BLOCK '',A12,
     1      26H'' IN TABLES OR XSM FILES '',A12,7H'' AND '',A12,
     2      16H'' IN DIRECTORY '',A12,2H'':)') NAMT,HNAME1,HNAME2,
     3      MYDIR(ILEV)
            DO 100 I=1,ILON1
            ABSEP=ABS(RRR1(I)-RRR2(I))
            IF(EPSMAX.LT.ABSEP) THEN
               EPSMAX=ABSEP
               INGRO=I
            ENDIF
            EPSAVG=EPSAVG+ABSEP
  100       CONTINUE
            EPSAVG=EPSAVG/REAL(ILON1)
            WRITE (6,'(/5H LEN=,I6,5X,7HEPSMAX=,E12.5,13H IN COMPONENT,
     1      I6/16X,7HEPSAVG=,E12.5)') ILON1,EPSMAX,INGRO,EPSAVG
         ELSE IF(ITY1.EQ.3) THEN
*           CHARACTER*4 DATA.
            CALL LCMGPD(IPLIS1,NAMT,PT_DATA1)
            CALL LCMGPD(IPLIS2,NAMT,PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, III1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, III2, (/ ILON2 /))
            DO 130 I=1,ILON1
              WRITE(CTMP1,'(A4)') III1(I)
              WRITE(CTMP2,'(A4)') III2(I)
              IF(CTMP1.NE.CTMP2) THEN
                WRITE(6,'(/39H INCONSISTENT CHARACTER DATA ON THE TWO,
     1          22H DIRECTORIES. RECORD='',A12,13H'' DIRECTORY '',
     2          A12,8H'' DATA='',A4,3H'' '',A4,1H'')') NAMT,MYDIR(ILEV),
     3          CTMP1,CTMP2
                GO TO 180
              ENDIF
  130       CONTINUE
         ELSE IF(ITY1.EQ.4) THEN
*           COMPARE THE TWO DOUBLE PRECISION BLOCKS.
            CALL LCMGPD(IPLIS1,NAMT,PT_DATA1)
            CALL LCMGPD(IPLIS2,NAMT,PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, DDD1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, DDD2, (/ ILON2 /))
            EPSMAX=0.0
            EPSAVG=0.0
            INGRO=0
            WRITE(6,'(/33H COMPARE DOUBLE PRECISION BLOCK '',A12,
     1      26H'' IN TABLES OR XSM FILES '',A12,7H'' AND '',A12,
     2      16H'' IN DIRECTORY '',A12,2H'':)') NAMT,HNAME1,HNAME2,
     3      MYDIR(ILEV)
            DO 150 I=1,ILON1
            ABSEP=REAL(ABS(DDD1(I)-DDD2(I)))
            IF(EPSMAX.LT.ABSEP) THEN
               EPSMAX=ABSEP
               INGRO=I
            ENDIF
            EPSAVG=EPSAVG+ABSEP
  150       CONTINUE
            EPSAVG=EPSAVG/REAL(ILON1)
            WRITE (6,'(/5H LEN=,I6,5X,7HEPSMAX=,E12.5,13H IN COMPONENT,
     1      I6/16X,7HEPSAVG=,E12.5)') ILON1,EPSMAX,INGRO,EPSAVG
         ELSE IF(ITY1.EQ.5) THEN
*           LOGICAL DATA.
            CALL LCMGPD(IPLIS1,NAMT,PT_DATA1)
            CALL LCMGPD(IPLIS2,NAMT,PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, LLL1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, LLL2, (/ ILON2 /))
            DO 160 I=1,ILON1
            IF(LLL1(I).NEQV.LLL2(I)) THEN
               WRITE(6,'(/40H INCONSISTENT LOGICAL DATA ON THE TWO DI,
     1         19HRECTORIES. RECORD='',A12,13H'' DIRECTORY='',A12,
     2         1H'')') NAMT,MYDIR(ILEV)
               GO TO 180
            ENDIF
  160       CONTINUE
         ELSE IF(ITY1.EQ.6) THEN
*           COMPARE THE TWO COMPLEX BLOCKS.
            CALL LCMGPD(IPLIS1,NAMT,PT_DATA1)
            CALL LCMGPD(IPLIS2,NAMT,PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, CCC1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, CCC2, (/ ILON2 /))
            EPSMAX=0.0
            EPSAVG=0.0
            INGRO=0
            WRITE(6,'(/32H COMPARE REAL OR COMPLEX BLOCK '',A12,
     1      26H'' IN TABLES OR XSM FILES '',A12,7H'' AND '',A12,
     2      16H'' IN DIRECTORY '',A12,2H'':)') NAMT,HNAME1,HNAME2,
     3      MYDIR(ILEV)
            DO 170 I=1,ILON1
            ABSEP=ABS(CCC1(I)-CCC2(I))
            IF(EPSMAX.LT.ABSEP) THEN
               EPSMAX=ABSEP
               INGRO=I
            ENDIF
            EPSAVG=EPSAVG+ABSEP
  170       CONTINUE
            EPSAVG=EPSAVG/REAL(ILON1)
            WRITE (6,'(/5H LEN=,I6,5X,7HEPSMAX=,E12.5,13H IN COMPONENT,
     1      I6/16X,7HEPSAVG=,E12.5)') ILON1,EPSMAX,INGRO,EPSAVG
         ELSE
            CALL XABORT('LCMSTA: INVALID DATA TYPE(1).')
         ENDIF
      ENDIF
  180 CALL LCMNXT(IPLIS1,NAMT)
      IF(NAMT.NE.FIRST(ILEV)) GO TO 15
      GO TO (185,185,370,370,380),IGO(ILEV)
*
  185 NAMT=PATH(ILEV)
      ILEV=ILEV-1
      IPLIS1=KDATA1(ILEV)
      IPLIS2=KDATA2(ILEV)
      CALL LCMNXT(IPLIS1,NAMT)
      IF(NAMT.NE.FIRST(ILEV)) GO TO 15
      GO TO (185,185,370,370,380),IGO(ILEV)
*
* LIST.
  190 IVEC(ILEV)=IVEC(ILEV)+1
      IF(IVEC(ILEV).GT.KJLON(ILEV)) THEN
         GO TO (185,185,370,370,380),IGO(ILEV)
      ENDIF
      CALL LCMLEL(KDATA1(ILEV),IVEC(ILEV),ILON1,ITY1)
      CALL LCMLEL(KDATA2(ILEV),IVEC(ILEV),ILON2,ITY2)
      IF((ILON1.NE.ILON2).OR.(ITY1.NE.ITY2)) THEN
         WRITE(6,'(/15H TWO COMPONENTS,I6,5H OF '',A12,7H'' AND '',
     1   A12,23H'' ARE OF UNEQUAL TYPE (,2I4,13H) OR LENGTH (,
     2   2I7,2H).)') IVEC(ILEV),HNAME1,HNAME2,ITY1,ITY2,ILON1,
     3   ILON2
         GO TO 190
      ENDIF
      IF((ILON1.NE.0).AND.(ITY1.EQ.0)) THEN
*        ASSOCIATIVE TABLE DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LCMSTA: TOO MANY DIRECTORY ',
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
            WRITE(HSMG,'(2A,A12,A)') 'LCMSTA: TOO MANY DIRECTORY ',
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
            DO 230 I=1,ILON1
            IF(III1(I).NE.III2(I)) THEN
               WRITE(6,'(/40H INCONSISTENT INTEGER DATA ON THE TWO DI,
     1         24HRECTORIES. LIST ELEMENT=,I5,1H.)') IVEC(ILEV)
               GO TO 190
            ENDIF
  230       CONTINUE
         ELSE IF((ITY1.EQ.2).OR.(ITY1.EQ.6)) THEN
*           COMPARE THE TWO SINGLE PRECISION BLOCKS.
            CALL LCMGPL(IPLIS1,IVEC(ILEV),PT_DATA1)
            CALL LCMGPL(IPLIS2,IVEC(ILEV),PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, RRR1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, RRR2, (/ ILON2 /))
            EPSMAX=0.0
            EPSAVG=0.0
            INGRO=0
            WRITE(6,'(/37H COMPARE REAL OR COMPLEX LIST ELEMENT,I5,
     1      25H IN TABLES OR XSM FILES '',A12,7H'' AND '',A12,2H'':)
     2      ') IVEC(ILEV),HNAME1,HNAME2
            DO 250 I=1,ILON1
            ABSEP=ABS(RRR1(I)-RRR2(I))
            IF(EPSMAX.LT.ABSEP) THEN
               EPSMAX=ABSEP
               INGRO=I
            ENDIF
            EPSAVG=EPSAVG+ABSEP
  250       CONTINUE
            EPSAVG=EPSAVG/REAL(ILON1)
            WRITE (6,'(/5H LEN=,I6,5X,7HEPSMAX=,E12.5,13H IN COMPONENT,
     1      I6/16X,7HEPSAVG=,E12.5)') ILON1,EPSMAX,INGRO,EPSAVG
         ELSE IF(ITY1.EQ.3) THEN
*           CHARACTER*4 DATA.
            CALL LCMGPL(IPLIS1,IVEC(ILEV),PT_DATA1)
            CALL LCMGPL(IPLIS2,IVEC(ILEV),PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, III1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, III2, (/ ILON2 /))
            DO 280 I=1,ILON1
              WRITE(CTMP1,'(A4)') III1(I)
              WRITE(CTMP2,'(A4)') III2(I)
              IF(CTMP1.NE.CTMP2) THEN
                WRITE(6,'(/40H INCONSISTENT CHARACTER DATA ON THE TWO ,
     1          26HDIRECTORIES. LIST ELEMENT=,I5,8H'' DATA='',A4,
     2          3H'' '',A4,2H''.)') IVEC(ILEV),CTMP1,CTMP2
                GO TO 190
              ENDIF
  280       CONTINUE
         ELSE IF(ITY1.EQ.4) THEN
*           COMPARE THE TWO DOUBLE PRECISION BLOCKS.
            CALL LCMGPL(IPLIS1,IVEC(ILEV),PT_DATA1)
            CALL LCMGPL(IPLIS2,IVEC(ILEV),PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, DDD1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, DDD2, (/ ILON2 /))
            EPSMAX=0.0
            EPSAVG=0.0
            INGRO=0
            WRITE(6,'(/33H COMPARE DOUBLE PRECISION BLOCK '',A12,
     1      26H'' IN TABLES OR XSM FILES '',A12,7H'' AND '',A12,
     2      2H'':)') NAMT,HNAME1,HNAME2
            DO 300 I=1,ILON1
            ABSEP=REAL(ABS(DDD1(I)-DDD2(I)))
            IF(EPSMAX.LT.ABSEP) THEN
               EPSMAX=ABSEP
               INGRO=I
            ENDIF
            EPSAVG=EPSAVG+ABSEP
  300       CONTINUE
            EPSAVG=EPSAVG/REAL(ILON1)
            WRITE (6,'(/5H LEN=,I6,5X,7HEPSMAX=,E12.5,13H IN COMPONENT,
     1      I6/16X,7HEPSAVG=,E12.5)') ILON1,EPSMAX,INGRO,EPSAVG
         ELSE IF(ITY1.EQ.5) THEN
*           LOGICAL DATA.
            CALL LCMGPL(IPLIS1,IVEC(ILEV),PT_DATA1)
            CALL LCMGPL(IPLIS2,IVEC(ILEV),PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, LLL1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, LLL2, (/ ILON2 /))
            DO 340 I=1,ILON1
            IF(LLL1(I).NEQV.LLL2(I)) THEN
               WRITE(6,'(/40H INCONSISTENT LOGICAL DATA ON THE TWO DI,
     1         24HRECTORIES. LIST ELEMENT=,I5,1H.)') IVEC(ILEV)
               GO TO 190
            ENDIF
  340       CONTINUE
         ELSE IF((ITY1.EQ.2).OR.(ITY1.EQ.6)) THEN
*           COMPARE THE TWO COMPLEX BLOCKS.
            CALL LCMGPL(IPLIS1,IVEC(ILEV),PT_DATA1)
            CALL LCMGPL(IPLIS2,IVEC(ILEV),PT_DATA2)
            CALL C_F_POINTER(PT_DATA1, CCC1, (/ ILON1 /))
            CALL C_F_POINTER(PT_DATA2, CCC2, (/ ILON2 /))
            EPSMAX=0.0
            EPSAVG=0.0
            INGRO=0
            WRITE(6,'(/37H COMPARE REAL OR COMPLEX LIST ELEMENT,I5,
     1      25H IN TABLES OR XSM FILES '',A12,7H'' AND '',A12,2H'':)
     2      ') IVEC(ILEV),HNAME1,HNAME2
            DO 350 I=1,ILON1
            ABSEP=ABS(CCC1(I)-CCC2(I))
            IF(EPSMAX.LT.ABSEP) THEN
               EPSMAX=ABSEP
               INGRO=I
            ENDIF
            EPSAVG=EPSAVG+ABSEP
  350       CONTINUE
            EPSAVG=EPSAVG/REAL(ILON1)
            WRITE (6,'(/5H LEN=,I6,5X,7HEPSMAX=,E12.5,13H IN COMPONENT,
     1      I6/16X,7HEPSAVG=,E12.5)') ILON1,EPSMAX,INGRO,EPSAVG
         ELSE
            CALL XABORT('LCMSTA: INVALID DATA TYPE(2).')
         ENDIF
      ENDIF
      GO TO 190
*
  370 ILEV=ILEV-1
      IPLIS1=KDATA1(ILEV)
      IPLIS2=KDATA2(ILEV)
      GO TO 190
*
  380 RETURN
      END
