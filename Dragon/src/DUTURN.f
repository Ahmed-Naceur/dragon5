*DECK DUTURN
      SUBROUTINE DUTURN(IHEX,TURN,NCEL,TURND,NCELA,CELL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Provide orientation of cell in an assembly with symetry IHEX.
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): M. Ouisloumen
*
*Parameters: input
* IHEX    symmetry type.
* NCEL    number of cells in symmetric assembly.
* NCELA   number of cells in unfolded assembly.
* CELL    cell index in symmetric assembly.
* TURN    cell orientation in symmetric assembly.
*
*Parameters: output
* TURND   cell orientation in unfolded assembly.
*
*-----------------------------------------------------------------------
*
      INTEGER TAB(6),CELL(NCELA),TAB6(6),TAB9(6),TURN(NCEL),
     +                     TURND(NCELA),TAB12(6),TABR8(6),
     +                     TABA8(6),TABB8(6)
      LOGICAL LGR8,LGSA,LGSB,LGSA6
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NUM,NTURN,ITAB
      SAVE TAB,TAB6,TAB9,TAB12,TABR8,TABA8,TABB8
      DATA TAB,TAB6,TAB9,TAB12,TABR8,TABA8,TABB8
     +     /1,6,5,4,3,2,2,1,6,5,4,3,3,2,1,6,5,4,3,4,5,6,1,2,4,5,6,1,2,3
     +     ,3,2,1,6,5,4,6,5,4,3,2,1/
*
      IFONC(N,L)= 2+(N-1)*(L+3*(N-2))
      IFCOUR(N)=NINT( (4.+SQRT(1.+4.*FLOAT(N-1)/3.)
     +                 +SQRT(1.+4.*FLOAT(N-2)/3.))*.25)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NUM(NCEL),NTURN(NCEL),ITAB(NCEL))
*
      DO 10 I=2,NCEL
         IF(TURN(I).EQ.1.OR.TURN(I).EQ.9) THEN
           ITAB(I)=2
         ELSEIF(TURN(I).EQ.2.OR.TURN(I).EQ.10) THEN
           ITAB(I)=3
         ELSEIF(TURN(I).EQ.3.OR.TURN(I).EQ.11) THEN
           ITAB(I)=4
         ELSEIF(TURN(I).EQ.4.OR.TURN(I).EQ.12) THEN
           ITAB(I)=5
         ELSEIF(TURN(I).EQ.5.OR.TURN(I).EQ.7) THEN
           ITAB(I)=6
         ELSEIF(TURN(I).EQ.6.OR.TURN(I).EQ.8) THEN
           ITAB(I)=1
         ELSE
           CALL XABORT('DUTURN : INVALID ORIENTATION ')
         ENDIF
10    CONTINUE
*
      NCOUR=IFCOUR(NCELA)
      IF(IHEX.EQ.1) THEN
        GOTO 20
      ELSEIF(IHEX.EQ.2.OR.IHEX.EQ.3) THEN
        GOTO 40
      ELSEIF(IHEX.EQ.4) THEN
        GOTO 60
      ELSEIF(IHEX.EQ.5) THEN
        GOTO 80
      ELSEIF(IHEX.GT.5.AND.IHEX.LT.9) THEN
        GOTO 100
      ELSE
        CALL XABORT('DUTURN : INVALID TYPE OF GEOMETRY ')
      ENDIF
*
   20 CONTINUE
*
*       DUPLICATION DE L'ORIENTATION DANS LA SYMETRIE S30
*
      TURND(1)=TURN(1)
*
*       DUPLICATION DE LA 2EME COURONE
*
      ITURN=ITAB(2)
      TURND(2)=TURN(2)
      DO 25 I=3,7
         ITURN=ITURN+1
         IF(ITURN.GT.6)ITURN=ITURN-6
         IF(ITURN.EQ.1) THEN
           TURND(I)=8
           IF(TURN(CELL(I)).LE.6)TURND(I)=6
         ELSEIF(ITURN.EQ.2) THEN
           TURND(I)=1
           IF(TURN(CELL(I)).GT.6)TURND(I)=9
         ELSEIF(ITURN.EQ.3) THEN
           TURND(I)=2
           IF(TURN(CELL(I)).GT.6)TURND(I)=10
         ELSEIF(ITURN.EQ.4) THEN
           TURND(I)=3
           IF(TURN(CELL(I)).GT.6)TURND(I)=11
         ELSEIF(ITURN.EQ.5) THEN
           TURND(I)=4
           IF(TURN(CELL(I)).GT.6)TURND(I)=12
         ELSEIF(ITURN.EQ.6) THEN
           TURND(I)=5
           IF(TURN(CELL(I)).GT.6)TURND(I)=7
         ELSE
           CALL XABORT('DUTURN : TURN DUPLICATION ALGORITHME ERROR ')
         ENDIF
   25 CONTINUE
*
*      DUPLICATON DES AUTRES COURONES
*
      JCEL=3
      DO 30 IC=3,NCOUR
         NCS=INT(AINT((REAL(IC)+1.)/2.))
         NCEL1=IFONC(IC,0)
         KCEL=JCEL+NCS-1
         DO 32 IN=JCEL,KCEL
            NTURN(IN)=TURN(IN)
            NUM(IN)=ITAB(IN)
   32    CONTINUE
         LAUX=1
         TURND(NCEL1)=TURN(KCEL)
      DO 35 JROT=0,11
         NCEL2=NCEL1+NCS-1
         IF(MOD(IC,2).EQ.0) THEN
           IF(LAUX.EQ.0) THEN
             NCEL2=NCEL2+1
             LAUX=1
           ELSE
             LAUX=0
           ENDIF
         ENDIF
         IF(JROT.EQ.11)NCEL2=NCEL2-1
         DO 33 J=NCEL1+1,NCEL2
            ITURN=NUM(CELL(J))
            KTURN=JROT+TAB(ITURN)
            IF(KTURN.GT.12) KTURN=KTURN-12
            IF(KTURN.GT.6) KTURN=KTURN-6
            NUM(CELL(J))=KTURN
            IF(KTURN.EQ.1) THEN
              TURND(J)=6
              IF(NTURN(CELL(J)).LE.6)TURND(J)=8
            ELSEIF(KTURN.EQ.2) THEN
              TURND(J)=1
              IF(NTURN(CELL(J)).LE.6)TURND(J)=9
            ELSEIF(KTURN.EQ.3) THEN
              TURND(J)=2
              IF(NTURN(CELL(J)).LE.6)TURND(J)=10
            ELSEIF(KTURN.EQ.4) THEN
              TURND(J)=3
              IF(NTURN(CELL(J)).LE.6)TURND(J)=11
            ELSEIF(KTURN.EQ.5) THEN
              TURND(J)=4
              IF(NTURN(CELL(J)).LE.6)TURND(J)=12
            ELSEIF(KTURN.EQ.6) THEN
              TURND(J)=5
              IF(NTURN(CELL(J)).LE.6)TURND(J)=7
            ELSE
              CALL XABORT('DUTURN : INVALID ORIENTATION 2 ')
            ENDIF
            NTURN(CELL(J))=TURND(J)
   33    CONTINUE
         NCEL1=NCEL2
   35 CONTINUE
      JCEL=KCEL+1
   30 CONTINUE
      GO TO 200
*
   40 CONTINUE
*
*        DUPLICATION DE L'ORIENTATION DES GEOMETRIES SA60 ET SB60
*
      TURND(1)=TURN(1)
      JCEL=2
      LGSA6=IHEX.EQ.2
      DO 55 IC=2,NCOUR
         NCS=IC
         NCEL1=IFONC(IC,0)
         NCEL10=0
         IF(.NOT.LGSA6) THEN
             NCS=2*NINT(REAL(IC)/2.)-1
             NCEL10=NCEL1
             NCEL1=NCEL1+NINT(REAL(IC+1)/2.)-1
         ENDIF
         KCEL=JCEL+NCS-1
         DO 50 IN=JCEL,KCEL
            NTURN(IN)=TURN(IN)
            NUM(IN)=ITAB(IN)
   50    CONTINUE
         IF(LGSA6) THEN
           TURND(NCEL1)=TURN(KCEL)
         ELSE
           KKK=KCEL-NCEL1+NCEL10-1
           NCFIN=NCEL1
           IF(MOD(IC,2).EQ.0)NCFIN=NCEL1-1
           DO 555 IK=NCEL10,NCFIN
             KKK=KKK+1
             TURND(IK)=TURN(KKK)
  555      CONTINUE
         ENDIF
         DO 54 JROT=0,5
            NCEL2=NCEL1+NCS-1
            IF(JROT.EQ.5) THEN
*              NCEL2=NCEL2-1
               IF(.NOT.LGSA6) THEN
                 NCEL2=NCEL2-NINT(REAL(NCS)/2.)
               ELSE
                 NCEL2=NCEL2-1
               ENDIF
            ENDIF
            DO 52 J=NCEL1,NCEL2
               ITURN=NUM(CELL(J))
               KTURN=0
               IF(LGSA6) THEN
                 KTURN=TAB(ITURN)+4*JROT
               ELSE
                 KTURN=TAB6(ITURN)+2*JROT
               ENDIF
               IF(KTURN.GT.24)KTURN=KTURN-24
               IF(KTURN.GT.12)KTURN=KTURN-12
               IF(KTURN.GT.6) KTURN=KTURN-6
               IF(.NOT.LGSA6) NUM(CELL(J))=KTURN
               ITTD=0
               IF(KTURN.EQ.1) THEN
                  ITTD=6
                  IF(NTURN(CELL(J)).LE.6)ITTD=8
               ELSEIF(KTURN.EQ.2) THEN
                  ITTD=1
                  IF(NTURN(CELL(J)).LE.6)ITTD=9
               ELSEIF(KTURN.EQ.3) THEN
                  ITTD=2
                  IF(NTURN(CELL(J)).LE.6)ITTD=10
               ELSEIF(KTURN.EQ.4) THEN
                  ITTD=3
                  IF(NTURN(CELL(J)).LE.6)ITTD=11
               ELSEIF(KTURN.EQ.5) THEN
                  ITTD=4
                  IF(NTURN(CELL(J)).LE.6)ITTD=12
               ELSEIF(KTURN.EQ.6) THEN
                  ITTD=5
                  IF(NTURN(CELL(J)).LE.6)ITTD=7
               ELSE
                  CALL XABORT('DUTURN : INVALID ORIENTATION 3 ')
               ENDIF
               IF(J.EQ.NCEL1) THEN
                 IF(LGSA6)GOTO 51
                 IF(MOD(IC,2).NE.0) GOTO 51
               ENDIF
               TURND(J)=ITTD
   51          NTURN(CELL(J))=ITTD
   52      CONTINUE
           NCEL1=NCEL2
           IF(.NOT.LGSA6) THEN
             IF(MOD(IC,2).EQ.0)NCEL1=NCEL1+1
           ENDIF
   54    CONTINUE
         JCEL=KCEL+1
   55  CONTINUE
       GO TO 200
*
   60 CONTINUE
*
*             DUPLICATION DE L'ORIENTATION DE LA GEOMETRIE S90
*
      TURND(1)=TURN(1)
      JCEL=2
      DO 75 IC=2,NCOUR
         NCS=IC+INT(AINT(REAL((IC+1)/2)))-1
         NCEL1=IFONC(IC,1)
         KCEL=JCEL+NCS-1
         DO 70 IN=JCEL,KCEL
            NTURN(IN)=TURN(IN)
            NUM(IN)=ITAB(IN)
   70    CONTINUE
         NCEL0=IFONC(IC,0)
         KKK=KCEL-NCEL1+NCEL0
         DO 71 IK=NCEL0,NCEL1
            KKK=KKK+1
            TURND(IK)=TURN(KKK)
   71    CONTINUE
         DO 74 JROT=0,3
            NCEL2=NCEL1+NCS-1
            IF(JROT.EQ.3) NCEL2=NCEL1+INT(AINT(REAL((IC+1)/2)))-2
            DO 72 J=NCEL1,NCEL2
               ITURN=NUM(CELL(J))
               KTURN=TAB9(ITURN)+3*JROT
               IF(KTURN.GT.12)KTURN=KTURN-12
               IF(KTURN.GT.6) KTURN=KTURN-6
               NUM(CELL(J))=KTURN
               IF(KTURN.EQ.1) THEN
                 TURND(J)=6
                 IF(NTURN(CELL(J)).LE.6)TURND(J)=8
               ELSEIF(KTURN.EQ.2) THEN
                 TURND(J)=1
                 IF(NTURN(CELL(J)).LE.6)TURND(J)=9
               ELSEIF(KTURN.EQ.3) THEN
                 TURND(J)=2
                 IF(NTURN(CELL(J)).LE.6)TURND(J)=10
               ELSEIF(KTURN.EQ.4) THEN
                 TURND(J)=3
                 IF(NTURN(CELL(J)).LE.6)TURND(J)=11
               ELSEIF(KTURN.EQ.5) THEN
                 TURND(J)=4
                 IF(NTURN(CELL(J)).LE.6)TURND(J)=12
               ELSEIF(KTURN.EQ.6) THEN
                 TURND(J)=5
                 IF(NTURN(CELL(J)).LE.6)TURND(J)=7
               ELSE
                 CALL XABORT('DUTURN : INVALID ORIENTATION 4 ')
               ENDIF
               NTURN(CELL(J))=TURND(J)
   72       CONTINUE
            NCEL1=NCEL2
            IF(MOD(IC,2).EQ.0) THEN
            IF(JROT.EQ.0.OR.JROT.EQ.2) NCEL1=NCEL1+1
            ENDIF
   74    CONTINUE
         JCEL=KCEL+1
   75 CONTINUE
*
      GO TO 200
*
   80 CONTINUE
*
*        DUPLICATION DE L'ORIENTATION DE LA SYMETRIE R120
*
      TURND(1)=TURN(1)
      JCEL=2
      DO 95 IC=2,NCOUR
         NCS=2*(IC-1)
         NCEL1=IFONC(IC,1)
         NCEL0=IFONC(IC,0)
         KCEL=JCEL+NCS-1
         DO 90 IN=JCEL,KCEL
            NTURN(IN)=TURN(IN)
            NUM(IN)=ITAB(IN)
  90     CONTINUE
         KK=KCEL
         DO 91 I=NCEL1,NCEL0,-1
            TURND(I)=TURN(KK)
            KK=KK-1
  91     CONTINUE
         NCEL1=NCEL1+1
      DO 94 JROT=0,1
         NCEL2=NCEL1+NCS-1
         DO 92 J=NCEL1,NCEL2
            ITURN=NUM(CELL(J))
            KTURN=TAB12(ITURN)
            NUM(CELL(J))=KTURN
            IF(KTURN.EQ.1) THEN
              TURND(J)=6
              IF(NTURN(CELL(J)).GT.6)TURND(J)=8
            ELSEIF(KTURN.EQ.2) THEN
              TURND(J)=1
              IF(NTURN(CELL(J)).GT.6)TURND(J)=9
            ELSEIF(KTURN.EQ.3) THEN
              TURND(J)=2
              IF(NTURN(CELL(J)).GT.6)TURND(J)=10
            ELSEIF(KTURN.EQ.4) THEN
              TURND(J)=3
              IF(NTURN(CELL(J)).GT.6)TURND(J)=11
            ELSEIF(KTURN.EQ.5) THEN
              TURND(J)=4
              IF(NTURN(CELL(J)).GT.6)TURND(J)=12
            ELSEIF(KTURN.EQ.6) THEN
              TURND(J)=5
              IF(NTURN(CELL(J)).GT.6)TURND(J)=7
            ELSE
              CALL XABORT('DUTURN : INVALID ORIENTATION 5 ')
            ENDIF
            NTURN(CELL(J))=TURND(J)
   92    CONTINUE
         NCEL1=NCEL2+1
   94 CONTINUE
      NCC=NCEL2+1
      DO 93 L=KK,JCEL,-1
        TURND(NCC)=TURN(L)
        NCC=NCC+1
   93 CONTINUE
      JCEL=KCEL+1
   95 CONTINUE
*
      GO TO 200
*
  100 CONTINUE
*
*         DUPLICATION DE L'ORIENTATION DES SYMETRIES R180,SA180 ET SB180
*
      TURND(1)=TURN(1)
      LGR8=.FALSE.
      LGSA=.FALSE.
      LGSB=.FALSE.
      IF(IHEX.EQ.6) THEN
        LGR8=.TRUE.
      ELSEIF(IHEX.EQ.7) THEN
        LGSA=.TRUE.
      ELSEIF(IHEX.EQ.8) THEN
        LGSB=.TRUE.
      ENDIF
      JCEL=2
      DO 115 IC=2,NCOUR
         NCEL1=IFONC(IC,1)
         NCEL10=NCEL1
         NCEL0=IFONC(IC,0)
         NCS=0
         IF(LGR8) THEN
           NCS=3*(IC-1)
           NCEL1=NCEL1+1
         ELSEIF(LGSA) THEN
           NCS=3*IC-2
         ELSEIF(LGSB) THEN
           NCC=INT(AINT(REAL(IC+1)/2.))-1
           NCS=2*IC-1+2*NCC
           NCEL1=NCEL1+IC+NCC
           NCEL10=NCEL1-1
         ENDIF
         KCEL=JCEL+NCS-1
         DO 110 IN=JCEL,KCEL
            NTURN(IN)=TURN(IN)
            NUM(IN)=ITAB(IN)
  110    CONTINUE
         NCEL2=NCEL1+NCS-1
         IF(LGSB) THEN
           IF(MOD(IC,2).NE.0)NCEL2=NCEL2-2
         ENDIF
         KK=KCEL
         DO 111 IZ=NCEL10,NCEL0,-1
           TURND(IZ)=TURN(KK)
           KK=KK-1
  111    CONTINUE
         LL=NCEL2
         DO 112 IZ=JCEL,KK
           LL=LL+1
           TURND(LL)=TURN(IZ)
  112    CONTINUE
         DO 102 J=NCEL1,NCEL2
            ITURN=NUM(CELL(J))
            KTURN=0
            IF(LGR8) THEN
              KTURN=TABR8(ITURN)
            ELSEIF(LGSA) THEN
              KTURN=TABA8(ITURN)
            ELSEIF(LGSB) THEN
              KTURN=TABB8(ITURN)
            ENDIF
            IF(KTURN.EQ.1) THEN
              TURND(J)=6
              IF(LGR8) THEN
                IF(NTURN(CELL(J)).GT.6) TURND(J)=8
              ELSE
                IF(NTURN(CELL(J)).LE.6) TURND(J)=8
              ENDIF
            ELSEIF(KTURN.EQ.2) THEN
              TURND(J)=1
              IF(LGR8) THEN
                IF(NTURN(CELL(J)).GT.6) TURND(J)=9
              ELSE
                IF(NTURN(CELL(J)).LE.6) TURND(J)=9
              ENDIF
            ELSEIF(KTURN.EQ.3) THEN
              TURND(J)=2
              IF(LGR8) THEN
                IF(NTURN(CELL(J)).GT.6) TURND(J)=10
              ELSE
                IF(NTURN(CELL(J)).LE.6) TURND(J)=10
              ENDIF
            ELSEIF(KTURN.EQ.4) THEN
              TURND(J)=3
              IF(LGR8) THEN
                IF(NTURN(CELL(J)).GT.6) TURND(J)=11
              ELSE
                IF(NTURN(CELL(J)).LE.6) TURND(J)=11
              ENDIF
            ELSEIF(KTURN.EQ.5) THEN
              TURND(J)=4
              IF(LGR8) THEN
                IF(NTURN(CELL(J)).GT.6) TURND(J)=12
              ELSE
                IF(NTURN(CELL(J)).LE.6) TURND(J)=12
              ENDIF
            ELSEIF(KTURN.EQ.6) THEN
              TURND(J)=5
              IF(LGR8) THEN
                IF(NTURN(CELL(J)).GT.6) TURND(J)=7
              ELSE
                IF(NTURN(CELL(J)).LE.6) TURND(J)=7
              ENDIF
            ELSE
              CALL XABORT('DUTURN : INVALID ORIENTATION 6 ')
            ENDIF
  102    CONTINUE
         JCEL=KCEL+1
  115 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  200 DEALLOCATE(ITAB,NTURN,NUM)
      RETURN
      END
