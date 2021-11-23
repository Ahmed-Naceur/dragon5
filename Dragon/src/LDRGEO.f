*DECK LDRGEO
      LOGICAL FUNCTION LDRGEO(IPGEOM,GEON1,GEON2,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compare two sub-geometries stored on LCM (do not compare mixture
* numbers).
*
*Copyright:
* Copyright (C) 1002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPGEOM  pointer to the geometry (L_GEOM signature).
* GEON1   name of the first sub-geometry.
* GEON2   name of the second sub-geometry.
* IMPX    print flag (impx=0 for no print).
*
*Parameters: output
* LDRGEO  equality flag (=.true. if the two geometries are identical).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPGEOM
      INTEGER IMPX
      CHARACTER GEON1*12,GEON2*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXLEV=50)
      TYPE(C_PTR) IPLIS1,IPLIS2,KDATA1(MAXLEV),KDATA2(MAXLEV)
      CHARACTER NAMT*12,GEON3*12,GEON4*12,HNAME*12,NAMMY1*12,NAMMY2*12,
     1 CTMP1*4,CTMP2*4,PATH(MAXLEV)*12,FIRST(MAXLEV)*12,HSMG*131
      INTEGER IGO(MAXLEV)
      LOGICAL EMPTY,LCM
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDATA1,IDATA2
      REAL, ALLOCATABLE, DIMENSION(:) :: RDATA1,RDATA2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DDATA1,DDATA2
*
      CALL LCMLEN(IPGEOM,GEON1,ILON1,ITY1)
      IF(ILON1.EQ.0) CALL XABORT('LDRGEO: UNKNOWN GEOMETRY (1).')
      CALL LCMLEN(IPGEOM,GEON1,ILON1,ITY1)
      IF(ILON1.EQ.0) CALL XABORT('LDRGEO: UNKNOWN GEOMETRY (2).')
      GEON3=GEON1
      GEON4=GEON2
      DO 120 IORDER=1,2
      IPLIS1=LCMGID(IPGEOM,GEON3)
      IPLIS2=LCMGID(IPGEOM,GEON4)
      LDRGEO=.TRUE.
      ILEV=1
      KDATA1(1)=IPLIS1
      KDATA2(1)=IPLIS2
      IGO(1)=3
*
* ASSOCIATIVE TABLE.
   10 CALL LCMINF(IPLIS1,HNAME,NAMMY1,EMPTY,ILONG,LCM)
      CALL LCMINF(IPLIS2,HNAME,NAMMY2,EMPTY,ILONG,LCM)
      IF(EMPTY) GO TO (100,100,110),IGO(ILEV)
      NAMT=' '
      CALL LCMNXT(IPLIS1,NAMT)
*
      FIRST(ILEV)=NAMT
   20 CALL LCMLEN(IPLIS1,NAMT,ILON1,ITY1)
      CALL LCMLEN(IPLIS2,NAMT,ILON2,ITY2)
      IF((ILON1.NE.ILON2).OR.(ITY1.NE.ITY2)) THEN
         LDRGEO=.FALSE.
         IF(IMPX.GT.0) WRITE (6,130) GEON3,GEON4,NAMT
         RETURN
      ENDIF
      IF(ITY1.EQ.0) THEN
*        ASSOCIATIVE TABLE DATA.
         ILEV=ILEV+1
         IF(ILEV.GT.MAXLEV) THEN
            WRITE(HSMG,'(2A,A12,A)') 'LDRGEO: TOO MANY DIRECTORY ',
     1      'LEVELS ON ''',HNAME,'''.'
            CALL XABORT(HSMG)
         ENDIF
         KDATA1(ILEV)=LCMGID(IPLIS1,NAMT)
         KDATA2(ILEV)=LCMGID(IPLIS2,NAMT)
         PATH(ILEV)=NAMT
         IPLIS1=KDATA1(ILEV)
         IPLIS2=KDATA2(ILEV)
         IGO(ILEV)=1
         GO TO 10
      ELSE IF(ITY1.LE.6) THEN
         IF(ITY1.EQ.1) THEN
*           INTEGER DATA.
            ALLOCATE(IDATA1(ILON1),IDATA2(ILON2))
            CALL LCMGET(IPLIS1,NAMT,IDATA1)
            CALL LCMGET(IPLIS2,NAMT,IDATA2)
            IF((NAMT.NE.'MIX').AND.(NAMT.NE.'STATE-VECTOR')) THEN
               DO 40 I=1,ILON1
               LDRGEO=LDRGEO.AND.(IDATA1(I).EQ.IDATA2(I))
   40          CONTINUE
            ELSE IF(NAMT.EQ.'STATE-VECTOR') THEN
               DO 50 I=1,6
               LDRGEO=LDRGEO.AND.(IDATA1(I).EQ.IDATA2(I))
   50          CONTINUE
            ENDIF
            DEALLOCATE(IDATA2,IDATA1)
         ELSE IF(ITY1.EQ.2) THEN
*           SINGLE PRECISION DATA.
            ALLOCATE(RDATA1(ILON1),RDATA2(ILON2))
            CALL LCMGET(IPLIS1,NAMT,RDATA1)
            CALL LCMGET(IPLIS2,NAMT,RDATA2)
            ZMAX=0.0
            DO 60 I=1,ILON1
            ZMAX=MAX(ZMAX,ABS(RDATA1(I)),ABS(RDATA2(I)))
   60       CONTINUE
            IF(ZMAX.EQ.0.0) ZMAX=1.0
            DO 70 I=1,ILON1
            EPS=ABS(RDATA1(I)-RDATA2(I))/ZMAX
            LDRGEO=LDRGEO.AND.(EPS.LT.1.0E-4)
   70       CONTINUE
            DEALLOCATE(RDATA2,RDATA1)
         ELSE IF(ITY1.EQ.3) THEN
*           CHARACTER*4 DATA.
            ALLOCATE(IDATA1(ILON1),IDATA2(ILON2))
            CALL LCMGET(IPLIS1,NAMT,IDATA1)
            CALL LCMGET(IPLIS2,NAMT,IDATA2)
            DO 80 I=1,ILON1
            WRITE(CTMP1,'(A4)') IDATA1(I)
            WRITE(CTMP2,'(A4)') IDATA2(I)
            LDRGEO=LDRGEO.AND.(CTMP1.EQ.CTMP2)
   80       CONTINUE
            DEALLOCATE(IDATA2,IDATA1)
         ELSE IF(ITY1.EQ.4) THEN
*           DOUBLE PRECISION DATA.
            ALLOCATE(DDATA1(ILON1),DDATA2(ILON2))
            CALL LCMGET(IPLIS1,NAMT,DDATA1)
            CALL LCMGET(IPLIS2,NAMT,DDATA2)
            ZMAX=0.0
            DO 85 I=1,ILON1
            ZMAX=MAX(ZMAX,REAL(ABS(DDATA1(I))),REAL(ABS(DDATA2(I))))
   85       CONTINUE
            IF(ZMAX.EQ.0.0) ZMAX=1.0
            DO 90 I=1,ILON1
            EPS=ABS(REAL(DDATA1(I)-DDATA2(I)))/ZMAX
            LDRGEO=LDRGEO.AND.(EPS.LT.1.0E-4)
   90       CONTINUE
            DEALLOCATE(DDATA2,DDATA1)
         ELSE
            CALL XABORT('LDRGEO: INVALID DATA TYPE.')
         ENDIF
         IF(.NOT.LDRGEO) THEN
            LDRGEO=.FALSE.
            IF(IMPX.GT.0) WRITE (6,130) GEON3,GEON4,NAMT
            RETURN
         ENDIF
      ENDIF
      CALL LCMNXT(IPLIS1,NAMT)
      IF(NAMT.NE.FIRST(ILEV)) GO TO 20
      GO TO (100,100,110),IGO(ILEV)
*
  100 NAMT=PATH(ILEV)
      ILEV=ILEV-1
      IPLIS1=KDATA1(ILEV)
      IPLIS2=KDATA2(ILEV)
      CALL LCMNXT(IPLIS1,NAMT)
      IF(NAMT.NE.FIRST(ILEV)) GO TO 20
      GO TO (100,100,110),IGO(ILEV)
  110 GEON3=GEON2
      GEON4=GEON1
  120 CONTINUE
      RETURN
*
  130 FORMAT (/34H LDRGEO: COMPARISON OF GEOMETRIES ,A12,5H AND ,A12,
     1 16H --- LCM BLOCK ',A12,15H' IS DIFFERENT.)
      END
