*DECK LIBLEX
      SUBROUTINE LIBLEX(NELE,PNTE,ELMT,NOTX,TERP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute factors for Lagrangian interpolation.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and G. Marleau
*
*Parameters: input
* NELE    number of elements in table.
* PNTE    extrapolation point.
* ELMT    values of elements in table.
* NOTX    order of interpolation:
*         >0 order NOTX in square root of PNTE;
*         <0 order -NOTX in PNTE.
*
*Parameters: output
* TERP    extrapolation factor.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NELE,NOTX
      REAL PNTE,ELMT(NELE)
      DOUBLE PRECISION TERP(NELE)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION FACTOR,TERM,SQPNTE
*
      I1=1
      I2=NELE
      NORD=ABS(NOTX)
      CALL XDDSET(TERP,NELE,0.0D0)
      DO 100 I=1,NELE
        IF(PNTE.EQ.ELMT(I)) THEN
           TERP(I)=1.0D0
           RETURN
        ELSE IF(PNTE.GT.ELMT(I)) THEN
           I1=MAX(I-(NORD-1)/2,I1)
        ELSE IF(PNTE.LT.ELMT(I)) THEN
           I2=MIN(I+(NORD-1)/2,I2)
        ENDIF
 100  CONTINUE
*
      IF(I1.EQ.1) THEN
        I2=MIN(NELE,NORD+1)
      ELSE IF(I2.EQ.NELE) THEN
        I1=MAX(1,NELE-NORD)
      ENDIF
      FACTOR=1.0D0
      IF(NOTX.LT.0) THEN
        SQPNTE=DBLE(PNTE)
        DO 110 I=I1,I2
          IF(ABS(PNTE-ELMT(I)).LE.1.0E-5*ABS(PNTE)) THEN
            TERP(I)=1.0D0
            RETURN
          ENDIF
          FACTOR=FACTOR*( SQPNTE-DBLE(ELMT(I)) )
 110    CONTINUE
        DO 120 I=I1,I2
          TERM=FACTOR/( SQPNTE-DBLE(ELMT(I)) )
          DO 130 J=I1,I2
            IF(I.NE.J) TERM=TERM/
     >               ( DBLE(ELMT(I))-DBLE(ELMT(J)) )
 130      CONTINUE
          TERP(I)=TERM
 120    CONTINUE
      ELSE
        SQPNTE=DBLE(SQRT(PNTE))
        DO 160 I=I1,I2
          IF(ABS(PNTE-ELMT(I)).LE.1.0E-5*ABS(PNTE)) THEN
            TERP(I)=1.0D0
            RETURN
          ENDIF
          FACTOR=FACTOR*( SQPNTE-DBLE(SQRT(ELMT(I))) )
 160    CONTINUE
        DO 170 I=I1,I2
          TERM=FACTOR/( SQPNTE-DBLE(SQRT(ELMT(I))) )
          DO 180 J=I1,I2
            IF(I.NE.J) TERM=TERM/
     >               ( DBLE(SQRT(ELMT(I)))-DBLE(SQRT(ELMT(J))) )
 180      CONTINUE
          TERP(I)=TERM
 170    CONTINUE
      ENDIF
      RETURN
      END
