*DECK AFERF
      FUNCTION AFERF (XX)
C
C-----------------------------------------------------------------------
C
C ERROR FUNCTION.
C
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C----
C FUNCTION
C----
      REAL  AFERF
C----
C VARIABLES
C----
      REAL     A0,A1,A2,A3,A4,A5,A6,B0,B1,B2,B3,
     &         B4,B5,B6,B7,D1,D2,X,XX,Z,Y
      INTEGER  IS
C----
      DATA A0 / 0.12837912 /, A1 / -0.37612384 /, A2 / 0.11281417 /
     1 ,   A3 / -0.26782287E-01 /, A4 / 0.50805667E-02 /
     2 ,   A5 / -0.72514300E-03 /, A6 / 0.58167427E-04 /
     3 ,   B0 /  0.39141751E-02 /, B1 / -0.17582889E-01 /
     4 ,   B2 /  0.35873539E-01 /, B3 / -0.42869095E-01 /
     5 ,   B4 /  0.32161925E-01 /, B5 / -0.11846550E-01 /
     6 ,   B6 /  0.30705572E-02 /, B7 /  0.59813837E-02 /
     7 ,   D1,D2 / 1.317, 2.040001 /
C
      X=XX
      IF (X.GE.0.) THEN
         IS=2
      ELSE
         IS=1
         X=-X
      ENDIF
      IF (X.GT.D1) GO TO 8
      Z=X*X
      Y=(((((A6*Z+A5)*Z+A4)*Z+A3)*Z+A2)*Z+A1)*Z+A0
      Y=X*Y+X
      GO TO 6
   8  IF (X.GT.D2) GO TO 10
      X=X-D2
      Y=((((((B7*X+B6)*X+B5)*X+B4)*X+B3)*X+B2)*X+B1)*X+B0
      GO TO 7
  10  Y=1.
      GO TO 6
   7  Y=1.-Y
   6  GO TO (17,18),IS
  17  AFERF=-Y
      GO TO 30
  18  AFERF=Y
  30  RETURN
      END
