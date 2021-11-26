*DECK ALQUAR
      SUBROUTINE ALQUAR(A,ROOTS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute the roots of the real quartic polynomial defined as
* A(1)+A(2)*Z + ... + A(5)*Z**4.
* NOTE: It is assumed that A(5) is non-zero. No test is made here.
*
*Author(s): A. H. Morris, W. L. Davis, A. Miller, and R. L. Carmichael
*
*Parameters: input
* A      polynomial coefficients
*
*Parameters: output
* ROOTS  complex roots
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      DOUBLE PRECISION, INTENT(IN) :: A(5)
      COMPLEX*16, INTENT(OUT) :: ROOTS(4)
*----
*  LOCAL VARIABLES
*----
      INTEGER :: I, J
      INTEGER,DIMENSION(1) :: K
      DOUBLE PRECISION :: B, B2, C, D, E, H, P, Q, R, T, WORK, U, V,
     > V1, V2, X, XX(3), Y, BQ, CQ, DQ, AA, BB
      COMPLEX*16 :: W, AAA, BBB, SQ1, TEST, SQRTM3
      DOUBLE PRECISION,DIMENSION(4) :: TEMP
      PARAMETER (SQRTM3=(0.0,1.73205080756888))
*
      IF(A(1)==0.0) THEN
        IF(A(2).EQ.0.0) THEN
           CQ=A(4)/A(5)
           DQ=A(3)/A(5)
           AAA=CQ*CQ-4.0D0*DQ
           AAA=SQRT(AAA)
           ROOTS(1)=0.0
           ROOTS(2)=0.0
           ROOTS(3)=-0.5D0*(CQ+AAA)
           ROOTS(4)=-0.5D0*(CQ-AAA)
        ELSE
           BQ=A(4)/A(5)
           CQ=A(3)/A(5)
           DQ=A(2)/A(5)
           AA=(3.0D0*CQ-BQ**2)/3.0D0
           BB=(2.0D0*BQ**3-9.0D0*BQ*CQ+27.0D0*DQ)/27.0D0
           SQ1=BB**2/4.0D0+AA**3/27.0D0
           TEST=BB/2.0D0-SQRT(SQ1)
           IF(DBLE(TEST).EQ.0.0) THEN
              AAA=0.0D0
           ELSE IF(DBLE(TEST).GT.0.0) THEN
             AAA=-(TEST)**(1.0D0/3.0D0)
           ELSE
             AAA=(-TEST)**(1.0D0/3.0D0)
           ENDIF
           TEST=BB/2.0D0+SQRT(SQ1)
           IF(DBLE(TEST).EQ.0.0) THEN
              BBB=0.0D0
           ELSE IF(DBLE(TEST).GT.0.0) THEN
              BBB=-(TEST)**(1.0D0/3.0D0)
           ELSE
              BBB=(-TEST)**(1.0D0/3.0D0)
           ENDIF
           ROOTS(1)=0.0
           ROOTS(2)=AAA+BBB-BQ/3.0D0
           ROOTS(3)=-(AAA+BBB)/2.0D0+(AAA-BBB)*SQRTM3/2.0D0-BQ/3.0D0
           ROOTS(4)=-(AAA+BBB)/2.0D0-(AAA-BBB)*SQRTM3/2.0D0-BQ/3.0D0
        ENDIF
        RETURN
      ENDIF
*----
*  Solve a quartic equation
*----
      B = A(4)/(4.0D0*A(5))
      C = A(3)/A(5)
      D = A(2)/A(5)
      E = A(1)/A(5)
      B2 = B*B

      P = 0.5D0*(C - 6.0D0*B2)
      Q = D - 2.0D0*B*(C - 4.0D0*B2)
      R = B2*(C - 3.0D0*B2) - B*D + E
*----
*  Solve the resolvent cubic equation. the cubic has at least one
*  nonnegative real root.  if W1, W2, W3 are the roots of the cubic
*  then the roots of the original equation are
*     ROOTS = -B + CSQRT(W1) + CSQRT(W2) + CSQRT(W3)
*  where the signs of the square roots are chosen so
*  that CSQRT(W1) * CSQRT(W2) * CSQRT(W3) = -Q/8.
*----
      TEMP(1) = -Q*Q/64.0D0
      TEMP(2) = 0.25D0*(P*P - R)
      TEMP(3) =  P
      TEMP(4) = 1.0D0
      BQ=TEMP(3)
      CQ=TEMP(2)
      DQ=TEMP(1)
      AA=(3.0D0*CQ-BQ**2)/3.0D0
      BB=(2.0D0*BQ**3-9.0D0*BQ*CQ+27.0D0*DQ)/27.0D0
      SQ1=BB**2/4.0D0+AA**3/27.0D0
      TEST=BB/2.0D0-SQRT(SQ1)
      IF(DBLE(TEST).EQ.0.0) THEN
         AAA=0.0D0
      ELSE IF(DBLE(TEST).GT.0.0) THEN
         AAA=-(TEST)**(1.0D0/3.0D0)
      ELSE
         AAA=(-TEST)**(1.0D0/3.0D0)
      ENDIF
      TEST=BB/2.0D0+SQRT(SQ1)
      IF(DBLE(TEST).EQ.0.0) THEN
         BBB=0.0D0
      ELSE IF(DBLE(TEST).GT.0.0) THEN
         BBB=-(TEST)**(1.0D0/3.0D0)
      ELSE
         BBB=(-TEST)**(1.0D0/3.0D0)
      ENDIF
      ROOTS(1)=AAA+BBB-BQ/3.0D0
      ROOTS(2)=-(AAA+BBB)/2.0D0+(AAA-BBB)*SQRTM3/2.0D0-BQ/3.0D0
      ROOTS(3)=-(AAA+BBB)/2.0D0-(AAA-BBB)*SQRTM3/2.0D0-BQ/3.0D0
      IF(AIMAG(ROOTS(2)).NE.0.0D0) GO TO 60
*----
*  The resolvent cubic has only real roots.
*  Reorder the roots in increasing order.
*----
      XX(1) = DBLE(ROOTS(1))
      XX(2) = DBLE(ROOTS(2))
      XX(3) = DBLE(ROOTS(3))
      DO 25 J=2,3
      X=XX(J)
      DO 10 I=J-1,1,-1
      IF(XX(I).LE.X) GOTO 20
      XX(I+1)=XX(I)
   10 CONTINUE
      I=0
   20 XX(I+1)=X
   25 CONTINUE

      U = 0.0D0
      IF(XX(3).GT.0.0D0) U = SQRT(XX(3))
      IF(XX(2).LE.0.0D0) GO TO 41
      IF(XX(1).GE.0.0D0) GO TO 30
      IF(ABS(XX(1)).GT.XX(2)) GO TO 40
      XX(1) = 0.0D0

   30 XX(1) = SQRT(XX(1))
      XX(2) = SQRT(XX(2))
      IF(Q.GT.0.0D0) XX(1) = -XX(1)
      TEMP(1) = (( XX(1) + XX(2)) + U) - B
      TEMP(2) = ((-XX(1) - XX(2)) + U) - B
      TEMP(3) = (( XX(1) - XX(2)) - U) - B
      TEMP(4) = ((-XX(1) + XX(2)) - U) - B

      DO J=1,3
        K=MINLOC(TEMP(J:))
        IF(J.NE.K(1)) THEN
          WORK = TEMP(J)
          TEMP(J) = TEMP(K(1))
          TEMP(K(1)) = WORK
        ENDIF
      ENDDO

      IF(ABS(TEMP(1)).GE.0.1D0*ABS(TEMP(4))) GO TO 31
      T = TEMP(2)*TEMP(3)*TEMP(4)
      IF(T.NE.0.0D0) TEMP(1) = E/T
   31 ROOTS(1) = CMPLX(TEMP(1), 0.0D0, KIND=KIND(ROOTS))
      ROOTS(2) = CMPLX(TEMP(2), 0.0D0, KIND=KIND(ROOTS))
      ROOTS(3) = CMPLX(TEMP(3), 0.0D0, KIND=KIND(ROOTS))
      ROOTS(4) = CMPLX(TEMP(4), 0.0D0, KIND=KIND(ROOTS))
      RETURN

   40 V1 = SQRT(ABS(XX(1)))
      V2 = 0.0D0
      GO TO 50
   41 V1 = SQRT(ABS(XX(1)))
      V2 = SQRT(ABS(XX(2)))
      IF(Q < 0.0D0) U = -U

   50 X = -U - B
      Y = V1 - V2
      ROOTS(1) = CMPLX(X, Y, KIND=KIND(ROOTS))
      ROOTS(2) = CMPLX(X,-Y, KIND=KIND(ROOTS))
      X =  U - B
      Y = V1 + V2
      ROOTS(3) = CMPLX(X, Y, KIND=KIND(ROOTS))
      ROOTS(4) = CMPLX(X,-Y, KIND=KIND(ROOTS))
      RETURN
*----
*  The resolvent cubic has complex roots.
*----
   60 T = DBLE(ROOTS(1))
      X = 0.0D0
      IF(T < 0.0D0) THEN
        GO TO 61
      ELSE IF(T.EQ.0.0D0) THEN
        GO TO 70
      ELSE
        GO TO 62
      ENDIF
   61 H = ABS(DBLE(ROOTS(2))) + ABS(AIMAG(ROOTS(2)))
      IF(ABS(T).LE.H) GO TO 70
      GO TO 80
   62 X = SQRT(T)
      IF(Q.GT.0.0D0) X = -X

   70 W = SQRT(ROOTS(2))
      U = 2.0D0*DBLE(W)
      V = 2.0D0*ABS(AIMAG(W))
      T =  X - B
      XX(1) = T + U
      XX(2) = T - U
      IF(ABS(XX(1)).LE.ABS(XX(2))) GO TO 71
      T = XX(1)
      XX(1) = XX(2)
      XX(2) = T
   71 U = -X - B
      H = U*U + V*V
      IF(XX(1)*XX(1) < 0.01D0*MIN(XX(2)*XX(2),H)) XX(1) = E/(XX(2)*H)
      ROOTS(1) = CMPLX(XX(1), 0.0D0, KIND=KIND(ROOTS))
      ROOTS(2) = CMPLX(XX(2), 0.0D0, KIND=KIND(ROOTS))
      ROOTS(3) = CMPLX(U, V, KIND=KIND(ROOTS))
      ROOTS(4) = CMPLX(U,-V, KIND=KIND(ROOTS))
      RETURN

   80 V = SQRT(ABS(T))
      ROOTS(1) = CMPLX(-B, V, KIND=KIND(ROOTS))
      ROOTS(2) = CMPLX(-B,-V, KIND=KIND(ROOTS))
      ROOTS(3) = ROOTS(1)
      ROOTS(4) = ROOTS(2)
      RETURN
      END
