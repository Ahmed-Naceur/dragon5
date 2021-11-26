*DECK DETFIL
      SUBROUTINE DETFIL(Y,X,XL,TC,DT,N)
*
*----------------------------------------------------------------------
*
*Purpose:
* Filters the n values of x vector for TC seconds.
* Formulations are taken from the expression of the different
* equation of a filter with a linear variation of x in DT time
* step.
*
*Author(s): 
* xxx
*
*Parameters:
* Y   real variables Y(I) at previous time
* X   real variables to filtered
* XL  real variables X(I) at previous time
* TC  filter time constant
* DT  time step between two calculations
* N   dimension of the vectors X,XL,Y
*
*--------------------------------------------------------------------
*
      IMPLICIT NONE
      INTEGER  N,I
      REAL     Y(N),X(N),XL(N),TC,DT,AA,A,B,C
*
*     COMPUTE PONDERATION FACTORS FOR Y,X ET XL
*
      AA =  - DT / TC
      A  = EXP ( AA )
      B  = 1. - A
      C  = 1. - B * TC/DT
*
*     COMPUTE NEW Y
*
      DO 10 I=1,N
*
        Y(I)  = A * Y(I) + ( B - C ) * XL(I) + C * X(I)
        XL(I) = X(I)
*
 10   CONTINUE
      RETURN
      END
