*DECK ALROOT
      SUBROUTINE ALROOT(A,M,ROOTS,LFAIL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* find the roots of a polynomial.
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
* A       polynomial coefficients.
* M       polynomial order.
*
*Parameters: output
* ROOTS   complex roots.
* LFAIL   flag set to .true. in case of failure.
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER M
      DOUBLE PRECISION A(M+1)
      COMPLEX*16 ROOTS(M)
      LOGICAL LFAIL
*----
*  LOCAL VARIABLES
*----
      COMPLEX*16 AAA,BBB,SQ1,TEST,SQRTM3
      PARAMETER (SQRTM3=(0.0,1.73205080756888))
      DOUBLE PRECISION EPS
      PARAMETER (EPS=1.D-6,MAXM=101)
      COMPLEX*16 AD(MAXM),X,B,C
*
      LFAIL=.FALSE.
      IF(M+1.GT.MAXM) CALL XABORT('ALROOT: INSUFFICIENT STORAGE.')
      IF(A(M+1).EQ.0.0D0) CALL XABORT('ALROOT: INVALID COEFFICIENT.')
      IF(M.EQ.1) THEN
         ROOTS(1)=-A(1)/A(2)
      ELSE IF(M.EQ.2) THEN
         CQ=A(2)/A(3)
         DQ=A(1)/A(3)
         AAA=CQ*CQ-4.0D0*DQ
         AAA=SQRT(AAA)
         ROOTS(1)=-0.5D0*(CQ+AAA)
         ROOTS(2)=-0.5D0*(CQ-AAA)
      ELSE IF(M.EQ.3) THEN
         IF(A(1).EQ.0.0) THEN
            CQ=A(3)/A(4)
            DQ=A(2)/A(4)
            AAA=CQ*CQ-4.0D0*DQ
            AAA=SQRT(AAA)
            ROOTS(1)=0.0
            ROOTS(2)=-0.5D0*(CQ+AAA)
            ROOTS(3)=-0.5D0*(CQ-AAA)
         ELSE
            BQ=A(3)/A(4)
            CQ=A(2)/A(4)
            DQ=A(1)/A(4)
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
         ENDIF
      ELSE IF(M.EQ.4) THEN
         CALL ALQUAR(A,ROOTS)
      ELSE
         DO 10 J=1,M+1
         AD(J)=CMPLX(A(J),0.D0,KIND=KIND(AD))
   10    CONTINUE
         DO 25 J=M,1,-1
         X=CMPLX(0.D0,0.D0,KIND=KIND(X))
         CALL ALGUER(AD,J,X,ITS,LFAIL)
         IF(LFAIL) RETURN
         IF(ABS(DIMAG(X)).LE.2.D0*EPS**2*ABS(DBLE(X)))
     1   X=CMPLX(DBLE(X),0.D0,KIND=KIND(X))
         ROOTS(J)=X
         B=AD(J+1)
         DO 20 JJ=J,1,-1
         C=AD(JJ)
         AD(JJ)=B
         B=X*B+C
   20    CONTINUE
   25    CONTINUE
         DO 30 J=1,M+1
         AD(J)=CMPLX(A(J),0.D0,KIND=KIND(AD))
   30    CONTINUE
         DO 40 J=1,M
         CALL ALGUER(AD,M,ROOTS(J),ITS,LFAIL)
         IF(LFAIL) RETURN
   40    CONTINUE
      ENDIF
*
      DO 70 J=2,M
      X=ROOTS(J)
      DO 50 I=J-1,1,-1
      IF(DBLE(ROOTS(I)).LE.DBLE(X)) GOTO 60
      ROOTS(I+1)=ROOTS(I)
   50 CONTINUE
      I=0
   60 ROOTS(I+1)=X
   70 CONTINUE
      RETURN
      END
