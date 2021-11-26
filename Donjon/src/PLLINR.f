*DECK PLLINR
      SUBROUTINE PLLINR(N0,M1,COUT,APLUS,BPLUS,BINF,BSUP,XOBJ,EPS,IMPR,
     > IERR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Minimizes a linear problem using LEMKE pivot.
* PLLINR = Linear Programmation LINeaR problem resolution
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert and R. Chambon
*
*Parameters: input
* N0      number of control variables.
* M1      number of constraints.
* COUT    costs of control variables.
* APLUS   coefficient matrix for the linear constraints.
* BPLUS   right hand sides corresponding to the coefficient matrix.
* BINF    lower bounds of control variables.
* BSUP    upper bounds of control variables.
* EPS     tolerence used for pivoting.
* IMPR    print flag.
*
*Parameters: ouput
* XOBJ    control variables.
* IERR    return code (=0: normal completion).
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER       N0,M1,IMPR,IERR
      DOUBLE PRECISION  BPLUS(M1+1),BINF(N0),BSUP(N0),XOBJ(N0),EPS,
     >                  APLUS(M1+1,N0),COUT(N0)
*----
*  LOCAL VARIABLES
*----
      INTEGER       N,NP1,NP2,I,J,II,IR
      DOUBLE PRECISION  XVAL
      CHARACTER*4   ROW(7)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IROW,ICOL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: P
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IROW(2*N0+M1),ICOL(2*N0+M1+1))
      ALLOCATE(P(2*N0+M1,2*N0+M1+2))
*----
*  SET-UP AND SOLVE THE PARAMETRIC COMPLEMENTARITY PROBLEM.
*----
      N = 2*N0 + M1
      NP1 = N + 1
      NP2 = N + 2
*
      DO 10 I=1,N
         DO 15 J=1,NP2
            P(I,J)=0.0D0
  15     CONTINUE
  10  CONTINUE
*
      DO 20 I=1,N0
         P(I,NP2) = COUT(I)
         P(N0+M1+I,NP2) = BSUP(I) - BINF(I)
         P(I,N0+M1+I) =  1.0D0
         P(N0+M1+I,I) = -1.0D0
*
         DO 25 J=1,M1
            P(I,N0+J)=APLUS(J,I)
  25     CONTINUE
  20  CONTINUE
*
      DO 40 I=1,M1
         XVAL = 0.0D0
         DO 30 J=1,N0
            XVAL = XVAL + APLUS(I,J)*BINF(J)
  30     CONTINUE
         P(N0+I,NP2) = BPLUS(I) - XVAL
*
         DO 45 J=1,N0
            P(N0+I,J)=-APLUS(I,J)
  45     CONTINUE
  40  CONTINUE
*
      DO 50 I=1,N
         IROW(I)  =  I
         ICOL(I)  = -I
         P(I,NP1) = 1.0D0
  50  CONTINUE
      ICOL(NP1) = -NP1
*
      CALL PLLEMK(N,NP2,EPS,IMPR,P,IROW,ICOL,IERR)
*
      IF (IERR.GE.1) THEN
         WRITE(6,1000) IERR
         GO TO 500
      ENDIF
*
      IF (IMPR.GE.3) THEN
         WRITE(6,2000)
*
         DO 70 I=1,N,7
            II = MIN0(I+6,N)
            DO 60 J=I,II
               IF (IROW(J).LT.0) THEN
                  WRITE (ROW(J-I+1),'(1HX,I3.3)') (-IROW(J))
               ELSE
                  WRITE (ROW(J-I+1),'(1HY,I3.3)') IROW(J)
               ENDIF
   60       CONTINUE
            WRITE(6,3000) (ROW(J-I+1),P(J,NP2),J=I,II)
   70    CONTINUE
      ENDIF
*
      DO 80 I=1,N0
         XOBJ(I) = BINF(I)
   80 CONTINUE
*
      DO 90 I=1,N
         IR = -IROW(I)
         IF ((IR.GT.0).AND.(IR.LE.N0)) THEN
            XOBJ(IR)=BINF(IR)+P(I,NP2)
         ENDIF
   90 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  500 DEALLOCATE(P)
      DEALLOCATE(ICOL,IROW)
      RETURN
*
 1000 FORMAT(//,5X,'PLLINR: FAILURE OF THE LINEAR COMPLEMENTARITY SO',
     > 'LUTION (IERR=',I5,').')
 2000 FORMAT(//,5X,'SOLUTION OF THE LINEAR COMPLEMENTARITY PROBLEM :',
     > '*** X: KUHN-TUCKER MULTIPLIERS ;',5X,
     > '*** Y: SLACK VARIABLES ',/)
 3000 FORMAT(7(1X,A4,'=',E12.5),/)
      END
