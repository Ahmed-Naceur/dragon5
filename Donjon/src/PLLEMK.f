*DECK PLLEMK
      SUBROUTINE PLLEMK(N,M,EPS,IMPR,P,IROW,ICOL,IERR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* This subroutine solve the parametric linear complementary problem.
* PLLEMK = Linear Programmation LEMKe
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
* N       number of control variables.
* M       number of constraints.
* EPS     tolerence used for pivoting.
* IMPR    print flag.
* P       coefficient matrix.
* IROW    permutation vector for row elements.
* ICOL    permutation vector for column elements.
*
*Parameters: ouput
* P       coefficient matrix.
* IROW    permutation vector for row elements.
* ICOL    permutation vector for column elements.
* IERR    return code (=0: normal completion).
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER     M,N,IERR,IMPR,IROW(N),ICOL(N+1)
      DOUBLE PRECISION  EPS,P(N,M)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION  WWW,S
      INTEGER     MNOP,NP1,IC,L,I,J,K,JJ,J1,LGAR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PLJ
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PLJ(M))
*
      MNOP = N + M
      NP1  = N + 1
      IC = 0
*
      S = 0.0D0
      L = 0
*
      DO 50 I=1,N
         P(I,NP1) = 1.0D0
         IF(P(I,M).GE.0.0D0) GO TO 50
         WWW = P(I,M)
         IF(WWW.LE.S) THEN
            L = I
            S = WWW
         ENDIF
   50 CONTINUE
*
      IF(L.EQ.0) GO TO 290
      LGAR = L
*
      K = NP1
*
   90 IC = IC + 1
      IF(IC.GT.MNOP) THEN
         IERR = 1
         GOTO 300
      ENDIF
*
      IF(IMPR.GE.5) WRITE (6,1000) L,K
      IF(IMPR.GE.6) THEN
         DO 110 J=1,M,12
            JJ = MIN0(J+11,M)
            WRITE (6,3000) (J1,J1=J,JJ)
            DO 100 I=1,N
               WRITE (6,4000) I,(P(I,J1),J1=J,JJ)
  100       CONTINUE
            IF(JJ.NE.M) WRITE (6,5000)
  110    CONTINUE
      ENDIF
*
      CALL PLPIVT(N,M,L,K,P,IROW,ICOL)
*
      IF(ICOL(K).EQ.-NP1) GOTO 290
      DO 120 J=1,NP1
         IF(J.EQ.K) GOTO 120
         IF(IABS(ICOL(K)).EQ.IABS(ICOL(J))) THEN
            K = J
            GOTO 130
         ENDIF
  120 CONTINUE
*
      IERR = 2
      GOTO 300
*
  130 DO 140 I=1,N
         IF(P(I,K)/P(I,M).LT.-EPS) THEN
            PLJ(I) = -P(I,M)/P(I,K)
         ELSE
            PLJ(I) = 1.0D50
         ENDIF
  140 CONTINUE
*
      S = PLJ(LGAR)
      L = LGAR
*
      IF(IMPR.GE.7) THEN
         WRITE (6,*) 'K=',K
         WRITE (6,6000) L,(PLJ(I),I=1,N)
      ENDIF
*
      DO 280 I=1,N
      WWW = PLJ(I)
      IF((ABS(WWW-S).GT.ABS(S)*EPS).AND.(WWW.LT.S)) THEN
         S = WWW
         L = I
      ENDIF
  280 CONTINUE
*
      IF(S.EQ.1.0D50) THEN
         IERR = 3
         GOTO 300
      ENDIF
*
      GOTO 90
  290 IERR = 0
*
  300 IF(IMPR.GE.6) THEN
         WRITE (6,2000) IC
         DO 315 J=1,M,12
            JJ = MIN0(J+11,M)
            WRITE (6,3000) (J1,J1=J,JJ)
            DO 310 I=1,N
               WRITE (6,4000) I,(P(I,J1),J1=J,JJ)
  310       CONTINUE
            IF(JJ.NE.M) WRITE (6,5000)
  315    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(PLJ)
      RETURN
*
1000  FORMAT (//,10X,'+ + + + +   MATRIX P(N,N+3) :',5X,'PIVOT = (',
     >           I4,' ,',I4,' )'/)
2000  FORMAT (//,10X,'NUMBER OF PIVOTS =',I5,
     >        //,10X,'+ + + + +   FINAL MATRIX P(N,N+3) :'/)
3000  FORMAT (5X,12I12)
4000  FORMAT (1X,I4,1P,10E12.4/(5X,1P,10E12.4))
5000  FORMAT (//)
6000  FORMAT (//,10X,'LGAR =',I4/(1X,1P,8E12.4))
      END
