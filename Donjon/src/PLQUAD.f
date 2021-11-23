*DECK PLQUAD
      SUBROUTINE PLQUAD(N0,M1,APLUS,BPLUS,PDG,XDROIT,COUT,XOBJ,EPS,
     >                IMPR,IERR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Minimize a linear problem with a quadratic constraint using a
* parametric complementarity principle.
* PLQUAD = Linear Programmation with QUADratic constraint
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
* APLUS   coefficient matrix for the linear constraints.
* BPLUS   right hand sides corresponding to the coefficient matrix.
* PDG     weights assigned to control variables in the quadratic
*         constraint.
* XDROIT  quadratic constraint radius squared.
* COUT    costs of control variables.
* XOBJ    control variables.
* EPS     tolerence used for pivoting.
* IMPR    print flag.
*
*Parameters: ouput
* IERR    return code (=0: normal completion).
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER     N0,M1,IERR,IMPR
      DOUBLE PRECISION BPLUS(M1+1),PDG(N0),XOBJ(N0),EPS,XDROIT,
     > APLUS(M1+1,N0),COUT(N0)
*----
*  LOCAL VARIABLES
*----
      CHARACTER*4 ROW(7)
      DOUBLE PRECISION  PVAL,POLY0,POLY1,POLY2,XVALIR,X,OBJ,DISCRI,
     >            XROOT1,XROOT2,XVAL,XTAUU,XVALL,XVALC,OBJLIN
      INTEGER     N,NP1,NP2,NP3,I,J,K,IS,JS,IROWIS,IR,IROWR,JR,IKIT,II
      DOUBLE PRECISION  XTAU,XTAUL,UI,XMIN,XVALU
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IROW,ICOL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: U,V,WRK
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: P
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IROW(M1+1),ICOL(M1+2))
      ALLOCATE(U(M1+1),V(M1+1))
      ALLOCATE(P(M1+1,M1+4),WRK(N0))
*
      N     = M1 + 1
      NP1   = N + 1
      NP2   = N + 2
      NP3   = N + 3
*----
*  STEP 2: SET-UP AND SOLVE THE PARAMETRIC COMPLEMENTARITY PROBLEM.
*----
      DO I=1,N
         DO J=1,N0
           WRK(J) = APLUS(I,J)/PDG(J)
         ENDDO
         DO K=1,N
            PVAL = 0.0D0
            DO J=1,N0
               PVAL = PVAL + WRK(J)*APLUS(K,J)
            ENDDO
            P(I,K) = PVAL
        ENDDO
      ENDDO
*
      DO I=1,N
         IROW(I) =  I
         ICOL(I) = -I
         P(I,NP1) = 1.0D0
         P(I,NP2) = 0.0D0
         P(I,NP3) = BPLUS(I)
      ENDDO
*
      ICOL(NP1) = -NP1
      P(N,NP2)  = 1.0D0
*
      CALL PLLEMK(N,NP3,EPS,IMPR,P,IROW,ICOL,IERR)
*
      IF (IERR.GE.1) THEN
         WRITE(6,1000) IERR
         GO TO 500
      ENDIF
*
      XTAU  = 0.0
      XTAUL = 0.0
      OBJLIN = BPLUS(N)
*----
*  COMPUTE VECTOR T=(NU,PI)=U+XTAU*V
*----
  110 POLY0 = 0.0D0
      POLY1 = 0.0D0
      POLY2 = 0.0D0
*
      DO 120 I=1,N
         IR = -IROW(I)
         IF (IR.GT.0) THEN
            U(IR) = P(I,NP3)
            V(IR) = P(I,NP2)
            POLY0 = POLY0 - P(I,NP3)*BPLUS(IR)
            POLY1 = POLY1 - P(I,NP2)*BPLUS(IR)
         ELSE
            U(-IR) = 0.0
            V(-IR) = 0.0
         ENDIF
         IF (IR.EQ.N) THEN
            POLY1 = POLY1 - P(I,NP3)
            POLY2 = (-P(I,NP2))
         ENDIF
  120 CONTINUE
*
      IF (IMPR.GE.3) THEN
         DO 121 I=1,N0
            XOBJ(I) = 0.0
  121    CONTINUE
*
         DO 123 I=1,N
            UI = U(I) + XTAUL*V(I)
            IF (UI.EQ.0.0) GO TO 123
            DO 122 J=1,N0
               XOBJ(J) = XOBJ(J) - UI*APLUS(I,J)/PDG(J)
  122       CONTINUE
  123    CONTINUE
*
         X = 0.0D0
         OBJ = 0.0D0
         DO 126 J=1,N0
            X = X + PDG(J)*XOBJ(J)*XOBJ(J)
            OBJ = OBJ + XOBJ(J)*COUT(J)
  126    CONTINUE
         WRITE(6,2000) OBJ,POLY0,X,POLY1,XTAUL,POLY2,(XOBJ(J),J=1,N0)
      ENDIF
      IF ((XTAU.EQ.0.0).AND.(POLY0.LE.XDROIT)) GO TO 230
*----
*  STEP 3
*----
      DO 130 I=1,N
         IF(P(I,NP2).LT.-EPS) GO TO 140
  130 CONTINUE
      GO TO 215
*----
*  STEP 4
*----
  140 XTAUU = 1.0E+25
*
      IR = 0
      DO 150 K=I,N
         IF(P(K,NP2).GE.-EPS) GO TO 150
         XVAL = -P(K,NP3)/P(K,NP2)
         IF(XVAL.GT.XTAUU) GO TO 150
         XTAUU = XVAL
         IR = K
  150 CONTINUE
*
      XVALU = (POLY2*XTAUU + POLY1)*XTAUU + POLY0
*----
*  STEP 5
*----
      IF(XVALU.LE.XDROIT) GO TO 215
      IROWR = IABS(IROW(IR))
      JR=0
      DO 160 K=1,NP1
         IF(IABS(ICOL(K)).EQ.IROWR) THEN
            JR=K
            GO TO 170
         ENDIF
  160 CONTINUE
      IERR = 5
      GO TO 500
*
  170 XTAUL = XTAUU
      XVALL = XVALU
      IF(P(IR,JR).LE.EPS) GO TO 180
      CALL PLPIVT(N,NP3,IR,JR,P,IROW,ICOL)
      GO TO 110
*
  180 XMIN=1.0E+25
*
      XVALIR = P(IR,NP3)/P(IR,NP2)
*
      DO 190 I=1,N
         IF(P(I,JR).GE.-EPS) GO TO 190
         XVAL = -1.0D0/P(I,JR)*(P(I,NP3) - P(I,NP2)*XVALIR)
         IF(XVAL.GE.XMIN) GO TO 190
         XMIN = XVAL
         IS = I
  190 CONTINUE
*
      IF (XMIN.EQ.1.0E+25) THEN
         IERR = 6
         GO TO 500
      ENDIF
*
      IROWIS=IABS(IROW(IS))
      DO 200 JS=1,N
         IF(IABS(ICOL(JS)).EQ.IROWIS) GO TO 210
  200 CONTINUE
*
  210 CALL PLPIVT(N,NP3,IR,JS,P,IROW,ICOL)
      CALL PLPIVT(N,NP3,IS,JR,P,IROW,ICOL)
      GO TO 110
*----
*  STEP 6
*----
  215 IKIT = 0
*
  216 XTAU  = (XTAUL + XTAUU)/2.0
      IKIT = IKIT + 1
      IF (IKIT.GT.50) GOTO 217
      XVALC = ((POLY2*XTAU + POLY1)*XTAU + POLY0)/XDROIT
      IF (IMPR.GE.3) THEN
         WRITE(6,5000) XTAUL,XTAUU,XTAU,POLY0,POLY1,POLY2,XDROIT,XVALC
      ENDIF
      IF (XVALC.GT.1.0) GO TO 220
      IF (XVALC.GE.0.99999) GO TO 230
      XTAUU = XTAU
      GO TO 216
  220 XTAUL = XTAU
      GO TO 216
*----
*  STEP 6
*----
  217 XTAU  = (XTAUL + XTAUU)/2.0
      XVALC = ((POLY2*XTAU + POLY1)*XTAU + POLY0)/XDROIT
      IF (IMPR.GE.3) THEN
         WRITE(6,5000) XTAUL,XTAUU,XTAU,POLY0,POLY1,POLY2,XDROIT,XVALC
      ENDIF
*
      IF (POLY0.EQ.0.0) THEN
         IF (POLY1.EQ.0.0) THEN
            IF (POLY2.EQ.0.0) THEN
               WRITE(6,6000) POLY0,POLY1,POLY2,XDROIT
               IERR = 7
               GO TO 500
            ELSE
               IF (POLY2.LT.0.0) THEN
                  WRITE(6,6000) POLY0,POLY1,POLY2,XDROIT
                  IERR = 7
                  GO TO 500
               ENDIF
               XTAU = SQRT(XDROIT/POLY2)
            ENDIF
         ELSE IF (POLY2.EQ.0.0) THEN
            XTAU = XDROIT/POLY1
         ELSE
            DISCRI = POLY1*POLY1 + 4.*POLY2*XDROIT
            IF (DISCRI.LT.0.0) THEN
               WRITE(6,6000) POLY0,POLY1,POLY2,XDROIT
               IERR = 7
               GO TO 500
            ENDIF
            XROOT1 = -POLY1 + SQRT(DISCRI)
            XROOT2 = -POLY1 - SQRT(DISCRI)
            XTAU = MAX(XROOT1,XROOT2)
            IF (XTAU.LE.0.0) THEN
               WRITE(6,6000) POLY0,POLY1,POLY2,XDROIT
               IERR = 7
               GO TO 500
            ENDIF
            XTAU = XTAU/(2.*POLY2)
         ENDIF
      ELSE IF (POLY1.EQ.0.0) THEN
         IF (POLY2.EQ.0.0) THEN
            IF ((POLY0.LT.(XDROIT-EPS)).OR.(POLY0.GT.XDROIT+EPS)) THEN
               WRITE(6,6000) POLY0,POLY1,POLY2,XDROIT
               IERR = 7
               GO TO 500
            ENDIF
         ELSE
            DISCRI = XDROIT-POLY0
            IF (DISCRI.LT.0.0) THEN
               WRITE(6,6000) POLY0,POLY1,POLY2,XDROIT
               IERR = 7
               GO TO 500
            ENDIF
            XTAU = SQRT(DISCRI/POLY2)
         ENDIF
      ELSE IF (POLY2.EQ.0.0) THEN
         XTAU = (XDROIT-POLY0)/POLY1
      ELSE
         DISCRI = POLY1*POLY1 - 4.*POLY2*(POLY0-XDROIT)
         IF (DISCRI.LT.0.0) THEN
            WRITE(6,6000) POLY0,POLY1,POLY2,XDROIT
            IERR = 7
            GO TO 500
         ENDIF
         XROOT1 = -POLY1 + SQRT(DISCRI)
         XROOT2 = -POLY1 - SQRT(DISCRI)
         XTAU = MAX(XROOT1,XROOT2)
         IF (XTAU.LE.0.0) THEN
            WRITE(6,6000) POLY0,POLY1,POLY2,XDROIT
            IERR = 7
            GO TO 500
         ENDIF
         XTAU = XTAU/(2.*POLY2)
      ENDIF
*
      IF (IMPR.GE.3) THEN
         WRITE(6,5000) XTAUL,XTAUU,XTAU,POLY0,POLY1,POLY2,XDROIT,XVALC
      ENDIF
*
      IF (ABS(XTAU).GT.XTAUU) THEN
         XTAU = XTAUU
      ENDIF
*----
*  END OF THE ALGORITHM. COMPUTE THE CONTROL VARIABLES.
*----
  230 XVALC=(POLY2*XTAU+POLY1)*XTAU+POLY0
*
      IF ((IMPR.GE.3).AND.(XVALC.NE.1.0)) THEN
         WRITE(6,5000) XTAUL,XTAUU,XTAU,POLY0,POLY1,POLY2,XDROIT,XVALC
      ENDIF
*
      IF (IMPR.GE.2) THEN
         WRITE(6,3000) XTAU,XVALC
         DO 255 I=1,N,7
            II = MIN0(I+6,N)
            DO 250 J=I,II
               IF (IROW(J).LT.0) THEN
                  WRITE (ROW(J-I+1),'(1HX,I3.3)') (-IROW(J))
               ELSE
                  WRITE (ROW(J-I+1),'(1HY,I3.3)') IROW(J)
               ENDIF
*
  250       CONTINUE
            WRITE(6,4000) (ROW(J-I+1),P(J,NP3)+XTAU*P(J,NP2),J=I,II)
  255    CONTINUE
      ENDIF
      IERR = 0
*
      CALL XDDSET(XOBJ,N0,0.0D0)
      DO 280 I=1,N
         UI = U(I) + XTAU*V(I)
         IF (UI.EQ.0.0) GO TO 280
         DO 270 J=1,N0
            XOBJ(J) = XOBJ(J) - UI*APLUS(I,J)/PDG(J)
  270    CONTINUE
  280 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  500 DEALLOCATE(WRK,P)
      DEALLOCATE(V,U)
      DEALLOCATE(ICOL,IROW)
      RETURN
*
 1000 FORMAT(//,5X,'PLQUAD: FAILURE OF THE PARAMETRIC LINEAR COMPLEME',
     > 'NTARITY SOLUTION (IERR=',I5,').')
 2000 FORMAT(//,5X,'SOLUTION AFTER PIVOTING : ',
     >        /,5X,'OBJECTIVE FUNCTION     = ',1P,E12.5,
     >        /,5X,'POLY0                  = ',1P,E12.5,
     >        /,5X,'QUADRATIC CONSTRAINT   = ',1P,E12.5,
     >        /,5X,'POLY1                  = ',1P,E12.5,
     >        /,5X,'XTAU PARAMETER         = ',1P,E12.5,
     >        /,5X,'POLY2                  = ',1P,E12.5,
     >        /,5X,'CONTROL VARIABLES      = ',/,(5X,1P,10E12.4))
 3000 FORMAT(//,5X,'SOLUTION OF THE PARAMETRIC LINEAR COMPLEMENTARITY',
     >        ' PROBLEM :','*** X: KUHN-TUCKER MULTIPLIERS ;',
     >          5X,'*** Y: SLACK VARIABLES ',/,
     >        /,5X,'TAU                    = ',1P,E12.5,
     >        /,5X,'QUADRATIC CONSTRAINT   = ',1P,E12.5,/)
 4000 FORMAT(7(1X,A4,'=',E12.5),/)
 5000 FORMAT(   8X,'XTAUL',7X,'XTAUU',7X,'XTAU ',7X,
     >             'POLY0',7X,'POLY1',7X,'POLY2',7X,
     >             'XDROIT',6X,'XVALC',/,
     >          5X,1P,8E12.5)
 6000 FORMAT(   8X,'POLY0',7X,'POLY1',7X,'POLY2',7X,
     >             'XDROIT'/5X,1P,4E12.5)
      END
