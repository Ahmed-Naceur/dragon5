*DECK PLMAP1
      SUBROUTINE PLMAP1(N0,M0,APLUS,PDG,BPLUS,INPLUS,XDROIT,COUT,OBJ,
     > XOBJ,IMTHD,IMPR,IERR,BINF,BSUP,SCALE,PX,RX,DELTA,BGAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solves a linear optimization problem with quadratic constraint using
* the method of approximation programming (MAP). 
* PLMAP1 = Linear Programmation MAP1
*
*Reference:
*  R.E. Griffith and R.A. Stewart, 'A non-linear programming technique
*  for the optimization of continuous processing systems', Management
*  Science, Vol. 7, NO. 4, 379 (1961).
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
* M0      number of constraints.
* APLUS   coefficient matrix for the linear constraints.
* PDG     weights assigned to control variables in the quadratic
*         constraint.
* BPLUS   right hand sides corresponding to the coefficient matrix.
* INPLUS  constraint relations (=-1 for .GE.; =0 for .EQ.; =1 for .LE.).
* XDROIT  quadratic constraint radius squared.
* COUT    costs of control variables.
* OBJ     objective function.
* XOBJ    control variables.
* IMTHD   type of solution (=1: SIMPLEX/LEMKE; =3: MAP).
* IMPR    print flag.
*
*Parameters: ouput
* IERR    return code (=0: normal completion).
*
*Parameters: scratch
* BINF
* BSUP
* SCALE
* PX
* RX
* DELTA
* BGAR
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N0,M0,INPLUS(M0+1),IMTHD,IMPR,IERR
      DOUBLE PRECISION PDG(N0),BPLUS(M0+2),XDROIT,XOBJ(N0),BINF(N0),
     > BSUP(N0),SCALE(N0),PX(N0),RX(N0),DELTA(N0),BGAR(M0+1),
     > APLUS(M0+2,N0),COUT(N0),OBJ
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION  DELF,X,DUMY,ZMAX,XVAL,SCAL,DELX,TEMP,ERR,
     > CONT,EPSIR,EPS,EPSS
      INTEGER     ITER,ITMAX,I,J,M
      CHARACTER   CLNAME*6
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ZMAXC,BGAR0
*----
*  DATA STATEMENTS
*----
      DATA     CLNAME /'PLMAP'/
*
      EPS   = 0.0001D0
      EPSIR =EPS
      ITMAX = 100
*----
*  CONTROL-VARIABLE SCALING
*----
      DO 10 J=1,N0
         SCAL = SQRT(XDROIT/PDG(J))
         SCALE(J) = SCAL
         COUT(J) = COUT(J)*SCAL
*
         PDG(J)  = 1.0D0
         BINF(J) = 0.0D0
         BSUP(J) = 0.0D0
*
         DO 20 I=1,M0
            APLUS(I,J) = APLUS(I,J)*SCAL
   20    CONTINUE
*
  10  CONTINUE
*----
*  PRINT TABLES AFTER SCALING OF CONTROL VARIABLES
*----
      IF(IMPR.GE.5) THEN
        CALL PLNTAB(COUT,APLUS,INPLUS,BPLUS,PDG,BINF,BSUP,
     >              N0,M0,CLNAME//' AFTER SCALING OF CONTROL VARIABLES')
      ENDIF
*
      XDROIT=1.0D0
*----
*  CONSTRAINT SCALING
*----
      ALLOCATE(ZMAXC(M0))
      DO 30 I=1,M0
         ZMAX = ABS(BPLUS(I))
*
         DO 40 J=1,N0
            ZMAX = MAX(ZMAX,ABS(APLUS(I,J)))
   40    CONTINUE
         BGAR(I) = BPLUS(I)/ZMAX
*
         DO 42 J=1,N0
            APLUS(I,J) = APLUS(I,J)/ZMAX
   42    CONTINUE
         ZMAXC(I) = ZMAX
  30  CONTINUE
*----
*  COST SCALING
*----
      ZMAX = 0.0D0
      DO 45 J=1,N0
         ZMAX = MAX(ZMAX,ABS(COUT(J)))
   45 CONTINUE
      DO 50 J=1,N0
         COUT(J) = COUT(J)/ZMAX
   50 CONTINUE
*----
*  PRINT TABLES AFTER SCALING OF COSTS AND CONSTRAINTS
*----
      IF(IMPR.GE.5) THEN
        CALL PLNTAB(COUT,APLUS,INPLUS,BGAR,PDG,BINF,BSUP,N0,M0,
     >              CLNAME//' AFTER SCALING OF COSTS AND CONSTRAINTS')
      ENDIF
      ALLOCATE(BGAR0(M0))
      DO 52 I=1,M0
         BGAR0(I) = BGAR(I)
   52 CONTINUE
*----
*  INITIAL ESTIMATES
*----
      DELX = SQRT(XDROIT)
      EPSS = EPS*DELX
*
      DO 55 I=1,N0
         DELTA(I) = DELX/10.0
         RX(I)    = 0.0
   55 CONTINUE
      TEMP = DELX/SQRT(REAL(N0))/10
*----
*  MAP ITERATIONS
*----
      ITER = 0
   60 ITER = ITER + 1
      CONT = 0.0
      DO 70 I=1,M0+1
         BGAR(I) = BGAR0(I)
   70 CONTINUE
*----
*  CONTROL VARIABLE BOUNDS
*----
      DO 90 I=1,N0
         IF(ITER.EQ.1) THEN
!            XOBJ(I) = EPSIR*10.0
            XOBJ(I) = 0.0
            BINF(I) = -TEMP
            BSUP(I) = TEMP
         ELSE
            BINF(I) = -DELTA(I)
            BSUP(I) = DELTA(I)
         ENDIF
*----
*  LINEARIZATION OF THE QUADRATIC CONSTRAINT
*----
         CONT = CONT + XOBJ(I)**2
         APLUS(M0+1,I) = 2.0*XOBJ(I)
         DO 95 J=1,M0
            BGAR(J) = BGAR(J) - APLUS(J,I)*XOBJ(I)
   95    CONTINUE
  90  CONTINUE
*
      INPLUS(M0+1) = 1
      BGAR(M0+1)   = XDROIT - CONT
      M = M0 + 1
*----
*  REORGANIZE TABLES FOR SIMPLEX
*----
      DO 120 I=1,M
         DUMY = 0.0D0
*
         DO 100 J=1,N0
            DUMY = DUMY + APLUS(I,J)*BINF(J)
  100    CONTINUE
*
         BGAR(I) = BGAR(I) - DUMY
         IF(BGAR(I).GE.0.0) GOTO 120
*
         DO 110 J=1,N0
            APLUS(I,J) = -APLUS(I,J)
  110    CONTINUE
*
         BGAR(I)   = -BGAR(I)
         BGAR0(I)= -BGAR0(I)
         BPLUS(I)  = -BPLUS(I)
         INPLUS(I) = -INPLUS(I)
*
  120 CONTINUE
*
      DO 130 J=1,N0
         BSUP(J) = BSUP(J) - BINF(J)
         BINF(J) = 0.0
  130 CONTINUE
*----
*  PRINT SIMPLEX TABLES
*----
      IF(IMPR.GE.5) THEN
        CALL PLNTAB(COUT,APLUS,INPLUS,BGAR ,XOBJ  ,BINF  ,BSUP,N0,M0,
     >              CLNAME//' AFTER REORGANIZATION FOR SIMPLEX')
      ENDIF
*----
*  SOLUTION OF A LINEAR PROGRAMMING PROBLEM USING THE SIMPLEX
*----
      CALL PLSPLX(N0,M,M0+2,1,COUT,APLUS,BGAR,INPLUS,BINF,BSUP,PX,
     >          DELF,EPSS,IMTHD,IMPR,IERR)
*
      DO 140 I=1,N0
         IF(ITER.EQ.1) THEN
            PX(I) = PX(I) - TEMP
         ELSE
            PX(I) = PX(I) - DELTA(I)
         ENDIF
  140 CONTINUE
*----
*  SOLUTION OF CURRENT ITERATION
*----
      IF(IMPR.GE.2) THEN
         IF(((ITER.GE.1).AND.(IMPR.LE.2)).OR.(IMPR.GE.3)) THEN
            WRITE (6,1000)
         ENDIF
         WRITE (6,2000) ITER,DELF,(PX(I),I=1,N0)
      ENDIF
*----
*  DEGENERESCENCE OR EPS TOO SMALL
*----
      IF(IERR.EQ.1) THEN
         WRITE(6,3000) ITER
         IERR = 3
         RETURN
*----
*  NO SOLUTION IF ITER=1
*----
      ELSE IF(IERR.EQ.2) THEN
         IF(IMPR.GE.1) WRITE(6,4000) ITER
         IF(ITER.GE.ITMAX) RETURN
      ENDIF
*
      ERR = 0.0
      DO 160 I=1,N0
*
         IF((RX(I)*PX(I).LT.0.0).AND.(IERR.EQ.0)) THEN
            DELTA(I) = DELTA(I)*0.5
         ENDIF
*
         RX(I)   = PX(I)
         XOBJ(I) = XOBJ(I) + PX(I)
         ERR = ERR + PX(I)**2
  160 CONTINUE
*
      ERR  = SQRT(ERR)
      EPSS = EPS*DELX/10.0
*
      IF(IMPR.GE.1) THEN
         WRITE(6,2000) ITER,DELF,(XOBJ(I),I=1,N0)
         WRITE(6,2000) ITER,0.0,(DELTA(I),I=1,N0)
      ENDIF
*
      IF(ERR.LE.EPSS) THEN
         IERR = 0
         GOTO 170
      ENDIF
*
      IF(ITER.GE.ITMAX) THEN
         IERR = 5
         WRITE (6,5000) ITER
         RETURN
      ENDIF
      GO TO 60
*----
*  RESCALE BACK AND PRINT THE SOLUTION
*----
  170 DO 175 J=1,N0
         SCAL    = SCALE(J)
         COUT(J) = COUT(J)*ZMAX/SCAL
         XOBJ(J) = XOBJ(J)*SCAL
         PDG(J)  = XDROIT/SCAL**2
*
         DO 177 I=1,M0
            APLUS(I,J) = APLUS(I,J)/SCAL
  177    CONTINUE
  175 CONTINUE
*
      X   = 0.0D0
      OBJ = 0.0D0
      DO 180 J=1,N0
         X   = X   + PDG(J)*XOBJ(J)*XOBJ(J)
         OBJ = OBJ + XOBJ(J)*COUT(J)
  180 CONTINUE
*
      IF(IMPR.GE.1) THEN
         WRITE (6,6000) OBJ,X,(XOBJ(J),J=1,N0)
         IF(M0.GT.0) WRITE (6,7000)
*
         DO 190 I=1,M0
            XVAL = BPLUS(I)
            DO 185 J=1,N0
               XVAL = XVAL - APLUS(I,J)*XOBJ(J)*ZMAXC(I)/SCALE(J)
  185       CONTINUE
            WRITE (6,8000) I,XVAL
  190    CONTINUE
      ENDIF
      DEALLOCATE(ZMAXC,BGAR0)
      RETURN
*
1000  FORMAT(/,5X,'ITERATION',8X,'DELF',5X,'CONTROL VARIABLES')
2000  FORMAT(5X,I6,5X,8E12.4,/,(28X,5E12.4))
3000  FORMAT(5X,I6,5X,'DEGENERESCENCE OR EPS TOO SMALL')
4000  FORMAT(5X,I6,5X,'NO SOLUTION')
5000  FORMAT(5X,I6,5X,'MAXIMUM ITERATION REACHED')
6000  FORMAT(//,5X,'FINAL SOLUTION (MAP1-SIMPLEX) ',
     >        /,5X,'------------------------',
     >        /,5X,'OBJECTIVE FUNCTION    : ',1P,E12.5,
     >        /,5X,'QUADRATIC CONSTRAINT  : ',1P,E12.5,
     >        /,5X,'CONTROL VARIABLES     : ',/,(10X,10E12.4))
7000  FORMAT(//,5X,'CONSTRAINT DEVIATIONS : ',/)
8000  FORMAT(2X,I3,'...',2X,1P,D12.4)
      END
