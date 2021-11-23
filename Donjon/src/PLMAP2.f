*DECK PLMAP2
      SUBROUTINE PLMAP2(N0,M0,APLUS,PDG,BPLUS,INPLUS,XDROIT,COUT,OBJ,
     > XOBJ,IMTHD,EPSIM,IMPR,IERR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solves a linear optimization problem with quadratic constraints using
* the method of LEMKE. 
* PLMAP2 = Linear Programmation MAP2
*
*Reference:
*    J. A. Ferland, 'A linear programming problem with an additional
*    quadratic constraint solved by parametric linear complementarity',
*    Publication number 497, Departement d'informatique et de recherche
*    operationnelle, Universite de Montreal, January 1984.
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
* IMTHD   type of solution (=1: SIMPLEX/LEMKE; =2: LEMKE/LEMKE).
* EPSIM   tolerence used for inner linear SIMPLEX calculation.
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
      INTEGER N0,M0,INPLUS(M0+1),IMTHD,IMPR,IERR
      DOUBLE PRECISION PDG(N0),BPLUS(M0+2),XDROIT,XOBJ(N0),EPSIM,
     > APLUS(M0+2,(M0+1)+N0),COUT(N0),OBJ
*----
*  LOCAL VARIABLES
*----
      CHARACTER   CLNAME*6
      DOUBLE PRECISION  X,ZMAX,XVAL,SCAL,EPS,FACTOR
      INTEGER     I,J,M0NEW
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BINF,BSUP,SCALE
*----
*  DATA STATEMENTS
*----
      DATA CLNAME /'PLMAP2'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(BINF(N0),BSUP(N0),SCALE(N0))
*
      EPS=EPSIM
*----
*  CONTROL-VARIABLE SCALING
*----
      FACTOR=MAX(XDROIT,EPSIM)
      DO 10 J=1,N0
         SCAL = SQRT(FACTOR/PDG(J))
         SCALE(J) = SCAL
         COUT(J) = COUT(J)*SCAL
*
         PDG(J)  = 1.0
         BINF(J) = 0.0
         BSUP(J) = 0.0
*
         DO 20 I=1,M0
            APLUS(I,J) = APLUS(I,J)*SCAL
   20    CONTINUE
*
   10 CONTINUE
*----
*  PRINT TABLES AFTER SCALING OF CONTROL VARIABLES
*----
      IF(IMPR.GE.5) THEN
         CALL PLNTAB(COUT,APLUS,INPLUS,BPLUS,PDG,BINF,BSUP,
     >               N0,M0,CLNAME//' AFTER SCALING OF VARIABLES')
      ENDIF
*
      XDROIT = XDROIT/FACTOR
*----
*  CONSTRAINT SCALING
*----
      DO 30 I=1,M0
         ZMAX = ABS(BPLUS(I))
*
         DO 40 J=1,N0
            ZMAX = MAX(ZMAX,ABS(APLUS(I,J)))
   40    CONTINUE
         BPLUS(I) = BPLUS(I)/ZMAX
*
         DO 42 J=1,N0
            APLUS(I,J) = APLUS(I,J)/ZMAX
   42    CONTINUE
   30 CONTINUE
*----
*  COST SCALING
*----
      ZMAX = 0.0D0
      DO 45 J=1,N0
         ZMAX = MAX(ZMAX,ABS(COUT(J)))
   45 CONTINUE
*
      DO 50 J=1,N0
         COUT(J) = COUT(J)/ZMAX
   50 CONTINUE
*----
*  STEP 1
*----
      M0NEW = M0 + 1
      DO 55 I=1,N0
         BINF(I) = -SQRT(XDROIT)
         BSUP(I) =  SQRT(XDROIT)
         APLUS(M0NEW,I) = 0.0D0
   55 CONTINUE
      BPLUS(M0NEW) = 0.0D0
*----
*  PRINT TABLES AFTER SCALING OF COSTS AND CONSTRAINTS
*----
      IF(IMPR.GE.5) THEN
         CALL PLNTAB(COUT,APLUS,INPLUS,BPLUS,PDG,BINF,BSUP,N0,M0,
     >               CLNAME//' AFTER SCALING OF COSTS AND CONSTRAINTS')
      ENDIF
*
      IF(IMTHD.EQ.1) THEN
*----
*  SOLUTION OF A LINEAR OPTIMIZATION PROBLEM USING THE SIMPLEX METHOD
*----
         CALL PLSPLX(N0,M0,M0+2,1,COUT,APLUS,BPLUS,INPLUS,BINF,BSUP,
     >             XOBJ,OBJ,EPS,IMTHD,IMPR,IERR)
*
         DO 70 I=1,M0
            IF(INPLUS(I).EQ.-1) THEN
               DO 60 J=1,N0
                  APLUS(I,J) = -APLUS(I,J)
   60          CONTINUE
               BPLUS(I)   = -BPLUS(I)
               INPLUS(I)  = 1
            ELSE IF(INPLUS(I).EQ.0) THEN
               DO 65 J=1,N0
                  APLUS(M0NEW,J) = APLUS(M0NEW,J) - APLUS(I,J)
   65          CONTINUE
               BPLUS(M0NEW)   = BPLUS(M0NEW)   - BPLUS(I)
            ENDIF
   70    CONTINUE
      ELSE
*----
*  SOLUTION OF A LINEAR OPTIMIZATION PROBLEM USING THE LINEAR LEMKE
*  METHOD
*----
         DO 90 I=1,M0
            IF(INPLUS(I).EQ.-1) THEN
               DO 75 J=1,N0
                  APLUS(I,J) = -APLUS(I,J)
   75          CONTINUE
               BPLUS(I)  = -BPLUS(I)
               INPLUS(I) = 1
            ELSE IF(INPLUS(I).EQ.0) THEN
               DO 80 J=1,N0
                  APLUS(M0NEW,J) = APLUS(M0NEW,J) - APLUS(I,J)
   80          CONTINUE
               BPLUS(M0NEW) = BPLUS(M0NEW) - BPLUS(I)
            ENDIF
   90    CONTINUE
         CALL PLLINR(N0,M0NEW,COUT,APLUS,BPLUS,BINF,BSUP,XOBJ,EPS,IMPR,
     >   IERR)
      ENDIF
*
      IF(IERR.GE.1) THEN
         WRITE (6,6000) IERR
         GO TO 500
      ENDIF
*
      X = 0.0D0
      OBJ = 0.0D0
      DO 100 J=1,N0
         X = X + PDG(J)*XOBJ(J)*XOBJ(J)
         OBJ = OBJ + COUT(J)*XOBJ(J)
  100 CONTINUE
      IF(IMPR.GE.2) THEN
         IF(IMTHD.EQ.1) THEN
            WRITE (6,1000) OBJ,X,(XOBJ(I),I=1,N0)
         ELSE IF(IMTHD.EQ.2) THEN
            WRITE (6,1500) OBJ,X,(XOBJ(I),I=1,N0)
         ENDIF
      ENDIF
*
      IF(IMPR.GE.5) THEN
         WRITE(6,*) 'AFTER LINEAR OPTIMIZATION'
         WRITE(6,*) 'XOBJ   ',(XOBJ(J),J=1,N0)
         WRITE(6,*) 'PDG    ',(PDG(J),J=1,N0)
         WRITE(6,*) 'OBJ    ',OBJ
         WRITE(6,*) 'X      ',X
         WRITE(6,*) 'XDROIT ',XDROIT
      ENDIF
*----
*  SOLUTION OF A LINEAR OPTIMIZATION PROBLEM WITH A QUADRATIC CONSTRAINT
*  USING THE GENERAL LEMKE METHOD
*----
      IF(X.GT.XDROIT) THEN
        DO J=1,N0
          APLUS(M0NEW+1,J) = COUT(J)
        ENDDO
        BPLUS(M0NEW+1) = OBJ
*
        CALL PLQUAD(N0,M0NEW,APLUS,BPLUS,PDG,XDROIT,COUT,XOBJ,EPS,IMPR,
     >  IERR)
*
        IF(IERR.GE.1) THEN
          WRITE(6,2000) IERR
          IERR = IERR + 10
          GO TO 500
        ENDIF
      ENDIF
*----
*  RESCALE BACK AND PRINT THE SOLUTION
*----
      DO 170 J=1,N0
         SCAL = SCALE(J)
         COUT(J) = COUT(J)*ZMAX/SCAL
         XOBJ(J) = XOBJ(J)*SCAL
         PDG(J)  = FACTOR/SCAL**2
*
         DO 175 I=1,M0
            APLUS(I,J) = APLUS(I,J)/SCAL
  175    CONTINUE
  170 CONTINUE
*----
*  COMPUTE THE NEW OPTIMAL POINT
*----
      X   = 0.0D0
      OBJ = 0.0D0
      DO 180 J=1,N0
         X   = X   + PDG(J)*XOBJ(J)*XOBJ(J)
         OBJ = OBJ + XOBJ(J)*COUT(J)
  180 CONTINUE
*
      IF(IMPR.GE.1) THEN
         WRITE (6,3000) OBJ,X,(XOBJ(J),J=1,N0)
         WRITE (6,4000)
*
         DO 190 I=1,M0
            XVAL = BPLUS(I)
            DO 185 J=1,N0
               XVAL = XVAL - APLUS(I,J)*XOBJ(J)
  185       CONTINUE
            WRITE (6,5000) I,XVAL
  190    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  500 DEALLOCATE(SCALE,BSUP,BINF)
      RETURN
*
1000  FORMAT(//,5X,'SOLUTION WITHOUT QUADRATIC CONSTRAINT (SIMPLEX) :',
     >        /,5X,'------------------------------------------------',
     >        /,5X,'OBJECTIVE FUNCTION    : ',1P,E12.5,
     >        /,5X,'QUADRATIC CONSTRAINT  : ',1P,E12.5,
     >        /,5X,'CONTROL VARIABLES     : ',/,(10X,10E12.4))
1500  FORMAT(//,5X,'SOLUTION WITHOUT QUADRATIC CONSTRAINT (LINR) :',
     >        /,5X,'---------------------------------------------',
     >        /,5X,'OBJECTIVE FUNCTION    : ',1P,E12.5,
     >        /,5X,'QUADRATIC CONSTRAINT  : ',1P,E12.5,
     >        /,5X,'CONTROL VARIABLES     : ',/,(10X,10E12.4))
2000  FORMAT(//,5X,'PLMAP2: ECHEC DU MODULE QUADR  IERR = ',I2)
3000  FORMAT(//,5X,'FINAL SOLUTION        :',
     >        /,5X,'---------------------',
     >        /,5X,'OBJECTIVE FUNCTION    : ',1P,E12.5,
     >        /,5X,'QUADRATIC CONSTRAINT  : ',1P,E12.5,
     >        /,5X,'CONTROL VARIABLES     : ',/,(10X,10E12.4))
4000  FORMAT(//,5X,'CONSTRAINT DEVIATIONS :',/)
5000  FORMAT(2X,I3,'...',2X,1P,E12.4)
6000  FORMAT(//,5X,'PLMAP2: FAILURE OF LINEAR ALGORITHM (IERR=',I5,')')
      END
