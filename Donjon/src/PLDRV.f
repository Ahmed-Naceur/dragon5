*DECK PLDRV
      SUBROUTINE PLDRV(IPOPT,N0,NCST,M0,MINMAX,IMTHD,FCOST,XOBJ,PDG,
     >           GRAD,INEGAL,CONTR,DINF,DSUP,XDROIT,EPSIM,IMPR,IERR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Prepares the different matrices for the resolution of a linear
* optimisation problem with a quadratic constraint. The different
* available methods are the MAP technique, the lemke, the augmented-
* Lagrangian and the penalty function.
* PLDRV = Linear Programmation DRiVer for options
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
* IPOPT   pointer to the L_OPTIMIZE object.
* N0      number of control variables.
* NCST    number of constraints.
* M0      number of constraints plus the number of lower/upper bounds
*         intercepting the quadratic constraint.
* MINMAX  type of optimization (=-1: minimize; =1: maximize).
* IMTHD   type of solution (=1: SIMPLEX/LEMKE; =2: LEMKE/LEMKE;
*         =3: MAP; =4: Augmented Lagragian; =5: External penalty
*         funnction).
* FCOST   objective function.
* XOBJ    control variables.
* PDG     weights assigned to control variables in the quadratic
*         constraint.
* GRAD    linearized gradients (GRAD(:,1) are control variable costs
*         and GRAD(:,2:NCST+1) are linear constraint coefficients).
* INEGAL  constraint relations (=-1 for .GE.; =0 for .EQ.; =1 for .LE.).
* CONTR   constraint right hand sides.
* DINF    lower bounds of control variables.
* DSUP    upper bounds of control variables.
* XDROIT  quadratic constraint radius squared.
* EPSIM   tolerence used for inner linear LEMKE or SIMPLEX calculation.
* IMPR    print flag.
*
*Parameters: ouput
* IERR    return code (=0: normal completion).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPOPT
      INTEGER N0,NCST,M0,MINMAX,IMTHD,INEGAL(NCST),IMPR,IERR
      DOUBLE PRECISION FCOST,XOBJ(N0),PDG(N0),GRAD(N0,NCST+1),
     >     CONTR(NCST),DINF(N0),DSUP(N0),XDROIT,EPSIM
*----
*  LOCAL VARIABLES
*----
      INTEGER      ME,MI,I,J,NPM
      DOUBLE PRECISION  XX
      CHARACTER    CLNAME*6
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INPLUS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BPLUS,GF
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: APLUS
*----
*  DATA STATEMENTS
*----
      DATA         CLNAME /'PLDRV'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INPLUS(M0+1))
      ALLOCATE(BPLUS(M0+2),GF(N0),APLUS(M0+2,(M0+1)+N0))
*----
*  SET COSTS OF CONTROL VARIABLES
*----
      DO 10 I=1,N0
         GF(I)=GRAD(I,1)*REAL(MINMAX)
   10 CONTINUE
*----
*  ORGANIZE TABLES APLUS AND BPLUS FOR EQUALITY CONSTRAINTS
*----
      ME=0
      DO 30 I=1,NCST
         IF(INEGAL(I).EQ.0) THEN
            ME = ME + 1
            DO 20 J=1,N0
               APLUS(ME,J) = GRAD(J,I+1)
   20       CONTINUE
            BPLUS(ME)  = CONTR(I)
            INPLUS(ME) = 0
         ENDIF
   30 CONTINUE
*----
*  ORGANIZE TABLES APLUS AND BPLUS FOR INEQUALITY CONSTRAINTS
*----
      MI=0
      DO 50 I=1,NCST
         IF(INEGAL(I).NE.0) THEN
            MI = MI + 1
            DO 40 J=1,N0
               APLUS(ME+MI,J) = GRAD(J,I+1)
   40       CONTINUE
            BPLUS(ME+MI)  = CONTR(I)
            INPLUS(ME+MI) = INEGAL(I)
         ENDIF
   50 CONTINUE
*----
*  ORGANIZE TABLES APLUS AND BPLUS FOR CONTROL-VARIABLE BOUNDS
*----
      DO 80 I=1,N0
         XX = SQRT(XDROIT/PDG(I))
         IF(DINF(I).GT.-XX) THEN
            MI = MI + 1
            DO 60 J=1,N0
               APLUS(ME+MI,J) = 0.0D0
   60       CONTINUE
            APLUS(ME+MI,I) = 1.0D0
            BPLUS(ME+MI)   = DINF(I)
            INPLUS(ME+MI)  = -1
         ENDIF
         IF(DSUP(I).LT.XX) THEN
            MI = MI + 1
            DO 70 J=1,N0
               APLUS(ME+MI,J) = 0.0D0
   70       CONTINUE
            APLUS(ME+MI,I) = 1.0D0
            BPLUS(ME+MI)   = DSUP(I)
            INPLUS(ME+MI)  = 1
         ENDIF
   80 CONTINUE
*
      DO 90 J=1,N0
         APLUS(M0+1,J) = 0.0D0
   90 CONTINUE
      BPLUS(M0+1)   = 0.0D0
      INPLUS(M0+1)  = 0
*
      IF(M0.NE.ME+MI) THEN
         WRITE (6,1000) M0,ME,MI
         CALL XABORT('PLDRV: M0 AND ME+MI ARE NOT THE SAME')
      ENDIF
*----
*  PRINT THE QUASILINEAR PROBLEM
*----
      IF(IMPR.GE.5) THEN
         CALL PLNTAB(GF,APLUS,INPLUS,BPLUS,PDG,DINF,DSUP,N0,M0,
     >             CLNAME)
      ENDIF
*----
*  LEMKE METHOD
*----
      NPM=(M0+1)+N0
      IF(IMTHD.LE.2) THEN
         CALL PLMAP2(N0,M0,APLUS,PDG,BPLUS,INPLUS,XDROIT,GF,FCOST,XOBJ,
     1   IMTHD,EPSIM,IMPR,IERR)
*----
*  MAP
*----
      ELSE IF(IMTHD.EQ.3) THEN
         CALL PLMAP1(N0,M0,APLUS,PDG,BPLUS,INPLUS,XDROIT,GF,FCOST,
     1   XOBJ,IMTHD,IMPR,IERR)
*----
*  AUGMENTED LAGRANGIAN
*----
      ELSE IF(IMTHD.EQ.4) THEN
         CALL PLLA(IPOPT,N0,M0,APLUS,PDG,BPLUS,INPLUS,XDROIT,
     1   FCOST,GF,XOBJ,IMPR,EPSIM,NCST,GRAD,CONTR,INEGAL,IERR)
*----
*  EXTERNAL PENALTY METHOD
*----
      ELSE IF(IMTHD.EQ.5) THEN
         CALL PLPNLT(IPOPT,N0,M0,APLUS,PDG,BPLUS,INPLUS,XDROIT,
     1   FCOST,GF,XOBJ,IMPR,EPSIM,NCST,GRAD,CONTR,INEGAL,IERR)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(APLUS,GF,BPLUS)
      DEALLOCATE(INPLUS)
      RETURN
*
1000  FORMAT(/' PLDRV: INCONSISTENCY BETWEEN M0 AND ME+MI ',3I5)
      END
