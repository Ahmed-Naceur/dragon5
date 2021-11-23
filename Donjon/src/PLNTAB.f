*DECK PLNTAB
      SUBROUTINE PLNTAB(GF,APLUS,INPLUS,BPLUS,XITK,XINF,XSUP,NDEC,M0,
     >                  SRCNAM)
*----------------------------------------------------------------------*
*                                                                      *
*Purpose:
* Print the arrays of the linear optimization problem. 
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* R. Chambon
*
*Parameters: input
* GF      costs of control variables.
* APLUS   coefficient matrix for the linear constraints.
* INPLUS  constraint relations (=-1 for .GE.; =0 for .EQ.; =1 for .LE.).
* BPLUS   right hand sides corresponding to the coefficient matrix.
* XITK    weights assigned to control variables in the quadratic
*         constraint.
* XINF    lower bounds of control variables.
* XSUP    upper bounds of control variables.
* NDEC    number of control variables.
* M0      number of constraints plus the number of lower/upper bounds
*         intercepting the quadratic constraint.
* SRCNAM  character text to print.
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER       M0
      DOUBLE PRECISION BPLUS(M0+2),XITK(NDEC),XINF(NDEC),XSUP(NDEC),
     >                 GF(NDEC),APLUS(M0+2,NDEC)
      INTEGER       INPLUS(M0+1)
      CHARACTER*(*) SRCNAM
*----
*  LOCAL VARIABLES
*----
      CHARACTER*2   CTYPES(-1:1)
      CHARACTER*80  FMT
*
      DATA CTYPES / '>=',' =','<=' /
*
      WRITE(6,1000) SRCNAM
*
      IF (NDEC.GT.8) THEN
         RETURN
      ELSE
         FMT = '(1P,XXE13.5,5X,A3,5X,1P,E13.5)'
         NVAL = NDEC
      ENDIF
*
      IDX = INDEX(FMT,'X')
      WRITE(FMT(IDX:IDX+1),'(I2.2)') NVAL
*----
*  PRINT CONTROL-VARIABLE COSTS
*----
      WRITE(6,2000) (I,I=1,NDEC)
      WRITE(6,3000) (GF(I),I=1,NDEC)
*----
*  PRINT COEFFICIENT MATRIX
*----
      IF(M0.GT.0) THEN
        WRITE(6,4000)
        DO 10 J=1,M0
          WRITE(6,FMT) (APLUS(J,I),I=1,NDEC),CTYPES(INPLUS(J)),BPLUS(J)
  10    CONTINUE
      ENDIF
*
      WRITE(6,5000) (XINF(I),I=1,NDEC)
      WRITE(6,6000) (XSUP(I),I=1,NDEC)
      WRITE(6,7000) (XITK(I),I=1,NDEC)
      RETURN
*
1000  FORMAT(//,5X,'PRINT LINEARIZED OPTIMIZATION PROBLEM IN ',A,/)
2000  FORMAT( /,5X,'COST(NDEC)',//,(10(5X,I3,5X)),//)
3000  FORMAT((1P,10E13.5))
4000  FORMAT( /,5X,'APLUS(M0,NDEC)',35X,'INPLUS(M0)',35X,'BPLUS(M0)',/)
5000  FORMAT( /,5X,'XINF(NDEC)  ',//,(1P,10E13.5))
6000  FORMAT( /,5X,'XSUP(NDEC)  ',//,(1P,10E13.5))
7000  FORMAT( /,5X,'WEIGHT(NDEC)',//,(1P,10E13.5))
      END
