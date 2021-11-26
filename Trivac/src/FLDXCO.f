*DECK FLDXCO
      SUBROUTINE FLDXCO(IPFLUX,L4,NUN,VECT,LMPR,B)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compare two solutions and print the logarithm of error.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPFLUX  L_FLUX pointer to the solution.
* L4      order of matrix systems.
* NUN     number of unknowns in each energy group.
* VECT    unknown vector.
* LMPR    logarithm print flag (.true. to print the logarithm value).
*
*Parameters: output
* B      base 10 logarithm of the error.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPFLUX
      INTEGER L4,NUN
      REAL VECT(NUN),B
      LOGICAL LMPR
*----
*  LOCAL VARIABLES
*----
      REAL, DIMENSION(:), ALLOCATABLE :: REF
*
      CALL LCMLEN(IPFLUX,'REF',ILONG,ITYLCM)
      IF(ILONG.EQ.0) RETURN
      IF(ILONG.NE.NUN) CALL XABORT('FLDXCO: INVALID LENGTH FOR REF.')
      ALLOCATE(REF(ILONG))
      CALL LCMGET(IPFLUX,'REF',REF)
      IN=0
      ERR1=0.0
      DO 5 I=1,L4
      IF(ABS(REF(I)).GT.ERR1) THEN
         IN=I
         ERR1=ABS(REF(I))
      ENDIF
    5 CONTINUE
      WEIGHT=REF(IN)/VECT(IN)
      ERR2=0.0
      DO 10 I=1,L4
      ERR2=AMAX1(ERR2,ABS(REF(I)-VECT(I)*WEIGHT))
   10 CONTINUE
      DEALLOCATE(REF)
      A=ERR2/ERR1
      IF(A.GT.0.0) THEN
         B=LOG10(A)
      ELSE
         B=-5.0
      ENDIF
      IF(LMPR) WRITE (6,20) A,B
      RETURN
*
   20 FORMAT (7H ERROR=,1P,E10.2,5X,11HLOG(ERROR)=,E10.2)
      END
