*DECK ALPLGN
      DOUBLE PRECISION FUNCTION ALPLGN(L,M,X)
*
*-----------------------------------------------------------------------
*
*Purpose:
* return the Ferrer definition of the associated Legendre function.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* L       main index
* M       secondary index
* X       direction cosine
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER L,M
      DOUBLE PRECISION X
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPX
*
      IF(M.LT.0) THEN
        CALL XABORT('ALPLGN: BAD ARGUMENTS (1)')
      ELSE IF(M.GT.L) THEN
        CALL XABORT('ALPLGN: BAD ARGUMENTS (2)')
      ELSE IF(ABS(X).GT.1.0D0) THEN
        CALL XABORT('ALPLGN: BAD ARGUMENTS (3)')
      ENDIF
      PMM=1.0D0
      IF(M.GT.0) THEN
        ALLOCATE(IPX(2*M))
        DO I=1,2*M
          IPX(I)=I
        ENDDO
        PMM=PRODUCT(IPX,MASK=MOD(IPX,2)==1)*SQRT((1.0D0-X)*(1.0D0+X))**M
        DEALLOCATE(IPX)
      ENDIF
      IF(L.EQ.M) THEN
        ALPLGN=PMM
      ELSE
        PMMP1=(2*M+1)*X*PMM
        IF(L.EQ.M+1) THEN
          ALPLGN=PMMP1
        ELSE
          PLL=0.0D0
          DO LL=M+2,L
            PLL=((2*LL-1)*X*PMMP1-(LL+M-1)*PMM)/(LL-M)
            PMM=PMMP1
            PMMP1=PLL
          ENDDO
          ALPLGN=PLL
        ENDIF
      ENDIF
      RETURN
      END
