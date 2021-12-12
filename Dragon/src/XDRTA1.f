*DECK XDRTA1
      SUBROUTINE XDRTA1(IPTRK,LBIC,LEXP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the tabulated functions required by the flux solution.
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
* IPTRK   pointer to the tracking information (L_TRACK signature).
* LBIC    compute the bickley function tables.
* LEXP    compute the bickley function tables.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPTRK
      LOGICAL LBIC,LEXP
*----
*  VARIABLES FOR BICKLEY FUNCTIONS
*----
      DOUBLE PRECISION DX,DEX
      PARAMETER (DX=0.02D0,NBX=600,NBICK=5,DEX=1.D0/512.D0,NBEX=7936)
      INTEGER MLOG(NBICK)
      DATA (MLOG(JJ),JJ=1,NBICK)/30,15,0,0,0/
*----
*  INITIALIZE COMMON BLOCKS FOR BICKLEY FUNCTIONS
*----
      IF(LBIC) THEN
         CALL XDRKIN(IPTRK,DX,NBX,MLOG)
      ENDIF
*----
*  INITIALIZE COMMON BLOCK FOR EXPONENTIAL FUNCTION
*----
      IF(LEXP) THEN
         CALL XDREXP(IPTRK,DEX,NBEX)
      ENDIF
*
      RETURN
      END
