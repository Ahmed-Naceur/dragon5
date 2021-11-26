*DECK RENDEG
      FUNCTION RENDEG(N,LC,IM,MCU,NODE,MASK)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the degree of a node in matrix graph in MSR format.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* N       order of the system.
* LC      size of MCU.
* IM      
* MCU     connection matrices which defines the graph of the ACA matrix.
* NODE    node considered.
* MASK    mask for node to be considered in this search.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*---
* SUBROUTINE ARGUMENTS
*---
      INTEGER RENDEG
      INTEGER N,LC,IM(N+1),MCU(LC),NODE,MASK(N)
*---
* LOCAL VARIABLES
*---
      INTEGER J,NJ
*
      RENDEG=0 
      DO J=IM(NODE)+1,IM(NODE+1) 
         NJ=MCU(J)
         IF (NJ.GT.0) THEN
            IF (MASK(NJ).EQ.1) RENDEG=RENDEG+1
         ENDIF
      ENDDO
*
      RETURN
      END
