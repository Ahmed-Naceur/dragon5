*DECK RENLST
      SUBROUTINE RENLST(N,LC,NFIRST,IM,MCU,TYPOR,NLEV,LEV,LEVPT,MASK)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Level-set traversal method of the graph of a matrix stored
* in MSR format.
*
*Reference
* Y. Saad, "Iterative Methods for Sparse Linear Systems",
* PWS Publishing Company, Boston, 1996
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
* NFIRST  starting node.
* IM      
* MCU     connection matrices which defines the graph of the ACA matrix.
* TYPOR   type of level traversal
*         0 : Breadth First Search
*         1 : Cuthill-McKee ordering
*
* Parameters: output
* NLEV    number of level in the last level-set traversal.
* LEV
* LEVPT   level data structure of the last level-set traversal.
*         LEV(LEVPT(I):LEVPT(I+1)-1) : nodes in level i.
* MASK    mask for node to be considered in this search.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*---
* SUBROUTINE ARGUMENTS
*---
      INTEGER N,LC,NFIRST,IM(N+1),MCU(LC),TYPOR,NLEV,LEV(N),
     1 LEVPT(N+1),MASK(N)
*---
* LOCAL VARIABLES
*---
      INTEGER IDEB,IEND,NEWEND,I,NODE,J,NJ,RENDEG
      INTEGER, DIMENSION(:), ALLOCATABLE :: DEG
*
      ALLOCATE(DEG(N))
      CALL XDISET(MASK,N,1)
      NLEV=1
      IDEB=1
      IEND=1
      LEVPT(NLEV)=IDEB
      LEV(1)=NFIRST
      MASK(NFIRST)=0
*
      DO WHILE (IEND.LT.N)
*     visit neighboring nodes of nodes LEV(IDEB^in:IEND^in)  
         NEWEND=IEND
         IF (TYPOR.EQ.1) THEN
*        Cuthill McKee ordering
*        find the degrees for this level
            DO I=IDEB,IEND
               NODE=LEV(I)
               DEG(I-IDEB+1)=RENDEG(N,LC,IM,MCU,NODE,MASK)
            ENDDO
*           sort this level by increasing degrees
            CALL RENINS((IEND-IDEB+1),LEV(IDEB),DEG)
         ENDIF
         DO I=IDEB,IEND
            NODE=LEV(I)
            DO J=IM(NODE)+1,IM(NODE+1)
               NJ=MCU(J)
               if (NJ.GT.0) THEN
                  if (MASK(NJ).EQ.1) THEN
                     NEWEND=NEWEND+1
                     MASK(NJ)=0
                     LEV(NEWEND)=NJ
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         IF (NEWEND.EQ.IEND) 
     1       CALL XABORT('RENLST: INCOHERENT MATRIX GRAPH')
         IDEB=IEND+1
         IEND=NEWEND
*        unmarked neighbors are added in LEV(IDEB^out:IEND^out) 
*        where IDEB^out=IEND^in + 1 
*              IEND^out=IEND^in + number of unmarked neighbors found
*        start new level
         NLEV=NLEV+1
         LEVPT(NLEV)=IEND+1
      ENDDO
      NLEV=NLEV-1
*
      DEALLOCATE(DEG)
*
      RETURN
      END
