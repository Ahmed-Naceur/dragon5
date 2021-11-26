*DECK RENUM
      SUBROUTINE RENUM(N,LC,NFIRST,IM,MCU,TYPOR1,TYPOR2,NLEV,LEV,LEVPT,
     1           IPERM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Renumbering of the unknowns for ilu0 decomposition of the a matrix.
* level-set traversal method which starts by finding a pseudo-peripheral
* node of the graph of the matrix.
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
* IM      
* MCU     connection matrices which defines the graph of the ACA matrix.
* TYPOR1  type of level traversal for pseudo-peripheral node.
* TYPOR2  type of level traversal for the last level-set traversal.
*         0 : Breadth First Search
*         1 : Cuthill-McKee ordering
*
*Parameters: input/output
* NFIRST  starting node for the traversal.
*
* Parameters: output
* NLEV    number of level in the last level-set traversal.
* LEV
* LEVPT   level data structure of the last level-set traversal.
*         LEV(LEVPT(I):LEVPT(I+1)-1) : nodes in level i.
* IPERM   permutation array: IPERM(I) : new index of node I.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*---
* SUBROUTINE ARGUMENTS
*---
      INTEGER N,LC,NFIRST,IM(N+1),MCU(LC),TYPOR1,TYPOR2,NLEV,LEV(N),
     1 LEVPT(N+1),IPERM(N)
*---
* LOCAL VARIABLES
*---
      INTEGER DELTA,J,MINDEG,DEG,PPN
      INTEGER RENDEG
      LOGICAL NOTDONE
*---
* PSEUDO-PERIPHERAL NODE SEARCH (IPERM IS USED AS A WORK ARRAY)
*---
      DELTA=0
      NOTDONE=.TRUE.
      DO WHILE (NOTDONE)
*     Level-set traversal (Breadth First Search or Cuthill Mc-Kee) from node NFIRST
      CALL RENLST(N,LC,NFIRST,IM,MCU,TYPOR1,NLEV,LEV,LEVPT,IPERM)
      IF (NLEV.GT.DELTA) THEN
         MINDEG=N+1
*        scan last level of the previous expansion to find a node (PPN) 
*        with minimum degree
         DO J=LEVPT(NLEV),LEVPT(NLEV+1)-1
            PPN=LEV(J)
            DEG=RENDEG(N,LC,IM,MCU,PPN,IPERM)
            IF (DEG.LT.MINDEG) THEN
                MINDEG=DEG
                NFIRST=PPN
            ENDIF
         ENDDO
        DELTA=NLEV
      ELSE
        NOTDONE=.FALSE.
      ENDIF
      ENDDO
*---
* LEVEL-SET TRAVERSAL FROM NODE NFIRST, A PSEUDO-PERIPHERAL NODE
*---
      IF (TYPOR1.NE.TYPOR2) THEN
*     Level-set traversal (Breadth First Search or Cuthill Mc-Kee) from node NFIRST
         CALL RENLST(N,LC,NFIRST,IM,MCU,TYPOR2,NLEV,LEV,LEVPT,IPERM)
      ENDIF
*---
* FORM IPERM ARRAY BY REVERSING THE ORDERING DEFINED BY THE LEV ARRAY
*---
      DO J=1,N
         IPERM(J)=LEV(N-J+1)
      ENDDO
*
      RETURN
      END
