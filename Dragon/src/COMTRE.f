*DECK COMTRE
      LOGICAL FUNCTION COMTRE(NPAR,NVP,ARB,DEB,MUPLET,IPARM,I0,II,JJ,
     1 LAST)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the index of the corresponding elementary calculation in the
* global parameter tree for a value of the tuple MUPLET. If the
* elementary calculation exists, set COMTRE=.true. otherwise, set the
* indices in the tree where the new calculation must be introduced.
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
* NPAR    number of global parameters.
* NVP     number of nodes in the global parameter tree.
* ARB     array arbval of the global parameters tree.
* DEB     array debarb of the global parameters tree.
* MUPLET  tuple of indices associated to each global parameter of the
*         elementary calculation.
*
*Parameters: output
* IPARM   index of the parameter not corresponding to a node.
* I0      index in DEB of the first element corresponding to
*         parameter iparm.
* II      index of the elementary calculation corresponding to the
*         tuple muplet (if exists). Otherwise, index in DEB of the
*         element that will contain the new elementary calculation.
* JJ      if the node has not been found, index in DEB of the
*         element corresponding to the next node.
* LAST    completion flag (=.true. if the node has not been found).
*         If LAST=.true., a node will be added at the end of the tree.
* COMTRE  index of the corresponding elementary calculation.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NPAR,NVP,ARB(NVP),DEB(NVP+1),MUPLET(NPAR),IPARM,I0,II,JJ
      LOGICAL LAST
*
      IPARM=NPAR
      II=1
      I0=1
      DO 30 IPAR=1,NPAR
         I0=DEB(I0)
         DO 10 I=DEB(II),DEB(II+1)-1
            IF(MUPLET(IPAR).EQ.ARB(I))THEN
               II=I
               GO TO 30
            ELSEIF(MUPLET(IPAR).LT.ARB(I))THEN
               JJ=I
               LAST=.FALSE.
               GO TO 20
            ENDIF
   10    CONTINUE
         JJ=DEB(II+1)
         LAST=JJ.EQ.DEB(I0)
   20    IPARM=IPAR
         COMTRE=.FALSE.
         RETURN
   30 CONTINUE
      II=DEB(II+1)
      COMTRE=.TRUE.
*
      RETURN
      END
