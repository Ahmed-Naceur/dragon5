*DECK COMARB
      SUBROUTINE COMARB(NPAR,NVPO,NVPN,OLDDEB,OLDARB,LGNEW,MUPLET,NCAL,
     1 NEWDEB,NEWARB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Add a node to the parameter tree.
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
* NPAR    number of parameters.
* NVPO    original number of nodes in the parameter tree.
* NVPN    new number of nodes in the parameter tree.
* OLDDEB  original array DEBARB of the parameters tree.
* OLDARB  original array ARBVAL of the parameters tree.
* LGNEW   new parameter flag (=.true. if the I-th parameter has changed
*         in the new elementary calculation).
* MUPLET  tuple of indices associated to each parameter of the
*         elementary calculation.
*
*Parameters: input/output
* NCAL    index of the last elementary calculation on input and
*         index of the new elementary calculation at output (value
*         is incremented by 1).
*
*Parameters: output
* NEWDEB  new array DEBARB of the parameters tree.
* NEWARB  new array ARBVAL of the parameters tree.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      LOGICAL LGNEW(NPAR)
      INTEGER NPAR,NVPO,NVPN,OLDDEB(NVPO+1),OLDARB(NVPO),MUPLET(NPAR),
     1 NCAL,NEWDEB(NVPN+1),NEWARB(NVPN)
*----
*  LOCAL VARIABLES
*----
      LOGICAL LAST,DUMMY,COMTRE
*----
*  Change addresses of parameter values if new value added.
*----
      I0=1
      DO 10 IPAR=1,NPAR
         I0=OLDDEB(I0)
         I1=OLDDEB(I0)-1
         IF(LGNEW(IPAR))THEN
            DO 5 I=I0,I1
               IF(OLDARB(I).GE.MUPLET(IPAR))OLDARB(I)=OLDARB(I)+1
    5       CONTINUE
         ENDIF
   10 CONTINUE
*----
*  Find point where new branch is to be added and copy the
*  unchanged part.
*----
      DUMMY=COMTRE(NPAR,NVPO,OLDARB,OLDDEB,MUPLET,ISTART,I0,II,JJ,LAST)
*
      I1=OLDDEB(I0)-1
      J0=II
      JX=JJ-OLDDEB(J0)+1
*
      DO 15 I=1,J0
         NEWDEB(I)=OLDDEB(I)
  15  CONTINUE
      DO 20 I=1,I0-1
         NEWARB(I)=OLDARB(I)
  20  CONTINUE
*----
*  Modified addresses, shifted in the array, and those of
*  inserted branch.
*  Computation of the address where the part of array with calc.
*  identifiers starts.
*----
      INCR=0
      DO 35 I=J0+1,NVPO-NCAL
         IF(I.EQ.JJ)THEN
            NEWDEB(I+INCR)=OLDDEB(I)+INCR+1
            INCR=INCR+1
            JX=JJ
            JJ=OLDDEB(JJ)
         ENDIF
         NEWDEB(I+INCR)=OLDDEB(I)+INCR+1
   35 CONTINUE
*----
*  Especial treatement if new added point is the rightmost point
*  in the tree.
*----
      IF(LAST)THEN
         IF(ISTART.LT.NPAR)THEN
            NEWDEB(NVPO+1-NCAL+INCR)=OLDDEB(NVPO+1-NCAL)+INCR+1
            INCR=INCR+1
         ENDIF
         JJ=NVPO+2
      ELSE
         IF(ISTART.EQ.NPAR)THEN
            JJ=OLDDEB(J0)+JX
         ELSE
            JJ=OLDDEB(JX)+1
         ENDIF
      ENDIF
*----
*  Address of next nonexisting point used do get dimension at the end
*----
      NEWDEB(NVPO+1-NCAL+INCR)=OLDDEB(NVPO+1-NCAL)+INCR+1
*----
*  Part of the NEWDEB array containing calculation numbers.
*----
      DO 37 I=NVPO+2-NCAL,JJ-1
         NEWDEB(I+INCR)=OLDDEB(I)
   37 CONTINUE
      NCAL=NCAL+1
      NEWDEB(JJ+INCR)=NCAL
      INCR=INCR+1
      DO 39 I=JJ,NVPO+1
         NEWDEB(I+INCR)=OLDDEB(I)
   39 CONTINUE
*----
*  Shifted copy for parameter values.
*  Computing the address for new added value.
*----
      DO 45 I=OLDDEB(II),OLDDEB(II+1)-1
         IF(MUPLET(ISTART).LT.OLDARB(I))THEN
            II=I
            GO TO 46
         ENDIF
   45 CONTINUE
      II=OLDDEB(II+1)
   46 CONTINUE
*
      INCR=0
      DO 70 IPAR=ISTART,NPAR
*
         DO 55 I=I0,II-1
            NEWARB(I+INCR)=OLDARB(I)
   55    CONTINUE
*
         NEWARB(II+INCR)=MUPLET(IPAR)
         INCR=INCR+1
*
         DO 65 I=II,I1
            NEWARB(I+INCR)=OLDARB(I)
   65    CONTINUE
*
         IF(IPAR.NE.NPAR)THEN
            II=OLDDEB(II)
            I0=OLDDEB(I0)
            I1=OLDDEB(I0)-1
         ENDIF
*
   70 CONTINUE
*
      RETURN
      END
