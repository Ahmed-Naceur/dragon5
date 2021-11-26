*DECK XCGBCM
      SUBROUTINE XCGBCM(IPTRK,NSOUT,NCODE,MATRT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Built boundary condition matrix for reflection and transmission.
*
*Copyright:
* Copyright (C) 1998 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G.Marleau
*
*Parameters: input
* IPTRK   pointer to the tracking file.
* NSOUT   number of outer surface.
* NCODE   albedo type.
*
*Parameters: output
* MATRT   BC matrix for reflection/transmission.
*
*----------------------------------------------------------------------
*
      USE GANLIB
      PARAMETER (NMCOD=6)
      TYPE(C_PTR) IPTRK
      INTEGER NSOUT,NCODE(NMCOD),MATRT(NSOUT),ISOUT
*----
*  INITIALIZE MATRT TO REFLECTION
*----
      DO 100 ISOUT=1,NSOUT
        MATRT(ISOUT)=ISOUT
 100  CONTINUE
*----
*  FOR CARTESIAN CELL LOOK AT PERIODIC BOUNDARY CONDITIONS
*  AND SET TRANSMISSION MATRIX
*----
      IF(NSOUT.EQ.4) THEN
        IF((NCODE(1) .EQ. 4) .AND. (NCODE(2) .EQ.4)) THEN
          MATRT(1)=3
          MATRT(3)=1
        ENDIF
        IF((NCODE(3) .EQ. 4) .AND. (NCODE(4) .EQ.4)) THEN
          MATRT(2)=4
          MATRT(4)=2
        ENDIF
      ENDIF
      CALL LCMPUT(IPTRK,'BC-REFL+TRAN',NSOUT,1,MATRT)
      RETURN
      END
