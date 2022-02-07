*DECK TRIKAX
      SUBROUTINE TRIKAX (IDIM,NCODE,XXX,YYY,ZZZ,LX,LY,LZ,IAXIS,CENTER)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculates the center of the external cylinder outside elements.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IDIM    number of dimensions.
* XXX     Cartesian coordinates of the domain along the X-axis.
* YYY     Cartesian coordinates of the domain along the Y-axis.
* ZZZ     Cartesian coordinates of the domain along the Z-axis.
* LX      number of parallelepipeds along the X-axis after mesh-
*         splitting.
* LY      number of parallelepipeds along the Y-axis.
* LZ      number of parallelepipeds along the Z-axis.
* NCODE   boundary condition relative to each side of the domain.
*
*Parameters: output
* CENTER  coordinates for center of cylinder.
* IAXIS   orientation of the cylinder axis: = 0 no cylinder at all;
*         = 1,2,3 axis of the cylinder.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IDIM,NCODE(6),LX,LY,LZ,IAXIS
      REAL XXX(LX+1),YYY(LY+1),ZZZ(LZ+1),CENTER(3)
*----
*  LOCAL VARIABLES
*----
      INTEGER IFC(3)
*
      IAXIS = 0
      DO 10 IC= 1,3
         CENTER(IC)= 0.0
         IFC(IC)= 0
   10 CONTINUE
      IF( IDIM.GE.2 )THEN
*----
* "X" AXIS STUDY
*----
         IF( NCODE(1).EQ.20.OR.NCODE(2).EQ.20 ) THEN
*        THERE IS AT LEAST ONE "X" CIRCULAR B.C.
            IFC(1)= 1
            IF( NCODE(1).EQ.20.AND.NCODE(2).EQ.20 )THEN
*           THERE IS TWO "X" CIRCULAR B.C.
               CENTER(1)= 0.5 * (XXX(LX+1) + XXX(1))
*              TAKE THE "X" CENTER AT THE MIDDLE OF ALL ELEMENTS
            ELSEIF( NCODE(1).EQ.5.OR.NCODE(2).EQ.5 )THEN
*           THERE IS ONE "X" SYMMETRIC B.C.
               IF( NCODE(1).EQ.5 )THEN
*              "X -" SYMMETRIC B.C.
                  CENTER(1)= 0.5 * (XXX(2) + XXX(1))
*                 TAKE THE "X" CENTER AT THE MIDDLE OF FIRST ELEMENT
               ELSE
*              "X +" SYMMETRIC B.C.
                  CENTER(1)= 0.5 * (XXX(LX+1) + XXX(LX))
*                 TAKE THE "X" CENTER AT THE MIDDLE OF LAST  ELEMENT
               ENDIF
            ELSE
*           ALL OTHER CASES
               IF( NCODE(1).EQ.20 )THEN
*              "X -" CIRCULAR B.C.
                  CENTER(1)= XXX(LX+1)
*                 TAKE THE "X" CENTER AT THE END    OF LAST  ELEMENT
               ELSE
*              "X +" SYMMETRIC B.C.
                  CENTER(1)= XXX(1)
*                 TAKE THE "X" CENTER AT THE BEGIN  OF FIRST ELEMENT
               ENDIF
            ENDIF
         ENDIF
*----
* "Y" AXIS STUDY
*----
         IF( NCODE(3).EQ.20.OR.NCODE(4).EQ.20 ) THEN
            IFC(2)= 1
*        THERE IS AT LEAST ONE "Y" CIRCULAR B.C.
            IF( NCODE(3).EQ.20.AND.NCODE(4).EQ.20 )THEN
*           THERE IS TWO "Y" CIRCULAR B.C.
               CENTER(2)= 0.5 * (YYY(LY+1) + YYY(1))
*              TAKE THE "Y" CENTER AT THE MIDDLE OF ALL ELEMENTS
            ELSEIF( NCODE(3).EQ.5.OR.NCODE(4).EQ.5 )THEN
*           THERE IS ONE "Y" SYMMETRIC B.C.
               IF( NCODE(3).EQ.5 )THEN
*              "Y -" SYMMETRIC B.C.
                  CENTER(2)= 0.5 * (YYY(2) + YYY(1))
*                 TAKE THE "Y" CENTER AT THE MIDDLE OF FIRST ELEMENT
               ELSE
*              "Y +" SYMMETRIC B.C.
                  CENTER(2)= 0.5 * (YYY(LY+1) + YYY(LY))
*                 TAKE THE "Y" CENTER AT THE MIDDLE OF LAST  ELEMENT
               ENDIF
            ELSE
*           ALL OTHER CASES
               IF( NCODE(3).EQ.20 )THEN
*              "Y -" CIRCULAR B.C.
                  CENTER(2)= YYY(LY+1)
*                 TAKE THE "Y" CENTER AT THE END    OF LAST  ELEMENT
               ELSE
*              "Y +" SYMMETRIC B.C.
                  CENTER(2)= YYY(1)
*                 TAKE THE "Y" CENTER AT THE BEGIN  OF FIRST ELEMENT
               ENDIF
            ENDIF
         ENDIF
         IF( IDIM.EQ.2 )THEN
            NONC  = IFC(1) + IFC(2)
            IF( NONC.GT.0 )THEN
               IAXIS = 3
            ENDIF
         ELSE
*----
* "Z" AXIS STUDY
*----
            IF( NCODE(5).EQ.20.OR.NCODE(6).EQ.20 ) THEN
*           THERE IS AT LEAST ONE "Y" CIRCULAR B.C.
               IFC(3)= 1
               IF( NCODE(5).EQ.20.AND.NCODE(6).EQ.20 )THEN
*              THERE IS TWO "Z" CIRCULAR B.C.
                  CENTER(3)= 0.5 * (ZZZ(LZ+1) + ZZZ(1))
*                 TAKE THE "Z" CENTER AT THE MIDDLE OF ALL ELEMENTS
               ELSEIF( NCODE(5).EQ.5.OR.NCODE(6).EQ.5 )THEN
*              THERE IS ONE "Z" SYMMETRIC B.C.
                  IF( NCODE(5).EQ.5 )THEN
*                 "Z -" SYMMETRIC B.C.
                     CENTER(3)= 0.5 * (ZZZ(2) + ZZZ(1))
*                    TAKE THE "Z" CENTER AT THE MIDDLE OF FIRST ELEMENT
                  ELSE
*                 "Z +" SYMMETRIC B.C.
                     CENTER(3)= 0.5 * (ZZZ(LZ+1) + ZZZ(LZ))
*                    TAKE THE "Z" CENTER AT THE MIDDLE OF LAST  ELEMENT
                  ENDIF
               ELSE
*              ALL OTHER CASES
                  IF( NCODE(5).EQ.20 )THEN
*                 "Z -" CIRCULAR B.C.
                     CENTER(3)= ZZZ(LZ+1)
*                    TAKE THE "Z" CENTER AT THE END    OF LAST  ELEMENT
                  ELSE
*                 "Z +" SYMMETRIC B.C.
                     CENTER(3)= ZZZ(1)
*                    TAKE THE "Z" CENTER AT THE BEGIN  OF FIRST ELEMENT
                  ENDIF
               ENDIF
            ENDIF
*
*           DETERMINE PRINCIPAL AXIS
            NONC= IFC(1) + IFC(2) + IFC(3)
            IF( NONC.GT.0 )THEN
               IF( NONC.EQ.2 )THEN
                  IF( IFC(1).EQ.0 ) IAXIS = 1
                  IF( IFC(2).EQ.0 ) IAXIS = 2
                  IF( IFC(3).EQ.0 ) IAXIS = 3
               ELSE
                  WRITE(6,1000)
                  CALL XABORT('TRIKAX: ALGORITHM FAILURE.')
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      RETURN
 1000 FORMAT(/1X,'*** NOT POSSIBLE TO DETERMINE THE PRINCIPAL AXIS'
     1       /1X,'***'
     2       /1X,'***  N O   C Y L I N D R I C   A L B E D O S'
     3       /1X,'***')
      END
