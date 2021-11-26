*DECK LIBLAG
      SUBROUTINE LIBLAG(NEF,XE,GE,XI,GS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Lagrange interpolation in a table of points using the APOLLO2 recipe.
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
* NEF     number of points.
* XE      x-values.
* GE      f(x)-values.
* XI      interpolating x-value.
*
*Parameters: output
* GS      interpolated value f(XI).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NEF
      REAL XE(NEF),GE(NEF),XI,GS
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NINT=3)
      REAL WEIGHT(NINT)
*----
*  CONSTANT FUNCTION DEFINED BY A SINGLE POINT
*----
      IF(NEF.EQ.1) THEN
         GS=GE(1)
         RETURN
      ENDIF
*
      IORD=MIN(NINT,NEF)
*----
*  LOCATE FIRST POINT TO THE RIGHT OF XI (PT IXP)
*----
      DO 30 IXP=1,NEF
       IF(ABS(XE(IXP)-XI).LE.1.0E-5*ABS(XI)) THEN
          GS=GE(IXP)
          RETURN
       ELSE IF(XE(IXP).GT.XI) THEN
          IMIN=IXP-(IORD+1)/2
          IMAX=IMIN+IORD-1
          GO TO 40
       ENDIF
   30 CONTINUE
*----
*  XI IS TO THE RIGHT OF EVERY POINT XE
*----
      IMAX=NEF
      IMIN=NEF-IORD+1
      GO TO 70
*
   40 IF(IMIN.LT.1) THEN
         IMIN=1
         IMAX=IORD
      ELSE IF(IMAX.GT.NEF) THEN
         IMIN=NEF-IORD+1
         IMAX=NEF
      ENDIF
*
   70 I0=IMIN-1
      DO 90 I=IMIN,IMAX
      PP=1.0
      DO 80 J=IMIN,IMAX
       IF(I.NE.J) PP=PP*((XI-XE(J))/(XE(I)-XE(J)))
   80 CONTINUE
      WEIGHT(I+1-IMIN)=PP
   90 CONTINUE
*
      GS=0.0
      DO 110 I=1,IORD
      GS=GS+WEIGHT(I)*GE(I+I0)
  110 CONTINUE
      DO 120 I=1,IORD
      I1=I+I0
      IF(XE(I1).GT.XI) THEN
         IF(I1-1.GT.0) THEN
            YMIN=MIN(GE(I1-1),GE(I1))
            YMAX=MAX(GE(I1-1),GE(I1))
            IF((GS.GT.YMAX).OR.(GS.LT.YMIN)) THEN
               GS=GE(I1-1)+(GE(I1)-GE(I1-1))*
     1         (XI-XE(I1-1))/(XE(I1)-XE(I1-1))
            ENDIF
         ELSE
            GS=GE(1)+(GE(2)-GE(1))*(XI-XE(1))/(XE(2)-XE(1))
            IF(GS.LE.0.) GS=GE(1)
         ENDIF
         RETURN
      ENDIF
  120 CONTINUE
      RETURN
      END
