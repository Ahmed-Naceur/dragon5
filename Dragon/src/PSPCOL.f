*DECK PSPCOL
      SUBROUTINE PSPCOL(ITCOL,NCOL,ICOL,RGB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Pick a color number from a N-color set.
*
*Copyright:
* Copyright (C) 1999 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* ITCOL   type of color set:
*         = 1 gray;
*         = 2 rgb;
*         = 3 cmyk;
*         = 4 hsb.
* NCOL    maximum number of color in set.
* ICOL    requested color number.
*
*Parameters: input
* RGB     color intensity:
*         for gray use only RGB(1);
*         for rgb use only RGB(1),RGB(2),RGB(3);
*         for cmyk use all;
*         for hsb use only RGB(1),RGB(2),RGB(3).
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='PSPCOL')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          ITCOL,NCOL,ICOL
      REAL             RGB(4)
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IDC,JCOL
      REAL             DELCOL,DELSAT,DELBLK
*----
*  LOCAL VARIABLES
*----
      IF(ITCOL .EQ. 4) THEN
        RGB(4)=0.0
        IF(ICOL .LE. 0 ) THEN
          RGB(1)=0.0
          RGB(2)=0.0
          RGB(3)=1.0
        ELSE
          DELCOL=0.6667/FLOAT(NCOL-1)
          DELSAT=0.5/FLOAT(NCOL-1)
          DELBLK=0.5/FLOAT(NCOL-1)
          JCOL=ICOL-1
          RGB(1)=0.6667-DELCOL*FLOAT(JCOL)
          RGB(2)=0.5+DELSAT*FLOAT(JCOL)
          RGB(3)=0.5+DELBLK*FLOAT(JCOL)
        ENDIF
      ELSE IF(ITCOL .EQ. 3) THEN
        RGB(4)=0.0
        IF(ICOL .LE. 0 ) THEN
          RGB(1)=0.0
          RGB(2)=0.0
          RGB(3)=0.0
        ELSE
          IF     (NCOL .LE.       8) THEN
            IDC=2
          ELSE IF(NCOL .LE.      64) THEN
            IDC=4
          ELSE IF(NCOL .LE.     512) THEN
            IDC=8
          ELSE IF(NCOL .LE.    4096) THEN
            IDC=16
          ELSE IF(NCOL .LE.   32768) THEN
            IDC=32
          ELSE IF(NCOL .LE.  262144) THEN
            IDC=64
          ELSE
            IDC=128
          ENDIF
          JCOL=ICOL-1
          DELCOL=1.0/FLOAT(IDC)
          RGB(1)=1.0-DELCOL*FLOAT(MOD(JCOL,IDC))
          JCOL=JCOL/IDC
          RGB(2)=1.0-DELCOL*FLOAT(MOD(JCOL,IDC))
          JCOL=JCOL/IDC
          RGB(3)=1.0-DELCOL*FLOAT(MOD(JCOL,IDC))
        ENDIF
      ELSE IF(ITCOL .EQ. 2) THEN
        RGB(4)=0.0
        IF(ICOL .LE. 0 ) THEN
          RGB(1)=1.0
          RGB(2)=1.0
          RGB(3)=1.0
        ELSE
          IF     (NCOL .LE.       8) THEN
            IDC=2
          ELSE IF(NCOL .LE.      64) THEN
            IDC=4
          ELSE IF(NCOL .LE.     512) THEN
            IDC=8
          ELSE IF(NCOL .LE.    4096) THEN
            IDC=16
          ELSE IF(NCOL .LE.   32768) THEN
            IDC=32
          ELSE IF(NCOL .LE.  262144) THEN
            IDC=64
          ELSE
            IDC=128
          ENDIF
          JCOL=ICOL-1
          DELCOL=1.0/FLOAT(IDC)
          RGB(1)=1.0-DELCOL*FLOAT(MOD(JCOL,IDC)+1)
          JCOL=JCOL/IDC
          RGB(2)=1.0-DELCOL*FLOAT(MOD(JCOL,IDC)+1)
          JCOL=JCOL/IDC
          RGB(3)=1.0-DELCOL*FLOAT(MOD(JCOL,IDC)+1)
        ENDIF
      ELSE
        IF(ICOL .LE. 0 ) THEN
          RGB(1)=0.0
          RGB(2)=0.0
          RGB(3)=0.0
        ELSE
          IF     (NCOL .LE.       8) THEN
            IDC=8
          ELSE IF(NCOL .LE.      64) THEN
            IDC=64
          ELSE IF(NCOL .LE.     512) THEN
            IDC=512
          ELSE IF(NCOL .LE.    4096) THEN
            IDC=4096
          ELSE IF(NCOL .LE.   32768) THEN
            IDC=32768
          ELSE
            IDC=262144
          ENDIF
          JCOL=ICOL-1
          DELCOL=1.0/FLOAT(IDC)
          RGB(1)=1.0-DELCOL*FLOAT(MOD(JCOL,IDC))
          RGB(2)=RGB(1)
          RGB(3)=RGB(1)
        ENDIF
      ENDIF
      RETURN
      END
