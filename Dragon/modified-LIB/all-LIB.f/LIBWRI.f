*DECK LIBWRI
      SUBROUTINE LIBWRI(NTMPR,NDILR,TMPISO,DILISO,TMPT,DILT,REST,RIT,
     >                  XSOUT,XSCOR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Resonance integral temperature and dilution interpolation.
*
*Copyright:
* Copyright (C) 1997 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): 
* G. Marleau
*
*Parameters: input
* NTMPR   number of temperature tables.             
* NDILR   number of dilution tables.                
* TMPISO  temperature of isotope.                   
* DILISO  dilution of isotope.                      
* TMPT    sqrt(temperature) in table.               
* DILT    sqrt(dilution) in table.                  
* REST    resonance rates input.                    
*
*Parameters: output
* XSOUT   resonance integrals.               
* XSCOR   resonance integrals correction.              
*
*Parameters: scratch
* RIT     dummy vector.                             
*
*-----------------------------------------------------------------------
*
      IMPLICIT   NONE
      REAL       EPSRI,SQDILI
      PARAMETER (EPSRI=0.0005,SQDILI=1.0E+5)
      INTEGER    NTMPR,NDILR
      REAL       TMPISO,DILISO,TMPT(NTMPR),DILT(NDILR),
     >           REST(NDILR,NTMPR),RIT(NDILR),XSOUT,XSCOR
      INTEGER    IDIL,NDILE,JDEPT,JFINT,IR,JD
      REAL       SQTD,SQTT,ALPHA,AIKINT,ASLOPE
*----
*  SIMPLE LINEAR INTERFOLATION IN SQRT(TMP)
*----
      SQTD=SQRT(DILISO)
      XSOUT=0.0
      DO 110 IDIL=NDILR,1,-1
        IF(DILT(IDIL).LT.SQDILI) THEN
          NDILE=IDIL
          GO TO 115
        ENDIF
 110  CONTINUE
      RETURN
 115  CONTINUE
      SQTT=SQRT(TMPISO)
      IF(NTMPR.EQ.1) THEN
        JDEPT=1
        JFINT=1
        ALPHA=0.0
      ELSE
        JDEPT=1
        DO 100 IR=1,NTMPR-1
          IF(SQTT.GE.TMPT(IR)) JDEPT=IR
 100    CONTINUE
        JFINT=JDEPT+1
        ALPHA=(SQTT-TMPT(JDEPT))/(TMPT(JFINT)-TMPT(JDEPT))
      ENDIF
      DO 120 JD=1,NDILR
        RIT(JD)=(1.-ALPHA)*REST(JD,JDEPT)+
     >           ALPHA*REST(JD,JFINT)
 120  CONTINUE
      IF(SQTD .GT. DILT(NDILE)) THEN
*----
*  INTERPOLATE LINEARLY BETWEEN LAST DILUTION IN TABLE
*  AND INFINITE DILUTION
*----
        ASLOPE=(DILT(NDILE)/SQTD)**2
        XSOUT=ASLOPE*RIT(NDILE)+(1.0-ASLOPE)*RIT(NDILR)
      ELSE
*----
*  AIKINT INTERPOLATION FOR DILUTION
*----
        XSOUT=AIKINT(SQTD,DILT,RIT,NDILE,EPSRI)
      ENDIF
      XSCOR=XSCOR+XSOUT
      RETURN
      END
