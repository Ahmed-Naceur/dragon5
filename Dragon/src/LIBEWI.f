*DECK LIBEWI
      SUBROUTINE LIBEWI(CFILNA,NEL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Initialize dimensions for depletion data on WIMS-AECL
* format library.
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
* CFILNA  file name.
*
*Parameters: output
* NEL     number of isotopes on library.
*
*Comments:
*   WIMS-AECL library parameters
*   MAXISO : max. nb. of iso = 246                
*   MLDEP  : maximum number of reaction per       
*            isotope = MAXISO +4
*   LPZ    : length of parameter array = 9   
*   LMASTB : length of mst tab = MAXISO+9         
*   LMASIN : length of mst idx = LMASTB-4         
*   LGENTB : length of gen tab = 6                
*   LGENIN : length of gen idx = LGENTB
*   MASTER : master index array                   
*   GENINX : general index array
*   NPZ    : list of main parameters              
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  WIMS-AECL LIBRARY PARAMETERS
*----
      INTEGER      LRIND,IACTO,IACTC,IUTYPE,MAXISO,MLDEP,LPZ,
     1             MAXTEM,LMASTB,LMASIN,LGENTB,LGENIN
      PARAMETER   (LRIND=256,IACTO=2,IACTC=1,IUTYPE=4,MAXISO=246,
     1             MLDEP=MAXISO+4,LPZ=9,MAXTEM=20,LMASTB=MAXISO+9,
     2             LMASIN=LMASTB-4,LGENTB=6,LGENIN=LGENTB)
      INTEGER      MASTER(LMASTB),GENINX(LGENTB),NPZ(LPZ)
*----
* EXTERNAL FUNCTIONS
*----
      INTEGER      KDROPN
*----
* LOCAL VARIABLES
*----
      INTEGER      NEL,IUNIT
      CHARACTER    CFILNA*8
*----
*  OPEN WIMS-AECL LIBRARY
*  READ INDEX AND GENERAL DIMENSIONING NPZ
*----
      IUNIT=KDROPN(CFILNA,IACTO,IUTYPE,LRIND)
      IF(IUNIT.LE.0) CALL XABORT('LIBEWI: WIMS-AECL LIBRARY '//
     >    CFILNA//' CANNOT BE OPENED FOR DEPLETION')
      CALL OPNIND(IUNIT,MASTER,LMASTB)
      CALL REDIND(IUNIT,MASTER,LMASIN,GENINX,LGENTB,1)
      CALL REDIND(IUNIT,GENINX,LGENIN,NPZ,LPZ,1)
      NEL=NPZ(1)
      IF(NEL.GT.MAXISO) CALL XABORT('LIBEWI: TOO MANY ISOTOPES '//
     >    'ON WIMS-AECL LIBRARY'//CFILNA)
*----
*  CLOSE WIMS-AECL LIBRARY AND
*  RETURN
*----
      CALL CLSIND(IUNIT)
*----
*  RETURN
*----
      RETURN
      END
