*DECK LIBENI
      SUBROUTINE LIBENI(CFILNA,NEL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Initialize dimensions for depletion data for WIMS-D4 format library.
*
*Copyright:
* Copyright (C) 1997 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau.
*
*Parameters: input
* CFILNA  WIMS file name.
*
*Parameters: output
* NEL     number of isotopes on library.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
* PARAMETERS
*----
      INTEGER      LRIND,IACTO,IACTC,IUTYPE,LPZ
      PARAMETER   (LRIND=0,IACTO=2,IACTC=1,IUTYPE=2,LPZ=8)
*----
* EXTERNAL FUNCTIONS
*----
      INTEGER      KDRCLS,KDROPN
*----
* LOCAL VARIABLES
*-----
      INTEGER      NEL,IUNIT,II,IERR
      CHARACTER    CFILNA*8
*----
*  WIMS-D4 LIBRARY PARAMETERS
*----
      INTEGER      NPZ(LPZ)
*----
*  OPEN WIMS-D4 LIBRARY
*  READ GENERAL DIMENSIONING
*----
      IUNIT=KDROPN(CFILNA,IACTO,IUTYPE,LRIND)
      IF(IUNIT.LE.0) CALL XABORT('LIBENI: WIMS-D4 LIBRARY '//
     >    CFILNA//' CANNOT BE OPENED FOR DEPLETION')
      READ(IUNIT) (NPZ(II),II=1,LPZ)
      NEL=NPZ(1)
      IERR=KDRCLS(IUNIT,IACTC)
      IF(IERR.LT.0)
     >  CALL XABORT('LIBENI: WIMS-D4 LIBRARY '//
     >    CFILNA//' CANNOT BE CLOSED')
*----
*  RETURN
*----
      RETURN
      END
