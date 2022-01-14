*DECK LIBND0
      SUBROUTINE LIBND0 (NAMFIL,NGRO,IPENER)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover energy group information from a NDAS library.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
*
*Author(s): A. Hebert
*
*Parameters: input
* NAMFIL  name of the NDAS file.
*
*Parameters: output
* NGRO    number of energy groups.
* IPENER  pointer of the energy mesh limit array.
*
*Reference:
* Copyright (C) from NDAS Atomic Energy of Canada Limited utility (2006)
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE FSDF
      IMPLICIT NONE
*----
*  Subroutine arguments
*----
      INTEGER NGRO
      CHARACTER NAMFIL*(*)
      TYPE(C_PTR) IPENER
*----
*  Local variables
*----
      INTEGER IERR,HEADER(16)
      REAL, POINTER, DIMENSION(:) :: ENERG
*
      CALL XSDOPN(NAMFIL,IERR)
      IF(IERR.NE.0) CALL XABORT('LIBND0: XSDOPN could not open Library'
     >  //' files')
      CALL XSDBLD(6001,HEADER,IERR)
      IF(IERR.NE.0) CALL XABORT('LIBND0: XSDBLD could not read library'
     > //' parameters')
      NGRO=HEADER(2)
      IPENER=LCMARA(NGRO+1)
      CALL C_F_POINTER(IPENER,ENERG,(/ NGRO+1 /))
      CALL XSDBLD(5019,ENERG,IERR)
      IF(IERR.NE.0) CALL XABORT('LIBND0: XSDBLD could not read energy '
     >  //'group limits')
      CALL XSDCL()
      RETURN
      END
