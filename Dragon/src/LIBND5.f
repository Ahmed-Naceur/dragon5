*DECK LIBND5
      SUBROUTINE LIBND5(CFILNA,NEL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Initialize dimensions for depletion data with NDAS library files.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
*
*Author(s): A. Hebert
*
*Parameters: input
* CFILNA  NDAS file name.
*
*Parameters: output
* NEL     number of isotopes on library.
*
*Reference:
* Copyright (C) from NDAS Atomic Energy of Canada Limited utility (2006)
*
*-----------------------------------------------------------------------
*
      USE FSDF
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER CFILNA*(*)
      INTEGER NEL
*----
*  Local variables
*----
      CHARACTER TEXT12*12
      INTEGER IND,IERR,HEADER(16),IHEAD(200)
*----
*  PROBE AND OPEN THE NDAS FILE.
*----
      CALL XSDOPN(CFILNA,IERR)
      IF(IERR.NE.0) THEN
         TEXT12=CFILNA
         CALL XABORT('LIBND5: NDAS library '//TEXT12//' cannot be'//
     >   ' opened')
      ENDIF
      CALL XSDBLD(6001,HEADER,IERR)
      IF(IERR.NE.0) CALL XABORT('LIBND5: XSDBLD could not read library'
     > //' parameters')
      NEL=0
      DO IND=1,HEADER(1)
*       Load nuclide header
        CALL XSDISO(7000,6001,IND,IHEAD,IERR)
        IF(IHEAD(1).NE.0) NEL=NEL+1
      ENDDO
      CALL XSDCL()
      RETURN
      END
