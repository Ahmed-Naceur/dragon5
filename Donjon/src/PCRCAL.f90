INTEGER FUNCTION PCRCAL(NPAR,NCAL,MUPLET,MUBASE) RESULT(ICAL)
!
!-----------------------------------------------------------------------
!
!Purpose:
! find the position of an elementary calculation in a PMAXS file.
!
!Copyright:
! Copyright (C) 2019 Ecole Polytechnique de Montreal
!
!Author(s): A. Hebert
!
!Parameters: input
! NPAR    number of parameters.
! NCAL    number of elementary calculations in the PMAXS file.
! MUPLET  tuple used to identify an elementary calculation.
!
!Parameters: output
! ICAL    position of the elementary calculation (=0 if does not exist;
!         =-1 if more than one exists).
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  !----
  !  SUBROUTINE ARGUMENTS
  !----
  INTEGER NPAR,NCAL,MUPLET(NPAR),MUBASE(NPAR,NCAL)
  !----
  !  LOCAL VARIABLES
  !----
  INTEGER I,J,NFIND
  !
  ICAL=0
  NFIND=0
  DO I=1,NCAL
    DO J=1,NPAR
      IF(MUPLET(J).NE.MUBASE(J,I)) GO TO 10
    ENDDO
    ICAL=I
    NFIND=NFIND+1
    10 CONTINUE
  ENDDO
  IF(NFIND.GT.1) ICAL=-1
END FUNCTION PCRCAL
