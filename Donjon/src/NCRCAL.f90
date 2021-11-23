RECURSIVE INTEGER FUNCTION NCRCAL(II,NVP,NPTOT,DEBARB,ARBVAL,MUPLET) RESULT(ICAL)
!
!-----------------------------------------------------------------------
!
!Purpose:
! find the position of an elementary calculation in the multicompo or
! in the Saphyb.
!
!Copyright:
! Copyright (C) 2012 Ecole Polytechnique de Montreal
!
!Author(s): 
! A. Hebert
!
!Parameters: input
! II      position in DEBARB. Must be set to 1 in the first call.
! NVP     number of nodes in the parameter tree.
! NPTOT   number of parameters.
! DEBARB  tree information
! ARBVAL  tree information
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
  INTEGER IKEEP, I, JICAL, NBOK
  INTEGER II,NVP,NPTOT,DEBARB(NVP+1),ARBVAL(NVP),MUPLET(NPTOT)
  !
  IF(NPTOT==0) THEN
    ICAL=DEBARB(II+1)
    RETURN
  ENDIF
  NBOK=0
  IKEEP=0
  DO I=DEBARB(II),DEBARB(II+1)-1
    IF((MUPLET(1)==0).OR.(MUPLET(1)==ARBVAL(I))) THEN
      JICAL=NCRCAL(I,NVP,NPTOT-1,DEBARB,ARBVAL,MUPLET(2))
      IF(JICAL > 0) THEN
        IKEEP=JICAL
        NBOK=NBOK+1
      ELSE IF(JICAL==-1) THEN
        NBOK=2
      ENDIF
    ENDIF
  ENDDO
  IF(NBOK > 1) THEN
    ! Many elementary calculation exist for this tuple.
    ICAL=-1
  ELSE IF(NBOK==0) THEN
    ! No elementary calculation exists for this tuple.
    ICAL=0
  ELSE
    ICAL=IKEEP
  ENDIF
END FUNCTION NCRCAL
