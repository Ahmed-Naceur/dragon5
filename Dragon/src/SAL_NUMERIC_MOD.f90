!
!---------------------------------------------------------------------
!
!Purpose:
! Support module for numerical functions.
!
!Copyright:
! Copyright (C) 2001 Ecole Polytechnique de Montreal.
!
!Author(s):
! X. Warin
!
!---------------------------------------------------------------------
!
MODULE SAL_NUMERIC_MOD
  
  USE PRECISION_AND_KINDS, ONLY : PDB

CONTAINS
  !
  FUNCTION SALACO(COSANG,Y)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! computes angle in radians for given cosinus and y component
    !
    !Parameters: input
    ! COSANG   cosinus of angle
    ! Y        component (to give sign)
    !
    !Parameters: output
    ! SALACO   angle in radiants
    !
    !---------------------------------------------------------------------
    !
    USE PRECISION_AND_KINDS, ONLY : PDB,PI,TWOPI
    !**
    REAL(PDB) :: SALACO
    REAL(PDB),INTENT(IN) :: COSANG,Y
    !*****
    IF(ABS(COSANG).LT.1.0_PDB) THEN
       SALACO=ACOS(COSANG)
    ELSEIF(COSANG.GE.1.0_PDB) THEN
       SALACO=0.0_PDB
    ELSE
       SALACO=PI
    ENDIF
    IF(Y.LT.0.0_PDB) SALACO=TWOPI-SALACO
    !
  END FUNCTION SALACO
  !
  SUBROUTINE SAL141(TYPE,RPAR,X,Y,IEND)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! computes coordinates of end of an element
    !
    !Parameters: input
    ! TYPE    type of element  1 (segment)  3 (arc of circle)
    ! RPAR    floating-point descriptors of the element
    ! IEND    = 1 (end is origin of element)
    !           2 (end is end of the element)
    !
    !Parameters: output
    ! X       abscissa coordinates of end
    ! Y       ordinate coordinates of end
    !
    !---------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER,   INTENT(IN)               :: TYPE,IEND
    REAL(PDB), INTENT(OUT)              :: X,Y
    REAL(PDB), INTENT(IN), DIMENSION(:) :: RPAR
    !     DIMENSION        RPAR(*)
    !****
    REAL(PDB) :: THETA,R
    !****
    X=RPAR(1)
    Y=RPAR(2)
    IF(TYPE.EQ.1)THEN
       !        segment
       IF(IEND.EQ.2)THEN
          X=X+RPAR(3)
          Y=Y+RPAR(4)
       ENDIF
    ELSEIF(TYPE.LE.3)THEN
       !        arc of circle
       IF(IEND.EQ.1)THEN
          THETA=RPAR(4)
       ELSE
          THETA=RPAR(5)
       ENDIF
       R=RPAR(3)
       X=X+R*COS(THETA)
       Y=Y+R*SIN(THETA)
    ELSE
       CALL XABORT('SAL141: not implemented')
    ENDIF
    !
  END SUBROUTINE SAL141

END MODULE SAL_NUMERIC_MOD
