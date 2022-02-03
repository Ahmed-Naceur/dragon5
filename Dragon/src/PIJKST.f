*DECK PIJKST
      SUBROUTINE PIJKST(IMPX,NREGIO,PIJSYM,PIJKS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Evaluate PIJK*=PIJ**(-1)*PIJK.
*
*Copyright:
* Copyright (C) 1994 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau, I. Petrovic
*
*Parameters: input
* IMPX    print/check flag. 
* NREGIO  number of regions considered.
* PIJSYM  group condensed reduce/symmetric PIJ.
*
*Parameters: output
* PIJKS   group condensed PIJK*.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
* PARAMETERS
*----
      INTEGER     IUNOUT
      PARAMETER  (IUNOUT=6)
*----
* INTERNAL FUNCTIONS
*----
      INTEGER     INDPOS
*----
* LOCAL VARIABLES
*----
      INTEGER     IMPX,NREGIO,IDIR,I,J,INDPIJ,IERROR
      REAL        PIJSYM(*),PIJKS(NREGIO,NREGIO,3)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: PIJSCT
*
*----- INTRINSIC FUNCTION FOR POSITION IN CONDENSE PIJ MATRIX
*
      INDPOS(I,J)=MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PIJSCT(NREGIO,2*NREGIO))
*----
*  FILL SYSTEM MATRIX WITH PIJ
*----
      DO 100 IDIR=1,3
        DO 110 I=1,NREGIO
          DO 120 J=1,NREGIO
            INDPIJ=INDPOS(I,J)
            PIJSCT(I,J)=DBLE(PIJSYM(INDPIJ))
            PIJSCT(J,NREGIO+I)=DBLE(PIJKS(I,J,IDIR))
 120      CONTINUE
 110    CONTINUE
        CALL ALSBD(NREGIO,NREGIO,PIJSCT,IERROR,NREGIO)
        IF(IERROR.NE.0) CALL XABORT('PIJKST: SINGULAR MATRIX.')
        DO 130 I=1,NREGIO
          DO 140 J=1,NREGIO
            PIJKS(I,J,IDIR)=REAL(PIJSCT(I,NREGIO+J))
 140      CONTINUE
 130    CONTINUE
        IF (IMPX.GE.8) THEN
          WRITE(IUNOUT,6000) (J,J=1,NREGIO)
          DO 150 I=1,NREGIO
            WRITE(IUNOUT,6001) I,(PIJKS(I,J,IDIR),J=1,NREGIO)
 150      CONTINUE
          WRITE(IUNOUT,'(//)')
        ENDIF
 100  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(PIJSCT)
      RETURN
*
 6000 FORMAT (//'COLLISION PROBAB. MATRIX PIJK*=((PIJ)**(-1))*PIJK:'//
     1 (11X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,
     2 I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,
     3 2HJ=,I4,:,5X,2HJ=,I4))
 6001 FORMAT (3H I=,I4,2H: ,1P,11E11.3/(9X,11E11.3))
      END

