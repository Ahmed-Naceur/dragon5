*DECK MSTGET
      SUBROUTINE MSTGET(IPSTR,IPRINT,IBEG,IEND,IINC,NAME)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Retrieve data from an existing block and put them into input variables.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* IPSTR   structure address.
* IPRINT  level of print index.
* IBEG    index of the first element.
* IEND    index of the last element.
* IINC    index increment between two consecutive elements.
* NAME    block name.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) :: IPSTR
      INTEGER :: IPRINT,IBEG,IEND,IINC
      CHARACTER(LEN=12) :: NAME  
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER :: IOUT=6
      INTEGER :: INDIC,NITMA,NARA,ITYP,SIZE,II,JJ
      REAL :: FLOTT
      DOUBLE PRECISION :: DFLOTT
      CHARACTER(LEN=12) :: TEXT12
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IARA
      REAL, ALLOCATABLE, DIMENSION(:) :: ARA
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DARA
*----
*  RETRIEVING BLOCK
*----
      CALL LCMLEN(IPSTR,NAME,SIZE,ITYP)
      IF (SIZE.LE.0) THEN
         CALL LCMLIB(IPSTR)
         CALL XABORT('MSTGET: INVALID BLOCK '//NAME//'.')
      ENDIF
      NARA=0
      IF (ITYP.EQ.1) THEN
         NARA=SIZE
         ALLOCATE(IARA(NARA))
         CALL LCMGET(IPSTR,NAME,IARA)
      ELSEIF (ITYP.EQ.2) THEN
         NARA=SIZE
         ALLOCATE(ARA(NARA))
         CALL LCMGET(IPSTR,NAME,ARA)
      ELSEIF (ITYP.EQ.3) THEN
         NARA=SIZE/3
         ALLOCATE(IARA(3*NARA))
         CALL LCMGET(IPSTR,NAME,IARA)
      ELSEIF (ITYP.EQ.4) THEN
         NARA=SIZE
         ALLOCATE(DARA(NARA))
         CALL LCMGET(IPSTR,NAME,DARA)
      ELSE 
         CALL XABORT('MSTGET: UNSUPPORTED TYPE')
      ENDIF 
      IF (IEND.GT.NARA) CALL XABORT('MSTGET: INCOMPATIBLE SIZE')
      IF (IPRINT.GT.2) 
     1 WRITE(IOUT,*) 'MSTGET: RETRIEVING DATA FROM '//NAME//' BLOCK'
*     PUT USER REQUESTED DATA IN INPUT VARIABLES
      DO II=IBEG,IEND,IINC
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF (-INDIC.NE.ITYP)
     1      CALL XABORT('MSTGET: INVALID VARIABLE TYPE.')
         IF (ITYP.EQ.1) THEN
            NITMA=IARA(II)
            DEALLOCATE(IARA)
         ELSEIF (ITYP.EQ.2) THEN
            FLOTT=ARA(II)
            DEALLOCATE(ARA)
         ELSEIF (ITYP.EQ.3) THEN
            WRITE(TEXT12,'(3A4)') (IARA(3*(II-1)+JJ),JJ=0,2)
            NITMA=12
            DEALLOCATE(IARA)
         ELSEIF (ITYP.EQ.4) THEN
            DFLOTT=DARA(II)
            DEALLOCATE(DARA)
         ENDIF
         CALL REDPUT(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
      ENDDO
      RETURN
      END
