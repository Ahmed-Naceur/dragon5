*DECK MSTPUT
      SUBROUTINE MSTPUT(IPSTR,IPRINT,IBEG,IEND,IINC,NAME,NBLOCK,BLNAM,
     1           BLTYP,BLLEN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Create or update a block in a structure from user input data.
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
* NBLOCK  number of existing block in the directory.
* BLNAM   names of these blocks.
* BLTYP   types of these blocks.
* BLLEN   lengths of these blocks.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) :: IPSTR
      INTEGER :: IPRINT,IBEG,IEND,IINC,NBLOCK,BLTYP(NBLOCK+1),
     1           BLLEN(NBLOCK+1)
      CHARACTER(LEN=12) :: BLNAM(NBLOCK+1),NAME 
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER :: IOUT=6
      INTEGER, PARAMETER :: NSTATE=40
      INTEGER :: INDIC,NITMA,ISTATE(NSTATE),II,JJ,ITYPO,NELEO,ITYP,
     1           NARA,SIZE
      REAL :: FLOTT
      DOUBLE PRECISION :: DFLOTT
      CHARACTER(LEN=12) :: TEXT12,WHITE12
      LOGICAL :: EXIST
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IARA
      REAL, ALLOCATABLE, DIMENSION(:) :: ARA
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DARA
*
*     BLOCK NAME
      EXIST=.FALSE.
      NELEO=0
*     SPECIAL CASE OF STATE-VECTOR MODIFICATION
      IF (NAME.EQ.'STATE-VECTOR') THEN
         IF (IEND.GT.40)
     1      CALL XABORT('MSTPUT: STATE-VECTOR SIZE IS LIMITED TO 40.')
         IF (IEND.EQ.40)
     2      CALL XABORT('MSTPUT: 40th STATE-VECTOR ELEMENT SHOULD'//
     3      ' NOT BE MODIFIED.')
         ITYPO=1
         NELEO=NSTATE
         EXIST=.TRUE.
      ENDIF
      IF (NBLOCK.NE.0) THEN
*     IS THIS BLOCK ALREADY PART OF THE STRUCTURE ?
         DO II=1,NBLOCK
            IF(BLNAM(II).EQ.NAME) THEN
               ITYPO=BLTYP(II)
               IF (ITYPO.EQ.0)
     1       CALL XABORT('MSTPUT: '//NAME//' IS AN EXISTING DIRECTORY.')
               NELEO=BLLEN(II)
               EXIST=.TRUE.
               GOTO 20
            ENDIF
         ENDDO
 20      CONTINUE
      ENDIF
      IF (IPRINT.GT.2) THEN
      IF (EXIST) THEN
*     YES: IT WILL BE UPDATED
         WRITE(IOUT,*) 'MSTPUT: BLOCK '//NAME//' IN UPDATE MODE'
      ELSE
*     NO: IT WILL BE CREATED
         WRITE(IOUT,*) 'MSTPUT: BLOCK '//NAME//' IN CREATION MODE'
      ENDIF
      ENDIF
*     FIRST ELEMENT
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      ITYP=INDIC
      IF ((EXIST).AND.(ITYP.NE.ITYPO))
     1   CALL XABORT('MSTPUT: HETEROGENEOUS BLOCK NOT SUPPORTED')
      IF (INDIC.GT.4)
     1   CALL XABORT('MSTPUT: UNSUPPORTED TYPE(1)')
      NARA=MAX(IEND,NELEO)
*     ALLOCATE MEMORY
      IF (ITYP.EQ.1) THEN
         SIZE=NARA
         ALLOCATE(IARA(NARA))
      ELSEIF (ITYP.EQ.2) THEN
         SIZE=NARA
         ALLOCATE(ARA(NARA))
      ELSEIF (ITYP.EQ.3) THEN
         SIZE=3*NARA
         ALLOCATE(IARA(3*NARA))
      ELSEIF (ITYP.EQ.4) THEN
         SIZE=NARA
         ALLOCATE(DARA(NARA))
      ENDIF
*     INITIALIZE BLOCK
      IF (EXIST) THEN
         IF (ITYP.EQ.1) THEN
            CALL LCMGET(IPSTR,NAME,IARA)
         ELSEIF (ITYP.EQ.2) THEN
            CALL LCMGET(IPSTR,NAME,ARA)
         ELSEIF (ITYP.EQ.3) THEN
            CALL LCMGET(IPSTR,NAME,IARA)
         ELSEIF (ITYP.EQ.4) THEN
            CALL LCMGET(IPSTR,NAME,DARA)
         ENDIF
      ELSE
         IF (ITYP.EQ.1) THEN
            IARA(:IEND)=0
         ELSEIF (ITYP.EQ.2) THEN
            ARA(:IEND)=0.0
         ELSEIF (ITYP.EQ.3) THEN
            WHITE12=' '
            DO II=1,IEND
               READ(WHITE12,'(3A4)') (IARA(3*(II-1)+JJ),JJ=0,2)
            ENDDO
         ELSEIF (ITYP.EQ.4) THEN
            DARA(:IEND)=0.D0
         ENDIF
      ENDIF
*     RETRIEVE USER'S INPUT VALUES
      DO II=IBEG,IEND,IINC
         IF (INDIC.NE.ITYP)
     1      CALL XABORT('MSTPUT: HETEROGENEOUS BLOCK NOT SUPPORTED')
         IF (INDIC.EQ.1) THEN
            IARA(II)=NITMA
         ELSEIF (INDIC.EQ.2) THEN
            ARA(II)=FLOTT
         ELSEIF (INDIC.EQ.3) THEN
            READ(TEXT12,'(3A4)') (IARA(3*(II-1)+JJ),JJ=0,2)
         ELSEIF (INDIC.EQ.4) THEN
            DARA(II)=DFLOTT
         ELSE
            CALL XABORT('MSTPUT: UNSUPPORTED TYPE(2)')
         ENDIF
         IF (II.LT.IEND) CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      ENDDO
      IF (INDIC.EQ.1) THEN
         CALL LCMPUT(IPSTR,NAME,SIZE,ITYP,IARA)
         DEALLOCATE(IARA)
      ELSEIF (INDIC.EQ.2) THEN
         CALL LCMPUT(IPSTR,NAME,SIZE,ITYP,ARA)
         DEALLOCATE(ARA)
      ELSEIF (INDIC.EQ.3) THEN
         CALL LCMPUT(IPSTR,NAME,SIZE,ITYP,IARA)
         DEALLOCATE(IARA)
      ELSEIF (INDIC.EQ.4) THEN
         CALL LCMPUT(IPSTR,NAME,SIZE,ITYP,DARA)
         DEALLOCATE(DARA)
      ENDIF
*----
*  UPDATE NB. BLOCKS, REC-NAMES, REC-TYPES, REC-LENGTHS IN STATE-VECTOR
*  IF REQUIRED
*----
      IF (.NOT.EXIST) THEN
         CALL LCMGET(IPSTR,'STATE-VECTOR',ISTATE)
         ISTATE(40)=ISTATE(40)+1
         BLNAM(NBLOCK+1)=NAME
         BLTYP(NBLOCK+1)=ITYP
         BLLEN(NBLOCK+1)=NARA
         CALL LCMPUT(IPSTR,'STATE-VECTOR',NSTATE,1,ISTATE)
         CALL LCMPTC(IPSTR,'REC-NAMES',12,NBLOCK+1,BLNAM)
         CALL LCMPUT(IPSTR,'REC-TYPES',(NBLOCK+1),1,BLTYP)
         CALL LCMPUT(IPSTR,'REC-LENGTHS',(NBLOCK+1),1,BLLEN)
      ENDIF
      RETURN
      END
