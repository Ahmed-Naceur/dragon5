*DECK MSTCPB
      SUBROUTINE MSTCPB(IPSTR,IPSTR2,IPRINT,IBEG,IEND,IINC,NAME,NAME2,
     1                  NBLOCK,BLNAM,BLTYP,BLLEN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Copy some elements from a structure's block to another.
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
* IPSTR   address of the structure from which the information is
*         retrieved.
* IPSTR2  destination structure address.
* IPRINT  level of print index.
* IBEG    index of the first element.
* IEND    index of the last element.
* IINC    index increment between two consecutive elements.
* NAME    name of the block from which the information is retrieved.
* NAME2   destination block name.
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
      TYPE(C_PTR) :: IPSTR,IPSTR2
      INTEGER :: IPRINT,IBEG,IEND,IINC,NBLOCK,BLTYP(NBLOCK+1),
     1           BLLEN(NBLOCK+1)
      CHARACTER(LEN=12) :: BLNAM(NBLOCK+1),NAME,NAME2  
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER :: IOUT=6,NSTATE=40
      INTEGER :: ISTATE(NSTATE),NARA,ITYP,SIZE,ITYPO,NELEO,NARA2,II,JJ,
     1           SIZE2
      CHARACTER(LEN=12) :: WHITE12
      LOGICAL :: EXIST
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IARA,IARA2
      REAL, ALLOCATABLE, DIMENSION(:) :: ARA,ARA2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DARA,DARA2
*----
*  RETRIEVING BLOCK TO BE COPIED
*----
      CALL LCMLEN(IPSTR,NAME,SIZE,ITYP)
      IF (SIZE.LE.0) THEN
         CALL LCMLIB(IPSTR)
         CALL XABORT('MSTCPB: INVALID BLOCK '//NAME//'.')
      ENDIF
      NARA=0
      IF (ITYP.EQ.1) THEN
         NARA=SIZE
         ALLOCATE(IARA(NARA))
         CALL LCMGET(IPSTR,NAME,IARA)
      ELSEIF (ITYP.LE.2) THEN
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
         CALL XABORT('MSTCPB: UNSUPPORTED TYPE')
      ENDIF
      IF (IEND.GT.NARA) CALL XABORT('MSTCPB: INCOMPATIBLE SIZE')
*     DOES THIS BLOCK ALREADY EXIST IN THE DESTINATION STRUCTURE ?
      EXIST=.FALSE.
      NELEO=0
*     SPECIAL CASE OF STATE-VECTOR MODIFICATION
      IF (NAME2.EQ.'STATE-VECTOR') THEN
         IF (IEND.GT.40)
     1      CALL XABORT('MSTCPM: STATE-VECTOR SIZE IS LIMITED TO 40.')
         IF (IEND.EQ.40)
     2      CALL XABORT('MSTCPM: 40th STATE-VECTOR ELEMENT SHOULD'//
     3      ' NOT BE MODIFIED.')
         ITYPO=1
         NELEO=NSTATE
         EXIST=.TRUE.
      ENDIF
      IF (NBLOCK.NE.0) THEN
         DO II=1,NBLOCK
            IF(BLNAM(II).EQ.NAME2) THEN
               ITYPO=BLTYP(II)
               IF (ITYPO.EQ.0)
     1           CALL XABORT('MSTCPM: '//NAME2//
     2                       ' IS AN EXISTING DIRECTORY.')
               IF (ITYPO.NE.ITYP)
     1           CALL XABORT('MSTCPM: INCOMPATIBLE TYPES')
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
         WRITE(IOUT,*) 'MSTCPB: BLOCK '//NAME//' IN UPDATE MODE'
      ELSE
*     NO: IT WILL BE CREATED
         WRITE(IOUT,*) 'MSTCPB: BLOCK '//NAME//' IN CREATION MODE'
      ENDIF
      ENDIF      
      NARA2=MAX(NARA,NELEO)
*     ALLOCATE MEMORY
      IF (ITYP.EQ.1) THEN
         SIZE2=NARA2
         ALLOCATE(IARA2(NARA2))
      ELSEIF (ITYP.EQ.2) THEN
         SIZE2=NARA2
         ALLOCATE(ARA2(NARA2))
      ELSEIF (ITYP.EQ.3) THEN
         SIZE2=3*NARA2
         ALLOCATE(IARA2(3*NARA2))
      ELSEIF (ITYP.EQ.4) THEN
         SIZE2=NARA2
         ALLOCATE(DARA2(NARA2))
      ENDIF
*     INITIALIZE BLOCK
      IF (EXIST) THEN
         IF (ITYP.EQ.1) THEN
            CALL LCMGET(IPSTR2,NAME2,IARA2)
         ELSEIF (ITYP.EQ.2) THEN
            CALL LCMGET(IPSTR2,NAME2,ARA2)
         ELSEIF (ITYP.EQ.3) THEN
            CALL LCMGET(IPSTR2,NAME2,IARA2)
         ELSEIF (ITYP.EQ.4) THEN
            CALL LCMGET(IPSTR2,NAME2,DARA2)
         ENDIF
      ELSE
         IF (ITYP.EQ.1) THEN
            IARA2(:NARA)=0
         ELSEIF (ITYP.EQ.2) THEN
            ARA2(:NARA)=0
         ELSEIF (ITYP.EQ.3) THEN
            WHITE12=' '
            DO II=1,NARA
               READ(WHITE12,'(3A4)') (IARA2(3*(II-1)+JJ),JJ=0,2)
            ENDDO
         ELSEIF (ITYP.EQ.4) THEN
            DARA2(:NARA)=0.D0
         ENDIF
      ENDIF
*     COPY ACTION
      DO II=IBEG,IEND,IINC
         IF (ITYP.EQ.1) THEN
            IARA2(II)=IARA(II)
         ELSEIF (ITYP.EQ.2) THEN
            ARA2(II)=ARA(II)
         ELSEIF (ITYP.EQ.3) THEN
            DO JJ=0,2
               IARA2(3*(II-1)+JJ)=IARA(3*(II-1)+JJ)
            ENDDO
         ELSEIF (ITYP.EQ.4) THEN
            DARA2(II)=DARA(II)
         ENDIF
      ENDDO
      IF (ITYP.EQ.1) THEN
         CALL LCMPUT(IPSTR2,NAME2,SIZE2,ITYP,IARA2)
         DEALLOCATE(IARA2)
         DEALLOCATE(IARA)
      ELSEIF (ITYP.EQ.2) THEN
         CALL LCMPUT(IPSTR2,NAME2,SIZE2,ITYP,ARA2)
         DEALLOCATE(ARA2)
         DEALLOCATE(ARA)
      ELSEIF (ITYP.EQ.3) THEN
         CALL LCMPUT(IPSTR2,NAME2,SIZE2,ITYP,IARA2)
         DEALLOCATE(IARA2)
         DEALLOCATE(IARA)
      ELSEIF (ITYP.EQ.4) THEN
         CALL LCMPUT(IPSTR2,NAME2,SIZE2,ITYP,DARA2)
         DEALLOCATE(DARA2)
         DEALLOCATE(DARA)
      ENDIF
*----
*  UPDATE NB. BLOCKS, REC-NAMES, REC-TYPES, REC-LENGTHS IN STATE-VECTOR
*  IF REQUIRED
*----
      IF (.NOT.EXIST) THEN
         CALL LCMGET(IPSTR2,'STATE-VECTOR',ISTATE)
         ISTATE(40)=ISTATE(40)+1
         BLNAM(NBLOCK+1)=NAME
         BLTYP(NBLOCK+1)=ITYP
         BLLEN(NBLOCK+1)=NARA
         CALL LCMPUT(IPSTR2,'STATE-VECTOR',NSTATE,1,ISTATE)
         CALL LCMPTC(IPSTR2,'REC-NAMES',12,NBLOCK+1,BLNAM)
         CALL LCMPUT(IPSTR2,'REC-TYPES',(NBLOCK+1),1,BLTYP)
         CALL LCMPUT(IPSTR2,'REC-LENGTHS',(NBLOCK+1),1,BLLEN)
      ENDIF
      RETURN
      END
