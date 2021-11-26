*DECK MSTANP
      SUBROUTINE MSTANP(NENTRY,IENTRY,JENTRY,KENTRY,PATH,IPSTR,ACSTR,
     1           TYSTR,NBDIR,DIRS,ROOT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Analyse user defined path.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): R. Le Tellier
*
*Parameters: input
* NENTRY  number of LCM objects or files used by the operator.
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
* PATH    user defined path.
*
*Parameters: output
* IPSTR   structure address.
* ACSTR   structure access.
* TYSTR   structure type.
* NBDIR   number of directories/blocks in PATH.
* DIRS    array of the directories/blocks names.
* ROOT    flag to know if the path is relative or absolute.
*
*-----------------------------------------------------------------------
*
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NENTRY                  
      CHARACTER DIRS(37)*12,PATH*72
      INTEGER IENTRY(NENTRY),JENTRY(NENTRY),ACSTR,TYSTR,NBDIR
      TYPE(C_PTR) KENTRY(NENTRY),IPSTR
      LOGICAL ROOT
*----
*  LOCAL VARIABLES
*----
      INTEGER II,IBEG,IEND,I
      CHARACTER*12 OBJNAM,MYDIR
*
*     SEARCH FOR OBJECT NAME
      OBJNAM=' '
      IBEG=0
      DO II=1,72
         IF (PATH(II:II).EQ.':') THEN
            OBJNAM=PATH(1:II-1)
            IBEG=II+1
            GOTO 10
         ENDIF
      ENDDO
 10   CONTINUE
      IF (OBJNAM.EQ.' ') THEN
         IBEG=1
      ELSE
         READ(OBJNAM,'(I12)') I
         IF ((I.GT.NENTRY).OR.(I.LT.1)) THEN
           CALL XABORT('MSTANP: INCORRECT OBJECT INDEX')
         ENDIF
         ACSTR=JENTRY(I)
         TYSTR=IENTRY(I)
         IPSTR=KENTRY(I)
         GOTO 15
      ENDIF
 15   CONTINUE
*     FIND THE HIERCHICAL DIRECTORIES STRUCTURE
      IF (PATH(IBEG:IBEG).EQ.'/') THEN
         ROOT=.TRUE.
      ELSE
         ROOT=.FALSE.
      ENDIF
      NBDIR=0
      DO II=IBEG,72
         IF ((PATH(II:II).EQ.'/').OR.
     1       (PATH(II:II).EQ.' ')) THEN
            IEND=II-1
            IF ((IBEG.LE.IEND).AND.(PATH(IBEG:IEND).NE.'.')) THEN
               NBDIR=NBDIR+1
               MYDIR=PATH(IBEG:IEND)
               DIRS(NBDIR)=MYDIR
            ENDIF
            IBEG=II+1
            IF (PATH(II:II).EQ.' ') GOTO 20
         ENDIF
      ENDDO
 20   CONTINUE
*
      RETURN
      END
