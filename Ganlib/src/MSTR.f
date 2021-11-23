*DECK MSTR
      SUBROUTINE MSTR(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Manage user-defined structures.
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
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): modification type(VECTOR);
*         HENTRY(2): read-only type(VECTOR).
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER     NENTRY                  
      CHARACTER   HENTRY(NENTRY)*12
      INTEGER     IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR) KENTRY(NENTRY)
*----
*  LOCAL VARIABLES
*----
      INTEGER I,ACSTR,TYSTR,ACSTRO,TYSTRO,ACSTR2,TYSTR2
      INTEGER, PARAMETER :: IOUT=6
      INTEGER, PARAMETER :: NSTATE=40
      INTEGER INDIC,NITMA,ISTATE(NSTATE),IPRINT,NBLOCK,ILEN,ITYP,NBDIR
      REAL FLOTT
      DOUBLE PRECISION DFLOTT
      CHARACTER TEXT4*4,NAME*12,PATH*72,DIRS(37)*12,NAME2*12,TEXT12*12
      LOGICAL ROOT
      TYPE(C_PTR) IPSTR,IPSTRO,IPSTR2
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: BLNAM
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BLTYP,BLLEN
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.LT.1)
     1  CALL XABORT('MSTR: AT LEAST ONE PARAMETER EXPECTED.')
      DO I=1,NENTRY
         IF(IENTRY(I).GT.2) 
     1     CALL XABORT('MSTR: LINKED LIST OR XSM EXPECTED.')
      ENDDO
*----
*  CREATE STATE-VECTOR/SIGNATURE FOR STRUCTURES IN CREATION MODE
*  VERIFY IF IT EXISTS FOR STRUCTURES IN MODIFICATION MODE
*----
      DO I=1,NENTRY
         ACSTR=JENTRY(I)
         IPSTR=KENTRY(I)
         IF (ACSTR.EQ.0) THEN
*        THE STRUCTURE IS CREATED:
*        ASSIGN IT A DEFAULT SIGNATURE AND A STATE-VECTOR
            NAME='VECTOR'
            CALL LCMPTC(IPSTR,'SIGNATURE',12,1,NAME)
            ISTATE(:NSTATE)=0
            CALL LCMPUT(IPSTR,'STATE-VECTOR',NSTATE,1,ISTATE)
         ELSEIF (ACSTR.EQ.1) THEN
            CALL LCMLEN(IPSTR,'SIGNATURE',ILEN,ITYP)
            IF ((ILEN.NE.3).OR.(ITYP.NE.3))
     1         CALL XABORT('MSTR: INVALID SIGNATURE FOR '//HENTRY(I)
     2                                                  //'.')
            CALL LCMLEN(IPSTR,'STATE-VECTOR',ILEN,ITYP)
            IF ((ILEN.NE.40).OR.(ITYP.NE.1))
     1          CALL XABORT('MSTR: INVALID STATE-VECTOR FOR '//HENTRY(I)
     2                                                      //'.')
         ENDIF
      ENDDO
*---
* PROCESS USER'S INPUT 
*---
      IPRINT=0
      ACSTR=JENTRY(1)
      TYSTR=IENTRY(1)
      IPSTR=KENTRY(1)
      ACSTR2=ACSTR
      TYSTR2=TYSTR
      IPSTR2=IPSTR
 50   CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.EQ.10) GO TO 60
      IF(INDIC.NE.3) CALL XABORT('MSTR: CHARACTER DATA EXPECTED(1).')
      IF(TEXT4.EQ.'EDIT') THEN
*     MODULE EDITION LEVEL
         CALL REDGET(INDIC,IPRINT,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MSTR: INTEGER DATA EXPECTED(1).')
      ELSEIF(TEXT4.EQ.'TYPE') THEN
*     USER DEFINED TYPE FOR THE STRUCTURE
         CALL REDGET(INDIC,NITMA,FLOTT,PATH,DFLOTT)
         IF(INDIC.NE.3) 
     1   CALL XABORT('MSTR: CHARACTER DATA EXPECTED(2).')
         CALL MSTANP(NENTRY,IENTRY,JENTRY,KENTRY,PATH,IPSTR,
     1               ACSTR,TYSTR,NBDIR,DIRS,ROOT)
         IF (NBDIR.NE.1) CALL XABORT('MSTR: INVALID TYPE ENTRY.')
         CALL LCMPTC(IPSTR,'SIGNATURE',12,1,DIRS(1))
      ELSEIF((TEXT4.EQ.'PUT').OR.
     1       (TEXT4.EQ.'GET').OR.
     2       (TEXT4.EQ.'CP')) THEN
*     PUT, GET OR CP ACTION
         CALL REDGET(INDIC,NELEM,FLOTT,TEXT12,DFLOTT)
*        NUMBER OF ELEMENTS
         IF(INDIC.NE.1) CALL XABORT('MSTR: INTEGER DATA EXPECTED(2).')
         IBEG=1
         IEND=NELEM
         IINC=1
         CALL REDGET(INDIC,NITMA,FLOTT,PATH,DFLOTT)
         IF(INDIC.EQ.1) THEN
*        STARTING INDEX
            IBEG=NITMA
            CALL REDGET(INDIC,NITMA,FLOTT,PATH,DFLOTT)
            IF(INDIC.EQ.1) THEN
*           INCREMENT
               IINC=NITMA
            ELSE
               GOTO 10
            ENDIF
            CALL REDGET(INDIC,NITMA,FLOTT,PATH,DFLOTT)
         ENDIF
 10      CONTINUE
         IEND=IBEG+(NELEM-1)*IINC
         IF (INDIC.NE.3)
     1        CALL XABORT('MSTR: CHARACTER DATA EXPECTED(3).')
         IPSTRO=IPSTR
         ACSTRO=ACSTR
         TYSTRO=TYSTR
*        ANALYSE USER'S PATH
         CALL MSTANP(NENTRY,IENTRY,JENTRY,KENTRY,PATH,IPSTR,
     1        ACSTR,TYSTR,NBDIR,DIRS,ROOT)
*        GO TO REQUESTED DIRECTORY         
         IF (NBDIR.GT.1) THEN
            CALL MSTMOV(IPSTR,ACSTR,IPRINT,NBDIR-1,DIRS,ROOT)
         ENDIF
         NAME=DIRS(NBDIR)
         IF (TEXT4.EQ.'PUT') THEN
*        CREATING OR UPDATING DATA IN A BLOCK
            IF (ACSTR.EQ.2)
     1      CALL XABORT('MSTR: PUT NOT PERMITTED IN READ-ONLY MODE')
            CALL LCMGET(IPSTR,'STATE-VECTOR',ISTATE)
            NBLOCK=ISTATE(40)
            ALLOCATE(BLNAM(NBLOCK+1),BLTYP(NBLOCK+1),BLLEN(NBLOCK+1))
            IF (NBLOCK.GT.0) THEN
               CALL LCMGTC(IPSTR,'REC-NAMES',12,NBLOCK+1,BLNAM)
               CALL LCMGET(IPSTR,'REC-TYPES',BLTYP)
               CALL LCMGET(IPSTR,'REC-LENGTHS',BLLEN)
            ENDIF
            CALL MSTPUT(IPSTR,IPRINT,IBEG,IEND,IINC,NAME,NBLOCK,BLNAM,
     1           BLTYP,BLLEN)
            DEALLOCATE(BLLEN,BLTYP,BLNAM)
         ELSEIF(TEXT4.EQ.'GET') THEN
*        RETRIEVING DATA FROM A BLOCK
            CALL MSTGET(IPSTR,IPRINT,IBEG,IEND,IINC,NAME)
         ELSEIF(TEXT4.EQ.'CP') THEN
*        COPYING A BLOCK FROM ONE PLACE TO ANOTHER
            CALL REDGET(INDIC,NITMA,FLOTT,PATH,DFLOTT)
            IF (INDIC.NE.3)
     1        CALL XABORT('MSTR: CHARACTER DATA EXPECTED(4).')
*           ANALYSE USER'S PATH
            CALL MSTANP(NENTRY,IENTRY,JENTRY,KENTRY,PATH,IPSTR2,
     1           ACSTR2,TYSTR2,NBDIR,DIRS,ROOT)    
            IF (ACSTR2.EQ.2)
     1           CALL XABORT('MSTR: CP NOT PERMITTED IN READ-ONLY MODE')
*           GO TO REQUESTED DIRECTORY         
            IF (NBDIR.GT.1) THEN
               CALL MSTMOV(IPSTR2,ACSTR2,IPRINT,NBDIR-1,DIRS,ROOT)
            ENDIF
            NAME2=DIRS(NBDIR)
            CALL LCMGET(IPSTR2,'STATE-VECTOR',ISTATE)
            NBLOCK=ISTATE(40)
            ALLOCATE(BLNAM(NBLOCK+1),BLTYP(NBLOCK+1),BLLEN(NBLOCK+1))
            IF (NBLOCK.GT.0) THEN
               CALL LCMGTC(IPSTR2,'REC-NAMES',12,NBLOCK+1,BLNAM)
               CALL LCMGET(IPSTR2,'REC-TYPES',BLTYP)
               CALL LCMGET(IPSTR2,'REC-LENGTHS',BLLEN)
            ENDIF
            CALL MSTCPB(IPSTR,IPSTR2,IPRINT,IBEG,IEND,IINC,NAME,NAME2,
     1           NBLOCK,BLNAM,BLTYP,BLLEN)
            DEALLOCATE(BLLEN,BLTYP,BLNAM)
         ENDIF
         IPSTR=IPSTRO
         ACSTR=ACSTRO
         TYSTR=TYSTRO
      ELSEIF(TEXT4.EQ.'CD') THEN
*     CHANGING DIRECTORY
         CALL REDGET(INDIC,NITMA,FLOTT,PATH,DFLOTT)
         IF(INDIC.NE.3) 
     1   CALL XABORT('MSTR: CHARACTER DATA EXPECTED(5).')
*        ANALYSE USER'S PATH
         CALL MSTANP(NENTRY,IENTRY,JENTRY,KENTRY,PATH,IPSTR,
     1        ACSTR,TYSTR,NBDIR,DIRS,ROOT)
*        GO TO REQUESTED DIRECTORY   
         CALL MSTMOV(IPSTR,ACSTR,IPRINT,NBDIR,DIRS,ROOT)
      ELSEIF(TEXT4.EQ.';') THEN
         GOTO 60
      ELSE
         CALL XABORT('MSTR: '//TEXT4//' IS AN INVALID KEY WORD.')
      ENDIF
      GOTO 50
*
 60   CONTINUE
*
      RETURN
      END
