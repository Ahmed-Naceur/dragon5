*DECK MSTMOV
      SUBROUTINE MSTMOV(IPSTR,ACSTR,IPRINT,NBDIR,DIRS,ROOT)
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
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* IPSTR   structure address.
* ACSTR   structure access.
* IPRINT  level of print index.
* NBDIR   number of successive directories.
* DIRS    array of directories names.
* ROOT    flag to know if the path is relative or absolute.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPRINT,ACSTR,NBDIR
      CHARACTER DIRS(NBDIR)*12
      TYPE(C_PTR) IPSTR
      LOGICAL ROOT
*----
*  LOCAL VARIABLES
*----
      INTEGER I,NBLOCK
      INTEGER, PARAMETER :: IOUT=6
      INTEGER, PARAMETER :: NSTATE=40
      INTEGER ISTATE(NSTATE)
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: BLNAM
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BLTYP,BLLEN
*
      IF (ROOT) THEN
*        FIRST OF ALL, GOING TO ROOT DIR IF THE PATH IS ABSOLUTE
         IF (IPRINT.GT.2) WRITE(IOUT,*) 'MSTMOV: GOING TO ROOT DIR'
         CALL LCMSIX(IPSTR,' ',0)
      ENDIF
      DO I=1,NBDIR
*     ENTERING SUCCESSIVE DIRECTORIES
*     (A DIRECTORY IS CREATED IF IT DOES NOT EXIST
*      AND THE STRUCTURE IS NOT IN READ-ONLY MODE)
         CALL LCMGET(IPSTR,'STATE-VECTOR',ISTATE)
         NBLOCK=ISTATE(40)
         ALLOCATE(BLNAM(NBLOCK+1),BLTYP(NBLOCK+1),BLLEN(NBLOCK+1))
         IF (NBLOCK.GT.0) THEN
            CALL LCMGTC(IPSTR,'REC-NAMES',12,NBLOCK+1,BLNAM)
            CALL LCMGET(IPSTR,'REC-TYPES',BLTYP)
            CALL LCMGET(IPSTR,'REC-LENGTHS',BLLEN)
         ENDIF
         CALL MSTCDI(IPSTR,ACSTR,IPRINT,NBLOCK,DIRS(I),BLNAM,BLTYP,
     1   BLLEN)
         DEALLOCATE(BLLEN,BLTYP,BLNAM)
      ENDDO
      RETURN
      END
