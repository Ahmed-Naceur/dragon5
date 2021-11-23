*DECK LCMGCD
      SUBROUTINE LCMGCD(IPLIST,NAMP,ILONG,HDATA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Copy a character variable from a table into memory.
*
*Copyright:
* Copyright (C) 2000 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIST  address of the table.
* NAMP    character*12 name of the existing block.
* ILONG   dimension of the character variable.
*
*Parameters: output
* HDATA   character variable.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIST
      INTEGER ILONG
      CHARACTER*(*) NAMP,HDATA(ILONG)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPLIST
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDATA
*
      ILE1=(LEN(HDATA(1))+3)/4
      CALL LCMLEN(IPLIST,NAMP,JLONG,ITYLCM)
      IF(ITYLCM.NE.10) CALL XABORT('LCMGCD: LIST EXPECTED.')
      IF(JLONG.GT.ILONG) CALL XABORT('LCMGCD: HDATA OVERFLOW.')
      JPLIST=LCMGID(IPLIST,NAMP)
      DO ISET=1,JLONG
         CALL LCMLEL(JPLIST,ISET,ILE2,ITYLCM)
         IF(ITYLCM.NE.3) CALL XABORT('LCMGCD: CHARACTER EXPECTED.')
         ALLOCATE(IDATA(ILE2))
         CALL LCMGDL(JPLIST,ISET,IDATA)
         HDATA(ISET)=' '
         WRITE(HDATA(ISET),'(100A4)') (IDATA(I),I=1,MIN(ILE1,ILE2))
         DEALLOCATE(IDATA)
      ENDDO
      DO ISET=JLONG+1,ILONG
         HDATA(ISET)=' '
      ENDDO
      RETURN
      END
