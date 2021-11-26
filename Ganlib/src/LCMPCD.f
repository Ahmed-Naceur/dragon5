*DECK LCMPCD
      SUBROUTINE LCMPCD(IPLIST,NAMP,ILONG,HDATA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Copy a character variable from memory into a table.
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
* NAMP    character*12 name of the block.
* ILONG   dimension of the character variable.
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
      CHARACTER*(*) NAMP,HDATA(*)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPLIST
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDATA
*
      ILEN=(LEN(HDATA(1))+3)/4
      IF(ILEN.GT.100) CALL XABORT('LCMPCD: ILEN OVERFLOW.')
      ALLOCATE(IDATA(ILEN))
      JPLIST=LCMLID(IPLIST,NAMP,ILONG)
      DO ISET=1,ILONG
         READ(HDATA(ISET),'(100A4)') (IDATA(I),I=1,ILEN)
         CALL LCMPDL(JPLIST,ISET,ILEN,3,IDATA)
      ENDDO
      DEALLOCATE(IDATA)
      RETURN
      END
