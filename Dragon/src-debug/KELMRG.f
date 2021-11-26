*DECK KELMRG
      FUNCTION    KELMRG(IPGEOM, NSURO, NVOLO, IDLGEO, MATGEO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Merge zones for a heterogeneous block.
*
*Copyright:
* Copyright (C) 1990 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPGEOM  pointer to the geometry. 
* NSURO   number of surfaces for a specific geometry.
* NVOLO   number of zones for a specific geometry.
* IDLGEO  specific position for a geometry.
*
*Parameters: output
* MATGEO  numbering of zones and surfaces for all geometries.
* KELMRG  number of surfaces and zones renumbered.
*
*-----------------------------------------------------------------------
*
      USE                 GANLIB
      IMPLICIT            NONE
      LOGICAL             SWONCE
      TYPE(C_PTR)         IPGEOM
      INTEGER             KELMRG,NSURO,NVOLO,IDLGEO,MATGEO(*)
      INTEGER             IOUT, IND, MATMIN, MATMAX, ITYLCM, IMRG, JMRG,
     >                    ILEN, I
      PARAMETER         ( IOUT=6 )
*
      IND(I)= IDLGEO + I
      CALL LCMLEN(IPGEOM, 'MERGE', ILEN, ITYLCM)
      IF( ILEN.EQ.0 )THEN
         KELMRG= NVOLO - NSURO + 1
      ELSE
         IF( ILEN.GT.NVOLO )
     >      CALL XABORT('KELMRG: MERGING HAS TOO MANY ZONES' )
         CALL LCMGET(IPGEOM, 'MERGE', MATGEO(IND(1)) )
         MATMIN= 100000000
         MATMAX=-100000000
         DO 10 IMRG= 1, ILEN
            IF( MATGEO(IND(IMRG)).LT.MATMIN) MATMIN= MATGEO(IND(IMRG))
            IF( MATGEO(IND(IMRG)).GT.MATMAX) MATMAX= MATGEO(IND(IMRG))
   10    CONTINUE
         IF( MATMIN.NE.1 )
     >      CALL XABORT('KELMRG: NO FIRST MERGING ZONE' )
         DO 30 JMRG= MATMIN, MATMAX
            SWONCE= .FALSE.
            DO 20 IMRG= 1, ILEN
               SWONCE= SWONCE.OR.(MATGEO(IND(IMRG)).EQ.JMRG)
   20       CONTINUE
            IF( .NOT.SWONCE )THEN
               WRITE(IOUT,*) 'WHERE IS MERGE REGION NO.', JMRG
               CALL XABORT('KELMRG: ERROR IN MERGE NUMBERING' )
            ENDIF
   30    CONTINUE
         KELMRG= MATMAX - NSURO + 1
      ENDIF
*
      RETURN
      END
