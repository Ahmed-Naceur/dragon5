*DECK SOURISO
      SUBROUTINE SOURISO(IPTRK,IPGEOM,NREG,NG,NANIS,NUNS,NDIM,NSOUR,
     1   ISOUR,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the fixed isotropic sources.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): C. Bienvenue 
*
*Parameters: input
*
*Parameters: output
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPGEOM
      INTEGER NREG,NG,NUNS,NDIM,NSOUR,NANIS
      REAL XMIN(NSOUR),XMAX(NSOUR),YMIN(NSOUR),YMAX(NSOUR),ZMIN(NSOUR),
     1 ZMAX(NSOUR),ISOUR(NG),SUNKNO(NUNS,NG)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE),MESH_LEN(3),SPLIT_LEN(3)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER,ALLOCATABLE,DIMENSION(:) :: MESHX,MESHY,MESHZ,SPLITX,
     1 SPLITY,SPLITZ
*----
*  RECOVER TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      IELEM=ISTATE(8)
      ISCAT=ISTATE(16)
*----
*  RECOVER GEOMETRY INFORMATION
*----
      IF(NDIM.EQ.1) THEN
      CALL LCMLEN(IPGEOM,'MESHX',MESH_LEN(1),ITYLCM)
      CALL LCMLEN(IPGEOM,'SPLITX',SPLIT_LEN(1),ITYLCM)
      ALLOCATE(MESHX(MESH_LEN(1)),SPLITX(SPLIT_LEN(1)))
      CALL LCMGET(IPGEOM,'MESHX',MESHX)
      CALL LCMGET(IPGEOM,'SPLITX',SPLITX)
      ELSE IF(NDIM.EQ.2) THEN

      ELSE IF(NDIM.EQ.3) THEN

      ENDIF
*----
*  CAS 1D
*----

      IF(NDIM.EQ.1) THEN
      NSCT=ISCAT
      DO 30 IG=1,NG 
      DO 20 IR=1,NREG
      DO 10 N=1,NSOUR
*      IF(XMIN(N).LE.XXX(IR).AND.XMAX(N).GE.XXX(IR)) THEN
*      IND=(IR-1)*NSCT*IELEM+1
*      SUNKNO(IND,IG)=SUNKNO(IND,IG)+ISOUR(NG)
*      ENDIF
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE     
      ENDIF










      END
