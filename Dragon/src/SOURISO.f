*DECK SOURISO
      SUBROUTINE SOURISO(IPTRK,IPGEOM,NREG,NG,NANIS,NUNS,NDIM,NSOUR,
     1 ISOUR,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,XXX,YYY,ZZZ,MESH_LEN,SUNKNO)
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
      INTEGER NREG,NG,NUNS,NDIM,NSOUR,NANIS,MESH_LEN(3)
      REAL XMIN(NSOUR),XMAX(NSOUR),YMIN(NSOUR),YMAX(NSOUR),ZMIN(NSOUR),
     1 ZMAX(NSOUR),ISOUR(NG),SUNKNO(NUNS,NG),XXX(MESH_LEN(1)),
     2 YYY(MESH_LEN(2)),ZZZ(MESH_LEN(3))        
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE),SPLIT_LEN(3),XP(NREG)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER,ALLOCATABLE,DIMENSION(:) :: SPLITX,SPLITY,SPLITZ
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
      CALL LCMLEN(IPGEOM,'SPLITX',SPLIT_LEN(1),ITYLCM)
      ALLOCATE(SPLITX(SPLIT_LEN(1)))
      CALL LCMGET(IPGEOM,'SPLITX',SPLITX)
      ELSE IF(NDIM.EQ.2) THEN

      ELSE IF(NDIM.EQ.3) THEN

      ENDIF
*----
*  CAS 1D
*----

      IF(NDIM.EQ.1) THEN

      K=1
      DO I=1,SPLIT_LEN(1)
      DO J=1,SPLITX(I)
      XP(K)=I
      K=K+1
      ENDDO
      ENDDO

      NSCT=ISCAT
      DO 30 IG=1,NG 
      DO 20 IR=1,NREG
      STEP=(XXX(XP(IR)+1)-XXX(XP(IR)))/SPLITX(XP(IR))
      IF(XP(IR).EQ.1) THEN
        X=XXX(XP(IR))+0.5*STEP+STEP*(IR-1)
      ELSE
        X=XXX(XP(IR))+0.5*STEP+STEP*(IR-SPLITX(XP(IR)-1)-1)
      ENDIF
      DO 10 N=1,NSOUR
      IF(XMIN(N).LE.X.AND.XMAX(N).GE.X) THEN
      IND=(IR-1)*NSCT*IELEM+1
      SUNKNO(IND,IG)=SUNKNO(IND,IG)+ISOUR(IG)
      ENDIF
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE     
      ENDIF

      RETURN
      END
