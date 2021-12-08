*DECK SOURISO
      SUBROUTINE SOURMONO(IPTRK,IPGEOM,NREG,LX,LY,LZ,NG,NUNS,NDIM,
     1 NSOUR,ISOUR,ISISOUR,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,XXX,YYY,ZZZ,
     2 MESH_LEN,BSINFO,BS,MAXL,NORM,DIR,MONOP,NBS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute moments of fixed isotropic sources.
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
      INTEGER NREG,LX,LY,LZ,NG,NUNS,NDIM,NSOUR,MESH_LEN(3),MONOP,
     1 ISISOUR(NG),BSINFO(3,NBS)
      REAL XMIN(NSOUR),XMAX(NSOUR),YMIN(NSOUR),YMAX(NSOUR),ZMIN(NSOUR),
     1 ZMAX(NSOUR),ISOUR(NG),SUNKNO(NUNS,NG),XXX(MESH_LEN(1)),
     2 YYY(MESH_LEN(2)),ZZZ(MESH_LEN(3)),NORM,DIR(3),BS(MAXL,NBS)
*----
*  LOCAL VARIABLES
*----

      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE),SPLIT_LEN(3),XP(NREG),YP(NREG),ZP(NREG),
     1 IPN
      REAL X(LX),Y(LX),Z(LX),IPVAL
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR) U_PTR,DU_PTR,DE_PTR,DZ_PTR
      INTEGER,ALLOCATABLE,DIMENSION(:) :: SPLITX,SPLITY,SPLITZ
      REAL, POINTER, DIMENSION(:) :: U,DU,DE,DZ
*----
*  RECOVER TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NLF=ISTATE(15)

      IF(NDIM.EQ.1) THEN
      CALL LCMGPD(IPTRK,'U',U_PTR)
      CALL C_F_POINTER(U_PTR,U,(/ NLF /))
      ELSE IF(NDIM.EQ.2) THEN
      CALL LCMLEN(IPTRK,'DU',NPQ,ITYLCM)
      CALL LCMGPD(IPTRK,'DU',DU_PTR)
      CALL LCMGPD(IPTRK,'DE',DE_PTR)
      CALL C_F_POINTER(DU_PTR,DU,(/ NPQ /))
      CALL C_F_POINTER(DE_PTR,DE,(/ NPQ /))
      ELSE IF(NDIM.EQ.3) THEN
      CALL LCMLEN(IPTRK,'DU',NPQ,ITYLCM)
      CALL LCMGPD(IPTRK,'DU',DU_PTR)
      CALL LCMGPD(IPTRK,'DE',DE_PTR)
      CALL LCMGPD(IPTRK,'DZ',DZ_PTR)
      CALL C_F_POINTER(DU_PTR,DU,(/ NPQ /))
      CALL C_F_POINTER(DE_PTR,DE,(/ NPQ /))
      CALL C_F_POINTER(DZ_PTR,DZ,(/ NPQ /))
      ENDIF

*----
*  RECOVER GEOMETRY INFORMATION
*----
      IF(NDIM.EQ.2) THEN
      IF(MONOP.NE.1.OR.MONOP.NE.-1) THEN     
        CALL LCMLEN(IPGEOM,'SPLITX',SPLIT_LEN(1),ITYLCM)
        ALLOCATE(SPLITX(SPLIT_LEN(1)))
        CALL LCMGET(IPGEOM,'SPLITX',SPLITX)
      ELSE IF(MONOP.NE.2.OR.MONOP.NE.-2) THEN
        CALL LCMLEN(IPGEOM,'SPLITY',SPLIT_LEN(2),ITYLCM)
        ALLOCATE(SPLITY(SPLIT_LEN(2)))
        CALL LCMGET(IPGEOM,'SPLITY',SPLITY)
      ENDIF
      ELSE IF(NDIM.EQ.3) THEN
      IF(MONOP.NE.1.OR.MONOP.NE.-1) THEN     
        CALL LCMLEN(IPGEOM,'SPLITX',SPLIT_LEN(1),ITYLCM)
        ALLOCATE(SPLITX(SPLIT_LEN(1)))
        CALL LCMGET(IPGEOM,'SPLITX',SPLITX)
      ELSE IF(MONOP.NE.2.OR.MONOP.NE.-2) THEN
        CALL LCMLEN(IPGEOM,'SPLITY',SPLIT_LEN(2),ITYLCM)
        ALLOCATE(SPLITY(SPLIT_LEN(2)))
        CALL LCMGET(IPGEOM,'SPLITY',SPLITY)
      ELSE IF(MONOP.NE.3.OR.MONOP.NE.-3) THEN
        CALL LCMLEN(IPGEOM,'SPLITZ',SPLIT_LEN(3),ITYLCM)
        ALLOCATE(SPLITZ(SPLIT_LEN(3)))
        CALL LCMGET(IPGEOM,'SPLITZ',SPLITZ)
      ENDIF 
      ENDIF

*----
*  1D CARTESIAN CASE
*----

      IF(NDIM.EQ.1) THEN

      ! Define the corresponding discrete ordinate  
      IPN=1
      IPVAL=(U(1)-DIR(1))**2
      DO IP=2,NLF
        IF((U(IP)-DIR(1))**2.LT.IPVAL) THEN
        IPN=IP
        IPVAL=(U(IP)-DIR(1))**2
        ENDIF
      ENDDO
 
      ! Save source information in BS variables
      IND=1
      DO IG=1,NG
        IF(ISISOUR(IG).EQ.1) THEN
        BSINFO(1,IND)=IG
        BSINFO(2,IND)=IPN
        BSINFO(3,IND)=MONOP
        BS(1,IND)=ISOUR(IG)
        IND=IND+1
        ENDIF
      ENDDO 
      
      ! Normalization
      NORM=SUM(ISOUR)

*----
*  2D CARTESIAN CASE
*----

      ELSE IF(NDIM.EQ.2) THEN

      CALL XABORT('SOUR: 2D CARTESIAN BOUNDARY SOURCES NOT IMPLEMENTED.'
     1 )

*----
*  3D CARTESIAN CASE
*----

      ELSE IF(NDIM.EQ.3) THEN

      CALL XABORT('SOUR: 3D CARTESIAN BOUNDARY SOURCES NOT IMPLEMENTED.'
     1 )

      ENDIF

      RETURN
      END
