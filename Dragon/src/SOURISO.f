*DECK SOURISO
      SUBROUTINE SOURISO(IPTRK,IPGEOM,NREG,LX,LY,LZ,NG,NUNS,NDIM,
     1 NSOUR,ISOUR,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,XXX,YYY,ZZZ,MESH_LEN,
     2 SUNKNO)
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
      INTEGER NREG,LX,LY,LZ,NG,NUNS,NDIM,NSOUR,MESH_LEN(3)
      REAL XMIN(NSOUR),XMAX(NSOUR),YMIN(NSOUR),YMAX(NSOUR),ZMIN(NSOUR),
     1 ZMAX(NSOUR),ISOUR(NG),SUNKNO(NUNS,NG),XXX(MESH_LEN(1)),
     2 YYY(MESH_LEN(2)),ZZZ(MESH_LEN(3))        
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE),SPLIT_LEN(3),XP(NREG),YP(NREG),ZP(NREG)
      REAL X(LX),Y(LX),Z(LX)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER,ALLOCATABLE,DIMENSION(:) :: SPLITX,SPLITY,SPLITZ
*----
*  RECOVER TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      ITYPE=ISTATE(6)
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
      CALL LCMLEN(IPGEOM,'SPLITX',SPLIT_LEN(1),ITYLCM)
      CALL LCMLEN(IPGEOM,'SPLITY',SPLIT_LEN(2),ITYLCM)
      CALL LCMLEN(IPGEOM,'SPLITZ',SPLIT_LEN(3),ITYLCM)
      ALLOCATE(SPLITX(SPLIT_LEN(1)),SPLITY(SPLIT_LEN(2)),
     1 SPLITZ(SPLIT_LEN(3)))
      CALL LCMGET(IPGEOM,'SPLITX',SPLITX)
      CALL LCMGET(IPGEOM,'SPLITY',SPLITY)
      CALL LCMGET(IPGEOM,'SPLITZ',SPLITZ)
      ENDIF
*----
*  1D CARTESIAN CASE
*----

      IF(NDIM.EQ.1) THEN

      ! CALCULATE X-COORDINATES OF EACH VOXELS
      K=1
      DO I=1,SPLIT_LEN(1)
      DO J=1,SPLITX(I)
      XP(K)=I
      K=K+1
      ENDDO
      ENDDO

      DO 10 IX=1,LX
      STEPX=(XXX(XP(IX)+1)-XXX(XP(IX)))/SPLITX(XP(IX))
      IF(XP(IX).EQ.1) THEN
        X(IX)=XXX(XP(IX))+0.5*STEPX+STEPX*(IX-1)
      ELSE
        X(IX)=XXX(XP(IX))+0.5*STEPX+STEPX*(IX-SPLITX(XP(IX)-1)-1)
      ENDIF
   10 CONTINUE

      ! CALCULATE THE SOURCE DENSITY
      NSCT=ISCAT
      DO 40 IX=1,LX
      IR=IX
      DO 30 N=1,NSOUR
      IF(XMIN(N).LE.X(IX).AND.XMAX(N).GE.X(IX)) THEN
      IND=(IR-1)*NSCT*IELEM+1
      DO 20 IG=1,NG
      SUNKNO(IND,IG)=SUNKNO(IND,IG)+ISOUR(IG)
   20 CONTINUE
      ENDIF
   30 CONTINUE
   40 CONTINUE

*----
*  2D CARTESIAN CASE 
*----

      ELSE IF(NDIM.EQ.2) THEN
      
      ! CALCULATE XY-COORDINATES OF EACH VOXELS
     
      K=1
      DO I=1,SPLIT_LEN(1)
      DO J=1,SPLITX(I)
      XP(K)=I
      K=K+1
      ENDDO
      ENDDO
      K=1
      DO I=1,SPLIT_LEN(2)
      DO J=1,SPLITY(I)
      YP(K)=I
      K=K+1
      ENDDO
      ENDDO

      DO 100 IX=1,LX
      STEPX=(XXX(XP(IX)+1)-XXX(XP(IX)))/SPLITX(XP(IX))
      IF(XP(IX).EQ.1) THEN
        X(IX)=XXX(XP(IX))+0.5*STEPX+STEPX*(IX-1)
      ELSE
        X(IX)=XXX(XP(IX))+0.5*STEPX+STEPX*(IX-SPLITX(XP(IX)-1)-1)
      ENDIF
  100 CONTINUE

      DO 110 IY=1,LY
      STEPY=(YYY(YP(IY)+1)-YYY(YP(IY)))/SPLITY(YP(IY))
      IF(YP(IY).EQ.1) THEN
        Y(IY)=YYY(YP(IY))+0.5*STEPY+STEPY*(IY-1)
      ELSE
        Y(IY)=YYY(YP(IY))+0.5*STEPY+STEPY*(IY-SPLITY(YP(IY)-1)-1)
      ENDIF
  110 CONTINUE

      ! CALCULATE THE SOURCE DENSITY
      NSCT=ISCAT*(ISCAT+1)/2
      DO 150 IY=1,LY
      DO 140 IX=1,LX
      IR=IX+(IY-1)*LX
      DO 130 N=1,NSOUR
      IF(XMIN(N).LE.X(IX).AND.XMAX(N).GE.X(IX).AND.
     1   YMIN(N).LE.Y(IY).AND.YMAX(N).GE.Y(IY)) THEN
      IND=(IR-1)*NSCT*IELEM*IELEM+1
      DO 120 IG=1,NG
      SUNKNO(IND,IG)=SUNKNO(IND,IG)+ISOUR(IG)
  120 CONTINUE
      ENDIF
  130 CONTINUE
  140 CONTINUE
  150 CONTINUE     

*----
*  3D CARTESIAN CASE
*----

      ELSE IF(NDIM.EQ.3) THEN
      
      ! CALCULATE XYZ-COORDINATES OF EACH VOXELS 
      K=1
      DO I=1,SPLIT_LEN(1)
      DO J=1,SPLITX(I)
      XP(K)=I
      K=K+1
      ENDDO
      ENDDO
      K=1
      DO I=1,SPLIT_LEN(2)
      DO J=1,SPLITY(I)
      YP(K)=I
      K=K+1
      ENDDO
      ENDDO
      K=1
      DO I=1,SPLIT_LEN(3)
      DO J=1,SPLITZ(I)
      ZP(K)=I
      K=K+1
      ENDDO
      ENDDO

      DO 200 IX=1,LX
      STEPX=(XXX(XP(IX)+1)-XXX(XP(IX)))/SPLITX(XP(IX))
      IF(XP(IX).EQ.1) THEN
        X(IX)=XXX(XP(IX))+0.5*STEPX+STEPX*(IX-1)
      ELSE
        X(IX)=XXX(XP(IX))+0.5*STEPX+STEPX*(IX-SPLITX(XP(IX)-1)-1)
      ENDIF
  200 CONTINUE

      DO 210 IY=1,LY
      STEPY=(YYY(YP(IY)+1)-YYY(YP(IY)))/SPLITY(YP(IY))
      IF(YP(IY).EQ.1) THEN
        Y(IY)=YYY(YP(IY))+0.5*STEPY+STEPY*(IY-1)
      ELSE
        Y(IY)=YYY(YP(IY))+0.5*STEPY+STEPY*(IY-SPLITY(YP(IY)-1)-1)
      ENDIF
  210 CONTINUE

      DO 220 IZ=1,LZ   
      STEPZ=(ZZZ(ZP(IZ)+1)-ZZZ(ZP(IZ)))/SPLITZ(ZP(IZ))
      IF(ZP(IZ).EQ.1) THEN
        Z(IZ)=ZZZ(ZP(IZ))+0.5*STEPZ+STEPZ*(IZ-1)
      ELSE
        Z(IZ)=ZZZ(ZP(IZ))+0.5*STEPZ+STEPZ*(IZ-SPLITZ(ZP(IZ)-1)-1)
      ENDIF        
  220 CONTINUE

      ! CALCULATE THE SOURCE DENSITY
      NSCT=(ISCAT)**2
      DO 270 IZ=1,LZ
      DO 260 IY=1,LY
      DO 250 IX=1,LX
      IR=IX+(IY-1)*LX+(IZ-1)*LX*LY
      DO 240 N=1,NSOUR
      IF(XMIN(N).LE.X(IX).AND.XMAX(N).GE.X(IX).AND.
     1   YMIN(N).LE.Y(IY).AND.YMAX(N).GE.Y(IY).AND.
     2   ZMIN(N).LE.Z(IZ).AND.ZMAX(N).GE.Z(IZ)) THEN
      IND=(IR-1)*NSCT*IELEM*IELEM*IELEM+1
      DO 230 IG=1,NG
      SUNKNO(IND,IG)=SUNKNO(IND,IG)+ISOUR(IG)
  230 CONTINUE
      ENDIF
  240 CONTINUE
  250 CONTINUE
  260 CONTINUE
  270 CONTINUE     

      ELSE
      CALL XABORT('SOUR: INVALID GEOMETRY, ONLY 1D, 2D AND 3D CARTESIAN'
     1 //' GEOMETRY ARE ACTUALLY IMPLEMENTED.')
      ENDIF


      RETURN
      END
