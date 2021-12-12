*DECK NXTACG
      SUBROUTINE NXTACG(IPGEO ,IPTRK ,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To analyze an assembly of cells containing
* clusters using the new EXCELL tracking procedure.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPGEO   pointer to the GEOMETRY data structure in
*         read only mode.
* IPTRK   pointer to the TRACKING data structure in
*         update or creation mode.
* IPRINT  print level.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPGEO,IPTRK
      INTEGER          IPRINT
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTACG')
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
*----
*  Local variables
*----
      INTEGER          NDIM,ITYPBC,ITYGEO,IDIRG
      INTEGER          ISTATG(NSTATE),ISTATT(NSTATE),IEDIMG(NSTATE)
      INTEGER          ILEAK,IDIAG,ISAXIS(3),NBOCEL,NBUCEL,
     >                 NOCELL(3),NUCELL(3),MAXCEL
      INTEGER          IDIR,IMCELL,ILCELL,MAXPIN,MAXMSP,
     >                 MAXMSH,MAXREG,NBTCLS
      INTEGER          NFSUR,NFREG,MXGSUR,MXGREG,IANIS,NBUNK
      INTEGER          NEREG,NESUR
      CHARACTER        NAMREC*12
*----
*  Update for Hexagonal geometry
*----
      INTEGER          IHSYM
*----
*  Update for prismatic geometry
*----
      INTEGER          IPRISM
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDIRR,ITURN,MERGE,KEYMRG,
     > MATRT,MATRTN,KEYFLX
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IDFEX,IDFRT,ITSYM,IUNFLD,
     > NAGGEO,MATALB,MATCOD
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLUME
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SURVOL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DGMESH
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DCMESH
*----
*  Processing starts:
*  print routine opening header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 10) WRITE(IOUT,6000) NAMSBR
      IF(IPRINT .GE. 1) WRITE(IOUT,6010)
*----
*  Get state vectors
*----
      CALL XDISET(ISTATG,NSTATE,0)
      CALL LCMGET(IPGEO,'STATE-VECTOR',ISTATG)
      CALL XDISET(ISTATT,NSTATE,0)
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATT)
*----
*  Get geometry state vector and test if geometry
*  is valid for NXT.
*  Valid geometries are:
*  Cartesian Boundary (ITYPBC=0):
*    CAR2D (5), CAR3D (7),
*    CARCEL (20), CARCELX (21), CARCELY (22), CARCELZ (23)
*  Annular Boundary (ITYPBC=1):
*    TUBE (3), TUBEX (10), TUBEY (11), TUBEZ (6)
*  Hexagonal Boundary (ITYPBC=2):
*    HEX (8) , HEXZ (9), HEXT (12) , HEXTZ (13)
*  Initialize
*  NDIM  :  number of dimensions for problem
*  IDIRG  : first direction to process
*  ITYPBC : type of boundary
*----
      NDIM=2
      IDIRG=3
      ITYPBC=0
      ITYGEO=ISTATG(1)
      IF(ITYGEO .EQ. 5) THEN
        IDIRG=1
      ELSE IF(ITYGEO .EQ. 7) THEN
        NDIM=3
        IDIRG=1
      ELSE IF(ITYGEO .EQ. 20) THEN
        IDIRG=1
      ELSE IF(ITYGEO .EQ. 21) THEN
        NDIM=3
        IDIRG=2
      ELSE IF(ITYGEO .EQ. 22) THEN
        NDIM=3
      ELSE IF(ITYGEO .EQ. 23) THEN
        NDIM=3
        IDIRG=1
      ELSE IF(ITYGEO .EQ. 3) THEN
        IDIRG=1
        ITYPBC=1
      ELSE IF(ITYGEO .EQ. 6) THEN
        NDIM=3
        IDIRG=1
        ITYPBC=1
      ELSE IF(ITYGEO .EQ. 10) THEN
        NDIM=3
        IDIRG=2
        ITYPBC=1
      ELSE IF(ITYGEO .EQ. 11) THEN
        NDIM=3
        ITYPBC=1
      ELSE IF(ITYGEO .EQ. 8) THEN
        IDIRG=1
        ITYPBC=2
      ELSE IF(ITYGEO .EQ. 9) THEN
        NDIM=3
        IDIRG=1
        ITYPBC=2
      ELSE IF(ITYGEO .EQ. 12) THEN
        IDIRG=1
        ITYPBC=2
      ELSE IF(ITYGEO .EQ. 13) THEN
        NDIM=3
        IDIRG=1
        ITYPBC=2
      ELSE IF(ITYGEO .EQ. 26) THEN
        IDIRG=1
        ITYPBC=2
      ELSE IF(ITYGEO .EQ. 27) THEN
        NDIM=3
        IDIRG=1
        ITYPBC=2
      ELSE
        CALL XABORT(NAMSBR//
     >  ': Module cannot analyze this geometry')
      ENDIF
      IMCELL=ISTATG(8)
      IF(IMCELL .EQ. 0) THEN
*----
*  Pure geometry
*  Create assembly of one cell
*----
        NOCELL(1)=1
        NOCELL(2)=1
        NOCELL(3)=1
        IF(NDIM .LT. 3) NOCELL(3)=0
        ILCELL=0
      ELSE
*----
*  Assembly
*----
        NOCELL(1)=ISTATG(3)
        NOCELL(2)=ISTATG(4)
        NOCELL(3)=ISTATG(5)
        ILCELL=1
      ENDIF
      IPRISM=ISTATT(39)
*----
*  Read and analyze boundary conditions
*----
      CALL NXTBCG(IPGEO ,IPTRK ,IPRINT,NDIM  ,ITYPBC,IDIRG ,
     >            IDIAG ,ISAXIS,IHSYM ,ILEAK ,IPRISM)
      IF(IDIAG .NE. 0) THEN
        IF(NOCELL(1) .NE. NOCELL(2)) CALL XABORT(NAMSBR//
     >  ': DIAG requires symmetric X-Y mesh.')
      ENDIF
      CALL XDISET(IEDIMG,NSTATE,0)
      IEDIMG(1)=NDIM
      IEDIMG(2)=ITYPBC
      IEDIMG(3)=IDIRG
*----
*  Compute global mesh in each direction and
*  total number of cells after unfolding of geometry
*----
      IF(ITYPBC .EQ. 0) THEN
* Full Cartesian
        IF(IDIAG .EQ. 0) THEN
          NBOCEL=1
          DO IDIR=1,3
            NBOCEL=NBOCEL*MAX(NOCELL(IDIR),1)
          ENDDO
        ELSE
          NBOCEL=MAX(NOCELL(1),1)
          NBOCEL=(NBOCEL*(NBOCEL+1))/2
          NBOCEL=NBOCEL*MAX(NOCELL(3),1)
        ENDIF
        NBUCEL=1
        MAXCEL=1
        DO IDIR =1,3
          IF(ABS(ISAXIS(IDIR)) .EQ. 1) THEN
            NUCELL(IDIR)=2*NOCELL(IDIR)-1
          ELSE IF (ABS(ISAXIS(IDIR)) .EQ. 2) THEN
            NUCELL(IDIR)=2*NOCELL(IDIR)
          ELSE
            NUCELL(IDIR)=NOCELL(IDIR)
          ENDIF
          NBUCEL=NBUCEL*MAX(NUCELL(IDIR),1)
          MAXCEL=MAX(MAXCEL,NUCELL(IDIR))
        ENDDO
      ELSE IF(ITYPBC .EQ. 1) THEN
* Annular or cylindrical
        IDIR=MOD(IDIRG+1,3)+1
        NBOCEL=MAX(NOCELL(IDIR),1)
        NBUCEL=NBOCEL
        MAXCEL=1
        DO IDIR =1,3
          NUCELL(IDIR)=NOCELL(IDIR)
          MAXCEL=MAX(MAXCEL,NUCELL(IDIR))
        ENDDO
      ELSE IF(ITYPBC .EQ. 2) THEN
* Hexagons
        NBOCEL=MAX(NOCELL(3),1)*MAX(NOCELL(1),1)
        NBUCEL=NBOCEL
        MAXCEL=1
        DO IDIR =1,3
          NUCELL(IDIR)=NOCELL(IDIR)
          MAXCEL=MAX(MAXCEL,NUCELL(IDIR))
        ENDDO
      ENDIF
*----
*  Create Array for testing symmetry and unfolding the
*  geometry according to symmetries
*----
      IEDIMG(4)=NBOCEL
      IEDIMG(5)=NBUCEL
      IEDIMG(6)=IDIAG
      IEDIMG(7)=ISAXIS(1)
      IEDIMG(8)=ISAXIS(2)
      IEDIMG(9)=ISAXIS(3)
      ALLOCATE(IDFEX(11,NBOCEL),IDFRT(8,NBOCEL),ITSYM(4,NBOCEL),
     > IUNFLD(2,NBUCEL))
      IF(ITYPBC .EQ. 0) THEN
        CALL NXTCUA(IPRINT,NDIM  ,IDIAG ,ISAXIS,
     >              NBOCEL,NBUCEL,NOCELL,NUCELL,ITSYM ,
     >              IDFEX ,IDFRT ,IUNFLD)
      ELSE IF(ITYPBC .EQ. 1) THEN
        CALL XABORT(NAMSBR//
     >  ': Annular boundary not programmed yet')
*        CALL NXTAUA(IPRINT,NDIM  ,IDIAG ,ISAXIS,
*     >              NBOCEL,NBUCEL,NOCELL,NUCELL,ITSYM ,
*     >              IDFEX ,IDFRT ,IUNFLD)
      ELSE IF(ITYPBC .EQ. 2) THEN
        CALL NXTHUA(IPRINT,NDIM  ,IHSYM ,ISAXIS,
     >              NBOCEL,NBUCEL,NOCELL,NUCELL,ITSYM ,
     >              IDFEX ,IDFRT ,IUNFLD)
      ENDIF
*----
*  Save cell unfolding and internal symmetrization
*  vectors
*----
      CALL LCMSIX(IPTRK,'NXTRecords  ',1)
      NAMREC='G00000001CUF'
      CALL LCMPUT(IPTRK,NAMREC,2*NBUCEL,1,IUNFLD)
      NAMREC='G00000001CIS'
      CALL LCMPUT(IPTRK,NAMREC,4*NBOCEL,1,ITSYM)
      NAMREC='G00000001CFE'
      CALL LCMPUT(IPTRK,NAMREC,11*NBOCEL,1,IDFEX)
      IEDIMG(10)=NOCELL(1)
      IEDIMG(11)=NOCELL(2)
      IEDIMG(12)=NOCELL(3)
      IEDIMG(13)=NUCELL(1)
      IEDIMG(14)=NUCELL(2)
      IEDIMG(15)=NUCELL(3)
      IEDIMG(16)=MAXCEL
*----
*  Get maximum dimensions and geometry limits
*  Test if assembly can be built and save global mesh for geometry.
*  Test if cluster pins are valid.
*----
      ALLOCATE(NAGGEO(3,NBOCEL),IDIRR(NBOCEL),ITURN(NBOCEL),
     > MERGE(NBOCEL))
      ALLOCATE(DCMESH(3,2,NBOCEL))
      ALLOCATE(DGMESH((MAXCEL+1),3))
      CALL NXTGMD(IPGEO ,IPTRK ,IPRINT,ITYPBC,ILCELL,NSTATE,
     >            NBOCEL,NBUCEL,MAXCEL,NUCELL,IUNFLD,IEDIMG,
     >            NAGGEO,ITURN ,MERGE ,IDIRR ,DCMESH,DGMESH)
      DEALLOCATE(DGMESH)
*----
*  Allocate memory to read global information
*  for each geometry and subgeometry
*----
      MAXMSH=IEDIMG(16)
      MAXREG=IEDIMG(17)
      NBTCLS=IEDIMG(18)
      MAXPIN=IEDIMG(19)
      MAXMSP=MAX(IEDIMG(20),IEDIMG(16),1)
*----
*  Create multicell description of geometry
*----
      CALL NXTMCD(IPGEO ,IPTRK ,IPRINT,NDIM  ,ILCELL,NBOCEL,
     >            MAXMSH,MAXREG,MAXPIN,NBTCLS,ITSYM ,IDFEX ,
     >            DCMESH,NAGGEO,ITURN ,IDIRR ,
     >            NFSUR ,NFREG ,MXGSUR,MXGREG)
      IEDIMG(22)=NFSUR
      IEDIMG(23)=NFREG
      IEDIMG(24)=MXGSUR
      IEDIMG(25)=MXGREG
      CALL LCMPUT(IPTRK,'G00000001DIM',NSTATE,1,IEDIMG)
*----
*  Compute surfaces and volumes and define MATALB
*----
      ALLOCATE(KEYMRG(-NFSUR:NFREG),MATALB(-NFSUR:NFREG,2))
      ALLOCATE(SURVOL(-NFSUR:NFREG))
      CALL NXTCVS(IPTRK ,IPRINT,NDIM  ,ITYPBC,NBOCEL,NFSUR ,
     >            NFREG ,MXGSUR,MXGREG,KEYMRG,MATALB,SURVOL)
*----
*  Create BC-REFL+TRAN vector
*----
      ALLOCATE(MATRT(NFSUR))
      CALL NXTBRT(IPTRK ,IPRINT,NDIM  ,ITYPBC,ISAXIS,NBOCEL,
     >            MAXMSP,MAXPIN,NFSUR ,MXGSUR,MXGREG,IDFRT ,
     >            MATRT )
      CALL LCMSIX(IPTRK,'NXTRecords  ',2)
*----
*  Save global tracking information
*  Including KEYMRG
*----
      ALLOCATE(MATCOD(NFREG,2),KEYFLX(NFREG),VOLUME(NFREG),
     > MATRTN(NFSUR))
      CALL NXTAGM(IPRINT,NFSUR ,NFREG ,NEREG ,NESUR ,KEYMRG,
     >            MATALB,MATRT ,SURVOL,KEYFLX,MATCOD,MATRTN,
     >            VOLUME)
      CALL LCMPUT(IPTRK,'BC-REFL+TRAN',NESUR,1,MATRTN)
      CALL LCMPUT(IPTRK,'MATCOD      ',NEREG,1,MATCOD(1,1))
      CALL LCMPUT(IPTRK,'HOMMATCOD   ',NEREG,1,MATCOD(1,2))
      CALL LCMPUT(IPTRK,'VOLUME      ',NEREG,2,VOLUME)
      CALL LCMPUT(IPTRK,'KEYFLX      ',NEREG,1,KEYFLX)
      DEALLOCATE(MATRTN,VOLUME,KEYFLX,MATCOD)
      DEALLOCATE(MATRT)
      DEALLOCATE(SURVOL,MATALB,KEYMRG)
      NBUNK=0
      DO IANIS=0,ISTATT(6)-1
        NBUNK=NBUNK+2*IANIS+1
      ENDDO
      NBUNK=NEREG*NBUNK
      ISTATT(1)=NEREG
      ISTATT(2)=NBUNK
      ISTATT(3)=ILEAK
      ISTATT(4)=ISTATG(7)
      ISTATT(5)=NESUR
      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,ISTATT)
*----
*  Release memory
*----
      DEALLOCATE(DCMESH,MERGE,ITURN,IDIRR,NAGGEO)
      DEALLOCATE(IUNFLD,ITSYM,IDFRT,IDFEX)
*----
*  Processing finished:
*  print routine closing header if required
*  and return
*----
      IF(IPRINT .GE. 1)  WRITE(IOUT,6011) NFREG,NEREG,NFSUR,NESUR
      IF(IPRINT .GE. 10) WRITE(IOUT,6001) NAMSBR
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(/' Analysis of geometry begins'/)
 6011 FORMAT(' Number of regions before merge =',I10/
     >       ' Number of regions after merge  =',I10/
     >       ' Number of surfaces before merge=',I10/
     >       ' Number of surfaces after merge =',I10//
     >       ' Analysis of geometry completed:'/
     >       )
      END
