*DECK MCTLDC
      SUBROUTINE MCTLDC(IPTRK,IPRINT,NDIM,MAXMSH,MXGSUR,MXGREG,NBIND,
     1           ICEL,INDX,ADDREC,DGMESH,IDIRC,MESHC,NSURC,NREGC,NTPIN,
     2           CELLPO,DCMESH,INDEX,ITPIN,DRAPIN,ID)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Load cell contents.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* IPTRK   pointer to the TRACKING data structure.
* IPRINT  print level.
* NDIM    problem dimensions.
* MAXMSH  maximum number of elements in MESH array.
* MXGSUR  maximum number of surfaces for any geometry.
* MXGREG  maximum number of region for any geometry.
* NBIND   first dimension of INDX.
* ICEL    requested cell index.
* INDX    position index in the geometry structure.
* ADDREC  name of additional requested record.
* DGMESH  meshing vector for global geometry.
*
*Parameters: output
* IDIRC   cylinders orientation.
* MESHC   cell meshes size.
* NSURC   number of surfaces for the cell.
* NREGC   number of regions for the cell.
* NTPIN   number of pins within the cell.
* CELLPO  cell global coordinates.
* DCMESH  cell meshing vector.
* INDEX   cell index vector.
* ID      additional requested record.
* ITPIN   undefined.
* DRAPIN  undefined.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER IPRINT,NDIM,MAXMSH,MXGSUR,MXGREG,NBIND,ICEL,INDX(NBIND),
     1 IDIRC,MESHC(4),NSURC,NREGC,NTPIN,INDEX(5,-MXGSUR:MXGREG),
     2 ITPIN(3,NTPIN),ID(MXGREG)
      DOUBLE PRECISION DGMESH(-1:MAXMSH,4),CELLPO(3,2),
     1 DCMESH(-1:MAXMSH,4),DRAPIN(-1:4,NTPIN)
      CHARACTER ADDREC*3
*----
*  LOCAL VARIABLES
*----
      INTEGER NSTATE,IOUT
      PARAMETER(NSTATE=40,IOUT=6)
      INTEGER ESTATE(NSTATE)
      INTEGER IDIR,ITYPG
      CHARACTER NAMCEL*9,NAMREC*12
      CHARACTER CDIR(4)*1
      DATA CDIR /'X','Y','Z','R'/
*----
*  LOAD CELL RECORDS
*----
      IF(IPRINT.GT.50) THEN
        WRITE(6,'(/21H MCTLDC: PROCESS CELL,I6)') ICEL
      ENDIF
      WRITE(NAMCEL,'(A1,I8.8)') 'C',ICEL
      NAMREC=NAMCEL//'DIM'
      CALL XDISET(ESTATE,NSTATE,0)
      CALL LCMGET(IPTRK,NAMREC,ESTATE)      
      ITYPG=ESTATE(1)
      MESHC(1)=ESTATE(3)
      MESHC(2)=ESTATE(4)
      MESHC(3)=ESTATE(5)
      MESHC(4)=ESTATE(2)
      NREGC=ESTATE(39)!8
      NSURC=ESTATE(40)!9
      NTPIN=ESTATE(18)
      NAMREC=NAMCEL//ADDREC
      CALL LCMGET(IPTRK,NAMREC,ID)
      NAMREC=NAMCEL//'VSC'
      CALL LCMGET(IPTRK,NAMREC,INDEX(1,-NSURC))
      DO IDIR=1,4
        NAMREC=NAMCEL//'SM'//CDIR(IDIR)
        IF(MESHC(IDIR) .GT. 0) THEN
          CALL LCMGET(IPTRK,NAMREC,DCMESH(-1,IDIR))
        ENDIF
      ENDDO
      IF(NTPIN .GT .0) THEN
        NAMREC=NAMCEL//'PIN'
        CALL LCMGET(IPTRK,NAMREC,DRAPIN)
        NAMREC=NAMCEL//'PNT'
        CALL LCMGET(IPTRK,NAMREC,ITPIN)
      ENDIF
      CELLPO(1,2)=DGMESH(INDX(1),1)
      CELLPO(1,1)=DGMESH(INDX(1)-1,1)
      CELLPO(2,2)=DGMESH(INDX(2),2)
      CELLPO(2,1)=DGMESH(INDX(2)-1,2)
      IF(NDIM .EQ. 3) THEN
         CELLPO(3,2)=DGMESH(INDX(3),3)
         CELLPO(3,1)=DGMESH(INDX(3)-1,3) 
      ENDIF      
      IF(ITYPG .EQ. 20 .OR. ITYPG .EQ. 21 .OR. 
     >   ITYPG .EQ. 22 .OR. ITYPG .EQ. 23) THEN
        IF(ITYPG .EQ. 21 ) THEN
          IDIRC=1
        ELSE IF(ITYPG .EQ. 22) THEN
          IDIRC=2
        ELSE
          IDIRC=3
        ENDIF
      ELSE
         IDIRC=0
      ENDIF
*
      RETURN
      END
