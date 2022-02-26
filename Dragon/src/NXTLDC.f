*DECK NXTLDC
      SUBROUTINE NXTLDC(IPTRK,MAXMSH,ICEL,IDIRC,MESHC,NSURC,NREGC,
     1           NTPIN,DCMESH,INDEX,IDREG,IDSUR,ITPIN,DRAPIN)
*-----------------------------------------------------------------------
*
*Purpose:
* Load cell contents.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): R. Le Tellier
*
*Parameters: input
* IPTRK   pointer to the TRACKING data structure.
* MAXMSH  maximum number of elements in MESH array.
* ICEL    requested cell index.
*
*Parameters: output
* IDIRC   cylinders orientation.
* MESHC   cell meshes size.
* NSURC   number of surfaces for the cell.
* NREGC   number of regions for the cell.
* NTPIN   number of pins within the cell.
* DCMESH  cell meshing vector.
* INDEX   cell index vector.
* IDREG   region index array.
* IDSUR   surface index array.
* ITPIN   pin integer descriptor.
* DRAPIN  pin double descriptor.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER MAXMSH,ICEL,IDIRC,MESHC(4),NSURC,NREGC,NTPIN,
     1 INDEX(5,-NSURC:NREGC),IDREG(NREGC),IDSUR(NSURC),ITPIN(3,NTPIN)
      DOUBLE PRECISION DCMESH(-1:MAXMSH,4),DRAPIN(-1:4,NTPIN)
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
      WRITE(NAMCEL,'(A1,I8.8)') 'C',ICEL
      NAMREC=NAMCEL//'DIM'
      CALL XDISET(ESTATE,NSTATE,0)
      CALL LCMGET(IPTRK,NAMREC,ESTATE)
      ITYPG=ESTATE(1)
      MESHC(1)=ESTATE(3)
      MESHC(2)=ESTATE(4)
      MESHC(3)=ESTATE(5)
      MESHC(4)=ESTATE(2)
      NREGC=ESTATE(8)
      NSURC=ESTATE(9)
      NTPIN=ESTATE(18)
      NAMREC=NAMCEL//'RID'
      CALL LCMGET(IPTRK,NAMREC,IDREG)
      NAMREC=NAMCEL//'SID'
      CALL LCMGET(IPTRK,NAMREC,IDSUR)
      NAMREC=NAMCEL//'VSI'
      CALL LCMGET(IPTRK,NAMREC,INDEX)
      DO IDIR=1,4
        NAMREC=NAMCEL//'SM'//CDIR(IDIR)
        IF(MESHC(IDIR).GT.0) THEN
          CALL LCMGET(IPTRK,NAMREC,DCMESH(-1,IDIR))
        ENDIF
      ENDDO
      IF(NTPIN.GT.0) THEN
        NAMREC=NAMCEL//'PIN'
        CALL LCMGET(IPTRK,NAMREC,DRAPIN)
        NAMREC=NAMCEL//'PNT'
        CALL LCMGET(IPTRK,NAMREC,ITPIN)
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
