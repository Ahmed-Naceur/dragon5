*DECK NXTLDP
      SUBROUTINE NXTLDP(IPTRK,MAXMSH,IPIN,MESHP,NSURP,NREGP,DPMESH,
     1                  INDEX,IDREG,IDSUR)
*-----------------------------------------------------------------------
*
*Purpose:
* Load pin contents.
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
* IPIN    requested pin index.
*
*Parameters: output
* MESHP   pin meshes size.
* NSURP   number of surfaces for the pin.
* NREGP   number of regions for the pin.
* DPMESH  pin meshing vector.
* INDEX   pin index vector.
* IDREG   region index array.
* IDSUR   surface index array.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER MAXMSH,IPIN,MESHP(4),NSURP,NREGP,
     1 INDEX(5,-NSURP:NREGP),IDREG(NREGP),IDSUR(NSURP)
      DOUBLE PRECISION DPMESH(-1:MAXMSH,4)
*----
*  LOCAL VARIABLES
*----
      INTEGER NSTATE,IOUT
      PARAMETER(NSTATE=40,IOUT=6)
      INTEGER ESTATE(NSTATE)
      INTEGER IDIR
      CHARACTER NAMPIN*9,NAMREC*12
      CHARACTER CDIR(4)*1
      DATA CDIR /'X','Y','Z','R'/
*----
*  LOAD PIN RECORDS
*----
      WRITE(NAMPIN,'(A1,I8.8)') 'P',IPIN
      NAMREC=NAMPIN//'DIM'
      CALL XDISET(ESTATE,NSTATE,0)
      CALL LCMGET(IPTRK,NAMREC,ESTATE)
      MESHP(1)=ESTATE(3)
      MESHP(2)=ESTATE(4)
      MESHP(3)=ESTATE(5)
      MESHP(4)=ESTATE(2)
      NREGP=ESTATE(8)
      NSURP=ESTATE(9)
      NAMREC=NAMPIN//'RID'
      CALL LCMGET(IPTRK,NAMREC,IDREG)
      NAMREC=NAMPIN//'SID'
      CALL LCMGET(IPTRK,NAMREC,IDSUR)
      NAMREC=NAMPIN//'VSI'
      CALL LCMGET(IPTRK,NAMREC,INDEX)
      DO IDIR=1,4
        NAMREC=NAMPIN//'SM'//CDIR(IDIR)
        IF(MESHP(IDIR) .GT. 0) THEN
          CALL LCMGET(IPTRK,NAMREC,DPMESH(-1,IDIR))
        ENDIF
      ENDDO
*
      RETURN
      END
