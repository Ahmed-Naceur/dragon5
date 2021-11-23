*DECK MCTLDP
      SUBROUTINE MCTLDP(IPTRK,IPRINT,MAXMSH,MXGSUR,MXGREG,ITPIN,ADDREC,
     1           MESHP,NSURP,NREGP,DPMESH,INDEX,ID)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Load pin contents.
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
* MAXMSH  maximum number of elements in MESH array.
* MXGSUR  maximum number of surfaces for any geometry.
* MXGREG  maximum number of region for any geometry.
* ITPIN   pin index.
* ADDREC  name of additional requested record.
*
*Parameters: output
* MESHP   pin meshes size.
* NSURP   number of surfaces for the pin.
* NREGP   number of regions for the pin.
* DPMESH  pin meshing vector.
* INDEX   pin index vector.
* ID      additional requested record.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER IPRINT,MAXMSH,MXGSUR,MXGREG,ITPIN,MESHP(4),NSURP,NREGP,
     1 INDEX(5,-MXGSUR:MXGREG),ID(MXGREG)
      DOUBLE PRECISION DPMESH(-1:MAXMSH,4)
      CHARACTER ADDREC*3
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
      IF(IPRINT.GT.50) THEN
        WRITE(6,'(/20H MCTLDP: PROCESS PIN,I6)') ITPIN
      ENDIF
      WRITE(NAMPIN,'(A1,I8.8)') 'P',ITPIN
      NAMREC=NAMPIN//'DIM'
      CALL XDISET(ESTATE,NSTATE,0)
      CALL LCMGET(IPTRK,NAMREC,ESTATE)
      MESHP(1)=ESTATE(3)
      MESHP(2)=ESTATE(4)
      MESHP(3)=ESTATE(5)
      MESHP(4)=ESTATE(2)
      NREGP=ESTATE(39)!8
      NSURP=ESTATE(40)!9
      NAMREC=NAMPIN//ADDREC
      CALL LCMGET(IPTRK,NAMREC,ID)
      NAMREC=NAMPIN//'VSC'
      CALL LCMGET(IPTRK,NAMREC,INDEX(1,-NSURP))
      DO IDIR=1,4
         NAMREC=NAMPIN//'SM'//CDIR(IDIR)
         IF(MESHP(IDIR) .GT. 0) THEN
            CALL LCMGET(IPTRK,NAMREC,DPMESH(-1,IDIR))
         ENDIF
      ENDDO
*
      RETURN
      END
