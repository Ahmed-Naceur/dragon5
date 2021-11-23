*DECK NXTMCA
      SUBROUTINE NXTMCA(IPTRK)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Add MC: specific geometry analysis info to NXTRecords.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): Romain Le Tellier
*
*Parameters: input
* IPTRK   pointer to the Tracking data structure.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
*----
*  LOCAL VARIABLES
*----
      INTEGER NSTATE
      PARAMETER(NSTATE=40)
      INTEGER NFREG,NMIX,NFSUR,NDIM,NBUCEL,NUCELL(3),MAXREG,NBTCLS,
     1 MAXPIN,MAXMSP,MAXRSP,MXGSUR,MXGREG,NUNK
      INTEGER GSTATE(NSTATE),ESTATE(NSTATE)
      CHARACTER NAMREC*12,CDIR(4)*1
      DATA CDIR /'X','Y','Z','R'/
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IUNFLD
*----
*  RECOVER SOME BASIC NXT GEOMETRY ANALYSIS INFO AND ALLOCATE RELATED
*  MEMORY
*----     
      CALL XDISET(GSTATE,NSTATE,0)
      CALL LCMGET(IPTRK,'STATE-VECTOR',GSTATE)
      NFREG    =GSTATE( 1)
      NMIX     =GSTATE( 4)
      NFSUR    =GSTATE( 5)
      IF (GSTATE(7).NE.4)
     1 CALL XABORT('NXTMCA: ONLY NXT: GEOMETRY ANALYSIS IS PERMITTED')
      CALL LCMSIX(IPTRK,'NXTRecords',1)
      CALL XDISET(ESTATE,NSTATE,0)
      CALL LCMGET(IPTRK,'G00000001DIM',ESTATE)
      NDIM     =ESTATE( 1)
      NBUCEL   =ESTATE( 5)
      NUCELL(1)=ESTATE(13)
      NUCELL(2)=ESTATE(14)
      NUCELL(3)=ESTATE(15)
      MAXREG   =ESTATE(17)
      NBTCLS   =ESTATE(18)
      MAXPIN   =ESTATE(19)
      MAXMSP   =ESTATE(20)
      MAXRSP   =ESTATE(21)
      IF (NFSUR.NE.ESTATE(22))
     1 CALL XABORT('NXTMCA: INCONSISTENT NUMBER OF OUTER SURFACES')
      IF (NFREG.NE.ESTATE(23))
     1 CALL XABORT('NXTMCA: INCONSISTENT NUMBER OF REGIONS')
      MXGSUR   =ESTATE(24)
      MXGREG   =ESTATE(25)
      NUNK=NFSUR+NFREG+1
*     cell index and orientation for the cells filling the geometry
      ALLOCATE(IUNFLD(2*NBUCEL))
      NAMREC='G00000001CUF'
      CALL LCMGET(IPTRK,NAMREC,IUNFLD)
*----
*  ADD MCA: SPECIFIC GEOMETRY ANALYSIS INFO TO NXTRecords
*----
      CALL NXTMCB(IPTRK,NUCELL,MXGSUR,MXGREG,MAXPIN,IUNFLD)
*----
*  RELEASE MEMORY
*----
      DEALLOCATE(IUNFLD)
*     
      CALL LCMSIX(IPTRK,' ',2)
      
      RETURN
      END
