*DECK MCTPSP
      SUBROUTINE MCTPSP(IPTRK,POS,IREG,IEV)
*-----------------------------------------------------------------------
*
*Purpose:
* Store position and region index in TRACKING table for PSP display.
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
* POS     point global coordinates.
* IREG    region/surface index.
* IEV     event index.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER IREG,IEV
      DOUBLE PRECISION POS(3)
*----
*  LOCAL VARIABLES
*----
      INTEGER ILONG,ITYLCM,NPOINT
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: REGI,EVENT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: COORD
*
      CALL LCMSIX(IPTRK,'MCpoints',1) 
      CALL LCMLEN(IPTRK,'REGI',ILONG,ITYLCM)
      NPOINT=ILONG+1
      ALLOCATE(COORD(3,NPOINT),REGI(NPOINT),EVENT(NPOINT))
      IF (ILONG.GT.0) THEN
         CALL LCMGET(IPTRK,'COORD',COORD)
         CALL LCMGET(IPTRK,'REGI',REGI)
         CALL LCMGET(IPTRK,'EVENT',EVENT)
      ENDIF
      COORD(1,NPOINT)=POS(1)
      COORD(2,NPOINT)=POS(2)
      COORD(3,NPOINT)=POS(3)
      REGI(NPOINT)=IREG
      EVENT(NPOINT)=IEV
      CALL LCMPUT(IPTRK,'COORD',3*NPOINT,4,COORD)
      CALL LCMPUT(IPTRK,'REGI',NPOINT,1,REGI)
      CALL LCMPUT(IPTRK,'EVENT',NPOINT,1,EVENT)
      DEALLOCATE(EVENT,REGI,COORD)
      CALL LCMSIX(IPTRK,' ',2)
*
      RETURN
      END
