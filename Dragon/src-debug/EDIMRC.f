*DECK EDIMRC
      SUBROUTINE EDIMRC(IPTRK ,IPRINT  ,NREGIO, NMERGE, IMERGE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To generate the region merging index for homogenisation
* per CELL for NXT treated geometry
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau.
*
*Parameters: input
* IPTRK   pointer to the TRACKING data structure.
* IPRINT  print level.
* NREGIO  number of regions.
*
*Parameters: output
* NMERGE  final number of merged regions.
* IMERGE  merged region index.
*
*-----------------------------------------------------------------------
*
      USE         GANLIB
      IMPLICIT    NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR) IPTRK
      INTEGER     IPRINT
      INTEGER     NREGIO
      INTEGER     NMERGE
      INTEGER     IMERGE(NREGIO)
*----
*  Local parameters
*----
      CHARACTER   NAMSBR*6
      PARAMETER  (NAMSBR='EDIMRC')
      INTEGER     NSTATE
      PARAMETER  (NSTATE=40)
      INTEGER     ILCMUP,ILCMDN
      PARAMETER  (ILCMUP=1,ILCMDN=2)
*----
*  Local variables
*----
      INTEGER     ISTATE(NSTATE),IEDIMG(NSTATE)
      CHARACTER   HSIGN*12
      INTEGER     NDIM,NFREG,NFSUR,NNC,NBUCEL,NUCELL(3),MAXREG
*----
*  Test if valid tracking data structure
*  EXCELL with type 4 tracking
*----
      CALL LCMGTC(IPTRK,'SIGNATURE',12,1,HSIGN)
*----
*  TEST IF GEOMETRY OR EXCELL TRACK DATA STRUCTURE
*----
      IF(HSIGN .NE. 'L_TRACK     ') CALL XABORT(NAMSBR//
     >': Invalid data structure for merge by cell')
      CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,HSIGN)
      IF((HSIGN .NE. 'EXCELL') .AND. (HSIGN .NE. 'MCCG')) THEN
        CALL XABORT(NAMSBR//': Invalid tracking for merge by cell')
      ENDIF
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      IF(ISTATE(7) .NE. 4) CALL XABORT(NAMSBR//
     >': Only NXT tracking permitted for merge by cell')
      IF(ISTATE(40) .EQ. 1) CALL XABORT(NAMSBR//
     >': Double heterogeneity (Bihet) not implemented')
      CALL LCMSIX(IPTRK,'NXTRecords  ',ILCMUP)
      CALL LCMGET(IPTRK,'G00000001DIM',IEDIMG)
      NDIM     =IEDIMG( 1)
      NNC      =IEDIMG( 4)
      NBUCEL   =IEDIMG( 5)
      NUCELL(1)=IEDIMG(13)
      NUCELL(2)=IEDIMG(14)
      NUCELL(3)=IEDIMG(15)
      NFSUR    =IEDIMG(22)
      NFREG    =IEDIMG(23)
      MAXREG   =IEDIMG(25)
      CALL EDIMCN(IPTRK ,IPRINT,NDIM  ,NUCELL,NBUCEL,MAXREG,
     >            NFREG,NFSUR,NNC,NREGIO,NMERGE,IMERGE)
      CALL LCMSIX(IPTRK,'NXTRecords  ',ILCMDN)
      RETURN
      END
