*DECK NXTCVM
      SUBROUTINE NXTCVM(IFTRK ,IPRINT,NFREG ,NFSUR ,NEREG ,NESUR ,
     >                  MATALB,SURVOL,KEYMRG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To compress VOLSUR and MATALB according to KEYMRG
* and save on IFTRK.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IFTRK   pointer to the TRACKING file in creation mode.
* IPRINT  print level.
* NFREG   number of regions (geometry).
* NFSUR   number of surfaces (geometry).
* NEREG   number of regions (compress).
* NESUR   number of surfaces (compress).
* MATALB  global mixture/albedo identification vector (geometry).
* SURVOL  global surface volume vector (geometry).
* KEYMRG  index array for surface and volume renumbering.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IFTRK,IPRINT
      INTEGER          NFREG,NFSUR,NEREG,NESUR
      INTEGER          MATALB(-NFSUR:NFREG),KEYMRG(-NFSUR:NFREG)
      DOUBLE PRECISION SURVOL(-NFSUR:NFREG)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTCVM')
*----
*  Local variables
*----
      INTEGER          IREG,JREG,IMIX,ITST,JJ
      DOUBLE PRECISION DVR
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ALBMAT
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLSUR
*----
*  Scratch storage allocation
*   ALBMAT  global mixture/albedo identification vector (compress).
*   VOLSUR  global surface volume vector (compress).
*----
      ALLOCATE(ALBMAT(-NESUR:NEREG),VOLSUR(-NESUR:NEREG))
*----
*  Processing starts:
*  print routine opening header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,*) 'Surface Merge',NFSUR
        WRITE(IOUT,'(5I10)') (KEYMRG(JREG),JREG=-1,-NFSUR,-1)
        WRITE(IOUT,*) 'Region Merge',NFREG
        WRITE(IOUT,'(5I10)') (KEYMRG(JREG),JREG=1,NFREG)
      ENDIF
*----
*  Compress regions
*----
      ALBMAT(0)=0
      DVR=0.0D0
      IMIX=0
      VOLSUR(0)=0.0
      DO IREG=1,NEREG
        ITST=-1
        DO JREG=1,NFREG
*          write(6,*) 'Merging regions',IREG,JREG,KEYMRG(JREG)
          IF(KEYMRG(JREG) .EQ. IREG) THEN
            IF(ITST .EQ. -1) THEN
              IMIX=MATALB(JREG)
              DVR=SURVOL(JREG)
              ITST=1
            ELSE
              IF(IMIX .NE. MATALB(JREG) ) CALL XABORT(NAMSBR//
     >': Merging region with different mixtures not permitted')
              DVR=DVR+SURVOL(JREG)
            ENDIF
          ENDIF
        ENDDO
        IF(ITST .EQ. -1) CALL XABORT(NAMSBR//
     >': One merge region not defined')
        VOLSUR(IREG)=REAL(DVR)
        ALBMAT(IREG)=IMIX
      ENDDO
*----
*  Compress surfaces
*----
      DO IREG=-1,-NESUR,-1
        ITST=-1
        DO JREG=-1,-NFSUR,-1
*          write(6,*) 'Merging surfaces',IREG,JREG,KEYMRG(JREG)
          IF(KEYMRG(JREG) .EQ. IREG) THEN
            IF(ITST .EQ. -1) THEN
              IMIX=MATALB(JREG)
              DVR=SURVOL(JREG)
              ITST=1
            ELSE
              IF(IMIX .NE. MATALB(JREG) ) CALL XABORT(NAMSBR//
     >': Merging surfaces with different albedos not permitted')
              DVR=DVR+SURVOL(JREG)
            ENDIF
          ENDIF
        ENDDO
        IF(ITST .EQ. -1) CALL XABORT(NAMSBR//
     >': One merge surface not defined')
        VOLSUR(IREG)=REAL(DVR/4.0D0)
        ALBMAT(IREG)=IMIX
      ENDDO
      WRITE(IFTRK) (VOLSUR(JJ),JJ=-NESUR,NEREG)
      WRITE(IFTRK) (ALBMAT(JJ),JJ=-NESUR,NEREG)
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(VOLSUR,ALBMAT)
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
