*DECK AXGGEO
      SUBROUTINE AXGGEO(IPGEOM,IPTRKM,IPRINT,GEONAM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Generate temporary tracking file to be used by PSPTRK.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy and G. Marleau
*
*Parameters: input
* IPGEOM  geometry data structures pointer
* IPTRKM  tracking data structures pointer         
* IPRINT  print level                              
* GEONAM  geometry name
*
*----
*
      USE          GANLIB
      IMPLICIT     NONE
      INTEGER      IOUT,NSTATE
      CHARACTER    NAMSBR*6
      PARAMETER   (IOUT=6,NSTATE=40,
     >             NAMSBR='AXGGEO')
*----
*  ROUTINE PARAMETERS
*----
      TYPE(C_PTR)  IPGEOM,IPTRKM
      INTEGER      IPRINT
      CHARACTER    GEONAM*12
*----
*  LOCAL PARAMETERS
*---- 
      INTEGER     ISTATE(NSTATE)
      INTEGER     ITYPEG,ITGEO
      CHARACTER   HSIGN*12
      INTEGER     NV,NS,NSOUT,NREG,NUNK,ICODE(6)
      REAL        EXTKOP(NSTATE)
      INTEGER     ITROP,MAXMIX,IREG,ISYMM
      INTEGER     IUEXP,KDROPN,KDRCLS,IRC 
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEYMRG,MATALB,MATMRG
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLSUR,VOLMRG
*----
*  STORE SIGNATURE AND TRACK TYPE ON IPTRKM
*----
      HSIGN='L_TRACK     '
      CALL LCMPTC(IPTRKM,'SIGNATURE',12,1,HSIGN)
      HSIGN='EXCELL      '
      CALL LCMPTC(IPTRKM,'TRACK-TYPE',12,1,HSIGN)
*----
*  ANALYZE GEOMETRY
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      ITYPEG= ISTATE(1)
      ITROP = 0
      IF(ITYPEG .EQ. 3 .OR. ITYPEG .EQ. 6 ) THEN
        ITGEO= 1
      ELSE IF(ITYPEG .EQ.  8 .OR. ITYPEG .EQ.  9 .OR.
     >        ITYPEG .EQ. 24 .OR. ITYPEG .EQ. 25      ) THEN
        ITGEO= 2
      ELSE IF(ITYPEG .EQ.  5 .OR. ITYPEG .EQ.  7 .OR.
     >        ITYPEG .EQ. 20 .OR. ITYPEG .EQ. 21 .OR.
     >        ITYPEG .EQ. 22 .OR. ITYPEG .EQ. 23      ) THEN
        ITGEO= 3
      ELSE
        ITGEO= 0
      ENDIF
      IF(ISTATE(13) .GE. 1) THEN
*----
*  CLUSTER GEOMETRY
*----     
        ISYMM=1
        CALL AXGXCW(IPGEOM ,IPTRKM,IPRINT,GEONAM,ISYMM )
        ITROP=3
      ELSE IF(ITGEO .EQ. 2 ) THEN
*----
*  HEXAGONAL 2D GEOMETRIES
*----
*          CALL AXGXHX(IPGEOM ,IPTRKM,IPRINT,GEONAM)
        ITROP=2
      ELSE IF(ITGEO .EQ. 3 ) THEN
*----
*  CARTESIAN 2D/3D ASSEMBLIES
*  CALL XELPRP TO GET GEOMETRY DIMENSIONING INFORMATION
*----
        CALL AXGXEL(IPGEOM ,IPTRKM,IPRINT,GEONAM)
        ITROP=1
      ELSE
        CALL XABORT(NAMSBR//': INVALID TYPE OF GEOMETRY')
      ENDIF
      CALL LCMGET(IPTRKM,'ICODE       ',ICODE)
      CALL LCMSIX(IPTRKM,'EXCELL      ',1)
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPTRKM,'STATE-VECTOR',ISTATE) 
      NV=ISTATE(3)
      NS=ISTATE(2)
      NUNK=NV+NS+1
      ALLOCATE(KEYMRG(NUNK),MATALB(NUNK),VOLSUR(NUNK))
      CALL LCMGET(IPTRKM,'KEYMRG      ',KEYMRG)
      CALL LCMGET(IPTRKM,'MATALB      ',MATALB)
      CALL LCMGET(IPTRKM,'VOLSUR      ',VOLSUR)
      CALL LCMSIX(IPTRKM,'EXCELL      ',2)
      ALLOCATE(MATMRG(NUNK),VOLMRG(NUNK))
      CALL XELCMP(NS    ,NV    ,
     >            VOLSUR,MATALB,KEYMRG,
     >            NSOUT ,NREG  ,VOLMRG,MATMRG,
     >            ITGEO ,ICODE )
      MAXMIX=0
      DO 100 IREG=1,NREG
        KEYMRG(IREG+NSOUT+1)= IREG
        MAXMIX=MAX(MAXMIX,MATMRG(IREG+NSOUT+1))
 100  CONTINUE
      CALL LCMPUT(IPTRKM,'MATCOD',NREG,1,MATMRG(NSOUT+2))
      CALL LCMPUT(IPTRKM,'VOLUME',NREG,2,VOLMRG(NSOUT+2))
      CALL LCMPUT(IPTRKM,'KEYFLX',NREG,1,KEYMRG(NSOUT+2))
      CALL XDRSET(EXTKOP,NSTATE,0.0)
      CALL LCMPUT(IPTRKM,'EXCELTRACKOP',NSTATE,2,EXTKOP)
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NREG
      ISTATE(2)=NREG
      ISTATE(4)=MAXMIX
      ISTATE(5)=NSOUT
      ISTATE(7)=ITROP
      ISTATE(8)=-1
      CALL LCMPUT(IPTRKM,'STATE-VECTOR',NSTATE,1,ISTATE)
      DEALLOCATE(VOLMRG,MATMRG,VOLSUR,MATALB,KEYMRG)
*----
*  IF IPRINT >= 20
*  EXPORT TEMPORARY TRACKING FILE
*----
      IF(IPRINT .GE. 10) THEN 
        IUEXP=KDROPN('AXGGEOEXPTRK',0,3,0,0)
        CALL LCMEXP(IPTRKM,IPRINT,IUEXP,2,1)
        IRC=KDRCLS(IUEXP,1)
      ENDIF
      RETURN
      END
