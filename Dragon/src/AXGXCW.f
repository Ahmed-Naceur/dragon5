*DECK AXGXCW
      SUBROUTINE AXGXCW(IPGEOM,IPTRKM,IPRINT,GEONAM,ISYMM )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Analyze XEL geometry WIMS-AECL type tracking with XCWTRK module.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPGEOM  geometry data structure pointer.
* IPTRKM  tracking data structure pointer.         
* IPRINT  print level.                              
* GEONAM  geometry name.
* ISYMM   geometry symmetry.
*
*-----------------------------------------------------------------------
*
      USE          GANLIB
      IMPLICIT     NONE
      INTEGER      IOUT,NALB,MREGIO,NSTATE
      PARAMETER   (IOUT=6,NALB=6,MREGIO=100000,NSTATE=40)
*----
*  ROUTINE PARAMETERS
*----
      TYPE(C_PTR)  IPGEOM,IPTRKM
      INTEGER      IPRINT,ISYMM
      CHARACTER*12 GEONAM
*----
*  INTEGER ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEYMRG,MATALB,NRINFO,NRODS,
     > NRODR,NXRS,NXRI,MATRT 
*----
*  REAL ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLSUR,RAN,RODS,RODR
*----
*  LOCAL VARIABLES
*---- 
      LOGICAL      ILK
      INTEGER      NCODE(NALB),ICODE(NALB)
      REAL         ZCODE(NALB),ALBEDO(NALB) 
      INTEGER      ISTATE(NSTATE)
      INTEGER      NDIM  ,NSUR  ,NVOL  ,MAXJ  ,IROT  ,NBAN  ,
     >             MNAN  ,NRT   ,MSROD ,MAROD ,NSURF ,NSURX ,
     >             NMAT  ,NUNK
      REAL         RADMIN,COTE  
*----
*  SET POSITION VECTOR AND READ ISTATE
*----
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/26H AXGXCW: PROCESS GEOMETRY ,A12)') GEONAM
      ENDIF
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      NDIM=2
      IF(ISTATE(1).EQ.3) THEN
        NSUR=1
      ELSE IF(ISTATE(1).EQ.20) THEN
        NSUR=4
      ELSE IF(ISTATE(1).EQ.24) THEN
        NSUR=6
      ENDIF
      MAXJ=1
      IROT=0
      CALL XCGDIM(IPGEOM,MREGIO,NSUR  ,IROT  ,ISYMM ,MAXJ  ,
     >            NVOL  ,NBAN  ,MNAN  ,NRT   ,MSROD ,MAROD ,
     >            NSURF )
*----
*  CHECK FOR SYMMETRY
*----
      NSURX=NSUR
      IF(ISYMM.GT.1) THEN
        IF(NSURX.EQ.4) THEN
          IROT=-ISYMM-400
        ELSE IF(NSURX.EQ.6) THEN
          IROT=-ISYMM-600
        ELSE
          IROT=-ISYMM-100
        ENDIF
        NSUR=1
      ENDIF
*----
*  ALLOCATE MEMORY FOR PROCESSING GEOMETRY INFORMATION
*----
      ALLOCATE(KEYMRG(NSUR+NVOL+1),MATALB(NSUR+NVOL+1),NRINFO(2*MNAN),
     > NRODS(3*NRT),NRODR(NRT),NXRS(NRT),NXRI(NRT*NBAN))
      ALLOCATE(VOLSUR(NSUR+NVOL+1),RAN(NBAN),RODS(2*NRT),
     > RODR(MSROD*NRT))
*
      CALL XCGGEO(IPGEOM,IROT  ,NSUR  ,NVOL  ,NBAN  ,MNAN  ,
     >            NRT   ,MSROD ,IPRINT,ILK   ,NMAT  ,RAN   ,
     >            NRODS ,RODS  ,NRODR ,RODR  ,NRINFO,MATALB,
     >            VOLSUR,COTE  ,RADMIN,NCODE ,ICODE ,ZCODE ,
     >            ALBEDO,KEYMRG,NXRS  ,NXRI)
*----
*  BUILD BOUNDARY CONDITION MATRIX FOR REFLECTION AND TRANSMISSION
*----
      ALLOCATE(MATRT(NSUR))
      CALL XCGBCM(IPTRKM,NSUR  ,NCODE ,MATRT )
*----
*  SAVE TRACKING FOR CLUSTER GEOMETRY
*----
      CALL XDISET(ISTATE,NSTATE,0)
      NUNK=NVOL+NSUR+1
      ISTATE(1)=NDIM
      ISTATE(2)=NSUR
      ISTATE(3)=NVOL
      ISTATE(4)=NSURX
      ISTATE(5)=NBAN
      ISTATE(6)=NUNK
      ISTATE(7)=NRT
      ISTATE(8)=MSROD
      ISTATE(9)=MAROD
      ISTATE(10)=MNAN
      CALL LCMSIX(IPTRKM,'EXCELL      ',1)
      CALL LCMPUT(IPTRKM,'STATE-VECTOR',NSTATE   ,1,ISTATE)
      CALL LCMPUT(IPTRKM,'RAN         ',NBAN     ,2,RAN   )
      IF(NSURX .EQ. 4)
     >CALL LCMPUT(IPTRKM,'COTE        ',1        ,2,COTE  )
      CALL LCMPUT(IPTRKM,'RADMIN      ',1        ,2,RADMIN)
      CALL LCMPUT(IPTRKM,'NRODS       ',3*NRT    ,1,NRODS )
      CALL LCMPUT(IPTRKM,'RODS        ',2*NRT    ,2,RODS  )
      CALL LCMPUT(IPTRKM,'NRODR       ',NRT      ,1,NRODR )
      CALL LCMPUT(IPTRKM,'RODR        ',MSROD*NRT,2,RODR  )
      CALL LCMPUT(IPTRKM,'NRINFO      ',2*NBAN   ,1,NRINFO)
      CALL LCMPUT(IPTRKM,'NXRI        ',NRT*NBAN ,1,NXRI  )
      CALL LCMPUT(IPTRKM,'NXRS        ',NRT      ,1,NXRS  )
      CALL LCMPUT(IPTRKM,'KEYMRG      ',NUNK     ,1,KEYMRG)
      CALL LCMPUT(IPTRKM,'MATALB      ',NUNK     ,1,MATALB)
      CALL LCMPUT(IPTRKM,'VOLSUR      ',NUNK     ,2,VOLSUR)
      CALL LCMSIX(IPTRKM,'EXCELL      ',2)
      CALL LCMPUT(IPTRKM,'ALBEDO      ',6        ,2,ALBEDO)
      CALL LCMPUT(IPTRKM,'ICODE       ',6        ,1,ICODE )
      CALL LCMPUT(IPTRKM,'NCODE       ',6        ,1,NCODE )
*----
*  RELEASE BLOCKS FOR GEOMETRY
*----
      DEALLOCATE(MATRT)
      DEALLOCATE(RODR,RODS,RAN,VOLSUR)
      DEALLOCATE(NXRI,NXRS,NRODR,NRODS,NRINFO,MATALB,KEYMRG)
      RETURN
      END
