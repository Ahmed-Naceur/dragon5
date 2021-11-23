*DECK NXTIND
      SUBROUTINE NXTIND(IX,IY,IZ,NFSUR,NFREG,MXGSUR,MXGREG,MAXMSH,NZP,
     1                  NUCELZ,MESHCZM,MESHC,NSURC,NREGC,INDEX,IDREG,
     2                  IDSUR,N2REG,N2SUR,IND2T3,REGI,NZC,IDZ,LSTORE,
     3                  I2SURC,N2REGC,II,JJ,LL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Locate the different regions corresponding to the same projection along
* an axis for a set of cells/pins.
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
* IX      first direction perpendicular to the projection axis.
* IY      second direction perpendicular to the projection axis.
* IZ      projection axis.
* NFSUR   number of outer surfaces in the 3D geometry.
* NFREG   number of regions in the 3D geometry.
* MXGSUR  maximum number of surfaces for any sub-geometry of 
*         the 3D geometry.
* MXGREG  maximum number of regions for any sub-geometry of 
*         the 3D geometry.
* MAXMSH  maximum dimension of any mesh in any sub-geometry of 
*         the 3D geometry.
* NZP     total number of plans in the 3D geometry.
* NUCELZ  number of cells/pins along the projection axis.
* MESHCZM maximum number of meshes along the projection axis within 
*         any cell/pin.
* MESHC   cells/pins meshes size.
* NSURC   number of surfaces for the cells/pins.
* NREGC   number of regions for the cells/pins.
* LSTORE  2D cell/pin storage flag.
* II      x index to locate.
* JJ      y index to locate.
* LL      r index to locate.
*
*Parameters: input/output
* INDEX   cells/pins index vector for 3D cells/pins and corresponding
*         2D cell/pin.
* IDSUR   surface index array for 3D cells/pins and corresponding 2D 
*         cell/pin.
* IDREG   region index array for 3D cells/pins and corresponding 2D 
*         cell/pin.
* N2REG   number of regions in the projected 2D geometry.
* N2SUR   number of outer surfaces in the projected 2D geometry.
* IND2T3  mapping index between the 2D projected geometries 
*         (plan by plan) and the initial 3D geometry.
* REGI    region sweeping flag array.
* NZC     array containing the number of meshes alon the projection 
*         axis for each cell/pin.
* I2SURC  initial/final outer surface position in surface index array 
*         for corresponding 2D cell/pin.
* N2REGC  initial/final outer surface position in region index array 
*         for corresponding 2D cell/pin.
*
*Parameters: temporary storage
* IDZ     work vector.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IX,IY,IZ,NFSUR,NFREG,MXGSUR,MXGREG,MAXMSH,NZP,
     1 NUCELZ,MESHCZM,MESHC(4,NUCELZ),NSURC(NUCELZ),NREGC(NUCELZ),
     2 INDEX(5,-MXGSUR:MXGREG,0:NUCELZ),IDSUR(MXGSUR,0:NUCELZ),
     3 IDREG(MXGREG,0:NUCELZ),N2REG,N2SUR,
     4 IND2T3(-NFSUR:NFREG,0:NUCELZ*MAXMSH+1),REGI(-NFSUR:NFREG),
     5 NZC(NUCELZ),IDZ(0:MESHCZM+1,NUCELZ),I2SURC,N2REGC,II,JJ,LL
      LOGICAL LSTORE
*----
*  LOCAL VARIABLES
*----
      INTEGER K,KK,MESHCZ,ITYP,ISUR,IREG
*
      DO K=1,NUCELZ
         MESHCZ=MESHC(IZ,K)
         CALL XDISET(IDZ(0,K),MESHCZ+2,0)
         CALL NXTFID(IX,IY,IZ,NSURC(K),NREGC(K),MESHCZ,II,JJ,LL,
     1        INDEX(1,-MXGSUR,K),IDSUR(1,K),IDREG(1,K),IDZ(0,K),ITYP)
      ENDDO
      IF (ITYP.EQ.-1) THEN
*     lateral surface found
         IF (REGI(IDZ(1,1)).EQ.0) THEN
            N2SUR=N2SUR-1
            DO 20 K=1,NUCELZ
            MESHCZ=MESHC(IZ,K)
            DO 10 KK=1,MESHCZ
               ISUR=IDZ(KK,K)
               IND2T3(N2SUR,NZC(K)+KK)=ISUR
               REGI(ISUR)=N2SUR
 10         CONTINUE
 20         CONTINUE
         ENDIF
      ELSEIF (ITYP.EQ.1) THEN
*     region and bottom/top surface found
         IF (REGI(IDZ(1,1)).EQ.0) THEN
            N2REG=N2REG+1
*           region
            DO 40 K=1,NUCELZ
            MESHCZ=MESHC(IZ,K)
            DO 30 KK=1,MESHCZ
               IREG=IDZ(KK,K)
               IND2T3(N2REG,NZC(K)+KK)=IREG
               REGI(IREG)=N2REG
 30         CONTINUE
 40         CONTINUE
*           top surface
            MESHCZ=MESHC(IZ,NUCELZ)
            ISUR=IDZ(MESHCZ+1,NUCELZ)
            IND2T3(N2REG,NZP+1)=ISUR
            REGI(ISUR)=N2REG
*           bottom surface
            ISUR=IDZ(0,1)
            IND2T3(N2REG,0)=ISUR   
            REGI(ISUR)=N2REG
         ENDIF           
      ENDIF
      IF (LSTORE) THEN !.AND.(ITYP.NE.0)) THEN
*     STORE THE CORRESPONDING 2D CELL CONTENTS
         IF (ITYP.LT.0) THEN
*        a surface
            I2SURC=I2SURC+1
            IF (ITYP.EQ.-1) THEN
               IDSUR(-I2SURC,0)=ABS(REGI(IDZ(1,1)))
            ELSE
               IDSUR(-I2SURC,0)=0
            ENDIF
            INDEX(1,I2SURC,0)=II
            INDEX(2,I2SURC,0)=JJ
            INDEX(3,I2SURC,0)=1
            INDEX(4,I2SURC,0)=LL
            INDEX(5,I2SURC,0)=0
         ELSE
*        a region
            N2REGC=N2REGC+1
            IF (ITYP.EQ.1) THEN
               IDREG(N2REGC,0)=ABS(REGI(IDZ(1,1)))
            ELSE
               IDREG(N2REGC,0)=0
            ENDIF
            INDEX(1,N2REGC,0)=II
            INDEX(2,N2REGC,0)=JJ
            INDEX(3,N2REGC,0)=1
            INDEX(4,N2REGC,0)=LL
            INDEX(5,N2REGC,0)=0
         ENDIF
      ENDIF
*
      RETURN
      END
