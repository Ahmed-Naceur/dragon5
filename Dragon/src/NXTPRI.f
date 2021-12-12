*DECK NXTPRI
      SUBROUTINE NXTPRI(IPTRK,JPTRK,IX,IY,IZ,NFREG,NFSUR,MAXMSH,NUCELL,
     1                  MXGSUR,MXGREG,INDEX,IDSUR,IDREG,MESHC,NSURC,
     2                  NREGC,IDIRC,NZP,N2REG,N2SUR,IND2T3,REGI,
     3                  DEPS,DCMESH,ZCORD,LFIRST,LSTORE,ILEV,IEL,N2EL,
     4                  N2SURC,N2REGC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Check compatibility of cells/pins along the projection axis for a 
* 3D prismatic treatment and generate corresponding 2D cell/pin.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): R. Le Tellier
*
*Parameters: input
* IPTRK   pointer to the NXT 3D geometry analysis.
* JPTRK   pointer to the NXT 2D projected geometry analysis.
* IX      first direction perpendicular to the projection axis.
* IY      second direction perpendicular to the projection axis.
* IZ      projection axis.
* NFREG   number of regions in the 3D geometry.
* NFSUR   number of outer surfaces in the 3D geometry.
* MAXMSH  maximum dimension of any mesh in any sub-geometry of the 3D
*         geometry.
* NUCELL  number of cells along the three axis in the 3D geometry.
* MXGSUR  maximum number of surfaces for any sub-geometry of the 3D
*         geometry.
* MXGREG  maximum number of regions for any sub-geometry of the 3D
*         geometry.
* INDEX   cells/pins index vector.
* IDSUR   surface index array.
* IDREG   region index array.
* MESHC   cells/pins meshes size.
* NSURC   number of surfaces for the cells/pins.
* NREGC   number of regions for the cells/pins.
* IDIRC   cylinders orientations.
* DEPS    comparison criterion for double precision values.
* DCMESH  cells/pins meshing vector.
* LFIRST  first cell/pin processed flag.
* LSTORE  2D cell/pin storage flag.
* ILEV    geometry level 1:cells / 2:pins.
* IEL     index of the first 3D cell/pin considered.
* N2EL    index of the corresponding 2D cell.
* 
*Parameters: input/output
* NZP     number of plans in the 3D prismatic geometry.
* N2REG   number of regions in the projected 2D geometry.
* N2SUR   number of outer surfaces in the projected 2D geometry.
* IND2T3  mapping index between the 2D projected geometries (plan by
*         plan) and the initial 3D geometry.
* REGI    region sweeping flag array.
* ZCORD   coordinates of the different plans of the 3D prismatic
*         geometry.
* N2SURC  number of outer surfaces in the 2D corresponding cell/pin.
* N2REGC  number of outer regions in the 2D corresponding cell/pin.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,JPTRK
      INTEGER IX,IY,IZ,NFREG,NFSUR,MAXMSH,NUCELL(3),MXGSUR,
     1 MXGREG,INDEX(5,-MXGSUR:MXGREG,0:NUCELL(IZ)),
     2 IDSUR(MXGSUR,0:NUCELL(IZ)),IDREG(MXGREG,0:NUCELL(IZ)),
     3 MESHC(4,NUCELL(IZ)),NSURC(NUCELL(IZ)),NREGC(NUCELL(IZ)),
     4 IDIRC(NUCELL(IZ)),NZP,N2REG,N2SUR,
     7 IND2T3(-NFSUR:NFREG,0:NUCELL(IZ)*MAXMSH+1),REGI(-NFSUR:NFREG),
     8 ILEV,IEL,N2EL,N2SURC,N2REGC
      DOUBLE PRECISION DEPS,DCMESH(-1:MAXMSH,4,0:NUCELL(IZ)),
     1 ZCORD(0:MAXMSH)
      LOGICAL LFIRST,LSTORE
*----
*  LOCAL VARIABLES
*----
      INTEGER K,JJ,MESHR(4),MESHCZM,I2SURC,II,LL
      DOUBLE PRECISION DELZ
      LOGICAL XDDCOM
      CHARACTER NAMCEL*9,NAMREC*12,NAMCE2*9
      CHARACTER CDIR(4)*1,CLEV(2)*1
      DATA CDIR /'X','Y','Z','R'/
      DATA CLEV /'C','P'/
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NZC,IDZ
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SM
*----
*  Scratch storage allocation
*----
      ALLOCATE(NZC(NUCELL(IZ)))
*----
*  TEST THE CELL/PIN DIMENSIONS COMPATIBILITY
*----
      IF ((IDIRC(1).NE.IZ).AND.(IDIRC(1).NE.0)) 
     1 CALL XABORT('NXTPRI: INVALID PRISMATIC GEOMETRY (IDIRC).')
      MESHR(1)=MESHC(IX,1)
      MESHR(2)=MESHC(IY,1)
      MESHR(4)=MESHC(4,1)
      NZC(1)=0
      MESHCZM=MESHC(IZ,1)
      DO K=2,NUCELL(IZ)
         IF ((IDIRC(K).NE.IZ).AND.(IDIRC(K).NE.0)) THEN
           WRITE(6,*) IDIRC(K),IZ
           CALL XABORT('NXTPRI: INVALID PRISMATIC GEOMETRY (IDIRC).')
         ENDIF
         IF (MESHC(IX,K).NE.MESHR(1)) THEN
           WRITE(6,*)  MESHC(IX,K),MESHR(1)
           CALL XABORT('NXTPRI: INVALID PRISMATIC GEOMETRY (MESHC(1)).')
         ENDIF
         IF (MESHC(IY,K).NE.MESHR(2)) THEN
           WRITE(6,*)  MESHC(IY,K),MESHR(2)
           CALL XABORT('NXTPRI: INVALID PRISMATIC GEOMETRY (MESHC(2)).')
         ENDIF
         IF (MESHC(4,K).NE.MESHR(4)) THEN
           WRITE(6,*)  MESHC(4,K),MESHR(4)
           CALL XABORT('NXTPRI: INVALID PRISMATIC GEOMETRY (MESHC(4)).')
         ENDIF
         NZC(K)=NZC(K-1)+MESHC(IZ,K)
         MESHCZM=MAX(MESHCZM,MESHC(IZ,K))
      ENDDO
*----
*  TEST THE CELL/PIN MESHES COMPATIBILITY
*----
      DO K=2,NUCELL(IZ)
         DO JJ=-1,MESHR(1)
            IF (.NOT.XDDCOM(DCMESH(JJ,IX,K),DCMESH(JJ,IX,1),DEPS)) THEN
             WRITE(6,*) IX,DCMESH(JJ,IX,K),DCMESH(JJ,IX,1)
             CALL XABORT('NXTPRI: INVALID PRISMATIC GEOMETRY (DCMESH).')
            ENDIF
         ENDDO
         DO JJ=-1,MESHR(2)
            IF (.NOT.XDDCOM(DCMESH(JJ,IY,K),DCMESH(JJ,IY,1),DEPS)) THEN
             WRITE(6,*) IY,DCMESH(JJ,IY,K),DCMESH(JJ,IY,1)
             CALL XABORT('NXTPRI: INVALID PRISMATIC GEOMETRY (DCMESH).')
          ENDIF
         ENDDO
         IF (MESHR(4).GT.0) THEN
         DO JJ=-1,MESHR(4)
            IF (.NOT.XDDCOM(DCMESH(JJ,4,K),DCMESH(JJ,4,1),DEPS)) THEN
             WRITE(6,*) 4,DCMESH(JJ,4,K),DCMESH(JJ,4,1)
             CALL XABORT('NXTPRI: INVALID PRISMATIC GEOMETRY (DCMESH).')
            ENDIF
         ENDDO
         ENDIF
      ENDDO
      IF (LFIRST) THEN
*----
*  CALCULATE GLOBAL MESH ALONG THE PROJECTION AXIS
*----
         NZP=0
         ZCORD(0)=0.0
         DO 20 K=1,NUCELL(IZ)
         DO 10 JJ=1,MESHC(IZ,K)
            NZP=NZP+1
            ZCORD(NZP)=ZCORD(NZP-1)+(DCMESH(JJ,IZ,K)-DCMESH(JJ-1,IZ,K))
 10      CONTINUE
 20      CONTINUE
      ELSE
*----
*  TEST THE COMPATIBILITY OF THE GLOBAL MESH ALONG THE PROJECTION AXIS
*----
         NZP=0
         DELZ=0.D0
         DO 40 K=1,NUCELL(IZ)
         DO 30 JJ=1,MESHC(IZ,K)
            NZP=NZP+1
            DELZ=DELZ+(DCMESH(JJ,IZ,K)-DCMESH(JJ-1,IZ,K))
            IF (.NOT.XDDCOM(DELZ,ZCORD(NZP),DEPS))
     1       CALL XABORT('NXTPRI: INVALID PRISMATIC GEOMETRY (ZCORD).')
            
 30      CONTINUE 
 40      CONTINUE 
      ENDIF
      IF (LSTORE) THEN
*----
*  STORE THE CORRESPONDING 2D CELL/PIN CONTENTS
*----
         WRITE(NAMCEL,'(A1,I8.8)') CLEV(ILEV),IEL
         WRITE(NAMCE2,'(A1,I8.8)') CLEV(ILEV),N2EL
         !write(*,*) 'storing ',NAMCE2,MESHR(1),MESHR(2),MESHR(4)
         ALLOCATE(SM(MESHR(1)+2))
           NAMREC=NAMCEL//'SM'//CDIR(IX)
           CALL LCMGET(IPTRK,NAMREC,SM)
           NAMREC=NAMCE2//'SM'//CDIR(1)
           CALL LCMPUT(JPTRK,NAMREC,(MESHR(1)+2),4,SM)
         DEALLOCATE(SM)
         ALLOCATE(SM(MESHR(2)+2))
           NAMREC=NAMCEL//'SM'//CDIR(IY)
           CALL LCMGET(IPTRK,NAMREC,SM)
           NAMREC=NAMCE2//'SM'//CDIR(2)
           CALL LCMPUT(JPTRK,NAMREC,(MESHR(2)+2),4,SM)
         DEALLOCATE(SM)
         IF (MESHR(4).GT.0) THEN
            ALLOCATE(SM(MESHR(4)+2))
              NAMREC=NAMCEL//'SM'//CDIR(4)
              CALL LCMGET(IPTRK,NAMREC,SM)
              NAMREC=NAMCE2//'SM'//CDIR(4)
              CALL LCMPUT(JPTRK,NAMREC,(MESHR(4)+2),4,SM)
            DEALLOCATE(SM)
         ENDIF
      ENDIF
*----
*  UPDATE IND2T3 FOR THIS SET OF CELLS/PINS
*----
      ALLOCATE(IDZ(NUCELL(IZ)*(MESHCZM+3)))
      N2SURC=(NSURC(1)-2*NREGC(1)/MESHC(IZ,1))/MESHC(IZ,1)
      I2SURC=-N2SURC-1
      N2REGC=0
      IF (ILEV.EQ.2) THEN
*     R+ surface     
         CALL NXTIND(IX,IY,IZ,NFSUR,NFREG,MXGSUR,MXGREG,MAXMSH,NZP,
     1        NUCELL(IZ),MESHCZM,MESHC,NSURC,NREGC,INDEX,IDREG,IDSUR,
     2        N2REG,N2SUR,IND2T3,REGI,NZC,IDZ,LSTORE,I2SURC,
     3        N2REGC,0,0,-2)
      ELSE
*     Y- Y+ surfaces
      DO 60 JJ=-2,-1
      DO 50 II=MESHC(IX,1),1,-1
         LL=0
         CALL NXTIND(IX,IY,IZ,NFSUR,NFREG,MXGSUR,MXGREG,MAXMSH,NZP,
     1        NUCELL(IZ),MESHCZM,MESHC,NSURC,NREGC,INDEX,IDREG,IDSUR,
     2        N2REG,N2SUR,IND2T3,REGI,NZC,IDZ,LSTORE,I2SURC,
     3        N2REGC,II,JJ,LL)
 50   CONTINUE
 60   CONTINUE
*     X- X+ surfaces
      DO 80 II=-2,-1
      DO 70 JJ=MESHC(IY,1),1,-1
         LL=0
         CALL NXTIND(IX,IY,IZ,NFSUR,NFREG,MXGSUR,MXGREG,MAXMSH,NZP,
     1        NUCELL(IZ),MESHCZM,MESHC,NSURC,NREGC,INDEX,IDREG,IDSUR,
     2        N2REG,N2SUR,IND2T3,REGI,NZC,IDZ,LSTORE,I2SURC,
     3        N2REGC,II,JJ,LL)
 70   CONTINUE
 80   CONTINUE
      ENDIF
*     regions
      DO 100 JJ=1,MESHC(IY,1)
      DO 90 II=1,MESHC(IX,1)
      DO LL=1,MESHC(4,1)
         CALL NXTIND(IX,IY,IZ,NFSUR,NFREG,MXGSUR,MXGREG,MAXMSH,NZP,
     1        NUCELL(IZ),MESHCZM,MESHC,NSURC,NREGC,INDEX,IDREG,IDSUR,
     2        N2REG,N2SUR,IND2T3,REGI,NZC,IDZ,LSTORE,I2SURC,
     3        N2REGC,II,JJ,LL)
      ENDDO
      IF (ILEV.EQ.1) THEN
         LL=0
         CALL NXTIND(IX,IY,IZ,NFSUR,NFREG,MXGSUR,MXGREG,MAXMSH,NZP,
     1        NUCELL(IZ),MESHCZM,MESHC,NSURC,NREGC,INDEX,IDREG,IDSUR,
     2        N2REG,N2SUR,IND2T3,REGI,NZC,IDZ,LSTORE,I2SURC,
     3        N2REGC,II,JJ,LL)  
      ENDIF
 90   CONTINUE
100   CONTINUE
      DEALLOCATE(IDZ)
      IF (LSTORE) THEN
*----
*  STORE THE CORRESPONDING 2D CELL/PIN CONTENTS
*----
         IF ((I2SURC.NE.-1).OR.(N2REGC.NE.NREGC(1)/MESHC(IZ,1))) THEN
            WRITE(6,*) I2SURC,-1,' OR', N2REGC,NREGC(1)/MESHC(IZ,1)
            CALL XABORT('NXTPRI: INVALID NUMBER OF SURFACES/REGIONS')
         ENDIF
         WRITE(NAMCEL,'(A1,I8.8)') CLEV(ILEV),N2EL
         NAMREC=NAMCEL//'VSI'
         INDEX(1,0,0)=0
         INDEX(2,0,0)=0
         INDEX(3,0,0)=0
         INDEX(4,0,0)=0         
         INDEX(5,0,0)=0
         CALL LCMPUT(JPTRK,NAMREC,5*(N2SURC+N2REGC+1),1,
     1                            INDEX(1,-N2SURC,0))
         NAMREC=NAMCEL//'RID'
         CALL LCMPUT(JPTRK,NAMREC,N2REGC,1,IDREG(1,0))
         NAMREC=NAMCEL//'SID'
         CALL LCMPUT(JPTRK,NAMREC,N2SURC,1,IDSUR(1,0))
      ENDIF
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(NZC)
      RETURN
      END
