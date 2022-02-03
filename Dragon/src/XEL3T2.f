*DECK XEL3T2
      SUBROUTINE XEL3T2(IX,IY,IZ,LDIM,N3MS,N3MR,N3RS,LMESH,NZP,N2MS,
     1                  N2MR,N2RS,N3S,N3R,NFI,MINDIM,MAXDIM,REMESH,
     2                  VOLSUR,MATALB,KEYMRG,INDEX,MAX2,MIN2,ICOR2,
     3                  REM2,VOL2,MAT2,KEY2,IND2,IND2T3,MATMRG,VOLMRG,
     4                  ZCOR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Create 2D projection (EXCELT geometry analysis) of a 3D prismatic
* geometry.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
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
* LDIM    dimension of MINDIM,MAXDIM,MAX2,MIN2 arrays.
* N3MS    maximum number of outer surfaces for the 3D geometry.
* N3MR    maximum number of regions for the 3D geometry.
* N3RS    second dimension of INDEX array.
* LMESH   dimension of REMESH array.
* NZP     number of plan in the 3D prismatic geometry.
* N2MS    maximum number of outer surfaces for the 2D projected
*         geometry.
* N2MR    maximum number of regions for the 2D projected geometry.
* N2RS    second dimension of IND2 array.
* MINDIM  min index values for the different axes of the 3D geometry.
* MAXDIM  max index values for the different axes of the 3D geometry.
* REMESH  different meshes of the 3D geometry.
* VOLSUR  outer surfaces and volumes for the 3D geometry.
* MATALB  albedo and material indexes for the 3D geometry.
* KEYMRG  merging index for the 3D geometry.
* INDEX   location index for the 3D geometry.
*
*Parameters: output
* N3S     number of outer surfaces for the 3D geometry after merging.
* N3R     number of regions for the 3D geometry after merging.
* NFI     effective dimension for MATMRG AND VOLMRG arrays. 
* MAX2    min index values for the different axes of the 2D projected
*         geometry.
* MIN2    max index values for the different axes of the 2D projected
*         geometry.
* ICOR2   undefined.
* REM2    different meshes of the 2D projected geometry.
* VOL2    outer surfaces and volumes for the 2D projected geometry.
* MAT2    albedo and material indexes for the 2D projected geometry.
* KEY2    merging index for the 2D projected geometry.
* IND2    location index for the 2D projected geometry.
* IND2T3  mapping index between the 2D projected geometries (plan by
*         plan) and the initial 3D geometry.
* MATMRG  albedo and material indexes for the 3D geometry after merging.
* VOLMRG  outer surfaces and volumes for the 3D geometry after merging.
* ZCOR    coordinates of the different plans of the 3D prismatic
*         geometry.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IX,IY,IZ,LDIM,N3MS,N3MR,N3RS,LMESH,NZP,N2MS,N2MR,N2RS,
     1 N3S,N3R,NFI,MINDIM(LDIM),MAXDIM(LDIM),MATALB(-N3MS:N3MR),
     2 KEYMRG(-N3MS:N3MR),INDEX(4,N3RS),MAX2(LDIM),MIN2(LDIM),
     3 ICOR2(LDIM),MAT2(-N2MS:N2MR),KEY2(-N2MS:N2MR),IND2(4,N2MR),
     4 IND2T3(-N2MS:N2MR,0:NZP+1),MATMRG(N3RS)
      REAL REMESH(LMESH),VOLSUR(-N3MS:N3MR),REM2(LMESH),
     1 VOL2(-N2MS:N2MR),VOLMRG(N3RS),ZCOR(0:NZP)
*----
*  LOCAL VARIABLES
*----
      INTEGER II,IR,IS,ZPL,IMRG,IPOS,ITYP,ZPLB,IPOSB,IIB,IPOS2,IMRG2
      REAL DELZ
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITEMP,ITEMP2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ITEMP(LMESH),ITEMP2(LMESH))
*---
* CREATE MATMRG AND VOLMRG ARRAYS FOR 3D GEOMETRY
*---
      CALL XDRSET(VOLMRG,N3RS,0.0)
      N3S=0
      DO IS=-N3MS,-1
         N3S=MIN(KEYMRG(IS),N3S)
      ENDDO
      N3S=-N3S
      DO IS=-N3MS,-1
         VOLMRG(N3S+1+KEYMRG(IS))=VOLMRG(N3S+1+KEYMRG(IS))+VOLSUR(IS)
         MATMRG(N3S+1+KEYMRG(IS))=MATALB(IS)
      ENDDO
      N3R=0
      MATMRG(N3S+1)=0
      DO IR=1,N3MR
         N3R=MAX(KEYMRG(IR),N3R)
         VOLMRG(N3S+1+KEYMRG(IR))=VOLMRG(N3S+1+KEYMRG(IR))+VOLSUR(IR)
         MATMRG(N3S+1+KEYMRG(IR))=MATALB(IR)
      ENDDO
      NFI=N3S+N3R+1
*---
* CREATE PERMUTATION ARRAY FOR THE CHANGE OF COORDINATES (IX,IY,IZ) ->
* (1,2,3)
*---
      IPOS=0
      DO II=MINDIM(IX),MAXDIM(IX)
         IPOS=IPOS+1
         ITEMP(IPOS)=II
      ENDDO
      DO II=MINDIM(IY),MAXDIM(IY)
         IPOS=IPOS+1
         ITEMP(IPOS)=II
      ENDDO
      DO II=MINDIM(IZ),MAXDIM(IZ)
         IPOS=IPOS+1
         ITEMP(IPOS)=II
      ENDDO
      DO II=4,LDIM
         DO IPOSB=MINDIM(II)-2,MAXDIM(II)
            IPOS=IPOS+1
            ITEMP(IPOS)=IPOSB
         ENDDO
      ENDDO
      DO II=1,LMESH
         ITEMP2(ITEMP(II))=II
      ENDDO
*---
* CREATE MAXDIM, MINDIM, REMESH, FOR 2D GEOMETRY
*---
      MIN2(1)=1
      MAX2(1)=MAXDIM(IX)-MINDIM(IX)+MIN2(1)
      ICOR2(1)=1
      DO II=MIN2(1),MAX2(1)
         REM2(II)=REMESH(ITEMP(II))
      ENDDO
      MIN2(2)=MAX2(1)+1
      MAX2(2)=MAXDIM(IY)-MINDIM(IY)+MIN2(2)
      ICOR2(2)=2
      DO II=MIN2(2),MAX2(2)
         REM2(II)=REMESH(ITEMP(II))
      ENDDO
      MIN2(3)=MAX2(2)+1
      MAX2(3)=MAXDIM(IZ)-MINDIM(IZ)+MIN2(3)
      ICOR2(3)=3
      DO II=MIN2(3),MAX2(3)
         REM2(II)=REMESH(ITEMP(II))
      ENDDO
      DO II=4,LDIM
         MIN2(II)=MINDIM(II)
         MAX2(II)=MAXDIM(II)
         ICOR2(II)=3
         DO IPOS=MIN2(II)-2,MAX2(II)
            REM2(IPOS)=REMESH(ITEMP(IPOS))
         ENDDO
      ENDDO
*---
* CREATE Z-COORDINATES ARRAY
*---
      DO II=MINDIM(IZ),MAXDIM(IZ)
         ZCOR(II-MINDIM(IZ))=REMESH(II)-REMESH(MINDIM(IZ))
      ENDDO
*---
* CREATE INDEX FOR 2D GEOMETRY (FROM FIRST Z-PLAN) AND MAPPING INDEX
* BETWEEN 2D AND 3D
*---      
      CALL XDISET(IND2T3(-N2MS,0),N2RS*(NZP+2),0)
*
      IR=0
      IS=-N3MS-1
      DO 10 II=1,N3RS
         IF ((INDEX(1,II).EQ.0).AND.
     1       (INDEX(2,II).EQ.0).AND.
     2       (INDEX(3,II).EQ.0)) GOTO 10
         ZPL=INDEX(IZ,II)-MINDIM(IZ)+1
         IPOS=0
*        what is the element we have encountered?
*        find if this (ix,iy,it) INDEX position has already been
*        encountered in another iz-plan
         IF ((ZPL.EQ.0).OR.(ZPL.EQ.(NZP+1))) THEN
            ITYP=-1
            IF (ZPL.EQ.0) THEN
*           It is a bottom surface
               ZPLB=NZP+1
            ELSE
*           It is a top surface
               ZPLB=0
            ENDIF
*           scan for a matching top/bottom surface
            DO IPOSB=1,N2MR
            IF (IND2T3(IPOSB,ZPLB).NE.0) THEN
               IIB=N3MS+1+IND2T3(IPOSB,ZPLB)
               IF ((INDEX(IX,IIB).EQ.INDEX(IX,II)).AND.
     1             (INDEX(IY,IIB).EQ.INDEX(IY,II)).AND.
     1             (INDEX(4,IIB).EQ.INDEX(4,II))) THEN
                  IPOS=IPOSB
                  GOTO 31
               ENDIF
            ENDIF
            ENDDO
*           scan for a matching region
            DO 21 ZPLB=1,NZP
            DO 20 IPOSB=1,N2MR
            IF (IND2T3(IPOSB,ZPLB).NE.0) THEN
               IIB=N3MS+1+IND2T3(IPOSB,ZPLB)
               IF ((INDEX(IX,IIB).EQ.INDEX(IX,II)).AND.
     1             (INDEX(IY,IIB).EQ.INDEX(IY,II)).AND.
     1             (INDEX(4,IIB).EQ.INDEX(4,II))) THEN
                  IPOS=IPOSB
                  GOTO 31
               ENDIF
            ENDIF
 20         CONTINUE
 21         CONTINUE
*           find an empty space
            DO 41 IPOSB=1,N2MR
               DO ZPLB=0,NZP+1
                  IF (IND2T3(IPOSB,ZPLB).NE.0) GOTO 41
               ENDDO
               IPOS=IPOSB
               GOTO 31
 41         CONTINUE
            CALL XABORT('XEL3T2: INCOMPATIBLE MESHES(1).')
 31         CONTINUE
         ELSE
            IF ((INDEX(IX,II).LT.MINDIM(IX)).OR.
     1          (INDEX(IY,II).LT.MINDIM(IY)).OR.
     2          (INDEX(IX,II).EQ.MAXDIM(IX)).OR.
     3          (INDEX(IY,II).EQ.MAXDIM(IY))) THEN
*           It is a lateral surface
               ITYP=-1
*              scan for a matching lateral surface
               DO 23 ZPLB=1,NZP
               DO 22 IPOSB=-N2MS,-1
               IF (IND2T3(IPOSB,ZPLB).NE.0) THEN
                  IIB=N3MS+1+IND2T3(IPOSB,ZPLB)
                  IF ((INDEX(IX,IIB).EQ.INDEX(IX,II)).AND.
     1                (INDEX(IY,IIB).EQ.INDEX(IY,II)).AND.
     1                (INDEX(4,IIB).EQ.INDEX(4,II))) THEN
                     IPOS=IPOSB
                     GOTO 32
                  ENDIF
               ENDIF
 22            CONTINUE
 23            CONTINUE
*              find an empty space
               DO 42 IPOSB=-N2MS,1
                  DO ZPLB=1,NZP
                     IF (IND2T3(IPOSB,ZPLB).NE.0) GOTO 42
                  ENDDO
                  IPOS=IPOSB
                  GOTO 32
 42            CONTINUE
               CALL XABORT('XEL3T2: INCOMPATIBLE MESHES(2).')
 32            CONTINUE
            ELSE
*           It is a region
               ITYP=1
*              scan for a matching top or bottom surface
               DO 25 ZPLB=0,NZP+1,NZP+1
               DO 24 IPOSB=1,N2MR
               IF (IND2T3(IPOSB,ZPLB).NE.0) THEN
                  IIB=N3MS+1+IND2T3(IPOSB,ZPLB)
                  IF ((INDEX(IX,IIB).EQ.INDEX(IX,II)).AND.
     1                (INDEX(IY,IIB).EQ.INDEX(IY,II)).AND.
     1                (INDEX(4,IIB).EQ.INDEX(4,II))) THEN
                     IPOS=IPOSB
                     GOTO 33
                  ENDIF
               ENDIF
 24            CONTINUE
 25            CONTINUE
*              scan for a matching region
               DO 27 ZPLB=1,NZP
               DO 26 IPOSB=1,N2MR
               IF (IND2T3(IPOSB,ZPLB).NE.0) THEN
                  IIB=IND2T3(IPOSB,ZPLB)+N3MS+1
                  IF ((INDEX(IX,IIB).EQ.INDEX(IX,II)).AND.
     1                (INDEX(IY,IIB).EQ.INDEX(IY,II)).AND.
     1                (INDEX(4,IIB).EQ.INDEX(4,II))) THEN
                     IPOS=IPOSB
                     GOTO 33
                  ENDIF
               ENDIF
 26            CONTINUE
 27            CONTINUE
*              find an empty space
               DO 43 IPOSB=1,N2MR
                  DO ZPLB=0,NZP+1
                     IF (IND2T3(IPOSB,ZPLB).NE.0) GOTO 43
                  ENDDO
                  IPOS=IPOSB
                  GOTO 33
 43            CONTINUE
               CALL XABORT('XEL3T2: INCOMPATIBLE MESHES(3).')
 33            CONTINUE
            ENDIF
         ENDIF
         IF (ITYP.EQ.-1) THEN
            IS=IS+1
            IMRG=IS
         ELSE
            IR=IR+1
            IMRG=IR
         ENDIF
         IND2T3(IPOS,ZPL)=IMRG
         DO ZPLB=0,NZP+1
            IF (IND2T3(IPOS,ZPLB).NE.0) THEN
               IIB=N3MS+1+IND2T3(IPOS,ZPLB)
               IF ((INDEX(IX,IIB).NE.INDEX(IX,II)).OR.
     1             (INDEX(IY,IIB).NE.INDEX(IY,II)).OR.
     1             (INDEX(4,IIB).NE.INDEX(4,II))) THEN
                  WRITE(6,*) ZPLB,IND2T3(IPOS,ZPLB),IIB
                  WRITE(6,*) INDEX(IX,IIB),INDEX(IX,II)
                  WRITE(6,*) INDEX(IY,IIB),INDEX(IY,II)
                  WRITE(6,*) INDEX(4,IIB),INDEX(4,II)
                  CALL XABORT('XEL3T2: PROJECTION FAILED (1).')
               ENDIF
            ENDIF
         ENDDO
 10   CONTINUE
      IF ((IR.NE.N3MR).OR.(-IS.NE.1)) THEN
         WRITE(6,*) N3MR,IR
         WRITE(6,*) 1,-IS
         CALL XABORT('XEL3T2: PROJECTION FAILED (2).')
      ENDIF
*---
* CREATE VOLSUR, MATALB, KEYMRG INDEX FOR 2D GEOMETRY (FROM FIRST
* Z-PLAN) APPLY KEYMRG TO MAPPING INDEX BETWEEN 2D AND 3D
*---
      DELZ=REMESH(MINDIM(IZ)+1)-REMESH(MINDIM(IZ))
      DO 60 IPOS=-N2MS,N2MR
         IMRG=IND2T3(IPOS,1)
         VOL2(IPOS)=VOLSUR(IMRG)/DELZ
         MAT2(IPOS)=MATALB(IMRG)
         KEY2(IPOS)=IPOS
         IPOS2=IPOS+N2MS+1
         IMRG2=IMRG+N3MS+1
         IF (INDEX(IX,IMRG2).LT.MINDIM(IX)) THEN
            IND2(1,IPOS2)=MIN2(1)-1
         ELSE
            IND2(1,IPOS2)=ITEMP2(INDEX(IX,IMRG2))
         ENDIF
         IF (INDEX(IY,IMRG2).LT.MINDIM(IY)) THEN
            IND2(2,IPOS2)=MIN2(2)-1
         ELSE
            IND2(2,IPOS2)=ITEMP2(INDEX(IY,IMRG2))
         ENDIF
         IND2(3,IPOS2)=ITEMP2(MINDIM(IZ))
         IF (INDEX(4,IMRG2).EQ.0) THEN
            IND2(4,IPOS2)=0
         ELSE
            IND2(4,IPOS2)=ITEMP2(INDEX(4,IMRG2))
         ENDIF
         DO 50 ZPL=0,NZP+1
            IMRG=IND2T3(IPOS,ZPL)
            IND2T3(IPOS,ZPL)=KEYMRG(IMRG)
 50      CONTINUE
 60   CONTINUE
      IND2(3,N2MS+1)=0
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ITEMP2,ITEMP)
      RETURN
      END
