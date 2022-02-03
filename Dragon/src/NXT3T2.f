*DECK NXT3T2
      SUBROUTINE NXT3T2(IPTRK,JPTRK,IX,IY,IZ,NFREG,NFSUR,NUNK,MAXMSH,
     1                  NUCELL,NBUCEL,MXGSUR,MXGREG,MAXPIN,MATALB,
     2                  SURVOL,IUNFLD,NZP,N2REG,N2SUR,N2CEL,N2PIN,
     3                  IND2T3,ZCORD,MATALB2,SURVOL2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Create 2D projection (NXT geometry analysis) of a 3D prismatic
* geometry.
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
* IPTRK   pointer to the NXT 3D geometry analysis.
* JPTRK   pointer to the NXT 2D projected geometry analysis.
* IX      first direction perpendicular to the projection axis.
* IY      second direction perpendicular to the projection axis.
* IZ      projection axis.
* NFREG   number of regions in the 3D geometry.
* NFSUR   number of outer surfaces in the 3D geometry.
* NUNK    NFREG+NFSUR+1.
* MAXMSH  maximum dimension of any mesh in any sub-geometry of the 3D
*         geometry.
* NUCELL  number of cells along the three axis in the 3D geometry.
* NBUCEL  total number of cells in the 3D geometry.
* MXGSUR  maximum number of surfaces for any sub-geometry of the 3D
*         geometry.
* MXGREG  maximum number of regions for any sub-geometry of the 3D
*         geometry.
* MAXPIN  maximum number of pins for any cell of the 3D geometry.
* MATALB  mixtures/albedos array for the 3D geometry.
* SURVOL  surfaces/volumes array for the 3D geometry.
* IUNFLD  assembly description array for the 3D geometry (*,*,*,*,0)
*                                   / projected 2D geometry (*,*,*,*,1).
*
*Parameters: output
* NZP     number of plans in the 3D prismatic geometry.
* N2REG   number of regions in the projected 2D geometry.
* N2SUR   number of outer surfaces in the projected 2D geometry.
* N2CEL   total number of cells in the projected 2D geometry.
* N2PIN   total number of pin descriptions in the projected 2D geometry.
* IND2T3  mapping index between the 2D projected geometries (plan by
*         plan) and the initial 3D geometry.
* ZCORD   coordinates of the different plans of the 3D prismatic
*         geometry.
* MATALB2 mixtures/albedos array for the projected 2D geometry.
* SURVOL2 surfaces/volumes array for the projected 2D geometry.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,JPTRK
      INTEGER IX,IY,IZ,NFREG,NFSUR,NUNK,MAXMSH,NUCELL(3),NBUCEL,MXGSUR,
     1 MXGREG,MAXPIN,MATALB(-NFSUR:NFREG),
     2 IUNFLD(2,NUCELL(1),NUCELL(2),NUCELL(3),0:1),NZP,N2REG,N2SUR,
     3 N2CEL,N2PIN,IND2T3(-NFSUR:NFREG,0:NUCELL(IZ)*MAXMSH+1),
     4 MATALB2(-NFSUR:NFREG)
      DOUBLE PRECISION SURVOL(-NFSUR:NFREG),ZCORD(0:MAXMSH),
     1 SURVOL2(-NFSUR:NFREG)
*----
*  LOCAL VARIABLES
*----
      INTEGER NSTATE
      DOUBLE PRECISION DEPS
      PARAMETER(NSTATE=40,DEPS=1.D-7)
      INTEGER ESTATE(NSTATE)
      INTEGER I,J,ITRN,ICEL,K,JJ,II,NTPINR,N2SURC,N2REGC,IPIN,N2SURP,
     1 N2REGP,NUNK2
      DOUBLE PRECISION DELZ,HPIN,APIN,RPIN,RADP
!!      CHARACTER SIZEX*5,FORM*30
      CHARACTER NAMCEL*9,NAMREC*12,NAMCE2*9
      LOGICAL LFIRST,XDDCOM,LPIN,LSTCEL,LSTPIN
      CHARACTER CDIR(4)*1
      DATA CDIR /'X','Y','Z','R'/
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NSURC,NREGC,IDIRC,NTPIN,
     > REGI,CELID,PINID
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IDSUR,IDREG,MESHC
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: INDEX,ITPIN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DCMESH,DRAPIN
*
!!!!      WRITE(SIZEX,*) NUCELL(1)
!!!!      FORM='('//SIZEX//'(I2,1X,1H(,I2,1H),2X))'
!!!!      WRITE(6,*) 'GLOBAL ASSEMBLY:'
!!!!      DO K=1,MAX(NUCELL(3),1)
!!!!         WRITE(6,*) K,' z-plan'
!!!!         DO J=NUCELL(2),1,-1
!!!!            WRITE(6,FORM) ((IUNFLD(JJ,I,J,K,0),JJ=1,2),I=1,NUCELL(1))
!!!!         ENDDO
!!!!      ENDDO
*----
*  Scratch storage allocation
*----
      ALLOCATE(INDEX(5,-MXGSUR:MXGREG,0:NUCELL(IZ)),
     >         IDSUR(MXGSUR,0:NUCELL(IZ)),IDREG(MXGREG,0:NUCELL(IZ)),
     >         MESHC(4,NUCELL(IZ)),NSURC(NUCELL(IZ)),NREGC(NUCELL(IZ)),
     >         IDIRC(NUCELL(IZ)),NTPIN(NUCELL(IZ)),
     >         ITPIN(3,MAXPIN,0:NUCELL(IZ)))
      ALLOCATE(DCMESH(-1:MAXMSH,4,0:NUCELL(IZ)),
     >         DRAPIN(-1:4,MAXPIN,0:NUCELL(IZ)),REGI(-NFSUR:NFREG),
     >         CELID(NBUCEL),PINID(NBUCEL*MAXPIN))
*
      CALL XDISET(REGI(-NFSUR),NUNK,0)
      CALL XDISET(CELID,NBUCEL,0)
      CALL XDISET(IND2T3(-NFSUR,0),NUNK*(MAXMSH*NUCELL(IZ)+2),0)
      N2SUR=0
      N2REG=0
      N2CEL=0
      N2PIN=0
      LFIRST=.TRUE.
      LSTCEL=.FALSE.
      DO 15 J=1,NUCELL(IY)
      DO 10 I=1,NUCELL(IX)
*----
*  LOOP OVER THE CELLS IN THE PLAN PERPENDICULAR TO THE PROJECTION AXIS
*----
*     ----
*     CELL LEVEL (1)
*     ----
      !write(*,*) 'CELL LEVEL (',I,J,' )'
      DO K=1,NUCELL(IZ)
         ICEL=IUNFLD(1,I,J,K,0)
         ITRN=IUNFLD(2,I,J,K,0)
         IF (ITRN.NE.IUNFLD(2,I,J,1,0))
     1    CALL XABORT('NXT3T2: INVALID PRISMATIC GEOMETRY (TURN).')
*        LOAD THE CONTENTS OF THE DIFFERENT CELLS (I,J,K=1,NUCELL(IZ)) 
         CALL NXTLDC(IPTRK,MAXMSH,ICEL,IDIRC(K),MESHC(1,K),NSURC(K),
     1        NREGC(K),NTPIN(K),DCMESH(-1,1,K),INDEX(1,-MXGSUR,K),
     2        IDREG(1,K),IDSUR(1,K),ITPIN(1,1,K),DRAPIN(-1,1,K))
         !write(*,*) 'loading cell',ICEL,MESHC(1,K),MESHC(2,K),MESHC(4,K)
         IF (K.EQ.1) THEN
            IF (CELID(ICEL).EQ.0) THEN
*           RECOVER DIM INFO FOR THE CORRESPONDING 2D CELL
               LSTCEL=.TRUE.
               N2CEL=N2CEL+1
               WRITE(NAMCEL,'(A1,I8.8)') 'C',ICEL
               !write(*,*) 'copying from ',NAMCEL
               NAMREC=NAMCEL//'DIM'
               CALL XDISET(ESTATE,NSTATE,0)
               CALL LCMGET(IPTRK,NAMREC,ESTATE)
               IF ((ESTATE(1).EQ.21).OR.
     1              (ESTATE(1).EQ.22).OR.
     2              (ESTATE(1).EQ.23)) THEN
                  ESTATE(1)=20
               ELSEIF(ESTATE(1).EQ.7) THEN
                  ESTATE(1)=5
               ENDIF
               ESTATE(5)=0
               ESTATE(6)=0
               ESTATE(12)=N2REG+1
               ESTATE(14)=N2SUR+1
               CELID(ICEL)=N2CEL
            ENDIF
            IUNFLD(1,I,J,1,1)=CELID(ICEL)
            IUNFLD(2,I,J,1,1)=ITRN
         ENDIF
      ENDDO
*     CHECK CELLS COMPATIBILITY, UPDATE IND2T3 FOR THIS SET OF CELLS
*     AND FILL-IN 2D CORRESPONDING CELL CONTENTS
      NTPINR=NTPIN(1)
      DO K=2,NUCELL(IZ)
         IF (NTPIN(K).NE.NTPINR) 
     1    CALL XABORT('NXT3T2: INVALID PRISMATIC GEOMETRY (NTPIN).')
      ENDDO
      CALL NXTPRI(IPTRK,JPTRK,IX,IY,IZ,NFREG,NFSUR,MAXMSH,NUCELL,MXGSUR,
     1     MXGREG,INDEX,IDSUR,IDREG,MESHC,NSURC,NREGC,IDIRC,NZP,N2REG,
     2     N2SUR,IND2T3,REGI,DEPS,DCMESH,ZCORD,LFIRST,LSTCEL,1,
     3     IUNFLD(1,I,J,1,0),N2CEL,N2SURC,N2REGC)
      IF (LSTCEL) THEN        
*     STORE 2D CELL CONTENTS: DIM ARRAY
         ESTATE(10)=N2REGC
         ESTATE(11)=N2SURC
         ESTATE(13)=N2REG
         ESTATE(15)=N2SUR
         !write(*,*) ESTATE(1),ESTATE(2)
         WRITE(NAMCE2,'(A1,I8.8)') 'C',N2CEL
         NAMREC=NAMCE2//'DIM'
         CALL LCMPUT(JPTRK,NAMREC,NSTATE,1,ESTATE)
      ENDIF
*     ----
*     PIN LEVEL (2)
*     ----
      CALL XDISET(PINID,NTPINR,0)
      DO II=1,NTPINR
         !write(*,*) 'PIN LEVEL ( ',II,')'
*        LOAD THE CONTENTS OF THE DIFFERENT PINS (II,K=1,NUCELL(IZ)) 
         IDIRC(1)=ABS(ITPIN(3,II,1))
         IPIN=ITPIN(2,II,1)
         HPIN=DRAPIN(IZ,II,1)
         DELZ=DCMESH(MESHC(IZ,1),IZ,1)-DCMESH(0,IZ,1)
         IF (.NOT.XDDCOM(DELZ,HPIN,DEPS))
     1    CALL XABORT('NXT3T2: INVALID PRISMATIC GEOMETRY (HPIN).')
         CALL NXTLDP(IPTRK,MAXMSH,IPIN,MESHC(1,1),NSURC(1),NREGC(1),
     1        DCMESH(-1,1,1),INDEX(1,-MXGSUR,1),IDREG(1,1),IDSUR(1,1))
         APIN=DRAPIN(-1,II,1)
         RPIN=DRAPIN( 0,II,1)
         RADP=DRAPIN( 4,II,1)
         LSTPIN=.FALSE.
         IF (PINID(IPIN).EQ.0) THEN
*        RECOVER DIM INFO FOR THE CORRESPONDING 2D PIN
            LSTPIN=.TRUE.
            N2PIN=N2PIN+1
            WRITE(NAMCEL,'(A1,I8.8)') 'P',IPIN
            !write(*,*) 'copying from ',NAMCEL
            NAMREC=NAMCEL//'DIM'
            CALL XDISET(ESTATE,NSTATE,0)
            CALL LCMGET(IPTRK,NAMREC,ESTATE)
            IF ((ESTATE(1).EQ.6).OR.
     1          (ESTATE(1).EQ.10).OR.
     2          (ESTATE(1).EQ.11)) THEN
               ESTATE(1)=3
            ENDIF
            ESTATE(5)=0
            ESTATE(6)=0
            ESTATE(12)=N2REG+1
            ESTATE(14)=N2SUR+1
            PINID(IPIN)=N2PIN
         ENDIF
         ITPIN(1,II,0)=ITPIN(1,II,1)
         ITPIN(2,II,0)=PINID(IPIN)
         ITPIN(3,II,0)=3
         DRAPIN(-1,II,0)=APIN
         DRAPIN( 0,II,0)=RPIN
         DRAPIN( 1,II,0)=0.D0
         DRAPIN( 2,II,0)=0.D0
         DRAPIN( 3,II,0)=1.D0
         DRAPIN( 4,II,0)=RADP
         DO K=2,NUCELL(IZ)
            DO JJ=1,NTPINR
               LPIN=.TRUE.
               LPIN=LPIN.AND.(XDDCOM(HPIN,DRAPIN(IZ,JJ,K),DEPS))
               LPIN=LPIN.AND.(XDDCOM(APIN,DRAPIN(-1,JJ,K),DEPS))
               LPIN=LPIN.AND.(XDDCOM(RPIN,DRAPIN( 0,JJ,K),DEPS))
               LPIN=LPIN.AND.(XDDCOM(RADP,DRAPIN( 4,JJ,K),DEPS))
               IF (LPIN) THEN
                  IPIN=ITPIN(2,JJ,K)
                  IDIRC(K)=ABS(ITPIN(3,JJ,K))
                  GOTO 20
               ENDIF
            ENDDO
            CALL XABORT('NXT3T2: INVALID PRISMATIC GEOMETRY (PIN).')
 20         CONTINUE
            CALL NXTLDP(IPTRK,MAXMSH,IPIN,MESHC(1,K),NSURC(K),NREGC(K),
     1           DCMESH(-1,1,K),INDEX(1,-MXGSUR,K),IDREG(1,K),
     2           IDSUR(1,K))
         ENDDO
*        CHECK PINS COMPATIBILITY AND UPDATE IND2T3 FOR THIS SET OF PINS
         CALL NXTPRI(IPTRK,JPTRK,IX,IY,IZ,NFREG,NFSUR,MAXMSH,NUCELL,
     1        MXGSUR,MXGREG,INDEX,IDSUR,IDREG,MESHC,NSURC,NREGC,IDIRC,
     2        NZP,N2REG,N2SUR,IND2T3,REGI,DEPS,DCMESH,ZCORD,LFIRST,
     3        LSTPIN,2,ITPIN(2,II,1),N2PIN,N2SURP,N2REGP)
         IF (LSTPIN) THEN        
*        STORE 2D PIN CONTENTS: DIM ARRAY
            ESTATE(10)=N2REGP
            ESTATE(11)=N2SURP
            ESTATE(13)=N2REG
            ESTATE(15)=N2SUR
            WRITE(NAMCE2,'(A1,I8.8)') 'P',N2PIN
            NAMREC=NAMCE2//'DIM'
            CALL LCMPUT(JPTRK,NAMREC,NSTATE,1,ESTATE)
         ENDIF
      ENDDO
      IF (LSTCEL) THEN
*     STORE 2D CELL CONTENTS: PIN RELATED
         IF (NTPINR.GT.0) THEN
            WRITE(NAMCE2,'(A1,I8.8)') 'C',N2CEL
            NAMREC=NAMCE2//'PIN'
            CALL LCMPUT(JPTRK,NAMREC,6*NTPINR,4,DRAPIN(-1,1,0))
            NAMREC=NAMCE2//'PNT'
            CALL LCMPUT(JPTRK,NAMREC,3*NTPINR,1,ITPIN(1,1,0))
         ENDIF
      ENDIF
      LSTCEL=.FALSE.
      LFIRST=.FALSE.
*----      
 10   CONTINUE
 15   CONTINUE
      N2SUR=-N2SUR
*----
*  FILL IN AND STORE MATALB AND SareaRvolume ARRAYS FOR THE 2D GEOMETRY
*----
      DELZ=ZCORD(1)
      DO I=-N2SUR,N2REG
         MATALB2(I)=MATALB(IND2T3(I,1))
         SURVOL2(I)=SURVOL(IND2T3(I,1))/DELZ
      ENDDO
      NUNK2=N2SUR+N2REG+1
      CALL LCMPUT(JPTRK,'MATALB      ',NUNK2,1,MATALB2(-N2SUR))
      CALL LCMPUT(JPTRK,'SAreaRvolume',NUNK2,4,SURVOL2(-N2SUR))
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(PINID,CELID,REGI,DRAPIN,DCMESH)
      DEALLOCATE(ITPIN,NTPIN,IDIRC,NREGC,NSURC,MESHC,IDREG,IDSUR,INDEX)
      RETURN
      END
