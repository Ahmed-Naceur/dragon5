*DECK NXTTCG
      SUBROUTINE NXTTCG(IPTRK ,IFTRK ,IPRINT,IGTRK ,NBSLIN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To track an assembly of cells containing
* clusters using the new EXCELL tracking procedure.
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
* IPTRK   pointer to the TRACKING data structure in
*         update or creation mode.
* IFTRK   pointer to the TRACKING file in
*         creation mode.
* IPRINT  print level.
* IGTRK   flag to generate the tracking file. In the case where:
*         IGTRK=1, the tracking is performed and
*         used to evaluate the track normalisation factor and the
*         tracking file is generated; when IGTRK=0, the tracking is
*         still performed and used to evaluate the
*         track normalisation factor but the tracking file is not
*         generated.
* NBSLIN  maximum number of segments in a single tracking line
*         computed by default in NXTTCG but limited to 100000
*         elements. This default value can be bypassed using
*         keyword NBSLIN.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPTRK
      INTEGER          IFTRK
      INTEGER          IPRINT,IGTRK
      INTEGER          NBSLIN
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTTCG')
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
      INTEGER          MXGAUS
      PARAMETER       (MXGAUS=64)
*----
*  Functions
*----
      INTEGER          KDROPN,IFTEMP,KDRCLS,ICLS
*----
*  Local variables
*----
      INTEGER          ISTATE(NSTATE),IEDIMG(NSTATE),ICODE(6)
      REAL             RSTATT(NSTATE),ALBEDO(6)
      INTEGER          RENO,LTRK,AZMOAQ,ISYMM,POLQUA,POLOAQ,AZMQUA,
     >                 AZMNBA,IMU
      DOUBLE PRECISION DENUSR,RCUTOF,DENLIN,SPACLN
      DOUBLE PRECISION WEIGHT
      DOUBLE PRECISION RADIUS,CENTER(3)
      INTEGER          NDIM,ITYPBC,IDIRG,NBOCEL,NBUCEL,IDIAG,
     >                 ISAXIS(3),NOCELL(3),NUCELL(3),MXMSH,MAXMSH,
     >                 MAXREG,NBTCLS,MAXPIN,MAXMSP,MAXRSP,NFSUR,
     >                 NFREG,MXGSUR,MXGREG,NUNK,NPLANE,NPOINT,
     >                 NTLINE,NBTDIR,MAXSUB,MAXSGL,NBDR
      INTEGER          IPER(3)
      INTEGER          JJ,NCOR,NQUAD,NANGL,NBANGL,LINMAX
      DOUBLE PRECISION DQUAD(4),ABSC(3,2),RCIRC
      CHARACTER        CTRK*4,COMENT*80,NAMREC*12
      REAL             XGSS(MXGAUS),WGSS(MXGAUS),XGSS0(2*MXGAUS),
     >                 WGSS0(2*MXGAUS)
      INTEGER          IFMT
      INTEGER          NEREG,NESUR,NUNKC,NNBANG
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:)   :: KEYMRG,MATALB
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IUNFLD,NBSANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)     :: SURVOL,DNSANG,
     >                                                   DDANG,DSNOR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: DGMESH,DDENWT,
     >                                                   DVNOR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DANGLT,DEPART
*----
*  Processing starts:
*  print routine opening header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
*----
*  Open temporary tracking file if required
*----
      IF(IGTRK .EQ. 1) THEN
        IFTEMP= KDROPN('DUMMYSQ',0,2,0)
        IF(IFTEMP .LE. 0) WRITE(IOUT,9010) NAMSBR
      ENDIF
*----
*  Get state vectors
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      CALL LCMGET(IPTRK,'EXCELTRACKOP',RSTATT)
      NEREG=ISTATE(1)
      NESUR=ISTATE(5)
      NUNKC=NEREG+NESUR+1
      RENO  =ISTATE(8)
      LTRK  =ISTATE(9)
      AZMOAQ=ISTATE(11)
      ISYMM =ISTATE(12)
      POLQUA=ISTATE(13)
      POLOAQ=ISTATE(14)
      AZMQUA=ISTATE(15)
      AZMNBA=ISTATE(16)
      IFMT=ISTATE(21)
      DENUSR=DBLE(RSTATT(2))
      RCUTOF=DBLE(RSTATT(3))
      DENLIN=DBLE(RSTATT(4))
      SPACLN=DBLE(RSTATT(5))
      WEIGHT=DBLE(RSTATT(6))
*----
*  Get main tracking records
*----
      CALL LCMGET(IPTRK,'ICODE       ',ICODE )
      CALL LCMGET(IPTRK,'ALBEDO      ',ALBEDO)
      CALL LCMSIX(IPTRK,'NXTRecords  ',1)
*----
*  Get general dimensioning vector for geometry tracking
*----
      CALL XDISET(IEDIMG,NSTATE,0)
      CALL LCMGET(IPTRK,'G00000001DIM',IEDIMG)
      NDIM     =IEDIMG( 1)
      ITYPBC   =IEDIMG( 2)
      IDIRG    =IEDIMG( 3)
      NBOCEL   =IEDIMG( 4)
      NBUCEL   =IEDIMG( 5)
      IDIAG    =IEDIMG( 6)
      ISAXIS(1)=IEDIMG( 7)
      ISAXIS(2)=IEDIMG( 8)
      ISAXIS(3)=IEDIMG( 9)
      NOCELL(1)=IEDIMG(10)
      NOCELL(2)=IEDIMG(11)
      NOCELL(3)=IEDIMG(12)
      NUCELL(1)=IEDIMG(13)
      NUCELL(2)=IEDIMG(14)
      NUCELL(3)=IEDIMG(15)
      MXMSH    =IEDIMG(16)
      MAXREG   =IEDIMG(17)
      NBTCLS   =IEDIMG(18)
      MAXPIN   =IEDIMG(19)
      MAXMSP   =IEDIMG(20)
      MAXRSP   =IEDIMG(21)
      NFSUR    =IEDIMG(22)
      NFREG    =IEDIMG(23)
      MXGSUR   =IEDIMG(24)
      MXGREG   =IEDIMG(25)
      NUNK=NFSUR+NFREG+1
      MAXMSH=MAX(MXMSH,MAXMSP,MAXREG)
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6011) NFREG,NEREG,NFSUR,NESUR
      ENDIF
      IF(ITYPBC .EQ. 0) THEN
*----
*  Define Cell for periodicity
*  Cartesian Boundary
*----
        IPER(1)=2
        IF(ABS(ISAXIS(1)) .EQ. 3) IPER(1)=1
        IPER(2)=2
        IF(ABS(ISAXIS(2)) .EQ. 3) IPER(2)=1
        IPER(3)=2
        IF(ABS(ISAXIS(3)) .EQ. 3) IPER(3)=1
*----
*  Use intrinsic geometry symmetries
*  to simplify tracking unless
*  NOSY tracking option activated
*----
        IF(ISYMM .NE. 0) THEN
          ISYMM=0
          IF(ABS(ISAXIS(1)) .EQ. 1 .OR. ABS(ISAXIS(1)) .EQ. 2) THEN
*----
*  X SYMMETRY
*----
            ISYMM=2
          ENDIF
          IF(ABS(ISAXIS(2)) .EQ. 1 .OR. ABS(ISAXIS(2)) .EQ. 2) THEN
*----
*  Y SYMMETRY
*----
            ISYMM=4+2*ISYMM
          ENDIF
          IF(NDIM .EQ. 3) THEN
            IF(ABS(ISAXIS(3)) .EQ. 1 .OR. ABS(ISAXIS(3)) .EQ. 2) THEN
*----
*  Z SYMMETRY
*----
              ISYMM=16+ISYMM
            ENDIF
          ENDIF
          IF(ISYMM .EQ. 0) ISYMM=1
        ENDIF
      ENDIF
*----
*  Get cell description of geometry
*----
      ALLOCATE(IUNFLD(2,NBUCEL))
      NAMREC='G00000001CUF'
      CALL LCMGET(IPTRK,NAMREC,IUNFLD)
      ALLOCATE(KEYMRG(-NFSUR:NFREG),MATALB(-NFSUR:NFREG),
     >         SURVOL(-NFSUR:NFREG))
      CALL LCMGET(IPTRK,'KEYMRG      ',KEYMRG)
      CALL LCMGET(IPTRK,'MATALB      ',MATALB)
      CALL LCMGET(IPTRK,'SAreaRvolume',SURVOL)
*----
*  Read global mesh for geometry
*----
      ALLOCATE(DGMESH(-1:MAXMSH,4))
      CALL NXTXYZ(IPTRK ,IPRINT,NDIM  ,ITYPBC,MAXMSH,NUCELL,
     >            ABSC,DGMESH)
*----
*  Verify tracking parameters and compute number of angles
*  associated with angular order
*  and spatial quadrature parameters
*  1. Isotropic tracking
*----
      NCOR= 1
      NPLANE=1
      IF(LTRK .EQ. 0) THEN
        IF(NDIM .EQ. 2) THEN
          NQUAD=2
          DQUAD(1)=DONE
          DQUAD(2)=DONE
          NANGL=AZMOAQ
          NBANGL=NANGL
          IF(ISYMM .EQ. 2 .OR. ISYMM .EQ. 8) THEN
            DQUAD(1)=DONE/DTWO
            DQUAD(2)=DZERO
          ENDIF
          ALLOCATE(DANGLT(NDIM,NQUAD,NBANGL),DDENWT(NBANGL,NQUAD))
          CALL NXTQAS(IPRINT,NDIM  ,AZMQUA,NANGL ,NQUAD ,NBANGL,
     >                DQUAD ,DANGLT,DDENWT)
        ELSE IF(NDIM .EQ. 3) THEN
          NQUAD=4
          NANGL=AZMOAQ
          IF(MOD(AZMOAQ,2) .EQ. 1)THEN
            WRITE(IOUT,9000) NAMSBR,NANGL,NANGL+1
            NANGL=NANGL+1
          ENDIF
          IF(AZMQUA .EQ. 1) THEN
            IF(AZMOAQ .GT. 16) THEN
              WRITE(IOUT,9001) NAMSBR
              AZMQUA=4
            ENDIF
            NBANGL=(NANGL * (NANGL+2)) / 8
          ENDIF
          IF(AZMQUA .EQ. 4) THEN
            IF(NANGL .GT. 44) THEN
              WRITE(IOUT,9004) NAMSBR,NANGL
              NANGL=44
            ENDIF
            NBANGL=3*(NANGL * (NANGL+2)) / 8
          ELSE IF(AZMQUA .EQ. 5) THEN
            IF(NANGL .GT. 44) THEN
              WRITE(IOUT,9004) NAMSBR,NANGL
              NANGL=44
            ENDIF
            NBANGL=3*(NANGL*NANGL)/2
          ELSE IF(AZMQUA .EQ. 6) THEN
            IF(NANGL .GT. 20) THEN
              WRITE(IOUT,9005) NAMSBR,NANGL
              NANGL=20
            ENDIF
            NBANGL=(NANGL * (NANGL+2)) / 8
          ELSE IF(AZMQUA .EQ. 7) THEN
            IF(NANGL .GT. 74) THEN
              WRITE(IOUT,9006) NAMSBR,NANGL
              NANGL=74
            ENDIF
            NBANGL=3*(NANGL * (NANGL+2)) / 8
          ENDIF
          AZMOAQ=NANGL
          DQUAD(1)=DONE
          DQUAD(2)=DONE
          DQUAD(3)=DONE
          DQUAD(4)=DONE
          IF( ISYMM .EQ. 8 .OR. ISYMM .EQ. 24 )THEN
            DQUAD(1)=DQUAD(1)/(DTWO*DTWO)
            DQUAD(2)=DZERO
            DQUAD(3)=DZERO
            DQUAD(4)=DZERO
          ELSE IF( ISYMM .EQ. 4 .OR. ISYMM .EQ. 20 )THEN
            DQUAD(1)=DQUAD(1)/DTWO
            DQUAD(2)=DZERO
            DQUAD(3)=DQUAD(3)/DTWO
            DQUAD(4)=DZERO
          ELSE IF( ISYMM .EQ. 2 .OR. ISYMM .EQ. 18 )THEN
            DQUAD(1)=DQUAD(1)/DTWO
            DQUAD(2)=DQUAD(2)/DTWO
            DQUAD(3)=DZERO
            DQUAD(4)=DZERO
          ENDIF
          ALLOCATE(DANGLT(NDIM,NQUAD,NBANGL),DDENWT(NQUAD,NBANGL))
          CALL NXTQAS(IPRINT,NDIM  ,AZMQUA,NANGL ,NQUAD ,NBANGL,
     >                DQUAD ,DANGLT,DDENWT)
        ENDIF
        LINMAX=NBUCEL*
     >         ((MAXPIN+1)*(2*MXGREG+2)+MXGSUR+16)
*----
*  Select standard spatial tracking parameters
*----
        CALL NXTQSS(IPRINT,NDIM  ,ITYPBC,MAXMSH,NUCELL,DENUSR,
     >              DGMESH,NPLANE,NPOINT,DENLIN,SPACLN,
     >              WEIGHT,RADIUS,CENTER)
        NNBANG=NBANGL
      ELSE IF( LTRK.EQ.1 )THEN
*----
*  2. Specular tracking
*----
        NPOINT=0
        NQUAD=1
        IF(NDIM .EQ. 3) CALL XABORT(NAMSBR//
     >  ': TSPC option not valid for 3-D geometries')
        IF(AZMOAQ .GT. 24) THEN
          WRITE(IOUT,9002) NAMSBR,AZMOAQ,24,30
          AZMOAQ=30
        ELSE IF(AZMOAQ .GT. 20) THEN
          IF(AZMOAQ .NE. 24) THEN
            WRITE(IOUT,9003) NAMSBR,AZMOAQ,20,24,24
            AZMOAQ=24
          ENDIF
        ELSE IF(AZMOAQ .GT. 18) THEN
          IF(AZMOAQ .NE. 20) THEN
            WRITE(IOUT,9003) NAMSBR,AZMOAQ,18,20,20
            AZMOAQ=20
          ENDIF
        ELSE IF(AZMOAQ .GT. 14) THEN
          IF(AZMOAQ .NE. 18) THEN
            WRITE(IOUT,9003) NAMSBR,AZMOAQ,14,18,18
            AZMOAQ=18
          ENDIF
        ELSE IF(AZMOAQ .GT. 12) THEN
          IF(AZMOAQ .NE. 14) THEN
            WRITE(IOUT,9003) NAMSBR,AZMOAQ,12,14,14
            AZMOAQ=14
          ENDIF
        ELSE IF(AZMOAQ .GT. 8) THEN
          IF(AZMOAQ .NE. 12) THEN
            WRITE(IOUT,9003) NAMSBR,AZMOAQ,8,12,12
            AZMOAQ=12
          ENDIF
        ELSE IF(AZMOAQ .GE. 0) THEN
          IF(AZMOAQ .NE. 8) THEN
            WRITE(IOUT,9003) NAMSBR,AZMOAQ,-1,8,8
            AZMOAQ=8
          ENDIF
        ENDIF
        NBANGL=AZMOAQ
        NANGL =AZMOAQ
        ALLOCATE(NBSANG(5,NBANGL))
        ALLOCATE(DANGLT(NDIM,NBANGL,4),DDENWT(NBANGL,4),
     >           DNSANG(NBANGL),DDANG(NBANGL))
        LINMAX=8*NANGL*NBUCEL*
     >         ((MAXPIN+1)*(2*MXGREG+2)+MXGSUR+16)
        RCIRC=SQRT(ABSC(1,1)**2+ABSC(2,1)**2)
        ABSC(1,1)= ABSC(1,1)/RCIRC
        ABSC(2,1)= ABSC(2,1)/RCIRC
        CALL NXTQAC(IPRINT,NDIM  ,NANGL ,NBANGL,DENUSR,ABSC  ,
     >              RCIRC ,AZMQUA,IPER  ,DANGLT,DDENWT,DNSANG,
     >              NBSANG,DDANG)
        DEALLOCATE(DDANG)
        DO JJ=1,NBANGL
          DANGLT(1,NBANGL-JJ+1,2)=-DANGLT(1,JJ,1)
          DANGLT(2,NBANGL-JJ+1,2)=DANGLT(2,JJ,1)
          DDENWT(NBANGL-JJ+1,2)=DDENWT(JJ,1)
        ENDDO
        DO JJ=1,NBANGL
          DANGLT(1,NBANGL-JJ+1,4)=DANGLT(1,JJ,1)
          DANGLT(2,NBANGL-JJ+1,4)=-DANGLT(2,JJ,1)
          DDENWT(NBANGL-JJ+1,4)=DDENWT(JJ,1)
          DANGLT(1,NBANGL-JJ+1,3)=DANGLT(1,JJ,2)
          DANGLT(2,NBANGL-JJ+1,3)=-DANGLT(2,JJ,2)
          DDENWT(NBANGL-JJ+1,3)=DDENWT(JJ,2)
        ENDDO
        NNBANG=4*NBANGL
        CALL LCMPUT(IPTRK,'TrackingDirc',NDIM*NNBANG,4,DANGLT)
        CALL LCMPUT(IPTRK,'TrackingTrkW',NNBANG,4,DDENWT)
        CALL LCMPUT(IPTRK,'TrackingSpaD',NBANGL,4,DNSANG)
*----
*  Select cyclic spatial tracking parameters
*----
        ALLOCATE(DEPART(NDIM,2,NBANGL))
        CALL NXTQSC(IPRINT,NDIM  ,NBANGL,MAXMSH,NUCELL,
     >              DGMESH,DANGLT,DDENWT,
     >              DNSANG,NBSANG,DEPART)
        CALL LCMPUT(IPTRK,'TrackingNbST',5*NBANGL,1,NBSANG)
        CALL LCMPUT(IPTRK,'TrackinDepT',NDIM*2*NBANGL,4,DEPART)
        IF(POLQUA .GT. 0) THEN
          IF(POLOAQ .LE. 0) THEN
            IF(ISYMM.LE.1) THEN
              POLOAQ=(AZMOAQ+1)/2
            ELSE IF(ISYMM.GE.2) THEN
              POLOAQ=AZMOAQ
            ENDIF
          ENDIF
          IF( POLQUA.EQ.0 )THEN
            CALL ALGPT(2*POLOAQ, -1.0, 1.0, XGSS0, WGSS0)
            DO IMU=1,POLOAQ
              XGSS(POLOAQ-IMU+1)=XGSS0(IMU)
              WGSS(POLOAQ-IMU+1)=WGSS0(IMU)
            ENDDO
          ELSE
            CALL ALCACT(POLQUA, POLOAQ, XGSS, WGSS)
          ENDIF
          CALL LCMPUT(IPTRK,'POLAR MU    ',POLOAQ,2,XGSS)
          CALL LCMPUT(IPTRK,'POLAR WEIGHT',POLOAQ,2,WGSS)
        ENDIF
      ENDIF
      RSTATT(4)=REAL(DENLIN)
      RSTATT(5)=REAL(SPACLN)
      RSTATT(6)=REAL(WEIGHT)
      RSTATT(7)=REAL(RADIUS)
      RSTATT(8)=REAL(CENTER(1))
      RSTATT(9)=REAL(CENTER(2))
      IF(NDIM .EQ. 3) RSTATT(10)=REAL(CENTER(3))
*----
*  Track
*----
      LINMAX=MIN(LINMAX,NBSLIN)
      IF(IPRINT .GE. 10) WRITE(IOUT,6010) LINMAX
      NBDR=1
      IF(RENO .EQ. -1) THEN
        NBDR=NQUAD*NNBANG+1
      ENDIF
      ALLOCATE(DVNOR(NFREG,NBDR),DSNOR(NFSUR))
      NBTDIR=0
      IF(LTRK .EQ. 0) THEN
*----
*  Standard (isotropic) tracking (white boundary conditions)
*----
        IF(IPRINT .GE. 1) THEN
          IF(NDIM .EQ.2 ) THEN
            WRITE(IOUT,6030) NBANGL*NQUAD,NPOINT
          ELSE
            WRITE(IOUT,6030) NBANGL*NQUAD,NPOINT*NPOINT*NPLANE
          ENDIF
        ENDIF
        MAXSUB=1
        CALL NXTTLS(IPTRK ,IFTEMP,IPRINT,IGTRK ,NDIM  ,MAXMSH,
     >              NFSUR ,NFREG ,NUCELL,NBUCEL,NBANGL,NQUAD ,
     >              NPLANE,NPOINT,LINMAX,MXGSUR,MXGREG,RENO  ,
     >              MAXPIN,NBTDIR,NBDR  ,ITYPBC,IFMT  ,
     >              RCUTOF,SPACLN,WEIGHT,RADIUS,CENTER,
     >              IUNFLD,SURVOL,DGMESH,DANGLT,DDENWT,
     >              MAXSGL,NTLINE,DVNOR ,DSNOR)
      ELSE
*----
*  Cyclic (specular) tracking (mirror like boundary conditions)
*----
        CALL NXTTLC(IPTRK ,IFTEMP,IPRINT,IGTRK ,NDIM  ,MAXMSH,
     >              NFSUR ,NFREG ,NUCELL,NBUCEL,NBANGL,
     >              LINMAX,MXGSUR,MXGREG,RENO  ,NBDR  ,ITYPBC,
     >              IFMT  ,MAXPIN,AZMQUA,IPER  ,IUNFLD,SURVOL,
     >              DGMESH,DANGLT,DDENWT,DNSANG,NBSANG,DEPART,
     >              MAXSUB,MAXSGL,NTLINE,DVNOR ,DSNOR)
      ENDIF
*----
*  Save track normalisation vector
*----
      CALL LCMPUT(IPTRK,'VTNormalize ',NFREG,4,DVNOR)
      IF(NBDR .GT. 1) THEN
        CALL LCMPUT(IPTRK,'VTNormalizeD',NFREG*(NBDR-1),4,
     >              DVNOR(1,2))
      ENDIF
      CALL LCMSIX(IPTRK,'NXTRecords  ',2)
      ISTATE(12)=ISYMM
      ISTATE(14)=POLOAQ
      ISTATE(17)=NPOINT
      ISTATE(18)=LINMAX
      ISTATE(19)=NTLINE
      ISTATE(20)=NBTDIR
      ISTATE(21)=NQUAD*NNBANG
      ISTATE(22)=NPLANE
      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPTRK,'EXCELTRACKOP',NSTATE,2,RSTATT)
*----
*  Renormalize tracks if required and transfer to final tracking file
*----
      IF(IGTRK .EQ. 1) THEN
        CTRK  = '$TRK'
        WRITE(IFTRK) CTRK,5,NTLINE,IFMT
        COMENT='CREATOR     : DRAGON'
        WRITE(IFTRK) COMENT
        COMENT='MODULE      : NXTTCG'
        WRITE(IFTRK) COMENT
        COMENT='TYPE        : CARTESIAN'
        WRITE(IFTRK) COMENT
        IF(RENO .EQ. -1) THEN
          COMENT='TRKNOR      : Directional '
        ELSE IF(RENO .EQ. 0) THEN
          COMENT='TRKNOR      : Global      '
        ELSE
          COMENT='TRKNOR      : Off         '
        ENDIF
        WRITE(IFTRK) COMENT
        IF(IFMT .EQ. 1) THEN
          COMENT='OPTION      : Extended    '
          WRITE(IFTRK) COMENT
        ELSE
          COMENT='OPTION      : Short       '
          WRITE(IFTRK) COMENT
        ENDIF
*----
*  Compress VOLSUR and MATALB according to KEYMRG
*  and save on IFTRK
*----
        WRITE(IFTRK) NDIM,LTRK,NEREG,NESUR,6,NCOR,NQUAD*NNBANG,MAXSUB,
     >  MAXSGL
        CALL NXTCVM(IFTRK,IPRINT,NFREG,NFSUR,NEREG,NESUR,MATALB,SURVOL,
     >              KEYMRG)
        WRITE(IFTRK) ( ICODE(JJ),JJ=1,6)
        WRITE(IFTRK) (ALBEDO(JJ),JJ=1,6)
        CALL NXTSQD(IFTRK,IPRINT,NDIM,NQUAD,NNBANG,DANGLT,DDENWT)
        REWIND IFTEMP
        CALL NXTTNS(IFTRK ,IFTEMP,IPRINT,RENO  ,NFSUR ,NFREG ,
     >              NDIM  ,MAXSUB,MAXSGL,NTLINE,NBDR  ,IFMT  ,
     >              KEYMRG,DVNOR)
*----
*  Close temporary tracking file if required
*----
        ICLS=KDRCLS(IFTEMP,2)
        IF(ICLS .NE. 0) WRITE(IOUT,9011) NAMSBR
      ENDIF
*----
*  Deallocate memory
*----
      DEALLOCATE(DSNOR,DVNOR)
      IF(LTRK .EQ. 1) DEALLOCATE(DEPART,NBSANG,DNSANG)
      DEALLOCATE(DDENWT,DANGLT)
      DEALLOCATE(DGMESH,SURVOL,MATALB,KEYMRG,IUNFLD)
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 1)  WRITE(IOUT,6012)
      IF(IPRINT .GE. 10) WRITE(IOUT,6001) NAMSBR
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(' Maximum length of a line =',I10)
 6011 FORMAT(/' Tracking of geometry begins:'//
     >       ' Number of regions before merge =',I10/
     >       ' Number of regions after merge  =',I10/
     >       ' Number of surfaces before merge=',I10/
     >       ' Number of surfaces after merge =',I10)
 6012 FORMAT(/' Tracking of geometry completed'/)
 6030 FORMAT(' Number of directions for tracking = ',I10/
     >       ' Number of lines per direction     = ',I10)
 9000 FORMAT(' ***** Warning in ',A6,' *****'/
     >       '       Order of quadrature must be even '/
     >       '       Order ',I10,' replaced by ',I10)
 9001 FORMAT(' ***** Warning in ',A6,' *****'/
     >       '       Order of quadrature for EQW limited to 16'/
     >       '       Use PNTN quadrature instead')
 9002 FORMAT(' ***** Warning in ',A6,' *****'/
     >       '       Number of speculular angles requested :',I10/
     >       '       For values > ',I10,' use ',I10)
 9003 FORMAT(' ***** Warning in ',A6,' *****'/
     >       '       Number of speculular angles requested :',I10/
     >       '       For values > ',I10,' and < ',I10,' use ',I10)
 9004 FORMAT(' ***** Warning in ',A6,' *****'/
     >       '       Order of quadrature for PNTN limited to 44'/
     >       '       Order ',I10,' replaced by 44')
 9005 FORMAT(' ***** Warning in ',A6,' *****'/
     >       '       Order of quadrature for LSN limited to 20'/
     >       '       Order ',I10,' replaced by 20')
 9006 FORMAT(' ***** Warning in ',A6,' *****'/
     >       '       Order of quadrature for QRN limited to 74'/
     >       '       Order ',I10,' replaced by 74')
 9010 FORMAT(' ***** Warning in ',A6,' *****'/
     >       '       Impossible to open temporary tracking file ')
 9011 FORMAT(' ***** Warning in ',A6,' *****'/
     >       '       Impossible to close temporary tracking file ')
      END
