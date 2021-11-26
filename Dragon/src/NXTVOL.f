*DECK NXTVOL
      SUBROUTINE NXTVOL(IPTRK ,IPRINT,MAXMSS,ITYPG ,IDIRC ,IGEO  ,
     >                  ILEV  ,NM    ,NREG  ,NSUR  ,NREGN ,NSURN ,
     >                  MAXPIN,NBPIN ,ITPIN ,DRAPIN,IDREG ,IDSUR ,
     >                  DAMESH,INDXSR,NAREG )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute regional volumes.
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
* IPTRK   pointer to the TRACKING data structure.
* IPRINT  intermediate printing level for output.
* MAXMSS  maximum number of elements in MESH array after split.
* ITYPG   type of geometry.
* IDIRC   direction of cell (1 for XYZ, 2 for YZX and 3 for ZXY).
*         Note that for CAR3D without pins IDIRC=1 while for
*         for CAR3D with pins IDIRC specified by pins direction.
* IGEO    geometry number.
* ILEV    geometry level.
* NM      mesh size in all directions ($X$, $Y$, $Z$ and $R$).
* NREG    maximum number of regions in splitted geometry.
* NSUR    maximum number of surfaces in splitted geometry.
* NREGN   number of regions in splitted geometry after symmetry.
* NSURN   number of surfaces in splitted geometry after symmetry.
* MAXPIN  maximum number of pins.
* NBPIN   number of pins.
* ITPIN   pins identification.
* DRAPIN  pins position.
* IDREG   region identifier after symmetry.
* IDSUR   surface identifier after symmetry.
* DAMESH  final mesh description for geometry.
*
*Parameters: input/output
* NAREG   last region number considered.
*
*Parameters: output
* INDXSR  local indexing of surfaces/regions.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPTRK
      INTEGER          IPRINT,MAXMSS,ITYPG,IGEO,ILEV,NM(4),
     >                 NREG,NSUR,NREGN,NSURN
      INTEGER          MAXPIN,NBPIN
      DOUBLE PRECISION DRAPIN(-1:4,MAXPIN)
      INTEGER          IDREG(NREG),IDSUR(NSUR),ITPIN(3,MAXPIN)
      DOUBLE PRECISION DAMESH(-1:MAXMSS,4)
      INTEGER          NAREG
      INTEGER          INDXSR(5,-NSUR:NREG)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTVOL')
      DOUBLE PRECISION DCUTOF,DCUTOS,DZERO,DONE
      PARAMETER       (DCUTOF=1.0D-8,DCUTOS=1.0D-6,DZERO=0.0D0,
     >                 DONE=1.0D0)
*----
*  Local variables
*----
      INTEGER          NDIM,IDIRC,IDIRCX,NBSUR,NBREG
      CHARACTER        NAMREC*12
      INTEGER          IREG,IDV,ISUR,IDS,INV,INS,LSTREG
      DOUBLE PRECISION VMAX,SMAX
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INREN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SURVOL,SVT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: POSTRI
*----
*  Data
*----
      CHARACTER        CLEV(2)*1
      SAVE             CLEV
      DATA             CLEV /'C','P'/
*----
*  Scratch storage allocation
*   SURVOL  area/volume of regions.
*   SVT     temporary area/volume of regions.
*   INREN   temporary vector for new region/surfaces identification.
*----
      ALLOCATE(INREN(-NSURN:NREGN))
      ALLOCATE(SURVOL(-NSUR:NREG),SVT(-NSURN:NREGN))
      CALL XDDSET(SURVOL,NSUR+NREG+1,DZERO)
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010)
        WRITE(IOUT,6014)  (IDREG(IDV),IDV=1,NREG)
        WRITE(IOUT,6012)
        WRITE(IOUT,6014)  (IDSUR(IDS),IDS=1,NSUR)
      ENDIF
      NDIM=3
      IF(ITYPG .EQ.  5 .OR. ITYPG .EQ. 7) THEN
        IF(ITYPG .EQ. 5) NDIM=2
        CALL NXTVCA(IPRINT,NDIM  ,IDIRC ,MAXMSS,NSUR  ,NREG  ,
     >              NM    ,DAMESH,NBSUR ,NBREG ,INDXSR,SURVOL)
        IF(NBPIN .GT. 0) THEN
*----
*  Remove pin contributions
*----
          CALL NXTPCA(IPRINT,NDIM  ,IDIRC ,MAXMSS,NSUR  ,NREG  ,
     >                NM    ,DAMESH,NBPIN ,ITPIN ,DRAPIN,
     >                NBSUR ,NBREG ,INDXSR,SURVOL)
        ENDIF
      ELSE IF(ITYPG .EQ.  3 .OR. ITYPG .EQ.  6) THEN
        IF(ITYPG .EQ. 3) NDIM=2
        IDIRCX=-IDIRC
        DAMESH(0,1)=-DAMESH(-1,1)-DAMESH(NM(4),4)
        DAMESH(NM(1),1)=-DAMESH(-1,1)+DAMESH(NM(4),4)
        DAMESH(0,2)=-DAMESH(-1,2)-DAMESH(NM(4),4)
        DAMESH(NM(2),2)=-DAMESH(-1,2)+DAMESH(NM(4),4)
        CALL NXTVCC(IPRINT,NDIM  ,IDIRCX,MAXMSS,NSUR  ,NREG ,
     >              NM    ,DAMESH,NBSUR ,NBREG ,INDXSR,SURVOL)
        IF(NBPIN .GT. 0) THEN
*----
*  Remove pin contributions
*----
          CALL NXTPCC(IPRINT,NDIM  ,IDIRCX,MAXMSS,NSUR  ,NREG  ,
     >                NM    ,DAMESH,NBPIN ,ITPIN ,DRAPIN,
     >                NBSUR ,NBREG ,INDXSR,SURVOL)
        ENDIF
      ELSE IF(ITYPG .EQ. 10 ) THEN
        NDIM=3
        IDIRCX=-IDIRC
        DAMESH(0,2)=-DAMESH(-1,2)-DAMESH(NM(4),4)
        DAMESH(NM(2),2)=-DAMESH(-1,2)+DAMESH(NM(4),4)
        DAMESH(0,3)=-DAMESH(-1,3)-DAMESH(NM(4),4)
        DAMESH(NM(3),3)=-DAMESH(-1,3)+DAMESH(NM(4),4)
        CALL NXTVCC(IPRINT,NDIM  ,IDIRCX,MAXMSS,NSUR  ,NREG  ,
     >              NM    ,DAMESH,NBSUR ,NBREG ,INDXSR,SURVOL)
        IF(NBPIN .GT. 0) THEN
*----
*  Remove pin contributions
*----
          CALL NXTPCC(IPRINT,NDIM  ,IDIRCX,MAXMSS,NSUR  ,NREG  ,
     >                NM    ,DAMESH,NBPIN ,ITPIN ,DRAPIN,
     >                NBSUR ,NBREG ,INDXSR,SURVOL)
        ENDIF
      ELSE IF(ITYPG .EQ. 11 ) THEN
        NDIM=3
        IDIRCX=-IDIRC
        DAMESH(0,3)=-DAMESH(-1,3)-DAMESH(NM(4),4)
        DAMESH(NM(3),3)=-DAMESH(-1,3)+DAMESH(NM(4),4)
        DAMESH(0,1)=-DAMESH(-1,1)-DAMESH(NM(4),4)
        DAMESH(NM(1),1)=-DAMESH(-1,1)+DAMESH(NM(4),4)
        CALL NXTVCC(IPRINT,NDIM  ,IDIRCX,MAXMSS,NSUR  ,NREG  ,
     >              NM    ,DAMESH,NBSUR ,NBREG ,INDXSR,SURVOL)
        IF(NBPIN .GT. 0) THEN
*----
*  Remove pin contributions
*----
          CALL NXTPCC(IPRINT,NDIM  ,IDIRCX,MAXMSS,NSUR  ,NREG  ,
     >                NM    ,DAMESH,NBPIN ,ITPIN ,DRAPIN,
     >                NBSUR ,NBREG ,INDXSR,SURVOL)
        ENDIF
      ELSE IF(ITYPG .EQ. 12 .OR. ITYPG .EQ. 13 ) THEN
        IF(ITYPG .EQ. 12) NDIM=2
        CALL NXTVHT(IPRINT,NDIM  ,MAXMSS,NSUR  ,NREG  ,
     >              NM    ,DAMESH,NBSUR ,NBREG ,INDXSR,SURVOL)
        IF(NBPIN .GT. 0) THEN
*----
*  Remove pin contributions
*----
          ALLOCATE(POSTRI(2,3,MAXMSS*MAXMSS,6))
          CALL NXTTLO(IPRINT,MAXMSS,NM    ,DAMESH,POSTRI)
          CALL NXTPHT(IPRINT,NDIM  ,IDIRC ,MAXMSS,NSUR  ,NREG  ,
     >                NM    ,DAMESH,NBPIN ,ITPIN ,DRAPIN,
     >                NBSUR ,NBREG ,INDXSR,SURVOL,POSTRI)
          DEALLOCATE(POSTRI)
        ENDIF
      ELSE IF(ITYPG .EQ. 20 .OR. ITYPG .EQ. 21 .OR.
     >        ITYPG .EQ. 22 .OR. ITYPG .EQ. 23) THEN
        IF(ITYPG .EQ. 20) NDIM=2
        CALL NXTVCC(IPRINT,NDIM  ,IDIRC ,MAXMSS,NSUR  ,NREG  ,
     >              NM    ,DAMESH,NBSUR ,NBREG ,INDXSR,SURVOL)
        IF(NBPIN .GT. 0) THEN
*----
*  Remove pin contributions
*----
          CALL NXTPCC(IPRINT,NDIM  ,IDIRC ,MAXMSS,NSUR  ,NREG  ,
     >                NM    ,DAMESH,NBPIN ,ITPIN ,DRAPIN,
     >                NBSUR ,NBREG ,INDXSR,SURVOL)
        ENDIF
      ELSE IF(ITYPG .EQ. 26 .OR. ITYPG .EQ. 27 ) THEN
        IF(ITYPG .EQ. 26) NDIM=2
        ALLOCATE(POSTRI(2,3,MAXMSS*MAXMSS,6))
        CALL NXTTLO(IPRINT,MAXMSS,NM    ,DAMESH,POSTRI)
        CALL NXTVHC(IPRINT,NDIM  ,MAXMSS,NSUR  ,NREG  ,
     >              NM    ,DAMESH,NBSUR ,NBREG ,INDXSR,SURVOL,
     >              POSTRI)
        IF(NBPIN .GT. 0) THEN
*----
*  Remove pin contributions
*----
          CALL NXTPHC(IPRINT,NDIM  ,MAXMSS,NSUR  ,NREG  ,
     >                NM    ,DAMESH,NBPIN ,ITPIN ,DRAPIN,
     >                NBSUR ,NBREG ,INDXSR,SURVOL,POSTRI)
        ENDIF
        DEALLOCATE(POSTRI)
      ENDIF
*----
*  Save surface and region identification on IPTRK
*----
      CALL XDDSET(SVT,NREGN+NSURN+1,DZERO)
      VMAX=0.0D0
      DO IREG=1,NBREG
        VMAX=MAX(VMAX,SURVOL(IREG))
        IDV=ABS(IDREG(IREG))
        IF(IDV .GT. NREGN) CALL XABORT(NAMSBR//
     >  ': Number of regions insufficient')
        IF(IDV .NE. 0) THEN
          SVT(IDV)=SVT(IDV)+SURVOL(IREG)
        ENDIF
      ENDDO
      SMAX=0.0D0
      DO ISUR=1,NBSUR
        SMAX=MAX(SMAX,SURVOL(-ISUR))
        IDS=ABS(IDSUR(ISUR))
        IF(IDS .GT. NSURN) CALL XABORT(NAMSBR//
     >  ': Number of surfaces insufficient')
        IF(IDS .NE. 0) THEN
          SVT(-IDS)=SVT(-IDS)+SURVOL(-ISUR)
        ENDIF
      ENDDO
*----
*  Remove region/surfaces with 0 volumes
*----
      INV=0
      INREN(0)=0
      DO IDV=1,NREGN
        IF(SVT(IDV)/VMAX .GT. DCUTOF) THEN
          INV=INV+1
          INREN(IDV)=INV
        ELSE
          INREN(IDV)=0
        ENDIF
      ENDDO
      LSTREG=INV+NAREG
      DO IREG=1,NBREG
        IDV=IDREG(IREG)
        INV=INREN(ABS(IDV))
        IF(INV .NE. 0) INV=INV+NAREG
        IF(IDV .LT. 0) THEN
          IDREG(IREG)=-INV
        ELSE
          IDREG(IREG)=INV
        ENDIF
      ENDDO
      INS=0
      DO IDS=1,NSURN
        IF(SVT(-IDS)/SMAX .GT. DCUTOS) THEN
          INS=INS+1
          INREN(-IDS)=INS
        ELSE
          INREN(-IDS)=0
        ENDIF
      ENDDO
      DO ISUR=1,NBSUR
        IDS=IDSUR(ISUR)
        INS=INREN(-ABS(IDS))
        IF(INS .NE. 0) INS=INS
        IF(IDS .LT. 0) THEN
          IDSUR(ISUR)=-INS
        ELSE
          IDSUR(ISUR)=INS
        ENDIF
      ENDDO
      NAREG=LSTREG
      WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'VSE'
      CALL LCMPUT(IPTRK,NAMREC,(NBSUR+NBREG+1),4,SURVOL)
      WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'VSI'
      CALL LCMPUT(IPTRK,NAMREC,(NBSUR+NBREG+1)*5,1,INDXSR)
      WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'RID'
      CALL LCMPUT(IPTRK,NAMREC,NBREG,1,IDREG)
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6011)
        WRITE(IOUT,6014)  (IDREG(IDV),IDV=1,NREG)
        WRITE(IOUT,6013)
        WRITE(IOUT,6014)  (IDSUR(IDS),IDS=1,NBSUR)
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(SVT,SURVOL)
      DEALLOCATE(INREN)
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(' Original regions ID')
 6011 FORMAT(' Final regions ID')
 6012 FORMAT(' Original surfaces ID')
 6013 FORMAT(' Final surfaces ID')
 6014 FORMAT(5I15)
      END
