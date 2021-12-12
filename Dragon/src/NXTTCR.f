*DECK NXTTCR
      SUBROUTINE NXTTCR(IPTRK ,IPRINT,ICEL  ,ITRN  ,IFSUR ,
     >                  NDIM  ,MAXMSH,LINMAX,MXGSUR,MXGREG,
     >                  MAXPIN,CELLPO,TRKLIM,TRKORI,ANGLES,
     >                  IBLIN ,IELIN ,NUMERO,DLENGT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To track a cell rotated according to its explicit
* position in the assembly.
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
* IPRINT  print level.
* ICEL    cell number.
* ITRN    cell rotation.
* IFSUR   surface treatment with:
*         IFSUR=-1 or +1 surfaces located
*         at the beginning or end of the track are considered;
*         IFSUR=0 no surface is considered.
* NDIM    problem dimensions.
* MAXMSH  maximum number of elements in MESH array.
* LINMAX  maximum number of segments in a track.
* MXGSUR  maximum number of surfaces for any geometry.
* MXGREG  maximum number of region for any geometry.
* MAXPIN  maximum number of pins in a cell.
* CELLPO  global cell position in space.
* TRKLIM  beginning and end of track in this cell.
* TRKORI  track origin.
* ANGLES  track direction.
* IBLIN   track line starting point.
*
*Parameters: output
* IELIN   track line ending point.
* NUMERO  region/surface identification number
*         for segment.
* DLENGT  spatial location of each line segment.
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
      INTEGER          IPRINT,ICEL,ITRN,IFSUR,
     >                 NDIM,MAXMSH,LINMAX,MXGSUR,MXGREG,MAXPIN
      DOUBLE PRECISION CELLPO(3,2),TRKLIM(2),TRKORI(NDIM),
     >                 ANGLES(NDIM)
      INTEGER          IBLIN,IELIN,NUMERO(LINMAX)
      DOUBLE PRECISION DLENGT(LINMAX)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTTCR')
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI
      INTEGER          NXTLCA,NXTLHA,NXTLHT,IRLA,NXTLCY,IRLCY
*----
*  Local variables
*----
      INTEGER          ITST,ILEV,IEDIMC(NSTATE),IEDIMP(NSTATE),
     >                 ITYPG,MESHC(5),MXSLIN,NTPIN,IPIN,MESHP(5),
     >                 IDIR,IDIRC,IDIRP,MXM1,NREGC,NSURC,NREGP,NSURP,
     >                 ICOMB,ITL,IOVER,MESHSP(5)
      INTEGER          NBCOR(2,3),NBSINT(3)
      DOUBLE PRECISION TRKORT(3),ANGROT(3),TRKORR(3),
     >                 TRKORP(3),COSDIR(3),
     >                 PINPOS(-1:1,5),ROTAX,COSR,SINR,
     >                 TRKROF(3),ANGROF(3),ENDNEW
      CHARACTER        NAMCEL*9,NAMREC*12
      INTEGER          IPRLOC,IKLIN,KFSUR
*----
*  Allocatable arrays
*   IDSUR   local surface identifier.
*   IDREG   local region identifier.
*   ITPIN   pin type identifier.
*   DCMESH  meshing vector for geometries.
*   DRAPIN  pin position identifier.
*   ICSINT  identification of spatial position for each
*           line segment in cell description of geometry.
*   DCSINT  position of each intersection point for each
*           line segment in cell description of geometry.
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDSUR,IDREG
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ITPIN
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: ICSINT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DCMESH,DRAPIN,
     > DCSINT
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: VSI
*----
*  Data
*----
      CHARACTER        CDIR(4)*1
      SAVE             CDIR
      CHARACTER        CLEV(2)*1
      SAVE             CLEV
      DATA             CDIR /'X','Y','Z','R'/
      DATA             CLEV /'C','P'/
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IPRLOC=IPRINT
      IF(IPRLOC .GT. 1000) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      ALLOCATE(IDSUR(MXGSUR),IDREG(MXGREG),ITPIN(3,MAXPIN+1),
     > ICSINT(0:5,4*(MXGREG+4),3))
      ALLOCATE(DCMESH(-1:MAXMSH,5),DRAPIN(-1:4,MAXPIN+1),
     > DCSINT(4*(MXGREG+4),3))
      PI=XDRCST('Pi',' ')
      ITST=1
      ILEV=1
      KFSUR=IFSUR
*----
*  Read cell information
*----
      WRITE(NAMCEL,'(A1,I8.8)') CLEV(ILEV),ICEL
      NAMREC=NAMCEL//'DIM'
      CALL XDISET(IEDIMC,NSTATE,0)
      CALL XDDSET(DCMESH(-1,1),(MAXMSH+2)*5,DZERO)
      CALL LCMGET(IPTRK,NAMREC,IEDIMC)
      ITYPG=IEDIMC(1)
      MESHC(1)=IEDIMC(3)
      MESHC(2)=IEDIMC(4)
      MESHC(3)=IEDIMC(5)
      MESHC(4)=IEDIMC(2)
      MESHC(5)=IEDIMC(3)
      NREGC=IEDIMC(8)
      NSURC=IEDIMC(9)
      MXSLIN=4*(MXGREG+4)
      NTPIN=IEDIMC(18)
      NAMREC=NAMCEL//'RID'
      CALL LCMGET(IPTRK,NAMREC,IDREG)
      NAMREC=NAMCEL//'SID'
      CALL LCMGET(IPTRK,NAMREC,IDSUR)
      DO IDIR=1,4
        NAMREC=NAMCEL//'SM'//CDIR(IDIR)
        IF(MESHC(IDIR) .GT. 0) THEN
          CALL LCMGET(IPTRK,NAMREC,DCMESH(-1,IDIR))
        ENDIf
      ENDDO
      IF(NTPIN .GT. 0) THEN
        NAMREC=NAMCEL//'PIN'
        CALL LCMGET(IPTRK,NAMREC,DRAPIN)
        NAMREC=NAMCEL//'PNT'
        CALL LCMGET(IPTRK,NAMREC,ITPIN)
      ENDIF
*----
*  Translate TRKORI in such a way that it is now defined with
*  respect to the cell center.
*----
      DO IDIR=1,NDIM
        TRKORT(IDIR)=TRKORI(IDIR)
     >              -(CELLPO(IDIR,1)+CELLPO(IDIR,2))/DTWO
      ENDDO
*----
*  Rotate tracking line according to ITRN provided
*----
      IF(IPRLOC .GT. 1000) THEN
        WRITE(IOUT,6030) ICEL,ITRN
      ENDIF
      CALL XDDSET(TRKORR,3,DZERO)
      CALL XDDSET(ANGROT,3,DZERO)
      CALL NXTRTL(IPRLOC,NDIM  ,ITRN  ,TRKORT,ANGLES,
     >            TRKORR,ANGROT)
      ITL=1
      IF(ITYPG .EQ.  5 .OR. ITYPG .EQ.  7 .OR.
     >   ITYPG .EQ. 20 .OR. ITYPG .EQ. 21 .OR.
     >   ITYPG .EQ. 22 .OR. ITYPG .EQ. 23) THEN
*----
*  Track Cartesian cell
*----
        IRLA=NXTLCA(IPRLOC   ,ITST  ,NDIM  ,MAXMSH,MXSLIN,
     >              MESHC ,TRKORR,ANGROT,DCMESH,
     >              NBCOR(1,ITL)   ,NBSINT(ITL)    ,
     >              ICSINT(0,1,ITL),DCSINT(1,ITL)  )
*----
*  If required, track annular mesh
*  and combine with Cartesian mesh.
*----
        IRLCY=0
        IDIRC=MOD(IEDIMC(20)+1,3)+1
        IF(ITYPG .EQ. 20 .OR. ITYPG .EQ. 21 .OR.
     >     ITYPG .EQ. 22 .OR. ITYPG .EQ. 23) THEN
          ITL=2
          IRLCY=NXTLCY(IPRLOC,ITST  ,NDIM  ,MAXMSH,MXSLIN,
     >                 MESHC ,TRKORR,ANGROT,DCMESH,IDIRC ,
     >                 NBCOR(1,ITL)   ,NBSINT(ITL)    ,
     >                 ICSINT(0,1,ITL),DCSINT(1,ITL)  )
          IF(IRLCY .EQ. 2) THEN
            ICOMB=1
            CALL NXTLCU(IPRLOC,MXSLIN,ICOMB ,
     >                  NBCOR ,NBSINT,ICSINT,DCSINT)
            ITL=3
          ELSE
            ITL=1
          ENDIF
        ENDIF
      ELSE IF(ITYPG .EQ.  8 .OR. ITYPG .EQ.  9 .OR.
     >        ITYPG .EQ. 24 .OR. ITYPG .EQ. 25) THEN
*----
*  Track hexagonal cell
*----
        IRLA=NXTLHA(IPRLOC   ,ITST  ,NDIM  ,MAXMSH,MXSLIN,
     >              MESHC ,TRKORR,ANGROT,DCMESH,
     >              NBCOR(1,ITL)   ,NBSINT(ITL)    ,
     >              ICSINT(0,1,ITL),DCSINT(1,ITL)  )
*----
*  If required, track annular mesh
*  and combine with Cartesian mesh.
*----
        IRLCY=0
        IF(ITYPG .EQ. 24 .OR. ITYPG .EQ. 25) THEN
          ITL=2
          IDIRC=3
          IRLCY=NXTLCY(IPRLOC,ITST  ,NDIM  ,MAXMSH,MXSLIN,
     >                 MESHC ,TRKORR,ANGROT,DCMESH,IDIRC ,
     >                 NBCOR(1,ITL)   ,NBSINT(ITL)    ,
     >                 ICSINT(0,1,ITL),DCSINT(1,ITL)  )
          IF(IRLCY .EQ. 2) THEN
            ICOMB=1
            CALL NXTLCU(IPRLOC,MXSLIN,ICOMB ,
     >                  NBCOR ,NBSINT,ICSINT,DCSINT)
            ITL=3
          ELSE
            ITL=1
          ENDIF
        ENDIF
      ELSE IF(ITYPG .EQ. 12 .OR. ITYPG .EQ. 13 .OR.
     >        ITYPG .EQ. 26 .OR. ITYPG .EQ. 27 ) THEN
*----
*  Track hexagonal cell with triangular mesh
*----
        MESHC(1)=2*MESHC(1)
        IRLA=NXTLHT(IPRLOC   ,ITST  ,NDIM  ,MAXMSH,MXSLIN,
     >              MESHC ,TRKORR,ANGROT,DCMESH,TRKLIM,
     >              NBCOR(1,ITL)   ,NBSINT(ITL)    ,
     >              ICSINT(0,1,ITL),DCSINT(1,ITL)  )
        IRLCY=0
        IF(ITYPG .EQ. 26 .OR. ITYPG .EQ. 27) THEN
          ITL=2
          IDIRC=3
          IRLCY=NXTLCY(IPRLOC,ITST  ,NDIM  ,MAXMSH,MXSLIN,
     >                 MESHC ,TRKORR,ANGROT,DCMESH,IDIRC ,
     >                 NBCOR(1,ITL)   ,NBSINT(ITL)    ,
     >                 ICSINT(0,1,ITL),DCSINT(1,ITL)  )
          IF(IRLCY .EQ. 2) THEN
            ICOMB=1
            CALL NXTLCU(IPRLOC,MXSLIN,ICOMB ,
     >                  NBCOR ,NBSINT,ICSINT,DCSINT)
            ITL=3
          ELSE
            ITL=1
          ENDIF
        ENDIF
      ENDIF
*----
*  Add tracking lines at the end of current tracking vector
*----
      IOVER=0
      IF(ITYPG .EQ. 12 .OR. ITYPG .EQ. 13 .OR.
     >   ITYPG .EQ. 26 .OR. ITYPG .EQ. 27) THEN
        ALLOCATE(VSI(5,-NSURC:NREGC))
        NAMREC=NAMCEL//'VSI'
        CALL LCMGET(IPTRK,NAMREC,VSI)
        CALL NXTLRH(IPRLOC,IOVER ,ITYPG ,LINMAX,MXSLIN,
     >              NREGC ,NSURC ,MESHC ,IDSUR ,IDREG ,
     >              NBCOR(1,ITL)   ,NBSINT(ITL)    ,
     >              ICSINT(0,1,ITL),DCSINT(1,ITL)  ,
     >              IBLIN ,IELIN ,NUMERO,DLENGT,VSI)
        DEALLOCATE(VSI)
      ELSE
        CALL NXTLRS(IPRLOC,IOVER ,ITYPG ,LINMAX,MXSLIN,
     >              NREGC ,NSURC ,MESHC ,IDSUR ,IDREG ,
     >              NBCOR(1,ITL)   ,NBSINT(ITL)    ,
     >              ICSINT(0,1,ITL),DCSINT(1,ITL)  ,
     >              IBLIN ,IELIN ,NUMERO,DLENGT)
      ENDIF
      IF(IELIN .GT. LINMAX) THEN
        CALL XABORT(NAMSBR//
     >  ': Tracking vector dimensions too small')
      ENDIF
*----
*  Translate TRKORR to take into account OFFCEN for cylinders
*  and pins.
*----
      IF(NTPIN .GT. 0) THEN
        DO IDIR=1,NDIM
          TRKORR(IDIR)=TRKORR(IDIR)-DCMESH(-1,IDIR)
        ENDDO
*----
*  Scan PINS for intersection with LINE
*----
        ILEV=2
        MESHP(1)=1
        MESHP(2)=1
        MESHP(3)=1
        MESHP(4)=1
        MESHP(5)=1
        MXM1=1
        IOVER=1
        DO IPIN=1,NTPIN
          IF(IPRLOC .GT. 1000) THEN
            WRITE(IOUT,6010) IPIN
          ENDIF
*----
*  Translate TRKORT in such a way that it is now defined with
*  respect to the pin center.
*----
          CALL XDDSET(PINPOS,15,DZERO)
          IDIRP=ABS(ITPIN(3,IPIN))
          PINPOS(1,4)=DRAPIN(4,IPIN)
          IF(IDIRP .EQ. 3) THEN
*----
*  2-D or Z directed pins
*----
            COSDIR(1)=DRAPIN(0,IPIN)*COS(DRAPIN(-1,IPIN))
            COSDIR(2)=DRAPIN(0,IPIN)*SIN(DRAPIN(-1,IPIN))
            COSDIR(3)=DZERO
            PINPOS(0,3)=-DRAPIN(3,IPIN)/DTWO
            PINPOS(1,3)=DRAPIN(3,IPIN)/DTWO
            DO IDIR=1,NDIM
              TRKORP(IDIR)=TRKORR(IDIR)-COSDIR(IDIR)
            ENDDO
          ELSE IF(IDIRP .EQ. 2) THEN
*----
*  Y directed pins
*----
            COSDIR(3)=DRAPIN(0,IPIN)*COS(DRAPIN(-1,IPIN))
            COSDIR(1)=DRAPIN(0,IPIN)*SIN(DRAPIN(-1,IPIN))
            COSDIR(2)=DZERO
            PINPOS(0,2)=-DRAPIN(2,IPIN)/DTWO
            PINPOS(1,2)=DRAPIN(2,IPIN)/DTWO
            DO IDIR=1,NDIM
              TRKORP(IDIR)=TRKORR(IDIR)-COSDIR(IDIR)
            ENDDO
          ELSE
*----
*  X directed pins
*----
            COSDIR(2)=DRAPIN(0,IPIN)*COS(DRAPIN(-1,IPIN))
            COSDIR(3)=DRAPIN(0,IPIN)*SIN(DRAPIN(-1,IPIN))
            COSDIR(1)=DZERO
            PINPOS(0,1)=-DRAPIN(1,IPIN)/DTWO
            PINPOS(1,1)=DRAPIN(1,IPIN)/DTWO
            DO IDIR=1,NDIM
              TRKORP(IDIR)=TRKORR(IDIR)-COSDIR(IDIR)
            ENDDO
          ENDIF
          ITL=2
          IRLCY=NXTLCY(IPRLOC,ITST  ,NDIM  ,MXM1  ,MXSLIN,
     >                 MESHP ,TRKORP,ANGROT,PINPOS,IDIRP ,
     >                 NBCOR(1,ITL)   ,NBSINT(ITL)    ,
     >                 ICSINT(0,1,ITL),DCSINT(1,ITL)  )
          IF(IRLCY .EQ. 2) THEN
*----
*  Rotate geometry by (Pi/2-alpha)
*  for pin at alpha.
*----
            ROTAX=PI/DTWO-DRAPIN(-1,IPIN)
            COSR=COS(ROTAX)
            SINR=SIN(ROTAX)
            IF(IDIRP .EQ. 3) THEN
*----
*  2-D or Z directed pins
*----
              TRKROF(1)=TRKORP(1)*COSR-TRKORP(2)*SINR
              TRKROF(2)=TRKORP(1)*SINR+TRKORP(2)*COSR
              TRKROF(3)=TRKORP(3)
              ANGROF(1)=ANGROT(1)*COSR-ANGROT(2)*SINR
              ANGROF(2)=ANGROT(1)*SINR+ANGROT(2)*COSR
              ANGROF(3)=ANGROT(3)
            ELSE IF(IDIRP .EQ. 2) THEN
              TRKROF(3)=TRKORP(3)*COSR-TRKORP(1)*SINR
              TRKROF(1)=TRKORP(3)*SINR+TRKORP(1)*COSR
              TRKROF(2)=TRKORP(2)
              ANGROF(3)=ANGROT(3)*COSR-ANGROT(1)*SINR
              ANGROF(1)=ANGROT(3)*SINR+ANGROT(1)*COSR
              ANGROF(2)=ANGROT(2)
            ELSE
              TRKROF(2)=TRKORP(2)*COSR-TRKORP(3)*SINR
              TRKROF(3)=TRKORP(2)*SINR+TRKORP(3)*COSR
              TRKROF(1)=TRKORP(1)
              ANGROF(2)=ANGROT(2)*COSR-ANGROT(3)*SINR
              ANGROF(3)=ANGROT(2)*SINR+ANGROT(3)*COSR
              ANGROF(1)=ANGROT(1)
            ENDIF
*----
*  Read pin information
*----
            WRITE(NAMCEL,'(A1,I8.8)') CLEV(ILEV),ITPIN(2,IPIN)
            NAMREC=NAMCEL//'DIM'
            CALL XDISET(IEDIMP,NSTATE,0)
            CALL LCMGET(IPTRK,NAMREC,IEDIMP)
            ITYPG=IEDIMP(1)
            MESHSP(1)=IEDIMP(3)
            MESHSP(2)=IEDIMP(4)
            MESHSP(3)=IEDIMP(5)
            MESHSP(4)=IEDIMP(2)
            MESHSP(5)=IEDIMP(3)
            NREGP=IEDIMP(8)
            NSURP=IEDIMP(9)
            NAMREC=NAMCEL//'RID'
            CALL LCMGET(IPTRK,NAMREC,IDREG)
            NAMREC=NAMCEL//'SID'
            CALL LCMGET(IPTRK,NAMREC,IDSUR)
            DO IDIR=1,4
              NAMREC=NAMCEL//'SM'//CDIR(IDIR)
              IF(MESHSP(IDIR) .GT. 0) THEN
                CALL LCMGET(IPTRK,NAMREC,DCMESH(-1,IDIR))
              ENDIF
            ENDDO
*----
*  Track annular mesh
*----
            ITL=2
            IRLCY=NXTLCY(IPRLOC,ITST  ,NDIM  ,MAXMSH,MXSLIN,
     >                   MESHSP,TRKORP,ANGROT,DCMESH,IDIRP ,
     >                   NBCOR(1,ITL)   ,NBSINT(ITL)    ,
     >                   ICSINT(0,1,ITL),DCSINT(1,ITL)  )
*----
*  Track Cartesian mesh
*----
            ITL=1
*     >         +500
            IRLA=NXTLCA(IPRLOC   ,ITST  ,NDIM  ,MAXMSH,MXSLIN,
     >                   MESHSP,TRKROF,ANGROF,DCMESH,
     >                   NBCOR(1,ITL)   ,NBSINT(ITL)    ,
     >                   ICSINT(0,1,ITL),DCSINT(1,ITL)  )
            IF(IRLA .EQ. 2) THEN
              ICOMB=2
              CALL NXTLCU(IPRLOC,MXSLIN,ICOMB ,
     >                    NBCOR ,NBSINT,ICSINT,DCSINT)
              ITL=3
            ELSE
              ITL=2
            ENDIF
*----
*  Insert tracking lines inside the current tracking vector
*----
            CALL NXTLRS(IPRLOC,IOVER ,ITYPG ,LINMAX,MXSLIN,
     >                  NREGP ,NSURP ,MESHSP,IDSUR ,IDREG ,
     >                  NBCOR(1,ITL)   ,NBSINT(ITL)    ,
     >                  ICSINT(0,1,ITL),DCSINT(1,ITL)  ,
     >                  IBLIN ,IELIN ,NUMERO,DLENGT)
            IF(IELIN .GT. LINMAX) THEN
              CALL XABORT(NAMSBR//
     >        ': Tracking vector dimensions too small')
            ENDIF
          ENDIF
        ENDDO
      ENDIF
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      ENDNEW=DLENGT(1)
      IF(IPRLOC .GT. 1000) THEN
        WRITE(IOUT,6020) NUMERO(1),ENDNEW
        DO IKLIN=1,IELIN-1
          IF(NUMERO(IKLIN) .GT. 0) THEN
            ENDNEW=ENDNEW+DLENGT(IKLIN)
            WRITE(IOUT,6021) NUMERO(IKLIN),ENDNEW
          ENDIF
        ENDDO
        WRITE(IOUT,6022) NUMERO(IELIN),DLENGT(IELIN)
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      DEALLOCATE(DCSINT,DRAPIN,DCMESH)
      DEALLOCATE(ICSINT,ITPIN,IDREG,IDSUR)
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT('   Processing PIN = ',I8)
 6020 FORMAT('LinePos ={ {',I10,',',F19.10,'},')
 6021 FORMAT('           {',I10,',',F19.10,'},')
 6022 FORMAT('           {',I10,',',F19.10,'}}')
 6030 FORMAT('   Processing CELL = ',I8,' with turn =',I8)
      END
