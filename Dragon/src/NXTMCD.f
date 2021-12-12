*DECK NXTMCD
      SUBROUTINE NXTMCD(IPGEO ,IPTRK ,IPRINT,NDIM  ,ILCELL,NBOCEL,
     >                  MAXMSH,MAXREG,MAXPIN,NBTCLS,ITSYM ,IDFEX ,
     >                  DCMESH,NAGGEO,ITURN ,IDIRR ,NFSUR ,NFREG ,
     >                  MXGSUR,MXGREG)
*
*----------
*
*Purpose:
* Create a multicell description for the geometry and verify
* if intrinsic cell symmetry required by boundary conditions
* are satisfied.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau.
*
*Parameters: input
* IPGEO   pointer to the GEOMETRY data structure.
* IPTRK   pointer to the TRACKING data structure.
* IPRINT  print level.
* NDIM    problem dimensions.
* ILCELL  cell level.
* NBOCEL  number of cells in original geometry.
* MAXMSH  maximum number of elements in MESH array.
* MAXREG  maximum number of elements in MIX array.
* MAXPIN  maximum number of pins in clusters.
* NBTCLS  total number of cluster geometry.
* ITSYM   array to identify the symmetry to test for each original
*         cell where:
*         ITSYM(1,*) identify $X$ symmetry;
*         ITSYM(2,*) identify $Y$ symmetry;
*         ITSYM(3,*) identify $Z$ symmetry;
*         ITSYM(4,*) identify $X-Y$ symmetry.
*         A value of 0 indicate that the geometry does not need
*         to be verified while a value of 1 implies a verification
*         of the geometry.
* IDFEX   identify faces associated with external boundary for a
*         generating cell and number of times this cell is used. Here:
*         IDFEX( 1,*)  identify bottom $U$ hexagonal face;
*         IDFEX( 2,*)  identify top $U$ hexagonal face;
*         IDFEX( 3,*)  identify bottom $V$ hexagonal face;
*         IDFEX( 4,*)  identify top $V$ hexagonal face;
*         IDFEX( 5,*)  identify bottom $Z$ face;
*         IDFEX( 6,*)  identify top $Z$ face;
*         IDFEX( 7,*)  not used;
*         IDFEX( 8,*)  not used;
*         IDFEX( 9,*)  identify bottom $W$ hexagonal face;
*         IDFEX(10,*) identify top $W$ hexagonal face.
* DCMESH  global mesh for each cell.
* NAGGEO  names of generating geometries.
* ITURN   geometry turns.
* IDIRR   direction of cell (1 for XYZ, 2 for YZX and 3 for ZXY).
*         Note: for CAR3D without pins IDIRR=1
*         for CAR3D with pins IDIRR specified by pins direction.
*
*Parameters: output
* NFSUR   final number of surfaces.
* NFREG   final number of regions.
* MXGSUR  maximum number of surfaces for any geometry.
* MXGREG  maximum number of region for any geometry.
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
      TYPE(C_PTR)      IPGEO,IPTRK
*      INTEGER          IPGEO,IPTRK
      INTEGER          IPRINT,NDIM,
     >                 ILCELL,NBOCEL,
     >                 MAXMSH,MAXREG,MAXPIN,NBTCLS,
     >                 ITSYM(4,NBOCEL),IDFEX(0:10,NBOCEL)
      DOUBLE PRECISION DCMESH(3,2,NBOCEL)
      INTEGER          NAGGEO(3,NBOCEL),ITURN(NBOCEL),
     >                 IDIRR(NBOCEL),NFSUR,NFREG,MXGSUR,MXGREG
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTMCD')
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
      DOUBLE PRECISION DCUTOF,DZERO,DONE,DTWO
      PARAMETER       (DCUTOF=1.0D-8,DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI
*----
*  Local variables
*----
      INTEGER          IOFNAC,ICEL,ITC,ITYPG,ITRN,IDIRC,ITST,MAXMSS,
     >                 NBGCLS,ICLS,ITYPC,IGCLS,ILEV,
     >                 NMIX,NMIXS,NMIXC,NMIXCS,
     >                 NREG,NSUR,NREGS,NSURS,NREGC,NSURC,NREGCS,NSURCS,
     >                 NREGN,NSURN,NREGCN,NSURCN
      INTEGER          ISTATG(NSTATE),ISTATC(NSTATE),
     >                 IEDIMX(NSTATE),IEDIMP(NSTATE)
      INTEGER          NM(4),NMC(4),NMS(4),NMCS(4),
     >                 ITSYMC(4),IGSYM,IDSYM
      CHARACTER        NAMGG*12,NAMCL*12,NAMREC*12
      INTEGER          ILCMLN,ILCMTY
      INTEGER          NAREG,NASUR,NEREN
      INTEGER          NPIN,IP,IPINT,ICPIN,NPIR,IDGPP,NBPIN
      DOUBLE PRECISION DELTA,TWOPI,DROUT,DHPIN
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ITPIN,ISPLT,MIX
      REAL, ALLOCATABLE, DIMENSION(:) :: RMESH,RPIN,APIN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DAMESH,DRAPIN
*
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: MIXC,ISPLTR
      INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: NAGCLS,INDXSR
      INTEGER, ALLOCATABLE, DIMENSION(:)     :: IDSUR,IDREG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DAMESR
*----
*  Data
*----
      CHARACTER        CLEV(2)*1
      SAVE             CLEV
      DATA             CLEV /'C','P'/
*----
*  Scratch storage allocation
*   ITPIN   temporary array for reading pin type.
*   ISPLT   temporary array for reading cell mesh split.
*   MIX     temporary array for reading cell mixture (including HMIX).
*   RMESH   temporary vector for reading cell mesh array.
*   RPIN    temporary vector for reading pin positions.
*   APIN    temporary vector for reading pin angles.
*   DAMESH  temporary vector for storing global mesh array.
*   DRAPIN  temporary vector for storing global pin positions.
*----
      ALLOCATE(ITPIN(3,MAXPIN),ISPLT(MAXMSH,4),MIX(MAXREG,2))
      ALLOCATE(RMESH(0:MAXMSH),RPIN(MAXPIN),APIN(MAXPIN))
      ALLOCATE(DAMESH(-1:MAXMSH,4),DRAPIN(-1:4,MAXPIN))
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      PI=XDRCST('Pi',' ')
      TWOPI=DTWO*PI
      ALLOCATE(MIXC(MAXREG,2,2),NAGCLS(3,NBTCLS),ISPLTR(MAXMSH,4,2))
      MXGSUR=0
      MXGREG=0
      IOFNAC=0
      IGCLS=0
      ICPIN=0
      NAREG=0
      NASUR=0
      DO ICEL=1,NBOCEL
        IF(ILCELL .EQ. 1) THEN
          WRITE(NAMGG,'(3A4)') (NAGGEO(ITC,ICEL),ITC=1,3)
          CALL LCMSIX(IPGEO,NAMGG,1)
        ELSE
          NAMGG='/           '
        ENDIF
        IF(IPRINT .GE. 10) THEN
          WRITE(IOUT,6010) ICEL,NAMGG
        ENDIF
        CALL XDISET(ISTATG,NSTATE,0)
        CALL LCMGET(IPGEO,'STATE-VECTOR',ISTATG)
        ITYPG=ISTATG(1)
        IF(ITYPG .EQ. 20) THEN
          IF(ISTATG(2) .EQ. 0) ITYPG=5
        ELSE IF(ITYPG .EQ. 21 .OR. ITYPG .EQ. 22 .OR.
     >          ITYPG .EQ. 23) THEN
          IF(ISTATG(2) .EQ. 0) ITYPG=7
        ELSE IF(ITYPG .EQ. 26) THEN
          IF(ISTATG(2) .EQ. 0) ITYPG=12
        ELSE IF(ITYPG .EQ. 27) THEN
          IF(ISTATG(2) .EQ. 0) ITYPG=13
        ENDIF
        NM(4)=ISTATG(2)
        NM(1)=ISTATG(3)
        NM(2)=ISTATG(4)
        NM(3)=ISTATG(5)
        NMIX=ISTATG(6)
        ITRN=ITURN(ICEL)
        IDIRC=IDIRR(ICEL)
        NBGCLS=ISTATG(13)
        ILEV=1
        CALL XDISET(ITPIN,3*MAXPIN,0)
        ICPIN=0
        IGSYM=0
        DO IDSYM=1,4
          ITSYMC(IDSYM)=0
          IGSYM=IGSYM+ABS(ITSYM(IDSYM,ICEL))
        ENDDO
        IF(NBGCLS .NE. 0) THEN
*----
*  Get pin cluster geometry information
*----
          CALL XDDSET(DRAPIN,6*MAXPIN,DZERO)
          ILEV=2
          CALL LCMGET(IPGEO,'CLUSTER',NAGCLS(1,IOFNAC+1))
          DO ICLS=1,NBGCLS
            IGCLS=IGCLS+1
            WRITE(NAMCL,'(3A4)') (NAGCLS(ITC,IOFNAC+ICLS),ITC=1,3)
            CALL LCMSIX(IPGEO,NAMCL,1)
            IF(IPRINT .GE. 10) THEN
              WRITE(IOUT,6011) ICLS,NAMCL
            ENDIF
            CALL XDISET(ISTATC,NSTATE,0)
            CALL LCMGET(IPGEO,'STATE-VECTOR',ISTATC)
            ITYPC=ISTATC(1)
            IF(ITYPC .EQ. 20) THEN
              IF(ISTATC(2) .EQ. 0) ITYPC=5
            ELSE IF(ITYPC .EQ. 21 .OR. ITYPC .EQ. 22 .OR.
     >              ITYPC .EQ. 23) THEN
              IF(ISTATC(2) .EQ. 0) ITYPC=7
            ENDIF
*----
*  Cartesian or Annular pin
*----
            IPINT=0
            IF(ITYPC .EQ.  5 .OR. ITYPC .EQ.  7 .OR.
     >         ITYPC .EQ. 20 .OR. ITYPC .EQ. 23      ) THEN
              IPINT=-3
            ELSE IF(ITYPC .EQ.  3 .OR. ITYPC .EQ.  6 ) THEN
              IPINT=3
            ELSE IF(ITYPC .EQ. 21) THEN
              IPINT=-1
            ELSE IF(ITYPC .EQ. 10) THEN
              IPINT=1
            ELSE IF(ITYPC .EQ. 22) THEN
              IPINT=-2
            ELSE IF(ITYPC .EQ. 11 ) THEN
              IPINT=2
            ENDIF
            NMC(4)=ISTATC(2)
            NMC(1)=ISTATC(3)
            NMC(2)=ISTATC(4)
            NMC(3)=ISTATC(5)
            NMIXC=ISTATC(6)
            CALL XDDSET(DAMESH,(MAXMSH+2)*4,DZERO)
            CALL XDISET(ISPLT,MAXMSH*4,1)
            CALL XDISET(MIX,MAXREG*2,0)
            CALL XDISET(MIXC,MAXREG*2*2,0)
            CALL NXTEGI(IPGEO ,IPRINT,ITYPC ,MAXMSH,NMIXC ,NMC   ,
     >                  MAXMSS,NMCS  ,NREGC ,NREGCS,NSURC ,NSURCS,
     >                  MIX   ,ISPLT ,DAMESH,RMESH ,MIXC  )
            ALLOCATE(IDSUR(NSURCS),IDREG(NREGCS))
            CALL XDISET(IDREG,NREGCS,0)
            CALL XDISET(IDSUR,NSURCS,0)
            NEREN=NREGCS+NSURCS
            ALLOCATE(DAMESR((MAXMSS+2),4,2))
            CALL XDDSET(DAMESR,(MAXMSS+2)*4*2,DZERO)
            NMIXCS=NREGCS
            DROUT=DAMESH(NMC(4),4)
            IF(ABS(IPINT) .EQ. 3) THEN
*----
*  For pins rotate symmetry by pi/2 around z axis
*----
              IF(ITSYM(1,ICEL) .NE. 0) THEN
                ITSYMC(2)=ITSYM(1,ICEL)
              ENDIF
              IF(ITSYM(2,ICEL) .NE. 0) THEN
                ITSYMC(1)=ITSYM(2,ICEL)
              ENDIF
              ITSYMC(3)=ITSYM(3,ICEL)
              ITSYMC(4)=4*ITSYM(4,ICEL)
            ELSE IF(ABS(IPINT) .EQ. 2) THEN
*----
*  For pins rotate symmetry by pi/2  around y axis
*----
              IF(ITSYM(1,ICEL) .NE. 0) THEN
                ITSYMC(3)=ITSYM(1,ICEL)
              ENDIF
              IF(ITSYM(3,ICEL) .NE. 0) THEN
                ITSYMC(1)=ITSYM(3,ICEL)
              ENDIF
              ITSYMC(2)=ITSYM(2,ICEL)
              ITSYMC(4)=4*ITSYM(4,ICEL)
            ELSE IF(ABS(IPINT) .EQ. 1) THEN
*----
*  For pins rotate symmetry by pi/2  around x axis
*----
              IF(ITSYM(3,ICEL) .NE. 0) THEN
                ITSYMC(2)=ITSYM(3,ICEL)
              ENDIF
              IF(ITSYM(2,ICEL) .NE. 0) THEN
                ITSYMC(3)=ITSYM(2,ICEL)
              ENDIF
              ITSYMC(1)=ITSYM(1,ICEL)
              ITSYMC(4)=4*ITSYM(4,ICEL)
            ENDIF
*----
*  Test if pin cells satisfy symmetry options
*----
            IF(IGSYM .GT. 0) THEN
*----
*  Rotate pin cells and test for intrinsic symmetry
*----
              ITST=1
              CALL NXTRIS(IPRINT,ITYPC ,MAXMSH,NREGC ,ITRN  ,ITST  ,
     >                    ITSYMC,NMC   ,MIX   ,ISPLT ,DAMESH,
     >                    NMCS  ,MIXC  ,ISPLTR,
     >                    DAMESR)
            ELSE
*----
*  Rotate pin cells
*----
              ITST=0
              CALL NXTRIS(IPRINT,ITYPC ,MAXMSH,NREGC ,ITRN  ,ITST  ,
     >                    ITSYMC,NMC   ,MIX   ,ISPLT ,DAMESH,
     >                    NMCS  ,MIXC  ,ISPLTR,DAMESR)
            ENDIF
*----
*  Final mesh after unfolding
*----
            IF(MOD(ITRN,2) .EQ. 0) THEN
*----
*  First turn geometry type if required
*----
              IF(ITYPC .EQ. 10) THEN
                ITYPC=11
              ELSE IF(ITYPC .EQ.11) THEN
                ITYPC=10
              ELSE IF(ITYPC .EQ.21) THEN
                ITYPC=22
              ELSE IF(ITYPC .EQ.22) THEN
                ITYPC=21
              ENDIF
            ENDIF
*----
*  Then create new mesh
*----
            CALL NXTSGI(IPTRK ,IPRINT,MAXMSH,ITYPC ,IGCLS ,ILEV  ,
     >                  MAXMSS,NMIXC ,NMC   ,MIX   ,DAMESH,ISPLT ,
     >                  NMIXCS,NMCS  ,DAMESR,ITSYMC,
     >                  NREGCS,NSURCS,NREGCN,NSURCN,NEREN ,
     >                  IDREG ,IDSUR )
            ALLOCATE(INDXSR(5,(NREGCS+NSURCS+1)))
            CALL XDISET(INDXSR,5*(NREGCS+NSURCS+1),0)
            NBPIN=0
            MXGSUR=MAX(MXGSUR,NSURCS)
            MXGREG=MAX(MXGREG,NREGCS)
            CALL XDISET(IEDIMP,NSTATE,0)
            IEDIMP(1)=ITYPC
            IEDIMP(2)=NMCS(4)
            IEDIMP(3)=NMCS(1)
            IEDIMP(4)=NMCS(2)
            IEDIMP(5)=NMCS(3)
            IEDIMP(6)=NMIXCS
            IEDIMP(7)=ILEV
            IEDIMP(8)=NREGCS
            IEDIMP(9)=NSURCS
            IEDIMP(10)=NREGCN
            IEDIMP(11)=NSURCN
            IEDIMP(12)=NAREG+1
            CALL NXTVOL(IPTRK ,IPRINT,MAXMSS,ITYPC ,IDIRC ,IGCLS ,
     >                  ILEV  ,NMCS  ,NREGCS,NSURCS,NREGCN,NSURCN,
     >                  MAXPIN,NBPIN ,ITPIN ,DRAPIN, IDREG,IDSUR ,
     >                  DAMESR,INDXSR,NAREG )
            IEDIMP(13)=NAREG
            DEALLOCATE(DAMESR)
*----
*  Read pin locations
*----
            CALL LCMGET(IPGEO,'NPIN',NPIN)
            CALL LCMLEN(IPGEO,'RPIN',ILCMLN,ILCMTY)
            IF(ILCMLN .EQ. 0) THEN
               CALL XABORT(NAMSBR//
     >        ': RPIN is absent for pin cluster')
            ELSE IF(ILCMLN .EQ. 1) THEN
              CALL LCMGET(IPGEO,'RPIN',RPIN)
              DO IP=1,NPIN
                DRAPIN(0,ICPIN+IP)=DBLE(RPIN(1))
              ENDDO
            ELSE IF(ILCMLN .EQ. NPIN) THEN
              CALL LCMGET(IPGEO,'RPIN',RPIN)
              DO IP=1,NPIN
                DRAPIN(0,ICPIN+IP)=DBLE(RPIN(IP))
              ENDDO
            ELSE
              CALL XABORT(NAMSBR//
     >        ': Length of RPIN vector is invalid')
            ENDIF
            CALL LCMLEN(IPGEO,'APIN',ILCMLN,ILCMTY)
            IF(ILCMLN .EQ. 0) THEN
              DRAPIN(-1,ICPIN+1)=DZERO
              DELTA=(DTWO*PI)/DBLE(NPIN)
              DO IP=2,NPIN
                DRAPIN(-1,ICPIN+IP)=DRAPIN(-1,ICPIN+IP-1)+DELTA
              ENDDO
            ELSE IF(ILCMLN .EQ. 1) THEN
              CALL LCMGET(IPGEO,'APIN',APIN)
              DRAPIN(-1,ICPIN+1)=DBLE(APIN(1))
              DELTA=(DTWO*PI)/DBLE(NPIN)
              DO IP=2,NPIN
                DRAPIN(-1,ICPIN+IP)=DRAPIN(-1,ICPIN+IP-1)+DELTA
              ENDDO
            ELSE IF(ILCMLN .EQ. NPIN) THEN
              CALL LCMGET(IPGEO,'APIN',APIN)
              DO IP=1,NPIN
                DRAPIN(-1,ICPIN+IP)=DBLE(APIN(IP))
              ENDDO
            ELSE
              CALL XABORT(NAMSBR//
     >        ': Length of APIN vector is invalid')
            ENDIF
*----
*  Find pin height
*----
            IDGPP=MOD(IDIRC+1,3)+1
            IF(NDIM .EQ. 3) THEN
              DHPIN=DAMESH(NMC(IDGPP),IDGPP)-DAMESH(0,IDGPP)
            ELSE
              DHPIN=DONE
            ENDIF
*----
*  Replace angles in range 0 to 2*Pi
*----
            DO IP=1,NPIN
              IF(DRAPIN(-1,ICPIN+IP) .GE. DZERO) THEN
                NPIR=INT((DRAPIN(-1,ICPIN+IP)+DCUTOF)/TWOPI)
                DRAPIN(-1,ICPIN+IP)=DRAPIN(-1,ICPIN+IP)
     >                             -DBLE(NPIR)*TWOPI
              ELSE
                NPIR=INT((DRAPIN(-1,ICPIN+IP)-DCUTOF)/TWOPI)
                DRAPIN(-1,ICPIN+IP)=DRAPIN(-1,ICPIN+IP)
     >                             -DBLE(NPIR-1)*TWOPI
              ENDIF
              DRAPIN(IDGPP,ICPIN+IP)=DHPIN
              DRAPIN(4,ICPIN+IP)=DROUT
              ITPIN(1,ICPIN+IP)=ICEL
              ITPIN(2,ICPIN+IP)=IGCLS
              ITPIN(3,ICPIN+IP)=IPINT
            ENDDO
            IF(IGSYM .GT. 0) THEN
*----
*  Test pin in cluster for symmetry
*----
              CALL NXTTPS(IPRINT,NPIN  ,IDGPP ,ITSYM(1,ICEL),
     >                    DRAPIN(-1,ICPIN+1))
            ENDIF
*----
*  Renumber pin surfaces
*----
            IEDIMP(14)=NASUR+1
            CALL NXTRPS(IPTRK ,IPRINT,ITYPC ,IGCLS ,ILEV  ,
     >                  NREGCS,NSURCS,NSURCN,IDFEX(0,ICEL),
     >                  INDXSR,DHPIN ,DCMESH(1,1,ICEL),
     >                  NASUR ,IDSUR )
            DEALLOCATE(INDXSR)
            DEALLOCATE(IDSUR,IDREG)
            IEDIMP(15)=NASUR
            IEDIMP(16)=ICPIN+1
            IEDIMP(17)=NPIN
            IEDIMP(20)=IDIRC
            WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGCLS,'DIM'
            CALL LCMPUT(IPTRK,NAMREC,NSTATE,1,IEDIMP)
            ICPIN=ICPIN+NPIN
            CALL LCMSIX(IPGEO,NAMCL,2)
          ENDDO
*----
*  Save cluster information for this cell
*----
          ILEV=1
          IOFNAC=IOFNAC+NBGCLS
          WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),ICEL,'PIN'
          CALL LCMPUT(IPTRK,NAMREC,6*ICPIN,4,DRAPIN)
          WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),ICEL,'PNT'
          CALL LCMPUT(IPTRK,NAMREC,3*ICPIN,1,ITPIN)
        ENDIF
        CALL XDDSET(DAMESH,(MAXMSH+2)*4,DZERO)
        CALL XDISET(ISPLT,4*MAXMSH,1)
        CALL XDISET(MIX,MAXREG*2,0)
        CALL XDISET(MIXC,2*MAXREG*2,0)
        IF(IPRINT .GE. 10) THEN
          WRITE(IOUT,6012)
        ENDIF
*----
*  Get main geometry information
*----
        CALL NXTEGI(IPGEO ,IPRINT,ITYPG ,MAXMSH,NMIX  ,NM    ,
     >              MAXMSS,NMS   ,NREG  ,NREGS ,NSUR  ,NSURS ,
     >              MIX   ,ISPLT ,DAMESH,
     >              RMESH ,MIXC  )
        NMIXS=NREGS
        ALLOCATE(IDREG(NREGS),IDSUR(NSURS))
        CALL XDISET(IDREG,NREGS,0)
        CALL XDISET(IDSUR,NSURS,0)
        NEREN=NREGS+NSURS
        ALLOCATE(DAMESR((MAXMSS+2),4,2))
        CALL XDDSET(DAMESR,(MAXMSS+2)*4*2,DZERO)
        IF(ITYPG .EQ. 12 .OR. ITYPG .EQ. 13 .OR. 
     >     ITYPG .EQ. 26 .OR. ITYPG .EQ. 27) THEN
          IF(IGSYM .GT. 0) THEN
*----
*  Rotate global geometry and test for intrinsic symmetry
*----
            ITST=1
            CALL NXTRTS(IPRINT,ITYPG ,MAXMSH,NREG  ,ITRN  ,ITST,
     >                  ITSYM(1,ICEL),NM    ,MIX   ,ISPLT ,DAMESH,
     >                  NMS   ,MIXC  ,ISPLTR,DAMESR)
          ELSE
*----
*  Rotate global geometry
*----
            ITST=0
            CALL NXTRTS(IPRINT,ITYPG ,MAXMSH,NREG  ,ITRN  ,ITST,
     >                  ITSYM(1,ICEL),NM    ,MIX   ,ISPLT ,DAMESH,
     >                  NMS   ,MIXC  ,ISPLTR,DAMESR)
          ENDIF
        ELSE
          IF(IGSYM .GT. 0) THEN
*----
*  Rotate global geometry and test for intrinsic symmetry
*----
            ITST=1
            CALL NXTRIS(IPRINT,ITYPG ,MAXMSH,NREG  ,ITRN  ,ITST,
     >                  ITSYM(1,ICEL),NM    ,MIX   ,ISPLT ,DAMESH,
     >                  NMS   ,MIXC  ,ISPLTR,DAMESR)
          ELSE
*----
*  Rotate global geometry
*----
            ITST=0
            CALL NXTRIS(IPRINT,ITYPG ,MAXMSH,NREG  ,ITRN  ,ITST,
     >                  ITSYM(1,ICEL),NM    ,MIX   ,ISPLT ,DAMESH,
     >                  NMS   ,MIXC  ,ISPLTR,DAMESR)
          ENDIF
        ENDIF
*----
*  Final mesh after unfolding
*----
        IF(MOD(ITRN,2) .EQ. 0) THEN
*----
*  First turn geometry type if required
*----
          IF(ITYPG .EQ. 10) THEN
            ITYPG=11
          ELSE IF(ITYPG .EQ.11) THEN
            ITYPG=10
          ELSE IF(ITYPG .EQ.21) THEN
            ITYPG=22
          ELSE IF(ITYPG .EQ.22) THEN
            ITYPG=21
          ENDIF
        ENDIF
*----
*  Then create new mesh
*----
        IF(ITYPG .EQ. 12 .OR. ITYPG .EQ. 13 .OR. 
     >     ITYPG .EQ. 26 .OR. ITYPG .EQ. 27) THEN
          CALL NXTSGT(IPTRK ,IPRINT,MAXMSH,ITYPG ,ICEL  ,ILEV  ,
     >                MAXMSS,NMIX  ,NM    ,MIX   ,DAMESH,ISPLT ,
     >                NMIXS ,NMS   ,DAMESR,
     >                ITSYM(1,ICEL),
     >                NREGS ,NSURS ,NREGN ,NSURN ,NEREN ,
     >                IDREG ,IDSUR )
        ELSE
          CALL NXTSGI(IPTRK ,IPRINT,MAXMSH,ITYPG ,ICEL  ,ILEV  ,
     >                MAXMSS,NMIX  ,NM    ,MIX   ,DAMESH,ISPLT ,
     >                NMIXS ,NMS   ,DAMESR,
     >                ITSYM(1,ICEL),
     >                NREGS ,NSURS ,NREGN ,NSURN ,NEREN ,
     >                IDREG ,IDSUR )
        ENDIF
        ALLOCATE(INDXSR(5,(NREGS+NSURS+1)))
        CALL XDISET(INDXSR,5*(NREGS+NSURS+1),0)
        NBPIN=ICPIN
        MXGSUR=MAX(MXGSUR,NSURS)
        MXGREG=MAX(MXGREG,NREGS)
        CALL XDISET(IEDIMX,NSTATE,0)
        IEDIMX(1)=ITYPG
        IEDIMX(2)=NMS(4)
        IEDIMX(3)=NMS(1)
        IEDIMX(4)=NMS(2)
        IEDIMX(5)=NMS(3)
        IEDIMX(6)=NMIXS
        IEDIMX(7)=ILEV
        IEDIMX(8)=NREGS
        IEDIMX(9)=NSURS
        IEDIMX(10)=NREGN
        IEDIMX(11)=NSURN
        IEDIMX(12)=NAREG+1
        CALL NXTVOL(IPTRK ,IPRINT,MAXMSS,ITYPG ,IDIRC ,ICEL  ,
     >              ILEV  ,NMS   ,NREGS ,NSURS ,NREGN ,NSURN ,
     >              MAXPIN,NBPIN ,ITPIN ,DRAPIN,IDREG ,IDSUR ,
     >              DAMESR,INDXSR,NAREG )
*----
*  Renumber cell surfaces
*----
        IEDIMX(13)=NAREG
        IEDIMX(14)=NASUR+1
        CALL NXTRCS(IPTRK ,IPRINT,ICEL ,ILEV  ,
     >              NREGS ,NSURS ,NSURN ,IDFEX(0,ICEL),
     >              INDXSR,NASUR ,IDSUR )
        IEDIMX(15)=NASUR
        IEDIMX(16)=NBGCLS
        IEDIMX(17)=IGCLS-NBGCLS+1
        IEDIMX(18)=NBPIN
        IEDIMX(19)=IDFEX(0,ICEL)
        IEDIMX(20)=IDIRC
        DEALLOCATE(INDXSR)
        WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),ICEL,'DIM'
        CALL LCMPUT(IPTRK,NAMREC,NSTATE,1,IEDIMX)
        DEALLOCATE(DAMESR)
        DEALLOCATE(IDSUR,IDREG)
        IF(ILCELL .EQ. 1) THEN
          CALL LCMSIX(IPGEO,NAMGG,2)
        ENDIF
      ENDDO
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      DEALLOCATE(ISPLTR,NAGCLS,MIXC)
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      NFSUR=NASUR
      NFREG=NAREG
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(DRAPIN,DAMESH)
      DEALLOCATE(APIN,RPIN,RMESH)
      DEALLOCATE(MIX,ISPLT)
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(' Processing CELL    = ',I5.4,5X,A12)
 6011 FORMAT(' Processing CLUSTER = ',I5.4,5X,A12)
 6012 FORMAT(' Global cell processing ')
      END
