*DECK NXTGMD
      SUBROUTINE NXTGMD(IPGEO ,IPTRK ,IPRINT,ITYPBC,ILCELL,NEDIM,
     >                  NBOCEL,NBUCEL,MAXCEL,NUCELL,IUNFLD,
     >                  IEDIMG,NAGGEO,ITURN ,MERGE ,IDIRR ,
     >                  DCMESH,DGMESH)
*
*----------
*
*Purpose:
* To obtain general dimension vector for geometry.
* To evaluate global mesh for assembly.
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
* IPGEO   pointer to the GEOMETRY data structure.
* IPTRK   pointer to the TRACKING data structure.
* IPRINT  print level.
* ITYPBC  type of boundary conditions where:
*         ITYPBC=0 for geometry with all Cartesian boundaries;
*         ITYPBC=1 for geometry with one annular boundary;
*         ITYPBC=2 for geometry with one hexagonal boundary.
* ILCELL  cell level.
* NEDIM   number of elements for general dimension vector.
* NBOCEL  number of cells in original geometry.
* NBUCEL  number of cells in unfolded geometry.
* MAXCEL  maximum number of elements in mesh vector for
*         each directions.
* NUCELL  number of cell after unfolding in
*         $X$, $Y$ and $Z$ directions.
* IUNFLD  array to identify the generating cell (IUNFLD(1,*))
*         and the rotation associated with this region in space.
*
*Parameters: input/output
* IEDIMG  general dimension vector for problem where
* NAGGEO  geometry names.
* ITURN   geometry turns.
* MERGE   geometry merge.
* IDIRR   direction of cell (1 for XYZ, 2 for YZX and 3 for ZXY).
*         Note: for CAR3D without pins IDIRR=1
*         for CAR3D with pins IDIRR specified by pins direction.
* DCMESH  global mesh for each cell.
* DGMESH  meshing vector for global geometry.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*  \\\\
*  This routine is based on the XELDCL routine written by
*  R. Roy for the EXCELT: module. It contains an additional
*  level for cluster subgeometry analysis.
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
      INTEGER          IPRINT,ITYPBC,ILCELL,NEDIM,NBOCEL,NBUCEL,MAXCEL,
     >                 NUCELL(3),IUNFLD(2,NBUCEL)
      INTEGER          IEDIMG(NEDIM),NAGGEO(3,NBOCEL),
     >                 ITURN(NBOCEL),MERGE(NBOCEL),IDIRR(NBOCEL)
      DOUBLE PRECISION DCMESH(3,2,NBOCEL),DGMESH(0:MAXCEL,3)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTGMD')
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
      DOUBLE PRECISION DCUT
      PARAMETER       (DCUT=1.0D-6)
*----
*  Local variables
*----
      INTEGER          ISTATG(NSTATE)
      INTEGER          NDIM,NSDIM,NCDIM,IDIRC
      INTEGER          MAXMSH,MAXREG,NBTCLS,
     >                 MAXPIN,MAXMSP,MAXRSP,NBGCLS,
     >                 NTPIN,NPIN,NBSP,NBRSP,NZU
      INTEGER          ITC,IDIR,ICEL,ICN,ISPL,ICLS,
     >                 IX,IY,IZ,IPOS,INTRN,IRTRN
      INTEGER          ILCMLN,ILCMTY
      CHARACTER        NAMGG*12,NAMREC*12,NAMCL*12
      REAL             OFFCEN(3)
      DOUBLE PRECISION DRW(3),DNW(3)
      REAL             SIDEH
      DOUBLE PRECISION ARGS
      INTEGER          IH,IR,IS,ISS,NCH,NCR,NCS,NSS,ISPLTH
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IGEN,ISPLT
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NAGCEL,NAGCLS
      REAL, ALLOCATABLE, DIMENSION(:) :: RMESH,XMESH
*----
*  Data
*----
      CHARACTER        CDIR(4)*1
      SAVE             CDIR
      DATA             CDIR /'X','Y','Z','R'/
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      NDIM=IEDIMG(1)
      MAXMSH=IEDIMG(16)
      MAXREG=1
      NBTCLS=0
      MAXPIN=0
      MAXMSP=0
      MAXRSP=0
      CALL XDDSET(DGMESH,3*(MAXCEL+1),0.0D0)
      IF(ILCELL .EQ. 1) THEN
        ALLOCATE(IGEN(NBOCEL),NAGCEL(3,NBOCEL))
        CALL LCMGET(IPGEO,'CELL',NAGCEL)
        CALL LCMLEN(IPGEO,'MIX',ILCMLN,ILCMTY)
        IF(ILCMLN .GT. 0) THEN
          CALL LCMGET(IPGEO,'MIX',IGEN)
          DO ICEL=1,NBOCEL
            ICN=-IGEN(ICEL)
            DO ITC=1,3
              NAGGEO(ITC,ICEL)=NAGCEL(ITC,ICN)
            ENDDO
          ENDDO
        ELSE
          CALL LCMGET(IPGEO,'GENERATING',IGEN)
          DO ICEL=1,NBOCEL
            ICN=IGEN(ICEL)
            DO ITC=1,3
              NAGGEO(ITC,ICEL)=NAGCEL(ITC,ICN)
            ENDDO
          ENDDO
        ENDIF
        DEALLOCATE(NAGCEL,IGEN)
        CALL LCMLEN(IPGEO,'TURN',ILCMLN,ILCMTY)
        IF(ILCMLN .GT. 0) THEN
          CALL LCMGET(IPGEO,'TURN',ITURN)
        ELSE
          CALL XDISET(ITURN,NBOCEL,1)
        ENDIF
        CALL LCMLEN(IPGEO,'MERGE',ILCMLN,ILCMTY)
        IF(ILCMLN .GT. 0) THEN
*----
*  Treatment of MERGE not yet implemented
*----
          CALL XABORT(NAMSBR//
     >    ': the module NXT cannot process geometry MERGE option')
*          CALL LCMGET(IPGEO,'MERGE',MERGE)
        ELSE
          DO ICEL=1,NBOCEL
            MERGE(ICEL)=ICEL
          ENDDO
        ENDIF
      ELSE
        NAMGG='            '
        READ(NAMGG,'(3A4)') (NAGGEO(ITC,1),ITC=1,3)
        ITURN(1)=1
        MERGE(1)=1
      ENDIF
      DO ICEL=1,NBOCEL
        IF(ILCELL .EQ. 1) THEN
          WRITE(NAMGG,'(3A4)') (NAGGEO(ITC,ICEL),ITC=1,3)
          CALL LCMSIX(IPGEO,NAMGG,1)
        ENDIF
        CALL XDISET(ISTATG,NSTATE,0)
        CALL LCMGET(IPGEO,'STATE-VECTOR',ISTATG)
        IRTRN=ITURN(ICEL)
*----
*  Find cell dimension and direction
*  Default is CAR3D
*----
        IDIRR(ICEL)=0
        NSDIM=3
        IF(ISTATG(1) .EQ. 5  .OR.
     >     ISTATG(1) .EQ. 8  .OR.
     >     ISTATG(1) .EQ. 12 .OR.
     >     ISTATG(1) .EQ. 20 .OR.
     >     ISTATG(1) .EQ. 26 ) THEN
          IDIRR(ICEL)=1
          NSDIM=2
        ELSE IF(ISTATG(1) .EQ. 21) THEN
          IDIRR(ICEL)=2
        ELSE IF(ISTATG(1) .EQ. 22) THEN
          IDIRR(ICEL)=3
        ELSE IF(ISTATG(1) .EQ. 23 .OR.
     >     ISTATG(1) .EQ. 27 ) THEN
          IDIRR(ICEL)=1
        ENDIF
        IF(NDIM .NE. NSDIM) CALL XABORT(NAMSBR//
     >  ': Geometry and sub-geometry must have the same dimensions')
        IF(ISTATG(1) .EQ. 12 .OR.
     >     ISTATG(1) .EQ. 13 ) THEN
          MAXMSH=MAX(MAXMSH,ISTATG(2),2*ISTATG(3)+1,ISTATG(5))
          MAXREG=MAX(MAXREG,6*ISTATG(3)*ISTATG(3)*MAX(ISTATG(5),1))
          ALLOCATE(ISPLT(MAXMSH))
          MAXMSP=2*ISTATG(3)+1
          NBRSP=6*ISTATG(3)*ISTATG(3)
          NAMREC='SPLITH      '
          CALL LCMLEN(IPGEO,NAMREC,ILCMLN,ILCMTY)
          IF(ILCMLN .GT. 0) THEN
            CALL LCMGET(IPGEO,NAMREC,ISPLTH)
            ISPLTH=ISPLTH*ISTATG(3)
            MAXMSP=MAX(2*ISPLTH+1,MAXMSP)
            NBRSP=MAX(6*ISPLTH*ISPLTH,NBRSP)
          ENDIF
          IDIR=3
          NAMREC='SPLIT'//CDIR(IDIR)//'      '
          CALL LCMLEN(IPGEO,NAMREC,ILCMLN,ILCMTY)
          IF(ILCMLN .GT. 0) THEN
            CALL LCMGET(IPGEO,NAMREC,ISPLT)
            NBSP=0
            DO ISPL=1,ILCMLN
              NBSP=NBSP+ABS(ISPLT(ISPL))
            ENDDO
          ELSE
            NBSP=MAX(ISTATG(5),1)
          ENDIF
          NBRSP=NBRSP*NBSP
          MAXMSP=MAX(MAXMSP,NBSP+1)
          MAXRSP=MAX(MAXRSP,NBRSP)
        ELSE IF(ISTATG(1) .EQ. 26 .OR.
     >          ISTATG(1) .EQ. 27 ) THEN
          MAXMSH=MAX(MAXMSH,ISTATG(2)+1,2*ISTATG(3)+1,ISTATG(5))
          MAXREG=MAX(MAXREG,(ISTATG(2)+1)*(6*ISTATG(3)*ISTATG(3))
     >                      *MAX(ISTATG(5),1))
          ALLOCATE(ISPLT(MAXMSH))
          MAXMSP=2*ISTATG(3)+1
          NBRSP=6*ISTATG(3)*ISTATG(3)
          NAMREC='SPLITH      '
          CALL LCMLEN(IPGEO,NAMREC,ILCMLN,ILCMTY)
          IF(ILCMLN .GT. 0) THEN
            CALL LCMGET(IPGEO,NAMREC,ISPLTH)
            ISPLTH=ISPLTH*ISTATG(3)
            MAXMSP=MAX(2*ISPLTH+1,MAXMSP)
            NBRSP=MAX(6*ISPLTH*ISPLTH,NBRSP)
          ENDIF
          IDIR=3
          NAMREC='SPLIT'//CDIR(IDIR)//'      '
          CALL LCMLEN(IPGEO,NAMREC,ILCMLN,ILCMTY)
          IF(ILCMLN .GT. 0) THEN
            CALL LCMGET(IPGEO,NAMREC,ISPLT)
            NBSP=0
            DO ISPL=1,ILCMLN
              NBSP=NBSP+ABS(ISPLT(ISPL))
            ENDDO
          ELSE
            NBSP=MAX(ISTATG(5),1)
          ENDIF
          NBRSP=(ISTATG(2)+1)*NBRSP*NBSP
          MAXMSP=MAX(MAXMSP,NBSP+1)
          MAXRSP=MAX(MAXRSP,NBRSP)
        ELSE
          MAXMSH=MAX(MAXMSH,ISTATG(2),ISTATG(3),ISTATG(4),ISTATG(5))
          MAXREG=MAX(MAXREG,(ISTATG(2)+1)*
     >           MAX(ISTATG(3),1)*MAX(ISTATG(4),1)*MAX(ISTATG(5),1))
          ALLOCATE(ISPLT(MAXMSH))
          NBRSP=1
          DO IDIR=1,4
            NAMREC='SPLIT'//CDIR(IDIR)//'      '
            CALL LCMLEN(IPGEO,NAMREC,ILCMLN,ILCMTY)
            IF(ILCMLN .GT. 0) THEN
              CALL LCMGET(IPGEO,NAMREC,ISPLT)
              NBSP=0
              DO ISPL=1,ILCMLN
                NBSP=NBSP+ABS(ISPLT(ISPL))
              ENDDO
*            NBRSP=NBRSP*(NBSP+1)
              IF(IDIR .EQ. 4) THEN
                NBRSP=NBRSP*(NBSP+1)
              ELSE
                NBRSP=NBRSP*NBSP
              ENDIF
              MAXMSP=MAX(MAXMSP,NBSP+1)
            ENDIF
          ENDDO
          MAXRSP=MAX(MAXRSP,NBRSP)
        ENDIF
        DEALLOCATE(ISPLT)
*----
*  Get off center
*----
        INTRN=1
        CALL XDRSET(OFFCEN,3,0.0)
        CALL LCMLEN(IPGEO,'OFFCENTER   ',ILCMLN,ILCMTY)
        IF(ILCMLN .GT. 0)
     >  CALL LCMGET(IPGEO,'OFFCENTER   ',OFFCEN)
        DO IDIR=1,3
          DRW(IDIR)=OFFCEN(IDIR)
        ENDDO
        CALL NXTTRM(IRTRN,INTRN,DRW,DCMESH(1,2,ICEL))
        IF(ITYPBC .EQ. 2) THEN
*----
*  Get H mesh
*----
          ALLOCATE(RMESH(MAXMSH+1))
          CALL LCMGET(IPGEO,'SIDE        ',SIDEH)
          DRW(1)=DBLE(SIDEH)
*----
*  Get z-mesh
*----
          IDIR=3
          NAMREC='MESH'//CDIR(IDIR)//'       '
          CALL LCMLEN(IPGEO,NAMREC,ILCMLN,ILCMTY)
          IF(ILCMLN .GT. 0) THEN
            CALL LCMGET(IPGEO,NAMREC,RMESH)
            DRW(IDIR)=DBLE(RMESH(ILCMLN))-DBLE(RMESH(1))
          ENDIF
          DEALLOCATE(RMESH)
          CALL NXTTRM(IRTRN,INTRN,DRW,DCMESH(1,1,ICEL))
*----
*  Fill DGMESH according to IUNFLD
*----
          IPOS=0
          NZU=MAX(NUCELL(3),1)
          NCH=NUCELL(1)
          ARGS=DBLE(12*NCH-3)
          NCR=(NINT(SQRT(ARGS))+3)/6
          DO IZ=1,MAX(NUCELL(3),1)
            IH=0
            DO IR=1,NCR
              NCS=6
              NSS=MAX(1,IR-1)
              IF(IR.EQ.1) NCS=1
              DO IS=1,NCS
                DO ISS=1,NSS
                  IH=IH+1
                  IPOS=IPOS+1
*----
*  Locate crown (IR), sector (IS) and hexagon number in sector (ISS)
*----
                  IF(IUNFLD(1,IPOS) .EQ. ICEL) THEN
*----
*  cell is located at this position
*  Turn width according to symmetry
*----
                    INTRN=IUNFLD(2,IPOS)
                    CALL NXTTRM(IRTRN,INTRN,DRW,DNW)
                    IF(DGMESH(0,1) .EQ. 0.0D0) THEN
                      DGMESH(0,1)=DRW(1)
                      DGMESH(0,2)=DRW(1)
                    ELSE IF(ABS(DGMESH(0,1)-DRW(1)) .GT. DCUT) THEN
                      CALL XABORT(NAMSBR//
     >              ': H- Mesh in assembly is not uniform')
                    ENDIF
                    CALL NXTHCL(IPRINT,IR,IS,ISS,DRW(1),
     >                          DGMESH(IH,1),DGMESH(IH,2))
                    IF(DGMESH(IZ,3) .EQ. 0.0D0) THEN
                      DGMESH(IZ,3)=DNW(3)
                    ELSE IF(ABS(DGMESH(IZ,3)-DNW(3)) .GT. DCUT) THEN
                      CALL XABORT(NAMSBR//
     >                ': Z Mesh in assembly is not uniform')
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE
*----
*  Get Cartesian mesh
*----
          ALLOCATE(XMESH(MAXMSH+1))
          DO IDIR=1,3
            NAMREC='MESH'//CDIR(IDIR)//'       '
            CALL LCMLEN(IPGEO,NAMREC,ILCMLN,ILCMTY)
            IF(ILCMLN .GT. 0) THEN
              CALL LCMGET(IPGEO,NAMREC,XMESH)
              DRW(IDIR)=DBLE(XMESH(ILCMLN))-DBLE(XMESH(1))
            ENDIF
          ENDDO
          DEALLOCATE(XMESH)
          CALL NXTTRM(IRTRN,INTRN,DRW,DCMESH(1,1,ICEL))
*----
*  Fill DGMESH according to IUNFLD
*----
          IPOS=0
          NZU=MAX(NUCELL(3),1)
          DO IZ=1,MAX(NUCELL(3),1)
            DO IY=1,NUCELL(2)
              DO IX=1,NUCELL(1)
                IPOS=IPOS+1
                IF(IUNFLD(1,IPOS) .EQ. ICEL) THEN
*----
*  cell is located at this position
*  Turn width according to symmetry
*----
                  INTRN=IUNFLD(2,IPOS)
                  CALL NXTTRM(IRTRN,INTRN,DRW,DNW)
                  IF(DGMESH(IX,1) .EQ. 0.0D0) THEN
                    DGMESH(IX,1)=DNW(1)
                  ELSE IF(ABS(DGMESH(IX,1)-DNW(1)) .GT. DCUT) THEN
                    CALL XABORT(NAMSBR//
     >              ': X Mesh in assembly is not uniform')
                  ENDIF
                  IF(DGMESH(IY,2) .EQ. 0.0D0) THEN
                    DGMESH(IY,2)=DNW(2)
                  ELSE IF(ABS(DGMESH(IY,2)-DNW(2)) .GT. DCUT) THEN
                    CALL XABORT(NAMSBR//
     >              ': Y Mesh in assembly is not uniform')
                  ENDIF
                  IF(DGMESH(IZ,3) .EQ. 0.0D0) THEN
                    DGMESH(IZ,3)=DNW(3)
                  ELSE IF(ABS(DGMESH(IZ,3)-DNW(3)) .GT. DCUT) THEN
                    CALL XABORT(NAMSBR//
     >              ': Z Mesh in assembly is not uniform')
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF(ISTATG(13) .GT. 0) THEN
          NBGCLS=ISTATG(13)
          NBTCLS=NBTCLS+NBGCLS
          ALLOCATE(NAGCLS(3,NBGCLS))
          CALL LCMGET(IPGEO,'CLUSTER',NAGCLS)
          NTPIN=0
          DO ICLS=1,NBGCLS
            WRITE(NAMCL,'(3A4)') (NAGCLS(ITC,ICLS),ITC=1,3)
            CALL LCMSIX(IPGEO,NAMCL,1)
            CALL XDISET(ISTATG,NSTATE,0)
            CALL LCMGET(IPGEO,'STATE-VECTOR',ISTATG)
*----
*  Find cluster dimension and direction
*----
            IDIRC=1
            NCDIM=3
            IF(ISTATG(1) .EQ. 3) THEN
              NCDIM=2
            ELSE IF(ISTATG(1) .EQ. 10) THEN
              IDIRC=2
            ELSE IF(ISTATG(1) .EQ. 11) THEN
              IDIRC=3
            ELSE IF(ISTATG(1) .EQ. 6) THEN
              IDIRC=1
            ELSE
              CALL XABORT(NAMSBR//': Invalid pin geometry'//
     >' -- Only TUBE* permitted')
            ENDIF
            IF(NSDIM .NE. NCDIM) CALL XABORT(NAMSBR//
     >': Sub-geometry and clusters must have the same dimensions')
            IF(IDIRC .NE. IDIRR(ICEL)) THEN
              IF(IDIRR(ICEL) .NE. 0) CALL XABORT(NAMSBR//
     >': Pin and geometry direction are incompatible')
              IDIRR(ICEL)=IDIRC
            ENDIF
            MAXMSH=MAX(MAXMSH,ISTATG(2),ISTATG(3),ISTATG(4),ISTATG(5))
            MAXREG=MAX(MAXREG,ISTATG(6))
            CALL LCMGET(IPGEO,'NPIN',NPIN)
            NTPIN=NTPIN+NPIN
            ALLOCATE(ISPLT(MAXMSH))
            NBRSP=1
            DO IDIR=1,4
              NAMREC='SPLIT'//CDIR(IDIR)//'      '
              CALL LCMLEN(IPGEO,NAMREC,ILCMLN,ILCMTY)
              IF(ILCMLN .GT. 0) THEN
                CALL LCMGET(IPGEO,NAMREC,ISPLT)
                NBSP=0
                DO ISPL=1,ILCMLN
                  NBSP=NBSP+ABS(ISPLT(ISPL))
                ENDDO
                IF(IDIR .EQ. 4) THEN
                  NBRSP=NBRSP*(NBSP+1)
                ELSE
                  NBRSP=NBRSP*NBSP
                ENDIF
                MAXMSP=MAX(MAXMSP,NBSP+1)
              ENDIF
            ENDDO
            MAXRSP=MAX(MAXRSP,NBRSP)
            DEALLOCATE(ISPLT)
            CALL LCMSIX(IPGEO,NAMCL,2)
          ENDDO
          MAXPIN=MAX(MAXPIN,NTPIN)
*----
*  Now test if pin overlapp
*----
          CALL NXTTPO(IPGEO ,IPRINT,ITYPBC,NBGCLS,NTPIN ,MAXMSH,NCDIM ,
     >                IDIRR(ICEL) ,DRW   ,OFFCEN,NAGCLS)
          DEALLOCATE(NAGCLS)
        ENDIF
        IF(IDIRR(ICEL) .EQ. 0) IDIRR(ICEL)=1
        IF(ILCELL .EQ. 1) THEN
          CALL LCMSIX(IPGEO,NAMGG,2)
        ENDIF
      ENDDO
      IEDIMG(16)=MAXMSH
      IEDIMG(17)=MAXREG
      IEDIMG(18)=NBTCLS
      IEDIMG(19)=MAX(MAXPIN,1)
      IEDIMG(20)=MAXMSP+1
      IEDIMG(21)=MAXRSP
      IF(ITYPBC .EQ. 2) THEN
        IF(NUCELL(1) .GT. 0) THEN
          DO IDIR=1,2
            NAMREC='G00000001SM'//CDIR(IDIR)
            CALL LCMPUT(IPTRK,NAMREC,(NUCELL(1)+1),4,
     >                  DGMESH(0,IDIR))
          ENDDO
        ENDIF
        IDIR=3
        IF(NUCELL(IDIR) .GT. 0) THEN
          DO IX=1,NUCELL(IDIR)
            DGMESH(IX,IDIR)=DGMESH(IX-1,IDIR)+DGMESH(IX,IDIR)
          ENDDO
          NAMREC='G00000001SM'//CDIR(IDIR)
          CALL LCMPUT(IPTRK,NAMREC,(NUCELL(IDIR)+1),4,
     >                DGMESH(0,IDIR))
        ENDIF
      ELSE
        DO IDIR=1,3
          IF(NUCELL(IDIR) .GT. 0) THEN
            DO IX=1,NUCELL(IDIR)
              DGMESH(IX,IDIR)=DGMESH(IX-1,IDIR)+DGMESH(IX,IDIR)
            ENDDO
            NAMREC='G00000001SM'//CDIR(IDIR)
            CALL LCMPUT(IPTRK,NAMREC,(NUCELL(IDIR)+1),4,
     >                  DGMESH(0,IDIR))
          ENDIF
        ENDDO
      ENDIF
      CALL LCMPUT(IPTRK,'G00000001DIM',NEDIM,1,IEDIMG)
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
