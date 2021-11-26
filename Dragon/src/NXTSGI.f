*DECK NXTSGI
      SUBROUTINE NXTSGI(IPTRK ,IPRINT,MAXMSH,ITYPG ,IGEO  ,ILEV  ,
     >                  MAXMSS,NMIX  ,NM    ,MIX   ,DAMESH,ISPLT ,
     >                  NMIXS ,NMS   ,DAMESS,
     >                  ITSYM ,NREGS ,NSURS ,NREGN ,NSURN ,NEREN ,
     >                  IDREG ,IDSUR )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Discretize geometry according to splitting options.
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
* MAXMSH  maximum number of elements in MESH array.
* ITYPG   type of geometry.
* IGEO    geometry number.
* ILEV    geometry level.
* MAXMSS  maximum number of elements in MESH array after split.
* NMIX    number of elements in MIX array.
* NM      mesh size in all directions ($X$, $Y$, $Z$ and $R$).
* MIX     final mixture description for geometry (including MMIX).
* DAMESH  final mesh description for geometry.
* NMIXS   number of regional mixtures.
* ISPLT   final split desctiption for geometry.
* NMS     mesh size after splitting.
* ITSYM   flag for symmetries to test.
* NREGS   maximum number of regions in splitted geometry.
* NSURS   maximum number of surfaces in splitted geometry.
* NREGN   number of regions in splitted geometry after symmetry.
* NSURN   number of surfaces in splitted geometry after symmetry.
* NEREN   maximum number of elements in IREN.
*
*Parameters: input/output
* DAMESS  mesh description for rotated geometry.
* IDREG   region identifier after symmetry.
* IDSUR   surface identifier after symmetry.
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
      INTEGER          IPRINT,MAXMSH,ITYPG,IGEO,ILEV,MAXMSS,NMIX,NM(4)
      INTEGER          MIX(NMIX,2)
      DOUBLE PRECISION DAMESH(-1:MAXMSH,4)
      INTEGER          ISPLT(MAXMSH,4)
      INTEGER          NMIXS,NMS(4)
      DOUBLE PRECISION DAMESS(-1:MAXMSS,4,2)
      INTEGER          ITSYM(4),NREGS,NSURS,NREGN,NSURN,NEREN
      INTEGER          IDREG(NREGS),IDSUR(NSURS)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTSGI')
      DOUBLE PRECISION DZERO,DONE,DDC,DCAOF
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DDC=1.0D-3,
     >                 DCAOF=DONE+DDC)
*----
*  Local variables
*----
      INTEGER          IDIR,NMTMP,IGEN,IGENS,NBS,IS,
     >                 NX,NY,NZ,NR,IX,IY,IZ,IR,
     >                 NXS,NYS,NZS,NRS,IXS,IYS,IZS,IRS,
     >                 NXX,NYY,NZZ,NRR,IXX,IYY,IZZ,IRR,
     >                 NRX,NRY,NRZ,NRXS,NRYS,NRZS,NZSR,NYSR,NXSR,
     >                 NSXS,NSYS,NSZS,IMIX,IMIXS,ISBOT,ISTOP
      INTEGER          IOFZ,IOFYZ,IOFXYZ,IOSZ,IOSYZ,IOSXYZ,
     >                 IODZ,IODYZ,IODXYZ,IORZ,IORYZ,IORXYZ
      INTEGER          IRO,IRN,IRT,IDS,IDV,NDIM,IPLOC
      DOUBLE PRECISION DDD,DDI,DDO
      CHARACTER        NAMREC*12
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IREN
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: MIXS
*----
*  Data
*----
      CHARACTER        CDIR(1:4)*1,CLEV(2)*1
      SAVE             CDIR,CLEV
      DATA             CDIR /'X','Y','Z','R'/
      DATA             CLEV /'C','P'/
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      ALLOCATE(IREN(NEREN),MIXS(NMIXS,2))
      CALL XDISET(MIXS,NMIXS*2,0)
      IPLOC=IPRINT
      IF(IPLOC .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010) (NM(IDIR),IDIR=1,4)
        DO IDIR=1,4
          NMTMP=NM(IDIR)
          IF(NMTMP .GT. 0) THEN
            WRITE(IOUT,6011) 'MESH'//CDIR(IDIR)//' ='
            WRITE(IOUT,6012) (DAMESH(IX,IDIR),IX=-1,NMTMP)
            WRITE(IOUT,6011) 'SPLT'//CDIR(IDIR)//' ='
            WRITE(IOUT,6013) (ISPLT(IX,IDIR),IX=1,NMTMP)
          ENDIF
        ENDDO
      ENDIF
      NR=0
      NRS=0
*----
*  Cartesian mesh
*----
      DO IDIR=1,3
        DAMESS(-1,IDIR,1)=DAMESH(-1,IDIR)
        DAMESS(0,IDIR,1)=DAMESH(0,IDIR)
        IGENS=0
        DO IGEN=1,NM(IDIR)
          NBS=ABS(ISPLT(IGEN,IDIR))
          DDD=(DAMESH(IGEN,IDIR)-DAMESH(IGEN-1,IDIR))/DBLE(NBS)
          DO IS=1,NBS
            IGENS=IGENS+1
            DAMESS(IGENS,IDIR,1)=DAMESS(IGENS-1,IDIR,1)+DDD
          ENDDO
        ENDDO
      ENDDO
*----
* Radial mesh
*----
      IDIR=4
      DAMESS(-1,IDIR,1)=DAMESH(-1,IDIR)
      DAMESS(0,IDIR,1)=DAMESH(0,IDIR)
      IGENS=0
      DO IGEN=1,NM(IDIR)
        NBS=ISPLT(IGEN,IDIR)
        IF(NBS .LT. 0) THEN
          NBS=-NBS
          DDI=DAMESH(IGEN-1,IDIR)*DAMESH(IGEN-1,IDIR)
          DDO=DAMESH(IGEN,IDIR)*DAMESH(IGEN,IDIR)
          DDD=(DDO-DDI)/DBLE(NBS)
          DO IS=1,NBS
            IGENS=IGENS+1
            DDO=DDI+DDD
            DAMESS(IGENS,IDIR,1)=SQRT(DDO)
            DDI=DDO
          ENDDO
        ELSE
          DDD=(DAMESH(IGEN,IDIR)-DAMESH(IGEN-1,IDIR))/DBLE(NBS)
          DO IS=1,NBS
            IGENS=IGENS+1
            DAMESS(IGENS,IDIR,1)=DAMESS(IGENS-1,IDIR,1)+DDD
          ENDDO
        ENDIF
      ENDDO
      IF(IPLOC .GE. 100) THEN
        WRITE(IOUT,6020) (NMS(IDIR),IDIR=1,4)
        DO IDIR=1,4
          NMTMP=NMS(IDIR)
          IF(NMTMP .GT. 0) THEN
            WRITE(IOUT,6011) 'MESH'//CDIR(IDIR)//' ='
            WRITE(IOUT,6012) (DAMESS(IX,IDIR,1),IX=-1,NMTMP)
          ENDIF
        ENDDO
      ENDIF
      NX=MAX(1,NM(1))
      NXS=MAX(1,NMS(1))
      NY=MAX(1,NM(2))
      NYS=MAX(1,NMS(2))
      NZ=MAX(1,NM(3))
      NZS=MAX(1,NMS(3))
      NRX=1
      NRXS=1
      NRY=1
      NRYS=1
      NRZ=1
      NRZS=1
      NSXS=NXS
      NSYS=NYS
      NSZS=NZS
      NDIM=3
      IF(ITYPG .EQ. 3 .OR. ITYPG .EQ. 5 .OR. ITYPG .EQ. 20)
     >  NDIM=2
      IF(ITYPG .EQ.  5 .OR. ITYPG .EQ. 7) THEN
        NR=1
        NRS=1
      ELSE IF(ITYPG .EQ.  3 .OR. ITYPG .EQ.  6) THEN
        NR=NM(4)
        NRS=NMS(4)
        NRZ=NR
        NRZS=NRS
        NSZS=0
      ELSE IF(ITYPG .EQ. 10 ) THEN
        NR=NM(4)
        NRS=NMS(4)
        NRX=NR
        NRXS=NRS
        NSXS=0
      ELSE IF(ITYPG .EQ. 11 ) THEN
        NR=NM(4)
        NRS=NMS(4)
        NRY=NR
        NRYS=NRS
        NSYS=0
      ELSE IF(ITYPG .EQ. 20 .OR. ITYPG .EQ. 23 ) THEN
        NR=NM(4)+1
        NRS=NMS(4)+1
        NRZ=NR
        NRZS=NRS
        IF(DAMESH(-1,4) .NE. DZERO) NRXS=NRS
      ELSE IF(ITYPG .EQ. 21 ) THEN
        NR=NM(4)+1
        NRS=NMS(4)+1
        NRX=NR
        NRXS=NRS
      ELSE IF(ITYPG .EQ. 22) THEN
        NR=NM(4)+1
        NRS=NMS(4)+1
        NRY=NR
        NRYS=NRS
      ENDIF
      CALL XDISET(IDREG,NREGS,0)
      CALL XDISET(IDSUR,NSURS,0)
*----
*  Regions
*  Mixture and global numbering
*----
      IZS=0
      DO IZ=1,NZ
        NZZ=MAX(1,ISPLT(IZ,3))
        IOFZ=(IZ-1)*NY
        DO IZZ=1,NZZ
          IZS=IZS+1
          IOSZ=(IZS-1)*NYS
          IYS=0
          DO IY=1,NY
            NYY=MAX(1,ISPLT(IY,2))
            IOFYZ=(IOFZ+(IY-1))*NX
            DO  IYY=1,NYY
              IYS=IYS+1
              IOSYZ=(IOSZ+(IYS-1))*NXS
              IXS=0
              DO IX=1,NX
                NXX=MAX(1,ISPLT(IX,1))
                IOFXYZ=(IOFYZ+(IX-1))*NR
                DO IXX=1,NXX
                  IXS=IXS+1
                  IOSXYZ=(IOSYZ+(IXS-1))*NRS
                  IRS=0
                  DO IR=1,NM(4)
                    NRR=MAX(1,ABS(ISPLT(IR,4)))
                    IMIX=IOFXYZ+IR
                    DO IRR=1,NRR
                      IRS=IRS+1
                      IMIXS=IOSXYZ+IRS
                      MIXS(IMIXS,1)=MIX(IMIX,1)
                      MIXS(IMIXS,2)=MIX(IMIX,2)
                      IDREG(IMIXS)=IMIXS
                    ENDDO
                  ENDDO
                  DO IR=NM(4)+1,NR
                    NRR=1
                    IMIX=IOFXYZ+IR
                    DO IRR=1,NRR
                      IRS=IRS+1
                      IMIXS=IOSXYZ+IRS
                      MIXS(IMIXS,1)=MIX(IMIX,1)
                      MIXS(IMIXS,2)=MIX(IMIX,2)
                      IDREG(IMIXS)=IMIXS
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
*----
*  Surfaces (X- and X+)
*----
      ISBOT=0
      ISTOP=ISBOT+NSZS*NSYS*NRXS
      DO IZS=1,NSZS
        DO IYS=1,NSYS
          DO IRS=1,NRXS
            ISBOT=ISBOT+1
            IDSUR(ISBOT)=ISBOT
            ISTOP=ISTOP+1
            IDSUR(ISTOP)=ISTOP
          ENDDO
        ENDDO
      ENDDO
*----
*  Surfaces (Y- and Y+)
*----
      ISBOT=ISTOP
      ISTOP=ISBOT+NSZS*NSXS*NRYS
      DO IXS=1,NSXS
        DO IZS=1,NSZS
          DO IRS=1,NRYS
            ISBOT=ISBOT+1
            IDSUR(ISBOT)=ISBOT
            ISTOP=ISTOP+1
            IDSUR(ISTOP)=ISTOP
          ENDDO
        ENDDO
      ENDDO
      IF(NDIM .EQ. 3) THEN
*----
*  Surfaces (Z- and Z+)
*----
        ISBOT=ISTOP
        ISTOP=ISBOT+NSYS*NSXS*NRZS
        DO IYS=1,NSYS
          DO IXS=1,NSXS
            DO IRS=1,NRZS
              ISBOT=ISBOT+1
              IDSUR(ISBOT)=ISBOT
              ISTOP=ISTOP+1
              IDSUR(ISTOP)=ISTOP
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF(NSXS .EQ. 0 .OR. NSYS .EQ. 0 .OR. NSZS .EQ. 0) THEN
        ISBOT=ISTOP+1
        IDSUR(ISBOT)=ISBOT
      ENDIF
*----
*  For tubes reset outer Cartesian mesh limits in normal plane
*----
      IF(ITYPG .EQ.  3 .OR. ITYPG .EQ.  6) THEN
        DAMESS(0,1,1)=-DAMESS(-1,1,1)-DAMESS(NMS(4),4,1)*DCAOF
        DAMESS(NMS(1),1,1)=-DAMESS(-1,1,1)+DAMESS(NMS(4),4,1)*DCAOF
        DAMESS(0,2,1)=-DAMESS(-1,2,1)-DAMESS(NMS(4),4,1)*DCAOF
        DAMESS(NMS(2),2,1)=-DAMESS(-1,2,1)+DAMESS(NMS(4),4,1)*DCAOF
      ELSE IF(ITYPG .EQ.  10) THEN
        DAMESS(0,2,1)=-DAMESS(-1,2,1)-DAMESS(NMS(4),4,1)*DCAOF
        DAMESS(NMS(2),2,1)=-DAMESS(-1,2,1)+DAMESS(NMS(4),4,1)*DCAOF
        DAMESS(0,3,1)=-DAMESS(-1,3,1)-DAMESS(NMS(4),4,1)*DCAOF
        DAMESS(NMS(3),3,1)=-DAMESS(-1,3,1)+DAMESS(NMS(4),4,1)*DCAOF
      ELSE IF(ITYPG .EQ.  11) THEN
        DAMESS(0,3,1)=-DAMESS(-1,3,1)-DAMESS(NMS(4),4,1)*DCAOF
        DAMESS(NMS(3),3,1)=-DAMESS(-1,3,1)+DAMESS(NMS(4),4,1)*DCAOF
        DAMESS(0,1,1)=-DAMESS(-1,1,1)-DAMESS(NMS(4),4,1)*DCAOF
        DAMESS(NMS(1),1,1)=-DAMESS(-1,1,1)+DAMESS(NMS(4),4,1)*DCAOF
      ENDIF
*----
*  Save MESH and mixture information on IPTRK
*----
      DO IDIR=1,4
        WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'SM'//CDIR(IDIR)
        IF(NMS(IDIR) .GT. 0) THEN
          CALL LCMPUT(IPTRK,NAMREC,(NMS(IDIR)+2),4,DAMESS(-1,IDIR,1))
        ENDIF
      ENDDO
      WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'MIX'
      CALL LCMPUT(IPTRK,NAMREC,NMIXS,1,MIXS(1,1))
      WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'HOM'
      CALL LCMPUT(IPTRK,NAMREC,NMIXS,1,MIXS(1,2))
      IF(IPLOC .GE. 100) THEN
        WRITE(IOUT,'(A38)') 'Regions and surfaces before symmetries'
        WRITE(IOUT,6030)
        WRITE(IOUT,6034)  (IDREG(IDV),IDV=1,NREGS)
        WRITE(IOUT,6032)
        WRITE(IOUT,6034)  (IDSUR(IDS),IDS=1,NSURS)
        WRITE(IOUT,6035) ITSYM
      ENDIF
*----
*  Diagonal X=Y symmetry
*----
      IF(ABS(ITSYM(4)) .EQ. 1) THEN
*----
*  1- Eliminate reflected regions
*----
        DO IZS=1,NZS
          IODZ=(IZS-1)*NYS
          IORZ=(IZS-1)*NYS
          DO IYS=1,NYS
            DO IXS=IYS+1,NYS
              IODXYZ=((IODZ+(IYS-1))*NXS+(IXS-1))*NRS
              IORXYZ=((IODZ+(IXS-1))*NXS+(IYS-1))*NRS
              DO IRS=1,NRS
                IDREG(IORXYZ+IRS)=-ABS(IDREG(IODXYZ+IRS))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
*----
*  2- Eliminate reflected Y- surfaces
*----
        ISBOT=0
        ISTOP=2*NSZS*NSYS*NRYS
        DO IZS=1,NSZS
          DO IYS=1,NSYS
            IXS=IYS
*----
*  X face
*----
            IODYZ=((IZS-1)*NSYS+(IYS-1))*NRYS
*----
*  Y face
*----
            IORYZ=((IXS-1)*NSZS+(IZS-1))*NRYS
            DO IRS=1,NRYS
              IDSUR(ISTOP+IORYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
            ENDDO
          ENDDO
        ENDDO
*----
*  3- Eliminate reflected Y+ surfaces
*----
        ISBOT=NSZS*NSYS*NRYS
        ISTOP=ISBOT+2*NSZS*NSYS*NRYS
        DO IZS=1,NSZS
          DO IYS=1,NSYS
            IXS=IYS
*----
*  X face
*----
            IODYZ=((IZS-1)*NSYS+(IYS-1))*NRYS
*----
*  Y face
*----
            IORYZ=((IXS-1)*NSZS+(IZS-1))*NRYS
            DO IRS=1,NRYS
              IDSUR(ISTOP+IORYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
            ENDDO
          ENDDO
        ENDDO
        IF(NDIM .EQ. 3) THEN
*----
*  4- Eliminate reflected Z+ surfaces
*----
          ISBOT=4*NSZS*NSYS*NRYS
          ISTOP=ISBOT+NSYS*NSXS*NRZS
          DO IYS=1,NSYS
            DO IXS=IYS+1,NSXS
              IODYZ=((IYS-1)*NSXS+(IXS-1))*NRS
              IORYZ=((IXS-1)*NSXS+(IYS-1))*NRS
              DO IRS=1,NRZS
                IDSUR(ISBOT+IORYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
                IDSUR(ISTOP+IORYZ+IRS)=-ABS(IDSUR(ISTOP+IODYZ+IRS))
              ENDDO
            ENDDO
          ENDDO
        ENDIF
*----
*  Diagonal X=-Y symmetry
*----
      ELSE IF(ABS(ITSYM(4)) .EQ. 4) THEN
*----
*  1- Eliminate reflected regions
*----
        DO IZS=1,NZS
          IODZ=(IZS-1)*NYS
          IORZ=(IZS-1)*NYS
          DO IYS=1,NYS
            DO IXS=NYS-IYS,1,-1
              IODXYZ=((IODZ+(IYS-1))*NXS+(IXS-1))*NRS
              IORXYZ=((IODZ+(NYS-IXS))*NXS+(NYS-IYS))*NRS
              DO IRS=1,NRS
                IDREG(IORXYZ+IRS)=-ABS(IDREG(IODXYZ+IRS))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
*----
*  2- Eliminate reflected Y+ surfaces
*----
        ISBOT=0
        ISTOP=ISBOT+3*NSZS*NSYS*NRYS
        DO IZS=1,NSZS
          DO IYS=1,NSYS
            IXS=IYS
*----
*  X face
*----
            IODYZ=((IZS-1)*NSYS+(IYS-1))*NRYS
*----
*  Y face
*----
            IORYZ=((NSYS-IXS)*NSZS+(IZS-1))*NRYS
            DO IRS=1,NRYS
              IDSUR(ISTOP+IORYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
            ENDDO
          ENDDO
        ENDDO
*----
*  3- Eliminate reflected Y- surfaces
*----
        ISBOT=NSZS*NSYS*NRYS
        ISTOP=ISBOT+NSZS*NSYS*NRYS
        DO IZS=1,NSZS
          DO IYS=1,NSYS
            IXS=IYS
*----
*  X face
*----
            IODYZ=((IZS-1)*NSYS+(IYS-1))*NRYS
*----
*  Y face
*----
            IORYZ=((NYS-IXS)*NSZS+(IZS-1))*NRYS
            DO IRS=1,NRYS
              IDSUR(ISTOP+IORYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
            ENDDO
          ENDDO
        ENDDO
        IF(NDIM .EQ. 3) THEN
*----
*  4- Eliminate reflected Z+ surfaces
*----
          ISBOT=4*NSZS*NSYS*NRYS
          ISTOP=ISBOT+NSYS*NSXS*NRZS
          DO IYS=1,NSYS
            DO IXS=NYS-IYS,1,-1
              IODYZ=((IYS-1)*NSXS+(IXS-1))*NRS
              IORYZ=((NYS-IXS)*NSXS+(NYS-IYS))*NRS
              DO IRS=1,NRZS
                IDSUR(ISBOT+IORYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
                IDSUR(ISTOP+IORYZ+IRS)=-ABS(IDSUR(ISTOP+IODYZ+IRS))
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
*----
*  Axial Z+ or Z- symmetry
*----
      IF(ABS(ITSYM(3)) .EQ. 1) THEN
*----
*  1- Eliminate reflected regions
*----
        NZSR=NZS/2
        DO IZS=1,NZSR
          IODZ=(IZS-1)*NYS
          IORZ=(NZS-IZS)*NYS
          DO IYS=1,NYS
            IODYZ=(IODZ+(IYS-1))*NXS
            IORYZ=(IORZ+(IYS-1))*NXS
            DO IXS=1,NXS
              IODXYZ=(IODYZ+(IXS-1))*NRS
              IORXYZ=(IORYZ+(IXS-1))*NRS
              DO IRS=1,NRS
                IDREG(IORXYZ+IRS)=-ABS(IDREG(IODXYZ+IRS))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
*----
*  2- Eliminate reflected X- and X+ surfaces
*----
        ISBOT=0
        ISTOP=ISBOT+NSZS*NSYS*NRXS
        DO IZS=1,NSZS/2
          IODZ=(IZS-1)*NSYS
          IORZ=(NSZS-IZS)*NSYS
          DO IYS=1,NSYS
            IODYZ=(IODZ+(IYS-1))*NRXS
            IORYZ=(IORZ+(IYS-1))*NRXS
            DO IRS=1,NRXS
              IDSUR(ISBOT+IORYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
              IDSUR(ISTOP+IORYZ+IRS)=-ABS(IDSUR(ISTOP+IODYZ+IRS))
            ENDDO
          ENDDO
        ENDDO
*----
*  3- Eliminate reflected Y- and Y+ surfaces
*----
        ISBOT=ISTOP+NSZS*NSYS*NRXS
        ISTOP=ISBOT+NSZS*NSXS*NRYS
        DO IXS=1,NSXS
          IODZ=(IXS-1)*NSZS
          IORZ=(IXS-1)*NSZS
          DO IZS=1,NSZS/2
            IODYZ=(IODZ+(IZS-1))*NRYS
            IORYZ=(IORZ+(NSZS-IZS))*NRYS
            DO IRS=1,NRYS
              IDSUR(ISBOT+IORYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
              IDSUR(ISTOP+IORYZ+IRS)=-ABS(IDSUR(ISTOP+IODYZ+IRS))
            ENDDO
          ENDDO
        ENDDO
*----
*  4- Eliminate reflected Z- surfaces
*----
        ISBOT=ISTOP+NSZS*NSXS*NRYS
        ISTOP=ISBOT+NSYS*NSXS*NRZS
        DO IYS=1,NSYS
          IODZ=(IYS-1)*NSXS
          DO IXS=1,NSXS
            IODYZ=(IODZ+(IXS-1))*NRZS
            DO IRS=1,NRZS
              IDSUR(ISTOP+IODYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
            ENDDO
          ENDDO
        ENDDO
      ENDIF
*----
*  Axial Y+ or Y- symmetry
*----
      IF(ABS(ITSYM(2)) .EQ. 1) THEN
*----
*  1- Eliminate reflected regions
*----
        NYSR=NYS/2
        DO IZS=1,NZS
          IODZ=(IZS-1)*NYS
          IORZ=(IZS-1)*NYS
          DO IYS=1,NYSR
            IODYZ=(IODZ+(IYS-1))*NXS
            IORYZ=(IORZ+(NYS-IYS))*NXS
            DO IXS=1,NXS
              IODXYZ=(IODYZ+(IXS-1))*NRS
              IORXYZ=(IORYZ+(IXS-1))*NRS
              DO IRS=1,NRS
                IDREG(IORXYZ+IRS)=-ABS(IDREG(IODXYZ+IRS))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
*----
*  2- Eliminate reflected X- and X+ surfaces
*----
        ISBOT=0
        ISTOP=ISBOT+NSZS*NSYS*NRXS
        DO IZS=1,NSZS
          IODZ=(IZS-1)*NSYS
          IORZ=(IZS-1)*NSYS
          DO IYS=1,NSYS/2
            IODYZ=(IODZ+(IYS-1))*NRXS
            IORYZ=(IORZ+(NSYS-IYS))*NRXS
            DO IRS=1,NRXS
              IDSUR(ISBOT+IORYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
              IDSUR(ISTOP+IORYZ+IRS)=-ABS(IDSUR(ISTOP+IODYZ+IRS))
            ENDDO
          ENDDO
        ENDDO
*----
*  3- Eliminate reflected Y+ surfaces
*----
        ISBOT=ISTOP+NSZS*NSYS*NRXS
        ISTOP=ISBOT+NSZS*NSXS*NRYS
        DO IXS=1,NSXS
          IODZ=(IXS-1)*NSZS
          DO IZS=1,NSZS
            IODYZ=(IODZ+(IZS-1))*NRYS
            DO IRS=1,NRYS
              IDSUR(ISTOP+IODYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
            ENDDO
          ENDDO
        ENDDO
        IF(NDIM .EQ. 3) THEN
*----
*  4- Eliminate reflected Z- and Z+ surfaces
*----
          ISBOT=ISTOP+NSZS*NSXS*NRYS
          ISTOP=ISBOT+NSYS*NSXS*NRZS
          DO IYS=1,NSYS/2
            IODZ=(IYS-1)*NSXS
            IORZ=(NSYS-IYS)*NSXS
            DO IXS=1,NSXS
              IODYZ=(IODZ+(IXS-1))*NRZS
              IORYZ=(IORZ+(IXS-1))*NRZS
              DO IRS=1,NRZS
                IDSUR(ISBOT+IORYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
                IDSUR(ISTOP+IORYZ+IRS)=-ABS(IDSUR(ISTOP+IODYZ+IRS))
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
*----
*  Axial X+ or X- symmetry
*----
      IF(ABS(ITSYM(1)) .EQ. 1) THEN
*----
*  1- Eliminate reflected regions
*----
        NXSR=NXS/2
        DO IZS=1,NZS
          IODZ=(IZS-1)*NYS
          IORZ=(IZS-1)*NYS
          DO IYS=1,NYS
            IODYZ=(IODZ+(IYS-1))*NXS
            IORYZ=(IORZ+(IYS-1))*NXS
            DO IXS=1,NXSR
              IODXYZ=(IODYZ+(IXS-1))*NRS
              IORXYZ=(IORYZ+(NXS-IXS))*NRS
              DO IRS=1,NRS
                IDREG(IORXYZ+IRS)=-ABS(IDREG(IODXYZ+IRS))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
*----
*  2- Eliminate reflected X+ surfaces
*----
        ISBOT=0
        ISTOP=ISBOT+NSZS*NSYS*NRXS
        DO IZS=1,NSZS
          IODZ=(IZS-1)*NSYS
          DO IYS=1,NSYS
            IODYZ=(IODZ+(IYS-1))*NRXS
            DO IRS=1,NRXS
              IDSUR(ISTOP+IODYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
            ENDDO
          ENDDO
        ENDDO
*----
*  3- Eliminate reflected Y- and Y+ surfaces
*----
        ISBOT=ISTOP+NSZS*NSYS*NRXS
        ISTOP=ISBOT+NSZS*NSXS*NRYS
        DO IXS=1,NSXS/2
          IODZ=(IXS-1)*NSZS
          IORZ=(NSXS-IXS)*NSZS
          DO IZS=1,NSZS
            IODYZ=(IODZ+(IZS-1))*NRYS
            IORYZ=(IORZ+(IZS-1))*NRYS
            DO IRS=1,NRYS
              IDSUR(ISBOT+IORYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
              IDSUR(ISTOP+IORYZ+IRS)=-ABS(IDSUR(ISTOP+IODYZ+IRS))
            ENDDO
          ENDDO
        ENDDO
        IF(NDIM .EQ. 3) THEN
*----
*  4- Eliminate reflected Z- and Z+ surfaces
*----
          ISBOT=ISTOP+NSZS*NSXS*NRYS
          ISTOP=ISBOT+NSYS*NSXS*NRZS
          DO IYS=1,NSYS
            IODZ=(IYS-1)*NSXS
            IORZ=(IYS-1)*NSXS
            DO IXS=1,NSXS/2
              IODYZ=(IODZ+(IXS-1))*NRZS
              IORYZ=(IORZ+(NSXS-IXS))*NRZS
              DO IRS=1,NRZS
                IDSUR(ISBOT+IORYZ+IRS)=-ABS(IDSUR(ISBOT+IODYZ+IRS))
                IDSUR(ISTOP+IORYZ+IRS)=-ABS(IDSUR(ISTOP+IODYZ+IRS))
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
*----
*  For annular geometry eliminate Cartesian boundary
*----
      IF(ITYPG .EQ.  6) THEN
        NR=NM(4)
        NRS=NMS(4)
        NRZ=NR
        NRZS=NRS
      ELSE IF(ITYPG .EQ. 10 ) THEN
        NR=NM(4)
        NRS=NMS(4)
        NRX=NR
        NRXS=NRS
      ELSE IF(ITYPG .EQ. 11 ) THEN
        NR=NM(4)
        NRS=NMS(4)
        NRY=NR
        NRYS=NRS
      ENDIF
*----
*  Renumber regions
*----
      IRN=0
      DO IRO=1,NREGS
        IDV=IDREG(IRO)
        IF(IDV .GT. 0) THEN
          IRN=IRN+1
          IREN(IRN)=IDV
          DO IRT=1,IRN-1
            IF(IDV .LT. IREN(IRT)) THEN
              DO IRR=IRN-1,IRT,-1
                IREN(IRR+1)=IREN(IRR)
              ENDDO
              IREN(IRT)=IDV
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      NREGN=IRN
      DO IRN=1,NREGN
        IDV=IREN(IRN)
        DO IRT=1,NREGS
          IF(IDREG(IRT) .EQ. IDV ) THEN
            IDREG(IRT)=IRN
          ELSE IF(IDREG(IRT) .EQ. -IDV ) THEN
            IDREG(IRT)=-IRN
          ENDIF
        ENDDO
      ENDDO
*----
*  Renumber surfaces
*----
      IRN=0
      DO IRO=1,NSURS
        IDS=IDSUR(IRO)
        IF(IDS .GT. 0) THEN
          IRN=IRN+1
          IREN(IRN)=IDS
          DO IRT=1,IRN-1
            IF(IDS .LT. IREN(IRT)) THEN
              DO IRR=IRN-1,IRT,-1
                IREN(IRR+1)=IREN(IRR)
              ENDDO
              IREN(IRT)=IDS
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      NSURN=IRN
      DO IRN=1,NSURN
        IDS=IREN(IRN)
        DO IRT=1,NSURS
          IF(IDSUR(IRT) .EQ. IDS ) THEN
            IDSUR(IRT)=IRN
          ELSE IF(IDSUR(IRT) .EQ. -IDS ) THEN
            IDSUR(IRT)=-IRN
          ENDIF
        ENDDO
      ENDDO
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPLOC .GE. 100) THEN
        WRITE(IOUT,6010) (NMS(IDIR),IDIR=1,4)
        DO IDIR=1,4
          NMTMP=NMS(IDIR)
          IF(NMTMP .GT. 0) THEN
            WRITE(IOUT,6011) 'MESH'//CDIR(IDIR)//' ='
            WRITE(IOUT,6012) (DAMESS(IX,IDIR,1),IX=-1,NMTMP)
          ENDIF
        ENDDO
        WRITE(IOUT,'(A37)') 'Regions and surfaces after symmetries'
        WRITE(IOUT,6030)
        WRITE(IOUT,6034)  (IDREG(IDV),IDV=1,NREGS)
        WRITE(IOUT,6032)
        WRITE(IOUT,6034)  (IDSUR(IDS),IDS=1,NSURS)
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      DEALLOCATE(MIXS,IREN)
      RETURN
*----
*  FORMATS
*----
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(1X,' DIMENSIONS =',5I10/1X,' ORIGINAL MESH ')
 6011 FORMAT(1X,A7)
 6012 FORMAT(5F15.9)
 6013 FORMAT(5I15)
 6020 FORMAT(1X,' DIMENSIONS =',5I10/1X,' SPLITTED MESH ')
 6030 FORMAT(' Regions ID')
 6032 FORMAT(' Surfaces ID')
 6034 FORMAT(5I15)
 6035 FORMAT('Symmetries =',4I15)
      END
