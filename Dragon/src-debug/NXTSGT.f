*DECK NXTSGT
      SUBROUTINE NXTSGT(IPTRK ,IPRINT,MAXMSH,ITYPG ,IGEO  ,ILEV  ,
     >                  MAXMSS,NMIX  ,NM    ,MIX   ,DAMESH,ISPLT ,
     >                  NMIXS ,NMS   ,DAMESS,
     >                  ITSYM ,NREGS ,NSURS ,NREGN ,NSURN ,NEREN ,
     >                  IDREG ,IDSUR )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Discretize geometry according to splitting options for
* HEXT, HEXCELT, HEXTZ and HEXTCELZ geometry.
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
* ISPLT   final split desctiption for geometry.
* NMIXS   number of regional mixtures.
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
      PARAMETER       (IOUT=6,NAMSBR='NXTSGT')
      DOUBLE PRECISION DZERO,DONE,DDC,DCAOF
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DDC=1.0D-3,
     >                 DCAOF=DONE+DDC)
*----
*  Local variables
*----
      INTEGER          NX ,NZ ,NRTP ,NRP ,NSTP ,NSP ,
     >                 NXS,NZS,NRTPS,NRPS,NSTPS,NSPS,
     >                 NZMX1,NZSMX1,NBSX,NBSZ,NZSR,NSPLAN
      INTEGER          IPLOC,IX,IXS,IZ,IZS,III,JJJ,
     >                 ISECT,IDMIX
      INTEGER          IOFZ,IOFZS,IOFS,IOFSS
      INTEGER          IODZ,IODHZ,IORHZ,IORZ,IRN,IRO,IRR,IRT,
     >                 ISTOP,ISBOT,IDS,IDV
      DOUBLE PRECISION DELTH,DELTH0,DDD
      CHARACTER        NAMREC*12
      INTEGER          IDIRR,IR,NR,NRS
      INTEGER          IDTRI,IGEN,IGENS,IOFT,IOFTS,IS,ISPZ,ITID,ITIDS,
     >                 NBS,NRR,NTL
      DOUBLE PRECISION DDI,DDO
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IREN,MIXS
*----
*  Data
*----
      CHARACTER        CLEV(2)*1
      SAVE             CLEV
      DATA             CLEV /'C','P'/
*----
*  Static variables used
*  before split
*  NX    =  mesh in x (from 0 to SIDE)
*  NZ    =  mesh in z (from zmin to zmax)
*  NR    =  radial mesh
*  NRTP  =  number of region in a triangle for each z plane
*  NRP   =  number of region in a z plane
*  NSTP  =  number of surfaces in a triangle for each z plane
*  NSP   =  number of surfaces in a z plane
*  after split
*  NXS   =  mesh in x (from 0 to SIDE)
*  NZS   =  mesh in z (from zmin to zmax)
*  NRS   =  radial mesh
*  NRTPS =  number of region in a triangle for each z plane
*  NRPS  =  number of region in a z plane
*  NSTPS =  number of surfaces in a triangle for each z plane
*  NSPS  =  number of surfaces in a z plane
*----
      ALLOCATE(IREN(2,NEREN),MIXS(NMIXS,2))
      CALL XDISET(MIXS,NMIXS*2,0)
      NR=NM(4)+1
      NRS=NMS(4)+1
      IDIRR=-99
      IF(ITYPG .EQ.  12 .OR. ITYPG .EQ. 13) THEN
        IDIRR=0
      ELSE IF(ITYPG .EQ.  26 .OR. ITYPG .EQ. 27) THEN  
        IDIRR=4
      ENDIF
      NX=NM(1)
      NZ=NM(3)
      NRTP=NX**2
      NRP=6*NRTP*NR
      NSTP=2*NX-1
      NSP=6*NSTP
      NXS=NMS(1)
      NZS=NMS(3)
      NRTPS=NXS**2
      NRPS=6*NRTPS*NRS
      NSTPS=2*NXS-1
      NSPS=6*NSTPS
*      write(6,*) 'NX,NZ,NRTP,NRP,NSTP,NSP,NXS,NZS,NRTPS,NRPS,NSTPS,NSPS'
*     >,NX,NZ,NRTP,NRP,NSTP,NSP,NXS,NZS,NRTPS,NRPS,NSTPS,NSPS
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IPLOC=IPRINT
      IF(IPLOC .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6002) ITYPG
        WRITE(IOUT,6010) NX,NZ
        WRITE(IOUT,6011) 'MESHH ='
        WRITE(IOUT,6012) (DAMESH(IX,1),IX=-1,2*NX)
        WRITE(IOUT,6011) 'SPLTH ='
        WRITE(IOUT,6013) ISPLT(1,1)
        IF(NZ .GT. 0) THEN
          WRITE(IOUT,6011) 'MESHZ ='
          WRITE(IOUT,6012) (DAMESH(IZ,3),IZ=-1,NZ)
          WRITE(IOUT,6011) 'SPLTZ ='
          WRITE(IOUT,6013) (ISPLT(IZ,3),IZ=1,NZ)
        ENDIF
        IF(IDIRR .EQ. 4) THEN
          WRITE(IOUT,6011) 'MESHR ='
          WRITE(IOUT,6012) (DAMESH(IR,4),IR=-1,NR)
          WRITE(IOUT,6011) 'SPLTR ='
          WRITE(IOUT,6013) (ISPLT(IR,4),IR=1,NR)
        ENDIF
        WRITE(IOUT,6033)
        WRITE(IOUT,6034)  (MIX(III,1),III=1,NMIX)
        WRITE(IOUT,6036)
        WRITE(IOUT,6034)  (MIX(III,2),III=1,NMIX)
      ENDIF
*      write(6,*) 'NX,NZ,NXS,NZS',NX,NZ,NXS,NZS
*      write(6,*) 'NREGS,NSURS',NREGS,NSURS
*----
*  Triangular mesh
*----
      NBSX=ABS(ISPLT(1,1))
*----
*  Note : last crown may have a width different
*  from the other crowns
*  Splitting is with respect to width of internal crowns
*  unless only one crown considered
*----
      IF(NX .EQ. 1) THEN
        DELTH0=(DAMESH(1,1)-DAMESH(0,1))/DBLE(NBSX)
        DELTH=DELTH0
      ELSE
        DELTH0=DAMESH(1,1)-DAMESH(0,1)
        DELTH=DAMESH(2,1)-DAMESH(1,1)
        DDD=DELTH-DELTH0
        DELTH=DELTH/DBLE(NBSX)
        IF(DDD .GT. DELTH) THEN
          WRITE(IOUT,9000) DELTH0,DELTH,DDD
          CALL XABORT(NAMSBR//
     >      ': Invalid SPLIT for HEXT and HEXTZ geometries')
        ENDIF
        DELTH0=DELTH-DDD
      ENDIF
      DAMESS(-1,1,1)=DAMESH(-1,1)
      DAMESS(0,1,1)=DAMESH(0,1)
      DAMESS(1,1,1)=DAMESS(0,1,1)+DELTH0
      III=1
      DO IXS=1,NBSX-1
        III=III+1
        DAMESS(III,1,1)=DAMESS(III-1,1,1)+DELTH
      ENDDO
      DO JJJ=2,2*NX-1
        DO IXS=1,NBSX
          III=III+1
          DAMESS(III,1,1)=DAMESS(III-1,1,1)+DELTH
        ENDDO
      ENDDO
      DO IXS=1,NBSX-1
        III=III+1
        DAMESS(III,1,1)=DAMESS(III-1,1,1)+DELTH
      ENDDO
      III=III+1
      DAMESS(III,1,1)=DAMESS(III-1,1,1)+DELTH0
*----
*  Y mesh identical to X mesh except for off center
*----
      DO IXS=-1,2*NXS
        DAMESS(IXS,2,1)=DAMESS(IXS,1,1)
      ENDDO
*----
*  Z mesh
*----
      IF(NZ .EQ. 0) THEN
        DAMESS(-1,3,1)=DZERO
        DAMESS(0,3,1)=DZERO
      ELSE
        DAMESS(-1,3,1)=DAMESH(-1,3)
        DAMESS(0,3,1)=DAMESH(0,3)
        III=0
        DO IZ=1,NZ
          NBSZ=ABS(ISPLT(IZ,3))
          DDD=(DAMESH(IZ,3)-DAMESH(IZ-1,3))/DBLE(NBSZ)
          DO IZS=1,NBSZ
            III=III+1
            DAMESS(III,3,1)=DAMESS(III-1,3,1)+DDD
          ENDDO
        ENDDO
      ENDIF
*----
* Radial mesh
*----
      IF(IDIRR .EQ. 4) THEN
        DAMESS(-1,IDIRR,1)=DAMESH(-1,IDIRR)
        DAMESS(0,IDIRR,1)=DAMESH(0,IDIRR)
        IGENS=0
        DO IGEN=1,NM(IDIRR)
          NBS=ISPLT(IGEN,IDIRR)
          IF(NBS .LT. 0) THEN
            NBS=-NBS
            DDI=DAMESH(IGEN-1,IDIRR)*DAMESH(IGEN-1,IDIRR)
            DDO=DAMESH(IGEN,IDIRR)*DAMESH(IGEN,IDIRR)
            DDD=(DDO-DDI)/DBLE(NBS)
            DO IS=1,NBS
              IGENS=IGENS+1
              DDO=DDI+DDD
              DAMESS(IGENS,IDIRR,1)=SQRT(DDO)
              DDI=DDO
            ENDDO
          ELSE
            DDD=(DAMESH(IGEN,IDIRR)-DAMESH(IGEN-1,IDIRR))/DBLE(NBS)
            DO IS=1,NBS
              IGENS=IGENS+1
              DAMESS(IGENS,IDIRR,1)=DAMESS(IGENS-1,IDIRR,1)+DDD
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      IF(IPLOC .GE. 100) THEN
        WRITE(IOUT,6020) NXS,NZS
        WRITE(IOUT,6011) 'MESHH ='
        WRITE(IOUT,6012) (DAMESS(IX,1,1),IX=-1,2*NXS)
        IF(NZS .GT. 0) THEN
          WRITE(IOUT,6011) 'MESHZ ='
          WRITE(IOUT,6012) (DAMESS(IZ,3,1),IZ=-1,NZS)
        ENDIF
        IF(IDIRR .EQ. 4) THEN
          WRITE(IOUT,6011) 'MESHR ='
          WRITE(IOUT,6012) (DAMESS(IR,4,1),IR=-1,NR)
        ENDIF
      ENDIF
      CALL XDISET(IDREG,NREGS,0)
      CALL XDISET(IDSUR,NSURS,0)
*----
*  Loop over planes (3D)
*----
      NZMX1=MAX(1,NZ)
      NZSMX1=MAX(1,NZS)
*----
*  Regions
*  Mixture and global numbering
*----
      IZS=0
      DO IZ=1,NZMX1
        IF(NZ .EQ. 0) THEN
          NBSZ=1
        ELSE
          NBSZ=MAX(1,ISPLT(IZ,3))
        ENDIF
        DO ISPZ=1,NBSZ
          IOFZ=(IZ-1)*NRP
          IOFZS=IZS*NRPS
          IZS=IZS+1
*----
*  Loop over sectors
*----
          DO ISECT=1,6
            IOFS=(ISECT-1)*NRTP*NR
            IOFSS=(ISECT-1)*NRTPS*NRS
*----
*  Loop over lines
*----
            IOFT=IOFZ+IOFS
            IOFTS=IOFZS+IOFSS
            ITID=0
            ITIDS=0
            DO IX=1,NX
              NTL=2*IX-1
*----
*  For remaining lines
*  IREN(1,*) is original cell position on line
*  IREN(2,*) is number of times this cell is used on
*  each initial splitted mesh
*  (NBSX for right triangle and 0 for left triangle)
*----
*  Loop over all right triangles in line
*----
              DO III=2,2*IX-2,2
                IREN(1,III)=IOFT
                IREN(2,III)=NBSX
                IOFT=IOFT+NR
              ENDDO
*----
*  Loop over all left triangles in line
*----
              DO III=1,2*IX-1,2
                IREN(1,III)=IOFT
                IREN(2,III)=0
                IOFT=IOFT+NR
              ENDDO
*      write(6,'(A6,2X,3I10)')
*     >     ('IREN =',III,IREN(1,III),IREN(2,III),III=1,2*IX-1)
              DO IXS=1,NBSX
*----
*  Process all triangles on coarse line
*  Extract only new (discretized) right triangles on sub-line
*  In first step (IXS=1), only coarse right triangles considered
*----
                DO III=1,2*IX-1
                  IDTRI=IREN(1,III)
                  DO JJJ=1,IREN(2,III)
                    DO IR=1,NM(4)
                      NRR=MAX(1,ABS(ISPLT(IR,4)))
                      IDMIX=IDTRI+IR
                      DO IRR=1,NRR
                        IOFTS=IOFTS+1
*      write(6,*) 'SPLITTING R',IZ,ISPZ,ISECT,IX,IXS,IR,IRR,IOFTS
                        MIXS(IOFTS,1)=MIX(IDMIX,1)
                        MIXS(IOFTS,2)=MIX(IDMIX,2)
                        IDREG(IOFTS)=IOFTS
                      ENDDO
                    ENDDO
                    IDMIX=IDTRI+NM(4)+1
                    IOFTS=IOFTS+1
*      write(6,*) 'SPLITTING R',IZ,ISPZ,ISECT,IX,IXS,NM(4),NRR+1,IOFTS
                    MIXS(IOFTS,1)=MIX(IDMIX,1)
                    MIXS(IOFTS,2)=MIX(IDMIX,2)
                    IDREG(IOFTS)=IOFTS
                  ENDDO
                ENDDO
*----
*  Prepare for right triangles on sub-line
*  1) Coarse right triangles
*  2) Coarse left triangles
*----
                DO III=2,2*IX-2,2
                  IREN(2,III)=NBSX-1
                ENDDO
                DO III=1,2*IX-1,2
                  IREN(2,III)=IREN(2,III)+1
                ENDDO
*      write(6,'(A6,2X,3I10)')
*     >     ('IREN =',III,IREN(1,III),IREN(2,III),III=1,2*IX-1)
*----
*  Process again all triangles on coarse line
*  Extract only new (discretized) left triangles on sub-line
*  In last step (IXS=NBXS), only coarse left triangles considered
*----
                DO III=1,2*IX-1
                  IDTRI=IREN(1,III)
                  DO JJJ=1,IREN(2,III)
                    DO IR=1,NM(4)
                      NRR=MAX(1,ABS(ISPLT(IR,4)))
                      IDMIX=IDTRI+IR
                      DO IRR=1,NRR
                        IOFTS=IOFTS+1
*      write(6,*) 'SPLITTING L',IZ,ISPZ,ISECT,IX,IXS,IR,IRR,IOFTS
                        MIXS(IOFTS,1)=MIX(IDMIX,1)
                        MIXS(IOFTS,2)=MIX(IDMIX,2)
                        IDREG(IOFTS)=IOFTS
                      ENDDO
                    ENDDO
                    IDMIX=IDTRI+NM(4)+1
                    IOFTS=IOFTS+1
*      write(6,*) 'SPLITTING L',IZ,ISPZ,ISECT,IX,IXS,NM(4),NRR+1,IOFTS
                    MIXS(IOFTS,1)=MIX(IDMIX,1)
                    MIXS(IOFTS,2)=MIX(IDMIX,2)
                    IDREG(IOFTS)=IOFTS
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
*----
*  Save MESH and mixture information on IPTRK
*----
      IF(NXS .GT. 0) THEN
        WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'SMX'
        CALL LCMPUT(IPTRK,NAMREC,(2*NXS+2),4,DAMESS(-1,1,1))
        WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'SMY'
        CALL LCMPUT(IPTRK,NAMREC,(2*NXS+2),4,DAMESS(-1,2,1))
      ENDIF
      IF(NZS .GT. 0) THEN
        WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'SMZ'
        CALL LCMPUT(IPTRK,NAMREC,(NZS+2),4,DAMESS(-1,3,1))
      ENDIF
      IF(IDIRR .EQ. 4) THEN
        WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'SMR'
        CALL LCMPUT(IPTRK,NAMREC,(NMS(4)+2),4,DAMESS(-1,4,1))
      ENDIF
      WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'MIX'
      CALL LCMPUT(IPTRK,NAMREC,NMIXS,1,MIXS(1,1))
      WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'HOM'
      CALL LCMPUT(IPTRK,NAMREC,NMIXS,1,MIXS(1,2))
      IF(IPLOC .GE. 100) THEN
        WRITE(IOUT,6033)
        WRITE(IOUT,6034)  (MIXS(III,1),III=1,NREGS)
        WRITE(IOUT,6036)
        WRITE(IOUT,6034)  (MIXS(III,2),III=1,NREGS)
        WRITE(IOUT,6030) 'before symmetries   '
        WRITE(IOUT,6034)  (IDREG(III),III=1,NREGS)
      ENDIF
*----
*  Hexagnal surfaces
*----
      ISBOT=0
      DO IZS=1,NZSMX1
        DO ISECT=1,6
          DO III=1,NXS-1
            ISBOT=ISBOT+1
            IDSUR(ISBOT)=ISBOT
          ENDDO
          DO III=1,NXS
            ISBOT=ISBOT+1
            IDSUR(ISBOT)=ISBOT
          ENDDO
        ENDDO
      ENDDO
      IF(NZ .GT. 0) THEN
*----
*  Z- and Z+ surfaces
*----
        ISTOP=ISBOT+NRPS
        DO III=1,NRPS
          ISBOT=ISBOT+1
          IDSUR(ISBOT)=ISBOT
          ISTOP=ISTOP+1
          IDSUR(ISTOP)=ISTOP
        ENDDO
      ENDIF
      IF(IPLOC .GE. 100) THEN
        WRITE(IOUT,6032) 'before symmetries   '
        WRITE(IOUT,6034)  (IDSUR(III),III=1,NSURS)
        WRITE(IOUT,6035)   ITSYM
      ENDIF
*----
*  Z symmetry
*----
      IF(ABS(ITSYM(3)) .EQ. 1) THEN
*----
*  1- Eliminate reflected regions
*----
        NZSR=NZS/2
        DO IZS=1,NZSR
          IODZ=(IZS-1)*NRPS
          IORZ=(NZS-IZS)*NRPS
          DO III=1,NRPS
            IODHZ=IODZ+III
            IORHZ=IORZ+III
            IDREG(IORHZ)=-ABS(IDREG(IODHZ))
          ENDDO
        ENDDO
*----
*  2- Eliminate reflected hexagonal surfaces
*----
        NSPLAN=12*NXS-6
        DO IZS=1,NZSR
          IODZ=(IZS-1)*NSPLAN
          IORZ=(NZS-IZS)*NSPLAN
          DO III=1,NSPLAN
            IODHZ=IODZ+III
            IORHZ=IORZ+III
            IDSUR(IORHZ)=-ABS(IDSUR(IODHZ))
          ENDDO
        ENDDO
*----
*  4- Eliminate reflected Z- surfaces
*----
        ISBOT=NSPLAN*NZSR
        ISTOP=ISBOT+NRPS
        DO III=1,NRPS
          ISBOT=ISBOT+1
          ISTOP=ISTOP+1
          IDSUR(ISTOP)=-ABS(IDSUR(ISBOT))
        ENDDO
      ENDIF
*----
*  Renumber regions
*----
      IRN=0
      DO IRO=1,NREGS
        IDV=IDREG(IRO)
        IF(IDV .GT. 0) THEN
          IRN=IRN+1
          IREN(1,IRN)=IDV
          DO IRT=1,IRN-1
            IF(IDV .LT. IREN(1,IRT)) THEN
              DO IRR=IRN-1,IRT,-1
                IREN(1,IRR+1)=IREN(1,IRR)
              ENDDO
              IREN(1,IRT)=IDV
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      NREGN=IRN
      DO IRN=1,NREGN
        IDV=IREN(1,IRN)
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
          IREN(1,IRN)=IDS
          DO IRT=1,IRN-1
            IF(IDS .LT. IREN(1,IRT)) THEN
              DO IRR=IRN-1,IRT,-1
                IREN(1,IRR+1)=IREN(1,IRR)
              ENDDO
              IREN(1,IRT)=IDS
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      NSURN=IRN
      DO IRN=1,NSURN
        IDS=IREN(1,IRN)
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
        WRITE(IOUT,6010) NXS,NZS
        IF(NXS .GT. 0) THEN
          WRITE(IOUT,6011) 'MESHH ='
          WRITE(IOUT,6012) (DAMESS(IX,1,1),IX=-1,2*NXS)
        ENDIF
        IF(NZS .GT. 0) THEN
          WRITE(IOUT,6011) 'MESHZ ='
          WRITE(IOUT,6012) (DAMESS(IX,3,1),IX=-1,NZS)
        ENDIF
        WRITE(IOUT,6030) 'after symmetries    '
        WRITE(IOUT,6034) (IDREG(IDV),IDV=1,NREGS)
        WRITE(IOUT,6032) 'after symmetries    '
        WRITE(IOUT,6034) (IDSUR(IDS),IDS=1,NSURS)
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      DEALLOCATE(MIXS,IREN)
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6002 FORMAT(' Type of geometry = ',I5)
 6010 FORMAT(1X,'DIMENSIONS =',2I10/1X,'ORIGINAL MESH ')
 6011 FORMAT(1X,A7)
 6012 FORMAT(5F20.10)
 6013 FORMAT(5I20)
 6020 FORMAT(1X,'DIMENSIONS =',2I10/1X,'SPLITTED MESH ')
 6030 FORMAT(' Regions ID ',A20)
 6032 FORMAT(' Surfaces ID ',A20)
 6033 FORMAT(' Mixtures ')
 6034 FORMAT(5I15)
 6035 FORMAT('Symmetries =',4I15)
 6036 FORMAT(' Virtual mixtures ')
 9000 FORMAT(' Problem with split in HEXT or HEXTZ geometry'//
     >3F20.10)
      END
