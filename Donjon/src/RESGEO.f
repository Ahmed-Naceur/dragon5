*DECK RESGEO
      SUBROUTINE RESGEO(IPMAP,IPMTX,LX,LY,LZ,NFUEL,IMPX,IGEO,NX,NY,NZ,
     1 NCH,NB,NTOT,LNAP,IPCPO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Create and check the fuel-map geometry.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* E. Varin, D. Sekki and V. Descotes
*
*Update(s):
* R. Chambon 2014
*
*Parameters: input
* IPMAP  pointer to fuel-map information.
* IPMTX  pointer to matex information.
* LX     number of elements along x-axis in geometry.
* LY     number of elements along y-axis in geometry.
* LZ     number of elements along z-axis in geometry.
* NFUEL  number of fuel types.
* IMPX   printing index (=0 for no print).
* IGEO   type of geometry (=7 or =9)
*
*Parameters: output
* NX     number of elements along x-axis in fuel map.
* NY     number of elements along y-axis in fuel map.
* NZ     number of elements along z-axis in fuel map.
* NCH    number of reactor channels.
* NB     number of fuel bundles per channel.
* NTOT   total number of fuel bundles.
* LNAP   Flag to call NAP: module to unfold geometry at assembly level
* IPCPO  pointer to multicompo information
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP,IPMTX,IPCPO,IPGNW
      INTEGER LX,LY,LZ,NFUEL,IGEO,NX,NY,NZ,NCH,NB,NTOT
      LOGICAL LNAP
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40,IOUT=6,EPSI=1.0E-4)
      INTEGER ISTATE(NSTATE),JENT(1),IENT(1),JENT2(3),IENT2(3)
      TYPE(C_PTR) KENT(1),KENT2(3)
      REAL    GEOXX(LX+1),GEOYY(LY+1),GEOZZ(LZ+1),GEOSI,GMAPSI
      CHARACTER HENT(1)*12,HENT2(3)*12,TEXT*12
      DOUBLE PRECISION DFLOT
*----
*  ALLOCATABLE STATEMENTS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: GMAPX,GMAPY,GMAPZ
*----
*  FUEL-MAP GEOMETRY
*----
      IF(IMPX.GT.1)WRITE(IOUT,*)'** CREATING FUEL-MAP GEOMETRY **'
      CALL LCMSIX(IPMAP,'GEOMAP',1)
      NENT=1
      JENT(1)=0
      HENT(1)='GEOMAP'
      IENT(1)=1
      KENT(1)=IPMAP
      CALL GEOD(NENT,HENT,IENT,JENT,KENT)
      IF(IMPX.GT.3)CALL LCMLIB(IPMAP)
      IF(LNAP) THEN
*----
*  FUEL-MAP GEOMETRY UNFOLDING WITH NAP:
*----
        IF(.NOT.C_ASSOCIATED(IPCPO)) THEN
          CALL XABORT('RESGEO: COMPO LCM OBJECT MISSING AT RHS.')
        ENDIF
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@RESGEO: CHARACTER DATA EXPECTED.')
        IF(TEXT.NE.':::') CALL XABORT('@RESGEO: ::: keyword EXPECTED.')
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@RESGEO: CHARACTER DATA EXPECTED.')
        IF(TEXT.NE.'NAP:') CALL XABORT('@RESGEO: NAP: keyword '
     1   //'EXPECTED.')
        CALL LCMOP(IPGNW,'GEONEW',0,1,0)
        CALL LCMSIX(IPMAP,' ',0)
        CALL LCMSIX(IPMAP,'GEOMAP',1)
        NENT2=3
        JENT2(1)=0
        JENT2(2)=2
        JENT2(3)=2
        HENT2(1)='GEONEW'
        HENT2(1)='GEOOLD'
        HENT2(1)='COMPO'
        IENT2(1)=1
        IENT2(2)=1
        IENT2(3)=1
        KENT2(1)=IPGNW
        KENT2(2)=IPMAP
        KENT2(3)=IPCPO
        CALL NAP(NENT2,HENT2,IENT2,JENT2,KENT2)
        CALL LCMSIX(IPMAP,' ',0)
        IF(IMPX.GT.3)CALL LCMLIB(IPMAP)
        CALL LCMDEL(IPMAP,'GEOMAP')
        IF(IMPX.GT.3)CALL LCMLIB(IPMAP)
        CALL LCMSIX(IPMAP,'GEOMAP',1)
        IF(IMPX.GT.3)CALL LCMLIB(IPMAP)
        CALL LCMEQU(IPGNW,IPMAP)
        IF(IMPX.GT.3)CALL LCMLIB(IPMAP)
        CALL LCMCL(IPGNW,1)
      ENDIF
****
      CALL LCMSIX(IPMAP,' ',0)
      IF(IMPX.GT.3)CALL LCMLIB(IPMAP)
      CALL LCMSIX(IPMAP,'GEOMAP',1)
      IF(IMPX.GT.3)CALL LCMLIB(IPMAP)
****
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.IGEO) CALL XABORT('@RESGEO: THE GEOMETRY '
     1 // 'IN FUEL-MAP MUST HAVE THE SAME TYPE AS IN THE MATEX-OBJECT')
      IGEO=ISTATE(1)
      NX=ISTATE(3)
      NY=ISTATE(4)
      NZ=ISTATE(5)
      IF(IMPX.GT.2)WRITE(IOUT,*)'CHECKING FUEL-MAP GEOMETRY'
      IF((IGEO.NE.7).AND.(IGEO.NE.9))CALL XABORT('@RESGEO: ONLY '
     1 //'3D-CARTESIAN OR 3D-HEXAGONAL GEOMETRY ALLOWED.')
      IF((LX.LT.NX).OR.(LY.LT.NY).OR.(LZ.LT.NZ))
     1  CALL XABORT('@RESGEO: WRONG GEOMETRY DEFINITION.')
      IF(NZ.LT.NB)THEN
        WRITE(IOUT,*)'@RESGEO: FOUND NZ=',NZ,' LESS THAN NB=',NB
        CALL XABORT('@RESGEO: WRONG FUEL-MAP GEOMETRY DEFINITION.')
      ENDIF
*----
*  CHECK MESHX OR SIDE
*----
      IF(IGEO.EQ.7) THEN
       ALLOCATE(GMAPX(NX+1))
       CALL XDRSET(GMAPX,NX+1,0.)
       CALL LCMGET(IPMAP,'MESHX',GMAPX)
       CALL XDRSET(GEOXX,LX+1,0.)
       CALL LCMGET(IPMTX,'MESHX',GEOXX)
       DO 10 IMP=1,NX+1
       DO IGM=1,LX+1
         IF(ABS(GMAPX(IMP)-GEOXX(IGM)).LT.EPSI)THEN
           GEOXX(IGM)=GMAPX(IMP)
           GOTO 10
         ENDIF
       ENDDO
       WRITE(IOUT,*)'@RESGEO: MESHX IN L_MAP ',GMAPX(IMP)
       CALL XABORT('@RESGEO: UNABLE TO FIND THIS MESHX IN L_GEOM.')
   10  CONTINUE
       DEALLOCATE(GMAPX)
       CALL LCMPUT(IPMTX,'MESHX',LX+1,2,GEOXX)
      ELSE IF(IGEO.EQ.9) THEN
       NY=1
       CALL LCMGET(IPMAP,'SIDE',GMAPSI)
       CALL LCMGET(IPMTX,'SIDE',GEOSI)
       IF(ABS(GMAPSI-GEOSI).LT.EPSI)THEN
           GEOSI=GMAPSI
           GOTO 20
       ENDIF
       WRITE(IOUT,*)'@RESGEO: SIDE IN L_MAP ',GMAPSI
       CALL XABORT('@RESGEO: UNABLE TO FIND THIS SIDE IN L_GEOM.')
   20  CONTINUE
       CALL LCMPUT(IPMTX,'SIDE',1,2,GEOSI)
      ENDIF
*----
*  CHECK MESHY (ONLY IF 3D-CARTESIAN GEOMETRY)
*----
      IF(IGEO.EQ.7) THEN
       ALLOCATE(GMAPY(NY+1))
       CALL XDRSET(GMAPY,NY+1,0.)
       CALL LCMGET(IPMAP,'MESHY',GMAPY)
       CALL XDRSET(GEOYY,LY+1,0.)
       CALL LCMGET(IPMTX,'MESHY',GEOYY)
       DO 30 IMP=1,NY+1
       DO IGM=1,LY+1
         IF(ABS(GMAPY(IMP)-GEOYY(IGM)).LT.EPSI)THEN
           GEOYY(IGM)=GMAPY(IMP)
           GOTO 30
         ENDIF
       ENDDO
       WRITE(IOUT,*)'@RESGEO: MESHY IN FUEL MAP ',GMAPY(IMP)
       CALL XABORT('@RESGEO: UNABLE TO FIND THIS MESHY IN L_GEOM.')
   30  CONTINUE
       DEALLOCATE(GMAPY)
       CALL LCMPUT(IPMTX,'MESHY',LY+1,2,GEOYY)
      ELSE IF(IGEO.EQ.9) THEN
       GOTO 40
   40  CONTINUE
      ENDIF
*----
*  CHECK MESHZ
*----
      ALLOCATE(GMAPZ(NZ+1))
      CALL XDRSET(GMAPZ,NZ+1,0.)
      CALL LCMGET(IPMAP,'MESHZ',GMAPZ)
      CALL XDRSET(GEOZZ,LZ+1,0.)
      CALL LCMGET(IPMTX,'MESHZ',GEOZZ)
      DO 50 IMP=1,NZ+1
      DO IGM=1,LZ+1
        IF(ABS(GMAPZ(IMP)-GEOZZ(IGM)).LT.EPSI)THEN
          GEOZZ(IGM)=GMAPZ(IMP)
          GOTO 50
        ENDIF
      ENDDO
      WRITE(IOUT,*)'@RESGEO: MESHZ IN FUEL MAP ',GMAPZ(IMP)
      CALL XABORT('@RESGEO: UNABLE TO FIND  THIS MESHZ IN L_GEOM.')
   50 CONTINUE
      DEALLOCATE(GMAPZ)
      CALL LCMPUT(IPMTX,'MESHZ',LZ+1,2,GEOZZ)
*----
*  CHECK FUEL MIXTURES
*----
      CALL RESPFM(IPMAP,IPMTX,NX,NY,NZ,LX,LY,LZ,NFUEL,IMPX,IGEO,NCH,NB,
     1 NTOT)
      RETURN
      END
