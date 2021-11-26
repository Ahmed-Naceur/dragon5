*DECK RESPFM
      SUBROUTINE RESPFM(IPMAP,IPMTX,NX,NY,NZ,LX,LY,LZ,NFUEL,IMPX,IGEO,
     > NCH,NB,NTOT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* recover, check and store the fuel mixtures.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input/output
* IPMAP  pointer to fuel-map information.
* IPMTX  pointer to matex information.
* NX     number of elements along x-axis in fuel map.
* NY     number of elements along y-axis in fuel map.
* NZ     number of elements along z-axis in fuel map.
* LX     number of elements along x-axis in geometry.
* LY     number of elements along y-axis in geometry.
* LZ     number of elements along z-axis in geometry.
* NFUEL  number of fuel types.
* IMPX   printing index (=0 for no print).
* IGEO   type of geometry (=7 or =9)
*
*Parameters: output
* NCH    number of fuel channels.
* NB     number of fuel bundles per channel.
* NTOT   total number of fuel bundles.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP,IPMTX
      INTEGER NX,NY,NZ,LX,LY,LZ,NFUEL,IMPX,IGEO,NCH,NB,NTOT
      TYPE(C_PTR) JPMAP,KPMAP
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIX,FMIX,FTOT,IFLMIX
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MIX(NX*NY*NZ),FMIX(NFUEL),FTOT(NFUEL))
*----
*  COMPUTE NUMBER OF FUEL CHANNELS AND NUMBER OF FUEL BUNDLES
*----
      IF(IMPX.GT.2) WRITE(IOUT,*)'SETTING FUEL-MAP MIXTURES'
      IF((IGEO.NE.7).AND.(IGEO.NE.9))THEN
        CALL XABORT('@RESPFM: WRONG TYPE OF GEOMETRY, 3D-CARTESIAN OR '
     >  //'3D-HEXAGONAL GEOMETRY EXPECTED')
      ENDIF
      CALL XDISET(MIX,NX*NY*NZ,0)
      CALL LCMGET(IPMAP,'MIX',MIX)
      NB=0
      DO IZ=1,NZ
        DO I=1,NX*NY
          IEL=(IZ-1)*NX*NY+I
          IF(MIX(IEL).NE.0) GOTO 10
        ENDDO
        CYCLE
   10   NB=NB+1
      ENDDO
      NCH=0
      DO I=1,NX*NY
        DO IZ=1,NZ
          IEL=(IZ-1)*NX*NY+I
          IF(MIX(IEL).NE.0) GOTO 20
        ENDDO
        CYCLE
   20   NCH=NCH+1
      ENDDO
      IF(IMPX.GT.0) WRITE(6,100) NCH,NB
      ALLOCATE(IFLMIX(NCH*NB))
*----
*  COMPUTE FLMIX AND FTOT
*----
      CALL XDISET(FMIX,NFUEL,0)
      CALL LCMGET(IPMTX,'FMIX',FMIX)
      CALL XDISET(FTOT,NFUEL,0)
      CALL XDISET(IFLMIX,NCH*NB,0)
      NTOT=0
      IB=0
      DO 50 IZ=1,NZ
      DO I=1,NX*NY
         IEL=(IZ-1)*NX*NY+I
         IF(MIX(IEL).NE.0) GOTO 30
      ENDDO
      GO TO 50
   30 IB=IB+1
      IF(IB.GT.NB) CALL XABORT('@RESPFM: NB OVERFLOW.')
      ICH=0
      DO 40 I=1,NX*NY
      DO K=1,NZ
        IF(MIX((K-1)*NX*NY+I).NE.0) GOTO 35
      ENDDO
      GO TO 40
   35 IEL=(IZ-1)*NX*NY+I
      ICH=ICH+1
      IF(ICH.GT.NCH) CALL XABORT('@RESPFM: NCH OVERFLOW.')
      IFLMIX((IB-1)*NCH+ICH)=MIX(IEL)
      IF(MIX(IEL).EQ.0) GO TO 40
      DO IFUEL=1,NFUEL
        IF(MIX(IEL).EQ.FMIX(IFUEL))THEN
          FTOT(IFUEL)=FTOT(IFUEL)+1
          NTOT=NTOT+1
          IF(NTOT.GT.NCH*NB)THEN
            WRITE(IOUT,*)'@RESPFM: TOTAL NUMBER OF BUNDLES =',NCH*NB
            WRITE(IOUT,*)'@RESPFM: READ TOTAL FUEL MIXTURES ',NTOT
            CALL XABORT('@RESPFM: WRONG FUEL-MAP DEFINITION.')
          ENDIF
          GOTO 40
        ENDIF
      ENDDO
      WRITE(IOUT,*)'@RESPFM: READ FUEL MIXTURE NUMBER ',MIX(IEL)
      CALL XABORT('@RESPFM: WRONG FUEL MIXTURE NUMBER.')
   40 CONTINUE
   50 CONTINUE
      IF(IMPX.GT.0) WRITE(6,110) NTOT
*----
*  STORE FUEL MIXTURES
*----
      IF(IMPX.GT.2) WRITE(IOUT,*)'STORING FUEL MIXTURES'
*     FUEL DIRECTORIES
      CALL LCMSIX(IPMAP,' ',0)
      JPMAP=LCMLID(IPMAP,'FUEL',NFUEL)
      DO IFUEL=1,NFUEL
        KPMAP=LCMDIL(JPMAP,IFUEL)
        CALL LCMPUT(KPMAP,'MIX',1,1,FMIX(IFUEL))
        CALL LCMPUT(KPMAP,'TOT',1,1,FTOT(IFUEL))
      ENDDO
      CALL LCMPUT(IPMAP,'FLMIX',NCH*NB,1,IFLMIX)
      DEALLOCATE(IFLMIX)
*     RENUMBERING
      NMIX=0
      DO IEL=1,NX*NY*NZ
        IF(MIX(IEL).NE.0)THEN
          NMIX=NMIX+1
          MIX(IEL)=NMIX
        ENDIF
      ENDDO
      CALL LCMPUT(IPMAP,'BMIX',NX*NY*NZ,1,MIX)
*     UPDATE MATERIAL INDEX
      IF(IGEO.EQ.7)THEN
        CALL RESIND(IPMAP,IPMTX,NX,NY,NZ,LX,LY,LZ,MIX,NFUEL,IMPX)
      ELSE IF(IGEO.EQ.9)THEN
        CALL RESHID(IPMAP,IPMTX,NX,NZ,LX,LZ,MIX,NFUEL,IMPX)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FTOT,FMIX,MIX)
      RETURN
*
  100 FORMAT(/33H RESPFM: NUMBER OF FUEL CHANNELS=,I5/9X,10HNUMBER OF ,
     > 25HFUEL BUNDLES PER CHANNEL=,I5)
  110 FORMAT(9X,29HTOTAL NUMBER OF FUEL BUNDLES=,I8)
      END
