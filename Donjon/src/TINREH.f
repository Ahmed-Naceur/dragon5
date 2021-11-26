*DECK TINREH
      SUBROUTINE TINREH(IPRES,IPMIC,IPMIC2,IPMIC3,NCH,NK,NH,NZ,NREG,
     +         NAMCHA,NS,MS,WINT,MIX,IHN,BS,PS,ISFT,POW,MAXS,NSS,IND,
     +         IPRT,KRF,LMIC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Refuel a channel according to a refuelling mode in hexagonal
* geometry.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
*
*Author(s): 
* E. Varin, M. Guyot and A. Hebert
*
*Parameters: input/output
* IPRES  Adress of the map Linked_List or XSM file.
* IPMIC  Adress of the L_LIBRARY in creation mode.
* IPMIC2 Adress of the fuel-map L_LIBRARY in read-only mode.
* IPMIC3 Adress of the L_LIBRARY in read-only mode, containing new 
*        fuel properties.
* NCH    Number of channels
* NK     Number of bundles per channel
* NH     Number of hexagons in the plane
* NZ     Number axial planes
* NREG   Number of regions in fuel map geometry
* NAMCHA Name of the channel to refuel
* NS     Number of bundles inserted
* MS     Old maximum of shift + 1.
* MIX    Fuel map bundle index
* IHN    Name of the channel according to the hexagonal position
* POW    Power distribution.
* INDEX  Fuel type indice
* IND    Fuel type indice in the channel to refuel
* MAXS   Maximum number of power shift
* IPRT   Flag for printing level
* KRF    Type of refueling
* LMIC   =.true. for a micro-refueling
*
*Parameters:
* WINT   
* BS     
* PS     
* ISFT   
* NSS    
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPRES,IPMIC,IPMIC2,IPMIC3
      INTEGER     NCH,NK,NH,NZ,NS,NREG,ILONG,ITYP,IH,IPRT,MS,IS,MAXS,
     1            KS,NSS(NCH),NNS
      REAL        WINT(NCH,NK),BS(NCH,NK,MS),PS(NCH,NK,MS),POW(NCH,NK)
      INTEGER     MIX(NREG),IHN(2,NH),ISFT(NCH,NK),IND(*)
      LOGICAL     LMIC
      REAL        TMPDAY(3)
      CHARACTER   NAMCHA*8
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,IOUT=6)
      INTEGER   ISTATE(NSTATE),I,J
      CHARACTER HNAM*8,CS*2,HSMG*131
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NW,NW2,NWU,ISONA,ISOMI,ISHF,
     + ICHMAP
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDEX,IWORK
      REAL, ALLOCATABLE, DIMENSION(:) :: DENIS,NDENS
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WORK
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: WORKS
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASK,MASKL
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NW(NK),NW2(NK),INDEX(NCH,NK),IWORK(NK,2),ICHMAP(NH))
      ALLOCATE(WORK(NK,2),WORKS(NK,MS,2))
*----
*  RECOVER SHIFT VECTOR
*----
      CALL LCMLEN(IPRES,'ISHIFT',ILS,ITYLCM)
      IF(ILS.NE.0) THEN
        CALL LCMGET(IPRES,'ISHIFT',ISFT(1,1))
        DO 18 I=1,NK
         DO 17 J=1,NCH
            MAXS=MAX(MAXS,ISFT(J,I))
 17      CONTINUE
 18     CONTINUE
      ELSE
        MAXS=0
        DO 115 I=1,NK
          DO 15 J=1,NCH
            ISFT(J,I) = 0
 15       CONTINUE
115     CONTINUE
      ENDIF
      DO 1 I=1,NK
       DO 2 J=1,NCH
         WINT(J,I) = 0.0
         DO 3 K=1,MS
           BS(J,I,K)=0.0
           PS(J,I,K)=0.0
  3      CONTINUE
  2    CONTINUE
  1   CONTINUE
*----
*  RECOVER FUEL BURNUPS
*----
      CALL LCMLEN(IPRES,'BURN-INST',ILONG,ITYP)
      IF(ILONG.EQ.0) THEN
         CALL XABORT('@TINREH: INITIAL BURNUP REQUIRED')
      ENDIF
      CALL LCMGET(IPRES,'BURN-INST',WINT)
*----
*  RECOVER FUEL INDEX
*----
      CALL LCMLEN(IPRES,'FLMIX',ILONG,ITYP)
      IF(ILONG.NE.0) THEN
         CALL LCMGET(IPRES,'FLMIX',INDEX)
      ELSE
         CALL XABORT('@TINREH: FLMIX ARE REQUIRED')
      ENDIF

      IF(MAXS.GT.0) THEN
        DO 16 IS=1,MAXS
          WRITE (CS,'(I2)') IS
          CALL LCMGET(IPRES,'BSHIFT'//CS,BS(1,1,IS))
          CALL LCMGET(IPRES,'PSHIFT'//CS,PS(1,1,IS))
 16     CONTINUE
      ENDIF
*----
*  SET THE CHANNEL INDEX MAP
*----
      CALL LCMSIX(IPRES,' ',0)
      CALL LCMGET(IPRES,'BMIX',MIX)
      CALL XDISET(ICHMAP,NH,0)
      ICH=0
      DO 25 IH=1,NH
      DO 23 IZ=1,NZ
      IF(MIX((IZ-1)*NH+IH).NE.0) GO TO 24
  23  CONTINUE
      GO TO 25
  24  ICH=ICH+1
      ICHMAP(IH)=ICH
  25  CONTINUE
      IF(ICH.NE.NCH) CALL XABORT('@TINREH: INVALID NUMBER OF CHANNELS')
*----
*  SEARCH FOR THE CHANNEL NUMBER FROM ITS NAME
*----
      IH = 0
      DO 10 I=1,NH
        WRITE(HNAM,'(2A4)') IHN(1,I),IHN(2,I)
        IF (HNAM.EQ.NAMCHA) THEN
           IH = I
           GOTO 21
        ENDIF
  10  CONTINUE
      WRITE(HSMG,'(26H@TINREH: NO CHANNEL NAMED ,A8,12H IN FUELMAP.)')
     + NAMCHA
      CALL XABORT(HSMG)

  21  ICH=ICHMAP(IH)
      IF(ICH.EQ.0) THEN
        WRITE(6,'(13H @TINREH: IH=,I6)') IH
        WRITE(HSMG,'(23H@TINREH: CHANNEL NAMED ,A8,13H HAS NO FUEL.)')
     +  NAMCHA
        CALL XABORT(HSMG)
      ENDIF
      IF(NSS(ICH).NE.0) THEN
        IF(ABS(NSS(ICH)).NE.ABS(NS)) THEN
          WRITE(6,'(14H @TINREH: ICH=,I6,5H NSS=,I6,4H NS=,I6)') ICH,
     +    NSS(ICH),NS
          CALL XABORT('@TINREH: WRONG REFUELING SCHEME')
        ENDIF
        NS = NSS(ICH)
      ENDIF
      IF( IPRT.GT.3 )THEN
        WRITE(6,*) ' REFUELING CHANNEL ',NAMCHA,' IH ',IH
        WRITE(6,*) ' REFUELING CHANNEL ',ICH,' SCHEME ',NS
        WRITE(6,*) ' INITIAL BURNUP ',(WINT(ICH,I),I=1,NK)
      ENDIF

      NNS = ABS(NS)
      CALL TINFL(NNS,NW,NW2,NK)

      II=0
      DO 30 K=1,NK
         KK = K
         IF (NS.LT.0) THEN
            KK = NK - K + 1
         ENDIF
         KA = NW(K)
*----
*  INSERTION OF A NEW BUNDLE OR REPOSITIONNING
*----
         IF (KA.EQ.0) THEN
            II=II+1
            WORK(KK,1) = 0.0
            IWORK(KK,1)=0
            IF( KRF.EQ.1 )THEN
              IWORK(KK,2)=INDEX(ICH,KK)
            ELSE
              IWORK(KK,2)=IND(II)
            ENDIF
            IF(MAXS.GT.0) THEN
              DO 39 IS=1,MAXS
                WORKS(KK,IS,1) = 0.0
                WORKS(KK,IS,2) = 0.0
  39          CONTINUE
            ENDIF
         ELSE
            IF (NS.LT.0) THEN
               KA = NK - KA + 1
            ENDIF
            WORK(KK,1) = WINT(ICH,KA)
            WORK(KK,2) = POW(ICH,KA)
            IWORK(KK,1)= ISFT(ICH,KA)
            IWORK(KK,2)= INDEX(ICH,KA)
            IF(MAXS.GT.0) THEN
              DO 19 IS=1,MAXS
                WORKS(KK,IS,1) = BS(ICH,KA,IS)
                WORKS(KK,IS,2) = PS(ICH,KA,IS)
  19          CONTINUE
            ENDIF
         ENDIF
  30  CONTINUE

      DO 40 K=1,NK
        WINT(ICH,K) = WORK(K,1)
        POW(ICH,K) = WORK(K,2)
        ISFT(ICH,K) = IWORK(K,1)
        INDEX(ICH,K) = IWORK(K,2)
        IF(MAXS.GT.0) THEN
          DO 22 IS=1,MAXS
            BS(ICH,K,IS)=WORKS(K,IS,1)
            PS(ICH,K,IS)=WORKS(K,IS,2)
  22      CONTINUE
        ENDIF
        IF(WORK(K,1).NE.0.0) THEN
          KS=ISFT(ICH,K)+1
          BS(ICH,K,KS)=WINT(ICH,K)
          PS(ICH,K,KS)=WORK(K,2)
          ISFT(ICH,K)=KS
        ENDIF
  40  CONTINUE

      MAXS=0
      DO 112 I=1,NK
        DO 12 J=1,NCH
          MAXS=MAX(MAXS,ISFT(J,I))
  12    CONTINUE
 112  CONTINUE
*----
* CALL THE SUBROUTINE FOR A MICROSCOPIC REFUEL
*----
      IF(LMIC) THEN
        CALL XDISET(ISTATE,NSTATE,0)
        CALL LCMGET(IPMIC2,'STATE-VECTOR',ISTATE)
        NISO=ISTATE(2)
        NDEP=ISTATE(12)
        IF(NDEP.NE.NK*NCH) CALL XABORT('@TINREH: WRONG NUMBER OF '
     +  //'DEPLETING MIXTURES IN THE LIBRARY.')
        CALL XDISET(ISTATE,NSTATE,0)
        CALL LCMGET(IPMIC3,'STATE-VECTOR',ISTATE)
        NISO2=ISTATE(2)
        ALLOCATE(NWU(NK))
        DO I=1,NK
          IF(NS.GT.0) THEN
            NWU(I)=NW(I)
          ELSE
            NWU(I)=NW2(I)
          ENDIF
        ENDDO
        ALLOCATE(NDENS(NISO),ISHF(NK))
        CALL TINMIC(IPMIC,IPMIC2,IPMIC3,NK,NCH,NWU,ICH,NISO,NISO2,
     1  IWORK,ISHF,NDENS)
        CALL LCMPUT(IPMIC,'ISOTOPESDENS',NISO,2,NDENS)
*----
*  COMPUTE THE MACROSCOPIC X-SECTIONS
*----
        CALL LCMGET(IPMIC,'STATE-VECTOR',ISTATE)
        MAXMIX=ISTATE(1)
        NBISO=ISTATE(2)
        NGRP=ISTATE(3)
        ALLOCATE(MASK(MAXMIX),MASKL(NGRP))
        ALLOCATE(ISONA(3*NBISO),ISOMI(NBISO),DENIS(NBISO))
        CALL LCMGET(IPMIC,'ISOTOPESUSED',ISONA)
        CALL LCMGET(IPMIC,'ISOTOPESMIX',ISOMI)
        CALL LCMGET(IPMIC,'ISOTOPESDENS',DENIS)
        CALL XDLSET(MASK,MAXMIX,.FALSE.)
        CALL XDLSET(MASKL,NGRP,.TRUE.)
        DO 13 I=1,NBISO
          IBM=ISOMI(I)
          MASK(IBM)=.TRUE.
  13    CONTINUE
        ITSTMP=0
        TMPDAY(1)=0.0
        TMPDAY(2)=0.0
        TMPDAY(3)=0.0
*       COMPUTATION OF THE MACROSCOPIC XS
        CALL LIBMIX(IPMIC,MAXMIX,NGRP,NBISO,ISONA,ISOMI,DENIS,MASK,
     1  MASKL,ITSTMP,TMPDAY)
        DEALLOCATE(DENIS,ISOMI,ISONA,MASKL,MASK)
        DEALLOCATE(NWU,NDENS,ISHF)
      ENDIF

      IF( IPRT.GT.3 )THEN
        WRITE(6,*) ' SHIFTING BURNUP ',(WINT(ICH,I),I=1,NK)
      ENDIF

      CALL LCMSIX(IPRES,' ',0)
      IF(IPRT.GT.3) WRITE(6,*) ' REFUELLING TYPE DIRECT OR HOMOGENOUS'
      CALL LCMPUT(IPRES,'BURN-INST',NCH*NK,2,WINT(1,1))
      CALL LCMPUT(IPRES,'FLMIX',NCH*NK,1,INDEX(1,1))
      CALL LCMPUT(IPRES,'ISHIFT',NCH*NK,1,ISFT(1,1))

      IF(MAXS.GT.0) THEN
        DO 14 IS=1,MAXS
          WRITE (CS,'(I2)') IS
          CALL LCMPUT(IPRES,'BSHIFT'//CS,NCH*NK,2,BS(1,1,IS))
          CALL LCMPUT(IPRES,'PSHIFT'//CS,NCH*NK,2,PS(1,1,IS))
 14    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORKS,WORK,ICHMAP,IWORK,INDEX,NW2,NW)
      RETURN
      END
