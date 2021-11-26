*DECK TINSHH
      SUBROUTINE TINSHH(IPRES,NCH,NK,NH,NZ,NREG,MS,NAMCHA,NAMCH2,
     + WINT,MIX,BS,PS,ISFT,IHN,IPRT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute new burnup values per channel after shuffling of two
* channels in hexagonal geometry.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
*
*Author(s): 
* E. Varin, M. Guyot and A. Hebert
*
*Parameters: input/output
* IPRES  ßAdress of the map Linked_List or XSM file.
* NAMCHA Name of the channel to refuel
* NAMCH2 Name of the channel to refuel
* NS     Number of bundles inserted
* MIX    Fuel map bundle index
* MS     Maximum number of power shift
*
*Parameters:
* NCH      
* NK       
* NH       
* NZ       
* NREG     
* WINT     
* BS       
* PS       
* ISFT     
* IHN      
* IPRT     
*      
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPRES
      INTEGER     NCH,NK,NH,NZ,NREG,ILONG,ITYP,IPRT,ICH1,ICH2,ILS,
     1            ITYLCM,IS,MAXS,MS
      REAL        WINT(NCH,NK),BS(NCH,NK,MS),PS(NCH,NK,MS)
      INTEGER     MIX(NREG),IHN(2,NH),ISFT(NCH,NK)
      CHARACTER   NAMCHA*8,NAMCH2*8
*----
*  LOCAL VARIABLES
*----
      INTEGER ICH,I,J,IZ,IH
      CHARACTER HNAM*8,CS*2,HSMG*131
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ICHMAP
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDEX
      REAL, ALLOCATABLE, DIMENSION(:,:) :: POOL
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ICHMAP(NH),INDEX(NCH,NK),POOL(NCH,NK))
*----
*  RECOVER INFORMATIONS FROM FUEL MAP OBJECT
*----
      DO 1 I=1,NK
       DO 2 J=1,NCH
         WINT(J,I) = 0.0
         ISFT(J,I) = 0
         POOL(J,I) = 0.0
         DO 3 IS=1,MS
           BS(J,I,IS)=0.0
           PS(J,I,IS)=0.0
  3      CONTINUE
  2   CONTINUE
  1   CONTINUE
*----
*  RECOVER FUEL BURNUPS
*----
      CALL LCMLEN(IPRES,'BURN-INST',ILONG,ITYP)
      IF(ILONG.EQ.0) THEN
         CALL XABORT('@TINSHH: INITIAL BURNUP REQUIRED')
      ENDIF
      CALL LCMGET(IPRES,'BURN-INST',WINT(1,1))
*----
*  RECOVER FUEL INDEX
*----
      CALL LCMLEN(IPRES,'FLMIX',ILONG,ITYP)
      IF(ILONG.NE.0) THEN
         CALL LCMGET(IPRES,'FLMIX',INDEX)
      ELSE
         CALL XABORT('@TINSHH: FLMIX ARE REQUIRED')
      ENDIF
*----
*  RECOVER SHIFT VECTOR
*----
      MAXS=0
      CALL LCMLEN(IPRES,'ISHIFT',ILS,ITYLCM)
      IF(ILS.NE.0) THEN
        CALL LCMGET(IPRES,'ISHIFT',ISFT(1,1))
        DO 16 I=1,NK
         DO 15 J=1,NCH
            MAXS=MAX(MAXS,ISFT(J,I))
 15     CONTINUE
 16     CONTINUE
      ELSE
        MAXS=0
      ENDIF

      IF(MAXS.GT.0) THEN
        DO 17 IS=1,MAXS
          WRITE (CS,'(I2)') IS
          CALL LCMGET(IPRES,'BSHIFT'//CS,BS(1,1,IS))
          CALL LCMGET(IPRES,'PSHIFT'//CS,PS(1,1,IS))
 17     CONTINUE
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
      IF(ICH.NE.NCH) CALL XABORT('@TINSHH: INVALID NUMBER OF CHANNELS')
*----
*  SEARCH FOR CHANNEL NUMBER TO MOVE
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

  21  ICH1 = ICHMAP(IH)
      IF(ICH1.EQ.0) THEN
        WRITE(6,'(13H @TINSHH: IH=,I6)') IH
        WRITE(HSMG,'(23H@TINSHH: CHANNEL NAMED ,A4,13H HAS NO FUEL.)')
     +  NAMCHA
        CALL XABORT(HSMG)
      ENDIF
      IF(IPRT.GT.3) THEN
        WRITE(6,*)
        WRITE(6,*) ' SHUFFLING CHANNEL ',NAMCHA,ICH1
        WRITE(6,*) ' BEFORE ',NAMCHA,(WINT(ICH1,I),I=1,NK)
      ENDIF
*----
*  SEARCH FOR CHANNEL NUMBER WHERE TO MOVE
*----
      IF(NAMCH2.NE.'POOL') THEN
        IH = 0
        DO 40 I=1,NH
          WRITE(HNAM,'(2A4)') IHN(1,I),IHN(2,I)
          IF (HNAM.EQ.NAMCH2) THEN
             IH = I
             GOTO 41
          ENDIF
  40    CONTINUE
        WRITE(HSMG,'(26H@TINREH: NO CHANNEL NAMED ,A8,12H IN FUELMAP.)')
     +  NAMCHA
        CALL XABORT(HSMG)

  41    ICH2 = ICHMAP(IH)
        IF(ICH2.EQ.0) CALL XABORT('@TINSHH: WRONG CHANNEL NAME')
        IF(IPRT.GT.3) THEN
          WRITE(6,*)
          WRITE(6,*) ' SHUFFLING CHANNEL ',NAMCH2,ICH2
          WRITE(6,*) ' BEFORE ',NAMCH2,(WINT(ICH2,I),I=1,NK)
        ENDIF
*----
*  SHUFFLING
*----
        DO 50 I=1,NK
          IF(WINT(ICH2,I).NE.0.0) THEN
             WRITE(6,*) ' BURNUP ',WINT(ICH2,I)
             CALL XABORT('@TINSHH: WRONG POSITION TO SHUFFLE, '
     +          //'CHANNEL NOT EMPTY')
          ENDIF
          WINT(ICH2,I) = WINT(ICH1,I)
          WINT(ICH1,I) = 0.0
          ISFT(ICH2,I) = ISFT(ICH1,I)
          ISFT(ICH1,I) = 0
          INDEX(ICH2,I) = INDEX(ICH1,I)
          IF(MAXS.GT.0) THEN
             DO 56 IS=1,MAXS
               BS(ICH2,I,IS) = BS(ICH1,I,IS)
               PS(ICH2,I,IS) = PS(ICH1,I,IS)
               BS(ICH1,I,IS) = 0.0
               PS(ICH1,I,IS) = 0.0
  56         CONTINUE
          ENDIF
  50    CONTINUE
        IF(IPRT.GT.3) THEN
          WRITE(6,*)
          WRITE(6,*) ' AFTER ',NAMCH2,(WINT(ICH2,I),I=1,NK)
        ENDIF
      ELSE
        WRITE(6,*) ' CHANNEL TO POOL '
*----
*  RECOVER DISCHARGED FUEL BURNUPS
*----
        CALL LCMLEN(IPRES,'BURN-POOL',ILONG,ITYP)
        IF(ILONG.NE.0) THEN
           CALL LCMGET(IPRES,'BURN-POOL',POOL(1,1))
        ENDIF
        DO 51 I=1,NK
          POOL(ICH1,I) = WINT(ICH1,I)
          WINT(ICH1,I) = 0.0
  51    CONTINUE
        CALL LCMPUT(IPRES,'BURN-POOL',NCH*NK,2,POOL(1,1))
      ENDIF
      IF(IPRT.GT.3)
     +    WRITE(6,*) ' AFTER ',NAMCHA,(WINT(ICH1,I),I=1,NK)
      CALL LCMSIX(IPRES,' ',0)
      CALL LCMPUT(IPRES,'BURN-INST',NCH*NK,2,WINT(1,1))
      CALL LCMPUT(IPRES,'FLMIX',NCH*NK,1,INDEX(1,1))
      CALL LCMPUT(IPRES,'ISHIFT',NCH*NK,1,ISFT(1,1))
      IF(MAXS.GT.0) THEN
        DO 53 IS=1,MAXS
          WRITE (CS,'(I2)') IS
          CALL LCMPUT(IPRES,'BSHIFT'//CS,NCH*NK,2,BS(1,1,IS))
          CALL LCMPUT(IPRES,'PSHIFT'//CS,NCH*NK,2,PS(1,1,IS))
 53     CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(POOL,INDEX,ICHMAP)
      RETURN
      END
