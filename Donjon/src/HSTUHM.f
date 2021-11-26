*DECK HSTUHM
      SUBROUTINE HSTUHM(IPHST, IPMAP, IPRINT, MAXL,   NCHA,   NBUN,   
     >                  MAXI, POWER,  BURNP, IREFUS, REFUT, BUNLEN,
     >                  IDCELL, IDFUEL, PARAML, DENI)
*
*----------
*
*Purpose:
* Store bundle power and depletion time in History 
* Refuel channel by performing fuel shuffling. 
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal.
*
*Author(s): 
* G. Marleau, E. Varin
*
*Parameters:
* IPHST   address of the \dds{history} data structure.
* IPMAP   address of the \dds{map} data structure.
* IPRINT  print level.
* MAXL    maximum number of local parameters.                   
* NCHA    number of fuel channels.                   
* NBUN    number of bundles per channel.            
* MAXI    maximum number of isotopes.            
* NBFUEL  number of fueled channels.            
* DELTAT  last character string read.
* POWER   power for each fuel bundle in each channel.
* BURNP   burnup for each fuel bundle in each channel.
* IREFUS  refueling strategy for each channel.
*         refueling strategy for each channel.
*         A channel is refueled using a NBS bundle 
*         shift procedure if IREFUS(I)=NBS. 
*         In the case where NBS $>$ 0,
*         bundles 1 to NBUN-NBS are displaced to position NBS+1 to
*         NBUN while locations 1 to NBS are filled with new fuel. 
*         In the case where NBS $<$ 0,
*         bundles -NBS+1 to NBUN are displaced to position 1 to
*         NBUN+NBS while locations NBUN+NBS+1 to NBUN are filled 
*         with new fuel.
* REFUT   refueling time for each channel.
* BUNLEN  length (cm) of a bundle.            
*
*Parameters: input/output
* IDCELL  cell identifier for each fuel bundle in each channel.
* IDFUEL  fuel type identifier for each fuel bundle in each channel.
*
*Parameters: work
* PARAML  local parameters.
* DENI    isotopic concentrations.
*
*----------
*
      USE GANLIB
      IMPLICIT         NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPHST,IPMAP
      INTEGER          IPRINT,MAXL,NCHA,NBUN,MAXI
      INTEGER          NBFUEL
      REAL             DELTAT, BUNLEN
      REAL             POWER(NCHA,NBUN),BURNP(NCHA,NBUN)
      INTEGER          IREFUS(NCHA)
      REAL             REFUT(NCHA)
      INTEGER          IDCELL(NBUN,NCHA),IDFUEL(NBUN,NCHA)
      REAL             PARAML(0:MAXL,2)
      REAL             DENI(0:MAXI)
*----
*  LOCAL PARAMETERS
*----
      TYPE(C_PTR)      JPMAP,KPMAP
      INTEGER          IOUT,NTC,ILCMUP,ILCMDN
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NTC=3,ILCMUP=1,ILCMDN=2,
     >                 NAMSBR='HSTUHM')
*----
*  LOCAL VARIABLES
*----
      INTEGER          ILONG,ITYP
      INTEGER          IUPDC,IUPDB
      INTEGER          IOK,ICT,IFT
      REAL             FDEN(2),RWEIGHT,WEIGHT,TIME
      REAL             TIMPOW(2)
      CHARACTER        NAMP*12
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISHUFF
*----
*  SCRATCH STORAGE ALLOCATION
*   ISHUFF  fuel shuffling index for a channel.
*----
      ALLOCATE(ISHUFF(NBUN))
*
      NBFUEL=0
      CALL XDRSET(PARAML,(MAXL+1)*2,0.0)
      DELTAT = 0.0
      TIME=0.0
*----
*  Get information in IPMAP
*----   
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,7000) NAMSBR
      ENDIF
*----   
      CALL HSTGMA(IPMAP ,NCHA ,NBUN   ,DELTAT,
     >            POWER ,BURNP,IREFUS ,REFUT ,NBFUEL)
*----
      DO 10 IUPDC=1,NCHA
        DO 11 IUPDB=1,NBUN
*
        IF(IDCELL(IUPDB,IUPDC) .LE. 0) THEN
          IDCELL(IUPDB,IUPDC)= IUPDC + (IUPDB - 1)*NCHA
          IDFUEL(IUPDB,IUPDC)=1 
        ENDIF
        IFT=IDFUEL(IUPDB,IUPDC)
        WRITE(NAMP,'(A4,I8.8)') 'FUEL',IFT
        CALL LCMSIX(IPHST,NAMP,ILCMUP)
        CALL LCMSIX(IPHST,NAMP,ILCMDN)
*----
*   store power and time after refueling
*   for all fuel cells
*----
        ICT=IDCELL(IUPDB,IUPDC)
        WRITE(NAMP,'(A4,I8.8)') 'CELL',ICT
        CALL LCMSIX(IPHST,NAMP,ILCMUP)
*----
*  Get fuel density or weight
*----         
        IOK=-1              
        RWEIGHT= 1.
        CALL HSTGSD(IPHST ,MAXI  ,IOK   ,DENI  ,FDEN   )
        IF(IOK.LT.0) THEN
          JPMAP=LCMGID(IPMAP,'FUEL')
          KPMAP=LCMGIL(JPMAP,1)
          CALL LCMLEN(KPMAP,'WEIGHT',ILONG,ITYP)
          IF (ILONG .EQ.0) 
     +     CALL XABORT(NAMSBR//' FUEL WEIGHT MUST BE SPECIFIED IN MAP')
          CALL LCMGET(KPMAP,'WEIGHT',WEIGHT)
          RWEIGHT= 1./WEIGHT
        ELSEIF(IOK.EQ.0) THEN
          RWEIGHT=1000.0/(FDEN(2)*BUNLEN)
        ENDIF
        IF(IPRINT .GE. 10) THEN
           WRITE(IOUT,7003) NAMP,POWER(IUPDC,IUPDB),BURNP(IUPDC,IUPDB),
     +     WEIGHT
        ENDIF
        POWER(IUPDC,IUPDB)=POWER(IUPDC,IUPDB)*RWEIGHT
        IF(DELTAT.EQ.0.0) THEN
          TIME = BURNP(IUPDC,IUPDB)/POWER(IUPDC,IUPDB)
        ELSE
          TIME = DELTAT
        ENDIF
*----
*  Save local parameters from cell IB after refueling
*----         
        IOK=-2
        CALL HSTGSL(IPHST ,MAXL  ,IOK   ,TIMPOW,PARAML)
        IF(IOK .NE. 0) THEN 
          CALL XDRSET(PARAML,MAXL+1,0.0)
        ENDIF
*-----
        IOK=2
        TIMPOW(1)=TIME
        TIMPOW(2)=POWER(IUPDC,IUPDB)
        CALL HSTGSL(IPHST ,MAXL  ,IOK   ,TIMPOW,PARAML)
        CALL LCMSIX(IPHST,NAMP,ILCMDN)
        IF(IPRINT .GE. 10) THEN
           WRITE(IOUT,7002) NAMP, POWER(IUPDC,IUPDB), TIME
        ENDIF
 11     CONTINUE
 10   CONTINUE 
**
      IF(NBFUEL .GT. 0) THEN
        CALL HSTREF(IPHST ,IPRINT,MAXL  ,NCHA ,NBUN ,MAXI ,
     >              DELTAT, POWER ,IREFUS,REFUT,
     >              IDCELL,IDFUEL,PARAML, DENI ,ISHUFF)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ISHUFF)
*----
*  Return
*----
      RETURN
*----
* Format
*----
 7000 FORMAT(' ***** OUTPUT FROM ',A6,' *****')
 7002 FORMAT(' Fuel cell  ',A12, ' with ',F12.4,' kW/kg ',
     >       F10.2,' days ')
 7003 FORMAT(' Fuel cell  ',A12, ' with ',F12.4,' kW/kg ',F12.3,
     >       ' kWd/kg  ',F12.3,' kg ')
      END 
