*DECK HSTREF
      SUBROUTINE HSTREF(IPHST,  IPRINT, MAXL,   NCHA,   NBUN,   MAXI,
     >                  DELTAT, POWER,  IREFUS, REFUT,  IDCELL, IDFUEL,
     >                  PARAML, DENI,   ISHUFF)
*
*----------
*
*Purpose:
* Refuel channel by performing fuel shuffling. 
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal.
*
*Author(s): 
* G. Marleau, E. Varin
*
*Parameters: input
* IPHST   address of the \dds{history} data structure.
* IPRINT  print level.
* MAXL    maximum number of local parameters.                   
* NCHA    number of fuel channels.                   
* NBUN    number of bundles per channel.            
* MAXI    maximum number of isotopes.            
* DELTAT  last character string read.
* POWER   burnup power for each fuel bundle in each channel.
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
*
*Parameters: input/output
* IDCELL  cell identifier for each fuel bundle in each channel.
* IDFUEL  fuel type identifier for each fuel bundle in each channel.
*
*Parameters: work
* PARAML  local parameters.
* DENI    isotopic concentrations.
* ISHUFF  fuel shuffling index for a channel.
*
*----------
*
      USE GANLIB
      IMPLICIT         NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPHST
      INTEGER          IPRINT,MAXL,NCHA,NBUN,MAXI
      REAL             DELTAT
      REAL             POWER(NCHA,NBUN)
      INTEGER          IREFUS(NCHA)
      REAL             REFUT(NCHA)
      INTEGER          IDCELL(NBUN,NCHA),IDFUEL(NBUN,NCHA)
      REAL             PARAML(0:MAXL,2)
      REAL             DENI(0:MAXI)
      INTEGER          ISHUFF(NBUN)
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IOUT,NTC,ILCMUP,ILCMDN
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NTC=3,ILCMUP=1,ILCMDN=2,
     >                 NAMSBR='HSTREF')
*----
*  LOCAL VARIABLES
*----
      INTEGER          IC,IB,IBS,IBO,ICT,IFT,IOK
      REAL             FDEN(2)
      REAL             TIMREF,TIMPOW(2)
      CHARACTER        NAMP*12 
*----
*  Take local paremeters after fueling
*  and store in local parameters before fueling
*  for all fuel cells
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,7000) NAMSBR
      ENDIF
      DO 100 IC=1,NCHA
        TIMREF=REFUT(IC)
        IBS=IREFUS(IC) 
        DO 110 IB=1,NBUN
          ICT=IDCELL(IB,IC)
          WRITE(NAMP,'(A4,I8.8)') 'CELL',ICT
          CALL LCMSIX(IPHST,NAMP,ILCMUP)
*----
*  Get local parameters from cell IB after refueling
*----
          IOK=-2              
          CALL HSTGSL(IPHST ,MAXL  ,IOK   ,TIMPOW,PARAML) 
          IF(IOK .NE. 0) THEN 
            CALL XDRSET(PARAML,MAXL+1,0.0)
          ENDIF
*----
*  Save local parameters from cell IB before refueling
*----         
          IOK=1
          TIMPOW(1)=TIMREF
          TIMPOW(2)=POWER(IC,IB)
          CALL HSTGSL(IPHST ,MAXL  ,IOK   ,TIMPOW,PARAML)
          IOK=2
          TIMPOW(1)=DELTAT-TIMREF
          TIMPOW(2)=POWER(IC,IB)
          CALL HSTGSL(IPHST ,MAXL  ,IOK   ,TIMPOW,PARAML)
          CALL LCMSIX(IPHST,NAMP,ILCMDN)
 110    CONTINUE
*----
*  Look for channel to refuel  
*  -> REFUT(IC) > 0.0
*  Refuel channel according to IREFUS(IC) bundle shift
*  IREFUS(IC) < 0  -> push bundles starting at I=NBUN side
*  IREFUS(IC) > 0  -> push bundles starting at I=1    side
*  For displaced fuel channels:
*     Change IDCELL to new cell identifier after displacement
*  For refuel channels 
*     Use IDCELL for channels removed from core and allocate
*     then to new fuel. 
*----
        IF(TIMREF .GT. 0.0) THEN
          IF(IPRINT .GE. 10) THEN
            WRITE(IOUT,7001) IC,IBS
          ENDIF
*----
*  Find ISHUFF(IB)=IBO
*  IBO > 0 is the position of the bundle IB before refueling
*  IBO < 0 is the free position availables for refueling
*----
          CALL XDISET(ISHUFF,NBUN,0)
          IF(IBS .GT. 0) THEN 
*----
*  push bundles starting at I=1 side
*  with +IBS > 0 bundle shifts
*  1) Displaced bundles :  position  1     -- NBUN-IBS 
*                       :  position  IBS+1 -- NBUN
*----
            IBO=0
            DO 120 IB=IBS+1,NBUN
              IBO=IBO+1
              ISHUFF(IB)=IBO
 120        CONTINUE 
*----
*  2) Inserted bundles  :  positions 1     -- IBS
*----
            IBO=NBUN-IBS
            DO 121 IB=1,IBS
              IBO=IBO+1
              ISHUFF(IB)=-IBO
 121        CONTINUE
          ELSE IF(IBS .LT. 0) THEN 
*----
*  push bundles starting at I=NBUN  side
*  with -IBS > 0 bundle shifts
*  1) Displaced bundles :  position  -IBS +1     -- NBUN
*                       :  position  1          -- NBUN+IBS
*----
            IBO=-IBS
            DO 130 IB=1,NBUN+IBS
              IBO=IBO+1
              ISHUFF(IB)=IBO
 130        CONTINUE
*----
*  2) Inserted bundles  :  positions NBUN+IBS+1 -- NBUN 
*----
            IBO=0
            DO 131 IB=NBUN+IBS+1,NBUN
              IBO=IBO+1
              ISHUFF(IB)=-IBO
 131        CONTINUE
          ENDIF
*----
*  treat refueling
*----
          DO 140 IB=1,NBUN
*----
*  Get local parameters from cell IB before refueling
*----
            ICT=IDCELL(IB,IC)
            WRITE(NAMP,'(A4,I8.8)') 'CELL',ICT
            CALL LCMSIX(IPHST,NAMP,ILCMUP)
            IOK=-1              
            CALL HSTGSL(IPHST ,MAXL  ,IOK   ,TIMPOW,PARAML) 
            CALL LCMSIX(IPHST,NAMP,ILCMDN)
*
            IBO=ISHUFF(IB)
            IF(IBO .GT. 0) THEN 
*----
*  Scan Displaced bundles
*  and save properties at old cell location
*----
              IF(IPRINT .GE. 10) THEN
                WRITE(IOUT,7010) IBO,IB
              ENDIF
*----
*  Save local parameters to cell IBO after refueling
*----         
              ICT=IDCELL(IBO,IC)
              WRITE(NAMP,'(A4,I8.8)') 'CELL',ICT
              CALL LCMSIX(IPHST,NAMP,ILCMUP)
              IOK=2
              TIMPOW(1)=DELTAT-TIMREF
              TIMPOW(2)=POWER(IC,IB)
              CALL HSTGSL(IPHST ,MAXL  ,IOK   ,TIMPOW,PARAML) 
              CALL LCMSIX(IPHST,NAMP,ILCMDN)
*----
*  Save in ISHUFF IDCELL for IBO
*----         
           ELSEIF(IBO .LT. 0) THEN
*----
*  Scan inserted fuel
*  and save properties at reused cell location
*----
              IF(IPRINT .GE. 10) THEN
                WRITE(IOUT,7011) IB
              ENDIF
              IBO=-IBO
*----
*  Get initial density for fuel type
*----         
              IFT=IDFUEL(IB,IC)
              WRITE(NAMP,'(A4,I8.8)') 'FUEL',IFT
              CALL LCMSIX(IPHST,NAMP,ILCMUP)
              IOK=-1              
              CALL HSTGSD(IPHST ,MAXI  ,IOK   ,DENI  ,FDEN   )
              CALL LCMSIX(IPHST,NAMP,ILCMDN)
*----
*  Save local parameters before and after refueling
*  from cell IBO before refueling
*  Save fuel density for fuel type
*----         
              ICT=IDCELL(IBO,IC)
              WRITE(NAMP,'(A4,I8.8)') 'CELL',ICT
              CALL LCMSIX(IPHST,NAMP,ILCMUP)
              IOK=1
              TIMPOW(1)=0.0
              TIMPOW(2)=POWER(IC,IB)
              CALL HSTGSL(IPHST ,MAXL  ,IOK   ,TIMPOW,PARAML) 
              IOK=2
              TIMPOW(1)=DELTAT-TIMREF
              TIMPOW(2)=POWER(IC,IB)
              CALL HSTGSL(IPHST ,MAXL  ,IOK   ,TIMPOW,PARAML) 
              IOK=2
              CALL HSTGSD(IPHST ,MAXI  ,IOK   ,DENI  ,FDEN  )
              CALL LCMSIX(IPHST,NAMP,ILCMDN)
*----
*  Save in ISHUFF IDCELL for IBO
*----         
            ENDIF
              ISHUFF(IB)=ICT
 140      CONTINUE
*----
*  Redefine IDCELL for new spatial location
*  of cells after refueling
*  Here assume that bundles are replaced
*  with fuels of the same type
*----
          DO 160 IB=1,NBUN
            IDCELL(IB,IC)=ISHUFF(IB)
 160      CONTINUE
        ENDIF
 100  CONTINUE 
*----
*  Save IDCELL and IDFUEL since they were updated
*----
      CALL LCMPUT(IPHST,'CELLID      ',NBUN*NCHA,1,IDCELL)
      CALL LCMPUT(IPHST,'FUELID      ',NBUN*NCHA,1,IDFUEL)
*----
*  Return
*----
      RETURN
*----
* Format
*----
 7000 FORMAT(' ***** OUTPUT FROM ',A6,' *****')
 7001 FORMAT(' Refueling channel ',I8, ' with ',I8,' bundle shifts')
 7010 FORMAT(10X,' Fuel bundle ',I8,' displaced to position ',I8)
 7011 FORMAT(10X,' Fresh fuel inserted at position ',I8)
      END 
