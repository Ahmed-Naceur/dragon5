*DECK HSTUHB
      SUBROUTINE HSTUHB(IPHST,  IPEVO,  IPRINT, MAXI,   NBBTS,  NCHA,
     >                  NBUN,   IUPDC,  IUPDB,  IDCELL, IDFUEL, DENI,
     >                  MAXL,   PARAML)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To update the HISTORY data structure using the information
* provided on the BURNUP data structure.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal.
*
*Author(s): 
* G. Marleau, E. Varin
*
*Parameters: input
* IPHST   address of the \dds{history} data structure.
* IPEVO   address of the \dds{burnup} data structure.
* IPRINT  print level.
* MAXI    maximum number of isotopes.            
* NBBTS   number of depletion steps.
* NCHA    number of fuel channels.                   
* NBUN    number of bundles per channel.
* IUPDC   number of the channel to analyze.                   
* IUPDB   number of the bundle to analyze.
* IDCELL  cell identifier for each fuel bundle in each channel.
* IDFUEL  fuel type identifier for each fuel bundle in each channel.
* IPHST   pointer to the HISTORY data structure
* IPEVO   pointer to the BURNUP data structure.
* IPRINT  print level. 
* MAXI    maximum number of isotopes.
* NBBTS   number of depletion steps.
* NCHA    number of fuel channels.
* NBUN    number of bundles per channels.
* IUPDC   channel number to process or update.
* IUPDB   bundle number to process or update.
* IDCELL  list of cell identifiers. 
* IDFUEL  list of fuel type identifiers. 
* MAXL    maximum number of local parameters.                   
*
*Parameters: work
* PARAML  local parameters.
* DENI    isotopic concentrations of the isotopes 
*         on the \dds{burnup} or \dds{history} structure.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT         NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPHST,IPEVO
      INTEGER          IPRINT,MAXI,NBBTS,MAXL
      INTEGER          NCHA,NBUN,IUPDC,IUPDB
      INTEGER          IDCELL(NBUN,NCHA),IDFUEL(NBUN,NCHA)
      REAL             DENI(0:MAXI) 
      REAL             PARAML(0:MAXL,2)
*----
*  LOCAL PARAMETERS
*  CDAY = conversion of days in 10^{8} seconds
*----
      INTEGER          IOUT
      INTEGER          ILCMUP,ILCMDN
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,ILCMUP=1,ILCMDN=2,NAMSBR='HSTUHB')
      REAL             CDAY,TIMPOW(2)
      PARAMETER       (CDAY=8.64E-4)
*----
*  LOCAL VARIABLES
*---- 
      INTEGER          ILCMLN,ILCMTY 
      CHARACTER        NAMTIM*12,NAMP*12 
      INTEGER          IFT,ICT,INEWF,INEWC
      INTEGER          ITS,ISO,IOK
      REAL             BITH(3),BITB(3)
      REAL             FDENC(2),FDENF(2),FDENB(2)
      REAL             REVOL(5)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIXIH,MIXIB
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NAMIH,NAMIB
      REAL, ALLOCATABLE, DIMENSION(:) :: DEPLT
*----
*  SCRATCH STORAGE ALLOCATION
*   NAMIH   name of isotopes on the \dds{history} structure.
*   MIXIH   mixture number associated with the isotopes 
*           on the \dds{history} structure.
*   NAMIB   name of isotopes on the \dds{burnup} structure.
*   MIXIB   mixture number associated with the isotopes 
*           on the \dds{burnup} structure.
*   DEPLT   time associated with each depletion step 
*           on the \dds{burnup} structure.
*----
      ALLOCATE(NAMIH(3,0:MAXI),MIXIH(0:MAXI),NAMIB(3,0:MAXI),
     > MIXIB(0:MAXI),DEPLT(0:NBBTS))
*----
*  Initialize test flags
*  INEWF -> new fuel type flag
*           = 0 fuel type does not exists/create it
*           = 1 fuel exists but does not contain isotopes
*           = 2 fuel type exists and contains isotopes 
*  INEWC -> new cell type flag
*           = 0 cell type does not exists/create it
*           = 1 cell type exists but isotopes densities missing
*           = 2 cell type exists and contains isotope densities
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      INEWF=2
      INEWC=2
      CALL XDRSET(DENI,MAXI+1,0.0)
      CALL XDRSET(FDENC,2,0.0)
      CALL XDRSET(FDENF,2,0.0)
      CALL XDRSET(FDENB,2,0.0)
      CALL XDRSET(BITH,3,0.0)
      CALL XDRSET(BITB,3,0.0)
      CALL XDRSET(PARAML,(MAXL+1)*2,0.0)
*----
*  Read HISTORY information for cell specified
*----
      IF(IUPDC .GT. 0 .AND. IUPDB .GT. 0) THEN
*----
*  Read isotope names and mixtures on FUEL TYPE
*  if available
*----
        IFT=IDFUEL(IUPDB,IUPDC)
        WRITE(NAMP,'(A4,I8.8)') 'FUEL',IFT 
        CALL LCMLEN(IPHST,NAMP,ILCMLN,ILCMTY)
        IF(ILCMLN .EQ. 0) THEN
          INEWF=0
        ELSE
          CALL LCMSIX(IPHST,NAMP,ILCMUP)
          CALL LCMLEN(IPHST,'ISOTOPESUSED',ILCMLN,ILCMTY)
          IF(ILCMLN .GT. 0 .AND. ILCMLN .LT. 3*MAXI*4) THEN
            CALL LCMGET(IPHST,'ISOTOPESUSED',NAMIH(1,1)) 
            CALL LCMGET(IPHST,'ISOTOPESMIX',MIXIH(1))
            CALL LCMGET(IPHST,'FUELDEN-INIT',FDENF)
          ELSE
            INEWF=1
          ENDIF
          CALL LCMSIX(IPHST,NAMP,ILCMDN)
        ENDIF
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6001) 'FUEL TYPE',IFT
          IF(INEWF .EQ. 2) THEN
            IF(IPRINT .GE. 100) THEN
              WRITE(IOUT,6010) 
              WRITE(IOUT,6011) 
     >        (NAMIH(1,ISO),NAMIH(2,ISO),NAMIH(3,ISO),ISO=1,MAXI)
              WRITE(IOUT,6020) 
              WRITE(IOUT,6021) 
     >        (MIXIH(ISO),ISO=1,MAXI)
            ENDIF
          ENDIF
        ENDIF 
*----
*  Read isotope densities on CELL TYPE
*  if available
*----
        ICT=IDCELL(IUPDB,IUPDC)
        WRITE(NAMP,'(A4,I8.8)') 'CELL',ICT
        CALL LCMLEN(IPHST,NAMP,ILCMLN,ILCMTY)
        IF(ILCMLN .EQ. 0) THEN
          INEWC=0
        ELSE
          CALL LCMSIX(IPHST,NAMP,ILCMUP)
          IOK=-1
          CALL HSTGSD(IPHST ,MAXI  ,IOK   ,DENI  ,FDENC )
          INEWC=1
          IF(IOK .EQ. 0) THEN
            INEWC=2
            CALL LCMGET(IPHST,'DEPL-PARAM  ',BITH)
          ENDIF
          CALL LCMSIX(IPHST,NAMP,ILCMDN)
        ENDIF
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6001) 'CELL TYPE',ICT
          IF(INEWF .EQ. 2) THEN
            IF(IPRINT .GE. 100) THEN
              WRITE(IOUT,6100)  
              WRITE(IOUT,6110) (DENI(ISO),ISO=1,MAXI)
            ENDIF
          ENDIF
        ENDIF 
*----
*  Read isotopes names and mixtures on BURNUP
*----
        CALL LCMGET(IPEVO,'ISOTOPESUSED',NAMIB(1,1))
        CALL LCMGET(IPEVO,'ISOTOPESMIX ',MIXIB(1))
        CALL LCMGET(IPEVO,'FUELDEN-INIT',FDENB)
*----
*  Test for coherence of isotopes names and mixture
*  between HISTORY and BURNUP if fuel type contains
*  isotopes description
*----
        IF(INEWF .EQ. 2) THEN
          DO 100 ISO=1,MAXI
            IF(NAMIH(ISO,1) .NE. NAMIB(ISO,1) .OR. 
     >         NAMIH(ISO,2) .NE. NAMIB(ISO,2) .OR. 
     >         NAMIH(ISO,3) .NE. NAMIB(ISO,3) .OR. 
     >         MIXIH(ISO)   .NE. MIXIB(ISO)   ) THEN
              CALL XABORT(NAMSBR//
     >        ': Isotopes on HISTORY and BURNUP not coherent')
            ENDIF 
 100      CONTINUE
          IF(FDENF(1) .NE. FDENB(1) .OR.
     >       FDENF(2) .NE. FDENB(2) ) THEN
            CALL XABORT(NAMSBR//
     >      ': Fuel DENSITY on HISTORY and BURNUP not coherent')
          ENDIF
        ENDIF
*----
*  Read calculation types on BURNUP
*----     
        CALL LCMGET(IPEVO,'EVOLUTION-R ',REVOL)
        CALL XDRSET(DEPLT,NBBTS+1,0.0)
        CALL LCMGET(IPEVO,'DEPL-TIMES  ',DEPLT(1)) 
*----
*  Read initial burnup information (FOR FUEL TYPE)
*  and save
*----
        ITS=1 
        IF(INEWF .NE. 2 ) THEN
          BITB(1)=DEPLT(ITS)/CDAY
          IF(BITB(1) .EQ. 0.0) THEN
            WRITE(NAMTIM,'(A8,I4.4)') 'DEPL-DAT',ITS
            CALL LCMSIX(IPEVO,NAMTIM,ILCMUP)
            CALL LCMGET(IPEVO,'ISOTOPESDENS',DENI(1))
            CALL LCMGET(IPEVO,'BURNUP-IRRAD',BITB(2))
            CALL LCMSIX(IPEVO,NAMTIM,ILCMDN)
          ELSE
            CALL XABORT(NAMSBR//
     >      ': Initial DENSITY on BURNUP required') 
          ENDIF
*----
*  Save isotopes names and mixtures for FUEL type
*----
          WRITE(NAMP,'(A4,I8.8)') 'FUEL',IFT 
          CALL LCMSIX(IPHST,NAMP,ILCMUP)
          CALL LCMPUT(IPHST,'ISOTOPESUSED',3*MAXI,3,NAMIB(1,1)) 
          CALL LCMPUT(IPHST,'ISOTOPESMIX',MAXI  ,1,MIXIB(1))
          CALL LCMPUT(IPHST,'FUELDEN-INIT',2     ,2,FDENB)
          CALL LCMPUT(IPHST,'ISOTOPESDENS',MAXI  ,2,DENI(1))
          CALL LCMSIX(IPHST,NAMP,ILCMDN)
          IF(IPRINT .GE. 100) THEN
            WRITE(IOUT,6010) 
            WRITE(IOUT,6011) 
     >        (NAMIB(1,ISO),NAMIB(2,ISO),NAMIB(3,ISO),ISO=1,MAXI)
            WRITE(IOUT,6020) 
            WRITE(IOUT,6021) 
     >        (MIXIB(ISO),ISO=1,MAXI)
          ENDIF
        ELSE
          BITB(1)=DEPLT(ITS)/CDAY
          WRITE(NAMTIM,'(A8,I4.4)') 'DEPL-DAT',ITS
          CALL LCMSIX(IPEVO,NAMTIM,ILCMUP)
          CALL LCMGET(IPEVO,'BURNUP-IRRAD',BITB(2))
          CALL LCMSIX(IPEVO,NAMTIM,ILCMDN)
*----
*  Test if initial BURNUP coherent with old history
*----
          IF(INEWC .EQ. 2 ) THEN
            IF(BITB(1) .NE. BITH(1) .OR.
     >         BITB(2) .NE. BITH(2) .OR.
     >         BITB(3) .NE. BITH(3) ) THEN
              WRITE(IOUT,6200) BITH(1)
            ENDIF
          ENDIF
        ENDIF
        ITS=NBBTS 
        BITB(1)=DEPLT(ITS)/CDAY
        WRITE(NAMTIM,'(A8,I4.4)') 'DEPL-DAT',ITS
        CALL LCMSIX(IPEVO,NAMTIM,ILCMUP)
        CALL LCMGET(IPEVO,'ISOTOPESDENS',DENI(1))
        CALL LCMGET(IPEVO,'BURNUP-IRRAD',BITB(2))
        CALL LCMSIX(IPEVO,NAMTIM,ILCMDN) 
*----
*  Save power desnity and depletion time in History
*   Modif EV 04/11/09
*----
        WRITE(NAMP,'(A4,I8.8)') 'CELL',ICT
        CALL LCMSIX(IPHST,NAMP,ILCMUP)
        IOK=-2
        CALL HSTGSL(IPHST ,MAXL  ,IOK   ,TIMPOW,PARAML(0,1)) 
        IF(IOK .NE. 0) THEN 
           CALL XDRSET(PARAML,MAXL+1,0.0)
        ENDIF
        IOK=2
        TIMPOW(1)= DEPLT(NBBTS)/CDAY
        TIMPOW(2)= REVOL(5)
        CALL HSTGSL(IPHST ,MAXL  ,IOK   ,TIMPOW,PARAML(0,1)) 
*----
*  Save last densities on BURNUP
*----
        IOK=2
        CALL HSTGSD(IPHST ,MAXI  ,IOK   ,DENI  ,FDENB )
        CALL LCMPUT(IPHST,'DEPL-PARAM  ',3,2,BITB)
        CALL LCMSIX(IPHST,NAMP,ILCMDN)
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6101)  
          WRITE(IOUT,6110) (DENI(ISO),ISO=1,MAXI)
        ENDIF
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DEPLT,MIXIB,NAMIB,MIXIH,NAMIH)
*----
*  Return
*----
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(' ****** OUTPUT FROM ',A6)
 6001 FORMAT(' Contents of ',A9,1X,I8) 
 6010 FORMAT(' NAME OF ISOTOPES ')
 6011 FORMAT(10(3A4,2X))
 6020 FORMAT(' MIXTURE OF ISOTOPES ')
 6021 FORMAT(10(I12,2X))
 6100 FORMAT(' INITIAL DENSITIES')
 6101 FORMAT(' FINAL DENSITIES')
 6110 FORMAT(1P,10E14.7)
 6200 FORMAT(' Update cell densities with no chronological burnup'/
     +       ' Old time ',F6.2,' days  should be zero.'/
     +       ' Possible errors or restart case')
      END 
