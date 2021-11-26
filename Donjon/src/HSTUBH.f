*DECK HSTUBH                                 
      SUBROUTINE HSTUBH(IPEVO,  IPHST,  IPRINT, MAXI,   NBBTS,   NCHA,
     >                  NBUN,   IUPDC,  IUPDB,  IDCELL, IDFUEL,  DENI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To update the BURNUP data structure using the information
* provided on the HISTORY data structure.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal.
*
*Author(s): 
* G. Marleau
*
*Parameters: input
* IPEVO   address of the \dds{burnup} data structure.
* IPHST   address of the \dds{history} data structure.
* IPRINT  print level.
* MAXI    maximum number of isotopes.            
* NBBTS   number of depletion steps.
* NCHA    number of fuel channels.                   
* NBUN    number of bundles per channel.
* IUPDC   number of the channel to analyze.                   
* IUPDB   number of the bundle to analyze.
* IDCELL  cell identifier for each fuel bundle in each channel.
* IDFUEL  fuel type identifier for each fuel bundle in each channel.
*
*Parameters: output
* NAMIH   name of isotopes on the \dds{history} 
*         or \dds{burnup} structure.
* MIXIH   mixture number associated with the isotopes 
*         on the \dds{history} or \dds{burnup} structure.
* DENI    isotopic concentrations of the isotopes 
*         on the \dds{history} or \dds{burnup} structure.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT         NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPHST,IPEVO
      INTEGER          IPRINT,MAXI,NBBTS
      INTEGER          NCHA,NBUN,IUPDC,IUPDB
      INTEGER          IDCELL(NBUN,NCHA),IDFUEL(NBUN,NCHA)
      REAL             DENI(0:MAXI) 
*----
*  LOCAL PARAMETERS
*  CDAY = conversion of days in 10^{8} seconds
*----
      INTEGER          IOUT
      INTEGER          ILCMUP,ILCMDN
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,ILCMUP=1,ILCMDN=2,NAMSBR='HSTUBH')
      REAL             CDAY
      PARAMETER       (CDAY=8.64E-4)
*----
*  LOCAL VARIABLES
*---- 
      INTEGER          ILCMLN,ILCMTY 
      CHARACTER        NAMTIM*12,NAMP*12 
      INTEGER          IFT,ICT
      INTEGER          ITS,ISO,IOK
      REAL             BITH(3)
      REAL             FDENC(2)
      REAL             FLXNOR,DELTA(2)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIXIH
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NAMIH
      REAL, ALLOCATABLE, DIMENSION(:) :: DEPLT
*----
*  SCRATCH STORAGE ALLOCATION
*   NAMIH   name of isotopes on the \dds{history} structure.
*   MIXIH   mixture number associated with the isotopes 
*           on the \dds{history} structure.
*   DEPLT   time associated with each depletion step 
*           on the \dds{burnup} structure.
*----
      ALLOCATE(NAMIH(3,0:MAXI),MIXIH(0:MAXI),DEPLT(0:NBBTS))
*----
*  Read HISTORY information for cell specified
*  1) Read fuel type information
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      IFT=IDFUEL(IUPDB,IUPDC)
      WRITE(NAMP,'(A4,I8.8)') 'FUEL',IFT 
      CALL LCMLEN(IPHST,NAMP,ILCMLN,ILCMTY)
      IF(ILCMLN .EQ. 0) CALL XABORT(NAMSBR//
     >':/ Fuel type absent -- BURNUP creation impossible')
      CALL LCMSIX(IPHST,NAMP,ILCMUP)
      CALL LCMLEN(IPHST,'ISOTOPESUSED',ILCMLN,ILCMTY)
      IF(ILCMLN .GT. 0 .AND. ILCMLN .LT. 3*MAXI*4) THEN
        CALL LCMGET(IPHST,'ISOTOPESUSED',NAMIH(1,1)) 
        CALL LCMGET(IPHST,'ISOTOPESMIX',MIXIH(1))
      ELSE
        CALL XABORT(NAMSBR//
     >  ':/ Isotopes are absent -- BURNUP creation impossible')
      ENDIF
      CALL LCMSIX(IPHST,NAMP,ILCMDN)
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6001) 'FUEL TYPE',IFT
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6010) 
          WRITE(IOUT,6011) 
     >    (NAMIH(1,ISO),NAMIH(2,ISO),NAMIH(3,ISO),ISO=1,MAXI)
          WRITE(IOUT,6020) 
          WRITE(IOUT,6021) 
     >    (MIXIH(ISO),ISO=1,MAXI)
        ENDIF
      ENDIF 
*----
*  2) Real cell type information
*----
      ICT=IDCELL(IUPDB,IUPDC)
      WRITE(NAMP,'(A4,I8.8)') 'CELL',ICT
      CALL LCMLEN(IPHST,NAMP,ILCMLN,ILCMTY)
      IF(ILCMLN .EQ. 0) CALL XABORT(NAMSBR//
     >':/ Cell type absent -- BURNUP creation impossible')
      CALL LCMSIX(IPHST,NAMP,ILCMUP)
      IOK=-1
      CALL HSTGSD(IPHST ,MAXI  ,IOK   ,DENI  ,FDENC )
      IF(IOK .NE. 0) CALL XABORT(NAMSBR//
     >':/ Densities are absent -- BURNUP creation impossible')
      CALL LCMGET(IPHST,'DEPL-PARAM  ',BITH)
      CALL LCMSIX(IPHST,NAMP,ILCMDN)
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6001) 'CELL TYPE',ICT
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6100)  
          WRITE(IOUT,6110) (DENI(ISO),ISO=1,MAXI)
        ENDIF
      ENDIF
*----
*  Save isotopes names and mixtures on BURNUP
*----
      CALL LCMPUT(IPEVO,'ISOTOPESUSED',3*MAXI,3,NAMIH(1,1))
      CALL LCMPUT(IPEVO,'ISOTOPESMIX ',MAXI  ,1,MIXIH(1))
      CALL LCMPUT(IPEVO,'FUELDEN-INIT',2     ,2,FDENC)
*----
*  Save current burnup information as initial time step
*----
      FLXNOR=0.0
      DELTA(1)=0.0
      DELTA(2)=0.0
      ITS=0 
      DEPLT(ITS)=BITH(1)*CDAY
      CALL LCMPUT(IPEVO,'DEPL-TIMES  ',1     ,2,DEPLT(ITS))
      WRITE(NAMTIM,'(A8,I4.4)') 'DEPL-DAT',ITS+1
      CALL LCMSIX(IPEVO,NAMTIM,ILCMUP)
      CALL LCMPUT(IPEVO,'ISOTOPESDENS',MAXI,2,DENI(1))
      CALL LCMPUT(IPEVO,'FLUX-NORM   ',   1,2,FLXNOR)
      CALL LCMPUT(IPEVO,'DELTA       ',   2,2,DELTA)
      CALL LCMPUT(IPEVO,'BURNUP-IRRAD',   2,2,BITH(2))
      CALL LCMSIX(IPEVO,NAMTIM,ILCMDN)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DEPLT,MIXIH,NAMIH)
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
 6110 FORMAT(1P,10E14.7)
      END 
