*DECK XENLIB
      SUBROUTINE XENLIB(IPLIB,MAXMIX,NMIX,NBISO,NGRP,XEN)

*
*-----------------------------------------------------------------------
*
*Purpose:
* Update the macroscopic cross sections thanks to the Xenon distribution
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
*
*Author(s): 
* M. Guyot
*
*Parameters: input/output
* IPLIB   adress of the L_LIBRARY
* MAXMIX  maximum number of mixtures in the library
* NMIX    number of mixtures present in the library
* NBISO   number of isotopes
* NGRP    number of energy groups
* XEN     xenon concentrations in each bundle
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER MAXMIX,NMIX,NBISO,NGRP
      REAL XEN(NMIX)
*----
*  LOCAL VARIABLES
*----
      INTEGER IMIX,ISO
      REAL TMPDAY(3)
      CHARACTER TEXT*8
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIX
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NAME,USED
      REAL, ALLOCATABLE, DIMENSION(:) :: DENS
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASK,MASKL
*----
*  SCRATCH STORAGE ALLOCATION
*   MIX     'ISOTOPESMIX'
*   NAME    'ISOTOPESNAME'
*   USED    'ISOTOPESUSED'
*   DENS    'ISOTOPESDENS' updated
*----
      ALLOCATE(MIX(NBISO),NAME(3,NBISO),USED(3,NBISO),DENS(NBISO))
*----
*  RECOVER INFORMATION
*----
      CALL LCMGET(IPLIB,'ISOTOPESMIX',MIX)
      CALL LCMGET(IPLIB,'ISOTOPERNAME',NAME)
      CALL LCMGET(IPLIB,'ISOTOPESUSED',USED)
      CALL LCMGET(IPLIB,'ISOTOPESDENS',DENS)
*----
*  PERFORM CALCULATION
*----
      IMIX=0
      DO 10 ISO=1,NBISO
        WRITE(TEXT,'(2A4)') (NAME(I,ISO),I=1,2)
        IF(TEXT.EQ.'Xe135   ') THEN
          IMIX=IMIX+1
          DENS(ISO)=XEN(IMIX)
        ENDIF
   10 CONTINUE
       
      IF(IMIX.NE.NMIX) CALL XABORT('@XENLIB: Xe135 SHOULD BE EXTRACTED '
     1 //'IN ALL MIXTURES .')

      CALL LCMPUT(IPLIB,'ISOTOPESDENS',NBISO,2,DENS)
*----
*  UPDATE MACROSCOPIC XS
*----
      ALLOCATE(MASK(MAXMIX),MASKL(NGRP))
      CALL XDLSET(MASK,MAXMIX,.FALSE.)
      CALL XDLSET(MASKL,NGRP,.TRUE.)
      DO 20 I=1,NBISO
       IBM=MIX(I)
       MASK(IBM)=.TRUE.
   20 CONTINUE
      ITSTMP=0
      TMPDAY(1)=0.0
      TMPDAY(2)=0.0
      TMPDAY(3)=0.0
*----
*  CALL THE DRAGON SUBROUTINE FOR THE COMPUTATION OF THE MACROSCOPIC XS
*----
      CALL LIBMIX(IPLIB,MAXMIX,NGRP,NBISO,USED,MIX,DENS,MASK,MASKL,
     1 ITSTMP,TMPDAY)
      DEALLOCATE(MASKL,MASK)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DENS,USED,NAME,MIX)
      RETURN
      END
