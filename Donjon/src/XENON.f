*DECK XENON
      SUBROUTINE XENON(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Computing the Xenon distribution
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): 
* M. Guyot
*
*Parameters: input
* NENTRY  number of data structures transfered to this module.
* HENTRY  name of the data structures.
* IENTRY  data structure type where:
*         IENTRY=1 for LCM memory object;
*         IENTRY=2 for XSM file;
*         IENTRY=3 for sequential binary file;
*         IENTRY=4 for sequential ASCII file.
* JENTRY  access permission for the data structure where:
*         JENTRY=0 for a data structure in creation mode;
*         JENTRY=1 for a data structure in modifications mode;
*         JENTRY=2 for a data structure in read-only mode.
* KENTRY  data structure pointer.
*
*Comments:
* The XENON: module specification is:
* MICROLIB := XENON: MICROLIB [ POWER ] :: (descxenon) ;
* where
*   MICROLIB : name of a \emph{library} object, that will be updated by the 
*     XENON : module. The Xenon should be extracted in this library for the use 
*     of this module.
*   POWER : name of a \emph{power} object containing the bundle fluxes, 
*     previously computed by the FLPOW: module. The fluxes should be normalized 
*     to the reactor power.
*   (descxenon) : structure describing the input data to the XENON: module.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40,IOUT=6)
      CHARACTER HSIGN*12,TEXT*12
      INTEGER ISTATE(NSTATE),ITYP,NITMA
      REAL FLOT
      DOUBLE PRECISION DFLOT
      LOGICAL LINI
      TYPE(C_PTR) IPLIB,IPPOW
      REAL, ALLOCATABLE, DIMENSION(:) :: XEN
*----
*  PARAMETER VALIDATION
*----
      IPLIB=C_NULL_PTR
      IPPOW=C_NULL_PTR
      IF((NENTRY.NE.1).AND.(NENTRY.NE.2)) 
     1  CALL XABORT('@XENON: 1 OR 2 PARAMETERS EXPECTED.')
      DO I=1,NENTRY
         IF((IENTRY(I).NE.1).AND.(IENTRY(I).NE.2))
     1   CALL XABORT('@XENON: LCM OBJECT EXPECTED AT LHS')
      ENDDO
      IF(JENTRY(1).NE.1)CALL XABORT('@XENON: MODIFICATION MODE EXPECTED'
     1 //' FOR L_LIBRARY.')
      IF(NENTRY.EQ.2) THEN
        IF(JENTRY(2).NE.2)CALL XABORT('@XENON: READ-ONLY MODE EXPECTED'
     1   //' FOR L_POWER AT LHS.')
      ENDIF
      DO IEN=1,NENTRY
        CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
*       L_LIBRARY
        IF(HSIGN.EQ.'L_LIBRARY')THEN
          IPLIB=KENTRY(IEN)
*       L_POWER
        ELSEIF(HSIGN.EQ.'L_POWER')THEN
          IPPOW=KENTRY(IEN)
        ELSE
          TEXT=HENTRY(IEN)
          CALL XABORT('@XENON: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1    '. L_LIBRARY OR L_POWER EXPECTED.')
        ENDIF
      ENDDO
*----
*  RECOVER INFORMATION
*----
*     L_LIBRARY
      CALL LCMSIX(IPLIB,' ',0)
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      MAXMIX=ISTATE(1)
      NBISO=ISTATE(2)
      NGRP=ISTATE(3)
      NMIX=ISTATE(14)
*     L_POWER
      IF(C_ASSOCIATED(IPPOW)) THEN
        CALL LCMSIX(IPPOW,' ',0)
        CALL XDISET(ISTATE,NSTATE,0)
        CALL LCMGET(IPPOW,'STATE-VECTOR',ISTATE)
        IF(ISTATE(1).NE.NGRP)CALL XABORT('@XENON: DIFFERENT NGR'
     1   //'P NUMBER IN L_LIBRARY AND L_POWER OBJECT.')
        NCH=ISTATE(6)
        NB=ISTATE(7)
        IF(NCH*NB.NE.NMIX)CALL XABORT('@XENON: DIFFERENT '
     1   //'MIXTURE NUMBER IN L_LIBRARY AND L_POWER OBJECT.')
      ENDIF
*----
*  READ INPUT DATA
*----
      IPRT=0
      LINI=.FALSE.
*     READ KEYWORD
   10 CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@XENON: CHARACTER DATA EXPECTED(1).')
      IF(TEXT.EQ.'EDIT')THEN
*       PRINTING INDEX
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@XENON: INTEGER DATA EXPECTED.')
        IPRT=MAX(0,NITMA)
        GOTO 10
      ELSEIF(TEXT.EQ.'INIT')THEN
        LINI=.TRUE.
        GOTO 10
      ELSEIF(TEXT.EQ.';')THEN
        GOTO 20
      ELSE
*       KEYWORD DOES NOT MATCH
        CALL XABORT('@XENON: WRONG KEYWORD:'//TEXT//'.')
      ENDIF

   20 IF((.NOT.C_ASSOCIATED(IPPOW)).AND.(.NOT.LINI)) THEN
         CALL XABORT('@XENON: L_POWER OBJECT REQUIRED .')
      ENDIF
      ALLOCATE(XEN(NMIX))
*----
*  COMPUTE THE VALUE OF THE XENON CONCENTRATIONS
*----
      IF(.NOT.LINI) THEN
        CALL XENCAL(IPLIB,IPPOW,NB,NCH,NGRP,NMIX,NBISO,XEN)
      ELSE
        CALL XDRSET(XEN,NMIX,0.0)
      ENDIF
*----
*  PUT THE CONCENTRATIONS IN THE LIBRARY AND COMPUTE NEW XS
*----
      CALL XENLIB(IPLIB,MAXMIX,NMIX,NBISO,NGRP,XEN)
      DEALLOCATE(XEN)
      RETURN
      END
