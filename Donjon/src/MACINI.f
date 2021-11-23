*DECK MACINI
      SUBROUTINE MACINI(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Construct a new macrolib that will contain one mixture number per
* material region; fuel-map macrolib is required.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* J. Koclas, E. Varin, D. Sekki
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
* The MACINI: module specification is:
* MACRO2 MATEX := MACINI: MATEX MACRO [ MACFL  ] :: [ EDIT iprint ] [ FUEL ] ;
* where
*   MACRO2 : name of the extended \emph{macrolib} to be created by the module.
*   MATEX  : name of the \emph{matex} object containing an extended material 
*     index over the reactor geometry. MATEX must be specified in the 
*     modification mode; it will store the recovered h-factors per each fuel 
*     region.
*   MACRO  : name of a \emph{macrolib}, created using either MAC:, CRE:, NCR: 
*     or AFM: module, for the evolution-independent material properties 
*     (see structure (desccre1) or refer to the DRAGON user guide).
*   MACFL  : name of a fuel-map \emph{macrolib}, created using either CRE:, 
*     NCR: or AFM: module, for the interpolated fuel properties (see structure
*     (desccre2) or refer to the DRAGON user guide).
*   EDIT   : keyword used to set iprint.
*   iprint : integer index used to control the printing on screen: = 0 for
*     no print; = 1 for minimum printing; larger values produce increasing 
*     amounts of output. The default value is iprint = 1.
*   FUEL   : keyword used to indicate that MACRO is a fuel-map \emph{macrolib}
*     in case where only two RHS objects are defined. By default, MACRO contains
*     evolution-independent cross sections.
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
      INTEGER ISTATE(NSTATE)
      DOUBLE PRECISION DFLOT
      LOGICAL LMAP,LWD1,LWD2
      TYPE(C_PTR) IPMAC,IPMTX,IPMAC1,IPMAC2,JPMAC,KPMAC
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIX
      REAL, ALLOCATABLE, DIMENSION(:) :: HFAC,WDLA
*----
*  PARAMETER VALIDATION
*----
      IF((NENTRY.LE.2).OR.(NENTRY.GE.5)) 
     1  CALL XABORT('@MACINI: 3 OR 4 PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))CALL XABORT('@MACINI:'
     1 //' LCM OBJECT EXPECTED AT LHS.')
      IF(JENTRY(1).NE.0)CALL XABORT('@MACINI: CREATE MODE EXPECTED'
     1 //' FOR L_MACROLIB AT LHS.')
      IPMAC=KENTRY(1)
*     L_MATEX
      IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2))CALL XABORT('@MACINI:'
     1 //' LCM OBJECT EXPECTED AT RHS.')
      IF(JENTRY(2).NE.1)CALL XABORT('@MACINI: MODIFICATION MODE EX'
     1 //'PECTED FOR L_MATEX OBJECT.')
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_MATEX')THEN
        TEXT=HENTRY(2)
        CALL XABORT('@MACINI: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1  '. L_MATEX EXPECTED AT RHS.')
      ENDIF
      IPMTX=KENTRY(2)
      DO IEN=3,NENTRY
      IF((IENTRY(IEN).NE.1).AND.(IENTRY(IEN).NE.2))CALL XABORT('@MACIN'
     1 //'I: LCM OBJECT EXPECTED AT RHS.')
      IF(JENTRY(IEN).NE.2)CALL XABORT('@MACINI: READ-ONLY MODE EXPEC'
     1 //'TED FOR THE LCM OBJECTS AT RHS.')
      ENDDO
*     L_MACROLIB(1)
      CALL LCMGTC(KENTRY(3),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_MACROLIB')THEN
        TEXT=HENTRY(3)
        CALL XABORT('@MACINI: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1  '. FIRST MACROLIB EXPECTED AT RHS.')
      ENDIF
      IPMAC1=KENTRY(3)
*     L_MACROLIB(2)
      IF(NENTRY.EQ.4) THEN
        CALL LCMGTC(KENTRY(4),'SIGNATURE',12,1,HSIGN)
      	IF(HSIGN.NE.'L_MACROLIB')THEN
          TEXT=HENTRY(4)
          CALL XABORT('@MACINI: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1    '. FUEL-MAP MACROLIB EXPECTED AT RHS.')
      	ENDIF
      	IPMAC2=KENTRY(4)
      ELSE
      	IPMAC2=C_NULL_PTR
      ENDIF
*----
*  RECOVER INFORMATION
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPMAC1,'STATE-VECTOR',ISTATE)
*     MACROLIB(1)-INFO
      NGRP=ISTATE(1)
      NMIX1=ISTATE(2)
      NL=ISTATE(3)
      NF1=ISTATE(4)
      NDEL1=ISTATE(7)
      NDEL2=0
      LEAK=ISTATE(9)
      NW1=ISTATE(10)
*     MACROLIB(2)-INFO
      NF2=1
      IF(NENTRY.EQ.4) THEN
      	CALL XDISET(ISTATE,NSTATE,0)
      	CALL LCMGET(IPMAC2,'STATE-VECTOR',ISTATE)
      	NMIX2=ISTATE(2)
        NDEL2=ISTATE(7)
        NL=MAX(ISTATE(3),NL)
        NW2=ISTATE(10)
        NF2=ISTATE(4)
        IF((NF2.NE.NF1).AND.(NF1.GT.1)) THEN
          WRITE(IOUT,*)'MACROLIB=',HENTRY(1),' NF=',NF1,' (0 EXPECTED)'
          WRITE(IOUT,*)'MACROLIB=',HENTRY(2),' NF=',NF2
          CALL XABORT('@MACINI: INCONSISTENT NUMBER OF FISSILE ISOTOPE'
     1    //'S.')
        ENDIF
        IF(ISTATE(1).NE.NGRP)CALL XABORT('@MACINI: DIFFERENT NGRP'
     1 //' NUMBER IN TWO MACROLIB OBJECTS.')
      	IF(ISTATE(3).NE.NL)CALL XABORT('@MACINI: INCONSISTENT NL '
     1 //'NUMBER IN TWO MACROLIB OBJECTS.')
      	IF(ISTATE(9).NE.LEAK)CALL XABORT('@MACINI: DIFFERENT LEAK'
     1 //' NUMBER IN TWO MACROLIB OBJECTS.')
      ENDIF
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPMTX,'STATE-VECTOR',ISTATE)
*     MATEX-INFO
      NMIX=ISTATE(2)
      NTOT=ISTATE(5)
      ALLOCATE(MIX(NTOT))
      CALL XDISET(MIX,NTOT,0)
      CALL LCMGET(IPMTX,'MAT',MIX)
      IMPX=1
      LMAP=.FALSE.
   10 CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.EQ.10)GOTO 20
      IF(ITYP.NE.3)CALL XABORT('@MACINI: CHARACTER DATA EXPECTED.')
      IF(TEXT.EQ.'EDIT')THEN
*       READ PRINTING INDEX
        CALL REDGET(ITYP,IMPX,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@MACINI: INTEGER FOR EDIT EXPECTED.')
      ELSE IF(TEXT.EQ.'FUEL')THEN
*       ASSUME FUEL-MAP MACROLIB
        IF(NENTRY.NE.3) CALL XABORT ('@MACINI: 3 PARAMETERS EXPECTED.')
        LMAP=.TRUE.
      ELSE IF(TEXT.EQ.';')THEN
        GO TO 20
      ELSE
        CALL XABORT('@MACINI: FINAL ; EXPECTED.')
      ENDIF
      GO TO 10
*----
*  NEW MACROLIB CREATION
*----
   20 IF(IMPX.GT.1)WRITE(IOUT,*)'NUMBER OF ENERGY GROUPS  ',NGRP
      IF(IMPX.GT.1)WRITE(IOUT,*)'TOTAL NUMBER OF MIXTURES ',NMIX
*     DO NOT INCLUDE FUEL PROPERTIES
      IF(IMPX.GT.0)WRITE(IOUT,*)'** TREATING FIRST MACROLIB **'
      CALL MACCRE(IPMAC1,IPMAC,NL,NW1,NF1,NGRP,NMIX1,NMIX,NTOT,MIX,LMAP,
     1 IMPX)
      IF(IMPX.GT.1)CALL LCMLIB(IPMAC)
*     INCLUDE FUEL PROPERTIES
      IF(NENTRY.EQ.4) THEN
      	LMAP=.TRUE.
      	IF(IMPX.GT.0)WRITE(IOUT,*)'** TREATING FUEL-MAP MACROLIB **'
      	CALL MACCRE(IPMAC2,IPMAC,NL,NW2,NF2,NGRP,NMIX2,NMIX,NTOT,MIX,
     1  LMAP,IMPX)
      ENDIF
      DEALLOCATE(MIX)
*----
*  RECOVER LAMBDA-D
*----
      CALL LCMLEN(IPMAC1,'LAMBDA-D',LENGTH,ITYLCM)
      LWD1=(LENGTH.EQ.NDEL1).AND.(NDEL1.GT.0)
      LWD2=.FALSE.
      IF(NENTRY.EQ.4) THEN
        CALL LCMLEN(IPMAC2,'LAMBDA-D',LENGTH,ITYLCM)
        LWD2=(LENGTH.EQ.NDEL2).AND.(NDEL2.GT.0)
      ENDIF
      NDEL=0
      IF(LWD1) THEN
        NDEL=NDEL1
        ALLOCATE(WDLA(NDEL))
        CALL LCMGET(IPMAC1,'LAMBDA-D',WDLA)
        CALL LCMPUT(IPMAC,'LAMBDA-D',NDEL,2,WDLA)
        DEALLOCATE(WDLA)
      ELSE IF(LWD2) THEN
        NDEL=NDEL2
        ALLOCATE(WDLA(NDEL))
        CALL LCMGET(IPMAC2,'LAMBDA-D',WDLA)
        CALL LCMPUT(IPMAC,'LAMBDA-D',NDEL,2,WDLA)
        DEALLOCATE(WDLA)
      ENDIF
*----
*  STATE-VECTOR
*----
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NGRP
      ISTATE(2)=NMIX
      ISTATE(3)=NL
      IF(NENTRY.EQ.3) THEN
        ISTATE(4)=NF1
      ELSE IF(NENTRY.EQ.4) THEN
        ISTATE(4)=NF2
      ENDIF
      ISTATE(7)=NDEL
      ISTATE(9)=LEAK
      CALL LCMPUT(IPMAC,'STATE-VECTOR',NSTATE,1,ISTATE)
      HSIGN='L_MACROLIB'
      CALL LCMPTC(IPMAC,'SIGNATURE',12,1,HSIGN)
      IF(IMPX.GT.0)CALL LCMLIB(IPMAC)
*----
*  RECOVER H-FACTOR AND SAVE ON L_MATEX
*----
      ALLOCATE(HFAC(NMIX*NGRP))
      JPMAC=LCMGID(IPMAC,'GROUP')
      DO JGR=1,NGRP
      KPMAC=LCMGIL(JPMAC,JGR)
      CALL LCMLEN(KPMAC,'H-FACTOR',LENGT,ITYP)
      IF(LENGT.NE.NMIX)CALL XABORT('@MACINI: UNABLE TO FIND H'
     1 //'-FACTOR BLOCK DATA IN THE NEW MACROLIB.')
      CALL LCMGET(KPMAC,'H-FACTOR',HFAC((JGR-1)*NMIX+1))
      ENDDO
      CALL LCMPUT(IPMTX,'H-FACTOR',NMIX*NGRP,2,HFAC)
      DEALLOCATE(HFAC)
      RETURN
      END
