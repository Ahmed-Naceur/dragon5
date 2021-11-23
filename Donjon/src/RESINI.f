*DECK RESINI
      SUBROUTINE RESINI(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Construct or modify a fuel-map object.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki and V. Descotes
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
* The RESINI: module specifications are:
* Option 1
* FLMAP MATEX := RESINI: MATEX [COMPO] :: (descresini1) ;
* Option 2
* FLMAP := RESINI: FLMAP [FLMAP2] :: (descresini2) ;
* where
*   FLMAP : name of the \emph{resini} object that will contain the fuel-lattice 
*     information. If FLMAP appears on both LHS and RHS, it will be updated; 
*     otherwise, it is created.
* 
*   MATEX : name of the \emph{matex} object specified in the modification mode. 
*     MATEX is required only when FLMAP is created.
*   COMPO : name of the \emph{multicompo} data structure (L\_COMPO signature) 
*     where the detailed subregion geometry at assembly level is stored.
*   FLMAP2 : name of the \emph{resini} object that contains the fuel-lattice 
*     information to recover from.
*   (descresini1) : structure describing the main input data to
*     the RESINI: module. Note that this input data is mandatory and
*     must be specified only when FLMAP is created.
*   (descresini2) : structure describing the input data for global and local 
*     parameters. This data is permitted to be modified in the subsequent calls 
*     to the RESINI: module.
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
      CHARACTER TEXT*12,HSIGN*12,HSIGN2*12
      INTEGER ISTATE(NSTATE),IGST(NSTATE)
      LOGICAL LNEW,LCPO,LMAP2
      TYPE(C_PTR) IPMTX,IPMAP,JPMAP,IPCPO,IPMP2
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.GT.3)CALL XABORT('@RESINI: 2 or 3 PARAMETERS ALLOWED.')
      LCPO=.FALSE.
      IPCPO=C_NULL_PTR
      IPMP2=C_NULL_PTR
      IF(IENTRY(1).GT.2) CALL XABORT('@RESINI: INVALID FIRST PARAMETER'
     1 //' TYPE.')
      LNEW=.TRUE.
      LMAP2=.FALSE.
      HSIGN2=' '
      IF(NENTRY.GE.2) CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN2)
      IF((NENTRY.EQ.1).OR.(HSIGN2.EQ.'L_MAP'))THEN
        IF(JENTRY(1).NE.1) CALL XABORT('@RESINI: OBJECT IN MODIFICATIO'
     1  //'N MODE EXPECTED.')
        CALL LCMGTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_MAP')THEN
          TEXT=HENTRY(1)
          CALL XABORT('@RESINI: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1    '. L_MAP EXPECTED.')
        ENDIF
        IF(JENTRY(1).NE.1)CALL XABORT('@RESINI: MODIFICATION MODE EX'
     1   //'PECTED FOR THE FUEL-MAP OBJECT.')
        LNEW=.FALSE.
        IF(HSIGN2.EQ.'L_MAP') THEN
          LMAP2=.TRUE.
          IPMP2=KENTRY(2)
        ENDIF
      ELSE
        IF(HSIGN2.NE.'L_MATEX')THEN
          TEXT=HENTRY(2)
          CALL XABORT('@RESINI: SIGNATURE OF '//TEXT//' IS '//HSIGN2//
     1    '. L_MATEX EXPECTED.')
        ENDIF
        IF(JENTRY(2).NE.1)CALL XABORT('@RESINI: MODIFICATION MODE EX'
     1   //'PECTED FOR THE MATEX OBJECT.')
        HSIGN='L_MAP'
        CALL LCMPTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
        IPMTX=KENTRY(2)
        IF(NENTRY.EQ.3) THEN
          LCPO=.TRUE.
          CALL LCMGTC(KENTRY(3),'SIGNATURE',12,1,HSIGN)
          IF(HSIGN.NE.'L_MULTICOMPO')THEN
            TEXT=HENTRY(3)
            CALL XABORT('@RESINI: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1      '. L_MULTICOMPO EXPECTED.')
          ENDIF
          IPCPO=KENTRY(3)
        ENDIF
      ENDIF
      IPMAP=KENTRY(1)
*----
*  RECOVER INFORMATION
*----
      IMPX=1
      CALL XDISET(ISTATE,NSTATE,0)
      IF(LNEW)THEN
        NPARM=0
        CALL LCMGET(IPMTX,'STATE-VECTOR',ISTATE)
        IGEO=ISTATE(6)
        IF((IGEO.NE.7).AND.(IGEO.NE.9))CALL XABORT('@RESINI: ONLY'
     1  //' 3D-CARTESIAN OR 3D-HEXAGONAL GEOMETRY ALLOWED.')
        NGRP=ISTATE(1)
        NFUEL=ISTATE(4)
        LX=ISTATE(8)
        LY=ISTATE(9)
        LZ=ISTATE(10)
*       MAIN INPUT
        CALL RESDRV(IPMAP,IPMTX,NFUEL,LX,LY,LZ,IMPX,IGEO,NCH,NB,NTOT,
     1  NCOMB,NSIMS,NASB,NAX,NAY,NIS,IPCPO)
        CALL XDISET(ISTATE,NSTATE,0)
        ISTATE(1)=NB
        ISTATE(2)=NCH
        ISTATE(3)=NCOMB
        ISTATE(4)=NGRP
        ISTATE(12)=IGEO
        ISTATE(7)=NFUEL
        ISTATE(8)=NPARM
        ISTATE(9)=NTOT
        ISTATE(13)=NSIMS
        ISTATE(14)=NASB
        ISTATE(15)=NAX
        ISTATE(16)=NAY
        ISTATE(18)=NIS
      ELSE
        CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
        NB=ISTATE(1)
        NCH=ISTATE(2)
        NCOMB=ISTATE(3)
        IGEO=ISTATE(12)
        NFUEL=ISTATE(7)
        NPARM=ISTATE(8)
        NTOT=ISTATE(9)
        NSIMS=ISTATE(13)
        NASB=ISTATE(14)
        NAX=ISTATE(15)
        NAY=ISTATE(16)
        NIS=ISTATE(18)
      ENDIF
      CALL XDISET(IGST,NSTATE,0)
      JPMAP=LCMGID(IPMAP,'GEOMAP')
      CALL LCMGET(JPMAP,'STATE-VECTOR',IGST)
      NX=IGST(3)
      NY=IGST(4)
      NZ=IGST(5)
      IF(IGEO.EQ.9) NY=1
*     INPUT OF PARAMETERS
      CALL RESPAR(IPMAP,NCH,NB,NFUEL,NCOMB,NPARM,NX,NY,NZ,NSTATE,
     1 ISTATE,IMPX,NASB,LMAP2,IPMP2)
      CALL LCMPUT(IPMAP,'STATE-VECTOR',NSTATE,1,ISTATE)
      IF(IMPX.GT.0)WRITE(IOUT,100) IMPX,(ISTATE(I),I=1,9),ISTATE(12),
     1 ISTATE(13),ISTATE(18)
      IF(IMPX.GT.5)CALL LCMLIB(IPMAP)
      RETURN
*
  100 FORMAT(/8H OPTIONS/8H -------/
     1 7H IMPX  ,I6,30H   (0=NO PRINT/1=SHORT/2=MORE)/
     2 7H NB    ,I6,39H   (NUMBER OF FUEL BUNDLES PER CHANNEL)/
     3 7H NCH   ,I6,28H   (NUMBER OF FUEL CHANNELS)/
     4 7H NCOMB ,I6,31H   (NUMBER OF COMBUSTION ZONES)/
     5 7H NGRP  ,I6,28H   (NUMBER OF ENERGY GROUPS)/
     6 7H INTER ,I6,26H   (TYPE OF INTERPOLATION)/
     7 7H ISHIFT,I6,28H   (NUMBER OF BUNDLE SHIFTS)/
     8 7H NFUEL ,I6,25H   (NUMBER OF FUEL TYPES)/
     9 7H NPARM ,I6,25H   (NUMBER OF PARAMETERS)/
     1 7H NTOT  ,I6,33H   (TOTAL NUMBER OF FUEL BUNDLES)/
     2 7H IGEO  ,I6,28H   (7=CARTESIAN/9=HEXAGONAL)/
     3 7H NSIMS ,I6,35H   (ASSEMBLY LAYOUT IN SIM: MODULE)/
     4 7H NIS   ,I6,38H   (NUMBER OF PARTICULARIZED ISOTOPES))
      END
