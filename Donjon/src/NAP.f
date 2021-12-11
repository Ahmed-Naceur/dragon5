*DECK NAP
      SUBROUTINE NAP(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* - Construct an 'enriched' multicompo with additional information  
* needed by Pin Power Reconstruction.
* - Performed the Pin Power Reconstruction
* - Split geometry from homogeneous to heterogeneous assemblies
*     Note : this function is also called directly from the RESINI: module
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal.
*
*Author(s): 
* R. Chambon
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
* The NAP: calling specifications are:
* Option 1: 
* COMPO := NAP: COMPO TRKNAM FLUNAM :: (descnap1) ;
* Option 2: 
* MAP := NAP: MAP TRKNAM FLUNAM MATEX MACRES :: (descnap2) ;
* Option 3: 
* GEONEW := NAP: GEOOLD COMPO :: (descnap3) ;
* where
*   COMPO : name of the \emph{multicompo} data structure (L\_COMPO signature) 
*     where the detailed subregion properties will be stored.
*   TRKNAM : name of the read-only \emph{tracking} data structure 
*     (L\_TRACK signature) containing the tracking. 
*   FLUNAM : name of the read-only \emph{fluxunk} data structure 
*     (L\_FLUX signature) containing a transport solution.
*   MAP : name of the \emph{map} data structure (L\_MAP signature) containing 
*     fuel regions description, global and local parameter information (burnup, 
*     fuel/coolant temperatures, coolant density, etc). A previous call to the 
*     FLPOW: module is highly recommended prior to the pin-power reconstruction 
*     to normalize the flux and compute the assembly power. If not, the 
*     pin-power reconstruction  will be normalized using the whole core power 
*     instead of a normalization for each assembly.
*   MATEX : name of the read-only \emph{matex} data structure 
*     (L\_MATEX signature). The object corresponds to the heterogeneously 
*     splited geometry.
*   MACRES : name of the read-only \emph{macrolib} data structure 
*     (L\_MACROLIB signature) containing a cross section for the fuel. The 
*     \emph{macrolib} data structure must have been created with a 
*     \emph{multicompo} data structure with pin level properties (transport 
*     flux, H-factor, infinite domain diffusion flux).
*   GEONEW : name of the created \emph{geometry} data structure 
*     (L\_GEOM signature) containing the detailed core geometry definition at 
*     heterogeneous assembly level.
*   GEOOLD : name of the read-only \emph{geometry} data structure 
*     (L\_GEOM signature) containing the core geometry definition with 
*     homogeneous assembly (only 1 mesh per assembly mandatory).
*   (descnap1) : structure containing the input data to this module to compute 
*     additional properties for subregions
*   (descnap2) : structure containing the input data to this module to perform 
*     pin power reconstruction
*   (descnap3) : structure containing the input data to this module to 
*     automatically define the core geometry with heterogeneous assembly
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      INTEGER IOUT,MAXPAR,MAXLIN,MAXVAL,NSTATE,MAXADD
      REAL REPS
      PARAMETER (REPS=1.0E-4,IOUT=6,MAXPAR=50,MAXLIN=50,MAXVAL=200,
     1 NSTATE=40,MAXADD=10)
      TYPE(C_PTR) IPCPO,IPFLU,IPTRK,IPMAP,IPMTX,IPGNW,IPGOD,IPMPP,IPMAC
      CHARACTER TEXT*12,HSIGN*12
      INTEGER KCHAR(3)
      INTEGER IEN,I
      LOGICAL ldebug

      IPMAP=C_NULL_PTR
      IPMTX=C_NULL_PTR
      IPGNW=C_NULL_PTR
      IPGOD=C_NULL_PTR
      IPCPO=C_NULL_PTR
      IPMPP=C_NULL_PTR
      IPMAC=C_NULL_PTR
      
      ldebug=.false.
      if(ldebug)write(6,*) 'NAP begin debug'
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.LE.2)CALL XABORT('@NAP: AT LEAST 3 PARAMETERS'
     > //' EXPECTED.')
      
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))CALL XABORT('@NAP'
     1 //': LCM OBJECT EXPECTED AT LHS.')
* NAPGEO
      if(ldebug)write(6,*) 'NAP begin NAPGEO'
      IF(JENTRY(1).EQ.0) THEN 
        IPGNW=KENTRY(1)
        DO IEN=2,3
          IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))CALL XABORT('@NAP'
     1     //': LCM OBJECT EXPECTED AT RHS.')
          CALL LCMGET(KENTRY(IEN),'SIGNATURE',KCHAR)
          WRITE(HSIGN,'(3A4)') (KCHAR(I),I=1,3)
          IF(HSIGN.EQ.'L_GEOM')THEN
            IPGOD=KENTRY(IEN)
          ELSEIF(HSIGN.EQ.'L_MULTICOMPO')THEN
            IPCPO=KENTRY(IEN)
          ELSE
              CALL XABORT('@NAP: COMPO OR GEOM OBJECT EXPECTED.')
          ENDIF
        ENDDO
        GOTO 3000
      ENDIF
* NAPCPO + NAPPPR
      if(ldebug)write(6,*) 'NAP begin NAPCPO + NAPPPR'
      CALL LCMGET(KENTRY(1),'SIGNATURE',KCHAR)
      WRITE(HSIGN,'(3A4)') (KCHAR(I),I=1,3)
      IF(HSIGN.EQ.'L_MULTICOMPO')THEN
        IPCPO=KENTRY(1)
      ELSEIF(HSIGN.EQ.'L_MAP')THEN
        IPMAP=KENTRY(1)
      ELSE
        CALL XABORT('@NAP: L_MULTICOMPO or L_MAP EXPECTED.')
      ENDIF
      DO 5 IEN=2,3
        IF((IENTRY(IEN).NE.1).AND.(IENTRY(IEN).NE.2))CALL XABORT('@N'
     1 //'AP: LCM OBJECT EXPECTED AT RHS.')
        IF(JENTRY(IEN).NE.2)CALL XABORT('@NAP: LCM OBJECT IN READ-ON'
     1 //'LY MODE EXPECTED AT RHS.')
        CALL LCMGET(KENTRY(IEN),'SIGNATURE',KCHAR)
        WRITE(HSIGN,'(3A4)') (KCHAR(I),I=1,3)
        IF(HSIGN.EQ.'L_FLUX')THEN
          IPFLU=KENTRY(IEN)
        ELSEIF(HSIGN.EQ.'L_TRACK')THEN
          IPTRK=KENTRY(IEN)
        ELSE
            CALL XABORT('@NAP: FLUX OR TRACKING OBJECT EXPECTED.')
        ENDIF
    5 CONTINUE
      IF(NENTRY.EQ.3) GOTO 1000
* NAPPPR
      if(ldebug)write(6,*) 'NAP begin NAPPPR'
      DO 7 IEN=4,NENTRY
        IF((IENTRY(IEN).NE.1).AND.(IENTRY(IEN).NE.2))CALL XABORT('@N'
     1   //'AP: LCM OBJECT EXPECTED AT RHS.')
        IF(JENTRY(IEN).NE.2)CALL XABORT('@NAP: LCM OBJECT IN READ-ON'
     1   //'LY MODE EXPECTED AT RHS.')
        CALL LCMGET(KENTRY(IEN),'SIGNATURE',KCHAR)
        WRITE(HSIGN,'(3A4)') (KCHAR(I),I=1,3)
C        IF(HSIGN.EQ.'L_MAP')THEN
C          IPMPP=KENTRY(IEN)
C        ELSEIF((HSIGN.EQ.'L_MATEX'))THEN
        IF((HSIGN.EQ.'L_MATEX'))THEN
          IPMTX=KENTRY(IEN)
        ELSEIF((HSIGN.EQ.'L_MACROLIB'))THEN
          IPMAC=KENTRY(IEN)
        ELSE
          TEXT=HENTRY(IEN)
          CALL XABORT('@NAP: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1    '. L_MATEX or L_MACROLIB EXPECTED.')
        ENDIF
    7 CONTINUE
      GOTO 2000
*----
*  enriched L_MULTICOMPO  computation
*----
 1000 CALL NAPCPO(IPCPO,IPTRK,IPFLU,NSTATE)
      GOTO 9000
*----
*  Pin Power Reconstruction
*----
 2000 CALL NAPPPR(IPMAP,IPTRK,IPFLU,IPMTX,IPMAC,NSTATE)
      GOTO 9000
*----
*  Automatic geometry unfolding
*----
 3000 CALL NAPGEO(IPGNW,IPGOD,IPCPO,NSTATE)
      GOTO 9000
*----
* END
*----
*
 9000 RETURN
      END
