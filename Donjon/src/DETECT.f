*DECK DETECT
      SUBROUTINE DETECT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
*  This module compute detectors readings
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): 
* E. Varin, M. Guyot
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
* The DETECT: module specifications are:
* DETEC := DETECT: DETEC FLUX  TRACK GEOM :: (descdetect) ;
* where
*   DETEC : name of the \emph{detect} containing the detector positions and 
*     responses. 
*   FLUX  : name of the \emph{flux} containing the flux solution computed by
*     the FLUD: or FLPOW: modules. To obtain a correct result, the best is to
*     use a normalized flux, coming from the FLPOW: module. In this case, the 
*     fluxes are normalized to the reactor power.
*   TRACK : name of the \emph{track} containing the TRIVAC tracking.
*   GEOM  : name of the \emph{geometry} containing the mesh-splitting 
*     geometry created by the USPLIT: or GEO: modules.
*   (descdetect) : structure containing the data to module DETECT:.
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
      INTEGER NSTATE,IOUT
      PARAMETER (NSTATE=40,IOUT=6)
      REAL    FLOT,DT,VNORM
      TYPE(C_PTR) IPFLU,JPFLUX,IPTRK,IPGEO,IPDET
      INTEGER ISTATE(NSTATE),NEL,NUN,
     1        PARAM(NSTATE),I,IPRT,ITYP,NITMA,KC,NX,NY,NZ,NXP1,
     2        NYP1,NZP1,NGRP,IGR,GEOTYP,ILONG,ITYLCM,IUN
      LOGICAL LTRK,LFLU,LGEO,LDET,LHEX,LNORM,LSIMEX,LPARAB
      CHARACTER HSIGN*12,TEXT*12
      DOUBLE PRECISION DFLOT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEYF
      REAL, ALLOCATABLE, DIMENSION(:) :: MESHX,MESHY,MESHZ,FUNKN
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FLU
*----
*  PARAMETERS VALIDATION
*----  
      IF(NENTRY.LE.3) CALL XABORT('@DETECT: FOUR PARAMETER EXPECTED.')
      LTRK = .FALSE.
      LFLU = .FALSE.
      LGEO = .FALSE.
      LDET = .FALSE.
      IPFLU = C_NULL_PTR
      IPTRK = C_NULL_PTR
      IPGEO = C_NULL_PTR
      IPDET = C_NULL_PTR
      DO 10 I=1,NENTRY
        IF((IENTRY(I).EQ.1).OR.(IENTRY(I).EQ.2)) THEN
          TEXT=HENTRY(I)
          CALL LCMSIX(KENTRY(I),' ',0)
          CALL LCMGTC(KENTRY(I),'SIGNATURE',12,1,HSIGN)
          IF (HSIGN.EQ.'L_DETECT') THEN
            IPDET=KENTRY(I)
            LDET = .TRUE.
            IF(JENTRY(I).NE.1) CALL XABORT('@DET'
     +      //'ECT: MODIFICATION MODE EXPECTED FOR OBJECT'//HSIGN//'.')
          ELSEIF (HSIGN.EQ.'L_GEOM') THEN
            IPGEO=KENTRY(I)
            LGEO = .TRUE.
            IF(JENTRY(I).NE.2) CALL XABORT('@DET'
     +      //'ECT: READ-ONLY MODE EXPECTED FOR OBJECT'//HSIGN//'.')
         ELSEIF (HSIGN.EQ.'L_TRACK') THEN
           IF (.NOT.LTRK) THEN
             IPTRK=KENTRY(I)
             LTRK = .TRUE.
           IF(JENTRY(I).NE.2) CALL XABORT('@DET'
     +      //'ECT: READ-ONLY MODE EXPECTED FOR OBJECT'//HSIGN//'.')
           ELSE
             CALL XABORT('@DETECT: ONLY ONE L_TRACK FILE IS REQUIRED')
           ENDIF
         ELSEIF ((HSIGN.EQ.'L_FLUX').AND.(.NOT.LFLU)) THEN
           IPFLU=KENTRY(I)
           LFLU = .TRUE.
           IF(JENTRY(I).NE.2) CALL XABORT('@DET'
     +      //'ECT: READ-ONLY MODE EXPECTED FOR OBJECT'//HSIGN//'.')
           ELSE
             CALL XABORT('@DETECT: ONLY ONE L_FLUX FILE IS REQUIRED')
           ENDIF
         ELSE
           CALL XABORT('@DETECT: INVALIV OBJECT='//TEXT)
         ENDIF
  10  CONTINUE
      IF (.NOT.(LFLU.AND.LGEO.AND.LTRK.AND.LDET))
     +  CALL XABORT('@DETECT: MISSING OBJECTS IN CALL')
*----
*  READ DATA
*----  
      IPRT = 1
      LHEX = .FALSE.
      LNORM = .FALSE.
      LSIMEX = .FALSE.
      LPARAB = .TRUE.
      DT = 0.0
      KC = 0

 15   CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.EQ.3) THEN
        IF (TEXT.EQ.'EDIT') THEN
          CALL REDGET(ITYP,IPRT,FLOT,TEXT,DFLOT)
          IF (ITYP.NE.1)
     +      CALL XABORT('@DETECT: INTEGER DATA EXPECTED(1)')
        ELSEIF (TEXT.EQ.'TIME') THEN
          CALL REDGET(ITYP,NITMA,DT,TEXT,DFLOT)
          IF (ITYP.NE.2)
     +      CALL XABORT('@DETECT: REAL DATA EXPECTED(1)')
        ELSEIF (TEXT.EQ.'REF') THEN
          CALL REDGET(ITYP,KC,FLOT,TEXT,DFLOT)
          IF (ITYP.NE.1)
     +      CALL XABORT('@DETECT: INTEGER DATA EXPECTED(2)')
        ELSEIF (TEXT.EQ.'SIMEX') THEN
            LSIMEX = .TRUE.
        ELSEIF (TEXT.EQ.'SPLINE') THEN
            IF(.NOT.LSIMEX) CALL XABORT('@DETECT: WRONG KEYWORD, '
     +               //' SIMEX REQUIRED')
            LPARAB = .FALSE.
        ELSEIF (TEXT.EQ.'PARAB') THEN
            IF(.NOT.LSIMEX) CALL XABORT('@DETECT: WRONG KEYWORD, '
     +               //' SIMEX REQUIRED')
            LPARAB = .TRUE.
        ELSEIF (TEXT.EQ.'NORM') THEN
          LNORM = .TRUE.
          CALL REDGET(ITYP,NITMA,VNORM,TEXT,DFLOT)
          IF (ITYP.NE.2)
     +      CALL XABORT('@DETECT: REAL DATA EXPECTED(3)')
          IF( VNORM.EQ.0.0 )CALL XABORT('@DETECT: ILLEGAL VALUE '
     +                 // 'OF NORM')
        ELSEIF (TEXT.EQ.';') THEN
          GOTO 20
        ELSE
          CALL XABORT('@DETECT: CONTROLLED TYPE EXPECTED'//TEXT)
        ENDIF
      ELSE
        CALL XABORT('@DETECT: CHARACTER DATA EXPECTED(1)')
      ENDIF
      GOTO 15
*----
*  RECOVER L_GEOM INFORMATION
*----  
 20   IF(DT.EQ.0.0) CALL XABORT('@DETECT: TIME NOT SET')
      IF(LSIMEX.AND.LNORM) CALL XABORT('@DETECT: WRONG ASSOCIATION '
     +  //' SIMEX INT AND NORMALIZATION')
      CALL LCMGET(IPDET,'STATE-VECTOR',PARAM)
      CALL LCMGET(IPGEO,'STATE-VECTOR',ISTATE)
      GEOTYP = ISTATE(1)
      IF(PARAM(3).EQ.1) LHEX = .TRUE.
      IF(LSIMEX.AND.GEOTYP.NE.7)
     +  CALL XABORT('@DETECT: SIMEX INTERPOLATION ONLY FOR 3D '
     +        //'CARTESIAN')
      IF((LHEX.AND.(GEOTYP.LT.8)).OR.(.NOT.LHEX.AND.(GEOTYP.GE.8)))
     +   CALL XABORT('@DETECT: INCOMPATIBLE DETECT WITH GEOMETRY')
      IF(GEOTYP.LT.5.OR.GEOTYP.EQ.6)
     +  CALL XABORT('@DETECT: GEOMETRY TYPE NOT SUPPORTED IN DETECT')
      NX = ISTATE(3)
      NY = ISTATE(4)
      IF(NY.EQ.0) NY=1
      NZ = ISTATE(5)
      IF(NZ.EQ.0) NZ=1
      NXP1 = NX+1
      NYP1 = NY+1
      NZP1 = NZ+1
      ALLOCATE(MESHX(NXP1),MESHY(NYP1),MESHZ(NZP1))
      IF((GEOTYP.EQ.7).OR.(GEOTYP.EQ.5)) THEN
        CALL LCMGET(IPGEO,'MESHX',MESHX)
        CALL LCMGET(IPGEO,'MESHY',MESHY)
      ELSE
        MESHY(1)=0.
        MESHY(2)=1.
        MESHX(1)=0.
        MESHX(2)=1.
      ENDIF
      IF(GEOTYP.EQ.9.OR.GEOTYP.EQ.7)THEN
       CALL LCMGET(IPGEO,'MESHZ',MESHZ)
      ELSE IF(GEOTYP.EQ.5.OR.GEOTYP.EQ.8)THEN
        MESHZ(1)=0.
        MESHZ(2)=1.
      ENDIF
*----
*  RECOVER L_TRACK INFORMATION
*----  
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NEL = ISTATE(1)
      NUN = ISTATE(2)
      ALLOCATE(KEYF(NEL))
      CALL LCMGET(IPTRK,'KEYFLX',KEYF)
      CALL LCMGET(IPDET,'STATE-VECTOR',ISTATE)
      NGRP = ISTATE(1)
*----
*  RECOVER L_FLUX INFORMATION
*----  
      CALL LCMGET(IPFLU,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.NGRP)CALL XABORT('@DETECT: NUMBER OF ENERGY '
     +  //'GROUPS INCOMPATIBLE BETWEEN FLUX AND DETECT')
       ALLOCATE(FLU(NUN,NGRP))
       CALL LCMSIX(IPFLU,' ',0)
       JPFLUX=LCMGID(IPFLU,'FLUX')
       CALL LCMLEL(JPFLUX,1,ILONG,ITYLCM)
       ALLOCATE(FUNKN(ILONG))
       DO 30 IGR=1,NGRP
        CALL LCMGDL(JPFLUX,IGR,FUNKN)
        DO 25 IUN=1,NUN
          FLU(IUN,IGR)=FUNKN(IUN)
   25   CONTINUE
   30  CONTINUE
       DEALLOCATE(FUNKN)
*----
*  CALL DRIVER
*----  
      CALL DETCDRV(IPDET,NGRP,NEL,NUN,NX,NY,NZ,MESHX,MESHY,MESHZ,KEYF,
     + FLU,IPRT,KC,DT,LHEX,LSIMEX,LNORM,VNORM,LPARAB)
*----
*  RELEASE MEMORY
*----  
      DEALLOCATE(FLU,KEYF,MESHX,MESHY,MESHZ)
      RETURN
      END
