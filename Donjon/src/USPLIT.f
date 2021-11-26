*DECK USPLIT
      SUBROUTINE USPLIT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Provide a link between the material index and reactor geometry;
* create a matex object.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* J. Koclas, D. Sekki, V. Descotes
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
* The USPLIT: module specification is: 
* GEOM MATEX := USPLIT: { GEOM | GEOMOLD } :: (desclink) ;
* where
*   GEOM : name of a \emph{geometry} object. This object is defined in creation 
*     (appears only on LHS) or modification (appears on both LHS and RHS) 
*     mode. An existing geometry previously created in the GEO: module is 
*     modified. Only 3-D Cartesian or 3-D Hexagonal reactor geometries are 
*     allowed.
*   MATEX name of a \emph{matex} object to be created by the module.
*   GEOMOLD : name of a \emph{geometry} object previously created in the GEO: 
*     module. This object must be specified in read-only mode (appears only on 
*     RHS). It is copied into GEOM at the beginning of USPLIT: module. Only 3-D 
*     Cartesian or 3-D Hexagonal reactor geometries are allowed.
*   (desclink) : structure describing the input data to the USPLIT: module.
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
      CHARACTER TEXT*12,HSIGN*12
      INTEGER ISTATE(NSTATE),NCODE(6),ICODE(6)
      REAL ZCODE(6)
      DOUBLE PRECISION DFLOT
      TYPE(C_PTR) IPGEO,IPMTX
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISPLX,ISPLY,ISPLZ,MAT,INDX,
     1 IRMIX,IFMIX,MIXA
      REAL, ALLOCATABLE, DIMENSION(:) :: XX,YY,ZZ
      LOGICAL LASBLY
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.LT.2)CALL XABORT('@USPLIT: 2 PARAMETERS EXPECTED.')
      IPGEO=KENTRY(1)
      IPMTX=KENTRY(2)
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))CALL XABORT('@USPLIT:'
     1 //' LCM OBJECT EXPECTED AT LHS.')
      IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2))CALL XABORT('@USPLIT:'
     1 //' LCM OBJECT EXPECTED AT LHS.')
      IF(JENTRY(1).EQ.0) THEN
        IF(NENTRY.LT.3)CALL XABORT('@USPLIT: 3 PARAMETERS EXPECTED.')
        IF((IENTRY(3).NE.1).AND.(IENTRY(3).NE.2))CALL XABORT('@USPLIT:'
     1  //' LCM OBJECT EXPECTED AT RHS.')
        IF(JENTRY(3).NE.2)CALL XABORT('@USPLIT: READ-ONLY MODE EXPECTE'
     1  //'D FOR L_GEOM.')
        CALL LCMGTC(KENTRY(3),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_GEOM')THEN
          TEXT=HENTRY(3)
          CALL XABORT('@USPLIT: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1    '. L_GEOM EXPECTED.')
        ENDIF
        CALL LCMEQU(KENTRY(3),IPGEO)
      ELSE IF(JENTRY(1).EQ.1) THEN
        CALL LCMGTC(IPGEO,'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_GEOM')THEN
          TEXT=HENTRY(1)
          CALL XABORT('@USPLIT: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1    '. L_GEOM EXPECTED.')
        ENDIF
      ELSE
        TEXT=HENTRY(1)
        CALL XABORT('@USPLIT: CREATE OR MODFICATION MODE EXPECTED FOR '
     1  //TEXT//'.')
      ENDIF
      IF(JENTRY(2).NE.0)CALL XABORT('@USPLIT: CREATE MODE EXPECTED FOR'
     1 //' L_MATEX.')
      LASBLY=.FALSE.
*----
*  RECOVER STATE-VECTOR INFORMATION
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPGEO,'STATE-VECTOR',ISTATE)
      IGEO=ISTATE(1)
      NMAT=ISTATE(7)
      CALL LCMLEN(IPGEO,'MIX-ASBLY',NITMA,ITYP)
      IF(NITMA.EQ.0) THEN
        NMIXA=0
      ELSE
        NMIXA=NITMA/2
      ENDIF
      CALL LCMLEN(IPGEO,'A-NMIXP',NITMA,ITYP)
      IF(NITMA.NE.1) THEN
        NMIXP=0
      ELSE
        CALL LCMGET(IPGEO,'A-NMIXP',NMIXP)
      ENDIF
      IF((IGEO.NE.7).AND.(IGEO.NE.9))CALL XABORT('@USPLIT: ONLY '
     1 //' 3D-CARTESIAN OR 3D-HEXAGONAL GEOMETRY ALLOWED.')
*----
*  READ INFORMATION
*----
      IMPX=1
      NREFL=0
      NFUEL=0
   10 CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@USPLIT: CHARACTER DATA EXPECTED.')
      IF(TEXT.EQ.'EDIT') THEN
*       READ PRINTING INDEX
        CALL REDGET(ITYP,IMPX,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@USPLIT: INTEGER FOR EDIT EXPECTED.')
      ELSE IF(TEXT.EQ.'NGRP') THEN
*       NUMBER OF ENERGY GROUPS
        CALL REDGET(ITYP,NGRP,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@USPLIT: INTEGER DATA EXPECTED(1).')
        IF(NGRP.LT.1)CALL XABORT('@USPLIT: INVALID NUMBER FOR NGRP.')
      ELSE IF(TEXT.EQ.'MAXR') THEN
*       MAXIMUM NUMBER OF REGIONS
        CALL REDGET(ITYP,MAXR,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@USPLIT: INTEGER DATA EXPECTED(2).')
        IF(MAXR.LT.1)CALL XABORT('@USPLIT: INVALID NUMBER FOR MAXR.')
      ELSE IF(TEXT.EQ.'NMIX') THEN
*       MAXIMUM NUMBER OF REGIONS
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@USPLIT: INTEGER DATA EXPECTED(3).')
        IF(NITMA.LT.NMAT)CALL XABORT('@USPLIT: INVALID NMIX < NMAT.')
        NMAT=NITMA
      ELSE IF(TEXT.EQ.'NREFL') THEN
*       NUMBER OF REFLECTOR TYPES
        CALL REDGET(ITYP,NREFL,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@USPLIT: INTEGER DATA EXPECTED(4).')
        IF(NREFL.LT.1)CALL XABORT('@USPLIT: REFLECTOR NOT DEFINED.')
        IF(NREFL.GT.NMAT-1)CALL XABORT('@USPLIT: WRONG NUMBER OF '
     1  //'REFLECTOR TYPES.')
*       REFLECTOR MIXTURES
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(TEXT.NE.'RMIX')CALL XABORT('@USPLIT: KEYWORD RMIX EXPECTED.')
        ALLOCATE(IRMIX(NREFL))
        DO I=1,NREFL
          CALL REDGET(ITYP,IRMIX(I),FLOT,TEXT,DFLOT)
          IF(ITYP.NE.1)CALL XABORT('@USPLIT: INTEGER DATA EXPECTED(5).')
          IF(IRMIX(I).LT.1)CALL XABORT('@USPLIT: INVALID RMIX'
     1    //' NUMBER < 1.')
          IF(IRMIX(I).GT.NMAT)CALL XABORT('@USPLIT: INVALID RMIX NUMBE'
     1    //'R > NBMIX.')
        ENDDO
      ELSE IF(TEXT.EQ.'NFUEL') THEN
*       NUMBER OF FUEL TYPES
        CALL REDGET(ITYP,NFUEL,FLOT,TEXT,DFLOT)
        IF(ITYP.EQ.1) THEN
*         general definition of fuel mixture
          IF(NFUEL.LT.1)CALL XABORT('@USPLIT: FUEL NOT DEFINED.')
          IF(NREFL+NFUEL.NE.NMAT)THEN
            WRITE(IOUT,*)'@USPLIT: NREFL:',NREFL,', NFUEL:',NFUEL 
            WRITE(IOUT,*)'@USPLIT: TOTAL NUMBER OF MATERIALS ',NMAT 
            CALL XABORT('@USPLIT: WRONG NUMBER OF REFLECTOR OR FUEL TY'
     1      //'PES.')
          ENDIF
*         FUEL MIXTURES
          CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
          IF(TEXT.NE.'FMIX')CALL XABORT('@USPLIT: KEYWORD FMIX EXPECTE'
     1    //'D.')
          ALLOCATE(IFMIX(NFUEL))
          DO I=1,NFUEL
            CALL REDGET(ITYP,IFMIX(I),FLOT,TEXT,DFLOT)
            IF(ITYP.NE.1)CALL XABORT('@USPLIT: INTEGER DATA EXPECTED(6'
     1      //').')
            IF(IFMIX(I).LT.1)CALL XABORT('@USPLIT: INVALID FMIX NUMBER'
     1      //' < 1.')
            IF(IFMIX(I).GT.NMAT)CALL XABORT('@USPLIT: INVALID FMIX NUM'
     1      //'BER > NBMIX.')
          ENDDO
        ELSEIF((ITYP.EQ.3).AND.(TEXT.EQ.'ASBLY')) THEN
*         automatic definition of fuel mixture from unfolded geometry
*         by NAP:
          LASBLY=.TRUE.
          IF(NMIXA.EQ.0)CALL XABORT('@USPLIT: No assembly previously '
     1    //'defined (NMIXA=0).')
          IF(NMIXP.EQ.0)CALL XABORT('@USPLIT: No assembly previously '
     1    //'defined (NMIXP=0).')
          NFUEL=NMIXA*(NMIXP+1)
          ALLOCATE(MIXA(2*NMIXA))
          CALL LCMGET(IPGEO,'MIX-ASBLY',MIXA)
          ALLOCATE(IFMIX(NFUEL))
          DO I=1,NMIXA
            IFMIX((I-1)*(NMIXP+1)+1)=MIXA(I)
            DO J=1,NMIXP
              IFMIX((I-1)*(NMIXP+1)+1+J)=MIXA(I+NMIXA)+J-1
            ENDDO
          ENDDO
        ELSE
          CALL XABORT('@USPLIT: INTEGER DATA or ASBLY keyword is EXPEC'
     1    //'TED.')
        ENDIF
      ELSE IF(TEXT.EQ.';') THEN
        GO TO 20
      ELSE
        CALL XABORT('@USPLIT: FINAL ; EXPECTED.')
      ENDIF
      GO TO 10
*----
*  RECOVER INFORMATION
*----
   20 IMPX0=MAX(0,IMPX-1)
      ALLOCATE(ISPLX(MAXR),ISPLY(MAXR),ISPLZ(MAXR),MAT(MAXR))
      ALLOCATE(XX(MAXR+1),YY(MAXR+1),ZZ(MAXR+1))
      CALL READ3D(MAXR,MAXR,MAXR,MAXR,IPGEO,IHEX,IR,ILK,SIDE,XX,YY,ZZ,
     1 IMPX0,LX,LY,LZ,MAT,NEL,NCODE,ICODE,ZCODE,ISPLX,ISPLY,ISPLZ,ISPLH,
     2 ISPLL)
      IF((NEL.NE.LX*LY*LZ).AND.(IHEX.EQ.0))CALL XABORT('@USPLIT: WRONG'
     1 // ' GEOMETRY.')
      IF((NEL.NE.LX*LZ).AND.(IHEX.NE.0))CALL XABORT('@USPLIT: WRONG' 
     1 // ' HEXAGONAL GEOMETRY, WRONG NUMBER OF ELEMENTS.')
      DEALLOCATE(ISPLZ,ISPLY,ISPLX)
      IF(IHEX.EQ.0) THEN
        CALL LCMPUT(IPMTX,'MESHX',LX+1,2,XX)
        CALL LCMPUT(IPMTX,'MESHY',LY+1,2,YY)
        CALL LCMPUT(IPMTX,'MESHZ',LZ+1,2,ZZ)
        CALL LCMPUT(IPGEO,'MESHX',LX+1,2,XX)
        CALL LCMPUT(IPGEO,'MESHY',LY+1,2,YY)
        CALL LCMPUT(IPGEO,'MESHZ',LZ+1,2,ZZ)
        DEALLOCATE(ZZ,YY,XX)
      ELSE
        CALL LCMPUT(IPMTX,'SIDE',1,2,SIDE)
        CALL LCMPUT(IPMTX,'MESHZ',LZ+1,2,ZZ)
        CALL LCMPUT(IPGEO,'SIDE',1,2,SIDE)
        CALL LCMPUT(IPGEO,'MESHZ',LZ+1,2,ZZ)
        DEALLOCATE(ZZ,YY,XX)
        LY=0
      ENDIF
*     NEW MATERIAL INDEX
      ALLOCATE(INDX(NEL))
      IF(NREFL.EQ.0) ALLOCATE(IRMIX(1))
      IF(NFUEL.EQ.0) ALLOCATE(IFMIX(1))
      CALL USPMIX(IPMTX,NEL,NREFL,NFUEL,MAT,IRMIX,IFMIX,INDX,NMIX)
*     MODIFY GEOMETRY
      ISTATE(3)=LX
      ISTATE(4)=LY
      ISTATE(5)=LZ
      ISTATE(6)=NEL
      ISTATE(7)=NMIX
      ISTATE(11)=0
      CALL LCMPUT(IPGEO,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPGEO,'MIX',NEL,1,INDX)
      DEALLOCATE(IFMIX,IRMIX,INDX,MAT)
      IF(LASBLY) DEALLOCATE(MIXA)
*----
*  STATE-VECTOR FOR MATEX
*----
      NTOT=NEL
      HSIGN='L_MATEX'
      CALL LCMPTC(IPMTX,'SIGNATURE',12,1,HSIGN)
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NGRP
      ISTATE(2)=NMIX
      ISTATE(3)=NREFL
      ISTATE(4)=NFUEL
      ISTATE(5)=NTOT
      ISTATE(6)=IGEO
      ISTATE(7)=NEL
      ISTATE(8)=LX
      ISTATE(9)=LY
      ISTATE(10)=LZ
      CALL LCMPUT(IPMTX,'STATE-VECTOR',NSTATE,1,ISTATE)
      IF(IMPX.GT.0) CALL LCMLIB(IPMTX)
      RETURN
      END
