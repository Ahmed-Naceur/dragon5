*DECK SOUR
      SUBROUTINE SOUR(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Define and set fixed external sources.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): C. Bienvenue
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): creation or modification type(L_SOURCE);
*         HENTRY(2): read-only type(L_MACROLIB);
*         HENTRY(3): read-only type(L_TRACKING);
*         HENTRY(4): read-only type(L_GEOM);
*         HENTRY(5+): optional read-only type(L_SOURCE).
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
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
      PARAMETER (IOUT=6,NSTATE=40,MAXPAR=10)
      TYPE(C_PTR) IPSOUR,JPSOUR,KPSOUR,IPMAC,IPTRK,IPGEOM,IPSOUR2,
     1 JPSOUR2,KPSOUR2,KPSOUR_INFO,KPSOUR2_INFO
      CHARACTER HSIGN*12,TEXT12*12
      INTEGER ISTATE(NSTATE),MESH_LEN(3),GN,NBS1,NBS2
      DOUBLE PRECISION DFLOTT
      REAL DIR(3),NORM,BSNORM
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISISOUR
      REAL, ALLOCATABLE, DIMENSION(:) :: XXX,YYY,ZZZ,ISOUR,
     1 XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,SOURTEMP
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SUNKNO,BS1,BS2,BS
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: BSINFO1,BSINFO2,BSINFO
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.LT.4) CALL XABORT('SOUR:AT LEAST FOUR PARAMETER'
     1 //' EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('SOUR: LI'
     1 //'NKED LIST OR XSM FILE EXPECTED AT LHS.')
      IF(JENTRY(1).NE.0.AND.JENTRY(1).NE.1) CALL XABORT('SOUR: ENTRY IN'
     1 //' CREATE OR MODIFICATION MODE EXPECTED.')
      IPSOUR=KENTRY(1)
      DO I=2,NENTRY
        IF((IENTRY(I).NE.1).AND.(IENTRY(I).NE.2)) CALL XABORT('SOUR: '
     1  //'LINKED LIST OR XSM FILE EXPECTED AT RHS.')
        IF(JENTRY(I).NE.2) CALL XABORT('SOUR: ENTRY IN READ-ONLY MODE'
     1  //' EXPECTED.')
      ENDDO
*----
*  RECOVER MACROLIB INFORMATION
*----
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.EQ.'L_MACROLIB') THEN
        IPMAC=KENTRY(2)
      ELSE IF(HSIGN.EQ.'L_LIBRARY') THEN
        IPMAC=LCMGID(KENTRY(2),'MACROLIB')
      ELSE
        TEXT12=HENTRY(2)
        CALL XABORT('SOUR: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1  '. L_MACROLIB OR L_LIBRARY EXPECTED')
      ENDIF
      CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
      NG=ISTATE(1)
      NMAT=ISTATE(2)
      NANIS=ISTATE(3)-1
      IADJ=ISTATE(13)

*----
*  RECOVER TRACKING INFORMATION
*----
      IPTRK=KENTRY(3)
      CALL LCMGTC(IPTRK,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_TRACK') THEN
        TEXT12=HENTRY(3)
        CALL XABORT('SOUR: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1  '. L_TRACK EXPECTED.')
      ENDIF
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NREG=ISTATE(1)
      NUNS=ISTATE(2)
      NDIM=ISTATE(9)
      LX=ISTATE(12)
      LY=ISTATE(13)
      LZ=ISTATE(14)
      NLF=ISTATE(15)
      IF(ISTATE(4).NE.NMAT) CALL XABORT('SOUR: INVALID NUMBER OF MI'
     1 //'XTURES.')

*----
*  RECOVER GEOMETRY INFORMATION
*----
      IPGEOM=KENTRY(4)
      CALL LCMGTC(IPGEOM,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_GEOM') THEN
        TEXT12=HENTRY(4)
        CALL XABORT('SOUR: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1  '. L_GEOM EXPECTED.')
      ENDIF

      IF(NDIM.EQ.1) THEN
      CALL LCMLEN(IPGEOM,'MESHX',MESH_LEN(1),ITYLCM)
      ALLOCATE(XXX(MESH_LEN(1)))
      CALL LCMGET(IPGEOM,'MESHX',XXX)
      ELSE IF(NDIM.EQ.2) THEN
      CALL LCMLEN(IPGEOM,'MESHX',MESH_LEN(1),ITYLCM)
      CALL LCMLEN(IPGEOM,'MESHY',MESH_LEN(2),ITYLCM)
      ALLOCATE(XXX(MESH_LEN(1)),YYY(MESH_LEN(2)))
      CALL LCMGET(IPGEOM,'MESHX',XXX)
      CALL LCMGET(IPGEOM,'MESHY',YYY)
      ELSE IF(NDIM.EQ.3) THEN
      CALL LCMLEN(IPGEOM,'MESHX',MESH_LEN(1),ITYLCM)
      CALL LCMLEN(IPGEOM,'MESHY',MESH_LEN(2),ITYLCM)
      CALL LCMLEN(IPGEOM,'MESHZ',MESH_LEN(3),ITYLCM)
      ALLOCATE(XXX(MESH_LEN(1)),YYY(MESH_LEN(2)),ZZZ(MESH_LEN(3)))
      CALL LCMGET(IPGEOM,'MESHX',XXX)
      CALL LCMGET(IPGEOM,'MESHY',YYY)
      CALL LCMGET(IPGEOM,'MESHZ',ZZZ)
      ENDIF

*----
* INITIALIZATION FLUX AND RECOVER INPUT SOURCE INFORMATION IF AVAILABLE
*----

      ALLOCATE(SUNKNO(NUNS,NG),SOURTEMP(NUNS))
      SUNKNO(:NUNS,:NG)=0.0
      NBS1=0
      BSNORM=-1.0

      IF(NDIM.EQ.1) THEN
        MAXL=1
      ELSE IF(NDIM.EQ.2) THEN
        MAXL=MAX(LX,LY)
      ELSE
        MAXL=MAX(LX*LY,LX*LZ,LY*LZ)
      ENDIF

      IF(JENTRY(1).EQ.1) THEN
      IPSOUR2=KENTRY(1)
      CALL LCMGTC(IPSOUR2,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_SOURCE') THEN
        TEXT12=HENTRY(1)
        CALL XABORT('SOUR: SIGNATURE OF'//TEXT12//' IS '//HSIGN//
     1  '. L_SOURCE EXPECTED.')
      ENDIF
      ! VOLUMIC SOURCES
      IF(IADJ.EQ.0) THEN
        JPSOUR2=LCMGID(IPSOUR2,'DSOUR') 
      ELSE IF(IADJ.EQ.1) THEN
        JPSOUR2=LCMGID(IPSOUR2,'ASOUR') 
      ENDIF
      KPSOUR2=LCMGIL(JPSOUR2,1)
      DO IG=1,NG
        CALL LCMGDL(KPSOUR2,IG,SOURTEMP(:NUNS))
        SUNKNO(:NUNS,IG)=SUNKNO(:NUNS,IG)+SOURTEMP(:NUNS)  
      ENDDO
      ! BOUNDARY SOURCES
      CALL LCMLEN(IPSOUR2,'NBS',ISBS,ITYLCM)
      IF(ISBS.GT.0) THEN
      CALL LCMGET(IPSOUR2,'NBS',NBS1)
      ALLOCATE(BS1(MAXL,NBS1),BSINFO1(3,NBS1))
      JPSOUR2=LCMGID(IPSOUR2,'DBSOUR')
      KPSOUR2_INFO=LCMGIL(JPSOUR2,1)
      KPSOUR2=LCMGIL(JPSOUR2,2)
      DO N=1,NBS1
        CALL LCMGDL(KPSOUR2_INFO,N,BSINFO1(1,N))
        CALL LCMGDL(KPSOUR2,N,BS1(1,N))
      ENDDO
      CALL LCMLEN(IPSOUR2,'NORM-FS',ILEN,ITYLCM)
      IF(ILEN.NE.0) CALL LCMGET(IPSOUR2,'NORM-FS',BSNORM)
      ENDIF
      ENDIF

      IF(.NOT.ALLOCATED(BS1)) ALLOCATE(BS1(MAXL,0))
      IF(.NOT.ALLOCATED(BSINFO1)) ALLOCATE(BSINFO1(3,0))
      
      IF(NENTRY.GE.5) THEN
      DO 5 NSOURIN=5,(NENTRY)
      IPSOUR2=KENTRY(NSOURIN)
      CALL LCMGTC(IPSOUR2,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_SOURCE') THEN
        TEXT12=HENTRY(NSOURIN)
        CALL XABORT('SOUR: SIGNATURE OF'//TEXT12//' IS '//HSIGN//
     1  '. L_SOURCE EXPECTED.')
      ENDIF
      ! VOLUMIC SOURCES
      IF(IADJ.EQ.0) THEN
        JPSOUR2=LCMGID(IPSOUR2,'DSOUR') 
      ELSE IF(IADJ.EQ.1) THEN
        JPSOUR2=LCMGID(IPSOUR2,'ASOUR') 
      ENDIF
      KPSOUR2=LCMGIL(JPSOUR2,1)
      DO IG=1,NG
        CALL LCMGDL(KPSOUR2,IG,SOURTEMP(:NUNS))
        SUNKNO(:NUNS,IG)=SUNKNO(:NUNS,IG)+SOURTEMP(:NUNS)
      ENDDO
      ! BOUNDARY SOURCES
      CALL LCMLEN(IPSOUR2,'NBS',ILEN,ITYLCM)
      IF(ILEN.GT.0) CALL XABORT('MULTIPLE BOUNDARY SOURCE AS INPUT NOT'
     1 //'IMPLEMENTED.')
    5 CONTINUE
      ENDIF

*----
*  READ THE INPUT DATA
*----

! PARAMETERS VALUES
      STYPE=-1
      NSOUR=-1
      ALLOCATE(ISOUR(NG),ISISOUR(NG))
      ISISOUR=0    
      ISOUR=0.0
      ISIDEF=0
      ISXDEF=0
      ISYDEF=0
      ISZDEF=0
      MONOP=0
      ISDIRDEF=0

  10  CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
  11  IF(INDIC.NE.3) CALL XABORT('SOUR: CHARACTER DATA EXPECTED.')
! ISOTROPIC SOURCE KEYWORD
      IF(TEXT12.EQ.'ISO') THEN
        IF(STYPE.NE.-1) CALL XABORT('SOUR: ONLY ONE SOURCE TYPE ALLOWED'
     1  //'.')
        STYPE=0
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('SOUR: INTEGER DATA EXPECTED'
     1  //' (NUMBER OF SOURCES).')
        NSOUR=NITMA
        IF(NSOUR.LE.0) CALL XABORT('SOUR: INVALID NUMBER OF SOURCES.')
        ALLOCATE(XMIN(NSOUR),XMAX(NSOUR),YMIN(NSOUR),YMAX(NSOUR),
     1  ZMIN(NSOUR),ZMAX(NSOUR))
! MONODIRECTIONNAL BOUNDARY SOURCE KEYWORD
      ELSE IF(TEXT12.EQ.'MONO') THEN
        IF(STYPE.NE.-1) CALL XABORT('SOUR: ONLY ONE SOURCE TYPE ALLOWED'
     1  //'.')
        STYPE=1
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('SOUR: INTEGER DATA EXPECTED'
     1  //' (NUMBER OF SOURCES).')
        NSOUR=NITMA
        IF(NSOUR.LE.0) CALL XABORT('SOUR: INVALID NUMBER OF SOURCES.')
        ALLOCATE(XMIN(NSOUR),XMAX(NSOUR),YMIN(NSOUR),YMAX(NSOUR),
     1  ZMIN(NSOUR),ZMAX(NSOUR))
! SOURCE INTENSITY PER ENERGY GROUP KEYWORD
      ELSE IF(TEXT12.EQ.'INTG') THEN
        ISIDEF=1
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.EQ.2) THEN
          ISOUR(1)=FLOTT
          IF(ISOUR(1).LT.0.0) CALL XABORT('SOUR: INVALID INTENSITY.')
          DO I=2,NG
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
            ISOUR(I)=FLOTT
            IF(ISOUR(I).LT.0.0) CALL XABORT('SOUR: INVALID INTENSITY.')
          ENDDO
          ISISOUR(1:NG)=1
        ELSE IF(INDIC.EQ.1) THEN
          GN=NITMA
          IF(GN.GT.NG.OR.GN.LT.1) CALL XABORT('SOUR: INVALID GROUP'
     1    //' NUMBER.')
          ISISOUR(GN)=1
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
          ISOUR(GN)=FLOTT      
          IF(ISOUR(I).LT.0.0) CALL XABORT('SOUR: INVALID INTENSITY.')
  12      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          IF(INDIC.EQ.1) THEN
            GN=NITMA
            IF(GN.GT.NG.OR.GN.LT.1) CALL XABORT('SOUR: INVALID GROUP'
     1      //' NUMBER.')
            CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
            ISOUR(GN)=FLOTT      
            IF(ISOUR(I).LT.0.0) CALL XABORT('SOUR: INVALID INTENSITY.')
            GO TO 12
          ELSE
            GO TO 11
          ENDIF
        ELSE
          CALL XABORT('SOUR: REAL OR INTEGER DATA EXPECTED.')
        ENDIF
! (X,Y,Z) LIMITS KEYWORDS
      ELSE IF(TEXT12.EQ.'XLIM') THEN
        DO I=1,NSOUR
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        XMIN(I)=FLOTT
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        XMAX(I)=FLOTT
        ISXDEF=1
        ENDDO
      ELSE IF(TEXT12.EQ.'YLIM') THEN       
        IF(NDIM.LT.2) CALL XABORT('SOUR: INVALID USE OF YLIM, DIM < 2 ')
        DO I=1,NSOUR
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        YMIN(I)=FLOTT
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        YMAX(I)=FLOTT
        ISYDEF=1
        ENDDO
      ELSE IF(TEXT12.EQ.'ZLIM') THEN
        IF(NDIM.LT.3) CALL XABORT('SOUR: INVALID USE OF YLIM, DIM < 3 ')
        DO I=1,NSOUR
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        ZMIN(I)=FLOTT
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        ZMAX(I)=FLOTT
        ISZDEF=1
        ENDDO
! BOUNDARY SOURCE LOCATION KEYWORD
      ELSE IF(TEXT12.EQ.'X-') THEN
        IF(MONOP.NE.0) CALL XABORT('SOUR: BOUNDARY SOURCE ALREADY'
     1  //'DEFINED')
        MONOP=-1
      ELSE IF(TEXT12.EQ.'X+') THEN
        IF(MONOP.NE.0) CALL XABORT('SOUR: BOUNDARY SOURCE ALREADY'
     1  //'DEFINED')
        MONOP=1
      ELSE IF(TEXT12.EQ.'Y-') THEN
        IF(NDIM.LT.2) CALL XABORT('SOUR: INVALID USE OF YLIM, DIM < 2 ')
        IF(MONOP.NE.0) CALL XABORT('SOUR: BOUNDARY SOURCE ALREADY'
     1  //'DEFINED')
        MONOP=-2
      ELSE IF(TEXT12.EQ.'Y+') THEN
        IF(NDIM.LT.2) CALL XABORT('SOUR: INVALID USE OF YLIM, DIM < 2 ')
        IF(MONOP.NE.0) CALL XABORT('SOUR: BOUNDARY SOURCE ALREADY'
     1  //'DEFINED')
        MONOP=2
      ELSE IF(TEXT12.EQ.'Z-') THEN
        IF(NDIM.LT.3) CALL XABORT('SOUR: INVALID USE OF YLIM, DIM < 3 ')
        IF(MONOP.NE.0) CALL XABORT('SOUR: BOUNDARY SOURCE ALREADY'
     1  //'DEFINED')
        MONOP=-3
      ELSE IF(TEXT12.EQ.'Z+') THEN
        IF(NDIM.LT.3) CALL XABORT('SOUR: INVALID USE OF YLIM, DIM < 3 ')
        IF(MONOP.NE.0) CALL XABORT('SOUR: BOUNDARY SOURCE ALREADY'
     1  //'DEFINED')
        MONOP=3
! MONODIRECTIONAL SOURCE DIRCTION (2 ANGLES)
      ELSE IF(TEXT12.EQ.'DIR') THEN
        ISDIRDEF=1 
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        TEMP1=FLOTT 
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        TEMP2=FLOTT 
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.EQ.2) THEN
          DIR(1)=TEMP1 !MU
          DIR(2)=TEMP2 !ETA
          DIR(3)=FLOTT !XI
          IF(DIR(1).GT.1.0.OR.DIR(1).LT.-1.0.OR.DIR(2).GT.1.0.OR.
     1    DIR(2).LT.-1.0.OR.DIR(3).GT.1.0.OR.DIR(3).LT.-1.0.OR.
     2    DIR(1)**2+DIR(2)**2+DIR(3)**2.GT.1.0) CALL XABORT('SOUR:'
     3    //' INVALID DIRECTION COSINES VALUES.')
        ELSE
          IF(TEMP1.GT.180.0.OR.TEMP1.LT.0.0.OR.TEMP2.GT.90.0.OR.
     1    TEMP2.LT.0.0) CALL XABORT('SOUR: INVALID POLAR OR AZIMUTAL'
     2    //' ANGLE.')
          DIR(1)=COS(TEMP1*3.1416/180)                   !MU
          DIR(2)=SQRT(1-DIR(1)**2)*COS(TEMP2*3.1416/180) !ETA
          DIR(3)=SQRT(1-DIR(1)**2)*SIN(TEMP2*3.1416/180) !XI
          GO TO 11
        ENDIF
! END OF INPUT DATA
      ELSE IF(TEXT12.EQ.';') THEN
        GO TO 20
! ERROR IF KEYWORD IS UNRECONIZED
      ELSE
        CALL XABORT('SOUR: '//TEXT12//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 10
  20  CONTINUE


! OTHER VERIFICATIONS

      IF(STYPE.EQ.-1) CALL XABORT('SOUR: SOURCE TYPE NOT DEFINED.')
      IF(ISIDEF.EQ.0) CALL XABORT('SOUR: SOURCE INTENSITY NOT DEFINED.')
      IF(MONOP.EQ.0.AND.STYPE.EQ.1) CALL XABORT('SOUR: BOUNDARY SOURCE'
     1 //' POSITION NOT DEFINED.')
      IF(ISDIRDEF.EQ.0.AND.STYPE.EQ.1) CALL XABORT('SOUR:' 
     1 //' MONODIRECTIONAL SOURCE DIRECTION NOT DEFINED.')

      ! X-BOUNDARY
      IF(ISXDEF.EQ.1) THEN
      IF(MONOP.EQ.-1.OR.MONOP.EQ.1) CALL XABORT('SOUR: BOUNDARY SOURCE'
     1 //' X- OR X+ : XLIM SHOULD NOT BE DEFINED.')
      DO I=1,NSOUR
        IF(XMIN(I).LT.XXX(1).OR.XMAX(I).GT.XXX(MESH_LEN(1)))
     1  CALL XABORT('SOUR: SOURCE X-LIMITS OUTSIDE THE GEOMETRY.')
        IF(XMIN(I).GE.XMAX(I)) CALL XABORT('SOUR: SOURCE X-LIMITS'
     1  //' INVALID.')
      ENDDO
      ELSE
        IF(STYPE.EQ.0) CALL XABORT('SOUR: SOURCE X-LIMITS NOT DEFINED.')
      ENDIF

      ! Y-BOUNDARY
      IF(ISYDEF.EQ.1) THEN
      IF(MONOP.EQ.-2.OR.MONOP.EQ.2) CALL XABORT('SOUR: BOUNDARY SOURCE'
     1 //' Y- OR Y+ : YLIM SHOULD NOT BE DEFINED.')
      DO I=1,NSOUR
        IF(YMIN(I).LT.YYY(1).OR.YMAX(I).GT.YYY(MESH_LEN(2)))
     1  CALL XABORT('SOUR: SOURCE Y-LIMITS OUTSIDE THE GEOMETRY.')
        IF(YMIN(I).GE.YMAX(I)) CALL XABORT('SOUR: SOURCE Y-LIMITS'
     1  //' INVALID.')
      ENDDO
      ELSE
        IF(NDIM.GE.2.AND.STYPE.EQ.0) CALL XABORT('SOUR: SOURCE Y-LIMITS'
     1  //' NOT DEFINED.')
        IF(NDIM.GE.2.AND.STYPE.EQ.1.AND.(MONOP.NE.2.OR.MONOP.NE.-2))
     1  CALL XABORT('SOUR: SOURCE Y-LIMITS NOT DEFINED.')
      ENDIF

      ! Z-BOUNDARY
      IF(ISZDEF.EQ.1) THEN
      IF(MONOP.EQ.-3.OR.MONOP.EQ.3) CALL XABORT('SOUR: BOUNDARY SOURCE'
     1 //' Z- OR Z+ : ZLIM SHOULD NOT BE DEFINED.')
      DO I=1,NSOUR
        IF(ZMIN(I).LT.ZZZ(1).OR.ZMAX(I).GT.ZZZ(MESH_LEN(3)))
     1  CALL XABORT('SOUR: SOURCE Z-LIMITS OUTSIDE THE GEOMETRY.')
        IF(ZMIN(I).GE.ZMAX(I)) CALL XABORT('SOUR: SOURCE Z-LIMITS'
     1  //' INVALID.')
      ENDDO
      ELSE
        IF(NDIM.GE.3.AND.STYPE.EQ.0) CALL XABORT('SOUR: SOURCE Z-LIMITS'
     1  //' NOT DEFINED.')
        IF(NDIM.GE.3.AND.STYPE.EQ.1.AND.(MONOP.NE.3.OR.MONOP.NE.-3))
     1  CALL XABORT('SOUR: SOURCE Z-LIMITS NOT DEFINED.')
      ENDIF

*----
*  COMPUTE THE FIXED EXTERNAL SOURCES
*----

      NBS2=SUM(ISISOUR)*NSOUR
      ALLOCATE(BSINFO2(3,NBS2),BS2(MAXL,NBS2))

      IF(STYPE.EQ.0) THEN
         CALL SOURISO(IPTRK,IPGEOM,NREG,LX,LY,LZ,NG,NUNS,NDIM,
     1   NSOUR,ISOUR,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,XXX,YYY,ZZZ,
     2   MESH_LEN,SUNKNO,NORM)
      ELSE IF(STYPE.EQ.1) THEN
         CALL SOURMONO(IPTRK,IPGEOM,NREG,LX,LY,LZ,NG,NUNS,NDIM,
     1   NSOUR,ISOUR,ISISOUR,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,XXX,YYY,ZZZ,
     2   MESH_LEN,BSINFO2,BS2,MAXL,NORM,DIR,MONOP,NBS2)
      ELSE
        CALL XABORT('SOUR: SOURCE TYPE EXPECTED.')
      ENDIF

*----
*  SAVE THE FIXED EXTERNAL SOURCES AND NORM ON LCM
*----

      ! SAVE VOLUMETRIC SOURCE
      IOF=1
      NDIR=0
      NCST=0
      IF(IADJ.EQ.0) THEN
        NDIR=1
        JPSOUR=LCMLID(IPSOUR,'DSOUR',NDIR)
      ELSE IF(IADJ.EQ.1) THEN
        NCST=1
        JPSOUR=LCMLID(IPSOUR,'ASOUR',NCST)
      ENDIF
      KPSOUR=LCMLIL(JPSOUR,IOF,NG)
      DO IG=1,NG
        CALL LCMPDL(KPSOUR,IG,NUNS,2,SUNKNO(1,IG))
      ENDDO
      DEALLOCATE(SUNKNO)

      ! SAVE BOUNDARY SOURCE
      IF(STYPE.EQ.0.AND.ISBS.EQ.1) THEN
      JPSOUR=LCMLID(IPSOUR,'DBSOUR',2)
      KPSOUR_INFO=LCMLIL(JPSOUR,1,NBS1)
      KPSOUR=LCMLIL(JPSOUR,2,NBS1)
      DO N=1,NBS1
        CALL LCMPDL(KPSOUR_INFO,N,3,1,BSINFO1(1,N))
        CALL LCMPDL(KPSOUR,N,MAXL,2,BS1(1,N))
      ENDDO
      CALL LCMPUT(IPSOUR,'NBS',1,1,NBS1)      
      DEALLOCATE(BSINFO1,BS1)
      ELSE IF(STYPE.EQ.1.AND.ISBS.EQ.0) THEN
      JPSOUR=LCMLID(IPSOUR,'DBSOUR',2)
      KPSOUR_INFO=LCMLIL(JPSOUR,1,NBS2)
      KPSOUR=LCMLIL(JPSOUR,2,NBS2)
      DO N=1,NBS2
        CALL LCMPDL(KPSOUR_INFO,N,3,1,BSINFO2(1,N))
        CALL LCMPDL(KPSOUR,N,MAXL,2,BS2(1,N))
      ENDDO
      CALL LCMPUT(IPSOUR,'NBS',1,1,NBS2)
      CALL LCMLEN(IPSOUR,'NBS',ILEN,ITYLCM)
      DEALLOCATE(BSINFO2,BS2)
      ELSE IF(STYPE.EQ.1.AND.ISBS.EQ.1) THEN
      NBS=NBS1+NBS2
      ALLOCATE(BSINFO(3,NBS),BS(MAXL,NBS))
      BSINFO(:3,:NBS1)=BSINFO1(:3,:NBS1)
      BSINFO(:3,NBS1+1:NBS)=BSINFO2(:3,:NBS2)
      BS(:MAXL,:NBS1)=BS1(:MAXL,:NBS1)
      BS(:MAXL,NBS1+1:NBS)=BS2(:MAXL,:NBS2)
      JPSOUR=LCMLID(IPSOUR,'DBSOUR',2)
      KPSOUR_INFO=LCMLIL(JPSOUR,1,NBS)
      KPSOUR=LCMLIL(JPSOUR,2,NBS)
      DO N=1,NBS
        CALL LCMPDL(KPSOUR_INFO,N,3,1,BSINFO(1,N))
        CALL LCMPDL(KPSOUR,N,MAXL,2,BS(1,N))
      ENDDO
      CALL LCMPUT(IPSOUR,'NBS',1,1,NBS)
      DEALLOCATE(BSINFO,BSINFO1,BSINFO2,BS,BS1,BS2)      
      ENDIF

      ! Save the source normalization factor
      IF(BSNORM.LT.0) THEN
        CALL LCMPUT(IPSOUR,'NORM-FS',1,2,NORM)
      ELSE
        CALL LCMPUT(IPSOUR,'NORM-FS',1,2,NORM+BSNORM)
      ENDIF
*----
*  SAVE THE SIGNATURE AND STATE VECTOR
*----

      HSIGN='L_SOURCE'
      CALL LCMPTC(IPSOUR,'SIGNATURE',12,1,HSIGN)
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NG
      ISTATE(2)=NUNS
      ISTATE(3)=NDIR
      ISTATE(4)=NCST
      CALL LCMPUT(IPSOUR,'STATE-VECTOR',NSTATE,1,ISTATE)
      RETURN
*
      END
