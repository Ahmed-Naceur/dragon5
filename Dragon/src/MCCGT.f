*DECK MCCGT
      SUBROUTINE MCCGT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Adapt EXCELL tracking to MCCG requirements.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and R. Le Tellier
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1) modification type(L_TRACK);
*         HENTRY(2) sequential binary tracking file;
*         HENTRY(3) read-only type(L_GEOM).
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
      TYPE(C_PTR) IPTRK,IPGEO
      INTEGER NSTATE,MXGAUS,IOUT,IBCV
      PARAMETER (NSTATE=40,MXGAUS=64,IOUT=6,IBCV=-7)
      INTEGER ITRK,IFTR,IGEO,IFTRAK,J,NCOMNT,NBTR,ICOM,
     1 NDIM,ISPEC,N2REG,N2SOU,NALBG,NCOR,NANGL,MXSEG,NREG,NSOU,NANIS,
     2 ISYMM,IMPX,LCACT,NMU,MAXI,IAAC,ISCR,KRYL,IDIFC,ILEXA,ILEXF,INDIC,
     3 NITMA,NZP,N2RS,DIMKEYF,TYPOR1,TYPOR2,LTMT,STIS,LMXMCU,TRTY,PACA,
     4 SSYM,H,IMU,NFI,LMCU,N3MAX,ILINE,IANGL,N2SEG,NSEG,LMCU0,NLONG,
     5 LPS,NFIRST,NLEV,IK,JK,K,ILAST,IH,KJ,IPOS,IJEND,IJ,NFUNL,NUN,IA,
     6 IR,IKEY,ICUR,NPJJM,IFORW,II,I,IBIHET,IQUA10,IR2,NREG2,IFMT,NSUB,
     7 MXSUB,NMOD,NLIN,IE
      INTEGER ISOU,IDIM,IDIR
      REAL EPSI,HDD,TMUIM,FACSYM,DELU,FLOTT,DUM
      DOUBLE PRECISION WEI2D,DFLOTT,CMU
      CHARACTER TEXT4*4,TEXT12*12,TITLE*72,HSIGN*12,CFTRAK*12,
     1 COMNT(10)*80
      LOGICAL LPRISM,ACFLAG,LACA,LSCR,CYCLIC,LBIHET
      INTEGER IGP(NSTATE),KTITL(18),NCODE(6),IGB(8)
      REAL ZREAL(4),ZMU(MXGAUS),WZMU(MXGAUS),XMU(MXGAUS),ALBEDO(6),
     1 EXTKOP(NSTATE),XMU0(2*MXGAUS),WZMU0(2*MXGAUS)
      DOUBLE PRECISION CMUV(MXGAUS),CMUIV(MXGAUS),SMUV(MXGAUS),
     1 SMUIV(MXGAUS),TMUV(MXGAUS),TMUIV(MXGAUS)
*----
*  ALLOCATABLE STATEMENTS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INDREG,NZON,NZONA,ITEMP,
     1 MCUW,MCUI,NRSEG,KANGL,INOM3D,KM,MCU,IM,IS,JS,IPI,INVPI,LEV,LEVPT,
     2 KMROR,MCUROR,IMROR,JU,IWORK,IM0,MCU0,KEYFLX,KEYCUR,KEYANI,MAT
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISGNR
      REAL, ALLOCATABLE, DIMENSION(:) :: ZZ,VV,RTEMP,VA,VOL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DENSTY,SEGLEN,T2D,
     1 H3D,SURFD,VNUM
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: CAZ,XSIXYZ
*----
*  DATA STATEMENTS
*----
      INTEGER FACMCU(3)
      DATA FACMCU / 2,8,12 /
*----
*  PARAMETER VALIDATION
*----
      ITRK=0
      IFTR=0
      IGEO=0
      IF(NENTRY.LE.1) CALL XABORT('MCCGT: two PARAMETERS EXPECTED.')
*     tracking table in modification mode
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))
     1  CALL XABORT('MCCGT: LINKED LIST EXPECTED AT LHS.')
      ITRK=1
      IF(JENTRY(ITRK).NE.1) 
     1  CALL XABORT('MCCGT: ENTRY IN MODIFICATION MODE EXPECTED(1).')
      IF(IENTRY(2).EQ.3) THEN
*     tracking file in read-only mode
         IFTR=2
         IF(JENTRY(IFTR).NE.2) 
     1     CALL XABORT('MCCGT: ENTRY IN READ-ONLY MODE EXPECTED(1).')
      ELSE
         CALL XABORT('MCCGT: INVALID OR MISSING ENTRY(1)')
      ENDIF
      IF(NENTRY.GE.3) THEN
         IF(IENTRY(3).LE.2) THEN
*        geometry table in read-only mode
            IGEO=3
            IF (JENTRY(IGEO).NE.2) 
     1        CALL XABORT('MCCGT: ENTRY IN READ-ONLY MODE EXPECTED(2).')
         ELSE
            CALL XABORT('MCCGT: INVALID OR MISSING ENTRY(2)')
         ENDIF
      ENDIF
*
      IPTRK=KENTRY(ITRK)
      IFTRAK=FILUNIT(KENTRY(IFTR))
      IF(IGEO.NE.0) THEN
        IPGEO=KENTRY(IGEO)
      ELSE
        IPGEO=C_NULL_PTR
      ENDIF
*      
      CALL LCMGTC(IPTRK,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_TRACK') THEN
         TEXT12=HENTRY(ITRK)
         CALL XABORT('MCCGT: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_TRACK EXPECTED.')
      ENDIF
      CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,HSIGN)
      IF(HSIGN.NE.'EXCELL') THEN
         TEXT12=HENTRY(ITRK)
         CALL XABORT('MCCGT: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. EXCELL EXPECTED.')
      ENDIF
*----
*  RECOVER GEOMETRY
*----
      IF(C_ASSOCIATED(IPGEO)) THEN
        CALL LCMGTC(IPGEO,'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_GEOM') THEN
           TEXT12=HENTRY(IGEO)
           CALL XABORT('MCCGT: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1          '. L_GEOM EXPECTED.')
        ENDIF
        TEXT12=HENTRY(IGEO)
        CALL LCMPTC(IPTRK,'LINK.GEOM',12,1,TEXT12)
      ENDIF
*----
*  RECOVER SEQUENTIAL BINARY TRACKING FILE CHARACTERISTICS
*----
      CFTRAK=HENTRY(IFTR)
      CALL LCMPTC(IPTRK,'LINK.FTRACK',12,1,CFTRAK)
      REWIND IFTRAK
      READ(IFTRAK) TEXT4,NCOMNT,NBTR,IFMT
      DO ICOM=1,NCOMNT
         READ(IFTRAK) COMNT(ICOM)
      ENDDO
      READ(IFTRAK) NDIM,ISPEC,N2REG,N2SOU,NALBG,NCOR,NANGL,MXSUB,MXSEG
      IF((NDIM.NE.2).AND.(NDIM.NE.3)) 
     & CALL XABORT('2D OR 3D EXCELT TRACKING EXPECTED')
*----
*  RECOVER TRACKING STATE-VECTOR AND USER INPUT INFORMATION
*----
      CALL XDISET(IGP,NSTATE,0) 
      CALL LCMGET(IPTRK,'STATE-VECTOR',IGP)
      CALL LCMGET(IPTRK,'ALBEDO',ALBEDO)
      NREG=IGP(1)
      NSOU=IGP(5)
      NANIS=IGP(6)
      TRTY=IGP(9)
      CYCLIC=(TRTY.EQ.1)
      ISYMM=IGP(12)
*
      IMPX=1
      LCACT=IGP(13)
      NMU=IGP(14)
      LBIHET=(IGP(40).NE.0)
      MAXI=20
      IAAC=1
      ISCR=0
      KRYL=10
      IDIFC=0
      EPSI=1.0E-5
      HDD=0.0
      PACA=3
      ILEXA=0
      LTMT=0
      ILEXF=0
      STIS=0
      LMXMCU=0
      IFORW=0
      NFUNL=1
      NLIN=1
      DELU=0.0
      FACSYM=0.0
      IF(NANIS.LE.4) STIS=1
*----
*  PROCESS DOUBLE HETEROGENEITY (BIHET) DATA (IF AVAILABLE)
*----
      IF(LBIHET) THEN
         IF(.NOT.C_ASSOCIATED(IPGEO)) CALL XABORT('MCCGT: NO RHS GEOME'
     >   //'TRY DEFINED.')
         CALL LCMSIX(IPTRK,'BIHET',1)
         CALL LCMGET(IPTRK,'PARAM',IGB)
         IR2=IGB(2)
         NREG2=IGB(3)
         IBIHET=IGB(6)
         IQUA10=IGB(8)
         ALLOCATE(MAT(NREG),VOL(NREG))
         CALL LCMGET(IPTRK,'IBI',MAT)
         CALL LCMGET(IPTRK,'VOLUME',VOL)
         CALL LCMSIX(IPTRK,' ',2)
         CALL LCMPUT(IPTRK,'MATCOD',NREG,1,MAT)
         CALL LCMPUT(IPTRK,'VOLUME',NREG,2,VOL)
         DEALLOCATE(VOL,MAT)
         IGP(1)=NREG2
         IGP(2)=IGP(2)-(NREG-NREG2)
         IGP(4)=IR2
         CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,IGP)
         NREG=NREG2
      ENDIF
*
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
   20 IF(INDIC.EQ.10) GO TO 30
      IF(INDIC.NE.3) CALL XABORT('MCCGT: CHARACTER DATA EXPECTED(1).')
      IF(TEXT4.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MCCGT: INTEGER DATA EXPECTED(1).')
      ELSE IF(TEXT4.EQ.'GAUS') THEN
         LCACT=-1
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) GO TO 20
         NMU=NITMA
      ELSE IF(TEXT4.EQ.'DGAU') THEN
         LCACT=0
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) GO TO 20
         NMU=NITMA
      ELSE IF(TEXT4.EQ.'CACA') THEN
         LCACT=1
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) GO TO 20
         NMU=NITMA
      ELSE IF(TEXT4.EQ.'CACB') THEN
         LCACT=2
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) GO TO 20
         NMU=NITMA
      ELSE IF(TEXT4.EQ.'LCMD') THEN
         LCACT=3
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) GO TO 20
         NMU=NITMA
      ELSE IF(TEXT4.EQ.'OPP1') THEN
         LCACT=4
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) GO TO 20
         NMU=NITMA
      ELSE IF(TEXT4.EQ.'OGAU') THEN
         LCACT=5
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) GO TO 20
         NMU=NITMA
      ELSE IF(TEXT4.EQ.'EPSI') THEN
*        CONVERGENCE CRITERION FOR INNER ITERATIONS.
         CALL REDGET(INDIC,NITMA,EPSI,TEXT4,DFLOTT)
         IF(INDIC.NE.2) CALL XABORT('MCCGT: REAL DATA EXPECTED(1).')
      ELSE IF(TEXT4.EQ.'SCR') THEN
*        SCR ACCELERATION FLAG.
         CALL REDGET(INDIC,ISCR,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MCCGT: INTEGER DATA EXPECTED(2).')
      ELSE IF(TEXT4.EQ.'KRYL') THEN
*        GMRES FLAG.
         CALL REDGET(INDIC,KRYL,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MCCGT: INTEGER DATA EXPECTED(3).')
      ELSE IF(TEXT4.EQ.'AAC') THEN
*        ACA ACCELERATION FLAG.
         CALL REDGET(INDIC,IAAC,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MCCGT: INTEGER DATA EXPECTED(4).')
      ELSE IF(TEXT4.EQ.'BICG') THEN
*        ACA SYSTEM RESOLUTION TYPE (obsolete because it is the only option!)
         IF(IAAC.EQ.0) CALL XABORT('MCCGT: BICG ONLY IF ACA IS ON.')
      ELSE IF(TEXT4.EQ.'ILU0') THEN
*        ILU0 PRECONDITIONER FOR SOLVING ACA SYSTEM
         PACA=3
         IF(IAAC.EQ.0) CALL XABORT('MCCGT: ILU0 ONLY IF ACA IS ON.')
      ELSE IF(TEXT4.EQ.'DIAG') THEN
*        DIAGONAL PRECONDITIONER FOR SOLVING ACA SYSTEM
         PACA=1
         IF(IAAC.EQ.0) CALL XABORT('MCCGT: DIAG ONLY IF ACA IS ON.')
      ELSE IF(TEXT4.EQ.'FULL') THEN
*        FULL MATRIX PRECONDITIONER FOR SOLVING ACA SYSTEM
         PACA=2
         IF(IAAC.EQ.0) CALL XABORT('MCCGT: FULL ONLY IF ACA IS ON.')
      ELSE IF(TEXT4.EQ.'NONE') THEN
*        NO PRECONDITIONER FOR SOLVING ACA SYSTEM
         PACA=0
         IF(IAAC.EQ.0) CALL XABORT('MCCGT: NONE ONLY IF ACA IS ON.')
      ELSE IF(TEXT4.EQ.'TMT') THEN
*        TO USE A TRACK MERGING TECHNIQUE IN ACA CALCULATION
         LTMT=1
         IF(IAAC.EQ.0) CALL XABORT('MCCGT: LTMT ONLY IF ACA IS ON.')
      ELSE IF(TEXT4.EQ.'LEXA') THEN
*        TO FORCE EXACT EXPONENTIALS IN PRECONDITIONER CALCULATIONS
         ILEXA=1
      ELSE IF(TEXT4.EQ.'DIFC') THEN
*        TRANSPORT/DIFFUSION SOLUTION FLAG.
         IDIFC=1
         IAAC=1
      ELSE IF(TEXT4.EQ.'MCU') THEN
*        MAXIMUM DIMENSION OF MCU FOR MEMORY ALLOCATION.
         CALL REDGET(INDIC,LMXMCU,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MCCGT: INTEGER DATA EXPECTED(5).')
      ELSE IF(TEXT4.EQ.'HDD') THEN
*        SELECTION OD STEP CHARACTERISTICS METHOD.
         CALL REDGET(INDIC,NITMA,HDD,TEXT4,DFLOTT)
         IF(INDIC.NE.2) CALL XABORT('MCCGT: REAL DATA EXPECTED(2).')
      ELSE IF(TEXT4.EQ.'STIS') THEN
*        'SOURCE TERM ISOLATION' FLAG
         CALL REDGET(INDIC,STIS,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MCCGT: INTEGER DATA EXPECTED(6).')
         IF(ABS(STIS).GT.1) THEN
            CALL XABORT('MCCGT: STIS MUST BE SET TO -1, 0 OR 1.')
         ENDIF
      ELSE IF(TEXT4.EQ.'LEXF') THEN
*        TO FORCE EXACT EXPONENTIALS IN FLUX CALCULATIONS
         ILEXF=1
      ELSE IF(TEXT4.EQ.'ADJ') THEN
*        ADJOINT FLUX CALCULATION
         IFORW=1
      ELSE IF(TEXT4.EQ.'MAXI') THEN
*        MAXIMUM NUMBER OF INNER ITERATIONS.
         CALL REDGET(INDIC,MAXI,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MCCGT: INTEGER DATA EXPECTED(7).')
      ELSE IF(TEXT4.EQ.'SC') THEN
*        STEP CHARACTERISTICS OR DD0 SCHEME.
         NLIN=1
      ELSE IF(TEXT4.EQ.'LDC') THEN
*        LINEAR DISCONTINUOUS CHARACTERISTICS OR DD1 SCHEME.
         NLIN=3
         STIS=0
      ELSE IF(TEXT4.EQ.';') THEN
         GO TO 30
      ELSE 
         CALL XABORT('MCCGT: '//TEXT4//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 10
*
   30 IF(NMU.EQ.0) THEN
         IF(ISPEC.EQ.0) THEN
            IF(ISYMM.LE.1) THEN
               NMU=(NANGL+1)/2
            ELSE IF(ISYMM.GE.2) THEN
               NMU=NANGL
            ENDIF
         ELSE IF(ISPEC.EQ.1) THEN
            IF(ISYMM.LE.1) THEN
               NMU=(NANGL/4+1)/2
            ELSE IF(ISYMM.GE.2) THEN
               NMU=NANGL/4
            ENDIF
         ENDIF
      ENDIF
      IGP(13)=LCACT
      LACA=(IAAC.GT.0)
      LSCR=((ISCR.GT.0).AND.(.NOT.CYCLIC))
      ZREAL(1)=EPSI
      ZREAL(2)=HDD
      ACFLAG=((LSCR).OR.(LACA))
      NZP=IGP(39)
      LPRISM=(NZP.NE.0)
      IF(LPRISM) THEN
*     3D PRISMATIC GEOMETRY
         CALL LCMGET(IPTRK,'NCODE',NCODE)
         IF(NCODE(6).EQ.30) THEN
            IF(NCODE(5).EQ.30) THEN
*           Z- and Z+ surfaces symmetry
               SSYM=2
               FACSYM=0.0
            ELSE
*           Z+ symmetry
               SSYM=1
               FACSYM=1.0
            ENDIF
         ELSE
            SSYM=0
            FACSYM=0.0
         ENDIF
         N2RS=N2SOU+N2REG+1
         ALLOCATE(ZZ(NZP+1),INDREG(N2RS*(NZP+2)))
         CALL LCMSIX(IPTRK,'PROJECTION',1)
         CALL LCMGET(IPTRK,'IND2T3',INDREG)
         CALL LCMGET(IPTRK,'ZCOORD',ZZ)
         CALL LCMSIX(IPTRK,'PROJECTION',2)
         CALL LCMGET(IPTRK,'EXCELTRACKOP',EXTKOP)
         DELU=EXTKOP(40)
         ZREAL(3)=DELU
      ELSE
         ZREAL(3)=0.0
      ENDIF 
*
      IF(IMPX.GT.1) THEN
         CALL LCMGET(IPTRK,'TITLE',KTITL)
         WRITE(TITLE,'(3A4)') (KTITL(J),J=1,3)
         WRITE(IOUT,100) TITLE
         IF(LPRISM) WRITE(IOUT,*) '3D PRISMATIC EXTENDED TRACKING'
      ENDIF
*---
* CALCULATE POLAR QUADRATURE IF REQUIRED
*--- 
      TMUIM=0.0
      IF(NDIM.EQ.2) THEN
         IF(LCACT.EQ.-1) THEN   
            CALL ALGPT ( 2*NMU, -1.0, 1.0, XMU0, WZMU0)
            DO IMU=1,NMU
               XMU(NMU-IMU+1)=XMU0(IMU)
               WZMU(NMU-IMU+1)=WZMU0(IMU)
            ENDDO
         ELSE IF(LCACT.EQ.0) THEN   
            CALL ALGPT ( NMU, 0.0, 1.0, XMU, WZMU)
         ELSE
            IF(LCACT.GE.3) THEN
               IF(NMU.GT.4) NMU=4
               IF(NMU.LT.2) NMU=2
            ENDIF
            CALL ALCACT( LCACT, NMU, XMU, WZMU)
         ENDIF
         IF(LPRISM) THEN
            DO IMU=1,NMU
               ZMU(IMU)=1.0
               WZMU(IMU)=0.5*WZMU(IMU)
               CMU=DBLE(XMU(IMU))
               CMUV(IMU)=CMU
               CMUIV(IMU)=1.D0/CMU
               SMUV(IMU)=SQRT(1.D0-CMU**2)
               SMUIV(IMU)=1.D0/SMUV(IMU)
               TMUV(IMU)=SMUV(IMU)/CMU
               TMUIV(IMU)=1.D0/TMUV(IMU)
               TMUIM=MAX(REAL(TMUIM),REAL(TMUIV(IMU)))
            ENDDO
         ELSE
            DO IMU=1,NMU
               DUM=SQRT(1.-XMU(IMU)*XMU(IMU))
               ZMU(IMU)=1./DUM
               WZMU(IMU)=WZMU(IMU)*DUM
            ENDDO
         ENDIF
         IF(IMPX.GT.3) THEN
            CALL PRINAM('XMU   ',XMU,NMU)
            CALL PRINAM('ZMU   ',ZMU,NMU)
            CALL PRINAM('WZMU  ',WZMU,NMU)
         ENDIF
      ELSE ! NDIM.EQ.3
         NMU=1
         WZMU(1)=1.0
         XMU(1)=0.0
         ZMU(1)=1.0
      ENDIF
      IF(LPRISM) THEN
         CALL LCMSIX(IPTRK,'PROJECTION',1)
         CALL LCMPUT(IPTRK,'CMU',NMU,4,CMUV)
         CALL LCMPUT(IPTRK,'CMUI',NMU,4,CMUIV)
         CALL LCMPUT(IPTRK,'SMU',NMU,4,SMUV)
         CALL LCMPUT(IPTRK,'SMUI',NMU,4,SMUIV)
         CALL LCMPUT(IPTRK,'TMU',NMU,4,TMUV)
         CALL LCMPUT(IPTRK,'TMUI',NMU,4,TMUIV)
         CALL LCMSIX(IPTRK,'PROJECTION',2)
      ENDIF
      CALL LCMPUT(IPTRK,'ZMU$MCCG',NMU,2,ZMU)
      CALL LCMPUT(IPTRK,'XMU$MCCG',NMU,2,XMU)
      CALL LCMPUT(IPTRK,'WZMU$MCCG',NMU,2,WZMU)

      NFI=NREG+NSOU
      ALLOCATE(VV(NFI+1),NZON(NFI+1),ITEMP(NSOU),RTEMP(NSOU))
      IF(ACFLAG) ALLOCATE(NZONA(NFI+1))
*---
* RECOVER VOLUME AND MATALB ARRAYS
*---
      IF(LPRISM) THEN
*     3D PRISMATIC GEOMETRY
         READ(IFTRAK)
         READ(IFTRAK)
         CALL LCMSIX(IPTRK,'PROJECTION',1)
         CALL LCMGET(IPTRK,'VOLSUR',VV)
         CALL LCMGET(IPTRK,'MATALB',NZON)
         CALL LCMSIX(IPTRK,'PROJECTION',2)
      ELSE
*     REGULAR 2D OR 3D GEOMETRY
         READ(IFTRAK) (VV(NSOU+II+1),II=-NSOU,NREG)
         READ(IFTRAK) (NZON(NSOU+II+1),II=-NSOU,NREG)
      ENDIF
*---
* REORDER VOLUME AND MATALB ARRAYS (SURFACE -j BECOMES NV+j)
*     [S(-NS) S(-NS+1) ... S(-1) 0 V(1) ... V(NV-1) V(NV)]
*     ->[V(1)... V(NV-1) V(NV) S(-1) ... S(-NS+1) S(-NS)]
*---
      DO I=1,NSOU
         ITEMP(I)=NZON(NSOU-I+1)
         RTEMP(I)=VV(NSOU-I+1)
      ENDDO
      DO I=1,NREG
         VV(I)=VV(NSOU+I+1)
         NZON(I)=NZON(NSOU+I+1)
         IF(ACFLAG) NZONA(I)=NZON(I)
      ENDDO
      DO I=1,NSOU
         NZON(NREG+I)=ITEMP(I)
         IF(ACFLAG) THEN
            IF(ALBEDO(-ITEMP(I)).GT.0.0) THEN
               NZONA(NREG+I)=ITEMP(I)
            ELSE
               NZONA(NREG+I)=IBCV
            ENDIF
         ENDIF
         VV(NREG+I)=RTEMP(I)
      ENDDO
      DEALLOCATE(RTEMP,ITEMP)
*
      ALLOCATE(DENSTY(NANGL),CAZ(NDIM,NANGL))
      READ(IFTRAK)
      READ(IFTRAK)
      READ(IFTRAK) ((CAZ(IDIM,II),IDIM=1,NDIM),II=1,NANGL)
      READ(IFTRAK) (DENSTY(II),II=1,NANGL)
      IF(LPRISM) NDIM=3
*---
* CALCULATE NUMERICAL SURFACES FOR NON CYCLIC TRACKING
* IF AAC OR SCR USED, CALCULATE CONNECTION MATRICES
*---
      IF(ACFLAG) THEN
         IF(LMXMCU.EQ.0) LMXMCU=FACMCU(NDIM)*NFI
         ALLOCATE(MCUW(LMXMCU),MCUI(LMXMCU))
         CALL XDISET(MCUW,LMXMCU,0)
         CALL XDISET(MCUI,LMXMCU,0)
      ENDIF
      ALLOCATE(SEGLEN(MXSEG),NRSEG(MXSEG),KANGL(MXSUB))
      ALLOCATE(SURFD(NSOU),XSIXYZ(NSOU,3))
      CALL XDDSET(SURFD,NSOU,0.D0)
      CALL XDDSET(XSIXYZ,NSOU*3,0.D0)
      LMCU=NFI
      IF(LPRISM) THEN
*       3D PRISMATIC GEOMETRY: 3D TRACKS ARE RECONSTRUCTED
        ALLOCATE(VNUM(2*NREG*NMU*NANGL))
        CALL XDDSET(VNUM,2*NREG*NANGL*NMU,0.D0)
        ALLOCATE(T2D(MXSEG))
        N3MAX=(INT(FACSYM)+1)*MXSEG*(NZP+2)
        IF(SSYM.LT.2) THEN
          ALLOCATE(INOM3D(N3MAX),H3D(N3MAX))
        ELSE
          TMUIM=TMUIM/ZZ(NZP+1)
        ENDIF
        DO ILINE=1,NBTR
          READ(IFTRAK) NSUB,N2SEG,WEI2D,(KANGL(II),II=1,NSUB),
     1        (NRSEG(II),II=1,N2SEG),(SEGLEN(II),II=1,N2SEG)
          IF(NSUB.GT.MXSUB) CALL XABORT('MCCGT: MXSUB OVERFLOW.')
          IANGL=KANGL(1)
          IF(N2SEG.GT.0) THEN
            T2D(1)=0.0D0
            DO II=1,N2SEG-1
               T2D(II+1)=T2D(II)+SEGLEN(II+1)
            ENDDO
            IF(SSYM.EQ.2) THEN
               FACSYM=MAX(TMUIM*REAL(T2D(N2SEG)),FACSYM)
               ALLOCATE(INOM3D((INT(FACSYM)+1)*N3MAX),
     1                     H3D((INT(FACSYM)+1)*N3MAX))
            ENDIF
!!!!            IF(N2SEG-2.GE.3) then
!!!!            do IZP=0,NZP
!!!!            IF(IZP.lt.NZP) write(8,900) T2D,T2D,
!!!!     1                                   ZZ(IZP+1),ZZ(IZP+2)
!!!!            do II=1,N2SEG-2
!!!!            write(8,900) T2D(II),T2D(II+1),
!!!!     1                   ZZ(IZP+1),ZZ(IZP+1)
!!!!            IF(IZP.lt.NZP) write(8,900) T2D(II+1),T2D(II+1)
!!!!     1                                  ,ZZ(IZP+1),ZZ(IZP+2)
!!!!            enddo
!!!!            enddo
!!!! 900        FORMAT(
!!!!     1  7H line([,E16.8,1H,,E16.8,3H],[,E16.8,1H,,E16.8,2H],,
!!!!     2 26H 'Color','r','Marker','o')/)
            CALL MCGPTV(N2SOU,N2REG,NZP,SSYM,NREG,NSOU,N2SEG,N2SEG-2,
     1           NANGL,NMU,LMCU,LMXMCU,IANGL,INDREG,NRSEG,MCUW,MCUI,ZZ,
     2           T2D,WEI2D,CMUV,CMUIV,SMUV,SMUIV,TMUV,TMUIV,WZMU,DELU,
     3           INOM3D,H3D,SURFD,VNUM,ACFLAG)
!!!!            endif
            IF(SSYM.EQ.2) DEALLOCATE(H3D,INOM3D)
          ENDIF
        ENDDO
        IF(SSYM.LT.2) DEALLOCATE(H3D,INOM3D)
        DEALLOCATE(T2D)
        CALL MCGPTN(IMPX,NREG,NSOU,NANGL,NMU,VV,VNUM,SURFD,DENSTY,WZMU)
        IF(IMPX.GT.4) CALL PRINDM('VNORF',VNUM,2*NREG*NANGL*NMU)
        CALL LCMSIX(IPTRK,'PROJECTION',1)
        CALL LCMPUT(IPTRK,'VNORF',2*NREG*NANGL*NMU,4,VNUM)
        CALL LCMSIX(IPTRK,'PROJECTION',2)
        DEALLOCATE(VNUM)
      ELSE
*       REGULAR 2D OR 3D GEOMETRY
        DO ILINE=1,NBTR
          READ(IFTRAK) NSUB,NSEG,WEI2D,(KANGL(II),II=1,NSUB),
     1        (NRSEG(II),II=1,NSEG),(SEGLEN(II),II=1,NSEG)
          IF(NSUB.GT.MXSUB) CALL XABORT('MCCGT: MXSUB OVERFLOW.')
          IANGL=KANGL(1)
          IF(NSEG.GT.0) THEN
            CALL MCGDTV(NDIM,NFI,NREG,NSOU,NSEG,NMU,LMCU,LMXMCU,
     1           NZONA,NRSEG,MCUW,MCUI,WEI2D,SEGLEN,WZMU,SURFD,
     2           CYCLIC,ACFLAG,ZMU,XSIXYZ,CAZ(1,IANGL))
          ENDIF
        ENDDO
        IF(.NOT.CYCLIC) THEN
          CALL XDRSDB(NSOU,VV(NREG+1),SURFD,1)
          DO IDIR=1,3
            DO ISOU=1,NSOU
              XSIXYZ(ISOU,IDIR)=XSIXYZ(ISOU,IDIR)/SURFD(ISOU)
            ENDDO
          ENDDO
          CALL LCMPUT(IPTRK,'XSI$MCCG',NSOU*3,4,XSIXYZ)
        ENDIF
      ENDIF
*
      DEALLOCATE(XSIXYZ,SURFD)
      DEALLOCATE(KANGL,NRSEG,SEGLEN,CAZ,DENSTY)
      IF(LPRISM) DEALLOCATE(INDREG,ZZ)
*---
* CREATE CONNECTION MATRICES IN KM/MCU FORMAT
*---
      LMCU0=0
      IF(ACFLAG) THEN
*        KM(i) is the number of non-diagonal element on row i
*        MCU gives the column indexes.
         ALLOCATE(KM(NFI),MCU(LMXMCU))
         CALL MCGREC(NFI,KM,MCUW,MCUI,MCU,LMCU,LMXMCU,0)
         DEALLOCATE(MCUI,MCUW)
         NLONG=NFI
         IF(CYCLIC) THEN
*        if cyclic tracking, only the volume related data are stored
            LMCU=0
            DO I=1,NREG
               LMCU=LMCU+KM(I)
            ENDDO
            NLONG=NREG
            IF(NSOU.EQ.0) NFI=NREG+1
         ENDIF
         ALLOCATE(IM(NLONG+1))
*        construct IM
*        containing number of sparse matrix elements in rows before I
*        so location of J-th non-zero in row # I is IM(I)+J i.e.
*        IM(K+1)=sum_i=1^K KM(i) where K in [1,NLONG]
         IM=0
         DO I=1,NLONG
            IM(I+1)=IM(I)+KM(I)
         ENDDO
         IF(LSCR) THEN
*---
* SCR ACCELERATION : CREATE INDEX FOR THE SURFACES NEIGHBORS
*---
            LPS=IM(NFI+1)-IM(NREG+1)
            ALLOCATE(IS(NSOU+1),JS(LPS))
            LPS=0
            DO I=NREG+1,NFI
               IS(I-NREG)=LPS
               DO J=IM(I)+1,IM(I+1)
                  IF(MCU(J).GT.0) THEN
                     LPS=LPS+1
                     JS(LPS)=MCU(J)
                  ENDIF
               ENDDO
            ENDDO
            IS(NSOU+1)=LPS
         ENDIF
         IF(LACA) THEN
*---
* ACA ACCELERATION : 
*---
            ALLOCATE(IPI(NFI),INVPI(NFI))
            IF(PACA.GE.2) THEN
               LMCU0=LMCU
               ALLOCATE(LEV(NLONG),LEVPT(NLONG+1),KMROR(NLONG),
     1         MCUROR(LMCU),IMROR(NLONG+1),JU(NLONG),VA(NFI))
               IF(PACA.EQ.3) ALLOCATE(IWORK(NFI),IM0(NFI+1),MCU0(LMCU0))
*              construct IPI permutation : old_index=IPI(new_index) or
*              F_new=F_old(IPI) reordering of the unknowns of the
*              corrective system for ilu0 preconditioner.
               NFIRST=1
               TYPOR1=0
               TYPOR2=0
               CALL RENUM(NLONG,LMCU,NFIRST,IM,MCU,TYPOR1,TYPOR2,NLEV,
     1         LEV,LEVPT,IPI)
               IF(CYCLIC) THEN
                  DO I=NLONG+1,NFI
                     IPI(I)=I
                  ENDDO
               ENDIF
*              reorder everything according to IPI
*              construct INVPI permutation : new_index=INVPI(old_index)
*              or F_old=F_new(INVPI)
               DO I=1,NFI
                  J=IPI(I)
                  INVPI(J)=I
               ENDDO
               DO I=1,NLONG
                  J=IPI(I)               
                  KMROR(I)=KM(J)
               ENDDO
               IMROR=0
               DO I=1,NLONG
                  IMROR(I+1)=IMROR(I)+KMROR(I)
               ENDDO
               DO I=1,NLONG
                  J=IPI(I)
                  IK=IMROR(I)
                  DO JK=IM(J)+1,IM(J+1)
                     IK=IK+1
                     IF(MCU(JK).GT.0) THEN
                        MCUROR(IK)=INVPI(MCU(JK))
                     ELSE
                        MCUROR(IK)=MCU(JK)
                     ENDIF
                  ENDDO
               ENDDO
*              sort each line by increasing column index
               DO I=1,NLONG
                  K=IMROR(I)+1
                  CALL SORTIN(KMROR(I),MCUROR(K))
               ENDDO
               DO I=1,NFI
                  J=IPI(I)
                  NZONA(I)=NZON(J)
                  VA(I)=VV(J)
                  IF(J.GT.NREG) THEN
                     IF(ALBEDO(-NZON(J)).EQ.0.0) NZONA(I)=IBCV
                  ENDIF
               ENDDO
               CALL XDISET(JU,NLONG,0)
               IF(PACA.EQ.3) THEN
                  CALL XDISET(IM0,(NLONG+1),LMCU)
                  CALL XDISET(IWORK,NLONG,0)
                  ILAST=0
                  LMCU0=0
               ENDIF
               DO 50 I=1,NLONG
*              construct JU (and IM0/MCU0 for optimized storage)
*              MCUROR(JU(i):IMROR(i+1)) corresponds to the upper triangular part of line i.
*              MCUROR(IMROR(i)+1:JU(i)-1) correspond to the lower triangular part of line i.
                  DO IH=IMROR(I)+1,IMROR(I+1)
                     H=MCUROR(IH)
                     IF(H.GT.0) THEN
                        IF((H.GT.I).AND.(JU(I).EQ.0)) JU(I)=IH
                        IF(PACA.EQ.3) IWORK(H)=IH
                     ENDIF
                  ENDDO
                  IF(JU(I).EQ.0) JU(I)=IMROR(I+1)+1
                  IF(PACA.EQ.3) THEN
                  DO IK=IMROR(I)+1,JU(I)-1
                     K=MCUROR(IK)
                     IF(K.GT.0) THEN
                     DO KJ=JU(K),IMROR(K+1)
                        J=MCUROR(KJ)
                        IF(IWORK(J).GT.0) THEN
                           IPOS=0
                           IJEND=MIN(IM0(I+1),LMCU0)
                           DO IJ=IM0(I)+1,IJEND
                           IF(MCU0(IJ).EQ.J) THEN
                              IPOS=IJ
                              GOTO 40
                           ENDIF
                           ENDDO
   40                      CONTINUE
                           IF(IPOS.EQ.0) THEN
                           IF(ILAST.NE.I) THEN
                           CALL XDISET(IM0(ILAST+1),(I-ILAST),LMCU0)
                           ILAST=I
                           ENDIF
                           LMCU0=LMCU0+1
                           MCU0(LMCU0)=J
                           ENDIF
                        ENDIF
                     ENDDO
                     ENDIF
                  ENDDO
                  DO IH=IMROR(I)+1,IMROR(I+1)
                     H=MCUROR(IH)
                     IF(H.GT.0) IWORK(H)=0
                  ENDDO
                  ENDIF
   50          CONTINUE
               IF(PACA.EQ.3) THEN
               IF(LMCU0.EQ.0) THEN
                  PACA=4 ! SPECIAL CASE WHEN THERE IS NO EXTRA-STORAGE FOR ILU0-ACA
               ELSE
                  CALL XDISET(IM0(ILAST+1),(NLONG+1-ILAST),LMCU0)
               ENDIF
               ENDIF
            ELSE
               DO I=1,NFI
                  IPI(I)=I
                  INVPI(I)=I
               ENDDO
            ENDIF
         ENDIF
      ENDIF
      IF(CYCLIC.AND.(NSOU.EQ.0)) THEN
         NZON(NREG+1)=-1
         NZONA(NREG+1)=-1
      ENDIF
      IF(IMPX.GT.3) THEN
         CALL PRINIM('MATALB',NZON,NFI)
         CALL PRINAM('VOLSUR',VV,NFI)
         IF(ACFLAG) THEN
            CALL PRINIM('MATALA',NZONA,NFI)
            WRITE(IOUT,'(16H MCGREC : LMCU =,I6)') LMCU
            CALL PRINIM('KM    ',KM,NLONG)
            CALL PRINIM('MCU   ',MCU,LMCU)
            IF((LACA).AND.(PACA.GE.2)) THEN
               CALL PRINIM('IPERM ',INVPI,NFI)
               CALL PRINIM('KMROR ',KMROR,NLONG)
               CALL PRINIM('MCUROR',MCUROR,LMCU)
               CALL PRINIM('JU    ',JU,NLONG)
               IF(PACA.GE.3) THEN
                  WRITE(IOUT,'(16H MCCGT : LMCU0 =,I6)') LMCU0
                  IF(LMCU0.GT.0) THEN
                     CALL PRINIM('IM0    ',IM0,NLONG+1)
                     CALL PRINIM('MCU0   ',MCU0,LMCU0)
                  ENDIF
               ENDIF
            ENDIF
            IF(LSCR) THEN
               CALL PRINIM('IS    ',IS,NSOU+1)
               CALL PRINIM('JS    ',JS,LPS)
            ENDIF
         ENDIF
      ENDIF
*     
      CALL LCMPUT(IPTRK,'NZON$MCCG',NFI,1,NZON)
      CALL LCMPUT(IPTRK,'MATCOD',NREG,1,NZON)
      CALL LCMPUT(IPTRK,'V$MCCG',NFI,2,VV)
      CALL LCMPUT(IPTRK,'VOLUME',NREG,2,VV)
      IF(ACFLAG) THEN
         IF(LACA) THEN
            CALL LCMPUT(IPTRK,'NZONA$MCCG',NFI,1,NZONA)
            IF(PACA.GE.2) THEN
               CALL LCMPUT(IPTRK,'VA$MCCG',NFI,2,VA)
               CALL LCMPUT(IPTRK,'KM$MCCG',NLONG,1,KMROR)
               CALL LCMPUT(IPTRK,'IM$MCCG',NLONG+1,1,IMROR)
               CALL LCMPUT(IPTRK,'MCU$MCCG',LMCU,1,MCUROR)
               CALL LCMPUT(IPTRK,'JU$MCCG',NLONG,1,JU)
               IF(PACA.EQ.3) THEN
                  CALL LCMPUT(IPTRK,'IM0$MCCG',NLONG+1,1,IM0)
                  CALL LCMPUT(IPTRK,'MCU0$MCCG',LMCU0,1,MCU0)    
               ENDIF
               IF(PACA.GE.3) DEALLOCATE(MCU0,IM0,IWORK)
               DEALLOCATE(VA,JU,IMROR,MCUROR,KMROR,LEVPT,LEV)
            ELSE
               CALL LCMPUT(IPTRK,'KM$MCCG',NLONG,1,KM)
               CALL LCMPUT(IPTRK,'IM$MCCG',NLONG+1,1,IM)
               CALL LCMPUT(IPTRK,'MCU$MCCG',LMCU,1,MCU) 
               CALL LCMPUT(IPTRK,'VA$MCCG',NLONG,2,VV)           
            ENDIF
            CALL LCMPUT(IPTRK,'INVPI$MCCG',NFI,1,INVPI)
            CALL LCMPUT(IPTRK,'PI$MCCG',NLONG,1,IPI)
            DEALLOCATE(INVPI,IPI)
         ENDIF
         IF(LSCR) THEN
            CALL LCMPUT(IPTRK,'IS$MCCG',NSOU+1,1,IS)
            CALL LCMPUT(IPTRK,'JS$MCCG',LPS,1,JS)
            DEALLOCATE(JS,IS)
         ENDIF
         DEALLOCATE(IM,MCU,KM,NZONA)
      ENDIF
      DEALLOCATE(NZON,VV)  
      IF(.NOT.LACA) LMCU=0
      IF(.NOT.LSCR) LPS=0
*---
* MODIFY KEYFLX FOR ANISOTROPIC SCATTERING
* CREATE KEYCUR
*---
      IF(NDIM.EQ.1) THEN
         NFUNL=NANIS
         NMOD=2
      ELSE IF(NDIM.EQ.2) THEN
         NFUNL=NANIS*(NANIS+1)/2
         NMOD=4
      ELSE ! NDIM.EQ.3
         NFUNL=NANIS*NANIS
         NMOD=8
      ENDIF
      DIMKEYF=NREG*NLIN*NFUNL

      IGP(2)=DIMKEYF
      TEXT12='MCCG'
      CALL LCMPTC(IPTRK,'TRACK-TYPE',12,1,TEXT12)
*     non-cyclic tracking -> MCCG used (else MOCC)
      IF(.NOT.CYCLIC) IGP(2)=IGP(2)+IGP(5)

      NUN=IGP(2)
      ALLOCATE(KEYFLX(DIMKEYF))
      IF(NLIN.EQ.1) THEN
         DO 65 IA=1,NFUNL
            DO 60 IR=1,NREG
               KEYFLX((IA-1)*NREG+IR)=(IA-1)*NREG+IR
 60         CONTINUE
 65      CONTINUE
      ELSE IF(NLIN.EQ.3) THEN
         DO 72 IA=1,NFUNL
            DO 71 IE=1,3
               DO 70 IR=1,NREG
                  KEYFLX((IA-1)*3*NREG+(IE-1)*NREG+IR)
     1                         =(IA-1)*3*NREG+(IE-1)*NREG+IR
 70            CONTINUE
 71         CONTINUE
 72      CONTINUE
      ENDIF
      IF(.NOT.CYCLIC) THEN
         ALLOCATE(KEYCUR(NSOU))
         IKEY=1
         ICUR=0
         DO I=1,NUN
            IF((KEYFLX(IKEY).NE.I).OR.(IKEY.GT.DIMKEYF)) THEN
               ICUR=ICUR+1
               IF(ICUR.GT.NSOU)
     1         CALL XABORT('MCCGT: INCORRECT NUMBER OF UNKNOWNS')
               KEYCUR(ICUR)=I
            ELSE
               IKEY=IKEY+1
            ENDIF
         ENDDO
         CALL LCMPUT(IPTRK,'KEYCUR$MCCG',NSOU,1,KEYCUR)
         DEALLOCATE(KEYCUR)      
      ENDIF
      CALL LCMPUT(IPTRK,'KEYFLX$ANIS',DIMKEYF,1,KEYFLX)
      CALL LCMPUT(IPTRK,'KEYFLX',NREG,1,KEYFLX(:NREG))
      DEALLOCATE(KEYFLX)
*---
* GENERATE ALL SIGNS FOR SPHERICAL HARMONICS
*---
      ALLOCATE(ISGNR(NMOD,NFUNL),KEYANI(NFUNL))
      CALL MOCIK3(NANIS-1,NFUNL,NMOD,ISGNR,KEYANI) 
      DEALLOCATE(ISGNR)
*---
* GENERATE INDEX FOR PJJ(NU'->NU) STORAGE
* IF 'SOURCE TERM ISOLATION' OPTION IS ON
*---
      IF((STIS.NE.0).OR.(ISCR.GT.0)) THEN
         CALL MCGPJJ(IPTRK,IMPX,NDIM,NANIS,NFUNL,NPJJM,KEYANI)
      ELSE
         NPJJM=1
      ENDIF
      DEALLOCATE(KEYANI)
*     
      IGP(14)=NMU
      IGP(16)=NDIM
      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,IGP)
*---
* GENERATE MCCG-STATE AND REAL-PARAM VECTORS
*---
      CALL XDISET(IGP,NSTATE,0)
      IGP(1)=LCACT
      IGP(2)=NMU
      IGP(3)=KRYL
      IGP(4)=IDIFC
      IGP(5)=MXSEG
      IGP(6)=LMCU
      IGP(7)=IAAC
      IGP(8)=ISCR
      IGP(9)=LPS
      IGP(10)=PACA
      IGP(11)=ILEXA
      IGP(12)=ILEXF
      IGP(13)=MAXI
      IGP(14)=LTMT
      IGP(15)=STIS
      IGP(16)=NPJJM
      IGP(17)=LMCU0
      IGP(18)=IFORW
      IGP(19)=NFUNL
      IGP(20)=NLIN
      CALL LCMPUT(IPTRK,'MCCG-STATE',NSTATE,1,IGP)
      ZREAL(4)=FACSYM
      CALL LCMPUT(IPTRK,'REAL-PARAM',4,2,ZREAL)
*
      IF(IMPX.GT.1) THEN
         CALL LCMGET(IPTRK,'MCCG-STATE',IGP)
         WRITE(IOUT,120) (IGP(I),I=1,11)
         WRITE(IOUT,130) (IGP(I),I=12,20)
         CALL LCMGET(IPTRK,'REAL-PARAM',ZREAL)
         WRITE(IOUT,140) (ZREAL(I),I=1,4)
      ENDIF
*----
*  PROCESS DOUBLE HETEROGENEITY (BIHET) DATA (IF AVAILABLE)
*----
      IF(LBIHET) THEN
         CALL LCMGET(IPTRK,'EXCELTRACKOP',EXTKOP)
         CALL XDRTBH(IPGEO,IPTRK,IQUA10,IBIHET,IMPX,EXTKOP(39))
      ENDIF
*
      IF(IMPX.GT.2) CALL LCMLIB(IPTRK)
      RETURN
*
  100 FORMAT(/
     1 44H MM      MM  CCCCC   CCCCC   GGGGG  TTTTTTTT/
     2 44H MMM    MMM CCCCCCC CCCCCCC GGGGGGG TTTTTTTT/
     4 41H MMMM  MMMM CC   CC CC   CC GG         TT/
     5 41H MM  MM  MM CC      CC      GG  GGG    TT/
     6 41H MM      MM CC      CC      GG  GGG    TT/
     7 41H MM      MM CC   CC CC   CC GG   GG    TT/
     8 41H MM      MM CCCCCCC CCCCCCC GGGGGGG    TT/
     9 41H MM      MM  CCCCC   CCCCC   GGGGG     TT/
     1 17H TRACKING TITLE: ,A72/)
  120 FORMAT(/
     1 55H STATE VECTOR RELATED TO THE METHOD OF CHARACTERISTICS:/
     2 7H LCACT ,I9,29H   (TYPE OF POLAR QUADRATURE)/
     1 7H NMU   ,I9,48H   (ORDER OF THE POLAR QUADRATURE IN 2D/1 IN 3D)/
     5 7H KRYL  ,I9,48H   (<0 Bi-CGSTAB SCHEME USED /0=KRYLOV SCHEMES N,
     6 30HOT USED/ >0=GMRES SCHEME USED)/
     7 7H IDIFC ,I9,39H   (0=TRANSPORT/1=CDD SOLUTION OF FLUX)/
     8 7H NMAX  ,I9,42H   (MAXIMUM NUMBER OF ELEMENTS IN A TRACK)/
     9 7H LMCU  ,I9,42H   (DIMENSION OF MCU FOR ACA ACCELERATION)/
     3 7H IAAC  ,I9,48H   (0=NO ACCELERATION/1=CDD ACCELERATION OF INNE,
     4 13HR ITERATIONS)/
     2 7H SCR   ,I9,48H   (0=NO ACCELERATION/1=SCR ACCELERATION OF INNE,
     3 13HR ITERATIONS)/
     4 7H LPS   ,I9,42H   (DIMENSION OF PSJ FOR SCR ACCELERATION)/
     5 7H PACA  ,I9,48H   (PRECONDITIONER FOR SOLVING THE ACA SYSTEM WI,
     6 38HTH BICGSTAB (>2=ILU0, 1=DIAG, 0=NONE))/
     7 7H LEXA  ,I9,48H   (1=FORCE EXACT EXPONENTIAL USAGE IN PRECONDIT,
     8 18HIONER CALCULATION))
  130 FORMAT(
     1 7H LEXF  ,I9,48H   (1=FORCE EXACT EXPONENTIAL USAGE IN FLUX CALC,
     2 8HULATION)/
     3 7H MAXI  ,I9,39H   (MAXIMUM NUMBER OF INNER ITERATIONS)/
     4 7H LTMT  ,I9,48H   (TO USE TRACK MERGING FOR ACA SYSTEM CALCULAT,
     5 4HION)/
     6 7H STIS  ,I9,48H   (1=SOURCE TERM ISOLATION FOR FLUX INTEGRATION,
     7 1H)/
     8 7H NPJJM ,I9,48H   (NUMBER OF PJJ MODES TO STORE FOR STIS OPTION,
     9 1H)/
     1 7H LMCU0 ,I9,48H   (DIMENSION OF MCU0 FOR ILU0-ACA ACCELERATION)/
     2 7H IFORW ,I9,40H   (0/1=DIRECT/ADJOINT FLUX CALCULATION)/
     3 7H NFUNL ,I9,45H   (NUMBER OF SPHERICAL HARMONICS COMPONENTS)/
     4 7H NLIN  ,I9,43H   (1/3=SC OR DD0 SCHEME/LDC OR DD1 SCHEME))
  140 FORMAT(/
     1 12H REAL PARAM:/
     2 7H EPSI  ,1P,E12.4,33H   (TOLERANCE ON INNER ITERATION)/
     3 7H HDD   ,1P,E12.4,41H   (0.0=STEP CHARACTERISTICS SOLUTION/>0.,
     4 32H0=DIAMOND DIFFERENCING SOLUTION)/
     5 7H DELU  ,1P,E12.4,42H   (TRACK SPACING FOR 3D PRISMATIC GEOMETR,
     6 2HY)/
     7 7H FACSYM,1P,E12.4,42H   (TRACKING SYMMETRY FACTOR FOR MAXIMUM T,
     8 38HRACK LENGTH FOR 3D PRISMATIC GEOMETRY)/)
      END
