*DECK NAPPPR
      SUBROUTINE NAPPPR(IPMAP,IPTRK,IPFLU,IPMTX,IPMAC,NSTATE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform the Pin Power Reconstruction for core with
* heterogeneous mixture
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal.
*
*Author(s): 
* R. Chambon (EPM) and R. Nguyen Van Ho (URANUS)
*
*Parameters: input/output
* IPMAP   LCM object address of Map.
* IPTRK   LCM object address of Tracking.
* IPFLU   LCM object address of Flux.
* IPMTX   LCM object address of Matex.
* IPMAC   LCM object address of Macrolib of the fuel.
* NSTATE  length of the state vector
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NSTATE
      TYPE(C_PTR)  IPMAP,IPTRK,IPFLU,IPMTX,IPMAC
*----
*  LOCAL VARIABLES
*----
      INTEGER IOUT,NGPT
      REAL REPS
      PARAMETER (REPS=1.0E-4,IOUT=6,NGPT=2)
      TYPE(C_PTR) JPFLU,JPMAP,KPMAP
      INTEGER INDIC,NITMA,LENGTH,NBPIN
      CHARACTER TEXT*12
      REAL FLOT
      DOUBLE PRECISION DFLOT
      INTEGER ISTATE(NSTATE),IMPX,IMETH
      INTEGER NXP,NYP,NXD,NYD,NZD,NAX,NAY,
     1  NASS,NCOMB,NG,NASS2,NREG,NXM,NYM,NZM,
     2  NXT,NYT,NZT,NXDA,NYDA,NZDA,NCH,NZASS,NPIN,IFX,
     3  NUN,IEL,NMIX,NAMIX,NGFF
      CHARACTER LABEL*8
      CHARACTER TFDINF*12
      INTEGER I,J,K,IP,JP,I1,I2,J1,J2,K1,K2,IASS,ICH,IG,IM,JM,ID,JD,KM,
     1 IAX,JAX,IGP,JGP,KGP,IMIX,IPIN,ICHX,IDIM,LC,L4,MAXKN,MKN,ITYLCM,
     2 ITYPE
      REAL POW,FACT,POWTOT,POWASS,DX,DY,DZ,FPD,FQ,PMAX,
     1  HOTPINPOW,PINPOW,FXY,VTOT
      REAL ZGKSIX(NGPT),ZGKSIY(NGPT),ZGKSIZ(NGPT),WGKSIX(NGPT),
     1  WGKSIY(NGPT),WGKSIZ(NGPT),X(NGPT),Y(NGPT),Z(NGPT),
     2  FLUGP(NGPT,NGPT,NGPT)
      REAL E(25)
      LOGICAL LSPX,LSPY,LSPZ,LCH,LPOW,LNOINT,LDEBUG
*----
*  ALLOCATABLE ARRAYS
*----
      CHARACTER*4, ALLOCATABLE, DIMENSION(:) :: NAMX,NAMY
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBAX,IBAX,BMIXP,AZONE,
     1  ACOMB,KN
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: CODEA
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: KEYFLX,BMIX,MAT
      REAL, ALLOCATABLE, DIMENSION(:) :: MXP,MYP,MXD,MYD,MZD,MXM,
     1  MYM,MZM,MXDA,MYDA,MZDA,FLXD,VOLM,FXYZ,PLINMAXZ,FXYASS,
     2  FACTASS,PWASS,PWASS2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FXTD,FYTD,BUNDPW
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: HFA,HFM,FDINFM,
     1  FTINFM,AXPOW,VPIN
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: FLXDA,VOL
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLXP,HF,FDINF,FTINF
*
      IMPX=0
      FACT=1.0
      LSPX=.FALSE.
      LSPY=.FALSE.
      LSPZ=.FALSE.
      LPOW=.FALSE.
      LNOINT=.FALSE.
      NZASS=0
      NPIN=0
      NBPIN=0
      IFX=0
      POW=1.0
      FQ=0.0
      FXY=0.0
      HOTPINPOW=0.0
      PINPOW=0.0
      PMAX=0.0
      VTOT=0.0
      LDEBUG=.false.
* Read mandatory keywords
      if(LDEBUG)write(6,*) 'NAPPPR begin debug'
* [EDIT] PPR
      CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
      IF(INDIC.NE.3) CALL XABORT('NAPPPR: character data expected.')
      IF(TEXT.EQ.'EDIT') THEN
        CALL REDGET(INDIC,IMPX,FLOT,TEXT,DFLOT)
        IF(INDIC.NE.1) CALL XABORT('NAPPPR: inteGEr data expected.')
        CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
        IF(INDIC.NE.3) CALL XABORT('NAPPPR: character data expected.')
      ENDIF
      IF(TEXT.NE.'PPR') CALL XABORT('NAPPPR: ''PPR'' keyword '//
     1  'expected.')
!* NPIN
!      CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
!      IF(INDIC.NE.3) CALL XABORT('NAPPPR: character data expected.')
!      IF(TEXT.NE.'NPIN') CALL XABORT('NAPPPR: ''NPIN'' keyword '//
!     1  'expected.')
!      CALL REDGET(INDIC,NPIN,FLOT,TEXT,DFLOT)
!      IF(INDIC.NE.1) CALL XABORT('NAPPPR: NPIN inteGEr expected.')
!      NXP=NPIN
!      NYP=NPIN
* NZASS
      CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
      IF(INDIC.NE.3) CALL XABORT('NAPPPR: character data expected.')
      IF(TEXT.NE.'NZASS') CALL XABORT('NAPPPR: ''NZASS'' keyword '//
     1  'expected.')
      CALL REDGET(INDIC,NZASS,FLOT,TEXT,DFLOT)
      IF(INDIC.NE.1) CALL XABORT('NAPPPR: NZASS inteGEr expected.')
C* SPIN + SGAP
C      CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
C      IF(INDIC.NE.3) CALL XABORT('NAPPPR: character data expected.')
C      IF(TEXT.NE.'DIM') CALL XABORT('NAPPPR: ''NZASS'' keyword '//
C     1  'expected.')
C      CALL REDGET(INDIC,NITMA,SPIN,TEXT,DFLOT)
C      IF(INDIC.NE.2) CALL XABORT('NAPPPR: SPIN real expected.')
C      CALL REDGET(INDIC,NITMA,SGAP,TEXT,DFLOT)
C      IF(INDIC.NE.2) CALL XABORT('NAPPPR: SGAP real expected.')
C      
* GEt core GEometry description in matex
      IF(IMPX.GE.100)WRITE(6,*) 'debug:GEt matex info'
      CALL LCMGET(IPMTX,'STATE-VECTOR',ISTATE)
      NG=ISTATE(1)
      NXD=ISTATE(8)
      NYD=ISTATE(9)
      NZD=ISTATE(10)
      ALLOCATE(MXD(NXD+1),MYD(NYD+1),MZD(NZD+1))
      CALL LCMGET(IPMTX,'MESHX',MXD)
      CALL LCMGET(IPMTX,'MESHY',MYD)
      CALL LCMGET(IPMTX,'MESHZ',MZD)
* GEt KEYFLX
      IF(IMPX.GE.100)WRITE(6,*) 'debug:GEt track info'
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NREG=ISTATE(1)
      NUN=ISTATE(2)
      ITYPE=ISTATE(6)
      IEL=ISTATE(9)
      L4=ISTATE(11)
      ICHX=ISTATE(12)
      NXT=ISTATE(14)
      NYT=ISTATE(15)
      NZT=ISTATE(16)
      IDIM=1
      IF((ITYPE.EQ.5).OR.(ITYPE.EQ.6).OR.(ITYPE.EQ.8)) IDIM=2
      IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) IDIM=3
      IF((NXD.NE.NXT).OR.(NYD.NE.NYT).OR.(NZD.NE.NZT)) CALL XABORT
     1  ('NAPPPR: dimension do not match between MATEX and TRACKING')
      ALLOCATE(KEYFLX(NXT,NYT,NZT),MAT(NXT,NYT,NZT))
      CALL LCMGET(IPTRK,'KEYFLX',KEYFLX)
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      ALLOCATE(FLXD(NUN))
* GEt assembly GEometry in map
      IF(IMPX.GE.100)WRITE(6,*) 'debug:GEt map info'
      CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
      NCH=ISTATE(2)
      NASS=ISTATE(14)
      NAX=ISTATE(15)
      NAY=ISTATE(16)
      ALLOCATE(AZONE(NCH))
      ALLOCATE(NAMX(NAX),NAMY(NAY))
      CALL LCMGET(IPMAP,'A-ZONE',AZONE)
      CALL LCMGTC(IPMAP,'AXNAME',4,NAX,NAMX)
      CALL LCMGTC(IPMAP,'AYNAME',4,NAY,NAMY)
      CALL LCMSIX(IPMAP,'GEOMAP',1)
      CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
      NXM=ISTATE(3)
      NYM=ISTATE(4)
      NZM=ISTATE(5)
      ALLOCATE(MXM(NXM+1),MYM(NYM+1),MZM(NZM+1))
      ALLOCATE(NBAX(NAY),IBAX(NAY))
      ALLOCATE(BMIX(NXM,NYM,NZM))
      CALL LCMGET(IPMAP,'MESHX',MXM)
      CALL LCMGET(IPMAP,'MESHY',MYM)
      CALL LCMGET(IPMAP,'MESHZ',MZM)
      CALL LCMLEN(IPMAP,'A-NX',LENGTH,INDIC)
      IF(LENGTH.NE.NAY) CALL XABORT('NAPPPR: Number of assembly along'
     1  //'Y direction do not match between MAP and embedded GEometry')
      CALL LCMGET(IPMAP,'A-NX',NBAX)
      CALL LCMGET(IPMAP,'A-IBX',IBAX)
      CALL LCMSIX(IPMAP,'GEOMAP',2)
      CALL LCMGET(IPMAP,'BMIX',BMIX)
C* GEt data in pin by pin assembly GEometry 
C      IF(IMPX.GE.100)WRITE(6,*) 'debug:GEt map pinBypin info'
C      CALL LCMGET(IPMPP,'STATE-VECTOR',ISTATE)
C      NCHP=ISTATE(2)
C      NASSP=ISTATE(14)
C*     total number of fuel bundles = tot. nb. of .XS 
C      NXS=ISTATE(9)
C      IF(NASS.NE.NASSP)CALL XABORT('NAPPPR: number of assembly do not '
C     1  //'match between unfolded GEometries')
C      ALLOCATE(AZONEP(NCHP))
C      CALL LCMGET(IPMPP,'A-ZONE',AZONEP)
C      CALL LCMSIX(IPMPP,'GEOMAP',1)
C      CALL LCMGET(IPMPP,'STATE-VECTOR',ISTATE)
C      NXMP=ISTATE(3)
C      NYMP=ISTATE(4)
C      NZMP=ISTATE(5)
C      CALL LCMSIX(IPMPP,'GEOMAP',2)
C      ALLOCATE(BMIXP(NXMP,NYMP,NZMP))
C      CALL LCMGET(IPMPP,'BMIX',BMIXP)
      IF(IMPX.GE.5) THEN
        WRITE(6,*) 'MATEX dimension   (het):',NXD,NYD,NZD
        WRITE(6,*) 'TRACKING dimension(het):',NXT,NYT,NZT
        WRITE(6,*) 'MAP dimension     (het):',NXM,NYM,NZM
      ENDIF
* Read remaining input file
      NCOMB=0
      ALLOCATE(ACOMB(NASS))
      IF(IMPX.GE.100)WRITE(6,*) 'debug: beg read input'
      IMETH=0
    5 CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
        IF(TEXT.EQ.'METH') THEN
        CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
        IF(INDIC.NE.3) CALL XABORT('NAPPPR: character data expected.')
        IF(TEXT.EQ.'GPPR') THEN
          IMETH=1
          CALL REDGET(INDIC,IFX,FLOT,TEXT,DFLOT)
          IF(INDIC.NE.1) CALL XABORT('NAPPPR: inteGEr data expected.')
          WRITE(TFDINF,500) IFX
        ELSE
          CALL XABORT('NAPPPR: '//TEXT//' is a wrong method keyword.')
        ENDIF
        GOTO 5
      ELSEIF(TEXT.EQ.'POWER') THEN
        LPOW=.TRUE.
        CALL REDGET(INDIC,NITMA,POW,TEXT,DFLOT)
        IF(INDIC.NE.2) CALL XABORT('NAPPPR: POWER real expected.')
        GOTO 5
      ELSEIF(TEXT.EQ.';') THEN
        GOTO 50
      ELSE
        CALL XABORT('NAPPPR: '//TEXT//' is a wrong keyword.')
      ENDIF
*-----------------------------
   50 CONTINUE
      IF(IMPX.GE.100)WRITE(6,*) 'debug: computation begin'
* Compute mesh X and Y for a pin-by-pin assembly
      CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
      NMIX=ISTATE(2)
      NAMIX=NMIX/NASS/NZASS
      IF(IMPX.GE.1) WRITE(6,*) 'Number of Mix per assembly per plane'//
     1 ' NAMIX = ',NAMIX
      NGFF=ISTATE(16)
      IF(NGFF.EQ.0) CALL XABORT('NAPPPR: NGFF.NE.0 expected.')
      CALL LCMSIX(IPMAC,'GFF',1)
      CALL LCMSIX(IPMAC,'GFF-GEOM',1)
      CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
      NXP=ISTATE(3)
      NYP=ISTATE(4)
      NPIN=NXP
      ALLOCATE(MXP(NXP+1),MYP(NYP+1))
      CALL LCMGET(IPMAC,'MESHX',MXP)
      CALL LCMGET(IPMAC,'MESHY',MYP)
      DO I=2,NXP+1
        MXP(I)=MXP(I)-MXP(1)
      ENDDO
      MXP(1)=0.0
      DO I=2,NYP+1
        MYP(I)=MYP(I)-MYP(1)
      ENDDO
      MYP(1)=0.0
      CALL LCMSIX(IPMAC,'GFF-GEOM',2)
      CALL LCMSIX(IPMAC,'GFF',2)
* Compute IX-,IX+,IY-,IY+,IZ-,IZ+ for each assembly in core GEometry
      IF(IMPX.GE.100) WRITE(6,*) 'debug PPR:IX-,IX+,IY-,IY+,IZ-,IZ+'
      ALLOCATE(CODEA(NASS,6))
      CALL XDISET(CODEA,NASS*6,0)
      ICH=0
      I1=0
      I2=0
      NASS2=0
      DO IASS=1,NASS
        CODEA(IASS,1)=NXD+1
        CODEA(IASS,2)=0
        CODEA(IASS,3)=NYD+1
        CODEA(IASS,4)=0
        CODEA(IASS,5)=NZD+1
        CODEA(IASS,6)=0
      ENDDO
      DO 150 JM=1,NYM
      DO 130 IM=1,NXM
      LCH=.TRUE.
      IASS=0
      DO 100 KM=1,NZM
      IF(BMIX(IM,JM,KM).NE.0) THEN
        IF(LCH) THEN
          ICH=ICH+1
          LCH=.FALSE.
          IASS=AZONE(ICH)
          NASS2=MAX(NASS2,IASS)
          DO I=1,NXD+1
            IF(MXD(I).EQ.MXM(IM)) I1=I
            IF(MXD(I).EQ.MXM(IM+1)) I2=I
          ENDDO
          CODEA(IASS,1)=MIN(I1,CODEA(IASS,1))
          CODEA(IASS,2)=MAX(I2,CODEA(IASS,2))
          DO I=1,NYD+1
            IF(MYD(I).EQ.MYM(JM)) I1=I
            IF(MYD(I).EQ.MYM(JM+1)) I2=I
          ENDDO
          CODEA(IASS,3)=MIN(I1,CODEA(IASS,3))
          CODEA(IASS,4)=MAX(I2,CODEA(IASS,4))
          DO I=1,NZD+1
            IF(MZD(I).EQ.MZM(KM)) I1=I
            IF(MZD(I).EQ.MZM(KM+1)) I2=I
          ENDDO
          CODEA(IASS,5)=MIN(I1,CODEA(IASS,5))
          CODEA(IASS,6)=MAX(I2,CODEA(IASS,6))
        ELSE
          DO I=2,NZD+1
            IF(MZD(I).EQ.MZM(KM+1)) I2=I
          ENDDO
          CODEA(IASS,6)=MAX(I2,CODEA(IASS,6))
        ENDIF
      ENDIF
  100 CONTINUE
  130 CONTINUE
  150 CONTINUE
      IF(IMPX.GE.10) THEN
        WRITE(6,*) 'Position of the assemblies in the core'
        WRITE(6,*) 'IX-,IX+,IY-,IY+,IZ-,IZ+'
        do iass=1,nass
          WRITE(6,*) 'Assembly #',iass,':',(CODEA(iass,i),i=1,6)
        ENDDO
      ENDIF
      IF(NASS2.NE.NASS)CALL XABORT('NAPPPR: number of assembly do not'
     1 //' match: NASS2.NE.NASS')
* For all assembly perform PPR 
      ALLOCATE(FLXP(NXP,NYP,NZASS,NG,NASS))
      ALLOCATE(AXPOW(NXP,NYP,NASS))
      ALLOCATE(VPIN(NXP,NYP,NASS))
      ALLOCATE(FXYASS(NASS))
      ALLOCATE(FACTASS(NASS))
      ALLOCATE(PWASS(NASS),PWASS2(NASS))
      ALLOCATE(BUNDPW(NASS,NZASS))
      IF(.NOT.LPOW) THEN
        CALL LCMGET(IPMAP,'BUND-PW',BUNDPW)
      ENDIF
      DO IASS=1,NASS
        PWASS(IASS)=0.0
        DO IP=1,NXP
          DO JP=1,NYP
            AXPOW(IP,JP,IASS)=0.0
            VPIN(IP,JP,IASS)=0.0
          ENDDO
        ENDDO
        FXYASS(IASS)=0.0
        IF(.NOT.LPOW) THEN
          DO K=1,NZASS
            PWASS(IASS)=PWASS(IASS)+BUNDPW(IASS,K)
          ENDDO
        ENDIF
      ENDDO
      DO IASS=1,NASS
*   GEt flux at core GEometry level for assembly only
      I1=CODEA(IASS,1)
      I2=CODEA(IASS,2)
      J1=CODEA(IASS,3)
      J2=CODEA(IASS,4)
      K1=CODEA(IASS,5)
      K2=CODEA(IASS,6)
      NXDA=I2-I1
      NYDA=J2-J1
      NZDA=K2-K1
      ALLOCATE(FLXDA(NXDA,NYDA,NZDA,NG))
      CALL XDRSET(FLXDA,NXDA*NYDA*NZDA*NG,0.0)
      IF(NZDA.NE.NZASS) CALL XABORT('NAPPPR: incoherent number of mesh' 
     1 //' in Z direction for an assembly: NZDA.NE.NZASS')
      JPFLU=LCMGID(IPFLU,'FLUX')
      IF((LNOINT).OR.(IMPX.GE.0)) THEN
        
        DO IG=1,NG
          CALL LCMGDL(JPFLU,IG,FLXD)
        DO I=I1,I2-1
        DO J=J1,J2-1
        DO K=K1,K2-1
          FLXDA(I-I1+1,J-J1+1,K-K1+1,IG)=FLXD(KEYFLX(I,J,K))
        ENDDO
C       end K
        ENDDO
C       end J
        ENDDO
C       end I
        ENDDO
C       end IG
      ENDIF
      ALLOCATE(MXDA(NXDA+1),MYDA(NYDA+1),MZDA(NZDA+1))
      DO I=I1,I2
        MXDA(I-I1+1)=MXD(I)-MXD(I1)+MXP(1)
      ENDDO
      IF(ABS(MXDA(NXDA+1)-MXDA(1)-MXP(NXP+1)+MXP(1)).GT.0.0001) THEN
        WRITE(6,*) 'Assembly Transport and Core meshX do not match:'// 
     1   'Transport=',MXP(NXP+1)-MXP(1),'Core=',MXDA(NXDA+1)-MXDA(1)
        CALL XABORT('Sizes do not match')
      ENDIF
      DO J=J1,J2
        MYDA(J-J1+1)=MYD(J)-MYD(J1)+MYP(1)
      ENDDO
      IF(ABS(MYDA(NYDA+1)-MYDA(1)-MYP(NYP+1)+MYP(1)).GT.0.0001) THEN
        WRITE(6,*) 'Assembly Transport and Core meshY do not match:'// 
     1   'Transport=',MYP(NYP+1)-MYP(1),'Core=',MYDA(NYDA+1)-MYDA(1)
        CALL XABORT('Sizes do not match')
      ENDIF
      DO K=K1,K2
        MZDA(K-K1+1)=MZD(K)
      ENDDO
      IF(IMPX.GE.10) THEN
        WRITE(6,*) 'Coarse Flux and mesh at assembly level'
        WRITE(6,*) 'Mesh X:',(MXDA(I),I=1,NXDA+1)
        WRITE(6,*) 'Mesh Y:',(MYDA(I),I=1,NYDA+1)
        WRITE(6,*) 'Mesh Z:',(MZDA(I),I=1,NZDA+1)
        WRITE(6,*) 'Flux:'
        DO IG=1,NG
        WRITE(6,*) 'Group #',IG
        DO K=1,NZDA
        WRITE(6,*) 'Plan #',K
        DO J=1,NYDA
        WRITE(6,*) (FLXDA(I,J,K,IG),I=1,NXDA)
        ENDDO
        ENDDO
        ENDDO
      ENDIF
*   project flux at assembly level
      ALLOCATE(FXTD(NXP,NXDA),FYTD(NYP,NYDA))
      CALL XDRSET(FXTD,NXP*NXDA,0.0)
      CALL XDRSET(FYTD,NYP*NYDA,0.0)
*   compute fraction of the transport volumes occupied by diffusion volumes  
      CALL NAPFTD(NXP,MXP,NXDA,MXDA,FXTD)
      CALL NAPFTD(NYP,MYP,NYDA,MYDA,FYTD)
!      DO IP=1,NXP
!      DXP=MXP(IP+1)-MXP(IP)
!      DO ID=1,NXDA
!        IF((MXDA(ID).LE.MXP(IP)).AND.(MXDA(ID+1).GE.MXP(IP+1))) THEN
!          FXTD(IP,ID)=1.0
!        ELSEIF ((MXDA(ID).LE.MXP(IP)).AND.(MXDA(ID+1).GT.MXP(IP))) THEN
!          FXTD(IP,ID)=(MXDA(ID+1)-MXP(IP))/DXP
!        ELSEIF ((MXDA(ID).GE.MXP(IP)).AND.
!     1          (MXDA(ID+1).LE.MXP(IP+1))) THEN
!          FXTD(IP,ID)=(MXDA(ID+1)-MXDA(ID))/DXP
!        ELSEIF ((MXDA(ID).LT.MXP(IP+1)).AND.
!     1          (MXDA(ID+1).GE.MXP(IP+1))) THEN
!          FXTD(IP,ID)=(MXP(IP+1)-MXDA(ID))/DXP
!        ENDIF
!      ENDDO
!      ENDDO
*
!      DO IP=1,NYP
!      DYP=MYP(IP+1)-MYP(IP)
!      DO ID=1,NYDA
!        IF((MYDA(ID).LE.MYP(IP)).AND.(MYDA(ID+1).GE.MYP(IP+1))) THEN
!          FYTD(IP,ID)=1.0
!        ELSEIF ((MYDA(ID).LE.MYP(IP)).AND.(MYDA(ID+1).GT.MYP(IP))) THEN
!          FYTD(IP,ID)=(MYDA(ID+1)-MYP(IP))/DYP
!        ELSEIF ((MYDA(ID).GE.MYP(IP)).AND.
!     1          (MYDA(ID+1).LE.MYP(IP+1))) THEN
!          FYTD(IP,ID)=(MYDA(ID+1)-MYDA(ID))/DYP
!        ELSEIF ((MYDA(ID).LT.MYP(IP+1)).AND.
!     1          (MYDA(ID+1).GE.MYP(IP+1))) THEN
!          FYTD(IP,ID)=(MYP(IP+1)-MYDA(ID))/DYP
!        ENDIF
!      ENDDO
!      ENDDO
!     adds up all fluxes
      if(LDEBUG)write(6,*)'NXP,NYP',NXP,NYP
      DO IG=1,NG
      IF(.NOT.LNOINT) CALL LCMGDL(JPFLU,IG,FLXD)
      DO K=1,NZASS
      DO IP=1,NXP
      DO JP=1,NYP
        FLXP(IP,JP,K,IG,IASS)=0.0
        DO ID=1,NXDA
        DO JD=1,NYDA
          IF(LNOINT) THEN
* No interpolation: use averaGE flux
            FLXP(IP,JP,K,IG,IASS)=FLXP(IP,JP,K,IG,IASS)
     1       +FLXDA(ID,JD,K,IG)*FXTD(IP,ID)*FYTD(JP,JD)
* Interpolate flux with polynomial representation
*   (only if pin and macro region have a non-nul intersection)
          ELSEIF(FXTD(IP,ID)*FYTD(JP,JD).NE.0.0) THEN
*     indent removed
*     compute gauss points and weights
      CALL ALGPT(NGPT,MAX(MXP(IP),MXDA(ID)),MIN(MXP(IP+1),MXDA(ID+1)),
     1   ZGKSIX,WGKSIX)
      DX=MIN(MXP(IP+1),MXDA(ID+1))-MAX(MXP(IP),MXDA(ID))
      CALL ALGPT(NGPT,MAX(MYP(JP),MYDA(JD)),MIN(MYP(JP+1),MYDA(JD+1)),
     1   ZGKSIY,WGKSIY)
      DY=MIN(MYP(JP+1),MYDA(JD+1))-MAX(MYP(JP),MYDA(JD))
      CALL ALGPT(NGPT,MZDA(K),MZDA(K+1),ZGKSIZ,WGKSIZ)
      DZ=MZDA(K+1)-MZDA(K)
      IF(IMPX.GE.10) then
        WRITE(6,*) 'IP,JP:',IP,JP,FXTD(IP,ID),'ID,JD:',ID,JD,FYTD(JP,JD)
        WRITE(6,*) 'Gauss point ZGWG:',(ZGKSIX(I),I=1,NGPT),
     1   (WGKSIX(I),I=1,NGPT),'DX',DX
        WRITE(6,*) 'Gauss point ZGWG:',(ZGKSIY(I),I=1,NGPT),
     1   (WGKSIY(I),I=1,NGPT),'DY',DY
        WRITE(6,*) 'Gauss point ZGWG:',(ZGKSIZ(I),I=1,NGPT),
     1   (WGKSIZ(I),I=1,NGPT),'DZ',DZ
      ENDIF
      
*     interpolate flux
      FPD=0.0
      DO IGP=1,NGPT
        X(IGP)=MXD(I1)-MXP(1)+ZGKSIX(IGP)
      ENDDO
      DO JGP=1,NGPT
        Y(JGP)=MYD(J1)-MYP(1)+ZGKSIY(JGP)
      ENDDO
      DO KGP=1,NGPT
        Z(KGP)=ZGKSIZ(KGP)
      ENDDO
      IF(IMPX.GE.10) then
        WRITE(6,*) 'Gauss point X:',(X(I),I=1,NGPT)
        WRITE(6,*) 'Gauss point Y:',(Y(I),I=1,NGPT)
        WRITE(6,*) 'Gauss point Z:',(Z(I),I=1,NGPT)
      ENDIF
      IF(ICHX.EQ.1) THEN
*       Variational collocation method
        CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
        MKN=MAXKN/(NXD*NYD*NZD)
        ALLOCATE(KN(MAXKN))
        CALL LCMGET(IPTRK,'KN',KN)
        CALL LCMSIX(IPTRK,'BIVCOL',1)
        CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
        CALL LCMGET(IPTRK,'E',E)
        CALL LCMSIX(IPTRK,' ',2)
        CALL VALUE2(LC,MKN,NXD,NYD,NZD,L4,X,Y,Z,MXD,MYD,MZD,
     1  FLXD,MAT,KN,NGPT,NGPT,NGPT,E,FLUGP)
        DEALLOCATE(KN)
      ELSE IF(ICHX.EQ.2) THEN
*       Raviart-Thomas finite element method
        CALL VALUE4(IEL,NUN,NXD,NYD,NZD,X,Y,Z,MXD,MYD,MZD,
     1  FLXD,MAT,KEYFLX,NGPT,NGPT,NGPT,FLUGP)
      ELSE IF(ICHX.EQ.3) THEN
*       Nodal collocation method (MCFD)
        CALL VALUE1(IDIM,NXD,NYD,NZD,L4,X,Y,Z,MXD,MYD,MZD,
     1  FLXD,MAT,IEL,NGPT,NGPT,NGPT,FLUGP)
      ELSE
        CALL XABORT('NAPPPR: INTERPOLATION NOT IMPLEMENTED.')
      ENDIF
      IF(IMPX.GE.10) then
        WRITE(6,*) 'Gauss flux values:'
        DO KGP=1,NGPT
          WRITE(6,*) 'KGP=:',KGP
        DO JGP=1,NGPT
          WRITE(6,*) (FLUGP(IGP,JGP,KGP),IGP=1,NGPT)
        ENDDO
        ENDDO
      ENDIF
*     integrate flux (gauss method)
      DO IGP=1,NGPT
      DO JGP=1,NGPT
      DO KGP=1,NGPT
        FPD=FPD+FLUGP(IGP,JGP,KGP)*WGKSIX(IGP)*WGKSIY(JGP)*WGKSIZ(KGP)
      ENDDO
      ENDDO
      ENDDO
*        GEt averaGE flux
      FPD=FPD/DX/DY/DZ
      if(LDEBUG)write(6,*)'FLXP,FPD,FXTD,FYTD',FLXP(IP,JP,K,IG,IASS),
     1                     FPD,FXTD(IP,ID),FYTD(JP,JD)
   
      FLXP(IP,JP,K,IG,IASS)=FLXP(IP,JP,K,IG,IASS)
     1                     +FPD*FXTD(IP,ID)*FYTD(JP,JD)
      if(LDEBUG)write(6,*)'FLXP after',FLXP(IP,JP,K,IG,IASS)
*     indent back
          ENDIF
        ENDDO
        ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
*
      DEALLOCATE(FXTD,FYTD)
      DEALLOCATE(MXDA,MYDA,MZDA)
      DEALLOCATE(FLXDA)
      IF(IMPX.GE.100)WRITE(6,*) 'debug PPR:projection flux for one '
     1  //'assem end'
!     end of DO IASS=1,NASS
      ENDDO
      IF(IMPX.GE.100)WRITE(6,*) 'debug PPR:projection flux for all '
     1  //'assem end'
* GPPR
      IF(IMETH.EQ.1) THEN
        IF(IMPX.GE.100)WRITE(6,*) 'debug PPR:',TFDINF
*   GEt Volume, phi^t,inf_p and phi^d,inf_m,p from macrolib of fuel
*      Note: if homoGEneous (normal PPR), m=1
        ALLOCATE(VOLM(NGFF),HFM(NMIX,NGFF,NG),
     1   FTINFM(NMIX,NGFF,NG),FDINFM(NMIX,NGFF,NG))
        ALLOCATE(VOL(NPIN,NPIN,NZASS,NASS))
        ALLOCATE(HF(NPIN,NPIN,NZASS,NG,NASS))
        ALLOCATE(FTINF(NPIN,NPIN,NZASS,NG,NASS))
        ALLOCATE(FDINF(NPIN,NPIN,NZASS,NG,NASS))
        ALLOCATE(BMIXP(NPIN*NPIN))
        CALL XDRSET(VOL,NPIN*NPIN*NZASS*NASS,0.0)
        CALL XDRSET(HF,NPIN*NPIN*NZASS*NG*NASS,0.0)
        CALL XDRSET(FTINF,NPIN*NPIN*NZASS*NG*NASS,0.0)
        CALL XDRSET(FDINF,NPIN*NPIN*NZASS*NG*NASS,0.0)

        if(LDEBUG) call LCMLIB(IPMAC)
        CALL LCMSIX(IPMAC,'GFF',1)
        if(LDEBUG) call LCMLIB(IPMAC)
        CALL LCMGET(IPMAC,'VOLUME',VOLM)
        CALL LCMGET(IPMAC,'H-FACTOR',HFM)
        CALL LCMGET(IPMAC,'NWT0',FTINFM)
        CALL LCMGET(IPMAC,TFDINF,FDINFM)
        CALL LCMSIX(IPMAC,'GFF-GEOM',1)
        CALL LCMGET(IPMAC,'MIX',BMIXP)
        CALL LCMSIX(IPMAC,'GFF-GEOM',2)
        CALL LCMSIX(IPMAC,'GFF',2)

        DO IG=1,NG
        
        DO IASS=1,NASS
          K1=CODEA(IASS,5)
          DO K=1,NZASS
!  NAMIX = 1 for homogeneous assembly 
!        > 1 for heterogeneous assembly
!             Note that all values of HFM are identical 
!             for all the mix in a specific assembly
          IMIX=(IASS-1+(K-1)*NASS)*NAMIX+1
          DO J=1,NPIN
            DO I=1,NPIN
            IPIN=I+(J-1)*NPIN
            HF(I,J,K,IG,IASS)=HFM(IMIX,BMIXP(IPIN),IG)
            FTINF(I,J,K,IG,IASS)=FTINFM(IMIX,BMIXP(IPIN),IG)
            FDINF(I,J,K,IG,IASS)=FDINFM(IMIX,BMIXP(IPIN),IG)
            VOL(I,J,K,IASS)=(MXP(I+1)-MXP(I))*(MYP(J+1)-MYP(J))
     3         *(MZM(K1+K)-MZM(K1+K-1))
            ENDDO
          ENDDO
          ENDDO
!       end of DO IASS=1,NASS
        ENDDO
!       end of DO IG=1,NG
        ENDDO
        IF(IMPX.GE.6) then
        DO iass=1,nass
          WRITE(6,*) 'XS for assembly #',IASS
          DO k=1,nzass
           WRITE(6,*) 'Plane #',K
          DO ig=1,ng
           WRITE(6,*) 'group #',ig
           WRITE(6,*) 'HF #'
          DO j=1,npin
           WRITE(6,*) (HF(I,J,K,ig,iass),I=1,NPIN)
          ENDDO
          WRITE(6,*) 'FTINF #'
          DO j=1,npin
           WRITE(6,*) (FTINF(I,J,K,ig,iass),I=1,NPIN)
          ENDDO
          WRITE(6,*) 'FLXP #'
          DO j=1,npin
           WRITE(6,*) (FLXP(I,J,K,ig,iass),I=1,NPIN)
          ENDDO
          WRITE(6,*) 'FDINF #'
          DO j=1,npin
           WRITE(6,*) (FDINF(I,J,K,ig,iass),I=1,NPIN)
          ENDDO
!           end of do ig=1,ng
          ENDDO
          WRITE(6,*) 'VOL #'
          DO j=1,npin
           WRITE(6,*) (VOL(I,J,K,iass),I=1,NPIN)
          ENDDO
          ENDDO
        ENDDO
        ENDIF
* Print and save reaction rates
        IF(IMPX.GE.100)WRITE(6,*) 'debug: Print and save reaction rates'
        POWTOT=0.0
        DO IASS=1,NASS
          PWASS2(IASS)=0.0
          K1=CODEA(IASS,5)
          DO K=1,NZASS
            DO J=1,NPIN
            DO I=1,NPIN
              VTOT=VTOT+VOL(I,J,K,IASS)
            DO IG=1,NG
              POWTOT=POWTOT+HF(I,J,K,IG,IASS)*FTINF(I,J,K,IG,IASS)
     1         *FLXP(I,J,K,IG,IASS)/FDINF(I,J,K,IG,IASS)
     2         *VOL(I,J,K,IASS)
              PWASS2(IASS)=PWASS2(IASS)
     1         +HF(I,J,K,IG,IASS)*FTINF(I,J,K,IG,IASS)
     1         *FLXP(I,J,K,IG,IASS)/FDINF(I,J,K,IG,IASS)
     2         *VOL(I,J,K,IASS)
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDDO
        IF(IMPX.GE.2) WRITE(6,*) 'POWTOT:',POWTOT
        IF(LPOW) THEN
          DO IASS=1,NASS
            FACTASS(IASS)=POW/POWTOT
          ENDDO
        ELSE
          DO IASS=1,NASS
            FACTASS(IASS)=PWASS(IASS)/PWASS2(IASS)
          ENDDO
        ENDIF
        IF(IMPX.GE.2) WRITE(6,*) 'FACTASS:',(FACTASS(I),I=1,NASS)
        ALLOCATE(HFA(NPIN,NPIN,NZASS))
        ALLOCATE(FXYZ(NZASS))
        ALLOCATE(PLINMAXZ(NZASS))
        DO K=1,NZASS
          FXYZ(K)=0.0
          PLINMAXZ(K)=0.0
        ENDDO
        JPMAP=LCMLID(IPMAP,'ASSEMBLY',NASS)
        IAX=0
        JAX=1
        DO IASS=1,NASS
        K1=CODEA(IASS,5)
        IAX=IAX+1
        IF(IAX.GT.NBAX(JAX)) THEN
          IAX=1
          JAX=JAX+1
        ENDIF
        WRITE(LABEL,'(A4,A4)') NAMY(JAX),NAMX(IBAX(JAX)+IAX-1)
        IF(IMPX.GE.5) THEN
          WRITE(6,*) 'Reaction rates for assembly #',IASS,' Label:',
     1     LABEL
        ENDIF
        DO K=1,NZASS
        IF(IMPX.GE.5) WRITE(6,*) 'Plane #',K
        DO J=1,NPIN
        DO I=1,NPIN
        HFA(I,J,K)=0.0
        DO IG=1,NG
          HFA(I,J,K)=HFA(I,J,K)+HF(I,J,K,IG,IASS)*FTINF(I,J,K,IG,IASS)
     1        *FLXP(I,J,K,IG,IASS)/FDINF(I,J,K,IG,IASS)
     2        *FACTASS(IASS)
     2        *VOL(I,J,K,IASS)
          IF((PLINMAXZ(K)*(MZM(K1+K)-MZM(K1+K-1))).LT.HFA(I,J,K)) THEN
            PLINMAXZ(K)=HFA(I,J,K)/(MZM(K1+K)-MZM(K1+K-1))
          ENDIF
!       end of DO IG=1,NG
        ENDDO
!       end of I=1,NPIN
        ENDDO
          IF(IMPX.GE.5) WRITE(6,*) (HFA(I,J,K),I=1,NPIN)
!       end of J=1,NPIN
        ENDDO
!       end of DO K=1,NZASS
        ENDDO
*
        KPMAP=LCMDIL(JPMAP,IASS)
        CALL LCMPTC(KPMAP,'LABEL',8,1,LABEL)
        CALL LCMPUT(KPMAP,'PIN-POWER',NPIN*NPIN*NZASS,2,HFA)
        CALL LCMPUT(KPMAP,'FLUX',NPIN*NPIN*NZASS*NG,2,
     1              FLXP(1,1,1,1,IASS))
        POWASS=0.0
        DO I=1,NPIN
        DO J=1,NPIN
        DO K=1,NZASS
          POWASS=POWASS+HFA(I,J,K)!power of the assembly iass
          VPIN(I,J,IASS)=VPIN(I,J,IASS)+VOL(I,J,K,IASS)
          !pin volume
        ENDDO
        ENDDO
        ENDDO
        DO I=1,NPIN
        DO J=1,NPIN
        DO K=1,NZASS
          AXPOW(I,J,IASS)=HFA(I,J,K)
     2    +AXPOW(I,J,IASS)
          !AXPOW:axially integrated pin power per pin
          !normalized to the pin mean power
          IF(PMAX.LT.HFA(I,J,K)) THEN
            PMAX=HFA(I,J,K)!maximal 3D local power
          ENDIF
        ENDDO
        AXPOW(I,J,IASS)=AXPOW(I,J,IASS)/(POWASS/NPIN/NPIN)
        ENDDO
        ENDDO
*
        IF(IMPX.GE.2) WRITE(6,*) 'Power of assembly #',IASS,":",POWASS
        DO I=1,NPIN
          DO J=1,NPIN
            IF(IMPX.GE.2) THEN
              WRITE(6,*) 'AXPOW for assembly #',IASS
              NBPIN=NBPIN+1
              WRITE(6,*) 'ASS:',IASS,'PIN #',NBPIN,":",AXPOW(I,J,IASS)
            ENDIF
            PINPOW=AXPOW(I,J,IASS)*VPIN(I,J,IASS)
            IF(HOTPINPOW.LT.PINPOW) THEN
              HOTPINPOW=PINPOW
        !power of the hot pin normalized to the pin mean power
            ENDIF
            IF(FXYASS(IASS).LT.AXPOW(I,J,IASS)) THEN
              FXYASS(IASS)=AXPOW(I,J,IASS)
            ENDIF
          ENDDO
        ENDDO
        NBPIN=0
*
        IF(IMPX.GE.1) THEN
          WRITE(6,*) 'Fxy for assembly #',IASS,":",FXYASS(IASS)
        ENDIF
        CALL LCMPUT(KPMAP,'ASS-POWER',1,2,POWASS)
!       end of DO IASS=1,NASS
        ENDDO
!     end of IF(IMETH.EQ.1) THEN
      ENDIF
*
      FQ=PMAX
      FXY=HOTPINPOW
*
      IF(IMPX.GE.0) THEN
        WRITE(6,*) 'FQ=',FQ
        WRITE(6,*) 'FXY=',FXY
      DO K=1,NZASS
        FXYZ(K)=PLINMAXZ(K)
        IF(IMPX.GE.0) WRITE(6,*) 'Plane #',K,'---> FXYZ(Z)=',FXYZ(K)
      ENDDO
      ENDIF
      CALL LCMPUT(IPMAP,'FQ',1,2,FQ)
      CALL LCMPUT(IPMAP,'FXY',1,2,FXY)
      CALL LCMPUT(IPMAP,'FXYZ',NZASS,2,FXYZ)
      CALL LCMPUT(IPMAP,'FXYASS',IASS,2,FXYASS)
      CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
      ISTATE(17)=NZASS
      CALL LCMPUT(IPMAP,'STATE-VECTOR',NSTATE,1,ISTATE)
!
      IF(IMPX.GE.100)WRITE(6,*) 'debug: beging deallacate'

      DEALLOCATE(FLXP,FLXD)
      DEALLOCATE(MXP,MYP)

      DEALLOCATE(CODEA)
      IF(IMETH.EQ.1) THEN 
        DEALLOCATE(VOLM,HFM,FTINFM,FDINFM)
        DEALLOCATE(VOL,HF,FTINF,FDINF,AXPOW,FXYASS)
        DEALLOCATE(HFA,FXYZ,PLINMAXZ,VPIN)
        DEALLOCATE(FACTASS,PWASS,PWASS2)
        DEALLOCATE(BUNDPW)
      ENDIF
      DEALLOCATE(ACOMB)
      DEALLOCATE(AZONE)
      DEALLOCATE(NAMX,NAMY)
      DEALLOCATE(BMIX,BMIXP)
      DEALLOCATE(MXM,MYM,MZM)
      DEALLOCATE(NBAX,IBAX)
      DEALLOCATE(KEYFLX,MAT)
      DEALLOCATE(MXD,MYD,MZD)

      RETURN
  500 FORMAT(5HFINF_,I3.3,4H    )
      END
