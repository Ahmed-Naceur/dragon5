*DECK NAPCPO
      SUBROUTINE NAPCPO(IPCPO,IPTRK,IPFLU,NSTATE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Construct an 'enriched' multicompo with additional information  
* needed by Pin Power Reconstruction.
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal.
*
*Author(s): 
* R. Chambon
*
*Parameters: input/output
* IPCPO   LCM object address of Multicompo.
* IPTRK   LCM object address of Tracking.
* IPFLU   LCM object address of Flux.
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
      TYPE(C_PTR)  IPCPO,IPTRK,IPFLU
*----
*  LOCAL VARIABLES
*----
      INTEGER NCRCAL,NGPT
      INTEGER IOUT,MAXPAR,MAXLIN,MAXVAL,MAXADD,MAXIFX
      REAL REPS
      PARAMETER (REPS=1.0E-4,IOUT=6,MAXPAR=50,MAXLIN=50,MAXVAL=200,
     1 MAXADD=10,MAXIFX=5,NGPT=2)
      CHARACTER PARKEY(MAXPAR)*12,PARFMT(MAXPAR)*8,RECNAM*12,
     1  COMMEN(MAXLIN)*80,PARKEL(MAXPAR)*12,VALH(MAXPAR)*12
      TYPE(C_PTR) JPCPO,KPCPO,LPCPO,JPMIC,JPFLU
      CHARACTER TEXT*12,HSMG*131,DIRHOM*12,VCHAR(MAXVAL)*12,HVECT*8
      INTEGER ISTATE(NSTATE),NPAR,NLOC,IMPX,IEL,NFDI,FINF(MAXIFX)
      INTEGER IPAR,IBMOLD,IFX
      REAL FLOT
      DOUBLE PRECISION DFLOT
      INTEGER VALI(MAXPAR),NVALUE(MAXPAR),VINTE(MAXVAL),
     1 MUPLET(2*MAXPAR)
      REAL VALR(2*MAXPAR,2),VREAL(MAXVAL),NVPO(2),PTR,PDF,PDF2
      REAL ZGKSIX(NGPT),ZGKSIY(NGPT),WGKSIX(NGPT), WGKSIY(NGPT),
     1 FLUGP(NGPT,NGPT),FPD,DX,DY
      INTEGER NGFF,NXP,NYP,ITYPGP,NMIXP,NMIL,NG,NCOMLI,MAXNVP,STYPP,
     1 NMCAL
      INTEGER I,J,ICAL,INDIC,ITYLCM,IX,IY,LENGTH,NITMA,IREG,IREGP,IG,
     1 IMIXP,ID,JD,IGP,JGP,IP,JP,J1,ICHX,IDIM,LC,L4,MAXKN,MKN
      INTEGER NREG,NUN,NXD,NYD,ITYPGD,NREGP
      REAL E(25)
      LOGICAL LNOINT,FLAG
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIXP
      REAL, ALLOCATABLE, DIMENSION(:) :: MXP,MYP,KN
      REAL, ALLOCATABLE, DIMENSION(:) :: MXD,MYD,XX,YY
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEYFLX,MATCOD,IXPD,
     1  IYPD,JDEBAR,JARBVA
      REAL, ALLOCATABLE, DIMENSION(:) :: FLXD
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FLXP,FT,FXTD,FYTD
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: FLAGMX
*
      IMPX=0
      CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
      IF(INDIC.NE.3) CALL XABORT('NAPCPO: character data expected.')
      IF(TEXT.EQ.'EDIT') THEN
        CALL REDGET(INDIC,IMPX,FLOT,TEXT,DFLOT)
        IF(INDIC.NE.1) CALL XABORT('NAPCPO: integer data expected.')
        CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
        IF(INDIC.NE.3) CALL XABORT('NAPCPO: character data expected.')
      ENDIF
      IF(TEXT.NE.'PROJECTION') CALL XABORT('NAPCPO: ''PROJECTION'' '//
     1  'EXPECTED.')
      CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
      IF(TEXT.NE.'STEP') CALL XABORT('NAPCPO: ''STEP'' '//
     1  'EXPECTED.')
      CALL REDGET(INDIC,NITMA,FLOT,DIRHOM,DFLOT)
      IF(INDIC.NE.3) CALL XABORT('NAPCPO: character data expected.')
      CALL LCMSIX(IPCPO,DIRHOM,1)
*      
      IFX=1
      LNOINT=.FALSE.
      CALL XDISET(FINF,MAXIFX,-1)
*----
*  RECOVER TABLE-OF-CONTENT INFORMATION FOR THE COMPO.
*----
      CALL LCMGET(IPCPO,'STATE-VECTOR',ISTATE)
      NMIL=ISTATE(1)
      NG=ISTATE(2)
      NMCAL=ISTATE(4)
      NPAR=ISTATE(5)
      NLOC=ISTATE(6)
      NCOMLI=ISTATE(10)
      NGFF=ISTATE(14)
      IF(NGFF.EQ.0) CALL XABORT('NAPCPO: NO GFF INFO IN MULTICOMPO.')
      CALL LCMGTC(IPCPO,'COMMENT',80,NCOMLI,COMMEN)
      IF(NPAR.GT.0)THEN
        CALL LCMSIX(IPCPO,'GLOBAL',1)
        CALL LCMGTC(IPCPO,'PARKEY',12,NPAR,PARKEY)
        CALL LCMGTC(IPCPO,'PARFMT',8,NPAR,PARFMT)
        CALL LCMGET(IPCPO,'NVALUE',NVALUE)
        IF(IMPX.GT.10)THEN
          DO IPAR=1,NPAR
            WRITE(RECNAM,'(''pval'',I8.8)') IPAR
            IF(PARFMT(IPAR).EQ.'INTEGER') THEN
              CALL LCMGET(IPCPO,RECNAM,VINTE)
              WRITE(IOUT,'(13H NAPCPO: KEY=,A,18H TABULATED POINTS=,
     1        1P,6I12/(43X,6I12))') PARKEY(IPAR),(VINTE(I),I=1,
     2        NVALUE(IPAR))
            ELSE IF(PARFMT(IPAR).EQ.'REAL') THEN
              CALL LCMGET(IPCPO,RECNAM,VREAL)
              WRITE(IOUT,'(13H NAPCPO: KEY=,A,18H TABULATED POINTS=,
     1        1P,6E12.4/(43X,6E12.4))') PARKEY(IPAR),(VREAL(I),I=1,
     2        NVALUE(IPAR))
            ELSE IF(PARFMT(IPAR).EQ.'STRING') THEN
              CALL LCMGTC(IPCPO,RECNAM,12,NVALUE(IPAR),VCHAR)
              WRITE(IOUT,'(13H NAPCPO: KEY=,A,18H TABULATED POINTS=,
     1        1P,6A12/(43X,6A12))') PARKEY(IPAR),(VCHAR(I),I=1,
     2        NVALUE(IPAR))
            ENDIF
          ENDDO
        ENDIF
        CALL LCMSIX(IPCPO,' ',2)
      ENDIF
      IF(NLOC.GT.0)THEN
        CALL LCMSIX(IPCPO,'LOCAL',1)
        CALL LCMGTC(IPCPO,'PARKEY',12,NLOC,PARKEL)
        CALL LCMSIX(IPCPO,' ',2)
        JPCPO=LCMGID(IPCPO,'MIXTURES')
        DO IBMOLD=1,NMIL
          KPCPO=LCMGIL(JPCPO,IBMOLD)
          LPCPO=LCMGID(KPCPO,'TREE')
          CALL LCMGET(LPCPO,'NVALUE',NVALUE)
          IF(IMPX.GT.10)THEN
            WRITE(IOUT,'(17H NAPCPO: MIXTURE=,I6)') IBMOLD
            DO IPAR=1,NLOC
              WRITE(RECNAM,'(''pval'',I8.8)') IPAR
              CALL LCMGET(LPCPO,RECNAM,VREAL)
              WRITE(IOUT,'(13H NAPCPO: KEY=,A,18H TABULATED POINTS=,
     1        1P,6E12.4/(43X,6E12.4))') PARKEL(IPAR),(VREAL(I),I=1,
     2        NVALUE(IPAR))
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      IF(IMPX.GT.10)WRITE(IOUT,'(1X,A)') (COMMEN(I),I=1,NCOMLI)
*----
*  READ (INTERP_DATA) AND SET VALI, VALR AND VALH PARAMETERS
*  CORRESPONDING TO THE INTERPOLATION POINT. FILL MUPLET FOR
*  PARAMETERS.
*----
      CALL XDISET(MUPLET,NPAR+NLOC,0)
 1020 CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
      IF(INDIC.NE.3) CALL XABORT('NAPCPO: character data expected.')
      IF(TEXT.EQ.'SET') THEN
         CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
         IF(INDIC.NE.3) CALL XABORT('NAPCPO: character data expected.')
         DO 50 I=1,NPAR
         IF(TEXT.EQ.PARKEY(I)) THEN
            IPAR=I
            GO TO 60
         ENDIF
   50    CONTINUE
         GO TO 100
   60    LPCPO=LCMGID(IPCPO,'GLOBAL')
         CALL LCMGET(LPCPO,'NVALUE',NVALUE)
         IF(NVALUE(IPAR).GT.MAXVAL) CALL XABORT('NAPCPO: MAXVAL OVERFL'
     1   //'OW.')
         WRITE(RECNAM,'(''pval'',I8.8)') IPAR
         CALL LCMLEN(LPCPO,RECNAM,LENGTH,ITYLCM)
         IF(LENGTH.EQ.0) THEN
            WRITE(HSMG,'(25HNAPCPO: GLOBAL PARAMETER ,A,9H NOT SET.)')
     1      PARKEY(IPAR)
            CALL XABORT(HSMG)
         ENDIF
         IF(PARFMT(IPAR).EQ.'INTEGER') THEN
            CALL REDGET(INDIC,VALI(IPAR),FLOT,TEXT,DFLOT)
            IF(INDIC.NE.1) CALL XABORT('NAPCPO: integer data expected.')
            CALL LCMGET(LPCPO,RECNAM,VINTE)
            DO J=1,NVALUE(IPAR)
              IF(VALI(IPAR).EQ.VINTE(J)) THEN
                MUPLET(IPAR)=J
*               MUTYPE(IPAR)=ITYPGD
                GO TO 1020
              ENDIF
            ENDDO
            WRITE(HSMG,'(26HNAPCPO: INTEGER PARAMETER ,A,9H WITH VAL,
     1      2HUE,I5,29H NOT FOUND IN COMPO DATABASE.)') PARKEY(IPAR),
     2      VALI(IPAR)
            CALL XABORT(HSMG)
         ELSE IF(PARFMT(IPAR).EQ.'REAL') THEN
            CALL REDGET(INDIC,NITMA,VALR(IPAR,1),TEXT,DFLOT)
            IF(INDIC.NE.2) CALL XABORT('NAPCPO: real data expected.')
            VALR(IPAR,2)=VALR(IPAR,1)
!            CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
            CALL LCMGET(LPCPO,RECNAM,VREAL)
            DO J=1,NVALUE(IPAR)
              IF(ABS(VALR(IPAR,1)-VREAL(J)).LE.REPS*ABS(VREAL(J))) THEN
                MUPLET(IPAR)=J
*               MUTYPE(IPAR)=ITYPGD
               GO TO 1020
              ENDIF
            ENDDO
            WRITE(HSMG,'(23HNAPCPO: REAL PARAMETER ,A,9H WITH VAL,
     1      2HUE,I5,29H NOT FOUND IN COMPO DATABASE.)') PARKEY(IPAR),
     2      VALR(IPAR,1)
            CALL XABORT(HSMG)
         ELSE IF(PARFMT(IPAR).EQ.'STRING') THEN
            CALL REDGET(INDIC,NITMA,FLOT,VALH(IPAR),DFLOT)
            IF(INDIC.NE.3) CALL XABORT('NAPCPO: STRING DATA EXPECTED.')
            CALL LCMGTC(LPCPO,RECNAM,12,NVALUE(IPAR),VCHAR)
            DO J=1,NVALUE(IPAR)
              IF(VALH(IPAR).EQ.VCHAR(J)) THEN
                MUPLET(IPAR)=J
*               MUTYPE(IPAR)=ITYPGD
                GO TO 1020
              ENDIF
            ENDDO
            WRITE(HSMG,'(25HNAPCPO: STRING PARAMETER ,A,10H WITH VALU,
     1      2HE ,A12,29H NOT FOUND IN COMPO DATABASE.)') PARKEY(IPAR),
     2      VALH(IPAR)
            CALL XABORT(HSMG)
         ENDIF
  100    DO 110 I=1,NLOC
         IF(TEXT.EQ.PARKEL(I)) THEN
            IPAR=NPAR+I
            GO TO 120
         ENDIF
  110    CONTINUE
         CALL XABORT('NAPCPO: PARAMETER '//TEXT//' NOT FOUND.')
  120    JPCPO=LCMGID(IPCPO,'MIXTURES')
         IBMOLD=1
         KPCPO=LCMGIL(JPCPO,IBMOLD)
         LPCPO=LCMGID(KPCPO,'TREE')
         CALL LCMGET(LPCPO,'NVALUE',NVALUE)
         CALL REDGET(INDIC,NITMA,VALR(IPAR,1),TEXT,DFLOT)
         IF(INDIC.NE.2) CALL XABORT('NAPCPO: real data expected.')
         VALR(IPAR,2)=VALR(IPAR,1)
         WRITE(RECNAM,'(''pval'',I8.8)') IPAR-NPAR
         CALL LCMLEN(LPCPO,RECNAM,LENGTH,ITYLCM)
         IF(LENGTH.EQ.0) THEN
            WRITE(HSMG,'(24HNAPCPO: LOCAL PARAMETER ,A,9H NOT SET.)')
     1      PARKEL(IPAR-NPAR)
            CALL XABORT(HSMG)
         ELSE IF(LENGTH.GT.MAXVAL) THEN
            CALL XABORT('NAPCPO: MAXVAL OVERFLOW.')
         ENDIF
         CALL LCMGET(LPCPO,RECNAM,VREAL)
         DO J=1,NVALUE(IPAR-NPAR)
           IF(ABS(VALR(IPAR,1)-VREAL(J)).LE.REPS*ABS(VREAL(J))) THEN
             MUPLET(IPAR)=J
*            MUTYPE(IPAR)=ITYPGD
             GO TO 1020
           ENDIF
         ENDDO
         WRITE(HSMG,'(26HNAPCPO: INTEGER PARAMETER ,A,9H WITH VAL,
     1   2HUE,I5,29H NOT FOUND IN COMPO DATABASE.)') PARKEY(IPAR),
     2   VALI(IPAR)
         CALL XABORT(HSMG)
      ELSEIF(TEXT.EQ.'IFX') THEN
        CALL REDGET(INDIC,IFX,FLOT,TEXT,DFLOT)
        IF(INDIC.NE.1) CALL XABORT('NAPCPO: integer data expected.')
        GO TO 1020
      ELSEIF(TEXT.EQ.'NOINTP') THEN
        LNOINT=.TRUE.
        GO TO 1020
      ELSEIF(TEXT.EQ.'INTERP') THEN
        LNOINT=.FALSE.
        GO TO 1020
      ELSEIF(TEXT.EQ.';') THEN
        GOTO 200
      ENDIF
      CALL XABORT('NAPCPO: '//TEXT//' is a wrong keyword')
*
  200 CONTINUE
      JPCPO=LCMGID(IPCPO,'MIXTURES')
      IBMOLD=1
      KPCPO=LCMGIL(JPCPO,IBMOLD)
      LPCPO=LCMGID(KPCPO,'TREE')
      CALL LCMGET(LPCPO,'NVP',NVPO)
      CALL LCMLEN(LPCPO,'ARBVAL',MAXNVP,ITYLCM)
      IF(NVPO(1).GT.MAXNVP) CALL XABORT('NAPCPO: NVP OVERFLOW.')
      ALLOCATE(JDEBAR(MAXNVP+1),JARBVA(MAXNVP))
      CALL LCMGET(LPCPO,'DEBARB',JDEBAR)
      CALL LCMGET(LPCPO,'ARBVAL',JARBVA)
      IF(IMPX.GE.20) THEN
        WRITE(6,*) 'MUPLET: ',(MUPLET(I),I=1,NPAR+NLOC)
      ENDIF
      ICAL=NCRCAL(1,NVPO(1),NPAR+NLOC,JDEBAR,JARBVA,MUPLET)
      IF(IMPX.GE.2) THEN
        WRITE(6,*) 'Performing projection for calculation: ',ICAL
      ENDIF
*
      LPCPO=LCMGID(KPCPO,'CALCULATIONS')
      JPMIC=LCMGIL(LPCPO,ICAL)
      CALL LCMGET(JPMIC,'STATE-VECTOR',ISTATE)
      CALL LCMSIX(JPMIC,'MACROLIB    ',1)
      CALL LCMSIX(JPMIC,'GFF         ',1)
      CALL LCMSIX(JPMIC,'GFF-GEOM    ',1)
C get dimension in geometry from L_MULTICOMPO
      CALL LCMGET(JPMIC,'STATE-VECTOR',ISTATE)
      ITYPGP=ISTATE(1)
      STYPP=ISTATE(11)
      IF(ITYPGP.NE.5) CALL XABORT('NAPCPO: CAR2D geometry type '
     1 //'expected in L_MULTICOMPO.')
      IF(STYPP.NE.0) CALL XABORT('NAPCPO: No split in geometry expected'
     1 //' in L_MULTICOMPO.')
      NXP=ISTATE(3)
      NYP=ISTATE(4)
      NREGP=ISTATE(6)
      NMIXP=ISTATE(7)
      IF(NMIXP.NE.NGFF) CALL XABORT('NAPCPO: INVALID GFF-GEOM.')
      ALLOCATE(MXP(NXP+1),MYP(NYP+1))
      ALLOCATE(IXPD(NXP+1),IYPD(NYP+1))
      ALLOCATE(MIXP(NREGP))
      CALL LCMGET(JPMIC,'MESHX',MXP)
      CALL LCMGET(JPMIC,'MESHY',MYP)
      CALL LCMGET(JPMIC,'MIX',MIXP)
      CALL LCMSIX(JPMIC,'GFF-GEOM    ',2)
      CALL XDISET(IXPD,NXP,0)
      CALL XDISET(IYPD,NYP,0)
C get dimension in geometry from L_TRACK
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NREG=ISTATE(1)
      NUN=ISTATE(2)
      ITYPGD=ISTATE(6)
      IEL=ISTATE(9)
      IF(ITYPGD.NE.5) CALL XABORT('NAPCPO: CAR2D geometry type expected'
     1 //' in L_TRACKING.')
      IEL=ISTATE(9)
      L4=ISTATE(11)
      ICHX=ISTATE(12)
      NXD=ISTATE(14)
      NYD=ISTATE(15)
      IDIM=2
      IF(NREG.NE.NXD*NYD) CALL XABORT('NAPCPO: No Splitting allowed in '
     1 //'CAR2D geometry type from L_TRACK.')
C     compute X and Y mesh from L_TRACK
      ALLOCATE(MXD(NXD+1),MYD(NYD+1))
      ALLOCATE(XX(NREG),YY(NREG))
      CALL LCMGET(IPTRK,'XX',XX)
      CALL LCMGET(IPTRK,'YY',YY)
      MXD(1)=MXP(1)
      DO I=1,NXD
        MXD(I+1)=MXD(I)+XX(I)
      ENDDO
      MYD(1)=MYP(1)
      DO I=1,NYD
        MYD(I+1)=MYD(I)+YY((I-1)*NXD+1)
      ENDDO
      if(IMPX.ge.10) then
        WRITE(6,*) 'Respective mesh (Diffusion vs. Transport):'
        WRITE(6,*) ' X direction :'
        WRITE(6,*) 'MXD:',(MXD(I),I=1,NXD+1)
        WRITE(6,*) 'MXP:',(MXP(I),I=1,NXP+1)
        WRITE(6,*) ' Y direction :'
        WRITE(6,*) 'MYD:',(MYD(I),I=1,NYD+1)
        WRITE(6,*) 'MYP:',(MYP(I),I=1,NYP+1)
      endif
      IF((ABS(MXD(NXD+1)-MXP(NXP+1)).GE.1E-3).OR.
     1   (ABS(MXD(NXD+1)-MXP(NXP+1)).GE.1E-3)) CALL XABORT('NAPCPO: '
     2 //'Diffusion and transport geometries total size mismach')
      ALLOCATE(FXTD(NXP,NXD),FYTD(NYP,NYD))
      CALL XDRSET(FXTD,NXP*NXD,0.0)
      CALL XDRSET(FYTD,NYP*NYD,0.0)
      CALL NAPFTD(NXP,MXP,NXD,MXD,FXTD)
      CALL NAPFTD(NYP,MYP,NYD,MYD,FYTD)
      IF(LNOINT) THEN
C     verify that both meshes match 
      J1=1
      DO I=2,NXD+1
        FLAG=.TRUE.
        DO J=J1,NXP+1
          IF(MXP(J).LT.MXD(I)) THEN
            IXPD(J)=I-1
          ENDIF
          IF(ABS(MXD(I)-MXP(J)).LE.ABS(1E-5*MXP(J))) THEN
            FLAG=.FALSE.
            IXPD(J)=I
            J1=J+1
          ENDIF
        ENDDO
        IF(FLAG) CALL XABORT('NAPCPO: a X mesh in L_TRACK does not '
     1 //'match the CAR2D geometry imbedded in L_MULTICOMPO.')
      ENDDO
      J1=1
      DO I=2,NYD+1
        FLAG=.TRUE.
        DO J=J1,NYP+1
          IF(MYP(J).LT.MYD(I)) THEN
            IYPD(J)=I-1
          ENDIF
          IF(ABS(MYD(I)-MYP(J)).LE.ABS(1E-5*MYP(J))) THEN
            FLAG=.FALSE.
            IYPD(J)=I
            J1=J+1
          ENDIF
        ENDDO
        IF(FLAG) CALL XABORT('NAPCPO: a Y mesh in L_TRACK does not '
     1  //'match the CAR2D geometry imbedded in L_MULTICOMPO.')
      ENDDO
      ENDIF
C     project flux 
      ALLOCATE(KEYFLX(NREG),MATCOD(NREG))
      CALL LCMGET(IPTRK,'KEYFLX',KEYFLX)
      CALL LCMGET(IPTRK,'MATCOD',MATCOD)
      ALLOCATE(FLXD(NUN),FLXP(NMIXP,NG))
      ALLOCATE(FLAGMX(NMIXP))
      JPFLU=LCMGID(IPFLU,'FLUX')
      DO IG=1,NG
        CALL LCMGDL(JPFLU,IG,FLXD)
        DO IP=1,NXP
          DO JP=1,NYP
            IREGP=IP+(JP-1)*NXP
            IF(IREGP.GT.NREGP) CALL XABORT('NAPCPO: NREGP OVERFLOW(1).')
            IMIXP=MIXP(IREGP)
            FLXP(IMIXP,IG)=0.0
            IF(LNOINT) THEN
*     integrated projected flux FLXP
              IREG=IXPD(IP)+(IYPD(JP)-1)*NXD
              FLXP(IMIXP,IG)=FLXD(KEYFLX(IREG))
            ELSE
*     interpolated projected flux FLXP
              DO ID=1,NXD
                DO JD=1,NYD
                  IF(FXTD(IP,ID)*FYTD(JP,JD).NE.0.0) THEN
* -----
      CALL ALGPT(NGPT,MAX(MXP(IP),MXD(ID)),MIN(MXP(IP+1),MXD(ID+1)),
     1   ZGKSIX,WGKSIX)
      DX=MIN(MXP(IP+1),MXD(ID+1))-MAX(MXP(IP),MXD(ID))
      CALL ALGPT(NGPT,MAX(MYP(JP),MYD(JD)),MIN(MYP(JP+1),MYD(JD+1)),
     1   ZGKSIY,WGKSIY)
      DY=MIN(MYP(JP+1),MYD(JD+1))-MAX(MYP(JP),MYD(JD))
      IF(IMPX.GE.5) THEN
        WRITE(6,*) 'IP,JP:',IP,JP,FXTD(IP,ID),'ID,JD:',ID,JD,FYTD(JP,JD)
        WRITE(6,*) 'Gauss point ZGWG:',(ZGKSIX(I),I=1,NGPT),
     1   (WGKSIX(I),I=1,NGPT),'DX',DX
        WRITE(6,*) 'Gauss point ZGWG:',(ZGKSIY(I),I=1,NGPT),
     1   (WGKSIY(I),I=1,NGPT),'DY',DY
      ENDIF
      FPD=0.0
*     interpolate flux
      IF(ICHX.EQ.1) THEN
*       Variational collocation method
        CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
        MKN=MAXKN/(NXD*NYD)
        ALLOCATE(KN(MAXKN))
        CALL LCMGET(IPTRK,'KN',KN)
        CALL LCMSIX(IPTRK,'BIVCOL',1)
        CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
        CALL LCMGET(IPTRK,'E',E)
        CALL LCMSIX(IPTRK,' ',2)
        CALL VALU2B(LC,MKN,NXD,NYD,L4,ZGKSIX,ZGKSIY,MXD,MYD,FLXD,MATCOD,
     1  KN,NGPT,NGPT,E,FLUGP)
      ELSE IF(ICHX.EQ.2) THEN
*       Raviart-Thomas finite element method
        CALL VALU4B(IEL,NUN,NXD,NYD,ZGKSIX,ZGKSIY,MXD,MYD,FLXD,MATCOD,
     1  KEYFLX,NGPT,NGPT,FLUGP)
      ELSE IF(ICHX.EQ.3) THEN
*       Nodal collocation method (MCFD)
        CALL VALU1B(IDIM,NXD,NYD,L4,ZGKSIX,ZGKSIY,MXD,MYD,FLXD,MATCOD,
     1  IEL,NGPT,NGPT,FLUGP)
      ELSE
        CALL XABORT('NAPCPO: INTERPOLATION NOT IMPLEMENTED.')
      ENDIF
      IF(IMPX.GE.5) THEN
        WRITE(6,*) 'Gauss flux values:'
        do JGP=1,NGPT
          WRITE(6,*) (FLUGP(IGP,JGP),IGP=1,NGPT)
        ENDDO
      ENDIF
*     integrate flux (gauss method)
      DO IGP=1,NGPT
      DO JGP=1,NGPT
        FPD=FPD+FLUGP(IGP,JGP)*WGKSIX(IGP)*WGKSIY(JGP)
      ENDDO
      ENDDO
*        get average flux
      FPD=FPD/DX/DY
      FLXP(IMIXP,IG)=FLXP(IMIXP,IG)+FPD*FXTD(IP,ID)*FYTD(JP,JD)
* -----
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C     flux normalization
C         get data from transport calculations
      ALLOCATE(FT(NMIXP,NG))
        CALL LCMGET(JPMIC,'NWT0',FT)
C           group by group
      DO IG=1,NG
C         compute average flux from transport calculations
      PTR=0.0
      IREGP=0
      DO IY=1,NYP
      DO IX=1,NXP
        IREGP=IREGP+1
        IF(IREGP.GT.NREGP) CALL XABORT('NAPCPO: NREGP OVERFLOW(2).')
        IMIXP=MIXP(IREGP)
        PTR=PTR+FT(IMIXP,IG)*(MXP(IX+1)-MXP(IX))*(MYP(IY+1)-MYP(IY))
      ENDDO
      ENDDO
C         compute average flux  with projected diffusion flux
      PDF=0.0
      IREGP=0
      DO IY=1,NYP
      DO IX=1,NXP
        IREGP=IREGP+1
        IF(IREGP.GT.NREGP) CALL XABORT('NAPCPO: NREGP OVERFLOW(3).')
        IMIXP=MIXP(IREGP)
        PDF=PDF+FLXP(IMIXP,IG)*(MXP(IX+1)-MXP(IX))*(MYP(IY+1)-MYP(IY))
      ENDDO
      ENDDO
C         renormalize flux
      DO IMIXP=1,NMIXP
        FLXP(IMIXP,IG)=FLXP(IMIXP,IG)/PDF*PTR
      ENDDO
C
      IF(IMPX.GT.5) THEN
        PDF2=0.0
        IREGP=0
        DO IY=1,NYP
        DO IX=1,NXP
          IREGP=IREGP+1
          IF(IREGP.GT.NREGP) CALL XABORT('NAPCPO: NREGP OVERFLOW(4).')
          IMIXP=MIXP(IREGP)
          PDF2=PDF2+FLXP(IMIXP,IG)
     1          *(MXP(IX+1)-MXP(IX))*(MYP(IY+1)-MYP(IY))
        ENDDO
        ENDDO
        WRITE(6,*)'NAPCPO: transport power:',PTR
        WRITE(6,*)'NAPCPO: diffusion power (before normalization):',PDF
        WRITE(6,*)'NAPCPO: diffusion power (after normalization):',PDF2
        IREGP=0
        WRITE(6,*) 'NAPCPO: FLXP/FT: group #',IG
        DO IY=1,NYP
          WRITE(6,*) (FLXP(MIXP(IREGP+I),IG)
     1                /FT(MIXP(IREGP+I),IG),I=1,NXP)
          IREGP=IREGP+NYP
        ENDDO
      ENDIF
C     verify that all mixtures have a projected flux 
        DO IMIXP=1,NMIXP
          IF(FLXP(IMIXP,IG).EQ.0.0) THEN
            WRITE(HSMG,'(42HNAPCPO: no projected flux for mixture and ,
     1      6Hgroup=,2I6,1H.)') IMIXP,IG
            CALL XABORT(HSMG)
          ENDIF
        ENDDO
C     end DO IG=1,NG
      ENDDO
C     save projected flux in L_MULTICOMPO for each original mixture
      DO I=1,NMIL
        JPCPO=LCMGID(IPCPO,'MIXTURES')
        KPCPO=LCMGIL(JPCPO,I)
        LPCPO=LCMGID(KPCPO,'CALCULATIONS')
        JPMIC=LCMGIL(LPCPO,ICAL)
        CALL LCMSIX(JPMIC,'MACROLIB    ',1)
        CALL LCMSIX(JPMIC,'GFF         ',1)
        CALL LCMLEN(JPMIC,'FINF_NUMBER ',NFDI,ITYLCM)
        IF(NFDI+1.GT.MAXIFX) CALL XABORT('NAPCPO: MAXIFX OVERFLOW.')
        IF(NFDI.GT.0) CALL LCMGET(JPMIC,'FINF_NUMBER ',FINF)
        FINF(NFDI+1)=IFX
        WRITE(HVECT,500) IFX
        CALL LCMPUT(JPMIC,'FINF_NUMBER ',NFDI+1,1,FINF)
        IF(IMPX.GE.10) THEN
          WRITE(6,'(17H NAPCPO: MIXTURE=,I5,8H RECORD ,A8,1H=)') I,HVECT
          DO IG=1,NG
            WRITE(6,'(7H GROUP=,I5/(1X,1P,12E13.4))') IG,
     1      (FLXP(IMIXP,IG),IMIXP=1,NMIXP)
          ENDDO
        ENDIF
        CALL LCMPUT(JPMIC,HVECT,NMIXP*NG,2,FLXP)
        CALL LCMSIX(JPMIC,'GFF         ',2)
        CALL LCMSIX(JPMIC,'*MAC*RES    ',2)
      ENDDO
      DEALLOCATE(FT)
      DEALLOCATE(FLAGMX)
      DEALLOCATE(FLXD,FLXP)
      DEALLOCATE(FXTD,FYTD)
      DEALLOCATE(KEYFLX,MATCOD)
      DEALLOCATE(MXD,MYD)
      DEALLOCATE(XX,YY)
      DEALLOCATE(MIXP)
      DEALLOCATE(MXP,MYP)
      DEALLOCATE(IXPD,IYPD)
      DEALLOCATE(JDEBAR,JARBVA)
      RETURN
*
  500 FORMAT(5HFINF_,I3.3)
      END
