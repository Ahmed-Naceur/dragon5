*DECK EXCELP
      SUBROUTINE EXCELP(  IPTRK, IFTRAK, IPRNTP, NREGIO,  NBMIX,   NANI,
     >                   MATCOD, VOLUME, NRENOR, XSSIGT, XSSIGW, NELPIJ,
     >                    IPIJK,    PIJ, LEAKSW,  NSBG,   NPSYS,   NPST,
     >                   PROBKS, TITREC, NALBP,   ALBP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the collision probabilities for EXCELL. All surfaces
* will disappear from the system using external boundary conditions.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  unit of the sequential binary tracking file.
* IPRNTP  print flag (equal to zero for no print).
* NREGIO  total number of merged blocks for which specific values
*         of the neutron flux and reactions rates are required.
* NBMIX   number of mixtures (NBMIX=max(MATCOD(i))).
* NANI    number of Legendre orders.
* MATCOD  index number of the mixture type assigned to each volume.
* VOLUME  volumes.
* NRENOR  normalization scheme for PIJ matrices.
* XSSIGT  total macroscopic cross sections ordered by mixture.
* XSSIGW  P0 within-group scattering macroscopic cross sections
*         ordered by mixture.
* NELPIJ  number of elements in symmetrized pij matrix.
* IPIJK   pij option (=1 pij, =4 pijk).
* LEAKSW  leakage flag (=.true. if neutron leakage through external
*         boundary is present).
* NSBG    number of energy groups.
* NPSYS   non-converged energy group indices.
* NPST    first dimension of matrix PROBKS.
* TITREC  title.
* NALBP   number of multigroup physical albedos.
* ALBP    multigroup physical albedos.
*
*
*Parameters: output
* PIJ     reduced and symmetrized collision probabilities.
* PROBKS  directional collision probabilities.
*
*-----------------------------------------------------------------------
*--------+---------------- R O U T I N E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                                *
*--------+-------------------------------------------------------------*
* Boundary conditions
*   PIJABC / TO ELIMINATE SURFACES USING B.C. OF THE SYSTEM
*   PIJAAA / TO ELIMINATE SURFACES FOR PIJKS USING B.C. OF THE SYSTEM
* CP INtegration
*   PIJI2D / TO INTEGRATE CP IN 2D GEOMETRIES (ISOTROPIC B.C.)
*   PIJI3D / TO INTEGRATE CP IN 3D GEOMETRIES (ISOTROPIC B.C.)
*   PIJS2D / TO INTEGRATE CP IN 2D GEOMETRIES (SPECULAR  B.C.)
*   PIJS3D / TO INTEGRATE CP IN 3D GEOMETRIES (SPECULAR  B.C.)
* CP Normalisation
*   PIJRDG / TO RENORMALIZE CP USING DIAGONAL COEFFICIENTS
*   PIJRGL / TO RENORMALIZE CP USING GELBARD HOMOGENEOUS SCHEME
*   PIJRNL / TO RENORMALIZE CP USING NON-LINEAR FACTORS
*   PIJRHL / TO RENORMALIZE CP USING HELIOS METHOD
* Various functions
*   PIJWPR / TO PRINT CP MATRICES IN SUM FORMAT
*   PIJSMD / TO EVALUATE SCATTERING-MODIFIED CP MATRIX
*   PIJCMP / COMPRESS CP MATRIX TO SYMETRIC FORMAT
*   PIJD2S / CHARGE PROBKS MATRICES IN THE DRAGON SQUARE FORMAT
*   PIJD2R / CHARGE PIJ MATRICES IN THE DRAGON SYMMETRIZED FORMAT
*   PIJKST / COMPUTE PIJK* MATRICES
* Inline tracking
*   NXTTGC / TRACK CYCLIC NXT LINE IN GEOMETRY
*   NXTTGS / TRACK STANDARD NXT LINE IN GEOMETRY
*   NXTXYZ / READ GEOMETRY LIMITS
*--------+-------------------------------------------------------------*
*
      USE               GANLIB
      IMPLICIT          NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER         TITREC*72
      LOGICAL           LEAKSW
      TYPE(C_PTR)       IPTRK
      INTEGER           IFTRAK, IPRNTP, NREGIO,  NBMIX,   NANI,
     >                  MATCOD(NREGIO),NRENOR, NELPIJ,  IPIJK,   NSBG,
     >                  NPSYS(NSBG), NPST,NALBP
      REAL              VOLUME(NREGIO), XSSIGT(0:NBMIX,NSBG),
     >                  XSSIGW(0:NBMIX,NANI,NSBG),
     >                  PIJ(NELPIJ,IPIJK,NSBG), PROBKS(NPST,NSBG),
     >                  ALBP(NALBP,NSBG)
*----
*  LOCAL VARIABLES
*----
      INTEGER           IOUT, ICPALL, ICPEND, MXGAUS, NSTATE
      PARAMETER       ( IOUT=6, ICPALL=4, ICPEND=3, MXGAUS=64,
     >                  NSTATE=40 )
      CHARACTER         NAMSBR*6
      PARAMETER       ( NAMSBR='EXCELP')
      INTEGER           MKI1, MKI2, MKI3, MKI4, MKI5
      PARAMETER       (MKI1=600,MKI2=600,MKI3=600,MKI4=600,MKI5=600)
      INTEGER           ILONG,ITYPE,ISTATE(NSTATE),ICODE(6)
      INTEGER           NPROB,N2PROB,ISBG,KSBG,ITYPBC
      REAL              ALBEDO(6),EXTKOP(NSTATE),CUTOF,RCUTOF,ASCRP,
     >                  YGSS,XGSS(MXGAUS),WGSS(MXGAUS),WGSSX(MXGAUS),
     >                  FACT,ALBG(6)
      LOGICAL           SWNZBC, SWVOID, LPIJK, LSKIP
      CHARACTER         CTRKT*4, COMENT*80
      DOUBLE PRECISION  WEIGHT,DANG0,DASCRP
*
      INTEGER           JJ,MSYM,IU,IL,ISOUT,IIN,I,J,IBM,IOP,INDPIJ,IJKS,
     >                  NALLOC,ITRAK,IANG,NBSEG,IC,IPRT,ISPEC,IUN,KSPEC,
     >                  LOPT,MXSEG,NALBG,NANGL,NCOMNT,NCOR,NCORT,NDIM,
     >                  NGSS,NNREG,NREG,NSCRP,NSOUT,NTRK,NUNKNO,IVV,
     >                  JGSS,JUN,IFMT,MXSUB,NSUB,ISA
*----
*  Variables for NXT: inline tracking
*----
      INTEGER           ILCMUP,ILCMDN
      PARAMETER        (ILCMUP=1,ILCMDN=2)
      DOUBLE PRECISION  DZERO,DONE,DTWO
      PARAMETER        (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
      INTEGER           IEDIMG(NSTATE),NPOINT,NBUCEL,MXMSH,MAXPIN,
     >                  MXGSUR,MXGREG,MAXMSH,NPLANE,NUCELL(3),IANGL
      CHARACTER         NAMREC*12
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MATALB
      REAL, ALLOCATABLE, DIMENSION(:,:) :: VOLSUR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NUMERO
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: LENGTH,DSV
      REAL, ALLOCATABLE, TARGET, DIMENSION(:,:) :: SIGTAL,SIGT00
      REAL, POINTER, DIMENSION(:,:) :: SIGT
*-- NXT TRACKING
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IUNFLD
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DVNOR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DGMESH,DANGLT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DWGTRK
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DORITR
*-- Temporary arrays
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MATRT
      REAL, ALLOCATABLE, DIMENSION(:) :: LOPATH,FFACT
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SIGANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PSST,PSVT,
     >                                               STAYIN,GOSOUT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DPROB,DPROBX,
     > PCSCT
*----
*  Common blocks for Bickley functions
*----
      INTEGER          L1, L2, L3, L4, L5
      REAL             PAS1,XLIM1,PAS2,XLIM2,PAS3,XLIM3,
     >                 PAS4,XLIM4,PAS5,XLIM5,BI1,BI2,BI3,BI4,BI5
      COMMON /BICKL1/  BI1(0:MKI1,3),PAS1,XLIM1,L1
      COMMON /BICKL2/  BI2(0:MKI2,3),PAS2,XLIM2,L2
      COMMON /BICKL3/  BI3(0:MKI3,3),PAS3,XLIM3,L3
      COMMON /BICKL4/  BI4(0:MKI4,3),PAS4,XLIM4,L4
      COMMON /BICKL5/  BI5(0:MKI5,3),PAS5,XLIM5,L5
      DOUBLE PRECISION ABSC(3,2)
*----
*  INTRINSIC FUNCTION FOR POSITION IN CONDENSE PIJ MATRIX
*----
      INTEGER INDPOS
      INDPOS(I,J)=MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
*----
* RECOVER EXCELL SPECIFIC TRACKING INFORMATION.
*             ALBEDO: SURFACE ALBEDOS (REAL(6))
*             KSPEC : KIND OF PIJ INTEGRATION (0:ISOTROPE,1:SPECULAR)
*             CUTOF : MFP CUTOFF FOR SPECULAR INTEGRATION
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      KSPEC=ISTATE(10)
      CALL LCMGET(IPTRK,'EXCELTRACKOP',EXTKOP)
      CUTOF=EXTKOP(1)
      CALL LCMGET(IPTRK,'ICODE',ICODE)
      CALL LCMGET(IPTRK,'ALBEDO',ALBG)
*
      IPRT  = IPRNTP
      IF( IPRT.GE.ICPEND ) WRITE(IOUT,'(1X,A72//)') TITREC
      NPLANE = 1
      IF(IFTRAK .NE. 0) THEN
        REWIND IFTRAK
        READ(IFTRAK) CTRKT,NCOMNT,NTRK,IFMT
        IF( CTRKT .NE.'$TRK' .OR.
     >      NCOMNT.LT.0      .OR.
     >      NTRK  .EQ.0          ) CALL XABORT(NAMSBR//
     >      ': Invalid tracking file')
        DO IC= 1,NCOMNT
           READ(IFTRAK) COMENT
        ENDDO
        READ(IFTRAK) NDIM,ISPEC,NREG,NSOUT,NALBG,NCOR,NANGL,MXSUB,MXSEG
        IF(NREG.NE.NREGIO )THEN
           CALL XABORT(NAMSBR//': TRACKING FILE HAS INVALID # OF ZONES')
        ENDIF
        NCORT=NCOR
      ELSE
        IF(ISTATE(7) .NE. 4) CALL XABORT(NAMSBR//
     >  ': Tracking file required unless NXT: tracking provided')
        NREG=ISTATE(1)
        NSOUT=ISTATE(5)
        ISPEC=ISTATE(9)
        NPOINT=ISTATE(17)
        MXSEG=ISTATE(18)
        NANGL=ISTATE(20)
        NPLANE=ISTATE(22)
        CALL LCMSIX(IPTRK,'NXTRecords  ',ILCMUP)
        CALL LCMGET(IPTRK,'G00000001DIM',IEDIMG)
        NDIM=IEDIMG(1)
        ITYPBC=IEDIMG( 2)
        NBUCEL=IEDIMG( 5)
        NUCELL(1)=IEDIMG(13)
        NUCELL(2)=IEDIMG(14)
        NUCELL(3)=IEDIMG(15)
        MXMSH=IEDIMG(16)
        MAXPIN=IEDIMG(19)
        MXGSUR=IEDIMG(24)
        MXGREG=IEDIMG(25)
        NCOR=1
        NCORT=NCOR
        NTRK=NANGL*NPLANE*NPOINT**(NDIM-1)
        IF(MXSEG .LE. 1) THEN
          IF(ISPEC .EQ. 0) THEN
            MXSEG=NBUCEL*
     >           ((MAXPIN+1)*(2*MXGREG+2)+MXGSUR+16)
          ELSE
            MXSEG=8*NANGL*NBUCEL*
     >            ((MAXPIN+1)*(2*MXGREG+2)+MXGSUR+16)
          ENDIF
        ENDIF
        MAXMSH=MAX(MXMSH,IEDIMG(17),IEDIMG(20))
      ENDIF
      NNREG = NREGIO*NREGIO
      NUNKNO= NREG+NSOUT+1
      ALLOCATE(MATALB(-NSOUT:NREG),VOLSUR(-NSOUT:NREG,NSBG))
      IF(IFTRAK .NE. 0) THEN
        READ(IFTRAK) (VOLSUR(JUN,1),JUN=-NSOUT,NREG)
        READ(IFTRAK) (MATALB(JUN),JUN=-NSOUT,NREG)
        READ(IFTRAK) ( NSCRP,JUN=1,NALBG)
        READ(IFTRAK) ( ASCRP,JUN=1,NALBG)
        READ(IFTRAK) DANG0,(DASCRP,IUN=2,NDIM),
     >               ((DASCRP,IUN=1,NDIM),JUN=2,NANGL)
        READ(IFTRAK) (DASCRP,JUN=1,NANGL)
      ELSE
        CALL LCMGET(IPTRK,'MATALB      ',MATALB)
        ALLOCATE(DSV(-NSOUT:NREG))
        CALL LCMGET(IPTRK,'SAreaRvolume',DSV)
        DO JJ=-NSOUT,0
          VOLSUR(JJ,1)=0.25*REAL(DSV(JJ))
        ENDDO
        DO JJ=1,NREG
          VOLSUR(JJ,1)=REAL(DSV(JJ))
        ENDDO
*----
*  Allocate memory for NXT tracking
*----
        ALLOCATE(DGMESH(-1:MAXMSH,4))
        CALL NXTXYZ(IPTRK,IPRNTP,NDIM,ITYPBC,MAXMSH,NUCELL,ABSC,DGMESH)
        ALLOCATE(IUNFLD(2,NBUCEL))
        ALLOCATE(DANGLT(NDIM,NANGL),DORITR(NDIM*(NDIM+1),NPLANE,NANGL),
     >  DWGTRK(NANGL),DVNOR(NREG))
        NAMREC='G00000001CUF'
        CALL LCMGET(IPTRK,NAMREC,IUNFLD)
        CALL LCMGET(IPTRK,'TrackingDirc',DANGLT)
        CALL LCMGET(IPTRK,'TrackingOrig',DORITR)
        CALL LCMGET(IPTRK,'TrackingWgtD',DWGTRK)
        CALL LCMGET(IPTRK,'VTNormalize ',DVNOR)
      ENDIF
      ALLOCATE(SIGTAL(-NSOUT:NREG,NSBG),SIGT00(-NSOUT:NREG,NSBG))
*
      SWNZBC= .FALSE.
      SWVOID= .FALSE.
      LPIJK= IPIJK.EQ.4
*----
*  PREPARE FOR MULTIGROUP CALCULATION
*----
      DO ISBG=2,NSBG
        DO IUN= -NSOUT, NREG
          VOLSUR(IUN,ISBG)=VOLSUR(IUN,1)
        ENDDO
      ENDDO
      DO ISBG=1,NSBG
        IF(NPSYS(ISBG).NE.0) THEN
          DO ISA=1,6
            ALBEDO(ISA)=ALBG(ISA)
          ENDDO
          IF(NALBP .GT. 0) THEN
            DO ISA=1,6
              IF(ICODE(ISA).GT.0) ALBEDO(ISA)=ALBP(ICODE(ISA),ISBG)
            ENDDO
          ENDIF
          DO IUN= -NSOUT, -1
            SIGT00(IUN,ISBG)= 0.0
            SIGTAL(IUN,ISBG)= ALBEDO(-MATALB(IUN))
            SWNZBC= SWNZBC.OR.(SIGTAL(IUN,ISBG).NE.0.0)
          ENDDO
          IUN=0
          SIGT00(IUN,ISBG)= 0.0
          SIGTAL(IUN,ISBG)= 0.0
          DO IUN= 1, NREG
            SIGT00(IUN,ISBG)= XSSIGT(MATCOD(IUN),ISBG)
            SIGTAL(IUN,ISBG)= XSSIGT(MATCOD(IUN),ISBG)
            IF( SIGTAL(IUN,ISBG) .EQ. 0.0 )THEN
              SWVOID= .TRUE.
            ELSE
              VOLSUR(IUN,ISBG)=VOLSUR(IUN,ISBG)*SIGTAL(IUN,ISBG)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
*----
*  CHOOSE ISOTROPIC OR SPECULAR B.C.
*----
      IF( KSPEC.EQ.0 )THEN
         SIGT => SIGT00
      ELSE
         SIGT => SIGTAL
      ENDIF
*
      NPROB = (NUNKNO*(NUNKNO+1))/2
      N2PROB= NUNKNO*NUNKNO
      IF(IPRNTP .GT. 1) THEN
         NALLOC=(2*N2PROB*NSBG)
         IF(LPIJK) NALLOC=NALLOC+(2*N2PROB*NSBG)
         WRITE(IOUT,6000) NALLOC/256
      ENDIF
      ALLOCATE(DPROB(N2PROB,NSBG))
      CALL XDDSET(DPROB,N2PROB*NSBG,0.0D0)
      IF(LPIJK)THEN
        ALLOCATE(DPROBX(N2PROB,NSBG))
        CALL XDDSET(DPROBX,N2PROB*NSBG,0.0D0)
      ELSE
        ALLOCATE(DPROBX(1,1))
      ENDIF
      IF(IPRNTP.GT.1) WRITE(IOUT,6001)
*
      IF(IPRNTP .GE. 10) WRITE(IOUT,6010) MXSEG
      ALLOCATE(NUMERO(MXSEG),LENGTH(MXSEG))
      IF( ISPEC.EQ.0 )THEN
*----
*  Standard tracking
*----
         IF( NDIM.EQ.2 )THEN
           ALLOCATE(LOPATH(MXSEG))
           DO ITRAK= 1, NTRK
             IF(IFTRAK .NE. 0) THEN
*----
*  Read tracks from file
*----
               READ(IFTRAK)  NSUB,NBSEG,WEIGHT,IANG,
     >         (NUMERO(IL),IL=1,NBSEG),(LENGTH(IL),IL=1,NBSEG)
               IF(NSUB.NE.1) CALL XABORT('EXCELP: NSUB.NE.1.')
             ELSE
*----
*  Generate selected track
*----
               CALL NXTTGS(IPTRK ,IPRNTP,NDIM  ,NANGL ,NPOINT,NTRK  ,
     >                     ITRAK ,MAXMSH,NSOUT ,NREG  ,NUCELL,NBUCEL,
     >                     MXGSUR,MXGREG,MAXPIN,MXSEG ,ITYPBC,IUNFLD,
     >                     MATALB,DSV   ,DGMESH,DANGLT,DVNOR ,DWGTRK,
     >                     DORITR,NBSEG ,NCORT ,WEIGHT,NUMERO,LENGTH)
               IF(NCORT .EQ. 0) GO TO 1005
             ENDIF
             CALL PIJI2D(NREG,NSOUT,NBSEG,NCOR,NSBG,SWVOID,SIGT,NPSYS,
     >                   WEIGHT,LENGTH,NUMERO,LOPATH,DPROB,
     >                   MKI1,BI1,PAS1,L1,
     >                   MKI2,BI2,PAS2,XLIM2,L2,
     >                   MKI3,BI3,PAS3,XLIM3)
             IF(LPIJK)THEN
               CALL PIJI2D(NREG,NSOUT,NBSEG,NCOR,NSBG,SWVOID,SIGT,NPSYS,
     >                     WEIGHT,LENGTH,NUMERO,LOPATH,DPROBX,
     >                     MKI3,BI3,PAS3,L3,
     >                     MKI4,BI4,PAS4,XLIM4,L4,
     >                     MKI5,BI5,PAS5,XLIM5)
             ENDIF
 1005        CONTINUE
           ENDDO
           DEALLOCATE(LOPATH)
         ELSE
           ALLOCATE(STAYIN(MXSEG),GOSOUT(MXSEG))
           DO ITRAK= 1, NTRK
             IF(IFTRAK .NE. 0) THEN
               READ(IFTRAK)  NSUB,NBSEG,WEIGHT,IANG,
     >         (NUMERO(IL),IL=1,NBSEG),(LENGTH(IL),IL=1,NBSEG)
               IF(NSUB.NE.1) CALL XABORT('EXCELP: NSUB.NE.1.')
             ELSE
               CALL NXTTGS(IPTRK ,IPRNTP,NDIM  ,NANGL ,NPOINT,NTRK  ,
     >                     ITRAK ,MAXMSH,NSOUT ,NREG  ,NUCELL,NBUCEL,
     >                     MXGSUR,MXGREG,MAXPIN,MXSEG ,ITYPBC,IUNFLD,
     >                     MATALB,DSV   ,DGMESH,DANGLT,DVNOR ,DWGTRK,
     >                     DORITR,NBSEG ,NCORT ,WEIGHT,NUMERO,LENGTH)
               IF(NCORT .EQ. 0) GO TO 1015
             ENDIF
             CALL PIJI3D(NREG,NSOUT,NBSEG,NCOR,NSBG,SWVOID,SIGT,NPSYS,
     >       WEIGHT,LENGTH,NUMERO,STAYIN,GOSOUT,DPROB)
 1015        CONTINUE
           ENDDO
           DEALLOCATE(GOSOUT,STAYIN)
           IF(LPIJK) CALL XABORT(NAMSBR//': 3D PIJK NOT SUPPORTED')
         ENDIF
      ELSEIF( ISPEC.EQ.1 )THEN
*----
*  CYCLIC TRACKING
*----
         RCUTOF= CUTOF
         IF( NDIM.EQ.2 )THEN
            IF( DANG0.EQ. 0.0D0 )THEN
               NGSS= NANGL/8
            ELSE
               NGSS= (NANGL/4+1)/2
            ENDIF
            CALL ALGPT( NGSS,0.0,1.0,XGSS,WGSS)
            ALLOCATE(SIGANG(NGSS,-NSOUT:NREG,NSBG),STAYIN(NGSS*MXSEG),
     >      GOSOUT(NGSS*MXSEG))
            DO JGSS= 1, NGSS
              YGSS= SQRT(1.0 - XGSS(JGSS)**2)
              WGSS(JGSS)= WGSS(JGSS) * YGSS
              XGSS(JGSS)= 1.0/YGSS
              WGSSX(JGSS)= WGSS(JGSS) / (XGSS(JGSS)**2)
              DO ISBG=1,NSBG
                IF(NPSYS(ISBG).NE.0) THEN
                  DO IUN= -NSOUT,NREG
                    IF( MATALB(IUN).LE.0 )THEN
                      SIGANG(JGSS,IUN,ISBG)= SIGT(IUN,ISBG)
                    ELSE
                      SIGANG(JGSS,IUN,ISBG)= SIGT(IUN,ISBG)*XGSS(JGSS)
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
*----
*  Loop over tracks
*  then loop over groups
*----
            DO ITRAK= 1, NTRK
              IF(IFTRAK .NE. 0) THEN
                READ(IFTRAK)  NSUB,NBSEG,WEIGHT,(IANGL,IL=1,NSUB),
     >          (NUMERO(IL),IL=1,NBSEG),(LENGTH(IL),IL= 1,NBSEG)
              ELSE
               CALL NXTTGC(IPTRK ,IPRNTP,NDIM  ,NANGL ,NPOINT,NTRK  ,
     >                     ITRAK ,MAXMSH,NSOUT ,NREG  ,NUCELL,NBUCEL,
     >                     MXGSUR,MXGREG,MAXPIN,MXSEG ,ITYPBC,IUNFLD,
     >                     MATALB,DSV   ,DGMESH,DANGLT,DVNOR ,DWGTRK,
     >                     DORITR,NBSEG ,NCORT ,WEIGHT,NUMERO,LENGTH)
                IF(NCORT .EQ. 0) GO TO 1025
              ENDIF
              CALL PIJS2D(NREG,NSOUT,NBSEG,NSBG,WEIGHT,RCUTOF,NGSS,
     >        SIGANG,XGSS,WGSS,NPSYS,LENGTH,NUMERO,STAYIN,GOSOUT,DPROB)
              IF(LPIJK)THEN
*               X-DIRECTION  PROBABILITIES CALCULATIONS ( PX=PY )
                CALL PIJS2D(NREG,NSOUT,NBSEG,NSBG,WEIGHT,RCUTOF,NGSS,
     >          SIGANG,XGSS,WGSSX,NPSYS,LENGTH,NUMERO,STAYIN,GOSOUT,
     >          DPROBX)
              ENDIF
 1025         CONTINUE
            ENDDO
            DEALLOCATE(GOSOUT,STAYIN,SIGANG)
         ELSE
            ALLOCATE(STAYIN(MXSEG),GOSOUT(MXSEG))
            DO ITRAK= 1, NTRK
              IF(IFTRAK .NE. 0) THEN
                READ(IFTRAK)  NSUB,NBSEG,WEIGHT,(IANGL,IL=1,NSUB),
     >          (NUMERO(IL),IL=1,NBSEG),(LENGTH(IL),IL=1,NBSEG)
              ELSE
               CALL NXTTGC(IPTRK ,IPRNTP,NDIM  ,NANGL ,NPOINT,NTRK  ,
     >                     ITRAK ,MAXMSH,NSOUT ,NREG  ,NUCELL,NBUCEL,
     >                     MXGSUR,MXGREG,MAXPIN,MXSEG ,ITYPBC,IUNFLD,
     >                     MATALB,DSV   ,DGMESH,DANGLT,DVNOR ,DWGTRK,
     >                     DORITR,NBSEG ,NCORT ,WEIGHT,NUMERO,LENGTH)
                IF(NBSEG .EQ. 0) GO TO 1035
              ENDIF
              CALL PIJS3D(NREG,NSOUT,NBSEG,NSBG,WEIGHT,RCUTOF,SIGT,
     >        NPSYS,LENGTH,NUMERO,STAYIN,GOSOUT,DPROBX)
 1035         CONTINUE
            ENDDO
            DEALLOCATE(GOSOUT,STAYIN)
         ENDIF
      ENDIF
      IF(IFTRAK .EQ. 0) THEN
        DEALLOCATE(DVNOR,DWGTRK,DORITR,DANGLT)
        DEALLOCATE(IUNFLD)
        DEALLOCATE(DGMESH,DSV)
        CALL LCMSIX(IPTRK,'NXTRecords  ',ILCMDN)
      ENDIF
*
      DO 2050 ISBG=1,NSBG
         IF(NPSYS(ISBG).EQ.0) GO TO 2050
         KSBG=(ISBG-1)*NUNKNO
         CALL PIJCMP(NREG,NSOUT,NCOR,DPROB(1,ISBG),
     >               VOLSUR(-NSOUT,ISBG),.FALSE.,DPROB(1,ISBG))
         IF(LPIJK)THEN
           CALL PIJCMP(NREG,NSOUT,NCOR,DPROBX(1,ISBG),
     >                 VOLSUR(-NSOUT,ISBG),.TRUE.,DPROBX(1,ISBG))
         ENDIF
 2050 CONTINUE
      DEALLOCATE(LENGTH,NUMERO)
*----
*  RENORMALIZE ALL ISOTROPIC PROBS WITH VARIOUS OPTIONS
*----
      DO 2060 ISBG=1,NSBG
      IF(NPSYS(ISBG).EQ.0) GO TO 2060
      IF( KSPEC.EQ.0 )THEN
         IF( NRENOR.EQ.1 )THEN
*
*           NORMALIZATION USING GELBARD SCHEME
            CALL PIJRGL(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                  DPROB(1,ISBG))
            IF(LPIJK) CALL PIJRGL(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                            DPROBX(1,ISBG))
         ELSEIF( NRENOR.EQ.2 )THEN
*
*           NORMALIZATION WORKING ON DIAGONAL COEFFICIENTS
            CALL PIJRDG(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                  DPROB(1,ISBG))
            IF(LPIJK) CALL PIJRDG(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                            DPROBX(1,ISBG))
         ELSEIF( NRENOR.EQ.3 )THEN
*
*           NORMALIZATION WORKING ON WEIGHT FACTORS TO KEEP DIAG = 0.0
            CALL PIJRNL(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                  DPROB(1,ISBG))
            IF(LPIJK) CALL PIJRNL(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                            DPROBX(1,ISBG))
         ELSEIF( NRENOR .EQ. 4 )THEN  ! ATTENTION
*
*           NORMALIZATION WORKING ON WEIGHT FACTORS ADDITIVE (HELIOS)
            CALL PIJRHL(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                  DPROB(1,ISBG))
            IF(LPIJK) CALL PIJRHL(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                            DPROBX(1,ISBG))
         ENDIF
         IF( IPRT.GE.ICPALL )THEN
            LOPT= -1
            MSYM=1
            WRITE(IOUT,'(1H )')
            WRITE(IOUT,'(35H   COLLISION PROBABILITIES OUTPUT: ,
     >                   35H *BEFORE* ALBEDO REDUCTION          )')
            CALL PIJWPR(LOPT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                  DPROB(1,ISBG),VOLSUR(1,ISBG),MSYM)
*
            IF(LPIJK)THEN
              WRITE(IOUT,'(35H   X-DIRECT. COLL. PROBAB. OUTPUT: ,
     >                     35H *BEFORE* ALBEDO REDUCTION          )')
              CALL PIJWPR(LOPT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                   DPROBX(1,ISBG),VOLSUR(1,ISBG),
     >                   MSYM)
            ENDIF
*
         ENDIF
      ENDIF
 2060 CONTINUE
      IF(LPIJK)THEN
         DO 2070 ISBG=1,NSBG
         IF(NPSYS(ISBG).EQ.0) GO TO 2070
         CALL PIJD2S(NREG,NSOUT,DPROBX(1,ISBG),PROBKS(1,ISBG))
 2070    CONTINUE
      ENDIF
      IF( KSPEC.EQ.0 )THEN
*----
*  ELIMINATION OF SURFACES FOR PIJ
*----
         IF( SWNZBC )THEN
            ALLOCATE(PSST(NSOUT*NSOUT),PSVT(NSOUT*NREG),MATRT(NSOUT))
            CALL LCMLEN(IPTRK,'BC-REFL+TRAN',ILONG,ITYPE)
            IF(ILONG.EQ.NSOUT) THEN
              CALL LCMGET(IPTRK,'BC-REFL+TRAN',MATRT)
            ELSE
               WRITE(IOUT,9000) NAMSBR
               DO 130 ISOUT=1,NSOUT
                 MATRT(ISOUT)=ISOUT
 130           CONTINUE
            ENDIF
            DO 2080 ISBG=1,NSBG
              IF(NPSYS(ISBG).EQ.0) GO TO 2080
              CALL PIJABC(NREG,NSOUT,NPROB,SIGTAL(-NSOUT,ISBG),MATRT,
     >                    DPROB(1,ISBG),PSST,PSVT)
*----
*  ELIMINATION OF SURFACES FOR PIJX AND CREATION OF PIJXX
*----
            IF(LPIJK)THEN
               CALL PIJAAA(NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                     DPROBX(1,ISBG),PSVT,PROBKS(1,ISBG))
               CALL PIJABC(NREG,NSOUT,NPROB,SIGTAL(-NSOUT,ISBG),MATRT,
     >                     DPROBX(1,ISBG),PSST,PSVT)
            ENDIF
 2080     CONTINUE
*
            DEALLOCATE(MATRT,PSVT,PSST)
         ENDIF
      ENDIF
*
      ALLOCATE(FFACT(NREG))
      DO 2090 ISBG=1,NSBG
      IF(NPSYS(ISBG).EQ.0) GO TO 2090
      IF( IPRT.GE.ICPEND )THEN
         LOPT= +1
         MSYM=1
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'(35H   COLLISION PROBABILITIES OUTPUT: ,
     >                35H *AFTER* ALBEDO REDUCTION          )')
         CALL PIJWPR(LOPT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >               DPROB(1,ISBG),VOLSUR(1,ISBG),MSYM)
*
         IF(LPIJK)THEN
           WRITE(IOUT,'(35H   X-DIRECT. COLL. PROBAB. OUTPUT: ,
     >                  35H *AFTER* ALBEDO REDUCTION          )')
           CALL PIJWPR(LOPT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                 DPROBX(1,ISBG),VOLSUR(1,ISBG),MSYM)
           WRITE(IOUT,'(35H0 X-DIRECT. COLL. PROBAB." OUTPUT: ,
     >                  35H PIJX"=PIJX+PISX*(1/(1-PSS))*PSJ   )')
           MSYM=0
           CALL PIJWPR(LOPT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                 DPROBX(1,ISBG),VOLSUR(1,ISBG),MSYM)
         ENDIF
*
      ENDIF
*----
*  CHARGE PIJ MATRIX IN THE DRAGON SYMMETRIZED FORMAT
*----
      DO 160 IIN=1,NREG
         IF(SIGTAL(IIN,ISBG).EQ.0.0) THEN
            FFACT(IIN)=1.0
         ELSE
            FFACT(IIN)=1.0/SIGTAL(IIN,ISBG)
         ENDIF
  160 CONTINUE
      CALL PIJD2R(NREG,NSOUT,DPROB(1,ISBG),FFACT,.FALSE.,NELPIJ,
     >            N2PROB,PIJ(1,1,ISBG))
*----
*  CHARGE PIJX AND PIJY MATRICES IN THE DRAGON SYMMETRIZED FORMAT
*  ( PIJX=PIJY ), AND PIJZ CALCULATION ( PIJZ=3*PIJ-PIJX-PIJY )
*  AND THE SAME FOR FULL MATRICES OF PIJX", PIJY" AND PIJZ"
*----
      IF(LPIJK)THEN
        CALL PIJD2R(NREG,NSOUT,DPROBX(1,ISBG),FFACT,.TRUE.,NELPIJ,
     >              N2PROB,PIJ(1,2,ISBG))
        IVV=0
        DO 181 IUN=1,NREG
          IU=IUN
          IL=(IUN-1)*NREG+1
          DO 191 JUN=1,IUN
            IVV=IVV+1
            PROBKS(IL,ISBG)=1.5*PROBKS(IL,ISBG)*FFACT(IUN)*FFACT(JUN)
            IF(IL.NE.IU)PROBKS(IU,ISBG)=1.5*PROBKS(IU,ISBG)*
     >                      FFACT(IUN)*FFACT(JUN)
            PIJ(IVV,3,ISBG)=PIJ(IVV,2,ISBG)
            PROBKS(NNREG+IL,ISBG)=PROBKS(IL,ISBG)
            PROBKS(NNREG+IU,ISBG)=PROBKS(IU,ISBG)
            PIJ(IVV,4,ISBG)=3*PIJ(IVV,1,ISBG)-PIJ(IVV,2,ISBG)
     >                                       -PIJ(IVV,3,ISBG)
            PROBKS(2*NNREG+IL,ISBG)=3*PIJ(IVV,1,ISBG)
     >      -PROBKS(IL,ISBG)-PROBKS(NNREG+IL,ISBG)
            PROBKS(2*NNREG+IU,ISBG)=3*PIJ(IVV,1,ISBG)
     >      -PROBKS(IU,ISBG)-PROBKS(NNREG+IU,ISBG)
            IU=IUN+JUN*NREG
            IL=IL+1
  191     CONTINUE
  181   CONTINUE
*----
*  COMPUTE PIJ**(-1)*PIJK*
*----
        CALL PIJKST(IPRNTP,NREGIO,PIJ(1,1,ISBG),PROBKS(1,ISBG))
      ENDIF
 2090 CONTINUE
      DEALLOCATE(FFACT,DPROBX)
*
      DEALLOCATE(DPROB,SIGT00,SIGTAL,MATALB,VOLSUR)
*----
*  CHECK IF SCATTERING REDUCTION IS REQUIRED
*----
      ALLOCATE(PCSCT(NREGIO,2*NREGIO))
      DO 3000 ISBG=1,NSBG
      IF(NPSYS(ISBG).EQ.0) GO TO 3000
      LSKIP=.TRUE.
      DO 200 IBM=1,NBMIX
        LSKIP=LSKIP.AND.(XSSIGW(IBM,1,ISBG).EQ.0.0)
  200 CONTINUE
*----
*  COMPUTE THE SCATTERING-REDUCED CP MATRICES
*----
      IOP=1
      IF(.NOT.LSKIP) THEN
        CALL PIJSMD(IPRNTP,NBMIX,NREGIO,MATCOD,VOLUME,XSSIGW(0,1,ISBG),
     >              XSSIGT(0,ISBG),LEAKSW,PIJ(1,1,ISBG),PCSCT,IOP)
        DO 220 I=1,NREGIO
          FACT=VOLUME(I)
          DO 210 J=1,NREGIO
            INDPIJ=INDPOS(I,J)
            PIJ(INDPIJ,1,ISBG)=REAL(PCSCT(I,J))*FACT
  210     CONTINUE
  220   CONTINUE
      ENDIF
*-------
      IF(IPIJK.EQ.4) THEN
        IOP=4
        IF(.NOT.LSKIP) THEN
*         P1 SCATTERING REDUCTION OF THE DIRECTIONNAL CP MATRICES.
          IF(NANI.LT.2) CALL XABORT('EXCELP: ANISOTROPIC SCAT MISSING.')
          DO 250 IJKS=1,3
            CALL PIJSMD(IPRNTP,NBMIX,NREGIO,MATCOD,VOLUME,
     >                  XSSIGW(0,2,ISBG),XSSIGT(0,ISBG),LEAKSW,
     >                  PIJ(1,IJKS+1,ISBG),PCSCT,IOP)
            DO 240 I=1,NREGIO
              FACT=VOLUME(I)
              DO 230 J=1,NREGIO
                INDPIJ=INDPOS(I,J)
                PIJ(INDPIJ,IJKS+1,ISBG)=REAL(PCSCT(I,J))*FACT
  230         CONTINUE
  240       CONTINUE
  250     CONTINUE
        ENDIF
      ENDIF
 3000 CONTINUE
      DEALLOCATE(PCSCT)
      RETURN
*
 6010 FORMAT(' Maximum length of a line =',I10)
 6000 FORMAT(' *** SPACE REQUIRED FOR CP MATRICES = ',I10,' K ***')
 6001 FORMAT(' *** CP MATRICES ALLOCATED            ',10X,'   ***')
 9000 FORMAT(1X,A6,': *** WARNING *** '/
     >       ' REFLECTION/TRANSMISSION MATRIX MISSING'/
     >       ' USE IDENTITY REFLECTION MATRIX')
      END
