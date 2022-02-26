*DECK TLM
      SUBROUTINE TLM(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Create Matlab procedure to trace the integration lines
* generated with the NXT tracking module of DRAGON.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* NENTRY  number of data structures transfered to this module.
* HENTRY  name of the data structures.
* IENTRY  data structure type where:
*         =1 for LCM memory object;
*         =2 for XSM file;
*         =3 for sequential binary file;
*         =4 for sequential ASCII file.
* JENTRY  access permission for the data structure where:
*         =0 for a data structure in creation mode;
*         =1 for a data structure in modifications mode;
*         =2 for a data structure in read-only mode.
* KENTRY  data structure pointer.
*
*Comments:
* Instructions for the use of the TLM: module:
*   M-file.m := TLM: [ M-file.m ] VOLTRK TRKFIL  ::
*    [ EDIT [ iprint ] ]
*    [ NTPO   nplots   ]
*    (TLMget) ;
*   where
*     M-file.m : SEQ_ASCII file containing Matlab instructions.
*     VOLTRK   : read-only tracking data structure
*                (signature L_TRACK).
*     TRKFIL   : read-only sequential binary tracking file.
*     EDIT     : keyword to specify print level.
*     iprint   : print level. By default, iprint=1.
*     NTPO     : keyword to specify number of plots generated
*                by this execution.
*     nplots   : number of plots. By default, nplots=1.
*     (TLMget) : Processing options to select types of plots.
*                (read from input using the TLMGET routine).
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: and TLM: Modules,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='TLM   ')
      INTEGER          ILCMUP,ILCMDN
      PARAMETER       (ILCMUP=1,ILCMDN=2)
      INTEGER          NSTATE,NIPLP
      PARAMETER       (NSTATE=40,NIPLP=6)
      DOUBLE PRECISION DZERO,DONE
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0)
*----
*  Variables for input via REDGET
*----
      INTEGER          ITYPLU,INTLIR
      CHARACTER        CARLIR*72,CARLST*72
      REAL             REALIR
      DOUBLE PRECISION DBLLIR
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MATALB
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IPLP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DPLP,DGMESH,DANGLT,
     > DVNOR
*----
*  Local functions
*----
      INTEGER          TLMVPL
      INTEGER          IVALID
*----
*  Local variables
*----
      TYPE(C_PTR)      IPTRK
      INTEGER          IMTRK,IFTRK,IMFTRK,IMMAT,IPMAT
      INTEGER          IEN
      CHARACTER        HSIGN*12,CMAT*4
      INTEGER          ISTATT(NSTATE),IEDIMG(NSTATE)
      CHARACTER        TITLE*72
      INTEGER          IPRINT,NPLOTS
      CHARACTER        CTYPE*4,COMNT*80
      INTEGER          NCOMNT,NBTR,ICOM
      INTEGER          NDIM,ISPEC,NREG,NSOUT,NALBG,NCOR,NANGL,IFMT,
     >                 MXSUB,MXSEG
      INTEGER          NBUCEL,NUCELL(3),MAXMSH,MAXMSP,MXGREG,MAXMDH
      INTEGER          II,KK,ITRKT,IRENOT,NBDR,NSUR,IPLOT
      INTEGER          NSKTRK
      INTEGER          ITGEO
      DOUBLE PRECISION XYZL(2,3)
      LOGICAL          LMIX
*----
*  Validate entry parameters
*----
      IF(NENTRY .NE. 3) CALL XABORT(NAMSBR//
     >  ': Three data structures required')
      IPTRK=C_NULL_PTR
      IMTRK=0
      IFTRK=0
      IMFTRK=0
      IEN=1
      IMMAT=JENTRY(IEN)
      IPMAT=FILUNIT(KENTRY(IEN))
      IF(IENTRY(IEN) .NE. 4 ) CALL XABORT(NAMSBR//
     >  ': Matlab .m file is not an ASCII file')
      IF(IMMAT .NE. 0   .AND.
     >   IMMAT .NE. 1 ) CALL XABORT(NAMSBR//
     >  ': Matlab .m file not in create or update mode')
*----
*  Scan data structure to determine type and mode
*----
      DO IEN=2,NENTRY
        IF(IENTRY(IEN) .EQ. 1 .OR. IENTRY(IEN) .EQ. 2) THEN
          IF(JENTRY(IEN) .EQ. 2) THEN
            CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
            IF(HSIGN .EQ. 'L_TRACK') THEN
              IPTRK=KENTRY(IEN)
              IMTRK=IEN
              CALL LCMGTC(KENTRY(IEN),'TRACK-TYPE',12,1,HSIGN)
              IF((HSIGN .NE. 'EXCELL').AND.(HSIGN .NE. 'MCCG')) THEN
                 CALL XABORT(NAMSBR//
     >           ': Tracking data structure type is invalid')
              ENDIF
            ELSE
              CALL XABORT(NAMSBR//
     >        ': Invalid signature for '//HENTRY(IEN))
            ENDIF
          ELSE
            CALL XABORT(NAMSBR//
     >      ': Tracking data structure not in read-only mode')
          ENDIF
        ELSE IF(IENTRY(IEN) .EQ. 3) THEN
          IF(JENTRY(IEN) .NE. 2) CALL XABORT(NAMSBR//
     >        ': Tracking file not in read-only mode')
          IFTRK=FILUNIT(KENTRY(IEN))
          IMFTRK=IEN
        ELSE
          CALL XABORT(NAMSBR//
     >    ': Invalid data structure format for '//HENTRY(IEN))
        ENDIF
      ENDDO
      IF(IMFTRK .EQ. 0) CALL XABORT(NAMSBR//
     >': No tracking file available')
      IF(IMTRK .EQ. 0) CALL XABORT(NAMSBR//
     >': No Tracking data structure available')
*----
*  Recover EDIT level and number of track processing option
*  [ EDIT [ iprint ] ]
*  [ NTPO [ nplots ] ]
*  by default iprint=1 and nplots=1.
*----
      IPRINT=1
      NPLOTS=1
      LMIX=.FALSE.
 1010 CONTINUE
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
 1011   CONTINUE
        IF(ITYPLU .EQ. 10) GO TO 1015
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': Read error -- Character variable expected')
        IF(CARLIR(1:4) .EQ. ';') THEN
          GO TO 1015
        ELSE IF(CARLIR(1:4) .EQ. 'EDIT') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) GO TO 1011
          IPRINT=INTLIR
        ELSE IF(CARLIR(1:7) .EQ. 'MIXTURE') THEN
          LMIX=.TRUE.
        ELSE IF(CARLIR(1:4) .EQ. 'NTPO') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) GO TO 1011
          NPLOTS=INTLIR
        ELSE
          GO TO 1015
        ENDIF
        GO TO 1010
 1015 CONTINUE
*----
*  Get Matlab plot options
*----
      CARLST=CARLIR
      ALLOCATE(IPLP(NIPLP,NPLOTS),DPLP(4*NPLOTS))
      CALL TLMGET(IPRINT,NPLOTS,NDIM,CARLST,IPLP,DPLP)
*----
*  Read tracking file parameters
*----
      READ(IFTRK) CTYPE,NCOMNT,NBTR,IFMT
      IF(CTYPE .NE. '$TRK') CALL XABORT(NAMSBR//
     >': Binary file is not a valid NXT: tracking file')
      IF(IFMT .NE. 1) CALL XABORT(NAMSBR//': IFMT.NE.1')
      ITRKT=1
      IRENOT=1
      DO ICOM=1,NCOMNT
        READ(IFTRK) COMNT
        IF(COMNT(1:12) .EQ. 'TRKNOR      ' ) THEN
          IF(COMNT(15:26) .EQ. 'Directional ' ) THEN
            IRENOT=-1
          ELSE IF(COMNT(15:26) .EQ. 'Global      ' ) THEN
            IRENOT=0
          ENDIF
        ELSE IF(COMNT(1:12) .EQ. 'OPTION      ' ) THEN
          IF(COMNT(15:26) .EQ. 'Extended    ' ) THEN
           ITRKT=0
          ENDIF
        ENDIF
      ENDDO
      IF(ITRKT .NE .0) CALL XABORT(NAMSBR//
     >': Insufficient information on tracking file'//
     >' Use EDIT -1000 in NXT:')
      READ(IFTRK) NDIM,ISPEC,NREG,NSOUT,NALBG,NCOR,NANGL,MXSUB,MXSEG
      READ(IFTRK) KK
      ALLOCATE(MATALB(-NSOUT:NREG))
      READ(IFTRK) (MATALB(II),II=-NSOUT,NREG)
      DO II=1,4
        READ(IFTRK) KK
      ENDDO
      NSKTRK=NCOMNT+8
      REWIND(IFTRK)
*----
*  Initialize tracking parameters to 0
*----
      CALL XDISET(ISTATT,NSTATE,0)
      CALL XDISET(IEDIMG,NSTATE,0)
*----
*  Read state vectors
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATT)
      CALL LCMGTC(IPTRK,'TITLE',72,1,TITLE)
      IF(ISTATT( 7) .NE. 4) CALL XABORT(NAMSBR//
     >': Tracking data structure incompatible with current module')
      NSUR=ISTATT(5)
      IRENOT=ISTATT(8)
      MXSEG=MAX(MXSEG,ISTATT(18))
      CALL LCMSIX(IPTRK,'NXTRecords  ',ILCMUP)
      CALL LCMGET(IPTRK,'G00000001DIM',IEDIMG)
      IF(ISTATT( 1) .NE. NREG  .OR.
     >   ISTATT( 5) .NE. NSOUT .OR.
     >   ISTATT(21) .NE. NANGL .OR.
     >   IEDIMG( 1) .NE. NDIM       ) THEN
        WRITE(IOUT,9000) ISTATT( 1),NREG ,ISTATT( 5),NSOUT,
     >                   ISTATT(21),NANGL,
     >                   IEDIMG( 1),NDIM
        CALL XABORT(NAMSBR//
     >': Tracking data structure and file do mot match')
      ENDIF
      ITGEO=ABS(IEDIMG( 2))
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6000) NAMSBR,NREG,NSOUT,NANGL,NDIM,ITGEO
      ENDIF
      NBUCEL=IEDIMG( 5)
      NUCELL(1)=IEDIMG(13)
      NUCELL(2)=IEDIMG(14)
      NUCELL(3)=IEDIMG(15)
      MAXMSH   =IEDIMG(16)
      MAXMSP   =IEDIMG(20)
      MXGREG   =IEDIMG(25)
      MAXMDH=MAX(MAXMSH,MAXMSP,MXGREG)
      ALLOCATE(DGMESH((MAXMDH+2)*4))
      CALL TLMGEO(IPTRK,IPMAT,IPRINT,ITGEO,MAXMDH,NDIM,NUCELL,
     >            DGMESH,XYZL)
      ALLOCATE(DANGLT(NDIM*NANGL))
      CALL LCMGET(IPTRK,'TrackingDirc',DANGLT)
      NBDR=1
      IF(IRENOT .EQ. -1) NBDR=NBDR+NANGL
      ALLOCATE(DVNOR(NREG*NBDR))
      CALL XDDSET(DVNOR,NREG*NBDR,DONE)
      IF(IRENOT .EQ. -1) THEN
        CALL LCMGET(IPTRK,'VTNormalizeD',DVNOR(NREG+1))
      ELSE IF(IRENOT .EQ. 0) THEN
        CALL LCMGET(IPTRK,'VTNormalize ',DVNOR)
      ENDIF
      CALL LCMSIX(IPTRK,'NXTRecords  ',ILCMDN)
*----
*  Read IPMAT to end-of-file and
*  insert pause if in update mode
*----
      IF(IMMAT .EQ. 1) THEN
 1000   CONTINUE
          READ(IPMAT,'(A4)',END=1005) CMAT
        GO TO 1000
 1005   CONTINUE
        WRITE(IPMAT,7010)
      ENDIF
*----
*  Write execution comments on IPMAT
*----
      WRITE(IPMAT,7000) NAMSBR,HENTRY(IMTRK),HENTRY(IMFTRK),TITLE
*----
*  Loop over PLOTS
*----
      DO IPLOT=1,NPLOTS
        IF(ABS(IPLP(1,IPLOT)) .EQ. 1) THEN
*----
*  POINTS
*----
          CALL TLMPNT(IPMAT ,IFTRK ,IPRINT,NSKTRK,NBTR  ,NDIM  ,
     >                NREG  ,NSUR  ,MXSUB ,MXSEG ,NANGL ,NBDR  ,
     >                NPLOTS,IPLOT ,IPLP  ,DANGLT,DVNOR)
        ELSE IF(ABS(IPLP(1,IPLOT)) .EQ. 2) THEN
*----
*  DIRECTIONS
*----
          CALL TLMDIR(IPMAT ,IFTRK ,IPRINT,ISPEC, NSKTRK,NBTR  ,
     >                NDIM  ,NSOUT, NREG  ,MXSUB ,MXSEG ,NANGL ,
     >                NBDR  ,NPLOTS,IPLOT ,IPLP  ,DANGLT,DVNOR ,
     >                MATALB,LMIX  )
         ELSE IF(ABS(IPLP(1,IPLOT)) .EQ. 3) THEN
*----
*  PLANA
*  Test if plane is valid
*----
          IVALID=TLMVPL(NDIM,NANGL,NPLOTS,IPLOT,IPLP,DPLP,DANGLT,XYZL)
          IF(IVALID . GE. 0) THEN
            CALL TLMPLA(IPMAT ,IFTRK ,IPRINT,NSKTRK,NBTR  ,NDIM  ,
     >                  NREG  ,MXSUB ,MXSEG ,NANGL ,NBDR  ,
     >                  NPLOTS,IPLOT ,IPLP  ,DPLP  ,DANGLT,DVNOR)
          ENDIF
        ELSE IF(ABS(IPLP(1,IPLOT)) .EQ. 4) THEN
*----
*  PLANP
*  Test if plane is valid
*----
          IVALID=TLMVPL(NDIM,NANGL,NPLOTS,IPLOT,IPLP,
     >                  DPLP,DANGLT,XYZL)
          IF(IVALID . GE. 0) THEN
            CALL TLMPLP(IPMAT ,IFTRK ,IPRINT,NSKTRK,NBTR  ,NDIM  ,
     >                  NREG  ,MXSUB ,MXSEG ,NANGL ,NBDR  ,
     >                  NPLOTS,IPLOT ,IPLP  ,DPLP  ,DANGLT,DVNOR)
          ENDIF
        ELSE IF(ABS(IPLP(1,IPLOT)) .EQ. 5) THEN
*----
*  REGION
*----
          CALL TLMREG(IPMAT ,IFTRK ,IPRINT,NSKTRK,NBTR  ,NDIM  ,
     >                NSOUT ,NREG  ,MXSUB ,MXSEG ,NANGL ,NBDR  ,
     >                NPLOTS,IPLOT ,IPLP  ,DANGLT,DVNOR ,MATALB,
     >                LMIX  )
        ENDIF
      ENDDO
*----
*  Release memory
*----
      DEALLOCATE(DPLP,IPLP,MATALB,DVNOR,DANGLT,DGMESH)
*----
*  Processing finished, return
*----
      RETURN
*----
*  Matlab .m file format
*----
 6000 FORMAT('(* Output from --',A6,'-- follows '/
     >       '   NREG =',I10/
     >       '   NSUR =',I10/
     >       '   NANGL=',I10/
     >       '   NDIM =',I10/
     >       '   ITGEO=',I10,10X,'*)')
 7000 FORMAT('%'/
     >       '% File generated using    : ',6X,A6/
     >       '% Tracking structure name : ',A12/
     >       '% Tracking file name      : ',A12/
     >       '% Title                   : ',A72/
     >       '%')
 7010 FORMAT('pause ;')
 9000 FORMAT('NREG =',2I10/
     >       'NSUR =',2I10/
     >       'NANGL=',2I10/
     >       'NDIM =',2I10)
      END
