*DECK XCWTRK
      SUBROUTINE XCWTRK(IPTRK ,IPGEOM,GEONAM,IDISP ,IFTEMP,IPRT  ,
     >                  NDIM  ,ITOPT ,NVOL  ,NSUR  ,NANGL ,ISYMM ,
     >                  DENS  ,PCORN ,MXSUB ,MXSEG ,ICODE ,TITREC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Analyse cluster geometry and perform specular or isotropic 
* traking if required.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPTRK   pointer to the excell tracking.
* IPGEOM  pointer to the geometry.
* GEONAM  geometry name.
* IFTEMP  temporary tracking file.
* IPRT    print option.
* TITREC  title for execution.
*
*Parameters: input/output
* IDISP   tracking file disposition:
*         = -2 no traking - only analyse geometry
*              then abort (option halt);
*         = -1 modify tracking file;
*         =  0 old tracking file;
*         =  1 new tracking file.
*
*Parameters: output
* NDIM    number of physical dimensions.
* ITOPT   tracking option:
*         = 0 finite;   
*         = 1 cyclic.
* NVOL    number of physical regions.
* NSUR    number of outer surface.
* NANGL   number of angles.
* ISYMM   symmetry factor.
* DENS    track density.
* PCORN   corner proximity.
* MXSUB   maximum number of subtracks.
* MXSEG   maximum segment length.
* ICODE   albedo associated with face.
*
*-------------------------    XCWTRK    -------------------------------
*
      USE          GANLIB
      IMPLICIT     NONE
      INTEGER      IOUT,NALB,NSTATE
      CHARACTER    NAMSBR*6
      PARAMETER   (IOUT=6,NALB=6,NSTATE=40,
     >             NAMSBR='XCWTRK')
*----
*  ROUTINE PARAMETERS
*----

      TYPE(C_PTR)  IPTRK,IPGEOM
      INTEGER      IDISP ,IFTEMP,IPRT  ,NDIM  ,ITOPT ,NVOL  ,NSUR  ,
     >             NANGL ,ISYMM ,MXSUB ,MXSEG ,ICODE(NALB)
      REAL         DENS  ,PCORN
      CHARACTER    GEONAM*12,TITREC*7
*----
*  REDGET VARIABLES
*----
      INTEGER     ITYPLU,INTLIR
      CHARACTER   CARLIR*12
      REAL        REALIR  
      DOUBLE PRECISION DBLLIR
*----
*  LOCAL VARIABLES
*---- 
      LOGICAL      SWZERO
      CHARACTER    COMENT*80
      INTEGER      NCODE(NALB),IMS(NALB)
      REAL         ALBEDO(NALB) 
      INTEGER      ISTATE(NSTATE)
      REAL         EXTKOP(NSTATE)
      INTEGER      ILENGT,ITYLCM,NANGR ,NCOMNT,NCOR  ,NALBG,
     >             MSROD ,MAROD ,MNAN  ,NRT   ,NSURX ,NBAN ,
     >             NUNK  ,JJ    ,IHS
      REAL         COTE  ,RADMIN 
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEYMRG,MATALB,NRINFO,NRODS,
     > NRODR,NXRI
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLSUR,RAN,RODS,RODR
*----
*  DEFAULT TRACKING OPTIONS:
*----
      PCORN=0.0
      CALL XDISET(ISTATE,NSTATE,0)
      CALL XDRSET(EXTKOP,NSTATE,0.0)
      CALL LCMLEN(IPTRK,'STATE-VECTOR',ILENGT,ITYLCM)
      IF(ILENGT .LE. 0 .OR. ILENGT .GT. NSTATE) THEN
        ITOPT=0
        NANGR=15
        ISYMM=1
        DENS=0.0
      ELSE
        CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
        CALL LCMGET(IPTRK,'EXCELTRACKOP',EXTKOP)
        ITOPT=ISTATE(9)
        NANGR=ISTATE(11)
        ISYMM=ISTATE(12)
        DENS=EXTKOP(2)
      ENDIF
*----
*  READ THE NEW TRACKING OPTIONS.
*----
      IF(IDISP .LE. 0) GO TO 200
 100  CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': CHARACTER DATA EXPECTED.')
      IF((CARLIR .EQ. 'TISO') .OR. (CARLIR .EQ. 'TSPC')) THEN
        IF(CARLIR .EQ. 'TSPC') THEN
          ITOPT=1
          SWZERO=.TRUE.
        ELSE
          ITOPT=0
        ENDIF
*----
*  2-D QUADRATURE PARAMETERS (ANGLE AND SPACE).
*----
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .EQ. 3) THEN
          IF(ITOPT .EQ. 1 .AND. CARLIR .EQ. 'MEDI') THEN
            SWZERO=.FALSE.
          ELSE
            CALL XABORT('XCWTRK: *MEDI* KEYWORD EXPECTED.')
          ENDIF
          CALL REDGET(ITYPLU,NANGR,REALIR,CARLIR,DBLLIR)
        ENDIF
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >    ': INTEGER DATA EXPECTED.') 
        NANGR=INTLIR
        IF(NANGR.LT.2) CALL XABORT(NAMSBR//
     >    ': THE NUMBER OF ANGLES MUST BE > 1.')
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >    ': REAL DATA EXPECTED.')
        DENS=REALIR
      ELSE IF(CARLIR .EQ. 'HALT') THEN
*----
*  NO TRACKING OPTION
*----
        IDISP=-2
      ELSE IF(CARLIR .EQ. 'SYMM') THEN
*----
*  SYMMETRY FACTOR
*----
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >    ': INTEGER DATA EXPECTED.')
        ISYMM=INTLIR
      ELSE IF(CARLIR .EQ. ';') THEN
        NANGL=NANGR
        GO TO 200
      ELSE
        CALL XABORT(NAMSBR//': INVALID KEYWORD.')
      ENDIF
      GO TO 100
 200  CONTINUE
*----
*  Set NANGL for specular tracking to a valid value
*----
      IF(ITOPT .EQ. 1) THEN
        NANGL=MIN(30,NANGL)
        IF(NANGL .GT. 24) THEN
          NANGL = 30
        ELSE IF(NANGL .GT. 20) THEN
          NANGL = 24
        ELSE IF(NANGL .GT. 18) THEN
          NANGL = 20
        ELSE IF(NANGL .GT. 14) THEN
          NANGL = 18
        ELSE IF(NANGL .GT. 12) THEN
          NANGL = 14
        ELSE IF(NANGL .GT. 8) THEN
          NANGL = 12
        ELSE
          NANGL = 8
        ENDIF
        ISYMM=1
      ENDIF
*----
*  SAVE EXCELL SPECIFIC TRACKING INFORMATION.
*----
      ISTATE(1)=NVOL
      ISTATE(5)=NSUR
      ISTATE(9)=ITOPT
      ISTATE(11)=NANGR
      ISTATE(12)=ISYMM
      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,ISTATE)
      EXTKOP(2)=DENS
      CALL LCMPUT(IPTRK,'EXCELTRACKOP',NSTATE,2,EXTKOP)
*----
*  ANALYZE GEOMETRY AND STORE DESCRIPTION ON TRACKING STRUCTURE
*---- 
      CALL AXGXCW(IPGEOM,IPTRK ,IPRT  ,GEONAM,ISYMM )
*----
*  READ TRACKING STRUCTURE
*     KEYMRG   : INTEGER MERGE VECTOR
*     VOLSUR   : REAL VOLUME-SURFACE VECTOR
*     MATALB   : INTEGER MATERIAL-FACE VECTOR
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMSIX(IPTRK,'EXCELL      ',1)
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE       )
      NDIM  = ISTATE(1)
      NSUR  = ISTATE(2)
      NVOL  = ISTATE(3)
      NSURX = ISTATE(4)
      NBAN  = ISTATE(5)
      NUNK  = ISTATE(6)
      NRT   = ISTATE(7)
      MSROD = ISTATE(8)
      MAROD = ISTATE(9)
      MNAN  = ISTATE(10)
      ALLOCATE(KEYMRG(NUNK),VOLSUR(NUNK),MATALB(NUNK))
      ALLOCATE(NRINFO(2*MNAN),NRODS(3*NRT),NRODR(NRT),NXRI(NRT*NBAN))
      ALLOCATE(RAN(NBAN),RODS(2*NRT),RODR(MSROD*NRT))
      CALL LCMGET(IPTRK,'RAN         ',RAN   )
      IF(NSURX .EQ. 4)
     >CALL LCMGET(IPTRK,'COTE        ',COTE  )
      CALL LCMGET(IPTRK,'RADMIN      ',RADMIN)
      CALL LCMGET(IPTRK,'NRODS       ',NRODS )
      CALL LCMGET(IPTRK,'RODS        ',RODS  )
      CALL LCMGET(IPTRK,'NRODR       ',NRODR )
      CALL LCMGET(IPTRK,'RODR        ',RODR  )
      CALL LCMGET(IPTRK,'NRINFO      ',NRINFO)
      CALL LCMGET(IPTRK,'NXRI        ',NXRI  )
      CALL LCMGET(IPTRK,'KEYMRG      ',KEYMRG)
      CALL LCMGET(IPTRK,'MATALB      ',MATALB)
      CALL LCMGET(IPTRK,'VOLSUR      ',VOLSUR)
      CALL LCMSIX(IPTRK,'EXCELL      ',2)
      CALL LCMGET(IPTRK,'ALBEDO      ',ALBEDO)
      CALL LCMGET(IPTRK,'ICODE       ',ICODE )
      CALL LCMGET(IPTRK,'NCODE       ',NCODE )
      IF(ISYMM.GT.1) THEN
        DO 110 IHS=1,NALB
          IMS(IHS)=1
 110    CONTINUE 
      ELSE 
        DO 111 IHS=1,NALB
          IMS(IHS)=IHS
 111    CONTINUE
      ENDIF
      IF(IDISP .EQ. 1) THEN
        MXSUB=1
        MXSEG=4*(NBAN+1+NRT*MSROD*MAROD)
        IF(ITOPT .EQ. 1) THEN
          MXSUB=4*NANGL
          MXSEG=16*NANGL*MXSEG
        ENDIF
        NCOMNT=5
        NCOR=1
        NALBG=NALB
        WRITE(IFTEMP) '$TRK',NCOMNT,0,0
        COMENT='CREATOR  : DRAGON'
        WRITE(IFTEMP) COMENT
        COMENT='MODULE   : XCWTRK'
        WRITE(IFTEMP) COMENT
        COMENT='TYPE     : CLUSTER'
        WRITE(IFTEMP) COMENT
        COMENT='GEOMETRY : '//GEONAM
        WRITE(IFTEMP) COMENT
        COMENT=TITREC
        WRITE(IFTEMP) COMENT
        IF(ITOPT .EQ. 1) THEN
          WRITE(IFTEMP) NDIM,ITOPT,NVOL,NSUR,NALBG,NCOR,4*NANGL,MXSUB,
     >    MXSEG
        ELSE
          WRITE(IFTEMP) NDIM,ITOPT,NVOL,NSUR,NALBG,NCOR,NANGL,MXSUB,
     >    MXSEG
        ENDIF
        WRITE(IFTEMP) (VOLSUR(JJ),JJ=1,1+NSUR+NVOL)
        WRITE(IFTEMP) (MATALB(JJ),JJ=1,1+NSUR+NVOL)
        WRITE(IFTEMP) (ICODE(JJ),JJ=1,NALBG)
        WRITE(IFTEMP) (ALBEDO(JJ),JJ=1,NALBG)
*----
*  SET DEFAULT TRACKING DENSITY
*----
        IF(DENS .EQ. 0.0) DENS=5.0/RADMIN
        IF(ITOPT .EQ. 1) THEN
*----
*  SPECULAR TRACKING
*----
          CALL XCWSCL(NDIM  ,NSURX ,NVOL  ,NBAN  ,NRT   ,MSROD ,MAROD ,
     >                NANGL ,DENS  ,IFTEMP,IPRT  ,NCODE ,SWZERO,NRINFO,
     >                RAN   ,COTE  ,NRODS ,RODS  ,NRODR ,RODR  ,MXSUB ,
     >                MXSEG ,NXRI  ,IMS   )
          NANGL=4*NANGL
        ELSE
*----
*  ISOTROPIC TRACKING
*----
          CALL XCWICL(NDIM  ,NSURX ,NVOL  ,NBAN  ,NRT   ,MSROD ,MAROD ,
     >                NANGL ,DENS  ,ISYMM ,IFTEMP,IPRT  ,NRINFO,RAN   ,
     >                COTE  ,NRODS ,RODS  ,NRODR ,RODR  ,MXSEG ,NXRI  ,
     >                IMS)
        ENDIF
      ENDIF
*----
*  RELEASE BLOCKS FOR GEOMETRY
*----
      DEALLOCATE(RODR,RODS,RAN)
      DEALLOCATE(NXRI,NRODR,NRODS,NRINFO)
      DEALLOCATE(MATALB,VOLSUR,KEYMRG)
      RETURN
      END
