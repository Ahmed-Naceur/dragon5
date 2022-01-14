*DECK XELNTR
      SUBROUTINE XELNTR(   NDIM,  IFOLD, IFTRAK,   NORE,  LMERG,
     >                     IPRT,     NS,     NV,  VOLIN,  MATIN,
     >                    MRGIN,  NSOUT,  NVOUT, VOLOUT, MATOUT,
     >                   CUTOFX,  ITGEO,  ICODE, ALBEDO,  NANGL,
     >                      KIN,    LIN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute renormalized tracks to obtain true volume values. The file 
* IFOLD contains the old tracks while the file IFTRAK will contain the 
* normalized tracks.
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* NDIM    number of dimensions.                  
* IFOLD   unnormalized tracking file number (at input).     
* IFTRAK  normalized tracking file number (at output).    
* NORE    integer flag for normalization:               
*         -1 normalize to volume (angle dependent);  
*          0 normalize to volume (angle independent);
*          1 do not normalize.                       
* LMERG   second integer flag for normalization:        
*          0 preserve volumes of fine regions;       
*          1 preserve volumes of merged regions .    
* IPRT    intermediate printing level for prinout.     
* NS      number of surfaces before merging.                
* NV      number of zones before merging.                   
* VOLIN   volumes and surfaces before merging.           
* MATIN   material numbers before merging.                  
* MRGIN   merging index.                               
* NSOUT   number of surfaces after  merging.                
* NVOUT   number of zones after  merging.                   
* VOLOUT  volumes and surfaces after  merging.           
* MATOUT  material numbers before merging.                  
* CUTOFX  cutoff factor.                               
* ITGEO   kind of geometry.                            
* ICODE   index of boundary conditions.                
* NANGL   number of angles to renormalize tracks by angle.  
* KIN     max. number of subtracks in a single track.       
* LIN     max. number of track segments in a single track.  
*
*Parameters: output
* ALBEDO  geometric albedos on external faces.         
*
*-----------------------------------------------------------------------
*
      IMPLICIT           NONE
*
      INTEGER            NDIM,IFOLD,IFTRAK,NORE,LMERG,IPRT,NS,NV,
     >                   MATIN(-NS:NV),MRGIN(-NS:NV),NSOUT,NVOUT,
     >                   MATOUT(-NSOUT:NVOUT),ITGEO,ICODE(6),NANGL,
     >                   KIN,LIN
      INTEGER            IANG,IC,IL,IP,IR,ISPEC,ITRAK,IVS,IVSC,JR,
     >                   LINE,MNSUR,MXSUB,MXSEG,MXVOL,NANG2,NBTRK,
     >                   NCOMNT,IOUT,JL,NALBG,NCOR,NSCRP,NSURC,NSURM,
     >                   NVOLC,NVOLM,NTMP,NSREN,NVREN,IPREN,NCSEG,
     >                   NOLDS,NNEWS,IFMT,II,NSUB,IND,IREG,IVSMAX(2)
      REAL               ERRCUR
      REAL               VOLIN(-NS:NV),VOLOUT(-NSOUT:NVOUT),
     >                   ASCRP,ALBEDO(6),ERRSUR,ERRVOL,CUTOFX,VOLMIN,
     >                   TMPERR(10),ERRVM,ERRSM
      DOUBLE PRECISION   APRSUR,APRVOL,TOTVOL,TOTSUR,ZERO,ONE,TWO,
     >                   FOUR,HALF,QUART,HUND,PI,FACVOL,FACSUR,VOLREF,
     >                   AVGREN,WEIGHT,RCUT,DASCRP
      LOGICAL            LNEW
      CHARACTER          CTRK*4, COMENT*80, CORIEN(0:3,-6:0)*4
      PARAMETER        ( PI=3.14159265358979323846D0, IOUT=6,
     >                   ZERO=0.D0, ONE=1.D0, TWO=2.D0, FOUR=4.D0,
     >                   HUND=1.D2, HALF=0.5D0, QUART=0.25D0 )
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NRSEG,KANGL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DENSTY
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PATH
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: ANGLES
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: VOLTRK
      DATA         ((CORIEN(JR,IR),IR=-6,0),JR=0,3)
     >       / ' 6  ',' 5  ',' 4  ',' 3  ',' 2  ',' 1  ','    ',
     >         ' Z+ ',' Z- ','****','****',' R+ ','****','    ',
     >         ' Z+ ',' Z- ','****','****','****','HBC ','    ',
     >         ' Z+ ',' Z- ',' Y+ ',' Y- ',' X+ ',' X- ','    '/
*----
*  SCRATCH STORAGE ALLOCATION
*   VOLTRK: volumes & surfaces as computed by tracking.
*   ANGLES: x,y,z components of angles.
*   DENSTY: weights by angle.
*   PATH  : relative path of each segment in a track.
*   NRSEG : material identification in a track.
*----
      ALLOCATE(KANGL(KIN),NRSEG(LIN))
      ALLOCATE(ANGLES(3,NANGL),DENSTY(NANGL),PATH(LIN))
      ALLOCATE(VOLTRK(-NS:NV,0:NANGL))
*
*     READ FIRST RECORDS OF THE TRACKING FILE
      READ (IFOLD ) CTRK,NCOMNT,NSCRP,NSCRP
      DO 5 IC= 1, NCOMNT
         READ (IFOLD ) COMENT
    5 CONTINUE
      READ (IFOLD ) NSCRP,ISPEC,NSCRP,NSCRP,NALBG,NCOR,NSCRP,NSCRP,NSCRP
      IF( NALBG.LE.0.OR.NALBG.GT.6 )THEN
         CALL XABORT('XELNTR: NALBG.GT.6.OR.NALBG.LE.0'//
     >               ' ON TRACKING FILE')
      ENDIF
      READ (IFOLD )  (ASCRP,IR=-NS,NV)
      READ (IFOLD )  (NSCRP,IR=-NS,NV)
      READ (IFOLD )  (NSCRP,     IR=1,NALBG)
      READ (IFOLD )  (ALBEDO(IR),IR=1,NALBG)
      READ (IFOLD ) ((ANGLES(IR,JR),IR=1,NDIM),JR=1,NANGL)
      READ (IFOLD )  (DENSTY(JR),JR=1,NANGL)
*
      FACSUR= 0.0D0
      FACVOL= TWO
      IF( ISPEC.EQ.0 )THEN
         IF( NDIM.EQ.2 )THEN
            FACSUR= QUART*PI
         ELSEIF( NDIM.EQ.3 )THEN
            FACSUR= ONE
         ENDIF
      ELSEIF( ISPEC.EQ.1 )THEN
         IF( NDIM.EQ.2 )THEN
            FACSUR= HALF*PI
         ELSEIF( NDIM.EQ.3 )THEN
            FACSUR= ONE
         ENDIF
      ENDIF
*
*     INITIALIZE NORMALIZED FACTORS
*     NORE = -1 -> ANGLE DEPENDENT NORMALIZATION
*          =  0 -> ANGLE INDEPENDENT NORMALIZATION
*          =  1 -> NO NORMALIZE BUT FIND TRACK ERROR ON MERGED VOLUME
*     LMERG = 0 -> NORMALIZATION PRESERVE FINE VOLUMES
*           = 1 -> NORMALIZATION PRESERVE MERGED VOLUMES FROM KEYMRG
      NSREN=NSOUT
      NVREN=NVOUT
      IF(LMERG.EQ.0) THEN
        NSREN=NS
        NVREN=NV
      ENDIF
      DO 10 IVS= -NSREN, NVREN
        DO 11 IANG= 0, NANGL
           VOLTRK(IVS,IANG)= ZERO
   11 CONTINUE
   10 CONTINUE
*
*     COMPUTE CUTOFF FOR LINE RELATIVE TO MERGED VOLUME
      VOLMIN=VOLOUT(1)
      DO 12 IVS= 2, NVOUT
        VOLMIN= MIN(VOLMIN,VOLOUT(IVS))
   12 CONTINUE
      RCUT= VOLMIN*CUTOFX
      IF( IPRT.GE.999 )THEN
         WRITE(IOUT,'(11X,A32,F20.15)')
     >      'CUTOFF FACTOR FOR LINES    =    ',RCUT
      ENDIF
*
*     LOOP OVER TRACKING (UNNORMALIZED TRACKS)
*     AND COMPUTE VOLUME OF TRACK
      NBTRK= 0
      MXSUB= 0
      MXSEG= 0
      IF( IPRT .GE. 999 ) THEN
        WRITE(IOUT,'(A22)') 'INITIAL line segments ' 
      ENDIF
   20 CONTINUE
      READ(IFOLD,END=40) NSUB,LINE,WEIGHT,
     >                   (KANGL(II),II=1,NSUB),
     >                   (NRSEG(IL),IL=1,LINE),
     >                   (PATH(IL),IL=1,LINE)
      IF( NSUB.GT.KIN )THEN
         WRITE(IOUT,'(11X,A32,I20)')
     >   'NUMBER OF SUBTRACKS IN LINE =    ',NSUB
         WRITE(IOUT,'(11X,A32,I20)')
     >   'MAXIMUM NUMBER OF SUBTRACKS =    ',KIN
         CALL XABORT( 'XELNTR: TRACKING FILE CORRUPTED' )
      ELSE IF( LINE.GT.LIN )THEN
         WRITE(IOUT,'(11X,A32,I20)')
     >   'NUMBER OF ELEMENTS IN LINE =    ',LINE
         WRITE(IOUT,'(11X,A32,I20)')
     >   'MAXIMUM NUMBER OF ELEMENTS =    ',LIN
         CALL XABORT( 'XELNTR: TRACKING FILE CORRUPTED' )
      ENDIF
      MXSUB= MAX(MXSUB,NSUB)
      MXSEG= MAX(MXSEG,LINE)
      IF(ISPEC.EQ.1) THEN
*       ANGULAR-DEPENDENT INTEGRATION OF A CYCLIC MULTI-TRACK
        IND=0
        LNEW=.TRUE.
        DO 31 IL=1,LINE
          IREG=NRSEG(IL)
          IF(LMERG.EQ.0) THEN
            IPREN=IREG
          ELSE
            IPREN=MRGIN(IREG)
          ENDIF
          IF(IREG.GT.0) THEN
            IF(LNEW) THEN
              IND=IND+1
              IF(IND.GT.NSUB) CALL XABORT('XELNTR: NSUB OVERFLOW')
              LNEW=.FALSE.
            ENDIF
            IANG=KANGL(1)
            VOLTRK(IPREN,IANG)= VOLTRK(IPREN,IANG) + WEIGHT*PATH(IL)
          ELSE
            LNEW=.TRUE.
            IANG=KANGL(1)
            VOLTRK(IPREN,IANG)= VOLTRK(IPREN,IANG) + WEIGHT*PATH(IL)
          ENDIF
   31   CONTINUE
        IF(IND.NE.NSUB) CALL XABORT('XELNTR: ALGORITHM FAILURE')
      ELSE
        DO 32 IL  = 1, LINE
          IF(LMERG.EQ.0) THEN
            IPREN=NRSEG(IL)
          ELSE
            IPREN=MRGIN(NRSEG(IL))
          ENDIF
          IANG=KANGL(1)
          VOLTRK(IPREN,IANG)= VOLTRK(IPREN,IANG) + WEIGHT*PATH(IL)
   32   CONTINUE
      ENDIF
      NBTRK= NBTRK + 1
      IF( IPRT .GE. 999 ) THEN
        WRITE(IOUT,6100) NBTRK,IANG,LINE,WEIGHT,WEIGHT/DENSTY(IANG) 
        WRITE(IOUT,6101)
     >       (PATH(JL),NRSEG(JL),MRGIN(NRSEG(JL)),JL=1,LINE)
      ENDIF
      GO TO 20
*
*     COMPUTE TRACK NORMALIZATION FACTOR
   40 CONTINUE
      DO 47 IVS=  -NSREN, NVREN
        IF(LMERG.EQ.0) THEN
          VOLREF=DBLE(VOLIN(IVS))
        ELSE
          VOLREF=DBLE(VOLOUT(IVS))
        ENDIF
        DO 46 IANG= 1, NANGL
           VOLTRK(IVS,0)= VOLTRK(IVS,0) + VOLTRK(IVS,IANG)
           VOLTRK(IVS,IANG)= VOLTRK(IVS,IANG)*DENSTY(IANG)
           IF( VOLTRK(IVS,IANG).NE.ZERO )THEN
              VOLTRK(IVS,IANG)= VOLREF/VOLTRK(IVS,IANG)
           ELSE
              VOLTRK(IVS,IANG)= ONE
           ENDIF
   46   CONTINUE
   47 CONTINUE
      IF(NORE .EQ. 0) THEN
        DO 48 IVS= 1, NVREN
          AVGREN=DBLE(VOLOUT(IVS))/(VOLTRK(IVS,0)*FACVOL)
          DO 44 IANG= 1, NANGL
            VOLTRK(IVS,IANG)=AVGREN
   44     CONTINUE
   48   CONTINUE
      ENDIF
*
*     COMPUTE ERRORS ON VOLUMES
      TOTSUR=ZERO
      APRSUR=ZERO
      TOTVOL=ZERO
      APRVOL=ZERO
      ERRSM=0.0
      ERRVM=0.0
      IVSMAX(1)=0
      IVSMAX(2)=0
      IVSC=0 
      DO 50 IVS= -NSREN, NVREN
         IF( LMERG.EQ.0 )THEN
            VOLREF=DBLE(VOLIN(IVS))
         ELSE
            VOLREF=DBLE(VOLOUT(IVS))
         ENDIF
         IF( VOLTRK(IVS,0).EQ.ZERO .AND. VOLREF.GT.ZERO )THEN
            IF( IVS.LT.0 )THEN
               IF( LMERG.EQ.0 )THEN
                  WRITE(IOUT,9010) IVS
               ELSE
                  WRITE(IOUT,9011) IVS
               ENDIF
            ELSEIF( IVS.GT.0 )THEN
               IF( LMERG.EQ.0 )THEN
                  WRITE(IOUT,9000) IVS
               ELSE
                  WRITE(IOUT,9001) IVS
               ENDIF
            ENDIF
            IVSC=IVS
         ENDIF
         IF( IVS.LT.0 )THEN
            VOLTRK(IVS,0)= FACSUR*VOLTRK(IVS,0)
            IF( VOLIN(IVS).NE.0.0 )THEN
               ERRCUR=REAL(100.0*ABS(1.0-VOLTRK(IVS,0)/VOLREF))
               IF(ERRCUR .GT. ERRSM) THEN
                 IVSMAX(1)=IVS  
                 ERRSM=ERRCUR
               ENDIF
            ENDIF
            TOTSUR=TOTSUR+VOLREF
            APRSUR=APRSUR+VOLTRK(IVS,0)
         ELSEIF( IVS.GT.0 )THEN
            VOLTRK(IVS,0)= FACVOL*VOLTRK(IVS,0)
            TOTVOL=TOTVOL+VOLREF
            APRVOL=APRVOL+VOLTRK(IVS,0)
            IF( VOLREF.NE.ZERO )THEN
               ERRCUR=REAL(100.0*ABS(1.0-VOLTRK(IVS,0)/VOLREF))
               IF(ERRCUR .GT. ERRVM) THEN
                 IVSMAX(2)=IVS  
                 ERRVM=ERRCUR
               ENDIF
            ENDIF
         ENDIF
   50 CONTINUE
      ERRSUR=100.*(1.0-REAL(APRSUR/TOTSUR))
      ERRVOL=100.*(1.0-REAL(APRVOL/TOTVOL))
*
*     CONSTRUCT THE NEW TRACKING FILE
      REWIND IFOLD
      READ (IFOLD ) CTRK,NSCRP,NSCRP,IFMT
      WRITE(IFTRAK) CTRK,NCOMNT,NBTRK,IFMT
      DO 55 IC= 1, NCOMNT
         READ (IFOLD ) COMENT
         WRITE(IFTRAK) COMENT
   55 CONTINUE
      READ (IFOLD )  (NSCRP,IR=1,9)
      WRITE(IFTRAK) NDIM,ISPEC,NVOUT,NSOUT,NALBG,NCOR,NANGL,MXSUB,MXSEG
      READ (IFOLD )  (ASCRP,     IR=-NS,NV)
      WRITE(IFTRAK)  (VOLOUT(IR),IR=-NSOUT,NVOUT)
      READ (IFOLD )  (NSCRP,     IR=-NS,NV)
      WRITE(IFTRAK)  (MATOUT(IR),IR=-NSOUT,NVOUT)
      READ (IFOLD )  (NSCRP,     IR=1,NALBG)
      WRITE(IFTRAK)  ( ICODE(IR),IR=1,NALBG)
      READ (IFOLD )  (ASCRP,     IR=1,NALBG)
      WRITE(IFTRAK)  (ALBEDO(IR),IR=1,NALBG)
      READ (IFOLD ) ((DASCRP,     IR=1,NDIM ),JR=1,NANGL)
      WRITE(IFTRAK) ((ANGLES(IR,JR),IR=1,NDIM),JR=1,NANGL)
      READ (IFOLD )  (DASCRP,     JR=1,NANGL)
      WRITE(IFTRAK)  (DENSTY(JR),JR=1,NANGL)
      IF( IPRT .GE. 999 ) THEN
        WRITE(IOUT,'(A22)') 'FINAL line segments   ' 
      ENDIF
      DO 70 ITRAK=1, NBTRK
         READ(IFOLD) NSUB,LINE,WEIGHT,
     >               (KANGL(II),II=1,NSUB),
     >               (NRSEG(IL),IL=1,LINE),
     >               (PATH(IL),IL=1,LINE)
         IF(NSUB.GT.MXSUB) CALL XABORT('XELNTR: MXSUB overflow.')
         IF(RCUT .GT. 0.0)THEN
            IL= 0
   23       CONTINUE
               IF(IL.EQ.LINE) GO TO 25
               IL= IL+1
               IF(PATH(IL).LT.RCUT)THEN
                  IF(IL.NE.LINE)THEN
                     DO 24 JL= IL+1, LINE
                        NRSEG(JL-1)= NRSEG(JL)
                        PATH(JL-1)=  PATH(JL)
   24                CONTINUE
                  ELSE
                     LINE= LINE-1
                     GO TO 25
                  ENDIF
                  LINE= LINE-1
                  IL= IL-1
               ENDIF
               GO TO 23
   25       CONTINUE
         ENDIF
*
*        RENORMALIZE TRACK LENGHTS
         IF((NORE.EQ.-1) .AND. (NSUB.GT.1) ) THEN
*          ANGULAR-DEPENDENT NORMALIZATION OF A CYCLIC MULTI-TRACK
           IND=0
           LNEW=.TRUE.
           DO 56 IL=1,LINE
             IREG=NRSEG(IL)
             IF(LMERG.EQ.0) THEN
                IPREN=IREG
             ELSE
                IPREN=MRGIN(IREG)
             ENDIF
             IF(IREG.GT.0) THEN
               IF(LNEW) THEN
                 IND=IND+1
                 IF(IND.GT.NSUB) CALL XABORT('XELNTR: NSUB overflow')
                 LNEW=.FALSE.
               ENDIF
               IANG=KANGL(1)
               PATH(IL)= VOLTRK(IPREN,IANG) * PATH(IL)
             ELSE
               LNEW=.TRUE.
             ENDIF
   56      CONTINUE
           IF(IND.NE.NSUB) CALL XABORT('XELNTR: Algorithm failure')
         ELSE IF(NORE.LT.1 ) THEN
            DO 60 IL  = 1, LINE
               IF( NRSEG(IL).GT.0 )THEN
                  IF(LMERG.EQ.0) THEN
                     IPREN=NRSEG(IL)
                  ELSE
                     IPREN=MRGIN(NRSEG(IL))
                  ENDIF
                  IANG=KANGL(1)
                  PATH(IL)= VOLTRK(IPREN,IANG) * PATH(IL)
               ENDIF
   60       CONTINUE
         ENDIF
*
*        CHANGE #ING ACCORDING TO MERGES
         DO 65 IL  = 1, LINE
             NRSEG(IL)= MRGIN(NRSEG(IL))
   65    CONTINUE
*
*        START MODIFICATIONS 98/06/02
*        COMPRESS TRACKING FILE FOR SUCCESSIVE IDENTICAL REGIONS
         NOLDS=NRSEG(1)
         NCSEG=1
         DO 66 IL = 2, LINE
            NNEWS=NRSEG(IL)
            IF( NNEWS.LT.0 .OR. NNEWS.NE.NOLDS )THEN
               NOLDS=NNEWS
               NCSEG=NCSEG+1
               NRSEG(NCSEG)=NRSEG(IL)
               PATH(NCSEG)=PATH(IL)
            ELSEIF( NNEWS.EQ.NOLDS )THEN
               PATH(NCSEG)=PATH(NCSEG)+PATH(IL)
            ENDIF
 66      CONTINUE
         WRITE(IFTRAK) NSUB,NCSEG,WEIGHT,
     >                (KANGL(II),II=1,NSUB),
     >                (NRSEG(IL),IL=1,NCSEG),
     >                (PATH(IL),IL=1,NCSEG)
         IF( IPRT.GE.999 ) THEN
              WRITE(IOUT,'(2H #,I8,12H WRITE IANG=,I8,5H LEN=,I8)')
     >                       ITRAK,            IANG,      NCSEG
              WRITE(IOUT,'(1P,(1X,E15.6,1X,I8))' )
     >                 (PATH(JL),NRSEG(JL),JL=1,NCSEG)
         ENDIF
*       END MODIFICATIONS 98/06/02
*
   70 CONTINUE
      IF(IPRT .GE. 5) THEN
         MNSUR = -NSREN
         MXVOL =  NVREN
         NSURC = -1
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,7000) ERRSUR,ERRSM
         WRITE(IOUT,7005) -IVSMAX(1),4.0*VOLTRK(IVSMAX(1),0)
         DO 80 IP  = 1, (9 - MNSUR) / 10
            NSURM= MAX( MNSUR, NSURC-9 )
            WRITE(IOUT,7100)(' FACE',-IR,IR=NSURC,NSURM,-1)
            IF(LMERG.EQ.0) THEN
              WRITE(IOUT,7110)
     >             (4.*VOLIN(IR),IR=NSURC,NSURM,-1)
              WRITE(IOUT,7111)
     >             (CORIEN(ITGEO,MATIN(IR)),IR=NSURC,NSURM,-1)
            ELSE
              WRITE(IOUT,7110)
     >             (4.*VOLOUT(IR),IR=NSURC,NSURM,-1)
              WRITE(IOUT,7111)
     >             (CORIEN(ITGEO,MATOUT(IR)),IR=NSURC,NSURM,-1)
            ENDIF
            WRITE(IOUT,7101) (FOUR*VOLTRK(IR,0),IR=NSURC,NSURM,-1)
            NTMP=0
            DO 81  IR=NSURC,NSURM,-1
              IF(LMERG.EQ.0) THEN
                VOLREF=DBLE(VOLIN(IR))
              ELSE
                VOLREF=DBLE(VOLOUT(IR))
              ENDIF
              NTMP=NTMP+1
              IF(VOLREF.NE.ZERO) THEN
                TMPERR(NTMP)=REAL(HUND-HUND*VOLTRK(IR,0)/VOLREF)
              ELSE
                TMPERR(NTMP)=0.0
              ENDIF
  81        CONTINUE
            WRITE(IOUT,7102) (TMPERR(IR),IR=1,NTMP)
            IF(LMERG.EQ.0) THEN
              WRITE(IOUT,7103) (' FACE',-MRGIN(IR),IR=NSURC,NSURM,-1)
            ENDIF
            WRITE(IOUT,7104)
            NSURC = NSURC - 10
   80    CONTINUE
         NVOLC= 1
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,7001) ERRVOL,ERRVM
         WRITE(IOUT,7006) IVSMAX(2),VOLTRK(IVSMAX(2),0)
         DO 90 IP  = 1, (9 + MXVOL) / 10
            NVOLM= MIN( MXVOL, NVOLC+9 )
            WRITE(IOUT,7100) (' ZONE',IR,IR=NVOLC,NVOLM)
            IF(LMERG.EQ.0) THEN
              WRITE(IOUT,7120)
     >             (VOLIN(IR),IR=NVOLC,NVOLM)
              WRITE(IOUT,7121)
     >             (' MIX ', MATIN(IR),IR=NVOLC,NVOLM)
            ELSE
              WRITE(IOUT,7120)
     >             (VOLOUT(IR),IR=NVOLC,NVOLM)
              WRITE(IOUT,7121)
     >             (' MIX ', MATOUT(IR),IR=NVOLC,NVOLM)
            ENDIF
            WRITE(IOUT,7101) (VOLTRK(IR,0),IR=NVOLC,NVOLM)
            NTMP=0
            DO 91  IR= NVOLC,NVOLM
              IF(LMERG.EQ.0) THEN
                VOLREF=DBLE(VOLIN(IR))
              ELSE
                VOLREF=DBLE(VOLOUT(IR))
              ENDIF
              NTMP=NTMP+1
              IF(VOLREF.NE.ZERO) THEN
                TMPERR(NTMP)=REAL(HUND-HUND*VOLTRK(IR,0)/VOLREF)
              ELSE
                TMPERR(NTMP)=0.0
              ENDIF
  91        CONTINUE
            WRITE(IOUT,7102) (TMPERR(IR),IR=1,NTMP)
            IF(LMERG.EQ.0) THEN
              WRITE(IOUT,7103) (' ZONE',MRGIN(IR),IR=NVOLC,NVOLM)
            ENDIF
            WRITE(IOUT,7104)
            NVOLC = NVOLC + 10
   90    CONTINUE
         IF(IPRT .GT. 5)THEN
           NVOLC= 1
           NANG2= NANGL+2
           WRITE(IOUT,'(1H )')
           IF( NORE.EQ.0 )THEN
              WRITE(IOUT,7004)
           ELSEIF(NORE.EQ.-1) THEN
              WRITE(IOUT,7002)
           ELSEIF(NORE.EQ.1) THEN
              WRITE(IOUT,7003)
           ENDIF
           DO 110 IP  = 1, (9 + MXVOL) / 10
             NVOLM= MIN( MXVOL, NVOLC+9 )
             WRITE(IOUT,7100) (' VOL ',IR,IR=NVOLC,NVOLM)
             IF(NORE .EQ. -2 ) THEN
               IANG=1
               WRITE(IOUT,7131)  
     >            (VOLTRK(IR,IANG),IR=NVOLC,NVOLM)
             ELSE
               DO 100 IANG= 1, NANGL
                 WRITE(IOUT,7130) IANG, 
     >            (VOLTRK(IR,IANG),IR=NVOLC,NVOLM)
  100          CONTINUE
             ENDIF
             WRITE(IOUT,7104)
             NVOLC = NVOLC + 10
  110      CONTINUE
         ENDIF
      ELSE IF(IPRT .GE. 1) THEN
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,7000) ERRSUR,ERRSM
         WRITE(IOUT,7005) -IVSMAX(1),4.0*VOLTRK(IVSMAX(1),0)
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,7001) ERRVOL,ERRVM
         WRITE(IOUT,7006) IVSMAX(2),VOLTRK(IVSMAX(2),0)
      ENDIF
      IF( IVSC.NE.0 )THEN
         WRITE(IOUT,9020)
         CALL XABORT( 'XELNTR: CHECK NUMBERING OR USE FINER TRACKING')
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(VOLTRK,NRSEG,KANGL,PATH,DENSTY,ANGLES)
      RETURN
*----
*  Formats
*----
 6100 FORMAT(1X,'TRACK # ',I8,3X,'ANGLE ',I8,
     >       1X,'WITH ',I8,1X,'SEGMENTS',
     >       3X,'TOTAL WEIGHT =',1P,E15.6,3X,'VOLUME WEIGHT =',E15.6)
 6101  FORMAT(1P,(1X,E15.6,1X,2I8))
 7000 FORMAT(/' TRACKING ERRORS ON SURFACE   AVERAGE ERROR: ',F10.4,
     >        ' % ',5X,'MAXIMUM ERROR: ',F10.4,' % ')
 7001 FORMAT( ' TRACKING ERRORS ON VOLUME    AVERAGE ERROR: ',F10.4,
     >        ' % ',5X,'MAXIMUM ERROR: ',F10.4,' % ')
 7002 FORMAT(/' ANGLE-BY-ANGLE RENORMALIZATION FACTORS: '/)
 7003 FORMAT(/' ANGLE-BY-ANGLE RENORMALIZATION FACTORS (NOT USED): '/)
 7004 FORMAT(/' GLOBAL RENORMALIZATION FACTORS: '/)
 7005 FORMAT(' MAXIMUM ERROR ON SURFACE=',I8,' WITH AREA  =',1P,E11.4)
 7006 FORMAT(' MAXIMUM ERROR IN REGION =',I8,' WITH VOLUME=',1P,E11.4)
 7100 FORMAT(10X,10(A5,I7))
 7101 FORMAT(' APPROX   ',1P,10E12.4)
 7102 FORMAT(' ERR(%)   ',10F12.5)
 7103 FORMAT(' MERGE TO ',10(A5,I7))
 7104 FORMAT(' ')
 7110 FORMAT(' SURFACE  ',1P,10E12.4)
 7111 FORMAT(' SIDE     ',10(A4,8X))
 7120 FORMAT(' VOLUME   ',1P,10E12.4)
 7121 FORMAT(' MIXTURE  ',10(A4,1X,I7))
 7130 FORMAT(' ANG ',I4,1X,1P,10E12.4)
 7131 FORMAT(10X,1P,10E12.4)
*
 9000 FORMAT(' *** WARNING - ORIGINAL VOLUME  # ',I10,' NOT TRACKED ')
 9001 FORMAT(' *** WARNING - MERGED VOLUME    # ',I10,' NOT TRACKED ')
 9010 FORMAT(' *** WARNING - ORIGINAL SURFACE # ',I10,' NOT TRACKED ')
 9011 FORMAT(' *** WARNING - MERGED SURFACE   # ',I10,' NOT TRACKED ')
 9020 FORMAT(' *** ERROR - ONE MERGED VOLUME OR SURFACE NOT TRACKED')
      END
