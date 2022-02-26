*DECK XELDRV
      SUBROUTINE XELDRV(IPTRK ,IPGEOM,IPRT  ,MAXPTS,NANIS ,NORE  ,
     >                  LMERG ,KSPEC ,KTOPT ,TITREC,CUTOFX,CFTRAK,
     >                  IFTRAK,IDISP ,ISYMM ,LCACT ,NMU   ,INSB  ,
     >                  LBIHET,LPRISM,IZ,DELU,FRTM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read a Dragon tracking file to compute pij, normalize a tracking 
* file to Dragon format and produce a new tracking file in Dragon 
* format.
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
* IPTRK   pointer to the excell tracking (L_TRACK).
* IPGEOM  pointer to the geometry (L_GEOM).
* IPRT    tracking print level.
* MAXPTS  number of zones according to user.
* NANIS   anisotropy of the solution.
* NORE    track normalization flag (<=0: yes; =1: no).
* LMERG   type of volume normalization.
* KSPEC   kind of pij integration (=0: isotr.; =1: spec.).
* KTOPT   tracking type option.
* TITREC  title for this case.
* CUTOFX  mfp cutoff for specular integration.
* CFTRAK  name of the sequential binary tracking file.
* IFTRAK  unit of the sequential binary tracking file.
* IDISP   mode of the sequential binary tracking file.
* LCACT   type of polar integration for the method of characteristics.
* NMU     number of polar angles for the method of characteristics.   
* ISYMM   symmetry factor.
* INSB    type of vectorization for the calculation of CP matrices.
* LBIHET  activation flag for the double heterogeneity option.
* LPRISM  flag for 3D prismatic geometry.
* IZ      projection axis for 3D prismatic geometry.
* DELU    user defined track spacing for 3D prismatic tracking.
* FRTM    minimum volume fraction of the grain in the representative 
*         volume for She-Liu-Shi model.
*
*-----------------------------------------------------------------------
*
      USE             GANLIB
      IMPLICIT        NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER       TITREC*72,CFTRAK*12
      TYPE(C_PTR)     IPTRK,IPGEOM
      INTEGER         IPRT,MAXPTS,NANIS,NORE,LMERG,KSPEC,KTOPT,IFTRAK,
     >                IDISP,ISYMM,LCACT,NMU,INSB,IZ
      REAL            CUTOFX,DELU,FRTM
      LOGICAL         LPRISM,LBIHET
*----
*  LOCAL VARIABLES
*----
      INTEGER         NSTATE
      PARAMETER      (NSTATE=40)
      INTEGER         NREG,NUNKNO,IUTYPE,ISTATE(NSTATE),IFTEMP,IFILE,
     >                ITYPE,ITOPT,ITROP,NCOMNT,NTRK,IREC,IC,IR,JR,IUNK,
     >                NSOUT,IDISPO,NDIM,NV,NS,NALBGO,NCOR,NANGL,MXSEG,
     >                NPRISM,NDIMO,NSO,NVO,NUNOLD,KDROPN,KDRCLS,NUNKNX,
     >                IOUT,ITGEO,IUSED(6),ICMAX,ICODE(6),ICOLD(6),
     >                NANGLO,MXSUB,MXSUBO,MXSEGO,ILONG,IFMT,I
      LOGICAL         LEAKSW, LELCHK, SWNOGE, SWCONS, EMPTY, LCM
      REAL            ALBEDO(6),ALBOLD(6),EXTKOP(NSTATE),
     >                ZERO,ONE,DENS,PCORN
      DOUBLE PRECISION DASCRP
      CHARACTER       GEONAM*12,CORIEN(0:3,6)*4,CUSED(0:1)*6,TEXT12*12,
     >                COMENT*80,CTRK*4
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MATCOD,KEYFLX,MATOLD,MATALB,
     > KEYMRG,MATMRG
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLUME,VOLOLD,VOLSUR,VOLMRG
*
      PARAMETER     ( IUTYPE=2, IOUT=6, ZERO=0.0, ONE=1.0 )
*
      DATA          (( CORIEN(JR,IR),IR=1,6),JR=0,3)
     >              / ' 1  ',' 2  ',' 3  ',' 4  ',' 5  ',' 6  ',
     >                '****',' R+ ','****','****',' Z- ',' Z+ ',
     >                'HBC ','****','****','****',' Z- ',' Z+ ',
     >                ' X- ',' X+ ',' Y- ',' Y+ ',' Z- ',' Z+ '/
      DATA           ( CUSED(JR),JR=0,1 ) / 'UNUSED','  USED' /
*----
*  SCRATCH STORAGE ALLOCATION
*   MATCOD  material numbers for zones in the supercell.
*   VOLUME  volumes of zones in the supercell.
*   KEYFLX  zone key for the unknown vectors (fluxes...).
*----
      ALLOCATE(MATCOD(MAXPTS),KEYFLX(MAXPTS))
      ALLOCATE(VOLUME(MAXPTS))
*
      NPRISM=0
      IF( IPRT.GT.0 ) WRITE(IOUT,'(1X,A72//)') TITREC
      SWNOGE= .NOT.C_ASSOCIATED(IPGEOM)
      IF(SWNOGE) THEN
         GEONAM=' '
         IF(INSB.EQ.2) CALL XABORT('XELDRV: GEOMETRY REQUESTED')
      ELSE
         CALL LCMINF(IPGEOM,GEONAM,TEXT12,EMPTY,ILONG,LCM)
      ENDIF
      NTRK  = 0
*
      IF( IPRT.GT.0 )THEN
         WRITE(IOUT,'(1H )')
         IF( SWNOGE )THEN
            WRITE(IOUT,'(27H ECHO = >>>   NO GEOMETRY     )')
         ELSE
            WRITE(IOUT,'(27H ECHO = >>> GEOMETRY NAME: ,A8)') GEONAM
         ENDIF
         WRITE(IOUT,'(27H ECHO = >>> TRACKING FILE: ,A8)') CFTRAK
         IF(     IDISP.EQ.-1 )THEN
            WRITE(IOUT,'(27H ECHO = >>>          DISP: ,A4)') 'MODT'
         ELSEIF( IDISP.EQ. 0 )THEN
            WRITE(IOUT,'(27H ECHO = >>>          DISP: ,A4)') 'OLDT'
         ELSEIF( IDISP.EQ.+1 )THEN
            WRITE(IOUT,'(27H ECHO = >>>          DISP: ,A4)') 'NEWT'
         ENDIF
         IF( NORE.EQ.-1 )THEN
            WRITE(IOUT,'(36H ECHO = >>> NORMALIZED-BY ANGLE     )')
         ELSEIF( NORE.EQ. 0 )THEN
            WRITE(IOUT,'(36H ECHO = >>> NORMALIZED-GLOBAL       )')
         ELSEIF( NORE.EQ.+1 )THEN
            WRITE(IOUT,'(36H ECHO = >>> NOT NORMALIZED          )')
         ENDIF
      ENDIF
*----
*  1)  REWIND TRACKING FILE --------------------------------------------
*----
      IF((IFTRAK.EQ.0).AND.(INSB.NE.2)) THEN
         CALL XABORT('XELDRV: NO SEQUENTIAL BINARY TRACKING FILE YET D'
     >   //'EFINED')
      ELSE IF((IFTRAK.NE.0).AND.(INSB.EQ.2)) THEN
         CALL XABORT('XELDRV: NO SEQUENTIAL BINARY TRACKING FILE EXPEC'
     >   //'TED WITH OPTION XCLL')
      ENDIF
      IF(INSB.NE.2) REWIND IFTRAK
*----
*  1.2) GET HEADER INFORMATIONS FROM *OLDT*/*MODT* FILES
*----
      NUNOLD= 0
      IFMT= 0
      IF( IDISP.LE.0 )THEN
         IREC= 1
         READ(IFTRAK,ERR=997) CTRK,NCOMNT,NTRK,IFMT
         DO 10 IC= 1, NCOMNT
            IREC= IREC+1
            READ (IFTRAK,ERR=997) COMENT
   10    CONTINUE
         IREC= IREC+1
         READ (IFTRAK,ERR=997) NDIMO,ITOPT,NVO,NSO,NALBGO,
     >                         NCOR,NANGL,MXSUB,MXSEG
         IF( NALBGO.LE.0.OR.NALBGO.GT.6 )THEN
            CALL XABORT('XELDRV: NALBG.GT.6.OR.NALBG.LE.0'//
     >                  ' ON TRACKING FILE')
         ENDIF
         NUNOLD= NVO+NSO+1
         ALLOCATE(VOLOLD(NUNOLD),MATOLD(NUNOLD))
         IREC= IREC+1
         READ (IFTRAK,ERR=997) (VOLOLD(IR),IR=1,NUNOLD)
         IREC= IREC+1
         READ (IFTRAK,ERR=997) (MATOLD(IR),IR=1,NUNOLD)
         IREC= IREC+1
         READ (IFTRAK,ERR=997) ( ICOLD(IR),IR=1,NALBGO)
         IREC= IREC+1
         READ (IFTRAK,ERR=997) (ALBOLD(IR),IR=1,NALBGO)
         IREC= IREC+1
         READ (IFTRAK,ERR=997) (DASCRP,IR=0,NDIMO*NANGL-1)
         IREC= IREC+1
         READ (IFTRAK,ERR=997) (DASCRP,IR=0,NANGL-1)
         REWIND IFTRAK
      ENDIF
*----
*  1.3) OPEN TEMP TRACKING FILE FOR *MODT*/*NEWT* FILES
*----
      ITROP = 0
      IF((IDISP.NE.0).AND.(INSB.NE.2)) THEN
         IFILE= KDROPN('DUMMYSQ',0,IUTYPE,0)
         IF( IFILE.LE.0 ) GO TO 998
         IF( IFILE.EQ.IFTRAK ) CALL XABORT('XELDRV: BAD TRACKING UNIT')
         IFTEMP = IFILE
         REWIND IFTEMP
         IF( IDISP.LT.0 )THEN
*
*           FOR *MODT* FILES, MAIN TRACKING IS COPIED ON TEMPORARY
            CALL XELCOP( IFTRAK, IFTEMP )
         ENDIF
      ENDIF
*----
*  2)  GET GEOMETRIC INFORMATIONS AND TRACK IF NECESSARY----------------
*----
      IF( SWNOGE )THEN
*----
*  2.1)  NO GEOMETRY, GET INFORMATIONS FROM *OLDT*/*MODT* FILE
*----
         IF(IDISP.GT.0) CALL XABORT('XELDRV: A RHS BINARY TRACKING FIL'
     >   //'E OR A RHS GEOMETRY MUST BE DEFINED.')
         NDIM=  NDIMO
         NSOUT= NSO
         NREG=  NVO
         NS=    NSO
         NV=    NVO
         NUNKNX= NUNOLD
         IF( NREG.GT.MAXPTS ) THEN
           WRITE(IOUT,'(28H ****** XELDRV ERROR  ******,/
     >      28H NUMBER OF REGION COMPUTED =,I10/
     >      28H MAXIMUM NUMBER OF REGION  =,I10)') NREG,MAXPTS
           CALL XABORT('XELDRV: MAXR TOO SMALL')
         ENDIF
         ALLOCATE(VOLSUR(NUNKNX),MATALB(NUNKNX),KEYMRG(NUNKNX))
         IUNK=  0
         DO 20 IR= -NS, NV
            VOLSUR(IUNK+1)= VOLOLD(IUNK+1)
            MATALB(IUNK+1)= MATOLD(IUNK+1)
            KEYMRG(IUNK+1)= IR
            IUNK= IUNK+1
   20    CONTINUE
         ALLOCATE(VOLMRG(NUNKNX),MATMRG(NUNKNX))
         DO 21 IUNK=1,NUNKNX
             VOLMRG(IUNK)= VOLOLD(IUNK)
             MATMRG(IUNK)= MATOLD(IUNK)
   21    CONTINUE
         DEALLOCATE(MATOLD,VOLOLD)
         DO 25 IR= 1, NALBGO
             ICODE(IR)=  ICOLD(IR)
   25    CONTINUE
         ITGEO= 0
      ELSE
*----
*  2.2)  THERE IS A GEOMETRY, GO TO EXCELL MODULES TO ANALYZE IT
*----
         CALL XDISET(ISTATE,NSTATE,0)
         CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
         ITYPE= ISTATE(1)
         IF( ITYPE.EQ.3.OR.ITYPE.EQ.6 )THEN
            ITGEO= 1
         ELSEIF( ITYPE.EQ. 8 .OR. ITYPE.EQ. 9 .OR.
     >           ITYPE.EQ.24 .OR. ITYPE.EQ.25 )THEN
            ITGEO= 2
         ELSEIF( ITYPE.EQ. 5 .OR. ITYPE.EQ. 7 .OR. ITYPE.EQ.20 .OR.
     >           ITYPE.EQ.21 .OR. ITYPE.EQ.22 .OR. ITYPE.EQ.23 )THEN
            ITGEO= 3
         ELSE
            ITGEO= 0
         ENDIF
         IDISPO= IDISP
         NANGLO= NANGL
         MXSUBO= MXSUB
         MXSEGO= MXSEG
         IF((INSB.EQ.2).AND.(ITGEO.NE.3)) THEN
            CALL XABORT('XELDRV: XCELL TRACKING NOT AVAILABLE.')
         ENDIF
         IF( ISTATE(13).GE.1 )THEN
           IF( ITYPE.EQ.3.OR.ITYPE.EQ.20.OR.ITYPE.EQ.24 )THEN
*----
*  2.2.1.1)     EXCELL DRIVER FOR CLUSTER SINGLE CELLS
*----
             CALL XCWTRK(IPTRK ,IPGEOM,GEONAM,IDISP ,IFTEMP,
     >                   IPRT  ,NDIM  ,ITOPT ,NV    ,NS    ,NANGL ,
     >                   ISYMM ,DENS  ,PCORN ,MXSUB ,MXSEG ,ICODE ,
     >                   TITREC)
             ITROP=3
           ELSE
             CALL XABORT('XELDRV: ONLY ONE-CELL TUBE/CARCEL/HEXCEL'//
     >                   ' CLUSTERS ARE AVAILABLE')
           ENDIF
         ELSEIF( ITGEO.EQ.2 )THEN
*----
*  2.2.1.2)  EXCELL DRIVER FOR HEXAGONAL 2D/3D ASSEMBLIES
*----
           CALL XHXTRK(IPTRK ,IPGEOM,GEONAM,IDISP,IFTEMP,
     >                 IPRT  ,NDIM  ,ITOPT ,NV   ,NS    ,NANGL ,
     >                 ISYMM ,DENS  ,PCORN ,MXSEG,ICODE ,TITREC)
           MXSUB=1
           ITROP=2
         ELSEIF( ITGEO.EQ.3 )THEN
*----
*  2.2.1.3)  EXCELL DRIVER FOR CARTESIAN 2D/3D ASSEMBLIES
*----
            CALL XELTRK(IPTRK ,IPGEOM,GEONAM,IDISP ,IFTEMP,
     >                  IPRT  ,NDIM  ,ITOPT ,NV    ,NS    ,NANGL ,
     >                  ISYMM ,DENS  ,PCORN ,MXSUB ,MXSEG ,ICODE ,
     >                  TITREC,INSB  ,IZ    ,LPRISM,NPRISM)
            ITROP=1
*----
*  For case with intrinsic symmetry
*  tracking performed on unfolded geometry assuming angular and spatial
*  symmetry
*  Normalization must be on global volume since only this option
*  makes sense.
*----  
            IF((ISYMM .GT. 1).AND.(LMERG .EQ. 0)) LMERG=1
         ELSE
            CALL XABORT('XELDRV: INVALID TYPE OF GEOMETRY')
         ENDIF
         IF((INSB.EQ.2).AND.(NDIM.NE.3)) THEN
            CALL XABORT('XELDRV: XCELL OPTION LIMITED TO 3D GEOMETRY.')
         ENDIF
*----
*  2.2.1.4)  RECOVER KEYMRG, MATALB AND VOLSUR
*----
         NUNKNX= NV+NS+1
         CALL LCMSIX(IPTRK,'EXCELL      ',1)
         CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
         IF(NUNKNX.NE.ISTATE(6)) CALL XABORT('XELDRV: INVALID NUNKNX.')
         ALLOCATE(VOLSUR(NUNKNX),MATALB(NUNKNX),KEYMRG(NUNKNX))
         CALL LCMGET(IPTRK,'KEYMRG      ',KEYMRG)
         CALL LCMGET(IPTRK,'MATALB      ',MATALB)
         CALL LCMGET(IPTRK,'VOLSUR      ',VOLSUR)
         CALL LCMSIX(IPTRK,'EXCELL      ',2)
*----
*  2.2.2)  MERGE SURFACES AND ZONES FOR THIS GEOMETRY
*----
         ALLOCATE(VOLMRG(NUNKNX),MATMRG(NUNKNX))
         CALL XELCMP( NS,   NV, VOLSUR, MATALB, KEYMRG,
     >               NSOUT,  NREG, VOLMRG, MATMRG, ITGEO, ICODE)
         IF( NREG.GT.MAXPTS ) THEN
           WRITE(IOUT,'(28H ****** XELDRV ERROR  ******,/
     >      28H NUMBER OF REGION COMPUTED =,I10/
     >      28H MAXIMUM NUMBER OF REGION  =,I10)') NREG,MAXPTS
           CALL XABORT('XELDRV: MAXR TOO SMALL')
         ENDIF
         IF((IPRT.GE.1).AND.(INSB.NE.2)) THEN
            WRITE(IOUT,6002) NANGL,ISYMM,CUTOFX,DENS,PCORN
         ENDIF
         IF(INSB.EQ.2) THEN
            IF(IDISP.NE.99) CALL XABORT('XELDRV: INCONSISTENT IDISP')
            CALL LCMGET(IPTRK,'ALBEDO',ALBOLD)
         ENDIF
*
         IF( IDISP.NE.IDISPO )THEN
           CALL XABORT('XELDRV: *HALT* OPTION REQUESTED '//
     >                 ' NO FURTHER CALCULATION IS POSSIBLE')
         ELSEIF( IDISP.GT.0 )THEN
           IF(INSB.NE.2) REWIND IFTEMP
         ELSE
           IF( NANGL.NE.NANGLO )THEN
             CALL XABORT('XELDRV: NOT POSSIBLE TO CHANGE '//
     >                   ' *NANGL* PARAMETER OF TRACKING FILE')
           ENDIF
           IF( MXSEG.NE.MXSEGO )THEN
             CALL XABORT('XELDRV: NOT POSSIBLE TO CHANGE '//
     >                   ' *MXSEG* PARAMETER OF TRACKING FILE')
           ENDIF
           IF( MXSUB.NE.MXSUBO )THEN
             CALL XABORT('XELDRV: NOT POSSIBLE TO CHANGE '//
     >                   ' *MXSUB* PARAMETER OF TRACKING FILE')
           ENDIF
         ENDIF
*----
*  2.2.3)  CHECK CONSISTENCY BETWEEN *MODT* FILE AND GEOMETRY BEFORE
*          MERGE
*----
         SWCONS=.FALSE.
         IF( IDISP.LT.0 )THEN
            SWCONS= LELCHK(NSO  ,NVO ,VOLOLD,MATOLD,ICOLD,
     >                     NS   ,NV  ,VOLSUR,MATALB,ICODE, 0)
            IF( IPRT.GT.0 )THEN
               WRITE(IOUT,'(1H )')
               IF( SWCONS )THEN
                  WRITE(IOUT,'(70H ECHO = >>>   CONSISTENCY BETWEEN '//
     >                     'TRACKING FILE AND UNMERGED GEOMETRY    /)')
               ELSE
                  WRITE(IOUT,'(70H ECHO = >>> INCONSISTENCY BETWEEN '//
     >                     'TRACKING FILE AND UNMERGED GEOMETRY    /)')
               ENDIF
            ENDIF
         ENDIF
*----
*  2.2.4)  CHECK CONSISTENCY BETWEEN *OLDT*/*MODT* FILE AND GEOMETRY
*          AFTER MERGE
*----
         IF( IDISP.EQ.0.OR.(IDISP.LT.0.AND.(.NOT.SWCONS)) )THEN
            IF( NDIMO.NE.NDIM )THEN
              CALL XABORT('XELDRV: DIMENSION (2-D/3-D) INCONSISTENCY')
            ENDIF
            SWCONS= LELCHK(NSO  ,NVO ,VOLOLD,MATOLD,ICOLD,
     >                     NSOUT,NREG,VOLMRG,MATMRG,ICODE,IPRT)
            IF( SWCONS )THEN
               IF( IPRT.GT.0 )THEN
                  WRITE(IOUT,'(1H )')
                  WRITE(IOUT,'(70H ECHO = >>>   CONSISTENCY BETWEEN '//
     >                     'TRACKING FILE AND   MERGED GEOMETRY    /)')
               ENDIF
            ELSE
               WRITE(IOUT,'(70H ECHO = >>> INCONSISTENCY BETWEEN '//
     >                  'TRACKING FILE AND   MERGED GEOMETRY    /)')
               CALL XABORT('XELDRV: INCONSISTENCY OF MERGED GEOMETRY '//
     >                    'WITH OLD TRACKING FILE ' )
            ENDIF
*
*           CONSISTENCY WITH MERGED GEOMETRY
*           COPY MERGED VOLUMES INTO VOLSUR ARRAY SINCE MERGE WAS DONE
            NS= NSOUT
            NV= NREG
            IUNK=  0
            DO 50 IR= -NS, NV
               VOLSUR(IUNK+1)= VOLOLD(IUNK+1)
               MATALB(IUNK+1)= MATOLD(IUNK+1)
               KEYMRG(IUNK+1)= IR
               IUNK= IUNK+1
   50       CONTINUE
         ENDIF
      ENDIF
*----
*  3)  NORMALIZE TEMPORARY FILE FOR *MODT*/*NEWT* FILES ----------------
*----
      IF((IDISP.NE.0).AND.(INSB.NE.2))THEN
*----
*  3.1)  WARNING IF THE FILE HAS *NTRK*.NE.0
*----
         IF( NTRK.NE.0 )THEN
            WRITE(IOUT,'(1H )')
            WRITE(IOUT,'(60H ECHO = >>> WARNING: TRACKING FILE'//
     >                 ' MAY ALREADY BE NORMALIZED             /)')
         ENDIF
         IF (LPRISM) THEN   
            CALL XELCTR(IFTEMP,IFTRAK,MXSUB,MXSEG,CUTOFX,ALBOLD)
         ELSE
            CALL XELNTR( NDIM, IFTEMP, IFTRAK, NORE, LMERG,
     >                   IPRT, NS, NV, VOLSUR, MATALB, KEYMRG,
     >                   NSOUT, NREG, VOLMRG, MATMRG, CUTOFX,
     >                   ITGEO, ICODE, ALBOLD, NANGL,  MXSUB,
     >                   MXSEG)
         ENDIF
         IFILE= IFTEMP
         IFTEMP= KDRCLS(IFTEMP,2)
         IF( IFTEMP.LT.0 ) GO TO 999
      ENDIF
*----
*  4)  CHARGE GEOMETRIC ALBEDOS & GET PHYSICAL ALBEDOS IF NECESSARY ----
*----
      ICMAX= 0
      DO 60 IR= 1, 6
        ALBEDO(IR)= ONE
        IUSED(IR)= 0
        ICMAX= MAX(ICMAX,ICODE(IR))
        IF( ICODE(IR).LT.0 ) ALBEDO(IR)= ALBOLD(-ICODE(IR))
   60 CONTINUE
      IF( ICMAX.GT.0 )THEN
         CALL XABORT('XELDRV: PHYSICAL ALBEDOS NOT IMPLEMENTED')
      ENDIF
      IF( KTOPT .EQ. -1) THEN
        KTOPT= ITOPT
      ENDIF
      IF(KSPEC .EQ. -1 ) THEN
        KSPEC= KTOPT
      ELSE
        KSPEC=MIN(KSPEC,KTOPT)
      ENDIF
*----
*  5)  STOCK INFORMATION (OUTPUT TO DRAGON DRIVER) ---------------------
*----
      DO 70 IR= 1, NREG
         KEYFLX(IR)= IR
         VOLUME(IR)= VOLMRG(IR+NSOUT+1)
         MATCOD(IR)= MATMRG(IR+NSOUT+1)
   70 CONTINUE
*
*     COMPUTE LEAKAGE SWITCH
      LEAKSW=.TRUE.
      DO 80 IR= -NSOUT, -1
         IUSED(-MATMRG(IR+NSOUT+1))= 1
         LEAKSW= LEAKSW .AND. ALBEDO(-MATMRG(IR+NSOUT+1)).EQ.ONE
   80 CONTINUE
      LEAKSW=.NOT.LEAKSW
      DEALLOCATE(MATMRG,VOLMRG)
      IF( (IDISP.LE.0).AND.(.NOT.SWNOGE) )THEN
         DEALLOCATE(MATOLD,VOLOLD)
      ENDIF
*
      IF( IPRT.GT.0 )THEN
         IF( IPRT.GT.1 )THEN
            WRITE(IOUT,'(8H  SIDE  ,2X,6(7X,A4))')
     >                                  (CORIEN(ITGEO,IR),IR=1,6)
            WRITE(IOUT,'(8H  GEOM #,2X,6(7X,I4.0))')
     >                                  (MAX(0,-ICODE(IR)),IR=1,6)
            WRITE(IOUT,'(8H  PHYS #,2X,6(7X,I4.0))')
     >                                  (MAX(0,ICODE(IR)),IR=1,6)
            WRITE(IOUT,'(8H  ALBEDO,2X,1P,6E11.4)')
     >                                  (ALBEDO(IR),IR=1,6)
            WRITE(IOUT,'(8H        ,2X,6(5X,A6))')
     >                                  (CUSED(IUSED(IR)),IR=1,6)
            WRITE(IOUT,'(1H )')
         ENDIF
         WRITE(IOUT,'(1H )')
         IF( KSPEC.EQ.0 )THEN
          WRITE(IOUT,'(40H ECHO = >>> ISOTROPIC CP CALCULATION       )')
         ELSEIF( KSPEC.EQ.1 )THEN
          WRITE(IOUT,'(40H ECHO = >>>  SPECULAR CP CALCULATION       )')
          IF( CUTOFX.EQ.ZERO )THEN
             WRITE(IOUT,'(27H ECHO = >>>    NO CUT-OFF     )')
          ELSE
             WRITE(IOUT,'(27H ECHO = >>> MFP.  CUT-OFF: ,1P,E11.4 )')
     >       CUTOFX
          ENDIF
         ENDIF
         WRITE(IOUT,'(28H ECHO = >>> NB. OF REGIONS: ,I5)') NREG
      ENDIF
*----
*  5.2) RELEASE SPACE ACCORDING TO INVERSE ORDER OF ALLOCATIONS
*----
      DEALLOCATE(KEYMRG,MATALB,VOLSUR)
*
*     SAVE GENERAL TRACKING INFORMATION.
      IF(NANIS.EQ.1) THEN
         NUNKNO= NREG
      ELSE
         IF( NDIM.EQ.1 )THEN
            NUNKNO= NANIS*NREG
         ELSE IF( NDIM.EQ.2 )THEN
            NUNKNO= ((NANIS*(NANIS+1))/2)*NREG
         ELSE
            NUNKNO= NANIS*NANIS*NREG
         ENDIF
      ENDIF
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      ISTATE(1)=NREG
      ISTATE(2)=NUNKNO
      IF(LEAKSW) THEN
         ISTATE(3)=0
      ELSE
         ISTATE(3)=1
      ENDIF
      IR=0
      DO 100 I=1,NREG
      IR=MAX(IR,MATCOD(I))
  100 CONTINUE
      ISTATE(4)=IR
      ISTATE(5)=NSOUT
      ISTATE(6)=NANIS
      ISTATE(7)=ITROP
      ISTATE(8)=NORE
      ISTATE(9)=KTOPT
      ISTATE(10)=KSPEC
      ISTATE(13)=LCACT
      ISTATE(14)=NMU
      ISTATE(16)=NDIM
      ISTATE(22)=INSB
      ISTATE(39)=NPRISM
      IF(LBIHET) ISTATE(40)=1
      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPTRK,'MATCOD',NREG,1,MATCOD)
      CALL LCMPUT(IPTRK,'VOLUME',NREG,2,VOLUME)
      CALL LCMPUT(IPTRK,'KEYFLX',NREG,1,KEYFLX)
      CALL LCMPUT(IPTRK,'ICODE',6,1,ICODE)
      CALL LCMPUT(IPTRK,'ALBEDO',6,2,ALBEDO)
      CALL LCMGET(IPTRK,'EXCELTRACKOP',EXTKOP)
      EXTKOP(39)=FRTM
      EXTKOP(40)=DELU
      CALL LCMPUT(IPTRK,'EXCELTRACKOP',NSTATE,2,EXTKOP)
      IF(IPRT.GE.1) THEN
        WRITE(IOUT,6000) IPRT,(ISTATE(IR),IR=1,12),ISTATE(22),
     >  ISTATE(16),ISTATE(39),ISTATE(40)
        WRITE(IOUT,6001)(EXTKOP(IR),IR=1,3)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(VOLUME)
      DEALLOCATE(KEYFLX,MATCOD)
      RETURN
*
 6000 FORMAT(' EXCELT PRINT LEVEL =',I8/
     > ' TRACK STATE-VECTOR'/' ------------------'/
     > ' NREG   =',I8,' (NUMBER OF REGIONS)'/
     > ' NUNKNO =',I8,' (NUMBER OF UNKNOWNS IN SYSTEM)'/
     > ' ILEAK  =',I8,' (LEAKAGE FLAG:0=PRESENT/1=ABSENT)'/
     > ' MAXMIX =',I8,' (NUMBER OF MIXTURES)'/
     > ' NSOUT  =',I8,' (NUMBER OF OUTER SURFACES)'/
     > ' NANIS  =',I8,' (FLUX ANISOTROPY ORDER)'/
     > ' ITROP  =',I8,' (GEOMETRY TYPE)'/
     > ' NORE   =',I8,' (TRACK NORMALIZATION OPTION)'/
     > ' KTOPT  =',I8,' (TYPE OF TRACKING:0=TISO/1=TSPC)'/
     > ' KSPEC  =',I8,' (TYPE OF BOUNDARY CONDITION)'/
     > ' NANGL  =',I8,' (USER-SUPPLIED NUMBER OF TRACKING ANGLES)'/
     > ' ISYMM  =',I8,' (USER-SUPPLIED TRACKING SYMMETRY FACTOR)'/
     > ' INSB   =',I8,' (TYPE OF VECTORIZATION:0=ONEG/1=ALLG/2=XCLL)'/
     > ' NDIM   =',I8,' (NUMBER OF GEOMETRIC DIMENSIONS)'/
     > ' NPRISM =',I8,' (NUMBER OF PLANS IN 3D PRISMATIC GEOMETRIES)'/
     > ' IBIHET =',I8,' (0/1=DOUBLE HETEROGENEITY IS NOT/IS ACTIVE)'/
     > ' -----------------'/)
 6001 FORMAT(
     > ' EXCELL TRACK OPTIONS '/
     > ' CUTOFX =',F20.8,' (CUTOFF FOR TRACK LENGTH)'/
     > ' DENS   =',F20.8,' (TRACK DENSITY)'/
     > ' PCORN  =',F20.8,' (CORNER DUPLICATION DISTANCE)'/
     > ' -----------------'/)
 6002 FORMAT(
     > ' RECOMPUTED PARAMETERS '/
     > ' NANGL  =',I10  ,' (NUMBER OF TRACKING ANGLES)'/
     > ' ISYMM  =',I10  ,' (TRACKING SYMMETRY FACTOR)'/
     > ' CUTOFX =',F10.5,' (CUTOFF FOR TRACK LENGTH)'/
     > ' DENS   =',F10.5,' (TRACK DENSITY)'/
     > ' PCORN  =',F10.5,' (CORNER DUPLICATION DISTANCE)'/
     > ' -----------------'/)
*
  997 WRITE(IOUT,'(31H ERROR= RECORD DESTROYED...        )')
      WRITE(IOUT,'(31H ERROR= UNABLE TO READ  RECORD ,I10)') IREC
      WRITE(IOUT,'(31H ERROR=              ON FILE FT,I2.2)') IFILE
      CALL XABORT('XELDRV: READ  TRACKING FILE FAILED' )
  998 WRITE(IOUT,'(31H ECHO = UNABLE TO OPEN  FILE FT,I4)') IFILE
      CALL XABORT('XELDRV: OPEN FAILED')
  999 WRITE(IOUT,'(31H ECHO = UNABLE TO CLOSE FILE FT,I4)') IFILE
      CALL XABORT('XELDRV: CLOSE FAILED')
*
      END
