*DECK XELDCL
      SUBROUTINE XELDCL( IPGEOM, GEONAM,   NDIM, MAXGRI, LCLSYM, NBLOCK,
     >                    NTYPO,    LL1,    LL2,   IPRT, NTOTCO, MAXRO ,
     >                   NGEOME,   NTYP,  NGIDL,  NTIDL,  NUNKO,  CELLG,
     >                    NSURO,  NVOLO, IDLDIM, IDLGEO, KEYTRN, KEYGEO,
     >                   IDLTYP, KEYTYP, MRGCEL, IDLBLK)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Associate all blocks of a problem to their block types and generate  
* almost all  useful integer values that will describe the problem.
*
*Copyright:
* Copyright (C) 1987 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPGEOM  pointer to the geometry (l_geom).             
* GEONAM  geometry name.                                
* NDIM    number of dimensions.                             
* MAXGRI  number of grid cell in X/Y/Z directions.           
* LCLSYM  symmetry flags (0: no,-1/+1: syme,,-2/+2: ssym). 
* NBLOCK  number of blocks.                                 
* NTYPO   old number of types.                              
* LL1     upper diag switch.                           
* LL2     lower diag switch.                           
* IPRT    intermediate printing level for output.      
*
*Parameters: output
* NTOTCO  tot number of cylinders in all geometries.         
* MAXRO   max number of words to stock meshes.              
* NGEOME  number of geometries.                             
* NTYP    new number of types.                              
* NGIDL   lenght of geometric numbering.               
* NTIDL   lenght of type numbering.                    
* NUNKO   old number of unknowns.                           
* CELLG   to keep cell geometry names.                 
* NSURO   number of surfaces of each geometry.              
* NVOLO   number of zones of each geometry.                 
* IDLDIM  position of each geoemtry in cylinders numbering. 
* IDLGEO  position of each geometry in the             
*         geometry numbering scheme.        
* KEYTRN  turn key for each block.                     
* KEYGEO  geometric key for each type.                 
* IDLTYP  position of each type in numbering scheme.   
* KEYTYP  type key for each block.                     
* MRGCEL  merging key of each block.                   
* IDLBLK  position of each block in numbering scheme.  
*
*-----------------------------------------------------------------------
*
      USE                GANLIB
      IMPLICIT           NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)        IPGEOM
      INTEGER            NDIM, NBLOCK, NTYPO, IPRT, NTOTCO,  MAXRO,
     >                   NGEOME, NTYP, NGIDL, NTIDL, NUNKO 
      CHARACTER          GEONAM*12
      INTEGER            MAXGRI(3),LCLSYM(3),CELLG(3*NBLOCK),
     >                   NSURO(NBLOCK),NVOLO(NBLOCK),IDLDIM(NBLOCK),
     >                   IDLGEO(NBLOCK),KEYTRN(NBLOCK),KEYGEO(NBLOCK),
     >                   IDLTYP(NBLOCK),KEYTYP(NBLOCK),MRGCEL(NBLOCK),
     >                   IDLBLK(NBLOCK)
      LOGICAL            LL1, LL2
*----
*  EXTERNAL FUNCTIONS
*----
      CHARACTER*2        AXGTRN
*----
*  LOCAL VARIABLES
*----
      INTEGER            IOUT,NSTATE,MAXTUR
      CHARACTER          NAMSBR*6
      PARAMETER         (IOUT=6,NSTATE=40,MAXTUR=12,NAMSBR='XELDCL')
*----
*  LOCAL PARAMETERS
*----
      CHARACTER          GEOC1*12,GEOC2*12,BLANC*8,CPLAN*8,GEOCV*12
      INTEGER            MINGRI(3),MEDGRI(3),ISTATE(NSTATE)
      LOGICAL            LLSYM
      INTEGER            IKG
      INTEGER            IBLK, I3, ISUB2, ITYP, JTYP, IX, IY, IZ,
     >                   NMERG1, NMERG2, NMERG3, IOFF, IOF1, IOF2,
     >                   IMERG1, NNCYL, NNSUR, NNVOL, NXC, NXM, IGEO,
     >                   IB1, IB2, IC1, IC2, IT1, IT2, IR, NLINP,
     >                   NTYP2, IOLTYP, NPROB, IDLPRB, IP, MAXREM
      INTEGER            KMESH,NXYZ,ITC
      INTEGER            IOT1,IOT2
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITGEOM,CELLT
*----
*  DATA STATEMENTS
*----
      DATA BLANC        / ' ' /
*----
*  SCRATCH STORAGE ALLOCATION
*   ITGEOM: turn by cell types
*   CELLT : cell type names
*----
      ALLOCATE(ITGEOM(NBLOCK),CELLT(3*NTYPO))
*----
*  INITIALIZE BLOCK INFORMATION
*----
      DO 10 IBLK= 1, NBLOCK
           NSURO(IBLK)= 0
           NVOLO(IBLK)= 0
          IDLGEO(IBLK)= 0
          KEYGEO(IBLK)= 0
          IDLTYP(IBLK)= 0
          KEYTYP(IBLK)= 0
          ITGEOM(IBLK)= 0
          KEYTRN(IBLK)= 0
          IDLBLK(IBLK)= 0
 10   CONTINUE
      LLSYM=.FALSE.
      DO 20 I3= 1, 3
         MINGRI(I3)= MAXGRI(I3)
         MEDGRI(I3)= 0
         IF( ABS(LCLSYM(I3)) .EQ. 1 )THEN
            MINGRI(I3)= (MAXGRI(I3)+1)/2
            LLSYM=.TRUE.
            IF( LCLSYM(I3).EQ.-1 )THEN
               MEDGRI(I3)= MINGRI(I3)-1
            ENDIF
         ELSE IF(ABS(LCLSYM(I3)) .EQ. 2 ) THEN
            MINGRI(I3)= MAXGRI(I3)/2
            LLSYM=.TRUE.
            IF( LCLSYM(I3).EQ.-2 )THEN
               MEDGRI(I3)= MINGRI(I3)
            ENDIF
         ELSE IF(LCLSYM(I3) .NE. 0) THEN  
            WRITE(IOUT,'(1H0,A8,4H -->,3(I8,1X))') 'LCLSYM', LCLSYM
            CALL XABORT(NAMSBR//': LCLSYM NOT WELL DEFINED' )
         ENDIF
 20   CONTINUE
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      NXYZ=MAX(ISTATE(3),ISTATE(4),ISTATE(5))
      ISUB2= ISTATE(9)
      KMESH=ISTATE(6)
      MAXRO = 0
      NTOTCO= 0 
      IF( ISUB2.GT.0 )THEN
         CALL LCMGET(IPGEOM,'CELL',CELLT)
         CALL LCMLEN(IPGEOM,'MIX', NMERG1, ITYP) 
         CALL LCMGET(IPGEOM,'MIX', KEYTYP)
         DO 30 IMERG1=1,NMERG1
           IF( KEYTYP(IMERG1).GT.0 )CALL XABORT(NAMSBR//': GENERATING '
     >     //'CELLS EXPECTED')
           KEYTYP(IMERG1)=-KEYTYP(IMERG1)
           IKG=KEYTYP(IMERG1)
           WRITE(GEOCV,'(3A4)')
     >     CELLT(3*IKG-2),CELLT(3*IKG-1),CELLT(3*IKG)
           CALL LCMSIX(IPGEOM,GEOCV,1)
           CALL XDISET(ISTATE,NSTATE,0)
           CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
           NXYZ=MAX(NXYZ,ISTATE(3),ISTATE(4),ISTATE(5))
           CALL LCMSIX(IPGEOM,GEOCV,2)
 30      CONTINUE
         CALL LCMLEN(IPGEOM,'MERGE', NMERG2, ITYP)
         IF( NMERG2.EQ.0 )THEN
            DO 100 IMERG1= 1, NMERG1
               MRGCEL(IMERG1)= IMERG1
 100        CONTINUE
         ELSEIF( NMERG2.EQ.NMERG1 )THEN
            CALL LCMGET(IPGEOM,'MERGE', MRGCEL)
         ELSE
            CALL XABORT(NAMSBR//': MERGES ARE INCOMPATIBLE' )
         ENDIF
*
         CALL LCMLEN(IPGEOM,'TURN', NMERG3, ITYP)
         IF( NMERG3.EQ.0 )THEN
            DO 110 IMERG1= 1, NMERG1
               ITGEOM(IMERG1)= 1
 110        CONTINUE
         ELSEIF( NMERG3.EQ.NMERG1 )THEN
            CALL LCMGET(IPGEOM,'TURN', ITGEOM)
            DO 120 IMERG1= 1, NMERG3
               IF( MOD(ITGEOM(IMERG1),MAXTUR).EQ.0.OR.
     >             MOD(ITGEOM(IMERG1),MAXTUR).GT.8 )
     >         CALL XABORT(NAMSBR//': INVALID TURNS (NO HEX CODES)' )
 120       CONTINUE
         ELSE
            CALL XABORT(NAMSBR//': TURNS ARE INCOMPATIBLE' )
         ENDIF
         IF(LL1 .OR. LL2) THEN
*----
*  Process diagonal symmetries
*----
           CALL AXGDIA( IPGEOM,   IPRT, NBLOCK,  NTYPO,   NXYZ, KMESH ,
     >                  GEONAM,    LL1,   LL2,  MINGRI,  CELLT, KEYTYP,
     >                  ITGEOM)
         ENDIF
         IF( LLSYM )THEN
*----
*  process x-x, y-y and z-z symmetry
*  1) Unfold geometry
*  2) Analyse symmetry
*----
           DO 300 IZ=MINGRI(3),1,-1
             IOF1=(IZ-1)*MINGRI(1)*MINGRI(2)
             IOF2=(IZ+MEDGRI(3)-1)*MAXGRI(1)*MAXGRI(2)
             DO 310 IY=MINGRI(2),1,-1
               DO 320 IX=MINGRI(1),1,-1
                 IOT2=IOF2+(IY+MEDGRI(2)-1)*MAXGRI(1)+IX+MEDGRI(1)
                 IOT1=IOF1+(IY-1)*MINGRI(1)+IX 
                 IF(IOT2 .NE. IOT1) THEN
                   IF(KEYTYP(IOT2) .NE. 0)THEN
                     CALL XABORT(NAMSBR//': PROBLEMS TO UNFOLD')
                   ELSE
                     KEYTYP(IOT2)=KEYTYP(IOT1)
                     KEYTYP(IOT1)= 0
                     ITGEOM(IOT2)=ITGEOM(IOT1)
                     ITGEOM(IOT1)= 0
                   ENDIF
                 ENDIF
  320          CONTINUE
  310        CONTINUE
  300      CONTINUE
           CALL AXGSYM( IPGEOM,   IPRT, NBLOCK,  NTYPO,   NXYZ,
     >                  GEONAM, LCLSYM, MINGRI, MAXGRI,  CELLT,
     >                  KEYTYP, ITGEOM)
         ENDIF
*----
*  FIND ALL DIFFERENT GEOMETRIES
*----
         IF( IPRT.GT.1 )THEN
            WRITE(IOUT,'(1H )')
            NXC= 1
            WRITE(IOUT,'(25H ===> CELL TYPES ARE:    /)')
            DO 400 IP= 1, (9+ISUB2)/10
              NXM= MIN( ISUB2, NXC+9 )
              WRITE(IOUT,'(1H ,10(I8.8,4X))')
     >        (IB1,IB1=NXC,NXM)
              WRITE(IOUT,'(1H ,30A4)')
     >        (CELLT(3*IB1-2),CELLT(3*IB1-1),CELLT(3*IB1),IB1=NXC,NXM)
              NXC= NXC + 10
  400       CONTINUE
*
*           PRINTING ASSEMBLY MAP
            CPLAN= BLANC
            NLINP= 3+(MAXGRI(2)+1)*((9+MAXGRI(1))/10+1)
            DO 410 IZ=1,MAXGRI(3)
              WRITE(IOUT,'(1H )')
              IF(NDIM.EQ.3)THEN
                 WRITE(CPLAN,'(4H (Z=,I3,1H))') IZ
              ENDIF
              WRITE(IOUT,'(/32H UNFOLD TYPE CELL MAP FOR PLANE ,A8)')
     >        CPLAN
              NXC= 1
              DO 415 IP  = 1, (9 + MAXGRI(1)) / 10
                NXM= MIN( MAXGRI(1), NXC+9 )
                WRITE(IOUT,'(1X,A8,1X,10(A4,I3,A4))')
     >             CPLAN, (' X= ',IR,' ROT',IR=NXC,NXM)
                NXC = NXC + 10
  415         CONTINUE
              WRITE(IOUT,'(1H )')
              DO 420 IY=1,MAXGRI(2)
              IOFF=((IZ-1)*MAXGRI(2)+(IY-1))*MAXGRI(1)
              NXC= 1
                DO 425 IP  = 1, (9 + MAXGRI(1)) / 10
                  NXM= MIN( MAXGRI(1), NXC+9 )
                  WRITE(IOUT,'(1X,A4,I3,2H=>,10(I7,1X,A2,1X))')
     >                ' Y= ',IY,(KEYTYP(IOFF+IR),
     >                          AXGTRN(ITGEOM(IOFF+IR)),IR=NXC,NXM)
                  NXC = NXC + 10
  425           CONTINUE
              WRITE(IOUT,'(1H )')
  420         CONTINUE
  410       CONTINUE
         ENDIF
         NGEOME= 0
         DO 40 IB1= 1, NBLOCK
            IC1= KEYTYP(IB1)
            IT1= ITGEOM(IB1)
            IF(IC1.LE.0.OR.IC1.GT.ISUB2 )THEN
               CALL XABORT(NAMSBR//': INVALID TYPE #')
            ENDIF
            NGEOME= NGEOME + 1
            CELLG(3*NGEOME-2)= CELLT(3*IC1-2)
            CELLG(3*NGEOME-1)= CELLT(3*IC1-1)
            CELLG(3*NGEOME  )= CELLT(3*IC1  )
            KEYTRN(NGEOME)= IT1
            WRITE( GEOC1(1: 4),'(A4)') CELLT(3*IC1-2)
            WRITE( GEOC1(5: 8),'(A4)') CELLT(3*IC1-1)
            WRITE( GEOC1(9:12),'(A4)') CELLT(3*IC1  )
*           SEARCH FOR SIMILAR GEOMETRIES IN PREVIOUS ONES
            IF( IB1.NE.1 )THEN
               DO 41 IB2= 1, IB1-1
                  IC2= KEYTYP(IB2)
                  IT2= ITGEOM(IB2)
                  IF( IT1.NE.IT2 ) GO TO 41
                  WRITE( GEOC2(1: 4),'(A4)') CELLT( 3*IC2-2 )
                  WRITE( GEOC2(5: 8),'(A4)') CELLT( 3*IC2-1 )
                  WRITE( GEOC2(9:12),'(A4)') CELLT( 3*IC2   )
                  IF( GEOC1.EQ.GEOC2 )THEN
                     KEYGEO(IB1)= KEYGEO(IB2)
                     NGEOME= NGEOME-1
                     GO TO 40
                  ENDIF
 41            CONTINUE
            ENDIF
*----
*  ANALYSE NEW GEOMETRY
*----
            CALL LCMSIX(IPGEOM,GEOC1,1)
            CALL XELPRC(IPGEOM,GEOC1,NDIM,NNCYL,NNSUR,NNVOL,MAXREM)
            IF( NNVOL.NE.0 )THEN
               NSURO(NGEOME)= -NNSUR
               NVOLO(NGEOME)=  NNVOL
               IDLDIM(NGEOME)= NTOTCO
               NTOTCO= NTOTCO + NNCYL + 3
               MAXRO= MAXRO + MAXREM
               IGEO= NGEOME
            ELSE
               NGEOME= NGEOME-1
               IGEO= -1
            ENDIF
            KEYGEO(IB1)= IGEO
            CALL LCMSIX(IPGEOM,' ',2)
 40      CONTINUE
         IF( IPRT.GT.1 )THEN
            WRITE(IOUT,'(1H )')
            NXC= 1
            WRITE(IOUT,'(25H ===> PHYSICAL CELLS ARE:    /)')
            DO 42 IP= 1, (9+NGEOME)/10
               NXM= MIN( NGEOME, NXC+9 )
               WRITE(IOUT,'(1H ,10(I8.8,4X))')
     >         (IB1,IB1=NXC,NXM)
               WRITE(IOUT,'(1H ,30A4)') (CELLG(3*IB1-2),
     >         CELLG(3*IB1-1),CELLG(3*IB1),IB1=NXC,NXM)
               WRITE(IOUT,'(1H ,10(A7,5X))')
     >         ('TURN '//AXGTRN(KEYTRN(IB1)),IB1=NXC,NXM)
               NXC= NXC + 10
 42          CONTINUE
         ENDIF
*----
*  RESTORE *KEYTYP* AND *KEYGEO* VALUES
*----
         NTYP2= NGEOME
         DO 43 IB1= 1, NBLOCK
            KEYTYP(IB1)= KEYGEO(IB1)
 43      CONTINUE
         DO 44 IC1= 1, NTYP2
            IF( KEYGEO(IC1).NE.-1 ) KEYGEO(IC1)= IC1
 44      CONTINUE
*----
*  DELETE ALL VIRTUAL CELLS
*----
         NTYP= 0
         DO 45 ITYP= 1, NTYP2
            IGEO= KEYGEO(NTYP+1)
            IF( IGEO.EQ.-1 )THEN
               DO 46 IBLK= 1, NBLOCK
                  IOLTYP= KEYTYP(IBLK)
                  IF( IOLTYP.EQ.NTYP+1 )THEN
                     KEYTYP(IBLK)= 0
                  ELSEIF( IOLTYP.GT.NTYP+1 )THEN
                     KEYTYP(IBLK)= KEYTYP(IBLK)-1
                  ENDIF
 46            CONTINUE
               DO 47 JTYP= NTYP+2, NTYP2
                  CELLG(3*JTYP-5)= CELLG(3*JTYP-2)
                  CELLG(3*JTYP-4)= CELLG(3*JTYP-1)
                  CELLG(3*JTYP-3)= CELLG(3*JTYP  )
                  KEYTRN(JTYP-1)= KEYTRN(JTYP)
                  KEYGEO(JTYP-1)= KEYGEO(JTYP)
 47            CONTINUE
            ELSE
               NTYP= NTYP+1
            ENDIF
 45      CONTINUE
      ELSE
*----
*  NO CELL IN THE GEOMETRY
*----    
         GEOCV='            '
         READ(GEOCV,'(3A4)') CELLT(1),CELLT(2),CELLT(3)         
         IF( NTYPO.NE.1 )
     >      CALL XABORT(NAMSBR//': INVALID GEOMETRY TYPE '//GEONAM)
         NGEOME=    1
         NTYP= NTYPO
         READ(GEONAM,'(3A4)') (CELLG(3*NGEOME+ITC),ITC=-2,0)
         KEYGEO(1)= 1
         KEYTYP(1)= 1
         MRGCEL(1)= 1
         ITGEOM(1)= 1
         KEYTRN(1)= 1
         IF( LLSYM )THEN
*----
*  process x-x, y-y and z-z symmetry
*  1) Unfold geometry
*  2) Analyse symmetry
*----
           DO 330 IZ=MINGRI(3),1,-1
             IOF1=(IZ-1)*MINGRI(1)*MINGRI(2)
             IOF2=(IZ+MEDGRI(3)-1)*MAXGRI(1)*MAXGRI(2)
             DO 340 IY=MINGRI(2),1,-1
               DO 350 IX=MINGRI(1),1,-1
                 IOT2=IOF2+(IY+MEDGRI(2)-1)*MAXGRI(1)+IX+MEDGRI(1)
                 IOT1=IOF1+(IY-1)*MINGRI(1)+IX 
                 IF(IOT2 .NE. IOT1) THEN
                   IF(KEYTYP(IOT2) .NE. 0)THEN
                     CALL XABORT(NAMSBR//': PROBLEMS TO UNFOLD')
                   ELSE
                     KEYTYP(IOT2)=KEYTYP(IOT1)
                     KEYTYP(IOT1)= 0
                     ITGEOM(IOT2)=ITGEOM(IOT1)
                     ITGEOM(IOT1)= 0
                   ENDIF
                 ENDIF
  350          CONTINUE
  340        CONTINUE
  330      CONTINUE
           CALL AXGSYM( IPGEOM,   IPRT, NBLOCK,  NTYPO,   NXYZ,
     >                  GEONAM, LCLSYM,  
     >                  MINGRI, MAXGRI,  CELLT,
     >                  KEYTYP, ITGEOM)
         ENDIF
         IF( IPRT.GT.1 )THEN
            IB1=1
            WRITE(IOUT,'(A32/)') ' ===> REFERENCE GEOMETRY IS:    '
            WRITE(IOUT,'(1X,A8/1X,A12)') '00000001',GEONAM
*----
*  PRINTING ASSEMBLY MAP
*----
            CPLAN= BLANC
            NLINP= 3+(MAXGRI(2)+1)*((9+MAXGRI(1))/10+1)
            DO 430 IZ=1,MAXGRI(3)
              WRITE(IOUT,'(1H )')
              IF(NDIM.EQ.3)THEN
                 WRITE(CPLAN,'(4H (Z=,I3,1H))') IZ
              ENDIF
              WRITE(IOUT,'(/32H UNFOLD TYPE CELL MAP FOR PLANE ,A8)')
     >        CPLAN
              NXC= 1
              DO 435 IP  = 1, (9 + MAXGRI(1)) / 10
                NXM= MIN( MAXGRI(1), NXC+9 )
                WRITE(IOUT,'(1X,A8,1X,10(A4,I3,A4))')
     >             CPLAN, (' X= ',IR,' ROT',IR=NXC,NXM)
                NXC = NXC + 10
  435         CONTINUE
              WRITE(IOUT,'(1H )')
              DO 440 IY=1,MAXGRI(2)
                IOFF=((IZ-1)*MAXGRI(2)+(IY-1))*MAXGRI(1)
                NXC= 1
                DO 445 IP  = 1, (9 + MAXGRI(1)) / 10
                   NXM= MIN( MAXGRI(1), NXC+9 )
                   WRITE(IOUT,'(1X,A4,I3,2H=>,10(I7,1X,A2,1X))')
     >                 ' Y= ',IY,(KEYTYP(IOFF+IR),
     >                            AXGTRN(ITGEOM(IOFF+IR)),IR=NXC,NXM)
                   NXC = NXC + 10
  445           CONTINUE
                WRITE(IOUT,'(1H )')
  440         CONTINUE
  430       CONTINUE
         ENDIF
         NGEOME= 0
         DO 50 IB1= 1, NBLOCK
           IC1= KEYTYP(IB1)
           IT1= ITGEOM(IB1)
           IF(IC1 .LE. 0) THEN
             CALL XABORT(NAMSBR//': INVALID TYPE #')
           ENDIF
           NGEOME= NGEOME + 1
           READ(GEONAM,'(3A4)') (CELLG(3*NGEOME+ITC),ITC=-2,0)
           KEYTRN(NGEOME)= IT1
           IF( IB1.NE.1 )THEN
             DO 51 IB2= 1, IB1-1
               IC2= KEYTYP(IB2)
               IT2= ITGEOM(IB2)
               IF( IT1.NE.IT2 ) GO TO 51
               KEYGEO(IB1)= KEYGEO(IB2)
               NGEOME= NGEOME-1
               GO TO 50
 51          CONTINUE
           ENDIF
*          ANALYSE GEOMETRY
           CALL XELPRC(IPGEOM,GEONAM,NDIM,NNCYL,NNSUR,NNVOL,MAXREM)
           IF( NNVOL.NE.0 )THEN
             NSURO(NGEOME)= -NNSUR
             NVOLO(NGEOME)=  NNVOL
             IDLDIM(NGEOME)= NTOTCO
             NTOTCO= NTOTCO + NNCYL + 3
             MAXRO= MAXRO + MAXREM
             IGEO= NGEOME
           ELSE
             NGEOME= NGEOME-1
             IGEO= -1
           ENDIF
           KEYGEO(IB1)= IGEO
 50      CONTINUE
         IF( IPRT.GT.1 )THEN
           WRITE(IOUT,'(1H )')
           NXC= 1
           WRITE(IOUT,'(25H ===> PHYSICAL CELLS ARE:    /)')
           DO 52 IP= 1, (9+NGEOME)/10
             NXM= MIN( NGEOME, NXC+9 )
             WRITE(IOUT,'(1H ,10(I8.8,4X))')
     >       (IB1,IB1=NXC,NXM)
             WRITE(IOUT,'(1H ,30A4)') 
     >         ((CELLG(3*IB1+ITC),ITC=-2,0),IB1=NXC,NXM)
             WRITE(IOUT,'(1H ,10(A7,5X))')
     >       ('TURN '//AXGTRN(KEYTRN(IB1)),IB1=NXC,NXM)
             NXC= NXC + 10
 52        CONTINUE
         ENDIF
*----
*  RESTORE *KEYTYP* AND *KEYGEO* VALUES
*----
         NTYP= NGEOME
         DO 53 IB1= 1, NBLOCK
            KEYTYP(IB1)= KEYGEO(IB1)
 53      CONTINUE
         DO 54 IC1= 1, NTYP
            KEYGEO(IC1)= IC1
 54      CONTINUE
      ENDIF
      IF( IPRT.GT.1 )THEN
         WRITE(IOUT,'(/35H ONE TRACKING FILE TO BE ATTACHED    /'//
     >              '1H ,12X,14H UNDER NAME : ,A12 )') GEONAM
      ENDIF
*----
*  DEFINITION OF INDEX VALUES, TO LOOK AT THE DOMAIN
*----
       NGIDL= 0
      DO 210 IGEO= 1, NGEOME
          IF( NSURO(IGEO).GE.0 )
     >       CALL XABORT(NAMSBR//': GEOMETRY NOT FOUND')
          IF( NVOLO(IGEO).LE.0 )
     >       CALL XABORT(NAMSBR//': GEOMETRY NOT FOUND')
          IDLGEO(IGEO)=  NGIDL - NSURO(IGEO) + 1
           NGIDL= NVOLO(IGEO) + IDLGEO(IGEO)
  210 CONTINUE
      NTIDL = 0
      NPROB = 0
      DO 220 ITYP= 1, NTYP
            IGEO= KEYGEO(ITYP)
          IF( IGEO.LE.0 )
     >       CALL XABORT(NAMSBR//': BLOC NOT FOUND')
          IDLTYP(ITYP)= NTIDL - NSURO(IGEO) + 1
          IDLPRB= NPROB  + (1-NSURO(IGEO))*(2-NSURO(IGEO))/2
          NTIDL = NVOLO(IGEO) + IDLTYP(ITYP)
          NPROB = NVOLO(IGEO)*(NVOLO(IGEO)-2*NSURO(IGEO)+3)/2+IDLPRB
  220 CONTINUE
      NUNKO= 0
      DO 230 IBLK= 1, NBLOCK
            ITYP= KEYTYP(IBLK)
          IF( ITYP.LT.0 )
     >       CALL XABORT(NAMSBR//': CELL NOT FOUND')
          IF( ITYP.EQ.0 )GO TO 230
            IGEO= KEYGEO(ITYP)
          IDLBLK(IBLK)= NUNKO - NSURO(IGEO) + 1
          NUNKO= NVOLO(IGEO) + IDLBLK(IBLK)
  230 CONTINUE
      IF( IPRT.GT.10 )THEN
*----
*  PRINTING INDEX VECTORS
*----
          WRITE(IOUT,'(1H )')
          WRITE(IOUT,'(1H0,A6,4H -->,I8)') 'MAXRO', MAXRO
          WRITE(IOUT,'(1H0,A6,4H -->,I8)') 'NTOTCO', NTOTCO
          WRITE(IOUT,'(1H0,A6,4H -->,3(I8,1X))') 'LCLSYM', LCLSYM
          WRITE(IOUT,'(1H0,A8,4H    ,5(A8,2X))') '  GEOM #',
     >    '   NSURO', '   NVOLO', '  IDLGEO', '  IDLDIM', '  KEYTRN'
          DO 250 IGEO= 1, NGEOME
              WRITE(IOUT,'(1H ,I8,4H -->,5(I8,2X))') IGEO,
     >   NSURO(IGEO),NVOLO(IGEO),IDLGEO(IGEO),IDLDIM(IGEO),KEYTRN(IGEO)
  250     CONTINUE
          WRITE(IOUT,'(1H )')
          WRITE(IOUT,'(1H0,A8,4H    ,2(A8,2X))') '  BLOC #',
     >    '  KEYGEO', '  IDLTYP'
          DO 260 ITYP= 1, NTYP
              WRITE(IOUT,'(1H ,I8,4H -->,2(I8,2X))') ITYP,
     >        KEYGEO(ITYP), IDLTYP(ITYP)
  260     CONTINUE
          WRITE(IOUT,'(1H )')
          WRITE(IOUT,'(1H0,A8,4H    ,3(A8,2X))') '  CELL #',
     >    '  KEYTYP', '  ITGEOM',  '  IDLBLK'
          DO 270 IBLK= 1, NBLOCK
              WRITE(IOUT,'(1H ,I8,4H -->,3(I8,2X))') IBLK,
     >        KEYTYP(IBLK), ITGEOM(IBLK), IDLBLK(IBLK)
  270     CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(CELLT,ITGEOM)
*----
*  RETURN
*----
      RETURN 
      END
