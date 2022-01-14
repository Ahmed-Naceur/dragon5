*DECK XELTRP
      SUBROUTINE XELTRP( IPGEOM,  NGIDL,   NDIM, NGEOME, L1CELL,
     >                   NTOTCO, NEXTGE,  MAXRO,   IPRT,  CELLG,
     >                    NSURO,  NVOLO, IDLDIM, IDLGEO, KEYTRN,
     >                    MAXDO,  MINDO, ICORDO, RMESHO, IDLREM,
     >                   INDEXO,  VOLSO, MATGEO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Prepare tracking by producing the required numbering and calculate 
* volumes and surfaces. 
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
* NGIDL   lenght of geometric numbering.              
* NDIM    number of dimensions (2 or 3).                   
* NGEOME  number of geometries.                            
* L1CELL  to indicate if there is just 1 cell.        
* NEXTGE  rectangular(0)/circular(1) boundary.        
* NTOTCO  tot number of cylinders in all geometries.  
* MAXRO   max number of real mesh values in RMESHO. 
* IPRT    intermediate printing level for output.     
* CELLG   to keep geomety names.                      
* NSURO   number of surfaces of each geometry.             
* NVOLO   number of zones of each geometry.                
* IDLDIM  position of each geometry in cylinder numbering. 
* IDLGEO  position of each geometry in the            
*         geometry numbering scheme.       
* KEYTRN  turn number of each geometry.                    
*
*Parameters: input
* MAXDO   max index values for all axes (rect/cyl).    
* MINDO   min index values for all axes (rect/cyl).    
* ICORDO  principal axes direction (X/Y/Z) for meshes. 
* RMESHO  real mesh values (rect/cyl).                 
* IDLREM  position of mesh values per geometry.        
* INDEXO  index for search in RMESHO.                
* VOLSO   volumes & surfaces for each geometry.        
* MATGEO  material numbers corresponding to geometries.     
*
*-----------------------------------------------------------------------
*
      USE               GANLIB
      IMPLICIT          NONE
*
      TYPE(C_PTR)       IPGEOM 
      INTEGER           NGIDL, NDIM, NGEOME, NTOTCO, NEXTGE, MAXRO, IPRT
      INTEGER           MAXDO(NTOTCO), MINDO(NTOTCO),   ICORDO(NTOTCO),
     >                  MATGEO(NGIDL), CELLG(3*NGEOME),
     >                  NSURO(NGEOME),  NVOLO(NGEOME), IDLDIM(NGEOME),
     >                  IDLGEO(NGEOME), IDLREM(NGEOME), KEYTRN(NGEOME),
     >                  INDEXO(4,NGIDL)
      REAL              RMESHO(MAXRO), VOLSO(NGIDL)
*
      INTEGER           NSTATE, IOUT, MAXTUR
      PARAMETER       ( NSTATE=40, IOUT=6, MAXTUR=12 )
      INTEGER           ISTATE(NSTATE)
      INTEGER           NTOTRM, NGEO, NTC, ITURN, NC, NCPC, NVSP1,
     >                  NO, NSYM, MAXC, KELRNG, KELMRG, KELSYM
      LOGICAL           L1CELL
      CHARACTER         CNAMEG*12, CTURN(2*MAXTUR)*2
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEYSYM
*----
*  DATA STATEMENTS
*----
      DATA       CTURN / ' A',' B',' C',' D',' E',' F',' G',' H',
     >                   ' I',' J',' K',' L',
     >                   '-A','-B','-C','-D','-E','-F','-G','-H',
     >                   '-I','-J','-K','-L' /
*----
*  SCRATCH STORAGE ALLOCATION
*   KEYSYM: symmetry key giving the symmetric surface
*----
      ALLOCATE(KEYSYM(NGIDL))
*
*     LOOP OVER ALL GEOMETRIES
      NTOTRM= 0
      DO 90 NGEO= 1, NGEOME
         NTC= IDLDIM(NGEO)+1
         ITURN= KEYTRN(NGEO)
         WRITE( CNAMEG( 1: 4),'(A4)') CELLG(3*NGEO-2)
         WRITE( CNAMEG( 5: 8),'(A4)') CELLG(3*NGEO-1)
         WRITE( CNAMEG( 9:12),'(A4)') CELLG(3*NGEO  )
         IF( .NOT.L1CELL ) CALL LCMSIX(IPGEOM, CNAMEG, 1)
         CALL XDISET(ISTATE,NSTATE,0)
         CALL LCMGET(IPGEOM, 'STATE-VECTOR', ISTATE)
         IF( ISTATE(1).GE.20.OR.ISTATE(1).EQ.3.OR.ISTATE(1).EQ.6 )THEN
            NC= 1
         ELSE
            NC= 0
         ENDIF
         IF( IPRT.GT.1 )THEN
            WRITE(IOUT,'(1H )')
            IF    ( NC.EQ.0 )THEN
              WRITE(IOUT,'(/27H NUMBERING PHYSICAL CELL # ,I8/6H  >>> ,
     >                     A12,6H /ROT ,A2,13H GEOMETRY <<<,
     >                     13H    (WITH  NO,11H CYLINDER ) )')
     >                                 NGEO,        CNAMEG,CTURN(ITURN)
            ELSEIF( NC.EQ.1 )THEN
              WRITE(IOUT,'(/27H NUMBERING PHYSICAL CELL # ,I8/6H  >>> ,
     >                     A12,6H /ROT ,A2,13H GEOMETRY <<<,
     >                     13H    (WITH ONE,11H CYLINDER ) )')
     >                                 NGEO,        CNAMEG,CTURN(ITURN)
            ELSE
              WRITE(IOUT,'(/27H NUMBERING PHYSICAL CELL # ,I8/6H  >>> ,
     >                     A12,6H /ROT ,A2,13H GEOMETRY <<<,
     >                     10H    (WITH ,I3,11H CYLINDERS) )')
     >                               NGEO, CNAMEG, CTURN(ITURN), NC
            ENDIF
         ENDIF
         NCPC  = NC + 3
         NVSP1 = NVOLO(NGEO) - NSURO(NGEO) + 1
*
*        LOOKING TO THE GEOMETRY
         CALL XELGRD( IPGEOM, IPRT, NDIM, NEXTGE, ITURN,
     >                MAXRO-NTOTRM, MAXC, RMESHO(NTOTRM+1),
     >                MINDO(NTC), MAXDO(NTC), ICORDO(NTC))
*
*        RENUMBER
         NO=   KELRNG(IPRT, NDIM, NEXTGE, NCPC,
     >                MINDO(NTC), MAXDO(NTC), ICORDO(NTC),
     >                NSURO(NGEO), NVOLO(NGEO), IDLGEO(NGEO),
     >                MAXC, RMESHO(NTOTRM+1), MATGEO, VOLSO, INDEXO)
*
*        MERGE
         NO= KELMRG(IPGEOM,NSURO(NGEO),NVOLO(NGEO),IDLGEO(NGEO),MATGEO)
         IF( NO.NE.NVSP1 )THEN
            IF( IPRT.GT.1 )THEN
               WRITE(IOUT,'(1H )')
               WRITE(IOUT,'(22H     MERGE INTO   >>> ,I8,
     >                  13H  ZONES   <<<)')
     >                       NO+NSURO(NGEO)-1
            ENDIF
         ENDIF
*
*        ESTABLISH NECESSARY SYMMETRIES
         NSYM= KELSYM( IPRT, NDIM, MAXDO(NTC), NSURO(NGEO), NVOLO(NGEO),
     >                 IDLGEO(NGEO), INDEXO, MATGEO,KEYSYM)
*
*        COMPUTE VOLUMES
         CALL XELVOL( IPRT, NDIM, NEXTGE, NCPC,
     >                MINDO(NTC), MAXDO(NTC), ICORDO(NTC),
     >                NSURO(NGEO), NVOLO(NGEO), IDLGEO(NGEO),INDEXO,
     >                MAXC, RMESHO(NTOTRM+1), MATGEO, VOLSO )
         IDLREM(NGEO)= NTOTRM
         NTOTRM= NTOTRM + MAXC
         IF( .NOT.L1CELL ) CALL LCMSIX(IPGEOM, ' ', 2 )
   90 CONTINUE
      IF( NTOTRM.GT.MAXRO )THEN
         CALL XABORT( 'XELTRP : INCREASE MAXREM => SEE DEVELOPPER')
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(KEYSYM)
*
      RETURN
      END
