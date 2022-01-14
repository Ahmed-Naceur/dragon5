*DECK XELEDC
      SUBROUTINE XELEDC(   NDIM, MAXGRI, NGEOME, NTOTCO, NTYPES,
     >                   NBLOCK,  NUNKO,
     >                    NSURO,  NVOLO,  MINDO,  MAXDO,
     >                   ICORDO, IDLDIM, KEYGEO,
     >                   KEYTYP, IDLBLK, KEYINT,
     >                   NTOTCL,   MAXR,   NSUR,   NVOL, KEYCYL )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Associate all blocks of a problem to only one geometry and generate
* the 4 useful integer values that will describe the problem
* in its exact geometric description.    
*
*Copyright:
* Copyright (C) 1990 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* NDIM    number of dimensions.                             
* MAXGRI  number of grid cell in x/y/z directions.          
* NGEOME  number of geometries.                             
* NTOTCO  tot number of cylinders in all geometries.         
* NTYPES  number of types.                                  
* NBLOCK  number of blocks.                                 
* NUNKO   number of unknowns.                               
* NSURO   number of surfaces of each geometry.              
* NVOLO   number of zones of each geometry.                 
* MINDO   min index in the remesh array.               
* MAXDO   min index in the remesh array.               
* ICORDO  coordinate for remesh array.               
* IDLDIM  position of each geoemtry in cylinders numbering. 
* KEYGEO  geometric key for each type.                 
* KEYTYP  type key for each block.                     
* IDLBLK  position of each block in numbering scheme.  
* KEYINT  numbering of cell interfaces.                     
*
*Parameters: input
* NTOTCL  tot number of cylinders in exact geometry.        
* MAXR    lenght to stock real abscissae.              
* NSUR    number of surfaces of exact geometry (negative).  
* NVOL    number of zones of exact geometry.                
* KEYCYL  index of cylinders by block.                 
*
*-----------------------------------------------------------------------
*
      IMPLICIT     NONE
*
      INTEGER      NDIM, NGEOME, NTOTCO, NTYPES, NBLOCK, NUNKO,
     >             NTOTCL, MAXR, NSUR, NVOL
      INTEGER      MAXGRI(3),      NSURO(NTYPES),  NVOLO(NTYPES),
     >             MINDO(NTOTCO),  MAXDO(NTOTCO),  ICORDO(NTOTCO),
     >             IDLDIM(NTYPES), KEYGEO(NTYPES),
     >             KEYTYP(NBLOCK), IDLBLK(NBLOCK), KEYCYL(NBLOCK),
     >             KEYINT( NUNKO)
*
      INTEGER      ICUR(3), IBLK, N, ICX, ITYP, IGEO, IDLD, MDMIN,
     >             NP1, NP2, IP1, IP2, IP3, NC, NSUX, NVOX, IVX
      INTEGER      NUMBLK, I, K
*
      NUMBLK(I,K)= I + IDLBLK(K)
*
      DO 5 IBLK= 1, NBLOCK
         KEYCYL(IBLK)= 0
    5 CONTINUE
*
*     DETERMINE: NTOTCL & MAXR
*.1)  RECONSTRUCT CARTESIAN MESH
      MAXR= 0
      NTOTCL= 3
      ICUR(1)= 1
      ICUR(2)= 1
      ICUR(3)= 1
      DO 30 N= 1, 3
*
*        SCANNING CELLS ON THE AXIS #N
         DO 20 ICX= 1, MAXGRI(N)
            ICUR(N)= ICX
            IF( NDIM.EQ.2 )THEN
               IBLK= MAXGRI(1) * (ICUR(2) - 1) + ICUR(1)
            ELSE
               IBLK= MAXGRI(1)*(MAXGRI(2)*ICUR(3)+ICUR(2)-
     >                          MAXGRI(2))+ICUR(1)-MAXGRI(1)
            ENDIF
            ITYP= KEYTYP(IBLK)
            IF( ITYP.EQ.0 ) GO TO 20
            IGEO= KEYGEO(ITYP)
            IDLD= IDLDIM(IGEO)
            MAXR= MAXR + (MAXDO(IDLD+N)-MINDO(IDLD+N))
   20    CONTINUE
         ICUR(N)= 1
         MAXR= MAXR+1
   30 CONTINUE
*
*.2)  RECONSTRUCT INFORMATIONS FOR CYLINDRICAL MESH
      IF( NDIM.EQ.2 )THEN
         MDMIN= 3
      ELSE
         MDMIN= 1
      ENDIF
      DO 130 N= MDMIN, 3
         ICUR(N)= 1
         NP1= MOD(N  ,3) + 1
         NP2= MOD(N+1,3) + 1
         DO 120 IP2= 1, MAXGRI(NP2)
         DO 110 IP1= 1, MAXGRI(NP1)
            ICUR(NP1)= IP1
            ICUR(NP2)= IP2
            IF( NDIM.EQ.2 )THEN
               IBLK= MAXGRI(1) * (ICUR(2) - 1) + ICUR(1)
            ELSE
               IBLK= MAXGRI(1)*(MAXGRI(2)*ICUR(3)+ICUR(2)-
     >                          MAXGRI(2))+ICUR(1)-MAXGRI(1)
            ENDIF
            ITYP= KEYTYP(IBLK)
            IF( ITYP.EQ.0 ) GO TO 110
            IGEO= KEYGEO(ITYP)
            IDLD= IDLDIM(IGEO)
            IF( IGEO.NE.NGEOME )THEN
               NC= IDLDIM(IGEO+1)-IDLD-3
            ELSE
               NC= NTOTCO-IDLD-3
            ENDIF
            IF( NC.EQ.1 )THEN
               IF( ICORDO(IDLD+4).EQ.N )THEN
                  NTOTCL= NTOTCL+1
                  MAXR= MAXR + 3 + (MAXDO(IDLD+4)-MINDO(IDLD+4))
                  DO 105 IP3= 1, MAXGRI(N)
                     ICUR(N)= IP3
                     IF( NDIM.EQ.2 )THEN
                        IBLK= MAXGRI(1) * (ICUR(2) - 1) + ICUR(1)
                     ELSE
                        IBLK= MAXGRI(1)*(MAXGRI(2)*ICUR(3)+ICUR(2)-
     >                                   MAXGRI(2))+ICUR(1)-MAXGRI(1)
                     ENDIF
                     KEYCYL(IBLK)= NTOTCL
  105             CONTINUE
                  ICUR(N)= 1
               ENDIF
            ENDIF
  110    CONTINUE
  120    CONTINUE
  130 CONTINUE
*
*     DETERMINE: NSUR & NVOL
      NSUR= 0
      NVOL= 0
      DO 230 IBLK= 1,NBLOCK
         ITYP= KEYTYP(IBLK)
         IF( ITYP.EQ.0 ) THEN
            CALL XABORT( '*** XELEDC: EXACT VOID CELL NOT ALLOWED')
         ENDIF
         IGEO= KEYGEO(ITYP)
         NSUX= NSURO(IGEO)
         NVOX= NVOLO(IGEO)
         DO 220 IVX= NSUX, NVOX
            IF( IVX.LT.0 )THEN
               IF( KEYINT(NUMBLK(IVX,IBLK)).EQ.0 ) NSUR= NSUR-1
            ELSEIF( IVX.GT.0 )THEN
               NVOL= NVOL + 1
            ENDIF
  220    CONTINUE
  230 CONTINUE
*
      RETURN
      END
