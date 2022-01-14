*DECK XELETR
      SUBROUTINE XELETR(   IPRT,   NDIM, MAXGRI, NGEOME, NTOTCO, NTYPES,
     >                    NTIDL, NBLOCK,   NSUR,   NVOL, NTOTCL,  NUNKO,
     >                    NSURO,  NVOLO,  MINDO,  MAXDO, ICORDO, IDLDIM,
     >                   IDLGEO, KEYGEO, IDLTYP, KEYTYP, IDLBLK, KEYCYL,
     >                   RMESHO, IDLREM, INDEXO,  VOLSO, MATGEO, KEYINT,
     >                   MATTYP, REMESH, MINDIM, MAXDIM,  ICORD, VOLSUR,
     >                   KEYMRG,  INDEX, INCELL, MATALB,  NSURC,  NVOLC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Prepare tracking by producing the required numbering and recalculate
* mesh for an exact geometry treatment.
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
* IPRT    intermediate printing level for output.      
* NDIM    number of dimensions (2 or 3).                    
* MAXGRI  number of blocks in X/Y/Z directions.             
* NGEOME  number of geometries.                             
* NTOTCO  tot number of cylinders in all geometries.   
* NTYPES  number of cell types.  
* NTIDL   lenght of type numbering.                    
* NBLOCK  number of blocks.                                 
* NSUR    number of surfaces.                               
* NVOL    number of zones.                                  
* NTOTCL  tot number of cylinders in exact geometry.        
* NUNKO   old number of unknowns.                           
* NSURO   number of surfaces of each geometry.              
* NVOLO   number of zones of each geometry.                 
* MINDO   min index values for all axes (rect/cyl).    
* MAXDO   max index values for all axes (rect/cyl).    
* ICORDO  principal axes direction (X/Y/Z) for meshes. 
* IDLDIM  position of each geometry in cylinder numbering.  
* IDLGEO  position of each geometry in the             
*         geometry numbering scheme.        
* KEYGEO  geometric key for each type.                 
* IDLTYP  position of each type in numbering scheme.   
* KEYTYP  type key for each block.                     
* IDLBLK  position of each block in numbering scheme.  
* KEYCYL  index of cylinders by block.                 
* RMESHO  real mesh values (rect/cyl).                 
* IDLREM  position of mesh values per geometry.        
* INDEXO  index for search in 'rmesho'.                
* VOLSO   volumes and surfaces for each geometry.        
* MATGEO  material numbers corresponding to geometries.     
* KEYINT  interface key (giving the connected surface).
* MATTYP  material numbers for zones of every type.         
*
*Parameters: output
* REMESH  real mesh values (rect/cyl).                 
* MINDIM  min index values for all axes (rect/cyl).    
* MAXDIM  max index values for all axes (rect/cyl).    
* ICORD   principal axes direction (X/Y/Z) for meshes. 
* VOLSUR  volume-surface vector of exact geometry.     
* KEYMRG  merging vector of exact geometry.     
* INDEX   numbering of surfaces and zones.                    
* INCELL  block numbering.                                
* MATALB  material types.                              
* NSURC   number of compressed surfaces.                    
* NVOLC   number of compressed zones.                       
*
*-----------------------------------------------------------------------
*
      IMPLICIT           NONE
*
      INTEGER              IPRT,   NDIM, NGEOME, NTOTCO, NTYPES,
     >                    NTIDL, NBLOCK,   NSUR,   NVOL, NTOTCL,  NUNKO,
     >                    NSURC,  NVOLC
      INTEGER            MAXGRI(3),
     >                    MAXDO(NTOTCO),  MINDO(NTOTCO), ICORDO(NTOTCO),
     >                    NSURO(NGEOME),  NVOLO(NGEOME), IDLDIM(NGEOME),
     >                   IDLGEO(NGEOME), IDLREM(NGEOME), KEYGEO(NTYPES),
     >                   IDLTYP(NTYPES),
     >                   KEYTYP(NBLOCK), IDLBLK(NBLOCK), KEYCYL(NBLOCK),
     >                   INDEXO(4,*), MATGEO(*),
     >                   KEYINT(NUNKO), MATTYP(NTIDL),
     >                   MINDIM(NTOTCL), MAXDIM(NTOTCL), ICORD(NTOTCL),
     >                   INDEX(4,*), KEYMRG(*), MATALB(*), INCELL(*)
      REAL               RMESHO(*), REMESH(*), VOLSO(*), VOLSUR(*)
*
      INTEGER            ICUR(4)
      INTEGER            NUNK, IDLGE2, IG2, I4, N, ICX, IREM, I, J, K,
     >                   IBLK, ITYP, IGEO, IDLD, IDLR, MINABS, MAXABS,
     >                   J1, NTOTCX, MDMIN, NP1, NP2, IP1, IP2, IP3,
     >                   IOLD, ISU2, IVO2, ICREM, MINP1, MINP2, MINC,
     >                   ICYL, IDLTYX, IDLGEX, NO, NC, IMYG, IVSN,
     >                   IVS, IKREM
      REAL               RSTART, RMINUS, XP1, XP2
      CHARACTER          TEMESH(4)*8
      INTEGER            IOUT
      PARAMETER        ( IOUT=6 )
      INTEGER            NUMBLK, KL
      DATA        TEMESH / 'X', 'Y', 'Z', 'C' /
*
      NUMBLK(I,K)= I + IDLBLK(K)
*
*     INITIALIZE: NO INTERFACE & PUT INDEXES TO 0.
      NUNK = NVOL + 1 - NSUR
      IDLGE2=       1 - NSUR
      DO 5 IG2= 1, NUNK
         VOLSUR(IG2)= 0.0
         KEYMRG(IG2)= 0
         MATALB(IG2)= 0
         INCELL(IG2)= 0
         DO 4 I4= 1,4
            INDEX(I4,IG2)= 0
    4    CONTINUE
    5 CONTINUE
*
      IF( IPRT.GE.1 )THEN
          WRITE(IOUT,'(1H )')
          WRITE(IOUT,'(/24H  ====> GLOBAL MESHING   )')
      ENDIF
*
*     RECONSTRUCT CARTESIAN MESH
      J= 0
      ICUR(1)= 1
      ICUR(2)= 1
      ICUR(3)= 1
      IDLD=0
      DO 30 N= 1, 3
         RSTART= 0.0
         ICORD(N)= N
         MINDIM(N)= J+1
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
            IGEO= KEYGEO(ITYP)
            IDLD= IDLDIM(IGEO)
            IDLR= IDLREM(IGEO)
            MINABS= MINDO(IDLD+N)
            MAXABS= MAXDO(IDLD+N)
            RMINUS= RSTART - RMESHO(IDLR+MINABS)
            DO 10 IREM= MINABS, MAXABS-1
               J= J+1
               REMESH(J)= RMESHO(IDLR+IREM)+RMINUS
   10       CONTINUE
            ICUR(N)= 1
            RSTART= RMESHO(IDLR+MAXABS)+RMINUS
   20    CONTINUE
         J= J+1
         REMESH(J)= RSTART
         MAXDIM(N)= J
         IF( IPRT.GE.1.AND.N.LE.NDIM )THEN
            WRITE(IOUT,'(8X,A1,14H-COORDINATES: /(9X,5(1X,F13.6)))')
     >               TEMESH(N), (REMESH(J1),J1=MINDIM(N),MAXDIM(N))
         ENDIF
   30 CONTINUE
      NTOTCX= 3
*
*     RECONSTRUCT CYLINDRICAL MESH
      IF( NDIM.EQ.2 )THEN
         MDMIN= 3
      ELSE
         MDMIN= 1
      ENDIF
      DO 130 N= MDMIN, 3
         ICUR(N)= 1
         NP1= MOD(N  ,3) + 1
         NP2= MOD(N+1,3) + 1
*
*        (XP1,XP2) ARE COORDINATES AT BEGINNING OF BLOCK (IP1,IP2)
         XP2= 0.0
         DO 120 IP2= 1, MAXGRI(NP2)
         XP1= 0.0
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
            IF( ITYP.EQ.0 ) GO TO 105
            IGEO= KEYGEO(ITYP)
            IDLD= IDLDIM(IGEO)
            IDLR= IDLREM(IGEO)
            IF( IGEO.NE.NGEOME )THEN
               NC= IDLDIM(IGEO+1)-IDLD-3
            ELSE
               NC= NTOTCO-IDLD-3
            ENDIF
            IF( NC.EQ.1 )THEN
              IF( ICORDO(IDLD+4).EQ.N )THEN
                NTOTCX= NTOTCX+1
                IF( NTOTCX.GT.NTOTCL )
     >             CALL XABORT( '** XELETR: TOO MANY CYLINDERS' )
                MINP1 = MINDO(IDLD+NP1)
                MINP2 = MINDO(IDLD+NP2)
                MINC  = MINDO(IDLD+4)
                ICORD(NTOTCX)= ICORDO(IDLD+4)
*
*               RECENTER CYLINDERS
                REMESH(J+1)= RMESHO(IDLR+MINC-2)-RMESHO(IDLR+MINP1)+XP1
                REMESH(J+2)= RMESHO(IDLR+MINC-1)-RMESHO(IDLR+MINP2)+XP2
                J= J+2
                MINDIM(NTOTCX)= J+1
                DO 95 IREM= MINC, MAXDO(IDLD+4)
                   J=J+1
                   REMESH(J)= RMESHO(IDLR+IREM)
   95           CONTINUE
                MAXDIM(NTOTCX)= J
                IF( IPRT.GE.1 )THEN
                   WRITE(IOUT,'(13H        CELL(,I8,1H,,I8,1H,,I8,1H),
     >                     3H  (,A1,1H,,A1,10H)- CENTRE: ,
     >                     2H (,2(1X,F13.6),1H) )')
     >                      ICUR(1), ICUR(2), ICUR(3),
     >                      TEMESH(MOD(ICORD(NTOTCX)   ,3)+1),
     >                      TEMESH(MOD(ICORD(NTOTCX)+1,3)+1),
     >                      REMESH(MINDIM(NTOTCX)-2),
     >                      REMESH(MINDIM(NTOTCX)-1)
                   IF( NDIM.EQ.3 )THEN
                     WRITE(IOUT,'(24X,A1,8H-RADII: /(25X,5(1X,F13.6)))')
     >               TEMESH(ICORD(NTOTCX)),
     >               (SQRT(REMESH(J1)),J1=MINDIM(NTOTCX),MAXDIM(NTOTCX))
                   ELSE
                     WRITE(IOUT,'(26X,7HRADII: /(26X,5(1X,F13.6)))')
     >               (SQRT(REMESH(J1)),J1=MINDIM(NTOTCX),MAXDIM(NTOTCX))
                   ENDIF
                ENDIF
              ENDIF
            ENDIF
  105       CONTINUE
            IOLD= ICUR(NP2)
            ICUR(NP2)=   1
            IF( NDIM.EQ.2 )THEN
               IBLK= MAXGRI(1) * (ICUR(2) - 1) + ICUR(1)
            ELSE
               IBLK= MAXGRI(1)*(MAXGRI(2)*ICUR(3)+ICUR(2)-
     >                          MAXGRI(2))+ICUR(1)-MAXGRI(1)
            ENDIF
            ITYP= KEYTYP(IBLK)
            IGEO= KEYGEO(ITYP)
            IDLD= IDLDIM(IGEO)
            IDLR= IDLREM(IGEO)
            MINABS= MINDO(IDLD+NP1)
            MAXABS= MAXDO(IDLD+NP1)
            XP1= XP1 + (RMESHO(IDLR+MAXABS)-RMESHO(IDLR+MINABS))
            ICUR(NP2)= IOLD
  110    CONTINUE
            ICUR(NP1)=   1
            IF( NDIM.EQ.2 )THEN
               IBLK= MAXGRI(1) * (ICUR(2) - 1) + ICUR(1)
            ELSE
               IBLK= MAXGRI(1)*(MAXGRI(2)*ICUR(3)+ICUR(2)-
     >                          MAXGRI(2))+ICUR(1)-MAXGRI(1)
            ENDIF
            ITYP= KEYTYP(IBLK)
            IGEO= KEYGEO(ITYP)
            IDLD= IDLDIM(IGEO)
            IDLR= IDLREM(IGEO)
            MINABS= MINDO(IDLD+NP2)
            MAXABS= MAXDO(IDLD+NP2)
            XP2= XP2 + (RMESHO(IDLR+MAXABS)-RMESHO(IDLR+MINABS))
  120    CONTINUE
  130 CONTINUE
*
*     REESTABLISH INDEXING OF ALL UNKNOWNS
*     NOW, *ICUR()* IS THE INCREMENT FOR CARTESIAN CELL MESHING
      ISU2= 0
      IVO2= 0
      ICREM = 0
      ICUR(3)= 0
      DO 230 IP3= 1,MAXGRI(3)
      ICUR(2)= 0
      DO 220 IP2= 1,MAXGRI(2)
      ICUR(1)= 0
      DO 210 IP1= 1,MAXGRI(1)
         IF( NDIM.EQ.2 )THEN
            IBLK= MAXGRI(1)*(IP2-1)+IP1
         ELSE
            IBLK= MAXGRI(1)*(MAXGRI(2)*IP3+IP2-MAXGRI(2))+IP1-MAXGRI(1)
         ENDIF
         ITYP=   KEYTYP(IBLK)
         ICYL=   KEYCYL(IBLK)
         IGEO=   KEYGEO(ITYP)
         IDLTYX= IDLTYP(ITYP)
         IDLD=   IDLDIM(IGEO)
         IDLGEX= IDLGEO(IGEO)
         IKREM = ICREM
         DO 200 IVS= 1, NVOLO(IGEO)
            NO= NUMBLK(IVS, IBLK)
            IMYG=0
            IF( KEYINT(NO).NE.0 ) GO TO 200
            IMYG=MATGEO(IDLGEX+IVS)
            IVO2= IVO2 + 1
            IVSN= IVO2
            IF( IMYG.GE.0 )THEN
               IF( IVO2.GT.NVOL )
     >           CALL XABORT( '** XELETR: TOO MANY ZONES' )
               KEYMRG( IDLGE2+IVSN)= IMYG+ICREM
               VOLSUR( IDLGE2+IVSN)= VOLSO( IDLGEX+IVS)
               MATALB( IDLGE2+IVSN)= MATTYP( IDLTYX+IVS)
               INDEX(1,IDLGE2+IVSN)= INDEXO(1,IDLGEX+IVS)+ICUR(1)
     >                               + (MINDIM(1)-MINDO(IDLD+1))
               INDEX(2,IDLGE2+IVSN)= INDEXO(2,IDLGEX+IVS)+ICUR(2)
     >                               + (MINDIM(2)-MINDO(IDLD+2))
               INDEX(3,IDLGE2+IVSN)= INDEXO(3,IDLGEX+IVS)+ICUR(3)
     >                             + (MINDIM(3)-MINDO(IDLD+3))
               INDEX(4,IDLGE2+IVSN)= 0
               INCELL( IDLGE2+IVSN)= IBLK
               IF( ICYL.NE.0 )THEN
                  IF( INDEXO(4,IDLGEX+IVS).NE.MAXDO(IDLD+4) )THEN
*                    IF WE ARE INSIDE THE CYLINDER:
                     INDEX(4,IDLGE2+IVSN)=  INDEXO(4,IDLGEX+IVS)
     >                                   + (MINDIM(ICYL)-MINDO(IDLD+4))
                  ENDIF
               ENDIF
              IKREM=IKREM+1
            ELSE
               KEYMRG( IDLGE2+IVSN)= 0
               INCELL( IDLGE2+IVSN)= IBLK
            ENDIF
  200    CONTINUE
         ICREM=IKREM
         DO 400 IVS= -1,NSURO(IGEO),-1
            NO= NUMBLK(IVS, IBLK)
            IMYG=0
            IF( KEYINT(NO).NE.0 ) GO TO 400
            IMYG=MATGEO(IDLGEX+IVS)
            IF( IMYG.LT.0 )THEN
               ISU2= ISU2 - 1
               IVSN= ISU2
               IF( ISU2.LT. NSUR )
     >            CALL XABORT( '** XELETR: TOO MANY SURFACES' )
               KEYMRG( IDLGE2+IVSN)= IVSN
               VOLSUR( IDLGE2+IVSN)= VOLSO( IDLGEX+IVS)
               MATALB( IDLGE2+IVSN)= MATTYP( IDLTYX+IVS)
               INDEX(1,IDLGE2+IVSN)= INDEXO(1,IDLGEX+IVS)+ICUR(1)
     >                               + (MINDIM(1)-MINDO(IDLD+1))
               INDEX(2,IDLGE2+IVSN)= INDEXO(2,IDLGEX+IVS)+ICUR(2)
     >                               + (MINDIM(2)-MINDO(IDLD+2))
               INDEX(3,IDLGE2+IVSN)= INDEXO(3,IDLGEX+IVS)+ICUR(3)
     >                             + (MINDIM(3)-MINDO(IDLD+3))
               INDEX(4,IDLGE2+IVSN)= 0
               INCELL( IDLGE2+IVSN)= IBLK
               IF( ICYL.NE.0 )THEN
                  IF( INDEXO(4,IDLGEX+IVS).NE.MAXDO(IDLD+4) )THEN
*                    IF WE ARE INSIDE THE CYLINDER:
                       INDEX(4,IDLGE2+IVSN)=  INDEXO(4,IDLGEX+IVS)
     >                                   + (MINDIM(ICYL)-MINDO(IDLD+4))
                  ENDIF
               ENDIF
            ENDIF
  400    CONTINUE
         ICUR(1)= ICUR(1) + (MAXDO(IDLD+1)-MINDO(IDLD+1))
  210 CONTINUE
         ICUR(2)= ICUR(2) + (MAXDO(IDLD+2)-MINDO(IDLD+2))
  220 CONTINUE
         ICUR(3)= ICUR(3) + (MAXDO(IDLD+3)-MINDO(IDLD+3))
  230 CONTINUE
*----
*  REMOVE ZONES AND SURFACES WITH VANISHING VOLSUR
*---- 
      IVS=0
      DO 410 ISU2=NSUR,-1
        IF(VOLSUR(ISU2-NSUR+1) .GT. 0.0) THEN
          IVS=IVS+1
          VOLSUR(IVS)=VOLSUR(ISU2-NSUR+1)
          MATALB(IVS)=MATALB(ISU2-NSUR+1)
          KEYMRG(IVS)=KEYMRG(ISU2-NSUR+1)
          INCELL(IVS)=INCELL(ISU2-NSUR+1)
          DO 411 J1=1,4
            INDEX(J1,IVS)=INDEX(J1,ISU2-NSUR+1)
 411      CONTINUE
        ENDIF
 410  CONTINUE
      NSURC=-IVS
      IVS=IVS+1                 
      VOLSUR(IVS)=0.0
      MATALB(IVS)=0
      KEYMRG(IVS)=0
      INCELL(IVS)=0
      DO 420 J1=1,4
        INDEX(J1,IVS)=0
 420  CONTINUE
      DO 430 IVO2=1,NVOL
        IF(VOLSUR(IVO2-NSUR+1) .GT. 0.0) THEN
          IVS=IVS+1               
          VOLSUR(IVS)=VOLSUR(IVO2-NSUR+1)
          MATALB(IVS)=MATALB(IVO2-NSUR+1)
          KEYMRG(IVS)=KEYMRG(IVO2-NSUR+1)
          INCELL(IVS)=INCELL(IVO2-NSUR+1)
          DO 431 J1=1,4
            INDEX(J1,IVS)=INDEX(J1,IVO2-NSUR+1)
 431      CONTINUE
        ENDIF
 430  CONTINUE
      NVOLC=IVS+NSURC-1
      KL=1-NSURC
      IF( IPRT.GE.5 )THEN
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'(/13H RENUMBERING ,I8,13H VOLUMES AND ,'//
     >           'I8,10H SURFACES.)') NVOL,-NSUR
         WRITE(IOUT,'(20H CARTESIAN MESH      ,7HMINDIM=,3I8)')
     >                  (MINDIM(J1),J1=1,3)
         WRITE(IOUT,'(20X                     ,7HMAXDIM=,3I8)')
     >                  (MAXDIM(J1),J1=1,3)
         IF( NTOTCL.GT.3 )THEN
            DO 540 J1= 4, NTOTCL
               WRITE(IOUT,'(10H CYLINDER ,I8,6X     ,7HMINDIM=,15X,I8)')
     >                                 J1-3,MINDIM(J1)
               WRITE(IOUT,'(20X                     ,7HMAXDIM=,15X,I8)')
     >                                      MAXDIM(J1)
  540       CONTINUE
         ENDIF
         DO 550 IVS= NSURC, NVOLC
             IF( IVS.LT.0 )THEN
                IF( KEYMRG(IVS+KL) .EQ. 0 )THEN
                   WRITE(IOUT,'(8H KEYMRG(,I8,2H)=,I8,
     >             7H INDEX=,4I8,7H BLOCK=,I8,9H SURFACE=,F20.7,
     >             17H ABSENT FROM CELL)')
     >             IVS,KEYMRG(IVS+KL),(INDEX(J1,IVS+KL),J1=1,4),
     >             INCELL(IVS+KL),4.*VOLSUR(IVS+KL)
                ELSE
                   WRITE(IOUT,'(8H KEYMRG(,I8,2H)=,I8,
     >             7H INDEX=,4I8,7H BLOCK=,I8,9H SURFACE=,F20.7)')
     >             IVS,KEYMRG(IVS+KL),(INDEX(J1,IVS+KL),J1=1,4),
     >             INCELL(IVS+KL),4.*VOLSUR(IVS+KL)
                ENDIF
             ELSE
                IF( KEYMRG(IVS+KL) .EQ. 0 )THEN
                   WRITE(IOUT,'(8H KEYMRG(,I8,2H)=,I8,
     >             7H INDEX=,4I8,7H BLOCK=,I8,9H VOLUME= ,F20.7,
     >             17H ABSENT FROM CELL)')
     >             IVS,KEYMRG(IVS+KL),(INDEX(J1,IVS+KL),J1=1,4),
     >             INCELL(IVS+KL),VOLSUR(IVS+KL)
                ELSE
                   WRITE(IOUT,'(8H KEYMRG(,I8,2H)=,I8,
     >             7H INDEX=,4I8,7H BLOCK=,I8,9H VOLUME= ,F20.7)')
     >             IVS,KEYMRG(IVS+KL),(INDEX(J1,IVS+KL),J1=1,4),
     >             INCELL(IVS+KL),VOLSUR(IVS+KL)
                ENDIF
             ENDIF
  550    CONTINUE                              
      ENDIF
      RETURN
      END
