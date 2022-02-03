*DECK XELBIN
      SUBROUTINE XELBIN( IPGEOM,   NDIM, NGEOME, L1CELL, NTYPES,  NGIDL,
     >                    NTIDL, NBLOCK, MAXGRI,  NUNKO,   IPRT,  CELLG,
     >                    NSURO,  NVOLO, IDLGEO, MATGEO, KEYGEO, IDLTYP,
     >                   IDLBLK, KEYTYP, MATTYP, KEYINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Identify every zone of every type  to its material and 
* interface all internal surfaces for cells present in the supercell. 
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
* NDIM    number of dimensions (2 or 3).                    
* NGEOME  number of geometries.                             
* L1CELL  .true. if only one cell.                     
* NTYPES  number of types.                                  
* NGIDL   lenght of geometric numbering.               
* NTIDL   lenght of type numbering.                    
* NBLOCK  number of blocks.                                 
* MAXGRI  number of cells along each axis.                  
* NUNKO   old number of unknowns.                           
* IPRT    intermediate printing level for output.      
* CELLG   to keep geometry  names.                     
* NSURO   number of surfaces of each geometry.              
* NVOLO   number of zones of each geometry.                 
* IDLGEO  position of each geometry in the             
*         geometry numbering scheme.        
* MATGEO  material numbers corresponding to geometries.     
* KEYGEO  geometric key for each type.                 
* IDLTYP  position of each type in numbering scheme.   
* IDLBLK  position of each block in numbering scheme.  
* KEYTYP  type key for each block.                     
*
*Parameters: output
* MATTYP  material numbers for zones of every type.         
* KEYINT  interface key (giving the connected surface).
*
*-----------------------------------------------------------------------
*
      USE                GANLIB
      IMPLICIT           NONE
*
      TYPE(C_PTR)        IPGEOM 
      INTEGER            NDIM, NGEOME, NTYPES, NGIDL, NTIDL, NBLOCK,
     >                   NUNKO, IPRT
      INTEGER             NSURO(NGEOME),  NVOLO(NGEOME), IDLGEO(NGEOME),
     >                   MATGEO( NGIDL), KEYGEO(NTYPES), IDLTYP(NTYPES),
     >                   MATTYP( NTIDL), KEYTYP(NBLOCK), IDLBLK(NBLOCK),
     >                   KEYINT(NUNKO ), MAXGRI(NDIM)  , CELLG(3*NTYPES)
*
      INTEGER            ILO(3,2), NO(2), KTYP(2),
     >                     KMAT(2),    KSUR(2), KABSO(2), KSID(2),
     >                   ICOORD(3),   NCODE(6)
      CHARACTER          GEOCEL*12, TEDATA*12, TEMESH(4)*7
      LOGICAL            SWKILL, L1CELL, LL1, LL2
      INTEGER            NSTATE, IOUT, MAXSPL
      PARAMETER        ( NSTATE=40, IOUT=6, MAXSPL=100 )
      INTEGER            ISTATE(NSTATE),ISPLT(MAXSPL)
      INTEGER            NUMGEO, NUMTYP, NUMBLK, I, K
      INTEGER            NBMD
      INTEGER            IMYT, IUNK, ITYP, IMYG, IGEO, NSUX, NVOX, ICYL,
     >                   ICX, ICY, ICZ, LR, LX, LY, LZ, KOLD, ITYPG,
     >                   ISUR, IX, IY, IZ, IOFF, KNEW, ILEN, ITYLCM,
     >                   ISX, ISY, ISZ, ISR, KIOFX, KIOFY, KIOFZ,
     >                   J0, J1, J2, JC, JR, IP0, IP1, IP2, N, NP1, NP2,
     >                   K0, K1, K2, K3, KR, IBLK, ISUX
      EQUIVALENCE      ( ICOORD(1),LX ),(ICOORD(2),LY),(ICOORD(3),LZ )
      DATA        TEMESH / 'X', 'Y', 'Z', 'R'/
*
      NUMGEO(I,K)= I + IDLGEO(K)
      NUMTYP(I,K)= I + IDLTYP(K)
      NUMBLK(I,K)= I + IDLBLK(K)
*
      SWKILL= .FALSE.
      LL1= .FALSE.
      LL2= .FALSE.
      DO 10 IMYT= 1,  NTIDL
          MATTYP(IMYT)=  0
   10 CONTINUE
      DO 20 IUNK= 1, NUNKO
          KEYINT(IUNK)=  0
   20 CONTINUE
      DO 40 ITYP= 1, NTYPES
         IGEO  = KEYGEO( ITYP   )
         NVOX  = NVOLO( IGEO )
         NSUX  = NSURO( IGEO ) 
         IF( .NOT.L1CELL )THEN
            WRITE(GEOCEL( 1: 4), '(A4)') CELLG(3*ITYP-2)
            WRITE(GEOCEL( 5: 8), '(A4)') CELLG(3*ITYP-1)
            WRITE(GEOCEL( 9:12), '(A4)') CELLG(3*ITYP  )
            CALL LCMSIX(IPGEOM, GEOCEL, 1)
         ELSE
            CALL LCMGET(IPGEOM,'NCODE', NCODE)
            LL1=((NCODE(2).EQ.3).AND.(NCODE(3).EQ.3))
            LL2=((NCODE(1).EQ.3).AND.(NCODE(4).EQ.3))
         ENDIF
         CALL XDISET(ISTATE,NSTATE,0)
         CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
         ITYPG= ISTATE(1)
         IF( ITYPG.EQ.20) THEN
*           FOR *CARCEL* GEOMETRIES
            ICYL= 1
            ICX=  1
            ICY=  2
            ICZ=  3
         ELSEIF(ITYPG.EQ.3.OR.ITYPG.EQ.6 )THEN
*           FOR *CARCEL*, *TUBE* OR *TUBEZ* GEOMETRIES
            ICYL= 1
            ICX=  1
            ICY=  2
            ICZ=  3
            IF( LL1.OR.LL2 )THEN
               CALL XABORT( 'XELBIN: DIAGONAL SYMETRIES NOT POSSIBLE')
            ENDIF
         ELSEIF( ITYPG.GT.20 )THEN
*           FOR *CARCELX*, *CARCELY* OR *CARCELZ*
            ICYL= 1
            ICZ= ITYPG-20
            ICX= MOD(ICZ  , 3) + 1
            ICY= MOD(ICZ+1, 3) + 1
         ELSE
*           FOR *CAR2D* OR *CAR3D*
            ICYL= 0
            ICX=  1
            ICY=  2
            ICZ=  3
         ENDIF
         LR=    ISTATE(2)
         LX=    MAX(1,ISTATE(3))
         LY=    MAX(1,ISTATE(4))
         LZ=    MAX(1,ISTATE(5))
         KOLD=  ISTATE(6)
         DO 30 ISUR= NSUX, -1
            MATTYP(NUMTYP(ISUR,ITYP))= MATGEO(NUMGEO(ISUR,IGEO))
   30    CONTINUE
*
*        GET MIXTURE NUMBERS
         CALL LCMLEN(IPGEOM, 'MIX', ILEN, ITYLCM)
         IF( ILEN.NE.KOLD )THEN
            WRITE(IOUT,*) 'LENGHT(MIX)=  ',ILEN
            WRITE(IOUT,*) '# OF VOLUMES= ',KOLD
            CALL LCMLIB(IPGEOM)
            CALL XABORT( 'XELBIN: INVALID NUMBER OF MIXTURES')
         ENDIF
         CALL LCMGET(IPGEOM,'MIX',MATTYP(NUMTYP(1,ITYP)))
*
*        IN THE CASE OF DIAGONAL SYMMETRY IN 'ONE-CELL'
*        CAR2D AND CAR3D  GEOMETRY UNFOLD MIXTURES 
*
      IF(ITYPG .LT. 20) THEN
         K3=ISTATE(6)
         NBMD=(LZ*LY*(LX+1))/2
         IF(K3 .EQ. NBMD) THEN
*----
* MIXTURE ENTERED IN DIAGONAL FORM
*----
         IF( LL1 )THEN
            DO 72 IZ=LZ,1,-1
            IOFF=(IZ-1)*LX*LY
            DO 71 IY=LY,1,-1
            DO 60 IX=LX,IY+1,-1
            MATTYP(NUMTYP(IOFF+(IY-1)*LX+IX,ITYP))=
     >                    MATTYP(NUMTYP(IOFF+(IX-1)*LY+IY,ITYP))
   60       CONTINUE
            DO 70 IX=IY,1,-1
            MATTYP(NUMTYP(IOFF+(IY-1)*LX+IX,ITYP))=
     >                    MATTYP(NUMTYP(K3,ITYP))
            K3=K3-1
   70       CONTINUE
   71       CONTINUE
   72       CONTINUE
            KOLD= LX*LY*LZ
         ELSEIF( LL2 )THEN
            DO 82 IZ=LZ,1,-1
            IOFF=(IZ-1)*LX*LY
            DO 81 IY=LY,1,-1
            DO 80 IX=LX,IY,-1
            MATTYP(NUMTYP(IOFF+(IY-1)*LX+IX,ITYP))=
     >                    MATTYP(NUMTYP(K3,ITYP))
            K3=K3-1
   80       CONTINUE
   81       CONTINUE
   82       CONTINUE
            DO 92 IZ=1,LZ
            IOFF=(IZ-1)*LX*LY
            DO 91 IY=1,LY
            DO 90 IX=1,IY-1
            MATTYP(NUMTYP(IOFF+(IY-1)*LX+IX,ITYP))=
     >               MATTYP(NUMTYP(IOFF+(IX-1)*LY+IY,ITYP))
   90       CONTINUE
   91       CONTINUE
   92       CONTINUE
            KOLD= LX*LY*LZ
         ENDIF
         ENDIF
      ENDIF
*
*        FOR THE PARTICULAR CASE OF *TUBE* OR *TUBEZ* GEOMETRIES
         IF( ITYPG.EQ.3.OR.ITYPG.EQ.6 )THEN
            DO 39 IZ= 1, LZ
               MATTYP(NUMTYP(KOLD+IZ,ITYP))= -2
   39       CONTINUE
            KOLD= KOLD+LZ
         ENDIF
*
*        FILL UP MATTYP ACCORDING TO SPLITTING VALUES.
         KNEW= NVOX
         ISR= 0
         ISX= 0
         ISY= 0
         ISZ= 0
         DO 308 K0= ICOORD(ICZ),1,-1
            KIOFZ= KOLD
            TEDATA= 'SPLIT'//TEMESH(ICZ)
            CALL LCMLEN(IPGEOM,TEDATA,ILEN,ITYLCM)
            IF( ILEN.GT.MAXSPL )THEN
               CALL XABORT('XELBIN: SPLIT OVERFLOW ('//TEDATA//')')
            ELSEIF( ILEN.EQ.0 )THEN
               ISZ= 1
            ELSE
               CALL LCMGET(IPGEOM,TEDATA,ISPLT)
               ISZ= ISPLT(K0)
            ENDIF
         DO 307 J0=ISZ,1,-1
            KOLD= KIOFZ
         DO 306 K1= ICOORD(ICY),1,-1
            KIOFY= KOLD
            TEDATA= 'SPLIT'//TEMESH(ICY)
            CALL LCMLEN(IPGEOM,TEDATA,ILEN,ITYLCM)
            IF( ILEN.GT.MAXSPL )THEN
               CALL XABORT('XELBIN: SPLIT OVERFLOW ('//TEDATA//')')
            ELSEIF( ILEN.EQ.0 )THEN
               ISY= 1
            ELSE
               CALL LCMGET(IPGEOM,TEDATA,ISPLT)
               ISY= ISPLT(K1)
            ENDIF
         DO 305 J1=ISY,1,-1
            KOLD= KIOFY
         DO 304 K2= ICOORD(ICX),1,-1
            KIOFX= KOLD
            TEDATA= 'SPLIT'//TEMESH(ICX)
            CALL LCMLEN(IPGEOM,TEDATA,ILEN,ITYLCM)
            IF( ILEN.GT.MAXSPL )THEN
               CALL XABORT('XELBIN: SPLIT OVERFLOW ('//TEDATA//')')
            ELSEIF( ILEN.EQ.0 )THEN
               ISX= 1
            ELSE
               CALL LCMGET(IPGEOM,TEDATA,ISPLT)
               ISX= ISPLT(K2)
            ENDIF
         DO 303 J2=ISX,1,-1
            KOLD= KIOFX
*           FOR RECTANGULAR OUTER REGIONS.
            IMYT= MATTYP(NUMTYP(KOLD,ITYP))
            MATTYP(NUMTYP(KNEW,ITYP))= IMYT
            KNEW= KNEW-1
            KOLD= KOLD-1
            IF( ICYL.EQ.1 )THEN
*              FOR CYLINDRICAL INNER REGIONS.
               DO 302 KR= LR,1,-1
                  TEDATA= 'SPLIT'//TEMESH(4)
                  CALL LCMLEN(IPGEOM,TEDATA,ILEN,ITYLCM)
                  IF( ILEN.GT.MAXSPL )THEN
                    CALL XABORT('XELBIN: SPLIT OVERFLOW ('//TEDATA//')')
                  ELSEIF( ILEN.EQ.0 )THEN
                     ISR= 1
                  ELSE
                     CALL LCMGET(IPGEOM,TEDATA,ISPLT)
                     ISR= ABS(ISPLT(KR))
                  ENDIF
                  IMYT= MATTYP(NUMTYP(KOLD,ITYP))
               DO 301 JR=ISR,1,-1
                  MATTYP(NUMTYP(KNEW,ITYP))= IMYT
                  KNEW= KNEW-1
  301          CONTINUE
                  KOLD= KOLD-1
  302          CONTINUE
            ENDIF
  303    CONTINUE
  304    CONTINUE
  305    CONTINUE
  306    CONTINUE
  307    CONTINUE
  308    CONTINUE
         IF( KNEW.NE.0 )THEN
            WRITE(IOUT,*) 'XELBIN: KNEW.NE.0 = PROBLEM WITH SPLITTING'
            SWKILL= .TRUE.
         ENDIF
         IF( KOLD.NE.0 )THEN
            WRITE(IOUT,*) 'XELBIN: KOLD.NE.0 = PROBLEM WITH SPLITTING'
            SWKILL= .TRUE.
         ENDIF
*
         IF( .NOT.L1CELL ) CALL LCMSIX(IPGEOM, ' ',    2)
   40 CONTINUE
*
*     RECOMPOSE INTERNAL SURFACES COUPLING (INTERFACES)
*        THIS ASSUMES THAT AN ORDERING OF SURFACES IS DONE
*        BECAUSE:  SIDE-BY-SIDE INTERFACES
*                  ARE SUPPOSED IN INCREASING POSITION.
      DO 220 N= 1, NDIM
*
*        DEFINITION OF THE SIDE NUMBER TO COUPLE.
         KSID(1)=  -2*N
         KSID(2)= (-2*N) + 1
         NP1   = MOD(N  ,NDIM) + 1
         IF( NDIM.EQ.3 )THEN
            NP2   = MOD(N+1,NDIM) + 1
            DO 112 IP1= 1, MAXGRI(NP1)
               ILO(NP1,1)= IP1
               ILO(NP1,2)= IP1
            DO 111 IP2= 1, MAXGRI(NP2)
               ILO(NP2,1)= IP2
               ILO(NP2,2)= IP2
            DO 110 IP0= 1, MAXGRI(N)-1
               ILO(N  ,1)= IP0
               ILO(N  ,2)= IP0 + 1
               DO 100  JC= 1, 2
                  NO(JC)= MAXGRI(1)*(MAXGRI(2)*ILO(3,JC)+ILO(2,JC)-
     >                               MAXGRI(2))+ILO(1,JC)-MAXGRI(1)
                  KTYP(JC)= KEYTYP( NO(JC) )
                  IF( KTYP(JC).EQ.0 ) GO TO 110
                  IGEO  = KEYGEO( KTYP(JC) )
*                 SEARCH FROM THE END
                  KSUR(JC)= NSURO(IGEO)
                  KMAT(JC)= MATTYP( NUMTYP(KSUR(JC),KTYP(JC)) )
  100          CONTINUE
*
*              ORDERING INTERFACING OF THE TWO BLOCKS.
  101          CONTINUE
                  IF( KMAT(1).EQ.KSID(1).AND.KMAT(2).EQ.KSID(2) )THEN
                     IF( KSUR(1).EQ.0 .OR. KSUR(2).EQ.0 ) GO TO 109
                     KABSO(1)= NUMBLK( KSUR(1),NO(1) )
                     KABSO(2)= NUMBLK( KSUR(2),NO(2) )
                     KEYINT( KABSO(1) )= KABSO(2)
                     KEYINT( KABSO(2) )= KABSO(1)
                     KSUR(1)= KSUR(1)+1
                     KSUR(2)= KSUR(2)+1
                  ELSE
                     IF( KMAT(1).NE.KSID(1) ) KSUR(1)= KSUR(1)+1
                     IF( KMAT(2).NE.KSID(2) ) KSUR(2)= KSUR(2)+1
                  ENDIF
                  IF( KSUR(1).NE.0 )THEN
                     KMAT(1)= MATTYP( NUMTYP(KSUR(1),KTYP(1)) )
                  ELSE
                     KMAT(1)= KSID(1)
                  ENDIF
                  IF( KSUR(2).NE.0 )THEN
                     KMAT(2)= MATTYP( NUMTYP(KSUR(2),KTYP(2)) )
                  ELSE
                     KMAT(2)= KSID(2)
                  ENDIF
               GO TO 101
  109          IF( KSUR(1).NE.0 .OR. KSUR(2).NE.0 )THEN
                  WRITE(IOUT,'(1H ,I8,4H OF ,I8,5H <=> ,I8,4H OF ,I8)')
     >                          KSUR(1),  NO(1),     KSUR(2),  NO(2)
                  SWKILL=.TRUE.
               ENDIF
  110       CONTINUE
  111       CONTINUE
  112       CONTINUE
         ELSEIF( NDIM.EQ.2 )THEN
            DO 215 IP1= 1, MAXGRI(NP1)
               ILO(NP1,1)= IP1
               ILO(NP1,2)= IP1
            DO 210 IP0= 1, MAXGRI(N)-1
               ILO(N  ,1)= IP0
               ILO(N  ,2)= IP0 + 1
               DO 200  JC= 1, 2
                  NO(JC)= MAXGRI(1) * (ILO(2,JC) - 1) + ILO(1,JC)
                  KTYP(JC)= KEYTYP( NO(JC) )
                  IF( KTYP(JC).EQ.0 ) GO TO 210
                  IGEO  = KEYGEO( KTYP(JC) )
*                 SEARCH FROM THE END
                  KSUR(JC)= NSURO(IGEO)
                  KMAT(JC)= MATTYP( NUMTYP(KSUR(JC),KTYP(JC)) ) 
  200          CONTINUE
*
*              ORDERING INTERFACING OF THE TWO BLOCKS.
  201          CONTINUE
                  IF( KMAT(1).EQ.KSID(1).AND.KMAT(2).EQ.KSID(2) )THEN
                     IF( KSUR(1).EQ.0 .OR. KSUR(2).EQ.0 ) GO TO 209
                     KABSO(1)= NUMBLK( KSUR(1),NO(1) )
                     KABSO(2)= NUMBLK( KSUR(2),NO(2) )
                     KEYINT( KABSO(1) )= KABSO(2)
                     KEYINT( KABSO(2) )= KABSO(1)
                     KSUR(1)= KSUR(1)+1
                     KSUR(2)= KSUR(2)+1
                  ELSE
                     IF( KMAT(1).NE.KSID(1) ) KSUR(1)= KSUR(1)+1
                     IF( KMAT(2).NE.KSID(2) ) KSUR(2)= KSUR(2)+1
                  ENDIF
                  IF( KSUR(1).NE.0 )THEN
                     KMAT(1)= MATTYP( NUMTYP(KSUR(1),KTYP(1)) )
                  ELSE
                     KMAT(1)= KSID(1)
                  ENDIF
                  IF( KSUR(2).NE.0 )THEN
                     KMAT(2)= MATTYP( NUMTYP(KSUR(2),KTYP(2)) )
                  ELSE
                     KMAT(2)= KSID(2)
                  ENDIF
               GO TO 201
  209          IF( KSUR(1).NE.0 .OR. KSUR(2).NE.0 )THEN
                  WRITE(IOUT,'(1H ,I8,4H OF ,I8,5H <=> ,I8,4H OF ,I8)')
     >                          KSUR(1),  NO(1),     KSUR(2),  NO(2)
                  SWKILL=.TRUE.
               ENDIF
  210       CONTINUE
  215       CONTINUE
         ELSE
            CALL LCMLIB(IPGEOM)
            CALL XABORT( 'XELBIN: *** FALSE NDIM VALUE')
         ENDIF
  220 CONTINUE
*
      IF( IPRT.GE.100 .OR. SWKILL )THEN
         IUNK= 0
         WRITE(IOUT,'(/40H       KEYINT       COUPLE    MATERIAL  )')
         DO 250 IBLK= 1, NBLOCK
            ITYP= KEYTYP(IBLK)
            IGEO= KEYGEO(ITYP)
            NVOX= NVOLO(IGEO)
            NSUX= NSURO(IGEO)
            DO 240 ISUX= NSUX, NVOX
               IUNK= IUNK+1
               IMYT= MATTYP( NUMTYP(ISUX,ITYP) )
               IF( ISUX.LT.0 )THEN
                  WRITE(IOUT,
     >            '(5H SUR(,I8,5H) => ,I8,4H OF ,    I8)')
     >                      IUNK,      KEYINT(IUNK), IMYT
               ELSEIF( ISUX.GT.0 )THEN
                  IMYG= MATGEO( NUMGEO(ISUX,IGEO) )
                  WRITE(IOUT,
     >            '(5H VOL(,I8,5H) => ,I8,4H OF ,    I8,1H(,I8,1H))')
     >                      IUNK,      KEYINT(IUNK), IMYT,  IMYG
               ENDIF
  240       CONTINUE
            WRITE(IOUT,'(/1X)')
  250    CONTINUE
      ENDIF
      IF( SWKILL ) CALL XABORT( 'XELBIN: IMPOSSIBLE TO INTERFACE')
*
      RETURN
      END
