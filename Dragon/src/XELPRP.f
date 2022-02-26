*DECK XELPRP
      SUBROUTINE XELPRP(IPGEOM, GEONAM,  NDIM,  NTYPO, NBLOCK,  NBMIX,
     >                  MAXGRI, ALBEDO, ICODE,  NCODE, LCLSYM, LCLTRA,
     >                  MRGSUR, LEAKSW,   LL1,   LL2,  L1CELL, NEXTGE,
     >                  IFCSYM, IPRT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Reads the geometry and check if the geometry 
* is acceptable for EXCELL.
*
*Copyright:
* Copyright (C) 1989 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPGEOM  pointer to the geometry (L_GEOM).             
* GEONAM  geometry name.                               
* IPRT    printing level.                              
*
*Parameters: output
* NDIM    number of dimensions.
* NTYPO   number of types.                                    
* NBLOCK  number of blocks.                                  
* NBMIX   number of mixtures.                                
* MAXGRI  grid dimensions (NX*NY*NZ).                   
* ALBEDO  geometric albedos on the six faces.          
* ICODE   index for boundary conditions.               
* NCODE   type of boundary conditions.                 
* LCLSYM  symmetry flags (0: no; -1/+1: syme; -2/+2: ssym).
* LCLTRA  translation flags (0: no; +1: tra).              
* MRGSUR  similarity between faces.                     
* LEAKSW  leakage switch.                               
* LL1     diagonal symmetry (2,3).                      
* LL2     diagonal symmetry (1,4).                      
* L1CELL  to indicate that there is only 1 cell.        
* NEXTGE  rectangular(0)/circular(1) boundary.          
* IFCSYM  number of symmetry in full assembly (1,2,3,4,5).   
*
*-----------------------------------------------------------------------
*
      USE          GANLIB
      IMPLICIT     NONE
*
      TYPE(C_PTR)  IPGEOM
      INTEGER      NDIM,  NTYPO, NBLOCK, NBMIX, NEXTGE, IFCSYM, IPRT
      INTEGER      MAXGRI(3),LCLSYM(3),LCLTRA(3),
     >             NCODE(6),ICODE(6),MRGSUR(-6:-1)
      LOGICAL      LEAKSW,LL1,LL2,L1CELL
      REAL         ALBEDO(6)
*
      INTEGER      NLCM, NIXS, NSTATE, IOUT
      PARAMETER  ( NLCM=26, NIXS=8, NSTATE=40, IOUT=6 )
      INTEGER      LNLCM(NLCM),INVLCM(NIXS),
     >             ISTATE(NSTATE),JCODE(6)
      REAL         ZCODE(6)
      LOGICAL      SWALBE(6)
      CHARACTER    LCMNM(NLCM)*12, GEONAM*12, CORIEN(-6:0)*4
      INTEGER      ILCM, IDIR, IIXS, ILONG, ITPLCM, ISUR, ITYPE,
     >             LREG, ISUB1, ISUB2, ISPLIT, ITRAN, I2, IAL
*
      DATA         CORIEN
     >             / ' Z+ ',' Z- ',' Y+ ',' Y- ',' X+ ',' X- ','    ' /
      DATA INVLCM/  6, 12, 16, 17, 18, 20, 21, 22 /
      DATA LCMNM /  'MIX',  'MESHX',  'MESHY',   'MESHZ',  'RADIUS',
     >             'SIDE', 'SPLITX', 'SPLITY',  'SPLITZ',  'SPLITR',
     >             'CELL',  'COORD',  'MERGE',    'TURN', 'CLUSTER',
     >             'NPIN',   'RPIN',   'APIN',   'BIHET',  'POURCE', 
     >           'PROCEL',   'IHEX',   'NCODE',   'ZCODE',  'ICODE',
     >           'CENTER'/
*
      IFCSYM=    1
      DO 10 ILCM= 1, NLCM
         CALL LCMLEN(IPGEOM, LCMNM(ILCM), LNLCM(ILCM), ITPLCM )
   10 CONTINUE
      IFCSYM=    1
      DO 11 IDIR=1,3
        LCLSYM(IDIR)=0
        LCLTRA(IDIR)=0
 11   CONTINUE
*
*     ELIMINATES THE INVALID OPTIONS
      DO 20 IIXS= 1, NIXS
        IF( LNLCM(INVLCM(IIXS)).NE.0 )
     >     CALL XABORT( 'XELPRP:*'//GEONAM//'* IS '//
     >                  'NOT A VALID GEOMETRY FOR EXCELL'//
     >                  ' (LCM BLOCK *'//LCMNM(INVLCM(IIXS))//'*)')
   20 CONTINUE
      CALL LCMLEN(IPGEOM,'STATE-VECTOR',ILONG,ITPLCM)
      IF( ILONG.LE.0 .OR. ILONG .GT. NSTATE )
     >   CALL XABORT( 'XELPRP: GEOMETRY HAS INVALID STATE VECTOR')
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      DO 35 ISUR= 1, 6
         SWALBE( ISUR)= .FALSE.
         ALBEDO( ISUR)= 1.0
         MRGSUR(-ISUR)= -ISUR
         ICODE ( ISUR)= -ISUR
   35 CONTINUE
*
      ITYPE=  ISTATE(1)
      LREG=   ISTATE(6)
      NBMIX=  ISTATE(7)
      ISUB1=  ISTATE(8)
      ISUB2=  ISTATE(9)
      ISPLIT= ISTATE(11)
      NEXTGE= 0
*
      IF( ISUB1.NE.0 )THEN
*
*        MANY CELLS
         L1CELL= .FALSE.
         MAXGRI(1)= MAX(1,ISTATE(3))
         MAXGRI(2)= MAX(1,ISTATE(4))
         MAXGRI(3)= MAX(1,ISTATE(5))
         NTYPO= ISUB2
         IF( ITYPE.EQ.5 )THEN
            NDIM= 2
            SWALBE(1)=.TRUE.
            SWALBE(2)=.TRUE.
            SWALBE(3)=.TRUE.
            SWALBE(4)=.TRUE.
            ICODE (5)= 0
            ICODE (6)= 0
         ELSEIF( ITYPE.EQ.7 )THEN
            NDIM= 3
            SWALBE(1)=.TRUE.
            SWALBE(2)=.TRUE.
            SWALBE(3)=.TRUE.
            SWALBE(4)=.TRUE.
            SWALBE(5)=.TRUE.
            SWALBE(6)=.TRUE.
         ELSE
            CALL XABORT( 'XELPRP: INVALID GEOMETRY FOR EXCELL')
         ENDIF
      ELSE
*
*        JUST ONE CELL
         L1CELL= .TRUE.
         MAXGRI(1)= 1
         MAXGRI(2)= 1
         MAXGRI(3)= 1
         NTYPO= 1
         IF( ITYPE.EQ. 3 .OR. ITYPE.EQ. 5 .OR.
     >       ITYPE.EQ.20 )THEN
            NDIM= 2
            IF( ITYPE.EQ.3 )THEN
               NEXTGE= 1
               ICODE (1)= 0
               SWALBE(2)=.TRUE.
               ICODE (3)= 0
               ICODE (4)= 0
               ICODE (5)= 0
               ICODE (6)= 0
            ELSE
               SWALBE(1)=.TRUE.
               SWALBE(2)=.TRUE.
               SWALBE(3)=.TRUE.
               SWALBE(4)=.TRUE.
               ICODE (5)= 0
               ICODE (6)= 0
            ENDIF
         ELSEIF( ITYPE.EQ. 6 .OR. ITYPE.EQ. 7 .OR.
     >           ITYPE.EQ.21 .OR. ITYPE.EQ.22 .OR. ITYPE.EQ.23 )THEN
            NDIM= 3
            IF( ITYPE.EQ.6 )THEN
               NEXTGE= 1
               ICODE (1)= 0
               SWALBE(2)=.TRUE.
               ICODE (3)= 0
               ICODE (4)= 0
               SWALBE(5)=.TRUE.
               SWALBE(6)=.TRUE.
            ELSE
               SWALBE(1)=.TRUE.
               SWALBE(2)=.TRUE.
               SWALBE(3)=.TRUE.
               SWALBE(4)=.TRUE.
               SWALBE(5)=.TRUE.
               SWALBE(6)=.TRUE.
            ENDIF
         ELSE
            CALL XABORT( 'XELPRP: INVALID GEOMETRY FOR EXCELL')
         ENDIF
      ENDIF
*
*     RECOVERS B.C.
      CALL LCMGET(IPGEOM,'NCODE',NCODE)
      CALL LCMGET(IPGEOM,'ZCODE',ZCODE)
      CALL LCMGET(IPGEOM,'ICODE',JCODE)
*
*     TREATMENT OF DIAGONAL B.C.
      LL1= .FALSE.
      LL2= .FALSE.
      ITRAN=0
      I2=0
      DO 50 IAL=1, 6
         IF( .NOT.SWALBE(IAL) ) GO TO 50
         IF( JCODE(IAL).NE.0 )THEN
            IF( ICODE(IAL).EQ.0 )THEN
               CALL XABORT('XELPRP: INVALID BOUNDARY CONDITION.')
            ENDIF
            ICODE(IAL)= JCODE(IAL)
            ZCODE(IAL)= 1.0
         ELSEIF( NCODE(IAL).EQ.0 )THEN
            CALL XABORT('XELPRP: A BOUNDARY CONDITION IS MISSING.')
         ENDIF
         IF( NCODE(IAL).EQ.2 )THEN
            ZCODE(IAL)= 1.0
         ELSEIF( NCODE(IAL).EQ.3 )THEN
            I2=I2+1
         ELSEIF( NCODE(IAL).EQ.4 )THEN
            ITRAN=ITRAN+1
            ZCODE(IAL)= 1.0
         ELSEIF( NCODE(IAL).EQ.6 )THEN
            NCODE(IAL)= 1
         ELSEIF( NCODE(IAL) .EQ. 7 .OR.
     >           NCODE(IAL) .EQ. 8 .OR.
     >           NCODE(IAL) .EQ. 9 .OR.
     >           NCODE(IAL) .GE. 11 )THEN
            CALL XABORT('XELPRP: INVALID B.C. FOR EXCELL')
         ENDIF
   50 CONTINUE
*
*     DIAGONAL  B.C.
      IF( I2.GT.0 )THEN
         IF( I2.NE.2 )
     >     CALL XABORT('XELPRP: NO MORE THAN 2 DIAGONAL CONDITIONS')
         IF( MAXGRI(1).NE.MAXGRI(2))
     >      CALL XABORT('XELPRP: LX=LY WITH A DIAGONAL SYMMETRY.')
         LL1=((NCODE(2).EQ.3).AND.(NCODE(3).EQ.3))
         LL2=((NCODE(1).EQ.3).AND.(NCODE(4).EQ.3))
         IFCSYM= IFCSYM+1
         IF( LL1 )THEN
            NCODE(2)=  NCODE(4)
            NCODE(3)=  NCODE(1)
            ICODE(2)=  ICODE(4)
            ICODE(3)=  ICODE(1)
            MRGSUR(-2)= -4
            MRGSUR(-3)= -1
            ZCODE(2)=  ZCODE(4)
            ZCODE(3)=  ZCODE(1)
         ELSEIF( LL2 )THEN
            NCODE(1)=  NCODE(3)
            NCODE(4)=  NCODE(2)
            ICODE(1)=  ICODE(3)
            ICODE(4)=  ICODE(2)
            MRGSUR(-1)= -3
            MRGSUR(-4)= -2
            ZCODE(1)=  ZCODE(3)
            ZCODE(4)=  ZCODE(2)
         ELSE
            CALL XABORT('XELPRP: THE DIAGONAL CONDITIONS '//
     >                  'X+: DIAG Y-: DIAG AND '//
     >                  'X-: DIAG Y+: DIAG ARE THE ONLY PERMITTED.')
         ENDIF
      ENDIF
*
*     TRANSLATION BC (PERIODIC CELL)
*     ONLY PAIRS PERMITTED:
*       1) X- TRAN X+ TRAN
*       2) Y- TRAN Y+ TRAN
*       3) Z- TRAN Z+ TRAN
      IF( ITRAN.GT.0 )THEN
         IF( MOD(ITRAN,2).EQ.1 )THEN
            CALL XABORT('XELPRP: TRANSLATION SYMETRIES COME IN PAIRS')
         ENDIF
         DO 45 IAL=1,6,2
            IF(SWALBE(IAL)) THEN
               IF( NCODE(IAL).EQ.4 .AND. NCODE(IAL+1).EQ.4 )THEN
                  LCLTRA((IAL+1)/2)=1
                  MRGSUR(-IAL  )=-IAL-1
                  MRGSUR(-IAL-1)=-IAL
                  ITRAN=ITRAN-2
               ENDIF
            ENDIF
 45      CONTINUE
         IF( ITRAN.NE.0 )THEN
            CALL XABORT('XELPRP: WRONG PAIRS OF TRANSLATION SYMETRIES')
         ENDIF
      ENDIF
*
*     SYMMETRIC B.C.
      DO 40 IAL= 1, 6
         IF( .NOT.SWALBE(IAL) ) GO TO 40
         ALBEDO( IAL)= ZCODE(IAL)
         IF( NCODE(IAL).EQ.5 )THEN
            MAXGRI((IAL+1)/2)= 2*MAXGRI((IAL+1)/2)-1
            IF( LCLSYM((IAL+1)/2).NE.0 )THEN
               CALL XABORT('XELPRP: 2 SYMMETRIES ON SAME FACE')
            ELSE
               IFCSYM= IFCSYM+1
               IF( MOD(IAL,2).EQ.0 )THEN
                  LCLSYM((IAL+1)/2)= +1
                  MRGSUR(-IAL)=  MRGSUR(-IAL+1)
                  ALBEDO( IAL)=  ZCODE(IAL-1)
                  ICODE ( IAL)=  ICODE(IAL-1)
               ELSE
                  LCLSYM((IAL+1)/2)= -1
                  MRGSUR(-IAL)=  MRGSUR(-IAL-1)
                  ALBEDO( IAL)=  ZCODE(IAL+1)
                  ICODE ( IAL)=  ICODE(IAL+1)
               ENDIF
            ENDIF
         ELSE IF( NCODE(IAL).EQ.10 )THEN
            MAXGRI((IAL+1)/2)= 2*MAXGRI((IAL+1)/2)
            IF( LCLSYM((IAL+1)/2).NE.0 )THEN
               CALL XABORT('XELPRP: 2 SYMMETRIES ON SAME FACE')
            ELSE
               IFCSYM= IFCSYM+1
               IF( MOD(IAL,2).EQ.0 )THEN
                  LCLSYM((IAL+1)/2)= +2
                  MRGSUR(-IAL)=  MRGSUR(-IAL+1)
                  ALBEDO( IAL)=  ZCODE(IAL-1)
                  ICODE ( IAL)=  ICODE(IAL-1)
               ELSE
                  LCLSYM((IAL+1)/2)= -2
                  MRGSUR(-IAL)=  MRGSUR(-IAL-1)
                  ALBEDO( IAL)=  ZCODE(IAL+1)
                  ICODE ( IAL)=  ICODE(IAL+1)
               ENDIF
            ENDIF
         ENDIF
   40 CONTINUE
*
      NBLOCK= MAXGRI(1)*MAXGRI(2)*MAXGRI(3)
      LEAKSW= .TRUE.
      DO 60 ISUR= 1, 6
         LEAKSW= LEAKSW .AND. ALBEDO( ISUR).EQ.1.0
   60 CONTINUE
      LEAKSW= .NOT. LEAKSW
      IF( IPRT.GT.2 )THEN
         IF( LEAKSW )THEN
          WRITE(IOUT,6000)
     >              (100.*(1.0-ALBEDO(IAL)), IAL= 1,6)
         ELSE
           WRITE(IOUT,6001)
         ENDIF
         WRITE(IOUT,6100)
     >           (CORIEN(MRGSUR(IAL)), IAL=-1,-6,-1)
      ENDIF
      IF( NEXTGE.NE.0 )THEN
           CALL XABORT( 'XELPRP:*'//GEONAM//'* IS '//
     >                  'A TUBE/TUBEZ GEOMETRY (NOT AVAILABLE)')
      ENDIF
*
      RETURN
 6000 FORMAT(/1X,'*** ONLY FOR GEOMETRIC ALBEDOS ***'       
     >       /1X,'PERCENT LEAKAGE X-: ',F5.1,'% X+: ',F5.1,'%'
     >       /1X,'(FULL UNFOLD    Y-: ',F5.1,'% Y+: ',F5.1,'%'
     >       /1X,' ASSEMBLY)      Z-: ',F5.1,'% Z+: ',F5.1,'%'//)
 6001 FORMAT(/1X,'*** ONLY FOR GEOMETRIC ALBEDOS ***'       
     >       /1X,'*** NO LEAKAGE ON THE ASSEMBLY ***'//)
 6100 FORMAT(/1X,'SIMILAR FACES   X-: ',A5,2X,'X+: ',A5
     >       /1X,'(FULL UNFOLD    Y-: ',A5,2X,'Y+: ',A5
     >       /1X,' ASSEMBLY)      Z-: ',A5,2X,'Z+: ',A5//)

      END
