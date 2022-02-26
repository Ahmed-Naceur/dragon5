*DECK LDRASS
      LOGICAL FUNCTION LDRASS(IPGEOM,IPRT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Reads the geometry on LCM and check compatibility for cell assemblies.
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
* IPGEOM  pointer to the geometry LCM object (L_GEOM).
* IPRT    print flag (iprt=0: no print).
*
*Parameters: output
* LDRASS  checking flag: =.true. if everything was 'ok' with assembly;
*         =.false. if nothing was checked.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)  IPGEOM
      INTEGER      IPRT
*----
*  LOCAL VARIABLES
*----
      PARAMETER  ( IOUT=6, NLCM=26, NSTATE=40, NIXS=12, NIST=1 )
      CHARACTER*12 LCMNM(NLCM), GEONAM, TEXT12
      INTEGER      LNLCM(NLCM),INVLCM(NIXS),INVSTA(NIST),ISTATE(NSTATE),
     >             MAXGRI(3),NCODE(6)
      LOGICAL      LL1,LL2,LCELL,LDRCEL,EMPTY,LCM
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEYTYP,KMERGE,KTURN,KEYCEL,
     > KCHECK
*----
*  DATA STATEMENTS
*----
      DATA INVLCM/ 6,  7,  8,  9, 10,     12,     16, 17, 18,
     >            20, 21, 22 /
      DATA INVSTA/11 /
      DATA LCMNM /   'MIX',  'MESHX',  'MESHY',  'MESHZ', 'RADIUS',
     >              'SIDE', 'SPLITX', 'SPLITY', 'SPLITZ', 'SPLITR',
     >              'CELL',  'COORD',  'MERGE',   'TURN','CLUSTER',
     >              'NPIN',   'RPIN',   'APIN',  'BIHET', 'POURCE',
     >            'PROCEL',   'IHEX',  'STATE',  'NCODE',  'ZCODE',
     >             'ICODE'/
*
*     I: CELL TYPE
*     J: CELL TURN (=J.EQ.1.OR.J.EQ.2)
      IADR( IAXIS, I1,J1, I2,J2)= IAXIS
     >     + 3*(J1-1+2*(J2-1+2*(I1-1+NTYPES*(I2-1))))
*
      CALL LCMINF(IPGEOM,GEONAM,TEXT12,EMPTY,ILONG,LCM)
      LDRASS= .TRUE.
      IF( IPRT.GT.1 )THEN
         WRITE(IOUT,'(/20H CHECKING ASSEMBLY:   ,A12)') GEONAM
      ENDIF
      DO 10 ILCM= 1, NLCM
         CALL LCMLEN(IPGEOM, LCMNM(ILCM), LNLCM(ILCM), ITPLCM )
   10 CONTINUE
*----
*  ELIMINATES OPTIONS NOT CHECKED BY THE ROUTINE
*----
      DO 20 IIXS= 1, NIXS
        IF( LNLCM(INVLCM(IIXS)).NE.0 )THEN
           LDRASS= .FALSE.
           GO TO 9999
        ENDIF
   20 CONTINUE
      CALL LCMLEN(IPGEOM,'STATE-VECTOR',ILEN,ITPLCM)
      IF( ITPLCM.NE.1 )THEN
         LDRASS= .FALSE.
         CALL XABORT( 'LDRASS: STATE VECTOR IS NOT AN INTEGER')
      ENDIF
      IF( ILEN.NE.NSTATE )THEN
         LDRASS= .FALSE.
         CALL XABORT( 'LDRASS: GEOMETRY HAS INVALID STATE VECTOR')
      ENDIF
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
*----
*  ELIMINATES THE INVALID OPTIONS
*----
      DO 30 IIST= 1, NIST
        IF( ISTATE(INVSTA(IIST)).NE.0 )THEN
           LDRASS= .FALSE.
           GO TO 9999
        ENDIF
   30 CONTINUE
*
      ITYPE=  ISTATE(1)
*----
*  FIX DIMENSION OF CELL ASSEMBLY
*----
      IF( ITYPE.EQ.5 )THEN
         NDIM= 2
      ELSEIF( ITYPE.EQ.7 )THEN
         NDIM= 3
      ELSE
         LDRASS= .FALSE.
         GO TO 9999
      ENDIF
      LREG=   ISTATE(6)
      NBMIX=  ISTATE(7)
      ISUB1=  ISTATE(8)
      ISUB2=  ISTATE(9)
*
      IF( ISUB1.NE.0 )THEN
         MAXGRI(1)= MAX(1,ISTATE(3))
         MAXGRI(2)= MAX(1,ISTATE(4))
         MAXGRI(3)= MAX(1,ISTATE(5))
         NBLOCK= MAXGRI(1)*MAXGRI(2)*MAXGRI(3)
         IF( NBLOCK.EQ.1 )THEN
*
*           JUST ONE CELL
            LDRASS=.FALSE.
            GO TO 9999
         ENDIF
*
*        MANY CELLS
         NTYPES= ISUB2
         IF( IPRT.GT.1 )THEN
            WRITE(IOUT,'(6H      ,I1,13H-D ASSEMBLY:   ,3I4)')
     >                         NDIM,               (MAXGRI(I),I=1,NDIM)
         ENDIF
      ELSE
*
*        JUST ONE CELL
         LDRASS=.FALSE.
         GO TO 9999
      ENDIF
*----
*  RECOVERS BOUNDARY CONDITIONS.
*----
      CALL LCMLEN(IPGEOM,'NCODE', ILEN, ITPLCM)
      IF( ITPLCM.NE.1 )
     >   CALL XABORT('LDRASS: THE NCODE BLOCK IS '//
     >               'NOT ADEQUATELY DEFINED')
      IF( ILEN.NE.6 )
     >   CALL XABORT('LDRASS: THE NCODE BLOCK HAS '//
     >               'INCORRECT DIMENSION')
      CALL LCMGET(IPGEOM,'NCODE',NCODE)
      IF( NDIM.EQ.2.AND.(NCODE(5).NE.0.OR.NCODE(6).NE.0) )
     >   CALL XABORT('LDRASS: 3-D NCODE VALUES FOR A 2-D ASSEMBLY...')
      NOCELL= NBLOCK
      NDIAG=0
      DO 40 IAL= 1, 2*NDIM
         IF( NCODE(IAL).EQ.0 )
     >      CALL XABORT('LDRASS: A BOUNDARY CONDITION IS MISSING.')
         IF( NCODE(IAL).EQ.3 )THEN
            NDIAG=NDIAG+1
         ENDIF
   40 CONTINUE
*
      LL1= .FALSE.
      LL2= .FALSE.
      MXDIAG=0
      NOC1=0
      NOC2=0
      NOCO=0
      IF( NDIAG.GT.0 )THEN
         IF( NDIAG.NE.2 )
     >      CALL XABORT('LDRASS: NO MORE THAN 2 DIAGONAL CONDITIONS')
         IF( MAXGRI(1).NE.MAXGRI(2))
     >      CALL XABORT('LDRASS: LX=LY WITH A DIAGONAL SYMMETRY.')
         MXDIAG= MAXGRI(1)
         NOC1=((MXDIAG+1)*MXDIAG)/2
         NOC2=MXDIAG*MXDIAG
         LL1=((NCODE(2).EQ.3).AND.(NCODE(3).EQ.3))
         LL2=((NCODE(1).EQ.3).AND.(NCODE(4).EQ.3))
         IF( LL1 )THEN
            NCODE(2)= NCODE(1)
            NCODE(3)= NCODE(4)
         ELSEIF( LL2 )THEN
            NCODE(1)= NCODE(2)
            NCODE(4)= NCODE(3)
         ELSE
            CALL XABORT('LDRASS: THE DIAGONAL CONDITIONS '//
     >                  'X+: DIAG Y-: DIAG OR  '//
     >                  'X-: DIAG Y+: DIAG ARE THE ONLY PERMITTED.')
         ENDIF
         NOCO=NOC1
         IF(NOC2*MAXGRI(3) .EQ. LREG) NOCO=NOC2 
         NOCELL=NOCO*MAXGRI(3) 
      ENDIF
      IF( NOCELL.GT.LREG )THEN
         CALL XABORT('LDRASS: # OF CELLS OF ASSEMBLY TOO LARGE...')
      ELSEIF( NOCELL.LT.LREG )THEN
         CALL XABORT('LDRASS: # OF CELLS OF ASSEMBLY TOO SMALL...')
      ENDIF
*----
*  CHECK IF GENERATING CELL VECTOR IS 'OK'
*----
      CALL LCMLEN(IPGEOM,'MIX', KCELL, ITPLCM)
      IF( ITPLCM.NE.1 )
     >   CALL XABORT('LDRASS: THE MIX BLOCK WITH CELLS IS '//
     >               'NOT ADEQUATELY DEFINED')
      IF( KCELL.GT.NOCELL )THEN
         CALL XABORT('LDRASS: THE ASSEMBLY HAS TOO MANY CELLS...')
      ELSEIF( KCELL.LT.NOCELL )THEN
         CALL XABORT('LDRASS: THE ASSEMBLY HAS NOT ENOUGH CELLS...')
      ENDIF
      ALLOCATE(KEYTYP(NBLOCK))
      CALL LCMGET(IPGEOM,'MIX', KEYTYP)
      DO 41 I=1,KCELL
      KEYTYP(I)=-KEYTYP(I)
   41 CONTINUE
*----
*  CHECK IF 'MERGE' ARE CORRECTLY DEFINED
*----
      CALL LCMLEN(IPGEOM,'MERGE', KCELL, ITPLCM)
      IF( KCELL.NE.0 )THEN
         IF( ITPLCM.NE.1 )
     >      CALL XABORT('LDRASS: THE MERGE BLOCK IS '//
     >                  'NOT ADEQUATELY DEFINED')
         IF( KCELL.GT.NOCELL )THEN
            CALL XABORT('LDRASS: THE ASSEMBLY HAS TOO MANY MERGE...')
         ELSEIF( KCELL.LT.NOCELL )THEN
            CALL XABORT('LDRASS: THE ASSEMBLY HAS NOT ENOUGH MERGE...')
         ENDIF
         ALLOCATE(KMERGE(NOCELL))
         CALL LCMGET(IPGEOM,'MERGE', KMERGE)
         DO 46 IM= 1, NOCELL
         IDKT= KMERGE(IM)
         IF( IDKT.GT.NOCELL )
     >      CALL XABORT( 'LDRASS: MERGE NUMBER > TOTAL # OF CELLS')
         DO 45 JM= IM+1, NOCELL
            JDKT= KMERGE(JM)
            IF( JDKT.EQ.IDKT )THEN
               IF( KEYTYP(IDKT).NE.KEYTYP(JDKT) )THEN
                  CALL XABORT( 'LDRASS: MERGE NUMBERS ARE NOT '//
     >                         'CONSISTENT WITH GEOMETRIC DEFINITION')
               ENDIF
            ENDIF
   45    CONTINUE
   46    CONTINUE
         DEALLOCATE(KMERGE)
      ENDIF
      ALLOCATE(KTURN(NBLOCK))
      DO 47 IT= 1, NOCELL
         KTURN(IT)= 1
   47 CONTINUE
      CALL LCMLEN(IPGEOM,'TURN', KCELL, ITPLCM)
      IF( KCELL.NE.0 )THEN
         IF( ITPLCM.NE.1 )
     >      CALL XABORT('LDRASS: THE MERGE BLOCK IS '//
     >                  'NOT ADEQUATELY DEFINED')
         IF( KCELL.GT.NOCELL )THEN
            CALL XABORT('LDRASS: THE ASSEMBLY HAS TOO MANY TURN...')
         ELSEIF( KCELL.LT.NOCELL )THEN
            CALL XABORT('LDRASS: THE ASSEMBLY HAS NOT ENOUGH TURN...')
         ENDIF
         CALL LCMGET(IPGEOM,'TURN', KTURN)
      ENDIF
*----
*  EFFECTIVE MOD2 FOR TURN IN CELLS
*----
      DO 48 IT= 1, NOCELL
         KTURN(IT)= MOD( KTURN(IT)+1,2 )+1
   48 CONTINUE
*----
*  CHECK IF CELL NAMES ARE 'OK'
*----
      CALL LCMLEN(IPGEOM,'CELL', KTYPES, ITPLCM)
      IF( ITPLCM.NE.3 )
     >   CALL XABORT('LDRASS: THE CELL NAMES ARE NOT STORED IN'
     >   //' CHARACTER*4')
      IF( KTYPES.GT.3*NTYPES )THEN
         CALL XABORT('LDRASS: THE ASSEMBLY HAS TOO MANY CELL NAMES')
      ELSEIF( KTYPES.LT.3*NTYPES )THEN
         CALL XABORT('LDRASS: THE ASSEMBLY HAS NOT ENOUGH CELLS NAMES')
      ENDIF
      ALLOCATE(KEYCEL(3*NTYPES),KCHECK(12*NTYPES*NTYPES))
      DO 50 IT= 1,12*NTYPES*NTYPES
         KCHECK(IT)= 0
   50 CONTINUE
      CALL LCMGET(IPGEOM,'CELL', KEYCEL)
*----
*  FILL UP "KEYTYP" ARRAY IN THE CASE OF DIAGONAL SYMMETRY
*----
      IF( LL1 )THEN
         IF(NOCO .EQ. NOC1) THEN
*----
*  LOCATE DIAGONAL ELEMENTS IN THEIR RESPECTIVE PLANES
*  WHILE UNFOLDING
*----
         K=LREG
         DO 72 IZ=MAXGRI(3),1,-1
         IOFF=(IZ-1)*MXDIAG*MXDIAG
         DO 71 IY=MXDIAG,1,-1
         DO 60 IX=MXDIAG,IY+1,-1 
         KEYTYP(IOFF+(IY-1)*MXDIAG+IX)=KEYTYP(IOFF+(IX-1)*MXDIAG+IY)
         KTURN(IOFF+(IY-1)*MXDIAG+IX)=KTURN(IOFF+(IX-1)*MXDIAG+IY)
   60    CONTINUE
         DO 70 IX=IY,1,-1
         KEYTYP(IOFF+(IY-1)*MXDIAG+IX)=KEYTYP(K)
         KTURN(IOFF+(IY-1)*MXDIAG+IX)=KTURN(K)
         K=K-1
   70    CONTINUE
   71    CONTINUE
   72    CONTINUE
         DO 77 IZ=MAXGRI(3),1,-1
         IOFF=(IZ-1)*MXDIAG*MXDIAG
         DO 76 IY=1,MXDIAG
         DO 75 IX=IY+1,MXDIAG
         KTURN(IOFF+(IY-1)*MXDIAG+IX)=
     >      MOD(KTURN(IOFF+(IY-1)*MXDIAG+IX),2) + 1
   75    CONTINUE
   76    CONTINUE
   77    CONTINUE
         IF (K.NE.0)
     >     CALL XABORT( 'LDRASS: UNABLE TO UNFOLD '//
     >                  'X+: DIAG Y-: DIAG ASSEMBLY...')
         ENDIF
      ELSEIF( LL2 )THEN
         IF(NOCO .EQ. NOC1) THEN
*----
*  LOCATE DIAGONAL ELEMENTS IN THEIR RESPECTIVE PLANES
*----
         K=LREG
         DO 82 IZ=MAXGRI(3),1,-1
         IOFF=(IZ-1)*MXDIAG*MXDIAG
         DO 81 IY=MXDIAG,1,-1
         DO 80 IX=MXDIAG,IY,-1
         KEYTYP(IOFF+(IY-1)*MXDIAG+IX)=KEYTYP(K)
         KTURN(IOFF+(IY-1)*MXDIAG+IX)=KTURN(K)
         K=K-1
   80    CONTINUE
   81    CONTINUE
   82    CONTINUE
*----
*  UNFOLD DIAGONAL ELEMENTS FOR EACH PLANE
*----
         DO 92 IZ=1,MAXGRI(3)
         IOFF=(IZ-1)*MXDIAG*MXDIAG
         DO 91 IY=1,MXDIAG
         DO 90 IX=1,IY-1
         KEYTYP(IOFF+(IY-1)*MXDIAG+IX)=KEYTYP(IOFF+(IX-1)*MXDIAG+IY)
         KTURN(IOFF+(IY-1)*MXDIAG+IX)=KTURN(IOFF+(IX-1)*MXDIAG+IY)
   90    CONTINUE
   91    CONTINUE
   92    CONTINUE
         IF (K.NE.0)
     >     CALL XABORT( 'LDRASS: UNABLE TO UNFOLD '//
     >                  'X-: DIAG Y+: DIAG ASSEMBLY...')
         DO 97 IZ=MAXGRI(3),1,-1
         IOFF=(IZ-1)*MXDIAG*MXDIAG
         DO 96 IY=1,MXDIAG
         DO 95 IX=1,IY-1
         KTURN(IOFF+(IY-1)*MXDIAG+IX)=
     >      MOD(KTURN(IOFF+(IY-1)*MXDIAG+IX),2) + 1
   95    CONTINUE
   96    CONTINUE
   97    CONTINUE
         ENDIF
      ENDIF
      IF(IPRT .GE. 10) THEN
        WRITE(IOUT,6100) 
        DO 600 IZ=MAXGRI(3),1,-1
          DO 601 IY=MAXGRI(2),1,-1 
            IOFF=(IZ-1)*MAXGRI(2)*MAXGRI(1)+(IY-1)*MAXGRI(1)
            WRITE(IOUT,6110) (KEYTYP(IOFF+IX),IX=1,MAXGRI(1))
 601      CONTINUE
          WRITE(IOUT,6111) 
 600    CONTINUE
        WRITE(IOUT,6101)
        DO 610 IZ=MAXGRI(3),1,-1
          DO 611 IY=MAXGRI(2),1,-1 
            IOFF=(IZ-1)*MAXGRI(2)*MAXGRI(1)+(IY-1)*MAXGRI(1)
            WRITE(IOUT,6110) (KTURN(IOFF+IX),IX=1,MAXGRI(1)) 
 611      CONTINUE
          WRITE(IOUT,6111) 
 610    CONTINUE
      ENDIF
*----
*  TRANSLATION B.C.: CHECK BEGIN-TO-END CONNEXIONS
*----
      DO 100 IC= 1, NDIM
         IF( NCODE(2*IC-1).EQ.4 )THEN
            IF( NCODE(2*IC).NE.4 )
     >      CALL XABORT( 'LDRASS: TRANSLATION B.C. IS NOT WELL DEFINED')
         ENDIF
  100 CONTINUE
*----
*  CHECK CELL INTERFACES
*----
      IOFF1= 0
      DO 112 IZ= 1, MAXGRI(3)
      DO 111 IY= 1, MAXGRI(2)
      DO 110 IX= 1, MAXGRI(1)
         IOFF1= IOFF1+1
         IT1= KEYTYP(IOFF1)
         JT1= KTURN(IOFF1)
         IF(IX.NE.MAXGRI(1).OR.(IX.EQ.MAXGRI(1).AND.NCODE(1).EQ.4))THEN
            IF(IX.NE.MAXGRI(1))THEN
               IOFF2= IOFF1 + 1
            ELSE
               IOFF2= IOFF1 + 1 - MAXGRI(1)
            ENDIF
            IT2= KEYTYP(IOFF2)
            JT2= KTURN(IOFF2)
            IF( KCHECK(IADR(1,IT1,JT1,IT2,JT2)).EQ.0 )THEN
               LCELL= LDRCEL(IPGEOM, IT1,JT1, IT2,JT2, KEYCEL,
     >                       NTYPES,  1,  NDIM, IPRT)
               KCHECK(IADR(1,IT1,JT1,IT2,JT2))= 1
            ENDIF
         ENDIF
         IF(IY.NE.MAXGRI(2).OR.(IY.EQ.MAXGRI(2).AND.NCODE(3).EQ.4))THEN
            IF( IY.NE.MAXGRI(2) )THEN
               IOFF2= IOFF1 + MAXGRI(1)
            ELSE
               IOFF2= IOFF1 +(1-MAXGRI(2)) * MAXGRI(1)
            ENDIF
            IT2= KEYTYP(IOFF2)
            JT2= KTURN(IOFF2)
            IF( KCHECK(IADR(2,IT1,JT1,IT2,JT2)).EQ.0 )THEN
               LCELL= LDRCEL(IPGEOM, IT1,JT1, IT2,JT2, KEYCEL,
     >                       NTYPES, 2,  NDIM, IPRT)
               KCHECK(IADR(2,IT1,JT1,IT2,JT2))= 1
            ENDIF
         ENDIF
         IF(IZ.NE.MAXGRI(3).OR.(IZ.EQ.MAXGRI(3).AND.NCODE(5).EQ.4))THEN
            IF( IZ.NE.MAXGRI(3) )THEN
               IOFF2= IOFF1 + MAXGRI(1)*MAXGRI(2)
            ELSE
               IOFF2= IOFF1 +(1-MAXGRI(3)) * MAXGRI(1) * MAXGRI(2)
            ENDIF
            IT2= KEYTYP(IOFF2)
            JT2= KTURN(IOFF2)
            IF( KCHECK(IADR(3,IT1,JT1,IT2,JT2)).EQ.0 )THEN
               LCELL= LDRCEL(IPGEOM, IT1,JT1, IT2,JT2, KEYCEL,
     >                       NTYPES, 3,  NDIM, IPRT)
               KCHECK(IADR(3,IT1,JT1,IT2,JT2))= 1
            ENDIF
         ENDIF
  110 CONTINUE
  111 CONTINUE
  112 CONTINUE
      DEALLOCATE(KCHECK,KEYCEL,KTURN,KEYTYP)
 9999 IF( IPRT.GT.0 )THEN
        IF( LDRASS )THEN
          WRITE(IOUT,'(1H ,A12,25H ASSEMBLY IS **OK**         )') GEONAM
        ELSE
          WRITE(IOUT,'(1H ,A12,25H WAS NOT ASSEMBLY CHECKED   )') GEONAM
        ENDIF
      ENDIF
      RETURN
*----
*  FORMATS
*----
 6100 FORMAT(/' Assembly by cell number' )
 6101 FORMAT(/' Assembly turns' )
 6110 FORMAT(40(1X,I5))
 6111 FORMAT(' ')
      END
