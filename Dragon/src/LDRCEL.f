*DECK LDRCEL
      LOGICAL FUNCTION LDRCEL(IPGEOM,   IT1,    JT1,    IT2,    JT2,
     >                         CELLT, NTYPES, IAXIS,   NDIM,   IPRT )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Check if two cells can be connected in an assembly of cells.
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
* IPGEOM  pointer to the geometry LCM object (L_GEOM signature).
* IT1     type of the  first cell.
* JT1     turn index of the  first cell.
* IT2     type of the second cell.
* JT2     turn index of the second cell.
* CELLT   to keep cell type names.
* NTYPES  number of types.
* IAXIS   axis of the connexion.
* NDIM    number of dimensions.
* IPRT    intermediate printing level for output.
*
*Parameters: output
* LDRCEL  checking flag: =.true. if cells do connect adequately;
*         =.false. if they do not connect adequately.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPGEOM
      INTEGER          IT1,JT1,IT2,JT2,CELLT(3*NTYPES),NTYPES,IAXIS,
     >                 NDIM,IPRT
*----
*  LOCAL VARIABLES
*----
      PARAMETER      ( IOUT=6,NSTATE=40 )
      INTEGER          ISTAT1(NSTATE),ISTAT2(NSTATE),ILEN(3), JLEN(3)
      CHARACTER        GEOC1*12,GEOC2*12,CAXIS(3)*1,GEOCT1*18,GEOCT2*18
      LOGICAL          LDRGEO
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: RIADD1,RIADD2,RIADD3,RJADD1,
     >                                   RJADD2,RJADD3
*----
*  DATA STATEMENTS
*----
      DATA     CAXIS   / 'X', 'Y', 'Z' /
*
      IST2=0
      IF( IAXIS.GT.NDIM )THEN
         LDRCEL= .FALSE.
         CALL XABORT( 'LDRCEL: # OF DIMENSIONS INCOMPATIBLE' )
      ELSE
         LDRCEL= .TRUE.
      ENDIF
      ICX1 = MOD(IAXIS  , NDIM) + 1
      ICX2 = ICX1
      ICY1 = MOD(IAXIS+1, NDIM) + 1
      ICY2 = ICY1
*
*     XY-ROTATES IF NECESSARY FROM JT1 AND JT2
      IF( JT1.EQ.2 )THEN
         IF( ICX1.EQ.2 )THEN
            ICX1= 1
         ELSEIF( ICX1.EQ.1 )THEN
            ICX1= 2
         ENDIF
         IF( ICY1.EQ.2 )THEN
            ICY1= 1
         ELSEIF( ICY1.EQ.1 )THEN
            ICY1= 2
         ENDIF
      ENDIF
      IF( JT2.EQ.2 )THEN
         IF( ICX2.EQ.2 )THEN
            ICX2= 1
         ELSEIF( ICX2.EQ.1 )THEN
            ICX2= 2
         ENDIF
         IF( ICY2.EQ.2 )THEN
            ICY2= 1
         ELSEIF( ICY2.EQ.1 )THEN
            ICY2= 2
         ENDIF
      ENDIF
      N1= 0
      N2= 0
      WRITE( GEOC1(1: 4),'(A4)') CELLT(3*IT1-2)
      WRITE( GEOC1(5: 8),'(A4)') CELLT(3*IT1-1)
      WRITE( GEOC1(9:12),'(A4)') CELLT(3*IT1    )
      IF( IT1.NE.IT2 )THEN
         WRITE( GEOC2(1: 4),'(A4)') CELLT(3*IT2-2 )
         WRITE( GEOC2(5: 8),'(A4)') CELLT(3*IT2-1 )
         WRITE( GEOC2(9:12),'(A4)') CELLT(3*IT2     )
         IF( JT1.EQ.1 )THEN
            GEOCT1= GEOC1//'      '
         ELSE
            GEOCT1= GEOC1//'/TURN '
         ENDIF
         IF( JT2.EQ.1 )THEN
            GEOCT2= GEOC2//'      '
         ELSE
            GEOCT2= GEOC2//'/TURN '
         ENDIF
         IF( IPRT.GT.1 )THEN
            WRITE(IOUT,'(1H ,A1,17H-CELL CONNEXION:  ,2A18)')
     >                    CAXIS(IAXIS),            GEOCT1, GEOCT2
         ENDIF
         IF( .NOT.LDRGEO(IPGEOM, GEOC1, GEOC2, 0) )THEN
*
*           ANALYSE GEOMETRIES
            CALL LCMSIX(IPGEOM, GEOC1, 1)
            CALL XDRSET(ISTAT1,NSTATE,0)
            CALL LCMGET(IPGEOM, 'STATE-VECTOR', ISTAT1)
            IST1= ISTAT1(1)
            IF( IST1.EQ.0 ) THEN
               CALL LCMSIX(IPGEOM, ' ', 2)
               RETURN
            ENDIF
            IF( NDIM.EQ.2 )THEN
               IF( IST1.EQ.5 .OR. IST1.EQ.20 )THEN
                  N1= 1
                  CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(ICX1), ILEN(1), IT)
                  ALLOCATE(RIADD1(ILEN(1)))
                  CALL LCMGET(IPGEOM, 'MESH'//CAXIS(ICX1), RIADD1)
               ELSE
                  CALL XABORT( 'LDRCEL: ONE CELL HAS INCORRECT DIM')
               ENDIF
            ELSEIF( NDIM.EQ.3 )THEN
               IF( IST1.EQ.7  .OR. IST1.EQ.21 .OR.
     >             IST1.EQ.22 .OR. IST1.EQ.23 )THEN
                  N1=2
                  CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(ICX1), ILEN(1), IT)
                  ALLOCATE(RIADD1(ILEN(1)))
                  CALL LCMGET(IPGEOM, 'MESH'//CAXIS(ICX1), RIADD1)
                  CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(ICY1), ILEN(2), IT)
                  ALLOCATE(RIADD2(ILEN(2)))
                  CALL LCMGET(IPGEOM, 'MESH'//CAXIS(ICY1), RIADD2)
                  IF( IAXIS.EQ.IST1-20 )THEN
                     N1=3
                     CALL LCMLEN(IPGEOM, 'RADIUS', ILEN(3), IT)
                     ALLOCATE(RIADD3(ILEN(3)))
                     CALL LCMGET(IPGEOM, 'RADIUS', RIADD3)
                  ENDIF
               ELSE
                  CALL XABORT( 'LDRCEL: ONE CELL HAS INCORRECT DIM')
               ENDIF
            ENDIF
            CALL LCMSIX(IPGEOM, ' ',   2)
*
*           SCAN THE SECOND GEOMETRY
            CALL LCMSIX(IPGEOM, GEOC2, 1)
            CALL XDRSET(ISTAT2,NSTATE,0)
            CALL LCMGET(IPGEOM, 'STATE-VECTOR', ISTAT2)
            IST2= ISTAT2(1)
            IF( IST2.EQ.0 ) THEN
               CALL LCMSIX(IPGEOM, ' ', 2)
               RETURN
            ENDIF
            IF( NDIM.EQ.2 )THEN
               IF( IST2.EQ.5 .OR. IST2.EQ.20 )THEN
                  N2=1
                  CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(ICX2), JLEN(1), IT)
                  ALLOCATE(RJADD1(JLEN(1)))
                  CALL LCMGET(IPGEOM, 'MESH'//CAXIS(ICX2), RJADD1)
               ELSE
                  CALL XABORT( 'LDRCEL: ONE CELL HAS INCORRECT DIM')
               ENDIF
            ELSEIF( NDIM.EQ.3 )THEN
               IF( IST2.EQ.7  .OR. IST2.EQ.21 .OR.
     >             IST2.EQ.22 .OR. IST2.EQ.23 )THEN
                  N2=2
                  CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(ICX2), JLEN(1), IT)
                  ALLOCATE(RJADD1(JLEN(1)))
                  CALL LCMGET(IPGEOM, 'MESH'//CAXIS(ICX2), RJADD1)
                  CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(ICY2), JLEN(2), IT)
                  ALLOCATE(RJADD2(JLEN(2)))
                  CALL LCMGET(IPGEOM, 'MESH'//CAXIS(ICY2), RJADD2)
                  IF( IAXIS.EQ.IST2-20 )THEN
                     N2=3
                     CALL LCMLEN(IPGEOM, 'RADIUS', JLEN(3), IT)
                     ALLOCATE(RJADD3(JLEN(3)))
                     CALL LCMGET(IPGEOM, 'RADIUS', RJADD3)
                  ENDIF
               ELSE
                  CALL XABORT( 'LDRCEL: ONE CELL HAS INCORRECT DIM')
               ENDIF
            ENDIF
            CALL LCMSIX(IPGEOM, ' ',   2)
         ELSE
*
*           GEOMETRIES ARE SIMILAR, CHECK FOR XY-ROTATIONAL INVARIANCE
            IF( JT1.NE.JT2 )THEN
               CALL LCMSIX(IPGEOM, GEOC1, 1)
               CALL XDRSET(ISTAT1,NSTATE,0)
               CALL LCMGET(IPGEOM, 'STATE-VECTOR', ISTAT1)
               IST1= ISTAT1(1)
               IF( IST1.EQ.0 ) THEN
                  CALL LCMSIX(IPGEOM, ' ', 2)
                  RETURN
               ENDIF
               IF( JT1.EQ.1 )THEN
                  GEOCT1= GEOC1//'      '
               ELSE
                  GEOCT1= GEOC1//'/TURN '
               ENDIF
               IF( JT2.EQ.1 )THEN
                  GEOCT2= GEOC1//'      '
               ELSE
                  GEOCT2= GEOC1//'/TURN '
               ENDIF
               IF( IPRT.GT.1 )THEN
                  WRITE(IOUT,'(1H ,A1,17H-CELL CONNEXION:  ,2A18)')
     >                          CAXIS(IAXIS),            GEOCT1, GEOCT2
               ENDIF
               IF( IST1.EQ.5 .OR. IST1.EQ.20 )THEN
                  N1= 1
                  N2= 1
                  CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(1), ILEN(1), IT)
                  ALLOCATE(RIADD1(ILEN(1)))
                  CALL LCMGET(IPGEOM, 'MESH'//CAXIS(1), RIADD1)
                  CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(2), JLEN(1), IT)
                  ALLOCATE(RJADD1(JLEN(1)))
                  CALL LCMGET(IPGEOM, 'MESH'//CAXIS(2), RJADD1) 
               ELSE IF( IST1.EQ.7  .OR. IST1.EQ.21 .OR.
     >             IST1.EQ.22 .OR. IST1.EQ.23 )THEN
                  N1=2
                  CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(ICX1), ILEN(1), IT)
                  ALLOCATE(RIADD1(ILEN(1)))
                  CALL LCMGET(IPGEOM, 'MESH'//CAXIS(ICX1), RIADD1)
                  CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(ICY1), ILEN(2), IT)
                  ALLOCATE(RIADD2(ILEN(2)))
                  CALL LCMGET(IPGEOM, 'MESH'//CAXIS(ICY1), RIADD2)
                  IF( IAXIS.EQ.IST1-20 )THEN
                     N1=3
                     CALL LCMLEN(IPGEOM, 'RADIUS', ILEN(3), IT)
                     ALLOCATE(RIADD3(ILEN(3)))
                     CALL LCMGET(IPGEOM, 'RADIUS', RIADD3)
                  ENDIF
                  N2=2
                  CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(ICX2), JLEN(1), IT)
                  ALLOCATE(RJADD1(JLEN(1)))
                  CALL LCMGET(IPGEOM, 'MESH'//CAXIS(ICX2), RJADD1)
                  CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(ICY2), JLEN(2), IT)
                  ALLOCATE(RJADD2(JLEN(2)))
                  CALL LCMGET(IPGEOM, 'MESH'//CAXIS(ICY2), RJADD2)
                  IF( IAXIS.EQ.IST2-20 )THEN
                     N2=3
                     CALL LCMLEN(IPGEOM, 'RADIUS', JLEN(3), IT)
                     ALLOCATE(RJADD3(JLEN(3)))
                     CALL LCMGET(IPGEOM, 'RADIUS', RJADD3)
                  ENDIF
               ELSE
                  CALL XABORT( 'LDRCEL: ONE CELL HAS INCORRECT DIM')
               ENDIF
               CALL LCMSIX(IPGEOM, ' ',   2)
            ENDIF
         ENDIF
      ELSE
*
*     GEOMETRY ARE OF THE SAME TYPE, CHECK FOR XY-ROTATIONAL INVARIANCE
         IF( JT1.NE.JT2 )THEN
            CALL LCMSIX(IPGEOM, GEOC1, 1)
            CALL XDRSET(ISTAT1,NSTATE,0)
            CALL LCMGET(IPGEOM, 'STATE-VECTOR', ISTAT1)
            IST1= ISTAT1(1)
            IF( IST1.EQ.0 ) THEN
               CALL LCMSIX(IPGEOM, ' ', 2)
               RETURN
            ENDIF
            IF( JT1.EQ.1 )THEN
               GEOCT1= GEOC1//'      '
            ELSE
               GEOCT1= GEOC1//'/TURN '
            ENDIF
            IF( JT2.EQ.1 )THEN
               GEOCT2= GEOC1//'      '
            ELSE
               GEOCT2= GEOC1//'/TURN '
            ENDIF
            IF( IPRT.GT.1 )THEN
               WRITE(IOUT,'(1H ,A1,17H-CELL CONNEXION:  ,2A18)')
     >                       CAXIS(IAXIS),            GEOCT1, GEOCT2
            ENDIF
            IF( IST1.EQ.5 .OR. IST1.EQ.20 )THEN
               N1= 1
               N2= 1
               CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(1), ILEN(1), IT)
               ALLOCATE(RIADD1(ILEN(1)))
               CALL LCMGET(IPGEOM, 'MESH'//CAXIS(1), RIADD1)
               CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(2), JLEN(1), IT)
               ALLOCATE(RJADD1(JLEN(1)))
               CALL LCMGET(IPGEOM, 'MESH'//CAXIS(2), RJADD1)
            ELSE IF( IST1.EQ.7  .OR. IST1.EQ.21 .OR.
     >          IST1.EQ.22 .OR. IST1.EQ.23 )THEN
               N1=2
               CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(ICX1), ILEN(1), IT)
               ALLOCATE(RIADD1(ILEN(1)))
               CALL LCMGET(IPGEOM, 'MESH'//CAXIS(ICX1), RIADD1)
               CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(ICY1), ILEN(2), IT)
               ALLOCATE(RIADD2(ILEN(2)))
               CALL LCMGET(IPGEOM, 'MESH'//CAXIS(ICY1), RIADD2)
               IF( IAXIS.EQ.IST1-20 )THEN
                  N1=3
                  CALL LCMLEN(IPGEOM, 'RADIUS', ILEN(3), IT)
                  ALLOCATE(RIADD3(ILEN(3)))
                  CALL LCMGET(IPGEOM, 'RADIUS', RIADD3)
               ENDIF
               N2=2
               CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(ICX2), JLEN(1), IT)
               ALLOCATE(RJADD1(JLEN(1)))
               CALL LCMGET(IPGEOM, 'MESH'//CAXIS(ICX2), RJADD1)
               CALL LCMLEN(IPGEOM, 'MESH'//CAXIS(ICY2), JLEN(2), IT)
               ALLOCATE(RJADD2(JLEN(2)))
               CALL LCMGET(IPGEOM, 'MESH'//CAXIS(ICY2), RJADD2)
               IF( IAXIS.EQ.IST2-20 )THEN
                  N2=3
                  CALL LCMLEN(IPGEOM, 'RADIUS', JLEN(3), IT)
                  ALLOCATE(RJADD3(JLEN(3)))
                  CALL LCMGET(IPGEOM, 'RADIUS', RJADD3)
               ENDIF
            ELSE
               CALL XABORT( 'LDRCEL: ONE CELL HAS INCORRECT DIM')
            ENDIF
            CALL LCMSIX(IPGEOM, ' ',   2)
         ENDIF
      ENDIF
*
*     CHECK VECTORS THAT HAVE BEEN EXTRACTED AND RELEASE SPACE...
      IF( N1.EQ.N2 )THEN
         DO 40 I= 1, N1
            IF( ILEN(I).NE.JLEN(I) )
     >         CALL XABORT( 'LDRCEL: INCORRECT MESHING DIM')
            IF(I.EQ.1) THEN
               DO 10 J= 1, ILEN(1)
               IF( RIADD1(J).NE.RJADD1(J))
     >            CALL XABORT( 'LDRCEL: INCORRECT MESHING(1)')
   10          CONTINUE
               DEALLOCATE(RJADD1,RIADD1)
            ELSE IF(I.EQ.2) THEN
               DO 20 J= 1, ILEN(2)
               IF( RIADD2(J).NE.RJADD2(J))
     >            CALL XABORT( 'LDRCEL: INCORRECT MESHING(2)')
   20          CONTINUE
               DEALLOCATE(RJADD2,RIADD2)
            ELSE IF(I.EQ.3) THEN
               DO 30 J= 1, ILEN(3)
               IF( RIADD3(J).NE.RJADD3(J))
     >            CALL XABORT( 'LDRCEL: INCORRECT MESHING(3)')
   30          CONTINUE
               DEALLOCATE(RJADD3,RIADD3)
            ENDIF
   40    CONTINUE
      ELSE
         CALL XABORT( 'LDRCEL: TYPES ARE INCOMPATIBLE')
      ENDIF
      RETURN
      END
