*DECK AXGDIA
      SUBROUTINE AXGDIA( IPGEOM,   IPRT, NBLOCK,  NTYPO,   NXYZ,  KMESH,
     >                   GEONAM,    LL1,   LL2,   MINGRI, CELLT, KEYTYP,
     >                   ITGEOM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Unfold assembly or cell according to diagonal $x-y$ symmetry 
* and verify if the symmetry is valid.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy and G. Marleau
*
*Parameters: input
* IPGEOM  pointer to the reference geometry data structure.
* IPRT    intermediate printing level for output.
* NBLOCK  number of block in geometry.
* NTYPO   number of types in geometry.
* NXYZ    maximum mesh size in directions $x$, $y$ and $z$.
* KMESH   number of mesh intervals in the geometry.
* GEONAM  name of the reference geometry.
* LL1     flag that is .TRUE. when the diagonal symmetry
*         is applied to surfaces X+ and Y- 
*         (upper diagonal symmetry).
* LL2     flag that is .TRUE. when the diagonal symmetry
*         is applied to surfaces X- and Y+
*         (lower diagonal symmetry).
* MINGRI  minimum grid cell in $x$, $y$ and $z$ directions.
* CELLT   cell type name.
*
*Parameters: input/output
* KEYTYP  type key for each block.
* ITGEOM  turn key associated with each cell type.
*
*External functions
* LELCSY  to verify if a geometry possesses the required internal
*         symmetry.
* AXGTRS  to modify current turn according to required internal
*         symmetry.
* AXGTRN  to associate a DRAGON turn name to a specific turn key. 
*
*-----------------------------------------------------------------------
*
      USE                GANLIB
      IMPLICIT           NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)        IPGEOM
      INTEGER            IPRT,NBLOCK,NTYPO,NXYZ,KMESH
      CHARACTER          GEONAM*12
      LOGICAL            LL1,LL2
      INTEGER            MINGRI(3),CELLT(3*NTYPO)
      INTEGER            KEYTYP(NBLOCK),ITGEOM(NBLOCK)
*----
*  EXTERNAL FUNCTIONS
*----
      LOGICAL            LELCSY
      INTEGER            AXGTRS
      CHARACTER          AXGTRN*2
*----
*  LOCAL VARIABLES
*----
      INTEGER            IOUT
      CHARACTER          NAMSBR*6
      PARAMETER         (IOUT=6,NAMSBR='AXGDIA')
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISPLT,ISPLT1
      REAL, ALLOCATABLE, DIMENSION(:) :: MESH,MESH1
*----
*  LOCAL PARAMETERS
*----
      INTEGER            KML,IX,IY,IZ,IOFF,IOF1,IOF2
      INTEGER            IKG,IKT(2)
      LOGICAL            VALSYM
      CHARACTER          GEOCV*12
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ISPLT(NXYZ),ISPLT1(3*3*NXYZ))
      ALLOCATE(MESH(NXYZ+1),MESH1(2*3*3*(NXYZ+1)))
*----
*  ANALYSE LL1 SYMMETRY (UPPER DIAGONAL SYMMETRY)
*----
      KML=KMESH 
      IF( LL1 )THEN
        DO 100 IZ=MINGRI(3),1,-1
          IOFF=(IZ-1)*MINGRI(1)*MINGRI(2)
          DO 110 IY=MINGRI(2),1,-1
            DO 120 IX=MINGRI(1),IY+1,-1
              KEYTYP(IOFF+(IY-1)*MINGRI(1)+IX)=
     >                      KEYTYP(IOFF+(IX-1)*MINGRI(2)+IY)
              ITGEOM(IOFF+(IY-1)*MINGRI(1)+IX)=
     >                  AXGTRS(ITGEOM(IOFF+(IX-1)*MINGRI(2)+IY),3)
 120        CONTINUE
            DO 130 IX=IY,1,-1
              KEYTYP(IOFF+(IY-1)*MINGRI(1)+IX)=KEYTYP(KML)
              ITGEOM(IOFF+(IY-1)*MINGRI(1)+IX)=ITGEOM(KML)
              IOF1=KML
              IOF2=IOFF+(IY-1)*MINGRI(1)+IX
              IF(IX .EQ. IY) THEN
                IKG=KEYTYP(IOF1)
                IKT(1)=ITGEOM(IOF1)
                IKT(2)=AXGTRS(IKT(1),3)
                WRITE(GEOCV,'(3A4)')
     >          CELLT(3*IKG-2),CELLT(3*IKG-1),CELLT(3*IKG)
                IF(GEOCV .EQ. '            ') THEN
                  IF(IPRT .GT. 10)        
     >              WRITE(IOUT,8000) NAMSBR,'X-Y',
     >              GEONAM,AXGTRN(IKT(1)),AXGTRN(IKT(2))
                ELSE
                  IF(IPRT .GT. 10)        
     >              WRITE(IOUT,8000) NAMSBR,'X-Y',
     >              GEOCV,AXGTRN(IKT(1)),AXGTRN(IKT(2))
                ENDIF
                VALSYM=LELCSY(IPGEOM,IPRT,GEONAM,GEOCV,NXYZ,IKT,
     >                        MESH,ISPLT,MESH1,ISPLT1)
                IF(.NOT. VALSYM) THEN 
                   WRITE(IOUT,8001) 'X-Y',GEOCV,
     >               AXGTRN(IKT(1)),AXGTRN(IKT(2))
                   CALL XABORT(NAMSBR//': INVALID SYMMETRY FOR CELL')
                ENDIF
              ENDIF
              KML=KML-1
 130        CONTINUE  
 110      CONTINUE  
 100    CONTINUE  
      ELSE IF( LL2 )THEN
*----
*  ANALYSE LL2 SYMMETRY (LOWER DIAGONAL SYMMETRY)
*----
        DO 200 IZ=MINGRI(3),1,-1
          IOFF=(IZ-1)*MINGRI(1)*MINGRI(2)
          DO 210 IY=MINGRI(2),1,-1
            DO 220 IX=MINGRI(1),IY,-1
              KEYTYP(IOFF+(IY-1)*MINGRI(1)+IX)=KEYTYP(KML)
              ITGEOM(IOFF+(IY-1)*MINGRI(1)+IX)=ITGEOM(KML)
              IOF1=KML
              IOF2=IOFF+(IY-1)*MINGRI(1)+IX
              IF(IX .EQ. IY) THEN
                IKG=KEYTYP(IOF1)
                IKT(1)=ITGEOM(IOF1)
                IKT(2)=AXGTRS(IKT(1),3)
                WRITE(GEOCV,'(3A4)')
     >            CELLT(3*IKG-2),CELLT(3*IKG-1),CELLT(3*IKG)       
                IF(GEOCV .EQ. '            ') THEN 
                  IF(IPRT .GT. 10)        
     >              WRITE(IOUT,8000) NAMSBR,'X-Y',
     >              GEONAM,AXGTRN(IKT(1)),AXGTRN(IKT(2))
                ELSE
                  IF(IPRT .GT. 10)        
     >              WRITE(IOUT,8000) NAMSBR,'X-Y',
     >              GEOCV,AXGTRN(IKT(1)),AXGTRN(IKT(2))
                ENDIF
                VALSYM=LELCSY(IPGEOM,IPRT,GEONAM,GEOCV,NXYZ,IKT,
     >                        MESH,ISPLT,MESH1,ISPLT1)
                IF(.NOT. VALSYM) THEN 
                  WRITE(IOUT,8001) 'X-Y',GEOCV,
     >            AXGTRN(IKT(1)),AXGTRN(IKT(2))
                  CALL XABORT(NAMSBR//': INVALID SYMMETRY FOR CELL')
                ENDIF
              ENDIF
              KML=KML-1                             
 220        CONTINUE
 210      CONTINUE
 200    CONTINUE
        DO 230 IZ=1,MINGRI(3)
          IOFF=(IZ-1)*MINGRI(1)*MINGRI(2)
          DO 240 IY=1,MINGRI(2)
            DO 250 IX=1,IY-1
              KEYTYP(IOFF+(IY-1)*MINGRI(1)+IX)=
     >                 KEYTYP(IOFF+(IX-1)*MINGRI(2)+IY)
              ITGEOM(IOFF+(IY-1)*MINGRI(1)+IX)=
     >             AXGTRS(ITGEOM(IOFF+(IX-1)*MINGRI(2)+IY),3)
 250        CONTINUE
 240      CONTINUE
 230    CONTINUE
      ENDIF
      IF(KML .NE. 0) CALL XABORT(NAMSBR//': DATA ERROR')
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(MESH1,MESH)
      DEALLOCATE(ISPLT1,ISPLT)
*----
*  RETURN
*----
      RETURN 
*----
*  FORMAT
*----
 8000 FORMAT(1X,A6,' NOW TESTING SYMMETRY ',A3,' FOR ',
     > A12,1X,'WITH ROTATION',1X,A2,' AND ',A2)
 8001 FORMAT('            INVALID SYMMETRY ',A3,' FOR ',
     > A12,1X,'WITH ROTATION',1X,A2,' AND ',A2)
      END
