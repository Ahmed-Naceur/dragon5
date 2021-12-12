*DECK AXGSYM
      SUBROUTINE AXGSYM( IPGEOM,   IPRT, NBLOCK,  NTYPO,   NXYZ,
     >                   GEONAM, LCLSYM, MINGRI, MAXGRI,  CELLT,
     >                   KEYTYP, ITGEOM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Unfold assembly or cell according to center cell symmetry in
* $x$, $y$ or $z$ and verify if the symmetry is valid.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPGEOM  pointer to the reference geometry data structure.
* IPRT    intermediate printing level for output.
* NBLOCK  number of block in geometry.
* NTYPO   number of types in geometry.
* NXYZ    maximum mesh size in directions $x$, $y$ and $z$.
* GEONAM  name of the reference geometry.
* LCLSYM  flag that is set to 1 when the $x$ (LCLSYM(1)),
*         $y$ (LCLSYM(2)) and/or $z$ (LCLSYM(3)) 
*         symmetries are required. 
* MINGRI  minimum grid cell in $x$, $y$ and $z$ directions.
* MAXGRI  maximum grid cell in $x$, $y$ and $z$ directions.
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
      INTEGER            IPRT,NBLOCK,NTYPO,NXYZ
      CHARACTER          GEONAM*12
      INTEGER            LCLSYM(3)
      INTEGER            MINGRI(3),MAXGRI(3),CELLT(3*NTYPO)
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
      PARAMETER         (IOUT=6,NAMSBR='AXGSYM')
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISPLT,ISPLT1
      REAL, ALLOCATABLE, DIMENSION(:) :: MESH,MESH1
*----
*  LOCAL PARAMETERS
*----
      INTEGER            IX,IY,IZ,IOF1,IOF2
      INTEGER            IKOF1,IKOF2,ITOF1,ITOF2
      INTEGER            IKG,IKT(2)
      LOGICAL            VALSYM
      CHARACTER          GEOCV*12 
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ISPLT(NXYZ),ISPLT1(3*3*NXYZ))
      ALLOCATE(MESH(NXYZ+1),MESH1(2*3*3*(NXYZ+1)))
*----
*  APPLY SYMMETRY IN Z
*----
      IF( LCLSYM(3) .NE. 0) THEN
        IF(IPRT .GT. 10) THEN
          WRITE(IOUT,8000) NAMSBR,'Z-Z'
        ENDIF
        DO 200 IZ=1,MINGRI(3)
          DO 210 IY=1,MAXGRI(2)
            DO 220 IX=1,MAXGRI(1)
              IOF1=((IZ-1)*MAXGRI(2)+(IY-1))*MAXGRI(1)+IX
              IOF2=((MAXGRI(3)-IZ)*MAXGRI(2)+(IY-1))*MAXGRI(1)+IX 
              IKOF1=KEYTYP(IOF1)
              IKOF2=KEYTYP(IOF2)
              ITOF1=ITGEOM(IOF1)
              ITOF2=ITGEOM(IOF2)
              IF(IPRT .GT. 10) THEN
                WRITE(IOUT,8010) IZ,IY,IX,
     >           IOF1,IOF2,IKOF1,IKOF2,ITOF1,ITOF2
              ENDIF
              IF( IKOF1 .NE. IKOF2) THEN
                IF( IKOF1 .GT. IKOF2) THEN 
                  IKOF2= IKOF1
                  ITOF2= AXGTRS(ITOF1,4) 
                  KEYTYP(IOF2)= IKOF2
                  ITGEOM(IOF2)= ITOF2
                ELSE
                  IKOF1= IKOF2
                  ITOF1= AXGTRS(ITOF2,4) 
                  KEYTYP(IOF1)= IKOF1
                  ITGEOM(IOF1)= ITOF1
                ENDIF
              ENDIF
              IF(IKOF1 .GT. 0) THEN
                IF(IZ .EQ. (MAXGRI(3)+1-IZ)) THEN
                  IKG=IKOF1
                  IKT(1)=ITOF1
                  IKT(2)=AXGTRS(IKT(1),4)
                  IF(IKG .GT. 0) THEN 
                    WRITE(GEOCV,'(3A4)')
     >              CELLT(3*IKG-2),CELLT(3*IKG-1),CELLT(3*IKG)
                    IF(GEOCV .EQ. '            ') THEN       
                      IF(IPRT .GT. 10) WRITE(IOUT,8001) 
     >                GEONAM,AXGTRN(IKT(1)),AXGTRN(IKT(2))
                    ELSE
                      IF(IPRT .GT. 10) WRITE(IOUT,8001) 
     >                GEOCV,AXGTRN(IKT(1)),AXGTRN(IKT(2))
                    ENDIF
                    VALSYM=LELCSY(IPGEOM,IPRT,GEONAM,GEOCV,NXYZ,IKT,
     >                            MESH,ISPLT,MESH1,ISPLT1)
                    IF(.NOT. VALSYM) THEN 
                      WRITE(IOUT,8002) 'Z-Z',GEOCV,
     >                AXGTRN(IKT(1)),AXGTRN(IKT(2))
                      CALL XABORT(NAMSBR//
     >                ': INVALID Z SYMMETRY FOR CELL')
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
 220        CONTINUE
 210      CONTINUE
 200    CONTINUE
      ENDIF
*----
*  APPLY SYMMETRY IN Y
*----
      IF( LCLSYM(2).NE.0)THEN
        IF(IPRT .GT. 10) THEN
          WRITE(IOUT,8000) NAMSBR,'Y-Y'
        ENDIF
        DO 300 IZ=1,MAXGRI(3)
          DO 310 IY=1,MINGRI(2)
            DO 320 IX=1,MAXGRI(1)
              IOF1=((IZ-1)*MAXGRI(2)+(IY-1))*MAXGRI(1)+IX
              IOF2=((IZ-1)*MAXGRI(2)+(MAXGRI(2)-IY))*MAXGRI(1)+IX
              IKOF1=KEYTYP(IOF1)
              IKOF2=KEYTYP(IOF2)
              ITOF1=ITGEOM(IOF1)
              ITOF2=ITGEOM(IOF2)
              IF(IPRT .GT. 10) THEN
                WRITE(IOUT,8010) IZ,IY,IX,
     >           IOF1,IOF2,IKOF1,IKOF2,ITOF1,ITOF2
              ENDIF
              IF( IKOF1 .NE. IKOF2) THEN
                IF( IKOF1 .GT. IKOF2) THEN 
                  IKOF2= IKOF1
                  ITOF2= AXGTRS(ITOF1,2) 
                  KEYTYP(IOF2)= IKOF2
                  ITGEOM(IOF2)= ITOF2
                ELSE
                  IKOF1= IKOF2
                  ITOF1= AXGTRS(ITOF2,2) 
                  KEYTYP(IOF1)= IKOF1
                  ITGEOM(IOF1)= ITOF1
                ENDIF
              ENDIF
              IF(IKOF1 .GT. 0) THEN
                IF(IY .EQ. (MAXGRI(2)+1-IY) ) THEN
                  IKG=IKOF1
                  IKT(1)=ITOF1
                  IKT(2)=AXGTRS(IKT(1),2)
                  IF(IKG .GT. 0) THEN 
                    WRITE(GEOCV,'(3A4)')
     >              CELLT(3*IKG-2),CELLT(3*IKG-1),CELLT(3*IKG)       
                    IF(GEOCV .EQ. '            ') THEN       
                      IF(IPRT .GT. 10)  WRITE(IOUT,8001) 
     >                  GEONAM,AXGTRN(IKT(1)),AXGTRN(IKT(2))
                    ELSE
                      IF(IPRT .GT. 10)  WRITE(IOUT,8001) 
     >                  GEOCV,AXGTRN(IKT(1)),AXGTRN(IKT(2))
                    ENDIF
                    VALSYM=LELCSY(IPGEOM,IPRT,GEONAM,GEOCV,NXYZ,IKT,
     >                            MESH,ISPLT,MESH1,ISPLT1)
                    IF(.NOT. VALSYM) THEN 
                      WRITE(IOUT,8002) 'Y-Y',GEOCV,
     >                AXGTRN(IKT(1)),AXGTRN(IKT(2))
                      CALL XABORT(NAMSBR//
     >                ': INVALID Y SYMMETRY FOR CELL')
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
 320        CONTINUE
 310      CONTINUE
 300    CONTINUE
      ENDIF
*----
*  APPLY SYMMETRY IN X
*----
      IF( LCLSYM(1).NE.0)THEN
        IF(IPRT .GT. 10) THEN
          WRITE(IOUT,8000) NAMSBR,'X-X'
        ENDIF
        DO 400 IZ=1,MAXGRI(3)
          DO 410 IY=1,MAXGRI(2)
            DO 420 IX=1,MINGRI(1)
              IOF1=((IZ-1)*MAXGRI(2)+(IY-1))*MAXGRI(1)+IX
              IOF2=((IZ-1)*MAXGRI(2)+(IY-1))*MAXGRI(1)+MAXGRI(1)+1-IX
              IKOF1=KEYTYP(IOF1)
              IKOF2=KEYTYP(IOF2)
              ITOF1=ITGEOM(IOF1)
              ITOF2=ITGEOM(IOF2)
              IF(IPRT .GT. 10) THEN
                WRITE(IOUT,8010) IZ,IY,IX,
     >           IOF1,IOF2,IKOF1,IKOF2,ITOF1,ITOF2
              ENDIF
              IF( IKOF1 .NE. IKOF2) THEN
                IF( IKOF1 .GT. IKOF2) THEN 
                  IKOF2= IKOF1
                  ITOF2= AXGTRS(ITOF1,1) 
                  KEYTYP(IOF2)= IKOF2
                  ITGEOM(IOF2)= ITOF2
                ELSE
                  IKOF1= IKOF2
                  ITOF1= AXGTRS(ITOF2,1) 
                  KEYTYP(IOF1)= IKOF1
                  ITGEOM(IOF1)= ITOF1
                ENDIF
              ENDIF
              IF(IKOF1 .GT. 0) THEN
                IF(IX .EQ. (MAXGRI(1)+1-IX)) THEN
                  IKG=IKOF1
                  IKT(1)=ITOF1
                  IKT(2)=AXGTRS(IKT(1),1)
                  IF(IKG .GT. 0) THEN 
                    WRITE(GEOCV,'(3A4)')
     >              CELLT(3*IKG-2),CELLT(3*IKG-1),CELLT(3*IKG)       
                    IF(GEOCV .EQ. '            ') THEN
                      IF(IPRT .GT. 10) WRITE(IOUT,8001)       
     >                GEONAM,AXGTRN(IKT(1)),AXGTRN(IKT(2))
                    ELSE
                      IF(IPRT .GT. 10) WRITE(IOUT,8001)       
     >                GEOCV,AXGTRN(IKT(1)),AXGTRN(IKT(2))
                    ENDIF
                    VALSYM=LELCSY(IPGEOM,IPRT,GEONAM,GEOCV,NXYZ,IKT,
     >                            MESH,ISPLT,MESH1,ISPLT1)
                    IF(.NOT. VALSYM) THEN 
                      WRITE(IOUT,8002) 'X-X',GEOCV,
     >                AXGTRN(IKT(1)),AXGTRN(IKT(2))
                      CALL XABORT(NAMSBR//
     >                ': INVALID X SYMMETRY FOR CELL')
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
 420        CONTINUE
 410      CONTINUE
 400    CONTINUE
      ENDIF
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
 8000 FORMAT(1X,A6,' NOW TESTING SYMMETRY ',A3)
 8001 FORMAT(7X,A12,1X,'WITH ROTATION',1X,A2,' AND ',A2)
 8002 FORMAT('            INVALID SYMMETRY ',A3,' FOR ',
     > A12,1X,'WITH ROTATION',1X,A2,' AND ',A2)
 8010 FORMAT(1X,'IZ=',I6,1X,'IY=',I6,1X,'IX=',I6/
     >       1X,'IOF1 =',I6,1X,'IOF2 =',I6,
     >       1X,'KOF1 =',I6,1X,'KOF2 =',I6,
     >       1X,'TOF1 =',I6,1X,'TOF2 =',I6)
      END
