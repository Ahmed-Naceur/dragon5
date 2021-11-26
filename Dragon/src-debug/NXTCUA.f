*DECK NXTCUA
      SUBROUTINE NXTCUA(IPRINT,NDIM  ,IDIAG ,ISAXIS,
     >                  NBOCEL,NBUCEL,NOCELL,NUCELL,
     >                  ITSYM ,IDFEX ,IDFRT ,IUNFLD)
*
*----------
*
*Purpose:
* To create the array for testing the geometry in
* a Cartesian assembly for internal symmetries and unfolding
* the assembly according to the symmetries.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau.
*
*Parameters: input
* IPRINT  print level.
* NDIM    problem dimensions.
* IDIAG   the diagonal symmetry flag where:
*         =-1 indicates X- Y+ DIAG symmetry;
*         = 1 indicates X+ Y- DIAG symmetry;
*         = 0 indicates no DIAG symmetry.
* ISAXIS  symmetry vector for each direction.
* NBOCEL  number of cells in original geometry.
* NBUCEL  number of cells in unfolded geometry.
* NOCELL  number of cell before unfolding in
*         $X$, $Y$ and $Z$ directions.
* NUCELL  number of cell after unfolding in
*         $X$, $Y$ and $Z$ directions.
*
*Parameters: output
* ITSYM   array to identify the symmetry to test for each original
*         cell where:
*         ITSYM(1,*) identify $X$ symmetry;
*         ITSYM(2,*) identify $Y$ symmetry;
*         ITSYM(3,*) identify $Z$ symmetry;
*         ITSYM(4,*) identify $X-Y$ symmetry.
*         A value of 0 indicate that the geometry does not need
*         to be verified while a value of 1 implies a verification
*         of the geometry.
* IDFEX   identify faces associated with external boundary for a
*         generating cell and number of times this cell is used. Here:
*         IDFEX( 1,*)  identify bottom $X$ hexagonal face;
*         IDFEX( 2,*)  identify top $X$ hexagonal face;
*         IDFEX( 3,*)  identify bottom $Y$ hexagonal face;
*         IDFEX( 4,*)  identify top $Y$ hexagonal face;
*         IDFEX( 5,*)  identify bottom $Z$ face;
*         IDFEX( 6,*)  identify top $Z$ face;
*         IDFEX( 7,*)  not used;
*         IDFEX( 8,*)  not used;
*         IDFEX( 9,*)  not used;
*         IDFEX(10,*)  not used.
* IDFRT   identify reflection/transmission faces.
* IUNFLD  array to identify the generating cell (IUNFLD(1,*))
*         and the rotation associated with this region in space.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,NDIM,IDIAG,ISAXIS(3)
      INTEGER          NBUCEL,NOCELL(3),NUCELL(3)
      INTEGER          NBOCEL,ITSYM(4,NBOCEL),IDFEX(0:10,NBOCEL),
     >                 IDFRT(8,NBOCEL),IUNFLD(2,NBUCEL)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTCUA')
*----
*  Functions
*----
      INTEGER          NXTTRS
*----
*  Local variables
*----
      INTEGER          IDIR,NSCELL(3),NFCELL(3),
     >                 NSUC(3),
     >                 IGEN,IGENT,IX,IY,IZ,ILOCD,
     >                 ILOCR,IOFYZ,IOFYZR,IOFZ,IOFZR
*----
*  Data
*----
      CHARACTER*2      CTRN(24)
      SAVE             CTRN
      DATA             CTRN
     > /'+A','+B','+C','+D','+E','+F','+G','+H','+I','+J','+K','+L',
     >  '-A','-B','-C','-D','-E','-F','-G','-H','-I','-J','-K','-L'/
*----
*  Processing starts:
*  print routine opening header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      CALL XDISET(IDFEX ,11*NBOCEL,0)
      CALL XDISET(IDFRT ,8*NBOCEL,0)
      CALL XDISET(ITSYM ,4*NBOCEL,0)
*----
*  Prepare direction control vector for
*  original assembly
*----
      DO IDIR=1,3
        NSCELL(IDIR)=MAX(1,NOCELL(IDIR))
        NSUC(IDIR)=MAX(1,NUCELL(IDIR))
        IF(ISAXIS(IDIR) .EQ. -2) THEN
          NFCELL(IDIR)=NOCELL(IDIR)
        ELSE IF(ISAXIS(IDIR) .EQ. -1) THEN
          NFCELL(IDIR)=NOCELL(IDIR)-1
        ELSE
          NFCELL(IDIR)=0
        ENDIF
*        WRITE(6,*) 'NSCELL =',IDIR,NSCELL(IDIR),NFCELL(IDIR)
      ENDDO
*----
*  Position original cell and process diagonal
*  symmetry
*----
      IGEN=0
      IF(IDIAG .EQ. -1) THEN
*----
*  Process X- Y+ diagonal symmetry
*----
        DO IZ=0,NSCELL(3)-1
          IOFZ=(IZ+NFCELL(3))*NUCELL(2)
          DO IY=0,NSCELL(2)-1
            IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
            IX=IY
            IGEN=IGEN+1
            IF(IGEN .GT. NBOCEL) CALL XABORT(NAMSBR//
     >      ': Cell number exceeds number of cells permitted')
            ILOCD=IOFYZ+IX+NFCELL(1)+1
            IUNFLD(1,ILOCD)=IGEN
            IUNFLD(2,ILOCD)=1
            ITSYM(4,IGEN)=IDIAG
            DO IX=IY+1,NSCELL(1)-1
              IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
              IOFYZR=(IX+NFCELL(2)+IOFZ)*NUCELL(1)
              IGEN=IGEN+1
              IF(IGEN .GT. NBOCEL) CALL XABORT(NAMSBR//
     >        ': Cell number exceeds number of cells permitted')
              ILOCD=IOFYZ+IX+NFCELL(1)+1
              ILOCR=IOFYZR+IY+NFCELL(1)+1
              IUNFLD(1,ILOCD)=IGEN
              IUNFLD(2,ILOCD)=1
              IUNFLD(1,ILOCR)=IGEN
              IUNFLD(2,ILOCR)=NXTTRS(1,2)
            ENDDO
          ENDDO
        ENDDO
*----
*  Identify cells to tests for X
*  reflection symmetry
*----
        IF(ISAXIS(1) .EQ. 1) THEN
          IX=NSCELL(1)-1
          DO IZ=0,NSCELL(3)-1
            IOFZ=(IZ+NFCELL(3))*NUCELL(2)
            DO IY=0,NSCELL(2)-1
              IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
              ILOCD=IOFYZ+IX+NFCELL(1)+1
              IGEN=IUNFLD(1,ILOCD)
              ITSYM(1,IGEN)=ISAXIS(1)
            ENDDO
          ENDDO
        ENDIF
*----
*  Identify cells to tests for Y
*  reflection symmetry
*----
        IF(ISAXIS(2) .EQ. -1) THEN
          IY=0
          DO IZ=0,NSCELL(3)-1
            IOFZ=(IZ+NFCELL(3))*NUCELL(2)
            IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
            DO IX=0,NSCELL(1)-1
              ILOCD=IOFYZ+IX+NFCELL(1)+1
              IGEN=IUNFLD(1,ILOCD)
              ITSYM(2,IGEN)=ISAXIS(2)
            ENDDO
          ENDDO
        ENDIF
*----
*  Identify cells to tests for Z
*  reflection symmetry
*----
        IF(ABS(ISAXIS(3)) .EQ. 1) THEN
          IF(ISAXIS(3) .EQ. -1) THEN
            IZ=0
          ELSE IF(ISAXIS(3) .EQ. 1) THEN
            IZ=NSCELL(3)-1
          ENDIF
          IOFZ=(IZ+NFCELL(3))*NUCELL(2)
          DO IY=0,NSCELL(2)-1
            IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
            DO IX=IY,NSCELL(1)-1
              IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
              ILOCD=IOFYZ+IX+NFCELL(1)+1
              IGEN=IUNFLD(1,ILOCD)
              ITSYM(3,IGEN)=ISAXIS(3)
            ENDDO
          ENDDO
        ENDIF
      ELSE IF(IDIAG .EQ. 1) THEN
*----
*  Process X+ Y- diagonal symmetry
*----
        DO IZ=0,NSCELL(3)-1
          IOFZ=(IZ+NFCELL(3))*NUCELL(2)
          DO IY=0,NSCELL(2)-1
            IGENT=IGEN+IY+2
            IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
            IX=IY
            IGEN=IGEN+1
            IF(IGEN .GT. NBOCEL) CALL XABORT(NAMSBR//
     >      ': Cell number exceeds number of cells permitted')
            IGENT=IGENT-1
            ILOCD=IOFYZ+IX+NFCELL(1)+1
            IUNFLD(1,ILOCD)=IGENT
            IUNFLD(2,ILOCD)=1
            ITSYM(4,IGENT)=IDIAG
            DO IX=IY-1,0,-1
              IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
              IOFYZR=(IX+NFCELL(2)+IOFZ)*NUCELL(1)
              IGEN=IGEN+1
              IF(IGEN .GT. NBOCEL) CALL XABORT(NAMSBR//
     >        ': Cell number exceeds number of cells permitted')
              IGENT=IGENT-1
              ILOCD=IOFYZ+IX+NFCELL(1)+1
              ILOCR=IOFYZR+IY+NFCELL(1)+1
              IUNFLD(1,ILOCD)=IGENT
              IUNFLD(2,ILOCD)=1
              IUNFLD(1,ILOCR)=IGENT
              IUNFLD(2,ILOCR)=NXTTRS(1,2)
            ENDDO
          ENDDO
        ENDDO
*----
*  Identify cells to tests for X
*  reflection symmetry
*----
        IF(ISAXIS(1) .EQ. -1) THEN
          IX=0
          DO IZ=0,NSCELL(3)-1
            IOFZ=(IZ+NFCELL(3))*NUCELL(2)
            DO IY=0,NSCELL(2)-1
              IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
              ILOCD=IOFYZ+IX+NFCELL(1)+1
              IGEN=IUNFLD(1,ILOCD)
              ITSYM(1,IGEN)=ISAXIS(1)
            ENDDO
          ENDDO
        ENDIF
*----
*  Identify cells to tests for Y
*  reflection symmetry
*----
        IF(ISAXIS(2) .EQ. 1) THEN
          IY=NSCELL(2)-1
          DO IZ=0,NSCELL(3)-1
            IOFZ=(IZ+NFCELL(3))*NUCELL(2)
            IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
            DO IX=IY,0,-1
              IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
              ILOCD=IOFYZ+IX+NFCELL(1)+1
              IGEN=IUNFLD(1,ILOCD)
              ITSYM(2,IGEN)=ISAXIS(2)
            ENDDO
          ENDDO
        ENDIF
*----
*  Identify cells to tests for Z
*  reflection symmetry
*----
        IF(ABS(ISAXIS(3)) .EQ. 1) THEN
          IF(ISAXIS(3) .EQ. -1) THEN
            IZ=0
          ELSE IF(ISAXIS(3) .EQ. 1) THEN
            IZ=NSCELL(3)-1
          ENDIF
          IOFZ=(IZ+NFCELL(3))*NUCELL(2)
          DO IY=0,NSCELL(2)-1
            IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
            DO IX=IY,0,-1
              IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
              ILOCD=IOFYZ+IX+NFCELL(1)+1
              IGEN=IUNFLD(1,ILOCD)
              ITSYM(3,IGEN)=ISAXIS(3)
            ENDDO
          ENDDO
        ENDIF
      ELSE
        DO IZ=0,NSCELL(3)-1
          IOFZ=(IZ+NFCELL(3))*NUCELL(2)
          DO IY=0,NSCELL(2)-1
            IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
            DO IX=0,NSCELL(1)-1
              IGEN=IGEN+1
              IF(IGEN .GT. NBOCEL) CALL XABORT(NAMSBR//
     >        ': Cell number exceeds number of cells permitted')
              ILOCD=IOFYZ+IX+NFCELL(1)+1
              IUNFLD(1,ILOCD)=IGEN
              IUNFLD(2,ILOCD)=1
            ENDDO
          ENDDO
        ENDDO
*----
*  Identify cells to tests for X
*  reflection symmetry
*----
        IF(ABS(ISAXIS(1)) .EQ. 1) THEN
          IF(ISAXIS(1) .EQ. -1) THEN
            IX=0
          ELSE IF(ISAXIS(1) .EQ. 1) THEN
            IX=NSCELL(1)-1
          ENDIF
          DO IZ=0,NSCELL(3)-1
            IOFZ=(IZ+NFCELL(3))*NUCELL(2)
            DO IY=0,NSCELL(2)-1
              IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
              ILOCD=IOFYZ+IX+NFCELL(1)+1
              IGEN=IUNFLD(1,ILOCD)
              ITSYM(1,IGEN)=ISAXIS(1)
            ENDDO
          ENDDO
        ENDIF
*----
*  Identify cells to tests for Y
*  reflection symmetry
*----
        IF(ABS(ISAXIS(2)) .EQ. 1) THEN
          IF(ISAXIS(2) .EQ. -1) THEN
            IY=0
          ELSE IF(ISAXIS(2) .EQ. 1) THEN
            IY=NSCELL(2)-1
          ENDIF
          DO IZ=0,NSCELL(3)-1
            IOFZ=(IZ+NFCELL(3))*NUCELL(2)
            IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
            DO IX=0,NSCELL(1)-1
              ILOCD=IOFYZ+IX+NFCELL(1)+1
              IGEN=IUNFLD(1,ILOCD)
              ITSYM(2,IGEN)=ISAXIS(2)
            ENDDO
          ENDDO
        ENDIF
*----
*  Identify cells to tests for Z
*  reflection symmetry
*----
        IF(ABS(ISAXIS(3)) .EQ. 1) THEN
          IF(ISAXIS(3) .EQ. -1) THEN
            IZ=0
          ELSE IF(ISAXIS(3) .EQ. 1) THEN
            IZ=NSCELL(3)-1
          ENDIF
          IOFZ=(IZ+NFCELL(3))*NUCELL(2)
          DO IY=0,NSCELL(2)-1
            IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
            DO IX=0,NSCELL(1)-1
              ILOCD=IOFYZ+IX+NFCELL(1)+1
              IGEN=IUNFLD(1,ILOCD)
              ITSYM(3,IGEN)=ISAXIS(3)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      DO ILOCD=1,NBOCEL
        IF(ITSYM(4,ILOCD) .EQ. -1) THEN
          IF(ITSYM(2,ILOCD) .EQ. -1) THEN
            ITSYM(1,ILOCD)=1
          ELSE IF(ITSYM(1,ILOCD) .EQ. 1) THEN
            ITSYM(2,ILOCD)=-1
          ENDIF
        ELSE IF(ITSYM(4,ILOCD) .EQ. 1) THEN
          IF(ITSYM(2,ILOCD) .EQ. 1) THEN
            ITSYM(1,ILOCD)=-1
          ELSE IF(ITSYM(1,ILOCD) .EQ. -1) THEN
            ITSYM(2,ILOCD)=1
          ENDIF
        ENDIF
      ENDDO
      IF(ISAXIS(1) .NE. 0) THEN
        IF(ISAXIS(1) .EQ. -2) THEN
*----
*  SSYM X-
*  Fill position IXR=1,NSCELL(1) with cells at
*  position IXD=NUCELL(1)-IXR+1
*----
          DO IZ=0,NSCELL(3)-1
            IOFZ=(IZ+NFCELL(3))*NUCELL(2)
            DO IY=0,NSCELL(2)-1
              IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
              DO IX=1,NSCELL(1)
                ILOCD=IOFYZ+NUCELL(1)+1-IX
                ILOCR=IOFYZ+IX
                IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
                IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),1)
              ENDDO
            ENDDO
          ENDDO
        ELSE IF(ISAXIS(1) .EQ. -1) THEN
*----
*  SYME X-
*  Fill position IXR=1,NSCELL(1)-1 with cells at
*  position IXD=NUCELL(1)-IXR+1
*  set test flag for IX=NSCELL(1)
*----
          DO IZ=0,NSCELL(3)-1
            IOFZ=(IZ+NFCELL(3))*NUCELL(2)
            DO IY=0,NSCELL(2)-1
              IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
              DO IX=1,NSCELL(1)-1
                ILOCD=IOFYZ+NUCELL(1)+1-IX
                ILOCR=IOFYZ+IX
                IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
                IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),1)
              ENDDO
            ENDDO
          ENDDO
        ELSE IF(ISAXIS(1) .EQ. 1) THEN
*----
*  SYME X+
*  Fill position IXR=NUCELL(1)-IXD+1 with cell
*  at position IXD=1,NSCELL(1)-1
*  set test flag for IX=NSCELL(1)
*----
          DO IZ=0,NSCELL(3)-1
            IOFZ=(IZ+NFCELL(3))*NUCELL(2)
            DO IY=0,NSCELL(2)-1
              IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
              DO IX=1,NSCELL(1)-1
                ILOCD=IOFYZ+IX
                ILOCR=IOFYZ+NUCELL(1)-IX+1
                IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
                IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),1)
              ENDDO
            ENDDO
          ENDDO
        ELSE IF(ISAXIS(1) .EQ. 2) THEN
*----
*  SSYM X+
*  Fill position IXR=NUCELL(1)-IXD+1 with cell
*  at position IXD=1,NSCELL(1)
*----
          DO IZ=0,NSCELL(3)-1
            IOFZ=(IZ+NFCELL(3))*NUCELL(2)
            DO IY=0,NSCELL(2)-1
              IOFYZ=(IY+NFCELL(2)+IOFZ)*NUCELL(1)
              DO IX=1,NSCELL(1)
                ILOCD=IOFYZ+IX
                ILOCR=IOFYZ+NUCELL(1)-IX+1
                IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
                IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),1)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        NSCELL(1)=MAX(1,NUCELL(1))
      ENDIF
      IF(ISAXIS(2) .NE. 0) THEN
        IF(ISAXIS(2) .EQ. -2) THEN
*----
*  SSYM Y-
*  Fill position IYR=1,NSCELL(2) with cells at
*  position IYD=NUCELL(2)-IYR+1
*----
          DO IZ=0,NSCELL(3)-1
            IOFZ=(IZ+NFCELL(3))*NUCELL(2)
            DO IY=1,NSCELL(2)
              IOFYZ=(NUCELL(2)-IY+IOFZ)*NUCELL(1)
              IOFYZR=(IY-1+IOFZ)*NUCELL(1)
              DO IX=1,NSCELL(1)
                ILOCD=IOFYZ+IX
                ILOCR=IOFYZR+IX
                IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
                IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),3)
              ENDDO
            ENDDO
          ENDDO
        ELSE IF(ISAXIS(2) .EQ. -1) THEN
*----
*  SYME Y-
*  Fill position IYR=1,NSCELL(2)-1 with cells at
*  position IYD=NUCELL(2)-IYR+1
*  set test flag for IY=NSCELL(2)
*----
          DO IZ=0,NSCELL(3)-1
            IOFZ=(IZ+NFCELL(3))*NUCELL(2)
            DO IY=1,NSCELL(2)-1
              IOFYZ=(NUCELL(2)-IY+IOFZ)*NUCELL(1)
              IOFYZR=(IY-1+IOFZ)*NUCELL(1)
              DO IX=1,NSCELL(1)
                ILOCD=IOFYZ+IX
                ILOCR=IOFYZR+IX
                IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
                IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),3)
              ENDDO
            ENDDO
          ENDDO
        ELSE IF(ISAXIS(2) .EQ. 1) THEN
*----
*  SYME Y+
*  Fill position IYR=NUCELL(2)-IYD+1 with cell
*  at position IYD=1,NSCELL(2)-1
*  set test flag for IY=NSCELL(2)
*----
          DO IZ=0,NSCELL(3)-1
            IOFZ=(IZ+NFCELL(3))*NUCELL(2)
            DO IY=1,NSCELL(2)-1
              IOFYZR=(NUCELL(2)-IY+IOFZ)*NUCELL(1)
              IOFYZ=(IY-1+IOFZ)*NUCELL(1)
              DO IX=1,NSCELL(1)
                ILOCD=IOFYZ+IX
                ILOCR=IOFYZR+IX
                IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
                IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),3)
              ENDDO
            ENDDO
          ENDDO
        ELSE IF(ISAXIS(2) .EQ. 2) THEN
*----
*  SSYM Y+
*  Fill position IYR=NUCELL(2)-IYD+1 with cell
*  at position IYD=1,NSCELL(2)
*----
          DO IZ=0,NSCELL(3)-1
            IOFZ=(IZ+NFCELL(3))*NUCELL(2)
            DO IY=1,NSCELL(2)
              IOFYZR=(NUCELL(2)-IY+IOFZ)*NUCELL(1)
              IOFYZ=(IY-1+IOFZ)*NUCELL(1)
              DO IX=1,NSCELL(1)
                ILOCD=IOFYZ+IX
                ILOCR=IOFYZR+IX
                IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
                IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),3)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        NSCELL(2)=MAX(1,NUCELL(2))
      ENDIF
      IF(ISAXIS(3) .NE. 0) THEN
        IF(ISAXIS(3) .EQ. -2) THEN
*----
*  SSYM Z-
*  Fill position IZR=1,NSCELL(3) with cells at
*  position IZD=NSUC(3)-IZR+1
*----
          DO IZ=1,NSCELL(3)
            IOFZR=(IZ-1)*NUCELL(2)
            IOFZ=(NSUC(3)-IZ)*NUCELL(2)
            DO IY=1,NSCELL(2)
              IOFYZ=(IY-1+IOFZ)*NUCELL(1)
              IOFYZR=(IY-1+IOFZR)*NUCELL(1)
              DO IX=1,NSCELL(1)
                ILOCD=IOFYZ+IX
                ILOCR=IOFYZR+IX
                IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
                IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),-1)
              ENDDO
            ENDDO
          ENDDO
        ELSE IF(ISAXIS(3) .EQ. -1) THEN
*----
*  SYME Z-
*  Fill position IZR=1,NSCELL(3)-1 with cells at
*  position IZD=NSUC(3)-IZR+1
*  set test flag for IZ=NSCELL(3)
*----
          DO IZ=1,NSCELL(3)-1
            IOFZR=(IZ-1)*NUCELL(2)
            IOFZ=(NSUC(3)-IZ)*NUCELL(2)
            DO IY=1,NSCELL(2)
              IOFYZ=(IY-1+IOFZ)*NUCELL(1)
              IOFYZR=(IY-1+IOFZR)*NUCELL(1)
              DO IX=1,NSCELL(1)
                ILOCD=IOFYZ+IX
                ILOCR=IOFYZR+IX
                IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
                IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),-1)
              ENDDO
            ENDDO
          ENDDO
        ELSE IF(ISAXIS(3) .EQ. 1) THEN
*----
*  SYME Z+
*  Fill position IZR=NSUC(3)-IZD+1 with cell
*  at position IZD=1,NSCELL(3)-1
*  set test flag for IZ=NSCELL(3)
*----
          DO IZ=1,NSCELL(3)-1
            IOFZ=(IZ-1)*NUCELL(2)
            IOFZR=(NSUC(3)-IZ)*NUCELL(2)
            DO IY=1,NSCELL(2)
              IOFYZ=(IY-1+IOFZ)*NUCELL(1)
              IOFYZR=(IY-1+IOFZR)*NUCELL(1)
              DO IX=1,NSCELL(1)
                ILOCD=IOFYZ+IX
                ILOCR=IOFYZR+IX
                IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
                IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),-1)
              ENDDO
            ENDDO
          ENDDO
        ELSE IF(ISAXIS(3) .EQ. 2) THEN
*----
*  SSYM Z+
*  Fill position IZR=NSUC(3)-IZD+1 with cell
*  at position IZD=1,NSCELL(3)
*----
          DO IZ=1,NSCELL(3)
            IOFZ=(IZ-1)*NUCELL(2)
            IOFZR=(NSUC(3)-IZ)*NUCELL(2)
            DO IY=1,NSCELL(2)
              IOFYZ=(IY-1+IOFZ)*NUCELL(1)
              IOFYZR=(IY-1+IOFZR)*NUCELL(1)
              DO IX=1,NSCELL(1)
                ILOCD=IOFYZ+IX
                ILOCR=IOFYZR+IX
                IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
                IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),-1)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        NSCELL(3)=MAX(1,NSUC(3))
      ENDIF
*----
*  Localize external faces
*  1. X- (1) AND X+ (2)
*----
      DO IZ=1,NSUC(3)
        DO IY=1,NUCELL(2)
*          write(6,*) ' X faces =',IY,IZ
          ILOCD=NUCELL(1)*(IY-1+NUCELL(2)*(IZ-1))+1
          IF(IUNFLD(2,ILOCD) .EQ. 1) THEN
            IGEN=IUNFLD(1,ILOCD)
            IDFEX(1,IGEN)=1
            IDFRT(1,IGEN)=IGEN
          ENDIF
*          write(6,*) ' X - =',ILOCD,IUNFLD(2,ILOCD),IDFEX(1,IGEN)
          ILOCD=NUCELL(1)*(IY-1+NUCELL(2)*(IZ-1))+NUCELL(1)
          IF(IUNFLD(2,ILOCD) .EQ. 1) THEN
            IGEN=IUNFLD(1,ILOCD)
            IDFEX(2,IGEN)=1
            IDFRT(2,IGEN)=IGEN
          ENDIF
*          write(6,*) ' X + =',ILOCD,IUNFLD(2,ILOCD),IDFEX(2,IGEN)
        ENDDO
      ENDDO
*----
*  2. Y- (3) Y+ (4)
*----
      DO IZ=1,NSUC(3)
        DO IX=1,NUCELL(1)
*          write(6,*) ' Y faces =',IX,IZ
          ILOCD=NUCELL(1)*NUCELL(2)*(IZ-1)+IX
          IF(IUNFLD(2,ILOCD) .EQ. 1) THEN
            IGEN=IUNFLD(1,ILOCD)
            IDFEX(3,IGEN)=1
            IDFRT(3,IGEN)=IGEN
          ENDIF
*          write(6,*) ' Y - =',ILOCD,IUNFLD(2,ILOCD),IDFEX(3,IGEN)
          ILOCD=NUCELL(1)*(NUCELL(2)-1+NUCELL(2)*(IZ-1))+IX
          IF(IUNFLD(2,ILOCD) .EQ. 1) THEN
            IGEN=IUNFLD(1,ILOCD)
            IDFEX(4,IGEN)=1
            IDFRT(4,IGEN)=IGEN
          ENDIF
*          write(6,*) ' Y + =',ILOCD,IUNFLD(2,ILOCD),IDFEX(4,IGEN)
        ENDDO
      ENDDO
*----
*  3. Z- (5) Z+ (6)
*----
      DO IY=1,NUCELL(2)
        DO IX=1,NUCELL(1)
*          write(6,*) ' Z faces =',IX,IY
          ILOCD=NUCELL(1)*(IY-1)+IX
          IF(IUNFLD(2,ILOCD) .EQ. 1) THEN
            IGEN=IUNFLD(1,ILOCD)
            IDFEX(5,IGEN)=1
            IDFRT(5,IGEN)=IGEN
          ENDIF
*          write(6,*) ' Z - =',ILOCD,IUNFLD(2,ILOCD),IDFEX(5,IGEN)
          ILOCD=NUCELL(1)*(IY-1+NUCELL(2)*(NSUC(3)-1))+IX
          IF(IUNFLD(2,ILOCD) .EQ. 1) THEN
            IGEN=IUNFLD(1,ILOCD)
            IDFEX(6,IGEN)=1
            IDFRT(6,IGEN)=IGEN
          ENDIF
*          write(6,*) ' Z + =',ILOCD,IUNFLD(2,ILOCD),IDFEX(6,IGEN)
        ENDDO
      ENDDO
*----
*  For translation BC, find translated cell
*  1. X translation
*----
      IF(ISAXIS(1) .EQ. 3) THEN
        DO IGEN=1,NBOCEL
          IF(IDFEX(1,IGEN) .EQ. 1) THEN
            DO IZ=1,NSUC(3)
              DO IY=1,NUCELL(2)
                ILOCD=NUCELL(1)*(IY-1+NUCELL(2)*(IZ-1))+1
                IF(IUNFLD(2,ILOCD) .EQ. 1) THEN
                  IF(IUNFLD(1,ILOCD) .EQ. IGEN) THEN
                    ILOCD=NUCELL(1)*(IY-1+NUCELL(2)*(IZ-1))+NUCELL(1)
                    IGENT=IUNFLD(1,ILOCD)
                    IDFRT(1,IGEN)=IGENT
                    IDFRT(2,IGENT)=IGEN
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF
*----
*  2. Y translation
*----
      IF(ISAXIS(2) .EQ. 3) THEN
        DO IGEN=1,NBOCEL
          IF(IDFEX(3,IGEN) .EQ. 1) THEN
            DO IZ=1,NSUC(3)
              DO IX=1,NUCELL(1)
                ILOCD=NUCELL(1)*NUCELL(2)*(IZ-1)+IX
                IF(IUNFLD(2,ILOCD) .EQ. 1) THEN
                  IF(IUNFLD(1,ILOCD) .EQ. IGEN) THEN
                    ILOCD=NUCELL(1)*(NUCELL(2)-1+NUCELL(2)*(IZ-1))+IX
                    IGENT=IUNFLD(1,ILOCD)
                    IDFRT(3,IGEN)=IGENT
                    IDFRT(4,IGENT)=IGEN
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF
*----
*  3. Z translation
*----
      IF(ISAXIS(3) .EQ. 3) THEN
        DO IGEN=1,NBOCEL
          IF(IDFEX(5,IGEN) .EQ. 1) THEN
            DO IY=1,NUCELL(2)
              DO IX=1,NUCELL(1)
                ILOCD=NUCELL(1)*(IY-1)+IX
                IF(IUNFLD(2,ILOCD) .EQ. 1) THEN
                  IF(IUNFLD(1,ILOCD) .EQ. IGEN) THEN
                    ILOCD=NUCELL(1)*(IY-1+NUCELL(2)*(NSUC(3)-1))+IX
                    IGENT=IUNFLD(1,ILOCD)
                    IDFRT(5,IGEN)=IGENT
                    IDFRT(6,IGENT)=IGEN
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF
*----
*  Analyze translation boundary
*----
*----
*  Compute the number of times each cell appears
*  after unfolding
*----
      DO IGEN=1,NBOCEL
        DO ILOCD=1,NBUCEL
          IF(ABS(IUNFLD(1,ILOCD)) .EQ. IGEN) THEN
            IDFEX(0,IGEN)=IDFEX(0,IGEN)+1
          ENDIF
        ENDDO
      ENDDO
*----
*  For 2-D cases reset components 5, 6 of IDFEX to 0
*----
      IF(NDIM .EQ. 2) THEN
        DO IGEN=1,NBOCEL
          IDFEX(5,IGEN)=0
          IDFEX(6,IGEN)=0
        ENDDO
      ENDIF
*----
*  Processing finished:
*  print routine output and closing header if required
*  and return
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6002)
        DO IZ=1,NSUC(3)
          IF(NDIM .EQ. 3) THEN
            WRITE(IOUT,6003) IZ
          ENDIF
          WRITE(IOUT,6004) (IX,IX=1,NUCELL(1))
          WRITE(IOUT,6012)
          DO IY=NUCELL(2),1,-1
            ILOCD=((IZ-1)*NUCELL(2)+(IY-1))*NUCELL(1)
            WRITE(IOUT,6005) IY,
     >        (IUNFLD(1,IX),IX=ILOCD+1,ILOCD+NUCELL(1))
          ENDDO
        ENDDO
        WRITE(IOUT,6008)
        DO IZ=1,NSUC(3)
          IF(NDIM .EQ. 3) THEN
            WRITE(IOUT,6003) IZ
          ENDIF
          WRITE(IOUT,6004) (IX,IX=1,NUCELL(1))
          WRITE(IOUT,6012)
          DO IY=NUCELL(2),1,-1
            ILOCD=((IZ-1)*NUCELL(2)+(IY-1))*NUCELL(1)
            WRITE(IOUT,6013) IY,
     >        (CTRN(IUNFLD(2,IX)),IX=ILOCD+1,ILOCD+NUCELL(1))
          ENDDO
        ENDDO
        WRITE(IOUT,6006)
        DO ILOCD=1,NBOCEL
          WRITE(IOUT,6007)
     >    ILOCD,(ITSYM(IX,ILOCD),IX=1,4)
        ENDDO
        WRITE(IOUT,6009)
        DO ILOCD=1,NBOCEL
          IF(NDIM .EQ. 3) THEN
            WRITE(IOUT,6007)
     >      ILOCD,(IDFEX(IX,ILOCD),IX=1,6),IDFEX(0,ILOCD)
          ELSE
            WRITE(IOUT,6011)
     >      ILOCD,(IDFEX(IX,ILOCD),IX=1,4),IDFEX(0,ILOCD)
          ENDIF
        ENDDO
        WRITE(IOUT,6010)
        DO ILOCD=1,NBOCEL
          WRITE(IOUT,6007)
     >    ILOCD,(IDFRT(IX,ILOCD),IX=1,2*NDIM)
        ENDDO
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6002 FORMAT(' Cells in assembly')
 6003 FORMAT(' Plane IZ =',I6)
 6004 FORMAT(' IX=',8X,24(I7,1X))
 6005 FORMAT(' IY=',I7,1X,24(I7.7,1X))
 6006 FORMAT(/' Symmetrized cell     X    Y    Z    D')
 6007 FORMAT(' Cell ',I7.7,5X,20I5)
 6008 FORMAT(/' Cell rotations in assembly')
 6009 FORMAT(/' External faces      -X   +X   -Y   +Y   -Z   +Z   ND')
 6010 FORMAT(/' Coupled  faces      -X   +X   -Y   +Y   -Z   +Z  ')
 6011 FORMAT(' Cell ',I7.7,5X,4I5,10X,14I5)
 6012 FORMAT('X/Y')
 6013 FORMAT(' IY=',I7,1X,24(5X,A2,1X))
      END
