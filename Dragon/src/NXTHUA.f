*DECK NXTHUA
      SUBROUTINE NXTHUA(IPRINT,NDIM  ,IHSYM ,ISAXIS,
     >                  NBOCEL,NBUCEL,NOCELL,NUCELL,
     >                  ITSYM ,IDFEX ,IDFRT ,IUNFLD)
*
*----------
*
*Purpose:
* To create the array for testing the geometry in
* an hexagonal assembly for internal symmetries and unfolding
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
* IHSYM   hexagonal symmetry option where:
*         = 0 geometry is not hexagonal;
*         = 1 for S30;
*         = 2 for SA60;
*         = 3 for SB60;
*         = 4 for S90;
*         = 5 for R120;
*         = 6 for R180;
*         = 7 for SA180;
*         = 8 for SB180;
*         = 9 for COMPLETE;
*         =10 for R60.
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
*         ITSYM(1,*) identify hexagonal symmetry;
*         ITSYM(2,*) not used;
*         ITSYM(3,*) identify $Z$ symmetry;
*         ITSYM(4,*) not used.
*         A value of 0 indicate that the geometry does not need
*         to be verified while a value of 1 implies a verification
*         of the geometry.
* IDFEX   identify faces associated with external boundary for a
*         generating cell and number of times this cell is used. Here:
*         IDFEX( 1,*)  identify bottom $U$ hexagonal face;
*         IDFEX( 2,*)  identify top $U$ hexagonal face;
*         IDFEX( 3,*)  identify bottom $V$ hexagonal face;
*         IDFEX( 4,*)  identify top $V$ hexagonal face;
*         IDFEX( 5,*)  identify bottom $Z$ face;
*         IDFEX( 6,*)  identify top $Z$ face;
*         IDFEX( 7,*)  not used;
*         IDFEX( 8,*)  not used;
*         IDFEX (9,*)  identify bottom $W$ hexagonal face;
*         IDFEX(10,*)  identify top $W$ hexagonal face.
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
      INTEGER          IPRINT,NDIM,IHSYM,ISAXIS(3)
      INTEGER          NBUCEL,NOCELL(3),NUCELL(3)
      INTEGER          NBOCEL,ITSYM(4,NBOCEL),IDFEX(0:10,NBOCEL),
     >                 IDFRT(8,NBOCEL),IUNFLD(2,NBUCEL)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTHUA')
*----
*  Functions
*----
      INTEGER          NXTTRS
*----
*  Local variables
*----
      INTEGER          IDIR,NSCELL(3),NFCELL(3),
     >                 NSUC(3),ID1,ID2,ID3,ISECT,
     >                 IGEN,IGENT,IX,IZ,ILOCD,
     >                 ILOCR,IOFZ,IOFZR,NCR,NCILC,IDD,ICR
      DOUBLE PRECISION ARGS
*----
*  Data
*----
      CHARACTER*2      CTRN(24)
      INTEGER          IDSEC(6)
      SAVE             CTRN,IDSEC
      DATA             CTRN
     > /'+A','+B','+C','+D','+E','+F','+G','+H','+I','+J','+K','+L',
     >  '-A','-B','-C','-D','-E','-F','-G','-H','-I','-J','-K','-L'/
      DATA             IDSEC
     > /4,9,1,3,10,2/
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
      IDIR=1
      NSCELL(IDIR)=MAX(1,NOCELL(IDIR))
      NSUC(IDIR)=MAX(1,NUCELL(IDIR))
      NFCELL(IDIR)=0
      IDIR=3
      NSCELL(IDIR)=MAX(1,NOCELL(IDIR))
      NSUC(IDIR)=MAX(1,NUCELL(IDIR))
      IF(ISAXIS(IDIR) .EQ. -2) THEN
        NFCELL(IDIR)=NOCELL(IDIR)
      ELSE IF(ISAXIS(IDIR) .EQ. -1) THEN
        NFCELL(IDIR)=NOCELL(IDIR)-1
      ELSE
        NFCELL(IDIR)=0
      ENDIF
      IGEN=0
      IF(IHSYM .EQ. 9) THEN
*----
*  Process complete cell
*----
        DO IZ=0,NSCELL(3)-1
          IOFZ=(IZ+NFCELL(3))*NUCELL(1)
          DO IX=0,NSCELL(1)-1
            IGEN=IGEN+1
            IF(IGEN .GT. NBOCEL) CALL XABORT(NAMSBR//
     >        ': Cell number exceeds number of cells permitted')
            ILOCD=IOFZ+IX+NFCELL(1)+1
            IUNFLD(1,ILOCD)=IGEN
            IUNFLD(2,ILOCD)=1
          ENDDO
        ENDDO
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
          IOFZ=(IZ+NFCELL(3))*NUCELL(1)
          DO IX=0,NSCELL(1)-1
            ILOCD=IOFZ+IX+NFCELL(1)+1
            IGEN=IUNFLD(1,ILOCD)
            ITSYM(3,IGEN)=ISAXIS(3)
          ENDDO
        ENDIF
      ENDIF
      IF(ISAXIS(3) .NE. 0) THEN
        IF(ISAXIS(3) .EQ. -2) THEN
*----
*  SSYM Z-
*  Fill position IZR=1,NSCELL(3) with cells at
*  position IZD=NSUC(3)-IZR+1
*----
          DO IZ=1,NSCELL(3)
            IOFZR=(IZ-1)*NUCELL(1)
            IOFZ=(NSUC(3)-IZ)*NUCELL(1)
            DO IX=1,NSCELL(1)
              ILOCD=IOFZ+IX
              ILOCR=IOFZR+IX
              IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
              IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),-1)
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
            IOFZR=(IZ-1)*NUCELL(1)
            IOFZ=(NSUC(3)-IZ)*NUCELL(1)
            DO IX=1,NSCELL(1)
              ILOCD=IOFZ+IX
              ILOCR=IOFZR+IX
              IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
              IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),-1)
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
            IOFZ=(IZ-1)*NUCELL(1)
            IOFZR=(NSUC(3)-IZ)*NUCELL(1)
            DO IX=1,NSCELL(1)
              ILOCD=IOFZ+IX
              ILOCR=IOFZR+IX
              IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
              IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),-1)
            ENDDO
          ENDDO
        ELSE IF(ISAXIS(3) .EQ. 2) THEN
*----
*  SSYM Z+
*  Fill position IZR=NSUC(3)-IZD+1 with cell
*  at position IZD=1,NSCELL(3)
*----
          DO IZ=1,NSCELL(3)
            IOFZ=(IZ-1)*NUCELL(1)
            IOFZR=(NSUC(3)-IZ)*NUCELL(1)
            DO IX=1,NSCELL(1)
              ILOCD=IOFZ+IX
              ILOCR=IOFZR+IX
              IUNFLD(1,ILOCR)=IUNFLD(1,ILOCD)
              IUNFLD(2,ILOCR)=NXTTRS(IUNFLD(2,ILOCD),-1)
            ENDDO
          ENDDO
        ENDIF
        NSCELL(3)=MAX(1,NSUC(3))
      ENDIF
*----
*  1. Localize external faces
*  Find number of crown : NCR
*  Find number of cell inside last crown : NCILC (no external faces)
*  Only last crown has external faces
*  Face notation:
*  W=-1 /\ V=-2
* U=-1 |  | U=-2
*  V=-1 \/ W=-2
*  -/+ U -> IDFEX(1,*),IDFEX(2,*)
*  -/+ V -> IDFEX(3,*),IDFEX(4,*)
*  -/+ Z -> IDFEX(5,*),IDFEX(6,*)
*  -/+ W -> IDFEX(9,*),IDFEX(10,*)
*----
      ARGS=DBLE(12*NUCELL(1)-3)
      NCR=(NINT(SQRT(ARGS))+3)/6
      IF(NCR .EQ. 1) THEN
        NCILC=1
      ELSE
        NCILC=6*(NCR-1)
      ENDIF
      IF(NCILC .EQ. 1) THEN
        ILOCD=0
        DO IZ=1,NSUC(3)
          ILOCD=ILOCD+1
          IGEN=IUNFLD(1,ILOCD)
          IDFEX(1,IGEN)=1
          IDFEX(2,IGEN)=1
          IDFEX(3,IGEN)=1
          IDFEX(4,IGEN)=1
          IDFEX(9,IGEN)=1
          IDFEX(10,IGEN)=1
          IDFRT(1,IGEN)=IGEN
          IDFRT(2,IGEN)=IGEN
        ENDDO
      ELSE
        DO IZ=1,NSUC(3)
          ILOCD=NUCELL(1)*IZ-NCILC
*----
*  Scan over all sectors
*----
          ID1=10
          ID2=2
          DO ISECT=1,6
            ILOCD=ILOCD+1
            IGEN=IUNFLD(1,ILOCD)
            ID3=IDSEC(ISECT)
            IDFEX(ID1,IGEN)=1
            IDFEX(ID2,IGEN)=1
            IDFEX(ID3,IGEN)=1
            IDFRT(1,IGEN)=IGEN
            IDFRT(2,IGEN)=IGEN
            IDD=IDD+1
            DO ICR=1,NCR-2
              ILOCD=ILOCD+1
              IGEN=IUNFLD(1,ILOCD)
              IDFEX(ID2,IGEN)=1
              IDFEX(ID3,IGEN)=1
              IDFRT(1,IGEN)=IGEN
              IDFRT(2,IGEN)=IGEN
            ENDDO
            ID1=ID2
            ID2=ID3
          ENDDO
        ENDDO
      ENDIF
*----
*  2. Z Faces
*  IDFEX(5,*) for Z-
*  IDFEX(6,*) FOR Z+
*----
      DO IX=1,NUCELL(1)
        ILOCD=IX
        IF(IUNFLD(2,ILOCD) .EQ. 1) THEN
          IGEN=IUNFLD(1,ILOCD)
          IDFEX(5,IGEN)=1
          IDFRT(5,IGEN)=IGEN
        ENDIF
        ILOCD=NUCELL(1)*(NSUC(3)-1)+IX
        IF(IUNFLD(2,ILOCD) .EQ. 1) THEN
          IGEN=IUNFLD(1,ILOCD)
          IDFEX(6,IGEN)=1
          IDFRT(6,IGEN)=IGEN
        ENDIF
      ENDDO
*----
*  Process. Z translation
*----
      IF(ISAXIS(3) .EQ. 3) THEN
        DO IGEN=1,NBOCEL
          IF(IDFEX(5,IGEN) .EQ. 1) THEN
            DO IX=1,NUCELL(1)
              ILOCD=IX
              IF(IUNFLD(2,ILOCD) .EQ. 1) THEN
                IF(IUNFLD(1,ILOCD) .EQ. IGEN) THEN
                  ILOCD=NUCELL(1)*(NSUC(3)-1)+IX
                  IGENT=IUNFLD(1,ILOCD)
                  IDFRT(5,IGEN)=IGENT
                  IDFRT(6,IGENT)=IGEN
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
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
*  print routine output and
*  closing header if required
*  and return
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6002)
        IF(NDIM .EQ. 3) THEN
          DO IZ=1,NSUC(3)
            WRITE(IOUT,6003) IZ
            ILOCD=(IZ-1)*NUCELL(1)
            WRITE(IOUT,6005)
     >          (IUNFLD(1,IX),IX=ILOCD+1,ILOCD+NUCELL(1))
          ENDDO
        ELSE
          WRITE(IOUT,6004)
          WRITE(IOUT,6005)
     >    (IUNFLD(1,IX),IX=1,NUCELL(1))
        ENDIF
        WRITE(IOUT,6008)
        IF(NDIM .EQ. 3) THEN
          DO IZ=1,NSUC(3)
            WRITE(IOUT,6003) IZ
            ILOCD=(IZ-1)*NUCELL(1)
            WRITE(IOUT,6011)
     >        (CTRN(IUNFLD(2,IX)),IX=ILOCD+1,ILOCD+NUCELL(1))
          ENDDO
        ELSE
          WRITE(IOUT,6004)
          WRITE(IOUT,6011)
     >    (CTRN(IUNFLD(2,IX)),IX=1,NUCELL(1))
        ENDIF
        WRITE(IOUT,6006)
        DO ILOCD=1,NBOCEL
          WRITE(IOUT,6007)
     >    ILOCD,(ITSYM(IX,ILOCD),IX=1,4)
        ENDDO
        WRITE(IOUT,6009)
        DO ILOCD=1,NBOCEL
          WRITE(IOUT,6007) ILOCD,
     >    (IDFEX(IX,ILOCD),IX=1,4),(IDFEX(IX,ILOCD),IX=9,10),
     >    (IDFEX(IX,ILOCD),IX=5,6),IDFEX(0,ILOCD)
        ENDDO
        WRITE(IOUT,6010)
        DO ILOCD=1,NBOCEL
          WRITE(IOUT,6007) ILOCD,IDFRT(1,ILOCD),IDFRT(2,ILOCD),
     >    IDFRT(5,ILOCD),IDFRT(6,ILOCD)
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
 6003 FORMAT(' Hexagons for plane IZ =',I6)
 6004 FORMAT(' Hexagons ')
 6005 FORMAT(24(I7.7,1X))
 6006 FORMAT(/' Symmetrized cell     X    Y    Z    D')
 6007 FORMAT(' Cell ',I7.7,5X,20I5)
 6008 FORMAT(/' Cell rotations in assembly')
 6009 FORMAT(/' External faces   ',
     >        '   -U   +U   -V   +V   -W   +W   -Z   +Z   ND')
 6010 FORMAT(/' Coupled  faces      FH   LH    -Z   +Z  ')
 6011 FORMAT(24(5X,A2,1X))
      END
