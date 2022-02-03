*DECK SYBUP0
      SUBROUTINE SYBUP0(ZZR,ZZI,NSURF,NREG,SIGT,TRONC,A,B,IMPX,VOL,PIJ,
     1 PIS,PSS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the one-group collision, DP-0 leakage and DP-0 transmission
* probabilities in a Cartesian or hexagonal non-sectorized cell.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* ZZR     real tracking elements.
* ZZI     integer tracking elements.
* NSURF   number of surfaces (4 or 6).
* NREG    number of regions.
* SIGT    total macroscopic cross section.
* TRONC   voided block criterion.
* A       Cartesian dimension of the cell along the X axis or side of
*         the hexagon.
* B       Cartesian dimension of the cell along the Y axis.
* IMPX    print flag.
*
*Parameters: output
* VOL     volumes.
* PIJ     cellwise reduced collision probability matrices.
* PIS     cellwise reduced escape probability matrices.
* PSS     cellwise reduced transmission probability matrices.
*         PSS(i,j) is the probability from surface i to surface j.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   ZZI(*),NSURF,NREG,IMPX
      REAL      ZZR(*),SIGT(NREG),TRONC,A,B,VOL(NREG),PIJ(NREG,NREG),
     1          PIS(NSURF,NREG),PSS(NSURF,NSURF)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (SIGVID=1.0E-10)
      REAL      SURF(6)
      REAL, ALLOCATABLE, DIMENSION(:) :: G,PIJS
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LGFULL
*----
*  INLINE FUNCTIONS
*----
      INDPOS(I,J)=MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(G(NREG+6),LGFULL(NREG))
*----
*  CHECK FOR VOIDED REGIONS
*----
      DO 10 IR=1,NREG
      VOL(IR)=ZZR(IR)
      IF(VOL(IR).GT.0.) THEN
         DR=SQRT(VOL(IR))
      ELSE
         DR=0.0
      ENDIF
      LGFULL(IR)=(SIGT(IR)*DR).GT.TRONC
      IF(SIGT(IR).LE.SIGVID) SIGT(IR)=SIGVID
   10 CONTINUE
*----
*  COMPUTE SYMMETRIZED CP MATRIX
*----
      IOFI=ZZI(1)
      IOFR=ZZI(2)
      ICARE=ZZI(3)
      MNA=ZZI(4)
      ALLOCATE(PIJS((NREG+NSURF)*(NREG+NSURF+1)/2))
      CALL SYBUQV(ZZR(IOFR),ZZI(IOFI),NSURF,NREG,SIGT,MNA,LGFULL,PIJS)
*----
*  APPLY SYMMETRIES
*----
      IF(NSURF.EQ.4) THEN
         SURF(1)=0.25*B
         SURF(2)=0.25*B
         SURF(3)=0.25*A
         SURF(4)=0.25*A
      ELSE
         DO 20 JC=1,6
         SURF(JC)=0.25*A
   20    CONTINUE
      ENDIF
      DO 30 I=1,NSURF
      G(I)=SURF(I)
   30 CONTINUE
      IF(ICARE.EQ.1) THEN
*        RECTANGULAR CELL.
         PIJS(2)=2.0*PIJS(2)
         PIJS(5)=0.5*PIJS(5)
         PIJS(9)=2.0*PIJS(9)
         PIJS(4)=PIJS(5)
         PIJS(7)=PIJS(5)
         PIJS(8)=PIJS(5)
         IOF=11
         DO 50 I=1,NREG
         G(4+I)=SIGT(I)*VOL(I)
         SUM1=PIJS(IOF)+PIJS(IOF+1)
         SUM2=PIJS(IOF+2)+PIJS(IOF+3)
         PIJS(IOF)=SUM1
         PIJS(IOF+1)=SUM1
         PIJS(IOF+2)=SUM2
         PIJS(IOF+3)=SUM2
         DO 40 J=4,3+I
         PIJS(IOF+J)=2.0*PIJS(IOF+J)
   40    CONTINUE
         IOF=IOF+4+I
   50    CONTINUE
      ELSE IF(ICARE.EQ.2) THEN
*        SQUARE CELL.
         PIJS(9)=2.0*PIJS(9)
         PIJS(2)=PIJS(9)
         PIJS(4)=PIJS(5)
         PIJS(7)=PIJS(5)
         PIJS(8)=PIJS(5)
         IOF=11
         DO 80 I=1,NREG
         G(4+I)=SIGT(I)*VOL(I)
         SUM=PIJS(IOF)+PIJS(IOF+1)+PIJS(IOF+2)+PIJS(IOF+3)
         DO 60 J=0,3
         PIJS(IOF+J)=SUM
   60    CONTINUE
         DO 70 J=4,3+I
         PIJS(IOF+J)=4.0*PIJS(IOF+J)
   70    CONTINUE
         IOF=IOF+4+I
   80    CONTINUE
      ELSE IF(ICARE.EQ.3) THEN
*        HEXAGONAL CELL.
         PIJS(12)=2.0*PIJS(12)
         PIJS(7)=PIJS(12)
         PIJS(18)=PIJS(12)
         PIJS(2)=PIJS(20)
         PIJS(5)=PIJS(20)
         PIJS(9)=PIJS(20)
         PIJS(14)=PIJS(20)
         PIJS(16)=PIJS(20)
         PIJS(4)=PIJS(11)
         PIJS(8)=PIJS(11)
         PIJS(13)=PIJS(11)
         PIJS(17)=PIJS(11)
         PIJS(19)=PIJS(11)
         IOF=22
         DO 120 I=1,NREG
         G(6+I)=SIGT(I)*VOL(I)
         SUM=0.0
         DO 90 J=0,5
         SUM=SUM+PIJS(IOF+J)
   90    CONTINUE
         DO 100 J=0,5
         PIJS(IOF+J)=SUM
  100    CONTINUE
         DO 110 J=6,5+I
         PIJS(IOF+J)=6.0*PIJS(IOF+J)
  110    CONTINUE
         IOF=IOF+6+I
  120    CONTINUE
      ENDIF
*----
*  FIRST APPLY THE ORTHONORMALIZATION FACTOR
*----
      DO 130 I=1,(NSURF+NREG)*(NSURF+NREG+1)/2
      PIJS(I)=PIJS(I)*ZZR(IOFR)*ZZR(IOFR)
  130 CONTINUE
*----
*  VILLARINO-STAMM'LER NORMALIZATION
*----
      CALL SYBRHL(IMPX,NSURF,NREG,G,PIJS)
*----
*  REFORMAT PIJ, PIS AND PSS MATRICES
*----
      DO 150 I=1,NSURF
      DO 140 J=1,NSURF
      PSS(I,J)=PIJS(INDPOS(I,J))/SURF(I)
  140 CONTINUE
  150 CONTINUE
      DO 170 I=1,NSURF
      DO 160 J=1,NREG
      PIS(I,J)=PIJS(INDPOS(I,NSURF+J))/VOL(J)/SIGT(J)
  160 CONTINUE
  170 CONTINUE
      DO 190 I=1,NREG
      DO 180 J=1,NREG
      PIJ(I,J)=PIJS(INDPOS(NSURF+I,NSURF+J))/(VOL(I)*SIGT(I)*SIGT(J))
  180 CONTINUE
  190 CONTINUE
      DEALLOCATE(PIJS)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(LGFULL,G)
      RETURN
      END
