*DECK SYB7TS
      SUBROUTINE SYB7TS(NA,NRD,NSECT,LSECT,NREG,HSIDE,RAYRE,ILIGN,IQW,
     1 DELR,LFAIRE,VOL,NZR,ZZR,NZI,ZZI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the tracking information related to an hexagonal sectorized
* heterogeneous cell.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NA      number of angles in (0,$\\pi$/2).
* NRD     one plus the number of tubes in the cell.
* NSECT   number of sectors.
* LSECT   type of sectorization:
*         =-999  no sectorization / processed as a sectorized cell;
*         =-101  X-type sectorization of the coolant;
*         =-1    X-type sectorization of the cell.
* NREG    number of regions in the cell.
* HSIDE   length of the hexagon side.
* RAYRE   radius of each cylinder.
* ILIGN   tracking print flag (=1 to print the tracking).
* IQW     equal weight quadrature flag (=1 to use equal weight
*         quadratures in angle and space).
* DELR    half distance between the tracks.
* LFAIRE  tracking calculation flag (=.FALSE. only compute the number
*         of tracks).
*
*Parameters: output
* VOL     volumes.
* NZR     number of real elements in vector ZZR.
* ZZR     real tracking information.
* NZI     number of integer elements in vector ZZI.
* ZZI     integer tracking information.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NA,NRD,NSECT,LSECT,NREG,ILIGN,IQW,NZR,NZI,ZZI(*)
      REAL HSIDE,RAYRE(NRD-1),DELR,VOL(NREG),ZZR(*)
      LOGICAL LFAIRE
*----
*  LOCAL VARIABLES
*----
      PARAMETER(DXMIN=1.E-3,PIO2=1.570796327,SQRT32=1.7320508075689/2.0)
      REAL ZA(64),WA(64)
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NUMREG
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLINT,XCOTE
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NUMREG(NSECT,NRD))
      ALLOCATE(VOLINT(NRD),XCOTE(NRD))
*
      IF(NA.GT.64) CALL XABORT('SYB7TS: NA IS GREATER THAN 64.')
      IF(RAYRE(NRD-1).GT.HSIDE) CALL XABORT('SYB7TS: A RADIUS IS GREAT'
     1  //'ER THAN THE HEXAGON SIDE LENGTH.')
      IF(IQW.EQ.0) THEN
*        GAUSS-LEGENDRE INTEGRATION POINTS.
         CALL ALGPT(NA,-1.0,1.0,ZA,WA)
      ELSE
*        EQUAL WEIGHT INTEGRATION POINTS.
         DO 10 I=1,NA
         ZA(I)=(2.0*REAL(I)-1.0)/REAL(NA)-1.0
         WA(I)=2.0/REAL(NA)
   10    CONTINUE
      ENDIF
*----
*  COMPUTE THE VOLUMES AND NUMREG
*----
      CALL SYB7VO(NRD,HSIDE,RAYRE,VOLINT)
      IND=0
      DO 50 I=1,NRD-1
      IF(ABS(LSECT).GT.100) THEN
         IND=IND+1
         DO 30 ISEC=1,NSECT
         NUMREG(ISEC,I)=IND
   30    CONTINUE
      ELSE IF(LSECT.EQ.-1) THEN
         NUMREG(1,I)=IND+5
         NUMREG(2,I)=IND+6
         NUMREG(3,I)=IND+1
         NUMREG(4,I)=IND+2
         NUMREG(5,I)=IND+3
         NUMREG(6,I)=IND+4
         IND=IND+6
      ELSE
         DO 40 ISEC=1,NSECT
         IND=IND+1
         NUMREG(ISEC,I)=IND
   40    CONTINUE
      ENDIF
   50 CONTINUE
      IF(LSECT.EQ.-999) THEN
         IND=IND+1
         DO 60 ISEC=1,NSECT
         NUMREG(ISEC,I)=IND
   60    CONTINUE
      ELSE IF((LSECT.EQ.-1).OR.(LSECT.EQ.-101)) THEN
         NUMREG(1,I)=IND+5
         NUMREG(2,I)=IND+6
         NUMREG(3,I)=IND+1
         NUMREG(4,I)=IND+2
         NUMREG(5,I)=IND+3
         NUMREG(6,I)=IND+4
         IND=IND+6
      ELSE
         DO 70 ISEC=1,NSECT
         IND=IND+1
         NUMREG(ISEC,I)=IND
   70    CONTINUE
      ENDIF
      DO 80 I=1,NREG
      VOL(I)=0.0
   80 CONTINUE
      DO 95 IR=1,NRD
      DO 90 IS=1,NSECT
      IND=NUMREG(IS,IR)
      VOL(IND)=VOL(IND)+VOLINT(IR)/6.0
   90 CONTINUE
   95 CONTINUE
*----
*  INTERSECTION OF THE HEXAGON SIDE WITH THE TUBES
*----
      HAUTEU=HSIDE*SQRT32
      H2=HAUTEU*HAUTEU
      DO 100 MRE=NRD-1,1,-1
      XI=RAYRE(MRE)*RAYRE(MRE)-H2
      IF(XI.GT.0.0) THEN
         XCOTE(MRE)=SQRT(XI)
      ELSE
         JMINR=MRE+1
         GO TO 110
      ENDIF
  100 CONTINUE
      JMINR=1
*
  110 NXMIN=999999999
      NXMAX=0
      CALL SYB7T0(NA,NRD,HSIDE,RAYRE,JMINR,XCOTE,LFAIRE,DELR,IQW,
     1 WA,ZA,NXMIN,NXMAX,MZRS,ZZR(1),MZIS,ZZI(3))
*
      IF(LFAIRE) THEN
*        SET ZZI(1:2) AND COMPUTE THE NUMERICAL ORTHONORMALIZATION
*        FACTORS.
         ZZI(1)=MZIS+3
         ZZI(2)=MZRS+1
         ZN1=0.0
         ZN2=0.0
         ZN3=0.0
         DO 120 IA=1,NA
         PHI=0.5*PIO2*(ZA(IA)+1.0)
         SI=SIN(PHI)
         ZN1=ZN1+SI*WA(IA)
         ZN2=ZN2+SI*SI*WA(IA)
         ZN3=ZN3+SI*SI*SI*WA(IA)
  120    CONTINUE
         ZN1=0.5*ZN1*PIO2
         ZN2=0.5*ZN2*PIO2
         ZN3=0.5*ZN3*PIO2
         ZZR(MZRS+1)=1.0/SQRT(ZN1)
         ZZR(MZRS+2)=1.0/SQRT(0.75*ZN3-0.7205061948*ZN2*ZN2/ZN1)
         ZZR(MZRS+3)=ZZR(MZRS+2)*0.8488263632*ZN2/ZN1
         ZZR(MZRS+4)=2.0/SQRT(3.0*(ZN1-ZN3))
         IF(ILIGN.GT.0) WRITE (6,210) (ZZR(MZRS+I),I=1,4)
*
*        UNFOLD THE TRACKS.
         IZI=MZIS+2
         IZR=MZRS+4
         DO 140 ISYM=-1,1,2
         DO 130 IFAC=3,8
         MZIR=MZIS
         MZRR=MZRS
         CALL SYB7TR(NA,NRD,MZIS,MZRS,IFAC,ISYM,NUMREG,ZZI(3),ZZR(1),
     1   MZIR,MZRR,ZZI(IZI+1),ZZR(IZR+1))
         IZI=IZI+MZIR
         IZR=IZR+MZRR
  130    CONTINUE
  140    CONTINUE
         NZI=IZI
         NZR=IZR
      ELSE
         NZI=13*MZIS+3
         NZR=13*MZRS+5
      ENDIF
*
      IF((ILIGN.GT.0).AND.(.NOT.LFAIRE)) THEN
         WRITE(6,200) NA,NRD,NSECT,HSIDE,DXMIN,DELR,NZI,NZR,NXMIN,
     1   NXMAX
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XCOTE,VOLINT)
      DEALLOCATE(NUMREG)
      RETURN
*
  200 FORMAT(/49H SYB7TS: TRACKING OF A SECTORIZED HEXAGONAL CELL./
     1 7H NA    ,I8,29H   (NUMBER OF ANGLES IN PI/2)/
     2 7H NRD   ,I8,22H   (1+NUMBER OF TUBES)/
     3 7H NSECT ,I8,22H   (NUMBER OF SECTORS)/
     4 7H HSIDE ,1P,E8.1,17H   (HEXAGON SIDE)/
     5 7H DXMIN ,1P,E8.1,24H   (GEOMETRICAL EPSILON)/
     6 7H DELR  ,1P,E8.1,37H   (HALF DISTANCE BETWEEN THE TRACKS)/
     7 7H NZI   ,I8,40H   (NUMBER OF INTEGER TRACKING ELEMENTS)/
     8 7H NZR   ,I8,37H   (NUMBER OF REAL TRACKING ELEMENTS)/
     9 7H NXMIN ,I8,37H   (MINIMUM NB. OF TRACKS PER REGION)/
     1 7H NXMAX ,I8,37H   (MAXIMUM NB. OF TRACKS PER REGION))
  210 FORMAT (/47H SYB7TS: NUMERICAL ORTHONORMALIZATION FACTORS =,1P,
     1 4E12.4/)
      END
