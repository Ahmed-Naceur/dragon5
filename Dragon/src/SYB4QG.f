*DECK SYB4QG
      SUBROUTINE SYB4QG (IMPX,NCURR,MNA4,NRD,NSECT,LSECT,NREG,ZZR,ZZI,
     1 A,B,RAYRE,SIGTR,TRONC,VOL,PIJ,PVS,PSS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the one-group collision, leakage and transmission 
* probabilities in a Cartesian sectorized cell.
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
* IMPX    print parameter (equal to zero for no print).
* NCURR   type of interface currents (=1 DP-0; =3 DP-1).
* MNA4    number of angles in (0,$\\pi$/2).
* NRD     one plus the number of tubes in the cell.
* NSECT   number of sectors.
* LSECT   type of sectorization:
*         =-999 no sectorization / processed as a sectorized cell;
*         =-101 X-type sectorization of the coolant;
*         =-1   X-type sectorization of the cell;
*         =101  +-type sectorization of the coolant;
*         =1    +-type sectorization of the cell;
*         =102  + and X-type sectorization of the coolant;
*         =2    + and X-type sectorization of the cell.
* NREG    number of regions.
* ZZR     real tracking elements.
* ZZI     integer tracking elements.
* A       size of the external X side.
* B       size of the external Y side.
* RAYRE   radius of the tubes.
* SIGTR   total macroscopic cross section.
* TRONC   voided block criterion.
*
*Parameters: output
* VOL     volumes.
* PIJ     volume to volume reduced probability.
* PVS     volume to surface probabilities:
*         XINF surface 1;  XSUP surface 2;
*         YINF surface 3;  YSUP surface 4.
* PSS     surface to surface probabilities in the following order:
*         PSS(i,j) is the probability from surface i to surface j.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,NCURR,MNA4,NRD,NSECT,LSECT,NREG,ZZI(*)
      REAL ZZR(*),A,B,RAYRE(NRD-1),SIGTR(NREG),TRONC,VOL(NREG),
     1 PIJ(NREG,NREG),PVS(4*NCURR,NREG),PSS(4*NCURR,4*NCURR)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (SIGVID=1.0E-10,NSURFQ=4)
      INTEGER IPER(3)
      REAL QSS(54)
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NUMREG
      REAL, ALLOCATABLE, DIMENSION(:) :: WORKIJ,G
      REAL, ALLOCATABLE, DIMENSION(:,:) :: VOLINT
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: PSIX
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LGFULL
*----
*  DATA STATEMENT AND INLINE FUNCTIONS
*----
      SAVE IPER
      DATA IPER/1,3,2/
      INC(IC,IH)=(IC-1)*NCURR+IPER(IH)
      INQ(IH,JH,IS)=(IS-1)*NCURR*NCURR+(IH-1)*NCURR+JH
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NUMREG(NSECT,NRD))
      ALLOCATE(VOLINT(NSECT,NRD),WORKIJ(0:(NREG+4)*(NREG+5)/2-1),
     1 PSIX(0:3,NCURR,NREG),G(NREG+4))
      ALLOCATE(LGFULL(NREG))
*----
*  COMPUTE THE VOLUMES
*----
      CALL SYB4VO(NSECT,NRD,A,B,RAYRE,VOLINT)
      IND=0
      DO 30 I=1,NRD-1
      IF(ABS(LSECT).GT.100) THEN
         IND=IND+1
         DO 10 ISEC=1,NSECT
         NUMREG(ISEC,I)=IND
   10    CONTINUE
      ELSE IF(LSECT.EQ.-1) THEN
         NUMREG(1,I)=IND+4
         NUMREG(2,I)=IND+1
         NUMREG(3,I)=IND+1
         NUMREG(4,I)=IND+2
         NUMREG(5,I)=IND+2
         NUMREG(6,I)=IND+3
         NUMREG(7,I)=IND+3
         NUMREG(8,I)=IND+4
         IND=IND+4
      ELSE
         DO 20 ISEC=1,NSECT
         IND=IND+1
         NUMREG(ISEC,I)=IND
   20    CONTINUE
      ENDIF
   30 CONTINUE
      IF(LSECT.EQ.-999) THEN
         IND=IND+1
         DO 40 ISEC=1,NSECT
         NUMREG(ISEC,I)=IND
   40    CONTINUE
      ELSE IF((LSECT.EQ.-1).OR.(LSECT.EQ.-101)) THEN
         NUMREG(1,I)=IND+4
         NUMREG(2,I)=IND+1
         NUMREG(3,I)=IND+1
         NUMREG(4,I)=IND+2
         NUMREG(5,I)=IND+2
         NUMREG(6,I)=IND+3
         NUMREG(7,I)=IND+3
         NUMREG(8,I)=IND+4
         IND=IND+4
      ELSE
         DO 50 ISEC=1,NSECT
         IND=IND+1
         NUMREG(ISEC,I)=IND
   50    CONTINUE
      ENDIF
      DO 60 I=1,NREG
      VOL(I)=0.0
   60 CONTINUE
      DO 75 IR=1,NRD
      DO 70 IS=1,NSECT
      IND=NUMREG(IS,IR)
      VOL(IND)=VOL(IND)+VOLINT(IS,IR)
   70 CONTINUE
   75 CONTINUE
*----
*  CHECH FOR VOIDED REGIONS
*----
      DO 80 IR=1,NREG
      IF(VOL(IR) .GT. 0.) THEN
         DR=SQRT(VOL(IR))
      ELSE
         DR=0.0
      ENDIF
      LGFULL(IR)=(SIGTR(IR)*DR).GT.TRONC
      IF(SIGTR(IR).LE.SIGVID) SIGTR(IR)=SIGVID
   80 CONTINUE
*----
*  COMPUTE COLLISION, DP-0 ESCAPE AND DP-0 TRANSMISSION PROBABILITIES
*----
      MZIS=ZZI(1)
      MZRS=ZZI(2)
      CALL SYBUQV(ZZR(MZRS),ZZI(MZIS),NSURFQ,NREG,SIGTR,MNA4,LGFULL,
     1 WORKIJ)
*----
*  STAMM'LER RENORMALIZATION
*----
      G(1)=A/4.0
      G(2)=B/4.0
      G(3)=A/4.0
      G(4)=B/4.0
      DO 100 IR=1,NREG
      G(4+IR)=SIGTR(IR)*VOL(IR)
  100 CONTINUE
*     FIRST APPLY THE ORTHONORMALIZATION FACTOR:
      DO 105 I=0,(NSURFQ+NREG)*(NSURFQ+NREG+1)/2-1
      WORKIJ(I)=WORKIJ(I)*ZZR(MZRS)*ZZR(MZRS)
  105 CONTINUE
*
*     THEN PERFORM STAMM'LER NORMALIZATION:
      CALL SYBRHL(IMPX,NSURFQ,NREG,G,WORKIJ)
*
      IIJ=NSURFQ*(NSURFQ+1)/2-1
      DO 120 JR=1,NREG
      IIJ=IIJ+NSURFQ
      DO 110 IR=1,JR-1
      AUX=WORKIJ(IIJ+IR)/(SIGTR(IR)*SIGTR(JR))
      PIJ(IR,JR)=AUX/VOL(IR)
      PIJ(JR,IR)=AUX/VOL(JR)
  110 CONTINUE
      IIJ=IIJ+JR
      AUX=WORKIJ(IIJ)/(SIGTR(JR)*SIGTR(JR))
      PIJ(JR,JR)=AUX/VOL(JR)
  120 CONTINUE
*----
*  PIS AND PSS CALCULATION
*----
      IF(NCURR.GT.1) THEN
*        PERFORM A DP-1 CALCULATION USING THE TRACKING.
         CALL SYBUQ0(ZZR(MZRS),ZZI(MZIS),NSURFQ,NREG,SIGTR,MNA4,
     1   LGFULL,PSIX(0,1,1),QSS)
*
         DO 132 JS=0,NSURFQ-1
         DO 131 IH=1,NCURR
         DO 130 IR=1,NREG
         ZNOR=G(JS+1)+G(NSURFQ+IR)
         PSIX(JS,IH,IR)=ZNOR*PSIX(JS,IH,IR)/SIGTR(IR)/VOL(IR)
  130    CONTINUE
  131    CONTINUE
  132    CONTINUE
         IIQ=1
         DO 142 JS=0,NSURFQ-1
         DO 141 IS=0,JS-1
         ZNOR=G(IS+1)+G(JS+1)
         DO 140 IH=1,NCURR*NCURR
         QSS(IIQ)=ZNOR*QSS(IIQ)
         IIQ=IIQ+1
  140    CONTINUE
  141    CONTINUE
  142    CONTINUE
      ELSE
*        RECOVER PSI AND PSS INFORMATION FROM DP-0 PIJ CALCULATION.
         IIQ=1
         IIJ=0
         DO 160 JS=0,NSURFQ-1
         DO 150 IS=0,JS-1
         QSS(IIQ)=4.0*WORKIJ(IIJ)
         IIQ=IIQ+NCURR*NCURR
         IIJ=IIJ+1
  150    CONTINUE
         IIJ=IIJ+1
  160    CONTINUE
         IIJ=NSURFQ*(NSURFQ+1)/2
         DO 180 IR=1,NREG
         DO 170 JS=0,NSURFQ-1
         PSIX(JS,1,IR)=WORKIJ(IIJ+JS)/SIGTR(IR)/VOL(IR)
  170    CONTINUE
         IIJ=IIJ+NSURFQ+IR
  180    CONTINUE
      ENDIF
*----
*  LOAD THE EURYDICE CP ARRAYS
*----
      DO 191 I=1,NREG
      DO 190 IH=1,NCURR
      PVS(INC(1,IH),I)=PSIX(3,IH,I)
      PVS(INC(2,IH),I)=PSIX(1,IH,I)
      PVS(INC(3,IH),I)=PSIX(0,IH,I)
      PVS(INC(4,IH),I)=PSIX(2,IH,I)
  190 CONTINUE
  191 CONTINUE
      DO 201 I=1,4*NCURR
      DO 200 J=1,4*NCURR
      PSS(I,J)=0.0
  200 CONTINUE
  201 CONTINUE
      DO 220 IH=1,NCURR
      DO 210 JH=1,NCURR
      PSS(INC(2,IH),INC(1,JH))=QSS(INQ(IH,JH,5))/B
      PSS(INC(3,IH),INC(1,JH))=QSS(INQ(JH,IH,4))/A
      PSS(INC(4,IH),INC(1,JH))=QSS(INQ(JH,IH,6))/A
      PSS(INC(1,IH),INC(2,JH))=QSS(INQ(IH,JH,5))/B
      PSS(INC(3,IH),INC(2,JH))=QSS(INQ(JH,IH,1))/A
      PSS(INC(4,IH),INC(2,JH))=QSS(INQ(IH,JH,3))/A
      PSS(INC(1,IH),INC(3,JH))=QSS(INQ(IH,JH,4))/B
      PSS(INC(2,IH),INC(3,JH))=QSS(INQ(IH,JH,1))/B
      PSS(INC(4,IH),INC(3,JH))=QSS(INQ(IH,JH,2))/A
      PSS(INC(1,IH),INC(4,JH))=QSS(INQ(IH,JH,6))/B
      PSS(INC(2,IH),INC(4,JH))=QSS(INQ(JH,IH,3))/B
      PSS(INC(3,IH),INC(4,JH))=QSS(INQ(IH,JH,2))/A
  210 CONTINUE
  220 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(LGFULL)
      DEALLOCATE(G,PSIX,WORKIJ,VOLINT)
      DEALLOCATE(NUMREG)
      RETURN
      END
