*DECK SYB7QG
      SUBROUTINE SYB7QG (IMPX,NCURR,MNA4,NRD,NSECT,LSECT,NREG,ZZR,ZZI,
     1 HSIDE,RAYRE,SIGTR,TRONC,VOL,PIJ,PVS,PSS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the one-group collision, leakage and transmission
* probabilities in a hexagonal sectorized cell.
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
* NCURR   type of interface current approximation:
*         =1 DP-0; =3 DP-1 interface currents.
* MNA4    number of angles in (0,$\\pi$/2).
* NRD     one plus the number of tubes in the cell.
* NSECT   number of sectors.
* LSECT   type of sectorization:
*         =-999  no sectorization / processed as a sectorized cell;
*         =-101  X-type sectorization of the coolant;
*         =-1    X-type sectorization of the cell.
* NREG    number of regions.
* ZZR     real tracking elements.
* ZZI     integer tracking elements.
* HSIDE   length of the hexagon sides.
* RAYRE   radius of the tubes.
* SIGTR   total macroscopic cross section.
* TRONC   voided block criterion.
*
*Parameters: output
* VOL     volumes.
* PIJ     volume to volume reduced probability.
* PVS     volume to surface probabilities.
* PSS     surface to surface probabilities in the following order:
*         PSS(i,j) is the probability from surface i to surface j.
*
*Comments:
*  hexagone surface identification.
*                     side 2
*                    xxxxxxxx
*                   x        x
*          side 3  x          x  side 1
*                 x            x
*                  x          x
*           side 4  x        x  side 6
*                    xxxxxxxx
*                     side 5
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,NCURR,MNA4,NRD,NSECT,LSECT,NREG,ZZI(*)
      REAL ZZR(*),HSIDE,RAYRE(NRD-1),SIGTR(NREG),TRONC,VOL(NREG),
     1 PIJ(NREG,NREG),PVS(6*NCURR,NREG),PSS(6*NCURR,6*NCURR)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (SIGVID=1.0E-10,NSURFQ=6)
      INTEGER IPER(3)
      REAL QSS(135)
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NUMREG
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLINT,WORKIJ,G
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
      ALLOCATE(VOLINT(NRD),WORKIJ(0:(NREG+6)*(NREG+7)/2-1),
     1 PSIX(0:5,NCURR,NREG),G(NREG+6))
      ALLOCATE(LGFULL(NREG))
*----
*  COMPUTE THE VOLUMES
*----
      CALL SYB7VO(NRD,HSIDE,RAYRE,VOLINT)
      IND=0
      DO 30 I=1,NRD-1
      IF(ABS(LSECT).GT.100) THEN
         IND=IND+1
         DO 10 ISEC=1,NSECT
         NUMREG(ISEC,I)=IND
   10    CONTINUE
      ELSE IF(LSECT.EQ.-1) THEN
         NUMREG(1,I)=IND+5
         NUMREG(2,I)=IND+6
         NUMREG(3,I)=IND+1
         NUMREG(4,I)=IND+2
         NUMREG(5,I)=IND+3
         NUMREG(6,I)=IND+4
         IND=IND+6
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
         NUMREG(1,I)=IND+5
         NUMREG(2,I)=IND+6
         NUMREG(3,I)=IND+1
         NUMREG(4,I)=IND+2
         NUMREG(5,I)=IND+3
         NUMREG(6,I)=IND+4
         IND=IND+6
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
      VOL(IND)=VOL(IND)+VOLINT(IR)/6.0
   70 CONTINUE
   75 CONTINUE
*----
*  CHECK FOR VOIDED REGIONS
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
      DO 90 IR=1,NSURFQ
      G(IR)=HSIDE/4.0
   90 CONTINUE
      DO 100 IR=1,NREG
      G(6+IR)=SIGTR(IR)*VOL(IR)
  100 CONTINUE
*     FIRST APPLY THE ORTHONORMALIZATION FACTOR:
      DO 110 I=0,(NSURFQ+NREG)*(NSURFQ+NREG+1)/2-1
      WORKIJ(I)=WORKIJ(I)*ZZR(MZRS)*ZZR(MZRS)
  110 CONTINUE
*
*     THEN PERFORM STAMM'LER NORMALIZATION:
      CALL SYBRHL(IMPX,NSURFQ,NREG,G,WORKIJ)
*
      IIJ=NSURFQ*(NSURFQ+1)/2-1
      DO 130 JR=1,NREG
      IIJ=IIJ+NSURFQ
      DO 120 IR=1,JR-1
      AUX=WORKIJ(IIJ+IR)/(SIGTR(IR)*SIGTR(JR))
      PIJ(IR,JR)=AUX/VOL(IR)
      PIJ(JR,IR)=AUX/VOL(JR)
  120 CONTINUE
      IIJ=IIJ+JR
      AUX=WORKIJ(IIJ)/(SIGTR(JR)*SIGTR(JR))
      PIJ(JR,JR)=AUX/VOL(JR)
  130 CONTINUE
*----
*  PIS AND PSS CALCULATION
*----
      IF(NCURR.GT.1) THEN
*        PERFORM A DP-1 CALCULATION USING THE TRACKING.
         CALL SYBUQ0(ZZR(MZRS),ZZI(MZIS),NSURFQ,NREG,SIGTR,MNA4,
     1   LGFULL,PSIX(0,1,1),QSS)
*
         DO 160 JS=0,NSURFQ-1
         DO 150 IH=1,NCURR
         DO 140 IR=1,NREG
         ZNOR=G(JS+1)+G(NSURFQ+IR)
         PSIX(JS,IH,IR)=ZNOR*PSIX(JS,IH,IR)/SIGTR(IR)/VOL(IR)
  140    CONTINUE
  150    CONTINUE
  160    CONTINUE
         IIQ=1
         DO 190 JS=0,NSURFQ-1
         DO 180 IS=0,JS-1
         ZNOR=G(IS+1)+G(JS+1)
         DO 170 IH=1,NCURR*NCURR
         QSS(IIQ)=ZNOR*QSS(IIQ)
         IIQ=IIQ+1
  170    CONTINUE
  180    CONTINUE
  190    CONTINUE
         DO 210 IIQ=1,135,9
         DO 200 IIS=1,8,2
         QSS(IIQ+IIS)=-QSS(IIQ+IIS)
  200    CONTINUE
  210    CONTINUE
      ELSE
*        RECOVER PSI AND PSS INFORMATION FROM DP-0 PIJ CALCULATION.
         IIQ=1
         IIJ=0
         DO 230 JS=0,NSURFQ-1
         DO 220 IS=0,JS-1
         QSS(IIQ)=4.0*WORKIJ(IIJ)
         IIQ=IIQ+NCURR*NCURR
         IIJ=IIJ+1
  220    CONTINUE
         IIJ=IIJ+1
  230    CONTINUE
         IIJ=NSURFQ*(NSURFQ+1)/2
         DO 250 IR=1,NREG
         DO 240 JS=0,NSURFQ-1
         PSIX(JS,1,IR)=WORKIJ(IIJ+JS)/SIGTR(IR)/VOL(IR)
  240    CONTINUE
         IIJ=IIJ+NSURFQ+IR
  250    CONTINUE
      ENDIF
*----
*  LOAD THE EURYDICE CP ARRAYS
*----
      DO 270 I=1,NREG
      DO 260 IH=1,NCURR
      SGN=1.0
      IF(IH.EQ.2) SGN=-1.0
      PVS(INC(1,IH),I)=SGN*PSIX(2,IH,I)
      PVS(INC(2,IH),I)=SGN*PSIX(3,IH,I)
      PVS(INC(3,IH),I)=SGN*PSIX(4,IH,I)
      PVS(INC(4,IH),I)=SGN*PSIX(5,IH,I)
      PVS(INC(5,IH),I)=SGN*PSIX(0,IH,I)
      PVS(INC(6,IH),I)=SGN*PSIX(1,IH,I)
  260 CONTINUE
  270 CONTINUE
      DO 290 I=1,6*NCURR
      DO 280 J=1,6*NCURR
      PSS(I,J)=0.0
  280 CONTINUE
  290 CONTINUE
      DO 310 IH=1,NCURR
      DO 300 JH=1,NCURR
      PSS(INC(2,IH),INC(1,JH))=QSS(INQ(IH,JH,6))/HSIDE
      PSS(INC(3,IH),INC(1,JH))=QSS(INQ(IH,JH,9))/HSIDE
      PSS(INC(4,IH),INC(1,JH))=QSS(INQ(JH,IH,13))/HSIDE
      PSS(INC(5,IH),INC(1,JH))=QSS(INQ(JH,IH,2))/HSIDE
      PSS(INC(6,IH),INC(1,JH))=QSS(INQ(JH,IH,3))/HSIDE
      PSS(INC(1,IH),INC(2,JH))=QSS(INQ(JH,IH,6))/HSIDE
      PSS(INC(3,IH),INC(2,JH))=QSS(INQ(IH,JH,10))/HSIDE
      PSS(INC(4,IH),INC(2,JH))=QSS(INQ(IH,JH,14))/HSIDE
      PSS(INC(5,IH),INC(2,JH))=QSS(INQ(JH,IH,4))/HSIDE
      PSS(INC(6,IH),INC(2,JH))=QSS(INQ(JH,IH,5))/HSIDE
      PSS(INC(1,IH),INC(3,JH))=QSS(INQ(JH,IH,9))/HSIDE
      PSS(INC(2,IH),INC(3,JH))=QSS(INQ(JH,IH,10))/HSIDE
      PSS(INC(4,IH),INC(3,JH))=QSS(INQ(IH,JH,15))/HSIDE
      PSS(INC(5,IH),INC(3,JH))=QSS(INQ(JH,IH,7))/HSIDE
      PSS(INC(6,IH),INC(3,JH))=QSS(INQ(JH,IH,8))/HSIDE
      PSS(INC(1,IH),INC(4,JH))=QSS(INQ(JH,IH,13))/HSIDE
      PSS(INC(2,IH),INC(4,JH))=QSS(INQ(JH,IH,14))/HSIDE
      PSS(INC(3,IH),INC(4,JH))=QSS(INQ(JH,IH,15))/HSIDE
      PSS(INC(5,IH),INC(4,JH))=QSS(INQ(JH,IH,11))/HSIDE
      PSS(INC(6,IH),INC(4,JH))=QSS(INQ(JH,IH,12))/HSIDE
      PSS(INC(1,IH),INC(5,JH))=QSS(INQ(IH,JH,2))/HSIDE
      PSS(INC(2,IH),INC(5,JH))=QSS(INQ(JH,IH,4))/HSIDE
      PSS(INC(3,IH),INC(5,JH))=QSS(INQ(IH,JH,7))/HSIDE
      PSS(INC(4,IH),INC(5,JH))=QSS(INQ(IH,JH,11))/HSIDE
      PSS(INC(6,IH),INC(5,JH))=QSS(INQ(IH,JH,1))/HSIDE
      PSS(INC(1,IH),INC(6,JH))=QSS(INQ(IH,JH,3))/HSIDE
      PSS(INC(2,IH),INC(6,JH))=QSS(INQ(IH,JH,5))/HSIDE
      PSS(INC(3,IH),INC(6,JH))=QSS(INQ(JH,IH,8))/HSIDE
      PSS(INC(4,IH),INC(6,JH))=QSS(INQ(IH,JH,12))/HSIDE
      PSS(INC(5,IH),INC(6,JH))=QSS(INQ(JH,IH,1))/HSIDE
  300 CONTINUE
  310 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(LGFULL)
      DEALLOCATE(G,PSIX,WORKIJ,VOLINT)
      DEALLOCATE(NUMREG)
      RETURN
      END
