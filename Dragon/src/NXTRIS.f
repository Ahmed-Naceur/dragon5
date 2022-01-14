*DECK NXTRIS
      SUBROUTINE NXTRIS(IPRINT,ITYPG ,MAXMSH,NREG  ,ITRN  ,ITST  ,
     >                  ITSYM ,NM    ,MIX   ,ISPLT ,DAMESH,
     >                  NMS   ,MIXR  ,ISPLTR,DAMESR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Rotate geometry according to reference turn and test, if required,
* in such a way that it satisfies intrinsic symmetries.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPRINT  intermediate printing level for output.
* ITYPG   geometry type.
* MAXMSH  maximum number of elements in MESH array.
* NREG    number of elements in MIX array.
* ITRN    geometry original turn number.
* ITST    flag for testing symmetry.
* ITSYM   flag for symmetries to test.
*
*Parameters: input/output
* NM      mesh size in all directions ($X$, $Y$, $Z$ and $R$).
* MIX     final mixture description for geometry (including HMIX).
* ISPLT   final split desctiption for geometry.
* DAMESH  final mesh description for geometry.
* NMS     mesh size after splitting.
*
*Parameters: temporary storage
* MIXR    mixture description for rotated geometry (including HMIX).
* ISPLTR  split desctiption for rotated geometry.
* DAMESR  mesh description for rotated geometry.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,ITYPG,MAXMSH,NREG,ITRN,ITST,
     >                 NM(4),ITSYM(4)
      INTEGER          ISPLT(0:MAXMSH-1,4),MIX(0:NREG-1,2)
      DOUBLE PRECISION DAMESH(-1:MAXMSH,4)
      INTEGER          NMS(4),ISPLTR(0:MAXMSH-1,4,2),
     >                 MIXR(0:NREG-1,2,2)
      DOUBLE PRECISION DAMESR(-1:MAXMSH,4,2)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTRIS')
      DOUBLE PRECISION DCUTOF,DZERO,DONE,DTWO
      PARAMETER       (DCUTOF=1.0D-6,DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Functions
*----
      INTEGER          NXTTRS
*----
*  Local variables
*----
      INTEGER          NR,NX,NY,NZ,ITM(4,2),NPG,IPG,IG,ICT,ITG,
     >                 IDIR,IKT,IDMI,ITMI,IX,IY,IZ,IR,NRP1,NMR,
     >                 NMT(4),NMTS(4),NMTMP
      DOUBLE PRECISION DDD
*----
*  Data
*----
      CHARACTER        CDIR(1:4)*1
      SAVE             CDIR
      DATA             CDIR /'X','Y','Z','R'/
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
*----
*  Turn reference geometry (IPG=1)
*  and symmetric geometries (IPG=2,3,4,5)
*----
      ICT=0
      NX=NM(1)
      NY=NM(2)
      NZ=MAX(NM(3),1)
      NR=NM(4)
      NRP1=NR+1
      NMR=NR
      IF(ITYPG .EQ.  3 .OR. ITYPG .EQ.  6 .OR.
     >   ITYPG .EQ. 10 .OR. ITYPG .EQ. 11 ) THEN
        NRP1=NR
        NMR=NR-1
      ENDIF
      ITM(3,1)=3
      ITM(3,2)=3
      ITM(4,1)=4
      ITM(4,2)=4
      NPG=1
      IF(ITST .EQ. 1) NPG=5
      DO IPG=1,NPG
        IF(IPG .EQ. 1) THEN
          IG=1
          ICT=ITRN
          DO IX=0,NR-1
            DAMESR(IX,4,IG)=DAMESH(IX,4)
            ISPLTR(IX,4,IG)=ISPLT(IX,4)
          ENDDO
          DAMESR(NR,4,IG)=DAMESH(NR,4)
        ELSE
          IG=2
          ITG=IPG-1
          IF(ABS(ITSYM(ITG)) .GE. 1) THEN
*----
*  Symmetry is valid
*  Determine final turn after applying symmetry on
*  current turn
*----
            IF(ITG .EQ. 1) THEN
*----
*  Symmetry in X
*----
              ICT=NXTTRS(ITRN,1)
            ELSE IF(ITG .EQ. 2) THEN
*----
*  Symmetry in Y
*----
              ICT=NXTTRS(ITRN,3)
            ELSE IF(ITG .EQ. 3) THEN
*----
*  Symmetry in Z
*----
              ICT=NXTTRS(ITRN,-1)
            ELSE IF(ITG .EQ. 4) THEN
*----
*  Symmetry in X=Y or X=-Y
*----
              IF(ABS(ITSYM(ITG)) .EQ. 1) THEN
                ICT=NXTTRS(ITRN,2)
              ELSE
                ICT=NXTTRS(ITRN,4)
              ENDIF
            ENDIF
          ELSE
*----
*  No need to test the geometry for this
*  intrinsic symmetry.
*----
            GO TO 1005
          ENDIF
        ENDIF
        IF(ICT .GT. 12 ) THEN
          IKT=12-ICT
        ELSE
          IKT=ICT
        ENDIF
        DO IX=0,NR-1
          DAMESR(IX,4,IG)=DAMESH(IX,4)
          ISPLTR(IX,4,IG)=ISPLT(IX,4)
        ENDDO
        DAMESR(NR,4,IG)=DAMESH(NR,4)
        IF(IKT .LT. 0) THEN
          DAMESR(-1,3,IG)=-DAMESH(-1,3)
          DAMESR(-1,4,IG)=-DAMESH(-1,4)
        ELSE
          DAMESR(-1,3,IG)=DAMESH(-1,3)
          DAMESR(-1,4,IG)=DAMESH(-1,4) 
        ENDIF
        IF (ABS(IKT) .EQ. 1) THEN
          ITM(1,IG)=1
          ITM(2,IG)=2
          DO 100 IX=0,NX-1
            DAMESR(IX,1,IG)=DAMESH(IX+1,1)-DAMESH(IX,1)
            ISPLTR(IX,1,IG)=ISPLT(IX,1)
 100      CONTINUE
          DAMESR(-1,1,IG)=DAMESH(-1,1)
          DO 110 IY=0,NY-1
            DAMESR(IY,2,IG)=DAMESH(IY+1,2)-DAMESH(IY,2)
            ISPLTR(IY,2,IG)=ISPLT(IY,2)
 110      CONTINUE
          DAMESR(-1,2,IG)=DAMESH(-1,2)
          IF(IKT .LT. 0) THEN
            DO 120 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(NZ-IZ,3)-DAMESH(NZ-IZ-1,3)
              ISPLTR(IZ,3,IG)=ISPLT(NZ-IZ-1,3)
              ITMI=IZ*NX*NY*NRP1
              IDMI=(NZ-IZ-1)*NX*NY*NRP1
              DO 121 IY=0,NY-1
                DO 122 IX=0,NX-1
                  DO 123 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 123              CONTINUE
 122            CONTINUE
 121          CONTINUE
 120        CONTINUE
          ELSE
            DO 130 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(IZ+1,3)-DAMESH(IZ,3)
              ISPLTR(IZ,3,IG)=ISPLT(IZ,3)
              ITMI=IZ*NX*NY*NRP1
              IDMI=ITMI
              DO 131 IY=0,NY-1
                DO 132 IX=0,NX-1
                  DO 133 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 133              CONTINUE
 132            CONTINUE
 131          CONTINUE
 130        CONTINUE
          ENDIF
        ELSE IF(ABS(IKT) .EQ. 2) THEN
*----
*  ROTATION OF PI/2
*----
          ITM(1,IG)=2
          ITM(2,IG)=1
          DO 200 IX=0,NY-1
            DAMESR(IX,1,IG)=DAMESH(IX+1,2)-DAMESH(IX,2)
            ISPLTR(IX,1,IG)=ISPLT(IX,2)
 200      CONTINUE
          DAMESR(-1,1,IG)=DAMESH(-1,2)
          DO 210 IY=0,NX-1
            DAMESR(IY,2,IG)=DAMESH(NX-IY,1)-DAMESH(NX-IY-1,1)
            ISPLTR(IY,2,IG)=ISPLT(NX-IY-1,1)
 210      CONTINUE
          DAMESR(-1,2,IG)=-DAMESH(-1,1)
          IF(IKT .LT. 0) THEN
            DO 220 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(NZ-IZ,3)-DAMESH(NZ-IZ-1,3)
              ISPLTR(IZ,3,IG)=ISPLT(NZ-IZ-1,3)
              DO 221 IY=0,NX-1
                DO 222 IX=0,NY-1
                  ITMI=IZ*NX*NY*NRP1+IY*NY*NRP1+
     >                 IX*NRP1
                  IDMI=(NZ-IZ-1)*NX*NY*NRP1+(NY-IX-1)*NX*NRP1+
     >                 IY*NRP1
                  DO 223 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 223              CONTINUE
 222            CONTINUE
 221          CONTINUE
 220        CONTINUE
          ELSE
            DO 230 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(IZ+1,3)-DAMESH(IZ,3)
              ISPLTR(IZ,3,IG)=ISPLT(IZ,3)
              DO 231 IY=0,NX-1
                DO 232 IX=0,NY-1
                  ITMI=IZ*NX*NY*NRP1+IY*NY*NRP1+
     >                 IX*NRP1
                  IDMI=IZ*NX*NY*NRP1+IX*NX*NRP1+
     >                 (NX-IY-1)*NRP1
                  DO 233 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 233              CONTINUE
 232            CONTINUE
 231          CONTINUE
 230        CONTINUE
          ENDIF
        ELSE IF(ABS(IKT) .EQ. 3) THEN
*----
*  ROTATION OF PI
*----
          ITM(1,IG)=1
          ITM(2,IG)=2
          DO 300 IX=0,NX-1
            DAMESR(IX,1,IG)=DAMESH(NX-IX,1)-DAMESH(NX-IX-1,1)
            ISPLTR(IX,1,IG)=ISPLT(NX-IX-1,1)
 300      CONTINUE
          DAMESR(-1,1,IG)=-DAMESH(-1,1)
          DO 310 IY=0,NY-1
            DAMESR(IY,2,IG)=DAMESH(NY-IY,2)-DAMESH(NY-IY-1,2)
            ISPLTR(IY,2,IG)=ISPLT(NY-IY-1,2)
 310      CONTINUE
          DAMESR(-1,2,IG)=-DAMESH(-1,2)
          IF(IKT .LT. 0) THEN
            DO 320 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(NZ-IZ,3)-DAMESH(NZ-IZ-1,3)
              ISPLTR(IZ,3,IG)=ISPLT(NZ-IZ-1,3)
              DO 321 IY=0,NY-1
                DO 322 IX=0,NX-1
                  ITMI=IZ*NX*NY*NRP1+IY*NX*NRP1+
     >                 IX*NRP1
                  IDMI=(NZ-IZ-1)*NX*NY*NRP1+(NY-IY-1)*NX*NRP1+
     >                 (NX-IX-1)*NRP1
                  DO 323 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 323              CONTINUE
 322            CONTINUE
 321          CONTINUE
 320        CONTINUE
          ELSE
            DO 330 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(IZ+1,3)-DAMESH(IZ,3)
              ISPLTR(IZ,3,IG)=ISPLT(IZ,3)
              DO 331 IY=0,NY-1
                DO 332 IX=0,NX-1
                  ITMI=IZ*NX*NY*NRP1+IY*NX*NRP1+
     >                 IX*NRP1
                  IDMI=IZ*NX*NY*NRP1+(NY-IY-1)*NX*NRP1+
     >                 (NX-IX-1)*NRP1
                  DO 333 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 333              CONTINUE
 332            CONTINUE
 331          CONTINUE
 330        CONTINUE
          ENDIF
        ELSE IF(ABS(IKT) .EQ. 4) THEN
*----
*  ROTATION OF 3*PI/2
*----
          ITM(1,IG)=2
          ITM(2,IG)=1
          DO 400 IX=0,NY-1
            DAMESR(IX,1,IG)=DAMESH(NY-IX,2)-DAMESH(NY-IX-1,2)
            ISPLTR(IX,1,IG)=ISPLT(NY-IX-1,2)
 400      CONTINUE
          DAMESR(-1,1,IG)=-DAMESH(-1,2)
          DO 410 IY=0,NX-1
            DAMESR(IY,2,IG)=DAMESH(IY+1,1)-DAMESH(IY,1)
            ISPLTR(IY,2,IG)=ISPLT(IY,1)
 410      CONTINUE
          DAMESR(-1,2,IG)=DAMESH(-1,1)
          IF(IKT .LT. 0) THEN
            DO 420 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(NZ-IZ,3)-DAMESH(NZ-IZ-1,3)
              ISPLTR(IZ,3,IG)=ISPLT(NZ-IZ-1,3)
              DO 421 IY=0,NX-1
                DO 422 IX=0,NY-1
                  ITMI=IZ*NX*NY*NRP1+IY*NY*NRP1+
     >                 IX*NRP1
                  IDMI=(NZ-IZ-1)*NX*NY*NRP1+(NY-IX-1)*NX*NRP1+
     >                 IY*NRP1
                  DO 423 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 423              CONTINUE
 422            CONTINUE
 421          CONTINUE
 420        CONTINUE
          ELSE
            DO 430 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(IZ+1,3)-DAMESH(IZ,3)
              ISPLTR(IZ,3,IG)=ISPLT(IZ,3)
              DO 431 IY=0,NX-1
                DO 432 IX=0,NY-1
                  ITMI=IZ*NX*NY*NRP1+IY*NY*NRP1+
     >                 IX*NRP1
                  IDMI=IZ*NX*NY*NRP1+(NY-IX-1)*NX*NRP1+
     >                 IY*NRP1
                  DO 433 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 433              CONTINUE
 432            CONTINUE
 431          CONTINUE
 430        CONTINUE
          ENDIF
        ELSE IF(ABS(IKT) .EQ. 5) THEN
*----
*  REFLECTION WITH RESPECT TO AXIS  // TO Y
*----
          ITM(1,IG)=1
          ITM(2,IG)=2
          DO 500 IX=0,NX-1
            DAMESR(IX,1,IG)=DAMESH(NX-IX,1)-DAMESH(NX-IX-1,1)
            ISPLTR(IX,1,IG)=ISPLT(NX-IX-1,1)
 500      CONTINUE
          DAMESR(-1,1,IG)=-DAMESH(-1,1)
          DO 510 IY=0,NY-1
            DAMESR(IY,2,IG)=DAMESH(IY+1,2)-DAMESH(IY,2)
            ISPLTR(IY,2,IG)=ISPLT(IY,2)
 510      CONTINUE
          DAMESR(-1,2,IG)=DAMESH(-1,2)
          IF(IKT .LT. 0) THEN
            DO 520 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(NZ-IZ,3)-DAMESH(NZ-IZ-1,3)
              ISPLTR(IZ,3,IG)=ISPLT(NZ-IZ-1,3)
              DO 521 IY=0,NY-1
                DO 522 IX=0,NX-1
                  ITMI=IZ*NX*NY*NRP1+IY*NX*NRP1+
     >                 IX*NRP1
                  IDMI=(NZ-IZ-1)*NX*NY*NRP1+IY*NX*NRP1+
     >                 (NX-IX-1)*NRP1
                  DO 523 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 523              CONTINUE
 522            CONTINUE
 521          CONTINUE
 520        CONTINUE
          ELSE
            DO 530 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(IZ+1,3)-DAMESH(IZ,3)
              ISPLTR(IZ,3,IG)=ISPLT(IZ,3)
              DO 531 IY=0,NY-1
                DO 532 IX=0,NX-1
                  ITMI=IZ*NX*NY*NRP1+IY*NX*NRP1+
     >                 IX*NRP1
                  IDMI=IZ*NX*NY*NRP1+IY*NX*NRP1+
     >                 (NX-IX-1)*NRP1
                  DO 533 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 533              CONTINUE
 532            CONTINUE
 531          CONTINUE
 530        CONTINUE
          ENDIF
        ELSE IF(ABS(IKT) .EQ. 6) THEN
*----
*  ROTATION OF PI/2 FOLLOWED BY
*  REFLECTION WITH RESPECT TO AXIS  // TO Y
*----
          ITM(1,IG)=2
          ITM(2,IG)=1
          DO 600 IX=0,NY-1
            DAMESR(IX,1,IG)=DAMESH(IX+1,2)-DAMESH(IX,2)
            ISPLTR(IX,1,IG)=ISPLT(IX,2)
 600      CONTINUE
          DAMESR(-1,1,IG)=DAMESH(-1,2)
          DO 610 IY=0,NX-1
            DAMESR(IY,2,IG)=DAMESH(IY+1,1)-DAMESH(IY,1)
            ISPLTR(IY,2,IG)=ISPLT(IY,1)
 610      CONTINUE
          DAMESR(-1,2,IG)=DAMESH(-1,1)
          IF(IKT .LT. 0) THEN
            DO 620 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(NZ-IZ,3)-DAMESH(NZ-IZ-1,3)
              ISPLTR(IZ,3,IG)=ISPLT(NZ-IZ-1,3)
              DO 621 IY=0,NX-1
                DO 622 IX=0,NY-1
                  ITMI=IZ*NX*NY*NRP1+IY*NY*NRP1+
     >                 IX*NRP1
                  IDMI=(NZ-IZ-1)*NX*NY*NRP1+IX*NX*NRP1+
     >                 IY*NRP1
                  DO 623 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 623              CONTINUE
 622            CONTINUE
 621          CONTINUE
 620        CONTINUE
          ELSE
            DO 630 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(IZ+1,3)-DAMESH(IZ,3)
              ISPLTR(IZ,3,IG)=ISPLT(IZ,3)
              DO 631 IY=0,NX-1
                DO 632 IX=0,NY-1
                  ITMI=IZ*NX*NY*NRP1+IY*NY*NRP1+
     >                 IX*NRP1
                  IDMI=IZ*NX*NY*NRP1+IX*NX*NRP1+
     >                 IY*NRP1
                  DO 633 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 633              CONTINUE
 632            CONTINUE
 631          CONTINUE
 630        CONTINUE
          ENDIF
        ELSE IF(ABS(IKT) .EQ. 7) THEN
*----
*  REFLECTION WITH RESPECT TO AXIS // TO X
*----
          ITM(1,IG)=1
          ITM(2,IG)=2
          DO 700 IX=0,NX-1
            DAMESR(IX,1,IG)=DAMESH(IX+1,1)-DAMESH(IX,1)
            ISPLTR(IX,1,IG)=ISPLT(IX,1)
 700      CONTINUE
          DAMESR(-1,1,IG)=DAMESH(-1,1)
          DO 710 IY=0,NY-1
            DAMESR(IY,2,IG)=DAMESH(NY-IY,2)-DAMESH(NY-IY-1,2)
            ISPLTR(IY,2,IG)=ISPLT(NY-IY-1,2)
 710      CONTINUE
          DAMESR(-1,2,IG)=-DAMESH(-1,2)
          IF(IKT .LT. 0) THEN
            DO 720 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(NZ-IZ,3)-DAMESH(NZ-IZ-1,3)
              ISPLTR(IZ,3,IG)=ISPLT(NZ-IZ-1,3)
              DO 721 IY=0,NY-1
                ITMI=IZ*NX*NY*NRP1+IY*NX*NRP1
                IDMI=(NZ-IZ-1)*NX*NY*NRP1+(NY-IY-1)*NX*NRP1
                DO 722 IX=0,NX-1
                  DO 723 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 723              CONTINUE
 722            CONTINUE
 721          CONTINUE
 720        CONTINUE
          ELSE
            DO 730 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(IZ+1,3)-DAMESH(IZ,3)
              ISPLTR(IZ,3,IG)=ISPLT(IZ,3)
              DO 731 IY=0,NY-1
                ITMI=IZ*NX*NY*NRP1+IY*NX*NRP1
                IDMI=IZ*NX*NY*NRP1+(NY-IY-1)*NX*NRP1
                DO 732 IX=0,NX-1
                  DO 733 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 733              CONTINUE
 732            CONTINUE
 731          CONTINUE
 730        CONTINUE
          ENDIF
        ELSE IF(ABS(IKT) .EQ. 8) THEN
*----
*  ROTATION OF PI/2 FOLLOWED BY
*  REFLECTION WITH RESPECT TO AXIS // TO X
*----
          ITM(1,IG)=2
          ITM(2,IG)=1
          DO 800 IX=0,NY-1
            DAMESR(IX,1,IG)=DAMESH(NY-IX,2)-DAMESH(NY-IX-1,2)
            ISPLTR(IX,1,IG)=ISPLT(NY-IX-1,2)
 800      CONTINUE
          DAMESR(-1,1,IG)=-DAMESH(-1,2)
          DO 810 IY=0,NX-1
            DAMESR(IY,2,IG)=DAMESH(NX-IY,1)-DAMESH(NX-IY-1,1)
            ISPLTR(IY,2,IG)=ISPLT(NX-IY-1,1)
 810      CONTINUE
          DAMESR(-1,2,IG)=-DAMESH(-1,1)
          IF(IKT .LT. 0) THEN
            DO 820 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(NZ-IZ,3)-DAMESH(NZ-IZ-1,3)
              ISPLTR(IZ,3,IG)=ISPLT(NZ-IZ-1,3)
              DO 821 IY=0,NX-1
                DO 822 IX=0,NY-1
                  ITMI=IZ*NX*NY*NRP1+IY*NY*NRP1+
     >                 IX*NRP1
                  IDMI=(NZ-IZ-1)*NX*NY*NRP1+(NY-IX-1)*NX*NRP1+
     >                 (NX-IY-1)*NRP1
                  DO 823 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 823              CONTINUE
 822            CONTINUE
 821          CONTINUE
 820        CONTINUE
          ELSE
            DO 830 IZ=0,NZ-1
              DAMESR(IZ,3,IG)=DAMESH(IZ+1,3)-DAMESH(IZ,3)
              ISPLTR(IZ,3,IG)=ISPLT(IZ,3)
              DO 831 IY=0,NX-1
                DO 832 IX=0,NY-1
                  ITMI=IZ*NX*NY*NRP1+IY*NY*NRP1+
     >                 IX*NRP1
                  IDMI=IZ*NX*NY*NRP1+(NY-IX-1)*NX*NRP1+
     >                 (NX-IY-1)*NRP1
                  DO 833 IR=0,NMR
                    MIXR(ITMI,IG,1)=MIX(IDMI,1)
                    MIXR(ITMI,IG,2)=MIX(IDMI,2)
                    ITMI=ITMI+1
                    IDMI=IDMI+1
 833              CONTINUE
 832            CONTINUE
 831          CONTINUE
 830        CONTINUE
          ENDIF
        ENDIF
        IF(IPRINT .GE. 100) THEN
*----
*  Print turned mesh if required
*----
          WRITE(IOUT,6010) (NM(ITM(IDIR,IG)),IDIR=1,3),NREG
          DO IDIR=1,4
            NMTMP=NM(ITM(IDIR,IG))
            IF(NMTMP .GT. 0) THEN
              WRITE(IOUT,6011) 'MESH'//CDIR(IDIR)//' ='
              WRITE(IOUT,6012) (DAMESR(IX,IDIR,IG),IX=-1,NMTMP)
              WRITE(IOUT,6011) 'SPLT'//CDIR(IDIR)//' ='
              WRITE(IOUT,6013) (ISPLTR(IX-1,IDIR,IG),IX=1,NMTMP)
            ENDIF
          ENDDO
          WRITE(IOUT,6011) 'MIX   ='
          WRITE(IOUT,6013) (MIXR(IX,IG,1),IX=0,NREG-1)
          WRITE(IOUT,6011) 'HMIX  ='
          WRITE(IOUT,6013) (MIXR(IX,IG,2),IX=0,NREG-1)
        ENDIF
        IF(IPG .GT. 1) THEN
*----
*  COMPARE GEOMETRY
*  1- MESH AND SPLIT IN X, Y AND Z
*  2- MIXTURES
*  3- OFFCENTER
*----
          DO 900 IDIR=1,3
            NMTMP=NM(ITM(IDIR,1))
            IF(NMTMP .NE. NM(ITM(IDIR,2))) CALL XABORT(NAMSBR//
     >      ': Symmetry invalid with this mesh')
            DO 910 IX=0,NMTMP-1
              DDD=ABS(DAMESR(IX,IDIR,1)-DAMESR(IX,IDIR,2))
              IF(DDD .GT. DCUTOF) CALL XABORT(NAMSBR//
     >        ': Symmetry invalid with this mesh')
              IF(ISPLTR(IX,IDIR,1) .NE. ISPLTR(IX,IDIR,2) )
     >        CALL XABORT(NAMSBR//
     >        ': Symmetry invalid with this split')
 910        CONTINUE
 900      CONTINUE
          DO 920 IX=0,NREG-1
            IF(MIXR(IX,1,1) .NE. MIXR(IX,2,1) ) CALL XABORT(NAMSBR//
     >      ': Symmetry invalid with this mixture')
            IF(MIXR(IX,1,2) .NE. MIXR(IX,2,2) ) CALL XABORT(NAMSBR//
     >      ': Symmetry invalid with this merging mixture')
 920      CONTINUE
          IF(DAMESR(-1,1,1) .NE. DAMESR(-1,1,2) .OR.
     >       DAMESR(-1,2,1) .NE. DAMESR(-1,2,2) .OR.
     >       DAMESR(-1,3,1) .NE. DAMESR(-1,3,2) ) CALL XABORT(NAMSBR//
     >    ': Symmetry invalid with this off center')
        ELSE
*----
*  Reset reference geometry for turn
*----
          DO IX=0,NR-1
            DAMESH(IX,4)=DAMESR(IX,4,IG)
            ISPLT(IX,4)=ISPLTR(IX,4,IG)
          ENDDO
          DAMESH(NR,4)=DAMESR(NR,4,IG)
          DAMESH(-1,4)=DAMESR(-1,4,IG)
*----
*  Find splitted mesh dimensions
*----
          DO 930 IDIR=1,4
            NMTMP=NM(ITM(IDIR,1))
            NMT(IDIR)=NMTMP
            NMTS(IDIR)=0
            DO 931 IX=0,NMTMP-1
              NMTS(IDIR)=NMTS(IDIR)+ABS(ISPLTR(IX,IDIR,1))
 931        CONTINUE
            IF(NMTS(IDIR) .NE. NMS(ITM(IDIR,1))) CALL XABORT(NAMSBR//
     >       ': Global symmetry invalid with this split')
 930      CONTINUE
        ENDIF
 1005   CONTINUE
      ENDDO
*----
*  Reset final mesh (center+original turn)
*----
      DO IDIR=1,3
        NMTMP=NMT(IDIR)
        CALL XDDSET(DAMESH(-1,IDIR),NM(IDIR)+2,DZERO)
        NM(IDIR)=NMTMP
        DDD=DZERO
        DO IX=0,NMTMP-1
          DDD=DDD+DAMESR(IX,IDIR,1)
        ENDDO
        DDD=DDD/DTWO
        DAMESH(-1,IDIR)=DAMESR(-1,IDIR,1)
        DAMESH(0,IDIR)=-DDD
        DO IX=1,NMTMP
          DAMESH(IX,IDIR)=DAMESH(IX-1,IDIR)+DAMESR(IX-1,IDIR,1)
        ENDDO
        DO IX=0,NMTMP
          ISPLT(IX,IDIR)=ISPLTR(IX,IDIR,1)
        ENDDO
      ENDDO
      DO IDIR=1,4
        NMTMP=NM(IDIR)
        NMS(IDIR)=0
        DO IX=0,NMTMP-1
          NMS(IDIR)=NMS(IDIR)+ABS(ISPLT(IX,IDIR))
        ENDDO
      ENDDO
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  FORMATS
*----
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(1X,' DIMENSIONS =',5I10/1X,' ORIGINAL MESH ')
 6011 FORMAT(1X,A7)
 6012 FORMAT(5F15.9)
 6013 FORMAT(5I15)
      END
