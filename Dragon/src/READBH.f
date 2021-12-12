*DECK READBH
      SUBROUTINE READBH (MAXPTS,IPGEOM,IR,IR2,NREG,NREG2,MAT,VOL,NG,
     1 NSMAX,MICRO,NS,IBI,RS,FRACT,VOLK,IMPX,IDIL,MIXGR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the input data for the double heterogeneity option (Bihet).
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
* MAXPTS  allocated storage for arrays of dimension NREG.
* IPGEOM  pointer to the geometry LCM object (L_GEOM signature).
* IR2     number of ordinary and composite mixtures (a mixture that
*         include a micro structure).
* NREG2   number of volumes in the macro geometry.
* NG      number of different kind of micro structures. A kind of
*         micro structure is characterized by the radius of its
*         tubes or shells. All the micro volumes of the same kind
*         should own the same nuclear properties in a given volume
*         of the macro geometry.
* NSMAX   maximum number of volumes (tubes or shells) in each kind
*         of micro structure).
* IBI     type of composite mixture in each volume of the macro
*         geometry. If IBI(IKK) is greater than IR, the volume IKK
*         contains a micro structure.
* IMPX    print flag (equal to zero for no print).
*
*Parameters: input/output
* VOL     volumes of the macro geometry on input and
*         volumes of the composite geometry at output.
*
*Parameters: output
* IR      number of ordinary mixtures.
* NREG    number of volumes in the composite geometry.
* MAT     index-number of the mixture type assigned to each volume
*         of the composite geometry.
* MICRO   type of micro volumes: =3 cylinder; =4 sphere.
* NS      number of tubes or shells in each kind of micro structure.
* RS      radius of the micro volumes.
* FRACT   volumic fractions of each type of micro volumes in each
*         ordinary or composite mixture.
* VOLK    volumic fractions of the tubes or shells in the micro volumes.
* IDIL    elementary mixture indices in the diluent of the composite
*         mixtures.
* MIXGR   elementary mixture indices in the micro structures.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPGEOM
      INTEGER MAXPTS,IR,IR2,NREG,NREG2,MAT(MAXPTS),NG,NSMAX,MICRO,
     1 NS(NG),IBI(NREG2),IMPX,IDIL(IR2),MIXGR(NSMAX,NG,IR2)
      REAL VOL(MAXPTS),RS(NSMAX+1,NG),FRACT(NG,IR2),VOLK(NG,NSMAX)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      CHARACTER GEONAM*12,HSMG*131,TEXT12*12
      LOGICAL EMPTY,LCM
      INTEGER ISTATE(NSTATE)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MILIEU
      REAL, ALLOCATABLE, DIMENSION(:) :: FTEMP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MILIEU(IR2))
      ALLOCATE(FTEMP(IR2))
*
      CALL XDISET(IDIL,IR2,0)
      CALL XDISET(MIXGR,NSMAX*NG*IR2,0)
      CALL LCMINF(IPGEOM,GEONAM,TEXT12,EMPTY,ILONG,LCM)
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      NMILG=ISTATE(3)
      MICRO=ISTATE(5)
      IF(NMILG.GT.IR2) CALL XABORT('READBH: INVALID VALUE FOR IR2.')
      CALL LCMGET(IPGEOM,'NS',NS)
      CALL LCMGET(IPGEOM,'RS',RS)
      CALL LCMGET(IPGEOM,'FRACT',FRACT)
      CALL LCMGET(IPGEOM,'MILIE',MILIEU)
      CALL LCMGET(IPGEOM,'MIXDIL',IDIL)
      CALL LCMGET(IPGEOM,'MIXGR',MIXGR)
*
      DO 30 J=1,NG
      FACT=RS(NS(J)+1,J)**(MICRO-1)
      W=0.0
      DO 20 K=1,NS(J)
      ZZZ=RS(K+1,J)**(MICRO-1)/FACT
      VOLK(J,K)=ZZZ-W
      W=ZZZ
   20 CONTINUE
   30 CONTINUE
*
      IND1=IR2
      DO 40 I=1,NMILG
      IF(MILIEU(I).GT.IR2) THEN
         WRITE (HSMG,390) MILIEU(I),IR2
         CALL XABORT(HSMG)
      ENDIF
      IND1=MIN(IND1,MILIEU(I))
   40 CONTINUE
*----
*  SET-UP THE NEW VOLUMES
*----
      NREG=NREG2
      DO 90 IKK=1,NREG2
      MAT(IKK)=IBI(IKK)
      IF(IBI(IKK).GE.IND1) THEN
         IND=0
         DO 50 I=1,NMILG
         IF(MILIEU(I).EQ.IBI(IKK)) IND=I
   50    CONTINUE
         IF(IND.EQ.0) THEN
            WRITE(HSMG,'(29HREADBH: A COMPOSITE MIXTURE (,I5,7H) IS NO,
     1      10HT DEFINED.)') IBI(IKK)
            CALL XABORT(HSMG)
         ENDIF
         DILF=1.0
         DO 60 J=1,NG
         DILF=DILF-FRACT(J,IND)
   60    CONTINUE
         VHET=VOL(IKK)
         MAT(IKK)=IDIL(IND)
         VOL(IKK)=VHET*DILF
         DO 80 J=1,NG
         FRT=FRACT(J,IND)
         IF(FRT.GT.0.00001) THEN
            FACT=RS(NS(J)+1,J)**(MICRO-1)
            W=0.0
            DO 70 K=1,NS(J)
            ZZZ=RS(K+1,J)**(MICRO-1)/FACT
            NREG=NREG+1
            MAT(NREG)=MIXGR(K,J,IND)
            VOL(NREG)=VHET*FRT*(ZZZ-W)
            W=ZZZ
   70       CONTINUE
          ENDIF
   80     CONTINUE
      ENDIF
   90 CONTINUE
      IF(NREG.GT.MAXPTS) CALL XABORT('READBH: MAXPTS IS TOO SMALL.')
      IR=0
      DO 100 I=1,NREG
      IR=MAX(IR,MAT(I))
  100 CONTINUE
      IF(IR+1.GT.IR2) CALL XABORT('READBH: INVALID MIX NUMBERS.')
      DO IND=1,IR2-IR
        IF(IDIL(IND).EQ.0) THEN
          WRITE(HSMG,'(15HREADBH: MIXTURE,I5,22H IS NOT USED IN THE GE,
     1    7HOMETRY.)') IR+IND-1
          CALL XABORT(HSMG)
        ENDIF
      ENDDO
      DO 135 J=1,NG
      DO 110 IND=1,NMILG
      FTEMP(IND)=FRACT(J,IND)
  110 CONTINUE
      DO 120 IND=1,IR2
      FRACT(J,IND)=0.0
  120 CONTINUE
      DO 130 IND=1,NMILG
      FRACT(J,MILIEU(IND))=FTEMP(IND)
  130 CONTINUE
  135 CONTINUE
*
      IF(IMPX.GE.1) THEN
         WRITE (6,300) GEONAM
         IF(MICRO.EQ.3) THEN
            WRITE (6,'(44H THE MICRO STRUCTURE IS MADE OF TUBES OR CYL,
     1      7HINDERS./)')
         ELSE IF(MICRO.EQ.4) THEN
            WRITE (6,'(44H THE MICRO STRUCTURE IS MADE OF SPHERES OR S,
     1      16HPHERICAL SHELLS./)')
         ENDIF
         WRITE (6,360) NREG2,NG,NSMAX
         WRITE (6,370) IR+1,IR2
         DO 140 J=1,NG
         WRITE (6,310) J
         WRITE (6,320) (RS(K,J),K=1,NS(J)+1)
         WRITE (6,330)
         WRITE (6,320) (FRACT(J,IBI(IKK)),IKK=1,NREG2)
  140    CONTINUE
         WRITE (6,'(///)')
         WRITE (6,400) NREG,MAXPTS,IR
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FTEMP)
      DEALLOCATE(MILIEU)
      RETURN
*
  300 FORMAT (///50H BIHET: INTRODUCTION OF A MICRO STRUCTURE IN THE M,
     1 26HACRO GEOMETRY LOCATED IN ',A12,2H'./)
  310 FORMAT (//23H MICRO STRUCTURE NUMBER,I4//20H RADIUS OF THE MICRO,
     1 17H TUBES OR SHELLS:)
  320 FORMAT (1X,1P,10E12.5)
  330 FORMAT (/53H VOLUMIC CONCENTRATIONS OF THE MICRO STRUCTURE IN EAC,
     1 31HH VOLUME OF THE MACRO GEOMETRY:)
  360 FORMAT (/42H NUMBER OF VOLUMES IN THE MACRO GEOMETRY =,I6/
     1         38H NUMBER OF KINDS OF MICRO STRUCTURES =,I6/
     2         49H MAXIMUM NUMBER OF VOLUMES IN A MICRO STRUCTURE =,I6/)
  370 FORMAT (/51H THE INDEX-NUMBERS OF THE MIXTURES WITH A MICRO STR,
     1 21HUCTURE VARIES BETWEEN,I6,4H AND,I6,1H.)
  390 FORMAT (34HREADBH: THE INPUT MIXTURE NUMBER (,I6,12H) IS GREATER,
     1 10H THAN IR (,I6,2H ))
  400 FORMAT (/20H NUMBER OF VOLUMES =,I6,5X,22HAVAILABLE STORAGE: MAX,
     1 5HPTS =,I6/21H NUMBER OF MIXTURES =,I6/)
      END
