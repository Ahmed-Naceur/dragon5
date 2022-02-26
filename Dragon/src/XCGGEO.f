*DECK XCGGEO
      SUBROUTINE XCGGEO(IPGEOM,IROT,NSOUT,NVOL,NBAN,MNAN,NRT,MSROD,
     >                  IPRT,ILK,NMAT,RAN,NRODS,RODS,NRODR,RODR,NRINFO,
     >                  MATALB,VOLSUR,COTE,RADMIN,NCODE,ICODE,ZCODE,
     >                  ALBEDO,KEYMRG,NXRS,NXRI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read and analyse 2-D cluster geometry.
*
*Copyright:
* Copyright (C) 1990 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPGEOM  pointer to the geometry.
* IROT    type of pij reconstruction:
*          <0 cp calculations with symmetries;
*          =0 cp calculations;
*          =1 direct jpm reconstruction;
*          =2 rot2 type reconstruction.
* NSOUT   number of outer surface.
* NVOL    maximum number of regions.
* NBAN    number of concentric regions.
* MNAN    maximum number of radius to read.
* NRT     number of rod types.
* MSROD   maximum number of subrods per rods.
* IPRT    impression level.
*
*Parameters: output
* ILK     leakage flag. ILK=.TRUE. if neutron leakage through
*         external boundary is present.
* NMAT    total number of materials.
* RAN     radius of annular regions.
* NRODS   integer description of rod of a given type:
*         NRODS(1,IRT) = number of rod;
*         NRODS(2,IRT) = number of subrods in rod;
*         NRODS(3,IRT) = first concentric region.
* RODS    real description of rod of a given type:
*         RODS(1,IRT) = rod center radius;
*         RODS(2,IRT) = angular position of first rod.
* NRODR   subrod region.
* RODR    subrod radius.
* NRINFO  annular region content.
*         NRINFO(1,IAN) = new region number.
*         NRINFO(2,IAN): = +I cluster number        (all);
*                       = 1000000+I cluster number cut (in);
*                       = 2000000+I cluster number cut (part);
*                       = 3000000+I cluster number cut (out);
*                       = 0 no cluster associated;
*                       = -I cluster at center     (all).
* MATALB  albedo-material of regions.
* VOLSUR  surface/4-volume of regions.
* COTE    additional side length for rectangle.
* RADMIN  minimum radius of region.
* NCODE   albedo type.
* ICODE   albedo number associated with face.
* ZCODE   albedo zcode vector.
* ALBEDO  albedo.
* KEYMRG  region-surface merge vector.
* NXRS    integer description of rod of a given type
*         last concentric region.
* NXRI    annular region content multi-rod.
*
*----------------------------------------------------------------------
*
      USE        GANLIB
      IMPLICIT   NONE
      INTEGER    IOUT,NSTATE,NMCOD
      REAL       PI,THSQ3
      PARAMETER (IOUT=6,NSTATE=40,NMCOD=6,PI=3.1415926535898,
     >           THSQ3=2.598076212)
      CHARACTER  NAMSBR*6
      PARAMETER (NAMSBR='XCGGEO')
*-----
*  ROUTINE PARAMETERS
*----
      TYPE(C_PTR) IPGEOM
      LOGICAL    ILK,EMPTY,LCM
      INTEGER    IROT,NSOUT,NVOL,NBAN,MNAN,NRT,MSROD,IPRT,
     >           NMAT,NRODS(3,NRT),NRODR(NRT),NRINFO(2,NBAN),
     >           MATALB(-NSOUT:NVOL),NCODE(NMCOD),ICODE(NMCOD),
     >           KEYMRG(-NSOUT:NVOL),NXRS(NRT),NXRI(NRT,NBAN)
      REAL       RAN(NBAN),RODS(2,NRT),RODR(MSROD,NRT),
     >           VOLSUR(-NSOUT:NVOL),COTE,RADMIN,ALBEDO(NMCOD),
     >           ZCODE(NMCOD)
*----
*  LOCAL VARIABLES
*----
      INTEGER    ISTATE(NSTATE)
      CHARACTER  GEONAM*12,TEXT12*12,CMSG*131
      INTEGER    IRT,IAN,IS,IC,ITRAN,I,NRANN,NRRANN,NSPLIT,ISA,ILSTP,
     >           ISPL,ISURW,NTAN,IPOS,ITYPE,IM,ISR,IZRT,JAN,JRT,KRT,
     >           ILR,JSUR,JSW,ISV,ILSTR,JPRT,LRT,IREG,ILONG
      REAL       RADL,RADN,VFIN,DELV,XTOP,XBOT,VOLI,VOLROD,VOLF,
     >           VOLIS,XNROD,VOLFS,VANSPI,VRPSPI,VRDSPI,XINT,
     >           YINT,ANGR,ANGA,VRGOU1,VRGIN1
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MATANN,ISPLIT,JGEOM
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: MATROD
      REAL, ALLOCATABLE, DIMENSION(:) :: RAD
      REAL, ALLOCATABLE, DIMENSION(:,:) :: VRGIO
*----
*  SCRATCH STORAGE ALLOCATION
*   MATANN  : TYPE OF MATERIAL FOR ANNULAR REGIONS        I(NBAN)
*   MATROD  : TYPE OF MATERIAL FOR EACH SUBROD            I(MSROD,NRT)
*   ISPLIT  : SPLITTING VECTOR FOR RODS                   I(NBAN)
*   RAD     : RADIUS VECTOR                               R(MNAN)
*   VRGIO   : DIVIDED ROD VOLUME                          R(2,NRT)
*           : 2 - INSIDE REGION
*           : 1 - OUTSIDE REGION
*----
      ALLOCATE(MATANN(NBAN),MATROD(MSROD,NRT),ISPLIT(NBAN))
      ALLOCATE(RAD(MNAN),VRGIO(2,NRT))
*----
*  INITIALIZE NRINFO, NXRI AND NXRS TO 0
*----
      DO 3 IRT=1,NRT
        NXRS(IRT)=0
        NRODR(IRT)=0
 3    CONTINUE
      DO 4 IAN=1,NBAN
        NRINFO(1,IAN)=0
        NRINFO(2,IAN)=0
        DO 5 IRT=1,NRT
          NXRI(IRT,IAN)=0
 5      CONTINUE
 4    CONTINUE
      DO 6 IS=-NSOUT,NVOL
        KEYMRG(IS)=IS
 6    CONTINUE
      VOLSUR(0)=0.0
      MATALB(0)=0
*----
*  READ GEOMETRY INFORMATIONS
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
*----
*  RECOVER THE BOUNDARY CONDITIONS.
*----
      CALL LCMGET(IPGEOM,'NCODE',NCODE)
      CALL LCMGET(IPGEOM,'ZCODE',ALBEDO)
      CALL LCMGET(IPGEOM,'ICODE',ICODE)
      DO 7 IC=1,NMCOD
        ZCODE(IC)=ALBEDO(IC)
        IF(ICODE(IC).NE.0) CALL XABORT(NAMSBR//
     >': MACROLIB DEFINED ALBEDOS ARE NOT IMPLEMENTED.')
 7    CONTINUE
      ITRAN=0
      DO 100 I=1,NMCOD
        IF ((NCODE(I) .EQ. 3) .OR. (NCODE(I) .EQ. 5) .OR.
     >      (NCODE(I) .GE. 7)) THEN
          CALL XABORT(NAMSBR//': INVALID TYPE OF B.C.')
        ELSE IF(NCODE(I).EQ.2) THEN
          ZCODE(I)=1.0
          ALBEDO(I)=1.0
        ELSE IF(NCODE(I).EQ.4) THEN
          ITRAN=ITRAN+1
          ZCODE(I)=1.0
          ALBEDO(I)=1.0
        ELSE IF(NCODE(I).EQ.6) THEN
           NCODE(I)=1
        ENDIF
 100  CONTINUE
      IF(NSOUT.EQ.1.AND.IROT.GT.-400) THEN
        MATALB(-1)=-2
        IF (NCODE(2).EQ.0)
     >    CALL XABORT(NAMSBR//': ANNULAR BOUNDARY CONDITION MISSING.')
        IF(ITRAN.NE.0) THEN
           NCODE(2)=2
        ENDIF
        IF(ICODE(2).EQ.0) ICODE(2)=-2
        ILK=( (NCODE(2).EQ.1) .OR. (ZCODE(2).NE.1.0) )
      ELSE IF(NSOUT.EQ.6.OR.IROT.LT.-600) THEN
        IF(ITRAN.NE.0) THEN
           NCODE(1)=2
        ENDIF
        IF(IROT.LT.0) THEN
          MATALB(-1)=-1
        ELSE
          MATALB(-1)=-1
          DO 101 IS=2,6
            ZCODE(IS)=ZCODE(1)
            MATALB(-IS)=-1
 101      CONTINUE
        ENDIF
        IF (NCODE(1).EQ.0) CALL XABORT(NAMSBR//
     >    ': HEXAGONAL BOUNDARY CONDITION MISSING.')
        IF(ICODE(1).EQ.0) ICODE(1)=-1
        ILK=( (NCODE(1).EQ.1) .OR. (ZCODE(1).NE.1.0) )
      ELSE
        IF(IROT.LT.0) THEN
          IF(ITRAN.NE.0) CALL XABORT(NAMSBR//
     >    ': CARTESIAN SYMMETRY NO TRANSLATION BOUNDARY CONDITIONS')
          IF(ZCODE(1).NE.ZCODE(2).OR.ZCODE(1).NE.ZCODE(3).OR.
     >       ZCODE(1).NE.ZCODE(4)) CALL XABORT(NAMSBR//
     >    ': CARTESIAN SYMMETRY REQUIRES '//
     >    ' IDENTICAL BOUNDARY CONDITION IN ALL DIRECTIONS.')
          MATALB(-1)=-1
          IF (NCODE(1).EQ.0) CALL XABORT(NAMSBR//
     >      ': CARTESIAN BOUNDARY CONDITION MISSING.')
          IF(ICODE(1).EQ.0) ICODE(1)=-1
          ILK=( (NCODE(1).EQ.1) .OR. (ZCODE(1).NE.1.0) )
        ELSE
          MATALB(-1)=-2
          MATALB(-2)=-4
          MATALB(-3)=-1
          MATALB(-4)=-3
          ZCODE(5)=ZCODE(1)
          ZCODE(1)=ZCODE(2)
          ZCODE(2)=ZCODE(4)
          ZCODE(4)=ZCODE(3)
          ZCODE(3)=ZCODE(5)
          ILK=.FALSE.
          DO 102 IS=1,NSOUT
            IF (NCODE(IS).EQ.0) CALL XABORT(NAMSBR//
     >      ': RECTANGLE BOUNDARY CONDITION MISSING.')
            IF(.NOT. ILK) THEN
              IF( (NCODE(IS).EQ.1) .OR. (ZCODE(IS).NE.1.0) ) THEN
                ILK=.TRUE.
              ENDIF
            ENDIF
            IF(ICODE(IS).EQ.0) ICODE(IS)=-IS
 102      CONTINUE
          IF(ITRAN .GT. 0) THEN
            IF(MOD(ITRAN,2) .EQ. 1) CALL XABORT(NAMSBR//
     >        ': TRANSLATION SYMMETRIES COME IN PAIRS')
            IF((NCODE(1) .EQ. 4) .AND. (NCODE(2) .EQ. 4)) THEN
              ITRAN=ITRAN-2
            ENDIF
            IF((NCODE(3) .EQ. 4) .AND. (NCODE(4) .EQ. 4)) THEN
              ITRAN=ITRAN-2
            ENDIF
            IF(ITRAN .NE. 0) CALL XABORT(NAMSBR//
     >        ': WRONG PAIRS OF TRANSLATION SYMMETRIES')
          ENDIF
        ENDIF
      ENDIF
*----
*  RECOVER THE MIXTURE FOR ANNULAR REGIONS
*----
      NRANN=ISTATE(6)
      CALL LCMGET(IPGEOM,'MIX',MATANN)
      NMAT=0
      DO 110 I=1,NRANN
        NMAT=MAX(NMAT,MATANN(I))
 110  CONTINUE
*----
*  RECOVER THE MESH COORDINATES
*----
      IF((IROT.LT.-400).OR.(NSOUT.GT.1)) THEN
        NRRANN=NRANN-1
        MATANN(NBAN)=MATANN(NRANN)
      ELSE
        NRRANN=NRANN
      ENDIF
      CALL LCMGET(IPGEOM,'RADIUS',RAD)
      IF(ISTATE(11).EQ.1) THEN
*----
*  SPLIT ANNULUS WHEN REQUIRED
*----
        CALL LCMLEN(IPGEOM,'SPLITR',ILONG,ITYPE)
        IF(ILONG.GT.NBAN) CALL XABORT(NAMSBR//': SPLITR OVERFLOW')
        CALL LCMGET(IPGEOM,'SPLITR',ISPLIT)
        NSPLIT=0
        DO 145 ISA=1,NRRANN
          NSPLIT=NSPLIT+ABS(ISPLIT(ISA))
 145    CONTINUE
        ILSTP=NSPLIT
        RADL=RAD(NRRANN+1)
        DO 155 ISA=NRRANN,1,-1
          RADN=RAD(ISA)
          RAN(ILSTP)=RADL
          MATANN(ILSTP)=MATANN(ISA)
          IF(ISPLIT(ISA).LT.0) THEN
*----
*  ANNULUS EQUAL VOLUMES SPLIT
*----
            VFIN=RADL*RADL
            DELV=(VFIN-RADN*RADN)/FLOAT(ABS(ISPLIT(ISA)))
            DO 165 ISPL=ABS(ISPLIT(ISA))-1,1,-1
              ILSTP=ILSTP-1
              VFIN=VFIN-DELV
              RAN(ILSTP)=SQRT(VFIN)
              MATANN(ILSTP)=MATANN(ISA)
 165        CONTINUE
          ELSE IF(ISPLIT(ISA).GT.0) THEN
*----
*  ANNULUS EQUAL TICKNESS SPLIT
*----
            VFIN=RADL
            DELV=(VFIN-RADN)/FLOAT(ISPLIT(ISA))
            DO 175 ISPL=ISPLIT(ISA)-1,1,-1
              ILSTP=ILSTP-1
              VFIN=VFIN-DELV
              RAN(ILSTP)=VFIN
              MATANN(ILSTP)=MATANN(ISA)
 175        CONTINUE
          ELSE
            CALL XABORT(NAMSBR//': A SPLIT OF 0 IS INVALID')
          ENDIF
          RADL=RADN
          ILSTP=ILSTP-1
 155    CONTINUE
      ELSE
        DO 20 IAN=1,NRRANN
          RAN(IAN)=RAD(IAN+1)
 20     CONTINUE
      ENDIF
      RADMIN=RAN(1)
      NTAN=NBAN
      IF(NSOUT.EQ.1.AND.IROT.GT.-400) THEN
        VOLSUR(-1)=0.5*PI*RAN(NBAN)
      ELSE IF(NSOUT.EQ.6.OR.IROT.LT.-600) THEN
        CALL LCMGET(IPGEOM,'SIDE',RAN(NBAN))
        NTAN=NBAN-1
        IF(IROT.LT.0) THEN
          VOLSUR(-1)=1.5*RAN(NBAN)
        ELSE
          VOLSUR(-1)=0.25*RAN(NBAN)
          DO 30 ISURW=-2,-6,-1
            VOLSUR(ISURW)=VOLSUR(-1)
 30       CONTINUE
        ENDIF
      ELSE
        CALL LCMGET(IPGEOM,'MESHX',RAD(1))
        CALL LCMGET(IPGEOM,'MESHY',RAD(3))
        RAN(NBAN)=RAD(2)-RAD(1)
        COTE=RAD(4)-RAD(3)
        NTAN=NBAN-1
        IF(IROT.LT.0) THEN
          IF(RAN(NBAN).NE.COTE) CALL XABORT(NAMSBR//
     >    ': CARTESIAN SYMMETRY REQUIRES  SQUARE CELL.')
          VOLSUR(-1)=COTE
        ELSE
          VOLSUR(-1)=0.25*COTE
          VOLSUR(-2)=0.25*RAN(NBAN)
          VOLSUR(-3)=VOLSUR(-1)
          VOLSUR(-4)=VOLSUR(-2)
        ENDIF
      ENDIF
*----
*  READ CLUSTER GEOMETRY AND ANALYSE
*----
      ALLOCATE(JGEOM(3*NRT))
      IPOS=1
      CALL LCMGET(IPGEOM,'CLUSTER',JGEOM)
*----
*  READ ROD DESCRIPTION AND SAVE
*----
      DO 120 IRT=1,NRT
        WRITE(TEXT12(1:4),'(A4)')  JGEOM(IPOS)
        WRITE(TEXT12(5:8),'(A4)')  JGEOM(IPOS+1)
        WRITE(TEXT12(9:12),'(A4)') JGEOM(IPOS+2)
        IPOS=IPOS+3
        CALL LCMSIX(IPGEOM,TEXT12,1)
        CALL XDISET(ISTATE,NSTATE,0)
        CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
        CALL LCMGET(IPGEOM,'MIX',MATROD(1,IRT))
        CALL LCMLEN(IPGEOM,'RADIUS',NRODS(2,IRT),ITYPE)
        CALL LCMGET(IPGEOM,'NPIN',NRODS(1,IRT))
        CALL LCMGET(IPGEOM,'RPIN',RODS(1,IRT))
        CALL LCMGET(IPGEOM,'APIN',RODS(2,IRT))
        NRODS(2,IRT)=NRODS(2,IRT)-1
        CALL LCMGET(IPGEOM,'RADIUS',RAD)
        DO 121 IM=1,NRODS(2,IRT)
          RODR(IM,IRT)=RAD(IM+1)
          NMAT=MAX(NMAT,MATROD(IM,IRT))
 121    CONTINUE
        IF(ISTATE(11).EQ.1) THEN
*----
*  SPLIT RODS WHEN REQUIRED
*----
          CALL LCMLEN(IPGEOM,'SPLITR',ILONG,ITYPE)
          IF(ILONG.GT.NBAN) CALL XABORT(NAMSBR//': SPLITR OVERFLOW')
          CALL LCMGET(IPGEOM,'SPLITR',ISPLIT)
          NSPLIT=0
          DO 140 ISR=1,NRODS(2,IRT)
            NSPLIT=NSPLIT+ABS(ISPLIT(ISR))
 140      CONTINUE
          ILSTP=NSPLIT
          RADL=RODR(NRODS(2,IRT),IRT)
          DO 150 ISR=NRODS(2,IRT),1,-1
            IF(ISR.EQ.1) THEN
              RADN=0.0
            ELSE
              RADN=RODR(ISR-1,IRT)
            ENDIF
            IF(ISPLIT(ISR).LT.0) THEN
*----
*  RODS EQUAL VOLUMES SPLIT
*----
              RODR(ILSTP,IRT)=RADL
              MATROD(ILSTP,IRT)=MATROD(ISR,IRT)
              VFIN=RADL*RADL
              DELV=(VFIN-RADN*RADN)/FLOAT(ABS(ISPLIT(ISR)))
              DO 160 ISPL=ABS(ISPLIT(ISR))-1,1,-1
                ILSTP=ILSTP-1
                VFIN=VFIN-DELV
                RODR(ILSTP,IRT)=SQRT(VFIN)
                MATROD(ILSTP,IRT)=MATROD(ISR,IRT)
 160          CONTINUE
            ELSE IF(ISPLIT(ISR).GT.0) THEN
*----
*  RODS EQUAL TICKNESS SPLIT
*----
              RODR(ILSTP,IRT)=RADL
              MATROD(ILSTP,IRT)=MATROD(ISR,IRT)
              VFIN=RADL
              DELV=(VFIN-RADN)/FLOAT(ISPLIT(ISR))
              DO 170 ISPL=ISPLIT(ISR)-1,1,-1
                ILSTP=ILSTP-1
                VFIN=VFIN-DELV
                RODR(ILSTP,IRT)=VFIN
                MATROD(ILSTP,IRT)=MATROD(ISR,IRT)
 170          CONTINUE
            ELSE
              CALL XABORT(NAMSBR//': A SPLIT OF 0 IS INVALID')
            ENDIF
            RADL=RADN
            ILSTP=ILSTP-1
 150      CONTINUE
          NRODS(2,IRT)=NSPLIT
        ENDIF
        RADMIN=MIN(RADMIN,RODR(1,IRT))
        CALL LCMSIX(IPGEOM,' ',2)
 120  CONTINUE
*----
*  CHECK ROD GEOMETRY AND REORDER IF NECESSARY
*----
      CALL XCGROD(NRT,MSROD,NRODS,RODS,MATROD,RODR)
*----
*  LOCALIZE ROD POSITION WITH RESPECT TO ANNULUS
*----
      IZRT=0
      DO 122 IRT=1,NRT
        XTOP=RODS(1,IRT)+RODR(NRODS(2,IRT),IRT)
        XBOT=RODS(1,IRT)-RODR(NRODS(2,IRT),IRT)
        IF((XBOT.LT.0.0).AND.(NRODS(1,IRT).GT.1)) THEN
          CALL XABORT(NAMSBR//': OVERLAPPING RODS')
        ELSE IF(RODS(1,IRT).EQ.0.0) THEN
          IF(RODR(NRODS(2,IRT),IRT).LE.RAN(1)) THEN
            NRODS(3,IRT)=-1
            NXRS(IRT)=-1
            NXRI(IRT,1)=-1
            NRINFO(2,1)=-IRT
          ELSE
            CALL XABORT(NAMSBR//': CENTRAL ROD OVERLAPP WITH ANNULUS')
          ENDIF
        ELSE
*----
*  SEARCH IN ANNULUS SINCE RODS MAY NOT BE LOCATED IN
*  SQUARE OR HEXAGONAL CROWN WHERE NTAN=NBAN-1
*----
          JAN=0
          KRT=0
          DO 130 IAN=1,NTAN
            JAN=IAN
            IF(XBOT.LE.RAN(IAN)) THEN
              NRODS(3,IRT)=IAN
              NXRS(IRT)=IAN
              NRINFO(2,IAN)=IRT
              DO 134 JRT=1,NRT
                KRT=JRT
                IF(NXRI(KRT,IAN).EQ.0) THEN
                  NXRI(KRT,IAN)=IRT
                  IZRT=MAX(IZRT,KRT)
                  GO TO 131
                ENDIF
 134          CONTINUE
            ENDIF
 130      CONTINUE
          WRITE(CMSG,9001) NAMSBR,IRT
          CALL XABORT(CMSG)
 131      CONTINUE
          IF(XTOP.GT.RAN(JAN)) THEN
            NXRI(KRT,IAN)=IRT+1000000
            DO 132 IAN=JAN+1,NTAN
              IF(XTOP.LE.RAN(IAN)) THEN
                NXRS(IRT)=IAN
                NRINFO(2,IAN)=IRT
                DO 135 JRT=1,NRT
                  KRT=JRT
                  IF(NXRI(KRT,IAN).EQ.0) THEN
                    NXRI(KRT,IAN)=IRT+3000000
                    IZRT=MAX(IZRT,KRT)
                    GO TO 133
                  ENDIF
 135            CONTINUE
              ELSE
                NXRS(IRT)=IAN
                NRINFO(2,IAN)=IRT
                DO 136 JRT=1,NRT
                  KRT=JRT
                  IF(NXRI(KRT,IAN).EQ.0) THEN
                    NXRI(KRT,IAN)=IRT+2000000
                    IZRT=MAX(IZRT,KRT)
                    GO TO 137
                  ENDIF
 136            CONTINUE
 137            CONTINUE
              ENDIF
 132        CONTINUE
            WRITE(CMSG,9001) NAMSBR,IRT
            CALL XABORT(CMSG)
 133        CONTINUE
          ENDIF
        ENDIF
*----
*  GEOMETRY CANNOT BE TRACKED BY JPM
*---
        IF(IROT.GT.0.AND.IZRT.GT.1) CALL XABORT(NAMSBR//
     >  ': ROD OVERLAPP -- JPM CAN NOT TRACK THIS GEOMETRY')
 122  CONTINUE
*----
*  CHECK FOR VALID CLUSTER IN JPM TRACKING
*----
      IF(IROT.GT.0) THEN
        DO  180 IAN=1,NTAN
          ILR=NRINFO(2,IAN)
          IF(ILR.GT.0) THEN
            IF(NXRI(1,IAN).NE.ILR) CALL XABORT(NAMSBR//
     >      ': ANNULUS OVERLAP PIN -- JPM CAN NOT TRACK THIS GEOMETRY')
          ENDIF
 180    CONTINUE
      ENDIF
      DEALLOCATE(JGEOM)
      IF(IPRT.GT.2) THEN
        WRITE(IOUT,6010)
        DO 600 IAN=1,NTAN
          IF((NRINFO(2,IAN).EQ.0).OR.
     >       (NRINFO(2,IAN).EQ.NXRI(1,IAN))) THEN
            WRITE(IOUT,6013) IAN,NRINFO(2,IAN),
     >                   RAN(IAN),MATANN(IAN)
          ELSE
            DO 601 IRT=1,NRT
              IF(NXRI(IRT,IAN).EQ.0) GO TO 602
              WRITE(IOUT,6013) IAN,NXRI(IRT,IAN),
     >                     RAN(IAN),MATANN(IAN)
 601        CONTINUE
          ENDIF
 602      CONTINUE
 600    CONTINUE
        IF(NSOUT.EQ.6.OR.IROT.LT.-600) THEN
          WRITE(IOUT,6030) NBAN,NRINFO(2,NBAN),RAN(NBAN),MATANN(NBAN)
        ELSE IF(NSOUT.EQ.4.OR.IROT.LT.-400) THEN
          WRITE(IOUT,6040) NBAN,NRINFO(2,NBAN),RAN(NBAN),
     >                       COTE,MATANN(NBAN)
        ENDIF
        WRITE(IOUT,6021) (4.0*VOLSUR(JSUR),JSUR=-1,-NSOUT,-1)
        WRITE(IOUT,6022) (ICODE(-MATALB(JSUR)),JSUR=-1,-NSOUT,-1)
        WRITE(IOUT,6023) (-JSW,ALBEDO(JSW),JSW=1,NMCOD)
        WRITE(IOUT,6011) (IRT,NRODS(1,IRT),NRODS(2,IRT),NRODS(3,IRT),
     >               NXRS(IRT),RODS(1,IRT),RODS(2,IRT),IRT=1,NRT)
        WRITE(IOUT,6012) ((IRT,ISR,RODR(ISR,IRT),MATROD(ISR,IRT),
     >                ISR=1,NRODS(2,IRT)),IRT=1,NRT)
      ENDIF
*----
*  FILL IN VOLSUR AND MATALB VECTORS
*----
      VOLI=0.0
      IPOS=0
      VOLFS=0.0
      DO 200 IAN=1,NTAN
        VOLROD=0.0
        VOLF=PI*RAN(IAN)*RAN(IAN)
        IF(NRINFO(2,IAN).NE.0) THEN
          IF(NRINFO(2,IAN).EQ.NXRI(1,IAN)) THEN
            VOLIS=0.0
            IRT=ABS(NRINFO(2,IAN))
            XNROD=FLOAT(NRODS(1,IRT))
            DO 202 ISV=1,NRODS(2,IRT)
              IPOS=IPOS+1
              VOLFS=PI*RODR(ISV,IRT)*RODR(ISV,IRT)*XNROD
              VOLSUR(IPOS)=VOLFS-VOLIS
              MATALB(IPOS)=MATROD(ISV,IRT)
              VOLIS=VOLFS
 202        CONTINUE
            NRODR(IRT)=IPOS
            VOLROD=VOLROD+VOLFS
          ELSE
            DO 210 IRT=1,NRT
              JRT=ABS(NXRI(IRT,IAN))
              IF(JRT.LT.1000000.AND.JRT.GT.0) THEN
                XNROD=FLOAT(NRODS(1,JRT))
                VOLIS=0.0
                ILSTR=NRODS(2,JRT)
                DO 211 ISV=1,ILSTR
                  IPOS=IPOS+1
                  VOLFS=PI*RODR(ISV,JRT)*RODR(ISV,JRT)*XNROD
                  VOLSUR(IPOS)=(VOLFS-VOLIS)
                  MATALB(IPOS)=MATROD(ISV,JRT)
                  VOLIS=VOLFS
 211            CONTINUE
                NRODR(JRT)=IPOS
                VOLROD=VOLROD+VOLFS
              ELSE IF(JRT.GT.0) THEN
*----
*  ANNULUS INTERSECT RODS
*  1) FIND X (XINT) AND Y (YINT) INTERSECTION
*     XINT=(RAN**2+RPIN**2-RODR**2)/(2*RPIN)
*     YINT=SQRT(RAN**2-XINT**2)
*  2) FIND OPENNING ANGLE FOR VOLUME LIMITED BY
*     ANNULUS (ANGA) AND ROD (ANGR)
*     ANGA=ACOS(XINT/RAN)
*     ANGR=ACOS((XINT-RPIN)/RODR)
*  3) EVALUATE VOLUME
*     VRDOUT=ANGR*RODR**2-YINT*(XINT-RPIN)
*     VANIN=ANGA*RAN**2-YINT*XINT
*     VRGOUT=VRDOUT-VANIN
*           =ANGR*RODR**2-ANGA*RAN**2+YINT*RPIN
*     VRGIN=PI*RODR*RODR-VRGOUT
*----
                JPRT=JRT/1000000
                JRT=MOD(JRT,1000000)
                ILSTR=NRODS(2,JRT)
                XNROD=FLOAT(NRODS(1,JRT))
                IF(JPRT.EQ.1) THEN
                  VANSPI=RAN(IAN)*RAN(IAN)
                  VRPSPI=RODS(1,JRT)*RODS(1,JRT)
                  VRDSPI=RODR(ILSTR,JRT)*RODR(ILSTR,JRT)
                  XINT=(VANSPI+VRPSPI-VRDSPI)/(2*RODS(1,JRT))
                  YINT=SQRT(VANSPI-XINT*XINT)
                  ANGR=ACOS((XINT-RODS(1,JRT))/RODR(ILSTR,JRT))
                  ANGA=ACOS(XINT/RAN(IAN))
                  VRGIO(1,JRT)=(ANGR*VRDSPI-ANGA*VANSPI)
     >                       +YINT*RODS(1,JRT)
                  VRGIO(2,JRT)=PI*VRDSPI-VRGIO(1,JRT)
*----
*  FIRST ANNULUS CROSSING ROD
*  COMPUTE ROD VOLUME AND ROD REGION NUMBER
*----
                  VOLIS=0.0
                  DO 212 ISV=1,ILSTR
                    IPOS=IPOS+1
                    VOLFS=PI*RODR(ISV,JRT)*RODR(ISV,JRT)*XNROD
                    VOLSUR(IPOS)=(VOLFS-VOLIS)
                    MATALB(IPOS)=MATROD(ISV,JRT)
                    VOLIS=VOLFS
 212              CONTINUE
                  NRODR(JRT)=IPOS
                  VOLROD=VOLROD+XNROD*VRGIO(2,JRT)
                ELSE IF(JPRT.EQ.2) THEN
*----
*  ROD OVERLAPP THIS ANNULUS AND PRECEEDING ANNULUS
*----
                  VANSPI=RAN(IAN)*RAN(IAN)
                  VRPSPI=RODS(1,JRT)*RODS(1,JRT)
                  VRDSPI=RODR(ILSTR,JRT)*RODR(ILSTR,JRT)
                  XINT=(VANSPI+VRPSPI-VRDSPI)/(2*RODS(1,JRT))
                  YINT=SQRT(VANSPI-XINT*XINT)
                  ANGR=ACOS((XINT-RODS(1,JRT))/RODR(ILSTR,JRT))
                  ANGA=ACOS(XINT/RAN(IAN))
                  VRGOU1=ANGR*VRDSPI-ANGA*VANSPI
     >                       +YINT*RODS(1,JRT)
                  VRGIN1=PI*VRDSPI-VRGOU1
                  VOLROD=VOLROD+XNROD*(VRGIN1-VRGIO(2,JRT))
                  VRGIO(1,JRT)=VRGOU1
                  VRGIO(2,JRT)=VRGIN1
                ELSE
*----
*  LAST ANNULUS CROSSING ROD
*----
                  VOLROD=VOLROD+XNROD*VRGIO(1,JRT)
                ENDIF
              ENDIF
 210        CONTINUE
          ENDIF
        ENDIF
        IPOS=IPOS+1
        VOLSUR(IPOS)=VOLF-VOLI-VOLROD
        MATALB(IPOS)=MATANN(IAN)
        NRINFO(1,IAN)=IPOS
        VOLI=VOLF
 200  CONTINUE
*----
*  FINAL REGION ANALYSIS FOR RECTANGLE AND HEXAGONE
*----
      IF(NSOUT.EQ.6.OR.IROT.LT.-600) THEN
        IPOS=IPOS+1
        MATALB(IPOS)=MATANN(NBAN)
        NRINFO(1,NBAN)=IPOS
        VOLF=THSQ3*RAN(NBAN)*RAN(NBAN)
        VOLSUR(IPOS)=VOLF-VOLI
      ELSE IF(NSOUT.EQ.4.OR.IROT.LT.-400) THEN
        IPOS=IPOS+1
        MATALB(IPOS)=MATANN(NBAN)
        NRINFO(1,NBAN)=IPOS
        VOLF=RAN(NBAN)*COTE
        VOLSUR(IPOS)=VOLF-VOLI
      ENDIF
*----
*  PRINT GEOMETRY INFORMATION IF REQUIRED
*----
      IF(IPRT.GT.0) THEN
        CALL LCMINF(IPGEOM,GEONAM,TEXT12,EMPTY,ILONG,LCM)
        IF(NSOUT.EQ.6.OR.IROT.LT.-600) THEN
          WRITE(IOUT,'(/31H 2-D HEXAGONAL CLUSTER GEOMETRY,
     >    21H BASED ON GEOMETRY : ,A12,1H./)') GEONAM
        ELSE IF(NSOUT.EQ.4.OR.IROT.LT.-400) THEN
          WRITE(IOUT,'(/28H 2-D SQUARE CLUSTER GEOMETRY,
     >    21H BASED ON GEOMETRY : ,A12,1H./)') GEONAM
        ELSE
          WRITE(IOUT,'(/33H 2-D CYLINDRICAL CLUSTER GEOMETRY,
     >    21H BASED ON GEOMETRY : ,A12,1H./)') GEONAM
        ENDIF
        IF (.NOT.ILK) WRITE(IOUT,'(17H INFINITE DOMAIN.)')
      ENDIF
*----
*  PRINT REGION VOLUME AND MATERIAL INFORMATION WHEN REQUIRED
*----
      IF(IPRT.GT.2) THEN
        WRITE(IOUT,6000)
        IREG=0
        DO 610 IAN=1,NTAN
          IREG=IREG+1
          IF(NRINFO(2,IAN).EQ.0) THEN
            WRITE(IOUT,6001) IAN,IREG,MATALB(IREG),VOLSUR(IREG)
          ELSE
            IF(NRINFO(2,IAN).EQ.NXRI(1,IAN)) THEN
              IRT=ABS(NRINFO(2,IAN))
              IF(IRT.LT.2000000) THEN
                LRT=MOD(IRT,1000000)
                DO 612 ISV=1,NRODS(2,LRT)
                  WRITE(IOUT,6002) ISV,IREG,
     >                   MATALB(IREG),VOLSUR(IREG)
                  IREG=IREG+1
 612            CONTINUE
              ENDIF
            ELSE
              DO 613 JRT=1,NRT
                KRT=ABS(NXRI(JRT,IAN))
                IF((KRT.LT.2000000).AND.(KRT.GE.1)) THEN
                  LRT=MOD(KRT,1000000)
                  DO 614 ISV=1,NRODS(2,LRT)
                    WRITE(IOUT,6002) ISV,IREG,
     >                MATALB(IREG),VOLSUR(IREG)
                    IREG=IREG+1
 614              CONTINUE
                ENDIF
 613          CONTINUE
            ENDIF
            WRITE(IOUT,6001) IAN,IREG,MATALB(IREG),VOLSUR(IREG)
          ENDIF
 610    CONTINUE
*----
*  LAST REGION FOR SQUARE AND HEXAGONES
*----
        IF(NSOUT.EQ.6.OR.IROT.LT.-600) THEN
          IREG=IREG+1
          WRITE(IOUT,6001) IAN,IREG,MATALB(IREG),VOLSUR(IREG)
        ELSE IF(NSOUT.EQ.4.OR.IROT.LT.-400) THEN
          IREG=IREG+1
          WRITE(IOUT,6001) IAN,IREG,MATALB(IREG),VOLSUR(IREG)
        ENDIF
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(VRGIO,RAD)
      DEALLOCATE(ISPLIT,MATROD,MATANN)
      RETURN
*----
*  GEOMETRY DESCRIPTION FORMATS
*----
 6000 FORMAT(//1X,'CLUSTER GEOMETRICAL DESCRIPTION.'/
     >1X,'ANN',2X,'ROD',2X,'REG',9X,'MATERIAL',7X,'VOLUME')
 6001 FORMAT(1X,I3,6X,I4,3X,I10,1P,5X,E15.7)
 6002 FORMAT(5X,I4,1X,I4,3X,I10,1P,5X,E15.7)
 6010 FORMAT(1X,'ANNULAR REGIONS DESCRIPTION'/
     >4X,'ANNULUS',5X,'ROD ARRAY',8X,'OUTER RADIUS',6X,'MIXTURE')
 6011 FORMAT(1X,'ROD CLUSTER DESCRIPTION'/
     >2X,'ROD ARRAY',5X,'NRODS',5X,'NSUBR',7X,'AND',7X,'ANF',8X,
     >'PITCH RADIUS',5X,'FIRST ROD ANGLE'/
     >(1X,5I10,5X,E15.7,5X,E15.7))
 6012 FORMAT(1X,'SUBROD DESCRIPTION'/
     >8X,'IRT',7X,'ISR',8X,'OUTER RADIUS',6X,'MIXTURE',1P/
     >(1X,2I10,5X,E15.7,1X,I10))
 6013 FORMAT(1P,(1X,I10,4X,I10,5X,E15.7,1X,I10))
 6021 FORMAT(1X,'OUTER SURFACE DESCRIPTION'/1P,6(5X,E15.7))
 6022 FORMAT(1X,'OUTER SURFACE ICODES     '/1P,6(5X,I10,5X))
 6023 FORMAT(1X,'GEOMETRICAL ALBEDOS      '/1P,6(2X,I3,E15.7))
 6040 FORMAT(1X,'RECTANGULAR REGION DESCRIPTION'/
     >2X,'RECTANGLE',5X,'ROD ARRAY',
     >8X,'X SIDE WIDTH',8X,'Y SIDE WIDTH',8X,'MIXTURE',1P/
     >(1X,I10,4X,I10,5X,E15.7,5X,E15.7,5X,I10))
 6030 FORMAT(1X,'HEXAGONAL REGIONS DESCRIPTION'/
     >3X,'HEXAHONE',5X,'ROD ARRAY',
     >10X,'SIDE WIDTH',8X,'MIXTURE',1P/
     >(1X,I10,4X,I10,5X,E15.7,5X,I10))
*----
*  ERROR MESSAGE FORMAT
*----
 9001 FORMAT(A6,': ROD TYPE ',I10,5X,'NOT INSIDE CLUSTER')
      END
