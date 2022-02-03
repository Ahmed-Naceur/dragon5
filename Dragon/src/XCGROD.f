*DECK XCGROD
      SUBROUTINE XCGROD(NRT,MSROD,NRODS,RODS,MATROD,RODR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Check geometry and reorder rod clusters if necessary.
*
*Copyright:
* Copyright (C) 1994 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* NRT     number of rod types.
* MSROD   maximum number of subrods per rods.
*
*Parameters: input/output
* NRODS   integer description of rod of a given type:
*         NRODS(1,IRT) = number of rod;
*         NRODS(2,IRT) = number of subrods in rod;
*         NRODS(3,IRT) = first concentric region.
* RODS    real description of rod of a given type:
*         RODS(1,IRT) = rod center radius;
*         RODS(2,IRT) = angular position of first rod.
* MATROD  type of material for each subrod.
* RODR    subrod radius.
*
*----------------------------------------------------------------------
*
      INTEGER    IOUT
      REAL       PI
      PARAMETER (IOUT=6,PI=3.1415926535898)
      INTEGER    NRT,NRODS(3,NRT),MATROD(MSROD,NRT)
      REAL       RODS(2,NRT),RODR(MSROD,NRT)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IORD
*----
*  SCRATCH STORAGE ALLOCATION
*   IORD    : NEW ROD CLUSTER ORDER                       I(NRT)
*----
      ALLOCATE(IORD(NRT))
*----
*  CLASSIFY ROD CLUSTER BY INCREASING DISTANCE OF CENTER AND ANGLE
*----
      DO 100 IRT=1,NRT
        IORD(IRT)=IRT
 100  CONTINUE
      DO 110 IRT=2,NRT
        REFR=RODS(1,IRT)
        REFA=RODS(2,IRT)
        IPOS=IORD(IRT)
        DO 111 JRT=IRT-1,1,-1
          KRT=JRT
          IF(RODS(1,JRT).GT.REFR) THEN
            RODS(1,JRT+1)=RODS(1,JRT)
            RODS(2,JRT+1)=RODS(2,JRT)
            IORD(JRT+1)=IORD(JRT)
          ELSE IF(RODS(1,JRT).EQ.REFR) THEN
            IPOS=-IPOS
            GO TO 112
          ELSE
            GO TO 112
          ENDIF
 111    CONTINUE
        KRT=0
 112    CONTINUE
        RODS(1,KRT+1)=REFR
        RODS(2,KRT+1)=REFA
        IORD(KRT+1)=IPOS
        IF(IPOS.LT.0) THEN
          DO 113 JRT=KRT,1,-1
            LRT=JRT
            IF((RODS(2,JRT).GT.REFA).AND.
     >         (RODS(1,JRT).EQ.REFR)) THEN
              RODS(1,JRT+1)=RODS(1,JRT)
              RODS(2,JRT+1)=RODS(2,JRT)
              IORD(JRT+1)=IORD(JRT)
            ELSE
              GO TO 114
            ENDIF
 113      CONTINUE
          LRT=0
 114      CONTINUE
          RODS(1,LRT+1)=REFR
          RODS(2,LRT+1)=REFA
          IORD(LRT+1)=-IPOS
        ENDIF
 110  CONTINUE
*----
*  REORDER REMAINING VECTORS NRODS,MATROD,RODR
*----
      DO 140 IRT=1,NRT
        JRT=IORD(IRT)
        IF(JRT.NE.IRT) THEN
          DO 141 IX=1,3
            NNR=NRODS(IX,IRT)
            NRODS(IX,IRT)=NRODS(IX,JRT)
            NRODS(IX,JRT)=NNR
 141      CONTINUE
          DO 142 IS=1,MSROD
            MATT=MATROD(IS,IRT)
            MATROD(IS,IRT)=MATROD(IS,JRT)
            MATROD(IS,JRT)=MATT
            RROD=RODR(IS,IRT)
            RODR(IS,IRT)=RODR(IS,JRT)
            RODR(IS,JRT)=RROD
 142      CONTINUE
          DO 143 KRT=IRT+1,NRT
            IF(IORD(KRT).EQ.IRT) THEN
              IORD(KRT)=JRT
              IORD(IRT)=IRT
              GO TO 144
            ENDIF
 143      CONTINUE
 144      CONTINUE
        ENDIF
 140  CONTINUE
*----
*  FIND IF ROD OVERLAPP
*----
      DO 150 IRT=1,NRT
        NRDB=NRODS(1,IRT)
        NSBRB=NRODS(2,IRT)
        RODRB=RODR(NSBRB,IRT)
        RODRB2=RODRB*RODRB
        RDPB=RODS(1,IRT)
        XBOT=RDPB-RODRB
        DANGB=2.*PI/FLOAT(NRDB)
        ANGB=RODS(2,IRT)
*----
*  CHECK FOR ROD OVERLAPP INSIDE EACH CLUSTER
*----
        IF(NRDB.GT.1) THEN
          IF(RODRB.GT.RDPB) THEN
            WRITE(IOUT,'(1X,24HROD OVERLAP IN CLUSTER =,I10)') IRT
            CALL XABORT('XCGROD: ROD OVERLAP IN A CLUSTER')
          ELSE
            ANGMIN=2.*ASIN(RODRB/RDPB)
            IF(DANGB.LE.ANGMIN) THEN
              WRITE(IOUT,'(1X,24HROD OVERLAP IN CLUSTER =,I10)') IRT
              CALL XABORT('XCGROD: ROD OVERLAP IN A CLUSTER')
            ENDIF
          ENDIF
        ENDIF
*----
*  CHECK FOR ROD OVERLAPP BETWEEN DIFFERENT CLUSTERS
*----
        DO 151 JRT=IRT-1,1,-1
          NRDT=NRODS(1,JRT)
          NSBRT=NRODS(2,JRT)
          RODRT=RODR(NSBRT,JRT)
          RODRT2=RODRT*RODRT
          RDPT=RODS(1,JRT)
          XTOP=RDPT+RODRT
          DANGT=2.*PI/FLOAT(NRDT)
          ANGT=RODS(2,JRT)
*----
*  NO OVERLAPP
*----
          IF(XTOP.LT.XBOT) GO TO 152
*----
*  SOME OVERLAPP POSSIBLE TEST FOR INTERSECTION
*----
          ANG1=ANGB
          DO 160 IA1=1,NRDB
*----
*  FIND POSITION OF ROD (X0,Y0)
*----
            X01=RDPB*COS(ANG1)
            Y01=RDPB*SIN(ANG1)
            RRX=RODRB2-X01*X01
            RRY=RODRB2-Y01*Y01
            XY=X01*Y01
            RR1=(RRX-Y01*Y01)
            ANG2=ANGT
            DO 161 IA2=1,NRDT
              X02=RDPT*COS(ANG2)
              Y02=RDPT*SIN(ANG2)
              RR2=(RODRT2-X02*X02-Y02*Y02)
*----
*  CHECK FOR ROD INSIDE ROD
*----
              DELX=X02-X01
              DELY=Y02-Y01
              DIST=SQRT(DELX**2+DELY**2)
              IF(DIST.LT.RODRT+RODRB) THEN
                WRITE(IOUT,'(1X,25HROD OVERLAP IN CLUSTERS =,2I10)')
     >                 IRT,JRT
                CALL XABORT('XCGROD: ROD OVERLAP IN 2 CLUSTERS')
              ENDIF
*----
*  FIND IF CIRCLES
*  (X-X01)**2+(Y-Y01)**2=RODRB*2
*  (X-X02)**2+(Y-Y02)**2=RODRT*2
*  INTERSECT
*----
              IF(X02.NE.X01) THEN
                CCR=1./DELX
                BBR=-DELY*CCR
                AAR=0.5*CCR*(RR1-RR2)
                ARGSQ=AAR*(2.*X01-2.*BBR*Y01-AAR)
     >               +BBR*(BBR*RRY+2.*XY)+RRX
              ELSE
                CCR=1./DELY
                BBR=-DELX*CCR
                AAR=0.5*CCR*(RR1-RR2)
                ARGSQ=AAR*(2.*Y01-2.*BBR*X01-AAR)
     >               +BBR*(BBR*RRX+2.*XY)+RRY
              ENDIF
              IF(ARGSQ.GE.0.0) THEN
                WRITE(IOUT,'(1X,25HROD OVERLAP IN CLUSTERS =,2I10)')
     >                   IRT,JRT
                CALL XABORT('XCGROD: ROD OVERLAP IN 2 CLUSTERS')
              ENDIF
              ANG2=ANG2+DANGT
 161        CONTINUE
            ANG1=ANG1+DANGB
 160      CONTINUE
 151    CONTINUE
 152    CONTINUE
 150  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IORD)
*----
*  RETURN
*----
      RETURN
      END
