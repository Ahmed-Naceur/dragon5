*DECK PSPXCG
      SUBROUTINE PSPXCG(IPRINT,ISPSP ,ICOLR ,NBAN  ,NRT   ,MSROD ,
     >                  NSURX ,NSUR  ,NVOL  ,COTE  ,
     >                  RAN   ,NRODS ,RODS  ,RODR  ,NRINFO,NRODR ,
     >                  NXRI  ,KEYMRG,COLREG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Graphics for 2-D cluster geometry.
*
*Copyright:
* Copyright (C) 1999 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPRINT  print level.
* ISPSP   psp file unit.
* ICOLR   color set used:
*         = -4 fill hsb with no-contour;
*         = -3 fill cmyk with no-contour;
*         = -2 fill rgb with no-contour;
*         = -1 fill bw with no-contour;
*         =  0 no fill contour only;
*         =  1 fill bw and contour;
*         =  2 fill rgb and contour;
*         =  3 fill cmyk and contour;
*         =  4 fill hsb and contour.
* NBAN    number of concentric regions.
* NRT     number of rod types.
* MSROD   maximum number of subrods per rod.
* NSURX   number of surfaces.
* NSUR    number of surfaces.
* NVOL    number of regions.
* COTE    Y dimension for rectangle.
* RAN     radius/lattice side of region.
* NRODS   integer description of rod type:
*         NRODS(1,IRT) = number of rod;
*         NRODS(2,IRT) = number of subrods in rod;
*         NRODS(3,IRT) = associated annulus.
* RODS    description of rod of a given type:
*         RODS(1,IRT) = rod center radius;
*         RODS(2,IRT) = angle position of one rod.
* RODR    subrod radii.
* NRINFO  annular region content:
*         NRINFO(1,IAN) = new region number;
*         NRINFO(2,IAN) = associated cluster;
*                       = 0 no cluster.
* NRODR   subrod region.
* NXRI    annular region content multi-rod.
* KEYMRG  merge index.
* COLREG  region color.
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
      INTEGER          IOUT,NPTS
      REAL             PI,DIMX,DIMY,WLINE
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NPTS=6,PI=3.1415926535897932,
     >                 DIMX=3.5,DIMY=3.5,WLINE=0.002,NAMSBR='PSPXCG')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          IPRINT,ISPSP,ICOLR,NBAN,NRT,MSROD,NSURX,
     >                 NSUR,NVOL
      INTEGER          NRODS(3,NRT),NRINFO(2,NBAN),NRODR(NRT),
     >                 NXRI(NRT,NBAN),KEYMRG(NSUR:NVOL)
      REAL             COTE,RAN(NBAN),RODS(2,NRT),
     >                 RODR(MSROD,NRT),COLREG(4,NVOL)
*----
*  LOCAL PARAMETERS
*----
      INTEGER          ICOL,ICONT,IVOL,IMRG,NTAN,IPT,IRT,
     >                 NPROD,NINRD,IROD,ISBR,IAN,NSEG,KRT,JRT
      REAL             XYPOS(2,NPTS),RADEQ,FACT,ANGD,ANGR(2),
     >                 DANGR,RPIN,RROD,XINT,ANGA,
     >                 WLFAC
      INTEGER          KFS,KFR,KSS,KSR
*----
*  INITIALIZE
*    ICOL  FOR COLOR   (NONE, BW, RGB)
*    ICONT FOR CONTOUR (WITH OR WITHOUT CONTOUR)
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      KFS=0
      KFR=0
      KSS=0
      KSR=0
      ICONT=1
      WLFAC=1.0
      ICOL=ABS(ICOLR)
      IF(ICOLR .EQ. 0) THEN
        WLFAC=2.5
      ELSE IF(ICOLR .LT. 0) THEN
        ICONT=0
      ELSE
        KFS=1
        KSR=1
      ENDIF
*----
*  LOCATE PEN AT CENTER OF CELL
*  DETERMINE DIMENSION OF GRAPH USING CELL LIMIT
*  FOR HEXAGONAL CELL PRINT HEXAGONAL REGION
*  FOR CARTESIAN CELL PRINT CARTESIAN REGION
*----
      XYPOS(1,1)=DIMX
      XYPOS(2,1)=DIMY
      CALL PSMOVE(ISPSP,XYPOS,-3)
      IF(NSURX.EQ.6) THEN
        RADEQ=RAN(NBAN)
        FACT=DIMX/RADEQ
        RADEQ=DIMX
        NTAN=NBAN-1
*----
*  POSITION OF POINTS DEFINING THE HEXAGONAL SHAPE TO FILL
*----
        ANGD=0.0
        DO 100 IPT=1,NSURX
          XYPOS(1,IPT)=COS(ANGD)*RADEQ
          XYPOS(2,IPT)=SIN(ANGD)*RADEQ
          ANGD=ANGD+PI/3.0
 100    CONTINUE
        IVOL=NRINFO(1,NBAN)
        IMRG=KEYMRG(IVOL)
*----
*  FILL IF REQUIRED
*----
        CALL PSDREG(ISPSP,NSURX,XYPOS)
        IF(ICOL. GT. 0) THEN
          CALL PSFILL(ISPSP,ICOL,COLREG(1,IVOL),KFS,KFR)
        ENDIF
*----
*  STROKE CONTOUR IF REQUIRED
*----
        IF(ICONT .EQ. 1) THEN
          CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
        ENDIF
      ELSE IF(NSURX.EQ.4) THEN
        RADEQ=0.5*MAX(RAN(NBAN),COTE)
        FACT=DIMX/RADEQ
        NTAN=NBAN-1
        XYPOS(1,1)=FACT*RAN(NBAN)/2
        XYPOS(2,1)=FACT*COTE/2
        XYPOS(1,2)=-XYPOS(1,1)
        XYPOS(2,2)=XYPOS(2,1)
        XYPOS(1,3)=XYPOS(1,2)
        XYPOS(2,3)=-XYPOS(2,2)
        XYPOS(1,4)=XYPOS(1,1)
        XYPOS(2,4)=XYPOS(2,3)
        IVOL=NRINFO(1,NBAN)
        IMRG=KEYMRG(IVOL)
*----
*  FILL IF REQUIRED
*----
        CALL PSDREG(ISPSP,NSURX,XYPOS)
        IF(ICOL. GT. 0) THEN
          CALL PSFILL(ISPSP,ICOL,COLREG(1,IVOL),KFS,KFR)
        ENDIF
*----
*  STROKE CONTOUR IF REQUIRED
*----
        IF(ICONT .EQ. 1) THEN
          CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
        ENDIF
      ELSE
        FACT=DIMX/RAN(NBAN)
        NTAN=NBAN
      ENDIF
*----
*  ANNULAR REGIONS
*----
      DO 110 IAN=NTAN,1,-1
        RADEQ=FACT*RAN(IAN)
        XYPOS(1,1)=0.0
        XYPOS(2,1)=0.0
        IVOL=NRINFO(1,IAN)
        IMRG=KEYMRG(IVOL)
*----
*  FILL IF REQUIRED
*----
        IF(ICOL. GT. 0) THEN
          CALL PSDCIR(ISPSP,XYPOS,RADEQ)
          CALL PSFILL(ISPSP,ICOL,COLREG(1,IVOL),0,0)
        ENDIF
*----
*  STROKE CONTOUR IF REQUIRED
*----
        IF(ICONT .EQ. 1) THEN
          IF(NRINFO(2,IAN) .NE. 0) THEN
            NSEG=0
            DO 111 KRT=NRINFO(2,IAN),1,-1
              JRT=NXRI(KRT,IAN)
              IF(JRT .GT. 1000000 .AND. JRT .LT. 3000000) THEN
                IRT=MOD(JRT,1000000)
                NSEG=NSEG+1
*----
*  IF ANNULAR REGION CUT BY PINS
*  DRAW ARC SEGMENT
*----
                NPROD=NRODS(1,IRT)
                NINRD=NRODS(2,IRT)
                DANGR=2.*PI/FLOAT(NPROD)
                ANGD=RODS(2,IRT)
                RROD=FACT*RODR(NINRD,IRT)
                RPIN=FACT*RODS(1,IRT)
*----
*  ANNULUS INTERSECT RODS
*  1) FIND X (XINT) AND Y (YINT) INTERSECTION
*     XINT=(RADEQ**2+RPIN**2-RROD**2)/(2*RPIN)
*     YINT=SQRT(RAN**2-XINT**2)
*  2) FIND OPENNING ANGLE FOR VOLUME LIMITED BY
*     ANNULUS (ANGA)
*     ANGA=ACOS(XINT/RADEQ)
*----
                XINT=(RADEQ**2+RPIN**2-RROD**2)
     >               /(2.0*RPIN)
                ANGA=ACOS(XINT/RADEQ)
                DO 112 IROD=1,NPROD
                  ANGR(1)=180.0*(ANGD+ANGA)/PI
                  ANGD=ANGD+DANGR
                  ANGR(2)=180.0*(ANGD-ANGA)/PI
                  CALL PSLINW(ISPSP,WLFAC*WLINE)
                  CALL PSSARC(ISPSP,XYPOS,RADEQ,ANGR)
 112            CONTINUE
              ENDIF
 111        CONTINUE
            IF(NSEG .EQ. 0) THEN
              CALL PSDCIR(ISPSP,XYPOS,RADEQ)
              CALL PSSTRK(ISPSP,WLFAC*WLINE,0,0)
            ENDIF
          ELSE
*----
*  IF ANNULAR REGION NOT CUT BY PINS
*  STROKE CIRCLES
*----
            CALL PSDCIR(ISPSP,XYPOS,RADEQ)
            CALL PSSTRK(ISPSP,WLFAC*WLINE,0,0)
          ENDIF
        ENDIF
 110  CONTINUE
*----
*  ROD CLUSTER
*----
      DO 120 IRT=NRT,1,-1
        NPROD=NRODS(1,IRT)
        NINRD=NRODS(2,IRT)
        DANGR=2.*PI/FLOAT(NPROD)
        ANGD=RODS(2,IRT)
        RPIN=FACT*RODS(1,IRT)
        DO 121 IROD=1,NPROD
          XYPOS(1,1)=RPIN*COS(ANGD)
          XYPOS(2,1)=RPIN*SIN(ANGD)
          ANGD=ANGD+DANGR
          DO 122 ISBR=NINRD,1,-1
            IVOL=NRODR(IRT)-NINRD+ISBR
            IMRG=KEYMRG(IVOL)
            RADEQ=FACT*RODR(ISBR,IRT)
*----
*  FILL IF REQUIRED
*----
            CALL PSDCIR(ISPSP,XYPOS,RADEQ)
            IF(ICOL. GT. 0) THEN
              CALL PSFILL(ISPSP,ICOL,COLREG(1,IVOL),KFS,KFR)
            ENDIF
*----
* STROKE IF REQUIRED
*----
            IF(ICONT .EQ. 1) THEN
              CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
            ENDIF
 122      CONTINUE
 121    CONTINUE
 120  CONTINUE
      XYPOS(1,1)=-DIMX
      XYPOS(2,1)=-DIMY
      CALL PSMOVE(ISPSP,XYPOS,-3)
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
