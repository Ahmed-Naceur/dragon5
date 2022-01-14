*DECK PSPXEL
      SUBROUTINE PSPXEL(IPRINT,ISPSP ,ICOLR ,NDIM  ,NSUR  ,NVOL  ,
     >                  NTOTCL,MAXR  ,
     >                  MINDIM,MAXDIM,KEYMRG,INDEX ,REMESH,COLREG)
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
* NDIM    number of dimensions.
* NSUR    number of surfaces.
* NVOL    number of regions.
* NTOTCL  number of cylinders.
* MAXR    dimension of REMESH vector.
* MINDIM  min index values for axes.
* MAXDIM  max index values for axes.
* KEYMRG  merge index.
* INDEX   numbering of surfaces and zones.
* REMESH  meshing.
* COLREG  region color.
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
      INTEGER          IOUT,NPTS,MXDIM,NXY,NINT
      CHARACTER        NAMSBR*6
      REAL             PI,DIMX,DIMY,WLINE
      PARAMETER       (IOUT=6,NPTS=4,MXDIM=3,NXY=2,NINT=16,
     >                 PI=3.1415926535897932,
     >                 DIMX=3.5,DIMY=3.5,WLINE=0.002,NAMSBR='PSPXEL')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          IPRINT,ISPSP,ICOLR,NDIM,
     >                 NSUR,NVOL,NTOTCL,MAXR
      INTEGER          MINDIM(NTOTCL),MAXDIM(NTOTCL),
     >                 KEYMRG(NSUR:NVOL),INDEX(4,NSUR:NVOL)
      REAL             REMESH(MAXR),COLREG(4,NVOL)
*----
*  LOCAL PARAMETERS
*----
      INTEGER          ICOL,ICONT,IDIR,IVOL,IMRG,
     >                 IX,IY,IR,ICL,NSEG,IORDER(NINT)
      REAL             WLFAC,RCIRC,OFFDIR(MXDIM),XYPOS(NXY,NPTS),
     >                 FACT,CENTER(NXY),RADANG(NXY,NINT)
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
*  COMPUTE THE CIRCUMSCRIBED RADIUS
*  THE COORDINATE FOR THE TRUE CENTER OF THE CELL
*----
      RCIRC= 0.0
      DO 100 IDIR=1,NDIM
        OFFDIR(IDIR)=0.5
     >              *(REMESH(MAXDIM(IDIR))+REMESH(MINDIM(IDIR)))
        RCIRC=MAX(RCIRC,
     >            0.5*(REMESH(MAXDIM(IDIR))-REMESH(MINDIM(IDIR))))
 100  CONTINUE
*----
*  LOCATE PEN AT CENTER OF CELL
*  DETERMINE DIMENSION OF GRAPH USING CELL LIMIT
*  FOR HEXAGONAL CELL PRINT HEXAGONAL REGION
*  FOR CARTESIAN CELL PRINT CARTESIAN REGION
*----
      XYPOS(1,1)=DIMX
      XYPOS(2,1)=DIMY
      CALL PSMOVE(ISPSP,XYPOS,-3)
      FACT=DIMX/RCIRC
*----
*  SCAN ALL REGIONS AND LOCATE POSITION
*  REGION NUMBER FROM INSIDE ANNULUS
*  TO EXTERIOR CARTESIAN
*----
      DO 110 IVOL=NVOL,1,-1
        IMRG=KEYMRG(IVOL)
        IF(IMRG .NE. 0) THEN
*----
*  CARTESIAN CELL POSITION IN X AND Y
*----
          IX=INDEX(1,IVOL)
          IY=INDEX(2,IVOL)
          XYPOS(1,1)=FACT*(REMESH(IX)-OFFDIR(1))
          XYPOS(2,1)=FACT*(REMESH(IY)-OFFDIR(2))
          XYPOS(1,2)=FACT*(REMESH(IX+1)-OFFDIR(1))
          XYPOS(2,2)=XYPOS(2,1)
          XYPOS(1,3)=XYPOS(1,2)
          XYPOS(2,3)=FACT*(REMESH(IY+1)-OFFDIR(2))
          XYPOS(1,4)=XYPOS(1,1)
          XYPOS(2,4)=XYPOS(2,3)
          IF(INDEX(4,IVOL) .EQ. 0) THEN
*----
*  CARTESIAN POSITION GEOMETRY LOCATED
*  COLOR AND TRACE IT
*----
            CALL PSDREG(ISPSP,NPTS,XYPOS)
            IF(ICOL. GT. 0) THEN
              CALL PSFILL(ISPSP,ICOL,COLREG(1,IVOL),KFS,KFR)
            ENDIF
            IF(ICONT.EQ.1) THEN
              CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
            ENDIF
          ELSE
*----
*  CARTESIAN GEOMETRY CONTAINS ANNULAR SUBDIVISION
*  DETERMINE WHICH ANNULUS
*----
            DO 111 ICL=4,NTOTCL
              IR=INDEX(4,IVOL)
              IF( IR .GE. MINDIM(ICL)-1 .AND.
     >            IR .LT. MAXDIM(ICL)      ) THEN
*----
*  ANNULUS IS DETERMINED
*  LOCATE ANNULAR/CARTESIAN AND ORDER CARTESIAN POINTS
*  FOR GEOMETRY TRACING
*----
                CENTER(1)=FACT*(REMESH(MINDIM(ICL)-2)-OFFDIR(1))
                CENTER(2)=FACT*(REMESH(MINDIM(ICL)-1)-OFFDIR(2))
                RCIRC=FACT*SQRT(REMESH(IR+1))
                CALL PSPRAI(NINT,NPTS,XYPOS,CENTER,RCIRC,
     >                      NSEG,IORDER,RADANG)
*----
*  COLOR AND TRACE RESULT
*----
                CALL PSMOVE(ISPSP,CENTER,-3)
                CALL PSDRAI(ISPSP,NSEG,IORDER,CENTER,RADANG)
                IF(ICOL. GT. 0) THEN
                  CALL PSFILL(ISPSP,ICOL,COLREG(1,IVOL),KFS,KFR)
                ENDIF
                IF(ICONT.EQ.1) THEN
                  CALL PSSTRK(ISPSP,WLFAC*WLINE,KSS,KSR)
                ENDIF
                CENTER(1)=-CENTER(1)
                CENTER(2)=-CENTER(2)
                CALL PSMOVE(ISPSP,CENTER,-3)
                GO TO 115
              ENDIF
 111        CONTINUE
 115        CONTINUE
          ENDIF
        ENDIF
 110  CONTINUE
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
