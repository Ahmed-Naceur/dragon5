*DECK TLMPLP
      SUBROUTINE TLMPLP(IPMAT ,IFTRK ,IPRINT,NSKTRK,NBTR  ,NDIM  ,
     >                  NREG  ,MXSUB ,MXSEG ,NANGL ,NBDR  ,
     >                  NPLOTS,IPLOT ,IPLP  ,DPLP  ,DANGLT,DVNOR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To generate the Matlab instruction for drawing the
* intersections between the lines and an plane (3D) or a line (2D)
* normal to the line directon.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* C. Plamondon, G. Marleau
*
*Parameters: input
* IPMAT   pointer to Matlab-m file.
* IFTRK   pointer to the TRACKING file.
* IPRINT  print level.
* NSKTRK  number of records to skip on track file before tracking
*         lines can be extracted.
* NBTR    numbre of tracks.
* NDIM    number of dimensions for problem.
* NREG    number of regions for problem.
* MXSUB   maximum number of subtracks in a line.
* MXSEG   maximum number of segments in a line.
* NANGL   number of direction for tracking.
* NBDR    number of direction for volume normalization.
* NPLOTS  number of plots.
* IPLOT   plot number being processed.
* IPLP    integer plot parameters.
* DPLP    real plot parameters.
* DANGLT  track directions.
* DVNOR   track normalization factor for regional volumes.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPMAT,IFTRK
      INTEGER          IPRINT,NSKTRK,NBTR,NDIM,NREG,MXSUB,MXSEG,
     >                 NANGL,NBDR,NPLOTS,IPLOT
      INTEGER          IPLP(6,NPLOTS)
      DOUBLE PRECISION DPLP(4,NPLOTS),DANGLT(NDIM,NANGL),
     >                 DVNOR(NREG,NBDR)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='TLMPLP')
      DOUBLE PRECISION DZERO,DONE
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0)
*----
*  Local variables for tracking file
*----
      INTEGER          ILINE,IDUM,NBSEG,NTLINE,ISEG,
     >                 IPLANE,IPTA2,IPTA3,NSUB,ISUB,II
      DOUBLE PRECISION WEIGHT
*----
*  Other local variables
*----
      INTEGER          IDIR,IREG,ILREG,ISV,IDL,IPL,IFDL,ILDL
      CHARACTER        TITLE*66
      DOUBLE PRECISION DXYZ(3),FLEN
      DOUBLE PRECISION NUMER,DENOM,DINT,DINP(2),DPLPT(4)
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NUMERO,KANGL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: LENGTH
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: TORIG
*----
*  Data
*----
      CHARACTER        ACOL(0:6)*2
      DATA             ACOL /'b.','g.','r.','c.','m.','y.','k.'/
*----
*  Scratch storage allocation
*   NUMERO  region/surface identification number for segment.
*   LENGTH  segment length.
*----
      ALLOCATE(NUMERO(MXSEG),LENGTH(MXSEG))
      ALLOCATE(KANGL(MXSUB),TORIG(NDIM,MXSUB))
*----
*  Processing starts:
*  print routine opening header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
*----
* Print IPMAT header
*----
      IDL=IPLP(2,IPLOT)
      IPL=IPLP(3,IPLOT)
      IF(IDL .EQ. 0) THEN
        IFDL=1
        ILDL=NANGL
      ELSE
        IFDL=IDL
        ILDL=IDL
      ENDIF
      DO IDL=IFDL,ILDL
      DO IDIR=1,NDIM
        DPLPT(IDIR)=DANGLT(IDIR,IDL)
      ENDDO
      DO IDIR=NDIM+1,3
        DPLPT(IDIR)=DZERO
      ENDDO
      DPLPT(4)=DPLP(4,IPLOT)
      IF(NDIM .EQ. 2) THEN
        WRITE(TITLE,'(A22,1P,4(1X,A2,F8.2))')
     >  'Lines crossing plane= ',
     >  'A=',DPLPT(1),'B=',DPLPT(2),'D=',DPLPT(4)
        WRITE(IPMAT,7000) NAMSBR,TITLE
      ELSE IF(NDIM .EQ. 3) THEN
        WRITE(TITLE,'(A22,1P,4(1X,A2,F8.2))')
     >  'Lines crossing plane= ',
     >  'A=',DPLPT(1),'B=',DPLPT(2),'C=',DPLPT(3),'D=',DPLPT(4)
        WRITE(IPMAT,7000) NAMSBR,TITLE
      ENDIF
*----
*  Print matlab instructions for line segment
*----
      DO IREG=1,NREG
        IF(IPRINT .GE. 10) THEN
          WRITE(IOUT,6002) IREG
        ENDIF
        DO ILINE=1,NSKTRK
          READ(IFTRK) IDUM
        ENDDO
        ILREG=0
*----
*  Scan over lines
*----
        DO ILINE=1,NBTR
          READ(IFTRK) NSUB,NBSEG,WEIGHT,
     >                (KANGL(II),II=1,NSUB),
     >                (NUMERO(ISEG),ISEG=1,NBSEG),
     >                (LENGTH(ISEG),ISEG=1,NBSEG),
     >                NTLINE,IPLANE,IPTA2,IPTA3,
     >                ((TORIG(IDIR,ISUB),IDIR=1,NDIM),ISUB=1,NSUB)
          IF(NSUB.NE.1) CALL XABORT(NAMSBR//
     >': Cyclic tracks not permitted for option PLANP')
*----
*  Verify for valid plane and angle
*----
          IF(KANGL(1) .EQ. IDL .AND.
     >      (IPL .EQ. IPLANE .OR. IPL .EQ. 0) ) THEN
*----
*  find location of first intersection between line and geometry
*----
            ISUB=1
            DO IDIR=1,NDIM
              DXYZ(IDIR)=TORIG(IDIR,ISUB)
            ENDDO
*----
*  For surface A*X+B*Y+C*Z=D
*  and line
*  X=X0+R*DANGLT(1)
*  Y=Y0+R*DANGLT(2)
*  Y=Z0+R*DANGLT(3)
*  Intersection must satisfy:
*  A*R*DANGLT(1)+B*R*DANGLT(2)+C*R*DANGLT(3)=D
*  R=(D-A*X0-B*Y0-C*Z0)/(A*DANGLT(1)+B*DANGLT(2)+C*DANGLT(3))
*----
            NUMER=DPLPT(4)
            DENOM=DZERO
            DO IDIR=1,NDIM
              DENOM=DENOM+DPLPT(IDIR)*DANGLT(IDIR,KANGL(1))
              NUMER=NUMER-DPLPT(IDIR)*DXYZ(IDIR)
            ENDDO
            IF(DENOM .NE. DZERO) THEN
*----
*  DINT=R
*----
              DINT=NUMER/DENOM
              DINP(1)=DZERO
              DO ISEG=1,NBSEG
                ISV=NUMERO(ISEG)
                IF(ISV .GT. 0) THEN
                  FLEN=LENGTH(ISEG)/DVNOR(ISV,1)
                  IF(NBDR .GT. 1)
     >            FLEN=LENGTH(ISEG)/DVNOR(ISV,KANGL(1)+1)
                  DINP(2)=DINP(1)+FLEN
                  IF(DINP(1) .LE. DINT .AND. DINT .LE. DINP(2)) THEN
                    DO IDIR=1,NDIM
                      DXYZ(IDIR)=DXYZ(IDIR)+DANGLT(IDIR,KANGL(1))*DINT
                    ENDDO
                    IF(ISV .EQ. IREG) THEN
                      ILREG=ILREG+1
                      IF(ILREG .EQ. 1) THEN
                        WRITE(IPMAT,7002)
                      ENDIF
                      WRITE(IPMAT,7004)
     >                  (DXYZ(IDIR),IDIR=1,NDIM),ILINE
                      GO TO 100
                    ENDIF
                  ENDIF
                  DINP(1)=DINP(2)
                ENDIF
              ENDDO
 100          CONTINUE
            ENDIF
          ENDIF
        ENDDO
*----
*  Write Matlab commands to print points
*----
        IF(ILREG .GE. 1) THEN
          WRITE(IPMAT,7003)
          IF(NDIM .EQ. 2) THEN
            WRITE(IPMAT,7010) ACOL(MOD(IREG,7))
          ELSE
            WRITE(IPMAT,7011) ACOL(MOD(IREG,7))
          ENDIF
          WRITE(IPMAT,7090)
          IF(IPLP(1,IPLOT) .GT. 0) WRITE(IPMAT,7091)
        ENDIF
        REWIND IFTRK
      ENDDO
      ENDDO
      WRITE(IPMAT,7093)
      IF(IPLP(1,IPLOT) .LT. 0) WRITE(IPMAT,7091)
*----
*  Processing finished, return
*----
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(KANGL,TORIG)
      DEALLOCATE(LENGTH,NUMERO)
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6002 FORMAT(' Processing lines for region = ',I8)
*----
*  Matlab .m file format
*----
 7000 FORMAT('%'/'% Output from ',A6/'%'/
     >'%figure;'/'hold on;'/7Htitle(',A66,3H');/
     >12Hxlabel('x');/12Hylabel('y');)
 7002 FORMAT('TLMSurfacePoints=[')
 7003 FORMAT(12X,'];')
 7004 FORMAT(3F18.10,2X,I10)
 7010 FORMAT('plot(TLMSurfacePoints(:,1),',
     >            'TLMSurfacePoints(:,2),',1H',A2,3H');)
 7011 FORMAT('plot3(TLMSurfacePoints(:,1),',
     >             'TLMSurfacePoints(:,2),',
     >             'TLMSurfacePoints(:,3),',1H',A2,3H');)
 7090 FORMAT('clear TLMSurfacePoints TLMcolorset;')
 7091 FORMAT('pause ;')
 7093 FORMAT('hold off;')
      END
