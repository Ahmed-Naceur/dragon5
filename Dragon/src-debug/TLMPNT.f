*DECK TLMPNT
      SUBROUTINE TLMPNT(IPMAT ,IFTRK ,IPRINT,NSKTRK,NBTR  ,NDIM  ,
     >                  NREG  ,NSUR  ,MXSUB ,MXSEG ,NANGL ,NBDR  ,
     >                  NPLOTS,IPLOT ,IPLP  ,DANGLT,DVNOR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To generate the Matlab instruction for drawing the
* external surface intersection points.
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
* NSUR    number of outer surfaces for problem.
* MXSUB   maximum number of subtracks in a line.
* MXSEG   maximum number of segments in a line.
* NANGL   number of direction for tracking.
* NBDR    number of direction for volume normalization.
* NPLOTS  number of plots.
* IPLOT   plot number being processed.
* IPLP    integer plot parameters.
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
      INTEGER          IPRINT,NSKTRK,NBTR,NDIM,NREG,NSUR,MXSUB,MXSEG,
     >                 NANGL,NBDR,NPLOTS,IPLOT
      INTEGER          IPLP(6,NPLOTS)
      DOUBLE PRECISION DANGLT(NDIM,NANGL),
     >                 DVNOR(NREG,NBDR)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='TLMPNT')
*----
*  Local variables for tracking file
*----
      INTEGER          ILINE,IDUM,NBSEG,NTLINE,ISEG,KSEG,
     >                 IPLANE,IPTA2,IPTA3,NSUB,ISUB,II
      DOUBLE PRECISION WEIGHT
*----
*  Other local variables
*----
      INTEGER          ISUR,IFACE,IDIR,ISV,IENTER
      DOUBLE PRECISION DXYZ(3),FLEN
      CHARACTER        TITLE*36
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
      WRITE(TITLE,'(A18,18X)') 'Points on surfaces'
      WRITE(IPMAT,7000) NAMSBR,TITLE
*----
* Identify points associated with each surface
*----
      DO ISUR=1,NSUR
        IF(IPRINT .GE. 10) THEN
          WRITE(IOUT,6002) ISUR
        ENDIF
        DO ILINE=1,NSKTRK
          READ(IFTRK) IDUM
        ENDDO
        IFACE=0
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
*----
*  Find line segment location
*----
          ISUB=0
          IENTER=-1
          DO ISEG=1,NBSEG
            ISV=NUMERO(ISEG)
            IF(ISV .EQ. -ISUR) THEN
              IF(IENTER .EQ. -1) THEN
                ISUB=ISUB+1
                IF(ISUB.GT.NSUB) THEN
                  WRITE(IOUT,9000) ILINE
                  WRITE(IOUT,9001)
     >            (NUMERO(KSEG),LENGTH(KSEG),KSEG=1,NBSEG)
                  CALL XABORT(NAMSBR//': Invalid tracking line')
                ENDIF
                DO IDIR=1,NDIM
                  DXYZ(IDIR)=TORIG(IDIR,ISUB)
                ENDDO
              ENDIF
              IENTER=-IENTER
              IFACE=IFACE+1
              IF(IFACE .EQ. 1) THEN
                WRITE(IPMAT,7002) ISUR
              ENDIF
              WRITE(IPMAT,7004) (DXYZ(IDIR),IDIR=1,NDIM),ILINE
            ELSE IF(ISV .GT. 0) THEN
              IF(NBDR.EQ.1) THEN
                FLEN=LENGTH(ISEG)/DVNOR(ISV,1)
              ELSE
                FLEN=LENGTH(ISEG)/DVNOR(ISV,KANGL(1)+1)
              ENDIF
              DO IDIR=1,NDIM
                DXYZ(IDIR)=DXYZ(IDIR)+
     >          DANGLT(IDIR,KANGL(1))*FLEN
              ENDDO
            ENDIF
          ENDDO
        ENDDO
*----
*  Write Matlab commands to print points
*----
        IF(IFACE .GE. 1) THEN
          WRITE(IPMAT,7003)
          IF(NDIM .EQ. 2) THEN
            WRITE(IPMAT,7010) ISUR,ACOL(MOD(ISUR,7))
          ELSE
            WRITE(IPMAT,7011) ISUR,ACOL(MOD(ISUR,7))
          ENDIF
          WRITE(IPMAT,7090)
          IF(IPLP(1,IPLOT) .GT. 0) WRITE(IPMAT,7091)
        ENDIF
        REWIND IFTRK
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
 6002 FORMAT(' Processing points for surface = ',I8)
*----
*  Matlab .m file format
*----
 7000 FORMAT('%'/'% Output from ',A6/'%'/'%hold on;'/
     >7Htitle(',A36,3H');)
 7002 FORMAT('% Points for Surface ',I10/
     >       'TLMSurfacePoints=[')
 7003 FORMAT(12X,'];')
 7004 FORMAT(3F18.10,2X,I10)
 7010 FORMAT('% Plot surface ',I10/
     >       'plot(TLMSurfacePoints(:,1),',
     >            'TLMSurfacePoints(:,2),',1H',A2,3H');)
 7011 FORMAT('% Plot surface ',I10/
     >       'plot3(TLMSurfacePoints(:,1),',
     >             'TLMSurfacePoints(:,2),',
     >             'TLMSurfacePoints(:,3),',1H',A2,3H');)
 7090 FORMAT('clear TLMSurfacePoints;')
 7091 FORMAT('pause ;')
 7093 FORMAT('hold off ;')
*----
*  Errors
*----
 9000 FORMAT(' ***** Error **** '/
     >       ' Number of track cycles exceeded for line ', I10)
 9001 FORMAT(1P,4(1X,I10,E20.10))
      END
