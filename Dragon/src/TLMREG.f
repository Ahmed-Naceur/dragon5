*DECK TLMREG
      SUBROUTINE TLMREG(IPMAT ,IFTRK ,IPRINT,NSKTRK,NBTR  ,NDIM  ,
     >                  NSOUT ,NREG  ,MXSUB ,MXSEG ,NANGL ,NBDR  ,
     >                  NPLOTS,IPLOT, IPLP  ,DANGLT,DVNOR ,MATALB,
     >                  LMIX  )
*
*-----------------------------------------------------------------------
*
*Purpose:
* To generate the Matlab instruction for drawing the
* lines for the region selected.
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
* NSOUT   number of outer surfaces for problem.
* MXSUB   maximum number of subtracks in a line.
* MXSEG   maximum number of segments in a line.
* NANGL   number of direction for tracking.
* NBDR    number of direction for volume normalization.
* NPLOTS  number of plots.
* IPLOT   plot number being processed.
* IPLP    integer plot parameters.
* DANGLT  track directions.
* DVNOR   track normalization factor for regional volumes.
* MATALB  surface direction and region material identification array.
* LMIX    flag set to .TRUE. to draw mixture lines.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPMAT,IFTRK
      INTEGER          IPRINT,NSKTRK,NBTR,NDIM,NSOUT,NREG,MXSUB,MXSEG,
     >                 NANGL,NBDR,NPLOTS,IPLOT
      INTEGER          IPLP(6,NPLOTS),MATALB(-NSOUT:NREG)
      DOUBLE PRECISION DANGLT(NDIM,NANGL),DVNOR(NREG,NBDR)
      LOGICAL          LMIX
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='TLMREG')
*----
*  Local variables for tracking file
*----
      INTEGER          ILINE,IDUM,NBSEG,NTLINE,ISEG,KSEG,
     >                 IPLANE,IPTA2,IPTA3,NSUB,ISUB,II
      DOUBLE PRECISION WEIGHT
*----
*  Other local variables
*----
      INTEGER          IREG,ILREG,IDIR,ISV,IENTER,ITRACE,IPM
      DOUBLE PRECISION DXYZ(3,2),FLEN
      CHARACTER        TITLE*36
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NUMERO,KANGL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: LENGTH
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: TORIG
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
      IF(IPRINT .GE. 1) WRITE(IOUT,6000) NAMSBR
*----
* Print IPMAT header
*----
      IREG=IPLP(2,IPLOT)
      WRITE(TITLE,'(A4,I4)') 'Reg=',IREG
      IF(LMIX) THEN
        WRITE(IPMAT,7000) NAMSBR,TITLE,MAXVAL(MATALB(1:NREG))+1
      ELSE
        WRITE(IPMAT,7000) NAMSBR,TITLE,NREG
      ENDIF
*----
*  Print matlab instructions for line segment
*----
      IF(IPRINT .GE. 10) WRITE(IOUT,6002) IREG
      DO ILINE=1,NSKTRK
        READ(IFTRK) IDUM
      ENDDO
      ILREG=0
*----
*  Scan over lines
*----
      DO ILINE=1,NBTR
        READ(IFTRK) NSUB,NBSEG,WEIGHT,
     >              (KANGL(II),II=1,NSUB),
     >              (NUMERO(ISEG),ISEG=1,NBSEG),
     >              (LENGTH(ISEG),ISEG=1,NBSEG),
     >              NTLINE,IPLANE,IPTA2,IPTA3,
     >              ((TORIG(IDIR,ISUB),IDIR=1,NDIM),ISUB=1,NSUB)
*----
*  Find line segment location
*----
        ISUB=0
        IENTER=-1
        ITRACE=1
        DO ISEG=1,NBSEG
          ISV=NUMERO(ISEG)
          IF(ISV .GT. 0) THEN
            FLEN=LENGTH(ISEG)/DVNOR(ISV,1)
            IF(NBDR .GT. 1) FLEN=LENGTH(ISEG)/DVNOR(ISV,KANGL(1)+1)
            DO IDIR=1,NDIM
              DXYZ(IDIR,2)=DXYZ(IDIR,1)+
     >        DANGLT(IDIR,KANGL(1))*FLEN
            ENDDO
            IF(IREG .EQ. ISV .AND. ITRACE .EQ. 1) THEN
              ILREG=ILREG+1
              IF(ILREG .EQ. 1) THEN
                WRITE(IPMAT,7002)
              ENDIF
              WRITE(IPMAT,7004)
     >        ((DXYZ(IDIR,IPM),IPM=1,2),IDIR=1,NDIM)
            ENDIF
            DO IDIR=1,NDIM
              DXYZ(IDIR,1)=DXYZ(IDIR,2)
            ENDDO
          ELSE
            IF(IENTER .EQ. -1) THEN
              ISUB=ISUB+1
              IF(ISUB.GT.NSUB) THEN
                WRITE(IOUT,9000) ILINE
                WRITE(IOUT,9001)
     >          (NUMERO(KSEG),LENGTH(KSEG),KSEG=1,NBSEG)
                CALL XABORT(NAMSBR//': Invalid tracking line')
              ENDIF
              DO IDIR=1,NDIM
                DXYZ(IDIR,1)=TORIG(IDIR,ISUB)
              ENDDO
            ENDIF
            IENTER=-IENTER
          ENDIF
        ENDDO
      ENDDO
*----
*  Write Matlab commands to trace lines
*----
      IF(ILREG .GE. 1) THEN
        WRITE(IPMAT,7003)
        IF(LMIX) THEN
          IF(NDIM .EQ. 2) THEN
            WRITE(IPMAT,7012) MATALB(IREG)+1
          ELSE
            WRITE(IPMAT,7013) MATALB(IREG)+1
          ENDIF
        ELSE
          IF(NDIM .EQ. 2) THEN
            WRITE(IPMAT,7012) IREG
          ELSE
            WRITE(IPMAT,7013) IREG
          ENDIF
        ENDIF
*----
*  Change colour for next region
*----
        WRITE(IPMAT,7090)
        IF(IREG .NE. NREG) THEN
          IF(IPLP(1,IPLOT) .GT. 0) WRITE(IPMAT,7091)
        ENDIF
      ENDIF
      REWIND IFTRK
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
 7000 FORMAT('%'/'% Output from ',A6/'%'
     >/7Htitle(',A36,3H');/9Hxcol=jet(,i5,2H);)
 7002 FORMAT('TLMIntegrationLines=[')
 7003 FORMAT(12X,'];')
 7004 FORMAT(6F18.10)
 7012 FORMAT('[m,n]=size(TLMIntegrationLines);'/
     > 'for i=1:m'/
     > '  TLMcolorset=line([TLMIntegrationLines(i,1),',
     >            'TLMIntegrationLines(i,2)],',
     >           '[TLMIntegrationLines(i,3),',
     >           'TLMIntegrationLines(i,4)]);'/
     > '  set(TLMcolorset,',8H'Color',,'xcol(',i5,',:));'/
     > 'end;')
 7013 FORMAT('[m,n]=size(TLMIntegrationLines);'/
     > 'for i=1:m'/
     > '  TLMcolorset=line([TLMIntegrationLines(i,1),',
     >            'TLMIntegrationLines(i,2)],',
     >           '[TLMIntegrationLines(i,3),',
     >            'TLMIntegrationLines(i,4)],',
     >           '[TLMIntegrationLines(i,5),',
     >            'TLMIntegrationLines(i,6)]);'/
     > '  set(TLMcolorset,',8H'Color',,'xcol(',i5,',:));'/
     > 'end;')
 7090 FORMAT('clear TLMIntegrationLines TLMcolorset ;')
 7091 FORMAT('pause ;')
*----
*  Errors
*----
 9000 FORMAT(' ***** Error **** '/
     >       ' Number of track cycles exceeded for line ', I10)
 9001 FORMAT(1P,4(1X,I10,E20.10))
      END
