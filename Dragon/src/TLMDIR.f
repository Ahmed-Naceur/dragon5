*DECK TLMDIR
      SUBROUTINE TLMDIR(IPMAT ,IFTRK ,IPRINT,ISPEC,NSKTRK,NBTR  ,
     >                  NDIM  ,NSOUT, NREG  ,MXSUB,MXSEG ,NANGL ,
     >                  NBDR  ,NPLOTS,IPLOT ,IPLP ,DANGLT,DVNOR ,
     >                  MATALB,LMIX  )
*
*-----------------------------------------------------------------------
*
*Purpose:
* To generate the Matlab instruction for drawing the
* lines for the directions selected.
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
* ISPEC   type of boundary (=0 TISO; =1 TSPC).
* NSKTRK  number of records to skip on track file before tracking
*         lines can be extracted.
* NBTR    number of tracks.
* NDIM    number of dimensions for problem.
* NSOUT   number of surfaces for problem.
* NREG    number of regions for problem.
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
* LMIX    flag set to .true. to draw mixture lines.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPMAT,IFTRK
      INTEGER          IPRINT,ISPEC,NSKTRK,NBTR,NDIM,NSOUT,NREG,MXSUB,
     >                 MXSEG,NANGL,NBDR,NPLOTS,IPLOT
      INTEGER          IPLP(6,NPLOTS),MATALB(-NSOUT:NREG)
      DOUBLE PRECISION DANGLT(NDIM,NANGL),DVNOR(NREG,NBDR)
      LOGICAL          LMIX
*----
*  Local parameters
*----
      INTEGER          IOUT
      LOGICAL          LNEW
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='TLMDIR')
*----
*  Local variables for tracking file
*----
      INTEGER          ILINE,IDUM,NSUB,NBSEG,NTLINE,ISEG,IPLANE,IPTA2,
     >                 IPTA3
      DOUBLE PRECISION WEIGHT
*----
*  Other local variables
*----
      INTEGER          IREG,ILREG,ILSUR,IDIR,ISV,IDL,IPL,IU,IV,IPM,
     >                 ITRACE,II,ISUB
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
      ALLOCATE(NUMERO(MXSEG),LENGTH(MXSEG),KANGL(MXSUB),
     > TORIG(NDIM,MXSUB))
*----
*  Processing starts:
*  print routine opening header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 1) WRITE(IOUT,6000) NAMSBR
*----
* Print IPMAT header
*----
      IDL=IPLP(2,IPLOT)
      IPL=IPLP(3,IPLOT)
      IU=IPLP(4,IPLOT)
      IV=IPLP(5,IPLOT)
      WRITE(TITLE,'(A4,I4,A7,I1,A4,I6,A4,I6)')
     >'Dir=',IDL,' Plane=',IPL,' IU=',IU,' IV=',IV
      IF(LMIX) THEN
        WRITE(IPMAT,7000) NAMSBR,TITLE,MAXVAL(MATALB(1:NREG))+1
      ELSE
        WRITE(IPMAT,7000) NAMSBR,TITLE,NREG
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
          IF(ISPEC.EQ.0) THEN
            IF(NSUB.NE.1) CALL XABORT(NAMSBR//': NSUB.NE.1')
            DO IDIR=1,NDIM
              DXYZ(IDIR,1)=TORIG(IDIR,1)
            ENDDO
            ITRACE=0
            IF(IDL .EQ. 0 .OR. IDL .EQ. KANGL(1) ) THEN
              IF(NDIM .EQ. 2) THEN
                ITRACE=1
              ELSE
                IF(IPL .EQ. 0 .OR. IPL .EQ. IPLANE) THEN
                  IF(IU .GT. 0) THEN
                    IF(IU .EQ. IPTA2) ITRACE=1
                  ELSE IF(IV .GT. 0) THEN
                    IF(IV .EQ. IPTA3) ITRACE=1
                  ELSE
                    ITRACE=1
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
            DO ISEG=1,NBSEG
              ISV=NUMERO(ISEG)
              IF(ISV .GT. 0) THEN
                IF(NBDR.EQ.1) THEN
                  FLEN=LENGTH(ISEG)/DVNOR(ISV,1)
                ELSE
                  FLEN=LENGTH(ISEG)/DVNOR(ISV,KANGL(1)+1)
                ENDIF
                DO IDIR=1,NDIM
                  DXYZ(IDIR,2)=DXYZ(IDIR,1)+DANGLT(IDIR,KANGL(1))*FLEN
                ENDDO
                IF(IREG .EQ. ISV .AND. ITRACE .EQ. 1) THEN
                  ILREG=ILREG+1
                  IF(ILREG .EQ. 1) WRITE(IPMAT,7002)
                  WRITE(IPMAT,7004)
     >            ((DXYZ(IDIR,IPM),IPM=1,2),IDIR=1,NDIM),
     >             FLOAT(ILINE),FLOAT(ISEG),LENGTH(ISEG)
                ENDIF
                DO IDIR=1,NDIM
                  DXYZ(IDIR,1)=DXYZ(IDIR,2)
                ENDDO
              ENDIF
            ENDDO
          ELSE IF(ISPEC.EQ.1) THEN
            ISUB=0
            LNEW=.TRUE.
            DO ISEG=1,NBSEG
              ISV=NUMERO(ISEG)
              IF(ISV.GT.0) THEN
                IF(LNEW) THEN
                  ISUB=ISUB+1
                  IF(ISUB.GT.NSUB) CALL XABORT(NAMSBR//
     >                             ': NSUB overflow')
                  DO IDIR=1,NDIM
                    DXYZ(IDIR,1)=TORIG(IDIR,ISUB)
                  ENDDO
                  LNEW=.FALSE.
                ENDIF
                IF(NBDR.EQ.1) THEN
                  FLEN=LENGTH(ISEG)/DVNOR(ISV,1)
                ELSE
                  FLEN=LENGTH(ISEG)/DVNOR(ISV,KANGL(ISUB)+1)
                ENDIF
                DO IDIR=1,NDIM
                 DXYZ(IDIR,2)=DXYZ(IDIR,1)+DANGLT(IDIR,KANGL(ISUB))*FLEN
                ENDDO
                IF(IREG .EQ. ISV) THEN
                  ILREG=ILREG+1
                  IF(ILREG .EQ. 1) WRITE(IPMAT,7002)
                  IF(IDL .EQ. 0 .OR. IDL .EQ. KANGL(ISUB)) THEN
                    WRITE(IPMAT,7004)
     >              ((DXYZ(IDIR,IPM),IPM=1,2),IDIR=1,NDIM),
     >               FLOAT(ILINE),FLOAT(ISEG),LENGTH(ISEG)
                    ENDIF
                ENDIF
                DO IDIR=1,NDIM
                  DXYZ(IDIR,1)=DXYZ(IDIR,2)
                ENDDO
              ELSE
                LNEW=.TRUE.
              ENDIF
            ENDDO
            IF(ISUB.NE.NSUB) CALL XABORT(NAMSBR//': Algorithm failure')
          ENDIF
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
      ENDDO
      IF(IPLP(6,IPLOT) .EQ. 0) WRITE(IPMAT,7093)
      WRITE(IPMAT,7091)
*----
*  Print matlab instructions for surface points
*----
      IF(IPLP(6,IPLOT) .EQ. 1) THEN
        DO ILINE=1,NSKTRK
          READ(IFTRK) IDUM
        ENDDO
        ILSUR=0
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
          IF(ISPEC.EQ.0) THEN
            IF(NSUB.NE.1) CALL XABORT(NAMSBR//': NSUB.NE.1')
            DO IDIR=1,NDIM
              DXYZ(IDIR,1)=TORIG(IDIR,1)
            ENDDO
            ITRACE=0
            IF(IDL .EQ. 0 .OR. IDL .EQ. KANGL(1) ) THEN
              IF(NDIM .EQ. 2) THEN
                ITRACE=1
              ELSE
                IF(IPL .EQ. 0 .OR. IPL .EQ. IPLANE) THEN
                  IF(IU .GT. 0) THEN
                    IF(IU .EQ. IPTA2) ITRACE=1
                  ELSE IF(IV .GT. 0) THEN
                    IF(IV .EQ. IPTA3) ITRACE=1
                  ELSE
                    ITRACE=1
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
            DO ISEG=1,NBSEG
              ISV=NUMERO(ISEG)
              IF(ISV .LT. 0) THEN
                IF(ITRACE .EQ. 1) THEN
                  ILSUR=ILSUR+1
                  IF(ILSUR .EQ. 1) WRITE(IPMAT,7005)
                  WRITE(IPMAT,7004)
     >            (DXYZ(IDIR,1),IDIR=1,NDIM),
     >             FLOAT(ILINE),FLOAT(ISEG),LENGTH(ISEG)
                ENDIF
              ELSE IF(ISV .GT. 0) THEN
                IF(NBDR.EQ.1) THEN
                  FLEN=LENGTH(ISEG)/DVNOR(ISV,1)
                ELSE
                  FLEN=LENGTH(ISEG)/DVNOR(ISV,KANGL(1)+1)
                ENDIF
                DO IDIR=1,NDIM
                  DXYZ(IDIR,1)=DXYZ(IDIR,1)+DANGLT(IDIR,KANGL(1))*FLEN
                ENDDO
              ENDIF
            ENDDO
          ELSE IF(ISPEC.EQ.1) THEN
            ISUB=0
            LNEW=.TRUE.
            DO ISEG=1,NBSEG
              ISV=NUMERO(ISEG)
              IF(ISV.GT.0) THEN
                IF(LNEW) THEN
                  ISUB=ISUB+1
                  IF(ISUB.GT.NSUB) CALL XABORT(NAMSBR//
     >                             ': NSUB overflow')
                  DO IDIR=1,NDIM
                    DXYZ(IDIR,1)=TORIG(IDIR,ISUB)
                  ENDDO
                  LNEW=.FALSE.
                ENDIF
                IF(NBDR.EQ.1) THEN
                  FLEN=LENGTH(ISEG)/DVNOR(ISV,1)
                ELSE
                  FLEN=LENGTH(ISEG)/DVNOR(ISV,KANGL(ISUB)+1)
                ENDIF
                DO IDIR=1,NDIM
                 DXYZ(IDIR,1)=DXYZ(IDIR,1)+DANGLT(IDIR,KANGL(ISUB))*FLEN
                ENDDO
              ELSE
                LNEW=.TRUE.
              ENDIF
            ENDDO
            IF(ISUB.NE.NSUB) CALL XABORT(NAMSBR//': Algorithm failure')
          ENDIF
        ENDDO
*----
*  Write Matlab commands to trace lines
*----
        WRITE(IPMAT,7003)
        IF(NDIM .EQ. 2) THEN
          WRITE(IPMAT,7010)
        ELSE
          WRITE(IPMAT,7011)
        ENDIF
        WRITE(IPMAT,7092)
        WRITE(IPMAT,7093)
        WRITE(IPMAT,7091)
        REWIND IFTRK
      ENDIF
*----
*  Processing finished, return
*----
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(TORIG,KANGL,LENGTH,NUMERO)
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
     >'hold on;'/7Htitle(',A36,3H');/9Hxcol=jet(,i5,2H);)
 7002 FORMAT('TLMIntegrationLines=[')
 7003 FORMAT(12X,'];')
 7004 FORMAT(9F18.10)
 7005 FORMAT('TLMSurfacePoints=[')
 7010 FORMAT('plot(TLMSurfacePoints(:,1),',
     >            'TLMSurfacePoints(:,2),',6H'k.');)
 7011 FORMAT('plot3(TLMSurfacePoints(:,1),',
     >             'TLMSurfacePoints(:,2),',
     >             'TLMSurfacePoints(:,3),',6H'k.');)
 7012 FORMAT('[m,n]=size(TLMIntegrationLines);'/
     > 'idreg=',I5,';'/
     > 'for i=1:m'/
     > '  TLMcolorset=line([TLMIntegrationLines(i,1),',
     >            'TLMIntegrationLines(i,2)],',
     >           '[TLMIntegrationLines(i,3),',
     >           'TLMIntegrationLines(i,4)]);'/
     > '  set(TLMcolorset,',8H'Color',,'xcol(idreg,:));'/
     > 'end;')
 7013 FORMAT('[m,n]=size(TLMIntegrationLines);'/
     > 'idreg=',I5,';'/
     > 'for i=1:m'/
     > '  TLMcolorset=line([TLMIntegrationLines(i,1),',
     >            'TLMIntegrationLines(i,2)],',
     >           '[TLMIntegrationLines(i,3),',
     >            'TLMIntegrationLines(i,4)],',
     >           '[TLMIntegrationLines(i,5),',
     >            'TLMIntegrationLines(i,6)]);'/
     > '  set(TLMcolorset,',8H'Color',,'xcol(idreg,:));'/
     > 'end;')
 7090 FORMAT('clear TLMIntegrationLines TLMcolorset ;')
 7091 FORMAT('pause ;')
 7092 FORMAT('clear TLMSurfacePoints ;')
 7093 FORMAT('hold off;')
      END
