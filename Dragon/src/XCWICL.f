*DECK XCWICL
      SUBROUTINE XCWICL(  NDIM, NSURX,  NVOL,  NBAN,   NRT, MSROD,
     >                   MAROD, NANGL,  DENS, ISYMM,IFTEMP,  IPRT,
     >                  NRINFO,   RAN,  COTE, NRODS,  RODS, NRODR,
     >                    RODR, MXSEG,  NXRI,   IMS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform isotropic tracking for 2-d cluster geometry.
*
*Copyright:
* Copyright (C) 1994 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G.Marleau
*
*Parameters: input
* NDIM    dimension of problem.
* NSURX   number of initial outer surfaces.
* NVOL    total number of regions.
* NBAN    number of concentric regions.
* NRT     number of rod types.
* MSROD   maximum number of subrod per rods.
* MAROD   maximum number of rod in any cluster.
* NANGL   number of integration angles.
* DENS    minimum parallel line trak density.
* ISYMM   integration symmetry factor.
* IFTEMP  temporary tracking file unit.
* IPRT    print level.
* NRINFO  type of concentric region:
*         NRINFO(1,IAN) = new region number;
*         NRINFO(2,IAN) = associated cluster;
*                       = 0 no cluster.
* RAN     radius/lattice side of region.
* COTE    y dimension for rectangle.
* NRODS   integer description of rod type:
*         NRODS(1,IRT) = number of rod;
*         NRODS(2,IRT) = number of subrods in rod;
*         NRODS(3,IRT) = associated annulus.
* RODS    description of rod of a given type:
*         RODS(1,IRT) = rod center radius;
*         RODS(2,IRT) = angle position of one rod.
* NRODR   subrod region.
* RODR    subrod radius.
* MXSEG   current maximum track length.
* NXRI    annular region content multi-rod.
* IMS     surface merge.
*
*----------------------------------------------------------------------
*
      PARAMETER (IUNOUT=6,PI=3.1415926535897932,SQ3=1.7320508075688773)
      INTEGER    NDIM,NSURX,NVOL,NBAN,NRT,MSROD,MAROD,NANGL,
     >           ISYMM,IFTEMP,IPRT,NXRI(NRT,NBAN),NRINFO(2,NBAN),
     >           NRODS(3,NRT),NRODR(NRT),MXSEG,INDS(2),IMS(6)
      LOGICAL    LINTER
      REAL       DENS,RAN(NBAN),COTE,RODS(2,NRT),RODR(MSROD,NRT),
     >           XPOS(2)
      DOUBLE PRECISION DCSA(2),SIDE(2),TRKPOS(2,2),ROTPOS(2,2),
     >                 DRADC,WEIGHT
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NRSEG,NNSEG
      REAL, ALLOCATABLE, DIMENSION(:) :: ATOP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DENSTY,SEGLEN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: ANGLES
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: RODP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NRSEG(MXSEG),NNSEG(MXSEG))
      ALLOCATE(ANGLES(NDIM,NANGL),DENSTY(NANGL),SEGLEN(MXSEG),
     > RODP(2,MAROD,NRT),ATOP(NRT))
*
      IF(IPRT.GE.1) THEN
        WRITE(IUNOUT,'(//1X,A20)') 'ISOTROPIC TRACKING  '
      ENDIF
*----
*  DETERMINE INTEGRATION LIMITS FOR CLUSTER REGIONS
*----
      IF(NSURX.EQ.6) THEN
        RADEQ=RAN(NBAN)
        NTAN=NBAN-1
      ELSE IF(NSURX.EQ.4) THEN
        SIDE(1)=DBLE(RAN(NBAN))
        SIDE(2)=DBLE(COTE)
        RADEQ=0.5*SQRT(RAN(NBAN)*RAN(NBAN)+COTE*COTE)
        NTAN=NBAN-1
      ELSE
        RADEQ=RAN(NBAN)
        NTAN=NBAN
      ENDIF
      IF(ISYMM.GT.1) THEN
        DANGI=4.0*PI/FLOAT(NANGL*ISYMM)
      ELSE
        DANGI=2.0*PI/FLOAT(NANGL)
      ENDIF
      NPLINE=INT(RADEQ*DENS+1.0)
      NPLINE=NPLINE+MOD(NPLINE+1,2)
      DRADI=RADEQ/FLOAT(NPLINE)
      IF(IPRT.GT.0) THEN
        WRITE(IUNOUT,6010) NVOL,NSURX,NBAN,NRT
        WRITE(IUNOUT,6011)
        WRITE(IUNOUT,6012) (II,NRODS(1,II),NRODS(2,II),
     >                      NRODS(3,II),II=1,NRT)
        WRITE(IUNOUT,6000) NANGL,DENS,NPLINE,1.0/DRADI,ISYMM
      ENDIF
      ANGD=-0.5*DANGI
      RADD=-0.5*DRADI
      WEIGHT=DRADI/DBLE(NANGL)
      DO 5 IANGL=1,NANGL
        ANGXX=ANGD+DANGI*FLOAT(IANGL)
        ANGLES(1,IANGL)=COS(ANGXX)
        ANGLES(2,IANGL)=SIN(ANGXX)
        DENSTY(IANGL)=REAL(2*NANGL)
 5    CONTINUE
      WRITE(IFTEMP) ((ANGLES(II,JJ),II=1,NDIM),JJ=1,NANGL)
      WRITE(IFTEMP) (DENSTY(JJ),JJ=1,NANGL)
*----
*  NUMBER OF RODS BETWEEN ORIGIN AND ROD 1
*----
      DO 90 IRT=1,NRT
        IF(NRODS(3,IRT).GT.0) THEN
          NBROD=NRODS(2,IRT)
          DANGR=2.*PI/FLOAT(NRODS(1,IRT))
          IF(RODR(NBROD,IRT).GT.RODS(1,IRT)) THEN
            ATOP(IRT)=0.0
          ELSE
            ATOP(IRT)=(RODS(2,IRT)
     >               +ASIN(RODR(NBROD,IRT)/RODS(1,IRT)))/DANGR
          ENDIF
        ENDIF
 90   CONTINUE
*----
*  SWEEP THROUGH TRACK ANGLES
*----
      DO 100 IANG=1,NANGL
        ANGD=ANGD+DANGI
        DCSA(1)=COS(DBLE(ANGD))
        DCSA(2)=SIN(DBLE(ANGD))
*----
*  LOCALIZE RODS WITH RESPECT TO TRAKING ANGLE
*  RODP(1,IRD,IRT)= X POSITION OF CENTER
*  RODP(2,IRD,IRT)= Y POSITION OF CENTER
*----
        DO 110 IRT=1,NRT
          IF(NRODS(3,IRT).GT.0) THEN
            DANGR=2.*PI/FLOAT(NRODS(1,IRT))
*----
*  NUMBER OF RODS BETWEEN FIRST ROD AND Y=0 TRACK
*----
            ANGC=(ANGD/DANGR)-ATOP(IRT)
            IF(ANGC.GT.0.0) THEN
              IRDEP=INT(ANGC+0.9999)
            ELSE
              IRDEP=INT(ANGC)
            ENDIF
            ANGC=RODS(2,IRT)-ANGD+IRDEP*DANGR
*----
*  STORE POSITION OF NRODS+1 RODS STARTING WITH FIRST
*  ROD ABOVE OR ON Y=0 TRACK
*----
            DO 120 IRD=1,NRODS(1,IRT)
              RODP(1,IRD,IRT)=RODS(1,IRT)*COS(ANGC)
              RODP(2,IRD,IRT)=RODS(1,IRT)*SIN(ANGC)
              ANGC=ANGC+DANGR
 120        CONTINUE
          ENDIF
 110    CONTINUE
        RADC=RADD
        DO 130 IRAD=1,NPLINE
*----
*  INITIALIZE REGION POSITION VECTOR
*----
          DO 135 ISEG=1,MXSEG
            NRSEG(ISEG)=0
            NNSEG(ISEG)=0
 135      CONTINUE
          RADC=RADC+DRADI
          RADC2=RADC*RADC
          DRADC=DBLE(RADC)
          NLSEG=MXSEG
          NFSEG=0
          NRIN=0
          IF(NSURX.EQ.6) THEN
            CALL XCWHEX(ANGD,RADC,RAN(NBAN),LINTER,XPOS,INDS,IMS)
          ELSE IF(NSURX.EQ.4) THEN
            TRKPOS(1,1)=-DRADC*DCSA(2)
            TRKPOS(2,1)=DRADC*DCSA(1)
            CALL XCWREC(DCSA,SIDE,TRKPOS,LINTER,ROTPOS,INDS,IMS)
            XPOS(1)=REAL(ROTPOS(1,1))
            XPOS(2)=REAL(ROTPOS(1,2))
          ELSE
            LINTER=.FALSE.
            INDS(1)=1
            INDS(2)=1
          ENDIF
          IF(LINTER) THEN
            NRSEG(NLSEG)=NRIN
            NNSEG(NFSEG+1)=NRIN
            SEGLEN(NLSEG)=XPOS(2)
            NLSEG=NLSEG-1
            NRIN=NRINFO(1,NBAN)
            NFSEG=NFSEG+1
            NRSEG(NFSEG)=NRIN
            NNSEG(NLSEG+1)=NRIN
            SEGLEN(NFSEG)=XPOS(1)
          ENDIF
*----
*  TRACK INSIDE ANNULAR REGIONS
*----
          DO 140 IAN=NTAN,1,-1
            IF(RADC.GT.RAN(IAN)) GO TO 141
*----
*  LINE INTERSECT ANNULUS IAN
*----
            XPOS(2)=SQRT(RAN(IAN)*RAN(IAN)-RADC2)
            XPOS(1)=-XPOS(2)
            NRSEG(NLSEG)=NRIN
            NNSEG(NFSEG+1)=NRIN
            SEGLEN(NLSEG)=XPOS(2)
            NLSEG=NLSEG-1
            NRIN=NRINFO(1,IAN)
            NFSEG=NFSEG+1
            NRSEG(NFSEG)=NRIN
            NNSEG(NLSEG+1)=NRIN
            SEGLEN(NFSEG)=XPOS(1)
            IF(NRINFO(2,IAN).NE.0) THEN
*----
*  TRACK INSIDE RODS
*----
              DO 146 KRT=1,NRT
                JRT=NXRI(KRT,IAN)
                IF((JRT.GT.3000000).OR.
     >            ((JRT.GT.0).AND.(JRT.LT.1000000)) ) THEN
                  LRT=MOD(JRT,1000000)
                  CALL XCWROD(NRIN,NRODS(1,LRT),NRODR(LRT),
     >                        RODR(1,LRT),RODP(1,1,LRT),DRADC,
     >                        NFSEG,NLSEG,SEGLEN,NRSEG,NNSEG)
                ELSE IF(JRT.EQ.0) THEN
                  GO TO 147
                ENDIF
 146          CONTINUE
 147          CONTINUE
              DO 143 KRT=1,NRT
                JRT=NXRI(KRT,IAN)
                IF(JRT.LT.0) THEN
                  IRT=-JRT
                  NXTR=NRODR(IRT)
                  DO 144 IRD=NRODS(2,IRT),1,-1
                    IF(RADC.GT.RODR(IRD,IRT)) GO TO 141
*----
*  LINE INTERSECT CENTERED ROD IRD
*----
                    XPOS(2)=SQRT(RODR(IRD,IRT)*RODR(IRD,IRT)-RADC2)
                    XPOS(1)=-XPOS(2)
                    NRSEG(NLSEG)=NRIN
                    NNSEG(NFSEG+1)=NRIN
                    SEGLEN(NLSEG)=XPOS(2)
                    NLSEG=NLSEG-1
                    NRIN=NXTR
                    NXTR=NXTR-1
                    NFSEG=NFSEG+1
                    NRSEG(NFSEG)=NRIN
                    NNSEG(NLSEG+1)=NRIN
                    SEGLEN(NFSEG)=XPOS(1)
 144              CONTINUE
                  GO TO 141
                ENDIF
 143          CONTINUE
            ENDIF
 140      CONTINUE
 141      CONTINUE
*----
*  COMPRESS AND SORT TRACK VECTOR
*----
          IF(IPRT.GE.20) THEN
            WRITE(IUNOUT,6020) IANG,ANGD,IRAD,RADC
          ENDIF
          CALL XCWSRT(IPRT,MXSEG,SEGLEN,NRSEG,NNSEG,NTSEG)
          NSEG=NTSEG
          IF(IPRT.GE.20) THEN
            WRITE(IUNOUT,6002) NSEG,-INDS(1),-INDS(2)
            WRITE(IUNOUT,6021) (SEGLEN(IIJJ),NRSEG(IIJJ),IIJJ=1,NSEG+1)
          ENDIF
*----
*  CONVERT SEGMENT DIVISION TO SEGMENT LENGTH
*----
          DO 160 ISEG=1,NSEG
            SEGLEN(ISEG)=SEGLEN(ISEG+1)-SEGLEN(ISEG)
 160      CONTINUE
          IF(NSEG+2.GT.MXSEG) THEN
            WRITE(IUNOUT,6023) NSEG,MXSEG
            WRITE(IUNOUT,6021) (SEGLEN(IIJJ),NRSEG(IIJJ),IIJJ=1,NSEG)
            CALL XABORT('XCWICL: NUMBER OF SEGMENT GREATER THAN'//
     >                  ' MAXUMUM ALLOWED')
          ENDIF
          IF(NSEG.GT.0) THEN
            WRITE(IFTEMP) 1,NSEG+2,WEIGHT,IANG,
     >                    -INDS(1),(NRSEG(JSEG),JSEG=1,NSEG),-INDS(2),
     >                    1.0D0,(SEGLEN(JSEG),JSEG=1,NSEG),1.0D0
          ENDIF
          IF(IPRT.GE.30) THEN
            WRITE(IUNOUT,6022) NSEG,-INDS(1),-INDS(2)
            WRITE(IUNOUT,6021) (SEGLEN(IIJJ),NRSEG(IIJJ),IIJJ=1,NSEG)
          ENDIF
 130    CONTINUE
 100  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ATOP,RODP,SEGLEN,DENSTY,ANGLES)
      DEALLOCATE(NNSEG,NRSEG)
      RETURN
*----
*  FORMATS
*----
 6000 FORMAT(1X,'INTEGRATION PARAMETERS',/
     >       1X,'        NUMBER OF ANGLES =',I10,/
     >       1X,'  MINIMUM TRACK DENSITY  =',1P,E15.7,/
     >       1X,'NUMBER OF PARALLEL LINES =',I10,/
     >       1X,'EFFECTIVE TRACK DENSITY  =',1P,E15.7,/
     >       1X,'         SYMMETRY FACTOR =',I10)
 6002 FORMAT(' FINAL TRACK POSITION WITH NUMBER OF SEGMENTS = ',I10/
     >       ' FIRST SURFACE INTERSECTED = ',I10,5X,
     >       '  LAST SURFACE INTERSECTED = ',I10)
 6010 FORMAT(1X,'          TOTAL NUMBER OF REGIONS =',I10/
     >       1X,'       NUMBER OF INITIAL SURFACES =',I10/
     >       1X,'        NUMBER OF ANNULAR REGIONS =',I10/
     >       1X,'             NUMBER OF RODS TYPES =',I10)
 6011 FORMAT(1X,'  ROD TYPE',10X,'  NB. RODS',10X,
     >       'NB. SUBROD',10X,'IN ANNULUS')
 6012 FORMAT((1X,I10,10X,I10,10X,I10,10X,I10))
 6020 FORMAT(//1X,' TRACKING INFORMATION'/
     >       1X,' ANGD(',I5,')=',F15.7/
     >       1X,' RADC(',I5,')=',F15.7/
     >       1X,' INTERSECTION AND REGION FOLLOWING')
 6021 FORMAT(4(5X,F15.7,I10))
 6022 FORMAT(' FINAL TRACKING LENGTH WITH NUMBER OF SEGMENTS = ',I10/
     >       ' FIRST SURFACE INTERSECTED = ',I10,4X,
     >       '  LAST SURFACE INTERSECTED = ',I10)
 6023 FORMAT(1X,' NUMBER OF SEGMENTS  ',I10,5X,'ALLOWED =',I10)
      END
