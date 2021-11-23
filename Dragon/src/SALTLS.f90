!
!-----------------------------------------------------------------------
!
!Purpose:
! To generate the standard tracking lines (isotropic tracking) for a
! geometry using the SALT algorithm.
!
!Copyright:
! Copyright (C) 2014 Ecole Polytechnique de Montreal.
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version
!
!Author(s):
! A. Hebert
!
!Parameters: input
! IFTEMP  pointer to a temporary tracking data structure in creation
!         mode.
! IPRINT  print level.
! IGTRK   flag to generate the tracking file. In the case where IGTRK=1,
!         the tracking is performed and used to evaluate the track
!         normalisation factor and the tracking file is generated.
!         When IGTRK=0, the tracking is still performed and used to
!         evaluate the track normalisation factor but the tracking file
!         is not generated.
! NFREG   number of regions.
! NBANGL  number of angles.
! NQUAD   number of quarter (in 2-D).
! RENO    track normalisation option. A value RENO=-1 implies
!         a direction dependent normalization of the tracks for 
!         the volume. A value reno=0, implies a global normalisation.
! NBDR    number of directions for track normalization (no normalization
!         when reno=1).
! IFMT    tracking file format:
!         IFMT=0 for short file;
!         IFMT=1 long file required by TLM: module.
! DENUSR  user defined track density.
! DANGLT  angle cosines.
! DDENWT  angular density for each angle.
!
!Parameters: output
! NBTDIR  number of tracks directions considered.
! MAXSGL  maximum number of segments in a line.
! NTLINE  total number of lines generated.
! DVNOR   ratio of analytic to tracked volume.
!
!-----------------------------------------------------------------------
!
SUBROUTINE SALTLS(IFTEMP,IPRINT,IGTRK,NFREG,NBANGL,NQUAD,RENO,NBDR,IFMT,DENUSR, &
                  DANGLT,DDENWT,NBTDIR,MAXSGL,NTLINE,DVNOR)
  USE PRECISION_AND_KINDS, ONLY : PDB,SMALL,PI,TWOPI,HALFPI,INFINITY
  USE SAL_GEOMETRY_MOD,    ONLY : GG
  USE SAL_TRACKING_TYPES,  ONLY : NMAX2,MINLEN,NNN,ITRAC2,RTRAC2,DPIECE,CNT,CNT0,NBTRAC,LGMORE, &
                                  LGMORE,IERR,DD0,EX0,EY0,DELX,AX,AY,DINIT,ANGTAB,ELMTAB,NB_MAX
  USE SAL_AUX_MOD,         ONLY : SAL231,SAL232,SAL235
  USE SAL_TRAJECTORY_MOD,  ONLY : SALTRA
  IMPLICIT         NONE
  !----
  !  subroutine arguments
  !----
  INTEGER     :: IFTEMP,IPRINT,IGTRK,NFREG,NBANGL,NQUAD,RENO,NBDR,IFMT,MAXSGL,NBTDIR,NTLINE,OK
  REAL(PDB)   :: DENUSR,DANGLT(2,NQUAD*NBANGL),DDENWT(NQUAD,NBANGL),DVNOR(NFREG,NBDR)
  !----
  !  local parameters
  !----
  INTEGER                 :: II0,IPHI,MQ,MQUAD,NPIECE,INODE,IANGL,II,IQUAD,NREG, &
                             FIRST,ICURR,JCURR,LASTI,NTSEG,IAVERR
  LOGICAL                 :: LGON
  REAL                    :: EPS0,X
  REAL(PDB)               :: WT,WR,DEL0,DEL2,DNEW,DOLD,ANGLE,DELM,KEEP_DELX, &
                             MOV_DELX,SUMM,NORM,XFACT,EPSILON_PDB,DCERR,DAVERR, &
                             DMVERR,DSVERR,TORIG(2)
  REAL(PDB), DIMENSION(2) :: THETA0
  REAL(PDB), ALLOCATABLE, DIMENSION(:)  :: VOL_MRG,VOLN,SURFN,CURRN
  REAL(PDB), ALLOCATABLE, DIMENSION(:,:) :: FACNRM ! aux for normalisation
  REAL, PARAMETER :: EPS3  = 1E-3
  INTEGER, PARAMETER :: FOUT =6
  !----
  !  recompute weights
  !----
  NBTDIR=0
  SUMM=0._PDB
  DO IANGL=1,NBANGL
     DO IQUAD=1,NQUAD
      IF(DDENWT(IQUAD,IANGL).GT.0._PDB) NBTDIR=NBTDIR+1
      SUMM=SUMM+DDENWT(IQUAD,IANGL)
    ENDDO
  ENDDO
  NORM=0.5_PDB/SUMM
  NB_MAX=1
  ALLOCATE(ANGTAB(2),ELMTAB(2))
  !----
  !  isotropic tracking loop
  !----
  !     define quadrature
  !     get limit intervals for radial quadrature
  !     minimum radial interval to contain one trajectory
  EPSILON_PDB=SQRT(EPSILON(X))
  EPS0=REAL(10.*EPSILON_PDB)
  DEL0=REAL(EPS3,PDB)/DENUSR
  DEL2=-1._PDB
  !     start a set of tracks
  !     CNT0     = address of beginning of trajectory - 1
  !     CNT      = address for trajectory data - 1
  !     NBTRAC   = number of trajectories in a record
  NBTRAC=0
  MAXSGL=0
  IPHI=0
  ALLOCATE(VOLN(GG%NB_NODE),SURFN(GG%NB_SURF2), &
           CURRN(GG%NB_SURF2),FACNRM(GG%NB_NODE,NBANGL*NQUAD),STAT=OK)
  IF(OK/=0) CALL XABORT('SALTLS: NOT ENOUGH MEMORY R')
  FACNRM(:GG%NB_NODE,:NBANGL*NQUAD)=0._PDB
  DO IANGL=1,NBANGL
     DO IQUAD=1,NQUAD
        IPHI=IPHI+1
        ! KEEP IPHI
        EX0=DANGLT(1,(IANGL-1)*NQUAD+IQUAD)
        EY0=DANGLT(2,(IANGL-1)*NQUAD+IQUAD)
        ANGLE=DACOS(EX0)
        ! get theta- and theta+ from angle (to be used in SAL235 to
        ! decide whether to include projections of tangents to arcs):
        THETA0(1)=ANGLE+HALFPI
        THETA0(2)=THETA0(1)+PI
        IF(THETA0(2)>TWOPI)THETA0(2)=THETA0(2)-TWOPI
        ! get projection of points onto axis orthogonal to tracking:
        ! only first and last points from macro perimeter
        CALL SAL235(NPIECE,THETA0,EX0,EY0,GG%IPAR,GG%RPAR,GG%PERIM_MAC2,GG%NPERIM_MAC2)
        ! integrate on each piece
        DOLD=DPIECE(1)
        DNEW=DPIECE(2)
        DELM=DNEW-DOLD
        IF(DELM>DEL0)THEN
           ! compute nber of intervals for step =< 1/denusr
           MQUAD=1+INT(DELM*DENUSR)
           DELM=DELM/MQUAD
           DO MQ=1,MQUAD
              DELX=DOLD+0.5_PDB*DELM
              KEEP_DELX=DELX
              WR=REAL(DELM)
              ! initialize entering distance
              DD0=-INFINITY
              LGON=.TRUE.
              DO WHILE (LGON)
                 CNT0=0
                 CNT=CNT0+NNN
                 IERR=0
                 MOV_DELX=0.
                 DO WHILE (IERR==0)
                    ! compute one trajectory:
                    AX=DELX*EY0
                    AY=-DELX*EX0
                    CALL SALTRA(DANGLT,GG%NPERIM_MAC2,GG%PERIM_MAC2,GG%ISURF2_ELEM,GG%IPAR, &
                                GG%RPAR,GG%PPERIM_NODE,GG%PERIM_NODE)
                    IF(IERR==0)THEN
                       ! if IERR=0,trajectory has entered into the element joint, moves
                       ! DELX -> DELX+EPSILON_PDB*N
                       MOV_DELX=MOV_DELX+EPS0
                       DELX=KEEP_DELX+MOV_DELX
                       CNT=CNT0+NNN
                    ENDIF
                 ENDDO
                 ! a trajectory has been completed: store angle order nber, total weight and wr
                 ITRAC2(CNT0+4)=IPHI
                 RTRAC2(CNT0+7)=DDENWT(IQUAD,IANGL)*WR*NORM
                 RTRAC2(CNT0+8)=WR
                 II0=CNT0+1
                 IF(IPRINT > 3) CALL SAL231(RTRAC2(II0:),ITRAC2(II0:),DELX,EX0,EY0,ANGLE)
                 DO II=1,ITRAC2(II0)
                   IF(RTRAC2(II0+II+NNN-1) <= 0.0) CALL XABORT('SALTLS: INVALID SEGMENT LENGTH')
                 ENDDO
                 ! compute volumes
                 CALL SAL232(ITRAC2(II0:),RTRAC2(II0:),FACNRM,SURFN,CURRN)
                 ! next line
                 IF(CNT+NNN+MINLEN>=NMAX2) CALL XABORT('SAL230: BUFFER OVERFLOW')
                 NBTRAC=NBTRAC+1
                 WT=RTRAC2(II0+7-1)
                 WR=RTRAC2(II0+8-1)
                 FIRST=1
                 LASTI=ITRAC2(II0)
                 IF(GG%NB_SURF2/=0)THEN
                   ICURR=ITRAC2(II0+5-1)
                   JCURR=ITRAC2(II0+6-1)
                 ELSE
                   ICURR=0
                   JCURR=0
                 ENDIF
                 NTSEG=LASTI-FIRST+3
                 MAXSGL=MAX(MAXSGL,NTSEG)
                 IF(IGTRK == 1) THEN
                   IF(IFMT == 1) THEN
                     TORIG(1)=AX+DANGLT(1,IPHI)*DINIT ; TORIG(2)=AY+DANGLT(2,IPHI)*DINIT ;
                     WRITE(IFTEMP) 1,NTSEG,WT,IPHI, &
                     -ICURR,(GG%NUM_MERGE(ITRAC2(II0+II-1)),II=FIRST+NNN,LASTI+NNN),-JCURR, &
                     0.5D0,(RTRAC2(II0+II-1),II=FIRST+NNN,LASTI+NNN),0.5D0, &
                     NBTRAC,1,MQ,1,TORIG(1),TORIG(2)
                   ELSE
                     WRITE(IFTEMP) 1,NTSEG,WT,IPHI, &
                     -ICURR,(GG%NUM_MERGE(ITRAC2(II0+II-1)),II=FIRST+NNN,LASTI+NNN),-JCURR, &
                     0.5D0,(RTRAC2(II0+II-1),II=FIRST+NNN,LASTI+NNN),0.5D0
                   ENDIF
                   IF(IPRINT>5) WRITE(FOUT,'(2X,''TRAJ# DELX '',''IERR = '',I6,3X,1P,E12.4,I5)') &
                   NBTRAC,DELX,IERR
                 ENDIF
                 LGON=LGMORE
              ENDDO
              ! end of one interval: move into beginning of next
             DOLD=DOLD+DELM
           ENDDO
           ! end of one piece
        ENDIF
        ! end of trajectories for this angle
        IF(IPRINT > 5) THEN
           WRITE(FOUT,'(''  ANGLE EX EY NTRA = '',1P,3E12.4,I8,/)') ANGLE,EX0,EY0,NBTRAC
        ENDIF
     ENDDO
  ENDDO
  NTLINE=NBTRAC
  !----
  !  Compute merged normalization factors
  !----
  IF(RENO/=1) THEN
     DO INODE=1,GG%NB_NODE
        VOLN(INODE)=0._PDB
        DO IANGL=1,NBANGL
           DO IQUAD=1,NQUAD
             VOLN(INODE)=VOLN(INODE)+2._PDB*FACNRM(INODE,IANGL+(IQUAD-1)*NBANGL)* &
                         DDENWT(IQUAD,IANGL)*NORM
           ENDDO
        ENDDO
     ENDDO
     DMVERR=0.0D0
     DSVERR=0.0D0
     DAVERR=0.0D0
     IAVERR=0
     DO INODE=1,GG%NB_NODE
        DO IPHI=1,NBANGL*NQUAD
           IF(RENO==0) THEN
              IF(ABS(VOLN(INODE))>SMALL) THEN
                 FACNRM(INODE,IPHI)=GG%VOL_NODE(INODE)/VOLN(INODE)
              ENDIF
           ELSE IF(RENO==-1) THEN
              IF(ABS(FACNRM(INODE,IPHI))>SMALL) THEN
                 FACNRM(INODE,IPHI)=GG%VOL_NODE(INODE)/FACNRM(INODE,IPHI)
              ENDIF
           ENDIF
        ENDDO
        IF(ABS(VOLN(INODE))>SMALL) THEN
          IAVERR=IAVERR+1
          DCERR=100.0D0*(1.0D0-GG%VOL_NODE(INODE)/VOLN(INODE))
          DMVERR=MAX(DMVERR,ABS(DCERR))
          DSVERR=DSVERR+DCERR*DCERR
          DAVERR=DAVERR+DCERR
        ENDIF
     ENDDO
     DSVERR=SQRT(DSVERR/DBLE(IAVERR))
     DAVERR=DAVERR/DBLE(IAVERR)
     IF(IPRINT > 0) WRITE(FOUT,6005) DSVERR,DMVERR,DAVERR
     IF(IPRINT > 5) THEN
        DO IPHI=1,NBANGL*NQUAD
           WRITE(*,*) 'iphi=',IPHI
           WRITE(*,*) 'facnrm : ',(FACNRM(INODE,IPHI),INODE=1,MIN(10,GG%NB_NODE))
        ENDDO
     ENDIF
     NREG=MAXVAL(GG%NUM_MERGE)
     ALLOCATE(VOL_MRG(NREG))
     DVNOR(:,:)=0._PDB ; VOL_MRG(:)=0._PDB
     DO INODE=1,GG%NB_NODE
       II=GG%NUM_MERGE(INODE)
       IF(II > NREG) CALL XABORT('bug')
       IF(II <= 0) CYCLE
       IF(RENO==0) THEN
           DVNOR(II,1)=DVNOR(II,1)+FACNRM(INODE,1)*GG%VOL_NODE(INODE)
       ELSE IF(RENO==-1) THEN
         DO IPHI=1,NBANGL*NQUAD
           XFACT=FACNRM(INODE,IPHI)
           DVNOR(II,IPHI+1)=DVNOR(II,IPHI+1)+XFACT*GG%VOL_NODE(INODE)
         ENDDO
       ENDIF
       VOL_MRG(II)=VOL_MRG(II)+GG%VOL_NODE(INODE)
     ENDDO
     DO IPHI=1,NBDR
       DVNOR(:,IPHI)=DVNOR(:,IPHI)/VOL_MRG(:)
     ENDDO
     DEALLOCATE(VOL_MRG)
     IF(IPRINT > 4) THEN
        DO IPHI=1,NBDR
           WRITE(*,*) 'iphi=',IPHI
           WRITE(*,*) 'dvnor : ',(DVNOR(INODE,IPHI),INODE=1,MIN(10,NFREG))
        ENDDO
     ENDIF
  ENDIF
  DEALLOCATE(FACNRM,CURRN,SURFN,VOLN)
  DEALLOCATE(ELMTAB,ANGTAB)
  RETURN
  6005 FORMAT(' SALTLS: Global RMS, maximum and average errors (%) ', &
  'on region volumes :',3(2X,F10.5))
END SUBROUTINE SALTLS
