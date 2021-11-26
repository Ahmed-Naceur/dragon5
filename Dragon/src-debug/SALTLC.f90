!
!-----------------------------------------------------------------------
!
!Purpose:
! To generate the cyclic tracking lines (specular tracking) for a
! geometry using the SALT algorithm.
!
!Copyright:
! Copyright (C) 2015 Ecole Polytechnique de Montreal.
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version
!
!Author(s):
! A. Hebert
!
!Parameters: input
! IFTEMP  pointer to a temporary TRACKING file in update or creation
!         mode.
! IPRINT  print level.
! IGTRK   flag to generate the tracking file. In the case where
!         IGTRK=1, the tracking is performed and
!         used to evaluate the track normalisation factor and the
!         tracking file is generated. When IGTRK=0, the tracking is
!         still performed and used to evaluate the
!         track normalisation factor but the tracking file is not
!         generated.
! NDIM    problem dimensions.
! NFREG   number of regions.
! NBANGL  number of angles.
! RENO    track normalisation option. A value RENO=-1 implies
!         a direction dependent normalization of the tracks
!         for the volume while a value RENO=0, implies
!         a global normalisation.
! NBDR    number of directions for track normalization.
! IFMT    tracking file format:
!         IFMT=0 for short file;
!         IFMT=1 long file required by TLM: module.
! DENUSR  user defined track density.
! DANGLT  angle cosines.
! DDENWT  angular density for each angle.
! NBSANG  number of subtracks for each angles.
!
!Parameters: output
! MAXSUB  maximum number of subtracks in a line.
! MAXSGL  maximum number of segments in a line.
! NTLINE  total number of lines generated.
! DVNOR   ratio of analytic to tracked volume.
!
!-----------------------------------------------------------------------
!
SUBROUTINE SALTLC(IFTEMP,IPRINT,IGTRK,NDIM,NFREG,NBANGL,RENO,NBDR,IFMT,DENUSR,DANGLT, &
                  DDENWT,NBSANG,MAXSUB,MAXSGL,NTLINE,DVNOR)
  USE PRECISION_AND_KINDS, ONLY : PDB,INFINITY,SMALL
  USE SAL_GEOMETRY_MOD,    ONLY : GG
  USE SAL_GEOMETRY_TYPES,  ONLY : ANGGEO,TYPGEO
  USE SAL_TRACKING_TYPES,  ONLY : NMAX2,MINLEN,NNN,ITRAC2,RTRAC2,DELR,CNT,CNT0,NBTRAC,IERR,DD0,EX0, &
                                  EY0,DELX,DINIT,EX,ANGTAB,ELMTAB,TORIG,DNEW,N_AXIS,NB_TOT,NB_MAX
  USE SAL_AUX_MOD,         ONLY : SAL231,SAL232,SAL237,SAL220_1
  USE SAL_TRAJECTORY_MOD,  ONLY : SALTRA
  IMPLICIT NONE
  !----
  !  Subroutine arguments
  !----
  INTEGER :: IFTEMP,IPRINT,IGTRK,NDIM,NFREG,NBANGL,MAXSUB,MAXSGL,NTLINE,RENO,NBDR,IFMT, &
             NBSANG(5,NBANGL)
  REAL(PDB) :: DENUSR,DANGLT(NDIM,4*NBANGL),DDENWT(4*NBANGL),DVNOR(NFREG,NBDR)
  !----
  INTEGER :: AXIS(2),PIECE,NPIECE,IANGL,II0,II,KEEP_NAXIS,NN,MM,IPHI,P1,P2,NREG,INODE, &
             OK,ICURR,JCURR,IJK1,ISURF,ITR,ITRS,JPHI,NA,NAOLD,NEST,NSEG,NTRACK,NTSEG,IAVERR, &
             NMAX3
  REAL(PDB) :: WR,WT,PROJTAB(6),KEEP_DELX,MOV_DELX,XFACT,ANGLE,COSSURF,COSX,DCERR,DAVERR, &
               DMVERR,DSVERR,EPS0
  REAL :: X
  REAL(PDB), ALLOCATABLE, DIMENSION(:)  :: VOL_MRG,VOLN
  REAL(PDB), ALLOCATABLE, DIMENSION(:,:) :: FACNRM ! aux for normalisation
  INTEGER, PARAMETER :: FOUT =6
  INTEGER,   ALLOCATABLE,    DIMENSION(:)   :: ITRACK_TMP
  REAL(PDB), ALLOCATABLE,    DIMENSION(:)   :: RTRACK_TMP,COSA
  INTEGER,   POINTER, DIMENSION(:) :: ITRAC3
  REAL(PDB), POINTER, DIMENSION(:) :: RTRAC3
  !----
  !  DATA STATEMENTS
  !----
  INTEGER   SURFMAT(4,4),ITRANS(4)
  DATA      SURFMAT/ 0,3,0,4, &
                     2,0,4,1, &
                     0,1,0,2, &
                     1,1,3,0 /
  DATA      ITRANS/ 4,3,2,1 / 
  SAVE      SURFMAT,ITRANS
  !----
  CALL SAL220_1(ANGGEO)
  !
  DELR=1.0D0/DENUSR
  NB_MAX=2*MAXVAL(NBSANG(1,:)+NBSANG(2,:))
  IF(TYPGEO.EQ.7) NB_MAX=2*NB_MAX
  ALLOCATE(ANGTAB(2*NB_MAX),ELMTAB(2*NB_MAX),TORIG(2,NB_MAX),STAT=OK)
  IF(OK.NE.0) CALL XABORT('SALTLC: Not enough memory ird')
  MAXSUB=0
  MAXSGL=0
  NBTRAC=0
  CNT0=0
  CNT=CNT0+NNN
  ALLOCATE(VOLN(GG%NB_NODE),FACNRM(GG%NB_NODE,2*NBANGL),STAT=OK)
  IF(OK/=0) CALL XABORT('SALTLS: NOT ENOUGH MEMORY R')
  FACNRM(:GG%NB_NODE,:2*NBANGL)=0._PDB
  ! initialize entering distance
  DD0=-INFINITY
  EPS0=REAL(10.*SQRT(EPSILON(X)))
  DO IANGL=1,NBANGL
     ! keep IANGL
     NN=NBSANG(2,IANGL) ; MM=NBSANG(1,IANGL)
     EX0=DANGLT(1,IANGL)
     EY0=DANGLT(2,IANGL)
     ! projections of geometry outline onto the two rotation or symmetry axis
     P1=GG%PPERIM_MAC2(3); P2=GG%PPERIM_MAC2(4)-1
     CALL SAL237(EX0,EY0,MM,NN,PROJTAB,AXIS)
     ! track in direction IANGL and its relative directions
     NPIECE=AXIS(1)
     ! tracking vector
     DO PIECE=1,NPIECE
       ! axis nber, position on the axis
       N_AXIS=AXIS(2)
       DELX=PROJTAB(3)+PROJTAB(5)*(PIECE-0.5)
       KEEP_NAXIS=N_AXIS
       KEEP_DELX=DELX
       ! radial weight of the track
       WR=PROJTAB(6)
       IERR=0
       MOV_DELX=0.
       DO WHILE (IERR==0)
          ! compute one trajectory:
          ! (1) compute one trajectory
          ! (2) keep entering points if we have more trajectories
          CALL SALTRA(DANGLT(:,:2*NBANGL),GG%NPERIM_MAC2,GG%PERIM_MAC2,GG%ISURF2_ELEM,GG%IPAR, &
                      GG%RPAR,GG%PPERIM_NODE,GG%PERIM_NODE,GG%IBC2_ELEM,GG%IDATA_BC2,GG%BCDATA, &
                      GG%PPERIM_MAC2,GG%DIST_AXIS)
          IF(IERR==1) DINIT=DNEW
          IF(IERR==0) THEN
             ! if ierr=0,trajectory has entered into the seam of element joint,
             ! moves delx -> delx+epsilon_pdb*n
             MOV_DELX=MOV_DELX+EPS0
             DELX=KEEP_DELX+MOV_DELX
             N_AXIS=KEEP_NAXIS
             CNT=CNT0+NNN
          ENDIF
       ENDDO
       ! if ierr=-1,tracking have not re-entering point=no trajectory
       IF(IERR/=-1) THEN
          ! a trajectory has been completed: store angle order nber, total weight and WR
          ANGLE=DACOS(EX)
          ITRAC2(CNT0+4)=IANGL
          RTRAC2(CNT0+7)=0.5_PDB*WR/DDENWT(IANGL)
          RTRAC2(CNT0+8)=WR
          II0=CNT0+1
          IF(IPRINT > 3) CALL SAL231(RTRAC2(II0:),ITRAC2(II0:),DELX,EX0,EY0,ANGLE)
          DO II=1,ITRAC2(II0)
            IF(RTRAC2(II0+II+NNN-1) <= 0.0) CALL XABORT('SALTLC: INVALID SEGMENT LENGTH')
          ENDDO
          ! compute volumes
          CALL SAL232(ITRAC2(II0:),RTRAC2(II0:),FACNRM)
          ! next line
          IF(CNT+NNN+MINLEN>=NMAX2) THEN
             NMAX3=CNT+NNN+MINLEN+1000
             ALLOCATE(ITRAC3(2*NMAX3),RTRAC3(NMAX3),STAT=OK)
             IF(OK/=0) CALL XABORT('SALTLC: NMAX2 overflow.')
             RTRAC3(:NMAX2)=RTRAC2(:NMAX2)
             ITRAC3(:2*NMAX2)=ITRAC2(:2*NMAX2)
             DEALLOCATE(RTRAC2,ITRAC2)
             RTRAC2=>RTRAC3
             ITRAC2=>ITRAC3
             NMAX2=CNT+NNN+MINLEN+1000
          ENDIF
          NBTRAC=NBTRAC+1
          !
          ! total weight and space weight
          COSSURF=RTRAC2(II0+1)
          WT=RTRAC2(II0-1+7)
          NTRACK=ITRAC2(II0)
          NB_TOT=ITRAC2(II0+1)
          IF(NB_TOT > NB_MAX) CALL XABORT('SALTLC: NB_MAX overflow.')
          !
          ! identify entering end leaving surfaces
          ICURR=0 ; JCURR=0 ;
          IJK1=NNN+NTRACK
          NTSEG=NTRACK+2*NB_TOT
          ALLOCATE(COSA(4*NBANGL))
          DO II=1,2*NBANGL
            COSA(II)=DANGLT(1,II)
          ENDDO
          DO II=1,2*NBANGL
            COSA(4*NBANGL-II+1)=COSA(II)
          ENDDO
          JPHI=ITRAC2(II0-1+IJK1+2)
          COSX=COSA(JPHI)
          IF(ABS(COSX-COSSURF) < 1E-6) THEN
            ICURR=1
          ELSEIF(ABS(SQRT(1-COSX*COSX)-COSSURF) < 1E-6) THEN
            ICURR=2
          ELSE
            WRITE(*,*) 'COSSURF,COSX,SQRT(1-COSX*COSX) :', &
                        COSSURF,COSX,SQRT(1-COSX*COSX) 
            CALL XABORT('SALTLC: problem 1')
          ENDIF
          DEALLOCATE(COSA)
          !
          MAXSGL=MAX(MAXSGL,NTSEG)
          MAXSUB=MAX(MAXSUB,NB_TOT)
          ALLOCATE(ITRACK_TMP(NTSEG),RTRACK_TMP(NTSEG))
          !
          ! loop over sub-trajectories
          NSEG=0
          NAOLD=0
          ITR=NNN
          DO ITRS=1,NB_TOT
            IJK1=IJK1+1
            NEST=ITRAC2(II0-1+IJK1)
            IJK1=IJK1+1
            JPHI=ITRAC2(II0-1+IJK1)
            IF(IPRINT > 5) WRITE(FOUT,'(I6,"*",6X,I6,5X,I6)') ITRS,NEST,JPHI
            NA=(JPHI-1)/NBANGL+1
            IF(ITRS == 1) THEN
              ITRACK_TMP(1)=-ICURR
              RTRACK_TMP(1)=0.5
              ISURF=ICURR
              IF(IPRINT > 5) WRITE(FOUT,*) '  -> surface',ISURF
            ELSE
              IF(NA == NAOLD) THEN ! specular translation
                ISURF=ITRANS(ISURF)
              ELSE
                ISURF=SURFMAT(NA,NAOLD)
              ENDIF
              IF(IPRINT > 5) THEN
                WRITE(FOUT,*)'NA=',NA,'NAOLD=',NAOLD,' ISURF=',ISURF
                IF(ISURF == 0) CALL XABORT('SALTLC: symmetry not implemented.')
                WRITE(FOUT,*) '  -> surface',ISURF
              ENDIF
            ENDIF
            IF(ISURF /= 0) THEN
               NSEG=NSEG+1
               ITRACK_TMP(NSEG)=-ISURF
               RTRACK_TMP(NSEG)=0.5
               IF(ITRS > 1) THEN
                 NSEG=NSEG+1
                 ITRACK_TMP(NSEG)=-ISURF
                 RTRACK_TMP(NSEG)=0.5
               ENDIF
            ENDIF
            DO II=1,NEST
              ITR=ITR+1
              IF(IPRINT > 5) THEN
                 WRITE(FOUT,'(5X,I6,"*",3(1P,I6,1P,E10.2,I7))') II, &
                 ITRAC2(II0-1+ITR),RTRAC2(II0-1+ITR),ITRAC2(II0-1+ITR+NMAX2)
              ENDIF
              NSEG=NSEG+1
              ITRACK_TMP(NSEG)=GG%NUM_MERGE(ITRAC2(II0-1+ITR))
              RTRACK_TMP(NSEG)=RTRAC2(II0-1+ITR)
            ENDDO
            NAOLD=NA
          ENDDO
          ! case of a geometry with specular reflective condition which is not a
          ! rectangular -> anisotropy treatment not supported
          IF(JCURR == 0) JCURR=1
          NSEG=NSEG+1
          IF(NSEG /= NTSEG) CALL XABORT('SALTLC: NTSEG inconsistency')
          ITRACK_TMP(NSEG)=-JCURR
          RTRACK_TMP(NSEG)=0.5
          IF(IPRINT > 5) THEN
            WRITE(FOUT,*) 'SALTLC: EXCELT entry with',NSEG,'segments:'
            WRITE(FOUT,*) (ITRACK_TMP(II),II=1,NSEG)
            WRITE(FOUT,*) (RTRACK_TMP(II),II=1,NSEG)
          ENDIF
          IF(IGTRK == 1) THEN
            IF(IFMT == 1) THEN
              WRITE(IFTEMP) NB_TOT,NTSEG,WT,(ANGTAB(II),II=2,2*NB_TOT,2),(ITRACK_TMP(II),II=1,NTSEG), &
              (RTRACK_TMP(II),II=1,NTSEG),NBTRAC,1,1,1,((TORIG(II0,II),II0=1,NDIM),II=1,NB_TOT)
            ELSE
              WRITE(IFTEMP) NB_TOT,NTSEG,WT,(ANGTAB(II),II=2,2*NB_TOT,2),(ITRACK_TMP(II),II=1,NTSEG), &
              (RTRACK_TMP(II),II=1,NTSEG)
            ENDIF
          ENDIF
          !
          DEALLOCATE(RTRACK_TMP,ITRACK_TMP)
          CNT0=CNT
          IF(CNT0+NNN+MINLEN >= NMAX2) THEN
            CALL XABORT('SALTLC: NMAX2 overflow(2)')
          ELSE
            CNT=CNT0+NNN
          ENDIF
       ENDIF
     ENDDO
     !  end of trajectories for this angle
     IF(IPRINT > 5) THEN
       WRITE(FOUT,'(''SALTLC:  ANGLE EX EY = '',1P,3E12.4/)') ANGLE,EX0,EY0
     ENDIF
  ENDDO
  NTLINE=NBTRAC
  !----
  !  Compute merged normalization factors
  !----
  IF(RENO/=1) THEN
     DO INODE=1,GG%NB_NODE
        VOLN(INODE)=0._PDB
        DO IANGL=1,NBANGL
           VOLN(INODE)=VOLN(INODE)+(FACNRM(INODE,IANGL)+FACNRM(INODE,NBANGL+IANGL)) &
           /DDENWT(IANGL)
        ENDDO
     ENDDO
     DMVERR=0.0D0
     DSVERR=0.0D0
     DAVERR=0.0D0
     IAVERR=0
     DO INODE=1,GG%NB_NODE
        DO IPHI=1,2*NBANGL
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
        DO IPHI=1,2*NBANGL
           WRITE(*,*) 'iphi=',IPHI
           WRITE(*,*) 'facnrm : ',(FACNRM(INODE,IPHI),INODE=1,MIN(10,GG%NB_NODE))
        ENDDO
     ENDIF
     NREG=MAXVAL(GG%NUM_MERGE)
     ALLOCATE(VOL_MRG(NREG))
     DVNOR(:,:)=0._PDB ; VOL_MRG(:)=0._PDB
     DO INODE=1,GG%NB_NODE
       II=GG%NUM_MERGE(INODE)
       IF(II <= 0) CYCLE
       IF(RENO==0) THEN
           DVNOR(II,1)=DVNOR(II,1)+FACNRM(INODE,1)*GG%VOL_NODE(INODE)
       ELSE IF(RENO==-1) THEN
         DO IPHI=1,2*NBANGL
           XFACT=FACNRM(INODE,IPHI)
           DVNOR(II,IPHI+1)=DVNOR(II,IPHI+1)+XFACT*GG%VOL_NODE(INODE)
           DVNOR(II,2*NBANGL+IPHI+1)=DVNOR(II,IPHI+1)
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
  !----
  !  Scratch storage deallocation
  !----
  DEALLOCATE(FACNRM,VOLN)
  DEALLOCATE(TORIG,ELMTAB,ANGTAB)
  RETURN
  6005 FORMAT(' SALTLC: Global RMS, maximum and average errors (%) ', &
  'on region volumes :',3(2X,F10.5))
END SUBROUTINE SALTLC
