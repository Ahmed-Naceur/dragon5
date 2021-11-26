!
!---------------------------------------------------------------------
!
!Purpose:
! Support module to compute a single track.
!
!Copyright:
! Copyright (C) 2001 Ecole Polytechnique de Montreal.
!
!Author(s):
! X. Warin
!
!---------------------------------------------------------------------
!
MODULE SAL_TRAJECTORY_MOD

  USE PRECISION_AND_KINDS, ONLY : PDB,TWOPI,SMALL,PI
  USE SAL_TRACKING_TYPES,  ONLY : NBER,NBINTE,LGOK,COSINE,AX,AY,EX,EY,ALPHA,F0,AT,BT, &
       LGTYPE,R,D0
  USE SAL_NUMERIC_MOD,     ONLY : SALACO

CONTAINS

  SUBROUTINE SALTRA(DANGLT,NPERIM_MAC2,PERIM_MAC2,ISURF2_ELEM,IPAR,RPAR,PPERIM, &
  PERIM,IBC2_ELEM,IDATA_BC2,BCDATA,PPERIM_MAC2,DIST_AXIS)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! computes intersection of trajectory (T): R=A+D*E with a mesh composed
    ! of nodes and elements
    !
    !                   begin
    !                     |
    !    first entrance (LGMORE=f), get entering point
    !                     |
    !    get boundary condition at entering point and treat it
    !                     |
    !    if there is a left domain,do left tracking,inverse the trajectory
    !                     |
    !    if LGGEO1=t,test if there is a left reentering point,keep it
    !                     |
    !            do basic tracking
    !                     |
    !    if LGGEO1=t,test if there is a right reentering point,keep it
    !                     |
    !                    end
    !
    !Parameters: input
    ! DANGLT       angle cosines
    ! NPERIM_MAC2  number of elements composing perimeter of domain
    ! PERIM_MAC2   elements composing perimeter of domain
    ! ISURF2_ELEM  relative 2D surf number per elem
    ! IPAR         integer geometry descriptors
    ! RPAR         REAL(PDB) geometry descriptors
    ! PPERIM       array pointer to elements in the perimeter of nodes
    ! PERIM        array of elements in perimeter of nodes
    !
    !Parameters: input (optional data for cyclic tracking)
    ! IBC2_ELEM    relative 2D boundary condition indices per element
    ! IDATA_BC2    position of data per 2D boundary condition
    ! BCDATA       table of boundary condition descriptor
    ! PPERIM_MAC2  pointer to 'perim' and 'dist_axis':
    !              PPERIM_MAC2(1): beginning of elements on axis 1
    !              PPERIM_MAC2(2): beginning of elements on axis 2
    !              PPERIM_MAC2(3): beginning of elements not on axis
    ! DIST_AXIS    distance of points on this axis to the center (0,0)
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES,  ONLY : ISPEC
    USE SAL_TRACKING_TYPES, ONLY : NNN,NMAX2,ITRAC2,ANGTAB,ELMTAB,CNT,CNT0,NB_TOT,DNEW,DINIT, &
                               NNEW,LNEW,IERR,LGMORE,DD0,NTRACK,EPS1,EX0,EY0,EX,EY,LGOK,IPART, &
                               DELX,N_AXIS
    IMPLICIT NONE
    REAL(PDB), INTENT(IN), DIMENSION(:,:)  :: DANGLT
    INTEGER,   INTENT(IN)                  :: NPERIM_MAC2
    INTEGER,   INTENT(IN), DIMENSION(:)    :: PERIM_MAC2,PPERIM,PERIM,ISURF2_ELEM
    INTEGER,   INTENT(IN), DIMENSION(:,:)  :: IPAR
    REAL(PDB), INTENT(IN), DIMENSION(:,:)  :: RPAR
    INTEGER,   INTENT(IN), DIMENSION(:), OPTIONAL :: IBC2_ELEM,IDATA_BC2,PPERIM_MAC2
    REAL(PDB), INTENT(IN), DIMENSION(:), OPTIONAL :: DIST_AXIS
    REAL(PDB), INTENT(IN), DIMENSION(:,:), OPTIONAL  :: BCDATA
    !***
    INTEGER   :: IOUT,LEN,P1,P2,IPHI
    LOGICAL   :: LGLEFT 
    REAL(PDB) :: RADIA
    !***
    !     initiate nber of sub-trajectories
    NB_TOT=0
    IF(ISPEC == 0) THEN
      !*    compute entering point for a basic trajectory
      CALL SAL240_3(PERIM_MAC2,NPERIM_MAC2,IPAR,RPAR)
    ELSE
      EX=EX0; EY=EY0
      LGOK=.FALSE.
      !     initiate all elements status to untreated
      IPART(1,:)=-1
      !     radia of the axis
      IF(PPERIM_MAC2(N_AXIS+1)==PPERIM_MAC2(N_AXIS)) THEN
        RADIA=0.
      ELSE
        RADIA=DIST_AXIS(PPERIM_MAC2(N_AXIS+1)-1)
      ENDIF
      IF(DELX<RADIA) THEN
        !        point 'DELX' is on one of the elements on the axis
        P1=PPERIM_MAC2(N_AXIS); P2=PPERIM_MAC2(N_AXIS+1)-1
        CALL SAL241_2(P2-P1+1,PERIM_MAC2(P1:P2),DIST_AXIS(P1:P2),IPAR)
        LGOK=.TRUE.
      ELSE
        WRITE(*,*) 'PPERIM_MAC2(N_AXIS+1),PPERIM_MAC2(N_AXIS) :',PPERIM_MAC2(N_AXIS+1),PPERIM_MAC2(N_AXIS)
        WRITE(*,*) 'DIST_AXIS(PPERIM_MAC2(N_AXIS+1)-1) :',DIST_AXIS(PPERIM_MAC2(N_AXIS+1)-1)
        WRITE(*,*) 'DELX :',DELX
        CALL XABORT('SAL240_1: Cant find entering point')
      ENDIF
    ENDIF
    DINIT=DNEW
    !*    treat boundary condition on entering point
    CALL SAL245(ISURF2_ELEM,IPAR,RPAR,IOUT,LNEW,NNEW,DNEW,LGLEFT)
    !*    abort if there is a left domain
    IF(LGLEFT) CALL XABORT('SALTRA: LGLEFT True')
    !*    track the basic domain
    IF(ISPEC == 0) THEN
      IPHI=ANGLE_TO_NUMBER(EX0,EY0,DANGLT)
      CALL SAL240_4(PERIM,PPERIM,IPAR,RPAR,IPHI,ISURF2_ELEM,IOUT)
    ELSE
      CALL SAL240_4_2(DANGLT,PPERIM,PERIM,IPAR,RPAR,IDATA_BC2, &
         BCDATA,PERIM_MAC2,PPERIM_MAC2,DIST_AXIS,IBC2_ELEM)
    ENDIF
    !     trajectory entered into the element joint
    IF(IERR==0) RETURN
    !     change DD0
    IF(LGMORE) DD0=DNEW+1.001*EPS1
    !*    put NB_TOT, ANGTAB to ITRAC2
    !     total length of the trajectory
    NTRACK=CNT-(CNT0+NNN)
    ITRAC2(CNT0+1)=NTRACK
    !     total nber of the sub-trajectories
    ITRAC2(CNT0+2)=NB_TOT
    !     put ANGTAB
    LEN=2*NB_TOT
    IF((CNT+LEN)>=NMAX2) THEN
       CALL XABORT('SALTRA: Buffer overflow')
    ELSE
       ITRAC2(CNT+1:CNT+LEN)=ANGTAB(1:LEN)
       ITRAC2(CNT+NMAX2+1:CNT+NMAX2+LEN)=ELMTAB(1:LEN)
       CNT=CNT+LEN
    ENDIF
    !
  END SUBROUTINE SALTRA

  SUBROUTINE SAL241_2(NPERIM,PERIM,DIST_AXIS,IPAR)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! compute an entering point on the axial element
    !
    !Parameters: input
    ! NPERIM          = number of elements on this axis
    ! PERIM           = elements on this axis in perimeter
    ! DIST_AXIS       = distance of points on this axis to the center (0,0)
    ! IPAR            = integer descriptor of elements
    !
    !---------------------------------------------------------------------
    !
    USE SAL_TRACKING_TYPES, ONLY : IPART,N_AXIS,DNEW,DELX,NNEW,LNEW,COSINE,AX, &
                                   AY,HX,HY,BX,BY,EX,EY
    INTEGER,   INTENT(IN)                  :: NPERIM
    INTEGER,   INTENT(IN), DIMENSION(:)    :: PERIM
    REAL(PDB), INTENT(IN), DIMENSION(:)    :: DIST_AXIS
    INTEGER,   INTENT(IN), DIMENSION(:,:)  :: IPAR
    INTEGER    :: I,J
    !***
    LNEW=0
    !*    compute crossed element
    DO I=1,NPERIM
       IF(DELX<=DIST_AXIS(I)) THEN
          LNEW=PERIM(I)
          EXIT
       ENDIF
    ENDDO
    IF(LNEW==0) CALL XABORT('SAL241_2: Error of distances on the axis')
    !*    get entered node
    NNEW=IPAR(2,LNEW)
    IF(NNEW<0) NNEW=IPAR(3,LNEW)
    IF(NNEW<0) CALL XABORT('SAL241_2: Error of element data')
    !*    compute DNEW at entering point
    DNEW=DELX*(EX*HX(N_AXIS)+EY*HY(N_AXIS))
    IF(N_AXIS>2) DNEW=DNEW+BX(N_AXIS)*EX+BY(N_AXIS)*EY
    !*    compute COSINE
    COSINE=ABS(HX(N_AXIS)*EY-HY(N_AXIS)*EX)
    !*    set all elements in this axis to be 'treated (0)'
    !     others are 'untreated'
    IPART(1,:)=-1
    DO I=1,NPERIM
       J=PERIM(I)
       IPART(1,J)=0
    ENDDO
    !*    initial point
    AX=BX(N_AXIS)+DELX*HX(N_AXIS)-DNEW*EX
    AY=BY(N_AXIS)+DELX*HY(N_AXIS)-DNEW*EY
    !
  END SUBROUTINE SAL241_2
  !
  SUBROUTINE SAL240_3(PERIM_MAC2,NPERIM_MAC2,IPAR,RPAR)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! compute entering point for a basic trajectory
    !
    !Parameters: input
    ! PERIM_MAC2   elements composing perimeter of domain
    ! NPERIM_MAC2  number of elements composing perimeter of domain
    ! IPAR         integer geometry descriptors
    ! RPAR         floating point geometry descriptors
    !
    !---------------------------------------------------------------------
    !
    USE SAL_TRACKING_TYPES, ONLY : IPART,EX0,EY0,EX,EY,DD0,LGOK,LGMORE,NBER, &
         DNEW,NNEW,LNEW,DINIT
    IMPLICIT NONE
    INTEGER,   INTENT(IN)                  :: NPERIM_MAC2
    INTEGER,   INTENT(IN), DIMENSION(:)    :: PERIM_MAC2
    REAL(PDB), INTENT(IN), DIMENSION(:,:)  :: RPAR
    INTEGER,   INTENT(IN), DIMENSION(:,:)  :: IPAR
    !***
    EX=EX0; EY=EY0
    !     enter trajectory with vectors a and e and initial distance dd0
    !     initiate all elements status to untreated
    IPART(1,:)=-1
    CALL SAL241(PERIM_MAC2,NPERIM_MAC2,IPAR,RPAR,DD0,-100,DNEW,NNEW,LNEW)
    DINIT=DNEW
    !     if we have sevaral entering points
    LGMORE=NBER>3
    !*    if not succeed
    IF(.NOT.LGOK) CALL XABORT('SAL240_3: Could not enter domain')
    !
  END SUBROUTINE SAL240_3
  !
  SUBROUTINE SAL245(ISURF2_ELEM,IPAR,RPAR,IOUT,LOLD,NOLD,DOLD,LGLEFT)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! boundary condition treatment at entering point
    !
    !Parameters: input
    ! ISURF2_ELEM  relative 2D surf nber per elem
    ! IPAR         integer geometry descriptors
    ! RPAR         floating point geometry descriptors
    !
    !Parameters: input/output
    ! LOLD         entering point
    ! NOLD         node index
    ! DOLD         distance
    ! LGLEFT       flag set to TRUE if there is a left domain
    !
    !Parameters: output
    ! IOUT         (to load outside surface) = 5 (left trajectory)
    !                                          6 (right trajectory)
    !
    !---------------------------------------------------------------------
    !
    USE SAL_TRACKING_TYPES, ONLY : ITRAC2,RTRAC2,CNT0,COSINE,EX,EY
    USE SAL_GEOMETRY_TYPES, ONLY : G_BC_TYPE,ISPEC
    IMPLICIT NONE
    INTEGER,   INTENT(IN), DIMENSION(:,:)  :: IPAR
    INTEGER,   INTENT(IN), DIMENSION(:)    :: ISURF2_ELEM
    REAL(PDB), INTENT(IN), DIMENSION(:,:)  :: RPAR
    INTEGER,   INTENT(INOUT)   :: LOLD,NOLD
    REAL(PDB), INTENT(INOUT)   :: DOLD
    LOGICAL, INTENT(INOUT)     :: LGLEFT
    INTEGER, INTENT(OUT)       :: IOUT
    !***
    INTEGER   :: BCIN,SURF
    !***
    !     initiate surface info
    ITRAC2(CNT0+5)=0
    ITRAC2(CNT0+6)=0
    !*    get boundary condition at entering surface
    BCIN=IPAR(2,LOLD)
    IF(NOLD==BCIN)BCIN=IPAR(3,LOLD)
    !     add entering cosine anyhow (for characteristics)
    RTRAC2(CNT0+2)=COSINE
    LGLEFT=.FALSE.
    IF(ISPEC == 1) RETURN
    IF(BCIN>=G_BC_TYPE(0).OR.BCIN==G_BC_TYPE(-1))THEN
       !*       trajectory enters through vacuum:
       !         - compute angle, trajectory-normal and store cos and sin
       !           (inverse vector E in order to pass outgoing direction)
       CALL SAL247_1(RPAR(:,LOLD),IPAR(:,LOLD),DOLD,RTRAC2(CNT0+3),RTRAC2(CNT0+5),-EX,-EY)
       !*       get index of the entering surface
       SURF=ISURF2_ELEM(LOLD)
       IF(SURF/=0) THEN
          !           store 2D surface index
          ITRAC2(CNT0+5)=SURF
       END IF
       IOUT=6
    ELSE
       CALL XABORT('SAL245: Reversed direction')
    ENDIF
    !
  END SUBROUTINE SAL245
  !
  SUBROUTINE SAL240_4(PERIM,PPERIM,IPAR,RPAR,IPHI,ISURF2_ELEM,IOUT)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! track a trajectory until leaving the domain
    !         begin: having (nnew,dnew,lnew)
    !               |
    !     (1)write horizontal angle information
    !               |
    !     (2)compute successive crossed nodes
    !               |
    !     (3)boundary condition:
    !        reentering: get reentering point,go to (1)
    !        go out: compute surface number,end
    !
    !Parameters: input
    ! PERIM        array of elements in perimeter of nodes
    ! PPERIM       array pointer to elements in the perimeter of nodes
    ! IPAR         integer geometry descriptors
    ! RPAR         floating point geometry descriptors
    ! IPHI         angular index of the track
    ! ISURF2_ELEM  relative 2D surf number per elem
    ! IOUT         (to load outside surface) = 5 (left trajectory)
    !
    !---------------------------------------------------------------------
    !
    USE SAL_TRACKING_TYPES, ONLY : NNN,NMAX2,ITRAC2,RTRAC2,ANGTAB,ELMTAB,PRTIND,CNT, &
         CNT0,EPS1,EX,EY,LGOK,NB_TOT,NB_MAX,DNEW,NNEW,LNEW,IERR
    USE SAL_GEOMETRY_TYPES, ONLY : G_BC_TYPE
    IMPLICIT NONE
    ! IN VARIABLE
    INTEGER,   INTENT(IN)                  :: IPHI
    INTEGER,   INTENT(IN), DIMENSION(:)    :: PPERIM,PERIM
    INTEGER,   INTENT(IN), DIMENSION(:,:)  :: IPAR
    REAL(PDB), INTENT(IN), DIMENSION(:,:)  :: RPAR
    INTEGER,   INTENT(IN), DIMENSION(:)    :: ISURF2_ELEM
    INTEGER,   INTENT(IN)                  :: IOUT
    !***
    INTEGER    :: NOLD,LOLD,P1,P2,CNT1,SURF
    REAL(PDB)  :: DOLD,LENGTH
    LOGICAL    :: LGON
    INTEGER, PARAMETER :: FOUT =6
    !***
    LGON = .TRUE.
    !     initiate counter
    CNT1=CNT
    EXTERIOR : DO WHILE(LGON)
       NB_TOT=NB_TOT+1
       IF(NB_TOT > NB_MAX) CALL XABORT('SAL240_4: NB_TOT overflow')
       !        keep horizontal phi nber in angtab
       IF(NB_TOT>1) CALL XABORT('SAL240_4: Angtab overflow')
       ANGTAB(2*NB_TOT)=IPHI
       ELMTAB(2*NB_TOT-1)=LNEW ; ELMTAB(2*NB_TOT)=0 ;
       !
       !*       track a sub-trajectory         
       INTERIOR: DO WHILE(NNEW>0)
          !           UPDATE DATA TO COMPUTE NEXT NODE:
          DOLD=DNEW
          LOLD=LNEW
          NOLD=NNEW
          !           crossing NODE NOLD
          !           input: trajectory (T):R=A+D*E => A = (AX,AY), E = (EX,EY)
          !           DOLD = D at last intersection
          !           NOLD = NODE just entered
          P1=PPERIM(NOLD); P2=PPERIM(NOLD+1)-1
          CALL SAL241(PERIM(P1:P2),P2-P1+1,IPAR,RPAR,DOLD,NOLD,DNEW,NNEW,LNEW)
          !           at return from SAL241:
          !           DNEW        = D at point exiting node
          !           COSINE      = cosine of trajectory with exiting normal
          !           NNEW        = new node entered
          !           LNEW        = element crossed when exiting node
          !           NBER        = number of intersections with perimeter
          !           LGOK        = .TRUE. if trajectory exits the node
          IF(.NOT.LGOK) THEN
             IERR=0
             IF(PRTIND>0)WRITE(FOUT,'("SAL240_4 ==> couldnt exit node ",I5)') NOLD
             RETURN
          ENDIF
          !           store data
          LENGTH=DNEW-DOLD
          !           store new length in track arrays
          IF(LENGTH.LT.EPS1) CYCLE
          CNT=CNT+1
          IF(CNT>=NMAX2) CALL XABORT('SAL240_4: NMAX2 overflow')
          RTRAC2(CNT)=LENGTH
          ITRAC2(CNT+NMAX2)=LNEW
          ITRAC2(CNT)=NOLD
       END DO INTERIOR
       !
       !*       exiting motif and analyzing bc condition
       LGON =(NNEW<=G_BC_TYPE(1)).AND.(NNEW>=G_BC_TYPE(5))
       !        STORE NBER OF REGIONS TO ANGTAB
       ANGTAB(2*NB_TOT-1)=CNT-CNT1
       CNT1=CNT
       IF(LGON) THEN
          CALL XABORT('SAL240_4: Lgon is true')
       ELSE
          !           compute leaving surface
          IF(NNEW<=0)THEN
             !        we got a surface: end of trajectory
             !        vacuum bd condition. put a marker (surface number) to allow
             !        psi and pss computation
             !        store exiting surface, cosphi and sinphi
             !        compute angle trajectory-normal
             CALL SAL247_1(RPAR(:,LNEW),IPAR(:,LNEW),DNEW, &
                  RTRAC2(CNT0+IOUT-2),RTRAC2(CNT0+IOUT),EX,EY)
             !        store 2D surface-cone nber
             SURF=ISURF2_ELEM(LNEW)
             IF(SURF/=0) ITRAC2(CNT0+IOUT)=SURF
          ENDIF
       ENDIF
       !
    ENDDO EXTERIOR
    !     set success flag
    IERR=1
    !
  END SUBROUTINE SAL240_4
  !
  SUBROUTINE SAL240_4_2(DANGLT,PPERIM,PERIM,IPAR,RPAR,IDATA_BC2,BCDATA, &
       PERIM_MAC2,PPERIM_MAC2,DIST_AXIS,IBC2_ELEM)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! track a cyclic trajectory until the track length equal to the predefined
    ! length of this trajectory
    !      begin: having (nnew,dnew,lnew)
    !               |
    !     (1)write horizontal angle information
    !               |
    !     (2)compute successive crossed nodes
    !               |
    !     (3)if abs(total_length*length_inv_cycl-1)>eps1 :
    !        yes:continue, get reentering point,go to (1)
    !        no: end of tracking
    !
    !Parameters: input
    ! DANGLT       angle cosines
    ! PERIM        array of elements in perimeter of nodes
    ! PPERIM       array pointer to elements in the perimeter of nodes
    ! IPAR         integer geometry descriptors
    ! RPAR         real geometry descriptors
    ! IDATA_BC2    position of bc data per 2D boundary conditions
    ! PERIM_MAC2   elements composing perimeter of domain
    ! PPERIM_MAC2  pointer to 'PERIM' and 'DIST_AXIS':
    !                  PPERIM_MAC2(1):beginning of elements on axis 1
    !                  PPERIM_MAC2(2):beginning of elements on axis 2
    !                  PPERIM_MAC2(3):beginning of elements not on axis
    ! DIST_AXIS    distance of points on this axis to the center (0,0)
    ! BCDATA       table of bc descriptor
    ! IBC2_ELEM    relative 2D bc nber per elem
    !
    !Parameters: output
    ! NB_TOT       total number of sub-trajectories
    ! ANGTAB       table of {N_K,ANGLE_K} (1<K<NB_TOT)
    !              N_K: nber of regions in kth sub-trajectory
    !              ANGLE_K: angle nber of kth sub-trajectory
    !
    !---------------------------------------------------------------------
    !
    USE SAL_TRACKING_TYPES, ONLY : NNN,NMAX2,ITRAC2,RTRAC2,ANGTAB,ELMTAB,N_AXIS,CNT, &
         EX,EY,AX,AY,LGOK,LENGTH_INV_CYCL,NB_TOT,NB_MAX,DNEW,NNEW,LNEW,IERR,EPS1,TORIG
    USE SAL_GEOMETRY_TYPES,      ONLY : G_BC_TYPE
    IMPLICIT NONE
    ! in variable
    !************
    REAL(PDB), INTENT(IN), DIMENSION(:,:)  :: DANGLT
    INTEGER,   INTENT(IN), DIMENSION(:)    :: PPERIM,PERIM
    INTEGER,   INTENT(IN), DIMENSION(:,:)  :: IPAR
    REAL(PDB), INTENT(IN), DIMENSION(:,:)  :: RPAR
    REAL(PDB), INTENT(IN), DIMENSION(:,:)  :: BCDATA
    INTEGER,   INTENT(IN), DIMENSION(:)    :: PPERIM_MAC2,IBC2_ELEM
    INTEGER,   INTENT(IN), DIMENSION(:)    :: PERIM_MAC2
    REAL(PDB), INTENT(IN), DIMENSION(:),OPTIONAL :: DIST_AXIS
    INTEGER,   INTENT(IN), DIMENSION(:),OPTIONAL :: IDATA_BC2
    ! local variable
    INTEGER    :: NOLD,LOLD,P1,P2,CNT1,IDATA,OK
    REAL(PDB)  :: DOLD,LENGTH,LENGTH_TOT
    LOGICAL    :: LGON
    INTEGER,   POINTER, DIMENSION(:) :: ITRAC3
    REAL(PDB), POINTER, DIMENSION(:) :: RTRAC3
    INTEGER, PARAMETER :: FOUT =6
    !***
    lgon=.true.
    !     initiate counter
    CNT1=CNT
    LENGTH_TOT=0.
    EXTERIOR : DO WHILE(LGON)
       NB_TOT=NB_TOT+1
       IF(NB_TOT > NB_MAX) CALL XABORT('SAL240_4_2: NB_TOT overflow')
       !        set horizontal angle index in angtab
       ANGTAB(2*NB_TOT)=ANGLE_TO_NUMBER(EX,EY,DANGLT)
       ELMTAB(2*NB_TOT-1)=LNEW ; ELMTAB(2*NB_TOT)=0 ;
       TORIG(1,NB_TOT)=AX+DNEW*EX ; TORIG(2,NB_TOT)=AY+DNEW*EY ;
       !
       !*       track a sub-trajectory         
       INTERIOR: DO WHILE(NNEW.GT.0)
          !           update data to compute next node:
          DOLD=DNEW
          LOLD=LNEW
          NOLD=NNEW
          !           crossing node NOLD
          !           input: trajectory (t):r=a+d*e => a = (ax,ay), e = (ex,ey)
          !           DOLD = d at last intersection
          !           NOLD = node just entered
          P1=PPERIM(NOLD); P2=PPERIM(NOLD+1)-1
          CALL SAL241(PERIM(P1:P2),P2-P1+1,IPAR,RPAR,DOLD,NOLD,DNEW,NNEW,LNEW)
          !           at return from SAL241:
          !           DNEW        = d at point exiting node
          !           COSINE      = cosine of trajectory with exiting normal
          !           NNEW        = new node entered
          !           LNEW        = element crossed when exiting node
          !           NBER        = nber of intersections with perimeter
          !           LGOK        = .true. if trajectory exits the node
          IF(.NOT.LGOK) THEN
             IERR=0
             WRITE(FOUT,'("SAL240_4_2 ==> couldnt exit node ",I5)') NOLD
             RETURN
          ENDIF
          !           store data
          LENGTH=DNEW-DOLD
          !           add to total length
          LENGTH_TOT=LENGTH_TOT+LENGTH
          IF(LENGTH.LT.EPS1) CYCLE
          CNT=CNT+1
          IF(CNT>=NMAX2) THEN
             ALLOCATE(ITRAC3(4*NMAX2),RTRAC3(2*NMAX2),STAT=OK)
             IF(OK/=0) CALL XABORT('SAL240_4_2: NMAX2 overflow.')
             RTRAC3(:NMAX2)=RTRAC2(:NMAX2)
             ITRAC3(:2*NMAX2)=ITRAC2(:2*NMAX2)
             DEALLOCATE(RTRAC2,ITRAC2)
             RTRAC2=>RTRAC3
             ITRAC2=>ITRAC3
             NMAX2=2*NMAX2
          ENDIF
          RTRAC2(CNT)=LENGTH
          ITRAC2(CNT+NMAX2)=LNEW
          ITRAC2(CNT)=NOLD
       ENDDO INTERIOR
       !
       !*       exiting motif and analyzing bc condition
       LGON=NNEW<=G_BC_TYPE(1).AND.NNEW>=G_BC_TYPE(5).AND.(ABS(LENGTH_TOT*LENGTH_INV_CYCL-1.)>EPS1)
       !        store nber of regions to angtab
       ANGTAB(2*NB_TOT-1)=CNT-CNT1
       CNT1=CNT
       IF(LGON)THEN
           !     treat boundary condition and get new entering point
           IF(PRESENT(IDATA_BC2).AND.PRESENT(DIST_AXIS)) THEN
               IDATA=IDATA_BC2(IBC2_ELEM(LNEW))
               !     treat bondary condition
               CALL SAL247_3(BCDATA(:,IDATA))
               !     at return from SAL247_3:
               !        AX AY        = entering point at boundary
               !        EX EY        = new direction
               !        compute entering point
               IF(NNEW==G_BC_TYPE(3).OR.NNEW==G_BC_TYPE(4).OR.NNEW==G_BC_TYPE(2)) THEN
                 P1=PPERIM_MAC2(N_AXIS); P2=PPERIM_MAC2(N_AXIS+1)-1
                 !        re-entering point is on the axis
                 CALL SAL241_2(P2-P1+1,PERIM_MAC2(P1:P2),DIST_AXIS(P1:P2),IPAR)
               ENDIF
            ELSE
               CALL XABORT('SAL240_4_2: missing IDATA_BC2 or DIST_AXIS argument')
            ENDIF
       ENDIF
       !
    ENDDO EXTERIOR
    !     set success flag
    IERR=1
    !
  END SUBROUTINE SAL240_4_2
  !
  SUBROUTINE SAL247_3(BCDATA)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! treatment of boundary conditions for the cyclic case
    !
    !Parameters: input
    ! BCDATA  boundary condition descriptor
    !
    !---------------------------------------------------------------------
    !
    USE SAL_TRACKING_TYPES, ONLY : NNEW,DNEW,AX,AY,EX,EY,DELX,N_AXIS
    USE SAL_GEOMETRY_TYPES, ONLY : TYPGEO
    IMPLICIT NONE
    REAL(PDB), INTENT(IN), DIMENSION(:) :: BCDATA
    !***
    REAL(PDB) :: AUX,ACX,ACY,COSTHE,SINTHE,CX,CY
    !***
    !     boundary point:
    AX=AX+DNEW*EX
    AY=AY+DNEW*EY
    SELECT CASE (NNEW)
       CASE(-4)
       !           symmetry with respect an axis (get axis data)
       !           axis: r=c+t*f (f unit vector of angle theta)
       CX=BCDATA(1)
       CY=BCDATA(2)
       COSTHE=BCDATA(3)
       SINTHE=BCDATA(4)
       ACX=AX-CX
       ACY=AY-CY
       AUX=2._PDB*(ACX*COSTHE+ACY*SINTHE)
       AX=AUX*COSTHE-ACX+CX
       AY=AUX*SINTHE-ACY+CY
       DELX=SQRT(AX*AX+AY*AY)
       AUX=2._PDB*(EX*COSTHE+EY*SINTHE)
       EX=AUX*COSTHE-EX
       EY=AUX*SINTHE-EY
       SELECT CASE(TYPGEO)
          CASE(6)
          IF(BCDATA(5)>0.) THEN
             !                 vertical axes
             IF(BCDATA(1)>0) THEN
                N_AXIS=4
             ELSE
                N_AXIS=2
             ENDIF
             DELX=AY
          ELSE
             !                 horizontal axes
             IF(BCDATA(2)>0) THEN
                N_AXIS=3
             ELSE
                N_AXIS=1
             ENDIF
             DELX=AX
          ENDIF
          CASE(7) 
          IF (BCDATA(1)>0) THEN
             !                 cx>0: axis 3 (vertical axis)
             N_AXIS=3
             DELX=AY
          ELSEIF (BCDATA(5)>0.) THEN
             !                 cx=0, angle>0: axis 2 (axis of angle pi/4)
             N_AXIS=2
             DELX=SQRT(AX*AX+AY*AY)
          ELSE
             !                 cx=0, angle=0: axis 1 (axis x)
             N_AXIS=1
             DELX=SQRT(AX*AX+AY*AY)
          ENDIF
          CASE DEFAULT
          CALL XABORT('SAL247_3: option not available(1)')
       END SELECT
       CASE DEFAULT
       CALL XABORT('SAL247_3: option not available(2)')
    END SELECT
    !
  END SUBROUTINE SAL247_3
  !
  INTEGER FUNCTION ANGLE_TO_NUMBER(EX,EY,DANGLT)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! search order number of a horizontal angle in the angular quadrature
    ! formula set
    !
    !Parameters: input
    ! EX       angle cosine
    ! EY       angle sine
    ! DANGLT   angle cosines table
    !
    !---------------------------------------------------------------------
    !
    IMPLICIT NONE
    REAL(PDB), INTENT(IN)               :: EX,EY
    REAL(PDB), INTENT(IN), DIMENSION(:,:) :: DANGLT
    !***
    INTEGER   :: I,NPHI
    REAL(PDB) :: EXREF,EYREF
    !***
    NPHI=SIZE(DANGLT,2)
    ANGLE_TO_NUMBER=0
    DO I=1,NPHI
       EXREF=DANGLT(1,I) ; EYREF=DANGLT(2,I) ;
       IF((ABS(EX-EXREF)<1.E-3).AND.(ABS(EY-EYREF)<1.E-3)) THEN
         ANGLE_TO_NUMBER=I
         GO TO 10
       ELSE IF((ABS(EX-EXREF)<1.E-3).AND.(ABS(EY+EYREF)<1.E-3)) THEN
         ANGLE_TO_NUMBER=2*NPHI-I+1
         GO TO 10
       ENDIF
    ENDDO
    CALL XABORT('ANGLE_TO_NUMBER: FAILURE')
    10 RETURN
    !
  END FUNCTION ANGLE_TO_NUMBER
  !
  SUBROUTINE SAL247_1(RPAR,IPAR,D,SINPHI,COSPHI,EX,EY)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! computes COSPHI and SINPHI for the intersection of the trajectory
    ! with element LNEW of descriptors RPAR and IPAR. SINPHI is
    ! computed with respect to directions outgoing with the trajectory
    !
    !Parameters: input
    ! RPAR         floating point geometry descriptors
    ! IPAR         integer geometry descriptors
    ! D            distance measured along the trajectory
    ! EX           first exiting unit vector in along trajectory
    ! EY           second exiting unit vector in along trajectory
    !
    !Parameters: output
    ! SINPHI       cosine at intersection
    ! COSPHI       sine at intersection
    !
    !---------------------------------------------------------------------
    !
    USE SAL_TRACKING_TYPES, ONLY : AX,AY
    IMPLICIT NONE
    INTEGER,   INTENT(IN), DIMENSION(:) :: IPAR
    REAL(PDB), INTENT(IN), DIMENSION(:) :: RPAR
    REAL(PDB), INTENT(IN)               :: D,EX,EY
    REAL(PDB), INTENT(OUT)              :: SINPHI,COSPHI
    !>    AX AY            = components of origin of trajectory
    !***
    INTEGER   :: TYPE
    REAL(PDB) :: NX,NY
    !***
    NX=0._PDB
    NY=0._PDB
    TYPE=IPAR(1)
    SELECT CASE (TYPE)
       CASE (1)
       !        TYPE=1=>  segment (s): R=C+T*F with T in (0,1)
       !        RPAR(1),RPAR(2) = C = (CX,CY)
       !        RPAR(3),RPAR(4) = F = (FX,FY)
       !        components of normal
       NX=RPAR(4)/RPAR(5)
       NY=-RPAR(3)/RPAR(5)
       CASE(2,3)
       !        TYPE=2,3=> arc of circle
       NX=(AX+D*EX-RPAR(1))/RPAR(3)
       NY=(AY+D*EY-RPAR(2))/RPAR(3)
    END SELECT
    SINPHI=NX*EY-NY*EX
    IF((EX*NX+EY*NY)<0._PDB) SINPHI=-SINPHI
    COSPHI=SQRT(1._PDB-SINPHI*SINPHI)
    !
  END SUBROUTINE SAL247_1
  !
  SUBROUTINE SAL241(PERIM,NPERIM,IPAR,RPAR,DOLD,NOLD,DNEW,NNEW,LNEW)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! computes intersection of trajectory (t): R=A+D*E with a perimeter
    ! composed of the elements given in array perim
    !
    !Parameters: input
    ! PERIM        elements in the perimeter of a node
    ! NPERIM       number of elements in the perimeter
    ! RPAR         floating point geometry descriptors
    ! IPAR         integer geometry descriptors
    ! DOLD         value of D at current position
    ! NOLD         current entered node
    !
    !Parameters: output
    ! DNEW         value of D at the intersection
    ! NNEW         node index that is enter after intersection
    ! LNEW         element index that is intersected
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES,  ONLY : NRPAR,NIPAR
    USE SAL_TRACKING_TYPES,  ONLY : NRPART,NIPART,RPART,IPART,EPS1
    !***
    IMPLICIT NONE
    INTEGER,   INTENT(IN), DIMENSION(:)   :: PERIM
    INTEGER,   INTENT(IN), DIMENSION(:,:) :: IPAR
    REAL(PDB), INTENT(IN), DIMENSION(:,:) :: RPAR
    INTEGER,   INTENT(IN)                 :: NPERIM,NOLD
    INTEGER,   INTENT(OUT)                :: NNEW,LNEW
    REAL(PDB), INTENT(IN)                 :: DOLD
    REAL(PDB), INTENT(OUT)                :: DNEW
    !***
    REAL(PDB) :: D
    INTEGER   :: I,L,INDEX,N,TYPE
    !     INFTY is used to initialize search for minimum distance
    REAL(PDB), PARAMETER :: INFTY=1.E+10
    INTEGER, PARAMETER :: FOUT =6
    !****
    !     initialize distance for intersection
    DNEW=INFTY
    NBER=0
    !     intersection:
    DO I=1,NPERIM
       !        get order nber of element in the perimeter
       L=PERIM(I)
       NBINTE=IPART(1,L)
       IF(NBINTE<0)THEN
          ! compute and store intersections
          TYPE=IPAR(1,L)
          IF(TYPE==1)THEN
             !              segment
             CALL SAL242(RPAR(:,L),IPAR(:,L),RPART(:,L),IPART(2:,L))
          ELSEIF(TYPE<=3)THEN
             !              arc of circle or circle
             CALL SAL243(RPAR(:,L),IPAR(:,L),RPART(:,L),IPART(2:,L))
          ELSE
             CALL XABORT('SAL241: Not implemented')
          ENDIF
          IPART(1,L)=NBINTE
       ENDIF
       !
       IF(NBINTE/=0)THEN
          DO INDEX=1,NBINTE
             D=RPART(INDEX,L)
             N=IPART(1+INDEX,L)
             ! analyzes feasability of intersection to eliminate conflicts due
             ! to concavities and crossing of mesh points
             IF(N/=NOLD.AND.D>(DOLD-EPS1))THEN
                NBER=NBER+1
                IF(D<(DNEW-EPS1))THEN
                   DNEW=D
                   COSINE=RPART(INDEX+2,L)
                   NNEW=N
                   LNEW=L
                ELSEIF(D<(DNEW+EPS1))THEN
                   ! case of two close intersections : liquidate smaller
                   NBER=NBER-1
                   IF(D>DNEW)THEN
                      DNEW=D
                      COSINE=RPART(INDEX+2,L)
                      NNEW=N
                      LNEW=L
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    !
    LGOK=DNEW/=INFTY
    IF(LGOK)THEN
       ! eliminate intersection
       IPART(1,LNEW)=IPART(1,LNEW)-1
       IF(IPART(1,LNEW)==1.AND.RPART(1,LNEW)==DNEW)THEN
          ! FOR AN ELEMENT WITH TWO INTERSECTIONS, MOVE 2ND
          ! INTERSECTION INTO FIRST IF FIRST HAS BEEN TAKEN
          RPART(1,LNEW)=RPART(2,LNEW)
          RPART(3,LNEW)=RPART(4,LNEW)
          IPART(2,LNEW)=IPART(3,LNEW)
       ENDIF
    ELSE
       ! print out problem
       IF(NOLD/=-100) THEN
          WRITE(FOUT,*)'Problem in SAL241: NOLD, DOLD, NBINTE ',NOLD,DOLD,NBINTE
       ENDIF
    ENDIF
  END SUBROUTINE SAL241
  !
  SUBROUTINE SAL242(RPAR,IPAR,D,NODE)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! analysis of the intersection of the trajectory (T):R=A+D*E and
    ! a segment
    !
    !Parameters: input
    ! RPAR         floating point geometry descriptors
    ! IPAR         integer geometry descriptors
    !
    !Parameters: output
    ! D            value of D at intersection
    ! NODE         order nber of node entered after intersection
    !
    !---------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER,   INTENT(IN),  DIMENSION(:) :: IPAR
    INTEGER,   INTENT(OUT), DIMENSION(:) :: NODE
    REAL(PDB), INTENT(IN),  DIMENSION(:) :: RPAR
    REAL(PDB), INTENT(OUT), DIMENSION(:) :: D
    !     DIMENSION        RPAR(*),IPAR(*),NODE(*),D(*)
    !***
    REAL(PDB) :: CAX,CAY,FX,FY,A,DELTAM,DELTA
    REAL(PDB), PARAMETER :: EPS2=0.
    !***
    !     TYPE=1=>  segment (S): R=C+T*F with T in (0,1)
    !     RPAR(1),RPAR(2) = C = (CX,CY)
    !     RPAR(3),RPAR(4) = F = (FX,FY)
    !     components of vector F
    FX=RPAR(3)
    FY=RPAR(4)
    !     DELTA=F X E
    DELTA=FX*EY-FY*EX
    IF(DELTA/=0._PDB)THEN
       !        components of vector CA=C-A
       CAX=RPAR(1)-AX
       CAY=RPAR(2)-AY
       !        A=AC X E
       A=CAY*EX-CAX*EY
       DELTAM=DELTA-A
       IF(DELTAM>(-EPS2).AND.A>(-EPS2))THEN
          !           crossing into + halfspace
          NODE(1)=IPAR(3)
       ELSEIF(DELTAM<=EPS2.AND.A<=EPS2)THEN
          !           crossing into - halfspace
          NODE(1)=IPAR(2)
       ELSE
          !           out-of-range crossing
          NBINTE=0
          RETURN
       ENDIF
       !        compute distance to intersection along trajectory
       !        D = (AC X F)/DELTA
       D(1)=(CAY*FX-CAX*FY)/DELTA
       D(3)=ABS(DELTA)
       NBINTE=1
    ELSE
       !        A/=0 => out-of-range crossing (infinity)
       !        A==0 => trajectory coincides with segment
       !        in any case neglect intersection
       NBINTE=0
    ENDIF
    !
  END SUBROUTINE SAL242
  !
  SUBROUTINE SAL243(RPAR,IPAR,D,NODE)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! analysis of the intersection of the trajectory (T):R=A+D*E and
    ! a circle or arc of circle
    !
    !Parameters: input
    ! RPAR         floating point geometry descriptors
    ! IPAR         integer geometry descriptors
    !
    !Parameters: output
    ! D            value of D at intersection
    ! NODE         order nber of node entered after intersection
    !
    !---------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER,   INTENT(IN),  DIMENSION(:) :: IPAR
    INTEGER,   INTENT(OUT), DIMENSION(:) :: NODE
    REAL(PDB), INTENT(IN),  DIMENSION(:) :: RPAR
    REAL(PDB), INTENT(OUT), DIMENSION(:) :: D
    !     DIMENSION        RPAR(*),IPAR(*),NODE(*),D(*)
    !***
    INTEGER   :: TYPE,I
    REAL(PDB) :: CAX,CAY,RAD,RAD2,RHOMI2,DMIN,THETA,THETA1,THETA2, &
         COSTHE,DELTA
    REAL(PDB), PARAMETER :: EPS2=0._PDB
    !***
    !     TYPE=2,3=> arc of circle (C): R=C+R*F(THETA),
    !                THETA in (THETA1,THETA2)
    !     RPAR(1),RPAR(2) = C = (CX,CY)
    !     RPAR(3)         = R = RADIUS
    !     RPAR(4),RPAR(5) = (THETA1,THETA2) in (0,2PI) with THETA1<THETA2
    !                        and THETA1<0 if arc crosses THETA=0
    !     LGTYPE = .TRUE.  => THETA1 > 0
    !            = .FALSE. => THETA1 < 0
    !     components of vector CA=C-A
    CAX=RPAR(1)-AX
    CAY=RPAR(2)-AY
    TYPE=IPAR(1)
    !     value of R2
    RAD=RPAR(3)
    RAD2=RAD**2
    !     RHOMI2=(CA X E)**2
    RHOMI2=(CAX*EY-CAY*EX)**2
    IF(RAD2>=RHOMI2)THEN
       DELTA=SQRT(RAD2-RHOMI2)
       !        tangent point = two very close points (to avoid infinite loop)
       if(delta==0._pdb) delta=small
       !        DMIN=CA*E
       DMIN=CAX*EX+CAY*EY
       !        D(1) D(2) = min and max distances for the two intersections
       D(1)=DMIN-DELTA
       D(2)=DMIN+DELTA
       D(3)=DELTA/RAD
       D(4)=D(3)
       IF(TYPE==2)THEN
          !           full circle. both intersections are possible
          NODE(1)=IPAR(2)
          NODE(2)=IPAR(3)
          NBINTE=2
       ELSE
          !           analysis for arc of circle:
          THETA1=RPAR(4)
          THETA2=RPAR(5)
          NBINTE=0
          !           compute angles for closest and farthest intersections
          !           and check feasability
          DO I=1,2
             COSTHE=(D(I)*EX-CAX)/RAD
             THETA=SALACO(COSTHE,D(I)*EY-CAY)
             IF( (((THETA1-THETA)<EPS2).AND.((THETA-THETA2)<EPS2)).OR. &
                 ((THETA-THETA2)<(EPS2-TWOPI)) )THEN
                NBINTE=NBINTE+1
                IF(NBINTE/=I)D(NBINTE)=D(I)
                NODE(NBINTE)=IPAR(I+1)
             END IF
          ENDDO
       ENDIF
    ELSE
       NBINTE=0
    ENDIF
    !
  END SUBROUTINE SAL243
  !
END MODULE SAL_TRAJECTORY_MOD
