!
!---------------------------------------------------------------------
!
!Purpose:
! To store common variables for the tracking calculation in SALT:
! module.
!
!Copyright:
! Copyright (C) 2001 Ecole Polytechnique de Montreal.
!
!Author(s):
! X. Warin
!
!---------------------------------------------------------------------
!
MODULE SAL_TRACKING_TYPES
  USE PRECISION_AND_KINDS, ONLY : PDB
  IMPLICIT NONE
  !       NBINT  = max number of intersections trajectory-element
  !       NIPART = number of intersections + one index = effective nber of
  !                intersections to be analyzed
  !       NRPART = 2 * nbe of intersections
  !       NNN    = length of header in a trajectory
  !       MINLEN = smallest length for a trajectory before emptying buffer
  !*      general angular quadrature formula
  !       LGSN   = .true. : use sn formulae
  !       ND_TURN= in case of symmetrical geometry, number of directions obtained
  !                from an angle by application of symmetry conditions.
  !                ND_TURN=1 in case of TISO tracking.
  !                
  INTEGER, PARAMETER :: NBINT=2, NIPART=1+NBINT, NRPART=2*NBINT, NNN=8, MINLEN=20
  INTEGER            :: NMAX2
  INTEGER            :: PRTIND
  !*      tracking data buffer:
  !       integers
  !       itrac2(nmax or 2*nmax) = integer tracking array
  !       itrac3(nmax or 2*nmax) = integer tracking array
  !
  !       *integer descriptors in itrac2 or itrac3 :
  !           1 = address of last data
  !           2 = total nber of sub-trajectories
  !           3 = 
  !           4 = phi for trajectory
  !          2d pca :             2d rac :
  !           5 = left cone        5 = left surface
  !           6 = right cone       6 = left horizontal ray
  !           7 = left phi (2D)    7 = right surface
  !           8 = right phi (2D)   8 = right horizontal ray
  !
  !       reals
  !       RTRAC2(NMAX OR 2*NMAX) = real tracking array (2d, with 3d geo)
  !       RTRAC3(NMAX OR 2*NMAX) = real tracking array
  !
  !       *real descriptors:
  !           1 =
  !           2 = cos phi entering basic
  !           3 = sin phi left surface
  !           4 = sin phi right surface
  !           5 = cos phi left surface
  !           6 = cos phi right surface
  !           7 = total weight (delr*wphi for 2D)
  !           8 = radial weight (DELR for 2D)
  !
  !      *contents of each piece of trajectory descriptor arrays:
  !        ITRACK: (small buffer)+(i1 i2 i3 ... in)*nt+(n_1 a_1 ... n_i a_i)
  !        RTRACK: (small buffer)+(l1 l2 l3 ... ln)*nt+(free space)
  !       where  NT = total subpieces of trajectories
  !            N_I = nber of regions intersected by horizontal angle a_i
  !            A_I = horizontal angle number for ith (1<i<nt) subpiece trajectory
  !        ITRACK(K) = ik, k'th intersected region
  !        RTRACK(K) = lk, chord length in the k'th intersected region
  !
  !       integers
  !       IPART(NIPART,MXELEM)   = to store integer intersection data
  !       ANGTAB(2)              = to store (1.number of regions crossed,2.angle nber)
  !                                for each sub-trajectory
  !       ELMTAB(2)              = to store the first element entered for each
  !       ANGTAB(2*ND_TURN)      = to store (1.nber of regions crossed,2.angle nber)
  !                                for each sub-trajectory
  !       ELMTAB(2*ND_TURN)      = to store the first element entered for each
  !                                sub-trajectory,for the dsa method
  !                                sub-trajectory,for the dsa method
  !       reals*8
  !       RPART(NRPART,MXELEM)   = to store real intersection data
  !       DPIECE(2    )          = to store projection of discontinuities
  !                                over line  orthogonal to tracking
  !                                 (with gauss radial quadrature)
  !       RPAR_MORE(5,2)         = to store re-entering points
  !                                RPAR_MORE(1:2,I) = (AX,AY) ith re-entering point
  !                                RPAR_MORE(3:4,I) = (EX,EY) ith re-entering track direction
  !                                RPAR_MORE(5,I)   = d,distance to (ax,ay)
  !       EPS1                   = epsilon to compute intersections element-trajectory
  !       DELR                   = delta r displacement for numerical quadrature
  !       LENGTH_INV_CYCL        = inversion of length of a cyclic trajectory
  !       TORIG(2,ND_TURN)       = to store the entering point for each sub-trajectory
  !
  INTEGER,   ALLOCATABLE, DIMENSION(:)     :: ANGTAB,ELMTAB
  INTEGER,   ALLOCATABLE, DIMENSION(:,:)   :: IPART
  INTEGER,   POINTER, DIMENSION(:)         :: ITRAC2
  REAL(PDB), POINTER, DIMENSION(:)         :: RTRAC2
  REAL(PDB), DIMENSION(2)     :: DPIECE
  REAL(PDB), ALLOCATABLE, DIMENSION(:,:)   :: RPART,TORIG
  !
  INTEGER   :: NBTRAC,CNT0,CNT
  !       NB_TOT     = total number of sub-trajectories in a trajectory
  !       NB_MAX     = maximum number of sub-trajectories in a trajectory
  !       length_inv_cycl= inversion of length of a cyclic trajectory
  INTEGER   :: NITER,NBER,NBINTE,NTRACK,NNEW,LNEW,IERR,NB_TOT,NB_MAX,N_AXIS,MXSEG
  LOGICAL   :: LGTYPE,LGOK,LGMORE
  REAL(PDB) :: DELMIN,DD0,COSINE,DNEW,DELX,EX0,EY0,AX,AY,EX,EY,R, &
       ALPHA,F0,D0,AT,BT,DINIT,EPS1,DELR,LENGTH_INV_CYCL,HX(6),HY(6), &
       BX(6),BY(6)
END MODULE SAL_TRACKING_TYPES
