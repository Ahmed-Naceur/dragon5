!
!---------------------------------------------------------------------
!
!Purpose:
! Support module used to create a 2D geometry.
!
!Copyright:
! Copyright (C) 2001 Ecole Polytechnique de Montreal.
!
!Author(s):
! X. Warin
!
!---------------------------------------------------------------------
!
MODULE SAL_GEOMETRY_MOD

  USE SAL_GEOMETRY_TYPES
  USE PRECISION_AND_KINDS, ONLY : PDB, PI,TWOPI,HALFPI
  USE SAL_NUMERIC_MOD,    ONLY : SAL141
  USE SALGET_FUNS_MOD
  TYPE(T_G_BASIC), POINTER :: GG

CONTAINS

  SUBROUTINE SAL100()
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! perform input data and first allocation for geometry OBJECT T_G_BASIC
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES, ONLY : TYPGEO,NBFOLD,EPS,ANGGEO, &
         INDEX,KNDEX,PREC,ISPEC,LGSPEC
    USE SAL_TRACKING_TYPES, ONLY : PRTIND
    !****
    IMPLICIT NONE
    !*****
    ! local variable
    ! **************
    INTEGER, PARAMETER :: N_DATAIN=25, N_DATARE=20
    INTEGER, DIMENSION (N_DATAIN) :: DATAIN
    REAL,    DIMENSION (N_DATARE) :: DATARE
    INTEGER :: NBFLUX
    INTEGER, DIMENSION(:), ALLOCATABLE :: IFLUX
    INTEGER :: OK
    INTEGER, PARAMETER :: FOUT =6
    !*****
    !     TYPGEO = type of geometry:
    !              0 = normal
    !              1 = geometry with symmetries of two axis of angle pi/n,n>0
    !              2 = geometry with rotation of angle 2*pi/n,n>1
    !              5 = rectangular geometry with translation on all sides
    !              6 = rectangular geometry with symmetry on all sides
    !              7 = 1/8 assembly with symmetries on all sides
    !              8 = hexagonal geometry with symmetries on all sides
    !                  (reduced to a triangle)
    !              9 = hexagonal geometry with translations on all sides
    !     NBFOLD = n in angle definition of rotation or symmetry geometry
    !     NB_NODE = number of nodes
    !     NBELEM = number of elements
    !
    CALL SALGET(DATAIN,6,F_GEO,FOUT0,'dimensions for geometry')
    TYPGEO=DATAIN(1)
    NBFOLD=DATAIN(2)
    GG%NB_NODE=DATAIN(3)
    GG%NB_ELEM=DATAIN(4)
    NBFLUX=DATAIN(6)
    SELECT CASE (TYPGEO)
       CASE(1)
       ANGGEO=TWOPI/NBFOLD
       CASE(2)
       ANGGEO=TWOPI/NBFOLD
       CASE(5:6)
       ANGGEO=HALFPI
       CASE(7)
       ANGGEO=PI*0.25
       CASE(8:9)
       ANGGEO=PI/3.
    END SELECT
    IF(PRTIND >= 1) THEN
      WRITE(FOUT,*) 'SAL100: TYPGEO=',TYPGEO,' NBFOLD=',NBFOLD
    ENDIF
    LGSPEC=(TYPGEO/=0).AND.(NBFOLD==0)
    IF(LGSPEC) THEN
      IF(ISPEC==0) THEN
        WRITE(*,*) 'SAL100: TYPGEO=',TYPGEO,' NBFOLD=',NBFOLD
        CALL XABORT('SAL100: TISO option is incompatible with the surfacic file')
      ENDIF
    ELSE
      IF(ISPEC==1) THEN
        WRITE(*,*) 'SAL100: TYPGEO=',TYPGEO,' NBFOLD=',NBFOLD
        CALL XABORT('SAL100: TSPC option is incompatible with the surfacic file')
      ENDIF
    ENDIF
    !
    !*    read printing indexes for general domain data and topological deformations
    !     INDEX  =  to print general domain data
    !     KNDEX  =  to print motions of topological adjustment
    !     PREC   =  if 0 then read RPAR & BCDATA with e20.12)
    CALL SALGET(DATAIN,3,F_GEO,FOUT0,'index kndex prec')
    INDEX=DATAIN(1)  
    KNDEX=DATAIN(2) 
    PREC=DATAIN(3) 
    !
    !*    read epsilons for topological deformations
    !     EPS    = if the distance of two ends of elements < eps, 
    !              they will be united to one point
    CALL SALGET(DATARE,1,F_GEO,FOUT0,'eps')
    EPS=DATARE(1)
    !
    IF(PRTIND >= 1) THEN
        WRITE(FOUT,'(//,5X,''domain checkout:'',/, &
        & 5X,''elements are in contact for distance < '',1P,E12.4,//)') EPS
    ENDIF
    !
    ALLOCATE(GG%NUM_MERGE(GG%NB_NODE),STAT =OK)
    IF(OK /= 0) CALL XABORT('SAL100: failure to allocate NB_NODE')
    CALL SALGET(GG%NUM_MERGE,GG%NB_NODE,F_GEO,FOUT0,'FLUX INDEX PER NODE')
    IF(MAXVAL(GG%NUM_MERGE) /= NBFLUX) CALL XABORT('SAL100: inconsistent NBFLUX')
    CALL SALGET(GG%NAME_GEOM,F_GEO,FOUT0,'NAMES OF MACROS')
    ALLOCATE(IFLUX(NBFLUX))
    CALL SALGET(IFLUX,NBFLUX,F_GEO,FOUT0,'macro order number per flux region')
    DEALLOCATE(IFLUX)
    !*    do the work (SAL100_2 is called here!):
    CALL SAL110(GG)
    !
  END SUBROUTINE SAL100
  !
  SUBROUTINE SAL110(GG)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! constructs geometrical domain
    !
    !Parameters: input/output
    ! GG    geometry descriptor
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES, ONLY : TYPGEO,NBFOLD,NIPAR,NRPAR,ALLSUR,NANIS,ISPEC
    USE SAL_TRACKING_TYPES, ONLY : PRTIND
    !****
    IMPLICIT NONE
    ! in variable
    ! ************
    TYPE(T_G_BASIC), INTENT(INOUT) :: GG

    ! local variable
    ! ***************
    INTEGER :: ELEM,OK,I,TYPE
    CHARACTER(LEN=4) :: HSYM
    REAL(PDB),PARAMETER :: CONV=PI/180._PDB
    INTEGER, PARAMETER :: FOUT =6
    !****
    ! allocate node arrays
    HSYM=' '
    IF((TYPGEO==0).OR.((TYPGEO==6).AND.(NBFOLD==0))) THEN
      CONTINUE
    ELSE IF(((TYPGEO==1).AND.(NBFOLD==8)).OR.((TYPGEO==7).AND.(NBFOLD==0))) THEN
      HSYM='EIGH'
    ELSE IF(((TYPGEO==1).AND.(NBFOLD==4)).OR.((TYPGEO==3).AND.(NBFOLD==0))) THEN
      HSYM='SYME'
    ELSE IF((TYPGEO==1).AND.(NBFOLD==2)) THEN
      HSYM='SY1D'
    ELSE
      WRITE(*,*) "TYPGEO=",TYPGEO," NBFOLD=",NBFOLD
      CALL XABORT('SAL110: non supported type of symmetry')
    ENDIF
    ALLOCATE (GG%IPAR(NIPAR,GG%NB_ELEM), GG%RPAR(NRPAR,GG%NB_ELEM), STAT=OK)
    IF(OK/=0) CALL XABORT('SAL110: not enough memory I,R')

    !*    read surfacic file
    CALL SALINP(GG)
    !
    !*    unite the ends of elements, redefine elements
    CALL SAL128(GG%RPAR,GG%IPAR,GG%NB_ELEM)
    !
    IF(ISPEC==0) THEN
      !*    unfold domain
      IF((LBCDIAG.AND.(TYPGEO==0)).OR.(LBCDIAG.AND.(TYPGEO==7))) THEN
        IF(PRTIND>0) WRITE(*,*) "SAL110: DIAG unfold"
        IF(NANIS>1) CALL XABORT('SAL110: unfold unsupported with NANIS>1')
        CALL SALFOLD('DIAG',GG)
      ELSE IF(HSYM=='EIGH') THEN
        IF(PRTIND>0) WRITE(*,*) "SAL110: DIAG + SYME unfold"
        IF(NANIS>1) CALL XABORT('SAL110: unfold unsupported with NANIS>1')
        CALL SALFOLD('DIAG',GG)
        CALL SALFOLD('SYMX',GG)
        CALL SALFOLD('SYMY',GG)
      ELSE IF(HSYM=='SYME') THEN
        IF(PRTIND>0) WRITE(*,*) "SAL110: SYME unfold"
        IF(NANIS>1) CALL XABORT('SAL110: unfold unsupported with NANIS>1')
        CALL SALFOLD('SYMX',GG)
        CALL SALFOLD('SYMY',GG)
      ELSE IF(HSYM=='SY1D') THEN
        IF(PRTIND>0) WRITE(*,*) "SAL110: SY1D unfold"
        IF(NANIS>1) CALL XABORT('SAL110: unfold unsupported with NANIS>1')
        CALL SALFOLD('SYMX',GG)
      ELSE
        IF(PRTIND>0) WRITE(*,*) "SAL110: no unfold"
      ENDIF
      IF((TYPGEO/=0).AND.(NBFOLD==0)) THEN
        TYPGEO=6
      ELSE
        TYPGEO=0; NBFOLD=0;
      ENDIF
    ENDIF
    IF(PRTIND>0) WRITE(FOUT,*) 'SAL110: after unfolding -- NB_ELEM=',GG%NB_ELEM

    IF(PRTIND>5) THEN
      !*    print surfacic file
      WRITE(FOUT,'(5H--cut,75(1H-))')
      WRITE(FOUT,'(5HBEGIN)')
      WRITE(FOUT,'(42H* typgeo nbfold nbnode nbelem nbmacr nbreg)')
      WRITE(FOUT,'(6I7)') TYPGEO,NBFOLD,GG%NB_NODE,GG%NB_ELEM,1,GG%NB_NODE
      WRITE(FOUT,'(20H* index  kndex  prec)')
      WRITE(FOUT,'(4I7)') 0,0,1
      WRITE(FOUT,'(18H* eps         eps0)')
      WRITE(FOUT,'(1P,2E18.9)') 1.0E-03,1.0E-05
      WRITE(FOUT,'(20H* num_of_region/mesh)')
      WRITE(FOUT,'(10I7)') (I,I=1,GG%NB_NODE)
      WRITE(FOUT,'(18H* name_of_geometry)')
      WRITE(FOUT,'(7H   GEOM)')
      WRITE(FOUT,'(35H* macro_order_index_per_flux_region)')
      WRITE(FOUT,'(10I7)') (1,I=1,GG%NB_NODE)
      DO ELEM=1,GG%NB_ELEM
        TYPE=GG%IPAR(1,ELEM)
        WRITE(FOUT,'(7h elem =,I6)') ELEM
        WRITE(FOUT,'(22H*type    node-   node+)')
        WRITE(FOUT,'(3I6)') (GG%IPAR(I,ELEM),I=1,3)
        WRITE(FOUT,'(63H*cx            cy            ex_or_R       ey_or_theta1  theta2)')
        IF(TYPE<=2) THEN
          WRITE(FOUT,'(1P,5E18.9)') (GG%RPAR(I,ELEM),I=1,5)
        ELSE IF(TYPE==3) THEN
          WRITE(FOUT,'(1P,5E18.9)') (GG%RPAR(I,ELEM),I=1,3),GG%RPAR(4,ELEM)/CONV, &
                                 (GG%RPAR(5,ELEM)-GG%RPAR(4,ELEM))/CONV
        ENDIF
      ENDDO
      WRITE(FOUT,'(40H*defaul  nbbcda  allsur  divsur  ndivsur)')
      WRITE(FOUT,'(1P,5I8)') GG%DEFAUL,GG%NBBCDA,ALLSUR,0,0
      WRITE(FOUT,'(17H*albedo  deltasur)')
      WRITE(FOUT,'(1P,2E18.9)') GG%ALBEDO,0.0
      DO ELEM=1,GG%NBBCDA
        WRITE(FOUT,'(37H particular boundary condition number,i12)') ELEM
        WRITE(FOUT,'(13H*type    nber)')
        WRITE(FOUT,'(1P,2I8)') GG%BCDATAREAD(ELEM)%SALTYPE,GG%BCDATAREAD(ELEM)%NBER
        WRITE(FOUT,'(14H*elems(1,nber))')
        WRITE(FOUT,'(1P,10I8)') (GG%BCDATAREAD(ELEM)%ELEMNB(I),I=1,GG%BCDATAREAD(ELEM)%NBER)
        IF(GG%BCDATAREAD(ELEM)%SALTYPE==0) THEN
          WRITE(FOUT,'(7H*albedo)')
          WRITE(FOUT,'(1P,E18.9)') GG%BCDATAREAD(ELEM)%BCDATA(1)
        ELSE
          WRITE(FOUT,'(22H*cx      cy      angle)')
          WRITE(FOUT,'(1P,3E18.9)') (GG%BCDATAREAD(ELEM)%BCDATA(I),I=1,3)
        ENDIF
      ENDDO
    ENDIF

    ! allocate media and element arrays 
    ALLOCATE (GG%VOL_NODE(GG%NB_NODE),GG%PPERIM_NODE(GG%NB_NODE+1),GG%IBC2_ELEM(GG%NB_ELEM), &
         GG%ISURF2_ELEM(GG%NB_ELEM),GG%MED(GG%NB_NODE), STAT=OK)
    IF(OK/=0) CALL XABORT('SAL110: not enough memory VOL')

    !*    2D boundary conditions and macro contacts:
    !     - defines NB_BC2, NBSUR2
    !     - defines surface strctures for each 2D macro
    !     - defines perimeter structure for each 2D macro
    !     - read 2D boundary conditions
    CALL SAL130(GG)
    !
    !*    topological check
    CALL SAL140(GG%NB_NODE,GG%RPAR,GG%IPAR,GG%PPERIM_NODE,GG%PERIM_NODE)
    !
    !*    volumes, surfaces, put local nbers in node, and read media:
    CALL SAL160(GG)
    IF(PRTIND>5) THEN
      WRITE(FOUT,'(12H* mil(nbreg))')
      WRITE(FOUT,'(10I7)') (GG%MED(I),I=1,GG%NB_NODE)
      WRITE(FOUT,'(3HEND)')
      WRITE(FOUT,'(5H--cut,75(1H-))')
    ENDIF
    !
    !*    printout basic domain
    CALL SAL170(GG)
    !
  END SUBROUTINE SAL110
  !
  SUBROUTINE SALINP(GG)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! read surfacic file
    !
    !Parameters: input/output
    ! GG    geometry descriptor
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES, ONLY : ALLSUR,PREC
    !****
    IMPLICIT NONE
    ! in variable
    ! ************
    TYPE(T_G_BASIC), INTENT(INOUT) :: GG
    !****
    INTEGER, PARAMETER :: N_DATAIN=25
    INTEGER, DIMENSION (N_DATAIN) :: DATAIN
    INTEGER   :: ELEM,I,TYPE,NBER
    INTEGER, PARAMETER, DIMENSION(0:4) :: READ_BC_LEN=(/1,1,2,3,3/)
    INTEGER, PARAMETER :: FOUT =6
    !
    !*    read element data
    DO ELEM=1,GG%NB_ELEM
       CALL SAL126(GG%RPAR(:,ELEM),GG%IPAR(:,ELEM))
    ENDDO
    !*    read
    !     DEFAUL = default bc condition
    !               0 = all external surfaces are isotropically reflected
    !                   with the same ALBEDO = ALBED0
    !               1 = all external surfaces are reflected
    !     NBBCDA = number of groups of bc conditions
    !
    CALL SALGET(DATAIN,3,F_GEO,FOUT0,'general bc data')
    GG%DEFAUL=DATAIN(1)
    GG%NBBCDA=DATAIN(2)
    ALLSUR=DATAIN(3)
    !     we can define only two default bc's
    IF(GG%DEFAUL>4.OR.GG%DEFAUL<0) THEN
         WRITE(FOUT,'(8H defaul=,I5)') GG%DEFAUL
         CALL XABORT('SAL131: wrong default bc type')
    ENDIF
    !*    read albedo : defaul bcdata
    CALL SALGET(GG%ALBEDO,F_GEO,FOUT0,PREC,'GENERAL ALBEDO')
    !
    !*    read detailed bcdata if required (motions)
    LBCDIAG=.FALSE.
    IF(GG%NBBCDA>0)THEN
       ALLOCATE(GG%BCDATAREAD(GG%NBBCDA))
       DO I=1,GG%NBBCDA
          !           read type of bc and nber of elements affected
          !         - TYPE of bc = 0 ~ 5
          !                  0 (vacuum + albedo), 1 (specular reflexion),
          !                  2 (translation),     3 (rotation),
          !                  4 (axial symmetry),  5 (central symmetry)
          !         - NBER = number of elements affected
          !
          CALL SALGET(DATAIN,2,F_GEO,FOUT0,'SPECIFIC BC: TYPE NBER')
          TYPE=DATAIN(1)
          NBER=DATAIN(2)
          GG%BCDATAREAD(I)%SALTYPE=TYPE
          GG%BCDATAREAD(I)%NBER=NBER
          IF(TYPE>5.OR.TYPE<0) CALL XABORT('SAL131: wrong bc type')
          IF(NBER>GG%NB_ELEM) CALL XABORT('SAL131: bc def exceeds nber of elements')
          !
          !           read order numbers of elements affected
          ALLOCATE(GG%BCDATAREAD(I)%ELEMNB(NBER))
          CALL SALGET(GG%BCDATAREAD(I)%ELEMNB,NBER,F_GEO,FOUT0,'BC ELEMENTS')
          !           read bc motion
          CALL SALGET(GG%BCDATAREAD(I)%BCDATA,READ_BC_LEN(TYPE),F_GEO,FOUT0,PREC, &
               'data for specific bc condition')
          LBCDIAG=LBCDIAG.OR.((GG%BCDATAREAD(I)%BCDATA(1)==0._PDB).AND.(GG%BCDATAREAD(I)%BCDATA(2)==0._PDB) &
                         .AND.(GG%BCDATAREAD(I)%BCDATA(3)==45._PDB))
       ENDDO
    ENDIF
  END SUBROUTINE SALINP
  !
  SUBROUTINE SAL126(RPAR,IPAR)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! reads one element data
    !
    !Parameters: input/output
    ! RPAR  floating point geometry descriptors
    ! IPAR  integer geometry descriptors
    !
    !---------------------------------------------------------------------
    !
    !****
    USE SAL_GEOMETRY_TYPES,  ONLY : PREC,G_ELE_TYPE
    !****
    IMPLICIT NONE
    ! in variable
    INTEGER,   INTENT(OUT), DIMENSION(:) :: IPAR
    REAL(PDB), INTENT(INOUT), DIMENSION(:) :: RPAR
    !****
    ! local variable
    INTEGER   :: TYPE,NBER,IAUX
    REAL(PDB) :: PHI1,PHI2,ANGMAX,DELPHI
    !****
    !     IPAR(1)          = 1 (segment), 2 (circle), 3(arc of circle)
    !     IPAR(2)          = order nber of node in - side
    !     IPAR(3)          = order nber of node in + side
    !
    !*    read integer element data
    CALL SALGET(IPAR,3,F_GEO,FOUT0,'integer descriptors')
    TYPE=IPAR(1)
    !
    !     ** segment:
    !        R = C+T*F  with T in (0,1)
    !        RPAR(1~5)  = CX CY FX FY F
    !
    !     ** arc of circle:
    !        R = C+R*F(THETA) with
    !            THETA in (THETA1 < THETA2)
    !            with THETA1 in (0,360)
    !        RPAR(1~5) = CX CY R THETA1 THETA2 (in degrees)
    !
    !*    read real element data

    SELECT CASE (TYPE)
       CASE (G_ELE_TYPE(1))
       NBER=4
       CASE (G_ELE_TYPE(2))
       NBER=3
       CASE (G_ELE_TYPE(3))
       NBER=5
       ANGMAX=360._PDB
    CASE DEFAULT
       WRITE(FOUT0,'(1X,''==> SAL126: unknown type '',I3)') TYPE
       CALL XABORT('SAL126: unknown element type')
    END SELECT

    CALL SALGET(RPAR,NBER,F_GEO,FOUT0,PREC,'real descriptors')

    SELECT CASE (TYPE)
       CASE (G_ELE_TYPE(1))
       !        segment: compute length
       RPAR(5)=SQRT(RPAR(3)*RPAR(3)+RPAR(4)*RPAR(4))
       RPAR(6)=0._PDB
       CASE (G_ELE_TYPE(2))
       !        full circle: set angles
       RPAR(4)=0._PDB
       RPAR(5)=TWOPI
       RPAR(6)=0._PDB
       CASE (G_ELE_TYPE(3),G_ELE_TYPE(4))
       !        check angles
       PHI1=RPAR(NBER-1)
       DELPHI=RPAR(NBER)
       !        order angles in increasing values:
       IF(DELPHI>0._PDB)THEN
          IF(DELPHI>ANGMAX)THEN
             WRITE(FOUT0,'(1X,''==> SAL126: DELPHI = '',1P,E12.4, &
                  &'' > '',1P,E12.4,''  FOR TYPE'',I3)')DELPHI,ANGMAX,TYPE
             CALL XABORT('SAL126: invalid value of delphi')
          ENDIF
          PHI2=PHI1+DELPHI
       ELSE
          IF(DELPHI<-ANGMAX)THEN
             WRITE(FOUT0,'(1X,''==> SAL126: DELPHI = '',1P,E12.4, &
                  &'' < '',1P,E12.4,''  FOR TYPE'',I3)')DELPHI,-ANGMAX,TYPE
             CALL XABORT('SAL126: invalid value of delphi')
          ENDIF
          PHI2=PHI1
          PHI1=PHI1+DELPHI
       ENDIF
       IF(TYPE==G_ELE_TYPE(3))THEN
          !           arc of circle: put phi1 within 0 and 360.
          IF(PHI1>360._PDB)THEN
             IAUX=INT(PHI1/360._PDB)
             DELPHI=360._PDB*IAUX
             PHI2=PHI2-DELPHI
             PHI1=PHI1-DELPHI
          ELSEIF(PHI1<0._PDB)THEN
             IAUX=INT((-PHI1+1.D-7)/360._PDB)+1
             DELPHI=360._PDB*IAUX
             PHI2=PHI2+DELPHI
             PHI1=PHI1+DELPHI
          ENDIF
          RPAR(6)=0._PDB
       ELSE
          CALL XABORT('SAL126: unsupported option')
       ENDIF
       !        convert to radians
       RPAR(NBER-1)=PHI1*(TWOPI/360._PDB)
       RPAR(NBER)=PHI2*(TWOPI/360._PDB)
    END SELECT
    !
  END SUBROUTINE SAL126
  !
  SUBROUTINE SAL128(RPAR,IPAR,NB_ELEM)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! create library of ends of elements and redefine element discriptors
    !
    !Parameters: input/output
    ! RPAR     floating point geometry descriptors
    ! IPAR     integer geometry descriptors
    ! NB_ELEM  number of surfacic elements
    !
    !---------------------------------------------------------------------
    !
    !****
    USE SAL_GEOMETRY_TYPES, ONLY : G_ELE_TYPE,EPS,KNDEX,TYPGEO
    ! in variable
    !************
    IMPLICIT NONE
    INTEGER,   INTENT(IN), DIMENSION(:,:)     :: IPAR
    REAL(PDB), INTENT(INOUT), DIMENSION(:,:)  :: RPAR
    INTEGER , INTENT(IN) :: NB_ELEM
    ! local variable
    ! ***************
    REAL(PDB), DIMENSION (:,:), ALLOCATABLE :: POINT    ! coordinates of the ends of elements
    INTEGER, DIMENSION (:,:), ALLOCATABLE  :: ELEM_NEW  ! point nbers of two ends of elements
    !****
    !     NN               = total number of points
    !     NBP              = count of number of ends having the same coordinates
    INTEGER   :: TYPE,NN,ELEM,I,P,POS,I1,I2,OK
    REAL(PDB) :: X(2),Y(2),D,DX,DY
    REAL      :: EPS2
    LOGICAL   :: LGONE
    INTEGER, DIMENSION(NB_ELEM)     :: NBP
    INTEGER, PARAMETER :: FOUT =6
    !****
    !     - POINT(2,NB_ELEM)    = coordinates of the ends of elements
    !     - ELEM_NEW(2,NB_ELEM) = point nbers of two ends of elements
    !
    ALLOCATE(POINT(2,NB_ELEM),ELEM_NEW(2,NB_ELEM),STAT=OK)
    IF(OK/=0) CALL XABORT('SAL100_3_1: not enough memory I,CH')

    EPS2=EPS*EPS
    IF(KNDEX/=0)WRITE(FOUT,'(//,&
         &" reunite the ends of elements of distance less than ",E13.4)') EPS
    !
    !*    nn counts the nber of points
    NBP(1:NB_ELEM)=0
    NN=0
    DO ELEM=1,NB_ELEM
       TYPE=IPAR(1,ELEM)
       IF(TYPE/=G_ELE_TYPE(2)) THEN
          DO I=1,2
             !              get coordinates of end
             CALL SAL141(TYPE,RPAR(:,ELEM),X(I),Y(I),I)
             !              find if there is a defined neighbouring point
             LGONE=.FALSE.
             DO P=1,NN
                DX=POINT(1,P)-X(I)
                DY=POINT(2,P)-Y(I)
                D=DX*DX+DY*DY
                IF(D<EPS2) THEN
                   LGONE=.TRUE.
                   POS=P
                   EXIT
                ENDIF
             ENDDO
             IF(LGONE) THEN
                !                 choose average value
                X(I)=(X(I)+POINT(1,POS)*NBP(POS))/REAL(NBP(POS)+1)
                Y(I)=(Y(I)+POINT(2,POS)*NBP(POS))/REAL(NBP(POS)+1)
             ELSE
                !                 add points
                NN=NN+1
                IF(NN > NB_ELEM) CALL XABORT('SAL128: point overflow')
                POS=NN
             ENDIF
             POINT(1,POS)=X(I)
             POINT(2,POS)=Y(I)
             NBP(POS)=NBP(POS)+1
             !              define the end of element
             ELEM_NEW(I,ELEM)=POS
          ENDDO
       ELSE
          ELEM_NEW(1:2,ELEM)=0
       ENDIF
    ENDDO
    !*    adjust points on the axes
    SELECT CASE(TYPGEO)
       CASE(1:2)
       CALL SAL128_1(NN,POINT,EPS,EPS2)
       CASE(5:6)
       CALL SAL128_3(NN,POINT,EPS,EPS2)
       CASE(7)
       CALL SAL128_4(NN,POINT,EPS,EPS2)
       CASE(8) 
       CALL SAL128_5(NN,POINT,EPS,EPS2)
       CASE(9) 
       CALL SAL128_6(NN,POINT,EPS,EPS2)
    END SELECT

    !*    redefine the elements in rpar
    IF(KNDEX/=0)WRITE(FOUT,'(/," redefine elements : ",/)')
    DO ELEM=1,NB_ELEM
       TYPE=IPAR(1,ELEM)
       IF(TYPE/=G_ELE_TYPE(2)) THEN
          I1=ELEM_NEW(1,ELEM)
          I2=ELEM_NEW(2,ELEM)
          CALL SAL129(POINT(1,I1),POINT(2,I1),POINT(1,I2),POINT(2,I2), &
               TYPE,RPAR(:,ELEM),ELEM,KNDEX,FOUT)
       ENDIF
    ENDDO
    !
    DEALLOCATE(POINT,ELEM_NEW)
    !
  END SUBROUTINE SAL128
  !
  SUBROUTINE SAL128_1(NN,POINT,EPS,EPS2)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! process geometry with rotation or symmetry (TYPGEO=1 & 2)
    ! adjusts points on the axes;  in case of rotation boundary,adjust the
    ! element lengths on the axis
    !
    !Parameters: input
    ! NN        total number of points
    ! EPS       when distance of points to axes less than EPS, displace
    !           the points onto the axes
    ! EPS2      variable set to EPS*EPS
    !
    !Parameters: input/output
    ! POINT     coordonats of points
    !
    !---------------------------------------------------------------------
    !
    USE PRECISION_AND_KINDS,    ONLY : PDB
    USE SAL_GEOMETRY_TYPES,     ONLY : ANGGEO,TYPGEO,KNDEX
    !****
    IMPLICIT NONE
    INTEGER,   INTENT(IN)     :: NN
    REAL,      INTENT(IN)     :: EPS,EPS2
    REAL(PDB), INTENT(INOUT), DIMENSION(:,:)  :: POINT
    !****
    INTEGER   :: NP,NAXIS,I,J,K,M,P,IPOINT(2,NN),BP(2),IP(2)
    LOGICAL   :: LGADJUST
    REAL(PDB) :: D,AUX,EX(2),EY(2),DIS(2),AX(2),AY(2),D_AXIS(2,NN)
    INTEGER, PARAMETER :: FOUT =6
    !****
    LGADJUST=TYPGEO==2
    !     nber of axes
    NAXIS=2
    !     directions of axes and their distances to (0,0)
    EX(1)=1.; EY(1)=0.; DIS(1)=0.
    EX(2)=COS(ANGGEO); EY(2)=SIN(ANGGEO); DIS(2)=0.
    !     number of intersections
    NP=1
    AX(1)=0.; AY(1)=0.
    !     beginning point of the axes
    IF(LGADJUST) THEN
       BP(1)=1; BP(2)=1; IP=0
    ENDIF
    ITER0:DO P=1,NN
       D=(POINT(1,P)-AX(1))**2+(POINT(2,P)-AY(1))**2
       IF(D.LT.EPS2) THEN
          !           move point to the center of the ordinates
          POINT(1,P)=AX(1) ; POINT(2,P)=AY(1)
          CYCLE ITER0
       ENDIF
       DO I=1,NAXIS
          !           distance to the axis
          D=POINT(1,P)*EY(I)-POINT(2,P)*EX(I)-DIS(I)
          !           when aux>0,the point is near the half axis
          AUX=POINT(1,P)*EX(I)+POINT(2,P)*EY(I)
          IF(ABS(D).LT.EPS.AND.AUX>=0) THEN
             !              move point to the axis
             POINT(1,P)=POINT(1,P)-D*EY(I)
             POINT(2,P)=POINT(2,P)+D*EX(I)
             IF(LGADJUST) THEN
                !                 compute distance to the original of the axis
                D=SQRT((POINT(1,P)-AX(BP(I)))**2+(POINT(2,P)-AY(BP(I)))**2)
                IP(I)=IP(I)+1
                D_AXIS(I,IP(I))=D
                IPOINT(I,IP(I))=P
             ENDIF
             CYCLE ITER0
          ENDIF
       ENDDO
    ENDDO ITER0
    !
    IF(LGADJUST) THEN
       I=1
       IF(IP(I)/=IP(I+1)) CALL XABORT('SAL128_1: axial points nber not the same,axis')
       !        sort the 'D_AXIS' table
       DO J=I,I+1
          DO P=1,IP(J)
             DO K=P+1,IP(J)
                IF(D_AXIS(J,P)>D_AXIS(J,K)) THEN
                   D=D_AXIS(J,P)
                   D_AXIS(J,P)=D_AXIS(J,K)
                   D_AXIS(J,K)=D
                   M=IPOINT(J,P)
                   IPOINT(J,P)=IPOINT(J,K)
                   IPOINT(J,K)=M
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       DO P=1,IP(I)
          IF(ABS(D_AXIS(I,P)-D_AXIS(I+1,P))>EPS) THEN
             IF(KNDEX/=0) &
                  WRITE(FOUT,'(" warning: too great axial length difference",&
                  & 2(/,2X,"axis = ",I3," point = ",I3," d = ",E13.6))')&
                  I,P,D_AXIS(I,P),I+1,P,D_AXIS(I+1,P)
          ENDIF
          D=(D_AXIS(I,P)+D_AXIS(I+1,P))*0.5
          DO J=I,I+1
             K=IPOINT(J,P)
             POINT(1,K)=D*EX(J)+AX(BP(J))
             POINT(2,K)=D*EY(J)+AY(BP(J))
          ENDDO
       ENDDO
    ENDIF
    !
  END SUBROUTINE SAL128_1
  !
  SUBROUTINE SAL128_3(NN,POINT,EPS,EPS2)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! process retangular geometry with translation or symmetry boundary condition
    ! (TYPGEO=5,6): adjusts points on the axes, computes the lengths of the
    ! retangle sides; in case of with translation,element lengths on the opposite
    ! axes will be adjusted to be the same
    !
    !Parameters: input
    ! NN        total number of points
    ! EPS       when distance of points to axes less than EPS, displace
    !           the points onto the axes
    ! EPS2      variable set to EPS*EPS
    !
    !Parameters: input/output
    ! POINT     coordonats of points
    !
    !---------------------------------------------------------------------
    !
    USE PRECISION_AND_KINDS,    ONLY : PDB
    USE SAL_GEOMETRY_TYPES,     ONLY : TYPGEO,KNDEX,LX=>LENGTHX,LY=>LENGTHY
    !****
    IMPLICIT NONE
    INTEGER,   INTENT(IN)     :: NN
    REAL,      INTENT(IN)     :: EPS,EPS2
    REAL(PDB), INTENT(INOUT), DIMENSION(:,:)  :: POINT
    !     DIMENSION        POINT(2,NN)
    !****
    INTEGER   :: NP,NAXIS,I,J,K,M,P,IPOINT(4,NN),BP(4),IP(4)
    LOGICAL   :: LGADJUST
    REAL(PDB) :: D,AUX,EX(4),EY(4),DIS(4),AX(4),AY(4),D_AXIS(4,NN)
    INTEGER, PARAMETER :: FOUT =6
    !****
    !     flag to adjust the element lengths on the opposite axes
    LGADJUST=TYPGEO==5
    !     compute sides of the rectangle
    LX=0.; LY=0.
    DO P=1,NN
       IF(ABS(POINT(1,P))<EPS) THEN
          IF(LY<POINT(2,P)) LY=POINT(2,P)
       ENDIF
       IF(ABS(POINT(2,P))<EPS) THEN
          IF(LX<POINT(1,P)) LX=POINT(1,P)
       ENDIF
    ENDDO
    !*    axis directions and their distance to (0,0)
    NAXIS=4
    EX(1)=1.; EY(1)=0.; DIS(1)=0.
    EX(2)=1.; EY(2)=0.; DIS(2)=-LY
    EX(3)=0.; EY(3)=1.; DIS(3)=0.
    EX(4)=0.; EY(4)=1.; DIS(4)=LX
    NP=4
    AX(1)=0.; AY(1)=0.
    AX(2)=0.; AY(2)=LY
    AX(3)=LX; AY(3)=0.
    AX(4)=LX; AY(4)=LY
    !     beginning point of the axes
    IF(LGADJUST) THEN
       !     axis numbering:
       !              axis 2
       !           ************
       !           *          *
       !    axis 3 *          * axis 4
       !           ************
       !              axis 1
       !
       BP(1)=1; BP(2)=2; BP(3)=1; BP(4)=3
       IP=0
    ENDIF
    ITER0:DO P=1,NN
       DO I=1,NP
          !           if it is the 4 corners of the rectangle
          D=(POINT(1,P)-AX(I))**2+(POINT(2,P)-AY(I))**2
          IF(D.LT.EPS2) THEN
             !              move point to the corners
             POINT(1,P)=AX(I) ; POINT(2,P)=AY(I)
             IF(LGADJUST) THEN
                !                 put distance to the axis origin
                SELECT CASE (I)
                   CASE (2)
                   !                    end of axis 3,keep its distance and number
                   IP(3)=IP(3)+1
                   D_AXIS(3,IP(3))=AY(I)
                   IPOINT(3,IP(3))=P
                   CASE (3)
                   !                    end of axis 1,keep its distance and number
                   IP(1)=IP(1)+1
                   D_AXIS(1,IP(1))=AX(I)
                   IPOINT(1,IP(1))=P
                   CASE (4)
                   !                    end of axes 2&4,keep its distance and number
                   IP(2)=IP(2)+1
                   D_AXIS(2,IP(2))=AX(I)
                   IPOINT(2,IP(2))=P
                   IP(4)=IP(4)+1
                   D_AXIS(4,IP(4))=AY(I)
                   IPOINT(4,IP(4))=P
                END SELECT
             ENDIF
             CYCLE ITER0
          ENDIF
       ENDDO
       DO I=1,NAXIS
          !           distance to the axis
          D=POINT(1,P)*EY(I)-POINT(2,P)*EX(I)-DIS(I)
          !           when aux>0,the point is near the half axis
          AUX=POINT(1,P)*EX(I)+POINT(2,P)*EY(I)
          IF(ABS(D).LT.EPS.AND.AUX>=0) THEN
             !              move point to the axis
             POINT(1,P)=POINT(1,P)-D*EY(I)
             POINT(2,P)=POINT(2,P)+D*EX(I)
             IF(LGADJUST) THEN
                !                 compute distance to the axis origin
                D=SQRT((POINT(1,P)-AX(BP(I)))**2+(POINT(2,P)-AY(BP(I)))**2)
                IP(I)=IP(I)+1
                D_AXIS(I,IP(I))=D
                IPOINT(I,IP(I))=P
             ENDIF
             CYCLE ITER0
          ENDIF
       ENDDO
    ENDDO ITER0
    !
    IF(LGADJUST) THEN
       DO I=1,NAXIS,2
          IF(IP(I)/=IP(I+1)) &
               CALL XABORT('SAL128_3: axial points nber not the same,axis')
          !           sort the 'd_axis' table
          DO J=I,I+1
             DO P=1,IP(J)
                DO K=P+1,IP(J)
                   IF(D_AXIS(J,P)>D_AXIS(J,K)) THEN
                      D=D_AXIS(J,P)
                      D_AXIS(J,P)=D_AXIS(J,K)
                      D_AXIS(J,K)=D
                      M=IPOINT(J,P)
                      IPOINT(J,P)=IPOINT(J,K)
                      IPOINT(J,K)=M
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          DO P=1,IP(I)
             IF(ABS(D_AXIS(I,P)-D_AXIS(I+1,P))>EPS) THEN
                IF(KNDEX/=0) &
                     WRITE(FOUT,'(" warning: too great axial length difference",&
                     & 2(/,2X,"axis = ",I3," point = ",I3," d = ",E13.6))')&
                     I,P,D_AXIS(I,P),I+1,P,D_AXIS(I+1,P)
             ENDIF
             D=(D_AXIS(I,P)+D_AXIS(I+1,P))*0.5
             DO J=I,I+1
                K=IPOINT(J,P)
                POINT(1,K)=D*EX(J)+AX(BP(J))
                POINT(2,K)=D*EY(J)+AY(BP(J))
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    !
  END SUBROUTINE SAL128_3
  !
  SUBROUTINE SAL128_4(NN,POINT,EPS,EPS2)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! process a 1/8 square geometry (typgeo=7): adjusts points on the axes,
    ! computes the square side length
    !
    !Parameters: input
    ! NN        total number of points
    ! EPS       when distance of points to axes less than EPS, displace
    !           the points onto the axes
    ! EPS2      variable set to EPS*EPS
    !
    !Parameters: input/output
    ! POINT     coordonats of points
    !
    !---------------------------------------------------------------------
    !
    USE PRECISION_AND_KINDS,    ONLY : PDB
    USE SAL_GEOMETRY_TYPES,     ONLY : ANGGEO,KNDEX,LX=>LENGTHX,LY=>LENGTHY
    !****
    IMPLICIT NONE
    INTEGER,   INTENT(IN)     :: NN
    REAL,      INTENT(IN)     :: EPS,EPS2
    REAL(PDB), INTENT(INOUT), DIMENSION(:,:) :: POINT
    !****
    INTEGER   :: NP,NAXIS,I,P
    REAL(PDB) :: D,AUX,EX(3),EY(3),DIS(3),AX(3),AY(3)
    INTEGER, PARAMETER :: FOUT =6
    !****
    !     compute the square sides
    LX=0.; LY=0.
    DO P=1,NN
       IF(LY<POINT(2,P)) LY=POINT(2,P)
       IF(LX<POINT(1,P)) LX=POINT(1,P)
    ENDDO
    IF(KNDEX>0.AND.ABS(LY-LX)>EPS) &
         WRITE(FOUT,'(" warning: the square sides are not the same!", &
         &  " lx = ",E12.4," ly = ",E12.4)')LX,LY
    LX=(LX+LY)*0.5
    LY=LX
    !*    adjust points on the axes
    NAXIS=3
    !     axis directions
    EX(1)=1.; EY(1)=0.; DIS(1)=0.
    EX(2)=COS(ANGGEO); EY(2)=SIN(ANGGEO); DIS(2)=0.
    EX(3)=0.; EY(3)=1.; DIS(3)=LX
    NP=3
    !     axis ends
    AX(1)=0.; AY(1)=0.
    AX(2)=LX; AY(2)=0.
    AX(3)=LX; AY(3)=LY
    ITER0:DO P=1,NN
       DO I=1,NP
          D=(POINT(1,P)-AX(I))**2+(POINT(2,P)-AY(I))**2
          IF(D.LT.EPS2) THEN
             !              move point to the axis origin
             POINT(1,P)=AX(I) ; POINT(2,P)=AY(I)
             CYCLE ITER0
          ENDIF
       ENDDO
       DO I=1,NAXIS
          !           distance to the axis
          D=POINT(1,P)*EY(I)-POINT(2,P)*EX(I)-DIS(I)
          !           when aux>0,the point is near the half axis
          AUX=POINT(1,P)*EX(I)+POINT(2,P)*EY(I)
          IF(ABS(D).LT.EPS.AND.AUX>=0) THEN
             !        move point to the axis
             POINT(1,P)=POINT(1,P)-D*EY(I)
             POINT(2,P)=POINT(2,P)+D*EX(I)
             CYCLE ITER0
          ENDIF
       ENDDO
    ENDDO ITER0
    !
  END SUBROUTINE SAL128_4
  !
  SUBROUTINE SAL128_5(NN,POINT,EPS,EPS2)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! process a right trianglar geometry (TYPGEO=8): adjusts points on the axes,
    ! computes the triangular side length
    !
    !Parameters: input
    ! NN        total number of points
    ! EPS       when distance of points to axes less than EPS, displace
    !           the points onto the axes
    ! EPS2      variable set to EPS*EPS
    !
    !Parameters: input/output
    ! POINT     coordonats of points
    !
    !---------------------------------------------------------------------
    !
    USE PRECISION_AND_KINDS,    ONLY : PDB
    USE SAL_GEOMETRY_TYPES,     ONLY : ANGGEO,KNDEX,LX=>LENGTHX,LY=>LENGTHY
    !****
    IMPLICIT NONE
    INTEGER,   INTENT(IN)     :: NN
    REAL,      INTENT(IN)     :: EPS,EPS2
    REAL(PDB), INTENT(INOUT), DIMENSION(:,:)  :: POINT
    !****
    INTEGER   :: NP,NAXIS,I,P
    REAL(PDB) :: D,AUX,EX(3),EY(3),DIS(3),AX(3),AY(3)
    INTEGER, PARAMETER :: FOUT =6
    !****
    !     compute the triangular side
    LX=0.; LY=0.
    DO P=1,NN
       IF(ABS(POINT(2,P))<EPS) THEN
          IF(LX<POINT(1,P)) LX=POINT(1,P)
       ENDIF
       IF(LY<POINT(2,P)) LY=POINT(2,P)
    ENDDO
    IF(KNDEX>0.AND.ABS(LY-LX*SIN(ANGGEO))>EPS) &
         WRITE(FOUT,'(" warning: h is different from a*sin(anggeo)!", &
         &  " a = ",E12.4," h = ",E12.4)')LX,LY
    !     adjust h according a
    LY=LX*SIN(ANGGEO)
    !*    adjust points on the axes
    NAXIS=3
    !     axis directions
    EX(1)=1.; EY(1)=0.; DIS(1)=0.
    EX(2)=COS(ANGGEO); EY(2)=SIN(ANGGEO); DIS(2)=0.
    EX(3)=COS(ANGGEO*2.); EY(3)=SIN(ANGGEO*2.); DIS(3)=LY
    NP=3
    !     axis ends
    AX(1)=0.; AY(1)=0.
    AX(2)=LX; AY(2)=0.
    AX(3)=LX*0.5; AY(3)=LY
    ITER0:DO P=1,NN
       DO I=1,NP
          D=(POINT(1,P)-AX(I))**2+(POINT(2,P)-AY(I))**2
          IF(D.LT.EPS2) THEN
             !              move point to the axis origin
             POINT(1,P)=AX(I) ; POINT(2,P)=AY(I)
             CYCLE ITER0
          ENDIF
       ENDDO
       DO I=1,NAXIS
          !           distance to the axis
          D=POINT(1,P)*EY(I)-POINT(2,P)*EX(I)-DIS(I)
          !           when aux>0,the point is near the half axis
          AUX=POINT(1,P)*EX(I)+POINT(2,P)*EY(I)
          IF(ABS(D).LT.EPS.AND.AUX>=0) THEN
             !              move point to the axis
             POINT(1,P)=POINT(1,P)-D*EY(I)
             POINT(2,P)=POINT(2,P)+D*EX(I)
             CYCLE ITER0
          ENDIF
       ENDDO
    ENDDO ITER0
    !
  END SUBROUTINE SAL128_5
  !
  SUBROUTINE SAL128_6(NN,POINT,EPS,EPS2)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! process a hexagonal geometry with translations on all sides (typgeo=9):
    ! adjusts points on the axes, computes the lengths of hexagon sides;
    ! element lengths on the opposite axes will be adjusted to be the same
    !
    !Parameters: input
    ! NN        total number of points
    ! EPS       when distance of points to axes less than EPS, displace
    !           the points onto the axes
    ! EPS2      variable set to EPS*EPS
    !
    !Parameters: input/output
    ! POINT     coordonats of points
    !
    !---------------------------------------------------------------------
    !
    USE PRECISION_AND_KINDS,    ONLY : PDB
    USE SAL_GEOMETRY_TYPES,     ONLY : ANGGEO,KNDEX,LX=>LENGTHX,LY=>LENGTHY
    !****
    IMPLICIT NONE
    INTEGER,   INTENT(IN)     :: NN
    REAL,      INTENT(IN)     :: EPS,EPS2
    REAL(PDB), INTENT(INOUT), DIMENSION(:,:)  :: POINT
    !****
    INTEGER   :: NP,NAXIS,I,J,K,M,P,IPOINT(6,NN),BP(6),IP(6)
    REAL(PDB) :: D,AUX,EX(6),EY(6),DIS(6),AX(6),AY(6),D_AXIS(6,NN)
    INTEGER, PARAMETER :: FOUT =6
    !****
    !     compute sides of the hexagon: LX=SIDE LENGTH,LY=LX*SIN(A)
    LX=0.; LY=0.
    DO P=1,NN
       IF(ABS(POINT(2,P))<EPS) THEN
          IF(LX<POINT(1,P)) LX=POINT(1,P)
       ENDIF
       IF(LY<POINT(2,P)) LY=POINT(2,P)
    ENDDO
    IF(KNDEX>0.AND.ABS(LY-LX*SIN(ANGGEO))>EPS) &
         WRITE(FOUT,'(" warning: h is different from A*SIN(ANGGEO)!", &
         & " a = ",E12.4," h = ",E12.4)')LX,LY
    !*    axis directions and their distance to (0,0)
    NAXIS=6
    EX(1)=1.; EY(1)=0.; DIS(1)=LY
    EX(2)=1.; EY(2)=0.; DIS(2)=-LY
    EX(3)=COS(ANGGEO); EY(3)=SIN(ANGGEO); DIS(3)=LY
    EX(4)=EX(3);       EY(4)=EY(3);       DIS(4)=-LY
    EX(5)=COS(2.*ANGGEO); EY(5)=SIN(2.*ANGGEO); DIS(5)=LY
    EX(6)=EX(5);          EY(6)=EY(5);          DIS(6)=-LY
    NP=6
    AX(1)=-LX*0.5; AY(1)=-LY
    AX(2)=LX*0.5;  AY(2)=-LY
    AX(3)=LX;      AY(3)=0.
    AX(4)=LX*0.5;  AY(4)=LY
    AX(5)=-LX*0.5; AY(5)=LY
    AX(6)=-LX;     AY(6)=0.
    !     axis numbering:
    !              axis 2
    !              *****
    !      axis 4 *     * axis 5
    !            *   *   *
    !      axis 6 *     * axis 3
    !              *****
    !              axis 1
    !
    BP(1)=1; BP(2)=5; BP(3)=2; BP(4)=6; BP(5)=3; BP(6)=1
    IP=0
    ITER0:DO P=1,NN
       DO I=1,NP
          !           if it is the 6 corners of the hexagon
          D=(POINT(1,P)-AX(I))**2+(POINT(2,P)-AY(I))**2
          IF(D.LT.EPS2) THEN
             !              move point to the corner
             POINT(1,P)=AX(I) ; POINT(2,P)=AY(I)
             !              put distance to the axis origin
             J=0
             SELECT CASE (I)
                CASE (1)
                CYCLE ITER0
                CASE (2)
                !                 END OF AXIS 1
                J=1
                CASE (3)
                !                 END OF AXIS 3
                J=3
                CASE (4)
                !                 END OF AXES 2&5
                J=2
                CASE (5)
                !                 END OF AXIS 4
                J=4
                CASE (6)
                !                 END OF AXIS 6
                J=6
             END SELECT
             !              keep its distance and number
10           IP(J)=IP(J)+1
             D_AXIS(J,IP(J))=LX
             IPOINT(J,IP(J))=P
             IF(J==2) THEN
                J=5
                GOTO 10
             ENDIF
             CYCLE ITER0
          ENDIF
       ENDDO
       DO I=1,NAXIS
          !           distance to the axis
          D=POINT(1,P)*EY(I)-POINT(2,P)*EX(I)-DIS(I)
          !           when aux>0,the point is near the half axis
          AUX=(POINT(1,P)-AX(BP(I)))*EX(I)+(POINT(2,P)-AY(BP(I)))*EY(I)
          IF(ABS(D).LT.EPS.AND.AUX>=0) THEN
             !              move point to the axis
             POINT(1,P)=POINT(1,P)-D*EY(I)
             POINT(2,P)=POINT(2,P)+D*EX(I)
             !              compute distance to the axis origin
             D=SQRT((POINT(1,P)-AX(BP(I)))**2+(POINT(2,P)-AY(BP(I)))**2)
             IP(I)=IP(I)+1
             D_AXIS(I,IP(I))=D
             IPOINT(I,IP(I))=P
             CYCLE ITER0
          ENDIF
       ENDDO
    ENDDO ITER0
    !
    do i=1,naxis,2
       IF(IP(I)/=IP(I+1)) &
            CALL XABORT('SAL128_6: axial points nber not the same,axis')
       !        sort the 'd_axis' table
       DO J=I,I+1
          DO P=1,IP(J)
             DO K=P+1,IP(J)
                IF(D_AXIS(J,P)>D_AXIS(J,K)) THEN
                   D=D_AXIS(J,P)
                   D_AXIS(J,P)=D_AXIS(J,K)
                   D_AXIS(J,K)=D
                   M=IPOINT(J,P)
                   IPOINT(J,P)=IPOINT(J,K)
                   IPOINT(J,K)=M
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       DO P=1,IP(I)
          IF(ABS(D_AXIS(I,P)-D_AXIS(I+1,P))>EPS) THEN
             IF(KNDEX/=0) &
                  WRITE(FOUT,'(" warning: too great axial length difference",&
                  & 2(/,2X,"axis = ",I3," point = ",I3," d = ",E13.6))')&
                  I,P,D_AXIS(I,P),I+1,P,D_AXIS(I+1,P)
          ENDIF
          D=(D_AXIS(I,P)+D_AXIS(I+1,P))*0.5
          DO J=I,I+1
             K=IPOINT(J,P)
             POINT(1,K)=D*EX(J)+AX(BP(J))
             POINT(2,K)=D*EY(J)+AY(BP(J))
          ENDDO
       ENDDO
    ENDDO
    !
  END SUBROUTINE SAL128_6

  SUBROUTINE SAL129(X1,Y1,X2,Y2,TYPE,RPAR,II,KNDEX,FOUT)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! recompute element data by giving two ends of the element
    !
    !Parameters: input
    ! X1      new position in X for the end 1
    ! Y1      new position in Y for the end 1
    ! X2      new position in X for the end 2
    ! Y2      new position in Y for the end 2
    ! TYPE    type of element  1 (segment)  3 (arc of circle)
    ! II      element to be changed
    ! KNDEX   if not 0 print out data when motion>EPS
    ! FOUT    printout file
    !
    !Parameters: input/output
    ! RPAR    descriptors of the element
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES,     ONLY : NRPAR,EPS,G_ELE_TYPE,G_ELE_LEN
    USE SAL_NUMERIC_MOD, ONLY : SALACO
    !****
    IMPLICIT NONE
    ! in variable
    ! ***********
    INTEGER, INTENT(IN)                    :: TYPE,KNDEX,FOUT,II
    REAL(PDB), INTENT(IN)                  :: X1,Y1,X2,Y2
    REAL(PDB), INTENT(INOUT), DIMENSION(:) :: RPAR
    !****

    ! local variable
    ! **************
    REAL(PDB) :: CX,CY,DELX,DELY,R,THETA,SX,SY,VX,VY,V2,DIST,MOV1,MOV2
    INTEGER   :: I
    REAL(PDB) :: RPAR_OLD(NRPAR)
    !****
    !     keep old parameters:
    RPAR_OLD(1:NRPAR)=RPAR(1:NRPAR)
    !
    MOV1=0._PDB; MOV2=0._PDB;
    SELECT CASE (TYPE)
       CASE (G_ELE_TYPE(1))
       !*       segment :
       !        compute the motions of two ends
       DELX=X1-RPAR(1)
       DELY=Y1-RPAR(2)
       MOV1=SQRT(DELX*DELX+DELY*DELY)
       DELX=X2-(RPAR(1)+RPAR(3))
       DELY=Y2-(RPAR(2)+RPAR(4))
       MOV2=SQRT(DELX*DELX+DELY*DELY)
       !        set new ends, compute direction
       RPAR(1)=X1
       RPAR(2)=Y1
       RPAR(3)=X2-X1
       RPAR(4)=Y2-Y1
       RPAR(5)=SQRT(RPAR(3)*RPAR(3)+RPAR(4)*RPAR(4))
       CASE(G_ELE_TYPE(3))
       !*       arc of circle :
       !        get old radius
       CX=RPAR(1); CY=RPAR(2)
       !        compute the motion of two ends
       R=RPAR(3)
       THETA=RPAR(4)
       DELX=CX+R*COS(THETA)-X1
       DELY=CY+R*SIN(THETA)-Y1
       MOV1=SQRT(DELX*DELX+DELY*DELY)
       THETA=RPAR(5)
       DELX=CX+R*COS(THETA)-X2
       DELY=CY+R*SIN(THETA)-Y2
       MOV2=SQRT(DELX*DELX+DELY*DELY)
       !        compute new radius, new center and new angles
       SX=X1-CX
       SY=Y1-CY
       !        get demi-vector from (x1,y1) to (x2,y2)
       VX=(X2-X1)/2.0
       VY=(Y2-Y1)/2.0
       V2=VX*VX+VY*VY
       !        get center closest to old center, change radius
       !        compute new r
       R=SQRT(V2+(SX*VY-SY*VX)**2/V2)
       RPAR(3)=R
       !        compute new center
       DIST=1.0_PDB+(VX*SX+VY*SY)/V2
       !        get new center and compute angles
       CX=CX+VX*DIST
       CY=CY+VY*DIST
       RPAR(1)=CX
       RPAR(2)=CY
       THETA=SALACO((X1-CX)/R,Y1-CY)
       RPAR(4)=THETA
       THETA=SALACO((X2-CX)/R,Y2-CY)
       RPAR(5)=THETA
       !        if phi1 > phi2 put phi2 equal to phi2+twopi
       IF(RPAR(4)>RPAR(5)) RPAR(5)=RPAR(5)+TWOPI
    END SELECT
    IF(KNDEX/=0) THEN
       IF(MOV1>EPS.OR.MOV2>EPS) THEN
          WRITE(FOUT,'(/," old element ",I4," : ",6(1P,E13.4))') &
               II,(RPAR_OLD(I),I=1,G_ELE_LEN(TYPE))
          WRITE(FOUT,'(" new element ",I4," : ",6(1P,E13.4))') &
               II,(RPAR(I),I=1,G_ELE_LEN(TYPE))
       ENDIF
    ENDIF
    !
  END SUBROUTINE SAL129
  !
  SUBROUTINE SAL130(GG)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! perform domain definition and set boundary conditions
    !   - computes node perimeters per 2d macro
    !   - reads boundary conditions
    !   - computes external perimeters 
    !
    !Parameters: input/output
    ! GG    geometry descriptor
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES,      ONLY : G_BC_TYPE,TYPGEO
    !****
    IMPLICIT NONE
    ! in variable
    ! ************
    TYPE(T_G_BASIC), INTENT(INOUT):: GG
    ! local variable
    !****
    INTEGER :: NN,OK
    INTEGER, DIMENSION(GG%NB_ELEM*2) :: AUX_ARR
    !****
    !*    compute node perimeters for the macro
    CALL SAL130_2(GG%NB_ELEM,GG%NB_NODE,GG%IPAR,GG%PPERIM_NODE, &
         GG%PERIM_NODE,AUX_ARR)
    !
    !*    - compute number of bc's per 2D macro,
    !       NB_BC2 counts total nber of 2D bc's
    !     - compute IBC2_ELEM, keep relative 2D bc nber to elements
    !     - allocation : 2D bc structures
    !                    2D perimeter structure for a macro
    !     - get list of elements in 2d macro perimeter

    CALL SAL130_4(GG%NB_ELEM,NN,GG%IPAR,GG%IBC2_ELEM,AUX_ARR)
    
    GG%NB_BC2=NN

    ALLOCATE (GG%TYPE_BC2(NN),GG%IDATA_BC2(NN),GG%PERIM_MAC2(NN), &
         GG%PPERIM_MAC2(7),STAT=OK)
    IF(OK/=0) CALL XABORT('SAL130: not enough memory I,R')
    GG%PERIM_MAC2(1:NN)=AUX_ARR(1:NN)
    GG%NPERIM_MAC2=NN
    !
    !*    read boundary conditions:
    CALL SAL131(GG)
    !
    !     - define IELEM_SURF2
    IF(GG%NB_SURF2>0) THEN
       ! allocate
       ALLOCATE(GG%SURF2(GG%NB_SURF2),STAT = OK)
       IF(OK /= 0) CALL XABORT('SAL130: not enough memory I,R')
       CALL SAL130_6(GG%NB_SURF2,GG%IBC2_SURF2,GG%PERIM_MAC2,GG%IELEM_SURF2)
    ENDIF
    !
    !*    construct perimeter structures for rotative or symmetrical geometry
    SELECT CASE(TYPGEO)
       CASE(1:2)
       CALL SAL130_8(GG%NPERIM_MAC2,GG%PERIM_MAC2,GG%PPERIM_MAC2, &
            GG%DIST_AXIS,GG%IBC2_ELEM,GG%TYPE_BC2,GG%IDATA_BC2,GG%BCDATA,GG%IPAR,GG%RPAR)
       CASE(5:9)
       CALL SAL130_10(GG%NPERIM_MAC2,GG%PERIM_MAC2,GG%PPERIM_MAC2, &
            GG%DIST_AXIS,GG%IBC2_ELEM,GG%TYPE_BC2,GG%IDATA_BC2,GG%BCDATA,GG%IPAR,GG%RPAR)
    END SELECT
    !
  END SUBROUTINE SAL130
  !
  SUBROUTINE SAL130_2(NB_ELEM,NB_NODE,IPAR,PPERIM,PERIM,LIST)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! compute node perimeters for one 2D macro
    !
    !Parameters: input
    ! NB_ELEM   number of elements
    ! NB_NODE   number of nodes
    ! IPAR      integer descriptors for elements
    !
    !Parameters: output
    ! PPERIM    array pointer to list of elements in perimeter per node
    ! PERIM     elements in perimeters per node
    ! LIST      temporary array
    !
    !---------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN)                    :: NB_ELEM,NB_NODE
    INTEGER, INTENT(IN),   DIMENSION(:,:)  :: IPAR
    INTEGER, INTENT(OUT),  DIMENSION(:)    :: LIST
    INTEGER, INTENT(OUT),  DIMENSION(:)    :: PPERIM
    INTEGER, POINTER,      DIMENSION(:)    :: PERIM
    !****
    INTEGER :: NT,NN,NODE,ELEM,OK
    !****
    NT=0
    PPERIM(1)=1
    DO NODE=1,NB_NODE
       DO ELEM=1, NB_ELEM
          IF(IPAR(2,ELEM)==NODE.OR.IPAR(3,ELEM)==NODE) THEN
             NT=NT+1
             LIST(NT)=ELEM
          ENDIF
       ENDDO
       NN=NT+1-PPERIM(NODE)
       IF(NN>0) THEN
          PPERIM(NODE+1)=NT+1
       ELSE
          CALL XABORT('SAL130_2: node without perimeter')
       ENDIF
    ENDDO
    IF(NT>0) THEN
       ALLOCATE (PERIM(NT), STAT=OK)
       IF(OK/=0) CALL XABORT('SAL130_2: not enough memory I')
       PERIM(1:NT)=LIST(1:NT)
    ENDIF
    !
  END SUBROUTINE SAL130_2
  !
  SUBROUTINE SAL130_4(NB_ELEM,NB_BC,IPAR,IBC2_ELEM,LIST_BC)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! count number of bc's in a 2D macro
    !
    !Parameters: input
    ! NB_ELEM   number of elements
    !
    !Parameters: input/output
    ! IPAR      integer descriptors for elements
    !
    !Parameters: output
    ! NB_BC     nber of bc's
    ! IBC2_ELEM relative 2D bc index per element
    ! LIST_BC   list of elements in boundary
    !
    !---------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN)                    :: NB_ELEM
    INTEGER, INTENT(OUT)                   :: NB_BC
    INTEGER, INTENT(INOUT), DIMENSION(:,:) :: IPAR
    INTEGER, INTENT(OUT),  DIMENSION(:)    :: IBC2_ELEM,LIST_BC
    !****
    INTEGER :: ELEM
    !****
    !     initiation
    IBC2_ELEM=0
    NB_BC=0
    DO ELEM=1, NB_ELEM
       IF(IPAR(2,ELEM)<=0.AND.IPAR(3,ELEM)<=0) THEN
          CALL XABORT('SAL130_4: two boundaries for element')
       ELSE IF(IPAR(2,ELEM)<=0.OR.IPAR(3,ELEM)<=0) THEN
          NB_BC=NB_BC+1
          LIST_BC(NB_BC)=ELEM
          IF(IBC2_ELEM(ELEM)/=0) THEN
             CALL XABORT('SAL130_4: two surfaces to element')
          ENDIF
          IBC2_ELEM(ELEM)=NB_BC
          !           set mark for the macro connection surface :
       ENDIF
    ENDDO
    !
  END SUBROUTINE SAL130_4
  !
  SUBROUTINE SAL130_6(NSURF,IBC2_SURF2,IELEM_BC2,IELEM_SURF2)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! get element order indices per surface
    !
    !Parameters: input
    ! NSURF       number of surfaces
    ! IBC2_SURF2  2D bc order index per surface
    ! IELEM_BC2   element order index per bc
    !
    !Parameters: output
    ! IELEM_SURF2 element order index per surface
    !
    !---------------------------------------------------------------------
    !
    !****
    USE SAL_GEOMETRY_TYPES, ONLY : G_BC_TYPE
    !****
    IMPLICIT NONE
    INTEGER, INTENT(IN)                    :: NSURF
    INTEGER, INTENT(IN),   DIMENSION(:)    :: IBC2_SURF2,IELEM_BC2
    INTEGER, INTENT(OUT),  DIMENSION(:)    :: IELEM_SURF2
    !****
    INTEGER :: SURF,IBC,ELEM
    !****
    IF(NSURF > 0) THEN
       DO SURF=1, NSURF
          !           get relative bc order number
          IBC=IBC2_SURF2(SURF)
          !           get relative element
          ELEM=IELEM_BC2(IBC)
          !           define ielem_surf2
          IELEM_SURF2(SURF)=ELEM
       ENDDO
    ENDIF
    !
  END SUBROUTINE SAL130_6
  !
  SUBROUTINE SAL130_8(NPERIM,PERIM,PPERIM,DIST_AXIS,IBC2_ELEM,TYPE_BC2,IDATA_BC2, &
       BCDATA,IPAR,RPAR)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! calculate PERIM_MAC2,DIST_AXIS for rotative or symmetrical geometry
    ! (TYPGEO=1 & 2)
    !
    !Parameters: input
    ! NPERIM     number of elements in perimeter
    ! IBC2_ELEM  2D bc order number per element
    ! TYPE_BC2   type of boundary conditions per 2D bc
    ! IDATA_BC2  position in the 'bcdata' table for each 2D bc
    ! BCDATA     table of boundary conditions data
    ! IPAR       integer element descriptor table
    ! RPAR       floating point element descriptor table
    !
    !Parameters: input/output
    ! PERIM      list of elements in perimeter,in return it will be in following
    !            order: (elems on axis 1)+(elems on axis 2)+(other elems)
    !
    !Parameters: output
    ! PPERIM     pointers to the table 'perim':
    !            (1): first elem on axis 1 (2): first elem on axis 2
    !            (3): first of other elems (4): NPERIM + 1
    ! DIST_AXIS  distance of points on axis to the center (0,0),in order of:
    !            (distances of points on axis 1)+(distances of points on axis 2)
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES,     ONLY : ANGGEO,EPS,G_BC_TYPE
    USE SAL_NUMERIC_MOD, ONLY : SAL141
    !****
    IMPLICIT NONE
    INTEGER, INTENT(IN)                  :: NPERIM
    INTEGER, INTENT(IN),    DIMENSION(:) :: IBC2_ELEM,TYPE_BC2,IDATA_BC2
    INTEGER, INTENT(INOUT), DIMENSION(:) :: PERIM
    INTEGER, INTENT(OUT),   DIMENSION(:) :: PPERIM
    REAL(PDB), POINTER,     DIMENSION(:) :: DIST_AXIS
    INTEGER, INTENT(IN),  DIMENSION(:,:) :: IPAR
    REAL(PDB), INTENT(INOUT),DIMENSION(:,:) :: RPAR
    REAL(PDB), INTENT(IN),DIMENSION(:,:) :: BCDATA
    !****
    !     LIST_ELEMS = table of elements,elements on axis are in increasing order
    !                  of distance to (0,0):
    !                  1=elements on axis 1; 2=elements on axis 2; 3=other elements
    !     AUX_DIST   = max distance of element ends to the beginnings of the axes:
    !                  1=distances on axis 1; 2=distances on axis 2
    INTEGER,   DIMENSION(NPERIM,3) :: LIST_ELEMS
    REAL(PDB), DIMENSION(NPERIM,2) :: AUX_DIST
    INTEGER   :: I,J,K,M,ELEM,TYPBC,IBC,IDATA,NBE(3),OK,NAXES
    REAL(PDB) :: ANGLE,X,Y,D
    !****
    NAXES=2
    !*    calculate number of elements on axis 1 & 2
    !     and their distances to the beginning of the axes
    NBE=0
    DO I=1,NPERIM
       ELEM=PERIM(I)
       IBC=IBC2_ELEM(ELEM)
       TYPBC=TYPE_BC2(IBC)
       IDATA=IDATA_BC2(IBC)
       !        default is the elements not on axis
       M=NAXES+1
       IF(TYPBC==G_BC_TYPE(3)) THEN
          !*          rotation:
          ANGLE=BCDATA(5,IDATA)
          IF(ABS(ANGLE-ANGGEO)<EPS) THEN
             !              element is on the axis 1 (angle=anggeo)
             M=1
          ELSEIF(ABS(ANGLE+ANGGEO)<EPS) THEN
             !              element is on the axis 2 (angle=-anggeo)
             M=2
          ENDIF
       ELSEIF(TYPBC==G_BC_TYPE(4)) THEN
          !*          symmetry:
          ANGLE=BCDATA(5,IDATA)
          IF(ABS(ANGLE)<EPS) THEN
             !              element is on the axis 1 (angle=0)
             M=1
          ELSEIF(ABS(ANGLE-ANGGEO)<EPS) THEN
             !              element is on the axis 2 (angle=anggeo)
             M=2
          ENDIF
       ENDIF
       NBE(M)=NBE(M)+1
       LIST_ELEMS(NBE(M),M)=ELEM
       IF(M/=NAXES+1) THEN
          !           sort the element list according their distance to (0,0)
          D=0.
          DO K=1,2
             CALL SAL141(IPAR(1,ELEM),RPAR(:,ELEM),X,Y,K)
             D=MAX(D,SQRT(X*X+Y*Y))
          ENDDO
          AUX_DIST(NBE(M),M)=D
          IF(NBE(M)>1) THEN
             DO K=1,NBE(M)-1
                IF(AUX_DIST(NBE(M),M)<AUX_DIST(K,M)) THEN
                   !                    insert the last element to position k
                   D=AUX_DIST(NBE(M),M)
                   ELEM=LIST_ELEMS(NBE(M),M)
                   DO J=NBE(M),K+1,-1
                      AUX_DIST(J,M)=AUX_DIST(J-1,M)
                      LIST_ELEMS(J,M)=LIST_ELEMS(J-1,M)
                   ENDDO
                   AUX_DIST(K,M)=D
                   LIST_ELEMS(K,M)=ELEM
                   EXIT
                ENDIF
             ENDDO
          ENDIF
       ENDIF
    ENDDO
    !     nber of blocks in 'perim'
    M=NAXES+1
    !*    organize elements as:
    !       (elems on axis 1)+(elems on axis 2)+(other elems)
    PPERIM(:)=0
    PPERIM(1)=1
    DO I=1,M
       PPERIM(I+1)=PPERIM(I)+NBE(I)
    ENDDO
    DO I=1,M
       PERIM(PPERIM(I):PPERIM(I+1)-1)=LIST_ELEMS(1:NBE(I),I)
    ENDDO
    !*    allocate and set dist_axis
    IF(PPERIM(NAXES+1)-1>0) THEN
       ALLOCATE(DIST_AXIS(PPERIM(NAXES+1)-1),STAT=OK)
       IF(OK.NE.0) &
            CALL XABORT('SAL130_8: not enough memory R')
       DO I=1,NAXES
          DIST_AXIS(PPERIM(I):PPERIM(I+1)-1)=AUX_DIST(1:NBE(I),I)
       ENDDO
    ELSE
       NULLIFY(DIST_AXIS)
    ENDIF
    !
  end subroutine sal130_8
  !
  SUBROUTINE SAL130_10(NPERIM,PERIM,PPERIM,DIST_AXIS,IBC2_ELEM,TYPE_BC2,IDATA_BC2, &
       BCDATA,IPAR,RPAR)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! calculate perim_mac2,dist_axis for the cyclical geometries:
    ! type 5&6: retangle with translations or symmetries on all sides
    ! type 7  : 1/8 assembly with symmetries on all sides
    ! type 8  : equilateral triangle with symmetries on all sides
    ! type 9  : hexagon with symmetries on all sides
    !
    !Parameters: input
    ! NPERIM     number of elements in perimeter
    ! IBC2_ELEM  2D bc order number per element
    ! TYPE_BC2   type of boundary conditions per 2D bc
    ! IDATA_BC2  position in the 'bcdata' table for each 2D bc
    ! BCDATA     table of boundary conditions data
    ! IPAR       integer element descriptor table
    ! RPAR       floating point element descriptor table
    !
    !Parameters: input/output
    ! PERIM      list of elements in perimeter,in return it will be in following
    !            order: (elems on axis 1)+(elems on axis 2)+(other elems)
    !
    !Parameters: output
    ! PPERIM     pointers to the table 'perim':
    ! DIST_AXIS  distance of points on axis to the axis origin
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES,     ONLY : TYPGEO,ANGGEO,EPS,LX=>LENGTHX,LY=>LENGTHY,G_BC_TYPE
    USE SAL_NUMERIC_MOD, ONLY : SAL141
    !****
    IMPLICIT NONE
    INTEGER, INTENT(IN)                  :: NPERIM
    INTEGER, INTENT(IN),    DIMENSION(:) :: IBC2_ELEM,TYPE_BC2,IDATA_BC2
    INTEGER, INTENT(INOUT), DIMENSION(:) :: PERIM
    INTEGER, INTENT(OUT),   DIMENSION(:) :: PPERIM
    REAL(PDB), POINTER,     DIMENSION(:) :: DIST_AXIS
    INTEGER, INTENT(IN),  DIMENSION(:,:) :: IPAR
    REAL(PDB), INTENT(INOUT),DIMENSION(:,:) :: RPAR
    REAL(PDB), INTENT(IN),DIMENSION(:,:) :: BCDATA
    !****
    !     LIST_ELEMS = table of elements,elements on axis are in increasing order
    !                  of distance to axis origin:
    !                  1=elements on axis 1; 2=elements on axis 2; ...
    !     AUX_DIST   = max distance of element ends to axis origin:
    !                  1=distances on axis 1; 2=distances on axis 2; ...
    INTEGER,   DIMENSION(NPERIM,6) :: LIST_ELEMS
    REAL(PDB), DIMENSION(NPERIM,6) :: AUX_DIST
    INTEGER   :: I,J,K,M,ELEM,TYPBC,IBC,IDATA,NBE(6),OK,NAXES
    REAL(PDB) :: ANGLE,X,Y,D
    !****
    NAXES=0
    IF(TYPGEO==5.OR.TYPGEO==6) THEN
       NAXES=4
    ELSEIF(TYPGEO==7) THEN
       NAXES=3
    ELSEIF(TYPGEO==8) THEN
       NAXES=3
    ELSEIF(TYPGEO==9) THEN
       NAXES=6
    ENDIF
    !*    calculate number of elements on axes,
    !     and their distances to the origin of the axis
    NBE=0
    DO I=1,NPERIM
       ELEM=PERIM(I)
       IBC=IBC2_ELEM(ELEM)
       TYPBC=TYPE_BC2(IBC)
       IDATA=IDATA_BC2(IBC)
       !        default is the elements not on axis
       M=NAXES+1
       SELECT CASE (TYPGEO)
          CASE (5)
          !*          rectangle with translations:
          !           axes definition:
          !                     axis 3
          !                  ************
          !                  *          *
          !           axis 2 *          * axis 4
          !                  ************
          !                     axis 1
          IF(TYPBC==G_BC_TYPE(2)) THEN
             IF(BCDATA(2,IDATA)>0) THEN
                M=1
             ELSEIF(BCDATA(2,IDATA)<0) THEN
                M=3
             ELSEIF(BCDATA(1,IDATA)>0) THEN
                M=2
             ELSEIF(BCDATA(1,IDATA)<0) THEN
                M=4
             ENDIF
          ELSE
             CALL XABORT('SAL130_10: wrong boundary condition for TYPGEO=5')
          ENDIF
          CASE (6)
          !*          rectangle with symmetries:
          !           axes definition:
          !                     axis 3
          !                  ************
          !                  *          *
          !           axis 2 *          * axis 4
          !                  ************
          !                     axis 1
          IF(TYPBC==G_BC_TYPE(4)) THEN
             ANGLE=BCDATA(5,IDATA)
             IF(ABS(ANGLE)<EPS) THEN
                IF(BCDATA(2,IDATA)>0) THEN
                   !                    cy>0: element is on the axis 3 (angle=0)
                   M=3
                ELSE
                   !                    cy=0: element is on the axis 1 (angle=0)
                   M=1
                ENDIF
             ELSEIF(ABS(ANGLE-ANGGEO)<EPS) THEN
                IF(BCDATA(1,IDATA)>0) THEN
                   !                    cx>0: element is on the axis 4 (angle=anggeo)
                   M=4
                ELSE
                   !                    cx=0: element is on the axis 2 (angle=anggeo)
                   M=2
                ENDIF
             ENDIF
          ELSE
             CALL XABORT('SAL130_10: wrong boundary condition for TYPGEO=6')
          ENDIF
          CASE (7:8)
          !*          triangles with symmetries:
          !           typgeo=7 axes definition:
          !                        *
          !               axis 2 * * axis 3
          !                    *   *
          !                  *******
          !                   axis 1
          !           typgeo=8 axes definition:
          !                        *
          !              axis 2  *   *  axis 3
          !                    *       *
          !                  *************
          !                      axis 1
          IF(TYPBC==G_BC_TYPE(4)) THEN
             ANGLE=BCDATA(5,IDATA)
             IF(ABS(ANGLE)<EPS) THEN
                !                 element is on the axis 1 (angle=0)
                M=1
             ELSEIF(ABS(ANGLE-ANGGEO)<EPS) THEN
                !                 element is on the axis 2 (angle=anggeo)
                M=2
             ELSEIF(ABS(ANGLE-HALFPI)<EPS) THEN
                !                 typgeo=7:element is on the axis 3 (angle=pi/2)
                M=3
             ELSEIF(ABS(ANGLE-2.*ANGGEO)<EPS) THEN
                !                 typgeo=8:element is on the axis 3 (angle=2*anggeo)
                M=3
             ELSE
                CALL XABORT('SAL130_10: boundary condition data error in element for TYPGEO=7,8:')
             ENDIF
          ELSE
             CALL XABORT('SAL130_10: wrong boundary condition for TYPGEO=7,8')
          ENDIF
          CASE (9)
          !*          hexagon with translations:
          !           axes definition:
          !                                axis 4 (0,-2ly)
          !                                  *****
          !              axis 3 (3/2lx,-ly) *     * axis 5 (-3/2lx,-ly)
          !                                *       *
          !              axis 2 (3/2lx, ly) *     * axis 6 (-3/2lx, ly)
          !                                  *****
          !                                axis 1 (0, 2ly)
          !           origins of axes:
          !           axis 1&2: (-lx/2,-ly)
          !           axis 3:   (-lx  ,  0)
          !           axis 4:   (-lx/2, ly)
          !           axis 5:   ( lx  ,  0)
          !           axis 6:   ( lx/2,-ly)
          IF(TYPBC==G_BC_TYPE(2)) THEN
             IF(ABS(BCDATA(1,IDATA))<EPS) THEN
                IF(BCDATA(2,IDATA)>0) THEN
                   M=1
                ELSE
                   M=4
                ENDIF
             ELSEIF(BCDATA(1,IDATA)>0.AND.BCDATA(2,IDATA)>0) THEN
                M=2
             ELSEIF(BCDATA(1,IDATA)>0.AND.BCDATA(2,IDATA)<0) THEN
                M=3
             ELSEIF(BCDATA(1,IDATA)<0.AND.BCDATA(2,IDATA)>0) THEN
                M=6
             ELSEIF(BCDATA(1,IDATA)<0.AND.BCDATA(2,IDATA)<0) THEN
                M=5
             ENDIF
          ELSE
             CALL XABORT('SAL130_10: wrong boundary condition for TYPGEO=9')
          ENDIF
       END SELECT
       IF(M==NAXES+1) CALL XABORT('SAL130_10: element not on axes')
       NBE(M)=NBE(M)+1
       LIST_ELEMS(NBE(M),M)=ELEM
       !        sort the element list according their distance to
       !        the origins of the axes
       D=0.
       SELECT CASE (TYPGEO)
          CASE (5:6)
          SELECT CASE (M)
             CASE (1:2)
             DO K=1,2
                CALL SAL141(IPAR(1,ELEM),RPAR(:,ELEM),X,Y,K)
                D=MAX(D,SQRT(X*X+Y*Y))
             ENDDO
             CASE (3)
             DO K=1,2
                CALL SAL141(IPAR(1,ELEM),RPAR(:,ELEM),X,Y,K)
                D=MAX(D,X)
             ENDDO
             CASE (4)
             DO K=1,2
                CALL SAL141(IPAR(1,ELEM),RPAR(:,ELEM),X,Y,K)
                D=MAX(D,Y)
             ENDDO
          END SELECT
          CASE (7)
          SELECT CASE (M)
             CASE (1:2)
             DO K=1,2
                CALL SAL141(IPAR(1,ELEM),RPAR(:,ELEM),X,Y,K)
                D=MAX(D,SQRT(X*X+Y*Y))
             ENDDO
             CASE (3)
             DO K=1,2
                CALL SAL141(IPAR(1,ELEM),RPAR(:,ELEM),X,Y,K)
                D=MAX(D,Y)
             ENDDO
          END SELECT
          CASE (8)
          SELECT CASE (M)
             CASE (1:2)
             DO K=1,2
                CALL SAL141(IPAR(1,ELEM),RPAR(:,ELEM),X,Y,K)
                D=MAX(D,SQRT(X*X+Y*Y))
             ENDDO
             CASE (3)
             DO K=1,2
                CALL SAL141(IPAR(1,ELEM),RPAR(:,ELEM),X,Y,K)
                !                 axis origin is (lx,0)
                D=MAX(D,SQRT((X-LX)*(X-LX)+Y*Y))
             ENDDO
          END SELECT
          CASE (9)
          !           origins of axes:
          !           axis 1&2: (-lx/2,-ly)
          !           axis 3:   (-lx  ,  0)
          !           axis 4:   (-lx/2, ly)
          !           axis 5:   ( lx  ,  0)
          !           axis 6:   ( lx/2,-ly)
          SELECT CASE (M)
             CASE (1:2)
             DO K=1,2
                CALL SAL141(IPAR(1,ELEM),RPAR(:,ELEM),X,Y,K)
                D=MAX(D,SQRT((X+LX*0.5)*(X+LX*0.5)+(Y+LY)*(Y+LY)))
             ENDDO
             CASE (3)
             DO K=1,2
                CALL SAL141(IPAR(1,ELEM),RPAR(:,ELEM),X,Y,K)
                D=MAX(D,SQRT((X+LX)*(X+LX)+Y*Y))
             ENDDO
             CASE (4)
             DO K=1,2
                CALL SAL141(IPAR(1,ELEM),RPAR(:,ELEM),X,Y,K)
                D=MAX(D,SQRT((X+LX*0.5)*(X+LX*0.5)+(Y-LY)*(Y-LY)))
             ENDDO
             CASE (5)
             DO K=1,2
                CALL SAL141(IPAR(1,ELEM),RPAR(:,ELEM),X,Y,K)
                D=MAX(D,SQRT((X-LX)*(X-LX)+Y*Y))
             ENDDO
             CASE (6)
             DO K=1,2
                CALL SAL141(IPAR(1,ELEM),RPAR(:,ELEM),X,Y,K)
                D=MAX(D,SQRT((X-LX*0.5)*(X-LX*0.5)+(Y+LY)*(Y+LY)))
             ENDDO
          END SELECT
       END SELECT
       AUX_DIST(NBE(M),M)=D
       IF(NBE(M)>1) THEN
          DO K=1,NBE(M)-1
             IF(AUX_DIST(NBE(M),M)<AUX_DIST(K,M)) THEN
                !                 insert the last element to position k
                D=AUX_DIST(NBE(M),M)
                ELEM=LIST_ELEMS(NBE(M),M)
                DO J=NBE(M),K+1,-1
                   AUX_DIST(J,M)=AUX_DIST(J-1,M)
                   LIST_ELEMS(J,M)=LIST_ELEMS(J-1,M)
                ENDDO
                AUX_DIST(K,M)=D
                LIST_ELEMS(K,M)=ELEM
                EXIT
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    !     nber of blocks in 'perim'
    M=NAXES
    !*    organize elements as:
    !       (elems on axis 1)+(elems on axis 2)+...
    PPERIM(1)=1
    DO I=1,M
       PPERIM(I+1)=PPERIM(I)+NBE(I)
    ENDDO
    DO I=1,M
       PERIM(PPERIM(I):PPERIM(I+1)-1)=LIST_ELEMS(1:NBE(I),I)
    ENDDO
    !*    allocate and set dist_axis
    IF(PPERIM(NAXES+1)-1>0) THEN
       ALLOCATE(DIST_AXIS(PPERIM(NAXES+1)-1),STAT=OK)
       IF(OK.NE.0) CALL XABORT('SAL130_10: not enough memory R')
       DO I=1,NAXES
          DIST_AXIS(PPERIM(I):PPERIM(I+1)-1)=AUX_DIST(1:NBE(I),I)
       ENDDO
    ELSE
       NULLIFY(DIST_AXIS)
    ENDIF
    !
  END SUBROUTINE SAL130_10
  !
  SUBROUTINE SAL131(GG)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! modifies GG%IPAR and constructs bcdata (modifies GG%IPAR and
    ! constructs BCDATA)
    !
    !Parameters: input/output
    ! GG    geometry descriptor
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES,      ONLY : G_BC_MAX_LEN,G_BC_LEN,G_BC_TYPE,NIPAR,EPS, &
                                    ANGGEO,TYPGEO,LENGTHX,LENGTHY,ALLSUR
    !****
    IMPLICIT NONE
    ! in variable
    ! ************
    TYPE(T_G_BASIC), INTENT(INOUT) :: GG
    !****
    LOGICAL   :: LGBC,LGALLS
    INTEGER   :: TYPE,II,NN,NBER,I,J,ITBC,IELEM
    REAL(PDB) :: ANGLE
    INTEGER   :: ELEM,IBC,OK,IDATA(2)
    INTEGER, POINTER, DIMENSION(:) :: TYPE_BC2,ISURF2_ELEM,IELEM_BC2
    !*****
    !**   TMP_BCDATA        = description of motions at the boundary
    !**   IPAR              = pointer to geometry descriptors
    INTEGER, DIMENSION(GG%NB_ELEM) :: AUX_ARR
    REAL(PDB), DIMENSION(G_BC_MAX_LEN,GG%NB_ELEM+3) :: TMP_BCDATA
    !*****
    !*    initialization
    IF(TYPGEO==1.OR.TYPGEO==2) IDATA(:)=0
    !*    BCDATA for surfaces of type G_BC_TYPE(-1)
    !
    LGALLS=ALLSUR/=0
    ITBC=1      ! the first bc data
    TMP_BCDATA(:,:) = 0._PDB
    TMP_BCDATA(1,ITBC)=1._PDB
    ITBC=ITBC+1
    NN=1
    TMP_BCDATA(NN,ITBC)=GG%ALBEDO
    !*    treat approximate boundary condictions
    IF(LGALLS)THEN
       !        ALL BC'S PRODUIT SURFACES
       IF(GG%DEFAUL==1)THEN
          !           SPECULAR REFLEXION -> ISOTROPIC REFLEXION WITH ALBEDO=1
          GG%DEFAUL=0 ; TMP_BCDATA(NN,ITBC)=1.
       ENDIF
    ENDIF
    !
    !*    put default value in all bc elements:
    CALL SAL131_2(GG%NB_ELEM,GG%DEFAUL,GG%IPAR,GG%IBC2_ELEM,GG%TYPE_BC2,GG%IDATA_BC2)
    !
    IF(GG%NBBCDA>0)THEN
       DO I=1,GG%NBBCDA
          ITBC=ITBC+1
          TMP_BCDATA(:,ITBC)=GG%BCDATAREAD(I)%BCDATA(:)
          TYPE=GG%BCDATAREAD(I)%SALTYPE
          SELECT CASE (TYPE)
             CASE(1)
             IF(LGALLS)THEN
                !                 specular reflexion -> isotropic reflexion with albedo=1
                TYPE=0 ; TMP_BCDATA(1,ITBC)=1.
             ENDIF
             CASE(2)
             SELECT CASE (TYPGEO)
                CASE (5) 
                !                 adjust translation data according to the sides of rectangle
                IF(TMP_BCDATA(1,ITBC)/=0) &
                     TMP_BCDATA(1,ITBC)=SIGN(LENGTHX,TMP_BCDATA(1,ITBC))
                IF(TMP_BCDATA(2,ITBC)/=0) &
                     TMP_BCDATA(2,ITBC)=SIGN(LENGTHY,TMP_BCDATA(2,ITBC))
                CASE(9)
                !                adjust translation data according to the hexagonal sides
                IF(ABS(TMP_BCDATA(1,ITBC))<EPS) THEN
                   !                    axes 1 & 4
                   TMP_BCDATA(2,ITBC)=SIGN(2.*LENGTHY,TMP_BCDATA(2,ITBC))
                ELSE
                   TMP_BCDATA(1,ITBC)=SIGN(LENGTHX*1.5,TMP_BCDATA(1,ITBC))
                   TMP_BCDATA(2,ITBC)=SIGN(LENGTHY,TMP_BCDATA(2,ITBC))
                ENDIF
             END SELECT
             IF(LGALLS)THEN
                !                 translation -> isotropic translation with albedo=1
                TYPE=TYPE+4
                !                 perform array shift
                TMP_BCDATA(:,ITBC)=CSHIFT(TMP_BCDATA(:,ITBC),1)
                !                 albedo=1
                TMP_BCDATA(1,ITBC)=1._PDB
             ENDIF
             CASE(3) 
             !              cases of rotation:
             !              read angle, compute cos and sin
             ANGLE=TMP_BCDATA(3,ITBC)*(TWOPI/360._PDB)
             IF(TYPGEO==1.OR.TYPGEO==2) THEN
                IF(ABS(ANGLE-ANGGEO)<EPS) THEN
                   !                    for the axis 1:keep anggeo,cos(anggeo),sin(anggeo)
                   ANGLE=ANGGEO
                   IF(IDATA(1)==0) IDATA(1)=ITBC
                ELSEIF(ABS(ANGLE+ANGGEO)<EPS) THEN
                   !                    for the axis 2:keep -anggeo,cos(-anggeo),sin(-anggeo)
                   ANGLE=-ANGGEO
                   IF(IDATA(2)==0) IDATA(2)=ITBC
                ELSE
                   CALL XABORT('SAL131: error in angle of rotative axis')
                ENDIF
             ENDIF
             TMP_BCDATA(3,ITBC)=COS(ANGLE)
             TMP_BCDATA(4,ITBC)=SIN(ANGLE)
             TMP_BCDATA(5,ITBC)=ANGLE
             IF(LGALLS)THEN
                !                 rotation -> isotropic rotation with albedo=1
                !                 axial symmetry -> isotropic axial symmetry with albedo=1
                TYPE=TYPE+4
                !                 perform array shift
                TMP_BCDATA(:,ITBC)=CSHIFT(TMP_BCDATA(:,ITBC),1)
                !                 albedo=1
                TMP_BCDATA(1,ITBC)=1._PDB
             ENDIF
             CASE(4)
             !              cases of symmetry:
             !              read center and angle, compute cos and sin
             ANGLE=TMP_BCDATA(3,ITBC)*(TWOPI/360._PDB)
             SELECT CASE (TYPGEO)
                CASE(1,2)
                IF(ABS(ANGLE-0.)<EPS) THEN
                   ANGLE=0.
                   IF(IDATA(1)==0) IDATA(1)=ITBC
                ELSEIF(ABS(ANGLE-ANGGEO)<EPS) THEN
                   ANGLE=ANGGEO
                   IF(IDATA(2)==0) IDATA(2)=ITBC
                ELSE
                   WRITE(*,*) ' itbc=',ITBC,' angle=',ANGLE,' anggeo=',ANGGEO
                   CALL XABORT('SAL131: error in angle of symmetry axis')
                ENDIF
                CASE(6:8)
                !                 adjust translation data according to the rectangle/triangle sides
                IF(TMP_BCDATA(1,ITBC)/=0) &
                     TMP_BCDATA(1,ITBC)=SIGN(LENGTHX,TMP_BCDATA(1,ITBC))
                IF(TMP_BCDATA(2,ITBC)/=0) &
                     TMP_BCDATA(2,ITBC)=SIGN(LENGTHY,TMP_BCDATA(2,ITBC))
             END SELECT
             TMP_BCDATA(3,ITBC)=COS(ANGLE)
             TMP_BCDATA(4,ITBC)=SIN(ANGLE)
             TMP_BCDATA(5,ITBC)=ANGLE
             IF(LGALLS)THEN
                !                 rotation -> isotropic rotation with albedo=1
                !                 axial symmetry -> isotropic axial symmetry with albedo=1
                TYPE=TYPE+4
                !                 perform array shift
                TMP_BCDATA(:,ITBC)=CSHIFT(TMP_BCDATA(:,ITBC),1)
                !                 albedo=1
                TMP_BCDATA(1,ITBC)=1._PDB
             ENDIF
          END SELECT
          !
          !           modify notation for bd conditions
          NBER=GG%BCDATAREAD(I)%NBER
          DO J=1,NBER
             ELEM=GG%BCDATAREAD(I)%ELEMNB(J)
             IF(ELEM>GG%NB_ELEM.OR.ELEM<=0) CALL XABORT('SAL131: unknown bc element')
             !        get local surface nber
             IBC=GG%IBC2_ELEM(ELEM)
             LGBC=GG%IPAR(2,ELEM)<=0
             II=0
             IF(LGBC)THEN
                II=2
             ELSE
                LGBC=GG%IPAR(3,ELEM)<=0
                IF(LGBC)II=3
             ENDIF
             IF(.NOT.LGBC) THEN
                WRITE(*,*) 'elem :',ELEM
                WRITE(*,*) 'GG%IPAR(:,ELEM) :',GG%IPAR(:,ELEM)
                CALL XABORT('SAL131: wrong bc element')
             ENDIF
             !              put bc type
             GG%IPAR(II,ELEM)=G_BC_TYPE(TYPE)
             GG%TYPE_BC2(IBC)=G_BC_TYPE(TYPE)
             !              put bc data position :
             GG%IDATA_BC2(IBC)=ITBC
          ENDDO
       ENDDO
    ENDIF
    !
    !*    - set BCDATA position for surfaces of type G_BC_TYPE(-1)
    !     - compute the nber of surfaces (type -1,0,-12,-13,-14,-15) : nbsur2
    !     - allocate structures for the surfaces
    !     - compute surf_mac2
    GG%NB_SURF2=0
    TYPE_BC2=>GG%TYPE_BC2
    ISURF2_ELEM=>GG%ISURF2_ELEM
    IELEM_BC2=>GG%PERIM_MAC2
    NN=0
    DO IBC=1,GG%NB_BC2
       IF(TYPE_BC2(IBC)==G_BC_TYPE(-1)) THEN
          !              macro contact surfaces : set bcdata position to 1
          GG%IDATA_BC2(IBC)=1
       ENDIF
       !           relative element nber
       IELEM=IELEM_BC2(IBC)
       !           count 2D surfaces number
       IF(TYPE_BC2(IBC)==G_BC_TYPE(-1) .OR. &
            TYPE_BC2(IBC)==G_BC_TYPE(0)) THEN
          NN=NN+1
          GG%NB_SURF2=GG%NB_SURF2+1
          AUX_ARR(NN)=IBC
          ISURF2_ELEM(IELEM)=NN
       ELSE
          ISURF2_ELEM(IELEM)=0
       ENDIF
    ENDDO
    GG%NB_SURF2=NN
    IF(NN>0) THEN
       ALLOCATE (GG%IBC2_SURF2(NN),GG%IELEM_SURF2(NN),STAT=OK)
       IF(OK/=0) CALL XABORT('SAL131: NOT ENOUGH MEMORY I,R')
       GG%IBC2_SURF2(1:NN)=AUX_ARR(1:NN)
    ELSE
       NULLIFY(GG%IBC2_SURF2,GG%IELEM_SURF2)
    ENDIF
    !
    !*    allocate idata_axis
    IF(TYPGEO==1.OR.TYPGEO==2) THEN
       DO I=1,2
          IF(IDATA(I)==0) THEN
             !              there is no elements on this axis,add a bcdata for this axis
             ITBC=ITBC+1
             IF(TYPGEO==1) THEN
                !                 symmetry
                IF(I==1) THEN
                   ANGLE=0
                ELSE
                   ANGLE=ANGGEO
                ENDIF
             ELSEIF(TYPGEO==2) THEN
                !                 rotation
                IF(I==1) THEN
                   ANGLE=ANGGEO
                ELSE
                   ANGLE=-ANGGEO
                ENDIF
             ENDIF
             TMP_BCDATA(1,ITBC)=0.
             TMP_BCDATA(2,ITBC)=0.
             TMP_BCDATA(3,ITBC)=COS(ANGLE)
             TMP_BCDATA(4,ITBC)=SIN(ANGLE)
             TMP_BCDATA(5,ITBC)=ANGLE
             IDATA(I)=ITBC
          ENDIF
       ENDDO
       ALLOCATE(GG%IDATA_AXIS(2), STAT=OK)
       IF(OK.NE.0) CALL XABORT('SAL131: not enough memory I')
       GG%IDATA_AXIS=IDATA
    ENDIF
    !*    allocate bcdata
    ALLOCATE (GG%BCDATA(G_BC_MAX_LEN,ITBC), STAT=OK)
    ! MISE A ZERO
    IF(OK/=0) CALL XABORT('SAL131: not enough memory R')
    GG%BCDATA(:,1:ITBC)=TMP_BCDATA(:,1:ITBC)
    GG%NALBG=ITBC
    !
  END SUBROUTINE SAL131
  !
  SUBROUTINE SAL131_2(NB_ELEM,DEFAUL,IPAR,IBC2_ELEM,TYPE_BC2,IDATA_BC2)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! put default bc type and BCDATA to a 2D macro
    !
    !Parameters: input
    ! NB_ELEM    nber of elements
    ! DEFAUL     default bc type
    ! IBC2_ELEM  relative 2D bc number per element
    !
    !Parameters: input/output
    ! IPAR       integer descriptors for elements
    !
    !Parameters: input
    ! TYPE_BC2    2D boundary condiction type
    ! IDATA_BC2   position in bcdata par 2D bc
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES, ONLY : G_BC_TYPE
    !****
    IMPLICIT NONE
    INTEGER, INTENT(IN)                    :: NB_ELEM,DEFAUL
    INTEGER, INTENT(INOUT), DIMENSION(:,:) :: IPAR
    INTEGER, INTENT(IN),   DIMENSION(:)    :: IBC2_ELEM
    INTEGER, INTENT(OUT),  DIMENSION(:)    :: TYPE_BC2
    INTEGER, INTENT(OUT),  DIMENSION(:)    :: IDATA_BC2
    !****
    INTEGER :: ELEM,II,IBC
    LOGICAL :: LGBC
    !****
    !     initiation
    TYPE_BC2=G_BC_TYPE(-1)
    !
    DO ELEM=1,NB_ELEM
       LGBC=IPAR(2,ELEM)==0
       II=0
       IF(LGBC) THEN
          II=2
          IF(IPAR(3,ELEM)<=0) CALL XABORT('SAL131_2: element with 2 bc''s')
       ELSE
          LGBC=IPAR(3,ELEM)==0
          IF(LGBC) II=3
       ENDIF
       IF(LGBC) THEN
          !           put bc type value to ipar
          IPAR(II,ELEM)=G_BC_TYPE(DEFAUL)
          !           put bc type value to bc type structure
          IBC=IBC2_ELEM(ELEM)
          IF(IBC==0) CALL XABORT('SAL131_2: surf-element relation error')
          IF(TYPE_BC2(IBC)/=G_BC_TYPE(-1)) CALL XABORT('SAL131_2: two elements to a surface')
          TYPE_BC2(IBC)=G_BC_TYPE(DEFAUL)
          !           put position of "defaul" in bcdata table (always the second position)
          IDATA_BC2(IBC)=2
       ENDIF
    ENDDO
    !
  END SUBROUTINE SAL131_2
  !
  SUBROUTINE SAL140(NB_NODE,RPAR,IPAR,PPERIM,PERIM)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! checks domain topology in a 2D macro
    !
    !Parameters: input
    ! NB_NODE    number of nodes in macro
    !
    !Parameters: input/output
    ! RPAR       floating point geometry descriptors
    ! IPAR       integer descriptors for elements
    ! PPERIM     pointer to list of elements in perimeter
    ! PERIM      list of perimeters
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES, ONLY : NRPAR,NIPAR
    !****
    IMPLICIT NONE
    INTEGER,   INTENT(IN)   :: NB_NODE
    INTEGER,   INTENT(IN),    DIMENSION (:)   :: PPERIM
    INTEGER,   INTENT(INOUT), DIMENSION (:,:) :: IPAR
    INTEGER,   INTENT(INOUT), DIMENSION (:)   :: PERIM
    REAL(PDB), INTENT(INOUT), DIMENSION (:,:) :: RPAR
    !****
    INTEGER   :: ELEM,NODE,ID,L,FIRST,LAST,NEXT,NLOOP, &
         NEWEND,KEEP,IEND,I1,I2
    REAL(PDB) :: X,Y,XNEW,YNEW,DIST
    LOGICAL   :: LGOPEN,LGON,LGERR
    INTEGER, PARAMETER      :: MXKEEP = 20
    REAL, DIMENSION(MXKEEP) :: KEEPD
    INTEGER, PARAMETER :: FOUT =6
    REAL(PDB),PARAMETER :: EPS2=1.0E-7_PDB
    !*****
    NEXT=0
    !*    checks topology of each node
    DO NODE=1,NB_NODE
       I1=PPERIM(NODE)
       !        ID counts elements in the perimeter already treated
       ID=I1-1
       I2=PPERIM(NODE+1)-1
       IF(I2<I1) CALL XABORT('SAL140: node without perimeter')
       !        mark elements in perimeter of node as untreated
       !        (except full circle: ipar=2)
       DO L=I1,I2
          ELEM=PERIM(L)
          IF(IPAR(1,ELEM)==2)THEN
             ID=ID+1
             IF(ID<L)THEN
                ! move circles at beginning
                PERIM(L)=PERIM(ID)
                PERIM(ID)=ELEM
             ENDIF
          ENDIF
       ENDDO
       ID=ID+1
       !        nloop counts nber of loops for node
       NLOOP=ID-1
       !        treat all elements in perimeter until they have been all done:
       DO WHILE (ID<I2)
          !           an element in the perimeter that has not been done
          FIRST=PERIM(ID)
          NLOOP=NLOOP+1
          !           get first end for this element:
          IEND=1
          !           determine the loop defined by element last:
          LGOPEN=.TRUE.
          LAST=FIRST
          DO WHILE(LGOPEN)
             !              get coordinates x and y of end of element last
             CALL SAL141(IPAR(1,LAST),RPAR(:,LAST),X,Y,IEND)
             !              get next element in the perimeter in touch with
             !              last element at its end iend
             LGON=.TRUE.
             L=ID+1
             KEEP=0
             DO WHILE(LGON.AND.L<=I2)
                NEXT=PERIM(L)
                !                 check whether next is really 'next' to last
                CALL SAL142(X,Y,XNEW,YNEW,IPAR(1,NEXT),RPAR(:,NEXT),NEWEND,EPS2,DIST)
                IF(KEEP<MXKEEP)KEEP=KEEP+1
                KEEPD(KEEP)=REAL(DIST)
                LGON=NEWEND<0
                IF(LGON)L=L+1
             ENDDO
             LGOPEN=.NOT.LGON
             IF(LGOPEN)THEN
                !                 replace last element in order in the perimeter
                ID=ID+1
                IF(ID<L)THEN
                   PERIM(L)=PERIM(ID)
                   PERIM(ID)=NEXT
                ENDIF
             ELSE
                LGERR=.FALSE.
                !                 try to close with first element of loop:
                CALL SAL142(X,Y,XNEW,YNEW,IPAR(1,FIRST),RPAR(:,FIRST),NEWEND,EPS2,DIST)
                LGOPEN=NEWEND<0
                IF(LGOPEN)THEN
                   !                    fatal error: cannot close loop for this node
                   WRITE(FOUT,'(//,''==> cannot close node '',i5, &
                        &''   AT ELEMENT'',I5,/)')NODE,LAST
                   LGERR=.TRUE.
                ELSEIF(FIRST==LAST)THEN
                   WRITE(FOUT,'(//,''==> node '',i5,''    with '', &
                        &''isolated element'',i5,/)')NODE,LAST
                   LGERR=.TRUE.
                ENDIF
                IF(LGERR) CALL XABORT('SAL140: node not closed')
                NEXT=FIRST
             ENDIF
             !              define last = next element and proceed
             LAST=NEXT
             IEND=3-NEWEND
          ENDDO
          ID=ID+1
       ENDDO
    ENDDO
    !
  END SUBROUTINE SAL140
  !
  SUBROUTINE SAL142(X,Y,XNEW,YNEW,TYPE,RPAR,IEND,EPS2,DIST)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! checks whether element is very close to point (X,Y)
    !
    !Parameters: input
    ! X        abscissa coordinate
    ! Y        ordinate coordinate
    ! TYPE     type of element  1 (segment), 3 (arc of circle)
    ! EPS2     criterium for closeness
    !
    !Parameters: input/output
    ! XNEW     abscissa coordinate of end of element close to point
    ! YNEW     ordinate coordinate of end of element close to point
    ! RPAR     floating point geometry descriptors
    ! IEND     =  1 (beginning is close to point dist < EPS2)
    !             2 (end is close to point DIST < EPS2)
    !            -1 (beginning is close to point dist > EPS2)
    !            -2 (end is close to point DIST > EPS2)
    ! DIST     distance from end of element to point (X,Y)
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES, ONLY : G_ELE_TYPE
    !****
    IMPLICIT NONE
    INTEGER,   INTENT(IN)               :: TYPE
    INTEGER,   INTENT(OUT)              :: IEND
    REAL(PDB), INTENT(IN)               :: X,Y,EPS2
    REAL(PDB), INTENT(OUT)              :: XNEW,YNEW,DIST
    REAL(PDB), INTENT(IN), DIMENSION(:) :: RPAR
    !     DIMENSION        RPAR(NRPAR)
    !****
    REAL(PDB) :: CX,CY,THETA,R,DIST2
    !****
    !*    function giving distance between two points:
    REAL(PDB) :: FUNC,X1,Y1,X2,Y2
    FUNC(X1,Y1,X2,Y2)=(X1-X2)**2+(Y1-Y2)**2
    !****
    CX=RPAR(1)
    CY=RPAR(2)
    DIST2=0._PDB
    IF(TYPE==G_ELE_TYPE(1))THEN
       !*       segment
       XNEW=CX
       YNEW=CY
       DIST=FUNC(XNEW,YNEW,X,Y)
       IF(DIST<=EPS2)THEN
          IEND=1
          RETURN
       ELSE
          XNEW=CX+RPAR(3)
          YNEW=CY+RPAR(4)
          DIST2=FUNC(XNEW,YNEW,X,Y)
       ENDIF
    ELSEIF(TYPE<=G_ELE_TYPE(3))THEN
       !*       arc of circle
       R=RPAR(3)
       THETA=RPAR(4)
       XNEW=CX+R*COS(THETA)
       YNEW=CY+R*SIN(THETA)
       DIST=FUNC(XNEW,YNEW,X,Y)
       IF(DIST<=EPS2)THEN
          IEND=1
          RETURN
       ELSE
          THETA=RPAR(5)
          XNEW=CX+R*COS(THETA)
          YNEW=CY+R*SIN(THETA)
          DIST2=FUNC(XNEW,YNEW,X,Y)
       ENDIF
    ELSE
       CALL XABORT('SAL142: not implemented')
    ENDIF
    !
    IF(DIST2<=EPS2)THEN
       DIST=DIST2
       IEND=2
       RETURN
    ELSE
       IF(DIST<=DIST2)THEN
          IEND=-1
       ELSE
          DIST=DIST2
          IEND=-2
       ENDIF
    ENDIF
    !
  END SUBROUTINE SAL142
  !
  SUBROUTINE SAL160(GG)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! analyse domain definition: 2D volumes, surfaces
    !  - compute node volumes
    !  - compute areas of 2d surfaces
    !  - read medium data
    !
    !Parameters: input/output
    ! GG    geometry descriptor
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES, ONLY : G_BC_TYPE,NBMED
    !****
    IMPLICIT NONE
    TYPE(T_G_BASIC), INTENT(INOUT)  :: GG
    !****
    INTEGER :: INB,IMED
    INTEGER, DIMENSION(GG%NB_NODE) :: DATAIN
    !****
    !     SUBROUTINE SAL160_2(NB_ELEM,IPAR,RPAR,VOL2,SURF2,ISURF2_ELEM, NB_SURF2)
    CALL SAL160_2(GG%NB_ELEM,GG%IPAR,GG%RPAR,GG%VOL_NODE,GG%SURF2,GG%ISURF2_ELEM, &
         GG%NB_SURF2)
    !
    !*    read medium per region
    CALL SALGET(DATAIN,GG%NB_NODE,F_GEO,FOUT0,'media per node')
	
    ! number of media fixed to maximum of datain
    NBMED = MAXVAL(DATAIN(1:GG%NB_NODE))
    !
    !*    check and define med for code regions
    DO INB=1,GG%NB_NODE
       IMED=DATAIN(INB)
       IF(IMED>NBMED.OR.IMED<0)THEN
          WRITE(*,*) 'medium : ',IMED
          WRITE(*,*) 'inb, nbmed, nbnode : ',INB,NBMED,GG%NB_NODE
          CALL XABORT('SAL160: wrong medium in a region')
       ENDIF
       GG%MED(INB)=IMED
    ENDDO
    !
  END SUBROUTINE SAL160
  !
  SUBROUTINE SAL160_2(NB_ELEM,IPAR,RPAR,VOL2,SURF2,ISURF2_ELEM,NB_SURF2)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! compute 2D volumes and surfaces
    !  - compute node volumes
    !  - compute 2D surface areas
    !
    !Parameters: input
    ! NB_ELEM      number of elements
    ! NB_NODE      number of nodes
    ! IPAR         integer descriptors for elements
    ! RPAR         floating point descriptors for elements
    ! ISURF2_ELEM  2D surface nber per elem
    ! NB_SURF2     number of 2D surface
    !
    !Parameters: output
    ! VOL2         2D volumes of node 
    ! SURF2        2D areas of node
    !
    !---------------------------------------------------------------------
    !
    USE PRECISION_AND_KINDS, ONLY : PDB
    USE SAL_GEOMETRY_TYPES, ONLY : G_BC_TYPE
    !****
    IMPLICIT NONE
    INTEGER, INTENT(IN)                    :: NB_ELEM
    INTEGER, INTENT(IN),   DIMENSION(:,:)  :: IPAR
    REAL(PDB), INTENT(IN), DIMENSION(:,:)  :: RPAR
    REAL(PDB), INTENT(OUT),  DIMENSION(:)  :: VOL2,SURF2
    INTEGER, INTENT(IN),   DIMENSION(:)    :: ISURF2_ELEM
    INTEGER, INTENT(IN)                    :: NB_SURF2
    !****
    INTEGER :: ELEM,NODEBC,NODE,ISURF
    LOGICAL :: LGBC
    REAL(PDB) :: AUX
    !****
    !     initiation
    VOL2=0._PDB
    IF(NB_SURF2 > 0) SURF2=0._PDB
    DO ELEM=1,NB_ELEM
       !*       compute volume of node and add to volume of region (local)
       CALL SAL161(IPAR(1,ELEM),RPAR(:,ELEM),AUX)
       NODEBC=IPAR(2,ELEM)
       LGBC=NODEBC<=0
       IF(.NOT.LGBC) THEN
          VOL2(NODEBC)=VOL2(NODEBC)+AUX
       ENDIF
       NODE=IPAR(3,ELEM)
       IF(NODE>0)THEN
          VOL2(NODE)=VOL2(NODE)-AUX
       ELSE
          IF(LGBC) CALL XABORT('SAL160_2: isolated element')
          LGBC=.TRUE.
          NODEBC=NODE
       ENDIF
       IF(LGBC) THEN
          IF(NODEBC==G_BC_TYPE(-1).OR.NODEBC==G_BC_TYPE(0))THEN
             !  compute external surface (loaded in total order number)
             ISURF=ISURF2_ELEM(ELEM)
             IF(ISURF==0) CALL XABORT('SAL160_2: wrong relation of element - surf ')
             CALL SAL162(IPAR(1,ELEM),RPAR(:,ELEM),SURF2(ISURF))
          ENDIF
       ENDIF
    ENDDO
    !
  END SUBROUTINE SAL160_2
  !
  SUBROUTINE SAL161(TYPE,RPAR,VOL2)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! computes 'volume' between an element and the x axis
    !
    !Parameters: input
    ! TYPE         type of element
    ! RPAR         floating point descriptors for elements
    !
    !Parameters: output
    ! VOL2         '2D area' between the element and the horizontal axis
    !
    !---------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER,   INTENT(IN)                 :: TYPE
    REAL(PDB), INTENT(IN),  DIMENSION (:) :: RPAR
    REAL(PDB), INTENT(OUT)                :: VOL2
    !     DIMENSION        RPAR(*)
    !****
    REAL(PDB) :: YC,EX,EY,R,PHI1,PHI2,COS1,COS2
    !****
    !     volume is added to node- and substracted from node+
    YC=RPAR(2)
    IF(TYPE==1)THEN
       !        segment:
       EX=RPAR(3)
       EY=RPAR(4)
       VOL2=EX*(YC+EY/2.)
    ELSEIF(TYPE==2)THEN
       !        whole circle
       VOL2=PI*RPAR(3)*RPAR(3)
    ELSEIF(TYPE==3)THEN
       !        arc of circle:
       R=RPAR(3)
       PHI1=RPAR(4)
       PHI2=RPAR(5)
       COS1=COS(PHI1)
       COS2=COS(PHI2)
       VOL2=R*(YC*(COS1-COS2)+(R/2.0)*(PHI2-PHI1+ &
            & (SIN(PHI1)*COS1-SIN(PHI2)*COS2)))
    ELSE
       CALL XABORT('SAL161: not implemented')
    ENDIF
    !
  END SUBROUTINE SAL161
  !
  SUBROUTINE SAL162(TYPE,RPAR,SURF2)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! computes surface of a boundary element
    !
    !Parameters: input
    ! TYPE         type of element
    ! RPAR         floating point descriptors for elements
    !
    !Parameters: output
    ! SURF2        '2D length' of element
    !
    !---------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER,   INTENT(IN)                :: TYPE
    REAL(PDB), INTENT(IN), DIMENSION (:) :: RPAR
    REAL(PDB), INTENT(OUT)               :: SURF2
    !****
    IF(TYPE==1)THEN
       !        segment
       SURF2=RPAR(5)
    ELSEIF(TYPE<=3)THEN
       !        arc of circle
       SURF2=(RPAR(5)-RPAR(4))*RPAR(3)
    ELSE
       CALL XABORT('SAL162: not implemented')
    ENDIF
    !
  END SUBROUTINE SAL162
  !
  SUBROUTINE SAL170(GG)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! prints out volume, surface and medium info for 2D macros
    !
    !Parameters: input/output
    ! GG    geometry descriptor
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES,      ONLY : G_BC_TYPE
    USE SAL_TRACKING_TYPES,      ONLY : PRTIND
    !****
    IMPLICIT NONE
    TYPE(T_G_BASIC), INTENT(INOUT)  :: GG
    INTEGER, PARAMETER :: FOUT =6

    !****
     IF(PRTIND == 1) THEN
         WRITE(FOUT,'(//,10x,''2D geometry'',/,10X,11(''=''),//, &
         &I8,''  regions'',/, &
         &I8,''  external surfaces'',//)') GG%NB_NODE,GG%NB_SURF2
     ENDIF
    !
  END SUBROUTINE SAL170
  !
  SUBROUTINE SALSYM(X1,Y1,X2,Y2,X0,Y0,X4,Y4)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! set the symmetric point (x4,y4) of (x0,y0) relative to symmetry axis
    ! (x1,y1)->(x2,y2)
    !
    !Parameters: input
    ! X1   abscissa coordinate of a point on the symmetry axis
    ! Y1   ordinate coordinate of a point on the symmetry axis
    ! X2   abscissa coordinate of a point on the symmetry axis
    ! Y2   ordinate coordinate of a point on the symmetry axis
    ! X0   abscissa coordinate of the point to mirror
    ! Y0   ordinate coordinate of the point to mirror
    !
    !Parameters: output
    ! X4   abscissa coordinate of the symmetric point
    ! Y4   ordinate coordinate of the symmetric point
    !
    !---------------------------------------------------------------------
    !
    IMPLICIT NONE
    REAL(PDB),INTENT(IN)  :: X1,Y1,X2,Y2,X0,Y0
    REAL(PDB),INTENT(OUT) :: X4,Y4

    REAL(PDB) :: A,B,C,DEN,X3,Y3
    A=Y2-Y1; B=X1-X2; C=X2*Y1-X1*Y2;
    DEN=A*A+B*B;
    IF(DEN==0._PDB) CALL XABORT('SALSYM: division by zero')
    X3=(B*(B*X0-A*Y0)-A*C)/DEN;
    Y3=(A*(-B*X0+A*Y0)-B*C)/DEN;
    X4=X0+2._PDB*(X3-X0); Y4=Y0+2._PDB*(Y3-Y0);
  END SUBROUTINE SALSYM
  !
  SUBROUTINE SALFOLD(HSYM,GG)
    ! unfold the domain
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! unfold the domain
    !
    !Parameters: input
    ! HSYM: type of symmetry: 'DIAG': diagonal symmetry; 'SYMX': symmetry
    !       relative to X axis; 'SYMY': symmetry relative to Y axis
    !
    !Parameters: input/output
    ! GG    geometry descriptor
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES, ONLY : NIPAR,NRPAR,ANGGEO,LENGTHX,LENGTHY
    IMPLICIT NONE
    CHARACTER(LEN=4),INTENT(IN)  :: HSYM
    TYPE(T_G_BASIC), INTENT(INOUT) :: GG
    !
    INTEGER, POINTER, DIMENSION(:,:) :: TMP_IPAR
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: I2
    INTEGER, ALLOCATABLE, DIMENSION(:) :: I1
    REAL(PDB), POINTER, DIMENSION(:,:) :: TMP_RPAR
    TYPE(T_SALBCDATA), POINTER, DIMENSION(:) :: TMP_BCDATAREAD
    INTEGER :: ELEM,TYPE,OK,TMP_NB_ELEM,TMP_NBBCDA,I,IBC,INDBC,TMP_IBC,IAUX
    REAL(PDB),PARAMETER :: EPS=1.0E-5_PDB
    REAL(PDB) :: X1,X2,X4,Y1,Y2,Y4,AXIS_X1,AXIS_X2,AXIS_Y1,AXIS_Y2, &
                 DX4,DY4,RAD,THETA1,THETA2,X1B,Y1B,X4B,Y4B,XMIN,YMIN, &
                 XMAX,YMAX,PHI1,PHI2,DELPHI
    !
    ! define symmetry axis
    IF(HSYM=='SYMX') THEN
      AXIS_X1=0._PDB; AXIS_X2=100._PDB; AXIS_Y1=0._PDB; AXIS_Y2=0._PDB;
    ELSE IF(HSYM=='SYMY') THEN
      AXIS_X1=0._PDB; AXIS_X2=0._PDB; AXIS_Y1=0._PDB; AXIS_Y2=100._PDB;
    ELSE IF(HSYM=='DIAG') THEN
      AXIS_X1=0._PDB; AXIS_X2=100._PDB; AXIS_Y1=0._PDB; AXIS_Y2=100._PDB;
      ANGGEO=2.0*ANGGEO
    ELSE
      CALL XABORT('SALFOLD: invalid type of symmetry axis')
    ENDIF
    XMIN=1.E10_PDB; YMIN=1.E10_PDB; XMAX=-1.E10_PDB; YMAX=-1.E10_PDB;
    !
    ! allocate new surfacic element containers
    ALLOCATE(TMP_IPAR(NIPAR,2*GG%NB_ELEM), TMP_RPAR(NRPAR,2*GG%NB_ELEM), &
             I2(2,GG%NB_ELEM), STAT=OK)
    IF(OK/=0) CALL XABORT('SALFOLD: not enough memory')
    TMP_IPAR(:,:)=0; TMP_RPAR(:,:)=0._PDB;
    !
    ! loop over old elements
    TMP_NB_ELEM=0
    THETA1=0._PDB; THETA2=0._PDB;
    DO ELEM=1,GG%NB_ELEM
      I2(:,ELEM)=0
      TYPE=GG%IPAR(1,ELEM)
      IF(TYPE==1) THEN
        X1=GG%RPAR(1,ELEM); Y1=GG%RPAR(2,ELEM);
        XMIN=MIN(XMIN,X1); YMIN=MIN(YMIN,Y1); XMAX=MAX(XMAX,X1); YMAX=MAX(YMAX,Y1);
        X2=X1+GG%RPAR(3,ELEM); Y2=Y1+GG%RPAR(4,ELEM);
        CALL SALSYM(AXIS_X1,AXIS_Y1,AXIS_X2,AXIS_Y2,X1,Y1,X4,Y4)
        CALL SALSYM(AXIS_X1,AXIS_Y1,AXIS_X2,AXIS_Y2,X2,Y2,DX4,DY4)
        IF(ABS(X1-X4)<EPS .AND. ABS(Y1-Y4)<EPS .AND. &
           ABS(X2-DX4)<EPS .AND. ABS(Y2-DY4)<EPS) CYCLE
      ELSE IF(TYPE==2) THEN
        X1=GG%RPAR(1,ELEM); Y1=GG%RPAR(2,ELEM); RAD=GG%RPAR(3,ELEM)
        CALL SALSYM(AXIS_X1,AXIS_Y1,AXIS_X2,AXIS_Y2,X1,Y1,X4,Y4)
        THETA1=0._PDB; THETA2=0._PDB;
      ELSE IF(TYPE==3) THEN
        X1=GG%RPAR(1,ELEM); Y1=GG%RPAR(2,ELEM); RAD=GG%RPAR(3,ELEM)
        CALL SALSYM(AXIS_X1,AXIS_Y1,AXIS_X2,AXIS_Y2,X1,Y1,X4,Y4)
        X1B=X1+RAD*COS(GG%RPAR(4,ELEM)); Y1B=Y1+RAD*SIN(GG%RPAR(4,ELEM));
        CALL SALSYM(AXIS_X1,AXIS_Y1,AXIS_X2,AXIS_Y2,X1B,Y1B,X4B,Y4B)
        IF((ABS(X4B-X4)<EPS*ABS(RAD)).AND.(Y4B-Y4 > 0._PDB)) THEN
          THETA1=PI/2._PDB
        ELSE IF((ABS(X4B-X4)<EPS*ABS(RAD)).AND.(Y4B-Y4 < 0._PDB)) THEN
          THETA1=3._PDB*PI/2._PDB
        ELSE IF(X4B-X4 > 0._PDB) THEN
          THETA1=ATAN((Y4B-Y4)/(X4B-X4))
        ELSE
          THETA1=ATAN((Y4B-Y4)/(X4B-X4))+PI
        ENDIF
        X1B=X1+RAD*COS(GG%RPAR(5,ELEM)); Y1B=Y1+RAD*SIN(GG%RPAR(5,ELEM));
        CALL SALSYM(AXIS_X1,AXIS_Y1,AXIS_X2,AXIS_Y2,X1B,Y1B,X4B,Y4B)
        IF((ABS(X4B-X4)<EPS*ABS(RAD)).AND.(Y4B-Y4 > 0._PDB)) THEN
          THETA2=PI/2._PDB
        ELSE IF((ABS(X4B-X4)<EPS*ABS(RAD)).AND.(Y4B-Y4 < 0._PDB)) THEN
          THETA2=3._PDB*PI/2._PDB
        ELSE IF(X4B-X4 > 0._PDB) THEN
          THETA2=ATAN((Y4B-Y4)/(X4B-X4))
        ELSE
          THETA2=ATAN((Y4B-Y4)/(X4B-X4))+PI
        ENDIF
      ELSE
        WRITE(*,*) " elem=",ELEM," type=",TYPE
        CALL XABORT('SALFOLD: invalid type of surfacic element')
      ENDIF
      TMP_NB_ELEM=TMP_NB_ELEM+1
      IF(TMP_NB_ELEM>2*GG%NB_ELEM) CALL XABORT('SALFOLD: tmp_nb_elem overflow(1)')
      TMP_IPAR(:,TMP_NB_ELEM)=GG%IPAR(:,ELEM)
      TMP_RPAR(:,TMP_NB_ELEM)=GG%RPAR(:,ELEM)
      I2(1,ELEM)=TMP_NB_ELEM
      TMP_NB_ELEM=TMP_NB_ELEM+1
      IF(TMP_NB_ELEM>2*GG%NB_ELEM) CALL XABORT('SALFOLD: tmp_nb_elem overflow(2)')
      I2(2,ELEM)=TMP_NB_ELEM
      TMP_RPAR(1,TMP_NB_ELEM)=X4; TMP_RPAR(2,TMP_NB_ELEM)=Y4;
      TMP_IPAR(1,TMP_NB_ELEM)=TYPE;
      IF(TYPE==1) THEN
        TMP_RPAR(3,TMP_NB_ELEM)=DX4-X4; TMP_RPAR(4,TMP_NB_ELEM)=DY4-Y4;
        TMP_RPAR(5,TMP_NB_ELEM)=SQRT(TMP_RPAR(3,TMP_NB_ELEM)**2+TMP_RPAR(4,TMP_NB_ELEM)**2)
        XMIN=MIN(XMIN,X4); YMIN=MIN(YMIN,Y4); XMAX=MAX(XMAX,X4); YMAX=MAX(YMAX,Y4);
        TMP_IPAR(2,TMP_NB_ELEM)=GG%IPAR(3,ELEM); TMP_IPAR(3,TMP_NB_ELEM)=GG%IPAR(2,ELEM);
      ELSE IF((TYPE==2).OR.(TYPE==3)) THEN
        TMP_RPAR(3,TMP_NB_ELEM)=GG%RPAR(3,ELEM) ! RADIUS
        IF(THETA2>THETA1) THETA1=THETA1+2._PDB*PI
        PHI1=THETA2; DELPHI=THETA1-THETA2;
        IF(DELPHI>0._PDB)THEN
          PHI2=PHI1+DELPHI
        ELSE
          PHI2=PHI1
          PHI1=PHI1+DELPHI
        ENDIF
        IF(TYPE==3)THEN
          ! arc of circle: put phi1 within 0 and 2*pi
          IF(PHI1>2._PDB*PI)THEN
             IAUX=INT(PHI1/(2._PDB*PI))
             DELPHI=(2._PDB*PI)*IAUX
             PHI2=PHI2-DELPHI
             PHI1=PHI1-DELPHI
          ELSEIF(PHI1<0._PDB)THEN
             IAUX=INT((-PHI1+1.D-7)/(2._PDB*PI))+1
             DELPHI=(2._PDB*PI)*IAUX
             PHI2=PHI2+DELPHI
             PHI1=PHI1+DELPHI
          ENDIF
        ENDIF
        TMP_RPAR(4,TMP_NB_ELEM)=PHI1; TMP_RPAR(5,TMP_NB_ELEM)=PHI2; ! ANGLES
        TMP_IPAR(2,TMP_NB_ELEM)=GG%IPAR(2,ELEM); TMP_IPAR(3,TMP_NB_ELEM)=GG%IPAR(3,ELEM);
      ENDIF
      TMP_RPAR(6,TMP_NB_ELEM)=0._PDB
    ENDDO
    DEALLOCATE(GG%IPAR,GG%RPAR)
    GG%IPAR=>TMP_IPAR; GG%RPAR=>TMP_RPAR;
    GG%NB_ELEM=TMP_NB_ELEM
    !
    ! translate the domain
    DO ELEM=1,GG%NB_ELEM
      GG%RPAR(1,ELEM)=GG%RPAR(1,ELEM)-XMIN
      GG%RPAR(2,ELEM)=GG%RPAR(2,ELEM)-YMIN
    ENDDO
    LENGTHX=XMAX-XMIN ; LENGTHY=YMAX-YMIN ;
    !
    ! allocate new boundary condition containers
    ALLOCATE(I1(GG%NBBCDA))
    I1(:)=0
    !
    ! loop over boundary conditions
    TMP_NBBCDA=0
    DO IBC=1,GG%NBBCDA
      DO I=1,GG%BCDATAREAD(IBC)%NBER
        INDBC=GG%BCDATAREAD(IBC)%ELEMNB(I)
        IF(INDBC==0) CYCLE
        IF(I2(1,INDBC)/=0) GO TO 10
      ENDDO
      CYCLE
      10 TMP_NBBCDA=TMP_NBBCDA+1
      I1(IBC)=TMP_NBBCDA
    ENDDO
    ALLOCATE(TMP_BCDATAREAD(2*TMP_NBBCDA))
    DO IBC=1,GG%NBBCDA
      TMP_IBC=I1(IBC)
      IF(TMP_IBC==0) CYCLE
      IF(GG%BCDATAREAD(IBC)%SALTYPE==0) THEN
        ALLOCATE(TMP_BCDATAREAD(TMP_IBC)%ELEMNB(2*GG%BCDATAREAD(IBC)%NBER))
        DO I=1,GG%BCDATAREAD(IBC)%NBER
          INDBC=GG%BCDATAREAD(IBC)%ELEMNB(I)
          TMP_BCDATAREAD(TMP_IBC)%ELEMNB(2*I-1)=I2(1,INDBC)
          TMP_BCDATAREAD(TMP_IBC)%ELEMNB(2*I)=I2(2,INDBC)
        ENDDO
        TMP_BCDATAREAD(TMP_IBC)%SALTYPE=GG%BCDATAREAD(IBC)%SALTYPE
        TMP_BCDATAREAD(TMP_IBC)%NBER=2*GG%BCDATAREAD(IBC)%NBER
        TMP_BCDATAREAD(TMP_IBC)%BCDATA=GG%BCDATAREAD(IBC)%BCDATA
      ELSE
        ALLOCATE(TMP_BCDATAREAD(TMP_IBC)%ELEMNB(GG%BCDATAREAD(IBC)%NBER))
        DO I=1,GG%BCDATAREAD(IBC)%NBER
          INDBC=GG%BCDATAREAD(IBC)%ELEMNB(I)
          TMP_BCDATAREAD(TMP_IBC)%ELEMNB(I)=I2(1,INDBC)
        ENDDO
        TMP_BCDATAREAD(TMP_IBC)%SALTYPE=GG%BCDATAREAD(IBC)%SALTYPE
        TMP_BCDATAREAD(TMP_IBC)%NBER=GG%BCDATAREAD(IBC)%NBER
        TMP_BCDATAREAD(TMP_IBC)%BCDATA=GG%BCDATAREAD(IBC)%BCDATA
        !
        TMP_NBBCDA=TMP_NBBCDA+1
        ALLOCATE(TMP_BCDATAREAD(TMP_NBBCDA)%ELEMNB(GG%BCDATAREAD(IBC)%NBER))
        DO I=1,GG%BCDATAREAD(IBC)%NBER
          INDBC=GG%BCDATAREAD(IBC)%ELEMNB(I)
          TMP_BCDATAREAD(TMP_NBBCDA)%ELEMNB(I)=I2(2,INDBC)
        ENDDO
        TMP_BCDATAREAD(TMP_NBBCDA)%SALTYPE=GG%BCDATAREAD(IBC)%SALTYPE
        TMP_BCDATAREAD(TMP_NBBCDA)%NBER=GG%BCDATAREAD(IBC)%NBER
        X1=GG%BCDATAREAD(IBC)%BCDATA(1)
        Y1=GG%BCDATAREAD(IBC)%BCDATA(2)
        IF(GG%BCDATAREAD(IBC)%SALTYPE >= 2) THEN
          CALL SALSYM(AXIS_X1,AXIS_Y1,AXIS_X2,AXIS_Y2,X1,Y1,X4,Y4)
          TMP_BCDATAREAD(TMP_NBBCDA)%BCDATA(1)=X4
          TMP_BCDATAREAD(TMP_NBBCDA)%BCDATA(2)=Y4
        ELSE
          ! X1 is the albedo
          TMP_BCDATAREAD(TMP_NBBCDA)%BCDATA(1)=X1
          TMP_BCDATAREAD(TMP_NBBCDA)%BCDATA(2)=Y1
        ENDIF
        IF((HSYM=='DIAG').AND.GG%BCDATAREAD(IBC)%BCDATA(3).EQ.0.0) THEN
          TMP_BCDATAREAD(TMP_NBBCDA)%BCDATA(3)=90.0
        ELSE IF((HSYM=='DIAG').AND.GG%BCDATAREAD(IBC)%BCDATA(3).EQ.90.0) THEN
          TMP_BCDATAREAD(TMP_NBBCDA)%BCDATA(3)=0.0
        ELSE
          TMP_BCDATAREAD(TMP_NBBCDA)%BCDATA(3)=GG%BCDATAREAD(IBC)%BCDATA(3)
        ENDIF
      ENDIF
      DEALLOCATE(GG%BCDATAREAD(IBC)%ELEMNB)
    ENDDO
    DEALLOCATE(GG%BCDATAREAD)
    GG%BCDATAREAD=>TMP_BCDATAREAD
    GG%NBBCDA=TMP_NBBCDA
    !
    DEALLOCATE(I1,I2)
  END SUBROUTINE SALFOLD
  !
END MODULE SAL_GEOMETRY_MOD
