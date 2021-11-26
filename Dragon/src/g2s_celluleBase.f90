!
!-----------------------------------------------------------------------
!
!Purpose:
! Creation of an array of type(t_celluleBase) structures.
!
!Copyright:
! Copyright (C) 2001 Ecole Polytechnique de Montreal.
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s):
! G. Civario (CS-SI)
!
!Comments:
! Il s'agit de toutes les cellules terminales de la geometrie (=celles qui
! ne comportent pas de sous-cellules).
! La creation se fait a l'aide de la routine recursive buildCellsBase.
! Elle explore toute l'arboressence de l'objet PyLCM donne en entree, et appelle
!  a chaque niveau la routine writeCellBase. 
! Celle-ci teste si la cellule est terminale, et cree une entree dans le tableau
! des cellules de base si c'est le cas.
! La cellule creee est ensuite completee si besoin, et sa coherence est
! verifiee
! \\\\
! variable globale
!  - tabCelluleBase : tableau des cellules de base
! \\\\
! fonctions du module
!  - initializeTabCelluleBase : mise a zero du tableau
!  - destoyTabCelluleBase : liberation de la memoire
!  - decale : decale les valeurs d'un tableau pour les faire demarer a 0
!  - createCB : constructeur d'une cellule de base
!  - destroyCB : destructeur d'une cellule de base
!  - createCluster : constructeur d'un cluster
!  - sortClusterTab : trie un tableau de clusters
!  - verrifieCB : verification de la coherence d'une cellule de base
!  - exploiteSplit : exploitation de la donnee split dans la cellule
!  - buildCellsBase : fonction recursive de creation du tableau de cellules
!  - writeCellBase : creation d'une cellule
!
!-----------------------------------------------------------------------
!
module celluleBase

  use constType
  use GANLIB
  use segArc, only : alloc_ok

  implicit none

  !cluster
  type t_cluster
     character*12                          :: name        !nom
     integer                               :: nbrPin      !nombre de crayons
     double precision                      :: radiusOfPin !rayon de la couronne
     double precision                      :: angleOfPin  !angle du 1er crayon
     double precision,dimension(:),pointer :: radius      !rayons des anneaux
     integer,dimension(:),pointer          :: mix         !milieux des crayons
  end type t_cluster

  !cellule generique de plus bas niveau hierarchique dans les gigognes
  type t_celluleBase
     character*12                          :: name      !nom
     integer,dimension(40)                 :: sv        !state vector
     integer,dimension(:),allocatable          :: mix       !milieux
     integer,dimension(:),allocatable          :: merge     !regroupements
     double precision,dimension(:),allocatable :: radius    !rayons
     double precision,dimension(3)         :: offcenter !x , y et z
     double precision,dimension(:),allocatable :: meshx     !en commancant a 0.0
     double precision,dimension(:),allocatable :: meshy     !en commancant a 0.0
     double precision                      :: side      !pour tri et hex
     integer,dimension(:),allocatable          :: splitr    !>0  rayon; <0  surface
     integer,dimension(:),allocatable          :: splitx    !>0
     integer,dimension(:),allocatable          :: splity    !>0
     type(t_cluster),dimension(:),pointer  :: cluster   !les clusters
     !tableau de bool qui donne la presence ou non de chacun des 10 champs
     !avec en plus l'indice 0 disant si la cellule est viable (complete)
     logical,dimension(0:12)               :: ok
  end type t_celluleBase

  !parametres pour le vecteur ok
  integer,parameter :: n_sv=1 , n_mix=2 , n_radius=3 , n_offcenter=4 , &
                     & n_meshx=5 , n_meshy=6 , n_side=7 , n_splitr=8 , &
                     & n_splitx=9 , n_splity=10 , n_cluster=11 , n_merge=12 , &
                     & n_tot=0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  declaration d'une variable tableau globale
  !!
  type(t_celluleBase),dimension(:),allocatable :: tabCelluleBase
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  subroutine initializeTabCelluleBase(dimTabCelluleBase)
    integer,intent(in) :: dimTabCelluleBase

    allocate(tabCelluleBase(dimTabCelluleBase),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: initializeTabCelluleBase(1) => allocation pb")
  end subroutine initializeTabCelluleBase
  
  subroutine destroyTabCelluleBase(szB)
    integer,intent(in) :: szB
    integer :: i

    do i = 1,szB
       call destroyCB(tabCelluleBase(i))
    end do
    deallocate(tabCelluleBase)
   end subroutine destroyTabCelluleBase

  !fait demarer un tableau de reels a 0.0
  subroutine decale(tab)
!    double precision,dimension(:),pointer :: tab
    double precision,dimension(:) :: tab
    integer :: i, lg
    
    lg = size(tab)
    do i = lg,1,-1
       tab(i) = tab(i) - tab(1)
    end do
  end subroutine decale
  
  !remplisage d'une cellule de base
  subroutine createCB(cell,name,ip)
    type(t_celluleBase),intent(out)  :: cell
    character*12,intent(in)          :: name
    type(c_ptr),intent(in)           :: ip

    integer :: lg,typ,i,lgm
    real,dimension(:),allocatable         :: tmpTabReal
    character*12,dimension(:),allocatable :: clusterName

    cell%name = name
    
    call LCMGET(ip,'STATE-VECTOR',cell%sv)
    cell%ok(n_sv) = .true.

    call LCMLEN(ip,'MIX         ',lgm,typ)
    cell%ok(n_mix) = (lgm/=0)
    if(cell%ok(n_mix)) then
       allocate(cell%mix(lgm),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: createCB(1) => allocation pb")
       call LCMGET(ip,'MIX         ',cell%mix)
    else
       call XABORT("G2S: no mix in the cellule " // name)
    end if

    call LCMLEN(ip,'MERGE       ',lg,typ)
    cell%ok(n_merge) = (lg/=0)
    if(cell%ok(n_merge)) then
       if(lg/=lgm) call XABORT("G2S: bad dimension for merge in the &
            &cellule " // name)
       allocate(cell%merge(lg),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: createCB(2) => allocation pb")
       call LCMGET(ip,'MERGE       ',cell%merge)
    else
       allocate(cell%merge(lgm),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: createCB(3) => allocation pb")
       cell%merge(:lgm) = (/(i,i=1,lgm)/)
    end if

   call LCMLEN(ip,'RADIUS      ',lg,typ)
    cell%ok(n_radius) = (lg/=0)
    if(cell%ok(n_radius)) then
       allocate(cell%radius(lg),tmpTabReal(lg),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: createCB(4) => allocation pb")
       call LCMGET(ip,'RADIUS      ',tmpTabReal)
       cell%radius(:lg)=tmpTabReal(:lg)
       deallocate(tmpTabReal)
    else if(cell%sv(1)==G_Carcel.or.cell%sv(1)==G_Hexcel) then
       !cas carcel 0 ou hexcel 0
       allocate(cell%radius(1))
       cell%radius(1)=0.d0
    else 
       !nullify(cell%radius)
    end if
    
    call LCMLEN(ip,'OFFCENTER   ',lg,typ)
    cell%ok(n_offcenter) = (lg/=0)
    if(cell%ok(n_offcenter)) then
       allocate(tmpTabReal(lg),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: createCB(5) => allocation pb")
       call LCMGET(ip,'OFFCENTER   ',tmpTabReal)
       cell%offcenter=tmpTabReal
       deallocate(tmpTabReal)
    else
       cell%offcenter=0.d0
    end if
    
    call LCMLEN(ip,'MESHX       ',lg,typ)
    cell%ok(n_meshx) = (lg/=0)
    if(cell%ok(n_meshx)) then
       allocate(cell%meshx(lg),tmpTabReal(lg),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: createCB(6) => allocation pb")
       call LCMGET(ip,'MESHX       ',tmpTabReal)
       cell%meshx(:lg)=tmpTabReal(:lg)
       deallocate(tmpTabReal)
       call decale(cell%meshx) !decalage pour demarer a 0.0
    else
       !nullify(cell%meshx)
    end if
    
    call LCMLEN(ip,'MESHY       ',lg,typ)
    cell%ok(n_meshy) = (lg/=0)
    if(cell%ok(n_meshy)) then
       allocate(cell%meshy(lg),tmpTabReal(lg),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: createCB(7) => allocation pb")
       call LCMGET(ip,'MESHY       ',tmpTabReal)
       cell%meshy(:lg)=tmpTabReal(:lg)
       deallocate(tmpTabReal)
       call decale(cell%meshy) !decalage pour demarer a 0.0
    else
       !nullify(cell%meshy)
    end if

    call LCMLEN(ip,'SIDE        ',lg,typ)
    cell%ok(n_side) = (lg/=0)
    if(cell%ok(n_side)) then
       allocate(tmpTabReal(lg),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: createCB(8) => allocation pb")
       call LCMGET(ip,'SIDE        ',tmpTabReal)
       cell%side=tmpTabReal(1)
       deallocate(tmpTabReal)
    else
       cell%side=0.d0
    end if

    call LCMLEN(ip,'SPLITR      ',lg,typ)
    cell%ok(n_splitr) = (lg/=0)
    if(cell%ok(n_splitr)) then
       allocate(cell%splitr(lg),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: createCB(9) => allocation pb")
       call LCMGET(ip,'SPLITR      ',cell%splitr)
    else
       !nullify(cell%splitr)
    end if

    call LCMLEN(ip,'SPLITX      ',lg,typ)
    cell%ok(n_splitx) = (lg/=0)
    if(cell%ok(n_splitx)) then
       allocate(cell%splitx(lg),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: createCB(10) => allocation pb")
       call LCMGET(ip,'SPLITX      ',cell%splitx)
    else if(cell%sv(1)==G_Car2d.or.cell%sv(1)==G_Carcel) then
       allocate(cell%splitx(cell%sv(3)),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: createCB(11) => allocation pb")
       cell%splitx = 1
    else
       allocate(cell%splitx(1))
       cell%splitx(1) = 1
    end if

    call LCMLEN(ip,'SPLITY      ',lg,typ)
    cell%ok(n_splity) = (lg/=0)
    if(cell%ok(n_splity)) then
       allocate(cell%splity(lg),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: createCB(12) => allocation pb")
       call LCMGET(ip,'SPLITY      ',cell%splity)
    else if(cell%sv(1)==G_Car2d.or.cell%sv(1)==G_Carcel) then
       allocate(cell%splity(cell%sv(4)),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: createCB(13) => allocation pb")
       cell%splity = 1
    else
       allocate(cell%splity(1))
       cell%splity(1) = 1
    end if

    call LCMLEN(ip,'CLUSTER     ',lg,typ)
    cell%ok(n_cluster) = (lg/=0)
    if(cell%ok(n_cluster)) then
       lg = lg/3
       allocate(cell%cluster(lg),clusterName(lg),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: createCB(14) => allocation pb")
       call LCMGTC(ip,'CLUSTER     ',12,lg,clusterName)
       do i = 1,lg
          cell%cluster(i) = createCluster(ip,clusterName(i))
       end do
       !on trie le tableau des cluster par rayon croissant
       call sortClusterTab(cell%cluster)
       deallocate(clusterName)
    else
       nullify(cell%cluster)
    end if

    call verrifieCB(cell)
    call exploiteSplit(cell)
  end subroutine createCB

  subroutine destroyCB(cell)
    type(t_celluleBase),intent(inout) :: cell
    integer :: i

    if(allocated(cell%mix)) deallocate(cell%mix)
    if(allocated(cell%merge)) deallocate(cell%merge)
    if(allocated(cell%radius)) deallocate(cell%radius)
    if(allocated(cell%meshx)) deallocate(cell%meshx)
    if(allocated(cell%meshy)) deallocate(cell%meshy)
    if(allocated(cell%splitr)) deallocate(cell%splitr)
    if(allocated(cell%splitx)) deallocate(cell%splitx)
    if(allocated(cell%splity)) deallocate(cell%splity)
    if(associated(cell%cluster)) then
       do i = 1,size(cell%cluster)
          deallocate(cell%cluster(i)%mix,cell%cluster(i)%radius)
          !nullify(cluster%mix,cell%cluster(i)%radius)
       end do
       deallocate(cell%cluster)
    end if
  end subroutine destroyCB

  function createCluster(cellBIp,clusterName)
    type(c_ptr),intent(in)      :: cellBIp
    character*12,intent(in) :: clusterName
    type(t_cluster)         :: createCluster

    type(c_ptr) :: ip
    integer :: lg,ty
    real    :: rpin,apin
    real,dimension(:),allocatable :: radius

    ip = cellBIp
    call LCMSIX(ip,clusterName,1)

    createCluster%name = clusterName

    call LCMGET(ip,'NPIN        ',createCluster%nbrPin)
    
    call LCMGET(ip,'RPIN        ',rpin)
    createCluster%radiusOfPin = rpin

    call LCMGET(ip,'APIN        ',apin)
    createCluster%angleOfPin = apin

    call LCMLEN(ip,'MIX         ',lg,ty)
    allocate(radius(lg+1),createCluster%radius(lg+1),createCluster%mix(lg),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: createCluster => allocation pb")
    call LCMGET(ip,'RADIUS      ',radius)
    createCluster%radius(:lg+1) = radius(:lg+1)
    deallocate(radius)

    call LCMGET(ip,'MIX         ',createCluster%mix)
  end function createCluster

  subroutine sortClusterTab(clusterPrt)
    type(t_cluster),dimension(:),pointer :: clusterPrt
    !trie le tableau des clusters par rayon croissant des couronnes
    type(t_cluster) :: tmpCluster
    integer         :: sz,indMax,i
    logical         :: trie
 
    if(.not. associated(clusterPrt)) return !pas de cluster
    sz = size(clusterPrt)
    do indMax = sz,2,-1
       trie = .true.
       do i = 1,indMax-1
          if(clusterPrt(i+1)%radiusOfPin<clusterPrt(i)%radiusOfPin) then
             tmpCluster=clusterPrt(i+1)
             clusterPrt(i+1)=clusterPrt(i)
             clusterPrt(i)=tmpCluster
             trie = .false.
          end if
       end do
       if(trie) exit
    end do
  end subroutine sortClusterTab

  !verification de la validite d'une cellule en fonction de son vecteur ok
  !et mise a jour du ok(0)
  subroutine verrifieCB(cell)
    type(t_celluleBase),intent(inout) :: cell

    integer :: typCell,sectori,sectorj,nx,ny

    typCell = cell%sv(1)
    sectori = cell%sv(14)
    sectorj = cell%sv(15)
    select case(typCell)
    case(G_Virtual)
       call XABORT("G2S: Type of geometry VIRTUAL not supported")
    case(G_Car1d)
       cell%ok(n_tot) = cell%ok(n_sv) .and. cell%ok(n_mix) &
            & .and. cell%ok(n_meshx) &
            & .and. .not. &
            & ( cell%ok(n_radius) .or. cell%ok(n_offcenter) &
            & .or. cell%ok(n_meshy) .or. cell%ok(n_side) &
            & .or. cell%ok(n_splitr) .or. cell%ok(n_splity) )
    case(G_Tube)
       cell%ok(n_tot) = cell%ok(n_sv) .and. cell%ok(n_mix) &
            & .and. cell%ok(n_radius) &
            & .and. .not. &
            & ( cell%ok(n_offcenter) .or. cell%ok(n_meshx) &
            & .or. cell%ok(n_meshy) .or. cell%ok(n_side) &
            & .or. cell%ok(n_splitx) .or. cell%ok(n_splity) )
       if(sectori/=S_not) call XABORT("G2S: SECT not allowed for TUBE")
    case(G_Car2d)
       cell%ok(n_tot) = cell%ok(n_sv) .and. cell%ok(n_mix) &
            & .and. cell%ok(n_meshx) .and. cell%ok(n_meshy) &
            & .and. .not. &
            & ( cell%ok(n_radius) .or. cell%ok(n_offcenter) &
            & .or. cell%ok(n_side) .or. cell%ok(n_splitr) )
       if(sectori/=S_not) call XABORT("G2S: SECT not allowed for CAR2D")
    case(G_Hex)
       cell%ok(n_tot) = cell%ok(n_sv) .and. cell%ok(n_mix) &
            & .and. cell%ok(n_side) &
            & .and. .not. &
            & ( cell%ok(n_meshx) .or. cell%ok(n_meshy) &
            & .or. cell%ok(n_splitr) .or. cell%ok(n_splitx) &
            & .or. cell%ok(n_splity) .or. cell%ok(n_radius) )
       if((sectori/=S_not).and.(size(cell%mix)/=6)) &
            call XABORT("G2S: wrong len for MIX with SECT")
    case(G_Tri)
       cell%ok(n_tot) = cell%ok(n_sv) .and. cell%ok(n_mix) &
            & .and. cell%ok(n_side) &
            & .and. .not. cell%ok(n_radius)
       if(sectori/=S_not) call XABORT("G2S: SECT not allowed for TRI2D")
    case(G_Carcel)
       cell%ok(n_tot) = cell%ok(n_sv) .and. cell%ok(n_mix) &
            & .and. cell%ok(n_radius) .and. cell%ok(n_meshx) &
            & .and. cell%ok(n_meshy) &
            & .and. .not. &
            & ( cell%ok(n_side) )
       select case(sectori)
       case(S_not)
          nx = size(cell%meshx)-1
          ny = size(cell%meshy)-1
          if(size(cell%mix)/=nx*ny*size(cell%radius)) &
               call XABORT("G2S: wrong len for MIX")
       case(S_X_tot,S_T_tot)
          if(size(cell%mix)/=4*size(cell%radius)-3*sectorj) &
               call XABORT("G2S: wrong len for MIX with SECT")
       case(S_TX_tot,S_TXS_tot)
          if(size(cell%mix)/=8*size(cell%radius)-7*sectorj) &
               call XABORT("G2S: wrong len for MIX with SECT")
       case(S_WM_tot)
          if(size(cell%mix)/=8*size(cell%radius)-7*sectorj+4) &
               call XABORT("G2S: wrong len for MIX with SECT")
       end select
    case(G_Hexcel)
       cell%ok(n_tot) = cell%ok(n_sv) .and. cell%ok(n_mix) &
            & .and. cell%ok(n_radius) .and. cell%ok(n_side) &
            & .and. .not. &
            & ( cell%ok(n_meshx) .or. cell%ok(n_meshy) &
            & .or. cell%ok(n_splitx) .or. cell%ok(n_splity) )
       select case(sectori)
       case(S_not)
           if(size(cell%mix)/=size(cell%radius)) &
               call XABORT("G2S: wrong len for MIX")
       case(S_X_tot)
          if(size(cell%mix)/=6*size(cell%radius)-5*sectorj) &
               call XABORT("G2S: wrong len for MIX with SECT")
       end select
    case default
       call XABORT("G2S: the type of geometry is not allowed")
    end select
  end subroutine verrifieCB

  subroutine exploiteSplit(cell)
    type(t_celluleBase),intent(inout) :: cell
    !modification des vecteurs radius, mix et merge si splitr active

    double precision,dimension(:),allocatable :: tmpRad
    integer,dimension(:),allocatable          :: tmpMix,tmpMrg
    double precision                          :: interval
    integer                                   :: longueur,i,j,k,sect,valS,sectori,sectorj,typecel,longueur2, &
                                                 k2,nsectext,nsectint,ip,nx,ny

    if(.not. cell%ok(n_splitr)) return
    if((size(cell%splitr)+1)/=size(cell%radius)) call XABORT("G2S : wrong length for SPLITR")
    longueur=0
    sectori = cell%sv(14)
    sectorj = cell%sv(15)
    typecel = cell%sv(1)
    do i = 1,size(cell%splitr)
       longueur=longueur+abs(cell%splitr(i))
    end do
    nx=0
    ny=0
    nsectint=0
    nsectext=0
    select case(typecel)
    case(G_Carcel)
       nx=size(cell%meshx)-1
       ny=size(cell%meshy)-1
       if((nx*ny.gt.1).and.sectori.ne.S_not) call XABORT("G2S : CARCEL lr lx>1 ly>1 is incompatinle with SECT")
       select case(sectori)
       case(S_not)
          nsectint=1
          nsectext=1
       case(S_X_tot,S_T_tot)
          nsectint=4
          if(sectorj /= 0) nsectint=1
          nsectext=4
       case(S_TX_tot,S_TXS_tot)
          nsectint=8
          if(sectorj /= 0) nsectint=1
          nsectext=8
       case(S_WM_tot)
          nsectint=12
          if(sectorj /= 0) nsectint=1
          nsectext=12
       end select
    case(G_Hexcel)
       select case(sectori)
       case(S_not)
          nsectint=1
          nsectext=1
       case(S_X_tot)
          nsectint=6
          if(sectorj /= 0) nsectint=1
          nsectext=6
       end select
    end select
    if(nx*ny.gt.1) then
       longueur2=nx*ny*(longueur+1)
    else
       longueur2=nsectint*longueur+nsectext
    endif
    allocate(tmpRad(longueur+1),tmpMix(longueur2),tmpMrg(longueur2),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: exploitSplit(1) => allocation pb")
    tmpRad(1)=0.d0
    k = 1
    k2= 1
    do i = 1,size(cell%splitr)
       valS = cell%splitr(i)
       if(valS==0) then
          call XABORT("G2S : SPLITR may not be null")
       else if(valS>0) then
          !on coupe en sous-rayons egaux -> moyenne arthmetique des rayons
          interval = (cell%radius(i+1)-cell%radius(i))/valS
          do j = 1,valS
             if(nx*ny.gt.1) then
                do ip=1,nx*ny
                   tmpMix((ip-1)*(longueur+1)+k) = cell%mix((ip-1)*(size(cell%radius))+i)
                   tmpMrg((ip-1)*(longueur+1)+k) = (ip-1)*(longueur+1)+k
                end do
             else
                do sect=1,nsectint
                   tmpMix(k2) = cell%mix(nsectint*(i-1)+sect)
                   !! BEFORE tmpMrg(k) = cell%merge(i)
                   tmpMrg(k2) = k2
                   k2= k2+ 1
                end do
             endif
             k = k + 1
             tmpRad(k) = tmpRad(k-1) + interval
          end do
       else
          !on coupe en sous-surfaces egales -> moyenne geometrique des rayons
          valS = abs(valS)
          interval = (cell%radius(i+1)**2-cell%radius(i)**2)/valS
          do j = 1,valS
             if(nx*ny.gt.1) then
                do ip=1,nx*ny
                   tmpMix((ip-1)*(longueur+1)+k) = cell%mix((ip-1)*(size(cell%radius))+i)
                   tmpMrg(k) = (ip-1)*(longueur+1)+k
                end do
             else
                do sect=1,nsectint
                   tmpMix(k2) = cell%mix(nsectint*(i-1)+sect)
                   !! BEFORE tmpMrg(k) = cell%merge(i)
                   tmpMrg(k2) = k2
                   k2= k2+ 1
                end do
             endif
             k = k + 1
             tmpRad(k) = sqrt(tmpRad(k-1)**2 + interval)
          end do
       end if
    end do
    if(nx*ny.gt.1) then
       do ip=1,nx*ny
          tmpMix(ip*(longueur+1)) = cell%mix(ip*(size(cell%radius)))
          tmpMrg(ip*(longueur+1)) = ip*(longueur+1)
       end do
    else
       do sect=1,nsectext
          tmpMix(k2) = cell%mix(nsectint*size(cell%splitr)+sect)
          !!BEFORE tmpMrg(longueur+1) = cell%merge(size(cell%merge))
          tmpMrg(k2) = k2
          k2= k2+ 1
       end do
    endif
    deallocate(cell%radius) ; allocate(cell%radius(longueur+1),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: exploiteSplit(2) => allocation pb")
    deallocate(cell%mix) ; allocate(cell%mix(longueur2),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: exploiteSplit(3) => allocation pb")
    deallocate(cell%merge) ; allocate(cell%merge(longueur2),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: exploiteSplit(4) => allocation pb")
    cell%radius(:longueur+1) = tmpRad(:longueur+1)
    cell%mix(:longueur2) = tmpMix(:longueur2) ; cell%merge(:longueur2) = tmpMrg(:longueur2)
    deallocate(tmpRad,tmpMix,tmpMrg)
  end subroutine exploiteSplit
  
  recursive subroutine buildCellsBase(ip,sz,dirname)
    type(c_ptr),intent(inout):: ip
    integer,intent(inout)    :: sz
    character*12 ,intent(in) :: dirname

    character*12             :: namp, savename, subdirname
    integer                  :: type, long

    integer,parameter        :: dir=0

    call writeCellBase(ip,sz,dirname)
    namp = ' '
    call LCMNXT(ip,namp)
    savename = namp
    do 
       if(namp == ' ') exit
       call LCMLEN(ip,namp,long,type)
       if(type == dir) then
          if(namp /= 'BIHET') then
             subdirname = namp
             call LCMSIX(ip,namp,1)
             call buildCellsBase(ip,sz,subdirname)
             call LCMSIX(ip,namp,2)
          end if
       end if
       call LCMNXT(ip,namp)
       if(namp == savename) exit
    end do
  end subroutine buildCellsBase

  subroutine writeCellBase(ip,sz,dirname)
    type(c_ptr),intent(in)   :: ip
    integer,intent(inout)    :: sz
    character*12 ,intent(in) :: dirname

    integer,dimension(40)    :: st !state vector
    integer                  :: type, long, i
    logical                  :: toCreate

    call LCMLEN(ip,'STATE-VECTOR',long,type)
    if(long==0) then
       !on est dans la partie nouvelle (donnees ajoutees par pretaitement
       !python => on sort de la subroutine
       return
    end if
    call LCMGET(ip,'STATE-VECTOR',st)
    call LCMLEN(ip,'NPIN        ',long,type) ! pour tester si c'est un cluster
    if(st(8)==0 .and. long==0) then !pas de sous-cellules
       toCreate = .true.
       do i = 1,sz
          if(tabCelluleBase(i)%name==dirname) toCreate = .false.
       end do
       if(toCreate) then
          sz = sz + 1
!!$          allocate(tabCelluleBase(sz)%p,stat=alloc_ok)
!!$          if (alloc_ok /= 0) call XABORT("G2S: writeCellBase => allocation pb")
          call createCB(tabCelluleBase(sz),dirname,ip)
       end if
    end if
  end subroutine writeCellBase

end module celluleBase
