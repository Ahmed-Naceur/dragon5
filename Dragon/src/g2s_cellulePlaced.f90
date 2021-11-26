!
!-----------------------------------------------------------------------
!
!Purpose:
! Creation of an array of type(t_cellulePlaced) structures.
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
! Les structures definies de type "t_cellulePlaced" possedent comme champs :
!  -une reference a une cellule de base (indice dans le tableau celluleBase)
!  -la position (x,y) du centre de la cellule
!  -l'orientation de la cellule (turn). La definition de ce turn depend du
!   type de la geometrie envisagee (les turns n'ont pas la meme signification)
! Le tableau est cree par la routine chargerNewData, qui en plus  complete
! les donnees des cellules de base. Pour cela, elle va lire les donnees crees 
! par les routines de pretraitement.
! \\\\
! variable globale:
!  - tabCellulePlaced : tableau des cellules placees
!  - geomTyp : type de la geometrie exterieure
! \\\\
! fonctions:
!  - initializeTabCellulePlaced : mise a zero du tableau
!  - destroyTabCellulePlaced : liberation de la memoire
!  - splitCells : eclatement des cellules en elements geometriques simples
!  - chargerNewData : recuperations des donnees supplementaires crees a
!                     l'etape de pretraitement
!  - litDonneesSup : lectures des donnees supplementaires
!
!-----------------------------------------------------------------------
!
module cellulePlaced
    use cast
    use celluleBase
    use constType
    use construire
    use GANLIB
    use segArc
    
    implicit none

  !cellule "prete a l'emploi" i.e une reference a une celluleBase
  !plus la position du centre de la cellule, et son orientation
  !par rapport a la celluleBase de reference
  type t_cellulePlaced
     integer          :: indice  !indice de la celluleBase dans le tableau
     double precision :: xcenter !abscice du centre
     double precision :: ycenter !ordonnee du centre
     integer          :: turn    !rotation par rapport a la celluleBase
     integer,dimension(:), allocatable :: gig !gigogne de la cellule
     integer,dimension(:), allocatable :: mrg !gigogne equivalente de la cellule
  end type t_cellulePlaced

  integer,parameter :: dimTabCellulePlaced=10000

  !variable globale de type tableau de cellulePlaced
  type(t_cellulePlaced),dimension(:), allocatable :: tabCellulePlaced

  !variable globale donnant le type de la geometrie envisagee
  ! (rectangle, hexagonale, triangulaire, tubulaire)
  integer,save :: geomTyp

  integer,parameter :: RecTyp=1 , HexTyp=2 , TriaTyp=3 , TubeTyp=4

contains

  subroutine initializeTabCellulePlaced()

    allocate(tabCellulePlaced(dimTabCellulePlaced),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: initializeTabCellulePlaced => allocation pb")
  end subroutine initializeTabCellulePlaced

  subroutine destroyTabCellulePlaced(szP)
    integer,intent(in) :: szP
    integer :: i

    do i = 1,szP
       deallocate(tabCellulePlaced(i)%gig,tabCellulePlaced(i)%mrg)
    end do
    deallocate(tabCellulePlaced)
   end subroutine destroyTabCellulePlaced

  subroutine splitCells(szP,szSA)
    integer,intent(in)    :: szP
    integer,intent(inout) :: szSA
    
    integer               :: i,j,s,keepSzSA,ip,ix,iy,lmx,lmy
    double precision      :: sx,sy,sxt,syt,offcx,offcy,cx,cy,cxloc,cyloc,offcxt,offcyt
    type(t_cellulePlaced) :: tcp
    type(t_celluleBase)   :: tcb
    integer,dimension(:),allocatable :: mix,psplitx,psplity

    !nullify(mix)
    if(szP > dimTabCellulePlaced) call XABORT('splitCells: dimTabCellulePlaced overflow.')
    do i = 1,szP
       tcp = tabCellulePlaced(i)
       tcb = tabCelluleBase(tcp%indice)
       keepSzSA = szSA
       select case(tcb%sv(1))
       case(G_Car2d)
          sx = tcb%meshx(size(tcb%meshx))
          sy = tcb%meshy(size(tcb%meshy))
          s = size(tcb%mix)
          allocate(mix(s))
          mix = (/(j,j=1,s)/)
          call construit_car2d(tcp%xcenter,tcp%ycenter,sx,sy,tcp%turn,&
               & tcb%meshx/sx,tcb%meshy/sy,tcb%splitx,tcb%splity,mix,szSA)
          deallocate(mix)
          !nullify(mix)
       case(G_Carcel)
          sxt = tcb%meshx(size(tcb%meshx))
          syt = tcb%meshy(size(tcb%meshy))
          s = size(tcb%radius)
          lmx = size(tcb%meshx)
          lmy = size(tcb%meshy)
          if ((lmx.eq.2).and.(lmy.eq.2)) then
             allocate(mix(s),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: splitCells(2) => allocation pb")
             mix = (/(j,j=1,s)/)
             call construit_carcel(tcp%xcenter,tcp%ycenter,sxt,syt,tcp%turn,&
                  & tcb%radius,tcb%offcenter(1),tcb%offcenter(2),tcb%splitx,&
                  tcb%splity,mix,tcb%sv(14),tcb%sv(15),tcb%cluster,szSA)
             deallocate(mix)
            ! nullify(mix)
         else
!         AFTER : CARCEL lr lx ly with lx>0 ly>0       
             ip = 0
             allocate(mix(s),psplitx(lmx-1),psplity(lmy-1),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: splitCells(3) => allocation pb")
             do ix = 2, lmx
                do iy = 2, lmy
                   ip = ip+1
                   cxloc = 0.5D0*(tcb%meshx(ix)+tcb%meshx(ix-1))
                   cyloc = 0.5D0*(tcb%meshy(iy)+tcb%meshy(iy-1))
                   sx = tcb%meshx(ix)-tcb%meshx(ix-1)
                   sy = tcb%meshy(iy)-tcb%meshy(iy-1)
                   offcx = 0.5D0*sxt-cxloc
                   offcy = 0.5D0*syt-cyloc
                   cx = tcp%xcenter-offcx
                   cy = tcp%ycenter-offcy
                   offcxt = offcx+tcb%offcenter(1)
                   offcyt = offcy+tcb%offcenter(2)
                   psplitx = (/(tcb%splitx(j),j=ix-1,lmx-1)/)
                   psplity = (/(tcb%splity(j),j=iy-1,lmy-1)/)
                   mix = (/(j,j=(ip-1)*s+1,ip*s)/)
                   call construit_carcel(cx,cy,sx,sy,tcp%turn,&
                        & tcb%radius,offcxt,offcyt,psplitx,&
                        & psplity,mix,tcb%sv(14),tcb%sv(15),tcb%cluster,szSA)
                end do
             end do
             deallocate(mix,psplitx,psplity)
             !nullify(mix,psplitx,psplity)
          endif
       case(G_Hex)
          call construit_hexhom(tcp%xcenter,tcp%ycenter,tcb%side,1,szSA,tcb%sv(14))
       case(G_Hexcel)
          s = size(tcb%radius)
          allocate(mix(s),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: splitCells(4) => allocation pb")
          mix = (/(j,j=1,s)/)
          call construit_hexcel(tcp%xcenter,tcp%ycenter,tcb%side,tcp%turn,&
               & tcb%radius,tcb%offcenter(1),tcb%offcenter(2),mix,    &
               & tcb%sv(14),tcb%sv(15),tcb%cluster,szSA)
          deallocate(mix)
          !nullify(mix)
       case(G_Tri)
          allocate(mix(1))
          mix = (/1/)
          call construit_tri2d(tcp%xcenter,tcp%ycenter,tcb%side,tcp%turn,&
               & tcb%sv(3),mix,szSA)
          deallocate(mix)
          !nullify(mix)
       case(G_Tube)
          s = size(tcb%mix)
          allocate(mix(s+1),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: splitCells(5) => allocation pb")
          mix(1:s) = (/(j,j=1,s)/)
          mix(s+1) = fooMix
          call construit_tube(0.d0,0.d0,tcb%radius,mix,tcb%cluster,szSA)
          deallocate(mix)
          !nullify(mix)
       case default
          call XABORT("G2S: splitCells --> Type of geometry not supported")
       end select
       !ajout du numero de la cellulePlaced dont sont issus les segArcs 
       do j = keepSzSA+1,szSA
          tabSegArc(j)%indCellPg = i
          tabSegArc(j)%indCellPd = i
       end do
    end do
  end subroutine splitCells
  
  !en sortie, toutes les cellules de base ont tous leurs
  !champs remplis
  subroutine chargerNewData(geoIp,szB,szP)
    type(c_ptr),intent(in)    :: geoIp
    integer,intent(inout) :: szB,szP

    integer  :: i
    type(c_ptr) :: ip

    ip = geoIp
    call LCMSIX(ip,'NEW-DATA    ',1)
    do i = 1,szB
       call litDonneesSup(ip,tabCelluleBase(i),i,szP)
    end do
    call LCMSIX(ip,'NEW-DATA    ',2)
  end subroutine chargerNewData

  subroutine litDonneesSup(ip,cellB,ind,szP)
    type(c_ptr),intent(inout)         :: ip
    integer,intent(inout)             :: szP 
    integer,intent(in)                :: ind
    type(t_celluleBase),intent(inout) :: cellB

    integer                           :: i,lg,typ,dimGig
    character*12                      :: posName,number,mrgName
    real,dimension(2)                 :: sidexy
    real,dimension(:),allocatable     :: cx,cy
    integer,dimension(:),allocatable  :: tu


    if (cellB%name/='/           ') then
       call LCMSIX(ip,cellB%name,1) !on entre dans le repertoire
       !lecture des dimension de la cellule
       call LCMLEN(ip,'SIDEXY      ',lg,typ)
       if (lg==0) then !pas de donnees supplementaires
          call LCMSIX(ip,' ',2)
          return
       else if (lg==2) then !c'est un rectangle
          call LCMGET(ip,'SIDEXY      ',sidexy)
          if (.not. cellB%ok(n_meshx)) then
             allocate(cellB%meshx(2))
             cellB%meshx(1) = 0.d0
             cellB%meshx(2) = sidexy(1)
          endif
          if (.not. cellB%ok(n_meshy)) then
             allocate(cellB%meshy(2))
             cellB%meshy(1) = 0.d0
             cellB%meshy(2) = sidexy(2)
          endif
       else !c'est un triangle ou un hexagone
          if (.not. cellB%ok(n_side)) then
             call LCMGET(ip,'SIDEXY      ',sidexy)
             cellB%side = sidexy(1)
          end if
       end if

       call LCMLEN(ip,'TURN        ',lg,typ)
       allocate(cx(lg),cy(lg),tu(lg),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: litDonneesSup(2) => allocation pb")
       call LCMGET(ip,'COORDX      ',cx)
       call LCMGET(ip,'COORDY      ',cy)
       call LCMGET(ip,'TURN        ',tu)
       do i = 1,lg
          szP = szP + 1
          tabCellulePlaced(szP)%indice = ind
          tabCellulePlaced(szP)%xcenter = cx(i)
          tabCellulePlaced(szP)%ycenter = cy(i)
          tabCellulePlaced(szP)%turn = tu(i)
          number = i2s(i)
          posName = 'POS' // number(:9)
          call LCMLEN(ip,posName,dimGig,typ)
          allocate(tabCellulePlaced(szP)%gig(dimGig),stat=alloc_ok)
          if (alloc_ok /= 0) then
             write(6,*) "litDonneesSup: szP=",szP," dimGig=",dimGig
             call XABORT("G2S: litDonneesSup(4) => allocation pb")
          endif
          call LCMGET(ip,posName,tabCellulePlaced(szP)%gig)
          mrgName = 'MRG' // number(:9)
          allocate(tabCellulePlaced(szP)%mrg(dimGig),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: litDonneesSup(5) => allocation pb")
          call LCMGET(ip,mrgName,tabCellulePlaced(szP)%mrg)
       end do
       deallocate(cx,cy,tu)
       !on sort du repertoire
       call LCMSIX(ip,cellB%name,2)
    else !une seule cellule
       szP = 1
       tabCellulePlaced(szP)%indice = 1
       tabCellulePlaced(szP)%xcenter = 0.d0
       tabCellulePlaced(szP)%ycenter = 0.d0
       tabCellulePlaced(szP)%turn = 1
       allocate (tabCellulePlaced(szP)%gig(1))
       tabCellulePlaced(szP)%gig = 1
       allocate (tabCellulePlaced(szP)%mrg(1))
       tabCellulePlaced(szP)%mrg = 1
    end if
  end subroutine litDonneesSup

end module cellulePlaced
