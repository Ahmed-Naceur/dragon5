!
!-----------------------------------------------------------------------
!
!Purpose:
! Generate the region numerotation of the geometry.
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
! La numerotation est utilisee pour le jeux de donnees SALOME. Il faut donner un
! numero a chaque region definie par des arcs, des segments et/ou des cercles.
! Les limites exterieures de la geometrie sont caracterisees par un indice zero
! pour le node.
! \\\\
! Le principe de l'algorithme utilise est le suivant :
!  - creation d'un tableau comportant tous les points de jonctions entre des
!    elements geometriques
!  - rangement des elements geometriques concernes dans le tableau, par ordre
!    croissant d'angle d'incidence de leur tangente au point concerne
!  - attribution arbitraire d'un numero de domaine entre les elements de ce
!    tableau, avec verification de coherence, et decalage d'indice dans le
!    cas de regions deja numerotees
!  - prise en compte des regions annulaires (caracterisees par aucune
!    intersection avec d'autre elements geometriques), et numerotation
!    arbitraire de ces anneaux
!  - renumerotation globale de l'ensemble des regions selon l'ordre suivant :
!    ymin -> ymax ; xmin -> xmax ; amax -> amin (angle d'incidence de la
!    tangente a la region au point considere). La renumerotation se fait
!    suivant l'ordre des point minimaux de chacunes des regions considerees.
! \\\\
! Finalement on obtient une numerotation de ce type :
! \\\\     ________
!     | 4 |   /|
!     |___|2 / |
!     | 1 | / 3|
!     |___|/___|
! \\\\
! variables globales:
!  - tabPtNode : tableau des nodes non-annulaires
!  - tabCercNode : tableau des nodes annulaires
!  - tabPtMN : tableau des points minimaux de chaque node
! \\\\
! fonctions:
!  - createNodes : fonction d'entree du module
!  - associatePoints : creation du tableau tabPtNode
!  - addInLine : ajout d'un element dans le tableau tabPtNode
!  - createPoint : construction d'un point
!  - calculAngleDep : calcul l'angle d'incidence de depart ou d'arrivee d'un
!                     element en un point
!  - isEqualPt : test de l'egalite de deux points
!  - solveNodeSystem : resolution du systeme cree avec tabPtNode
!  - XetNodeAY : trouve ou positionne la valeur du node avant ou apres une
!                valeur d'angle d'incidence (X = g ou s ; Y = v ou p)
!  - getMixAY : trouve la valeur du mix avant ou apres une valeur d'angle
!               d'incidence (Y = v ou p)
!  - remplaceNodeEtDecale : decalage des valeurs de node affectees dans le cas
!                           ou une mauvaise valeur a ete preaffectee (cf algo)
!  - circularNodes : traite les nodes annulaires
!  - addCircel : ajout d'un cercle dans le tableau tabCercNode
!  - solveCircularNodes : resolution du systeme cree par tabCercNode
!  - disproj : calcule la distance entre un point et un segment (pour
!              determiner le milieux exterieur de l'anneau maximum)
!  - lessThanPtNd : fonction d'ordre sur les t_ptMinNode
!  - renumNodes : renumerotation des nodes dans l'ordre lexicographique
!  - givePtMin : donne le point minimal d'un node
!  - sortTabPtMN : trie le tableau tabPtMN
!  - getGigogneData : recupere les donnees relatives a la gigogne de provenance
!                     du node
!  - setIntIfPos : affectation d'une valeur entiere a une autre si elle est positive,
!                  et test de coherence si besoin est
!
!-----------------------------------------------------------------------
!
module ptNodes
  use boundCond
  use celluleBase
  use cellulePlaced
  use constUtiles
  use segArc

  implicit none

  !pour les segments et les arcs
  type t_ptNode
     type(t_point)                         :: position        !position
     integer,dimension(:), allocatable     :: listOfSA        !indice dans le
        !tableau des SA
     logical,dimension(:), allocatable     :: listOfDeparture !booleen
        !correspondant
     double precision,dimension(:),allocatable :: listOfDirection !angle de depart
        !du SA (eventuellement decale de +/-epsilon pour les pb de tangence)
     integer                               :: sz              !taille des
        !tableaux
  end type t_ptNode

  type(t_ptNode),dimension(:),allocatable,save :: tabPtNode

  !pour les cercles
  type t_cercNode
     type(t_point)                         :: centre       !position
     logical                               :: withSect     !dit si la cellule
        !envisagee est sectorisee
     integer,dimension(:),allocatable      :: listOfSA     !indice dans le
        !tableau des SA
     double precision,dimension(:),allocatable :: listOfRadius !rayons
     integer                               :: imin         !position du dernier cercle
     integer                               :: sz           !taille des tableaux
     integer                               :: indNodeExt   !numero de l'element
        !donnant le node exterieur
     logical                               :: coteGauche   !dit si le node
        !exterieur est celui de gauche de l'element designe par indNodeExt
  end type t_cercNode

  type(t_cercNode),dimension(:),allocatable,save :: tabCercNode

  type t_ptMinNode
     type(t_point)    :: ptMin   !point le plus bas et le plus a gauche
     double precision :: alfa    !angle de depart
     integer          :: indNode !numero du node avant renumerotation
  end type t_ptMinNode

  type(t_ptMinNode),dimension(:),allocatable,save :: tabPtMN

  type t_nodeGigSect
     integer  :: indTabCellPlac !indice de la cellulePlaced d'origine (gig)
     integer  :: ring           !indice de l'anneau ou se trouve le node
     integer  :: sect           !indice du secteur angulaire du node
     integer  :: neutronicMix   !milieux neutronique du node
     integer  :: merge          !merge du node
     integer  :: dragSector     !indice de type DRAGON du node
     logical  :: clust          !.true. si le node est un cercle de cluster
  end type t_nodeGigSect

  type(t_nodeGigSect),dimension(:),allocatable,save :: tabNodeGigSect

contains

  subroutine createNodes(szSA,dimTabCelluleBase,nbNode,merg)
    implicit none
    integer,intent(in)  :: szSA,dimTabCelluleBase
    integer,intent(out) :: nbNode,merg(dimTabCelluleBase)
    
    integer  :: i,nbPts,nbCers,nbNode_noclust

    !cas des Arcs et Segments du systeme :
    ! preparation du tableau

    allocate(tabPtNode(szSA*2),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: createNodes(1) => allocation pb")
    tabPtNode(1:szSA*2)%sz=0
    nbPts=0
    nbNode=0
    !creation d'un systeme comprenant une liste d'indices
    ! d'elements geometriques decrits dans le sens trigo
    ! tels que : chaque element de la liste represente
    ! les elements geometriques partant ou arrivant sur un
    ! point specifique (ainsi que leur sens depuis ce point :
    ! .true.=>depart , .false.=>arrivee)
    call associatePoints(szSA,nbPts)

    !resolution du systeme precedent
    call solveNodeSystem(szSA,nbPts,nbNode)

    ! nettoyage du tableau
    do i = 1,nbPts
       if(allocated(tabPtNode(i)%listOfSA))        &
            deallocate(tabPtNode(i)%listOfSA)
       if(allocated(tabPtNode(i)%listOfDeparture)) &
            deallocate(tabPtNode(i)%listOfDeparture)
       if(allocated(tabPtNode(i)%listOfDirection)) &
            deallocate(tabPtNode(i)%listOfDirection)
    enddo
    deallocate(tabPtNode)

    !prise en compte des cercles complets :
    ! preparation du tableau
    allocate(tabCercNode(szSA),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: createNodes(2) => allocation pb")
    do i = 1,szSA
       tabCercNode(i)%imin=0
       tabCercNode(i)%sz=0
       tabCercNode(i)%indNodeExt=-1
    enddo
    !on classe les cercles suivant leur centre, et par rayon croissant,
    ! en ajoutant dans la liste en derniere position un arc de cercle ou un
    ! segment le plus proche du plus grand cercle, pour donner une reference
    ! de Node
    nbCers = 0
    nbNode_noclust=nbNode
    call circularNodes(szSA,nbNode,nbCers)

    ! nettoyage du tableau
    do i = 1,nbCers
       if(allocated(tabCercNode(i)%listOfSA))      &
            deallocate(tabCercNode(i)%listOfSA)
       if(allocated(tabCercNode(i)%listOfRadius))  &
            deallocate(tabCercNode(i)%listOfRadius)
    enddo
    deallocate(tabCercNode)

    !renumerotation des nodes dans l'ordre lexicographique
    call renumNodes(szSA,nbNode,nbNode_noclust)

    !recuperation des informations sur les gigognes
    call getGigogneData(szSA,nbNode)
    if(nbNode.gt.dimTabCelluleBase) call XABORT('createNodes: merg overflow.')
    do i = 1,nbNode
      merg(i)=tabNodeGigSect(i)%merge
    enddo
    deallocate(tabNodeGigSect)
  end subroutine createNodes

  subroutine associatePoints(szSA,nbPts)
    implicit none
    integer,intent(in)    :: szSA
    integer,intent(inout) :: nbPts

    type(t_point)  :: pt
    type(t_segArc) :: sa
    integer        :: i,j,k
    logical        :: isNewPt,dep

    !programmation defensive
    integer :: taille_table_tabPtNode

    taille_table_tabPtNode = size(tabPtNode)
    
    !creation du tableau
    do i = 1,szSA
       sa = tabSegArc(i)
       if(sa%typ==tcer) cycle
       do k = 1,2
          dep=(k==1)
          pt = createPoint(i,dep)
          isNewPt = .true.
          do j = 1,nbPts
             isNewPt = .not.isEqualPt(pt,tabPtNode(j)%position)
             if(.not.isNewPt) then
                call addInLine(j,i,dep)
                exit
             endif
          enddo
          if(isNewPt) then
             nbPts = nbPts+1
             if (nbPts > taille_table_tabPtNode) &
                  call XABORT("G2S : memory problem in routine associatePoints")
             call addInLine(nbPts,i,dep)
          endif
       enddo
    enddo
 end subroutine associatePoints

  subroutine addInLine(indTab,indSA,dep)
    integer,intent(in) :: indTab,indSA
    logical,intent(in) :: dep
    !permet d'ajouter la reference a un segArc dans le tableau des noeuds
    ! (les angles des secteurs sont ordonnes croissants dans le sens trigo)
    integer,dimension(tabPtNode(indTab)%sz+1)          :: lSA
    logical,dimension(tabPtNode(indTab)%sz+1)          :: lDep
    double precision,dimension(tabPtNode(indTab)%sz+1) :: lDir
    integer                                   :: sz,i,j
    double precision                          :: angl
    logical                                   :: flag

    sz   = tabPtNode(indTab)%sz
    angl = calculAngleDep(indSA,dep)
    if(sz/=0) then
       flag = .true.
       do i = 1,sz
          if(angl > tabPtNode(indTab)%listOfDirection(i)) then
             lSA(i)  = tabPtNode(indTab)%listOfSA(i)
             lDep(i) = tabPtNode(indTab)%listOfDeparture(i)
             lDir(i) = tabPtNode(indTab)%listOfDirection(i)
          else
             lSA(i)  = indSa
             lDep(i) = dep
             lDir(i) = angl
             do j = i,sz
                lSA(j+1)  = tabPtNode(indTab)%listOfSA(j)
                lDep(j+1) = tabPtNode(indTab)%listOfDeparture(j)
                lDir(j+1) = tabPtNode(indTab)%listOfDirection(j)          
             enddo
             flag = .false.
             exit
          endif
       enddo
       if(flag) then
          lSA(sz+1)  = indSa
          lDep(sz+1) = dep
          lDir(sz+1) = angl
       endif
       deallocate( tabPtNode(indTab)%listOfSA &
               & , tabPtNode(indTab)%listOfDeparture &
               & , tabPtNode(indTab)%listOfDirection )
    else
       tabPtNode(indTab)%position = createPoint(indSA,dep)
       lSA(1)  = indSa
       lDep(1) = dep
       lDir(1) = angl
    endif
    allocate( tabPtNode(indTab)%listOfSA(sz+1)        &
          & , tabPtNode(indTab)%listOfDeparture(sz+1) &
          & , tabPtNode(indTab)%listOfDirection(sz+1) ,stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: addInLine(3) => allocation pb")
    tabPtNode(indTab)%listOfSA=lSA
    tabPtNode(indTab)%listOfDeparture=lDep
    tabPtNode(indTab)%listOfDirection=lDir
    tabPtNode(indTab)%sz=sz+1
  end subroutine addInLine

  function calculAngleDep(indSA,dep)
    integer,intent(in) :: indSA
    logical,intent(in) :: dep
    double precision   :: calculAngleDep

    type(t_segArc) :: sa

    calculAngleDep = 0.d0
    sa = tabSegArc(indSA)
    if(sa%typ==tseg) then
       if(dep) then
          calculAngleDep = calculeAngle(sa%x,sa%y,sa%dx,sa%dy)
       else
          calculAngleDep = calculeAngle(sa%dx,sa%dy,sa%x,sa%y)         
       endif
    else if(sa%typ==tarc) then
       !on modifie l'angle de +/- 5*epsilon pour bien prendre en compte
       ! les positions relatives avec les arcs dans les cas de tangence
       ! (le rapport de modification prend en compte le rayon de l'element
       ! pour gerer les cas de tangence entre cercles)
       if(dep) then
          calculAngleDep = angleNormal(sa%a+pi_2_c+muleps1*(2-tanh(sa%r))*epsilon)
       else
          calculAngleDep = angleNormal(sa%b-pi_2_c-muleps1*(2-tanh(sa%r))*epsilon)
       endif
    else
       call XABORT("G2S : internal error in function calculAngleDep")
    endif
  end function calculAngleDep

  function createPoint(indSA,dep)
    integer,intent(in) :: indSA
    logical,intent(in) :: dep
    type(t_point)      :: createPoint
    type(t_segArc)  :: sa

    createPoint%x=0.d0 ; createPoint%y=0.d0
    sa = tabSegArc(indSA)
    if(sa%typ==tseg) then
       if(dep) then
          createPoint%x=sa%x ; createPoint%y=sa%y
       else
          createPoint%x=sa%dx ; createPoint%y=sa%dy
       endif
    else if(sa%typ==tarc) then
       if(dep) then
          createPoint%x=sa%x+sa%r*cos(sa%a) ; createPoint%y=sa%y+sa%r*sin(sa%a)
       else
          createPoint%x=sa%x+sa%r*cos(sa%b) ; createPoint%y=sa%y+sa%r*sin(sa%b)
       endif
    else
       call XABORT("G2S : internal error in function createPoint")
    endif
  end function createPoint

  function isEqualPt(p1,p2)
    type(t_point),intent(in) :: p1,p2
    logical                  :: isEqualPt

    isEqualPt = isEqualConst(p1%x,p2%x) .and. isEqualConst(p1%y,p2%y)
  end function isEqualPt

  subroutine solveNodeSystem(szSA,nbPts,nbNode)
    integer,intent(in)  :: szSA,nbPts
    integer,intent(out) :: nbNode

    integer  :: i,j,indNodeAv,indNodeAp,mixAv,mixAp,sz
    ! correction of plain CAR2D bug by Alain Hebert (May 2016)
    integer  :: nbFile
    logical  :: isOpen
    logical,parameter :: drawMix = .true.
    real,parameter,dimension(2) :: zoomx = (/ 0.0, 1.0 /) ! no x zoom on postscript plot
    real,parameter,dimension(2) :: zoomy = (/ 0.0, 1.0 /) ! no y zoom on postscript plot

    nbNode = 0
    do i = 1,nbPts
       sz = tabPtNode(i)%sz
       do j = 1,sz
          mixAv=getMixAv(i,j) ;  mixAp=getMixAp(i,j)
          if(mixAv<0 .or. mixAp<0) then
             call setNodeAv(i,j,0) ; call setNodeAp(i,j,0)
             cycle
          endif
! correction of plain CAR2D bug by Alain Hebert (May 2016)
          if(mixAv/=mixAp) then
             nbFile = 50
             do
                nbFile = nbFile + 1
                inquire(nbFile,opened=isOpen)
                if(isOpen) cycle
               open(nbFile,file='errorMix.ps')
                exit
             enddo
             call drawSegArc(nbFile,szSA,.false.,drawMix,zoomx,zoomy)
             close(nbFile)
             write(*,*) 'i,j,mixAv,mixAp : ',i,j,mixAv,mixAp
             call XABORT("G2S: internal problem for mix values. See the file &
                  &errorMix.ps")
          endif
! CS-IV : fin de la mise en commentaires de Alain
          indNodeAv=getNodeAv(i,j) ; indNodeAp=getNodeAp(i,j)
          if(indNodeAv>0 .and. indNodeAp<0) then
             call setNodeAp(i,j,indNodeAv)
          else if(indNodeAp>0 .and. indNodeAv<0) then
             call setNodeAv(i,j,indNodeAp)
          else if(indNodeAp<0 .and. indNodeAv<0) then
             nbNode = nbNode + 1
             call setNodeAv(i,j,nbNode) ; call setNodeAp(i,j,nbNode)
          else if(indNodeAv/=indNodeAp) then
             !incoherence de numerotation (possible si plus de 3 cotes a
             ! un domaine par remplisage des coins opposes) => on decremente
             ! le nombre de nodes, et on rend la numerotation coherente
             ! (remplacement du plus grand des deux par le plus petit, et
             ! decalage de -1 des numeros superieurs au plus grand)
             nbNode = nbNode - 1
             call remplaceNodeEtDecale(max(indNodeAv,indNodeAp), &
                  & min(indNodeAv,indNodeAp),szSA)
          endif
       enddo
    enddo

  end subroutine solveNodeSystem

  function getNodeAv(indTPN,indLSA)
    integer,intent(in) :: indTPN,indLSA
    integer            :: getNodeAv
    integer :: ind

    if(indLSA==1) then
       ind = tabPtNode(indTPN)%sz
    else
       ind = indLSA-1
    endif
    if(tabPtNode(indTPN)%listOfDeparture(ind)) then
       getNodeAv=tabSegArc(tabPtNode(indTPN)%listOfSA(ind))%nodeg
    else
       getNodeAv=tabSegArc(tabPtNode(indTPN)%listOfSA(ind))%noded
    endif
  end function getNodeAv

  function getNodeAp(indTPN,indLSA)
    integer,intent(in) :: indTPN,indLSA
    integer            :: getNodeAp

    if(tabPtNode(indTPN)%listOfDeparture(indLSA)) then
       getNodeAp=tabSegArc(tabPtNode(indTPN)%listOfSA(indLSA))%noded
    else
       getNodeAp=tabSegArc(tabPtNode(indTPN)%listOfSA(indLSA))%nodeg
    endif
  end function getNodeAp

  subroutine setNodeAv(indTPN,indLSA,val)
    integer,intent(in) :: indTPN,indLSA,val

    integer :: ind

    if(indLSA==1) then
       ind = tabPtNode(indTPN)%sz
    else
       ind = indLSA-1
    endif
    if(tabPtNode(indTPN)%listOfDeparture(ind)) then
       tabSegArc(tabPtNode(indTPN)%listOfSA(ind))%nodeg=val
    else
       tabSegArc(tabPtNode(indTPN)%listOfSA(ind))%noded=val
    endif

  end subroutine setNodeAv

  subroutine setNodeAp(indTPN,indLSA,val)
    integer,intent(in) :: indTPN,indLSA,val

    if(tabPtNode(indTPN)%listOfDeparture(indLSA)) then
       tabSegArc(tabPtNode(indTPN)%listOfSA(indLSA))%noded=val
    else
       tabSegArc(tabPtNode(indTPN)%listOfSA(indLSA))%nodeg=val
    endif

  end subroutine setNodeAp

  function getMixAv(indTPN,indLSA)
    integer,intent(in) :: indTPN,indLSA
    integer            :: getMixAv

    integer :: ind

    if(indLSA==1) then
       ind = tabPtNode(indTPN)%sz
    else
       ind = indLSA-1
    endif
    if(tabPtNode(indTPN)%listOfDeparture(ind)) then
       getMixAv=tabSegArc(tabPtNode(indTPN)%listOfSA(ind))%mixg
    else
       getMixAv=tabSegArc(tabPtNode(indTPN)%listOfSA(ind))%mixd
    endif
  end function getMixAv

  function getMixAp(indTPN,indLSA)
    integer,intent(in) :: indTPN,indLSA
    integer            :: getMixAp

    if(tabPtNode(indTPN)%listOfDeparture(indLSA)) then
       getMixAp=tabSegArc(tabPtNode(indTPN)%listOfSA(indLSA))%mixd
    else
       getMixAp=tabSegArc(tabPtNode(indTPN)%listOfSA(indLSA))%mixg
    endif
  end function getMixAp

  subroutine remplaceNodeEtDecale(badNum,goodNum,szSA)
    integer,intent(in) :: badNum,goodNum,szSA

    integer :: i
    do i = 1,szSA
       if(tabSegArc(i)%nodeg == badNum) tabSegArc(i)%nodeg = goodNum
       if(tabSegArc(i)%noded == badNum) tabSegArc(i)%noded = goodNum
       if(tabSegArc(i)%nodeg > badNum) &
            & tabSegArc(i)%nodeg = tabSegArc(i)%nodeg - 1
       if(tabSegArc(i)%noded > badNum) &
            & tabSegArc(i)%noded = tabSegArc(i)%noded - 1
    enddo
  end subroutine remplaceNodeEtDecale

  subroutine circularNodes(szSA,nbNode,szTabCer)
    integer,intent(in)    :: szSA
    integer,intent(inout) :: nbNode,szTabCer

    type(t_point)  :: centre
    integer        :: i,j
    logical        :: isNewCentre

    !creation du systeme
    do i = 1,szSA
       if(tabSegArc(i)%typ == tseg) cycle
       centre%x=tabSegArc(i)%x ; centre%y=tabSegArc(i)%y
       isNewCentre = .true.
       do j = 1,szTabCer
          isNewCentre = .not.isEqualPt(centre,tabCercNode(j)%centre)
          if(.not.isNewCentre) then
             call addCircel(j,i)
             exit
          endif
       enddo
       if(isNewCentre) then
          szTabCer = szTabCer + 1
          tabCercNode(szTabCer)%centre = centre
          tabCercNode(szTabCer)%withSect = &
            (tabCelluleBase(tabCellulePlaced(tabSegArc(i)%indCellPg)%indice)%sv(14) &
             /= S_not)
          call addCircel(szTabCer,i)
       endif
    enddo
    !determination du node exterieur et resolution du syteme
    call solveCircularNodes(szTabCer,szSA,nbNode)
  end subroutine circularNodes

  subroutine addCircel(indTab,indSA)
    integer,intent(in)    :: indTab,indSA
    
    integer,dimension(:), allocatable         :: lSA
    double precision,dimension(:),allocatable :: lRad
    integer                               :: sz,i,j,imin
    double precision                      :: radius
    logical                               :: flag

    sz   = tabCercNode(indTab)%sz
    imin   = tabCercNode(indTab)%imin
    radius = tabSegArc(indSA)%r
    do i = 1,sz
       if(abs(radius-tabCercNode(indTab)%listOfRadius(i))<muleps1*epsilon) then
          if(tabCercNode(indTab)%listOfSA(i) == indSA) return
       endif
    enddo
    allocate(lSA(sz+1),lRad(sz+1),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: addCircel => allocation pb")
    if(sz/=0) then
       flag = .true.
       do i = 1,sz
          if(radius > tabCercNode(indTab)%listOfRadius(i)) then
             lSA(i)  = tabCercNode(indTab)%listOfSA(i)
             lRad(i) = tabCercNode(indTab)%listOfRadius(i)
          else
             lSA(i)  = indSa
             lRad(i) = radius
             do j = i,sz
                lSA(j+1)  = tabCercNode(indTab)%listOfSA(j)
                lRad(j+1) = tabCercNode(indTab)%listOfRadius(j)
             enddo
             flag = .false.
             exit
          endif
       enddo
       if(flag) then
          lSA(sz+1)  = indSa
          lRad(sz+1) = radius
       endif
       deallocate( tabCercNode(indTab)%listOfSA &
               & , tabCercNode(indTab)%listOfRadius )
    else
       lSA(1)  = indSa
       lRad(1) = radius
    endif
    allocate( tabCercNode(indTab)%listOfSA(sz+1)        &
          & , tabCercNode(indTab)%listOfRadius(sz+1) ,stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: addCircle1 => allocation pb")
    tabCercNode(indTab)%listOfSA(:sz+1) = lSA(:sz+1)
    tabCercNode(indTab)%listOfRadius(:sz+1) = lRad(:sz+1)
    tabCercNode(indTab)%sz=sz+1
    tabCercNode(indTab)%imin=imin
  end subroutine addCircel

  subroutine solveCircularNodes(szTabCer,szSA,nbNode)
    integer,intent(in)    :: szTabCer,szSA
    integer,intent(inout) :: nbNode

    integer          :: i,j,indCer,indNodeAv,indNodeAp
    double precision :: d,dist,angl,cx,cy,radius
    type(t_segArc)   :: sa,cer

    do j = 1,szTabCer
       indCer = 0
       do i=tabCercNode(j)%sz,1,-1
          tabCercNode(j)%imin = i
          indCer = tabCercNode(j)%listOfSA(i)
          radius = tabCercNode(j)%listOfRadius(i)
          if(tabSegArc(indCer)%typ == tcer) go to 10
       enddo
       tabCercNode(j)%imin = 0
       cycle
       10 cer = tabSegArc(indCer)
       d = infinity
       do i = 1,szSA
          if(i==indCer) cycle
          sa = tabSegArc(i)
          if(sa%typ == tseg) then
             dist = disproj(tabCercNode(j)%centre,sa) - radius
             if((dist < 0.0).or.(dist >= d)) cycle
             tabCercNode(j)%coteGauche = .not.estADroite(sa%x,sa%y,sa%dx,sa%dy,cer%x,cer%y)
          else if(sa%typ == tcer) then
             if(isEqualConst(sa%x,cer%x).and.isEqualConst(sa%y,cer%y)) cycle
             dist = sa%r - longVect(sa%x-cer%x,sa%y-cer%y) - radius
             if((dist < 0.0).or.(dist >= d)) cycle
             tabCercNode(j)%coteGauche = .true.
          else
             if(sa%b > sa%a) then
                angl=(sa%b+sa%a)*0.5
             else
                angl=(sa%b+sa%a)*0.5+pi_c
             endif
             cx = sa%x+cos(angl)*sa%r ; cy = sa%y+sin(angl)*sa%r
             dist = longVect(cer%x-cx,cer%y-cy) - radius
             if((dist < 0.0).or.(dist >= d)) cycle
             tabCercNode(j)%coteGauche = longVect(cer%x-sa%x,cer%y-sa%y) < sa%r
          endif
          d = dist
          tabCercNode(j)%indNodeExt = i
       enddo
    enddo
    !resolution du systeme
    do j = 1,szTabCer
       if(tabCercNode(j)%indNodeExt < 0) cycle
       if(tabSegArc(tabCercNode(j)%listOfSA(1))%nodeg==fooNode) then
          nbNode = nbNode + 1
          tabSegArc(tabCercNode(j)%listOfSA(1))%nodeg = nbNode
       endif
       if(tabCercNode(j)%imin > 0) then
          if(tabCercNode(j)%coteGauche) then
             tabSegArc(tabCercNode(j)%listOfSA(tabCercNode(j)%imin))%noded = &
                tabSegArc(tabCercNode(j)%indNodeExt)%nodeg
          else
             tabSegArc(tabCercNode(j)%listOfSA(tabCercNode(j)%imin))%noded = &
                tabSegArc(tabCercNode(j)%indNodeExt)%noded
          endif
       endif
       do i = 1,tabCercNode(j)%sz-1
          if(i == tabCercNode(j)%imin) cycle
          if(tabSegArc(tabCercNode(j)%listOfSA(i))%noded==fooNode .and. &
              tabSegArc(tabCercNode(j)%listOfSA(i+1))%nodeg==fooNode) then
             nbNode = nbNode + 1
             tabSegArc(tabCercNode(j)%listOfSA(i))%noded = nbNode
             tabSegArc(tabCercNode(j)%listOfSA(i+1))%nodeg = nbNode
          else if(tabSegArc(tabCercNode(j)%listOfSA(i))%noded==fooNode) then
             tabSegArc(tabCercNode(j)%listOfSA(i))%noded = &
                  tabSegArc(tabCercNode(j)%listOfSA(i+1))%nodeg
          else if(tabSegArc(tabCercNode(j)%listOfSA(i+1))%nodeg==fooNode)then
             tabSegArc(tabCercNode(j)%listOfSA(i+1))%nodeg = &
                  tabSegArc(tabCercNode(j)%listOfSA(i))%noded
          else if((tabSegArc(tabCercNode(j)%listOfSA(i))%noded /= &
               tabSegArc(tabCercNode(j)%listOfSA(i+1))%nodeg)) then
             !on fait le remplacement dans touts les cas si il n'y a pas
             !de sectorisation de la cellule, et seulement si l'element
             !geometrique est un cercle complet si il y a sectorisation
             if(.not.(tabSegArc(tabCercNode(j)%listOfSA(i))%typ/=tcer)) then
                nbNode = nbNode - 1
                indNodeAv = tabSegArc(tabCercNode(j)%listOfSA(i))%noded
                indNodeAp = tabSegArc(tabCercNode(j)%listOfSA(i+1))%nodeg
                call remplaceNodeEtDecale(max(indNodeAv,indNodeAp), &
                     & min(indNodeAv,indNodeAp),szSA)
             endif
          endif
       enddo
    enddo
    !validation du node exterieur
    do j = 1,szTabCer
      if(tabSegArc(tabCercNode(j)%listOfSA(tabCercNode(j)%sz))%noded == fooNode) then
          if(tabCercNode(j)%indNodeExt>0) then
             if(tabCercNode(j)%coteGauche) then
                tabSegArc(tabCercNode(j)%listOfSA(tabCercNode(j)%sz))%noded &
                     = tabSegArc(tabCercNode(j)%indNodeExt)%nodeg
             else
                tabSegArc(tabCercNode(j)%listOfSA(tabCercNode(j)%sz))%noded &
                     = tabSegArc(tabCercNode(j)%indNodeExt)%noded
             endif
          else
             tabSegArc(tabCercNode(j)%listOfSA(tabCercNode(j)%sz))%noded = 0
          endif
       endif
    enddo
    
  end subroutine solveCircularNodes

  function disproj(pt,sg)
    type(t_point),intent(in)  :: pt
    type(t_segArc),intent(in) :: sg
    double precision          :: disproj

    double precision  :: prjx,prjy

    disproj = distance(pt%x,pt%y,sg%x,sg%y,sg%dx,sg%dy,prjx,prjy)
    if(isIn(prjx,prjy,sg)==0) then
       disproj = min(longVect(sg%x-pt%x,sg%y-pt%y),&
            longVect(sg%dx-pt%x,sg%dy-pt%y)) + muleps2*epsilon 
    endif
    !pour privilegier les arcs par rapport aux segments
    disproj = disproj + muleps2*epsilon
  end function disproj

  function lessThanPtNd(first,second)
    type(t_ptMinNode),intent(in) :: first,second
    logical                      :: lessThanPtNd

    if(.not.(isEqualConst(first%ptMin%y,second%ptMin%y))) then
       lessThanPtNd = (first%ptMin%y < second%ptMin%y)
    else if(.not.(isEqualConst(first%ptMin%x,second%ptMin%x))) then
       lessThanPtNd = (first%ptMin%x < second%ptMin%x)
    else
       lessThanPtNd = (first%alfa > second%alfa) !ordre inverse sur les angles
    endif
  end function lessThanPtNd

  subroutine renumNodes(szSA,nbNode,nbNode_noclust)
    integer,intent(in) :: szSA,nbNode,nbNode_noclust

    type(t_ptMinNode) :: tmpPtMN
    type(t_segArc)    :: sa
    type(t_point)     :: toOrig
    integer           :: i,j,numNod
    integer,dimension(:),allocatable :: newNumNode

    allocate(tabPtMN(nbNode),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: renumNodes(1) => allocation pb")
    !initialisation du tableau
    do i = 1,nbNode
       tabPtMN(i)%ptMin%y = infinity
       tabPtMN(i)%ptMin%x = infinity
       tabPtMN(i)%alfa    = 0.d0
       tabPtMN(i)%indNode = i
    enddo

    !remplisage du tableau
    do i = 1,szSA
       sa = tabSegArc(i)
       tmpPtMN = givePtMin(sa)
       numNod = sa%nodeg !travail sur le node gauche
       do j = 1,2
          if(numNod>0) then
             if(lessThanPtNd(tmpPtMN,tabPtMN(numNod))) then
                tabPtMN(numNod)%ptMin = tmpPtMN%ptMin
                tabPtMN(numNod)%alfa  = tmpPtMN%alfa
             endif
          endif
          numNod = sa%noded !travail sur le node droit
       enddo
    enddo

    !triage du tableau
    call sortTabPtMN(nbNode)

    !recuperation du vecteur de translation des numeros de nodes
    allocate(newNumNode(nbNode),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: renumNodes(2) => allocation pb")
    do i = 1,nbNode
       newNumNode(tabPtMN(i)%indNode)=i
    enddo
    toOrig = tabPtMN(1)%ptMin

    !renumerotation des nodes
    do i = 1,szSA
       numNod = tabSegArc(i)%nodeg
       if(numNod>0) then
         tabSegArc(i)%nodeg = newNumNode(numNod)
         tabSegArc(i)%clusg = (numNod.gt.nbNode_noclust)
       endif
       numNod = tabSegArc(i)%noded
       if(numNod>0) then
         tabSegArc(i)%noded = newNumNode(numNod)
         tabSegArc(i)%clusd = (numNod.gt.nbNode_noclust)
       endif
    enddo

    deallocate(tabPtMN,newNumNode)

    !translation vers l'origine des elements
    do i = 1,szSA
       tabSegArc(i)%x=tabSegArc(i)%x - toOrig%x
       tabSegArc(i)%y=tabSegArc(i)%y - toOrig%y
       if(tabSegArc(i)%typ==tseg) then
          tabSegArc(i)%dx=tabSegArc(i)%dx - toOrig%x
          tabSegArc(i)%dy=tabSegArc(i)%dy - toOrig%y
       endif
    enddo
    !ecriture dans les donnees de condition aux limites du vecteur de
    !translation
    bCData%toOrig_xy(1) = toOrig%x ; bCData%toOrig_xy(2) = toOrig%y
  end subroutine renumNodes

  function givePtMin(sa)
    type(t_segArc),intent(in) :: sa
    type(t_ptMinNode)         :: givePtMin

    double precision :: ax,ay,bx,by

    givePtMin%ptMin%x=infinity
    givePtMin%ptMin%y=infinity
    givePtMin%alfa=0.d0
    givePtMin%indNode=-1 !ne doit pas servir
    select case(sa%typ)
    case(tseg)
       if(isEqualConst(sa%y,sa%dy)) then
          if(sa%x<sa%dx) then
             givePtMin%ptMin%x = sa%x  ; givePtMin%ptMin%y = sa%y
             givePtMin%alfa = calculeAngle(sa%x,sa%y,sa%dx,sa%dy)
          else
             givePtMin%ptMin%x = sa%dx ; givePtMin%ptMin%y = sa%dy
             givePtMin%alfa = calculeAngle(sa%dx,sa%dy,sa%x,sa%y)
          endif
       else if(sa%y<sa%dy) then
          givePtMin%ptMin%x = sa%x  ; givePtMin%ptMin%y = sa%y
          givePtMin%alfa = calculeAngle(sa%x,sa%y,sa%dx,sa%dy)
       else
          givePtMin%ptMin%x = sa%dx ; givePtMin%ptMin%y = sa%dy
          givePtMin%alfa = calculeAngle(sa%dx,sa%dy,sa%x,sa%y)
       endif
    case(tcer)
       givePtMin%ptMin%x = sa%x ; givePtMin%ptMin%y = sa%y - sa%r
       givePtMin%alfa = pi_c
    case(tarc)
       if(isAngleInArc(-pi_2_c,sa)/=0) then
          givePtMin%ptMin%x = sa%x ; givePtMin%ptMin%y = sa%y - sa%r
          givePtMin%alfa = pi_c
       else
          call extremitesArc(sa,ax,ay,bx,by)
          if(isEqualConst(ay,by)) then
             if(ax<bx) then
                givePtMin%ptMin%x = ax ; givePtMin%ptMin%y = ay
                givePtMin%alfa = sa%a+pi_2_c+muleps1*epsilon
             else
                givePtMin%ptMin%x = bx ; givePtMin%ptMin%y = by
                givePtMin%alfa = sa%b-pi_2_c-muleps1*epsilon
             endif
          else if(ay<by) then
             givePtMin%ptMin%x = ax ; givePtMin%ptMin%y = ay
             givePtMin%alfa = sa%a+pi_2_c+muleps1*epsilon
          else
             givePtMin%ptMin%x = bx ; givePtMin%ptMin%y = by
             givePtMin%alfa = sa%b-pi_2_c-muleps1*epsilon
          endif
       endif
    case default
       call XABORT("G2S: internal error in function givePtMin")
    end select
  end function givePtMin

  subroutine sortTabPtMN(sz)
    implicit none
    integer,intent(in) :: sz
    !triage du tableau de pointsMin dans l'odre lexicographique
    !(utilisation d'un bubble-sort)
    type(t_ptMinNode) :: tmp
    integer           :: indMax,i
    logical           :: trie

    do indMax = sz,2,-1
       trie = .true.
       do i = 1,indMax-1
          if(lessThanPtNd(tabPtMN(i+1),tabPtMN(i))) then
             tmp=tabPtMN(i+1) ; tabPtMN(i+1)=tabPtMN(i) ; tabPtMN(i)=tmp
             trie = .false.
          endif
       enddo
       if(trie) exit
    enddo
  end subroutine sortTabPtMN

  subroutine getGigogneData(szSA,nbNode)
    integer,intent(in) :: szSA,nbNode

    type(t_segArc)        :: sa
    type(t_cellulePlaced) :: tcp
    type(t_celluleBase)   :: tcb
    integer               :: i,j,numNod,typCell,sectori,sectorj,ds,cluster,lg,szM
    logical               :: cl
    integer,dimension(:),allocatable :: neutronicMix,mrg

    allocate(tabNodeGigSect(nbNode),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: getGigogneData(1) => allocation pb")
    !initialisation
    do i = 1,nbNode
       tabNodeGigSect(i)%indTabCellPlac = 0
       tabNodeGigSect(i)%ring = 0
       tabNodeGigSect(i)%sect = 0
       tabNodeGigSect(i)%neutronicMix = 0
       tabNodeGigSect(i)%merge = 0
       tabNodeGigSect(i)%dragSector = 0
       tabNodeGigSect(i)%clust = .false.
    enddo
    !remplissage
    do i = 1,szSA
       sa = tabSegArc(i)
       numNod = sa%nodeg !travail sur le node gauche
       if(numNod>0) then
          call setIntIfPos(tabNodeGigSect(numNod)%indTabCellPlac,sa%indCellPg)
          call setIntIfPos(tabNodeGigSect(numNod)%ring,sa%mixg)
          tabNodeGigSect(numNod)%clust=tabNodeGigSect(numNod)%clust.or.sa%clusg
          if(sa%typ==tarc) then
             !sur le node interieur, en cas de sectorisation exterieure,
             !il doit apparaitre une discontinuite de sectorisation.
             !=> on passe outre le test de coherence de sect
             tabNodeGigSect(numNod)%sect = sa%sectg
          else
            call setIntIfPos(tabNodeGigSect(numNod)%sect,sa%sectg)
          endif
        endif
       numNod = sa%noded !travail sur le node droit
       if(numNod>0) then
          call setIntIfPos(tabNodeGigSect(numNod)%indTabCellPlac,sa%indCellPd)
          call setIntIfPos(tabNodeGigSect(numNod)%ring,sa%mixd)
          tabNodeGigSect(numNod)%clust=tabNodeGigSect(numNod)%clust.or.sa%clusd
          call setIntIfPos(tabNodeGigSect(numNod)%sect,sa%sectd)
       endif
    enddo
    !calcul des secteurs et des milieux neutroniques
    do i = 1,nbNode
       tcp = tabCellulePlaced(tabNodeGigSect(i)%indTabCellPlac)
       tcb = tabCelluleBase(tcp%indice)
       typCell = tcb%sv(1)
       cluster = tcb%sv(13)
       sectori = tcb%sv(14)
       sectorj = tcb%sv(15)
       !creation des tableaux de reference sur les milieux et les merges
       !par ajout des donnees eventuelles sur les clusters
       if(cluster/=0 .and. sectori/=S_not) call XABORT("G2S: CLUSTER not&
            &allowed with SECT")
       if(cluster==0) then
          allocate(neutronicMix(size(tcb%mix)),mrg(size(tcb%merge)),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: getGigogneData(2) => allocation pb")
          if((typCell == G_Hex).and.(tcb%name == '/')) then
            if(size(tcb%mix) /= nbNode) call XABORT("G2S: getGigogneData=> invalid size")
            neutronicMix(1) = tcb%mix(tabNodeGigSect(i)%indTabCellPlac)
          else
            neutronicMix(:size(tcb%mix)) = tcb%mix
          endif
          mrg(:size(tcb%merge)) = tcb%merge
       else
          lg = size(tcb%mix)
          do j = 1,cluster
             lg = lg + size(tcb%cluster(j)%mix)
          enddo
          allocate(neutronicMix(lg),mrg(lg),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: getGigogneData(3) => allocation pb")
          !remplissage de neutronicMix
          lg = size(tcb%mix)
          neutronicMix(1:lg) = tcb%mix(1:lg)
          do j = 1,cluster
             szM = size(tcb%cluster(j)%mix)
             neutronicMix(lg+1:lg+szM) = tcb%cluster(j)%mix(1:szM)
             lg = lg + szM
          enddo
          !remplissage de mrg
          ! remplissage par defaut pour le merge des clusters
          mrg(1:lg) =  (/(j,j=1,lg)/)
          ! remplissage du debut avec le merge de la cellule
          lg = size(tcb%merge)
          mrg(1:lg) = tcb%merge(1:lg)
       endif
       !traitement
       ds = 0
       cl = .false.
       select case(typCell)
       case(G_Car2d,G_Tri,G_Tube)
          if(sectori/=S_not) call XABORT("G2S: no SECT allowed for CAR2D, &
               &TRI2D, or TUBE")
          ds = tabNodeGigSect(i)%ring
          cl = cluster/=0 .and. tabNodeGigSect(i)%clust
       case(G_Carcel)
          select case (sectori)
          case(S_not)
             ds = tabNodeGigSect(i)%ring
             cl = cluster/=0 .and. tabNodeGigSect(i)%clust
          case(S_X_TOT, S_T_TOT)
             if (tabNodeGigSect(i)%ring<= sectorj) then
                ds = sectorj
             else
                ds = (tabNodeGigSect(i)%ring -sectorj-1)*4   &
                     + sectorj                                 &
                     + tabNodeGigSect(i)%sect
             endif
          case(S_TX_TOT, S_TXS_TOT, S_WM_TOT)
             if (tabNodeGigSect(i)%ring<= sectorj) then
                ds =tabNodeGigSect(i)%ring
             else
                ds = (tabNodeGigSect(i)%ring -sectorj-1)*8   &
                     + sectorj                                 &
                     + tabNodeGigSect(i)%sect
             endif
          end select

       case(G_Hex)
          if(sectori==S_not) then
             ds = 1
          else
             ds = tabNodeGigSect(i)%sect
          endif

       case(G_Hexcel)
          if(sectori==S_not) then
             ds = tabNodeGigSect(i)%ring
             cl = cluster/=0 .and. tabNodeGigSect(i)%clust
          else if(sectori==S_X_tot) then
             if (tabNodeGigSect(i)%ring<= sectorj) then
                ds = sectorj
             else
                ds = (tabNodeGigSect(i)%ring-sectorj-1)*6   &
                     + sectorj                                 &
                     + tabNodeGigSect(i)%sect
             endif
          else
             call XABORT("G2S: value for SECT not allowed")
          endif
       case default
          call XABORT("G2S: internal error in subroutine getGigogneData")
       end select
       tabNodeGigSect(i)%neutronicMix = neutronicMix(ds)
       ! correction of MERGE bug by Alain Hebert (January 2016)
       ! tabNodeGigSect(i)%merge = mrg(ds)
       ! new correction of MERGE bug for clusters by Alain Hebert (November 2019)
       if(.not.cl) then
         tabNodeGigSect(i)%merge = i
       else
         tabNodeGigSect(i)%merge = nbNode+ds
       endif
       tabNodeGigSect(i)%dragSector = ds
       deallocate(neutronicMix,mrg)
    enddo
! CS-IV : visualisation pour debug
!    call PrintTabNodeGigSect(nbNode)

    !remontage des infos
    do i = 1,szSA
       numNod = tabSegArc(i)%nodeg !travail sur le node gauche
       if(numNod>0) &
            tabSegArc(i)%neutronicMixg = tabNodeGigSect(numNod)%neutronicMix
       numNod = tabSegArc(i)%noded !travail sur le node droit
       if(numNod>0) &
            tabSegArc(i)%neutronicMixd = tabNodeGigSect(numNod)%neutronicMix
    enddo
  end subroutine getGigogneData

  subroutine setIntIfPos (toModif,valToSet)
    implicit none
    integer,intent(inout) :: toModif
    integer,intent(in)    :: valToSet
    
    if(valToSet<=0) return !rien a faire
    if(toModif<=0) then
       toModif = valToSet
! correction of plain CAR2D bug by Alain Hebert (May 2016)
    else if(toModif/=valToSet) then
       write(6,*) " setIntIfPos: ",toModif,"/=",valToSet
       call XABORT("G2S: internal error in subroutine setIntIfPos")
    endif
  end subroutine setIntIfPos

  subroutine PrintTabPtNode(size)
    integer, intent(in) :: size
    integer :: i,j
    write(*,*) "Impression de TabPtNode de ",size," elements"
    do i=1, size
       write(*,10) i
       write(*,20) TabPtNode(i)%position
       write(*,30) TabPtNode(i)%sz
       do j = 1, TabPtNode(i)%sz
          write(*,40) TabPtNode(i)%listOfSA(j),TabPtNode(i)%listOfDeparture(j),&
          TabPtNode(i)%listOfDirection(j)
       end do
    enddo
10  format(("**** element ****", i6))
20  format("*----------- position = ", f13.6,";",f13.6)
30  format("*----------- size     = ", i6)
40  format("*----------- SA/Depart/dir = ", i5,"/",l5,"/",F13.6)
  end subroutine PrintTabPtNode

  subroutine PrintTabNodeGigSect(size)
    integer, intent(in) :: size
    integer :: i

    do i=1, size
       write(*,10) i
       write(*,20) tabNodeGigSect(i)%indtabcellplac
       write(*,30) tabNodeGigSect(i)%ring
       write(*,35) tabNodeGigSect(i)%clust
       write(*,40) tabNodeGigSect(i)%sect
       write(*,50) tabNodeGigSect(i)%neutronicmix
       write(*,60) tabNodeGigSect(i)%merge
       write(*,70) tabNodeGigSect(i)%dragsector
    end do
10  format("**** element ****", i6)
20  format("*+++++++++++ IndTabCellPlac = ", i6)
30  format("*+++++++++++ Ring           = ", i6)
35  format("*+++++++++++ Cluster ring   = ", l6)
40  format("*+++++++++++ Sect           = ", i6)
50  format("*+++++++++++ NeutronicMix   = ", i6)
60  format("*+++++++++++ Merge          = ", i6)
70  format("*+++++++++++ DragSector     = ", i6)
  end subroutine PrintTabNodeGigSect
end module ptNodes
