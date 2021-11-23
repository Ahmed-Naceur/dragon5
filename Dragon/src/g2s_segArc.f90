!
!-----------------------------------------------------------------------
!
!Purpose:
! Collect processing functions related to geometric elements.
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
! Une structure "t_segArc" est definie, avec tous les champs necessaires aux
! calculs. On stocke par ailleurs ces elements dans un tableau global de
! pointeurs "tabSegArc".
! \\\\
! variable globale:
!  - tabSegArc : le tableau de travail principal
! \\\\
! fonctions du module:
!  - initializeTabSegArc : mise a zero des tableaux
!  - destroyTabSegArc : liberation de la memoire
!  - createSeg : constructeur de segment
!  - createArc : constructeur de cercle ou d'arc de cercle
!  - giveOrigine : retourne les coordonnees de l'origine d'un segArc
!  - giveExtremite : retourne les coordonnees de l'extremite d'un segArc
!  - extremitesArc : determination des points extremites d'un arc
!  - estColieaire : test de la colinearite de 2 segments
!  - isAngleInArc : teste l'appartenance d'un angle au domaine d'un arc
!  - isIn : teste l'appartenance d'un point, sur le support d'un segment ou
!           d'un arc, au segment ou a l'arc
!  - isSameWay : teste si deux segments alignes ont meme orientation
!  - turnBackSide : retourne le segment donne en entree
!  - adjustMix : ajuste les milieux a droite et a gauche d'un segment, par
!                rapport a ceux d'un segment de reference le contenant
!  - interSgAr : calcule les points d'intersection entre un segment et un arc
!  - interSgSg : calcule les points d'intersection entre deux segments
!  - interSgSg_V2 : recherche les points d'itersection entre deux segments 
!                   (extremites inclues)
!  - interArcArc : calcule les points d'intersection entre deux arcs
!  - translateAndTurnTriSg : effectue une translation-rotation d'un cote
!                            de triangle (fonction tres specifique)
!  - giveExtremalsAngles : donne les angles d'incidence des extrmites d'un
!                          element, par rapport a un point
!  - drawSegArc : affiche le tableau des elements dans un fichier postscript
!  - coupeCercle : intersecte les couronnes et les clusters
!  - splitSegsForSect : coupe l'encadrement d'une cellule rectangulaire en cas
!                       de sectorisation
!  - majSectori : met a jour la donnee sur le secteur des elements
!  - addSegsAndClean : ajout de nouveaux segments lors de la mise cote a cote
!                      des differentes cellules, et elimination des doublons
!  - cutClusters : prise en compte des intersections couronnes/cluster, et
!                  elimination de elements internes aux clusters
!  - segWithSameCoord :  true if the transmitted segments have the same coordinates
!
!-----------------------------------------------------------------------
!
module segArc
  use constUtiles
  use constType
  use derivedPSPLOT

  implicit none

  interface printSegArc
     module procedure printSegArc1, printSegArc2
  end interface printSegArc

  type t_segArc
     integer           :: typ   !1=segment, 2=arc de cercle
     double precision  :: x,y   !origine si 1 , centre si 2
     double precision  :: dx,dy !extremite si 1
     double precision  :: r     !rayon si 2
     double precision  :: a,b   !angle debut et fin si 2 (entre -pi et pi)
     integer           :: mixg  =0!numero de couronne a gauche ou interieur
     integer           :: mixd  =0!numero de couronne a droite ou exterieur
     integer           :: nodeg =0!node a gauche ou interieur
     integer           :: noded =0!node a droite ou exterieur
     integer           :: indCellPg=0 !cellulePlaced d'origine a gauche
     integer           :: indCellPd=0 !cellulePlaced d'origine a droite
     integer           :: sectg=0 !numero du secteur a gauche dans la cellule
     integer           :: sectd=0 !numero du secteur a droite dans la cellule
     integer           :: neutronicMixg=0 !milieux neutronique a gauche
     integer           :: neutronicMixd=0 !milieux neutronique a droite
     logical           :: clusg=.false. ! circle belong to a cluster 
     logical           :: clusd=.false. ! circle belong to a cluster 
  end type t_segArc

  integer,parameter :: tseg=1 , tcer=2 , tarc=3 , fooMix=-99 , fooNode=-999 , &
       & tRec=1 , tHex=2 , tTri=3

  !variable globale de type tableau d'elements geometriques
  type(t_segArc), dimension(:), allocatable :: tabSegArc

  ! storage of  cutting straights and segments
  type(t_segArc), dimension(:), allocatable :: tabStrCut,tabSegCut

  type segArcArrayBis
     type(t_segArc) :: sa
     logical        :: keep
  end type segArcArrayBis

  type segArcArrayTer
     type(t_segArc) :: sa
     logical        :: keep
     logical        :: cl   !condition limite ou cluster selon les cas
  end type segArcArrayTer

  ! Programmation defensive
  integer :: alloc_ok

contains

  function createSeg(ox,oy,ex,ey,mg,md)
    double precision,intent(in) :: ox,oy,ex,ey
    integer,intent(in)          :: mg,md
    type(t_segArc)              :: createSeg

    createSeg%typ   = tseg
    createSeg%x     = ox ; createSeg%y     = oy
    createSeg%dx    = ex ; createSeg%dy    = ey
    createSeg%r     = 0.d0
    createSeg%a     = 0.d0      ; createSeg%b     = 0.d0
    createSeg%mixg  = mg        ; createSeg%mixd  = md
    createSeg%nodeg = fooNode   ; createSeg%noded = fooNode
    createSeg%indCellPg = 0     ; createSeg%indCellPd = 0
    createSeg%sectg = 0         ; createSeg%sectd = 0
    createSeg%neutronicMixg = 0 ; createSeg%neutronicMixd = 0
    createSeg%clusg = .false.   ; createSeg%clusd = .false.
  end function createSeg

  function copySegExceptOrigin(ox,oy,Seg)
    double precision,intent(in) :: ox,oy
    type(t_segArc),intent(in)   :: Seg
    type(t_segArc)              :: copySegExceptOrigin
    if (Seg%typ /= tseg) call XABORT("copySegExceptOrigin : inconsistency")
    copySegExceptOrigin    = Seg
    copySegExceptOrigin%x  = ox ; copySegExceptOrigin%y  = oy
  end function copySegExceptOrigin

  function copySegExceptEnd(ex,ey,Seg)
    double precision,intent(in) :: ex,ey
    type(t_segArc),intent(in)   :: Seg
    type(t_segArc)              :: copySegExceptEnd
    if (Seg%typ /= tseg) call XABORT("copySegExceptEnd : inconsistency")
    copySegExceptEnd    = Seg
    copySegExceptEnd%dx = ex ; copySegExceptEnd%dy  = ey
  end function copySegExceptEnd

  function copySegWithNewExtrems(ox,oy,ex,ey,seg)
    double precision, intent(in) :: ox, oy, ex, ey
    type(t_segArc),intent(in)   :: Seg
    type(t_segArc)              :: copySegWithNewExtrems
    if (Seg%typ /= tseg) call XABORT("copySegWithNewExtrems: inconsistency")
    copySegWithNewExtrems    = Seg
    copySegWithNewExtrems%x  = ox ; copySegWithNewExtrems%y  = oy
    copySegWithNewExtrems%dx = ex ; copySegWithNewExtrems%dy = ey
  end function copySegWithNewExtrems

  function copyUTurnSeg(seg)
    type(t_segArc), intent(in) :: seg
    type(t_segArc)             :: copyUTurnSeg
    if (Seg%typ /= tseg) call XABORT("copyUTurnSeg : inconsistency")
    copyUTurnSeg%typ = tseg
    copyUTurnSeg%x  = seg%dx    ; copyUTurnSeg%y  = seg%dy
    copyUTurnSeg%dx = seg%x     ; copyUTurnSeg%dy = seg%y
    copyUTurnSeg%r = seg%r ; copyUTurnSeg%a = seg%a ; copyUTurnSeg%b = seg%b
    copyUTurnSeg%mixg      = seg%mixd      ; copyUTurnSeg%mixd      = seg%mixg
    copyUTurnSeg%nodeg     = seg%noded     ; copyUTurnSeg%noded     = seg%nodeg
    copyUTurnSeg%indCellPg = seg%indCellPd ; copyUTurnSeg%indCellPd = seg%indCellPg 
    copyUTurnSeg%sectg     = seg%sectd     ; copyUTurnSeg%sectd     = seg%sectg
    copyUTurnSeg%neutronicMixg = seg%neutronicMixd
    copyUTurnSeg%neutronicMixd = seg%neutronicMixg
  end function copyUTurnSeg

  function createArc(cx,cy,r,a,b,mi,me)
    double precision,intent(in) :: cx,cy,r,a,b
    integer,intent(in)          :: mi,me
    type(t_segArc)              :: createArc

    createArc%a     = a ; createArc%b     = b
    if (isEqual(createArc%a,createArc%b)) then ; createArc%typ  = tcer
    else ; createArc%typ  = tarc 
    end if
    createArc%x     = cx   ; createArc%y     = cy
    createArc%dx    = 0.d0 ; createArc%dy    = 0.d0
    createArc%r     = r
    createArc%mixg  = mi        ; createArc%mixd  = me
    createArc%nodeg = fooNode   ; createArc%noded = fooNode
    createArc%indCellPg = 0     ; createArc%indCellPd = 0
    createArc%sectg = 0         ; createArc%sectd = 0
    createArc%neutronicMixg = 0 ; createArc%neutronicMixd = 0
  end function createArc

  function copyArcExceptOrigin(o,Arc)
    double precision,intent(in) :: o
    type(t_segArc),intent(in)   :: Arc
    type(t_segArc)              :: copyArcExceptOrigin
    if (Arc%typ /= tarc) call XABORT("copyArcExceptOrigin : inconsistency")
    copyArcExceptOrigin    = Arc
    copyArcExceptOrigin%a  = o
  end function copyArcExceptOrigin

  function copyArcExceptEnd(e,Arc)
    double precision,intent(in) :: e
    type(t_segArc),intent(in)   :: Arc
    type(t_segArc)              :: copyArcExceptEnd
    if (Arc%typ /= tarc) call XABORT("copyArcExceptEnd : inconsistency")
    copyArcExceptEnd   = Arc
    copyArcExceptEnd%b = e
  end function copyArcExceptEnd

  function copyArcWithNewAngles(o,e,Arc)
    double precision,intent(in) :: o, e
    type(t_segArc),intent(in)   :: Arc
    type(t_segArc)              :: copyArcWithNewAngles
    if (Arc%typ /= tarc) call XABORT("copyArcWithNewAngles : inconsistency")
    copyArcWithNewAngles   = Arc
    copyArcWithNewAngles%a = o
    copyArcWithNewAngles%b = e
  end function copyArcWithNewAngles

  function createArcFromCircle(o,e,Circ)
    double precision,intent(in) :: o, e
    type(t_segArc),intent(in)   :: Circ
    type(t_segArc)              :: createArcFromCircle
    if (Circ%typ /= tcer) call XABORT("createArcFromCircle : inconsistency")
    createArcFromCircle     = Circ
    createArcFromCircle%typ = tarc
    createArcFromCircle%a   = o
    createArcFromCircle%b   = e
  end function createArcFromCircle

  subroutine giveOrigine(sa,xx,yy)
    type(t_segArc),intent(in)    :: sa
    double precision,intent(out) :: xx,yy
    select case(sa%typ)
    case(tseg) ; xx = sa%x ; yy = sa%y
    case(tarc) ; xx = sa%x+sa%r*cos(sa%a) ; yy=sa%y+sa%r*sin(sa%a)
    case(tcer) ; xx = sa%x+sa%r ; yy = sa%y
    end select
  end subroutine giveOrigine

  subroutine giveExtremite(sa,xx,yy)
    type(t_segArc),intent(in)    :: sa
    double precision,intent(out) :: xx,yy
    select case(sa%typ)
    case(tseg) ; xx = sa%dx ; yy = sa%dy
    case(tarc) ; xx = sa%x+sa%r*cos(sa%b) ; yy=sa%y+sa%r*sin(sa%b)
    case(tcer) ; xx = sa%x+sa%r ; yy = sa%y
    end select
  end subroutine giveExtremite

  subroutine giveFourPointsOnCircle(sa,P1,P2,P3,P4)
    type(t_segArc), intent(in) :: sa
    type(t_point), intent(out) :: P1,P2,P3,P4
    if (sa%typ/=tcer) call XABORT("giveFourPointsOnCircle : bad use")
    P1 = t_point(sa%x+sa%r,sa%y)
    P2 = t_point(sa%x-sa%r,sa%y)
    P3 = t_point(sa%x,sa%y+sa%r)
    P4 = t_point(sa%x,sa%y-sa%r)

  end subroutine giveFourPointsOnCircle

  subroutine extremitesArc(ar,pt1x,pt1y,pt2x,pt2y)
    type(t_segArc),intent(in)    :: ar
    double precision,intent(out) :: pt1x,pt1y,pt2x,pt2y

    pt1x=ar%x+ar%r*cos(ar%a) ; pt1y=ar%y+ar%r*sin(ar%a)
    pt2x=ar%x+ar%r*cos(ar%b) ; pt2y=ar%y+ar%r*sin(ar%b)
  end subroutine extremitesArc

  subroutine MediumPointOnArc(Ar,P)
    type(t_segArc), intent(in) :: Ar
    type(t_point), intent(out) :: P
    double precision :: angle
    if (Ar%a < Ar %b) then
       angle = angleNormal((Ar%a + Ar%b) * 0.5)
    else
       angle = angleNormal((Ar%a + Ar%b ) *0.5 + pi_c)
    endif
    p%x = Ar%x + Ar%r * cos(angle)
    p%y = Ar%y + Ar%r * sin(angle)
  end subroutine MediumPointOnArc

  function estColineaire(sa1,sa2)
    type(t_segArc),intent(in) :: sa1,sa2
    logical                   :: estColineaire

    if (sa1%typ/=tseg .or. sa2%typ/=tseg) then
       estColineaire = .false.
    else
       estcolineaire = estcoli(sa1%dx-sa1%x,sa1%dy-sa1%y,sa2%dx-sa2%x,sa2%dy-sa2%y)
    end if
  end function estColineaire

  function isAngleInArc(a,arc)
    double precision,intent(in) :: a
    type(t_segArc),intent(in)   :: arc
    integer                     :: isAngleInArc
    double precision :: angl
    !dit si un angle est sur l'intervalle d'un arc : 0 -> pas dessus,
    ! 1 -> c'est l'origine, 2 -> entre les deux, 3 -> c'est l'extremite
    if (arc%typ==tcer) then
       isAngleInArc = 2
       return
    end if

    angl = angleNormal(a)
    if      (isEqualConst(arc%a,angl)) then ; isAngleInArc = 1
    else if (isEqualConst(arc%b,angl)) then ; isAngleInArc = 3 
    else if (arc%a<arc%b)              then
       if((arc%a<angl).and.(angl<arc%b)) then ; isAngleInArc = 2
       else                                   ; isAngleInArc = 0
       end if
    else
       if((arc%b<angl).and.(angl<arc%a)) then ; isAngleInArc = 0
       else                                   ; isAngleInArc = 2
       end if
    end if
  end function isAngleInArc

  function isIn(ptx,pty,sa)
    double precision,intent(in) :: ptx,pty
    type(t_segArc),intent(in)   :: sa
    integer                     :: isIn

    double precision :: angl

    !dit si un point est sur un segment ou un arc : 0 -> pas dessus,
    ! 1 -> c'est l'origine, 2 -> entre les deux, 3 -> c'est l'extremite
    !ATTENTION :le point doit etre deja sur la droite ou le cercle d'appui
    if (sa%typ==tseg) then
       if (isEqualConst(sa%x,ptx) .and. isEqualConst(sa%y,pty)) then
          isIn = 1
       else if (isEqualConst(sa%dx,ptx) .and. isEqualConst(sa%dy,pty)) then
          isIn = 3
       else if (((sa%x-ptx)*(sa%dx-ptx)+(sa%y-pty)*(sa%dy-pty))<0) then
          isIn = 2
       else
          isIn = 0
       end if
    else 
       angl = calculeAngle(sa%x,sa%y,ptx,pty)
       isIn = isAngleInArc(angl,sa)
    end if
  end function isIn

  function isIn_V2(ptx,pty,sa)
    double precision,intent(in) :: ptx,pty
    type(t_segArc),intent(in)   :: sa
    integer                     :: isIn_V2

    double precision :: angl, tmp

    !dit si un point est sur un segment ou un arc : 0 -> pas dessus,
    ! 1 -> c'est l'origine, 2 -> entre les deux, 3 -> c'est l'extremite
    !ATTENTION :le point doit etre deja sur la droite ou le cercle d'appui
    if (sa%typ==tseg) then
       if (isEqualConst(sa%x,ptx) .and. isEqualConst(sa%y,pty)) then
          isIn_V2 = 1
       else if (isEqualConst(sa%dx,ptx) .and. isEqualConst(sa%dy,pty)) then
          isIn_V2 = 3
       else 
          tmp = (sa%x-ptx)*(sa%dx-ptx)+(sa%y-pty)*(sa%dy-pty)
          if (tmp < dp_0) then
             isIn_V2 = 2
          else if (isEqualConst(tmp,dp_0)) then
             isIn_V2 = 4
          else
             isIn_V2 = 0
          end if
       end if
    elseif (sa%typ == tarc) then
       angl = calculeAngle(sa%x,sa%y,ptx,pty)
       isIn_V2 = isAngleInArc(angl,sa)
    else
       tmp = sqrt((sa%x-ptx)*(sa%x-ptx)+(sa%y-pty)*(sa%y-pty))
       if (isEqualConst(tmp,sa%r)) then
          isIn_V2 = 2
       else
          isIn_V2 = 0
       endif
    end if
  end function isIn_V2

  function isSameWay(sg1,sg2)
    type(t_segArc),intent(in) :: sg1,sg2
    logical                   :: isSameWay
    !dit si 2 segments alignes ont meme orientation

    isSameWay=(((sg1%dx-sg1%x)*(sg2%dx-sg2%x)+(sg1%dy-sg1%y)*(sg2%dy-sg2%y))>0)
  end function isSameWay

  subroutine AddSA(SA,szSA)
    type(t_segArc), intent(in) :: SA
    integer, intent(inout) :: szSA
    if (szSA == size(tabSegArc)) call XABORT("AddSA : memory problem")
    szSA = szSA + 1
    tabSegArc(szSA) = SA
  end subroutine AddSA

  subroutine Add2SA(SA1,SA2,szSA)
    type(t_segArc), intent(in) :: SA1,SA2
    integer, intent(inout) :: szSA
    if (szSA == size(tabSegArc)) call XABORT("Add2SA : memory problem")
    szSA = szSA + 2
    tabSegArc(szSA-1) = SA1
    tabSegArc(szSA)   = SA2
  end subroutine Add2SA

  function ReplaceSA(newSA,iSA)
    type(t_segArc), intent(in)    :: newSA
    integer,        intent(in)    :: iSA
    type(t_segArc)                :: ReplaceSA
    ReplaceSA      = newSA
    tabSegArc(iSA) = newSA
  end function ReplaceSA

  subroutine SlideSA(iSA,szSA)
    integer, intent(in)    :: iSA
    integer, intent(inout) :: szSA
    tabSegArc(iSA:szSA-1) = tabSegArc(iSA+1:szSA)
    szSA = szSA - 1
  end subroutine SlideSA

  subroutine AddSC(SC,szSC)
    type(t_segArc), intent(in) :: SC
    integer, intent(inout)     :: szSC
    if (szSC == size(tabSegCut)) then
       print*,szSC,size(tabSegCut)
       call XABORT("AddSC : memory problem")
    endif
    szSC = szSC + 1
    tabSegCut(szSC) = SC
  end subroutine AddSC

  subroutine Add2SC(SC1,sC2,szSC)
    type(t_segArc), intent(in) :: SC1,SC2
    integer, intent(inout) :: szSC
    if (szSC == size(tabSegCut)) call XABORT("Add2SC : memory problem")
    szSC = szSC + 2
    tabSegCut(szSC-1) = SC1
    tabSegCut(szSC)   = SC2
  end subroutine Add2SC

  function ReplaceSC(newSC,iSC)
    type(t_segArc), intent(in)    :: newSC
    integer,        intent(in)    :: iSC
    type(t_segArc)                :: ReplaceSC
    ReplaceSC      = newSC
    tabSegCut(iSC) = newSC
  end function ReplaceSC
  
  subroutine OverlaidSegmentsManagement(n1, n2, iSC, szSC, SC, iSA, szSA, SA)
    integer,        intent(in)    :: n1, n2
    integer,        intent(inout) :: iSC, szSC, iSA, szSA
    type(t_SegArc), intent(inout) :: SC, SA
    type(t_segArc)                :: SC1, SC2, SC3, SA1, SA2
    select case(n1)
    case(11)             ! SC and SA have same origin
       select case(n2)
       case(23)          
          SA1 = CopySegExceptOrigin(SC%dx,SC%dy,SA)
          SA = ReplaceSA(SA1,iSA)
          iSA = iSA + 1
       case(32)          ! end of SA internal to SC
          SC1 = CopySegExceptOrigin(SA%dx,SA%dy,SC)
          SC2 = CopySegExceptEnd(SA%dx,SA%dy,SC)
          call AddSC(SC1,szSC)
          SC = ReplaceSC(SC2,iSC)
          call SlideSA(iSA,szSA)
       case(33)          ! SC and SA have same end
          call SlideSA(iSA,szSA)
       case default
          call XABORT("OverlaidSegmentsManagement : impossible circumstance")
       end select
    case(12)             ! the origin of SC is internal to SA
       select case (n2)
       case(21)
          SC1 = CopySegExceptEnd(SA%x,SA%y,SC)
          SC2 = CopySegExceptOrigin(SA%x,SA%y,SC)
          SA1 = CopySegExceptOrigin(SC%x,SC%y,SA)
          SC = ReplaceSC(SC1,iSC)
          SA = ReplaceSA(SA1,iSA)
          call AddSC(SC2,szSC)
          iSA = iSA + 1
       case(23)
          SC1 = CopySegExceptEnd(SA%x,SA%y,SC)
          SC2 = CopySegExceptOrigin(SA%x,SA%y,SC)
          SC3 = CopySegWithNewExtrems(SA%x,SA%y,SA%dx,SA%dy,SC)
          SA1 = CopySegExceptOrigin(SC%dx,SC%dy,SA)
          SC = ReplaceSC(SC3,iSC)
          SA = ReplaceSA(SA1,iSA)
          call AddSC(SC2,szSC)
          iSa = iSA + 1
       case(31)
          SC1 = CopySegExceptOrigin(SA%x,SA%y,SC)
          SC2 = CopySegWithNewExtrems(SA%dx,SA%dy,SA%x,SA%y,SC)
          SC = ReplaceSC(SC2,iSC)
          call AddSC(SC1,szSC)
          call SlideSA(iSA,szSA)
       case(32)
          if (isSameWay(SC,SA)) then
             SC1 = CopySegExceptEnd(SA%x,SA%y,SC)
             SC2 = CopySegExceptOrigin(SA%dx,SA%dy,SC)
             SC3 = CopySegWithNewExtrems(SA%x,SA%y,SA%dx,SA%dy,SC)
             SC =  ReplaceSC(SC3,iSC)
           else
             SC1 = CopySegExceptEnd(SA%dx,SA%dy,SC)
             SC2 = CopySegExceptOrigin(SA%x,SA%y,SC)
             SC3 = CopySegWithNewExtrems(SA%dx,SA%dy,SA%x,SA%y,SC)
             SC = ReplaceSC(SC3,iSC)
          endif
          call Add2SC(SC1,SC2,szSC)
          call SlideSA(iSA,szSA)
       case(33)
          SC1 = CopySegExceptEnd(SA%x,SA%y,SC)
          SC2 = CopySegWithNewExtrems(SA%x,SA%y,SA%dx,SA%dy,SC)
          SC = ReplaceSC(SC2,iSC)
          call AddSC(SC1,szSC)
          call SlideSA(iSA,szSA)
       case default
          call XABORT("OverlaidSegmentsManagement : impossible circumstance")
       end select
    case(13)          ! the origin of SC is the end of SA
       select case(n2)
       case(21)
          SA1 = CopySegExceptOrigin(SC%x,SC%y,SA)
          call AddSA(SA1,szSA)
          iSA = iSA + 1
       case(31)
          SC1 = CopySegWithNewExtrems(SA%dx,SA%dy,SA%x,SA%y,SC)
          SC = ReplaceSC(SC1,iSC)
          call SlideSA(iSA,szSA)
       case(32)
          SC1 = CopySegExceptEnd(SA%dx,SA%dy,SC)
          SC2 = CopySegExceptOrigin(SA%dx,SA%dx,SC)
          SC = ReplaceSC(SC2,iSC)
          call AddSC(SC1,szSC)
          call SlideSA(iSa,szSA)
       case default
          print*,"Pb n1,n2",n1,N2
          call XABORT("OverlaidSegmentsManagement : impossible circumstance")
       end select    
    case(21)
       select case(n2)
       case(12)
          SC1 = CopySegExceptEnd(SA%x,SA%y,SC)
          SC2 = CopySegExceptOrigin(SA%x,SA%y,SC)
          SA1 = CopySegExceptOrigin(SC%x,SC%y,SA)
          SC = ReplaceSC(SC1,iSC)
          call AddSC(SC2,szSC)
          SA = ReplaceSA(SA1,iSA)
          iSA = iSA + 1
       case(13)
          SA1 = CopySegExceptOrigin(SC%x,SC%y,SC)
          SA = ReplaceSA(SA1,iSA)
          iSA = iSA + 1
       case(23)
          if (isSameWay(SC,SA)) then
             SA1 = CopySegExceptEnd(SC%x,SC%y,SA)
             SA2 = CopySegExceptOrigin(SC%dx,SC%dy,SA)
          else
             SA1 = CopySegExceptEnd(SC%Dx,SC%Dy,SA)
             SA2 = CopySegExceptOrigin(SC%x,SC%y,SA)  
         endif
         SA = ReplaceSA(SA1,iSA)
         call AddSC(SC2,szSC)
         iSA = iSA + 1
       case(32)
          SC1 = CopySegExceptEnd(SA%dx,SA%dy,SC)
          SC2 = CopySegExceptOrigin(SA%dx,SA%dy,SC)
          SA1 = CopySegExceptEnd(SC%x,SC%y,SA)
          SC = ReplaceSC(SC1,iSC)
          call AddSC(SC2,szSC)
          SA = ReplaceSA(SA1,iSA)
          iSA = iSA + 1
       case(33)
          SA1 = CopySegExceptEnd(SC%x,SC%y,SC)
          SA = ReplaceSA(SA1,iSA)
          iSA = iSA + 1
       case default
          print*,"Pb n1,n2",n1,N2
          call XABORT("OverlaidSegmentsManagement : impossible circumstance")
       end select
    case(23)
       select case (n2)
       case(11)
          SA1 = CopySegExceptOrigin(SC%dx,SC%dy,SA)
          SA = ReplaceSA(SA1,iSA)
          iSA = iSA + 1
       case(12)  
          SC1 = CopySegExceptEnd(SA%x,SA%y,SC)
          SC2 = CopySegExceptOrigin(SA%x,SA%y,SC)
          SC = ReplaceSC(SC1,iSC)
          call AddSC(SC2,szSC)
          SA1 = CopySegExceptOrigin(SC%dx,SC%dy,SA)
          SA = ReplaceSA(SA1,iSA)
          iSA = iSA + 1
       case(21)
          if (isSameWay(SC,SA)) then
             SA1 = CopySegExceptEnd(SC%x,SC%y,SA)
             SA2 = CopySegExceptOrigin(SC%dx,SC%dy,SA)
          else
             SA1 = CopySegExceptEnd(SC%dx,SC%dy,SA)
             SA2 = CopySegExceptOrigin(SC%x,SC%y,SA)   
          endif
          SA = ReplaceSA(SA1,iSA)
          call AddSA(SA2,szSA)
          iSA = iSA + 1
       case(31)
          SA1 = CopySegExceptEnd(SC%dx,SC%dy,SA)
          SA = ReplaceSA(SA1,iSA)
          iSA = iSA + 1
       case(32)
          SC1 = CopySegExceptEnd(SA%dx,SA%dy,SC)
          SC2 = CopySegExceptOrigin(SA%dx,SA%dy,SC)
          SA1 = CopySegExceptEnd(SC%dx,SC%dy,SA)
          SC = ReplaceSC(SC1,iSC)
          SA = ReplaceSA(SA1,iSA)
          call AddSC(SC2,szSC)
          iSA = iSA + 1
       case default
          print*,"Pb n1,n2",n1,N2
          call XABORT("OverlaidSegmentsManagement : impossible circumstance")
       end select          
    case(31)
       select case (n2)
       case(12)
          SC1 = CopySegExceptEnd(SA%x,SA%y,SC)
          SC2 = CopySegExceptOrigin(SA%x,SA%y,SC)
          SC = ReplaceSC(SC1,iSC)
          call AddSc(SC2,szSC)
          call SlideSA(iSA,szSA)
       case(13) 
          SC1 = CopySegWithNewExtrems(SA%dx,SA%dy,SA%x,SA%y,SC)
          SC = ReplaceSC(SC1,iSC)
          call SlideSA(iSA,szSA)
       case(23)
          SA1 = CopySegExceptEnd(SC%dx,SC%dy,SA)
          SA = ReplaceSA(SA1,iSA)
          iSA = iSA + 1
       case default
          print*,"Pb n1,n2",n1,N2
          call XABORT("OverlaidSegmentsManagement : impossible circumstance")
       end select          
    case(32)
       select case (n2)
       case(11)
          SC1 = CopySegExceptEnd(SA%dx,SA%dy,SC)
          SC2 = CopySegExceptOrigin(SA%dx,SA%dy,SC)
          SC = ReplaceSC(SC1,iSC)
          call AddSC(SC2,szSC)
          iSA = iSA + 1
       case(12)
          if (isSameWay(SA,SC)) then
             SC1 = CopySegExceptEnd(SA%x,SA%y,SC)
             SC2 = CopySegExceptOrigin(SA%dx,SA%dy,SC)
             SC3 = CopySegWithNewExtrems(SA%x,SA%y,SA%dx,SA%dy,SC)
             SC = ReplaceSC(SC3,iSC)
          else
             SC1 = CopySegExceptEnd(SA%dx,SA%dy,SC)
             SC2 = CopySegExceptOrigin(SA%x,SA%y,SC)
             SC3 = CopySegWithNewExtrems(SA%dx,SA%dy,SA%x,SA%y,SC)
             SC = ReplaceSC(SC3,iSC)
          end if
          call Add2SC(SC1,SC2,szSC)
          call SlideSA(iSA,szSA)
       case(13)
          SC1 = CopySegExceptEnd(SA%dx,SA%dy,SC)
          SC2 = CopySegExceptOrigin(SA%dx,SA%dy,SC)
          SC = ReplaceSC(SC1,iSC)
          call AddSC(SC2,szSC)
          call SlideSA(iSA,szSA)
       case(21)
          SC1 = CopySegExceptEnd(SA%dx,SA%dy,SC)
          SC2 = CopySegExceptOrigin(SA%dx,SA%dy,SC)
          SA1 = CopySegExceptOrigin(SC%x,SC%y,SA)
          SC = ReplaceSC(SC1,iSC)
          SA = ReplaceSA(SA1,iSA)
          call AddSC(SC2,szSC)
          iSA = iSA + 1
       case(23)
          SC1 = CopySegExceptEnd(SA%dx,SA%dy,SC)
          SC2 = CopySegExceptOrigin(SA%dx,SA%dy,SC)
          SA1 = CopySegExceptOrigin(SC%dx,SC%dy,SA)
          SC = ReplaceSC(SC1,iSC)
          SA = ReplaceSA(SA1,iSA)
          call AddSC(SC2,szSC)
          iSA = iSA + 1
       case default
          print*,"Pb n1,n2",n1,N2
          call XABORT("OverlaidSegmentsManagement : impossible circumstance")
       end select     
    case(33)
       select case (n2)
       case(11)
          SC1 = CopySegWithNewExtrems(SA%x,SA%y,SA%dx,SA%dy,SC)
          SC = ReplaceSC(SC1,iSC)
          call SlideSA(iSA,szSA)
       case(12)
          SC1 = CopySegExceptEnd(SA%x,SA%y,SC)
          SC2 = CopySegWithNewExtrems(SA%x,SA%y,SA%dx,SA%dy,SC)
          SC = ReplaceSC(SC2,iSC)
          call AddSC(SC1,szSC)
          call SlideSA(iSA,szSA)
       case(21)
          SA1 = CopySegExceptEnd(SC%x,SC%y,SA)
          SA = ReplaceSA(SA1,iSA)
          iSA = iSA + 1
       case default
          print*,"Pb n1,n2",n1,N2
          call XABORT("OverlaidSegmentsManagement : impossible circumstance")
       end select
    case default
       print*,"Pb n1,n2",n1,N2
       call XABORT("OverlaidSegmentsManaement : ??")
    end select
  end subroutine OverlaidSegmentsManagement

  function turnBackSide(sg)
    type(t_segArc),intent(inout) :: sg
    type(t_segArc)               :: turnBackSide
    !tourne le sens du segment en entree

    turnBackSide=createSeg(sg%dx,sg%dy,sg%x,sg%y,sg%mixd,sg%mixg)
    turnBackSide%nodeg=sg%noded
    turnBackSide%noded=sg%nodeg
    turnBackSide%indCellPg=sg%indCellPd
    turnBackSide%indCellPd=sg%indCellPg
    turnBackSide%sectg=sg%sectd
    turnBackSide%sectd=sg%sectg
    turnBackSide%neutronicMixg=sg%neutronicMixd
    turnBackSide%neutronicMixd=sg%neutronicMixg
    sg=turnBackSide
  end function turnBackSide

  subroutine adjustMix(sgToAdjust,sgOfRef)
    type(t_segArc),intent(inout) :: sgToAdjust
    type(t_segArc),intent(in)    :: sgOfRef
    !ajuste les milieux a droite et a gauche de sgToAdjust
    !par rapport a ceux de sgOfRef (sgToAdjust doit etre inclus dans sgOfRef)
    if (sgOfRef%mixg/=fooMix) then
       sgToAdjust%mixg          = sgOfRef%mixg
       sgToAdjust%indCellPg     = sgOfRef%indCellPg
       sgToAdjust%sectg         = sgOfRef%sectg
       sgToAdjust%neutronicMixg = sgOfRef%neutronicMixg
    end if
    if (sgOfRef%mixd/=fooMix) then
       sgToAdjust%mixd          = sgOfRef%mixd
       sgToAdjust%indCellPd     = sgOfRef%indCellPd
       sgToAdjust%sectd         = sgOfRef%sectd
       sgToAdjust%neutronicMixd = sgOfRef%neutronicMixd
    end if
  end subroutine adjustMix

  logical function interSgAr_V2(sg,ar, n1, n2, P1, P2)
    ! Return true only if there is at least an intersection between sg and ar.
    ! If there is an point of intersection or an tangency point, P1 is this point.
    ! If there are two intersection points, P1 and P2 are these points.
    ! N1 and n2 let us pinpoint the point of intersection on both segment and arc.
    ! If there is no point of intersection, the values for n1 and n2 are 0.
    type(t_segArc),intent(in)    :: sg,ar
    integer, intent(out)::  n1, n2
    type (t_point), intent(out) :: P1, P2

    double precision    :: sox, soy, sfx, sfy, cx, cy, ray, aox, aoy
    double precision    :: prjx, prjy, d, u, v, norm, tmp, Ix, Iy, Jx, Jy
    integer :: IinS,IinA,JinS,JinA

    ! Initialization
    n1 = 0                  ; n2 = 0
    P1 = t_point(dp_0,dp_0) ; P2 = t_point(dp_0,dp_0)
    InterSgAr_V2 = .false.
    ! Coordinates
    sox = sg%x ; sfx = sg%dx ; soy = sg%y ; sfy = sg%dy
    cx = ar%x ; cy = ar%y  ; ray = ar%r

    ! Special case : ar is a point
    if ((ar%typ == tarc.and.isEqualConst(angleNormal(ar%a),angleNormal(ar%b))))   then
       call giveOrigine(ar,aox,aoy)
       if (isIn(aox,aoy,sg) /= 0) then
          interSgAr_V2 =.true. ; n1 = isIn(aox,aoy,sg)+10 
          P1 = t_point(aox,aoy)
       endif
       goto 10
    endif
    ! Special case : sg is a point
    if (isEqualConst(sox,sfx).and.isEqualConst(soy,sfy)) then
       d = sqrt((sox-cx)*(sox-cx) + (soy-cy)*(soy-cy))
       if (isEqualConst(d,ray)) then
          interSgAr_V2 =.true. ; n1 = 1+isIn(sox,soy,ar)*10 
          P1 = t_point(sox,soy)
       endif
       goto 10
    endif
    ! Distance Segment-Arc
    d = distance(cx,cy,sox,soy,sfx,sfy,prjx,prjy)
    ! Distance analysis
    if (.not. isEqualConst(d,ar%r)) then
       if (d < ray) then
          ! no one, one or two distinct intersection points
          ! P1 before and P2 after on the straight for segment
          u = sfx - sox ; v = sfy - soy ; norm = u*u + v*v
          if (.not.isEqualConst(sqrt(norm),dp_0)) then
             tmp = sqrt((ray*ray - d*d)/norm)
             u = tmp*u ; v = tmp*v
             Ix=prjx-u ; Iy=prjy-v ; IinS=isIn_V2(Ix,Iy,sg) ; IinA=isIn(Ix,Iy,ar)
             Jx=prjx+u ; Jy=prjy+v ; JinS=isIn_V2(Jx,Jy,sg) ; JinA=isIn(Jx,Jy,ar)
             if (IinS/=0.and.IinA/=0) then
                interSgAr_V2 = .true. ; n1 = IinS + 10*IinA ; P1 = t_point(Ix,Iy)
                if (JinS/=0.and.JinA/=0) then
                   n2 = JinS + 10*JinA ; P2 = t_point(Jx,Jy)
                   goto 10
                endif
             else if (JinS/=0.and.JinA/=0) then
                interSgAr_V2 = .true. ; n1 = JinS + 10*JinA ; P1 = t_point(Jx,Jy)
                goto 10
             end if
          else
             call XABORT("interSgAr_V2 : seg was not detected as a point but its norm is zero")
          end if
       end if
    else
       ! Tangency point
       Ix=prjx ; Iy=prjy ; IinS=isIn(Ix,Iy,sg) ; IinA=isIn(Ix,Iy,ar)
       if ((IinS/=0.and.IinA/=0)) then
          interSgAr_V2=.true. ; n1 = 100 + 10 * IinA + IinS ; P1 = t_point(Ix,Iy)
       end if
    end if
    ! 
10  continue
    return
  end function interSgAr_V2


  function interSgAr(sg,ar,pt1x,pt1y,pt2x,pt2y)
    type(t_segArc),intent(in)    :: sg,ar
    double precision,intent(out) :: pt1x,pt1y,pt2x,pt2y
    integer                      :: interSgAr

    double precision    :: prjx,prjy,d
    double precision    :: u,v,tmp,Ix,Iy,Jx,Jy
    integer :: IinS,IinA,JinS,JinA
    !renvoie le nombre de points d'intersection entre sg et ar,
    !et leurs coordonnees eventuelles
    !Si un seul point d'intersection, renvoi 1 si pt==I et -1 si pt==J
    !Si l'arc et le segment sont tangents, renvoie 2, et pt1=pt2

    interSgAr = 0
    pt1x=0.d0 ; pt1y=0.d0 ; pt2x=1.d0 ; pt2y=1.d0
    d = distance(ar%x,ar%y,sg%x,sg%y,sg%dx,sg%dy,prjx,prjy)
    if (.not. isEqualConst(d,ar%r)) then
       if (d>ar%r) then
          interSgAr = 0
       else
          !2 points d'intersection distincts I et J avec
          !I avant J sur la droite d'appui du segment
          u=sg%dx-sg%x ; v=sg%dy-sg%y
          tmp=sqrt((ar%r*ar%r-d*d)/(u*u+v*v))
          u=tmp*u ; v=tmp*v
          Ix=prjx-u ; Iy=prjy-v ; IinS=isIn(Ix,Iy,sg) ; IinA=isIn(Ix,Iy,ar)
          Jx=prjx+u ; Jy=prjy+v ; JinS=isIn(Jx,Jy,sg) ; JinA=isIn(Jx,Jy,ar)
          if ((IinS/=0.and.IinA/=0).and.(IinS==2.or.IinA==2)) then
             interSgAr=1 ; pt1x=Ix ; pt1y=Iy
          end if
          if ((JinS/=0.and.JinA/=0).and.(JinS==2.or.JinA==2)) then
             if (interSgAr==1) then
                interSgAr=2 ; pt2x=Jx ; pt2y=Jy
             else
                interSgAr=-1 ; pt1x=Jx ; pt1y=Jy
             end if
          end if
       end if
    else
       Ix=prjx ; Iy=prjy ; IinS=isIn(Ix,Iy,sg) ; IinA=isIn(Ix,Iy,ar)
       if ((IinS/=0.and.IinA/=0).and.(IinS==2.or.IinA==2)) then
          interSgAr=2 ; pt1x=Ix ; pt1y=Iy ; pt2x=Ix ; pt2y=Iy
       end if
    end if
  end function interSgAr

  logical function InterSgSg_V2(sg1,sg2,n1, n2, P1,P2)
    type(t_segArc),intent(in)    :: sg1,sg2
    integer, intent(inout) :: n1, n2
    type(t_point), intent(out) :: P1, P2

    double precision :: xo1, xf1, xo2, xf2, yo1, yf1, yo2, yf2
    double precision :: min1, min2, max1, max2, m, prodvect
    double precision :: a1, a2, b1, b2, y

    ! Initialization
    n1 = 0                  ; n2 = 0
    P1 = t_point(dp_0,dp_0) ; P2 = t_point(dp_0,dp_0)
    InterSgSg_V2 = .false.

    ! Coordinates
    xo1=sg1%x ; xf1 =sg1%dx ; yo1 = sg1%y ; yf1 = sg1%dy
    xo2=sg2%x ; xf2 =sg2%dx ; yo2 = sg2%y ; yf2 = sg2%dy

    ! Very Special case : point instead of segment
    min1 = sqrt((xf1-xo1)*(xf1-xo1)+(yf1-yo1)*(yf1-yo1))
    min2 = sqrt((xf2-xo2)*(xf2-xo2)+(yf2-yo2)*(yf2-yo2))
    if (isEqualConst(min1,dp_0)) then
       prodvect = (xf2-xo2)*(yo1-yo2)-(yf2-yo2)*(xo1-xo2)
       if (isEqualConst(prodvect,dp_0) .and. (isIn(xo1,yo1,Sg2)/=0)) then
          InterSgSg_V2 = .true. ; n1 = 1 ; P1 = t_point(xo1,yo1)
       endif
       goto 10
    elseif(isEqualConst(min2,dp_0)) then
       prodvect = (xf1-xo1)*(yo2-yo1)-(yf1-yo1)*(xo2-xo1)
       if (isEqualConst(prodvect,dp_0) .and. (isIn(xo2,yo2,Sg1)/=0)) then
          InterSgSg_V2 = .true. ; n1 = 1 ; P1 = t_point(xo2,yo2)
       endif
       goto 10
    else

       ! Special case : two vertical segments
       if (IsEqualConst(xo1,xf1).and.IsEqualConst(xo2,xf2)) then
          if (IsEqualConst(xo1,xo2)) then
             min1 = min(yo1,yf1) ; max1 = max(yo1,yf1)
             min2 = min(yo2,yf2) ; max2 = max(yo2,yf2)
             if (isEqualConst(min1,max2).or.isEqualConst(min2,max1)) then
                InterSgSg_V2 = .true. ;  n1 = 1 
                if (isEqualConst(min1,max2)) P1 = t_point(xo1,min1)
                if (isEqualConst(max1,min2)) P1 = t_point(xo2,min2) 
             elseif (max2 < min1 .or. max1 < min2) then
                InterSgSg_V2 = .false.
             else
                InterSgSg_V2 = .true. 
                n1 = 1 ; P1 = t_point(xo1,max(min1,min2))
                n2 = 1 ; P2 = t_point(xo1,min(max1,max2))
             endif
          endif
       else
          ! Special case : only one vertical segment
          if (isEqualConst(xo1,xf1).and..not.isEqualConst(xo2,xf2)) then 
             a2 = (yf2 - yo2) / (xf2 - xo2) ; b2 = yo2 - a2 * xo2
             if ((isEqualConst(xo1,min(xo2,xf2)).or.xo1>min(xo2,xf2)) .and.&
                  (isEqualConst(xo2,max(xo2,xf2)).or.xo1<max(xo2,xf2)))then
                y = a2 * xo1 + b2
                if ((isEqualConst(y,min(yo1,yf1)).or.y>min(yo1,yf1)).and.&
                     (isEqualConst(y,max(yo1,yf1)).or.y<max(yo1,yf1)).and.&
                     (isEqualConst(y,min(yo2,yf2)).or.y>min(yo2,yf2)).and.&
                     (isEqualConst(y,max(yo2,yf2)).or.y<max(yo2,yf2))) then
                   InterSgSg_V2 = .true. ; n1 = 1 ; P1 = t_point(xo1,y)
                else
                   InterSgSg_V2 = .false.
                endif
             else
                InterSgSg_V2 = .false.
             endif
          elseif (isEqualConst(xo2,xf2).and..not.isEqualConst(xo1,xf1)) then 
             if ((isEqualConst(xo2,min(xo1,xf1)).or.xo2>min(xo1,xf1) ).and.&
                  (isEqualConst(xo2,max(xo1,xf1)).or.xo2<max(xo1,xf1))) then
                a1 = (yf1 - yo1) / (xf1 - xo1) ; b1 = yo1 - a1 * xo1
                y = a1 * xo2 + b1
                if ((isEqualConst(y,min(yo1,yf1)).or.y>min(yo1,yf1)).and.&
                     (isEqualConst(y,max(yo1,yf1)).or.y<max(yo1,yf1)).and.&
                     (isEqualConst(y,min(yo2,yf2)).or.y>min(yo2,yf2)).and.&
                     (isEqualConst(y,max(yo2,yf2)).or.y<max(yo2,yf2))) then
                   InterSgSg_V2 = .true. ; n1 = 1 ; P1 = t_point(xo2,y)
                else
                   InterSgSg_V2 = .false.
                endif
             else
                InterSgSg_V2 = .false.
             endif
          else
             ! general case
             a1 = (yf1 - yo1) / (xf1 - xo1) ; b1 = yo1 - a1 * xo1
             a2 = (yf2 - yo2) / (xf2 - xo2) ; b2 = yo2 - a2 * xo2
             min1 = min(xo1,xf1) ; max1 = max(xo1,xf1)
             min2 = min(xo2,xf2) ; max2 = max(xo2,xf2)
             ! parallel lines
             if (isEqualCOnst(a1,a2)) then
                if (.not. isEqualConst(b1,b2)) then
                   InterSgSg_V2 = .false.
                else
                   if (max2 < min1 .or. max1 < min2) then
                      InterSgSg_V2 = .false.
                   elseif(isEqualConst(min1,max2).or. &
                          isEqualCOnst(max1,min2)) then
                      InterSgSg_V2 = .true. ; n1 = 1
                      if (isEqualConst(min1,max2)) then
                         P1 = t_point(min1,a1*min1+b1)
                      else
                         P1 = t_point(min2,a2*min2+b2)
                      endif
                   else
                      InterSgSg_V2 = .true.
                      m = max(min1,min2)
                      n1 = 1 ; P1 = t_point(m, a1*m+b1)
                      m = min(max1,max2)
                      n2 = 1 ; P2 = t_point(m, a1*m+b1)
                   endif
                endif
                ! Secant lines
             else
                m = (b2-b1)/(a1-a2)
                if ( (isEqualConst(m,min1).or.m>min1) .and. &
                     (isEqualConst(m,max1).or.m<max1) .and. &
                     (isEqualConst(m,min2).or.m>min2) .and. &
                     (isEqualConst(m,max2).or.m<max2)) then
                   InterSgSg_V2 = .true.
                   n1 = 1; P1 = t_point(m, a1*m+b1)
                else
                   InterSgSg_V2 = .false.
                endif
             endif
          end if
       end if
    endif
    ! Points qualification
10  continue
    if (n1 /= 0) then
       xo1 = P1%x ; yo1 = P1%y
       n1 = isIn(xo1,yo1,Sg1) + 10 * isIn(xo1,yo1,Sg2)
    endif
    if (n2 /= 0) then       
       xo1 = P2%x ; yo1 = P2%y
       n2 = isIn(xo1,yo1,Sg1) + 10 * isIn(xo1,yo1,Sg2)
       ! bug correction made on May 17, 2019
       if((n1.EQ.13).and.(n2.EQ.32)) n2=21
    endif
  end function InterSgSg_V2

  function interSgSg(sg1,sg2,ix,iy)
    type(t_segArc),intent(in)    :: sg1,sg2
    double precision,intent(out) :: ix,iy
    logical                      :: interSgSg

    double precision :: A1,B1,C1,A2,B2,C2,delta
    ix = 0.d0 ; iy = 0.d0
    A1=sg1%y-sg1%dy ; B1=sg1%dx-sg1%x ; C1=sg1%x*sg1%dy-sg1%dx*sg1%y
    A2=sg2%y-sg2%dy ; B2=sg2%dx-sg2%x ; C2=sg2%x*sg2%dy-sg2%dx*sg2%y
    delta = A1*B2-A2*B1
    interSgSg = .not. isEqualConst(delta,0.d0)
    if (.not. interSgSg) return
    delta = 1./delta
    ix = (B1*C2-B2*C1)*delta ; iy = (A2*C1-A1*C2)*delta
    interSgSg = (isIn(ix,iy,sg1)==2 .and. isIn(ix,iy,sg2)==2)
  end function interSgSg

  function translateAndTurnTriSg(cx,cy,turn,sg)
    double precision,intent(in) :: cx,cy
    integer,intent(in)          :: turn
    type(t_segArc),intent(in)   :: sg
    type(t_segArc)              :: translateAndTurnTriSg

    translateAndTurnTriSg = sg
    select case(turn)
    case(1)
       translateAndTurnTriSg%x  =  sg%x  + cx
       translateAndTurnTriSg%y  =  sg%y  + cy
       translateAndTurnTriSg%dx =  sg%dx + cx
       translateAndTurnTriSg%dy =  sg%dy + cy
    case(2)
       translateAndTurnTriSg%x  =  5.d-1*sg%x  - sqrt3_2d*sg%y  + cx
       translateAndTurnTriSg%y  =  5.d-1*sg%y  + sqrt3_2d*sg%x  + cy
       translateAndTurnTriSg%dx =  5.d-1*sg%dx - sqrt3_2d*sg%dy + cx
       translateAndTurnTriSg%dy =  5.d-1*sg%dy + sqrt3_2d*sg%dx + cy       
    case(3)
       translateAndTurnTriSg%x  = -5.d-1*sg%x  - sqrt3_2d*sg%y  + cx
       translateAndTurnTriSg%y  = -5.d-1*sg%y  + sqrt3_2d*sg%x  + cy
       translateAndTurnTriSg%dx = -5.d-1*sg%dx - sqrt3_2d*sg%dy + cx
       translateAndTurnTriSg%dy = -5.d-1*sg%dy + sqrt3_2d*sg%dx + cy       
    case(4)
       translateAndTurnTriSg%x  = -sg%x  + cx
       translateAndTurnTriSg%y  = -sg%y  + cy
       translateAndTurnTriSg%dx = -sg%dx + cx
       translateAndTurnTriSg%dy = -sg%dy + cy
    case(5)
       translateAndTurnTriSg%x  = -5.d-1*sg%x  + sqrt3_2d*sg%y  + cx
       translateAndTurnTriSg%y  = -5.d-1*sg%y  - sqrt3_2d*sg%x  + cy
       translateAndTurnTriSg%dx = -5.d-1*sg%dx + sqrt3_2d*sg%dy + cx
       translateAndTurnTriSg%dy = -5.d-1*sg%dy - sqrt3_2d*sg%dx + cy       
    case(6)
       translateAndTurnTriSg%x  =  5.d-1*sg%x  + sqrt3_2d*sg%y  + cx
       translateAndTurnTriSg%y  =  5.d-1*sg%y  - sqrt3_2d*sg%x  + cy
       translateAndTurnTriSg%dx =  5.d-1*sg%dx + sqrt3_2d*sg%dy + cx
       translateAndTurnTriSg%dy =  5.d-1*sg%dy - sqrt3_2d*sg%dx + cy       
    case(7)
       translateAndTurnTriSg%mixg  =  sg%mixd
       translateAndTurnTriSg%mixd  =  sg%mixg
       translateAndTurnTriSg%x  =  sg%x  + cx
       translateAndTurnTriSg%y  = -sg%y  + cy
       translateAndTurnTriSg%dx =  sg%dx + cx
       translateAndTurnTriSg%dy = -sg%dy + cy
    case(8)
       translateAndTurnTriSg%mixg  =  sg%mixd
       translateAndTurnTriSg%mixd  =  sg%mixg
       translateAndTurnTriSg%x  =  5.d-1*sg%x  - sqrt3_2d*sg%y  + cx
       translateAndTurnTriSg%y  = -5.d-1*sg%y  - sqrt3_2d*sg%x  + cy
       translateAndTurnTriSg%dx =  5.d-1*sg%dx - sqrt3_2d*sg%dy + cx
       translateAndTurnTriSg%dy = -5.d-1*sg%dy - sqrt3_2d*sg%dx + cy       
    case(9)
       translateAndTurnTriSg%mixg  =  sg%mixd
       translateAndTurnTriSg%mixd  =  sg%mixg
       translateAndTurnTriSg%x  = -5.d-1*sg%x  - sqrt3_2d*sg%y  + cx
       translateAndTurnTriSg%y  =  5.d-1*sg%y  - sqrt3_2d*sg%x  + cy
       translateAndTurnTriSg%dx = -5.d-1*sg%dx - sqrt3_2d*sg%dy + cx
       translateAndTurnTriSg%dy =  5.d-1*sg%dy - sqrt3_2d*sg%dx + cy       
    case(10)
       translateAndTurnTriSg%mixg  =  sg%mixd
       translateAndTurnTriSg%mixd  =  sg%mixg
       translateAndTurnTriSg%x  = -sg%x  + cx
       translateAndTurnTriSg%y  =  sg%y  + cy
       translateAndTurnTriSg%dx = -sg%dx + cx
       translateAndTurnTriSg%dy =  sg%dy + cy
    case(11)
       translateAndTurnTriSg%mixg  =  sg%mixd
       translateAndTurnTriSg%mixd  =  sg%mixg
       translateAndTurnTriSg%x  = -5.d-1*sg%x  + sqrt3_2d*sg%y  + cx
       translateAndTurnTriSg%y  =  5.d-1*sg%y  + sqrt3_2d*sg%x  + cy
       translateAndTurnTriSg%dx = -5.d-1*sg%dx + sqrt3_2d*sg%dy + cx
       translateAndTurnTriSg%dy =  5.d-1*sg%dy + sqrt3_2d*sg%dx + cy       
    case(12)
       translateAndTurnTriSg%mixg  =  sg%mixd
       translateAndTurnTriSg%mixd  =  sg%mixg
       translateAndTurnTriSg%x  =  5.d-1*sg%x  + sqrt3_2d*sg%y  + cx
       translateAndTurnTriSg%y  = -5.d-1*sg%y  + sqrt3_2d*sg%x  + cy
       translateAndTurnTriSg%dx =  5.d-1*sg%dx + sqrt3_2d*sg%dy + cx
       translateAndTurnTriSg%dy = -5.d-1*sg%dy + sqrt3_2d*sg%dx + cy       
    end select
  end function translateAndTurnTriSg

  function interArcArc(clus,ring,aClus,bClus,aRing,bRing,revert)
    type(t_segArc),intent(in)    :: ring,clus
    double precision,intent(out) :: aClus,bClus,aRing,bRing
    logical,intent(out)          :: revert !les arcs sont dans le meme sens?
    integer                      :: interArcArc

    double precision :: d,maxD,minD,aC,aR,d2,rr2,rc2
    integer          :: inClus,inRing
    logical          :: aIn,bIn

    !donne le nombre de points d'intersection entre clus et ring
    !et les angles eventuels d'intersection dans clus et ring
    !avec pour reference l'arc clus (pour les point designes)
    interArcArc = 0
    aClus = 0.d0 ; bClus = 0.d0
    aRing = 0.d0 ; bRing = 0.d0

    d = longVect(ring%x-clus%x,ring%y-clus%y)
    maxD = ring%r + clus%r + epsilon
    minD = abs(ring%r-clus%r) - epsilon
    d2 = d**2 ; rr2 = ring%r**2 ; rc2 = clus%r**2
    revert = (rr2>(d2+rc2)).or.(rc2>(d2+rr2))

    if (d>(maxD+epsilon).or.d<(minD-epsilon)) return !pas d'intersection

    if (isEqualConst(d,maxD).or.isEqualConst(d,minD)) then
       !tangence entre les deux cercles d'appuis
       aClus  = calculeAngle(clus%x,clus%y,ring%x,ring%y)
       if (revert) aClus = angleNormal(aClus - pi_c)
       inRing = isIn(clus%x+clus%r*cos(aClus),clus%y+clus%r*sin(aClus),ring)
       aRing  = calculeAngle(ring%x,ring%y,clus%x,clus%y)
       inClus = isIn(ring%x+ring%r*cos(aRing),ring%y+ring%r*sin(aRing),clus)
       if (inRing==2 .and. inClus==2) interArcArc = 1 !dans les arcs
       return
    end if

    !il y a intersection reelle des deux cercles d'appuis
    ! l'angle a1 adjacent a r1 est determinable a l'aide de la formule
    !  cos(a1) = (d\B2+r1\B2-r2\B2)/(2*d*r1) (idem mutatis-mutandis pour a2)
    aC = acos((d2+rc2-rr2)/(2*d*clus%r))
    aClus = angleNormal(calculeAngle(clus%x,clus%y,ring%x,ring%y) - aC)
    bClus = angleNormal(calculeAngle(clus%x,clus%y,ring%x,ring%y) + aC)
    aR = acos((d2+rr2-rc2)/(2*d*ring%r))
    aRing = angleNormal(calculeAngle(ring%x,ring%y,clus%x,clus%y) - aR)
    bRing = angleNormal(calculeAngle(ring%x,ring%y,clus%x,clus%y) + aR)
    aIn = (isAngleInArc(aClus,clus)==2) .and. (isAngleInArc(aRing,ring)==2)
    bIn = (isAngleInArc(bClus,clus)==2) .and. (isAngleInArc(bRing,ring)==2)
    if (aIn .and. bIn) then
       interArcArc = 2
    else if (aIn) then
       interArcArc = 1
    else if (bIn) then
       interArcArc = 1
       aClus = bClus
       aRing = bRing
    end if
  end function interArcArc

  function giveExtremalsAngles(sa,cx,cy,ao,ae)
    type(t_segArc),intent(in)    :: sa
    double precision,intent(in)  :: cx,cy
    double precision,intent(out) :: ao,ae
    logical                      :: giveExtremalsAngles

    double precision :: ox,oy,ex,ey

    select case(sa%typ)
    case(tcer)
       ! on sort avec une valeur false
       ao = 0.d0 ; ae=0.d0 ; giveExtremalsAngles = .false.
       return
    case(tarc)
       call extremitesArc(sa,ox,oy,ex,ey)
    case(tseg)
       ox = sa%x ; oy = sa%y ; ex = sa%dx ; ey = sa%dy
    end select
    ao = angleNormal(calculeAngle(cx,cy,ox,oy))
    ae = angleNormal(calculeAngle(cx,cy,ex,ey))
    giveExtremalsAngles = .not.( (isEqualConst(cx,ox).and.isEqualConst(cy,oy))&
         .or.(isEqualConst(cx,ex).and.isEqualConst(cy,ey)) )
  end function giveExtremalsAngles

  subroutine drawSegArc(fileNbr,szSA,withNodes,drawMix,zoomx,zoomy)
    integer,intent(in) :: fileNbr,szSA
    logical,intent(in) :: withNodes,drawMix
    real,intent(in)    :: zoomx(2),zoomy(2)

    type(t_segArc) :: sa
    integer        :: i
    real           :: cx,cy,angl,lx,ly,tailleNbr,delx,dely

    g_psp_bBoxXmin = 1.e10
    g_psp_bBoxYmin = 1.e10
    g_psp_bBoxXmax = -1.e10
    g_psp_bBoxYmax = -1.e10
    !recuperation des donnees permettants le centrage
    do i = 1,szSA
       sa = tabSegArc(i)
       if (sa%typ==tseg) then
          g_psp_bBoxXmin = min(real(sa%x),real(sa%dx),g_psp_bBoxXmin)
          g_psp_bBoxYmin = min(real(sa%y),real(sa%dy),g_psp_bBoxYmin)
          g_psp_bBoxXmax = max(real(sa%x),real(sa%dx),g_psp_bBoxXmax)
          g_psp_bBoxYmax = max(real(sa%y),real(sa%dy),g_psp_bBoxYmax)
       else
          g_psp_bBoxXmin = min(real(sa%x-sa%r),g_psp_bBoxXmin)
          g_psp_bBoxYmin = min(real(sa%y-sa%r),g_psp_bBoxYmin)
          g_psp_bBoxXmax = max(real(sa%x+sa%r),g_psp_bBoxXmax)
          g_psp_bBoxYmax = max(real(sa%y+sa%r),g_psp_bBoxYmax)
       end if
    end do
    delx = g_psp_bBoxXmax - g_psp_bBoxXmin
    dely = g_psp_bBoxYmax - g_psp_bBoxYmin
    g_psp_bBoxXmax = g_psp_bBoxXmin + zoomx(2)*delx
    g_psp_bBoxXmin = g_psp_bBoxXmin + zoomx(1)*delx
    g_psp_bBoxYmax = g_psp_bBoxYmin + zoomy(2)*dely
    g_psp_bBoxYmin = g_psp_bBoxYmin + zoomy(1)*dely
    lx = 0.1 * ( g_psp_bBoxXmax - g_psp_bBoxXmin )
    ly = 0.1 * ( g_psp_bBoxYmax - g_psp_bBoxYmin )
    write(*,10) g_psp_bBoxXmin,g_psp_bBoxXmax,g_psp_bBoxYmin,g_psp_bBoxYmax
    10 format(' g2s_segArc: plot domain=(',1p,e14.7,':',e14.7,',',e14.7,':',e14.7,')')
    g_psp_bBoxXmin = g_psp_bBoxXmin - lx
    g_psp_bBoxYmin = g_psp_bBoxYmin - ly
    g_psp_bBoxXmax = g_psp_bBoxXmax + lx
    g_psp_bBoxYmax = g_psp_bBoxYmax + ly
    tailleNbr = min(595./(g_psp_bBoxXmax-g_psp_bBoxXmin), &
         842./(g_psp_bBoxYmax-g_psp_bBoxYmin))
    tailleNbr = 3.6 / tailleNbr
    !impression
    call psinit(fileNbr,.true.)
    do i = 1,szSA
       sa = tabSegArc(i)
       if (sa%typ==tseg) then
          call line(sa%x,sa%y,sa%dx,sa%dy)
          cx=real((sa%dx+sa%x)*0.5d0) ; cy=real((sa%dy+sa%y)*0.5d0)
          angl = real(calculeAngle(sa%x,sa%y,sa%dx,sa%dy)*rad2deg-90.d0)
          if (withNodes .and. drawMix) then
             call keknum(cx,cy,tailleNbr,real(sa%nodeg),angl,-1,2)
             call keknum(cx,cy,tailleNbr,real(sa%noded),angl,-1,0)
          else if (drawMix) then
             call keknum(cx,cy,tailleNbr,real(sa%neutronicMixg),angl,-1,2)
             call keknum(cx,cy,tailleNbr,real(sa%neutronicMixd),angl,-1,0)
          end if
       else if (sa%typ==tcer) then
          call CIRCLE(real(sa%x),real(sa%y),real(sa%r),.false.)
          cx=real(sa%x+sa%r) ; cy=real(sa%y)
          if (withNodes .and. drawMix) then
             call keknum(cx,cy,tailleNbr,real(sa%nodeg),0.,-1,2)
             call keknum(cx,cy,tailleNbr,real(sa%noded),0.,-1,0)
          else if (drawMix) then
             call keknum(cx,cy,tailleNbr,real(sa%neutronicMixg),0.,-1,2)
             call keknum(cx,cy,tailleNbr,real(sa%neutronicMixd),0.,-1,0)
          end if
       else
          call ARC(real(sa%x),real(sa%y),real(sa%r), &
               & real(sa%a*rad2deg),real(sa%b*rad2deg))
          if (sa%b>sa%a) then
             angl=real((sa%b+sa%a)*0.5d0)
          else
             angl=real((sa%b+sa%a)*0.5d0+pi_c)
          end if
          cx=real(sa%x+cos(angl)*sa%r) ; cy=real(sa%y+sin(angl)*sa%r)
          angl=real(angl*rad2deg)
          if (withNodes .and. drawMix) then
             call keknum(cx,cy,tailleNbr,real(sa%nodeg),angl,-1,2)
             call keknum(cx,cy,tailleNbr,real(sa%noded),angl,-1,0)
          else if (drawMix) then
             call keknum(cx,cy,tailleNbr,real(sa%neutronicMixg),angl,-1,2)
             call keknum(cx,cy,tailleNbr,real(sa%neutronicMixd),angl,-1,0)
          end if
       end if
    end do
    call PLOTND()
  end subroutine drawSegArc

  subroutine coupeCercle(indDeb,szSA,nbSeg,typ)
    integer,intent(in)           :: indDeb,nbSeg,typ
    integer,intent(inout)        :: szSA

    type(t_segArc)   :: sg,ar
    double precision :: pt1x,pt1y,pt2x,pt2y,angl1,angl2,midx,midy,tmp
    integer          :: i,j,k,tailleAv,tailleAp,nbPtInter
    type(segArcArrayBis),dimension(:),allocatable :: tmpTabSegArc

    ! pour programmation defensive
    integer          :: taille_table_temp 

    !RQ: indDeb est l'indice du dernier element du tableau qui
    ! ne doit pas etre pris en compte. nbSeg est le nombre
    ! de segments paralleles aux axes, et typ correspond au type (3,4,6,ou...)
    ! CS-IV : Correction des commentaires d'origine
    ! CS-IV : indDeb : indice du dernier element a ne pas prendre en compte 
    ! CS-IV :          (point de depart)
    ! CS-IV : nbSeg : nombre de SA a prendre en charge
    ! CS-IV : szSA : nombre total de SA dans le tableau

    !preparation du tableau de travail
    tailleAp = szSA-indDeb

    ! CS-IV : pourquoi ce 25 ?(on considere que la decoupe ne va pas 
    ! CS-IV : accroitre le nombre de segments de plus de 25 fois
    ! CS-IV : ce coef devrait etre reflechi et dependre du nbr de cercles et
    ! CS-IV : et de segments

    ! pour programmation defensive
    taille_table_temp = 25*tailleAp
    allocate(tmpTabSegArc(1:taille_table_temp),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: coupeCercle (1) => allocation pb")

    !copie des segments et arcs construits dans le tableau
    !temporaire
    ! CS-IV : et suppression dans le tableau initial
    tmpTabSegArc(1:tailleAp)%sa = tabSegArc(indDeb+1:indDeb+tailleAp)
    tmpTabSegArc(1:tailleAp)%keep = .true.

    szSa = indDeb
    ! on coupe tout les elements possible presents, en ne gardant pas ceux qui 
    ! sont a couper, mais seulement les morceaux obtenus
    ! pour i->segment, pour j->arc ou cercle
    do
       tailleAv = tailleAp
       i = 0
       do
          i = i+1
          if (.not. tmpTabSegArc(i)%keep) then
             if (i/=tailleAv ) cycle
             exit
          end if
          sg = tmpTabSegArc(i)%sa
          if (sg%typ/=tseg) then
             if (i/=tailleAv ) cycle
             exit
          end if
          ! CS-IV : le SA tmpTabSegArc(i) est un segment
          do j = 1,tailleAv
             if ((j==i).or.(.not. tmpTabSegArc(j)%keep)) cycle
             ar = tmpTabSegArc(j)%sa
             if (ar%typ==tseg) cycle
             ! CS-IV : le SA tmpTabSegArc(j) est un arc
             nbPtInter = interSgAr(sg,ar,pt1x,pt1y,pt2x,pt2y)
             if (nbPtInter==2) then
                !!deux points d'intersection : 
                ! CS-IV : on obtient 3 segments et 2 arcs
                tmpTabSegArc(i)%keep=.false. ; tmpTabSegArc(j)%keep=.false.
                ! => 3 segments (qui gardent les caracteristiques sect*)
                do k = 1,3
                   tailleAp=tailleAp+1
                   ! pour programmation defensive
                   if (tailleAP > taille_table_temp) &
                        call XABORT("G2S : memory problem in routine CoupeCercle (1)")
                   tmpTabSegArc(tailleAp)%sa=sg
                end do
                tmpTabSegArc(tailleAp-2)%sa%dx=pt1x
                tmpTabSegArc(tailleAp-2)%sa%dy=pt1y
                tmpTabSegArc(tailleAp-1)%sa%x=pt1x
                tmpTabSegArc(tailleAp-1)%sa%y=pt1y
                tmpTabSegArc(tailleAp-1)%sa%dx=pt2x
                tmpTabSegArc(tailleAp-1)%sa%dy=pt2y
                tmpTabSegArc(tailleAp)%sa%x=pt2x
                tmpTabSegArc(tailleAp)%sa%y=pt2y
                tmpTabSegArc(tailleAp-2)%sa%mixg=ar%mixd
                tmpTabSegArc(tailleAp-2)%sa%mixd=ar%mixd
                tmpTabSegArc(tailleAp-1)%sa%mixg=ar%mixg
                tmpTabSegArc(tailleAp-1)%sa%mixd=ar%mixg
                tmpTabSegArc(tailleAp)%sa%mixg=ar%mixd
                tmpTabSegArc(tailleAp)%sa%mixd=ar%mixd
                ! elimination d'un eventuel segment nul
                do k = 0,2
                   tmpTabSegArc(tailleAp-k)%keep= .not. &
                        & ( (isEqual( tmpTabSegArc(tailleAp-k)%sa%x , &
                        &             tmpTabSegArc(tailleAp-k)%sa%dx) ) .and. &
                        &   (isEqual( tmpTabSegArc(tailleAp-k)%sa%y , &
                        &             tmpTabSegArc(tailleAp-k)%sa%dy) ) )
                end do
                angl1=calculeAngle(ar%x,ar%y,pt1x,pt1y)
                angl2=calculeAngle(ar%x,ar%y,pt2x,pt2y)
                !travail sur les arcs
                if (ar%typ==tcer) then
                   !c'etait un cercle complet
                   if (isEqualAngl(angl1,angl2)) then !tangence
                      tailleAp=tailleAp+1
                      ! pour programmation defensive
                      if (tailleAP > taille_table_temp) &
                           call XABORT("G2S : memory problem in routine CoupeCercle (2)")
                      tmpTabSegArc(tailleAp)%sa=ar
                      tmpTabSegArc(tailleAp)%keep=.true.
                      tmpTabSegArc(tailleAp)%sa%typ=tarc
                      tmpTabSegArc(tailleAp)%sa%a=angl1
                      tmpTabSegArc(tailleAp)%sa%b=angl1
                      !positionnement des sect*
                      if (estAGauche(sg%x,sg%y,sg%dx,sg%dy,ar%x,ar%y)) then
                         !cercle completement a gauche de sg
                         tmpTabSegArc(tailleAp)%sa%sectg=sg%sectg
                         tmpTabSegArc(tailleAp)%sa%sectd=sg%sectg
                      else
                         !cercle completement a droite de sg
                         tmpTabSegArc(tailleAp)%sa%sectg=sg%sectd
                         tmpTabSegArc(tailleAp)%sa%sectd=sg%sectd
                      end if
                   else !pas tangence
                      !un arc de chaque cote de sg (1\B0 a drt, 2\B0 a gch)
                      do k = 1,2
                         tailleAp=tailleAp+1
                         ! pour programmation defensive
                         if (tailleAP > taille_table_temp) &
                              call XABORT("G2S : memory problem in routine CoupeCercle (3)")
                         tmpTabSegArc(tailleAp)%sa=ar
                         tmpTabSegArc(tailleAp)%keep=.true.
                         tmpTabSegArc(tailleAp)%sa%typ=tarc
                      end do
                      tmpTabSegArc(tailleAp-1)%sa%a=angl1
                      tmpTabSegArc(tailleAp-1)%sa%b=angl2
                      tmpTabSegArc(tailleAp-1)%sa%sectg=sg%sectd
                      tmpTabSegArc(tailleAp-1)%sa%sectd=sg%sectd
                      tmpTabSegArc(tailleAp)%sa%a=angl2
                      tmpTabSegArc(tailleAp)%sa%b=angl1
                      tmpTabSegArc(tailleAp)%sa%sectg=sg%sectg
                      tmpTabSegArc(tailleAp)%sa%sectd=sg%sectg
                   end if
                else !c'etait deja un arc de cercle
                   do k = 1,3
                      tailleAp=tailleAp+1
                      ! pour programmation defensive
                      if (tailleAP > taille_table_temp) &
                           call XABORT("G2S : memory problem in routine CoupeCercle (4)")
                      tmpTabSegArc(tailleAp)%sa=ar
                   end do
                   tmpTabSegArc(tailleAp-2)%sa%b=angl1
                   !on test si on est bien dans le bon sens
                   if (isIn(pt2x,pt2y,tmpTabSegArc(tailleAp-2)%sa)==2) then
                      !on est a l'envers (ouverture a droite)
                      tmpTabSegArc(tailleAp-2)%sa%b=angl2
                      tmpTabSegArc(tailleAp-1)%sa%a=angl2
                      tmpTabSegArc(tailleAp-1)%sa%b=angl1
                      tmpTabSegArc(tailleAp)%sa%a=angl1
                      !ajustement des sect*
                      tmpTabSegArc(tailleAp-2)%sa%sectg=sg%sectd
                      tmpTabSegArc(tailleAp-2)%sa%sectd=sg%sectd
                      tmpTabSegArc(tailleAp-1)%sa%sectg=sg%sectg
                      tmpTabSegArc(tailleAp-1)%sa%sectd=sg%sectg
                      tmpTabSegArc(tailleAp)%sa%sectg=sg%sectd
                      tmpTabSegArc(tailleAp)%sa%sectd=sg%sectd
                   else !on est a l'endroit (ouverture a gauche)
                      tmpTabSegArc(tailleAp-1)%sa%a=angl1
                      tmpTabSegArc(tailleAp-1)%sa%b=angl2
                      tmpTabSegArc(tailleAp)%sa%a=angl2
                      !ajustement des sect*
                      tmpTabSegArc(tailleAp-2)%sa%sectg=sg%sectg
                      tmpTabSegArc(tailleAp-2)%sa%sectd=sg%sectg
                      tmpTabSegArc(tailleAp-1)%sa%sectg=sg%sectd
                      tmpTabSegArc(tailleAp-1)%sa%sectd=sg%sectd
                      tmpTabSegArc(tailleAp)%sa%sectg=sg%sectg
                      tmpTabSegArc(tailleAp)%sa%sectd=sg%sectg
                   end if
                   do k = 0,2
                      tmpTabSegArc(tailleAp-k)%keep=.not. &
                           & ( isEqualAngl( tmpTabSegArc(tailleAp-k)%sa%a , &
                           &                tmpTabSegArc(tailleAp-k)%sa%b ) )
                   end do
                end if
                i=0
                exit
             else if (abs(nbPtInter)==1) then
                !! un seul point d'intersection
                ! CS-IV : on obtient 2 segments & 2 arcs
                tmpTabSegArc(i)%keep=.false.
                tmpTabSegArc(j)%keep=.false.
                do k = 1,2
                   tailleAp=tailleAp+1
                   ! pour programmation defensive
                   if (tailleAP > taille_table_temp) &
                        call XABORT("G2S : memory problem in routine CoupeCercle (5)")
                   tmpTabSegArc(tailleAp)%sa=sg
                end do
                tmpTabSegArc(tailleAp-1)%sa%dx=pt1x
                tmpTabSegArc(tailleAp-1)%sa%dy=pt1y
                tmpTabSegArc(tailleAp)%sa%x=pt1x
                tmpTabSegArc(tailleAp)%sa%y=pt1y
                if (nbPtInter==1) then
                   !intersection en I (origne de sg a l'exterieur de ar)
                   tmpTabSegArc(tailleAp-1)%sa%mixg=ar%mixd
                   tmpTabSegArc(tailleAp-1)%sa%mixd=ar%mixd
                   tmpTabSegArc(tailleAp)%sa%mixg=ar%mixg
                   tmpTabSegArc(tailleAp)%sa%mixd=ar%mixg
                else
                   !intersection en J (origne de sg a l'interieur de ar)
                   tmpTabSegArc(tailleAp-1)%sa%mixg=ar%mixg
                   tmpTabSegArc(tailleAp-1)%sa%mixd=ar%mixg
                   tmpTabSegArc(tailleAp)%sa%mixg=ar%mixd
                   tmpTabSegArc(tailleAp)%sa%mixd=ar%mixd
                end if
                do k = 0,1
                   tmpTabSegArc(tailleAp-k)%keep=.not. &
                        & ( ( isEqual( tmpTabSegArc(tailleAp-k)%sa%x , &
                        &              tmpTabSegArc(tailleAp-k)%sa%dx ) ) .and.&
                        &   ( isEqual( tmpTabSegArc(tailleAp-k)%sa%y , &
                        &              tmpTabSegArc(tailleAp-k)%sa%dy ) ) )
                end do
                angl1=calculeAngle(ar%x,ar%y,pt1x,pt1y)
                if (ar%typ==tcer) then
                   !ar etait un cercle
                   tailleAp=tailleAp+1
                   ! pour programmation defensive
                   if (tailleAP > taille_table_temp) &
                        call XABORT("G2S : memory problem in routine CoupeCercle (6)")
                   tmpTabSegArc(tailleAp)%sa=ar
                   tmpTabSegArc(tailleAp)%sa%typ=tarc
                   tmpTabSegArc(tailleAp)%sa%a=angl1
                   tmpTabSegArc(tailleAp)%sa%b=angl1
                   tmpTabSegArc(tailleAp)%keep=.true.
                   !pas d'ajustement des sect*, car l'arc n'est ni d'un cote
                   !ni de l'autre
                else
                   do k = 1,2
                      tailleAp=tailleAp+1
                      ! pour programmation defensive
                      if (tailleAP > taille_table_temp) &
                           call XABORT("G2S : memory problem in routine CoupeCercle (7)")
                      tmpTabSegArc(tailleAp)%sa=ar
                   end do
                   tmpTabSegArc(tailleAp-1)%sa%b=angl1
                   tmpTabSegArc(tailleAp)%sa%a=angl1
                   !ajustement des sect*
                   if (nbPtInter==1) then
                      !intersection en I (origne de sg a l'exterieur de ar)
                      tmpTabSegArc(tailleAp-1)%sa%sectg=sg%sectg
                      tmpTabSegArc(tailleAp-1)%sa%sectd=sg%sectg
                      tmpTabSegArc(tailleAp)%sa%sectg=sg%sectd
                      tmpTabSegArc(tailleAp)%sa%sectd=sg%sectd
                   else
                      !intersection en J (origne de sg a l'interieur de ar)
                      tmpTabSegArc(tailleAp-1)%sa%sectg=sg%sectd
                      tmpTabSegArc(tailleAp-1)%sa%sectd=sg%sectd
                      tmpTabSegArc(tailleAp)%sa%sectg=sg%sectg
                      tmpTabSegArc(tailleAp)%sa%sectd=sg%sectg
                   end if
                   do k = 0,1
                      tmpTabSegArc(tailleAp-k)%keep=.not. &
                           & ( isEqualAngl( tmpTabSegArc(tailleAp-k)%sa%a , &
                           &                tmpTabSegArc(tailleAp-k)%sa%b ) )
                   end do
                end if
                i=0
                exit
             end if
          end do
          if (i==tailleAv) exit
       end do
       if ((tailleAp==tailleAv) .and. (i==tailleAv)) exit
    end do

    !on elimine les elements en dehors des limites de la cellule, et on remet le
    !milieu exterieur a la cellule a fooMix
    !Remarque : les elements tmpTabSegArc(1:nbSeg)%sa sont les segments
    !initiaux paralleles aux bords du domaine. Ils peuvent etre a jeter
    !(keep==.false.), mais ils servent tout de meme de reference. Pour savoir
    !si ils sont au bord de la cellule, on teste si le milieu est fooMix.
    !Si oui, on elimine tous les elements qui sont de ce cote.
    if (typ==tRec .or. typ==tHex) then
       do i = 1,nbSeg
          sg = tmpTabSegArc(i)%sa
          do j = 1,tailleAp
             ar = tmpTabSegArc(j)%sa
             if (ar%typ==tseg) then
                if (.not. &
                     (estColi(sg%dx-sg%x,sg%dy-sg%y,ar%dx-ar%x,ar%dy-ar%y) &
                     .and. pointsAlignes(sg%x,sg%y,sg%dx,sg%dy,ar%x,ar%y)) &
                     ) cycle
                if (sg%mixd==fooMix) then
                   tmpTabSegArc(j)%sa%mixd=fooMix
                else if (sg%mixg==fooMix) then
                   tmpTabSegArc(j)%sa%mixg=fooMix
                end if
                cycle
             end if
             if ((j==i) .or. (.not. tmpTabSegArc(j)%keep)) cycle
             if (ar%b>ar%a) then
                angl1=(ar%b+ar%a)*0.5d0
             else
                angl1=(ar%b+ar%a)*0.5d0+pi_c
             end if
             !point milieu de l'arc
             midx=ar%x+ar%r*cos(angl1)
             midy=ar%y+ar%r*sin(angl1)
             !sinus de l'angle (OE,OM)
             tmp = sin( calculeAngle(sg%x,sg%y,midx,midy) - &
                  &     calculeAngle(sg%x,sg%y,sg%dx,sg%dy) )
             if (sg%mixd==fooMix) then !on garde a gauche
                tmpTabSegArc(j)%keep=(tmp>0.).or.(isEqualAngl(ar%a,ar%b).and. &
                     & estAGauche(sg%x,sg%y,sg%dx,sg%dy,ar%x,ar%y))
             else if (sg%mixg==fooMix) then  !on garde a droite
                tmpTabSegArc(j)%keep=(tmp<0.).or.(isEqualAngl(ar%a,ar%b).and. &
                     & estADroite(sg%x,sg%y,sg%dx,sg%dy,ar%x,ar%y))
             end if
          end do
       end do
    else
       call XABORT("G2S : type of geometrie not supported")
    end if

    !elimination des segments nuls eventuellement restants
    do i = 1,tailleAp
       sg = tmpTabSegArc(i)%sa
       if (sg%typ==tseg) tmpTabSegArc(i)%keep=tmpTabSegArc(i)%keep .and. &
            .not.(isEqual(tmpTabSegArc(i)%sa%x,tmpTabSegArc(i)%sa%dx) .and. &
            isEqual(tmpTabSegArc(i)%sa%y,tmpTabSegArc(i)%sa%dy))
    end do

    !on recupere dans le tableau general les elements restants
    do i = 1,tailleAp
       if (.not. tmpTabSegArc(i)%keep) cycle
       szSa = szSa+1
       tabSegArc(szSa) = tmpTabSegArc(i)%sa
    end do

    !on nettoie le tableau temporaire
    deallocate(tmpTabSegArc)
  end subroutine coupeCercle

  subroutine splitSegsForSect(indDeb,szSA)
    integer,intent(in)    :: indDeb
    integer,intent(inout) :: szSA

    type(t_segArc) :: sgi,sgj
    integer        :: i,j

    ! pour programmation defensive
    integer :: taille_table_tabSegArc

    !! permet de spliter les segments en bordure de domaine,
    !! losqu'ils sont interceptes par les segments delimitants
    !! la sectorisation. Si on a intersection, c'est toujours a
    !! l'interieur d'un segment de bordure, et au bout d'un segment
    !! de secteur.                      -->^----->          <-----^<--
    !! Dans ce cas, on split ainsi: (1)  i |  szSA  ou: (2)  szSA | i 

    ! pour programmation defensive
    taille_table_tabSegArc = size(tabSegArc)

    i = indDeb
    do
       i = i + 1
       if (i>szSA) exit
       sgi=tabSegArc(i)
       if (sgi%typ/=tseg) cycle
       do j = indDeb+1,szSA
          if (j==i) cycle
          sgj=tabSegArc(j)
          if (sgj%typ/=tseg) cycle
          if (PointsAlignes(sgi%x,sgi%y,sgi%dx,sgi%dy,sgj%dx,sgj%dy).and.&
               (isIn(sgj%dx,sgj%dy,sgi)==2).and. &
               .not.PointsAlignes(sgi%x,sgi%y,sgi%dx,sgi%dy,sgj%x,sgj%y)) then
             szSA=szSA+1
             ! pour programmation defensive
             if (szSA > taille_table_tabSegArc) &
                  call XABORT("G2S : memory problem in routine splitSegsForSect")
             tabSegArc(szSA) = sgi
             tabSegArc(i)%dx = sgj%dx   ; tabSegArc(i)%dy = sgj%dy
             tabSegArc(szSA)%x = sgj%dx ; tabSegArc(szSA)%y = sgj%dy

             if (estAGaucheStrict(sgj%x,sgj%y,sgj%dx,sgj%dy,sgi%x,sgi%y)) then
                !! cas (1)
                tabSegArc(i)%sectg = sgj%sectg
                tabSegArc(i)%sectd = sgj%sectg
                tabSegArc(szSA)%sectg = sgj%sectd
                tabSegArc(szSA)%sectd = sgj%sectd
             else
                !! cas (2)
                tabSegArc(i)%sectg = sgj%sectd
                tabSegArc(i)%sectd = sgj%sectd
                tabSegArc(szSA)%sectg = sgj%sectg
                tabSegArc(szSA)%sectd = sgj%sectg
             end if
             sgi = tabSegArc(i)
          end if
       end do
    end do
  end subroutine splitSegsForSect

  subroutine majSectori(indDeb,szSA,sectori,typGeo,cx,cy)
    double precision, parameter :: PI = 3.141592653589793d0
    integer,intent(in)          :: indDeb,szSA,sectori,typGeo
    double precision,intent(in) :: cx,cy
    !met a jour les secteurs pour les elements

    double precision :: ao,ae !angles origine et extremite
    integer          :: i,j,nbSect,oIn,eIn
    logical          :: goodAngles
    double precision,dimension(9) :: limit

    if ((typGeo/=tRec).and.(typGeo/=tHex)) &
         call XABORT("G2S: internal error with typGeo in subroutine majSectori")
    !creation des zones angulaires
    nbSect = 0
    select case(sectori)
    case(S_not)
       return !rien a faire
    case(S_X_tot)
       if (typGeo == tRec) then
          nbSect = 4
          limit(1) = pi_2_c*5.d-1 ; limit(2) = pi_2_c*1.5d0
          limit(3) = -limit(2)    ; limit(4) = -limit(1)
          limit(5) = limit(1)
       else
          nbSect = 6
          limit(1) = 0.d0 ; limit(2) = pi_3_c    ; limit(3) = pi_3_c*2.d0
          limit(4) = pi_c ; limit(5) = -limit(3) ; limit(6) = -limit(2)
          limit(7) = limit(1)
       end if
    case(S_T_tot)
       nbSect = 4
       limit(1) = 0.d0 ; limit(2) = pi_2_c
       limit(3) = pi_c ; limit(4) = -limit(2)
       limit(5) = limit(1)
    case(S_TX_tot)
       nbSect = 8
       limit(1) = 0.d0      ; limit(2) = pi_2_c*5.d-1
       limit(3) = pi_2_c    ; limit(4) = pi_2_c*1.5d0
       limit(5) = pi_c      ; limit(6) = -limit(4)
       limit(7) = -limit(3) ; limit(8) = -limit(2)
       limit(9) = limit(1)
    case(S_TXS_tot)
       nbsect = 8
       limit(1) = PI/8. 
       do i=1, 7 
          limit(i+1) = limit(i) + PI/4. 
       enddo
       limit(9) = limit(1) 
    case(S_WM_tot)
       nbsect = 8
       limit(1) = PI/8. 
       do i=1, 7 
          limit(i+1) = limit(i) + PI/4. 
       enddo
       limit(9) = limit(1) 
    case default
       call XABORT("G2S: internal error with sectori in subroutine majSectori")
    end select
    !traitement des elements
    do i = indDeb+1,szSA
       if (tabSegArc(i)%typ == tcer) cycle !on ne s'interesse pas aux cercles
       !on recupere les angles des bouts
       goodAngles = giveExtremalsAngles(tabSegArc(i),cx,cy,ao,ae)
       if (.not.goodAngles) cycle !une des extremites est le centre 
       !  => on ne s'en occupe pas
       do j = 1,nbSect
          !boucle sur les secteurs
          oIn = isAngleInInterval(ao,limit(j),limit(j+1))
          eIn = isAngleInInterval(ae,limit(j),limit(j+1))
          if ((oIn==0) .or. (eIn==0)) then
             !pas dans le secteur
             cycle !on passe au secteur suivant
          else if ((oIn==2) .or. (eIn==2)) then
             !on est completement dans le secteur
             tabSegArc(i)%sectg = j
             tabSegArc(i)%sectd = j
             exit !on sort de la boucle en j
          else if (oIn/=eIn) then
             !element occupant le secteur angulaire total => idem cas precedent
             tabSegArc(i)%sectg = j
             tabSegArc(i)%sectd = j
             exit !on sort de la boucle en j
          else
             !les 2 extremites sont sur le meme separateur
             ! => on ne s'en occupe pas
             exit !on sort de la boucle en j
          end if
       end do
    end do
  end subroutine majSectori

  subroutine addSegsAndClean(sizeSA)
    integer,intent(inout) :: sizeSA

    integer :: i,j,sizeTmp
    integer :: oj,ej !test l'appartenance a sgi des bouts de sgj
    type(t_segArc) :: sgi,sgj
    type(segArcArrayBis),dimension(:),allocatable :: tmpTab

    ! pour programmation defensive
    integer          :: taille_table_temp, taille_table_tabSegArc

    ! copie dans un tableau temporaire des segments et nettoyage du tableau
    ! global
    ! pour programmation defensive
    taille_table_temp = sizeSA*10
    taille_table_tabSegArc = size(tabSegArc)

    allocate(tmpTab(taille_table_temp),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: addSegsAndClean(1) => allocation pb")
    tmpTab(1:sizeSA)%sa=tabSegArc(1:sizeSA)
    tmpTab(1:sizeSA)%keep=.true.

    sizeTmp=sizeSA
    i=0
    do
       i = i+1
       if (i>sizeTmp) exit
       if (.not. tmpTab(i)%keep) cycle
       if (tmpTab(i)%sa%typ/=tseg) cycle
       sgi=tmpTab(i)%sa !c'est bien un segment a garder pour le moment
       j = 0
       do
          j = j+1
          if (j>sizeTmp) exit
          if (j==i) cycle
          if (.not. tmpTab(j)%keep) cycle
          if (tmpTab(j)%sa%typ/=tseg) cycle
          sgj=tmpTab(j)%sa
          if (.not. estColineaire(sgi,sgj)) cycle
          if (.not. pointsAlignes(sgi%x,sgi%y,sgi%dx,sgi%dy,sgj%x,sgj%y)) cycle
          !on renverse le sens de sgj si besoin
          if (.not. isSameWay(sgi,sgj)) sgj=turnBackSide(tmpTab(j)%sa)
          oj=isIn(sgj%x,sgj%y,sgi) ; ej=isIn(sgj%dx,sgj%dy,sgi)
          select case(ej) !seuls cas 0,2,3 utiles pour ej
          case(0) !seuls cas 0,1,2 utiles pour oj
             select case(oj)
             case(0) !si isIn(sgi%x,sgi%y,sgj)/=2 rien sinon, sgi dans sgj
                if (isIn(sgi%x,sgi%y,sgj)/=2) cycle
                call adjustMix(tmpTab(i)%sa,sgj)
                sizeTmp=sizeTmp+1 
                ! pour programmation defensive
                if (sizeTmp > taille_table_temp) &
                     call XABORT("G2S : memory problem in addSegsAndClean (1)")
                tmpTab(sizeTmp)%keep=.true. ; tmpTab(sizeTmp)%sa=sgj
                tmpTab(j)%sa%dx=sgi%x ; tmpTab(j)%sa%dy=sgi%y
                tmpTab(sizeTmp)%sa%x=sgi%dx ; tmpTab(sizeTmp)%sa%y=sgi%dy
             case(1) !origines confondues et sgi dans sgj
                call adjustMix(tmpTab(i)%sa,sgj)
                tmpTab(j)%sa%x=sgi%dx ; tmpTab(j)%sa%y=sgi%dy
             case(2) !sgi et sgj s'intersectent reellement (sgj au dessus)
                sizeTmp=sizeTmp+1
                ! pour programmation defensive
                if (sizeTmp > taille_table_temp) &
                     call XABORT("G2S : memory problem in addSegsAndClean (2)")
                tmpTab(sizeTmp)%keep=.true. ; tmpTab(sizeTmp)%sa=sgi
                tmpTab(i)%sa%dx=sgj%x ; tmpTab(i)%sa%dy=sgj%y
                tmpTab(sizeTmp)%sa%x=sgj%x ; tmpTab(sizeTmp)%sa%y=sgj%y
                tmpTab(sizeTmp)%sa%dx=sgi%dx ; tmpTab(sizeTmp)%sa%dy=sgi%dy
                call adjustMix(tmpTab(sizeTmp)%sa,sgj)
                tmpTab(j)%sa%x=sgi%dx ; tmpTab(j)%sa%y=sgi%dy
             end select
          case(2) !seuls cas 0,1,2 utiles pour oj
             select case(oj)
             case(0) !sgi et sgj s'intersectent reellement (sgj en dessous)
                sizeTmp=sizeTmp+1 
                ! pour programmation defensive
                if (sizeTmp > taille_table_temp) &
                     call XABORT("G2S : memory problem in addSegsAndClean (3)")
                tmpTab(sizeTmp)%keep=.true. ; tmpTab(sizeTmp)%sa=sgi
                tmpTab(i)%sa%x=sgj%dx ; tmpTab(i)%sa%y=sgj%dy
                tmpTab(sizeTmp)%sa%x=sgi%x ; tmpTab(sizeTmp)%sa%y=sgi%y
                tmpTab(sizeTmp)%sa%dx=sgj%dx ; tmpTab(sizeTmp)%sa%dy=sgj%dy
                call adjustMix(tmpTab(sizeTmp)%sa,sgj)
                tmpTab(j)%sa%dx=sgi%x ; tmpTab(j)%sa%dy=sgi%y
             case(1) !origines confondues et sgj dans sgi
                call adjustMix(tmpTab(j)%sa,sgi)
                tmpTab(i)%sa%x=sgj%dx ; tmpTab(i)%sa%y=sgj%dy
             case(2) !sgj dans sgi
                call adjustMix(tmpTab(j)%sa,sgi)
                sizeTmp=sizeTmp+1 
                ! pour programmation defensive
                if (sizeTmp > taille_table_temp) &
                     call XABORT("G2S : memory problem in addSegsAndClean (4)")
                tmpTab(sizeTmp)%keep=.true. ; tmpTab(sizeTmp)%sa=sgi
                tmpTab(i)%sa%dx=sgj%x ; tmpTab(i)%sa%dy=sgj%y
                tmpTab(sizeTmp)%sa%x=sgj%dx ; tmpTab(sizeTmp)%sa%y=sgj%dy
             end select
          case(3) !seuls cas 0,1,2 utiles pour oj
             select case(oj)
             case(0) !extremites confondues et sgi dans sgj
                call adjustMix(tmpTab(i)%sa,sgj)
                tmpTab(j)%sa%dx=sgi%x ; tmpTab(j)%sa%dy=sgi%y
             case(1) !sgi et sgj confondus
                call adjustMix(tmpTab(i)%sa,sgj)
                tmpTab(j)%keep=.false.
             case(2) !extremites confondues et sgj dans sgi
                call adjustMix(tmpTab(j)%sa,sgi)
                tmpTab(i)%sa%dx=sgj%x ; tmpTab(i)%sa%dy=sgj%y
             end select
          end select
       end do
    end do

    do i = 1,sizeTmp
       if (tmpTab(i)%keep) then
          if (tmpTab(i)%sa%typ==tseg) then
             if (isEqual(tmpTab(i)%sa%x,tmpTab(i)%sa%dx).and. &
                  isEqual(tmpTab(i)%sa%y,tmpTab(i)%sa%dy)) tmpTab(i)%keep=.false.
          end if
       end if
    end do

    !recopie dans le tableau global des elements a conserver
    sizeSA=0
    do i = 1,sizeTmp
       if (tmpTab(i)%keep) then
          sizeSA=sizeSA+1
          ! pour programmation defensive
          if (sizeSA > taille_table_tabSegArc) &
               call XABORT("G2S : memory problem in addSegsAndClean (5)")
          tabSegArc(sizeSA) = tmpTab(i)%sa
       end if
    end do
    deallocate(tmpTab)

  end subroutine addSegsAndClean

  subroutine cutClusters(nbRing,nbClus,szSa)
    integer,intent(in)    :: nbRing,nbClus
    integer,intent(inout) :: szSa

    integer          :: i,j,szTmp,nbInter
    double precision :: aClus,bClus,aRing,bRing
    logical          :: revert,eqA,eqB
    type(t_segArc)   :: ring,clus
    type(segArcArrayTer),dimension(:),allocatable :: tmpTabAr

    ! pour programmation defensive
    integer :: taille_table_temp

    taille_table_temp = nbRing*nbClus*4 !ca me semble suffisant (note de l'auteur )
    allocate(tmpTabAr(taille_table_temp),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: cutCluster(1) => allocation pb") 
    szTmp = 0

    !initialisation du tableau temporaire
    !recopie des cluster et des anneaux dans l'ordre inverse de leur
    !creation => du plus grand vers le plus petit
    do i = 1,nbClus
       szTmp = szTmp + 1
       tmpTabAr(szTmp)%keep = .true.
       tmpTabAr(szTmp)%cl = .true.
       tmpTabAr(szTmp)%sa = tabSegArc(szSa)
       szSa = szSa - 1
    end do
    do i = 1,nbRing
       szTmp = szTmp + 1
       tmpTabAr(szTmp)%keep = .true.
       tmpTabAr(szTmp)%cl = .false.
       tmpTabAr(szTmp)%sa = tabSegArc(szSa)
       szSa = szSa - 1
    end do

    i=0
    do
       i = i+1
       if (i>szTmp) exit
       if (.not. (tmpTabAr(i)%keep .and. tmpTabAr(i)%cl)) cycle
       clus = tmpTabAr(i)%sa
       j = 0
       do
          j = j+1
          if (j>szTmp) exit
          if ((.not. tmpTabAr(j)%keep) .or. tmpTabAr(j)%cl) cycle
          ring = tmpTabAr(j)%sa
          nbInter = interArcArc(clus,ring,aClus,bClus,aRing,bRing,revert)
          select case(nbInter)
          case(2)
             if (clus%typ==tcer) then
                !demi cercle interieur
                tmpTabAr(i)%sa%typ = tarc
                tmpTabAr(i)%sa%mixd = ring%mixg
                tmpTabAr(i)%sa%a = aClus
                tmpTabAr(i)%sa%b = bClus
                clus = tmpTabAr(i)%sa
                !demi cercle exterieur
                szTmp = szTmp + 1
                ! programmation defensive
                if (sztmp > taille_table_temp) &
                     call XABORT("G2S : memory problem in routine cutClusters (1)")
                tmpTabAr(szTmp)%keep = .true.
                tmpTabAr(szTmp)%cl = .true.
                tmpTabAr(szTmp)%sa = clus
                tmpTabAr(szTmp)%sa%mixd = ring%mixd
                tmpTabAr(szTmp)%sa%a = bClus
                tmpTabAr(szTmp)%sa%b = aClus                 
             else
                eqA = isEqualAngl(clus%a,aClus)
                eqB = isEqualAngl(clus%b,bClus)
                !debut
                if (.not.eqA) then
                   szTmp = szTmp + 1
                   ! programmation defensive
                   if (sztmp > taille_table_temp) &
                        call XABORT("G2S : memory problem in routine cutClusters (2)")
                   tmpTabAr(szTmp)%keep = .true.
                   tmpTabAr(szTmp)%cl = .true.
                   tmpTabAr(szTmp)%sa = clus
                   tmpTabAr(szTmp)%sa%mixd = ring%mixd
                   tmpTabAr(szTmp)%sa%b = aClus
                end if
                !fin
                if (.not.eqB) then
                   szTmp = szTmp + 1
                   ! programmation defensive
                   if (sztmp > taille_table_temp) &
                        call XABORT("G2S : memory problem in routine cutClusters (3)")
                   tmpTabAr(szTmp)%keep = .true.
                   tmpTabAr(szTmp)%cl = .true.
                   tmpTabAr(szTmp)%sa = clus
                   tmpTabAr(szTmp)%sa%mixd = ring%mixd
                   tmpTabAr(szTmp)%sa%a = bClus
                end if
                !milieu
                if ((.not.eqA).or.(.not.eqB)) then
                   tmpTabAr(i)%sa%a = aClus
                   tmpTabAr(i)%sa%b = bClus
                   tmpTabAr(i)%sa%mixd = ring%mixg
                   clus = tmpTabAr(i)%sa
                end if
             end if
             if (ring%typ==tcer) then
                tmpTabAr(j)%sa%typ = tarc
                tmpTabAr(j)%sa%a = bRing
                tmpTabAr(j)%sa%b = aRing
                ring = tmpTabAr(j)%sa
             else
                tmpTabAr(j)%keep = .false.
                eqA = isEqualAngl(ring%a,aRing)
                eqB = isEqualAngl(ring%b,bRing)
                !premier bout d'arc
                if (.not.eqA) then
                   szTmp = szTmp + 1
                   ! programmation defensive
                   if (sztmp > taille_table_temp) &
                        call XABORT("G2S : memory problem in routine cutClusters (4)")
                   tmpTabAr(szTmp)%keep = .true.
                   tmpTabAr(szTmp)%cl = .false.
                   tmpTabAr(szTmp)%sa = ring
                   tmpTabAr(szTmp)%sa%b = aRing
                end if
                !second bout d'arc
                if (.not.eqB) then
                   szTmp = szTmp + 1
                   ! programmation defensive
                   if (sztmp > taille_table_temp) &
                        call XABORT("G2S : memory problem in routine cutClusters (5)")
                   tmpTabAr(szTmp)%keep = .true.
                   tmpTabAr(szTmp)%cl = .false.
                   tmpTabAr(szTmp)%sa = ring
                   tmpTabAr(szTmp)%sa%a = bRing
                end if
             end if
          case(1)
             !a priori, tel que l'algo est fait, le cas ne se produit que si
             !le cluster est tangent a un anneau
             !travail sur l'anneau
             if (ring%typ==tcer) then
                tmpTabAr(j)%sa%typ = tarc
                tmpTabAr(j)%sa%a = aRing
                tmpTabAr(j)%sa%b = aRing
                ring = tmpTabAr(j)%sa
             else
                tmpTabAr(j)%keep = .false.
                eqA = isEqualAngl(ring%a,aRing)
                eqB = isEqualAngl(ring%b,aRing)
                !premier bout d'arc
                if (.not.eqA) then
                   szTmp = szTmp + 1
                   ! programmation defensive
                   if (sztmp > taille_table_temp) &
                        call XABORT("G2S : memory problem in routine cutClusters (6)")
                   tmpTabAr(szTmp)%keep = .true.
                   tmpTabAr(szTmp)%cl = .false.
                   tmpTabAr(szTmp)%sa = ring
                   tmpTabAr(szTmp)%sa%b = aRing
                end if
                !second bout d'arc
                if (.not.eqB) then
                   szTmp = szTmp + 1
                   ! programmation defensive
                   if (sztmp > taille_table_temp) &
                        call XABORT("G2S : memory problem in routine cutClusters (7)")
                   tmpTabAr(szTmp)%keep = .true.
                   tmpTabAr(szTmp)%cl = .false.
                   tmpTabAr(szTmp)%sa = ring
                   tmpTabAr(szTmp)%sa%a = aRing
                end if
             end if
             !travail sur le cluster
             if (clus%typ==tcer) then
                tmpTabAr(i)%sa%typ = tarc
                tmpTabAr(i)%sa%a = aClus
                tmpTabAr(i)%sa%b = aClus
                clus = tmpTabAr(i)%sa
                if (revert) then !cluster a l'exterieur de l'anneau
                   tmpTabAr(i)%sa%mixd = ring%mixd
                else              !cluster a l'interieur de l'anneau
                   tmpTabAr(i)%sa%mixd = ring%mixg
                end if
             else ! le cluster a deja ete coupe => tangence exterieur
                !   et derniere utilisation
                eqA = isEqualAngl(clus%a,aClus)
                eqB = isEqualAngl(clus%b,aClus)
                !debut
                if (.not.eqA) then
                   szTmp = szTmp + 1
                   ! programmation defensive
                   if (sztmp > taille_table_temp) &
                        call XABORT("G2S : memory problem in routine cutClusters (8)")
                   tmpTabAr(szTmp)%keep = .true.
                   tmpTabAr(szTmp)%cl = .true.
                   tmpTabAr(szTmp)%sa = clus
                   tmpTabAr(szTmp)%sa%mixd = ring%mixd
                   tmpTabAr(szTmp)%sa%b = aClus
                end if
                !fin
                if (.not.eqB) then
                   szTmp = szTmp + 1
                   ! programmation defensive
                   if (sztmp > taille_table_temp) &
                        call XABORT("G2S : memory problem in routine cutClusters (9)")
                   tmpTabAr(szTmp)%sa%a = aClus
                   tmpTabAr(szTmp)%sa = clus
                   tmpTabAr(szTmp)%sa%mixd = ring%mixd
                   tmpTabAr(szTmp)%sa%a = aClus
                end if
                !milieu
                if ((.not.eqA).or.(.not.eqB)) then
                   tmpTabAr(i)%keep = .false.
                   i = 0
                end if
             end if

          end select
       end do
    end do
    !recopie des resultats
    do i = 1,szTmp
       if (.not. tmpTabAr(i)%keep) cycle
       szSa = szSa + 1
       tabSegArc(szSa) = tmpTabAr(i)%sa
    end do
    deallocate(tmpTabAr)
  end subroutine cutClusters

  logical function segWithSameCoord(S1,S2,n)
    type(t_segArc), intent(in) :: S1, S2
    integer, intent(out) :: n

    segWithSameCoord = .false. ; n=0
    if (isEqualConst(S1%x,S2%x)) then
       if (isEqualConst(S1%dx,S2%dx)) then
          if (isEqualConst(S1%y,S2%y)) then
             if (isEqualConst(S1%dy,S2%dy)) then
                segWithSameCoord = .true. ; n=1
             endif
          endif
       endif
    elseif (isEqualConst(S1%x,S2%dx)) then
       if (isEqualConst(S1%dx,S2%x)) then
          if (isEqualConst(S1%y,S2%dy)) then
             if (isEqualConst(S1%dy,S2%y)) then
                segWithSameCoord = .true. ; n=-1
             endif
          endif
       endif
    endif
  end function segWithSameCoord

  subroutine PrintTabSegArc(szSA)
    integer, intent(in) :: szSA
    type(t_segArc) :: sa
    integer :: i

    write(*,10) szSA

    do i = 1, szSA
       sa = tabSegArc(i)
       if (sa%typ == tseg) then
          write(*,20) i, sa%x, sa%y, sa%dx, sa%dy
       elseif (sa%typ == tarc) then
          write(*,30) i, sa%x, sa%y, sa%r, sa%a, sa%b
       elseif (sa%typ == tcer) then
          write(*,35) i, sa%x, sa%y, sa%r
       else
          write(*,38) i, sa%typ, sa%x, sa%y
       endif
       write(*,40) sa%sectg, sa%sectd
       write(*,50) sa%mixg, sa%mixd
       write(*,60) sa%nodeg, sa%noded
       write(*,70) sa%indcellpg, sa%indcellpd
       write(*,80) sa%neutronicMixg,sa%neutronicMixd
    enddo

10  format("Number of SegArc :", i4)
20  format("N. SegArc ", i6," segment type : ",/, &
         "|-> origin/extrem : (", f7.4, ";", f7.4,")/(", f7.4, ";", f7.4,")")
30  format("N. SegArc ", i6," arc type:",/, &
         "|-> center/radius/angles : (", f7.4,";",f7.4,")/", f7.4,"/",f7.4,";", f7.4)
35  format("N. SegArc ", i6," circle type :",/, &
         "|-> center/radius : (", f7.4,";",f7.4,")/", f7.4)
38  format("N. SegArc ", i6," unknown type :",/, &
         "|-> type/origin/ : (", i6,"/",f7.4,";",f7.4,")")
40  format( "|-> sectg sectd  : ", i6,2x,i6)
50  format( "|-> mixg mixd  : ", i6,2x,i6)
60  format( "|-> nodeg noded  : ", i6,2x,i6)
70  format( "|-> IndCellPg IndCellPd  : ", i6,2x,i6)
80  format( "|-> neutronicMixg neutronicMixd  : ", i6,2x,i6)
  end subroutine PrintTabSegArc

  subroutine PrintSegArc1(i)
    integer, intent(in) ::i
    type(t_segArc) :: sa

    sa = tabSegArc(i)
    if (sa%typ == tseg) then
       write(*,20) i, sa%x, sa%y, sa%dx, sa%dy
    elseif (sa%typ == tarc) then
       write(*,30) i, sa%x, sa%y, sa%r, sa%a, sa%b
    elseif (sa%typ == tcer) then
       write(*,35) i, sa%x, sa%y, sa%r
    else
       write(*,38) i, sa%typ, sa%x, sa%y
    endif
    write(*,40) sa%sectg, sa%sectd
    write(*,50) sa%mixg, sa%mixd
    write(*,60) sa%nodeg, sa%noded
    write(*,70) sa%indcellpg, sa%indcellpd
    write(*,80) sa%neutronicMixg,sa%neutronicMixd

20  format("N. SegArc ", i6," segment type : ",/, &
         "|-> origin/extrem : (", f7.4, ";", f7.4,")/(", f7.4, ";", f7.4,")")
30  format("N. SegArc ", i6," arc type:",/, &
         "|-> center/radius/angles : (", f7.4,";",f7.4,")/", f7.4,"/",f7.4,";", f7.4)
35  format("N. SegArc ", i6," circle type :",/, &
         "|-> center/radius : (", f7.4,";",f7.4,")/", f7.4)
38  format("N. SegArc ", i6," unknown type :",/, &
         "|-> type/origin/ : (", i6,"/",f7.4,";",f7.4,")")
40  format( "|-> sectg sectd  : ", i6,2x,i6)
50  format( "|-> mixg mixd  : ", i6,2x,i6)
60  format( "|-> nodeg noded  : ", i6,2x,i6)
70  format( "|-> IndCellPg IndCellPd  : ", i6,2x,i6)
80  format( "|-> neutronicMixg neutronicMixd  : ", i6,2x,i6)

  end subroutine PrintSegArc1

  subroutine PrintSegArc2(sa)
    type(t_segArc), intent(in) :: sa

    if (sa%typ == tseg) then
       write(*,20) sa%x, sa%y, sa%dx, sa%dy
    elseif (sa%typ == tarc) then
       write(*,30) sa%x, sa%y, sa%r, sa%a, sa%b
    elseif (sa%typ == tcer) then
       write(*,35) sa%x, sa%y, sa%r
    else
       write(*,38) sa%typ, sa%x, sa%y
    endif
    write(*,40) sa%sectg, sa%sectd
    write(*,50) sa%mixg, sa%mixd
    write(*,60) sa%nodeg, sa%noded
    write(*,70) sa%indcellpg, sa%indcellpd
    write(*,80) sa%neutronicMixg,sa%neutronicMixd

20  format(" Segment type : ",/, &
         "|-> origin/extrem : (", f7.4, ";", f7.4,")/(", f7.4, ";", f7.4,")")
30  format(" Arc type:",/, &
         "|-> center/radius/angles : (", f7.4,";",f7.4,")/", f7.4,"/",f7.4,";", f7.4)
35  format("Circle type :",/, &
         "|-> center/radius : (", f7.4,";",f7.4,")/", f7.4)
38  format("unknown type :",/, &
         "|-> type/origin/ : (", i6,"/",f7.4,";",f7.4,")")
40  format( "|-> sectg sectd  : ", i6,2x,i6)
50  format( "|-> mixg mixd  : ", i6,2x,i6)
60  format( "|-> nodeg noded  : ", i6,2x,i6)
70  format( "|-> IndCellPg IndCellPd  : ", i6,2x,i6)
80  format( "|-> neutronicMixg neutronicMixd  : ", i6,2x,i6)

  end subroutine PrintSegArc2



end module segArc
