!
!-----------------------------------------------------------------------
!
!Purpose:
! Collect the construction functions for different data types used in
! G2S: module.
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
!  - construit_xxxxx : cree un xxxxx, avec xxxxx parmi car2d, carcel, hex,
!                      hexcel, tri2d ou tube
!
!-----------------------------------------------------------------------
!
module construire
  
  use celluleBase
  use constType
  use constUtiles
  use segArc
  
  implicit none
  
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: tg = 0.414213562373095d0

contains

  subroutine construit_car2d(cx,cy,sx,sy,turn,inMeshx,inMeshy,splitx,splity, &
       inMix,szSA)
    double precision,intent(in)              :: cx,cy,sx,sy
    double precision,dimension(:),intent(in) :: inMeshx,inMeshy !entre 0 et 1
    integer,intent(in)                       :: turn
    integer,dimension(:), intent(in)         :: splitx,splity,inMix
    integer,intent(inout)                    :: szSA

    double precision                          :: xmin,xmax,ymin,ymax,ssx,ssy,dd
    double precision,dimension(:),allocatable :: mmeshx,mmeshy,meshx,meshy
    integer,dimension(:),allocatable          :: mmix,mix,mixDiag
    integer                                   :: i,j,dimx,dimy,smix,tmpSzx
    integer                                   :: ind,lgx,lgy,tmpSzy
    logical                                   :: ldiag
    type(t_segArc)                            :: sg

    ! programmation defensive
    integer :: dimtabSegArc 
    dimtabSegArc = size(tabSegArc)

    !prise en compte du split => creation de meshx,meshy et mix a partir des
    !donnees d'entree

    lgx = size(inMeshx)-1 ; lgy = size(inMeshy)-1
    ldiag=(lgx == lgy).and.(lgx*(lgx+1)/2 == size(inMix))
    tmpSzx = sum(splitx(1:lgx)) ; tmpSzy = sum(splity(1:lgy)) ; smix = tmpSzx*tmpSzy
    allocate(meshx(tmpSzx+1),meshy(tmpSzy+1),mix(smix),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: construit_car2d(1) => allocation pb")

    !sur x
    meshx(1) = 0.d0
    ind = 1
    do i = 1,lgx
       dd = (inMeshx(i+1)-inMeshx(i)) / splitx(i)
       do j = 1,splitx(i)-1
          ind = ind + 1
          meshx(ind) = meshx(ind-1) + dd
       end do
       ind = ind + 1
       meshx(ind) = inMeshx(i+1)
    end do
    !sur y
    meshy(1) = 0.d0
    ind = 1
    do i = 1,lgy
       dd = (inMeshy(i+1)-inMeshy(i)) / splity(i)
       do j = 1,splity(i)-1
          ind = ind + 1
          meshy(ind) = meshy(ind-1) + dd
       end do
       ind = ind + 1
       meshy(ind) = inMeshy(i+1)
    end do
    !mix
    ! correction of plain CAR2D bug by Alain Hebert (May 2016)
    if(ldiag) then
       allocate(mixDiag(lgx*lgx),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: construit_car2d(2) => allocation pb")
       ind = 0
       do j = 1,lgx
          do i = j,lgx
             ind=ind+1
             if(ind > size(inMix)) call XABORT('overflow1')
             if(i+(j-1)*lgx > lgx*lgx) call XABORT('overflow2')
             if(j+(i-1)*lgx > lgx*lgx) call XABORT('overflow3')
             mixDiag(i+(j-1)*lgx)=inMix(ind)
             mixDiag(j+(i-1)*lgx)=inMix(ind)
          end do
       end do
    endif
    ind = 1
    do j = 1,lgy
       do i = 1,lgx
          if(ldiag) then
             mix(ind:ind+splitx(i)-1) = mixDiag(i+(j-1)*lgx)
          else
             mix(ind:ind+splitx(i)-1) = inMix(i+(j-1)*lgx)
          endif
          ind = ind + splitx(i)
       end do
       do i = 2,splity(j)
          mix(ind:ind+tmpSzx-1) = mix(ind-tmpSzx:ind-1)
          ind = ind + tmpSzx
       end do
    end do
    if(ldiag) deallocate(mixDiag)

    !creation des tableaux de travail pour les abscisses, les ordonnees et les
    !milieux, en fonction du turn (la construction des segments se faisant
    !ensuite toujours de la meme maniere que pour un turn = 1
    if(mod(turn,2) == 0) then
       ssx = sy ; ssy = sx ; dimx = tmpSzy ; dimy = tmpSzx
    else
       ssx = sx ; ssy = sy ; dimx = tmpSzx ; dimy = tmpSzy
    end if
    allocate(mmeshx(dimx+1),mmeshy(dimy+1),mmix(smix),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: construit_car2d(3) => allocation pb")
    select case(turn)
    case(1)
       mmeshx(:dimx+1) = meshx(:dimx+1)      ; mmeshy(:dimy+1) = meshy(:dimy+1)
       mmix(:smix) = mix(:smix)
    case(2)
       mmeshx(:dimx+1) = meshy(:dimx+1)
       do i = 1,dimy+1
          mmeshy(i) = 1.-meshx(dimy+2-i)
       end do
       mmix(:smix) = (/(mix(i:i+(dimx-1)*dimy:dimy),i=dimy,1,-1)/)
    case(3)
       do i = 1,dimx+1
          mmeshx(i) = 1.-meshx(dimx+2-i)
       end do
       do i = 1,dimy+1
          mmeshy(i) = 1.-meshy(dimy+2-i)
       end do
       mmix(:smix) = mix(dimx*dimy:1:-1)
    case(4)
       do i = 1,dimx+1
          mmeshx(i) = 1.-meshy(dimx+2-i)
       end do
       mmeshy(:dimy+1) = meshx(:dimy+1)
       mmix(:smix) = (/(mix(i:i-(dimx-1)*dimy:-dimy),i=1+(dimx-1)*dimy,dimx*dimy,1)/)
    case(5)
       do i = 1,dimx+1
          mmeshx(i) = 1.-meshx(dimx+2-i)
       end do
       mmeshy(:dimy+1) = meshy(:dimy+1)
       mmix(:smix) = (/(mix(i:i-(dimx-1):-1),i=dimx,dimx*dimy,dimx)/)
    case(6)
       mmeshx(:dimx+1) = meshy(:dimx+1)
       mmeshy(:dimy+1) = meshx(:dimy+1)
       mmix(:smix) = (/(mix(i:i+(dimx-1)*dimy:dimy),i=1,dimy,1)/)
    case(7)
       mmeshx(:dimx+1) = meshx(:dimx+1)
       do i = 1,dimy+1
          mmeshy(i) = 1.-meshy(dimy+2-i)
       end do
       mmix(:smix) = (/(mix(i:i+(dimx-1):1),i=1+(dimy-1)*dimx,1,-dimx)/)
    case(8)
       do i = 1,dimx+1
          mmeshx(i) = 1.-meshy(dimx+2-i)
       end do
       do i = 1,dimy+1
          mmeshy(i) = 1.-meshx(dimy+2-i)
       end do
       mmix(:smix) = (/(mix(i:i-(dimx-1)*dimy:-dimy),i=dimx*dimy,1+(dimx-1)*dimy,-1)/)
    end select
    deallocate(meshx,meshy,mix)

    !construction des segments :
    !elle se fait dans l'odre suivant:     ...       n-1   n
    !                                     ---->     ---->---->
    !                                     ....................
    !                                       1    2        ...
    !                                     ---->---->     ---->
    !puis:      |. :   |p
    !              :   |p-1
    !           |2 :
    !           |1 :   |.
    !Remarque: on distingue le cas particulier des bords, pour l'attribution
    !du milieux fooMix
    xmin = cx-0.5d0*ssx
    xmax = xmin+ssx
    ymin = cy-0.5d0*ssy
    ymax = ymin+ssy
    do i=1,dimx
       !cas particulier j=1
       szSA=szSA+1
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine construit_car2d (1)")
       sg=createSeg(xmin+mmeshx(i)*ssx, ymin, xmin+mmeshx(i+1)*ssx, ymin   , &
            & mmix(i)             , fooMix )
       tabSegArc(szSA)=sg
       do j=2,dimy
          szSA=szSA+1
          ! programmation defensive
          if (szSA > dimTabSegArc)&
               call XABORT("G2S: memory problem in routine construit_car2d (2)")
          sg=createSeg(xmin+mmeshx(i)*ssx  ,ymin+mmeshy(j)*ssy, &
               & xmin+mmeshx(i+1)*ssx,ymin+mmeshy(j)*ssy, &
               & mmix((j-1)*dimx+i)  ,mmix((j-2)*dimx+i) )
          tabSegArc(szSA)=sg
       end do
       !cas particulier j=dimy+1
       szSA=szSA+1
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine construit_car2d (3)")

       sg=createSeg(xmin+mmeshx(i)*ssx  ,ymax                , &
            & xmin+mmeshx(i+1)*ssx,ymax                , &
            & fooMix              ,mmix((dimy-1)*dimx+i))
       tabSegArc(szSA)=sg
    end do
    do j=1,dimy
       !cas particulier i=1
       szSA=szSA+1
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine construit_car2d (4)")
       sg=createSeg(xmin   ,ymin+mmeshy(j)*ssy   , &
            & xmin   ,ymin+mmeshy(j+1)*ssy , &
            & fooMix ,mmix((j-1)*dimx+1)    )
       tabSegArc(szSA)=sg
       do i=2,dimx
          szSA=szSA+1
          ! programmation defensive
          if (szSA > dimTabSegArc)&
               call XABORT("G2S: memory problem in routine construit_car2d (5)")
          sg=createSeg(xmin+mmeshx(i)*ssx  ,ymin+mmeshy(j)*ssy  , &
               & xmin+mmeshx(i)*ssx  ,ymin+mmeshy(j+1)*ssy, &
               & mmix((j-1)*dimx+i-1),mmix((j-1)*dimx+i)   )
          tabSegArc(szSA)=sg
       end do
       !cas particulier i=dimx+1
       szSA=szSA+1
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine construit_car2d (6)")
       sg=createSeg(xmax        ,ymin+mmeshy(j)*ssy  , &
            & xmax        ,ymin+mmeshy(j+1)*ssy, &
            & mmix(j*dimx),fooMix              )
       tabSegArc(szSA)=sg
    end do
    deallocate(mmeshx,mmeshy,mmix)
    !call  PrintTabSegArc(szSA)
  end subroutine construit_car2d

  subroutine construit_carcel(cx,cy,sx,sy,turn,radius,offcx,offcy, &
       & splitx,splity,mix,sectori,sectorj,cluster,szSA)
    double precision,intent(in)           :: cx,cy,sx,sy,offcx,offcy
    double precision,dimension(:),intent(in) :: radius
    integer,intent(in)                    :: turn,sectori,sectorj
    integer,dimension(:),intent(in)       :: splitx,splity,mix
    type(t_cluster),dimension(:),pointer  :: cluster
    integer,intent(inout)                 :: szSA

    integer                       :: i,j,keepSz,nbseg,tmp,nb
    double precision              :: ccx,ccy,d,dorig,dextr,tmpx,tmpy
    double precision             :: tgx, tgy
    double precision              :: pt1x,pt1y,pt2x,pt2y,ss,ssx,ssy
    double precision,dimension(8) :: coef,coef2
    integer,dimension(:),allocatable  :: tmpMix
    type(t_segArc)                :: sg,ar,tmpSg
    double precision              :: dist1, dist2

    ! programmation defensive
    integer :: dimtabSegArc
    dimtabSegArc = size(tabSegArc)
    keepSz=szSA
    select case (turn)
    case(1) ;       ccx=cx+offcx ; ccy=cy+offcy
    case(2) ;       ccx=cx+offcy ; ccy=cy-offcx
    case(3) ;       ccx=cx-offcx ; ccy=cy-offcy
    case(4) ;       ccx=cx-offcy ; ccy=cy+offcx
    case(5) ;       ccx=cx-offcx ; ccy=cy+offcy
    case(6) ;       ccx=cx+offcy ; ccy=cy+offcx
    case(7) ;       ccx=cx+offcx ; ccy=cy-offcy
    case(8) ;       ccx=cx-offcy ; ccy=cy-offcx
    end select

    if(mod(turn,2)==0) then ; ssx = sy ; ssy = sx
    else ; ssx = sx ; ssy = sy ;  end if

    !! creation du cadre
    allocate(tmpMix(1))
    tmpMix(1) = mix(size(mix))
    call construit_car2d(cx,cy,sx,sy,turn,(/0.d0,1.d0/),(/0.d0,1.d0/), &
         splitx,splity,tmpMix,szSA)
    deallocate(tmpMix)
    nbseg = 2*splitx(1)*splity(1) + splitx(1) + splity(1)
    if(nbseg/=4 .and. sectori/=S_not) call XABORT("G2S: mixing SPLIT and SECT&
         & not allowed")

    !! gestion des secteurs
    ss = min(sx,sy)
    if(sectori==s_not) then
       !! rien
    else if((sectori==S_X_tot)) then
       tmp = mix(size(mix))
       coef = 0.5d0 * (/1.d0,-1.d0,-1.d0, 1.d0,&
                        1.d0, 1.d0,-1.d0,-1.d0/) * ss
       !!   2 1
       !!    X
       !!   3 4
       if (sectorj==0) then
          do i = 1,4
             nbseg=nbseg+1
             szSA=szSA+1
             ! programmation defensive
             if (szSA > dimTabSegArc)&
                  call XABORT("G2S: memory problem in routine &
                  &construit_carcel (1)")
             tmpSg = createSeg(cx,cy,cx+coef(i),cy+coef(i+4),tmp,tmp)
             tmpSg%sectg = i
             tmpSg%sectd = mod(i+2,4)+1
             tabSegArc(szSA) = tmpSg
          end do
       else
          ar = createArc(ccx,ccy,radius(sectorj+1),0.d0,0.d0,tmp,tmp)
          do i = 1,4
             sg = createSeg(cx,cy,cx+coef(i),cy+coef(i+4),tmp,tmp)
             nb = abs(interSgAr(sg,ar,pt1x,pt1y,pt2x,pt2y))
             if(nb==1) then
                if(.not.(isEqualConst(pt1x,cx+coef(i)).and. &
                     isEqualConst(pt1y,cy+coef(i+4)))) then
                   nbseg=nbseg+1
                   szSA=szSA+1
                   ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine &
                        &construit_carcel (1)")
                   tmpSg = createSeg(pt1x,pt1y,cx+coef(i),cy+coef(i+4),tmp,tmp)
                   tmpSg%sectg = i
                   tmpSg%sectd = mod(i+2,4)+1
                   tabSegArc(szSA) = tmpSg
                end if
             else if(nb==2) then
                if(.not.(isEqualConst(pt2x,cx+coef(i)).and. &
                     isEqualConst(pt2y,cy+coef(i+4)))) then
                   nbseg=nbseg+1
                   szSA=szSA+1
                   ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine &
                        &construit_carcel (3)")
                   tmpSg = createSeg(pt2x,pt2y,cx+coef(i),cy+coef(i+4),tmp,tmp)
                   tmpSg%sectg = i
                   tmpSg%sectd = mod(i+2,4)+1
                   tabSegArc(szSA) = tmpSg
                end if
             end if
          end do
       end if
    else if(sectori==S_T_tot) then
       tmp=mix(size(mix))
       coef=0.5d0*(/ ssx,0.d0,-ssx,0.d0,&
            0.d0, ssy,0.d0,-ssy/)
       !!    2
       !!   3+1
       !!    4
       if (sectorj==0) then
          do i = 1,4
             nbseg=nbseg+1
             szSA=szSA+1
             ! programmation defensive
             if (szSA > dimTabSegArc)&
                  call XABORT("G2S: memory problem in routine &
                  &construit_carcel (4)")
             tmpSg = createSeg(cx,cy,cx+coef(i),cy+coef(i+4),tmp,tmp)
             tmpSg%sectg = i
             tmpSg%sectd = mod(i+2,4)+1
             tabSegArc(szSA) = tmpSg
          end do
       else
          ar = createArc(ccx,ccy,radius(sectorj+1),0.d0,0.d0,tmp,tmp)
          do i = 1,4
             sg = createSeg(cx,cy,cx+coef(i),cy+coef(i+4),tmp,tmp)
             nb = abs(interSgAr(sg,ar,pt1x,pt1y,pt2x,pt2y))
             if(nb==1) then
                if(.not.(isEqualConst(pt1x,cx+coef(i)).and. &
                     isEqualConst(pt1y,cy+coef(i+4)))) then
                   nbseg=nbseg+1
                    ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine &
                        &construit_carcel (5)")
                   szSA=szSA+1
                   tmpSg = createSeg(pt1x,pt1y,cx+coef(i),cy+coef(i+4),tmp,tmp)
                   tmpSg%sectg = i
                   tmpSg%sectd = mod(i+2,4)+1
                   tabSegArc(szSA) = tmpSg
                end if
             else if(nb==2) then
                if(.not.(isEqualConst(pt2x,cx+coef(i)).and. &
                     isEqualConst(pt2y,cy+coef(i+4)))) then
                   nbseg=nbseg+1
                   szSA=szSA+1
                   ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine &
                        &construit_carcel (6)")
                   tmpSg = createSeg(pt2x,pt2y,cx+coef(i),cy+coef(i+4),tmp,tmp)
                   tmpSg%sectg = i
                   tmpSg%sectd = mod(i+2,4)+1
                   tabSegArc(szSA) = tmpSg
                end if
             end if
          end do
       end if
    else if(sectori==S_TX_tot) then
       tmp=mix(size(mix))
       coef  = 0.5d0 * (/ ssx,ss,0.d0,-ss,-ssx,-ss,0.d0, ss/)
       coef2 = 0.5d0 * (/0.d0,ss, ssy, ss,0.d0,-ss,-ssy,-ss/)
       !!   432
       !!   5*1
       !!   678   
       if (sectorj==0) then
          do i = 1,8
             nbseg=nbseg+1
             szSA=szSA+1
             ! programmation defensive
             if (szSA > dimTabSegArc)&
                  call XABORT("G2S: memory problem in routine &
                  &construit_carcel (7)")
             tmpSg = createSeg(cx,cy,cx+coef(i),cy+coef2(i),tmp,tmp)
             tmpSg%sectg = i
             tmpSg%sectd = mod(i+6,8)+1
             tabSegArc(szSA) = tmpSg
          end do
       else 
          ar = createArc(ccx,ccy,radius(sectorj+1),0.d0,0.d0,tmp,tmp)
          do i = 1,8
             sg = createSeg(cx,cy,cx+coef(i),cy+coef2(i),tmp,tmp)
             nb = abs(interSgAr(sg,ar,pt1x,pt1y,pt2x,pt2y))
             if(nb==1) then
                if(.not.(isEqualConst(pt1x,cx+coef(i)).and. &
                     isEqualConst(pt1y,cy+coef2(i)))) then
                   nbseg=nbseg+1
                   szSA=szSA+1
                   ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine &
                        &construit_carcel (7)")
                   tmpSg = createSeg(pt1x,pt1y,cx+coef(i),cy+coef2(i),tmp,tmp)
                   tmpSg%sectg = i
                   tmpSg%sectd = mod(i+6,8)+1
                   tabSegArc(szSA) = tmpSg
                end if
             else if(nb==2) then
                if(.not.(isEqualConst(pt2x,cx+coef(i)).and. &
                     isEqualConst(pt2y,cy+coef2(i)))) then
                   nbseg=nbseg+1
                   szSA=szSA+1
                   ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine &
                        &construit_carcel (8)")
                   tmpSg = createSeg(pt2x,pt2y,cx+coef(i),cy+coef2(i),tmp,tmp)
                   tmpSg%sectg = i
                   tmpSg%sectd = mod(i+6,8)+1
                   tabSegArc(szSA) = tmpSg
                end if
             end if
          end do
       end if
    elseif (sectori == S_TXS_tot) then
       tgx = ssx * tg
       tgy = ssy * tg
       tmp = mix(size(mix))
       coef  = 0.5d0 *(/ ssx, tgy, -tgy, -ssx, -ssx, -tgy,  tgy,  ssx /)
       coef2 = 0.5d0 *(/ tgx, ssy, ssy,   tgx, -tgx, -ssy, -ssy, -tgx /) 
       if (sectorj==0) then
          do i = 1,8
             nbseg=nbseg+1 !
             szSA=szSA+1
             ! programmation defensive
             if (szSA > dimTabSegArc)&
                  call XABORT("G2S: memory problem in routine &
                  &construit_carcel (9)")
             tmpSg = createSeg(cx,cy,cx+coef(i),cy+coef2(i),tmp,tmp)
             tmpSg%sectg = i
             tmpSg%sectd = mod(i+6,8)+1
             tabSegArc(szSA) = tmpSg
          end do
       else 
          ar = createArc(ccx,ccy,radius(sectorj+1),0.d0,0.d0,tmp,tmp)
          do i = 1,8
             sg = createSeg(cx,cy,cx+coef(i),cy+coef2(i),tmp,tmp)
             nb = abs(interSgAr(sg,ar,pt1x,pt1y,pt2x,pt2y))
             if(nb==1) then
                if(.not.(isEqualConst(pt1x,cx+coef(i)).and. &
                     isEqualConst(pt1y,cy+coef2(i)))) then
                   nbseg=nbseg+1
                   szSA=szSA+1
                   ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine &
                        &construit_carcel (10)")
                   tmpSg = createSeg(pt1x,pt1y,cx+coef(i),cy+coef2(i),tmp,tmp)
                   tmpSg%sectg = i
                   tmpSg%sectd = mod(i+6,8)+1
                   tabSegArc(szSA) = tmpSg
                end if
             else if(nb==2) then
                if(.not.(isEqualConst(pt2x,cx+coef(i)).and. &
                     isEqualConst(pt2y,cy+coef2(i)))) then
                   nbseg=nbseg+1
                   szSA=szSA+1
                   ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine&
                        & construit_carcel (11)")
                   tmpSg = createSeg(pt2x,pt2y,cx+coef(i),cy+coef2(i),tmp,tmp)
                   tmpSg%sectg = i
                   tmpSg%sectd = mod(i+6,8)+1
                   tabSegArc(szSA) = tmpSg
                end if
             end if
          end do
      end if
    elseif (sectori == S_WM_tot) then
       tgx = ssx * tg
       tgy = ssy * tg
       tmp = mix(size(mix))
       coef  = 0.5d0 *(/ ssx, tgy, -tgy, -ssx, -ssx, -tgy,  tgy,  ssx /)
       coef2 = 0.5d0 *(/ tgx, ssy, ssy,   tgx, -tgx, -ssy, -ssy, -tgx /) 
       if (sectorj==0) then
          do i = 1,8
             nbseg=nbseg+1 !
             szSA=szSA+1
             ! programmation defensive
             if (szSA > dimTabSegArc)&
                  call XABORT("G2S: memory problem in routine &
                  &construit_carcel (12)")
             tmpSg = createSeg(cx,cy,cx+coef(i),cy+coef2(i),tmp,tmp)
             tmpSg%sectg = i
             tmpSg%sectd = mod(i+6,8)+1
             tabSegArc(szSA) = tmpSg
          end do
          do i = 1, 4
             nbseg=nbseg+1
             szSA=szSA+1
             ! programmation defensive
             if (szSA > dimTabSegArc)&
                  call XABORT("G2S: memory problem in routine &
                  &construit_carcel (13)")
             tmpSg = createSeg(cx+coef(2*i-1),cy+coef2(2*i-1),cx+coef(2*i),cy+coef2(2*i),tmp,tmp)
             tmpSg%sectg=i
             tmpSg%sectd= mod(i+6,8)+1
             tabSegArc(szSA) = tmpSg
          enddo
        else 
          ar = createArc(ccx,ccy,radius(sectorj+1),0.d0,0.d0,tmp,tmp)
          do i = 1,8
             sg = createSeg(cx,cy,cx+coef(i),cy+coef2(i),tmp,tmp)
             nb = abs(interSgAr(sg,ar,pt1x,pt1y,pt2x,pt2y))
             if(nb==1) then
                if(.not.(isEqualConst(pt1x,cx+coef(i)).and. &
                     isEqualConst(pt1y,cy+coef2(i)))) then
                   nbseg=nbseg+1
                   szSA=szSA+1
                   ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine &
                        &construit_carcel (14)")
                   tmpSg = createSeg(pt1x,pt1y,cx+coef(i),cy+coef2(i),tmp,tmp)
                   tmpSg%sectg = i
                   tmpSg%sectd = mod(i+6,8)+1
                   tabSegArc(szSA) = tmpSg
                end if
             else if(nb==2) then
                if(.not.(isEqualConst(pt2x,cx+coef(i)).and. &
                     isEqualConst(pt2y,cy+coef2(i)))) then
                   nbseg=nbseg+1
                   szSA=szSA+1
                   ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine &
                        &construit_carcel (15)")
                   tmpSg = createSeg(pt2x,pt2y,cx+coef(i),cy+coef2(i),tmp,tmp)
                   tmpSg%sectg = i
                   tmpSg%sectd = mod(i+6,8)+1
                   tabSegArc(szSA) = tmpSg
                end if
             end if
          end do
          ! creation des 4 segments fermant les ailes
          do i = 1, 4
             pt1x = coef(2*i-1)+cx
             pt1y = coef2(2*i-1)+cy
             pt2x = coef(2*i)+cx
             pt2y = coef2(2*i)+cy
             dist1 = (pt1x-ccx)*(pt1x-ccx)+(pt1y-ccy)*(pt1y-ccy)
             dist2 = (pt2x-ccx)*(pt2x-ccx)+(pt2y-ccy)*(pt2y-ccy)

             if (sqrt(dist1)>radius(sectorj+1).and.sqrt(dist2)>radius(sectorj)) then
                nbseg=nbseg+1
                szSA=szSA+1
                ! programmation defensive
                if (szSA > dimTabSegArc)&
                     call XABORT("G2S: memory problem in routine &
                     &construit_carcel (16)")
                tmpSg = createSeg(cx+coef(2*i-1),cy+coef2(2*i-1),cx+coef(2*i),cy+coef2(2*i),tmp,tmp)
                !tmpSg%sectg=i
                !tmpSg%sectd= mod(i+6,8)+1
                 tmpSg%sectg=mod(2*i+5,8)+1
                tmpSg%sectd=2*i
                tabSegArc(szSA) = tmpSg
             else
                call XABORT("G2S : Intersection troubles with WindMill")
             endif
          enddo
       end if
    else
       call XABORT("G2S : type of sectorisation not recognised")
    endif

    !prise en compte des milieux pour un segment completement a l'interieur
    !d'une zone annulaire (ie non coupe)
    ! CS-IV : La boucle s'applique sur les nbseg derniers segments crees
    ! CS-IV : ie sur les segments que l'on vient de faire (nbseg et szSA sont
    ! CS-IV : incrementes de la meme facon mais le premier part de 0 alors que 
    ! CS-IV : part du nombre de SA deja crees.
    ! CS-IV : donc le commentaire d'origine est faux, et devrait etre : 
    ! CS-IV : Prise en compte des milieux pour les segments crees.
     do i = 0,nbseg-1
       sg = tabSegArc(szSA-i)
       d = distance(ccx,ccy,sg%x,sg%y,sg%dx,sg%dy,tmpx,tmpy)
       tmpx=ccx-sg%x  ; tmpy=ccy-sg%y  ; dorig = sqrt((tmpx*tmpx)+(tmpy*tmpy))
       tmpx=ccx-sg%dx ; tmpy=ccy-sg%dy ; dextr = sqrt((tmpx*tmpx)+(tmpy*tmpy))
       do j = 1,size(radius)-1
          if((min(d,dorig,dextr)+epsilon>radius(j)).and. &
               & (max(d,dorig,dextr)-epsilon<radius(j+1))) then
             if(sg%mixg/=fooMix) tabSegArc(szSA-i)%mixg=mix(j)
             if(sg%mixd/=fooMix) tabSegArc(szSA-i)%mixd=mix(j)
          end if
       end do
    end do
    call construit_tube(ccx,ccy,radius,mix,cluster,szSA)
    call coupeCercle(keepSz,szSA,nbseg,tRec)
    
    if(sectori/=S_not) then
       call splitSegsForSect(keepSz,szSA)
       call majSectori(keepSz,szSA,sectori,tRec,cx,cy)
    end if
  end subroutine construit_carcel

  subroutine construit_hexhom(cx,cy,sd,mix,szSA,sectori)
    double precision,intent(in)          :: cx,cy,sd
    integer,intent(in)                   :: mix
    integer,intent(inout)                :: szSA
    integer,intent(in)                   :: sectori

    double precision,dimension(7)  :: xx,yy
    double precision               :: sqrt3_2S,S_2
    integer                        :: i
    integer,parameter,dimension(7) :: sect=(/1,2,3,4,5,6,1/)
    type(t_segArc)                 :: tmpSg

    ! programmation defensive
    integer :: dimtabSegArc
    dimtabSegArc = size(tabSegArc)

    S_2=5.d-1*sd ; sqrt3_2S=S_2*sqrt(3.d0)
    xx(1) = cx + S_2 ; yy(1) = cy + sqrt3_2S
    xx(2) = cx - S_2 ; yy(2) = yy(1)
    xx(3) = cx - sd  ; yy(3) = cy
    xx(4) = xx(2)    ; yy(4) = cy - sqrt3_2S
    xx(5) = xx(1)    ; yy(5) = yy(4)
    xx(6) = cx + sd  ; yy(6) = cy
    xx(7) = xx(1)    ; yy(7) = yy(1)
    do i = 1,6
       szSA=szSA+1
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine construit_hexhom (1)")
       tmpSg = createSeg(xx(i),yy(i),xx(i+1),yy(i+1),mix,fooMix)
       tabSegArc(szSA) = tmpSg
    end do
    if(sectori==S_X_tot) then
       do i = 1,6
          szSA=szSA+1
          ! programmation defensive
          if (szSA > dimTabSegArc)&
               call XABORT("G2S: memory problem in routine construit_hexhom (2)")
          tmpSg = createSeg(cx,cy,xx(i),yy(i),mix,mix)
          tmpSg%sectg = sect(i+1)
          tmpSg%sectd = sect(i)
          tabSegArc(szSA) = tmpSg
       end do
    else if(sectori/=S_not) then
       call XABORT("G2S : type of sectorisation not recognised for &
            &homogeneous hexagonal geometry")
    endif
  end subroutine construit_hexhom

  subroutine construit_hex(cx,cy,sd,mix,szSA,sectori,sectorj,radius,ocx,ocy)
    double precision,intent(in)          :: cx,cy,sd
    integer,intent(in)                   :: mix
    integer,intent(inout)                :: szSA
    integer,intent(in)                   :: sectori,sectorj
    double precision,dimension(:),intent(in) :: radius
    double precision,intent(in)          :: ocx,ocy

    double precision,dimension(7)  :: xx,yy
    double precision               :: sqrt3_2S,S_2
    double precision               :: pt1x,pt1y,pt2x,pt2y
    integer                        :: i,nb
    integer,parameter,dimension(7) :: sect=(/1,2,3,4,5,6,1/)
    type(t_segArc)                 :: sg,ar,tmpSg

    ! programmation defensive
    integer :: dimtabSegArc 
    dimtabSegArc = size(tabSegArc)

    S_2=5.d-1*sd ; sqrt3_2S=S_2*sqrt(3.d0)
    xx(1) = cx + S_2 ; yy(1) = cy + sqrt3_2S
    xx(2) = cx - S_2 ; yy(2) = yy(1)
    xx(3) = cx - sd  ; yy(3) = cy
    xx(4) = xx(2)    ; yy(4) = cy - sqrt3_2S
    xx(5) = xx(1)    ; yy(5) = yy(4)
    xx(6) = cx + sd  ; yy(6) = cy
    xx(7) = xx(1)    ; yy(7) = yy(1)
    do i = 1,6
       szSA=szSA+1
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine construit_hex (1)")
       tmpSg = createSeg(xx(i),yy(i),xx(i+1),yy(i+1),mix,fooMix)
       tabSegArc(szSA) = tmpSg
    end do
    if (sectori==S_X_tot) then
       if (sectorj==0) then
          do i = 1,6      
             szSA=szSA+1
             ! programmation defensive
             if (szSA > dimTabSegArc)&
                  call XABORT("G2S: memory problem in routine construit_hex (2)")
             tmpSg = createSeg(cx,cy,xx(i),yy(i),mix,mix)
             tmpSg%sectg = sect(i+1)
             tmpSg%sectd = sect(i)
             tabSegArc(szSA) = tmpSg
          end do
       else 
          ! creation du + grand cercle impermeable a la sectorisation
          ar = createArc(ocx,ocy,radius(sectorj+1),0.d0,0.d0,mix,mix)
          do i = 1,6
             sg=createSeg(cx,cy,xx(i),yy(i),mix,mix)
             nb=interSgAr(sg,ar,pt1x,pt1y,pt2x,pt2y)
             select case(nb)
             case(-1)
                if(.not.(isEqualConst(pt1x,xx(i)).and. &
                     isEqualConst(pt1y,yy(i)))) then
                   szSA=szSA+1
                   ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine &
                        &construit_hex (3)")
                   tmpSg = createSeg(pt1x,pt1y,xx(i),yy(i),mix,mix)
                   tmpSg%sectg = sect(i+1)
                   tmpSg%sectd = sect(i)
                   tabSegArc(szSA) = tmpSg
                end if
             case(0)
                if(.not.(isEqualConst(cx,xx(i)).and.isEqualConst(cy,yy(i)))) then
                   szSA=szSA+1
                   ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine &
                        &construit_hex (4)")
                   tmpSg = createSeg(cx,cy,xx(i),yy(i),mix,mix)
                   tmpSg%sectg = sect(i+1)
                   tmpSg%sectd = sect(i)
                   tabSegArc(szSA) = tmpSg           
                end if
             case(1)
                if(.not.(isEqualConst(cx,pt1x).and.isEqualConst(cy,pt1y))) then
                   szSA=szSA+1
                   ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine &
                        &construit_hex (5)")
                   tmpSg = createSeg(cx,cy,pt1x,pt1y,mix,mix)
                   tmpSg%sectg = sect(i+1)
                   tmpSg%sectd = sect(i)
                   tabSegArc(szSA) = tmpSg
                end if
             case(2)
                if(.not.(isEqualConst(cx,pt1x).and.isEqualConst(cy,pt1y))) then
                   szSA=szSA+1
                   ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine &
                        &construit_hex (6)")
                   tmpSg = createSeg(cx,cy,pt1x,pt1y,mix,mix)
                   tmpSg%sectg = sect(i+1)
                   tmpSg%sectd = sect(i)
                   tabSegArc(szSA) = tmpSg
                end if
                if(.not.(isEqualConst(pt2x,xx(i)).and. &
                     isEqualConst(pt2y,yy(i)))) then
                   szSA=szSA+1
                   ! programmation defensive
                   if (szSA > dimTabSegArc)&
                        call XABORT("G2S: memory problem in routine &
                        &construit_hex (7)")
                   tmpSg = createSeg(pt2x,pt2y,xx(i),yy(i),mix,mix)
                   tmpSg%sectg = sect(i+1)
                   tmpSg%sectd = sect(i)
                   tabSegArc(szSA) = tmpSg
                end if
             end select
          end do
       endif
    else if(sectori/=S_not) then
       call XABORT("G2S : type of sectorisation not recognised for &
            &hexagonal pincell geometry")
    endif
  end subroutine construit_hex

  subroutine construit_hexcel(cx,cy,sd,turn,radius,offcx,offcy, &
       & mix,secori,secorj,cluster,szSA)
    double precision,intent(in)           :: cx,cy,sd,offcx,offcy
    double precision,dimension(:),intent(in) :: radius
    integer,intent(in)                    :: turn,secori,secorj
    integer,dimension(:)                  :: mix
    type(t_cluster),dimension(:),pointer  :: cluster
    integer,intent(inout)                 :: szSA

    integer                       :: i,j,keepSz,nbseg
    double precision              :: tocx,tocy,tpi_3,costp,sintp
    double precision              :: d,tmpx,tmpy,dorig,dextr
    type(t_segArc)                :: sg

    keepSz=szSA
    if(turn<=6) then
       tpi_3=(turn-1)*pi_3_c ; costp=cos(tpi_3) ; sintp=sin(tpi_3)
       tocx = cx + costp*offcx - sintp*offcy
       tocy = cy + sintp*offcx + costp*offcy
    else
       tpi_3=(12-turn)*pi_3_c ; costp=cos(tpi_3) ; sintp=sin(tpi_3)
       tocx = cx + costp*offcx - sintp*offcy
       tocy = cy - sintp*offcx - costp*offcy
    end if
    call construit_hex(cx,cy,sd,mix(size(mix)),szSA,secori,secorj,radius,tocx,tocy)
    !prise en compte des milieux pour un segement completement a l'interieur
    !d'une zone annulaire (ie non coupee)
    do i = 0,szSA-keepSz-1
       sg = tabSegArc(szSA-i)
       d = distance(tocx,tocy,sg%x,sg%y,sg%dx,sg%dy,tmpx,tmpy)
       tmpx=tocx-sg%x  ; tmpy=tocy-sg%y  ; dorig=sqrt((tmpx*tmpx)+(tmpy*tmpy))
       tmpx=tocx-sg%dx ; tmpy=tocy-sg%dy ; dextr=sqrt((tmpx*tmpx)+(tmpy*tmpy))
       do j = 2,size(radius)-1
          if((min(d,dorig,dextr)>radius(j)).and. &
               & (max(d,dorig,dextr)<radius(j+1))) then
             if(sg%mixg/=fooMix) tabSegArc(szSA-i)%mixg=mix(j)
             if(sg%mixd/=fooMix) tabSegArc(szSA-i)%mixd=mix(j)
          end if
       end do
    end do
    nbseg=szSA-keepSz
    !creation des anneaux
    call construit_tube(tocx,tocy,radius,mix,cluster,szSA)
    !creation des anneaux
    call coupeCercle(keepSz,szSA,nbseg,tHex)
    !creation des anneaux
    if(secori/=S_not) call majSectori(keepSz,szSA,secori,tHex,cx,cy)
    !creation des anneaux
  end subroutine construit_hexcel

  subroutine construit_tri2d(cx,cy,sd,turn,split,mix,szSA)
    double precision,intent(in)  :: cx,cy,sd
    integer,intent(in)           :: turn,split
    integer,dimension(:),allocatable :: mix
    integer,intent(inout)        :: szSA

    integer          :: i,j,k,dimx,dimy,indMix,deltaIndMix
    type(t_segArc)   :: sg
    double precision :: deltax,deltay

    ! CS-IV : F_C_2
    !double precision,dimension(:),allocatable :: xx,yy
    double precision, dimension(2*split+1) :: xx
    double precision, dimension(split+1) :: yy

    ! programmation defensive
    integer :: dimtabSegArc 
    dimtabSegArc = size(tabSegArc)

    !on cree les triangles sans turn en (0,0), puis on les tourne
    !et on les translate

    !dimensionnement et remplissage des coordonnees
    ! x
    dimx = 2*split+1
    deltax = 5.d-1*sd/split
    ! CS-IV : F_C_2 allocate(xx(dimx))
    xx(1) = -5.d-1*sd
    do i = 2,dimx
       xx(i) = xx(i-1) + deltax
    end do
    ! y
    dimy = split+1
    ! CS-IV : F_C_2 allocate(yy(dimy))
    deltay = deltax*sqrt(3.d0)
    yy(1) = -25.d-2*sd*sqrt(3.d0)
    do j = 2,dimy
       yy(j) = yy(j-1) + deltay
    end do
    !creation des triangles
    !en trois etapes successives : 1=_ ; 2=\ ; 3=/

    !Premiere ligne
    ! Triangle de debut de ligne
    ! 1
    szSA=szSA+1
    ! programmation defensive
    if (szSA > dimTabSegArc)&
         call XABORT("G2S: memory problem in routine construit_tri2d (1)")
    sg = createSeg(xx(1),yy(1),xx(3),yy(1),mix(1),fooMix)
    tabSegArc(szSA) = translateAndTurnTriSg(cx,cy,turn,sg)
    if(split/=1) then
       ! 2
       szSA=szSA+1
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine construit_tri2d (2)")
       sg = createSeg(xx(3),yy(1),xx(2),yy(2),mix(1),mix(2))
       tabSegArc(szSA)=translateAndTurnTriSg(cx,cy,turn,sg)
    end if
    ! 3
    szSA=szSA+1
    ! programmation defensive
    if (szSA > dimTabSegArc)&
         call XABORT("G2S: memory problem in routine construit_tri2d (3)")
    sg = createSeg(xx(2),yy(2),xx(1),yy(1),mix(1),fooMix)
    tabSegArc(szSA)=translateAndTurnTriSg(cx,cy,turn,sg)
    ! Triangles de milieu de ligne
    do i = 2,split-1
       ! 1
       szSA=szSA+1
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine construit_tri2d (4)")
       sg = createSeg(xx(2*i-1),yy(1),xx(2*i+1),yy(1),mix(2*i-1),fooMix)
       tabSegArc(szSA)=translateAndTurnTriSg(cx,cy,turn,sg)
       ! 2
       szSA=szSA+1
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine construit_tri2d (5)")
       sg = createSeg(xx(2*i+1),yy(1),xx(2*i),yy(2),mix(2*i-1),mix(2*i))
       tabSegArc(szSA) = translateAndTurnTriSg(cx,cy,turn,sg)
       ! 3
       szSA=szSA+1
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine construit_tri2d (6)")
       sg = createSeg(xx(2*i),yy(2),xx(2*i-1),yy(1),mix(2*i-1),mix(2*i-2))
       tabSegArc(szSA) = translateAndTurnTriSg(cx,cy,turn,sg)
    end do
    ! Triangle de fin de ligne
    if(split/=1) then
       ! 1
       szSA=szSA+1
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine construit_tri2d (7)")
       sg = createSeg(xx(dimx-2),yy(1),xx(dimx),yy(1),mix(2*split-1),fooMix)
       tabSegArc(szSA)=translateAndTurnTriSg(cx,cy,turn,sg)
    end if
    ! 2
    szSA=szSA+1
    ! programmation defensive
    if (szSA > dimTabSegArc)&
         call XABORT("G2S: memory problem in routine construit_tri2d (8)")
    sg = createSeg(xx(dimx),yy(1),xx(dimx-1),yy(2),mix(2*split-1),fooMix)
    tabSegArc(szSA)=translateAndTurnTriSg(cx,cy,turn,sg)
    if(split/=1) then
       ! 3
       szSA=szSA+1
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine construit_tri2d (9)")
       sg = createSeg(xx(dimx-1),yy(2),xx(dimx-2),yy(1),mix(2*split-1), &
            & mix(2*split-2))
       tabSegArc(szSA)=translateAndTurnTriSg(cx,cy,turn,sg)
    end if

    !Autres lignes
    indMix = 2*split+1
    do j = 2,split
       k = j
       indMix = indMix - 1
       deltaIndMix = 2*(split-j+1)     
       do i = 1,dimy-j
          ! 1
          szSA=szSA+1
          ! programmation defensive
          if (szSA > dimTabSegArc)&
               call XABORT("G2S: memory problem in routine construit_tri2d (10)")
          sg = createSeg(xx(k),yy(j),xx(k+2),yy(j), &
               & mix(indMix),mix(indMix-deltaIndMix))
          tabSegArc(szSA)=translateAndTurnTriSg(cx,cy,turn,sg)
          ! 2
          szSA=szSA+1
          ! programmation defensive
          if (szSA > dimTabSegArc)&
               call XABORT("G2S: memory problem in routine construit_tri2d (11)")
          if(i==dimy-j) then
             sg = createSeg(xx(k+2),yy(j),xx(k+1),yy(j+1), &
                  & mix(indMix),fooMix)
          else
             sg = createSeg(xx(k+2),yy(j),xx(k+1),yy(j+1), &
                  & mix(indMix),mix(indMix+1))
          end if
          tabSegArc(szSA)=translateAndTurnTriSg(cx,cy,turn,sg)
          ! 3
          szSA=szSA+1
          ! programmation defensive
          if (szSA > dimTabSegArc)&
               call XABORT("G2S: memory problem in routine construit_tri2d (12)")
          if(i==1) then
             sg = createSeg(xx(k+1),yy(j+1),xx(k),yy(j), &
                  & mix(indMix),fooMix)
          else
             sg = createSeg(xx(k+1),yy(j+1),xx(k),yy(j), &
                  & mix(indMix),mix(indMix-1))
          end if
          tabSegArc(szSA)=translateAndTurnTriSg(cx,cy,turn,sg)
          k = k + 2
          indMix = indMix + 2
       end do
    end do
    ! CS-IV : F_C_2 deallocate(xx,yy)
  end subroutine construit_tri2d

  subroutine construit_tube(cx,cy,radius,mix,cluster,szSA)
    double precision,intent(in)           :: cx,cy
    double precision,dimension(:),intent(in) :: radius
    integer,dimension(:), intent(in)         :: mix
    type(t_cluster),dimension(:),pointer  :: cluster
    integer,intent(inout)                 :: szSA

    integer          :: i,j,k,sr,sc,nbp,keepSize,nbRing,nbClus,merge
    double precision :: ccx,ccy,alpha
    type(t_segArc)   :: ce
    integer,dimension(:),allocatable :: lastMix

    ! programmation defensive
    integer :: dimtabSegArc
    dimtabSegArc = size(tabSegArc)

    keepSize = szSa

    !anneaux
    sr = size(radius)
    do i = 2,sr
       szSa = szSa + 1
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine construit_tube (1)")
       ce = createArc(cx,cy,radius(i),0.d0,0.d0,mix(i-1),mix(i))
       tabSegArc(szSa) = ce
    end do

    nbRing = szSa - keepSize
    keepSize = szSa

    merge = size(mix) !servira a identifier les milieux dans les clusters
    !clusters
    if(.not. associated(cluster)) return !pas de cluster
    sc = size(cluster)
    if(sc==0) return !pas de cluster
    ! calcul du milieux exterieur de chaque cluster
    allocate(lastMix(sc))
    do i = 1,sc
       lastMix(i) = fooMix
       do j = 2,sr
          if(cluster(i)%radiusOfPin < radius(j)) then
             lastMix(i) = mix(j-1)
             exit
          end if
       end do
    end do
    ! creation des clusters par rayon croissant (pour que les cercles
    ! de rayon max rencontrent les anneaux en premier car la recopie
    ! inverse l'ordre)
    do i = 1,sc
       sr = size(cluster(i)%radius)
       nbp = cluster(i)%nbrPin
       do j = 1,nbp
          !          write(*,*) 'offset cluster ',j,' : ',cluster(i)%angleOfPin
          !          alpha = deg2rad*cluster(i)%angleOfPin + (j-1)*dpi_c/nbp
          alpha = cluster(i)%angleOfPin + (j-1)*dpi_c/nbp
          ccx = cx + cluster(i)%radiusOfPin * cos(alpha)
          ccy = cy + cluster(i)%radiusOfPin * sin(alpha)
          !tous les cercles sauf le dernier
          do k = 1,sr-2
             szSa = szSa + 1
             ! programmation defensive
             if (szSA > dimTabSegArc)&
                  call XABORT("G2S: memory problem in routine construit_tube (2)")
             ce = createArc(ccx,ccy,cluster(i)%radius(k+1),0.d0,0.d0, &
                  merge+k,merge+k+1)
             tabSegArc(szSa) = ce
          end do
          !cercle le plus grand
          szSa = szSa + 1
          ! programmation defensive
          if (szSA > dimTabSegArc)&
               call XABORT("G2S: memory problem in routine construit_tube (3)")
          ce = createArc(ccx,ccy,cluster(i)%radius(sr),0.d0,0.d0, &
               merge+sr-1,lastMix(i))
          tabSegArc(szSa) = ce        
       end do
       merge = merge + sr - 1
    end do
    deallocate(lastMix)

    nbClus = szSa - keepSize
    !coupage des anneaux par les clusters
    call cutClusters(nbRing,nbClus,szSa)
  end subroutine construit_tube

end module construire
