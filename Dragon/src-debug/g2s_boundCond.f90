!
!-----------------------------------------------------------------------
!
!Purpose:
! Process data relative to boundary conditions.
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
! Deux structures differentes ont ete definies en variables globales:
!  - bCData : pour stocker les donnees en entree du code
!  - SALbCData : pour stocker les donnees correspondantes a fournir dans
!                le jeux de donnees SAL genere
! Le fonction definies sont
!  - appliBoundariConditions : calcule les valeurs des champs des structures
!                              en appelant les fonctions specifiques a chaque
!                              type de geometrie
! \\\\
!  - initializebCData : mise a zero des stuctures
!  - destroybCData : liberation de la memoire
!  - appliBoundariConditionsForXXX : fonction specifique au type XXX avec
!                                    XXX = Rec, Hex, Tri ou Tub
!  - setBoundSide : applique une condition limite de contact aux elements
!  - setBoundCut : applique une condition limite avec reduction du domaine
!                  aux elements
!  - prepareSALBCData : prepare les donnees de conditions limites pour SAL
!
!-----------------------------------------------------------------------
!
module boundCond
  use cellulePlaced
  use constType
  use constUtiles
  use GANLIB
  use segArc

  implicit none  

  type t_bCData
     double precision,dimension(2) :: sidexy !longueur en xy de la gigogne
     !ext (pour rec et tri ST60 ou COMPLETE), ou cote d'un hexagone et cote
     !de l'assemblage (pour hexa et tri S30 et SA60)
     double precision,dimension(4) :: minmaxXY
     integer,dimension(6)          :: bc
     double precision,dimension(6) :: albedo
     integer,dimension(6)          :: albInd
     double precision,dimension(2) :: toOrig_xy !vecteur de translation du
     !centre vers le coins inferieur gauche
     integer                       :: iHex !type de geo si hexagone
     integer                       :: iTri !type de geo si triangle
  end type t_bCData

  type(t_bCData),save  :: bCData

  type t_SALbCData
     integer                      :: sunsetType
     integer                      :: SALtype
     integer                      :: nber
     integer,dimension(:),allocatable :: elemNb
     real                         :: albedo
     real                         :: tx,ty
     real                         :: cx,cy,angle
  end type t_SALbCData

  type(t_SALbCData),dimension(4),save  :: SALbCDataTab

contains

  subroutine initializebCData()
    SALbCDataTab(1:4)%nber = 0
  end subroutine initializebCData

  subroutine destroybCData()
    integer :: i
    do i = 1,4
       if (allocated(SALbCDataTab(i)%elemNb)) then
          deallocate(SALbCDataTab(i)%elemNb)
       end if
    end do
  end subroutine destroybCData

  subroutine appliBoundariConditions(ip,szSA,nbCLP)
    type(c_ptr),intent(in):: ip
    integer,intent(inout) :: szSA
    integer,intent(out)   :: nbCLP

    select case(geomTyp)
    case(RecTyp)
       !       write(*,*) 'entering appliBoundariConditionsForRec'
       call appliBoundariConditionsForRec(ip,szSA,nbCLP)
    case(HexTyp)
       !       write(*,*) 'entering appliBoundariConditionsForHex'
       call appliBoundariConditionsForHex(ip,szSA,nbCLP)
    case(TriaTyp)
       !       write(*,*) 'entering appliBoundariConditionsForTri'
       call appliBoundariConditionsForTri(ip,szSA,nbCLP)
    case(TubeTyp)
       !       write(*,*) 'entering appliBoundariConditionsForTub'
       call appliBoundariConditionsForTub(ip,nbCLP)
    end select
  end subroutine appliBoundariConditions

  subroutine appliBoundariConditionsForRec(geoIp,szSA,nbCLP)
    type(c_ptr),intent(in)    :: geoIp
    integer,intent(inout) :: szSA
    integer,intent(out)   :: nbCLP

    double precision              :: rminx,rminy,rmaxx,rmaxy
    double precision,dimension(4) :: x,y,xx,yy,cx,cy,cxx,cyy
    type(c_ptr)                   :: ip
    integer                       :: i,nbCut
    real,dimension(2)             :: tmpTab2
    real,dimension(4)             :: tmpTab4
    real,dimension(6)             :: tmpTab6
    type(t_segArc)                :: sg

    ! programmation defensive
    integer :: dimTabSegArc

    dimTabSegArc = size(tabSegArc)

    ip = geoIp
    !recuperations des donnees sur les conditions aux limites
    call LCMGET(ip,'NCODE       ',bCData%bc)
    call LCMGET(ip,'ZCODE       ',tmpTab6) ; bCData%albedo=tmpTab6
    call LCMGET(ip,'ICODE       ',bCData%albInd)
    !    write(*,*) 'NCODE :',bCData%bc
    !    write(*,*) 'ZCODE :',bCData%albedo
    !    write(*,*) 'ICODE :',bCData%albInd
    call LCMSIX(ip,'NEW-DATA    ',1)
    call LCMSIX(ip,'BOUND-DATA  ',1)
    call LCMGET(ip,'SIDEXY      ',tmpTab2) ; bCData%sidexy=tmpTab2
    call LCMGET(ip,'MINMAXXY    ',tmpTab4) ; bCData%minmaxXY=tmpTab4
    call LCMSIX(ip,'BOUND-DATA  ',2)
    call LCMSIX(ip,'NEW-DATA    ',2)

    !exploitation des donnees
    rmaxx = 0.5d0*bCData%sidexy(1) ; rminx = -rmaxx
    rmaxy = 0.5d0*bCData%sidexy(2) ; rminy = -rmaxy
    x(1) = rminx ; y(1) = rmaxy ; xx(1) = rminx ; yy(1) = rminy
    x(2) = rmaxx ; y(2) = rminy ; xx(2) = rmaxx ; yy(2) = rmaxy
    x(3) = rminx ; y(3) = rminy ; xx(3) = rmaxx ; yy(3) = rminy
    x(4) = rmaxx ; y(4) = rmaxy ; xx(4) = rminx ; yy(4) = rmaxy 
    cx(1)  = bCData%minmaxXY(1) ; cy(1)  = rmaxy
    cxx(1) = bCData%minmaxXY(1) ; cyy(1) = rminy
    cx(2)  = bCData%minmaxXY(3) ; cy(2)  = rminy
    cxx(2) = bCData%minmaxXY(3) ; cyy(2) = rmaxy
    cx(3)  = rminx              ; cy(3)  = bCData%minmaxXY(2)
    cxx(3) = rmaxx              ; cyy(3) = bCData%minmaxXY(2)
    cx(4)  = rmaxx              ; cy(4)  = bCData%minmaxXY(4)
    cxx(4) = rminx              ; cyy(4) = bCData%minmaxXY(4)

    !creation des conditions aux limites (pour une geometrie carre)
    !lorsqu'on coupe, on garde le cote gauche (sens trigo)
    nbCLP=0
    nbCut=0
    do i = 1,4
       select case(bCData%bc(i))
       case(B_Void,B_Refl,B_Ssym,B_Albe,B_Zero,B_Tran,B_Pi_2,B_Pi)
          !write(*,*) 'reflective boundary condition ',i
          if (bCData%bc(i)/=B_Void) nbCLP=nbCLP+1
          call setBoundSide(x(i),y(i),xx(i),yy(i),bCData%bc(i)+100*i,szSA)
       case(B_Diag)
          !write(*,*) 'diagonal symmetry ',i,"seg ",szsa+1
          if (i>=3) cycle
          nbCut=nbCut+1 ; nbCLP=nbCLP+1
          szSA=szSA+1 
          ! programmation defensive
          if (szSA > dimTabSegArc)&
               call XABORT("G2S: memory problem in routine &
               &appliBoundaryConditionsForRec (1)")
          if (i==1) then
             sg=createSeg(x(4),y(4),x(3),y(3),fooMix,-(bCData%bc(i)+100*i))
          else
             sg=createSeg(x(3),y(3),x(4),y(4),fooMix,-(bCData%bc(i)+100*i))
          end if
          tabSegArc(szSA)=sg
       case(B_Syme)
          !write(*,*) 'syme ', i,"seg ",szsa+1
          nbCut=nbCut+1 ; nbCLP=nbCLP+1
          szSA=szSA+1  
          ! programmation defensive
          if (szSA > dimTabSegArc)&
               call XABORT("G2S: memory problem in routine &
               &appliBoundaryConditionsForRec (1)")
          sg=createSeg(cx(i),cy(i),cxx(i),cyy(i),fooMix,-(bCData%bc(i)+100*i))
          tabSegArc(szSA)=sg
       end select
    end do
    if (nbCut/=0) then
       !write(*,*) 'entering setBoundCut :',nbCut,szSA
       call setBoundCut_V2(nbCut,szSA)
       !write(*,*) 'leaving setBoundCut :',nbCut,szSA
    end if
  end subroutine appliBoundariConditionsForRec

  subroutine setBoundSide(x,y,xx,yy,nbCL,szSA)
    double precision,intent(in) :: x,y,xx,yy
    integer,intent(in)          :: nbCL,szSA

    type(t_segArc) :: sa
    integer        :: i

    do i = 1,szSA
       sa = tabSegArc(i)
       if (sa%typ/=tseg) cycle
       if (.not. ( estColi(xx-x,yy-y,sa%dx-sa%x,sa%dy-sa%y) .and. &
            !               &   estAligne(x,y,xx,yy,sa%x,sa%y) ) ) cycle
            &   pointsAlignes(x,y,xx,yy,sa%x,sa%y) ) ) cycle
       if (sa%mixd==fooMix) then
          tabSegArc(i)%mixd=-nbCL
       else if (sa%mixg==fooMix) then
          tabSegArc(i)%mixg=-nbCL
       end if
    end do
  end subroutine setBoundSide

  subroutine setBoundCut_V2(nbCut,szSA)
    integer, intent(in)    :: nbCut
    integer, intent(inout) :: szSA

    ! segments management
    integer :: iSA, iSC, iStr, thisSA, thisOtherSA, szSC
    type(t_SegArc) :: SA, SA1, SA2, SA3, SC, SC1, SC2, SC3
    double precision :: interangle, interAngle2, RefAngle, RefOtherAngle
    double precision :: thisAngle, thisOtherAngle, angle1,angle2
    ! follow up for intersections
    logical :: L_Inter, test1, test2, test3, test4       
    type(t_point) :: P1, P2, P3, P4
    integer :: n1, n2
    double precision :: aox,aoy,afx,afy
    ! sort
    ! find minimal angle
    double precision :: SCAngle, SA_x,SA_y,SA_dx,SA_dy

    ! Part 1 :  Backup of cutting straights & updating
    allocate(tabStrCut(nbCut), tabSegCut(szSA),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: setBoundCut_V2 => allocation pb")
    tabStrCut(1:nbCut) = tabSegArc(szSA-nbCut+1:szSA+1)
    szSA = szSA - nbCut
    !
    ! Part2 : Adaptation for cutting Straights 
    B_A1: do iStr = 1, nbCut
       SC = tabStrCut(iStr)
       B_B1: do iSA = iStr+1, nbCut
          SA = tabStrCut(iSA)
          L_Inter = InterSgSg_V2(SC,SA,n1,n2,P1,P2)
          if (L_Inter.and.n2==0) then
             SA1 = copySegExceptEnd(P1%x,P1%y,SA)
             if (estAGauche(SC%x,SC%y,SC%dx,SC%dy,SA1%x,SA1%y)) then
                tabStrCut(iSA) = SA1
             else
                SA1 = copySegExceptOrigin(P1%x,P1%y,SA)
                tabStrCut(iSA) = SA1
             endif
          end if
       enddo B_B1
    enddo B_A1

    ! Part 3 : process for the current straight
    B_A2: do iStr = 1, nbCut  
       szSC = 1 
       tabSegCut(szSC) = tabStrCut(iStr)
       iSC = szSC 

       B_B2: do
          if (iSC > szSC)  exit B_B2
          SC = tabSegCut(iSC)
          iSA = 1
          B_C2: do
             if (iSA > szSA) exit B_C2
             SA = tabSegArc(iSA)
             if (SA%typ == tseg) then  
                ! SA is a segment
                L_Inter = InterSgSg_V2(SC,SA,n1,n2,P1,P2)
                if (.not.L_Inter) then
                   ! no intersection
                   if  (estAGauche(SC%x,SC%y,SC%dx,SC%dy,SA%x,SA%y).or.&
                        estAGauche(SC%x,SC%y,SC%dx,SC%dy,SA%dx,SA%dy)) then
                      iSA = iSA + 1
                   else
                      call SlideSA(iSA,szSA)
                   endif
                   cycle B_C2 
                elseif (L_Inter.and.(n2==0)) then
                   ! only one intersection point
                   if (n1 == 22) then
                      SC1 = createSeg(SC%x,SC%y,P1%x,P1%y,SC%mixg,SC%mixd)
                      SC2 = createSeg(P1%x,P1%y,SC%dx,SC%dy,SC%mixg,SC%mixd)
                      SC = ReplaceSC(SC1,iSC)
                      call AddSC(SC2,szSC)
                      SA1 = copySegExceptEnd(P1%x,P1%y,SA)
                      if (.not.estAGauche(SC%x,SC%y,SC%dx,SC%dy,SA1%x,SA1%y)) &
                         SA1 = copySegExceptOrigin(P1%x,P1%y,SA)
                      SA = ReplaceSA(SA1,iSA) ; iSA = iSA + 1

                   elseif ((n1 == 21).or.(n1 == 23)) then
                      ! the intersection point 
                      SA1 = copySegExceptEnd(P1%x,P1%y,SA)
                      if (.not.estAGauche(SC%x,SC%y,SC%dx,SC%dy,SA1%x,SA1%y))&
                           SA1 = copySegExceptOrigin(P1%x,P1%y,SA)
                      SA = ReplaceSA(SA1,iSA) ; iSA = iSA + 1

                   elseif ((n1==11).or.(n1==13).or.(n1==31).or.(n1==33)) then
                      if ( ((n1==11).and.&
                           estAGAuche(SC%x,SC%y,SC%dx,SC%dy,SA%dx,SA%dy)) .or. &
                           ((n1==13).and.&
                           estAGauche(SC%x,SC%y,SC%dx,SC%dy,SA%dx,SA%dy)) .or. &
                           ((n1==31).and.&
                           estAGauche(SC%x,SC%y,SC%dx,SC%dy,SA%x,SA%y  )) .or. &
                           ((n1==33).and.&
                           estAGauche(SC%x,SC%y,SC%dx,SC%dy,SA%x, SA%y ))) then
                         iSA = iSA + 1
                      else
                         call SlideSA(iSA,szSA)
                      endif

                   elseif ((n1 == 12).or.(n1 == 32)) then
                      SC1 = createSeg(SC%x,SC%y,P1%x,P1%y,SC%mixg,SC%mixd)
                      SC2 = createSeg(P1%x,P1%y,SC%dx,SC%dy,SC%mixg,SC%mixd)
                      SC = ReplaceSC(SC1,iSC)
                      call AddSC(SC2,szSC)
                      if (estAGauche(SC%x,SC%y,SC%dx,SC%dy,SA%x,SA%y)) then
                         iSA = iSA + 1
                      else
                         call SlideSA(iSA,szSA)
                      endif
                   else
                      call XABORT("SetBoundCut_V2 : inconsistency")
                   endif
                   cycle B_C2

                else  ! (L_Inter.and.n2 /= 0 )
                  call OverlaidSegmentsManagement(n1,n2,iSC,szSC,SC,iSA,szSA,SA)
                  cycle B_C2    
                endif
             else      
                ! Sa is an arc or a circle
                L_Inter = interSgAr_V2(SC,SA,n1,n2,P1,P2)
                if (SA%typ == tarc) then
                   call giveOrigine(SA,Aox,Aoy) ; call giveExtremite(SA,Afx,Afy)
                endif
                if (.not.L_Inter) then
                   if (SA%typ /= tcer) then
                      if (estAGauche(SC%x,SC%y,SC%dx,SC%dy,aox,aoy).or. &
                           estAGauche(SC%x,SC%y,SC%dx,SC%dy,afx,afy)) then
                         iSA = iSA + 1
                      else
                         call SlideSA(iSA,szSA)
                      endif
                   else 
                      ! SA is a circle
                      call giveFourPointsOnCircle(SA,P1,P2,P3,P4)
                      test1 = estAGauche(SC%x,SC%y,SC%dx,SC%dy,P1%x,P1%y)
                      test2 = estAGauche(SC%x,SC%y,SC%dx,SC%dy,P2%x,P2%y)
                      test3 = estAGauche(SC%x,SC%y,SC%dx,SC%dy,P3%x,P3%y)
                      test4 = estAGauche(SC%x,SC%y,SC%dx,SC%dy,P4%x,P4%y)
                      if  (test1.or.test2.or.test3.or.test4) then
                         ! circle will be cut by another SC
                         iSA = iSA + 1
                      else
                         call SlideSA(iSA,szSA)
                      endif
                   endif
                   cycle B_C2
                elseif (n2 == 0) then
                   if (n1>100) then
                      ! one intersection point which is a tangency point
                      if (estAGauche(SC%x,SC%y,SC%dx,SC%dy,aox,aoy)) then
                         SC1 = createSeg(SC%x,SC%y,P1%x,P1%y,SC%mixg,SC%mixd)
                         SC2 = createSeg(P1%x,P1%y,SC%x,SC%y,SC%mixg,SC%mixd)
                         SC = ReplaceSC(SC1,iSC)
                         call AddSC(SC2,szSC)
                         interAngle = calculeAngle(SA%x,SA%y,P1%x,P1%y)
                         if (SA%typ == tarc) then
                            SA1 = copyArcExceptEnd(interAngle,SA)
                            SA2 = copyArcExceptOrigin(interAngle,SA)
                            SA = ReplaceSA(SA1,iSA)
                            call AddSA(SA2,szSA)
                         else
                            interangle2 = anglenormal(interangle+pi_c)
                            SA1 = createArcFromCircle(interAngle,interAngle2,SA)
                            SA2 = createArcFromCircle(interAngle2,interAngle,SA)
                            SA = ReplaceSA(SA1,iSA)
                            call AddSA(SA2,szSA)
                         end if
                         iSA = iSA + 1
                      else
                         call SlideSA(iSA,szSA)
                      end if
                   elseif (n1<100) then
                      ! one intersection point P1
                      if (mod(n1,10) == 2) then
                         ! P1 is inside segment
                         SC1 = createSeg(SC%x,SC%y,P1%x,P1%y,SC%mixg,SC%mixd)
                         SC2 = createSeg(P1%x,P1%y,SC%dx,SC%dy,SC%mixg,SC%mixd)
                         SC = ReplaceSC(SC1,iSC)
                         call AddSC(SC2,szSC)
                         interAngle = calculeAngle(SA%x,SA%y,P1%x,P1%y)
                         if (n1/10 == 2) then
                            if (SA%typ == tarc) then
                               SA1 = copyArcExceptEnd(interAngle,SA)
                               if (.not.estAGauche &
                                    (SC%x,SC%y,SC%dx,SC%dy,aox,aoy)) &
                                   SA1 = copyArcExceptOrigin(interAngle,SA)
                               SA = ReplaceSA(SA1,iSA)
                            else
                               interangle2 = anglenormal(interangle+pi_c)
                               SA1 = createArcFromCircle &
                                    (interAngle,interAngle2,SA) 
                               call MediumPointOnArc(SA1,P1)
                               if (.not.estAGauche &
                                    (SC%x,SC%y,SC%dx,SC%dy,P1%x,P1%y)) &
                                    SA1 = createArcFromCircle &
                                          (interAngle2,interAngle,SA) 
                               SA = ReplaceSA(SA1,iSA)
                            endif
                         endif
                      else
                         ! P1 is the begin or the end of SC
                         if (n1/10 == 2) then
                            ! P1 is inside SA
                            interAngle = calculeAngle(SA%x,SA%y,P1%x,P1%y)
                            if (SA%typ == tarc) then
                               SA1 = copyArcExceptEnd(interAngle,SA)
                               if (.not.estAGauche &
                                        (SC%x,SC%y,SC%dx,SC%dy,aox,aoy)) &
                                  SA1 = copyArcExceptOrigin(interAngle,SA)
                            else
                               interangle2 = anglenormal(interangle+pi_c)
                               SA1 = createArcFromCircle &
                                     (interAngle,interAngle2,SA) 
                               call MediumPointOnArc(SA1,P1)
                               if (.not.estAGauche &
                                        (SC%x,SC%y,SC%dx,SC%dy,P1%x,P1%y)) &
                                  SA1 = createArcFromCircle &
                                        (interAngle2,interAngle,SA) 
                            endif
                            SA = ReplaceSA(SA1,iSA)
                            SC = ReplaceSC(SC1,iSC)
                            iSA = iSA + 1
                         else
                            ! P1 is at the begin or the end of SA
                            if (((n1/10==1).and.                               &
                                 estAGauche(SC%x,SC%y,SC%dx,SC%dy,afx,afy))    &
                                 .or.                                          &
                                 ((n1/10==3).and.                              &
                                 estAGauche(SC%x,SC%y,SC%dx,SC%dy,aox,aoy))) then
                               iSA =iSA + 1
                            else
                               call SlideSA(iSA,szSA)
                            endif
                         endif
                      endif
                   endif
                else 
                   ! two points of intersection
                   if ((n1==21.and.n2==23) .or.(n1==23.and.n2==21)) then
                      ! the intersection points are the extremities of SA
                      interangle  = calculeAngle(SA%x,SA%y,P1%x,P1%y)
                      interangle2 = calculeAngle(SA%x,SA%y,P2%x,P2%y)
                      angle1 = min(interangle,interangle2)
                      angle2 = max(interangle,interangle2)
                      if (SA%a > SA%b) then
                         SA1 = CopyArcExceptEnd(angle1,SA)
                         SA2 = CopyArcWithNewAngles(angle1,angle2,SA)
                         SA3 = CopyArcExceptOrigin(angle2,SA)
                      else
                         SA1 = CopyArcExceptEnd(angle2,SA)
                         SA2 = CopyArcWithNewAngles(angle2,angle1,SA)
                         SA3 = CopyArcExceptOrigin(angle1,SA)
                      endif
                      call MediumPointOnArc(SA2,P1)
                      if (estAGauche(SC%x,SC%y,SC%dx,SC%dy,P1%x,P1%y)) then
                         SA = ReplaceSA(SA2,iSA)
                      else
                         SA = ReplaceSA(SA1,iSA)
                         call AddSA(SA3,szSA)
                      end if
                      iSA = iSA + 1
                   elseif((n1==12.and.n2==32).or.(n1==32.and.n2==12)) then
                      SC1 = createSeg(SC%x,SC%y,P1%x,P1%y,SC%mixg,SC%mixd)
                      SC2 = createSeg(P1%x,P1%y,P2%x,P2%y,SC%mixg,SC%mixd)
                      SC3 = createSeg(P2%x,P2%y,SC%dx,SC%dy,SC%mixg,SC%mixd)
                      SC = ReplaceSC(SC1,iSC)
                      call Add2SC(SC2,SC3,szSC)
                      iSA = iSA + 1
                   elseif ((n1==11.and.n2==33).or. &
                           (n1==31.and.n2==13).or. &
                           (n1==13.and.n2==31)) then
                      ! P1 and P2 are at the begin or the end of SA and SC
                      ! nothing to do
                      iSA = iSA + 1
                   else
                      SC1 = createSeg(SC%x,SC%y,P1%x,P1%y,SC%mixg,SC%mixd)
                      SC2 = createSeg(P1%x,P1%y,P2%x,P2%y,SC%mixg,SC%mixd)
                      SC3 = createSeg(P2%x,P2%y,SC%dx,SC%dy,SC%mixg,SC%mixd)
                      if (SA%typ == tarc) then
                         interAngle = calculeAngle(SA%x,SA%y,P1%x,P1%y)
                         SA1 = copyArcExceptEnd(interAngle,SA)
                         interAngle = calculeAngle(SA%x,SA%y,P2%x,P2%y)
                         SA2 = copyArcWithNewAngles(SA1%b,interAngle,SA)
                         SA3 = copyArcExceptOrigin(SA2%a,SA)
                         ! if SA1 is on left side, SA3 is on left side too
                         ! if SA1 isn't on left side, SA2 is on left side
                         call giveOrigine(SA1,aox,aoy)
                         if (estAGauche(SC%x,SC%y,SC%dx,SC%dy,aox,aoy)) then
                            SA = ReplaceSA(SA1,iSA)
                            call AddSA(SA3,szSA)
                         else
                            SA = ReplaceSA(SA2,iSA)
                         endif
                      else
                         interAngle = calculeAngle(SA%x,SA%y,P1%x,P1%y)
                         interAngle2 = calculeAngle(SA%x,SA%y,P2%x,P2%y)
                         SA1 = createArcFromCircle(interAngle2,interAngle,SA)
                         call MediumPointOnArc(SA1,P3)
                         if (.not.estAGauche(SC%x,SC%y,SC%dx,SC%dy,P3%x,P3%y)) &
                              SA1 = createArcFromCircle &
                                    (interAngle,interAngle2,SA)
                         SA = ReplaceSA(SA1,iSA)
                      endif
                      SC = ReplaceSC(SC1,iSC)
                      call Add2SC(SC2,SC3,szSC)
                       iSA = iSA + 1
                   endif
                endif
                cycle B_C2
             endif
          enddo B_C2
          iSC = iSC + 1
       enddo B_B2
       
!!$       ! Part 3 : Sort
!!$       B_B3:do iSC = 1,szSC-1
!!$          P1%x = tabSegCut(iSC)%dx ; P1%y = tabSegCut(iSC)%dy
!!$          B_C3:do jSC = iSC+1,szSC
!!$             if (tabSegCut(jSC)%x == P1%x .and.tabSegCut(jSC)%y == P1%y) then
!!$                if(jSC==iSC+1)  then ; cycle B_B3
!!$                else
!!$                   SCBuffer = tabSegCut(iSC+1)
!!$                   tabSegCut(iSC+1) = tabSegCut(jSC)
!!$                   tabSegCut(jSC) = SCBuffer
!!$                endif
!!$             end if
!!$          end do B_C3
!!$       end do B_B3
       
       ! Part 4 : search for neighbourhood informations
       B_B4: do iSC = 1, szSC
          SC = tabSegCut(iSC)
          thisSA  = 0         ; thisOtherSA = 0
          RefAngle = infinity ; RefOtherAngle = infinity
          SCAngle = calculeAngle(SC%x,SC%y,SC%dx,SC%dy)
          ! Mix awarding
          B_C4: do iSA = 1, szSA
             SA = tabSegArc(iSA)
             if (SA%typ == tseg) then
                if (IsEqualConst(SC%x,SA%x).and.IsEqualConst(SC%y,SA%y)) then
                   thisAngle = calculeAngle(SA%x,SA%y,SA%dx,SA%dy) - SCAngle
                   if (thisAngle < RefAngle) then
                      thisSA = iSA ; RefAngle = thisAngle
                   endif
                elseif (IsEqualConst(SC%x,SA%dx).and.IsEqualConst(SC%y,SA%dy)) then
                   thisAngle = calculeAngle(SA%x,SA%y,SA%dx,SA%dy) + SCAngle
                   if (thisAngle < RefAngle) then
                      thisSA = iSA ; RefAngle = thisAngle
                   endif
                endif
             else if (SA%typ == tarc) then
                call giveOrigine(SA,SA_x,SA_y)
                call giveExtremite(SA,SA_dx,SA_dy)
                if (IsEqualConst(SC%x,SA_x) .and. IsEqualConst(SC%y,SA_y)) then
                   thisAngle = SA%a + pi_2_c - SCAngle
                   if (thisAngle < RefAngle) then
                      thisSA = iSA ; RefAngle = thisAngle
                   endif
                elseif (IsEqualConst(SC%x,SA_dx).and.IsEqualConst(SC%y,SA_dy)) then
                   thisAngle = SA%b - pi_2_c + SCAngle
                   if (thisAngle < RefAngle) then
                      thisSA = iSA ; RefAngle = thisAngle
                   endif
                endif
             endif
          enddo B_C4
          if (thisSA == 0) then
            call XABORT('unable to find element thisSA')
          endif
          if (tabSegArc(thisSA)%typ == tseg) then
             if (isEqualConst(tabSegArc(ThisSA)%x,SC%x).and.           &
                  isEqualConst(tabSegArc(ThisSA)%y,SC%y)) then
                SC%mixg      = tabSegArc(thisSA)%mixd
                SC%nodeg     = tabSegArc(thisSA)%noded
                SC%IndCellPg = tabSegArc(thisSA)%indCellPd 
             elseif(isEqualConst(tabSegArc(ThisSA)%dx,SC%x).and.       &
                  isEqualConst(tabSegArc(ThisSA)%dy,SC%y)) then
                SC%mixg = tabSegArc(thisSA)%mixg
                SC%nodeg     = tabSegArc(thisSA)%nodeg
                SC%IndCellPg = tabSegArc(thisSA)%indCellPg
              else 
                call XABORT("SetBoundCut_V2 : mix error (1)")
             endif
          else
             call giveOrigine(tabSegArc(thisSA),SA_x,SA_y)
             call giveExtremite(tabSegArc(thisSA),SA_dx,SA_dy)
             if (isEqualConst(SA_x,SC%x).and.isEqualConst(SA_y,SC%y)) then
                SC%mixg = tabSegArc(thisSA)%mixd
                SC%nodeg     = tabSegArc(thisSA)%noded
                SC%IndCellPg = tabSegArc(thisSA)%indCellPd 
              elseif (isEqualConst(SA_dx,SC%x).and.isEqualConst(SA_dy,SC%y)) then
                SC%mixg = tabSegArc(thisSA)%mixg
                SC%nodeg     = tabSegArc(thisSA)%nodeg
                SC%IndCellPg = tabSegArc(thisSA)%indCellPg
              else
                call XABORT("SetBoundCut_V2 : mix error (2)")
             endif
          endif

          ! Mix control
          SCAngle = calculeAngle(SC%dx,SC%dy,SC%x,SC%y)
          B_C5: do iSA = 1,szSA
             SA = tabSegArc(iSA)
             if (SA%typ == tseg) then
                if (IsEqualConst(SC%dx,SA%dx).and.IsEqualConst(SC%dy,SA%dy)) then
                   thisOtherAngle = calculeAngle(SA%x,SA%y,SA%dx,SA%dy) + SCAngle
                   if (thisOtherAngle < RefOtherAngle) then
                      thisOtherSA = iSA ; RefOtherAngle = thisOtherAngle
                   endif
                elseif (IsEqualConst(SC%dx,SA%x).and.IsEqualConst(SC%dy,SA%y)) then
                   thisOtherAngle = calculeAngle(SA%x,SA%y,SA%dx,SA%dy) - SCAngle
                   if (thisOtherAngle < RefOtherAngle) then
                      thisOtherSA = iSA ; RefOtherAngle = thisOtherAngle
                   endif
                endif
             elseif (SA%typ == tarc) then
                call giveOrigine(SA,SA_x,SA_y)
                call giveExtremite(SA,SA_dx,SA_dy)
                if (IsEqualConst(SC%dx,SA_dx).and.IsEqualConst(SC%dy,SA_dy)) then
                   thisOtherAngle = SA%B + pi_2_c - SCAngle
                   if (thisOtherAngle < RefOtherAngle) then
                      thisOtherSA = iSA ; RefOtherAngle = thisOtherAngle
                   endif
                elseif (IsEqualConst(SC%dx,SA_x).and.IsEqualConst(SC%dy,SA_y)) then
                   thisOtherAngle = SA%a - pi_2_c + SCAngle
                   if (thisOtherAngle < RefOtherAngle) then
                      thisOtherSA = iSA ; RefOtherAngle = thisOtherAngle
                   endif
                endif
             endif
          enddo B_C5
          if (thisOtherSA == 0) then
             ! we do not find a neigbour for SC
             call XABORT("SetBoundCut_V2 : incredible situation too")
          endif
          if (tabSegArc(thisOtherSA)%typ == tseg) then
             if (IsEqualConst(SC%dx,tabSegArc(thisOtherSA)%x).and. &
                  IsEqualConst(SC%dy,tabSegArc(thisOtherSA)%y)) then
                if (SC%mixg /= tabSegArc(thisOtherSA)%mixg) then
                   call XABORT("SetBoundCut_V2 : mix error (3)")
                endif
             elseif (IsEqualConst(SC%dx,tabSegArc(thisOtherSA)%dx).and. &
                  IsEqualConst(SC%dy,tabSegArc(thisOtherSA)%dy)) then
                if (SC%mixg /= tabSegArc(thisOtherSA)%mixd)  then
                   call XABORT("SetBoundCut_V2 : mix error (4)")
                end if
             else
                call XABORT("SetBoundCut_V2 : mix error (5)")
             endif
          else
             call giveOrigine(tabSegArc(thisOtherSA),SA_x,SA_y)
             call giveExtremite(tabSegArc(thisOtherSA),SA_dx,SA_dy)
             if (isEqualConst(SC%dx,SA_x).and.isEqualConst(SC%dy,SA_y)) then
                if (SC%mixg /= tabSegArc(thisOtherSA)%mixg) then
                   call XABORT("SetBoundCut_V2 : mix error (6)")
                end if
             elseif(isEqualCOnst(SC%dx,SA_dx).and. &
                    isEqualConst(SC%dy,SA_dy)) then
                if (SC%mixg /= tabSegArc(thisOtherSA)%mixd) then
                    call XABORT("SetBoundCut_V2 : mix error (7)")
                endif
             else
                call XABORT("SetBoundCut_V2 : mix error (8)")
             endif
          endif
          tabSegCut(iSC) = SC
       enddo B_B4
       if (szSA+szSC > size(TabSegArc)) &
            call XABORT("InterSgSg_V2 : memory problem (1)")
       ! Transfer
       tabSegArc(szSA+1:szSA+szSC) = tabSegCut(1:szSC)
       szSA = szSA + szSC
       szSC = 0
    enddo B_A2

    deallocate(tabSegCut, tabStrCut)
  end subroutine setBoundCut_V2

  subroutine setBoundCut(nbCut,szSA)
    integer,intent(in)    :: nbCut
    integer,intent(inout) :: szSA

    type(t_segArc)    :: sa,sgi,sg,ar
    type(segArcArrayTer),dimension(:),allocatable :: tmpTab
    double precision  :: intx,inty,pt1x,pt1y,pt2x,pt2y,angl
    integer           :: i,j,sizeTmp,nbPtInter,sztmp
    logical           :: oInSa,eInSa,oInSgi,eInSgi

    ! programmation defensive
    integer :: taille_table_tmpTab

    !on coupe le domaine de maniere symetrique
    !  et on garde les elements a gauche de l'axe de coupe
    !copie dans un tableau temporaire des segments
    !et nettoyage du tableau global

    sztmp = szSA*10
    ! programmation defensive 
    taille_table_tmpTab = sztmp
    allocate(tmpTab(sztmp),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: setBoundCut(1) => allocation pb")
    do i = 1,nbCut
       tmpTab(i)%sa=tabSegArc(szSA-nbCut+i)
       tmpTab(i)%keep=.true. ; tmpTab(i)%cl=.true.
    end do
    do i = 1,szSA-nbCut
       tmpTab(i+nbCut)%sa=tabSegArc(i)
       tmpTab(i+nbCut)%keep=.true. ; tmpTab(i+nbCut)%cl=.false.
    end do
    sizeTmp=szSA

    !coupage des segments   
    j=0
    do
       j = j+1
       if (j>sizeTmp) exit
       if (.not. tmpTab(j)%cl .or. .not. tmpTab(j)%keep) cycle
       sa = tmpTab(j)%sa
       if (isEqualConst(sa%x,sa%dx).and.isEqualConst(sa%y,sa%dy)) then
          tmpTab(j)%keep=.false.
          cycle
       end if
       i=0
       do
          i = i+1
          if (i>sizeTmp) exit
          if ((.not. tmpTab(i)%keep) .or. i==j) cycle
          if (tmpTab(i)%sa%typ/=tseg) cycle
          sgi=tmpTab(i)%sa
          !cas ou les segments de coupe rencontrent d'autres segments sans les
          ! couper. exemple: les coins des carres pour une symetrie diagonnale
          oInSa=pointsAlignes(sgi%x,sgi%y,sa%x,sa%y,sa%dx,sa%dy).and. &
               & (isIn(sgi%x,sgi%y,sa)==2)
          eInSa=pointsAlignes(sgi%dx,sgi%dy,sa%x,sa%y,sa%dx,sa%dy).and. &
               & (isIn(sgi%dx,sgi%dy,sa)==2)
          if (oInSa.and.(.not.eInSa)) then
             tmpTab(j)%keep=.false.
             sizeTmp=sizeTmp+1 
             ! programmation defensive
             if (sizeTmp > taille_table_tmpTab) &
                  call XABORT("G2S : memory problem in routine SetBoundCut")
             tmpTab(sizeTmp)%sa=sa ; tmpTab(sizeTmp)%cl=.true.
             tmpTab(sizeTmp)%keep=.true.
             if (estAGauche(sa%x,sa%y,sa%dx,sa%dy,sgi%dx,sgi%dy)) then
                if (sgi%mixd/=fooMix) tmpTab(sizeTmp)%sa%mixg = sgi%mixd
             else
                if (sgi%mixg/=fooMix) tmpTab(sizeTmp)%sa%mixg = sgi%mixg
             end if
             tmpTab(sizeTmp)%sa%x = sgi%x ; tmpTab(sizeTmp)%sa%y = sgi%y
             sizeTmp=sizeTmp+1 
             ! programmation defensive
             if (sizeTmp > taille_table_tmpTab) &
                  call XABORT("G2S : memory problem in routine SetBoundCut")
             tmpTab(sizeTmp)%sa=sa ; tmpTab(sizeTmp)%cl=.true.
             tmpTab(sizeTmp)%keep=.true.
             if (estAGauche(sa%x,sa%y,sa%dx,sa%dy,sgi%dx,sgi%dy)) then
                if (sgi%mixg/=fooMix) tmpTab(sizeTmp)%sa%mixg = sgi%mixg
             else
                if (sgi%mixd/=fooMix) tmpTab(sizeTmp)%sa%mixg = sgi%mixd
             end if
             tmpTab(sizeTmp)%sa%dx = sgi%x ; tmpTab(sizeTmp)%sa%dy = sgi%y
             exit
          else if ((.not.oInSa).and.eInSa) then
             tmpTab(j)%keep=.false.
             sizeTmp=sizeTmp+1 
             ! programmation defensive
             if (sizeTmp > taille_table_tmpTab) &
                  call XABORT("G2S : memory problem in routine SetBoundCut")
             tmpTab(sizeTmp)%sa=sa ; tmpTab(sizeTmp)%cl=.true.
             tmpTab(sizeTmp)%keep=.true.
             if (estAGauche(sa%x,sa%y,sa%dx,sa%dy,sgi%x,sgi%y)) then
                !write(*,*) "zaza sbc 5", i, j, sizetmp, sgi%mixG
                if (sgi%mixg/=fooMix) tmpTab(sizeTmp)%sa%mixg = sgi%mixg
             else
                !write(*,*) "zaza sbc 6" ,i, j, sizetmp, sgi%mixd
                if (sgi%mixd/=fooMix) tmpTab(sizeTmp)%sa%mixg = sgi%mixd
             end if
             tmpTab(sizeTmp)%sa%x = sgi%dx ; tmpTab(sizeTmp)%sa%y = sgi%dy
             !print*, tmpTab(sizeTmp)%p%x ,tmpTab(sizeTmp)%p%y ,tmptab(sizeTmp)%p%dx ,tmpTab(sizeTmp)%p%dy 
             !print*, tmpTab(sizeTmp)%p%mixg, tmpTab(sizeTmp)%p%mixd
             !print*, sa%x, sa%y, sa%dx, sa%dy
             !print*, sgi%x, sgi%y, sgi%dx, sgi%dy
              sizeTmp=sizeTmp+1 
             ! programmation defensive
             if (sizeTmp > taille_table_tmpTab) &
                  call XABORT("G2S : memory problem in routine SetBoundCut")
             tmpTab(sizeTmp)%sa=sa ; tmpTab(sizeTmp)%cl=.true.
             tmpTab(sizeTmp)%keep=.true.
             if (estAGauche(sa%x,sa%y,sa%dx,sa%dy,sgi%x,sgi%y)) then
                !write(*,*) "zaza sbc 7", i, j, sizetmp, sgi%mixd
                if (sgi%mixd/=fooMix) tmpTab(sizeTmp)%sa%mixg = sgi%mixd
             else
                !write(*,*) "zaza sbc 8", i, j, sizetmp, sgi%mixG
                if (sgi%mixg/=fooMix) tmpTab(sizeTmp)%sa%mixg = sgi%mixg
             end if
             tmpTab(sizeTmp)%sa%dx = sgi%dx ; tmpTab(sizeTmp)%sa%dy = sgi%dy
             !print*, tmpTab(sizeTmp)%p%x ,tmpTab(sizeTmp)%p%y ,tmpTab(sizeTmp)%p%dx ,tmpTab(sizeTmp)%p%dy 
             !print*, tmpTab(sizeTmp)%p%mixg, tmpTab(sizeTmp)%p%mixd
             !print*, sa%x, sa%y, sa%dx, sa%dy
             !print*, sgi%x, sgi%y, sgi%dx, sgi%dy
             exit
          end if
          !cas ou les segments a couper interceptent le segment de coupe a
          ! une de ses extremites. Exemple: bord du domaine pour la symetrie
          ! non diagonale sur un carre
          ! Remarque: pas d'exit , mais un cycle a la fin pour que le segment
          ! de coupe puisse etre encore utilise
          oInSgi=pointsAlignes(sa%x,sa%y,sgi%x,sgi%y,sgi%dx,sgi%dy).and. &
               & (isIn(sa%x,sa%y,sgi)==2)
          eInSgi=pointsAlignes(sa%dx,sa%dy,sgi%x,sgi%y,sgi%dx,sgi%dy).and. &
               & (isIn(sa%dx,sa%dy,sgi)==2)
          if (oInSgi.and.(.not.eInSgi)) then
             tmpTab(i)%keep=.false.
             sizeTmp=sizeTmp+1 
             ! programmation defensive
             if (sizeTmp > taille_table_tmpTab) &
                  call XABORT("G2S : memory problem in routine SetBoundCut")
             tmpTab(sizeTmp)%sa=sgi ; tmpTab(sizeTmp)%cl=tmpTab(i)%cl
             if (estAGauche(sa%x,sa%y,sa%dx,sa%dy,sgi%x,sgi%y)) then
                tmpTab(sizeTmp)%sa%dx=sa%x ;  tmpTab(sizeTmp)%sa%dy=sa%y
                !write(*,*) "zaza sbc 9", i, j, sizetmp, sgi%mixg
                if (sgi%mixg/=fooMix) tmpTab(j)%sa%mixg=sgi%mixg
             else
                tmpTab(sizeTmp)%sa%x=sa%x ;  tmpTab(sizeTmp)%sa%y=sa%y
                !write(*,*) "zaza sbc 10", i, j, sizetmp, sgi%mixd
                if (sgi%mixd/=fooMix) tmpTab(j)%sa%mixg=sgi%mixd
             end if
             tmpTab(sizeTmp)%keep=.true.
             cycle
          else if (eInSgi.and.(.not.oInSgi)) then
             tmpTab(i)%keep=.false.
             sizeTmp=sizeTmp+1 
             ! programmation defensive
             if (sizeTmp > taille_table_tmpTab) &
                  call XABORT("G2S : memory problem in routine SetBoundCut")
             tmpTab(sizeTmp)%sa=sgi ; tmpTab(sizeTmp)%cl=tmpTab(i)%cl
             if (estAGauche(sa%x,sa%y,sa%dx,sa%dy,sgi%x,sgi%y)) then
                tmpTab(sizeTmp)%sa%dx=sa%dx ;  tmpTab(sizeTmp)%sa%dy=sa%dy
                !write(*,*) "zaza sbc 11", i, j, sizetmp, sgi%mixd
                if (sgi%mixd/=fooMix) tmpTab(j)%sa%mixg=sgi%mixd
             else
                tmpTab(sizeTmp)%sa%x=sa%dx ;  tmpTab(sizeTmp)%sa%y=sa%dy
                !write(*,*) "zaza sbc 12", i, j, sizetmp, sgi%mixg
                if (sgi%mixg/=fooMix) tmpTab(j)%sa%mixg=sgi%mixg
             end if
             tmpTab(sizeTmp)%keep=.true.
             cycle
          end if
          !cas ou on coupe vraiment
          if ( interSgSg(sa,sgi,intx,inty)) then
             tmpTab(j)%keep=.false. ; tmpTab(i)%keep=.false.
             sizeTmp=sizeTmp+1 ; 
             tmpTab(sizeTmp)%sa=sgi ; tmpTab(sizeTmp)%cl=tmpTab(i)%cl
             tmpTab(sizeTmp)%sa%x=intx ; tmpTab(sizeTmp)%sa%y=inty
             sizeTmp=sizeTmp+1 
             ! programmation defensive
             if (sizeTmp > taille_table_tmpTab) &
                  call XABORT("G2S : memory problem in routine SetBoundCut")

             tmpTab(sizeTmp)%sa=sgi ; tmpTab(sizeTmp)%cl=tmpTab(i)%cl
             tmpTab(sizeTmp)%sa%dx=intx ; tmpTab(sizeTmp)%sa%dy=inty
             tmpTab(sizeTmp-1)%keep = estAGauche(sa%x,sa%y,sa%dx,sa%dy, &
                  & sgi%dx,sgi%dy)
             tmpTab(sizeTmp)%keep = .not.tmpTab(sizeTmp-1)%keep
             sizeTmp=sizeTmp+1 
             ! programmation defensive
             if (sizeTmp > taille_table_tmpTab) &
                  call XABORT("G2S : memory problem in routine SetBoundCut")
             tmpTab(sizeTmp)%sa=sa ; tmpTab(sizeTmp)%keep=.true.
             tmpTab(sizeTmp)%cl=.true.
             tmpTab(sizeTmp)%sa%x=intx ; tmpTab(sizeTmp)%sa%y=inty
             sizeTmp=sizeTmp+1 
             ! programmation defensive
             if (sizeTmp > taille_table_tmpTab) &
                  call XABORT("G2S : memory problem in routine SetBoundCut")

             tmpTab(sizeTmp)%sa=sa ; tmpTab(sizeTmp)%keep=.true.
             tmpTab(sizeTmp)%cl=.true.
             tmpTab(sizeTmp)%sa%dx=intx ; tmpTab(sizeTmp)%sa%dy=inty
             if (tmpTab(sizeTmp-2)%keep) then
                tmpTab(sizeTmp)%sa%mixg=sgi%mixd
                tmpTab(sizeTmp-1)%sa%mixg=sgi%mixg
                !write(*,*) "zaza sbc 13", i, j, sizetmp, sgi%mixd,sgi%mixg
             else
                tmpTab(sizeTmp)%sa%mixg=sgi%mixg
                tmpTab(sizeTmp-1)%sa%mixg=sgi%mixd
                !write(*,*) "zaza sbc 14" ,i, j, sizetmp, sgi%mixd,sgi%mixg
             end if
             exit
          end if
       end do
    end do
    do j = 1,nbCut
       sa = tmpTab(j)%sa
       do i = nbCut+1,sizeTmp
          if (.not. tmpTab(i)%keep) cycle
          sgi=tmpTab(i)%sa
          if (sgi%typ/=tseg) then
             !elimination des arcs et cercles du mauvais cote
             ! le centre de l'arc est a droite du segment de coupe
             ! => on l'elimine
             tmpTab(i)%keep=.not.estADroiteStrict(sa%x,sa%y,sa%dx,sa%dy, &
                  & sgi%x,sgi%y)
          else !elimination des segments du mauvais cote
             tmpTab(i)%keep=.not.(estADroiteStrict(sa%x,sa%y,sa%dx,sa%dy, &
                  & sgi%x,sgi%y) .or. &
                  & estADroiteStrict(sa%x,sa%y,sa%dx,sa%dy,sgi%dx,sgi%dy) )
          end if
       end do
    end do
    !coupage des arcs
    j=0
    do
       j = j+1
       if (j>sizeTmp) exit
       if ((.not. tmpTab(j)%cl).or.(.not. tmpTab(j)%keep)) cycle
       sg = tmpTab(j)%sa
       if (isEqual(tmpTab(j)%sa%x,tmpTab(j)%sa%dx) .and. &
            &  isEqual(tmpTab(j)%sa%y,tmpTab(j)%sa%dy)) then
          tmpTab(j)%keep=.false.
          cycle
       end if
       i=0
       do
          i = i+1
          if (i>sizeTmp) exit
          if (.not. tmpTab(i)%keep) cycle
          if (tmpTab(i)%sa%typ==tseg) cycle
          ar=tmpTab(i)%sa
          nbPtInter = abs(interSgAr(sg,ar,pt1x,pt1y,pt2x,pt2y))
          if (nbPtInter==2) then
             !on coupe le segment en trois morceaux
             tmpTab(j)%keep=.false.
             sizeTmp=sizeTmp+1 
             ! programmation defensive
             if (sizeTmp > taille_table_tmpTab) &
                  call XABORT("G2S : memory problem in routine SetBoundCut")
             !morceau 1
             tmpTab(sizeTmp)%sa=sg ; tmpTab(sizeTmp)%cl=.true.
             tmpTab(sizeTmp)%keep=.true.
             tmpTab(sizeTmp)%sa%dx=pt1x ; tmpTab(sizeTmp)%sa%dy=pt1y
             tmpTab(sizeTmp)%sa%mixg=ar%mixd
             sizeTmp=sizeTmp+1 
             ! programmation defensive
             if (sizeTmp > taille_table_tmpTab) &
                  call XABORT("G2S : memory problem in routine SetBoundCut")
             !morceau 2
             tmpTab(sizeTmp)%sa=sg ; tmpTab(sizeTmp)%cl=.true.
             tmpTab(sizeTmp)%keep=.true.
             tmpTab(sizeTmp)%sa%x=pt1x ; tmpTab(sizeTmp)%sa%y=pt1y
             tmpTab(sizeTmp)%sa%dx=pt2x ; tmpTab(sizeTmp)%sa%dy=pt2y
             tmpTab(sizeTmp)%sa%mixg=ar%mixg
             sizeTmp=sizeTmp+1 
             ! programmation defensive
             if (sizeTmp > taille_table_tmpTab) &
                  call XABORT("G2S : memory problem in routine SetBoundCut")
             !morceau 3
             tmpTab(sizeTmp)%sa=sg ; tmpTab(sizeTmp)%cl=.true.
             tmpTab(sizeTmp)%keep=.true.
             tmpTab(sizeTmp)%sa%x=pt2x ; tmpTab(sizeTmp)%sa%y=pt2y
             tmpTab(sizeTmp)%sa%mixg=ar%mixd
             !on coupe l'arc en 2 morceaux et on n'en garde qu'un
             tmpTab(i)%keep=.false.
             if (ar%typ==tcer) then
                sizeTmp=sizeTmp+1
                tmpTab(sizeTmp)%sa=ar ; tmpTab(sizeTmp)%cl=.false.
                tmpTab(sizeTmp)%keep=.true.
                tmpTab(sizeTmp)%sa%typ=tarc
                tmpTab(sizeTmp)%sa%a=calculeAngle(ar%x,ar%y,pt2x,pt2y)
                tmpTab(sizeTmp)%sa%b=calculeAngle(ar%x,ar%y,pt1x,pt1y)
             else
                ! il s'agit ici d'un arc centre sur le segment et le coupant 2
                ! fois, dans une geometrie symetrique... => a priori,
                !  cela ne doit jamais arriver
                call XABORT("G2S : Error, not symetrical geometry")
             end if
             exit
          else if (nbPtInter==1) then
             !on coupe le segment en deux morceaux
             tmpTab(j)%keep=.false.
             sizeTmp=sizeTmp+1 
             ! programmation defensive
             if (sizeTmp > taille_table_tmpTab) &
                  call XABORT("G2S : memory problem in routine SetBoundCut") 
             !morceau 1
             tmpTab(sizeTmp)%sa=sg ; tmpTab(sizeTmp)%cl=.true.
             tmpTab(sizeTmp)%keep=.true.
             tmpTab(sizeTmp)%sa%dx=pt1x ; tmpTab(sizeTmp)%sa%dy=pt1y
             sizeTmp=sizeTmp+1 
             ! programmation defensive
             if (sizeTmp > taille_table_tmpTab) &
                  call XABORT("G2S : memory problem in routine SetBoundCut")
             !morceau 2
             tmpTab(sizeTmp)%sa=sg ; tmpTab(sizeTmp)%cl=.true.
             tmpTab(sizeTmp)%keep=.true.
             tmpTab(sizeTmp)%sa%x=pt1x ; tmpTab(sizeTmp)%sa%y=pt1y
             !determination des milieux pour les segments selon l'ordre
             !relatif des points Osg, InterSgAr, Car... 
             !(en fait selon le sens relatif des segments OI et CI)
             if (((pt1x-sg%x)*(pt1x-ar%x)+(pt1y-sg%y)*(pt1y-ar%y))>0.) then
                !dans le meme sens => premier segment a l'interieur
                tmpTab(sizeTmp-1)%sa%mixg=ar%mixg
                tmpTab(sizeTmp)%sa%mixg=ar%mixd
             else !=> premier segment a l'exterieur
                tmpTab(sizeTmp-1)%sa%mixg=ar%mixd
                tmpTab(sizeTmp)%sa%mixg=ar%mixg
             end if
             !on coupe l'arc en 1 ou 2 morceaux et on n'en garde qu'un
             tmpTab(i)%keep=.false.
             if (ar%typ==tcer) then
                sizeTmp=sizeTmp+1 
                ! programmation defensive
                if (sizeTmp > taille_table_tmpTab) &
                     call XABORT("G2S : memory problem in routine SetBoundCut")
                tmpTab(sizeTmp)%sa=ar ; tmpTab(sizeTmp)%cl=.false.
                tmpTab(sizeTmp)%keep=.true.
                tmpTab(sizeTmp)%sa%typ=tarc
                angl=calculeAngle(ar%x,ar%y,pt1x,pt1y)
                if (isPI(angl)) then
                   tmpTab(sizeTmp)%sa%a=-pi_c
                   tmpTab(sizeTmp)%sa%b=pi_c
                else
                   tmpTab(sizeTmp)%sa%a=angl
                   tmpTab(sizeTmp)%sa%b=angl
                end if
             else
                sizeTmp=sizeTmp+1 
                tmpTab(sizeTmp)%sa=ar ; tmpTab(sizeTmp)%cl=.false.
                tmpTab(sizeTmp)%keep=.true.
                !pour savoir quelle partie on garde, on se sert du meme test
                !que pour les milieux
                if (((pt1x-sg%x)*(pt1x-ar%x)+(pt1y-sg%y)*(pt1y-ar%y))>0.) then
                   !dans le meme sens => on garde la deuxieme moitie de l'arc
                   tmpTab(sizeTmp)%sa%a=calculeAngle(ar%x,ar%y,pt1x,pt1y)
                else !=> on garde la premiere moitie de l'arc
                   tmpTab(sizeTmp)%sa%b=calculeAngle(ar%x,ar%y,pt1x,pt1y)
                end if
             end if
             exit
          else !pas d'intersection mais l'arc est peut-etre a eliminer
             if (ar%typ/=tarc) cycle
             call extremitesArc(ar,pt1x,pt1y,pt2x,pt2y)
             if ( ( estADroiteStrict(sg%x,sg%y,sg%dx,sg%dy,pt1x,pt1y) .or.   &
                  & estADroiteStrict(sg%x,sg%y,sg%dx,sg%dy,pt2x,pt2y) )      &
                  &                        .and.     .not.                   &
                  & ( estAGaucheStrict(sg%x,sg%y,sg%dx,sg%dy,pt1x,pt1y) .or. &
                  &   estAGaucheStrict(sg%x,sg%y,sg%dx,sg%dy,pt2x,pt2y) ) )  &
                  & tmpTab(i)%keep=.false.
          end if
       end do
    end do
    !elimination des segments nuls eventuellement restants
    do i = 1,sizeTmp
       sa = tmpTab(i)%sa
       if (sa%typ==tseg) tmpTab(i)%keep=tmpTab(i)%keep .and. &
            & .not.(isEqual(tmpTab(i)%sa%x,tmpTab(i)%sa%dx) .and.&
            &       isEqual(tmpTab(i)%sa%y,tmpTab(i)%sa%dy))
    end do
    !elimination des segments non de coupe se trouvant sur les segments de
    ! coupe (ca arrive pour les sectorisations)
    do j = 1,sizeTmp
       if (.not.(tmpTab(j)%cl .and. tmpTab(j)%keep)) cycle
       sa = tmpTab(j)%sa
       if (sa%typ/=tseg) cycle
       do i = 1,sizeTmp
          if (.not. tmpTab(i)%keep .or. tmpTab(i)%cl ) cycle
          sgi=tmpTab(i)%sa
          if (sgi%typ/=tseg) cycle
          if (.not. estColineaire(sa,sgi)) cycle
          if (.not. pointsAlignes(sa%x,sa%y,sa%dx,sa%dy,sgi%x,sgi%y)) cycle
          !on renverse le sens de sgi si besoin
          if (.not. isSameWay(sa,sgi)) sgi=turnBackSide(tmpTab(i)%sa)
          if (isEqualConst(sa%x,sgi%x).and.isEqualConst(sa%dx,sgi%dx).and. &
               & isEqualConst(sa%y,sgi%y).and.isEqualConst(sa%dy,sgi%dy)) then
             tmpTab(j)%sa%mixg=sgi%mixg
             tmpTab(i)%keep=.false.
             exit
          end if
       end do
    end do
    !recopie dans le tableau global des elements a conserver
    szSA=0
    ! programmation defensive
    if (sizeTmp > size(tabSegArc)) &
         call XABORT("G2S : memory problem in routine SetBoundCut (2)")
    do i = 1,sizeTmp
       if (tmpTab(i)%keep) then
          szSA=szSA+1
          tabSegArc(szSA)=tmpTab(i)%sa
       end if
    end do
    deallocate(tmpTab)
  end subroutine setBoundCut

  subroutine prepareSALBCData(szSA,nbCLP)
    integer,intent(in) :: szSA,nbCLP

    integer        :: i,j,k,sunsetType,tmpTabSize,indice
    logical        :: newBC
    type(t_segArc) :: sa
    integer,dimension(:),allocatable :: tmpTabElemNb

    !enregistrement des numeros des differentes conditions limites
!    write(*,*) 'nbCLP :',nbCLP
    do i = 1,nbCLP
       SALbCDataTab(i)%sunsetType = -1
       sunsetType = 0
       do j = 1,szSA
          sa = tabSegArc(j)
          if (sa%mixg==fooMix .or. sa%mixd==fooMix) then
             cycle
          else if (sa%mixg<0) then
             sunsetType = -sa%mixg
          else if (sa%mixd<0) then
             sunsetType = -sa%mixd
          else
             cycle
          end if
          newBC = .true.
          do k = 1,i
             newBC = (newBC .and. (sunsetType /= SALbCDataTab(k)%sunsetType) &
                  & .and. (mod(sunsetType,100) /= B_Void) )
          end do
          if (.not. newBC) cycle
          !on a un type de condition limite non encore pris en compte
          SALbCDataTab(i)%sunsetType = sunsetType
          exit
       end do
    end do
    !creation des donnees SAL pour chaque CL
    do i = 1,nbCLP
       SALbCDataTab(i)%nber = 0
       sunsetType = mod(SALbCDataTab(i)%sunsetType,100)
       indice = SALbCDataTab(i)%sunsetType / 100
       select case(geomTyp)
       case(RecTyp)
          select case(sunsetType)
          case(B_Refl,B_Ssym)
             SALbCDataTab(i)%albedo = 1.
             !SALbCDataTab(i)%SALtype = 1 !ne marche pas pour le moment
             ! => on passe par une symetrie, mais exterieure
             SALbCDataTab(i)%SALtype = 4
             select case(indice)
             case(1)
                SALbCDataTab(i)%cx = 0.
                SALbCDataTab(i)%cy = 0.
                SALbCDataTab(i)%angle = 90.
             case(2)
                SALbCDataTab(i)%cx = real(0.5*bCData%sidexy(1)-bCData%toOrig_xy(1))
                SALbCDataTab(i)%cy = 0.
                SALbCDataTab(i)%angle = 90.
             case(3)
                SALbCDataTab(i)%cx = 0.
                SALbCDataTab(i)%cy = 0.
                SALbCDataTab(i)%angle = 0.
             case(4)
                SALbCDataTab(i)%cx = 0.
                SALbCDataTab(i)%cy = real(0.5*bCData%sidexy(2)-bCData%toOrig_xy(2))
                SALbCDataTab(i)%angle = 0.
             end select
          case(B_Tran)
             SALbCDataTab(i)%SALtype = 2
             select case(indice)
             case(1)
                SALbCDataTab(i)%tx = real(bCData%sidexy(1))
                SALbCDataTab(i)%ty = 0.
             case(2)
                SALbCDataTab(i)%tx = real(-bCData%sidexy(1))
                SALbCDataTab(i)%ty = 0.
             case(3)
                SALbCDataTab(i)%tx = 0.
                SALbCDataTab(i)%ty = real(bCData%sidexy(2))
             case(4)
                SALbCDataTab(i)%tx = 0.
                SALbCDataTab(i)%ty = real(-bCData%sidexy(2))
             end select
          case(B_Diag)
             SALbCDataTab(i)%SALtype = 4
             SALbCDataTab(i)%cx = 0.
             SALbCDataTab(i)%cy = 0.
             SALbCDataTab(i)%angle = 45.          
          case(B_Syme)
             SALbCDataTab(i)%SALtype = 4
             select case(indice)
             case(1)
                SALbCDataTab(i)%cx = 0.
                SALbCDataTab(i)%cy = 0.
                SALbCDataTab(i)%angle = 90.
             case(2)
                SALbCDataTab(i)%cx = real(bCData%minmaxXY(3)-bCData%toOrig_xy(1))
                SALbCDataTab(i)%cy = 0.
                SALbCDataTab(i)%angle = 90.
             case(3)
                SALbCDataTab(i)%cx = 0.
                SALbCDataTab(i)%cy = 0.
                SALbCDataTab(i)%angle = 0.
             case(4)
                SALbCDataTab(i)%cx = 0.
                SALbCDataTab(i)%cy = real(bCData%minmaxXY(4)-bCData%toOrig_xy(2))
                SALbCDataTab(i)%angle = 0.
             end select
          case(B_Albe)
             SALbCDataTab(i)%SALtype = 0
             SALbCDataTab(i)%albedo = real(bCData%albedo(indice))
          case(B_Zero)
             SALbCDataTab(i)%SALtype = 0
             SALbCDataTab(i)%albedo = 0.
          case(B_Pi_2)
             SALbCDataTab(i)%SALtype = 3
             SALbCDataTab(i)%cx = 0.
             SALbCDataTab(i)%cy = 0.
             if (indice==1) then
                SALbCDataTab(i)%angle = -90.
             else if (indice==3) then
                SALbCDataTab(i)%angle = 90.
             end if
          case(B_Pi)
             SALbCDataTab(i)%SALtype = 3
             SALbCDataTab(i)%cx = real(-bCData%toOrig_xy(1))
             SALbCDataTab(i)%cy = real(-bCData%toOrig_xy(2))
             SALbCDataTab(i)%angle = 180.
          end select
       case(HexTyp)
          SALbCDataTab(i)%SALtype = 4
          select case(indice)
          case(1)
             if (bCData%iHex==H_S30) then
                SALbCDataTab(i)%cx    = 0.
                SALbCDataTab(i)%cy    = 0.
                SALbCDataTab(i)%angle = 0.
             else if (bCData%iHex==H_SA60) then
                SALbCDataTab(i)%cx    = real(-bCData%toOrig_xy(1))
                SALbCDataTab(i)%cy    = real(-bCData%toOrig_xy(2))
                SALbCDataTab(i)%angle = -30.
             else
                call XABORT("G2S: internal error in prepareSALBCData(1)")
             end if
          case(2)
             SALbCDataTab(i)%angle = 30.
             if (bCData%iHex==H_S30) then
                SALbCDataTab(i)%cx    = 0.
                SALbCDataTab(i)%cy    = 0.
             else if (bCData%iHex==H_SA60) then
                SALbCDataTab(i)%cx    = real(-bCData%toOrig_xy(1))
                SALbCDataTab(i)%cy    = real(-bCData%toOrig_xy(2))
             else
                call XABORT("G2S: internal error in prepareSALBCData(2)")
             end if
          case(3)
             SALbCDataTab(i)%cy    = 0.
             SALbCDataTab(i)%angle = 90.
             if (bCData%iHex==H_S30) then
                SALbCDataTab(i)%cx = real(.5*sqrt(3.)*bCData%sidexy(2))
             else if (bCData%iHex==H_SA60) then
                SALbCDataTab(i)%cx = real(.5*sqrt(3.)*bCData%sidexy(2) &
                     -bCData%toOrig_xy(1))
             else
                call XABORT("G2S: internal error in prepareSALBCData(3)")
             end if
          end select
       case(TriaTyp)
          select case(sunsetType)
          case(B_Refl)
             SALbCDataTab(i)%SALtype = 1
             SALbCDataTab(i)%albedo = 1.
          case(B_Albe)
             SALbCDataTab(i)%SALtype = 0
             SALbCDataTab(i)%albedo = real(bCData%albedo(indice))          
          case(B_Zero)
             SALbCDataTab(i)%SALtype = 0
             SALbCDataTab(i)%albedo = 0.
          case(B_Syme)
             SALbCDataTab(i)%SALtype = 4
             select case(indice)
             case(1)
                if (bCData%iTri==T_S30) then
                   SALbCDataTab(i)%cx    = 0.
                   SALbCDataTab(i)%cy    = 0.
                   SALbCDataTab(i)%angle = 0.
                else if (bCData%iTri==T_SA60) then
                   SALbCDataTab(i)%cx    = real(-bCData%toOrig_xy(1))
                   SALbCDataTab(i)%cy    = real(-bCData%toOrig_xy(2))
                   SALbCDataTab(i)%angle = -30.
                else
                   call XABORT("G2S: internal error in prepareSALBCData(4)")
                end if
             case(2)
                if (bCData%iTri==T_S30) then
                   SALbCDataTab(i)%cx    = 0.
                   SALbCDataTab(i)%cy    = 0.
                   SALbCDataTab(i)%angle = 30.
                else if (bCData%iTri==T_SA60) then
                   SALbCDataTab(i)%cx    = real(-bCData%toOrig_xy(1))
                   SALbCDataTab(i)%cy    = real(-bCData%toOrig_xy(2))
                   SALbCDataTab(i)%angle = 30.
                else
                   call XABORT("G2S: internal error in prepareSALBCData(5)")
                end if
             case(3)
                SALbCDataTab(i)%cx    = real(bCData%minmaxXY(3)-bCData%toOrig_xy(1))
                SALbCDataTab(i)%cy    = 0.
                SALbCDataTab(i)%angle = 90.
             end select
          end select
       end select
    end do
    !remplissage du tableau des indices des segArc concernes par la condition
    do i = 1,nbCLP
       sunsetType = SALbCDataTab(i)%sunsetType
       allocate(tmpTabElemNb(szSA),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: prepareSALBCData(1) => allocation pb")
       tmpTabSize = 0
       do j = 1,szSA
          sa = tabSegArc(j)
          if (sa%mixg==-sunsetType .or. sa%mixd==-sunsetType) then
             tmpTabSize = tmpTabSize + 1
             tmpTabElemNb(tmpTabSize) = j
          end if
       end do
       SALbCDataTab(i)%nber = tmpTabSize
       allocate(SALbCDataTab(i)%elemNb(tmpTabSize),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: prepareSALBCData(2) => allocation pb")
       do j = 1,tmpTabSize
          SALbCDataTab(i)%elemNb(j) = tmpTabElemNb(j)
       end do
       deallocate(tmpTabElemNb)
    end do
  end subroutine prepareSALBCData


  subroutine appliBoundariConditionsForHex(geoIp,szSA,nbCLP)
    type(c_ptr),intent(in):: geoIp
    integer,intent(inout) :: szSA
    integer,intent(out)   :: nbCLP

    real,dimension(6) :: tmpTab6
    type(c_ptr)       :: ip
    integer           :: nbCut
    double precision  :: fx,fy,tx,ty,sqrt3_2
    double precision  :: rayIntCell,rayIntCore,rayExtCell,rayExtCore
    type(t_segArc)    :: sg

    ! programmation defensive
    integer :: dimTabSegArc
    dimTabSegArc = size(tabSegArc)

    ip = geoIp
    call LCMGET(ip,'NCODE       ',bCData%bc)
    call LCMGET(ip,'ZCODE       ',tmpTab6) ; bCData%albedo=tmpTab6
    call LCMGET(ip,'ICODE       ',bCData%albInd)
    call LCMGET(ip,'IHEX        ',bCData%iHex)
    bCData%sidexy(:)=0.0d0
    sqrt3_2 = 5.d-1*sqrt(3.d0)
    nbCut = 0 ; nbCLP = 0
    rayExtCell = bCData%sidexy(1)
    rayIntCell = sqrt3_2*rayExtCell
    rayExtCore = bCData%sidexy(2)
    rayIntCore = sqrt3_2*rayExtCore
    select case(bCData%iHex)
    case(H_S30)
       !1
       nbCut=nbCut+1 ; nbCLP=nbCLP+1
       szSA=szSA+1           
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine &
            &appliBoundaryConditionsForHex")
       fx = -rayExtCell           ; fy = 0.d0
       tx = rayIntCore+rayExtCell ; ty = 0.d0
       sg = createSeg(fx,fy,tx,ty,fooMix,-100-B_Syme)
       tabSegArc(szSA) = sg
       !2
       nbCut=nbCut+1 ; nbCLP=nbCLP+1
       szSA=szSA+1            
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine &
            &appliBoundaryConditionsForHex")
       fx = rayIntCore+sqrt3_2*rayIntCell
       fy = 5.d-1*(rayExtCore+rayIntCell)
       tx = 0.d0 ; ty = 0.d0
       sg = createSeg(fx,fy,tx,ty,fooMix,-200-B_Syme)
       tabSegArc(szSA) = sg
       if (bCData%bc(1)==B_Syme) then
          !3
          nbCut=nbCut+1 ; nbCLP=nbCLP+1
          szSA=szSA+1            
          ! programmation defensive
          if (szSA > dimTabSegArc)&
               call XABORT("G2S: memory problem in routine &
               &appliBoundaryConditionsForHex")
          fx = rayIntCore ; fy = 0.d0
          tx = fx         ; ty = 5.d-1*rayExtCore
          sg = createSeg(fx,fy,tx,ty,fooMix,-300-bCData%bc(1))
          tabSegArc(szSA) = sg
       end if
    case(H_SA60)
       !1
       nbCut=nbCut+1 ; nbCLP=nbCLP+1
       szSA=szSA+1            
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine &
            &appliBoundaryConditionsForHex")
       fx = -sqrt3_2*rayIntCell
       fy = 5.d-1*rayIntCell
       tx = rayIntCore+sqrt3_2*rayIntCell
       ty = -5.d-1*(rayExtCore+rayIntCell)
       sg = createSeg(fx,fy,tx,ty,fooMix,-100-B_Syme)
       tabSegArc(szSA) = sg
       !2
       nbCut=nbCut+1 ; nbCLP=nbCLP+1
       szSA=szSA+1            
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine &
            &appliBoundaryConditionsForHex")
       fx = rayIntCore+sqrt3_2*rayIntCell
       fy = 5.d-1*(rayExtCore+rayIntCell)
       tx = 0.d0 ; ty = 0.d0
       sg = createSeg(fx,fy,tx,ty,fooMix,-200-B_Syme)
       tabSegArc(szSA) = sg
       if (bCData%bc(1)==B_Syme) then
          !3
          nbCut=nbCut+1 ; nbCLP=nbCLP+1
          szSA=szSA+1            
          ! programmation defensive
          if (szSA > dimTabSegArc)&
               call XABORT("G2S: memory problem in routine &
               &appliBoundaryConditionsForHex")
          fx = rayIntCore ; fy = -5.d-1*rayExtCore
          tx = fx         ; ty = -fy
          sg = createSeg(fx,fy,tx,ty,fooMix,-300-bCData%bc(1))
          tabSegArc(szSA) = sg
       end if
    case(H_SB60)
       !1
       nbCut=nbCut+1 ; nbCLP=nbCLP+1
       szSA=szSA+1            
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine &
            &appliBoundaryConditionsForHex")
       fx = -rayExtCell           ; fy = 0.d0
       tx = rayIntCore+rayExtCell ; ty = 0.d0
       sg = createSeg(fx,fy,tx,ty,fooMix,-100-B_Syme)
       tabSegArc(szSA) = sg
       !2
       nbCut=nbCut+1 ; nbCLP=nbCLP+1
       szSA=szSA+1            
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine &
            &appliBoundaryConditionsForHex")
       fx = 5.d-1*(rayIntCore+rayExtCell)
       fy = sqrt3_2*(rayIntCore+rayExtCell) 
       tx = 0.d0 ; ty = 0.d0
       sg = createSeg(fx,fy,tx,ty,fooMix,-200-B_Syme)
       tabSegArc(szSA) = sg       
    case(H_S90)
       !1
       nbCut=nbCut+1 ; nbCLP=nbCLP+1
       szSA=szSA+1            
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine &
            &appliBoundaryConditionsForHex")
       fx = -rayExtCell           ; fy = 0.d0
       tx = rayIntCore+rayExtCell ; ty = 0.d0
       sg = createSeg(fx,fy,tx,ty,fooMix,-100-B_Syme)
       tabSegArc(szSA) = sg
       !2
       nbCut=nbCut+1 ; nbCLP=nbCLP+1
       szSA=szSA+1            
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine &
            &appliBoundaryConditionsForHex")
       fx = 0.d0 ; fy = rayExtCore+rayIntCell
       tx = 0.d0 ; ty = 0.d0
       sg = createSeg(fx,fy,tx,ty,fooMix,-200-B_Syme)
       tabSegArc(szSA) = sg              
    case(H_R120)
       !!!appliquer les rotations adequat (pas simple)
    case(H_R180)
       !!!idem mais plus simple
    case(H_SA180)
       !1
       nbCut=nbCut+1 ; nbCLP=nbCLP+1
       szSA=szSA+1            
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine &
            &appliBoundaryConditionsForHex")
       fx = 0.d0 ; fy = rayExtCore+rayIntCell
       tx = 0.d0 ; ty = -fy
       sg = createSeg(fx,fy,tx,ty,fooMix,-100-B_Syme)
       tabSegArc(szSA) = sg                    
    case(H_SB180)
       !1
       nbCut=nbCut+1 ; nbCLP=nbCLP+1
       szSA=szSA+1            
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine &
            &appliBoundaryConditionsForHex")
       fx = -rayExtCore-rayIntCell ; fy = 0.d0
       tx = -fx                    ; ty = 0.d0
       sg = createSeg(fx,fy,tx,ty,fooMix,-100-B_Syme)
       tabSegArc(szSA) = sg       
    case(H_Complete)
    end select
    if (nbCut/=0) then
       call setBoundCut(nbCut,szSA)
    end if
  end subroutine appliBoundariConditionsForHex

  subroutine appliBoundariConditionsForTri(geoIp,szSA,nbCLP)
    type(c_ptr),intent(in):: geoIp
    integer,intent(inout) :: szSA
    integer,intent(out)   :: nbCLP

    real,dimension(2) :: tmpTab2
    real,dimension(4) :: tmpTab4
    real,dimension(6) :: tmpTab6
    type(c_ptr)       :: ip
    integer           :: nbCut,i
    double precision  :: dy
    type(t_segArc)    :: sg
    integer,dimension(40)         :: sv
    double precision,dimension(4) :: x,y
    
    ! programmation defensive
    integer :: dimTabSegArc 
    dimTabSegArc = size(tabSegArc)
    
    ip = geoIp
    call LCMGET(ip,'STATE-VECTOR',sv)
    call LCMGET(ip,'NCODE       ',bCData%bc)
    call LCMGET(ip,'ZCODE       ',tmpTab6) ; bCData%albedo=tmpTab6
    call LCMGET(ip,'ICODE       ',bCData%albInd)
    call LCMGET(ip,'ITRI        ',bCData%iTri)
    call LCMSIX(ip,'NEW-DATA    ',1)
    call LCMSIX(ip,'BOUND-DATA  ',1)
    call LCMGET(ip,'SIDEXY      ',tmpTab2) ; bCData%sidexy=tmpTab2
    call LCMGET(ip,'MINMAXXY    ',tmpTab4) ; bCData%minmaxXY=tmpTab4

    nbCut = 0
    nbCLP = 0
    select case(bCData%iTri)
    case(T_S30)
       !1 (on ne coupe pas car inutile)
       nbCLP=nbCLP+1
       call setBoundSide(0.d0,0.d0,bCData%sidexy(2),0.d0,100+B_Syme,szSA)
       !2
       nbCut=nbCut+1 ; nbCLP=nbCLP+1
       szSA=szSA+1            
       ! programmation defensive
       if (szSA > dimTabSegArc)&
            call XABORT("G2S: memory problem in routine &
            &appliBoundaryConditionsForTri")
       dy = bCData%sidexy(1) * sqrt3_2d
       if (estPaire(sv(4))) then
          x(1) = bCData%sidexy(2)
          y(1) = dy * sv(4)
       else
          x(1) = bCData%sidexy(2) - 2.5d-1*bCData%sidexy(1)
          y(1) = dy * (sv(4) - 5.d-1)
       end if
       sg = createSeg(x(1),y(1),0.d0,0.d0,fooMix,-200-B_Syme)
       tabSegArc(szSA) = sg
       if (bCData%bc(1)==B_Syme) then
          !3
          nbCut=nbCut+1 ; nbCLP=nbCLP+1
          szSA=szSA+1        ! programmation defensive
          if (szSA > dimTabSegArc)&
               call XABORT("G2S: memory problem in routine &
               &appliBoundaryConditionsForTri")

          if (estPaire(sv(4))) then
             x(2) = bCData%sidexy(2)
             y(2) = dy * sv(4)
          else
             x(2) = bCData%sidexy(2) - 5.d-1*bCData%sidexy(1)
             y(2) = dy * (sv(4) - 2.d0/3)
          end if
          bCData%minmaxXY(3)=x(2) !stockage de l'abscisse de l'axe de coupe
          sg = createSeg(x(2),0.d0,x(2),y(2),fooMix,-300-B_Syme)
          tabSegArc(szSA) = sg
       end if
    case(T_SA60)
       !1
       nbCut=nbCut+1 ; nbCLP=nbCLP+1
       szSA=szSA+1        
       ! programmation defensive
       if (szSA > dimTabSegArc) call XABORT("G2S: memory problem in routine &
               &appliBoundaryConditionsForTri")
       dy = bCData%sidexy(1) * sqrt3_2d
       if (estPaire(sv(4))) then
          x(1) = bCData%sidexy(2)
          y(1) = -dy * sv(4)
       else
          x(1) = bCData%sidexy(2) - 2.5d-1*bCData%sidexy(1)
          y(1) = -dy * (sv(4) - 5.d-1)
       end if
       sg = createSeg(0.d0,0.d0,x(1),y(1),fooMix,-100-B_Syme)
       tabSegArc(szSA) = sg
       !2
       nbCut=nbCut+1 ; nbCLP=nbCLP+1
       szSA=szSA+1        
       ! programmation defensive
       if (szSA > dimTabSegArc)&
              call XABORT("G2S: memory problem in routine &
               &appliBoundaryConditionsForTri")
       sg = createSeg(x(1),-y(1),0.d0,0.d0,fooMix,-200-B_Syme)
       tabSegArc(szSA) = sg
       if (bCData%bc(1)==B_Syme) then
          !3
          nbCut=nbCut+1 ; nbCLP=nbCLP+1
          szSA=szSA+1        ! programmation defensive
          if (szSA > dimTabSegArc)&
                 call XABORT("G2S: memory problem in routine &
               &appliBoundaryConditionsForTri")
          if (estPaire(sv(4))) then
             x(2) = bCData%sidexy(2)
             y(2) = dy * sv(4)
          else
             x(2) = bCData%sidexy(2) - 5.d-1*bCData%sidexy(1)
             y(2) = dy * (sv(4) - 2.d0/3)
          end if
          bCData%minmaxXY(3)=x(2) !stockage de l'abscisse de l'axe de coupe
          sg = createSeg(x(2),-y(2),x(2),y(2),fooMix,-300-B_Syme)
          tabSegArc(szSA) = sg
       end if
    case(T_ST60)
       x(1) = -5.d-1*bCData%sidexy(1) ; x(2) = 0. ; x(3) = -x(1) ; x(4) = x(1)
       y(1) = sqrt3_2d*x(1) ; y(2) = -y(1) ; y(3) = y(1) ; y(4) = y(1)
       do i = 1,3
          if (bCData%bc(i)/= B_Void) nbCLP=nbCLP+1
          call setBoundSide(x(i+1),y(i+1),x(i),y(i),bCData%bc(i)+100*i,szSA)
       end do
    case(T_Complete)
       y(1) = -5.d-1*bCData%sidexy(2) ; y(2) =  y(1)
       y(3) = -y(1)                   ; y(4) =  y(1)
       x(1) = -5.d-1*bCData%sidexy(1) ; x(2) = -x(1)
       x(3) =  x(2) + y(3)            ; x(4) =  x(1) + y(3)
       call setBoundSide(x(4),y(4),x(1),y(1),bCData%bc(1)+100,szSA)
       call setBoundSide(x(2),y(2),x(3),y(3),bCData%bc(2)+200,szSA)
       call setBoundSide(x(1),y(1),x(2),y(2),bCData%bc(3)+300,szSA)
       call setBoundSide(x(3),y(3),x(4),y(4),bCData%bc(4)+400,szSA)
    case default
       call XABORT("G2S: internal error in subroutine &
            &appliBoundariConditionsForTri")
    end select
    if (nbCut/=0) then
       call setBoundCut(nbCut,szSA)
    end if
  end subroutine appliBoundariConditionsForTri

  subroutine appliBoundariConditionsForTub(geoIp,nbCLP)
    type(c_ptr),intent(in):: geoIp
    integer,intent(out)   :: nbCLP

    real,dimension(6) :: tmpTab6

    nbCLP = 0
    call LCMGET(geoIp,'NCODE       ',bCData%bc)
    call LCMGET(geoIp,'ZCODE       ',tmpTab6) ; bCData%albedo=tmpTab6
    call LCMGET(geoIp,'ICODE       ',bCData%albInd)
  end subroutine appliBoundariConditionsForTub

end module boundCond
