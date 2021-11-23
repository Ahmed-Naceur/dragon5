!
!-----------------------------------------------------------------------
!
!Purpose:
! Generate a dataset for the Tripoli code.
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
! fonctions:
!  - initMC : initialisation des tableaux
!  - destroyMC : liberation des pointeurs
!  - createElemGene : constructeur d'un element generique
!  - isSameEG : teste l'equivalence de deux elements generiques
!  - isSameWayEG : teste si deux elements generiques ont meme orientation
!  - prepareMCData : preparation des donnes necessaires et remplissage
!                    des tableaux globaux
!  - generateTripoliFile : ecriture du fichier de donnees Tripoli4
!  - generateMCNPFile : ecriture du fichier de donnees MCNP
!  - generateSerpentFile : ecriture du fichier de donnees Serpent
!  - putOn80col : formattage sur 80 colonnes style MCNP, d'un buffer a afficher
!  - findParalleleWithTrans : trouve l'element geometrique a associer avec un
!                             autre, dans le cas d'une translation
!
!-----------------------------------------------------------------------
!
module monteCarlo
    use boundCond
    use cast
    use cellulePlaced
    use constType
    use constUtiles
    use segArc
    use generSAL
  
    implicit none

  type t_elemGene
     logical          :: isPlan !plan -> true      ; cylindre -> false
     double precision :: x,y    !     -> origine   ;          -> centre
     double precision :: dx,dy  !     -> extremite ;          -> 0.
     double precision :: r      !     -> 0.        ;          -> rayon
     integer          :: limite ! ==Tri_XXXX
  end type t_elemGene

  type(t_elemGene),dimension(:),allocatable,save :: tabEG

  type t_volume
     integer                      :: lg      !longueur des tableaux
     integer                      :: mix     !numero du milieu
     integer,dimension(:), allocatable :: indElem !indice des elements generiques
                                             !du contour
     logical,dimension(:), allocatable :: side    !true -> cote + ; false -> cote -
     integer,dimension(:), allocatable :: typCL   !condition limite de la face
  end type t_volume

  type(t_volume),dimension(:),allocatable,save   :: tabVolume

  integer,parameter :: Tri_Not=-1,Tri_Void=0,Tri_Refl=1,Tri_Trans=2,Tri_Cos=3

contains

  subroutine initMC(nbNode)
    integer,intent(in) :: nbNode

    allocate(tabVolume(nbNode),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: initMC => allocation pb")
    tabVolume(1:nbNode)%lg = 0
  end subroutine initMC

  subroutine destroyMC()
    integer :: i

    do i = 1,size(tabVolume)
       deallocate(tabVolume(i)%indElem)
       deallocate(tabVolume(i)%side)
       deallocate(tabVolume(i)%typCL)
    end do
    deallocate(tabEG,tabVolume)
  end subroutine destroyMC

  function createElemGene(sa,strDCL)
    type(t_segArc),intent(in)   :: sa
    character(len=4),intent(in) :: strDCL
    type(t_elemGene)            :: createElemGene

    integer :: sunsetCL

    createElemGene%isPlan = (sa%typ==tseg)
    createElemGene%x = sa%x
    createElemGene%y = sa%y
    if (createElemGene%isPlan) then
       createElemGene%dx = sa%dx
       createElemGene%dy = sa%dy
       createElemGene%r  = 0.d0
    else
       createElemGene%dx = 0.d0
       createElemGene%dy = 0.d0
       createElemGene%r  = sa%r
    end if
    createElemGene%limite = Tri_Not
    sunsetCL = minval((/sa%nodeg,sa%noded/))
    if (sunsetCL<=0) then
       sunsetCL = mod(-sunsetCL,100)
       select case(sunsetCL)
       case(B_Void,B_Zero)
          createElemGene%limite = Tri_Void
       case(B_Refl,B_Ssym,B_Syme,B_Diag)
          createElemGene%limite = Tri_Refl
       case(B_Tran)
          createElemGene%limite = Tri_Trans
       case(B_Albe)
          createElemGene%limite = Tri_Cos
       case(-fooMix,0)
          !il faut utiliser la CL par defaut
          if (strDCL=='ALBE') then
             createElemGene%limite = Tri_Cos
          else if (strDCL=='REFL') then
             createElemGene%limite = Tri_Refl
          else ! => VOID
             createElemGene%limite = Tri_Void
          end if
       case default
          call XABORT("G2MC: boundary condition not allowed")
       end select
    end if
  end function createElemGene

  function isSameEG(eg1,eg2)
    type(t_elemGene),intent(in) :: eg1,eg2
    logical                     :: isSameEG

    if (eg1%isPlan.neqv.eg2%isPlan) then
       isSameEG = .false.
    else if (eg1%isPlan) then
       isSameEG = pointsAlignes(eg1%x,eg1%y,eg1%dx,eg1%dy,eg2%x,eg2%y)  &
            .and. pointsAlignes(eg1%x,eg1%y,eg1%dx,eg1%dy,eg2%dx,eg2%dy)
    else
       isSameEG = isEqualConst(eg1%x,eg2%x) &
            .and. isEqualConst(eg1%y,eg2%y) &
            .and. isEqualConst(eg1%r,eg2%r)
    end if
    isSameEG = isSameEG .and. (eg1%limite==eg2%limite)
  end function isSameEG

  function isSameWayEG(eg1,eg2)
    type(t_elemGene),intent(in) :: eg1,eg2
    logical                     :: isSameWayEG

    if (eg1%isPlan .and. eg2%isPlan) then
       isSameWayEG = &
            (((eg1%dx-eg1%x)*(eg2%dx-eg2%x)+(eg1%dy-eg1%y)*(eg2%dy-eg2%y))>0)
    else
       isSameWayEG = .true.
    end if
  end function isSameWayEG

  subroutine prepareMCData(szSA,nbNode,szEG)
    integer,intent(in)    :: szSA,nbNode
    integer,intent(inout) :: szEG

    integer          :: i,j,volNbr,defautCl,lgMax,indFictive
    logical          :: found,toPut
    real             :: albedo
    character*4      :: strDCL
    double precision :: ptx,pty,rxx,ryy
    type(t_elemGene) :: eg
    integer,dimension(:),allocatable :: indEG,nbPoints
    logical,dimension(:),allocatable :: isGoodWay,withFictiveEG,inACircle
    double precision,dimension(:,:),allocatable :: xx,yy
    type(t_elemGene),dimension(:),allocatable   :: tmpTabEG

    !recuperation des donnees de conditions aux limites
    call calculDefaultCl(defautCl,albedo,strDCL)
    !preparation des donnees de remplissage des tableaux globaux
    allocate(indEG(szSA),isGoodWay(szSA),tmpTabEG(szSA+nbNode),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: prepareMCData(1) => allocation pb")
    do i = 1,szSA
       eg = createElemGene(tabSegArc(i),strDCL)
       found = .false.
       do j = 1,szEG
          found = isSameEG(eg,tmpTabEG(j))
          if (found) exit
       end do
       if (found) then
          indEG(i) = j
          isGoodWay(i) = isSameWayEG(eg,tmpTabEG(j))
       else
          szEG = szEG + 1
          tmpTabEG(szEG) = eg
          indEG(i) = szEG
          isGoodWay(i) = .true.
       end if
    end do
    !preparation des donnes permettant de savoir si un element fictif
    !englobant est necessaire
    allocate(inACircle(nbNode),withFictiveEG(nbNode),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: prepareMCData(2) => allocation pb")
    withFictiveEG(:nbNode) = .false. ; inACircle(:nbNode) = .false.
    !calcul des dimensions des pointeurs du tableau des volumes
    do i = 1,szSA
       !cote gauche
       volNbr = tabSegArc(i)%nodeg
       if (volNbr>0) then
          tabVolume(volNbr)%lg = tabVolume(volNbr)%lg + 1
          if (tabSegArc(i)%typ==tcer .and. .not. inACircle(volNbr)) then
             !le node est a l'interieur d'un cercle => on n'a pas besoin de
             !l'englober
             inACircle(volNbr) = .true.
             if (withFictiveEG(volNbr)) then
                !on avait commencer a englober fictivement => on laisse tomber
                tabVolume(volNbr)%lg = tabVolume(volNbr)%lg - 1
                withFictiveEG(volNbr) = .false.
             end if
          end if
       end if
       !cote droit
       volNbr = tabSegArc(i)%noded
       if (volNbr>0) then
          tabVolume(volNbr)%lg = tabVolume(volNbr)%lg + 1
          if (tabSegArc(i)%typ==tarc .and. .not. withFictiveEG(volNbr) &
              .and. .not. inACircle(volNbr)) then
             !le node est a l'exterieur d'un arc => on englobe si pas encore
             !fait
             tabVolume(volNbr)%lg = tabVolume(volNbr)%lg + 1
             withFictiveEG(volNbr) = .true.
          end if
       end if
    end do
    deallocate(inACircle)
    lgMax = maxval(tabVolume(:nbNode)%lg)
    !allocation des pointeurs du tableau des volumes
    do i = 1,nbNode
       allocate(tabVolume(i)%indElem(tabVolume(i)%lg),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: prepareMCData(3) => allocation pb")
       tabVolume(i)%indElem = 0
       allocate(tabVolume(i)%side(tabVolume(i)%lg),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: prepareMCData(4) => allocation pb")
       allocate(tabVolume(i)%typCL(tabVolume(i)%lg),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: prepareMCData(5) => allocation pb")
       tabVolume(i)%typCL = 0
       !remise a zero du compteur, pour gerer les doublons
       tabVolume(i)%lg = 0
    end do
    !remplissage des pointeurs, avec elimination des doublons
    allocate(xx(nbNode,lgMax),yy(nbNode,lgMax))
    allocate(nbPoints(nbNode),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: prepareMCData(6) => allocation pb")
    nbPoints(:nbNode) = 0 ; xx(:nbNode,:lgMax) = 0.d0 ; yy(:nbNode,:lgMax) = 0.d0
    do i = 1,szSA
       !cote gauche
       volNbr = tabSegArc(i)%nodeg
       if (volNbr>0) then
          toPut = ( count(tabVolume(volNbr)%indElem(:)==indEG(i)) == 0 )
          if (toPut) then
             tabVolume(volNbr)%mix = tabSegArc(i)%neutronicMixg
             tabVolume(volNbr)%lg = tabVolume(volNbr)%lg + 1
             tabVolume(volNbr)%indElem(tabVolume(volNbr)%lg) = indEG(i)
             tabVolume(volNbr)%typCL(tabVolume(volNbr)%lg) = &
                  tmpTabEG(indEG(i))%limite
             tabVolume(volNbr)%side(tabVolume(volNbr)%lg) = .not.isGoodWay(i)
          end if
          !preparation du calcul des coordonnees du cercle englobant
          if (withFictiveEG(volNbr).and.tabSegArc(i)%typ/=tcer) then
             call giveOrigine(tabSegArc(i),ptx,pty)
             nbPoints(volNbr) = nbPoints(volNbr) + 1
             xx(volNbr,nbPoints(volNbr)) = ptx
             yy(volNbr,nbPoints(volNbr)) = pty
          end if
       end if
       !cote droit
       volNbr = tabSegArc(i)%noded
       if (volNbr>0) then
          toPut = ( count(tabVolume(volNbr)%indElem(:)==indEG(i)) == 0 )
          if (toPut) then
             tabVolume(volNbr)%mix = tabSegArc(i)%neutronicMixd
             tabVolume(volNbr)%lg = tabVolume(volNbr)%lg + 1
             tabVolume(volNbr)%indElem(tabVolume(volNbr)%lg) = indEG(i)
             tabVolume(volNbr)%typCL(tabVolume(volNbr)%lg) = &
                  tmpTabEG(indEG(i))%limite
             tabVolume(volNbr)%side(tabVolume(volNbr)%lg) = isGoodWay(i)
          end if
          !preparation du calcul des coordonnees du cercle englobant
          if (withFictiveEG(volNbr).and.tabSegArc(i)%typ/=tcer) then
             call giveExtremite(tabSegArc(i),ptx,pty)
             nbPoints(volNbr) = nbPoints(volNbr) + 1
             xx(volNbr,nbPoints(volNbr)) = ptx
             yy(volNbr,nbPoints(volNbr)) = pty
          end if
       end if
    end do
    !calcul des cercles englobants
    do i = 1,nbNode
       if (.not.withFictiveEG(i)) cycle
       !creation du cercle
       eg%isPlan = .false.
       eg%x = sum(xx(i,:nbPoints(i)))/nbPoints(i)
       eg%y = sum(yy(i,:nbPoints(i)))/nbPoints(i)
       eg%dx = 0.d0 ; eg%dy = 0.d0 ; eg%r = 0.d0 ; eg%limite = Tri_Not
       do j = 1,nbPoints(i)
          rxx = xx(i,j)-eg%x ; ryy = yy(i,j)-eg%y
          eg%r = max(eg%r,longVect(rxx,ryy))
       end do
       !integration eventuelle dans le tableau des elements geometriques
       found = .false.
       do j = 1,szEG
          found = isSameEG(eg,tmpTabEG(j))
          if (found) exit
       end do
       if (found) then
          indFictive = j
       else
          szEG = szEG + 1
          tmpTabEG(szEG) = eg
          indFictive = szEG
       end if
       !ajout eventuel dans la liste des elements geometriques constitutifs
       !de la zone
       toPut = ( count(tabVolume(i)%indElem(:)==indFictive) == 0 )
       if (toPut) then
          tabVolume(i)%lg = tabVolume(i)%lg + 1
          tabVolume(i)%indElem(tabVolume(i)%lg) = indFictive
          tabVolume(i)%typCL(tabVolume(i)%lg) = Tri_Not
          tabVolume(i)%side(tabVolume(i)%lg) = .false.
       end if
    end do
    deallocate(nbPoints,xx,yy,indEG,isGoodWay,withFictiveEG)
    !recopie du tableau global des elements geometriques, et liberation
    !du temporaire
    allocate(tabEG(szEG),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: prepareMCData(7) => allocation pb")
    tabEG(1:szEG) = tmpTabEG(1:szEG)
    deallocate(tmpTabEG)
  end subroutine prepareMCData

  subroutine generateTripoliFile(fileTripoli,szSA,nbNode)
    integer,intent(in) :: fileTripoli,szSA,nbNode

    real    :: a,b,c,d,xx,yy,rr
    integer :: i,j,szEG,nbPlus,nbMoins,nbMix,nbCLTot,nbCL
    integer,dimension(:),allocatable :: sp,sm,listMix,volmix,clList

    szEG = 0
    call initMC(nbNode)
    call prepareMCData(szSA,nbNode,szEG)

    write(fileTripoli,'(/"GEOMETRY"/)')
    write(fileTripoli,'("//jeu de donnees geometriques TRIPOLI")')
    write(fileTripoli,'("//genere pour DRAGON par le module G2MC"/)')
    write(fileTripoli,'("TITRE geometrie DRAGON pour TRIPOLI"/)')
    !definition des surfaces
    write(fileTripoli,'("//definition des surfaces"/)')
    c = 0.
    do i = 1,szEG
       if (tabEG(i)%isPlan) then
          a = real(tabEG(i)%dy - tabEG(i)%y)
          b = real(tabEG(i)%x - tabEG(i)%dx)
          d = real(tabEG(i)%y*tabEG(i)%dx - tabEG(i)%x*tabEG(i)%dy)
          write(fileTripoli,*) " SURF ",i," PLAN ",a,b,c,d
       else
          xx = real(tabEG(i)%x)
          yy = real(tabEG(i)%y)
          rr = real(tabEG(i)%r)
          write(fileTripoli,*) " SURF ",i," CYLZ ",xx,yy,rr
       end if
    end do
    write(fileTripoli,*) " SURF ",szEG+1," PLANZ  0. //plan inferieur" 
    write(fileTripoli,*) " SURF ",szEG+2," PLANZ 30. //plan superieur"
    !definition des volumes
    write(fileTripoli,'(/"//definition des volumes"/)')
    do i = 1,nbNode
       nbPlus = count(tabVolume(i)%side(1:tabVolume(i)%lg))
       nbMoins = tabVolume(i)%lg - nbPlus
       allocate(sp(nbPlus),sm(nbMoins),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: generateTripoliFile(1) => allocation pb")
       sp = pack(tabVolume(i)%indElem(1:tabVolume(i)%lg), &
                 tabVolume(i)%side(1:tabVolume(i)%lg))
       sm = pack(tabVolume(i)%indElem(1:tabVolume(i)%lg), &
                 .not.tabVolume(i)%side(1:tabVolume(i)%lg))
       write(fileTripoli,*) " VOLU"
       write(fileTripoli,*) "  ",i
       write(fileTripoli,*) "   EQUA"
       write(fileTripoli,*) "     PLUS",nbPlus+1
       write(fileTripoli,*) "      ",szEG+1,sp(:) !plan inf en premier
       write(fileTripoli,*) "     MOINS",nbMoins+1
       write(fileTripoli,*) "      ",sm(:),szEG+2 !plan sup en dernier
       write(fileTripoli,*) " FINV"
       write(fileTripoli,*)
       deallocate(sp,sm)
    end do
    write(fileTripoli,'("FINGEOM")')
    
    !conditions aux limites
    write(fileTripoli,'(/"//conditions aux limites")')
    write(fileTripoli,'("LIMIT")')
    !calcul du nombre total de lignes de conditions limites particulieres
    nbCLTot = 2*nbNode
    do i = 1,nbNode
       nbCLTot = nbCLTot &
            + count(tabVolume(i)%typCL(1:tabVolume(i)%lg)==Tri_Refl) &
            + count(tabVolume(i)%typCL(1:tabVolume(i)%lg)==Tri_Trans) &
            + count(tabVolume(i)%typCL(1:tabVolume(i)%lg)==Tri_Cos)
    end do
    write(fileTripoli,*) nbCLTot
    !affichage des conditions limites
    do i = 1,nbNode
       !Reflexions
       write(fileTripoli,*) "  ",i," REFLECTION ",szEG+1
       write(fileTripoli,*) "  ",i," REFLECTION ",szEG+2
       nbCL = count(tabVolume(i)%typCL(1:tabVolume(i)%lg)==Tri_Refl)
       if (nbCL/=0) then
          allocate(clList(nbCL),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: generateTripoliFile(2) => allocation pb")
          clList = pack(tabVolume(i)%indElem(1:tabVolume(i)%lg), &
                        tabVolume(i)%typCL(1:tabVolume(i)%lg)==Tri_Refl)
          do j = 1,nbCL
             write(fileTripoli,*) "  ",i," REFLECTION ",clList(j)
          end do
          deallocate(clList)
       end if
       !Translations
       nbCL = count(tabVolume(i)%typCL(1:tabVolume(i)%lg)==Tri_Trans)
       if (nbCL/=0) then
          allocate(clList(nbCL),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: generateTripoliFile(3) => allocation pb")
          clList = pack(tabVolume(i)%indElem(1:tabVolume(i)%lg), &
                        tabVolume(i)%typCL(1:tabVolume(i)%lg)==Tri_Trans)
          do j = 1,nbCL
             write(fileTripoli,*) "  ",i," TRANSLATION ",clList(j)
          end do
          deallocate(clList)
       end if
       !Cosinus
       nbCL = count(tabVolume(i)%typCL(1:tabVolume(i)%lg)==Tri_Cos)
       if (nbCL/=0) then
          allocate(clList(nbCL),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: generateTripoliFile(4) => allocation pb")
          clList = pack(tabVolume(i)%indElem(1:tabVolume(i)%lg), &
                        tabVolume(i)%typCL(1:tabVolume(i)%lg)==Tri_Cos)
          do j = 1,nbCL
             write(fileTripoli,*) "  ",i," COSINUS ",clList(j)
          end do
          deallocate(clList)
       end if
    end do
    write(fileTripoli,'("FIN_LIMIT"/)')

    ! insertion des donnees neutroniques des materiaux
    ! le milieu ne correspond au materiaux nomme "MIX_n"
    write(fileTripoli,'("//insertion des donnees neutroniques")')
    write(fileTripoli,'("//des matreriaux")')
    write(fileTripoli,'("FILE")')
    write(fileTripoli,'("   mix_data"/)') !le nom est fixe, mais pourra etre
                                          !change si besoin  ######

    !definition du milieu par volume
    write(fileTripoli,'("//definition du milieu par volume")')
    write(fileTripoli,'("GEOMCOMP")')
    allocate(listMix(nbNode),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: generateTripoliFile(5) => allocation pb")
    listMix(:nbNode) = 0
    nbMix = 0
    do i = 1,nbNode
       if (count(listMix(:nbNode)==tabVolume(i)%mix)==0) then
          nbMix = nbMix + 1
          listMix(nbMix) = tabVolume(i)%mix
       end if
    end do
    do i = 1,nbMix
       allocate(volMix(count(tabVolume(:)%mix==listMix(i))),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: generateTripoliFile(6) => allocation pb")
       volMix = pack((/(j,j=1,nbNode)/),tabVolume(:)%mix==listMix(i))
       write(fileTripoli,*) " MIX_"//trim(i2s(listMix(i))),size(volMix)
       write(fileTripoli,*) "   ",volMix
       deallocate(volMix)
    end do
    deallocate(listMix)
    write(fileTripoli,'("FIN_GEOMCOMP"/)')

    call destroyMC()
  end subroutine generateTripoliFile

  subroutine generateMCNPFile(fileMCNP,szSA,nbNode)
    integer,intent(in) :: fileMCNP,szSA,nbNode

    integer :: i,j,k,szEG,ind,transInd,nbEGCL
    real    :: a,b,c,d,xx,yy,rr
    character(len=1)    :: sgn
    character(len=5)    :: nbr,mix
    character(len=6)    :: other,signedNbr
    character(len=12)   :: text12
    character(len=1000) :: buffer
    integer,dimension(:),allocatable :: tabEGCL

    szEG = 0
    call initMC(nbNode)
    call prepareMCData(szSA,nbNode,szEG)
    
    !title card
    write(fileMCNP,'("Jeu de donnees geometriques MCNP genere pour DRAGON &
         &par le module G2MC")')
    !message en commentaire
    write(fileMCNP,'("C Attention ! Les densites des milieux n''etant pas &
         &disponibles dans les")')
    write(fileMCNP,'("C jeux de donnees geometriques DRAGON, il faut passer &
         &une moulinette pour")')
    write(fileMCNP,'("C remplacer les ""DENSITE_I"" presents dans ce fichier, &
         &par leur valeur."/)')
    !cell card
    write(fileMCNP,'("C   cell card")')
    do i = 1,nbNode
       text12 = i2s(i)
       nbr = text12(:5)
       text12 = i2s(tabVolume(i)%mix)
       mix = text12(:5)
       buffer =nbr//" "//trim(mix)//" DENSITE_"//trim(mix)
       ind = len_trim(buffer)+1
       write(buffer(ind:),*) merge(tabVolume(i)%indElem(1:tabVolume(i)%lg), &
            -tabVolume(i)%indElem(1:tabVolume(i)%lg), &
            tabVolume(i)%side(1:tabVolume(i)%lg)) , &
            szEG+1 , -(szEG+2)
       call putOn80col(fileMCNP,buffer)
    end do
    text12=i2s(nbNode+1)
    nbr = text12(:5)
    buffer = nbr//" 0"
    ind = 9
    if (geomTyp==HexTyp.or.geomTyp==TriaTyp) then
       !on a des risques de geometrie non convexe
       ! => on prend les complementaires
       do i = 1,nbNode
          text12 = i2s(i)
          nbr = text12(:5)
          buffer(ind:) = "#"//nbr
          ind = ind + len_trim(nbr) + 2
       end do
    else
       !la geometrie est convexe
       ! on prend l'union des complementaires des cotes a condition limite
       nbEGCL = count(tabEG(:)%limite/=Tri_Not)
       allocate(tabEGCL(nbEGCL),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: generateTripoliFile(7) => allocation pb")
       tabEGCL = pack((/(i,i=1,szEG)/),tabEG(:)%limite/=Tri_Not)
       do i = 1,nbEGCL
          do j = 1,nbNode
             do k = 1,tabVolume(j)%lg
                if (tabVolume(j)%indElem(k)==tabEGCL(i)) then
                   if (tabVolume(j)%side(k)) tabEGCL(i) = -tabEGCL(i)
                   exit
                end if
             end do
          end do
       end do
       do i = 1,nbEGCL
          text12 = i2s(tabEGCL(i))
          signedNbr = text12(:6)
          buffer(ind:) = trim(signedNbr)//" :"
          ind = ind + len_trim(signedNbr) + 3
       end do
       buffer(ind:) = trim(i2s(-(szEG+1)))//" : "//i2s(szEG+2)
       deallocate(tabEGCL)
    end if
    call putOn80col(fileMCNP,buffer)

    !surface card
    write(fileMCNP,'(/"C   surface card")')
    c = 0.
    do i = 1,szEG
       text12 = i2s(i)
       nbr = text12(:5)
       sgn = ' '
       other = ' '
       select case(tabEG(i)%limite)
       case(Tri_Refl)
          sgn = '*'
       case(Tri_Trans)
          transInd = findParalleleWithTrans(i,szEG)
          text12=i2s(transInd)
          if (isSameWayEG(tabEG(i),tabEG(transInd))) then
             other = '-'//text12(:6)
          else
             other = text12(:6)
          end if
       case(Tri_Cos)
          sgn = '+'
       end select
       if (tabEG(i)%isPlan) then
          a = real(tabEG(i)%dy - tabEG(i)%y)
          b = real(tabEG(i)%x - tabEG(i)%dx)
          d = real(tabEG(i)%y*tabEG(i)%dx - tabEG(i)%x*tabEG(i)%dy)
          buffer = adjustl(sgn//nbr//trim(" "//other)//" P")
          ind = len_trim(buffer) + 1
          write(buffer(ind:),*) a,b,c,d
       else
          xx = real(tabEG(i)%x)
          yy = real(tabEG(i)%y)
          rr = real(tabEG(i)%r)
          buffer = adjustl(sgn//nbr//" C/Z")
          ind = len_trim(buffer) + 1
          write(buffer(ind:),*) xx,yy,rr
       end if
       call putOn80col(fileMCNP,buffer)
    end do
    text12 = i2s(szEG+1)
    nbr = text12(:5)
    buffer = "*"//nbr//" PZ  0. $plan inferieur"
    call putOn80col(fileMCNP,buffer)
    text12 = i2s(szEG+2)
    nbr = text12(:5)
    buffer = "*"//nbr//" PZ 30. $plan superieur"
    call putOn80col(fileMCNP,buffer)

    !data card

    call destroyMC()
  end subroutine generateMCNPFile

  subroutine generateSerpentFile(fileSerpent,szSA,nbNode)
    integer,intent(in) :: fileSerpent,szSA,nbNode

    integer :: i,j,k,szEG,ind,transInd,nbEGCL,nb2,isurf,jjj 
    real    :: a,b,c,d,xx,yy,rr,cuboid(6)
    character(len=1)    :: sgn
    character(len=5)    :: nbr,mix
    character(len=6)    :: other
    character(len=12)   :: text12
    character(len=1000) :: buffer
    integer,dimension(:),allocatable :: tabEGCL,tmpnod,tmpnod2
    logical :: lcuboid

    szEG = 0
    call initMC(nbNode)
    call prepareMCData(szSA,nbNode,szEG)

    !title card
    write(fileSerpent,'("% Jeu de donnees geometriques Serpent genere par le module&
         & G2MC de DRAGON")')
    !cell card
    nbEGCL = count(tabEG(:)%limite/=Tri_Not)
    allocate(tabEGCL(nbEGCL),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: generateSerpentFile(1) => allocation pb")
    tabEGCL = pack((/(i,i=1,szEG)/),tabEG(:)%limite/=Tri_Not)
    do i = 1,nbEGCL
       do j = 1,nbNode
          do k = 1,tabVolume(j)%lg
             if (tabVolume(j)%indElem(k)==tabEGCL(i)) then
                if (tabVolume(j)%side(k)) tabEGCL(i) = -tabEGCL(i)
                exit
             end if
          end do
       end do
    end do
    write(fileSerpent,'("%   cell card")')
    nb2=0
    cuboid=(/ 0.0, 0.0, 0.0, 0.0, -1.0E10, 1.0E10 /)
    do i = 1,nbNode
       allocate(tmpnod(tabVolume(i)%lg),tmpnod2(tabVolume(i)%lg),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: generateSerpentFile(2) => allocation pb")
       text12 = i2s(i)
       nbr = text12(:5)
       text12 = i2s(tabVolume(i)%mix)
       mix = text12(:5)
       buffer ="cell "//nbr//" 0 MIX_"//trim(mix)
       ind = len_trim(buffer)+1
       tmpnod(:tabVolume(i)%lg)=merge(-tabVolume(i)%indElem(1:tabVolume(i)%lg), &
            tabVolume(i)%indElem(1:tabVolume(i)%lg), &
            tabVolume(i)%side(1:tabVolume(i)%lg))
       jjj=0
       lcuboid=.true.
       do j = 1,tabVolume(i)%lg
         do k =1,nbEGCL
           if(tmpnod(j) == tabEGCL(k)) then
             isurf=abs(tabEGCL(k))
             if (.not.tabEG(isurf)%isPlan) call XABORT('g2s_generatingMC: unsupported BC')
             a = real(tabEG(isurf)%dy - tabEG(isurf)%y)
             b = real(tabEG(isurf)%x - tabEG(isurf)%dx)
             d = real(tabEG(isurf)%y*tabEG(isurf)%dx - tabEG(isurf)%x*tabEG(isurf)%dy)
             if ((a == 0.0).and.(tabEGCL(k) < 0)) then
               cuboid(3)=-d/b
             else if ((a == 0.0).and.(tabEGCL(k) > 0)) then
               cuboid(4)=-d/b
             else if ((b == 0.0).and.(tabEGCL(k) < 0)) then
               cuboid(1)=-d/a
             else if ((b == 0.0).and.(tabEGCL(k) > 0)) then
               cuboid(2)=-d/a
             endif
             if(lcuboid) then
               jjj=jjj+1
               tmpnod2(jjj)=-szEG-1
               lcuboid=.false.
             endif
             go to 10
           endif
         end do
         jjj=jjj+1
         tmpnod2(jjj)=tmpnod(j)
         if(tabEG(abs(tmpnod(j)))%isPlan) tmpnod2(jjj)=-tmpnod(j)
    10   continue
       end do
       if(jjj.gt.0) then
         write(buffer(ind:),*) tmpnod2(:jjj)
         write(fileSerpent,'(a)') trim(buffer)
       endif
       deallocate(tmpnod2,tmpnod)
    end do
    if (nbEGCL > 0) then
       text12 = i2s(nbNode+1)
       nbr = text12(:5)
       text12 = i2s(szEG+1)
       buffer ="cell "//nbr//" 0 outside         "//trim(text12)
       write(fileSerpent,'(a)') trim(buffer)
    endif

    !surface card
    write(fileSerpent,'(/"%   surface card")')
    c = 0.
    do i = 1,szEG
       do k = 1,nbEGCL
         if (i == abs(tabEGCL(k))) go to 20
       end do
       text12 = i2s(i)
       nbr = text12(:5)
       sgn = ' '
       other = ' '
       select case(tabEG(i)%limite)
       case(Tri_Refl)
          sgn = '*'
       case(Tri_Trans)
          transInd = findParalleleWithTrans(i,szEG)
          text12=i2s(transInd)
          if (isSameWayEG(tabEG(i),tabEG(transInd))) then
             other = '-'//text12(:6)
          else
             other = text12(:6)
          end if
       case(Tri_Cos)
          sgn = '+'
       end select
       if (tabEG(i)%isPlan) then
         a = real(tabEG(i)%dy - tabEG(i)%y)
         b = real(tabEG(i)%x - tabEG(i)%dx)
         d = -real(tabEG(i)%y*tabEG(i)%dx - tabEG(i)%x*tabEG(i)%dy)
         buffer = adjustl("surf "//sgn//nbr//trim(" "//other)//" plane")
         ind = len_trim(buffer) + 1
         write(buffer(ind:),*) a,b,c,d
       else
         xx = real(tabEG(i)%x)
         yy = real(tabEG(i)%y)
         rr = real(tabEG(i)%r)
         buffer = adjustl("surf "//sgn//nbr//" cyl")
         ind = len_trim(buffer) + 1
         write(buffer(ind:),*) xx,yy,rr
       end if
       write(fileSerpent,'(a)') trim(buffer)
   20  continue
    end do
    if (nbEGCL > 0) then
       text12 = i2s(szEG+1)
       nbr = text12(:5)
       buffer = adjustl("surf "//sgn//nbr//" cuboid ")
       ind = len_trim(buffer) + 1
       write(buffer(ind:),*) (cuboid(i),i=1,6)
       write(fileSerpent,'(a)') trim(buffer)
    endif
    deallocate(tabEGCL)
    !data card

    call destroyMC()
  end subroutine generateSerpentFile

  subroutine putOn80col(outFile,strIn)
    integer,intent(in)          :: outFile
    character(len=*),intent(in) :: strIn

    integer           :: lg,lastBlank
    character(len=15) :: frm
    character(len=80) :: truncStr
    character(len=len(strIn)) :: str

    str = strIn
    do
       lg = len_trim(str)
       if (lg <= 80) then
          frm = "(a"//i2s(lg)//")"
          write(outFile,frm) str
          return
       end if
       truncStr = str(1:78)
       do lastBlank = 78,1,-1
          if (truncStr(lastBlank:lastBlank)==' ') exit
       end do
       truncStr = trim(str(1:lastBlank-1))//' &'
       frm = "(a"//i2s(len_trim(truncStr))//")"
       write(outFile,frm) truncStr
       str = '     '//adjustl(str(lastBlank:))
    end do
  end subroutine putOn80col

  function findParalleleWithTrans(ind,szEG)
    integer,intent(in) :: ind,szEG
    integer            :: findParalleleWithTrans

    integer :: i

    findParalleleWithTrans = 0
    do i = 1,szEG
       if (i==ind) cycle
       if (.not.tabEG(i)%isPlan) cycle
       if (tabEG(i)%limite/=Tri_Trans) cycle
       if (.not.estColi(tabEG(i)%dx-tabEG(i)%x,tabEG(i)%dy-tabEG(i)%y, &
            tabEG(ind)%dx-tabEG(ind)%x,tabEG(ind)%dy-tabEG(ind)%y)) cycle
       findParalleleWithTrans = i
       return
    end do
    call XABORT("G2MC: error, no parallel side found for translation &
         &boundary condition")
  end function findParalleleWithTrans

end module monteCarlo
