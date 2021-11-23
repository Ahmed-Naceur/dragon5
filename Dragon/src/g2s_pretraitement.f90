!
!-----------------------------------------------------------------------
!
!Purpose:
! Preprocessing of geometric data recovered from a geometry LCM data
! structure.
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
! Ce fichier regroupe toutes les actions effectuees en amont du code, pour
! verifier et completer le jeux de donnees geometriques. Il utilise de maniere
! intensive l'api de PyLCM.
! Les fonctions de ce fichier ont ete initialement developpees en python, puis
! traduite en fortran90. C'est pourquoi la forme des algorithmes peut etre
! parfois deroutante... 
! \\\\
! Le module pretraitement definit quelques constantes, dont en particulier la
! precision admise en entree du code (pour la coherence des donnees).
! De plus, les matrices "composeRotRec" et "composeRotTri" permettent
! d'effectuer la composition de deux transformations definies selon le mode
! Dragon.
! \\\\
! Le seul pretraitement qui soit reellement delicat, est celui qui correspond 
! a une geometrie de type gigogne rectangulaire. Il est effectue par la routine
! "creerEtChargerNewDataRec", et est documente dans le source de cette routine.
!
!-----------------------------------------------------------------------
!
module pretraitement
    use boundcond
    use cast
    use celluleBase
    use cellulePlaced
    use constType
    use GANLIB

    implicit none

  real,parameter :: geomPrec = 1.e-5
  real,parameter :: sqrt_3f  = 1.732050807568

  !composeRotXXX(t1,t2)=t2ot1
  integer,dimension(8,8),parameter :: composeRotRec = &
       & reshape((/1,2,3,4,5,6,7,8, &
       &           2,3,4,1,8,5,6,7, &
       &           3,4,1,2,7,8,5,6, &
       &           4,1,2,3,6,7,8,5, &
       &           5,6,7,8,1,2,3,4, &
       &           6,7,8,5,4,1,2,3, &
       &           7,8,5,6,3,4,1,2, &
       &           8,5,6,7,2,3,4,1/) , (/8,8/))

  integer,dimension(12,12),parameter :: composeRotTri = &
       & reshape((/ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12, &
       &            2, 3, 4, 5, 6, 1, 8, 9,10,11,12, 7, &
       &            3, 4, 5, 6, 1, 2, 9,10,11,12, 7, 8, &
       &            4, 5, 6, 1, 2, 3,10,11,12, 7, 8, 9, &
       &            5, 6, 1, 2, 3, 4,11,12, 7, 8, 9,10, &
       &            6, 1, 2, 3, 4, 5,12, 7, 8, 9,10,11, &
       &            7,12,11,10, 9, 8, 1, 6, 5, 4, 3, 2, &
       &            8, 7,12,11,10, 9, 2, 1, 6, 5, 4, 3, &
       &            9, 8, 7,12,11,10, 3, 2, 1, 6, 5, 4, &
       &           10, 9, 8, 7,12,11, 4, 3, 2, 1, 6, 5, &
       &           11,10, 9, 8, 7,12, 5, 4, 3, 2, 1, 6, &
       &           12,11,10, 9, 8, 7, 6, 5, 4, 3, 2, 1/) , (/12,12/))

contains
  subroutine prepareData(geoIp,sizeB,sizeP,lgMaxGig)
    type(c_ptr),intent(in):: geoIp
    integer,intent(inout) :: sizeB,sizeP,lgMaxGig
    
    integer,dimension(40) :: st
    type(c_ptr)           :: ip
    integer               :: nbNewMix,i

    !premier completion de la structure geometrique pour transformer
    !les constructions melangeant MIX et CELL, en des constructions
    !uniformes de ce point de vue

    nbNewMix = 0
    call separateMixAndCell(geoIp,nbNewMix)

    ip = geoIp
    !creation des cellules de bases
    call buildCellsBase(ip,sizeB,'/           ')
    ip = geoIp

    call LCMGET(ip,'STATE-VECTOR',st)
    select case(st(1))
    case(G_Car2d,G_Carcel)
       geomTyp = RecTyp
       call creerEtChargerNewDataRec(ip,sizeB,sizeP)
    case(G_Hex,G_Hexcel)
       geomTyp = HexTyp
       call creerEtChargerNewDataHex(ip,sizeB,sizeP)
    case(G_Tri)
       geomTyp = TriaTyp
       call creerEtChargerNewDataTri(ip,sizeB,sizeP)
    case(G_Tube)
       geomTyp = TubeTyp
       if (sizeB/=1) call XABORT("G2S: more than one cellule in a cylindrical&
            & geometry")
       sizeP = 1
       tabCellulePlaced(1)%indice = 1
       tabCellulePlaced(1)%xcenter = 0.d0
       tabCellulePlaced(1)%ycenter = 0.d0
       tabCellulePlaced(1)%turn = 1
       allocate (tabCellulePlaced(1)%gig(1),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: prepareData(2) => allocation pb")
       tabCellulePlaced(1)%gig = 1
       allocate (tabCellulePlaced(1)%mrg(1),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: prepareData(3) => allocation pb")
       tabCellulePlaced(1)%mrg = 1
    case default
       call XABORT("G2S: Type of geometry not supported")
    end select

    lgMaxGig = 0
    do i = 1,sizeP
       lgMaxGig = max(lgMaxGig,size(tabCellulePlaced(i)%gig))
    end do
  end subroutine prepareData

  recursive subroutine separateMixAndCell(ipIn,nbNewMix)
    type(c_ptr),intent(in):: ipIn
    integer,intent(inout) :: nbNewMix

    type(c_ptr)  :: ip,jp
    integer      :: i,lgC,lg,ty,indOfCellToAdd,mixToPut,nbCellToAdd,lgC2
    character*12 :: newCellName,text12
    integer,dimension(40) :: sv,newSv
    integer,dimension(:),allocatable      :: mix
    character*12,dimension(:),allocatable :: cell

    ip = ipIn
    call LCMLEN(ip,'CELL        ',lgC,ty)
    lgC = lgC/3
    if (lgC/=0) then !il y a des sous-cellules
       allocate(cell(lgC),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: separateMixAndCell(1) => allocation pb")
       call LCMGTC(ip,'CELL        ',12,lgC,cell)
       do i = 1,lgC
          !traitement recursif sur les sous cellules de la cellule courante
          call LCMLEN(ip,cell(i),lg,ty)
          if (ty==0) then !la sous cellule est bien definie ici
             jp=LCMDID(ip,cell(i))
             call separateMixAndCell(jp,nbNewMix)
          end if
       end do
       deallocate(cell)
       !test sur la mixite de la cellule courante
       call LCMLEN(ip,'MIX         ',lg,ty)
       allocate(mix(lg),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: prepareData(2) => allocation pb")
       call LCMGET(ip,'MIX         ',mix)
       nbCellToAdd = count(mix(:lg)>0)
       if (nbCellToAdd/=0) then !la cellule courante est mixte MIX/CELL
          call LCMGET(ip,'STATE-VECTOR',sv)
          newSv(:) = 0
          newSv(1) = sv(1)
          newSv(2) = 0
          newSv(3) = 1
          newSv(4) = 1
          newSv(5) = 0
          newSv(6) = 1
          newSv(7) = 1
          allocate(cell(lg),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: prepareData(3) => allocation pb")
          cell(:lg)=' '
          call LCMLEN(ip,'CELL        ',lgC2,ty)
          lgC2 = lgC2/3
          call LCMGTC(ip,'CELL        ',12,lgC2,cell)
          indOfCellToAdd = lgC
          do i = 1,lg
             if (mix(i)>0) then
                !creation d'une nouvelle cellule
                nbNewMix = nbNewMix+1
                indOfCellToAdd = indOfCellToAdd+1
                mixToPut = mix(i)
                mix(i) = -indOfCellToAdd
                text12 = i2s(nbNewMix)
                newCellName = 'Nmix'//text12(:8)
                cell(indOfCellToAdd) = newCellName
                jp=LCMDID(ip,newCellName)
                call LCMPUT(jp,'MIX         ',1,1,mixToPut)
                newSv(7) = mixToPut
                call LCMPUT(jp,'STATE-VECTOR',40,1,newSv)
                call LCMPUT(jp,'NCODE       ',6,1,(/0,0,0,0,0,0/))
                call LCMPUT(jp,'ICODE       ',6,1,(/0,0,0,0,0,0/))
                text12='L_GEOM'
                call LCMPTC(jp,'SIGNATURE',12,1,text12)
                call LCMPUT(jp,'ZCODE       ',6,2,(/0.,0.,0.,0.,0.,0./))
             end if
          end do
          !modification du STATE-VECTOR, du MIX et du CELL de la
          !cellule courante
          sv(9) = sv(9)+nbCellToAdd
          call LCMPUT(ip,'STATE-VECTOR',40,1,sv)
          call LCMPUT(ip,'MIX         ',lg,1,mix)
          call LCMPTC(ip,'CELL        ',12,lg,cell)
          deallocate(cell)
       end if
       deallocate(mix)
    end if
  end subroutine separateMixAndCell

  subroutine creerEtChargerNewDataRec(geoIp,szB,szP)
    type(c_ptr),intent(in):: geoIp
    integer,intent(inout) :: szB,szP

    !effectue le pretraitement initialement ecrit en python.
    !celui-ci utilisait intensivement les dictionnaires python.
    !Ceux-ci seront donc remplaces par des objets LCM
    type(c_ptr) :: ip,sidesIp,lCentreIp,lMinMaxIp,centreCalculesIp,lcIp,dicoIp
    integer     :: szLc

    ip = geoIp
    !modification de la geometrie dans le cas diagonal, par ajout des
    !cellules manquantes dans la gigogne exterieure (ceci pour permettre
    !la mise en oeuvre du traitement general aux geometries rectangulaires)
    call traiteConditionDiagonale(ip)
    !creation de la liste des carres
    call LCMOP(lcIp,'lc          ',0,1,0)
    call creeListeCellules(lcIp,ip,szLc)
    !creation du dictionnaire de verification des donnees geometriques
    call LCMOP(dicoIp,'dico        ',0,1,0)
    call creeDico(dicoIp,lcIp,szLc)

    !creation de la liste des longueurs des cotes et resolution du systeme
    call LCMOP(sidesIp,'sides       ',0,1,0)
    call resoudDico(dicoIp,sidesIp)
    !cree et prepare la liste des centres des cellules de base,
    !avec leur gigogne d'origine
    !et destruction du dictionnaire
    call LCMOP(lCentreIp,'centres     ',0,1,0)
    call prepareCentres(lCentreIp,lcIp,sidesIp,dicoIp)
    call LCMCL(dicoIp,2)
    !preparation des longeurs de la gigogne exterieure, et des coordonnees
    !min et max de ses sous cellules
    call LCMOP(lMinMaxIp,'minAndMaxXY ',0,1,0)
    call prepareMinMax(lMinMaxIp,lCentreIp)
    !resolution du systeme donnant les coordonnees des centres des cellules
    !de base et stockage dans une nouvelle structure 
    call LCMOP(centreCalculesIp,'resCentres  ',0,1,0)
    call calculeCentre(centreCalculesIp,lCentreIp,0)
    !compilation des resultats et ajout dans la geometrie des nouvelles donnees
    call compileResutats(ip,centreCalculesIp,sidesIp,lMinMaxIp)
    call LCMCL(lcIp,2)
    call LCMCL(centreCalculesIp,2)
    call LCMCL(sidesIp,2)
    call LCMCL(lMinMaxIp,2)
    !utilisation des donnees (on accroche le traitement classique developpe
    !pour python)
    call chargerNewData(geoIp,szB,szP)
  end subroutine creerEtChargerNewDataRec

  subroutine traiteConditionDiagonale(geoIp)
    type(c_ptr),intent(in) :: geoIp

    integer :: i,j,k,pos,n,lgAv,lgAp,lg,ty
    integer,dimension(6)  :: ncode
    integer,dimension(40) :: sv
    integer,dimension(:),allocatable :: mixAv,mixAp,turnAv,turnAp
    integer,dimension(:),allocatable :: mergeAv,mergeAp

    call LCMGET(geoIp,'NCODE       ',ncode)
    if ((ncode(1)/=B_Diag) .and. (ncode(2)/=B_Diag)) &
         return !pas de condition DIAG
    call LCMGET(geoIp,'STATE-VECTOR',sv)
    if (sv(8)==0) return !pas de sous cellules
    call LCMLEN(geoIp,'MIX         ',lgAv,ty)
    n = (nint(sqrt(1.+8.*lgAv)) - 1) / 2
    lgAp = n*n

    allocate(mixAv(lgAv),turnAv(lgAv),mergeAv(lgAv),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: traiteConditionDiagonale(1) => allocation pb")
    allocate(mixAp(lgAp),turnAp(lgAp),mergeAp(lgAp),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: traiteConditionDiagonale(2) => allocation pb")
    call LCMGET(geoIp,'MIX         ',mixAv)
    call LCMLEN(geoIp,'TURN        ',lg,ty)
    if (lg == 0) then
       turnAv(:lgAv) = 1
    else
       call LCMGET(geoIp,'TURN        ',turnAv)
    end if
    call LCMLEN(geoIp,'MERGE       ',lg,ty)
    if (lg == 0) then
       mergeAV(:lgAv) = (/(i,i=1,lgAv)/)
    else
       call LCMGET(geoIp,'MERGE       ',mergeAV)
    end if

    pos = 0
    if ((ncode(1)==B_Diag) .and. (ncode(4)==B_Diag)) then      !triangle inf
       do i = 1,n
          do j = 1,i-1
             k = i + ((j-1)*(2*n-j))/2
             pos = pos + 1
             mixAp(pos)   = mixAv(k)
             mergeAP(pos) = mergeAV(k)
             turnAp(pos)  = composeRotRec(turnAv(k),6)
          end do
          do j = i,n
             k = j + ((i-1)*(2*n-i))/2
             pos = pos + 1
             mixAp(pos)   = mixAv(k)
             mergeAP(pos) = mergeAV(k)
             turnAp(pos)  = turnAv(k)
          end do
       end do
    else if ((ncode(2)==B_Diag) .and. (ncode(3)==B_Diag)) then !triangle sup
       do i = 1,n
          do j = 1,i
             k = j + (i*(i-1))/2
             pos = pos + 1
             mixAp(pos)   = mixAv(k)
             mergeAP(pos) = mergeAV(k)
             turnAp(pos)  = turnAv(k)
          end do
          do j = i+1,n
             k = i + (j*(j-1))/2
             pos = pos + 1
             mixAp(pos)   = mixAv(k)
             mergeAP(pos) = mergeAV(k)
             turnAp(pos)  = composeRotRec(turnAv(k),6)
          end do
       end do
    else
       call XABORT("G2S: internal error in routine traiteConditionDiagonale")
    end if

    call LCMPUT(geoIp,'MIX         ',pos,1,mixAp)
    call LCMPUT(geoIp,'TURN        ',pos,1,turnAp)
    call LCMPUT(geoIp,'MERGE       ',pos,1,mergeAp)
    deallocate(mixAv,turnAv,mixAp,turnAp,mergeAv,mergeAp)
  end subroutine traiteConditionDiagonale

  subroutine creeListeCellules(lcIp,geoIp,nbCarre)
    type(c_ptr),intent(in) :: lcIp,geoIp
    integer,intent(out) :: nbCarre

    integer  :: nbAV
    
    nbCarre = 0
    call rempliListeCellules(lcIp,geoIp,nbCarre)
    do
       nbAV = nbCarre
       call rempliListeCellules(lcIp,geoIp,nbCarre)
       if (nbCarre==nbAV) exit
    end do
  end subroutine creeListeCellules

  subroutine rempliListeCellules(lcIp,geoIp,nbCarre)
    type(c_ptr),intent(in):: lcIp,geoIp
    integer,intent(inout) :: nbCarre

    character*12 :: subName
    type(c_ptr)  :: subIp,jplist
    integer      :: i,j,long,typ,nbCell,long2
    character*12,dimension(:),allocatable :: dataCh12
    
    if (nbCarre==0) then
       nbCarre = 1
       jplist=LCMDID(lcIp,'root        ')
       call LCMEQU(geoIp,jplist)
    end if
    subName = ' '
    do i = 1,nbCarre
       call LCMNXT(lcIp,subName)
       subIp=LCMGID(lcIp,subName)
       call LCMLEN(subIp,'CELL        ',long,typ)
       if (long/=0) then
          nbCell = long/3
          allocate(dataCh12(nbCell),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: rempliListeCellules => allocation pb")
          call LCMGTC(subIp,'CELL        ',12,nbCell,dataCh12)
          do j = 1,nbCell
             call LCMLEN(subIp,dataCh12(j),long,typ)
             call LCMLEN(lcIp,dataCh12(j),long2,typ)
             if (long/=0.and.long2==0) then
                call LCMSIX(subIp,dataCh12(j),1)
                nbCarre = nbCarre + 1
                jplist=LCMDID(lcIp,dataCh12(j))
                call LCMEQU(subIp,jplist)
                call LCMSIX(subIp,dataCh12(j),2)
             end if
          end do
          deallocate(dataCh12)
       end if
    end do
  end subroutine rempliListeCellules

  subroutine creeDico(dicoIp,lcIp,szLc)
    type(c_ptr),intent(in) :: dicoIp,lcIp
    integer,intent(in) :: szLc
    
    character*12 :: carreName,sideName,tmpStr,arrayName,nbStr
    type(c_ptr)  :: cIp,dIp,dJp
    integer      :: i,j,k,long,typ,lenMix,ind,lignNbr,lenCell
    real         :: side
    integer,dimension(40)                 :: sv
    integer,dimension(:),allocatable      :: mix,turn,turnStr,merg,mergStr
    character*12,dimension(:),allocatable :: cell,tabCh
    
    !initialisation
    dIp = dicoIp
    carreName = ' '
    do i = 1,szLc
       call LCMNXT(lcIp,carreName)
       cIp=LCMGID(lcIp,carreName)
       sideName = 'x' // carreName(1:11)
       do j = 1,2
          call getSquareSide(cIp,j,side)
          dJp=LCMDID(dIp,sideName)
          call LCMPUT(dJp,'value',1,2,side)
          sideName(1:1)='y'
       end do
    end do
    !remplissage
    dIp = dicoIp
    carreName = ' '
    lignNbr = 0
    do k = 1,szLc
       !tavail sur la kieme cellule
       call LCMNXT(lcIp,carreName)
       cIp=LCMGID(lcIp,carreName)
       call LCMGET(cIp,'STATE-VECTOR',sv)
       lenMix = sv(3)*sv(4)
       if ((sv(1)/=G_Car2d).or.(sv(9)==0)) cycle !pas d'info a retirer
       call LCMLEN(cIp,'CELL        ',lenCell,typ)
       lenCell=lenCell/3
       allocate(mix(lenMix),turn(lenMix),cell(lenCell),merg(lenMix),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: creeDico(1) => allocation pb")
       call LCMLEN(cIp,'TURN        ',long,typ)
       if (long==0) then
          turn(:lenMix) = 1
       else
          call LCMGET(cIp,'TURN        ',turn)
       end if
       call LCMLEN(cIp,'MERGE       ',long,typ)
       if (long==0) then
          merg(:lenMix) = (/(i,i=1,lenMix)/)
       else
          call LCMGET(cIp,'MERGE       ',merg)
       end if
       call LCMGET(cIp,'MIX         ',mix)
       !on teste si tous les milieux sont bien negatifs
       do i = 1,lenMix
          if (mix(i)<0) cycle
          call XABORT("G2S: error, meltig MIX and CELL not supported")
       end do
       call LCMGTC(cIp,'CELL        ',12,lenCell,cell)
       !travail sur l'axe des x
       allocate(tabCh(2*sv(3)-1),turnStr(sv(3)),mergStr(sv(3)),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: creeDico(2) => allocation pb")
       tabCh = '+'
       sideName = 'x' // carreName(1:11)
       dJp=LCMDID(dIp,sideName)
       do j = 1,sv(4)
          do i = 1,sv(3)
             ind = i+(j-1)*sv(3)
             tmpStr = tourne('x',turn(ind)) // cell(-mix(ind))(1:11)
             tabCh(2*i-1) = tmpStr
             turnStr(i) = turn(ind)
             mergStr(i) = merg(ind)
          end do
          lignNbr = lignNbr + 1
          nbStr = i2s(lignNbr)
          arrayName = "array" // nbStr(:7)
          call LCMPTC(dJp,arrayName,12,2*sv(3)-1,tabCh)
          arrayName = "turns" // nbStr(:7)
          call LCMPUT(dJp,arrayName,sv(3),1,turnStr)
          arrayName = "merge" // nbStr(:7)
          call LCMPUT(dJp,arrayName,sv(3),1,mergStr)
          call exploiteStr(tabCh,cIp,sideName,dicoIp,lignNbr)
       end do
       deallocate(tabCh,turnStr,mergStr)
       !travail sur l'axe des y
       allocate(tabCh(2*sv(4)-1),turnStr(sv(4)),mergStr(sv(4)),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: creeDico(3) => allocation pb")
       tabCh = '+'
       sideName(1:1) = 'y'
       dJp=LCMDID(dIp,sideName)
       do i = 1,sv(3)
          do j = 1,sv(4)
             ind = i+(j-1)*sv(3)
             tmpStr = tourne('y',turn(ind)) // cell(-mix(ind))(1:11)
             tabCh(2*j-1) = tmpStr
             turnStr(j) = turn(ind)
             mergStr(j) = merg(ind)
          end do
          lignNbr = lignNbr + 1
          nbStr = i2s(lignNbr)
          arrayName = "array" // nbStr(:7)
          call LCMPTC(dJp,arrayName,12,2*sv(4)-1,tabCh)
          arrayName = "turns" // nbStr(:7)
          call LCMPUT(dJp,arrayName,sv(4),1,turnStr)
          arrayName = "merge" // nbStr(:7)
          call LCMPUT(dJp,arrayName,sv(4),1,mergStr)
          call exploiteStr(tabCh,cIp,sideName,dicoIp,lignNbr)
       end do
       deallocate(tabCh,turnStr,mergStr)
       !fin du travail sur la cellule
       deallocate(mix,turn,cell,merg)
    end do
  end subroutine creeDico
  
  function tourne(axe,turn)
    character,intent(in) :: axe
    integer,intent(in)   :: turn
    character            :: tourne 
    
    tourne = axe
    if (mod(turn,2)==1) return
    if (axe=='x') then
       tourne = 'y'
    else
       tourne = 'x'
    end if
  end function tourne

  subroutine exploiteStr(inStr,carreIp,sideName,dicoIp,lignNbr)
    character*12,dimension(:),intent(in) :: inStr
    type(c_ptr),intent(in)               :: carreIp,dicoIp
    character*12,intent(in)              :: sideName
    integer,intent(inout)                :: lignNbr

    character*12,dimension((size(inStr)+1)/2) :: tmpStr,workStr
    character*12,dimension(:),allocatable :: resStr
    real,dimension(:),allocatable :: mesh
    character*12,dimension(1) :: name2
    type(c_ptr)  :: workIp
    integer      :: i,j,k,nbOcc,szTS,szIS,nbRest,long,typ
    character*12 :: arrayName,name1,meshName,number,text12
    real         :: value,res

    workIp = dicoIp
    szIS = size(inStr)
    szTS = (szIS+1)/2
    !on enleve les "+" dans le tableau
    do i = 1,szTS
       tmpStr(i)=inStr(2*i-1)
    enddo
    !creation a partir de r=a+b+c+b+a, de a=r-b-c-b/2,
    ! b=r-a-c-a/2 et c=r-a-b-b-a
    do i = 1,szTS
       nbOcc = 0
       do j = 1,szTS
          if (tmpStr(j)==tmpStr(i)) then
             nbOcc = nbOcc + 1
          else if (j > nbOcc) then
             workStr(j-nbOcc) = tmpStr(j)
          end if
       end do
       nbRest = szTS - nbOcc
       if (nbOcc==1) then
          allocate(resStr(2*nbRest+1),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: exploiteStr(2) => allocation pb")
       else
          allocate(resStr(2*nbRest+3),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: exploiteStr(3) => allocation pb")
       end if
       resStr(1) = sideName
       do k = 1,nbRest
          resStr(2*k)   = '-'
          resStr(2*k+1) = workStr(k)
       end do
       long=2*nbRest+1
       if (nbOcc/=1) then
          resStr(2*nbRest+2) = '/'
          number = i2s(nbOcc)
          resStr(2*nbRest+3) = number
          long=2*nbRest+3
       end if
       lignNbr = lignNbr + 1
       text12 = i2s(lignNbr)
       arrayName = "expl" // text12(:8)
       call LCMSIX(workIp,tmpStr(i),1)
       call LCMPTC(workIp,arrayName,12,long,resStr)
       call LCMSIX(workIp,tmpStr(i),2)
       deallocate(resStr)
       !exploitation orthogonale de la chaine ('xa + xb' => ya=yb)
       name1 = tourne(tmpStr(i)(1:1),2) // tmpStr(i)(2:12)
       do k = 1,nbRest
          name2(1) = tourne(workStr(k)(1:1),2) // workStr(k)(2:12)
          lignNbr = lignNbr + 1
          arrayName = "orth" // text12(:8)
          call LCMSIX(workIp,name1,1)
          call LCMPTC(workIp,arrayName,12,1,name2)
          call LCMSIX(workIp,name1,2)
       end do
    end do
    !exploitation eventuelle du maillage en entree
    if (sideName(1:1)=='x') then
       meshName = 'MESHX'
    else if (sideName(1:1)=='y') then
       meshName = 'MESHY'
    else
       call XABORT("G2S: internal error in subroutine exploiteStr")
    end if
    call LCMLEN(carreIp,meshName,long,typ)
    if (long/=0) then
       allocate(mesh(long),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: exploiteStr(4) => allocation pb")
       call LCMGET(carreIp,meshName,mesh)
       do i = 1,szTS
          call getValueOf(tmpStr(i),workIp,value)
          res = mesh(i+1)-mesh(i)
          if ((value>0.) .and. (abs(value-res)>geomPrec)) then
             print *,"  ",meshName," -> ",tmpStr(i)
             call XABORT("G2S: error, incoherent geometry(1)")
          endif
          call setValueOf(tmpStr(i),workIp,res)
       end do
       deallocate(mesh)
    end if
  end subroutine exploiteStr

  subroutine getSquareSide(cIp,axis,side)
    type(c_ptr),intent(in) :: cIp
    integer,intent(in) :: axis
    real,intent(out)   :: side
    
    integer                       :: typ,long
    real,dimension(:),allocatable :: tr
    
    side = -1.
    if (axis==1) then !axe des x
       call LCMLEN(cIp,'MESHX       ',long,typ)
       if (long/=0) then
          allocate(tr(long),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: getSquareSide(1) => allocation pb")
          call LCMGET(cIp,'MESHX       ',tr)
          side = tr(long) - tr(1)
          deallocate(tr)
       end if
    else if (axis==2) then !axe des y
       call LCMLEN(cIp,'MESHY       ',long,typ)
       if (long/=0) then
          allocate(tr(long),stat=alloc_ok)
          if (alloc_ok /= 0) call XABORT("G2S: getSquareSide(2) => allocation pb")
          call LCMGET(cIp,'MESHY       ',tr)
          side = tr(long) - tr(1)
          deallocate(tr)
       end if
    else
       call XABORT("G2S: internal error in subroutine getSquareSide")
    end if
  end subroutine getSquareSide

  subroutine getValueOf(sideName,dicoIp,value)
    character*12,intent(in) :: sideName
    type(c_ptr),intent(in)      :: dicoIp
    real,intent(out)        :: value

    type(c_ptr) :: ip,jp

    ip = dicoIp
    jp = LCMGID(ip,sideName)
    call LCMGET(jp,'value',value)
  end subroutine getValueOf

  subroutine setValueOf(sideName,dicoIp,value)
    character*12,intent(in) :: sideName
    type(c_ptr),intent(in)      :: dicoIp
    real,intent(in)         :: value

    type(c_ptr) :: ip,jp

    ip = dicoIp
    jp = LCMDID(ip,sideName)
    call LCMPUT(jp,'value',1,2,value)
  end subroutine setValueOf

  subroutine resoudDico(dicoIp,resIp)
    type(c_ptr),intent(in) :: dicoIp,resIp

    logical      :: newVal,flag1,flag2,empty,lcm
    character*12 :: sideName,saveSideName,namlcm,nammy
    type(c_ptr)  :: ip,jp
    integer      :: ilong
    real         :: val

    !la resolution se fait en travaillant sur toutes les variables une a une
    !on resout les unes apres les autres toutes les equations regissant
    !chaque variable (si possible)
    ! on teste a chaque passe si on a trouve une nouvelle valeur
    ! si oui on continue, sinon (apres une nouvelle passe infructueuse), on
    ! arrete le traitement et on teste si tous les resultats sont bons
    !  oui -> on sort ; non -> on abandonne car manque de donnees

    flag1 = .true.
    do
       flag2 = flag1 !garde la valeur de flag1 a l'essai precedent
       flag1 = .false.
       sideName = ' '
       call LCMINF(dicoIp,namlcm,nammy,empty,ilong,lcm)
       if (.not. empty) then
          call LCMNXT(dicoIp,sideName)
          saveSideName = sideName
          do !boucle sur toutes les donnees du dico
             newVal = evaluateEquation(dicoIp,sideName)
             !test si newVal a ete vrai au moins une fois
             if (newVal) flag1 = .true.
             call LCMNXT(dicoIp,sideName)
             if (sideName == saveSideName) exit
          end do
       end if
       if (.not.(flag1 .or. flag2)) then
          !il ne s'est rien passe au cours des 2 tentatives
          !precedentes => on abandonne si une valeur est non resolue
          saveSideName = sideName
          do !boucle sur toutes les donnees du dico
             call getValueOf(sideName,dicoIp,val)
             if (val<0) then !une valeur est non resolue -> ABORT
                print *,"Curent values of geometry :"
                saveSideName = sideName
                do !boucle sur toutes les donnees du dico
                   call getValueOf(sideName,dicoIp,val)
                   print *,"  ",sidename," -> ",val
                   call LCMNXT(dicoIp,sideName)
                   if (sideName == saveSideName) exit
                end do
                call XABORT("G2S: not enought data in the geometry(1)")
             end if
             call LCMNXT(dicoIp,sideName)
             if (sideName == saveSideName) exit
          end do
          exit !toutes les valeurs sont bonnes
       end if
    end do
    !recopie des resultats
    ip = dicoIp
    sideName = ' '
    call LCMINF(ip,namlcm,nammy,empty,ilong,lcm)
    if (empty) return
    call LCMNXT(ip,sideName)
    saveSideName = sideName    
    do
       jp=LCMGID(ip,sideName)
       call LCMGET(jp,'value',val)
       call LCMPUT(resIp,sideName,1,2,val)
       call LCMNXT(ip,sideName)
       if (sideName == saveSideName) exit
    end do
  end subroutine resoudDico

  function evaluateEquation(dicoIp,sideName)
    type(c_ptr),intent(in)  :: dicoIp
    character*12,intent(in) :: sideName
    logical                 :: evaluateEquation
    !evaluation une a une des equations de dico relatives a sideName
    !retourne .true. si nouvelle affectation et .false. sinon

    real             :: res,val
    double precision :: dres
    type(c_ptr)      :: sideIp
    integer          :: i,long,typ,ilong
    logical          :: goodLine,empty,lcm
    character*12     :: eqName,saveEqName,namlcm,nammy
    character*12,dimension(:),allocatable :: str

    evaluateEquation = .false.
    sideIp = dicoIp
    call LCMSIX(sideIp,sideName,1)
    eqName = ' '
    call LCMINF(sideIp,namlcm,nammy,empty,ilong,lcm)
    if (empty) call XABORT("G2S: intenal error in data structure in &
         &function evaluateEquation")
    call LCMNXT(sideIp,eqName)
    saveEqName = eqName
    do !recherche de la premiere occurence d'une arrayXXX
       call LCMNXT(sideIp,eqName)
       if (eqName(1:5)=='array') exit
       if (eqName==saveEqName) then
          !pas de donnee interessante => on quitte
          return
       end if
    end do
    saveEqName = eqName
    do !boucle sur toutes les equations de sideName
       if (eqName(1:5)/='array') then
          call LCMNXT(sideIp,eqName)
          if (eqName==saveEqName) exit
          cycle
       end if
       !on travaille sur une arrayXXX
       goodLine = .true.
       call LCMLEN(sideIp,eqName,long,typ)
       long = long/3
       allocate(str(long),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("G2S: evaluateEquation => allocation pb")
       call LCMGTC(sideIp,eqName,12,long,str)
       call getValueOf(str(1),dicoIp,res)
       if (res<0) then
          goodLine = .false.
       end if
       dres = res
       do i = 2,long,2
          if (str(i)=='+') then
             call getValueOf(str(i+1),dicoIp,val)
             if (val<0) then 
                goodLine = .false.
                exit
             end if
             dres = dres + val
          else if (str(i)=='-') then
             call getValueOf(str(i+1),dicoIp,val)
             if (val<0) then 
                goodLine = .false.
                exit
             end if
             dres = dres - val             
          else if (str(i)=='/') then
             val = s2i(str(i+1))
             dres = dres / val
          else
             call XABORT("G2S: internal error in function evaluateEquation")
          end if
       end do
       res = real(dres)
       !si on sort normalement de la boucle, on teste la coherence
       !ou on assigne le resultat obtenu
       if (goodLine) then
          call getValueOf(sideName,dicoIp,val)
          if (val<0) then
             call setValueOf(sideName,dicoIp,res)
             evaluateEquation = .true. !on a resolu une nouvelle equation
          else if (abs(val-res)>geomPrec) then
             call XABORT("G2S: error, incoherent geometry(2)")
          end if
       end if
       deallocate(str)
       call LCMNXT(sideIp,eqName)
       if (eqName==saveEqName) exit 
    end do
  end function evaluateEquation

  subroutine prepareCentres(lCentreIpIn,lcIp,sidesIp,dicoIpIn)
    type(c_ptr),intent(in) :: lCentreIpIn,lcIp,sidesIp,dicoIpIn

    type(c_ptr)  :: tmpIp,coordCentIp,ip,dicoIp,lCentreIp
    integer      :: i,j,lg,lgc,ty,dimValToPut,nbArray,gigNum,ilong
    logical      :: empty,lcm
    character*12 :: namlcm,nammy
    character*12 :: carreName,saveCarreName,tmpName,eqName,saveEqName
    character*12 :: sideName,saveSideName,arrayName,saveArrayName
    character*12 :: turnName,mergeName
    real         :: saveVal,val
    real,dimension(:),allocatable         :: valToPut,xx,yy
    integer,dimension(40)                 :: sv
    integer,dimension(:),allocatable      :: turns,merg
    character*12,dimension(:),allocatable :: eqStr,arrayNameStr

    lCentreIp = lCentreIpIn
    dicoIp = dicoIpIn
    !initialisation de la liste et creation et initialisation
    !d'une liste temporaire a partir du dictionnaire utilise a
    !l'etape precedente, et d'un autre temporaire appelle coordCent
    call LCMOP(tmpIp,'lCtemp      ',0,1,0)
    carreName = ' '
    call LCMINF(lcIp,namlcm,nammy,empty,ilong,lcm)
    if (empty) call XABORT("G2S: intenal error in data structure in &
         &subroutine prepareCentres")
    call LCMNXT(lcIp,carreName)
    saveCarreName = carreName
    do
       !travail sur lCtemp
       call LCMSIX(tmpIp,carreName,1)
       tmpName = 'x' // carreName(1:11)
       call LCMSIX(dicoIp,tmpName,1)
       eqName = ' '
       call LCMNXT(dicoIp,eqName)
       saveEqName = eqName
       do
          if (eqName(1:5)=='array') then !l'equation est a conserver
             call putEquationIn(dicoIp,eqName,tmpIp)
          end if
          call LCMNXT(dicoIp,eqName)
          if (eqName==saveEqName) exit
       end do
       call LCMSIX(dicoIp,tmpName,2)
       call LCMSIX(tmpIp,carreName,2)
       !travail sur lCentre
       call LCMSIX(lCentreIp,carreName,1)
       call LCMSIX(lCentreIp,carreName,2)
       call LCMNXT(lcIp,carreName)
       if (carreName==saveCarreName) exit
    end do

    !creation d'un autre temporaire appelle coordCent
    call LCMOP(coordCentIp,'lCoordCent  ',0,1,0)
    sideName = ' '
    call LCMNXT(dicoIp,sideName)
    saveSideName = sideName
    do
       call LCMSIX(dicoIp,sideName,1)
       eqName = ' '
       call LCMNXT(dicoIp,eqName)
       saveEqName = eqName
       do
          if (eqName(1:5)=='array') then !l'equation est a exploiter
             call LCMLEN(dicoIp,eqName,lg,ty)
             lg = lg / 3
             dimValToPut = (lg+1) / 2
             allocate(eqStr(lg),valToPut(dimValToPut),stat=alloc_ok)
             if (alloc_ok /= 0) call XABORT("G2S: prepareCentres(1) => allocation pb")
             call LCMGTC(dicoIp,eqName,12,lg,eqStr)
             !calcul des coordonnees des sous cellules de la cellule
             !consideree et stockage de la valeur obtenu dans coordCent
             call LCMGET(sidesIp,sideName,val)
             valToPut(1) = -0.5 * val
             call LCMGET(sidesIp,eqStr(1),val)
             valToPut(1) = valToPut(1) + 0.5 * val
             do i = 2,dimValToPut
                saveVal =  val
                call LCMGET(sidesIp,eqStr(2*i-1),val)
                valToPut(i) = valToPut(i-1) + 0.5 * (val+saveVal)
             end do
             call LCMPUT(coordCentIp,sideName,dimValToPut,2,valToPut)
             deallocate(eqStr,valToPut)
             exit !on sort car on n'a pas a verifier la coherence
          end if
          call LCMNXT(dicoIp,eqName)
          if (eqName==saveEqName) exit
       end do
       call LCMSIX(dicoIp,sideName,2)
       call LCMNXT(dicoIp,sideName)
       if (sideName==saveSideName) exit
    end do

    carreName = ' '
    call LCMNXT(lcIp,carreName)
    saveCarreName = carreName
    do
       ip=LCMGID(lcIp,carreName)
       call LCMGET(ip,'STATE-VECTOR',sv)
       !recuperation des noms des equations utiles
       call LCMSIX(tmpIp,carreName,1)
       call LCMINF(tmpIp,namlcm,nammy,empty,ilong,lcm)
       if (.not. empty) then
          allocate(arrayNameStr(max(sv(3),sv(4))))
          nbArray = 0
          arrayName = ' '
          call LCMNXT(tmpIp,arrayName)
          saveArrayName = arrayName
          do
             if (arrayName(1:5)/='array') then
                call LCMNXT(tmpIp,arrayName)
                if (arrayName == saveArrayName) exit
                cycle
             end if
             nbArray = nbArray + 1
             arrayNameStr(nbArray) = arrayName
             call LCMNXT(tmpIp,arrayName)
             if (arrayName == saveArrayName) exit
          end do
          !trie des noms des equations
          call sortTabCh12(arrayNameStr,nbArray)
          !travail sur les equations
          gigNum = 0
          do i = 1,nbArray
             arrayName = arrayNameStr(i)
             call LCMLEN(tmpIp,arrayName,lgc,ty)
             lgc = lgc / 3
             allocate(eqStr(lgc))
             lg = (lgc+1) / 2
             allocate(turns(lg),merg(lg),xx(lg),yy(nbArray))
             call LCMGTC(tmpIp,arrayName,12,lgc,eqStr)
             turnName = 'turns' // arrayName(6:12)
             call LCMGET(tmpIp,turnName,turns)
             mergeName = 'merge' // arrayName(6:12)
             call LCMGET(tmpIp,mergeName,merg)
             do j = 1,lg
                gigNum = gigNum + 1
                sideName = 'x' // carreName(1:11)
                call LCMGET(coordCentIp,sideName,xx)
                sideName(1:1) = 'y'
                call LCMGET(coordCentIp,sideName,yy)
                sideName=eqStr(2*j-1)(2:12)
                !ajout du resultat dans lCentre
                call addInRes(lCentreIp,sideName,carreName, &
                     & xx(j),yy(i),turns(j),gigNum,merg(j))
             end do
             deallocate(eqStr,turns,merg,xx,yy)
          end do
          deallocate(arrayNameStr)
       end if

       call LCMSIX(tmpIp,carreName,2)
       call LCMNXT(lcIp,carreName)
       if (carreName==saveCarreName) exit
    end do
    !destruction du repertoire 'root' dans lCentre
    call LCMDEL(lCentreIp,'root        ')
    !fermeture des temporaires
    call LCMCL(coordCentIp,2)
    call LCMCL(tmpIp,2)
  end subroutine prepareCentres

  subroutine addInRes(ipIn,extName,intName,xx,yy,turn,gigNum,merg)
    type(c_ptr),intent(in)  :: ipIn
    integer,intent(in)      :: turn,gigNum,merg
    character*12,intent(in) :: extName,intName
    real,intent(in)         :: xx,yy

    real,dimension(:),allocatable    :: tabR
    integer,dimension(:),allocatable :: tabI
    character*12 :: reg,gigName
    integer      :: lg,ty
    type(c_ptr)  :: ip

    ip = ipIn
    call LCMSIX(ip,extName,1)
    call LCMSIX(ip,intName,1)
    reg = 'x'
    call LCMLEN(ip,reg,lg,ty)
    allocate(tabR(lg+1))
    !ajout de xx
    if (lg /= 0) call LCMGET(ip,reg,tabR)
    tabR(lg+1) = xx
    call LCMPUT(ip,reg,lg+1,2,tabR)
    !ajout de yy
    reg = 'y'
    if (lg /= 0) call LCMGET(ip,reg,tabR)
    tabR(lg+1) = yy
    call LCMPUT(ip,reg,lg+1,2,tabR)
    deallocate(tabR)
    !ajout de turn
    reg = 't'
    allocate(tabI(lg+1))
    if (lg /= 0) call LCMGET(ip,reg,tabI)
    tabI(lg+1) = turn
    call LCMPUT(ip,reg,lg+1,1,tabI)
    deallocate(tabI)
    !ajout de la gigogne et du merge
    gigName = i2s(lg+1)
    call LCMSIX(ip,gigName,1)
    reg = 'gig'
    call LCMPUT(ip,reg,1,1,gigNum)
    reg = 'mrg'
    call LCMPUT(ip,reg,1,1,merg)
    call LCMSIX(ip,gigName,2)
    call LCMSIX(ip,intName,2)
    call LCMSIX(ip,extName,2)
  end subroutine addInRes

  subroutine sortTabCh12(tab,sz)
    character*12,dimension(:),intent(inout) :: tab
    integer,intent(in)                      :: sz

    integer      :: indMax,i,s1,s2
    logical      :: trie
    character*12 :: tmp

    do indMax = sz,2,-1
       trie = .true.
       do i = 1,indMax-1
          s1 = s2i(tab(i)(6:12))
          s2 = s2i(tab(i+1)(6:12))
          if (s2<s1) then
             tmp=tab(i+1) ; tab(i+1)=tab(i) ; tab(i)=tmp
             trie = .false.
          end if
       end do
       if (trie) exit
    end do
  end subroutine sortTabCh12

  subroutine putEquationIn(dicoIp,eqName,lCentreIp)
    type(c_ptr),intent(in)   :: dicoIp,lCentreIp
    character*12,intent(in)  :: eqName

    character*12,dimension(:),allocatable :: str12
    integer,dimension(:),allocatable     :: turns,merg
    character*12 :: turnsName,mergeName
    integer      :: lg,ty,lgt

    call LCMLEN(dicoIp,eqName,lg,ty)
    if (lg/=0) then
       turnsName = 'turns' // eqName(6:12)
       mergeName = 'merge' // eqName(6:12)
       lg = lg/3
       lgt = (lg+1)/2
       allocate(str12(lg),turns(lgt),merg(lgt))
       call LCMGTC(dicoIp,eqName,12,lg,str12)
       call LCMGET(dicoIp,turnsName,turns)
       call LCMGET(dicoIp,mergeName,merg)
       call LCMPTC(lCentreIp,eqName,12,lg,str12)
       call LCMPUT(lCentreIp,turnsName,lgt,1,turns)
       call LCMPUT(lCentreIp,mergeName,lgt,1,merg)
       deallocate(str12,turns,merg)
    else
       call XABORT("G2S: internal error in subroutine putEquationIn")
    end if
  end subroutine putEquationIn

  subroutine prepareMinMax(lMinMaxIp,lCentreIp)
    type(c_ptr),intent(in) :: lMinMaxIp,lCentreIp

    character*12 :: carreName,saveCarreName,namlcm,nammy
    type(c_ptr)  :: ipLC
    integer      :: lg,ty,ilong
    logical      :: empty,lcm
    real         :: minX,maxX,minY,maxY
    real,dimension(:),allocatable :: xx,yy

    ipLC = lCentreIp
    minX = 0. ; maxX = 0. ; minY = 0. ; maxY = 0.
    carreName = ' '
    call LCMINF(ipLC,namlcm,nammy,empty,ilong,lcm)
    if (.not. empty) then !plus d'une seule cellule
       call LCMNXT(ipLC,carreName)
       saveCarreName = carreName
       do
          call LCMSIX(ipLC,carreName,1)
          call LCMLEN(ipLC,'root        ',lg,ty)
          if (lg/=0) then
             call LCMSIX(ipLC,'root        ',1)
             call LCMLEN(ipLC,'x           ',lg,ty)
             allocate(xx(lg))
             call LCMGET(ipLC,'x           ',xx)
             call LCMLEN(ipLC,'y           ',lg,ty)
             allocate(yy(lg))
             call LCMGET(ipLC,'y           ',yy)
             minX = min(minX,minval(xx))
             maxX = max(maxX,maxval(xx))
             minY = min(minY,minval(yy))
             maxY = max(maxY,maxval(yy))
             deallocate(xx,yy)
             call LCMSIX(ipLC,'root        ',2)
          end if
          call LCMSIX(ipLC,carreName,2)
          call LCMNXT(ipLC,carreName)
          if (carreName==saveCarreName) exit
       end do
    end if
    call LCMPUT(lMinMaxIp,'minX        ',1,2,minX)
    call LCMPUT(lMinMaxIp,'maxX        ',1,2,maxX)
    call LCMPUT(lMinMaxIp,'minY        ',1,2,minY)
    call LCMPUT(lMinMaxIp,'maxY        ',1,2,maxY)
  end subroutine prepareMinMax

  recursive subroutine calculeCentre(resIp,lcIp,deep)
    type(c_ptr),intent(inout) :: resIp
    type(c_ptr),intent(in)    :: lcIp
    integer,intent(in)    :: deep

    character*12 :: dirName,saveDirName,subDirName,saveSubDirName
    character*12 :: namlcm,nammy
    type(c_ptr)  :: valuesIp,tmpLcIp,cpResIp,tmpDic
    integer      :: i,lg,ty,nbGig,ilong
    logical      :: empty,lcm

    call LCMEQU(lcIp,resIp)
    tmpLcIp = lcIp
    dirName = ' '
    call LCMINF(tmpLcIp,namlcm,nammy,empty,ilong,lcm)
    if (empty) then
       !une seule cellule => traitement particulier
       call LCMSIX(resIp,'root        ',1)
       call LCMSIX(resIp,'root        ',1)
       call LCMPUT(resIp,'x           ',1,2,0.)
       call LCMPUT(resIp,'y           ',1,2,0.)
       call LCMPUT(resIp,'t           ',1,1,1)
       call LCMSIX(resIp,'1           ',1)
       call LCMPUT(resIp,'gig         ',1,1,1)
       call LCMPUT(resIp,'mrg         ',1,1,1)
       call LCMSIX(resIp,'1           ',2)
       call LCMSIX(resIp,'root        ',2)
       call LCMSIX(resIp,'root        ',2)
       return
    end if
    call LCMOP(valuesIp,'temporariObj',0,1,0)
    call LCMNXT(tmpLcIp,dirName)
    saveDirName = dirName
    nbGig = 0
    do
       call LCMSIX(tmpLcIp,dirName,1)
       subDirName = ' '
       call LCMINF(tmpLcIp,namlcm,nammy,empty,ilong,lcm)
       if (.not. empty) then
          call LCMNXT(tmpLcIp,subDirName)
          saveSubDirName = subDirName
          do
             call LCMLEN(valuesIp,subDirName,lg,ty)
             if (lg==0) then
                nbGig = nbGig + 1
                call LCMSIX(valuesIp,subDirName,1)
                call LCMSIX(valuesIp,subDirName,2)
             end if
             call LCMNXT(tmpLcIp,subDirName)
             if (subDirName==saveSubDirName) exit
          end do
       end if
       call LCMSIX(tmpLcIp,dirName,2)
       call LCMNXT(tmpLcIp,dirName)
       if (dirName==saveDirName) exit
    end do
    !test de sortie de la recurrence
    if (nbGig <= 1) then
       call LCMEQU(lcIp,resIp)
       call LCMCL(valuesIp,2)
       return
    end if
    !liste des cellules a garder
    dirName = ' '
    call LCMNXT(lcIp,dirName)
    saveDirName = dirName
    do
       call LCMLEN(valuesIp,dirName,lg,ty)
       if (lg==0) call LCMDEL(resIp,dirName)
       call LCMNXT(lcIp,dirName)
       if (dirName==saveDirName) exit
    end do

    dirName = ' '
    call LCMNXT(lcIp,dirName)
    saveDirName = dirName
    do
       call LCMLEN(valuesIp,dirName,lg,ty)
       if (lg==0) then
          call LCMSIX(tmpLcIp,dirName,1)
          subDirName = ' '
          call LCMNXT(tmpLcIp,subDirName)
          saveSubDirName = subDirName
          do
             call LCMSIX(tmpLcIp,subDirName,1)
             if (subDirName/='root') then
                call LCMLEN(tmpLcIp,'t           ',lg,ty)
                do i = 1,lg
                   call LCMOP(tmpDic,'tmpDic      ',0,1,0)
                   call translateDic(tmpDic,lcIp,dirName,subDirName,i)
                   call mergeDic(resIp,dirName,tmpDic)
                   call LCMCL(tmpDic,2)
                end do
             else
                call LCMOP(tmpDic,'tmpDic      ',0,1,0)
                call LCMSIX(tmpDic,subDirName,1)
                call LCMEQU(tmpLcIp,tmpDic) !copie uniquement de la branche
                                            !interresante (lc[dirname]['root'])
                call LCMSIX(tmpDic,subDirName,2)
                call mergeDic(resIp,dirName,tmpDic)
                call LCMCL(tmpDic,2)
             end if
             call LCMSIX(tmpLcIp,subDirName,2)
             call LCMNXT(tmpLcIp,subDirName)
             if (subDirName==saveSubDirName) exit
          end do
          call LCMSIX(tmpLcIp,dirName,2)
       end if
       call LCMNXT(lcIp,dirName)
       if (dirName==saveDirName) exit
    end do
    !iteration
    call LCMOP(cpResIp,'resCentres  ',0,1,0)
    call LCMEQU(resIp,cpResIp)
    call calculeCentre(resIp,cpResIp,deep+1)
    call LCMCL(cpResIp,2)
    !elimination des cellules non terminales si deep==0
    if (deep==0) then
       dirName = ' '
       call LCMNXT(valuesIp,dirName)
       saveDirName = dirName
       do
          call LCMLEN(resIp,dirName,lg,ty)
          if (lg/=0) call LCMDEL(resIp,dirName)
          call LCMNXT(valuesIp,dirName)
          if (dirName==saveDirName) exit
       end do
    end if
    call LCMCL(valuesIp,2)
  end subroutine calculeCentre

  subroutine translateDic(resIp,dicIp,dName,sDName,ind)
    type(c_ptr),intent(in)  :: resIp,dicIp
    integer,intent(in)      :: ind
    character*12,intent(in) :: dName,sDName
    !effectue la translation dans le repere de la cellule englobante
    !de la liste des cellules de dicIp[sDName], suivant les vecteurs de
    !dicIp[dName][sDName] (qui comprennent x(i),y(i),t(i),[`i`][gig] i=1,sz)

    character*12 :: dirName,saveDirName,gigName
    type(c_ptr)  :: tmpRes,tmpDic
    integer      :: i,lg,ty,tval,lggigI,lggigval
    real         :: xval,yval
    integer,dimension(:),allocatable :: tt,gigI,gigval,mrgI,mrgval
    real,dimension(:),allocatable    :: xx,yy

    tmpRes = resIp ; tmpDic = dicIp
    !recuperation des donnees de translation
    call LCMSIX(tmpDic,dName,1)
    call LCMSIX(tmpDic,sDName,1)
    call LCMLEN(tmpDic,'t           ',lg,ty)
    allocate(xx(lg),yy(lg),tt(lg))
    call LCMGET(tmpDic,'t           ',tt) ; tval = tt(ind)
    call LCMGET(tmpDic,'x           ',xx) ; xval = xx(ind)
    call LCMGET(tmpDic,'y           ',yy) ; yval = yy(ind)
    deallocate(xx,yy,tt)
    gigName = i2s(ind)
    call LCMSIX(tmpDic,gigName,1)
    call LCMLEN(tmpDic,'gig         ',lggigval,ty)
    allocate(gigval(lggigval),mrgval(lggigval))
    call LCMGET(tmpDic,'gig         ',gigval)
    call LCMGET(tmpDic,'mrg         ',mrgval)
    tmpDic = dicIp
    !translation
    call LCMSIX(tmpDic,sDName,1)
    dirName = ' '
    call LCMNXT(tmpDic,dirName)
    saveDirName = dirName
    do
       call LCMSIX(tmpDic,dirName,1)
       call LCMSIX(tmpRes,dirName,1)
       call LCMLEN(tmpDic,'t           ',lg,ty)
       allocate(xx(lg),yy(lg),tt(lg))
       call LCMGET(tmpDic,'t           ',tt)
       call LCMGET(tmpDic,'x           ',xx)
       call LCMGET(tmpDic,'y           ',yy)
       do i = 1,lg
          gigName = i2s(i)
          call LCMSIX(tmpDic,gigName,1)
          call LCMSIX(tmpRes,gigName,1)
          call LCMLEN(tmpDic,'gig         ',lggigI,ty)
          allocate(gigI(lggigI+lggigval),mrgI(lggigI+lggigval))
          call LCMGET(tmpDic,'gig         ',gigI)
          call LCMGET(tmpDic,'mrg         ',mrgI)
          call translateWithTurn(xx(i),yy(i),tt(i),xval,yval,tval)
          gigI(lggigI+1:lggigI+lggigval) = gigval
          call LCMPUT(tmpRes,'gig         ',lggigI+lggigval,1,gigI)
          mrgI(lggigI+1:lggigI+lggigval) = mrgval
          call LCMPUT(tmpRes,'mrg         ',lggigI+lggigval,1,mrgI)
          deallocate(gigI,mrgI)
          call LCMSIX(tmpRes,gigName,2)
          call LCMSIX(tmpDic,gigName,2)
       end do
       call LCMPUT(tmpRes,'t           ',lg,1,tt)
       call LCMPUT(tmpRes,'x           ',lg,2,xx)
       call LCMPUT(tmpRes,'y           ',lg,2,yy)       
       deallocate(xx,yy,tt)
       call LCMSIX(tmpRes,dirName,2)
       call LCMSIX(tmpDic,dirName,2)
       call LCMNXT(tmpDic,dirName)
       if (dirName==saveDirName) exit
    end do
    deallocate(gigval,mrgval)
  end subroutine translateDic

  subroutine translateWithTurn(xx,yy,tt,xval,yval,tval)
    real,intent(inout)    :: xx,yy
    integer,intent(inout) :: tt
    real,intent(in)       :: xval,yval
    integer,intent(in)    :: tval

    select case(geomTyp)
    case(RecTyp)
       select case (tt)
       case(1)
          xx=xx+xval ; yy=yy+yval
       case(2)
          xx=xx+yval ; yy=yy-xval
       case(3)
          xx=xx-xval ; yy=yy-yval
       case(4)
          xx=xx-yval ; yy=yy+xval
       case(5)
          xx=xx-xval ; yy=yy+yval
       case(6)
          xx=xx+yval ; yy=yy+xval
       case(7)
          xx=xx+xval ; yy=yy-yval
       case(8)
          xx=xx-yval ; yy=yy-xval
       end select
       tt = composeRotRec(tval,tt)
     case(TriaTyp)
       select case (tt)
       case(1)
          xx=xx+xval ; yy=yy+yval
       case(2)
          xx=xx+0.5*(xval-sqrt_3f*yval)
          yy=yy+0.5*(yval+sqrt_3f*xval)
       case(3)
          xx=xx+0.5*(-xval-sqrt_3f*yval)
          yy=yy+0.5*(-yval+sqrt_3f*xval)
       case(4)
          xx=xx-xval ; yy=yy-yval
       case(5)
          xx=xx+0.5*(-xval+sqrt_3f*yval)
          yy=yy+0.5*(-yval-sqrt_3f*xval)
       case(6)
          xx=xx+0.5*(xval+sqrt_3f*yval)
          yy=yy+0.5*(yval-sqrt_3f*xval)
       case(7)
          xx=xx+xval ; yy=yy-yval
       case(8)
          xx=xx+0.5*(xval-sqrt_3f*yval)
          yy=yy-0.5*(yval+sqrt_3f*xval)
       case(9)
          xx=xx+0.5*(-xval-sqrt_3f*yval)
          yy=yy-0.5*(-yval+sqrt_3f*xval)
       case(10)
          xx=xx-xval ; yy=yy+yval
       case(11)
          xx=xx+0.5*(-xval+sqrt_3f*yval)
          yy=yy+0.5*(yval+sqrt_3f*xval)
       case(12)
          xx=xx+0.5*(xval+sqrt_3f*yval)
          yy=yy-0.5*(yval-sqrt_3f*xval)
       end select
       tt = composeRotTri(tval,tt)
    case default
       call XABORT("G2S: internal error in subroutine translateWithTurn")
    end select
  end subroutine translateWithTurn

  subroutine mergeDic(resIp,dirName,otherDicIp)
    type(c_ptr),intent(in)  :: resIp,otherDicIp
    character*12,intent(in) :: dirName
    !concatene deux dictionnaires (res[dirName] et otherDic)
    !contenants des positions (x,y,t,(1:gig),..,(size(t):gig))

    character*12 :: sdName,nGig,namlcm,nammy
    type(c_ptr)  :: tmpResIp,tmpOtherIp
    integer      :: i,lg,ty,lgN,lgO,ilong
    logical      :: empty,lcm
    real,dimension(:),allocatable    :: nra,ora
    integer,dimension(:),allocatable :: nia,oia,mrg

    tmpResIp = resIp ; tmpOtherIp = otherDicIp
    call LCMSIX(tmpResIp,dirName,1)
    call LCMINF(tmpOtherIp,namlcm,nammy,empty,ilong,lcm)
    if (empty) return !###
    sdName = ' '
    call LCMNXT(tmpOtherIp,sdName)
    call LCMSIX(tmpOtherIp,sdName,1)
    call LCMSIX(tmpResIp,sdName,1)
    call LCMLEN(tmpResIp,'t           ',lgN,ty)
    if (lgN/=0) then !le repertoire existait deja => ajout de donnees
       call LCMLEN(tmpOtherIp,'t           ',lgO,ty)
       allocate(nra(lgO+lgN),ora(lgO),nia(lgO+lgN),oia(lgO))
       !t
       call LCMGET(tmpResIp,'t           ',nia)
       call LCMGET(tmpOtherIp,'t           ',oia)
       nia(lgN+1:lgO+lgN) = oia
       call LCMPUT(tmpResIp,'t           ',lgO+lgN,1,nia)
       !x
       call LCMGET(tmpResIp,'x           ',nra)
       call LCMGET(tmpOtherIp,'x           ',ora)
       nra(lgN+1:lgO+lgN) = ora
       call LCMPUT(tmpResIp,'x           ',lgO+lgN,2,nra)
       !y
       call LCMGET(tmpResIp,'y           ',nra)
       call LCMGET(tmpOtherIp,'y           ',ora)
       nra(lgN+1:lgO+lgN) = ora
       call LCMPUT(tmpResIp,'y           ',lgO+lgN,2,nra)
       deallocate(nra,ora,nia,oia)
       lg = lgO
       do i = 1,lg
          nGig = i2s(i)
          call LCMSIX(tmpOtherIp,nGig,1)
          call LCMLEN(tmpOtherIp,'gig         ',lgO,ty)
          allocate(oia(lgO),mrg(lgO))
          call LCMGET(tmpOtherIp,'gig         ',oia)
          call LCMGET(tmpOtherIp,'mrg         ',mrg)
          call LCMSIX(tmpOtherIp,nGig,2)
          nGig = i2s(i+lgN)
          call LCMSIX(tmpResIp,nGig,1)
          call LCMPUT(tmpResIp,'gig         ',lgO,1,oia)
          call LCMPUT(tmpResIp,'mrg         ',lgO,1,mrg)
          call LCMSIX(tmpResIp,nGig,2)
          deallocate(oia,mrg)
       end do
    else !le repertoire n'existait pas encore => affectation directe
       !call LCMEQU(tmpOtherIp,tmpResIp)
       call LCMLEN(tmpOtherIp,'t           ',lg,ty)
       allocate(ora(lg),oia(lg))
       !t
       call LCMGET(tmpOtherIp,'t           ',oia)
       call LCMPUT(tmpResIp,'t           ',lg,1,oia)
       !x
       call LCMGET(tmpOtherIp,'x           ',ora)
       call LCMPUT(tmpResIp,'x           ',lg,2,ora)
       !y
       call LCMGET(tmpOtherIp,'y           ',ora)
       call LCMPUT(tmpResIp,'y           ',lg,2,ora)
       deallocate(ora,oia)
       do i = 1,lg
          nGig = i2s(i)
          call LCMSIX(tmpOtherIp,nGig,1)
          call LCMLEN(tmpOtherIp,'gig         ',lgO,ty)
          allocate(oia(lgO),mrg(lgO))
          call LCMGET(tmpOtherIp,'gig         ',oia)
          call LCMGET(tmpOtherIp,'mrg         ',mrg)
          call LCMSIX(tmpOtherIp,nGig,2)
          call LCMSIX(tmpResIp,nGig,1)
          call LCMPUT(tmpResIp,'gig         ',lgO,1,oia)
          call LCMPUT(tmpResIp,'mrg         ',lgO,1,mrg)
          call LCMSIX(tmpResIp,nGig,2)
          deallocate(oia,mrg)
       end do
    end if
  end subroutine mergeDic

  subroutine compileResutats(geoIp,lCentreIp,lSideIp,boundDataIp)
    type(c_ptr),intent(in) :: geoIp,lCentreIp,lSideIp,boundDataIp

    type(c_ptr)  :: geo,centres
    integer      :: i,lg,ty,lgG
    character*12 :: carreName,saveCarreName,gigName,num,mrgName
    real,dimension(2)                :: sideXY
    real,dimension(4)                :: minmaxXY
    real,dimension(:),allocatable    :: xx,yy
    integer,dimension(:),allocatable :: tt,gigI,mrgI

    geo = geoIp ; centres = lCentreIp
    call LCMSIX(geo,'NEW-DATA    ',1)
    carreName = ' '
    call LCMNXT(centres,carreName)
    saveCarreName = carreName
    do
       call LCMSIX(geo,carreName,1)
       !SIDEXY
       call LCMGET(lSideIp,'x'//carreName,sideXY(1))
       call LCMGET(lSideIp,'y'//carreName,sideXY(2))
       if (geomTyp==RecTyp) then
          call LCMPUT(geo,'SIDEXY      ',2,2,sideXY)
       else
          call LCMPUT(geo,'SIDEXY      ',2,2,sideXY(1))
       end if
       !COORDX , COORDY , TURN , POSi
       call LCMSIX(centres,carreName,1)
       call LCMSIX(centres,'root        ',1)
       call LCMLEN(centres,'t           ',lg,ty)
       allocate(xx(lg),yy(lg),tt(lg))
       call LCMGET(centres,'x           ',xx)
       call LCMPUT(geo,'COORDX      ',lg,2,xx)
       call LCMGET(centres,'y           ',yy)
       call LCMPUT(geo,'COORDY      ',lg,2,yy)
       call LCMGET(centres,'t           ',tt)
       call LCMPUT(geo,'TURN        ',lg,1,tt)
       do i = 1,lg
          num = i2s(i)
          gigName = num
          call LCMSIX(centres,gigName,1)
          call LCMLEN(centres,'gig         ',lgG,ty)
          allocate(gigI(lgG),mrgI(lgG))
          call LCMGET(centres,'gig         ',gigI)
          call LCMGET(centres,'mrg         ',mrgI)
          call LCMSIX(centres,gigName,2)
          gigName = 'POS' // num(:9)
          call LCMPUT(geo,gigName,lgG,1,gigI)
          mrgName = 'MRG' // num(:9)
          call LCMPUT(geo,mrgName,lgG,1,mrgI)
          deallocate(gigI,mrgI)
       end do
       deallocate(xx,yy,tt)
       call LCMSIX(centres,'root        ',2)
       call LCMSIX(centres,carreName,2)
       !on cycle
       call LCMSIX(geo,carreName,2)
       call LCMNXT(centres,carreName)
       if (carreName==saveCarreName) exit
    end do
    !donnes sur les CL
    call LCMSIX(geo,'BOUND-DATA  ',1)
    !SIDEXY
    call LCMGET(lSideIp,'xroot       ',sideXY(1))
    call LCMGET(lSideIp,'yroot       ',sideXY(2))
    call LCMPUT(geo,'SIDEXY      ',2,2,sideXY)
    !MINMAXXY
    call LCMGET(boundDataIp,'minX        ',minmaxXY(1))
    call LCMGET(boundDataIp,'minY        ',minmaxXY(2))
    call LCMGET(boundDataIp,'maxX        ',minmaxXY(3))
    call LCMGET(boundDataIp,'maxY        ',minmaxXY(4))
    call LCMPUT(geo,'MINMAXXY    ',4,2,minmaxXY)
  end subroutine compileResutats

  subroutine creerEtChargerNewDataHex(geoIp,szB,szP)
    type(c_ptr),intent(in)    :: geoIp
    integer,intent(inout) :: szB,szP
    
    real,dimension(:),allocatable :: txx,tyy
    integer,dimension(:),allocatable :: turn,mix,merg
    character*12,dimension(:),allocatable :: cells
    real                  :: side
    double precision      :: sdX,sdY
    integer,dimension(40) :: sv
    type(c_ptr)           :: ip
    integer               :: nbH,lg,typ,iHex,i,j,k,l
    integer               :: a,b,aa,bb,nbC,ind,lgSide
    integer,dimension(6)  :: da,db
    real,dimension(4)     :: fooData
    
    ip = geoIp
    !ajout de donnees pour le cas d'une seule cellule
    call LCMSIX(ip,'NEW-DATA    ',1)
    call LCMSIX(ip,'BOUND-DATA  ',1)
    fooData = 0.
    call LCMPUT(ip,'SIDEXY      ',4,2,fooData)
    call LCMPUT(ip,'MINMAXXY    ',4,2,fooData)
    call LCMSIX(ip,'BOUND-DATA  ',2)
    call LCMSIX(ip,'NEW-DATA    ',2)
    !placement des cellules
    call LCMGET(ip,'STATE-VECTOR',sv)
    nbH=sv(3)
    nbC=sv(9)
    if ((sv(8)==0).and.(nbH==1)) then !pas de sous-cellules
       call chargerNewData(geoIp,szB,szP)
       return
    end if
    allocate(txx(nbH),tyy(nbH),turn(nbH),mix(nbH),merg(nbH))
    call LCMGET(ip,'MIX         ',mix)
    call LCMLEN(ip,'TURN        ',lg,typ)
    if (lg/=0) then
       call LCMGET(ip,'TURN        ',turn)
    else
       turn = 1
    end if
    call LCMLEN(ip,'MERGE       ',lg,typ)
    if (lg/=0) then
       call LCMGET(ip,'MERGE       ',merg)
    else
       merg(:nbH) = (/(i,i=1,nbH)/)
    end if
    call LCMGET(ip,'IHEX        ',iHex)
    lgSide = 0 
    select case(iHex)
    case(H_S30)
       a = -1 ; l = 0 ; i = 0
       do
          i = i + 1
          do j = 0,1
             a = a + 1 ; b = j
             do k = 1,i
                l = l + 1 ; txx(l) = a ; tyy(l) = b ; lgSide = 2*b ; b = b + 2
                if (l==nbH) goto 10
             end do
          end do
       end do
    case(H_SA60)
       l = 0 ; i = -1
       do
          i = i + 1 ; a = i ; b = -i
          do j = 0,i
             l = l + 1 ; txx(l) = a ; tyy(l) = b ; lgSide = 2*b ; b = b + 2
             if (l==nbH) goto 10
          end do
       end do
    case(H_SB60)
       a = -1 ; l = 0 ; i = 0
       do
          i = i + 1
          do j = 0,1
             a = a + 1 ; b = j
             do k = 0,i
                l = l + 1 ; txx(l) = a ; tyy(l) = b ; lgSide = 2*b ; b = b + 2
                if (l==nbH) goto 10
             end do
             aa = a ; bb = b
             do k = 0,i-1
                l = l + 1 ; aa = aa - 1 ; bb = bb + 1
                txx(l) = aa ; tyy(l) = bb
                if (l==nbH) goto 10
             end do
          end do
       end do
    case(H_S90)
       da = (/0,-1,0,0,0,0/) ; db = (/2,1,0,0,0,0/) ; l = 1 ; i = 0
       txx(l) = 0 ; tyy(l) = 0
       if (l==nbH) goto 10
       do
          i = i + 1 ; a = i ; b = -i
          do j = 1,2
             do k = 1,i
                a = a + da(j) ; b = b + db(j)
                if (b>=0) then
                   l = l + 1 ; txx(l) = a ; tyy(l) = b ; lgSide = b
                   if (l==nbH) goto 10
                end if
             end do
          end do
       end do
    case(H_R120)
       da = (/0,-1,0,0,0,0/) ; db = (/2,1,0,0,0,0/) ; l = 1 ; i = 0
       txx(l) = 0 ; tyy(l) = 0
       if (l==nbH) goto 10
       do
          i = i + 1 ; a = i ; b = -i
          do j = 1,2
             do k = 1,i
                a = a + da(j) ; b = b + db(j)
                l = l + 1 ; txx(l) = a ; tyy(l) = b ; lgSide = b
                if (l==nbH) goto 10
             end do
          end do
       end do
    case(H_R180)
       da = (/1,0,-1,0,0,0/) ; db = (/1,2,1,0,0,0/) ; l = 1 ; i = 0
       txx(l) = 0 ; tyy(l) = 0
       if (l==nbH) goto 10
       do
          i = i + 1 ; a = 0 ; b = -2*i
          do j = 1,3
             do k = 1,i
                a = a + da(j) ; b = b + db(j)
                l = l + 1 ; txx(l) = a ; tyy(l) = b ; lgSide = b
                if (l==nbH) goto 10
             end do
          end do
       end do
    case(H_SA180)
       da = (/1,0,-1,0,0,0/) ; db = (/1,2,1,0,0,0/) ; l = 1 ; i = 0
       txx(l) = 0 ; tyy(l) = 0
       if (l==nbH) goto 10
       do
          i = i + 1 ; a = 0 ; b = -2*i
          l = l + 1 ; txx(l) = a ; tyy(l) = b
          if (l==nbH) goto 10
          do j = 1,3
             do k = 1,i
                a = a + da(j) ; b = b + db(j)
                l = l + 1 ; txx(l) = a ; tyy(l) = b ; lgSide = b
                if (l==nbH) goto 10
             end do
          end do
       end do
    case(H_SB180)
       da = (/0,-1,-1,0,0,0/) ; db = (/2,1,-1,-2,0,0/) ; l = 1 ; i = 0
       txx(l) = 0 ; tyy(l) = 0
       if (l==nbH) goto 10
       do
          i = i + 1 ; a = i ; b = -i
          do j = 1,4
             do k = 1,i
                if (b>=0) then
                   l = l + 1 ; txx(l) = a ; tyy(l) = b ; lgSide = -2*a
                   if (l==nbH) goto 10
                end if
                a = a + da(j) ; b = b + db(j)
             end do
          end do
       end do
    case(H_Complete)
       da = (/-1,-1,0,1,1,0/) ; db = (/1,-1,-2,-1,1,2/) ; l = 1 ; i = 0
       txx(l) = 0 ; tyy(l) = 0
       if (l==nbH) goto 10
       do
          i = i + 1 ; a = i ; b = i
          do j = 1,6
             do k = 1,i
                l = l + 1 ; txx(l) = a ; tyy(l) = b ; lgSide = 2*a
                if (l==nbH) goto 10
                a = a + da(j) ; b = b + db(j)
             end do
          end do
       end do
    end select
10  continue !sortie de la boucle de remplissage des positions (bhaaa, un goto)
    call giveSide(ip,szB,side)
    sdX = side*1.5d0 ; sdY = side*0.5d0*sqrt(3.d0)
    bCData%sidexy(1) = side ; bCData%sidexy(2) = lgSide * sdY
    if (sv(8)==0) then !no sub-geometries
       do j = 1,nbH
          szP = szP + 1
          tabCellulePlaced(szP)%indice  = 1
          tabCellulePlaced(szP)%xcenter = txx(j)*sdX
          tabCellulePlaced(szP)%ycenter = tyy(j)*sdY
          tabCellulePlaced(szP)%turn    = turn(j)
          allocate(tabCellulePlaced(szP)%gig(1))
          tabCellulePlaced(szP)%gig(1) = j
          allocate(tabCellulePlaced(szP)%mrg(1))
          tabCellulePlaced(szP)%mrg(1) = merg(j)
       enddo
    else
       allocate(cells(nbC))
       call LCMGTC(ip,'CELL        ',12,nbC,cells)
       do i = 1,szB
          ind = 0
          do j = 1,szB
             if (cells(i)==tabCelluleBase(j)%name) then
                ind = j
                exit
             end if
          end do
          if (ind==0) call XABORT("G2S: internal error in function&
               & creerEtChargerNewDataHex")
          do j = 1,nbH
             if (mix(j)==-i) then
                !on va creer une cellule placee d'indice ind en position j
                szP = szP + 1
                tabCellulePlaced(szP)%indice  = ind
                tabCellulePlaced(szP)%xcenter = txx(j)*sdX
                tabCellulePlaced(szP)%ycenter = tyy(j)*sdY
                tabCellulePlaced(szP)%turn    = turn(j)
                allocate(tabCellulePlaced(szP)%gig(1))
                tabCellulePlaced(szP)%gig(1) = j
                allocate(tabCellulePlaced(szP)%mrg(1))
                tabCellulePlaced(szP)%mrg(1) = merg(j)
             end if
          end do
       end do
       deallocate(cells)
    endif
    deallocate(txx,tyy,turn,mix,merg)
  end subroutine creerEtChargerNewDataHex

  subroutine giveSide(ip,szB,side)
    type(c_ptr),intent(in) :: ip
    integer,intent(in) :: szB
    real,intent(inout) :: side
    
    integer          :: i,lg,typ
    logical          :: gotIt
    double precision :: dside
    
    gotIt = .false.
    call LCMLEN(ip,'SIDE        ',lg,typ)
    if (lg/=0) then
       call LCMGET(ip,'SIDE        ',side)
       gotIt = .true.
       dside = side
    end if
    do i = 1,szB
       if (tabCelluleBase(i)%ok(n_side)) then
          if (gotIt) then
             if (.not. isEqual(tabCelluleBase(i)%side,dside)) &
                  & call XABORT("G2S: Error in the value of argument SIDE&
                  & in the geometry")
          else
             gotIt = .true.
             dside = tabCelluleBase(i)%side
          end if
       else if (gotIt) then
          tabCelluleBase(i)%side = dside
       end if
    end do
    if (.not. gotIt) then
       call XABORT("G2S: Error in the value of argument SIDE in the geometry")
    end if
    do i = 1,szB
       tabCelluleBase(i)%ok(n_side) = .true.
       tabCelluleBase(i)%side = dside
    end do
    side = real(dside)
  end subroutine giveSide

  subroutine creerEtChargerNewDataTri(geoIp,szB,szP)
    type(c_ptr),intent(in):: geoIp
    integer,intent(inout) :: szB,szP
    
    type(c_ptr) :: ip,lcIp,dicoIp,sidesIp,lCentreIp,centreCalculesIp
    integer :: szLc
    
    ip = geoIp
    !creation de la liste des triangles
    call LCMOP(lcIp,'lc          ',0,1,0)
    call creeListeCellules(lcIp,ip,szLc)

    !creation d'une structure pour le calcul des longueurs des cotes
    call LCMOP(dicoIp,'dico        ',0,1,0)
    call creeStructure(dicoIp,lcIp,szLc)

    !creation de la liste des longueurs des cotes et resolution du systeme
    !et destruction de la structure
    call LCMOP(sidesIp,'sides       ',0,1,0)
    call resoudStructure(dicoIp,sidesIp)
    call LCMCL(dicoIp,2)

    !cree et prepare la liste des centres des cellules de base,
    !avec leur gigogne d'origine
    call LCMOP(lCentreIp,'centres     ',0,1,0)
    call prepareCentresTri(lCentreIp,lcIp,szLc,sidesIp)

    !resolution du systeme donnant les coordonnees des centres des cellules
    !de base et stockage dans une nouvelle structure 
    call LCMOP(centreCalculesIp,'resCentres  ',0,1,0)
    call calculeCentre(centreCalculesIp,lCentreIp,0)

    !compilation des resultats et ajout dans la geometrie des nouvelles donnees
    call compileResutatsTri(ip,lcIp,szLc,centreCalculesIp,sidesIp)
    call LCMCL(lcIp,2)
    call LCMCL(centreCalculesIp,2)
    call LCMCL(sidesIp,2)
    !utilisation des donnees (on accroche le traitement classique developpe
    !pour python)
    call chargerNewData(geoIp,szB,szP)
  end subroutine creerEtChargerNewDataTri

  subroutine creeStructure(dicoIp,lcIp,szLc)
    type(c_ptr),intent(in) :: dicoIp,lcIp
    integer,intent(in) :: szLc
    
    character*12 :: triName
    type(c_ptr)  :: triIp,dIp
    integer      :: i,j,k,lg,ty,lenMix,lignNbr
    real         :: side,value
    integer,dimension(40)                 :: sv
    integer,dimension(:),allocatable      :: mix
    character*12,dimension(:),allocatable :: cell
    
    !initialisation
    dIp = dicoIp
    triName = ' '
    do i = 1,szLc
       call LCMNXT(lcIp,triName)
       triIp=LCMGID(lcIp,triName)
       side = -1.
       call LCMLEN(triIp,'SIDE        ',lg,ty)
       if (lg/=0) call LCMGET(triIp,'SIDE        ',side)
       call LCMSIX(dIp,triName,1)
       call LCMPUT(dIp,'value',1,2,side)
       call LCMSIX(dIp,triName,2)
    end do
    !remplissage
    dIp = dicoIp
    triName = ' '
    lignNbr = 0
    do k = 1,szLc
       !tavail sur la kieme cellule
       call LCMNXT(lcIp,triName)
       triIp=LCMGID(lcIp,triName)
       call LCMGET(triIp,'STATE-VECTOR',sv)
       call LCMLEN(triIp,'MIX         ',lenMix,ty)
       if ((lenMix==1).or.(sv(9)==0)) cycle !pas d'info a retirer
       allocate(mix(lenMix),cell(sv(9)))
       call LCMGET(triIp,'MIX         ',mix)
       !on teste si tous les milieux sont bien negatifs
       do i = 1,lenMix
          if (mix(i)<0) cycle
          call XABORT("G2S: error, meltig MIX and CELL not supported")
       end do
       call LCMGTC(triIp,'CELL        ',12,sv(9),cell)
       !on traite les donnees
       !dans le triangle considere
       call LCMSIX(dIp,triName,1)
       value = 1.*sv(3)
       do i = 1,lenMix
          call LCMPUT(dIp,cell(-mix(i)),1,2,value)
       end do
       call LCMSIX(dIp,triName,2)
       !dans les autres
       value = 1./sv(3)
       do i = 1,sv(9)
          call LCMSIX(dIp,cell(i),1)
          do j = 1,sv(9)
             if (j==i) cycle
             call LCMPUT(dIp,cell(j),1,2,1.)
          end do
          call LCMPUT(dIp,triName,1,2,value)
          call LCMSIX(dIp,cell(i),2)
       end do
       deallocate(mix,cell)
    end do
  end subroutine creeStructure

  subroutine resoudStructure(dicoIp,resIp)
    type(c_ptr),intent(in) :: dicoIp,resIp

    type(c_ptr)  :: dIp
    logical      :: fini,flag1,flag2
    character*12 :: triName,saveTriName
    real         :: val

    fini = .false.
    flag1 = .true.
    do
       if (fini) exit
       fini = .true.
       flag2 = flag1 !garde la valeur de flag1 a l'essai precedent
       flag1 = .false.
       triName = ' '
       call LCMNXT(dicoIp,triName)
       saveTriName = triName
       do !boucle sur toutes les donnees du dico
          fini = fini .and. evaluateLine(dicoIp,triName)
          !test si fini a ete vrai au moins une fois
          if (fini) flag1 = .true.
          call LCMNXT(dicoIp,triName)
          if (triName == saveTriName) exit
       end do
       if (.not.(flag1 .or. flag2)) then
          !il ne s'est rien passe au cours des 2 tentatives
          !precedentes => on abandonne
          call XABORT("G2S: not enought data in the geometry(2)")
       end if
    end do
    !recopie des resultats
    dIp = dicoIp
    triName = ' '
    call LCMNXT(dIp,triName)
    saveTriName = triName  
    do
       call LCMSIX(dIp,triName,1)
       call LCMGET(dIp,'value',val)
       call LCMPUT(resIp,triName,1,2,val)
       call LCMSIX(dIp,triName,2)
       call LCMNXT(dIp,triName)
       if (triName == saveTriName) exit
    end do
  end subroutine resoudStructure

  function evaluateLine(dicoIp,triName)
    type(c_ptr),intent(in)  :: dicoIp
    character*12,intent(in) :: triName
    logical                 :: evaluateLine

    type(c_ptr)  :: triIp
    character*12 :: eqName,saveEqName
    real         :: factor,val,res

    evaluateLine = .true.
    triIp = dicoIp
    call LCMSIX(triIp,triName,1)
    eqName = ' '
    call LCMNXT(triIp,eqName)
    saveEqName = eqName
    do !boucle sur toutes les equations de triName
       if (eqName=='value') then
          call LCMNXT(triIp,eqName)
          if (eqName==saveEqName) exit
       end if
       call getValueOf(eqName,dicoIp,val)
       if (val<0.) then
          evaluateLine = .false.
       else
          call LCMGET(triIp,eqName,factor)
          res = val * factor
          call LCMGET(triIp,'value',val)
          if (val<0.) then
             call LCMPUT(triIp,'value',1,2,res)
          else if (abs(val-res)>geomPrec) then
             call XABORT("G2S: error, incoherent geometry(3)")
          end if
       end if
       call LCMNXT(triIp,eqName)
       if (eqName==saveEqName) exit 
    end do
  end function evaluateLine

  subroutine prepareCentresTri(lCentreIp,lcIp,szLc,sidesIp)
    type(c_ptr),intent(in) :: lCentreIp,lcIp,sidesIp
    integer,intent(in) :: szLc

    character*12 :: triName
    type(c_ptr)  :: triIp
    integer      :: i,j,k,lenMix,ind,iTri,lg,ty
    real         :: dsidex,dsidey,xi,yj
    integer,dimension(40) :: sv
    real,dimension(:),allocatable         :: xx,yy
    integer,dimension(:),allocatable      :: mix,turn,merg
    character*12,dimension(:),allocatable :: cell

    triName = ' '
    do k = 1,szLc
       !tavail sur la kieme cellule
       call LCMNXT(lcIp,triName)
       triIp=LCMGID(lcIp,triName)
       call LCMGET(triIp,'STATE-VECTOR',sv)
       call LCMLEN(triIp,'MIX         ',lenMix,ty)
       if ((lenMix==1).or.(sv(9)==0)) cycle !pas d'info a retirer
       allocate(mix(lenMix),turn(lenMix),merg(lenMix), &
            xx(lenMix),yy(lenMix),cell(sv(9)))
       call LCMGET(triIp,'MIX         ',mix)
       call LCMGET(triIp,'TURN        ',turn)
       call LCMLEN(triIp,'MERGE       ',lg,ty)
       if (lg/=0) then
          call LCMGET(triIp,'MERGE       ',merg)
       else
          merg = (/(i,i=1,lenMix)/)
       end if
       call LCMGTC(triIp,'CELL        ',12,sv(9),cell)
       call LCMGET(sidesIp,cell(1),dsidex)
       dsidex = 0.5*dsidex
       dsidey = sqrt_3f*dsidex
       iTri = T_ST60
       if (triName=='root') call LCMGET(triIp,'ITRI        ',iTri)
       select case(iTri)
       case(T_ST60) !triangle normal ou 'root en ST60'
          ind = 0
          yj = (1-sv(4))*0.5*dsidey
          do j = 1,sv(4)
             xi = (j-sv(3))*dsidex
             do i = 1,2*(sv(3)-j)+1
                ind = ind + 1
                xx(ind) = xi
                yy(ind) = yj
                xi = xi + dsidex
             end do
             yj = yj + dsidey
          end do
       case(T_S30)
          ind = 0
          yj = 0.5*dsidey
          do j = 1,sv(4)
             xi = (3*j-2)*dsidex
             do i = 1,sv(3)-3*(j-1)
                ind = ind + 1
                xx(ind) = xi
                yy(ind) = yj
                xi = xi + dsidex
             end do
             yj = yj + dsidey
          end do
       case(T_SA60)
          ind = 0
          yj = -(2*sv(4)-1)*0.5*dsidey
          do j = 1,sv(4)
             xi = (3*(sv(4)-j)+1)*dsidex
             do i = 1,3*j-2+mod(sv(3)-1,3)
                ind = ind + 1
                xx(ind) = xi
                yy(ind) = yj
                xi = xi + dsidex
             end do
             yj = yj + dsidey
          end do
          yj = 0.5*dsidey
          do j = 1,sv(4)
             xi = (3*j-2)*dsidex
             do i = 1,sv(3)-3*(j-1)
                ind = ind + 1
                xx(ind) = xi
                yy(ind) = yj
                xi = xi + dsidex
             end do
             yj = yj + dsidey
          end do
       case(T_Complete)
          ind = 0
          yj = (1-sv(4))*0.5*dsidey
          do j = 1,sv(4)
             xi = (j-sv(3))*dsidex
             do i = 1,2*sv(3)
                ind = ind + 1
                xx(ind) = xi
                yy(ind) = yj
                xi = xi + dsidex
             end do
             yj = yj + dsidey
          end do
       end select
       do i = 1,lenMix
          call addInRes(lCentreIp,cell(-mix(i)),triName,xx(i), &
               yy(i),turn(i),i,merg(i))
       end do
       deallocate(mix,turn,merg,xx,yy,cell)
    end do
  end subroutine prepareCentresTri

  subroutine compileResutatsTri(geoIp,lcIp,szLc,centreCalculesIp,sidesIp)
    type(c_ptr),intent(in) :: geoIp,lcIp,centreCalculesIp,sidesIp
    integer,intent(in) :: szLc

    type(c_ptr)  :: lMinMaxIp,sidexyIp
    integer      :: i,iTri
    character*12 :: cellName
    real         :: val
    integer,dimension(40) :: sv

    call LCMOP(lMinMaxIp,'minAndMaxXY ',0,1,0)
    call LCMPUT(lMinMaxIp,'minX        ',1,2,0.)
    call LCMPUT(lMinMaxIp,'maxX        ',1,2,0.)
    call LCMPUT(lMinMaxIp,'minY        ',1,2,0.)
    call LCMPUT(lMinMaxIp,'maxY        ',1,2,0.)
    
    call LCMOP(sidexyIp,'sidesXY     ',0,1,0)
    cellName = ' '
    do i = 1,szLc
       call LCMNXT(lcIp,cellName)
       call LCMGET(sidesIp,cellName,val)
       call LCMPUT(sidexyIp,'x'//cellName(1:11),1,2,val)
       call LCMPUT(sidexyIp,'y'//cellName(1:11),1,2,val)
       if (cellName=='root') then !modification de xroot et yroot si besoin
          call LCMGET(geoIp,'STATE-VECTOR',sv)
          call LCMGET(geoIp,'ITRI        ',iTri)
          select case(iTri)
          case(T_Complete)
             val=val*sv(4)/real(sv(3))
             call LCMPUT(sidexyIp,'y'//cellName(1:11),1,2,val)
          case(T_S30,T_SA60)
             val=val/real(sv(3))
             call LCMPUT(sidexyIp,'x'//cellName(1:11),1,2,val)
             val=val*int((sv(3)+1)/2)
             call LCMPUT(sidexyIp,'y'//cellName(1:11),1,2,val)
          end select
       end if   
    end do

    call compileResutats(geoIp,centreCalculesIp,sidexyIp,lMinMaxIp)

    call LCMCL(lMinMaxIp,2)
    call LCMCL(sidexyIp,2)
  end subroutine compileResutatsTri

end module pretraitement
