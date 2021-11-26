!
!-----------------------------------------------------------------------
!
!Purpose:
! Generate a dataset for use in a Monte Carlo code.
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
!Parameters: input/output
!  NENTRY : NUMBER OF LINKED LISTS AND FILES USED BY THE MODULE.
!  HENTRY : CHARACTER*12 NAME OF EACH LINKED LIST OR FILE.
!  IENTRY : =1 LINKED LIST; =2 XSM FILE; =3 SEQUENTIAL BINARY FILE;
!           =4 SEQUENTIAL ASCII FILE; =5 DIRECT ACCESS FILE.
!  JENTRY : =0 THE LINKED LIST OR FILE IS CREATED;
!           =1 THE LINKED LIST OR FILE IS OPEN FOR MODIFICATIONS;
!           =2 THE LINKED LIST OR FILE IS OPEN IN READ-ONLY MODE.
!  KENTRY : FILE UNIT NUMBER OR LINKED LIST ADDRESS.
!
!-----------------------------------------------------------------------
!
subroutine G2MC(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
  use SALGET_FUNS_MOD
  use celluleBase
  use cellulePlaced
  use boundCond
  use ptNodes
  use pretraitement
  use derivedPSPLOT
  use monteCarlo
  use track
  use segArc
  use GANLIB
  use generTabSegArc

  implicit none

  integer NENTRY
  integer IENTRY,JENTRY
  type(c_ptr) KENTRY
  character*12 HENTRY
  dimension IENTRY(*),JENTRY(*),KENTRY(*),HENTRY(*)


  integer,parameter :: dimTabCelluleBase = 20000
  integer,parameter :: dimTabSegArc = 100000

  type(c_ptr) :: ipGeo,ipGeo_1
  integer     :: ipMC,ipSal,ipPs,sizeB,sizeP,sizeSA,nbNode,nbCLP,nbFlux
  logical,parameter :: drawMix = .false. ! set to .true. to draw the mix numbers
  real,parameter,dimension(2) :: zoomx = (/ 0.0, 1.0 /) ! no x zoom on postscript plot
  real,parameter,dimension(2) :: zoomy = (/ 0.0, 1.0 /) ! no y zoom on postscript plot
  integer      :: lgMaxGig=0
  integer,dimension(10) :: datain
  integer,allocatable,dimension(:) :: merg
  character(len=12) :: name_geom

  if ((nentry == 2).and.(IENTRY(2) == 4)) then
     !generating Monte-Carlo file from Sal file
     ipMC   = FILUNIT(KENTRY(1)) ! Monte-Carlo file generated
     ipSal  = FILUNIT(KENTRY(2)) ! input Salome file (surfacic elements)
     ipPs   = -1         ! no postscript file
     ipGeo_1= c_null_ptr ! no geometry read
     ! check that second argumnet is file to write
     ! then the tracking object
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
        call XABORT('G2MC: a new ascii file expected at LHS for containing MC info')
     if ((IENTRY(2) /= 4) .or. (JENTRY(2) /= 2)) &
        call XABORT('G2MC: read-only ascii file expected at RHS with surfacic elements')
  else if ((nentry == 2).and.(IENTRY(2) <= 2)) then
     !generating Monte-Carlo file from LCM geometry
     ipPs   = -1        ! no postscript file
     ipMC   = FILUNIT(KENTRY(1)) ! Monte-Carlo file generated
     ipGeo_1= KENTRY(2) ! geometry read
     ! check that second argumnet is file to write
     ! then the tracking object
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
        call XABORT('G2MC: a new ascii file expected at LHS for containing MC info')
  else if ((nentry == 3).and.(IENTRY(1) == 4)) then
     !generating Sal file and ps file from LCM geometry
     ipMC   = FILUNIT(KENTRY(1)) ! Monte-Carlo file generated
     ipPs   = FILUNIT(KENTRY(2)) ! psfile generated
     g_psp_isEpsFile = (index(HENTRY(2),'.eps')/=0) !is it an eps file ?
     ipGeo_1= KENTRY(3) ! geometry read
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
        call XABORT('G2MC: a new ascii file expected at LHS for containing MC info')
     if ((IENTRY(2) /= 4) .or. (JENTRY(2) /= 0)) &
        call XABORT('G2MC: a new file was expected for the postscript file')
  else
     call XABORT('G2MC: you must provide 2 or 3 arguments')
  end if
  
  sizeB = 0   !cellules de base
  sizeP = 0   !cellules placees
  sizeSA = 0  !elements geometriques
  if(c_associated(ipGeo_1)) then
     ! copy the input geometric object
     call lcmop(ipGeo,'geom_copy',0,1,0)
     call lcmequ(ipGeo_1,ipGeo)

     !initialisation des differents tableaux
     call initializeData(dimTabCelluleBase,dimTabSegArc)

     !unfold the geometry
     call g2s_unfold(ipGeo,0)

     !pretraitement des donnees lues (remplace la partie python)
     !+completion des cellules de base et remplissage du tableau
     !des cellules placees
     call prepareData(ipGeo,sizeB,sizeP,lgMaxGig)

     !en sortie, toutes les cellules de base ont tous leurs
     !champs remplis, et le tableau des cellules placees est pret

     !eclatement des cellules
     call splitCells(sizeP,sizeSA)

     !creation de nouveaux segments aux interfaces des cellules
     !et elimination des doublons
     call addSegsAndClean(sizeSA)

     !prise en compte des conditions aux limites
     call appliBoundariConditions(ipGeo,sizeSA,nbCLP)

     !calcul des nodes delimites par les elements
     allocate(merg(dimTabCelluleBase),stat=alloc_ok)
     if (alloc_ok /= 0) call XABORT("G2S: g2s_g2mc(1) => allocation pb(1)")
     call createNodes(sizeSA,dimTabCelluleBase,nbNode,merg)
     if(sizeSA > dimTabSegArc) call XABORT('g2s_g2mc: sizeSA overflow')
  else
     if (JENTRY(nentry) == 0) call XABORT('G2M: an existing Salomon file is expected')
     !initialisation de TabSegArc
     call SALGET(datain,4,ipSal,0,'dimensions for geometry')
     nbNode=datain(3)
     sizeSA=datain(4)
     rewind(ipSal)
     allocate(tabSegArc(sizeSA))
     call initializebCData()  
     allocate(merg(nbNode),stat=alloc_ok)
     if (alloc_ok /= 0) call XABORT("G2S: g2s_g2mc => allocation pb")
     call generateTabSegArc(ipSal,sizeSA,nbNode,nbCLP,nbFlux,merg,name_geom)
  endif
  deallocate(merg)

  !impression des segArc charges
  if (ipPs /= -1) call drawSegArc(ipPs,sizeSA,.true.,drawMix,zoomx,zoomy)

  !creation du fichier de commande Monte-Carlo
  if (index(HENTRY(1),'.tp')/=0) then
     ! generate a Tripoli4 datafile
     call generateTripoliFile(ipMC,sizeSA,nbNode)
  else if (index(HENTRY(1),'.sp')/=0) then
     ! generate a Serpent datafile
     call generateSerpentFile(ipMC,sizeSA,nbNode)
  else
     ! generate a MCNP datafile
     call generateMCNPFile(ipMC,sizeSA,nbNode)
  end if

  print *,"  At end of G2MC:","    ",nbNode,"volumes"

  !liberation de la memoire allouee
  call destroyData(sizeB,sizeP)
end subroutine G2MC
