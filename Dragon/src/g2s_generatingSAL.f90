!
!-----------------------------------------------------------------------
!
!Purpose:
! Generate the surfacic geometry ascii file.
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
! Attention, comme SAL est toujours en phase de developpement, des
! modifications ont ete apportees a la routine "calculTypgeo" pour bien prendre
! en compte les condition limites de reflexion. En effet, celle-ci sont
! traitees comme des symetries axiales, dont l'axe est le bord de la geometrie.
! Cependant, une evolution previsible de SAL devrait rendre obsolette cette
! assimilation, et permettre une prise en compte directe des reflexion.
! \\
! fonctions:
!  - generateSALFile : creation du fichier de donnees SAL
!  - calculTypgeo : determination des donnees typgeo et nbfold
!  - calculDefaultCl : dertermination de la meilleur condition limite par
!                      defaut
!
!-----------------------------------------------------------------------
!
module generSAL
    use boundCond
    use cellulePlaced
    use constType
    use constUtiles
    use segArc

    implicit none
 
contains

  subroutine generateSALFile(fileNbr,szSA,nbNode,nbCLP,nbFlux,merg,name_geom)
   integer,intent(in) :: fileNbr,szSA,nbNode,nbCLP,nbFlux
    integer,dimension(nbNode),intent(in) :: merg
    character(len=*),intent(in) :: name_geom

    type(t_segArc)                   :: sa
    integer                          :: i,j,nmoins,nplus
    integer                          :: typgeo,nbfold,defautCl
    integer,dimension(:),allocatable :: milTab
    integer,dimension(:),allocatable :: tmpTab
    real                             :: cx,cy,tx,ty,delta,albedo
    character(len=4)                 :: strDCL

    rewind(fileNbr)
    !preparer les donnees de CL pour SAL
    call prepareSALBCData(szSA,nbCLP)

    write(fileNbr,'(a5)')   'BEGIN'

    write(fileNbr,'(/a24)')  'DEFINE DOMAINE'
    write(fileNbr,'(a24/)') '=============='

    write(fileNbr,'(a28/)') '1.main dimensions:'
    write(fileNbr,'(a32)')  '*typgeo  nbfold  nbnode  nbelem '
    call calculTypgeo(typgeo,nbfold)
    write(fileNbr,'(10'//formati//')') typgeo,nbfold,nbNode,szSA,1,nbFlux

    write(fileNbr,'(/a37)') '2.impression and precision:'
    write(fileNbr,'(/a21)') '*index   kndex   prec'
    write(fileNbr,'(10'//formati//')') 0,0,1

    write(fileNbr,'(/a40)') '3.precision of geometry data:'
    write(fileNbr,'(/a4)')  '*eps'
    write(fileNbr,'('//formatr//')') gSALeps

    write(fileNbr,'(/a58)') '4.flux region number per geometry region (mesh):'
    write(fileNbr,'(/a6)')  '*merge'
    write(fileNbr,'(10'//formati//')') merg

    write(fileNbr,'(/a29)') '5.name of geometry:'
    write(fileNbr,'(/a11)')  '*macro_name'
    write(fileNbr,'(4(3x,a12))') name_geom

    write(fileNbr,'(/a46)') '6.macro order number per flux region:'
    write(fileNbr,'(/a14)')  '*macro_indices'
    write(fileNbr,'(10'//formati//')') (1,i=1,nbFlux)

    write(fileNbr,'(/a57)') '7.read integer and real data for each elements:'
    do i = 1,szSA
       sa = tabSegArc(i)
       cx     = real(sa%x)
       cy     = real(sa%y)
       if (sa%typ==tseg) then
          nmoins = sa%noded
          nplus  = sa%nodeg
          tx     = real(sa%dx-sa%x)
          ty     = real(sa%dy-sa%y)
          delta  = 0.
       else
          nmoins = sa%nodeg
          nplus  = sa%noded
          tx     = real(sa%r)
          if (sa%typ==tarc) then
             ty     = real(sa%a*rad2deg)
             if (sa%b>sa%a) then
                delta  = real((sa%b-sa%a)*rad2deg)
             else
                delta  = real((sa%b-sa%a)*rad2deg+360.d0)
             end if
          else
             ty     = 0.
             delta  = 0.
          end if
       end if
       write(fileNbr,*)
       write(fileNbr,*) 'elem =',i 
       write(fileNbr,'(a22)')    '*type    node-   node+' 
       write(fileNbr,'(10'//formati//')')    tabSegArc(i)%typ,nmoins,nplus
       write(fileNbr,'(a63)')    '*cx            cy            ex or R       &
            &ey or theta1  theta2' 
       write(fileNbr,'(5'//formatr//')')  cx,cy,tx,ty,delta
    end do

    write(fileNbr,'(/a63)') '8.read integer and real data for boundary &
         &conditions:'
    !test de la condition aux limites par defaut en fonction du type de geo
    call calculDefaultCl(defautCl,albedo,strDCL)
    write(fileNbr,'(/a40)') '*defaul  nbbcda  allsur  divsur  ndivsur'
    write(fileNbr,'(10'//formati//')') defautCl,nbCLP,0,0,0
    write(fileNbr,'(/a24)') 'DEFAULT = ' // strDCL
    write(fileNbr,'(a24)')  '=============='
    write(fileNbr,'(/a17)') '*albedo  deltasur'
    write(fileNbr,'(5'//formatr//')') albedo,0.0
    do i = 1,nbCLP
       write(fileNbr,*)
       write(fileNbr,*) 'particular boundary condition number',i
       write(fileNbr,'(/a13)') '*type    nber'
       write(fileNbr,'(10'//formati//')') SALbCDataTab(i)%SALtype,SALbCDataTab(i)%nber
     allocate(tmpTab(SALbCDataTab(i)%nber),stat=alloc_ok)
     if (alloc_ok /= 0) call XABORT("G2S: generateSALFile(1) => allocation pb")
       do j = 1,SALbCDataTab(i)%nber
          tmpTab(j) = SALbCDataTab(i)%elemNb(j)
       end do
       write(fileNbr,'(/a14)') '*elems(1,nber)'
       write(fileNbr,'(10'//formati//')') tmpTab
       deallocate(tmpTab)
       select case(SALbCDataTab(i)%SALtype)
       case(0,1)
          write(fileNbr,'(/a7)')    '*albedo'
          write(fileNbr,'(5'//formatr//')') SALbCDataTab(i)%albedo
       case(2)
          write(fileNbr,'(/a11)') '*tx      ty'
          write(fileNbr,'(5'//formatr//')') SALbCDataTab(i)%tx,SALbCDataTab(i)%ty
       case(3,4)
          write(fileNbr,'(/a22)') '*cx      cy      angle'
          write(fileNbr,'(5'//formatr//')') SALbCDataTab(i)%cx,SALbCDataTab(i)%cy &
               & ,SALbCDataTab(i)%angle
       end select
    end do

  allocate(milTab(nbNode),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: generateSALFile(2) => allocation pb")
    do i = 1,szSA
       if (tabSegArc(i)%nodeg>0) &
            & milTab(tabSegArc(i)%nodeg) = tabSegArc(i)%neutronicMixg
       if (tabSegArc(i)%noded>0) &
            & milTab(tabSegArc(i)%noded) = tabSegArc(i)%neutronicMixd
    end do
    write(fileNbr,'(/a28)') '9.medium per node:'
    write(fileNbr,'(/a11)') '*mil(nbreg)'
    write(fileNbr,'(10'//formati//')') milTab
    deallocate(milTab)

    write(fileNbr,'(/a3)') 'END'
  end subroutine generateSALFile


  subroutine calculTypgeo(typgeo,nbfold)
    integer,intent(out) :: typgeo,nbfold

    integer :: i
    integer,dimension(4) :: bc

    typgeo = 0 ; nbfold = 0
    select case(geomTyp)
    case(RecTyp)
       bc=bCData%bc(1:4)
       if ( (bc(1)==B_Pi_2)                    .and. &
            (bc(2)==B_Syme .or. bc(2)==B_Refl .or. bc(2)==B_Ssym) .and. &
            (bc(3)==B_Pi_2)                    .and. &
            (bc(4)==B_Syme .or. bc(3)==B_Refl .or. bc(3)==B_Ssym)) then
          typgeo = 11
       else if (all(bc==(/B_Pi_2,B_Tran,B_Pi_2,B_Tran/))) then
          typgeo = 10
       else if ((bc(1)==B_Diag )                   .and. &
            (bc(2)==B_Syme .or. bc(2)==B_Refl .or. bc(2)==B_Ssym) .and. &
            (bc(3)==B_Syme .or. bc(3)==B_Refl .or. bc(3)==B_Ssym) .and. &
            (bc(4)==B_Diag                   )) then
          typgeo = 7
       else if ((bc(1)==B_Syme .or. bc(1)==B_Refl .or. bc(1)==B_Ssym) .and. &
            (bc(2)==B_Syme .or. bc(2)==B_Refl .or. bc(2)==B_Ssym) .and. &
            (bc(3)==B_Syme .or. bc(3)==B_Refl .or. bc(3)==B_Ssym) .and. &
            (bc(4)==B_Syme .or. bc(3)==B_Refl .or. bc(3)==B_Ssym)) then
          typgeo = 6
       else if (all(bc==(/B_Tran,B_Tran,B_Tran,B_Tran/))) then
          typgeo = 5
       else if (((bc(1)==B_Tran).and.(bc(2)==B_Tran)) &
            & .or.((bc(3)==B_Tran).and.(bc(4)==B_Tran))) then
          typgeo = 4
       else if ( ( ((bc(1)==B_Syme).or.(bc(1)==B_Refl).or.(bc(1)==B_Ssym))   &
            & .and.((bc(2)==B_Syme).or.(bc(2)==B_Refl).or.(bc(2)==B_Ssym)) ) &
            & .or.(((bc(3)==B_Syme).or.(bc(3)==B_Refl).or.(bc(3)==B_Ssym))   &
            & .and.((bc(4)==B_Syme).or.(bc(4)==B_Refl).or.(bc(4)==B_Ssym)) ) )then
          typgeo = 3
       else if ((bc(1)==B_Pi_2).and.(bc(3)==B_Pi_2)) then
          typgeo = 2 ; nbfold = 4
       else if (   ((bc(1)==B_Syme).or.(bc(1)==B_Refl).or.(bc(1)==B_Ssym))   &
            & .and.((bc(3)==B_Syme).or.(bc(3)==B_Refl).or.(bc(3)==B_Ssym)) ) then
          typgeo = 1 ; nbfold = 4
       else if ((bc(1)==B_Diag).and.(bc(3)==B_Syme.or.bc(3)==B_Refl.or.bc(3)==B_Ssym) &
            .and. (bc(4)==B_Diag)) then
          typgeo = 1 ; nbfold = 8
       else if ((bc(3)==B_Syme).or.(bc(3)==B_Refl).or.(bc(3)==B_Ssym)) then
          typgeo = 1 ; nbfold = 2
       end if
       if (typgeo==0 .and. all((/ &
            ( bc(i)==B_Refl.or.bc(i)==B_Ssym.or.bc(i)==B_Syme.or.bc(i)==B_Diag,i=1,4)/))) &
            call XABORT("G2S: Type of boundary conditions not supported by SAL yet(1)")
    case(HexTyp)
       if (bCData%iHex==H_S30) then
          typgeo = 1 ; nbfold = 12
       end if
       if (typgeo==0 .and. bCData%bc(1)==B_Refl .or. bCData%bc(1)==B_Syme) &
            call XABORT("G2S: Type of boundary conditions not supported by SAL yet(2)")
    case(TriaTyp)
       if (bCData%iTri==T_S30) then
          typgeo = 1 ; nbfold = 12
       end if
    case(TubeTyp)
       !rien de particulier
    end select
  end subroutine calculTypgeo

  subroutine calculDefaultCl(defautCl,albedo,strDCL)
    integer,intent(out)     :: defautCl
    real,intent(out)        :: albedo
    character*4,intent(out) :: strDCL

    defautCl = 0
    albedo   = 0.0
    strDCL   = 'VOID'
    select case(geomTyp)
    case(HexTyp) 
       select case(bCData%bc(1))
       case(B_Void,B_Syme)
          !rien a faire
       case(B_Albe)
          albedo = real(bCData%albedo(1))
          strDCL = 'ALBE'
       case(B_Refl)
          defautCl = 1
          strDCL = 'REFL'
       case default
          call XABORT("G2S: type of boundary condition not allowed")
       end select
    case(TriaTyp)
       if (bCData%iTri==T_S30 .or. bCData%iTri==T_SA60) then
          select case(bCData%bc(1))
          case(B_Void,B_Syme)
             !rien a faire
          case(B_Albe)
             albedo = real(bCData%albedo(1))
             strDCL = 'ALBE'
          case(B_Refl)
             defautCl = 1
             strDCL = 'REFL'
          case default
             call XABORT("G2S: type of boundary condition not allowed")
          end select
       end if
    case(TubeTyp)
       !ATTENTION, l'indice de la cl dans le jdd est 2 et non 1
       select case(bCData%bc(2))
       case(B_Void)
          !rien a faire
       case(B_Albe)
          albedo = real(bCData%albedo(2))
          strDCL = 'ALBE'
       case(B_Refl)
          defautCl = 1
          strDCL = 'REFL'
       case default
          call XABORT("G2S: type of boundary condition not allowed")
       end select
    end select
  end subroutine calculDefaultCl
end module generSAL
