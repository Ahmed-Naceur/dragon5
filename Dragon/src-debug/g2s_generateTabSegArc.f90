!
!-----------------------------------------------------------------------
!
!Purpose:
! Fill the TabSegArc structure with information recovered in a
! surfacic file.
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
!-----------------------------------------------------------------------
!
module generTabSegArc
  use constUtiles
  use segArc
  use boundCond
  use SALGET_FUNS_MOD
  use precision_and_kinds, only : pdb
  implicit none

contains

  subroutine generateTabSegArc(ipSal,sizeSA,nbNode,nbCLP,nbFlux,merg,name_geom)
    integer,intent(inout) :: ipSal
    integer,intent(in) :: nbNode,sizeSA
    integer,intent(out) :: nbCLP,nbFlux
    integer,dimension(nbNode),intent(out) :: merg
    character(len=12),intent(out) :: name_geom

    integer, parameter :: n_datain=25, n_datare=20
    integer, dimension (n_datain) :: datain
    real,    dimension (n_datare) :: datare
    real(pdb),    dimension (n_datare) :: datade
    integer, parameter :: fout0=6
    integer :: type,nber,prec,elem,i,nbMacro
    integer, parameter, dimension(0:4) :: read_bc_len=(/1,1,2,3,3/)
    ! internal : albedo
    ! vacuum surface : albedo
    ! specular reflexion : none
    ! translation :     tx ty (t=translation vector)
    ! rotation :        cx cy cos(theta) sin(theta) theta
    !                   (c= center,theta= axis angle)
    ! axial symmetry :  cx cy cos(theta) sin(theta) theta
    !                   (c= center,theta= axis angle)
    ! central symetry : cx cy (c= center)
    integer, allocatable, dimension(:) :: medium,iflux

    call SALGET(datain,6,ipSal,fout0,'dimensions for geometry')
    if(nbNode /= datain(3)) call XABORT('g2s_generateTabSegArc: nbNode error')
    if(sizeSA /= datain(4)) call XABORT('g2s_generateTabSegArc: sizeSA error')
    nbMacro=datain(5)
    nbFlux=datain(6)
    print *,'nbNode=',nbNode,' sizeSA=',sizeSA
    call SALGET(datain,3,ipSal,fout0,'index kndex prec')
    prec=datain(3)
    call SALGET(datare,1,ipSal,fout0,'eps')
    call SALGET(merg,nbNode,ipSal,fout0,'flux index per node')
    call SALGET(name_geom,ipSal,fout0,'names of macros')
  allocate(iflux(nbFlux),stat=alloc_ok)
  if (alloc_ok /= 0) call XABORT("G2S: generateTabSegArc(1) => allocation pb")
    call SALGET(iflux,nbFlux,ipSal,fout0,'macro order number per flux region')
    deallocate(iflux)
    do elem=1,sizeSA
       call SALGET(datain,3,ipSal,fout0,'integer descriptors')
       type=datain(1)
       tabSegArc(elem)%typ=type
       tabSegArc(elem)%noded=datain(2)
       tabSegArc(elem)%nodeg=datain(3)
       select case (type)
       case (1)
          nber=4
       case (2)
          nber=3
       case (3)
          nber=5
       case default
          write(fout0,'(1x,''==> sal126: unknown type '',i3)')type
          call xabort('g2_generateTabSegArc: unknown type')
       end select
       call SALGET(datade,nber,ipSal,fout0,prec,'real descriptors')
       tabSegArc(elem)%x=datade(1)
       tabSegArc(elem)%y=datade(2)
       select case (type)
       case (1)
          tabSegArc(elem)%dx=datade(1)+datade(3)
          tabSegArc(elem)%dy=datade(2)+datade(4)
       case (2)
          tabSegArc(elem)%r=datade(3)
          tabSegArc(elem)%a=0.0
          tabSegArc(elem)%b=0.0
       case (3)
          tabSegArc(elem)%r=datade(3)
          tabSegArc(elem)%a=datade(4)/rad2deg
          tabSegArc(elem)%b=(datade(4)+datade(5))/rad2deg
       end select
    enddo
    call SALGET(datain,3,ipSal,fout0,'general bc data')
    nbCLP=datain(2)
    call SALGET(datade(1),ipSal,fout0,prec,'general albedo')
    do i=1,nbCLP
       call SALGET(datain,2,ipSal,fout0,'specific bc: type nber')
       type=datain(1)
       nber=datain(2)
       SALbCDataTab(i)%SALtype=type
       SALbCDataTab(i)%nber=nber
    allocate(SALbCDataTab(i)%elemNb(nber),stat=alloc_ok)
    if (alloc_ok /= 0) call XABORT("G2S: generateTabSegArc(3) => allocation pb")
       call SALGET(SALbCDataTab(i)%elemNb,nber,ipSal,fout0,'bc elements')
       ! read bc motion
       call SALGET(datade,read_bc_len(type),ipSal,fout0,prec,'data for specific bc condition')
       select case(type)
       case(0,1)
          SALbCDataTab(i)%albedo=real(datade(1))
       case(2)
          SALbCDataTab(i)%tx=real(datade(1))
          SALbCDataTab(i)%ty=real(datade(2))
       case(3,4)
          SALbCDataTab(i)%cx=real(datade(1))
          SALbCDataTab(i)%cy=real(datade(2))
          SALbCDataTab(i)%angle=real(datade(3))
       end select
    enddo
  allocate(medium(nbNode),stat=alloc_ok)
  if (alloc_ok /= 0) call XABORT("G2S: generateTabSegArc(4) => allocation pb")
    call SALGET(medium,nbNode,ipSal,fout0,'media per node')
    do i=1,sizeSA
       if(tabSegArc(i)%nodeg>0) tabSegArc(i)%neutronicMixg=medium(tabSegArc(i)%nodeg)
       if(tabSegArc(i)%noded>0) tabSegArc(i)%neutronicMixd=medium(tabSegArc(i)%noded)
    enddo
    deallocate(medium)
    print *,'leaving g2s_generateTabSegArc'
  end subroutine generateTabSegArc
end module generTabSegArc
