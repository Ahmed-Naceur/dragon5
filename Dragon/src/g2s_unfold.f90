!
!-----------------------------------------------------------------------
!
!Purpose:
! Unfold the geometry.
!
!Copyright:
! Copyright (C) 2020 Ecole Polytechnique de Montreal.
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s):
! A. Hebert
!
!-----------------------------------------------------------------------
!
subroutine g2s_unfold(geoIp,impx)
  use GANLIB
  use constType
  type(c_ptr),intent(in):: geoIp
  integer,parameter :: nstate=40
  integer,dimension(nstate) :: st
  integer,allocatable,dimension(:) :: idp,ind1,ind2
  !
  call LCMGET(geoIp,'STATE-VECTOR',st)
  select case(st(1))
  case(G_Hex)
    call LCMLIB(geoIp)
    call LCMGET(geoIp,'IHEX        ',iHex)
    lxold=st(6)
    if((iHex /= 9).and.(lxold > 1)) then
      call LCMLEN(geoIp,'TURN        ',ilong,itylcm)
      if(ilong > 0) call XABORT('g2s_unfold: TURN not supported.')
      ! caution: HEXCEL cells are not rotated according to symmetries in BIVALL
      maxpts=12*lxold
      allocate(idp(maxpts))
      call BIVALL(maxpts,iHex,lxold,lx,idp)
      if(impx > 0) write(*,*) 'g2s_unfold: nb of hexagons=',lxold,'-->',lx
      allocate(ind1(lxold),ind2(lx))
      call LCMGET(geoIp,'MIX         ',ind1)
      do i=1,lx
        ind2(i)=ind1(idp(i))
      enddo
      call LCMPUT(geoIp,'MIX         ',lx,1,ind2)
      call LCMLEN(geoIp,'MERGE       ',ilong,itylcm)
      if(ilong > 0) then
        call LCMGET(geoIp,'MERGE       ',ind1)
        do i=1,lx
          ind2(i)=ind1(idp(i))
        enddo
        call LCMPUT(geoIp,'MERGE       ',lx,1,ind2)
      endif
      deallocate(ind2,ind1,idp)
      st(3)=lx
      st(6)=lx
      call LCMPUT(geoIp,'STATE-VECTOR',nstate,1,st)
      iHex=9
      call LCMPUT(geoIp,'IHEX        ',1,1,iHex)
    endif
  end select
end subroutine g2s_unfold
