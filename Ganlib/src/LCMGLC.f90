!
!-----------------------------------------------------------------------
!
!Purpose:
! get an array of character variables from a LCM list component.
!
!Copyright:
! Copyright (C) 2010 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version
!
!Author(s): A. Hebert
!
!Parameters: input
! iplcm   pointer to the LCM object.
! ipos    list index.
! leng    length of each character variable in the array carr.
! nlin    dimension of array carr.
!
!Parameters: output
! carr    array of character variables.
!
!-----------------------------------------------------------------------
!
subroutine LCMGLC(iplcm,ipos,leng,nlin,carr)
   use GANLIB
   !----
   !  Subroutine arguments
   !----
   type(c_ptr) iplcm
   integer :: ipos,leng,nlin
   character*(*) :: carr(nlin)
   !----
   !  Local variables
   !----
   character(len=13) :: fmt
   character(len=131) :: hsmg
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   allocate(ibase(nlin*(leng+3)/4))
   !----
   !  Read from LCM object
   !----
   call lcmlel(iplcm,ipos,ilong,itylcm)
   if(ilong == 0) then
      call LCMLIB(iplcm)
      write(hsmg,'(8hLCMGLC: ,i5,21h-th record not found.)') ipos
      call XABORT(hsmg)
   endif
   call LCMGDL(iplcm,ipos,ibase)
   !----
   !  Define format and dimensions for conversion
   !----
   n=leng/4
   m=mod(leng,4)
   if(n > 0)write(fmt(2:9),'(i6,''a4'')') n
   if(m > 0)then
      if(n.gt.0)fmt(10:10)=','
      write(fmt(11:12),'(''a'',i1)') m
      n=n+1
   endif
   !----
   !  Convert integers to strings
   !----
   do i=1,nlin
      write(carr(i)(1:leng),fmt) (ibase((i-1)*n+j),j=1,n)
   enddo
   deallocate(ibase)
end subroutine LCMGLC
