!
!-----------------------------------------------------------------------
!
!Purpose:
! get an array of character variables from a LCM table component.
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
! name    character name of the LCM node.
! leng    length of each character variable in the array carr.
! nlin    dimension of array carr.
!
!Parameters: output
! carr    array of character variables.
!
!-----------------------------------------------------------------------
!
subroutine LCMGTC(iplcm,name,leng,nlin,carr)
   use GANLIB
   !----
   !  Subroutine arguments
   !----
   type(c_ptr) iplcm
   integer :: leng,nlin
   character*(*) :: name,carr(nlin)
   !----
   !  Local variables
   !----
   character(len=13)  :: fmt
   character(len=12) :: text12
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   allocate(ibase(nlin*(leng+3)/4))
   !----
   !  Read from LCM object
   !----
   call LCMLEN(iplcm,name,ilong,itylcm)
   if(ilong == 0) then
      call LCMLIB(iplcm)
      text12=name
      call XABORT('LCMGTC: record '//text12//' not found.')
   endif
   call LCMGET(iplcm,name,ibase)
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
end subroutine LCMGTC
