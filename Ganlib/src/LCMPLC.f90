!
!-----------------------------------------------------------------------
!
!Purpose:
! put an array of character variables in a LCM list component.
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
! carr    array of character variables.
!
!-----------------------------------------------------------------------
!
subroutine LCMPLC(iplcm,ipos,leng,nlin,carr)
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
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   !----
   !  Define format and dimensions for conversion
   !----
   allocate(ibase(nlin*(leng+3)/4))
   n=leng/4
   m=mod(leng,4)
   if(n > 0)write(fmt(2:9),'(i6,''a4'')') n
   if(m > 0)then
      if(n.gt.0)fmt(10:10)=','
      write(fmt(11:12),'(''a'',i1)') m
      n=n+1
   endif
   !----
   !  Convert strings to integers
   !----
   do i=1,nlin
      read(carr(i)(1:leng),fmt) (ibase((i-1)*n+j),j=1,n)
   enddo
   !----
   !  Write to LCM object
   !----
   call LCMPDL(iplcm,ipos,n*nlin,3,ibase)
   deallocate(ibase)
end subroutine LCMPLC
