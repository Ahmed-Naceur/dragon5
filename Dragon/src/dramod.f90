!
!-----------------------------------------------------------------------
!
!Purpose:
! Dispatch to a calculation module in DRAGON. ANSI-C interoperable.
!
!Copyright:
! Copyright (C) 2009 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s): A. Hebert
!
!-----------------------------------------------------------------------
!
integer(c_int) function dramod(cmodul, nentry, hentry, ientry, jentry, &
               kentry, hparam_c) bind(c)
!
   use GANLIB
   implicit none
!----
!  subroutine arguments
!----
   character(kind=c_char), dimension(*) :: cmodul
   integer(c_int), value :: nentry 
   character(kind=c_char), dimension(13,*) :: hentry
   integer(c_int), dimension(nentry) :: ientry, jentry
   type(c_ptr), dimension(nentry) :: kentry
   character(kind=c_char), dimension(73,*) :: hparam_c
!----
!  local variables
!----
   integer :: i, ier
   character :: hmodul*12, hsmg*131, hparam*72
   character(len=12), allocatable :: hentry_f(:)
   type FIL_file_array
     type(FIL_file), pointer :: my_file
   end type FIL_file_array
   type(FIL_file_array), pointer :: my_file_array(:)
   integer, external :: KDRDRV
!
   allocate(hentry_f(nentry),my_file_array(nentry))
   call STRFIL(hmodul, cmodul)
   do i=1,nentry
      call STRFIL(hentry_f(i), hentry(1,i))
      if((ientry(i) >= 3).and.(ientry(i) <= 5)) then
!        open a Fortran file.
         call STRFIL(hparam, hparam_c(1,i))
         my_file_array(i)%my_file=>FILOPN(hparam,jentry(i),ientry(i)-1,0)
         if(.not.associated(my_file_array(i)%my_file)) then
            write(hsmg,'(29hdramod: unable to open file '',a12,2h''.)') hentry_f(i)
            call XABORT(hsmg)
         endif
         kentry(i)=c_loc(my_file_array(i)%my_file)
      endif
   enddo
!  ----------------------------------------------------------
   dramod=KDRDRV(hmodul,nentry,hentry_f,ientry,jentry,kentry)
!  ----------------------------------------------------------
   do i=1,nentry
      if(jentry(i) == -2) then
!        destroy a LCM object or a Fortran file.
         if(ientry(i) <= 2) then
            call LCMCL(kentry(i),2)
            kentry(i)=c_null_ptr
         else if((ientry(i) >= 3).and.(ientry(i) <= 5)) then
            ier=FILCLS(my_file_array(i)%my_file,2)
            if(ier < 0) then
               write(hsmg,'(32hdramod: unable to destroy file '',a12,2h''.)') hentry_f(i)
               call XABORT(hsmg)
            endif
            kentry(i)=c_null_ptr
         endif
      else
!        close a Fortran file.
         if((ientry(i) >= 3).and.(ientry(i) <= 5)) then
            ier=FILCLS(my_file_array(i)%my_file,1)
            if(ier < 0) then
               write(hsmg,'(30hdramod: unable to close file '',a12,2h''.)') hentry_f(i)
               call XABORT(hsmg)
            endif
         endif
      endif
   enddo
   deallocate(my_file_array,hentry_f)
   flush(6)
   return
end function dramod
