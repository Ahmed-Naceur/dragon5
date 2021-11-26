!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for lcm -- part 1.
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
module LCMAUX
contains
subroutine STRCUT(name1, name2)
   ! transform a Fortran string into a C null-terminated string
   use, intrinsic :: iso_c_binding
   character(kind=c_char), dimension(*) :: name1
   character(len=*) :: name2
   integer :: ilong
   interface
      subroutine strcut_c (s, ct, n) bind(c)
         use, intrinsic :: iso_c_binding
         character(kind=c_char), dimension(*) :: s, ct
         integer(c_int), value :: n
      end subroutine strcut_c
   end interface
   ilong=len(name2)
   call strcut_c(name1, name2, ilong)
end subroutine STRCUT
!
subroutine STRFIL(name1, name2)
   ! transform a C null-terminated string into a Fortran string
   use, intrinsic :: iso_c_binding
   character(len=*) :: name1
   character(kind=c_char), dimension(*) :: name2
   integer :: ilong
   interface
      subroutine strfil_c (s, ct, n) bind(c)
         use, intrinsic :: iso_c_binding
         character(kind=c_char), dimension(*) :: s, ct
         integer(c_int), value :: n
      end subroutine strfil_c
   end interface
   ilong=len(name1)
   call strfil_c(name1, name2, ilong)
end subroutine STRFIL
!
function LCMARA(ilong)
   ! allocate an array of length ilong and return a c_ptr pointer
   use, intrinsic :: iso_c_binding
   type(c_ptr) LCMARA
   integer :: ilong
   interface
      function setara_c (length) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) setara_c
         integer(c_int), value :: length
      end function setara_c
   end interface
   LCMARA=setara_c(ilong)
end function LCMARA
!
subroutine LCMDRD(ipdata)
   ! deallocate an array allocated by LCMARA
   use, intrinsic :: iso_c_binding
   type(c_ptr) ipdata
   interface
      subroutine rlsara_c (ipd) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr),value :: ipd
      end subroutine rlsara_c
   end interface
   call rlsara_c(ipdata)
end subroutine LCMDRD
!
subroutine LCMOP(iplist, name, imp, medium, impx)
   ! open a LCM object
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist
   character(len=*) :: name
   character(kind=c_char), dimension(73) :: name73
   integer imp, medium, impx
   interface
      subroutine lcmop_c (iplist, namp, imp, medium, impx) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: imp, medium, impx
      end subroutine lcmop_c
   end interface
   call STRCUT(name73, name)
   call lcmop_c(iplist, name73, imp, medium, impx)
end subroutine LCMOP
!
subroutine LCMPPD(iplist, name, ilong, itype, pt_data)
   ! store a record in an associative table via its c_ptr pointer
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   interface
      subroutine lcmppd_c (iplist, namp, ilong, itype, iofdum) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: iofdum
      end subroutine lcmppd_c
   end interface
   call STRCUT(name13, name)
   call lcmppd_c(iplist, name13, ilong, itype, pt_data)
   pt_data=c_null_ptr
end subroutine LCMPPD
!
subroutine LCMGPD(iplist, name, pt_data)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   interface
      subroutine lcmgpd_c (iplist, namp, iofdum) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr) :: iofdum
      end subroutine lcmgpd_c
   end interface
   call STRCUT(name13, name)
   call lcmgpd_c(iplist, name13, pt_data)
end subroutine LCMGPD
!
subroutine LCMLEN(iplist, name, ilong, itylcm)
   ! recover length and type of a record in an associative table
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist
   character(len=*),intent(in) :: name
   integer :: ilong, itylcm
   character(kind=c_char), dimension(13) :: name13
   interface
      subroutine lcmlen_c (iplist, namp, ilong, itylcm) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: ilong, itylcm
      end subroutine lcmlen_c
   end interface
   call STRCUT(name13, name)
   call lcmlen_c(iplist, name13, ilong, itylcm)
end subroutine LCMLEN
!
subroutine LCMINF(iplist, fnamlcm, fnammy, fempty, ilong, flcml)
   ! recover general info about a LCM object
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist
   character(len=*) :: fnamlcm, fnammy
   logical :: fempty, flcml
   integer :: empty, ilong, lcml, access
   character(kind=c_char), dimension(73) :: namlcm
   character(kind=c_char), dimension(13) :: nammy
   interface
      subroutine lcminf_c (iplist, namlcm, nammy, empty, ilong, lcml, access) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namlcm, nammy
         integer(c_int) :: empty, ilong, lcml, access
      end subroutine lcminf_c
   end interface
   call lcminf_c(iplist, namlcm, nammy, empty, ilong, lcml, access)
   call STRFIL(fnamlcm, namlcm)
   call STRFIL(fnammy, nammy)
   fempty=(empty == 1)
   flcml=(lcml == 1)
end subroutine LCMINF
!
subroutine LCMNXT(iplist, name)
   ! recover name of next record in an associative table
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist
   character(len=*) :: name
   character(kind=c_char), dimension(13) :: name13
   interface
      subroutine lcmnxt_c (iplist, namp) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
      end subroutine lcmnxt_c
   end interface
   call STRCUT(name13, name)
   call lcmnxt_c(iplist, name13)
   call STRFIL(name, name13)
end subroutine LCMNXT
!
subroutine LCMVAL(iplist, name)
   ! validate an associative table, starting from name
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   interface
      subroutine lcmval_c (iplist, namp) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
      end subroutine lcmval_c
   end interface
   call STRCUT(name13, name)
   call lcmval_c(iplist, name13)
end subroutine LCMVAL
!
subroutine LCMDEL(iplist, name)
   ! delete a record in an associative table
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   interface
      subroutine lcmdel_c (iplist, namp) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
      end subroutine lcmdel_c
   end interface
   call STRCUT(name13, name)
   call lcmdel_c(iplist, name13)
end subroutine LCMDEL
!
function LCMDID(iplist, name)
   ! create/access a daughter table in a parent table in modification mode
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: LCMDID,iplist
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   interface
      function lcmdid_c (iplist, namp) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: lcmdid_c,iplist
         character(kind=c_char), dimension(*) :: namp
      end function lcmdid_c
   end interface
   call STRCUT(name13, name)
   LCMDID = lcmdid_c(iplist, name13)
end function LCMDID
!
function LCMLID(iplist, name, ilong)
   ! create/access a daughter list in a parent table in modification mode
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: LCMLID,iplist
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong
   interface
      function lcmlid_c (iplist, namp, ilong) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: lcmlid_c,iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong
      end function lcmlid_c
   end interface
   call STRCUT(name13, name)
   LCMLID = lcmlid_c(iplist, name13, ilong)
end function LCMLID
!
function LCMDIL(iplist, ipos)
   ! create/access a daughter table in a parent list in modification mode
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: LCMDIL,iplist
   integer :: ipos
   interface
      function lcmdil_c (iplist, ipos) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: lcmdil_c,iplist
         integer(c_int), value :: ipos
      end function lcmdil_c
   end interface
   LCMDIL = lcmdil_c(iplist, ipos-1)
end function LCMDIL
!
function LCMLIL(iplist, ipos, ilong)
   ! create/access a daughter list in a parent list in modification mode
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: LCMLIL,iplist
   integer :: ipos, ilong
   interface
      function lcmlil_c (iplist, ipos, ilong) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: lcmlil_c,iplist
         integer(c_int), value :: ipos, ilong
      end function lcmlil_c
   end interface
   LCMLIL = lcmlil_c(iplist, ipos-1, ilong)
end function LCMLIL
!
function LCMGID(iplist, name)
   ! access a daughter table/list in a parent table in read-only mode
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: LCMGID,iplist
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   interface
      function lcmgid_c (iplist, namp) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: lcmgid_c,iplist
         character(kind=c_char), dimension(*) :: namp
      end function lcmgid_c
   end interface
   call STRCUT(name13, name)
   LCMGID = lcmgid_c(iplist, name13)
end function LCMGID
!
function LCMGIL(iplist, ipos)
   ! access a daughter table/list in a parent list in read-only mode
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: LCMGIL,iplist
   integer :: ipos
   interface
      function lcmgil_c (iplist, ipos) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: lcmgil_c,iplist
         integer(c_int), value :: ipos
      end function lcmgil_c
   end interface
   LCMGIL = lcmgil_c(iplist, ipos-1)
end function LCMGIL
!
subroutine LCMSIX(iplist, name, iact)
   ! create/access a daughter table in a parent table
   ! depreciated: better to use LCMDID
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist
   character(len=*),intent(in) :: name
   integer :: iact
   character(kind=c_char), dimension(13) :: name13
   interface
      subroutine lcmsix_c (iplist, namp, iact) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: iact
      end subroutine lcmsix_c
   end interface
   call STRCUT(name13, name)
   call lcmsix_c(iplist, name13, iact)
end subroutine LCMSIX
!
subroutine LCMPPL(iplist, ipos, ilong, itype, pt_data)
   ! store a record in an heterogeneous list via its c_ptr pointer
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   interface
      subroutine lcmppl_c (iplist, ipos, ilong, itype, iofdum) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: iofdum
      end subroutine lcmppl_c
   end interface
   call lcmppl_c(iplist, ipos-1, ilong, itype, pt_data)
   pt_data=c_null_ptr
end subroutine LCMPPL
!
subroutine LCMLEL(iplist, ipos, ilong, itylcm)
   ! recover length and type of a record in an heterogeneous list
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist
   integer :: ipos, ilong, itylcm
   interface
      subroutine lcmlel_c (iplist, ipos, ilong, itylcm) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         integer(c_int) :: ilong, itylcm
      end subroutine lcmlel_c
   end interface
   call lcmlel_c(iplist, ipos-1, ilong, itylcm)
end subroutine LCMLEL
!
subroutine LCMGPL(iplist, ipos, pt_data)
   ! recover a record from an heterogeneous list via its c_ptr pointer
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   interface
      subroutine lcmgpl_c (iplist, ipos, iofdum) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist, iofdum
         integer(c_int), value :: ipos
      end subroutine lcmgpl_c
   end interface
   call lcmgpl_c(iplist, ipos-1, pt_data)
end subroutine LCMGPL
!
subroutine LCMCL(iplist, iact)
   ! close a LCM object
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist
   integer :: iact
   interface
      subroutine lcmcl_c (iplist, iact) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: iact
      end subroutine lcmcl_c
   end interface
   call lcmcl_c(iplist, iact)
end subroutine LCMCL
!
subroutine LCMEQU(iplis1, iplis2)
   ! deep copy of a LCM object
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplis1, iplis2
   interface
      subroutine lcmequ_c (iplis1, iplis2) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplis1, iplis2
      end subroutine lcmequ_c
   end interface
   call lcmequ_c(iplis1, iplis2)
end subroutine LCMEQU
!
subroutine LCMEXP(iplist, impx, nunit, imode, idir)
   ! import/export a LCM object
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, file = c_null_ptr
   integer, intent(in) :: impx, nunit, imode, idir
   character(len=72) :: filename
   character(kind=c_char), dimension(73) :: filename_c
   integer(c_int) :: ier
   interface
      function fopen (filename_c, mode) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) fopen
         character(kind=c_char), dimension(*) :: filename_c, mode
      end function fopen
   end interface
   interface
      function fclose (file) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) fclose
         type(c_ptr), value :: file
      end function fclose
   end interface
   interface
      subroutine lcmexp_c (iplist, impx, file, imode, idir) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         type(c_ptr), value :: file
         integer(c_int), value :: impx, imode, idir
      end subroutine lcmexp_c
   end interface
   interface
      function stdfil_c (s) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) stdfil_c
         character(kind=c_char), dimension(*) :: s
      end function stdfil_c
   end interface
!
   if(nunit == 0) then
      file=c_null_ptr
   else if(nunit == 6) then
      file=stdfil_c("stdout"//c_null_char)
      flush(6)
   else
      inquire(nunit,name=filename)
      close(nunit,status='keep')
      call STRCUT(filename_c, filename)
      if(imode == 1) then
         if(idir == 1) then
            file=fopen(filename_c, "wb"//c_null_char)
         else
            file=fopen(filename_c, "rb"//c_null_char)
         endif
      else if(imode == 2) then
         if(idir == 1) then
            file=fopen(filename_c, "w"//c_null_char)
         else
            file=fopen(filename_c, "r"//c_null_char)
         endif
      endif
      if(.not.c_associated(file)) call XABORT('LCMEXP: UNABLE TO OPEN FILE '//filename(:44))
   endif
   call lcmexp_c(iplist, impx, file, imode, idir)
   if(nunit /= 6) then
      ier = fclose(file)
      if(ier /= 0) call XABORT('LCMEXP: UNABLE TO CLOSE FILE '//filename(:43))
      if(imode == 1) then
        open(nunit,file=filename,status='old',form='unformatted',position='append')
      else
        open(nunit,file=filename,status='old',position='append')
      endif
   endif
end subroutine LCMEXP
!
!-----------------------------------------------------------------------
! additionnal lcm subroutine specific to Version 3
! R. Chambon (based on Version 4)
!-----------------------------------------------------------------------
subroutine LCMEXPV3(iplist, impx, nunit, imode, idir)
   ! import/export a LCM object V3
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, file = c_null_ptr
   integer, intent(in) :: impx, nunit, imode, idir
   character(len=72) :: filename
   character(kind=c_char), dimension(73) :: filename_c
   integer(c_int) :: ier
   interface
      function fopen (filename_c, mode) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) fopen
         character(kind=c_char), dimension(*) :: filename_c, mode
      end function fopen
   end interface
   interface
      function fclose (file) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) fclose
         type(c_ptr), value :: file
      end function fclose
   end interface
   interface
      subroutine lcmexpv3_c (iplist, impx, file, imode, idir) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         type(c_ptr), value :: file
         integer(c_int), value :: impx, imode, idir
      end subroutine lcmexpv3_c
   end interface
   interface
      function stdfil_c (s) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) stdfil_c
         character(kind=c_char), dimension(*) :: s
      end function stdfil_c
   end interface
!
   if(nunit == 0) then
      file=c_null_ptr
   else if(nunit == 6) then
      file=stdfil_c("stdout"//c_null_char)
   else
      inquire(nunit,name=filename)
      close(nunit,status='keep')
      call STRCUT(filename_c, filename)
      if(imode == 1) then
         if(idir == 1) then
            file=fopen(filename_c, "wb"//c_null_char)
         else
            file=fopen(filename_c, "rb"//c_null_char)
         endif
      else if(imode == 2) then
         if(idir == 1) then
            file=fopen(filename_c, "w"//c_null_char)
         else
            file=fopen(filename_c, "r"//c_null_char)
         endif
      endif
      if(.not.c_associated(file)) call XABORT('LCMEXPV3: UNABLE TO OPEN FILE '//filename(:44))
   endif
   call lcmexpv3_c(iplist, impx, file, imode, idir)
   if(nunit /= 6) then
      ier = fclose(file)
      if(ier /= 0) call XABORT('LCMEXPV3: UNABLE TO CLOSE FILE '//filename(:43))
      if(imode == 1) then
        open(nunit,file=filename,status='old',form='unformatted',position='append')
      else
        open(nunit,file=filename,status='old',position='append')
      endif
   endif
end subroutine LCMEXPV3
end module LCMAUX
