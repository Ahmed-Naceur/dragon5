!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for kdi.
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
function KDIOP(name, iactio)
   ! open a KDI file
   use, intrinsic :: iso_c_binding
   use LCMAUX
   type(c_ptr) KDIOP
   character(len=*) :: name
   integer :: iactio
   character(kind=c_char), dimension(13) :: name13
   interface
      function kdiop_c (name_c, iactio) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) kdiop_c
         character(kind=c_char), dimension(*) :: name_c
         integer(c_int), value :: iactio
      end function kdiop_c
   end interface
   call STRCUT(name13, name)
   KDIOP=kdiop_c(name13, iactio)
end function KDIOP
!
subroutine KDIPUT(my_file, idata, iofset, length)
   ! store a data array in a KDI file at offset iofset
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: my_file, pt_data
   integer, target, dimension(*) :: idata
   integer :: iofset, length
   interface
      subroutine kdiput_c (my_file, idata, iofset, length) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: my_file, idata
         integer(c_int), value :: iofset, length
      end subroutine kdiput_c
   end interface
   pt_data=c_loc(idata)
   call kdiput_c(my_file, pt_data, iofset, length)
end subroutine KDIPUT
!
subroutine KDIGET(my_file, idata, iofset, length)
   ! read a data array from a KDI file at offset iofset
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: my_file, pt_data
   integer, target, dimension(*) :: idata
   integer :: iofset, length
   interface
      subroutine kdiget_c (my_file, idata, iofset, length) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: my_file, idata
         integer(c_int), value :: iofset, length
      end subroutine kdiget_c
   end interface
   pt_data=c_loc(idata)
   call kdiget_c(my_file, pt_data, iofset, length)
end subroutine KDIGET
!
function KDICL(my_file, istatu)
   ! close a KDI file
   use, intrinsic :: iso_c_binding
   integer(c_int) KDICL
   type(c_ptr) :: my_file
   integer :: istatu
   interface
      function kdicl_c (my_file, istatu) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) kdicl_c
         type(c_ptr), value :: my_file
         integer(c_int), value :: istatu
      end function kdicl_c
   end interface
   KDICL=kdicl_c(my_file, istatu)
end function KDICL
