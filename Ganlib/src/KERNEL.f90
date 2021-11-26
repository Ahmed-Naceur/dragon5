!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for CLE-2000. Call the CLE-2000 driver.
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
integer function KERNEL(dummod, iprint)
   use GANLIB
   implicit none
!----
!  subroutine arguments
!----
   integer :: iprint
   interface
      function dummod(cmodul, nentry, hentry, ientry, jentry, kentry, &
                     hparam_c) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) dummod
         character(kind=c_char), dimension(*) :: cmodul
         integer(c_int), value :: nentry 
         character(kind=c_char), dimension(13,*) :: hentry
         integer(c_int), dimension(nentry) :: ientry, jentry
         type(c_ptr), dimension(nentry) :: kentry
         character(kind=c_char), dimension(73,*) :: hparam_c
       end function dummod
   end interface
!----
!  local variables
!----
   interface
      function cle2000_c(ilevel, dummod_pt, filein, iprint, my_param) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) cle2000_c
         integer(c_int), value :: ilevel, iprint
         type(c_funptr), value :: dummod_pt
         character(kind=c_char), dimension(*) :: filein
         type(c_ptr), value :: my_param
      end function cle2000_c
   end interface
   interface
      function stdfil_c (s) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) stdfil_c
         character(kind=c_char), dimension(*) :: s
      end function stdfil_c
   end interface
   integer :: ilevel = 1
   type(c_funptr) :: dummod_pt
!----
!  call the driver
!----
   dummod_pt=c_funloc(dummod)
   KERNEL=cle2000_c(ilevel, dummod_pt, " "//c_null_char, iprint, c_null_ptr)
   return
end function KERNEL
