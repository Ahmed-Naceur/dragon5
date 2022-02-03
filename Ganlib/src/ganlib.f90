!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for Ganlib support. This module defines the
! interface prototypes of the Ganlib Fortran API, defines TYPE(C_PTR)
! and defines the external functions in the Ganlib API.
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
module GANLIB
   use FILMOD
   use LCMAUX
   use LCMMOD
   use OPNMOD
   use XDRMOD
   use, intrinsic :: iso_c_binding
   integer, parameter :: dp = kind(0.0d0)
   interface
      subroutine CUT(name1, name2, ilong)
         use, intrinsic :: iso_c_binding
         character(kind=c_char), dimension(*) :: name1
         character(len=*) :: name2
         integer :: ilong
      end subroutine CUT
   end interface
   interface
      subroutine FIL(name1, name2, ilong)
         use, intrinsic :: iso_c_binding
         character(len=*) :: name1
         character(kind=c_char), dimension(*) :: name2
         integer :: ilong
      end subroutine FIL
   end interface
   interface
      subroutine XABORT(msg)
         character(len=*) :: msg
      end subroutine XABORT
   end interface
   interface
      subroutine REDGET(ityp, nitma, flott, text, dflot)
         integer :: ityp, nitma
         real :: flott
         character(len=*) :: text
         double precision :: dflot
      end subroutine REDGET
   end interface
   interface
      subroutine REDPUT(ityp, nitma, flott, text, dflot)
         integer :: ityp, nitma
         real :: flott
         character(len=*) :: text
         double precision :: dflot
      end subroutine REDPUT
   end interface
   interface
      subroutine REDOPN(iinp1, iout1, nrec)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iinp1, iout1
         integer :: nrec
      end subroutine REDOPN
   end interface
   interface
      subroutine REDCLS(iinp1, iout1, nrec)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iinp1, iout1
         integer :: nrec
      end subroutine REDCLS
   end interface
   interface
      function KDIOP(name, iactio)
         use, intrinsic :: iso_c_binding
         type(c_ptr) KDIOP
         character(len=*) :: name
         integer :: iactio
      end function KDIOP
   end interface
   interface
      function KDICL(my_file, istatu)
         use, intrinsic :: iso_c_binding
         integer(c_int) KDICL
         type(c_ptr) :: my_file
         integer :: istatu
      end function KDICL
   end interface
   interface
      integer function GANDRV(hmodul, nentry, hentry, ientry, jentry, kentry)
         use, intrinsic :: iso_c_binding
         character(len=*), intent(in) :: hmodul
         integer, intent(in) :: nentry
         character(len=12), dimension(nentry), intent(in) :: hentry
         integer, dimension(nentry), intent(in) :: ientry, jentry
         type(c_ptr), dimension(nentry), intent(in) :: kentry
      end function GANDRV
   end interface
end module GANLIB
