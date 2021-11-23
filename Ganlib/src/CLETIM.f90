!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for CLETIM.
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
subroutine CLETIM(sec, nsec)
   ! abort execution
   use, intrinsic :: iso_c_binding
   integer :: sec, nsec
   interface
      subroutine cletim_c (sec, nsec) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) :: sec, nsec
      end subroutine cletim_c
   end interface
   call cletim_c(sec, nsec)
end subroutine CLETIM
