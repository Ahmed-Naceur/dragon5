!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for XABORT.
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
subroutine XABORT(msg)
   ! abort execution
   use, intrinsic :: iso_c_binding
   use LCMAUX
   character(len=*) :: msg
   character(kind=c_char), dimension(73) :: msg_c
   interface
      subroutine xabort_c (msg) bind(c)
         use, intrinsic :: iso_c_binding
         character(kind=c_char), dimension(*) :: msg
      end subroutine xabort_c
   end interface
   flush(6)
   call STRCUT(msg_c, msg)
   call xabort_c(msg_c)
end subroutine XABORT
