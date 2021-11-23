!
!-----------------------------------------------------------------------
!
!Purpose:
! recover user memory used
!
!Copyright:
! Copyright (C) 2019 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s): A. Hebert
!
!Parameters: output
! utime   allocated memory in bytes.
!
!-----------------------------------------------------------------------
!
subroutine KDRMEM(utime)
   use, intrinsic :: iso_c_binding
   double precision :: utime
   interface
      function getusage () bind(c)
         use, intrinsic :: iso_c_binding
         real(c_double) :: getusage
      end function getusage
   end interface
   utime=getusage()
end subroutine KDRMEM
