subroutine KDRCPU(cpusec)
!
!-----------------------------------------------------------------------
!
!Purpose:
! system clock support.
!
!Copyright:
! Copyright (C) 2002 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s): A. Hebert
!
!Parameters: output
!  cpusec : number of seconds elapsed since the first call to KDRCPU.
!
!-----------------------------------------------------------------------
!
   integer :: itloc,irate
   integer,save :: isave=0,itloc0
!
   if(isave==0) then
      call system_clock(count=itloc0)
      isave=1
   endif
   call system_clock(count=itloc,count_rate=irate)
   cpusec=real(itloc-itloc0)/real(irate)
   return
end subroutine KDRCPU
