integer function GANDRV(hmodul,nentry,hentry,ientry,jentry,kentry)
!
!-----------------------------------------------------------------------
!
!Purpose:
! standard utility operator driver for Ganlib.
!
!Copyright:
! Copyright (C) 2002 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version
!
!Author(s): A. Hebert
!
!Parameters: input/output
! hmodul  name of the operator.
! nentry  number of LCM objects or files used by the operator.
! hentry  name of each LCM object or file.
! ientry  type of each LCM object or file:
!         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
!         =4 sequential ascii file.
! jentry  access of each LCM object or file:
!         =0 the LCM object or file is created;
!         =1 the LCM object or file is open for modifications;
!         =2 the LCM object or file is open in read-only mode.
! kentry  LCM object address or file unit number.
!
!Parameters: output
! kdrstd  completion flag (=0: operator hmodul exists; =1: does not exists).
!
!-----------------------------------------------------------------------
!
   use, intrinsic :: iso_c_binding
!----
!  subroutine arguments
!----
   integer nentry
   character hmodul*(*),hentry(nentry)*12
   integer ientry(nentry),jentry(nentry)
   type(c_ptr) kentry(nentry)
!
   real tbeg,tend
   double precision dmemb,dmemd
!
   GANDRV=0
   call KDRCPU(tbeg)
   call KDRMEM(dmemb)
   if(hmodul == 'EQU:' )then
!     standard equality module.
      call DRVEQU(nentry,hentry,ientry,jentry,kentry)
   else if(hmodul == 'GREP:') then
!     standard grep module.
      call DRVGRP(nentry,hentry,ientry,jentry,kentry)
   else if(hmodul == 'UTL:') then
!     standard LCM/XSM utility module.
      call DRVUTL(nentry,hentry,ientry,jentry,kentry)
   else if(hmodul == 'ADD:') then
!     standard addition module.
      call DRVADD(nentry,hentry,ientry,jentry,kentry)
   else if(hmodul == 'MPX:') then
!     standard multiplication module.
      call DRVMPX(nentry,hentry,ientry,jentry,kentry)
   else if(hmodul == 'STAT:') then
!     standard compare module.
      call DRVSTA(nentry,hentry,ientry,jentry,kentry)
   else if(hmodul == 'BACKUP:') then
!     standard backup module.
      call DRVBAC(nentry,hentry,ientry,jentry,kentry)
   else if(hmodul == 'RECOVER:') then
!     standard recovery module.
      call DRVREC(nentry,hentry,ientry,jentry,kentry)
   else if(hmodul == 'FIND0:') then
!     standard module to find zero of a continuous function.
      call DRV000(nentry,hentry,ientry,jentry,kentry)
   else if(hmodul == 'MSTR:') then
!     manage user-defined structures.
      call MSTR(nentry,hentry,ientry,jentry,kentry)
   else if(hmodul == 'MODUL1:') then
!     user-defined module.
      call DRVMO1(nentry,hentry,ientry,jentry,kentry)
   else if(hmodul == 'ABORT:') then
!     requested abort.
      call XABORT('GANDRV: requested abort.')
#if defined(MPI)
   elseif(hmodul == 'DRVMPI:') then
!     initialize mpi.
      call DRVMPI(nentry,hentry,ientry,jentry,kentry)
   elseif(hmodul == 'SNDMPI:') then
!     export LCM or XSM using mpi.
      call SNDMPI(nentry,hentry,ientry,jentry,kentry)
#endif /* defined(MPI) */
   else
      GANDRV=1
   endif
   call KDRCPU(tend)
   call KDRMEM(dmemd)
   write(6,5000) hmodul,(tend-tbeg),real(dmemd-dmemb)
   return
!
   5000 format('-->>module ',a12,': time spent=',f13.3,' memory usage=',1p,e10.3)
end function GANDRV
