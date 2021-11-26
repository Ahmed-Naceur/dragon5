!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for fortran direct access in WIMS-AECL. Open,
! read or close indexed random file using fortran direct access files
! Subroutines:
!  OPNIND: open file and read master index
!  REDIND: read data on indexed file
!  CLSIND: close file
!
!Copyright:
! Copyright (C) 2020 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s): J. Donnelly
!
!Parameters: input/output
! IUNIT   read unit
! INDEX   index table (master index for OPNIND)
! LINDEX  length of index table (LINDEX >=(# of entries) + 1)
! DATA    data array to retreive from file
! NWORDS  lenght of data array to retreive from file
! KEY     location of data array in index
!
!Internal parameter description
! IOUT    output unit = 6
! NBLOCKS number of blocks per record = 256
! IOFSET  offset for index length = 65536
!
!-----------------------------------------------------------------------
!
module OPNMOD
   private
   integer iunitr
   public :: OPNIND, REDIND, CLSIND
   interface REDIND
      ! read data on indexed file
      module procedure REDIND_I1, REDIND_R1
   end interface
contains
   subroutine OPNIND(iunit,index,lindex)
   parameter (iout=6,nblock=256,iofset=65536)
   integer    iunit,index(lindex),lindex
   logical    exst, opnd
   character  dirst*7
   !----
   ! unit number must be > zero
   !----
   if(iunit.le.0) then
     write(iout,6000) iunit
     call XABORT('OPNIND: (readonly) illegal unit number')
   endif
   !----
   ! find out file status, and if unit already associated with
   ! an open file
   !----
   inquire(unit=iunit,exist=exst,opened=opnd,direct=dirst)
   if(.not.opnd) then
     !----
     ! file closed
     !----
     write(iout,6010) iunit
     call XABORT('OPNIND: (readonly) file not opened')
   endif
   if( exst .and. dirst .eq. 'no' ) then
     !----
     ! file already exists, but is not direct access
     !----
     write(iout,6020) iunit
     call XABORT('OPNIND: (readonly) file is nor direct access')
   endif
   if(.not.exst) then
     !----
     ! file does not exists
     !----
     write(iout,6030) iunit
     call XABORT('OPNIND: (readonly) file does not exists ')
   endif
   !----
   ! process the file master index
   !----
   iunitr=iunit
   irec = 1
   minw = 1
   40 continue
   maxw = min0( minw + nblock - 1 , lindex )
   ierr=0
   read(iunit,rec=irec,iostat=ierr) (index(i),i=minw,maxw)
   if(ierr.ne.0) then
      write(iout,6040) iunit,ierr
      call XABORT('OPNIND: read master record error')
   endif
   irec = irec + 1
   if( maxw .ne. lindex ) then
     minw = maxw + 1
     go to 40
   endif
   return
   !----
   !  format
   !----
   6000 format(' //// OPNIND: file, ',i5,' invalid')
   6010 format(' //// OPNIND: file, ',i5,' has not been opened with KDROPN')
   6020 format(' //// OPNIND: unit ',i5,' is not a direct access file')
   6030 format(' //// OPNIND: unit ',i5,' does not exists (readonly version)')
   6040 format(' //// OPNIND: error during reading of master index ', &
               'of unit ',i5,' status = ',i5)
   end subroutine OPNIND
   !----
   subroutine CLSIND(iunit)
   parameter (iout=6)
   integer    iunit
   if(iunitr.ne.iunit) then
     write(iout,6100) iunit
     call XABORT('CLSIND: file not opened by OPNIND')
   endif
   iunitr=kdrcls(iunit,1)
   if(iunitr.ne.0) then
     write(iout,6110) iunitr
     call XABORT('CLSIND: error in file closing')
   endif
   return
   !----
   ! format
   !----
   6100 format(' //// CLSIND: unit ',i5,' not found')
   6110 format(' //// CLSIND: error status =',i5,' from kdrcls')
   end subroutine CLSIND
   !----
   subroutine REDIND_I1(iunit,index,lindex,data,nwords,key)
   parameter (iout=6,nblock=256,iofset=65536)
   integer    iunit,index(lindex),lindex,nwords,key
   integer    data(nwords)
   !
   if(iunitr.ne.iunit) then
     write(iout,6200) iunit
     call XABORT('REDIND_I1: file not opened by OPNIND')
   endif
   !----
   ! validate key number
   !----
   if(key.le.0.or.key.ge.lindex) then
     write(iout,6210) iunit, key
     call XABORT('REDIND_I1: invalid key')
   endif
   !----
   ! key number valid, validate record number
   !----
   indr=index(key+1)
   if(indr.eq.0) then
     write(iout,6220) iunit, key
     call XABORT('REDIND_I1: invalid record number for key')
   endif
   !----
   ! validate record length
   !----
   lrecrd = (nwords-1)/nblock + 1
   loldrc = indr/iofset
   if(loldrc.lt.lrecrd) then
     write(iout,6230) iunit, key
     call XABORT('REDIND_I1: invalid record length')
   endif
   !----
   ! record found, read the data
   !----
   nrec = mod( indr, iofset )
   minw = 1
   50 continue
   maxw = min0( minw + nblock - 1 , nwords )
   ierr=0
   read(iunit,rec=nrec,iostat=ierr) (data(i),i=minw,maxw)
   if(ierr.ne.0) then
     write(iout,6240) iunit,ierr
     call XABORT('REDIND_I1: read record error')
   endif
   nrec = nrec + 1
   if( maxw .ne. nwords ) then
     minw = maxw + 1
     go to 50
   endif
   return
   !----
   ! format
   !----
   6200 format(' //// REDIND_I1: unit ',i5,' not found')
   6210 format(' //// REDIND_I1: invalid record number, unit ',i5,' key= ',i10)
   6220 format(' //// REDIND_I1: non-existant record, unit ',i5, &
               ' record key =',i10)
   6230 format(' //// REDIND_I1: data count exceeds record, unit ',i5, &
               ' record key =',i10)
   6240 format(' //// REDIND_I1: error during reading of record ', &
               'of unit ',i5,' status = ',i5)
   end subroutine REDIND_I1
   !----
   subroutine REDIND_R1(iunit,index,lindex,data,nwords,key)
   parameter (iout=6,nblock=256,iofset=65536)
   integer    iunit,index(lindex),lindex,nwords,key
   real       data(nwords)
   !
   if(iunitr.ne.iunit) then
     write(iout,6200) iunit
     call XABORT('REDIND_R1: file not opened by OPNIND')
   endif
   !----
   ! validate key number
   !----
   if(key.le.0.or.key.ge.lindex) then
     write(iout,6210) iunit, key
     call XABORT('REDIND_R1: invalid key')
   endif
   !----
   ! key number valid, validate record number
   !----
   indr=index(key+1)
   if(indr.eq.0) then
     write(iout,6220) iunit, key
     call XABORT('REDIND_R1: invalid record number for key')
   endif
   !----
   ! validate record length
   !----
   lrecrd = (nwords-1)/nblock + 1
   loldrc = indr/iofset
   if(loldrc.lt.lrecrd) then
     write(iout,6230) iunit, key
     call XABORT('REDIND_R1: invalid record length')
   endif
   !----
   ! record found, read the data
   !----
   nrec = mod( indr, iofset )
   minw = 1
   50 continue
   maxw = min0( minw + nblock - 1 , nwords )
   ierr=0
   read(iunit,rec=nrec,iostat=ierr) (data(i),i=minw,maxw)
   if(ierr.ne.0) then
     write(iout,6240) iunit,ierr
     call XABORT('REDIND_R1: read record error')
   endif
   nrec = nrec + 1
   if( maxw .ne. nwords ) then
     minw = maxw + 1
     go to 50
   endif
   return
   !----
   ! format
   !----
   6200 format(' //// REDIND_R1: unit ',i5,' not found')
   6210 format(' //// REDIND_R1: invalid record number, unit ',i5,' key= ',i10)
   6220 format(' //// REDIND_R1: non-existant record, unit ',i5, &
               ' record key =',i10)
   6230 format(' //// REDIND_R1: data count exceeds record, unit ',i5, &
               ' record key =',i10)
   6240 format(' //// REDIND_R1: error during reading of record ', &
               'of unit ',i5,' status = ',i5)
   end subroutine REDIND_R1
end module OPNMOD
