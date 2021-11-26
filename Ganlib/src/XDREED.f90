!
!-----------------------------------------------------------------------
!
!Purpose:
! Read or write CCCC format records.
! XDREED:  transfer data from CCCC file to array.
! XDRITE:  transfer data from array to CCCC file.
! XDRCLS:  close file and release unit number.
!
!Copyright:
! Copyright (C) 1991 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s): G. Marleau
!
!Parameters: input
! iucccc  file unit
! numrec  record number to read
! nwds    number of words to read
!
!Parameters: input/output
! array   location in central memory where information is to be stored
!
!Reference:
! R. D. O'Dell, 'Standard interface files and procedures for
! reactor physics codes, Version IV', Los Alamos National
! Laboratory Report LA-6941-MS (Sept. 1977).
!
!-----------------------------------------------------------------------
!
module XDRMOD
   private
   integer, parameter :: nbunit=99
   integer, save, dimension(nbunit) :: ipunit(1:nbunit)=0
   public :: XDREED, XDRITE, XDRCLS
contains
   subroutine XDREED(iucccc,numrec,array,nwds)
      integer, intent(in) :: iucccc,numrec,nwds
      real, intent(out) :: array(nwds)
      !
      if (numrec.eq.0) then
        call XABORT('XDREED: record number 0 cannot be read '// &
                    'from cccc file')
      endif
      if (nwds.eq.0) return
      if (numrec.eq.1) then
        rewind iucccc
      else
        nskip=numrec-ipunit(iucccc)-1
        if (nskip.gt.0) then
          do i=1,nskip
            read(iucccc) dum
          enddo
        else if (nskip.lt.0) then
          do i=1,-nskip
            backspace iucccc
          enddo
        endif
      endif
      read(iucccc) (array(jj),jj=1,nwds)
      ipunit(iucccc)=numrec
   end subroutine XDREED
   !
   subroutine XDRITE(iucccc,numrec,array,nwds)
      integer, intent(in) :: iucccc,numrec,nwds
      real, intent(in) :: array(nwds)
      !
      if (numrec.eq.0) then
        call XABORT('XDRITE: record number 0 cannot be written '// &
                    'on cccc file')
      endif
      if (nwds.eq.0) return
      if (numrec.eq.1) rewind iucccc
      write(iucccc) (array(jj),jj=1,nwds)
   end subroutine XDRITE
   !
   subroutine XDRCLS(iucccc)
      integer, intent(in) :: iucccc
      !
      ipunit(iucccc)=0
   end subroutine XDRCLS
end module XDRMOD
