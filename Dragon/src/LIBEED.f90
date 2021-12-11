!
!-----------------------------------------------------------------------
!
!Purpose:
! Read bcd-formatted MATXS format records.
! LIBEED:  transfer data from CCCC file to array.
! LIBCLS:  close file and release unit number.
!
!Copyright:
! Copyright (C) 2021 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s): A. Hebert
!
!Parameters: input
! iucccc  file unit
! numrec  record number to read
! nwds    number of words to read
!
!Parameters: output
! ra      location in central memory where information is to be stored
!
!-----------------------------------------------------------------------
!
module LIBEEDR
   use, intrinsic :: iso_c_binding
   private
   integer, parameter :: nbunit=99,ntoc=5000
   integer, save :: ipunit,nbrec,npart,ntype,nmat,nsub
   integer, save :: iucccc_old=0
   integer, save, dimension(ntoc) :: atoc(1:ntoc)=0
   double precision, allocatable, save, dimension(:) :: hmatn
   integer, allocatable, save, dimension(:) :: nsubm
   character(len=3), save, dimension(ntoc) :: btoc
   public :: LIBEED, LIBCLS
contains
   subroutine LIBEED(iucccc,numrec,ra,nwds)
      integer, intent(in) :: iucccc,numrec,nwds
      real, intent(out), target :: ra(nwds)
      character(len=72) :: text72
      character(len=131) :: hsmg
      !
      type(c_ptr) ra_ptr
      integer, pointer :: ia(:)
      double precision, pointer :: da(:)
      !
      if (numrec.eq.0) then
        call XABORT('LIBEED: record number 0 cannot be read '// &
                    'from cccc file')
      endif
      !
      if (iucccc.ne.iucccc_old) iucccc_old=0
      if (iucccc_old.eq.0) then
        iucccc_old=iucccc
        ipunit=0
        npart=0
        ntype=0
        nmat=0
        rewind iucccc
        i=0
        nbrec=0
        do
          i=i+1
          read(iucccc,'(A72)',end=10) text72
          if((text72(3:3).eq.'d').or.(text72(3:3).eq.'v')) then
            nbrec=nbrec+1
            if(nbrec.gt.ntoc) call XABORT('LIBEED: ntoc overflow.')
            atoc(nbrec)=i
            btoc(nbrec)=text72(1:3)
          endif
        enddo
        10 rewind iucccc
      endif
      !
      if (nwds.eq.0) return
      if (numrec.eq.1) then
        rewind iucccc
      else if(numrec.gt.nbrec) then
        call XABORT('LIBEED: nbrec overflow.')
      else
        nskip=atoc(numrec)-ipunit-1
        if (nskip.gt.0) then
          do i=1,nskip
            read(iucccc,'(a72)') text72
          enddo
        else if (nskip.lt.0) then
          do i=1,-nskip
            backspace iucccc
          enddo
        endif
      endif
      ra_ptr=c_loc(ra)
      call c_f_pointer(ra_ptr,ia,(/ nwds /))
      call c_f_pointer(ra_ptr,da,(/ nwds /))
      ipunit=atoc(numrec)
      if(btoc(numrec).eq.' 0v') then
        read(iucccc,'(4x,a8,1x,2a8,1x,i6)') (da(jj),jj=1,3),ia(7)                                       
      else if(btoc(numrec).eq.' 1d') then
        read(iucccc,'(6x,6i6)') (ia(jj),jj=1,nwds)
        npart=ia(1)
        ntype=ia(2)
        nmat=ia(4)
        allocate(hmatn(nmat),nsubm(nmat))
      else if(btoc(numrec).eq.' 2d') then
        read(iucccc,'(4x/(9a8))') (da(jj),jj=1,nwds/2)
        ipunit=ipunit+1+(nwds/2-1)/9
      else if(btoc(numrec).eq.' 3d') then
        ndr=npart+ntype+nmat
        nir=npart+2*ntype+2*nmat
        if(2*ndr+nir.ne.nwds) call XABORT('LIBEED: invalid nwds(1).')
        read(iucccc,'(8x,8a8:/(9a8))') (da(jj),jj=1,ndr)
        ipunit=ipunit+ndr/9
        read(iucccc,'(12i6)') (ia(2*ndr+i),i=1,nir)
        ipunit=ipunit+1+(nir-1)/12
        if(.not.allocated(hmatn)) call XABORT('LIBEED: hmatn not allocated.')
        hmatn(:nmat)=da(npart+ntype+1:ndr)
        nsubm(:nmat)=ia(2*ndr+npart+2*ntype+1:2*ndr+npart+2*ntype+nmat)
      else if(btoc(numrec).eq.' 6d') then
        ndr=nwds/4
        read(iucccc,'(8x,8a8:/(9a8))') (da(jj),jj=1,ndr)
        ipunit=ipunit+ndr/9
        read(iucccc,'(12i6)') (ia(2*ndr+i),i=1,2*ndr)
        ipunit=ipunit+1+(2*ndr-1)/12
      else if((btoc(numrec).eq.' 4d').or.(btoc(numrec).eq.' 7d').or. &
              (btoc(numrec).eq.' 9d').or.(btoc(numrec).eq.'10d')) then
        read(iucccc,'(12x,5e12.0:/(6e12.0))') (ra(jj),jj=1,nwds)
        ipunit=ipunit+nwds/6
      else if(btoc(numrec).eq.' 5d') then
        read(iucccc,'(4x,a8,e12.0)') da(1),ra(3)
        nsub=0
        do i=1,nmat
          if(hmatn(i).eq.da(1)) then
            nsub=nsubm(i)
            go to 20
          endif
        enddo
        write(hsmg,'(49HLIBEED: unable to find material control data for , &
        & a,1h.)') da(1)
        call XABORT(HSMG)
        20 do i=1,nsub
          ll=4+6*(i-1)
          read(iucccc,'(2e12.0,4i6)') ra(ll),ra(ll+1),ia(ll+2),ia(ll+3), &
          ia(ll+4),ia(ll+5)
          ipunit=ipunit+1
        enddo
        if(3+6*nsub.ne.nwds) call XABORT('LIBEED: invalid nwds(2).')
      else if(btoc(numrec).eq.' 8d') then
        read(iucccc,'(8x,a8/(12i6))') da(1),(ia(jj),jj=3,nwds)
        ipunit=ipunit+1+(nwds-3)/12
      else
        call XABORT('LIBEED: invalid record type.')
      endif
   end subroutine LIBEED
   !
   subroutine LIBCLS()
      if(allocated(hmatn)) deallocate(nsubm,hmatn)
      ipunit=0
      npart=0
      ntype=0
      nmat=0
      iucccc_old=0
   end subroutine LIBCLS
end module LIBEEDR
