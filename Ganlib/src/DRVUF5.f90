subroutine DRVUH5(nentry,hentry,ientry,jentry,kentry)
  !
  !-----------------------------------------------------------------------
  !
  !Purpose:
  ! standard utility module for HDF5 files.
  !
  !Copyright:
  ! Copyright (C) 2021 Ecole Polytechnique de Montreal
  ! This library is free software; you can redistribute it and/or
  ! modify it under the terms of the GNU Lesser General Public
  ! License as published by the Free Software Foundation; either
  ! version 2.1 of the License, or (at your option) any later version
  !
  !Author(s): A. Hebert
  !
  !Parameters: input/output
  ! NENTRY  number of LCM objects or files used by the operator.
  ! HENTRY  name of each LCM object or file:
  !         HENTRY(1): read-only or modification type(HDF5_FILE).
  ! IENTRY  type of each LCM object or file:
  !         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
  !         =4 sequential ascii file; =6 HDF5 file.
  ! JENTRY  access of each LCM object or file:
  !         =0 the LCM object or file is created;
  !         =1 the LCM object or file is open for modifications;
  !         =2 the LCM object or file is open in read-only mode.
  ! KENTRY  LCM object address or file unit number.
  !
  ! List of utility actions:
  !   DIR  : print the table of content.
  !   INFO : print information about a dataset.
  !   TEST : test if a group exists.
  !   IMPR : print a dataset.
  !   CREA : create a dataset.
  !   GREP : recover a single component in a dataset of rank 1.
  !
  !-----------------------------------------------------------------------
  !
  use hdf5_wrap
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env
  !----
  !  subroutine arguments
  !----
  integer :: nentry,ientry(nentry),jentry(nentry)
  type(c_ptr) :: kentry(nentry)
  character(len=12) :: hentry(nentry)
  !----
  !  local variables
  !----
  character :: text4*4,text12*12,text72*72,hsmg*131
  integer :: dimsr(5)
  integer :: rank,type,nbyte
  double precision :: dflott
  type(c_ptr) :: my_hdf5
  !----
  !  allocatable arrays
  !----
  integer, allocatable, dimension(:) :: nitmaV1
  integer, allocatable, dimension(:,:) :: nitmaV2
  character(len=32) :: text32V0
  character(len=64) :: text64V0
  character(len=32), allocatable, dimension(:) :: text32V1
  character(len=64), allocatable, dimension(:) :: text64V1
  real, allocatable, dimension(:) :: flottV1
  real, allocatable, dimension(:,:) :: flottV2
  double precision, allocatable, dimension(:) :: dflottV1
  double precision, allocatable, dimension(:,:) :: dflottV2
  !----
  !  parameter validation.
  !----
  if(nentry.eq.0) call XABORT('DRVUH5: parameter expected.')
  text12=hentry(1)
  ind=jentry(1)
  if(ientry(1).ne.6) call XABORT('DRVUH5: the utility module works on' &
  & //'ly for hdf5 files.')
  my_hdf5=kentry(1)
  !----
  !  perform some utility actions.
  !----
  10 call REDGET(indic,nitma,flott,text4,dflott)
  if(indic.ne.3) call XABORT('DRVUH5: character data expected.')
  20 if(text4.eq.'DIR') then
    ! print the group content
    flush(OUTPUT_UNIT)
    call REDGET(indprt,nitma,flott,text72,dflott)
    if(indprt.ne.3) call XABORT('DRVUH5: dataset name expected.')
    if((text72.eq.'INFO').or.(text72.eq.'TEST').or.(text72.eq.'IMPR').or. &
       (text72.eq.'CREA').or.(text72.eq.'DIR').or.(text72.eq.';')) then
      text4=text72(:4)
      call hdf5_list(my_hdf5,' ')
      go to 20
    endif
    call hdf5_list(my_hdf5,text72)
  else if(text4.eq.'INFO') then
    ! print a dataset.
    call REDGET(indprt,nitma,flott,text72,dflott)
    if(indprt.ne.3) call XABORT('DRVUH5: dataset name expected.')
    call hdf5_info(my_hdf5,text72,rank,type,nbyte,dimsr)
    write(6,'(/32h DRVUF5: information on dataset ,a,1h:)') text72
    write(6,*) 'rank=',rank,' type=',type,' nbyte=',nbyte,' dimsr=',dimsr(:rank)
  else if(text4.eq.'TEST') then
    ! test if a group exists.
    call REDGET(indprt,nitma,flott,text72,dflott)
    if(indprt.ne.3) call XABORT('DRVUH5: dataset name expected.')
    if (hdf5_group_exists(my_hdf5,text72)) then
      write(6,'(/15h DRVUF5: group ,a,8h exists.)') trim(text72)
    else
      write(6,'(/15h DRVUF5: group ,a,15h doesn''t exist.)') trim(text72)
    endif
  else if(text4.eq.'IMPR') then
    ! print a dataset.
    call REDGET(indprt,nitma,flott,text72,dflott)
    if(indprt.ne.3) call XABORT('DRVUH5: dataset name expected.')
    call hdf5_info(my_hdf5,text72,rank,type,nbyte,dimsr)
    write(6,'(/29h DRVUF5: printout of dataset ,a,1h:)') text72
    if(type.eq.1) then
      if((rank.eq.1).and.(dimsr(1).eq.1)) then
        call hdf5_read_data(my_hdf5,text72,nitma)
        write(6,'(4x,i12)') nitma
      else if(rank.eq.1) then
        call hdf5_read_data(my_hdf5,text72,nitmaV1)
        write(6,'(4x,10i12)') nitmaV1(:)
        deallocate(nitmaV1)
      else if(rank.eq.2) then
        call hdf5_read_data(my_hdf5,text72,nitmaV2)
        write(6,'(4x,10i12)') nitmaV2(:,:)
        deallocate(nitmaV2)
      else
        write(hsmg,100) type,rank
        call XABORT(hsmg)
      endif
    else if(type.eq.2) then
      if((rank.eq.1).and.(dimsr(1).eq.1)) then
        call hdf5_read_data(my_hdf5,text72,flott)
        write (6,'(1x,1p,e13.4)') flott
      else if(rank.eq.1) then
        call hdf5_read_data(my_hdf5,text72,flottV1)
        write (6,'(1x,1p,10e13.4)') flottV1(:)
        deallocate(flottV1)
      else if(rank.eq.2) then
        call hdf5_read_data(my_hdf5,text72,flottV2)
        write (6,'(1x,1p,10e13.4)') flottV2(:,:)
        deallocate(flottV2)
      else
        write(hsmg,100) type,rank
        call XABORT(hsmg)
      endif
    else if(type.eq.3) then
      if((rank.eq.1).and.(dimsr(1).eq.1)) then
        if(nbyte.le.32) then
          call hdf5_read_data(my_hdf5,text72,text32V0)
          write(6,'(4x,a)') text32V0
        else
          call hdf5_read_data(my_hdf5,text72,text64V0)
          write(6,'(4x,a)') text64V0
        endif
      else if(rank.eq.1) then
        if(nbyte.le.32) then
          call hdf5_read_data(my_hdf5,text72,text32V1)
          write(6,'(4x,5a32)') text32V1(:)
          deallocate(text32V1)
        else
          call hdf5_read_data(my_hdf5,text72,text64V1)
          write(6,'(4x,3a64)') text64V1(:)
          deallocate(text64V1)
        endif
      else
        write(hsmg,100) type,rank
        call XABORT(hsmg)
      endif
    else if(type.eq.4) then
      if((rank.eq.1).and.(dimsr(1).eq.1)) then
        call hdf5_read_data(my_hdf5,text72,dflott)
        write (6,'(1x,1p,d21.12)') dflott
      else if(rank.eq.1) then
        call hdf5_read_data(my_hdf5,text72,dflottV1)
        write (6,'(1x,1p,6d21.12)') dflottV1(:)
        deallocate(dflottV1)
      else if(rank.eq.2) then
        call hdf5_read_data(my_hdf5,text72,dflottV2)
        write (6,'(1x,1p,6d21.12)') dflottV2(:,:)
        deallocate(dflottV2)
      else
        write(hsmg,100) type,rank
        call XABORT(hsmg)
      endif
    else
      write(hsmg,100) type,rank
      call XABORT(hsmg)
    endif
  else if(text4.eq.'CREA') then
    if(ind.eq.2) call XABORT('DRVUF5: crea is a forbidden operation in read-only mode.')
    call REDGET(ntype,iset,flott,text72,dflott)
    indico=0
    ilong0=0
    if(ntype.eq.3) then
      call hdf5_info(my_hdf5,text72,rank,indico,nbyte,dimsr)
      if(rank.gt.1) call XABORT('DRVUF5: rank>1 forbidden.')
      ilong0=nbyte/4
    else
      call XABORT('DRVUF5: character data expected.')
    endif
    call REDGET(indic,nitma,flott,text4,dflott)
    if(indic.eq.1) then
      ilong2=nitma
    else if((indic.eq.3).and.(text4.eq.'=')) then
      call REDGET(indic,nitma,float,text4,dflott)
      if(indic.eq.1) then
        call hdf5_write_data(my_hdf5, text72, nitma)
      else if(indic.eq.2) then
        call hdf5_write_data(my_hdf5, text72, float)
      else if(indic.eq.3) then
        call hdf5_write_data(my_hdf5, text72, text4)
      else if(indic.eq.4) then
        call hdf5_write_data(my_hdf5, text72, dflott)
      else
        call XABORT('DRVUF5: invalid type.')
      endif
      go to 10
    else
      call XABORT('DRVUF5: integer data or = expected.')
    endif
    ilong1=1
    30 call REDGET(indic,ilong,flott,text4,dflott)
    if(indic.eq.1) then
      if(ilong0.eq.0) call XABORT('DRVUF5: lower index not expected.')
      ilong1=ilong2
      ilong2=ilong
      go to 30
    else if((indic.ne.3).or.(text4.ne.'=')) then
       call XABORT('DRVUF5: = sign expected.')
    endif
    call REDGET(indic,nitma,float,text4,dflott)
    if(indico.eq.99) then
       indico=indic
    else if(indic.ne.indico) then
       call XABORT('DRVUF5: inconsistent data type(1).')
    endif
    if(indic.eq.1) then
      if(ilong0.ne.0) then
        call hdf5_read_data(my_hdf5, text72, nitmaV1)
      else
        allocate(nitmaV1(ilong2))
      endif
      nitmaV1(ilong1)=nitma
    else if(indic.eq.2) then
      if(ilong0.ne.0) then
        call hdf5_read_data(my_hdf5, text72, flottV1)
      else
        allocate(flottV1(ilong2))
      endif
      flottV1(ilong1)=float
    else if(indic.eq.3) then
      if(ilong0.ne.0) then
        call hdf5_read_data(my_hdf5, text72, text32V1)
      else
        allocate(text32V1(ilong2))
      endif
      text32V1(ilong1)=text4
    else if(indic.eq.4) then
      if(ilong0.ne.0) then
        call hdf5_read_data(my_hdf5, text72, dflottV1)
      else
        allocate(dflottV1(ilong2))
      endif
      dflottV1(ilong1)=dflott
    endif
    do i=ilong1+1,ilong2
      call REDGET(indic,nitma,float,text4,dflott)
      if(indic.ne.indico) then
        call XABORT('DRVUF5: inconsistent data type(2).')
      else if(indic.eq.1) then
        nitmaV1(i)=nitma
      else if(indic.eq.2) then
        flottV1(i)=float
      else if(indic.eq.3) then
        text32V1(i)=text4
      else if(indic.eq.4) then
        dflottV1(i)=dflott
      endif
    enddo
    if(indico.eq.1) then
      call hdf5_write_data(my_hdf5, text72, nitmaV1)
      deallocate(nitmaV1)
    else if(indico.eq.2) then
      call hdf5_write_data(my_hdf5, text72, flottV1)
      deallocate(flottV1)
    else if(indico.eq.3) then
      call hdf5_write_data(my_hdf5, text72, text32V1)
      deallocate(text32V1)
    else if(indico.eq.4) then
      call hdf5_write_data(my_hdf5, text72, dflottV1)
      deallocate(dflottV1)
    endif
  else if(text4.eq.'GREP') then
    ! grep a single value in a rank 1 dataset.
    call REDGET(indprt,nitma,flott,text72,dflott)
    if(indprt.ne.3) call XABORT('DRVUH5: dataset name expected.')
    call hdf5_info(my_hdf5,text72,rank,type,nbyte,dimsr)
    if(rank.ne.1) call XABORT('DRVUH5: rank 1 dataset expected.')
    call REDGET(indic,index,flott,text12,dflott)
    if(indic.lt.0) then
      index=1
    else if(indic.eq.1) then
      call REDGET(indic,nitma,flott,text12,dflott)
      if(indic.ge.0) call XABORT('DRVUH5: >>...<< expected.')
    else
      call XABORT('DRVUH5: integer value or >>...<< expected.')
    endif
    write(6,'(/19h DRVUF5: grep value,i8,12h in dataset ,a,1h:)') index,text72
    if(index.gt.dimsr(1)) call XABORT('DRVUH5: index overflow.')
    indic=-indic
    if(indic.ne.type) then
      write(hsmg,'(33hDRVUH5: inconststent REDPUT type=,i2,14h dataset type=,i2,1h.)') &
      & indic,type
      call XABORT(hsmg)
    endif 
    if(type.eq.1) then
        call hdf5_read_data(my_hdf5,text72,nitmaV1)
        call REDPUT(indic,nitmav1(index),flott,text12,dflott)
        deallocate(nitmaV1)
    else if(type.eq.2) then
        call hdf5_read_data(my_hdf5,text72,flottV1)
        call REDPUT(indic,nitma,flottv1(index),text12,dflott)
        deallocate(flottV1)
    else if(type.eq.3) then
        if(nbyte.le.32) then
          call hdf5_read_data(my_hdf5,text72,text32V1)
          call REDPUT(indic,nitma,flott,text32v1(index),dflott)
          deallocate(text32V1)
        else
          call hdf5_read_data(my_hdf5,text72,text64V1)
          call REDPUT(indic,nitma,flott,text64v1(index),dflott)
          deallocate(text64V1)
        endif
    else if(type.eq.4) then
        call hdf5_read_data(my_hdf5,text72,dflottV1)
        call REDPUT(indic,nitma,flott,text12,dflottv1(index))
        deallocate(dflottV1)
    else
      write(hsmg,100) type,rank
      call XABORT(hsmg)
    endif
  else if(text4.eq.';') then
    return
  else
    write(hsmg,'(8hDRVUH5: ,a4,30h is an invalid utility action.)') text4
    call XABORT(hsmg)
  endif
  go to 10
  !
  100 format(12hDRVUF5: type,i3,9h and rank,i3,19h are not supported.)
end subroutine DRVUH5
