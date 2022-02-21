!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for hdf5.
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
!-----------------------------------------------------------------------
!
module hdf5_wrap
   use, intrinsic :: iso_c_binding
   private
   integer, parameter :: MAX_NAME = 1024
   public :: hdf5_open_file,hdf5_close_file,hdf5_get_dimensions,hdf5_get_shape, &
             hdf5_list,hdf5_info,hdf5_read_data,hdf5_write_data,hdf5_list_datasets, &
             hdf5_list_groups,hdf5_group_exists

   interface hdf5_read_data
     module procedure hdf5_read_data_0d_int4,  hdf5_read_data_1d_int4,  &
                      hdf5_read_data_2d_int4,  hdf5_read_data_0d_real4, &
                      hdf5_read_data_1d_real4, hdf5_read_data_2d_real4, &
                      hdf5_read_data_3d_real4, hdf5_read_data_4d_real4, &
                      hdf5_read_data_0d_real8, hdf5_read_data_1d_real8, &
                      hdf5_read_data_2d_real8, hdf5_read_data_3d_real8, &
                      hdf5_read_data_4d_real8, hdf5_read_data_0d_string,&
                      hdf5_read_data_1d_string
   end interface hdf5_read_data
   interface hdf5_write_data
     module procedure hdf5_write_data_0d_int4, hdf5_write_data_1d_int4,   &
                      hdf5_write_data_2d_int4, hdf5_write_data_0d_real4,  &
                      hdf5_write_data_1d_real4, hdf5_write_data_2d_real4, &
                      hdf5_write_data_3d_real4, hdf5_write_data_4d_real4, &
                      hdf5_write_data_0d_real8, hdf5_write_data_1d_real8, &
                      hdf5_write_data_2d_real8, hdf5_write_data_3d_real8, &
                      hdf5_write_data_4d_real8, hdf5_write_data_0d_string,&
                      hdf5_write_data_1d_string
   end interface hdf5_write_data

   character(len=131) :: hsmg
   character(kind=c_char), dimension(MAX_NAME) :: name1024

contains
subroutine STRCUT(name1, name2)
   ! transform a Fortran string into a C null-terminated string
   character(kind=c_char), dimension(*) :: name1
   character(len=*) :: name2
   integer :: ilong
   interface
      subroutine strcut_c (s, ct, n) bind(c)
         use, intrinsic :: iso_c_binding
         character(kind=c_char), dimension(*) :: s, ct
         integer(c_int), value :: n
      end subroutine strcut_c
   end interface
   ilong=len(trim(adjustl(name2)))
   call strcut_c(name1, trim(adjustl(name2)), ilong)
end subroutine STRCUT
!
subroutine STRFIL(name1, name2, nbyte)
   ! transform a C null-terminated string into a Fortran string
   character(len=*) :: name1
   character(kind=c_char), dimension(*) :: name2
   integer :: nbyte, ilong
   interface
      subroutine strfil_c (s, ct, n) bind(c)
         use, intrinsic :: iso_c_binding
         character(kind=c_char), dimension(*) :: s, ct
         integer(c_int), value :: n
      end subroutine strfil_c
   end interface
   ilong=len(name1)
   name1=' '
   call strfil_c(name1, name2, min(ilong,nbyte))
end subroutine STRFIL
!
subroutine STRFIL1D(name1, name2, nbyte)
   ! transform a C null-terminated string into a Fortran string array
   character(len=*),dimension(:) :: name1
   character(kind=c_char), dimension(*) :: name2
   integer :: nbyte, ilong
   character(len=MAX_NAME) :: text1024
   interface
      subroutine strfil_c (s, ct, n) bind(c)
         use, intrinsic :: iso_c_binding
         character(kind=c_char), dimension(*) :: s, ct
         integer(c_int), value :: n
      end subroutine strfil_c
   end interface
   ilong=len(name1)
   idim=size(name1,1)
   iof=1
   do i=1,idim
     name1(i)=' '
     call strfil_c(text1024, name2(iof), nbyte)
     name1(i)=text1024(:min(ilong,nbyte))
     iof=iof+nbyte
   enddo
end subroutine STRFIL1D
!
subroutine hdf5_open_file(fname, ifile, rdonly)
   ! open a HDF5 file
   character(len=*), intent(in)  :: fname
   type(c_ptr), intent(out) :: ifile
   logical, optional :: rdonly
   interface
      subroutine hdf5_open_file_c (fname, ifile, irdonly) bind(c)
         use, intrinsic :: iso_c_binding
         character(kind=c_char), dimension(*) :: fname
         type(c_ptr) :: ifile
         integer(c_int), value :: irdonly
      end subroutine hdf5_open_file_c
   end interface
   !
   irdonly=0
   if(present(rdonly)) then
     if(rdonly) irdonly=1
   endif
   call STRCUT(name1024, fname)
   call hdf5_open_file_c (name1024, ifile, irdonly)
end subroutine hdf5_open_file
!
subroutine hdf5_close_file(ifile)
   ! close a HDF5 file
   type(c_ptr),intent(in) :: ifile
   interface
      subroutine hdf5_close_file_c (ifile) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
      end subroutine hdf5_close_file_c
   end interface
   call hdf5_close_file_c(ifile)
end subroutine hdf5_close_file
!
subroutine hdf5_list(ifile, name)
   ! table of contents
   type(c_ptr) :: ifile
   character(len=*),intent(in) :: name
   !
   interface
      subroutine hdf5_list_c (ifile, namp) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
      end subroutine hdf5_list_c
   end interface
   call STRCUT(name1024, name)
   call hdf5_list_c(ifile, name1024)
end subroutine hdf5_list
!
subroutine hdf5_info(ifile, name, rank, type, nbyte, dimsr)
   ! find dataset info
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   integer,intent(out) :: rank, type, nbyte
   integer,target,dimension(5),intent(out) :: dimsr
   !
   type(c_ptr) :: pt_dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   dimsr(:5)=0
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
end subroutine hdf5_info
!
integer function hdf5_get_dimensions(ifile, name)
   ! find dataset rank (number of dimensions
   type(c_ptr),intent(in) :: ifile
   character(len=*), intent(in) :: name
   integer :: rank
   !
   interface
      subroutine hdf5_get_dimensions_c(ifile, namp, rank) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank
      end subroutine hdf5_get_dimensions_c
   end interface
   call STRCUT(name1024, name)
   call hdf5_get_dimensions_c(ifile, name1024, rank)
   hdf5_get_dimensions = rank
end function hdf5_get_dimensions
!
subroutine hdf5_get_shape(ifile, name, dimsr)
   ! find dataset shape
   type(c_ptr), intent(in) :: ifile
   character(len=*), intent(in) :: name
   integer, allocatable, target :: dimsr(:)
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr
   !
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   rank=hdf5_get_dimensions(ifile, name)
   allocate(dimsr(rank))
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
end subroutine hdf5_get_shape
!
subroutine hdf5_list_datasets(ifile, name, dsets)
   ! collect daughter dataset names in a group
   type(c_ptr), intent(in) :: ifile
   character(len=*), intent(in) :: name
   character(len=*), allocatable, dimension(:) :: dsets
   integer :: nbobj, ndsets
   !
   character(kind=c_char), allocatable, dimension(:) :: pt_strim
   interface
      subroutine hdf5_get_num_group_c(ifile, namp, nbobj) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: nbobj
      end subroutine hdf5_get_num_group_c
   end interface
   interface
      subroutine hdf5_list_datasets_c(ifile, namp, ndsets, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp, idata
         integer(c_int) :: ndsets
      end subroutine hdf5_list_datasets_c
   end interface
   !
   call STRCUT(name1024, name)
   call hdf5_get_num_group_c (ifile, name1024, nbobj)
   if (nbobj .eq. 0) return
   allocate(pt_strim(nbobj*MAX_NAME))
   call hdf5_list_datasets_c (ifile, name1024, ndsets, pt_strim)
   allocate(dsets(ndsets))
   call STRFIL1D(dsets,pt_strim,MAX_NAME)
   deallocate(pt_strim)
end subroutine hdf5_list_datasets
!
subroutine hdf5_list_groups(ifile, name, groups)
   ! collect daughter group names in a group
   type(c_ptr), intent(in) :: ifile
   character(len=*), intent(in) :: name
   character(len=*), allocatable, dimension(:) :: groups
   integer :: nbobj, ngroups
   !
   character(kind=c_char), allocatable, dimension(:) :: pt_strim
   interface
      subroutine hdf5_get_num_group_c(ifile, namp, nbobj) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: nbobj
      end subroutine hdf5_get_num_group_c
   end interface
   interface
      subroutine hdf5_list_groups_c(ifile, namp, ngroups, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp, idata
         integer(c_int) :: ngroups
      end subroutine hdf5_list_groups_c
   end interface
   !
   call STRCUT(name1024, name)
   call hdf5_get_num_group_c (ifile, name1024, nbobj)
   if (nbobj .eq. 0) return
   allocate(pt_strim(nbobj*MAX_NAME))
   call hdf5_list_groups_c (ifile, name1024, ngroups, pt_strim)
   allocate(groups(ngroups))
   call STRFIL1D(groups,pt_strim,MAX_NAME)
   deallocate(pt_strim)
end subroutine hdf5_list_groups
!
function hdf5_group_exists(ifile, name) result(lexist)
   ! test for existence of a group
   type(c_ptr), intent(in) :: ifile
   character(len=*), intent(in) :: name
   logical :: lexist
   interface
      function hdf5_group_exists_c (ifile, namp) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) :: hdf5_group_exists_c
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
      end function hdf5_group_exists_c
   end interface
   call STRCUT(name1024, name)
   lexist = (hdf5_group_exists_c(ifile, name1024) == 0)
end function hdf5_group_exists
!
subroutine hdf5_read_data_0d_int4(ifile, name, idata)
   ! read a rank 0 integer dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   integer, target :: idata
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr, pt_data
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_int_c(ifile, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine hdf5_read_data_int_c
   end interface
   pt_data=c_loc(idata)
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if((rank.ne.1).and.(dimsr(1).ne.1)) then
     write(hsmg,'(49hhdf5_read_data_0d_int4: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   call hdf5_read_data_int_c(ifile, name1024, pt_data)
end subroutine hdf5_read_data_0d_int4
!
subroutine hdf5_read_data_1d_int4(ifile, name, idata)
   ! read a rank 1 integer dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   integer, allocatable, dimension(:), target :: idata
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr, pt_data
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_int_c(ifile, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine hdf5_read_data_int_c
   end interface
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if(rank.ne.1) then
     write(hsmg,'(49hhdf5_read_data_1d_int4: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   allocate(idata(dimsr(1)))
   pt_data=c_loc(idata)
   call hdf5_read_data_int_c(ifile, name1024, pt_data)
end subroutine hdf5_read_data_1d_int4
!
subroutine hdf5_read_data_2d_int4(ifile, name, idata)
   ! read a rank 2 integer dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   integer, allocatable, dimension(:,:), target :: idata
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr, pt_data
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_int_c(ifile, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine hdf5_read_data_int_c
   end interface
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if(rank.ne.2) then
     write(hsmg,'(49hhdf5_read_data_2d_int4: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   allocate(idata(dimsr(1),dimsr(2)))
   pt_data=c_loc(idata)
   call hdf5_read_data_int_c(ifile, name1024, pt_data)
end subroutine hdf5_read_data_2d_int4
!
subroutine hdf5_read_data_0d_real4(ifile, name, rdata)
   ! read a rank 0 real4 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(4), target :: rdata
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr, pt_data
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_real4_c(ifile, namp, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: rdata
      end subroutine hdf5_read_data_real4_c
   end interface
   pt_data=c_loc(rdata)
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if((rank.ne.1).and.(dimsr(1).ne.1)) then
     write(hsmg,'(50hhdf5_read_data_0d_real4: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   call hdf5_read_data_real4_c(ifile, name1024, pt_data)
end subroutine hdf5_read_data_0d_real4
!
subroutine hdf5_read_data_1d_real4(ifile, name, rdata)
   ! read a rank 1 real4 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(4), allocatable, dimension(:), target :: rdata
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr, pt_data
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_real4_c(ifile, namp, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: rdata
      end subroutine hdf5_read_data_real4_c
   end interface
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if(rank.ne.1) then
     write(hsmg,'(50hhdf5_read_data_1d_real4: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   allocate(rdata(dimsr(1)))
   pt_data=c_loc(rdata)
   call hdf5_read_data_real4_c(ifile, name1024, pt_data)
end subroutine hdf5_read_data_1d_real4
!
subroutine hdf5_read_data_2d_real4(ifile, name, rdata)
   ! read a rank 2 real4 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(4), allocatable, dimension(:,:), target :: rdata
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr, pt_data
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_real4_c(ifile, namp, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: rdata
      end subroutine hdf5_read_data_real4_c
   end interface
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if(rank.ne.2) then
     write(hsmg,'(50hhdf5_read_data_2d_real4: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   allocate(rdata(dimsr(1),dimsr(2)))
   pt_data=c_loc(rdata)
   call hdf5_read_data_real4_c(ifile, name1024, pt_data)
end subroutine hdf5_read_data_2d_real4
!
subroutine hdf5_read_data_3d_real4(ifile, name, rdata)
   ! read a rank 3 real4 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(4), allocatable, dimension(:,:,:), target :: rdata
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr, pt_data
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_real4_c(ifile, namp, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: rdata
      end subroutine hdf5_read_data_real4_c
   end interface
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if(rank.ne.3) then
     write(hsmg,'(50hhdf5_read_data_3d_real4: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   allocate(rdata(dimsr(1),dimsr(2),dimsr(3)))
   pt_data=c_loc(rdata)
   call hdf5_read_data_real4_c(ifile, name1024, pt_data)
end subroutine hdf5_read_data_3d_real4
!
subroutine hdf5_read_data_4d_real4(ifile, name, rdata)
   ! read a rank 4 real4 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(4), allocatable, dimension(:,:,:,:), target :: rdata
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr, pt_data
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_real4_c(ifile, namp, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: rdata
      end subroutine hdf5_read_data_real4_c
   end interface
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if(rank.ne.4) then
     write(hsmg,'(50hhdf5_read_data_4d_real4: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   allocate(rdata(dimsr(1),dimsr(2),dimsr(3),dimsr(4)))
   pt_data=c_loc(rdata)
   call hdf5_read_data_real4_c(ifile, name1024, pt_data)
end subroutine hdf5_read_data_4d_real4
!
subroutine hdf5_read_data_0d_real8(ifile, name, rdata)
   ! read a rank 0 real8 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(8), target :: rdata
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr, pt_data
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_real8_c(ifile, namp, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: rdata
      end subroutine hdf5_read_data_real8_c
   end interface
   pt_data=c_loc(rdata)
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if((rank.ne.1).and.(dimsr(1).ne.1)) then
     write(hsmg,'(50hhdf5_read_data_0d_real8: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   call hdf5_read_data_real8_c(ifile, name1024, pt_data)
end subroutine hdf5_read_data_0d_real8
!
subroutine hdf5_read_data_1d_real8(ifile, name, rdata)
   ! read a rank 1 real8 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(8), allocatable, dimension(:), target :: rdata
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr, pt_data
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_real8_c(ifile, namp, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: rdata
      end subroutine hdf5_read_data_real8_c
   end interface
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if(rank.ne.1) then
     write(hsmg,'(50hhdf5_read_data_1d_real8: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   allocate(rdata(dimsr(1)))
   pt_data=c_loc(rdata)
   call hdf5_read_data_real8_c(ifile, name1024, pt_data)
end subroutine hdf5_read_data_1d_real8
!
subroutine hdf5_read_data_2d_real8(ifile, name, rdata)
   ! read a rank 2 real8 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(8), allocatable, dimension(:,:), target :: rdata
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr, pt_data
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_real8_c(ifile, namp, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: rdata
      end subroutine hdf5_read_data_real8_c
   end interface
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if(rank.ne.2) then
     write(hsmg,'(50hhdf5_read_data_2d_real8: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   allocate(rdata(dimsr(1),dimsr(2)))
   pt_data=c_loc(rdata)
   call hdf5_read_data_real8_c(ifile, name1024, pt_data)
end subroutine hdf5_read_data_2d_real8
!
subroutine hdf5_read_data_3d_real8(ifile, name, rdata)
   ! read a rank 3 real8 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(8), allocatable, dimension(:,:,:), target :: rdata
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr, pt_data
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_real8_c(ifile, namp, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: rdata
      end subroutine hdf5_read_data_real8_c
   end interface
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if(rank.ne.3) then
     write(hsmg,'(50hhdf5_read_data_3d_real8: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   allocate(rdata(dimsr(1),dimsr(2),dimsr(3)))
   pt_data=c_loc(rdata)
   call hdf5_read_data_real8_c(ifile, name1024, pt_data)
end subroutine hdf5_read_data_3d_real8
!
subroutine hdf5_read_data_4d_real8(ifile, name, rdata)
   ! read a rank 4 real8 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(8), allocatable, dimension(:,:,:,:), target :: rdata
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr, pt_data
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_real8_c(ifile, namp, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: rdata
      end subroutine hdf5_read_data_real8_c
   end interface
   pt_dimsr=c_loc(dimsr)
   call STRCUT(name1024, name)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if(rank.ne.4) then
     write(hsmg,'(50hhdf5_read_data_4d_real8: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   allocate(rdata(dimsr(1),dimsr(2),dimsr(3),dimsr(4)))
   pt_data=c_loc(rdata)
   call hdf5_read_data_real8_c(ifile, name1024, pt_data)
end subroutine hdf5_read_data_4d_real8
!
subroutine hdf5_read_data_0d_string(ifile, name, idata)
   ! read a rank 0 string dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*), intent(in) :: name
   character(len=*), target :: idata
   character(kind=c_char), allocatable, dimension(:) :: pt_data
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr
   !character(kind=c_char) :: pt_data
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_string_c(ifile, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp, idata
      end subroutine hdf5_read_data_string_c
   end interface
   call STRCUT(name1024, name)
   pt_dimsr=c_loc(dimsr)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if((rank.ne.1).and.(dimsr(1).ne.1)) then
     write(hsmg,'(51hhdf5_read_data_0d_string: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   allocate(pt_data(nbyte+1)) ! add one byte for null termination
   call hdf5_read_data_string_c(ifile, name1024, pt_data)
   call STRFIL(idata,pt_data,nbyte)
   deallocate(pt_data)
end subroutine hdf5_read_data_0d_string
!
subroutine hdf5_read_data_1d_string(ifile, name, idata)
   ! read a rank 1 string dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*), intent(in) :: name
   character(len=*), allocatable, dimension(:), target :: idata
   character(kind=c_char), allocatable, dimension(:) :: pt_data
   !
   integer :: rank, type, nbyte
   type(c_ptr) :: pt_dimsr
   integer,target,dimension(5) :: dimsr
   interface
      subroutine hdf5_info_c(ifile, namp, rank, type, nbyte, dimsr) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: rank, type, nbyte
         type(c_ptr), value :: dimsr
      end subroutine hdf5_info_c
   end interface
   interface
      subroutine hdf5_read_data_string_c(ifile, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp, idata
      end subroutine hdf5_read_data_string_c
   end interface
   call STRCUT(name1024, name)
   pt_dimsr=c_loc(dimsr)
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   if(rank.ne.1) then
     write(hsmg,'(51hhdf5_read_data_1d_string: invalid info for dataset ,a)')  name
     call XABORT(hsmg)
   endif
   call hdf5_info_c(ifile, name1024, rank, type, nbyte, pt_dimsr)
   allocate(idata(dimsr(1)))
   allocate(pt_data(dimsr(1)*nbyte+1)) ! add one byte for null termination
   call hdf5_read_data_string_c(ifile, name1024, pt_data)
   call STRFIL1D(idata,pt_data,nbyte)
   deallocate(pt_data)
end subroutine hdf5_read_data_1d_string
!
subroutine hdf5_write_data_0d_int4(ifile, name, idata)
   ! write a rank 1 integer dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   integer,intent(in), target :: idata
   !
   integer, parameter :: rank = 1
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   interface
      subroutine hdf5_write_data_int_c(ifile, namp, rank, dimsf, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank
         type(c_ptr), value :: dimsf, idata
      end subroutine hdf5_write_data_int_c
   end interface
   dimsr(1)=1
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(idata)
   call STRCUT(name1024, name)
   call hdf5_write_data_int_c(ifile, name1024, rank, pt_dimsr, pt_data)
end subroutine hdf5_write_data_0d_int4
!
subroutine hdf5_write_data_1d_int4(ifile, name, idata)
   ! write a rank 1 integer dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   integer,dimension(:),intent(in), target :: idata
   !
   integer, parameter :: rank = 1
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   integer, pointer :: idata_p
   interface
      subroutine hdf5_write_data_int_c(ifile, namp, rank, dimsf, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank
         type(c_ptr), value :: dimsf, idata
      end subroutine hdf5_write_data_int_c
   end interface
   dimsr(1)=size(idata,1)
   idata_p => idata(1)
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(idata_p)
   call STRCUT(name1024, name)
   call hdf5_write_data_int_c(ifile, name1024, rank, pt_dimsr, pt_data)
end subroutine hdf5_write_data_1d_int4
!
subroutine hdf5_write_data_2d_int4(ifile, name, idata)
   ! write a rank 2 integer dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   integer,dimension(:,:),intent(in), target :: idata
   !
   integer, parameter :: rank = 2
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   integer, pointer :: idata_p
   interface
      subroutine hdf5_write_data_int_c(ifile, namp, rank, dimsf, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank
         type(c_ptr), value :: dimsf, idata
      end subroutine hdf5_write_data_int_c
   end interface
   dimsr(2)=size(idata,1)
   dimsr(1)=size(idata,2)
   idata_p => idata(1,1)
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(idata_p)
   call STRCUT(name1024, name)
   call hdf5_write_data_int_c(ifile, name1024, rank, pt_dimsr, pt_data)
end subroutine hdf5_write_data_2d_int4
!
subroutine hdf5_write_data_0d_real4(ifile, name, rdata)
   ! write a rank 0 real4 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(4),intent(in), target :: rdata
   !
   integer, parameter :: rank = 1
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   interface
      subroutine hdf5_write_data_real4_c(ifile, namp, rank, dimsf, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank
         type(c_ptr), value :: dimsf, rdata
      end subroutine hdf5_write_data_real4_c
   end interface
   dimsr(1)=1
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(rdata)
   call STRCUT(name1024, name)
   call hdf5_write_data_real4_c(ifile, name1024, rank, pt_dimsr, pt_data)
end subroutine hdf5_write_data_0d_real4
!
subroutine hdf5_write_data_1d_real4(ifile, name, rdata)
   ! write a rank 1 real4 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(4),dimension(:),intent(in), target :: rdata
   !
   integer, parameter :: rank = 1
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   real, pointer :: rdata_p
   interface
      subroutine hdf5_write_data_real4_c(ifile, namp, rank, dimsf, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank
         type(c_ptr), value :: dimsf, rdata
      end subroutine hdf5_write_data_real4_c
   end interface
   dimsr(1)=size(rdata,1)
   rdata_p => rdata(1)
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(rdata_p)
   call STRCUT(name1024, name)
   call hdf5_write_data_real4_c(ifile, name1024, rank, pt_dimsr, pt_data)
end subroutine hdf5_write_data_1d_real4
!
subroutine hdf5_write_data_2d_real4(ifile, name, rdata)
   ! write a rank 2 real4 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(4),dimension(:,:),intent(in), target :: rdata
   !
   integer, parameter :: rank = 2
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   real, pointer :: rdata_p
   interface
      subroutine hdf5_write_data_real4_c(ifile, namp, rank, dimsf, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank
         type(c_ptr), value :: dimsf, rdata
      end subroutine hdf5_write_data_real4_c
   end interface
   dimsr(2)=size(rdata,1)
   dimsr(1)=size(rdata,2)
   rdata_p => rdata(1,1)
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(rdata_p)
   call STRCUT(name1024, name)
   call hdf5_write_data_real4_c(ifile, name1024, rank, pt_dimsr, pt_data)
end subroutine hdf5_write_data_2d_real4
!
subroutine hdf5_write_data_3d_real4(ifile, name, rdata)
   ! write a rank 3 real4 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(4),dimension(:,:,:),intent(in), target :: rdata
   !
   integer, parameter :: rank = 3
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   real, pointer :: rdata_p
   interface
      subroutine hdf5_write_data_real4_c(ifile, namp, rank, dimsf, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank
         type(c_ptr), value :: dimsf, rdata
      end subroutine hdf5_write_data_real4_c
   end interface
   dimsr(3)=size(rdata,1)
   dimsr(2)=size(rdata,2)
   dimsr(1)=size(rdata,3)
   rdata_p => rdata(1,1,1)
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(rdata_p)
   call STRCUT(name1024, name)
   call hdf5_write_data_real4_c(ifile, name1024, rank, pt_dimsr, pt_data)
end subroutine hdf5_write_data_3d_real4
!
subroutine hdf5_write_data_4d_real4(ifile, name, rdata)
   ! write a rank 4 real4 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(4),dimension(:,:,:,:),intent(in), target :: rdata
   !
   integer, parameter :: rank = 4
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   real, pointer :: rdata_p
   interface
      subroutine hdf5_write_data_real4_c(ifile, namp, rank, dimsf, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank
         type(c_ptr), value :: dimsf, rdata
      end subroutine hdf5_write_data_real4_c
   end interface
   dimsr(4)=size(rdata,1)
   dimsr(3)=size(rdata,2)
   dimsr(2)=size(rdata,3)
   dimsr(1)=size(rdata,4)
   rdata_p => rdata(1,1,1,1)
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(rdata_p)
   call STRCUT(name1024, name)
   call hdf5_write_data_real4_c(ifile, name1024, rank, pt_dimsr, pt_data)
end subroutine hdf5_write_data_4d_real4
!
subroutine hdf5_write_data_0d_real8(ifile, name, rdata)
   ! write a rank 0 real8 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(8),intent(in), target :: rdata
   !
   integer, parameter :: rank = 1
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   interface
      subroutine hdf5_write_data_real8_c(ifile, namp, rank, dimsf, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank
         type(c_ptr), value :: dimsf, rdata
      end subroutine hdf5_write_data_real8_c
   end interface
   dimsr(1)=1
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(rdata)
   call STRCUT(name1024, name)
   call hdf5_write_data_real8_c(ifile, name1024, rank, pt_dimsr, pt_data)
end subroutine hdf5_write_data_0d_real8
!
subroutine hdf5_write_data_1d_real8(ifile, name, rdata)
   ! write a rank 1 real8 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(8),dimension(:),intent(in), target :: rdata
   !
   integer, parameter :: rank = 1
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   real(8), pointer :: rdata_p
   interface
      subroutine hdf5_write_data_real8_c(ifile, namp, rank, dimsf, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank
         type(c_ptr), value :: dimsf, rdata
      end subroutine hdf5_write_data_real8_c
   end interface
   dimsr(1)=size(rdata,1)
   rdata_p => rdata(1)
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(rdata_p)
   call STRCUT(name1024, name)
   call hdf5_write_data_real8_c(ifile, name1024, rank, pt_dimsr, pt_data)
end subroutine hdf5_write_data_1d_real8
!
subroutine hdf5_write_data_2d_real8(ifile, name, rdata)
   ! write a rank 2 real8 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(8),dimension(:,:),intent(in), target :: rdata
   !
   integer, parameter :: rank = 2
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   real(8), pointer :: rdata_p
   interface
      subroutine hdf5_write_data_real8_c(ifile, namp, rank, dimsf, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank
         type(c_ptr), value :: dimsf, rdata
      end subroutine hdf5_write_data_real8_c
   end interface
   dimsr(2)=size(rdata,1)
   dimsr(1)=size(rdata,2)
   rdata_p => rdata(1,1)
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(rdata_p)
   call STRCUT(name1024, name)
   call hdf5_write_data_real8_c(ifile, name1024, rank, pt_dimsr, pt_data)
end subroutine hdf5_write_data_2d_real8
!
subroutine hdf5_write_data_3d_real8(ifile, name, rdata)
   ! write a rank 3 real8 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(8),dimension(:,:,:),intent(in), target :: rdata
   !
   integer, parameter :: rank = 3
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   real(8), pointer :: rdata_p
   interface
      subroutine hdf5_write_data_real8_c(ifile, namp, rank, dimsf, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank
         type(c_ptr), value :: dimsf, rdata
      end subroutine hdf5_write_data_real8_c
   end interface
   dimsr(3)=size(rdata,1)
   dimsr(2)=size(rdata,2)
   dimsr(1)=size(rdata,3)
   rdata_p => rdata(1,1,1)
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(rdata_p)
   call STRCUT(name1024, name)
   call hdf5_write_data_real8_c(ifile, name1024, rank, pt_dimsr, pt_data)
end subroutine hdf5_write_data_3d_real8
!
subroutine hdf5_write_data_4d_real8(ifile, name, rdata)
   ! write a rank 4 real8 dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   real(8),dimension(:,:,:,:),intent(in), target :: rdata
   !
   integer, parameter :: rank = 4
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   real(8), pointer :: rdata_p
   interface
      subroutine hdf5_write_data_real8_c(ifile, namp, rank, dimsf, rdata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank
         type(c_ptr), value :: dimsf, rdata
      end subroutine hdf5_write_data_real8_c
   end interface
   dimsr(4)=size(rdata,1)
   dimsr(3)=size(rdata,2)
   dimsr(2)=size(rdata,3)
   dimsr(1)=size(rdata,4)
   rdata_p => rdata(1,1,1,1)
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(rdata_p)
   call STRCUT(name1024, name)
   call hdf5_write_data_real8_c(ifile, name1024, rank, pt_dimsr, pt_data)
end subroutine hdf5_write_data_4d_real8
!
subroutine hdf5_write_data_0d_string(ifile, name, idata)
   ! write a rank 1 string dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   character(len=*),intent(in), target :: idata
   !
   integer, parameter :: rank = 1
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   character(len=1), pointer :: idata_p
   interface
      subroutine hdf5_write_data_string_c(ifile, namp, rank, lenstr, dimsf, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank, lenstr
         type(c_ptr), value :: dimsf, idata
      end subroutine hdf5_write_data_string_c
   end interface
   dimsr(1)=1
   idata_p => idata
   lenstr=len(idata)
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(idata_p)
   call STRCUT(name1024, name)
   call hdf5_write_data_string_c(ifile, name1024, rank, lenstr, pt_dimsr, pt_data)
end subroutine hdf5_write_data_0d_string
!
subroutine hdf5_write_data_1d_string(ifile, name, idata)
   ! write a rank 1 string dataset
   type(c_ptr),intent(in) :: ifile
   character(len=*),intent(in) :: name
   character(len=*),intent(in), dimension(:), target :: idata
   !
   integer, parameter :: rank = 1
   integer,dimension(rank),target :: dimsr
   type(c_ptr) :: pt_dimsr, pt_data
   character(len=1), pointer :: idata_p
   interface
      subroutine hdf5_write_data_string_c(ifile, namp, rank, lenstr, dimsf, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: ifile
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: rank, lenstr
         type(c_ptr), value :: dimsf, idata
      end subroutine hdf5_write_data_string_c
   end interface
   dimsr(1)=size(idata,1)
   idata_p => idata(1)
   lenstr=len(idata(1))
   pt_dimsr=c_loc(dimsr)
   pt_data=c_loc(idata_p)
   call STRCUT(name1024, name)
   call hdf5_write_data_string_c(ifile, name1024, rank, lenstr, pt_dimsr, pt_data)
end subroutine hdf5_write_data_1d_string
end module hdf5_wrap
