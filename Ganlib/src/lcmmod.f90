!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for lcm -- part 2.
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
module LCMMOD
   use LCMAUX
   private
   public :: LCMPUT, LCMGET, LCMPDL, LCMGDL
   interface LCMPUT
      ! store a record in an associative table
      MODULE PROCEDURE LCMPUT_I0, LCMPUT_R0, LCMPUT_D0, LCMPUT_L0, LCMPUT_C0, &
      LCMPUT_I1, LCMPUT_R1, LCMPUT_D1, LCMPUT_L1, LCMPUT_C1, &
      LCMPUT_I2, LCMPUT_R2, LCMPUT_D2, LCMPUT_L2, LCMPUT_C2, &
      LCMPUT_I3, LCMPUT_R3, LCMPUT_D3, LCMPUT_L3, LCMPUT_C3, &
      LCMPUT_I4, LCMPUT_R4, LCMPUT_D4, LCMPUT_L4, LCMPUT_C4
   end interface
   interface LCMGET
      ! recover a record from an associative table
      MODULE PROCEDURE LCMGET_I0, LCMGET_R0, LCMGET_D0, LCMGET_L0, LCMGET_C0, &
      LCMGET_I1, LCMGET_R1, LCMGET_D1, LCMGET_L1, LCMGET_C1, &
      LCMGET_I2, LCMGET_R2, LCMGET_D2, LCMGET_L2, LCMGET_C2, &
      LCMGET_I3, LCMGET_R3, LCMGET_D3, LCMGET_L3, LCMGET_C3, &
      LCMGET_I4, LCMGET_R4, LCMGET_D4, LCMGET_L4, LCMGET_C4
   end interface
   interface LCMPDL
      ! store a record in an heterogeneous list
      MODULE PROCEDURE LCMPDL_I0, LCMPDL_R0, LCMPDL_D0, LCMPDL_L0, LCMPDL_C0, &
      LCMPDL_I1, LCMPDL_R1, LCMPDL_D1, LCMPDL_L1, LCMPDL_C1, &
      LCMPDL_I2, LCMPDL_R2, LCMPDL_D2, LCMPDL_L2, LCMPDL_C2, &
      LCMPDL_I3, LCMPDL_R3, LCMPDL_D3, LCMPDL_L3, LCMPDL_C3
   end interface
   interface LCMGDL
      ! recover a record from an heterogeneous list
      MODULE PROCEDURE LCMGDL_I0, LCMGDL_R0, LCMGDL_D0, LCMGDL_L0, LCMGDL_C0, &
      LCMGDL_I1, LCMGDL_R1, LCMGDL_D1, LCMGDL_L1, LCMGDL_C1, &
      LCMGDL_I2, LCMGDL_R2, LCMGDL_D2, LCMGDL_L2, LCMGDL_C2, &
      LCMGDL_I3, LCMGDL_R3, LCMGDL_D3, LCMGDL_L3, LCMGDL_C3
   end interface
contains
subroutine LCMPUT_I0(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   integer, target :: idata
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if((itype.ne.1).and.(itype.ne.3)) call XABORT('LCMPUT_I0: type 1 or 3 expected.')
   pt_data=c_loc(idata)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_I0
!
subroutine LCMPUT_R0(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   real, target :: idata
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.2) call XABORT('LCMPUT_R0: type 2 expected.')
   pt_data=c_loc(idata)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_R0
!
subroutine LCMPUT_D0(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   double precision, target :: idata
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.4) call XABORT('LCMPUT_D0: type 4 expected.')
   pt_data=c_loc(idata)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_D0
!
subroutine LCMPUT_L0(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   logical, target :: idata
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.5) call XABORT('LCMPUT_L0: type 5 expected.')
   pt_data=c_loc(idata)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_L0
!
subroutine LCMPUT_C0(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   complex, target :: idata
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.6) call XABORT('LCMPUT_C0: type 6 expected.')
   pt_data=c_loc(idata)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_C0
!
subroutine LCMPUT_I1(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   integer, target, dimension(:) :: idata
   integer, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if((itype.ne.1).and.(itype.ne.3)) call XABORT('LCMPUT_I1: type 1 or 3 expected.')
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_I1
!
subroutine LCMPUT_R1(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   real, target, dimension(:) :: idata
   real, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.2) call XABORT('LCMPUT_R1: type 2 expected.')
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_R1
!
subroutine LCMPUT_D1(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   double precision, target, dimension(:) :: idata
   double precision, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.4) call XABORT('LCMPUT_D1: type 4 expected.')
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_D1
!
subroutine LCMPUT_L1(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   logical, target, dimension(:) :: idata
   logical, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.5) call XABORT('LCMPUT_L1: type 5 expected.')
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_L1
!
subroutine LCMPUT_C1(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   complex, target, dimension(:) :: idata
   complex, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.6) call XABORT('LCMPUT_C1: type 6 expected.')
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_C1
!
subroutine LCMPUT_I2(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   integer, target, dimension(:,:) :: idata
   integer, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if((itype.ne.1).and.(itype.ne.3)) call XABORT('LCMPUT_I2: type 1 or 3 expected.')
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_I2
!
subroutine LCMPUT_R2(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   real, target, dimension(:,:) :: idata
   real, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.2) call XABORT('LCMPUT_R2: type 2 expected.')
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_R2
!
subroutine LCMPUT_D2(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   double precision, target, dimension(:,:) :: idata
   double precision, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.4) call XABORT('LCMPUT_D2: type 4 expected.')
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_D2
!
subroutine LCMPUT_L2(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   logical, target, dimension(:,:) :: idata
   logical, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.5) call XABORT('LCMPUT_L2: type 5 expected.')
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_L2
!
subroutine LCMPUT_C2(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   complex, target, dimension(:,:) :: idata
   complex, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.6) call XABORT('LCMPUT_C2: type 6 expected.')
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_C2
!
subroutine LCMPUT_I3(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   integer, target, dimension(:,:,:) :: idata
   integer, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if((itype.ne.1).and.(itype.ne.3)) call XABORT('LCMPUT_I3: type 1 or 3 expected.')
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_I3
!
subroutine LCMPUT_R3(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   real, target, dimension(:,:,:) :: idata
   real, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.2) call XABORT('LCMPUT_R3: type 2 expected.')
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_R3
!
subroutine LCMPUT_D3(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   double precision, target, dimension(:,:,:) :: idata
   double precision, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.4) call XABORT('LCMPUT_D3: type 4 expected.')
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_D3
!
subroutine LCMPUT_L3(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   logical, target, dimension(:,:,:) :: idata
   logical, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.5) call XABORT('LCMPUT_L3: type 5 expected.')
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_L3
!
subroutine LCMPUT_C3(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   complex, target, dimension(:,:,:) :: idata
   complex, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.6) call XABORT('LCMPUT_C3: type 6 expected.')
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_C3
!
subroutine LCMPUT_I4(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   integer, target, dimension(:,:,:,:) :: idata
   integer, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if((itype.ne.1).and.(itype.ne.3)) call XABORT('LCMPUT_I4: type 1 or 3 expected.')
   idata_p => idata(1,1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_I4
!
subroutine LCMPUT_R4(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   real, target, dimension(:,:,:,:) :: idata
   real, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.2) call XABORT('LCMPUT_R4: type 2 expected.')
   idata_p => idata(1,1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_R4
!
subroutine LCMPUT_D4(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   double precision, target, dimension(:,:,:,:) :: idata
   double precision, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.4) call XABORT('LCMPUT_D4: type 4 expected.')
   idata_p => idata(1,1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_D4
!
subroutine LCMPUT_L4(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   logical, target, dimension(:,:,:,:) :: idata
   logical, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.5) call XABORT('LCMPUT_L4: type 5 expected.')
   idata_p => idata(1,1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_L4
!
subroutine LCMPUT_C4(iplist, name, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   complex, target, dimension(:,:,:,:) :: idata
   complex, pointer :: idata_p
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   if(itype.ne.6) call XABORT('LCMPUT_C4: type 6 expected.')
   idata_p => idata(1,1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT_C4
!
subroutine LCMGET_I0(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer, target :: idata
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   pt_data=c_loc(idata)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_I0
!
subroutine LCMGET_R0(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   real, target :: idata
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   pt_data=c_loc(idata)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_R0
!
subroutine LCMGET_D0(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   double precision, target :: idata
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   pt_data=c_loc(idata)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_D0
!
subroutine LCMGET_L0(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   logical, target :: idata
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   pt_data=c_loc(idata)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_L0
!
subroutine LCMGET_C0(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   complex, target :: idata
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   pt_data=c_loc(idata)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_C0
!
subroutine LCMGET_I1(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer, target, dimension(:) :: idata
   integer, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_I1
!
subroutine LCMGET_R1(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   real, target, dimension(:) :: idata
   real, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_R1
!
subroutine LCMGET_D1(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   double precision, target, dimension(:) :: idata
   double precision, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_D1
!
subroutine LCMGET_L1(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   logical, target, dimension(:) :: idata
   logical, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_L1
!
subroutine LCMGET_C1(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   complex, target, dimension(:) :: idata
   complex, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_C1
!
subroutine LCMGET_I2(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer, target, dimension(:,:) :: idata
   integer, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_I2
!
subroutine LCMGET_R2(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   real, target, dimension(:,:) :: idata
   real, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_R2
!
subroutine LCMGET_D2(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   double precision, target, dimension(:,:) :: idata
   double precision, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_D2
!
subroutine LCMGET_L2(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   logical, target, dimension(:,:) :: idata
   logical, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_L2
!
subroutine LCMGET_C2(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   complex, target, dimension(:,:) :: idata
   complex, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_C2
!
subroutine LCMGET_I3(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer, target, dimension(:,:,:) :: idata
   integer, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_I3
!
subroutine LCMGET_R3(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   real, target, dimension(:,:,:) :: idata
   real, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_R3
!
subroutine LCMGET_D3(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   double precision, target, dimension(:,:,:) :: idata
   double precision, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_D3
!
subroutine LCMGET_L3(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   logical, target, dimension(:,:,:) :: idata
   logical, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_L3
!
subroutine LCMGET_C3(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   complex, target, dimension(:,:,:) :: idata
   complex, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_C3
!
subroutine LCMGET_I4(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer, target, dimension(:,:,:,:) :: idata
   integer, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_I4
!
subroutine LCMGET_R4(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   real, target, dimension(:,:,:,:) :: idata
   real, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_R4
!
subroutine LCMGET_D4(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   double precision, target, dimension(:,:,:,:) :: idata
   double precision, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_D4
!
subroutine LCMGET_L4(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   logical, target, dimension(:,:,:,:) :: idata
   logical, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_L4
!
subroutine LCMGET_C4(iplist, name, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   complex, target, dimension(:,:,:,:) :: idata
   complex, pointer :: idata_p
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   idata_p => idata(1,1,1,1)
   pt_data=c_loc(idata_p)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET_C4
!
subroutine LCMPDL_I0(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   integer, target :: idata
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if((itype.ne.1).and.(itype.ne.3)) call XABORT('LCMPDL_I0: type 1 or 3 expected.')
   pt_data=c_loc(idata)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_I0
!
subroutine LCMPDL_R0(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   real, target :: idata
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.2) call XABORT('LCMPDL_R0: type 2 expected.')
   pt_data=c_loc(idata)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_R0
!
subroutine LCMPDL_D0(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   double precision, target :: idata
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.4) call XABORT('LCMPDL_D0: type 4 expected.')
   pt_data=c_loc(idata)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_D0
!
subroutine LCMPDL_L0(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   logical, target :: idata
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.5) call XABORT('LCMPDL_L0: type 5 expected.')
   pt_data=c_loc(idata)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_L0
!
subroutine LCMPDL_C0(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   complex, target :: idata
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.6) call XABORT('LCMPDL_C0: type 6 expected.')
   pt_data=c_loc(idata)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_C0
!
subroutine LCMPDL_I1(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   integer, target, dimension(:) :: idata
   integer, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if((itype.ne.1).and.(itype.ne.3)) call XABORT('LCMPDL_I1: type 1 or 3 expected.')
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_I1
!
subroutine LCMPDL_R1(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   real, target, dimension(:) :: idata
   real, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.2) call XABORT('LCMPDL_R1: type 2 expected.')
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_R1
!
subroutine LCMPDL_D1(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   double precision, target, dimension(:) :: idata
   double precision, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.4) call XABORT('LCMPDL_D1: type 4 expected.')
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_D1
!
subroutine LCMPDL_L1(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   logical, target, dimension(:) :: idata
   logical, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.5) call XABORT('LCMPDL_L1: type 5 expected.')
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_L1
!
subroutine LCMPDL_C1(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   complex, target, dimension(:) :: idata
   complex, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.6) call XABORT('LCMPDL_C1: type 6 expected.')
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_C1
!
subroutine LCMPDL_I2(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   integer, target, dimension(:,:) :: idata
   integer, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if((itype.ne.1).and.(itype.ne.3)) call XABORT('LCMPDL_I2: type 1 or 3 expected.')
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_I2
!
subroutine LCMPDL_R2(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   real, target, dimension(:,:) :: idata
   real, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.2) call XABORT('LCMPDL_R2: type 2 expected.')
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_R2
!
subroutine LCMPDL_D2(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   double precision, target, dimension(:,:) :: idata
   double precision, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.4) call XABORT('LCMPDL_D2: type 4 expected.')
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_D2
!
subroutine LCMPDL_L2(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   logical, target, dimension(:,:) :: idata
   logical, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.5) call XABORT('LCMPDL_L2: type 5 expected.')
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_L2
!
subroutine LCMPDL_C2(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   complex, target, dimension(:,:) :: idata
   complex, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.6) call XABORT('LCMPDL_C2: type 6 expected.')
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_C2
!
subroutine LCMPDL_I3(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   integer, target, dimension(:,:,:) :: idata
   integer, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if((itype.ne.1).and.(itype.ne.3)) call XABORT('LCMPDL_I3: type 1 or 3 expected.')
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_I3
!
subroutine LCMPDL_R3(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   real, target, dimension(:,:,:) :: idata
   real, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.2) call XABORT('LCMPDL_R3: type 2 expected.')
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_R3
!
subroutine LCMPDL_D3(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   double precision, target, dimension(:,:,:) :: idata
   double precision, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.4) call XABORT('LCMPDL_D3: type 4 expected.')
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_D3
!
subroutine LCMPDL_L3(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   logical, target, dimension(:,:,:) :: idata
   logical, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.5) call XABORT('LCMPDL_L3: type 5 expected.')
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_L3
!
subroutine LCMPDL_C3(iplist, ipos, ilong, itype, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   complex, target, dimension(:,:,:) :: idata
   complex, pointer :: idata_p
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   if(itype.ne.6) call XABORT('LCMPDL_C3: type 6 expected.')
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL_C3
!
subroutine LCMGDL_I0(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   integer, target :: idata
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   pt_data=c_loc(idata)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_I0
!
subroutine LCMGDL_R0(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   real, target :: idata
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   pt_data=c_loc(idata)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_R0
!
subroutine LCMGDL_D0(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   double precision, target :: idata
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   pt_data=c_loc(idata)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_D0
!
subroutine LCMGDL_L0(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   logical, target :: idata
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   pt_data=c_loc(idata)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_L0
!
subroutine LCMGDL_C0(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   complex, target :: idata
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   pt_data=c_loc(idata)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_C0
!
subroutine LCMGDL_I1(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   integer, target, dimension(:) :: idata
   integer, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_I1
!
subroutine LCMGDL_R1(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   real, target, dimension(:) :: idata
   real, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_R1
!
subroutine LCMGDL_D1(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   double precision, target, dimension(:) :: idata
   double precision, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_D1
!
subroutine LCMGDL_L1(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   logical, target, dimension(:) :: idata
   logical, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_L1
!
subroutine LCMGDL_C1(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   complex, target, dimension(:) :: idata
   complex, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_C1
!
subroutine LCMGDL_I2(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   integer, target, dimension(:,:) :: idata
   integer, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_I2
!
subroutine LCMGDL_R2(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   real, target, dimension(:,:) :: idata
   real, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_R2
!
subroutine LCMGDL_D2(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   double precision, target, dimension(:,:) :: idata
   double precision, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_D2
!
subroutine LCMGDL_L2(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   logical, target, dimension(:,:) :: idata
   logical, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_L2
!
subroutine LCMGDL_C2(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   complex, target, dimension(:,:) :: idata
   complex, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1,1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_C2
!
subroutine LCMGDL_I3(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   integer, target, dimension(:,:,:) :: idata
   integer, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_I3
!
subroutine LCMGDL_R3(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   real, target, dimension(:,:,:) :: idata
   real, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_R3
!
subroutine LCMGDL_D3(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   double precision, target, dimension(:,:,:) :: idata
   double precision, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_D3
!
subroutine LCMGDL_L3(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   logical, target, dimension(:,:,:) :: idata
   logical, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_L3
!
subroutine LCMGDL_C3(iplist, ipos, idata)
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   complex, target, dimension(:,:,:) :: idata
   complex, pointer :: idata_p
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   idata_p => idata(1,1,1)
   pt_data=c_loc(idata_p)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL_C3
end module LCMMOD
