!
!-----------------------------------------------------------------------
!
!Purpose:
! allocate and release file units associated to a given file name. Word
! addressable (KDI), sequential (formatted or not) and direct access
! (DA) files are permitted. These functions are Ganlib wrappers for the
! KDROPN and KDIOP utilities.
! FILOPN:  open file and allocate unit number. Allocate a unit number to
!          to file name. If unit is already opened, returns its address.
! FILCLS:  close file and release unit number.
! FILUNIT: recover Fortran file unit number.
! FILKDI:  recover KDI file c_ptr.
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
!Parameters: input
! iplcm   pointer to the LCM object.
! cuname  filename. if cuname=' ', use a default name
! iactio  action on file
!         =0 to allocate a new file
!         =1 to access and modify an existing file
!         =2 to access an existing file in read-only mode
!         =3 unknown
! iutype  file type
!         =1 KDI word addressable file
!         =2 sequential unformatted
!         =3 sequential formatted
!         =4 direct access (DA) unformatted file
! ldra    number of words in DA a record (required for iutype=4 only)
! file_pt type(c_ptr) address of file
! my_file type(FIL_file) address of file
! iactio  action on file
!         = 1 to keep the file;
!         = 2 to delete the file.
!
!Parameters: output
! FILOPN  type(FIL_file) address of file (successful allocation).
!          = NULL in case of allocation failure.
! FILCLS  error status
!         =  0 unit closed
!         = -1 error on close of kdi file
!         = -2 error on close of Fortran file
!         = -3 unknown file at close
! FILUNIT Fortran file unit number
! FILKDI  KDI file c_ptr
!
!-----------------------------------------------------------------------
!
module FILMOD
   use, intrinsic :: iso_c_binding
   type FIL_file
      integer :: unit
      type(c_ptr) :: kdi_file
   end type FIL_file
contains
   function FILOPN(cuname,iactio,iutype,lrda) result(my_file)
      use, intrinsic :: iso_c_binding
      character(len=*) :: cuname
      integer :: iactio,iutype,lrda
      !----
      !   local variables
      !----
      type(FIL_file), pointer :: my_file
      type(c_ptr) :: my_kdi_file
      type(c_ptr),external :: KDIOP
      integer,external :: KDROPN
      integer :: ret_val
      !----
      !  kdi (word addressable) file open
      !----
      if(iutype == 1) then
         my_kdi_file=KDIOP(crdnam,iactio)
         if(.not.c_associated(my_kdi_file)) go to 6000
         allocate(my_file)
         my_file%kdi_file=my_kdi_file
         my_file%unit=0
      else
      !----
      !  Fortran file open
      !----
         ret_val=KDROPN(cuname,iactio,iutype,lrda)
         if(ret_val <= 0) go to 6000
         allocate(my_file)
         my_file%kdi_file=c_null_ptr
         my_file%unit=ret_val
      endif
      return
      !----
      !  Error
      !----
      6000 NULLIFY(my_file)
      return
   end function FILOPN
   integer function FILCLS(my_file,iactio)
      use, intrinsic :: iso_c_binding
      type(FIL_file), pointer :: my_file
      integer :: iactio
      integer, parameter :: ndummy=4
      type(c_ptr) :: my_kdi_file
      integer,external :: KDICL,KDRCLS
      integer :: ret_val
      !
      itapno=my_file%unit
      my_kdi_file=my_file%kdi_file
      !----
      !  kdi (word addressable) file open
      !----
      if((itapno == 0).and.c_associated(my_kdi_file)) then
         iercod=KDICL(my_kdi_file,iactio)
         ret_val=-1
         if(iercod /= 0) go to 7000
      else if((itapno > 0).and..not.c_associated(my_kdi_file)) then
         ret_val=-2
         iercod=KDRCLS(itapno,iactio)
         if(iercod /= 0) go to 7000
      else
         ret_val=-3
         go to 7000
      endif
      deallocate(my_file)
      FILCLS=0
      return
      7000 FILCLS=ret_val
      return
   end function FILCLS
   integer function FILUNIT(file_pt)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in) :: file_pt
      type(FIL_file), pointer :: my_file
      !
      call c_f_pointer(file_pt,my_file)
      if(c_associated(my_file%kdi_file)) then
         FILUNIT = -1
         return
      endif
      FILUNIT=my_file%unit
      return
   end function FILUNIT
   function FILKDI(file_pt)
      use, intrinsic :: iso_c_binding
      type(c_ptr) FILKDI
      type(c_ptr), intent(in) :: file_pt
      type(FIL_file), pointer :: my_file
      !
      call c_f_pointer(file_pt,my_file)
      if(my_file%unit > 0) then
         FILKDI=c_null_ptr
         return
      endif
      FILKDI=my_file%kdi_file
      return
   end function FILKDI
end module
