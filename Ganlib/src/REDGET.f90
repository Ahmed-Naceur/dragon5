!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for CLE-2000. REDGET and REDPUT support
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
subroutine REDGET(ityp, nitma, flott, text, dflot)
   ! read a value from input deck
   use, intrinsic :: iso_c_binding
   use LCMAUX
   integer :: ityp, nitma
   real :: flott
   character(len=*) :: text
   double precision :: dflot
   character(kind=c_char), dimension(73) :: text_c
   interface
      subroutine redget_c (ityp, nitma, flott, text_c, dflot) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) :: ityp, nitma
         real(c_float) :: flott
         character(kind=c_char), dimension(*) :: text_c
         real(c_double) :: dflot
      end subroutine redget_c
   end interface
   call redget_c(ityp, nitma, flott, text_c, dflot)
   if(ityp == 3) call STRFIL(text, text_c)
end subroutine REDGET
!
subroutine REDPUT(ityp, nitma, flott, text, dflot)
   ! write a value into the input deck
   use, intrinsic :: iso_c_binding
   use LCMAUX
   integer :: ityp, nitma
   real :: flott
   character(len=*) :: text
   double precision :: dflot
   character(kind=c_char), dimension(73) :: text_c
   interface
      subroutine redput_c (ityp, nitma, flott, text_c, dflot) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) :: ityp, nitma
         real(c_float) :: flott
         character(kind=c_char), dimension(*) :: text_c
         real(c_double) :: dflot
      end subroutine redput_c
   end interface
   if(ityp == 3) call STRCUT(text_c, text)
   call redput_c(ityp, nitma, flott, text_c, dflot)
end subroutine REDPUT
!
subroutine REDOPN(iinp1, iout1, nrec)
   ! read a value from input deck
   use, intrinsic :: iso_c_binding
   use LCMAUX
   type(c_ptr) :: iinp1, file
   integer :: iout1, nrec
   character(len=72) :: filename
   character(kind=c_char), dimension(73) :: filename_c
   interface
      subroutine redopn_c (iinp1, file, filename_c, nrec) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: iinp1, file
         character(kind=c_char), dimension(*) :: filename_c
         integer(c_int), value :: nrec
      end subroutine redopn_c
   end interface
   interface
      function fopen (filename_c, mode) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) fopen
         character(kind=c_char), dimension(*) :: filename_c, mode
      end function fopen
   end interface
   interface
      function stdfil_c (s) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) stdfil_c
         character(kind=c_char) :: s
      end function stdfil_c
   end interface
   if(iout1 == 0) then
      file=c_null_ptr
      filename_c=c_null_char
   else if(iout1 == 6) then
      file=stdfil_c("stdout"//c_null_char)
      filename_c=c_null_char
   else
      inquire(iout1,name=filename)
      close(iout1,status='keep')
      call STRCUT(filename_c, filename)
      file=fopen(filename_c, "w"//c_null_char)
      if(.not.c_associated(file)) call XABORT('REDOPN: UNABLE TO OPEN FILE '//filename(:44))
   endif
   call redopn_c(iinp1, file, filename_c, nrec)
end subroutine REDOPN
!
subroutine REDCLS(iinp1, iout1, nrec)
   ! read a value from input deck
   use, intrinsic :: iso_c_binding
   use LCMAUX
   type(c_ptr) :: iinp1, file
   integer :: iout1, nrec, ier
   character(len=72) :: filename
   character(kind=c_char), dimension(73) :: filename_c
   interface
      subroutine redcls_c (iinp1, file, filename_c, nrec) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iinp1, file
         character(kind=c_char), dimension(*) :: filename_c
         integer(c_int) :: nrec
      end subroutine redcls_c
   end interface
   interface
      function fclose (file) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) fclose
         type(c_ptr), value :: file
      end function fclose
   end interface
   interface
      function stdfil_c (s) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) stdfil_c
         character(kind=c_char) :: s
      end function stdfil_c
   end interface
   call redcls_c(iinp1, file, filename_c, nrec)
   if(c_associated(file,c_null_ptr)) then
      iout1=0
   else if(c_associated(file,stdfil_c("stdout"//c_null_char))) then
      iout1=6
   else
      call STRFIL(filename, filename_c)
      ier=fclose(file)
      if(ier /= 0) call XABORT('REDOPN: UNABLE TO CLOSE FILE '//filename(:44))
      iout1=KDROPN(filename,1,3,0)
   endif
end subroutine REDCLS
