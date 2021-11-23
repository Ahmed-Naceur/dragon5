!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for the NDAS C API.
!
!Copyright:
! Copyright (C) 2012 Atomic Energy of Canada Limited and Ecole
! Polytechnique de Montreal.
!
!Author(s): A. Hebert
!
!-----------------------------------------------------------------------
!
module FSDF
   private
   public :: XSDOPN, XSDNAM, XSDBLD, XSDISO, XSDTHE, XSDRES, XSDTAB, XSDCL
   interface XSDBLD
      module procedure XSDBLD_I1, XSDBLD_R1
   end interface
   interface XSDISO
      module procedure XSDISO_I1, XSDISO_R1
   end interface
contains
subroutine XSDOPN(namfil, ierr)
   ! open the NDAS file
   use, intrinsic :: iso_c_binding
   use LCMAUX
   character(len=*) :: namfil
   integer ierr
   character(kind=c_char), dimension(73) :: name73
   interface 
      subroutine xsdopn_c(namp, ierr) bind(c)
         use, intrinsic :: iso_c_binding
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: ierr
      end subroutine xsdopn_c
   end interface
   call STRCUT(name73, namfil)
   call xsdopn_c(name73, ierr)
end subroutine XSDOPN
!
subroutine XSDNAM(iset, numericId, isonam, ierr)
   ! recover an isotope name from NDAS file
   use LCMAUX
   use, intrinsic :: iso_c_binding
   integer iset,numericId, ierr
   character(len=*) :: isonam
   character(kind=c_char), dimension(73) :: name73
   interface 
      subroutine xsdnam_c(iset, numericId, namp, ierr) bind(c)
         use, intrinsic :: iso_c_binding
         character(kind=c_char), dimension(*) :: namp
         integer(c_int) :: iset,numericId,ierr
      end subroutine xsdnam_c
   end interface
   call xsdnam_c(iset, numericId, name73, ierr)
   call STRFIL(isonam, name73)
end subroutine XSDNAM
!
subroutine XSDBLD_I1(item, where, ierr)
   ! recover a header or integer record from NDAS file
   use, intrinsic :: iso_c_binding
   integer item,ierr
   integer, dimension(*) :: where
   interface 
      subroutine xsdbld_c(item, where, ierr) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) :: item, ierr
         integer(c_int), dimension(*) :: where
      end subroutine xsdbld_c
   end interface
   call xsdbld_c(item, where, ierr)
end subroutine XSDBLD_I1
!
subroutine XSDBLD_R1(item, where, ierr)
   ! recover a header or real record from NDAS file
   use, intrinsic :: iso_c_binding
   integer item,ierr
   real, dimension(*) :: where
   interface 
      subroutine xsdbld_c(item, where, ierr) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) :: item, ierr
         real(c_float), dimension(*) :: where
      end subroutine xsdbld_c
   end interface
   call xsdbld_c(item, where, ierr)
end subroutine XSDBLD_R1
!
subroutine XSDISO_I1(groupRange, item, nuclideIndex, where, ierr)
   ! recover an integer header for an isotope
   use, intrinsic :: iso_c_binding
   integer groupRange, item, nuclideIndex, ierr
   integer, dimension(*) :: where
   interface 
      subroutine xsdiso_c(groupRange, item, nuclideIndex, where, ierr) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) groupRange, item, nuclideIndex, ierr
         integer(c_int), dimension(*) :: where
      end subroutine xsdiso_c
   end interface
   call xsdiso_c(groupRange, item, nuclideIndex, where, ierr)
end subroutine XSDISO_I1
!
subroutine XSDISO_R1(groupRange, item, nuclideIndex, where, ierr)
   ! recover a real header for an isotope
   use, intrinsic :: iso_c_binding
   integer groupRange, item, nuclideIndex, ierr
   real, dimension(*) :: where
   interface 
      subroutine xsdiso_c(groupRange, item, nuclideIndex, where, ierr) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) groupRange, item, nuclideIndex, ierr
         real(c_float), dimension(*) :: where
      end subroutine xsdiso_c
   end interface
   call xsdiso_c(groupRange, item, nuclideIndex, where, ierr)
end subroutine XSDISO_R1
!
subroutine XSDTHE(groupRange, item, nuclideIndex, index, where, ierr)
   ! recover a cross-section array
   use, intrinsic :: iso_c_binding
   integer groupRange, item, nuclideIndex, index, ierr
   real, dimension(*) :: where
   interface 
      subroutine xsdthe_c(groupRange, item, nuclideIndex, index, where, ierr) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) groupRange, item, nuclideIndex, index, ierr
         real(c_float), dimension(*) :: where
      end subroutine xsdthe_c
   end interface
   call xsdthe_c(groupRange, item, nuclideIndex, index, where, ierr)
end subroutine XSDTHE
!
subroutine XSDRES(nuclideIndex, where, ierr)
   ! recover a resonance information array
   use, intrinsic :: iso_c_binding
   integer nuclideIndex, ierr
   integer, dimension(*) :: where
   interface 
      subroutine xsdres_c(nuclideIndex, where, ierr) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) nuclideIndex, ierr
         integer(c_int), dimension(*) :: where
      end subroutine xsdres_c
   end interface
   call xsdres_c(nuclideIndex, where, ierr)
end subroutine XSDRES
!
subroutine XSDTAB(item, nuclideIndex, resGroup, where, ierr)
   ! recover a resonance cross-section array
   use, intrinsic :: iso_c_binding
   integer item, nuclideIndex, resGroup, ierr
   real, dimension(*) :: where
   interface 
      subroutine xsdtab_c(item, nuclideIndex, resGroup, where, ierr) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int) item, nuclideIndex, resGroup, ierr
         real(c_float), dimension(*) :: where
      end subroutine xsdtab_c
   end interface
   call xsdtab_c(item, nuclideIndex, resGroup, where, ierr)
end subroutine XSDTAB
!
subroutine XSDCL()
   ! close the NDAS file
   interface 
      subroutine xsdcl_c() bind(c)
      end subroutine xsdcl_c
   end interface
   call xsdcl_c()
end subroutine XSDCL
end module FSDF
