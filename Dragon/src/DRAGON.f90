!
!-----------------------------------------------------------------------
!
!Purpose:
! DRAGON main program.
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
program DRAGON
   use GANLIB
   implicit none
   integer, parameter :: iout=6
   character(len=131) :: hsmg
!----
! local storage
!----
   integer       :: iprint,ier
!----
! gan-2000 external functions
!----
   integer, external :: KERNEL
   interface
      integer(c_int) function dramod(cmodul, nentry, hentry, ientry, jentry, &
                     kentry, hparam_c) bind(c)
         use, intrinsic :: iso_c_binding
         character(kind=c_char), dimension(*) :: cmodul
         integer(c_int), value :: nentry 
         character(kind=c_char), dimension(13,*) :: hentry
         integer(c_int), dimension(nentry) :: ientry, jentry
         type(c_ptr), dimension(nentry) :: kentry
         character(kind=c_char), dimension(73,*) :: hparam_c
       end function dramod
    end interface
!----
!  variables for DRAGON version
!----
      integer           :: imvers
      character(len=64) :: date
      character(len=48) :: rev
      character(len=6), parameter :: namsbr='dragon'
!----
!  version information recovered from cvs
!----
   imvers=5
   call KDRVER(rev,date)
   write(iout,6000) namsbr,imvers,rev,date
   write(iout,6010) namsbr
!----
!  execute the cle-2000 driver
!----
   iprint=0
   ier=KERNEL(dramod,iprint)
   if( ier /= 0 )then
      write(hsmg,'(27hDRAGON: kernel error (code=,I5,2h).)') ier
      call XABORT(hsmg)
   endif
!----
!  all modules processed
!----
   write(iout,6030) namsbr,imvers,rev
   stop
!----
!  formats
!----
 6000 FORMAT( &
       ' DDDDDD   RRRRRR     AAA     GGGGG    OOOOO   NN     NN'/ &
       ' DDDDDDD  RRRRRRR   AAAAA   GGGGGGG  OOOOOOO  NNN    NN'/ &
       ' DD   DD  RR   RR  AA   AA  GG       OO   OO  NNNN   NN'/ &
       ' DD   DD  RRRRRR   AA   AA  GG  GGG  OO   OO  NN NN  NN'/ &
       ' DD   DD  RRRR     AAAAAAA  GG  GGG  OO   OO  NN  NN NN'/ &
       ' DD   DD  RR RR    AAAAAAA  GG   GG  OO   OO  NN   NNNN'/ &
       ' DDDDDDD  RR  RR   AA   AA  GGGGGGG  OOOOOOO  NN    NNN'/ &
       ' DDDDDD   RR   RR  AA   AA   GGGGG    OOOOO   NN     NN'// &
       ' VERSION ',A6,I2,2X,A,4X,A/ &
       ' GROUPE D''ANALYSE NUCLEAIRE'/ &
       ' ECOLE POLYTECHNIQUE DE MONTREAL'///)
 6010 FORMAT( &
       ' COPYRIGHT NOTICE FOR THIS VERSION OF ',A6,':'/ &
       ' --------------------------------------------'/ &
       ' Copyright (C) 2002 Ecole Polytechnique de Montreal '/ &
       ' This library is free software; you can redistribute it' / &
       ' and/or modify it under the terms of the GNU Lesser General' / &
       ' Public License as published by the Free Software Foundation;' / &
       ' either version 2.1 of the License, or (at your option) any' / &
       ' later version '////)
 6030 FORMAT(/1x,'normal end of execution for ',a6,i2,2x,a/ &
              1x,'check for warning in listing'/ &
              1x,'before assuming your run was successful')
end program DRAGON
