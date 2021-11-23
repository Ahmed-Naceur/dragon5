program GANMAIN
   use GANLIB
   implicit none
   character(len=131) :: hsmg
!
! local storage
   integer       :: iprint,ier
   integer           :: imvers
   character(len=64) :: date
   character(len=48) :: rev
   integer, parameter :: iout=6
   character(len=6), parameter :: namsbr='ganlib'
!
! gan-2000 external functions
   integer, external :: KERNEL
   interface
      integer(c_int) function ganmod(cmodul, nentry, hentry, ientry, jentry, &
                     kentry, hparam_c) bind(c)
         use, intrinsic :: iso_c_binding
         character(kind=c_char), dimension(*) :: cmodul
         integer(c_int), value :: nentry 
         character(kind=c_char), dimension(13,*) :: hentry
         integer(c_int), dimension(nentry) :: ientry, jentry
         type(c_ptr), dimension(nentry) :: kentry
         character(kind=c_char), dimension(73,*) :: hparam_c
       end function ganmod
    end interface
!----
!  version information recovered from cvs
!----
   imvers=5
   call KDRVER(rev,date)
!----
!  execute the cle-2000 driver
!----
   iprint=0
   ier=KERNEL(ganmod,iprint)
   if( ier /= 0 )then
      write(hsmg,'(28hGANMAIN: kernel error (code=,I5,2h).)') ier
      call XABORT(hsmg)
   endif
   write(iout,6030) namsbr,imvers,rev
   stop
   6030 format(/1x,'normal end of execution for ',a6,i2,2x,a/ &
                1x,'check for warning in listing'/ &
                1x,'before assuming your run was successful')
end program GANMAIN
