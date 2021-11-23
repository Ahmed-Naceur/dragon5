!
!---------------------------  KDROPN  ----------------------------------
!
! 1- programme statistics:
!     name     : KDROPN, KDRCLS
!     use      : allocate and release file units associated to a given
!                file name. word addressable (KDI), sequential (formatted
!                or not) and direct access (DA) files are permitted
!     modified : 91-01-24
!     author   : G. Marleau and A. Hebert
!
! 2- routine parameters:
!
! - function KDROPN(cuname,iactio,iutype,lrda)
!
!   open file and allocate unit number. allocate unit number to file
!   name if unit is already opened, returns unit number.
!
!   input
!     cuname   : filename                                     c*12
!                if cuname=' ', use a default name
!     iactio   : action on file                               i
!                =0 to allocate a new file
!                =1 to access and modify an existing file
!                =2 to access an existing file in
!                   read-only mode
!                =3 unknown
!     iutype   : file type                                    i
!                =1  KDI word addressable file
!                =2  sequential unformatted
!                =3  sequential formatted
!                =4  direct access (DA) unformatted file
!     ldra     : number of words in DA record                 i
!                required for  iutype = 4 only
!   output
!     KDROPN   : unit number/error status                     type(c_ptr)
!                address of file (successful allocation)
!                == NULL  allocation failure
!
!   error codes:
!                = -1  no more unit available
!                = -2  file type requested inconsistent with file type of this file
!                = -3  this file has been already opened
!                = -4  file name is reserved or too long (FT06F00, FT07F00)
!                = -5  illegal file type 1 < iutype < 4
!                = -6  (not used)
!                = -7  error on open of unformatted sequential file
!                = -8  error on open of formatted sequential file
!                = -9  error on open of direct access file
!                =-10  invalid number of word in direct access record
!                      lrda must be > 0
!
!---------------------------  KDROPN  ----------------------------------
!
integer function KDROPN(cuname,iactio,iutype,lrda)
!----
!   subroutine arguments
!----
   character(len=*) :: cuname
   integer :: iactio,iutype,lrda
!----
!   local variables
!----
   integer :: ret_val=0
   integer, parameter :: nbtape=99,nreser=2,ndummy=4
   character :: crdnam*72,cform*11,cstatu*12
   integer :: itapno
   logical :: lfilop
   character(len=8),save,dimension(ndummy) :: cdummy= &
      (/ 'DUMMYKDI','DUMMYSQ ','DUMMYCA ','DUMMYIN ' /)
   character(len=8),save,dimension(nreser) :: creser= &
      (/ 'FT05F001','FT06F001' /)
   character(len=22),save,dimension(ndummy) :: ctype= &
      (/ 'WORD ADDRESSABLE KDI  ','SEQUENTIAL UNFORMATTED', &
         'SEQUENTIAL CHARACTER  ','DIRECT ACCESS DA      ' /)
!----
!  check if iutype is valid
!----
   if((iutype > 4).or.(iutype <= 1)) then
      ret_val=-5
      go to 6000
   endif
!----
!  check if lrda is valid
!----
   if((iutype == 4).and.(lrda < 1)) then
      ret_val=-10
      go to 6000
   endif
!----
!  check if file name is more than 72 characters
!----
   luname= len(cuname)
   if(luname > 72) then
      ret_val=-4
      go to 6000
   endif
!----
!  check if file name not forbidden
!----
   if(luname < 8) go to 120
   do ireser=1,nreser
      if(cuname(:8) == creser(ireser)) then
         ret_val=-4
         go to 6000
      endif
   enddo
!----
!  check for dummy file name/allocate dummy file name if requested
!----
   do idummy=1,ndummy
      if(cuname(:8) == cdummy(idummy)) then
         if(idummy /= iutype) then
            ret_val=-2
            go to 6000
         endif
      endif
   enddo
   120 if(cuname == ' ') then
      crdnam=cdummy(iutype)
   else
      crdnam=cuname
   endif
!----
!  check if file opened/permitted
!----
   inquire(file=crdnam,opened=lfilop)
   if(lfilop) then
      ret_val = -3
      go to 6000
   endif
!----
!  look for never allocated unit location
!----
   do jboucl=nbtape,1,-1
      itapno=jboucl
      inquire(unit=itapno,opened=lfilop)
      if(.not.lfilop) go to 121
   enddo
!----
!  error - no unit number available
!----
   ret_val = -1
   go to 6000
!
   121 if(iutype == 2) then
!----
!  open sequential unformatted file
!----
      ret_val=-7
      cform='UNFORMATTED'
   else if(iutype == 3) then
!----
!  open sequential formatted file
!----
      ret_val=-8
      cform='FORMATTED'
   else if(iutype == 4) then
!----
!  open DA file
!----
      ret_val=-9
      cform='UNFORMATTED'
   endif
   if(iactio == 0) then
      cstatu='NEW'
   else if(iactio == 1) then
      cstatu='OLD'
   else if(iactio == 2) then
      cstatu='OLD'
   else
      cstatu='UNKNOWN'
   endif
   if((iutype == 4).and.(iactio == 2)) then
      idummy=0
      inquire(iolength=lrecl) (idummy,i=1,lrda)
      open(unit=itapno,file=crdnam,err=7000,iostat=iercod,form=cform, &
      access='DIRECT',recl=lrecl,status=cstatu,action='READ')
   else if(iutype == 4) then
      idummy=0
      inquire(iolength=lrecl) (idummy,i=1,lrda)
      open(unit=itapno,file=crdnam,err=7000,iostat=iercod,form=cform, &
      access='DIRECT',recl=lrecl,status=cstatu)
   else if(((iutype == 2).or.(iutype == 3)).and.(iactio == 0)) then
      open(unit=itapno,file=crdnam,err=7000,iostat=iercod,form=cform, &
      access='SEQUENTIAL',status=cstatu)
   else if(((iutype == 2).or.(iutype == 3)).and.(iactio == 1)) then
      open(unit=itapno,file=crdnam,err=7000,iostat=iercod,form=cform, &
      access='SEQUENTIAL',position='APPEND',status=cstatu)
   else if(((iutype == 2).or.(iutype == 3)).and.(iactio == 2)) then
      open(unit=itapno,file=crdnam,err=7000,iostat=iercod,form=cform, &
      access='SEQUENTIAL',status=cstatu,action='READ')
      rewind(itapno)
   endif
   KDROPN=itapno
   return
   6000 write(6,8000) crdnam,ctype(iutype),ret_val
   KDROPN=ret_val
   return
   7000 write(6,9000) crdnam,ctype(iutype),ret_val,iercod
   KDROPN=ret_val
   return
!----
!  error format
!----
   8000 format('1',5x,'ERROR IN OPENING OF FILE IN KDROPN'/ &
               6x,'FILE NAME = ',a/6x,'FILE TYPE = ',a22/ &
               6x,'ERROR CODE= ',i7)
   9000 format('1',5x,'ERROR IN OPENING OF FILE IN KDROPN'/ &
               6x,'FILE NAME = ',a/6x,'FILE TYPE = ',a22/ &
               6x,'ERROR CODE= ',i7,'   IERCOD = ',i7)
end function KDROPN
!
!-------------------------- KDRCLS -------------------------------------
!
! - function KDRCLS(itapno,iactio)
!
!   close file and release unit number.
!
!   input
!     itapno   : unit number                                  i
!                = 0    close all units
!                > 0    unit to close
!     iactio   : action on file                               i
!                = 1    to keep the file;
!                = 2    to delete the file.
!   output
!     KDRCLS   : error status                                 i
!                =  0  unit closed
!                = -2  file not opened
!                = -3  this file has been opened by routines other than KDROPN
!                = -4  file unit is reserved (5,6)
!                = -5  illegal unit number
!                = -6  (not used)
!                = -7  error on close of unformatted sequential file
!                = -8  error on close of formatted sequential file
!                = -9  error on close of DA file
!                =-10  invalid close action (iactio=1,2 permitted only)
!                =-11  type of file not supported
!
!-------------------------- KDRCLS -------------------------------------
!
integer function KDRCLS(itapno,iactio)
!----
!   subroutine arguments
!----
   integer :: itapno,iactio
!----
!   local variables
!----
   integer, parameter :: ndummy=4
   character(len=10) :: acc
   character(len=11) :: frm
   character(len=22),save,dimension(ndummy) :: ctype= &
    (/ '                      ','SEQUENTIAL UNFORMATTED', 'SEQUENTIAL CHARACTER  ', &
       'DIRECT ACCESS DA      ' /)
!
   integer, parameter :: nbtape=99
   integer :: ret_val=0,itapet=0
   logical :: lfilop,lnmd
   character (len=72) :: cuname
!----
!  invalid unit number
!----
   if((itapno <= 0).or.(itapno > nbtape)) then
      ret_val=-5
      go to 7000
   endif
   inquire(unit=itapno,opened=lfilop,named=lnmd)
   if((.not.lfilop).or.(.not.lnmd)) then
      ret_val=-2
      go to 7000
   endif
!----
!  close the file
!----
   inquire(unit=itapno,access=acc,form=frm)
   if((acc == 'SEQUENTIAL').and.(frm == 'UNFORMATTED')) then
      itapet=2
      ret_val=-7
   else if((acc == 'SEQUENTIAL').and.(frm == 'FORMATTED')) then
      itapet=3
      ret_val=-8
   else if((acc == 'DIRECT').and.(frm == 'UNFORMATTED')) then
      itapet=4
      ret_val=-9
   else
      ret_val=-11
      go to 7000
   endif
   if(iactio == 1) then
       close(itapno,iostat=iercod,status='KEEP',err=7000)
   else if(iactio == 2) then
       close(itapno,iostat=iercod,status='DELETE',err=7000)
   else
      ret_val=-10
      go to 7000
   endif
   KDRCLS=0
   return
   7000 inquire(unit=itapno,name=cuname)
   write(6,9000) itapno,cuname,ctype(itapet),ret_val,iercod
   KDRCLS=ret_val
   return
!----
!  error format
!----
   9000 format('1',5x,'ERROR IN CLOSE OF FILE IN KDROPN'/ &
        6x,'UNIT NB.  = ',i10/6x,'FILE NAME = ',a7/6x,'FILE TYPE = ',a22/ &
        6x,'ERROR CODE= ',i7,'   IERCOD = ',i7)
end function KDRCLS
