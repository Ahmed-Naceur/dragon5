*DECK ALINDX
      SUBROUTINE ALINDX(n,arr,indx)
*
*-----------------------------------------------------------------------
*
*Purpose:
* indexes an array arr(1:n) such that abs(arr(indx(j))) is in descending
* order for j=1:n.
* 
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de MontCOMPLEX(KIND=8)
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* n       size of the array a.
* a       array a.
*
*Parameters: output
* indx    indexing vector.
*
*-----------------------------------------------------------------------
*
*----
*  subroutine arguments
*----
      INTEGER n,indx(n)
      COMPLEX(KIND=8) arr(n)
*----
*  local variables
*----
      INTEGER nstack,i,indxt,ir,itemp,j,jstack,k,l
      COMPLEX(KIND=8) a
      INTEGER, parameter :: M=7
      INTEGER, allocatable, dimension(:) :: istack
*
      nstack=2*ceiling(log(real(n))/log(2.0))
      allocate(istack(nstack))
      do j=1,n
        indx(j)=j
      enddo
      jstack=0
      l=1
      ir=n
    1 if(ir-l.lt.M) then
        do j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do i=j-1,l,-1
            if(absc(arr(indx(i))).ge.absc(a)) goto 2
            indx(i+1)=indx(i)
          enddo
          i=l-1
    2     indx(i+1)=indxt
        enddo
        if(jstack.eq.0)go to 5
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(absc(arr(indx(l))).lt.absc(arr(indx(ir))))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(absc(arr(indx(l+1))).lt.absc(arr(indx(ir)))) then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(absc(arr(indx(l))).lt.absc(arr(indx(l+1)))) then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
    3   i=i+1
        if(absc(arr(indx(i))).gt.absc(a)) goto 3
    4   j=j-1
        if(absc(arr(indx(j))).lt.absc(a)) goto 4
        if(j.ge.i) then
          itemp=indx(i)
          indx(i)=indx(j)
          indx(j)=itemp
          goto 3
        endif
        indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.nstack) call XABORT('ALINDX: nstack too small.')
        if(ir-i+1.ge.j-l) then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
    5 deallocate(istack)
      return
*
      contains
      function absc(val) result(aval)
        ! definition of combined value for complex numbers
        COMPLEX(kind=8), intent(in) :: val
        REAL(kind=8) :: aval
        aval=100.0d0*real(val)+aimag(val)
      end function absc
      END
