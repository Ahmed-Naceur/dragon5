!
!-----------------------------------------------------------------------
!
!Purpose:
! GMRES(m) linear equation solver.
!
!Copyright:
! Copyright (c) 2019 Ecole Polytechnique de Montreal
! this library is free software; you can redistribute it and/or
! modify it under the terms of the gnu lesser general public
! license as published by the free software foundation; either
! version 2.1 of the license, or (at your option) any later version
!
!Author(s): A. Hebert
!
!Reference:
! based on a Matlab script by C. T. Kelley, July 10, 1994.
!
!Parameters: input
! B       fixed source
! atv     function pointer for the matrix-vector product returning
!         X+M*(B-A*X) where X is the unknown and B is the source.
!         The format for atv is "Y=atv(X,B,n,...)"
! n       order of matrix A
! ertol   iteration convergence criterion
! nstart  restarts the GMRES method every nstart iterations
! maxit   maximum number of GMRES iterations.
! impx    print parameter: =0: no print; =1: minimum printing.
! iptrk   L_TRACK pointer to the tracking information
! ipsys   L_SYSTEM pointer to system matrices
! ipflux  L_FLUX pointer to the solution
!
!Parameters: input/output
! X       initial estimate / solution of the linear system.
!
!Parameters: output
! iter    actual number of iterations
!
!----------------------------------------------------------------------------
!
subroutine FLDMRA(B,atv,n,ertol,nstart,maxit,impx,iptrk,ipsys,ipflux,X,iter)
  use GANLIB
  implicit real(kind=8) (a-h,o-z)
  !----
  !  subroutine arguments
  !----
  real(kind=8), dimension(n), intent(in) :: B
  integer, intent(in) :: nstart,maxit,impx
  real(kind=8), intent(in) :: ertol
  interface
    function atv(X,B,n,iptrk,ipsys,ipflux) result(Y)
      use GANLIB
      integer, intent(in) :: n
      real(kind=8), dimension(n), intent(in) :: X, B
      real(kind=8), dimension(n) :: Y
      type(c_ptr) iptrk,ipsys,ipflux
    end function atv
  end interface
  real(kind=8), dimension(n), intent(inout) :: X
  integer, intent(out) :: iter
  type(c_ptr) iptrk,ipsys,ipflux
  !----
  !  local variables
  !----
  integer, parameter :: iunout=6
  !----
  !  allocatable arays
  !----
  real(kind=8), allocatable, dimension(:) :: r,qq,g,c,s
  real(kind=8), allocatable, dimension(:,:) :: v,h
  !----
  !  scratch storage allocation
  !----
  allocate(v(n,nstart+1),g(nstart+1),h(nstart+1,nstart+1), &
  c(nstart+1),s(nstart+1))
  !----
  !  global GMRES(m) iteration.
  !----
  allocate(r(n),qq(n))
  eps1=ertol*sqrt(dot_product(B(:n),B(:n)))
  rho=1.0d10
  iter=0
  do while((rho > eps1).and.(iter < maxit))
    r(:)=atv(X,B,n,iptrk,ipsys,ipflux)-X(:)
    rho=sqrt(dot_product(r(:n),r(:n)))
    !----
    !  test for termination on entry
    !----
    if(rho < eps1) then
       deallocate(qq,r)
       go to 100
    endif
    !
    g(:nstart+1)=0.0d0
    h(:nstart,:nstart)=0.0d0
    v(:n,:nstart+1)=0.0d0
    c(:nstart+1)=0.0d0
    s(:nstart+1)=0.0d0
    g(1)=rho
    v(:n,1)=r(:n)/rho
    !----
    !  gmres(1) iteration
    !----
    k=0
    do while((rho > eps1).and.(k < nstart).and.(iter < maxit))
      k=k+1
      iter=iter+1
      if(impx > 2) write(iunout,200) iter,rho,eps1
      qq(:n)=0.0d0
      r(:)=atv(v(:,k),qq,n,iptrk,ipsys,ipflux)
      v(:n,k+1)=v(:n,k)-r(:n)
      !----
      !  modified Gram-Schmidt
      !----
      do j=1,k
        hr=dot_product(v(:n,j),v(:n,k+1))
        h(j,k)=hr
        v(:n,k+1)=v(:n,k+1)-hr*v(:n,j)
      enddo
      h(k+1,k)=sqrt(dot_product(v(:n,k+1),v(:n,k+1)))
      !----
      !  reorthogonalize
      !----
      do j=1,k
        hr=dot_product(v(:n,j),v(:n,k+1))
        h(j,k)=h(j,k)+hr
        v(:n,k+1)=v(:n,k+1)-hr*v(:n,j)
      enddo
      h(k+1,k)=sqrt(dot_product(v(:n,k+1),v(:n,k+1)))
      !----
      !  watch out for happy breakdown 
      !----
      if(h(k+1,k) /= 0.0) then
        v(:n,k+1)=v(:n,k+1)/h(k+1,k)
      endif
      !----
      !  form and store the information for the new Givens rotation
      !----
      do i=1,k-1
        w1=c(i)*h(i,k)-s(i)*h(i+1,k)
        w2=s(i)*h(i,k)+c(i)*h(i+1,k)
        h(i,k)=w1
        h(i+1,k)=w2
      enddo
      znu=sqrt(h(k,k)**2+h(k+1,k)**2)
      if(znu /= 0.0) then
        c(k)=h(k,k)/znu
        s(k)=-h(k+1,k)/znu
        h(k,k)=c(k)*h(k,k)-s(k)*h(k+1,k)
        h(k+1,k)=0.0d0
        w1=c(k)*g(k)-s(k)*g(k+1)
        w2=s(k)*g(k)+c(k)*g(k+1)
        g(k)=w1
        g(k+1)=w2
      endif
      !----
      !  update the residual norm
      !----
      rho=abs(g(k+1))
    enddo
    !----
    !  at this point either k > nstart or rho < eps1.
    !  it's time to compute x and cycle.
    !----
    h(:k,k+1)=g(:k)
    call ALSBD(k,1,h,ier,nstart+1)
    if(ier /= 0) call XABORT('FLDMRA: singular matrix.')
    do i=1,n
      X(i)=X(i)+dot_product(v(i,:k),h(:k,k+1))
    enddo
  enddo
  deallocate(qq,r)
  !----
  !  scratch storage deallocation
  !----
  100 deallocate(s,c,h,g,v)
  return
  !
  200 format(24h FLDMRA: outer iteration,i4,10h  L2 norm=,1p,e11.4, &
  6h eps1=,e11.4)
end subroutine FLDMRA
