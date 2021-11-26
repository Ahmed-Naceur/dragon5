!
!-----------------------------------------------------------------------
!
!Purpose:
! find the eigenvalues and corresponding eigenvectors of equation
! (a-eval)*evect=0 using the power method with the shifted Hessenberg
! QR algorithm.
!
!Copyright:
! Copyright (C) 2020 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s): A. Hebert
!
!Reference:
! G. E. Robles, "Implementing the QR algorithm for efficiently
! computing matrix eigenvalues and eigenvectors," Final Degree
! Dissertation in Mathematics, Universidad del Pais Vasco, Spain
! (2017).
!
!Parameters: input
! ndim    dimensioned column length of A.
! n       order of matrix A.
! A       input matrix.
! maxiter maximum number of iterations.
!
!Parameters: output
! iter    actual number of iterations.
! V       eigenvector matrix.
! D       eigenvalue diagonal matrix.
!
!-----------------------------------------------------------------------
!
subroutine ALHQR(ndim,n,A,maxiter,iter,V,D)
  implicit none
  !----
  !  Subroutine arguments
  !----
  integer, intent(in) :: ndim,n,maxiter
  integer, intent(out) :: iter
  real(kind=8), dimension(ndim,n), intent(in) :: A
  complex(kind=8), dimension(n,n), intent(out) :: V,D
  !----
  !  Local variables
  !----
  integer :: i,j,k,i1,i2,nset,ier,ii
  complex(kind=8) :: kappa,p,qq,r,mu,csum,denom,m(2,2)
  real(kind=8) :: sgn,tau,nu,normH,sum,AA(2,2),BB(2,2),CC(2,2)
  real(kind=8), parameter :: eps=epsilon(A)
  !----
  !  Allocatable arrays
  !----
  integer, allocatable, dimension(:) :: iset
  real(kind=8), allocatable, dimension(:) :: c
  complex(kind=8), allocatable, dimension(:) :: s,t,work1d
  real(kind=8), allocatable, dimension(:,:) :: VR
  complex(kind=8), allocatable, dimension(:,:) :: H,Q,work2d
  !----
  ! Perform Householder transformation to upper Hessenberg form
  !----
  allocate(H(n,n), Q(n,n), VR(n,n-2))
  H(:n,:n)=A(:n,:n)
  do k = 1,(n-2)
    VR(k+1:n,k) = real(H(k+1:n,k))
    sgn = sign(1.0d0,VR(k+1,k))
    VR(k+1,k) = VR(k+1,k) + sgn * sqdotv(VR(k+1:n,k))
    VR(k+1:n,k) = VR(k+1:n,k) / sqdotv(VR(k+1:n,k))
    H(k+1:n,k:n) = H(k+1:n,k:n) - 2.0d0 * matmul(reshape(VR(k+1:n,k),(/n-k, 1/)), &
                   reshape(matmul(VR(k+1:n,k),H(k+1:n,k:n)),(/1, n-k+1/)))
    H(:,k+1:n) = H(:,k+1:n) - 2.0d0 * matmul(reshape(matmul(H(:,k+1:n),VR(k+1:n,k)),(/n, 1/)), &
                   reshape(VR(k+1:n,k),(/1, n-k/)))
  enddo
  !----
  ! Construct Q matrix
  !----
  Q(:n,:n) = 0.0D0
  do j=1,n
    Q(j,j)=1.0d0
  enddo
  do j = (n-2),1,-1
    Q(j+1:n,:n) = Q(j+1:n,:n) - 2.0d0 * matmul(reshape(VR(j+1:n,j),(/n-j, 1/)), &
                  reshape(matmul(VR(j+1:n,j),Q((j+1):n,:)),(/1, n/)))
  enddo
  deallocate(VR)
  !----
  ! Perform Schur factorization
  !----
  i2 = n
  allocate(c(n),s(n),t(n))
  c(:n)=0.0d0; s(:n)=0.0d0; t(:n)=0.0d0;
  iter = 0
  do
    iter = iter + 1
    if(iter > maxiter) then
      call xabort('ALHQR: maximum number of iterations exceeded.')
    endif
    ! Check subdiagonal for near zeros, deflating points. Finds deflating rows
    ! on a complex Schur form matrix.
    i1 = i2
    normH = sqdotm(abs(H(:n,:n)))
    do
      if(i1 == 1) exit
      if(abs(H(i1,i1-1)) < eps*normH) then
        H(i1,i1-1) = 0.0d0
        if(i1 == i2) then
          i2 = i1 - 1; i1 = i1 - 1;
        else
          exit
        endif
      else
        i1 = i1 - 1
      endif
    enddo
    !----
    ! End the function if H is upper triangular
    !----
    if(i2 == 1) exit
    ! Compute Wilkinson shift
    kappa = H(i2,i2)
    sum = abs(H(i2-1,i2-1)) + abs(H(i2-1,i2)) + abs(H(i2,i2-1)) + abs(H(i2,i2))
    if(sum /= 0) then
      qq = (H(i2-1,i2)/sum)*(H(i2,i2-1)/sum)
      if(qq /= 0) then
        p = 0.5*((H(i2-1,i2-1)/sum) - (H(i2,i2)/sum))
        r = sqrt(p*p + qq);
        if( (real(p)*real(r) + imag(p)*imag(r)) < 0 ) then
          r = -r
        endif
        kappa = kappa - sum*(qq/(p+r))
      endif
    endif
    ! Apply shift to the element of the diagonal that is left out of the loop
    H(i1,i1) = H(i1,i1) - kappa
    do j = i1,i2-1 ! Loop reducing the matrix to triangular form
      ! Apply Givens rotation so that the subdiagonal is set to zero
      if(H(j+1,j) == 0) then
        c(j) = 1.0d0; s(j) = 0.0d0;
      elseif(H(j,j) == 0) then
        c(j) = 0.0d0; s(j) = 1; H(j,j) = H(j+1,j); H(j+1,j) = 0.0d0;
      else
        mu = H(j,j)/abs(H(j,j))
        tau = abs(real(H(j,j))) + abs(imag(H(j,j))) + abs(real(H(j+1,j))) &
                + abs(imag(H(j+1,j)))
        nu = tau*sqrt(abs(H(j,j)/tau)**2 + abs(H(j+1,j)/tau)**2)
        c(j) = abs(H(j,j))/nu
        s(j) = mu*conjg(H(j+1,j))/nu
        H(j,j) = nu*mu
        H(j+1,j) = 0.0d0
      endif
      ! Apply shift to diagonal
      H(j+1,j+1) = H(j+1,j+1) - kappa
      ! Modify the involved rows using a plane rotation 
      t(j+1:n) = c(j)*H(j,j+1:n) + s(j)*H(j+1,j+1:n)
      H(j+1,j+1:n) = c(j)*H(j+1,j+1:n) - conjg(s(j))*H(j,j+1:n)
      H(j,j+1:n) = t(j+1:n)
    enddo
    do k = i1,i2-1
      ! Loop applying the back multiplication using a plane rotation 
      t(1:k+1) = c(k)*H(1:k+1,k) + conjg(s(k))*H(1:k+1,k+1);
      H(1:k+1,k+1) = c(k)*H(1:k+1,k+1) - s(k)*H(1:k+1,k)
      H(1:k+1,k) = t(1:k+1)
      ! Accumulate transformations using a plane rotation 
      t(1:n) = c(k)*Q(1:n,k) + conjg(s(k))*Q(1:n,k+1)
      Q(1:n,k+1) = c(k)*Q(1:n,k+1) - s(k)*Q(1:n,k)
      Q(1:n,k) = t(1:n)
      H(k,k) = H(k,k) + kappa
    enddo
    H(i2,i2) = H(i2,i2) + kappa
  enddo
  deallocate(t,s,c)
  !----
  ! Construct the orthonormal basis
  !----
  V(:n,:n)=0.0d0
  D(:n,:n)=0.0d0
  do i=1,n
    V(i,i)=1.0d0
    D(i,i)=H(i,i)
  enddo
  do j=2,n
    do i=j-1,1,-1
      denom=H(i,i)-H(j,j)
      if(denom /= 0) then
        csum=0.0d0
        do k=i+1,j
          csum=csum+H(i,k)*V(k,j)
        enddo
        V(i,j)=V(i,j)-csum/denom
      endif
    enddo
  enddo
  V=matmul(Q,V)
  deallocate(Q,H)
  !----
  ! Sort and normalize the eigensolution
  !----
  allocate(iset(n),work1d(n),work2d(n,n))
  do i=1,n
    work1d(i) = D(i,i)
  enddo
  call ALINDX(n, work1d, iset)
  do i=1,n
    work1d(i) = D(iset(i),iset(i))
    work2d(:n,i) = V(:n,iset(i))
  enddo
  do i=1,n
    D(i,i)=work1d(i)
  enddo
  V(:n,:n) = work2d(:n,:n)
  deallocate(work2d,work1d)
  nset=0
  do i=1,n
    if(abs(imag(D(i,i))) > 1.0e-10) then
      nset=nset+1
      iset(nset)=i
    endif
  enddo
  do i=1,n
    ii=findlc(iset(:nset),i)
    if(mod(ii-1,2)+1.eq.1) then
      j=iset(ii+1)
      m=reshape( (/V(i,i), V(j,i), V(i,j), V(j,j)/), (/2, 2/) )
      m(:,1)=m(:,1)/sqdotv(abs(m(1:2,1)))
      m(:,2)=m(:,2)/sqdotv(abs(m(1:2,2)))
      AA=reshape( (/real(m(1,1))+real(m(2,1)), aimag(m(1,1))+aimag(m(2,1)), &
        -aimag(m(1,1))-aimag(m(2,1)), real(m(1,1))+real(m(2,1)) /), (/2, 2/) )
      BB=reshape( (/real(m(1,2))+real(m(2,2)), -aimag(m(1,2))-aimag(m(2,2)), &
        -aimag(m(1,2))-aimag(m(2,2)), -real(m(1,2))-real(m(2,2)) /), (/2, 2/) )
      call ALINVD(2,BB,2,ier)
      if(ier.ne.0) call xabort('ALHQR: singular matrix')
      CC=matmul(BB,AA)
      V(:,i)=V(:,i)*cmplx(CC(1,1),CC(1,2),kind=8)
    elseif (mod(ii-1,2)+1.eq.2) then
      j=iset(ii-1)
      if(abs(D(i,i)-conjg(D(j,j))) > 1.0e-10) then
        call xabort('ALHQR: pathological ordering')
      endif
      D(i,i)=conjg(D(j,j))
    else
      D(i,i)=real(D(i,i))
    endif
    V(:,i)=V(:,i)/sqdotv(abs(V(:,i)))
  enddo
  deallocate(iset)
  return

  contains
  function sqdotv(vec) result(vsum)
    ! function emulating the vectorial norm2 function in Fortran 2008
    real(kind=8), dimension(:), intent(in) :: vec
    real(kind=8) :: vsum
    vsum=sqrt(dot_product(vec(:),vec(:)))
  end function sqdotv
  function sqdotm(mat) result(vsum)
    ! function emulating the matrix norm2 function in Fortran 2008
    real(kind=8), dimension(:,:), intent(in) :: mat
    real(kind=8) :: vsum
    vsum=0.0d0
    do i=1,size(mat,1)
      do j=1,size(mat,2)
        vsum=vsum+mat(i,j)**2
      enddo
    enddo
    vsum=sqrt(vsum)
  end function sqdotm
  function findlc(iset,itest) result(ii)
    ! function emulating the findloc function in Fortran 2008
    integer, dimension(:), intent(in) :: iset
    integer, intent(in) :: itest
    integer :: ii
    ii=0
    do j=1,size(iset)
      if(iset(j) == itest) then
        ii=j
        exit
      endif
    enddo
  end function findlc
end subroutine ALHQR
