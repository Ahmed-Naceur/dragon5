!
!----------------------------------------------------------------------------
!
!Purpose:
! Find a few eigenvalues and eigenvectors for the standard eigenvalue problem
! A*x = lambda*x using the implicit restarted Arnoldi method (IRAM). 
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
! J. Baglama, "Augmented Block Householder Arnoldi Method,"
! Linear Algebra Appl., 429, Issue 10, 2315-2334 (2008).
!
!Parameters: input
! atv     function pointer for the matrix-vector product returning Ab
!         where b is input. The format for atv is "x=atv(b,n,blsz,iter,...)"
! n       order of matrix A.
! blsz    block size of the Arnoldi Hessenberg matrix (blsz=3 recommended)
! K_org   number of desired eigenvalues
! maxit   maximum number of iterations
! tol     tolerance used for convergence (tol=1.0d-6 recommended)
! impx    print parameter: =0: no print; =1: minimum printing.
! iptrk   L_TRACK pointer to the tracking information
! ipsys   L_SYSTEM pointer to system matrices
! ipflux  L_FLUX pointer to the solution
!
!Parameters: output
! iter    actual number of iterations
! V       eigenvector matrix
! D       eigenvalue diagonal matrix
!
!----------------------------------------------------------------------------
!
subroutine ALBEIGS(atv,n,blsz,K_org,maxit,tol,impx,iter,V,D,iptrk,ipsys,ipflux)
  use GANLIB
  implicit complex(kind=8)(a-h,o-z)
  !----
  !  Subroutine arguments
  !----
  interface
    function atv(b,n,blsz,iter,iptrk,ipsys,ipflux) result(x)
      use GANLIB
      integer, intent(in) :: n,blsz,iter
      complex(kind=8), dimension(n,blsz), intent(in) :: b
      complex(kind=8), dimension(n,blsz) :: x
      type(c_ptr) iptrk,ipsys,ipflux
    end function atv
  end interface
  integer, intent(in) :: n,blsz,K_org,maxit,impx
  real(kind=8), intent(in) :: tol
  integer, intent(out) :: iter
  complex(kind=8), dimension(n,K_org), intent(inout) :: V
  complex(kind=8), dimension(K_org,K_org), intent(out) :: D
  type(c_ptr) iptrk,ipsys,ipflux
  !----
  !  Local variables
  !----
  integer :: Hsz_n,Hsz_m,Tsz_m,H_sz1,H_sz2,adjust
  real(kind=8) :: tol2
  integer, parameter :: nbls_init=10 ! Maximum number of Arnoldi vectors.
  integer, parameter :: iunout=6
  real(kind=8), parameter :: eps=epsilon(tol)
  complex(kind=8), dimension(1,1) :: onebyone
  complex(kind=8), dimension(2,2) :: twobytwo
  logical :: eig_conj
  character(len=1) :: conv
  character(len=131) :: hsmg
  !----
  !  Allocatable arrays
  !----
  integer, allocatable, dimension(:) :: vadjust,iset
  real(kind=8), allocatable, dimension(:) :: residuals,tau
  real(kind=8), allocatable, dimension(:,:) :: work2d_r,w,QR,T_Schur
  complex(kind=8), allocatable, dimension(:) :: V_sign,work1d
  complex(kind=8), allocatable, dimension(:,:) :: T_blk,H_R,DD,Q,R,work2d
  complex(kind=8), pointer, dimension(:) :: DH_eig
  complex(kind=8), pointer, dimension(:,:) :: T,T_old,H,H_old,VC,VC_old,H_eig,H_eigv

  ! Resize Krylov subspace if blsz*nbls (i.e. number of Arnoldi vectors)
  ! is larger than n (i.e. the size of the matrix A).
  if(impx > 1) write(iunout,'(/23H ALBEIGS: IRAM solution)')
  nbls = nbls_init
  if(blsz*nbls >= n) then
    nbls = floor(real(n)/real(blsz))
    write(6,'(26h ALBEIGS: Changing nbls to,i5)') nbls
  endif

  ! Increase the number of desired values to help increase convergence.
  ! Set K_org+adjust(1) to be next multiple of blsz > K_org.
  allocate(vadjust(blsz))
  do i=1,blsz
    vadjust(i)=mod(K_org+i-1,blsz)
  enddo
  adjust=findlc(vadjust,0)
  deallocate(vadjust)
  K = K_org + adjust
  allocate(VC(n,blsz), T(blsz,blsz))

  ! Check for input errors in the structure array.
  if(K <= 0)    call XABORT('ALBEIGS: K must be a positive value.')
  if(K > n)     call XABORT('ALBEIGS: K is too large.')
  if(blsz <= 0) call XABORT('ALBEIGS: blsz must be a positive value.')
  if(blsz > K_org)  call XABORT('ALBEIGS: blsz <= K_org expected.')

  ! Automatically adjust Krylov subspace to accommodate larger values of K.
  if(blsz*(nbls-1) - blsz - K - 1 < 0) then
    nbls = ceiling((K+1)/blsz+2.1)
    write(6,'(26h ALBEIGS: Changing nbls to,i5)') nbls
  endif
  if(blsz*nbls >= n) then
    nbls = floor(n/blsz-0.1)
    write(6,'(26h ALBEIGS: Changing nbls to,i5)') nbls
  endif
  if(blsz*(nbls-1) - blsz - K - 1 < 0)  call XABORT('ALBEIGS: K is too large.')
  VC(:n,:blsz) = V(:n,:blsz) ! set initial eigenvector estimate
  
  tol2 = tol
  if(tol < eps) tol2 = eps ! Set tolerance to machine precision if tol < eps.

  allocate(tau(n),residuals(K_org))
  ! allocate adjustable T matrix in the WY representation of the Householder
  ! products.
  m = blsz     ! Current number of columns in the matrix VC.
  nullify(DH_eig, H_eig, H_eigv)
  allocate(H(blsz,0))
  iter = 0     ! Main loop iteration count.
  do
    iter = iter+1
    if(iter > maxit) then
      call XABORT('ALBEIGS: maximum number of IRAM iterations reached')
    endif
    allocate(T_blk(blsz,blsz)); T_blk(:blsz,:blsz)=0.0d0;
    ! Compute the block Householder Arnoldi decomposition.
    if(blsz*(nbls+1) > n) nbls = floor(real(n)/real(blsz))-1

    ! Begin of main iteration loop for the Augmented Block Householder Arnoldi
    ! decomposition.
    do
      if(m > blsz*(nbls+1)) exit
      allocate(work2d_r(n-m+blsz,blsz))
      work2d_r(:n-m+blsz,:blsz) = real(VC(m-blsz+1:n,m-blsz+1:m))
      call ALST2F(n-m+blsz,n-m+blsz,blsz,work2d_r,tau)
      VC(m-blsz+1:n,m-blsz+1:m) = work2d_r(:n-m+blsz,:blsz)
      deallocate(work2d_r)
      if(m > blsz) then
        H_old => H
        allocate(H(m,m-blsz)); H(:m,:m-blsz) = 0.0d0;
        H(:m-blsz,:m-2*blsz) = H_old(:,:)
        deallocate(H_old)
        H(:m-blsz,m-2*blsz+1:m-blsz) = VC(:m-blsz,m-blsz+1:m)
        do i=m-blsz+1,m
          H(i,i-blsz:m-blsz) = VC(i,i:m)
        enddo
      endif
      m2 = 2*m
      if(m2 > n) m2=n
      do j=1,blsz
        VC(1:m-blsz+j,m-blsz+j) = 0.0d0
      enddo
      do i=m-blsz+1,m
        VC(i,i)=1.0d0
      enddo
      onebyone=matmul(tcg(VC(:,m-blsz+1:m-blsz+1)),VC(:,m-blsz+1:m-blsz+1))
      T_blk(1,1) = -2.0d0/onebyone(1,1)
      do i=2,blsz
        T_blk(:i,i:i) = tcg(matmul(tcg(VC(:,m-blsz+i:m-blsz+i)),VC(:,m-blsz+1:m-blsz+i)))
        T_blk(i,i) = -2.0d0/T_blk(i,i)
        T_blk(:i-1,i) = T_blk(i,i)*matmul(T_blk(:i-1,:i-1),T_blk(:i-1,i))
      enddo

      ! Matrix T expansion
      if(m == blsz) then
        T(:m,:m) = T_blk(:m,:m)
      else
        T_old => T
        allocate(T(m,m))
        T(:m-blsz,:m-blsz) = T_old(:m-blsz,:m-blsz); T(m-blsz+1:m,:m-blsz) = 0.0d0;
        T(:m-blsz,m-blsz+1:m) = matmul(matmul(T_old,tcg(matmul(tcg(VC(:,m-blsz+1:m)),VC(:,1:m-blsz)))),T_blk)
        T(m-blsz+1:m,m-blsz+1:m) = T_blk(:blsz,:blsz)
        deallocate(T_old)
      endif

      ! Resize and reactualize VC
      if(m <= blsz*nbls) then
        VC_old => VC
        allocate(VC(n,m+blsz))
        do j=1,m
          VC(:n,j) = VC_old(:n,j)
        enddo
        deallocate(VC_old)
        VC(:,m+1:m+blsz) = matmul(VC(:,:m),matmul(T,tcg(VC(m-blsz+1:m,:m))))
        do i=1,blsz
          VC(i+m-blsz,i+m) = VC(i+m-blsz,i+m) + 1.0d0
        enddo
        VC(:,m+1:m+blsz) = atv(VC(:,m+1:m+blsz),n,blsz,iter,iptrk,ipsys,ipflux)
        allocate(work2d(n,blsz))
        work2d(:,:) = matmul(VC(:,:m),tcg(matmul(matmul(tcg(VC(:,m+1:m+blsz)),VC(:,:m)),T)))
        VC(:,m+1:m+blsz) = VC(:,m+1:m+blsz) + work2d(:,:)
        deallocate(work2d)
      endif
      m = m + blsz
    enddo
    deallocate(T_blk)

    ! Determine the size of the block Hessenberg matrix H. Possible truncation may occur
    ! if an invariant subspace has been found.
    Hsz_n = size(H,1); Hsz_m = size(H,2);
    
    ! Compute the eigenvalue decomposition of the block Hessenberg H(:Hsz_m,:).
    if(associated(H_eigv)) deallocate(H_eigv,H_eig,DH_eig)
    allocate(H_eigv(Hsz_m,Hsz_m), H_eig(Hsz_m,Hsz_m), DH_eig(Hsz_m))
    allocate(work2d_r(Hsz_m, Hsz_m))
    work2d_r(:Hsz_m,:Hsz_m)=real(H(:Hsz_m,:Hsz_m))
    call ALHQR(Hsz_m, Hsz_m, work2d_r, 200, jter, H_eigv, H_eig)
    deallocate(work2d_r)
    do i=1,Hsz_m
      DH_eig(i) = H_eig(i,i)
    enddo

    ! Check the accuracy of the computation of the eigenvalues of the
    ! Hessenberg matrix. This is used to monitor balancing.
    conv = 'F';  ! Boolean to determine if all desired eigenpairs have converged.

    ! Sort the eigenvalue and eigenvector arrays.
    allocate(iset(Hsz_m),work1d(Hsz_m),work2d(Hsz_m,Hsz_m))
    call ALINDX(Hsz_m, DH_eig(:), iset)
    do i=1,Hsz_m
      work1d(i) = DH_eig(iset(i))
      work2d(:Hsz_m,i) = H_eigv(:Hsz_m,iset(i))
    enddo
    DH_eig(:Hsz_m) = work1d(:Hsz_m); H_eigv(:Hsz_m,:Hsz_m) = work2d(:Hsz_m,:Hsz_m)
    deallocate(work2d,work1d,iset)

    ! Compute the residuals for the K_org Ritz values.
    residuals(:K_org) = sqrt(sum(abs(matmul(H(Hsz_n-blsz+1:Hsz_n,Hsz_m-blsz+1:Hsz_m), &
             H_eigv(Hsz_m-blsz+1:Hsz_m,:K_org)))**2, 1))
    if(impx > 1) write(iunout,200) iter,residuals(:K_org)
             
    ! Check for convergence.
    conv = 'T'
    do i=1,K_org
      if(residuals(i) >= tol2*abs(DH_eig(i))) conv = 'F'
    enddo

    ! Adjust K to include more vectors as the number of vectors converge.
    K =  K_org + adjust
    do i=1,K_org
      if(residuals(i) < eps*abs(DH_eig(i))) K = K+1
    enddo
    if(K > Hsz_m - 2*blsz-1) K = Hsz_m - 2*blsz-1

    ! Determine if K splits a conjugate pair. If so replace K with K + 1.
    if(aimag(DH_eig(K)) /= 0.0d0) then
      eig_conj = .true.
      if(K < Hsz_m) then
        if(abs(aimag(DH_eig(K)) + aimag(DH_eig(K+1))) < sqrt(eps)) then
          K = K + 1
          eig_conj = .false.
        endif
      endif
      if(K > 1 .and. eig_conj) then
        if(abs(aimag(DH_eig(K)) + aimag(DH_eig(K-1))) < sqrt(eps)) then
          eig_conj = .false.
        endif
      endif
      if(eig_conj) then
        write(hsmg,'(9h ALBEIGS:,i5,25h-th conjugate pair split.)') K
        call XABORT(hsmg)
      endif
    endif

    ! If all desired Ritz values converged then exit main loop.
    if(conv == 'T')  exit

    ! Compute the QR factorization of H_eigv(:,:K).
    allocate(iset(Hsz_m))
    nset=0
    do i=1,Hsz_m
      if(abs(aimag(DH_eig(i))) > 1.0d-10) then
        nset=nset+1
        iset(nset)=i
      endif
    enddo
    allocate(Q(Hsz_m,K))
    if(nset == 0) then
      Q(:,:) = H_eigv(:,:K)
    else
      ! Convert the complex eigenvectors of the eigenvalue decomposition of H
      ! to real vectors and convert the complex diagonal matrix to block diagonal.
      allocate(work2d(Hsz_m,Hsz_m)); work2d(:Hsz_m,:Hsz_m) = 0.0d0;
      do i=1,Hsz_m
        work2d(i,i) = 1.0d0
      enddo
      twobytwo(1,1) = cmplx(1.0d0, 0.0d0, kind=8); twobytwo(2,1) = cmplx(0.0d0, 1.0d0, kind=8);
      twobytwo(1,2) = cmplx(1.0d0, 0.0d0, kind=8); twobytwo(2,2) = -cmplx(0.0d0, 1.0d0, kind=8);
      do i=1,Hsz_m
        ii=findlc(iset(:nset),i)
        if(mod(ii-1,2)+1.eq.1) then
          if(conjg(DH_eig(i)) /= DH_eig(i+1)) call XABORT('ALBEIGS: invalid diagonal')
          work2d(i:i+1,i:i+1) = twobytwo;
        endif
      enddo
      call ALINVC(Hsz_m,work2d,Hsz_m,ier)
      if(ier /= 0) call XABORT('ALBEIGS: singular matrix(1)')
      Q(:,:) = matmul(H_eigv(:,:Hsz_m),work2d(:Hsz_m,:K))
      deallocate(work2d)
    endif
    deallocate(iset)
    allocate(work2d_r(Hsz_m,K))
    do i=1,Hsz_m
      work2d_r(i,:K) = real(Q(i,:K))
    enddo
    deallocate(Q)
    call ALST2F(Hsz_m,Hsz_m,K,work2d_r(:,:K),tau)
    allocate(QR(Hsz_m,K)); QR(:Hsz_m,:K) = 0.0d0;
    do i=1,K
      QR(i,i) = 1.0d0
    enddo
    do j = K,1,-1
      allocate(w(Hsz_m-j+1,1))
      w(:,:) = reshape((/1.0d0, work2d_r(j+1:Hsz_m,j)/), (/Hsz_m-j+1, 1/))
      QR(j:Hsz_m,:) = QR(j:Hsz_m,:)+tau(j)*matmul(w,matmul(transpose(w),QR(j:Hsz_m,:)))
      deallocate(w)
    enddo
    deallocate(work2d_r)

    ! The Schur matrix for H.
    allocate(T_Schur(K,K))
    T_Schur = matmul(matmul(transpose(QR),real(H(:Hsz_m,:))),QR)
    do i=3,K
      T_Schur(i,:i-2) = 0.0d0
    enddo

    ! Compute the starting vectors and the residual vectors from the Householder
    ! WY form. The starting vectors will be the first K Schur vectors and the
    ! residual vectors are stored as the last blsz vectors in the Householder WY form.
    Tsz_m = size(T,1)
    VC(:,Hsz_n-blsz+1:Hsz_n)= matmul(VC(:,:Tsz_m),matmul(T,tcg(VC(Tsz_m-blsz+1:Tsz_m,:Tsz_m))))
    do i=Tsz_m-blsz+1,Tsz_m-blsz+blsz
      VC(i,i) = VC(i,i) + 1.0d0
    enddo
    allocate(work2d(n,K))
    do j=1,K
      work2d(:,j) = matmul(VC(:,:Hsz_m),matmul(T(:Hsz_m,:Hsz_m),matmul(tcg(VC(:Hsz_m,:Hsz_m)),QR(:,j))))
    enddo
    do j=1,K
      VC(:,j) = work2d(:,j)
      VC(:Hsz_m,j) = QR(:Hsz_m,j) + VC(:Hsz_m,j)
    enddo
    deallocate(work2d)

    ! Set the size of the large matrix VC and move the residual vectors.
    m = K + 2*blsz; VC(:,K+1:K+blsz) = VC(:,Hsz_n-blsz+1:Hsz_n);

    ! Set the new starting vector(s) to be the desired vectors VC(:,:K) with the
    ! residual vectors VC(:,Hsz_n-blsz+1:Hsz_n). Place all vectors in the compact
    ! WY form of the Householder product. Compute the next set of vectors by
    ! computing A*VC(:,Hsz_n-blsz+1:Hsz_n) and store this in VC(:,Hsz_n+1:Hsz_n+blsz).
    m2=m-blsz
    allocate(R(m2,m2), V_sign(m2), DD(m2,m2))
    R(:m2,:m2) = 0.0d0; V_sign(:m2) = 1.0d0; DD(:m2,:m2) = 0.0d0;
    deallocate(T)
    allocate(T(m2,m2))
    T = VC(:m2,:m2)
    VC(:,m2+1:m) = atv(VC(:,m2-blsz+1:m2),n,blsz,iter,iptrk,ipsys,ipflux)
    do i =1,m2
      V_sign(i) = VC(i,i)/abs(VC(i,i))
      if(VC(i,i) == 0.0d0) V_sign(i)=1.0d0
      R(i,i) = -V_sign(i)
      Vdot = 1.0d0 + V_sign(i)*VC(i,i)  ! Dot product of Householder vectors.
      VC(i,i) = VC(i,i) + V_sign(i)     ! Reflection to the ith axis.
      DD(i,i) = 1.0d0/VC(i,i)           ! Used for scaling. Note: VC(i,i) >= 1.
      VC(:m2,i+1:m2) = VC(:m2,i+1:m2) - (V_sign(i)/Vdot)*matmul(VC(:m2,i:i),VC(i:i,i+1:m2))
    enddo
    VC(:m2,:m2) = matmul(VC(:m2,:m2),DD)
    deallocate(DD, V_sign)
    VC(:,m2+1:m) = matmul(VC(:,m2+1:m),R(m2-blsz+1:m2,m2-blsz+1:m2))
    T = matmul(T,R)
    do i=1,m2
      VC(i,i+1:m2) = 0.0d0
      T(i,i) = T(i,i) - 1.0d0
    enddo
    allocate(H_R(m2,m2))
    H_R(:m2,:m2) = tcg(VC(:m2,:m2))
    call ALINVC(m2,H_R,m2,ier)
    if(ier /= 0) call XABORT('ALBEIGS: singular matrix(2)')
    T(:m2,:m2) = matmul(T, H_R)
    H_R(:m2,:m2) = VC(:m2,:m2)
    call ALINVC(m2,H_R,m2,ier)
    if(ier /= 0) call XABORT('ALBEIGS: singular matrix(3)')
    T(:m2,:m2) = matmul(H_R, T)
    do i=2,m2
      T(i,:i-1) = 0.0d0
    enddo
    H_R(:m2,:m2) = matmul(T,tcg(VC(:m2,:m2)))
    call ALINVC(m2,H_R,m2,ier)
    if(ier /= 0) call XABORT('ALBEIGS: singular matrix(4)')
    allocate(work2d(n-m2,m2))
    do j=1,m2
      work2d(:n-m2,j) = matmul(VC(m2+1:n,:m2),matmul(R(:m2,:m2),H_R(:m2,j)))
    enddo
    do j=1,m2
      VC(m2+1:n,j) = work2d(:n-m2,j)
    enddo
    deallocate(work2d, H_R)
    VC(:,m2+1:m) = VC(:,m2+1:m) + matmul(VC(:,:m2),tcg(matmul(matmul(tcg(VC(:,m2+1:m)),VC(:,:m2)),T)))

    ! Compute the first K columns and K+blsz rows of the matrix H, used in augmenting.
    allocate(H_R(blsz,blsz))
    H_R(:blsz,:blsz) = H(Hsz_n-blsz+1:Hsz_n,Hsz_m-blsz+1:Hsz_m)
    deallocate(H)
    allocate(H(K+blsz,K)); H(:K+blsz,:K)=0.0d0;
    H(:K,:K) = matmul(R(:K,:K),matmul(T_Schur(:K,:K),R(:K,:K)))
    H(K+1:K+blsz,:K) = matmul(R(K+1:K+blsz,K+1:K+blsz),matmul(H_R, &
                           matmul(QR(Hsz_m-(blsz-1):Hsz_m,:K),R(:K,:K))))
    deallocate(T_Schur, H_R, R, QR)
  enddo
  deallocate(residuals,tau)

  ! Truncated eigenvalue and eigenvector arrays to include only desired eigenpairs.
  Tsz_m = size(T,1); H_sz1 = size(H_eigv,1); H_sz2 = size(H_eigv,2);
  do j=1,H_sz2
    VC(:,j) = matmul(VC(:,:Tsz_m),matmul(T,matmul(tcg(VC(:H_sz1,:Tsz_m)),H_eigv(:H_sz1,j))))
  enddo
  VC(:H_sz1,:H_sz2) = H_eigv + VC(:H_sz1,:H_sz2)

  ! Set the first K_org eigensolutions
  D(:K_org,:K_org)=0.0d0
  do i=1,K_org
    V(:,i) = VC(:,i)
    D(i,i) = DH_eig(i)
  enddo
  deallocate(H_eigv,H_eig,DH_eig,VC)
  return
  !
  200 format(25h ALBEIGS: outer iteration,i4,12h  residuals=,1p,10e12.4/(41x,10e12.4))

  contains
  function tcg(ac) result(bc)
    ! function emulating complex conjugate transpose in Matlab
    complex(kind=8), dimension(:,:), intent(in) :: ac
    complex(kind=8), dimension(size(ac,2),size(ac,1)) :: bc
    bc(:,:)=transpose(conjg(ac(:,:)))
  end function tcg
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
end subroutine ALBEIGS
