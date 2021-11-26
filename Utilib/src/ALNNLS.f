*DECK ALNNLS
      SUBROUTINE ALNNLS(A,M,N,MP,NP,B,X,RNORM,MODE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* implementation of the nonnegative least square algorithm.
*
*Author(s): A. Miller
*
*Parameters: input
* A       input rectangular matrix.
* M,N     first/second mathematical dimension of matrix A.
* MP,NP   first/second physical dimension of matrix A.
* B       source M-vector.
*
*Parameters: output
* A       product matrix, Q*A, where Q IS AN M x M orthogonal matrix
*         generated implicitly by this subroutine.
* B       product Q*B.
* X       solution N-vector.
* RNORM   Euclidean norm of the residual vector.
* MODE    success-failure flag (1: the solution has been computed
*         successfully; 2: the dimensions of the problem are
*         inconsistent, either m.LE.0 or n.LE.0; 3: iteration count
*         exceeded with more than 10*n iterations).
*
*Reference:
*  Chen, Donghui; Plemmons, Robert J. (2009). Nonnegativity constraints
*  in numerical analysis. Symposium on the Birth of Numerical Analysis.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  Subroutine arguments
*----
      INTEGER M, N, MP, NP, MODE
      DOUBLE PRECISION A(MP,NP), B(M), X(N), RNORM
*----
*  Local variables
*----
      INTEGER I, II, IP, ITER, ITMAX, IZ, IZ1, IZ2, IZMAX, J, JJ, JZ,
     >        L, MDA, NPP1, NSETP
      DOUBLE PRECISION DUMMY(1), ALPHA, ASAVE, CC, FACTOR, SM, SS, T,
     >                 TEMP, UNORM, UP, WMAX, ZTEST, XR, YR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INDX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: W, ZZ
      PARAMETER(FACTOR = 0.01D0)
*
      MODE = 1
      IF(M.LE.0 .OR. N.LE.0) THEN
        MODE = 2
        RETURN
      ENDIF
      ITER = 0
      ITMAX = 10*N
      ALLOCATE(INDX(N),W(N), ZZ(M))
*----
*  Initialize the arrays INDX() and X().
*----
      DO I = 1,N
        X(I) = 0.0D0
        INDX(I) = I
      ENDDO
      IZ2 = N
      IZ1 = 1
      NSETP = 0
      NPP1 = 1
*----
*                ******  MAIN LOOP BEGINS HERE  ******
*  Quit if all coefficients are already in the solution or if M cols of
*  A have been triangularized.
*----
   30 IF(IZ1.GT.IZ2 .OR. NSETP.GE.M) GO TO 350
*----
*  Compute components of the dual (negative gradient) vector W.
*----
      DO IZ = IZ1,IZ2
        J = INDX(IZ)
        W(J) = DOT_PRODUCT(A(NPP1:M,J), B(NPP1:M))
      ENDDO
*----
*  Find largest positive W(J).
*----
   60 IZMAX = 0
      WMAX = 0.0D0
      DO IZ = IZ1,IZ2
        J = INDX(IZ)
        IF(W(J).GT.WMAX) THEN
          WMAX = W(J)
          IZMAX = IZ
        ENDIF
      ENDDO
*----
*  If WMAX.le.0. go to termination. This indicates satisfaction of the
*  Kuhn-Tucker conditions.
*----
      IF(WMAX.LE.0.0D0) GO TO 350
      IZ = IZMAX
      J = INDX(IZ)
*----
*  The sign of W(J) is ok for J to be moved to set P. Begin the
*  transformation and check new diagonal element to avoid near linear
*  dependence.
*----
      ASAVE = A(NPP1,J)
      CALL ALH12(1, NPP1, NPP1+1, M, A(:,J), UP, DUMMY, 1, 1, 0)
      UNORM = 0.0D0
      IF(NSETP.NE.0) UNORM = SUM( A(1:NSETP,J)**2 )
      UNORM = SQRT(UNORM)
      IF(UNORM + ABS(A(NPP1,J))*FACTOR - UNORM .GT. 0.0D0) THEN
*----
*  Col J is sufficiently independent. Copy B into ZZ, update ZZ
*  and solve for ZTEST ( = proposed new value for X(J) ).
*----
        ZZ(1:M) = B(1:M)
        CALL ALH12(2, NPP1, NPP1+1, M, A(:,J), UP, ZZ, 1, 1, 1)
        ZTEST = ZZ(NPP1)/A(NPP1,J)
*----
*  See if ZTEST is positive.
*----
        IF(ZTEST.GT.0.0D0) GO TO 140
      ENDIF
*----
*  Reject J as a candidate to be moved from set Z to set P.  Restore
*  A(NPP1,J), set W(J) = 0., and loop back to test dual coeffs again.
*----
      A(NPP1,J) = ASAVE
      W(J) = 0.0D0
      GO TO 60
*----
*  The index  J = INDX(IZ)  has been selected to be moved from
*  set Z to set P. Update B, update indices, apply Householder
*  transformations to cols in new set Z, sero subdiagonal elts in
*  col J, set W(J) = 0.
*----
  140 B(1:M) = ZZ(1:M)
      INDX(IZ) = INDX(IZ1)
      INDX(IZ1) = J
      IZ1 = IZ1+1
      NSETP = NPP1
      NPP1 = NPP1+1
      MDA = SIZE(A,1)
      JJ = 0
      IF(IZ1.LE.IZ2) THEN
        DO JZ = IZ1,IZ2
          JJ = INDX(JZ)
          CALL ALH12(2, NSETP, NPP1, M, A(:,J), UP, A(:,JJ), 1, MDA, 1)
        ENDDO
      ENDIF
      IF(NSETP.NE.M)  A(NPP1:M,J) = 0.0D0
      W(J) = 0.0D0
*----
*  Solve the triangular system. Store the solution temporarily in ZZ().
*----
      DO L = 1, NSETP
        IP = NSETP+1-L
        IF(L .NE. 1) ZZ(1:IP) = ZZ(1:IP) - A(1:IP,JJ)*ZZ(IP+1)
        JJ = INDX(IP)
        ZZ(IP) = ZZ(IP) / A(IP,JJ)
      ENDDO
*----
*             ******  SECONDARY LOOP BEGINS HERE ******
*----
  210 ITER = ITER+1
      IF(ITER.GT.ITMAX) THEN
        MODE = 3
        WRITE (*,'(/A)') ' NNLS QUITTING ON ITERATION COUNT.'
        GO TO 350
      ENDIF
*----
*  See if all new constrained coeffs are feasible; if not compute ALPHA.
*----
      ALPHA = 2.0D0
      DO IP = 1,NSETP
        L = INDX(IP)
        IF(ZZ(IP).LE.0.0D0) THEN
          T = -X(L)/(ZZ(IP)-X(L))
          IF(ALPHA.GT.T) THEN
            ALPHA = T
            JJ = IP
          ENDIF
        ENDIF
      ENDDO
*----
*  If all new constrained coeffs are feasible then ALPHA will still be
*  equal to 2. If so exit from secondary loop to main loop.
*----
      IF(ALPHA == 2.0D0) GO TO 330
*----
*  Otherwise use ALPHA which will be between 0. and 1. to interpolate
*  between the old X and the new ZZ.
*----
      DO IP = 1,NSETP
        L = INDX(IP)
        X(L) = X(L) + ALPHA*(ZZ(IP)-X(L))
      ENDDO
*----
*  Modify A and B and the index arrays to move coefficient I from set
*  P to set Z.
*----
      I = INDX(JJ)
  260 X(I) = 0.0D0
*----
*  Compute.. matrix   (C, S) so that (C, S)(A) = (SQRT(A**2+B**2))
*                     (-S,C)         (-S,C)(B)   (   0          )
*  Compute SIG = SQRT(A**2+B**2)
*     SIG is computed last to allow for the possibility that SIG
*     may be in the same location as A OR B .
*----
      IF(JJ.NE.NSETP) THEN
        JJ = JJ+1
        DO J = JJ,NSETP
          II = INDX(J)
          INDX(J-1) = II
          IF(ABS(A(J-1,II)).GT.ABS(A(J,II))) THEN
            XR = A(J,II) / A(J-1,II)
            YR = SQRT(1.0D0 + XR**2)
            CC = SIGN(1.0D0/YR, A(J-1,II))
            SS = CC * XR
            A(J-1,II) = ABS(A(J-1,II)) * YR
          ELSE IF(A(J,II).NE.0.D0) THEN
            XR = A(J-1,II) / A(J,II)
            YR = SQRT(1.0D0 + XR**2)
            SS = SIGN(1.0D0/YR, A(J,II))
            CC = SS * XR
            A(J-1,II) = ABS(A(J,II)) * YR
          ELSE
            CC = 0.0D0
            SS = 1.0D0
          ENDIF
          A(J,II) = 0.0D0
          DO L = 1,N
            IF(L.NE.II) THEN
*----
*  Apply procedure G2 (CC,SS,A(J-1,L),A(J,L))
*----
              TEMP = A(J-1,L)
              A(J-1,L) = CC*TEMP + SS*A(J,L)
              A(J,L)   = -SS*TEMP + CC*A(J,L)
            ENDIF
          ENDDO
*----
*  Apply procedure G2 (CC,SS,B(J-1),B(J))
*----
          TEMP = B(J-1)
          B(J-1) = CC*TEMP + SS*B(J)
          B(J)   = -SS*TEMP + CC*B(J)
        ENDDO
      ENDIF
      NPP1 = NSETP
      NSETP = NSETP-1
      IZ1 = IZ1-1
      INDX(IZ1) = I
*----
*  See if the remaining coeffs in set P are feasible. They should be
*  because of the way ALPHA was determined. if any are infeasible it is
*  due to round-off error. any that are nonpositive will be set to 0.0d0
*  and moved from set P to set Z.
*----
      DO JJ = 1,NSETP
        I = INDX(JJ)
        IF(X(I).LE.0.0D0) GO TO 260
      ENDDO
*----
*  Copy B( ) into ZZ( ). Then solve again and loop back.
*----
      ZZ(1:M) = B(1:M)
      DO L = 1, NSETP
        IP = NSETP+1-L
        IF(L .NE. 1) ZZ(1:IP) = ZZ(1:IP) - A(1:IP,JJ)*ZZ(IP+1)
        JJ = INDX(IP)
        ZZ(IP) = ZZ(IP) / A(IP,JJ)
      ENDDO
      GO TO 210
*----
*          ******  End of secondary loop  ******
*----
  330 DO IP = 1,NSETP
        I = INDX(IP)
        X(I) = ZZ(IP)
      ENDDO
*----
*  All new coeffs are positive. Loop back to beginning.
*----
      GO TO 30
*----
*  Come to here for termination. Compute the norm of the final residual
*  vector.
*----
  350 SM = 0.0D0
      IF(NPP1.LE.M) THEN
        SM = SUM( B(NPP1:M)**2 )
      ELSE
        W(1:N) = 0.0D0
      ENDIF
      RNORM = SQRT(SM)
      DEALLOCATE(ZZ,W,INDX)
      RETURN
      END
