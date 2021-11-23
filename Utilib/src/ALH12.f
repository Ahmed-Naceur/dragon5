*DECK ALH12
      SUBROUTINE ALH12(MODE,LPIVOT,L1,M,U,UP,C,ICE,ICV,NCV)
*
*-----------------------------------------------------------------------
*
*Purpose:
* construction and/or application of a single Householder
* transformation (Q = I + U*(U**T)/B).
*
*Author(s): C. L. Lawson and R. J. Hanson
*
*Parameters: input
* MODE    algorithm flag (=1/2: Selects algorithm H1 to construct and
*         apply a Householder transformation / Algorithm H2 to apply a
*         previously constructed transformation).
* LPIVOT  index of the pivot element.
* L1,M    if L1.le.M, the transformation will be constructed to
*         zero elements indexed from L1 through M. If L1.gt.M,
*         the subroutine does an identity transformation.
* U       pivot vector (if MODE=1). At exit, contains quantities
*         defining the vector U of the Householder transformation.
*         On entry with MODE=2, U should contain information
*         previously computed with MODE=1.
* C       matrix which will be regarded as a set of vectors to which the
*         Householder transformation is to be applied.
* ICE     storage increment between elements of vectors in C.
* ICV     storage increment between vectors in C.
* NCV     number of vectors in C to be transformed. if NCV.le.0, no
*         operations will be done on C.
*
*Parameters: output
* U       quantitie defining the vector U of the Householder
*         transformation. These will not be  modified during the entry
*         with MODE=2.
* UP      quantity defining the vector U of the Householder
*         transformation. On entry with MODE = 2, UP should
*         contain information previously computed with MODE=1. These
*         will not be modified during the entry with MODE=2.
* C       set of transformed vectors.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  Subroutine arguments
*----
      INTEGER MODE, LPIVOT, L1, M, ICE, ICV, NCV
      DOUBLE PRECISION U(M),C(M),UP
*----
*  Local variables
*----
      INTEGER I, I2, I3, I4, INCR, J
      DOUBLE PRECISION B, CL, CLINV, SM
*
      IF(0.GE.LPIVOT .OR. LPIVOT.GE.L1 .OR. L1.GT.M) RETURN
      CL = ABS(U(LPIVOT))
      IF(MODE.NE.2) THEN
*----
*  Construct the transformation.
*----
        DO J = L1, M
          CL=MAX(ABS(U(J)),CL)
        ENDDO
        IF(CL.LE.0) RETURN
        CLINV = 1.0D0 / CL
        SM = (U(LPIVOT)*CLINV) ** 2 + SUM( (U(L1:M)*CLINV)**2 )
        CL = CL * SQRT(SM)
        IF(U(LPIVOT).GT.0) THEN
          CL = -CL
        ENDIF
        UP = U(LPIVOT) - CL
        U(LPIVOT) = CL
      ELSE
*----
*  Apply the transformation  I+U*(U**T)/B to C.
*----
        IF(CL.LE.0) RETURN
      ENDIF
      IF(NCV.LE.0) RETURN
      B = UP * U(LPIVOT)
*----
*  B must be nonpositive here. If B = 0., return.
*----
      IF(B.LT.0.0) THEN
        B = 1.0D0 / B
        I2 = 1 - ICV + ICE * (LPIVOT-1)
        INCR = ICE * (L1-LPIVOT)
        DO J = 1, NCV
          I2 = I2 + ICV
          I3 = I2 + INCR
          I4 = I3
          SM = C(I2) * UP
          DO I = L1, M
            SM = SM + C(I3) * U(I)
            I3 = I3 + ICE
          ENDDO
          IF(SM.NE.0) THEN
            SM = SM * B
            C(I2) = C(I2) + SM * UP
            DO I = L1, M
              C(I4) = C(I4) + SM * U(I)
              I4 = I4 + ICE
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END
