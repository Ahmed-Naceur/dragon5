*DECK FLDARN
      SUBROUTINE FLDARN (FLDATV,IPTRK,IPSYS,IPFLUX,LL4,NUN,NGRP,LMOD,
     1 IBLSZ,ADJ,IMPX,EPSOUT,MAXOUT,EVECT,FKEFFV)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of a multigroup eigenvalue system for the calculation of the
* LMOD first orthogonal harmonics of the diffusion or SPN equation.
* Use the implicit restarted Arnoldi method (IRAM).
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* FLDATV  function pointer for the multiplication of A^(-1)B times the
*         harmonic flux
* IPTRK   L_TRACK pointer to the BIVAC tracking information.
* IPSYS   L_SYSTEM pointer to system matrices.
* IPFLUX  L_FLUX pointer to the solution.
* LL4     order of the system matrices.
* NUN     number of unknowns in each energy group.
* NGRP    number of energy groups.
* LMOD    number of orthogonal harmonics to compute.
* IBLSZ   block size of the Arnoldi Hessenberg matrix.
* ADJ     adjoint calculation flag.
* IMPX    print parameter: =0: no print; =1: minimum printing.
* EPSOUT  convergence criteria for the flux.
* MAXOUT  maximum number of outer iterations.
* EVECT   initial estimate of the unknown vector.
*
*Parameters: output
* EVECT   converged unknown vector.
* FKEFFV  effective multiplication factor of each harmonic.
*
*Reference:
* J. BAGLAMA, "Augmented Block Householder Arnoldi Method,"
* Linear Algebra Appl., 429, Issue 10, 2315-2334 (2008).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS,IPFLUX
      INTEGER LL4,NUN,NGRP,LMOD,IBLSZ,IMPX,MAXOUT
      LOGICAL ADJ
      REAL EPSOUT
      COMPLEX EVECT(NUN,NGRP,LMOD),FKEFFV(LMOD)
*----
*  LOCAL VARIABLES
*----
      INTERFACE
        FUNCTION FLDATV(F,N,IBLSZ,ITER,IPTRK,IPSYS,IPFLUX) RESULT(X)
          USE GANLIB
          INTEGER, INTENT(IN) :: N,IBLSZ,ITER
          REAL(KIND=8), DIMENSION(N,IBLSZ), INTENT(IN) :: F
          REAL(KIND=8), DIMENSION(N,IBLSZ) :: X
          TYPE(C_PTR) IPTRK,IPSYS,IPFLUX
        END FUNCTION FLDATV
      END INTERFACE
      REAL TIME(2)
      REAL(KIND=8) DEPSOUT
      CHARACTER(LEN=8) TEXT8
      TYPE(C_PTR) JPFLUX,KPFLUX,MPFLUX
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR
      COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: V, D
*----
*  SCRATCH STORAGE ALLOCATION
*----
      N=LL4*NGRP
      ALLOCATE(V(N,LMOD),D(LMOD,LMOD),GAR(NUN))
*----
*  SET TIMER
*----
*     TIME(1) : CPU TIME FOR THE SOLUTION OF LINEAR SYSTEMS.
*     TIME(2) : CPU TIME FOR BILINEAR PRODUCT EVALUATIONS.
      TIME(1)=0.0
      TIME(2)=0.0
      CALL LCMPUT(IPFLUX,'CPU-TIME',2,2,TIME)
*----
*  FLUX INITIALIZATION
*----
      DO IMOD=1,LMOD
        V(:N,IMOD)=1.0D0
        V(1:MIN(IBLSZ,IMOD)-1,IMOD)=0.0D0
      ENDDO
      CALL LCMLEN(IPFLUX,'MODE',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
        DO IMOD=1,LMOD
          JPFLUX=LCMGID(IPFLUX,'MODE')
          CALL LCMLEL(JPFLUX,IMOD,ILONG,ITYLCM)
          IF(ILONG.EQ.0) CYCLE
          KPFLUX=LCMGIL(JPFLUX,IMOD)
          IF(ADJ) THEN
            CALL LCMLEN(KPFLUX,'AFLUX',LENA,ITYLCM)
            IF(LENA.EQ.0) CYCLE
            MPFLUX=LCMGID(KPFLUX,'AFLUX')
            DO IGR=1,NGRP
              IF(ITYLCM.EQ.2) THEN
                CALL LCMGDL(MPFLUX,IGR,GAR)
                EVECT(:NUN,IGR,IMOD)=GAR(:NUN)
              ELSE IF(ITYLCM.EQ.6) THEN
                CALL LCMGDL(MPFLUX,IGR,EVECT(1,IGR,IMOD))
              ENDIF
            ENDDO
          ELSE
            CALL LCMLEN(KPFLUX,'FLUX',LEND,ITYLCM)
            IF(LEND.EQ.0) CYCLE
            MPFLUX=LCMGID(KPFLUX,'FLUX')
            DO IGR=1,NGRP
              IF(ITYLCM.EQ.2) THEN
                CALL LCMGDL(MPFLUX,IGR,GAR)
                EVECT(:NUN,IGR,IMOD)=GAR(:NUN)
              ELSE IF(ITYLCM.EQ.6) THEN
                CALL LCMGDL(MPFLUX,IGR,EVECT(1,IGR,IMOD))
              ENDIF
            ENDDO
          ENDIF
          DO IGR=1,NGRP
            DO IUN=1,LL4
              IOF=(IGR-1)*LL4+IUN
              V(IOF,IMOD)=EVECT(IUN,IGR,IMOD)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
*----
*  CALL IRAM SOLVER
*----
      DEPSOUT=EPSOUT
      CALL ALBEIGS(FLDATV,N,IBLSZ,LMOD,MAXOUT,DEPSOUT,IMPX,ITER,V,D,
     1 IPTRK,IPSYS,IPFLUX)
      DO IMOD=1,LMOD
        FKEFFV(IMOD)=CMPLX(D(IMOD,IMOD),KIND=4)
        DO IGR=1,NGRP
          DO IUN=1,LL4
            IOF=(IGR-1)*LL4+IUN
            EVECT(IUN,IGR,IMOD)=CMPLX(V(IOF,IMOD),KIND=4)
          ENDDO
        ENDDO
      ENDDO
*----
*  PRINTOUTS
*----
      IF(IMPX.GE.1) THEN
        CALL LCMGET(IPFLUX,'CPU-TIME',TIME)
        WRITE (6,650) ITER,TIME(1),TIME(2),TIME(1)+TIME(2)
        WRITE (6,670) (FKEFFV(IMOD),IMOD=1,LMOD)
      ENDIF
      IF(IMPX.GE.3) THEN
        TEXT8=' DIRECT'
        IF(ADJ) TEXT8=' ADJOINT'
        DO IMOD=1,LMOD
          WRITE (6,'(/A8,13H HARMONIC NB.,I3/)') TEXT8,IMOD
          DO IGR=1,NGRP
            WRITE (6,680) IGR,(REAL(EVECT(I,IGR,IMOD)),I=1,LL4)
          ENDDO
        ENDDO
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAR,D,V)
      RETURN
*
  650 FORMAT(/31H FLDARN: CONVERGENCE OF IRAM IN,I5,11H ITERATIONS/
     1 9X,54HCPU TIME USED TO SOLVE THE TRIANGULAR LINEAR SYSTEMS =,
     2 F10.3/23X,34HTO COMPUTE THE BILINEAR PRODUCTS =,F10.3,20X,
     3 16HTOTAL CPU TIME =,F10.3)
  670 FORMAT(//21H FLDARN: EIGENVALUES:/(5X,1P,E17.10,3H + ,E17.10,1Hi))
  680 FORMAT(43H FLDARN: EIGENVECTOR CORRESPONDING TO GROUP,I4//
     1 (5X,1P,8E14.5))
      END
