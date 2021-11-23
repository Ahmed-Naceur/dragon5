*DECK FLDONE
      FUNCTION FLDONE(X,B,N,IPTRK,IPSYS,IPFLUX) RESULT(Y)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Computation of a single X+M*(B-A*X) iteration in TRIVAC.
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
* X       initial flux.
* B       fixed source.
* N       number of unknowns in the flux.
* IPTRK   L_TRACK pointer to the tracking information.
* IPSYS   L_SYSTEM pointer to system matrices.
* IPFLUX  L_FLUX pointer to the solution.
*
*Parameters: output
* Y       flux at the next iteration.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER, INTENT(IN) :: N
      REAL(KIND=8), DIMENSION(N), INTENT(IN) :: X, B
      REAL(KIND=8), DIMENSION(N) :: Y
      TYPE(C_PTR) IPTRK,IPSYS,IPFLUX
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE)
      CHARACTER*12 TEXT12
      REAL, DIMENSION(:), ALLOCATABLE :: WORK1,WORK2,GAR
*
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NLF=ISTATE(30)
      CALL LCMGET(IPSYS,'STATE-VECTOR',ISTATE)
      NGRP=ISTATE(1)
      LL4=ISTATE(2)
      ITY=ISTATE(4)
      NBMIX=ISTATE(7)
      NAN=ISTATE(8)
      IF(ITY.EQ.13) LL4=LL4*NLF/2 ! SPN cases
      CALL LCMGET(IPFLUX,'STATE-VECTOR',ISTATE)
      IGR=ISTATE(39)
      IF(LL4.NE.N) CALL XABORT('FLDONE: INCONSISTENT UNKNOWNS.')
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK1(LL4),WORK2(LL4))
      WRITE(TEXT12,'(1HA,2I3.3)') IGR,IGR
      WORK1(:LL4)=REAL(B(:LL4))
      WORK2(:LL4)=REAL(X(:LL4))
      IF(ITY.EQ.2) THEN
*       CLASSICAL TREATMENT
        ALLOCATE(GAR(LL4))
        CALL MTLDLM(TEXT12,IPTRK,IPSYS,LL4,ITY,WORK2,GAR)
        GAR(:LL4)=WORK1(:LL4)-GAR(:LL4)
        CALL MTLDLS(TEXT12,IPTRK,IPSYS,LL4,ITY,GAR)
        WORK2(:LL4)=WORK2(:LL4)+GAR(:LL4)
        DEALLOCATE(GAR)
      ELSE IF(ITY.EQ.3) THEN
*       THOMAS-RAVIART/DIFFUSION TRIVAC TRACKING.
        CALL FLDTRS(TEXT12,IPTRK,IPSYS,LL4,WORK1,WORK2,1)
      ELSE IF(ITY.EQ.13) THEN
*       THOMAS-RAVIART/SIMPLIFIED PN TRIVAC TRACKING.
        IF(NAN.EQ.0) CALL XABORT('FLDONE: SPN-ONLY ALGORITHM(2).')
        CALL FLDSPN(TEXT12,IPTRK,IPSYS,LL4,NBMIX,NAN,WORK1,WORK2,1)
      ELSE
        CALL XABORT('FLDONE: INVALID TYPE.')
      ENDIF
      Y(:LL4)=WORK2(:LL4)
      DEALLOCATE(WORK2,WORK1)
      RETURN
      END FUNCTION FLDONE
