*DECK FLDADI
      SUBROUTINE FLDADI (NAMP,IPTRK,IPSYS,LL4,ITY,F1,NADI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform NADI inner iterations with the ADI preconditionning method.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NAMP    name of the ADI-splitted matrix.
* IPTRK   L_TRACK pointer to the tracking information.
* IPSYS   L_SYSTEM pointer to system matrices.
* LL4     order of the matrix.
* ITY     type of coefficient matrix (2: classical Trivac;
*         3: Thomas-Raviart; 13: SPN/Thomas-Raviart).
* F1      source term of the linear system.
* NADI    number of inner ADI iterations.
*
*Parameters: output
* F1      solution of the linear system after NADI iterations.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS
      CHARACTER NAMP*12
      INTEGER LL4,ITY,NADI
      REAL F1(LL4)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ITP(NSTATE)
      REAL, DIMENSION(:), ALLOCATABLE :: S1,GAR
*
      ALLOCATE(S1(LL4))
      S1(:LL4)=F1(:LL4) ! SOURCE TERM
      F1(:LL4)=0.0
      IF(ITY.EQ.2) THEN
*       CLASSICAL TREATMENT
        ALLOCATE(GAR(LL4))
        DO IADI=1,NADI
          IF(IADI.EQ.1) THEN
            GAR(:LL4)=S1(:LL4)
          ELSE
            CALL MTLDLM(NAMP,IPTRK,IPSYS,LL4,ITY,F1,GAR)
            GAR(:LL4)=S1(:LL4)-GAR(:LL4)
          ENDIF
          CALL MTLDLS(NAMP,IPTRK,IPSYS,LL4,ITY,GAR)
          F1(:LL4)=F1(:LL4)+GAR(:LL4)
        ENDDO
        DEALLOCATE(GAR)
      ELSE IF(ITY.EQ.3) THEN
*       THOMAS-RAVIART/DIFFUSION TRIVAC TRACKING.
        CALL FLDTRS(NAMP,IPTRK,IPSYS,LL4,S1,F1,NADI)
      ELSE IF(ITY.EQ.13) THEN
*       THOMAS-RAVIART/SIMPLIFIED PN TRIVAC TRACKING.
        CALL LCMGET(IPSYS,'STATE-VECTOR',ITP)
        NBMIX=ITP(7)
        NAN=ITP(8)
        IF(NAN.EQ.0) CALL XABORT('FLDADI: SPN-ONLY ALGORITHM.')
        CALL FLDSPN(NAMP,IPTRK,IPSYS,LL4,NBMIX,NAN,S1,F1,NADI)
      ENDIF
      DEALLOCATE(S1)
      RETURN
      END
