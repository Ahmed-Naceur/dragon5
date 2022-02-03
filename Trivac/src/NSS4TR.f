*DECK NSS4TR
      SUBROUTINE NSS4TR(NEL,NMIX,MAT,XX,IQFR,QFR,DIFF,SIGR,FD,A11)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of leakage system matrices for the coarse mesh finite
* difference method.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NEL     number of nodes
* NMIX    number of mixtures
* MAT     node mixtures
* XX      node widths
* IQFR    boundary conditions
* QFR     albedo functions
* DIFF    diffusion coefficients
* SIGR    macroscopic removal cross sections
* FD      discontinuity factors
*
*Parameters: output
* A11     assembly matrix.
*
*-----------------------------------------------------------------------
*
      INTEGER NEL,NMIX,MAT(NEL),IQFR(6,NEL)
      REAL XX(NEL),QFR(6,NEL),DIFF(NMIX),SIGR(NMIX),FD(NMIX,2),
     1 A11(3*NEL,3*NEL)
*
      A11(:3*NEL,:3*NEL)=0.0
      ! WEIGHT RESIDUAL EQUATIONS:
      NUM1=0
      DO KEL=1,NEL
        IBM=MAT(KEL)
        DX2=XX(KEL)**2
        SIGG=SIGR(IBM)
        DIDD=DIFF(IBM)
        A11(NUM1+1,NUM1+1)=SIGG
        A11(NUM1+1,NUM1+3)=-2.0*DIDD/DX2
        A11(NUM1+2,NUM1+2)=SIGG/12.0
        A11(NUM1+3,NUM1+3)=SIGG/180.0
        NUM1=NUM1+3
      ENDDO
      ! continuity relations:
      NUM1=0
      DO KEL=1,NEL-1
        IBM=MAT(KEL)
        IBMP=MAT(KEL+1)
        DIDD=DIFF(IBM)
        DIDDP=DIFF(IBMP)
        NUM2=NUM1+3
        ! flux continuity:
        FDP=FD(IBM,2)
        FDM=FD(IBMP,1)
        A11(NUM1+2,NUM1+1)=FDP
        A11(NUM1+2,NUM1+2)=FDP/2.0
        A11(NUM1+2,NUM1+3)=FDP/6.0
        A11(NUM1+2,NUM2+1)=-FDM
        A11(NUM1+2,NUM2+2)=FDM/2.0
        A11(NUM1+2,NUM2+3)=-FDM/6.0
        ! current contunuity:
        A11(NUM1+3,NUM1+2)=DIDD/XX(KEL)
        A11(NUM1+3,NUM1+3)=DIDD/XX(KEL)
        A11(NUM1+3,NUM2+2)=-DIDDP/XX(KEL+1)
        A11(NUM1+3,NUM2+3)=DIDDP/XX(KEL+1)
        NUM1=NUM1+3
      ENDDO
      ! left boundary condition:
      IBM=MAT(1)
      IF((IQFR(1,1) == -1).OR.(IQFR(1,1) > 0)) THEN
        ! VOID
        AFACTOR=QFR(1,1)
        A11(NUM1+2,1)=AFACTOR
        A11(NUM1+2,2)=-(AFACTOR/2.0+DIFF(IBM)/XX(1))
        A11(NUM1+2,3)=(AFACTOR/6.0+DIFF(IBM)/XX(1))
      ELSE IF(IQFR(1,1) == -2) THEN
        ! REFL
        A11(NUM1+2,2)=1.0
        A11(NUM1+2,3)=-1.0
      ELSE IF(IQFR(1,1) == -3) THEN
        ! ZERO
        A11(NUM1+2,1)=1.0
        A11(NUM1+2,2)=-1.0/2.0
        A11(NUM1+2,3)=1.0/6.0
      ENDIF
      ! right boundary condition:
      IBM=MAT(NEL)
      IF((IQFR(2,NEL) == -1).OR.(IQFR(2,NEL) > 0)) THEN
        NUM2=3*(NEL-1)
        ! VOID
        AFACTOR=QFR(2,NEL)
        A11(NUM1+3,NUM2+1)=AFACTOR
        A11(NUM1+3,NUM2+2)=(AFACTOR/2.0+DIFF(IBM)/XX(NEL))
        A11(NUM1+3,NUM2+3)=(AFACTOR/6.0+DIFF(IBM)/XX(NEL))
      ELSE IF(IQFR(2,NEL) == -2) THEN
        NUM2=3*(NEL-1)
        ! REFL
        A11(NUM1+3,NUM2+2)=1.0
        A11(NUM1+3,NUM2+3)=1.0
      ELSE IF(IQFR(2,NEL) == -3) THEN
        NUM2=3*(NEL-1)
        ! ZERO
        A11(NUM1+3,NUM2+1)=1.0
        A11(NUM1+3,NUM2+2)=1.0/2.0
        A11(NUM1+3,NUM2+3)=1.0/6.0
      ENDIF
      END SUBROUTINE NSS4TR
