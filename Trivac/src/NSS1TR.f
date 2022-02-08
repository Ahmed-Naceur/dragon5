*DECK NSS1TR
      SUBROUTINE NSS1TR(ITRIAL,NEL,NMIX,MAT,XX,IQFR,QFR,DIFF,SIGR,FD,
     1 A11)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of leakage system matrices for the nodal expansion method.
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
* ITRIAL  type of base (=1: polynomial; =2: hyperbolic).
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
      INTEGER ITRIAL(NMIX),NEL,NMIX,MAT(NEL),IQFR(6,NEL)
      REAL XX(NEL),QFR(6,NEL),DIFF(NMIX),SIGR(NMIX),FD(NMIX,2),
     1 A11(5*NEL,5*NEL)
*
      A11(:5*NEL,:5*NEL)=0.0
      ! WEIGHT RESIDUAL EQUATIONS:
      NUM1=0
      DO KEL=1,NEL
        IBM=MAT(KEL)
        DX2=XX(KEL)**2
        SIGG=SIGR(IBM)
        DIDD=DIFF(IBM)
        ETA=XX(KEL)*SQRT(SIGG/DIDD)
        A11(NUM1+1,NUM1+1)=SIGG
        A11(NUM1+1,NUM1+3)=-6.0*DIDD/DX2
        A11(NUM1+2,NUM1+2)=SIGG/12.0
        A11(NUM1+3,NUM1+3)=SIGG/20.0
        IF(ITRIAL(IBM) == 1) THEN
          A11(NUM1+1,NUM1+5)=-2.0*DIDD/(5.0*DX2)
          A11(NUM1+2,NUM1+4)=-SIGG/120.0-DIDD/(2.0*DX2)
          A11(NUM1+3,NUM1+5)=-SIGG/700.0-DIDD/(5.0*DX2)
        ELSE
          ALP0=2.0*ETA*SINH(ETA/2.0)
          A11(NUM1+1,NUM1+5)=-DIDD*ALP0/DX2
        ENDIF
        NUM1=NUM1+5
      ENDDO
      ! continuity relations:
      NUM1=0
      DO KEL=1,NEL-1
        IBM=MAT(KEL)
        IBMP=MAT(KEL+1)
        DIDD=DIFF(IBM)
        DIDDP=DIFF(IBMP)
        ETA=XX(KEL)*SQRT(SIGR(IBM)/DIDD)
        ETAP=XX(KEL+1)*SQRT(SIGR(IBMP)/DIDDP)
        NUM2=NUM1+5
        ! flux continuity:
        FDP=FD(IBM,2)
        FDM=FD(IBMP,1)
        A11(NUM1+4,NUM1+1)=FDP
        A11(NUM1+4,NUM1+2)=FDP/2.0
        A11(NUM1+4,NUM1+3)=FDP/2.0
        A11(NUM1+4,NUM2+1)=-FDM
        A11(NUM1+4,NUM2+2)=FDM/2.0
        A11(NUM1+4,NUM2+3)=-FDM/2.0
        IF(ITRIAL(IBM) == 2) THEN
          ALP1=ETA*COSH(ETA/2.0)-2.0*SINH(ETA/2.0)
          A11(NUM1+4,NUM1+4)=FDP*SINH(ETA/2.0)
          A11(NUM1+4,NUM1+5)=FDP*ALP1/ETA
        ENDIF
        IF(ITRIAL(IBMP) == 2) THEN
          ALP1P=ETAP*COSH(ETAP/2.0)-2.0*SINH(ETAP/2.0)
          A11(NUM1+4,NUM2+4)=FDM*SINH(ETAP/2.0)
          A11(NUM1+4,NUM2+5)=-FDM*ALP1P/ETAP
        ENDIF
        ! current contunuity:
        A11(NUM1+5,NUM1+2)=DIDD/XX(KEL)
        A11(NUM1+5,NUM1+3)=3.0*DIDD/XX(KEL)
        A11(NUM1+5,NUM2+2)=-DIDDP/XX(KEL+1)
        A11(NUM1+5,NUM2+3)=3.0*DIDDP/XX(KEL+1)
        IF(ITRIAL(IBM) == 1) THEN
          A11(NUM1+5,NUM1+4)=DIDD/(2.0*XX(KEL))
          A11(NUM1+5,NUM1+5)=DIDD/(5.0*XX(KEL))
        ELSE
          A11(NUM1+5,NUM1+4)=(DIDD/XX(KEL))*ETA*COSH(ETA/2.0)
          A11(NUM1+5,NUM1+5)=(DIDD/XX(KEL))*ETA*SINH(ETA/2.0)
        ENDIF
        IF(ITRIAL(IBMP) == 1) THEN
          A11(NUM1+5,NUM2+4)=-DIDDP/(2.0*XX(KEL+1))
          A11(NUM1+5,NUM2+5)=DIDDP/(5.0*XX(KEL+1))
        ELSE
          A11(NUM1+5,NUM2+4)=-(DIDDP/XX(KEL+1))*ETAP*COSH(ETAP/2.0)
          A11(NUM1+5,NUM2+5)=(DIDDP/XX(KEL+1))*ETAP*SINH(ETAP/2.0)
        ENDIF
        NUM1=NUM1+5
      ENDDO
      ! left boundary condition:
      IBM=MAT(1)
      ETA=XX(1)*SQRT(SIGR(MAT(1))/DIFF(IBM))
      IF((IQFR(1,1) == -1).OR.(IQFR(1,1) > 0)) THEN
        ! VOID
        AFACTOR=QFR(1,1)
        A11(NUM1+4,1)=AFACTOR
        A11(NUM1+4,2)=-(AFACTOR/2.0+DIFF(IBM)/XX(1))
        A11(NUM1+4,3)=(AFACTOR/2.0+3.0*DIFF(IBM)/XX(1))
        IF(ITRIAL(IBM) == 1) THEN
          A11(NUM1+4,4)=-DIFF(IBM)/(2.0*XX(1))
          A11(NUM1+4,5)=DIFF(IBM)/(5.0*XX(1))
        ELSE
          ALP1=ETA*COSH(ETA/2.0)-2.0*SINH(ETA/2.0)
          A11(NUM1+4,4)=-(AFACTOR*SINH(ETA/2.0)+(DIFF(IBM)/XX(1))*
     1    ETA*COSH(ETA/2.0))
          A11(NUM1+4,5)=AFACTOR*ALP1/ETA+(DIFF(IBM)/XX(1))*ETA*
     1    SINH(ETA/2.0)
        ENDIF
      ELSE IF(IQFR(1,1) == -2) THEN
        ! REFL
        A11(NUM1+4,2)=1.0
        A11(NUM1+4,3)=-3.0
        IF(ITRIAL(IBM) == 1) THEN
          A11(NUM1+4,4)=1.0/2.0
          A11(NUM1+4,5)=-1.0/5.0
        ELSE
          A11(NUM1+4,4)=ETA*COSH(ETA/2.0)
          A11(NUM1+4,5)=-ETA*SINH(ETA/2.0)
        ENDIF
      ELSE IF(IQFR(1,1) == -3) THEN
        ! ZERO
        A11(NUM1+4,1)=1.0
        A11(NUM1+4,2)=-1.0/2.0
        A11(NUM1+4,3)=1.0/2.0
        IF(ITRIAL(IBM) == 2) THEN
          ALP1=ETA*COSH(ETA/2.0)-2.0*SINH(ETA/2.0)
          A11(NUM1+4,4)=-SINH(ETA/2.0)
          A11(NUM1+4,5)=ALP1/ETA
        ENDIF
      ENDIF
      ! right boundary condition:
      IBM=MAT(NEL)
      ETA=XX(NEL)*SQRT(SIGR(IBM)/DIFF(IBM))
      IF((IQFR(2,NEL) == -1).OR.(IQFR(2,NEL) > 0)) THEN
        NUM2=5*(NEL-1)
        ! VOID
        AFACTOR=QFR(2,NEL)
        A11(NUM1+5,NUM2+1)=AFACTOR
        A11(NUM1+5,NUM2+2)=(AFACTOR/2.0+DIFF(IBM)/XX(NEL))
        A11(NUM1+5,NUM2+3)=(AFACTOR/2.0+3.0*DIFF(IBM)/XX(NEL))
        IF(ITRIAL(IBM) == 1) THEN
          A11(NUM1+5,NUM2+4)=DIFF(IBM)/(2.0*XX(NEL))
          A11(NUM1+5,NUM2+5)=DIFF(IBM)/(5.0*XX(NEL))
        ELSE
          ALP1=ETA*COSH(ETA/2.0)-2.0*SINH(ETA/2.0)
          A11(NUM1+5,NUM2+4)=AFACTOR*SINH(ETA/2.0)+(DIFF(IBM)/
     1    XX(NEL))*ETA*COSH(ETA/2.0)
          A11(NUM1+5,NUM2+5)=AFACTOR*ALP1/ETA+(DIFF(IBM)/
     1    XX(NEL))*ETA*SINH(ETA/2.0)
        ENDIF
      ELSE IF(IQFR(2,NEL) == -2) THEN
        NUM2=5*(NEL-1)
        ! REFL
        A11(NUM1+5,NUM2+2)=1.0
        A11(NUM1+5,NUM2+3)=3.0
        IF(ITRIAL(IBM) == 1) THEN
          A11(NUM1+5,NUM2+4)=1.0/2.0
          A11(NUM1+5,NUM2+5)=1.0/5.0
        ELSE
          A11(NUM1+5,NUM2+4)=ETA*COSH(ETA/2.0)
          A11(NUM1+5,NUM2+5)=ETA*SINH(ETA/2.0)
        ENDIF
      ELSE IF(IQFR(2,NEL) == -3) THEN
        NUM2=5*(NEL-1)
        ! ZERO
        A11(NUM1+5,NUM2+1)=1.0
        A11(NUM1+5,NUM2+2)=1.0/2.0
        A11(NUM1+5,NUM2+3)=1.0/2.0
        IF(ITRIAL(IBM) == 2) THEN
          ALP1=ETA*COSH(ETA/2.0)-2.0*SINH(ETA/2.0)
          A11(NUM1+5,NUM2+4)=SINH(ETA/2.0)
          A11(NUM1+5,NUM2+5)=ALP1/ETA
        ENDIF
      ENDIF
      END SUBROUTINE NSS1TR
