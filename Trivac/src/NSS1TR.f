*DECK NSS1TR
      SUBROUTINE NSS1TR(ITRIAL,NEL,NMIX,MAT,XX,KN,QFR,DIFF,SIGR,FD,A11)
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
* DIFF    diffusion coefficients
* SIGR    macroscopic cross sections
* FD      discontinuity factors
*
*Parameters: output
* A11     assembly matrix.
*
*-----------------------------------------------------------------------
*
      INTEGER ITRIAL(NMIX),NEL,NMIX,MAT(NEL),KN(6,NEL)
      REAL XX(NEL),QFR(6,NEL),DIFF(NMIX),SIGR(NMIX),FD(NMIX,2),
     1 A11(5*NEL,5*NEL)
*
      A11(:5*NEL,:5*NEL)=0.0
      ! WEIGHT RESIDUAL EQUATIONS:
      NUM1=0
      DO KEL=1,NEL
        IBM=MAT(KEL)
        DX2=XX(KEL)**2
        SIGG=SIGR(MAT(KEL))
        DIDD=DIFF(MAT(KEL))
        ETA=XX(KEL)*SQRT(SIGG/DIDD)
        A11(NUM1+1,NUM1+1)=SIGG
        A11(NUM1+1,NUM1+3)=-2.0*DIDD/DX2
        A11(NUM1+2,NUM1+2)=SIGG/12.0
        A11(NUM1+3,NUM1+3)=SIGG/180.0
        IF (ITRIAL(IBM) == 1) THEN
          A11(NUM1+1,NUM1+5)=-2.0*DIDD/(5.0*DX2)
          A11(NUM1+2,NUM1+4)=-SIGG/120.0-DIDD/(2.0*DX2)
          A11(NUM1+3,NUM1+5)=-SIGG/2100.0-DIDD/(15.0*DX2)
        ELSE
          ALP0=2.0*ETA*SINH(ETA/2)
          ALP1=ETA*COSH(ETA/2)-2.0*SINH(ETA/2)
          ALP2=((12.0+ETA**2)*SINH(ETA/2)-6.0*ETA*COSH(ETA/2))/(3.0*ETA)
          A11(NUM1+1,NUM1+5)=-DIDD*ALP0/DX2
          A11(NUM1+2,NUM1+4)=(SIGG/(ETA**2)-DIDD/DX2)*ALP1
          A11(NUM1+3,NUM1+5)=(SIGG/(ETA**2)-DIDD/DX2)*ALP2
        ENDIF
        NUM1=NUM1+5
      ENDDO
      ! continuity relations:
      NUM1=0
      DO KEL=1,NEL-1
        IBM=MAT(KEL)
        IBMP=MAT(KEL+1)
        DIDD=DIFF(MAT(KEL))
        DIDDP=DIFF(MAT(KEL+1))
        ETA=XX(KEL)*SQRT(SIGR(MAT(KEL))/DIDD)
        ETAP=XX(KEL+1)*SQRT(SIGR(MAT(KEL+1))/DIDDP)
        NUM2=NUM1+5
        ! flux continuity:
        FDP=FD(MAT(KEL),2)
        FDM=FD(MAT(KEL+1),1)
        A11(NUM1+4,NUM1+1)=FDP
        A11(NUM1+4,NUM1+2)=FDP/2.0
        A11(NUM1+4,NUM1+3)=FDP/6.0
        A11(NUM1+4,NUM2+1)=-FDM
        A11(NUM1+4,NUM2+2)=FDM/2.0
        A11(NUM1+4,NUM2+3)=-FDM/6.0
        IF (ITRIAL(IBM) == 2) THEN
          ALP1=ETA*COSH(ETA/2)-2.0*SINH(ETA/2)
          A11(NUM1+4,NUM1+4)=FDP*SINH(ETA/2)
          A11(NUM1+4,NUM1+5)=FDP*ALP1/ETA
        ENDIF
        IF (ITRIAL(IBMP) == 2) THEN
          ALP1P=ETAP*COSH(ETAP/2)-2.0*SINH(ETAP/2)
          A11(NUM1+4,NUM2+4)=FDM*SINH(ETAP/2)
          A11(NUM1+4,NUM2+5)=-FDM*ALP1P/ETAP
        ENDIF
        ! current contunuity:
        A11(NUM1+5,NUM1+2)=DIDD/XX(KEL)
        A11(NUM1+5,NUM1+3)=DIDD/XX(KEL)
        A11(NUM1+5,NUM2+2)=-DIDDP/XX(KEL+1)
        A11(NUM1+5,NUM2+3)=DIDDP/XX(KEL+1)
        IF (ITRIAL(IBM) == 1) THEN
          A11(NUM1+5,NUM1+4)=DIDD/(2.0*XX(KEL))
          A11(NUM1+5,NUM1+5)=DIDD/(5.0*XX(KEL))
        ELSE
          A11(NUM1+5,NUM1+4)=(DIDD/XX(KEL))*ETA*COSH(ETA/2)
          A11(NUM1+5,NUM1+5)=(DIDD/XX(KEL))*ETA*SINH(ETA/2)
        ENDIF
        IF (ITRIAL(IBMP) == 1) THEN
          A11(NUM1+5,NUM2+4)=-DIDDP/(2.0*XX(KEL+1))
          A11(NUM1+5,NUM2+5)=DIDDP/(5.0*XX(KEL+1))
        ELSE
          A11(NUM1+5,NUM2+4)=-(DIDDP/XX(KEL+1))*ETAP*COSH(ETAP/2)
          A11(NUM1+5,NUM2+5)=(DIDDP/XX(KEL+1))*ETAP*SINH(ETAP/2)
        ENDIF
        NUM1=NUM1+5
      ENDDO
      ! left boundary condition:
      IBM=MAT(1)
      ETA=XX(1)*SQRT(SIGR(MAT(1))/DIFF(MAT(1)))
      IF (KN(1,1) == -1) THEN
        ! VOID
        AFACTOR=QFR(1,1)
        A11(NUM1+4,1)=AFACTOR
        A11(NUM1+4,2)=-(AFACTOR/2.0+DIFF(MAT(1))/XX(1))
        A11(NUM1+4,3)=(AFACTOR/6.0+DIFF(MAT(1))/XX(1))
        IF (ITRIAL(IBM) == 1) THEN
          A11(NUM1+4,4)=-DIFF(MAT(1))/(2.0*XX(1))
          A11(NUM1+4,5)=DIFF(MAT(1))/(5.0*XX(1))
        ELSE
          ALP1=ETA*COSH(ETA/2)-2.0*SINH(ETA/2)
          A11(NUM1+4,4)=-(AFACTOR*SINH(ETA/2)+(DIFF(MAT(1))/XX(1))*
     1    ETA*COSH(ETA/2))
          A11(NUM1+4,5)=AFACTOR*ALP1/ETA+(DIFF(MAT(1))/XX(1))*ETA*
     1    SINH(ETA/2)
        ENDIF
      ELSE IF (KN(1,1) == -2) THEN
        ! REFL
        A11(NUM1+4,2)=1.0
        A11(NUM1+4,3)=-1.0
        IF (ITRIAL(IBM) == 1) THEN
          A11(NUM1+4,4)=1.0/2.0
          A11(NUM1+4,5)=-1.0/5.0
        ELSE
          A11(NUM1+4,4)=ETA*COSH(ETA/2)
          A11(NUM1+4,5)=-ETA*SINH(ETA/2)
        ENDIF
      ELSE IF (KN(1,1) == -3) THEN
        ! ZERO
        A11(NUM1+4,1)=1.0
        A11(NUM1+4,2)=-1.0/2.0
        A11(NUM1+4,3)=1.0/6.0
        IF (ITRIAL(IBM) == 2) THEN
          ALP1=ETA*COSH(ETA/2)-2.0*SINH(ETA/2)
          A11(NUM1+4,4)=-SINH(ETA/2)
          A11(NUM1+4,5)=ALP1/ETA
        ENDIF
      ENDIF
      ! right boundary condition:
      IBM=MAT(NEL)
      ETA=XX(NEL)*SQRT(SIGR(MAT(NEL))/DIFF(MAT(NEL)))
      IF (KN(2,NEL) == -1) THEN
        NUM2=5*(NEL-1)
        ! VOID
        AFACTOR=QFR(2,NEL)
        A11(NUM1+5,NUM2+1)=AFACTOR
        A11(NUM1+5,NUM2+2)=(AFACTOR/2.0+DIFF(MAT(NEL))/XX(NEL))
        A11(NUM1+5,NUM2+3)=(AFACTOR/6.0+DIFF(MAT(NEL))/XX(NEL))
        IF (ITRIAL(IBM) == 1) THEN
          A11(NUM1+5,NUM2+4)=DIFF(MAT(NEL))/(2.0*XX(NEL))
          A11(NUM1+5,NUM2+5)=DIFF(MAT(NEL))/(5.0*XX(NEL))
        ELSE
          ALP1=ETA*COSH(ETA/2)-2.0*SINH(ETA/2)
          A11(NUM1+5,NUM2+4)=AFACTOR*SINH(ETA/2)+(DIFF(MAT(NEL))/
     1    XX(NEL))*ETA*COSH(ETA/2)
          A11(NUM1+5,NUM2+5)=AFACTOR*ALP1/ETA+(DIFF(MAT(NEL))/
     1    XX(NEL))*ETA*SINH(ETA/2)
        ENDIF
      ELSE IF (KN(2,NEL) == -2) THEN
        NUM2=5*(NEL-1)
        ! REFL
        A11(NUM1+5,NUM2+2)=1.0
        A11(NUM1+5,NUM2+3)=1.0
        IF (ITRIAL(IBM) == 1) THEN
          A11(NUM1+5,NUM2+4)=1.0/2.0
          A11(NUM1+5,NUM2+5)=1.0/5.0
        ELSE
          A11(NUM1+5,NUM2+4)=ETA*COSH(ETA/2)
          A11(NUM1+5,NUM2+5)=ETA*SINH(ETA/2)
        ENDIF
      ELSE IF (KN(2,NEL) == -3) THEN
        NUM2=5*(NEL-1)
        ! ZERO
        A11(NUM1+5,NUM2+1)=1.0
        A11(NUM1+5,NUM2+2)=1.0/2.0
        A11(NUM1+5,NUM2+3)=1.0/6.0
        IF (ITRIAL(IBM) == 2) THEN
          ALP1=ETA*COSH(ETA/2)-2.0*SINH(ETA/2)
          A11(NUM1+5,NUM2+4)=SINH(ETA/2)
          A11(NUM1+5,NUM2+5)=ALP1/ETA
        ENDIF
      ENDIF
      END SUBROUTINE NSS1TR
