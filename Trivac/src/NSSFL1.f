*DECK NSSFL1
      SUBROUTINE NSSFL1(NUN,NG,LX1,NMIX,NALB,ITRIAL,EPSOUT,MAXOUT,MAT,
     1 XX,IQFR,QFR,DIFF,SIGR,CHI,SIGF,SCAT,BETA,FD,IPRINT,EVAL,EVECT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Flux calculation for the nodal expansion method in Cartesian 1D
* geometry.
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
* NUN     number of unknowns (=4*LX1+1).
* NG      number of energy groups.
* LX1     number of nodes in the nodal calculation.
* NMIX    number of mixtures in the nodal calculation.
* NALB    number of physical albedos.
* ITRIAL  type of expansion functions in the nodal calculation
*         (=1: polynomial; =2: hyperbolic).
* EPSOUT  convergence epsilon for the power method.
* MAXOUT  maximum number of iterations for the power method.
* MAT     material mixtures.
* XX      mesh spacings.
* IQFR    boundary condition information.
* QFR     albedo function information.
* DIFF    diffusion coefficients
* SIGR    removal cross sections.
* CHI     fission spectra.
* SIGF    nu times fission cross section.
* SCAT    scattering cross section.
* BETA    albedos.
* FD      discontinuity factors.
* IPRINT  edition flag.
*
*Parameters: output
* EVAL    effective multiplication factor
* EVECT   neutron flux
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NUN,NG,LX1,NMIX,NALB,ITRIAL(NMIX,NG),MAXOUT,IPRINT,
     1 MAT(LX1),IQFR(6,LX1)
      REAL EPSOUT,XX(LX1),QFR(6,LX1),DIFF(NMIX,NG),SIGR(NMIX,NG),
     1 CHI(NMIX,NG),SIGF(NMIX,NG),SCAT(NMIX,NG,NG),BETA(NALB,NG,NG),
     2 FD(NMIX,2,NG,NG),EVAL,EVECT(NUN,NG)
*----
*  LOCAL VARIABLES
*----
      INTEGER DIM
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK
      REAL, ALLOCATABLE, DIMENSION(:,:) :: A,B,AI,A11,QFR2,FUNKN
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      DIM=5*LX1
      ALLOCATE(FUNKN(DIM,NG),A(DIM*NG,DIM*NG),B(DIM*NG,DIM*NG))
      ALLOCATE(WORK(NMIX),A11(DIM,DIM),QFR2(6,LX1))
*----
*  COMPUTE NODAL SOLUTION
*----
      DIM=5*LX1
      A(:DIM*NG,:DIM*NG)=0.0
      B(:DIM*NG,:DIM*NG)=0.0
      QFR2(:6,:LX1)=0.0
      DO J=1,NG
        IOF1=(J-1)*DIM
        DO I=1,NG
          DO IQW=1,2
            DO IEL=1,LX1
              IALB=IQFR(IQW,IEL)
              IF(IALB.GT.0) THEN
                IF(IALB.GT.NALB) CALL XABORT('NSSFL1: BETA OVERFLOW.')
                QFR2(IQW,IEL)=QFR(IQW,IEL)*ALB(BETA(IALB,I,J))
              ELSE IF(I == J) THEN
                QFR2(IQW,IEL)=QFR(IQW,IEL)
              ELSE
                QFR2(IQW,IEL)=0.0
              ENDIF
            ENDDO
          ENDDO
          DO IBM=1,NMIX
            WORK(IBM)=CHI(IBM,I)*SIGF(IBM,J)
          ENDDO
          IOF2=(I-1)*DIM
          IF(I == J) THEN
            CALL NSS1TR(ITRIAL(1,J),LX1,NMIX,MAT,XX,IQFR,QFR2,DIFF(:,I),
     1      SIGR(:,I),FD(:,:,I,J),A11)
            A(IOF1+1:IOF1+DIM,IOF1+1:IOF1+DIM)=A11(:,:)
          ELSE
            CALL NSS2TR(ITRIAL(1,J),LX1,NMIX,MAT,XX,IQFR,QFR2,DIFF(:,J),
     1      SIGR(:,J),SCAT(:,I,J),FD(:,:,I,J),A11)
            A(IOF2+1:IOF2+DIM,IOF1+1:IOF1+DIM)=-A11(:,:)
          ENDIF
          CALL NSS3TR(ITRIAL(1,J),LX1,NMIX,MAT,XX,DIFF(:,J),SIGR(:,J),
     1    WORK(:),A11)
          B(IOF2+1:IOF2+DIM,IOF1+1:IOF1+DIM)=A11(:,:)
        ENDDO
      ENDDO
      DEALLOCATE(QFR2,A11,WORK)
*----
*  SOLVE EIGENVALUE MATRIX SYSTEM
*----
      CALL ALINV(DIM*NG,A,DIM*NG,IER)
      IF(IER.NE.0) CALL XABORT('NSSFL1: SINGULAR MATRIX')
      ALLOCATE(AI(DIM*NG,DIM*NG))
      AI(:DIM*NG,:DIM*NG)=MATMUL(A(:DIM*NG,:DIM*NG),B(:DIM*NG,:DIM*NG))
      FUNKN(:,:)=0.0
      NUM1=0
      DO IEL=1,LX1
        FUNKN(NUM1+1,:)=1.0
        NUM1=NUM1+5
      ENDDO
      CALL AL1EIG(DIM*NG,AI,EPSOUT,MAXOUT,ITER,FUNKN,EVAL,IPRINT)
      IF(IPRINT.GT.0) WRITE(6,10) EVAL,ITER
      DEALLOCATE(AI,B,A)
*----
*  NORMALIZE THE FLUX
*----
      FLMAX=0.0
      DO IG=1,NG
        NUM1=0
        DO IEL=1,LX1
          IF(ABS(FUNKN(NUM1+1,IG)).GT.ABS(FLMAX)) FLMAX=FUNKN(NUM1+1,IG)
          NUM1=NUM1+5
        ENDDO
      ENDDO
      FUNKN(:,:)=FUNKN(:,:)/FLMAX
*----
*  COMPUTE INTERFACE FLUXES AND CURRENTS
*----
      IOF1=LX1
      IOF2=2*LX1
      IOF3=3*LX1
      IF(IOF3+LX1+1.NE.NUN) CALL XABORT('NSSFL1: NUN ERROR.')
      DO IG=1,NG
        DO KEL=1,LX1
          IBM=MAT(KEL)
          IOF=(KEL-1)*5
          EVECT(KEL,IG)=FUNKN(IOF+1,IG)
          EVECT(IOF1+KEL,IG)=FUNKN(IOF+1,IG)+0.5*(-FUNKN(IOF+2,IG)+
     1    FUNKN(IOF+3,IG))
          EVECT(IOF2+KEL,IG)=FUNKN(IOF+1,IG)+0.5*(FUNKN(IOF+2,IG)+
     1    FUNKN(IOF+3,IG))
          IF(ITRIAL(IBM,IG).EQ.1) THEN
            EVECT(IOF3+KEL,IG)=-(DIFF(IBM,IG)/XX(KEL))*(FUNKN(IOF+2,IG)-
     1      3.0*FUNKN(IOF+3,IG)+FUNKN(IOF+4,IG)/2.0-FUNKN(IOF+5,IG)/5.0)
          ELSE
            ETA=XX(KEL)*SQRT(SIGR(IBM,IG)/DIFF(IBM,IG))
            ALP1=ETA*COSH(ETA/2.0)-2.0*SINH(ETA/2.0)
            EVECT(IOF1+KEL,IG)=EVECT(IOF1+KEL,IG)-FUNKN(IOF+4,IG)*
     1      SINH(ETA/2.0)+FUNKN(IOF+5,IG)*ALP1/ETA
            EVECT(IOF2+KEL,IG)=EVECT(IOF2+KEL,IG)+FUNKN(IOF+4,IG)*
     1      SINH(ETA/2.0)+FUNKN(IOF+5,IG)*ALP1/ETA
            EVECT(IOF3+KEL,IG)=-(DIFF(IBM,IG)/XX(KEL))*(FUNKN(IOF+2,IG)-
     1      3.0*FUNKN(IOF+3,IG)+FUNKN(IOF+4,IG)*ETA*COSH(ETA/2.0)-
     2      FUNKN(IOF+5,IG)*ETA*SINH(ETA/2.0))
          ENDIF
        ENDDO
        IBM=MAT(LX1)
        IOF=(LX1-1)*5
        IF(ITRIAL(IBM,IG).EQ.1) THEN
          EVECT(IOF3+LX1+1,IG)=-(DIFF(IBM,IG)/XX(LX1))*(FUNKN(IOF+2,IG)+
     1    3.0*FUNKN(IOF+3,IG)+FUNKN(IOF+4,IG)/2.0+FUNKN(IOF+5,IG)/5.0)
        ELSE
          ETA=XX(LX1)*SQRT(SIGR(IBM,IG)/DIFF(IBM,IG))
          EVECT(IOF3+LX1+1,IG)=-(DIFF(IBM,IG)/XX(LX1))*(FUNKN(IOF+2,IG)+
     1    3.0*FUNKN(IOF+3,IG)+FUNKN(IOF+4,IG)*ETA*COSH(ETA/2.0)+
     2    FUNKN(IOF+5,IG)*ETA*SINH(ETA/2.0))
        ENDIF
        IF(IPRINT.GT.0) THEN
          WRITE(6,'(/33H NSSFL1: AVERAGED FLUXES IN GROUP,I5)') IG
          WRITE(6,'(1P,10e12.4)') (EVECT(I,IG),I=1,LX1)
          WRITE(6,'(/39H NSSFL1: SURFACIC NET CURRENTS IN GROUP,I5)') IG
          WRITE(6,'(1P,10e12.4)') (EVECT(IOF3+I,IG),I=1,LX1+1)
        ENDIF
      ENDDO
      DEALLOCATE(FUNKN)
      RETURN
*
   10 FORMAT(14H NSSFL1: KEFF=,F11.8,12H OBTAINED IN,I5,11H ITERATIONS)
      END
