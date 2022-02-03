*DECK NSSFL2
      SUBROUTINE NSSFL2(NUN,NG,LX1,NMIX,NALB,EPSOUT,MAXOUT,MAT,
     1 XX,IQFR,QFR,DIFF,SIGR,CHI,SIGF,SCAT,BETA,FD,IPRINT,EVAL,EVECT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Flux calculation for the coarse mesh finite differences method in
* Cartesian 1D geometry.
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
* NUN     number of unknowns (=3*LX1).
* NG      number of energy groups.
* LX1     number of nodes in the nodal calculation.
* NMIX    number of mixtures in the nodal calculation.
* NALB    number of physical albedos.
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
      INTEGER NUN,NG,LX1,NMIX,NALB,MAXOUT,IPRINT,MAT(LX1),IQFR(6,LX1)
      REAL EPSOUT,XX(LX1),QFR(6,LX1),DIFF(NMIX,NG),SIGR(NMIX,NG),
     1 CHI(NMIX,NG),SIGF(NMIX,NG),SCAT(NMIX,NG,NG),BETA(NALB,NG,NG),
     2 FD(NMIX,2,NG,NG),EVAL,EVECT(NUN,NG)
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:,:) :: A,B,AI,A11,QFR2,COUR
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*----
* SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(COUR(LX1+1,NG),A(NUN*NG,NUN*NG),B(NUN*NG,NUN*NG))
      ALLOCATE(A11(NUN,NUN),QFR2(6,LX1))
*----
*  COMPUTE NODAL SOLUTION
*----
      A(:NUN*NG,:NUN*NG)=0.0
      B(:NUN*NG,:NUN*NG)=0.0
      QFR2(:6,:LX1)=0.0
      DO J=1,NG
        IOF1=(J-1)*NUN
        DO I=1,NG
          DO IQW=1,2
            DO IEL=1,LX1
              IALB=IQFR(IQW,IEL)
              IF(IALB.GT.0) THEN
                IF(IALB.GT.NALB) CALL XABORT('NSSFL2: BETA OVERFLOW.')
                QFR2(IQW,IEL)=QFR(IQW,IEL)*ALB(BETA(IALB,I,J))
              ELSE IF(I == J) THEN
                QFR2(IQW,IEL)=QFR(IQW,IEL)
              ELSE
                QFR2(IQW,IEL)=0.0
              ENDIF
            ENDDO
          ENDDO
          IOF2=(I-1)*NUN
          IF(I == J) THEN
            CALL NSS4TR(LX1,NMIX,MAT,XX,IQFR,QFR2,DIFF(:,I),SIGR(:,I),
     1      FD(:,:,I,J),A11)
            A(IOF1+1:IOF1+NUN,IOF1+1:IOF1+NUN)=A11(:,:)
          ELSE
            CALL NSS5TR(LX1,NMIX,MAT,IQFR,QFR2,SCAT(:,I,J),FD(:,:,I,J),
     1      A11)
            A(IOF2+1:IOF2+NUN,IOF1+1:IOF1+NUN)=-A11(:,:)
          ENDIF
          B(IOF2+1:IOF2+NUN,IOF1+1:IOF1+NUN)=0.0
          NUM1=0
          DO IEL=1,LX1
            IBM=MAT(IEL)
            B(IOF2+NUM1+1,IOF1+NUM1+1)=CHI(IBM,I)*SIGF(IBM,J)
            NUM1=NUM1+3
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(QFR2,A11)
*----
*  SOLVE EIGENVALUE MATRIX SYSTEM
*----
      EPS=1.0E-7
      CALL ALINV(NUN*NG,A,NUN*NG,IER)
      IF(IER.NE.0) CALL XABORT('NSSFL2: SINGULAR MATRIX')
      ALLOCATE(AI(NUN*NG,NUN*NG))
      AI(:NUN*NG,:NUN*NG)=MATMUL(A(:NUN*NG,:NUN*NG),B(:NUN*NG,:NUN*NG))
      EVECT(:,:)=0.0
      NUM1=0
      DO IEL=1,LX1
        EVECT(NUM1+1,:)=1.0
        NUM1=NUM1+3
      ENDDO
      CALL AL1EIG(NUN*NG,AI,EPSOUT,MAXOUT,ITER,EVECT,EVAL,IPRINT)
      IF(IPRINT.GT.0) WRITE(6,10) EVAL,ITER
      DEALLOCATE(AI,B,A)
*----
*  NORMALIZE THE FLUX
*----
      FLMAX=0.0
      DO IG=1,NG
        NUM1=0
        DO IEL=1,LX1
          IF(ABS(EVECT(NUM1+1,IG)).GT.ABS(FLMAX)) FLMAX=EVECT(NUM1+1,IG)
          NUM1=NUM1+3
        ENDDO
      ENDDO
      EVECT(:,:)=EVECT(:,:)/FLMAX
*----
*  COMPUTE CURRENTS
*----
      DO IG=1,NG
        DO KEL=1,LX1
          IBM=MAT(KEL)
          IOF=(KEL-1)*3
          COUR(KEL,IG)=-(DIFF(IBM,IG)/XX(KEL))*(EVECT(IOF+2,IG)-
     1    3.0*EVECT(IOF+3,IG))
        ENDDO
        IBM=MAT(LX1)
        IOF=(LX1-1)*3
        COUR(LX1+1,IG)=-(DIFF(IBM,IG)/XX(LX1))*(EVECT(IOF+2,IG)+
     1  3.0*EVECT(IOF+3,IG))
        IF(IPRINT.GT.0) THEN
          WRITE(6,'(/33H NSSFL2: AVERAGED FLUXES IN GROUP,I5)') IG
          WRITE(6,'(1P,10e12.4)') (EVECT((I-1)*3+1,IG),I=1,LX1)
          WRITE(6,'(/39H NSSFL2: SURFACIC NET CURRENTS IN GROUP,I5)') IG
          WRITE(6,'(1P,10e12.4)') (COUR(I,IG),I=1,LX1+1)
        ENDIF
      ENDDO
      DEALLOCATE(COUR)
      RETURN
*
   10 FORMAT(14H NSSFL2: KEFF=,F11.8,12H OBTAINED IN,I5,11H ITERATIONS)
      END
