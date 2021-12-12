*DECK XDRH13
      SUBROUTINE XDRH13 (IR1,NMILG,NG,NSMAX,IQUAD,FRTM,NS,IDIL,MIXGR,
     1 RS,FRACT,SIGMA,SIGMS,P1I,P1DI,P1KI,SIGA1)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the reduced collision probabilities for the She-Liu-Shi
* double heterogeneity model (part 1).
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Chambon
*
*Parameters: input
* IR1     number of elementary mixtures in the domain.
* NMILG   number of composite mixtures in the domain.
* NG      number of different kind of micro structures. A kind of
*         micro structure is characterized by the radius of its
*         micro volumes. All the micro volumes of the same kind
*         should own the same nuclear properties in a given macro
*         volume.
* NSMAX   maximum number of volumes (tubes or shells) in each kind of
*         micro structure.
* IQUAD   quadrature parameter for the treatment of the micro volumes.
*         if IQUAD < 0, lines with regular interval are applied.
* FRTM    minimum volume fraction of the grain in the representative 
*         volume for She-Liu-Shi models.
* NS      number of volumes in each kind of micro structure.
* IDIL    elementary mixture indices in the diluent of the composite
*         mixtures.
* MIXGR   elementary mixture indices in the micro structures.
* RS      radius of the micro volumes.
* FRACT   volumic fractions of the micro volumes.
*
*Parameters: input/output
* SIGMA   total macroscopic cross sections in each mixture of the
*         composite geometry.
* SIGMS   scattering macroscopic cross sections in each mixture of the
*         composite geometry.
*
*Parameters: output
* P1I     non collision probability in subvolume with 1 grain type.    
* P1KI    escape probability from layer k in subvolume 
*         with 1 grain type.
* P1DI    escape probability from matrix in subvolume 
*         with 1 grain type.
* SIGA1   output cross sections.
*
*Reference:
* D. She, Z. Liu, and L. Shi, An Equivalent Homogenization Method for
* Treating the Stochastic Media, Nucl. Sci. Eng., 185,  351-360 (2018)
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR1,NMILG,NG,NSMAX,IQUAD,NS(NG),IDIL(NMILG),
     1 MIXGR(NSMAX,NG,NMILG)
      REAL FRTM,RS(NSMAX+1,NG),FRACT(NG,IR1+NMILG),SIGMA(0:IR1+NMILG),
     1 SIGMS(0:IR1+NMILG),P1I(NG,NMILG),P1DI(NG,NMILG),
     2 P1KI(NSMAX,NG,NMILG),SIGA1(NG,NMILG)
*----
*  LOCAL VARIABLES
*----
      INTEGER NR,IBM,J,K,N,M,MIL
      REAL DILF,DX,DXFACT,DRMIN,X,EP,EPI1,EPI2,LM,LGAR,SIGMA1,SIGMS1,
     > FRT,RMAX,SIGT,XI,FRTT
      DOUBLE PRECISION P1,P1D
      DOUBLE PRECISION P1K(NSMAX)
      REAL RGAR(NSMAX+2),LR(NSMAX+1)
      REAL, ALLOCATABLE, DIMENSION(:) :: SIG
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SIG(NSMAX))
*----
*  COMPUTE THE EQUIVALENT TOTAL AND SCATTERING CROSS SECTIONS 
*  IN COMPOSITE REGIONS
*----
      DXFACT=100.0
      IF(IQUAD.LT.0) DXFACT=REAL(-IQUAD)
      DO 180 IBM=1,NMILG
        MIL=IR1+IBM
        SIGT=SIGMA(IDIL(IBM))
        DILF=1.0
        SIGMA(IR1+IBM)=0.0
        SIGMS(IR1+IBM)=0.0
        DO 10 J=1,NG
          DILF=DILF-FRACT(J,MIL)
   10   CONTINUE
        DO 130 J=1,NG
          FRT=FRACT(J,MIL)
          IF(FRT.LE.0.00001) GO TO 130
          NR=NS(J)
          DRMIN=RS(NS(J)+1,J)
          RGAR(1)=0.0
          DO 15 K=2,NS(J)+1
            DRMIN=MIN(DRMIN,RS(K,J)-RS(K-1,J))
            RGAR(K)=RS(K,J)
   15     CONTINUE
          DO 20 K=1,NS(J)
            P1K(K)=0.0D0
            SIG(K)=SIGMA(MIXGR(K,J,IBM))
   20     CONTINUE
          FRTT=1.0-DILF
*         FRT too small -> additional ring of matrix
          IF((1.0-DILF).LT.FRTM) THEN
            NR=NR+1
            RGAR(NR+1)=RGAR(NR)*(FRTM/FRT)**(1.0/3.0)   
            SIG(NR)=SIGT
            FRTT=FRTM
          ENDIF
          RMAX=RGAR(NR+1)
          LGAR=4.0/3.0*RMAX/FRTT
          P1=0.0D0
          P1D=0.0D0
          DX=DRMIN/DXFACT
          XI=-0.5
          K=1
* integral over radius to compute collision prob.
   30     XI=XI+1.0
          X=DX*XI
          IF (X.GT.RGAR(K+1)) K=K+1
          IF (K.GT.NR) GO TO 100
*    Ref 1): Eq 13-17 
*     compute segment lengths 
          LM=LGAR/2.0
          DO 40 N=1,NR
            IF (N.LT.K) THEN
              LR(N)=0.0D0
            ELSEIF (N.EQ.K) THEN
              LR(N)=(RGAR(N+1)**2.0 - X**2.0)**0.5
            ELSE 
              LR(N)=(RGAR(N+1)**2.0 - X**2.0)**0.5
     1              -(RGAR(N)**2.0 - X**2.0)**0.5
            ENDIF
            LM=LM-LR(N)  
   40     CONTINUE
*     Ref 1): Eq 18-19
          EP=2*SIGT*LM
          DO 50 N=1,NR
            EP=EP+2*LR(N)*SIG(N)
   50     CONTINUE
          P1=P1+X*DX*EXP(-EP)
          DO 70 N=K,NR 
            EPI1=SIGT*LM
            EPI2=SIGT*LM
            DO 60 M=1,NR
              IF (M.LT.N) THEN
                EPI2=EPI2+2*LR(M)*SIG(M)
              ELSEIF (M.EQ.N) THEN
                EPI2=EPI2+LR(M)*SIG(M)
              ELSE 
                EPI1=EPI1+LR(M)*SIG(M)
                EPI2=EPI2+LR(M)*SIG(M)
              ENDIF
   60       CONTINUE  
            P1K(N)=P1K(N)+X*DX*(EXP(-EPI1)+EXP(-EPI2))*
     1              (1.0D0-EXP(-LR(N)*SIG(N)))
   70     CONTINUE
          GO TO 30
  100     CONTINUE
          P1=P1*2/RMAX**2.0
          P1I(J,IBM)=REAL(P1)
          P1D=1.0D0-P1
          DO 110 K=1,NS(J)
            P1K(K)=P1K(K)*2/RMAX**2.0
            P1KI(K,J,IBM)=REAL(P1K(K))
*      collision prob. conservation, Ref 1): Eq 4
            P1D=P1D-P1K(K)
  110      CONTINUE
          P1DI(J,IBM)=REAL(P1D)
*     Ref 1): Eq 5
          SIGMA1=REAL(-LOG(P1)/LGAR)
          SIGA1(J,IBM)=SIGMA1
          SIGMS1=REAL(P1D/(1.0-P1)*SIGMA1/SIGT*SIGMS(IDIL(IBM)))
          DO 120 K=1,NS(J)
            SIGMS1=REAL(SIGMS1+P1K(K)/(1.0-P1)*SIGMA1/SIG(K)*
     1               SIGMS(MIXGR(K,J,IBM)))
  120      CONTINUE
*     Ref 1): Eq 26
          SIGMA(IR1+IBM)=SIGMA(IR1+IBM)+SIGMA1*FRT/(1-DILF)
          SIGMS(IR1+IBM)=SIGMS(IR1+IBM)+SIGMS1*FRT/(1-DILF)
  130     CONTINUE
  180   CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SIG)
      RETURN
      END
