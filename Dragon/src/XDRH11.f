*DECK XDRH11
      SUBROUTINE XDRH11 (IR1,NMILG,NG,NSMAX,MICRO,IQUAD,NS,IDIL,MIXGR,
     1 RS,FRACT,VOLK,SIGMA,SIGMS,NCO,RRRR,QKOLD,QKDEL,PKL,COEF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the reduced collision probabilities for the Sanchez-
* Pomraning double heterogeneity model (part 1).
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
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
* MICRO   type of micro volumes (=3 cylinder; =4 sphere).
* IQUAD   quadrature parameter for the treatment of the micro volumes.
* NS      number of volumes in each kind of micro structure.
* IDIL    elementary mixture indices in the diluent of the composite
*         mixtures.
* MIXGR   elementary mixture indices in the micro structures.
* RS      radius of the micro volumes.
* FRACT   volumic fractions of the micro volumes.
* VOLK    volumic fractions of the tubes or shells in the micro volumes.
*
*Parameters: input/output
* SIGMA   total macroscopic cross sections in each mixture of the
*         composite geometry.
* SIGMS   scattering macroscopic cross sections in each mixture of the
*         composite geometry.
*
*Parameters: output
* NCO     number of volumes in each composite mixture.
* RRRR    information used by XDRH20, XDRH23, XDRH30 and XDRH33.
* QKOLD   information used by XDRH20, XDRH23, XDRH30 and XDRH33.
* QKDEL   information used by XDRH20, XDRH23, XDRH30 and XDRH33.
* PKL     information used by XDRH20, XDRH23, XDRH30 and XDRH33.
* COEF    information used by XDRH20, XDRH23, XDRH30 and XDRH33.
*
*References:
* R. Sanchez and G. C. Pomraning, A Statistical Analysis of the Double
* Heterogeneity Problem, Ann. Nucl. Energy, 18, 371-395 (1991).
* \\\\
* R. Sanchez and E. Masiello, Treatment of the Double Heterogeneity
* with the Method of Characteristics", PHYSOR 2002, Seoul, Korea (2002).
* \\\\
* R. Sanchez, Renormalized Treatment of the Double Heterogeneity with
* the Method of Characteristics, PHYSOR 2004, Chicago, USA (2004).
* 
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR1,NMILG,NG,NSMAX,MICRO,IQUAD,NS(NG),IDIL(NMILG),
     1 MIXGR(NSMAX,NG,NMILG),NCO(NMILG)
      REAL RS(NSMAX+1,NG),FRACT(NG,IR1+NMILG),VOLK(NG,NSMAX),
     1 SIGMA(0:IR1+NMILG),SIGMS(0:IR1+NMILG),RRRR(NMILG),
     2 QKOLD(NG,NSMAX,NMILG),QKDEL(NG,NSMAX,NMILG),
     3 PKL(NG,NSMAX,NSMAX,NMILG)
      DOUBLE PRECISION COEF(1+NG*NSMAX,1+NG*NSMAX,NMILG)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(EPS1=1.0E-5)
      DOUBLE PRECISION DP0,DP1,DP1OLD,QKD,CHORD,CHORDK,DDOT
      REAL, ALLOCATABLE, DIMENSION(:) :: SIG,ZZ
      REAL, ALLOCATABLE, DIMENSION(:,:) :: QKN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: RHS
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(QKN(NSMAX,NSMAX),RHS(1+NG*NSMAX))
*----
*  COMPUTE THE EQUIVALENT TOTAL CROSS SECTIONS IN COMPOSITE REGIONS
*----
      CALL XDRSET(QKOLD,NG*NSMAX*NMILG,0.0)
      CALL XDRSET(QKDEL,NG*NSMAX*NMILG,0.0)
      CALL XDRSET(PKL,NG*NSMAX*NSMAX*NMILG,0.0)
      IPOW=MICRO-1
      DO 110 IBM=1,NMILG
      MIL=IR1+IBM
      DILF=1.0
      DO 10 J=1,NG
      DILF=DILF-FRACT(J,MIL)
   10 CONTINUE
      DP1OLD=0.0D0
      ITER=0
      DO
         ITER=ITER+1
         IF(ITER.GT.100) CALL XABORT('XDRH11: CONVERGENCE FAILURE.')
         DP1=DILF*SIGMA(IDIL(IBM))
         DO 80 J=1,NG
         FRT=FRACT(J,MIL)
         IF(FRT.LE.0.00001) GO TO 80
         ALLOCATE(SIG(NS(J)))
         CHORD=0.0D0
         DO 20 K=1,NS(J)
         SIG(K)=REAL(SIGMA(MIXGR(K,J,IBM))-DP1OLD)
         SIG(K)=MAX(0.0,SIG(K))
         CHORD=CHORD+(RS(K+1,J)**IPOW-RS(K,J)**IPOW)*SIG(K)
   20    CONTINUE
         CHORD=4.0D0*CHORD/(REAL(IPOW)*RS(NS(J)+1,J)**(IPOW-1))
         ALLOCATE(ZZ(1+IQUAD*((NS(J)*(5+NS(J)))/2)))
         IF(MICRO.EQ.3) THEN
            CALL SYBT1D(NS(J),RS(1,J),.FALSE.,IQUAD,ZZ)
            CALL SYBALC(NS(J),NSMAX,RS(1,J),SIG,IQUAD,0.0,ZZ,QKN)
         ELSE IF(MICRO.EQ.4) THEN
            CALL SYBT1D(NS(J),RS(1,J),.TRUE.,IQUAD,ZZ)
            CALL SYBALS(NS(J),NSMAX,RS(1,J),SIG,IQUAD,0.0,ZZ,QKN)
         ENDIF
         DEALLOCATE(ZZ)
         IF(CHORD.GE.1.0E4) THEN
            DO 25 K=1,NS(J)
            CHORDK=4.0D0*(RS(K+1,J)**IPOW-RS(K,J)**IPOW)/(REAL(IPOW)
     1      *RS(NS(J)+1,J)**(IPOW-1))
            QKDEL(J,K,IBM)=REAL(CHORDK/CHORD)
            DP1=DP1+FRT*VOLK(J,K)*QKDEL(J,K,IBM)*SIG(K)
   25       CONTINUE
         ELSE
            DO 40 K=1,NS(J)
            QKD=1.0D0
            DO 30 N=1,NS(J)
            QKD=QKD-QKN(K,N)*SIG(N)
   30       CONTINUE
            QKDEL(J,K,IBM)=REAL(QKD)
            DP1=DP1+FRT*VOLK(J,K)*QKDEL(J,K,IBM)*SIG(K)
   40       CONTINUE
         ENDIF
         IF(ITER.EQ.1) THEN
            DO 50 K=1,NS(J)
            QKOLD(J,K,IBM)=QKDEL(J,K,IBM)
   50       CONTINUE
            IF(CHORD.GE.1.0E4) THEN
               DO 60 K=1,NS(J)
               PKL(J,K,K,IBM)=1.0/SIGMA(MIXGR(K,J,IBM))
   60          CONTINUE
            ELSE
               DO 75 K=1,NS(J)
               DO 70 N=1,NS(J)
               PKL(J,K,N,IBM)=QKN(K,N)
   70          CONTINUE
   75          CONTINUE
            ENDIF
         ENDIF
         DEALLOCATE(SIG)
   80    CONTINUE
         IF(ABS(DP1OLD-DP1/DILF).LE.EPS1*ABS(DP1)) EXIT
         DP1OLD=DP1/DILF
      ENDDO
      RRRR(IBM)=DILF
      SIGMIN=REAL(DP1)/DILF
      DO 100 J=1,NG
      FRT=FRACT(J,MIL)
      IF(FRT.LE.0.00001) GO TO 100
      DO 90 K=1,NS(J)
      RRRR(IBM)=RRRR(IBM)+FRT*VOLK(J,K)*QKDEL(J,K,IBM)
      SIGMIN=MIN(SIGMIN,SIGMA(MIXGR(K,J,IBM)))
   90 CONTINUE
  100 CONTINUE
      IF((SIGMIN*(1.0+EPS1).LT.DP1/DILF).AND.(MICRO.EQ.3)) THEN
         CALL XABORT('XDRH11: SANCHEZ-POMRANING MODEL FAILURE.')
      ENDIF
      SIGMA(IR1+IBM)=REAL(DP1)/DILF
  110 CONTINUE
*----
*  COMPUTE THE EQUIVALENT SCATTERING CROSS SECTIONS IN COMPOSITE REGIONS
*----
      CALL XDDSET(COEF,NMILG*(1+NG*NSMAX)**2,0.0D0)
      DO 170 IBM=1,NMILG
      MIL=IR1+IBM
      NCO(IBM)=1
      DILF=1.0
      DP0=0.0D0
      DO 130 J=1,NG
      FRT=FRACT(J,MIL)
      DILF=DILF-FRT
      IF(FRT.LE.0.00001) GO TO 130
      DO 120 K=1,NS(J)
      DP0=DP0+FRT*VOLK(J,K)*QKOLD(J,K,IBM)*SIGMA(MIXGR(K,J,IBM))
  120 CONTINUE
  130 CONTINUE
      DP0=DP0+DILF*SIGMA(IDIL(IBM))
      COEF(1,1,IBM)=1.0D0
      RHS(1)=DILF*SIGMS(IDIL(IBM))/DP0
      IND2=1
      DO 160 J=1,NG
      FRT=FRACT(J,MIL)
      IF(FRT.LE.0.00001) GO TO 160
      DO 150 K=1,NS(J)
      NCO(IBM)=NCO(IBM)+1
      COEF(1,IND2+K,IBM)=-FRT*VOLK(J,K)*QKOLD(J,K,IBM)*
     1 SIGMS(MIXGR(K,J,IBM))/DP0
      COEF(IND2+K,IND2+K,IBM)=1.0D0
      DO 140 N=1,NS(J)
      COEF(IND2+K,IND2+N,IBM)=COEF(IND2+K,IND2+N,IBM)
     1 -PKL(J,K,N,IBM)*SIGMS(MIXGR(N,J,IBM))
  140 CONTINUE
      COEF(IND2+K,1,IBM)=-(QKOLD(J,K,IBM)-QKDEL(J,K,IBM))
      RHS(IND2+K)=QKDEL(J,K,IBM)
  150 CONTINUE
      IND2=IND2+NS(J)
  160 CONTINUE
      CALL ALINVD(NCO(IBM),COEF(1,1,IBM),1+NG*NSMAX,IER)
      IF(IER.NE.0) CALL XABORT('XDRH11: SINGULAR MATRIX.')
      DP0=DDOT(NCO(IBM),COEF(1,1,IBM),1+NG*NSMAX,RHS,1)
      SIGMS(IR1+IBM)=REAL(DP0)*SIGMA(IR1+IBM)
  170 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(RHS,QKN)
      RETURN
      END
