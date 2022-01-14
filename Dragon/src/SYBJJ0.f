*DECK SYBJJ0
      SUBROUTINE SYBJJ0 (IPAS,NSUPCE,NPIJ,NUNKNO,EPSJ,FUNKNO,SUNKNO,
     1 IMPX,ISTAT,NMC,PROCEL,PIJW,PISW,PSJW,PSSW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the neutron flux and interface currents in a do-it-yourself
* geometry using the current iteration method.
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
* IPAS    total number of regions.
* NSUPCE  number of cells.
* NPIJ    length of cellwise scattering-reduced collision probability
*         matrices.
* EPSJ    stopping criterion for flux-current iterations.
* NUNKNO  total number of unknowns in vectors SUNKNO and FUNKNO.
* SUNKNO  input source vector.
* IMPX    print flag (equal to 0 for no print).
* ISTAT   statistical approximation flag (set with ISTAT=1).
* NMC     offset of the first volume in each cell.
* PROCEL  user supplied geometrical matrix.
* PIJW    cellwise scattering-reduced collision probability matrices.
* PISW    cellwise scattering-reduced escape probability matrices.
* PSJW    cellwise scattering-reduced collision probability matrices
*         for incoming neutrons.
* PSSW    cellwise scattering-reduced transmission probability matrices.
*
*Parameters: input/output
* FUNKNO  unknown vector.
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPAS,NSUPCE,NPIJ,NUNKNO,IMPX,ISTAT,NMC(NSUPCE+1)
      REAL EPSJ,FUNKNO(NUNKNO),SUNKNO(NUNKNO),PROCEL(NSUPCE,NSUPCE),
     1 PIJW(NPIJ),PISW(IPAS),PSJW(IPAS),PSSW(NSUPCE)
*----
*  LOCAL VARIABLES
*----
      REAL PIJ
      LOGICAL LOGTES
      PARAMETER (MAXIT=400,LACCFC=2,ICL1=3,ICL2=3)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, DIMENSION(:), POINTER :: INDPIJ
      DOUBLE PRECISION, DIMENSION(:), POINTER :: CIT0
      DOUBLE PRECISION, DIMENSION(:,:), POINTER :: CITR,AITR
      DOUBLE PRECISION, DIMENSION(:), POINTER :: WCURR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDPIJ(NSUPCE))
      ALLOCATE(CITR(3,NSUPCE),CIT0(NSUPCE),AITR(2,NSUPCE))
      ALLOCATE(WCURR(NSUPCE))
*
      IPIJ=0
      DO 10 JKG=1,NSUPCE
      J2=NMC(JKG+1)-NMC(JKG)
      INDPIJ(JKG)=IPIJ
      IPIJ=IPIJ+J2*J2
   10 CONTINUE
*----
*  PROCESS STATISTICAL APPROXIMATION
*----
      IF(ISTAT.NE.0) THEN
         X1=0.0D0
         DO 20 IKK=1,NSUPCE
         X1=X1+PSSW(IKK)*PROCEL(1,IKK)
   20    CONTINUE
         X1=1.0D0/(1.0D0-X1)
         SSS=0.0D0
         DO 35 IKK=1,NSUPCE
         I1P=NMC(IKK)
         I2=NMC(IKK+1)-I1P
         DO 30 I=1,I2
         SSS=SSS+PROCEL(1,IKK)*X1*PSJW(I1P+I)*SUNKNO(I1P+I)
   30    CONTINUE
   35    CONTINUE
         IT3=1
         DO 40 IKK=1,NSUPCE
         CITR(IT3,IKK)=SSS
   40    CONTINUE
         GO TO 190
      ENDIF
*----
*  COMPUTE PSJW * Q(*) CONTRIBUTION
*----
      DO 52 IKK=1,NSUPCE
      CIT0(IKK)=0.0D0
      CITR(1,IKK)=FUNKNO(IPAS+IKK)
      DO 51 JKK=1,NSUPCE
      I1P=NMC(JKK)
      I2=NMC(JKK+1)-I1P
      DO 50 I=1,I2
      CIT0(IKK)=CIT0(IKK)+PROCEL(IKK,JKK)*PSJW(I1P+I)*SUNKNO(I1P+I)
   50 CONTINUE
   51 CONTINUE
   52 CONTINUE
*----
*  COMPUTE NORMALIZATION VECTOR WCURR
*----
      DO 65 JKK=1,NSUPCE
      WCURR(JKK)=1.0D0
      DO 60 IKK=1,NSUPCE
      WCURR(JKK)=WCURR(JKK)-PROCEL(IKK,JKK)*PSSW(JKK)
   60 CONTINUE
   65 CONTINUE
*
      ISTART=1
      TEST=0.0D0
      ITER=0
   70 ITER=ITER+1
      IF(ITER.GT.MAXIT) THEN
         WRITE(6,'(/47H SYBJJ0: *** WARNING *** MAXIMUM NUMBER OF ITER,
     1   15HATIONS REACHED.)')
         GO TO 190
      ENDIF
      IT3=MOD(ITER,3)+1
      IT2=MOD(ITER-1,3)+1
      IT1=MOD(ITER-2,3)+1
      DO 80 I=1,NSUPCE
      CITR(IT3,I)=CIT0(I)
   80 CONTINUE
*----
*  COMPUTE PSSW * J(-) CONTRIBUTION
*----
      DO 95 IKK=1,NSUPCE
      DO 90 JKK=1,NSUPCE
      PSS=PROCEL(IKK,JKK)*PSSW(JKK)
      CITR(IT3,IKK)=CITR(IT3,IKK)+PSS*CITR(IT2,JKK)
   90 CONTINUE
   95 CONTINUE
*----
*  NORMALIZATION
*----
      S1=0.0D0
      S2=0.0D0
      DO 100 I=1,NSUPCE
      S1=S1+WCURR(I)*CITR(IT3,I)
      S2=S2+CIT0(I)
  100 CONTINUE
      ZNORM=S2/S1
      IF(ZNORM.LT.0.0D0) ZNORM=1.0D0
      DO 110 I=1,NSUPCE
      CITR(IT3,I)=CITR(IT3,I)*ZNORM
  110 CONTINUE
*----
*  ONE/TWO PARAMETER ACCELERATION
*----
      ALP=1.0D0
      BET=0.0D0
      LOGTES=(1+MOD(ITER-ISTART,ICL1+ICL2).GT.ICL1)
      IF(LOGTES) THEN
         DO 135 IKK=1,NSUPCE
         AITR(1,IKK)=CITR(IT3,IKK)-CITR(IT2,IKK)
         AITR(2,IKK)=CITR(IT2,IKK)-CITR(IT1,IKK)
         DO 130 JKK=1,NSUPCE
         PSS=PROCEL(IKK,JKK)*PSSW(JKK)
         AITR(1,IKK)=AITR(1,IKK)-PSS*(CITR(IT3,JKK)-CITR(IT2,JKK))
         AITR(2,IKK)=AITR(2,IKK)-PSS*(CITR(IT2,JKK)-CITR(IT1,JKK))
  130    CONTINUE
  135    CONTINUE
         IF((LACCFC.EQ.1).OR.(MOD(ITER-ISTART,ICL1+ICL2).EQ.ICL1)) THEN
            S1=0.0D0
            S2=0.0D0
            DO 140 I=1,NSUPCE
            S1=S1+(CITR(IT3,I)-CITR(IT2,I))*AITR(1,I)
            S2=S2+AITR(1,I)*AITR(1,I)
  140       CONTINUE
            IF(S2.EQ.0.0D0) THEN
               ISTART=ITER+1
            ELSE
               ALP=S1/S2
               IF(ALP.LE.0.0D0) THEN
                  ISTART=ITER+1
                  ALP=1.0D0
               ENDIF
            ENDIF
            DO 150 I=1,NSUPCE
            CITR(IT3,I)=CITR(IT2,I)+ALP*(CITR(IT3,I)-CITR(IT2,I))
  150       CONTINUE
         ELSE IF(LACCFC.EQ.2) THEN
            S1=0.0D0
            S2=0.0D0
            S3=0.0D0
            S4=0.0D0
            S5=0.0D0
            DO 160 I=1,NSUPCE
            S1=S1+(CITR(IT3,I)-CITR(IT2,I))*AITR(1,I)
            S2=S2+AITR(1,I)*AITR(1,I)
            S3=S3+(CITR(IT3,I)-CITR(IT2,I))*AITR(2,I)
            S4=S4+AITR(1,I)*AITR(2,I)
            S5=S5+AITR(2,I)*AITR(2,I)
  160       CONTINUE
            DET=S2*S5-S4*S4
            IF(DET.EQ.0.0D0) THEN
               ISTART=ITER+1
            ELSE
               ALP=(S5*S1-S4*S3)/DET
               BET=(S2*S3-S4*S1)/DET
               IF(ALP.LE.0.0D0) THEN
                  ISTART=ITER+1
                  ALP=1.0D0
                  BET=0.0D0
               ENDIF
            ENDIF
            DO 170 I=1,NSUPCE
            CITR(IT3,I)=CITR(IT2,I)+ALP*(CITR(IT3,I)-CITR(IT2,I))+
     1      BET*(CITR(IT2,I)-CITR(IT1,I))
  170       CONTINUE
         ENDIF
      ENDIF
*----
*  CHECK THE CONVERGENCE ERROR
*----
      ERR1=0.0D0
      ERR2=0.0D0
      DO 180 I=1,NSUPCE
      ERR1=MAX(ERR1,ABS(CITR(IT3,I)-CITR(IT2,I)))
      ERR2=MAX(ERR2,ABS(CITR(IT3,I)))
  180 CONTINUE
      IF(IMPX.GT.3) WRITE(6,'(30H SYBJJ0: CURRENT ITERATION NB.,I4,
     1 7H ERROR=,1P,E10.3,5H OVER,E10.3,15H NORMALIZATION=,E10.3,
     2 14H ACCELERATION=,2E11.3,1H.)') ITER,ERR1,ERR2,ZNORM,ALP,
     3 BET/ALP
      IF(ITER.EQ.1) TEST=ERR1/ERR2
      IF((ITER.GT.20).AND.(ERR1/ERR2.GT.TEST)) CALL XABORT('SYBJJ0: '
     1 //'CONVERGENCE FAILURE.')
      IF(LOGTES.OR.(ERR1.GT.EPSJ*ERR2)) GO TO 70
      IF(IMPX.GT.2) WRITE(6,'(37H SYBJJ0: CURRENT CONVERGENCE AT ITERA,
     1 8HTION NB.,I4,7H ERROR=,1P,E10.3,5H OVER,E10.3,1H.)') ITER,ERR1,
     2 ERR2
*
  190 DO 200 I=1,IPAS
      FUNKNO(I)=0.0
  200 CONTINUE
      DO 210 I=1,NSUPCE
      FUNKNO(IPAS+I)=REAL(CITR(IT3,I))
  210 CONTINUE
*----
*  COMPUTE ( PISW * J(-) ) + ( PIJW * Q(*) ) CONTRIBUTION
*----
      DO 240 IKK=1,NSUPCE
      I1P=NMC(IKK)
      I2=NMC(IKK+1)-I1P
      DO 230 J=1,I2
      FUNKNO(I1P+J)=FUNKNO(I1P+J)+PISW(I1P+J)*FUNKNO(IPAS+IKK)
      DO 220 I=1,I2
      PIJ=PIJW(INDPIJ(IKK)+(I-1)*I2+J)
      FUNKNO(I1P+J)=FUNKNO(I1P+J)+PIJ*SUNKNO(I1P+I)
  220 CONTINUE
  230 CONTINUE
  240 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WCURR)
      DEALLOCATE(AITR,CIT0,CITR)
      DEALLOCATE(INDPIJ)
      RETURN
      END
