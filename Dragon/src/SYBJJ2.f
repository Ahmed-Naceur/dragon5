*DECK SYBJJ2
      SUBROUTINE SYBJJ2 (IPAS,NMCEL,NMERGE,NGEN,IJAT,NPIJ,NPIS,EPSJ,
     1 NUNKNO,FUNKNO,SUNKNO,IMPX,NCOUR,NMC,IFR,ALB,INUM,MIX,DVX,IGEN,
     2 PIJW,PISW,PSJW,PSSW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the neutron flux and interface currents in a 2-D Cartesian
* or hexagonal assembly using the current iteration method with
* Roth X 4, DP0 or DP1 approximation.
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
* NMCEL   total number of cells in the domain.
* NMERGE  total number of merged cells for which specific values
*         of the neutron flux and reactions rates are required.
*         Many cells with different position in the domain can
*         be merged before the neutron flux calculation if they
*         own the same generating cell (NMERGE.le.NMCEL).
* NGEN    total number of generating cells. A generating cell is
*         defined by its material and dimensions, irrespective of
*         its position in the domain (NGEN.le.NMERGE).
* IJAT    total number of distinct out-currents.
* NPIJ    size of cellwise scattering-reduced collision probability
*         matrices.
* NPIS    size of cellwise scattering-reduced escape probability
*         matrices.
* EPSJ    stopping criterion for flux-current iterations.
* NUNKNO  total number of unknowns in vectors SUNKNO and FUNKNO.
* SUNKNO  input source vector.
* IMPX    print flag (equal to 0 for no print).
* NCOUR   number of incoming currents (=4 Cartesian lattice;
*         =6 hexagonal lattice).
* NMC     offset of the first volume in each generating cell.
* IFR     index-number of in-currents.
* ALB     transmission/albedo associated with each in-current.
* INUM    index-number of the merged cell associated to each cell.
* MIX     index-number of out-currents.
* DVX     weight associated with each out-current.
*         Note: IFR, ALB, MIX and DVX contains information to rebuild
*         the geometrical 'A' matrix.
* IGEN    index-number of the generating cell associated with each
*         merged cell.
* PIJW    cellwise scattering-reduced collision probability matrices.
* PISW    cellwise scattering-reduced escape probability matrices.
* PSJW    cellwise scattering-reduced collision probability matrices
*         for incoming neutrons.
* PSSW    cellwise scattering-reduced transmission probability
*         matrices.
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
      INTEGER IPAS,NMCEL,NMERGE,NGEN,IJAT,NPIJ,NPIS,NUNKNO,IMPX,NCOUR,
     1 NMC(NGEN+1),IFR(NCOUR*NMCEL),INUM(NMCEL),MIX(NCOUR*NMERGE),
     2 IGEN(NMERGE)
      REAL EPSJ,FUNKNO(NUNKNO),SUNKNO(NUNKNO),ALB(NCOUR*NMCEL),
     1 DVX(NCOUR*NMERGE),PIJW(NPIJ),PISW(NCOUR*NPIS),PSJW(NCOUR*NPIS),
     2 PSSW(NGEN*NCOUR*NCOUR)
*----
*  LOCAL VARIABLES
*----
      REAL PIJ,PIS
      LOGICAL LOGTES
      PARAMETER (MAXIT=400,LACCFC=2,ICL1=3,ICL2=3)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, DIMENSION(:), POINTER :: INDPIJ,INDNMC
      DOUBLE PRECISION, DIMENSION(:), POINTER :: CIT0
      DOUBLE PRECISION, DIMENSION(:,:), POINTER :: CITR,AITR
      DOUBLE PRECISION, DIMENSION(:), POINTER :: WCURR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDPIJ(NGEN),INDNMC(NMERGE))
      ALLOCATE(CITR(3,IJAT),CIT0(IJAT),AITR(2,IJAT))
      ALLOCATE(WCURR(IJAT))
*
      IPIJ=0
      DO 10 JKG=1,NGEN
      J2=NMC(JKG+1)-NMC(JKG)
      INDPIJ(JKG)=IPIJ
      IPIJ=IPIJ+J2*J2
   10 CONTINUE
      KNMC=0
      DO 20 JKK=1,NMERGE
      JKG=IGEN(JKK)
      J2=NMC(JKG+1)-NMC(JKG)
      INDNMC(JKK)=KNMC
      KNMC=KNMC+J2
   20 CONTINUE
*
      DO 30 I=1,IJAT
      WCURR(I)=1.0D0
      CIT0(I)=0.0D0
      CITR(1,I)=FUNKNO(IPAS+I)
   30 CONTINUE
*----
*  COMPUTE PSJW * Q(*) CONTRIBUTION
*----
      DO 42 IKK=1,NMERGE
      IKG=IGEN(IKK)
      I1P=NMC(IKG)
      I2=NMC(IKG+1)-I1P
      IT=NCOUR*(IKK-1)
      KNMC=INDNMC(IKK)
      DO 41 I=1,I2
      DO 40 IC=1,NCOUR
      JCC=MIX(IT+IC)
      PBJ=PSJW(I1P*NCOUR+(I-1)*NCOUR+IC)
      CIT0(JCC)=CIT0(JCC)+PBJ*DVX(IT+IC)*SUNKNO(KNMC+I)
   40 CONTINUE
   41 CONTINUE
   42 CONTINUE
*----
*  COMPUTE NORMALIZATION VECTOR WCURR
*----
      DO 52 ICEL=1,NMCEL
      IKK=INUM(ICEL)
      IT=NCOUR*(IKK-1)
      IS=NCOUR*(ICEL-1)
      IKG=IGEN(IKK)
      IPSS=(IKG-1)*NCOUR*NCOUR
      DO 51 JC=1,NCOUR
      J1=IFR(IS+JC)
      DO 50 IC=1,NCOUR
      J2=MIX(IT+IC)
      PSS=PSSW(IPSS+(JC-1)*NCOUR+IC)
      WCURR(J1)=WCURR(J1)-PSS*ALB(IS+JC)*DVX(IT+IC)
   50 CONTINUE
   51 CONTINUE
   52 CONTINUE
*
      ISTART=1
      TEST=0.0D0
      ITER=0
   70 ITER=ITER+1
      IF(ITER.GT.MAXIT) THEN
         WRITE(6,'(/47H SYBJJ2: *** WARNING *** MAXIMUM NUMBER OF ITER,
     1   15HATIONS REACHED.)')
         GO TO 190
      ENDIF
      IT3=MOD(ITER,3)+1
      IT2=MOD(ITER-1,3)+1
      IT1=MOD(ITER-2,3)+1
      DO 80 I=1,IJAT
      CITR(IT3,I)=CIT0(I)
   80 CONTINUE
*----
*  COMPUTE PSSW * J(-) CONTRIBUTION
*----
      DO 92 ICEL=1,NMCEL
      IKK=INUM(ICEL)
      IT=NCOUR*(IKK-1)
      IS=NCOUR*(ICEL-1)
      IKG=IGEN(IKK)
      IPSS=(IKG-1)*NCOUR*NCOUR
      DO 91 JC=1,NCOUR
      J1=IFR(IS+JC)
      DO 90 IC=1,NCOUR
      J2=MIX(IT+IC)
      PSS=PSSW(IPSS+(JC-1)*NCOUR+IC)
      CITR(IT3,J2)=CITR(IT3,J2)+PSS*ALB(IS+JC)*DVX(IT+IC)*CITR(IT2,J1)
   90 CONTINUE
   91 CONTINUE
   92 CONTINUE
*----
*  NORMALIZATION
*----
      S1=0.0D0
      S2=0.0D0
      DO 100 I=1,IJAT
      S1=S1+WCURR(I)*CITR(IT3,I)
      S2=S2+CIT0(I)
  100 CONTINUE
      ZNORM=S2/S1
      IF(ZNORM.LT.0.0D0) ZNORM=1.0D0
      DO 110 I=1,IJAT
      CITR(IT3,I)=CITR(IT3,I)*ZNORM
  110 CONTINUE
*----
*  ONE/TWO PARAMETER ACCELERATION
*----
      ALP=1.0D0
      BET=0.0D0
      LOGTES=(1+MOD(ITER-ISTART,ICL1+ICL2).GT.ICL1)
      IF(LOGTES) THEN
         DO 120 I=1,IJAT
         AITR(1,I)=CITR(IT3,I)-CITR(IT2,I)
         AITR(2,I)=CITR(IT2,I)-CITR(IT1,I)
  120    CONTINUE
         DO 132 ICEL=1,NMCEL
         IKK=INUM(ICEL)
         IT=NCOUR*(IKK-1)
         IS=NCOUR*(ICEL-1)
         IKG=IGEN(IKK)
         IPSS=(IKG-1)*NCOUR*NCOUR
         DO 131 JC=1,NCOUR
         J1=IFR(IS+JC)
         DO 130 IC=1,NCOUR
         J2=MIX(IT+IC)
         PSS=PSSW(IPSS+(JC-1)*NCOUR+IC)*ALB(IS+JC)*DVX(IT+IC)
         AITR(1,J2)=AITR(1,J2)-PSS*(CITR(IT3,J1)-CITR(IT2,J1))
         AITR(2,J2)=AITR(2,J2)-PSS*(CITR(IT2,J1)-CITR(IT1,J1))
  130    CONTINUE
  131    CONTINUE
  132    CONTINUE
         IF((LACCFC.EQ.1).OR.(MOD(ITER-ISTART,ICL1+ICL2).EQ.ICL1)) THEN
            S1=0.0D0
            S2=0.0D0
            DO 140 I=1,IJAT
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
            DO 150 I=1,IJAT
            CITR(IT3,I)=CITR(IT2,I)+ALP*(CITR(IT3,I)-CITR(IT2,I))
  150       CONTINUE
         ELSE IF(LACCFC.EQ.2) THEN
            S1=0.0D0
            S2=0.0D0
            S3=0.0D0
            S4=0.0D0
            S5=0.0D0
            DO 160 I=1,IJAT
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
            DO 170 I=1,IJAT
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
      DO 180 I=1,IJAT
      ERR1=MAX(ERR1,ABS(CITR(IT3,I)-CITR(IT2,I)))
      ERR2=MAX(ERR2,ABS(CITR(IT3,I)))
  180 CONTINUE
      IF(IMPX.GT.3) WRITE(6,'(30H SYBJJ2: CURRENT ITERATION NB.,I4,
     1 7H ERROR=,1P,E10.3,5H OVER,E10.3,15H NORMALIZATION=,E10.3,
     2 14H ACCELERATION=,2E11.3,1H.)') ITER,ERR1,ERR2,ZNORM,ALP,
     3 BET/ALP
      IF(ITER.EQ.1) TEST=ERR1/ERR2
      IF((ITER.GT.20).AND.(ERR1/ERR2.GT.TEST)) THEN
         WRITE(6,'(/45H SYBJJ2: *** WARNING *** CONVERGENCE DIFFICUL,
     1   5HTIES.)')
         GO TO 190
      ENDIF
      IF(LOGTES.OR.(ERR1.GT.EPSJ*ERR2)) GO TO 70
      IF(IMPX.GT.2) WRITE(6,'(37H SYBJJ2: CURRENT CONVERGENCE AT ITERA,
     1 8HTION NB.,I4,7H ERROR=,1P,E10.3,5H OVER,E10.3,1H.)') ITER,ERR1,
     2 ERR2
*
  190 DO 200 I=1,IPAS
      FUNKNO(I)=0.0
  200 CONTINUE
      DO 210 I=1,IJAT
      FUNKNO(IPAS+I)=REAL(CITR(IT3,I))
  210 CONTINUE
*----
*  COMPUTE PISW * J(-) CONTRIBUTION
*----
      DO 240 ICEL=1,NMCEL
      IKK=INUM(ICEL)
      IS=NCOUR*(ICEL-1)
      IKG=IGEN(IKK)
      I1P=NMC(IKG)
      I2=NMC(IKG+1)-I1P
      KNMC=INDNMC(IKK)
      DO 230 J=1,I2
      DO 220 JC=1,NCOUR
      J1=IFR(IS+JC)
      PIS=PISW(I1P*NCOUR+(JC-1)*I2+J)
      FUNKNO(KNMC+J)=FUNKNO(KNMC+J)+PIS*ALB(IS+JC)*FUNKNO(IPAS+J1)
  220 CONTINUE
  230 CONTINUE
  240 CONTINUE
*----
*  COMPUTE PIJW * Q(*) CONTRIBUTION
*----
      DO 270 IKK=1,NMERGE
      IKG=IGEN(IKK)
      I2=NMC(IKG+1)-NMC(IKG)
      KNMC=INDNMC(IKK)
      DO 260 I=1,I2
      DO 250 J=1,I2
      PIJ=PIJW(INDPIJ(IKG)+(I-1)*I2+J)
      FUNKNO(KNMC+J)=FUNKNO(KNMC+J)+PIJ*SUNKNO(KNMC+I)
  250 CONTINUE
  260 CONTINUE
  270 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WCURR)
      DEALLOCATE(AITR,CIT0,CITR)
      DEALLOCATE(INDNMC,INDPIJ)
      RETURN
      END
