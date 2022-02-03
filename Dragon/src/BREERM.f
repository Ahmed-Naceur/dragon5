*DECK BREERM
      SUBROUTINE BREERM(IPMAC1,NC,NG,LX1,NMIX1,ITRIAL,IMIX,ICODE,ISPH,
     1 ZKEFF,B2,ENER,VOL1,FLX1,DC1,TOT1,CHI1,SIGF1,SCAT1,JXM,JXP,FHETXM,
     2 FHETXP,ADF1,NGET,ADFREF,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Implement the 1D ERM-NEM reflector model.
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
* IPMAC1  nodal macrolib.
* NC      number of sn macrolibs.
* NG      number of energy groups.
* LX1     number of nodes in the reflector model.
* NMIX1   number of mixtures in the nodal calculation.
* ITRIAL  type of expansion functions in the nodal calculation.
*         (=1: polynomial; =2: hyperbolic).
* IMIX    mix index of each node.
* ICODE   physical albedo index on each side of the domain.
* ISPH    SPH flag (=0: use discontinuity factors; =1: use SPH factors).
* ZKEFF   effective multiplication factor.
* B2      buckling.
* ENER    energy limits.
* VOL1    volumes.
* FLX1    averaged fluxes
* DC1     diffusion coefficients.
* TOT1    total cross sections.
* CHI1    fission spectra.
* SIGF1   nu*fission cross sections.
* SCAT1   scattering P0 cross sections.
* JXM     left boundary currents.
* JXP     right boundary currents.
* FHETXM  left boundary fluxes.
* FHETXP  right boundary fluxes.
* ADF1    assembly discontinuity factors from macrolib.
* NGET    type of NGET normalization if discontinuity factors
*         (=0: simple; =1: imposed ADF on fuel assembly; =2: recover
*         fuel assembly ADF from input macrolib).
* ADFREF  imposed ADF values on fuel assembly side.
* IPRINT  edition flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC1
      INTEGER NC,NG,LX1,NMIX1,ITRIAL(NG),IMIX(LX1),ICODE(2),ISPH,NGET,
     1 IPRINT
      REAL ZKEFF(NC),B2(NC),ENER(NG+1),VOL1(NMIX1,NC),FLX1(NMIX1,NG,NC),
     1 DC1(NMIX1,NG,NC),TOT1(NMIX1,NG,NC),CHI1(NMIX1,NG,NC),
     2 SIGF1(NMIX1,NG,NC),SCAT1(NMIX1,NG,NG,NC),JXM(NMIX1,NG,NC),
     3 JXP(NMIX1,NG,NC),FHETXM(NMIX1,NG,NC),FHETXP(NMIX1,NG,NC),
     4 ADF1(NMIX1,NG,NC),ADFREF(NG)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE)
      REAL SX(5),AA11(4,4),A11(5,5),Q(5)
      CHARACTER(LEN=8) HADF(2)
      LOGICAL LFISS
      TYPE(C_PTR) JPMAC1,KPMAC1
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK,ETA,VOL
      REAL, ALLOCATABLE, DIMENSION(:,:) :: AB,ALPHA,FLX,DC,TOT,CHI,SIGF,
     1 ADF,AFACTOR,BETA
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: FDXM,FDXP,SCAT
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: TAU,B,X
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: WORK2
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: FHOMM,FHOMP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ETA(NG),ALPHA(5,NG),FDXM(NMIX1,NG,NG),FDXP(NMIX1,NG,NG),
     1 AFACTOR(NG,NG),BETA(NG,NG),FHOMM(NC,NG,NMIX1),FHOMP(NC,NG,NMIX1))
      ALLOCATE(VOL(NMIX1),FLX(NMIX1,NG),DC(NMIX1,NG),TOT(NMIX1,NG),
     1 CHI(NMIX1,NG),SIGF(NMIX1,NG),SCAT(NMIX1,NG,NG),ADF(NMIX1,NG))
*----
*  AVERAGE THE OUTPUT NODAL MACROLIB
*----
      VOL(:)=0.0
      FLX(:,:)=0.0
      DC(:,:)=0.0
      TOT(:,:)=0.0
      CHI(:,:)=0.0
      SIGF(:,:)=0.0
      SCAT(:,:,:)=0.0
      ADF(:,:)=0.0
      DO IC=1,NC
        DO IBM=1,NMIX1
          VOL(IBM)=VOL(IBM)+VOL1(IBM,IC)
          DO IGR=1,NG
            FLX(IBM,IGR)=FLX(IBM,IGR)+FLX1(IBM,IGR,IC)
            DC(IBM,IGR)=DC(IBM,IGR)+DC1(IBM,IGR,IC)
            TOT(IBM,IGR)=TOT(IBM,IGR)+TOT1(IBM,IGR,IC)
            CHI(IBM,IGR)=CHI(IBM,IGR)+CHI1(IBM,IGR,IC)
            SIGF(IBM,IGR)=SIGF(IBM,IGR)+SIGF1(IBM,IGR,IC)
            DO JGR=1,NG
              SCAT(IBM,IGR,JGR)=SCAT(IBM,IGR,JGR)+SCAT1(IBM,IGR,JGR,IC)
            ENDDO
            ADF(IBM,IGR)=ADF(IBM,IGR)+ADF1(IBM,IGR,IC)
          ENDDO
        ENDDO
      ENDDO
      VOL(:)=VOL(:)/REAL(NC)
      FLX(:,:)=FLX(:,:)/REAL(NC)
      DC(:,:)=DC(:,:)/REAL(NC)
      TOT(:,:)=TOT(:,:)/REAL(NC)
      CHI(:,:)=CHI(:,:)/REAL(NC)
      SIGF(:,:)=SIGF(:,:)/REAL(NC)
      SCAT(:,:,:)=SCAT(:,:,:)/REAL(NC)
      ADF(:,:)=ADF(:,:)/REAL(NC)
*----
*  LOOP OVER CASES
*----
      IF(ISPH.EQ.1) CALL XABORT('BREERM: SPH OPTION NOT IMPLEMENTED.')
      J_FUEL=0
      DO IC=1,NC
*----
*  SET AND SOLVE NODAL SYSTEM
*----
        J_FUEL=0
        DO J=1,LX1
          IBM=IMIX(J)
          IF(IBM.EQ.0) CYCLE
          LFISS=.FALSE.
          DO IGR=1,NG
            IF(SIGF(IBM,IGR).GT.0.0) LFISS=.TRUE.
          ENDDO
          IF(LFISS) THEN
            J_FUEL=J
            ALLOCATE(AB(4*NG,4*NG+1))
            DO IGR=1,NG
              DIFF=DC1(IBM,IGR,IC)
              SIGR=TOT1(IBM,IGR,IC)+B2(IC)*DIFF-SCAT1(IBM,IGR,IGR,IC)
              ETA(IGR)=VOL1(IBM,IC)*SQRT(SIGR/DIFF)
              DO JGR=1,NG
                IF(JGR.EQ.IGR) THEN
                  SIGR=SIGR-CHI1(IBM,IGR,IC)*SIGF1(IBM,IGR,IC)/ZKEFF(IC)
                  CALL BRESS1(ITRIAL(IGR),VOL1(IBM,IC),DIFF,SIGR,AA11)
                ELSE
                  SIGT=-SCAT1(IBM,JGR,IGR,IC)-CHI1(IBM,JGR,IC)*
     1            SIGF1(IBM,IGR,IC)/ZKEFF(IC)
                  CALL BRESS2(ITRIAL(IGR),VOL1(IBM,IC),DIFF,SIGR,SIGT,
     1            AA11)
                ENDIF
                DO K1=1,4
                  DO K2=1,4
                    AB((JGR-1)*4+K1,(IGR-1)*4+K2)=AA11(K1,K2)
                  ENDDO
                ENDDO
              ENDDO
              SX = (/0.0,0.0,0.0,JXM(IBM,IGR,IC),JXP(IBM,IGR,IC)/) 
              DO K1=1,4
                AB((IGR-1)*4+K1,4*NG+1)=SX(K1+1)
              ENDDO
            ENDDO
            CALL ALSB(4*NG,1,AB,IER,4*NG)
            IF (IER.NE.0) CALL XABORT('BREERM: ALBS FAILURE(1).')
            DO IGR=1,NG
              ALPHA(1,IGR)=FLX1(IBM,IGR,IC)
              DO I=1,4
                ALPHA(I+1,IGR)=AB((IGR-1)*4+I,4*NG+1)
              ENDDO
            ENDDO
            DEALLOCATE(AB)
          ELSE
*           use averaged cross section values for the reflector
            ALLOCATE(AB(5*NG,5*NG+1))
            DO IGR=1,NG
              DIFF=DC(IBM,IGR)
              SIGR=TOT(IBM,IGR)+B2(IC)*DIFF-SCAT(IBM,IGR,IGR)
              ETA(IGR)=VOL(IBM)*SQRT(SIGR/DIFF)
              DO JGR=1,NG
                IF(JGR.EQ.IGR) THEN
                  CALL BRESS3(ITRIAL(IGR),VOL(IBM),DIFF,SIGR,A11)
                ELSE
                  SIGT=-SCAT(IBM,JGR,IGR)
                  CALL BRESS4(ITRIAL(IGR),VOL(IBM),DIFF,SIGR,SIGT,
     1            A11)
                ENDIF
                DO K1=1,5
                  DO K2=1,5
                    AB((JGR-1)*5+K1,(IGR-1)*5+K2)=A11(K1,K2)
                  ENDDO
                ENDDO
              ENDDO
              SX = (/0.0,0.0,0.0,JXM(IBM,IGR,IC),JXP(IBM,IGR,IC)/) 
              DO K1=1,5
                AB((IGR-1)*5+K1,5*NG+1)=SX(K1)
              ENDDO
            ENDDO
            CALL ALSB(5*NG,1,AB,IER,5*NG)
            IF (IER.NE.0) CALL XABORT('BREERM: ALBS FAILURE(2).')
            DO IGR=1,NG
              DO I=1,5
                ALPHA(I,IGR)=AB((IGR-1)*5+I,5*NG+1)
              ENDDO
            ENDDO
            DEALLOCATE(AB)
          ENDIF
          IF(IPRINT.GT.1) THEN
            WRITE(6,'(/9H MIXTURE=,I5,6H CASE=,I3)') J,IC
            WRITE(6,20) 'ALPHA',ALPHA(:5,:NG)
          ENDIF
*----
*  COMPUTE NODAL SURFACE FLUXES
*----
          DO IGR=1,NG
            IF (ITRIAL(IGR) == 1) THEN
              Q(1) = ALPHA(3,IGR)/2.
              FHOMM(IC,IGR,IBM)=-ALPHA(2,IGR)/2.+ALPHA(1,IGR)+Q(1)
              FHOMP(IC,IGR,IBM)=ALPHA(2,IGR)/2.+ALPHA(1,IGR)+Q(1)
            ELSE
              Q(1) = ETA(IGR)/2.
              Q(2) = SINH(Q(1))
              Q(3) = ALPHA(3,IGR)/2.
              Q(4) = ALPHA(4,IGR)*Q(2)
              Q(5) = ALPHA(5,IGR)*(COSH(Q(1)) - (2*Q(2))/ETA(IGR))
              FHOMM(IC,IGR,IBM)=-ALPHA(2,IGR)/2.+ALPHA(1,IGR)+Q(3)-Q(4)+
     1        Q(5)
              FHOMP(IC,IGR,IBM)=ALPHA(2,IGR)/2.+ALPHA(1,IGR)+Q(3)+Q(4)+
     1        Q(5)
            ENDIF
          ENDDO
        ENDDO
        IF(IPRINT.GT.0) THEN
          WRITE(6,'(/39H BREERM: NODAL SURFACE FLUXES FOR CASE=,I5)') IC
          DO IBM=1,NMIX1
            WRITE(6,'(/9H MIXTURE=,I5)') IBM
            WRITE(6,20) 'FHOMM',FHOMM(IC,:NG,IBM)
            WRITE(6,20) 'FHOMP',FHOMP(IC,:NG,IBM)
          ENDDO
        ENDIF
*----
*  END OF LOOP OVER CASES
*----
      ENDDO
*----
*  COMPUTE DISCONTINUITY AND ALBEDO FACTORS
*----
      AFACTOR(:,:)=0.0
      DO IBM=1,NMIX1
        IF(NC.EQ.1) THEN
          ! DF-NEM approach
          FDXM(IBM,:,:)=0.0
          FDXP(IBM,:,:)=0.0
          DO IGR=1,NG
            FDXM(IBM,IGR,IGR)=FHETXM(IBM,IGR,1)/REAL(FHOMM(1,IGR,IBM))
            FDXP(IBM,IGR,IGR)=FHETXP(IBM,IGR,1)/REAL(FHOMP(1,IGR,IBM))
          ENDDO
          IF(IBM.EQ.NMIX1) THEN
            DO IGR=1,NG
              AFACTOR(IGR,IGR)=JXP(IBM,IGR,1)/REAL(FHOMP(1,IGR,IBM))
            ENDDO
          ENDIF
        ELSE IF(NC.LT.NG) THEN
          CALL XABORT('BREERM: DEGENERATE SYSTEM')
        ELSE IF(NC.EQ.NG) THEN
          ! ERM-NEM approach: linear system resolution
          ALLOCATE(WORK2(NC,2*NG))
          DO IGR=1,NG
            DO IC=1,NC
              WORK2(IC,IGR)=FHOMM(IC,IGR,IBM)
              WORK2(IC,NG+IGR)=FHETXM(IBM,IGR,IC)
            ENDDO
          ENDDO
          CALL ALSBD(NC,NG,WORK2,IER,NC)
          IF(IER.NE.0) CALL XABORT('BREERM: SINGULAR MATRIX(1).')
          DO IGR=1,NG
            DO IC=1,NC
              FDXM(IBM,IGR,IC)=REAL(WORK2(IC,NG+IGR))
            ENDDO
          ENDDO
          DO IGR=1,NG
            DO IC=1,NC
              WORK2(IC,IGR)=FHOMP(IC,IGR,IBM)
              WORK2(IC,NG+IGR)=FHETXP(IBM,IGR,IC)
            ENDDO
          ENDDO
          CALL ALSBD(NC,NG,WORK2,IER,NC)
          IF(IER.NE.0) CALL XABORT('BREERM: SINGULAR MATRIX(2).')
          DO IGR=1,NG
            DO IC=1,NC
              FDXP(IBM,IGR,IC)=REAL(WORK2(IC,NG+IGR))
            ENDDO
          ENDDO
          IF(IBM.EQ.NMIX1) THEN
            DO IGR=1,NG
              DO IC=1,NC
                WORK2(IC,IGR)=FHOMP(IC,IGR,IBM)
                WORK2(IC,NG+IGR)=JXP(IBM,IGR,IC)
              ENDDO
            ENDDO
            CALL ALSBD(NC,NG,WORK2,IER,NC)
            IF(IER.NE.0) CALL XABORT('BREERM: SINGULAR MATRIX(3).')
            DO IGR=1,NG
              DO JGR=1,NG
                AFACTOR(IGR,JGR)=REAL(WORK2(JGR,NG+IGR))
              ENDDO
            ENDDO
          ENDIF
          DEALLOCATE(WORK2)
        ELSE IF(NC.GE.NG) THEN
          ! ERM-NEM approach: pseudo inversion
          ALLOCATE(TAU(NG),B(NC),X(NG))
          CALL ALST2F(NC,NC,NG,FHOMM(1,1,IBM),TAU)
          DO IGR=1,NG
            B(:)=FHETXM(IBM,IGR,:)
            CALL ALST2S(NC,NC,NG,FHOMM(1,1,IBM),TAU,B,X)
            FDXM(IBM,IGR,:)=REAL(X(:))
          ENDDO
          CALL ALST2F(NC,NC,NG,FHOMP(1,1,IBM),TAU)
          DO IGR=1,NG
            B(:)=FHETXP(IBM,IGR,:)
            CALL ALST2S(NC,NC,NG,FHOMP(1,1,IBM),TAU,B,X)
            FDXP(IBM,IGR,:)=REAL(X(:))
          ENDDO
          IF(IBM.EQ.NMIX1) THEN
            DO IGR=1,NG
              B(:)=JXP(IBM,IGR,:)
              CALL ALST2S(NC,NC,NG,FHOMP(1,1,IBM),TAU,B,X)
              AFACTOR(IGR,:)=REAL(X(:))
            ENDDO
          ENDIF
          DEALLOCATE(X,B,TAU)
        ENDIF
      ENDDO
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/48H BREERM: DISCONTINUITY FACTORS BEFORE NORMALIZAT,
     1  3HION)')
        DO IBM=1,NMIX1
          WRITE(6,'(/9H MIXTURE=,I5)') IBM
          WRITE(6,20) 'FDXM',FDXM(IBM,:NG,:NG)
          WRITE(6,20) 'FDXP',FDXP(IBM,:NG,:NG)
        ENDDO
      ENDIF
*----
*  COMPUTE ALBEDOS
*----
      IF(ICODE(2).NE.0) THEN
        BETA(:,:)=0.0
        IF(ICODE(2).NE.0) THEN
          IBM=IMIX(LX1)
          DO IGR=1,NG
            DO JGR=1,NG
              BETA(IGR,JGR)=(1.0-2.0*AFACTOR(IGR,JGR))/(1.0+2.0*
     1        AFACTOR(IGR,JGR))
            ENDDO
          ENDDO
        ENDIF
        IF(IPRINT.GT.0) THEN
          WRITE(6,'(/16H BREERM: ALBEDOS)')
          WRITE(6,20) 'BETA',BETA(:NG,:NG)
        ENDIF
      ENDIF
*----
*  NGET NORMALIZATION OF THE DISCONTINUITY FACTORS
*----
      ALLOCATE(WORK2(NG,2*NG))
      DO J=1,LX1-1
        IBM=IMIX(J)
        IBMP=IMIX(J+1)
        DO IGR=1,NG
          DO JGR=1,NG
            WORK2(IGR,JGR)=FDXP(IBM,IGR,JGR)
            WORK2(IGR,NG+JGR)=FDXM(IBMP,IGR,JGR)
          ENDDO
        ENDDO
        CALL ALSBD(NG,NG,WORK2,IER,NG)
        IF(IER.NE.0) CALL XABORT('BREERM: SINGULAR MATRIX(3).')
        DO IGR=1,NG
          ! impose the adf on the fuel assembly side
          IF((J.EQ.J_FUEL).AND.(NGET.EQ.1)) THEN
            FNORM=ADFREF(IGR)
          ELSE IF((J.EQ.J_FUEL).AND.(NGET.EQ.2)) THEN
            FNORM=ADF(IBM,IGR)
          ELSE
            FNORM=FDXP(IBM,IGR,IGR)
          ENDIF
          FDXP(IBM,IGR,:)=0.0
          FDXP(IBM,IGR,IGR)=FNORM
          DO JGR=1,NG
            FDXM(IBMP,IGR,JGR)=REAL(WORK2(IGR,NG+JGR))*FNORM
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(WORK2)
      IF(J_FUEL.GT.0) THEN
        IBMM=0
        DO J=J_FUEL,1,-1
          IBM=IMIX(J)
          IF(J>1) IBMM=IMIX(J-1)
          DO IGR=1,NG
            FNORM=FDXP(IBM,IGR,IGR)/FDXM(IBM,IGR,IGR)
            DO JGR=1,NG
              IF(J>1) FDXP(IBMM,IGR,JGR)=FDXP(IBMM,IGR,JGR)*FNORM
              FDXM(IBM,IGR,JGR)=FDXM(IBM,IGR,JGR)*FNORM
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IBMP=0
      DO J=J_FUEL+1,LX1
        IBM=IMIX(J)
        IF(J<LX1) IBMP=IMIX(J+1)
        DO IGR=1,NG
          FNORM=FDXM(IBM,IGR,IGR)/FDXP(IBM,IGR,IGR)
          DO JGR=1,NG
            IF(J<LX1) FDXM(IBMP,IGR,JGR)=FDXM(IBMP,IGR,JGR)*FNORM
            FDXP(IBM,IGR,JGR)=FDXP(IBM,IGR,JGR)*FNORM
          ENDDO
        ENDDO
      ENDDO
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/48H BREERM: DISCONTINUITY FACTORS AFTER NGET NORMAL,
     1  7HIZATION)')
        DO IBM=1,NMIX1
          WRITE(6,'(/9H MIXTURE=,I5)') IBM
          WRITE(6,20) 'FDXM',FDXM(IBM,:NG,:NG)
          WRITE(6,20) 'FDXP',FDXP(IBM,:NG,:NG)
        ENDDO
      ENDIF
*----
*  SAVE THE OUTPUT NODAL MACROLIB
*----
      ALLOCATE(IJJ(NMIX1),NJJ(NMIX1),IPOS(NMIX1),WORK(NMIX1*NG))
      ISTATE(:)=0
      ISTATE(1)=NG
      ISTATE(2)=NMIX1
      ISTATE(3)=1
      IF(J_FUEL.GT.0) ISTATE(4)=1
      IF(ICODE(2).NE.0) ISTATE(8)=1  ! physical matrix albedo info
      ISTATE(9)=1  ! diffusion coefficient information
      IF(ISPH.EQ.0) ISTATE(12)=4 ! discontinuity factor information
      CALL LCMPUT(IPMAC1,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPMAC1,'ENERGY',NG+1,2,ENER)
      CALL LCMPUT(IPMAC1,'VOLUME',NMIX1,2,VOL)
      IF(ICODE(2).NE.0) CALL LCMPUT(IPMAC1,'ALBEDO',NG*NG,2,BETA)
      IF(ISPH.EQ.0) THEN
        CALL LCMSIX(IPMAC1,'ADF',1)
          NTYPE=2
          HADF(1)='ERM_M'
          HADF(2)='ERM_P'
          CALL LCMPTC(IPMAC1,'HADF',8,NTYPE,HADF)
          CALL LCMPUT(IPMAC1,HADF(1),NMIX1*NG*NG,2,FDXM)
          CALL LCMPUT(IPMAC1,HADF(2),NMIX1*NG*NG,2,FDXP)
        CALL LCMSIX(IPMAC1,' ',2)
      ENDIF
      JPMAC1=LCMLID(IPMAC1,'GROUP',NG)
      DO IGR=1,NG
        KPMAC1=LCMDIL(JPMAC1,IGR)
        DO IBM=1,NMIX1
          WORK(IBM)=VOL(IBM)*FLX(IBM,IGR)
        ENDDO
        CALL LCMPUT(KPMAC1,'FLUX-INTG',NMIX1,2,WORK)
        CALL LCMPUT(KPMAC1,'NTOT0',NMIX1,2,TOT(:,IGR))
        CALL LCMPUT(KPMAC1,'DIFF',NMIX1,2,DC(:,IGR))
        DO IBM=1,NMIX1
          WORK(IBM)=SCAT(IBM,IGR,IGR)
        ENDDO
        CALL LCMPUT(KPMAC1,'SIGW00',NMIX1,2,WORK)
        CALL LCMPUT(KPMAC1,'CHI',NMIX1,2,CHI(:,IGR))
        CALL LCMPUT(KPMAC1,'NUSIGF',NMIX1,2,SIGF(:,IGR))
        IPOSDE=0
        DO IBM=1,NMIX1
          J2=IGR
          J1=IGR
          DO JGR=1,NG
            IF(SCAT(IBM,IGR,JGR).NE.0.0) THEN
              J2=MAX(J2,JGR)
              J1=MIN(J1,JGR)
            ENDIF
          ENDDO
          NJJ(IBM)=J2-J1+1
          IJJ(IBM)=J2
          IPOS(IBM)=IPOSDE+1
          DO JGR=J2,J1,-1
            IPOSDE=IPOSDE+1
            IF(IPOSDE.GT.NG*NMIX1) CALL XABORT('BREERM: SCAT OVERFLOW.')
            WORK(IPOSDE)=SCAT(IBM,IGR,JGR)
          ENDDO
        ENDDO
        CALL LCMPUT(KPMAC1,'SCAT00',IPOSDE,2,WORK)
        CALL LCMPUT(KPMAC1,'NJJS00',NMIX1,1,NJJ)
        CALL LCMPUT(KPMAC1,'IJJS00',NMIX1,1,IJJ)
        CALL LCMPUT(KPMAC1,'IPOS00',NMIX1,1,IPOS)
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK,IPOS,NJJ,IJJ,ADF,SCAT,SIGF,CHI,TOT,DC,FLX,VOL)
      DEALLOCATE(FHOMP,FHOMM,BETA,AFACTOR,FDXP,FDXM,ALPHA,ETA)
      RETURN
   20 FORMAT(1X,A9,1P,10E12.4,/(10X,10E12.4))
      END
