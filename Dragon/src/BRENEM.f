*DECK BRENEM
      SUBROUTINE BRENEM(IPMAC1,NG,LX1,NMIX1,ITRIAL,IMIX,ICODE,ISPH,
     1 ZKEFF,B2,ENER,VOL1,FLX1,DC1,TOT1,CHI1,SIGF1,SCAT1,JXM,JXP,FHETXM,
     2 FHETXP,ADF1,NGET,ADFREF,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Implement the 1D FD-NEM reflector model.
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
      INTEGER NG,LX1,NMIX1,ITRIAL(NG),IMIX(LX1),ICODE(2),ISPH,NGET,
     1 IPRINT
      REAL ZKEFF,B2,ENER(NG+1),VOL1(NMIX1),FLX1(NMIX1,NG),DC1(NMIX1,NG),
     1 TOT1(NMIX1,NG),CHI1(NMIX1,NG),SIGF1(NMIX1,NG),
     2 SCAT1(NMIX1,NG,NG),JXM(NMIX1,NG),JXP(NMIX1,NG),FHETXM(NMIX1,NG),
     3 FHETXP(NMIX1,NG),ADF1(NMIX1,NG),ADFREF(NG)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE)
      REAL SX(4),A11(4,4),Q(5)
      CHARACTER HADF*8
      TYPE(C_PTR) JPMAC1,KPMAC1
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK,ETA,AFACTOR,BETA
      REAL, ALLOCATABLE, DIMENSION(:,:) :: AB,ALPHA,FDXM,FDXP
*----
*  LOOP OVER EQUIVALENT REFLECTOR NODES
*----
      ALLOCATE(ETA(NG),AB(4*NG,4*NG+1),ALPHA(4,NG),FDXM(NMIX1,NG),
     1 FDXP(NMIX1,NG),AFACTOR(NG),BETA(NG))
*----
*  SET AND SOLVE NODAL SYSTEM
*----
      J_FUEL=0
      DO J=1,LX1
        IBM=IMIX(J)
        IF(IBM.EQ.0) CYCLE
        DO IGR=1,NG
          IF(SIGF1(IBM,IGR).GT.0.0) J_FUEL=J
          DIFF=DC1(IBM,IGR)
          SIGR=TOT1(IBM,IGR)+B2*DIFF-SCAT1(IBM,IGR,IGR)
          ETA(IGR)=VOL1(IBM)*SQRT(SIGR/DIFF)
          DO JGR=1,NG
            IF(JGR.EQ.IGR) THEN
             SIGR=SIGR-CHI1(IBM,IGR)*SIGF1(IBM,IGR)/ZKEFF
             CALL BRESS1(ITRIAL(IGR),VOL1(IBM),DIFF,SIGR,A11)
            ELSE
             SIGT=-SCAT1(IBM,JGR,IGR)-CHI1(IBM,JGR)*SIGF1(IBM,IGR)/ZKEFF
             CALL BRESS2(ITRIAL(IGR),VOL1(IBM),DIFF,SIGR,SIGT,A11)
            ENDIF
            DO K1=1,4
              DO K2=1,4
                AB((JGR-1)*4+K1,(IGR-1)*4+K2)=A11(K1,K2)
             ENDDO
            ENDDO
          ENDDO
          SX = (/0.0,0.0,JXM(IBM,IGR),JXP(IBM,IGR)/) 
          DO K1 =1,4
            AB((IGR-1)*4+K1,4*NG+1)=SX(K1)
          ENDDO
        ENDDO
        CALL ALSB(4*NG,1,AB,IER,4*NG)
        IF (IER.NE.0) CALL XABORT('BRENEM: ALBS FAILURE')
        DO IGR=1,NG
          DO I=1,4
            ALPHA(I,IGR)=AB((IGR-1)*4+I,4*NG+1)
          ENDDO
        ENDDO
        IF(IPRINT.GT.1) THEN
          WRITE(6,'(/9H MIXTURE=,I5)') J
          WRITE(6,20) 'ALPHA',ALPHA(:4,:NG)
        ENDIF
*----
*  COMPUTE DISCONTINUITY FACTORS
*----
        DO IGR=1,NG
          IF (ITRIAL(IGR) == 1) THEN
            Q(1) = ALPHA(2,IGR)/2.
            FHOMM=-ALPHA(1,IGR)/2.+FLX1(IBM,IGR)+Q(1)
            FHOMP=ALPHA(1,IGR)/2.+FLX1(IBM,IGR)+Q(1)
          ELSE
            Q(1) = ETA(IGR)/2.
            Q(2) = SINH(Q(1))
            Q(3) = ALPHA(2,IGR)/2.
            Q(4) = ALPHA(3,IGR)*Q(2)
            Q(5) = ALPHA(4,IGR)*(COSH(Q(1)) - (2*Q(2))/ETA(IGR))
            FHOMM=-ALPHA(1,IGR)/2.+FLX1(IBM,IGR)+Q(3)-Q(4)+Q(5)
            FHOMP=ALPHA(1,IGR)/2.+FLX1(IBM,IGR)+Q(3)+Q(4)+Q(5)
          ENDIF
          FDXM(IBM,IGR)=FHETXM(IBM,IGR)/FHOMM
          FDXP(IBM,IGR)=FHETXP(IBM,IGR)/FHOMP
        ENDDO
      ENDDO
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/48H BRENEM: DISCONTINUITY FACTORS BEFORE NORMALIZAT,
     1  3HION)')
        DO IBM=1,NMIX1
          WRITE(6,'(/9H MIXTURE=,I5)') IBM
          WRITE(6,20) 'FDXM',FDXM(IBM,:NG)
          WRITE(6,20) 'FDXP',FDXP(IBM,:NG)
        ENDDO
      ENDIF
*----
*  COMPUTE ALBEDOS
*----
      IF(ICODE(2).NE.0) THEN
        BETA(:)=0.0
        IF(ICODE(2).NE.0) THEN
          IBM=IMIX(LX1)
          DO IGR=1,NG
            AFACTOR(IGR)=FDXP(IBM,IGR)*JXP(IBM,IGR)/FHETXP(IBM,IGR)
            BETA(IGR)=(1.0-2.0*AFACTOR(IGR))/(1.0+2.0*AFACTOR(IGR))
          ENDDO
        ENDIF
        IF(IPRINT.GT.0) THEN
          WRITE(6,'(/16H BRENEM: ALBEDOS)')
          WRITE(6,20) 'BETA',BETA(:NG)
        ENDIF
      ENDIF
*----
*  NGET NORMALIZATION OF THE DISCONTINUITY FACTORS
*----
      IF(J_FUEL.GT.0) THEN
        IF(NGET.GT.0) THEN
          IBMP=0
          IBM=IMIX(J_FUEL)
          IF(J_FUEL<LX1) IBMP=IMIX(J_FUEL+1)
          DO IGR=1,NG
             ! impose the adf on the fuel assembly side
             IF(NGET.EQ.1) THEN
               FNORM=ADFREF(IGR)/FDXP(IBM,IGR)
             ELSE
               FNORM=ADF1(IBM,IGR)/FDXP(IBM,IGR)
             ENDIF
             FDXP(IBM,IGR)=FDXP(IBM,IGR)*FNORM
             IF(J_FUEL<LX1) FDXM(IBMP,IGR)=FDXM(IBMP,IGR)*FNORM
          ENDDO
        ENDIF
        IBMM=0
        DO J=J_FUEL,1,-1
          IBM=IMIX(J)
          IF(J>1) IBMM=IMIX(J-1)
          DO IGR=1,NG
            IF(J>1) FDXP(IBMM,IGR)=FDXP(IBMM,IGR)*FDXP(IBM,IGR)/
     1      FDXM(IBM,IGR)
            FDXM(IBM,IGR)=FDXP(IBM,IGR)
          ENDDO
        ENDDO
      ENDIF
      IBMP=0
      DO J=J_FUEL+1,LX1
        IBM=IMIX(J)
        IF(J<LX1) IBMP=IMIX(J+1)
        DO IGR=1,NG
          IF(J<LX1) FDXM(IBMP,IGR)=FDXM(IBMP,IGR)*FDXM(IBM,IGR)/
     1    FDXP(IBM,IGR)
          FDXP(IBM,IGR)=FDXM(IBM,IGR)
        ENDDO
      ENDDO
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/48H BRENEM: DISCONTINUITY FACTORS AFTER NGET NORMAL,
     1  7HIZATION)')
        DO IBM=1,NMIX1
          WRITE(6,'(/9H MIXTURE=,I5)') IBM
          WRITE(6,20) 'FDX',FDXM(IBM,:NG)
        ENDDO
      ENDIF
*----
*  APPLY SPH FACTORS
*----
      IF(ISPH.EQ.1) THEN
        DO IGR=1,NG
          DO J=1,LX1
            IBM=IMIX(J)
            TOT1(IBM,IGR)=TOT1(IBM,IGR)/FDXM(IBM,IGR)
            DC1(IBM,IGR)=DC1(IBM,IGR)/FDXM(IBM,IGR)
            SIGF1(IBM,IGR)=SIGF1(IBM,IGR)/FDXM(IBM,IGR)
            DO JGR=1,NG
              SCAT1(IBM,IGR,JGR)=SCAT1(IBM,IGR,JGR)/FDXM(IBM,JGR)
            ENDDO
          ENDDO
        ENDDO
        IF(ICODE(2).NE.0) THEN
          BETA(:)=0.0
          IF(ICODE(2).NE.0) THEN
            IBM=IMIX(LX1)
            DO IGR=1,NG
              AFACTOR(IGR)=AFACTOR(IGR)/FDXM(IBM,IGR)
              BETA(IGR)=(1.0-2.0*AFACTOR(IGR))/(1.0+2.0*AFACTOR(IGR))
            ENDDO
          ENDIF
          IF(IPRINT.GT.0) THEN
            WRITE(6,'(/30H BRENEM: SPH CORRECTED ALBEDOS)')
            WRITE(6,20) 'BETA',BETA(:NG)
          ENDIF
        ENDIF
      ENDIF
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/31H BRENEM: DIFFUSION COEFFICIENTS)')
        DO IBM=1,NMIX1
          WRITE(6,'(/9H MIXTURE=,I5)') IBM
          WRITE(6,20) 'DIFF',DC1(IBM,:NG)
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
      IF(ICODE(2).NE.0) ISTATE(8)=1  ! physical albedo information
      ISTATE(9)=1  ! diffusion coefficient information
      IF(ISPH.EQ.0) ISTATE(12)=3 ! discontinuity factor information
      CALL LCMPUT(IPMAC1,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPMAC1,'ENERGY',NG+1,2,ENER)
      CALL LCMPUT(IPMAC1,'VOLUME',NMIX1,2,VOL1)
      IF(ICODE(2).NE.0) CALL LCMPUT(IPMAC1,'ALBEDO',NG,2,BETA)
      IF(ISPH.EQ.0) THEN
        CALL LCMSIX(IPMAC1,'ADF',1)
          NTYPE=1
          HADF='FD_B'
          CALL LCMPTC(IPMAC1,'HADF',8,1,HADF)
          CALL LCMPUT(IPMAC1,HADF,NMIX1*NG,2,FDXM)
        CALL LCMSIX(IPMAC1,' ',2)
      ENDIF
      JPMAC1=LCMLID(IPMAC1,'GROUP',NG)
      DO IGR=1,NG
        KPMAC1=LCMDIL(JPMAC1,IGR)
        DO IBM=1,NMIX1
          WORK(IBM)=VOL1(IBM)*FLX1(IBM,IGR)
        ENDDO
        CALL LCMPUT(KPMAC1,'FLUX-INTG',NMIX1,2,WORK)
        CALL LCMPUT(KPMAC1,'NTOT0',NMIX1,2,TOT1(:,IGR))
        CALL LCMPUT(KPMAC1,'DIFF',NMIX1,2,DC1(:,IGR))
        DO IBM=1,NMIX1
          WORK(IBM)=SCAT1(IBM,IGR,IGR)
        ENDDO
        CALL LCMPUT(KPMAC1,'SIGW00',NMIX1,2,WORK)
        CALL LCMPUT(KPMAC1,'CHI',NMIX1,2,CHI1(:,IGR))
        CALL LCMPUT(KPMAC1,'NUSIGF',NMIX1,2,SIGF1(:,IGR))
        IF(ISPH.EQ.1) THEN
          DO IBM=1,NMIX1
            WORK(IBM)=1.0/FDXM(IBM,IGR)
          ENDDO
          CALL LCMPUT(KPMAC1,'NSPH',NMIX1,2,WORK)
        ENDIF
        IPOSDE=0
        DO IBM=1,NMIX1
          J2=IGR
          J1=IGR
          DO JGR=1,NG
            IF(SCAT1(IBM,IGR,JGR).NE.0.0) THEN
              J2=MAX(J2,JGR)
              J1=MIN(J1,JGR)
            ENDIF
          ENDDO
          NJJ(IBM)=J2-J1+1
          IJJ(IBM)=J2
          IPOS(IBM)=IPOSDE+1
          DO JGR=J2,J1,-1
            IPOSDE=IPOSDE+1
            IF(IPOSDE.GT.NG*NMIX1) CALL XABORT('BRENEM: SCAT OVERFLOW.')
            WORK(IPOSDE)=SCAT1(IBM,IGR,JGR)
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
      DEALLOCATE(WORK,IPOS,NJJ,IJJ,BETA,AFACTOR,FDXP,FDXM,ALPHA,AB,ETA)
      RETURN
   20 FORMAT(1X,A9,1P,10E12.4,/(10X,10E12.4))
      END
