*DECK DUO004
      SUBROUTINE DUO004(IPLIB,IPRINT,NMIX,NISOT,NGRP,LENER,ZKEFF1,
     > ZKEFF2,RHSI1,RHSI2,LHSI1,LHSI2,FLUXI2,AFLUXI1)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print the detail of isotopic delta-rho discrepancies between two
* calculations
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIB   microlib.
* IPRINT  print parameter.
* NMIX    number of mixtures.
* NISOT   number of isotopes.
* NGRP    number of energy groups.
* LENER   energy group analysis flag.
* ZKEFF1  effective multiplication factor of the first calculation.
* ZKEFF2  effective multiplication factor of the second calculation.
* RHSI1   absorption macroscopic cross-section matrix for the first
*         calculation.
* RHSI2   absorption macroscopic cross-section matrix for the second
*         calculation.
* LHSI1   production macroscopic cross-section matrix for the first
*         calculation.
* LHSI2   production macroscopic cross-section matrix for the second
*         calculation.
* FLUXI2  flux for the second calculation.
* AFLUXI1 adjoint flux for the first calculation.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER IPRINT,NMIX,NISOT,NGRP
      LOGICAL LENER
      REAL ZKEFF1,ZKEFF2,RHSI1(NGRP,NGRP,NISOT+NMIX),
     > RHSI2(NGRP,NGRP,NISOT+NMIX),LHSI1(NGRP,NGRP,NISOT+NMIX),
     > LHSI2(NGRP,NGRP,NISOT+NMIX),FLUXI2(NGRP,NISOT+NMIX),
     > AFLUXI1(NGRP,NISOT+NMIX)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION ZNUM,ZDEN,RHO12
      CHARACTER TEXT12*12
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IHUSED
      REAL, ALLOCATABLE, DIMENSION(:) :: RHO1,RHO2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RHO
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IHUSED(3,NISOT))
      ALLOCATE(RHO(NGRP,NISOT+NMIX),RHO1(NGRP),RHO2(NISOT+NMIX))
*----
*  RAYLEIGH RATIO FOR THE FIRST SYSTEM
*----
      IF(IPRINT.GT.1) THEN
        ZNUM=0.0D0
        ZDEN=0.0D0
        DO ISOT=1,NISOT+NMIX
          DO IGR=1,NGRP
            DO JGR=1,NGRP
              ZNUM=ZNUM+LHSI1(JGR,IGR,ISOT)*AFLUXI1(JGR,ISOT)*
     >        FLUXI2(IGR,ISOT)
              ZDEN=ZDEN+RHSI1(JGR,IGR,ISOT)*AFLUXI1(JGR,ISOT)*
     >        FLUXI2(IGR,ISOT)
            ENDDO
          ENDDO
        ENDDO
        WRITE(6,100) 1,ZNUM/ZDEN,ZKEFF1
      ENDIF
*----
*  RAYLEIGH RATIO FOR THE SECOND SYSTEM
*----
      IF(IPRINT.GT.1) THEN
        ZNUM=0.0D0
        ZDEN=0.0D0
        DO ISOT=1,NISOT+NMIX
          DO IGR=1,NGRP
            DO JGR=1,NGRP
              ZNUM=ZNUM+LHSI2(JGR,IGR,ISOT)*AFLUXI1(JGR,ISOT)*
     >        FLUXI2(IGR,ISOT)
              ZDEN=ZDEN+RHSI2(JGR,IGR,ISOT)*AFLUXI1(JGR,ISOT)*
     >        FLUXI2(IGR,ISOT)
            ENDDO
          ENDDO
        ENDDO
        WRITE(6,100) 2,ZNUM/ZDEN,ZKEFF2
      ENDIF
*----
*  PERTURBATIVE ANALYSIS WITH THE CLIO FORMULA
*----
      CALL LCMGET(IPLIB,'ISOTOPESUSED',IHUSED)
      CALL XDRSET(RHO,NGRP*(NISOT+NMIX),0.0)
      CALL XDRSET(RHO1,NGRP,0.0)
      CALL XDRSET(RHO2,NISOT+NMIX,0.0)
      ZDEN=0.0D0
      DO ISOT=1,NISOT+NMIX
        DO IGR=1,NGRP
          DO JGR=1,NGRP
            DRHS=(RHSI2(JGR,IGR,ISOT)-RHSI1(JGR,IGR,ISOT))
            DLHS=(LHSI2(JGR,IGR,ISOT)-LHSI1(JGR,IGR,ISOT))
            RHO(IGR,ISOT)=RHO(IGR,ISOT)+(DRHS-DLHS/ZKEFF1)*
     >      AFLUXI1(JGR,ISOT)*FLUXI2(IGR,ISOT)
            ZDEN=ZDEN+LHSI2(JGR,IGR,ISOT)*AFLUXI1(JGR,ISOT)*
     >      FLUXI2(IGR,ISOT)
          ENDDO
        ENDDO
      ENDDO
      RHO12=0.0D0
      DO ISOT=1,NISOT+NMIX
        DO IGR=1,NGRP
          RHO(IGR,ISOT)=RHO(IGR,ISOT)*1.0E5/REAL(ZDEN)
          RHO1(IGR)=RHO1(IGR)+RHO(IGR,ISOT)
          RHO2(ISOT)=RHO2(ISOT)+RHO(IGR,ISOT)
          RHO12=RHO12+RHO(IGR,ISOT)
        ENDDO
      ENDDO
      DELLAM=((1./ZKEFF2)-(1./ZKEFF1))*1.0E5
      DELTA=ABS(DELLAM-REAL(RHO12))
      IF(DELTA.GT.2.0) CALL XABORT('DUO004: FAILURE OF DUO: MODULE TO '
     > //'ANALYSE THE DELTA-RHO DISCREPANCY WITHIN 2 PCM.')
*----
*  PRINT DELTA-RHO
*----
      IF(LENER) THEN
        WRITE(6,'(/48H DUO004: DELTA-RHO ISOTOPIC-MULTIGROUP DISCREPAN,
     >  12HCIES IN PCM:)')
        WRITE(6,'(14X,9I12)') (IGR,IGR=1,NGRP)
        IF(IPRINT.GT.0) THEN
          DO ISOT=1,NISOT
            WRITE(TEXT12,'(3A4)') (IHUSED(I0,ISOT),I0=1,3)
            WRITE(6,110) TEXT12,(RHO(IGR,ISOT),IGR=1,NGRP)
          ENDDO
          DO IBM=1,NMIX
            WRITE(TEXT12,'(4HLEAK,I6.6)') IBM
            WRITE(6,110) TEXT12,(RHO(IGR,NISOT+IBM),IGR=1,NGRP)
          ENDDO
        ENDIF
        WRITE(6,120) (RHO1(IGR),IGR=1,NGRP)
      ENDIF
      WRITE(6,'(/49H DUO004: DELTA-RHO ISOTOPIC DISCREPANCIES IN PCM:)')
      DO ISOT=1,NISOT
        WRITE(TEXT12,'(3A4)') (IHUSED(I0,ISOT),I0=1,3)
        WRITE(6,110) TEXT12,RHO2(ISOT)
      ENDDO
      DO IBM=1,NMIX
        WRITE(TEXT12,'(4HLEAK,I6.6)') IBM
        WRITE(6,110) TEXT12,RHO2(NISOT+IBM)
      ENDDO
      WRITE(6,120) RHO12
      WRITE(6,'(14H  *** SUM *** ,F12.2,8H (EXACT))') DELLAM
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(RHO2,RHO1,RHO)
      DEALLOCATE(IHUSED)
      RETURN
*
  100 FORMAT(16H DUO004: SYSTEM=,I2,21H DIRECT-ADJOINT KEFF=,1P,E13.5,
     > 13H EXACT VALUE=,E13.5)
  110 FORMAT(1X,A12,1X,9F12.2/(14X,9F12.2))
  120 FORMAT(/14H  *** SUM *** ,9F12.2/(14X,9F12.2))
      END
