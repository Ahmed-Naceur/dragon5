*DECK DUO002
      SUBROUTINE DUO002(IPRINT,NMIX,NGRP,LENER,ZKEFF1,ZKEFF2,RHS1,RHS2,
     > LHS1,LHS2,FLUX2,AFLUX1)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print the detail of mixture delta-rho discrepancies between two
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
* IPRINT  print parameter.
* NMIX    number of mixtures.
* NGRP    number of energy groups.
* LENER   energy group analysis flag.
* ZKEFF1  effective multiplication factor of the first calculation.
* ZKEFF2  effective multiplication factor of the second calculation.
* RHS1    absorption macroscopic cross-section matrix for the first
*         calculation.
* RHS2    absorption macroscopic cross-section matrix for the second
*         calculation.
* LHS1   production macroscopic cross-section matrix for the first
*         calculation.
* LHS2    production macroscopic cross-section matrix for the second
*         calculation.
* FLUX2   flux for the second calculation.
* AFLUX1  adjoint flux for the first calculation.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPRINT,NMIX,NGRP
      LOGICAL LENER
      REAL ZKEFF1,ZKEFF2,RHS1(NGRP,NGRP,NMIX),RHS2(NGRP,NGRP,NMIX),
     > LHS1(NGRP,NGRP,NMIX),LHS2(NGRP,NGRP,NMIX),FLUX2(NGRP,NMIX),
     > AFLUX1(NGRP,NMIX)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION ZNUM,ZDEN,RHO12
      CHARACTER TEXT12*12
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: RHO1,RHO2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RHO
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(RHO(NGRP,NMIX),RHO1(NGRP),RHO2(NMIX))
*----
*  RAYLEIGH RATIO FOR THE FIRST SYSTEM
*----
      IF(IPRINT.GT.1) THEN
        ZNUM=0.0D0
        ZDEN=0.0D0
        DO IBM=1,NMIX
          DO IGR=1,NGRP
            DO JGR=1,NGRP
              ZNUM=ZNUM+LHS1(JGR,IGR,IBM)*AFLUX1(JGR,IBM)*FLUX2(IGR,IBM)
              ZDEN=ZDEN+RHS1(JGR,IGR,IBM)*AFLUX1(JGR,IBM)*FLUX2(IGR,IBM)
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
        DO IBM=1,NMIX
          DO IGR=1,NGRP
            DO JGR=1,NGRP
              ZNUM=ZNUM+LHS2(JGR,IGR,IBM)*AFLUX1(JGR,IBM)*FLUX2(IGR,IBM)
              ZDEN=ZDEN+RHS2(JGR,IGR,IBM)*AFLUX1(JGR,IBM)*FLUX2(IGR,IBM)
            ENDDO
          ENDDO
        ENDDO
        WRITE(6,100) 2,ZNUM/ZDEN,ZKEFF2
      ENDIF
*----
*  PERTURBATIVE ANALYSIS WITH THE CLIO FORMULA
*----
      CALL XDRSET(RHO,NGRP*(NMIX),0.0)
      CALL XDRSET(RHO1,NGRP,0.0)
      CALL XDRSET(RHO2,NMIX,0.0)
      ZDEN=0.0D0
      DO IBM=1,NMIX
        DO IGR=1,NGRP
          DO JGR=1,NGRP
            DRHS=(RHS2(JGR,IGR,IBM)-RHS1(JGR,IGR,IBM))
            DLHS=(LHS2(JGR,IGR,IBM)-LHS1(JGR,IGR,IBM))
            RHO(IGR,IBM)=RHO(IGR,IBM)+(DRHS-DLHS/ZKEFF1)*
     >      AFLUX1(JGR,IBM)*FLUX2(IGR,IBM)
            ZDEN=ZDEN+LHS2(JGR,IGR,IBM)*AFLUX1(JGR,IBM)*
     >      FLUX2(IGR,IBM)
          ENDDO
        ENDDO
      ENDDO
      RHO12=0.0D0
      DO IBM=1,NMIX
        DO IGR=1,NGRP
          RHO(IGR,IBM)=RHO(IGR,IBM)*1.0E5/REAL(ZDEN)
          RHO1(IGR)=RHO1(IGR)+RHO(IGR,IBM)
          RHO2(IBM)=RHO2(IBM)+RHO(IGR,IBM)
          RHO12=RHO12+RHO(IGR,IBM)
        ENDDO
      ENDDO
      DELLAM=((1./ZKEFF2)-(1./ZKEFF1))*1.0E5
      DELTA=ABS(DELLAM-REAL(RHO12))
      IF(DELTA.GT.2.0) CALL XABORT('DUO002: FAILURE OF DUO: MODULE TO '
     > //'ANALYSE THE DELTA-RHO DISCREPANCY WITHIN 2 PCM.')
*----
*  PRINT DELTA-RHO
*----
      IF(LENER) THEN
        WRITE(6,'(/47H DUO002: DELTA-RHO MIXTURE-MULTIGROUP DISCREPAN,
     >  12HCIES IN PCM:)')
        WRITE(6,'(14X,9I12)') (IGR,IGR=1,NGRP)
        IF(IPRINT.GT.0) THEN
          DO IBM=1,NMIX
            WRITE(TEXT12,'(4HMIXT,I6.6)') IBM
            WRITE(6,110) TEXT12,(RHO(IGR,IBM),IGR=1,NGRP)
          ENDDO
        ENDIF
        WRITE(6,120) (RHO1(IGR),IGR=1,NGRP)
      ENDIF
      WRITE(6,'(/48H DUO002: DELTA-RHO MIXTURE-DISCREPANCIES IN PCM:)')
      DO IBM=1,NMIX
        WRITE(TEXT12,'(4HMIXT,I6.6)') IBM
        WRITE(6,110) TEXT12,RHO2(IBM)
      ENDDO
      WRITE(6,120) RHO12
      WRITE(6,'(14H  *** SUM *** ,F12.2,8H (EXACT))') DELLAM
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(RHO2,RHO1,RHO)
      RETURN
*
  100 FORMAT(16H DUO002: SYSTEM=,I2,21H DIRECT-ADJOINT KEFF=,1P,E13.5,
     > 13H EXACT VALUE=,E13.5)
  110 FORMAT(1X,A12,1X,9F12.2/(14X,9F12.2))
  120 FORMAT(/14H  *** SUM *** ,9F12.2/(14X,9F12.2))
      END
