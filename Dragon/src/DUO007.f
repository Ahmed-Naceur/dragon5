*DECK DUO007
      SUBROUTINE DUO007(IPLIB,IPRINT,NISOT,NGRP,LENER,RHSI1,RHSI2,
     > LHSI2,FLUXI2,AFLUXI1,RHOREA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print the detail of isotopic delta-rho discrepancies between two
* calculations for a single reaction
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
* NISOT   number of isotopes.
* NGRP    number of energy groups.
* LENER   energy group analysis flag.
* RHSI1   absorption macroscopic cross-section matrix for the first
*         calculation.
* RHSI2   absorption macroscopic cross-section matrix for the second
*         calculation.
* LHSI2   production macroscopic cross-section matrix for the second
*         calculation.
* FLUXI2  flux for the second calculation.
* AFLUXI1 adjoint flux for the first calculation.
*
*Parameters: output
* RHOREA  total delta-rho for the reaction.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER IPRINT,NISOT,NGRP
      LOGICAL LENER
      REAL RHSI1(NGRP,NGRP,NISOT),RHSI2(NGRP,NGRP,NISOT),
     > LHSI2(NGRP,NGRP,NISOT),FLUXI2(NGRP,NISOT),AFLUXI1(NGRP,NISOT),
     > RHOREA
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION ZDEN,RHO12
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
      ALLOCATE(RHO(NGRP,NISOT),RHO1(NGRP),RHO2(NISOT))
*----
*  PERTURBATIVE ANALYSIS WITH THE CLIO FORMULA
*----
      CALL LCMGET(IPLIB,'ISOTOPESUSED',IHUSED)
      CALL XDRSET(RHO,NGRP*NISOT,0.0)
      CALL XDRSET(RHO1,NGRP,0.0)
      CALL XDRSET(RHO2,NISOT,0.0)
      ZDEN=0.0D0
      DO ISOT=1,NISOT
        DO IGR=1,NGRP
          DO JGR=1,NGRP
            DRHS=(RHSI2(JGR,IGR,ISOT)-RHSI1(JGR,IGR,ISOT))
            RHO(IGR,ISOT)=RHO(IGR,ISOT)+DRHS*AFLUXI1(JGR,ISOT)*
     >      FLUXI2(IGR,ISOT)
            ZDEN=ZDEN+LHSI2(JGR,IGR,ISOT)*AFLUXI1(JGR,ISOT)*
     >      FLUXI2(IGR,ISOT)
          ENDDO
        ENDDO
      ENDDO
      RHO12=0.0D0
      DO ISOT=1,NISOT
        DO IGR=1,NGRP
          RHO(IGR,ISOT)=RHO(IGR,ISOT)*1.0E5/REAL(ZDEN)
          RHO1(IGR)=RHO1(IGR)+RHO(IGR,ISOT)
          RHO2(ISOT)=RHO2(ISOT)+RHO(IGR,ISOT)
          RHO12=RHO12+RHO(IGR,ISOT)
        ENDDO
      ENDDO
*----
*  PRINT DELTA-RHO
*----
      IF(LENER) THEN
        WRITE(6,'(/48H DUO007: DELTA-RHO ISOTOPIC-MULTIGROUP DISCREPAN,
     >  12HCIES IN PCM:)')
        WRITE(6,'(14X,9I12)') (IGR,IGR=1,NGRP)
        IF(IPRINT.GT.0) THEN
          DO ISOT=1,NISOT
            WRITE(TEXT12,'(3A4)') (IHUSED(I0,ISOT),I0=1,3)
            WRITE(6,110) TEXT12,(RHO(IGR,ISOT),IGR=1,NGRP)
          ENDDO
        ENDIF
        WRITE(6,120) (RHO1(IGR),IGR=1,NGRP)
      ENDIF
      WRITE(6,'(/49H DUO007: DELTA-RHO ISOTOPIC DISCREPANCIES IN PCM:)')
      DO ISOT=1,NISOT
        WRITE(TEXT12,'(3A4)') (IHUSED(I0,ISOT),I0=1,3)
        WRITE(6,110) TEXT12,RHO2(ISOT)
      ENDDO
      WRITE(6,120) RHO12
      RHOREA=REAL(RHO12)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(RHO2,RHO1,RHO)
      DEALLOCATE(IHUSED)
      RETURN
*
  110 FORMAT(1X,A12,1X,9F12.2/(14X,9F12.2))
  120 FORMAT(/14H  *** SUM *** ,9F12.2/(14X,9F12.2))
      END
