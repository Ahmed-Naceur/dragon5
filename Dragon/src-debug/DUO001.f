*DECK DUO001
      SUBROUTINE DUO001(IPMAC,IPRINT,NMIX,NGRP,NFIS,IDIV,ZKEFF,RHS,LHS,
     > FLUX,AFLUX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Processing one of the two macrolibs and return mixture-dependent
* RHS and LHS matrices.
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
* IPMAC   macrolib.
* IPRINT  print parameter.
* NMIX    number of mixtures.
* NGRP    number of energy groups.
* NFIS    number of fissile isotopes.
* IDIV    type of divergence term processing (=0: no processing;
*         =1: direct processing; =2: adjoint processing;
*         =3: direct-adjoint processing).
*
*Parameters: output
* ZKEFF   effective multiplication factor.
* RHS     absorption macroscopic cross-section matrix.
* LHS     production macroscopic cross-section matrix.
* FLUX    integrated direct flux.
* AFLUX   integrated adjoint flux.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC
      INTEGER IPRINT,NMIX,NGRP,NFIS,IDIV
      REAL ZKEFF,RHS(NGRP,NGRP,NMIX),LHS(NGRP,NGRP,NMIX),
     > FLUX(NGRP,NMIX),AFLUX(NGRP,NMIX)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAC,KPMAC
      DOUBLE PRECISION SUM
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL,GAR,GAR2,DLK,ALK,V,W
      REAL, ALLOCATABLE, DIMENSION(:,:) :: NUF
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: CHI
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(NMIX),NJJ(NMIX),IPOS(NMIX))
      ALLOCATE(VOL(NMIX),GAR(NMIX),CHI(NMIX,NFIS,NGRP),NUF(NMIX,NFIS),
     > GAR2(NMIX*NGRP))
*----
*  COMPUTE THE RHS AND LHS MATRICES
*----
      CALL XDRSET(RHS,NGRP*NGRP*NMIX,0.0)
      CALL XDRSET(LHS,NGRP*NGRP*NMIX,0.0)
      CALL LCMGET(IPMAC,'K-EFFECTIVE',ZKEFF)
      IF(IPRINT.GT.1) WRITE(6,'(35H DUO001: EFFECTIVE MULTIPLICATION F,
     > 6HACTOR=,1P,E12.5)') ZKEFF
      CALL LCMGET(IPMAC,'VOLUME',VOL)
      JPMAC=LCMGID(IPMAC,'GROUP')
      DO IGR=1,NGRP
        KPMAC=LCMGIL(JPMAC,IGR)
        CALL LCMGET(KPMAC,'CHI',CHI(1,1,IGR))
      ENDDO
      DO IGR=1,NGRP
        KPMAC=LCMGIL(JPMAC,IGR)
        CALL LCMGET(KPMAC,'FLUX-INTG',GAR)
        DO IBM=1,NMIX
          FLUX(IGR,IBM)=GAR(IBM)/VOL(IBM)
        ENDDO
        CALL LCMLEN(KPMAC,'NWAT0',ILONG,ITYLCM)
        IF(ILONG.EQ.NMIX) THEN
          CALL LCMGET(KPMAC,'NWAT0',GAR)
          DO IBM=1,NMIX
            AFLUX(IGR,IBM)=GAR(IBM)
          ENDDO
        ELSE
          CALL XDRSET(AFLUX(1,IBM),NMIX,1.0)
        ENDIF
        CALL LCMGET(KPMAC,'NTOT0',GAR)
        CALL LCMGET(KPMAC,'SCAT00',GAR2)
        CALL LCMGET(KPMAC,'NJJS00',NJJ)
        CALL LCMGET(KPMAC,'IJJS00',IJJ)
        CALL LCMGET(KPMAC,'IPOS00',IPOS)
        DO IBM=1,NMIX
          IPOSDE=IPOS(IBM)
          DO JGR=IJJ(IBM),IJJ(IBM)-NJJ(IBM)+1,-1
            RHS(IGR,JGR,IBM)=RHS(IGR,JGR,IBM)-GAR2(IPOSDE) ! IGR <-- JGR
            IPOSDE=IPOSDE+1
          ENDDO
          RHS(IGR,IGR,IBM)=RHS(IGR,IGR,IBM)+GAR(IBM)
        ENDDO
        CALL LCMGET(KPMAC,'NUSIGF',NUF)
        DO IBM=1,NMIX
          DO IFIS=1,NFIS
            DO JGR=1,NGRP
              LHS(JGR,IGR,IBM)=LHS(JGR,IGR,IBM)+CHI(IBM,IFIS,JGR)*
     >        NUF(IBM,IFIS)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
*----
*  INTRODUCE THE DIRECT OR ADJOINT DIVERGENCE COMPONENT IN THE RHS
*  MATRIX
*----
      DO IBM=1,NMIX
        IF(IDIV.EQ.1) THEN
          DO JGR=1,NGRP
            SUM=0.0D0
            DO IGR=1,NGRP
              SUM=SUM+(RHS(JGR,IGR,IBM)-LHS(JGR,IGR,IBM)/ZKEFF)*
     >        FLUX(IGR,IBM)
            ENDDO
            RHS(JGR,JGR,IBM)=RHS(JGR,JGR,IBM)-REAL(SUM)/FLUX(JGR,IBM)
          ENDDO
        ELSE IF(IDIV.EQ.2) THEN
          DO IGR=1,NGRP
            SUM=0.0D0
            DO JGR=1,NGRP
              SUM=SUM+(RHS(JGR,IGR,IBM)-LHS(JGR,IGR,IBM)/ZKEFF)*
     >        AFLUX(JGR,IBM)
            ENDDO
            RHS(IGR,IGR,IBM)=RHS(IGR,IGR,IBM)-REAL(SUM)/AFLUX(IGR,IBM)
          ENDDO
        ELSE IF(IDIV.EQ.3) THEN
          ALLOCATE(DLK(NGRP),ALK(NGRP))
          DO JGR=1,NGRP
            SUM=0.0D0
            DO IGR=1,NGRP
              SUM=SUM+(RHS(JGR,IGR,IBM)-LHS(JGR,IGR,IBM)/ZKEFF)*
     >        FLUX(IGR,IBM)
            ENDDO
            DLK(JGR)=REAL(SUM)
          ENDDO
          DO IGR=1,NGRP
            SUM=0.0D0
            DO JGR=1,NGRP
              SUM=SUM+(RHS(JGR,IGR,IBM)-LHS(JGR,IGR,IBM)/ZKEFF)*
     >        AFLUX(JGR,IBM)
            ENDDO
            ALK(IGR)=REAL(SUM)
          ENDDO
          ALLOCATE(V(NGRP),W(NGRP))
          CALL DUO005(NGRP,DLK,ALK,FLUX(1,IBM),AFLUX(1,IBM),V,W)
          DO IGR=1,NGRP
            DO JGR=1,NGRP
              RHS(IGR,JGR,IBM)=RHS(IGR,JGR,IBM)-V(IGR)-W(JGR)
            ENDDO
          ENDDO
          DEALLOCATE(W,V,ALK,DLK)
        ENDIF
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAR2,NUF,CHI,GAR,VOL)
      DEALLOCATE(IPOS,NJJ,IJJ)
      RETURN
      END
