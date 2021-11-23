*DECK DUO003
      SUBROUTINE DUO003(IPLIB,IPRINT,NMIX,NISOT,NGRP,IDIV,ZKEFF,RHS,LHS,
     > FLUX,AFLUX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Processing one of the two microlibs and return isotope-dependent
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
* IPLIB   microlib.
* IPRINT  print parameter.
* NMIX    number of mixtures.
* NISOT   number of isotopes.
* NGRP    number of energy groups.
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
      TYPE(C_PTR) IPLIB
      INTEGER IPRINT,NMIX,NISOT,NGRP,IDIV
      REAL ZKEFF,RHS(NGRP,NGRP,NISOT+NMIX),LHS(NGRP,NGRP,NISOT+NMIX),
     > FLUX(NGRP,NISOT+NMIX),AFLUX(NGRP,NISOT+NMIX)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) KPLIB
      CHARACTER HSMG*131
      DOUBLE PRECISION SUM
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IMIX
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IHUSED
      REAL, ALLOCATABLE, DIMENSION(:) :: DENS,VOL,TOTAL,ZNUSF,CHI,SIGS,
     > DLK,ALK,V,W
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SCAT
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IHUSED(3,NISOT),IMIX(NISOT))
      ALLOCATE(DENS(NISOT),VOL(NISOT),TOTAL(NGRP),ZNUSF(NGRP),
     > CHI(NGRP),SCAT(NGRP,NGRP))
      ALLOCATE(IPISO(NISOT))
*----
*  FIND ISOTOPE POINTERS IN INPUT MICROLIB
*----
      CALL LIBIPS(IPLIB,NISOT,IPISO)
*----
*  COMPUTE THE RHS AND LHS MATRICES
*----
      CALL LCMGET(IPLIB,'ISOTOPESUSED',IHUSED)
      CALL LCMGET(IPLIB,'ISOTOPESMIX',IMIX)
      CALL LCMGET(IPLIB,'ISOTOPESDENS',DENS)
      CALL LCMGET(IPLIB,'ISOTOPESVOL',VOL)
      CALL LCMGET(IPLIB,'K-EFFECTIVE',ZKEFF)
      IF(IPRINT.GT.1) WRITE(6,'(35H DUO003: EFFECTIVE MULTIPLICATION F,
     > 6HACTOR=,1P,E12.5)') ZKEFF
      CALL XDRSET(RHS,NGRP*NGRP*(NISOT+NMIX),0.0)
      CALL XDRSET(LHS,NGRP*NGRP*(NISOT+NMIX),0.0)
      DO ISOT=1,NISOT
        IF(IPRINT.GT.4) WRITE(6,'(29H DUO003: PROCESSING ISOTOPE '',
     >  3A4,2H''.)') (IHUSED(I0,ISOT),I0=1,3)
        KPLIB=IPISO(ISOT) ! set ISOT-th isotope
        IF(.NOT.C_ASSOCIATED(KPLIB)) THEN
          WRITE(HSMG,'(18H DUO003: ISOTOPE '',3A4,7H'' (ISO=,I8,4H) IS,
     >    31H NOT AVAILABLE IN THE MICROLIB.)') (IHUSED(I0,ISOT),
     >    I0=1,3),ISOT
          CALL XABORT(HSMG)
        ENDIF
        CALL LCMGET(KPLIB,'NWT0',FLUX(1,ISOT))
        CALL LCMLEN(KPLIB,'NWAT0',ILON,ITYLCM)
        IF(ILON.NE.0) THEN
          CALL LCMGET(KPLIB,'NWAT0',AFLUX(1,ISOT))
        ELSE
          CALL XDRSET(AFLUX(1,ISOT),NGRP,1.0)
        ENDIF
        DO IGR=1,NGRP
          FLUX(IGR,ISOT)=FLUX(IGR,ISOT)*VOL(ISOT)
          AFLUX(IGR,ISOT)=AFLUX(IGR,ISOT)*VOL(ISOT)
        ENDDO
        CALL LCMGET(KPLIB,'NTOT0',TOTAL)
        CALL LCMLEN(KPLIB,'NUSIGF',ILON,ITYLCM)
        IF(ILON.GT.0) THEN
          CALL LCMGET(KPLIB,'NUSIGF',ZNUSF)
          CALL LCMGET(KPLIB,'CHI',CHI)
          DO IGR=1,NGRP
            DO JGR=1,NGRP
              LHS(JGR,IGR,ISOT)=LHS(JGR,IGR,ISOT)+DENS(ISOT)*CHI(JGR)*
     >        ZNUSF(IGR)
            ENDDO
          ENDDO
        ENDIF
        ALLOCATE(SIGS(NGRP))
        CALL XDRLGS(KPLIB,-1,IPRINT,0,0,1,NGRP,SIGS,SCAT,ITYPRO)
        DEALLOCATE(SIGS)
        DO IGR=1,NGRP
          DO JGR=1,NGRP
            RHS(JGR,IGR,ISOT)=RHS(JGR,IGR,ISOT)-DENS(ISOT)*SCAT(JGR,IGR)
          ENDDO
          RHS(IGR,IGR,ISOT)=RHS(IGR,IGR,ISOT)+DENS(ISOT)*TOTAL(IGR)
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
            FLUMIX=0.0
            AFLUMI=0.0
            DO ISOT=1,NISOT
              IF(IMIX(ISOT).EQ.IBM) THEN
                FLUMIX=FLUX(JGR,ISOT)
                AFLUMI=AFLUX(JGR,ISOT)
                DO IGR=1,NGRP
                  SUM=SUM+(RHS(JGR,IGR,ISOT)-LHS(JGR,IGR,ISOT)/ZKEFF)*
     >            FLUX(IGR,ISOT)
                ENDDO
              ENDIF
            ENDDO
            RHS(JGR,JGR,NISOT+IBM)=-REAL(SUM)/FLUMIX
            FLUX(JGR,NISOT+IBM)=FLUMIX
            AFLUX(JGR,NISOT+IBM)=AFLUMI
          ENDDO
        ELSE IF(IDIV.EQ.2) THEN
          DO IGR=1,NGRP
            SUM=0.0D0
            FLUMIX=0.0
            AFLUMI=0.0
            DO ISOT=1,NISOT
              IF(IMIX(ISOT).EQ.IBM) THEN
                FLUMIX=FLUX(IGR,ISOT)
                AFLUMI=AFLUX(IGR,ISOT)
                DO JGR=1,NGRP
                  SUM=SUM+(RHS(JGR,IGR,ISOT)-LHS(JGR,IGR,ISOT)/ZKEFF)*
     >            AFLUX(JGR,ISOT)
                ENDDO
              ENDIF
            ENDDO
            RHS(IGR,IGR,NISOT+IBM)=-REAL(SUM)/AFLUMI
            FLUX(IGR,NISOT+IBM)=FLUMIX
            AFLUX(IGR,NISOT+IBM)=AFLUMI
          ENDDO
        ELSE IF(IDIV.EQ.3) THEN
          ALLOCATE(DLK(NGRP),ALK(NGRP))
          DO JGR=1,NGRP
            SUM=0.0D0
            FLUMIX=0.0
            AFLUMI=0.0
            DO ISOT=1,NISOT
              IF(IMIX(ISOT).EQ.IBM) THEN
                FLUMIX=FLUX(JGR,ISOT)
                AFLUMI=AFLUX(JGR,ISOT)
                DO IGR=1,NGRP
                  SUM=SUM+(RHS(JGR,IGR,ISOT)-LHS(JGR,IGR,ISOT)/ZKEFF)*
     >            FLUX(IGR,ISOT)
                ENDDO
              ENDIF
            ENDDO
            DLK(JGR)=REAL(SUM)
            FLUX(JGR,NISOT+IBM)=FLUMIX
          ENDDO
          DO IGR=1,NGRP
            SUM=0.0D0
            FLUMIX=0.0
            AFLUMI=0.0
            DO ISOT=1,NISOT
              IF(IMIX(ISOT).EQ.IBM) THEN
                FLUMIX=FLUX(IGR,ISOT)
                AFLUMI=AFLUX(IGR,ISOT)
                DO JGR=1,NGRP
                  SUM=SUM+(RHS(JGR,IGR,ISOT)-LHS(JGR,IGR,ISOT)/ZKEFF)*
     >            AFLUX(JGR,ISOT)
                ENDDO
              ENDIF
            ENDDO
            ALK(IGR)=REAL(SUM)
            AFLUX(IGR,NISOT+IBM)=AFLUMI
          ENDDO
          ALLOCATE(V(NGRP),W(NGRP))
          CALL DUO005(NGRP,DLK,ALK,FLUX(1,NISOT+IBM),
     >    AFLUX(1,NISOT+IBM),V,W)
          DO IGR=1,NGRP
            DO JGR=1,NGRP
              RHS(IGR,JGR,NISOT+IBM)=RHS(IGR,JGR,NISOT+IBM)-
     >        V(IGR)-W(JGR)
            ENDDO
          ENDDO
          DEALLOCATE(W,V,ALK,DLK)
        ENDIF
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SCAT,CHI,ZNUSF,TOTAL,VOL,DENS)
      DEALLOCATE(IMIX,IHUSED)
      RETURN
      END
