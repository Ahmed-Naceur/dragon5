*DECK DUO006
      SUBROUTINE DUO006(IPLIB,IPRINT,NISOT,NGRP,HREAC,IDIV,RHS,
     > FLUX,AFLUX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Processing one of the two microlibs and return the RHS matrix for
* the single reaction HREAC.
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
* HREAC   character*8 reaction name of the reaction to process.
* IDIV    type of divergence term processing (=0: no processing;
*         =1: direct processing; =2: adjoint processing;
*         =3: direct-adjoint processing).
*
*Parameters: output
* RHS     macroscopic cross-section matrix.
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
      INTEGER IPRINT,NISOT,NGRP,IDIV
      CHARACTER HREAC*8
      REAL RHS(NGRP,NGRP,NISOT),FLUX(NGRP,NISOT),AFLUX(NGRP,NISOT)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) KPLIB
      CHARACTER HSMG*131
      DOUBLE PRECISION SUM
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IHUSED
      REAL, ALLOCATABLE, DIMENSION(:) :: DENS,VOL,VECTOR,CHI,SIGS,DLK,
     > ALK,V,W
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SCAT,RATE
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IHUSED(3,NISOT))
      ALLOCATE(DENS(NISOT),VOL(NISOT),VECTOR(NGRP),SCAT(NGRP,NGRP),
     > CHI(NGRP))
      ALLOCATE(IPISO(NISOT))
*----
*  FIND ISOTOPE POINTERS IN INPUT MICROLIB
*----
      CALL LIBIPS(IPLIB,NISOT,IPISO)
*----
*  COMPUTE THE RHS AND LHS MATRICES
*----
      CALL LCMGET(IPLIB,'ISOTOPESUSED',IHUSED)
      CALL LCMGET(IPLIB,'ISOTOPESDENS',DENS)
      CALL LCMGET(IPLIB,'ISOTOPESVOL',VOL)
      CALL LCMGET(IPLIB,'K-EFFECTIVE',ZKEFF)
      IF(IPRINT.GT.4) WRITE(6,'(35H DUO006: EFFECTIVE MULTIPLICATION F,
     > 6HACTOR=,1P,E12.5)') ZKEFF
      CALL XDRSET(RHS,NGRP*NGRP*NISOT,0.0)
      DO ISOT=1,NISOT
        IF(IPRINT.GT.4) WRITE(6,'(29H DUO006: PROCESSING ISOTOPE '',
     >  3A4,2H''.)') (IHUSED(I0,ISOT),I0=1,3)
        KPLIB=IPISO(ISOT) ! set ISOT-th isotope
        IF(.NOT.C_ASSOCIATED(KPLIB)) THEN
          WRITE(HSMG,'(18H DUO006: ISOTOPE '',3A4,7H'' (ISO=,I8,4H) IS,
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
        CALL LCMLEN(KPLIB,HREAC,ILONG,ITYLCM)
        IF((ILONG.EQ.0).AND.(HREAC.NE.'LEAK')) CYCLE
        IF(HREAC.EQ.'SCAT00') THEN
          ALLOCATE(SIGS(NGRP))
          CALL XDRLGS(KPLIB,-1,IPRINT,0,0,1,NGRP,SIGS,SCAT,ITYPRO)
          DEALLOCATE(SIGS)
          DO IGR=1,NGRP
            DO JGR=1,NGRP
              RHS(JGR,IGR,ISOT)=RHS(JGR,IGR,ISOT)+DENS(ISOT)*
     >        SCAT(JGR,IGR)
            ENDDO
          ENDDO
        ELSE IF((HREAC.EQ.'NUSIGF').OR.(HREAC.EQ.'CHI')) THEN
          CALL LCMGET(KPLIB,'NUSIGF',VECTOR)
          CALL LCMGET(KPLIB,'CHI',CHI)
          DO IGR=1,NGRP
            DO JGR=1,NGRP
              RHS(JGR,IGR,ISOT)=RHS(JGR,IGR,ISOT)+DENS(ISOT)*CHI(JGR)*
     >        VECTOR(IGR)
            ENDDO
          ENDDO
        ELSE IF(HREAC(:3).EQ.'NWT') THEN
          WRITE(HSMG,'(8HDUO006: ,A8,25H IS A FORBIDDEN REACTION.)')
     >    HREAC
          CALL XABORT(HSMG)
        ELSE IF(HREAC.EQ.'LEAK') THEN
          ALLOCATE(RATE(NGRP,NGRP))
          CALL XDRSET(RATE,NGRP*NGRP,0.0)
          CALL LCMLEN(KPLIB,'NUSIGF',ILON,ITYLCM)
          IF(ILON.GT.0) THEN
            CALL LCMGET(KPLIB,'NUSIGF',VECTOR)
            CALL LCMGET(KPLIB,'CHI',CHI)
            DO IGR=1,NGRP
              DO JGR=1,NGRP
                RATE(JGR,IGR)=RATE(JGR,IGR)-DENS(ISOT)*CHI(JGR)*
     >          VECTOR(IGR)/ZKEFF
              ENDDO
            ENDDO
          ENDIF
          ALLOCATE(SIGS(NGRP))
          CALL XDRLGS(KPLIB,-1,IPRINT,0,0,1,NGRP,SIGS,SCAT,ITYPRO)
          DEALLOCATE(SIGS)
          CALL LCMGET(KPLIB,'NTOT0',VECTOR)
          DO IGR=1,NGRP
            DO JGR=1,NGRP
              RATE(JGR,IGR)=RATE(JGR,IGR)-DENS(ISOT)*SCAT(JGR,IGR)
            ENDDO
            RATE(IGR,IGR)=RATE(IGR,IGR)+DENS(ISOT)*VECTOR(IGR)
          ENDDO
          IF(IDIV.EQ.1) THEN
            DO JGR=1,NGRP
              SUM=0.0D0
              DO IGR=1,NGRP
                SUM=SUM+RATE(JGR,IGR)*FLUX(IGR,ISOT)
              ENDDO
              RHS(JGR,JGR,ISOT)=-REAL(SUM)/FLUX(JGR,ISOT)
            ENDDO
          ELSE IF(IDIV.EQ.2) THEN
            DO IGR=1,NGRP
              SUM=0.0D0
              DO JGR=1,NGRP
                SUM=SUM+RATE(JGR,IGR)*AFLUX(JGR,ISOT)
              ENDDO
              RHS(IGR,IGR,ISOT)=-REAL(SUM)/AFLUX(IGR,ISOT)
            ENDDO
          ELSE IF(IDIV.EQ.3) THEN
            ALLOCATE(DLK(NGRP),ALK(NGRP))
            DO JGR=1,NGRP
              SUM=0.0D0
              DO IGR=1,NGRP
                SUM=SUM+RATE(JGR,IGR)*FLUX(IGR,ISOT)
              ENDDO
              DLK(JGR)=REAL(SUM)
            ENDDO
            DO IGR=1,NGRP
              SUM=0.0D0
              DO JGR=1,NGRP
                SUM=SUM+RATE(JGR,IGR)*AFLUX(JGR,ISOT)
              ENDDO
              ALK(IGR)=REAL(SUM)
            ENDDO
            ALLOCATE(V(NGRP),W(NGRP))
            CALL DUO005(NGRP,DLK,ALK,FLUX(1,ISOT),AFLUX(1,ISOT),V,W)
            DO IGR=1,NGRP
              DO JGR=1,NGRP
                RHS(IGR,JGR,ISOT)=RHS(IGR,JGR,ISOT)-V(IGR)-W(JGR)
              ENDDO
            ENDDO
            DEALLOCATE(W,V,ALK,DLK)
          ENDIF
          DEALLOCATE(RATE)
        ELSE
          CALL LCMGET(KPLIB,HREAC,VECTOR)
          DO IGR=1,NGRP
            RHS(IGR,IGR,ISOT)=RHS(IGR,IGR,ISOT)+DENS(ISOT)*VECTOR(IGR)
          ENDDO
        ENDIF
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IPISO)
      DEALLOCATE(CHI,SCAT,VECTOR,VOL,DENS)
      DEALLOCATE(IHUSED)
      RETURN
      END
