*DECK EDIWP1
      SUBROUTINE EDIWP1(IPFLUX,NW,NGROUP,NUN,NREGIO,NDIM,IADJ,NLIN,
     > NFUNL,NGCOND,NMERGE,KEYANI,VOLUME,IGCOND,IMERGE,FLUXES,AFLUXE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Evaluate the PN weighting spectra for an homogenization based on
* spherical harmonic moments of the flux.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPFLUX  pointer to the flux LCM object.
* NW      order of the spherical harmonic expansion for the flux.
* NGROUP  number of energy groups.
* NUN     number of unknowns in flux array.
* NREGIO  number of regions.
* NDIM    number of dimensions.
* IADJ    type of flux weighting:
*         =0: direct flux weighting;
*         =1: direct-adjoint flux weighting.
* NLIN    number of polynomial components in flux.
* NFUNL   number of spherical harmonic components in flux.
* NGCOND  number of merged energy groups.
* NMERGE  number of merged regions.
* KEYANI  position of spherical harmonic components in unknown vector.
* VOLUME  volumes.
* IGCOND  limit condensed groups.
* IMERGE  region merging matrix.
*
*Parameters: input/output
* FLUXES  weighting function for PN fluxes.
* AFLUXE  weighting function for PN adjoint fluxes.
*
*Reference:
* Jean-Francois Vidal et al., APOLLO3 homogenization techniques for
* transport core calculations - application to the ASTRID CFV core,
* Nuclear Engineering and Technology 49 (2017) 1379 - 1387.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPFLUX
      INTEGER    NW,NREGIO,NGROUP,NUN,NDIM,IADJ,NLIN,NFUNL,NGCOND,
     >           NMERGE,KEYANI(NREGIO,NLIN,NFUNL),IGCOND(NGCOND),
     >           IMERGE(NREGIO)
      REAL       VOLUME(NREGIO),FLUXES(NREGIO,NGROUP,NW),
     >           AFLUXE(NREGIO,NGROUP,NW)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPFLUX,JPFLUA
      DOUBLE PRECISION DVOL
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) ::  WORKF,WORKA
      REAL, ALLOCATABLE, DIMENSION(:,:) ::  FDEN,ADEN
      REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  FLUANI,AFLANI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::  SVOL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) ::  DFLUA,DAFLA
*----
*  INITIALIZATION
*----
      IF(NFUNL.EQ.1) CALL XABORT('EDIWP1: ANIS.GE.2 EXPECTED IN TRACKI'
     > //'NG.')
      JPFLUX=LCMGID(IPFLUX,'FLUX')
      ALLOCATE(WORKF(NUN))
      IF(IADJ.EQ.1) THEN
        JPFLUA=LCMGID(IPFLUX,'AFLUX')
        ALLOCATE(WORKA(NUN))
      ENDIF
*----
*  PROCESS TRIVIAL 1D CASE
*----
      IF(NDIM.EQ.1) THEN
        DO IL=1,NW
          DO IGR=1,NGROUP
            IF(IADJ.EQ.0) THEN
              CALL LCMGDL(JPFLUX,IGR,WORKF)
              DO IREG=1,NREGIO
                FLUXES(IREG,IGR,IL)=WORKF(KEYANI(IREG,1,IL+1))
                AFLUXE(IREG,IGR,IL)=1.0
              ENDDO
            ELSE IF(IADJ.EQ.1) THEN
              CALL LCMGDL(JPFLUX,IGR,WORKF)
              CALL LCMGDL(JPFLUA,IGR,WORKA)
              DO IREG=1,NREGIO
                FLUXES(IREG,IGR,IL)=WORKF(KEYANI(IREG,1,IL+1))
                AFLUXE(IREG,IGR,IL)=WORKA(KEYANI(IREG,1,IL+1))
              ENDDO
            ENDIF
            DO IREG=1,NREGIO
              FLUXES(IREG,IGR,IL)=MAX(ABS(FLUXES(IREG,IGR,IL)),1.0E-10)
              AFLUXE(IREG,IGR,IL)=MAX(ABS(AFLUXE(IREG,IGR,IL)),1.0E-10)
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(WORKF)
        IF(IADJ.EQ.1) DEALLOCATE(WORKA)
        RETURN
      ENDIF
*----
*  RECOVER PN MOMENTS OF THE FLUX
*----
      IOF=1
      DO IL=1,NW
        IOF0=IOF
        NTRM=1
        IF(NDIM.EQ.2) THEN
          NTRM=IL+1
        ELSE IF(NDIM.EQ.3) THEN
          NTRM=2*IL+1
        ENDIF
        ALLOCATE(FLUANI(NREGIO,NGROUP,NTRM),AFLANI(NREGIO,NGROUP,NTRM))
        DO IGR=1,NGROUP
          IF(IADJ.EQ.0) THEN
            CALL LCMGDL(JPFLUX,IGR,WORKF)
            DO IREG=1,NREGIO
              IOF=IOF0
              ID=0
              DO IM=-IL,IL
                IF((NDIM.EQ.2).AND.(MOD(IL+IM,2).EQ.1)) CYCLE
                IOF=IOF+1
                IF(IOF.GT.NFUNL) CALL XABORT('EDIWP1: KEYANI OVERFLOW.')
                ID=ID+1
                FLUANI(IREG,IGR,ID)=WORKF(KEYANI(IREG,1,IOF))
                AFLANI(IREG,IGR,ID)=1.0
              ENDDO
            ENDDO
          ELSE IF(IADJ.EQ.1) THEN
            CALL LCMGDL(JPFLUX,IGR,WORKF)
            CALL LCMGDL(JPFLUA,IGR,WORKA)
            DO IREG=1,NREGIO
              IOF=IOF0
              ID=0
              DO IM=-IL,IL
                IF((NDIM.EQ.2).AND.(MOD(IL+IM,2).EQ.1)) CYCLE
                IOF=IOF+1
                IF(IOF.GT.NFUNL) CALL XABORT('EDIWP1: KEYANI OVERFLOW.')
                ID=ID+1
                FLUANI(IREG,IGR,ID)=WORKF(KEYANI(IREG,1,IOF))
                AFLANI(IREG,IGR,ID)=WORKA(KEYANI(IREG,1,IOF))
              ENDDO
            ENDDO
          ENDIF
        ENDDO
*----
*  CONDENSATION AND HOMOGENIZATION OF SPHERICAL HARMONIC MOMENTS
*----
        ALLOCATE(DFLUA(NMERGE,NGCOND,NTRM),DAFLA(NMERGE,NGCOND,NTRM),
     >  SVOL(NMERGE))
        DFLUA(:NMERGE,:NGCOND,:NTRM)=0.0D0
        DAFLA(:NMERGE,:NGCOND,:NTRM)=0.0D0
        IGRFIN=0
        DO IGRC=1,NGCOND
          IGRDEB=IGRFIN+1
          IGRFIN=IGCOND(IGRC)
          DO IGR=IGRDEB,IGRFIN
            SVOL(:NMERGE)=0.0D0
            DO IREG=1,NREGIO
              IRA=IMERGE(IREG)
              IF(IRA.EQ.0) CYCLE
              DVOL=VOLUME(IREG)
              SVOL(IRA)=SVOL(IRA)+DVOL
              DO ID=1,NTRM
                DFLUA(IRA,IGRC,ID)=DFLUA(IRA,IGRC,ID)+
     >                             FLUANI(IREG,IGR,ID)*DVOL
                DAFLA(IRA,IGRC,ID)=DAFLA(IRA,IGRC,ID)+
     >                             AFLANI(IREG,IGR,ID)*DVOL
              ENDDO
            ENDDO
            DO IRA=1,NMERGE
              DO ID=1,NTRM
                DFLUA(IRA,IGRC,ID)=DFLUA(IRA,IGRC,ID)/SVOL(IRA)
                DAFLA(IRA,IGRC,ID)=DAFLA(IRA,IGRC,ID)/SVOL(IRA)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
*----
*  USE APOLLO3 FORMULA
*----
        ALLOCATE(FDEN(NREGIO,NGROUP),ADEN(NREGIO,NGROUP))
        FLUXES(:NREGIO,:NGROUP,IL)=0.0D0
        AFLUXE(:NREGIO,:NGROUP,IL)=0.0D0
        FDEN(:NREGIO,:NGROUP)=0.0D0
        ADEN(:NREGIO,:NGROUP)=0.0D0
        IGRFIN=0
        DO IGRC=1,NGCOND
          IGRDEB=IGRFIN+1
          IGRFIN=IGCOND(IGRC)
          DO IGR=IGRDEB,IGRFIN
            DO IREG=1,NREGIO
              IRA=IMERGE(IREG)
              IF(IRA.EQ.0) CYCLE
              DO ID=1,NTRM
                FLUXES(IREG,IGR,IL)=FLUXES(IREG,IGR,IL)+
     >                    REAL(FLUANI(IREG,IGR,ID)*DFLUA(IRA,IGRC,ID))
                AFLUXE(IREG,IGR,IL)=AFLUXE(IREG,IGR,IL)+
     >                    REAL(AFLANI(IREG,IGR,ID)*DAFLA(IRA,IGRC,ID))
                FDEN(IREG,IGR)=FDEN(IREG,IGR)+REAL(DFLUA(IRA,IGRC,ID))
                ADEN(IREG,IGR)=ADEN(IREG,IGR)+REAL(DAFLA(IRA,IGRC,ID))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO IGR=1,NGROUP
          DO IREG=1,NREGIO
            FLUXES(IREG,IGR,IL)=FLUXES(IREG,IGR,IL)/FDEN(IREG,IGR)
            AFLUXE(IREG,IGR,IL)=AFLUXE(IREG,IGR,IL)/ADEN(IREG,IGR)
            FLUXES(IREG,IGR,IL)=MAX(ABS(FLUXES(IREG,IGR,IL)),1.0E-10)
            AFLUXE(IREG,IGR,IL)=MAX(ABS(AFLUXE(IREG,IGR,IL)),1.0E-10)
          ENDDO
        ENDDO
        DEALLOCATE(ADEN,FDEN,SVOL,DAFLA,DFLUA,AFLANI,FLUANI)
      ENDDO
      DEALLOCATE(WORKF)
      IF(IADJ.EQ.1) DEALLOCATE(WORKA)
      RETURN
      END
