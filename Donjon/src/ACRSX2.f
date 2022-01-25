*DECK ACRSXS2
      SUBROUTINE ACRSX2(IPAPX,RECNAM,NREA,NGRP,NISOF,NISOP,NL,INDX,
     1 NOMREA,B2APEX,FACT,WEIGHT,SPH,FLUXS,IREAF,LPURE,LXS,XS,SIGS,
     2 SS2D,TAUXFI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the cross sections of an elementary calculation and single
* mixture in an Apex file and perform multiparameter interpolation.
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
* IPAPX   pointer to the Apex file.
* RECNAM  character identification of calculation.
* NREA    number of reactions in the Apex file.
* NGRP    number of energy groups.
* NISOF   number of fissile isotopes.
* NISOP   number of fission products.
* NL      maximum Legendre order (NL=1 is for isotropic scattering).
* INDX    position of isotopic set in current mixture (=-2: residual
*         set; -1: total set; >0 isotope index).
* NOMREA  names of reactions in the Apex file.
* B2APEX  buckling as recovered from the Apex file
* FACT   number density ratio for the isotope
* WEIGHT  interpolation weight
* SPH     SPH factors
* FLUXS   averaged flux
* IREAF   position of 'NUFI' reaction in NOMREA array
* LPURE   =.true. if the interpolation is a pure linear interpolation 
*         with TERP factors.
*
*Parameters: input/output
* LXS     existence flag of each reaction.
* XS      interpolated cross sections per reaction
* SIGS    interpolated scattering cross sections
* SS2D    interpolated scattering matrix
* TAUXFI  interpolated fission rate
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPAPX
      CHARACTER RECNAM*80
      INTEGER NREA,NGRP,NISOF,NISOP,NL,INDX,IREAF
      REAL B2APEX,FACT,WEIGHT,SPH(NGRP),FLUXS(NGRP),SS2D(NGRP,NGRP,NL),
     1 SIGS(NGRP,NL),XS(NGRP,NREA),TAUXFI
      LOGICAL LXS(NREA),LPURE
      CHARACTER NOMREA(NREA)*12
*----
*  LOCAL VARIABLES
*----
      INTEGER RANK,TYPE,NBYTE,DIMSR(5),IREA,IOF,IL,IGR,JGR
      REAL TAUXF,XSECT
      CHARACTER RECNAM2*80,RECNAM3*80
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK1D
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WORK2D,SIGSB,XSB
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: WORK3D,SS2DB
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: WORK4D
*----
*  FILL OUTPUT ARRAYS
*----
      ALLOCATE(SIGSB(NGRP,NL),SS2DB(NGRP,NGRP,NL),XSB(NGRP,NREA))
      SIGSB(:NGRP,:NL)=0.0
      SS2DB(:NGRP,:NGRP,:NL)=0.0
      XSB(:NGRP,:NREA)=0.0
      IOF=0
      IF(INDX.EQ.-2) THEN
        ! residual set
        RECNAM2=TRIM(RECNAM)//"mac/RESIDUAL/"
      ELSE IF(INDX.EQ.-1) THEN
        ! total set
        RECNAM2=TRIM(RECNAM)//"mac/TOTAL/"
      ELSE IF((INDX.GE.1).AND.(INDX.LE.NISOF)) THEN
        ! particularized fissile isotope set
        IOF=0
        RECNAM2=TRIM(RECNAM)//"mic/f.p./"
      ELSE IF((INDX.GE.NISOF+1).AND.(INDX.LE.NISOF+NISOP)) THEN
        ! particularized fission product set
        IOF=NISOF
        RECNAM2=TRIM(RECNAM)//"mic/fiss/"
      ELSE IF(INDX.GE.NISOF+NISOP+1) THEN
        ! particularized stable isotope set
        IOF=NISOF+NISOP
        RECNAM2=TRIM(RECNAM)//"mic/othe/"
      ENDIF
      LXS(:NREA)=.FALSE.
      DO IREA=1,NREA
        RECNAM3=TRIM(RECNAM2)//NOMREA(IREA)
        IF(NOMREA(IREA).EQ.'PROF') CYCLE
        CALL hdf5_info(IPAPX,RECNAM3,RANK,TYPE,NBYTE,DIMSR)
        IF(TYPE.NE.99) THEN
          LXS(IREA)=.TRUE.
          IF(NOMREA(IREA).EQ.'DIFF') THEN
            IF(INDX.LT.0) THEN
              CALL hdf5_read_data(IPAPX,RECNAM3,WORK2D)
              SIGSB(:,:)=WORK2D(:,:)
              DEALLOCATE(WORK2D)
            ELSE
              CALL hdf5_read_data(IPAPX,RECNAM3,WORK3D)
              SIGSB(:,:)=WORK3D(:,:,INDX-IOF)
              DEALLOCATE(WORK3D)
            ENDIF
          ELSE IF(NOMREA(IREA).EQ.'SCAT') THEN
            IF(INDX.LT.0) THEN
              CALL hdf5_read_data(IPAPX,RECNAM3,WORK3D)
              SS2DB(:,:,:)=WORK3D(:,:,:)
              DEALLOCATE(WORK3D)
            ELSE
              CALL hdf5_read_data(IPAPX,RECNAM3,WORK4D)
              SS2DB(:,:,:)=WORK4D(:,:,:,INDX-IOF)
              DEALLOCATE(WORK4D)
            ENDIF
            NL=SIZE(SS2DB,3)
            DO IL=2,NL
              SS2DB(:,:,IL)=SS2DB(:,:,IL)/REAL(2*IL-1)
            ENDDO
          ELSE
            IF(INDX.LT.0) THEN
              CALL hdf5_read_data(IPAPX,RECNAM3,WORK1D)
              XSB(:,IREA)=WORK1D(:)
              DEALLOCATE(WORK1D)
            ELSE
              CALL hdf5_read_data(IPAPX,RECNAM3,WORK2D)
              XSB(:,IREA)=WORK2D(:,INDX-IOF)
              DEALLOCATE(WORK2D)
            ENDIF
          ENDIF
        ENDIF
      ENDDO
*----
*  COMPUTE FISSION RATE FOR AN ELEMENTARY CALCULATION
*----
      TAUXF=0.0
      IF(.NOT.LPURE.AND.(IREAF.GT.0)) THEN
        DO IGR=1,NGRP
          TAUXF=TAUXF+XSB(IGR,IREAF)*FLUXS(IGR)
        ENDDO
        TAUXFI=TAUXFI+WEIGHT*FACT*TAUXF
      ENDIF
*----
*  WEIGHT MICROSCOPIC CROSS SECTION DATA IN AN INTERPOLATED MICROLIB
*----
      DO IGR=1,NGRP
        DO IREA=1,NREA
          IF(.NOT.LXS(IREA)) CYCLE
          IF(LPURE.AND.NOMREA(IREA).EQ.'CHI') THEN
            XS(IGR,IREA)=XS(IGR,IREA)+WEIGHT*XSB(IGR,IREA)
          ELSE IF(NOMREA(IREA).EQ.'CHI') THEN
            IF(IREAF.EQ.0) CALL XABORT('ACRSX2: IREAF=0.')
            XS(IGR,IREA)=XS(IGR,IREA)+WEIGHT*FACT*TAUXF*XSB(IGR,IREA)
          ELSE IF(NOMREA(IREA).EQ.'LEAK') THEN
            IF(B2APEX.NE.0.0) THEN
              XSECT=XSB(IGR,IREA)/B2APEX
              XS(IGR,IREA)=XS(IGR,IREA)+SPH(IGR)*FACT*WEIGHT*XSECT
            ENDIF
          ELSE
            XS(IGR,IREA)=XS(IGR,IREA)+FACT*SPH(IGR)*WEIGHT*XSB(IGR,IREA)
          ENDIF
        ENDDO
        DO IL=1,NL
          IF(MOD(IL,2).EQ.1) THEN
            SIGS(IGR,IL)=SIGS(IGR,IL)+FACT*SPH(IGR)*WEIGHT*SIGSB(IGR,IL)
          ELSE
            DO JGR=1,NGRP
              SIGS(IGR,IL)=SIGS(IGR,IL)+FACT*WEIGHT*SS2DB(JGR,IGR,IL)
     1        /SPH(JGR)
            ENDDO
          ENDIF
        ENDDO
        DO JGR=1,NGRP
          DO IL=1,NL
            IF(MOD(IL,2).EQ.1) THEN
              SS2D(IGR,JGR,IL)=SS2D(IGR,JGR,IL)+FACT*SPH(JGR)*WEIGHT*
     1        SS2DB(IGR,JGR,IL)
            ELSE
              SS2D(IGR,JGR,IL)=SS2D(IGR,JGR,IL)+FACT*WEIGHT*
     1        SS2DB(IGR,JGR,IL)/SPH(IGR)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(XSB,SS2DB,SIGSB)
      RETURN
      END
