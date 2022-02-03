*DECK ACRSXS
      SUBROUTINE ACRSXS(IPAPX,RECNAM,NREA,NGROUP,NISOF,NISOP,NL,INDX,
     1 NOMREA,SIGS,SS2D,XS,LXS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the cross sections of an elementary calculation and single
* mixture in an Apex file.
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
* NGROUP  number of energy groups.
* NISOF   number of fissile isotopes.
* NISOP   number of fission products.
* NL      maximum Legendre order (NL=1 is for isotropic scattering).
* INDX    position of isotopic set in current mixture (=-2: residual
*         set; -1: total set; >0 isotope index).
* NOMREA  names of reactions in the Apex file.
*
*Parameters: output
* SIGS    scattering cross sections.
* SS2D    complete scattering matrix.
* XS      cross sections per reaction.
* LXS     existence flag of each reaction.
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
      INTEGER NREA,NGROUP,NISOF,NISOP,NL,INDX
      REAL SS2D(NGROUP,NGROUP,NL),SIGS(NGROUP,NL),XS(NGROUP,NREA)
      LOGICAL LXS(NREA)
      CHARACTER NOMREA(NREA)*12
*----
*  LOCAL VARIABLES
*----
      INTEGER RANK,TYPE,NBYTE,DIMSR(5),IREA,IOF,IL
      CHARACTER RECNAM2*80,RECNAM3*80
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK1D
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WORK2D
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: WORK3D
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: WORK4D
*----
*  FILL OUTPUT ARRAYS
*----
      SIGS(:NGROUP,:NL)=0.0
      SS2D(:NGROUP,:NGROUP,:NL)=0.0
      XS(:NGROUP,:NREA)=0.0
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
              SIGS(:,:)=WORK2D(:,:)
              DEALLOCATE(WORK2D)
            ELSE
              CALL hdf5_read_data(IPAPX,RECNAM3,WORK3D)
              SIGS(:,:)=WORK3D(:,:,INDX-IOF)
              DEALLOCATE(WORK3D)
            ENDIF
          ELSE IF(NOMREA(IREA).EQ.'SCAT') THEN
            IF(INDX.LT.0) THEN
              CALL hdf5_read_data(IPAPX,RECNAM3,WORK3D)
              SS2D(:,:,:)=WORK3D(:,:,:)
              DEALLOCATE(WORK3D)
            ELSE
              CALL hdf5_read_data(IPAPX,RECNAM3,WORK4D)
              SS2D(:,:,:)=WORK4D(:,:,:,INDX-IOF)
              DEALLOCATE(WORK4D)
            ENDIF
            NL=SIZE(SS2D,3)
            DO IL=2,NL
              SS2D(:,:,IL)=SS2D(:,:,IL)/REAL(2*IL-1)
            ENDDO
          ELSE
            IF(INDX.LT.0) THEN
              CALL hdf5_read_data(IPAPX,RECNAM3,WORK1D)
              XS(:,IREA)=WORK1D(:)
              DEALLOCATE(WORK1D)
            ELSE
              CALL hdf5_read_data(IPAPX,RECNAM3,WORK2D)
              XS(:,IREA)=WORK2D(:,INDX-IOF)
              DEALLOCATE(WORK2D)
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      RETURN
      END
