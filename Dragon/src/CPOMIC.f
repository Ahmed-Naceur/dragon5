*DECK CPOMIC
      SUBROUTINE CPOMIC(IPCPO ,IPEDIT,IPRINT,NGROUP,NMERGE,NBMICR,
     >                  NL    ,IMRG  ,ISOR  ,NPROC ,ISOTMP,IDIMIX,
     >                  INDPRO,ITYPRO,DENCPO,DENTMP,DXSMIC,DSCMIC,
     >                  DMJCPO,DFLUX )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute microscopic cross sections for the Compo file.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau
*
*Parameters: input
* IPCPO   pointer to the compo (L_COMPO signature).
* IPEDIT  pointer to edit information (L_EDIT signature).
* IPRINT  print parameter. Equal to zero for no print.
* NGROUP  number of groups condensed .
* NMERGE  number of regions merged.
* NBMICR  maximum number of isotopes.
* NL      number of Legendre orders.
* IMRG    merge region indices.
* ISOR    Compo isotope number.
* NPROC   number of isotopes.
* ISOTMP  name of isotopes in EDIT.
* IDIMIX  isotopes identifier in each Compo material.
* INDPRO  identifier for xs processing.
* ITYPRO  identifier for xs processed.
* DENCPO  Compo isotopes concentration.
* DENTMP  Dragon isotopes concentration.
*
*Parameters: input/output
* DXSMIC  micro vector xs.
* DSCMIC  micro scat matrix xs.
* DFLUX   flux.
* DMJCPO  fission energy for macro.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT  NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPCPO,IPEDIT
      INTEGER          IPRINT,NGROUP,NMERGE,NBMICR,NL,IMRG,ISOR,NPROC,
     >                 ISOTMP(3,NBMICR),IDIMIX(NMERGE,NBMICR),
     >                 INDPRO(NPROC),ITYPRO(NPROC)
      REAL             DENTMP(NBMICR)
      DOUBLE PRECISION DENCPO(NBMICR),DXSMIC(NGROUP,NPROC),
     >                 DSCMIC(NGROUP,NGROUP,NL),DMJCPO(2,NBMICR),
     >                 DFLUX(NGROUP)
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: HFACT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: XSREC
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XSCAT
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IOUT,NDPROC,IORD
      PARAMETER       (IOUT=6,NDPROC=20)
      TYPE(C_PTR)      KPEDIT
      INTEGER          NBISOE,ISOE,IADDXS,ITC,IXSR,JXSR,KXSR,IL,
     >                 IGR,JGR,ILCMLN,ITYLCM
      CHARACTER        NAMISO*12
      DOUBLE PRECISION DSCPO,DNUFI,FACDEN,FACAXS
*----
*  SCRATCH STORAGE ALLOCATION
*   XSREC   fission energy for macro.
*   XSCAT   compress scattering.
*----
      ALLOCATE(XSREC(NGROUP,NPROC),XSCAT(NGROUP,NGROUP,NL),
     > IPISO(NBMICR))
*----
*  FIND ISOTOPE POINTERS IN INPUT MICROLIB
*----
      CALL LIBIPS(IPEDIT,NBMICR,IPISO)
*----
*  COMPUTE FINAL DENSITY OF ISOTOPE ISOR IN MIXTURE IMRG
*----
      NBISOE=0
      DO 100 ISOE=1,NBMICR
        IF(IDIMIX(IMRG,ISOE).EQ.ISOR) THEN
          NBISOE=NBISOE+1
          DENCPO(ISOR)=DENCPO(ISOR)+DENTMP(ISOE)
        ENDIF
 100  CONTINUE
      IF(DENCPO(ISOR).EQ.0.0D0) THEN
        FACDEN=1.0D0/DBLE(NBISOE)
        IADDXS=0
      ELSE
        FACDEN=1.0D0/DENCPO(ISOR)
        IADDXS=1
      ENDIF
*----
*  SCAN ALL ISOTOPES FOR THIS MERGE REGION AND
*  LOCATE ADD XS FOR THOSE ASSOCIATED WITH ISOR
*----
      DSCPO=0.0D0
      DO 110 ISOE=1,NBMICR
        IF(IDIMIX(IMRG,ISOE).EQ.ISOR) THEN
          IF(IADDXS.EQ.0) THEN
            FACAXS=FACDEN
          ELSE
            FACAXS=FACDEN*DBLE(DENTMP(ISOE))
          ENDIF
*----
*  READ MICRO XS
*----
          WRITE(NAMISO,'(3A4)') (ISOTMP(ITC,ISOE),ITC=1,3)
          IF(IPRINT.GE.10) WRITE(IOUT,6000) NAMISO
          KPEDIT=IPISO(ISOE) ! set ISOE-th isotope
          IORD=1
          CALL XDRLGS(KPEDIT,-1,IPRINT,0,NL-1,IORD,NGROUP,
     >                XSREC(1,NDPROC+1),XSCAT,ITYPRO(NDPROC+1))
          CALL CPOLGX(KPEDIT,-1,IPRINT,IORD,NGROUP,INDPRO,XSREC(1,1),
     >                ITYPRO)
*----
*  ADD MICRO XS TO CPO ISOTOPE
*----
          DO 120 IXSR=1,4
            IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
              DO 121 IGR=1,NGROUP
                DXSMIC(IGR,IXSR)=DXSMIC(IGR,IXSR)
     >            +FACAXS*DBLE(XSREC(IGR,IXSR))
 121          CONTINUE
            ENDIF
 120      CONTINUE
          IXSR=4
          IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
*----
*  COMPUTE TOTAL NEUTRON PRODUCED AND TOTAL ENERGY GENERATED
*----
            CALL LCMLEN(KPEDIT,'H-FACTOR',ILCMLN,ITYLCM)
            IF(ILCMLN.EQ.NGROUP) THEN
               ALLOCATE(HFACT(NGROUP))
               CALL LCMGET(KPEDIT,'H-FACTOR',HFACT)
               DO 122 IGR=1,NGROUP
                 DMJCPO(1,ISOR)=DMJCPO(1,ISOR)+DENTMP(ISOE)
     >               *DFLUX(IGR)*DBLE(HFACT(IGR))*1.0E18
                 DMJCPO(2,ISOR)=DMJCPO(2,ISOR)+DENTMP(ISOE)
     >               *DFLUX(IGR)*DBLE(XSREC(IGR,IXSR))
 122           CONTINUE
               DEALLOCATE(HFACT)
            ENDIF
          ENDIF
          IXSR=5
          IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
            JXSR=3
            DNUFI=0.0D0
            DO 130 IGR=1,NGROUP
              DNUFI=DNUFI+DFLUX(IGR)*DBLE(XSREC(IGR,JXSR))
 130        CONTINUE
            DNUFI=DNUFI*FACAXS
            DSCPO=DSCPO+DNUFI
            DO 140 IGR=1,NGROUP
              DXSMIC(IGR,IXSR)=DXSMIC(IGR,IXSR)
     >          +DNUFI*DBLE(XSREC(IGR,IXSR))
 140        CONTINUE
          ENDIF
          DO 150 IXSR=7,NDPROC
            IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
              DO 151 IGR=1,NGROUP
                DXSMIC(IGR,IXSR)=DXSMIC(IGR,IXSR)
     >            +FACAXS*DBLE(XSREC(IGR,IXSR))
 151          CONTINUE
            ENDIF
 150      CONTINUE
          IL=0
          DO 160 IXSR=NDPROC+1,NDPROC+NL
            IL=IL+1
            IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
              DO 170 IGR=1,NGROUP
                DXSMIC(IGR,IXSR)=DXSMIC(IGR,IXSR)
     >            +FACAXS*DBLE(XSREC(IGR,IXSR))
                DO 171 JGR=1,NGROUP
                  DSCMIC(IGR,JGR,IL)=DSCMIC(IGR,JGR,IL)
     >              +FACAXS*DBLE(XSCAT(IGR,JGR,IL))
 171            CONTINUE
 170          CONTINUE
            ENDIF
 160      CONTINUE
        ENDIF
 110  CONTINUE
*----
*   SAVE AVERAGE XS
*----
      CALL XDRSET(XSREC,NGROUP*NPROC,0.0)
      CALL XDRSET(XSCAT,NGROUP*NGROUP*NL,0.0)
      DO 180 IXSR=1,4
        IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
          DO 181 IGR=1,NGROUP
            XSREC(IGR,IXSR)=REAL(DXSMIC(IGR,IXSR))
 181      CONTINUE
        ENDIF
 180  CONTINUE
      IXSR=5
      IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
        IF(DSCPO.NE.0.0D0) THEN
          DO 190 IGR=1,NGROUP
            XSREC(IGR,IXSR)=REAL(DXSMIC(IGR,IXSR)/DSCPO)
 190      CONTINUE
        ENDIF
      ENDIF
      IXSR=6
      IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
        JXSR=4
        KXSR=3
        DO 191 IGR=1,NGROUP
          IF(DXSMIC(IGR,JXSR).GT.0.0) THEN
            XSREC(IGR,IXSR)=REAL(DXSMIC(IGR,KXSR)/DXSMIC(IGR,JXSR))
          ENDIF
 191    CONTINUE
      ENDIF
      DO 200 IXSR=7,NDPROC
        IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
          DO 201 IGR=1,NGROUP
            XSREC(IGR,IXSR)=REAL(DXSMIC(IGR,IXSR))
 201      CONTINUE
        ENDIF
 200  CONTINUE
      IL=0
      DO 210 IXSR=NDPROC+1,NDPROC+NL
        IL=IL+1
        IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
          DO 220 IGR=1,NGROUP
            XSREC(IGR,IXSR)=REAL(DXSMIC(IGR,IXSR))
            DO 221 JGR=1,NGROUP
              XSCAT(IGR,JGR,IL)=REAL(DSCMIC(IGR,JGR,IL))
 221        CONTINUE
 220      CONTINUE
        ENDIF
 210  CONTINUE
*----
*  SAVE CPO MICRO
*----
      IORD=1
      CALL XDRLGS(IPCPO,1,IPRINT,0,NL-1,IORD,NGROUP,XSREC(1,NDPROC+1),
     >            XSCAT,ITYPRO(NDPROC+1))
      CALL CPOLGX(IPCPO,1,IPRINT,IORD,NGROUP,INDPRO,XSREC(1,1),ITYPRO)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IPISO,XSCAT,XSREC)
      RETURN
*----
*  PRINT FORMAT
*----
 6000 FORMAT(' CPOMIC: INCLUDE ISOTOPE          = ',A12)
      END
