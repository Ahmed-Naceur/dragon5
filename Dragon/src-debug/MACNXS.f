*DECK MACNXS
      SUBROUTINE MACNXS(IPLIST,MAXFIS,NGROUP,NBMIX,NIFISS,NANISO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Normalization of macroscopic cross section information.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPLIST  LCM pointer to the macrolib.
* MAXFIS  set to max(1,NIFISS).
* NGROUP  number of energy groups.
* NBMIX   number of mixtures.
* NIFISS  number of fissile isotopes.
* NANISO  maximum Legendre order:
*         =1 isotropic collision;
*         =2 linearly anisotropic collision.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIST
      INTEGER    MAXFIS,NGROUP,NBMIX,NIFISS,NANISO
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPLIST,KPLIST
      CHARACTER  CANISO*2
*----
* ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INGSCT,IFGSCT
      REAL, ALLOCATABLE, DIMENSION(:) :: XSWORK,XSWOR2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: CHWORK
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: SCWORK
*----
*  SCRATCH STORAGE ALLOCATION
*   INGSCT  number of scattering group for cross sections.
*   IFGSCT  first scattering group for cross sections.
*----
      ALLOCATE(INGSCT(NBMIX),IFGSCT(NBMIX))
      ALLOCATE(XSWORK(NBMIX*NGROUP),XSWOR2(NBMIX*NIFISS),
     >         CHWORK(NBMIX,MAXFIS))
      ALLOCATE(SCWORK(NBMIX,NANISO,NGROUP))
*
      DO 100 IMIX=1,NBMIX
        DO 110 IAN=1,NANISO
          DO 120 IG=1,NGROUP
            SCWORK(IMIX,IAN,IG)=0.0D0
 120      CONTINUE
 110    CONTINUE
        DO 130 JFIS=1,NIFISS
          CHWORK(IMIX,JFIS)=0.0
 130    CONTINUE
 100  CONTINUE
      JPLIST=LCMGID(IPLIST,'GROUP')
      DO 140 IGR=1,NGROUP
        KPLIST=LCMGIL(JPLIST,IGR)
*----
*  COMPUTE SUM OF FISSION SPECTRUM.
*----
        CALL LCMLEN(KPLIST,'CHI',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          CALL LCMGET(KPLIST,'CHI',XSWOR2)
          DO 150 IFISS=1,NIFISS
            DO 160 IMAT=1,NBMIX
              CHWORK(IMAT,IFISS)=CHWORK(IMAT,IFISS)
     >                          +XSWOR2((IFISS-1)*NBMIX+IMAT)
 160        CONTINUE
 150      CONTINUE
        ENDIF
*----
*  SUM TRANSFER MATRICES OVER SECONDARY GROUPS.
*----
        DO 170 IANIS=1,NANISO
          WRITE(CANISO,'(I2.2)') IANIS-1
          CALL LCMLEN(KPLIST,'NJJS'//CANISO,ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPLIST,'NJJS'//CANISO,INGSCT)
            CALL LCMGET(KPLIST,'IJJS'//CANISO,IFGSCT)
            CALL LCMGET(KPLIST,'SCAT'//CANISO,XSWORK)
            IPO=0
            DO 180 IMAT=1,NBMIX
              IDG=IFGSCT(IMAT)
              IFG=IDG-INGSCT(IMAT)+1
              DO 190 JGR=IDG,IFG,-1
                IPO=IPO+1
                SCWORK(IMAT,IANIS,JGR)=SCWORK(IMAT,IANIS,JGR)
     >                                +XSWORK(IPO)
 190          CONTINUE
 180        CONTINUE
          ENDIF
 170    CONTINUE
 140  CONTINUE
*----
*  WRITE NORMALIZED X-S ON THE MACROLIB.
*----
      DO 200 IGR=1,NGROUP
        KPLIST=LCMGIL(JPLIST,IGR)
        CALL LCMLEN(KPLIST,'CHI',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          CALL LCMGET(KPLIST,'CHI',XSWOR2)
          DO 210 IFISS=1,NIFISS
            DO 220 IMAT=1,NBMIX
              IF(CHWORK(IMAT,IFISS).GT.0.5) XSWOR2((IFISS-1)*NBMIX+IMAT)
     >        =XSWOR2((IFISS-1)*NBMIX+IMAT)/CHWORK(IMAT,IFISS)
 220        CONTINUE
 210      CONTINUE
          CALL LCMPUT(KPLIST,'CHI',NBMIX*NIFISS,2,XSWOR2)
        ENDIF
        DO 230 IANIS=1,NANISO
          WRITE(CANISO,'(I2.2)') IANIS-1
          CALL LCMLEN(KPLIST,'SIGS'//CANISO,ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            DO 240 IMAT=1,NBMIX
              XSWORK(IMAT)=REAL(SCWORK(IMAT,IANIS,IGR))
 240        CONTINUE
            CALL LCMPUT(KPLIST,'SIGS'//CANISO,NBMIX,2,XSWORK)
          ENDIF
 230    CONTINUE
 200  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SCWORK)
      DEALLOCATE(CHWORK,XSWOR2,XSWORK)
      DEALLOCATE(IFGSCT,INGSCT)
      RETURN
      END
