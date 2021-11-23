*DECK MCTLIB
      SUBROUTINE MCTLIB(IPLIB,NMIX,NGRP,NL,NFM,NDEL,NED,NAMEAD,LN2N,
     <                  XSTOT,XSS,XSSNN,XSNUSI,XSCHI,XSN2N,XSN3N,XSEDI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover macroscopic cross-section information from the macrolib.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): B. Arsenault
*
*Parameters: input
* IPLIB   pointer to the LIBRARY data structure.
* NMIX    number of mixtures in the geometry.
* NGRP    number of energy groups.
* NL      number of Legendre orders required in the estimations
*         (NL=1 or higher).
* NFM     number of fissile isotopes.
* NDEL    number of delayed precursor groups.
* NED     number of extra edit vectors.
* NAMEAD  names of these extra edits.
* LN2N    N2N cross section recovery flag.
*
*Parameters: output
* XSTOT   total macroscopic cross sections for each mixture and energy
*         group.
* XSS     total scattering cross sections for each mixture and energy
*         group.
* XSSNN   in-group and out-of-group macroscopic transfert cross sections
*         for each mixture and energy group.
* XSNUSI  the values of Nu time the fission cross sections for each
*         isotope per mixture and energy group.
* XSCHI   the values of fission spectrum per isotope per mixture for
*         each energy group.
* XSN2N   N2N macroscopic cross sections for each mixture and energy
*         group.
* XSN3N   N3N macroscopic cross sections for each mixture and energy
*         group.
* XSEDI   extra edit cross sections for each mixture and energy group.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER NMIX,NGRP,NL,NFM,NDEL,NED,NAMEAD(2,NED)
      LOGICAL LN2N
      REAL    XSTOT(NMIX,NGRP),XSS(NMIX,NGRP,NL),XSN2N(NMIX,NGRP),
     <        XSN3N(NMIX,NGRP),XSSNN(NMIX,NGRP,NGRP,NL),
     <        XSNUSI(NMIX,NFM,NGRP,1+NDEL),XSCHI(NMIX,NFM,NGRP,1+NDEL),
     <        XSEDI(NMIX,NGRP,NED)
*----
*  LOCAL VARIABLES 
*----
      TYPE(C_PTR) JPMC,KPMC
      INTEGER   IGROUP,JGROUP,IMAT,IED,IPOS,IEN0,IENBR,ILONG,ITYLCM,
     <          IMAX,IL,IDEL
      DOUBLE PRECISION SUM
      CHARACTER TEXT12*12,CM*2
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJS00,NJJS00,IPOS00
      REAL, ALLOCATABLE, DIMENSION(:) :: SCAT
*----
*  ALLOCATE THE MEMORY THAT IS REQUIRED TO READ THE SCATTERING MATTRICES
*----
      ALLOCATE(IJJS00(NMIX),NJJS00(NMIX),IPOS00(NMIX))
*----
*  PROCESS THE CROSS SECTIONS FOR EACH ENERGY GROUP 
*  THIS IS THE MAIN LOOP
*----
      CALL XDRSET(XSSNN,NMIX*NGRP*NGRP*NL,0.0)
      JPMC = LCMGID(IPLIB,'GROUP')
      DO IGROUP=1,NGRP
        KPMC = LCMGIL(JPMC,IGROUP)
*----
*  READ THE TOTAL MACROSCOPIC CROSS SECTIONS
*----
        CALL LCMGET(KPMC,'NTOT0',XSTOT(1,IGROUP))
*----
*  READ THE TOTAL SCATTERING CROSS SECTIONS AND MATRICES
*----
        DO IL=1,NL
          WRITE(CM,'(I2.2)') IL-1
          CALL LCMGET(KPMC,'SIGS'//CM,XSS(1,IGROUP,IL))
          CALL LCMGET(KPMC,'IJJS'//CM,IJJS00)
          CALL LCMGET(KPMC,'NJJS'//CM,NJJS00)
          CALL LCMGET(KPMC,'IPOS'//CM,IPOS00)
          IMAX=0
          DO IMAT=1,NMIX
            IMAX=IMAX+NJJS00(IMAT)
          ENDDO
          ALLOCATE(SCAT(IMAX))
          CALL LCMGET(KPMC,'SCAT'//CM,SCAT)
          DO IMAT=1,NMIX
            IPOS=IPOS00(IMAT)
            IEN0=IJJS00(IMAT)
            IENBR=NJJS00(IMAT)
            DO WHILE (IENBR.GE.1) 
              XSSNN(IMAT,IGROUP,IEN0,IL)=SCAT(IPOS)
              IPOS=IPOS+1
              IENBR=IENBR-1
              IEN0=IEN0-1
            ENDDO
          ENDDO
          DEALLOCATE(SCAT)
        ENDDO
*----
*  RECOVER THE N2N MACROSCOPIC CROSS SECTIONS
*----
        IF(LN2N) THEN
          CALL LCMLEN(KPMC,'N2N',ILONG,ITYLCM)
          IF(ILONG.GT.0) THEN
            CALL LCMGET(KPMC,'N2N',XSN2N(1,IGROUP))
          ELSE
            CALL XDRSET(XSN2N(1,IGROUP),NMIX,0.0)
          ENDIF
          CALL LCMLEN(KPMC,'N3N',ILONG,ITYLCM)
          IF(ILONG.GT.0) THEN
            CALL LCMGET(KPMC,'N3N',XSN3N(1,IGROUP))
          ELSE
            CALL XDRSET(XSN3N(1,IGROUP),NMIX,0.0)
          ENDIF
          DO IMAT=1,NMIX
           XSS(IMAT,IGROUP,1)=XSS(IMAT,IGROUP,1)-2.0*XSN2N(IMAT,IGROUP)
     1                       -3.0*XSN3N(IMAT,IGROUP)
           IF(XSS(IMAT,IGROUP,1).LT.0.0) CALL XABORT('MCTLIB: BUG1')
           XSS(IMAT,IGROUP,1)=MIN(XSTOT(IMAT,IGROUP),XSS(IMAT,IGROUP,1))
          ENDDO
        ELSE
          CALL XDRSET(XSN2N(1,IGROUP),NMIX,0.0)
          CALL XDRSET(XSN3N(1,IGROUP),NMIX,0.0)
*         N2N CORRECTION IN UPPER ENERGY GROUPS
          DO IMAT=1,NMIX
           IF(XSS(IMAT,IGROUP,1).GT.XSTOT(IMAT,IGROUP)) THEN
              XSN2N(IMAT,IGROUP)=XSS(IMAT,IGROUP,1)-XSTOT(IMAT,IGROUP)
              XSS(IMAT,IGROUP,1)=2.0*XSTOT(IMAT,IGROUP)-
     1        XSS(IMAT,IGROUP,1)
            ENDIF
            IF(XSS(IMAT,IGROUP,1).LT.0.0) CALL XABORT('MCTLIB: BUG2')
          ENDDO
        ENDIF
*----
*  RECOVER FISSION INFORMATION
*----
        IF(NFM.GT.0) THEN
          CALL LCMGET(KPMC,'NUSIGF',XSNUSI(1,1,IGROUP,1))
          CALL LCMGET(KPMC,'CHI',XSCHI(1,1,IGROUP,1))
          DO IDEL=1,NDEL
            WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL
            CALL LCMGET(KPMC,TEXT12,XSNUSI(1,1,IGROUP,1+IDEL))
            WRITE(TEXT12,'(3HCHI,I2.2)') IDEL
            CALL LCMGET(KPMC,TEXT12,XSCHI(1,1,IGROUP,1+IDEL))
          ENDDO
        ENDIF
*----
*  RECOVER SPECIAL EDIT CROSS SECTIONS
*----
        DO IED=1,NED
          WRITE(TEXT12,'(2A4)') NAMEAD(1,IED),NAMEAD(2,IED)
          CALL LCMLEN(KPMC,TEXT12,ILONG,ITYLCM)
          IF(ILONG.GT.0) THEN
            CALL LCMGET(KPMC,TEXT12,XSEDI(1,IGROUP,IED))
          ELSE
            CALL XDRSET(XSEDI(1,IGROUP,IED),NMIX,0.0)
          ENDIF
        ENDDO
      ENDDO
*----
*  RELEASE THE TEMPORARY MEMORY ALLOCATION
*----
      DEALLOCATE(IPOS00,NJJS00,IJJS00)
*----
*  SCATTERING MATRIX NORMALIZATION
*----
      DO IL=1,NL
        DO IMAT=1,NMIX
          DO IGROUP=1,NGRP
            SUM=0.0D0
            DO JGROUP=1,NGRP
              SUM=SUM+XSSNN(IMAT,JGROUP,IGROUP,IL) ! JGROUP <-- IGROUP
            ENDDO
            IF(SUM.NE.0.0) THEN
              DO JGROUP=1,NGRP
               XSSNN(IMAT,JGROUP,IGROUP,IL)=XSSNN(IMAT,JGROUP,IGROUP,IL)
     1         *XSS(IMAT,IGROUP,IL)/REAL(SUM)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
