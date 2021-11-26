*DECK SCRISO
      SUBROUTINE SCRISO(IPLIB,NREA,NGRP,NL,NPRC,NOMREA,NWT0,XS,SIGS,
     > SS2D,TAUXFI,LXS,LAMB,CHIRS,BETAR,INVELS,INAME,LSTRD,LPURE,ITRANC,
     > IFISS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Store an isotopic data recovered from a Saphyb into a Microlib.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPLIB   address of the output microlib LCM object
* NREA    number of reactions in the Saphyb object
* NGRP    number of energy groups
* NL      maximum Legendre order (NL=1 is for isotropic scattering)
* NPRC    number of delayed neutron precursor groups
* NOMREA  names of reactions in the Saphyb
* NWT0    average flux
* XS      cross sections per reaction
* SIGS    scattering cross sections
* SS2D    complete scattering matrix
* TAUXFI  interpolated fission rate
* LXS     existence flag of each reaction
* LAMB    decay constants of the delayed neutron precursor groups
* CHIRS   delayed neutron emission spectrums
* BETAR   delayed neutron fractions
* INVELS  group-average of the inverse neutron velocity
* INAME   name of the isotope.
* LSTRD   flag set to .true. if B2=0.0.
* LPURE   =.true. if the interpolation is a pure linear interpolation 
*         with TERP factors.
*
*Parameters: output
* ITRANC  transport correction flag
* IFISS   fission flag
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER NREA,NGRP,NL,NPRC,INAME(2),ITRANC,IFISS
      REAL NWT0(NGRP),XS(NGRP,NREA),SIGS(NGRP,NL),SS2D(NGRP,NGRP,NL),
     > TAUXFI,LAMB(NPRC),CHIRS(NGRP,NPRC),BETAR(NPRC),INVELS(NGRP)
      LOGICAL LXS(NREA),LSTRD,LPURE
      CHARACTER NOMREA(NREA)*12
*----
*  LOCAL VARIABLES
*----
      INTEGER I0, IGFROM, IGMAX, IGMIN, IGR, IGTO, ILEG, IPRC, IREA,
     & NXSCMP
      LOGICAL LDIFF,LZERO
      CHARACTER TEXT12*12
      DOUBLE PRECISION XDRCST
      CHARACTER HCM(0:10)*2,NAMLEG*2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYPRO,NJJ,IJJ
      REAL, ALLOCATABLE, DIMENSION(:) :: STRD,WRK,XSSCMP
      DATA HCM  /'00','01','02','03','04','05','06','07','08','09','10'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(STRD(NGRP))
*----
*  BUILD MICROLIB
*----
      WRITE(TEXT12,'(2A4)') (INAME(I0),I0=1,2)
      CALL LCMPTC(IPLIB,'ALIAS',12,1,TEXT12)
      CALL LCMPUT(IPLIB,'NWT0',NGRP,2,NWT0)
      IF(NPRC.GT.0) THEN
        CALL LCMPUT(IPLIB,'LAMBDA-D',NPRC,2,LAMB)
        CALL LCMPUT(IPLIB,'OVERV',NGRP,2,INVELS)
      ENDIF
      ITRANC=0
      IFISS=0
      LDIFF=.FALSE.
      STRD(:NGRP)=0.0
      DO IREA=1,NREA
        IF(.NOT.LXS(IREA)) CYCLE
        LZERO=.TRUE.
        DO IGR=1,NGRP
          LZERO=LZERO.AND.(XS(IGR,IREA).EQ.0.0)
        ENDDO
        IF(LZERO) CYCLE
        IF(NOMREA(IREA).EQ.'TOTALE') THEN
          IF(LSTRD) THEN
            DO IGR=1,NGRP
              STRD(IGR)=STRD(IGR)+XS(IGR,IREA)
            ENDDO
          ENDIF
          CALL LCMPUT(IPLIB,'NTOT0',NGRP,2,XS(1,IREA))
        ELSE IF(NOMREA(IREA).EQ.'TOTALE P1') THEN
          CALL LCMPUT(IPLIB,'NTOT1',NGRP,2,XS(1,IREA))
        ELSE IF(NOMREA(IREA).EQ.'EXCESS') THEN
*         correct scattering XS with excess XS
          DO IGR=1,NGRP
            SIGS(IGR,1)=SIGS(IGR,1)+XS(IGR,IREA)
          ENDDO
          CALL LCMPUT(IPLIB,'N2N',NGRP,2,XS(1,IREA))
        ELSE IF(NOMREA(IREA).EQ.'FISSION') THEN
          CALL LCMPUT(IPLIB,'NFTOT',NGRP,2,XS(1,IREA))
        ELSE IF(NOMREA(IREA).EQ.'ABSORPTION') THEN
          CALL LCMPUT(IPLIB,'NG',NGRP,2,XS(1,IREA))
        ELSE IF(NOMREA(IREA).EQ.'SPECTRE') THEN
          IF(.NOT.LPURE) THEN
            DO IGR=1,NGRP
              IF(XS(IGR,IREA).NE.0.0) THEN
                XS(IGR,IREA)=XS(IGR,IREA)/TAUXFI
              ENDIF
            ENDDO
          ENDIF
          CALL LCMPUT(IPLIB,'CHI',NGRP,2,XS(1,IREA))
          DO IPRC=1,NPRC
            WRITE(TEXT12,'(A3,I2.2)') 'CHI',IPRC
            CALL LCMPUT(IPLIB,TEXT12,NGRP,2,CHIRS(1,IPRC))
          ENDDO
        ELSE IF(NOMREA(IREA).EQ.'NU*FISSION') THEN
          IFISS=1
          CALL LCMPUT(IPLIB,'NUSIGF',NGRP,2,XS(1,IREA))
          IF(NPRC.GT.0) THEN
            ALLOCATE(WRK(NGRP))
            DO IPRC=1,NPRC
              DO IGR=1,NGRP
                WRK(IGR)=XS(IGR,IREA)*BETAR(IPRC)
              ENDDO
              WRITE(TEXT12,'(A6,I2.2)') 'NUSIGF',IPRC
              CALL LCMPUT(IPLIB,TEXT12,NGRP,2,WRK)
            ENDDO
            DEALLOCATE(WRK)
          ENDIF
        ELSE IF(NOMREA(IREA).EQ.'ENERGIE') THEN
          ALLOCATE(WRK(NGRP))
          DO IGR=1,NGRP
            WRK(IGR)=XS(IGR,IREA)*1.0E6*REAL(XDRCST('eV','J'))
          ENDDO
          CALL LCMPUT(IPLIB,'H-FACTOR',NGRP,2,WRK)
          DEALLOCATE(WRK)
        ELSE IF(NOMREA(IREA).EQ.'SELF') THEN
          CALL LCMPUT(IPLIB,'SIGW00',NGRP,2,XS(1,IREA))
        ELSE IF(NOMREA(IREA).EQ.'TRANSP-CORR') THEN
          ITRANC=2
          IF(LSTRD) THEN
            DO IGR=1,NGRP
              STRD(IGR)=STRD(IGR)-XS(IGR,IREA)
            ENDDO
          ENDIF
          CALL LCMPUT(IPLIB,'TRANC',NGRP,2,XS(1,IREA))
        ELSE IF(NOMREA(IREA).EQ.'FUITES') THEN
          LDIFF=LSTRD
          IF(.NOT.LSTRD) THEN
            DO IGR=1,NGRP
              LDIFF=LDIFF.OR.(XS(IGR,IREA).NE.0.0)
              STRD(IGR)=XS(IGR,IREA)
            ENDDO
          ENDIF
        ELSE IF(NOMREA(IREA).EQ.'DIFFUSION') THEN
          CYCLE
        ELSE IF(NOMREA(IREA).EQ.'TRANSFERT') THEN
          CYCLE
        ELSE
          CALL LCMPUT(IPLIB,NOMREA(IREA),NGRP,2,XS(1,IREA))
        ENDIF
      ENDDO
      IF(LSTRD) THEN
        IF((ITRANC.EQ.0).AND.(NL.GT.1)) THEN
*         Apollo-type transport correction
          DO IGR=1,NGRP
            STRD(IGR)=STRD(IGR)-SIGS(IGR,2)
          ENDDO
        ENDIF
      ELSE
        DO IGR=1,NGRP
          STRD(IGR)=1.0/(3.0*STRD(IGR))
        ENDDO
      ENDIF
      IF((ITRANC.EQ.0).AND.(NL.GT.1)) THEN
*       Apollo-type transport correction
        ITRANC=2
        CALL LCMPUT(IPLIB,'TRANC',NGRP,2,SIGS(1,2))
      ENDIF
      IF(LDIFF) CALL LCMPUT(IPLIB,'STRD',NGRP,2,STRD)
*----
*  SAVE SCATTERING VECTORS AND MATRICES (DO NOT USE XDRLGS TO SAVE CPU
*  TIME)
*----
      ALLOCATE(NJJ(NGRP),IJJ(NGRP),XSSCMP(NGRP*NGRP),ITYPRO(NL))
      DO ILEG=1,NL
        IF(ILEG.LE.11) THEN
          NAMLEG=HCM(ILEG-1)
        ELSE
          WRITE(NAMLEG,'(I2.2)') ILEG-1
        ENDIF
        CALL LCMPUT(IPLIB,'SIGS'//NAMLEG,NGRP,2,SIGS(1,ILEG))
        NXSCMP=0
        DO IGTO=1,NGRP
          IGMIN=IGTO
          IGMAX=IGTO
          DO IGFROM=1,NGRP
            IF(SS2D(IGTO,IGFROM,ILEG).NE.0.0) THEN
              IGMIN=MIN(IGMIN,IGFROM)
              IGMAX=MAX(IGMAX,IGFROM)
            ENDIF
          ENDDO
          IJJ(IGTO)=IGMAX
          NJJ(IGTO)=IGMAX-IGMIN+1
          DO IGFROM=IGMAX,IGMIN,-1
            NXSCMP=NXSCMP+1
            XSSCMP(NXSCMP)=SS2D(IGTO,IGFROM,ILEG)
          ENDDO
        ENDDO
        CALL LCMPUT(IPLIB,'NJJS'//NAMLEG,NGRP,1,NJJ)
        CALL LCMPUT(IPLIB,'IJJS'//NAMLEG,NGRP,1,IJJ)
        CALL LCMPUT(IPLIB,'SCAT'//NAMLEG,NXSCMP,2,XSSCMP)
        ITYPRO(ILEG)=1
      ENDDO
      CALL LCMPUT(IPLIB,'SCAT-SAVED',NL,1,ITYPRO)
      DEALLOCATE(ITYPRO,XSSCMP,IJJ,NJJ)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(STRD)
      RETURN
      END
