*DECK DMASOU
      SUBROUTINE DMASOU(IPRINT,IPDMA,IPMAC,IPFLX,NG,NREG,NMIL,NL,
     1 NDEL,NED,NAMEAD,NUN,NMERGE,NGCOND,NCST,IMERGE,INDGRP,MATCOD,
     2 KEYFLX,VOL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the GPT sources corresponding to the gradient of a macrolib.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPRINT  print parameter
* IPDMA   pointer to the DMA data structure.
* IPMAC   pointer to the macrolib structure.
* IPFLX   pointer to the multigroup flux.
* NG      number of energy groups.
* NREG    number of regions.
* NMIL    number of material mixtures.
* NL      number of Legendre orders.
* NDEL    number of delayed precursors.
* NED     number of extra edit vectors.
* NAMEAD  names of these extra edits.
* NUN     number of unknowns per energy group.
* NMERGE  number of merged regions.
* NGCOND  number of condensed energy groups.
* NCST    number of DMA fixed sources.
* IMERGE  merging indices.
* INDGRP  condensation indices.
* MATCOD  material mixture indices per region.
* KEYFLX  position of averaged fluxes in unknown vector.
* VOL     volumes.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDMA,IPMAC,IPFLX
      INTEGER IPRINT,NG,NREG,NMIL,NL,NDEL,NED,NAMEAD(2,NED),NMERGE,
     1 NGCOND,NCST,IMERGE(NREG),INDGRP(NG),MATCOD(NREG),KEYFLX(NREG)
      REAL VOL(NREG)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSECT=4,EPSMAX=1.0E-7)
      TYPE(C_PTR) JPFLX,JPDMA,KPDMA,JPMAC,KPMAC
      CHARACTER TEXT12*12,TEXB12*12,HSECT(NSECT)*12,CM*2
      DOUBLE PRECISION WW,SUM,FUNC,ZN
*----
*  DATA STATEMENTS
*----
      DATA HSECT/'NTOT0','SIGS00','N2N','N3N'/
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) ::IKEP
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJS00,NJJS00,IPOS00
      REAL, ALLOCATABLE, DIMENSION(:) :: FLUX,SIGT,CHI,EPS,SCAT
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XSSNN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: GAR1,GAR2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) ::GAR3
      REAL, POINTER, DIMENSION(:) :: SUNK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IKEP(NG))
      ALLOCATE(FLUX(NUN),SIGT(NMIL),CHI(NMIL),XSSNN(NMIL,NG,NG),
     1 EPS(NCST))
      ALLOCATE(GAR1(NMERGE,NGCOND),GAR2(NMERGE,NGCOND),
     1 GAR3(NMERGE,NGCOND,NGCOND))
*
      IOF=0
      CALL XDRSET(EPS,NCST,0.0)
      JPFLX=LCMGID(IPFLX,'FLUX')
      JPDMA=LCMLID(IPDMA,'ASOUR',NCST)
*----
*  NWT0 INFORMATION
*----
      SUM=0.0D0
      CALL XDDSET(GAR2,NMERGE*NGCOND,0.0D0)
      DO IG=1,NG
        IGCND=INDGRP(IG)
        CALL LCMGDL(JPFLX,IG,FLUX)
        DO IR=1,NREG
          IF(KEYFLX(IR).EQ.0) CYCLE
          IMERG=IMERGE(IR)
          WW=FLUX(KEYFLX(IR))*VOL(IR)
          SUM=SUM+WW
          IF((IGCND.NE.0).AND.(IMERG.NE.0)) THEN            
            GAR2(IMERG,IGCND)=GAR2(IMERG,IGCND)+WW
          ENDIF
        ENDDO
      ENDDO
      DO IGCND=1,NGCOND
        DO IMERG=1,NMERGE
          IOF=IOF+1
          IF(IOF.GT.NCST) CALL XABORT('DMASOU: NCST OVERFLOW(1).')
          DO IG=1,NG
            IKEP(IG)=LCMARA(NUN)
            CALL C_F_POINTER(IKEP(IG),SUNK,(/ NUN /))
            CALL XDRSET(SUNK,NUN,0.0)
            DO IR=1,NREG
              IF(KEYFLX(IR).EQ.0) CYCLE
              IUNK=KEYFLX(IR)
              FUNC=VOL(IR)*GAR2(IMERG,IGCND)/SUM
              IF((IMERGE(IR).EQ.IMERG).AND.(INDGRP(IG).EQ.IGCND))
     1        THEN
                SUNK(IUNK)=REAL(FUNC/GAR2(IMERG,IGCND))
              ENDIF
              SUNK(IUNK)=REAL(SUNK(IUNK)-FUNC/SUM)
              ZN=SUM/FUNC
              EPS(IOF)=MAX(EPS(IOF),ABS(SUNK(IUNK)*REAL(ZN)))
            ENDDO
          ENDDO
          IF(EPS(IOF).GT.EPSMAX) THEN
            KPDMA=LCMLIL(JPDMA,IOF,NG)
            DO IG=1,NG
              CALL LCMPPL(KPDMA,IG,NUN,2,IKEP(IG))
            ENDDO
          ELSE
            DO IG=1,NG
              CALL LCMDRD(IKEP(IG))
            ENDDO
          ENDIF
        ENDDO
      ENDDO
*----
*  SET OF NSECT BASIC CROSS SECTIONS
*----
      JPMAC=LCMGID(IPMAC,'GROUP')
      DO ISECT=1,NSECT
        TEXT12=HSECT(ISECT)
        CALL XDDSET(GAR1,NMERGE*NGCOND,0.0D0)
        CALL XDDSET(GAR2,NMERGE*NGCOND,0.0D0)
        DO IG=1,NG
          KPMAC=LCMGIL(JPMAC,IG)
          CALL LCMLEN(KPMAC,TEXT12,ILONG,ITYLCM)
          IF(ILONG.GT.0) THEN
            CALL LCMGET(KPMAC,TEXT12,SIGT)
            IGCND=INDGRP(IG)
            CALL LCMGDL(JPFLX,IG,FLUX)
            DO IR=1,NREG
              IF(KEYFLX(IR).EQ.0) CYCLE
              IMERG=IMERGE(IR)
              WW=FLUX(KEYFLX(IR))*VOL(IR)
              IF((IGCND.NE.0).AND.(IMERG.NE.0)) THEN            
                IBM=MATCOD(IR)
                GAR1(IMERG,IGCND)=GAR1(IMERG,IGCND)+WW
                IF(IBM.GT.0) THEN
                  GAR2(IMERG,IGCND)=GAR2(IMERG,IGCND)+SIGT(IBM)*WW
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        DO IGCND=1,NGCOND
          DO IMERG=1,NMERGE
            IOF=IOF+1
            IF(IOF.GT.NCST) CALL XABORT('DMASOU: NCST OVERFLOW(2).')
            IF(GAR2(IMERG,IGCND).NE.0.0) THEN
              DO IG=1,NG
                KPMAC=LCMGIL(JPMAC,IG)
                CALL LCMLEN(KPMAC,TEXT12,ILONG,ITYLCM)
                IKEP(IG)=LCMARA(NUN)
                CALL C_F_POINTER(IKEP(IG),SUNK,(/ NUN /))
                CALL XDRSET(SUNK,NUN,0.0)
                IF(ILONG.GT.0) THEN
                  CALL LCMGET(KPMAC,TEXT12,SIGT)
                  DO IR=1,NREG
                    IF(KEYFLX(IR).EQ.0) CYCLE
                    IF((IMERGE(IR).EQ.IMERG).AND.(INDGRP(IG).EQ.IGCND))
     1              THEN
                      IUNK=KEYFLX(IR)
                      IBM=MATCOD(IR)
                      FUNC=VOL(IR)*GAR2(IMERG,IGCND)/GAR1(IMERG,IGCND)
                      IF(IBM.EQ.0) THEN
                        SUNK(IUNK)=REAL(-FUNC/GAR1(IMERG,IGCND))
                      ELSE
                        SUNK(IUNK)=REAL(FUNC*(SIGT(IBM)/
     1                  GAR2(IMERG,IGCND)-1.0D0/GAR1(IMERG,IGCND)))
                      ENDIF
                      ZN=SUM/FUNC
                      EPS(IOF)=MAX(EPS(IOF),ABS(SUNK(IUNK)*REAL(ZN)))
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
              IF(EPS(IOF).GT.EPSMAX) THEN
                KPDMA=LCMLIL(JPDMA,IOF,NG)
                DO IG=1,NG
                  CALL LCMPPL(KPDMA,IG,NUN,2,IKEP(IG))
                ENDDO
              ELSE
                DO IG=1,NG
                  CALL LCMDRD(IKEP(IG))
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
*----
*  SCATTERING CROSS SECTION INFORMATION
*----
      ALLOCATE(IJJS00(NMIL),NJJS00(NMIL),IPOS00(NMIL))
      DO IL=1,NL
        WRITE(CM,'(I2.2)') IL-1
        CALL XDRSET(XSSNN,NMIL*NG*NG,0.0)
        DO JG=1,NG
          KPMAC=LCMGIL(JPMAC,JG)
          CALL LCMGET(KPMAC,'IJJS'//CM,IJJS00)
          CALL LCMGET(KPMAC,'NJJS'//CM,NJJS00)
          CALL LCMGET(KPMAC,'IPOS'//CM,IPOS00)
          IMAX=0
          DO IBM=1,NMIL
            IMAX=IMAX+NJJS00(IBM)
          ENDDO
          ALLOCATE(SCAT(IMAX))
          CALL LCMGET(KPMAC,'SCAT'//CM,SCAT)
          DO IBM=1,NMIL
            IPOS=IPOS00(IBM)
            IG=IJJS00(IBM)
            IENBR=NJJS00(IBM)
            DO WHILE (IENBR.GE.1) 
              XSSNN(IBM,JG,IG)=SCAT(IPOS) ! JG <-- IG
              IPOS=IPOS+1
              IENBR=IENBR-1
              IG=IG-1
            ENDDO
          ENDDO
          DEALLOCATE(SCAT)
        ENDDO
        CALL XDDSET(GAR1,NMERGE*NGCOND,0.0D0)
        CALL XDDSET(GAR3,NMERGE*NGCOND*NGCOND,0.0D0)
        DO IG=1,NG
          IGCND=INDGRP(IG)
          CALL LCMGDL(JPFLX,IG,FLUX)
          DO IR=1,NREG
            IF(KEYFLX(IR).EQ.0) CYCLE
            IMERG=IMERGE(IR)
            WW=FLUX(KEYFLX(IR))*VOL(IR)
            IF((IGCND.NE.0).AND.(IMERG.NE.0)) THEN            
              IBM=MATCOD(IR)
              GAR1(IMERG,IGCND)=GAR1(IMERG,IGCND)+WW
              IF(IBM.GT.0) THEN
                DO JG=1,NG
                  JGCND=INDGRP(JG)
                  GAR3(IMERG,JGCND,IGCND)=GAR3(IMERG,JGCND,IGCND)
     1                                   +XSSNN(IBM,JG,IG)*WW
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        DO JGCND=1,NGCOND
          DO IGCND=1,NGCOND
            DO IMERG=1,NMERGE
              IOF=IOF+1
              IF(IOF.GT.NCST) CALL XABORT('DMASOU: NCST OVERFLOW(3).')
              IF(GAR3(IMERG,JGCND,IGCND).NE.0.0) THEN
                DO IG=1,NG
                  IKEP(IG)=LCMARA(NUN)
                  CALL C_F_POINTER(IKEP(IG),SUNK,(/ NUN /))
                  CALL XDRSET(SUNK,NUN,0.0)
                  DO JG=1,NG
                    DO IR=1,NREG
                      IF(KEYFLX(IR).EQ.0) CYCLE
                      IF((IMERGE(IR).EQ.IMERG).AND.(INDGRP(IG).EQ.IGCND)
     1                .AND.(INDGRP(JG).EQ.JGCND)) THEN
                        IBM=MATCOD(IR)
                        IF(IBM.NE.0) THEN
                          IUNK=KEYFLX(IR)
                          FUNC=VOL(IR)*GAR3(IMERG,JGCND,IGCND)/
     1                                 GAR1(IMERG,IGCND)
                          SUNK(IUNK)=REAL(SUNK(IUNK)+FUNC*
     1                    XSSNN(IBM,JG,IG)/GAR3(IMERG,JGCND,IGCND))
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDDO
                  DO IR=1,NREG
                    IF(KEYFLX(IR).EQ.0) CYCLE
                    IF((IMERGE(IR).EQ.IMERG).AND.(INDGRP(IG).EQ.IGCND))
     1              THEN
                      IUNK=KEYFLX(IR)
                      FUNC=VOL(IR)*GAR3(IMERG,JGCND,IGCND)/
     1                             GAR1(IMERG,IGCND)
                      SUNK(IUNK)=SUNK(IUNK)-REAL(FUNC/GAR1(IMERG,IGCND))
                      ZN=SUM/FUNC
                      EPS(IOF)=MAX(EPS(IOF),ABS(SUNK(IUNK)*REAL(ZN)))
                    ENDIF
                  ENDDO
                ENDDO
                IF(EPS(IOF).GT.EPSMAX) THEN
                  KPDMA=LCMLIL(JPDMA,IOF,NG)
                  DO IG=1,NG
                    CALL LCMPPL(KPDMA,IG,NUN,2,IKEP(IG))
                  ENDDO
                ELSE
                  DO IG=1,NG
                    CALL LCMDRD(IKEP(IG))
                  ENDDO
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(IPOS00,NJJS00,IJJS00)
*----
*  FISSION INFORMATION
*----
      DO IDEL=1,1+NDEL
        IF(IDEL.EQ.1) THEN
          TEXT12='NUSIGF'
          TEXB12='CHI'
        ELSE
          WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL-1
          WRITE(TEXB12,'(3HCHI,I2.2)') IDEL-1
        ENDIF
        CALL XDDSET(GAR1,NMERGE*NGCOND,0.0D0)
        CALL XDDSET(GAR2,NMERGE*NGCOND,0.0D0)
        CALL XDDSET(GAR3(1,1,1),NMERGE*NGCOND,0.0D0)
        DO IG=1,NG
          KPMAC=LCMGIL(JPMAC,IG)
          CALL LCMLEN(KPMAC,TEXT12,ILONG,ITYLCM)
          IF(ILONG.GT.0) THEN
            CALL LCMGET(KPMAC,TEXT12,SIGT)
            CALL LCMGET(KPMAC,TEXB12,CHI)
            IGCND=INDGRP(IG)
            CALL LCMGDL(JPFLX,IG,FLUX)
            DO IR=1,NREG
              IF(KEYFLX(IR).EQ.0) CYCLE
              IMERG=IMERGE(IR)
              WW=FLUX(KEYFLX(IR))*VOL(IR)
              IF((IGCND.NE.0).AND.(IMERG.NE.0)) THEN            
                IBM=MATCOD(IR)
                GAR1(IMERG,IGCND)=GAR1(IMERG,IGCND)+WW
                IF(IBM.GT.0) THEN
                  GAR2(IMERG,IGCND)=GAR2(IMERG,IGCND)+SIGT(IBM)*WW
                  GAR3(IMERG,IGCND,1)=GAR3(IMERG,IGCND,1)+CHI(IBM)*
     1            SIGT(IBM)*WW
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        DO IGCND=1,NGCOND
          DO IMERG=1,NMERGE
            IOF=IOF+1
            IF(IOF.GT.NCST) CALL XABORT('DMASOU: NCST OVERFLOW(4).')
            IF(GAR2(IMERG,IGCND).NE.0.0) THEN
              DO IG=1,NG
                KPMAC=LCMGIL(JPMAC,IG)
                CALL LCMLEN(KPMAC,TEXT12,ILONG,ITYLCM)
                IKEP(IG)=LCMARA(NUN)
                CALL C_F_POINTER(IKEP(IG),SUNK,(/ NUN /))
                CALL XDRSET(SUNK,NUN,0.0)
                IF(ILONG.GT.0) THEN
                  CALL LCMGET(KPMAC,TEXT12,SIGT)
                  DO IR=1,NREG
                    IF(KEYFLX(IR).EQ.0) CYCLE
                    IF((IMERGE(IR).EQ.IMERG).AND.(INDGRP(IG).EQ.IGCND))
     1              THEN
                      IUNK=KEYFLX(IR)
                      IBM=MATCOD(IR)
                      FUNC=VOL(IR)*GAR2(IMERG,IGCND)/GAR1(IMERG,IGCND)
                      IF(IBM.EQ.0) THEN
                        SUNK(IUNK)=REAL(-FUNC/GAR1(IMERG,IGCND))
                      ELSE
                        SUNK(IUNK)=REAL(FUNC*(SIGT(IBM)/
     1                  GAR2(IMERG,IGCND)-1.0D0/GAR1(IMERG,IGCND)))
                      ENDIF
                      ZN=SUM/FUNC
                      EPS(IOF)=MAX(EPS(IOF),ABS(SUNK(IUNK)*REAL(ZN)))
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
              IF(EPS(IOF).GT.EPSMAX) THEN
                KPDMA=LCMLIL(JPDMA,IOF,NG)
                DO IG=1,NG
                  CALL LCMPPL(KPDMA,IG,NUN,2,IKEP(IG))
                ENDDO
              ELSE
                DO IG=1,NG
                  CALL LCMDRD(IKEP(IG))
                ENDDO
              ENDIF
            ENDIF
            IOF=IOF+1
            IF(IOF.GT.NCST) CALL XABORT('DMASOU: NCST OVERFLOW(5).')
            IF(GAR3(IMERG,IGCND,1).NE.0.0) THEN
              DO IG=1,NG
                KPMAC=LCMGIL(JPMAC,IG)
                CALL LCMLEN(KPMAC,TEXT12,ILONG,ITYLCM)
                IKEP(IG)=LCMARA(NUN)
                CALL C_F_POINTER(IKEP(IG),SUNK,(/ NUN /))
                CALL XDRSET(SUNK,NUN,0.0)
                IF(ILONG.GT.0) THEN
                  CALL LCMGET(KPMAC,TEXT12,SIGT)
                  CALL LCMGET(KPMAC,TEXB12,CHI)
                  DO IR=1,NREG
                    IF(KEYFLX(IR).EQ.0) CYCLE
                    IF((IMERGE(IR).EQ.IMERG).AND.(INDGRP(IG).EQ.IGCND))
     1              THEN
                      IBM=MATCOD(IR)
                      IF(IBM.NE.0) THEN
                        IUNK=KEYFLX(IR)
                        FUNC=VOL(IR)*SIGT(IBM)*GAR3(IMERG,IGCND,1)/
     1                  GAR2(IMERG,IGCND)
                        SUNK(IUNK)=REAL(FUNC*(CHI(IBM)/
     1                  GAR3(IMERG,IGCND,1)-1.0D0/GAR2(IMERG,IGCND)))
                        ZN=SUM/FUNC
                        EPS(IOF)=MAX(EPS(IOF),ABS(SUNK(IUNK)*REAL(ZN)))
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
              IF(EPS(IOF).GT.EPSMAX) THEN
                KPDMA=LCMLIL(JPDMA,IOF,NG)
                DO IG=1,NG
                  CALL LCMPPL(KPDMA,IG,NUN,2,IKEP(IG))
                ENDDO
              ELSE
                DO IG=1,NG
                  CALL LCMDRD(IKEP(IG))
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
*----
*  ADDITIONAL CROSS SECTION INFORMATION
*----
      DO IED=1,NED
        WRITE(TEXT12,'(2A4)') NAMEAD(1,IED),NAMEAD(2,IED)
        CALL XDDSET(GAR1,NMERGE*NGCOND,0.0D0)
        CALL XDDSET(GAR2,NMERGE*NGCOND,0.0D0)
        DO IG=1,NG
          KPMAC=LCMGIL(JPMAC,IG)
          CALL LCMLEN(KPMAC,TEXT12,ILONG,ITYLCM)
          IF(ILONG.GT.0) THEN
            CALL LCMGET(KPMAC,TEXT12,SIGT)
            IGCND=INDGRP(IG)
            CALL LCMGDL(JPFLX,IG,FLUX)
            DO IR=1,NREG
              IF(KEYFLX(IR).EQ.0) CYCLE
              IMERG=IMERGE(IR)
              WW=FLUX(KEYFLX(IR))*VOL(IR)
              IF((IGCND.NE.0).AND.(IMERG.NE.0)) THEN            
                IBM=MATCOD(IR)
                GAR1(IMERG,IGCND)=GAR1(IMERG,IGCND)+WW
                IF(IBM.GT.0) THEN
                  GAR2(IMERG,IGCND)=GAR2(IMERG,IGCND)+SIGT(IBM)*WW
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        DO IGCND=1,NGCOND
          DO IMERG=1,NMERGE
            IOF=IOF+1
            IF(IOF.GT.NCST) CALL XABORT('DMASOU: NCST OVERFLOW(6).')
            IF(GAR2(IMERG,IGCND).NE.0.0) THEN
              DO IG=1,NG
                KPMAC=LCMGIL(JPMAC,IG)
                CALL LCMLEN(KPMAC,TEXT12,ILONG,ITYLCM)
                IKEP(IG)=LCMARA(NUN)
                CALL C_F_POINTER(IKEP(IG),SUNK,(/ NUN /))
                CALL XDRSET(SUNK,NUN,0.0)
                IF(ILONG.GT.0) THEN
                  CALL LCMGET(KPMAC,TEXT12,SIGT)
                  DO IR=1,NREG
                    IF(KEYFLX(IR).EQ.0) CYCLE
                    IF((IMERGE(IR).EQ.IMERG).AND.(INDGRP(IG).EQ.IGCND))
     1              THEN
                      IUNK=KEYFLX(IR)
                      IBM=MATCOD(IR)
                      FUNC=VOL(IR)*GAR2(IMERG,IGCND)/GAR1(IMERG,IGCND)
                      IF(IBM.EQ.0) THEN
                        SUNK(IUNK)=REAL(-FUNC/GAR1(IMERG,IGCND))
                      ELSE
                        SUNK(IUNK)=REAL(FUNC*(SIGT(IBM)/
     1                  GAR2(IMERG,IGCND)-1.0D0/GAR1(IMERG,IGCND)))
                      ENDIF
                      ZN=SUM/FUNC
                      EPS(IOF)=MAX(EPS(IOF),ABS(SUNK(IUNK)*REAL(ZN)))
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
              IF(EPS(IOF).GT.EPSMAX) THEN
                KPDMA=LCMLIL(JPDMA,IOF,NG)
                DO IG=1,NG
                  CALL LCMPPL(KPDMA,IG,NUN,2,IKEP(IG))
                ENDDO
              ELSE
                DO IG=1,NG
                  CALL LCMDRD(IKEP(IG))
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
*----
*  CHECK SOURCE ORTHOGONALITY
*----
      ALLOCATE(SUNK(NUN))
      DO IOF=1,NCST
        CALL LCMLEL(JPDMA,IOF,ILONG,ITYLCM)
        IF(ILONG.NE.0) THEN
          KPDMA=LCMGIL(JPDMA,IOF)
          SUM=0.0D0
          DO IG=1,NG
            CALL LCMGDL(KPDMA,IG,SUNK)
            CALL LCMGDL(JPFLX,IG,FLUX)
            DO IR=1,NREG
              IUNK=KEYFLX(IR)
              IF(IUNK.GT.0) SUM=SUM+SUNK(IUNK)*FLUX(IUNK)
            ENDDO
          ENDDO
          IF(IPRINT.GT.0) THEN
            WRITE(6,'(14H SOURCE INDEX=,I10,14H  DOT PRODUCT=,1P,E11.4,
     1      19H  SOURCE INTENSITY=,E11.4)') IOF,ABS(SUM),EPS(IOF)
          ENDIF
          IF(ABS(SUM).GT.1.0E-5) THEN
            WRITE(TEXT12,'(I10,2X)') IOF
            CALL XABORT('DMASOU: NON ORTHOGONAL SOURCE (IOF='//
     1      TEXT12(:10)//').')
          ENDIF
        ENDIF
      ENDDO
      DEALLOCATE(SUNK)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAR3,GAR2,GAR1)
      DEALLOCATE(EPS,XSSNN,CHI,SIGT,FLUX)
      DEALLOCATE(IKEP)
      RETURN
      END
