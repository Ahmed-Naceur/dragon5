*DECK SPHSAP
      SUBROUTINE SPHSAP(IPSAP,IPMAC,ICAL,IMPX,HEQUI,HMASL,NMIL,NGROUP,
     > SPH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Extract a Macrolib corresponding to an elementary calculation in a
* Saphyb.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPSAP   pointer to the Saphyb (L_SAPHYB signature).
* ICAL    index of the elementary calculation being considered.
* IMPX    print parameter (equal to zero for no print).
* HEQUI   keyword of SPH-factor set to be recovered.
* HMASL   keyword of MASL data set to be recovered.
* NMIL    number of mixtures in the elementary calculation.
* NGROUP  number of energy groups in the elementary calculation.
*
*Parameters: output
* IPMAC   pointer to the Macrolib (L_MACROLIB signature).
* SPH     SPH-factor set extracted from the Saphyb.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSAP,IPMAC
      INTEGER ICAL,IMPX,NMIL,NGROUP
      REAL SPH(NMIL,NGROUP)
      CHARACTER HEQUI*4,HMASL*4
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40,MAXREA=25,MAXMAC=2,MAXDIV=3,MAXLOC=10,
     1 MAXMIL=100)
      INTEGER ISTATE(NSTATE),DIMSAP(50)
      INTEGER IPOS(MAXMIL),NJJM(MAXMIL),IJJM(MAXMIL)
      REAL VALDIV(MAXDIV)
      LOGICAL LSTRD,LDIFF,LSPH,LMASL
      CHARACTER TEXT12*12,HSMG*131,NOMREA(MAXREA)*12,CM*2,
     1 IDVAL(MAXDIV)*4,LOCTYP(MAXLOC)*4,LOCKEY(MAXLOC)*4,
     2 NOMMIL(MAXMIL)*20,TEXT8*8,TEXT9*8
      TYPE(C_PTR) JPMAC,KPMAC
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IADRX,LOCAD,ISADRX,LENGDX,
     1 LENGDP,IDATA,IHEDI,TOTM,RESM,ISOTS,NOMISO(:)
      REAL, ALLOCATABLE, DIMENSION(:) :: ENER,XVOLM,FLUXS,RDATA,STR,WRK,
     1 DIFHO,FLUHO,NWT0,SIGS,XS,SS2D,SCAT,GAR,RVALO,SIGSB,SS2DB,XSB,
     2 CONCES,LAMB,CHIRS,BETAR,INVELS,SURF,FMASL
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGS0,SURFLX,WORK
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LXS
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HADF
*----
*  SCRATCH STORAGE ALLOCATION
*   SIGS0    P0 scattering cross sections.
*----
      ALLOCATE(SIGS0(NMIL,NGROUP),FMASL(NMIL))
      CALL XDRSET(FMASL,NMIL,0.0)
*----
*  RECOVER SAPHYB CHARACTERISTICS
*----
      CALL LCMLEN(IPSAP,'DIMSAP',ILENG,ITYLCM)
      IF(ILENG.EQ.0) CALL XABORT('SPHSAP: INVALID SAPHYB.')
      CALL LCMGET(IPSAP,'DIMSAP',DIMSAP)
      IF(NMIL.NE.DIMSAP(7)) THEN
         CALL XABORT('SPHSAP: INVALID VALUE OF NMIL.')
      ELSE IF(NGROUP.NE.DIMSAP(20)) THEN
         CALL XABORT('SPHSAP: INVALID VALUE OF NGROUP.')
      ENDIF
      NREA=DIMSAP(4)   ! number of reactions
      NISO=DIMSAP(5)   ! number of particularized isotopes
      NMAC=DIMSAP(6)   ! number of macroscopic sets
      NPARL=DIMSAP(11) ! number of local variables
      NADRX=DIMSAP(18) ! number of address sets
      NCALS=DIMSAP(19) ! number of elementary calculations in the Saphyb
      NPRC=DIMSAP(31)  ! number of delayed neutron precursor groups
      NISOTS=DIMSAP(32) ! number of isotopes in edition tables
      IF(IMPX.GT.1) THEN
        WRITE(6,'(30H SPHSAP: number of reactions =,I3)') NREA
        WRITE(6,'(44H SPHSAP: number of particularized isotopes =,I4)')
     1  NISO
        WRITE(6,'(37H SPHSAP: number of macroscopic sets =,I2)') NMAC
        WRITE(6,'(29H SPHSAP: number of mixtures =,I5)') NMIL
        WRITE(6,'(36H SPHSAP: number of local variables =,I4)') NPARL
        WRITE(6,'(33H SPHSAP: number of address sets =,I4)') NADRX
        WRITE(6,'(33H SPHSAP: number of calculations =,I5)') NCALS
        WRITE(6,'(34H SPHSAP: number of energy groups =,I4)') NGROUP
        WRITE(6,'(37H SPHSAP: number of precursor groups =,I4)') NPRC
        WRITE(6,'(46H SPHSAP: number of isotopes in output tables =,
     1  I4)') NISOTS
      ENDIF
      IF(NREA.GT.MAXREA) CALL XABORT('SPHSAP: MAXREA OVERFLOW')
      IF(NMAC.GT.MAXMAC) CALL XABORT('SPHSAP: MAXMAC OVERFLOW')
      INDX=NISO+NMAC
      IF(INDX.EQ.0) CALL XABORT('SPHSAP: NO CROSS SECTIONS FOUND.')
*----
*  RECOVER INFORMATION FROM constphysiq DIRECTORY.
*----
      ALLOCATE(ENER(NGROUP+1))
      CALL LCMSIX(IPSAP,'constphysiq',1)
      CALL LCMGET(IPSAP,'ENRGS',ENER)
      CALL LCMSIX(IPSAP,' ',2)
      DO IGR=1,NGROUP+1
        ENER(IGR)=ENER(IGR)/1.0E-6
      ENDDO
      CALL LCMPUT(IPMAC,'ENERGY',NGROUP+1,2,ENER)
      DEALLOCATE(ENER)
*----
*  RECOVER INFORMATION FROM contenu DIRECTORY.
*----
      ALLOCATE(TOTM(NMIL),RESM(NMIL))
      CALL LCMSIX(IPSAP,'contenu',1)
      IF(NREA.GT.0) THEN
        CALL LCMGTC(IPSAP,'NOMREA',12,NREA,NOMREA)
        IF(IMPX.GT.1) THEN
          WRITE(6,'(29H SPHSAP: Available reactions:/(1X,10A13))')
     1    (NOMREA(I),I=1,NREA)
        ENDIF
      ENDIF
      CALL LCMGET(IPSAP,'TOTMAC',TOTM)
      CALL LCMGET(IPSAP,'RESMAC',RESM)
      IF(NISO.GT.0) THEN
        ALLOCATE(NOMISO(NISO*2))
        CALL LCMGET(IPSAP,'NOMISO',NOMISO)
      ENDIF
      CALL LCMSIX(IPSAP,' ',2)
*----
*  RECOVER INFORMATION FROM adresses DIRECTORY.
*----
      NL=0
      IF(NADRX.GT.0) THEN
         ALLOCATE(IADRX((NREA+2)*(NISO+NMAC)*NADRX))
         CALL LCMSIX(IPSAP,'adresses',1)
         CALL LCMGET(IPSAP,'ADRX',IADRX)
         CALL LCMSIX(IPSAP,' ',2)
         DO IAD=1,NADRX
          DO ISO=1,NISO+NMAC
            IOF=(NREA+2)*(NISO+NMAC)*(IAD-1)+(NREA+2)*(ISO-1)+NREA+1
            NL=MAX(NL,IADRX(IOF))
            IOF=(NREA+2)*(NISO+NMAC)*(IAD-1)+(NREA+2)*(ISO-1)+NREA+2
            NL=MAX(NL,IADRX(IOF))
          ENDDO
         ENDDO
      ENDIF
      IF(IMPX.GT.1) THEN
        WRITE(6,'(36H SPHSAP: number of legendre orders =,I4)') NL
      ENDIF
*----
*  RECOVER INFORMATION FROM geom DIRECTORY.
*----
      NSURFD=0
      CALL LCMSIX(IPSAP,'geom',1)
      ALLOCATE(XVOLM(NMIL))
      CALL LCMGET(IPSAP,'XVOLMT',XVOLM)
      CALL LCMGTC(IPSAP,'NOMMIL',20,NMIL,NOMMIL)
      CALL LCMLEN(IPSAP,'outgeom',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
        CALL LCMSIX(IPSAP,'outgeom',1)
        CALL LCMLEN(IPSAP,'SURF',NSURFD,ITYLCM)
        IF(IMPX.GT.1) THEN
          WRITE(6,'(42H SPHSAP: number of discontinuity factors =,I4/)')
     1    NSURFD
        ENDIF
        CALL LCMSIX(IPSAP,' ',2)
      ENDIF
      ALLOCATE(SURFLX(NSURFD,NGROUP),SURF(NSURFD))
      IF(NSURFD.GT.0) THEN
        CALL LCMSIX(IPSAP,'outgeom',1)
        CALL LCMGET(IPSAP,'SURF',SURF)
        CALL LCMSIX(IPSAP,' ',2)
      ENDIF
      CALL LCMSIX(IPSAP,' ',2)
*----
*  RECOVER INFORMATION FROM caldir DIRECTORY.
*----
      WRITE(TEXT12,'(4Hcalc,I8)') ICAL
      CALL LCMLEN(IPSAP,TEXT12,ILENG,ITYLCM)
      IF(ILENG.EQ.0) THEN
         WRITE(HSMG,'(29HSPHSAP: MISSING CALCULATION '',A12,2H''.)')
     1   TEXT12
         CALL XABORT(HSMG)
      ENDIF
      CALL LCMSIX(IPSAP,TEXT12,1)
      CALL LCMSIX(IPSAP,'info',1)
      LSPH=.FALSE.
      LMASL=.FALSE.
      IF(NPARL.GT.0) THEN
        CALL LCMGET(IPSAP,'NLOC',NLOC)
        IF(NLOC.GT.MAXLOC) CALL XABORT('SPHSAP: MAXLOC OVERFLOW')
        CALL LCMGTC(IPSAP,'LOCTYP',4,NLOC,LOCTYP)
        CALL LCMGTC(IPSAP,'LOCKEY',4,NLOC,LOCKEY)
        ALLOCATE(LOCAD(NLOC+1))
        CALL LCMGET(IPSAP,'LOCADR',LOCAD)
        DO ILOC=1,NLOC
          LSPH=LSPH.OR.((LOCTYP(ILOC).EQ.'EQUI').AND.
     1                  (LOCKEY(ILOC).EQ.HEQUI))
          LMASL=LMASL.OR.((LOCTYP(ILOC).EQ.'MASL').AND.
     1                    (LOCKEY(ILOC).EQ.HMASL))
        ENDDO
      ENDIF
      IF((HEQUI.NE.' ').AND.(.NOT.LSPH)) THEN
        WRITE(HSMG,'(46HSPHSAP: UNABLE TO FIND A LOCAL PARAMETER SET O,
     1  25HF TYPE EQUI WITH KEYWORD ,A4,1H.)') HEQUI
        CALL XABORT(HSMG)
      ELSE IF((HMASL.NE.' ').AND.(.NOT.LMASL)) THEN
        WRITE(HSMG,'(46HSPHSAP: UNABLE TO FIND A LOCAL PARAMETER SET O,
     1  25HF TYPE MASL WITH KEYWORD ,A4,1H.)') HMASL
        CALL XABORT(HSMG)
      ENDIF
      ALLOCATE(ISADRX(NMIL),LENGDX(NMIL),LENGDP(NMIL))
      CALL LCMGET(IPSAP,'ISADRX',ISADRX)
      CALL LCMGET(IPSAP,'LENGDX',LENGDX)
      CALL LCMGET(IPSAP,'LENGDP',LENGDP)
      IF(NISOTS.GT.0) THEN
        ALLOCATE(ISOTS(NISOTS*2))
        CALL LCMGET(IPSAP,'ISOTS',ISOTS)
      ENDIF
      CALL LCMSIX(IPSAP,' ',2)
      CALL LCMSIX(IPSAP,'divers',1)
      CALL LCMGET(IPSAP,'NVDIV',NVDIV)
      B2=0.0
      IF(NVDIV.GT.0) THEN
        IF(NVDIV.GT.MAXDIV) CALL XABORT('SPHSAP: MAXDIV OVERFLOW.')
        CALL LCMGTC(IPSAP,'IDVAL',4,NVDIV,IDVAL)
        CALL LCMGET(IPSAP,'VALDIV',VALDIV)
        DO I=1,NVDIV
          IF(IMPX.GT.1) THEN
            WRITE(6,'(9H SPHSAP: ,I3,2X,A,1H=,1P,E13.5)') I,IDVAL(I),
     1      VALDIV(I)
          ENDIF
          IF(IDVAL(I).EQ.'KEFF') THEN
            CALL LCMPUT(IPMAC,'K-EFFECTIVE',1,2,VALDIV(I))
          ELSE IF(IDVAL(I).EQ.'KINF') THEN
            CALL LCMPUT(IPMAC,'K-INFINITY',1,2,VALDIV(I))
          ELSE IF(IDVAL(I).EQ.'B2') THEN
            B2=VALDIV(I)
            CALL LCMPUT(IPMAC,'B2  B1HOM',1,2,VALDIV(I))
          ENDIF
        ENDDO
      ENDIF
      CALL LCMSIX(IPSAP,' ',2)
*----
*  ALLOCATE MACROLIB WORKING ARRAYS.
*----
      ALLOCATE(LXS(NREA),NWT0(NGROUP*NMIL),SIGS(NGROUP*NMIL*NL),
     1 SS2D(NGROUP*NGROUP*NMIL*NL),XS(NGROUP*NMIL*NREA))
      CALL XDRSET(NWT0,NGROUP*NMIL,0.0)
      CALL XDRSET(SIGS,NGROUP*NMIL*NL,0.0)
      CALL XDRSET(SS2D,NGROUP*NGROUP*NMIL*NL,0.0)
      CALL XDRSET(XS,NGROUP*NMIL*NREA,0.0)
      CALL XDLSET(LXS,NREA,.FALSE.)
*----
*  ALLOCATE DELAYED NEUTRON WORKING ARRAYS.
*----
      ALLOCATE(LAMB(NPRC),CHIRS(NGROUP*NMIL*NPRC),BETAR(NMIL*NPRC),
     1 INVELS(NMIL*NGROUP))
      CALL XDRSET(LAMB,NPRC,0.0)
      CALL XDRSET(CHIRS,NGROUP*NMIL*NPRC,0.0)
      CALL XDRSET(BETAR,NMIL*NPRC,0.0)
      CALL XDRSET(INVELS,NMIL*NGROUP,0.0)
      CALL LCMSIX(IPSAP,'divers',1)
      CALL LCMLEN(IPSAP,'NPR',ILONG,ITYLCM)
      IF((NPRC.GT.0).AND.(ILONG.EQ.1)) THEN
        CALL LCMGET(IPSAP,'NPR',NPR)
        IF(NPR.NE.NPRC) CALL XABORT('SPHSAP: NPR INCONSISTENCY(1).')
        CALL LCMGET(IPSAP,'LAMBRS',LAMB)
        DO IBM=1,NMIL
          CALL LCMGET(IPSAP,'CHIRS',CHIRS((IBM-1)*NPRC*NGROUP+1))
          CALL LCMGET(IPSAP,'BETARS',BETAR((IBM-1)*NPRC+1))
          CALL LCMGET(IPSAP,'INVELS',INVELS((IBM-1)*NGROUP+1))
        ENDDO
      ENDIF
      CALL LCMSIX(IPSAP,' ',2)
*----
*  LOOP OVER SAPHYB MIXTURES.
*----
      IF(NADRX.EQ.0) CALL XABORT('SPHSAP: NO ADDRESS SETS AVAILABLE.')
      DO IBM=1,NMIL
        WRITE(TEXT12,'(4Hmili,I8)') IBM
        CALL LCMLEN(IPSAP,TEXT12,ILENG,ITYLCM)
        IF(ILENG.EQ.0) THEN
          WRITE(HSMG,'(29HSPHSAP: MISSING MIXTURE '',A12,2H''.)')
     1    TEXT12
          CALL XABORT(HSMG)
        ENDIF
        CALL LCMSIX(IPSAP,TEXT12,1)
        IMAC=TOTM(IBM)
        IRES=RESM(IBM)
        IAD=ISADRX(IBM)
        NDATAX=LENGDX(IBM)
        NDATAP=LENGDP(IBM)
        ALLOCATE(FLUXS(NGROUP),RDATA(NDATAX),IDATA(NDATAP))
        CALL LCMGET(IPSAP,'FLUXS',FLUXS)
        CALL LCMGET(IPSAP,'RDATAX',RDATA)
        CALL LCMGET(IPSAP,'IDATAP',IDATA)
        DO I=1,NGROUP
          J=(I-1)*NMIL+IBM
          NWT0(J)=NWT0(J)+FLUXS(I)
        ENDDO
        ALLOCATE(SIGSB(NGROUP*NL),SS2DB(NGROUP*NGROUP*NL),
     1  XSB(NGROUP*NREA))
        IF(IMAC.NE.0) THEN
          CALL SPHSXS(NREA,INDX,NADRX,NGROUP,NL,NDATAX,NDATAP,
     1    NISO+IMAC,IAD,IADRX,RDATA,IDATA,NOMREA,SIGSB,SS2DB,
     2    XSB,LXS)
          DO I=1,NGROUP*NL
            J=(I-1)*NMIL+IBM
            SIGS(J)=SIGS(J)+SIGSB(I)
          ENDDO
          DO I=1,NGROUP*NGROUP*NL
            J=(I-1)*NMIL+IBM
            SS2D(J)=SS2D(J)+SS2DB(I)
          ENDDO
          DO I=1,NGROUP*NREA
            J=(I-1)*NMIL+IBM
            XS(J)=XS(J)+XSB(I)
          ENDDO
        ELSE IF(NISO.NE.0) THEN
          IF(NISOTS.EQ.0) CALL XABORT('SPHSAP: MISSING CONCES INFO.')
          ALLOCATE(CONCES(NISOTS))
          CALL LCMGET(IPSAP,'CONCES',CONCES)
          DO ISO=1,NISO
            WRITE(TEXT8,'(2A4)') (NOMISO(2*(ISO-1)+I0),I0=1,2)
            ISOKEP=0
            DO IS2=1,NISOTS
              ISOKEP=IS2
              WRITE(TEXT9,'(2A4)') (ISOTS(2*(IS2-1)+I0),I0=1,2)
              IF(TEXT9.EQ.TEXT8) GO TO 10
            ENDDO
            CYCLE
   10       DEN=CONCES(ISOKEP)
            IF(DEN.NE.0.0) THEN
              CALL SPHSXS(NREA,INDX,NADRX,NGROUP,NL,NDATAX,NDATAP,ISO,
     1        IAD,IADRX,RDATA,IDATA,NOMREA,SIGSB,SS2DB,XSB,LXS)
              DO I=1,NGROUP*NL
                J=(I-1)*NMIL+IBM
                SIGS(J)=SIGS(J)+DEN*SIGSB(I)
              ENDDO
              DO I=1,NGROUP*NGROUP*NL
                J=(I-1)*NMIL+IBM
                SS2D(J)=SS2D(J)+DEN*SS2DB(I)
              ENDDO
              DO I=1,NGROUP*NREA
                J=(I-1)*NMIL+IBM
                XS(J)=XS(J)+DEN*XSB(I)
              ENDDO
            ENDIF
          ENDDO
          DEALLOCATE(CONCES)
          IF(IRES.NE.0) THEN
            CALL SPHSXS(NREA,INDX,NADRX,NGROUP,NL,NDATAX,NDATAP,
     1      NISO+IRES,IAD,IADRX,RDATA,IDATA,NOMREA,SIGSB,SS2DB,
     2      XSB,LXS)
            DO I=1,NGROUP*NL
              J=(I-1)*NMIL+IBM
              SIGS(J)=SIGS(J)+SIGSB(I)
            ENDDO
            DO I=1,NGROUP*NGROUP*NL
              J=(I-1)*NMIL+IBM
              SS2D(J)=SS2D(J)+SS2DB(I)
            ENDDO
            DO I=1,NGROUP*NREA
              J=(I-1)*NMIL+IBM
              XS(J)=XS(J)+XSB(I)
            ENDDO
          ENDIF
        ELSE
          CALL XABORT('SPHSAP: NO MACROSCOPIC SET.')
        ENDIF
        DEALLOCATE(XSB,SS2DB,SIGSB,IDATA,RDATA,FLUXS)
*
        IF(LSPH) THEN
          ALLOCATE(RVALO(LOCAD(NLOC+1)-1))
          CALL LCMGET(IPSAP,'RVALOC',RVALO)
          DO ILOC=1,NLOC
            IF((LOCTYP(ILOC).EQ.'EQUI').AND.(LOCKEY(ILOC).EQ.HEQUI))
     1      THEN
              IF(LOCAD(ILOC+1)-LOCAD(ILOC).NE.NGROUP) THEN
                CALL XABORT('SPHSAP: INVALID NUMBER OF COMPONENTS FOR '
     1          //'SPH FACTORS')
              ENDIF
              DO IGR=1,NGROUP
                SPH(IBM,IGR)=RVALO(LOCAD(ILOC)+IGR-1)
              ENDDO
            ENDIF
          ENDDO
          DEALLOCATE(RVALO)
        ELSE
          SPH(IBM,:NGROUP)=1.0
        ENDIF
        IF(LMASL) THEN
          ALLOCATE(RVALO(LOCAD(NLOC+1)-1))
          CALL LCMGET(IPSAP,'RVALOC',RVALO)
          DO ILOC=1,NLOC
            IF((LOCTYP(ILOC).EQ.'MASL').AND.(LOCKEY(ILOC).EQ.HMASL))
     1      THEN
              IF(LOCAD(ILOC+1)-LOCAD(ILOC).NE.1) THEN
                CALL XABORT('SPHSAP: INVALID NUMBER OF COMPONENTS FOR '
     1          //'MASL')
              ENDIF
              FMASL(IBM)=RVALO(LOCAD(ILOC))
            ENDIF
          ENDDO
          DEALLOCATE(RVALO)
        ENDIF
*
        CALL LCMLEN(IPSAP,'cinetique',ILONG,ITYLCM)
        IF((NPRC.GT.0).AND.(ILONG.NE.0)) THEN
          CALL LCMSIX(IPSAP,'cinetique',1)
          CALL LCMGET(IPSAP,'NPR',NPR)
          IF(NPR.NE.NPRC) CALL XABORT('SPHSAP: NPR INCONSISTENCY(2).')
          CALL LCMGET(IPSAP,'LAMBRS',LAMB)
          CALL LCMGET(IPSAP,'CHIRS',CHIRS((IBM-1)*NPRC*NGROUP+1))
          CALL LCMGET(IPSAP,'BETARS',BETAR((IBM-1)*NPRC+1))
          CALL LCMGET(IPSAP,'INVELS',INVELS((IBM-1)*NGROUP+1))
          CALL LCMSIX(IPSAP,' ',2)
        ENDIF
        CALL LCMSIX(IPSAP,' ',2)
*       END OF LOOP OVER SAPHYB MIXTURES
      ENDDO
      IF(NPARL.GT.0) DEALLOCATE(LOCAD)
*----
*  RECOVER DISCONTINUITY FACTOR INFORMATION
*----
      IDF=0
      IF(NSURFD.GT.0) THEN
        IDF=2
        CALL LCMSIX(IPSAP,'outflx',1)
        CALL LCMGET(IPSAP,'SURFLX',SURFLX)
        CALL LCMSIX(IPSAP,' ',2)
        CALL LCMSIX(IPMAC,'ADF',1)
        CALL LCMPUT(IPMAC,'NTYPE',1,1,NSURFD)
        ALLOCATE(HADF(NSURFD),FLUHO(NGROUP),WORK(NMIL,NGROUP))
        DO I=1,NSURFD
          WRITE(HADF(I),'(3HFD_,I5.5)') I
          DO IGR=1,NGROUP
            WORK(:,IGR)=SURFLX(I,IGR)/SURF(I)
          ENDDO
          CALL LCMPUT(IPMAC,HADF(I),NMIL*NGROUP,2,WORK)
        ENDDO
        CALL LCMPTC(IPMAC,'HADF',8,NSURFD,HADF)
        DEALLOCATE(WORK,FLUHO,HADF)
        CALL LCMSIX(IPMAC,' ',2)
      ENDIF
      DEALLOCATE(SURFLX,SURF)
      CALL LCMSIX(IPSAP,' ',2)
*----
*  IDENTIFY SPECIAL FLUX EDITS
*----
      ALLOCATE(IHEDI(2*NREA))
      NED=0
      DO IREA=1,NREA
        IF(.NOT.LXS(IREA)) CYCLE
        IF(NOMREA(IREA).EQ.'TOTALE') CYCLE
        IF(NOMREA(IREA).EQ.'TOTALE P1') CYCLE
        IF(NOMREA(IREA).EQ.'EXCESS') CYCLE
        IF(NOMREA(IREA).EQ.'FISSION') CYCLE
        IF(NOMREA(IREA).EQ.'SPECTRE') CYCLE
        IF(NOMREA(IREA).EQ.'NU*FISSION') CYCLE
        IF(NOMREA(IREA)(:7).EQ.'ENERGIE') CYCLE
        IF(NOMREA(IREA).EQ.'SELF') CYCLE
        IF(NOMREA(IREA).EQ.'TRANSP-CORR') CYCLE
        IF(NOMREA(IREA).EQ.'FUITES') CYCLE
        IF(NOMREA(IREA).EQ.'DIFFUSION') CYCLE
        IF(NOMREA(IREA).EQ.'TRANSFERT') CYCLE
        NED=NED+1
        READ(NOMREA(IREA),'(2A4)') IHEDI(2*NED-1),IHEDI(2*NED)
      ENDDO
*----
*  STORE MACROLIB.
*----
      CALL LCMPUT(IPMAC,'VOLUME',NMIL,2,XVOLM)
      IF(LMASL) CALL LCMPUT(IPMAC,'MASL',NMIL,2,FMASL)
      IF(NPRC.GT.0) CALL LCMPUT(IPMAC,'LAMBDA-D',NPRC,2,LAMB)
      IFISS=0
      ITRANC=0
      LSTRD=(B2.EQ.0.0)
      LDIFF=.FALSE.
      NW=0
      ALLOCATE(STR(NMIL),WRK(NMIL),DIFHO(NGROUP),FLUHO(NGROUP))
      CALL XDRSET(DIFHO,NGROUP,0.0)
      CALL XDRSET(FLUHO,NGROUP,0.0)
      CALL XDRSET(SIGS0,NMIL*NGROUP,0.0)
      JPMAC=LCMLID(IPMAC,'GROUP',NGROUP)
      DO IGR=1,NGROUP
        CALL XDRSET(STR,NMIL,0.0)
        KPMAC=LCMDIL(JPMAC,IGR)
        IOF=NMIL*(IGR-1)
        CALL LCMPUT(KPMAC,'FLUX-INTG',NMIL,2,NWT0(IOF+1))
        IF(LSPH) CALL LCMPUT(KPMAC,'NSPH',NMIL,2,SPH(1,IGR))
        IF(NPRC.GT.0) THEN
          DO IBM=1,NMIL
            WRK(IBM)=INVELS((IBM-1)*NGROUP+IGR)
          ENDDO
          CALL LCMPUT(KPMAC,'OVERV',NMIL,2,WRK)
        ENDIF
        DO IREA=1,NREA
          IF(.NOT.LXS(IREA)) CYCLE
          IF(NOMREA(IREA).EQ.'TOTALE') THEN
            IOF=NMIL*NGROUP*(IREA-1)+NMIL*(IGR-1)
            IF(LSTRD) THEN
              DO IBM=1,NMIL
                STR(IBM)=STR(IBM)+XS(IOF+IBM)
              ENDDO
            ENDIF
            CALL LCMPUT(KPMAC,'NTOT0',NMIL,2,XS(IOF+1))
          ELSE IF(NOMREA(IREA).EQ.'TOTALE P1') THEN
            NW=1
            IOF=NMIL*NGROUP*(IREA-1)+NMIL*(IGR-1)
            CALL LCMPUT(KPMAC,'NTOT1',NMIL,2,XS(IOF+1))
          ELSE IF(NOMREA(IREA).EQ.'EXCESS') THEN
            IOF=NMIL*NGROUP*(IREA-1)+NMIL*(IGR-1)
*           correct scattering XS with excess XS
            DO IBM=1,NMIL
              SIGS0(IBM,IGR)=SIGS0(IBM,IGR)+XS(IOF+IBM)
            ENDDO
            CALL LCMPUT(KPMAC,'N2N',NMIL,2,XS(IOF+1))
          ELSE IF(NOMREA(IREA).EQ.'FISSION') THEN
            IOF=NMIL*NGROUP*(IREA-1)+NMIL*(IGR-1)
            CALL LCMPUT(KPMAC,'NFTOT',NMIL,2,XS(IOF+1))
          ELSE IF(NOMREA(IREA).EQ.'SPECTRE') THEN
            IOF=NMIL*NGROUP*(IREA-1)+NMIL*(IGR-1)
            CALL LCMPUT(KPMAC,'CHI',NMIL,2,XS(IOF+1))
            DO IPRC=1,NPRC
              DO IBM=1,NMIL
                IOF=(IBM-1)*NGROUP*NPRC+(IPRC-1)*NGROUP+IGR
                WRK(IBM)=CHIRS(IOF)
              ENDDO
              WRITE(TEXT12,'(A3,I2.2)') 'CHI',IPRC
              CALL LCMPUT(KPMAC,TEXT12,NMIL,2,WRK)
            ENDDO
          ELSE IF(NOMREA(IREA).EQ.'NU*FISSION') THEN
            IFISS=1
            IOF=NMIL*NGROUP*(IREA-1)+NMIL*(IGR-1)
            CALL LCMPUT(KPMAC,'NUSIGF',NMIL,2,XS(IOF+1))
            DO IPRC=1,NPRC
              DO IBM=1,NMIL
                WRK(IBM)=XS(IOF+IBM)*BETAR((IBM-1)*NPRC+IPRC)
              ENDDO
              WRITE(TEXT12,'(A6,I2.2)') 'NUSIGF',IPRC
              CALL LCMPUT(KPMAC,TEXT12,NMIL,2,WRK)
            ENDDO
          ELSE IF(NOMREA(IREA).EQ.'ENERGIE') THEN
            IOF=NMIL*NGROUP*(IREA-1)+NMIL*(IGR-1)
            CALL LCMPUT(KPMAC,'H-FACTOR',NMIL,2,XS(IOF+1))
          ELSE IF(NOMREA(IREA).EQ.'SELF') THEN
            IOF=NMIL*NGROUP*(IREA-1)+NMIL*(IGR-1)
            CALL LCMPUT(KPMAC,'SIGW00',NMIL,2,XS(IOF+1))
          ELSE IF(NOMREA(IREA).EQ.'TRANSP-CORR') THEN
            ITRANC=2
            IOF=NMIL*NGROUP*(IREA-1)+NMIL*(IGR-1)
            IF(LSTRD) THEN
              DO IBM=1,NMIL
                STR(IBM)=STR(IBM)-XS(IOF+IBM)
              ENDDO
            ENDIF
            CALL LCMPUT(KPMAC,'TRANC',NMIL,2,XS(IOF+1))
          ELSE IF(NOMREA(IREA).EQ.'FUITES') THEN
            LDIFF=LSTRD
            IF(.NOT.LSTRD) THEN
              IOF=NMIL*NGROUP*(IREA-1)+NMIL*(IGR-1)
              DO IBM=1,NMIL
                LDIFF=LDIFF.OR.(XS(IOF+IBM).NE.0.0)
                STR(IBM)=XS(IOF+IBM)/B2
              ENDDO
            ENDIF
          ELSE IF(NOMREA(IREA).EQ.'DIFFUSION') THEN
            DO IL=1,NL
              WRITE(CM,'(I2.2)') IL-1
              IOF=NMIL*NGROUP*(IL-1)+NMIL*(IGR-1)
              IF(IL.EQ.1) THEN
                DO IBM=1,NMIL
                  SIGS0(IBM,IGR)=SIGS0(IBM,IGR)+SIGS(IOF+IBM)
                ENDDO
              ELSE
                CALL LCMPUT(KPMAC,'SIGS'//CM,NMIL,2,SIGS(IOF+1))
              ENDIF
            ENDDO
          ELSE IF(NOMREA(IREA).EQ.'TRANSFERT') THEN
            ALLOCATE(SCAT(NGROUP*NMIL),GAR(NMIL))
            DO IL=1,NL
              WRITE(CM,'(I2.2)') IL-1
              IPOSDE=0
              DO IBM=1,NMIL
                IPOS(IBM)=IPOSDE+1
                IGMIN=IGR
                IGMAX=IGR
                DO JGR=NGROUP,1,-1
                  IOF=NMIL*NGROUP*NGROUP*(IL-1)+NMIL*NGROUP*(JGR-1)+
     1                NMIL*(IGR-1)
                  IF(SS2D(IOF+1).NE.0.0) THEN
                    IGMIN=MIN(IGMIN,JGR)
                    IGMAX=MAX(IGMAX,JGR)
                  ENDIF
                ENDDO
                IJJM(IBM)=IGMAX
                NJJM(IBM)=IGMAX-IGMIN+1
                DO JGR=IGMAX,IGMIN,-1
                  IPOSDE=IPOSDE+1
                  IOF=NMIL*NGROUP*NGROUP*(IL-1)+NMIL*NGROUP*(JGR-1)+
     1                NMIL*(IGR-1)
                  SCAT(IPOSDE)=SS2D(IOF+IBM)
                ENDDO
                GAR(IBM)=SCAT(IPOS(IBM)+IJJM(IBM)-IGR)
              ENDDO
              CALL LCMPUT(KPMAC,'SCAT'//CM,IPOSDE,2,SCAT)
              CALL LCMPUT(KPMAC,'NJJS'//CM,NMIL,1,NJJM)
              CALL LCMPUT(KPMAC,'IJJS'//CM,NMIL,1,IJJM)
              CALL LCMPUT(KPMAC,'IPOS'//CM,NMIL,1,IPOS)
              CALL LCMPUT(KPMAC,'SIGW'//CM,NMIL,2,GAR)
            ENDDO
            DEALLOCATE(GAR,SCAT)
          ELSE
            IOF=NMIL*NGROUP*(IREA-1)+NMIL*(IGR-1)
            CALL LCMPUT(KPMAC,NOMREA(IREA),NMIL,2,XS(IOF+1))
          ENDIF
        ENDDO
        IF(LSTRD) THEN
          IF((ITRANC.EQ.0).AND.(NL.GT.1)) THEN
*           Apollo-type transport correction
            IOF=NMIL*(NGROUP+IGR-1)
            DO IBM=1,NMIL
              STR(IBM)=STR(IBM)-SIGS(IOF+IBM)
            ENDDO
          ENDIF
          DO IBM=1,NMIL
            STR(IBM)=1.0/(3.0*STR(IBM))
          ENDDO
          LDIFF=.TRUE.
        ENDIF
        IF((ITRANC.EQ.0).AND.(NL.GT.1)) THEN
*         Apollo-type transport correction
          IF(IGR.EQ.NGROUP) ITRANC=2
          IOF=NMIL*(NGROUP+IGR-1)
          CALL LCMPUT(KPMAC,'TRANC',NMIL,2,SIGS(IOF+1))
        ENDIF
        IF(LDIFF) THEN
          CALL LCMPUT(KPMAC,'DIFF',NMIL,2,STR)
          IOF=NMIL*(IGR-1)
          DO IBM=1,NMIL
            FLUHO(IGR)=FLUHO(IGR)+NWT0(IOF+IBM)
            DIFHO(IGR)=DIFHO(IGR)+NWT0(IOF+IBM)*STR(IBM)
          ENDDO
        ENDIF
      ENDDO
      IF(LDIFF) THEN
        DO IGR=1,NGROUP
          DIFHO(IGR)=DIFHO(IGR)/FLUHO(IGR)
        ENDDO
        CALL LCMPUT(IPMAC,'DIFFB1HOM',NGROUP,2,DIFHO)
      ENDIF
      DEALLOCATE(FLUHO,DIFHO,WRK,STR)
*----
*  RELEASE MEMORY
*----
      DEALLOCATE(INVELS,BETAR,CHIRS,LAMB,LXS,XS,SS2D,SIGS,NWT0,LENGDP,
     1 LENGDX,ISADRX,XVOLM)
      IF(NISOTS.GT.0) DEALLOCATE(ISOTS)
      IF(NADRX.GT.0) DEALLOCATE(IADRX)
      IF(NISO.GT.0) DEALLOCATE(NOMISO)
      DEALLOCATE(RESM,TOTM)
*----
*  SAVE SCATTERING P0 INFO
*----
      DO IGR=1,NGROUP
        KPMAC=LCMDIL(JPMAC,IGR)
        CALL LCMPUT(KPMAC,'SIGS00',NMIL,2,SIGS0(1,IGR))
      ENDDO
*----
*  WRITE STATE VECTOR
*----
      TEXT12='L_MACROLIB'
      CALL LCMPTC(IPMAC,'SIGNATURE',12,1,TEXT12)
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NGROUP
      ISTATE(2)=NMIL
      ISTATE(3)=NL ! 1+scattering anisotropy
      ISTATE(4)=IFISS
      ISTATE(5)=NED
      ISTATE(6)=ITRANC
      ISTATE(7)=NPRC
      IF(LDIFF) ISTATE(9)=1
      ISTATE(10)=NW
      ISTATE(12)=IDF
      CALL LCMPUT(IPMAC,'STATE-VECTOR',NSTATE,1,ISTATE)
      IF(NED.GT.0) CALL LCMPUT(IPMAC,'ADDXSNAME-P0',2*NED,3,IHEDI)
      DEALLOCATE(IHEDI)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FMASL,SIGS0)
      RETURN
      END
