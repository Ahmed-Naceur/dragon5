*DECK ACRCAL
      SUBROUTINE ACRCAL(IPAPX,IPMAC,ICAL,IMPX,NMIL,NGROUP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Extract a Macrolib corresponding to an elementary calculation in an
* Apex file
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPAPX   pointer to the Apex file.
* IPMAC   address of the output Macrolib LCM object.
* ICAL    index of the elementary calculation being considered.
* IMPX    print parameter (equal to zero for no print).
* NMIL    number of mixtures in the elementary calculation.
* NGROUP  number of energy groups in the elementary calculation.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPAPX,IPMAC
      INTEGER ICAL,IMPX,NMIL,NGROUP
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::NSTATE=40
      INTEGER, PARAMETER::MAXFRD=4
      REAL FLOTT, B2, DEN
      INTEGER I, J, K, I0, IBM, IFISS, IGMAX, IGMIN, IGR, IL, IMAC,
     & IPOSDE, IREA, IRES, ISO, ITRANC, JGR, NCALS, NED, NBISO, NISOTS,
     & NL, NBMAC, NPAR, NREA, NSURFD, NISOF, NISOP, NISOS, NBYTE, NLAM,
     & NVP, NPRC, RANK, TYPE, DIMSR(5)
      INTEGER ISTATE(NSTATE)
      LOGICAL LSTRD,LDIFF,LHFACT
      CHARACTER RECNAM*80,TEXT12*12,CM*2,TEXT8*8,HHAD(MAXFRD)*16
      TYPE(C_PTR) JPMAC,KPMAC
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITOTM,IRESM,IPOS,NJJM,IJJM
      REAL, ALLOCATABLE, DIMENSION(:) :: ENER,XVOLM,FLUXS,STR,WRK,
     1 DIFHO,FLUHO,SCAT,GAR,CONCES
      REAL, ALLOCATABLE, DIMENSION(:,:) :: NWT0, EFACT, XSB,SIGSB,SIGS0,
     1 ADF
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XS,SIGS,SS2DB
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: SS2D
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LXS
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HEDI,HADF,NOMISO,
     1 NOMMAC
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: NOMREA
      CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: NOMMIL
*----
*  SCRATCH STORAGE ALLOCATION
*   SIGS0    P0 scattering cross sections.
*----
      ALLOCATE(IPOS(NMIL),NJJM(NMIL),IJJM(NMIL),NOMMIL(NMIL))
      ALLOCATE(SIGS0(NMIL,NGROUP))
*----
*  RECOVER APEX FILE CHARACTERISTICS
*----
      I=0
      CALL ACRTOC(IPAPX,0,NLAM,NREA,NBISO,NBMAC,NMIL,NPAR,NVP,NISOF,
     1 NISOP,NISOS,NCALS,I,NISOTS,NSURFD,NPRC)
      IF(NGROUP.NE.I) CALL XABORT('ACRCAL: INVALID VALUE OF NGROUP.')
      IF(NBISO+NBMAC.EQ.0) CALL XABORT('ACRCAL: NO CROSS SECTIONS.')
*----
*  RECOVER INFORMATION FROM physconst GROUP.
*----
      CALL hdf5_read_data(IPAPX,"/physconst/ENRGS",ENER)
      DO IGR=1,NGROUP+1
        ENER(IGR)=ENER(IGR)/1.0E-6
      ENDDO
      CALL LCMPUT(IPMAC,'ENERGY',NGROUP+1,2,ENER)
      DEALLOCATE(ENER)
*----
*  RECOVER INFORMATION FROM explicit GROUP.
*----
      ALLOCATE(ITOTM(NMIL),IRESM(NMIL))
      ITOTM(:)=0
      IRESM(:)=0
      IF(NREA.GT.0) THEN
        CALL hdf5_read_data(IPAPX,"/explicit/REANAME",NOMREA)
        IF(IMPX.GT.1) THEN
          WRITE(IOUT,'(29H ACRCAL: Available reactions:/(1X,10A13))')
     1    (NOMREA(I),I=1,NREA)
        ENDIF
      ENDIF
      IF(NBISO.GT.0) THEN
        CALL hdf5_read_data(IPAPX,"/explicit/ISONAME",NOMISO)
      ENDIF
      IF(NBMAC.GT.0) THEN
        CALL hdf5_read_data(IPAPX,"/explicit/MACNAME",NOMMAC)
        DO I=1,NBMAC
          IF(NOMMAC(I).EQ.'TOTAL') ITOTM(:)=I
          IF(NOMMAC(I).EQ.'RESIDUAL') IRESM(:)=I
        ENDDO
        DEALLOCATE(NOMMAC)
      ENDIF
*----
*  RECOVER INFORMATION FROM miscellaneous GROUP
*----
      WRITE(RECNAM,'(4Hcalc,I8,15H/miscellaneous/)') ICAL
      CALL hdf5_info(IPAPX,TRIM(RECNAM)//"KEFF",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"KEFF",FLOTT)
        CALL LCMPUT(IPMAC,'K-EFFECTIVE',1,2,FLOTT)
      ENDIF
      CALL hdf5_info(IPAPX,TRIM(RECNAM)//"KINF",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"KINF",FLOTT)
        CALL LCMPUT(IPMAC,'K-INFINITY',1,2,FLOTT)
      ENDIF
      CALL hdf5_info(IPAPX,TRIM(RECNAM)//"B2",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"B2",B2)
        CALL LCMPUT(IPMAC,'B2  B1HOM',1,2,B2)
      ENDIF
      K=0
      CALL hdf5_info(IPAPX,TRIM(RECNAM)//"ADF",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        HHAD(K+1)='ADF'
        K=K+1
      ENDIF
      CALL hdf5_info(IPAPX,TRIM(RECNAM)//"CPDF",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        HHAD(K+1)='CPDF'
        K=K+1
      ENDIF
      CALL hdf5_info(IPAPX,TRIM(RECNAM)//"INTERNAL_ADF",RANK,TYPE,NBYTE,
     1 DIMSR)
      IF(TYPE.NE.99) THEN
        HHAD(K+1)='INTERNAL_ADF'
        K=K+1
      ENDIF
      CALL hdf5_info(IPAPX,TRIM(RECNAM)//"INTERNAL_CPDF",RANK,TYPE,
     1 NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        HHAD(K+1)='INTERNAL_CPDF'
        K=K+1
      ENDIF
      IF(4*K.NE.NSURFD) CALL XABORT('ACRCAL: INVALID ADF COUNT.')
      CALL LCMSIX(IPMAC,'ADF',1)
      CALL LCMPUT(IPMAC,'NTYPE',1,1,NSURFD)
      ALLOCATE(WRK(NGROUP),HADF(NSURFD))
      DO I=1,K
        CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//HHAD(I),ADF)
        DO I0=1,4
          IF(HHAD(I).EQ.'ADF') THEN
            WRITE(TEXT8,'(3HADF,I1)') I0
          ELSE IF(HHAD(I).EQ.'CPDF') THEN
            WRITE(TEXT8,'(4HCPDF,I1)') I0
          ELSE IF(HHAD(I).EQ.'INTERNAL_ADF') THEN
            WRITE(TEXT8,'(6HIN_ADF,I1)') I0
          ELSE IF(HHAD(I).EQ.'INTERNAL_CPDF') THEN
            WRITE(TEXT8,'(7HIN_CPDF,I1)') I0
          ENDIF
          HADF((I-1)*4+I0)=TEXT8
          WRK(:)=ADF(I0,:)
          CALL LCMPUT(IPMAC,TEXT8,NGROUP,2,WRK)
        ENDDO
        DEALLOCATE(ADF)
      ENDDO
      CALL LCMPTC(IPMAC,'HADF',8,NSURFD,HADF)
      DEALLOCATE(HADF,WRK)
      CALL LCMSIX(IPMAC,' ',2)
*----
*  FIND SCATTERING ANISOTROPY.
*----
      WRITE(RECNAM,'(4Hcalc,I8,4H/xs/)') ICAL
      CALL hdf5_info(IPAPX,TRIM(RECNAM)//"mac/TOTAL/DIFF",RANK,
     1 TYPE,NBYTE,DIMSR)
      IF(TYPE.EQ.99) CALL XABORT('ACRCAL: MISSING SCATTERING INFO.')
      NL=DIMSR(2)
      IF(IMPX.GT.1) THEN
        WRITE(IOUT,'(36H ACRCAL: number of Legendre orders =,I4)') NL
      ENDIF
*----
*  ALLOCATE MACROLIB WORKING ARRAYS.
*----
      ALLOCATE(LXS(NREA),NWT0(NMIL,NGROUP),EFACT(NMIL,NGROUP),
     1 XVOLM(NMIL),SIGS(NMIL,NGROUP,NL),SS2D(NMIL,NGROUP,NGROUP,NL),
     2 XS(NMIL,NGROUP,NREA))
      NWT0(:NMIL,:NGROUP)=0.0
      EFACT(:NMIL,:NGROUP)=0.0
      SIGS(:NMIL,:NGROUP,:NL)=0.0
      SS2D(:NMIL,:NGROUP,:NGROUP,:NL)=0.0
      XS(:NMIL,:NGROUP,:NREA)=0.0
      LXS(:NREA)=.FALSE.
*----
*  LOOP OVER APEX MIXTURES.
*----
      DO IBM=1,NMIL
*----
*  RECOVER MIXTURE VOLUMES AND FLUXES.
*----
        CALL hdf5_info(IPAPX,TRIM(RECNAM)//"MEDIA_VOLUME",RANK,
     1  TYPE,NBYTE,DIMSR)
        IF(TYPE.EQ.99) THEN
          XVOLM(IBM)=1.0
          WRITE(IOUT,'(44H ACRCAL: WARNING -- Record MEDIA_VOLUME is m,
     1    42Hissing in the Apex file. Volume set to 1.0)')
        ELSE
          CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"MEDIA_VOLUME",
     1    XVOLM(IBM))
        ENDIF
        CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"FLUX",FLUXS)
        DO I=1,NGROUP
          NWT0(IBM,I)=NWT0(IBM,I)+FLUXS(I)*XVOLM(IBM)
        ENDDO
        DEALLOCATE(FLUXS)
*----
*  RECOVER CROSS SECTIONS.
*----
        IMAC=ITOTM(IBM)
        IRES=IRESM(IBM)
        ALLOCATE(SIGSB(NGROUP,NL),SS2DB(NGROUP,NGROUP,NL),
     1  XSB(NGROUP,NREA))
        IF(IMAC.NE.0) THEN
          CALL ACRSXS(IPAPX,RECNAM,NREA,NGROUP,NISOF,NISOP,NL,-1,NOMREA,
     1    SIGSB,SS2DB,XSB,LXS)
          DO J=1,NL
            DO I=1,NGROUP
              SIGS(IBM,I,J)=SIGS(IBM,I,J)+SIGSB(I,J)
            ENDDO
          ENDDO
          DO K=1,NL
            DO J=1,NGROUP
              DO I=1,NGROUP
                SS2D(IBM,I,J,K)=SS2D(IBM,I,J,K)+SS2DB(I,J,K)
              ENDDO
            ENDDO
          ENDDO
          DO J=1,NREA
            DO I=1,NGROUP
              XS(IBM,I,J)=XS(IBM,I,J)+XSB(I,J)
            ENDDO
          ENDDO
        ELSE IF(NBISO.GT.0) THEN
          CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"mic/CONC",CONCES)
          DO ISO=1,NBISO
            DEN=CONCES(ISO)
            IF(DEN.NE.0.0) THEN
              CALL ACRSXS(IPAPX,RECNAM,NREA,NGROUP,NISOF,NISOP,NL,ISO,
     1        NOMREA,SIGSB,SS2DB,XSB,LXS)
              DO J=1,NL
                DO I=1,NGROUP
                  SIGS(IBM,I,J)=SIGS(IBM,I,J)+DEN*SIGSB(I,J)
                ENDDO
              ENDDO
              DO K=1,NL
                DO J=1,NGROUP
                  DO I=1,NGROUP
                    SS2D(IBM,I,J,K)=SS2D(IBM,I,J,K)+DEN*SS2DB(I,J,K)
                  ENDDO
                ENDDO
              ENDDO
              DO J=1,NREA
                DO I=1,NGROUP
                  XS(IBM,I,J)=XS(IBM,I,J)+DEN*XSB(I,J)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
          DEALLOCATE(CONCES)
          IF(IRES.NE.0) THEN
            CALL ACRSXS(IPAPX,RECNAM,NREA,NGROUP,NISOF,NISOP,NL,-2,
     1      NOMREA,SIGSB,SS2DB,XSB,LXS)
            DO J=1,NL
              DO I=1,NGROUP
                SIGS(IBM,I,J)=SIGS(IBM,I,J)+SIGSB(I,J)
              ENDDO
            ENDDO
            DO K=1,NL
              DO J=1,NGROUP
                DO I=1,NGROUP
                  SS2D(IBM,I,J,K)=SS2D(IBM,I,J,K)+SS2DB(I,J,K)
                ENDDO
              ENDDO
            ENDDO
            DO J=1,NREA
              DO I=1,NGROUP
                XS(IBM,I,J)=XS(IBM,I,J)+XSB(I,J)
              ENDDO
            ENDDO
          ENDIF
        ELSE
          CALL XABORT('ACRCAL: NO MACROSCOPIC SET.')
        ENDIF
        DEALLOCATE(XSB,SS2DB,SIGSB)
*       END OF LOOP OVER APEX MIXTURES
      ENDDO
*----
*  IDENTIFY SPECIAL FLUX EDITS
*----
      ALLOCATE(HEDI(NREA))
      NED=0
      DO IREA=1,NREA
        IF(NOMREA(IREA).EQ.'ABSO') THEN
          NED=NED+1
          HEDI(NED)=NOMREA(IREA)(:8)
          EXIT
        ENDIF
      ENDDO
*----
*  STORE MACROLIB.
*----
      CALL LCMPUT(IPMAC,'VOLUME',NMIL,2,XVOLM)
      IFISS=0
      ITRANC=0
      LSTRD=(B2.EQ.0.0)
      LDIFF=.FALSE.
      LHFACT=.FALSE.
      ALLOCATE(STR(NMIL),WRK(NMIL),DIFHO(NGROUP),FLUHO(NGROUP))
      DIFHO(:NGROUP)=0.0
      FLUHO(:NGROUP)=0.0
      SIGS0(:NMIL,:NGROUP)=0.0
      JPMAC=LCMLID(IPMAC,'GROUP',NGROUP)
      DO IGR=1,NGROUP
        STR(:NMIL)=0.0
        KPMAC=LCMDIL(JPMAC,IGR)
        CALL LCMPUT(KPMAC,'FLUX-INTG',NMIL,2,NWT0(1,IGR))
        DO IREA=1,NREA
          IF(.NOT.LXS(IREA)) CYCLE
          IF(NOMREA(IREA).EQ.'TOTA') THEN
            IF(LSTRD) THEN
              DO IBM=1,NMIL
                STR(IBM)=STR(IBM)+XS(IBM,IGR,IREA)
              ENDDO
            ENDIF
            CALL LCMPUT(KPMAC,'NTOT0',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'N2N') THEN
*           correct scattering XS with excess XS
            SIGS0(:,IGR)=SIGS0(:,IGR)+XS(:,IGR,IREA)
            CALL LCMPUT(KPMAC,'N2N',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'FISS') THEN
            CALL LCMPUT(KPMAC,'NFTOT',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'CHI') THEN
            CALL LCMPUT(KPMAC,'CHI',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'NUFI') THEN
            IFISS=1
            CALL LCMPUT(KPMAC,'NUSIGF',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'ENER') THEN
            LHFACT=.TRUE.
            EFACT(:,IGR)=EFACT(:,IGR)+XS(:,IGR,IREA)
          ELSE IF(NOMREA(IREA).EQ.'KAFI') THEN
            LHFACT=.TRUE.
            EFACT(:,IGR)=EFACT(:,IGR)+XS(:,IGR,IREA)
          ELSE IF(NOMREA(IREA).EQ.'EGAM') THEN
            LHFACT=.TRUE.
            EFACT(:,IGR)=EFACT(:,IGR)+XS(:,IGR,IREA)
          ELSE IF(NOMREA(IREA).EQ.'LEAK') THEN
            LDIFF=LSTRD
            IF(.NOT.LSTRD) THEN
              DO IBM=1,NMIL
                LDIFF=LDIFF.OR.(XS(IBM,IGR,IREA).NE.0.0)
                STR(IBM)=XS(IBM,IGR,IREA)/B2
              ENDDO
            ENDIF
          ELSE IF(NOMREA(IREA).EQ.'DIFF') THEN
            DO IL=1,NL
              WRITE(CM,'(I2.2)') IL-1
              IF(IL.EQ.1) THEN
                DO IBM=1,NMIL
                  SIGS0(IBM,IGR)=SIGS0(IBM,IGR)+SIGS(IBM,IGR,IL)
                ENDDO
              ELSE
                CALL LCMPUT(KPMAC,'SIGS'//CM,NMIL,2,SIGS(1,IGR,IL))
              ENDIF
            ENDDO
          ELSE IF(NOMREA(IREA).EQ.'SCAT') THEN
            ALLOCATE(SCAT(NGROUP*NMIL),GAR(NMIL))
            DO IL=1,NL
              WRITE(CM,'(I2.2)') IL-1
              IPOSDE=0
              DO IBM=1,NMIL
                IPOS(IBM)=IPOSDE+1
                IGMIN=IGR
                IGMAX=IGR
                DO JGR=NGROUP,1,-1
                  IF(SS2D(IBM,IGR,JGR,IL).NE.0.0) THEN
                    IGMIN=MIN(IGMIN,JGR)
                    IGMAX=MAX(IGMAX,JGR)
                  ENDIF
                ENDDO
                IJJM(IBM)=IGMAX
                NJJM(IBM)=IGMAX-IGMIN+1
                DO JGR=IGMAX,IGMIN,-1
                  IPOSDE=IPOSDE+1
                  SCAT(IPOSDE)=SS2D(IBM,IGR,JGR,IL)
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
            CALL LCMPUT(KPMAC,NOMREA(IREA),NMIL,2,XS(1,IGR,IREA))
          ENDIF
        ENDDO
        IF(LSTRD) THEN
          IF((ITRANC.EQ.0).AND.(NL.GT.1)) THEN
*           Apollo-type transport correction
            DO IBM=1,NMIL
              STR(IBM)=STR(IBM)-SIGS(IBM,IGR,2)
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
          CALL LCMPUT(KPMAC,'TRANC',NMIL,2,SIGS(1,IGR,2))
        ENDIF
        IF(LDIFF) THEN
          CALL LCMPUT(KPMAC,'DIFF',NMIL,2,STR)
          DO IBM=1,NMIL
            FLUHO(IGR)=FLUHO(IGR)+NWT0(IBM,IGR)
            DIFHO(IGR)=DIFHO(IGR)+NWT0(IBM,IGR)*STR(IBM)
          ENDDO
        ENDIF
        IF(LHFACT) CALL LCMPUT(KPMAC,'H-FACTOR',NMIL,2,EFACT(1,IGR))
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
      DEALLOCATE(LXS,XS,SS2D,SIGS,XVOLM,EFACT,NWT0)
      IF(NBISO.GT.0) DEALLOCATE(NOMISO)
      DEALLOCATE(NOMREA,IRESM,ITOTM)
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
      ISTATE(:NSTATE)=0
      ISTATE(1)=NGROUP
      ISTATE(2)=NMIL
      ISTATE(3)=NL ! 1+scattering anisotropy
      ISTATE(4)=IFISS
      ISTATE(5)=NED
      ISTATE(6)=ITRANC
      IF(LDIFF) ISTATE(9)=1
      IF(NSURFD.GT.0) ISTATE(12)=3 ! ADF/CPDF information
      CALL LCMPUT(IPMAC,'STATE-VECTOR',NSTATE,1,ISTATE)
      IF(NED.GT.0) CALL LCMPTC(IPMAC,'ADDXSNAME-P0',8,NED,HEDI)
      DEALLOCATE(HEDI)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SIGS0)
      DEALLOCATE(NOMMIL,IJJM,NJJM,IPOS)
      RETURN
      END
