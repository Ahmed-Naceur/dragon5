*DECK MACRDM
      SUBROUTINE MACRDM(IPMACR,IPRINT,IEN   ,NTOTMX,NGROUP,NANISO,
     >                  NBMIXF,NIFISF,NEDF  ,NDELF ,NREACD,NTREA ,
     >                  IMLOC ,NAMREA,NAMEDN,NUMPX ,IXSPRO,XSGEN ,
     >                  XSIGS ,XSSCAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read an old macrolib and transfer information to vectors for a new
* macrolib.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPMACR  pointer to structures.
* IPRINT  print level.
* IEN     macrolib index to process.
* NTOTMX  maximum number of mixtures.
* NGROUP  number of groups.
* NANISO  maximun scattering anisotropy.
* NBMIXF  final number of mixtures.
* NIFISF  final number fissile isotopes.
* NEDF    final number of aditional x-s.
* NDELF   final number of precursor groups.
* NREACD  number of default x-s.
* NTREA   total number of x-s types.
* IMLOC   mixture location.
* NAMREA  names of default x-s.
* NAMEDN  total number of x-s.
* NUMPX   correspondence between old and new 'NUSIGF' arrays.
*
*Parameters: output
* IXSPRO  flag for x-s processing.
* XSGEN   general x-s vector.
* XSIGS   scattering x-s vector.
* XSSCAT  general scattering matrix.
*
*-----------------------------------------------------------------------
*
      USE          GANLIB
      IMPLICIT     NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)  IPMACR
      INTEGER      IPRINT,IEN,NTOTMX,NGROUP,NANISO,NBMIXF,NIFISF,NEDF,
     >             NDELF,NREACD,NTREA,IMLOC(2,NTOTMX),NAMEDN(2,NEDF),
     >             NUMPX(NBMIXF,NIFISF),IXSPRO(NTREA+2*NANISO+1)
      REAL         XSGEN(NBMIXF,NTREA+2),XSIGS(NBMIXF,NANISO),
     >             XSSCAT(NGROUP,NBMIXF,NANISO)
      CHARACTER    NAMREA(NREACD)*12
*----
*  LOCAL VARIABLES
*----
      INTEGER      IOUT
      PARAMETER   (IOUT=6)
      INTEGER      IDEL,ILCMLN,ILCMTY,IMIX,IREA,IREAF,IREAP,IREAA,
     >             IFIS,IED,IANIS,IOMIX,NGF,IGD,IGF,IPOS,IGT,ITC
      CHARACTER    NAMADD*12,CANISO*2,CHID*12,NUSIGD*12
*----
* ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISCAT
      REAL, ALLOCATABLE, DIMENSION(:) :: XSTMP,SCTMP
*----
*  SCRATCH STORAGE ALLOCATION
*   XSTMP   temporary x-s vector
*   SCTMP   temporary scattering matrix
*   ISCAT   scattering pointer
*----
      ALLOCATE(ISCAT(NTOTMX,3))
      ALLOCATE(XSTMP(NTOTMX*(NIFISF+1)),SCTMP(NGROUP*NTOTMX))
*----
*  PRINT HEADER IF REQUIRED
*----
      IF(IPRINT.GE.10) WRITE(IOUT,6000)
*----
*  1) DEFAULT XS
*----
      DO 110 IREA=1,NREACD
        CALL LCMLEN(IPMACR,NAMREA(IREA),ILCMLN,ILCMTY)
        IF(ILCMLN.GT.0) THEN
          IF(IPRINT.GE.10) WRITE(IOUT,6010) NAMREA(IREA)
          IXSPRO(IREA)=1
          CALL LCMGET(IPMACR,NAMREA(IREA),XSTMP)
          DO 100 IMIX=1,NTOTMX
            IF(IMLOC(1,IMIX).EQ.IEN) THEN
              IOMIX=IMLOC(2,IMIX)
              XSGEN(IMIX,IREA)=XSTMP(IOMIX)
            ENDIF
 100      CONTINUE
        ENDIF
 110  CONTINUE
*----
*  2) NUSIGF AND CHI
*----
      IF(NIFISF.GT.0) THEN
        IREAF=NREACD
        IREAP=NREACD+NIFISF
        CALL LCMLEN(IPMACR,'NUSIGF',ILCMLN,ILCMTY)
        IF(ILCMLN.GT.0) THEN
          IF(IPRINT.GE.10) THEN
            WRITE(IOUT,6010) 'NUSIGF      '
            WRITE(IOUT,6010) 'CHI         '
          ENDIF
          CALL LCMGET(IPMACR,'NUSIGF',XSTMP)
          IXSPRO(IREAF+1)=1
          IXSPRO(IREAP+1)=1
          DO 130 IMIX=1,NTOTMX
            IF(IMLOC(1,IMIX).EQ.IEN) THEN
              DO 120 IFIS=1,NIFISF
                IOMIX=NUMPX(IMIX,IFIS)
                IF(IOMIX.NE.0) XSGEN(IMIX,IREAF+IFIS)=XSTMP(IOMIX)
 120          CONTINUE
            ENDIF
 130      CONTINUE
          CALL LCMGET(IPMACR,'CHI',XSTMP)
          DO 150 IMIX=1,NTOTMX
            IF(IMLOC(1,IMIX).EQ.IEN) THEN
              DO 140 IFIS=1,NIFISF
                IOMIX=NUMPX(IMIX,IFIS)
                IF(IOMIX.NE.0) XSGEN(IMIX,IREAP+IFIS)=XSTMP(IOMIX)
 140          CONTINUE
            ENDIF
 150      CONTINUE
        ENDIF
        DO 200 IDEL=1,NDELF
          IREAF=IREAF+2*NIFISF
          IREAP=IREAP+2*NIFISF
          WRITE(NUSIGD,'(A6,I2.2)') 'NUSIGF',IDEL
          WRITE(CHID,'(A3,I2.2)') 'CHI',IDEL
          CALL LCMLEN(IPMACR,NUSIGD,ILCMLN,ILCMTY)
          IF(ILCMLN.GT.0) THEN
            IF(IPRINT.GE.10) THEN
              WRITE(IOUT,6010) NUSIGD
              WRITE(IOUT,6010) CHID
            ENDIF
            CALL LCMGET(IPMACR,NUSIGD,XSTMP)
            IXSPRO(IREAF+1)=1
            IXSPRO(IREAP+1)=1
            DO 170 IMIX=1,NTOTMX
              IF(IMLOC(1,IMIX).EQ.IEN) THEN
                DO 160 IFIS=1,NIFISF
                  IOMIX=NUMPX(IMIX,IFIS)
                  IF(IOMIX.NE.0) XSGEN(IMIX,IREAF+IFIS)=XSTMP(IOMIX)
 160            CONTINUE
              ENDIF
 170        CONTINUE
            CALL LCMGET(IPMACR,CHID,XSTMP)
            DO 190 IMIX=1,NTOTMX
              IF(IMLOC(1,IMIX).EQ.IEN) THEN
                DO 180 IFIS=1,NIFISF
                  IOMIX=NUMPX(IMIX,IFIS)
                  IF(IOMIX.NE.0) XSGEN(IMIX,IREAP+IFIS)=XSTMP(IOMIX)
 180            CONTINUE
              ENDIF
 190        CONTINUE
          ENDIF
 200    CONTINUE
      ENDIF
*----
*  3) ADDITIONAL EDIT XS
*----
      IF(NEDF.GT.0) THEN
        IREAA=NREACD+2*NIFISF*(NDELF+1)
        DO 220 IED=1,NEDF
          WRITE(NAMADD,'(A4,A2)') (NAMEDN(ITC,IED),ITC=1,2)
          CALL LCMLEN(IPMACR,NAMADD,ILCMLN,ILCMTY)
          IF(ILCMLN.GT.0) THEN
            IF(IPRINT.GE.10) WRITE(IOUT,6010) NAMADD
            IXSPRO(IREAA+IED)=1
            CALL LCMGET(IPMACR,NAMADD,XSTMP)
            DO 210 IMIX=1,NTOTMX
              IF(IMLOC(1,IMIX).EQ.IEN) THEN
                IOMIX=IMLOC(2,IMIX)
                XSGEN(IMIX,IREAA+IED)=XSTMP(IOMIX)
              ENDIF
 210        CONTINUE
          ENDIF
 220    CONTINUE
      ENDIF
*----
*  5) SCATTERING XS
*----
      DO 250 IANIS=1,NANISO
        WRITE(CANISO,'(I2.2)') IANIS-1
        CALL LCMLEN(IPMACR,'SCAT'//CANISO,ILCMLN,ILCMTY)
        IF(ILCMLN.GT.0) THEN
          IXSPRO(NTREA+IANIS)=1
          IF(IPRINT.GE.10) WRITE(IOUT,6010) 'SCATTERING'//CANISO
*----
*  4.3) TREAT SCAT
*----
          CALL LCMGET(IPMACR,'IJJS'//CANISO,ISCAT(1,1))
          CALL LCMGET(IPMACR,'NJJS'//CANISO,ISCAT(1,2))
          CALL LCMGET(IPMACR,'IPOS'//CANISO,ISCAT(1,3))
          CALL LCMGET(IPMACR,'SCAT'//CANISO,SCTMP)
          DO 240 IMIX=1,NTOTMX
            IF(IMLOC(1,IMIX).EQ.IEN) THEN
              IOMIX=IMLOC(2,IMIX)
              NGF=ISCAT(IOMIX,2)
              IF(NGF.GT.0) THEN
                IGD=ISCAT(IOMIX,1)
                IGF=IGD-NGF+1
                IPOS=ISCAT(IOMIX,3)
                DO 230 IGT=IGD,IGF,-1
                  XSSCAT(IGT,IMIX,IANIS)=SCTMP(IPOS)
                  IPOS=IPOS+1
 230            CONTINUE
              ENDIF
            ENDIF
 240      CONTINUE
        ENDIF
        CALL LCMLEN(IPMACR,'SIGS'//CANISO,ILCMLN,ILCMTY)
        IF(ILCMLN.GT.0) THEN
          IXSPRO(NTREA+NANISO+IANIS)=1
          CALL LCMGET(IPMACR,'SIGS'//CANISO,XSTMP)
          DO 245 IMIX=1,NTOTMX
            IF(IMLOC(1,IMIX).EQ.IEN) THEN
              IOMIX=IMLOC(2,IMIX)
              XSIGS(IMIX,IANIS)=XSTMP(IOMIX)
            ENDIF
 245      CONTINUE
        ENDIF
 250  CONTINUE
      DEALLOCATE(SCTMP,XSTMP)
      DEALLOCATE(ISCAT)
*----
*  6) STOPPING POWER
*----
      CALL LCMLEN(IPMACR,'ESTOPW',ILCMLN,ILCMTY)
      IF(ILCMLN.GT.0) THEN
        IF(IPRINT.GE.10) WRITE(IOUT,6010) 'ESTOPW'
        ALLOCATE(XSTMP(ILCMLN))
        IXSPRO(NTREA+2*NANISO+1)=1
        CALL LCMGET(IPMACR,'ESTOPW',XSTMP)
        DO 260 IMIX=1,NTOTMX
          IF(IMLOC(1,IMIX).EQ.IEN) THEN
            IOMIX=IMLOC(2,IMIX)
            IF(IOMIX.GT.ILCMLN/2) CALL XABORT('MACRDM: XSTMP OVERFLOW.')
            XSGEN(IMIX,NTREA+1)=XSTMP(IOMIX)
            XSGEN(IMIX,NTREA+2)=XSTMP(ILCMLN/2+IOMIX)
          ENDIF
 260    CONTINUE
        DEALLOCATE(XSTMP)
      ENDIF
      RETURN
*----
*  EDIT FORMATS
*----
 6000 FORMAT(1X,'MACRDM - READING CROSS SECTIONS '/)
 6010 FORMAT(7X,   '   READING RECORD         : ',A12)
      END
