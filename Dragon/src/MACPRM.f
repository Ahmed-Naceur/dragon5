*DECK MACPRM
      SUBROUTINE MACPRM(IPMACR,IPRINT,NGROUP,NANISO,NBMIXF,NIFISF,
     >                  NEDF  ,NDELF ,NREACD,NTREA ,IGR   ,NAMREA,
     >                  NAMEDN,IXSPRO,XSGEN ,XSIGS ,XSSCAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Save a new macrolib created from old macrolibs.
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
* NGROUP  number of groups.
* NANISO  maximun scattering anisotropy.
* NBMIXF  final number of mixtures.
* NIFISF  final number fissile isotopes.
* NEDF    final number of aditional x-s.
* NDELF   final number of precursor groups.
* NREACD  number of default x-s.
* NTREA   total number of x-s types.
* IGR     current group being processed.
* NAMREA  name of default x-s.
* NAMEDN  name of edit x-s.
* IXSPRO  flag for cross section processing.
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
      INTEGER      IPRINT,NGROUP,NANISO,NBMIXF,NIFISF,NEDF,NDELF,NREACD,
     >             NTREA,IGR,NAMEDN(2,NEDF),IXSPRO(NTREA+2*NANISO+1)
      CHARACTER    NAMREA(NREACD)*12
      REAL         XSGEN(NBMIXF,NTREA+2),XSIGS(NBMIXF,NANISO),
     >             XSSCAT(NGROUP,NBMIXF,NANISO)
*----
*  LOCAL VARIABLES
*----
      INTEGER      IOUT
      PARAMETER   (IOUT=6)
      INTEGER      IDEL,IMIX,IREA,IREAF,IREAP,IREAA,IED,IANIS,IGD,IGF,
     >             JGR,ITC,NELEM,IOF
      CHARACTER    NAMADD*12,CANISO*2,CHID*12,NUSIGD*12
*----
* ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISCAT
      REAL, ALLOCATABLE, DIMENSION(:) :: SCTMP
*----
*  SCRATCH STORAGE ALLOCATION
*   SCTMP   temporary scattering matrix
*   ISCAT   scattering pointer
*----
      ALLOCATE(ISCAT(NBMIXF,3))
      ALLOCATE(SCTMP(NGROUP*NBMIXF))
*----
*  PRINT HEADER
*----
      IF(IPRINT.GE.10) WRITE(IOUT,6000)
*----
*  1) DEFAULT XS
*----
      DO 100 IREA=1,NREACD
        IF(IXSPRO(IREA).EQ.1) THEN
          IF(IPRINT.GE.10) WRITE(IOUT,6010) NAMREA(IREA)
          CALL LCMPUT(IPMACR,NAMREA(IREA),NBMIXF,2,XSGEN(1,IREA))
        ENDIF
 100  CONTINUE
*----
*  2) CHI AND NUSIGF
*----
      IF(NIFISF.GT.0) THEN
        IREAF=NREACD+1
        IREAP=NREACD+NIFISF+1
        IF(IXSPRO(IREAF).EQ.1) THEN
          IF(IPRINT.GE.10) WRITE(IOUT,6010) 'NUSIGF      '
          CALL LCMPUT(IPMACR,'NUSIGF',NBMIXF*NIFISF,2,XSGEN(1,IREAF))
        ENDIF
        IF(IXSPRO(IREAP).EQ.1) THEN
          IF(IPRINT.GE.10) WRITE(IOUT,6010) 'CHI         '
          CALL LCMPUT(IPMACR,'CHI',NBMIXF*NIFISF,2,XSGEN(1,IREAP))
        ENDIF
        DO 110 IDEL=1,NDELF
          WRITE(NUSIGD,'(A6,I2.2)') 'NUSIGF',IDEL
          IF(IPRINT.GE.10) WRITE(IOUT,6010) NUSIGD
          IOF=IREAF+2*IDEL*NIFISF
          CALL LCMPUT(IPMACR,NUSIGD,NBMIXF*NIFISF,2,XSGEN(1,IOF))
          WRITE(CHID,'(A3,I2.2)') 'CHI',IDEL
          IF(IPRINT.GE.10) WRITE(IOUT,6010) CHID
          IOF=IREAP+2*IDEL*NIFISF
          CALL LCMPUT(IPMACR,CHID,NBMIXF*NIFISF,2,XSGEN(1,IOF))
 110    CONTINUE
      ENDIF
*----
*  3) ADDITIONAL EDIT XS
*----
      IF(NEDF.GT.0) THEN
        IREAA=NREACD+2*NIFISF*(NDELF+1)
        DO 120 IED=1,NEDF
          IF(IXSPRO(IREAA+IED).EQ.1) THEN
            WRITE(NAMADD,'(A4,A2)') (NAMEDN(ITC,IED),ITC=1,2)
            IF(IPRINT.GE.10) WRITE(IOUT,6010) NAMADD
            CALL LCMPUT(IPMACR,NAMADD,NBMIXF,2,XSGEN(1,IREAA+IED))
          ENDIF
 120    CONTINUE
      ENDIF
*----
*  4) SCATTERING XS
*----
      DO 130 IANIS=1,NANISO
        IF(IXSPRO(NTREA+IANIS).EQ.1) THEN
          WRITE(CANISO,'(I2.2)') IANIS-1
          IF(IPRINT.GE.10) WRITE(IOUT,6010) 'SCATTERING'//CANISO
          CALL XDRSET(XSGEN,NBMIXF,0.0)
*----
*  4.3) TREAT SCAT
*----
          NELEM=0
          DO 140 IMIX=1,NBMIXF
            IGD=IGR
            IGF=IGR
            XSGEN(IMIX,1)=XSSCAT(IGR,IMIX,IANIS)
            DO 141 JGR=1,NGROUP
              IF(XSSCAT(JGR,IMIX,IANIS).NE.0.0) THEN
                IGD=MAX(IGD,JGR)
                IGF=MIN(IGF,JGR)
              ENDIF
 141        CONTINUE
            ISCAT(IMIX,1)=IGD
            ISCAT(IMIX,2)=IGD-IGF+1
            ISCAT(IMIX,3)=NELEM+1
            DO 142 JGR=IGD,IGF,-1
              NELEM=NELEM+1
              SCTMP(NELEM)=XSSCAT(JGR,IMIX,IANIS)
 142        CONTINUE
 140      CONTINUE
          CALL LCMPUT(IPMACR,'SIGW'//CANISO,NBMIXF,2,XSGEN(1,1))
          CALL LCMPUT(IPMACR,'IJJS'//CANISO,NBMIXF,1,ISCAT(1,1))
          CALL LCMPUT(IPMACR,'NJJS'//CANISO,NBMIXF,1,ISCAT(1,2))
          CALL LCMPUT(IPMACR,'IPOS'//CANISO,NBMIXF,1,ISCAT(1,3))
          CALL LCMPUT(IPMACR,'SCAT'//CANISO,NELEM,2,SCTMP)
        ENDIF
        IF(IXSPRO(NTREA+NANISO+IANIS).EQ.1) THEN
          CALL LCMPUT(IPMACR,'SIGS'//CANISO,NBMIXF,2,XSIGS(1,IANIS))
        ENDIF
 130  CONTINUE
      DEALLOCATE(SCTMP)
      DEALLOCATE(ISCAT)
*----
*  6) STOPPING POWER
*----
      IF(IXSPRO(NTREA+2*NANISO+1).EQ.1) THEN
        IF(IPRINT.GE.10) WRITE(IOUT,6010) 'ESTOPW'
        CALL LCMPUT(IPMACR,'ESTOPW',2*NBMIXF,2,XSGEN(1,NTREA+1))
      ENDIF
      RETURN
*----
*  EDIT FORMATS
*----
 6000 FORMAT(1X,'MACPRM - SAVING CROSS SECTIONS '/)
 6010 FORMAT(7X,   '   SAVING RECORD         : ',A12)
      END
