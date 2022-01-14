*DECK MACNFI
      SUBROUTINE MACNFI(IPMACR,IPRINT,IEN   ,NTOTMX,NGROUP,NIFISS,
     >                  NEDMAC,NBMIXF,NGROF ,NIFISF,NEDF  ,NDELF ,
     >                  NBMIXO,NIFISO,NEDO  ,NDELO ,IMLOC ,ENERGN,
     >                  NAMEDN,NUMFN ,NUMPX )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Update list of fissile isotopes from fissile isotopes on a specific
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
* NIFISS  maximum number fissile isotopes per mixture.
* NEDMAC  number of aditional edition x-s.
* NBMIXO  number of mixtures in IPMACR.
* NIFISO  number of fissile isotopes in IPMACR.
* NEDO    number of aditional x-s in IPMACR.
* NDELO   number of precursor groups in IPMACR.
* IMLOC   mixture location.
* NBMIXF  final number of mixtures.
*
*Parameters: input/output
* NGROF   number of groups tested.
* NIFISF  final number fissile isotopes.
* NEDF    final number of aditional x-s.
* NDELF   final number of precursor groups.
* ENERGN  final energy/lethargy vector.
* NAMEDN  final edit names.
* NUMFN   final 'FISSIONINDEX' record.
* NUMPX   correspondence between old and new 'NUSIGF' arrays.
*
*-----------------------------------------------------------------------
*
      USE          GANLIB
      IMPLICIT     NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)  IPMACR
      INTEGER      IPRINT,IEN ,NTOTMX,NGROUP,NIFISS,NEDMAC,
     >             NBMIXF,NGROF,NIFISF,NEDF,NDELF,NBMIXO,
     >             NIFISO,NEDO,NDELO,IMLOC(2,NTOTMX),
     >             NAMEDN(2,NEDMAC),NUMFN(NBMIXF,NIFISS),
     >             NUMPX(NBMIXF,NIFISS)
      REAL         ENERGN(2*NGROUP+1)
*----
*  LOCAL VARIABLES
*----
      INTEGER      IOUT
      PARAMETER   (IOUT=6)
      INTEGER      IGR,ILO,ILN,JLN,IMXN,IMAC,IMIX,ITC,ISOT,ILCMLN,
     >             ILCMTY,NGROO
*----
* ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NAMEDO,NUMFO
      REAL, ALLOCATABLE, DIMENSION(:) :: ENERGO
*----
*  SCRATCH STORAGE ALLOCATION
*   ENERGO  energy/lethargy vector in IPMACR
*   NAMEDO  edit names in IPMACR
*   NUMFO   'FISSIONINDEX' record in IPMACR
*----
      ALLOCATE(NAMEDO(2,NEDO),NUMFO(NBMIXO,NIFISO))
*----
*  PRINT HEADER IF REQUIRED
*----
      IF(IPRINT.GE.10) WRITE(IOUT,6000) IEN
*----
*  TEST FOR ENERGY
*----
      NGROO=0
      CALL LCMLEN(IPMACR,'ENERGY',ILCMLN,ILCMTY)
      IF(ILCMLN.GT.0) THEN
        NGROO=ILCMLN-1
        ALLOCATE(ENERGO(2*NGROO+1))
        IF(NGROF.GT.0) THEN
          CALL LCMGET(IPMACR,'ENERGY',ENERGO(1))
          DO IGR=1,NGROO
            ENERGO(NGROO+1+IGR)=LOG(ENERGO(IGR)/ENERGO(IGR+1))
          ENDDO
          DO IGR=1,2*NGROO+1
            IF(ENERGN(IGR).NE.ENERGO(IGR)) THEN
              WRITE(IOUT,9000) IEN
              GO TO 110
            ENDIF
          ENDDO
        ELSE
          CALL LCMGET(IPMACR,'ENERGY',ENERGN(1))
          DO IGR=1,NGROO
            ENERGN(NGROO+1+IGR)=LOG(ENERGN(IGR)/ENERGN(IGR+1))
          ENDDO
          NGROF=NGROO
        ENDIF
      ENDIF
*----
*  TEST FOR ADDITIONAL EDIT XS
*----
 110  IF(NEDO.GT.0) THEN
        CALL LCMGET(IPMACR,'ADDXSNAME-P0',NAMEDO)
        IF(IPRINT.GE.10) THEN
          WRITE(IOUT,6010) 'ADDXSNAME-P0'
          WRITE(IOUT,6011) ((NAMEDO(ITC,ILO),ITC=1,2),ILO=1,NEDO)
        ENDIF
        DO 140 ILO=1,NEDO
          DO 120 ILN=1,NEDF
            IF( NAMEDO(1,ILO) .EQ. NAMEDN(1,ILN) .AND.
     >          NAMEDO(1,ILO) .EQ. NAMEDN(1,ILN)      ) GO TO 130
 120      CONTINUE
          NEDF=NEDF+1
          NAMEDN(1,ILN)=NAMEDO(1,ILO)
          NAMEDN(2,ILN)=NAMEDO(2,ILO)
 130      CONTINUE
 140    CONTINUE
      ENDIF
*----
*  TEST FOR PRECURSOR GROUPS
*----
      IF(NDELO.GT.0) THEN
        IF(NDELF.EQ.0) THEN
          NDELF=NDELO
        ELSE IF(NDELF.NE.NDELO) THEN
          CALL XABORT('MACNFI: INVALID NUMBER OF PRECURSOR GROUPS.')
        ENDIF
      ENDIF
*----
*  TEST FOR FISSILE ISOTOPES NAMES
*  STORE IN NUMFN THE LOCATION OF CROSS SECTION IN OLD NUSIGF AND CHI
*----
      IF(NIFISO.GT.0) THEN
        CALL LCMLEN(IPMACR,'FISSIONINDEX',ILCMLN,ILCMTY)
        IF(ILCMLN.EQ.0) THEN
          IF(NIFISO.EQ.1) THEN
*          IF(NIFISF.GT.1) CALL XABORT('MACNFI: MISSING FISSIONINDEX RE'
*     >    //'CORD.')
          DO 145 IMXN=1,NBMIXF ! loop over new mixture indices
            IMAC=IMLOC(1,IMXN) ! old macrolib index
            IMIX=IMLOC(2,IMXN) ! old mixture index
            IF(IMAC.EQ.IEN) THEN
              NIFISF=1
              NUMFN(IMXN,1)=1
              NUMPX(IMXN,1)=IMIX
            ENDIF
 145        CONTINUE
            GO TO 190
          ENDIF
          CALL XDISET(NUMFO,NBMIXO*NIFISO,-1)
        ELSE
          IF(ILCMLN.GT.NBMIXO*NIFISO)
     >    CALL XABORT('MACNFI: FISSIONINDEX OVERFLOW,')
          CALL LCMGET(IPMACR,'FISSIONINDEX',NUMFO)
        ENDIF
        IF(IPRINT.GE.10) THEN
          WRITE(IOUT,6010) 'FISSIONINDEX'
          WRITE(IOUT,6012) ((NUMFO(ITC,ILO),ITC=1,NBMIXO),ILO=1,NIFISO)
        ENDIF
        DO 180 IMXN=1,NBMIXF ! loop over new mixture indices
          IMAC=IMLOC(1,IMXN) ! old macrolib index
          IMIX=IMLOC(2,IMXN) ! old mixture index
          IF(IMAC.EQ.IEN) THEN
            DO 170 ILO=1,NIFISO ! loop over old fissile isotopes
              ISOT=NUMFO(IMIX,ILO) ! a reference to the old microlib
              DO 150 JLN=1,NIFISF ! loop over new fissile isotopes
                IF(NUMFN(IMXN,JLN).EQ.ISOT) GO TO 170
 150          CONTINUE
              DO 160 JLN=1,NIFISF
                IF(NUMFN(IMXN,JLN).EQ.0) THEN
                  NUMFN(IMXN,JLN)=ISOT
                  NUMPX(IMXN,JLN)=(ILO-1)*NBMIXO+IMIX
                  GO TO 170
                ENDIF
 160          CONTINUE
              NIFISF=NIFISF+1
              IF(NIFISF.GT.NIFISS) CALL XABORT('MACNFI: NUMFN OVERFLOW')
              NUMFN(IMXN,NIFISF)=ISOT
              NUMPX(IMXN,NIFISF)=(ILO-1)*NBMIXO+IMIX
 170        CONTINUE
          ENDIF
 180    CONTINUE
 190    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      IF(NGROO.GT.0) DEALLOCATE(ENERGO)
      DEALLOCATE(NUMFO,NAMEDO)
      RETURN
*----
*  EDIT FORMATS
*----
 6000 FORMAT(1X,'MACNFI - PROCESSING MACROLIB : ',I12)
 6010 FORMAT(7X,   '   PRECESSING RECORD   : ',A12)
 6011 FORMAT(10(2A4,4X))
 6012 FORMAT(10(I8,4X))
*----
*  WARNING FORMATS
*----
 9000 FORMAT(' **** WARNING IN MACNFI FOR MACROLIB : ',I12/
     >       '      ENERGY GROUP STRUCTURE NOT COMPATIBLE'/
     >       ' **** CORRECTION:  USE LAST ENERGY STRUCTURE')
      END
