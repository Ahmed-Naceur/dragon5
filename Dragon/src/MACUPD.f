*DECK MACUPD
      SUBROUTINE MACUPD(NENTRY,KENTRY,IPRINT,NTOTMX,NBMIX ,NGROUP,
     >                  NANISO,NIFISS,NEDMAC,ITRANC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Update Dragon macrolib with other Dragon macrolib.
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
* NENTRY  number of structures.
* KENTRY  pointer to structures.
* IPRINT  print level.
* NTOTMX  maximum number of mixtures in input macrolibs.
* NBMIX   number of mixtures on output macrolib.
* NGROUP  number of groups.
* NANISO  maximun scattering anisotropy.
* NIFISS  number fissile isotopes per mixture.
* NEDMAC  number of aditional edition x-s.
* ITRANC  type of transport correction.
*
*-----------------------------------------------------------------------
*
      USE          GANLIB
      IMPLICIT     NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)  KENTRY(NENTRY)
      INTEGER      NENTRY,IPRINT,NTOTMX,NBMIX,NGROUP,NANISO,NIFISS,
     >             NEDMAC,ITRANC
*----
*  LOCAL VARIABLES
*----
      INTEGER      IOUT,NSTATE,MAXPAR
      PARAMETER   (IOUT=6,NSTATE=40,MAXPAR=10)
*----
*  INPUT
*----
      CHARACTER    CARLIR*12
      INTEGER      ITYPLU,INTLIR
      REAL         REALIR
      DOUBLE PRECISION DBLLIR
*----
*  LOCAL PARAMETERS
*----
      TYPE(C_PTR)  IPMACR
      INTEGER      ISTATE(NSTATE),ITEXT4,IMIX,IEN,NUMNEW,NUMOLD,NBMIXF,
     >             NIFISF,NGROF,NEDF,NDELF,NBMIXO,NIFISO,NEDO,NDELO,
     >             NGROO,ILCMLN,ILCMTY,ITC,NPART,I
      CHARACTER    TEXT4*4,HGROUP*12,HPART0*1,HPART(MAXPAR)*1
*----
* ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NAMEN,NUMFN,NUMPX
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IMLOC
      REAL, ALLOCATABLE, DIMENSION(:) :: ENERN
*----
*  SCRATCH STORAGE ALLOCATION
*   IMLOC   mixture location
*----
      ALLOCATE(IMLOC(2,NTOTMX))
*----
*  INITIALIZE IMLOC FOR MIXTURE ALREADY PRESENT ON OUTPUT MACROLIB
*----
      TEXT4='    '
      READ(TEXT4,'(A4)') ITEXT4
      CALL XDISET(IMLOC,2*NTOTMX,0)
      DO 100 IMIX=1,NBMIX
        IMLOC(1,IMIX)=1
        IMLOC(2,IMIX)=IMIX
 100  CONTINUE
      IF(IPRINT.GE.5) WRITE(IOUT,6000)
*----
*  READ INPUT DATA
*  TYPICAL FORMAT:
*  EDIT iprint
*  CTRA { ON | OFF }
*  MIX numnew numold [{ UPDL | OLDL }]
*    WHERE
*    iprint =  PRINT LEVEL
*    numnew =  NEW MIXTURE NUMBER
*    numold =  OLD MIXTURE NUMBER
*    UPDL   =  TAKE numold FORM LIBRARY TO UPDATE
*    OLDL   =  TAKE numold FORM OLD LIBRARY
*----
      ITYPLU = 3
      CARLIR = 'MIX'
 1000 CONTINUE
      IF(ITYPLU.NE.3) CALL XABORT('MACUPD: CHARACTER KEYWORD EXPECTED.')
*----
*  CHECK FOR STOP/RETURN
*----
      IF(CARLIR .EQ. ';') THEN
        GO TO 1005
      ELSE IF(CARLIR(1:3).EQ.'MIX') THEN
*----
*  READ MIX CARD
*----
        CALL REDGET(ITYPLU,NUMNEW,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('MACUPD: NEW MIXTURE NUMBER IS NOT'
     >  //' AN INTEGER.')
        CALL REDGET(ITYPLU,NUMOLD,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('MACUPD: OLD MIXTURE NUMBER IS NOT'
     >  //' AN INTEGER.')
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('MACUPD: STRUCTURE TYPE IS NOT CHA'
     >  //'RACTER.')
*---
*  TEST IF NUMNEW IS VALID
*----
        IF(NUMNEW.GT.NTOTMX) CALL XABORT('MACUPD: NEW MATERIAL NUMBER '
     >  //'IS TOO LARGE.')
*----
*  BY DEFAULT CARLIR IS ASSUMED TO BE OLDL
*  IF CARLIR IS UPDL STORE INFORMATION IN IMLOC
*  KEYWORD OLDL NOT PROCESSED IF PRESENT
*  PROCESS KEYWORD OLDL ACCORDING TO USER'S GUIDE IN THE
*  CASE WHERE IT IS PRESENT
*----
        IF(CARLIR(1:4).EQ.'UPDL') THEN
          IMLOC(1,NUMNEW)=1
          IMLOC(2,NUMNEW)=NUMOLD
        ELSE IF(CARLIR(1:4).EQ.'OLDL') THEN
          IMLOC(1,NUMNEW)=2
          IMLOC(2,NUMNEW)=NUMOLD
        ELSE
          IMLOC(1,NUMNEW)=2
          IMLOC(2,NUMNEW)=NUMOLD
          GO TO 1000
        ENDIF
      ELSE
        CALL XABORT('MACUPD: KEYWORD '//CARLIR//' NOT PERMITTED.')
      ENDIF
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      GO TO 1000
 1005 CONTINUE
*----
*  FIND TOTAL NUMBER OF MIXTURES CREATED
*----
      DO 120 IMIX=NTOTMX,1,-1
        IF(IMLOC(2,IMIX).NE.0) THEN
          NBMIXF=IMIX
          GO TO 125
        ENDIF
 120  CONTINUE
      CALL XABORT('MACUPD: NO MIXTURES FOUND.')
 125  CONTINUE
*----
*  TEST FOR ENERGY
*  FIND ADDITIONAL XS NAME
*  FIND TOTAL NUMBER OF FISSILE ISOTOPES AND THEIR NAME
*----
      ALLOCATE(NAMEN(2*NEDMAC),NUMFN(NBMIXF*NIFISS),
     > NUMPX(NBMIXF*NIFISS),ENERN(2*NGROUP+1))
      CALL XDISET(NAMEN,2*NEDMAC     ,ITEXT4)
*----
*  INITIALIZE VECTOR
*----
      CALL XDISET(NUMFN,NBMIXF*NIFISS,0)
      CALL XDISET(NUMPX,NBMIXF*NIFISS,0)
      NIFISF=0
      NGROF =0
      NEDF  =0
      NDELF =0
      DO 130 IEN=1,NENTRY
        IPMACR=KENTRY(IEN)
        DO 131 IMIX=1,NTOTMX
          IF(IMLOC(1,IMIX).EQ.IEN) THEN
            CALL XDISET(ISTATE,NSTATE,0)
            CALL LCMGET(IPMACR,'STATE-VECTOR',ISTATE)
            NBMIXO=ISTATE(2)
            CALL LCMLEN(IPMACR,'DELTAU',ILCMLN,ILCMTY)
            NGROO=ILCMLN
            NIFISO=ISTATE(4)
            NEDO=ISTATE(5)
            NDELO=ISTATE(7)
            CALL MACNFI(IPMACR,IPRINT,IEN   ,NTOTMX,NGROUP,NIFISS,
     >                  NEDMAC,NBMIXF,NGROF ,NIFISF,NEDF  ,NDELF ,
     >                  NBMIXO,NGROO ,NIFISO,NEDO  ,NDELO ,IMLOC ,
     >                  ENERN ,NAMEN ,NUMFN ,NUMPX )
            GO TO 135
          ENDIF
 131    CONTINUE
 135    CONTINUE
 130  CONTINUE
*----
*  SAVE ENERGY, ADDITIONAL XS NAME AND FISSILE ISOTOPES INFORMATION
*  ON FINAL MACROLIB
*----
      IPMACR=KENTRY(1)
      IF(NGROF.GT.0) THEN
        CALL LCMPUT(IPMACR,'ENERGY',NGROF+1,2,ENERN)
        CALL LCMPUT(IPMACR,'DELTAU',NGROF,2,ENERN(NGROF+2))
        IF(IPRINT.GE.5) THEN
          WRITE(IOUT,6010) 'ENERGY      '
          WRITE(IOUT,6013) (ENERN(ITC),ITC=1,NGROF+1)
          WRITE(IOUT,6010) 'DELTAU      '
          WRITE(IOUT,6013) (ENERN(ITC),ITC=NGROF+2,2*NGROF+1)
        ENDIF
      ENDIF
      IF(NEDF.GT.0) THEN
        CALL LCMPUT(IPMACR,'ADDXSNAME-P0',2*NEDF,3,NAMEN)
        IF(IPRINT.GE.5) THEN
          WRITE(IOUT,6010) 'ADDXSNAME-P0'
          WRITE(IOUT,6011) (NAMEN(ITC),ITC=1,2*NEDO)
        ENDIF
      ENDIF
      IF(NIFISF.GT.0) THEN
        IF(IPRINT.GE.5) THEN
          WRITE(IOUT,6010) 'FISSIONINDEX'
          WRITE(IOUT,6012) (NUMFN(ITC),ITC=1,NBMIXF*NIFISF)
        ENDIF
      ENDIF
      DEALLOCATE(ENERN)
*----
*  CROSS SECTION PROCESSING IN GROUP LIST DIRECTORY WITH UPDATE.
*----
      HGROUP='GROUP'
      CALL MACUPG(KENTRY,HGROUP,NENTRY,NIFISF,NDELF,NEDF,NGROUP,NBMIXF,
     > NIFISS,NANISO,NEDMAC,NTOTMX,ITRANC,IPRINT,NAMEN,NUMPX,IMLOC)
*----
*  RESET NUMFFN TO ONE FOR TERMS WHICH ARE NOT 0 AND SAVE
*----
      IF(NIFISF.GT.0) THEN
        CALL LCMPUT(IPMACR,'FISSIONINDEX',NBMIXF*NIFISF,1,NUMFN)
      ENDIF
*----
*  CCROSS SECTION PROCESSING IN COMPANION GROUP LIST DIRECTORY WITH
*  UPDATE.
*----
      NPART=0
      DO 140 IEN=1,NENTRY
        IPMACR=KENTRY(IEN)
           call lcmlib(IPMACR)
        CALL LCMGET(IPMACR,'STATE-VECTOR',ISTATE)
        IF(ISTATE(17).GT.0) THEN
          NPART=ISTATE(17)+1
          IF(NPART.GT.MAXPAR) CALL XABORT('MAXUPD: MAXPAR OVERFLOW.')
          CALL LCMGTC(IPMACR,'PARTICLE',1,1,HPART0)
          CALL LCMGTC(IPMACR,'PARTICLE-NAM',1,NPART,HPART)
          print *,'HPART=',HPART(:)
          GO TO 150
        ENDIF
  140 CONTINUE
  150 DO I=1,NPART
        IF(HPART(I).EQ.HPART0) CYCLE
        HGROUP='GROUP-'//HPART(I)
        CALL MACUPG(KENTRY,HGROUP,NENTRY,NIFISF,NDELF,NEDF,NGROUP,
     >  NBMIXF,NIFISS,NANISO,NEDMAC,NTOTMX,ITRANC,IPRINT,NAMEN,NUMPX,
     >  IMLOC)
      ENDDO
      DEALLOCATE(NUMPX,NUMFN,NAMEN)
      NBMIX=NBMIXF
      NIFISS=NIFISF
      NEDMAC=NEDF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IMLOC)
      RETURN
*----
*  EDIT FORMATS
*----
 6000 FORMAT(1X,'MACUPD - UPDATING  MACROLIB ')
 6010 FORMAT(7X,   '   PRECESSING RECORD   : ',A12)
 6011 FORMAT(10(2A4,4X))
 6012 FORMAT(10(I8,4X))
 6013 FORMAT(1P,8E15.7)
      END
