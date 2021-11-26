*DECK PSOUR
      SUBROUTINE PSOUR(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Set the fixed source from companion particles in a coupled-particle
* problem.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): creation type(L_SOURCE);
*         HENTRY(2): read-only type(L_MACROLIB) for the source;
*         HENTRY(3): read-only type(L_TRACKING);
*         HENTRY(4): read-only type(L_FLUX) for companion particle 1;
*         HENTRY(5): read-only type(L_FLUX) for companion particle 2;
*         ...
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IOUT=6,NSTATE=40,MAXPAR=10)
      TYPE(C_PTR) IPSOUR,JPSOUR,KPSOUR,IPTRK,IPMAC,JPMAC,KPMAC,JPFLX
      CHARACTER HSIGN*12,TEXT12*12,CMODUL*12,HPART*6,HSMG*131
      INTEGER ISTATE(NSTATE)
      DOUBLE PRECISION DFLOTT
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MATCOD
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FUNKNO,SUNKNO
      CHARACTER(LEN=1), ALLOCATABLE, DIMENSION(:) :: HPARTC
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: HPFLX
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPFLX
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.LT.4) CALL XABORT('PSOUR: FOUR PARAMETERS EXPECTED.')
      IF(MOD(NENTRY,2).EQ.1) CALL XABORT('PSOUR: EVEN NUMBER OF PARAME'
     1 //'TERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('PSOUR: LI'
     1 //'NKED LIST OR XSM FILE EXPECTED AT LHS.')
      IF(JENTRY(1).NE.0) CALL XABORT('PSOUR: ENTRY IN CREATE MODE EXPE'
     1 //'CTED.')
      IPSOUR=KENTRY(1)
      DO I=2,NENTRY
        IF((IENTRY(I).NE.1).AND.(IENTRY(I).NE.2)) CALL XABORT('PSOUR: '
     1  //'LINKED LIST OR XSM FILE EXPECTED AT RHS.')
        IF(JENTRY(I).NE.2) CALL XABORT('PSOUR: ENTRY IN READ-ONLY MODE'
     1  //' EXPECTED.')
      ENDDO
*----
*  RECOVER MACROLIB INFORMATION FOR THE SOURCE
*----
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.EQ.'L_MACROLIB') THEN
        IPMAC=KENTRY(2)
      ELSE IF(HSIGN.EQ.'L_LIBRARY') THEN
        IPMAC=LCMGID(KENTRY(2),'MACROLIB')
      ELSE
        TEXT12=HENTRY(2)
        CALL XABORT('PSOUR: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1  '. L_MACROLIB OR L_LIBRARY EXPECTED.')
      ENDIF
      CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
      NG2=ISTATE(1)
      NMAT=ISTATE(2)
      NANIS=ISTATE(3)-1
      IADJ=ISTATE(13)
      CALL LCMGTC(IPMAC,'PARTICLE',6,1,HPART)
*----
*  RECOVER TRACKING INFORMATION
*----
      IPTRK=KENTRY(3)
      CALL LCMGTC(IPTRK,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_TRACK') THEN
        TEXT12=HENTRY(3)
        CALL XABORT('PSOUR: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1  '. L_TRACK EXPECTED.')
      ENDIF
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NREG=ISTATE(1)
      NUNS=ISTATE(2)
      IF(ISTATE(4).NE.NMAT) CALL XABORT('PSOUR: INVALID NUMBER OF MI'
     1 //'XTURES.')
      ALLOCATE(MATCOD(NREG))
      CALL LCMGET(IPTRK,'MATCOD',MATCOD)
      CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,CMODUL)
*----
*  RECOVER COMPANION FLUX INFORMATION
*----
      NPART0=NENTRY-3
      ALLOCATE(IPFLX(NPART0),HPFLX(NPART0),HPARTC(NPART0))
      DO IPART=1,NPART0
*       Flux for particle IPART
        IPFLX(IPART)=KENTRY(3+IPART)
        HPFLX(IPART)=HENTRY(3+IPART)
        CALL LCMGTC(IPFLX(IPART),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_FLUX') THEN
           TEXT12=HENTRY(3+IPART)
           CALL XABORT('PSOUR: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1     '. L_FLUX EXPECTED.')
        ENDIF
      ENDDO
*----
*  READ THE INPUT DATA
*----
      IMPX=1
      NPART=0
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('PSOUR: CHARACTER DATA EXPECTED.')
      IF(TEXT12.EQ.'EDIT') THEN
*       READ THE PRINT INDEX.
        CALL REDGET(INDIC,IMPX,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('PSOUR: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT12.EQ.'PARTICLE') THEN
*       READ THE PARTICLE TYPE ('N', 'G', 'B', 'C', 'P')
        NPART=NPART+1
        IF(NPART.GT.NPART0) CALL XABORT('PSOUR: NPART0 OVERFLOW.')
        CALL REDGET(INDIC,NITMA,FLOTT,HPARTC(NPART),DFLOTT)
        IF(INDIC.NE.3) CALL XABORT('PSOUR: CHARACTER DATA EXPECTED.')
        TEXT12='GROUP-'//HPARTC(NPART)
        CALL LCMLEN(IPMAC,TEXT12,ILENG,ITYLCM)
        IF(ILENG.EQ.0) THEN
          WRITE(HSMG,'(14HPSOUR: RECORD ,A12,22H IS NOT AVAILABLE IN T,
     1    12HHE MACROLIB.)') TEXT12
          CALL LCMLIB(IPMAC)
          CALL XABORT(HSMG)
        ENDIF
      ELSE IF(TEXT12.EQ.';') THEN
        GO TO 20
      ELSE
        CALL XABORT('PSOUR: '//TEXT12//' IS AN INVALID KEYWORD.')
      ENDIF
      GO TO 10
*----
*  COMPUTE THE FIXED SOURCE
*----
   20 IF(NPART.NE.NPART0) THEN
        WRITE(HSMG,'(30HPSOUR: INVALID NUMBER OF RHS (,I2,10H) NUMBER O,
     1  23HF COMPANION PARTICLES =,I2,1H.)') NENTRY,NPART
        CALL XABORT(HSMG)
      ENDIF
      ALLOCATE(SUNKNO(NUNS,NG2))
      SUNKNO(:NUNS,:NG2)=0.0
      IF(IMPX.GT.0) WRITE(IOUT,100)
      DO IPART=1,NPART
        CALL LCMGET(IPFLX(IPART),'STATE-VECTOR',ISTATE)
        NG1=ISTATE(1)
        NUNF=ISTATE(2)
        ALLOCATE(FUNKNO(NUNF,NG1))
        IF(IMPX.GT.0) WRITE(IOUT,110) IPART,HPARTC(NPART),HPFLX
        TEXT12='GROUP-'//HPARTC(IPART)
        JPMAC=LCMGID(IPMAC,TEXT12)
        IF(IADJ.EQ.0) THEN
          JPFLX=LCMGID(IPFLX(IPART),'FLUX')
        ELSE IF(IADJ.EQ.1) THEN
          JPFLX=LCMGID(IPFLX(IPART),'AFLUX')
        ENDIF
        DO IG=1,NG1
          IF(IADJ.EQ.0) THEN
            CALL LCMGDL(JPFLX,IG,FUNKNO(1,IG))
          ELSE IF(IADJ.EQ.1) THEN
            CALL LCMGDL(JPFLX,NG1-IG+1,FUNKNO(1,IG))
          ENDIF
        ENDDO
        DO IG=1,NG2
          KPMAC=LCMGIL(JPMAC,IG)
          CALL PSOUSN(NUNF,NUNS,IG,IPTRK,KPMAC,NANIS,NREG,NMAT,NG1,NG2,
     1    MATCOD,FUNKNO,SUNKNO)
        ENDDO
        DEALLOCATE(FUNKNO)
      ENDDO
*----
*  SAVE THE FIXED SOURCE ON LCM
*----
      IOF=1
      NDIR=0
      NCST=0
      IF(IADJ.EQ.0) THEN
        NDIR=1
        JPSOUR=LCMLID(IPSOUR,'DSOUR',NDIR)
      ELSE IF(IADJ.EQ.1) THEN
        NCST=1
        JPSOUR=LCMLID(IPSOUR,'ASOUR',NCST)
      ENDIF
      KPSOUR=LCMLIL(JPSOUR,IOF,NG2)
      DO IG=1,NG2
        CALL LCMPDL(KPSOUR,IG,NUNS,2,SUNKNO(1,IG))
      ENDDO
      DEALLOCATE(SUNKNO,MATCOD,HPARTC,HPFLX,IPFLX)
*----
*  SAVE THE SIGNATURE AND STATE VECTOR
*----
      HSIGN='L_SOURCE'
      CALL LCMPTC(IPSOUR,'SIGNATURE',12,1,HSIGN)
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NG2
      ISTATE(2)=NUNS
      ISTATE(3)=NDIR
      ISTATE(4)=NCST
      IF(IMPX.GT.0) WRITE(6,120) (ISTATE(I),I=1,4)
      CALL LCMPUT(IPSOUR,'STATE-VECTOR',NSTATE,1,ISTATE)
      RETURN
*
  100 FORMAT(/44H PSOUR: TRANSITIONS FROM COMPANION PARTICLES/4X,
     1 14HPARTICLE......, 3X,13H PRIMARY_FLUX)
  110 FORMAT(1X,I4,A13,4X,A12)
  120 FORMAT(/8H OPTIONS/8H -------/
     1 7H NG    ,I8,28H   (NUMBER OF ENERGY GROUPS)/
     2 7H NUN   ,I8,40H   (NUMBER OF UNKNOWNS PER ENERGY GROUP)/
     3 7H NDIR  ,I8,35H   (NUMBER OF DIRECT FIXED SOURCES)/
     4 7H NCST  ,I8,36H   (NUMBER OF ADJOINT FIXED SOURCES))
      END
