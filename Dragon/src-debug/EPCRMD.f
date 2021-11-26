*DECK EPCRMD
      SUBROUTINE EPCRMD(NENTRY,KENTRY,IPRINT,NOPT,IOPT,CARRET)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for Error Propagation Module with option "RNDMicXS".
*
*Copyright:
* Copyright (C) 2009 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* NENTRY  number of data structures transfered to this module.
* KENTRY  data structure pointer.
* IPRINT  print level.
* NOPT    length of array IOPT containing the procession option flags.
* IOPT    processing option with:
*         IOPT(1)  type of processing (3 for current option);
*         IOPT(2)  entry number for L_EPC structure;
*         IOPT(3)  entry number for read only SEQ_ASCII or BINARY file;
*         IOPT(4)  entry number for update microlib.
* CARRET  last input option read.
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NENTRY
      TYPE(C_PTR)      KENTRY(NENTRY)
      INTEGER          IPRINT,NOPT,IOPT(NOPT)
      CHARACTER*12     CARRET
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='EPCRMD')
      INTEGER          ILCMUP,ILCMDN
      PARAMETER       (ILCMUP=1,ILCMDN=2)
      INTEGER          NSTATE,NTC
      PARAMETER       (NSTATE=40,NTC=3)
      INTEGER          NDPROC
      PARAMETER       (NDPROC=20)
*----
*  Input and output parameters
*----
      INTEGER          ITYPLU,INTLIR
      CHARACTER        CARLIR*12
      REAL             REALIR
      DOUBLE PRECISION DBLLIR
*----
*  Local variables
*----
      TYPE(C_PTR)      IPEPC,IPMIC
      INTEGER          IPCOV,IFMT,IEN,IPRTL,IXS
      INTEGER          ISTATE(NSTATE),ISTATM(NSTATE)
      CHARACTER*6      NAMDXS(NDPROC+1)
      CHARACTER*8      XSN
      INTEGER          NGR,NIS,NXS,NCV,ISOREC
      INTEGER          NMIXT,NBISO,NGROUP,NIFISS,NFI,NISL
      INTEGER          ITOTL,ISCAT
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NAMISO,NISOU,ISOMIX,IDVF,
     > IDMF,IDXS
*----
*  Read instruction for RNDMicXS
*----
      IPRTL=IPRINT
      IF(IPRTL .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      IF(CARRET .EQ. ';') GO TO 105
 100  CONTINUE
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      CONTINUE
      IF(ITYPLU .EQ. 10) GO TO 105
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': Read error -- Character variable expected')
      CARRET=CARLIR
      IF(CARLIR .EQ. ';') THEN
        GO TO 105
      ELSE
        CALL XABORT(NAMSBR//': Read error -- keyword'//CARLIR//
     >'for processing option is invalid.')
      ENDIF
      GO TO 100
 105  CONTINUE
*----
*  Output structure
*----
      IEN=ABS(IOPT(3))
      IPEPC=KENTRY(ABS(IOPT(2)))
      IPCOV=FILUNIT(KENTRY(IEN))
      IPMIC=KENTRY(ABS(IOPT(4)))
      IF(IOPT(3) .NE. 0) THEN
*----
*  Transfer Variance and covariance data on EPC data structure
*----
        IFMT=IOPT(3)/IEN
        IF(IFMT .GT. 0) THEN
          READ(IPCOV,'(3I8)') NGR,NIS,NXS
        ELSE
          READ(IPCOV) NGR,NIS,NXS
        ENDIF
        IF(IPRTL .GE. 10) THEN
          WRITE(IOUT,6010) NGR,NIS,NXS
        ENDIF
        NCV=(NGR*(NGR+1))/2
        CALL EPCRMV(IPEPC,IPCOV,IPRINT,IFMT,NGR,NIS,NXS,NCV)
        CALL XDISET(ISTATE,NSTATE,0)
        ISTATE(1)=IOPT(1)
        ISTATE(2)=NGR
        ISTATE(3)=NIS
        ISTATE(4)=NXS
        CALL LCMPUT(IPEPC,'STATE-VECTOR',NSTATE,1,ISTATE)
      ENDIF
      IF(IOPT(4) .NE. 0) THEN
*----
*  Get info for EPC data structure
*----
        CALL XDISET(ISTATE,NSTATE,0)
        CALL LCMGET(IPEPC,'STATE-VECTOR',ISTATE)
        IF(ISTATE(1) .NE. 3) CALL XABORT(NAMSBR//
     >': Invalid format for EPC data structure')
        NGR=ISTATE(2)
        NIS=ISTATE(3)
        NXS=ISTATE(4)
        ISOREC=ISTATE(5)
        NFI=ISTATE(6)
        NISL=ISTATE(7)
        IF(IPRTL .GE. 10) THEN
          WRITE(IOUT,6011) NGR,NIS,NXS
        ENDIF
        NCV=(NGR*(NGR+1))/2
*----
*  Get info for microlib data structure
*----
        CALL XDISET(ISTATM,NSTATE,0)
        CALL LCMGET(IPMIC,'STATE-VECTOR',ISTATM)
        NMIXT=ISTATM(1)
        NBISO=ISTATM(2)
        NGROUP=ISTATM(3)
        IF(NGROUP .NE. NGR) CALL XABORT(NAMSBR//
     >': Number of groups in MICROLIB and EPC incoherent.')
        CALL LCMSIX(IPMIC,'MACROLIB    ',ILCMUP)
        CALL XDISET(ISTATM,NSTATE,0)
        CALL LCMGET(IPMIC,'STATE-VECTOR',ISTATM)
        NIFISS=ISTATM(4)
        CALL LCMSIX(IPMIC,'MACROLIB    ',ILCMDN)
*----
*  Allocate memory
*----
        ALLOCATE(NAMISO(3*NIS),NISOU(3*NBISO),ISOMIX(NBISO),IDVF(2*NIS),
     >  IDMF(2*NBISO))
        CALL LCMGET(IPMIC,'ISOTOPESUSED',NISOU)
        CALL LCMGET(IPMIC,'ISOTOPESMIX ',ISOMIX)
        CALL LCMSIX(IPEPC,'XSVariances ',ILCMUP)
        CALL LCMGET(IPEPC,'NAMEISO     ',NAMISO)
*----
*  Variance XS
*----
        ALLOCATE(IDXS(2*NXS))
        CALL LCMGET(IPEPC,'NAMEXS      ',IDXS)
        ISCAT=0
        ITOTL=0
        DO IXS=1,NXS
          WRITE(XSN,'(2A4)') IDXS(2*IXS-1),IDXS(2*IXS)
          NAMDXS(IXS)=XSN(1:6)
          IF(NAMDXS(IXS).EQ.'TOTAL ') THEN
            ITOTL=IXS
          ELSE IF(NAMDXS(IXS).EQ.'SIGS00') THEN
            ISCAT=IXS
          ENDIF
        ENDDO
        DEALLOCATE(IDXS)
        IF(ISOREC .EQ. 0) THEN
          CALL EPCRMI(IPMIC,IPRINT,NIS,NBISO,NMIXT,NIFISS,
     >                NAMISO,NISOU,IDVF,IDMF)
          CALL LCMPUT(IPEPC,'INDEXISOV',2*NIS,1,IDVF)
          CALL LCMPUT(IPEPC,'INDEXISOM',2*NBISO,1,IDMF)
          ISOREC=1
          NFI=NIFISS
          NISL=NBISO
          ISTATE(5)=ISOREC
          ISTATE(6)=NFI
          ISTATE(7)=NISL
          CALL LCMPUT(IPEPC,'STATE-VECTOR',NSTATE,1,ISTATE)
        ELSE
          IF(NFI .NE. NIFISS) CALL XABORT(NAMSBR//
     >': Number of fissile isotopes in MICROLIB and EPC incoherent.')
          IF(NISL .NE. NBISO) CALL XABORT(NAMSBR//
     >': Number of isotopes in MICROLIB and EPC incoherent.')
          CALL LCMGET(IPEPC,'INDEXISOV',IDVF)
          CALL LCMGET(IPEPC,'INDEXISOM',IDMF)
        ENDIF
*----
*  Update macrolib from random error distribution
*----
        CALL EPCRMU(IPEPC,IPMIC,IPRINT,NGR,NIS,NXS,NCV,
     >              NBISO,NMIXT,NIFISS,ITOTL,ISCAT,NAMDXS,
     >              NAMISO,NISOU,ISOMIX,IDVF,IDMF)
        DEALLOCATE(IDMF,IDVF,ISOMIX,NISOU,NAMISO)
        CALL LCMSIX(IPEPC,'XSVariances ',ILCMDN)
      ENDIF
      IF(IPRTL .GE. 10) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT('  Creating EPC data structure with:'/
     >       '  NGROUP  = ',I8/
     >       '  NISO    = ',I8/
     >       '  NXS     = ',I8)
 6011 FORMAT('  Updating MICROLIB from EPC with:'/
     >       '  NGROUP  = ',I8/
     >       '  NISO    = ',I8/
     >       '  NXS     = ',I8)
      END
