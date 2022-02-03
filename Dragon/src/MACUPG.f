*DECK MACUPG
      SUBROUTINE MACUPG(KENTRY,HGROUP,NENTRY,NIFISF,NDELF,NEDF,NGROUP,
     > NBMIXF,NIFISS,NANISO,NEDMAC,NTOTMX,ITRANC,IPRINT,NAMEN,NUMPX,
     > IMLOC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Cross section processing in GROUP list directory with update.
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
* KENTRY  array of macrolib structures.
* HGROUP  character*12 name of GROUP list directory.
* NENTRY  number of structures.
* NIFISF  maximum number fissile isotopes per mixture.
* NDELF   final number of precursor groups.
* NEDF    final number of aditional x-s.
* NGROUP  number of groups.
* NBMIXF  final number of mixtures.
* NIFISS  number fissile isotopes per mixture.
* NANISO  maximun scattering anisotropy.
* NEDMAC  number of aditional edition x-s.
* NTOTMX  maximum number of mixtures in input macrolibs.
* ITRANC  type of transport correction.
* IPRINT  print level.
* NAMEN   total number of x-s.
* NUMPX   correspondence between old and new 'NUSIGF' arrays.
* IMLOC   mixture location.
*
*-----------------------------------------------------------------------
*
      USE          GANLIB
      IMPLICIT     NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER HGROUP*12
      INTEGER NENTRY,NIFISF,NDELF,NEDF,NGROUP,NBMIXF,NIFISS,NANISO,
     > NEDMAC,NTOTMX,ITRANC,IPRINT,NAMEN(2*NEDMAC),NUMPX(NBMIXF*NIFISS),
     > IMLOC(2,NTOTMX)
      TYPE(C_PTR) KENTRY(NENTRY)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPMACR,JPMACR,KPMACR
      INTEGER IOUT,NREACD,NTREA,IGR,IEN,IMIX,IPRG
      PARAMETER (IOUT=6,NREACD=16)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IXSPRO
      REAL, ALLOCATABLE, DIMENSION(:) :: XSGEN,XSIGS,SCAT
*----
*  DATA
*----
      CHARACTER    NAMREA(NREACD)*12
      SAVE         NAMREA
      DATA         NAMREA
     >  /'NTOT0       ','NTOT1       ','TRANC       ','FIXE        ',
     >   'DIFF        ','DIFFX       ','DIFFY       ','DIFFZ       ',
     >   'NTOT0       ','NTOT1       ','NSPH        ','H-FACTOR    ',
     >   'C-FACTOR    ','OVERV       ','FLUX-INTG   ','FLUX-INTG-P1'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      NTREA=NREACD+2*NIFISF*(1+NDELF)+NEDF
      ALLOCATE(XSGEN(NBMIXF*(NTREA+2)),XSIGS(NBMIXF*NANISO),
     > IXSPRO(NTREA+2*NANISO+1),SCAT(NGROUP*NBMIXF*NANISO*2))
      CALL XDISET(IXSPRO,NTREA+2*NANISO+1,0)
      IPRG=IPRINT
      DO 140 IGR=1,NGROUP
        CALL XDRSET(XSGEN,NBMIXF*(NTREA+2),0.0)
        CALL XDRSET(XSIGS,NBMIXF*NANISO,0.0)
        CALL XDRSET(SCAT,NGROUP*NBMIXF*NANISO,0.0)
        DO 150 IEN=1,NENTRY
          IPMACR=KENTRY(IEN)
          IF(IEN.EQ.1) THEN
*           IPMACR IS OPEN IN CREATION/MODIFICATION MODE
            JPMACR=LCMLID(IPMACR,HGROUP,NGROUP)
            KPMACR=LCMDIL(JPMACR,IGR)
          ELSE
*           IPMACR IS OPEN IN READ-ONLY MODE
            JPMACR=LCMGID(IPMACR,HGROUP)
            KPMACR=LCMGIL(JPMACR,IGR)
          ENDIF
          IF(IPRG.GE.10) WRITE(IOUT,6020) IEN
          DO 151 IMIX=1,NTOTMX
            IF(IMLOC(1,IMIX).EQ.IEN) THEN
              CALL MACRDM(KPMACR,IPRG  ,IEN   ,NTOTMX,NGROUP,NANISO,
     >                    NBMIXF,NIFISF,NEDF  ,NDELF ,NREACD,NTREA ,
     >                    IMLOC ,NAMREA,NAMEN ,NUMPX ,IXSPRO,XSGEN ,
     >                    XSIGS ,SCAT  )
              GO TO 155
            ENDIF
 151      CONTINUE
 155      CONTINUE
 150    CONTINUE
        IPMACR=KENTRY(1)
        JPMACR=LCMLID(IPMACR,HGROUP,NGROUP)
        KPMACR=LCMDIL(JPMACR,IGR)
*----
*  FOR TRANC OFF DO NOT SAVE TRANSPORT CORRECTION
*----
        IF(ITRANC.EQ.0) IXSPRO(2)=0
        CALL MACPRM(KPMACR,IPRG  ,NGROUP,NANISO,NBMIXF,NIFISF,
     >              NEDF  ,NDELF ,NREACD,NTREA ,IGR   ,NAMREA,
     >              NAMEN ,IXSPRO,XSGEN ,XSIGS ,SCAT  )
        IPRG=0
 140  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SCAT,IXSPRO,XSIGS,XSGEN)
      RETURN
 6020 FORMAT(' MACUPG: PROCESSING MACROLIB',I12)
      END
