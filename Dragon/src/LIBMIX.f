*DECK LIBMIX
      SUBROUTINE LIBMIX(IPLIB,NBMIX,NGROUP,NBISO,ISONAM,MIX,DEN,MASK,
     1 MASKL,ITSTMP,TMPDAY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transformation of the isotope ordered microscopic cross sections to
* group ordered macroscopic cross sections (part 1).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* NBMIX   number of material mixtures.
* NGROUP  number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* ISONAM  names of microlib isotopes.
* MIX     mixture number of each isotope (can be zero).
* DEN     density of each isotope.
* MASK    mixture mask (=.true. if a mixture is to be made).
* MASKL   group mask (=.true. if an energy group is to be treated).
* ITSTMP  type of cross section perturbation (=0 perturbation
*         forbidden; =1 perturbation not used even if present;
*         =2 perturbation used if present).
* TMPDAY  time stamp in day/burnup/irradiation.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER NBMIX,NGROUP,NBISO,ISONAM(3,NBISO),MIX(NBISO),ITSTMP
      REAL DEN(NBISO),TMPDAY(3)
      LOGICAL MASK(NBMIX),MASKL(NGROUP)
*----
*  LOCAL VARIABLES
*----
      INTEGER NBLK,NSTATE
      PARAMETER (NBLK=50,NSTATE=40)
      LOGICAL LSAME,LSTOPW
      INTEGER ISTATE(NSTATE),I,IPROB,ITRANC,LENGTH,ITYLCM,MAXNFI,NBESP,
     1 NDEL,NED,NESP,NFISSI,NL,NPART
      CHARACTER TEXT12*12,HPRT1*1
      REAL OLDTIM(3)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JNED
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  RECOVER SOME LIBRARY PARAMETERS.
*----
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      NL=ISTATE(4)
      ITRANC=ISTATE(5)
      IPROB=ISTATE(6)
      NED=ISTATE(13)
      NBESP=ISTATE(16)
      NDEL=ISTATE(19)
      NPART=ISTATE(26)
!      PRINT *,NPART
!      CALL XABORT('AHMED, LIBMIX.f')
      ALLOCATE(JNED(2*NED))
      IF(NED.GT.0) CALL LCMGET(IPLIB,'ADDXSNAME-P0',JNED)
*----
*  LOOK FOR OLD LIBRARY DATA
*----
      CALL LCMLEN(IPLIB,'MACROLIB',LENGTH,ITYLCM)
      IF(LENGTH.EQ.-1) THEN
        CALL LCMSIX(IPLIB,'MACROLIB',1)
        CALL LCMGTC(IPLIB,'SIGNATURE',12,1,TEXT12)
        IF(TEXT12.NE.'L_MACROLIB') THEN
          CALL XABORT('LIBMIX: INVALID SIGNATURE ON THE MACROLIB.')
        ENDIF
        CALL LCMLEN(IPLIB,'TIMESTAMP',LENGTH,ITYLCM)
        IF((LENGTH.GT.0).AND.(LENGTH.LE.3)) THEN
          CALL LCMGET(IPLIB,'TIMESTAMP',OLDTIM)
          IF(ITSTMP.EQ.0) THEN
            TMPDAY(1)=OLDTIM(1)
            TMPDAY(2)=OLDTIM(2)
            TMPDAY(3)=OLDTIM(3)
          ENDIF
        ENDIF
        CALL LCMSIX(IPLIB,' ',2)
      ENDIF
*----
*  SET THE LCM MICROLIB ISOTOPEWISE DIRECTORIES.
*----
      ALLOCATE(IPISO(NBISO))
      CALL LIBIPS(IPLIB,NBISO,IPISO)
*----
*  TRANSPOSE THE MICROSCOPIC CROSS SECTIONS TO ADJOINT ORDERING.
*----
      IF(IPROB.EQ.1) THEN
         CALL LIBADJ (IPLIB,NGROUP,NBISO,NL,NDEL,NBESP,IPISO,NED,JNED)
      ENDIF
*----
*  SET MULTIPLE FISSION SPECTRA INFORMATION.
*----
      IF(NBESP.EQ.0) THEN
        NESP=1
      ELSE
        NESP=NBESP
      ENDIF
*----
*  COMPUTE THE MAXIMUM NUMBER OF FISSIONABLE ISOTOPES IN A MIXTURE.
*----
      DO 20 I=1,NBISO
      IF(MIX(I).GT.NBMIX) CALL XABORT('LIBMIX: NBMIX OVERFLOW.')
   20 CONTINUE
      MAXNFI=MIN(NBISO,200)
      CALL LIBNFI (IPLIB,NGROUP,NBISO,NBMIX,NDEL,NESP,IPISO,MIX,MAXNFI,
     1 NFISSI,LSAME)
*----
*  BUILD THE MACROSCOPIC CROSS SECTIONS.
*----
      CALL LIBDEN (IPLIB,NGROUP,NBISO,NBMIX,NL,NDEL,NESP,ISONAM,IPISO,
     1 MIX,DEN,MASK,MASKL,NED,JNED,ITRANC,NFISSI,NPART,LSAME,ITSTMP,
     2 TMPDAY)
*----
* RECOVER STOPPING POWERS.
*----
      LSTOPW=.FALSE.
      IF(NPART.GT.0) THEN
         CALL LCMGTC(IPLIB,'PARTICLE',1,1,HPRT1)
         LSTOPW=((HPRT1.EQ.'B').OR.(HPRT1.EQ.'C'))
      ENDIF

      !PRINT *,HPRT1
      !PRINT *,LSTOPW
      !CALL XABORT('AHMED,LIBMIX.f')
      !LSTOPW=.TRUE. !AHMED
      
      IF(LSTOPW) THEN
         CALL LIBEST (IPLIB,NGROUP,NBISO,NBMIX,IPISO,MIX,DEN,MASK,MASKL,
     1   NED,JNED,ITSTMP,TMPDAY)
      ENDIF
*----
*  TRANSPOSE THE MICROSCOPIC CROSS SECTIONS BACK TO FORWARD ORDERING.
*----
      IF(IPROB.EQ.1) THEN
         CALL LIBADJ (IPLIB,NGROUP,NBISO,NL,NDEL,NBESP,IPISO,NED,JNED)
      ENDIF
      DEALLOCATE(IPISO,JNED)
      RETURN
      END
