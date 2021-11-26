*DECK SIMLIB
      SUBROUTINE SIMLIB(IMPX,MODE,KPMAP,IPLIB,NTOT,NIS,IFMIX,HFOLLO,
     > RFOLLO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Put/get number densities of particularized isotopes in the microlib
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IMPX    print parameter.
* MODE    transfert mode (=1: get from KPMAP; =2: put to KPMAP).
* KPMAP   HCYCLE subdirectory in the fuelmap.
* IPLIB   pointer to the microlib.
* NTOT    number of fuel bundles.
* NIS     number of particularized isotopes.
* IFMIX   fuel mixture assigned to each fuel bundle.
* HFOLLO  character*8 names of the particularized isotopes.
*
*Parameters: input/output
* RFOLLO  number densities of the particularized isotopes.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KPMAP,IPLIB
      INTEGER IMPX,MODE,NIS,NTOT,IFMIX(NTOT)
      REAL RFOLLO(NTOT,NIS)
      CHARACTER*8 HFOLLO(NIS)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE)
      CHARACTER*12 HCYCL
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IMIX,IVB
      REAL, ALLOCATABLE, DIMENSION(:) :: DENS
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: HUSE
*
      IF(.NOT.C_ASSOCIATED(IPLIB)) THEN
        CALL XABORT('SIMLIB: MICROLIB LCM OBJECT MISSING AT RHS.')
      ENDIF
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      NBMIX=ISTATE(1)
      NBISO=ISTATE(2)
      IF(NTOT.GT.NBMIX) CALL XABORT('SIMLIB: NBMIX OVERFLOW.')
      ALLOCATE(HUSE(NBISO),DENS(NBISO),IMIX(NBISO),IVB(NBMIX))
      CALL LCMGTC(IPLIB,'ISOTOPESUSED',12,NBISO,HUSE)
      CALL LCMGET(IPLIB,'ISOTOPESMIX',IMIX)
      CALL LCMGET(IPLIB,'ISOTOPESDENS',DENS)
      CALL XDISET(IVB,NBMIX,0)
      IBM=0
      DO ITOT=1,NTOT
        IF(IFMIX(ITOT).EQ.0) CYCLE
        IBM=IBM+1
        IVB(IBM)=ITOT
      ENDDO
      CALL LCMGTC(KPMAP,'ALIAS',12,1,HCYCL)
      IF(MODE.EQ.1) THEN
*       recover number densities from KCYCLE directory
        IF(IMPX.GE.0) WRITE(6,'(/34H SIMLIB: recover number densities ,
     >  5Hfrom ,A,11H directory.)') HCYCL
        CALL LCMGET(KPMAP,'FOLLOW',RFOLLO)
        DO ISO=1,NBISO
          IBM=IMIX(ISO)
          ITOT=IVB(IBM)
          IF(ITOT.EQ.0) CALL XABORT('SIMLIB: MISSING FUEL BUNDLE(1).')
          DO JSO=1,NIS
            IF(HUSE(ISO)(:8).EQ.HFOLLO(JSO)) THEN
              DENS(ISO)=RFOLLO(ITOT,JSO)
              GO TO 10
            ENDIF
          ENDDO
   10     CONTINUE
        ENDDO
        CALL LCMPUT(IPLIB,'ISOTOPESDENS',NBISO,2,DENS)
      ELSE IF(MODE.EQ.2) THEN
*       put number densities in KCYCLE directory
        IF(IMPX.GE.0) WRITE(6,'(/33H SIMLIB: put number densities in ,
     >  A,11H directory.)') HCYCL
        CALL XDRSET(RFOLLO,NTOT*NIS,0.0)
        DO ISO=1,NBISO
          IBM=IMIX(ISO)
          ITOT=IVB(IBM)
          IF(ITOT.EQ.0) CALL XABORT('SIMLIB: MISSING FUEL BUNDLE(2).')
          DO JSO=1,NIS
            IF(HUSE(ISO)(:8).EQ.HFOLLO(JSO)) THEN
              RFOLLO(ITOT,JSO)=DENS(ISO)
              GO TO 20
            ENDIF
          ENDDO
   20     CONTINUE
        ENDDO
        CALL LCMPUT(KPMAP,'FOLLOW',NTOT*NIS,2,RFOLLO)
      ELSE
        CALL XABORT('SIMLIB: INVALID VALUE OF MODE.')
      ENDIF
      DEALLOCATE(IVB,IMIX,DENS,HUSE)
      RETURN
      END
