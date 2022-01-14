*DECK LIBEST
      SUBROUTINE LIBEST (IPLIB,NGROUP,NBISO,NBMIX,IPISO,MIX,DEN,
     1 MASK,MASKL,NED,NAMEAD,ITSTMP,TMPDAY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the stopping powers from microscopic cross-section library
* and generate record 'ESTOPW' in the group ordered macroscopic
* cross-section library.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and A. Naceur
*
*Parameters: input
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* NGROUP  number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NBMIX   number of mixtures present in the calculation domain.
* IPISO   pointer array towards microlib isotopes.
* MIX     mixture number of each isotope (can be zero).
* DEN     density of each isotope.
* MASK    mixture mask (=.true. if a mixture is to be made).
* MASKL   group mask (=.true. if an energy group is to be treated).
* NED     number of extra edit vectors.
* NAMEAD  names of these extra edits.
* ITSTMP  type of cross section perturbation (=0: perturbation
*         forbidden; =1: perturbation not used even if present;
*         =2: perturbation used if present).
* TMPDAY  time stamp in day/burnup/irradiation.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
      INTEGER NGROUP,NBISO,NBMIX,MIX(NBISO),NED,NAMEAD(2,NED),ITSTMP
      REAL DEN(NBISO),TMPDAY(3)
      REAL DENMIXR, DENMIX(NBMIX)!AHMED
      LOGICAL MASK(NBMIX),MASKL(NGROUP)
      LOGICAL STERNCORR
*----
*  LOCAL VARIABLES
*----
      INTEGER NSTATE,IOUT
      PARAMETER (NSTATE=40,IOUT=6)
      CHARACTER CV*12,HSMG*131,TEXT12*12,HPRT1*1,NORD(3)*4
      LOGICAL EXIST
      INTEGER IDATA(NSTATE),I0,LLL,IBM,IED,NXSPER,ISOT,ILONG,
     1 ITYLCM,IXSPER,NGROUPS,ICV
      REAL TMPPER(2,3),TIMFCT,DENISO,XTF
      TYPE(C_PTR) JPLIB,KPLIB
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KGAS
      REAL, ALLOCATABLE, DIMENSION(:) :: GA1,ENER
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ESTOP,GAF
      REAL, ALLOCATABLE, DIMENSION(:,:) :: DENMAT
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPGRP
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: ISONRF
*----
*  DATA STATEMENTS
*----
      DATA NORD/'    ',' LIN',' QUA'/
*----
*  SCRATCH STORAGE ALLOCATION
*   IPGRP   LCM pointers of the macrolib groupwise directories.
*----
      ALLOCATE(KGAS(NBMIX))
      ALLOCATE(ESTOP(NBMIX,NGROUP+1,2),GAF(NBMIX,NGROUP+1,2))
      ALLOCATE(DENMAT(NBMIX,NGROUP+1)) !AHMED
      ALLOCATE(IPGRP(NGROUP))
*----
*  SET CROSS SECTION PERTURBATION INFORMATION.
*----
      NXSPER=1
      TIMFCT=0.0
      CALL LCMLEN(IPLIB,'TIMESPER',ILONG,ITYLCM)
      IF((ILONG.GE.2).AND.(ILONG.LE.6)) THEN
        IF(ITSTMP.EQ.0) THEN
          CALL XABORT('LIBDEN: XS PERTURBATION FORBIDDEN.')
        ELSE IF(ITSTMP.EQ.2) THEN
          CALL LCMGET(IPLIB,'TIMESPER',TMPPER)
          TIMFCT=TMPDAY(1)-TMPPER(1,1)
          XTF=TIMFCT/TMPPER(2,1)
          IF(XTF.NE.0.0) NXSPER=2
          IF(XTF.LT.0.0) THEN
            WRITE(IOUT,6000) TMPPER(1,1),TMPDAY(1)
          ELSE IF(XTF.GT.1.0) THEN
            WRITE(IOUT,6001) TMPPER(1,1)+TMPPER(2,1),TMPDAY(1)
          ENDIF
        ENDIF
      ENDIF
*----
*  VALIDATE MACROLIB SIGNATURE AND STATE-VECTOR.
*----
      CALL LCMSIX(IPLIB,'MACROLIB',1)
      CALL LCMGTC(IPLIB,'SIGNATURE',12,1,TEXT12)
      IF(TEXT12.NE.'L_MACROLIB') THEN
         CALL XABORT('LIBEST: INVALID SIGNATURE ON THE MACROLIB.')
      ENDIF
      CALL LCMGTC(IPLIB,'PARTICLE',1,1,HPRT1)
      IF((HPRT1.NE.'B').AND.(HPRT1.NE.'C')) THEN
         CALL XABORT('LIBEST: INVALID PARTICLE TYPE. B OR C EXPECTED.')
      ENDIF
      CALL LCMGET(IPLIB,'STATE-VECTOR',IDATA)
      IF(IDATA(1).NE.NGROUP) THEN
         WRITE(HSMG,'(37HLIBEST: EXISTING MACROLIB HAS NGROUP=,I4,
     1   25H NEW MACROLIB HAS NGROUP=,I4,1H.)') IDATA(1),NGROUP
         CALL XABORT(HSMG)
      ELSE IF(IDATA(2).NE.NBMIX) THEN
         WRITE(HSMG,'(36HLIBEST: EXISTING MACROLIB HAS NBMIX=,I4,
     1   24H NEW MACROLIB HAS NBMIX=,I4,1H.)') IDATA(2),NBMIX
         CALL XABORT(HSMG)
      ENDIF
      CALL LCMSIX(IPLIB,' ',2)
*----
*  SELECT NUMBER OF GROUPS TO PROCESS
*----
      NGROUPS=0
      DO LLL=1,NGROUP
         IF(MASKL(LLL)) NGROUPS=NGROUPS+1
      ENDDO
      IF(NGROUPS.EQ.0) GO TO 50
*----
*  SET THE LCM MACROLIB GROUPWISE AND MICROLIB ISOTOPEWISE DIRECTORIES
*----
      CALL LCMSIX(IPLIB,'MACROLIB',1)
      JPLIB=LCMLID(IPLIB,'GROUP',NGROUP)
      DO LLL=1,NGROUP
         IPGRP(LLL)=LCMDIL(JPLIB,LLL)
      ENDDO
      CALL LCMSIX(IPLIB,' ',2)
*----
*  PROCESS THE STOPPING POWERS.
*----
      ESTOP(:NBMIX,:NGROUP+1,:2)=0.0
      EXIST=.FALSE.
      ALLOCATE(GA1(NGROUP+1))
      DO 40 IED=1,NED
      WRITE(CV,'(2A4)') (NAMEAD(I0,IED),I0=1,2)
      IF((CV(:3).EQ.'BST').OR.(CV(:3).EQ.'CST')) THEN
         ICV=0
         IF(CV(2:4).EQ.'STC') THEN
            ICV=1
         ELSE IF(CV(2:4).EQ.'STR') THEN
            ICV=2
         ELSE
            CALL XABORT('LIBEST: BSTC, BSTR, CSTC OR CSTR EXPECTED.')
         ENDIF
         DO 30 IBM=1,NBMIX
         IF(MASK(IBM)) THEN
            DO 20 ISOT=1,NBISO
            IF((MIX(ISOT).NE.IBM).OR.(DEN(ISOT).EQ.0.0)) GO TO 20
            JPLIB=IPISO(ISOT)
            IF(.NOT.C_ASSOCIATED(JPLIB)) GO TO 20
*-
            DENISO=DEN(ISOT)
            DO IXSPER=1,NXSPER
               CALL LCMLEN(JPLIB,CV(:8)//NORD(IXSPER),ILONG,ITYLCM)
               IF(ILONG.EQ.0) GO TO 10
               EXIST=.TRUE.
               GA1(:NGROUP+1)=0.0
               CALL LCMGET(JPLIB,CV(:8)//NORD(IXSPER),GA1)
               DO LLL=1,NGROUP+1
                  ESTOP(IBM,LLL,ICV)=ESTOP(IBM,LLL,ICV)+GA1(LLL)*DENISO
               ENDDO
               DENISO=DENISO*TIMFCT
            ENDDO
*-
   10       CONTINUE
   20       CONTINUE
         ENDIF
   30    CONTINUE
      ENDIF
   40 CONTINUE
      DEALLOCATE(GA1)
*----
*  APPLY STERNHEIMER DENSITY CORRECTION
*----
      STERNCORR=.TRUE.
      IF (STERNCORR) THEN 
         ALLOCATE(ISONRF(NBISO),ENER(NGROUP+1))
         CALL LCMGTC(IPLIB,'ISOTOPERNAME',12,NBISO,ISONRF)
         CALL LCMGET(IPLIB,'ENERGY',ENER)
         CALL LCMGET(IPLIB,'MIXTUREGAS',KGAS)
         CALL LIBSDC(NBMIX,NGROUP,NBISO,ISONRF,MIX,DEN,MASK,ENER,KGAS,
     1   DENMAT)
         DO IBM=1,NBMIX
            DO LLL=1,NGROUP+1
               PRINT *,ENER(LLL)
               !PRINT *,"ESTOP(IBM,LLL,1)=",ESTOP(IBM,LLL,1)
               !PRINT *,"DENMAT(IBM,LLL )=",DENMAT(IBM,LLL)
               ESTOP(IBM,LLL,1)=ESTOP(IBM,LLL,1)-DENMAT(IBM,LLL) !MeV/cm
               !PRINT *,"ESTOP(IBM,LLL,1)=",ESTOP(IBM,LLL,1)
            ENDDO
         ENDDO
         DEALLOCATE(ENER,ISONRF)
      ENDIF
*----
*  SAVE STOPPING POWERS IN THE MACROLIB
*----
      IF(EXIST) THEN
         DO LLL=1,NGROUP+1
            GAF(:NBMIX,LLL,1)=ESTOP(:NBMIX,LLL,1)+ESTOP(:NBMIX,LLL,2)
         ENDDO
         DO LLL=1,NGROUP
         IF(MASKL(LLL)) THEN
            KPLIB=IPGRP(LLL)
            CALL LCMLEN(KPLIB,'ESTOPW',ILONG,ITYLCM)
            IF(ILONG.GT.0) THEN
               CALL XDRSET(GAF(1,LLL,2),NBMIX,0.0)
               CALL XDRSET(GAF(1,LLL+1,2),NBMIX,0.0)
               CALL LCMGET(KPLIB,'ESTOPW',GAF(1,LLL,2))
               DO IBM=1,NBMIX
                  IF(.NOT.MASK(IBM)) THEN
                     GAF(IBM,LLL:LLL+1,1)=GAF(IBM,LLL:LLL+1,2)
                  ENDIF
               ENDDO
            ENDIF
            CALL LCMPUT(KPLIB,'ESTOPW',NBMIX*2,2,GAF(1,LLL:LLL+1,1))
            !CALL LCMPUT(KPLIB,'STERNHEIMER',NBMIX*1,2,DENMAT(1,LLL))
            ENDIF
         ENDDO
      ENDIF

!AHMED
*----
*  RECOVER AD MIXTURES DENSITIES
*----
      DO IBM=1,NBMIX
         CALL LIBCON(IPLIB,IBM,NBISO,MIX,DEN,DENMIXR,2)
         DENMIX(IBM)=DENMIXR !g/cm3
      ENDDO
      CALL LCMPUT(KPLIB,'DENMIXTURES',NBMIX,2,DENMIX)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
   50 DEALLOCATE(IPGRP)
      DEALLOCATE(GAF,ESTOP)
      DEALLOCATE(KGAS)
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(' WARNING IN LIBEST FOR PERTURBATION'/
     >       ' EXTRAPOLATION BELOW PRETURBATION TABLES'/
     >       ' INITIAL TIME       = ',F15.6,' DAYS'/
     >       ' EXTRAPOLATION TIME = ',F15.6,' DAYS')
 6001 FORMAT(' WARNING IN LIBEST FOR PERTURBATION'/
     >       ' EXTRAPOLATION ABOVE PRETURBATION TABLES'/
     >       ' FINAL TIME         = ',F15.6,' DAYS'/
     >       ' EXTRAPOLATION TIME = ',F15.6,' DAYS')
      END
