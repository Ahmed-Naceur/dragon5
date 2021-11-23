*DECK PCRONE
      SUBROUTINE PCRONE(IMPX,ICAL,IPMIC,NCAL,NGRP,XS_CALC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Copy an elementary calculation of the PMAXS file into a microlib.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert and D. Calic
*
*Parameters: input
* IMPX    print parameter.
* ICAL    index of the elementary calculation.
* IPMIC   address of the microlib.
* NCAL    number of elementary calculations in the PMAXS file.
* NGRP    number of energy groups.
* XS_CALC pointers towards PMAXS elementary calculations.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE PCRDATA
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,ICAL,NCAL,NGRP
      TYPE(C_PTR) IPMIC
      TYPE(XSBLOCK_ITEM),TARGET :: XS_CALC(NCAL)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMIC,KPMIC
      INTEGER NSTATE
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE),ITYPR(1)
      REAL DENS(3)
      DOUBLE PRECISION DELTA
      CHARACTER(LEN=8) :: HVECT(3)
      CHARACTER(LEN=12) :: HNAME,HISONA(3)
      TYPE(XSBLOCK_TYPE),POINTER :: XSONE,XSREF
      TYPE(TH_INDEP_VAR),POINTER :: TIVONE
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: SIG1,GAR,ENERGY
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIG2
*----
*  SET SIGNATURE AND STATE VECTOR
*----
      OPEN(UNIT=39,FILE='PCRONE.txt',STATUS='UNKNOWN')

      NED=1
      HVECT(1)='H-FACTOR'
      IF(NXST.GE.7) THEN
         NED=2
         HVECT(2)='NFTOT'
      ENDIF
      IF(NXST.EQ.8) THEN
         NED=3
         HVECT(3)='DETEC'
      ENDIF
      NBISO=1 ! number of isotopes
      IF(NXST.GT.4) NBISO=3 ! include Xe and Sm
      HNAME='L_LIBRARY'
      CALL LCMPTC(IPMIC,'SIGNATURE',12,1,HNAME)
      ISTATE(:)=0
      ISTATE(1)=1
      ISTATE(2)=NBISO
      ISTATE(3)=NGRP
      ISTATE(4)=1 ! isotropic scattering
      ISTATE(13)=NED ! number of additional edits
      ISTATE(19)=NDLAY ! number of delayed neutron groups
      CALL LCMPUT(IPMIC,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPTC(IPMIC,'ADDXSNAME-P0',8,NED,HVECT)
      JPMIC=LCMLID(IPMIC,'ISOTOPESLIST',NBISO)
*
      ALLOCATE(SIG1(NGRP),SIG2(NGRP,NGRP))
      XSONE=>XS_CALC(ICAL)%XS
      XSREF=>XS_CALC(XS_CALC(ICAL)%IBURN)%XS
      WRITE(39,*)XS_CALC(ICAL)%IBURN
      TIVONE=>XS_CALC(ICAL)%TIV
      KPMIC=LCMDIL(JPMIC,1) ! step up isot 1
      HISONA(1)='*MAC*RES'
      DENS(1)=1.0
      CALL LCMPTC(KPMIC,'ALIAS',12,1,HISONA(1))
*----
*  PROCESS VECTORIAL CROSS SECTIONS
*----
*    1   2   3   4   5   6   7
*   xtr,xab,xnf,xkf,xfi,xxe,xsm
*----
*  RUN ELEMENTARY CALC IN CASE THE PMAXS FILE IS GIVEN AS DERIVATIVES
*----
      IF(derivatives) THEN
        DELTA=XS_CALC(ICAL)%DELTA
        DO IG=1,NGRP
          SIG1(IG)=REAL(XSREF%sig(IG,1)+DELTA*XSONE%sig(IG,1))
        ENDDO
        CALL LCMPUT(KPMIC,'STRD',NGRP,2,SIG1)
        DO IG=1,NGRP
          SIG1(IG)=REAL(XSREF%sig(IG,2)+SUM(XSREF%sct(IG,:))+DELTA*
     1    XSONE%sig(IG,2)+SUM(XSONE%sct(IG,:)))
        ENDDO
        CALL LCMPUT(KPMIC,'NTOT0',NGRP,2,SIG1)
        DO IG=1,NGRP
          SIG1(IG)=REAL(XSREF%sig(IG,3)+DELTA*XSONE%sig(IG,3))
          WRITE(39,*)SIG1(IG)
        ENDDO
        CALL LCMPUT(KPMIC,'NUSIGF',NGRP,2,SIG1)
        DO IG=1,NGRP
          SIG1(IG)=REAL(XSREF%sig(IG,4)+DELTA*XSONE%sig(IG,4))
        ENDDO
        CALL LCMPUT(KPMIC,'H-FACTOR',NGRP,2,SIG1)
        IF(NXST.GT.4) THEN
          DO IG=1,NGRP
            SIG1(IG)=REAL(XSREF%sig(IG,7)+DELTA*XSONE%sig(IG,7))
          ENDDO
          CALL LCMPUT(KPMIC,'NFTOT',NGRP,2,SIG1)
          SIG1(:NGRP)=0.0
          SIG1(1)=1.0
          CALL LCMPUT(KPMIC,'CHI',NGRP,2,SIG1)
        ENDIF
        IF(NXST.EQ.8) THEN
          DO IG=1,NGRP
            SIG1(IG)=REAL(XSREF%sig(IG,8)+DELTA*XSONE%sig(IG,8))
          ENDDO
          CALL LCMPUT(KPMIC,'DETEC',NGRP,2,SIG1)
        ENDIF
        IF(lamb) THEN
          ALLOCATE(GAR(ELAM-BLAM+1))
          GAR(:ELAM-BLAM+1)=REAL(TIVONE%kinp(BLAM:ELAM))
          CALL LCMPUT(KPMIC,'LAMBDA',ELAM-BLAM+1,2,GAR)
          DEALLOCATE(GAR)
        ENDIF
        DO j=1,iXSTI
          k=iTIV(j)
          IF(k.GT.0) THEN
            DO IG=1,NGRP
              SIG1(IG)=REAL(TIVONE%sig(IG,k))
            ENDDO
            IF(j.EQ.1) THEN
              CALL LCMPUT(KPMIC,'CHI',NGRP,2,SIG1)
            ELSE IF(j.EQ.2) THEN
              CALL LCMPUT(KPMIC,'CHID',NGRP,2,SIG1)
            ELSE IF(j.EQ.3) THEN
              CALL LCMPUT(KPMIC,'INVEL',NGRP,2,SIG1)
            ENDIF
          ENDIF
        ENDDO
*----
*  PROCESS SCATTERING INFORMATION (JG --> IG)
*----
        DO IG=1,NGRP
          SIG1(IG)=REAL(SUM(XSREF%sct(IG,:))+DELTA*SUM(XSONE%sct(IG,:)))
          DO JG=1,NGRP
            SIG2(IG,JG)=REAL(XSREF%sct(JG,IG)+DELTA*XSONE%sct(JG,IG))
          ENDDO
        ENDDO
        ITYPR(1)=1
        CALL XDRLGS(KPMIC,1,IMPX,0,0,1,NGRP,SIG1,SIG2,ITYPR)
        IF(IMPX.GT.5) CALL LCMLIB(KPMIC)
*----
*  PROCESS Xe and Sm
*----
        IF(NXST.GT.4) THEN
          KPMIC=LCMDIL(JPMIC,2) ! step up isot 2
          HISONA(2)='Xe135'
          DENS(2)=0.0
          CALL LCMPTC(KPMIC,'ALIAS',12,1,HISONA(2))
          DO IG=1,NGRP
            SIG1(IG)=REAL(XSREF%sig(IG,5)+DELTA*XSONE%sig(IG,5))
          ENDDO
          CALL LCMPUT(KPMIC,'NTOT0',NGRP,2,SIG1)
          CALL LCMPUT(KPMIC,'STRD',NGRP,2,SIG1)
          DO IG=1,NGRP
            SIG1(IG)=0.0
            DO JG=1,NGRP
              SIG2(IG,JG)=0.0
            ENDDO
          ENDDO
          ITYPR(1)=1
          CALL XDRLGS(KPMIC,1,IMPX,0,0,1,NGRP,SIG1,SIG2,ITYPR)
          KPMIC=LCMDIL(JPMIC,3) ! step up isot 3
          HISONA(3)='Sm149'
          DENS(3)=0.0
          CALL LCMPTC(KPMIC,'ALIAS',12,1,HISONA(3))
          DO IG=1,NGRP
            SIG1(IG)=REAL(XSREF%sig(IG,6)+DELTA*XSONE%sig(IG,6))
          ENDDO
          CALL LCMPUT(KPMIC,'NTOT0',NGRP,2,SIG1)
          CALL LCMPUT(KPMIC,'STRD',NGRP,2,SIG1)
          DO IG=1,NGRP
            SIG1(IG)=0.0
            DO JG=1,NGRP
              SIG2(IG,JG)=0.0
            ENDDO
          ENDDO
          ITYPR(1)=1
          CALL XDRLGS(KPMIC,1,IMPX,0,0,1,NGRP,SIG1,SIG2,ITYPR)
        ENDIF
        CALL LCMPTC(IPMIC,'ISOTOPESUSED',12,NBISO,HISONA)
        CALL LCMPUT(IPMIC,'ISOTOPESDENS',NBISO,2,DENS)
        DEALLOCATE(SIG2,SIG1)
      ELSE
*----
*  RUN ELEMENTARY CALC IN CASE THE PMAXS FILE IS GIVEN AS RAW CROSS
*  SECTIONS
*----
        DELTA=XS_CALC(ICAL)%DELTA
        DO IG=1,NGRP
          SIG1(IG)=REAL(XSONE%sig(IG,1))
        ENDDO
        CALL LCMPUT(KPMIC,'STRD',NGRP,2,SIG1)
        DO IG=1,NGRP
          SIG1(IG)=REAL(XSONE%sig(IG,2)+SUM(XSONE%sct(IG,:)))
        ENDDO
        CALL LCMPUT(KPMIC,'NTOT0',NGRP,2,SIG1)
        DO IG=1,NGRP
          SIG1(IG)=REAL(XSONE%sig(IG,3))
          WRITE(39,*)XSONE%sig(IG,3)
        ENDDO
        CALL LCMPUT(KPMIC,'NUSIGF',NGRP,2,SIG1)
        DO IG=1,NGRP
          SIG1(IG)=REAL(XSONE%sig(IG,4))
        ENDDO
        CALL LCMPUT(KPMIC,'H-FACTOR',NGRP,2,SIG1)
        IF(NXST.GT.4) THEN
          DO IG=1,NGRP
            SIG1(IG)=REAL(XSONE%sig(IG,7))
          ENDDO
          CALL LCMPUT(KPMIC,'NFTOT',NGRP,2,SIG1)
          SIG1(:NGRP)=0.0
          SIG1(1)=1.0
          CALL LCMPUT(KPMIC,'CHI',NGRP,2,SIG1)
        ENDIF
        IF(NXST.EQ.8) THEN
          DO IG=1,NGRP
            SIG1(IG)=REAL(XSONE%sig(IG,8))
          ENDDO
          CALL LCMPUT(KPMIC,'DETEC',NGRP,2,SIG1)
        ENDIF
        IF(lamb) THEN
          ALLOCATE(GAR(ELAM-BLAM+1))
          GAR(:ELAM-BLAM+1)=REAL(TIVONE%kinp(BLAM:ELAM))
          CALL LCMPUT(KPMIC,'LAMBDA',ELAM-BLAM+1,2,GAR)
          DEALLOCATE(GAR)
        ENDIF
        DO j=1,iXSTI
          k=iTIV(j)
          IF(k.GT.0) THEN
            DO IG=1,NGRP
              SIG1(IG)=REAL(TIVONE%sig(IG,k))
            ENDDO
            IF(j.EQ.1) THEN
              CALL LCMPUT(KPMIC,'CHI',NGRP,2,SIG1)
            ELSE IF(j.EQ.2) THEN
              CALL LCMPUT(KPMIC,'CHID',NGRP,2,SIG1)
            ELSE IF(j.EQ.3) THEN
              CALL LCMPUT(KPMIC,'INVEL',NGRP,2,SIG1)
            ENDIF
          ENDIF
        ENDDO
*----
*  PROCESS SCATTERING INFORMATION
*----
        DO IG=1,NGRP
           SIG1(IG)=REAL(SUM(XSONE%sct(IG,:)))
          DO JG=1,NGRP
            SIG2(IG,JG)=REAL(XSONE%sct(JG,IG)) ! JG --> IG
          ENDDO
        ENDDO
        ITYPR(1)=1
        CALL XDRLGS(KPMIC,1,IMPX,0,0,1,NGRP,SIG1,SIG2,ITYPR)
        IF(IMPX.GT.5) CALL LCMLIB(KPMIC)
*----
*  PROCESS Xe and Sm
*----
        IF(NXST.GT.4) THEN
          KPMIC=LCMDIL(JPMIC,2) ! step up isot 2
          HISONA(2)='Xe135'
          DENS(2)=0.0
          CALL LCMPTC(KPMIC,'ALIAS',12,1,HISONA(2))
          DO IG=1,NGRP
            SIG1(IG)=REAL(XSONE%sig(IG,5))
          ENDDO
          CALL LCMPUT(KPMIC,'NTOT0',NGRP,2,SIG1)
          CALL LCMPUT(KPMIC,'STRD',NGRP,2,SIG1)
          DO IG=1,NGRP
            SIG1(IG)=0.0
            DO JG=1,NGRP
              SIG2(IG,JG)=0.0
            ENDDO
          ENDDO
          ITYPR(1)=1
          CALL XDRLGS(KPMIC,1,IMPX,0,0,1,NGRP,SIG1,SIG2,ITYPR)
          KPMIC=LCMDIL(JPMIC,3) ! step up isot 3
          HISONA(3)='Sm149'
          DENS(3)=0.0
          CALL LCMPTC(KPMIC,'ALIAS',12,1,HISONA(3))
          DO IG=1,NGRP
            SIG1(IG)=REAL(XSONE%sig(IG,6))
          ENDDO
          CALL LCMPUT(KPMIC,'NTOT0',NGRP,2,SIG1)
          CALL LCMPUT(KPMIC,'STRD',NGRP,2,SIG1)
          DO IG=1,NGRP
            SIG1(IG)=0.0
            DO JG=1,NGRP
              SIG2(IG,JG)=0.0
            ENDDO
          ENDDO
          ITYPR(1)=1
          CALL XDRLGS(KPMIC,1,IMPX,0,0,1,NGRP,SIG1,SIG2,ITYPR)
        ENDIF
        CALL LCMPTC(IPMIC,'ISOTOPESUSED',12,NBISO,HISONA)
        CALL LCMPUT(IPMIC,'ISOTOPESDENS',NBISO,2,DENS)
        DEALLOCATE(SIG2,SIG1)
      ENDIF
*----
*  SET ENERGY MESH
*----
      ALLOCATE(ENERGY(NGRP+1))
      IF(NGRP.EQ.2) THEN
        ENERGY(:)=(/ 1.964E7,  6.25E-1,  1.1E-4 /)
      ELSE
        CALL XABORT('PCRONE: UNKNOWN ENERGY MESH')
      ENDIF
      CALL LCMPUT(IPMIC,'ENERGY',NGRP+1,2,ENERGY)
      DEALLOCATE(ENERGY)
      RETURN
      END
