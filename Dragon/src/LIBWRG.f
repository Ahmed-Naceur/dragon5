*DECK LIBWRG
      SUBROUTINE LIBWRG(IUNIT,NTYP,NGR,NRTOT,MAXTEM,MAXDIL,NSRES,RID,
     >                  NTM,NDI,RTMP,RDIL,RESI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read resonance information from WIMS-D4 library.
*
*Copyright:
* Copyright (C) 1997 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): 
* G. Marleau
*
*Parameters: input
* IUNIT   WIMS-D4 read unit.                    
* NTYP    number of resonance tables per isotopes.
* NGR     number of resonance groups.
* NRTOT   number of resonance sets.
* MAXTEM  max nb temperature.
* MAXDIL  max nb dilutions.
* NSRES   nb of resonance set.
* RID     resonance id.
* NTM     number of temperatures.
* NDI     number of dilutions.                  
* RTMP    resonance temperature.
* RDIL    resonance dilution.
* RESI    resonance integrals.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
* PARAMETERS
*----
      INTEGER    IOUT
      PARAMETER (IOUT=6)
*----
* INTERFACE PARAMETERS
*----
      INTEGER    IUNIT,NTYP,NGR,NRTOT,MAXTEM,MAXDIL
      INTEGER    NTM(NTYP,NRTOT,NGR),NDI(NTYP,NRTOT,NGR)
*
      REAL       RID(NRTOT),RTMP(MAXTEM,NTYP,NRTOT,NGR),
     1           RDIL(MAXDIL,NTYP,NRTOT,NGR),
     2           RESI(MAXDIL,MAXTEM,NTYP,NRTOT,NGR)
*----
* LOCAL VARIABLES
*----
      INTEGER    IGR,NSRES,ISRES,IPREV,IRS,M1,M2,IT,ID,ISR,ITYP,
     1           NTIS
      REAL       XIDR,ENDR
*----
* ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: TMPT,DILT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: REST
*----
*  SCRATCH STORAGE ALLOCATION
*     TMPT   : TEMPERATURE
*     DILT   : DILUTION
*     REST   : RESONANCE INTEGRALS
*----
      ALLOCATE(TMPT(MAXTEM),DILT(MAXDIL),REST(MAXDIL,MAXTEM))
*----
*  SCAN OVER RESONANCE GROUPS
*----
      NSRES=0
      ISRES=0
      DO 100 IGR=1,NGR
        IPREV=0
*----
*  SCAN OVER RESONANCE SETS + 1
*  AND READ RESONANCE INFO
*----
        DO 110 IRS=1,NTYP*NRTOT+1
          READ(IUNIT) XIDR,M1,M2,
     >     (TMPT(IT),IT=1,M1),(DILT(ID),ID=1,M2),
     >    ((REST(ID,IT),ID=1,M2),IT=1,M1)
          IF(XIDR.EQ.0.0) GO TO 115
          IF((M1.EQ.0).AND.(M2.EQ.0)) GO TO 110
          DO 120 ISR=1,NSRES
            IF(XIDR.EQ.RID(ISR)) THEN
              ISRES=ISR
              GO TO 125
            ENDIF
 120      CONTINUE
          NSRES=NSRES+1
          IF(NSRES.GT.NRTOT) THEN
            CALL XABORT('LIBWRG: TO MANY RESONANCE SET')
          ENDIF
          ISRES=NSRES
          IPREV=0
          RID(ISRES)=XIDR
 125      CONTINUE
          IF(ISRES.NE.IPREV) THEN
            ITYP=1
            IPREV=ISRES
          ELSE IF((ISRES.EQ.IPREV).AND.(ITYP.EQ.1)) THEN
            ITYP=2
          ELSE IF((ISRES.EQ.IPREV).AND.(ITYP.EQ.2)) THEN
            ITYP=3
            IPREV=0
          ENDIF
          NTIS=NTM(ITYP,ISRES,IGR)
          IF(NTIS.GT.0) THEN
            WRITE(IOUT,9000) IGR,ISRES,ITYP,XIDR
            CALL XABORT('LIBWRG: DUPLICATE RESONANCE SET')
          ENDIF
*----
*  SAVE RESONANCE INFORMATION FOR THIS SET
*----
          NTM(ITYP,ISRES,IGR)=M1
          NDI(ITYP,ISRES,IGR)=M2
          DO 130 IT=1,M1
            RTMP(IT,ITYP,ISRES,IGR)=TMPT(IT)
 130      CONTINUE
          DO 131 ID=1,M2
            RDIL(ID,ITYP,ISRES,IGR)=DILT(ID)
 131      CONTINUE
          DO 140 IT=1,M1
            DO 141 ID=1,M2
              RESI(ID,IT,ITYP,ISRES,IGR)=REST(ID,IT)
 141        CONTINUE
 140      CONTINUE
 110    CONTINUE
 115    CONTINUE
        IF(NTYP.EQ.2) READ(IUNIT) ENDR
 100  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(REST,DILT,TMPT)
      RETURN
*----
*  FORMAT
*----
 9000 FORMAT(' LIBWRG ERROR - WIMS-D4 DUPLICATE RESONANCE SET'/
     >       ' RESONANCE GROUP = ',I10/
     >       '   RESONANCE SET = ',I10/
     >       '   INTEGRAL TYPE = ',I10/
     >       '    RESONANCE ID = ',F20.5)
      END
