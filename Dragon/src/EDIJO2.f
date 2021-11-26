*DECK EDIJO2
      SUBROUTINE EDIJO2(IPMAC2,IPTRK1,IPFLUX,IPRINT,NGCOND,IGCOND)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover ALBS information from last component of unknown array for use
* with SPH equivalence techniques. MCCG compatible version.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMAC2  pointer to condensed macrolib information (L_MACROLIB
*         signature) built by EDI:.
* IPTRK1  pointer to the reference tracking object.
* IPFLUX  pointer to the reference solution (L_FLUX signature).
* IPRINT  print index.
* NGCOND  number of condensed groups.
* IGCOND  limit of condensed groups.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC2,IPTRK1,IPFLUX
      INTEGER IPRINT,NGCOND,IGCOND(NGCOND)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      TYPE(C_PTR) JPFLUX
      INTEGER ISTATE(NSTATE)
      REAL ALBEDO(6)
      CHARACTER CDOOR*12
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NZON,KEYSUR
      REAL, ALLOCATABLE, DIMENSION(:) :: WORKD,VOLSUR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: OUTG
*----
*  RECOVER FLUX OBJECT INFORMATION
*----
      CALL LCMGET(IPFLUX,'STATE-VECTOR',ISTATE)
      NUNKNO=ISTATE(2)
      ILEAK=ISTATE(7)
*----
*  RECOVER TRACKING INFORMATION
*----
      CALL LCMGTC(IPTRK1,'TRACK-TYPE',12,1,CDOOR)
      CALL LCMGET(IPTRK1,'STATE-VECTOR',ISTATE)
      IF(CDOOR.NE.'MCCG') CALL XABORT('EDIJO2: MCCG TRACKING EXPECTED.')
      NBVOL=ISTATE(1)
      NBSUR=ISTATE(5)
      IF(NBSUR.EQ.0) CALL XABORT('EDIJO2: NO BOUNDARY LEAKAGE.')
      ALLOCATE(VOLSUR(NBVOL+NBSUR),NZON(NBVOL+NBSUR),KEYSUR(NBSUR))
      CALL LCMGET(IPTRK1,'V$MCCG',VOLSUR)
      CALL LCMGET(IPTRK1,'NZON$MCCG',NZON)
      CALL LCMGET(IPTRK1,'KEYCUR$MCCG',KEYSUR)
      CALL LCMGET(IPTRK1,'ALBEDO',ALBEDO)
*----
*  COMPUTE THE OUTGOING CURRENT
*----
      ALLOCATE(OUTG(NGCOND,2))
      IGRFIN=0
      CALL LCMSIX(IPMAC2,'ADF',1)
      ALLOCATE(WORKD(NUNKNO))
      DO 30 IGRCD=1,NGCOND
      OUTG(IGRCD,:2)=0.0
      IGRDEB=IGRFIN+1
      IGRFIN=IGCOND(IGRCD)
      CALL LCMLEN(IPFLUX,'FLUX',ILON,ITYLCM)
      IF(ILON.EQ.0) CALL XABORT('EDIJO2: MISSING FLUX INFO(1).')
      JPFLUX=LCMGID(IPFLUX,'FLUX')
      DO 20 IGR=IGRDEB,IGRFIN
      CALL LCMLEL(JPFLUX,IGR,ILONG,ITYLCM)
      IF(ILONG.NE.NUNKNO) CALL XABORT('EDIJO2: MISSING FLUX INFO(2).')
      CALL LCMGDL(JPFLUX,IGR,WORKD)
      DO 10 IS=1,NBSUR
      IUN=KEYSUR(IS)
      IF(IUN.EQ.0) GO TO 10
      IAL=-NZON(NBVOL+IS)
      OUTG(IGRCD,1)=OUTG(IGRCD,1)+WORKD(IUN)*VOLSUR(NBVOL+IS)
      OUTG(IGRCD,2)=OUTG(IGRCD,2)+WORKD(IUN)*VOLSUR(NBVOL+IS)*
     > ALBEDO(IAL)
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE
      DEALLOCATE(WORKD)
      CALL LCMPUT(IPMAC2,'ALBS00',NGCOND*2,2,OUTG)
      IF(IPRINT.GT.3) THEN
         WRITE(6,900) (OUTG(IGR,1),IGR=1,NGCOND)
         WRITE(6,910) (OUTG(IGR,2),IGR=1,NGCOND)
         WRITE(6,'(/)')
      ENDIF
      CALL LCMSIX(IPMAC2,' ',2)
      DEALLOCATE(KEYSUR,NZON,VOLSUR)
      DEALLOCATE(OUTG)
      RETURN
*
  900 FORMAT(/49H EDIJO2: OUT-CURRENTS (4J-/S) PER MACRO-GROUP ARE/
     > (1X,1P,10E13.5))
  910 FORMAT(/49H EDIJO2:  IN-CURRENTS (4J+/S) PER MACRO-GROUP ARE/
     > (1X,1P,10E13.5))
      END
