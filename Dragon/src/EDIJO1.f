*DECK EDIJO1
      SUBROUTINE EDIJO1(IPMAC2,IPTRK1,IPFLUX,IPRINT,NGCOND,IGCOND)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover ALBS information from last component of unknown array for use
* with SPH equivalence techniques. SYBILF compatible version. SYBILF is
* activated with ARM keyword in ASM: module.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
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
      INTEGER ISTATE(NSTATE),IPAR(16)
      CHARACTER CDOOR*12
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR) IFR_PTR,ALB_PTR,SUR_PTR,INUM_PTR,MIX_PTR
      INTEGER, POINTER, DIMENSION(:) :: IFR,INUM,MIX
      REAL, POINTER, DIMENSION(:) :: ALB,SUR
      REAL, ALLOCATABLE, DIMENSION(:) :: WORKD
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
      IF(CDOOR.NE.'SYBIL') CALL XABORT('EDIJO2: SYBIL TRACKING EXPECTE'
     > //'D.')
      CALL LCMGET(IPTRK1,'STATE-VECTOR',ISTATE)
      NREG=ISTATE(1)
      ITG=ISTATE(6)
      IF(ITG.NE.4) CALL XABORT('EDIJO1: JOUT OPTION FORBIDDEN.')
      NUNKNO=ISTATE(2)+ISTATE(9)
      CALL LCMSIX(IPTRK1,'EURYDICE',1)
      CALL LCMGET(IPTRK1,'PARAM',IPAR)
      IHEX=IPAR(1)
      MULTC=IPAR(2)
      NMCEL=IPAR(4)
      NMERGE=IPAR(5)
      NCOUR=4
      IF(IHEX.NE.0) NCOUR=6
      IF(MULTC.EQ.4) NCOUR=3*NCOUR
      CALL LCMGPD(IPTRK1,'IFR',IFR_PTR)
      CALL LCMGPD(IPTRK1,'ALB',ALB_PTR)
      CALL LCMGPD(IPTRK1,'SUR',SUR_PTR)
      CALL LCMGPD(IPTRK1,'INUM',INUM_PTR)
      CALL LCMGPD(IPTRK1,'MIX',MIX_PTR)
      CALL LCMSIX(IPTRK1,' ',2)
*
      CALL C_F_POINTER(IFR_PTR,IFR,(/ NCOUR*NMCEL /))
      CALL C_F_POINTER(ALB_PTR,ALB,(/ NCOUR*NMCEL /))
      CALL C_F_POINTER(SUR_PTR,SUR,(/ NCOUR*NMCEL /))
      CALL C_F_POINTER(INUM_PTR,INUM,(/ NMCEL /))
      CALL C_F_POINTER(MIX_PTR,MIX,(/ NCOUR*NMERGE /))
*----
*  COMPUTE THE OUTGOING CURRENT
*----
      ALLOCATE(OUTG(NGCOND,2))
      IGRFIN=0
      CALL LCMSIX(IPMAC2,'ADF',1)
      DO 70 IGRCD=1,NGCOND
      OUTG(IGRCD,:2)=0.0
      IGRDEB=IGRFIN+1
      IGRFIN=IGCOND(IGRCD)
      CALL LCMLEN(IPFLUX,'FLUX',ILON,ITYLCM)
      IF(ILON.EQ.0) CALL XABORT('EDIJO1: MISSING FLUX INFO(1).')
      JPFLUX=LCMGID(IPFLUX,'FLUX')
      DO 60 IGR=IGRDEB,IGRFIN
      CALL LCMLEL(JPFLUX,IGR,ILCMLN,ITYLCM)
      IF(ILCMLN.EQ.0) CALL XABORT('EDIJO1: MISSING FLUX INFO(2).')
      IF(ILEAK.LE.4) THEN
        IF(ILCMLN.NE.NUNKNO) CALL XABORT('EDIJO1: ARM KEYWORD MUST B'
     1  //'E SET IN ASM: MODULE(1).')
        ALLOCATE(WORKD(NUNKNO))
      ELSE IF(ILEAK.EQ.5) THEN
        IF(ILCMLN.NE.2*NUNKNO) CALL XABORT('EDIJO1: ARM KEYWORD MUST'
     1  //' BE SET IN ASM: MODULE(2).')
        ALLOCATE(WORKD(2*NUNKNO))
      ELSE
        CALL XABORT('EDIJO1: INVALID TYPE OF LEAKAGE.')
      ENDIF
      CALL LCMGDL(JPFLUX,IGR,WORKD)
      OUTC1=0.0
      OUTC2=0.0
      SURT=0.0
      IF(MULTC.EQ.1) THEN
        DO 20 ICEL=1,NMCEL
        IKK=INUM(ICEL)
        IF(IKK.EQ.0) GO TO 20
        IT0=NCOUR*(ICEL-1)
        DO 10 JC=1,NCOUR
        IF((IKK.EQ.IFR(IT0+JC)).AND.(SUR(IT0+JC).NE.0.0)) THEN
          J1=IFR(IT0+JC)
          OUTC1=OUTC1+WORKD(NREG+J1)*SUR(IT0+JC)
          OUTC2=OUTC2+WORKD(NREG+J1)*SUR(IT0+JC)*ALB(IT0+JC)
          SURT=SURT+SUR(IT0+JC)
        ENDIF
   10   CONTINUE
   20   CONTINUE
      ELSE
        ISTR=1
        IF((NCOUR.EQ.12).OR.(NCOUR.EQ.18)) ISTR=3
        DO 50 ICEL=1,NMCEL
        IKK=INUM(ICEL)
        IF(IKK.EQ.0) GO TO 50
        IT0=NCOUR*(ICEL-1)
        IT1=NCOUR*(IKK-1)
        DO 40 JC=1,NCOUR,ISTR
        IF((MIX(IT1+JC).EQ.IFR(IT0+JC)).AND.(SUR(IT0+JC).NE.0.0)) THEN
          J1=IFR(IT0+JC)
          OUTC1=OUTC1+WORKD(NREG+J1)*SUR(IT0+JC)
          OUTC2=OUTC2+WORKD(NREG+J1)*SUR(IT0+JC)*ALB(IT0+JC)
          SURT=SURT+SUR(IT0+JC)
        ENDIF
   40   CONTINUE
   50   CONTINUE
      ENDIF
      OUTG(IGRCD,1)=OUTG(IGRCD,1)+OUTC1/SURT
      OUTG(IGRCD,2)=OUTG(IGRCD,2)+OUTC2/SURT
      DEALLOCATE(WORKD)
   60 CONTINUE
   70 CONTINUE
      CALL LCMPUT(IPMAC2,'ALBS00',NGCOND*2,2,OUTG)
      IF(IPRINT.GT.3) THEN
         WRITE(6,900) (OUTG(IGR,1),IGR=1,NGCOND)
         WRITE(6,910) (OUTG(IGR,2),IGR=1,NGCOND)
         WRITE(6,'(/)')
      ENDIF
      CALL LCMSIX(IPMAC2,' ',2)
      DEALLOCATE(OUTG)
      RETURN
*
  900 FORMAT(/49H EDIJO1: OUT-CURRENTS (4J-/S) PER MACRO-GROUP ARE/
     > (1X,1P,10E13.5))
  910 FORMAT(/49H EDIJO1:  IN-CURRENTS (4J+/S) PER MACRO-GROUP ARE/
     > (1X,1P,10E13.5))
      END
