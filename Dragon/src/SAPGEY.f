*DECK SAPGEY
      SUBROUTINE SAPGEY(IPSAP,IPEDIT,NISO,NMAC,NG,NMIL,NBISO,ICAL,NDFI,
     1 NISFS,NISPS,NISYS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To recover the fission yields of an elementary calculation.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPSAP   pointer to the Saphyb.
* IPEDIT  pointer to the edition object (L_EDIT signature).
* NISO    number of particularized isotopes.
* NMAC    number of macros.
* NG      number of condensed energy groups.
* NMIL    number of mixtures in the Saphyb.
* NBISO   number of isotopes in the condensed microlib of the edition
*         object. A given isotope may appear in many mixtures.
* ICAL    index of the current elementary calculation.
* NDFI    number of fissile isotopes producing fission products in
*         the edition object.
* NISFS   number of particularized fissile isotopes.
* NISPS   number of particularized fission products.
* NISYS   number of particularized fissile isotopes, fission products
*         and macros.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSAP,IPEDIT
      INTEGER NISO,NMAC,NG,NMIL,NBISO,ICAL,NDFI,NISFS,NISPS,NISYS
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXISO=800)
      TYPE(C_PTR) JPEDIT,KPEDIT
      CHARACTER TEXT8*8,TEXT12*12,DIRNAM*12,NOMISO(MAXISO)*8
      LOGICAL LGIMF
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIX,ITYPE,PIFI,ADRY
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISONAM
      REAL, ALLOCATABLE, DIMENSION(:) :: DEN,PYIELD,SIG,PFIRA
      REAL, ALLOCATABLE, DIMENSION(:,:) :: YLDS,FLUXES
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) ::  IPISO
*----
*  SCRATCH STORAGE ALLOCATION
*   PFIRA   fission rate.
*   ADRY    offset in YLDS array for fissile isotopes (positive) and
*           fission products (negative).
*----
      ALLOCATE(ISONAM(3,NBISO),MIX(NBISO),ITYPE(NBISO),PIFI(NDFI))
      ALLOCATE(YLDS(NISFS+NMAC,NISPS),DEN(NBISO),PYIELD(NDFI),
     1 FLUXES(NMIL,NG),SIG(NG),PFIRA(NBISO),ADRY(NISYS))
      ALLOCATE(IPISO(NBISO))
*----
*  RECOVER INFORMATION FROM THE 'CONTENU' DIRECTORY.
*----
      CALL LCMSIX(IPSAP,'contenu',1)
      IF(NISO.GT.0) CALL LCMGTC(IPSAP,'NOMISO',8,NISO,NOMISO)
      CALL LCMSIX(IPSAP,' ',2)
*----
*  MOVE TO THE 'CALC' DIRECTORY.
*----
      WRITE(DIRNAM,'(''calc'',I8)') ICAL
      CALL LCMSIX(IPSAP,DIRNAM,1)
*
      CALL LCMGET(IPEDIT,'ISOTOPESUSED',ISONAM)
      CALL LCMGET(IPEDIT,'ISOTOPESMIX',MIX)
      CALL LCMGET(IPEDIT,'ISOTOPESDENS',DEN)
      CALL LCMGET(IPEDIT,'ISOTOPESTYPE',ITYPE)
      CALL LIBIPS(IPEDIT,NBISO,IPISO)
*----
*  COMPUTE ARRAY ADRY.
*----
      ISF=0
      ISP=0
      CALL XDISET(ADRY,NISYS,0)
      DO 30 ISO=1,NISO
      DO 10 IBISO=1,NBISO
      WRITE(TEXT8,'(2A4)') (ISONAM(I0,IBISO),I0=1,2)
      IF(NOMISO(ISO).EQ.TEXT8) GO TO 20
   10 CONTINUE
      GO TO 30
   20 IF(ITYPE(IBISO).EQ.2) THEN
         ISF=ISF+1
         ADRY(ISO)=ISF
      ELSEIF(ITYPE(IBISO).EQ.3) THEN
         ISP=ISP+1
         ADRY(ISO)=-ISP
      ENDIF
   30 CONTINUE
      IF(NMAC.GT.0) ADRY(NISO+1)=ISF+1
      LGIMF=NISYS.GT.NISO
      IMF=0
      IF(LGIMF) IMF=ADRY(NISO+1)
*----
*  RECOVER THE NEUTRON FLUX.
*----
      CALL LCMSIX(IPEDIT,'MACROLIB',1)
      JPEDIT=LCMGID(IPEDIT,'GROUP')
      DO 40 IGR=1,NG
      KPEDIT=LCMGIL(JPEDIT,IGR)
      CALL LCMGET(KPEDIT,'FLUX-INTG',FLUXES(1,IGR))
   40 CONTINUE
      CALL LCMSIX(IPEDIT,' ',2)
*----
*  RECOVER THE FISSION RATES.
*----
      DO 65 IBISO=1,NBISO
      GAR=0.0
      IF(MIX(IBISO).EQ.0) GO TO 60
      KPEDIT=IPISO(IBISO)
      CALL LCMLEN(KPEDIT,'NFTOT',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
         CALL LCMGET(KPEDIT,'NFTOT',SIG)
         DO 50 IGR=1,NG
         GAR=GAR+FLUXES(MIX(IBISO),IGR)*DEN(IBISO)*SIG(IGR)
   50    CONTINUE
      ENDIF
   60 PFIRA(IBISO)=GAR
   65 CONTINUE
*----
*  LOOP OVER SAPHYB MIXTURES TO RECOVER THE FISSION YIELDS.
*----
      DO 140 IMIL=1,NMIL
      DO 75 IFP=1,NISPS
      DO 70 IFI=1,NISFS+NMAC
      YLDS(IFI,IFP)=0.0
   70 CONTINUE
   75 CONTINUE
      DO 130 IBISO=1,NBISO
      IF(MIX(IBISO).EQ.IMIL) THEN
         WRITE(TEXT12,'(3A4)') (ISONAM(I0,IBISO),I0=1,3)
         DO 80 ISO=1,NISO
         IISO=ISO
         IF(NOMISO(ISO).EQ.TEXT12(:8)) GO TO 90
   80    CONTINUE
         GO TO 130
   90    KPEDIT=IPISO(IBISO)
*
*        RECOVER THE FISSION YIELDS.
         CALL LCMLEN(KPEDIT,'PYIELD',ILONG,ITYLCM)
         IF((ILONG.GT.0).AND.(ILONG.EQ.NDFI)) THEN
            CALL LCMGET(KPEDIT,'PIFI',PIFI)
            CALL LCMGET(KPEDIT,'PYIELD',PYIELD)
         ELSE
            GO TO 130
         ENDIF
         IFP=-ADRY(IISO)
         IF(IFP.GT.0) THEN
*           Particular fission product found.
*           If exists in medium, find position in microlib
*           and search all fissiles.
            YLDW=0.0
            DO 120 IDFI=1,NDFI
            JBISO=PIFI(IDFI)
            IF(JBISO.EQ.0) GO TO 120
            IF(JBISO.GT.NBISO) CALL XABORT('SAPGEY: MIX OVERFLOW.')
            IF(MIX(JBISO).NE.IMIL) CALL XABORT('SAPGEY: BAD MIX.')
            WRITE(TEXT8,'(3A4)') (ISONAM(I0,JBISO),I0=1,2)
            DO 100 JSO=1,NISO
            JISO=JSO
            IF(NOMISO(JSO).EQ.TEXT8) GO TO 110
  100       CONTINUE
*           Mother isotope is in residual macro.
            YLDW=YLDW+PFIRA(JBISO)
            IF(IMF.EQ.0) CALL XABORT('SAPGEY: LGIMF IS FALSE.')
            YLDS(IMF,IFP)=YLDS(IMF,IFP)+PYIELD(IDFI)*PFIRA(JBISO)
            GO TO 120
*
*           Yield for selected isotopes.
  110       IFI=ADRY(JISO)
            IF(IFI.LE.0) CALL XABORT('SAPGEY: BAD ADRY.')
            YLDS(IFI,IFP)=PYIELD(IDFI)
  120       CONTINUE
            IF(LGIMF) THEN
               IF(YLDW.NE.0.0) YLDS(IMF,IFP)=YLDS(IMF,IFP)/YLDW
            ENDIF
         ENDIF
      ENDIF
  130 CONTINUE
*----
*  STORE INFORMATION IN THE MIXTURE DIRECTORY.
*----
      WRITE(DIRNAM,'(''mili'',I8)') IMIL
      CALL LCMSIX(IPSAP,DIRNAM,1)
      CALL LCMPUT(IPSAP,'YLDS',(NISFS+NMAC)*NISPS,2,YLDS)
      CALL LCMSIX(IPSAP,' ',2)
  140 CONTINUE
*
      CALL LCMSIX(IPSAP,'info',1)
      CALL LCMPUT(IPSAP,'ADRY',NISYS,1,ADRY)
      CALL LCMSIX(IPSAP,' ',2)
*----
*  MOVE TO THE SAPHYB ROOT DIRECTORY.
*----
      CALL LCMSIX(IPSAP,' ',2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IPISO)
      DEALLOCATE(ADRY,PFIRA,SIG,FLUXES,PYIELD,DEN,YLDS)
      DEALLOCATE(PIFI,ITYPE,MIX,ISONAM)
      RETURN
      END
