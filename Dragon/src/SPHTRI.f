*DECK SPHTRI
      SUBROUTINE SPHTRI(IPTRK2,JPSYS,NPSYS,NREG,NUN,NMERGE,NALBP,NGCOND,
     1 NBMIX,MAT,VOL,KEY,MERG,SUNMER,FLXMER,SPH,IEX,FUNKNO,SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Diffusion theory calculation with TRIVAC over the macro-geometry.
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
* IPTRK2  pointer to the TRIVAC tracking of the macro-geometry
*         (L_TRACK signature).
* JPSYS   pointer to the 'GROUP' entry in the system LCM object.
* NPSYS   group masks.
* NREG    number of macro-regions (in the macro calculation).
* NUN     number of unknowns in the macro-calculation.
* NMERGE  number of merged regions.
* NALBP   number of physical albedos.
* NGCOND  number of condensed groups.
* NBMIX   number of macro-mixtures.
* MAT     mixture index per macro-region.
* VOL     volume of macro-regions.
* KEY     position of the flux components associated with each volume.
* MERG    index of merged macro-regions per macro-mixture.
* SUNMER  incoming source (scattering+fission) cross sections.
* FLXMER  flux estimate per mixture.
* SPH     SPH factors.
* IEX     iteration index.
* FUNKNO  neutron flux.
*
*Parameters: output
* SUNKNO  neutron sources.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK2,JPSYS
      INTEGER NPSYS(NGCOND),NREG,NUN,NMERGE,NALBP,NGCOND,NBMIX,
     1 MAT(NREG),KEY(NREG),MERG(NBMIX),IEX
      REAL VOL(NREG),SUNMER(NMERGE,NGCOND,NGCOND),FLXMER(NMERGE,NGCOND),
     1 SPH(NMERGE+NALBP,NGCOND),FUNKNO(NUN,NGCOND),SUNKNO(NUN,NGCOND)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (NSTATE=40)
      INTEGER     IPAR(NSTATE)
      LOGICAL     CYLIND,LHEX
      CHARACTER   HSMG*131
      TYPE(C_PTR) KPSYS
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KN,IQFR,IPERT
      REAL, ALLOCATABLE, DIMENSION(:) :: XX,YY,ZZ,DD,QFR,FRZ,GAMMA
*----
*  RECOVER TRIVAC SPECIFIC TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK2,'STATE-VECTOR',IPAR)
      ITYPE=IPAR(6)
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
      LHEX=(ITYPE.EQ.8).OR.(ITYPE.EQ.9)
      IELEM=IPAR(9)
      ICOL=IPAR(10)
      LL4=IPAR(11)
      ICHX=IPAR(12)
      ISPLH=IPAR(13)
      LX=IPAR(14)
      LY=IPAR(15)
      LZ=IPAR(16)
      ISCAT=IPAR(32)
      IF(ISCAT.GT.1) THEN
         WRITE(HSMG,'(45HSPHTRI: MACRO-CALCULATION WITH ANISOTROPIC SC,
     1   58HATTERING CURRENTLY NOT IMPLEMENTED; USE SCAT 1 KEYWORD IN ,
     2   13HTRIVAT: DATA.)')
         CALL XABORT(HSMG)
      ENDIF
      CALL LCMLEN(IPTRK2,'KN',MAXKN,ITYLCM)
      CALL LCMLEN(IPTRK2,'QFR',MAXQF,ITYLCM)
      ALLOCATE(XX(LX*LY*LZ),YY(LX*LY*LZ),ZZ(LX*LY*LZ),DD(LX*LY*LZ),
     1 KN(MAXKN),QFR(MAXQF),IQFR(MAXQF))
      IF(LHEX) THEN
         CALL LCMGET(IPTRK2,'SIDE',SIDE)
      ELSE
         CALL LCMGET(IPTRK2,'XX',XX)
         CALL LCMGET(IPTRK2,'YY',YY)
         CALL LCMGET(IPTRK2,'DD',DD)
      ENDIF
      CALL LCMGET(IPTRK2,'ZZ',ZZ)
      CALL LCMGET(IPTRK2,'KN',KN)
      CALL LCMGET(IPTRK2,'QFR',QFR)
      CALL LCMGET(IPTRK2,'IQFR',IQFR)
*----
*  LOOP OVER MACRO ENERGY GROUPS.
*----
      DO 30 IGR=1,NGCOND
      CALL XDRSET(SUNKNO(1,IGR),NUN,0.0)
      IF(NPSYS(IGR).EQ.0) GO TO 30
      KPSYS=LCMDIL(JPSYS,IGR)
*----
*  APPLY PHYSICAL ALBEDOS.
*----
      CALL LCMLEN(KPSYS,'ALBEDO-FU',NALBP,ITYLCM)
      IF(NALBP.GT.0) THEN
         ALLOCATE(GAMMA(NALBP))
         CALL LCMGET(KPSYS,'ALBEDO-FU',GAMMA)
         DO IQW=1,MAXQF
            IALB=IQFR(IQW)
            IF(IALB.NE.0) QFR(IQW)=QFR(IQW)*GAMMA(IALB)
         ENDDO
         DEALLOCATE(GAMMA)
      ENDIF
*----
*  COMPUTE THE NEUTRON SOURCES IN GROUP IGR.
*----
      COUR=0.0
      DO 10 I=1,NREG
      IF(MAT(I).GT.NBMIX) CALL XABORT('SPHTRI: NBMIX OVERFLOW.')
   10 CONTINUE
      IF((IELEM.GT.0).AND.(ICHX.EQ.2).AND.(.NOT.LHEX)) THEN
*        NON-HEXAGONAL GEOMETRY WITH THOMAS-RAVIART DISCRETIZATION
         CALL TRIGSO(IEX,NGCOND,NMERGE,NALBP,FLXMER,SPH,SUNMER(1,1,IGR),
     1   CYLIND,IELEM,NREG,NUN,NBMIX,MAT,XX,DD,KN,QFR,VOL,MERG,COUR,
     2   FUNKNO,SUNKNO(1,IGR))
      ELSE IF((IELEM.GT.0).AND.(ICHX.EQ.2).AND.LHEX) THEN
*        HEXAGONAL GEOMETRY WITH THOMAS-RAVIART-SCHNEIDER DISCRETIZATION
         NBLOS=LX*LZ/3
         ALLOCATE(IPERT(NBLOS),FRZ(NBLOS))
         CALL LCMGET(IPTRK2,'IPERT',IPERT)
         CALL LCMGET(IPTRK2,'FRZ',FRZ)
         CALL TRIGSR(IEX,NGCOND,NMERGE,NALBP,FLXMER,SPH,SUNMER(1,1,IGR),
     1   IELEM,NBLOS,NUN,NBMIX,SIDE,ZZ,FRZ,MAT,IPERT,KN,QFR,MERG,COUR,
     2   FUNKNO,SUNKNO(1,IGR))
         DEALLOCATE(FRZ,IPERT)
      ELSE
         CALL XABORT('SPHTRI: THIS DISCRETIZATION IS NOT AVAILABLE.')
      ENDIF
*----
*  DIVISION OF THE EVEN-PARITY SOURCES BY THE VOLUMES.
*----
      DO 20 I=1,NREG
      IF(MAT(I).EQ.0) GO TO 20
      JND1=KEY(I)
      SUNKNO(JND1,IGR)=SUNKNO(JND1,IGR)/VOL(I)
   20 CONTINUE
   30 CONTINUE
*----
*  RELEASE TRIVAC SPECIFIC TRACKING INFORMATION.
*----
      DEALLOCATE(IQFR,QFR,KN,DD,ZZ,YY,XX)
      RETURN
      END
