*DECK SPHBIV
      SUBROUTINE SPHBIV(IPTRK2,JPSYS,NPSYS,NREG,NUN,NMERGE,NALBP,NGCOND,
     1 NBMIX,MAT,VOL,KEY,MERG,SUNMER,FLXMER,SPH,IEX,FUNKNO,SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Diffusion theory calculation with BIVAC over the macro-geometry.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPTRK2  pointer to the tracking of the macro-geometry
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
      LOGICAL     CYLIND 
      CHARACTER   HSMG*131
      TYPE(C_PTR) KPSYS
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KN,IQFR,IPERT
      REAL, ALLOCATABLE, DIMENSION(:) :: XX,YY,DD,QFR,R,RS,T,TS,RH,RT,
     1 DIFF,GAMMA
*----
*  RECOVER BIVAC SPECIFIC TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK2,'STATE-VECTOR',IPAR)
      ITYPE=IPAR(6)
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
      IELEM=IPAR(8)
      ICOL=IPAR(9)
      ISPLH=IPAR(10)
      LL4=IPAR(11)
      LX=IPAR(12)
      LY=IPAR(13)
      ISCAT=IPAR(16)
      IF(ISCAT.GT.1) THEN
         WRITE(HSMG,'(45HSPHBIV: MACRO-CALCULATION WITH ANISOTROPIC SC,
     1   58HATTERING CURRENTLY NOT IMPLEMENTED; USE SCAT 1 KEYWORD IN ,
     2   13HBIVACT: DATA.)')
         CALL XABORT(HSMG)
      ENDIF
      CALL LCMLEN(IPTRK2,'KN',MAXKN,ITYLCM)
      CALL LCMLEN(IPTRK2,'QFR',MAXQF,ITYLCM)
      ALLOCATE(XX(LX*LY),YY(LX*LY),DD(LX*LY),KN(MAXKN),QFR(MAXQF),
     1 IQFR(MAXQF))
      IF(ITYPE.EQ.8) THEN
         CALL LCMGET(IPTRK2,'SIDE',SIDE)
      ELSE
         CALL LCMGET(IPTRK2,'XX',XX)
         CALL LCMGET(IPTRK2,'YY',YY)
         CALL LCMGET(IPTRK2,'DD',DD)
      ENDIF
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
      IF(MAT(I).GT.NBMIX) CALL XABORT('SPHBIV: NBMIX OVERFLOW.')
   10 CONTINUE
      IF((IELEM.LT.0).AND.(ITYPE.NE.8)) THEN
*        NON-HEXAGONAL GEOMETRY WITH PRIMAL DIFFUSION DISCRETIZATION
         CALL LCMSIX(IPTRK2,'BIVCOL',1)
         CALL LCMLEN(IPTRK2,'T',LC,ITYLCM)
         ALLOCATE(R(LC*LC),RS(LC*LC),T(LC),TS(LC))
         CALL LCMGET(IPTRK2,'R',R)
         CALL LCMGET(IPTRK2,'RS',RS)
         CALL LCMGET(IPTRK2,'T',T)
         CALL LCMGET(IPTRK2,'TS',TS)
         CALL LCMSIX(IPTRK2,' ',2)
         CALL BIVFSO(IEX,MAXKN,NGCOND,NMERGE,NALBP,FLXMER,SPH,
     1   SUNMER(1,1,IGR),CYLIND,NREG,NUN,NBMIX,XX,DD,MAT,KN,QFR,VOL,
     2   MERG,COUR,FUNKNO,LC,T,TS,R,RS,SUNKNO(1,IGR))
         DEALLOCATE(TS,T,RS,R)
      ELSE IF((IELEM.GT.0).AND.(ITYPE.NE.8)) THEN
*        NON-HEXAGONAL GEOMETRY WITH THOMAS-RAVIART DISCRETIZATION
         CALL BIVGSO(IEX,NGCOND,NMERGE,NALBP,FLXMER,SPH,SUNMER(1,1,IGR),
     1   CYLIND,IELEM,NREG,NUN,NBMIX,MAT,XX,DD,KN,
     2   QFR,VOL,MERG,COUR,FUNKNO,SUNKNO(1,IGR))
      ELSE IF((IELEM.LT.0).AND.(ITYPE.EQ.8)) THEN
*        HEXAGONAL GEOMETRY WITH PRIMAL DIFFUSION DISCRETIZATION
         CALL LCMSIX(IPTRK2,'BIVCOL',1)
         ALLOCATE(T(2),RH(36),RT(9))
         CALL LCMGET(IPTRK2,'T',T)
         CALL LCMGET(IPTRK2,'RH',RH)
         CALL LCMGET(IPTRK2,'RT',RT)
         CALL LCMSIX(IPTRK2,' ',2)
         IF(ISPLH.EQ.1) THEN
            NELEM=MAXKN/7
         ELSE
            NELEM=MAXKN/4
         ENDIF
         CALL BIVFSH(IEX,MAXKN,MAXQF,NGCOND,NMERGE,NALBP,FLXMER,SPH,
     1   SUNMER(1,1,IGR),NREG,NUN,ISPLH,NELEM,NBMIX,SIDE,MAT,
     2   KN,QFR,VOL,MERG,COUR,FUNKNO,T,RH,RT,SUNKNO(1,IGR))
         DEALLOCATE(T,RT,RH)
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.8).AND.(ICOL.EQ.4)) THEN
*        HEXAGONAL GEOMETRY WITH MCFD DISCRETIZATION
         ALLOCATE(DIFF(NBMIX+1))
         CALL LCMGET(KPSYS,'DRAGON-DIFF',DIFF)
         CALL BIVGSH(IEX,MAXKN,MAXQF,NGCOND,NMERGE,NALBP,FLXMER,SPH,
     1   SUNMER(1,1,IGR),NREG,NUN,LL4,ISPLH,NBMIX,SIDE,MAT,KN,QFR,
     2   VOL,DIFF(2),MERG,COUR,FUNKNO,SUNKNO(1,IGR))
         DEALLOCATE(DIFF)
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.8)) THEN
*        HEXAGONAL GEOMETRY WITH THOMAS-RAVIART-SCHNEIDER DISCRETIZATION
         NBLOS=LX/3
         ALLOCATE(IPERT(NBLOS))
         CALL LCMGET(IPTRK2,'IPERT',IPERT)
         CALL BIVGSR(IEX,NGCOND,NMERGE,NALBP,FLXMER,SPH,SUNMER(1,1,IGR),
     1   IELEM,NBLOS,NUN,NBMIX,SIDE,MAT,IPERT,KN,QFR,MERG,COUR,FUNKNO,
     2   SUNKNO(1,IGR))
         DEALLOCATE(IPERT)
      ELSE
         CALL XABORT('SPHBIV: UNKNOWN TYPE OF GEOMETRY(1).')
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
*  RELEASE BIVAC SPECIFIC TRACKING INFORMATION.
*----
      DEALLOCATE(IQFR,QFR,KN,DD,YY,XX)
      RETURN
      END
