*DECK BIVGSH
      SUBROUTINE BIVGSH (IEX,MAXKN,MAXQF,NGCOND,NMERGE,NALBP,FLXMER,
     1 SPH,SUNMER,NREG,NUN,LL4,ISPLH,NBMIX,SIDE,MAT,KN,QFR,VOL,DIFF0,
     2 MERG,COUR,FUNKNO,SOURCE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Source term calculation for mesh centered finite differences in
* hexagonal geometry.
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
* IEX     iteration number.
* MAXKN   dimension of array KN.
* MAXQF   dimension of array QFR.
* NGCOND  number of groups condensed.
* NMERGE  number of merged regions.
* NALBP   number of physical albedos.
* FLXMER  flux estimate per mixture.
* SPH     SPH factors.
* SUNMER  incoming source (scattering+fission) cross sections.
* NREG    number of hexagons in BIVAC.
* NUN     number of unknowns per group in BIVAC.
* LL4     order of the system matrices in BIVAC. Equal to the number
*         of finite elements (hexagons or triangles) excluding the
*         virtual elements.
* ISPLH   hexagonal geometry flag:
*         =1: hexagonal elements; >1: triangular elements
* NBMIX   number of macro-mixtures.
* SIDE    side of the hexagons.
* MAT     mixture index per hexagon.
* KN      element-ordered unknown list.
* QFR     element-ordered information.
* VOL     volume of hexagons.
* DIFF0   diffusion coefficients.
* MERG    index of merged hexagons per macro-mixture.
* COUR    four times the incoming current per unit surface.
* FUNKNO  previously calculated fluxes.
*
*Parameters: output
* SOURCE  fission and diffusion sources.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IEX,MAXKN,MAXQF,NGCOND,NMERGE,NALBP,NREG,NUN,LL4,ISPLH,
     1 NBMIX,MAT(NREG),KN(MAXKN),MERG(NBMIX)
      REAL FLXMER(NMERGE,NGCOND),SPH(NMERGE+NALBP,NGCOND),
     1 SUNMER(NMERGE,NGCOND),SIDE,QFR(MAXQF),VOL(NREG),DIFF0(NBMIX),
     2 COUR,FUNKNO(NUN,NGCOND),SOURCE(NREG)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION A,DHARM
      DHARM(X1,X2,DIF1,DIF2)=2.0D0*DIF1*DIF2/(X1*DIF2+X2*DIF1)
*
      IF(ISPLH.EQ.1) THEN
         DS=SQRT(3.0)*SIDE
         FACT=2.0/(3.0*DS)
         NSURF=6
      ELSE
         DS=SIDE/(SQRT(3.0)*REAL(ISPLH-1))
         FACT=4.0/(3.0*DS)
         NSURF=3
      ENDIF
*
      CALL XDRSET(SOURCE,LL4,0.0)
*----
*  INCOMING CURRENT SOURCE.
*----
      IF(COUR.NE.0.0) THEN
         NUM1=0
         DO 50 IND1=1,LL4
         KHEX=KN(NUM1+NSURF+1)
         IF(VOL(KHEX).EQ.0.0) GO TO 40
         L=MAT(KHEX)
         VOL0=QFR(NUM1+NSURF+1)
         SIDEB=FACT*VOL0
         DO 30 IC=1,NSURF
         IF(KN(NUM1+IC).EQ.-1) THEN
            A=DHARM(DS,DS,DIFF0(L),0.5*DS*QFR(NUM1+IC))
            SOURCE(IND1)=SOURCE(IND1)+REAL(A)*COUR*SIDEB
         ENDIF
   30    CONTINUE
   40    NUM1=NUM1+NSURF+1
   50    CONTINUE
      ENDIF
*----
*  DIFFUSION AND FISSION SOURCE.
*----
      IF(IEX.EQ.1) GO TO 90
      DO 80 JGR=1,NGCOND
      NUM1=0
      DO 70 IND1=1,LL4
      KHEX=KN(NUM1+NSURF+1)
      IF(VOL(KHEX).EQ.0.0) GO TO 60
      L=MAT(KHEX)
      VOL0=QFR(NUM1+NSURF+1)
      GARS=SUNMER(MERG(L),JGR)*SPH(MERG(L),JGR)
      SOURCE(IND1)=SOURCE(IND1)+FUNKNO(IND1,JGR)*VOL0*GARS
   60 NUM1=NUM1+NSURF+1
   70 CONTINUE
   80 CONTINUE
      RETURN
*----
*  INITIAL ITERATION.
*----
   90 NUM1=0
      DO 120 IND1=1,LL4
      KHEX=KN(NUM1+NSURF+1)
      IF(VOL(KHEX).EQ.0.0) GO TO 110
      L=MAT(KHEX)
      VOL0=QFR(NUM1+NSURF+1)
      PV=0.0
      DO 100 JGR=1,NGCOND
      PV=PV+SUNMER(MERG(L),JGR)*FLXMER(MERG(L),JGR)
  100 CONTINUE
      SOURCE(IND1)=SOURCE(IND1)+VOL0*PV
  110 NUM1=NUM1+NSURF+1
  120 CONTINUE
      RETURN
      END
