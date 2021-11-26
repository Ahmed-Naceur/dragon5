*DECK SYB005
      SUBROUTINE SYB005 (NGEN,NPIJ,NPIS,NRAYRE,SIGT2,SIGW2,IMPX,NCOUR,
     1 IQUAD,XX,YY,LSECT,NMC,NMCR,RAYRE,MAIL,IZMAIL,RZMAIL,PIJW,PISW,
     2 PSJW,PSSW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the cellwise scattering-reduced collision, escape and
* transmission probabilities in a 2-D Cartesian or hexagonal assembly
* with DP-1 approximation.
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
* NGEN    total number of generating cells.
* NPIJ    length of cellwise scattering-reduced collision probability
*         matrices.
* NPIS    length of cellwise scattering-reduced escape probability
*         matrices (NPIS=NMC(NGEN+1)).
* NRAYRE  size of array rayre (NRAYRE=NMCR(NGEN+1)).
* SIGT2   total macroscopic cross sections.
* SIGW2   P0 within-group scattering macroscopic cross sections.
* IMPX    print flag (equal to 0 for no print).
* NCOUR   number of currents surrounding the cells (=12 Cartesian
*         lattice; =18 hexagonal lattice).
* IQUAD   quadrature parameters.
* XX      X-thickness of the generating cells.
* YY      Y-thickness of the generating cells.
* LSECT   type of sectorization.
* NMC     offset of the first volume in each generating cell.
* NMCR    offset of the first radius in each generating cell.
* RAYRE   radius of the tubes in each generating cell.
* MAIL    offset of the first tracking information in each generating
*         cell.
* IZMAIL  integer tracking information.
* RZMAIL  real tracking information.
*
*Parameters: output
* PIJW    cellwise scattering-reduced collision probability matrices.
* PISW    cellwise scattering-reduced escape probability matrices.
* PSJW    cellwise scattering-reduced collision probability matrices
*         for incoming neutrons.
* PSSW    cellwise scattering-reduced transmission probability
*         matrices.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGEN,NPIJ,NPIS,NRAYRE,IMPX,NCOUR,IQUAD(4),LSECT(NGEN),
     1 NMC(NGEN+1),NMCR(NGEN+1),MAIL(2,NGEN),IZMAIL(*)
      REAL SIGT2(NPIS),SIGW2(NPIS),XX(NGEN),YY(NGEN),RAYRE(NRAYRE),
     1 RZMAIL(*),PIJW(NPIJ),PISW(NCOUR*NPIS),PSJW(NCOUR*NPIS),
     2 PSSW(NGEN*NCOUR*NCOUR)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (PI=3.141592654,SMALL=5.0E-3,SQRT3=1.732050807568877)
      LOGICAL LSKIP
      REAL PSS(324),SURFA(6),ALPA(64),PWA(64)
      DOUBLE PRECISION PPP(16)
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL,WORK
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PIS,PSJ,PP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PIS(NCOUR,NPIS),PSJ(NCOUR,NPIS))
*
      IPIJ=0
      IPIS=0
      IPSS=0
      DO 240 JKG=1,NGEN
      J1=NMC(JKG)
      J2=NMC(JKG+1)-J1
      J1R=NMCR(JKG)
      J2R=NMCR(JKG+1)-J1R
      ALLOCATE(PP(J2,J2),VOL(J2))
*----
*  COMPUTE THE REDUCED COLLISION PROBABILITY MATRIX
*----
      A=XX(JKG)
      B=YY(JKG)
      IF((NCOUR.EQ.12).AND.(LSECT(JKG).NE.0)) THEN
*        SECTORIZED CARTESIAN CELL.
         IB1=MAIL(1,JKG)
         IB2=MAIL(2,JKG)
         IF(LSECT(JKG).EQ.-999) THEN
            NSECT=4
         ELSE IF((LSECT(JKG).EQ.-1).OR.(LSECT(JKG).EQ.-101)) THEN
            NSECT=8
         ELSE
            NSECT=4*MOD(ABS(LSECT(JKG)),100)
         ENDIF
         MNA4=4*IQUAD(1)
         CALL SYB4QG(IMPX,3,MNA4,J2R,NSECT,LSECT(JKG),J2,RZMAIL(IB2),
     1   IZMAIL(IB1),A,B,RAYRE(J1R+2),SIGT2(J1+1),SMALL,VOL,PP,PIS,PSS)
      ELSE IF(LSECT(JKG).NE.0) THEN
*        SECTORIZED HEXAGONAL CELL.
         IB1=MAIL(1,JKG)
         IB2=MAIL(2,JKG)
         NSECT=6
         MNA4=12*IQUAD(1)
         CALL SYB7QG(IMPX,3,MNA4,J2R,NSECT,LSECT(JKG),J2,RZMAIL(IB2),
     1   IZMAIL(IB1),A,RAYRE(J1R+2),SIGT2(J1+1),SMALL,VOL,PP,PIS,PSS)
      ELSE IF((NCOUR.EQ.12).AND.(J2.EQ.1)) THEN
         CALL ALGPT(IQUAD(3),-1.0,1.0,ALPA,PWA)
         CALL RECT2(IQUAD(3),A,B,SIGT2(J1+1),SMALL,PP,PIS,PSS,ALPA,PWA)
         VOL(1)=A*B
      ELSE IF(J2.EQ.1) THEN
         CALL ALGPT(IQUAD(3),-1.0,1.0,ALPA,PWA)
         CALL XHX2D1(IQUAD(3),ALPA,PWA,A,SIGT2(J1+1),SMALL,PP,PIS,PSS,
     1   PPP)
         VOL(1)=1.5*SQRT3*A*A
      ELSE
*        NON-SECTORIZED CARTESIAN OR HEXAGONAL CELL.
         IB1=MAIL(1,JKG)
         IB2=MAIL(2,JKG)
         CALL SYBUP1(RZMAIL(IB2),IZMAIL(IB1),NCOUR/3,J2,SIGT2(J1+1),
     1   SMALL,A,B,IMPX,VOL,PP,PIS,PSS)
      ENDIF
      IF(IMPX.GE.8) CALL SYBPRX(NCOUR,J2,JKG,SIGT2(J1+1),PP,PIS,PSS)
*----
*  COMPUTE THE REDUCED COLLISION PROBABILITY MATRIX FOR INCOMING
*  NEUTRONS
*----
      DO 65 I=1,J2
      IF(NCOUR.EQ.12) THEN
         SURFA(1)=0.25*B
         SURFA(2)=0.25*B
         SURFA(3)=0.25*A
         SURFA(4)=0.25*A
      ELSE
         DO 50 JC=1,6
         SURFA(JC)=0.25*A
   50    CONTINUE
      ENDIF
      DO 60 JC=1,NCOUR
      PSJ(JC,I)=PIS(JC,I)*VOL(I)/SURFA(1+(JC-1)/3)
   60 CONTINUE
   65 CONTINUE
      DEALLOCATE(VOL)
*----
*  CHECK IF SCATTERING REDUCTION IS REQUIRED
*----
      LSKIP=.TRUE.
      DO 70 I=1,J2
      LSKIP=LSKIP.AND.(SIGW2(J1+I).EQ.0.0)
   70 CONTINUE
*----
*  SCATTERING REDUCTION IF LSKIP=.FALSE.
*----
      IF(LSKIP) THEN
*        DO NOT PERFORM SCATTERING REDUCTION.
         DO 85 I=1,J2
         DO 80 J=1,J2
         PIJW(IPIJ+(J-1)*J2+I)=PP(I,J)
   80    CONTINUE
   85    CONTINUE
         DO 95 I=1,J2
         DO 90 JC=1,NCOUR
         PISW(IPIS+(JC-1)*J2+I)=PIS(JC,I)
         PSJW(IPIS+(I-1)*NCOUR+JC)=PSJ(JC,I)
   90    CONTINUE
   95    CONTINUE
         DO 105 IC=1,NCOUR
         DO 100 JC=1,NCOUR
         PSSW(IPSS+(JC-1)*NCOUR+IC)=PSS((JC-1)*NCOUR+IC)
  100    CONTINUE
  105    CONTINUE
      ELSE
*        COMPUTE THE SCATTERING-REDUCED COLLISION AND ESCAPE MATRICES.
         DO 120 I=1,J2
         DO 110 J=1,J2
         PIJW(IPIJ+(J-1)*J2+I)=-PP(I,J)*SIGW2(J1+J)
  110    CONTINUE
         PIJW(IPIJ+(I-1)*J2+I)=1.0+PIJW(IPIJ+(I-1)*J2+I)
  120    CONTINUE
         CALL ALINV(J2,PIJW(IPIJ+1),J2,IER)
         IF(IER.NE.0) CALL XABORT('SYB005: SINGULAR MATRIX.')
         ALLOCATE(WORK(J2))
         DO 175 I=1,J2
         DO 130 K=1,J2
         WORK(K)=PIJW(IPIJ+(K-1)*J2+I)
  130    CONTINUE
         DO 150 J=1,J2
         WGAR=0.0
         DO 140 K=1,J2
         WGAR=WGAR+WORK(K)*PP(K,J)
  140    CONTINUE
         PIJW(IPIJ+(J-1)*J2+I)=WGAR
  150    CONTINUE
         DO 170 JC=1,NCOUR
         WGAR=0.0
         DO 160 K=1,J2
         WGAR=WGAR+WORK(K)*PIS(JC,K)
  160    CONTINUE
         PISW(IPIS+(JC-1)*J2+I)=WGAR
  170    CONTINUE
  175    CONTINUE
         DEALLOCATE(WORK)
*
*        COMPUTE THE SCATTERING-REDUCED COLLISION PROBABILITY MATRIX
*        FOR INCOMING NEUTRONS.
         DO 200 IC=1,NCOUR
         DO 190 J=1,J2
         WGAR=PSJ(IC,J)
         DO 180 K=1,J2
         WGAR=WGAR+PSJ(IC,K)*SIGW2(J1+K)*PIJW(IPIJ+(J-1)*J2+K)
  180    CONTINUE
         PSJW(IPIS+(J-1)*NCOUR+IC)=WGAR
  190    CONTINUE
  200    CONTINUE
*
*        COMPUTE THE SCATTERING-REDUCED TRANSMISSION PROBABILITY MATRIX.
         DO 230 IC=1,NCOUR
         DO 220 JC=1,NCOUR
         WGAR=PSS((JC-1)*NCOUR+IC)
         DO 210 K=1,J2
         WGAR=WGAR+PSJ(IC,K)*SIGW2(J1+K)*PISW(IPIS+(JC-1)*J2+K)
  210    CONTINUE
         PSSW(IPSS+(JC-1)*NCOUR+IC)=WGAR
  220    CONTINUE
  230    CONTINUE
      ENDIF
      DEALLOCATE(PP)
      IPIJ=IPIJ+J2*J2
      IPIS=IPIS+J2*NCOUR
      IPSS=IPSS+NCOUR*NCOUR
  240 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(PSJ,PIS)
      RETURN
      END
