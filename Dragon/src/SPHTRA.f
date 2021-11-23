*DECK SPHTRA
      SUBROUTINE SPHTRA(JPSYS,IEX,NPSYS,KSPH,NREG,NUN,NMERGE,NALBP,
     1 NGCOND,SUNMER,FLXMER,NBMIX,MAT,VOL,KEY,MERG,SPH,SIGW,SIGT,
     2 COURIN,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transport calculation over the macro-geometry using the collision
* probability technique. Use the Bell factor acceleration strategy.
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
* JPSYS   pointer to the 'GROUP' directory in the system LCM object.
* IEX     iteration number.
* NPSYS   group masks.
* KSPH    type of SPH factor normalization:
*         <0 asymptotic normalization;
*         =1 average flux normalization;
*         =2 Selengut normalization;
*         =3 generalized Selengut normalization (EDF type);
*         =4 Selengut normalization with surface leakage.
* NREG    number of macro-regions (in the macro calculation).
* NUN     number of unknowns per group in macro-calculation.
* NMERGE  number of merged regions.
* NALBP   number of physical albedos.
* NGCOND  number of condensed groups.
* SUNMER  incoming source (scattering+fission) cross sections.
* FLXMER  flux estimate per mixture.
* NBMIX   number of material mixtures.
* MAT     mixture index per macro-region.
* VOL     volume of macro-regions.
* KEY     position of the flux components associated with each volume.
* MERG    index of merged regions.
* SPH     SPH factors.
* SIGW    transport correction.
* SIGT    macroscopic total cross section.
* COURIN  averaged flux if KSPH=1. Equal to 4 times the incoming current
*         per unit surface if KSPH=2 or 3.
*
*Parameters: output
* FUNKNO  neutron flux.
*
*Reference(s):
* P. Blanc-Tranchant, A. Santamarina, G. Willermoz and A. Hebert,
* "Definition and Validation of a 2-D Transport Scheme for PWR Control
* Rod Clusters", paper presented at the Int. Conf. on Mathematics and
* Computation, Reactor Physics and Environmental Analysis in Nuclear
* Applications, Madrid, Spain, September 27-30, 1999.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) JPSYS
      INTEGER IEX,NPSYS(NGCOND),KSPH,NREG,NUN,NMERGE,NALBP,NGCOND,NBMIX,
     1 MAT(NREG),KEY(NREG),MERG(NBMIX)
      REAL SUNMER(NMERGE,NGCOND,NGCOND),FLXMER(NMERGE,NGCOND),VOL(NREG),
     1 SPH(NMERGE+NALBP,NGCOND),SIGW(NMERGE,NGCOND),SIGT(NMERGE,NGCOND),
     2 COURIN(NGCOND),FUNKNO(NUN,NGCOND)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) KPSYS
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGMA,SUNKNO
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PIJ
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: WORK2,WORK3
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SIGMA(0:NBMIX),SUNKNO(NREG),PIJ(NREG,NREG))
      ALLOCATE(WORK2(NREG+1,NREG+1),WORK3(NREG,NREG+1))
*----
*  GLOBAL SOURCE FOR THE BELL FACTOR METHOD.
*----
      DO 100 IGR=1,NGCOND
      IF(NPSYS(IGR).EQ.0) GO TO 100
      CALL XDRSET(SUNKNO,NREG,0.0)
      IF(IEX.EQ.1) THEN
        DO 20 IREG=1,NREG
        IMAT=MAT(IREG)
        IMERG=MERG(IMAT)
        IF(IMAT.EQ.0) GO TO 20
        IF(VOL(IREG).EQ.0.0) GO TO 20
        SUM=-(SIGT(IMERG,IGR)-SIGW(IMERG,IGR))*FLXMER(IMERG,IGR)
        DO 10 JGR=1,NGCOND
        SUM=SUM+SUNMER(IMERG,JGR,IGR)*FLXMER(IMERG,JGR)
   10   CONTINUE
        SUNKNO(IREG)=SUM
   20   CONTINUE
      ELSE
        DO 30 IREG=1,NREG
        IMAT=MAT(IREG)
        IMERG=MERG(IMAT)
        IF(IMAT.EQ.0) GO TO 30
        IF(VOL(IREG).EQ.0.0) GO TO 30
        GARS=-(SIGT(IMERG,IGR)-SIGW(IMERG,IGR))*SPH(IMERG,IGR)
        SUM=FUNKNO(KEY(IREG),IGR)*GARS
        DO 25 JGR=1,NGCOND
        GARS=SUNMER(IMERG,JGR,IGR)*SPH(IMERG,JGR)
        SUM=SUM+FUNKNO(KEY(IREG),JGR)*GARS
   25   CONTINUE
        SUNKNO(IREG)=SUM
   30   CONTINUE
      ENDIF
*----
*  COMPUTE THE WORK2 MATRIX.
*----
      KPSYS=LCMGIL(JPSYS,IGR)
      CALL LCMGET(KPSYS,'DRAGON-TXSC',SIGMA)
      CALL LCMGET(KPSYS,'DRAGON-PCSCT',PIJ)
      DO 45 I=1,NREG
      WORK2(I,NREG+1)=0.0D0
      DO 40 J=1,NREG
      WORK2(I,NREG+1)=WORK2(I,NREG+1)+PIJ(I,J)*VOL(I)*SUNKNO(J)
      WORK2(I,J)=PIJ(I,J)*VOL(I)
   40 CONTINUE
   45 CONTINUE
*----
*  COMPUTE THE NEUTRON FLUXES.
*----
      IF(KSPH.LT.0) THEN
*        ASYMPTOTIC NORMALIZATION.
         VOLTOT=0.0
         DO 60 I=1,NREG
         IF(MAT(I).EQ.-KSPH) THEN
            VOLTOT=VOLTOT+VOL(I)
            WORK2(NREG+1,I)=VOL(I)
         ELSE
            WORK2(NREG+1,I)=0.0D0
         ENDIF
         DO 50 J=1,NREG
         JBM=MAT(J)
         WORK2(I,J)=-SIGMA(JBM)*WORK2(I,J)
   50    CONTINUE
         WORK2(I,I)=WORK2(I,I)+VOL(I)
   60    CONTINUE
         WORK2(NREG+1,NREG+1)=COURIN(IGR)*VOLTOT
      ELSE
*        INTEGRATED FLUX OR SELENGUT NORMALIZATION.
         VOLTOT=0.0
         DO 80 I=1,NREG
         VOLTOT=VOLTOT+VOL(I)
         WORK2(NREG+1,I)=VOL(I)
         DO 70 J=1,NREG
         JBM=MAT(J)
         WORK2(I,J)=-SIGMA(JBM)*WORK2(I,J)
   70    CONTINUE
         WORK2(I,I)=WORK2(I,I)+VOL(I)
   80    CONTINUE
         WORK2(NREG+1,NREG+1)=COURIN(IGR)*VOLTOT
      ENDIF
      CALL ALSVDF(WORK2,NREG+1,NREG,NREG+1,NREG,WORK3(1,NREG+1),
     1 WORK3)
      CALL ALSVDS(WORK2,WORK3(1,NREG+1),WORK3,NREG+1,NREG,NREG+1,
     1 NREG,WORK2(1,NREG+1),WORK2(1,NREG+1))
      CALL XDRSET(FUNKNO(1,IGR),NUN,0.0)
      DO 90 I=1,NREG
      FUNKNO(KEY(I),IGR)=REAL(WORK2(I,NREG+1))
   90 CONTINUE
  100 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK3,WORK2)
      DEALLOCATE(PIJ,SUNKNO,SIGMA)
      RETURN
      END
