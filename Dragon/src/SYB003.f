*DECK SYB003
      SUBROUTINE SYB003 (NGEN,NPIJ,NPIS,SIGT2,SIGW2,IMPX,NCOUR,IWIGN,
     1 IQUAD,XX,YY,NMC,RAYRE,MAIL,RZMAIL,PIJW,PISW,PSJW,PSSW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the cellwise scattering-reduced collision, escape and
* transmission probabilities in a 2-D Cartesian or hexagonal assembly
* with Roth x 4 or Roth x 6 approximation.
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
* SIGT2   total macroscopic cross sections.
* SIGW2   P0 within-group scattering macroscopic cross sections.
* IMPX    print flag (equal to 0 for no print).
* NCOUR   number of currents surrounding the cells (=4 Cartesian
*         lattice; =6 hexagonal lattice).
* IWIGN   type of cylinderization.
* IQUAD   quadrature parameters.
* XX      X-thickness of the generating cells.
* YY      Y-thickness of the generating cells.
* NMC     offset of the first volume in each generating cell.
* RAYRE   radius of the tubes in each generating cell.
* MAIL    offset of the first tracking information in each generating
*         cell.
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
      INTEGER NGEN,NPIJ,NPIS,IMPX,NCOUR,IWIGN,IQUAD(4),NMC(NGEN+1),
     1 MAIL(2,NGEN)
      REAL SIGT2(NPIS),SIGW2(NPIS),XX(NGEN),YY(NGEN),RAYRE(NPIS),
     1 RZMAIL(*),PIJW(NPIJ),PISW(NCOUR*NPIS),PSJW(NCOUR*NPIS),
     2 PSSW(NGEN*NCOUR*NCOUR)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (PI=3.141592654)
      LOGICAL LSKIP
      INTEGER ISLR(4,4),ISLH(6,6)
      REAL PBB(6),PSS(36)
      REAL, ALLOCATABLE, DIMENSION(:) :: RAYR2,WORK
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PIS,PSJ,PP
*----
*  DATA STATEMENTS
*----
      SAVE ISLR,ISLH
      DATA ISLR/0,4,1,1,4,0,1,1,3,3,0,2,3,3,2,0/
      DATA ISLH/0,1,2,3,2,1,1,0,1,2,3,2,2,1,0,1,2,3,
     1          3,2,1,0,1,2,2,3,2,1,0,1,1,2,3,2,1,0/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PIS(NCOUR,NPIS),PSJ(NCOUR,NPIS))
*
      MR=IQUAD(4)
      IPIJ=0
      IPIS=0
      IPSS=0
      DO 240 JKG=1,NGEN
      J1=NMC(JKG)
      J2=NMC(JKG+1)-J1
*----
*  CYLINDERIZATION OPTIONS 
*----
      A=XX(JKG)
      B=YY(JKG)
      IB=MAIL(2,JKG)
      RJ1=RAYRE(NMC(JKG+1))
      SCALE1=1.0
      SCALE2=1.0
      ROUT=0.0
      IF((NCOUR.EQ.4).AND.(IWIGN.EQ.1)) THEN
*        ASKEW CYLINDERIZATION CARTESIAN.
         RJ2=(A+B)/PI
         SCALE1=(A*B-PI*RJ1**2)/(PI*RJ2**2-PI*RJ1**2)
         ROUT=RJ2
      ELSE IF((NCOUR.EQ.4).AND.(IWIGN.EQ.2)) THEN
*        WIGNER CYLINDERIZATION CARTESIAN.
         ROUT=SQRT(A*B/PI)
      ELSE IF((NCOUR.EQ.4).AND.(IWIGN.EQ.3)) THEN
*        SANCHEZ CYLINDERIZATION CARTESIAN.
         SCALE2=SQRT(PI*A*B)/(A+B)
         ROUT=SQRT(A*B/PI)
      ELSE IF(IWIGN.EQ.1) THEN
*        ASKEW CYLINDERIZATION HEXAGONAL.
         RJ2=3.0*A/PI
         SCALE1=(1.5*SQRT(3.0)*A*A-PI*RJ1**2)/(PI*RJ2**2-PI*RJ1**2)
         ROUT=RJ2
      ELSE IF(IWIGN.EQ.2) THEN
*        WIGNER CYLINDERIZATION HEXAGONAL.
         ROUT=SQRT(1.5*SQRT(3.0)/PI)*A
      ELSE IF(IWIGN.EQ.3) THEN
*        SANCHEZ CYLINDERIZATION HEXAGONAL.
         SCALE2=SQRT(PI*SQRT(3.0)/6.0)
         ROUT=SQRT(1.5*SQRT(3.0)/PI)*A
      ENDIF
      IF(ROUT.LE.RJ1) CALL XABORT('SYB003: CYLINDERIZATION ERROR.')
*----
*  COMPUTE THE REDUCED COLLISION PROBABILITY MATRIX
*----
      SURFA=0.5*PI*ROUT
      ALLOCATE(PP(J2,J2),RAYR2(J2+1))
      DO 10 I=1,J2
      RAYR2(I)=RAYRE(J1+I)
   10 CONTINUE
      RAYR2(J2+1)=ROUT
      SIGT2(J1+J2)=SIGT2(J1+J2)*SCALE1
      CALL SYBALC(J2,J2,RAYR2,SIGT2(J1+1),MR,0.0,RZMAIL(IB),PP)
*
      PSSX=0.0
      RJ1=0.0
      DO 50 I=1,J2
      PISX=1.0
      RJ2=RAYR2(I+1)**2
      VV=PI*(RJ2-RJ1)
      DO 20 J=1,J2
      PISX=PISX-PP(I,J)*SIGT2(J1+J)
   20 CONTINUE
      PSSX=PSSX+PISX*SIGT2(J1+I)*VV/SURFA
      IF(NCOUR.EQ.4) THEN
         DEN1=2.0*(A+B)
         PIS(1,I)=PISX*B/DEN1
         PIS(2,I)=PISX*B/DEN1
         PIS(3,I)=PISX*A/DEN1
         PIS(4,I)=PISX*A/DEN1
      ELSE
         DO 30 IC=1,NCOUR
         PIS(IC,I)=PISX/6.0
   30    CONTINUE
      ENDIF
*----
*  COMPUTE THE REDUCED COLLISION PROBABILITY MATRIX FOR INCOMING
*  NEUTRONS
*----
      SURFA2=(0.5*PI*ROUT)/SCALE2
      DO 40 IC=1,NCOUR
      IF(I.LE.J2-1) THEN
         PSJ(IC,I)=PISX*VV/SURFA2
      ELSE IF(I.EQ.J2) THEN
         PSJ(IC,I)=PISX*VV*SCALE1/SURFA2
      ENDIF
   40 CONTINUE
      RJ1=RJ2
   50 CONTINUE
      DEALLOCATE(RAYR2)
*----
*  COMPUTE THE TRANSMISSION PROBABILITIES
*----
      PSSX=1.0-SCALE2*PSSX
      IF(NCOUR.EQ.4) THEN
         A=XX(JKG)
         B=YY(JKG)
         DEN1=MAX(2.0*B+A,2.0*A+B)
         PBB(1)=B/DEN1
         PBB(2)=MAX(A,2.0*A-B)/DEN1
         PBB(3)=A/DEN1
         PBB(4)=MAX(B,2.0*B-A)/DEN1
      ELSE
         PBB(1)=1.0/5.0
         PBB(2)=1.0/5.0
         PBB(3)=1.0/5.0
      ENDIF
      DO 65 JC=1,NCOUR
      DO 60 IC=1,NCOUR
      IF(NCOUR.EQ.4) THEN
         IB=ISLR(IC,JC)
      ELSE
         IB=ISLH(IC,JC)
      ENDIF
      IF(IB.GT.0) THEN
         PSS((JC-1)*NCOUR+IC)=PSSX*PBB(IB)
      ELSE
         PSS((JC-1)*NCOUR+IC)=0.0
      ENDIF
   60 CONTINUE
   65 CONTINUE
      IF(IMPX.GE.8) CALL SYBPRX(NCOUR,J2,JKG,SIGT2(J1+1),PP,PIS,PSS)
      SIGT2(J1+J2)=SIGT2(J1+J2)/SCALE1
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
         DO 80 J=1,J2-1
         PIJW(IPIJ+(J-1)*J2+I)=PP(I,J)
   80    CONTINUE
         PIJW(IPIJ+(J2-1)*J2+I)=PP(I,J2)*SCALE1
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
         DO 110 J=1,J2-1
         PIJW(IPIJ+(J-1)*J2+I)=-PP(I,J)*SIGW2(J1+J)
  110    CONTINUE
         PIJW(IPIJ+(J2-1)*J2+I)=-PP(I,J2)*SIGW2(J1+J2)*SCALE1
         PIJW(IPIJ+(I-1)*J2+I)=1.0+PIJW(IPIJ+(I-1)*J2+I)
  120    CONTINUE
         CALL ALINV(J2,PIJW(IPIJ+1),J2,IER)
         IF(IER.NE.0) CALL XABORT('SYB003: SINGULAR MATRIX.')
         ALLOCATE(WORK(J2))
         DO 175 I=1,J2
         DO 130 K=1,J2
         WORK(K)=PIJW(IPIJ+(K-1)*J2+I)
  130    CONTINUE
         DO 150 J=1,J2-1
         WGAR=0.0
         DO 140 K=1,J2
         WGAR=WGAR+WORK(K)*PP(K,J)
  140    CONTINUE
         PIJW(IPIJ+(J-1)*J2+I)=WGAR
  150    CONTINUE
         WGAR=0.0
         DO 155 K=1,J2
         WGAR=WGAR+WORK(K)*PP(K,J2)
  155    CONTINUE
         PIJW(IPIJ+(J2-1)*J2+I)=WGAR*SCALE1
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
