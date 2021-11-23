*DECK SYB001
      SUBROUTINE SYB001 (NREG,NSUPCE,NPIJ,SIGT,SIGW,IMPX,IQUAD,NMC,
     1 RAYRE,PIJW,PISW,PSJW,PSSW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the cellwise scattering-reduced collision, escape and
* transmission probabilities for the 'do-it-yourself' approach.
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
* NREG    total number of regions (NREG=NMC(NSUPCE+1)).
* NSUPCE  total number of cells.
* NPIJ    length of cellwise scattering-reduced collision probability
*         matrices.
* SIGT    total macroscopic cross sections.
* SIGW    P0 within-group scattering macroscopic cross sections.
* IMPX    print flag (equal to 0 for no print).
* IQUAD   quadrature parameter.
* NMC     offset of the first volume in each cell.
* RAYRE   radius of the tubes in each cell.
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
      INTEGER NSUPCE,NREG,NPIJ,IMPX,IQUAD,NMC(NSUPCE+1)
      REAL SIGT(NREG),SIGW(NREG),RAYRE(NREG+NSUPCE),PIJW(NPIJ),
     1 PISW(NREG),PSJW(NREG),PSSW(NSUPCE)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (PI=3.141592654)
      LOGICAL LSKIP
      REAL, ALLOCATABLE, DIMENSION(:) :: PIS,PSJ,ZTR,WORK
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PIS(NREG),PSJ(NREG))
*
      IPIJ=0
      DO 160 IKK=1,NSUPCE
      J1=NMC(IKK)
      J2=NMC(IKK+1)-J1
*----
*  COMPUTE THE REDUCED COLLISION PROBABILITY MATRIX
*----
      ALLOCATE(PP(J2,J2),ZTR(1+IQUAD*((J2*(5+J2))/2)))
      CALL SYBT1D(J2,RAYRE(J1+IKK),.FALSE.,IQUAD,ZTR)
      CALL SYBALC(J2,J2,RAYRE(J1+IKK),SIGT(J1+1),IQUAD,0.0,ZTR,PP)
      DEALLOCATE(ZTR)
      SURFA=2.0*PI*RAYRE(J1+J2+IKK)
      PSS=0.0
      RJN=0.0
      DO 20 I=1,J2
      PIS(I)=0.0
      DO 10 J=1,J2
      PIS(I)=PIS(I)+PP(I,J)*SIGT(J+J1)
   10 CONTINUE
      PIS(I)=1.0-PIS(I)
      RJN1=RAYRE(I+J1+IKK)**2
      PSJ(I)=4.0*PI*(RJN1-RJN)*PIS(I)/SURFA
      PSS=PSS+PSJ(I)*SIGT(I+J1)
      RJN=RJN1
   20 CONTINUE
      PSS=1.0-PSS
      IF(IMPX.GE.8) CALL SYBPRX (1,J2,IKK,SIGT(J1+1),PP,PIS,PSS)
*----
*  CHECK IF SCATTERING REDUCTION IS REQUIRED
*----
      LSKIP=.TRUE.
      DO 30 I=1,J2
      LSKIP=LSKIP.AND.(SIGW(J1+I).EQ.0.0)
   30 CONTINUE
*----
*  SCATTERING REDUCTION IF LSKIP=.FALSE.
*----
      IF(LSKIP) THEN
*        DO NOT PERFORM SCATTERING REDUCTION.
         DO 45 I=1,J2
         DO 40 J=1,J2
         PIJW(IPIJ+(J-1)*J2+I)=PP(I,J)
   40    CONTINUE
   45    CONTINUE
         DO 50 I=1,J2
         PISW(J1+I)=PIS(I)
         PSJW(J1+I)=PSJ(I)
   50    CONTINUE
         PSSW(IKK)=PSS
      ELSE
*        COMPUTE THE SCATTERING-REDUCED COLLISION AND ESCAPE MATRICES.
         DO 70 I=1,J2
         DO 60 J=1,J2
         PIJW(IPIJ+(J-1)*J2+I)=-PP(I,J)*SIGW(J1+J)
   60    CONTINUE
         PIJW(IPIJ+(I-1)*J2+I)=1.0+PIJW(IPIJ+(I-1)*J2+I)
   70    CONTINUE
         CALL ALINV(J2,PIJW(IPIJ+1),J2,IER)
         IF(IER.NE.0) CALL XABORT('SYB001: SINGULAR MATRIX.')
         ALLOCATE(WORK(J2))
         DO 120 I=1,J2
         DO 80 K=1,J2
         WORK(K)=PIJW(IPIJ+(K-1)*J2+I)
   80    CONTINUE
         DO 100 J=1,J2
         WGAR=0.0
         DO 90 K=1,J2
         WGAR=WGAR+WORK(K)*PP(K,J)
   90    CONTINUE
         PIJW(IPIJ+(J-1)*J2+I)=WGAR
  100    CONTINUE
         WGAR=0.0
         DO 110 K=1,J2
         WGAR=WGAR+WORK(K)*PIS(K)
  110    CONTINUE
         PISW(J1+I)=WGAR
  120    CONTINUE
         DEALLOCATE(WORK)
*
*        COMPUTE THE SCATTERING-REDUCED COLLISION PROBABILITY MATRIX
*        FOR INCOMING NEUTRONS.
         DO 140 J=1,J2
         WGAR=PSJ(J)
         DO 130 K=1,J2
         WGAR=WGAR+PSJ(K)*SIGW(J1+K)*PIJW(IPIJ+(J-1)*J2+K)
  130    CONTINUE
         PSJW(J1+J)=WGAR
  140    CONTINUE
*
*        COMPUTE THE SCATTERING-REDUCED TRANSMISSION PROBABILITY MATRIX.
         WGAR=PSS
         DO 150 K=1,J2
         WGAR=WGAR+PSJ(K)*SIGW(J1+K)*PISW(J1+K)
  150    CONTINUE
         PSSW(IKK)=WGAR
      ENDIF
      DEALLOCATE(PP)
      IPIJ=IPIJ+J2*J2
  160 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(PSJ,PIS)
      RETURN
      END
