*DECK SYBRXE
      SUBROUTINE SYBRXE(IPAS,NPIJ,NSUPCE,RAYRE,SIGT,SIGW,P,IQUAD,ISTAT,
     1 NMC,PROCEL,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the reduced collision probabilities for the
* 'do-it-yourself' approach.
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
*  IPAS   number of volumes.
*  NPIJ   length of cellwise scattering-reduced collision probability
*         matrices.
*  NSUPCE number of cells.
*  RAYRE  radius of the tubes in each cell.
*  SIGT   total macroscopic cross sections.
*  SIGW   P0 within-group scattering macroscopic cross sections.
*  IQUAD  quadrature parameter.
*  ISTAT  istat=1 for the statistical approximation.
*  NMC    offset of the first volume in each cell.
*  PROCEL user supplied geometrical matrix.
*  IMPX   print flag (equal to 0 for no print).
*
*Parameters: output
*  P      reduced collision probabilities.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPAS,NPIJ,NSUPCE,IQUAD,ISTAT,NMC(NSUPCE+1),IMPX
      REAL RAYRE(NSUPCE+IPAS),SIGT(IPAS),SIGW(IPAS),P(IPAS,IPAS),
     1 PROCEL(NSUPCE,NSUPCE)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (PI=3.141592654)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IAPPAR
      REAL, ALLOCATABLE, DIMENSION(:) :: PIJW,PISW,PSJW,PSSW
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PSSB
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IAPPAR(IPAS))
      ALLOCATE(PIJW(NPIJ),PISW(IPAS),PSJW(IPAS),PSSB(NSUPCE,2*NSUPCE),
     1 PSSW(NSUPCE))
*
      CALL SYB001 (NMC(NSUPCE+1),NSUPCE,NPIJ,SIGT,SIGW,IMPX,IQUAD,NMC,
     1 RAYRE,PIJW,PISW,PSJW,PSSW)
*
      IPIJ=0
      DO 30 IKK=1,NSUPCE
      J1=NMC(IKK)
      J2=NMC(IKK+1)-J1
      DO 20 I=1,J2
      IAPPAR(I+J1)=IKK
      DO 10 J=1,J2
      P(J1+I,J1+J)=PIJW(IPIJ+(J-1)*J2+I)
   10 CONTINUE
   20 CONTINUE
      IPIJ=IPIJ+J2*J2
   30 CONTINUE
*----
*  COMPUTATION OF QIJ MATRIX
*----
      IF(ISTAT.EQ.0) THEN
         DO 50 I=1,NSUPCE
         DO 40 J=1,NSUPCE
         PSSB(I,J)=-PROCEL(I,J)*PSSW(J)
         PSSB(I,NSUPCE+J)=PROCEL(I,J)
   40    CONTINUE
         PSSB(I,I)=1.0+PSSB(I,I)
   50    CONTINUE
         CALL ALSB(NSUPCE,NSUPCE,PSSB,IER,NSUPCE)
         IF(IER.NE.0) CALL XABORT('SYBRXE: SINGULAR MATRIX.')
*        COMPUTATION OF PIJ MATRIX.
         DO 70 I=1,IPAS
         K=IAPPAR(I)
         DO 60 J=1,IPAS
         L=IAPPAR(J)
         XX=PISW(I)*PSSB(K,NSUPCE+L)*PSJW(J)
         IF(L.EQ.K) THEN
            P(I,J)=P(I,J)+XX
         ELSE
            P(I,J)=XX
         ENDIF
   60    CONTINUE
   70    CONTINUE
      ELSE
         X1=0.0
         DO 80 I=1,NSUPCE
         X1=X1+PSSW(I)*PROCEL(1,I)
   80    CONTINUE
         X1=1.0/(1.0-X1)
         DO 100 J=1,IPAS
         L=IAPPAR(J)
         ZZZ=PROCEL(1,L)*X1*PSJW(J)
         DO 90 I=1,IPAS
         K=IAPPAR(I)
         XX=PISW(I)*ZZZ
         IF(L.EQ.K) THEN
            P(I,J)=P(I,J)+XX
         ELSE
            P(I,J)=XX
         ENDIF
   90    CONTINUE
  100    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(PSSW,PSSB,PSJW,PISW,PIJW)
      DEALLOCATE(IAPPAR)
      RETURN
      END
