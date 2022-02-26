*DECK SYBRX2
      SUBROUTINE SYBRX2 (IPAS,NPIJ,NPIS,SIGT,SIGW,P,IMPX,NCOUR,IWIGN,
     1 NMCEL,NMERGE,NGEN,IQUAD,XX,YY,NMC,RAYRE,MAIL,RZMAIL,IFR,ALB,
     2 INUM,IGEN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the global scattering-reduced collision probabilities in a
* 2-D Cartesian or hexagonal assembly using the interface current
* method with Roth approximation.
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
* IPAS    total number of volumes.
* NPIJ    length of cellwise scattering-reduced collision probability
*         matrices.
* NPIS    length of cellwise scattering-reduced collision probability
*         matrices (NPIS=NMC(NGEN+1)).
* SIGT    total macroscopic cross sections.
* SIGW    within group scattering cross sections.
* IMPX    print flag (equal to 0 for no print).
* NCOUR   number of currents surrounding the cells (=4 Cartesian
*         lattice; =6 hexagonal lattice).
* IWIGN   type of cylinderization.
* IQUAD   quadrature parameters.
* NMCEL   total number of cells in the domain.
* IFR     index-number of in-currents.
* ALB     transmission/albedo associated with each in-current.
* NMERGE  total number of merged cells for which specific values
*         of the neutron flux and reactions rates are required.
*         Many cells with different position in the domain can
*         be merged before the neutron flux calculation if they
*         own the same generating cell. This allows some reduction
*         in cpu time and memory (NMERGE.le.NMCEL).
* INUM    index-number of the merged cell associated to each cell.
*         Note: IFR and ALB contains information to rebuild the
*         geometrical 'A' matrix.
* NGEN    total number of generating cells. A generating cell is
*         defined by its material and its position in the domain
*         (NGEN.le.NMERGE).
* XX      X-thickness of the generating cells.
* YY      Y-thickness of the generating cells.
* NMC     offset of the first volume in each generating cell.
* RAYRE   radius of the tubes in each generating cell.
* MAIL    offset of the first tracking information in each
*         generating cell.
* RZMAIL  real tracking information.
* IGEN    index-number of the generating cell associated with each
*         merged cell.
*
*Parameters: output
* P       reduced collision probabilities.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPAS,NPIJ,NPIS,IMPX,NCOUR,IWIGN,NMCEL,NMERGE,NGEN,
     1 IQUAD(4),NMC(NGEN+1),MAIL(2,NGEN),IFR(NCOUR*NMCEL),INUM(NMCEL),
     2 IGEN(NMERGE)
      REAL SIGT(IPAS),SIGW(IPAS),P(IPAS,IPAS),XX(NGEN),YY(NGEN),
     1 RAYRE(NPIS),RZMAIL(*),ALB(NCOUR*NMCEL)
*----
*  LOCAL VARIABLES
*----
      REAL PIBB(6)
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGT2,SIGW2,PIJW,PISW,PSJW,PSSW
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PSSB
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PSSB(NMERGE,2*NMERGE),SIGT2(IPAS),SIGW2(IPAS),PIJW(NPIJ),
     1 PISW(NPIS),PSJW(NPIS),PSSW(NGEN))
*
      DO 20 I=1,IPAS
      DO 10 J=1,IPAS
      P(I,J)=0.0
   10 CONTINUE
   20 CONTINUE
      I1=0
      DO 40 IKK=1,NMERGE
      IKG=IGEN(IKK)
      J1=NMC(IKG)
      I2=NMC(IKG+1)-J1
      DO 30 I=1,I2
      SIGT2(J1+I)=SIGT(I1+I)
      SIGW2(J1+I)=SIGW(I1+I)
   30 CONTINUE
      I1=I1+I2
   40 CONTINUE
*
      CALL SYB002 (NGEN,NPIJ,NPIS,SIGT2,SIGW2,IMPX,NCOUR,IWIGN,IQUAD,
     1 XX,YY,NMC,RAYRE,MAIL,RZMAIL,PIJW,PISW,PSJW,PSSW)
*
      IPIJ=0
      DO 80 JKG=1,NGEN
      J2=NMC(JKG+1)-NMC(JKG)
      I1=0
      DO 70 IKK=1,NMERGE
      IKG=IGEN(IKK)
      I2=NMC(IKG+1)-NMC(IKG)
      IF(IKG.EQ.JKG) THEN
         DO 60 J=1,J2
         DO 50 I=1,J2
         P(I1+I,I1+J)=PIJW(IPIJ+(J-1)*J2+I)
   50    CONTINUE
   60    CONTINUE
      ENDIF
      I1=I1+I2
   70 CONTINUE
      IPIJ=IPIJ+J2*J2
   80 CONTINUE
*----
*  COMPUTATION OF PSSB=A*(I-PSS*A)**-1
*----
      DO 100 I=1,NMERGE
      DO 90 J=1,NMERGE
      PSSB(I,J)=0.0
      PSSB(I,NMERGE+J)=0.0
   90 CONTINUE
      PSSB(I,I)=1.0
  100 CONTINUE
      DO 130 ICEL=1,NMCEL
      IKK=INUM(ICEL)
      IKG=IGEN(IKK)
      IS=NCOUR*(ICEL-1)
      IF(NCOUR.EQ.4) THEN
         A=XX(IKG)
         B=YY(IKG)
         DEN1=2.0*(A+B)
         PIBB(1)=B/DEN1
         PIBB(2)=B/DEN1
         PIBB(3)=A/DEN1
         PIBB(4)=A/DEN1
      ELSE
         DO 110 JC=1,NCOUR
         PIBB(JC)=1.0/6.0
  110    CONTINUE
      ENDIF
      ZZZ=PSSW(IKG)
      DO 120 JC=1,NCOUR
      J1=IFR(IS+JC)
      ALBEDO=PIBB(JC)*ALB(IS+JC)
      PSSB(J1,NMERGE+IKK)=PSSB(J1,NMERGE+IKK)+ALBEDO
      PSSB(J1,IKK)=PSSB(J1,IKK)-ZZZ*ALBEDO
  120 CONTINUE
  130 CONTINUE
      CALL ALSB(NMERGE,NMERGE,PSSB,IER,NMERGE)
      IF(IER.NE.0) CALL XABORT('SYBRX2: SINGULAR MATRIX.')
*----
*  COMPUTATION OF PIS*PSSB*PSJ
*----
      I1=0
      DO 170 IKK=1,NMERGE
      IKG=IGEN(IKK)
      I1P=NMC(IKG)
      I2=NMC(IKG+1)-I1P
      DO 160 I=1,I2
      ZZZ=PISW(I1P+I)
      J1=0
      DO 150 JKK=1,NMERGE
      JKG=IGEN(JKK)
      J1P=NMC(JKG)
      J2=NMC(JKG+1)-J1P
      DO 140 J=1,J2
      P(I1+I,J1+J)=P(I1+I,J1+J)+ZZZ*PSSB(JKK,NMERGE+IKK)*PSJW(J1P+J)
  140 CONTINUE
      J1=J1+J2
  150 CONTINUE
  160 CONTINUE
      I1=I1+I2
  170 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(PSSW,PSJW,PISW,PIJW,SIGW2,SIGT2,PSSB)
      RETURN
      END
