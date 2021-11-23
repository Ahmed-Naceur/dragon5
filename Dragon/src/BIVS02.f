*DECK BIVS02
      SUBROUTINE BIVS02(CYLIND,IELEM,ICOL,NREG,NUN,LL4,IIMAX,MAT,KN,BFR,
     1 VOL,MU,SOURCE,SYS,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* One-speed flux calculation in mixed-dual finite element approximation
* (Cartesian geometry).
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
* CYLIND  cylinderization flag (=.TRUE. for cylindrical geometry).
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*          =2 (parabolic); =3 (cubic); =4 (quartic).
* ICOL    type of quadrature: =1 (analytical integration);
*         =2 (Gauss-Lobatto); =3 (Gauss-Legendre).
* NREG    number of elements in BIVAC.
* NUN     dimension of array FUNKNO.
* LL4     number of unknowns per group in BIVAC.
* IIMAX   allocated dimension of array SYS.
* MAT     mixture index per region.
* KN      element-ordered unknown list.
* BFR     element-ordered surface fractions.
* VOL     volume of regions.
* MU      indices used with compressed diagonal storage mode matrix SYS.
* SOURCE  fission and diffusion sources.
* SYS     factorized system matrix.
*
*Parameters: output
* FUNKNO  neutron fluxes (surface-averaged flux are in position
*         FUNKNO(NUN)).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IELEM,ICOL,NREG,NUN,LL4,IIMAX,MAT(NREG),KN(5*NREG),
     1 MU(LL4)
      REAL BFR(4*NREG),VOL(NREG),SOURCE(LL4),SYS(IIMAX),FUNKNO(NUN)
      LOGICAL CYLIND
*
      IF((CYLIND).AND.((IELEM.GT.1).OR.(ICOL.NE.2)))
     1 CALL XABORT('BIVS02: TYPE OF DISCRETIZATION NOT IMPLEMENTED.')
*----
*  RESOLUTION.
*----
      DO 90 I=1,LL4
      FUNKNO(I)=SOURCE(I)
   90 CONTINUE
      CALL ALLDLS (LL4,MU,SYS,FUNKNO)
*----
*  CALCULATION OF THE SURFACE-AVERAGED FLUX.
*----
      FUNKNO(NUN)=0.0
      NUM1=0
      NUM2=0
      DO 160 K=1,NREG
      IF(MAT(K).EQ.0) GO TO 160
      IF(VOL(K).EQ.0.0) GO TO 150
      IF(BFR(NUM2+1).GT.0.0) THEN
         SG=1.0
         DO 110 I0=1,IELEM
         IND1=KN(NUM1+1)+I0-1
         FUNKNO(NUN)=FUNKNO(NUN)+SG*SQRT(REAL(2*I0-1))*FUNKNO(IND1)*
     1   BFR(NUM2+1)
         SG=-SG
  110    CONTINUE
      ENDIF
      IF(BFR(NUM2+2).GT.0.0) THEN
         DO 120 I0=1,IELEM
         IND1=KN(NUM1+1)+I0-1
         FUNKNO(NUN)=FUNKNO(NUN)+SQRT(REAL(2*I0-1))*FUNKNO(IND1)*
     1   BFR(NUM2+2)
  120    CONTINUE
      ENDIF
      IF(BFR(NUM2+3).GT.0.0) THEN
         SG=1.0
         DO 130 I0=1,IELEM
         IND1=KN(NUM1+1)+(I0-1)*IELEM
         FUNKNO(NUN)=FUNKNO(NUN)+SG*SQRT(REAL(2*I0-1))*FUNKNO(IND1)*
     1   BFR(NUM2+3)
         SG=-SG
  130    CONTINUE
      ENDIF
      IF(BFR(NUM2+4).GT.0.0) THEN
         DO 140 I0=1,IELEM
         IND1=KN(NUM1+1)+(I0-1)*IELEM
         FUNKNO(NUN)=FUNKNO(NUN)+SQRT(REAL(2*I0-1))*FUNKNO(IND1)*
     1   BFR(NUM2+4)
  140    CONTINUE
      ENDIF
  150 NUM1=NUM1+5
      NUM2=NUM2+4
  160 CONTINUE
      RETURN
      END
