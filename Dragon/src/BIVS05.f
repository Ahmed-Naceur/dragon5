*DECK BIVS05
      SUBROUTINE BIVS05(IELEM,NBLOS,NUN,LL4,IIMAX,IPERT,KN,BFR,MU,
     1 SOURCE,SYS,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* One-speed flux calculation in a Thomas-Raviart-Schneider (dual) finite
* element diffusion approximation (hexagonal geometry).
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*          =2 (parabolic); =3 (cubic).
* NBLOS   number of lozenges per direction, taking into account
*         mesh-splitting.
* NUN     dimension of array FUNKNO.
* LL4     number of unknowns per group in BIVAC.
* IIMAX   allocated dimension of array SYS.
* IPERT   mixture permutation index.
* KN      element-ordered unknown list.
* BFR     element-ordered surface fractions.
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
      INTEGER IELEM,NBLOS,NUN,LL4,IIMAX,IPERT(NBLOS),
     1 KN(NBLOS,4+6*IELEM*(IELEM+1)),MU(LL4)
      REAL BFR(NBLOS,6),SOURCE(LL4),SYS(IIMAX),FUNKNO(NUN)
*----
*  RESOLUTION.
*----
      DO 10 I=1,LL4
      FUNKNO(I)=SOURCE(I)
   10 CONTINUE
      CALL ALLDLS (LL4,MU,SYS,FUNKNO)
*----
*  CALCULATION OF THE SURFACE-AVERAGED FLUX.
*----
      FUNKNO(NUN)=0.0
      NUM=0
      DO 100 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 100
      NUM=NUM+1
      IF(BFR(NUM,1).GT.0.0) THEN
         SG=1.0
         DO 40 I0=1,IELEM
         IND1=KN(NUM,1)+I0-1
         FUNKNO(NUN)=FUNKNO(NUN)+SG*SQRT(REAL(2*I0-1))*FUNKNO(IND1)*
     1   BFR(NUM,1)
         SG=-SG
   40    CONTINUE
      ENDIF
      IF(BFR(NUM,2).GT.0.0) THEN
         DO 50 I0=1,IELEM
         IND1=KN(NUM,2)+(I0-1)*IELEM
         FUNKNO(NUN)=FUNKNO(NUN)+SQRT(REAL(2*I0-1))*FUNKNO(IND1)*
     1   BFR(NUM,2)
   50    CONTINUE
      ENDIF
      IF(BFR(NUM,3).GT.0.0) THEN
         SG=1.0
         DO 60 I0=1,IELEM
         IND1=KN(NUM,2)+I0-1
         FUNKNO(NUN)=FUNKNO(NUN)+SG*SQRT(REAL(2*I0-1))*FUNKNO(IND1)*
     1   BFR(NUM,3)
         SG=-SG
   60    CONTINUE
      ENDIF
      IF(BFR(NUM,4).GT.0.0) THEN
         DO 70 I0=1,IELEM
         IND1=KN(NUM,3)+(I0-1)*IELEM
         FUNKNO(NUN)=FUNKNO(NUN)+SQRT(REAL(2*I0-1))*FUNKNO(IND1)*
     1   BFR(NUM,4)
   70    CONTINUE
      ENDIF
      IF(BFR(NUM,5).GT.0.0) THEN
         SG=1.0
         DO 80 I0=1,IELEM
         IND1=KN(NUM,3)+I0-1
         FUNKNO(NUN)=FUNKNO(NUN)+SG*SQRT(REAL(2*I0-1))*FUNKNO(IND1)*
     1   BFR(NUM,5)
         SG=-SG
   80    CONTINUE
      ENDIF
      IF(BFR(NUM,6).GT.0.0) THEN
         DO 90 I0=1,IELEM
         IND1=KN(NUM,4)+(I0-1)*IELEM
         FUNKNO(NUN)=FUNKNO(NUN)+SQRT(REAL(2*I0-1))*FUNKNO(IND1)*
     1   BFR(NUM,6)
   90    CONTINUE
      ENDIF
  100 CONTINUE
      RETURN
      END
