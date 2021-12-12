*DECK BIVS04
      SUBROUTINE BIVS04(MAXKN,MAXQF,NREG,NUN,LL4,ISPLH,IIMAX,SIDE,KN,
     1 QFR,BFR,VOL,IDL,MU,SOURCE,SYS,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* One-speed flux calculation in mixed-dual finite element approximation
* (hexagonal geometry).
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
* MAXKN   dimension of array KN.
* MAXQF   dimension of arrays QFR and BFR.
* NREG    number of hexagons in BIVAC.
* NUN     dimension of array FUNKNO.
* LL4     number of unknowns per group in BIVAC. Equal to the number
*         of finite elements (hexagons or triangles) excluding the
*         virtual elements.
* ISPLH   type of hexagonal mesh-splitting:
*         =1: hexagonal elements; >1: triangular elements.
* IIMAX   allocated dimension of array SYS.
* SIDE    side of the hexagons.
* KN      element-ordered unknown list.
* QFR     element-ordered information.
* BFR     element-ordered surface fractions.
* VOL     volume of hexagons.
* IDL     position of the average flux component associated with
*         each volume.
* MU      indices used with the compressed diagonal storage mode matrix
*         SYS.
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
      INTEGER MAXKN,MAXQF,NREG,NUN,LL4,ISPLH,IIMAX,KN(MAXKN),IDL(NREG),
     1 MU(LL4)
      REAL SIDE,QFR(MAXQF),BFR(MAXQF),VOL(NREG),SOURCE(LL4),SYS(IIMAX),
     1 FUNKNO(NUN)
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
*----
*  RESOLUTION.
*----
      DO 10 I=1,LL4
      FUNKNO(I)=SOURCE(I)
   10 CONTINUE
      CALL ALLDLS (LL4,MU,SYS,FUNKNO)
*----
*  CALCULATION OF ELEMENT-AVERAGED FLUXES.
*----
      IF(ISPLH.GT.1) THEN
         NSURF=3
         DO 20 K=1,NREG
         IF(IDL(K).NE.0) FUNKNO(IDL(K))=0.0
   20    CONTINUE
         NUM1=0
         DO 40 IND1=1,LL4
         K=KN(NUM1+NSURF+1)
         IF(VOL(K).EQ.0.0) GO TO 30
         FUNKNO(IDL(K))=FUNKNO(IDL(K))+QFR(NUM1+NSURF+1)*FUNKNO(IND1)/
     1   VOL(K)
   30    NUM1=NUM1+NSURF+1
   40    CONTINUE
      ELSE
         NSURF=6
      ENDIF
*----
*  CALCULATION OF SURFACE-AVERAGED FLUX.
*----
      FUNKNO(NUN)=0.0
      NUM1=0
      DO 70 IND1=1,LL4
      K=KN(NUM1+NSURF+1)
      IF(VOL(K).EQ.0.0) GO TO 60
      DO 50 IC=1,NSURF
      IF(BFR(NUM1+IC).GT.0.0) THEN
         FUNKNO(NUN)=FUNKNO(NUN)+FUNKNO(IND1)*BFR(NUM1+IC)
      ENDIF
   50 CONTINUE
   60 NUM1=NUM1+NSURF+1
   70 CONTINUE
      RETURN
      END
