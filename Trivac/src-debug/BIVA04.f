*DECK BIVA04
      SUBROUTINE BIVA04(ITY,MAXKN,MAXQF,SGD,NREG,LL4,ISPLH,NBMIX,IIMAX,
     1 SIDE,MAT,KN,QFR,VOL,MU,SYS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of a within-group (leakage and removal) or out-of-group
* system matrix in mesh-centered finite-difference diffusion
* approximation (hexagonal geometry).
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
* ITY     type of assembly: =0: leakage-removal matrix assembly;
*         =1: cross section matrix assembly.
* MAXKN   dimension of array KN.
* MAXQF   dimension of array QFR.
* SGD     nuclear properties. SGD(:,1) and SGD(:,2) are diffusion
*         coefficients. SGD(:,3) are removal macroscopic cross sections.
* NREG    number of hexagons in BIVAC.
* LL4     number of unknowns per group in BIVAC. Equal to the number
*         of finite elements (hexagons or triangles) excluding the
*         virtual elements.
* ISPLH   type of hexagonal mesh-splitting:
*         =1: hexagonal elements; >1: triangular elements.
* NBMIX   number of macro-mixtures.
* IIMAX   allocated dimension of array SYS.
* SIDE    side of the hexagons.
* MAT     mixture index per hexagon.
* KN      element-ordered unknown list.
* QFR     element-ordered information.
* VOL     volume of hexagons.
* MU      indices used with the compressed diagonal storage mode matrix
*         SYS.
*
*Parameters: output
* SYS     system matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ITY,MAXKN,MAXQF,NREG,LL4,ISPLH,NBMIX,IIMAX,MAT(NREG),
     1 KN(MAXKN),MU(LL4)
      REAL SGD(NBMIX,3),SIDE,QFR(MAXQF),VOL(NREG),SYS(IIMAX)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION A1,DHARM,VAR1
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
*----
*  ASSEMBLY OF A SYSTEM MATRIX.
*----
      IF(ITY.EQ.0) THEN
*        LEAKAGE-REMOVAL SYSTEM MATRIX ASSEMBLY.
         NUM1=0
         DO 35 IND1=1,LL4
         KHEX=KN(NUM1+NSURF+1)
         IF(VOL(KHEX).EQ.0.0) GO TO 30
         L=MAT(KHEX)
         VOL0=QFR(NUM1+NSURF+1)
         SIDEB=FACT*VOL0
         VAR1=0.0D0
         KEY0=MU(IND1)-IND1
         DO 20 IX=1,NSURF
         IND2=KN(NUM1+IX)
         A1=0.0
         IF(IND2.GT.0) THEN
            LL=MAT(KN(IND2*(NSURF+1)))
            A1=DHARM(DS,DS,SGD(L,1),SGD(LL,1))*SIDEB
         ELSE IF(IND2.EQ.-1) THEN
            A1=DHARM(DS,DS,SGD(L,1),DS*QFR(NUM1+IX)/2.0)*SIDEB
         ELSE IF(IND2.EQ.-2) THEN
            A1=0.0D0
         ELSE IF(IND2.EQ.-3) THEN
            A1=2.0D0*DHARM(DS,DS,SGD(L,1),SGD(L,1))*SIDEB
         ENDIF
         VAR1=VAR1+A1
         IF(IND2.GT.0) THEN
            IF(IND2.LT.IND1) THEN
               KEY=KEY0+IND2
               SYS(KEY)=SYS(KEY)-REAL(A1)
            ENDIF
         ENDIF
   20    CONTINUE
         KEY=KEY0+IND1
         SYS(KEY)=SYS(KEY)+REAL(VAR1)+SGD(L,3)*VOL0
   30    NUM1=NUM1+NSURF+1
   35    CONTINUE
      ELSE
*        CROSS SECTION SYSTEM MATRIX ASSEMBLY.
         NUM1=0
         DO 45 IND1=1,LL4
         KHEX=KN(NUM1+NSURF+1)
         IF(VOL(KHEX).EQ.0.0) GO TO 40
         L=MAT(KHEX)
         KEY=MU(IND1)
         SYS(KEY)=SYS(KEY)+SGD(L,1)*QFR(NUM1+NSURF+1)
   40    NUM1=NUM1+NSURF+1
   45    CONTINUE
      ENDIF
      RETURN
      END
