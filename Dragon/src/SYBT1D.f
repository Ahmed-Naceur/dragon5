*DECK SYBT1D
      SUBROUTINE SYBT1D(NPIJ,RAD,LGSPH,NGAUSS,Z)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Tracking for pij calculation using the method of Kavenoky in annular
* or spherical geometry. The tracking is used by SYBALC.
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
* NPIJ    number of regions.
* RAD     radius of regions array.
* LGSPH   geometry flag (=.TRUE. for spherical geometry).
* NGAUSS  number of Gauss points.
*
*Parameters: output
* Z       tracking information.
*
*-----------------------------------------------------------------------
*
* REFERENCE:
* A. KAVENOKY, 'CALCUL ET UTILISATION DES PROBABILITES DE PREMIERE
* COLLISION POUR LES MILIEUX HETEROGENES A UNE DIMENSION: LES PROGRAMMES
* ALCOLL ET CORTINA', NOTE CEA-N-1077, COMMISSARIAT A L'ENERGIE
* ATOMIQUE, MARS 1969.
*----
*  SUBROUTINE ARGUMENTS
*----
      LOGICAL LGSPH
      INTEGER NPIJ,NGAUSS
      REAL RAD(0:NPIJ),Z(*)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(MAXGAU=64,PI=3.1415926535)
      REAL ALPR(MAXGAU),PWR(MAXGAU)
*
      CALL ALGJP(NGAUSS,ALPR,PWR)
      SUM=0.0
      IZ=1
      RIK1=RAD(0)
      RI=1.0/RAD(NPIJ)
      DO 60 IK=1,NPIJ
      RIK2=RAD(IK)
      RD=RIK2-RIK1
      DO 50 IL=1,NGAUSS
      R=RIK2-RD*ALPR(IL)**2
      R2=R*R
      IZ=IZ+2
      IZ1=IZ
*
*----
*  STORE INTEGRATION WEIGHT (TIMES 2) FOR THIS LINE
*----
      AUX=RD*PWR(IL)*4.0
      IF(LGSPH) AUX=AUX*R*PI
      Z(IZ)=AUX
      CT1=0.
      CT2=0.
      DO 40 I=IK,NPIJ
      CT2=SQRT(RAD(I)**2-R2)
      IZ=IZ+1
*----
*  STORE LENGTH OF PATH
*----
      Z(IZ)=CT2-CT1
      CT1=CT2
   40 CONTINUE
*----
*  STORE COS(PHI)*INTEGRATION WEIGHT ( TIMES 2 )
*----
      Z(IZ1-1)=Z(IZ1)*CT2*RI
      SUM=SUM+AUX/CT2
   50 CONTINUE
      RIK1=RIK2
   60 CONTINUE
      IZ=IZ+1
      IF(LGSPH) THEN
         Z(1)=SUM/(PI*2.0*RAD(NPIJ))
      ELSE
         Z(1)=SUM/PI
      ENDIF
      RETURN
      END
