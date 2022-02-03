*DECK SYBALS
      SUBROUTINE SYBALS(NPIJ,MAXPTS,RAYRE,SIG,NGAUSS,ALBEDO,Z,PIJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Pij calculation in 1D spherical geometry. The tracking is computed
* by subroutine SYBT1D.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NPIJ    number of regions.
* MAXPTS  first dimension of matrix PIJ.
* RAYRE   radius of regions array.
* SIG     total cross section array.
* NGAUSS  number of Gauss points.
* ALBEDO  outside albedo.
* Z       tracking information.
*
*Parameters: output
* PIJ     reduced collision probability matrix.
*
*Reference:
* A. Kavenoky, 'Calcul et utilisation des probabilites de premiere
* collision pour les milieux heterogenes a une dimension: Les programmes
* ALCOLL et CORTINA', note CEA-N-1077, Commissariat a l'energie
* atomique, mars 1969.
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NPIJ,MAXPTS,NGAUSS
      REAL RAYRE(NPIJ+1),SIG(NPIJ),PIJ(MAXPTS,NPIJ),ALBEDO,Z(*)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (PI=3.1415926535)
      LOGICAL LGEMPT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: AUXI
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(AUXI(NPIJ,3))
*----
*  TEST FOR VOIDED REGIONS
*----
      LGEMPT=.FALSE.
      VOLI=PI*RAYRE(1)**2
      DO 10 IP=1,NPIJ
      LGEMPT=LGEMPT.OR.(2.0*(RAYRE(IP+1)-RAYRE(IP))*SIG(IP).LE.0.004)
      AUXI(IP,1)=(4.0/3.0)*PI*RAYRE(IP+1)**3-VOLI
      AUXI(IP,2)=MAX(1.0E-10,SIG(IP))
      VOLI=(4.0/3.0)*PI*RAYRE(IP+1)**3
   10 CONTINUE
      SURF=4.0*PI*RAYRE(NPIJ+1)**2
      CALL XDRSET(PIJ,MAXPTS*NPIJ,0.0)
      IZ=1
      IF(.NOT.LGEMPT) THEN
*        NO VOIDED REGIONS DETECTED.
         DO 42 IX=1,NPIJ
         DO 41 I=1,NGAUSS
         IZ=IZ+2
         W=Z(IZ)
         DO 20 ITR=IX,NPIJ
         IZ=IZ+1
         AUXI(ITR,3)=AUXI(ITR,2)*Z(IZ)
   20    CONTINUE
         AUX0=2.0*AUXI(IX,3)
         EXP0=EXP(-AUX0)
         DII=AUX0-1.0+EXP0
         PIJ(IX,IX)=PIJ(IX,IX)+W*DII/AUXI(IX,2)**2
         TAU=AUX0
         TAU1J=0.0
         DO 40 IP=IX+1,NPIJ
         AUX1=AUXI(IP,3)
         EXP1=EXP(-TAU)
         EXP2=EXP(-AUX1)
         EXP3=EXP(-TAU1J)
         DII=AUX1-1.0+EXP2
         CII=EXP1*(1.0-2.0*EXP2+EXP2*EXP2)
         CIJ1=EXP3*(1.0-EXP0-EXP2+EXP0*EXP2)
         PIJ(IP,IP)=PIJ(IP,IP)+W*(2.0*DII+CII)/AUXI(IP,2)**2
         PIJ(IX,IP)=PIJ(IX,IP)+W*CIJ1/(AUXI(IX,2)*AUXI(IP,2))
         TAUIJ=0.0
         DO 30 JP=IP+1,NPIJ
         EXP4=EXP(-TAUIJ)
         EXP5=EXP(-AUXI(JP,3))
         CIJ2=EXP4*(1.0-EXP2-EXP5+EXP2*EXP5)
         CIJ3=EXP1*EXP2*EXP4*(1.0-EXP2-EXP5+EXP2*EXP5)
         PIJ(IP,JP)=PIJ(IP,JP)+W*(CIJ2+CIJ3)/(AUXI(IP,2)*AUXI(JP,2))
         TAUIJ=TAUIJ+AUXI(JP,3)
   30    CONTINUE
         TAU=TAU+2.0*AUX1
         TAU1J=TAU1J+AUX1
   40    CONTINUE
   41    CONTINUE
   42    CONTINUE
      ELSE
         DO 72 IX=1,NPIJ
         DO 71 I=1,NGAUSS
         IZ=IZ+2
         W=Z(IZ)
         DO 50 ITR=IX,NPIJ
         IZ=IZ+1
         AUXI(ITR,3)=AUXI(ITR,2)*Z(IZ)
   50    CONTINUE
         CALL SYB43C(DII,2.0*AUXI(IX,3))
         PIJ(IX,IX)=PIJ(IX,IX)+W*DII/AUXI(IX,2)**2
         TAU=2.0*AUXI(IX,3)
         TAU1J=0.0
         DO 70 IP=IX+1,NPIJ
         CALL SYB43C(DII,AUXI(IP,3))
         CALL SYB41C(CII,TAU,AUXI(IP,3),AUXI(IP,3))
         CALL SYB41C(CIJ1,TAU1J,2.0*AUXI(IX,3),AUXI(IP,3))
         PIJ(IP,IP)=PIJ(IP,IP)+W*(2.0*DII+CII)/AUXI(IP,2)**2
         PIJ(IX,IP)=PIJ(IX,IP)+W*CIJ1/(AUXI(IX,2)*AUXI(IP,2))
         TAUIJ=0.0
         DO 60 JP=IP+1,NPIJ
         CALL SYB41C(CIJ2,TAUIJ,AUXI(IP,3),AUXI(JP,3))
         CALL SYB41C(CIJ3,TAUIJ+TAU+AUXI(IP,3),AUXI(IP,3),AUXI(JP,3))
         PIJ(IP,JP)=PIJ(IP,JP)+W*(CIJ2+CIJ3)/(AUXI(IP,2)*AUXI(JP,2))
         TAUIJ=TAUIJ+AUXI(JP,3)
   60    CONTINUE
         TAU=TAU+2*AUXI(IP,3)
         TAU1J=TAU1J+AUXI(IP,3)
   70    CONTINUE
   71    CONTINUE
   72    CONTINUE
      ENDIF
*
      DO 85 I=1,NPIJ
      DO 80 J=I,NPIJ
      VAL=PIJ(I,J)
      PIJ(I,J)=VAL/AUXI(I,1)
      PIJ(J,I)=VAL/AUXI(J,1)
   80 CONTINUE
   85 CONTINUE
*----
*  COMPUTING REFLECTED PROBABILITIES ASSUMING WHITE BOUNDARY CONDITION.
*----
      IF(ALBEDO.NE.0.0) THEN
         PSS=1.0
         DO 100 IK=1,NPIJ
         AUXI(IK,3)=1.0
         DO 90 JK=1,NPIJ
         AUXI(IK,3)=AUXI(IK,3)-PIJ(IK,JK)*AUXI(JK,2)
   90    CONTINUE
         PSS=PSS-4.0*AUXI(IK,1)*AUXI(IK,2)*AUXI(IK,3)/SURF
  100    CONTINUE
         AUX0=ALBEDO/(1.0-ALBEDO*PSS)
         DO 120 JK=1,NPIJ
         AUX1=AUX0*(4.0*AUXI(JK,1)/SURF)*AUXI(JK,3)
         DO 110 IK=1,NPIJ
         PIJ(IK,JK)=PIJ(IK,JK)+AUXI(IK,3)*AUX1
  110    CONTINUE
  120    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(AUXI)
      RETURN
      END
