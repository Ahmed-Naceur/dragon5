*DECK SYBALC
      SUBROUTINE SYBALC(NPIJ,MAXPTS,RAYRE,SIG,NGAUSS,ALBEDO,Z,PIJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Pij calculation in 1D cylindrical geometry. The tracking is computed
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
      PARAMETER (MKI3=600,PI=3.1415926535,ZI30=0.785398164)
      LOGICAL LGEMPT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: AUXI
      COMMON /BICKL3/BI3(0:MKI3),BI31(0:MKI3),BI32(0:MKI3),PAS3,XLIM3,L3
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
      AUXI(IP,1)=PI*RAYRE(IP+1)**2-VOLI
      AUXI(IP,2)=MAX(1.0E-10,SIG(IP))
      VOLI=PI*RAYRE(IP+1)**2
   10 CONTINUE
      SURF=2.0*PI*RAYRE(NPIJ+1)
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
         DII=AUX0-ZI30+TABKI(3,AUX0)
         PIJ(IX,IX)=PIJ(IX,IX)+W*DII/AUXI(IX,2)**2
         TAU=AUX0
         TAU1J=0.0
         DO 40 IP=IX+1,NPIJ
         AUX1=AUXI(IP,3)
         AUX2=TAU+2.0*AUXI(IP,3)
         DII=AUX1-ZI30+TABKI(3,AUX1)
         TAB0=TABKI(3,TAU)
         TAB1=TABKI(3,TAU+AUX1)
         TAB2=TABKI(3,AUX2)
         CII=TAB0-2.0*TAB1+TAB2
         TAB0=TABKI(3,TAU1J)
         TAB1=TABKI(3,TAU1J+AUX0)
         TAB2=TABKI(3,TAU1J+AUX1)
         TAB3=TABKI(3,TAU1J+AUX0+AUX1)
         CIJ1=TAB0-TAB1-TAB2+TAB3
         PIJ(IP,IP)=PIJ(IP,IP)+W*(2.0*DII+CII)/AUXI(IP,2)**2
         PIJ(IX,IP)=PIJ(IX,IP)+W*CIJ1/(AUXI(IX,2)*AUXI(IP,2))
         TAUIJ=0.0
         DO 30 JP=IP+1,NPIJ
         AUX3=AUXI(JP,3)
         IF(TAUIJ+AUX1+AUX3.GE.XLIM3) THEN
            TAB0=TABKI(3,TAUIJ)
            TAB1=TABKI(3,TAUIJ+AUX1)
            TAB2=TABKI(3,TAUIJ+AUX3)
            TAB3=TABKI(3,TAUIJ+AUX1+AUX3)
         ELSE
            K=NINT(TAUIJ*PAS3)
            TAB0=BI3(K)+TAUIJ*(BI31(K)+TAUIJ*BI32(K))
            TAUX=TAUIJ+AUX1
            K=NINT(TAUX*PAS3)
            TAB1=BI3(K)+TAUX*(BI31(K)+TAUX*BI32(K))
            TAUX=TAUIJ+AUX3
            K=NINT(TAUX*PAS3)
            TAB2=BI3(K)+TAUX*(BI31(K)+TAUX*BI32(K))
            TAUX=TAUIJ+AUX1+AUX3
            K=NINT(TAUX*PAS3)
            TAB3=BI3(K)+TAUX*(BI31(K)+TAUX*BI32(K))
         ENDIF
         CIJ2=TAB0-TAB1-TAB2+TAB3
         IF(TAUIJ+AUX2+AUX3.GE.XLIM3) THEN
            TAB0=TABKI(3,TAUIJ+TAU+AUX1)
            TAB1=TABKI(3,TAUIJ+AUX2)
            TAB2=TABKI(3,TAUIJ+TAU+AUX1+AUX3)
            TAB3=TABKI(3,TAUIJ+AUX2+AUX3)
         ELSE
            TAUX=TAUIJ+TAU+AUX1
            K=NINT(TAUX*PAS3)
            TAB0=BI3(K)+TAUX*(BI31(K)+TAUX*BI32(K))
            TAUX=TAUIJ+AUX2
            K=NINT(TAUX*PAS3)
            TAB1=BI3(K)+TAUX*(BI31(K)+TAUX*BI32(K))
            TAUX=TAUIJ+TAU+AUX1+AUX3
            K=NINT(TAUX*PAS3)
            TAB2=BI3(K)+TAUX*(BI31(K)+TAUX*BI32(K))
            TAUX=TAUIJ+AUX2+AUX3
            K=NINT(TAUX*PAS3)
            TAB3=BI3(K)+TAUX*(BI31(K)+TAUX*BI32(K))
         ENDIF
         CIJ3=TAB0-TAB1-TAB2+TAB3
         PIJ(IP,JP)=PIJ(IP,JP)+W*(CIJ2+CIJ3)/(AUXI(IP,2)*AUXI(JP,2))
         TAUIJ=TAUIJ+AUX3
   30    CONTINUE
         TAU=AUX2
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
         CALL SYB33C(DII,2.0*AUXI(IX,3))
         PIJ(IX,IX)=PIJ(IX,IX)+W*DII/AUXI(IX,2)**2
         TAU=2.0*AUXI(IX,3)
         TAU1J=0.0
         DO 70 IP=IX+1,NPIJ
         CALL SYB33C(DII,AUXI(IP,3))
         CALL SYB31C(CII,TAU,AUXI(IP,3),AUXI(IP,3))
         CALL SYB31C(CIJ1,TAU1J,2.0*AUXI(IX,3),AUXI(IP,3))
         PIJ(IP,IP)=PIJ(IP,IP)+W*(2.0*DII+CII)/AUXI(IP,2)**2
         PIJ(IX,IP)=PIJ(IX,IP)+W*CIJ1/(AUXI(IX,2)*AUXI(IP,2))
         TAUIJ=0.0
         DO 60 JP=IP+1,NPIJ
         CALL SYB31C(CIJ2,TAUIJ,AUXI(IP,3),AUXI(JP,3))
         CALL SYB31C(CIJ3,TAUIJ+TAU+AUXI(IP,3),AUXI(IP,3),AUXI(JP,3))
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
