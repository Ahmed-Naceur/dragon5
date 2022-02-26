*DECK SYBALP
      SUBROUTINE SYBALP(NPIJ,MAXPTS,Y,SIG,NCOD,ALB,PIJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Pij calculation using the method of Kavenoky in 1D slab geometry.
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
* Y       abscissa array.
* SIG     total cross section array.
* NCOD    left and right type of boundary conditions (=1 void;
*         =2 refl; =4 tran).
* ALB     left and right albedos.
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
      INTEGER NPIJ,MAXPTS,NCOD(2)
      REAL Y(NPIJ+1),SIG(NPIJ),ALB(2),PIJ(MAXPTS,NPIJ)
*----
*  LOCAL VARIABLES
*----
      CHARACTER BC*8
      REAL, ALLOCATABLE, DIMENSION(:,:) :: AUXI,F2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(AUXI(2*NPIJ,3),F2(2*NPIJ,2*NPIJ))
*----
*  SET THE BOUNDARY CONDITIONS
*----
      IF((NCOD(1).EQ.1).AND.(NCOD(2).EQ.1)) THEN
         BC='VOID'
         IF((ALB(1).NE.0.0).OR.(ALB(2).NE.0.0)) BC='ALBE'
      ELSE IF((NCOD(1).EQ.4).AND.(NCOD(2).EQ.4)) THEN
         BC='TRAN'
      ELSE IF((ALB(1).EQ.0.0).AND.(ALB(2).EQ.0.0)) THEN
         BC='VOID'
      ELSE
         BC='ALBE'
      ENDIF
*----
*  COMPUTE THE PIJ MATRIX
*----
      CALL XDRSET(PIJ,MAXPTS*NPIJ,0.0)
      IF(BC.EQ.'VOID') THEN
         DO 10 IP=1,NPIJ
         AUXI(IP,1)=Y(IP+1)-Y(IP)
         AUXI(IP,3)=MAX(1.0E-10,AUXI(IP,1)*SIG(IP))
         AUXI(IP,2)=AUXI(IP,3)/AUXI(IP,1)
   10    CONTINUE
         DO 25 IP=1,NPIJ
         CALL SYBRII(RIIP,1.0,0.0,AUXI(IP,3))
         PIJ(IP,IP)=RIIP/AUXI(IP,2)**2
         TAU0=0.0
         DO 20 JP=IP+1,NPIJ
         CALL SYBRIJ(RIJP,1.0,TAU0,AUXI(IP,3),AUXI(JP,3))
         PIJ(IP,JP)=RIJP/(AUXI(IP,2)*AUXI(JP,2))
         TAU0=TAU0+AUXI(JP,3)
   20    CONTINUE
   25    CONTINUE
         DO 35 IP=1,NPIJ
         DO 30 JP=IP,NPIJ
         PIJ(JP,IP)=PIJ(IP,JP)
   30    CONTINUE
   35    CONTINUE
      ELSE IF(BC.EQ.'TRAN') THEN
         DO 40 IP=1,NPIJ
         AUXI(IP,1)=Y(IP+1)-Y(IP)
         AUXI(IP,3)=MAX(1.0E-10,AUXI(IP,1)*SIG(IP))
         AUXI(IP,2)=AUXI(IP,3)/AUXI(IP,1)
   40    CONTINUE
         TAUCEL=0.0
         DO 50 IP=1,NPIJ
         TAUCEL=TAUCEL+AUXI(IP,3)
   50    CONTINUE
         M=-1
   60    M=M+1
         IF(M.GT.100) CALL XABORT('SYBALP: UNABLE TO CONVERGE(1).')
         CALL XDRSET(F2,4*NPIJ*NPIJ,0.0)
         SMALL=0.0
         DO 75 IP=1,NPIJ
         CALL SYBRII(RIIP,1.0,M*TAUCEL,AUXI(IP,3))
         CALL SYBRII(RIIM,-1.0,(M+1)*TAUCEL,AUXI(IP,3))
         F2(IP,IP)=F2(IP,IP)+(RIIP+RIIM)/AUXI(IP,2)**2
         SMALL=MAX(SMALL,ABS(F2(IP,IP)*AUXI(IP,2)))
         TAU0=0.0
         DO 70 JP=IP+1,NPIJ
         CALL SYBRIJ(RIJP,1.0,M*TAUCEL+TAU0,AUXI(IP,3),AUXI(JP,3))
         CALL SYBRIJ(RIJM,-1.0,(M+1)*TAUCEL-TAU0,AUXI(IP,3),AUXI(JP,3))
         F2(IP,JP)=F2(IP,JP)+(RIJP+RIJM)/(AUXI(IP,2)*AUXI(JP,2))
         TAU0=TAU0+AUXI(JP,3)
         SMALL=MAX(SMALL,ABS(F2(IP,JP)*AUXI(JP,2)))
   70    CONTINUE
   75    CONTINUE
         DO 85 IP=1,NPIJ
         DO 80 JP=IP,NPIJ
         PIJ(IP,JP)=PIJ(IP,JP)+F2(IP,JP)
   80    CONTINUE
   85    CONTINUE
         IF(SMALL.LE.1.0E-6) GO TO 90
         GO TO 60
   90    DO 105 IP=1,NPIJ
         DO 100 JP=IP,NPIJ
         PIJ(JP,IP)=PIJ(IP,JP)
  100    CONTINUE
  105    CONTINUE
      ELSE IF(BC.EQ.'ALBE') THEN
         TAUCEL=0.0
         DO 110 IP=1,NPIJ
         AUXI(IP,1)=Y(IP+1)-Y(IP)
         AUXI(IP,3)=MAX(1.0E-10,AUXI(IP,1)*SIG(IP))
         AUXI(IP,2)=AUXI(IP,3)/AUXI(IP,1)
         AUXI(2*NPIJ-IP+1,1)=AUXI(IP,1)
         AUXI(2*NPIJ-IP+1,2)=AUXI(IP,2)
         AUXI(2*NPIJ-IP+1,3)=AUXI(IP,3)
         TAUCEL=TAUCEL+2.0*AUXI(IP,3)
  110    CONTINUE
         CALL XDRSET(F2,4*NPIJ*NPIJ,0.0)
         DO 125 IP=1,2*NPIJ
         CALL SYBRII(RIIP,1.0,0.0,AUXI(IP,3))
         F2(IP,IP)=F2(IP,IP)+RIIP/AUXI(IP,2)**2
         TAU0=0.0
         DO 120 JP=IP+1,2*NPIJ
         CALL SYBRIJ(RIJP,1.0,TAU0,AUXI(IP,3),AUXI(JP,3))
         F2(IP,JP)=F2(IP,JP)+RIJP/(AUXI(IP,2)*AUXI(JP,2))
         F2(JP,IP)=F2(JP,IP)+RIJP/(AUXI(IP,2)*AUXI(JP,2))
         TAU0=TAU0+AUXI(JP,3)
  120    CONTINUE
  125    CONTINUE
         DO 135 IP=1,NPIJ
         DO 130 JP=1,NPIJ
         PIJ(IP,JP)=PIJ(IP,JP)+F2(IP,JP)+ALB(2)*F2(2*NPIJ+1-IP,JP)
  130    CONTINUE
  135    CONTINUE
         M=0
  140    M=M+1
         IF(M.GT.100) CALL XABORT('UNABLE TO CONVERGE(2).')
         CALL XDRSET(F2,4*NPIJ*NPIJ,0.0)
         SMALL=0.0
         DO 155 IP=1,2*NPIJ
         CALL SYBRII(RIIP,1.0,M*TAUCEL,AUXI(IP,3))
         F2(IP,IP)=F2(IP,IP)+ALB(1)**M*RIIP/AUXI(IP,2)**2
         SMALL=MAX(SMALL,ABS(F2(IP,IP)*AUXI(IP,2)))
         TAU0=0.0
         DO 150 JP=IP+1,2*NPIJ
         CALL SYBRIJ(RIJP,1.0,M*TAUCEL+TAU0,AUXI(IP,3),AUXI(JP,3))
         F2(IP,JP)=F2(IP,JP)+ALB(1)**M*RIJP/(AUXI(IP,2)*AUXI(JP,2))
         CALL SYBRIJ(RIJM,-1.0,M*TAUCEL-TAU0,AUXI(JP,3),AUXI(IP,3))
         F2(JP,IP)=F2(JP,IP)+ALB(1)**M*RIJM/(AUXI(IP,2)*AUXI(JP,2))
         TAU0=TAU0+AUXI(JP,3)
         SMALL=MAX(SMALL,ABS(F2(IP,JP)*AUXI(JP,2)))
         SMALL=MAX(SMALL,ABS(F2(JP,IP)*AUXI(IP,2)))
  150    CONTINUE
  155    CONTINUE
         DO 165 IP=1,NPIJ
         DO 160 JP=1,NPIJ
         PIJ(IP,JP)=PIJ(IP,JP)+ALB(2)**M*F2(IP,JP)+ALB(2)**(M-1)
     1   *F2(2*NPIJ+1-IP,JP)
  160    CONTINUE
  165    CONTINUE
         CALL XDRSET(F2,4*NPIJ*NPIJ,0.0)
         DO 175 IP=1,NPIJ
         CALL SYBRII(RIIM,-1.0,M*TAUCEL,AUXI(IP,3))
         F2(IP,IP)=F2(IP,IP)+ALB(1)**M*RIIM/AUXI(IP,2)**2
         TAU0=0.0
         DO 170 JP=IP+1,2*NPIJ
         CALL SYBRIJ(RIJM,-1.0,M*TAUCEL-TAU0,AUXI(IP,3),AUXI(JP,3))
         F2(IP,JP)=F2(IP,JP)+ALB(1)**M*RIJM/(AUXI(IP,2)*AUXI(JP,2))
         CALL SYBRIJ(RIJP,1.0,M*TAUCEL+TAU0,AUXI(JP,3),AUXI(IP,3))
         F2(JP,IP)=F2(JP,IP)+ALB(1)**M*RIJP/(AUXI(IP,2)*AUXI(JP,2))
         TAU0=TAU0+AUXI(JP,3)
         SMALL=MAX(SMALL,ABS(F2(IP,JP)*AUXI(JP,2)))
         SMALL=MAX(SMALL,ABS(F2(JP,IP)*AUXI(IP,2)))
  170    CONTINUE
  175    CONTINUE
         DO 185 IP=1,NPIJ
         DO 180 JP=1,NPIJ
         PIJ(IP,JP)=PIJ(IP,JP)+ALB(2)**M*F2(IP,JP)+ALB(2)**(M+1)*
     1   F2(2*NPIJ+1-IP,JP)
  180    CONTINUE
  185    CONTINUE
         IF(SMALL.LE.1.0E-6) GO TO 190
         GO TO 140
      ENDIF
  190 DO 210 IP=1,NPIJ
      DO 200 JP=1,NPIJ
      PIJ(IP,JP)=PIJ(IP,JP)/AUXI(IP,1)
  200 CONTINUE
  210 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(F2,AUXI)
      RETURN
      END
