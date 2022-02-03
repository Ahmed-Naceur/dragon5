*DECK TRIFLV
      SUBROUTINE TRIFLV(KPSYS,INCONV,NGIND,IPTRK,IMPX,NGEFF,NREG,NUN,
     1 MAT,VOL,KEYFLX,FUNKNO,SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve N-group transport equation for fluxes using the diffusion
* approximation or simplified PN method in TRIVAC.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* KPSYS   pointer to the assembly matrices. KPSYS is an array of
*         directories.
* INCONV  energy group convergence flag (set to .FALSE. if converged).
* NGIND   energy group indices assign to the NGEFF set.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IMPX    print flag (equal to zero for no print).
* NGEFF   number of energy groups processed in parallel.
* NREG    total number of regions for which specific values of the
*         neutron flux and reactions rates are required.
* NUN     total number of unknowns in vectors SUNKNO and FUNKNO.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  position of averaged flux elements in FUNKNO vector.
* SUNKNO  input source vector.
*
*Parameters: input/output
* FUNKNO  unknown vector.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KPSYS(NGEFF),IPTRK
      INTEGER     NGEFF,NGIND(NGEFF),IMPX,NREG,NUN,MAT(NREG),
     1            KEYFLX(NREG)
      LOGICAL     INCONV(NGEFF)
      REAL        VOL(NREG),FUNKNO(NUN,NGEFF),SUNKNO(NUN,NGEFF)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,NSTATE=40,MAXIT=20,EPSINR=1.0E-5,ICL1=3,
     1            ICL2=3)
      DOUBLE PRECISION F1,F2,R1,R2,DMU
      INTEGER     IPAR(NSTATE)
      CHARACTER   NAMP*12
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR,OLD1,OLD2
*----
*  RECOVER TRIVAC SPECIFIC PARAMETERS.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
      ITY=2
      IELEM=ABS(IPAR(9))
      LL4=IPAR(11)
      ISPLH=IPAR(13)
      LX=IPAR(14)
      LZ=IPAR(16)
      NLF=IPAR(30)
      IF(IPAR(12).EQ.2) ITY=3
      IF((NLF.GT.0).AND.(ITY.GE.3)) ITY=10+ITY
      IF((ITY.EQ.11).OR.(ITY.EQ.13)) LL4=LL4*NLF/2
      NADI=IPAR(33)
*----
*  MAIN LOOP OVER ENERGY GROUPS.
*----
      ALLOCATE(GAR(NUN),OLD1(NUN),OLD2(NUN))
      DO 130 II=1,NGEFF
      IF(.NOT.INCONV(II)) GO TO 130
      IF(IMPX.GT.1) WRITE(IUNOUT,'(/25H TRIFLV: PROCESSING GROUP,I5,
     1 6H WITH ,A,1H.)') NGIND(II),'TRIVAC'
*----
*  MULTIPLICATION OF THE SOURCES BY THE VOLUMES.
*----
      DO 10 K=1,NREG
      IF(MAT(K).EQ.0) GO TO 10
      JND1=KEYFLX(K)
      SUNKNO(JND1,II)=SUNKNO(JND1,II)*VOL(K)
   10 CONTINUE
*----
*  SOLVE FOR THE FLUXES. USE EQUATION (C.24) IN IGE-281.
*----
      NAMP='A001001'
      CALL XDRSET(OLD2,NUN,0.0)
      TEST=0.0
      ITER=0
   20 ITER=ITER+1
      IF(ITER.GT.MAXIT) THEN
      WRITE(IUNOUT,'(46H TRIFLV: MAXIMUM NUMBER OF ONE-SPEED ITERATION,
     1   9H REACHED.)')
         GO TO 110
      ENDIF
      DO 30 I=1,NUN
      OLD1(I)=OLD2(I)
      OLD2(I)=FUNKNO(I,II)
   30 CONTINUE
      CALL MTLDLM(NAMP,IPTRK,KPSYS(II),LL4,ITY,FUNKNO(1,II),GAR)
      DO 40 I=1,NUN
      GAR(I)=SUNKNO(I,II)-GAR(I)
   40 CONTINUE
      CALL FLDADI(NAMP,IPTRK,KPSYS(II),LL4,ITY,GAR,NADI)
      DO 50 I=1,NUN
      FUNKNO(I,II)=FUNKNO(I,II)+GAR(I)
   50 CONTINUE
*----
*  VARIATIONAL ACCELERATION.
*----
      DMU=1.0D0
      IF(MOD(ITER-1,ICL1+ICL2).GE.ICL1) THEN
         F1=0.0D0
         F2=0.0D0
         DO  80 I=1,NUN
         R1=OLD2(I)-OLD1(I)
         R2=FUNKNO(I,II)-OLD2(I)
         F1=F1+R1*(R2-R1)
         F2=F2+(R2-R1)*(R2-R1)
   80    CONTINUE
         DMU=-F1/F2
         IF(DMU.GT.0.0) THEN
            DO  90 I=1,NUN
            FUNKNO(I,II)=OLD2(I)+REAL(DMU)*(FUNKNO(I,II)-OLD2(I))
            OLD2(I)=OLD1(I)+REAL(DMU)*(OLD2(I)-OLD1(I))
   90       CONTINUE
         ENDIF
      ENDIF      
*----
*  CALCULATE ERROR AND TEST FOR CONVERGENCE.
*----
      AAA=0.0
      BBB=0.0
      DO 100 I=1,NREG
      IF(KEYFLX(I).EQ.0) GO TO 100
      AAA=MAX(AAA,ABS(FUNKNO(KEYFLX(I),II)-OLD2(KEYFLX(I))))
      BBB=MAX(BBB,ABS(FUNKNO(KEYFLX(I),II)))
  100 CONTINUE
      IF(IMPX.GT.2) WRITE(IUNOUT,300) ITER,AAA,BBB,DMU
      IF(AAA.LE.EPSINR*BBB) GO TO 110
      IF(ITER.EQ.1) TEST=AAA
      IF((ITER.GE.10).AND.(AAA.GT.TEST)) THEN
         WRITE(IUNOUT,'(43H TRIFLV: UNABLE TO CONVERGE ONE-SPEED ITERA,
     1   6HTIONS.)')
         GO TO 110
      ENDIF
      GO TO 20
*----
*  DIVISION OF THE SOURCES BY THE VOLUMES.
*----
  110 DO 120 K=1,NREG
      IF(MAT(K).EQ.0) GO TO 120
      JND1=KEYFLX(K)
      SUNKNO(JND1,II)=SUNKNO(JND1,II)/VOL(K)
  120 CONTINUE
*----
* END OF LOOP OVER ENERGY GROUPS
*----
  130 CONTINUE
      DEALLOCATE(OLD2,OLD1,GAR)
      RETURN
*
  300 FORMAT(28H TRIFLV: ONE-SPEED ITERATION,I3,8H  ERROR=,1P,E11.4,
     1 5H OVER,E11.4,22H  ACCELERATION FACTOR=,0P,F7.3)
      END
