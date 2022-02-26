*DECK PNFLV
      SUBROUTINE PNFLV(KPSYS,INCONV,NGIND,IPTRK,IMPX,NGEFF,NREG,NBMIX,
     1 NUN,MAT,VOL,KEYFLX,FUNKNO,SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve N-group transport equation for fluxes using the spherical
* harmonics (PN) method in BIVAC.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
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
* INCONV  energy group convergence flag (set to .false. if converged).
* NGIND   energy group indices assign to the NGEFF set.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IMPX    print flag (equal to zero for no print).
* NGEFF   number of energy groups processed in parallel.
* NREG    total number of regions for which specific values of the
*         neutron flux and reactions rates are required.
* NBMIX   number of mixtures.
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
      INTEGER     NGEFF,NGIND(NGEFF),IMPX,NREG,NBMIX,NUN,MAT(NREG),
     1            KEYFLX(NREG)
      LOGICAL     INCONV(NGEFF)
      REAL        VOL(NREG),FUNKNO(NUN,NGEFF),SUNKNO(NUN,NGEFF)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,NSTATE=40,MAXIT=20,EPSINR=1.0E-5,ICL1=3,
     1            ICL2=3)
      INTEGER     IPAR(NSTATE)
      DOUBLE PRECISION F1,F2,R1,R2,DMU
      CHARACTER   NAMP*12
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MU,KN,IPERT
      REAL, ALLOCATABLE, DIMENSION(:) :: QFR,SGDI,RR,VV,SYS,OLD1,OLD2,
     1 XX,YY
*----
*  RECOVER PN SPECIFIC PARAMETERS.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
      IF(NREG.NE.IPAR(1)) CALL XABORT('PNFLV: INVALID VALUE OF NREG.')
      IF(NUN.NE.IPAR(2)) CALL XABORT('PNFLV: INVALID VALUE OF NUN.')
      ITYPE=IPAR(6)
      IELEM=IPAR(8)
      ICOL=IPAR(9)
      ISPLH=IPAR(10)
      L4=IPAR(11)
      LX=IPAR(12)
      NLF=IPAR(14)
      ISPN=IPAR(15)
      ISCAT=IPAR(16)
      NVD=IPAR(17)
*----
*  RECOVER TRACKING INFORMATION.
*----
      ALLOCATE(MU(L4))
      CALL LCMGET(IPTRK,'MU',MU)
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      CALL LCMLEN(IPTRK,'QFR',MAXQF,ITYLCM)
      ALLOCATE(KN(MAXKN),QFR(MAXQF))
      CALL LCMGET(IPTRK,'KN',KN)
      CALL LCMGET(IPTRK,'QFR',QFR)
      IIMAX=MU(L4)*(NLF/2)
      ALLOCATE(SYS(IIMAX))
*----
*  RECOVER THE FINITE ELEMENT UNIT STIFFNESS MATRIX.
*----
      CALL LCMSIX(IPTRK,'BIVCOL',1)
      CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
      ALLOCATE(RR(LC*LC),VV(LC*(LC-1)))
      CALL LCMGET(IPTRK,'R',RR)
      CALL LCMGET(IPTRK,'V',VV)
      CALL LCMSIX(IPTRK,' ',2)
*----
*  MAIN LOOP OVER ENERGY GROUPS.
*----
      ALLOCATE(OLD1(NUN),OLD2(NUN))
      DO 140 II=1,NGEFF
      IF(.NOT.INCONV(II)) GO TO 140
      IF(IMPX.GT.1) WRITE(IUNOUT,'(/24H PNFLV: PROCESSING GROUP,I5,
     1 6H WITH ,A,1H.)') NGIND(II),'BIVAC/PN'
*----
*  RECOVER CROSS SECTIONS.
*----
      CALL LCMGET(KPSYS(II),'STATE-VECTOR',IPAR)
      NAN=IPAR(8)
      ALLOCATE(SGDI(NBMIX*NAN))
      DO 10 IL=0,NAN-1
      WRITE(NAMP,'(4HSCAI,I2.2,6H001001)') IL
      CALL LCMGET(KPSYS(II),NAMP,SGDI(IL*NBMIX+1))
   10 CONTINUE
*----
*  MULTIPLICATION OF THE EVEN-PARITY SOURCES BY THE VOLUMES.
*----
      DO 20 K=1,NREG
      IF(MAT(K).EQ.0) GO TO 20
      JND1=KEYFLX(K)
      SUNKNO(JND1,II)=SUNKNO(JND1,II)*VOL(K)
   20 CONTINUE
*----
*  INNER ITERATION LOOP FOR ONE-GROUP TRANSPORT EQUATION.
*----
      CALL LCMGET(KPSYS(II),'IA001001',SYS)
      CALL XDRSET(OLD2,NUN,0.0)
      TEST=0.0
      ITER=0
   30 ITER=ITER+1
      IF(ITER.GT.MAXIT) THEN
         WRITE(IUNOUT,'(43H PNFLV: MAXIMUM NUMBER OF ONE-SPEED ITERATI,
     1   11HON REACHED.)')
         GO TO 100
      ENDIF
      DO 40 I=1,NUN
      OLD1(I)=OLD2(I)
      OLD2(I)=FUNKNO(I,II)
   40 CONTINUE
      IF((ITYPE.EQ.2).OR.((ITYPE.EQ.5).AND.(ISPN.EQ.1))) THEN
         ALLOCATE(XX(NREG),YY(NREG))
         CALL LCMGET(IPTRK,'XX',XX)
         CALL LCMGET(IPTRK,'YY',YY)
         CALL PNFL2E(NREG,IELEM,ICOL,XX,YY,MAT,VOL,NBMIX,NLF,NVD,
     1   NAN,SGDI,L4,KN,QFR,MU,IIMAX,LC,RR,VV,SYS,SUNKNO(1,II),
     2   FUNKNO(1,II),1)
         DEALLOCATE(YY,XX)
      ELSE IF((ITYPE.EQ.2).OR.((ITYPE.EQ.8).AND.(ISPN.EQ.1))) THEN
         NBLOS=LX/3
         ALLOCATE(IPERT(NBLOS))
         CALL LCMGET(IPTRK,'SIDE',SIDE)
         CALL LCMGET(IPTRK,'IPERT',IPERT)
         CALL PNFH2E(IELEM,ICOL,NBLOS,SIDE,NLF,NVD,L4,IPERT,KN,QFR,MU,
     1   IIMAX,LC,VV,SYS,SUNKNO(1,II),FUNKNO(1,II),1)
         DEALLOCATE(IPERT)
      ELSE
         CALL XABORT('PNFLV: TYPE OF DISCRETIZATION NOT IMPLEMENTED.')
      ENDIF
      IF(NLF.EQ.2) GO TO 100
*----
*  VARIATIONAL ACCELERATION.
*----
      DMU=1.0D0
      IF(MOD(ITER-1,ICL1+ICL2).GE.ICL1) THEN
         F1=0.0
         F2=0.0
         DO  70 I=1,NUN
         R1=OLD2(I)-OLD1(I)
         R2=FUNKNO(I,II)-OLD2(I)
         F1=F1+R1*(R2-R1)
         F2=F2+(R2-R1)*(R2-R1)
   70    CONTINUE
         DMU=-F1/F2
         IF(DMU.GT.0.0) THEN
            DO  80 I=1,NUN
            FUNKNO(I,II)=OLD2(I)+REAL(DMU)*(FUNKNO(I,II)-OLD2(I))
            OLD2(I)=OLD1(I)+REAL(DMU)*(OLD2(I)-OLD1(I))
   80       CONTINUE
         ENDIF
      ENDIF
*---  CALCULATE ERROR AND TEST FOR CONVERGENCE.
      AAA=0.0
      BBB=0.0
      DO 90 I=1,NREG
      IF(KEYFLX(I).EQ.0) GO TO 90
      AAA=MAX(AAA,ABS(FUNKNO(KEYFLX(I),II)-OLD2(KEYFLX(I))))
      BBB=MAX(BBB,ABS(FUNKNO(KEYFLX(I),II)))
   90 CONTINUE
      IF(IMPX.GT.2) WRITE(IUNOUT,300) ITER,AAA,BBB,DMU
      IF(AAA.LE.EPSINR*BBB) GO TO 100
      IF(ITER.EQ.1) TEST=AAA
      IF((ITER.GE.10).AND.(AAA.GT.TEST)) THEN
         WRITE(IUNOUT,'(43H PNFLV: UNABLE TO CONVERGE ONE-SPEED ITERAT,
     1   5HIONS.)')
         GO TO 100
      ENDIF
      GO TO 30
  100 DEALLOCATE(SGDI)
*----
*  DIVISION OF THE EVEN-PARITY SOURCES BY THE VOLUMES.
*----
      DO 130 K=1,NREG
      IF(MAT(K).EQ.0) GO TO 130
      JND1=KEYFLX(K)
      SUNKNO(JND1,II)=SUNKNO(JND1,II)/VOL(K)
  130 CONTINUE
*----
*  END OF LOOP OVER ENERGY GROUPS.
*----
  140 CONTINUE
      DEALLOCATE(OLD2,OLD1,SYS,VV,RR,QFR,KN,MU)
      IF((IMPX.GT.0).AND.(NLF.GT.2)) WRITE(IUNOUT,'(15H PNFLV: NUMBER ,
     1 24HOF ONE-SPEED ITERATIONS=,I5,1H.)') ITER
      RETURN
*
  300 FORMAT(27H PNFLV: ONE-SPEED ITERATION,I3,8H  ERROR=,1P,E11.4,
     1 5H OVER,E11.4,22H  ACCELERATION FACTOR=,0P,F7.3)
      END
