*DECK SNF
      SUBROUTINE SNF (IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,IDIR,NREG,
     1 NBMIX,NUN,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,TITR,DCUTOFF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve N-group transport equation for fluxes using the discrete
* ordinates (SN) method.
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
* IPSYS   pointer to the assembly LCM object (L_PIJ signature). IPSYS is
*         a list of directories.
* NPSYS   index array pointing to the IPSYS list component corresponding
*         to each energy group. Set to zero if a group is not to be
*         processed. Usually, NPSYS(I)=I.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  not used.
* IMPX    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* IDIR    not used.
* NREG    total number of regions for which specific values of the
*         neutron flux and reactions rates are required.
* NBMIX   number of mixtures.
* NUN     total number of unknowns in vectors SUNKNO and FUNKNO.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  position of averaged flux elements in FUNKNO vector.
* SUNKNO  input source vector.
* TITR    title.
*
*Parameters: input/output
* FUNKNO  unknown vector.
* DCUTOFF energy deposition under the energy cutoff
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER   TITR*72
      TYPE(C_PTR) IPSYS,IPTRK
      INTEGER     NGRP,NPSYS(NGRP),IFTRAK,IMPX,IDIR,NREG,NBMIX,NUN,
     1            MAT(NREG),KEYFLX(NREG)
      REAL        VOL(NREG),FUNKNO(NUN,NGRP),SUNKNO(NUN,NGRP),
     1            DCUTOFF(NREG)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,NSTATE=40)
      TYPE(C_PTR) JPSYS
      INTEGER     IPAR(NSTATE)
      LOGICAL     LBIHET,LIVO
      DOUBLE PRECISION F1,F2,R1,R2,DMU
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INGIND
      REAL, ALLOCATABLE, DIMENSION(:) :: FGAR,TEST
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SUNKN,OLD1,OLD2
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: INCONV
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: KPSYS
*----
*  RECOVER SN SPECIFIC PARAMETERS
*----
      IF(IMPX.GT.2) THEN
        WRITE(IUNOUT,'(//6H SNF: ,A72)') TITR
        CALL KDRCPU(TK1)
      ENDIF
      IF(IDIR.NE.0) CALL XABORT('SNF: EXPECTING IDIR=0')
      IF(IFTRAK.NE.0) CALL XABORT('SNF: EXPECTING IFTRAK=0')
      CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
      IF(NREG.NE.IPAR(1)) CALL XABORT('SNF: INVALID NREG ON LCM.')
      IF(NUN.NE.IPAR(2)) CALL XABORT('SNF: INVALID NUN ON LCM.')
      NSTART=IPAR(20)
      MAXIT=IPAR(22)
      LIVO=(IPAR(23).EQ.1)
      ICL1=IPAR(24)
      ICL2=IPAR(25)
      IBFP=IPAR(31)
      LBIHET=IPAR(40).NE.0
      CALL LCMGET(IPTRK,'EPSI',EPSINR)
      IF(IMPX.GT.3) THEN
        ALLOCATE(FGAR(NREG))
        DO IG=1,NGRP
          IF(NPSYS(IG).EQ.0) CYCLE
          CALL XDRSET(FGAR,NREG,0.0)
          DO I=1,NREG
            IF(KEYFLX(I).NE.0) FGAR(I)=SUNKNO(KEYFLX(I),IG)
          ENDDO
          WRITE(IUNOUT,'(/33H N E U T R O N    S O U R C E S (,I5,
     1    3H ):)') IG
          WRITE(IUNOUT,'(1P,6(5X,E15.7))') (FGAR(I),I=1,NREG)
        ENDDO
        DEALLOCATE(FGAR)
      ENDIF
*----
* DOUBLE HETEROGENEITY TREATMENT
*----
      NBMIXG=0
      NREGAR=0
      CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
      LBIHET=IPAR(40).NE.0
      IF(LBIHET) THEN
         ALLOCATE(SUNKN(NUN,NGRP))
         NBMIXG=NBMIX
         NREGAR=NREG
         DO IG=1,NGRP
           IOFSET=NPSYS(IG)
           IF(IOFSET.NE.0) THEN
             JPSYS=LCMGIL(IPSYS,IOFSET)
             SUNKN(:NUN,IG)=SUNKNO(:NUN,IG)
             CALL DOORFB2(JPSYS,IPTRK,IMPX,NBMIX,NREG,NUN,KEYFLX,
     1       NBMIX2,NREG2,SUNKNO(1,IG))
           ENDIF
         ENDDO
         NBMIX=NBMIX2
         NREG=NREG2
      ENDIF
*---
* DETERMINE THE NUMBER OF GROUPS TO BE PROCESSED
* RECOVER FLUXES FROM A PREVIOUS SELF-SHIELDING CALCULATION IF AVAILABLE
*---
      NGEFF=0
      DO IG=1,NGRP
         IOFSET=NPSYS(IG)
         IF(IOFSET.NE.0) THEN
            NGEFF=NGEFF+1
            JPSYS=LCMGIL(IPSYS,IOFSET)
            CALL LCMLEN(JPSYS,'FUNKNO$USS',ILENG,ITYLCM)
            IF(ILENG.EQ.NUN) THEN
               CALL LCMGET(JPSYS,'FUNKNO$USS',FUNKNO(1,IG))
            ENDIF
         ENDIF
      ENDDO
      IF(NGEFF.EQ.0) GO TO 100
*---
* RECOVER POINTERS TO EACH GROUP PROPERTIES
* CREATE AN INDEX FOR THE GROUPS TO BE PROCESSED
*---
      ALLOCATE(INGIND(NGEFF),KPSYS(NGEFF))
      II=1
      DO IG=1,NGRP
         IOFSET=NPSYS(IG)
         IF(IOFSET.NE.0) THEN
            INGIND(II)=IG
            IF(LBIHET) THEN
               JPSYS=LCMGIL(IPSYS,IOFSET)
               KPSYS(II)=LCMGID(JPSYS,'BIHET')
            ELSE
               KPSYS(II)=LCMGIL(IPSYS,IOFSET)
            ENDIF
            II=II+1
         ENDIF
      ENDDO
*
      IF(NSTART.GT.0) THEN
*----
*  GMRES(M) INNER ITERATION LOOP FOR ONE-SPEED TRANSPORT EQUATION
*----
         CALL SNGMRE (KPSYS,INGIND,IPTRK,IMPX,NGRP,NGEFF,NREG,NBMIX,
     1   NUN,NSTART,MAXIT,EPSINR,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,DCUTOFF)
      ELSE
*----
*  LIVOLANT INNER ITERATION LOOP FOR ONE-SPEED TRANSPORT EQUATION
*----
         ALLOCATE(INCONV(NGEFF),OLD1(NUN,NGEFF),OLD2(NUN,NGEFF),
     1   TEST(NGEFF))
         LNCONV=NGEFF
         CALL XDLSET(INCONV,NGEFF,.TRUE.)
         CALL XDRSET(OLD2,NUN*NGEFF,0.0)
         CALL XDRSET(TEST,NGEFF,0.0)
*
         ITER=0
   10    ITER=ITER+1
         IF(ITER.GT.MAXIT) THEN
            WRITE(IUNOUT,'(40H SNF: MAXIMUM NUMBER OF ONE-SPEED ITERAT,
     1      12HION REACHED.)')
            GO TO 70
         ENDIF
*
         DO II=1,NGEFF
           IG=INGIND(II)
           DO I=1,NUN
             OLD1(I,II)=OLD2(I,II)
             OLD2(I,II)=FUNKNO(I,IG)
           ENDDO
         ENDDO
*----
*  UPDATE THE FIXED SOURCE AND COMPUTE THE FLUX
*----
         CALL SNFLUX(KPSYS,INCONV,INGIND,IPTRK,IMPX,NGRP,NGEFF,NREG,
     1   NBMIX,NUN,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,ITER,DCUTOFF)
*----
*  LOOP OVER ENERGY GROUPS
*----
         DO 60 II=1,NGEFF
         IF(INCONV(II)) THEN
           IG=INGIND(II)
*----
*  VARIATIONAL ACCELERATION. LIVOLANT INNER ITERATION LOOP FOR ONE-GROUP
*  TRANSPORT EQUATION.
*----
           DMU=1.0D0
           IF(LIVO.AND.(MOD(ITER-1,ICL1+ICL2).GE.ICL1)) THEN
              F1=0.0
              F2=0.0
              DO 30 I=1,NUN
              R1=OLD2(I,II)-OLD1(I,II)
              R2=FUNKNO(I,IG)-OLD2(I,II)
              F1=F1+R1*(R2-R1)
              F2=F2+(R2-R1)*(R2-R1)
   30         CONTINUE
              DMU=-F1/F2
              IF(DMU.GT.0.0) THEN
                RDMU=REAL(DMU)
                DO 40 I=1,NUN
                FUNKNO(I,IG)=OLD2(I,II)+RDMU*(FUNKNO(I,IG)-OLD2(I,II))
                OLD2(I,II)=OLD1(I,II)+RDMU*(OLD2(I,II)-OLD1(I,II))
   40           CONTINUE
              ENDIF
           ENDIF
*----
*  CALCULATE ERROR AND TEST FOR CONVERGENCE
*----
           AAA=0.0
           BBB=0.0
           DO 50 I=1,NREG
           IF(KEYFLX(I).EQ.0) GO TO 50
           AAA=MAX(AAA,ABS(FUNKNO(KEYFLX(I),IG)-OLD2(KEYFLX(I),II)))
           BBB=MAX(BBB,ABS(FUNKNO(KEYFLX(I),IG)))
   50      CONTINUE
           IF(IMPX.GT.2) WRITE(IUNOUT,300) IG,ITER,AAA,BBB,AAA/BBB,DMU
           IF(IMPX.GT.5) THEN
             ALLOCATE(FGAR(NREG))
             CALL XDRSET(FGAR,NREG,0.0)
             DO I=1,NREG
              IF(KEYFLX(I).NE.0) FGAR(I)=FUNKNO(KEYFLX(I),IG)
             ENDDO
             WRITE(IUNOUT,'(//33H N E U T R O N    F L U X E S   :)')
             WRITE(IUNOUT,'(1P,6(5X,E15.7))') (FGAR(I),I=1,NREG)
             DEALLOCATE(FGAR)
           ENDIF
           IF(AAA.LE.0.1*EPSINR*BBB) THEN
             LNCONV=LNCONV-1
             INCONV(II)=.FALSE.
           ENDIF
           IF(ITER.EQ.1) TEST(II)=AAA
           IF((ITER.GE.10).AND.(AAA.GT.TEST(II))) THEN
             WRITE(IUNOUT,'(39H SNF: UNABLE TO CONVERGE ONE-SPEED ITER,
     1       15HATIONS IN GROUP,I5,1H.)') IG
             LNCONV=LNCONV-1
             INCONV(II)=.FALSE.
           ENDIF
         ENDIF
   60    CONTINUE
         IF(LNCONV.EQ.0) GO TO 70
         GO TO 10
*----
* CONVERGENCE OF ONE-SPEED ITERATIONS IN ALL NGEFF GROUPS
*----
   70    IF(IMPX.GT.1) WRITE(IUNOUT,'(29H SNF: NUMBER OF ONE-SPEED ITE,
     1   8HRATIONS=,I5,1H.)') ITER
         DEALLOCATE(TEST,OLD2,OLD1,INCONV)
      ENDIF
      DEALLOCATE(KPSYS,INGIND)
*----
* DOUBLE HETEROGENEITY TREATMENT
*----
  100 IF(LBIHET) THEN
         DO IG=1,NGRP
           IOFSET=NPSYS(IG)
           IF(IOFSET.NE.0) THEN
             SUNKNO(:NUN,IG)=SUNKN(:NUN,IG)
             JPSYS=LCMGIL(IPSYS,IOFSET)
             CALL DOORFB3(JPSYS,IPTRK,IMPX,NBMIXG,NREGAR,NUN,KEYFLX,
     1       SUNKNO(1,IG),FUNKNO(1,IG))
           ENDIF
         ENDDO
         NBMIX=NBMIXG
         NREG=NREGAR
         DEALLOCATE(SUNKN)
      ENDIF
*
      IF(IMPX.GT.3) THEN
        ALLOCATE(FGAR(NREG))
        DO IG=1,NGRP
          IF(NPSYS(IG).EQ.0) CYCLE
          CALL XDRSET(FGAR,NREG,0.0)
          DO I=1,NREG
            IF(KEYFLX(I).NE.0) FGAR(I)=FUNKNO(KEYFLX(I),IG)
          ENDDO
          WRITE(IUNOUT,'(//31H N E U T R O N    F L U X E S (,I5,
     1    3H ):)') IG
          WRITE(IUNOUT,'(1P,6(5X,E15.7))') (FGAR(I),I=1,NREG)
        ENDDO
        DEALLOCATE(FGAR)
      ENDIF
      IF(IMPX.GT.2) THEN
        CALL KDRCPU(TK2)
        WRITE(IUNOUT,'(15H SNF: CPU TIME=,1P,E11.3,8H SECOND./)')
     1  TK2-TK1
      ENDIF
      RETURN
*
  300 FORMAT(11H SNF: GROUP,I5,20H ONE-SPEED ITERATION,I4,8H  ERROR=,
     1 1P,E11.4,5H OVER,E11.4,5H PREC,E12.4,22H  ACCELERATION FACTOR=,
     2 0P,F7.3)
      END
