*DECK SNF
      SUBROUTINE SNF(KPSYS,IPTRK,IFTRAK,IMPX,NGEFF,NGIND,IDIR,NREG,
     1 NBMIX,NUN,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,TITR,NBS,KPSOU1,KPSOU2,
     2 FLUXC)
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
* KPSYS   pointer to the assembly LCM object (L_PIJ signature). KPSYS is
*         an array of directories.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  not used.
* IMPX    print flag (equal to zero for no print).
* NGEFF   number of energy groups processed in parallel.
* NGIND   energy group indices assign to the NGEFF set.
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
* NBS
* KPSOU1
* KPSOU2
*
*Parameters: input/output
* FUNKNO  unknown vector.
* FLUXC   flux at the cutoff energy.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER   TITR*72
      TYPE(C_PTR) KPSYS(NGEFF),IPTRK,KPSOU1(NGEFF),KPSOU2(NGEFF)
      INTEGER     NGEFF,NGIND(NGEFF),IFTRAK,IMPX,IDIR,NREG,NBMIX,NUN,
     1            MAT(NREG),KEYFLX(NREG),NBS(NGEFF)
      REAL        VOL(NREG),FUNKNO(NUN,NGEFF),SUNKNO(NUN,NGEFF)
      REAL,OPTIONAL :: FLUXC(NREG)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,NSTATE=40)
      INTEGER     IPAR(NSTATE)
      LOGICAL     LIVO
      DOUBLE PRECISION F1,F2,R1,R2,DMU
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: FGAR,TEST
      REAL, ALLOCATABLE, DIMENSION(:,:) :: OLD1,OLD2
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: INCONV
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
      NSTART=IPAR(20)
      MAXIT=IPAR(22)
      LIVO=(IPAR(23).EQ.1)
      ICL1=IPAR(24)
      ICL2=IPAR(25)
      IBFP=IPAR(31)
      CALL LCMGET(IPTRK,'EPSI',EPSINR)
      IF(IMPX.GT.3) THEN
        ALLOCATE(FGAR(NREG))
        DO II=1,NGEFF
          CALL XDRSET(FGAR,NREG,0.0)
          DO I=1,NREG
            IF(KEYFLX(I).NE.0) FGAR(I)=SUNKNO(KEYFLX(I),II)
          ENDDO
          WRITE(IUNOUT,'(/33H N E U T R O N    S O U R C E S (,I5,
     1    3H ):)') NGIND(II)
          WRITE(IUNOUT,'(1P,6(5X,E15.7))') (FGAR(I),I=1,NREG)
        ENDDO
        DEALLOCATE(FGAR)
      ENDIF
*
      IF(NSTART.GT.0) THEN
*----
*  GMRES(M) INNER ITERATION LOOP FOR ONE-SPEED TRANSPORT EQUATION
*----
         CALL SNGMRE (KPSYS,NGIND,IPTRK,IMPX,NGEFF,NREG,NBMIX,NUN,
     1   NSTART,MAXIT,EPSINR,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,NBS,KPSOU1,
     2   KPSOU2,FLUXC)
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
           DO I=1,NUN
             OLD1(I,II)=OLD2(I,II)
             OLD2(I,II)=FUNKNO(I,II)
           ENDDO
         ENDDO
*----
*  UPDATE THE FIXED SOURCE AND COMPUTE THE FLUX
*----
         CALL SNFLUX(KPSYS,INCONV,NGIND,IPTRK,IMPX,NGEFF,NREG,
     1   NBMIX,NUN,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,ITER,NBS,KPSOU1,
     2   KPSOU2,FLUXC)
*----
*  LOOP OVER ENERGY GROUPS
*----
         DO 60 II=1,NGEFF
         IF(INCONV(II)) THEN
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
              R2=FUNKNO(I,II)-OLD2(I,II)
              F1=F1+R1*(R2-R1)
              F2=F2+(R2-R1)*(R2-R1)
   30         CONTINUE
              DMU=-F1/F2
              IF(DMU.GT.0.0) THEN
                RDMU=REAL(DMU)
                DO 40 I=1,NUN
                FUNKNO(I,II)=OLD2(I,II)+RDMU*(FUNKNO(I,II)-OLD2(I,II))
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
           AAA=MAX(AAA,ABS(FUNKNO(KEYFLX(I),II)-OLD2(KEYFLX(I),II)))
           BBB=MAX(BBB,ABS(FUNKNO(KEYFLX(I),II)))
   50      CONTINUE
           IF(IMPX.GT.2) WRITE(IUNOUT,300) NGIND(II),ITER,AAA,BBB,
     1     AAA/BBB,DMU
           IF(IMPX.GT.5) THEN
             ALLOCATE(FGAR(NREG))
             CALL XDRSET(FGAR,NREG,0.0)
             DO I=1,NREG
              IF(KEYFLX(I).NE.0) FGAR(I)=FUNKNO(KEYFLX(I),II)
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
     1       15HATIONS IN GROUP,I5,1H.)') NGIND(II)
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
