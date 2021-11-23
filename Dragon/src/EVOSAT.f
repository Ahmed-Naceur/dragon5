*DECK EVOSAT
      SUBROUTINE EVOSAT(IMPX,MAXA,MAXB,MAXY,LOGY,NSAT,NVAR,KSAT,YST1,
     1 YSAT,MU1,IMA,NSUPF,NFISS,IDIRAC,KFISS,YSF,ADPL,BDPL,NSUPFG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Lumping of the depletion matrix, fission yields, sources and initial
* conditions to take into account the saturation of depleting nuclides.
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
* IMPX    print parameter.
* MAXA    first dimension of matrices ADPL and AGAR.
* MAXB    first dimension of matrices BDPL, IMA and MU1.
* MAXY    second dimension of matrix YSF.
* LOGY    number of passes through EVOSAT:
*         first pass updates YSAT and YST1;
*         second pass does not update YSAT and YST1.
* NSAT    number of saturating nuclides.
* NVAR    number of nuclides in the complete depletion chain.
* KSAT    position in chain of the saturating nuclides.
* NFISS   number of fissile isotopes producing fission products.
* IDIRAC  saturation model flag (=1 to use Dirac function contributions
*         in the saturating nuclide number densities).
* MU1     position of each diagonal element in vector ADPL.
* IMA     position of the first non-zero column element in vector ADPL.
* NSUPF   number of depleting fission products.
* KFISS   position in chain of the fissile isotopes.
* YSF     product of the fission yields and fission rates.
* ADPL    depletion matrix.
* BDPL    depletion source.
*
*Parameters: input/output
* YST1    number densities for all isotopes as input and of
*         the non-saturated isotopes as output.
*
*Parameters: output
* NSUPFG  number of lumped depleting fission products.
* YSAT    number densities of the saturating isotopes.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,MAXA,MAXB,MAXY,LOGY,NSAT,NVAR,KSAT(NSAT),MU1(MAXB),
     1 IMA(MAXB),NSUPF,NFISS,IDIRAC,KFISS(NFISS),NSUPFG
      REAL YST1(NVAR),YSAT(NSAT),YSF(NFISS,MAXY,LOGY),ADPL(MAXA,LOGY),
     1 BDPL(MAXB,LOGY)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(EPS=1.0E-5)
      CHARACTER HSMG*131
      LOGICAL LTEST
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEV,MGAR,IGAR
      REAL, ALLOCATABLE, DIMENSION(:) :: YSTG,AGAR,BGAR,GAR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: A22,YSFG
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: A21,A12
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(KEV(NVAR),MGAR(NVAR-NSAT),IGAR(NVAR-NSAT))
      ALLOCATE(YSTG(NVAR-NSAT),A22(NSAT,NSAT),A21(NSAT,NVAR-NSAT,LOGY),
     1 A12(NVAR-NSAT,NSAT,LOGY),AGAR(MAXA),BGAR(NVAR-NSAT),
     2 YSFG(NFISS,NSUPF),GAR(NSAT))
*
      NSUPL=NVAR-NSUPF
      I0=0
      DO 40 I=1,NVAR
      DO 10 II=1,NSAT
      IF(I.EQ.KSAT(II)) GO TO 20
   10 CONTINUE
      I0=I0+1
      KEV(I)=I0
      GO TO 40
   20 DO 25 L=1,LOGY
      IF(ADPL(MU1(I),L).EQ.0.0) CALL XABORT('EVOSAT: ZERO DIAGONAL COM'
     1 //'PONENT FOR A SATURATING ISOTOPE.')
   25 CONTINUE
      DO 30 II=1,NFISS
      IF(I.EQ.KFISS(II)) CALL XABORT('EVOSAT: A FISSILE ISOTOPE IS SAT'
     1 //'URATING.')
   30 CONTINUE
      KEV(I)=0
   40 CONTINUE
      DO 50 I=1,NFISS
      KFISS(I)=KEV(KFISS(I))
   50 CONTINUE
*----
*  FIRST LOOP OVER LOGY
*----
      DO 275 L=1,LOGY
*----
*  COMPUTE MATRICES A22**-1, A21, AND A12
*----
      DO 90 II=1,NSAT
      I=KSAT(II)
      IMAM1=0
      IF(I.GT.1) IMAM1=IMA(I-1)
      DO 60 JJ=1,NSAT
      J=KSAT(JJ)
      IF((J.LE.I).AND.(J.GT.I+IMAM1-MU1(I))) THEN
         A22(II,JJ)=ADPL(MU1(I)-I+J,L)
      ELSE IF((I.LE.J).AND.(I.GE.J-IMA(J)+MU1(J))) THEN
         A22(II,JJ)=ADPL(MU1(J)+J-I,L)
      ELSE
         A22(II,JJ)=0.0
      ENDIF
   60 CONTINUE
      JMAM1=0
      DO 75 J=1,NVAR
      J0=KEV(J)
      IF(J0.EQ.0) GO TO 70
      IF((J.LE.I).AND.(J.GT.I+IMAM1-MU1(I))) THEN
         A21(II,J0,L)=ADPL(MU1(I)-I+J,L)
      ELSE IF((I.LE.J).AND.(I.GE.J-IMA(J)+MU1(J))) THEN
         A21(II,J0,L)=ADPL(MU1(J)+J-I,L)
      ELSE
         A21(II,J0,L)=0.0
      ENDIF
      IF((I.LE.J).AND.(I.GT.J+JMAM1-MU1(J))) THEN
         A12(J0,II,L)=ADPL(MU1(J)-J+I,L)
      ELSE IF((J.LE.I).AND.(J.GE.I-IMA(I)+MU1(I))) THEN
         A12(J0,II,L)=ADPL(MU1(I)+I-J,L)
      ELSE
         A12(J0,II,L)=0.0
      ENDIF
   70 JMAM1=IMA(J)
   75 CONTINUE
      IF(I.GT.NSUPL) THEN
         DO 80 K=1,NFISS
         A21(II,KFISS(K),L)=A21(II,KFISS(K),L)+YSF(K,I-NSUPL,L)
   80    CONTINUE
      ENDIF
   90 CONTINUE
      CALL ALINV(NSAT,A22,NSAT,IER)
      IF(IER.NE.0) CALL XABORT('EVOSAT: SINGULAR MATRIX.')
*----
*  COMPUTE VECTOR YSTG ANT YSAT
*----
      IF(L.EQ.1) THEN
*        BEGINNING-OF-STAGE DIRAC DELTA CONTRIBUTIONS:
         DO 100 I=1,NSAT
         YSAT(I)=YST1(KSAT(I))
  100    CONTINUE
         DO 110 I=1,NVAR
         IF(KEV(I).GT.0) YSTG(KEV(I))=YST1(I)
  110    CONTINUE
         IF(IDIRAC.EQ.0) THEN
            DO 125 I=1,NSAT
            GAR(I)=BDPL(KSAT(I),L)
            DO 120 J=1,NVAR-NSAT
            GAR(I)=GAR(I)+A21(I,J,L)*YSTG(J)
  120       CONTINUE
  125       CONTINUE
            DO 135 I=1,NSAT
            YSAT(I)=0.0
            DO 130 J=1,NSAT
            YSAT(I)=YSAT(I)-A22(I,J)*GAR(J)
  130       CONTINUE
  135       CONTINUE
            GO TO 220
         ENDIF
         ITER=0
  140    ITER=ITER+1
         IF(ITER.GT.50) CALL XABORT('EVOSAT: CONVERGENCE FAILURE.')
         DO 155 I=1,NSAT
         GAR(I)=BDPL(KSAT(I),L)
         DO 150 J=1,NVAR-NSAT
         GAR(I)=GAR(I)+A21(I,J,L)*YSTG(J)
  150    CONTINUE
  155    CONTINUE
         ERR1=0.0
         ERR2=0.0
         DO 170 I=1,NSAT
         ZCOMP=YSAT(I)
         YSAT(I)=0.0
         DO 160 J=1,NSAT
         YSAT(I)=YSAT(I)-A22(I,J)*GAR(J)
  160    CONTINUE
         ERR1=MAX(ERR1,ABS(ZCOMP-YSAT(I)))
         ERR2=MAX(ERR2,ABS(YSAT(I)))
  170    CONTINUE
         DO 185 I=1,NSAT
         GAR(I)=0.0
         DO 180 J=1,NSAT
         GAR(I)=GAR(I)-A22(I,J)*(YST1(KSAT(J))-YSAT(J))
  180    CONTINUE
  185    CONTINUE
         DO 190 I=1,NVAR
         IF(KEV(I).GT.0) YSTG(KEV(I))=YST1(I)
  190    CONTINUE
         DO 210 I=1,NVAR-NSAT
         DO 200 J=1,NSAT
         YSTG(I)=YSTG(I)+A12(I,J,L)*GAR(J)
  200    CONTINUE
         ERR2=MAX(ERR2,ABS(YSTG(I)))
  210    CONTINUE
         IF(ERR1.LE.EPS*ERR2) GO TO 220
         GO TO 140
      ENDIF
*----
*  COMPUTE MATRICES A21 AND BGAR
*----
  220 DO 235 I=1,NSAT
      GAR(I)=0.0
      DO 230 J=1,NSAT
      GAR(I)=GAR(I)-A22(I,J)*BDPL(KSAT(J),L)
  230 CONTINUE
  235 CONTINUE
      CALL XDRSET(BGAR,NVAR-NSAT,0.0)
      DO 240 I=1,NVAR
      IF(KEV(I).GT.0) BGAR(KEV(I))=BDPL(I,L)
  240 CONTINUE
      DO 255 I=1,NVAR-NSAT
      DO 250 J=1,NSAT
      BGAR(I)=BGAR(I)+A12(I,J,L)*GAR(J)
  250 CONTINUE
  255 CONTINUE
      DO 272 J=1,NVAR-NSAT
      BDPL(J,L)=BGAR(J)
      IF(L.EQ.1) YST1(J)=YSTG(J)
      DO 260 K=1,NSAT
      GAR(K)=A21(K,J,L)
  260 CONTINUE
      DO 271 I=1,NSAT
      A21(I,J,L)=0.0
      DO 270 K=1,NSAT
      A21(I,J,L)=A21(I,J,L)+A22(I,K)*GAR(K)
  270 CONTINUE
  271 CONTINUE
  272 CONTINUE
*
  275 CONTINUE
*----
*  DETERMINE THE PROFILE PATTERN OF THE LUMPED DEPLETION MATRIX.
*----
      NSUPLG=NSUPL
      DO 280 I=1,NVAR
      IF((KEV(I).EQ.0).AND.(I.LE.NSUPL)) NSUPLG=NSUPLG-1
  280 CONTINUE
      NSUPFG=NVAR-NSAT-NSUPLG
      CALL XDISET(MGAR,NVAR-NSAT,1)
      CALL XDISET(IGAR,NVAR-NSAT,1)
      IMAM1=0
      DO 305 I=1,NVAR
      IKEV=KEV(I)
      IF(IKEV.EQ.0) GO TO 300
      DO 290 J=1,NVAR
      JKEV=KEV(J)
      IF(JKEV.EQ.0) GO TO 290
      IF((J.LE.I).AND.(J.GT.I+IMAM1-MU1(I))) THEN
         MGAR(IKEV)=MAX(MGAR(IKEV),IKEV-JKEV+1)
      ELSE IF((I.LE.J).AND.(I.GE.J-IMA(J)+MU1(J))) THEN
         IGAR(JKEV)=MAX(IGAR(JKEV),JKEV-IKEV+1)
      ENDIF
  290 CONTINUE
  300 IMAM1=IMA(I)
  305 CONTINUE
      DO 335 J=1,NVAR-NSAT
      JIFI=0
      DO 310 IFI=1,NFISS
      IF(J.EQ.KFISS(IFI)) JIFI=IFI
  310 CONTINUE
      DO 330 I=1,NVAR-NSAT
      IF((I.GT.NSUPLG).AND.(JIFI.GT.0)) GO TO 330
      LTEST=.FALSE.
      DO 325 L=1,LOGY
      DO 320 K=1,NSAT
      LTEST=LTEST.OR.(A12(I,K,L)*A21(K,J,L).NE.0.0)
  320 CONTINUE
  325 CONTINUE
      IF(LTEST.AND.(J.LE.I)) THEN
         MGAR(I)=MAX(MGAR(I),I-J+1)
      ELSE IF(LTEST) THEN
         IGAR(J)=MAX(IGAR(J),J-I+1)
      ENDIF
  330 CONTINUE
  335 CONTINUE
      II=0
      DO 340 I=1,NVAR-NSAT
      II=II+MGAR(I)
      MGAR(I)=II
      II=II+IGAR(I)-1
      IGAR(I)=II
  340 CONTINUE
      IF(IMPX.GT.8) WRITE(6,'(/27H EVOSAT: REAL SIZE OF ADPL=,I9,3H AL,
     1 13HLOCATED SIZE=,I9,1H.)') IGAR(NVAR-NSAT),MAXA
      IF(IGAR(NVAR-NSAT).GT.MAXA) THEN
         WRITE(HSMG,'(24HEVOSAT: IGAR(NVAR-NSAT)=,I6,6H MAXA=,I6)')
     1   IGAR(NVAR-NSAT),MAXA
         CALL XABORT(HSMG)
      ENDIF
*----
*  SECOND LOOP OVER LOGY
*----
      DO 540 L=1,LOGY
*----
*  COMPUTE MATRIX AGAR AND YIELDS YSFG.
*----
      CALL XDRSET(AGAR,IGAR(NVAR-NSAT),0.0)
      IMAM1=0
      DO 445 I=1,NVAR
      IKEV=KEV(I)
      IF(IKEV.EQ.0) GO TO 440
      DO 420 J=1,NVAR
      JKEV=KEV(J)
      IF(JKEV.EQ.0) GO TO 420
      IF((J.LE.I).AND.(J.GT.I+IMAM1-MU1(I))) THEN
         AGAR(MGAR(IKEV)-IKEV+JKEV)=ADPL(MU1(I)-I+J,L)
      ELSE IF((I.LE.J).AND.(I.GE.J-IMA(J)+MU1(J))) THEN
         AGAR(MGAR(JKEV)+JKEV-IKEV)=ADPL(MU1(J)+J-I,L)
      ENDIF
  420 CONTINUE
      IF(I.GT.NSUPL) THEN
         DO 430 K=1,NFISS
         YSFG(K,IKEV-NSUPLG)=YSF(K,I-NSUPL,L)
  430    CONTINUE
      ENDIF
  440 IMAM1=IMA(I)
  445 CONTINUE
      DO 495 J=1,NVAR-NSAT
      JIFI=0
      DO 450 IFI=1,NFISS
      IF(J.EQ.KFISS(IFI)) JIFI=IFI
  450 CONTINUE
      IMAM1=0
      DO 490 I=1,NVAR-NSAT
      IF((I.GT.NSUPLG).AND.(JIFI.GT.0)) GO TO 480
      IF((J.LE.I).AND.(J.GT.I+IMAM1-MGAR(I))) THEN
         DO 460 K=1,NSAT
         AGAR(MGAR(I)-I+J)=AGAR(MGAR(I)-I+J)-A12(I,K,L)*A21(K,J,L)
  460    CONTINUE
      ELSE IF((I.LE.J).AND.(I.GE.J-IGAR(J)+MGAR(J))) THEN
         DO 470 K=1,NSAT
         AGAR(MGAR(J)+J-I)=AGAR(MGAR(J)+J-I)-A12(I,K,L)*A21(K,J,L)
  470    CONTINUE
      ENDIF
  480 IMAM1=IGAR(I)
  490 CONTINUE
  495 CONTINUE
      DO 510 I=NSUPLG+1,NVAR-NSAT
      DO 505 IFI=1,NFISS
      J=KFISS(IFI)
      DO 500 K=1,NSAT
      YSFG(IFI,I-NSUPLG)=YSFG(IFI,I-NSUPLG)-A12(I,K,L)*A21(K,J,L)
  500 CONTINUE
  505 CONTINUE
  510 CONTINUE
*----
*  REPLACE THE ORIGINAL INFORMATION WITH THE LUMPED ONE
*----
      DO 520 I=1,IGAR(NVAR-NSAT)
      ADPL(I,L)=AGAR(I)
  520 CONTINUE
      DO 535 I=1,NFISS
      DO 530 J=1,NSUPFG
      YSF(I,J,L)=YSFG(I,J)
  530 CONTINUE
  535 CONTINUE
  540 CONTINUE
      DO 550 I=1,NVAR-NSAT
      IMA(I)=IGAR(I)
      MU1(I)=MGAR(I)
  550 CONTINUE
      RETURN
      END
