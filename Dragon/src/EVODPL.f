*DECK EVODPL
      SUBROUTINE EVODPL(IMPX,YDPL,NVAR,XT,EPS1,EXPMAX,H1,ITYPE,IDIRAC,
     1 IEVOL2,MU1,IMA,MAXA,NSUPF,NFISS,KFISS,YSF,ADPL,BDPL,ICHAIN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multi-purpose driver for solving the isotopic depletion equations,
* taking into account the saturation phenomena.
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
*Parameters: input/output
* IMPX    print flag (equal to zero for no print).
* YDPL    initial/final number densities.
* NVAR    number of nuclides in the complete depletion chain.
* XT      initial and final value of the independent variable.
* EPS1    required accuracy for the ODE solver.
* EXPMAX  saturation limit. A nuclide is saturating if
*         -ADPL(MU1(I))*(XT(2)-XT(1)).GT.EXPMAX. Suggested value:
*         EXPMAX=80.0. EXPMAX=0.0 means that the saturation model is
*         not used.
* H1      guessed first stepsize.
* ITYPE   type of ODE solution:
*         =1 fifth-order Runge-Kutta method;
*         =2 fourth-order Kaps-Rentrop method.
* IDIRAC  saturation model flag (=1 to use Dirac function contributions
*         in the saturating nuclide number densities.
* IEVOL2  flag making an isotope non-depleting:
*         =1 to force an isotope to be non-depleting;
*         =2 to force an isotope to be depleting;
*         =3 to force an isotope to be at saturation.
* MU1     position of each diagonal element in matrix ADPL.
* IMA     position of the first non-zero column element in matrix ADPL.
* MAXA    first dimension of matrix ADPL.
* NSUPF   number of depleting fission products.
* NFISS   number of fissile isotopes producing fission products.
* KFISS   position in chain of the fissile isotopes.
* YSF     initial/final product of the fission yields and fission
*         rates.
* ADPL    initial/final depletion matrix.
* BDPL    initial/final depletion source.
* ICHAIN  name of the isotopes in the depletion chain.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,NVAR,ITYPE,IDIRAC,IEVOL2(NVAR),MU1(NVAR),IMA(NVAR),
     1 MAXA,NSUPF,NFISS,KFISS(NFISS),ICHAIN(2,NVAR)
      REAL YDPL(NVAR,2),XT(2),EPS1,EXPMAX,H1,YSF(NFISS,NSUPF,2),
     1 ADPL(MAXA,2),BDPL(NVAR,2)
*----
*  LOCAL VARIABLES
*----
      LOGICAL LSAT
      CHARACTER*2 SHOW(120,120)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KSAT,IPERM,MU12,IMA2,KFIS2
      REAL, ALLOCATABLE, DIMENSION(:) :: YST1,YSAT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: ADPL2,BDPL2,BDPL3
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: YSF2,YSF3
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(KSAT(NVAR),IPERM(NVAR),MU12(NVAR),IMA2(NVAR),
     1 KFIS2(NFISS))
      ALLOCATE(YST1(NVAR),YSAT(NVAR),ADPL2(MAXA,2),BDPL2(NVAR,2),
     1 YSF2(NFISS,NSUPF,2))
*----
*  COMPUTE THE LUMPING INDEX VECTOR IPERM
*----
      DO 10 I=1,NVAR
      IPERM(I)=I
   10 CONTINUE
      NTER=0
   20 NTER=NTER+1
      INDSAT=0
      DO 30 I=1,NVAR
      IF(IPERM(I).GE.0) THEN
         LSAT=(IEVOL2(I).EQ.3).AND.(EXPMAX.GT.0.0)
         IF(EXPMAX.GT.0.0) THEN
            LSAT=LSAT.OR.((ABS(ADPL(MU1(I),1)*(XT(2)-XT(1))).GT.EXPMAX)
     >      .AND.(ABS(ADPL(MU1(I),2)*(XT(2)-XT(1))).GT.EXPMAX))
         ENDIF
         IF(LSAT) THEN
            DO 25 II=1,NFISS
            IF(I.EQ.KFISS(II)) GO TO 30
   25       CONTINUE
            IPERM(I)=0
            IF(INDSAT.EQ.0) THEN
               IF(IMPX.GT.5) WRITE(6,'(17H EVODPL: ISOTOPE ,2A4,
     1         18H IS SATURATING(1).)') ICHAIN(1,I),ICHAIN(2,I)
               IPERM(I)=-NTER
               INDSAT=I
            ENDIF
         ENDIF
      ENDIF
   30 CONTINUE
      IF(INDSAT.EQ.0) GO TO 60
      DO 50 I=INDSAT+1,NVAR
      JMN=I-MU1(I)+IMA(I-1)+1
      IMN=I-IMA(I)+MU1(I)
      DO 40 J=MIN(JMN,IMN),I-1
      IF((IPERM(I).EQ.0).AND.(IPERM(J).EQ.-NTER)) THEN
         IF(IMPX.GT.5) WRITE(6,'(17H EVODPL: ISOTOPE ,2A4,
     1   18H IS SATURATING(2).)') ICHAIN(1,I),ICHAIN(2,I)
         IPERM(I)=-NTER
         GO TO 50
      ENDIF
   40 CONTINUE
   50 CONTINUE
      GO TO 20
   60 NTER=NTER-1
      N=0
      DO 70 I=1,NVAR
      IF(IPERM(I).GT.0) THEN
         N=N+1
         IPERM(I)=N
      ENDIF
   70 CONTINUE
      IF(IMPX.GT.3) THEN
         WRITE(6,400) NVAR,XT(1),XT(2),EPS1,H1,ITYPE,NTER,NVAR-N,NFISS,
     1   NSUPF,(IPERM(I),I=1,NVAR)
         WRITE(6,410) (YDPL(I,1),I=1,NVAR)
      ENDIF
      IF(IMPX.GT.5) THEN
         NVARM=MIN(NVAR,120)
         WRITE (6,'(//34H EVODPL: DEPLETION MATRIX PROFILE:/)')
         DO 85 I=1,NVARM
         DO 80 J=1,NVARM
         SHOW(I,J)=' '
   80    CONTINUE
   85    CONTINUE
         IMAM1=0
         DO 120 I=1,NVARM
         DO 90 J=I-MU1(I)+IMAM1+1,I-1
         SHOW(I,J)='*'
   90    CONTINUE
         DO 100 J=I-IMA(I)+MU1(I),I-1
         SHOW(J,I)='*'
  100    CONTINUE
         IF(I.GT.NVAR-NSUPF) THEN
            DO 110 K=1,NFISS
            KFI=KFISS(K)
            IF((KFI.GT.0).AND.(KFI.LE.120)) SHOW(I,KFI)='-'
  110       CONTINUE
         ENDIF
         SHOW(I,I)='+'
         IMAM1=IMA(I)
  120    CONTINUE
         DO 130 I=1,NVARM
         WRITE (6,'(1X,I4,1X,2A4,1X,120A2)') I,ICHAIN(1,I),ICHAIN(2,I),
     1   (SHOW(I,J),J=1,NVARM)
  130    CONTINUE
         IF(NVAR.GT.120) 
     >   WRITE(6,'(34H MATRIX TRUNCATED TO 120 ELEMENTS.)')
         IF(IMPX .GE. 1000) THEN
           write(6,'(A)') 'ORIGINAL DEPLETION SYSTEM'
           write(6,'(3I10)') NVAR,NFISS,NSUPF
           write(6,'(A6)') 'MU1   '
           write(6,'(20I5)') (MU1(I),I=1,NVAR)
           write(6,'(A6)') 'IMA   '
           write(6,'(20I5)') (IMA(I),I=1,NVAR)
           write(6,'(A6)') 'ADPL1 '
           write(6,'(1P,5E20.12)') (ADPL(I,1),I=1,IMA(NVAR))
           write(6,'(A6)') 'BDPL1 '
           write(6,'(1P,5E20.12)') (BDPL(I,1),I=1,NVAR)
           write(6,'(A6)') 'KFISS '
           write(6,'(20I5)') (KFISS(K),K=1,NFISS)
           write(6,'(A6)') 'YSF1  '
           write(6,'(1P,5E20.12)') ((YSF(I,J,1),I=1,NFISS),J=1,NSUPF)
         ENDIF
      ENDIF
*----
*  LUMPING OF THE DEPLETION MATRICES
*----
      DO 135 IFI=1,NFISS
      KFIS2(IFI)=KFISS(IFI)
  135 CONTINUE
      DO 140 I=1,NVAR
      YST1(I)=YDPL(I,1)
      MU12(I)=MU1(I)
      IMA2(I)=IMA(I)
  140 CONTINUE
      DO 162 L=1,2
      DO 145 I=1,NVAR
      BDPL2(I,L)=BDPL(I,L)
  145 CONTINUE
      DO 150 I=1,IMA(NVAR)
      ADPL2(I,L)=ADPL(I,L)
  150 CONTINUE
      DO 161 I=1,NFISS
      DO 160 J=1,NSUPF
      YSF2(I,J,L)=YSF(I,J,L)
  160 CONTINUE
  161 CONTINUE
  162 CONTINUE
      NVAR2=NVAR
      NSUPF2=NSUPF
      DO 180 ITER=1,NTER
      I0=0
      NSAT=0
      DO 170 I=1,NVAR
      IF((IPERM(I).GT.0).OR.(IPERM(I).LT.-ITER)) THEN
        I0=I0+1
      ELSE IF(IPERM(I).EQ.-ITER) THEN
        I0=I0+1
        NSAT=NSAT+1
        KSAT(NSAT)=I0
      ENDIF
  170 CONTINUE
      IF(I0.NE.NVAR2) CALL XABORT('EVODPL: ALGORITHM FAILURE 1.')
      MAXB=NVAR
      MAXY=NSUPF
      CALL EVOSAT(IMPX,MAXA,MAXB,MAXY,2,NSAT,NVAR2,KSAT,YST1,YSAT,MU12,
     1 IMA2,NSUPF2,NFISS,IDIRAC,KFIS2,YSF2(1,1,1),ADPL2(1,1),BDPL2(1,1),
     2 NSUPF3)
      NVAR2=NVAR2-NSAT
      NSUPF2=NSUPF3
      NSAT=0
      I0=0
      DO 175 I=1,NVAR
      IF((IPERM(I).GT.0).OR.(IPERM(I).LT.-ITER)) THEN
         I0=I0+1
         YDPL(I,1)=YST1(I0)
      ELSE IF(IPERM(I).EQ.-ITER) THEN
         NSAT=NSAT+1
         YDPL(I,1)=YSAT(NSAT)
      ENDIF
  175 CONTINUE
  180 CONTINUE
      IF(IMPX.GT.4) WRITE(6,420) (YDPL(I,1),I=1,NVAR)
*----
*  SOLUTION OF THE LUMPED DEPLETION SYSTEM
*----
      DO 185 I=1,NVAR
      YDPL(I,2)=YDPL(I,1)
  185 CONTINUE
      IF(NVAR2.EQ.0) GO TO 315
      DO 190 I=1,NVAR2
      FACT=(BDPL2(I,2)-BDPL2(I,1))/(XT(2)-XT(1))
      BDPL2(I,1)=BDPL2(I,1)-FACT*XT(1)
      BDPL2(I,2)=FACT
  190 CONTINUE
      DO 200 I=1,IMA2(NVAR2)
      FACT=(ADPL2(I,2)-ADPL2(I,1))/(XT(2)-XT(1))
      ADPL2(I,1)=ADPL2(I,1)-FACT*XT(1)
      ADPL2(I,2)=FACT
  200 CONTINUE
      DO 215 I=1,NFISS
      DO 210 J=1,NSUPF2
      FACT=(YSF2(I,J,2)-YSF2(I,J,1))/(XT(2)-XT(1))
      YSF2(I,J,1)=YSF2(I,J,1)-FACT*XT(1)
      YSF2(I,J,2)=FACT
  210 CONTINUE
  215 CONTINUE
      IF(IMPX.GT.4) THEN
         WRITE(6,430) NSUPF2
         WRITE(6,440) (YST1(I),I=1,NVAR2)
      ENDIF
      IF(IMPX.GT.5) THEN
         NVARM=MIN(NVAR2,120)
         WRITE (6,'(//41H EVODPL: LUMPED DEPLETION MATRIX PROFILE:/)')
         DO 225 I=1,NVARM
         DO 220 J=1,NVARM
         SHOW(I,J)=' '
  220    CONTINUE
  225    CONTINUE
         IMAM1=0
         DO 260 I=1,NVARM
         DO 230 J=I-MU12(I)+IMAM1+1,I-1
         SHOW(I,J)='*'
  230    CONTINUE
         DO 240 J=I-IMA2(I)+MU12(I),I-1
         SHOW(J,I)='*'
  240    CONTINUE
         IF(I.GT.NVAR2-NSUPF2) THEN
            DO 250 K=1,NFISS
            KFI=KFIS2(K)
            IF((KFI.GT.0).AND.(KFI.LE.60)) SHOW(I,KFI)='-'
  250       CONTINUE
         ENDIF
         SHOW(I,I)='+'
         IMAM1=IMA2(I)
  260    CONTINUE
         DO 270 I=1,NVARM
         WRITE (6,'(1X,I4,1X,2A4,1X,120A2)') I,ICHAIN(1,I),ICHAIN(2,I),
     1   (SHOW(I,J),J=1,NVARM)
  270    CONTINUE
         IF(NVAR.GT.120) 
     >   WRITE(6,'(34H MATRIX TRUNCATED TO 120 ELEMENTS.)')
         IF(IMPX .GE. 1000) THEN
           write(6,'(A)') 'LUMPED DEPLETION SYSTEM'
           write(6,'(3I10)') NVAR2,NFISS,NSUPF2
           write(6,'(A6)') 'MU1   '
           write(6,'(20I5)') (MU12(I),I=1,NVAR2)
           write(6,'(A6)') 'IMA   '
           write(6,'(20I5)') (IMA2(I),I=1,NVAR2)
           write(6,'(A6)') 'ADPL2 '
           write(6,'(1P,5E20.12)') (ADPL2(I,1),I=1,IMA2(NVAR2))
           write(6,'(A6)') 'BDPL2 '
           write(6,'(1P,5E20.12)') (BDPL2(I,1),I=1,NVAR2)
           write(6,'(A6)') 'KFISS '
           write(6,'(20I5)') (KFIS2(K),K=1,NFISS)
           write(6,'(A6)') 'YSF1  '
           write(6,'(1P,5E20.12)') ((YSF2(I,J,1),I=1,NFISS),J=1,NSUPF2)
         ENDIF
      ENDIF
      ALLOCATE(BDPL3(NVAR2,2),YSF3(NFISS,NSUPF2,2))
      DO 280 I=1,NVAR2
      BDPL3(I,1)=BDPL2(I,1)
      BDPL3(I,2)=BDPL2(I,2)
  280 CONTINUE
      DO 295 I=1,NFISS
      DO 290 J=1,NSUPF2
      YSF3(I,J,1)=YSF2(I,J,1)
      YSF3(I,J,2)=YSF2(I,J,2)
  290 CONTINUE
  295 CONTINUE
      CALL EVOODE(YST1,NVAR2,XT(1),XT(2),EPS1,H1,NOK,NBAD,ITYPE,MU12,
     1 IMA2,MAXA,NSUPF2,NFISS,KFIS2,YSF3,ADPL2,BDPL3)
      DEALLOCATE(YSF3,BDPL3)
      IF(IMPX.GT.4) THEN
         WRITE(6,450) (YST1(I),I=1,NVAR2)
         IF(ITYPE.LE.2) WRITE(6,'(13H EVODPL: NOK=,I5,6H NBAD=,I5)')
     1   NOK,NBAD
      ENDIF
      DO 310 I=1,NVAR
      IF(IPERM(I).GT.0) YDPL(I,2)=YST1(IPERM(I))
  310 CONTINUE
*----
*  COMPUTE NUMBER DENSITIES OF THE SATURATED ISOTOPES
*----
  315 IF(NTER.EQ.0) GO TO 370
      DO 320 I=1,NVAR
      YST1(I)=YDPL(I,2)
      BDPL2(I,2)=BDPL(I,2)
      MU12(I)=MU1(I)
      IMA2(I)=IMA(I)
  320 CONTINUE
      DO 330 I=1,IMA(NVAR)
      ADPL2(I,2)=ADPL(I,2)
  330 CONTINUE
      DO 345 I=1,NFISS
      KFIS2(I)=KFISS(I)
      DO 340 J=1,NSUPF
      YSF2(I,J,2)=YSF(I,J,2)
  340 CONTINUE
  345 CONTINUE
      NVAR2=NVAR
      NSUPF2=NSUPF
      DO 365 ITER=1,NTER
      I0=0
      NSAT=0
      DO 350 I=1,NVAR
      IF((IPERM(I).GT.0).OR.(IPERM(I).LT.-ITER)) THEN
        I0=I0+1
      ELSE IF(IPERM(I).EQ.-ITER) THEN
        I0=I0+1
        NSAT=NSAT+1
        KSAT(NSAT)=I0
      ENDIF
  350 CONTINUE
      IF(I0.NE.NVAR2) CALL XABORT('EVODPL: ALGORITHM FAILURE 2.')
      MAXB=NVAR
      MAXY=NSUPF
      CALL EVOSAT(IMPX,MAXA,MAXB,MAXY,1,NSAT,NVAR2,KSAT,YST1,YSAT,MU12,
     1 IMA2,NSUPF2,NFISS,IDIRAC,KFIS2,YSF2(1,1,2),ADPL2(1,2),BDPL2(1,2),
     2 NSUPF3)
      IF(IMPX.GT.4) WRITE(6,425) ITER,(YSAT(I),I=1,NSAT)
      NVAR2=NVAR2-NSAT
      NSUPF2=NSUPF3
      NSAT=0
      I0=0
      DO 360 I=1,NVAR
      IF((IPERM(I).GT.0).OR.(IPERM(I).LT.-ITER)) THEN
        I0=I0+1
        YDPL(I,2)=YST1(I0)
      ELSE IF(IPERM(I).EQ.-ITER) THEN
        NSAT=NSAT+1
        YDPL(I,2)=YSAT(NSAT)
      ENDIF
  360 CONTINUE
  365 CONTINUE
  370 IF(IMPX.GT.3) WRITE(6,460) (YDPL(I,2),I=1,NVAR)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(YSF2,BDPL2,ADPL2,YSAT,YST1)
      DEALLOCATE(KFIS2,IMA2,MU12,IPERM,KSAT)
      RETURN
*
  400 FORMAT(//45H EVODPL: SOLUTION OF THE DEPLETION EQUATIONS.//14X,
     1 25HTOTAL NUMBER OF NUCLIDES=,I5/26X,13HINITIAL TIME=,1P,E12.4/
     2 28X,11HFINAL TIME=,E12.4/15X,24HACCURACY FOR ODE SOLVER=,E12.4/
     3 16X,23HGUESSED FIRST STEPSIZE=,E12.4,0P/22X,17HTYPE OF SOLUTION=,
     4 I3/39H NUMBER OF GROUP OF SATURATED NUCLIDES=,I5/10X,
     5 29HNUMBER OF SATURATED NUCLIDES=,I5/12X,19HNUMBER OF FISSILE N,
     6 8HUCLIDES=,I5/12X,27HNUMBER OF FISSION PRODUCTS=,I5//
     7 22H LUMPING INDEX VECTOR:/(1X,20I5))
  410 FORMAT(/48H EVODPL: INITIAL VALUES OF THE DEPLETION SYSTEM:/
     1 (1X,1P,10E12.4))
  420 FORMAT(/53H EVODPL: SATURATED INITIAL CONDITIONS OF THE DEPLETIO,
     1 9HN SYSTEM:/(1X,1P,10E12.4))
  425 FORMAT(/51H EVODPL: FINAL VALUES OF THE SATURATED NUCLIDES IN ,
     1 9HGROUP NO.,I5//(1X,1P,10E12.4))
  430 FORMAT(/42H NUMBER OF NON-SATURATED FISSION PRODUCTS=,I5)
  440 FORMAT(/55H EVODPL: INITIAL VALUES OF THE LUMPED DEPLETION SYSTEM:
     1 /(1X,1P,10E12.4))
  450 FORMAT(/53H EVODPL: ODE SOLUTION OF THE LUMPED DEPLETION SYSTEM:/
     1 (1X,1P,10E12.4))
  460 FORMAT(/42H EVODPL: SOLUTION OF THE DEPLETION SYSTEM:/
     1 (1X,1P,10E12.4))
      END
