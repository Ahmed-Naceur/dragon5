*DECK EVOMU1
      SUBROUTINE EVOMU1(IMPX,NVAR,NREAC,LP,XT,LCOOL,NPAR,KPAR,DCR,SIG1,
     1 SIG2,EXPMAX,IEVOLB,MU1,IMA,MAXA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Determination of the profile of the depletion matrix (not taking into
* account the fission yields).
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
* IMPX    print flag (equal to zero for no print).
* NVAR    number of nuclides in the complete depletion chain.
* NREAC   one plus the number of neutron-induced depletion reactions.
* LP      lumping index (set to zero to get rid of unused isotopes).
* XT      initial and final value of the independent variable.
* LCOOL   out-of-core depletion flag (LCOOL=.true. to set flag).
* NPAR    maximum number of parent nuclides in the depletion chain.
* KPAR    position in chain of the parent nuclide and type of
*         reaction.
* DCR     sum of radioactive decay constants in 10**-8/s.
* SIG1    initial reaction rates:
*         SIG1(I,1) fission reaction rate for nuclide I;
*         SIG1(I,2) gamma reaction rate for nuclide I;
*         SIG1(I,3) N2N reaction rate for nuclide I;
*         ...;
*         SIG1(I,NREAC) neutron-induced energy released for nuclide I;
*         SIG1(I,NREAC+1) decay energy released for nuclide I.
* SIG2    final reaction rates.
* EXPMAX  saturation limit. A nuclide is saturating if
*         -ADPL(MU1(I))*(XT(2)-XT(1)).GT.EXPMAX. Suggested value:
*         EXPMAX=80.0.
* IEVOLB  flag making an isotope non-depleting:
*         =1 to force an isotope to be non-depleting;
*         =2 to force an isotope to be depleting;
*         =3 to force an isotope to be at saturation.
*
*Parameters: output
* MU1     position of each diagonal element in vector ADPL.
* IMA     position of the first non-zero column element in vector ADPL.
* MAXA    maximum number of terms in ADPL, taking into account
*         saturation.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      LOGICAL LCOOL
      INTEGER IMPX,NVAR,NREAC,LP(NVAR),NPAR,KPAR(NPAR,NVAR),
     1 IEVOLB(NVAR),MU1(NVAR+1),IMA(NVAR+1),MAXA
      REAL XT(2),DCR(NVAR),SIG1(NVAR+1,NREAC+1),SIG2(NVAR+1,NREAC+1),
     1 EXPMAX
*----
*  LOCAL VARIABLES
*----
      LOGICAL LSAT
      CHARACTER*2 SHOW(65,65)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPERM
      REAL, ALLOCATABLE, DIMENSION(:) :: DIAG
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IPERM(NVAR),DIAG(NVAR))
*
      DO 10 IS=1,NVAR
      IF(LP(IS).EQ.0) GO TO 10
      MU1(LP(IS))=1
      IMA(LP(IS))=1
   10 CONTINUE
      NVAR2=0
      DO 30 IS=1,NVAR
      IF(LP(IS).EQ.0) GO TO 30
      NVAR2=MAX(NVAR2,LP(IS))
      DO 20 K=1,NPAR
      IF(KPAR(K,IS).EQ.0) GO TO 30
      IF(KPAR(K,IS).EQ.2) GO TO 20
      JS=KPAR(K,IS)/100
      KT=KPAR(K,IS)-JS*100
      IF((LCOOL.AND.(KT.GE.2)).OR.(LP(JS).EQ.0)) GO TO 20
      MU1(LP(IS))=MAX(MU1(LP(IS)),LP(IS)-LP(JS)+1)
      IMA(LP(JS))=MAX(IMA(LP(JS)),LP(JS)-LP(IS)+1)
   20 CONTINUE
   30 CONTINUE
*
      DO 40 IS=1,NVAR
      DIAG(IS)=0.0
      IF(LP(IS).EQ.0) GO TO 40
      SIGE1=0.0
      SIGE2=0.0
      IF(.NOT.LCOOL) THEN
         DO 35 IREAC=1,NREAC-1
         SIGE1=SIGE1+SIG1(IS,IREAC)
         SIGE2=SIGE2+SIG2(IS,IREAC)
   35    CONTINUE
      ENDIF
      DIAG(IS)=(MIN(SIGE1,SIGE2)+DCR(IS))*(XT(2)-XT(1))
   40 CONTINUE
*----
*  COMPUTE THE LUMPING INDEX VECTOR IPERM
*----
      DO 50 I=1,NVAR2
      IPERM(I)=I
   50 CONTINUE
      NTER=0
   60 NTER=NTER+1
      INDSAT=0
      DO 70 IS=1,NVAR
      IF(LP(IS).GT.0) THEN
         LSAT=(IEVOLB(IS).EQ.3).AND.(EXPMAX.GT.0.0)
         LSAT=LSAT.OR.((EXPMAX.GT.0.0).AND.(DIAG(IS).GT.EXPMAX))
         IF((IPERM(LP(IS)).GE.0).AND.LSAT) THEN
            IPERM(LP(IS))=0
            IF(INDSAT.EQ.0) THEN
               IPERM(LP(IS))=-NTER
               INDSAT=LP(IS)
            ENDIF
         ENDIF
      ENDIF
   70 CONTINUE
      IF(INDSAT.EQ.0) GO TO 100
      DO 90 I=INDSAT+1,NVAR2
      DO 80 J=MIN(I-MU1(I)+1,I-IMA(I)+1),I-1
      IF((IPERM(I).EQ.0).AND.(IPERM(J).EQ.-NTER)) THEN
         IPERM(I)=-NTER
         GO TO 90
      ENDIF
   80 CONTINUE
   90 CONTINUE
      GO TO 60
  100 NTER=NTER-1
      N=0
      DO 110 I=1,NVAR2
      IF(IPERM(I).GT.0) THEN
         N=N+1
         IPERM(I)=N
      ENDIF
  110 CONTINUE
*
      MAXA=0
      DO 175 ITER=1,NTER
      JMIN=NVAR2
      JMAX=1
      DO 130 I=1,NVAR2
      IF(-IPERM(I).NE.ITER) GO TO 130
      DO 120 J=1,NVAR2
      IF((J.LE.I).AND.(J.GT.I-MU1(I))) THEN
         JMIN=MIN(JMIN,J)
         JMAX=MAX(JMAX,J)
      ELSE IF((I.LE.J).AND.(I.GT.J-IMA(J))) THEN
         JMIN=MIN(JMIN,J)
         JMAX=MAX(JMAX,J)
      ENDIF
  120 CONTINUE
  130 CONTINUE
      IMIN=NVAR2
      IMAX=1
      DO 145 I=1,NVAR2
      DO 140 J=1,NVAR2
      IF(-IPERM(J).NE.ITER) GO TO 140
      IF((J.LE.I).AND.(J.GT.I-MU1(I))) THEN
         IMIN=MIN(IMIN,I)
         IMAX=MAX(IMAX,I)
      ELSE IF((I.LE.J).AND.(I.GT.J-IMA(J))) THEN
         IMIN=MIN(IMIN,I)
         IMAX=MAX(IMAX,I)
      ENDIF
  140 CONTINUE
  145 CONTINUE
      DO 170 I=IMIN,IMAX
      IF(-IPERM(I).EQ.ITER) GO TO 170
      DO 160 J=JMIN,JMAX
      IF(-IPERM(J).EQ.ITER) GO TO 160
      IF((J.LE.I).AND.(J.GT.I-MU1(I))) GO TO 160
      IF((I.LE.J).AND.(I.GE.J-IMA(J))) GO TO 160
      MAXA=MAXA+1
  160 CONTINUE
  170 CONTINUE
  175 CONTINUE
*
      MAXROW=1
      MAXCOL=1
      II=0
      DO 180 I=1,NVAR2
      MAXROW=MAX(MAXROW,MU1(I))
      MAXCOL=MAX(MAXCOL,IMA(I))
      II=II+MU1(I)
      MU1(I)=II
      II=II+IMA(I)-1
      IMA(I)=II
  180 CONTINUE
      MAXA=MAXA+IMA(NVAR2)
      IF(IMPX.GT.3) WRITE (6,'(/34H EVOMU1: MAXIMUM ROW PROFILE WIDTH,
     1 4X,1H=,I5/9X,30HMAXIMUM COLUMN PROFILE WIDTH =,I5)') MAXROW,
     2 MAXCOL
      IF(IMPX.GT.9) THEN
         NVARM=MIN(NVAR2,65)
         WRITE (6,'(//34H EVOMU1: DEPLETION MATRIX PROFILE:/)')
         DO 195 I=1,NVARM
         DO 190 J=1,NVARM
         SHOW(I,J)=' '
  190    CONTINUE
  195    CONTINUE
         IMAM1=0
         DO 220 I=1,NVARM
         DO 200 J=I-MU1(I)+IMAM1+1,I
         IF(IPERM(I).GT.0) THEN
            SHOW(I,J)='*'
         ELSE
            SHOW(I,J)='+'
         ENDIF
  200    CONTINUE
         DO 210 J=I-IMA(I)+MU1(I),I
         IF(IPERM(J).GT.0) THEN
            SHOW(J,I)='*'
         ELSE
            SHOW(J,I)='+'
         ENDIF
  210    CONTINUE
         IMAM1=IMA(I)
  220    CONTINUE
         DO 230 I=1,NVARM
         WRITE (6,'(1X,65A2)') (SHOW(I,J),J=1,NVARM)
  230    CONTINUE
         IF(NVAR2.GT.65) WRITE(6,'(18H MATRIX TRUNCATED.)')
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DIAG,IPERM)
      RETURN
      END
