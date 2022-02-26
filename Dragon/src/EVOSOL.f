*DECK EVOSOL
      SUBROUTINE EVOSOL(IMPX,LCOOL,NVAR,NREAC,NDFP,NPAR,NFISS,XT,EPS1,
     1 EXPMAX,H1,ITYPE,IDIRAC,DCR,KPAR,BPAR,KFISS,KPF,YIELD,LP,IEVOLB,
     2 SIG1,SIG2,NVAR2,NFISS2,NSUPF2,MU1,IMA,MAXA,YDPL,ICHAIN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Put the depletion matrix system in sparse storage mode for a single
* depleting mixture. Solve this system between times XT(1) and XT(2).
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
* LCOOL   out-of-core depletion flag (LCOOL=.true. to set flag).
* NVAR    number of depleting nuclides.
* NREAC   one plus the number of neutron-induced depletion reactions.
* NDFP    number of direct fission products (fission fragments).
* NPAR    maximum number of parent nuclides in the depletion chain.
* NFISS   number of fissile isotopes producing fission products.
* XT      initial and final time (independent variable).
* EPS1    required accuracy for the ODE solver.
* EXPMAX  saturation limit. A nuclide is saturating if
*         -ADPL(MU1(I))*(XT(2)-XT(1)).GT.EXPMAX. Suggested value:
*         EXPMAX=80.0.
* H1      guessed first stepsize.
* ITYPE   type of ODE solution:
*         =1 fifth-order Runge-Kutta method;
*         =2 fourth-order Kaps-Rentrop method.
* IDIRAC  saturation model flag (=1 to use Dirac function contributions
*         in the saturating nuclide number densities).
* DCR     sum of radioactive decay constants in 10**-8/s
* KPAR    position in chain of the parent nuclide and type of
*         reaction.
* BPAR    branching ratio for neutron induced reactions.
* KFISS   position in chain of the fissile isotopes.
* KPF     position in chain of the direct fission products (fission
*         fragments).
* YIELD   fission yields.
* LP      index vector used to remove unused isotopes from the
*         depletion chain.
* IEVOLB  flag making an isotope non-depleting:
*         =1 to force an isotope to be non-depleting;
*         =2 to force an isotope to be depleting;
*         =3 to force an isotope to be at saturation.
* SIG1    initial reaction rates for nuclide I:
*         SIG1(I,1) fission reaction rate;
*         SIG1(I,2) gamma reaction rate;
*         SIG1(I,3) N2N reaction rate;
*         ...;
*         SIG1(I,NREAC) neutron-induced energy;
*         SIG1(I,NREAC+1) decay energy released.
* SIG2    final reaction rates.
* NVAR2   number of used isotopes in the depletion chain, where
*         NVAR2=max(LP(I)).
* NFISS2  number of used fissile isotopes producing fission products.
* NSUPF2  number of used fission products.
* MU1     position of each diagonal element in matrix ADPL.
* IMA     position of the first non-zero column element in matrix ADPL.
* MAXA    first dimension of matrix ADPL.
* YDPL    initial/final number density of each depleting isotope.
*         YDPL(NVAR+1,2) is the stage burnup increment.
* ICHAIN  name of the used isotopes in the depletion chain.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      LOGICAL LCOOL
      INTEGER IMPX,NVAR,NREAC,NDFP,NPAR,NFISS,ITYPE,IDIRAC,
     1 KPAR(NPAR,NVAR),KFISS(NFISS),KPF(NDFP),LP(NVAR),IEVOLB(NVAR),
     2 NVAR2,NFISS2,NSUPF2,MU1(NVAR2+1),IMA(NVAR2+1),MAXA,
     3 ICHAIN(2,NVAR2+1)
      REAL XT(2),EPS1,EXPMAX,H1,DCR(NVAR),BPAR(NPAR,NVAR),
     1 YIELD(NFISS,NDFP),SIG1(NVAR+1,NREAC+1),SIG2(NVAR+1,NREAC+1),
     2 YDPL(NVAR+1,2)
*----
*  LOCAL VARIABLES
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LQ,KFISS2,IEVOL2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: ADPL,BDPL,YDPL2
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: YSF
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(LQ(NFISS),KFISS2(NFISS2),IEVOL2(NVAR2+1))
      ALLOCATE(ADPL(MAXA,2),BDPL(NVAR2+1,2),YSF(NFISS2,NSUPF2+1,2),
     1 YDPL2(NVAR2+1,2))
*----
*  COMPUTE LQ AND KFISS2
*----
      I0=0
      DO 10 I=1,NFISS
      LQ(I)=0
      IF(KFISS(I).EQ.0) GO TO 10
      IF((LP(KFISS(I)).EQ.0).OR.(LP(KFISS(I)).GT.NVAR2)) GO TO 10
      I0=I0+1
      IF(I0.GT.NFISS2) CALL XABORT('EVOSOL: NFISS2 TOO SMALL.')
      KFISS2(I0)=LP(KFISS(I))
      LQ(I)=I0
   10 CONTINUE
*----
*  BUILD THE SPARSE DEPLETION MATRICES ADPL, BDPL AND YSF
*----
      DO 80 IP=1,2
      DO 20 I=1,IMA(NVAR2+1)
      ADPL(I,IP)=0.0
   20 CONTINUE
      BDPL(NVAR2+1,IP)=0.0
      DO 40 IS=1,NVAR
      IF((LP(IS).EQ.0).OR.(LP(IS).GT.NVAR2)) GO TO 40
      BDPL(LP(IS),IP)=0.0
      SIGE=0.0
      IF((.NOT.LCOOL).AND.(IP.EQ.1)) THEN
         DO 25 IREAC=1,NREAC-1
         SIGE=SIGE+SIG1(IS,IREAC)
   25    CONTINUE
      ELSE IF((.NOT.LCOOL).AND.(IP.EQ.2)) THEN
         DO 26 IREAC=1,NREAC-1
         SIGE=SIGE+SIG2(IS,IREAC)
   26    CONTINUE
      ENDIF
      ADPL(MU1(LP(IS)),IP)=-SIGE-DCR(IS)
      DO 30 IPAR=1,NPAR
      IF(KPAR(IPAR,IS).EQ.0) GO TO 40
      IF(KPAR(IPAR,IS).EQ.2) GO TO 30
      JS=KPAR(IPAR,IS)/100
      KT=KPAR(IPAR,IS)-JS*100
      IF((LCOOL.AND.(KT.GE.2)).OR.(LP(JS).EQ.0).OR.(LP(JS).GT.NVAR2))
     1 GO TO 30
      SIGE=0.0
      IF(KT.EQ.1) THEN
         SIGE=BPAR(IPAR,IS)*DCR(JS)
      ELSE IF((KT.GE.3).AND.(IP.EQ.1)) THEN
         SIGE=BPAR(IPAR,IS)*SIG1(JS,KT-1)
      ELSE IF((KT.GE.3).AND.(IP.EQ.2)) THEN
         SIGE=BPAR(IPAR,IS)*SIG2(JS,KT-1)
      ELSE
         CALL XABORT('EVOSOL: UNKNOWN REACTION.')
      ENDIF
      IF(LP(JS).LE.LP(IS)) THEN
         ADPL(MU1(LP(IS))-LP(IS)+LP(JS),IP)=
     1   ADPL(MU1(LP(IS))-LP(IS)+LP(JS),IP)+SIGE
      ELSE
         ADPL(MU1(LP(JS))+LP(JS)-LP(IS),IP)=
     1   ADPL(MU1(LP(JS))+LP(JS)-LP(IS),IP)+SIGE
      ENDIF
   30 CONTINUE
   40 CONTINUE
      IF(LCOOL) THEN
         DO 51 IFIS=1,NFISS2
         DO 50 ISUPF=1,NSUPF2+1
         YSF(IFIS,ISUPF,IP)=0.0
   50    CONTINUE
   51    CONTINUE
      ELSE
*        ADD ONE EQUATION TO COMPUTE THE BURNUP.
         DO 55 JS=1,NVAR
         IF((LP(JS).EQ.0).OR.(LP(JS).GT.NVAR2)) GO TO 55
         IF(IP.EQ.1) THEN
            ADPL(MU1(NVAR2+1)-(NVAR2+1)+LP(JS),IP)=SIG1(JS,NREAC)+
     &      SIG1(JS,NREAC+1)
         ELSE IF(IP.EQ.2) THEN
            ADPL(MU1(NVAR2+1)-(NVAR2+1)+LP(JS),IP)=SIG2(JS,NREAC)+
     &      SIG2(JS,NREAC+1)
         ENDIF
   55    CONTINUE
*
*        ADD THE FISSION YIELD CONTRIBUTIONS.
         DO 70 IFIS=1,NFISS
         IF(LQ(IFIS).EQ.0) GO TO 70
         DO 56 ISUPF=1,NSUPF2+1
         YSF(LQ(IFIS),ISUPF,IP)=0.0
   56    CONTINUE
         IF(NSUPF2.EQ.0) GO TO 70
         DO 60 ISUPF=1,NDFP
         LPP=LP(KPF(ISUPF))
         IF((LPP.EQ.0).OR.(LPP.GT.NVAR2)) GO TO 60
         IF(LPP+NSUPF2.LE.NVAR2) CALL XABORT('EVOSOL: FAILURE.')
         IF(IP.EQ.1) THEN
            YSF(LQ(IFIS),LPP-NVAR2+NSUPF2,IP)=YIELD(IFIS,ISUPF)*
     1      SIG1(KFISS(IFIS),1)
         ELSE IF(IP.EQ.2) THEN
            YSF(LQ(IFIS),LPP-NVAR2+NSUPF2,IP)=YIELD(IFIS,ISUPF)*
     1      SIG2(KFISS(IFIS),1)
         ENDIF
   60    CONTINUE
   70    CONTINUE
      ENDIF
   80 CONTINUE
*----
*  SOLVE THE DEPLETION SYSTEM. EQUATION NVAR2+1 IS USED TO COMPUTE
*  THE BURNUP
*----
      CALL XDRSET(YDPL2(1,1),NVAR2+1,0.0)
      CALL XDISET(IEVOL2,NVAR2+1,0)
      DO 90 IS=1,NVAR
      IF((LP(IS).EQ.0).OR.(LP(IS).GT.NVAR2)) GO TO 90
      YDPL2(LP(IS),1)=YDPL(IS,1)
      IEVOL2(LP(IS))=IEVOLB(IS)
   90 CONTINUE
      CALL EVODPL(IMPX,YDPL2(1,1),NVAR2+1,XT,EPS1,EXPMAX,H1,ITYPE,
     1 IDIRAC,IEVOL2,MU1,IMA,MAXA,NSUPF2+1,NFISS2,KFISS2,YSF,ADPL,
     2 BDPL,ICHAIN)
      YDPL(NVAR+1,2)=YDPL2(NVAR2+1,2)
      DO 100 IS=NVAR,1,-1
      IF((LP(IS).EQ.0).OR.(LP(IS).GT.NVAR2)) THEN
         YDPL(IS,2)=YDPL(IS,1)
      ELSE
         YDPL(IS,1)=YDPL2(LP(IS),1)
         YDPL(IS,2)=MAX(0.0,YDPL2(LP(IS),2))
      ENDIF
  100 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(YDPL2,YSF,BDPL,ADPL)
      DEALLOCATE(IEVOL2,KFISS2,LQ)
      RETURN
      END
