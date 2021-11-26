*DECK LNSR
      SUBROUTINE LNSR(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* single iteration for the line optimization of the objective function.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* NENTRY  number of data structures transfered to this module.
* HENTRY  name of the data structures.
* IENTRY  data structure type where:
*         IENTRY=1 for LCM memory object;
*         IENTRY=2 for XSM file;
*         IENTRY=3 for sequential binary file;
*         IENTRY=4 for sequential ASCII file.
* JENTRY  access permission for the data structure where:
*         JENTRY=0 for a data structure in creation mode;
*         JENTRY=1 for a data structure in modifications mode;
*         JENTRY=2 for a data structure in read-only mode.
* KENTRY  data structure pointer.
*
*Comments:
* The calling specifications are:
* OPTIM := LNSR: OPTIM :: (lnsr\_data) ;
* where
*   OPTIM : name of the \emph{optimize} object (L\_OPTIMIZE signature) 
*     containing the optimization informations. Object OPTIM must appear on 
*     both LHS and RHS to be able to update the previous values.
*   (lnsr\_data) : structure containing the data to the module LNSR.
*
*Reference:
*  L. Armijo, "Minimization of functions having Lipschitz continuous
*  first partial derivatives," Pacific journal of mathematics, Vol. 16,
*  No. 1, 1-3, 1966.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,MAXINT=30)
      TYPE(C_PTR) IPGRAD
      CHARACTER TEXT12*12,HSIGN*12
      INTEGER ISTATE(NSTATE),CNVTST,DNVTST,hist_nr
      DOUBLE PRECISION OPTPRR(NSTATE)
      REAL FLOTT
      DOUBLE PRECISION DFLOTT,SR,DSAVE(3)
      PARAMETER(XI=0.5D0,WIDTH=0.5D0) ! Armijo parameters
*----
*  ALLOCATABLE ARRAYS
*----
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: XP,GP,V,Y,YGG,GGY,
     1 FF,UD,GAMMA
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: X,P,G,XMIN,XMAX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: AA,GG,DFF,TDFF,
     1 SS,YY
*----
*  PARAMETER VALIDATION.
*----
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('LNSR: LCM'
     1 //' OBJECT EXPECTED AT LHS.')
      IF(JENTRY(1).NE.1) CALL XABORT('LNSR: OBJECT IN MODIFICATION MOD'
     1  //'E EXPECTED.')
      IPGRAD=KENTRY(1)
      CALL LCMGTC(IPGRAD,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_OPTIMIZE') THEN
         TEXT12=HENTRY(1)
         CALL XABORT('LNSR: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_OPTIMIZE EXPECTED.')
      ENDIF
      CNVTST=-1
      ICONV =0
      DNVTST=-1
*----
*  READ INPUT PARAMETERS
*----
      CALL LCMGET(IPGRAD,'STATE-VECTOR',ISTATE)
      IF((ISTATE(2).NE.0).AND.(ISTATE(8).NE.4)) THEN
        CALL XABORT('LNSR: CONSTRAINTS NOT IMPLEMENTED.')
      ENDIF
      NVAR  =ISTATE(1)
      NFUNC =ISTATE(2)+1
      IOPT  =ISTATE(3)
      ICONV =ISTATE(4)
      IF((IOPT.NE.1).AND.(IOPT.NE.-1)) CALL XABORT('LNSR: IOPT not equ'
     1  //'al to 1 or -1')
      IEXT  =ISTATE(5)
      IF(IEXT.EQ.0) IEXT=1
      IEDSTP=ISTATE(6)
      IHESS =ISTATE(7)
      IMETH =ISTATE(8)
      ISTEP =ISTATE(10)
      JCONV =ISTATE(11)
      MAXEXT=ISTATE(12)
      NSTART=ISTATE(13)
      CALL LCMGET(IPGRAD,'OPT-PARAM-R',OPTPRR)
      SR=OPTPRR(1)
      EPS1=OPTPRR(2)
      EPS2=OPTPRR(3)
      EPS3=OPTPRR(4)
      IPICK=0
      hist_nr=10
      IPRINT=1
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF(INDIC.EQ.10) GO TO 20
      IF(INDIC.NE.3) CALL XABORT('LNSR: CHARACTER DATA EXPECTED(1).')
   15 IF(TEXT12.EQ.'EDIT') THEN
        CALL REDGET(INDIC,IPRINT,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('LNSR: INTEGER DATA EXPECTED(1).')
      ELSE IF(TEXT12.EQ.'MINIMIZE') THEN
        IOPT=1
      ELSE IF(TEXT12.EQ.'MAXIMIZE') THEN
        IOPT=-1
      ELSE IF(TEXT12.EQ.'OUT-STEP-LIM') THEN
*       Set maximum step for line optimization.
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.EQ.2) THEN
          SR=FLOTT
        ELSE IF(INDIC.EQ.4) THEN
          SR=DFLOTT
        ELSE
          CALL XABORT('LNSR: REAL OR DOUBLE PRECISION VALUE EXPECTED.')
        ENDIF
      ELSE IF(TEXT12.EQ.'INN-STEP-EPS') THEN
*       Set the tolerence used for line optimization.
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.EQ.2) THEN
          EPS3=FLOTT
        ELSE IF(INDIC.EQ.4) THEN
          EPS3=DFLOTT
        ELSE
          CALL XABORT('LNSR: REAL OR DOUBLE PRECISION VALUE EXPECTED.')
        ENDIF
      ELSE IF(TEXT12.EQ.'OUT-STEP-EPS') THEN
*       Set the tolerence used for external iterations.
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.EQ.2) THEN
          EPS2=FLOTT
        ELSE IF(INDIC.EQ.4) THEN
          EPS2=DFLOTT
        ELSE
          CALL XABORT('LNSR: REAL OR DOUBLE PRECISION VALUE EXPECTED.')
        ENDIF
      ELSE IF(TEXT12.EQ.'OUT-ITER-MAX') THEN
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('LNSR: INTEGER DATA EXPECTED(2).')
        IF(MAXEXT.EQ.0) MAXEXT=NITMA
      ELSE IF(TEXT12.EQ.'OUT-RESTART') THEN
        CALL REDGET(INDIC,NSTART,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('LNSR: INTEGER DATA EXPECTED(3).')
      ELSE IF(TEXT12.EQ.'SD') THEN
        IHESS=0
      ELSE IF(TEXT12.EQ.'CG') THEN
        IHESS=1
      ELSE IF(TEXT12.EQ.'BFGS') THEN
        IHESS=2
      ELSE IF(TEXT12.EQ.'LBFGS') THEN
        IHESS=3
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.EQ.1) THEN
*         hist_nr: number of corrections stored in LBFGS method
          hist_nr=NITMA
        ELSE IF(INDIC.EQ.3) THEN
          GO TO 15
        ELSE
          CALL XABORT('LNSR: INTEGER OR CHARACTER VALUE EXPECTED.')
        ENDIF
      ELSE IF(TEXT12.EQ.'NEWT') THEN
        IHESS=4
      ELSE IF(TEXT12.EQ.'INN-CONV-TST') THEN
*       Internal convergence test
        IPICK=1
        GO TO 20
      ELSE IF(TEXT12(:1).EQ.';') THEN
        GO TO 20
      ELSE 
        CALL XABORT('LNSR: '//TEXT12//' IS AN INVALID KEYWORD.')
      ENDIF
      GO TO 10
*----
*  RECOVER INFORMATION FROM OPTIM OBJECT
*----
   20 ISTEP=ISTEP+1
      ALLOCATE(X(NVAR),P(NVAR),G(NVAR),XMIN(NVAR),XMAX(NVAR))
      IF(IMETH.EQ.4) THEN
        ALLOCATE(FF(NFUNC))
        CALL LCMGET(IPGRAD,'FOBJ-CST-VAL',FF)
        F=DOT_PRODUCT(FF(:NFUNC),FF(:NFUNC))
        DEALLOCATE(FF)
      ELSE
        CALL LCMGET(IPGRAD,'FOBJ-CST-VAL',F)
        CALL LCMGET(IPGRAD,'GRADIENT',G)
      ENDIF
      CALL LCMGET(IPGRAD,'VAR-VALUE',X)
      CALL LCMGET(IPGRAD,'VAR-VAL-MIN',XMIN)
      CALL LCMGET(IPGRAD,'VAR-VAL-MAX',XMAX)
      CALL LCMLEN(IPGRAD,'LNSR-INFO',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
        CALL LCMGET(IPGRAD,'LNSR-INFO',DSAVE)
        SLOPE=DSAVE(1)
        ALAM=DSAVE(2)
        GNORM=DSAVE(3)
      ELSE
        SLOPE=0.0D0
        ALAM=0.0D0
        GNORM=0.0D0
      ENDIF
*----
*  SET THE DIRECTION AND INITIALIZATION OF THE LINE SEARCH
*----
      IF(ISTEP.EQ.1) THEN
        IF(IPRINT.GT.0) WRITE(6,100) IEXT,F
        IF(IHESS.EQ.0) THEN
*         Steepest descent
          P(:NVAR)=-G(:NVAR)
        ELSE IF(IHESS.EQ.1) THEN
          IF(IEXT.EQ.1) THEN
*           Steepest descent
            P(:NVAR)=-G(:NVAR)
            GNORM=DOT_PRODUCT(G(:NVAR),G(:NVAR))/REAL(NVAR)
          ELSE
*           Conjugate gradient
            GNORMP=GNORM
            CALL LCMGET(IPGRAD,'DIRECTION',P)
            GNORM=DOT_PRODUCT(G(:NVAR),G(:NVAR))/REAL(NVAR)
            P(:NVAR)=-G(:NVAR)+(GNORM/GNORMP)*P(:NVAR)
          ENDIF
        ELSE IF(IHESS.EQ.2) THEN
*         BFGS
          IF(IEXT.EQ.1) THEN
            ALLOCATE(GG(NVAR,NVAR))
            GG(:NVAR,:NVAR)=0.0D0
            DO I=1,NVAR
              GG(I,I)=1.0D0
            ENDDO
*           Steepest descent
            P(:NVAR)=-G(:NVAR)
          ELSE
            ALLOCATE(V(NVAR),Y(NVAR),XP(NVAR),GP(NVAR))
            CALL LCMSIX(IPGRAD,'OLD-VALUE',1)
              CALL LCMGET(IPGRAD,'VAR-VALUE',XP)
              CALL LCMGET(IPGRAD,'GRADIENT',GP)
            CALL LCMSIX(IPGRAD,' ',2)
            V(:NVAR)=X(:NVAR)-XP(:NVAR)
            Y(:NVAR)=G(:NVAR)-GP(:NVAR)
            SVY=DOT_PRODUCT(V(:NVAR),Y(:NVAR))
            IF(SVY.EQ.0.0D0) CALL XABORT('LNSR: DIVIDE CHECK IN BFGS.')
            DEALLOCATE(GP,XP)
            ALLOCATE(GG(NVAR,NVAR),GGY(NVAR),YGG(NVAR),AA(NVAR,NVAR))
            CALL LCMGET(IPGRAD,'HESSIAN',GG)
            SVYI=1.0D0/SVY
            DO I=1,NVAR
              TMP1=0.0D0
              TMP2=0.0D0
              DO J=1,NVAR
                AA(J,I)=V(J)*V(I)*SVYI
                TMP1=TMP1+GG(I,J)*Y(J)
                TMP2=TMP2+Y(J)*GG(J,I)
              ENDDO
              GGY(I)=TMP1
              YGG(I)=TMP2
            ENDDO
            B=1.0D0
            DO I=1,NVAR
              B=B+Y(I)*GGY(I)*SVYI
            ENDDO
            AA(:NVAR,:NVAR)=AA(:NVAR,:NVAR)*B
            DO J=1,NVAR
              DO I=1,NVAR
                AA(I,J)=AA(I,J)-(V(I)*YGG(J)+GGY(I)*V(J))*SVYI
              ENDDO
            ENDDO
            GG(:NVAR,:NVAR)=GG(:NVAR,:NVAR)+AA(:NVAR,:NVAR)
            P(:NVAR)= 0.0D0
            DO I=1,NVAR
              P(:NVAR)=P(:NVAR)-GG(:NVAR,I)*G(I)
            ENDDO
            DEALLOCATE(AA,YGG,GGY,Y,V)
          ENDIF
          CALL LCMPUT(IPGRAD,'HESSIAN',NVAR*NVAR,4,GG)
          DEALLOCATE(GG)
        ELSE IF(IHESS.EQ.3) THEN
*         Limited memory BFGS
          ALLOCATE(SS(NVAR,hist_nr),YY(NVAR,hist_nr))
          P(:NVAR)=G(:NVAR)
          IF(IEXT.EQ.1) THEN
            CALL XDDSET(SS,NVAR*hist_nr,0.0D0)
            CALL XDDSET(YY,NVAR*hist_nr,0.0D0)
          ELSE
*           quasi-Newton search
            ALLOCATE(GAMMA(hist_nr),XP(NVAR),GP(NVAR))
            CALL LCMGET(IPGRAD,'LBFGS-S',SS)
            CALL LCMGET(IPGRAD,'LBFGS-Y',YY)
            CALL LCMSIX(IPGRAD,'OLD-VALUE',1)
              CALL LCMGET(IPGRAD,'VAR-VALUE',XP)
              CALL LCMGET(IPGRAD,'GRADIENT',GP)
            CALL LCMSIX(IPGRAD,' ',2)
            J=MOD(IEXT-1,hist_nr)+1
            SS(:NVAR,J)=X(:NVAR)-XP(:NVAR)
            YY(:NVAR,J)=G(:NVAR)-GP(:NVAR)
            SVY=DOT_PRODUCT(SS(:NVAR,J),YY(:NVAR,J))
            IF(SVY.EQ.0.0D0) CALL XABORT('LNSR: DIVIDE CHECK IN LBFGS.')
            DEALLOCATE(GP,XP)
            IBOUND=MIN(IEXT-1,hist_nr)
            DO IB=IBOUND,1,-1
              J=MOD(IEXT+IB-IBOUND-1,hist_nr)+1
              TAU=DOT_PRODUCT(SS(:NVAR,J),YY(:NVAR,J))
              GAMMA(IB)=DOT_PRODUCT(SS(:NVAR,J),P(:NVAR))/TAU
              P(:NVAR)=P(:NVAR)-GAMMA(IB)*YY(:NVAR,J)
            ENDDO
            DO IB=1,IBOUND
              J=MOD(IEXT+IB-IBOUND-1,hist_nr)+1
              TAU=DOT_PRODUCT(SS(:NVAR,J),YY(:NVAR,J))
              BETA=DOT_PRODUCT(YY(:NVAR,J),P(:NVAR))/TAU
              P(:NVAR)=P(:NVAR)+(GAMMA(IB)-BETA)*SS(:NVAR,J)
            ENDDO
            DEALLOCATE(GAMMA)
          ENDIF
          CALL LCMPUT(IPGRAD,'LBFGS-S',NVAR*hist_nr,4,SS)
          CALL LCMPUT(IPGRAD,'LBFGS-Y',NVAR*hist_nr,4,YY)
          DEALLOCATE(YY,SS)
          P(:NVAR)=-P(:NVAR)
        ELSE IF(IHESS.EQ.4) THEN
*         Newton method for unconstrained optimization
          ALLOCATE(FF(NFUNC),DFF(NVAR,NFUNC),TDFF(NFUNC,NVAR),
     1    UD(NVAR))
          CALL LCMGET(IPGRAD,'FOBJ-CST-VAL',FF)
          CALL LCMGET(IPGRAD,'GRADIENT',DFF)
          G(:NVAR)=2.0D0*MATMUL(DFF,FF)
          TDFF=TRANSPOSE(DFF)
          CALL ALST2F(NFUNC,NFUNC,NVAR,TDFF,UD)
          CALL ALST2S(NFUNC,NFUNC,NVAR,TDFF,UD,FF,P)
          P(:NVAR)=-P(:NVAR)
          DEALLOCATE(UD,TDFF,DFF,FF)
        ENDIF
        GNORM=DOT_PRODUCT(G(:NVAR),G(:NVAR))/REAL(NVAR)
        PABS=SQRT(DOT_PRODUCT(P(:NVAR),P(:NVAR)))
        P(:NVAR)=P(:NVAR)*SR/PABS ! stepsize normalization
        SLOPE=DOT_PRODUCT(G(:NVAR),P(:NVAR))
        ALAM=1.0D0
        IF(IOPT.EQ.-1) F=-F
        FOLD=F
        CALL LCMPUT(IPGRAD,'DIRECTION',NVAR,4,P)
        CALL LCMSIX(IPGRAD,'OLD-VALUE',1)
          CALL LCMPUT(IPGRAD,'VAR-VALUE',NVAR,4,X)
          CALL LCMPUT(IPGRAD,'FOBJ-CST-VAL',1,4,FOLD)
        CALL LCMSIX(IPGRAD,' ',2)
        GO TO 30
      ELSE
*       recover values at beginning of line search
        CALL LCMGET(IPGRAD,'DIRECTION',P)
        CALL LCMSIX(IPGRAD,'OLD-VALUE',1)
          CALL LCMGET(IPGRAD,'VAR-VALUE',X)
          CALL LCMGET(IPGRAD,'FOBJ-CST-VAL',FOLD)
          IF(IOPT.EQ.-1) FOLD=-FOLD
        CALL LCMSIX(IPGRAD,' ',2)
      ENDIF
*----
*  SINGLE INNER ITERATION OF THE LINE OPTIMIZATION
*----
      IF(IOPT.EQ.-1) F=-F
      IF(F.LE.FOLD+XI*ALAM*SLOPE) THEN
*       Armijo condition
        JCONV =1
        GO TO 40
      ELSE IF(ISTEP.GT.MAXINT) THEN
        JCONV =2
        GO TO 40
      ENDIF
      ALAM=ALAM*WIDTH
   30 X(:NVAR)=X(:NVAR)+ALAM*P(:NVAR)
      DO I=1,NVAR
        X(I)=MAX(XMIN(I),MIN(XMAX(I),X(I)))
      ENDDO
      CALL LCMPUT(IPGRAD,'VAR-VALUE',NVAR,4,X)
   40 DEALLOCATE(XMAX,XMIN,G,P,X)
      IF(IPRINT.GT.0) WRITE(6,110) IEXT,ISTEP,ALAM,F,JCONV
      IF(IPRINT.GT.2) THEN
        ALLOCATE(X(NVAR),P(NVAR))
        CALL LCMGET(IPGRAD,'DIRECTION',P)
        CALL LCMGET(IPGRAD,'VAR-VALUE',X)
        WRITE(6,120) '    LINE SEARCH DIRECTION',(P(I),I=1,NVAR)
        WRITE(6,120) 'OUTPUT DECISION VARIABLES',(X(I),I=1,NVAR)
        DEALLOCATE(P,X)
      ENDIF
*----
*  TEST FOR EXTERNAL ITERATION CONVERGENCE
*----
      IF(JCONV.GE.1) THEN
        DNVTST=1
        TEST2=ABS(F-FOLD)
        IF(GNORM.LT.0.01*EPS2) THEN
          IF(IPRINT.GT.0) PRINT *,'>>> OUTER CONVERGED WRT GNORM'
          CNVTST=1
          ICONV =1
        ELSE IF((TEST2.LT.EPS2).AND.(ISTEP.GT.1)) THEN
          IF(IPRINT.GT.0) PRINT *,'>>> OUTER CONVERGED WRT F-FOLD'
          CNVTST=1
          ICONV =1
        ELSE IF(IEXT.GE.MAXEXT) THEN
          IF(IPRINT.GT.0) PRINT *,'>>> OUTER REACHES MAXIMUM ITERATION'
          CNVTST=1
          ICONV =1
        ENDIF
        IF(IPRINT.GT.0) WRITE(6,130) IEXT,ABS(ALAM),GNORM,TEST2,EPS2
*----
*  RESTART CG OR BFGS HESSIAN MATRIX CALCULATION
*----
        IF((NSTART.NE.0).AND.(IEXT.GE.NSTART)) THEN
          IEXT=0
          MAXEXT=MAXEXT-NSTART
        ENDIF
*----
*  SAVE OLD GRADIENT
*----
        ALLOCATE(G(NVAR),P(NVAR))
        CALL LCMGET(IPGRAD,'GRADIENT',G)
        CALL LCMGET(IPGRAD,'DIRECTION',P)
        CALL LCMSIX(IPGRAD,'OLD-VALUE',1)
          CALL LCMPUT(IPGRAD,'GRADIENT',NVAR,4,G)
        CALL LCMSIX(IPGRAD,' ',2)
        DEALLOCATE(P,G)
        IEXT=IEXT+1
      ENDIF
*----
*  SAVE THE STATE VECTORS
*----
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NVAR
      ISTATE(3)=IOPT
      ISTATE(4)=ICONV
      ISTATE(5)=IEXT
      ISTATE(6)=IEDSTP
      ISTATE(7)=IHESS
      ISTATE(8)=IMETH
      ISTATE(10)=ISTEP
      ISTATE(11)=JCONV
      ISTATE(12)=MAXEXT
      ISTATE(13)=NSTART
      IF(IPRINT.GT.0) WRITE(6,140) (ISTATE(I),I=1,13)
      CALL LCMPUT(IPGRAD,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL XDDSET(OPTPRR,NSTATE,0.0D0)
      OPTPRR(1)=SR
      OPTPRR(2)=EPS1
      OPTPRR(3)=EPS2
      OPTPRR(4)=EPS3
      IF(IPRINT.GT.0) WRITE(6,150) (OPTPRR(I),I=1,4)
      CALL LCMPUT(IPGRAD,'OPT-PARAM-R',NSTATE,4,OPTPRR)
      DSAVE(1)=SLOPE
      DSAVE(2)=ALAM
      DSAVE(3)=GNORM
      CALL LCMPUT(IPGRAD,'LNSR-INFO',3,4,DSAVE)
      IF(IPRINT.GT.2) CALL LCMLIB(IPGRAD)
*----
*  RECOVER THE CONVERGENCE FLAGS AND SAVE IT IN A CLE-2000 VARIABLE
*----
      IF(IPICK.EQ.1) THEN
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.-5) CALL XABORT('LNSR: OUTPUT LOGICAL EXPECTED.')
        INDIC=5
        CALL REDPUT(INDIC,DNVTST,FLOTT,TEXT12,DFLOTT)
   50   CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.3) CALL XABORT('LNSR: CHARACTER DATA EXPECTED(2).')
        IF(TEXT12.EQ.'OUT-CONV-TST') THEN
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
          IF(INDIC.NE.-5) CALL XABORT('LNSR: OUTPUT LOGICAL EXPECTED.')
          INDIC=5
          CALL REDPUT(INDIC,CNVTST,FLOTT,TEXT12,DFLOTT)
          GO TO 50
        ELSE IF (TEXT12.EQ.';') THEN
          RETURN
        ELSE
          CALL XABORT('LNSR: ; CHARACTER EXPECTED.')
        ENDIF      
      ENDIF      
      RETURN
*
  100 FORMAT(/14H LNSR: ##ITER=,I8,20H OBJECTIVE FUNCTION=,1P,E14.6)
  110 FORMAT(/21H LNSR: EXTERNAL ITER=,I5,18H LINE SEARCH ITER=,I4,
     1 7H ALPHA=,1P,E17.10,20H OBJECTIVE FUNCTION=,E17.10,6H CONV=,I2)
  120 FORMAT(/7H LNSR: ,A25,1H=,1P,8E12.4/(33X,8E12.4))
  130 FORMAT(/26H LNSR: EXTERNAL ITERATION=,I4,12H ACCURACIES=,1P,
     1 3E12.4,6H EPS2=,E12.4)
  140 FORMAT(/8H OPTIONS/8H -------/
     1 7H NVAR  ,I8,32H   (NUMBER OF CONTROL VARIABLES)/
     2 7H NCST  ,I8,26H   (NUMBER OF CONSTRAINTS)/
     3 7H IOPT  ,I8,37H   (=1/-1: MINIMIZATION/MAXIMIZATION)/
     4 7H ICONV ,I8,43H   (=0/1: EXTERNAL NOT CONVERGED/CONVERGED)/
     5 7H IEXT  ,I8,32H   (INDEX OF EXTERNAL ITERATION)/
     6 7H IEDSTP,I8,13H   (NOT USED)/
     7 7H IHESS ,I8,46H   (=0/1/2/3/4: STEEPEST/CG/BFGS/LBFGS/NEWTON)/
     8 7H ISEARC,I8,35H   (=0/1/2: NO SEARCH/OPTEX/NEWTON)/
     9 7H IMETH ,I8,13H   (NOT USED)/
     1 7H ISTEP ,I8,35H   (INDEX OF LINE SEARCH ITERATION)/
     2 7H JCONV ,I8,48H   (=0/1/2: LINE SEARCH NOT CONVERGED/CONVERGED)/
     3 7H MAXEXT,I8,42H   (MAXIMUM NUMBER OF EXTERNAL ITERATIONS)/
     4 7H NSTART,I8,37H   (EXTERNAL ITERATION RESTART CYCLE))
  150 FORMAT(/
     1 12H REAL PARAM:,1P/12H -----------/
     2 7H SR    ,D12.4,33H   (MAXIMUM LINE SEARCH STEPSIZE)/
     3 7H EPS1  ,D12.4,13H   (NOT USED)/
     4 7H EPS2  ,D12.4,31H   (EXTERNAL CONVERGENCE LIMIT)/
     5 7H EPS3  ,D12.4,31H   (INTERNAL CONVERGENCE LIMIT))
      END
