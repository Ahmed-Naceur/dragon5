*DECK PLSPLX
      SUBROUTINE PLSPLX(N,M,MAXM,MINMAX,COUT,APLUS,B,INEGAL,BINF,BSUP,
     > XOBJ,F,EPS,IMTHD,IMPR,IERR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solves a linear problem using the revisited simplex method. 
* PLSLPX = Linear Programmation SimPLeX
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert and R. Chambon
*
*Parameters: input/ouput
* N       number of control variables.
* M       number of constraints.
* MAXM    first dimension of matrix APLUS.
* MINMAX  type of optimization (=-1: minimize; =1: maximize).
* COUT    costs of control variables.
* APLUS   coefficient matrix for the linear constraints.
* B       right hand sides corresponding to the coefficient matrix.
* INEGAL  constraint relations (=-1 for .GE.; =0 for .EQ.; =1 for .LE.).
* BINF    lower bounds of control variables.
* BSUP    upper bounds of control variables.
* EPS     tolerence used for SIMPLEX calculation.
* IMTHD   type of solution (=1: SIMPLEX/LEMKE; =3: MAP).
* IMPR    print flag.
*
*Parameters: ouput
* XOBJ    control variables.
* F       objective function.
* IERR    return code (=0: normal completion).
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N,M,MAXM,MINMAX,INEGAL(MAXM),IMTHD,IMPR,IERR
      DOUBLE PRECISION B(M+1),BINF(N),BSUP(N),XOBJ(N),EPS,COUT(N),
     > APLUS(MAXM,N+M),F
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION   DUMY,DELTA,DELTAM,CSTE,XIJ,XIKI,BBAR,PM,QM
      INTEGER  I,J,MP1,MP2,IDIMC,JJ,IJK,J1,IPHASE,IND,K,IP,IQ,L,LL,
     >         LARTF
      LOGICAL  LARTIF
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IVARS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: CRE,BORNE,XIK,BGAR,
     > GSUP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: BASE0,AGAR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IVARS(M+2))
      ALLOCATE(CRE(N+M),BORNE(N+M),XIK(M+2),BGAR(M+2),GSUP(N+M))
      ALLOCATE(BASE0(M+2,M+2),AGAR(2,N+M))
*
      CALL XDDSET(CRE,N+M,0.0D0)
      LARTIF = .FALSE.
      IF(M.GT.MAXM) THEN
         IERR = 9
         GO TO 500
      ENDIF
      MP1 = M+1
      MP2 = M+2
*----
*  ORGANIZATION OF THE SIMPLEX TABLES
*----
      DO 3 I=1,M
         DUMY = 0.0D0
         DO 1 J=1,N
            DUMY = DUMY + APLUS(I,J)*BINF(J)
    1    CONTINUE
*
         BGAR(I) = B(I)-DUMY
         IF(BGAR(I).GE.0.0) GOTO 3
         DO 2 J=1,N
            APLUS(I,J) = -APLUS(I,J)
    2    CONTINUE
*
         BGAR(I)   = -BGAR(I)
         B(I)      = -B(I)
         INEGAL(I) = -INEGAL(I)
    3 CONTINUE
*
      CSTE = 0.0D0
      DO 4 J=1,N
         CSTE      = CSTE + COUT(J)*BINF(J)
         GSUP(J)   = BSUP(J) - BINF(J)
         AGAR(1,J) = REAL(MINMAX)*COUT(J)
    4 CONTINUE
*
      DO 8 J=1,N
         AGAR(2,J) = 0.0D0
         DO 7 I=1,M
            AGAR(2,J) = AGAR(2,J) - APLUS(I,J)
    7    CONTINUE
    8 CONTINUE
*
      BGAR(MP1) = 0.0
      BGAR(MP2) = 0.0
*
      DO 9 I=1,M
         BGAR(MP2) = BGAR(MP2) - BGAR(I)
    9 CONTINUE
*----
*  SLACK VARIABLES
*----
      IDIMC = N
      DO 13 I=1,M
         IF(INEGAL(I).EQ.0) GOTO 13
         IDIMC = IDIMC + 1
         GSUP(IDIMC) = 1.0E+25
         DO 12 J=1,M
            APLUS(J,IDIMC) = 0.0D0
   12    CONTINUE
*
         APLUS(I,IDIMC) = REAL(INEGAL(I))
         AGAR(1,IDIMC)  = 0.0D0
         AGAR(2,IDIMC)  = -REAL(INEGAL(I))
   13 CONTINUE
*
      DO 15 I=1,MP2
         IVARS(I) = IDIMC + I
         DO 14 J=1,MP2
            BASE0(I,J) = 0.0D0
   14    CONTINUE
         BASE0(I,I) = 1.0D0
   15 CONTINUE
*
      IF(IMPR.GE.6) THEN
         WRITE (6,800)
         WRITE (6,840) N,M,IDIMC
         DO 510 J=1,N+M,12
            JJ = MIN0(J+11,N+M)
            WRITE (6,841) (J1,J1=J,JJ)
            DO 505 I=1,M
               WRITE (6,805) I,(APLUS(I,J1),J1=J,JJ)
  505       CONTINUE
            IF(JJ.NE.N+M) WRITE (6,803)
  510    CONTINUE
*
         WRITE (6,810) (I,COUT(I),I=1,N)
         WRITE (6,820) (I,BGAR(I),I=1,MP2)
         WRITE (6,830) (I,INEGAL(I),I=1,M)
         WRITE (6,835) (I,GSUP(I),I=1,IDIMC)
      ENDIF
*----
*  END OF PHASE 1 IF BGAR(M+2)=0 +- EPS
*----
      IERR = 0
      IND=0
      LARTF=0
      DO 11 I=1,IDIMC
         IF(GSUP(I).EQ.0.0) GSUP(I)=1.0E-6
         BORNE(I)=0.0
   11 CONTINUE
*
      IPHASE = 1
*
   10 IF((IPHASE.EQ.2) .AND. (.NOT.LARTIF)) GOTO 21
      IF((ABS(BGAR(MP2)).LT.EPS) .OR. (LARTIF)) THEN
        IPHASE=2
        LARTIF=.FALSE.
        IJK=1
   99   IF(.NOT.((IJK.GT.IDIMC) .OR. (LARTIF))) THEN
          IF(IVARS(IJK).GT.IDIMC) THEN
            LARTIF=.TRUE.
            LARTF=IJK
          ENDIF
          IJK=IJK+1
          GOTO 99
        ENDIF
      ENDIF
*----
*  RETURN IF NO BASE SWITCH HAS BEEN DONE
*----
   20 IND=MP2
   21 IF(IPHASE.EQ.2) IND=MP2-1
      IF(IMPR.GE.6) THEN
         WRITE (6,845) IPHASE
         WRITE (6,843)
         WRITE (6,850)
         DO 520 J=1,MP2,12
         JJ=MIN0(J+11,MP2)
         WRITE (6,841) (J1,J1=J,JJ)
         DO 515 I=1,MP2
         WRITE (6,805) I,(BASE0(I,J1),J1=J,JJ)
  515    CONTINUE
         IF(JJ.NE.MP2) WRITE (6,803)
  520    CONTINUE
         WRITE (6,820) (I,BGAR(I),I=1,MP2)
         WRITE (6,855) (I,IVARS(I),I=1,MP2)
         WRITE (6,860) (I,BORNE(I),I=1,IDIMC)
      ENDIF
*----
*  SELECTION OF THE BASE K VECTOR TO PICK IN
*----
      K=0
      DELTAM=0.0D0
      DO 100 J=1,IDIMC
      DELTA=0.0D0
      DO 30 JJ=1,MP2
      IF(IVARS(JJ).EQ.J) GOTO 90
   30 CONTINUE
      DELTA=BASE0(IND,MP1)*AGAR(1,J)+BASE0(IND,MP2)*AGAR(2,J)
      DO 35 I=1,M
        DELTA=DELTA+BASE0(IND,I)*APLUS(I,J)
   35 CONTINUE
      IF(IMPR.GE.7) WRITE (6,865) J,DELTA
*
      IF(BORNE(J).EQ.GSUP(J)) GOTO 45
      IF(DELTA.GE.0.0) GOTO 90
      IF(K.EQ.0) GOTO 40
   39 IF(DELTA.GE.DELTAM) GOTO 90
   40 K=J
      DELTAM=DELTA
      GOTO 90
   45 IF(DELTA.LE.0.0) GOTO 90
      DELTA=-DELTA
      IF(K.EQ.0) GOTO 40
      IF(K.GT.0) GOTO 39
*
   90 CRE(J)=DELTA
  100 CONTINUE
*
      IF(IMPR.GE.7) THEN
         WRITE (6,870) (I,CRE(I),I=1,N+M)
         WRITE (6,875) K
      ENDIF
*
      IF(K.NE.0) GOTO 125
      DO 105 JJ=1,M
      IF(IVARS(JJ).GT.IDIMC) GOTO 110
  105 CONTINUE
      IF(IPHASE.EQ.2) GOTO 115
      IERR=1
      GO TO 500
  110 IF(IPHASE.EQ.2) GOTO 115
      IF(IMTHD.EQ.1) THEN
*        no solution.
         IERR=2
         GO TO 500
      ELSEIF(IMTHD.EQ.3) THEN
*        no feasible solution in the reduced domain. Go to PHASE 2 with
*        modified constraints.      
         IERR=2
         GO TO 115
      ENDIF
*
  115 DO 999 I=1,MP2
      DO 950 J=1,IDIMC
      XIJ=0.0D0
      DO 920 JJ=1,M
      IF(J.EQ.IVARS(JJ)) GOTO 950
  920 CONTINUE
      IF(BORNE(J).NE.GSUP(J)) GOTO 950
      XIJ=BASE0(I,MP1)*AGAR(1,J)+BASE0(I,MP2)*AGAR(2,J)
      DO 925 JJ=1,M
        XIJ=XIJ+BASE0(I,JJ)*APLUS(JJ,J)
  925 CONTINUE
      BGAR(I)=BGAR(I)-XIJ*GSUP(J)
  950 CONTINUE
  999 CONTINUE
      DO 400 I=1,N
      DO 390 J=1,M
      IF(IVARS(J).NE.I) GO TO 390
      XOBJ(I)=BINF(I)+BGAR(J)
      GO TO 400
  390 CONTINUE
      XOBJ(I)=BINF(I)+BORNE(I)
  400 CONTINUE
      DO 410 I=1,IDIMC
        CRE(I)=CRE(I)*REAL(MINMAX)
  410 CONTINUE
      F=REAL(-MINMAX)*BGAR(MP1)+CSTE
      IF(IMPR.GE.5) THEN
         WRITE (6,890) F,(I,XOBJ(I),I=1,N)
         WRITE (6,820) (I,BGAR(I),I=1,MP2)
         WRITE (6,855) (I,IVARS(I),I=1,MP2)
         WRITE (6,860) (I,BORNE(I),I=1,IDIMC)
      ENDIF
      GO TO 500
*----
*  SELECTION OF THE VECTOR TO PICK OUT
*----
  125 IF(IMPR.GE.7) WRITE (6,880)
      IP=0
      IQ=0
*
      IF(K.GT.N+M) THEN
         WRITE(6,'(1X,A)') 'PLSPLX: AGAR OVERFLOW.'
         IERR=10
         GO TO 500
      ENDIF
      DO 200 I=1,M
*
      BBAR=0.0D0
      DO 150 J=1,IDIMC
      IF(BORNE(K).EQ.GSUP(K) .AND. J.EQ.K) GOTO 150
      DO 130 JJ=1,M
      IF(J.EQ.IVARS(JJ)) GOTO 150
  130 CONTINUE
      IF(BORNE(J).NE.GSUP(J)) GOTO 150
      XIJ=BASE0(I,MP1)*AGAR(1,J)+BASE0(I,MP2)*AGAR(2,J)
      DO 135 JJ=1,M
        XIJ=XIJ+BASE0(I,JJ)*APLUS(JJ,J)
  135 CONTINUE
      BBAR=BBAR-GSUP(J)*XIJ
      IF(IMPR.GE.7) THEN
         WRITE (6,885) I,J,XIJ,BBAR
      ENDIF
  150 CONTINUE
      BBAR=BBAR+BGAR(I)
      XIKI=BASE0(I,MP1)*AGAR(1,K)+BASE0(I,MP2)*AGAR(2,K)
      DO 155 JJ=1,M
        XIKI=XIKI+BASE0(I,JJ)*APLUS(JJ,K)
  155 CONTINUE
      XIK(I)=XIKI
*----
*  CASE I: XOBJ(K)=0.0
*----
      IF(BORNE(K).EQ.GSUP(K)) GOTO 175
      IF(ABS(XIKI).LT.EPS) GOTO 200
      IF(XIKI.GT.0.0) GOTO 165
      IF(LARTIF.AND.I.EQ.LARTF) GOTO 165
*----
*  TEST FOR ARTIFICIAL VARIABLES
*----
      IF(IVARS(I).GT.IDIMC) GOTO 200
      IF(IQ.EQ.0) GOTO 160
      IF(QM.LE.((GSUP(IVARS(I))-BBAR)/(-XIKI))) GOTO 200
  160 CONTINUE
      IQ=I
      QM=(GSUP(IVARS(I))-BBAR)/(-XIKI)
      GOTO 200
  165 CONTINUE
      IF(IP.EQ.0) GOTO 170
      IF(PM.LE.(BBAR/XIKI)) GOTO 200
  170 IP=I
      PM=BBAR/XIKI
      GOTO 200
*----
*  CASE II: XOBJ(K) = UPPER BOUND
*----
  175 IF(ABS(XIKI).LT.EPS) GOTO 200
      IF(XIKI.GT.0.0) GOTO 185
      IF(IQ.EQ.0) GOTO 180
      IF(QM.GE.(BBAR/XIKI)) GOTO 200
  180 IQ=I
      QM=BBAR/XIKI
      GOTO 200
*----
*  TEST FOR ARTIFICIAL VARIABLES
*----
  185 IF(IVARS(I).GT.IDIMC) GOTO 200
      IF(IP.EQ.0) GOTO 190
      IF(PM.GE.((BBAR-GSUP(IVARS(I)))/XIKI)) GOTO 200
  190 IP=I
      PM=(BBAR-GSUP(IVARS(I)))/XIKI
  200 CONTINUE
*
      IF(IMPR.GE.7) WRITE (6,894) IQ,QM,IP,PM
*
      XIK(MP1)=BASE0(MP1,MP1)*AGAR(1,K)+BASE0(MP1,MP2)*AGAR(2,K)
      XIK(MP2)=BASE0(MP2,MP1)*AGAR(1,K)+BASE0(MP2,MP2)*AGAR(2,K)
      DO 204 JJ=1,M
      XIK(MP1)=XIK(MP1)+BASE0(MP1,JJ)*APLUS(JJ,K)
        XIK(MP2)=XIK(MP2)+BASE0(MP2,JJ)*APLUS(JJ,K)
  204 CONTINUE
*
      IF(BORNE(K).EQ.GSUP(K)) GOTO 250
      IF(IP.EQ.0) GOTO 205
      IF(IQ.EQ.0) GOTO 220
      IF(QM.LE.PM) GOTO 210
      GOTO 220
  205 IF(IQ.EQ.0) GOTO 211
  210 IF(QM.LE.GSUP(K)) GOTO 215
  211 BORNE(K)=GSUP(K)
      GOTO 20
  215 L=IQ
      LL=IVARS(L)
      IF(LL.LE.IDIMC) BORNE(LL)=GSUP(LL)
      IVARS(L)=K
      IF(IMPR.GE.7) WRITE (6,895) L,LL,IVARS(L)
      GOTO 300
  220 IF(PM.GT.GSUP(K)) GOTO 211
      L=IP
      LL=IVARS(L)
      IF(LL.LE.IDIMC) BORNE(LL)=0.0
      IVARS(L)=K
      IF(IMPR.GE.7) WRITE (6,895) L,LL,IVARS(L)
      GOTO 300
*
  250 IF(IP.EQ.0) GOTO 255
      IF(IQ.EQ.0) GOTO 265
      IF(QM.GE.PM) GOTO 256
      GOTO 265
  255 IF(IQ.EQ.0) GOTO 257
  256 IF(QM.GE.0.0) GOTO 260
  257 BORNE(K)=0.0
      GOTO 20
  260 L=IQ
      LL=IVARS(L)
      IF(LL.LE.IDIMC) BORNE(LL)=0.0
      IVARS(L)=K
      IF(IMPR.GE.7) WRITE (6,895) L,LL,IVARS(L)
      GOTO 300
  265 IF(PM.LE.0.0) GOTO 257
      L=IP
      LL=IVARS(L)
      IF(LL.LE.IDIMC) BORNE(LL)=GSUP(LL)
      IVARS(L)=K
      IF(IMPR.GE.7) WRITE (6,895) L,LL,IVARS(L)
*
*     PIVOTAGE SUR XIK
  300 IF(IMPR.GE.7) WRITE (6,905) L,XIK(L)
      DO 325 I=1,MP2
      IF(I.EQ.L) GOTO 325
      BGAR(I)=BGAR(I)-BGAR(L)/XIK(L)*XIK(I)
  325 CONTINUE
      BGAR(L)=BGAR(L)/XIK(L)
*
      DO 375 I=1,MP2
      IF(I.EQ.L) GOTO 375
      DO 335 J=1,MP2
      BASE0(I,J)=BASE0(I,J)-(BASE0(L,J)/XIK(L))*XIK(I)
  335 CONTINUE
  375 CONTINUE
      DO 380 J=1,MP2
        BASE0(L,J)=BASE0(L,J)/XIK(L)
  380 CONTINUE
*
      GOTO 10
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  500 DEALLOCATE(AGAR,BASE0)
      DEALLOCATE(GSUP,BGAR,XIK,BORNE,CRE)
      DEALLOCATE(IVARS)
      RETURN
*
*     .................................................................
*     ...
*     ...     F O R M A T S
*     ...
*     .................................................................
*
  800 FORMAT (///25X,35H- - - -   SIMPLEX ALGORITHM   - - -)
  803 FORMAT (//10X,8HCONTN.../)
  805 FORMAT (1X,I4,1P,12E10.2)
  810 FORMAT (//10X,28H+ + + + +    TABLE COUT(N) :/
     1 3(/6(1X,5HCOUT(,I4,4H) = ,1P,E10.3,:)))
  820 FORMAT (//10X,30H+ + + + +    TABLE BGAR(M+2) :/
     1 3(/6(1X,5HBGAR(,I4,4H) = ,1P,E10.3,:)))
  830 FORMAT (//10X,30H+ + + + +    TABLE INEGAL(M) :/
     1 3(/6(1X,5HINEG(,I4,4H) = ,I4,:)))
  835 FORMAT (//10X,32H+ + + + +    TABLE GSUP(IDIMC) :/
     1 5(/6(1X,5HGSUP(,I4,4H) = ,1P,E10.3,:)))
  840 FORMAT (//34H-- NUMBER OF CONTROL VARIABLES (N):,I4/
     1 10X,29H-- NOMBER OF CONSTRAINTS (M):,I4/
     2 10X,38H-- NUMBER OF SLACK VARIABLES (IDIMC) :,I4//10X,
     3 32H+ + + + +   TABLE APLUS(M,N+M) :/)
  841 FORMAT (1X,12I10)
  843 FORMAT (/10X,34H-- SELECTION OF VECTOR TO PICK IN.)
  845 FORMAT (///25X,42H- - - - -   START OF ITERATION   - - - - -//
     1 10X,17H-- PHASE NUMBER :,I4/)
  850 FORMAT (//10X,34H+ + + + +   TABLE BASE0(M+2,M+2) : /)
  855 FORMAT (//10X,31H+ + + + +    TABLE IVARS(M+2) :/
     1 2(/6(1X,5HIVAR(,I4,3H) =,I5,:)))
  860 FORMAT (//10X,33H+ + + + +    TABLE BORNE(IDIMC) :/
     1 5(/6(1X,5HBORN(,I4,3H) =,1P,E11.3,:)))
  865 FORMAT (/10X,21H-- CANDIDATE VECTOR :,I4,2X,7HVALUE :,1P,E10.3)
  870 FORMAT (//10X,29H+ + + + +    TABLE CRE(N+M) :/
     1 3(/6(1X,4HCRE(,I4,3H) =,1P,E11.3,:)))
  875 FORMAT (/10X,19H-- SELECTED PIVOT :,I4)
  880 FORMAT (/10X,35H-- SELECTION OF VECTOR TO PICK OUT.)
  885 FORMAT(//10X,3HI =,I4,5X,3HJ =,I4,5X,5HXIJ =,1P,E11.3,5X,
     1 6HBBAR =,E11.3)
  890 FORMAT (//10X,28HOPTIMAL OBJECTIVE FUNCTION =,1P,E12.4//10X,
     1 28H+ + + + +    TABLE XOBJ(N) :/5(/6(1X,5HXOBJ(,I4,3H) =,1P,
     2 E11.3,:)))
  894 FORMAT (10X,27H-- SELECTED VARIABLES ATE :/15X,4HIQ =,I4,5X,
     1 4HQM =,1P,E11.3,10X,4HIP =,I4,5X,4HPM =,E11.3)
  895 FORMAT (10X,3HL =,I4,5X,4HLL =,I4,5X,10HIVARS(L) =,I5)
  905 FORMAT (10X,25H-- START OF PIVOTING; L =,I4,5X,8HXIK(L) =,
     1 1P,E11.3)
      END
