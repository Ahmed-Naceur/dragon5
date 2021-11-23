*DECK TABEN
      FUNCTION TABEN(L,X)
*
*-----------------------------------------------------------------------
*
* COMPUTE AN EXPONENTIAL INTEGRAL OF ORDER L.
*
*----------------------------------------------------------- R. ROY ----
*
      IMPLICIT  NONE
      INTEGER   IOUT,NA,NB,NC
      CHARACTER NAMSBR*6
      PARAMETER(IOUT=6,NA=6,NB=4,NC=4,NAMSBR='TABEN ')
C----
C  ROUTINE PARAMETERS
C-----
      INTEGER   L
      REAL      X
C----
C  FUNCTION TYPE
C----
      REAL      TABEN
C----
C  LOCAL PARAMETERS
C----
      INTEGER   ITERM
      REAL      EX,P,R,S
C----
C  DATA
C----
      REAL      A(6),B(4),C(4)
      SAVE      A,B,C
      DATA      A
     >  /-.5772156649,.99999193,-.24991055,.05519968,
     >   -.00976004,.00107857/
      DATA      B
     >  /8.5733287401,18.0590169730,8.634760825,.2677737343/
      DATA      C
     >  /9.5733223454,25.6329561486,21.0996530827,3.9584969228/
      IF(X .LE. 0.)THEN
        IF(L .LE. 1)THEN
          TABEN= 1.0E20
        ELSE
          TABEN=1.0/REAL(L-1)
        ENDIF
      ELSE
        IF (X .LT. 50.) THEN
          EX=EXP(-X)
        ELSE
          EX=0.0
        ENDIF
        IF (L .EQ. 0) THEN
          TABEN=EX/X
        ELSE
          IF (X .LE. 1.0) THEN
            P=A(1)+X*(A(2)+X*(A(3)+X*(A(4)+X*(A(5)+X*A(6)))))
            TABEN=P-LOG(X)
          ELSE
            R=B(4)+X*(B(3)+X*(B(2)+X*(B(1)+X)))
            S=C(4)+X*(C(3)+X*(C(2)+X*(C(1)+X)))
            TABEN=R/S*EX/X
          ENDIF
          DO 100 ITERM=1,L-1
            TABEN=(EX-X*TABEN)/REAL(ITERM)
 100      CONTINUE
        ENDIF
      ENDIF
      RETURN
      END
