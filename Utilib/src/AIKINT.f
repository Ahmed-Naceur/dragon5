*DECK AIKINT
       FUNCTION AIKINT(Z,X,Y,N,EPS)
C----
C  REVISION HISTORY
C     - CREATED 80-AUG,  BY E.G. LONG
C     - CHECK THAT X-VALUES ARE MONOTONICALLY INCREASING,
C     - REVISED FROM AELIB 1985 OCT 16 BY P CARLSON
C     - CHANGED TO RELATIVE EPS AUGUST 1986 JV DONNELLY
C     - IMPLEMENTED IN DRAGON AUGUST 1993 (G. MARLEAU)
C  ABNORMAL TERMINATION
C      1. PARAMETER CHECKS FAIL-
C          A. X NOT MONOTONICALLY INCREASING
C      2. REQUESTED ACCURACY NOT OBTAINED
C  PARAMETERS
C     Z      : INTERPOLATION POINT
C     X      : INTERPOLATION POINT TABLULATION
C     Y      : TABULATED FUNCTION AT POINTS X
C     N      : NUMBER OF TABULATED POINTS
C     EPS    : TABULATION ERROR PERMITTED
C----
      SAVE       ICALLB,ICALLA
      PARAMETER (IOUT=6)
      LOGICAL    RELOK,ABSOK
      REAL       X(N),Y(N),W(10),F(10)
      DATA       ICALLB,ICALLA /0,0/
C----
C  CHECK FOR MONOTONICALLY INCREASING TABLE
C----
      IL=0
      MC=0
      IF(N.GE.2) THEN
        DO 120 I = 2, N
          IF (X(I) .LE. X(I-1)) THEN
            WRITE(IOUT,2000)
            CALL XABORT('AIKINT: ILLEGAL FORMAT TABLE')
          ENDIF
 120    CONTINUE
      ENDIF
C----
C  HANDLE INTEROPOLATION WITH N<=2
C----
      IF(N.LE.2) THEN
        IF(N.LE.0) AIKINT=0.0
        IF(N.EQ.1) AIKINT=Y(1)
        IF(N.EQ.2) AIKINT=((Z-X(2))*Y(1)-(Z-X(1))*Y(2))
     >                     /(X(1) - X(2))
C----
C  RETURN FROM AIKINT
C----
        RETURN
      ENDIF
C----
C  INTERPOLATION FOR N>2
C  CHECK IF Z=X(1)
C----
      IX=1
      IF(X(1).NE.Z) GO TO 20
C----
C Z IS AT ONE OF X(I)
C----
 21   CONTINUE
      AIKINT=Y(IX)
C----
C  RETURN FROM AIKINT
C----
      RETURN
C
C-----------------------------------------------------------
 20   CONTINUE
      IF(X(IX).GT.Z) THEN
C----
C  EXTRAPOLATION BEFORE TABLE BEGINS
C----
        IF(ICALLB.EQ.0) WRITE(IOUT,2001)
        ICALLB=ICALLB+1
        IU = 1
        JA = 1
        K = 2
        NA=MIN0(N,10)
        DO 40 I=1,NA
          W(I) = X(I)
          F(I) = Y(I)
 40     CONTINUE
        GO TO 600
      ELSE
C----
C  CHECK IF Z=X(N)
C----
        IX=N
        IF(X(IX).EQ.Z) GO TO 21
        IF(X(IX).LT.Z) THEN
C----
C  EXTRAPOLATION BEYOND END OF TABLE
C----
          IF(ICALLA.EQ.0) WRITE(IOUT,2002)
C          ICALLA=ICALLA+1
C          IL=N
C          JA=2
C          K=2
C          NA=MIN0(N,10)
C          DO 41 I=1,NA
C            J = N - I + 1
C            W(I) = X(J)
C            F(I) = Y(J)
C 41       CONTINUE
C          GO TO 600
           AIKINT=Y(N)
           RETURN
        ENDIF
      ENDIF
C----
C  Z IS WITHIN X(I), FIND WHERE Z LIES IN THE TABLE
C----
      IL = 1
      IU = N
      JA = 3
      A = N
      M = INT(ALOG(FLOAT(N))/0.693147) + 2
      DO 11 I = 1, M
        IR = IU - IL
        IF( IR .EQ. 1 ) GO TO 30
        IX = IL + IR/2
        IF( X(IX) .EQ. Z ) GO TO 21
        IF( X(IX) .GT. Z ) THEN
          IU = IX
        ELSE
          IL = IX
        ENDIF
 11   CONTINUE
      IU = N
      IL = N - 1
 30   CONTINUE
      K=0
      MC=3
C----
C  FINDING NEAREST ARGUMENTS TO Z
C  IF POSSIBLE, THE ARGUMENTS ARE CHOSEN IN PAIRS, SO THAT THEY ARE
C  ON THE SAME SIDE OF Z.  THE FIRST LINEAR CROSS MEANS IS
C  CALCULATED USING THE CLOSEST PAIR.  SUBSEQUENT LINEAR CROSS
C  MEANS ARE CALCULATED USING FIRST THE CLOSEST ARGUMENT TO Z OF
C  THE NEXT PAIR AND THEN THE OTHER ARGUMENT
C  OTHERWISE THE ARGUMENTS ARE CHOSEN IN ORDER OF CLOSENESS.
C----
      NA=MIN0(N,10)
 601  IF( K .EQ. NA ) GO TO 501
      I450=1
      I400=1
      IF( IU .EQ. (N+1) ) THEN
        I450=0
        I400=0
      ELSE IF( IL .EQ. 0) THEN
        I450=0
      ELSE IF(MC-2 .EQ. -1) THEN
        I450=0
        I400=0
      ELSE IF(MC-2 .EQ.  0) THEN
        I450=0
      ENDIF
* I450
      IF( I450 .EQ. 1) THEN
        MC=0
        D1=ABS(X(IU)-Z)
        D2=ABS(X(IL)-Z)
        IF( D1 .GT. D2 ) THEN
          I400=0
        ENDIF
      ENDIF
* I400
      IF( I400 .EQ. 1) THEN
        MC = MC + 1
        K = K + 1
        F(K) = Y(IU)
        W(K) = X(IU)
        IU = IU + 1
      ELSE
        MC = MC + 2
        K = K + 1
        F(K) = Y(IL)
        W(K) = X(IL)
        IL = IL - 1
      ENDIF
      IF(K .LT. 2) GO TO 601
C----
C  EVALUATION OF POSSIBLE ANSWERS
C----
 600  CONTINUE
      KA = K - 1
      DO 200 I= 1, KA
        F(K) = ( ( Z - W(K) )*F(I) -
     >         ( Z - W(I) )*F(K) )/( W(I) - W(K) )
 200  CONTINUE
C----
C  TEST FOR CONVERGENCE OF INTERPOLATION
C----
      IF( F(KA) .EQ. 0.0 ) THEN
        IF( ABS(F(K)) .LT. EPS ) THEN
          AIKINT = F(K)
          RETURN
        ENDIF
      ELSE
        IF( ABS( 1. - F(K)/F(KA) ) .LT. EPS) THEN
          AIKINT = F(K)
          RETURN
        ENDIF
      ENDIF
C----
C  NOT CONVERGED YET, TRY NEXT ORDER IF POSSIBLE
C----
      IF( JA .GT. 2 ) GO TO 601
      IF( K .EQ. NA ) GO TO 501
      K = K + 1
      GO TO 600
C----
C  REQUESTED ACCURACY WAS NOT OBTAINED
C  REVERT TO LEAST DIVERGENT INTERPOLATION
C  BASED ON RELATVE AND ABSOLUTE CONVERGENCE
C----
 501  CONTINUE
      KR = 1
      CONR = 100.
      KA = 1
      CONA = 100.*ABS(F(1))
      ABSOK = .TRUE.
      RELOK = .TRUE.
      DO 10 I = 2, NA
        IF( F(I-1) .NE. 0.0 ) THEN
          CON = ABS( 1.0 - F(I)/F(I-1) )
        ELSE
          CON = 1.0
        ENDIF
        IF( RELOK .AND. CON .LT. CONR ) THEN
          KR = I
          CONR = CON
        ELSE
          RELOK = .FALSE.
        ENDIF
        CON = ABS( F(I) - F(I-1) )
        IF( ABSOK .AND. CON .LT. CONA ) THEN
          KA = I
          CONA = CON
        ELSE
          ABSOK = .FALSE.
        ENDIF
 10   CONTINUE
      IF( KR .GT. KA ) THEN
        KK = KR
        CONMIN = CONR
      ELSE
        KK = KA
        CONMIN = CONA
        IF( F(KK) .NE. 0.0 ) CONMIN = CONMIN/F(KK)
      ENDIF
      AIKINT = F(KK)
      RETURN
C----
C  FORMAT
C----
 2000 FORMAT(' ---   ERROR IN AIKINT:  X-VALUES ARE',
     >       ' NOT MONOTONICALLY INCREASING  ---')
 2001 FORMAT(' ---   WARNING FROM AIKINT: EXTRAPOLATION',
     >       ' BEFORE TABLE BEGINS AT LEAST ONCE  ---')
 2002 FORMAT(' ---   WARNING FROM AIKINT: EXTRAPOLATION',
     >       ' BEYOND END OF TABLE  AT LEAST ONCE  ---')
      END
