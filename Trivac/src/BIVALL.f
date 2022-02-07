*DECK BIVALL
      SUBROUTINE BIVALL (MAXPTS,IHEX,NH,NTH,ITAB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Unfold any hexagonal geometry to produce a complete domain.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Benaboud
*
*Parameters: input
* MAXPTS  maximum number of hexagons.
* IHEX    type of symmetry:
*         =1: S30;   =2: SA60;   =3: SB60;   =4: S90;   =5: R120;
*         =6: R180;  =7: SA180;  =8: SB180;  =9: COMPLETE.
* NH      total number of hexagons in the partial hexagonal geometry.
*
*Parameters: output
* NTH     total number of hexagons in the complete geometry.
* ITAB    correspondance table.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXPTS,IHEX,NH,NTH,ITAB(MAXPTS)
*----
*  LOCAL VARIABLES
*----
      LOGICAL LPAIR
      CHARACTER TEXT4*4
      INTEGER NP(7)
      INTEGER, DIMENSION(:), ALLOCATABLE :: J1,J2,J3,K1,K2,K3,K4
*
      NC=0
      IF((IHEX.EQ.1).OR.(IHEX.EQ.10)) THEN
         VI = 2.* SQRT(REAL(NH)) - 1.
         VP = SQRT(REAL(4*NH+1)) - 1.
         IF(AINT(VI).EQ.VI) THEN
            NC = INT(VI)
         ELSE IF(AINT(VP).EQ.VP) THEN
            NC = INT(VP)
         ELSE
            CALL XABORT('BIVALL: INVALID NUMBER OF HEXAGONS (1).')
         ENDIF
      ELSE IF((IHEX.EQ.2).OR.(IHEX.EQ.11)) THEN
         VA = (SQRT(REAL(8*NH+1)) - 1.)/2.
         IF(AINT(VA).EQ.VA) THEN
            NC = INT(VA)
         ELSE
            CALL XABORT('BIVALL: INVALID NUMBER OF HEXAGONS (2).')
         ENDIF
      ELSE IF(IHEX.EQ.3) THEN
         VI = SQRT(REAL(2*NH-1))
         VP = SQRT(REAL(2*NH))
         IF(AINT(VI).EQ.VI) THEN
            NC = INT(VI)
         ELSE IF(AINT(VP).EQ.VP) THEN
            NC = INT(VP)
         ELSE
            CALL XABORT('BIVALL: INVALID NUMBER OF HEXAGONS (3).')
         ENDIF
      ELSE IF(IHEX.EQ.4) THEN
         VI = SQRT(REAL((4*NH-1)/3))
         VP = SQRT(REAL(4*NH/3))
         IF(AINT(VI).EQ.VI) THEN
            NC = INT(VI)
         ELSE IF(AINT(VP).EQ.VP) THEN
            NC = INT(VP)
         ELSE
            CALL XABORT('BIVALL: INVALID NUMBER OF HEXAGONS (4).')
         ENDIF
      ELSE IF(IHEX.EQ.5) THEN
         VA = (SQRT(REAL(4*(NH-1)+1)) + 1.)/2.
         IF(AINT(VA).EQ.VA) THEN
            NC = INT(VA)
         ELSE
            CALL XABORT('BIVALL: INVALID NUMBER OF HEXAGONS (5).')
         ENDIF
      ELSE IF(IHEX.EQ.6) THEN
         VA = (SQRT(REAL(8*(NH-1)/3+1)) + 1)/2
         IF(AINT(VA).EQ.VA) THEN
            NC = INT(VA)
         ELSE
            CALL XABORT('BIVALL: INVALID NUMBER OF HEXAGONS (6).')
         ENDIF
      ELSE IF(IHEX.EQ.7) THEN
         VA = (SQRT(REAL(24*NH+1)) + 1.)/6.
         IF(AINT(VA).EQ.VA) THEN
            NC = INT(VA)
         ELSE
            CALL XABORT('BIVALL: INVALID NUMBER OF HEXAGONS (7).')
         ENDIF
      ELSE IF(IHEX.EQ.8) THEN
         VI = (1.+SQRT(REAL(3*(2*NH-1)+1)))/3.
         VP = (1.+SQRT(REAL(6*NH+1)))/3.
         IF(AINT(VI).EQ.VI) THEN
            NC = INT(VI)
         ELSE IF(AINT(VP).EQ.VP) THEN
            NC = INT(VP)
         ELSE
            CALL XABORT('BIVALL: INVALID NUMBER OF HEXAGONS (8).')
         ENDIF
      ELSE IF(IHEX.EQ.9) THEN
         VA = (SQRT(REAL((4*NH-1)/3)) + 1.)/2.
         IF(AINT(VA).EQ.VA) THEN
            NC = INT(VA)
         ELSE
            CALL XABORT('BIVALL: INVALID NUMBER OF HEXAGONS (9).')
         ENDIF
      ELSE
         WRITE(TEXT4,'(I4)') IHEX
         CALL XABORT('BIVALL: INVALID TYPE OF SYMMETRY (IHEX='//TEXT4//
     1   ').')
      ENDIF
      NTH = 1 + 3 * NC * (NC - 1)
      IF(NTH.GT.MAXPTS) CALL XABORT('BIVALL: MAXPTS OVERFLOW.')
      ITAB(1) = 1
      ALLOCATE(J1(NC+2),J2(NC+2),J3(NC+2),K1(NC+2),K2(NC+2),K3(NC+2),
     1 K4(NC+2))
      J1(1) = 1
      J2(1) = 1
      J3(1) = 1
      K1(1) = 1
      K2(1) = 1
      K3(1) = 1
      DO 10 L = 2,NC+1
         J1(L) = (L-1)*6
         J3(L) = 1+3*L*(L-1)
         J2(L) = 1+J3(L-1)
 10   CONTINUE
*
      IF((IHEX.EQ.1).OR.(IHEX.EQ.10)) THEN
         IL=0
         DO 20 L = 1,NC+1,2
            K1(L) = 1 + IL
            K1(L+1)   = 1 + IL
            IL = IL+1
 20      CONTINUE
         DO 30 L = 2,NC+1
            K2(L) = K2(L-1) + K1(L-1)
 30      CONTINUE
         IL=0
         DO 40 L = 1,NC+1,2
            K3(L) = K2(L) + IL
            K3(L+1)   = K2(L+1)   + IL
            IL = IL+1
 40      CONTINUE
      ELSE IF((IHEX.EQ.2).OR.(IHEX.EQ.11)) THEN
         K1(2) = 2
         DO 50 L = 2,NC+1
            K2(L) = K2(L-1) + L
            K1(L+1)   = K1(L) + L
 50      CONTINUE
      ELSE IF(IHEX.EQ.3) THEN
         K1(2) = 2
         DO 60 L = 1,NC+1
            K1(L+1) = K1(L) + L
 60      CONTINUE
         IL=0
         DO 70 L  = 1,NC+1,2
            K4(L) = 1 + IL
            K4(L+1)   = 1 + IL
            IL = IL + 2
 70      CONTINUE
         DO 80 L  = 2,NC+1
            K2(L) = K2(L-1) + K4(L-1)
            K3(L) = K3(L-1) + K4(L)
 80      CONTINUE
      ELSE IF(IHEX.EQ.4) THEN
         IL=0
         DO 90 L  = 1,NC+1,2
            K4(L) = L + IL
            K4(L+1)   = L + IL + 1
            IL = IL + 1
 90      CONTINUE
         DO 100 L  = 2,NC+1
            K1(L) = K1(L-1) + K4(L-1)
            K3(L) = K3(L-1) + K4(L)
100      CONTINUE
         IL=0
         DO 110 L = 1,NC+1,2
            K2(L) = K1(L) + IL
            K2(L+1)   = K1(L+1)   + IL
            IL = IL+1
110      CONTINUE
      ELSE IF(IHEX.EQ.5) THEN
         DO 120 L  = 2,NC+1
            K2(L) = 2 * (L-1)
            K1(L) = K1(L-1) + K2(L)
120      CONTINUE
      ELSE IF(IHEX.EQ.6) THEN
         DO 130 L  = 2,NC+1
            K2(L) = 3 * (L-1)
            K1(L) = K1(L-1) + K2(L)
130      CONTINUE
      ELSE IF(IHEX.EQ.7) THEN
         DO 140 L  = 2,NC+1
            K2(L) = 3 + K2(L-1)
            K1(L) = K1(L-1) + K2(L)
140      CONTINUE
      ELSE IF(IHEX.EQ.8) THEN
         IL = 1
         IF = 1
         DO 150 L  = 2,NC+1,2
            K2(L) = 3 * (L-1)
            K2(L+1)   = 3 * L + 1
150      CONTINUE
         DO 160 L  = 2,NC+1
            IL = IL + K2(L)
            IF = IF + K2(L-1)
            K1(L) = (IF + IL) / 2
160      CONTINUE
      ENDIF
*
      DO 300 N = 2,NTH
         I=0
         J=0
         DO 170 I0 = 2,NC
            IF((N.GE.J2(I0)).AND.(N.LE.J3(I0))) THEN
               I=I0
               GO TO 180
            ENDIF
170      CONTINUE
         IF(I.EQ.0) CALL XABORT('BIVALL: ALGORITHM FAILURE(1).')
180      DO 190 K = 1,6
            NP(K) = J2(I) + (K - 1) * (I - 1)
190      CONTINUE
         NP(7) = J3(I)
         COURS2 = REAL(I)/2.
         LPAIR = (AINT(COURS2).EQ.COURS2)
*
      IF((IHEX.EQ.1).OR.(IHEX.EQ.10)) THEN
         IF(N.LE.7) THEN
            ITAB(N) = 2
            GO TO 300
         ENDIF
         DO 200 L = 1,6
            IF((N.GE.NP(L)).AND.(N.LT.NP(L+1))) J = L
200      CONTINUE
         IF(N.EQ.NP(7)) J = 6
         IF(J.EQ.0) CALL XABORT('BIVALL: ALGORITHM FAILURE(2).')
         IC = 0
         IF(J.EQ.6) IC = 1
         N12 = (NP(J) + NP(J+1)+IC)/2
         N13 = N12 + 1
         IF(N.EQ.NP(J)) THEN
            ITAB(N) = K3(I)
         ELSE IF(N.EQ.NP(7)) THEN
            ITAB(N) = K3(I) - 1
         ELSE IF((N.GT.NP(J)).AND.(N.LT.N12)) THEN
            ITAB(N) = K3(I) - (N - NP(J))
         ELSE IF((N.EQ.N12).OR.((N.EQ.N13).AND.LPAIR)) THEN
            ITAB(N) = K2(I)
         ELSE IF((N.EQ.N13).AND.(.NOT.LPAIR)) THEN
            ITAB(N) = K3(I) - (NP(J+1) + IC - N)
         ELSE IF((N.GT.N13).AND.(N.LT.NP(J+1))) THEN
            ITAB(N) = K3(I) - (NP(J+1) + IC - N)
         ENDIF
*
      ELSE IF((IHEX.EQ.2).OR.(IHEX.EQ.11)) THEN
         DO 210 L = 1,6,2
            IF((N.GE.NP(L)).AND.(N.LT.NP(L+2))) J = L
210      CONTINUE
         IF(N.EQ.NP(7)) J = 5
         IF(J.EQ.0) CALL XABORT('BIVALL: ALGORITHM FAILURE(3).')
         IF(N.EQ.NP(J)) THEN
            ITAB(N) = K2(I)
         ELSE IF(N.EQ.NP(7)) THEN
            ITAB(N) = K2(I) - 1
         ELSE IF((N.GT.NP(J)).AND.(N.LT.NP(J+1))) THEN
            ITAB(N) = K2(I) - (N - NP(J))
         ELSE IF(N.EQ.NP(J+1)) THEN
            ITAB(N) = K1(I)
         ELSE IF((N.GT.NP(J+1)).AND.(N.LT.NP(J+2))) THEN
            ITAB(N) = K1(I) + (N - NP(J+1))
         ENDIF
*
      ELSE IF(IHEX.EQ.3) THEN
         IF(N.LE.7) THEN
            ITAB(N) = 2
            GO TO 300
         ENDIF
         DO 220 L = 1,6,2
            IF((N.GE.NP(L)).AND.(N.LT.NP(L+2))) J = L
220      CONTINUE
         IF(N.EQ.NP(7)) J = 5
         IF(J.EQ.0) CALL XABORT('BIVALL: ALGORITHM FAILURE(4).')
         IC = 0
         IF(J.EQ.5) IC = 1
         N12 = (NP(J) + NP(J+1))/2
         N13 = N12 + 1
         N14 = (NP(J+1) + NP(J+2)+IC)/2
         N15 = N14 + 1
         IF((N.EQ.NP(J)).OR.(N.EQ.NP(J+1))) THEN
            ITAB(N) = K1(I)
         ELSE IF(N.EQ.NP(7)) THEN
            ITAB(N) = K1(I) - 1
         ELSE IF((N.GT.NP(J)).AND.(N.LT.N12)) THEN
            ITAB(N) = K1(I) + (N - NP(J))
         ELSE IF((N.EQ.N12).OR.((N.EQ.N13).AND.LPAIR)) THEN
            ITAB(N) = K3(I)
         ELSE IF((N.EQ.N13).AND.(.NOT.LPAIR)) THEN
            ITAB(N) = K3(I) - 1
         ELSE IF((N.GT.N13).AND.(N.LT.NP(J+1))) THEN
            ITAB(N) = K1(I) + (NP(J+1) - N)
         ELSE IF((N.GT.NP(J+1)).AND.(N.LT.N14)) THEN
            ITAB(N) = K1(I) - (N - NP(J+1))
         ELSE IF((N.EQ.N14).OR.((N.EQ.N15).AND.LPAIR)) THEN
            ITAB(N) = K2(I)
         ELSE IF((N.EQ.N15).AND.(.NOT.LPAIR)) THEN
            ITAB(N) = K2(I) + 1
         ELSE IF((N.GT.N15).AND.(N.LT.NP(J+2))) THEN
            ITAB(N) = K1(I) - (NP(J+2) + IC - N)
         ENDIF
*
      ELSE IF(IHEX.EQ.4) THEN
         IF(N.EQ.7) THEN
            ITAB(N) = 2
            GO TO 300
         ENDIF
         DO 230 L = 1,6,3
            IF((N.GE.NP(L)).AND.(N.LT.NP(L+3))) J = L
230      CONTINUE
         IF(N.EQ.NP(7)) J = 4
         IF(J.EQ.0) CALL XABORT('BIVALL: ALGORITHM FAILURE(5).')
         IC = 0
         IF(J.EQ.4) IC = 1
         N12 = (NP(J+2) + NP(J+3)+IC)/2
         N13 = N12 + 1
         IF((N.EQ.NP(J)).OR.(N.EQ.NP(J+2))) THEN
            ITAB(N) = K2(I)
         ELSE IF(N.EQ.NP(7)) THEN
            ITAB(N) = K2(I) - 1
         ELSE IF((N.GT.NP(J)).AND.(N.LT.NP(J+1))) THEN
            ITAB(N) = K2(I) + (N - NP(J))
         ELSE IF(N.EQ.NP(J+1)) THEN
            ITAB(N) = K3(I)
         ELSE IF((N.GT.NP(J+1)).AND.(N.LE.N12).AND.(N.NE.NP(J+2))) THEN
            ITAB(N) = K2(I) - (N - NP(J+2))
         ELSE IF((N.EQ.N13).AND.(.NOT.LPAIR)) THEN
            ITAB(N) = K2(I) - (NP(J+3) + IC - N)
         ELSE IF((N.EQ.N13).AND.LPAIR) THEN
            ITAB(N) = K1(I)
         ELSE IF((N.GT.N13).AND.(N.LT.NP(J+3))) THEN
            ITAB(N) = K2(I) - (NP(J+3) + IC - N)
         ENDIF
*
      ELSE IF(IHEX.EQ.5) THEN
         IF(N.EQ.7) THEN
            ITAB(N) = 3
            GO TO 300
         ELSE IF((N.EQ.11).OR.(N.EQ.15).OR.(N.EQ.19)) THEN
            ITAB(N) = 4
            GO TO 300
         ENDIF
         DO 240 L = 1,6,2
            IF((N.GE.NP(L)).AND.(N.LT.NP(L+2))) J = L
240      CONTINUE
         IF(N.EQ.NP(7)) J = 5
         IF(J.EQ.0) CALL XABORT('BIVALL: ALGORITHM FAILURE(6).')
         IC = 0
         IF(J.EQ.5) IC = 1
         IF((N.GE.NP(J)).AND.(N.LE.NP(J+1))) THEN
            ITAB(N) = K1(I) - (NP(J+1) - N)
         ELSE IF((N.GT.NP(J+1)).AND.(N.LE.NP(J+2)+IC)) THEN
            ITAB(N) = K1(I) -(2*(NP(J+2)+IC)-NP(J+1)- N)
         ENDIF
*
      ELSE IF(IHEX.EQ.6) THEN
         DO 250 L = 1,6,3
            IF((N.GE.NP(L)).AND.(N.LT.NP(L+3))) J = L
250      CONTINUE
         IF(N.EQ.NP(7)) J = 4
         IF(J.EQ.0) CALL XABORT('BIVALL: ALGORITHM FAILURE(7).')
         IC = 0
         IF(J.EQ.4) IC = 1
         IF((N.GE.NP(J)).AND.(N.LE.NP(J+1))) THEN
            ITAB(N) = K1(I) - (NP(J+1) - N)
         ELSE IF((N.GT.NP(J+1)).AND.(N.LE.NP(J+2))) THEN
            ITAB(N) =  K1(I) - 2*(NP(J+2)-NP(J+1))-(NP(J+2)-N)
         ELSE IF((N.GT.NP(J+2)).AND.(N.LE.NP(J+3)+IC)) THEN
            ITAB(N) = K1(I)-(2*(NP(J+3)+IC)-NP(J+2)-N)
         ENDIF
*
      ELSE IF(IHEX.EQ.7) THEN
         IF((N.GE.NP(1)).AND.(N.LE.NP(2))) THEN
            ITAB(N) = K1(I) - (NP(2) - N)
         ELSE IF((N.GT.NP(2)).AND.(N.LE.NP(3))) THEN
            ITAB(N) = K1(I) - (NP(3) - NP(2)) + (NP(3) - N)
         ELSE IF((N.GT.NP(3)).AND.(N.LE.NP(4))) THEN
            ITAB(N) = K1(I) - (NP(4) - NP(2)) + (NP(4) - N)
         ELSE IF((N.GT.NP(4)).AND.(N.LE.NP(5))) THEN
            ITAB(N) = K1(I) - (NP(5) - NP(2)) + (NP(5) - N)
         ELSE IF((N.GT.NP(5)).AND.(N.LE.NP(6))) THEN
            ITAB(N) = K1(I) - (NP(4) - NP(2)) - (NP(6) - N)
         ELSE IF((N.GT.NP(6)).AND.(N.LE.NP(7)+1)) THEN
            ITAB(N) = K1(I) - (NP(3) - NP(2)) - (NP(7) + 1 - N)
         ENDIF
*
      ELSE IF(IHEX.EQ.8) THEN
         N12 = (NP(3) +  NP(4)) / 2
         N13 = (NP(6) +  NP(7) + 1) / 2
         IF((N.GE.NP(1)).AND.(N.LE.N12)) THEN
            ITAB(N) = K1(I) - (NP(2) - N)
         ELSE IF((N.GT.N12).AND.(N.LE.N13)) THEN
            ITAB(N) = K1(I) + (NP(5) - N)
         ELSE IF((N.GT.N13).AND.(N.LE.NP(7)+1)) THEN
            ITAB(N) = K1(I) - (NP(6) - NP(5)) - (NP(7) + 1 - N)
         ENDIF
*
      ELSE IF(IHEX.EQ.9) THEN
         ITAB(N) = N
      ENDIF
300   CONTINUE
      DEALLOCATE(K4,K3,K2,K1,J3,J2,J1)
      RETURN
      END
