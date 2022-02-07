*DECK NEIGH1
      SUBROUTINE NEIGH1 (NC,N,K,M,POIDS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the index of a neighbour hexagon for a given symmetry.
* The following SUBROUTINE points are available:
* NEIGH1: S30 symmetry;              NEIGH2: SA60 symmetry;
* NEIGH3: SB60 symmetry;             NEIGH4: S90 symmetry;
* NEIGH5: R120 symmetry;             NEIGH6: R180 symmetry;
* NEIGH7: SA180 symmetry;            NEIGH8: SB180 symmetry;
* NEIGH9: complete assembly;         NEIG10: S30 symmetry with HBC SYME;
* NEIG11: SA60 symmetry with HBC SYME.
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
* NC      total number of hexagonal crowns.
* N       index of the considered hexagon.
* K       index of the side.
* POIDS   weight of the hexagon.
*
*Parameters: output
* M       index of the neighbour hexagon (=n: reflection on side k;
*         .LT.0: axial symmetry or rotation).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   NC,N,K,M
      REAL      POIDS
*----
*  LOCAL VARIABLES
*----
      LOGICAL   EVEN
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBA,FSTA,LSTA
*
      ALLOCATE(NBA(NC+2),FSTA(NC+2),LSTA(NC+2))
      EVEN=.TRUE.
      CALL XDISET(NBA,NC+2,0)
      CALL XDISET(FSTA,NC+2,0)
      CALL XDISET(LSTA,NC+2,0)
      FSTA(1) = 1
      IL=0
      DO 1 L = 1,NC+1,2
         NBA(L)   = 1+IL
         NBA(L+1) = 1+IL
         IL = IL+1
 1    CONTINUE
      DO 2 L = 2,NC+1
         FSTA(L) = FSTA(L-1)+NBA(L-1)
 2    CONTINUE
      IL=0
      DO 3 L = 1,NC+1,2
         LSTA(L)   = FSTA(L)+IL
         LSTA(L+1) = FSTA(L+1)+IL
         IL = IL+1
 3    CONTINUE
*
      I=1
      IF (N.GT.1) THEN
         I=0
         DO 4 I0 = 1,NC
            IF ((N.GE.FSTA(I0)).AND.(N.LE.LSTA(I0))) THEN
               I=I0
               GO TO 5
            ENDIF
 4       CONTINUE
 5       IF (I.EQ.0) CALL XABORT('NEIGH1: ALGORITHM FAILURE.')
      ENDIF
*
      N1 = FSTA(I)
      N2 = LSTA(I)
      EVEN = MOD(I,2).EQ.0
*
      IF (K.EQ.1) THEN
*
         IF (N.EQ.1) THEN
            M = 2
         ELSE IF (EVEN) THEN
            M = N+NBA(I)+1
         ELSE IF (.NOT.EVEN) THEN
            M = N+NBA(I)
         ENDIF
*
      ELSE IF (K.EQ.2) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF (N.EQ.N2) THEN
            M = -(LSTA(I+1)-1)
         ELSE
            M =  N+1
         ENDIF
*
      ELSE IF (K.EQ.3) THEN
*
         IF ((N.EQ.1).OR.(N.EQ.2)) THEN
            M = -2
         ELSE IF (N.EQ.N2) THEN
            M = -(N-1)
         ELSE IF (EVEN) THEN
            M =  N-NBA(I-1)+1
         ELSE IF (.NOT.EVEN) THEN
            M =  N-NBA(I-1)
         ENDIF
*
      ELSE IF (K.EQ.4) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -FSTA(I-1)
         ELSE IF (EVEN) THEN
            M =  N-NBA(I-1)
         ELSE IF (.NOT.EVEN) THEN
            M =  N-NBA(I-1)-1
         ENDIF
*
      ELSE IF (K.EQ.5) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF ((N.EQ.N1).AND.EVEN) THEN
            M =  N
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -(N+1)
         ELSE
            M =  N-1
         ENDIF
*
      ELSE IF (K.EQ.6) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -FSTA(I+1)
         ELSE IF (EVEN) THEN
            M =  N+NBA(I)
         ELSE IF (.NOT.EVEN) THEN
            M =  N+NBA(I)-1
         ENDIF
*
      ENDIF
*
      IF (N.EQ.1) THEN
         POIDS = 1./12.
      ELSE IF (N.EQ.N2) THEN
         POIDS = 0.5
      ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
         POIDS = 0.5
      ELSE
         POIDS = 1.
      ENDIF
      DEALLOCATE(NBA,FSTA,LSTA)
      RETURN
      END
*
      SUBROUTINE NEIGH2 (NC,N,K,M,POIDS)
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   NC,N,K,M
      REAL      POIDS
*----
*  LOCAL VARIABLES
*----
      LOGICAL   EVEN
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBA,FSTA,LSTA
*
      ALLOCATE(NBA(NC+2),FSTA(NC+2),LSTA(NC+2))
      EVEN=.TRUE.
      CALL XDISET(NBA,NC+2,0)
      CALL XDISET(FSTA,NC+2,0)
      CALL XDISET(LSTA,NC+2,0)
      NBA(1)  = 1
      LSTA(1) = 1
      FSTA(1) = 1
      FSTA(2) = 2
      DO 7 L = 2,NC+1
         NBA(L)  = L
         LSTA(L) = L+LSTA(L-1)
         FSTA(L+1) = L+FSTA(L)
 7    CONTINUE
*
      I=0
      DO 8 I0 = 1,NC
         IF ((N.GE.FSTA(I0)).AND.(N.LE.LSTA(I0))) THEN
            I=I0
            GO TO 9
         ENDIF
 8    CONTINUE
      IF (I.EQ.0) CALL XABORT('NEIGH2: ALGORITHM FAILURE.')
*
 9    N1 = FSTA(I)
      N2 = LSTA(I)
*
      IF (K.EQ.1) THEN
*
         M =  N+NBA(I)+1
*
      ELSE IF (K.EQ.2) THEN
*
         IF (N.EQ.N2) THEN
            M = -(N+NBA(I))
         ELSE
            M =  N+1
         ENDIF
*
      ELSE IF (K.EQ.3) THEN
*
         IF (N.EQ.1) THEN
            M = -3
         ELSE IF (N.EQ.N2) THEN
            M = -(N-1)
         ELSE
            M =  N-NBA(I-1)
         ENDIF
*
      ELSE IF (K.EQ.4) THEN
*
         IF (N.EQ.N1) THEN
            M = -(N+1)
         ELSE
            M =  N-NBA(I-1)-1
         ENDIF
*
      ELSE IF (K.EQ.5) THEN
*
         IF (N.EQ.N1) THEN
            M = -(N+NBA(I)+1)
         ELSE
            M =  N-1
         ENDIF
*
      ELSE IF (K.EQ.6) THEN
*
         M =  N+NBA(I)
*
      ENDIF
*
      IF (N.EQ.1) THEN
         POIDS = 1./6.
      ELSE IF ((N.EQ.N1).OR.(N.EQ.N2)) THEN
         POIDS = 0.5
      ELSE
         POIDS = 1.
      ENDIF
      DEALLOCATE(NBA,FSTA,LSTA)
      RETURN
      END
*
      SUBROUTINE NEIGH3 (NC,N,K,M,POIDS)
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   NC,N,K,M
      REAL      POIDS
*----
*  LOCAL VARIABLES
*----
      LOGICAL   EVEN
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBA,FSTA,LSTA
*
      ALLOCATE(NBA(NC+2),FSTA(NC+2),LSTA(NC+2))
      EVEN=.TRUE.
      CALL XDISET(NBA,NC+2,0)
      CALL XDISET(FSTA,NC+2,0)
      CALL XDISET(LSTA,NC+2,0)
      LSTA(1) = 1
      FSTA(1) = 1
      IL=0
      DO 10 L = 1,NC+1,2
         NBA(L)   = 1+IL
         NBA(L+1) = 1+IL
         IL = IL+2
 10   CONTINUE
      DO 11 L = 2,NC+1
         FSTA(L) = FSTA(L-1)+NBA(L-1)
        LSTA(L) = NBA(L)+LSTA(L-1)
 11   CONTINUE
*
      I=0
      N1=0
      N2=0
      N3=0
      IF (N.EQ.1) GOTO 14
      DO 12 I0 = 1,NC
         IF ((N.GE.FSTA(I0)).AND.(N.LE.LSTA(I0))) THEN
            I=I0
            GO TO 13
         ENDIF
 12   CONTINUE
      IF (I.EQ.0) CALL XABORT('NEIGH3: ALGORITHM FAILURE.')
*
 13   N1 = FSTA(I)
      N2 = (FSTA(I)+LSTA(I))/2
      N3 = LSTA(I)
      EVEN = MOD(I,2).EQ.0
*
 14   IF (K.EQ.1) THEN
*
         IF (N.EQ.1) THEN
            M =  2
         ELSE IF ((N.GE.N1).AND.(N.LE.N3)) THEN
            M =  N+I
         ENDIF
*
      ELSE IF (K.EQ.2) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF (N.EQ.2) THEN
            M =  5
         ELSE IF ((N.GE.N1).AND.(N.LT.N2)) THEN
            M =  N+1
         ELSE IF ((N.GE.N2).AND.(N.LT.N3)) THEN
            M =  N+I+1
         ELSE IF ((N.EQ.N3).AND.EVEN) THEN
            M =  LSTA(I+1)
         ELSE IF ((N.EQ.N3).AND.(.NOT.EVEN)) THEN
            M = -LSTA(I+1)
         ENDIF
*
      ELSE IF (K.EQ.3) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF (N.EQ.2) THEN
            M =  2
         ELSE IF ((N.GE.N1).AND.(N.LT.N2)) THEN
            M =  N+2-I
         ELSE IF ((N.GE.N2).AND.(N.LT.N3)) THEN
            M =  N+1
         ELSE IF ((N.EQ.N3).AND.EVEN) THEN
            M =  N
         ELSE IF ((N.EQ.N3).AND.(.NOT.EVEN)) THEN
            M = -(N-1)
         ENDIF
*
      ELSE IF (K.EQ.4) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF ((N.EQ.N1).AND.EVEN) THEN
            M =  N+1-I
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -(N+2-I)
         ELSE IF ((N.GT.N1).AND.(N.LT.N3)) THEN
            M =  N+1-I
         ELSE IF ((N.EQ.N3).AND.EVEN) THEN
            M =  N+1-I
         ELSE IF ((N.EQ.N3).AND.(.NOT.EVEN)) THEN
            M = -(N-I)
         ENDIF
*
      ELSE IF (K.EQ.5) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF ((N.EQ.N1).AND.EVEN) THEN
            M =  N
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -(N+1)
         ELSE IF ((N.GT.N1).AND.(N.LE.N2)) THEN
            M =  N-1
         ELSE IF ((N.GT.N2).AND.(N.LE.N3)) THEN
            M =  N-I
         ENDIF
*
      ELSE IF (K.EQ.6) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF ((N.EQ.N1).AND.EVEN) THEN
            M =  N+I-1
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -(N+I)
         ELSE IF ((N.GT.N1).AND.(N.LE.N2)) THEN
            M =  N+I-1
         ELSE IF ((N.GT.N2).AND.(N.LE.N3)) THEN
            M =  N-1
         ENDIF
*
      ENDIF
*
      IF (N.EQ.1) THEN
         POIDS = 1./6.
      ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
         POIDS = 0.5
      ELSE IF ((N.EQ.N3).AND.(.NOT.EVEN)) THEN
         POIDS = 0.5
      ELSE
         POIDS = 1.
      ENDIF
      DEALLOCATE(NBA,FSTA,LSTA)
      RETURN
      END
*
      SUBROUTINE NEIGH4 (NC,N,K,M,POIDS)
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   NC,N,K,M
      REAL      POIDS
*----
*  LOCAL VARIABLES
*----
      LOGICAL   EVEN
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBA,FSTA,LSTA,INTA
*
      ALLOCATE(NBA(NC+2),FSTA(NC+2),LSTA(NC+2),INTA(NC+2))
      EVEN=.TRUE.
      CALL XDISET(NBA,NC+2,0)
      CALL XDISET(FSTA,NC+2,0)
      CALL XDISET(LSTA,NC+2,0)
      CALL XDISET(INTA,NC+2,0)
      LSTA(1) = 1
      FSTA(1) = 1
      IL=0
      DO 15 L = 1,NC+1,2
         NBA(L)   = L+IL
         NBA(L+1) = L+1+IL
         IL = IL+1
 15   CONTINUE
      DO 16 L = 2,NC+1
         FSTA(L) = FSTA(L-1)+NBA(L-1)
         LSTA(L) = NBA(L)+LSTA(L-1)
 16   CONTINUE
      IL=0
      DO 17 L = 1,NC+1,2
         INTA(L)   = FSTA(L)+IL
         INTA(L+1) = FSTA(L+1)+IL
         IL = IL+1
 17   CONTINUE
*
      I=0
      N1=0
      N2=0
      N3=0
      IF (N.EQ.1) GOTO 20
      DO 18 I0 = 1,NC
         IF ((N.GE.FSTA(I0)).AND.(N.LE.LSTA(I0))) THEN
            I=I0
            GO TO 19
         ENDIF
 18   CONTINUE
      IF (I.EQ.0) CALL XABORT('NEIGH4: ALGORITHM FAILURE.')
*
 19   N1 = FSTA(I)
      N2 = INTA(I)
      N3 = LSTA(I)
      EVEN = MOD(I,2).EQ.0
*
 20   IF (K.EQ.1) THEN
*
         IF (N.EQ.1) THEN
            M =  2
         ELSE IF ((N.GE.N1).AND.(N.LE.N3).AND.EVEN) THEN
            M =  N+NBA(I)+1
         ELSE IF ((N.GE.N1).AND.(N.LE.N3).AND.(.NOT.EVEN)) THEN
            M =  N+NBA(I)
         ENDIF
*
      ELSE IF (K.EQ.2) THEN
*
         IF (N.EQ.1) THEN
            M =  3
         ELSE IF (N.EQ.2) THEN
            M =  6
         ELSE IF (N.EQ.N1) THEN
            M =  N+1
         ELSE IF ((N.GT.N1).AND.(N.LT.N2)) THEN
            M =  N+1
         ELSE IF ((N.GE.N2).AND.(N.LE.N3)) THEN
            M =  N+NBA(I+1)
         ENDIF
*
      ELSE IF (K.EQ.3) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF (N.EQ.2) THEN
            M =  3
         ELSE IF ((N.EQ.N1).AND.EVEN) THEN
            M =  FSTA(I-1)+1
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M =  FSTA(I-1)
         ELSE IF ((N.GT.N1).AND.(N.LT.N2).AND.EVEN) THEN
            M =  N-(I-1)-(INTA(I-1)-FSTA(I-1))+1
         ELSE IF ((N.GT.N1).AND.(N.LT.N2).AND.(.NOT.EVEN)) THEN
            M =  N-(I-1)-(INTA(I-1)-FSTA(I-1))
         ELSE IF ((N.GE.N2).AND.(N.LT.N3)) THEN
            M =  N+1
         ELSE IF (N.EQ.N3) THEN
            M = -(LSTA(I+1)-1)
         ENDIF
*
      ELSE IF (K.EQ.4) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF ((N.EQ.N1).AND.EVEN) THEN
            M =  FSTA(I-1)
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -FSTA(I-1)
         ELSE IF ((N.GT.N1).AND.(N.LT.N3).AND.EVEN) THEN
            M =  N-NBA(I-1)
         ELSE IF ((N.GT.N1).AND.(N.LT.N3).AND.(.NOT.EVEN)) THEN
            M =  N-NBA(I-1)-1
         ELSE IF (N.EQ.N3) THEN
            M = -(N-1)
         ENDIF
*
      ELSE IF (K.EQ.5) THEN
*
         IF (N.EQ.1) THEN
            M = -3
         ELSE IF ((N.EQ.N1).AND.EVEN) THEN
            M =  N
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -(N+1)
         ELSE IF ((N.GT.N1).AND.(N.LE.N2)) THEN
            M =  N-1
         ELSE IF ((N.GT.N2).AND.(N.LE.N3)) THEN
            M =  N-NBA(I)
         ENDIF
*
      ELSE IF (K.EQ.6) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF ((N.EQ.N1).AND.EVEN) THEN
            M =  FSTA(I+1)
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -FSTA(I+1)
         ELSE IF ((N.GT.N1).AND.(N.LT.N2).AND.EVEN) THEN
            M =  N+I+INTA(I)-FSTA(I)
         ELSE IF ((N.GT.N1).AND.(N.LT.N2).AND.(.NOT.EVEN)) THEN
            M =  N+I+INTA(I)-FSTA(I)-1
         ELSE IF (N.EQ.N2) THEN
            M =  INTA(I+1)-1
         ELSE IF ((N.GT.N2).AND.(N.LE.N3)) THEN
            M =  N-1
         ENDIF
*
      ENDIF
*
      IF (N.EQ.1) THEN
         POIDS = 0.25
      ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
         POIDS = 0.5
      ELSE IF (N.EQ.N3) THEN
         POIDS = 0.5
      ELSE
         POIDS = 1.
      ENDIF
      DEALLOCATE(NBA,FSTA,LSTA,INTA)
      RETURN
      END
*
      SUBROUTINE NEIGH5 (NC,N,K,M,POIDS)
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   NC,N,K,M
      REAL      POIDS
*----
*  LOCAL VARIABLES
*----
      LOGICAL   EVEN
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBA,FSTA,LSTA
*
      ALLOCATE(NBA(NC+2),FSTA(NC+2),LSTA(NC+2))
      EVEN=.TRUE.
      CALL XDISET(NBA,NC+2,0)
      CALL XDISET(FSTA,NC+2,0)
      CALL XDISET(LSTA,NC+2,0)
      NBA (1) = 1
      LSTA(1) = 1
      FSTA(1) = 1
      DO 21 L = 2,NC+1
         NBA(L)  = 2*(L-1)
         LSTA(L) = NBA(L)+LSTA(L-1)
         FSTA(L) = NBA(L-1)+FSTA(L-1)
 21   CONTINUE
*
      I=0
      DO 22 I0 = 1,NC
         IF ((N.GE.FSTA(I0)).AND.(N.LE.LSTA(I0))) THEN
            I=I0
            GO TO 23
         ENDIF
 22   CONTINUE
      IF (I.EQ.0) CALL XABORT('NEIGH5: ALGORITHM FAILURE.')
*
 23   N1 = FSTA(I)
      N2 = FSTA(I) + (I-2)
      N3 = LSTA(I)
*
      IF (K.EQ.1) THEN
*
         IF (N.EQ.1) THEN
            M =  2
         ELSE IF ((N.GE.N1).AND.(N.LT.N2).AND.(N.NE.1)) THEN
            M =  N+NBA(I+1)-1
         ELSE IF ((N.GE.N2).AND.(N.LT.N3)) THEN
            M =  N+NBA(I)+1
         ELSE IF (N.EQ.N3) THEN
            M =  N+NBA(I+1)-1
         ENDIF
*
      ELSE IF (K.EQ.2) THEN
*
         IF (N.EQ.1) THEN
            M =  3
         ELSE IF (N.EQ.2) THEN
            M =  6
         ELSE IF ((N.GE.N1).AND.(N.LT.N2)) THEN
            M =  N+1
         ELSE IF ((N.GE.N2).AND.(N.LT.N3)) THEN
            M =  N+NBA(I)+2
         ELSE IF (N.EQ.N3) THEN
            M =  N+NBA(I+1)
         ENDIF
*
      ELSE IF (K.EQ.3) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF (N.EQ.2) THEN
            M =  3
         ELSE IF ((N.GE.N1).AND.(N.LT.N2).AND.(N.NE.1)) THEN
            M =  N-NBA(I-1)
         ELSE IF ((N.GE.N2).AND.(N.LT.N3)) THEN
            M =  N+1
         ELSE IF (N.EQ.N3) THEN
            M = -(N+1)
         ENDIF
*
      ELSE IF (K.EQ.4) THEN
*
         IF (N.EQ.1) THEN
            M = -3
         ELSE IF (N.EQ.2) THEN
            M =  1
         ELSE IF (N.EQ.3) THEN
            M = -2
         ELSE IF ((N.EQ.N1).AND.(N.NE.1)) THEN
            M = -LSTA(I-1)
         ELSE IF ((N.GT.N1).AND.(N.LT.N3)) THEN
            M =  N-NBA(I-1)-1
         ELSE IF ((N.EQ.N3).AND.(N.NE.3)) THEN
            M = -(N-NBA(I-1)-1)
         ENDIF
*
      ELSE IF (K.EQ.5) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF (N.EQ.2) THEN
            M = -3
         ELSE IF ((N.EQ.N1).AND.(N.NE.1)) THEN
            M = -LSTA(I)
         ELSE IF ((N.GT.N1).AND.(N.LE.N2).AND.(N.NE.3)) THEN
            M =  N-1
         ELSE IF ((N.GT.N2).AND.(N.LT.N3)) THEN
            M =  N-NBA(I-1)-2
         ELSE IF (N.EQ.N3) THEN
            M =  N-NBA(I)
         ENDIF
*
      ELSE IF (K.EQ.6) THEN
*
         IF (N.EQ.1) THEN
            M = -3
         ELSE IF (N.EQ.3) THEN
            M =  2
         ELSE IF ((N.EQ.N1).AND.(N.NE.1)) THEN
            M =  FSTA(I+1)
         ELSE IF ((N.GT.N1).AND.(N.LT.N2)) THEN
            M =  N+NBA(I+1)-2
         ELSE IF (N.EQ.N2) THEN
            M =  N+NBA(I)
         ELSE IF ((N.GT.N2).AND.(N.LE.N3)) THEN
            M =  N-1
         ENDIF
*
      ENDIF
*
      IF (N.EQ.1) THEN
         POIDS = 1./3.
      ELSE
         POIDS = 1.
      ENDIF
      DEALLOCATE(NBA,FSTA,LSTA)
      RETURN
      END
*
      SUBROUTINE NEIGH6 (NC,N,K,M,POIDS)
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   NC,N,K,M
      REAL      POIDS
*----
*  LOCAL VARIABLES
*----
      LOGICAL   EVEN
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBA,FSTA,LSTA
*
      ALLOCATE(NBA(NC+2),FSTA(NC+2),LSTA(NC+2))
      EVEN=.TRUE.
      CALL XDISET(NBA,NC+2,0)
      CALL XDISET(FSTA,NC+2,0)
      CALL XDISET(LSTA,NC+2,0)
      NBA (1) = 1
      LSTA(1) = 1
      FSTA(1) = 1
      DO 24 L = 2,NC+1
         NBA(L)  = 3*(L-1)
         LSTA(L) = NBA(L)+LSTA(L-1)
         FSTA(L) = NBA(L-1)+FSTA(L-1)
 24   CONTINUE
*
      I=0
      N1=0
      N2=0
      N3=0
      N4=0
      IF (N.EQ.1) GOTO 27
      DO 25 I0 = 1,NC
         IF ((N.GE.FSTA(I0)).AND.(N.LE.LSTA(I0))) THEN
            I=I0
            GO TO 26
         ENDIF
 25   CONTINUE
      IF (I.EQ.0) CALL XABORT('NEIGH6: ALGORITHM FAILURE.')
*
 26   N1 = FSTA(I)
      N2 = LSTA(I-1) +   (I-1)
      N3 = LSTA(I-1) + 2*(I-1)
      N4 = LSTA(I)
*
 27   IF (K.EQ.1) THEN
*
         IF (N.EQ.1) THEN
            M =  3
         ELSE IF (N.EQ.2) THEN
            M =  7
         ELSE IF ((N.GE.N1).AND.(N.LT.N2)) THEN
            M =  N+1
         ELSE IF ((N.GE.N2).AND.(N.LE.N4)) THEN
            M =  N+NBA(I)+2
         ENDIF
*
      ELSE IF (K.EQ.2) THEN
*
         IF (N.EQ.1) THEN
            M =  4
         ELSE IF (N.EQ.2) THEN
            M =  3
         ELSE IF ((N.GE.N1).AND.(N.LT.N2)) THEN
            M =  N-NBA(I-1)
         ELSE IF ((N.GE.N2).AND.(N.LT.N3)) THEN
            M =  N+1
         ELSE IF ((N.GE.N3).AND.(N.LE.N4)) THEN
            M =  N+NBA(I+1)
         ENDIF
*
      ELSE IF (K.EQ.3) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF (N.EQ.2) THEN
            M =  1
         ELSE IF (N.EQ.N1) THEN
            M = -(N-1)
         ELSE IF ((N.GT.N1).AND.(N.LT.N3)) THEN
            M =  N-NBA(I-1)-1
         ELSE IF ((N.GE.N3).AND.(N.LT.N4)) THEN
            M =  N+1
         ELSE IF (N.EQ.N4) THEN
            M = -(N+1)
         ENDIF
*
      ELSE IF (K.EQ.4) THEN
*
         IF (N.EQ.1) THEN
            M = -3
         ELSE IF (N.EQ.2) THEN
            M = -4
         ELSE IF (N.EQ.3) THEN
            M =  1
         ELSE IF (N.EQ.N1) THEN
            M = -LSTA(I)
         ELSE IF ((N.GT.N1).AND.(N.LE.N2)) THEN
            M =  N-1
         ELSE IF ((N.GT.N2).AND.(N.LT.N4)) THEN
            M =  N-NBA(I)+1
         ELSE IF (N.EQ.N4) THEN
            M = -(N-NBA(I)+1)
         ENDIF
*
      ELSE IF (K.EQ.5) THEN
*
         IF (N.EQ.1) THEN
            M = -4
         ELSE IF (N.EQ.3) THEN
            M =  2
         ELSE IF ((N.GE.N1).AND.(N.LE.N2)) THEN
            M =  N+NBA(I)
         ELSE IF ((N.GT.N2).AND.(N.LE.N3)) THEN
            M =  N-1
         ELSE IF ((N.GT.N3).AND.(N.LE.N4)) THEN
            M =  N-NBA(I)
         ENDIF
*
      ELSE IF (K.EQ.6) THEN
*
         IF (N.EQ.1) THEN
            M =  2
         ELSE IF ((N.GE.N1).AND.(N.LE.N3)) THEN
            M =  N+NBA(I)+1
         ELSE IF ((N.GT.N3).AND.(N.LE.N4)) THEN
            M =  N-1
         ENDIF
*
      ENDIF
*
      IF (N.EQ.1) THEN
         POIDS = 1./2.
      ELSE
         POIDS = 1.
      ENDIF
      DEALLOCATE(NBA,FSTA,LSTA)
      RETURN
      END
*
      SUBROUTINE NEIGH7 (NC,N,K,M,POIDS)
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   NC,N,K,M
      REAL      POIDS
*----
*  LOCAL VARIABLES
*----
      LOGICAL   EVEN
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBA,FSTA,LSTA
*
      ALLOCATE(NBA(NC+2),FSTA(NC+2),LSTA(NC+2))
      EVEN=.TRUE.
      CALL XDISET(NBA,NC+2,0)
      CALL XDISET(FSTA,NC+2,0)
      CALL XDISET(LSTA,NC+2,0)
      NBA (1) = 1
      LSTA(1) = 1
      FSTA(1) = 1
      DO 28 L = 2,NC+1
         NBA(L)  = 3+NBA(L-1)
         LSTA(L) = NBA(L)+LSTA(L-1)
         FSTA(L) = 1+LSTA(L-1)
 28   CONTINUE
*
      I=0
      DO 29 I0 = 1,NC
         IF ((N.GE.FSTA(I0)).AND.(N.LE.LSTA(I0))) THEN
            I=I0
            GO TO 30
         ENDIF
 29   CONTINUE
      IF (I.EQ.0) CALL XABORT('NEIGH7: ALGORITHM FAILURE.')
*
 30   N1 = FSTA(I)
      N2 = FSTA(I) + (I-1)
      N3 = FSTA(I) + 2*(I-1)
      N4 = LSTA(I)
*
      IF (K.EQ.1) THEN
*
         IF (N.EQ.1) THEN
            M =  4
         ELSE IF ((N.GE.N1).AND.(N.LT.N2).AND.(N.NE.1)) THEN
            M =  N+1
         ELSE IF ((N.GE.N2).AND.(N.LT.N3)) THEN
            M =  NBA(I)+N+2
         ELSE IF ((N.GE.N3).AND.(N.LE.N4)) THEN
            M =  NBA(I+1)+N-1
         ENDIF
*
      ELSE IF (K.EQ.2) THEN
*
         IF (N.EQ.1) THEN
            M =  5
         ELSE IF ((N.GE.N1).AND.(N.LT.N2).AND.(N.NE.1)) THEN
            M =  N-NBA(I-1)
         ELSE IF ((N.GE.N2).AND.(N.LT.N3)) THEN
            M =  N+1
         ELSE IF ((N.GE.N3).AND.(N.LE.N4)) THEN
            M =  N+NBA(I+1)
         ENDIF
*
      ELSE IF (K.EQ.3) THEN
*
         IF (N.EQ.1) THEN
            M = -4
         ELSE IF ((N.EQ.N1).AND.(N.NE.1)) THEN
            M = -(N+1)
         ELSE IF ((N.GT.N1).AND.(N.LT.N3)) THEN
            M =  N-NBA(I-1)-1
         ELSE IF ((N.GE.N3).AND.(N.LT.N4)) THEN
            M =  N+1
         ELSE IF (N.EQ.N4) THEN
            M = -(N+NBA(I+1)-1)
         ENDIF
*
      ELSE IF (K.EQ.4) THEN
*
         IF (N.EQ.1) THEN
            M = -3
         ELSE IF ((N.EQ.N1).AND.(N.NE.1)) THEN
            M = -(FSTA(I+1)+1)
         ELSE IF ((N.GT.N1).AND.(N.LE.N2)) THEN
            M =  N-1
         ELSE IF ((N.GT.N2).AND.(N.LT.N3)) THEN
            M =  N-NBA(I-1)-2
         ELSE IF ((N.GE.N3).AND.(N.LT.N4)) THEN
            M =  N-NBA(I)+1
         ELSE IF (N.EQ.N4) THEN
            M = -(N-1)
         ENDIF
*
      ELSE IF (K.EQ.5) THEN
*
         IF (N.EQ.N1) THEN
            M =  FSTA(I+1)
         ELSE IF ((N.GT.N1).AND.(N.LE.N2)) THEN
            M =  N+NBA(I)
         ELSE IF ((N.GT.N2).AND.(N.LE.N3)) THEN
            M =  N-1
         ELSE IF ((N.GT.N3).AND.(N.LE.N4)) THEN
            M =  N-NBA(I)
         ENDIF
*
      ELSE IF (K.EQ.6) THEN
*
         IF (N.EQ.N1) THEN
            M =  FSTA(I+1)+1
         ELSE IF ((N.GT.N1).AND.(N.LT.N3)) THEN
            M =  N+NBA(I)+1
         ELSE IF (N.EQ.N3) THEN
            M =  N+NBA(I+1)-2
         ELSE IF ((N.GT.N3).AND.(N.LE.N4)) THEN
            M =  N-1
         ENDIF
*
      ENDIF
*
      IF ((N.EQ.N1).OR.(N.EQ.N4)) THEN
         POIDS = 0.5
      ELSE
         POIDS = 1.
      ENDIF
      DEALLOCATE(NBA,FSTA,LSTA)
      RETURN
      END
*
      SUBROUTINE NEIGH8 (NC,N,K,M,POIDS)
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   NC,N,K,M
      REAL      POIDS
*----
*  LOCAL VARIABLES
*----
      LOGICAL   EVEN
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBA,FSTA,LSTA
*
      ALLOCATE(NBA(NC+2),FSTA(NC+2),LSTA(NC+2))
      EVEN=.TRUE.
      CALL XDISET(NBA,NC+2,0)
      CALL XDISET(FSTA,NC+2,0)
      CALL XDISET(LSTA,NC+2,0)
      NBA (1) = 1
      LSTA(1) = 1
      FSTA(1) = 1
      DO 31 L = 2,NC+1,2
         NBA(L)  = 3*(L-1)
         NBA(L+1)  = 3*L+1
 31   CONTINUE
      DO 32 L = 2,NC+1
         LSTA(L) = NBA(L)+LSTA(L-1)
         FSTA(L) = NBA(L-1)+FSTA(L-1)
 32   CONTINUE
*
      I=0
      N1=0
      N2=0
      N3=0
      N4=0
      N5=0
      IF (N.EQ.1) GOTO 35
      DO 33 I0 = 1,NC
         IF ((N.GE.FSTA(I0)).AND.(N.LE.LSTA(I0))) THEN
            I=I0
            GO TO 34
         ENDIF
 33   CONTINUE
      IF (I.EQ.0) CALL XABORT('NEIGH8: ALGORITHM FAILURE.')
*
 34   N1 = FSTA(I)
      N2 = (FSTA(I) + LSTA(I))/2 - (I-1)
      N3 = (FSTA(I) + LSTA(I))/2
      N4 = (FSTA(I) + LSTA(I))/2 + (I-1)
      N5 = LSTA(I)
      EVEN = MOD(I,2).EQ.0
*
 35   IF (K.EQ.1) THEN
*
         IF (N.EQ.1) THEN
            M =  2
         ELSE IF ((N.GE.N1).AND.(N.LE.N3)) THEN
            M =  N+3*I-2
         ELSE IF ((N.GT.N3).AND.(N.LE.N4)) THEN
            M =  N-1
         ELSE IF ((N.GT.N4).AND.(N.LE.N5)) THEN
            M =  N-(3*I-2)
         ENDIF
*
      ELSE IF (K.EQ.2) THEN
*
         IF (N.EQ.1) THEN
            M =  3
         ELSE IF (N.EQ.2) THEN
            M =  7
         ELSE IF ((N.GE.N1).AND.(N.LT.N2)) THEN
            M =  N+1
         ELSE IF ((N.GE.N2).AND.(N.LE.N4)) THEN
            M =  N+3*I-1
         ELSE IF ((N.GT.N4).AND.(N.LE.N5)) THEN
            M =  N-1
         ENDIF
*
      ELSE IF (K.EQ.3) THEN
*
         IF (N.EQ.1) THEN
            M =  4
         ELSE IF (N.EQ.2) THEN
            M =  3
         ELSE IF ((N.GE.N1).AND.(N.LT.N2)) THEN
            M =  N-3*(I-2)
         ELSE IF ((N.GE.N2).AND.(N.LT.N3)) THEN
            M =  N+1
         ELSE IF ((N.GE.N3).AND.(N.LE.N5)) THEN
            M =  N+3*I
         ENDIF
*
      ELSE IF (K.EQ.4) THEN
*
         IF (N.EQ.1) THEN
            M = -4
         ELSE IF ((N.EQ.N1).AND.EVEN) THEN
            M =  FSTA(I-1)
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -FSTA(I-1)
         ELSE IF ((N.GT.N1).AND.(N.LT.N3)) THEN
            M =  N-3*I+5
         ELSE IF ((N.GE.N3).AND.(N.LT.N4)) THEN
            M =  N+1
         ELSE IF ((N.GE.N4).AND.(N.LT.N5)) THEN
            M =  N+3*I+1
         ELSE IF ((N.EQ.N5).AND.EVEN) THEN
            M =  LSTA(I+1)
         ELSE IF ((N.EQ.N5).AND.(.NOT.EVEN)) THEN
            M = -LSTA(I+1)
         ENDIF
*
      ELSE IF (K.EQ.5) THEN
*
         IF (N.EQ.1) THEN
            M = -3
         ELSE IF ((N.EQ.N1).AND.EVEN) THEN
            M =  N
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -(N+1)
         ELSE IF ((N.GT.N1).AND.(N.LE.N2)) THEN
            M =  N-1
         ELSE IF ((N.GT.N2).AND.(N.LT.N4)) THEN
            M =  N-3*I+4
         ELSE IF ((N.GE.N4).AND.(N.LT.N5)) THEN
            M =  N+1
         ELSE IF ((N.EQ.N5).AND.EVEN) THEN
            M =  N
         ELSE IF ((N.EQ.N5).AND.(.NOT.EVEN)) THEN
            M = -(N-1)
         ENDIF
*
      ELSE IF (K.EQ.6) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF ((N.EQ.N1).AND.EVEN) THEN
            M =  FSTA(I+1)
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -FSTA(I+1)
         ELSE IF ((N.GT.N1).AND.(N.LE.N2)) THEN
            M =  N+3*(I-1)
         ELSE IF ((N.GT.N2).AND.(N.LE.N3)) THEN
            M =  N-1
         ELSE IF ((N.GT.N3).AND.(N.LT.N5)) THEN
            M =  N-3*(I-1)
         ELSE IF ((N.EQ.N5).AND.EVEN) THEN
            M =  LSTA(I-1)
         ELSE IF ((N.EQ.N5).AND.(.NOT.EVEN)) THEN
            M = -LSTA(I-1)
         ENDIF
*
      ENDIF
*
      IF (N.EQ.1) THEN
         POIDS = 0.5
      ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
         POIDS = 0.5
      ELSE IF ((N.EQ.N5).AND.(.NOT.EVEN)) THEN
         POIDS = 0.5
      ELSE
         POIDS = 1.
      ENDIF
      DEALLOCATE(NBA,FSTA,LSTA)
      RETURN
      END
*
      SUBROUTINE NEIGH9 (NC,N,K,M,POIDS)
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   NC,N,K,M
      REAL      POIDS
*----
*  LOCAL VARIABLES
*----
      LOGICAL   EVEN
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBA,FSTA,LSTA
*
      ALLOCATE(NBA(NC+2),FSTA(NC+2),LSTA(NC+2))
      EVEN=.TRUE.
      CALL XDISET(NBA,NC+2,0)
      CALL XDISET(FSTA,NC+2,0)
      CALL XDISET(LSTA,NC+2,0)
      POIDS = 1.
      NBA(1)  = 1
      LSTA(1) = 1
      FSTA(1) = 1
      DO 36 L = 2,NC+1
         NBA(L)  = (L-1)*6
         LSTA(L) = 1+3*L*(L-1)
         FSTA(L) = 1+LSTA(L-1)
 36   CONTINUE
*
      I=0
      IF (N.EQ.1) THEN
         M =  K+1
         RETURN
      ELSE IF(N.GT.1) THEN
         DO 37 I0 = 2,NC
            IF ((N.GE.FSTA(I0)).AND.(N.LE.LSTA(I0))) THEN
               I=I0
               GO TO 38
            ENDIF
 37      CONTINUE
         IF (I.EQ.0) CALL XABORT('NEIGH9: ALGORITHM FAILURE.')
      ENDIF
*
 38   N1 = FSTA(I)
      N2 = FSTA(I) + (I-1)
      N3 = FSTA(I) + 2*(I-1)
      N4 = FSTA(I) + 3*(I-1)
      N5 = FSTA(I) + 4*(I-1)
      N6 = FSTA(I) + 5*(I-1)
      N7 = LSTA(I)
*
      IF (K.EQ.1) THEN
*
         IF (N.EQ.N1) THEN
            M =  FSTA(I+1)
         ELSE IF ((N.GT.N1).AND.(N.LT.N2)) THEN
            M =  N+NBA(I)
         ELSE IF (N.EQ.N2) THEN
            M =  FSTA(I+1)+I-1
         ELSE IF ((N.GT.N2).AND.(N.LE.N3)) THEN
            M =  N-1
         ELSE IF ((N.GT.N3).AND.(N.LT.N4)) THEN
            M =  N-NBA(I-1)-3
         ELSE IF (N.EQ.5) THEN
            M =  1
         ELSE IF (N.EQ.N4) THEN
            M =  FSTA(I-1)+3*(I-2)
         ELSE IF ((N.GT.N4).AND.(N.LT.N5)) THEN
            M =  N-NBA(I-1)-3
         ELSE IF ((N.GE.N5).AND.(N.LT.N6)) THEN
            M =  N+1
         ELSE IF (N.EQ.7) THEN
            M =  19
         ELSE IF ((N.GE.N6).AND.(N.LE.N7)) THEN
            M =  N+NBA(I)+6
         ENDIF
*
      ELSE IF (K.EQ.2) THEN
*
         IF (N.EQ.N1) THEN
            M =  FSTA(I+1)+1
         ELSE IF ((N.GT.N1).AND.(N.LT.N2)) THEN
            M =  N+NBA(I)+1
         ELSE IF (N.EQ.N2) THEN
            M =  FSTA(I+1)+I
         ELSE IF ((N.GT.N2).AND.(N.LE.N3)) THEN
            M =  N+NBA(I)+1
         ELSE IF ((N.GT.N3).AND.(N.LE.N4)) THEN
            M =  N-1
         ELSE IF ((N.GT.N4).AND.(N.LT.N5)) THEN
            M =  N-NBA(I-1)-4
         ELSE IF (N.EQ.6) THEN
            M =  1
         ELSE IF (N.EQ.N5) THEN
            M =  FSTA(I-1)+4*(I-2)
         ELSE IF ((N.GT.N5).AND.(N.LT.N6)) THEN
            M =  N-NBA(I-1)-4
         ELSE IF ((N.GE.N6).AND.(N.LT.N7)) THEN
            M =  N+1
         ELSE IF (N.EQ.N7) THEN
            M =  FSTA(I)
         ENDIF
*
      ELSE IF (K.EQ.3) THEN
*
         IF ((N.GE.N1).AND.(N.LT.N2)) THEN
            M =  N+1
         ELSE IF (N.EQ.N2) THEN
            M =  FSTA(I+1)+I+1
         ELSE IF ((N.GT.N2).AND.(N.LE.N4)) THEN
            M =  N+NBA(I)+2
         ELSE IF ((N.GT.N4).AND.(N.LE.N5)) THEN
            M =  N-1
         ELSE IF (N.EQ.7) THEN
            M =  1
         ELSE IF ((N.GT.N5).AND.(N.LT.N7)) THEN
            M =  N-NBA(I-1)-5
         ELSE IF (N.EQ.N7) THEN
            M =  FSTA(I-1)
         ENDIF
*
      ELSE IF (K.EQ.4) THEN
*
         IF (N.EQ.2) THEN
            M =  1
         ELSE IF ((N.GE.N1).AND.(N.LT.N2)) THEN
            M =  N-NBA(I-1)
         ELSE IF ((N.GE.N2).AND.(N.LT.N3)) THEN
            M =  N+1
         ELSE IF ((N.GE.N3).AND.(N.LE.N5)) THEN
            M =  N+NBA(I)+3
         ELSE IF ((N.GT.N5).AND.(N.LE.N6)) THEN
            M =  N-1
         ELSE IF ((N.GT.N6).AND.(N.LT.N7)) THEN
            M =  N-NBA(I-1)-6
         ELSE IF (N.EQ.N7) THEN
            M =  LSTA(I-1)
         ENDIF
*
      ELSE IF (K.EQ.5) THEN
*
         IF (N.EQ.N1) THEN
            M =  LSTA(I)
         ELSE IF (N.EQ.3) THEN
            M =  1
         ELSE IF (N.EQ.7) THEN
            M =  17
         ELSE IF ((N.GT.N1).AND.(N.LT.N3)) THEN
            M =  N-NBA(I-1)-1
         ELSE IF ((N.GE.N3).AND.(N.LT.N4)) THEN
            M =  N+1
         ELSE IF ((N.GE.N4).AND.(N.LE.N6)) THEN
            M =  N+NBA(I)+4
         ELSE IF ((N.GT.N6).AND.(N.LE.N7)) THEN
            M =  N-1
         ENDIF
*
      ELSE IF (K.EQ.6) THEN
*
         IF (N.EQ.N1) THEN
            M =  LSTA(I+1)
         ELSE IF ((N.GT.N1).AND.(N.LE.N2)) THEN
            M =  N-1
         ELSE IF ((N.GT.N2).AND.(N.LT.N4)) THEN
            M =  N-NBA(I-1)-2
         ELSE IF ((N.GE.N4).AND.(N.LT.N5)) THEN
            M =  N+1
         ELSE IF ((N.GE.N5).AND.(N.LE.N7)) THEN
            M =  N+NBA(I)+5
         ENDIF
      ENDIF
      DEALLOCATE(NBA,FSTA,LSTA)
      RETURN
      END
*
      SUBROUTINE NEIG10 (NC,N,K,M,POIDS)
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   NC,N,K,M
      REAL      POIDS
*----
*  LOCAL VARIABLES
*----
      LOGICAL   EVEN,OUTER
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBA,FSTA,LSTA
*
      ALLOCATE(NBA(NC+2),FSTA(NC+2),LSTA(NC+2))
      EVEN=.TRUE.
      CALL XDISET(NBA,NC+2,0)
      CALL XDISET(FSTA,NC+2,0)
      CALL XDISET(LSTA,NC+2,0)
      FSTA(1) = 1
      IL=0
      DO 39 L = 1,NC+1,2
         NBA(L)   = 1+IL
         NBA(L+1) = 1+IL
         IL = IL+1
 39   CONTINUE
      DO 40 L = 2,NC+1
         FSTA(L) = FSTA(L-1)+NBA(L-1)
 40   CONTINUE
      IL=0
      DO 41 L = 1,NC+1,2
         LSTA(L)   = FSTA(L)+IL
         LSTA(L+1) = FSTA(L+1)+IL
         IL = IL+1
 41   CONTINUE
*
      I=1
      IF (N.GT.1) THEN
         I=0
         DO 42 I0 = 1,NC
            IF ((N.GE.FSTA(I0)).AND.(N.LE.LSTA(I0))) THEN
               I=I0
               GO TO 43
            ENDIF
 42      CONTINUE
 43      IF (I.EQ.0) CALL XABORT('NEIG10: ALGORITHM FAILURE.')
      ENDIF
*
      N1 = FSTA(I)
      N2 = LSTA(I)
      EVEN = MOD(I,2).EQ.0
      OUTER = I.EQ.NC
*
      IF (K.EQ.1) THEN
*
         IF (N.EQ.1) THEN
            M =  2
         ELSE IF (OUTER.AND.(N.EQ.2)) THEN
            M = -2
         ELSE IF (OUTER.AND.(N.EQ.N2)) THEN
            M = -(N-1)
         ELSE IF (OUTER.AND.EVEN) THEN
            M = -(N-NBA(I-1)+1)
         ELSE IF (OUTER.AND.(.NOT.EVEN)) THEN
            M = -(N-NBA(I-1))
         ELSE IF (EVEN) THEN
            M =  N+NBA(I)+1
         ELSE IF (.NOT.EVEN) THEN
            M =  N+NBA(I)
         ENDIF
*
      ELSE IF (K.EQ.2) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF (OUTER.AND.(N.EQ.N2)) THEN
            M = -LSTA(I-1)
         ELSE IF (N.EQ.N2) THEN
            M = -(LSTA(I+1)-1)
         ELSE
            M =  N+1
         ENDIF
*
      ELSE IF (K.EQ.3) THEN
*
         IF ((N.EQ.1).OR.(N.EQ.2)) THEN
            M = -2
         ELSE IF (N.EQ.N2) THEN
            M = -(N-1)
         ELSE IF (EVEN) THEN
            M =  N-NBA(I-1)+1
         ELSE IF (.NOT.EVEN) THEN
            M =  N-NBA(I-1)
         ENDIF
*
      ELSE IF (K.EQ.4) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -FSTA(I-1)
         ELSE IF (EVEN) THEN
            M =  N-NBA(I-1)
         ELSE IF (.NOT.EVEN) THEN
            M =  N-NBA(I-1)-1
         ENDIF
*
      ELSE IF (K.EQ.5) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF ((N.EQ.N1).AND.EVEN) THEN
            M =  N
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -(N+1)
         ELSE
            M =  N-1
         ENDIF
*
      ELSE IF (K.EQ.6) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF (OUTER.AND.(N.EQ.N1)) THEN
            M = -FSTA(I-1)
         ELSE IF (OUTER.AND.EVEN) THEN
            M = -(N-NBA(I-1))
         ELSE IF (OUTER.AND.(.NOT.EVEN)) THEN
            M = -(N-NBA(I-1)-1)
         ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
            M = -FSTA(I+1)
         ELSE IF (EVEN) THEN
            M =  N+NBA(I)
         ELSE IF (.NOT.EVEN) THEN
            M =  N+NBA(I)-1
         ENDIF
*
      ENDIF
*
      IF (N.EQ.1) THEN
         POIDS = 1./12.
      ELSE IF (OUTER.AND.(N.EQ.N2)) THEN
         POIDS = 1./6.
      ELSE IF (OUTER.AND.(N.EQ.N1).AND.(.NOT.EVEN)) THEN
         POIDS = 0.25
      ELSE IF (OUTER.OR.(N.EQ.N2)) THEN
         POIDS = 0.5
      ELSE IF ((N.EQ.N1).AND.(.NOT.EVEN)) THEN
         POIDS = 0.5
      ELSE
         POIDS = 1.
      ENDIF
      DEALLOCATE(NBA,FSTA,LSTA)
      RETURN
      END
*
      SUBROUTINE NEIG11 (NC,N,K,M,POIDS)
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   NC,N,K,M
      REAL      POIDS
*----
*  LOCAL VARIABLES
*----
      LOGICAL   EVEN,OUTER
      INTEGER, DIMENSION(:), ALLOCATABLE :: NBA,FSTA,LSTA
*
      ALLOCATE(NBA(NC+2),FSTA(NC+2),LSTA(NC+2))
      EVEN=.TRUE.
      CALL XDISET(NBA,NC+2,0)
      CALL XDISET(FSTA,NC+2,0)
      CALL XDISET(LSTA,NC+2,0)
      NBA(1)  = 1
      LSTA(1) = 1
      FSTA(1) = 1
      FSTA(2) = 2
      DO 45 L = 2,NC+1
         NBA(L)  = L
         LSTA(L) = L+LSTA(L-1)
         FSTA(L+1) = L+FSTA(L)
 45   CONTINUE
*
      I=0
      DO 46 I0 = 1,NC
         IF ((N.GE.FSTA(I0)).AND.(N.LE.LSTA(I0))) THEN
            I=I0
            GO TO 47
         ENDIF
 46   CONTINUE
      IF (I.EQ.0) CALL XABORT('NEIG11: ALGORITHM FAILURE.')
*
 47   N1 = FSTA(I)
      N2 = LSTA(I)
      OUTER = I.EQ.NC
*
      IF (K.EQ.1) THEN
*
         IF (N.EQ.1) THEN
            M =  3
         ELSE IF (OUTER.AND.(N.EQ.N2)) THEN
            M = -(N-1)
         ELSE IF (OUTER) THEN
            M = -(N-NBA(I-1))
         ELSE
            M =  N+NBA(I)+1
         ENDIF
*
      ELSE IF (K.EQ.2) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF (OUTER.AND.(N.EQ.N2)) THEN
            M = -(N-NBA(I-1)-1)
         ELSE IF (N.EQ.N2) THEN
            M = -(N+NBA(I))
         ELSE
            M =  N+1
         ENDIF
*
      ELSE IF (K.EQ.3) THEN
*
         IF (N.EQ.1) THEN
            M = -3
         ELSE IF (N.EQ.N2) THEN
            M = -(N-1)
         ELSE
            M =  N-NBA(I-1)
         ENDIF
*
      ELSE IF (K.EQ.4) THEN
*
         IF (N.EQ.1) THEN
            M = -2
         ELSE IF (N.EQ.N1) THEN
            M = -(N+1)
         ELSE
            M =  N-NBA(I-1)-1
         ENDIF
*
      ELSE IF (K.EQ.5) THEN
*
         IF (N.EQ.1) THEN
            M = -3
         ELSE IF (OUTER.AND.(N.EQ.N1)) THEN
            M = -(N-NBA(I-1))
         ELSE IF (N.EQ.N1) THEN
            M = -(N+NBA(I)+1)
         ELSE
            M =  N-1
         ENDIF
*
      ELSE IF (K.EQ.6) THEN
*
         IF (N.EQ.1) THEN
            M =  2
         ELSE IF (OUTER.AND.(N.EQ.N1)) THEN
            M = -(N+1)
         ELSE IF (OUTER) THEN
            M = -(N-NBA(I-1)-1)
         ELSE
            M =  N+NBA(I)
         ENDIF
*
      ENDIF
*
      IF (N.EQ.1) THEN
         POIDS = 1./6.
      ELSE IF (OUTER.AND.((N.EQ.N1).OR.(N.EQ.N2))) THEN
         POIDS = 1./6.
      ELSE IF (OUTER.OR.(N.EQ.N1).OR.(N.EQ.N2)) THEN
         POIDS = 0.5
      ELSE
         POIDS = 1.
      ENDIF
      DEALLOCATE(NBA,FSTA,LSTA)
      RETURN
      END
