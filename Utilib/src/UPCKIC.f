*DECK UPCKIC
      SUBROUTINE UPCKIC(IV,VC,N)
C----
C  UNPACK CHARACTERS FROM AN INTEGER ARRAY 4 AT A TIME
C----
      IMPLICIT   NONE
      INTEGER    NANSI
      PARAMETER (NANSI=256)
      INTEGER    N,I,J,K,II,KK,NUM,ND,NR
      INTEGER    IV(*)
      CHARACTER  VC(N)*8
C----
C  UNPACK ANSI CHARACTER FROM INTEGER VECTOR
C----
      II = 1
      DO 10 I = 1, N
        KK = 1
        DO 20 J = 1, 2
          NUM = IV(II)
          II = II + 1
          DO 30 K = 1, 4
            ND = NUM/NANSI
            NR = NUM -NANSI*ND
            NUM = ND
            VC(I)(KK:KK)=CHAR(NR)
            KK = KK + 1
 30       CONTINUE
 20     CONTINUE
 10   CONTINUE
      RETURN
      END
