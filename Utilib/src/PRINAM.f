*DECK PRINAM
      SUBROUTINE PRINAM(T,A,N)
*
*-----------------------------------------------------------------------
*
* PRINT A REAL, INTEGER OR DOUBLE PRECISION ARRAY.
*
* INPUT PARAMETER:
*  T  : CHARACTER*6 NAME OF THE ARRAY.
*  A  : REAL ARRAY TO PRINT. DIMENSION A(N)
*  IA : INTEGER ARRAY TO PRINT. DIMENSION IA(N)
*  DA : DOUBLE PRECISION ARRAY TO PRINT. DIMENSION DA(N)
*
*-----------------------------------------------------------------------
*
      CHARACTER*6 T
      INTEGER N
      REAL A(N)
      PARAMETER (MAX5=5,MAX10=10,IOUT=6)
*
      MN=MIN0(N,MAX5)
      WRITE (IOUT,10) T,(A(I),I=1,MN)
   10 FORMAT(1X,20('-')/1X,A6,1X,1P,5E13.6)
      IF(N.LT.MAX5) GOTO 30
      MN=MN+1
      WRITE (IOUT,20) (A(I),I=MN,N)
   20 FORMAT(8X,1P,5E13.6)
   30 CONTINUE
      RETURN
      END SUBROUTINE PRINAM
C ------------------ P R I N I M
      SUBROUTINE PRINIM(T,IA,N)
      CHARACTER*6 T
      INTEGER N,IA(N)
      PARAMETER (MAX5=5,MAX10=10,IOUT=6)
*
      MN=MIN0(N,MAX10)
      WRITE (IOUT,40) T,(IA(I),I=1,MN)
   40 FORMAT(1X,20('-')/1X,A6,1X,10I6)
      IF(N.LT.MAX10) GOTO 60
      MN=MN+1
      WRITE (IOUT,50) (IA(I),I=MN,N)
   50 FORMAT(8X,10I6)
   60 CONTINUE
      RETURN
      END SUBROUTINE PRINIM
C ------------------ P R I N D M
      SUBROUTINE PRINDM(T,DA,N)
      CHARACTER*6 T
      INTEGER N
      DOUBLE PRECISION DA(N)
      PARAMETER (MAX5=5,MAX10=10,IOUT=6)
*
      MN=MIN0(N,MAX5)
      WRITE (IOUT,70) T,(DA(I),I=1,MN)
   70 FORMAT(1X,20('-')/1X,A6,1X,1P,5D13.6)
      IF(N.LT.MAX5) GOTO 90
      MN=MN+1
      WRITE (IOUT,80) (DA(I),I=MN,N)
   80 FORMAT(8X,1P,5D13.6)
   90 CONTINUE
      RETURN
      END SUBROUTINE PRINDM
