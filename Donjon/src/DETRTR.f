*DECK DETRTR
      SUBROUTINE DETRTR(DA,A,IA,A1,A2,A3,II1,II2,II3)
*
*----------------------------------------------------------------------
*Purpose:
* Obtain the coordinates of a point where the interpolation is
* performed
*
*Author(s):
* ???
*
*Parameters: 
* DA      
* A       
* IA      
* A1      
* A2      
* A3      
* II1     
* II2     
* II3     
*
*----------------------------------------------------------------------
*
      DIMENSION A(*)
      CHARACTER*6 CLNAME
*
      CLNAME = 'SORTR '
      DIF1 = 1000000.
      DIF2 = 1000001.
      DIF3 = 1000002.
      II1  = 1000000
      II2  = 1000001
      II3  = 1000002
*
      DO 10 II=1,IA
         DIF = ABS(DA-A(II))
         IF ( DIF .LE. DIF1 ) THEN
            DIF3 = DIF2
            DIF2 = DIF1
            DIF1 = DIF
            II3  = II2
            II2  = II1
            II1  = II
         ELSE IF ( DIF .LE. DIF2 ) THEN
            DIF3 = DIF2
            DIF2 = DIF
            II3  = II2
            II2  = II
         ELSE IF ( DIF .LE. DIF3 ) THEN
            DIF3 = DIF
            II3  = II
         ENDIF
  10  CONTINUE
      A1 = A(II1)
      A2 = A(II2)
      A3 = A(II3)
      RETURN
      END
