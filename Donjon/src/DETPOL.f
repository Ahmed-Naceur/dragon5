*DECK DETPOL
      REAL FUNCTION DETPOL(VECT,IXX,JJJ,K0,I1,I2,I3,X1,X2,X3,X)
*
*----------------------------------------------------------------------
*Purpose:
* Function performing the parabolic interpolation at X.
*
*Author(s): 
* M. Beaudet
*
*Parameters: 
* DETPOL   
* VECT     
* IXX      
* JJJ      
* K0       
* I1       
* I2       
* I3       
* X1       
* X2       
* X3       
* X        
*
*----------------------------------------------------------------------
*
      INTEGER*4    IXX(*)
      REAL*4       VECT(*)
*
      CHARACTER*6  CLNAME
      DATA CLNAME /'INTPOL'/
*
      IJK1 = IXX(JJJ+K0+I1)
      IJK2 = IXX(JJJ+K0+I2)
      IJK3 = IXX(JJJ+K0+I3)
*
      IZERO = 0
*
      IF (IJK1.LE.0) IZERO = IZERO + 1
      IF (IJK2.LE.0) IZERO = IZERO + 1
      IF (IJK3.LE.0) IZERO = IZERO + 1
*
      IF (IZERO.GE.2) CALL XABORT('DETPOL: INVALID VALUE OF INDICES')
*
      IF (IJK1.LE.0) THEN
         A2   = VECT(IJK2)
         A3   = VECT(IJK3)
         CALL DETLIN(X2,X3,A2,A3,BE,CE)
         AH = 0.0
*
      ELSE IF (IJK2.LE.0) THEN
         A1   = VECT(IJK1)
         A3   = VECT(IJK3)
         CALL DETLIN(X1,X3,A1,A3,BE,CE)
         AH = 0.0
*
      ELSE IF (IJK3.LE.0) THEN
         A1   = VECT(IJK1)
         A2   = VECT(IJK2)
         CALL DETLIN(X1,X2,A1,A2,BE,CE)
         AH = 0.0
*
      ELSE
         A1   = VECT(IJK1)
         A2   = VECT(IJK2)
         A3   = VECT(IJK3)
         CALL DETPAR(X1,X2,X3,A1,A2,A3,AH,BE,CE)
      ENDIF
*
      DETPOL = AH*X*X + BE*X + CE
*
      RETURN
      END
