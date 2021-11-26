*DECK DETLIN
      SUBROUTINE DETLIN(X1,X2,Y1,Y2,BS,CS)
*
*----------------------------------------------------------------------
*Purpose:
* Routine calculating the linear coefficient needed for a linear
* interpolation Y = BS*X + CS
*
*Author(s): 
* M. Beaudet
*
*Parameters: 
* X1      
* X2      
* Y1      
* Y2      
* BS      
* CS      
*
*----------------------------------------------------------------------
*
      BS   = (Y1-Y2)/(X1-X2)
      CS   = Y1-BS*X1
      RETURN
      END
