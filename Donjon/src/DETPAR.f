*DECK DETPAR
      SUBROUTINE DETPAR(X1,X2,X3,Y1,Y2,Y3,AS,BS,CS)
*
*----------------------------------------------------------------------
*Purpose:
* Routine calculating the parabolic coefficients needed for
* a parabolic interpolation Y = AS*X*X + BS*X + CS 
*
*Author(s): 
* M. Beaudet
*
*Parameters: 
* X1      
* X2      
* X3      
* Y1      
* Y2      
* Y3      
* AS      
* BS      
* CS      
*
*----------------------------------------------------------------------
*
      ANUM = Y1*(X2-X3)+Y3*(X1-X2)+Y2*(X3-X1)
      ADEN = (X1-X2)*(X1-X3)*(X2-X3)
      AS   = ANUM/ADEN
      BS   = (Y2-Y3-AS*(X2*X2-X3*X3))/(X2-X3)
      CS   = Y1-BS*X1-AS*X1*X1
      RETURN
      END
