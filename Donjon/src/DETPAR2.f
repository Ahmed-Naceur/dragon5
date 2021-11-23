*DECK DETPAR2
      SUBROUTINE DETPAR2(V1,V2,V3,U1,U2,U3,AS,BS,CS)
*
*----------------------------------------------------------------------
*Purpose: routine de HQSIMEX
*
*Author(s): 
* M. Beaudet
*
*Parameters: 
* V1     
* V2     
* V3     
* U1     
* U2     
* U3     
* AS     
* BS     
* CS     
*
*----------------------------------------------------------------------
*
      CHARACTER*6 CLNAME
      CLNAME = 'PAR'
      ANUM = U1*(V2-V3)+U3*(V1-V2)+U2*(V3-V1)
      ADEN = (V1-V2)*(V1-V3)*(V2-V3)
      AS   = ANUM/ADEN
      BS   = (U2-U3-AS*(V2*V2-V3*V3))/(V2-V3)
      CS   = U1-BS*V1-AS*V1*V1
      RETURN
      END
