*DECK TINFL
      SUBROUTINE TINFL (NNS,NW,NW2,NK)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Produce the useful vector for refuelling, according
* to a given refuelling-scheme
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
*
*Author(s): 
* E. Varin, M. Guyot
*
*Parameters: input/output
* NNS    Number corresponding to the refuelling type
* NW     Vector corresponding to the refuelling type
*           NW > 0  : Position of the bundle before refuelling     
*           NW = 0  : Insertion of a new bundle
* NW2    Vector NW when the refueling is negative
* NK     Number of bundles per channel
*
*-----------------------------------------------------------------------
*
      IMPLICIT    NONE
*
      INTEGER     NNS,NK,NW(NK),NW2(NK)
      INTEGER     MODEID,I,IP,NL
*
      MODEID = NNS
*
*---- MODE DE RECHARGEMENT GENERALISE
*
      IF(MODEID.GT.NK) THEN
        WRITE(6,'(13H @TINFL: NNS=,I6,4H NK=,I6)') NNS,NK
        CALL XABORT('@TINFL: ONLY BI-DIRECTIONNAL REFUELING ')
      ELSE
*
*------- MODE DE RECHARGEMENT DIRECT
*
         CALL XDISET(NW2,NK,0)
         DO 10 I=1,MODEID
            NW(I) = 0
 10      CONTINUE
*
         IF(MODEID.NE.NK) THEN
            NL = NK - MODEID
            DO 20 I=1,NL
               IP = MODEID + I
               NW(IP) = I
               NW2(I)=I+NNS
  20        CONTINUE
          ENDIF
*
      ENDIF
      RETURN
      END
