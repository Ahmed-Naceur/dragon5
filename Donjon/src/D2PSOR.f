*DECK D2PSOR
      SUBROUTINE D2PSOR(TAB,DIMTAB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Sort D2PSOR state variable real array to match GENPMAXS order, in
* ascendent order
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* TAB     vector of rank index of state variables
* DIMTAB  dimension of TAB
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER DIMTAB
      REAL TAB(DIMTAB)
*----
*  LOCAL VARIABLES
*----
      REAL Rtmp
      INTEGER :: I, J

      DO I = 2, DIMTAB
         Rtmp = TAB(I)
         DO J = I-1, 1, -1
          IF (Rtmp < TAB(J)) THEN
           TAB(J+1) = TAB(J)
          ELSE
           EXIT
          ENDIF
         ENDDO
         TAB(J+1) = Rtmp
      ENDDO
      RETURN
      END
