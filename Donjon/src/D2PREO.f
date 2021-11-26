*DECK D2PREO
      SUBROUTINE D2PREO(IPDAT,VALPAR,IND,NPAR,NVAL,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* take into account the meanning of the control rod composition in
* Saphyb, attribute to each value of control rod the corresponding value
* in GENMAPXS formalism
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPDAT   address of the INFO data block
* VALPAR  vector of values for each state variable
* IND     index of the control rod parameter
* NPAR    number of state variables
* NVAL    number of values for control rod parameter
* IPRINT  control the printing on screen
*
*Parameters: 
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT
      INTEGER NVAL
      REAL VALPAR(NPAR,100)
*----
*  LOCAL VARIABLES
*----
      ! USER INPUT: MEANING OF BARR PARAMETERS : LOCATED AT INFO/CRDINF
      INTEGER CRDINF(NVAL)

      ! RECOVER CRDINF DATA BLOCK
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMGET(IPDAT,'BARR_INFO',CRDINF)

      DO I=1, NVAL
        VALPAR(IND,I)=CRDINF(I)
        IF (CRDINF(I)<0) THEN
          CALL XABORT('@D2PREO: CONTROL ROD COMPO MUST BE POSITIVE')
        ENDIF
      ENDDO
      ! ATTRIBUTION OF CRDINF TO THE BARR PARAMETERS

      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)

      ! EDIT THE LISTING FILE
      IF(IPRINT > 1)  THEN
        WRITE(6,*)
        WRITE(6,*) "******  CONTROL ROD COMPOSITION (IN SAPHYB) ******"
        WRITE(6,*)
        WRITE(6,*) "UNRODDED CROSS SECTIONS :",VALPAR(IND,1)
        DO J=2, NVAL
         WRITE(6,*) "RODDED COMPOSITION ",J-1,": ",VALPAR(IND,J)
        ENDDO
        WRITE(6,*)
      ENDIF
      END
