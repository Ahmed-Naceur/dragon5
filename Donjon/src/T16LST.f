*DECK T16LST
      SUBROUTINE T16LST(IFT16 )
*
*----
*
*Purpose:
*  Print records stored on tape16.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IFT16   tape16 file unit.
*
*----
*
      IMPLICIT         NONE
      INTEGER          IFT16
*----
*  T16 KEYS
*----
      CHARACTER        TKEY1*10,TKEY2*10
      INTEGER          NKEY,IOPT,NBE
*----
*  LOCAL VARIABLES
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='T16LST')
      INTEGER          IPRINT
*----
*  LIST TAPE16 RECORDS AFTER REWINDING
*----
      WRITE(IOUT,6000) NAMSBR
      IPRINT=10000
      IOPT=1
      NKEY=1
      TKEY1='          '
      TKEY2=TKEY1
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE)
      WRITE(IOUT,6001)
      RETURN
      
*----
*  PRINT FORMAT
*----
 6000 FORMAT( 1X, 'PRINTING CONTENTS OF TAPE16 FILE USING ',A6)
 6001 FORMAT( 1X, 'END OF TAPE16 FILE REACHED')
      END
