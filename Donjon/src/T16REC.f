*DECK T16REC
      SUBROUTINE T16REC(IFT16 ,IPRINT,INEXTR)
*
*----
*
*Purpose:
*  Locate next set of records on tape16.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IFT16   tape16 file unit.
* IPRINT  print level where:
*         =0 for no print; >=  1 print processing option.
* INEXTR  next record set to read.
*
*----
*
      IMPLICIT         NONE
      INTEGER          IFT16,IPRINT,INEXTR
*----
*  T16 PARAMETERS
*----
      INTEGER          MAXKEY
      PARAMETER       (MAXKEY=3)
      CHARACTER        TKEY1(MAXKEY)*10,TKEY2(MAXKEY)*10,
     >                 RKEY1*10,RKEY2*10
      INTEGER          NKEY,IOPT,NBE,NSKIPR,ISKIPR
*----
*  LOCAL VARIABLES
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='T16REC')
*----
*  REWIND AND SKIP FIRST INEXTR-1 SETS OF RECORDS
*----
      REWIND(IFT16)
      NSKIPR=INEXTR
      TKEY1(1)='MTR       '
      TKEY2(1)='FEWGROUPS '
      NKEY=1
      IOPT=-1
      DO ISKIPR=1,NSKIPR
        CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >              NBE   )
        IF(NBE .EQ. -1) THEN
          WRITE(IOUT,9000) NAMSBR,TKEY1(1),TKEY2(1),INEXTR
          CALL XABORT(NAMSBR//': INVALID RECORD NUMBER ON TAPE16')
        ENDIF
        READ(IFT16) RKEY1,RKEY2,NBE
      ENDDO
      RETURN
*----
*  ABORT FORMAT
*----
 9000 FORMAT(1X,A6,1X,7('*'),' ERROR ',7('*')/
     >       8X,I6,' TAPE16 RECORD WITH KEYS =',2(A10,2X),
     >       'NOT FOUND'/
     >       8X,21('*'))
      END
