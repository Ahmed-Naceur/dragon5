*DECK T16FND
      SUBROUTINE T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >                  NBELEM)
*
*----
*
*Purpose:
*  Find next record on tape16 identified by keys TKEY1 and TKEY2.
*
*Author(s): 
* G. Marleau
*
*Parameters: input
* IFT16   tape16 file unit.
* IPRINT  print level where:
*         <100 for no print; 
*         >=100 print record to read;
*         >=10000 print all record read.
* IOPT    processing option with:
*         =-1 start at current position and read to end of file with
*             no backspace before return;
*         =0 start at current position and read to end of file with
*             backspace before return;
*         =1 rewind before reading and read to end of file;
*         =2 start at current position, rewind, start at beginning of
*         file until end of file.
* NKEY    number of keys set to test:
*         =1 search for TKEY1(1),TKEY2(1) until end of file;
*         >1 search for TKEY1(1),TKEY2(1) until
*         (TKEY1(IK),TKEY2(IK),IK=2,NKEY) or end of file.
* TKEY1   primary key.
* TKEY2   secondary key.
*
*Parameters: output
* NBELEM  number of element found on record with:
*         <-1 record not found before alternative keys -NBELEM ;
*         =-1 record not found before end of files;
*         >=0 record found with NBELEM elements.
*
*----
*
      IMPLICIT         NONE
      INTEGER          IFT16,IPRINT,IOPT,NKEY,NBELEM
      CHARACTER        TKEY1(NKEY)*10,TKEY2(NKEY)*10
*----
*  LOCAL VARIABLES
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='T16FND')
      CHARACTER        RKEY1*10,RKEY2*10
      INTEGER          NBE,IEND,IKEY
*----
*  Print keys if required
*----
      IF(IPRINT .GE. 100) THEN
        IF(IPRINT .LT. 10000) THEN
          WRITE(6,6000) TKEY1(1),TKEY2(1)
        ENDIF
      ENDIF
*----
*  REWIND FILE FIRST IF IOPT=1
*----
      IEND=1
      IF(IOPT .EQ. 1) THEN
        REWIND(IFT16)
      ELSE IF (IOPT .EQ. 2) THEN
        IEND=0
      ENDIF
*----
*  LOOP FOR READ
*----
 100  CONTINUE
      READ(IFT16,END=105) RKEY1,RKEY2,NBE
      IF(IPRINT .GE. 10000) THEN
        WRITE(6,6003) RKEY1,RKEY2,NBE
      ENDIF
      IF(RKEY1 .EQ. TKEY1(1) .AND.
     >   RKEY2 .EQ. TKEY2(1)       ) THEN
*----
*  KEYS FOUND BACKSPACE AND RETURN
*----
        NBELEM=NBE
        IF(IOPT .GE. 0) BACKSPACE(IFT16)
        IF(IPRINT .GE. 100) THEN
          WRITE(6,6001) RKEY1,RKEY2,NBELEM
        ENDIF
        RETURN
      ELSE IF(NKEY .GE. 2) THEN
        DO IKEY=2,NKEY
          IF(RKEY1 .EQ. TKEY1(IKEY) .AND.
     >       RKEY2 .EQ. TKEY2(IKEY)       ) THEN
            NBELEM=-IKEY
            IF(IOPT .GE. 0) BACKSPACE(IFT16)
            IF(IPRINT .GE. 100) THEN
              WRITE(6,6004) RKEY1,RKEY2,NBE,
     >                      TKEY1(1),TKEY2(1)
            ENDIF
            RETURN
          ENDIF
        ENDDO
      ENDIF
*----
*  KEYS NOT FOUND READ NEXT RECORD
*----
      GO TO 100
*----
*  END OF FILE REACHED
*----
 105  CONTINUE
      IF(IEND .EQ. 0) THEN
*----
*  REWIND FILE AND CONTINUE READ
*----
        IEND=1
        REWIND(IFT16)
        GO TO 100
      ENDIF
*----
*  RECORD ABSENT, RETURN
*----
      NBELEM=-1
      IF(IPRINT .GE. 100) THEN
        IF(IPRINT .LT. 10000) THEN
          WRITE(6,6002) TKEY1(1),TKEY2(1)
        ENDIF
      ENDIF
      RETURN
*----
*  PRINT FORMAT
*----
 6000 FORMAT( 1X, 'FIND T16 RECORD = ',2(A10,2X))
 6001 FORMAT( 1X, '     T16 RECORD = ',2(A10,2X),I10,
     >        1X,'FOUND')
 6002 FORMAT( 1X, '     T16 RECORD = ',2(A10,2X),10X,
     >        1X,'NOT FOUND')
 6003 FORMAT(11X,'T16 RECORD READ = ',2(A10,2X),I10)
 6004 FORMAT( 1X,'T16 STOP RECORD = ',2(A10,2X),I10,
     >        1X,'FOUND BEFORE RECORD = ',2(A10,2X))
      END
