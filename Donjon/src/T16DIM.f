*DECK T16DIM
      SUBROUTINE T16DIM(IFT16 ,IPRINT,MXGRP ,SUBTIT,NEL   ,NG    ,
     >                  NGMTR ,NMATZ ,MTRMSH,NZONE ,NGREAC,NRCELA,
     >                  NRREGI,IFGMTR,IFGEDI)
*
*----
*
*Purpose:
*  Scan WIMS-AECL tape16 file for general dimensioning information.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IFT16   tape16 file unit.
* IPRINT  print level where:
*         =0 for no print; >=  1 print record to read;
*         >= 10 print all record read.
* MXGRP   maximum number or groups.
*
*Parameters: output
* SUBTIT  subtitle.
* NEL     number of isotopes on X-S library.
* NG      number of groups on X-S library.
* NGMTR   number of main transport groups.
* NMATZ   number of mixtures.
* MTRMSH  number of main transport mesh points.
* NZONE   number of zones.
* NGREAC  number of edit groups.
* NRCELA  number of CELLAV sets of records.
* NRREGI  number of REGION sets of records.
* IFGMTR  fewgroups for main transport.
* IFGEDI  fewgroups for edit.
*
*----
*
      IMPLICIT         NONE
      INTEGER          IFT16,IPRINT,MXGRP,NEL,NG,
     >                 NGMTR,NMATZ,MTRMSH,NZONE,
     >                 NGREAC,NRCELA,NRREGI
      INTEGER          IFGMTR(MXGRP),IFGEDI(MXGRP)
      CHARACTER        SUBTIT*240
*----
*  T16 KEYS
*----
      CHARACTER        CWVER*80,CLIBN*16,CASETL*128,
     >                 TKEY1*10,TKEY2*10,RKEY1*10,RKEY2*10,
     >                 WLEAK*10, WDIFF*10,WEDIT*10,BLANK*2
      INTEGER          NKEY,IOPT,NBE,NID,NJD,IR,JR
      REAL             RID
*----
*  LOCAL VARIABLES
*----
      INTEGER          IOUT,NFPR,NREGON,NM
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='T16DIM')
*----
*  READ GENERAL TAPE16 INFORMATION
*----
      IOPT=0
      NKEY=1
      SUBTIT=' '
      REWIND(IFT16)
*----
*  1) WIMS-AECL VERSION
*----
      TKEY1='PROCESSING'
      TKEY2='PROCESSING'
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1  ,TKEY2  ,
     >            NBE   )
      IF(NBE .GT. 0) THEN
        READ(IFT16) RKEY1,RKEY2,NBE,CWVER
      ELSE
        CALL XABORT(NAMSBR//': KEYS '//TKEY1//','//
     >              TKEY2//' NOT FOUND ON TAPE16')
      ENDIF
      SUBTIT(1:80)=CWVER
*----
*  2) LIBRARY NAME
*----
      TKEY1='PROCESSING'
      TKEY2='NDASTITLE '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1  ,TKEY2  ,
     >            NBE   )
      IF(NBE .GT. 0) THEN
        READ(IFT16) RKEY1,RKEY2,NBE, CLIBN
      ELSE
        CALL XABORT(NAMSBR//': KEYS '//TKEY1//','//
     >              TKEY2//' NOT FOUND ON TAPE16')
      ENDIF
      SUBTIT(81:104)=' ------ '//CLIBN
*----
*  3) CASE TITLE
*----
      TKEY1='TITLE     '
      TKEY2='CARD      '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1  ,TKEY2  ,
     >            NBE   )
      IF(NBE .GT. 0) THEN
        READ(IFT16) RKEY1,RKEY2,NBE,CASETL
      ELSE
        CALL XABORT(NAMSBR//': KEYS '//TKEY1//','//
     >              TKEY2//' NOT FOUND ON TAPE16')
      ENDIF
      SUBTIT(105:240)=' ------ '//CASETL
*----
*  4) WIMS CONSTANTS
*----
      TKEY1='WIMS      '
      TKEY2='CONSTANTS '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1  ,TKEY2  ,
     >            NBE   )
      IF(NBE .GT. 0) THEN
        READ(IFT16) RKEY1,RKEY2,NBE,NEL,NG,(NID,IR=1,8),NGMTR,
     >              (NID,IR=1,6),NMATZ,NM
      ELSE
        CALL XABORT(NAMSBR//': KEYS '//TKEY1//','//
     >              TKEY2//' NOT FOUND ON TAPE16')
      ENDIF
*----
*  5) MAIN TRANSPORT GROUPS
*----
      TKEY1='MTR       '
      TKEY2='FEWGROUPS '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1  ,TKEY2  ,
     >            NBE   )
      IF(NBE .GT. 0) THEN
        READ(IFT16) RKEY1,RKEY2,NBE,(IFGMTR(IR),IR=1,NGMTR)
      ELSE
        CALL XABORT(NAMSBR//': KEYS '//TKEY1//','//
     >              TKEY2//' NOT FOUND ON TAPE16')
      ENDIF
*----
*  6) DIMENSION OF TRANSPORT MESH
*     PRESENT ONLY IF MTRFLX KEY ACTIVATED
*----
      TKEY1='MTRFLX    '
      TKEY2='FLUX      '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1  ,TKEY2  ,
     >            NBE   )
      IF(NBE .GT. 0) THEN
        READ(IFT16) RKEY1,RKEY2,NBE,NID,MTRMSH
      ELSE
        REWIND(IFT16)
        MTRMSH=0
        IF(IPRINT .GE. 10)
     >  WRITE(IOUT,8000) NAMSBR,TKEY1,TKEY2,'MTRMSH',MTRMSH
      ENDIF
*----
*  7) NUMBER OF FUEL PIN RINGS
*     PRESENT ONLY FOR BURNUP CASES WITH CLUSTER GEOMETRY
*----
*----- A.ZH. THIS RECORD CAN HAVE A DIFFERENT INTERPRETATION-----
      TKEY1='CELLAV    '
      TKEY2='PINBURNUP '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1  ,TKEY2  ,
     >            NBE   )
      IF(NBE .GT. 0) THEN
        READ(IFT16) RKEY1,RKEY2,NBE
        NFPR=(NBE-1)/3
      ELSE
        REWIND(IFT16)
        NFPR=0
        IF(IPRINT .GE. 10)
     >  WRITE(IOUT,8000) NAMSBR,TKEY1,TKEY2,'NFPR  ',NFPR
      ENDIF   
*----
*  8) NUMBER OF ZONES
*----     
      REWIND(IFT16) 
      TKEY1='REGION    '
      TKEY2='DESCRIPTON'
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1  ,TKEY2  ,
     >            NBE   )      
      IF(NBE .GT. 0) THEN       
        READ(IFT16) RKEY1,RKEY2,NBE,NZONE
      ELSE       
        CALL XABORT(NAMSBR//': KEYS '//TKEY1//','//
     >              TKEY2//' NOT FOUND ON TAPE16')
      ENDIF      
*----
*  9) NUMBER OF EDIT REGIONS
*----
      TKEY1='REGION    '
      TKEY2='DIMENSIONS'
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1  ,TKEY2  ,
     >            NBE   )
      IF(NBE .GT. 0) THEN
        READ(IFT16) RKEY1,RKEY2,NBE,NREGON
      ELSE
        CALL XABORT(NAMSBR//': KEYS '//TKEY1//','//
     >              TKEY2//' NOT FOUND ON TAPE16')
      ENDIF
*----
* 10) NUMBER OF EDIT GROUPS
*     PRESENT ONLY IF REACTION KEY ACTIVATED
*----
      TKEY1='REACTION  '
      TKEY2='FLUX      '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1  ,TKEY2  ,
     >            NBE   )
      IF(NBE .GT. 0) THEN
         READ(IFT16) RKEY1,RKEY2,NBE,WLEAK,WDIFF,WEDIT,BLANK,
     >               (NID,IR=1,2),NGREAC,
     >               ((RID,IR=1,NZONE),JR=1,NG),
     >               (IFGEDI(IR),IR=1,NGREAC)
      ELSE
        NGREAC=0
        IF(IPRINT .GE. 10)
     >  WRITE(IOUT,8000) NAMSBR,TKEY1,TKEY2,'NGREAC',NGREAC
      ENDIF
*----
*  FIND THE NUMBER OF SETS OF CELLAV RECORDS
*  BASED ON THE PRESENCE OF CELLAV,NGROUP KEYS
*  ALSO TEST FOR NGMTR CONSISTENCY
*----  
      REWIND(IFT16)
      NRCELA=0
      TKEY1='CELLAV    '
      TKEY2='NGROUPS   '
 100  CONTINUE
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1  ,TKEY2  ,
     >            NBE   )
      IF(NBE .EQ. 1) THEN
        NRCELA=NRCELA+1
        READ(IFT16) RKEY1,RKEY2,NBE,NID
        IF(NID .NE. NGMTR) THEN
          WRITE(IOUT,9000) NAMSBR,NGMTR,NRCELA,NID
          CALL XABORT(NAMSBR//': INVALID CELLAV STRUCTURE')
        ENDIF
        GO TO 100
      ELSE IF(NBE .EQ. -1) THEN
        GO TO 105
      ELSE
        WRITE(IOUT,9001) NAMSBR,1,NBE
        CALL XABORT(NAMSBR//': INVALID CELLAV STRUCTURE')
      ENDIF
 105  CONTINUE
*----
*  FIND THE NUMBER OF SETS OF REGION RECORD NRREGI
*  BASED ON THE PRESENCE OF REGION,DESCRIPTON KEYS
*  ALSO TEST FOR NZONE, NGMTR AND NREGON CONSISTENCY
*----
      REWIND(IFT16)
      NRREGI=0
      TKEY1='REGION    '
      TKEY2='DESCRIPTON'
 110  CONTINUE
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1  ,TKEY2  ,
     >            NBE   )
      IF(NBE .GT. 0) THEN
        NRREGI=NRREGI+1
        READ(IFT16) RKEY1,RKEY2,NBE,NID
        IF(NID .NE. NZONE ) THEN
          WRITE(IOUT,9010) NAMSBR,NZONE,NRREGI,NID
          CALL XABORT(NAMSBR//': INVALID REGION STRUCTURE')
        ENDIF
        READ(IFT16) RKEY1,RKEY2,NBE,NID,NJD
        IF(NID .NE. NREGON ) THEN
          WRITE(IOUT,9010) NAMSBR,NREGON,NRREGI,NID
          CALL XABORT(NAMSBR//': INVALID REGION STRUCTURE')
        ENDIF
        IF(NJD .NE. NGMTR ) THEN
          WRITE(IOUT,9010) NAMSBR,NGMTR,NRREGI,NJD
          CALL XABORT(NAMSBR//': INVALID REGION STRUCTURE')
        ENDIF
        GO TO 110
      ELSE
        GO TO 115
      ENDIF
 115  CONTINUE
*----
*  PROCESS PRINT LEVEL
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR,
     >  SUBTIT(113:240),SUBTIT(1:80),SUBTIT(89:104)
        WRITE(IOUT,6010) NEL,NG,NGMTR,NMATZ,NM,MTRMSH,
     >                   NFPR,NZONE,NREGON,NGREAC,NRCELA,NRREGI
        WRITE(IOUT,6001)
      ENDIF
      RETURN
*----
*  PRINT FORMAT
*----
 6000 FORMAT(1X,5('*'),' OUTPUT FROM ',A6,1X,5('*')/
     >       6X,'CONTENTS OF TAPE16 FILE :'/A128/
     >       6X,'WIMS-AECL VERSION = ',A80/
     >       6X,'LIBRARY VERSION   = ',A16)
 6001 FORMAT(1X,30('*'))
 6010 FORMAT(6X,'DIMENSIONING DATA '/
     >       6X,'NEL    : NB. ISOTOPES             = ',I10/
     >       6X,'NG     : NB. GROUPS               = ',I10/
     >       6X,'NGMTR  : NB. MAIN TRANSPORT GROUP = ',I10/
     >       6X,'NMATZ  : NB. MIXTURES             = ',I10/
     >       6X,'NM     : NB. BURNABLE MATERIALS   = ',I10/
     >       6X,'MTRMSH : NB. TRANSPORT MESH POINTS= ',I10/
     >       6X,'NFPR   : NB. FUEL PIN RINGS       = ',I10/
     >       6X,'NZONE  : NB. ZONES                = ',I10/
     >       6X,'NREGON : NB. EDIT REGIONS         = ',I10/
     >       6X,'NGREAC : NB. EDIT GROUPS          = ',I10/
     >       6X,'NRCELA : NB. CELLAV  RECORDS      = ',I10/
     >       6X,'NRREGI : NB. REGION  RECORDS      = ',I10)
*----
*  WARNING FORMAT
*----
 8000 FORMAT(1X,A6,1X,6('*'),' WARNING ',6('*')/
     >       8X,'RECORD WITH KEYS ',2(A10,2X),'NOT FOUND'/
     >       8X,'USE DEFAULT VALUE FOR ',A6,' = ',I10/
     >       8X,21('*'))
*----
*  ABORT FORMAT
*----
 9000 FORMAT(1X,A6,1X,7('*'),' ERROR ',7('*')/
     >       8X,6X,' NUMBER OF MAIN TRANSPORT GROUP ',I10/
     >       8X,I6,' CELLAV NGROUPS RECORD GIVES    ',I10/
     >       8X,21('*'))
 9001 FORMAT(1X,A6,1X,7('*'),' ERROR ',7('*')/
     >       8X,' NB ELEMENT ALLOWED ON CELLAV NGROUPS  ',I10/
     >       8X,' NB ELEMENT READ ON CELLAV NGROUPS      ',I10/
     >       8X,21('*'))
 9010 FORMAT(1X,A6,1X,7('*'),' ERROR ',7('*')/
     >       8X,6X,' NUMBER OF ZONES     ',I10/
     >       8X,I6,' REGION RECORD ',I10,' GIVES ',I10/
     >       8X,21('*'))
      END
