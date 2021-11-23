*DECK XDRCAS
      SUBROUTINE XDRCAS(DIR,TEXT)
*
*-----------------------------------------------------------------------
*
* CONVERT A LOWER-CASE CHARACTER VARIABLE TO UPPER CASE OR
* UPPER CASE CHARACTER VARIABLE TO LOWER-CASE
*
* INPUT/OUTPUT VARIABLE:
*  DIR  : DIRECTION OF CONVERSION
*         ='LOWTOUP' FOR LOWER TO UPPER
*         ='UPTOLOW' FOR UPPER TO LOWER
*  TEXT : CHARACTER VARIABLE TO BE CONVERTED.
*
*-----------------------------------------------------------------------
*
      CHARACTER    DIR*(*),TEXT*(*)
C----
C  LOCAL PARAMETERS
C----
      PARAMETER   (NCAR=26)
      INTEGER      LENTEX,ITEX,ICAR
      CHARACTER    LOWCAS(NCAR)*1,UPCAS(NCAR)*1
      SAVE         LOWCAS,UPCAS
      DATA LOWCAS /'a','b','c','d','e','f','g','h','i','j','k','l','m',
     >             'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      DATA UPCAS  /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     >             'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
      LENTEX=LEN(TEXT)
      IF(DIR.EQ.'LOWTOUP') THEN
        DO 100 ITEX=1,LENTEX
          DO 110 ICAR=1,NCAR
            IF(TEXT(ITEX:ITEX).EQ.LOWCAS(ICAR)) THEN
              TEXT(ITEX:ITEX)=UPCAS(ICAR)
              GO TO 115
            ENDIF
 110      CONTINUE
 115      CONTINUE
 100    CONTINUE
      ELSE IF (DIR.EQ.'UPTOLOW') THEN
        DO 200 ITEX=1,LENTEX
          DO 210 ICAR=1,NCAR
            IF(TEXT(ITEX:ITEX).EQ.UPCAS(ICAR)) THEN
              TEXT(ITEX:ITEX)=LOWCAS(ICAR)
              GO TO 215
            ENDIF
 210      CONTINUE
 215      CONTINUE
 200    CONTINUE
      ENDIF
      RETURN
      END
