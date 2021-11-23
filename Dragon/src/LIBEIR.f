*DECK LIBEIR
      SUBROUTINE LIBEIR(MAXR,NEL,NMDEPL,ITNAM,ITZEA,KPAX,BPAX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read depletion data on input file.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert and G. Marleau.
*
*Parameters: input
* MAXR    number of reaction types.
* NEL     number of isotopes on library.
* NMDEPL  names of reactions:
*           NMDEPL(1)='DECAY'; NMDEPL(2)='NFTOT';
*           NMDEPL(3)='NG'   ; NMDEPL(4)='N2N';
*           etc.
*
*Parameters: output
* ITNAM   reactive isotope names in chain.
* ITZEA   6-digit nuclide identifier:
*         atomic number z*10000 (digits) + mass number a*10 +
*         energy state (0 = ground state, 1 = first state, etc.).
* KPAX    complete reaction type matrix.
* BPAX    complete branching ratio matrix.
*
*Comments:
*  INPUT FORMAT
*    CHAIN
*    [[ hnamson [ izea ]
*      [ [[ { DECAY    constant |
*            reaction [energy] } ]]  ]
*    [ { STABLE |
*       FROM  [[ { DECAY | reaction }
*          [[ yield hnampar ]] ]] } ]
*    ]]
*    ENDCHAIN
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXR,NEL,ITNAM(3,NEL),ITZEA(NEL),KPAX(NEL+MAXR,NEL)
      CHARACTER NMDEPL(MAXR)*8
      REAL BPAX(NEL+MAXR,NEL)
*----
*  INPUT FILE PARAMETERS
*----
      CHARACTER TEXT12*12
      INTEGER KNADPL(3)
      DOUBLE PRECISION DBLINP
*----
*  INTERNAL PARAMETERS
*   KFISSP : DRAGON FISSION PRODUCT FLAG = 2
*            POSITION OF NFTOT IN NMDEPL
*----
      INTEGER KFISSP
      PARAMETER (KFISSP=2)
      INTEGER INDIC,NITMA,NDEPL,IEL,JEL,IDEPL,INTG,IREAC,ISOT,JREL,JDEPL
      REAL FLOTT,RRAT
*----
*  READ LIST OF ISOTOPES AND PROPERTIES
*----
      TEXT12=' '
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
      IF(INDIC.NE.3.OR.TEXT12.NE.'CHAIN')
     >  CALL XABORT('LIBEIR: KEYWORD CHAIN MISSING')
      NDEPL=0
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
      DO 100 IEL=1,NEL
*----
*  EXIT IF ENDCHAIN READ
*----
        IF(TEXT12.EQ.'ENDCHAIN') GO TO 105
*----
*  ISOTOPE NAME READ
*  IF NAME ALREADY EXISTS SELECT ISOTOPE NUMBER
*  IF NAME NOT DEFINED ADD TO ISOTOPE LIST
*----
        IF(INDIC.NE.3)
     >    CALL XABORT('LIBEIR: ISOTOPE NAME HNAMSON MISSING')
        READ(TEXT12,'(3A4)') KNADPL(1),KNADPL(2),KNADPL(3)
        DO 110 JEL=1,NDEPL
          IF(KNADPL(1).EQ.ITNAM(1,JEL).AND.
     >       KNADPL(2).EQ.ITNAM(2,JEL).AND.
     >       KNADPL(3).EQ.ITNAM(3,JEL)     ) THEN
            IDEPL=JEL
            GO TO 115
          ENDIF
 110    CONTINUE
        NDEPL=NDEPL+1
        IF(NDEPL.GT.NEL)
     >    CALL XABORT('LIBEIR: TO MANY ISOTOPES')
        IDEPL=NDEPL
        ITNAM(1,IDEPL)=KNADPL(1)
        ITNAM(2,IDEPL)=KNADPL(2)
        ITNAM(3,IDEPL)=KNADPL(3)
 115    CONTINUE
*----
*  READ IZEA
*----
        CALL REDGET(INDIC,INTG,FLOTT,TEXT12,DBLINP)
        IF(INDIC.EQ.1) THEN
           ITZEA(IDEPL)=INTG
           CALL REDGET(INDIC,INTG,FLOTT,TEXT12,DBLINP)
        ELSE
           ITZEA(IDEPL)=0
        ENDIF
*----
*  LOOP OVER ALL PARAMETERS ASSOCIATED WITH SON ISOTOPES
*----
 120    CONTINUE
        IF(INDIC.NE.3) CALL XABORT('LIBEIR: REACTION TYPE EXPECTED FOR'
     >  //' ISOTOPE '//TEXT12)
*----
*  IF KEYWORD IS 'FROM' READ LIST OF PARENT NUCLIDES
*----
        IF(TEXT12.EQ.'FROM') THEN
*----
*  LOOP OVER ALL PARAMETERS ASSOCIATED WITH PARENT ISOTOPES
*----
          CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
 130      CONTINUE
          IF(INDIC.NE.3)
     >      CALL XABORT('LIBEIR: REACTION TYPE EXPECTED.')
          DO 140 IREAC=1,MAXR
            RRAT=1.0
*----
*  TEST IF KEYWORD IS A REACTION
*----
            IF(TEXT12.EQ.NMDEPL(IREAC)) THEN
*----
*  READ LIST OF YIELD AND PARENT ISOTOPES
*----
              DO 150 JEL=1,NEL
*----
*  IF YIELD ABSENT GO TO TEST FOR NEW REACTION TYPE
*----
                CALL REDGET(INDIC,ISOT,RRAT,TEXT12,DBLINP)
                IF(INDIC.NE.2) GO TO 130
                CALL REDGET(INDIC,ISOT,FLOTT,TEXT12,DBLINP)
                IF(INDIC.NE.3)
     >            CALL XABORT('LIBEIR: ISOTOPE NAME hnampar MISSING')
*----
*  ISOTOPE NAME READ
*  IF NAME ALREADY EXISTS SELECT ISOTOPE NUMBER
*  IF NAME NOT DEFINED ADD TO ISOTOPE LIST
*----
                READ(TEXT12,'(3A4)') KNADPL(1),KNADPL(2),KNADPL(3)
                DO 160 JREL=1,NDEPL
                  IF(KNADPL(1).EQ.ITNAM(1,JREL).AND.
     >               KNADPL(2).EQ.ITNAM(2,JREL).AND.
     >               KNADPL(3).EQ.ITNAM(3,JREL)     ) THEN
                    JDEPL=JREL
                    GO TO 165
                  ENDIF
 160            CONTINUE
                NDEPL=NDEPL+1
                IF(NDEPL.GT.NEL) CALL XABORT('LIBEIR: TO MANY ISOTOPES')
                JDEPL=NDEPL
                ITNAM(1,JDEPL)=KNADPL(1)
                ITNAM(2,JDEPL)=KNADPL(2)
                ITNAM(3,JDEPL)=KNADPL(3)
 165            CONTINUE
                KPAX(IDEPL,JDEPL)=IREAC
                BPAX(IDEPL,JDEPL)=RRAT
 150          CONTINUE
              CALL XABORT('LIBEIR: TO MANY PARENT ISOTOPES')
            ENDIF
 140      CONTINUE
        ELSE IF(TEXT12.EQ.'STABLE') THEN
          DO 141 IREAC=1,MAXR
          IF(KPAX(NEL+IREAC,IDEPL).NE.0) KPAX(NEL+IREAC,IDEPL)=-9999
 141      CONTINUE
          CALL REDGET(INDIC,INTG,FLOTT,TEXT12,DBLINP)
*----
*  READ NEXT KEYWORD FOR THIS ISOTOPE
*----
        ELSE
          DO 170 IREAC=1,MAXR
            RRAT=0.0
            IF(TEXT12.EQ.NMDEPL(IREAC)) THEN
              CALL REDGET(INDIC,ISOT,RRAT,TEXT12,DBLINP)
              IF(INDIC.EQ.1) THEN
                CALL XABORT('LIBEIR: INVALID INTEGER')
              ELSE IF(INDIC.EQ.2) THEN
                CALL REDGET(INDIC,INTG,FLOTT,TEXT12,DBLINP)
              ENDIF
              KPAX(NEL+IREAC,IDEPL)=1
              BPAX(NEL+IREAC,IDEPL)=RRAT
*----
*  READ NEXT KEYWORD FOR THIS ISOTOPE
*----
              GO TO 120
            ENDIF
 170      CONTINUE
        ENDIF
 100  CONTINUE
      IF(INDIC.NE.3.OR.TEXT12.NE.'ENDCHAIN')
     >  CALL XABORT('LIBEIR: KEYWORD ENDCHAIN MISSING')
 105  CONTINUE
*----
*  FIND FISSION PRODUCTS
*----
      DO 200 IEL=1,NDEPL
        DO 210 JEL=1,NDEPL
          IF(KPAX(JEL,IEL).EQ.KFISSP) KPAX(NEL+KFISSP,JEL)=-1
 210    CONTINUE
 200  CONTINUE
*----
*  RETURN FROM LIBEIR
*----
      RETURN
      END
