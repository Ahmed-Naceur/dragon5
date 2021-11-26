*DECK PCREIR
      SUBROUTINE PCREIR(NMDEPL,MD2,NEL,ITNAM,ITZEA,KPAX,BPAX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read depletion data on input file. Based on LIBEIR.f routine in
* DRAGON.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* NMDEPL  names of reactions:
*           NMDEPL(1)='DECAY'; NMDEPL(2)='NFTOT';
*           NMDEPL(3)='NG'   ; NMDEPL(4)='N2N';
*           etc
* MD2     dimension of arrays ITNAM, ITZEA, KPAX and BPAX
*
*Parameters: output
* NEL     number of particularized isotopes including macro
* ITNAM   reactive isotope names in chain
* ITZEA   6-digit nuclide identifier
*         atomic number z*10000 (digits) + mass number a*10 +
*         energy state (0 = ground state, 1 = first state, etc.)
* KPAX    complete reaction type matrix
* BPAX    complete branching ratio matrix
*
*-----------------------------------------------------------------------
*
*----
*  INPUT FORMAT
*----
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
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER, PARAMETER::MAXR=12
      INTEGER MD2,NEL,ITNAM(3,MD2),ITZEA(MD2),KPAX(MD2+MAXR,MD2)
      CHARACTER NMDEPL(MAXR)*8
      REAL BPAX(MD2+MAXR,MD2)
*----
*  INPUT FILE PARAMETERS
*----
      CHARACTER TEXT12*12
      INTEGER KNADPL(2)
      DOUBLE PRECISION DBLINP
*----
*  INTERNAL PARAMETERS
*   KFISSP : FISSION PRODUCT FLAG = 2 (POSITION OF NFTOT IN NMDEPL)
*----
      INTEGER KFISSP
      PARAMETER (KFISSP=2)
      INTEGER INDIC,NITMA,IEL,JEL,IDEPL,INTG,IREAC,ISOT,JREL,JDEPL
      REAL FLOTT,RRAT
*
      NEL=0
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
  105 IF(INDIC.NE.3) CALL XABORT('PCREIR: CHARACTER DATA EXPECTED')
*----
*  EXIT IF ENDCHAIN READ
*----
      IF(TEXT12.EQ.'ENDCHAIN') GO TO 190
*----
*  ISOTOPE NAME READ
*  IF NAME ALREADY EXISTS SELECT ISOTOPE NUMBER
*  IF NAME NOT DEFINED ADD TO ISOTOPE LIST
*----
      IDEPL=0
      READ(TEXT12,'(2A4)') KNADPL(1),KNADPL(2)
      DO 110 JEL=1,NEL
        IF(KNADPL(1).EQ.ITNAM(1,JEL).AND.
     >     KNADPL(2).EQ.ITNAM(2,JEL)) THEN
          IDEPL=JEL
          GO TO 115
        ENDIF
 110  CONTINUE
      NEL=NEL+1
      IF(NEL.GT.MD2) CALL XABORT('PCREIR: MD2 OVERFLOW(1).')
      IDEPL=NEL
      ITNAM(1,NEL)=KNADPL(1)
      ITNAM(2,NEL)=KNADPL(2)
*----
*  READ IZEA
*----
 115  CALL REDGET(INDIC,INTG,FLOTT,TEXT12,DBLINP)
      IF(INDIC.EQ.1) THEN
         ITZEA(IDEPL)=INTG
         CALL REDGET(INDIC,INTG,FLOTT,TEXT12,DBLINP)
      ELSE
         ITZEA(IDEPL)=0
      ENDIF
*----
*  LOOP OVER ALL PARAMETERS ASSOCIATED WITH SON ISOTOPES
*----
 120  IF(INDIC.NE.3) CALL XABORT('PCREIR: REACTION TYPE EXPECTED FOR'
     > //' ISOTOPE '//TEXT12)
*----
*  IF KEYWORD IS 'FROM' READ LIST OF PARENT NUCLIDES
*----
      IF(TEXT12.EQ.'FROM') THEN
*----
*  LOOP OVER ALL PARAMETERS ASSOCIATED WITH PARENT ISOTOPES
*----
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
 130    IF(INDIC.NE.3) CALL XABORT('PCREIR: REACTION TYPE EXPECTED.')
        DO 140 IREAC=1,MAXR
          RRAT=1.0
*----
*  TEST IF KEYWORD IS A REACTION
*----
          IF(TEXT12.EQ.NMDEPL(IREAC)) THEN
*----
*  READ LIST OF YIELD AND PARENT ISOTOPES
*----
            JDEPL=0
            DO 150 JEL=1,MD2
*----
*  IF YIELD ABSENT GO TO TEST FOR NEW REACTION TYPE
*----
              CALL REDGET(INDIC,ISOT,RRAT,TEXT12,DBLINP)
              IF(INDIC.NE.2) GO TO 130
              CALL REDGET(INDIC,ISOT,FLOTT,TEXT12,DBLINP)
              IF(INDIC.NE.3)
     >          CALL XABORT('PCREIR: ISOTOPE NAME hnampar MISSING')
*----
*  ISOTOPE NAME READ
*  IF NAME ALREADY EXISTS SELECT ISOTOPE NUMBER
*  IF NAME NOT DEFINED ADD TO ISOTOPE LIST
*----
              READ(TEXT12,'(2A4)') KNADPL(1),KNADPL(2)
              DO 160 JREL=1,MD2
                IF(KNADPL(1).EQ.ITNAM(1,JREL).AND.
     >             KNADPL(2).EQ.ITNAM(2,JREL)) THEN
                  JDEPL=JREL
                  GO TO 165
                ENDIF
 160          CONTINUE
              NEL=NEL+1
              IF(NEL.GT.MD2) CALL XABORT('PCREIR: MD2 OVERFLOW(2).')
              JDEPL=NEL
              ITNAM(1,NEL)=KNADPL(1)
              ITNAM(2,NEL)=KNADPL(2)
 165          KPAX(IDEPL,JDEPL)=IREAC
              BPAX(IDEPL,JDEPL)=RRAT
 150        CONTINUE
            CALL XABORT('PCREIR: TO MANY PARENT ISOTOPES')
          ENDIF
 140    CONTINUE
      ELSE IF(TEXT12.EQ.'STABLE') THEN
        DO 141 IREAC=1,MAXR
        IF(KPAX(MD2+IREAC,IDEPL).NE.0) KPAX(MD2+IREAC,IDEPL)=-9999
 141    CONTINUE
        DO 142 IEL=1,MD2
        KPAX(IDEPL,IEL)=0
 142    CONTINUE
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
              CALL XABORT('PCREIR: INVALID INTEGER')
            ELSE IF(INDIC.EQ.2) THEN
              CALL REDGET(INDIC,INTG,FLOTT,TEXT12,DBLINP)
            ENDIF
            KPAX(MD2+IREAC,IDEPL)=1
            BPAX(MD2+IREAC,IDEPL)=RRAT
*----
*  READ NEXT KEYWORD FOR THIS ISOTOPE
*----
            GO TO 120
          ENDIF
 170    CONTINUE
      ENDIF
      GO TO 105
*----
*  FIND FISSION PRODUCTS
*----
 190  DO 200 IEL=1,MD2
        DO 210 JEL=1,MD2
          IF(KPAX(JEL,IEL).EQ.KFISSP) KPAX(MD2+KFISSP,JEL)=-1
 210    CONTINUE
 200  CONTINUE
      IF(NEL.NE.MD2) CALL XABORT('PCREIR: INVALID VALUE OF MD2.')
*----
*  RETURN FROM PCREIR
*----
      RETURN
      END
