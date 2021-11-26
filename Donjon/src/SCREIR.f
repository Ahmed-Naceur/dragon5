*DECK SCREIR
      SUBROUTINE SCREIR(NMDEPL,MY1,MY2,MD1,NEL,NOMIS,ADRY,NVTOT,VTOT,
     > YLDS,DECAY,ITNAM,ITZEA,KPAX,BPAX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read depletion data on input file. Based on LIBEIR.f routine in
* DRAGON.
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal
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
* MY1     first dimension of matrix YLDS
* MY2     number of particularized fission products
* MD1     number of types of radioactive decay reactions
* NEL     number of particularized isotopes including macro
* NOMIS   names of isotopes in chain
* ADRY    indices of fissile isotopes (positive values) and fission
*         products (negative values) in array YLDS
* NVTOT   number of Saphyb calls
* VTOT    volume of updated core per Saphyb call
* YLDS    fission yields
* DECAY   radioactive decay constants
*
*Parameters: output
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
*      [ [[ reaction [energy] ]]  ]
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
      INTEGER MY1,MY2,MD1,NEL,ADRY(NEL),NVTOT,ITNAM(3,NEL),
     > ITZEA(NEL),KPAX(NEL+MAXR,NEL)
      CHARACTER NMDEPL(MAXR)*8,NOMIS(NEL)*8
      REAL BPAX(NEL+MAXR,NEL)
      DOUBLE PRECISION VTOT(NVTOT),YLDS(MY1,MY2,NVTOT),
     > DECAY(MD1,NEL,NVTOT)
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
      CHARACTER HSMG*131
      INTEGER INDIC,NITMA,I0,IEL,JEL,IDEPL,INTG,IREAC,ISOT,JREL,JDEPL,
     > IY1,IY2,IV
      REAL FLOTT,RRAT
      DOUBLE PRECISION ZN,ZD
*----
*  READ LIST OF ISOTOPES AND PROPERTIES
*----
      DO 70 IY1=1,MY1
        IEL=0
        DO I0=1,NEL
          IF(ADRY(I0).EQ.IY1) THEN
*           IEL is a fissile isotope
            IEL=I0
          ENDIF
        ENDDO
        IF(IEL.EQ.0) GO TO 70
        DO 60 IY2=1,MY2
          JEL=0
          DO I0=1,NEL
            IF(-ADRY(I0).EQ.IY2) THEN
*             JEL is a fission fragment
              JEL=I0
            ENDIF
          ENDDO
          IF(JEL.EQ.0) GO TO 60
          KPAX(JEL,IEL)=KFISSP
          ZN=0.0D0
          ZD=0.0D0
          DO IV=1,NVTOT
            ZN=ZN+YLDS(IY1,IY2,IV)*VTOT(IV)
            ZD=ZD+VTOT(IV)
          ENDDO
          BPAX(JEL,IEL)=REAL(ZN/ZD)
   60   CONTINUE
   70 CONTINUE
      DO 100 IEL=1,NEL
        TEXT12=' '
        TEXT12(:8)=NOMIS(IEL)
        READ(TEXT12,'(3A4)') (ITNAM(I0,IEL),I0=1,3)
        KPAX(NEL+1,IEL)=1
        BPAX(NEL+1,IEL)=0.0
        DO 80 I0=1,MD1
        ZN=0.0D0
        ZD=0.0D0
        DO IV=1,NVTOT
          ZN=ZN+DECAY(I0,IEL,IV)*VTOT(IV)
          ZD=ZD+VTOT(IV)
        ENDDO
        BPAX(NEL+1,IEL)=BPAX(NEL+1,IEL)+REAL(ZN/ZD)*1.0E8
   80   CONTINUE
  100 CONTINUE
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
  105 IF(INDIC.NE.3) CALL XABORT('SCREIR: CHARACTER DATA EXPECTED')
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
      WRITE(HSMG,'(16HSCREIR: ISOTOPE ,2A4,24H IS MISSING AMONG PARTIC,
     > 35HULARIZED ISOTOPES OF THE SAPHYB(1).)') KNADPL(1),KNADPL(2)
      CALL XABORT(HSMG)
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
 120  IF(INDIC.NE.3) CALL XABORT('SCREIR: REACTION TYPE EXPECTED FOR'
     > //' ISOTOPE '//TEXT12)
*----
*  IF KEYWORD IS 'FROM' READ LIST OF PARENT NUCLIDES
*----
      IF(TEXT12.EQ.'FROM') THEN
*----
*  LOOP OVER ALL PARAMETERS ASSOCIATED WITH PARENT ISOTOPES
*----
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
 130    IF(INDIC.NE.3) CALL XABORT('SCREIR: REACTION TYPE EXPECTED.')
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
            DO 150 JEL=1,NEL
*----
*  IF YIELD ABSENT GO TO TEST FOR NEW REACTION TYPE
*----
              CALL REDGET(INDIC,ISOT,RRAT,TEXT12,DBLINP)
              IF(INDIC.NE.2) GO TO 130
              CALL REDGET(INDIC,ISOT,FLOTT,TEXT12,DBLINP)
              IF(INDIC.NE.3)
     >          CALL XABORT('SCREIR: ISOTOPE NAME hnampar MISSING')
*----
*  ISOTOPE NAME READ
*  IF NAME ALREADY EXISTS SELECT ISOTOPE NUMBER
*  IF NAME NOT DEFINED ADD TO ISOTOPE LIST
*----
              READ(TEXT12,'(2A4)') KNADPL(1),KNADPL(2)
              DO 160 JREL=1,NEL
                IF(KNADPL(1).EQ.ITNAM(1,JREL).AND.
     >             KNADPL(2).EQ.ITNAM(2,JREL)) THEN
                  JDEPL=JREL
                  GO TO 165
                ENDIF
 160          CONTINUE
              WRITE(HSMG,'(16HSCREIR: ISOTOPE ,2A4,16H IS MISSING AMON,
     >        43HG PARTICULARIZED ISOTOPES OF THE SAPHYB(2).)')
     >        KNADPL(1),KNADPL(2)
              CALL XABORT(HSMG)
 165          KPAX(IDEPL,JDEPL)=IREAC
              BPAX(IDEPL,JDEPL)=RRAT
 150        CONTINUE
            CALL XABORT('SCREIR: TO MANY PARENT ISOTOPES')
          ENDIF
 140    CONTINUE
      ELSE IF(TEXT12.EQ.'STABLE') THEN
        DO 141 IREAC=1,MAXR
        IF(KPAX(NEL+IREAC,IDEPL).NE.0) KPAX(NEL+IREAC,IDEPL)=-9999
 141    CONTINUE
        DO 142 IEL=1,NEL
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
              CALL XABORT('SCREIR: INVALID INTEGER')
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
 170    CONTINUE
      ENDIF
      GO TO 105
*----
*  FIND FISSION PRODUCTS
*----
 190  DO 200 IEL=1,NEL
        DO 210 JEL=1,NEL
          IF(KPAX(JEL,IEL).EQ.KFISSP) KPAX(NEL+KFISSP,JEL)=-1
 210    CONTINUE
 200  CONTINUE
*----
*  RETURN FROM SCREIR
*----
      RETURN
      END
