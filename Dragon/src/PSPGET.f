*DECK PSPGET
      SUBROUTINE PSPGET(IPRINT,ITYPE,ICOLR,NGROUP,NGT,ICOND)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read PSP: module input data.
*
*Copyright:
* Copyright (C) 1999 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPRINT  print level.                           
* ITYPE   type of graphic:                        
*         =0 color per region number;
*         =1 color per material;
*         =2 color for flux (one group);
*         =3 color for flux (multigroup);
*         =4 color per material for homogenizatION (HMIX);
*         =5 color for mode (one group);
*         =6 color for mode (multigroup).
* ICOLR   color set used:                         
*         = -4 fill hsb with no-contour;
*         = -3 fill cmyk with no-contour;
*         = -2 fill rgb with no-contour;
*         = -1 fill bw with no-contour;
*         =  0 no fill contour only;
*         =  1 fill bw and contour;
*         =  2 fill rgb and contour;
*         =  3 fill cmyk and contour;
*         =  4 fill hsb and contour.
* NGROUP  number of groups for flux.
* NGT     number of condensed groups for flux.
* ICOND   upper group condensation limit.
*
*Comments:
*  Input instructions:
*     [ EDIT iprint ]
*     [ FILL  { NONE | GRAY | RGB | CMYK | HSB }  [ NOCONTOUR ]  ]
*     [ TYPE  { REGION | MIXTURE | FLUX | HMIX | 
*               MGFLUX (icond(i),i=1,ngt) }  ]       ;
*     DEFAULT:
*       IPRINT = 1 -> EDIT  1
*       ITYPE  = 0 -> PER REGION NUMBER
*       ICOLR  = 4 -> FILL HSB WITH CONTOUR
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='PSPGET')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          IPRINT,ITYPE,ICOLR,NGROUP,NGT
      INTEGER          ICOND(NGROUP)
*----
*  REDGET INPUT VARIABLES
*----
      INTEGER          ITYPLU,INTLIR
      CHARACTER        CARLIR*12
      REAL             REALIR
      DOUBLE PRECISION DBLLIR
*----
*  LOCAL PARAMETERS
*----
      INTEGER          ICOL,ITY,ICONT,IGT
*----
*  READ OPTIONS
*----
      IPRINT=1
      ICOL=4
      ITY=0
      ICONT=1
 100  CONTINUE
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
 101  CONTINUE
      IF(ITYPLU .EQ. 10) THEN
        GO TO 105
      ELSE IF(ITYPLU .NE. 3) THEN
        CALL XABORT(NAMSBR//': ERROR -> CHARACTER VARIABLE EXPECTED')
      ENDIF
      IF(CARLIR(1:1) .EQ. ';' ) THEN
        GO TO 105
      ELSE IF(CARLIR .EQ. 'EDIT' ) THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1 ) GO TO 101
        IPRINT=INTLIR
      ELSE IF(CARLIR(1:4) .EQ. 'FILL' ) THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3 ) GO TO 101
        IF(CARLIR .EQ. 'NONE') THEN
          ICOL=0
        ELSE IF(CARLIR .EQ. 'GRAY') THEN
          ICOL=1
        ELSE IF(CARLIR .EQ. 'RGB') THEN
          ICOL=2
        ELSE IF(CARLIR .EQ. 'CMYK') THEN
          ICOL=3
        ELSE IF(CARLIR .EQ. 'HSB') THEN
          ICOL=4
        ELSE
          CALL XABORT(NAMSBR//': ILEGAL FILL KEYWORD '//CARLIR//
     >    'KEYWORD EXPECTED: NONE, GRAY, RGB, CMYK, HSB')
        ENDIF
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3 ) GO TO 101
        IF(CARLIR(1:4) .EQ. 'NOCO') THEN
          ICONT=0
        ELSE
          GO TO 101
        ENDIF
      ELSE IF(CARLIR(1:4) .EQ. 'TYPE' ) THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(CARLIR(1:4) .EQ. 'REGI') THEN
          ITY=0
        ELSE IF(CARLIR(1:4) .EQ. 'MIXT') THEN
          ITY=1
        ELSE IF(CARLIR(1:4) .EQ. 'FLUX') THEN
          ITY=2
          NGT=1
          ICOND(NGT)=NGROUP
        ELSE IF(CARLIR(1:4) .EQ. 'MODE') THEN
          ITY=5
          NGT=1
          ICOND(NGT)=NGROUP
        ELSE IF(CARLIR(1:4) .EQ. 'MGFL') THEN
          ITY=3
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) THEN
            NGT=NGROUP
            DO IGT=1,NGT
              ICOND(IGT)=IGT
            ENDDO
            GO TO 101
          ENDIF
          NGT=0
          DO IGT=1,NGROUP
            NGT=NGT+1
            IF(INTLIR .LT. 1 .OR. INTLIR .GT. NGROUP)
     >CALL XABORT(NAMSBR//': illegal group condensation number')
            IF(IGT .GT. 1) THEN
              IF(INTLIR .LE. ICOND(IGT-1))
     >CALL XABORT(NAMSBR//': group numbers must be increasing')
            ENDIF
            ICOND(IGT)=INTLIR
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU .NE. 1 ) THEN
              IF(ICOND(IGT) .NE. NGROUP) THEN
                NGT=NGT+1
                ICOND(NGT)=NGROUP
               ENDIF
              GO TO 101
            ENDIF
          ENDDO
        ELSE IF(CARLIR(1:4) .EQ. 'MGMD') THEN
          ITY=6
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) THEN
            NGT=NGROUP
            DO IGT=1,NGT
              ICOND(IGT)=IGT
            ENDDO
            GO TO 101
          ENDIF
          NGT=0
          DO IGT=1,NGROUP
            NGT=NGT+1
            IF(INTLIR .LT. 1 .OR. INTLIR .GT. NGROUP)
     >CALL XABORT(NAMSBR//': illegal group condensation number')
            IF(IGT .GT. 1) THEN
              IF(INTLIR .LE. ICOND(IGT-1))
     >CALL XABORT(NAMSBR//': group numbers must be increasing')
            ENDIF
            ICOND(IGT)=INTLIR
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU .NE. 1 ) THEN
              IF(ICOND(IGT) .NE. NGROUP) THEN
                NGT=NGT+1
                ICOND(NGT)=NGROUP
               ENDIF
              GO TO 101
            ENDIF
          ENDDO
        ELSE IF(CARLIR(1:4) .EQ. 'HMIX') THEN
          ITY=4
        ELSE
          CALL XABORT(NAMSBR//': ILEGAL TYPE KEYWORD '//CARLIR//
     >    'KEYWORD EXPECTED: REGION, MIXTURE, FLUX, MGFLUX')
        ENDIF
      ELSE
*----
*  INVALID OPTION
*----
        CALL XABORT(NAMSBR//': ILEGAL MAIN KEYWORD '//CARLIR//
     >  'KEYWORD EXPECTED: FILL, TYPE, EDIT OR ; ')
      ENDIF
      GO TO 100
 105  CONTINUE
*----
*  TEST READ OPTIONS
*  IF FILL = NONE (ICOL = 0) IMPOSE CONTOUR
*----
      IF(ICONT .EQ. 0) THEN
        ICOLR=-ICOL
      ELSE
        ICOLR=ICOL
      ENDIF
      ITYPE=ITY
*----
*  PRINT ECHO OF PSP OPTIONS THAT WILL BE USED
*----
      IF(IPRINT .GE. 1 ) THEN
        WRITE(IOUT,6000) IPRINT,ICOL,ITY,ICONT
      ENDIF
*----
*  RETURN
*----
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(' ------  PSP EXECUTION OPTIONS --------'/
     >       ' PRINT LEVEL = ',I8                    /
     >       ' COLOR       = ',I8                    /
     >       ' TYPE        = ',I8                    /
     >       ' CONTOUR     = ',I8                    /
     >       ' --------------------------------------')
      END
