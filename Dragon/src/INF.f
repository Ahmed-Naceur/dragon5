*DECK INF
      SUBROUTINE INF(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for the dragon information module to recover cle-2000 values
* from the xs libraries.
*
*Copyright:
* Copyright (C) 1995 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): R. Roy
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file.
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      INTEGER  IOUT
      CHARACTER NAMSBR*6
      PARAMETER (IOUT=6,NAMSBR='INF   ')
      CHARACTER TEXT12*12,HNAMIS(3)*8,TEXTT(3)*12
      CHARACTER TEXT64*64,CFILNA*64
      LOGICAL LTMP,LPUR,LENR,LISO,LPRES
      INTEGER IPARM,I,IPISO(3),ITYPL
      INTEGER ITYP,NITMA,NOUT,NCARS
      INTEGER ITYPE,ILOOP,NBISO,IPRINT
      DOUBLE  PRECISION DFLOTT
      REAL    FLOTT,RBASE(3),AWR(3),PRES
      REAL    TEMPC,TEMPK,PURWGT,PURATM,ENRWGT,ENRATM,TOTMU
      IF(NENTRY.NE.0)THEN
        CALL XABORT(NAMSBR//': NO DATA STRUCTURE EXPECTED')
      ENDIF
      CFILNA=' '
      IPRINT= 1
      ITYPE=  2
      IPARM=  0
      NBISO=  0
      LTMP=.FALSE.
      LPRES=.FALSE.
      LPUR=.FALSE.
      LENR=.FALSE.
      LISO=.FALSE.
      NOUT= 1
      ITYPL=0
      ENRWGT=0.0
      PRES=0.0
      NCARS=0
      DO ILOOP=1,3
        TEXTT(ILOOP)='            '
      ENDDO
   20 CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
      IF( ITYP.NE.3 )CALL XABORT(NAMSBR//': CHARACTER DATA EXPECTED.')
      TEXTT(1)=TEXT12
      IF(TEXTT(1).EQ.';' )THEN
         GO TO 40
      ELSEIF(TEXTT(1).EQ.'EDIT' )THEN
         CALL REDGET(ITYP,IPRINT,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.1 )CALL XABORT(NAMSBR//': INTEGER EXPECTED.')
         DO I=1,NENTRY
           WRITE(IOUT,*) HENTRY(I),IENTRY(I),JENTRY(I)
           IF(IENTRY(I).LE.2) CALL LCMLIB(KENTRY(I))
         ENDDO
         GO TO 20
      ELSEIF(TEXTT(1).EQ.'LIB:' )THEN
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.3 ) CALL XABORT
     >     (NAMSBR//': LIBRARY TYPE NOT SPECIFIED FOLLOWING LIB:')
         IF(TEXT12.EQ.'WIMSAECL') THEN
           ITYPL=1
         ELSE IF(TEXT12.EQ.'WIMSD4') THEN
           ITYPL=2
         ELSE IF(TEXT12.EQ.'APLIB1') THEN
           ITYPL=3
         ELSE IF(TEXT12.EQ.'DRAGON') THEN
           ITYPL=4
         ELSE IF(TEXT12.EQ.'MATXS ') THEN
           ITYPL=5
         ELSE IF(TEXT12.EQ.'MATXS2') THEN
           ITYPL=6
         ELSE IF(TEXT12.EQ.'NDAS') THEN
           ITYPL=7
         ELSE IF(TEXT12.EQ.'WIMSE') THEN
           ITYPL=8
         ELSE
           CALL XABORT(NAMSBR//': ILLEGAL LIBRARY TYPE FOLLOWING LIB:')
         ENDIF
          CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.3.OR.TEXT12.NE.'FIL:' )
     >      CALL XABORT(NAMSBR//': *FIL:* EXPECTED.')
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT64,DFLOTT)
         IF( ITYP.NE.3 )CALL XABORT(NAMSBR//': LIBRARY NAME EXPECTED.')
         CFILNA= TEXT64
         GO TO 20
      ELSEIF(TEXTT(1).EQ.'TMP:' )THEN
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.2 )
     >      CALL XABORT(NAMSBR//': TEMPERATURE EXPECTED.')
         TEMPK = FLOTT
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.3)
     >      CALL XABORT(NAMSBR//': *C* OR *K* UNIT EXPECTED.')
         IF(TEXT12.EQ.'C')THEN
            TEMPK = TEMPK + 273.15
         ELSEIF( TEXT12.NE.'K' )THEN
            CALL XABORT(NAMSBR//': *C* OR *K* UNIT EXPECTED.')
         ENDIF
         TEMPC = TEMPK-273.15
         LTMP=.TRUE.
         GO TO 20
      ELSEIF(TEXTT(1).EQ.'PRES:' )THEN
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.2 )
     >      CALL XABORT(NAMSBR//': PRESSURE EXPECTED (Pa).')
         PRES = FLOTT
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.3)
     >      CALL XABORT(NAMSBR//
     >      ': *Pa*, *kPa*, *MPa* or *bar* UNITS EXPECTED.')
         IF( TEXT12.EQ.'kPa' ) THEN
            PRES = PRES*1000.0
         ELSE IF( TEXT12.EQ.'bar' ) THEN
            PRES = PRES*100000.0
         ELSE IF( TEXT12.EQ.'MPa' ) THEN
            PRES = PRES*1000000.0
         ELSE IF( TEXT12.NE.'Pa' ) THEN
            CALL XABORT(NAMSBR//
     >      ': *Pa*, *kPa*, *MPa* or *bar* UNITS EXPECTED.')
         ENDIF
         LPRES=.TRUE.
         GO TO 20
      ELSEIF(TEXTT(1).EQ.'PUR:' )THEN
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.2 )
     >      CALL XABORT(NAMSBR//': PURITY EXPECTED.')
         PURWGT = FLOTT
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.3)CALL XABORT(NAMSBR//': *ATM%* OR *WGT%* '
     >                           //'UNIT EXPECTED.')
         IF( TEXT12.EQ.'ATM%' )THEN
            PURATM= PURWGT
            PURWGT= 100.0/(1. + 0.8994866*(100./PURATM - 1.))
         ELSEIF( TEXT12.EQ.'WGT%' )THEN
            PURATM= 100.0/(1. + 1.1117435*(100./PURWGT - 1.))
         ELSE
            CALL XABORT(NAMSBR//': *ATM%* OR *WGT%* UNIT EXPECTED.')
         ENDIF
         LPUR=.TRUE.
         GO TO 20
      ELSEIF(TEXTT(1).EQ.'ENR:' )THEN
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.2 )
     >      CALL XABORT(NAMSBR//': ENRICHMENT EXPECTED.')
         ENRWGT = FLOTT
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.3)CALL XABORT(NAMSBR//': *ATM%* OR *WGT%* '
     >                           //'UNIT EXPECTED.')
         IF( TEXT12.EQ.'ATM%' )THEN
            ENRATM= ENRWGT
            ENRWGT= 100.0/(1. + 1.01279335*(100./ENRATM - 1.))
         ELSEIF( TEXT12.EQ.'WGT%' )THEN
            ENRATM= 100.0/(1. + 0.98736825*(100./ENRWGT - 1.))
         ELSE
            CALL XABORT(NAMSBR//': *ATM%* OR *WGT%* UNIT EXPECTED.')
         ENDIF
         LENR=.TRUE.
         GO TO 20
      ELSEIF(TEXTT(1).EQ.'ISO:' )THEN
         IF( NBISO.NE.0 )
     >      CALL XABORT(NAMSBR//': PREVIOUS ISOTOPES NOT USED.')
         CALL REDGET(ITYP,NBISO,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.1 )
     >      CALL XABORT(NAMSBR//': NUMBER OF ISOTOPES EXPECTED.')
         IF( NBISO.LE.0.OR.NBISO.GT.3 )
     >      CALL XABORT(NAMSBR//
     >      ': NB OF ISOTOPES MUST BE BETWEEN 1 AND 3.')
         DO I=1,NBISO
            CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
            IF( ITYP.NE.3)
     >         CALL XABORT(NAMSBR//': ISOTOPE NAME EXPECTED.')
            HNAMIS(I)= TEXT12(1:8)
         ENDDO
         IF(ITYPL.EQ.1)THEN
            CALL INFWIM(CFILNA,IPRINT,NBISO,HNAMIS,AWR)
         ELSE IF(ITYPL.EQ.2)THEN
            CALL INFWD4(CFILNA,4,IPRINT,NBISO,HNAMIS,AWR)
         ELSE IF(ITYPL.EQ.3)THEN
            CALL INFAPL(CFILNA,IPRINT,NBISO,HNAMIS,AWR)
         ELSE IF(ITYPL.EQ.4)THEN
            CALL INFDRA(CFILNA,IPRINT,NBISO,HNAMIS,AWR)
         ELSE IF(ITYPL.EQ.5)THEN
            CALL INFTR1(CFILNA,IPRINT,NBISO,HNAMIS,AWR)
         ELSE IF(ITYPL.EQ.6)THEN
            CALL INFTR2(CFILNA,IPRINT,NBISO,HNAMIS,AWR)
         ELSE IF(ITYPL.EQ.7)THEN
            CALL INFNDA(CFILNA,IPRINT,NBISO,HNAMIS,AWR)
         ELSE IF(ITYPL.EQ.8)THEN
            CALL INFWD4(CFILNA,5,IPRINT,NBISO,HNAMIS,AWR)
         ENDIF
         GO TO 20
      ELSEIF(TEXTT(1).EQ.'CALC' )THEN
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.3 )CALL XABORT(NAMSBR//': *DENS* EXPECTED.')
         TEXTT(2)=TEXT12
         IF(TEXTT(2).EQ.'DENS' ) THEN
            CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(ITYP.NE.3) CALL XABORT(NAMSBR//': *WATER* EXPECTED.')
            TEXTT(3)=TEXT12
            NCARS=3
            IF(TEXTT(3).EQ.'WATER' ) THEN
               IF( .NOT.LTMP )
     >            CALL XABORT(NAMSBR//': NO TEMPERATURE GIVEN.')
               IF( .NOT.LPUR )
     >            CALL XABORT(NAMSBR//': NO PURITY      GIVEN.')
               IF(LPRES) WRITE(IOUT,9000) NAMSBR
               CALL INFWAT(TEMPC,PURWGT,RBASE(1))
               NOUT= 1
            ELSEIF(TEXTT(3).EQ.'PWATER' ) THEN
               IF( .NOT.LTMP )
     >            CALL XABORT(NAMSBR//': NO TEMPERATURE GIVEN.')
               IF( .NOT.LPUR )
     >            CALL XABORT(NAMSBR//': NO PURITY      GIVEN.')
               IF( .NOT.LPRES) THEN
                 CALL INFPSA(IPRINT,TEMPK,PURWGT,PRES)
               ENDIF
               CALL INFWAN(TEMPK,PURWGT,PRES,RBASE(1))
               NOUT= 1
            ELSE
               CALL XABORT(NAMSBR//': *WATER* or *PWATER* EXPECTED.')
            ENDIF
         ELSEIF(TEXTT(2).EQ.'WGT%')THEN
            CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
            IF( ITYP.NE.3)
     >         CALL XABORT(NAMSBR//
     >         ': *D2O*, *H2O*, *UO2* OR *THO2* EXPECTED.')
            TEXTT(3)=TEXT12
            NCARS=4
            IF( NBISO.NE.3 )
     >         CALL XABORT(NAMSBR//': NB OF ISOTOPES MUST BE 3.')
            IPISO(1)= 0
            IPISO(2)= 0
            IPISO(3)= 0
            IF(TEXTT(3).EQ.'UO2' )THEN
               IF( .NOT.LENR )
     >            CALL XABORT(NAMSBR//': NO ENRICHMENT  GIVEN.')
               DO I=1,NBISO
                  IF( 15.8.LT.AWR(I).AND.AWR(I).LT.16.2 )THEN
                     IPISO(1)= I
                  ELSEIF( 234.8.LT.AWR(I).AND.AWR(I).LT.235.2 )THEN
                     IPISO(2)= I
                  ELSEIF( 237.8.LT.AWR(I).AND.AWR(I).LT.238.2 )THEN
                     IPISO(3)= I
                  ELSE
                     CALL XABORT(NAMSBR//': NOT A U5,U8 OR O ISOTOPE.')
                  ENDIF
               ENDDO
               IF( IPISO(1)*IPISO(2)*IPISO(3).EQ.0 )
     >            CALL XABORT(NAMSBR//': MISSING ONE OF TH2,U3 OR O.')
               RBASE(IPISO(2))=       ENRWGT
               RBASE(IPISO(3))= 100.- ENRWGT
               RBASE(IPISO(1))= 2.*AWR(IPISO(1))*
     > (RBASE(IPISO(2))/AWR(IPISO(2))+ RBASE(IPISO(3))/AWR(IPISO(3)))
               TOTMU=  RBASE(IPISO(1))+RBASE(IPISO(2))+RBASE(IPISO(3))
               RBASE(IPISO(1))= 100.*RBASE(IPISO(1))/TOTMU
               RBASE(IPISO(2))= 100.*RBASE(IPISO(2))/TOTMU
               RBASE(IPISO(3))= 100.*RBASE(IPISO(3))/TOTMU
             ELSEIF(TEXTT(3).EQ.'THO2' )THEN
               IF( .NOT.LENR )
     >            CALL XABORT(NAMSBR//': NO ENRICHMENT  GIVEN.')
               DO I=1,NBISO
                  IF( 15.8.LT.AWR(I).AND.AWR(I).LT.16.2 )THEN
                     IPISO(1)= I
                  ELSEIF( 232.8.LT.AWR(I).AND.AWR(I).LT.233.2 )THEN
                     IPISO(2)= I
                  ELSEIF( 231.8.LT.AWR(I).AND.AWR(I).LT.232.2 )THEN
                     IPISO(3)= I
                  ELSE
                     CALL XABORT(NAMSBR//
     >               ': NOT A TH2,U3 OR O ISOTOPE.')
                  ENDIF
               ENDDO
               IF( IPISO(1)*IPISO(2)*IPISO(3).EQ.0 )
     >            CALL XABORT(NAMSBR//': MISSING ONE OF TH2,U3 OR O.')
               RBASE(IPISO(2))=       ENRWGT
               RBASE(IPISO(3))= 100.- ENRWGT
               RBASE(IPISO(1))= 2.*AWR(IPISO(1))*
     > (RBASE(IPISO(2))/AWR(IPISO(2))+ RBASE(IPISO(3))/AWR(IPISO(3)))
               TOTMU=  RBASE(IPISO(1))+RBASE(IPISO(2))+RBASE(IPISO(3))
               RBASE(IPISO(1))= 100.*RBASE(IPISO(1))/TOTMU
               RBASE(IPISO(2))= 100.*RBASE(IPISO(2))/TOTMU
               RBASE(IPISO(3))= 100.*RBASE(IPISO(3))/TOTMU
           ELSEIF(TEXTT(3).EQ.'D2O' .OR. TEXTT(3).EQ.'H2O')THEN
               IF( .NOT.LPUR )
     >            CALL XABORT(NAMSBR//': NO PURITY      GIVEN.')
               DO I=1,NBISO
                  IF( 15.8.LT.AWR(I).AND.AWR(I).LT.16.2 )THEN
                     IPISO(1)= I
                  ELSEIF( 0.8.LT.AWR(I).AND.AWR(I).LT.1.2 )THEN
                     IPISO(2)= I
                  ELSEIF( 1.8.LT.AWR(I).AND.AWR(I).LT.2.2 )THEN
                     IPISO(3)= I
                  ELSE
                     CALL XABORT(NAMSBR//': NOT A H1,D2 OR O ISOTOPE.')
                  ENDIF
               ENDDO
               IF( IPISO(1)*IPISO(2)*IPISO(3).EQ.0 )
     >            CALL XABORT(NAMSBR//': MISSING ONE OF H1,D2 OR O.')
               RBASE(IPISO(2))= (100.-PURWGT)*2.*AWR(IPISO(2))/
     >                         (2.*AWR(IPISO(2))+AWR(IPISO(1)))
               RBASE(IPISO(3))=       PURWGT *2.*AWR(IPISO(3))/
     >                         (2.*AWR(IPISO(3))+AWR(IPISO(1)))
               RBASE(IPISO(1))= 100.-(RBASE(IPISO(2))+RBASE(IPISO(3)))
            ELSE
               CALL XABORT(NAMSBR//
     >         ': *D2O*, *H2O*, *UO2* OR *THO2* EXPECTED.')
            ENDIF
            NOUT= NBISO
            NBISO= 0
         ELSE
            CALL XABORT(NAMSBR//': *DENS* OR *WGT%* EXPECTED.')
         ENDIF
      ELSEIF(TEXTT(1).EQ.'GET' )THEN
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.3 )CALL XABORT(NAMSBR//': *MASS* EXPECTED.')
         TEXTT(2)=TEXT12
         NCARS=2
         IF(TEXTT(2).EQ.'MASS' ) THEN
            IF( NBISO.EQ.0 )
     >         CALL XABORT(NAMSBR//': ISOTOPE LIST NOT SPECIFIED.')
            NOUT= NBISO
            NBISO= 0
            DO ILOOP= 1, NOUT
               RBASE(ILOOP)= AWR(ILOOP)
            ENDDO
         ELSE
            CALL XABORT(NAMSBR//': *MASS* EXPECTED.')
         ENDIF
      ELSE
         CALL XABORT(NAMSBR//': '//TEXT12//' IS AN INVALID KEYWORD.')
      ENDIF
*----
*  PUT PARMS IN CLE-2000 REAL VARIABLES (WRITE MODE).
*----
      DO ILOOP= 1, NOUT
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.-ITYPE )THEN
            CALL XABORT(NAMSBR//': INVALID TYPE FOR OUTPUT VALUE')
         ELSEIF( IPRINT.GT.0 )THEN
            IF(NCARS .EQ. 2) THEN
              WRITE(IOUT,6000) NAMSBR,TEXTT(1),TEXTT(2),
     >        HNAMIS(ILOOP),RBASE(ILOOP)
            ELSE IF(NCARS .EQ. 3) THEN
              WRITE(IOUT,6001) NAMSBR,TEXTT(1),TEXTT(2),TEXTT(3),
     >        RBASE(ILOOP)
            ELSE IF(NCARS .EQ. 4) THEN
              WRITE(IOUT,6002) NAMSBR,TEXTT(1),TEXTT(2),TEXTT(3),
     >        HNAMIS(ILOOP),RBASE(ILOOP)
            ENDIF
         ENDIF
         CALL REDPUT(ITYPE,NITMA,RBASE(ILOOP),TEXT12,DFLOTT)
      ENDDO
      GO TO 20
   40 CONTINUE
      RETURN
*
 6000 FORMAT(A6,': ',2(A12,1X),'Isotope ',A8,' <- ',1P,E15.7)
 6001 FORMAT(A6,': ',3(A12,1X),' <- ',1P,E15.7)
 6002 FORMAT(A6,': ',3(A12,1X),'Isotope ',A8,' <- ',1P,E15.7)
 9000 FORMAT('*****  WARNING in ',A6,'*****'/
     >       ' Pressure is not used with option -WATER-'/
     >       ' For pressure dependence use option -PWATER-')
      END
