*DECK DRVGRP
      SUBROUTINE DRVGRP(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* standard grep module to recover cle-2000 values in a linked list or
* in an xsm file.
*
*Copyright:
* Copyright (C) 2000 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): read-only type(VECTOR).
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
      INTEGER   NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR) KENTRY(NENTRY)
      CHARACTER HENTRY(*)*12
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXT12*12,HLIST*12,HNAME*12,TEXT72*72,TEXTIN*72,NAMT*12
      CHARACTER CLOGBG(5)*12,HSMG*131
      TYPE(C_PTR) IPDATA,JPDATA
      INTEGER   ITYLCM,IACTIO,I,J,K,N
      INTEGER   ITYP,ITYP1,NOUT,NSTP
      INTEGER   IBCHAR,NBCHAR,NRESID,IOFSET,ISET
      INTEGER   ITYPE,IOUT,ILENG,ILENG2,IPRINT
      INTEGER   NITMA, INTGMX, INTMIN, INTMAX, INDMIN, INDMAX
      PARAMETER  ( INTGMX= 2147483646 )
      DOUBLE PRECISION DFLOTT, DBLEMX, DBLMIN, DBLMAX, DBLMNV
      PARAMETER  ( DBLEMX= 1.D+100 )
      REAL      FLOTT, REALMX, RELMIN, RELMAX, RELMNV
      PARAMETER  ( REALMX= 1.E+30 )
      INTEGER   ISEEME(2),ITRANS
      REAL      ASEEME(2),ATRANS
      DOUBLE PRECISION DSEEME
      EQUIVALENCE ( DSEEME, ISEEME, ASEEME )
      EQUIVALENCE ( ITRANS, ATRANS )
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDATA,JDATA
      DATA CLOGBG
     1  / 'INTEGER' , 'REAL'    , 'STRING'   , 'DOUBLE', 'LOGICAL' /
*----
*  PARAMETER VALIDATION.
*----
      IACTIO=0
      INDMIN= 0
      INDMAX= 0
      NBCHAR= 0
      IF( NENTRY.NE.1 )THEN
         CALL XABORT('DRVGRP: MORE THAN ONE ENTRY')
      ELSEIF( IENTRY(1).NE.1.AND.IENTRY(1).NE.2 )THEN
         CALL XABORT('DRVGRP: RHS LINKED LIST '
     >             //'OR XSM FILE PARAMETER EXPECTED.')
      ELSEIF( JENTRY(1).NE.2 )THEN
         CALL XABORT('DRVGRP: RHS PARAMETER IN '
     >             //'READ-ONLY MODE EXPECTED.')
      ENDIF
*
      HLIST=HENTRY(1)
      IPDATA=KENTRY(1)
      I= 1
      IPRINT= 0
   20 CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
   25 CONTINUE
      IF( ITYP.NE.3 )CALL XABORT('DRVGRP: CHARACTER DATA EXPECTED.')
      IF( TEXT12.EQ.'EDIT' )THEN
         CALL REDGET(ITYP,IPRINT,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.1 )CALL XABORT('DRVGRP: NO INTEGER AFTER *EDIT*.')
      ELSE IF(TEXT12.EQ.'STEP') THEN
*        CHANGE THE HIERARCHICAL LEVEL ON THE LCM OBJECT.
         JPDATA=C_NULL_PTR
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(ITYP.NE.3) CALL XABORT('DRVGRP: CHARACTER DATA EXPECTED.')
         IF(TEXT12.EQ.'UP') THEN
            CALL REDGET(ITYP,NITMA,FLOTT,NAMT,DFLOTT)
            IF(ITYP.NE.3) CALL XABORT('DRVGRP: DIR-NAME EXPECTED.')
            JPDATA=LCMGID(IPDATA,NAMT)
         ELSE IF(TEXT12.EQ.'AT') THEN
            CALL REDGET(ITYP,NITMA,FLOTT,NAMT,DFLOTT)
            IF(ITYP.NE.1) CALL XABORT('DRVGRP: INTEGER EXPECTED.')
            JPDATA=LCMGIL(IPDATA,NITMA)
         ELSE
            CALL XABORT('DRVGRP: *UP* OR *AT* EXPECTED.')
         ENDIF
         IPDATA=JPDATA
      ELSEIF( TEXT12.EQ.'GETVAL'.OR.TEXT12.EQ.'MAXVAL'.OR.
     >        TEXT12.EQ.'MINVAL'.OR.TEXT12.EQ.'INDMAX'.OR.
     >        TEXT12.EQ.'INDMIN'.OR.TEXT12.EQ.'MEAN'  .OR.
     >        TEXT12.EQ.'TYPE'  .OR.TEXT12.EQ.'LENGTH')THEN
         IF(     TEXT12.EQ.'GETVAL' )THEN
            IACTIO= 1
         ELSEIF( TEXT12.EQ.'MAXVAL' )THEN
            IACTIO= 2
         ELSEIF( TEXT12.EQ.'MINVAL' )THEN
            IACTIO= 3
         ELSEIF( TEXT12.EQ.'INDMAX' )THEN
            IACTIO= 4
         ELSEIF( TEXT12.EQ.'INDMIN' )THEN
            IACTIO= 5
         ELSEIF( TEXT12.EQ.'MEAN'   )THEN
            IACTIO= 6
         ELSEIF( TEXT12.EQ.'TYPE'   )THEN
            IACTIO= 7
         ELSEIF( TEXT12.EQ.'LENGTH' )THEN
            IACTIO= 8
         ENDIF
*
*        FIND BLOCK NAME
         CALL REDGET(ITYP1,ISET ,FLOTT,HNAME ,DFLOTT)
         IF(ITYP1.EQ.1) THEN
            CALL LCMLEL(IPDATA,ISET  ,ILENG2,ITYLCM)
         ELSE IF(ITYP1.EQ.3) THEN
            CALL LCMLEN(IPDATA,HNAME ,ILENG2,ITYLCM)
         ELSE
            CALL XABORT('DRVGRP: BLOCK-NAME OR LIST INDEX EXPECTED.')
         ENDIF
         IF((ITYLCM.EQ.4).OR.(ITYLCM.EQ.6)) ILENG2=2*ILENG2
         IF( IACTIO.EQ.7 ) THEN
            NOUT= 1
            NSTP= 1
            ITYPE= 1
            ALLOCATE(JDATA(NOUT))
            JDATA(1)= ITYLCM
            CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
            GO TO 310
         ELSE IF( IACTIO.EQ.8 ) THEN
            NOUT= 1
            NSTP= 1
            ITYPE= 1
            ALLOCATE(JDATA(NOUT))
            JDATA(1)= ILENG2
            CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
            GO TO 310
         ENDIF
         IF( ILENG2.EQ.0 )THEN
            CALL LCMLIB(IPDATA)
            CALL XABORT('DRVGRP: BLOCK *'//HNAME//'* IS NOT STORED IN *'
     >                 //HLIST//'*.')
         ELSE IF( ITYLCM.EQ.10 ) THEN
            CALL XABORT('DRVGRP: '//HNAME//' IS A LIST OF ARRAYS. USE A'
     >                 //' STEP UP KEYWORD TO ACCESS THE LIST.')
         ENDIF
         ALLOCATE(IDATA(ILENG2))
         ALLOCATE(JDATA(ILENG2))
         IF( ITYLCM.EQ.3 )THEN
            ILENG= ILENG2*4
         ELSE
            ILENG= ILENG2
         ENDIF
*
*        GET BLOCK
         IF(ITYP1.EQ.1) THEN
            CALL LCMGDL(IPDATA,ISET  ,IDATA)
         ELSE IF(ITYP1.EQ.3) THEN
            CALL LCMGET(IPDATA,HNAME ,IDATA)
         ENDIF
*
         CALL REDGET(ITYP,I    ,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.1.OR.I.LT.1 )
     >      CALL XABORT('DRVGRP: POSITIVE INDEX EXPECTED.')
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         J= I
         K= 1
         N= 1
         IF( ITYP.EQ.1 )THEN
            J= NITMA
            IF( J.LT.I )
     >         CALL XABORT('DRVGRP: SECOND INDEX EXPECTED GREATER '
     >                   //'OR EQUAL THAN FIRST.')
            CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
            IF( ITYP.EQ.1 )THEN
               K= NITMA
               IF( K.LT.1 )
     >            CALL XABORT('DRVGRP: POSITIVE THIRD INDEX EXPECTED.')
               IF( MOD(J-I,K).NE.0 )
     >            CALL XABORT('DRVGRP: THIRD INDEX EXPECTED TO BALANCE'
     >                      //' STEPS FROM FIRST TO SECOND INDEX.')
               CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
            ENDIF
            N= (J-I)/K + 1
            IF( N.LT.1 )
     >         CALL XABORT('DRVGRP: INCONSISTENT NUMBER OF WORDS.')
         ELSEIF( (ITYP.EQ.3).AND.(TEXT12.EQ.'*') )THEN
            J= ILENG
            CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
            IF( ITYP.EQ.1 )THEN
               K= NITMA
               IF( K.LT.1 )
     >            CALL XABORT('DRVGRP: POSITIVE THIRD INDEX EXPECTED.')
               IF( MOD(J-I,K).NE.0 )
     >            CALL XABORT('DRVGRP: THIRD INDEX EXPECTED TO BALANCE'
     >                      //' STEPS FROM FIRST TO SECOND INDEX.')
               CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
            ENDIF
            N= (J-I)/K + 1
            IF( N.LT.1 )
     >         CALL XABORT('DRVGRP: INCONSISTENT NUMBER OF WORDS.')
         ENDIF
         IF( TEXT12.EQ.'NVAL' )THEN
            IF( N.NE.1.OR.K.NE.1 )THEN
               CALL XABORT('DRVGRP: NVALUE ALREADY GIVEN FROM INDEX.')
            ENDIF
            CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
            IF( ITYP.EQ.1)THEN
               N= NITMA
               IF( N.LT.1 )
     >         CALL XABORT('DRVGRP: POSITIVE NVALUE EXPECTED.')
               J= I + N - 1
            ELSEIF( (ITYP.EQ.3).AND.(TEXT12.EQ.'*') )THEN
               J= ILENG
               N= ILENG - I + 1
            ELSE
               CALL XABORT('DRVGRP: NVAL IS FOLLOWED BY * OR INTEGER')
            ENDIF
            CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         ENDIF
         IF( J.GT.ILENG )THEN
            WRITE(HSMG,'(29HDRVGRP: THE VALUE OF INDEX2 (,I8,7H) IS GR,
     >      29HEATER THAN THE BLOCK LENGTH (,I8,2H).)') J,ILENG
            CALL XABORT(HSMG)
         ENDIF
         GO TO 30
      ELSEIF( TEXT12.EQ.';' )THEN
         GO TO 40
      ELSE
         CALL XABORT('DRVGRP: '//TEXT12//' IS AN INVALID KEYWORD.')
      ENDIF
      GO TO 20
*----
*  PROCESSING THE COMMAND NUMBER: IACTIO.
*----
   30 CONTINUE
      NSTP= 1
      NOUT= 0
      ITYPE= ITYLCM
      IF(     ITYLCM.EQ.1 )THEN
*
*        GREP INTEGER DATA
         INTMIN= INTGMX
         INTMAX=-INTGMX
         DO 301 IOUT=I,J,K
            NOUT= NOUT+1
            ITRANS= IDATA(IOUT)
            JDATA(NOUT)= ITRANS
            IF( ITRANS.GT.INTMAX )THEN
               INTMAX= ITRANS
               INDMAX= IOUT
            ENDIF
            IF( ITRANS.LT.INTMIN )THEN
               INTMIN= ITRANS
               INDMIN= IOUT
            ENDIF
  301    CONTINUE
         IF(     IACTIO.EQ.2 )THEN
            NOUT= 1
            JDATA(1)= INTMAX
         ELSEIF( IACTIO.EQ.3 )THEN
            NOUT= 1
            JDATA(1)= INTMIN
         ELSEIF( IACTIO.EQ.4 )THEN
            NOUT= 1
            JDATA(1)= INDMAX
            ITYPE= 1
         ELSEIF( IACTIO.EQ.5 )THEN
            NOUT= 1
            JDATA(1)= INDMIN
            ITYPE= 1
         ELSEIF( IACTIO.NE.1 )THEN
            CALL XABORT('DRVGRP: INVALID ACTION ON INTEGERS')
         ENDIF
      ELSEIF( ITYLCM.EQ.2 )THEN
*
*        GREP REAL DATA
         RELMIN= REALMX
         RELMAX=-REALMX
         RELMNV= 0.0
         DO 302 IOUT=I,J,K
            NOUT= NOUT+1
            ITRANS= IDATA(IOUT)
            JDATA(NOUT)= ITRANS
            IF( ATRANS.GT.RELMAX )THEN
               RELMAX= ATRANS
               INDMAX= IOUT
            ENDIF
            IF( ATRANS.LT.RELMIN )THEN
               RELMIN= ATRANS
               INDMIN= IOUT
            ENDIF
            RELMNV= RELMNV+ATRANS
  302    CONTINUE
         IF(     IACTIO.EQ.2 )THEN
            NOUT= 1
            ATRANS= RELMAX
            JDATA(1)= ITRANS
         ELSEIF( IACTIO.EQ.3 )THEN
            NOUT= 1
            ATRANS= RELMIN
            JDATA(1)= ITRANS
         ELSEIF( IACTIO.EQ.4 )THEN
            NOUT= 1
            JDATA(1)= INDMAX
            ITYPE= 1
         ELSEIF( IACTIO.EQ.5 )THEN
            NOUT= 1
            JDATA(1)= INDMIN
            ITYPE= 1
         ELSEIF( IACTIO.EQ.6 )THEN
            ATRANS= RELMNV/FLOAT(NOUT)
            NOUT= 1
            JDATA(1)= ITRANS
         ENDIF
      ELSEIF( ITYLCM.EQ.3 )THEN
         IF( IACTIO.NE.1 )THEN
            CALL XABORT('DRVGRP: INVALID ACTION ON STRING')
         ELSEIF( (J-I)/K.GT.71 )THEN
            CALL XABORT('DRVGRP: STRING HAS LENGTH .GT. 72')
         ENDIF
         TEXT72= ' '
         IOFSET= 0
         IF( ILENG.GE.72 )THEN
            DO 313 IBCHAR= 1, ILENG/72
               WRITE(TEXTIN,'(18A4)') (IDATA(IOFSET+IOUT),IOUT=1,18)
               DO 303 IOUT= I,J,K
                  IF( IOUT.GT.NBCHAR.AND.IOUT.LE.NBCHAR+72 )THEN
                     NOUT= NOUT+1
                     TEXT72(NOUT:NOUT)= TEXTIN(IOUT-NBCHAR:IOUT-NBCHAR)
                  ENDIF
  303          CONTINUE
               IOFSET= IOFSET+18
               NBCHAR= NBCHAR+72
  313       CONTINUE
         ENDIF
         NRESID= (ILENG-NBCHAR)/4
         IF( NRESID.GT.0 )THEN
            WRITE(TEXTIN,'(18A4)') (IDATA(IOFSET+IOUT),IOUT=1,NRESID)
            DO 323 IOUT= I,J,K
               IF( IOUT.GT.NBCHAR.AND.IOUT.LE.NBCHAR+NRESID*4 )THEN
                  NOUT= NOUT+1
                  TEXT72(NOUT:NOUT)= TEXTIN(IOUT-NBCHAR:IOUT-NBCHAR)
               ENDIF
  323       CONTINUE
         ENDIF
         NBCHAR= NOUT
         NOUT= 1
      ELSEIF( ITYLCM.EQ.4 )THEN
*
*        GREP DOUBLE PRECISION DATA
         I= I+I
         J= J+J
         K= K+K
         DBLMIN= DBLEMX
         DBLMAX=-DBLEMX
         DBLMNV= 0.0D0
         DO 304 IOUT=I,J,K
            NOUT= NOUT+1
            ISEEME(1)= IDATA(IOUT-1)
            ISEEME(2)= IDATA(IOUT)
            JDATA(2*NOUT-1)= ISEEME(1)
            JDATA(2*NOUT)= ISEEME(2)
            IF( DSEEME.GT.DBLMAX )THEN
               DBLMAX= DSEEME
               INDMAX= IOUT
            ENDIF
            IF( DSEEME.LT.DBLMIN )THEN
               DBLMIN= DSEEME
               INDMIN= IOUT
            ENDIF
            DBLMNV= DBLMNV+DSEEME
  304    CONTINUE
         IF(     IACTIO.EQ.2 )THEN
            NOUT= 1
            DSEEME= DBLMAX
            JDATA(1)= ISEEME(1)
            JDATA(2)= ISEEME(2)
         ELSEIF( IACTIO.EQ.3 )THEN
            NOUT= 1
            DSEEME= DBLMIN
            JDATA(1)= ISEEME(1)
            JDATA(2)= ISEEME(2)
         ELSEIF( IACTIO.EQ.4 )THEN
            NOUT= 1
            JDATA(1)= INDMAX
            ITYPE= 1
         ELSEIF( IACTIO.EQ.5 )THEN
            NOUT= 1
            JDATA(1)= INDMIN
            ITYPE= 1
         ELSEIF( IACTIO.EQ.6 )THEN
            DSEEME= DBLMNV/DBLE(NOUT)
            NOUT= 1
            JDATA(1)= ISEEME(1)
            JDATA(2)= ISEEME(2)
         ENDIF
         IF( ITYPE.EQ.4 )THEN
            NSTP= 2
            NOUT= 2*NOUT-1
         ENDIF
      ELSEIF( ITYLCM.EQ.5 )THEN
*
*        GREP LOGICAL DATA
         IF( IACTIO.NE.1 )THEN
            CALL XABORT('DRVGRP: INVALID ACTION ON LOGICALS')
         ENDIF
         DO 305 IOUT=I,J,K
            NOUT= NOUT+1
            JDATA(NOUT)= IDATA(IOUT)
  305    CONTINUE
      ELSE
         CALL XABORT('DRVGRP: INVALID DATA TYPE.')
      ENDIF
      DEALLOCATE(IDATA)
*----
*  PUT CLE-2000 PARMS IN CREATE OR READ/WRITE MODES.
*----
  310 DO 35 IOUT= 1, NOUT, NSTP
         IF( -ITYP.NE.ITYPE )THEN
            CALL XABORT('DRVGRP: NOT ENOUGH CLE-2000 PARAMETERS '
     >                //'TO CONTAIN ALL VALUES ASKED TO BE PICKED')
         ENDIF
         ITYP= ITYPE
         IF( ITYP.EQ.1.OR.ITYP.EQ.5 )THEN
            NITMA= JDATA(IOUT)
            IF( IPRINT.GT.0 )THEN
               IF( ITYP.EQ.1 )THEN
                  WRITE(6,*) CLOGBG(ITYP),TEXT12,'<-',NITMA
               ELSE
                  IF( NITMA.EQ.+1 )THEN
                     WRITE(6,*) CLOGBG(ITYP),TEXT12,'<- $True_L'
                  ELSEIF( NITMA.EQ.-1 )THEN
                     WRITE(6,*) CLOGBG(ITYP),TEXT12,'<- $False_L'
                  ELSE
                     WRITE(6,*) CLOGBG(ITYP),TEXT12,'<- ?_L'
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF( ITYP.EQ.2 )THEN
            ITRANS= JDATA(IOUT)
            FLOTT= ATRANS
            IF( IPRINT.GT.0 )THEN
               WRITE(6,*)CLOGBG(ITYP),TEXT12,'<-',FLOTT
            ENDIF
         ELSEIF( ITYP.EQ.3 )THEN
            NITMA=NBCHAR
            IF( IPRINT.GT.0 )THEN
               WRITE(6,*)CLOGBG(ITYP),TEXT12,'<-"',TEXT72(1:NBCHAR),'"'
            ENDIF
         ELSEIF( ITYP.EQ.4 )THEN
            ISEEME(1)= JDATA(IOUT)
            ISEEME(2)= JDATA(IOUT+1)
            DFLOTT= DSEEME
            IF( IPRINT.GT.0 )THEN
               WRITE(6,*)CLOGBG(ITYP),TEXT12,'<-',DFLOTT
            ENDIF
         ENDIF
         CALL REDPUT(ITYP,NITMA,FLOTT,TEXT72,DFLOTT)
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
   35 CONTINUE
      DEALLOCATE(JDATA)
      GO TO 25
*----
*  ENDING COMMANDS: CHECK UP/DOWN BALANCE AND REMAINING PARMS.
*----
   40 CONTINUE
      RETURN
      END
